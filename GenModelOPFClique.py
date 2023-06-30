from pyomo.environ import ConcreteModel, RangeSet, Var, Reals, ConstraintList, Objective, Param, minimize, SolverFactory, value
from RecupDataMat import parse_matpower_file
from RecupDataGOC import parse_GOC_file
import numpy as np
from ChordalNx import find_clique_tree_from_caseMatpower, find_clique_tree_from_caseGOC
from modelOPF import construct_Y_matrix
import time

# Cherches la liste des arrêtes en utilsant l'indices dans cliques pour chaque arrête du clique tree
def edge_ind_cliques(cliques, edges):
    liste_aretes = []
    for tuple_indices in edges:
        clique_indices, clique_indices_suivant = map(list, tuple_indices)
        arete = (cliques.index(clique_indices) + 1, cliques.index(clique_indices_suivant) + 1)
        liste_aretes.append(arete)
    return liste_aretes

# Cherche le nombre d'apparition d'un noeud dans une clique
def number_i_cliques(cliques, n_bus):
    liste_number = []
    for i in range(n_bus):
        j = 0
        for k in range(len(cliques)):
            if i in cliques[k]:
                j += 1
        liste_number.append(j)
    return liste_number

# Recherche une clique contenant les noeuds i et j
def clique_with_var1_var2(cliques,i,j):
    min = 1e10
    ind = None
    #print(i)
    #print(j)
    for k in range(len(cliques)):
        if i in cliques[k] and j in cliques[k] and len(cliques[k])<min:
            min = len(cliques[k])
            ind = k+1
    return ind 

# Recherche une clique contenant le noeuds i
def clique_with_var(cliques,i):
    min = 1e10
    ind = None
    for k in range(len(cliques)):
        if i in cliques[k] and len(cliques[k])<min:
            min = len(cliques[k])
            ind = k+1    
    return ind    


# Fonction qui crée à partir d'un problème OPF contenues dans ppc, d'une liste de cliques et du cliques tree associé un model définissant le problème OPF 
def GenSepModelOPFR(ppc, cliques, clique_tree_edges):

    # Récupération des données de ppc
    baseMVA = float(ppc['baseMVA'])
    bus = np.array(ppc['bus'])
    gen = np.array(ppc['gen'])
    branch = np.array(ppc['branch'])
    gencost = np.array(ppc['gencost'])
    
    #création de la matrice d'admittance
    Ybranch = construct_Y_matrix(branch)

    #Initialisation du model
    model = ConcreteModel()
    model.BaseMVA = baseMVA
    model.nb_bus = len(bus)
    model.nb_gen = len(gen)
    model.nb_branch = len(branch)
    model.nb_cost = len(gencost)
    model.nb_cliques = len(cliques)

    model.BUS = RangeSet(1, model.nb_bus)
    model.GEN = RangeSet(1, model.nb_gen)
    model.BRANCH = RangeSet(1, model.nb_branch)
    model.GENCOST = RangeSet(1, model.nb_cost)
    model.CLIQUES = RangeSet(1, model.nb_cliques)

    iBus = {}
    for i in model.BUS:
        iBus[bus[i-1,0]]= i

    # Paramètres des demandes de puissance active et réactive
    model.Pd = Param(model.BUS, initialize={i: bus[i-1, 2] for i in model.BUS})
    model.Qd = Param(model.BUS, initialize={i: bus[i-1, 3] for i in model.BUS})

    # Paramètres de limites de tension
    model.Vmax = Param(model.BUS, initialize={i: bus[i-1, 11] for i in model.BUS})
    model.Vmin = Param(model.BUS, initialize={i: bus[i-1, 12] for i in model.BUS})

    # Paramètre b pour l'équation d'équilibre de puissance réactive
    model.b = Param(model.BUS, initialize={i: bus[i-1, 5] for i in model.BUS})
    
    # Paramètre g pour l'équation d'équilibre de puissance active
    model.g = Param(model.BUS, initialize={i: bus[i-1, 4] for i in model.BUS})

    #Création des Variables du problème

    # Définir les variables pour chaque clique
    model.Vre = Var(model.CLIQUES, model.BUS, within=Reals, initialize=lambda model, j, i: bus[i-1, 7], name='Vre')
    model.Vim = Var(model.CLIQUES, model.BUS, within=Reals, initialize=0.0)

    # initialisation de la valeur objective
    model.val_obj = 0    

    #Contraintes liant les valeurs des memes noeuds dans différente cliques lié par le clique_tree
    model.similNode = ConstraintList()
    edges_ind = edge_ind_cliques(cliques,clique_tree_edges)
    for edge in edges_ind:
        common_elements = list(set(cliques[edge[0]-1]).intersection(cliques[edge[1]-1]))
        for i in common_elements:
            model.similNode.add(model.Vre[edge[0], i+1] == model.Vre[edge[1], i+1])
            model.similNode.add(model.Vim[edge[0], i+1] == model.Vim[edge[1], i+1])
            

    #Contraintes de puissance active pour chaque noeud en utilisant la forme rectangulaire
    model.ActivePowerBalance = ConstraintList()
    
    for i in model.BUS:
        
        i_is_Gen = False
        ind_Gen = None

        # demande de puissance active au noeud i
        aPd = model.Pd[i]
        
        #for k in model.BRANCH:
            #if iBus[branch[k-1,0]] == i:
                #print(bus[0])
                #print((i,branch[k-1,0], branch[k-1,1],bus[int(iBus[branch[k - 1, 0]]-1)]))


        # Définition de Porig et Pdest à partire de la matrice d'admittance et des tensions au différents noeuds en veillant à chercher des cliques qui contiennent les noeuds de la branche
        Porig = sum(
            ((Ybranch[k-1, 0, 0].real * model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]] - Ybranch[k-1, 0, 0].imag * model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]]
            + Ybranch[k-1, 0, 1].real * model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]] - Ybranch[k-1, 0, 1].imag * model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]])*model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]]
            +
            (Ybranch[k-1, 0, 0].imag * model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]] + Ybranch[k-1, 0, 0].real * model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]]
            + Ybranch[k-1, 0, 1].imag * model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]] + Ybranch[k-1, 0, 1].real * model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]])*model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]])
            for k in model.BRANCH if iBus[branch[k - 1, 0]] == i )
        
        Pdest = sum(((Ybranch[k-1, 1, 0].real * model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]] - Ybranch[k-1, 1, 0].imag * model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]]
            + Ybranch[k-1, 1, 1].real * model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]] - Ybranch[k-1, 1, 1].imag * model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]])*model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]]
            +
            (Ybranch[k-1, 1, 0].imag * model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]] + Ybranch[k-1, 1, 0].real * model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]]
            + Ybranch[k-1, 1, 1].imag * model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]] + Ybranch[k-1, 1, 1].real * model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]])*model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]]) 
            for k in model.BRANCH if iBus[branch[k - 1, 1]] == i )
        
        # on vérifie dans la boucle suivante si le noeud i fait partie des générateurs
        for k in model.GEN:
            if iBus[gen[k-1,0]] == i:
                i_is_Gen = True
                ind_Gen = k-1

        if i_is_Gen:
            
            # si le noeud est générateurs on définit la contraintes de flux de puissance actives au noeud et établit les bornes de la puissance active générée en ce noeud           
            #en veillant à chercher une clique contenant le noeud
            model.ActivePowerBalance.add(
                aPd+baseMVA*(Porig+Pdest)+ bus[i-1, 4] * (model.Vim[clique_with_var(cliques,i-1), i]**2 + model.Vre[clique_with_var(cliques,i-1), i]**2) <= gen[ind_Gen, 8]
                )
            model.ActivePowerBalance.add(
                aPd+baseMVA*(Porig+Pdest)+ bus[i-1, 4] * (model.Vim[clique_with_var(cliques,i-1), i]**2 + model.Vre[clique_with_var(cliques,i-1), i]**2) >= gen[ind_Gen, 9]
            )

            # Puis on ajoute [Cg * Pg] à la fonction objective on retirant le coefficient quadratique lorsqu'il existe
            if gencost[ind_Gen, 3] == 3:
                model.val_obj += (
                        gencost[ind_Gen, 5] * (aPd + baseMVA*(Porig+Pdest)+bus[i-1,4]*(model.Vim[clique_with_var(cliques,i-1), i]**2 + model.Vre[clique_with_var(cliques,i-1), i]**2))
                        + gencost[ind_Gen, 6] 
                )
            else:
                # Sinon la valeur de puisssance active générée au noeud i est nulle
                model.val_obj +=  (
                                gencost[ind_Gen, 4] * (aPd + baseMVA*(Porig+Pdest)+bus[i-1,4]*(model.Vim[clique_with_var(cliques,i-1), i]**2 + model.Vre[clique_with_var(cliques,i-1), i]**2))
                                + gencost[ind_Gen, 5] 
                            )

        else:
            # Sinon la valeur de puisssance active générée au noeud i est nulle

            # code qui permet de verifier que la contrainte n'est pas vide
            #Ajout = sum(1 for k in model.BRANCH if branch[k - 1, 0] == i )
            #Ajout += sum(1 for k in model.BRANCH if branch[k - 1, 1] == i )
            #if Ajout>0 or bus[i-1,4]>0.0:

            model.ActivePowerBalance.add(
                (aPd+baseMVA*(Porig+Pdest)+ bus[i-1, 4] * (model.Vim[clique_with_var(cliques,i-1), i]**2 + model.Vre[clique_with_var(cliques,i-1), i]**2)) == 0
            )
            
              
        

    # Contraintes de puissance réactive pour chaque noeud en utilisant la forme rectangulaire
    model.ReactivePowerBalance = ConstraintList()
    for i in model.BUS:
    
        i_is_Gen = False
        ind_Gen = None
        
        # demande de puissance réactive au noeud i       
        aQd = model.Qd[i]

        # Définition de Qorig et Qdest à partire de la matrice d'admittance et des tensions au différents noeuds en veillant à chercher des cliques qui contiennent les noeuds de la branche
        Qorig = sum(
            ((Ybranch[k-1, 0, 0].real * model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]] - Ybranch[k-1, 0, 0].imag * model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]]
            + Ybranch[k-1, 0, 1].real * model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]] - Ybranch[k-1, 0, 1].imag * model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]])*model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]]
            -
            (Ybranch[k-1, 0, 0].imag * model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]] + Ybranch[k-1, 0, 0].real * model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]]
            + Ybranch[k-1, 0, 1].imag * model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]] + Ybranch[k-1, 0, 1].real * model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]])*model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]])  
            for k in model.BRANCH if iBus[branch[k - 1, 0]] == i )
        Qdest = sum(
            ((Ybranch[k-1, 1, 0].real * model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]] - Ybranch[k-1, 1, 0].imag * model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]]
            + Ybranch[k-1, 1, 1].real * model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]] - Ybranch[k-1, 1, 1].imag * model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]])*model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]]
            -
            (Ybranch[k-1, 1, 0].imag * model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]] + Ybranch[k-1, 1, 0].real * model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,0]]]
            + Ybranch[k-1, 1, 1].imag * model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]] + Ybranch[k-1, 1, 1].real * model.Vim[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]])*model.Vre[clique_with_var1_var2(cliques,iBus[branch[k-1,0]]-1,iBus[branch[k-1,1]]-1), iBus[branch[k-1,1]]])
            for k in model.BRANCH if iBus[branch[k - 1, 1]] == i )

        # on vérifie dans la boucle suivante si le noeud i fait partie des générateurs
        for k in model.GEN:
            if iBus[gen[k-1,0]]+1 == i:
                i_is_Gen = True
                ind_Gen = k-1
        
        if i_is_Gen:
            # si le noeud est générateurs on définit la contraintes de flux de puissance réactives au noeud et établit les bornes de la puissance réactive générée en ce noeud
            # en veillant à chercher une clique contenant le noeud
            model.ReactivePowerBalance.add(
                aQd + baseMVA*(Qorig + Qdest)-(bus[i-1, 5]*(model.Vim[clique_with_var(cliques,i-1), i]**2 + model.Vre[clique_with_var(cliques,i-1), i]**2))   <= gen[ind_Gen, 3]
            )
            model.ReactivePowerBalance.add(
                aQd + baseMVA*(Qorig + Qdest)-(bus[i-1, 5]*(model.Vim[clique_with_var(cliques,i-1), i]**2 + model.Vre[clique_with_var(cliques,i-1), i]**2))   >= gen[ind_Gen, 4]
            )
        else:
            # Sinon la valeur de puisssance réactive générée au noeud i est nulle
            
            # Code qui permet de verifier que la contrainte n'est pas vide
            
            #Ajout = sum(1 for k in model.BRANCH if branch[k - 1, 0] == i )
            #Ajout += sum(1 for k in model.BRANCH if branch[k - 1, 1] == i )
            #if Ajout>0 or bus[i-1,5]>0.0:
            
            model.ReactivePowerBalance.add(
                (aQd + baseMVA*(Qorig + Qdest)-(bus[i-1, 5]*(model.Vim[clique_with_var(cliques,i-1), i]**2 + model.Vre[clique_with_var(cliques,i-1), i]**2))) >= 0
            )
            
             
    # Définition des bornes du module de la tension
    model.VoltageBounds = ConstraintList()
    for clique in model.CLIQUES:
        for i in cliques[clique-1]:
            model.VoltageBounds.add(model.Vmin[i+1] ** 2 <= (model.Vre[clique, i+1] ** 2 + model.Vim[clique, i+1] ** 2))
            model.VoltageBounds.add((model.Vre[clique, i+1] ** 2 + model.Vim[clique, i+1] ** 2) <= model.Vmax[i+1] ** 2)

    
    # Définition des Limites de courrants sur une branche pour chaque clique
    model.BranchCurrentBalance = ConstraintList()
    for clique in model.CLIQUES:
        for k in model.BRANCH:
            orig = iBus[branch[k - 1, 0]]
            dest = iBus[branch[k - 1, 1]]
            Imax = branch[k - 1, 5] / baseMVA

            if ((orig-1)in cliques[clique-1])and ((dest-1)in cliques[clique-1]) and Imax != 0:
                # Si une limite à été définit : on établit les limites du courrant en utilsant [I = Y.V]
                
                # Limite pour Iorig
                model.BranchCurrentBalance.add(
                    ((Ybranch[k-1, 0, 0].real * model.Vre[clique, orig] - Ybranch[k-1, 0, 0].imag * model.Vim[clique, orig]
                    + Ybranch[k-1, 0, 1].real * model.Vre[clique, dest] - Ybranch[k-1, 0, 1].imag * model.Vim[clique, dest])**2
                    +
                    (Ybranch[k-1, 0, 0].imag * model.Vre[clique, orig] + Ybranch[k-1, 0, 0].real * model.Vim[clique, orig]
                    + Ybranch[k-1, 0, 1].imag * model.Vre[clique, dest] + Ybranch[k-1, 0, 1].real * model.Vim[clique, dest])**2)
                    <= Imax**2
                )
                # Limite pour Idest
                model.BranchCurrentBalance.add(
                    ((Ybranch[k-1, 1, 0].real * model.Vre[clique, orig] - Ybranch[k-1, 1, 0].imag * model.Vim[clique, orig]
                    + Ybranch[k-1, 1, 1].real * model.Vre[clique, dest] - Ybranch[k-1, 1, 1].imag * model.Vim[clique, dest])**2
                    +
                    (Ybranch[k-1, 1, 0].imag * model.Vre[clique, orig] + Ybranch[k-1, 1, 0].real * model.Vim[clique, orig]
                    + Ybranch[k-1, 1, 1].imag * model.Vre[clique, dest] + Ybranch[k-1, 1, 1].real * model.Vim[clique, dest])**2)
                    <=Imax**2
                )
    # Définition de la fonction objective
    model.obj = Objective(expr=model.val_obj, sense=minimize)

    return model

def test_caseMatPower(filepath):
    start = time.time()
    case = parse_matpower_file(filepath)
    end = time.time()
    print(" le temps utilisé pour lire le fichier et crééle dictionnaire contenant les données est de : ", end-start)
    tot_t = end-start
    start = time.time()
    ch_mat, cliques, clique_tree_edges = find_clique_tree_from_caseMatpower(filepath)
    end = time.time()
    print(" le temps utilisé pour obtenir le graph chordale, la liste de clique et l'arbre  est de : ", end-start)
    tot_t += (end-start)
    #print(cliques," ", clique_tree_edges)
    #print("--------------------------------------clique trouver--------------------------------------------------")
    start = time.time()
    model = GenSepModelOPFR(case, cliques, clique_tree_edges)
    end = time.time()
    print(" le temps utilisé pour créé le model(cliques) est de : ", end-start)
    tot_t += (end-start)
    #print("-----------------------------------------model creer--------------------------------------------------")

    #print(cliques)
    # affichage du model (Variables , Param, Contraintes et fonction objective)
    #model.pprint()

    # Création de l'objet solveur
    solver = SolverFactory('ipopt', executable='/Users/djelelbenabdellah/anaconda3/envs/opf-env/bin/ipopt')

    # Définir le nombre maximum d'itération de Ipopt
    solver.options['max_iter'] = 150000

    # Résolution du problème d'optimisation sans affichage
    start = time.time()
    results = solver.solve(model)
    end = time.time()
    print(" le temps mis par le solveur pour résoudre le problème est de : ", end-start)
    tot_t += (end-start)
    print("le temps total est de : ", tot_t)
    # Affichage des valeurs des variables
    #for v in model.component_data_objects(Var, active=True):
    #    print(f"{v} = {value(v)}")

    # Affichage de la fonction objectif
    print(f"Fonction objective = {value(model.obj)}")

def test_caseGOC(filepath):
    start = time.time()
    case = parse_GOC_file(filepath)
    end = time.time()
    print(" le temps utilisé pour lire le fichier et crééle dictionnaire contenant les données est de : ", end-start)
    tot_t = (end-start)

    #print(case['bus'])
    start = time.time()
    ch_mat, cliques, clique_tree_edges = find_clique_tree_from_caseGOC(filepath)
    end = time.time()
    print(" le temps utilisé pour obtenir le graph chordale, la liste de clique et l'arbre  est de : ", end-start)
    tot_t += (end-start)

    #print("-------------------------------------- clique trouver ------------------------------------------------------------------")
    start = time.time()
    model = GenSepModelOPFR(case, cliques, clique_tree_edges)
    end = time.time()
    print(" le temps utilisé pour créé le model(cliques) est de : ", end-start)
    tot_t += (end-start)

    #print(cliques)
    #model.pprint()

    # Création de l'objet solveur
    solver = SolverFactory('ipopt', executable='/Users/djelelbenabdellah/anaconda3/envs/opf-env/bin/ipopt')

    # Définir le nombre maximum d'itération de Ipopt
    solver.options['max_iter'] = 150000

    # Résolution du problème d'optimisation sans affichage
    start = time.time()
    results = solver.solve(model)
    end = time.time()
    print(" le temps mis par le solveur pour résoudre le problème est de : ", end-start)
    tot_t += (end-start)
    print("le temps total est de : ", tot_t)

    # Affichage des valeurs des variables
    #for v in model.component_data_objects(Var, active=True):
    #    print(f"{v} = {value(v)}")

    # Affichage de la fonction objectif
    print(f"Fonction objective = {value(model.obj)}")

# Tests matpower
print("-----------------------------------------------------------------------------------------------------------------------")

print("\n case 9 :")
test_caseMatPower('cases/case9.m')
print("\n case 14 :")
test_caseMatPower('cases/case14.m')
print("\n case 30 :")
test_caseMatPower('cases/case30.m')
print("\n case 39 :")
test_caseMatPower('cases/case39.m')
print("\n case 57 :")
test_caseMatPower('cases/matpower/case57.m')
print("\n case 145 :")
test_caseMatPower('cases/matpower/case145.m')
#print("\n case 300 :")
#test_caseMatPower('cases/case300.m')

print("-----------------------------------------------------------------------------------------------------------------------")


# Test GOC
#print("test goc:")
#test_caseGOC('scenario/scenario_035/case.raw')
