from pyomo.environ import ConcreteModel, RangeSet, Var, Reals, ConstraintList, Objective, Param, minimize, SolverFactory, value
import pandapower.networks as pn
from RecupDataMat import parse_matpower_file
from RecupDataGOC import parse_GOC_file
import numpy as np
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
import time



# fonction qui définit la matrice d'admittance
def construct_Y_matrix(branch):
    # nombre de liens dans le modele
    num_br = len(branch)
    
    # initialisation de la matrice
    Y = np.zeros((num_br, 2, 2), dtype=np.complex128)
    
    # Boucle qui complète pour chacune des branches du problème la matrice d'admittance
    for i in range(num_br):
        
        yl = 1/(branch[i,2] + 1j*branch[i,3])
        bl = branch[i,4] / 2 
        tau_l = branch[i,8]
        theta_l = branch[i,9]

        if tau_l == 0:
            tau_l = 1

        Y[i, 0, 0] = (yl + 1j * bl) / (tau_l**2)
        Y[i, 1, 1] = (yl + 1j * bl )
        Y[i, 0, 1] = (-yl) / (tau_l * np.exp(-1j * theta_l))
        Y[i, 1, 0] = (-yl) / (tau_l * np.exp(1j * theta_l))

    return Y



# Fonction qui crée à partir d'un problème OPF de matpower chargé dans ppc un model pyomo D'oPF en nombre réels
def GenModelOPFR(ppc):

    # récupération des données du dictionnaire contenant les informations du model
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

    model.BUS = RangeSet(1, model.nb_bus)
    model.GEN = RangeSet(1, model.nb_gen)
    model.BRANCH = RangeSet(1, model.nb_branch)
    model.GENCOST = RangeSet(1, model.nb_cost)

    ind_Bus = {}
    for i in model.BUS:
        ind_Bus[bus[i-1,0]]= i

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
    #model.Vre = Var(model.BUS, within=Reals, initialize={i: bus[i-1, 7] for i in model.BUS})
    model.Vre = Var(model.BUS, within=Reals, initialize=1.0)
    model.Vim = Var(model.BUS, within=Reals, initialize=0.0)

    # initialisation de la valeur objective
    model.val_obj = 0
    

    #Contraintes de puissance active pour chaque noeud en utilisant la forme rectangulaire
    model.ActivePowerBalance = ConstraintList()
    for i in model.BUS:

        i_is_Gen = False
        ind_Gen = None

        # demande de puissance active au noeud i
        aPd = model.Pd[i]
        
        # Définition de Porig et Pdest à partire de la matrice d'admittance et des tensions au différents noeuds
        Porig = sum(
            ((Ybranch[k-1, 0, 0].real * model.Vre[ind_Bus[branch[k-1,0]]] - Ybranch[k-1, 0, 0].imag * model.Vim[ind_Bus[branch[k-1,0]]]
            + Ybranch[k-1, 0, 1].real * model.Vre[ind_Bus[branch[k-1,1]]] - Ybranch[k-1, 0, 1].imag * model.Vim[ind_Bus[branch[k-1,1]]])*model.Vre[ind_Bus[branch[k-1,0]]]
            +
            (Ybranch[k-1, 0, 0].imag * model.Vre[ind_Bus[branch[k-1,0]]] + Ybranch[k-1, 0, 0].real * model.Vim[ind_Bus[branch[k-1,0]]]
            + Ybranch[k-1, 0, 1].imag * model.Vre[ind_Bus[branch[k-1,1]]] + Ybranch[k-1, 0, 1].real * model.Vim[ind_Bus[branch[k-1,1]]])*model.Vim[ind_Bus[branch[k-1,0]]])
             for k in model.BRANCH if branch[k - 1, 0] == bus[i - 1, 0])
        
        Pdest = sum(((Ybranch[k-1, 1, 0].real * model.Vre[ind_Bus[branch[k-1,0]]] - Ybranch[k-1, 1, 0].imag * model.Vim[ind_Bus[branch[k-1,0]]]
              + Ybranch[k-1, 1, 1].real * model.Vre[ind_Bus[branch[k-1,1]]] - Ybranch[k-1, 1, 1].imag * model.Vim[ind_Bus[branch[k-1,1]]])*model.Vre[ind_Bus[branch[k-1,1]]]
            +
            (Ybranch[k-1, 1, 0].imag * model.Vre[ind_Bus[branch[k-1,0]]] + Ybranch[k-1, 1, 0].real * model.Vim[ind_Bus[branch[k-1,0]]]
             + Ybranch[k-1, 1, 1].imag * model.Vre[ind_Bus[branch[k-1,1]]] + Ybranch[k-1, 1, 1].real * model.Vim[ind_Bus[branch[k-1,1]]])*model.Vim[ind_Bus[branch[k-1,1]]]) 
            for k in model.BRANCH if branch[k - 1, 1] == bus[i - 1, 0])

        # on vérifie dans la boucle suivante si le noeud i fait partie des générateurs
        for k in model.GEN:
            if ind_Bus[gen[k-1,0]] == i:
                i_is_Gen = True 
                ind_Gen = k-1

        if i_is_Gen:

            # si le noeud est générateurs on définit la contraintes de flux de puissance actives au noeud et établit les bornes de la puissance active générée en ce noeud
            model.ActivePowerBalance.add(
                aPd+baseMVA*(Porig + Pdest) + bus[i - 1, 4] * (model.Vim[i] ** 2 + model.Vre[i] ** 2) <= gen[ind_Gen,8]
            )
            model.ActivePowerBalance.add(
                aPd+baseMVA*(Porig + Pdest)  + bus[i - 1, 4] * (model.Vim[i] ** 2 + model.Vre[i] ** 2) >= gen[ind_Gen,9] 
            )
            # Puis on ajoute [Cg * Pg] à la fonction objective on retirant le coefficient quadratique lorsqu'il existe
            if(gencost[ind_Gen,3] == 3):
                model.val_obj += (
                              gencost[ind_Gen, 5] * (aPd+baseMVA*(Porig + Pdest )+ bus[i - 1, 4] * (model.Vim[i] ** 2 + model.Vre[i] ** 2)) +
                              gencost[ind_Gen, 6])
            else:
                model.val_obj += (
                              gencost[ind_Gen, 4] * (aPd+baseMVA*(Porig + Pdest )+ bus[i - 1, 4] * (model.Vim[i] ** 2 + model.Vre[i] ** 2)) +
                              gencost[ind_Gen, 5])
        else:
            # Sinon la valeur de puisssance active générée au noeud i est nulle
            model.ActivePowerBalance.add(
                aPd+ baseMVA*(Porig + Pdest )+ bus[i - 1, 4] * (model.Vim[i] ** 2 + model.Vre[i] ** 2) == 0
            )

        

    # Contraintes de puissance réactive pour chaque noeud en utilisant la forme rectangulaire
    model.ReactivePowerBalance = ConstraintList()
    for i in model.BUS:
        
        i_is_Gen = False
        ind_Gen = None

        # demande de puissance réactive au noeud i
        aQd = model.Qd[i]

        # Définition de Qorig et Qdest à partire de la matrice d'admittance et des tensions au différents noeuds
        Qorig = sum(
            ((Ybranch[k-1, 0, 0].real * model.Vre[ind_Bus[branch[k-1,0]]] - Ybranch[k-1, 0, 0].imag * model.Vim[ind_Bus[branch[k-1,0]]]
            + Ybranch[k-1, 0, 1].real * model.Vre[ind_Bus[branch[k-1,1]]] - Ybranch[k-1, 0, 1].imag * model.Vim[ind_Bus[branch[k-1,1]]])*model.Vim[ind_Bus[branch[k-1,0]]]
            -
            (Ybranch[k-1, 0, 0].imag * model.Vre[ind_Bus[branch[k-1,0]]] + Ybranch[k-1, 0, 0].real * model.Vim[ind_Bus[branch[k-1,0]]]
            + Ybranch[k-1, 0, 1].imag * model.Vre[ind_Bus[branch[k-1,1]]] + Ybranch[k-1, 0, 1].real * model.Vim[ind_Bus[branch[k-1,1]]])*model.Vre[ind_Bus[branch[k-1,0]]])  
            for k in model.BRANCH if branch[k - 1, 0] == bus[i - 1, 0])
        Qdest = sum(
            ((Ybranch[k-1, 1, 0].real * model.Vre[ind_Bus[branch[k-1,0]]] - Ybranch[k-1, 1, 0].imag * model.Vim[ind_Bus[branch[k-1,0]]]
              + Ybranch[k-1, 1, 1].real * model.Vre[ind_Bus[branch[k-1,1]]] - Ybranch[k-1, 1, 1].imag * model.Vim[ind_Bus[branch[k-1,1]]])*model.Vim[ind_Bus[branch[k-1,1]]]
            -
            (Ybranch[k-1, 1, 0].imag * model.Vre[ind_Bus[branch[k-1,0]]] + Ybranch[k-1, 1, 0].real * model.Vim[ind_Bus[branch[k-1,0]]]
             + Ybranch[k-1, 1, 1].imag * model.Vre[ind_Bus[branch[k-1,1]]] + Ybranch[k-1, 1, 1].real * model.Vim[ind_Bus[branch[k-1,1]]])*model.Vre[ind_Bus[branch[k-1,1]]])
             for k in model.BRANCH if branch[k - 1, 1] == bus[i - 1, 0])
        
        # on vérifie dans la boucle suivante si le noeud i fait partie des générateurs
        for k in model.GEN:
            if ind_Bus[gen[k-1,0]] == i:
                i_is_Gen = True 
                ind_Gen = k-1
        
        if i_is_Gen:

            # si le noeud est générateurs on définit la contraintes de flux de puissance réactives au noeud et établit les bornes de la puissance réactive générée en ce noeud
            model.ReactivePowerBalance.add(
                aQd+baseMVA*( Qorig + Qdest)- (bus[i - 1, 5] * (model.Vim[i] ** 2 + model.Vre[i] ** 2)) <= gen[ind_Gen, 3]
            )
            model.ReactivePowerBalance.add(
                aQd+baseMVA*(Qorig + Qdest )- (bus[i - 1, 5] * (model.Vim[i] ** 2 + model.Vre[i] ** 2)) >= gen[ind_Gen, 4]
            )
        else:
            # Sinon la valeur de puisssance réactive générée au noeud i est nulle
            model.ReactivePowerBalance.add(
                aQd+baseMVA*(Qorig + Qdest )- (bus[i - 1, 5] * (model.Vim[i] ** 2 + model.Vre[i] ** 2)) == 0
            )
         

    # Définition des Bornes du module de la tension
    model.VoltageBounds = ConstraintList()
    for i in model.BUS:
        
        # Borne inférieur
        model.VoltageBounds.add(model.Vmin[i] ** 2 <= (model.Vre[i] ** 2 + model.Vim[i] ** 2))
        
        # Borne Supérieur
        model.VoltageBounds.add((model.Vre[i] ** 2 + model.Vim[i] ** 2) <= model.Vmax[i] ** 2)
        
    # Définition des Limites de courrants sur une branche
    model.BranchCurrentBalance = ConstraintList()
    for k in model.BRANCH:
        orig = branch[k - 1, 0]
        dest = branch[k - 1, 1]
        rateA = branch[k-1, 5]

        Imax = rateA / baseMVA
        


        if Imax > 1e-10:
            # Si une limite à été définit : on établit les limites du courrant en utilsant [I = Y.V]

            # Limite pour Iorig
            model.BranchCurrentBalance.add(
                ((Ybranch[k-1, 0, 0].real * model.Vre[ind_Bus[orig]] - Ybranch[k-1, 0, 0].imag * model.Vim[ind_Bus[orig]]
                + Ybranch[k-1, 0, 1].real * model.Vre[ind_Bus[dest]] - Ybranch[k-1, 0, 1].imag * model.Vim[ind_Bus[dest]])**2
                +
                (Ybranch[k-1, 0, 0].imag * model.Vre[ind_Bus[orig]] + Ybranch[k-1, 0, 0].real * model.Vim[ind_Bus[orig]]
                + Ybranch[k-1, 0, 1].imag * model.Vre[ind_Bus[dest]] + Ybranch[k-1, 0, 1].real * model.Vim[ind_Bus[dest]])**2)
                <= Imax**2
            )

            # limite pour Idest
            model.BranchCurrentBalance.add(
                ((Ybranch[k-1, 1, 0].real * model.Vre[ind_Bus[orig]] - Ybranch[k-1, 1, 0].imag * model.Vim[ind_Bus[orig]]
                + Ybranch[k-1, 1, 1].real * model.Vre[ind_Bus[dest]] - Ybranch[k-1, 1, 1].imag * model.Vim[ind_Bus[dest]])**2
                +
                (Ybranch[k-1, 1, 0].imag * model.Vre[ind_Bus[orig]] + Ybranch[k-1, 1, 0].real * model.Vim[ind_Bus[orig]]
                + Ybranch[k-1, 1, 1].imag * model.Vre[ind_Bus[dest]] + Ybranch[k-1, 1, 1].real * model.Vim[ind_Bus[dest]])**2)
                <=Imax**2
            )
     
    # Définition de la fonction objective  
    model.obj = Objective(expr=model.val_obj, sense=minimize)

    return model

def test_model_matpower(filepath):

    #récupération du problème OPF de matpower
    start = time.time()
    case = parse_matpower_file(filepath)
    end = time.time()
    print(" Le temps utilisé pour lire le fichier et crée le dictionnaire contenant les données est de : ", end-start)
    tot_t = end-start
    # Définition du model lié au problème récupéré
    start = time.time()
    model = GenModelOPFR(case)
    end = time.time()
    print(" le temps utilisé pour creer le model est de : ",end-start)
    tot_t += (end-start)

    # Affichage du model pour vérification
    #print("\n\n\n\n\n\n\n\n\nle modèle crée est le suivant :")
    #model.pprint()

    # Création de l'objet solveur
    solver = SolverFactory('ipopt', executable='/Users/djelelbenabdellah/anaconda3/envs/opf-env/bin/ipopt')
    #solver = SolverFactory('mosek')
    # Définir le nombre maximum d'itération de Ipopt
    #solver.options['max_iter'] = 100000
    

    # Résolution du problème d'optimisation sans affichage
    start = time.time()
    results = solver.solve(model)
    end = time.time()
    print(" le temps utilisé par le solveur pour résoudre le model est de : ", end-start)
    tot_t += (end-start)
    print("le temps utilisé pour réaliser ces trois etapes est de : ", tot_t)
    # Affichage des valeurs des variables
    #for v in model.component_data_objects(Var, active=True):
    #    print(f"{v} = {value(v)}")

    # Affichage de la fonction objectif
    print(f"Fonction objective = {value(model.obj)}")

    #ecriture du model dans le fichier 
    # Ouvrez un fichier texte en mode écriture
    with open('modele_pprint.txt', 'w') as f:
        # Redirigez la sortie du modèle vers le fichier
        model.pprint(ostream=f)

    ## Vérifiez que le solveur a convergé
    #print(results.solver.status)
    #print(results.solver.termination_condition)

def test_model_GOC(filepath):

    #récupération du problème OPF de matpower
    case = parse_GOC_file(filepath)

    # Définition du model lié au problème récupéré
    model = GenModelOPFR(case)

    # Affichage du model pour vérification
    #print("\n\n\n\n\n\n\n\n\nle modèle crée est le suivant :")
    #model.pprint()

    # Création de l'objet solveur
    solver = SolverFactory('ipopt', executable='/Users/djelelbenabdellah/anaconda3/envs/opf-env/bin/ipopt')

    # Définir le nombre maximum d'itération de Ipopt
    solver.options['max_iter'] = 100000

    # Résolution du problème d'optimisation sans affichage
    results = solver.solve(model)

    # Affichage des valeurs des variables
    #for v in model.component_data_objects(Var, active=True):
    #    print(f"{v} = {value(v)}")

    # Affichage de la fonction objectif
    print(f"Fonction objective = {value(model.obj)}")

    ## Vérifiez que le solveur a convergé
    #print(results.solver.status)
    #print(results.solver.termination_condition)



# Tests

# Tests Matpower
print("-----------------------------------------------------------------------------------------------------------------------")
print("\ncase 9 :")
test_model_matpower('./data/matpower/case9.m')
print("\ncase 14 :")
test_model_matpower('./data/matpower/case14.m')
print("\ncase 30 :")
test_model_matpower('./datat/matpower/case30.m')
print("\ncase 39 :")
test_model_matpower('./data/matpower/case39.m')
print("\ncase 57 :")
test_model_matpower('./data/matpower/case57.m')
print("\ncase 145 :")
test_model_matpower('./data/matpower/case145.m')
print("-----------------------------------------------------------------------------------------------------------------------")


#Tests GOC
#print("test goc sans cliques:")
#test_model_GOC('./data/GOC/scenario_009/case.raw')




