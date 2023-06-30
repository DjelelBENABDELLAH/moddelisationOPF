import numpy as np
import networkx as nx
from RecupDataMat import parse_matpower_file
from RecupDataGOC import parse_GOC_file
import random

# fonction qui permet d'obtenir la matrice d'adjacence de l'OPF à partir de la matrice bus et branch
def create_adjacency_matrix(branch, bus):
    N = len(bus) # Trouver le numéro de nœud maximum
    adjacency_matrix = np.zeros((N, N))  # Créer une matrice de zéros avec les dimensions appropriées
    bus = np.array(bus)
    branch = np.array(branch)
    iBus = {}
    for i in range(N):
        iBus[int(bus[i,0])]= i
    for b in branch:
        node_from = int(iBus[b[0]])  # Noeud d'origine
        node_to = int(iBus[b[1]])  # Noeud de destination

        # Définition d'une connexion entre les nœuds d'origine et de destination
        adjacency_matrix[node_to, node_from] = 1
        adjacency_matrix[node_from, node_to] = 1 

    return adjacency_matrix

def chordal_completion(adjacency_matrix):
    chordal_matrix = adjacency_matrix
    # Convert the adjacency matrix to a NetworkX graph
    graph = nx.Graph(adjacency_matrix)

    # Perform chordal completion
    i = random.randint(0,3)
    print(i)
    if i==0:
        chordal_graph, mapping = nx.complete_to_chordal_graph(graph)
        print(chordal_graph)

    if i==1:
        chordal_graph = nx.chordal_graph_cliques(graph)
        print(chordal_graph)

    if i == 2:
        chordal_graph = nx.chordal_graph_treewidth(graph)
        print(chordal_graph)

    if i == 3:
        chordal_graph = nx.chordal_cycle_graph
        print(chordal_graph)


    for (i,j) in chordal_graph.edges():
        if chordal_matrix[i,j]==0 :
            chordal_matrix[i,j]= 1
            chordal_matrix[j,i]= 1
    return chordal_matrix

def find_cliques_from_adjacency_matrix(adjacency_matrix):
    # Créez un graphe à partir de la matrice d'adjacence
    graph = nx.Graph(adjacency_matrix)

    # Recherche de toutes les cliques dans le graphe
    cliques = list(nx.find_cliques(graph))

    return cliques

def find_clique_tree_from_adacency_matrix(adjacency_matrix):
    chordal_matrix = chordal_completion(adjacency_matrix)
    #print('chordal matrix.')
    cliques = find_cliques_from_adjacency_matrix(chordal_matrix)
    #print('cliques')
    clique_tree = construct_max_clique_tree(chordal_matrix) 
    #print('clique tree')
    return chordal_matrix, cliques, clique_tree.edges

def find_clique_tree_from_caseMatpower(filename):
    opf = parse_matpower_file(filename)
    adjacency_matrix = create_adjacency_matrix(opf['branch'],opf['bus'])
    #print('adjacency matrix')
    return find_clique_tree_from_adacency_matrix(adjacency_matrix)

def find_clique_tree_from_caseGOC(filename):
    opf = parse_GOC_file(filename)
    adjacency_matrix = create_adjacency_matrix(opf['branch'],opf['bus'])
    #print('adjacency matrix')
    return find_clique_tree_from_adacency_matrix(adjacency_matrix)


def compute_similarity(clique1, clique2):
    return len(set(clique1) & set(clique2))

def construct_max_clique_tree(chordal_matrix):

    # Créer un graphe non orienté
    graph = nx.Graph()

    # rechercher les cliques du graph
    cliques = find_cliques_from_adjacency_matrix(chordal_matrix)

    # Ajouter les cliques en tant que nœuds du graphe
    for clique in cliques:
        graph.add_node(tuple(clique))

    # Ajouter les arêtes entre les cliques basées sur une mesure de similarité ou de poids
    for i in range(len(cliques)):
        for j in range(i + 1, len(cliques)):
            similarity = compute_similarity(cliques[i], cliques[j])  # Remplacez compute_similarity par votre mesure de similarité ou de poids
            graph.add_edge(tuple(cliques[i]), tuple(cliques[j]), weight=similarity)

    # Utiliser l'algorithme de Prim pour trouver l'arbre de poids couvrant maximal
    clique_tree = nx.maximum_spanning_tree(graph)


    return clique_tree



#test chordal completion

#adj_matrix = np.array([[0, 1, 1, 0], [1, 0, 0, 1], [1, 0, 0, 1], [0,1,1,0]]) # Création d'un cycle de taille 4
#print(chordal_completion(adj_matrix))

# test intermédiaire

opf = parse_matpower_file('cases/case9.m')
adj_matrix = create_adjacency_matrix(opf['branch'],opf['bus'])


#print(adj_matrix)
#print('-------------------------------------------')
print(chordal_completion(adj_matrix))

#ch_mat, cliques , clique_tree_edges = find_clique_tree_from_adacency_matrix(adj_matrix)

#print('-----------------------------------')
#print('ch_mat : ')
#print(ch_mat)
#print('-----------------------------------')
#print('cliques : ')
#print(cliques)
#print('-----------------------------------')
#print('clique_tree_edges :')
#print(clique_tree_edges)
#print('-----------------------------------')

# Tests matpower

# ch_mat, cliques, clique_tree_edges = find_clique_tree_from_caseMatpower('cases/case9.m')
# print('-----------------------------------')
# print('ch_mat : ')
# print(ch_mat)
# print('-----------------------------------')
# print('cliques : ')
# print(cliques)
# print('-----------------------------------')
# print('clique_tree_edges :')
# print(clique_tree_edges)
# print('-----------------------------------')

#Tests GOC

# ch_mat, cliques, clique_tree_edges = find_clique_tree_from_caseGOC('scenario/case.raw')
# print('-----------------------------------')
# print('ch_mat : ')
# print(ch_mat)
# print('-----------------------------------')
# print('cliques : ')
# print(cliques)
# print('-----------------------------------')
# print('clique_tree_edges :')
# print(clique_tree_edges)
# print('-----------------------------------')
