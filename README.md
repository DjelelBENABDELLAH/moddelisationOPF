
# Projet modelisation OPF (Optimal Power Flow):

    ## Description:

    Le projet modelOPF est une implémentation du problème de l'Optimal Power Flow (OPF) en utilisant le langage de programmation Python et la bibliothèque Pyomo.

    L'objectif de ce projet est de modéliser et résoudre le problème OPF, qui consiste à optimiser la production et la distribution d'électricité dans un réseau électrique tout en respectant des contraintes de capacité, de flux et de stabilité.
    En passant par une formulation semi-definite positive(SDP) qui a pour objectif d'être ensuite redéfini en plusieurs blocs selon les cliques de la version chordale du graph représentant l'OPF.

    ## Fonctionnalités :

    - Modélisation du problème OPF en valeur réélle à l'aide de Pyomo
    - Résolution du modèle OPF en utilisant des solveurs d'optimisation ("ipopt")
    - Séparation du modèle OPF en utilisant des cliques générées à partir de la version chordale du graphe obtenue à l'aide de la fonction fournit par networkX
    - Modélisation et résolution du modèle OPF séparé en utilisant les cliques et le clique-tree obtenue

    ## librairie requise :

    - Python
    - Pyomo
    - Solveur d'optimisation (dans le code actuel Ipopt mais peut être modifié)
    - NetworkX
    - numpy

    ## Installation :
    .....

    ## Utilisation : 
    
    Il existe dans le code des fonctions tests permettant l'utilisation du programme. Qui peuvent ensuite être adapté selon les preferences de l'utilisateur : 
      - ajout d'un fichier de sortie qui représente le model
      - affichage du model ou des valeurs objectives dans le terminal 


    


