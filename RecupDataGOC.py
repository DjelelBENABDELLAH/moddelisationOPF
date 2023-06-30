import json
import numpy as np

# Ajouter param JSON lorsque la facon de retrouver les valeurs C1 et C0 aura été trouver
def parse_GOC_file(filenameRAW):
    with open(filenameRAW, 'r') as f:
        ppc = {}
        first_line = f.readline().strip()
        elements = first_line.split(',')
        ppc['baseMVA'] = elements[1]
    # Ignorer les deux premières lignes
        f.readline()
        f.readline()

    # Récupération des donées contenues dans case.raw en utilisant les information contenue dans le pdf 'Challenge2_Problem_Formulation_20210531.pdf' et une visualisation du fichier
        # Récupération des données relatives au bus
        BusData = []
        for line in f:
            line = line.strip()
            if line == '0 / END OF BUS DATA BEGIN LOAD DATA':
                break
            elements = line.split(',')
            elements[0] = int(elements[0])  # Deuxième élément
            elements[2] = float(elements[2].strip("'"))  # Troisième élément
            BusData.append(elements)
        
        # Récupération des données de Puisssance demandé en chaque noeud
        LoadData = []
        for line in f:
            line = line.strip()
            if line == '0 / END OF LOAD DATA BEGIN FIXED SHUNT DATA':
                break
            elements = line.split(',')
            elements[0] = int(elements[0])
            elements[5] = float(elements[5])
            elements[6] = float(elements[6])
            LoadData.append(elements)

        FixedShuntData = []
        for line in f:
            line = line.strip()
            if line == '0 / END OF FIXED SHUNT DATA BEGIN GENERATOR DATA':
                break
            elements = line.split(',')
            elements[0]=int(elements[0])
            FixedShuntData.append(elements)
        
        GenData = []
        for line in f:
            line = line.strip()
            if line == '0 / END OF GENERATOR DATA BEGIN BRANCH DATA':
                break
            elements = line.split(',')
            GenData.append(elements)
        
        BranchData = []
        for line in f:
            line = line.strip()
            if line == '0 / END OF BRANCH DATA BEGIN TRANSFORMER DATA':
                break
            elements = line.split(',')
            BranchData.append(elements)
        
        TransformerData = []
        for line in f:
            line = line.strip()
            if line == '0 / END OF TRANSFORMER DATA BEGIN AREA DATA':
                break
            elements = line.split(',')
            TransformerData.append(elements)
        
        AreaData = []
        for line in f:
            line = line.strip()
            if line == '0 / END OF AREA DATA BEGIN TWO-TERMINAL DC DATA':
                break
            elements = line.split(',')
            AreaData.append(elements)
        
        TwoTerminalDCData = []
        for line in f:
            line = line.strip()
            if line == '0 / END OF TWO-TERMINAL DC DATA BEGIN VSC DC LINE DATA':
                break
            elements = line.split(',')
            TwoTerminalDCData.append(elements)
        
        VSCDCLineData = []
        for line in f:
            line = line.strip()
            if line == '0 / END OF VSC DC LINE DATA BEGIN IMPEDANCE CORRECTION DATA':
                break
            elements = line.split(',')
            VSCDCLineData.append(elements)
        
        ImpedanceCorrectionData = []
        for line in f:
            line = line.strip()
            if line == '0 / END OF IMPEDANCE CORRECTION DATA BEGIN MULTI-TERMINAL DC DATA':
                break
            elements = line.split(',')
            ImpedanceCorrectionData.append(elements)

        MultiTerminalDCData = []
        for line in f:
            line = line.strip()
            if line == '0 / END OF MULTI-TERMINAL DC DATA BEGIN MULTI-SECTION LINE DATA':
                break
            elements = line.split(',')
            MultiTerminalDCData.append(elements)
        
        MultiSectionLineData = []
        for line in f:
            line = line.strip()
            if line == '0 / END OF MULTI-SECTION LINE DATA BEGIN ZONE DATA':
                break
            elements = line.split(',')
            MultiSectionLineData.append(elements)

        ZoneData = []
        for line in f:
            line = line.strip()
            if line == '0 / END OF ZONE DATA BEGIN INTER-AREA TRANSFER DATA':
                break
            elements = line.split(',')
            ZoneData.append(elements)

        InterAreaTransferData = []
        for line in f:
            line = line.strip()
            if line == '0 / END OF INTER-AREA TRANSFER DATA BEGIN OWNER DATA':
                break
            elements = line.split(',')
            InterAreaTransferData.append(elements)
        
        OwnerData = []
        for line in f:
            line = line.strip()
            if line == '0 / END OF OWNER DATA BEGIN FACTS DEVICE DATA':
                break
            elements = line.split(',')
            OwnerData.append(elements)
        
        FactsDeviceData = []
        for line in f:
            line = line.strip()
            if line == '0 / END OF FACTS DEVICE DATA BEGIN SWITCHED SHUNT DATA':
                break
            elements = line.split(',')
            FactsDeviceData.append(elements)
        
        SwitchedShuntData = []
        for line in f:
            line = line.strip()
            if line == '0 / END OF SWITCHED SHUNT DATA BEGIN GNE DATA':
                break
            elements = line.split(',')
            SwitchedShuntData.append(elements)

        GneData = []
        for line in f:
            line = line.strip()
            if line == '0 / END OF GNE DATA BEGIN INDUCTION MACHINE DATA':
                break
            elements = line.split(',')
            GneData.append(elements)
        
        InductionMachineData = []
        for line in f:
            line = line.strip()
            if line == '0 / END OF INDUCTION MACHINE DATA':
                break
            elements = line.split(',')
            InductionMachineData.append(elements)

    # Tentatives de récupération des couts à partire de case.json mais données différentes de matpower nécessite une redéfinissions
    #with open(filenameJSON, "r") as f:
    #    donnees = json.load(f)
    #    genCostData = donnees['generators']
       
    # Création de la matrice contenant les données relatives aux bus
    Bus = []
    iBus = {}
    i = 0
    for row in BusData:
        Bus.append([row[0],0,0,0,0,0,int(row[4]),float(row[7]),float(row[8]),row[2],int(row[5]),float(row[9]),float(row[10])])
        iBus[row[0]] = i
        i += 1
    for row in LoadData:
        #print(row)
        Bus[iBus[row[0]]][2] = float(row[5])
        Bus[iBus[row[0]]][3] = float(row[6])
    for row in FixedShuntData:
        #print(row)
        Bus[iBus[row[0]]][4] = float(row[3])
        Bus[iBus[row[0]]][5] = float(row[4])

    # Création de des matrice contenant les données relatives aux générateurs et à leurs couts(à redéfinire)
    Gen = []
    GenCost = []
    for row in GenData:
        Gen.append([int(row[0]),float(row[2]),float(row[3]),float(row[4]),float(row[5]), 0, float(row[8]),int(row[14]), float(row[16]),float(row[17]), 0,0,0,0,0, 0,0,0,0,0, 0])
        GenCost.append([2, 0, 0, 3, 0, 1, 0])

    # Création de la première partie de la matrice branche contenant les données des branches qui ne sont pas des transformateurs
    Branch = []
    for row in BranchData:
        Branch.append([int(row[0]),int(row[1]),float(row[3]),float(row[4]), float(row[5]), float(row[6]), float(row[7]), float(row[8]), 0, 0, int(row[13]), 0, 0])

    # Ajout de la seconde partie de la matrice branche contenant les données des transformateurs
    # les données des transformateurs sont définis par bloc de 4 lignes
    for i in range(0, len(TransformerData), 4):
        row1 = TransformerData[i]
        row2 = TransformerData[i + 1]
        row3 = TransformerData[i + 2]
        row4 = TransformerData[i + 3]
        Branch.append([int(row1[0]),int(row1[1]),float(row2[0]),float(row2[1]),0,0,0,0, float(row3[8]), float(row3[2]), int(row1[11]), -360, 360])



    # Création du dictionnaire similaire à celui récupéré dans un fichier matpower
    ppc['bus'] = Bus
    ppc['gen'] = Gen
    ppc['branch'] = Branch
    ppc['gencost'] = GenCost

    return ppc

#test 
# ppc = parse_GOC_file('scenario/scenario_043/case.raw')
# bus = np.array(ppc['bus'])
# print(sum(bus[:,2]))
# print(sum(bus[:,3]))
# print(ppc['baseMVA'])
# print(ppc['gen'][0])
# print('-----------------------------------------------------------')
# print(ppc['gencost'][0])
# print('-----------------------------------------------------------')
# print(ppc['branch'][0])
# print('-----------------------------------------------------------')
# print(ppc['bus'][0])
# print('-----------------------------------------------------------')
# print(len(ppc['gen']))
# print(len(ppc['gencost']))
# print(len(ppc['branch']))
# print(len(ppc['bus']))
