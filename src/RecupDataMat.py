
import numpy as np

# fonction qui permet de récupérer les données contenue dans un fichier caseN.m de matpower dans un dictionnaire
def parse_matpower_file(filename):

    with open(filename, 'r') as file:
        content = file.read()

    # Initialize the data dictionary
    data = {}

    # Récupération de la donnée BaseMVA 
    baseMVA_start = content.index("mpc.baseMVA = ") + len("mpc.baseMVA = ")
    baseMVA_end = content.index(";", baseMVA_start)
    baseMVA = float(content[baseMVA_start:baseMVA_end])
    data['baseMVA'] = baseMVA

    # Récupération de la matrice gen 
    gen_start = content.index("mpc.gen = ") + len("mpc.gen = ")
    gen_end = content.index("];", gen_start) + 1
    gen_data = content[gen_start:gen_end]
    gen_lines = gen_data.splitlines()
    gen_lines = [line.replace(';', '') for line in gen_lines if not line.startswith('[') and not line.startswith(']') and not line.startswith(' ]')  and line.strip()]
    data['gen'] = np.array([list(map(float, line.split())) for line in gen_lines])

    # Récupértaion de la matrice gencost
    gencost_start = content.index("mpc.gencost = ") + len("mpc.gencost = ")
    gencost_end = content.index("];", gencost_start) + 1
    gencost_data = content[gencost_start:gencost_end]
    gencost_lines = gencost_data.splitlines()
    gencost_lines = [line.replace(';', '') for line in gencost_lines if not line.startswith('[') and not line.startswith(']') and not line.startswith(' ]') and line.strip()]
    data['gencost'] = np.array([list(map(float, line.split())) for line in gencost_lines])

    # Récupération de la matrice branch
    branch_start = content.index("mpc.branch = ") + len("mpc.branch = ")
    branch_end = content.index("];", branch_start) + 1
    branch_data = content[branch_start:branch_end]
    branch_lines = branch_data.splitlines()
    branch_lines = [line.replace(';', '') for line in branch_lines if not line.startswith('[') and not line.startswith(']') and not line.startswith(' ]') and line.strip()]
    data['branch'] = np.array([list(map(float, line.split())) for line in branch_lines])

    # Récupération de la matrice bus
    bus_start = content.index("mpc.bus = ") + len("mpc.bus = ")
    bus_end = content.index("];", bus_start) + 1
    bus_data = content[bus_start:bus_end]
    bus_lines = bus_data.splitlines()
    bus_lines = [line.replace(';', '') for line in bus_lines if not line.startswith('[') and not line.startswith(']')  and not line.startswith(' ]') and line.strip()]
    data['bus'] = np.array([list(map(float, line.split())) for line in bus_lines])

    return data

# Test

#filename = './data/matpower/case9.m'
#data = parse_matpower_file(filename)
#print(data['baseMVA'])
#print(data['gen'])
#print(data['gen'][0])

