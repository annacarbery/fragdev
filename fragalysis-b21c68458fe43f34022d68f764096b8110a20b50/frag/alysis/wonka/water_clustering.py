import os
from get_functions import _get_water_coords
from rdkit import Chem
from cluster_functions import cluster_dp
import json


DATA_DIRECTORY = os.path.abspath('data')
RESULTS_DIRECTORY = os.path.abspath('results')

for dir in os.listdir(DATA_DIRECTORY):
    if dir not in os.listdir(RESULTS_DIRECTORY):
        os.makedirs(os.path.join(RESULTS_DIRECTORY, dir))
    out_list = []
    id_list = []
    for file in os.listdir(os.path.join(DATA_DIRECTORY, dir)):
        pdb = open(os.path.join(DATA_DIRECTORY, dir, file)).readlines()
        waters = [x for x in pdb if x[17:20] == "HOH"]
        rd_waters = [Chem.MolFromPDBBlock(i) for i in waters]
        print(file, rd_waters)
        #coords = _get_water_coords(waters)
        for i in rd_waters:
            conf = i.GetConformer()
            for j in range(i.GetNumAtoms()):
                cp = conf.GetAtomPosition(j)
                out_list.append((cp.x, cp.y, cp.z))
                id_list.append(file)
    print(len(out_list))
    clusters = cluster_dp(out_list, 1.5, id_list)
    json.dump(clusters, open(os.path.join(RESULTS_DIRECTORY, dir, 'water_cluster.json'), 'w'))
