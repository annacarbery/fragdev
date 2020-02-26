import os
from frag.utils.rdkit_utils import _get_water_coords
from rdkit import Chem
from frag.alysis.run_clustering import cluster_dp


DATA_DIRECTORY = os.path.abspath('data')


for dir in os.listdir(DATA_DIRECTORY):
    out_list = []
    id_list = []
    for file in os.listdir(os.path.join(DATA_DIRECTORY, dir)):
        pdb = open(os.path.join(DATA_DIRECTORY, dir, file)).readlines()
        waters = [x for x in pdb if x[17:20] == "HOH"]
        rd_waters = [Chem.MolFromPDBBlock(i) for i in waters]
        coords = _get_water_coords(waters)
        for i in rd_waters:
            conf = i.GetConformer()
            for j in range(i.GetNumAtoms()):
                cp = conf.GetAtomPosition(j)
                out_list.append((cp.x, cp.y, cp.z))
                id_list.append(file)
    print(len(out_list))
    clusters = cluster_dp(out_list, 1.5, id_list)

for i in clusters:
    print(i, len(clusters[i]['mol_ids']), clusters[i])

