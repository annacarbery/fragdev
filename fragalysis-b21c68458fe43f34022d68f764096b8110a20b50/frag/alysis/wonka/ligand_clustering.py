import os
from cluster_functions import cluster_dp
from rdkit_functions import parse_ligand_ph4s
import json
from rdkit import Chem

PH4_LAMBDA = 2.0
C_OF_M_LAMBDA = 6.0


def _parse_pdb(data):
    return Chem.MolFromPDBFile(data)


def build_type_dict(mol_ph4_list, identifiers):
    type_dict = {}
    for i, mol in enumerate(mol_ph4_list):
        for ph4 in mol:
            x = ph4[0]
            y = ph4[1]
            z = ph4[2]
            ph4_type = ph4[3]
            if ph4_type in type_dict:
                type_dict[ph4_type]["coords"].append((x, y, z))
                type_dict[ph4_type]["mols"].append(identifiers[i])
            else:
                type_dict[ph4_type] = {"coords": [(x, y, z)], "mols": [identifiers[i]]}
    return type_dict


def run_lig_cluster(mols, identifiers):
    """
    Cluster a list of molecules
    :param mols: the input list of RDKit molecules
    :param identifiers: the corresponding list of identifiers
    :return:
    """
    # First we get the list of mols with their Ph4s
    print('parsing ligand ph4s')
    mol_ph4_list = parse_ligand_ph4s(mols)
    print('created ph4s')
    # Then we build a dict of type: coords: [coords list], mols: [mol_index]
    type_dict = build_type_dict(mol_ph4_list, identifiers)
    print('built type dictionary')
    # Then we cluster coords
    clusters = {}
    for ph4_type in type_dict:
        if ph4_type == "c_of_m":
            clusters[ph4_type] = cluster_dp(
                type_dict[ph4_type]["coords"],
                C_OF_M_LAMBDA,
                type_dict[ph4_type]["mols"],
            )
        else:
            clusters[ph4_type] = cluster_dp(
                type_dict[ph4_type]["coords"], PH4_LAMBDA, type_dict[ph4_type]["mols"]
            )
    print('clustered')
    return clusters


def cluster(DATA_DIRECTORY, RESULTS_DIRECTORY):
 #   for dir in os.listdir(DATA_DIRECTORY):
    dir = 'Mpro_allPdb_06-Mar-2020'
    if dir not in os.listdir(RESULTS_DIRECTORY):
        os.makedirs(os.path.join(RESULTS_DIRECTORY, dir))
    print('dir: ', dir)
    mols = []
    identifiers = []
    for file in os.listdir(os.path.join(DATA_DIRECTORY, dir)):
        print('file: ', file)
        mols.append(_parse_pdb(os.path.join(DATA_DIRECTORY, dir, file)))
        identifiers.append(str(file))
    print('out of loop')
    output = run_lig_cluster(mols, identifiers)
    print('wonka: ', output)
    json.dump(output, open(os.path.join(RESULTS_DIRECTORY, dir, 'ligand_cluster.json'), 'w'))


DATA_DIRECTORY = os.path.abspath('data')
RESULTS_DIRECTORY = os.path.abspath('results')

cluster(DATA_DIRECTORY, RESULTS_DIRECTORY)
