import os
import math
from rdkit.Chem import ChemicalFeatures
from rdkit import Chem
from frag.alysis.cluster import dp_means


def _get_water_coords(waters):
    """Helper function to get the coordinates from a load of waters."""
    rd_waters = Chem.MolFromPDBBlock("\n".join(waters))
    out_list = []
    if rd_waters is None:
        print("Warning - unable to parse waters.")
    if rd_waters is not None:
        # Check the waters exist
        conf = rd_waters.GetConformer()
        # Delete them for this protein
        for i in range(rd_waters.GetNumAtoms()):
            cp = conf.GetAtomPosition(i)
            if rd_waters.GetAtomWithIdx(i).GetSmarts() != "O":
                print("Warning - skipping a water")
                continue
            out_list.append((cp.x, cp.y, cp.z))
    return out_list

def map_cluster(dp_means_cluster, mol_id_list):
    """

    :param dp_means_cluster:
    :param mol_id_list:
    :return:
    """
    out_dict = {}
    for cluster in dp_means_cluster.clusters:
        out_dict[cluster] = {
            "centre_of_mass": dp_means_cluster.clusters[cluster],
            "mol_ids": [],
        }
    for i, cluster_id in enumerate(dp_means_cluster.dataClusterId):
        out_dict[cluster_id]["mol_ids"].append(mol_id_list[i])
    return out_dict

def cluster_dp(vect_list, lam, mol_list):
    """
    Perform a DP Means clustering.
    :param vect_list: a list of lists of vectors
    :param lam: the clustering parameters
    :param mol_list: the molecular identifers in the same order as the list of vectors
    :return: a dictionary of the form {cluster_id: {centre_of_mass: (x,y,z), mol_ids: [1,5,12]}}
    """
    return map_cluster(dp_means(vect_list, lam), mol_list)