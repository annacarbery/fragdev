import os
from frag.utils.rdkit_utils import _parse_pdb
from rdkit import Chem
from frag.alysis.run_clustering import cluster_dp
#from frag.utils.parser import parse_residues
from frag.alysis.models import Object, Owner
from frag.utils.rdkit_utils import (
    _parse_ligand_sdf,
    _get_c_of_mass,
    RDKitPh4,
    RDKitAtom,
    _get_water_coords,
    _get_waters,
    _get_res,
    _get_res_rmsds,
    _parse_pdb,
)


def parse_residues(input_pdbs, input_mol=None, max_dist=10.0):
    """
    Function to parse a series of PDB files of proteins.
    :param input_pdb: the input PDB files - with identifiers
    :return: a dict (Key Residue -> value list of molecules)
    """
    owner_list = []
    res_dict = {}
    for input_pdb in input_pdbs:
        # Loop through the residues
        mol = _parse_pdb(input_pdb)
        this_res_dict = _get_res(mol)
        for key in this_res_dict:
            if key in res_dict:
                res_dict[key].append(this_res_dict[key])
            else:
                res_dict[key] = [this_res_dict[key]]
    for res in res_dict:
        rmsd_coords = _get_res_rmsds(res_dict[res])
        out_l = []
        res = Object(rmsd_coords, res)
        out_l.append(res)
        owner = Owner(out_l, input_pdb)
        owner_list.append(owner)
    return owner_list


DATA_DIRECTORY = os.path.abspath('data')


for dir in os.listdir(DATA_DIRECTORY):
    pdb_files = []
    for file in os.listdir(os.path.join(DATA_DIRECTORY, dir)):
        pdb = os.path.join(DATA_DIRECTORY, dir, file)
        pdb_files.append(pdb)
    residues = parse_residues(pdb_files)
    for i in residues:
        print(i.object_list)
