from rdkit import Chem
from get_functions import _get_res, _get_res_rmsds
from models import Object, Owner

def _parse_ligand_sdf(input_file):
    """
    Function to parse a series of ligands - return RDKit mols.
    :param input_file: the file to parse
    :param input_type: the type of ligands
    :return: the molecules parsed
    """
    mols = Chem.SDMolSupplier(input_file)
    return mols


def _parse_mols(input_file, input_format):
    if input_format == "smi":
        return Chem.SmilesMolSupplier(input_file, delimiter=",")
    else:
        return Chem.SDMolSupplier(input_file)


def _parse_pdb(data):
    return Chem.MolFromPDBFile(data)


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
        if 'HOH' not in res and len(res_dict[res]) == 15 and '_A_' in res:
            rmsd_coords = _get_res_rmsds(res_dict[res])
            out_l = []
            res = Object(rmsd_coords, res)
            out_l.append(res)
            owner = Owner(out_l, input_pdb)
            owner_list.append(owner)
    return owner_list
