from rdkit import Chem
import math

def _get_c_of_mass(rdmol):
    """
    Get the unweighted centre of mass of an RDKit Molecule
    :param rdmol:
    :return:
    """
    atoms = rdmol.GetAtoms()
    conf = rdmol.GetConformer()
    x_coord = y_coord = z_coord = 0.0
    numatoms = 0.0
    # Assume all heavy atoms have the same mass
    for atom in atoms:
        if atom.GetAtomicNum() == 1 or atom.GetSmarts() == "[*]":
            continue
        numatoms += 1.0
        coords = conf.GetAtomPosition(atom.GetIdx())

        x_coord += float(coords.x)
        y_coord += float(coords.y)
        z_coord += float(coords.z)
    # Now we have all the coords -> we want to loop through
    if numatoms == 0:
        raise ValueError("No atoms in Molecules")
    return x_coord / numatoms, y_coord / numatoms, z_coord / numatoms


def _get_waters(input_lines):
    """
    Very strict but efficient way of getting waters
    :param input_lines: the lines of a PDB file we want to consider
    :return:
    """
    return [x for x in input_lines if x[17:20] == "HOH"]


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


def _get_file(file_path, output_format, file_counter):
    if output_format == "smi":
        return Chem.SmilesWriter(file_path + "_" + str(file_counter) + ".smi")
    else:
        return Chem.SDWriter(file_path + "_" + str(file_counter) + ".sdf")


def find_dist(mol_1_x, mol_1_y, mol_1_z, mol_2_x, mol_2_y, mol_2_z):
    """Function to find the square distance between two points in 3D
    Takes two len=3 tuples
    Returns a float"""
    return (
        pow((mol_1_x - mol_2_x), 2)
        + pow((mol_1_y - mol_2_y), 2)
        + pow((mol_1_z - mol_2_z), 2)
    )


def _get_mean_dist(res_one, res_two):
    tot_dist = 0.0
    num_matches = 0
    for atom in res_one:
        if atom in res_two:
            # Find the distance
            dist = find_dist(
                res_one[atom][0],
                res_one[atom][1],
                res_one[atom][2],
                res_two[atom][0],
                res_two[atom][1],
                res_two[atom][2],
            )
            tot_dist += dist
            num_matches += 1
    # Find the mean square distance
    return float(tot_dist) / float(num_matches)


def _get_res_rmsds(input_res_list):
    """
    Helper function to get the RMSDs for a list of Residues.
    :param input_res_list: The list of RDKit molecules of residues
    :return: a list of lists of RMSDS.
    """
    # Calculate RMSD from corresponding
    num_res = len(input_res_list)
    tot_res_rmsd_list = []
    for i in range(num_res):
        this_res_rmsd_list = []
        for j in range(num_res):
            if i == j:
                continue
            mean_dist = _get_mean_dist(input_res_list[i], input_res_list[j])
            # Append the root mean square distance
            root_mean_sqr = math.sqrt(mean_dist)
            this_res_rmsd_list.append(root_mean_sqr)
        tot_res_rmsd_list.append(this_res_rmsd_list)
    return tot_res_rmsd_list


def get_res_atom_name(atom, conf):
    """
    Get the information for each atom
    :param atom: the input atom
    :param conf: the molecular conformer
    :return: the unqiue residue level name, the atom name and the position
    """
    res_info = atom.GetPDBResidueInfo()
    atom_pos = conf.GetAtomPosition(atom.GetIdx())
    position = [atom_pos.x, atom_pos.y, atom_pos.z]
    identifiers = [
        res_info.GetResidueNumber(),
        res_info.GetChainId(),
        res_info.GetResidueName(),
        res_info.GetAltLoc(),
    ]
    unique_name = "_".join([str(x) for x in identifiers if x != " "])
    atom_name = res_info.GetName().strip()
    return unique_name, atom_name, position


def _get_res(mol):
    """
    Get a list of residues with RDKit mols (RDMol,RES_NAME)
    :param input_data:
    :return:
    """
    out_dict = {}
    conf = mol.GetConformer()
    atoms = mol.GetAtoms()
    for atom in atoms:
        # Get res_name
        res_name, atom_name, position = get_res_atom_name(atom, conf)
        if res_name in out_dict:
            out_dict[res_name][atom_name] = position
        else:
            out_dict[res_name] = {atom_name: position}
    return out_dict


PH4_LAMBDA = 1.0
C_OF_M_LAMBDA = 6.0
