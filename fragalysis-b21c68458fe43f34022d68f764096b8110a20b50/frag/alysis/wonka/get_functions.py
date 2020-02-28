from rdkit import Chem

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


def _get_file(file_path, output_format, file_counter):
    if output_format == "smi":
        return Chem.SmilesWriter(file_path + "_" + str(file_counter) + ".smi")
    else:
        return Chem.SDWriter(file_path + "_" + str(file_counter) + ".sdf")

