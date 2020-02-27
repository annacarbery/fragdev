import rdkit
from rdkit import Chem
import os
from rdkit.Chem import ChemicalFeatures
from get_functions import _get_c_of_mass

def parse_ligand_ph4s(input_mols):
    """
    Function to return a series of ligand based pharmacophores.
    :param input_mols: the RDKit molecules
    :return: the molecule based pharmacophores
    """
    rdkit_ph4 = RDKitPh4()
    rdkit_atom = RDKitAtom()
    output_pharma_list = []
    for mol in input_mols:
        if not mol:
            pharma_list = []
        else:
            pharma_list = rdkit_ph4.generate_ph4_for_mol(rdmol=mol)
            atom_list = rdkit_atom.generate_atoms_for_mol(mol)
            x, y, z = _get_c_of_mass(rdmol=mol)
            c_of_m_feat = (x, y, z, "c_of_m")
            pharma_list.append(c_of_m_feat)
            pharma_list.extend(atom_list)
            print(mol, 'added to output ph4 list')
        output_pharma_list.append(pharma_list)
    return output_pharma_list


class RDKitPh4(object):

    factory = None

    def __init__(self):
        # Generate the factory on init
        self.get_factory()

    def get_factory(self):
        """
        Generate the Ph4 feature factory
        :return:
        """
        if self.factory is None:
            this_dir, this_filename = os.path.split(__file__)
            data_path = os.path.join(this_dir, "RDKitPh4.fdef")
            self.factory = ChemicalFeatures.BuildFeatureFactory(data_path)
        return self.factory

    def generate_ph4_for_mol(self, rdmol):
        """
        Generate a pharmacophore from an input molecule and a feature factory.
        :param rdmol: the input RDKit molecule
        :param factory: the feature factory
        :return: a list of 4 tuples (x,y,z, feature)
        """
        feats = self.get_factory().GetFeaturesForMol(rdmol)
        return [
            (feat.GetPos().x, feat.GetPos().y, feat.GetPos().z, feat.GetType())
            for feat in feats
        ]


class RDKitAtom(object):

    def generate_atoms_for_mol(self, rdmol):
        """
        Generate the atoms from an input molecule and a feature factory.
        :param rdmol: the input RDKit molecule
        :return: a list of 4 tuples (x,y,z, atom description)
        """
        out_list = []
        conf = rdmol.GetConformer()
        for atom in rdmol.GetAtoms():
            atom_desc = self.get_atom_description(atom)
            atom_pos = conf.GetAtomPosition(atom.GetIdx())
            out_list.append((atom_pos.x, atom_pos.y, atom_pos.z, atom_desc))
        return out_list

    def get_atom_description(self, atom):
        """
        Generate a unique description of an atom
        :param atom: the input atom
        :return: a hash string of atomic number and  hybridization state
        """
        return "_".join(
            [str(x) for x in [atom.GetAtomicNum(), atom.GetHybridization()]]
        )

