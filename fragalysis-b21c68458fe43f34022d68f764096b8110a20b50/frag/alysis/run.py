from frag.alysis.run_clustering import run_lig_cluster
from frag.utils.rdkit_utils import _parse_pdb
import os

DATA_DIRECTORY = os.path.abspath('data')
print(DATA_DIRECTORY)

for dir in os.listdir(DATA_DIRECTORY):
    mols = []
    identifiers = []
    for file in os.listdir(os.path.join(DATA_DIRECTORY, dir)):
        mols.append(_parse_pdb(os.path.join(DATA_DIRECTORY, dir, file)))
        identifiers.append(str(file))
    print(run_lig_cluster(mols, identifiers))
