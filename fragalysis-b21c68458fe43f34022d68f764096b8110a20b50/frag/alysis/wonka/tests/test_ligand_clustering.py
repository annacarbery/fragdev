import unittest
import os
from frag.alysis.wonka.ligand_clustering import cluster
import json

DATA_DIRECTORY = os.path.abspath('data')
RESULTS_DIRECTORY = os.path.abspath('results')

class MyTestCase(unittest.TestCase):
    def test_something(self):
        cluster(DATA_DIRECTORY, RESULTS_DIRECTORY)
        ligand_clusters = json.load(
            open(os.path.join(RESULTS_DIRECTORY, 'ATAD_allPdb_28-Feb-2020', 'ligand_cluster.json'), 'r'))
        self.assertEqual(len(ligand_clusters['SingleAtomDonor']['0']['mol_ids']), 8)


if __name__ == '__main__':
    unittest.main()
