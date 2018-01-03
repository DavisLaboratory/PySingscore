import unittest, pandas
from singscore.singscore import score

class SingscoreTestCase(unittest.TestCase):

    def test_score_entrez_up_down_path(self):
        """
        test score() on entrez identifiers for sig and samples

        :return: .
        """

        up_path = 'all_sigs/entrez_up.txt'
        down_path = 'all_sigs/entrez_down.txt'

        sample = pandas.read_csv(open('all_data/validation_merge.txt', 'r'),
                                 header = 'infer', sep='\t')
        sample = sample.set_index(keys='GeneID')


        self.assertIsInstance(score(up_gene=up_path, down_gene=down_path,
                                    sample= sample))