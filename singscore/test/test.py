import unittest, pandas,os
from singscore.singscore import score

class SingscoreTestCase(unittest.TestCase):





    def test_score_entrez_up_down_path(self):
        """
        test score() on entrez identifiers for sig and samples

        :return: .
        """
        file_path = (os.path.dirname(__file__)) + '/all_data/entrez_sample.txt'
        f = open(file_path, 'r')
        up_path = (os.path.dirname(__file__)) +'/all_sigs/entrez_up.txt'
        down_path = (os.path.dirname(__file__)) + '/all_sigs/entrez_down.txt'

        sample = pandas.read_csv(f,header = 'infer', sep='\t')
        sample = sample.set_index(keys='GeneID')
        f.close()
        
        self.assertIsInstance(score(up_gene=up_path, down_gene=down_path,
                                    sample= sample), pandas.DataFrame)


if __name__ == '__main__':
    unittest.main()