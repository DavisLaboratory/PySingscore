import unittest, pandas,os, matplotlib
from singscore.singscore import score,rank, permutate, empiricalpval, \
    plotrankdist, nulldistribution, plotdispersion

class SingscoreTestCase(unittest.TestCase):


    def test_score_same_identifiers(self):
        """
        test score() when gene identifiers are the same

        :return: .
        """
        file_path = (os.path.dirname(__file__)) + \
                    '/test_data/entrez_sample.txt'
        f = open(file_path, 'r')
        sample = pandas.read_csv(f, header='infer', sep='\t')
        sample = sample.set_index(keys=sample.columns[0])
        f.close()
        # prepare signatures
        sig_path = (os.path.dirname(__file__)) + '/test_sigs/tgfb_upDown.txt'
        sig = open(sig_path, 'r')
        sigs = pandas.read_csv(sig, header=
        'infer', sep='\t')
        sig.close()
        # subset the data for up and down
        up = sigs[sigs['upDown'] == 'up']
        down = sigs[sigs['upDown'] == 'down']
        # get a list of ids
        up = list(up['EntrezID'])
        down = list(down['EntrezID'])

        self.assertIsInstance(score(up_gene=up, down_gene=down,
                                    sample= sample), pandas.DataFrame)

    def test_rank_same_identifiers(self):
        """
        test rank() when gene identifiers are the same

        :return: .
        """
        file_path = (os.path.dirname(__file__)) + \
                    '/test_data/entrez_sample.txt'
        f = open(file_path, 'r')
        sample = pandas.read_csv(f, header='infer', sep='\t')
        sample = sample.set_index(keys=sample.columns[0])
        f.close()
        # prepare signatures
        sig_path = (os.path.dirname(__file__)) +'/test_sigs/tgfb_upDown.txt'
        sig = open(sig_path, 'r')
        sigs = pandas.read_csv(sig, header=
        'infer', sep='\t')
        # subset the data for up and down
        up = sigs[sigs['upDown'] == 'up']
        down = sigs[sigs['upDown'] == 'down']
        # get a list of ids
        up = list(up['EntrezID'])
        down = list(down['EntrezID'])


        sig.close()
        self.assertIsInstance(rank(up_gene=up, down_gene=down,
                                    sample= sample), pandas.DataFrame)

    def test_permutation(self):
        """
        test permutations() using only 10 reps for speed.
        :return: .
        """
        file_path = (os.path.dirname(__file__)) + \
                    '/test_data/entrez_sample.txt'
        f = open(file_path, 'r')

        sample = pandas.read_csv(f, header='infer', sep='\t')
        f.close()
        self.assertIsInstance(permutate(sample=sample, n_up=50, n_down=50,
                                        reps=10),pandas.DataFrame)


    def test_pvalue(self):
        """
        test empirical_pval() using only 10 reps for speed.
        :return: .
        """
        file_path = (os.path.dirname(__file__)) + \
                    '/test_data/entrez_sample.txt'
        f = open(file_path, 'r')
        sample = pandas.read_csv(f, header='infer', sep='\t')
        sample = sample.set_index(keys=sample.columns[0])
        f.close()
        # prepare signatures
        sig_path = (os.path.dirname(__file__)) + '/test_sigs/tgfb_upDown.txt'
        sig = open(sig_path, 'r')
        sigs = pandas.read_csv(sig, header=
        'infer', sep='\t')
        sig.close()
        # subset the data for up and down
        up = sigs[sigs['upDown'] == 'up']
        down = sigs[sigs['upDown'] == 'down']
        # get a list of ids
        up = list(up['EntrezID'])
        down = list(down['EntrezID'])
        p = permutate(sample=sample, n_up=50, n_down=50,
                                        reps=10)
        scores = score(up_gene=up, down_gene=down,sample=sample)

        self.assertIsInstance(empiricalpval(permutations=p, score=scores),
                              pandas.DataFrame)

    def test_barcode(self):

        file_path = (os.path.dirname(__file__)) + \
                    '/test_data/entrez_sample.txt'
        f = open(file_path, 'r')
        sample = pandas.read_csv(f, header='infer', sep='\t')
        sample = sample.set_index(keys=sample.columns[0])
        f.close()
        # prepare signatures
        sig_path = (os.path.dirname(__file__)) + '/test_sigs/tgfb_upDown.txt'
        sig = open(sig_path, 'r')
        sigs = pandas.read_csv(sig, header=
        'infer', sep='\t')
        sig.close()
        # subset the data for up and down
        up = sigs[sigs['upDown'] == 'up']
        down = sigs[sigs['upDown'] == 'down']
        # get a list of ids
        up = list(up['EntrezID'])
        down = list(down['EntrezID'])
        r = rank(up_gene=up, down_gene=down,
                           sample=sample[['D_Ctrl_R1']],
                           norm_method='theoretical')

        self.assertIsInstance(plotrankdist(ranks=r, show=False),
                              matplotlib.figure.Figure)

    def test_dispersion_plot(self):

        file_path = (os.path.dirname(__file__)) + \
                    '/test_data/entrez_sample.txt'
        f = open(file_path, 'r')
        sample = pandas.read_csv(f, header='infer', sep='\t')
        sample = sample.set_index(keys=sample.columns[0])
        f.close()
        # prepare signatures
        sig_path = (os.path.dirname(__file__)) + '/test_sigs/tgfb_upDown.txt'
        sig = open(sig_path, 'r')
        sigs = pandas.read_csv(sig, header=
        'infer', sep='\t')
        sig.close()
        # subset the data for up and down
        up = sigs[sigs['upDown'] == 'up']
        down = sigs[sigs['upDown'] == 'down']
        # get a list of ids
        up = list(up['EntrezID'])
        down = list(down['EntrezID'])

        scores = score(up_gene=up, down_gene=down,
                                    sample=sample, full_data=True)
        self.assertIsInstance(plotdispersion(score=scores,show=False),
                              matplotlib.figure.Figure)

    def test_null_dist(self):

        file_path = (os.path.dirname(__file__)) + \
                    '/test_data/entrez_sample.txt'
        f = open(file_path, 'r')
        sample = pandas.read_csv(f, header='infer', sep='\t')
        sample = sample.set_index(keys=sample.columns[0])
        f.close()
        # prepare signatures
        sig_path = (os.path.dirname(__file__)) + '/test_sigs/tgfb_upDown.txt'
        sig = open(sig_path, 'r')
        sigs = pandas.read_csv(sig, header=
        'infer', sep='\t')
        sig.close()
        # subset the data for up and down
        up = sigs[sigs['upDown'] == 'up']
        down = sigs[sigs['upDown'] == 'down']
        # get a list of ids
        up = list(up['EntrezID'])
        down = list(down['EntrezID'])
        p = permutate(sample=sample[['D_Ctrl_R1']], n_up=50, n_down=50,
                      reps=10)
        scores = score(up_gene=up, down_gene=down, sample=sample[['D_Ctrl_R1']])

        self.assertIsInstance(nulldistribution(permutations=p, score=scores,
                                               show = False),
                              matplotlib.figure.Figure)


if __name__ == '__main__':
    unittest.main()