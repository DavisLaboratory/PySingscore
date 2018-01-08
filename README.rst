=========
singscore
=========
**A Python implementation of Foroutan et al (doi:10.1101/231217) simple sample scoring method for gene-set enrichment analysis**

    A package of functions that can be used to implement singscore.

    Future features will include interactivity in plots and plotting of a
    landscape of scores.

**Author:** Kristy Horan

**Brief Description:** This collection of methods will take a gene set or signature and a single sample (although multiple samples are also acceptable) and return a score reflective of the gene-set or signature enrichment for that single sample.

Install
-------
Ensure that you have git installed, in Terminal

::

    pip install git

Then
::

    pip install git+https://github.com/kristyhoran/singscore


Example workflow
----------------
The example data was originally used in Foroutan et al 2017 (DOI: 10
.1158/1541-7786.MCR-16-0313).
Data used here can be downloaded from https://github.com/kristyhoran/singscore/tree/master/singscore/test

Generate a dataframe and set index to be the Gene Identifier
column, in this case 'EntrezID'
::

    data = pandas.read_csv(open('test_data/entrez_sample.txt', 'r'), header =
    'infer', sep='\t')
    data = data.set_index(keys='EntrezID')


Get the signatures to be used

::

    sigs = pandas.read_csv(open('test_sigs/tgfb_upDown.txt', 'r'), header =
    'infer', sep = '\t')

This signature is in both up and down direction, so subset based on the
direction of the gene-set
::

    up = sigs[sigs['upDown'] == 'up']
    down = sigs[sigs['upDown'] == 'down']

and get a list of ids. Alternatively, a path to the up and/or down
signatures may be supplied. The target file should be a text file, with each
 identifier on a new line and with no headings or other characters.
::

    up = list(up['EntrezID'])
    down = list(down['EntrezID'])

scoring
::

    scored_data = score(up_gene=up, down_gene=down,sample=data,
                    norm_method='theoretical', full_data=True)


plot dispersion of scores
::

    plotdispersion(scored_data, ctrlstring='Ctrl', teststring='TGFb',
                testlabel='TGFb',colour_1='g', colour_2='b',show=True)


plotting distribution of ranks single sample, D_Ctrl_R1 -  first find ranks
::

    ranked_data = rank(up_gene=up, down_gene=down, sample=data[['D_Ctrl_R1']],
                      norm_method='theoretical')

    plotrankdist(ranks=ranked_data, colour_1='r', colour_2='b', show=True)

for significance test
::

    permd = permutate(data, n_up=193, n_down=108)
    pvals = empiricalpval(permutations=permd, score=scored_data)

    nulldistribution(permutations=permd, score=scored_data, nrows=2, ncols=5,
                  threshold=0.05, show=True)
