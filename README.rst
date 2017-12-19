=========
singscore
=========
**A Python implementation of Foroutan et al (doi:10.1101/231217) simple sample scoring method for gene-set enrichment analysis**

**Author:** Kristy Horan

**Brief Description:** This collection of methods will take a gene set or signature and a single sample (although multiple samples are also acceptable) and return a score reflective of the gene-set or signature enrichment for that single sample.

Methods
-------

score(up_gene, sample, down_gene = False,norm_method = 'standard',norm_down = 0,full_data= False)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

rank(up_gene, sample, down_gene = False,norm_method = 'standard')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

plotrankdist(ranks, nrows= 1, ncols = 1, counter = 0, t = False, colour_1 ='black', colour_2 = 'grey', singledir = True, output =False, show= True)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

plotdispersion(score, nrows = 1, ncols = 1, counter = 0, ctrlstring =False,teststring= False, testlabel = False, colour_1 ='grey', colour_2 = 'black', outpath = False, show = True)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^