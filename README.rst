=========
singscore
=========
**A Python implementation of Foroutan et al (doi:10.1101/231217) simple sample scoring method for gene-set enrichment analysis**

    A package of functions that can be used to implement singscore.

    Future features will include interactivity in plots and plotting of a
    landscape of scores.

**Author:** Kristy Horan

**Brief Description:** This collection of methods will take a gene set or signature and a single sample (although multiple samples are also acceptable) and return a score reflective of the gene-set or signature enrichment for that single sample.

Methods
-------

score(up_gene, sample, down_gene = False,norm_method = 'standard',norm_down = 0,full_data= False)

    This function will generate a score, using singscore method for each
    sample in a cohort. It may be used for either single direction signatures or both up and down. Gene identifiers used in signature must be the same as those used in the
    sample for example both Entrez or GeneSymbol

rank(up_gene, sample, down_gene = False,norm_method = 'standard')

    This function will generate a dataframe of ranks, using singscore method
    for each gene in each sample in a cohort. It may be used for either single
    direction signatures or both up and down. The difference between rank
    and score is that the rank simply shows where in the library the genes
    of a signature are placed high normalised rank means that the gene is
    highly expressed, whereas a low rank means the gene is lowly expressed.

    Gene identifiers used in signature must be the same as those used in the
    sample for example both Entrez or GeneSymbol


plotrankdist(ranks, nrows= 1, ncols = 1, counter = 0, t = False, colour_1 ='black', colour_2 = 'grey', singledir = True, output =False, show= True)

    Takes the dataframe output of rank and produces barcode plots. options to
    save figure to output path (if supplied) and show figure.
    Must supply nrows and ncols that are sutiable for the number of plots to
    generated.

plotdispersion(score, nrows = 1, ncols = 1, counter = 0, ctrlstring =False,teststring= False, testlabel = False, colour_1 ='grey', colour_2 = 'black', outpath = False, show = True)

    Takes a dataframe, the output of score(), with full_data set to True. It
    will plot the MAD of up, down and total scores if supplied or just the
    MAD of total score

permutate(sample, n_up, n_down = False, reps= 100, norm_method ='standard', rs_down = 0)

    Bootstrap a random population of scores for a given sample, that is
    dependent on a the number of genes in a signature. Take a sample and
    score it with randomly selected genes from a gene list. Returns a
    dataframe the length of the permutations (reps) desired.

empiricalpval(permutations, score)

    The empirical p value is the probability of observing a score greater
    than observed based on the permutation of the sample.

    p = (r + 1)/(m + 1), where r is the number of scores in the permutated
    data that is greater than the actual score and m is the number of
    permutations


nulldistribution(permutations, score,  nrows = 1, ncols = 1,counter = 0, outpath = False, show = True, color ='grey', threshold = False)

    Generate histogram/density plot of permutated data, with actual score
    indicated by a vertical line (blue) and a significance threshold
    indicated by red vertical line

