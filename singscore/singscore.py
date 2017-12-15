import pandas, sys, os, numpy, pandas, matplotlib, matplotlib.pyplot, itertools, \
    seaborn, scipy, scipy.stats
# class Singscore:


def getsignature(path):
    """

    :param path: path to the signature must have no header
    :return: a list containing all the genes in the signature
    """
    sig = open(path, 'rt')
    s = []
    for line in sig.readlines():
        if line.strip().isdigit():
            s.append(int(line.strip()))
        else:
            s.append(line)

    return s

def normalisation(norm_method, score_list, score, library_len):

    """

    :param norm_method: method of normalisation, standard or theoretical
    :param score_list: list of scores will each be normalised
    :param score: average score (average of list)
    :param library_len: lengthe of the library
    :return: a tuple, containing the normalised score and an array of each
    genes score normalised
    """

    if norm_method == 'standard':
        norm = score / library_len
        u = numpy.array(score_list) / library_len
    elif norm_method == 'theoretical':
        low_bound = (library_len + 1) / 2
        upper_bound = library_len - ((library_len - 1) / 2)
        u = ((numpy.array(score_list)) - low_bound) / (upper_bound - low_bound)
        norm = (score - low_bound) / (upper_bound - low_bound)

    return norm, u

def score(up_gene, sample, down_gene = False,
          norm_method = 'standard', norm_down = 0,
          full_data= False):
    """
    This function will generate a score, using singscore method for each
    sample in a cohort. It may be used for either single direction signatures
    or both up and down.
    Gene identifiers used in signature must be the same as those used in the
    sample for example both Entrez or GeneSymbol

    :param up_gene: can be either a path to a .txt file containing genes
    (if so there must be no header) OR a list of genes
    :param sample: a dataframe, with row index as gene id (must be the same
    as the the gene identifier in the signature)
    :param down_gene: can be either a path to a .txt file containing genes
    (if so there must be no header) OR a list of genes. Default in False,
    for use in single direction signatures
    :param norm_method: choose a normalisation method. Default is
    'standard', where score is simply divided by the library size. This is
    acceptable for most case. Theoretcial uses the theoretical minimum and
    maximum to normalise the score.
    :param norm_down: if down_gene is False, then this is the value used to
    calculate the total score (0)
    :param full_data: if True then a dataframe with scores (both up and
    down if down_gene != False) and dispersion
    will be returned, otherwise just the scores will be returned

    :return: a dataframe of scores with or without dispersion
    """

    data = pandas.DataFrame()
    # ============ get up and down signatures ===============
    # if up_gene and/or down_gene are a path to txt file then gene_path =
    # True and use getsignature function to open and generate a list of
    # genes

    if type(up_gene) is str:
        up_gene = getsignature(up_gene)
        if down_gene != False:
            down_gene = getsignature(down_gene)

    for i in sample.columns:
        # rank the genes -> Ties will be taken as the rank of the first
        # appearance
        up_sort = sample[i].rank(method='min', ascending=True)
        # su is a list to contain the score for each gene in the gene list
        su = []

        # for every gene in the list gene get the value at that
        # index/rowname (the gene) and the sample that is equal to i
        for j in up_gene:
            if j in up_sort.index:
                su.append(up_sort.get_value(j, i))
        # normalise the score for the number of genes in the signature
        score_up = numpy.mean(su)

        # normalisation
        norm_up = normalisation(norm_method= norm_method, library_len=len(
            sample.index), score_list=su, score = score_up)
        u = norm_up[1]
        norm_up = norm_up[0]

        # find dispersion
        median_up = numpy.median(u)
        mad_up = numpy.median(abs(u-median_up))


        # ==== repeat with down genes,flipping the data frame around
        if down_gene != False:
            # this is the standard for scoring, opposite to up
            down_sort = sample[i].rank(method='min', ascending=False)

            # sd is a list to contain the score for each gene in the
            # gene list
            sd = []
            # for every gene in the list gene get the value at that
            # index/rowname (the gene) and the sample that is equal to i
            for k in down_gene:
                if k in sample.index:
                    sd.append(down_sort.get_value(k,i))

            score_down = numpy.mean(sd)

            # normalisation
            norm_down = normalisation(norm_method=norm_method, library_len=len(
                sample.index), score_list=sd, score=score_down)
            d = norm_down[1]
            norm_down= norm_down[0]
            # find dispersion
            median_down = numpy.median(d)
            mad_down = numpy.median(abs(d-median_down))

        total_score = norm_up + norm_down
        # make the score dataframe
        if full_data == True and down_gene != False: # if all data is
            # wanted and there is a down gene list
            temp_df = pandas.DataFrame({'up_score': norm_up,
                                    'mad_up':mad_up,'down_score':norm_down,
                                    'mad_down':mad_down,
                                    'total score': total_score,
                                    'total_mad':mad_down + mad_up},
                                   index=[i])
        elif full_data == True and down_gene == False: # if all data is
            # wanted and there is only up gene list
            temp_df = pandas.DataFrame({'score': total_score,
                                        'mad': mad_up,
                                        },index=[i])
        else: # default, regardless of down gene list, just make total
            # score
            temp_df = pandas.DataFrame({'score':total_score},
                                       index=[i])

        if len(data.columns) == 0:
            data = temp_df

        else:
            data = data.append(temp_df)


    return data


def rank(up_gene, sample, down_gene = False, gene_path = True,
          norm_method = 'standard'):

    data = pandas.DataFrame()
    # ============ get up and down signatures ===============
    # if up_gene and/or down_gene are a path to txt file then gene_path =
    # True and use getsignature function to open and generate a list of
    # genes

    if gene_path:
        up_gene = getsignature(up_gene)
        if down_gene != False:
            down_gene = getsignature(down_gene)
    su = {}
    for i in sample.columns:
        # rank the genes -> Ties will be taken as the rank of the first
        # appearance
        up_sort = sample[i].rank(method='min', ascending=True)
        # su is a dictionary to contain the gene id, rank

        su[i] = []
        # for every gene in the list gene get the value at that
        # index/rowname (the gene) and the sample that is equal to i
        for j in up_gene:
            if j in up_sort.index:
                su[i].append((j, up_sort.get_value(j,i)))


        # ==== repeat with down genes
        if down_gene != False:
            # for ranking use the same direction as up
            down_sort = sample[i].rank(method='min', ascending=True)

            # sd is a dictionary to contain the score for each gene in the
            # gene list
            sd = {}
            sd[i] = []
            # for every gene in the list gene get the value at that
            # index/rowname (the gene) and the sample that is equal to i
            for k in down_gene:
                if k in sample.index:
                    sd[i].append((k, down_sort.get_value(k,i)))

    # dataframes of ranks for each gene in each sample and then normalise
    # ranks (0 to 1, based on library size) standard = simply divide by the
    # library size theoretical = use the theoretical minimum and maximum for
    up_ranks = pandas.DataFrame({s: pandas.Series({g: r for g,r in su[s]}) for
                                s in su})
    if norm_method == 'standard':
        up_ranks = up_ranks/len(sample.index)
    elif norm_method == 'theoretical':
        low_bound_up = (len(up_gene)+1)/2
        upper_bound_up = len(sample.index) - ((len(up_gene)-1)/2)
        up_ranks = (up_ranks- low_bound_up)/(upper_bound_up-low_bound_up)



    if down_gene != False:
        down_ranks = pandas.DataFrame({s: pandas.Series({g: r for g,r in sd[s]})
                                       for s in sd})
        if norm_method == 'standard':
            down_ranks = up_ranks/len(sample.index)
        elif norm_method == 'theoretical':
            low_bound_up = (len(down_gene)+1)/2
            upper_bound_up = len(sample.index) - ((len(down_gene)-1)/2)
            down_ranks = (down_ranks- low_bound_up)/(
                upper_bound_up-low_bound_up)

        down_ranks['up_or_down'] = 'down'
        up_ranks['up_or_down'] = 'up'

    if down_gene != False:
        ranks = up_ranks.append(down_ranks)
    else:
        ranks = up_ranks

    return (ranks)



up = '../all_sigs/kh_intersect_up.txt'
down = '../all_sigs/kh_intersect_down.txt'
sample = pandas.read_csv(open('../all_data/validation_merge.txt', 'r'),
                         header = 'infer', sep = '\t')
sample = sample.set_index(keys='GeneID')

# x = score(up_gene=up, down_gene=down, sample=sample)
# print(x)

# r = rank(up_gene=up,  down_gene= down, sample= sample)
# print(r)

if type(list(up)) is str:
    print('is a list')
else:
    print('not a list')