import pandas, sys, os, numpy, pandas, matplotlib, matplotlib.pyplot, itertools, \
    seaborn, scipy, scipy.stats
class Singscore:


    def getsignature(self, path):
        """

        :param path: path to the signature-needs to be as entrez id no header
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


    def score(self, up_gene, sample, down_gene = False, gene_path = True,
              norm_method = 'standard', ss = False, norm_down = 0,
              full_data= False):


        data = pandas.DataFrame()
        # ============ get up and down signatures ===============
        # if up_gene and/or down_gene are a path to txt file then gene_path =
        # True and use getsignature function to open and generate a list of
        # genes

        if gene_path:
            up_gene = self.getsignature(up_gene)
            if down_gene != False:
                down_gene = self.getsignature(down_gene)

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

            # normalise the score for library size
            # standard = simply divide by the library size
            # theoretical = use the theoretical minumun and maximum for
            # normalisation of library size
            if norm_method == 'standard':
                norm_up = score_up/len(up_sort)
                u = numpy.array(su)/len(up_sort)
            elif norm_method == 'theoretical':
                low_bound_up = (len(up_gene)+1)/2
                upper_bound_up = len(up_sort) - ((len(up_gene)-1)/2)
                u = ((numpy.array(su))- low_bound_up)/(upper_bound_up-low_bound_up)
                norm_up = (score_up-low_bound_up)/(upper_bound_up-low_bound_up)

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
                # normalise the score for library size
                # standard = simply divide by the library size
                # theoretical = use the theoretical minumun and maximum for
                # normalisation of library size
                if norm_method == 'standard':
                    norm_down = score_down/len(down_sort)
                    d = numpy.array(sd)/len(down_sort)
                if norm_method == 'theoretical':
                    low_bound_down = (len(down_gene) + 1) / 2
                    upper_bound_down = len(down_sort) - ((len(down_gene) - 1) / 2)
                    d = ((numpy.array(sd))-low_bound_down)/(upper_bound_down-low_bound_down)
                    norm_down = (score_down - low_bound_down) / (upper_bound_down-low_bound_down)
                # find dispersion
                median_down = numpy.median(d)
                mad_down = numpy.median(abs(d-median_down))

            total_score = norm_up + norm_down
            # make the score dataframe
            if full_data == True and down_gene != False:
                temp_df = pandas.DataFrame({'up_score': norm_up,
                                        'mad_up':mad_up,'down_score':norm_down,
                                        'mad_down':mad_down,
                                        'total score': total_score,
                                        'total_mad':mad_down + mad_up},
                                       index=[i])
            elif full_data == True and down_gene == False:
                temp_df = pandas.DataFrame({'score': total_score,
                                            'mad': mad_up,
                                            },index=[i])
            else:
                temp_df = pandas.DataFrame({'score':total_score},
                                           index=[i])

            if len(data.columns) == 0:
                # print('ok')
                data = temp_df

            else:
                # print('append')
                data = data.append(temp_df)


        return data


    def rank(self, up_gene, sample, down_gene = False, gene_path = True,
              norm_method = 'standard'):

        data = pandas.DataFrame()
        # ============ get up and down signatures ===============
        # if up_gene and/or down_gene are a path to txt file then gene_path =
        # True and use getsignature function to open and generate a list of
        # genes

        if gene_path:
            up_gene = self.getsignature(up_gene)
            if down_gene != False:
                down_gene = self.getsignature(down_gene)

        for i in sample.columns:
            # rank the genes -> Ties will be taken as the rank of the first
            # appearance
            up_sort = sample[i].rank(method='min', ascending=True)
            # su is a list to contain the gene id, rank
            su = {}
            su[i] = []
            # for every gene in the list gene get the value at that
            # index/rowname (the gene) and the sample that is equal to i
            for j in up_gene:
                if j in up_sort.index:
                    # su.append((j,up_sort.get_value(j, i)))
                    su[i] = su[i].append((j, up_sort.get_value(j,i)))
                    # normalise the rank for library size
            # standard = simply divide by the library size
            # theoretical = use the theoretical minumun and maximum for
            # normalisation of library size
            if norm_method == 'standard':
                u = numpy.array(su)/len(up_sort)
            elif norm_method == 'theoretical':
                low_bound_up = (len(up_gene)+1)/2
                upper_bound_up = len(up_sort) - ((len(up_gene)-1)/2)
                u = ((numpy.array(su))- low_bound_up)/(upper_bound_up-low_bound_up)

                        # ==== repeat with down genes,flipping the data frame around
            if down_gene != False:
                # this is the standard for scoring, opposite to up
                down_sort = sample[i].rank(method='min', ascending=True)

                # sd is a list to contain the score for each gene in the
                # gene list
                sd = []
                # for every gene in the list gene get the value at that
                # index/rowname (the gene) and the sample that is equal to i
                for k in down_gene:
                    if k in sample.index:
                        sd.append(down_sort.get_value(k,i))

                score_down = numpy.mean(sd)
                # normalise the score for library size
                # standard = simply divide by the library size
                # theoretical = use the theoretical minumun and maximum for
                # normalisation of library size
                if norm_method == 'standard':
                    norm_down = score_down/len(down_sort)
                    d = numpy.array(sd)/len(down_sort)
                if norm_method == 'theoretical':
                    low_bound_down = (len(down_gene) + 1) / 2
                    upper_bound_down = len(down_sort) - ((len(down_gene) - 1) / 2)
                    d = ((numpy.array(sd))-low_bound_down)/(upper_bound_down-low_bound_down)
                    norm_down = (score_down - low_bound_down) / (upper_bound_down-low_bound_down)
                # find dispersion
                median_down = numpy.median(d)
                mad_down = numpy.median(abs(d-median_down))

            total_score = norm_up + norm_down