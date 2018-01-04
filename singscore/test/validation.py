import pandas
from singscore.singscore import *

# compare score output to R script

data = pandas.read_csv(open('test_data/entrez_sample.txt', 'r'), header =
'infer', sep='\t')
data = data.set_index(keys='EntrezID')

# prepare signatures
sigs = pandas.read_csv(open('test_sigs/tgfb_upDown.txt', 'r'), header =
'infer', sep = '\t')
# subset the data for up and down
up = sigs[sigs['upDown'] == 'up']
down = sigs[sigs['upDown'] == 'down']
# get a list of ids
up = list(up['EntrezID'])
down = list(down['EntrezID'])

# scoring
scored_data = score(up_gene=up, down_gene=down,sample=data,
                    norm_method='theoretical', full_data=True)
print(scored_data)

# plot dispersion of scores
plotdispersion(scored_data, ctrlstring='Ctrl', teststring='TGFb',
               testlabel='TGFb',colour_1='g', colour_2='b',
               outpath='output/dispersion.pdf',show=False)


# plotting distribution of ranks single sample, D_Ctrl_R1 -  first find ranks
ranked_data = rank(up_gene=up, down_gene=down, sample=data[['D_Ctrl_R1']],
                   norm_method='theoretical')
# plot
plotrankdist(ranks=ranked_data, colour_1='r', colour_2='b',
             output='output/barcode.pdf', show=False)

# for significance test
permd = permutate(data, n_up=193, n_down=108)
pvals = empiricalpval(permutations=permd, score=scored_data)
print(pvals)
nulldistribution(permutations=permd, score=scored_data, nrows=2, ncols=5,
                 threshold=0.05, outpath='output/nulldist.pdf', show=False)
