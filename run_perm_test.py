import pandas as pd
from src.perm_test import *
from src.utils import *
import argparse
import matplotlib.pyplot as plt
import json
import os
from pathlib import Path

pth = "results/perm_test"
Path(pth).mkdir(parents=True, exist_ok=True)


parser = argparse.ArgumentParser()
parser.add_argument("-i","--infile", help="Input file in wordlist tsv format", default='data/nostratic.tsv')
parser.add_argument("-m","--method", choices=['p1-dolgo','turchin','sca','lexstat'], default='turchin', help="Method for computing similarity metrics")
parser.add_argument("-n","--nruns", type=int, default=100, help="Num of runs")
parser.add_argument("-p", "--plotfig", action='store_true', help="Plot bilateral significance (p < 0.05), only for complete nostratic")

args = parser.parse_args()

## load data
infile = args.infile
dataset_name = infile.split('/')[-1].split('.tsv')[0]
nos_df = pd.read_csv(infile,sep='\t')
nos_words = [''] + nos_df['IPA'].tolist()
langs = nos_df['DOCULECT'].unique().tolist()
concepts = nos_df['CONCEPT'].unique().tolist()

p1_dolgo, turchin_arr = create_arr(nos_words=nos_words)

method = args.method
nruns = args.nruns

compare_func = None
wl = None
cognates = None

word_mat = np.zeros((len(langs),len(concepts)), dtype=np.int32)
if method in ['p1-dolgo', 'turchin']:
    arr = None
    if method == 'p1-dolgo':
        arr = p1_dolgo
    else:
        arr = turchin_arr
    compare_func = lambda x,y: identity(x,y,arr)
    cognates = lambda x,y: get_cognates(x,y,word_mat,langs,nos_words,arr)
elif method in ['sca', 'lexstat']:
    wl = LexStat(infile, model= 'sca')
    if method == 'lexstat':
        wl.get_scorer(runs=1000, force=True)
    langs = wl.cols
    concepts = wl.rows
    compare_func = lambda x,y: lexstat_score(x,y,wl,method=method)
    cognates = lambda x,y: get_cognates_wl(x,y,word_mat,langs,wl,method)

for index, row in nos_df.iterrows():
    word_mat[langs.index(row['DOCULECT']),concepts.index(row['CONCEPT'])] = row['ID']


print(langs)
ret_dict = perm_cluster(word_mat, langs, len(concepts), compare_func=compare_func, num_runs=nruns)

mat, p_mat = ret_dict['bi_scores'], ret_dict['bi_pval']
mat *= (mat > 0).astype(np.int32)
mat = np.abs(1-mat)
for i in range(mat.shape[0]):
    for j in range(i,mat.shape[0]):
        if i==j:
            mat[i,j] = 0
            p_mat[i,j] = 0.0
        else:
            mat[i,j] = mat[j,i]
            p_mat[i,j] = p_mat[j,i]

most_similar = lambda x: sort_by_similarity(x, mat, langs, p_mat)
print('similarity:',round(ret_dict['last_score'],3),'p-value:',round(ret_dict['last_pval'],3))

prefix = 'nos'
if 'nostratic' in args.infile:
    labels = ['Georgian', 'Old_Kannada', 'Old_Telugu', 'Old_Tamil', 'Old_Malayalam', 'Greek_Anc', 'Armenian', 'Middle_Persian','Sanskrit', 'Pali', \
         'Old_Church_Slavonic', 'Old_Irish', 'Latin', 'Old_French', 'Old_High_German', 'Old_English', 'Old_Norse']
    short_labels = ['Ge', 'Ka', 'Te', 'Ta', 'Ma', 'Gr', 'Ar', 'Pe', 'Sa', 'Pa', \
         'CS', 'Ir', 'La', 'Fr', 'HG', 'En', 'No']
else: ## mayan_mixezoque_utoaztecan
    prefix = 'mmu'
    labels = ['HUASTEC', 'TOJOLABAL', 'CHUJ', 'JACALTEC', \
       'QANJOBAL_SANTA_EULALIA', 'ACATECO_SAN_MIGUEL_ACATAN', \
       'IXIL_CHAJUL', 'AGUACATEC', 'TECO_TECTITAN', 'MAM_NORTHERN', \
       'SIPAKAPENSE', 'SACAPULTECO_SACAPULAS_CENTRO', 'CENTRAL_QUICHE', \
       'SOUTHERN_CAKCHIQUEL_SAN_ANDRES_ITZAPA', \
       'TZUTUJIL_SAN_JUAN_LA_LAGUNA', 'POQOMCHI_WESTERN', \
       'POCOMAM_EASTERN', 'USPANTEKO', 'EASTERN_KEKCHI_CAHABON', \
       'TZELTAL_BACHAJON', 'CHOL_TUMBALA', 'CHORTI', 'ITZAJ', 'MOPAN', \
       'MAYA_YUCATAN', 'ZINACANTAN_TZOTZIL', 'CHONTAL_TABASCO', \
       'LACANDON', 'MOCHO', 'CHICOMUCELTEC', 'NORTH_HIGHLAND_MIXE', \
       'SOUTH_HIGHLAND_MIXE', 'LOWLAND_MIXE', 'SAYULA_POPOLUCA', \
       'OLUTA_POPOLUCA', 'TEXISTEPEC_ZOQUE', 'SOTEAPAN_ZOQUE', \
       'MARIA_CHIMALAPA', 'MIGUEL_CHIMALAPA', 'CHIAPAS_ZOQUE', \
       'ClassicalNahuatl', 'TetelcingoNahuatl', 'NorthPueblaNahuatl', \
       'MecayapanNahuat', 'PajapanNahuat', 'JalupaNahuat', 'Pipil', \
       'Pochutec', 'ProtoAztecan']
    short_labels = [str(i+1) for i in range(len(labels))]

if args.plotfig:
    inds = [langs.index(l) for l in labels]
    related = (p_mat[inds,:][:,inds] < 0.05)
    fig = plt.pcolormesh(related, edgecolor='black', linestyle=':', lw=1)
    ax = fig.axes
    ax.invert_yaxis()
    plt.yticks(np.arange(0.5,len(labels)), short_labels)
    plt.xticks(np.arange(0.5,len(labels)), short_labels)
    plt.savefig(f'results/{prefix}_{method}.png')

save_dict = {'sim_score': ret_dict['last_score'], 'p_val': ret_dict['last_pval'], 'p_val_ls': 1/nruns, 'tree': ret_dict['tree']}
with open(f"results/perm_test/{dataset_name}_{method}.json",'w') as fp:
    json.dump(save_dict, fp)