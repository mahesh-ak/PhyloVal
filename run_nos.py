import subprocess
import time
import os

data_files_sets = [['data/drav.tsv', 'data/ie.tsv', 'data/drav_ie.tsv', 'data/mayan.tsv'],\
              ['data/mixezoque.tsv','data/mayan_mixezoque.tsv', 'data/utoaztecan.tsv'], ['data/mayan_utoaztecan.tsv', 'data/monkhmer1_utoaztecan.tsv',\
              'data/afrasian_loloburmese.tsv'], ['data/monkhmer.tsv','data/munda.tsv', 'data/monkhmer_munda.tsv', 'data/monkhmer1_mayan.tsv']]
methods = ['p1-dolgo','turchin','lexstat', 'sca']
nruns = 100

for m in methods:
    if os.path.isfile(f"results/perm_test/nostratic_{m}.json"):
        continue

    if m in ['p1-dolgo', 'turchin']:
            nruns = 1000
    else:
            nruns = 100
    subprocess.Popen(f'python run_perm_test.py -i data/nostratic.tsv -m {m} -p -n {nruns}'.split()) 

for data_files in data_files_sets:
    skip = True
    for d in data_files:
        for m in methods:
            family = d.split('/')[-1].replace('.tsv','')
            if os.path.isfile(f"results/perm_test/{family}_{m}.json"):
                continue
            else:
                skip = False

            if m in ['p1-dolgo', 'turchin']:
                nruns = 1000
            else:
                nruns = 100
            subprocess.Popen(f'python run_perm_test.py -i {d} -m {m} -n {nruns}'.split())
    if not skip:
        time.sleep(7200) ## sleep for 2 hr before running next files set

