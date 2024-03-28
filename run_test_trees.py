import subprocess
import time
import os

test_files = ['data/test_aa.tsv', 'data/test_an.tsv', 'data/test_ie.tsv', 'data/test_pn.tsv', 'data/test_st.tsv']
methods = ['p1-dolgo','turchin','lexstat', 'sca']
nruns = 100

for d in test_files:
    skip = True
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
        time.sleep(5400) ## sleep for 1.5 hr before running next file