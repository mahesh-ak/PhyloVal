import os
import regex as re
from tqdm import tqdm
import numpy as np
from pathlib import Path
import json
import sys
from scipy.stats import ttest_ind

eps = 1e-8
model = 'dolgo'
base_freq = '+F'
if model == 'sca':
    base_freq= ''

PINV = 0.06
PINV0 = 0.01

def run_ml(file, simulate=False, num_sim=100, pinv= None):
    print(f"Running ML tree search on {file}...")
    cmd = f"./tools/iqtree-2.2.2.6-Linux/bin/iqtree2 -s {file} -T 5 -redo -quiet -m poisson{base_freq}"

    if pinv:
        cmd += "+I{"+f"{pinv}"+'}'
    if simulate:
        cmd += f" --alisim {file.replace('.fa',f'_sim')} --out-format fasta"# --num-alignments {num_sim}"

    os.system(cmd)
    
    log = open(f"{file}.log",'r').read()
    match =  re.search(r'Optimal log-likelihood: (-\d+\.\d+)', log)
    likelihood = float(match.group(1))
    return likelihood



def run_param_bootstrap(family, num_sim=20):
    src_pth = "processed/aligned"
    res_pth = f"results/lrt/{family}"
    Path(res_pth).mkdir(parents=True, exist_ok=True)

    os.system(f"cp {src_pth}/{family}_{model}.fa {res_pth}/")
    file = f"{res_pth}/{family}_{model}.fa"

    print("Inferring deltas on bootstrap replicates ....")

    print(f"Running ML tree search on bootstrap reps ...")
    delta_null = []
    delta_alt = []
    for i in tqdm(range(num_sim)):
        l0 = run_ml(file, pinv= PINV0, simulate= True)
        l1 = run_ml(file, pinv= PINV)
        d = 2*(l1-l0)
        delta_alt.append(d)

        file_sim = file.replace('.fa',f'_sim.fa')#_{i+1}.fa')

        l0 = run_ml(file_sim, pinv= PINV0)
        l1 = run_ml(file_sim, pinv= PINV)
        d0 = 2*(l1-l0)

        delta_null.append(d0)
        print(f"Run {i}: delta_alt={d}, delta_null={d0}")
    return delta_alt, delta_null

def LRT(family):
    d, deltas = run_param_bootstrap(family)

    pval = ttest_ind(d,deltas,equal_var=False, alternative='greater')[1]

    print(f"Delta: {round(np.mean(d),3)}, p-val: {pval}")
    with open(f"results/lrt_{family}.json",'w') as fp:
        json.dump({'deltas': deltas, 'observed': round(np.mean(d),3), 'observed_list': d, 'p-val': pval},fp)

if __name__=='__main__':
    if len(sys.argv) != 2:
        print("Usage: python src/lrt.py family")
    else:
        family = sys.argv[1]
        LRT(family)



