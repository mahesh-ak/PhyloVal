import json, subprocess
import regex as re
import pandas as pd
from ete3 import PhyloTree
from pathlib import Path

pth = "results/pred_trees"
Path(pth).mkdir(parents=True, exist_ok=True)



################################ MPT with LRT #################################################################################
families = ['drav', 'ie', 'drav_ie', 'nostratic', 'mayan', 'mixezoque', 'utoaztecan', 'mayan_mixezoque', 'mayan_utoaztecan', 'monkhmer', 'munda', 'monkhmer_munda', 'monkhmer1_utoaztecan' , 'monkhmer1_mayan', 'afrasian_loloburmese']
methods = ['p1-dolgo','turchin','lexstat', 'sca']

lines = ["Method\t" + '\t'.join(families)]
for m in methods:
    line = [m]
    for family in families:
        fname = f"results/perm_test/{family}_{m}.json"
        with open(fname,'r') as fp:
            scores = json.load(fp)
        p_val = scores['p_val']
        if p_val == 0.0:
            p_val = f"<{scores['p_val_ls']}" 
        else:
            p_val = f"{p_val}"
        line.append(f"{round(scores['sim_score'],3)} ({p_val})")
    lines.append('\t'.join(line))

line = ['lrt_pinv']
for family in families:
    fname = f"results/lrt_{family}.json"
    with open(fname,'r') as fp:
            scores = json.load(fp)
    p_val = scores['p-val']
    if p_val < 0.001:
        p_val = f"<0.001" 
    else:
        p_val = f"{p_val:.3f}"
    line.append(f"{round(scores['observed'],3)} ({p_val})")
lines.append('\t'.join(line))

with open("results/summary_lrt.tsv",'w') as fp:
    fp.write('\n'.join(lines))


############################ Tree benchmark ####################################
fnames = ['test_aa', 'test_an', 'test_ie', 'test_pn', 'test_st']
methods = ['p1-dolgo','turchin','lexstat','sca', 'ml-poisson+I+G2']

map_df = pd.read_csv("data/map_glottocode.txt", sep='\t')
map_dict = dict(zip(map_df.lang, map_df.glottcode))

for f in fnames:
    for m in methods:
        if m == 'ml-poisson+I+G2':
            tree = open(f"results/{f}_ml.tre",'r').read()
        else: 
            with open(f"results/perm_test/{f}_{m}.json",'r') as fp:
                res_dict = json.load(fp)
            tree = res_dict['tree']+';'
        tree = re.sub(r'(\d)*(\.)*(\d)*:(\d)+\.(\d)+','',tree)
        for k,v in map_dict.items():
            tree = tree.replace(f"({k},",f"({str(v)},")
            tree = tree.replace(f",{k})",f",{str(v)})")
            tree = tree.replace(f",{k},",f",{str(v)},")
        with open(f"results/pred_trees/{f}_{m}.tre",'w') as fp:
            fp.write(tree)

lines = ["Method\tAA\tAN\tIE\tPN\tST"]
for m in methods:
    line = [m]
    for f in fnames:
        t1 = PhyloTree(f"results/pred_trees/{f}_{m}.tre")
        t2 = PhyloTree(f"data/gold_trees/{f}.tre")
        leaves_2 = [nm for nm in t2.get_leaf_names()]
        node2labels = t1.get_cached_content(store_attr="name")
        t1.set_species_naming_function(lambda node: node.name)
        t1 = t1.collapse_lineage_specific_expansions()
        t1.prune(leaves_2)
        leaves_1 = [nm for nm in t1.get_leaf_names()]
        set_1 = set(leaves_1)
        set_2 = set(leaves_2)
        duplicates_1 = [l for l in set_1 if leaves_1.count(l) > 1]
        t1.write(format=9, outfile=f"results/pred_trees/{f}_{m}.tre")
    
        qdist = subprocess.check_output(["tools/qdist-src-2.0/qdist", f"results/pred_trees/{f}_{m}.tre", f"data/gold_trees/{f}.tre"])
        x=str(qdist).split("\\n")[1].split("\\t")
        gqd = float(x[4])/float(x[2])
        line.append(f"{round(gqd,3)}")
    lines.append('\t'.join(line))

with open("results/summary_test_gqd.tsv",'w') as fp:
    fp.write('\n'.join(lines))
