## Bayesian/ML phylogeny
## FASTA conversion
import pandas as pd
from lingpy import *
from tqdm import tqdm
import os
from typing import List, Union
from Bio import SeqIO
import sys
from pathlib import Path

model = 'dolgo'
unalgn_pth = "processed/unaligned_temp"
algn_pth = "processed/aligned_temp"
for pth in unalgn_pth, algn_pth, "processed/aligned":
    Path(pth).mkdir(parents=True, exist_ok=True)

sca = rc(model)
model_alpha = list(set([y for x,y in sca.converter.items()]))
protein_alpha = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
model2protein = {d: d for d in model_alpha}
if model == 'dolgo':
    spl = {'J': 'I', '+': 'Q', '_': 'L', '1': 'Y', '0': '-'}
elif model == 'sca':
    spl = {'J': 'G', 'U': 'V', 'O': 'A', '9': 'N', 'B': 'F', 'G': 'H', 
           '1': 'Q', '2': 'Q', '3': 'Q', '4': 'Q', '5': 'Q', '6': 'Q', '_': 'Q', '!': 'Q' ,'+': 'Q', '-':'-',
          '0':'-'}

for k, v in spl.items():
    model2protein[k] = v

    
def Model2protein(l):
    init = []
    if l[0] == 'V':
        init = ['V']
    return init + [model2protein[x] for x in l if x != 'V']


## Align and store individual meanings in algn_path
def Align(fname, algn_path= algn_pth):
    nos_df = pd.read_csv(f"data/{fname}.tsv", sep='\t')

    taxa = nos_df['DOCULECT'].unique().tolist()
    loci = nos_df['CONCEPT'].unique().tolist()

    if model.upper() in nos_df.columns:
        nos_df[model] = nos_df.apply(lambda x: ''.join(Model2protein(x[model.upper()].split())), axis=1) 
    else:
        nos_df[model] = nos_df.apply(lambda x: ''.join(Model2protein(tokens2class(ipa2tokens(x['IPA'].replace(' ','_')),model))), axis=1) 

    print("Converting meanings to FASTA ...")
    for l in tqdm(loci):
        mng = l.replace(' ','_').replace('(','').replace(')','').replace('/','_')
        lines = []
        for t in taxa:
            seq = nos_df[(nos_df["CONCEPT"] == l) & (nos_df["DOCULECT"] == t)][model].unique()
            if seq.shape[0] > 0:
                if len(seq[0]) == 1: ## Short form
                    continue
                lines.append(f"> gnl|taxa|{t}")
                lines.append(f"{seq[0]}")

            lines.append("")
        
        with open(f"processed/unaligned_temp/{mng}.fa",'w') as fp:
            fp.write('\n'.join(lines))

    files = os.listdir(unalgn_pth)

    os.system(f"rm {algn_path}/*")
    print(f"Aligning conceptes with ClustalW and saving to {algn_path} ...")
    for f in tqdm(files):
        os.system(f"./tools/clustalw-2.1/src/clustalw2 -ALIGN -MATRIX=ID -QUIET -INFILE={unalgn_pth}/{f} -OUTFILE={algn_pth}/{f} -OUTPUT=FASTA")
    os.system(f"rm {unalgn_pth}/*")
    algn_files = os.listdir(algn_pth)
    print('Issues:', set(files).difference(set(algn_files)))




## save alignments generated in algn_pth from fname with model in fasta and nexus format at processed/aligned/{fname}_{model}.*
def SaveAlign(fname, algn_path= algn_pth):

    nos_df = pd.read_csv(f"data/{fname}.tsv", sep='\t')
    taxa = nos_df['DOCULECT'].unique().tolist()

    algn_files = os.listdir(algn_path)
    concat = {lng: '' for lng in taxa}

    print(f"Combining aligned concepts from {algn_path} ...")
    for fi in tqdm(algn_files):
        records = SeqIO.parse(f"{algn_path}/{fi}","fasta")
        rec_dict = {}
        seq_len = 0
        for rec in records:
            seq = rec.seq
            seq_len = len(seq)
            langs = rec.id.split('|')[-1]
            rec_dict[langs] = seq
        for lng in concat:
            if lng in rec_dict:
                concat[lng] += rec_dict[lng]
            else:
                concat[lng] += '-'*seq_len

    lines = []
    for t, seq in concat.items():
        lines.append(f">{t}")
        lines.append(f"{seq}")
        lines.append(f"")

    with open(f"processed/aligned/{fname}_{model}.fa",'w') as fp:
        fp.write('\n'.join(lines))

    in_file = f"processed/aligned/{fname}_{model}.fa"
    out_file = f"processed/aligned/{fname}_{model}.nex"

    SeqIO.convert(in_file, "fasta", out_file, "nexus", "protein")

def AlignSave(fname):
    Align(fname)
    SaveAlign(fname)

if __name__=='__main__':
    if len(sys.argv) != 2:
        print("Usage: python src/align.py family")
    else:
        fname = sys.argv[1]
        AlignSave(fname)