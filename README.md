# Code for "A Likelihood Ratio Test of Genetic Relationship among Languages"

V.S.D.S. Mahesh Akavarapu, Arnab Bhattacharya

Runs on Linux 64-bit
The datafiles are present in `data/`. Pip requirements are in `requirements.txt`.
Note that all the test is run on various the families and methods in parallel, hence make sure that resources are sufficient.
Tested on 12 core CPU 3.2 GHz and 32 GB RAM
The programs re-run only if the relavent result files are missing. We provide all of them beforehand. If one wishes to re-run from scratch clear the directory `results/` and run the experiments. 

## Running LRT

To generate the results for LRT saved in `results/lrt_<family>.json`, run the following command (should take less an hour):

>`python run_lrt.py`

## Running Multilateral Permutation Tests

To generate the results for multilateral permutation tests, run the following command (Should take about 8 hrs):

>`python run_nos.py`

Individual families and methods can be tested by the running `run_perm_test.py`, run the following command for help:

>`python run_perm_test.py --help`

## Running Tree Construction Task

To generate the results for tree construction task, run the following command (should take about 8 hrs):

>`python run_test_trees.py`

This command, however does not generate ML-trees, they have to be manually generated on the alignment files generated by `src/align.py` at `processed/aligned/`,
on the families `test_aa`, `test_an`, `test_ie`, `test_pn` and `test_st`

For example, for family `test_aa`, run `python src/align.py test_aa` and then manually analyze maximum likelihood tree from `processed/aligned/test_aa_dolgo.fa`
by using MEGA11 GUI (install from [here](https://www.megasoftware.net/)) and the model POISSON+I+G2 by NNI search (fastest) and save the Newick tree as `results/test_aa_ml.tre`

## Tabulating the results

Once all the required files are generated, the results can be tabulated by running:

>`python tabulate_res.py`

The results are tabulated in `results/summary_<task>.tsv`. The primary results as well as those for evaluation of macro-families are generated together.ho

Analysis on Nostratic is performed manually using the `.fa` files from `processed/aligned/` on MEGA11 GUI with model Poisson+I.

