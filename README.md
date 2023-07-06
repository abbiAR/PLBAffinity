### Paper: Protein-ligand binding affinity prediction exploiting sequence constituent homology

1. CASF experiments:

1.1 DATA
For 2007, 2013 and 2016:
CASF_YEAR_aa.csv: amino acid count for binding site
CASF_YEAR_ligprop.csv: ligand features
CASF_YEAR_activities_test.csv: activities for train set
CASF_YEAR_activities_train.csv: activities for test set

For 2019:
CASF_2019_aa.csv: amino acid count for binding site
CASF_2019_ligprop.csv: ligand features
CASF_2019_activities.csv: activities

1.2 CODE
coreset.R - code for running CASF 2007, 2013 and 2016 with dedicated test sets (core sets)
not_coreset.R - code for running cross validation on CASF2019


2. CHEMBL experiments:

2.1 DATA
Experiment 1-3 csv files (ligand features and amino acid count for active site):
Number of ligand thresholds: None, 250 and 500
Experiment 4-6 csv files (ligand features and amino acid count for entire chain):
Number of ligand thresholds: None, 250 and 500

2.2 CODE
benchmark.R - code for running the experiments
