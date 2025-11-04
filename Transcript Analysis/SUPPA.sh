#!/bin/bash

# Script to run SUPPA2 for differential splicing analysis

# Step 1: Calculate PSI per event
suppa.py psiPerEvent -i "./model.ioe" -e "./suppa_expression.txt" -o events

# Step 2: Create subset files for the "60-70" group
cut -f 1-11 ./events.psi > ./60_70.psi
cut -f 1-11 ./suppa_expression.txt > ./60_70.tpm

# Step 3: Create subset files for the "30-40" group
cut -f 1,12-20 ./events.psi > ./30_40.psi
cut -f 1,12-20 ./suppa_expression.txt > ./30_40.tpm

# Step 4: Run differential splice analysis
# This step uses SUPPA to perform differential splicing analysis between two conditions.
# Before running this step, ensure that the headers and data in your input files (PSI and TPM) are correctly matched.
# Mismatched headers and data may cause SUPPA to fail.
suppa.py diffSplice --method empirical --input ./model.ioe --psi ./30_40.psi ./60_70.psi --tpm ./30_40.tpm ./60_70.tpm -gc -o dpsi

# Step 5: Filter the output to remove 'nan' values, keeping the header
awk -F'\t' 'NR==1 || ($2 != "nan") { print }' ./dpsi.dpsi > ./output2.txt

echo "Differential splicing analysis finished successfully!"
