#!/bin/bash

# Define an array of parameters for the R script
params=(
    # ""                          # basecase
    "_Blo"
    "_Bhi"
    "_Clo"
    "_Chi"
)

# Loop through each parameter and execute the R script
for param in "${params[@]}"; do
    echo "$param";
    R --slave --vanilla --args "$param" < 2_RRcalculations.R
done



