#!/bin/bash

plink2 --bfile filtered -indep-pairwise 1000 100 0.8 --out 02/pruned/pruned
