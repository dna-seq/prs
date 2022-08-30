#!/bin/bash
plink2 --pfile raw --score score.txt no-mean-imputation list-variants --out PRS
