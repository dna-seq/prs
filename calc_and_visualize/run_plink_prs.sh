#!/bin/bash
cut -f 2 pickard_allvars.bim | sort | uniq -d > 1.dups
# https://www.cog-genomics.org/plink/1.9/score
# By default, final scores are averages of valid per-allele scores. The 'sum' modifier causes sums to be reported instead
plink --bfile pickard_allvars --exclude 1.dups --score ../PGS001298_chrpos.raw sum --allow-extra-chr
