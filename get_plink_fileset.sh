#!/bin/bash
"""produce a binary fileset"""
plink2 --vcf antonkulaga.vcf.gz --make-bed --out ex_antonk --allow-extra-chr --max-alleles 2

