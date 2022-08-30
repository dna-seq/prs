#!/bin/bash
"""produce a binary fileset"""
plink2 --vcf /path/example.hg37.pickard.vcf.gz --make-bed --out pickard_allvars --allow-extra-chr --max-alleles 2 --set-all-var-ids @#
