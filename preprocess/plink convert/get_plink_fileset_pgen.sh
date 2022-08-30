#!/bin/bash
"""produce a binary fileset"""
plink2 --vcf  /path/example.hg37.crossmap.alt.vcf.gz --out crossmap_alt --allow-extra-chr --make-pgen --sort-vars
