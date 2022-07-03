#!/bin/bash
# allvars_dbsnp.vcf.gz can be downloaded from https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/
bcftools index allvars_dbsnp.vcf.gz
bcftools view example_38.vcf -Oz -o example_38.vcf.gz
bcftools index example_38.vcf.gz
bcftools annotate -c ID -a allvars_dbsnp.vcf.gz example_38.vcf.gz -o example_allvars_dbsnp.vcf.gz
bcftools annotate --set-id +'%CHROM\_%POS' example_allvars_dbsnp.vcf.gz -o example_allvars_dbsnp_chrompos.vcf.gz