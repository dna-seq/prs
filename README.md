## Scripts to calculate polygenic risk score
### Preprocessing data
#### Weights files
Script `preprocess pgs/preprocess_prs_file.py` for preprocessing weights files from [PGS Catalog](https://www.pgscatalog.org/), see --help for help with arguments.
#### VCF files
VCF files can be annotated with bcftools, source of annotation can be selected from [NCBI Human Variation Sets in VCF Format](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/) (ClinVar, dbSNP)<br>
See launch example in `preprocess/annotate_vcf.sh`
#### Plink files
Get plink binary files by running the script `plink/get_plink_fileset_bin.sh`.
If get into troubles with getting binary fileset, some intermediate steps may help. See `plink/get_plink_fileset_pgen.sh`, `plink/pgen_to_bed.sh`, or just try to fix initial .vcf file.
### Calculate polygenic risk score
- `plink/run_plink_prs.sh` <br>
The most preferable and simplest case: check duplicates and save them in `1.dups`, then run calculation of the score.
- python script `run_hail_prs.py` for launch [hail](https://hail.is/docs/0.2/guides/genetics.html#polygenic-score-calculation) <br>
Some java compatibility issues can easily arise.
- R script `run_Rprs.R` for launch [Rprs](https://github.com/statgen/Rprs)
