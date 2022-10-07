This is a repository for research purposes. 
If you want to know more about [Longevity Module](https://github.com/dna-seq/oakvar-longevity) and how you can work with it, please read [this](https://just-dna-seq.readthedocs.io/en/oakvar/) documentation.


## Calculate polygenic risk score
### 1. Preprocessing data
#### 1.1 Weights files
Check prs scoring files (for example, files from [PGS Catalog](https://www.pgscatalog.org/)) for errors:<br>
script `preprocess pgs/preprocess_prs_file.py`,<br> see --help for help with arguments,<br> 'PGS001833.txt' example in `data/`.
#### 1.2 VCF files
Check your vcf file for SNPs IDs - rsIDs. To get more info to your vcf:
- VCF files can be annotated with [oakvar](https://github.com/rkimoakbioinformatics/oakvar/)
- VCF files can be annotated with bcftools, source of annotation can be selected from [NCBI Human Variation Sets in VCF Format](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/) (ClinVar, dbSNP)<br>
See launch example in `preprocess/annotate_vcf.sh`

### 2. Calculate polygenic risk score
Check header of script `calc_and_visualize/calc_prs.py`: file paths to vcf and prs scoring file must be specified, to use the optional enable/disable some SNPs, the paths to these lists can be set.

### 3. Visualize
To draw bullet plot fill the data (custom_snp_number, overall_snp_number variables) in header `calc_and_visualize/draw_bullet_plot.py` and run.

### 3. Calculate PRS with plink
We can also perform additional check with plink.
#### 3.1 Plink
- Get plink binary files by running the script `preprocess/plink convert/get_plink_fileset_bin.sh`.<br>
If get into troubles with getting binary fileset, some intermediate steps may help. See `plink/get_plink_fileset_pgen.sh`, `plink/pgen_to_bed.sh`, or just try to fix initial .vcf file. <br>
- Run `calc_and_visualize/run_plink_prs.sh`
#### 3.2 Plink2
- Get plink pfiles by running the script `preprocess/plink convert/get_plink_fileset_pgen.sh`.
- Run `calc_and_visualize/run_plink2_prs.sh`
`plink/run_plink_prs.sh` <br>
