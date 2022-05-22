## Scripts to run polygenic risk score
### Preprocessing data
#### Weights files
Script `preprocess_prs_file.py` for preprocessing weights files from [PGS Catalog](https://www.pgscatalog.org/), see --help for help with arguments.
#### Plink files
Get plink binary files by running the script `get_plink_fileset.sh`
### Calculate polygenic risk score
- python script `run_hail_prs.py` for launch [hail](https://hail.is/docs/0.2/guides/genetics.html#polygenic-score-calculation)
- R script `run_Rprs.R` for launch [Rprs](https://github.com/statgen/Rprs)
