cut -f 2 pickard_allvars.bim | sort | uniq -d > 1.dups
plink --bfile pickard_allvars --exclude 1.dups --score ../PGS001298_chrpos.raw --allow-extra-chr
