library(Rprs)

setwd("/home/rstudio/data")
pgs_files = list.files(path = ".", pattern = 'PGS(.+).txt')
for (f in pgs_files)
{
  p <- prs_gw(f, "antonkulaga_37.vcf.gz", NULL)
  outputfile <- paste(getwd(), paste(paste('prs', unlist(strsplit(f, "\\."))[1], sep='_'), 'csv', sep="."), sep="/")
  write.table(p, outputfile, col.names = TRUE, row.names = FALSE, sep = "\t", quote = F)
}
