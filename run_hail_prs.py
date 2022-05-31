import hail as hl


mt = hl.import_plink(
    bed="data/public/prs/ex_antonk.bed", bim="data/public/prs/ex_antonk.bim", fam="data/public/prs/ex_antonk.fam",
    quant_pheno=True, missing='-9')

mt = hl.variant_qc(mt)

scores = hl.import_table('data/public/prs/data/PGS000706_edited.txt', delimiter=' ', key='rsid',
                         types={'score': hl.tfloat32})
mt = mt.annotate_rows(**scores[mt.rsid])
flip = hl.case().when(mt.allele == mt.alleles[0], True).when(mt.allele == mt.alleles[1], False).or_missing()
mt = mt.annotate_rows(flip=flip)
mt = mt.annotate_rows(prior=2 * hl.if_else(mt.flip, mt.variant_qc.AF[0], mt.variant_qc.AF[1]))
mt = mt.annotate_cols(prs=hl.agg.sum(mt.score * hl.coalesce(hl.if_else(mt.flip, 2 - mt.GT.n_alt_alleles(),
                                                                       mt.GT.n_alt_alleles()), mt.prior)))

