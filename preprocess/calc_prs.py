import numpy as np
import pandas as pd
import gzip

vcf_file_path = '/home/alina_grf/Documents/longevity/PRS/antonkulaga_38_all_rsid_chrompos.vcf.gz'
pgs_file_path = '/home/alina_grf/Documents/longevity/PRS/pgs/obesity/PGS001298_header.txt'
# nopred_path = '/home/alina_grf/Documents/longevity/PRS/pgs/obesity/launches/launch4_chrom_pos/plink.nopred'
nopred_path = ''
ploidy = 2


def get_vcf_header_line_names(vcf_path):
    with gzip.open(vcf_path, "rt") as ifile:
        for line in ifile:
            if line.startswith("#CHROM"):
                vcf_names = [x for x in line.split('\t')]
                break
    ifile.close()
    return vcf_names


def get_zygosity(df, idx_vcf):
    no_pred = list()
    if nopred_path:
        df_nopred = pd.read_csv(nopred_path,  names=['nosnp', 'rsID'], header=None, delimiter='\t')
        no_pred = (df_nopred.loc[:, 'rsID']).to_list()
    zygos_dict = dict()
    for idx in idx_vcf:
        zygo = (df.iloc[idx]['default\n'].split(',')[0]).split(':')[0]
        ref = df.iloc[idx]['REF']
        alt = df.iloc[idx]['ALT']
        id = df.iloc[idx]['ID']
        print("zygo", zygo, "ref", ref, "alt", alt, 'ID', df.iloc[idx]['ID'])
        if no_pred and id in no_pred:
            zygos_dict[df.iloc[idx]['ID']] = ('missing', False)
            print("nopred {}".format(id))
            continue
        if zygo == '0/0' or zygo == '1/1':
            zygos_dict[df.iloc[idx]['ID']] = ('homo', alt)
        elif zygo == '0/1':
            zygos_dict[df.iloc[idx]['ID']] = ('hetero', alt)
        else:
            zygos_dict[df.iloc[idx]['ID']] = ('missing', False)
            print("missing genotype")
    return zygos_dict


def calc_prs(vcf_path, pgs_path, header_line_names):
    prs_score_sum = 0
    non_missing_snp_count = 0
    df_pgs = pd.read_csv(pgs_path, delimiter="\t")
    chunks = 0
    with pd.read_csv(vcf_path, compression='gzip', comment='#', chunksize=1000, delim_whitespace=True,
                     header=None,
                     names=header_line_names) as vcf_reader:
        for df_chunk in vcf_reader:
            chunks += 1
            # print(df_chunk['#CHROM'], df_chunk['ID'], df_chunk['REF'], df_chunk['ALT'], df_chunk['default\n'])

            isin_pgs = df_pgs['rsID'].isin(df_chunk['ID'])
            isin_vcf = df_chunk['ID'].isin(df_pgs['rsID'])

            idx_match_pgs = df_pgs.index[isin_pgs == True].tolist()
            if idx_match_pgs:
                if chunks > 1:
                    idx_match_vcf = (df_chunk.index[isin_vcf == True] - 1000*chunks).tolist()
                else:
                    idx_match_vcf = df_chunk.index[isin_vcf == True].tolist()
                zygos_dict = get_zygosity(df_chunk, idx_match_vcf)
            else:
                continue
            for idx in idx_match_pgs:
                weight = df_pgs.iloc[idx]['WEIGHT_eff']
                rsid = df_pgs.iloc[idx]['rsID']
                effect_allele = df_pgs.iloc[idx]['EA']
                print('rsID:', rsid, "EA:", effect_allele, "OA:", df_pgs.iloc[idx]['OA'])
                if zygos_dict[rsid][1] != effect_allele:
                    print("warning", zygos_dict[rsid][1], effect_allele)
                    continue
                if zygos_dict[rsid][0] == 'hom':  # and zygos_dict[rsid][1] == effect_allele:
                    weight = 2 * weight
                # elif zygos_dict[rsid][0] == 'hetero' and zygos_dict[rsid][1] == effect_allele:
                #     print("hetero {}".format(rsid))
                elif zygos_dict[rsid][0] == 'missing':
                    # weight = 1.5 * weight
                    weight = 0
                if weight:
                    prs_score_sum += weight
                    non_missing_snp_count += 1
    if non_missing_snp_count:
        prs_score_avg = prs_score_sum / (ploidy*non_missing_snp_count)
        print("prs_score_avg", prs_score_avg, 'prs_score_sum', prs_score_sum, "non_missing_snp_count", non_missing_snp_count)
    else:
        print("no match")


names = get_vcf_header_line_names(vcf_file_path)
print("names", names)
calc_prs(vcf_file_path, pgs_file_path, names)
# for row in df.itertuples(index=True, name='Pandas'):
#     print(getattr(row, "Name"), getattr(row, "Percentage"))
