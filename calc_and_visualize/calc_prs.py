import pandas as pd
import gzip
import os

pwd = os.getcwd()
vcf_file_path = '/path_to/example.vcf'
sample_name = 'default'  # last column: column after 'FORMAT' in vcf
pgs_file_path = 'path_to/prs_scores.txt'
pre_snp_list_path = '/path_to/snp_list.txt'  # to get in analysis only snps from this list, leave blank '' if no list
nopred_path = '/path_to/no_pred.txt'  # to exclude from analysis snps from this list, leave blank '' if no list
ploidy = 2
pre_snps = set()
no_pred_snps = set()
processed_snps = set()


def get_vcf_header(vcf_path):
    try:
        with gzip.open(vcf_path, "rt") as ifile:
            for line in ifile:
                if line.startswith("#CHROM"):
                    vcf_names = [x for x in line.split('\t')]
                    break
    except gzip.BadGzipFile:
        with open(vcf_path, "rt") as ifile:
            for line in ifile:
                if line.startswith("#CHROM"):
                    vcf_names = [x for x in line.split('\t')]
                    break
    return vcf_names


def get_zygosity(df, idx_vcf):
    zygos_dict = dict()
    for idx in idx_vcf:
        zygo = (df.iloc[idx]['{}\n'.format(sample_name)].split(',')[0]).split(':')[0]
        ref = df.iloc[idx]['REF']
        alt = df.iloc[idx]['ALT']
        rsid = df.iloc[idx]['ID']
        print("zygo", zygo, "ref", ref, "alt", alt, 'ID', df.iloc[idx]['ID'])
        if no_pred_snps and rsid in no_pred_snps:
            zygos_dict[rsid] = ('missing', '')
            print("rsid {} from nopred list, insert missing genotype to skip in analysis".format(rsid))
            continue
        """
        - 0/0 : the sample is homozygous reference
        - 0/1 : the sample is heterozygous, carrying 1 copy of each of the REF and ALT alleles
        - 1/1 : the sample is homozygous alternate
        """
        if zygo == '0/0' or zygo == '0|0':
            zygos_dict[rsid] = ('homo', [ref])
        if zygo == '1/1' or zygo == '1|1':
            zygos_dict[rsid] = ('homo', [alt])
        if zygo == './.':
            zygos_dict[rsid] = ('missing', '')
        else:  # zygo == '0/1' or zygo == '0|1' or zygo == '1/2':
            zygos_dict[rsid] = ('hetero', [ref, alt])
    return zygos_dict


def read_vcf(vcf_path, header_line_names):
    try:
        with pd.read_csv(vcf_path, compression='gzip', comment='#', chunksize=1000, delim_whitespace=True,
                         header=None,
                         names=header_line_names) as vcf_reader:
            for df_chunk in vcf_reader:
                yield df_chunk
    except gzip.BadGzipFile:
        with pd.read_csv(vcf_path, comment='#', chunksize=1000, delim_whitespace=True,
                         header=None,
                         names=header_line_names) as vcf_reader:
            for df_chunk in vcf_reader:
                yield df_chunk


def calc_weight(zygos_dict, rsid, df_pgs, idx, score_sum, non_miss_snp_count):
    weight = df_pgs.iloc[idx]['effect_weight']
    if weight is None or weight == "":
        print('None weight for ', rsid)
        return
    if zygos_dict[rsid][0] == 'homo':
        weight = 2 * weight
    elif zygos_dict[rsid][0] == 'missing':
        # weight = 2 * frequency* weight
        print("missing genotype for", rsid)
        weight = 0
    if weight:
        score_sum += weight
        non_miss_snp_count += 1
        return score_sum, non_miss_snp_count


def calc_prs(vcf_path, pgs_path, header_line_names, pre_snp_list):
    score_sum = 0
    non_miss_snp_count = 0
    df_pgs = pd.read_csv(pgs_path, delimiter="\t")
    chunks = 0
    for df_chunk in read_vcf(vcf_path, header_line_names):
        chunks += 1
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
            rsid = df_pgs.iloc[idx]['rsID']
            vcf_alleles = zygos_dict[rsid][1]
            effect_allele = df_pgs.iloc[idx]['effect_allele']
            if pre_snp_list:
                if rsid in pre_snp_list:
                    try:
                        score_sum, non_miss_snp_count = calc_weight(zygos_dict, rsid, df_pgs, idx, score_sum,
                                                                    non_miss_snp_count)
                        processed_snps.add(rsid)
                    except TypeError:
                        pass
            else:
                for allele in vcf_alleles:
                    if effect_allele in allele:
                        try:
                            score_sum, non_miss_snp_count = calc_weight(zygos_dict, rsid, df_pgs, idx, score_sum,
                                                                        non_miss_snp_count)
                            processed_snps.add(rsid)
                            break
                        except TypeError:
                            pass
    if non_miss_snp_count:
        score_avg = score_sum / (ploidy * non_miss_snp_count)
        print('sum:', score_sum, ", count:", non_miss_snp_count, ', average:', score_avg)
        return non_miss_snp_count, score_sum, score_avg
    else:
        print("No match: None weights")


def write_result(non_miss_snp_count, score_sum, score_avg):
    prs_file_name = (pgs_file_path.split('/')[-1]).split('.')[0]
    score_res_path = os.path.join(pwd, '{}.{}'.format(prs_file_name, 'score'))
    df = pd.DataFrame({'snp_count': [non_miss_snp_count], 'score_sum': [score_sum], 'average': [score_avg]})
    print("score result recorded to {}".format(score_res_path))
    df.to_csv(score_res_path, '\t', index=False)
    processed_path = os.path.join(pwd, '{}.{}'.format(prs_file_name, 'vars'))
    df = pd.DataFrame(list(processed_snps))
    df.to_csv(processed_path, '\t', index=False, header=False)
    print("scoring snp list recorded to {}".format(processed_path))


if __name__ == "__main__":
    vcf_header = get_vcf_header(vcf_file_path)
    print("vcf_header", vcf_header)
    if pre_snp_list_path:
        with open(pre_snp_list_path, 'r') as f:
            for snp in f:
                pre_snps.add(snp.replace('\n', ''))
        print('snp_list of length: {}'.format(len(pre_snps)))
    if nopred_path:
        with open(nopred_path, 'r') as f:
            for snp in f:
                no_pred_snps.add(snp.replace('\n', ''))
        print('no_pred_snps of length: {}'.format(len(no_pred_snps)))
    try:
        non_missing_snp_count, prs_score_sum, prs_score_avg = calc_prs(vcf_file_path, pgs_file_path, vcf_header, pre_snps)
        write_result(non_missing_snp_count, prs_score_sum, prs_score_avg)
    except TypeError:
        print('No intersection of rsids in vcf and rsids in prs scoring file, or missing weights')