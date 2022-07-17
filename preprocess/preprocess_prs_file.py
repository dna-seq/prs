import argparse
import os.path
import csv

"""
Script for Unix-like OS
"""


def parse_dir(input_dir):
    for infile in os.scandir(input_dir):
        if infile.name.split('.')[-1] == 'txt':
            yield os.path.join(input_dir, infile.name)


def get_pgs_header_line_names(pgs_path):
    with open(pgs_path, "rt") as ifile:
        for line in ifile:
            if not line.startswith("CHROM"):
                pgs_names = [x for x in line.split('\t')]
                break
    ifile.close()
    return pgs_names


def get_edited_file_path(output_dir, file_path):
    path_splitted = file_path.split('/')
    name_pgs = path_splitted[-1].split('.')
    cleaned_file_name = '{}_edited.{}'.format(name_pgs[0], name_pgs[1])
    cleaned_file_path = os.path.join(output_dir, cleaned_file_name)
    return cleaned_file_path


def get_rsid_idx(splited_line):
    def split_line_casefold(x):
        return [xx.casefold() for xx in x]
    if 'rsid' in split_line_casefold(splited_line):
        index = split_line_casefold(splited_line).index('rsid')
        return index
    else:
        return None


def change_header(header):
    s1 = header.replace('chr_name', 'CHROM')
    s2 = s1.replace('chr_position', 'POS')
    s3 = s2.replace('effect_allele', 'EA')
    s4 = s3.replace('other_allele', 'OA')
    edited_header = s4.replace('effect_weight', 'WEIGHT_eff')
    return edited_header


def change_header_without_rsid(header):
    s1 = header.replace('chr_name', 'CHROM')
    s2 = s1.replace('chr_position', 'POS')
    s3 = s2.replace('effect_allele', 'EA')
    s4 = s3.replace('other_allele', 'OA')
    s5 = s4.replace('effect_weight', 'WEIGHT_eff')
    edited_header = 'rsID\t' + s5
    return edited_header


def move_empty_rsid(prs_file_path, output_dir, replace_flag):
    rs_id_list = list()
    sniffer = csv.Sniffer()
    edited_file_path = get_edited_file_path(output_dir, prs_file_path)
    if os.path.isfile(edited_file_path):
        raise IOError("Edited file path {} already exists, remove or rename it first".format(edited_file_path))
    i = -1
    without_rsid = False
    with open(prs_file_path, 'r') as f:
        with open(edited_file_path, 'a') as output:
            for line in f:
                if not line.startswith('#') and i == -1:
                    dialect = sniffer.sniff(line)
                    s = line.split(dialect.delimiter)
                    rsid_idx = get_rsid_idx(s)
                    if rsid_idx is not None:
                        new_header = change_header(line)
                        output.write(new_header)
                        print("new header recorded")
                        i += 1
                        continue
                    else:
                        print("no rsid column in {}, rsid column will be inserted and values will be replaced by "
                              "CHROM_POS ".format(prs_file_path))
                        new_header = change_header_without_rsid(line)
                        output.write(new_header)
                        print("new header recorded")
                        rsid_idx = 0
                        without_rsid = True
                        i += 1
                        continue
                if not line.startswith('#') and i >= 0:
                    s = line.split(dialect.delimiter)
                    if (not s[rsid_idx] or without_rsid) and replace_flag == 'y':
                        i += 1
                        print("new_header", repr(new_header))
                        chrom, pos = new_header.split('\t').index('CHROM') - 1, new_header.split('\t').index('POS') - 1
                        print("chrom, pos", chrom, pos)
                        print("s", s)
                        if s[chrom] and s[pos]:
                            s[rsid_idx] = '{}_{}'.format(s[chrom], s[pos])
                        else:
                            print("error file: no CHROM and POS")
                            return edited_file_path, i
                        edited_line = dialect.delimiter.join(s)
                        output.write(edited_line)
                    else:
                        output.write(line)
                    if s[rsid_idx] in rs_id_list:
                        print("warning: {} is duplicate".format(s[rsid_idx]))
                    rs_id_list.append(s[rsid_idx])
    return edited_file_path, i


def main(prs_dir, out_dir, replace_flag):
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    for prs_file_path in parse_dir(prs_dir):
        print("start checking file {}".format(prs_file_path))
        edited_file_path, counter = move_empty_rsid(prs_file_path, out_dir, replace_flag)
        if counter < 1:
            os.rename(edited_file_path, edited_file_path.replace('_edited', ''))
            print("{} doesn't require to be changed, original file with replaced header has been written to the output "
                  "folder".
                  format(prs_file_path))
        else:
            print("{} has been recorded, number of rsid replacements = {}".format(edited_file_path, counter - 1))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='Absolute path to the folder with PRS scoring file', nargs='?', required=True)
    parser.add_argument('--o', help='Absolute path to the folder with output files, if doesn\'t exist, '
                                    'will be auto created', nargs='?', required=True)
    parser.add_argument('--r', help='Replace unknown SNP with SHROM_POS', nargs='?', default='n')
    args = parser.parse_args()
    print("Options in effect:\n--i absolute path of folder with prs scoring file {}\n"
          "--o absolute path to the folder with output files {}".format(args.i, args.o))
    try:
        main(args.i, args.o, args.r)
    except BaseException as e:
        print("Unexpected error: {}".format(e))
        raise e
