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


def get_edited_file_path(output_dir, file_path):
    path_splitted = file_path.split('/')
    name_pgs = path_splitted[-1].split('.')
    cleaned_file_name = '{}_edited.{}'.format(name_pgs[0], name_pgs[1])
    cleaned_file_path = os.path.join(output_dir, cleaned_file_name)
    return cleaned_file_path


def get_rsid_idx(split_line):
    def split_line_casefold(x):
        return [xx.casefold() for xx in x]
    if 'rsid' in split_line_casefold(split_line):
        index = split_line_casefold(split_line).index('rsid')
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


def move_empty_rsid(prs_file_path, output_dir):
    sniffer = csv.Sniffer()
    edited_file_path = get_edited_file_path(output_dir, prs_file_path)
    print("edited_file_path", edited_file_path)
    if os.path.isfile(edited_file_path):
        raise IOError("Edited file path already exists, remove or rename it first")
    i = 0
    with open(prs_file_path, 'r') as f:
        with open(edited_file_path, 'a') as output:
            for line in f:
                if not line.startswith('#') and not i:
                    i += 1
                    dialect = sniffer.sniff(line)
                    s = line.split(dialect.delimiter)
                    new_header = change_header(line)
                    output.write(new_header)
                    rsid_idx = get_rsid_idx(s)
                    continue
                if not line.startswith('#') and i > 0:
                    s = line.split(dialect.delimiter)
                    if rsid_idx is not None:
                        if not s[rsid_idx]:
                            i += 1
                            s[rsid_idx] = 'unknown_{}'.format(i)
                            print('empty rsid has been changed to \'unknown_{}\' name'.format(i))
                            edited_line = dialect.delimiter.join(s)
                            output.write(edited_line)
                        else:
                            output.write(line)
    return edited_file_path, i


def main(prs_dir, out_dir):
    for prs_file_path in parse_dir(prs_dir):
        edited_file_path, i = move_empty_rsid(prs_file_path, out_dir)
        if not i:
            os.rename(edited_file_path, edited_file_path.replace('_edited', ''))
            print("{} doesn't require to be changed, original file has been written to the output folder".
                  format(prs_file_path))
        else:
            print("{} has been recorded".format(edited_file_path))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', help='Absolute path to the folder with PRS scoring file', nargs='?', required=True)
    parser.add_argument('--o', help='Absolute path to the folder with output files, if doesn\'t exist, '
                                    'will be auto created', nargs='?', required=True)
    args = parser.parse_args()
    print("Options in effect:\n--i absolute path of folder with prs scoring file {}\n"
          "--o absolute path to the folder with output files {}".format(args.i, args.o))
    try:
        main(args.i, args.o)
    except BaseException as e:
        print("Unexpected error: {}".format(e))
        raise e
