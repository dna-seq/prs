import gzip
import pandas as pd

gwas_statistics_file = '/path/'
rsid_path = '/path/'
result_file = '/path/'
df = pd.read_excel(rsid_path, header=1, dtype=str, sheet_name='Table Name')
N = 1216  # number of required rsid
required_rsids = df['PRS-4'].iloc[:N]  # for instance
presented_idx = list()
absent_rsid = []
absent_pos = []


def Diff(li1, li2):
    li_dif = [i for i in li1 + li2 if i not in li1 or i not in li2]
    return li_dif


def form_prs_file():
    dupl = 0
    i = 0
    unique_rsid = set()
    with gzip.open(gwas_statistics_file, 'rb') as f:
        with open(result_file, 'w+') as r:
            for line in f:
                line = line.decode("utf-8")
                if i == 0:
                    r.write(line)
                    print("write header:\n", line)
                    i += 1
                    continue
                line_splitted = line.split('\t')
                if line_splitted[0] in required_rsids.values:
                    presented_idx.append(required_rsids[required_rsids == line_splitted[0]].index[0])
                    r.write(line)
                    unique_rsid.add(line_splitted[0])
                    i += 1
                elif absent_pos and absent_rsid:
                    for pos, rs in zip(absent_pos, absent_rsid):
                        if line_splitted[0] == pos:
                            line_with_rsid = line.replace(line_splitted[0], rs)
                            r.write(line_with_rsid)
                            unique_rsid.add(rs)
                            i += 1
                diff = i - len(unique_rsid)
                if diff == 2 and dupl == 0:
                    dupl += 1
                    print("dupl {}: i={}".format(dupl, i), 'len(unique_rsid)=', len(unique_rsid), line_splitted[0])
                elif diff > 2:
                    dupl += 1
                    print("dupl {}: i={}".format(dupl, i), 'len(unique_rsid)=', len(unique_rsid), line_splitted[0])

    print("number of iterations over SNP items", i)
    print("unique iterations", len(unique_rsid))
    return presented_idx


def check_absent_items(presence_idx, required_number):
    li2 = list(range(required_number))
    li3 = Diff(presence_idx, li2)
    print("diff idx: ", li3)
    absent_rsids = list()
    for i in li3:
        absent_rsids.append(required_rsids[i])
    print("absent_rsids: ", absent_rsids)


def main():
    presented_indexes = form_prs_file()
    check_absent_items(presented_indexes, N)


main()