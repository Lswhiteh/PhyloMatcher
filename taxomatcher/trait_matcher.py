import pandas as pd
import csv
from tqdm import tqdm
import os
import subprocess


def file_len(fname):
    p = subprocess.Popen(
        ["wc", "-l", fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])


def main(traitfile, speciesfile, outfile, header):
    if "/" in outfile:
        os.makedirs(os.path.dirname(outfile), exist_ok=True)

    with open(speciesfile, "r") as r_obj:
        csv_reader = csv.reader(r_obj, delimiter="\t")
        spec_list = list(csv_reader)

    with open(traitfile, "r") as ifile:
        dialect = csv.Sniffer().sniff(ifile.readline())
        delim = dialect.delimiter

    print("[INFO] Getting length of file")
    f_len = file_len(traitfile)
    print("[INFO] Running")

    with open(outfile, "w") as ofile:
        with open(traitfile, "r") as tfile:
            reader = csv.reader(tfile, delimiter=delim, quoting=csv.QUOTE_NONE)
            if header:
                header = next(reader)
                ofile.write(",".join(header))

            for line in tqdm(
                reader,
                total=f_len,
                desc="Searching for matches",
            ):
                t_spec = line[0]
                for s in spec_list:
                    if t_spec in s:
                        line[0] = s[0]
                ofile.write(",".join(line) + "\n")
