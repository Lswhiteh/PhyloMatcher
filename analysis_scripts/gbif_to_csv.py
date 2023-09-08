import argparse
import csv
import subprocess
from tqdm import tqdm


def file_len(fname):
    p = subprocess.Popen(
        ["wc", "-l", fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])


def get_ua():
    ap = argparse.ArgumentParser(
        description="""Used to convert gbif ocurrence data to easily-parsable trait files. 
        Example usage: python gbif_to_csv.py -i ../data/test_trait.csv -o ../data/test_trait_out.csv -s species -c decimalLatitude decimalLongitude"""
    )
    ap.add_argument("-i", "--input-file", help="File to parse data from.")
    ap.add_argument("-o", "--output-file", help="File to write.")
    ap.add_argument(
        "-s",
        "--species-columns",
        help="Column(s) to use as the new index column (e.g. genera, species). Multiple names will be merged using underscores if provided.",
        nargs="+",
    )
    ap.add_argument(
        "-c",
        "--keep-columns",
        help="Column(s) to keep from the input file and output in the output file.",
        nargs="+",
    )
    return ap.parse_args()


def main():
    ua = get_ua()
    print("[INFO] Getting file length")
    flen = file_len(ua.input_file)
    print("[INFO] Running")
    with open(ua.input_file, "r") as ifile:
        dialect = csv.Sniffer().sniff(ifile.readline())
        delim = dialect.delimiter

    with open(ua.output_file, "w") as ofile:
        with open(ua.input_file, "r") as ifile:
            reader = csv.reader(ifile, delimiter=delim, quoting=csv.QUOTE_NONE)

            header = next(reader)
            tr_col = header.index("taxonRank")
            sp_idxs = [header.index(i) for i in ua.species_columns]
            c_idxs = [header.index(i) for i in ua.keep_columns]

            for line in tqdm(reader, total=flen):
                try:
                    # This is for gbif data only
                    if line[tr_col] != "SPECIES":
                        continue

                    sp_name = "_".join([line[i].replace(" ", "_") for i in sp_idxs])
                    k_cols = ",".join([line[i] for i in c_idxs])
                    if k_cols == "".join([","] * (len(c_idxs) - 1)):
                        continue

                    ofile.write(f"{sp_name},{k_cols}\n")

                except Exception as e:
                    print(f"Couldn't parse line because of {e}")


if __name__ == "__main__":
    main()
