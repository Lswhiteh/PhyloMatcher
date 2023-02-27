#!/usr/bin/env python

"""
A utility to match NCBI Taxonomy synonym names given a taxon name.

Resources:
- https://harryincupboard.blog/3
- https://biopython.org/docs/1.76/api/Bio.Entrez.html
- https://www.ncbi.nlm.nih.gov/books/NBK25497/
- https://www.ncbi.nlm.nih.gov/books/NBK25500/

Logan Whitehouse - lswhiteh@unc.edu
"""
import argparse
import multiprocessing as mp
from tqdm import tqdm
from pygbif import species


def get_ua():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-c",
        "--csv",
        dest="input_csv",
        required=True,
        help="CSV where first column is a list of target species names to look up. ",
    )
    ap.add_argument("-t", "--threads", dest="threads", required=False, default=4)
    return ap.parse_args()


def read_csv(csvfile):
    target_list = []
    with open(csvfile, "r") as ifile:
        for line in ifile.readlines():
            if "Tree_Sp_Name" in line:
                continue
            else:
                target_list.append(line.strip().split(",")[0])

    return target_list


def get_sp_id(sp):
    sp_dict = species.name_backbone(sp)
    if "speciesKey" in sp_dict:
        return sp_dict["speciesKey"], sp_dict["species"]
    else:
        return None, None


def get_synonyms(sp_key):
    nu = species.name_usage(sp_key, data="synonyms")["results"]
    if len(nu) > 0:
        return [i["canonicalName"] for i in nu if "canonicalName" in i]
    else:
        return []


def worker(sp):
    key, curr_name = get_sp_id(sp)
    if key:
        synonyms = get_synonyms(key)
        synonyms.insert(0, sp)
        synonyms.append(curr_name)
    else:
        synonyms = [sp]

    return synonyms


def main():
    ua = get_ua()

    sp_list = read_csv(ua.input_csv)
    cleaned_sp_list = [i.replace("_", " ") for i in sp_list]

    with mp.Pool(ua.threads) as p:
        synonyms = list(
            tqdm(
                p.imap(worker, cleaned_sp_list, chunksize=4),
                desc="[INFO] Fetching GBIF information",
                total=len(cleaned_sp_list),
            )
        )

    max_len = max([len(i) for i in synonyms])
    eq_headers = (
        ["Tree_Sp_Name"] + [f"Eq_{i}" for i in range(max_len - 1)] + ["Curr_Name"]
    )

    with open(f"output/gbif_output.tsv", "w") as ofile:
        ofile.write("\t".join(eq_headers) + "\n")
        for names in synonyms:
            ofile.write("\t".join([i.replace(" ", "_") for i in names]) + "\n")


if __name__ == "__main__":
    main()
