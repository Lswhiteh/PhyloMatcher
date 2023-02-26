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

from Bio import Entrez


def get_ua():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-e", "--email", dest="email", required=True
    )
    ap.add_argument(
        "-c",
        "--csv",
        dest="input_csv",
        required=True,
        help="CSV where first column is a list of target species names to look up. Will check other columns too. Ignores final column, assumes is notes.",
    )
    ap.add_argument("-t", "--threads", dest="threads", default=4, required=False)
    return ap.parse_args()


def get_tax_id(species):
    """Searches Entrez Taxonomy and grabs the official ID"""
    e_handle = Entrez.esearch(db="taxonomy", term=species, rettype="gb")
    record = Entrez.read(e_handle)["IdList"]
    e_handle.close()

    return record[0]


def get_tax_info(tax_id):
    search = Entrez.efetch(id=tax_id, db="taxonomy")
    taxinfo = Entrez.read(search)
    return taxinfo[0]


def get_other_names(tax_dict):
    synonyms = tax_dict["OtherNames"]["Synonym"]
    for name in tax_dict["OtherNames"]["Name"]:
        if name["ClassCDE"] == "authority":
            synonyms.append(" ".join(name["DispName"].split()[:2]))

    return list(set(synonyms))


def read_csv(csvfile):
    sp_list = []
    target_list = []
    with open(csvfile, "r") as ifile:
        for line in ifile.readlines():
            if "Tree_Sp_Name" in line:
                continue
            else:
                sp_list.extend([i for i in line.strip().split(",")[:-1] if i])
                target_list.append(line.strip().split(",")[0])

    return sp_list, target_list


def worker(sp_str):
    try:
        tax_id = get_tax_id(sp_str)
        tax_info = get_tax_info(tax_id)
        namelist = get_other_names(tax_info)
        if sp_str in namelist:
            namelist.remove(sp_str)

        namelist = [i.replace(" ", "_") for i in namelist]
        namelist.append(sp_str.replace(" ", "_"))

        return sorted(namelist)

    except Exception as e:
        return ("Failed", sp_str)


def main():
    ua = get_ua()
    Entrez.email = ua.email

    sp_list, target_list = read_csv(ua.input_csv)
    cleaned_sp_list = [i.replace("_", " ") for i in sp_list]

    print("[INFO] Searching NCBI")
    pool = mp.Pool(8)
    work_res = pool.imap(worker, cleaned_sp_list, chunksize=4)
    pool.close()

    fail_list = []
    name_res = []

    for res in work_res:
        if isinstance(res, tuple):
            fail_list.append(res[-1])
        else:
            name_res.append(res)

    # Remove redundant lists to allow for checking all csv entries
    name_res = set([tuple(i) for i in name_res])
    name_res = [list(i) for i in name_res]

    for n in name_res:
        for s in target_list:
            if s in n:
                n.remove(s)
                n.insert(0, s)

    max_len = max([len(i) for i in name_res])
    eq_headers = ["Tree_Sp_Name"] + [f"Eq_{i}" for i in range(max_len - 1)]

    with open(f"output/multi_taxomatcher_output.tsv", "w") as ofile:
        ofile.write("\t".join(eq_headers) + "\n")
        for names in name_res:
            ofile.write("\t".join(names) + "\n")

    with open("output/multi_taxomatcher_fails.tsv", "w") as failfile:
        for i in fail_list:
            failfile.write(i.replace(" ", "_") + "\n")

    print("[INFO] Done")


if __name__ == "__main__":
    main()
