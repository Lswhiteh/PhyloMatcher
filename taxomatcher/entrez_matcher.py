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
import os

from tqdm import tqdm
from Bio import Entrez

def get_ua():
    ap = argparse.ArgumentParser()
    ap.add_argument("-e", "--email", dest="email", required=True)
    ap.add_argument(
        "-c",
        "--csv",
        dest="input_csv",
        required=True,
        help="CSV where first column is a list of target species names to look up. Will check other columns too. Ignores final column, assumes is notes.",
    )
    return ap.parse_args()


def run_esearch(sp):
    try:
        e_handle = Entrez.esearch(db="taxonomy", term=sp, rettype="gb")
        record = Entrez.read(e_handle)["IdList"]
        e_handle.close()
        return record[0]
    except:
        raise Exception


def get_tax_id(species):
    """Searches Entrez Taxonomy and grabs the official ID"""
    try:
        sp_id = run_esearch(species)
        return sp_id
    except:
        # Check for misspellings
        if species[-2:] == "um":
            species = "us".join(species.rsplit("um", 1))
        elif species[-2:] == "us":
            species = "um".join(species.rsplit("us", 1))
        elif species[-2:] == "ea":
            species = "eus".join(species.rsplit("ea", 1))
        elif species[-2:] == "eus":
            species = "ea".join(species.rsplit("eus", 1))
        try:
            sp_id = run_esearch(species)
            return sp_id
        except:
            return None


def get_tax_info(tax_ids):
    fetch = Entrez.efetch(id=tax_ids, db="taxonomy")
    taxinfo = Entrez.read(fetch)
    return taxinfo


def get_other_names(tax_dict):
    synonyms = []
    if "OtherNames" in tax_dict:
        if "Synonym" in tax_dict["OtherNames"]:
            synonyms.extend(tax_dict["OtherNames"]["Synonym"])

        if "Name" in tax_dict["OtherNames"]:
            for name in tax_dict["OtherNames"]["Name"]:
                if name["ClassCDE"] == "authority":
                    synonyms.append(" ".join(name["DispName"].split()[:2]))

    return list(set(synonyms))


def read_csv(csvfile):
    sp_list = []
    target_list = []
    idx_list = []
    with open(csvfile, "r") as ifile:
        for idx, line in enumerate(ifile.readlines(), start=-1):
            if "Tree_Sp_Name" in line:
                continue
            else:
                _line = [i for i in line.strip().split(",")[:-1] if i]
                sp_list.extend(_line)
                target_list.append(line.strip().split(",")[0])
                idx_list.extend([idx] * len(_line))

    return sp_list, target_list, idx_list


def main(input_csv, outfile, user_email):
    Entrez.email = user_email
    os.makedirs(os.path.dirname(outfile), exist_ok=True)

    sp_list, target_list, idx_list = read_csv(input_csv)
    cleaned_sp_list = [i.replace("_", " ") for i in sp_list]

    name_res = []
    pass_list = []
    fail_list = []
    tax_ids = []
    for idx, sp_str in tqdm(
        enumerate(cleaned_sp_list), "[INFO] Fetching IDs", total=len(cleaned_sp_list)
    ):
        try:
            tax_id = get_tax_id(sp_str)
            if tax_id:
                tax_ids.append(get_tax_id(sp_str))
                pass_list.append(sp_str)
            else:
                fail_list.append(target_list[idx_list[idx]])

        except:
            fail_list.append(target_list[idx_list[idx]])

    # Filter out redundant hits

    unique_tax_ids = set(tax_ids)
    pass_list = [pass_list[i] for i in [tax_ids.index(j) for j in unique_tax_ids]]

    tax_info = get_tax_info(unique_tax_ids)

    for name, ti in tqdm(
        zip(pass_list, tax_info), "[INFO] Parsing XML data", total=len(pass_list)
    ):
        try:
            namelist = get_other_names(ti)
        except:
            fail_list.append(name)
            continue
        if name in namelist:
            namelist.remove(name)

        namelist = [i.replace(" ", "_") for i in namelist]
        namelist.append(name.replace(" ", "_"))
        name_res.append(namelist)

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

    with open(outfile, "w") as ofile:
        ofile.write("\t".join(eq_headers) + "\n")
        for names in name_res:
            ofile.write("\t".join(names) + "\n")

    with open(f"{outfile.split('.')[0]}_fails.csv", "w") as failfile:
        for i in set(fail_list):
            failfile.write(i.replace(" ", "_") + "\n")

    print("[INFO] Done")
