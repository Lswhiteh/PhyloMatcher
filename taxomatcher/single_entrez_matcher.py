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
from Bio import Entrez


def get_ua():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-e", "--email", dest="email", required=True
    )
    ap.add_argument(
        "-s",
        "--species",
        dest="species_name",
        default="Cyclophiops_semicarinatus",
        required=False,
    )
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


def main():
    ua = get_ua()
    Entrez.email = ua.email

    sp_str = ua.species_name.replace("_", " ")

    tax_id = get_tax_id(sp_str)
    tax_info = get_tax_info(tax_id)
    namelist = get_other_names(tax_info)
    if sp_str in namelist:
        namelist.remove(sp_str)

    namelist.insert(0, sp_str)
    namelist = [i.replace(" ", "_") for i in namelist]

    eq_headers = ["Tree_Sp_Name"] + [f"Eq_{i}" for i in range(len(namelist[1:]))]
    with open(f"output/{ua.species_name}_taxomatcher_output.tsv", "w") as ofile:
        ofile.write("\t".join(eq_headers) + "\n")
        ofile.write("\t".join(namelist))


if __name__ == "__main__":
    main()
