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
import os

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


def worker(sp, global_synonyms_dict):
    key, curr_name = get_sp_id(sp)
    if key:
        synonyms = get_synonyms(key)
        synonyms.insert(0, sp)
        synonyms.append(curr_name)
    else:
        synonyms = [sp]

    global_synonyms_dict[sp] = synonyms


def main(input_csv, outfile, threads):
    sp_list = read_csv(input_csv)
    cleaned_sp_list = [i.replace("_", " ") for i in sp_list]

    synonyms_dict = {}

    # List of tuples (species, synonyms dict), allows us to pass synonyms dict into worker
    species_args = [(species, synonyms_dict) for species in cleaned_sp_list]

    # Create shared dictionary
    manager = mp.Manager()
    synonyms_dict = manager.dict()

    # Use starmap to call worker function with multiple arguments
    with mp.Pool(threads) as p:
        p.starmap(worker, [(sp, synonyms_dict) for sp in cleaned_sp_list])

    # Convert shared dictionary to regular dictionary
    synonyms_dict = dict(synonyms_dict)

    # Convert shared dictionary to regular dictionary
    synonyms_dict = dict(synonyms_dict)

    print(synonyms_dict)

    # max_len = max([len(i) for i in synonyms])
    # eq_headers = (
    #     ["Tree_Sp_Name"] + [f"Eq_{i}" for i in range(max_len - 1)] + ["Curr_Name"]
    # )

    # os.makedirs("../output", exist_ok=True)
    # filename = f"../output/{run_name}_gbif_output.tsv"

    # if not os.path.isfile(filename):
    #     # if output file does not exist, create an empty file
    #     open(filename, 'a').close()
    max_len = max([len(i) for i in synonyms])
    eq_headers = (
        ["Tree_Sp_Name"] + [f"Eq_{i}" for i in range(max_len - 1)] + ["Curr_Name"]
    )
    
    os.makedirs(os.path.dirname(outfile), exist_ok=True)

    # with open(filename, "w") as ofile:
    #     ofile.write("\t".join(eq_headers) + "\n")
    #     for names in synonyms:
    #         ofile.write("\t".join([i.replace(" ", "_") for i in names]) + "\n")

main("../data/3_species.csv", "..", 4)