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
import regex as re

def get_sp_id(sp):
    pattern = fr"(?e)({sp}){{e<=1}}" 

    sp_dict = species.name_backbone(sp_key)
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

sp_key = get_sp_id("coelognathus enganensis")[0] # outputs tuple - [0] is ID
print(sp_key)


