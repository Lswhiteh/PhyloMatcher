# Taxomatcher
Python module to query the NCBI Taxonomy database for synonym names of target species.

## Installation
Required packages:
- Biopython

To use conda:
`conda env create -f taxomatch_env.yaml`

## Usage
Single lookup:
`python taxomatcher/single_entrez_matcher.py -e <email@email.com> -s <Species_Name>`

CSV lookup:
`python taxomatcher/multi_entrez_matcher.py -e <email@email.com> -c data/taxolist.csv`

Currently missing a lot of hits on the multi-matcher. Work in progress.
