# Taxomatcher
Python module to query the NCBI Taxonomy database for synonym names of target species.

## Installation
Required packages:
- Biopython

To use conda:
`conda env create -f taxomatch_env.yaml`

## Usage
Single lookup:
```$ python taxomatcher/single_entrez_matcher.py -h
usage: single_entrez_matcher.py [-h] -e EMAIL [-s SPECIES_NAME]

options:
  -h, --help            show this help message and exit
  -e EMAIL, --email EMAIL
  -s SPECIES_NAME, --species SPECIES_NAME
```

CSV lookup:
```$ python taxomatcher/multi_entrez_matcher.py -h
usage: multi_entrez_matcher.py [-h] -e EMAIL -c INPUT_CSV

options:
  -h, --help            show this help message and exit
  -e EMAIL, --email EMAIL
  -c INPUT_CSV, --csv INPUT_CSV
                        CSV where first column is a list of target species names to look up. Will check other columns too. Ignores final column,
                        assumes is notes.
```
