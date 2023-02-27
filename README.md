# Taxomatcher
Python module to query the NCBI Taxonomy database for synonym names of target species.

## Installation
Required packages:
- Biopython
- pygbif

To use conda:
`conda env create -f taxomatch_env.yaml`

## Input

### GBIF
CSV input file should be a single column of target species to look up (other columns will be ignored). Names can be either space or underscore ("_") separated, i.e. "Sphenodon_punctatus" is equivalent to "Sphenodon punctatus". 

E.g.
```
Sphenodon_punctatus
Gonyosoma_prasinus
```

### NCBI
Due to the specificity of Entrez results and the relative sparseness of taxonomy data a CSV intended for NCBI matching can have multiple columns, assuming the first column is the target species and the remaining columns are prior-known synonyms. All names will be flattened and searched to increase the chances of matches in the database. A single-column CSV file will also work, identically to the GBIF format.

Note: currently the multi-Entrez script assumes the final column in the CSV is for notes. If people decide this is worth fixing I will, but it seems like the GBIF approach is much better across the board.

E.g.
```
Sphenodon_punctatus,Hatteria_punctata
Gonyosoma_prasinus,Coluber_prasinus,Elaphe_prasina,Rhadinophis_prasinus,Rhadinophis_prasina
```

## Usage

### GBIF (recommended)
```$ python taxomatcher/multi_gbif_matcher.py -h
usage: multi_gbif_matcher.py [-h] -c INPUT_CSV [-t THREADS]

options:
  -h, --help            show this help message and exit
  -c INPUT_CSV, --csv INPUT_CSV
                        CSV where first column is a list of target species names to look up.
  -t THREADS, --threads THREADS
```

### NCBI (lot of misses, not recommended)
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

