# Taxomatcher

Python modules to query the GFIB or NCBI Taxonomy databases for synonym names of target species.

## Installation

Required packages:

- Biopython
- pygbif
- pandas
- tqdm

Easiest installation using a conda environment and pip:

```{bash}
conda create -n tm-env -c conda-forge python
conda activate tm-env
pip install taxomatcher
```

---

## Modes

### GBIF

CSV input file should be a single column of target species to look up (other columns will be ignored). Names can be either space or underscore ("_") separated, i.e. "Sphenodon_punctatus" is equivalent to "Sphenodon punctatus". 

E.g.
```
Sphenodon_punctatus
Gonyosoma_prasinus
```

Usage:

```{bash}
$ taxomatcher gbif -h
usage: taxomatcher gbif [-h] -i INPUT_CSV -o OUTFILE [-t THREADS]

options:
  -h, --help            show this help message and exit
  -i INPUT_CSV, --csv INPUT_CSV
                        CSV where first column is a list of target species names to look up.
  -o OUTFILE, --outfile OUTFILE
                        Path to output.
  -t THREADS, --threads THREADS
```

### NCBI

Due to the specificity of Entrez results and the relative sparseness of taxonomy data a CSV intended for NCBI matching can have multiple columns, assuming the first column is the target species and the remaining columns are prior-known synonyms. All names will be flattened and searched to increase the chances of matches in the database. A single-column CSV file will also work, identically to the GBIF format.

Note: currently the multi-Entrez script assumes the final column in the CSV is for notes. If people decide this is worth fixing I will, but it seems like the GBIF approach is much better across the board.

E.g.

```{text}
Sphenodon_punctatus,Hatteria_punctata
Gonyosoma_prasinus,Coluber_prasinus,Elaphe_prasina,Rhadinophis_prasinus,Rhadinophis_prasina
```

Usage:

```{bash}
$ taxomatcher ncbi -h
usage: taxomatcher ncbi [-h] -i INPUT_CSV -o OUTFILE -e EMAIL

options:
  -h, --help            show this help message and exit
  -i INPUT_CSV, --csv INPUT_CSV
                        CSV where first column is a list of target species names to look up.
  -o OUTFILE, --outfile OUTFILE
                        Path to output.
  -e EMAIL, --email EMAIL
```

### Matching trait files to taxomatcher output

Once you have found synonyms from tree output you can match that to any phenotypic/trait data you have. Similar to the other inputs the first column must contain the species names you're targeting, and they should be formatted identically to how the taxomatcher output looks (case-sensitive and separated by underscores.)

Usage:

```{bash}
$ taxomatcher trait -h
usage: taxomatcher trait [-h] -t TRAITFILE -s SPECIESFILE -o OUTFILE

options:
  -h, --help            show this help message and exit
  -t TRAITFILE, --traitfile TRAITFILE
                        CSV of trait values, first column must be species names.
  -s SPECIESFILE, --speciesfile SPECIESFILE
                        CSV of species synonyms output by the gbif or ncbi modules.
  -o OUTFILE, --outfile OUTFILE
                        Path to output.
```