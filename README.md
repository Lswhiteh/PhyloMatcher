# PhyloMatcher

Python modules to query the GFIB or NCBI Taxonomy databases for synonym names of target species.

## For use as GUI 

To begin, install the executable for your operating system from the [latest release](https://github.com/Lswhiteh/PhyloMatcher/releases). Run it, and select either PhyloMatcher (recommended) or the standalone Traitmatcher.

### PhyloMatcher 
PhyloMatcher is an end-to-end tool which turns a CSV containing a column of species names into a column of rows in |canonical name, synonym 1, synonym 2, ...| TSV format. It can then optionally feed this output into Traitmatcher if given a CSV file of traits (in rows of |Species_name, trait value 1, trait value 2, ...|), standardizing all species_names. See Help for examples.
To use:

1. Select PhyloMatcher on start.
2. Click "Browse" and select a Species CSV file from your computer.
   - To see an example, click "Help" and then the 'Example Sp_CSV' tab.
3. Click the newly appeared "Browse" button and select a place on your computer for PhyloMatcher to store its output.
4. Choose between searching from the GBIF or NCBI database. GBIF tends to produce significantly higher quality results and is recommended.
   - If choosing NCBI, input an email.
   - If choosing GBIF, optionally choose a thread count.
5. Click "Run".
7. A progress bar should appear and gradually fill. After PhyloMatcher is done, a popup should confirm that it is finished.
Optional:
8. Input a trait CSV file if you wish to run Traitmatcher
   - To see an example, click Help and then go the 'Example Trait_CSV' tab
9. Click run. A popup should appear confirming this.


### Standalone Traitmatcher 
Standalone Traitmatcher allows you to use just the Traitmatcher functionality, but will require you to provide a TSV file of synonyms in the same format as PhyloMatcher's output - e.i. canonical names first (case-sensitive, and using underscore seperation (ex: "Homo_sapiens")), followed by synonyms. See Help for examples.

To use: 
1. Select Standalone Traitmatcher on start
2. Select a trait CSV file from your computer using the first browse button. Specify whether or not it has a header (as is present in the example)
   - To see an example, click Help and then go the 'Example Trait_CSV' tab
3. Select a species TSV file from your computer using the second browse button
    - To see an example, click Help and then go the 'Example TSV for standalone Traitmatcher' tab
5. Select a place for output to be stored with the third button
6. Once these have been selected, click the run button that appears.


## Command Line Installation

Required packages:

- Biopython
- pygbif
- pandas
- tqdm

Easiest installation using a conda environment and pip:

```{bash}
conda create -n tm-env -c conda-forge python
conda activate tm-env
pip install phylomatcher
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
$ phylomatcher gbif -h
usage: phylomatcher gbif [-h] -i INPUT_CSV -o OUTFILE [-t THREADS]

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
$ phylomatcher ncbi -h
usage: phylomatcher ncbi [-h] -i INPUT_CSV -o OUTFILE -e EMAIL

options:
  -h, --help            show this help message and exit
  -i INPUT_CSV, --csv INPUT_CSV
                        CSV where first column is a list of target species names to look up.
  -o OUTFILE, --outfile OUTFILE
                        Path to output.
  -e EMAIL, --email EMAIL
```

### Matching trait files to phylomatcher output

Once you have found synonyms from tree output you can match that to any phenotypic/trait data you have. Similar to the other inputs the first column must contain the species names you're targeting, and they should be formatted identically to how the phylomatcher output looks (case-sensitive and separated by underscores.)

Usage:

```{bash}
$ phylomatcher trait -h
usage: phylomatcher trait [-h] -t TRAITFILE -s SPECIESFILE -o OUTFILE

options:
  -h, --help            show this help message and exit
  -t TRAITFILE, --traitfile TRAITFILE
                        CSV of trait values, first column must be species names.
  -s SPECIESFILE, --speciesfile SPECIESFILE
                        CSV of species synonyms output by the gbif or ncbi modules.
  -o OUTFILE, --outfile OUTFILE
                        Path to output.
```
