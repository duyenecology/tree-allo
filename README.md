# Tree allometry

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


Code repostitory to run the analysis and generate the manuscript for Nguyen and Katabuchi "Saturating allometric relationships reveal how wood density shapes global tree architecture"

## Reproduce the results

Codes (R and STAN) and workflow are managed with the R package `targets` (https://github.com/ropensci/targets).

### Running code on local

To run analysis:

```bash
# To install R packages for the first run
# Rscript -e "renv::restore()"
Rscript run.R
```

To generate the manuscript:

```bash
make
```

Requirements:

- cmdstan 2.34.1
- quarto
- latexdiff
- R (4.3.2)
	- renv (`renv::restore()` will install all the R packages)

### Running code in Apptainer (Linux)

To build the Apptainer container:

```bash
sudo apptainer build radian.sif radian.def
```

or you can download the container from XXX.

To run analysis or to render the manuscript:

```bash
# To install R packages for the first run
# apptainer exec --env RENV_PATHS_CACHE=/home/${USER}/renv \
		# --env RENV_PATHS_PREFIX_AUTO=TRUE \
 		# radian.sif Rscript run.R -e "renv::restore()"

> ./run.sh
1) tar_make() on local
2) tar_make() on Apptainer
3) Enter in the Apptainer container
4) make (render manuscript) on the Apptainer container
Enter number：
```

Requirements:

- Apptainer (or singularity)
- `renv` directory on your home
