# *LoCo* - _Local Correlation Analysis_

[![Tests](https://github.com/tstohn/LoCo/actions/workflows/tests_cpp.yml/badge.svg?branch=main)](https://github.com/tstohn/LoCo/actions/workflows/tests_cpp.yml)

<img src="https://github.com/tstohn/LoCo/raw/master/docs/LoCo5.png" width="200" />

LoCo (Local Correlation Analysis) detects locally structured correlation patterns from the continuous single-cell state space.
Input to LoCo is a single-cell features matrix and output are correlations between feature pairs. The reported correlations
change a lot across regions of the single-cell space (within the single-cell dataset) but vary only slighlty between close regions.
LoCo calcualtes p-values for the reported ranking scores of local correlations, reports sets of features that seemn to co-correlate 
with the correlated feature-pair and reports all grouped neighborhoods, correlations within them etc.

Rather than partitioning cells into discrete clusters, LoCo defines local neighborhoods based on similarity in the single-cell manifold. 
The neighborhood structure can be constructed using all measured
features or a biologically informed subset that defines the relevant cellular state space (cell-state features).
Within each local neighborhood, LoCo computes pairwise correlations between selected molecular
features, capturing context-specific correlation patterns. To identify correlations that reflect structured
biological variation rather than noise, we assess how these local correlations change across the space.
Specifically, LoCo create a neighborhood graph by connecting neighborhoods that are close in the
single-cell space and employs a Laplacian-based scoring approach that prioritizes correlations which
vary smoothly with respect to the neighborhood graph while exhibiting substantial variation across the
global single-cell space.

<img src="https://github.com/tstohn/LoCo/raw/master/docs/Overview.png" />


# Getting started:

LoCo is a cpp-tool that has to be compiled and can then be run from the command line in a terminal - an R-package is following soon. To compile the tool LoCo requires boost::program_options as well as nanoflann (https://github.com/jlblancoc/nanoflann), which is included as a header in this git-repo.
To download LoCo, install boost::program_options and compile the tool run the commands below:

```bash
  git clone https://github.com/tstohn/LoCo
  make install
  make loco
```

You will then find LoCo as an executable in the folder 'bin'. 
To see a description of the input parameters and how to use LoCo run 'bin/loco --help' from 'bin'. 
The only one compulsary parameter of LoCo is the input file:
  - the input file as a tsv file of features counts with cells in the rows and features in the columns (tab-seperated)
Nevertheless, it might make sense to set additional parameters like number of neighbourhoods, number of cells within a neighbourhood, etc.
For some examples you can have a look into the Makefile under 'make test' to see some examples of using loco.

# Output:

LoCo will create several files that can be used to analyze/ plot local correlation patterns in the data. Those files will state neighborhood-ids (aas the id of the anchor cell) and cell-ids for cells in the neighborhoods. All indices start from zero and index the row of the origional input file.
Among those the most important ones are:
  - LoCo_correlations.tsv: The first column contains the index for the neighborhood (this is the row index if the anchor cell around which the neighborhood was build) and one column for every found correlation pair.
  - LoCo_laplacian.tsv: This file contains all the laplacian scores/ p-values for found correlations.
  - LoCo_coord.tsv: This file contains the coordinates of all the enighborhoods. The coordinates are defined by the used features and are the counts for all the features of the anchor cells that were used to construct those neighborhoods.
  - LoCo_cells.tsv: The first row contains all neighborhood indices, then all rows below this one contains the cell indices of the cells that are part of this neighborhood (including the anchor cell).
