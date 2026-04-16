# *LoCo* - _Local Correlation Analysis_

[![Tests](https://github.com/tstohn/LoCo/actions/workflows/tests_cpp.yml/badge.svg?branch=master)](https://github.com/tstohn/LoCo/actions/workflows/tests_cpp.yml)

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
LoCo can be easily installed as an R-package or be build as a command line tool in cpp.

## 1.) R
### a) Install R-package loco

You can easily install LoCo from within R:
```R
remotes::install_github("https://github.com/tstohn/LoCo")
```

### b) Run

Once installed you can load LoCo in R and run it:
```R
library(loco)
locoResults <- run_loco("src/test/data_1.tsv", correlationCutoff= 0.5)
```

For the exact output format see below or in the manual.
After running it you can find the detected 'local' correlations in:
```r
locoResults$LaplacianScores
```

To the plot how these correlations change across the single-cell space you can plot the correlations per neighbourhood across any low dimensional representation of the data. We recommend to plot the correlations on the UMAP-space generated with the origional `RawData`. This way you can analyze the single-cell data with standard single-cell methods (clustering of data/ cell-tyep identification, etc.) and then plot in the same representation the distribution of correlations in neighbourhoods. Therefore, we simply plot neihgbourhoods at the coordinates of their anchor cells. You can do so by:

```r
# assign UMAP coordiantes to the results of loco
locoResults <- add_umap_coords(locoResults)
# plot the neighbourhoods into the UMAP of the RawData and
# color by the correlations of the correlation pair with the lowest laplacian score
# ( locoResults$Laplacian is ordered in increasing order by the laplacian score (lower score = more local correlations with strong global variation)
plot_local_correlation_map(locoResults, locoResults$Laplacian$FeaturePair[1])
```

### c) Output

LoCo returns a named list containing four distinct `data.frame` objects that desribe the neighbourhoods and the found local correlations.

### 1. `RawData`
Stores the processed input expression matrix in a **long-format** structure.
* **Cell Identification**: Includes a `cellID` column. If your input lacked IDs, LoCo automatically generates them in the format `C_<index>` (e.g., `C_0`, `C_1`).
* **Structure**: Each row represents a single cell, and columns represent the measured protein or gene features.

### 2. `LaplacianScores`
The primary results table summarizing feature pairs that show statistically significant local correlation patterns.
* **FeaturePair**: The names of the two features being compared (e.g., `feature1_feature2`).
* **LaplacianScore**: The calculated score used to rank the strength of the local relationship.
* **p_value**: A permutation-based significance value.
* **FeatureSet**: A comma-separated list of features that form larger "co-correlated" clusters.

### 3. `Correlations`
Provides all correlations between the filtered feature-pairs in all neighbourhoods.
* **CorrelationPair**: Matches the pairs found in the `LaplacianScores` table.
* **NeighbourhoodID**: The specific local group (e.g., `N_42`) where the correlation was calculated.
* **Correlation**: The actual correlation coefficient for that pair within that specific neighborhood.

### 4. `Neighbourhoods`
Describes the sampled neighbourhoods. Neighbourhoods are created by first sampling anchor cells (AnchorCellID) around which LoCo creates local neighbourhoods (AllCellIDs). The CellIDs are the same ones as used in `RawData` and NeighbourhoodIDs are the same as in `Correlations`.
* **NeighborhoodID**: Unique identifier for each local group in the form `N_<index>`.
* **AnchorCellID**: The ID of the "center" cell sampled to seed the neighborhood.
* **AllCellIDs**: A comma-separated list of all neighboring cells included in that specific local group (selected via KNN).


## 2.) CPP-tool
### a) Install cpp-tool
The C++ command-line tool depends on Boost.Program_options and zlib, which must be installed on your system before building. These dependencies are required for parsing command-line arguments and handling compressed data streams.
You can easily install boost/ zlib and build loco with following command:

```bash
  git clone https://github.com/tstohn/LoCo
  make install
  make loco
```

### b) Run

You will then find LoCo as an executable in the folder 'bin'. 
To see a description of the input parameters and how to use LoCo run 'bin/loco --help' from 'bin'. 
The only one compulsary parameter of LoCo is the input file:
  - the input file as a tsv file of features counts with cells in the rows and features in the columns (tab-seperated)
Nevertheless, it might make sense to set additional parameters like number of neighbourhoods, number of cells within a neighbourhood, etc.
For some examples you can have a look into the Makefile under 'make test' to see some examples of using loco.
The cpp package provides the same functionality as the R-package plus it can run LoCo with various granularities. Instead of one parameter for <neighborhoodSize> (-s / the number of cells within one neighbourhood) you can run LoCo with an array of <neighborhoodSize>, each of them generating one output to analyze correlations on many granularities (different scales from small to bigger neighbourhoods)

### c) Output:

LoCo will create several files that can be used to analyze/ plot local correlation patterns in the data. Those files will state neighborhood-ids (aas the id of the anchor cell) and cell-ids for cells in the neighborhoods. All indices start from zero and index the row of the origional input file.
Among those the most important ones are:
  - LoCo_correlations.tsv: The first column contains the index for the neighborhood (this is the row index if the anchor cell around which the neighborhood was build) and one column for every found correlation pair.
  - LoCo_laplacian.tsv: This file contains all the laplacian scores/ p-values for found correlations.
  - LoCo_coord.tsv: This file contains the coordinates of all the enighborhoods. The coordinates are defined by the used features and are the counts for all the features of the anchor cells that were used to construct those neighborhoods.
  - LoCo_cells.tsv: The first row contains all neighborhood indices, then all rows below this one contains the cell indices of the cells that are part of this neighborhood (including the anchor cell).
