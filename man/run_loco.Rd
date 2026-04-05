\name{run_loco}
\alias{run_loco}
\title{Run LoCo analysis on single-cell data}

\description{
Runs the LoCo algorithm to detect local correlations between features (e.g., genes)
across neighborhoods in single-cell data. The method identifies feature pairs whose
correlations vary across the cellular state space and reports correlation and slope
Laplacians, neighborhood-level statistics, and correlated feature sets.
}

\usage{
run_loco(
  inFile,
  outFile,
  prefix = "LoCo",
  del = "\t",
  col = FALSE,
  row = FALSE,
  zscore = TRUE,
  thread = 1,
  correlatedSetMode = 2,
  numberCorrelations = 0,
  cellStateGeneFile = "",
  correlationStateGeneFile = "",
  numNeighborhoods = 0,
  neighborhoodSizeString = "50",
  neighborhoodKNN = 5,
  correlationCutoff = 0.7,
  permutations = 100,
  minSetSize = 2,
  corrSetAbundance = 0.01
)
}

\arguments{
  \item{inFile}{Path to input file (TSV format), containing a cell-by-feature matrix.}

  \item{outFile}{Output directory where results will be written.}

  \item{prefix}{Prefix for output files.}

  \item{del}{Delimiter used in the input file (default: tab).}

  \item{col}{Logical; whether to explicitly parse column names.}

  \item{row}{Logical; whether to explicitly parse row names.}

  \item{zscore}{Logical; whether to z-score normalize the data (default: TRUE).}

  \item{thread}{Integer; number of threads to use.}

  \item{correlatedSetMode}{Integer; method for detecting correlated feature sets:
    \describe{
      \item{0}{Fully connected components (cliques)}
      \item{1}{Connected components}
      \item{2}{Default mode (recommended)}
    }
  }

  \item{numberCorrelations}{Integer; number of top correlations/slopes to retain.
  If greater than zero, only feature pairs with lowest Laplacian scores are written.}

  \item{cellStateGeneFile}{Optional file containing a comma-separated list of genes
  used to define the state space (neighborhood construction).}

  \item{correlationStateGeneFile}{Optional file containing a comma-separated list of genes
  used for correlation analysis.}

  \item{numNeighborhoods}{Integer; number of neighborhoods. If 0, it is automatically
  set based on the number of cells.}

  \item{neighborhoodSizeString}{Character; size of neighborhoods. Can be a single value
  (e.g., "50") or a vector-like string (e.g., "[10,50,100]") for multi-scale analysis.}

  \item{neighborhoodKNN}{Integer; number of nearest neighbors used to connect neighborhoods.}

  \item{correlationCutoff}{Numeric; minimum correlation threshold for feature connections.}

  \item{permutations}{Integer; number of permutations used to compute p-values
  for Laplacian scores.}

  \item{minSetSize}{Integer; minimum size of correlated feature sets to consider.}

  \item{corrSetAbundance}{Numeric; minimum fraction of neighborhoods in which
  a correlation must be observed to be retained.}
}

\details{
LoCo constructs a graph of neighborhoods in single-cell data and computes correlations
between features within each neighborhood. It identifies feature pairs whose correlations
change across the state space using Laplacian-based scoring.

The method first detects sets of correlated features based on a graph representation,
then evaluates pairwise relationships within these sets. Correlations that persist
across multiple connected neighborhoods are prioritized.

Outputs include:
\itemize{
  \item Laplacian scores for correlation and slope per feature pair
  \item Neighborhood-level correlation and slope matrices
  \item Coordinates of neighborhoods
  \item Filtered correlated feature sets
}
}

\value{
Called for its side effects. Writes multiple output files to \code{outFile}.
}

\examples{
\dontrun{
run_loco(
  inFile = "data.tsv",
  outFile = "results",
  prefix = "example",
  numNeighborhoods = 20,
  correlationCutoff = 0.7
)
}
}

\author{
Tim Stohn
}