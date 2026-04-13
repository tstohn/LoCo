#' Run LoCo analysis
#'
#' @param inFile Input file path.
#' @param del Delimiter used in input file.
#' @param col Whether feature names are in columns.
#' @param row Whether cell anmes are in rows.
#' @param zscore Apply z-score normalization.
#' @param thread Number of threads.
#' @param correlatedSetMode Mode for correlated sets (connected components, fully connected, connected by min x edges)
#' @param numberCorrelations Number of correlations to compute.
#' @param cellStateGeneFile Path to cell state gene file. Only these features are used to create cell neighbourhoods.
#' @param correlationStateGeneFile Path to correlation state gene file: only these features are used to calculate pairwise correlations.
#' @param neighborhoodSize Size of neighbourhoods.
#' @param neighborhoodKNN Number of nearest neighbours.
#' @param correlationCutoff Correlation threshold (minimum correlation between pairs to consider this correlation)
#' @param permutations Number of permutations for p-value calculation.
#' @param minSetSize Minimum set size.
#' @param corrSetAbundance Minimum abundance threshold (percentage of neighbourhoods that must contain a correlated pair to consider this pair).
#' @return A named list containing LoCo results with four elements:
#'
#' \describe{
#'   \item{RawData}{
#'     A data.frame containing the input expression matrix in long form,
#'     with an additional column \code{cellID}. If the raw data did not contain named cellIDs the new
#'     cellIDs are enumerated in the form 'C_<index>', where the first index starts from 0.
#'     Rows correspond to cells and columns correspond to features.
#'   }
#'
#'   \item{LaplacianScores}{
#'     A data.frame summarising statistically significant feature pairs
#'     based on Laplacian scoring. Contains:
#'     \itemize{
#'       \item \code{FeaturePair}: names of feature pairs as 'feature1_feature2'
#'       \item \code{LaplacianScore}: computed laplacian score
#'       \item \code{p_value}: permutation-based p-value
#'       \item \code{FeatureSet}: comma-separated list of features forming sets of co-correlated features
#'     }
#'   }
#'
#'   \item{CorrelationMatrix}{
#'     A data.frame of correlation values per neighbourhood.
#'     The first column \code{CorrelationPair} contains feature pairs (same name as in LaplacianScore),
#'     and each additional column corresponds to a neighbourhood ID,
#'     containing correlation values for that neighbourhood.
#'   }
#'
#'   \item{Neighbourhoods}{
#'     A data.frame describing all constructed neighbourhoods:
#'     \itemize{
#'       \item \code{NeighborhoodID}: unique ID of each neighbourhood in the form 'N_<index>' where the index starts from 0.
#'       \item \code{AnchorCellID}: central cell of the neighbourhood. LoCo starts by sampling random cells (anchor cells) and builds local neighbourhoods around them.
#'       \item \code{AllCellIDs}: comma-separated list of all cells in this neighbourhood. These were the closest cells to the anchor cell.
#'     }
#'   }
#' }
#'
#' @examples
#' # Load example dataset shipped with the package
#' infile <- system.file("example", "data_1.tsv",
#'                       package = "loco")
#'
#' # Run LoCo on a small example dataset: find all local correlations with a correlation above 0.4
#' res <- run_loco(
#'   inFile = infile,
#'   correlationCutoff = 0.4)
#'
#' # Inspect results
#' head(res$LaplacianScores)
#' @details
#' The result is a list of data.frame containing: 1.) the RawData, 
#' 2.) LaplacianScores: scored correlation pairs that vary across single-cell sapce, 
#' 3.) CorrelationMatrix: all correlations in all neighbourhoods,
#' 4.) Neighbourhoods: the definition of neighbourhoods by their anchor cell and all contained cells
#' @export
run_loco <- function(
  inFile,
  del = "\t",
  col = FALSE,
  row = FALSE,
  zscore = TRUE,
  thread = 1,
  correlatedSetMode = 2,
  numberCorrelations = 0,
  cellStateGeneFile = "",
  correlationStateGeneFile = "",
  neighborhoodSize = 50,
  neighborhoodKNN = 5,
  correlationCutoff = 0.7,
  permutations = 100,
  minSetSize = 2,
  corrSetAbundance = 0.01
) {

  # ---- checks ----
  if (!file.exists(inFile)) {
    stop("Input file does not exist: ", inFile)
  }

  if (!is.character(del) || nchar(del) != 1) {
    stop("`del` must be a single character")
  }

  if (!is.numeric(thread) || thread < 1) {
    stop("`thread` must be >= 1")
  }

  # ---- call C++ ----
  res <- run_loco_cpp(
    inFile,
    del,
    col,
    row,
    zscore,
    as.integer(thread),
    as.integer(correlatedSetMode),
    as.integer(numberCorrelations),
    cellStateGeneFile,
    correlationStateGeneFile,
    neighborhoodSize,
    as.integer(neighborhoodKNN),
    correlationCutoff,
    as.integer(permutations),
    as.integer(minSetSize),
    corrSetAbundance
  )

  return(res)
}

# TODO: add umap coords to res, when plot_n is called the first time these umap coords are filled

#' Plot correlation for features by neighbourhoods
#' Therefore create UMAP space first to then plot neighbourhoods into the same space
#' @import ggplot2
#' @export
plot_neighbourhood <- function(res, featureA = "x", featureB = "y") 
{


  p +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::theme_minimal()
}