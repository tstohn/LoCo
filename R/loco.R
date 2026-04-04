#' Run LoCo analysis
#'
#' @export
run_loco <- function(
  inFile,
  outFile = NULL,
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

  if (!is.null(outFile)) {
    warning("`outFile` is deprecated and ignored.")
  }

  # ---- call C++ ----
  res <- .Call(
    "_loco_run_loco",
    inFile,
    "",  # no longer used
    prefix,
    del,
    col,
    row,
    zscore,
    as.integer(thread),
    as.integer(correlatedSetMode),
    as.integer(numberCorrelations),
    cellStateGeneFile,
    correlationStateGeneFile,
    as.integer(numNeighborhoods),
    neighborhoodSizeString,
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