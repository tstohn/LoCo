#' Run LoCo analysis
#'
#' @export
run_loco <- function(
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
) {

  .Call(
    "_loco_run_loco",
    inFile,
    outFile,
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
}