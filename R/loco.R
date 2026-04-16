#' Run LoCo analysis
#'
#' @param inFile Input file path.
#' @param del Delimiter used in input file.
#' @param col Whether feature names are in columns.
#' @param row Whether cell names are in rows.
#' @param zscore Apply z-score normalization.
#' @param thread Number of threads.
#' @param correlatedSetMode Mode for correlated sets (connected components, fully connected, connected by min x edges)
#' @param numberCorrelations Number of correlations to compute.
#' @param cellStateGeneFile Path to cell state gene file. Only these features are used to create cell neighbourhoods.
#' @param correlationStateGeneFile Path to correlation state gene file: only these features are used to calculate pairwise correlations.
#' @param numberNeighbourhoods number of neighbourhoods to create. Per default (0), LoCo creates x = (totalCellNumber/ neighborhoodSize) neighbourhoods.
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
  numberNeighbourhoods = 0,
  neighborhoodSize = 100,
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

  if (!is.logical(zscore)) {
    stop("`zscore` must be TRUE or FALSE")
  }
  if (!is.logical(col)) {
    stop("`col` must be TRUE or FALSE")
  }
  if (!is.logical(row)) {
    stop("`row` must be TRUE or FALSE")
  }
  if (!is.numeric(correlatedSetMode) || correlatedSetMode < 1) {
    stop("`correlatedSetMode` must be >= 1")
  }
  if (!is.numeric(numberCorrelations) || numberCorrelations < 0) {
    stop("`numberCorrelations` must be >= 0")
  }
  if (!is.character(cellStateGeneFile)) {
    stop("`cellStateGeneFile` must be a character string")
  }
  if (!is.character(correlationStateGeneFile)) {
    stop("`correlationStateGeneFile` must be a character string")
  }
  if (!is.numeric(numberNeighbourhoods) || numberNeighbourhoods < 0) {
    stop("`numberNeighbourhoods` must be >= 0")
  }
  if (!is.numeric(neighborhoodSize) || neighborhoodSize < 0) {
    stop("`neighborhoodSize` must be >= 0")
  }
  if (!is.numeric(neighborhoodKNN) || neighborhoodKNN < 0) {
    stop("`neighborhoodKNN` must be >= 0")
  }
  if (!is.numeric(correlationCutoff) || correlationCutoff < 0 || correlationCutoff > 1) {
    stop("`correlationCutoff` must be between 0 and 1")
  }
  if (!is.numeric(permutations) || permutations < 0) {
    stop("`permutations` must be >= 0")
  }
  if (!is.numeric(minSetSize) || minSetSize < 0) {
    stop("`minSetSize` must be >= 0")
  }
  if (!is.numeric(corrSetAbundance) || corrSetAbundance < 0) {
    stop("`corrSetAbundance` must be >= 0")
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
    as.integer(numberNeighbourhoods),
    as.integer(neighborhoodSize),
    as.integer(neighborhoodKNN),
    correlationCutoff,
    as.integer(permutations),
    as.integer(minSetSize),
    corrSetAbundance
  )

  return(res)
}

#' To plot local correlations across neighbourhoods LoCo needs a space to plot those correlations in.
#' One possibility is to create the UMAP space of the raw data. Within this UMAP space you can inspect
#' origional raw features and additionally assign UMAP coordiantes to neighbourhoods by taking the central
#' cell of all neighbourhoods (the anchor cell) and assigning the UMAP coordinates of this cell to their
#' neighbourhood. After running add_umap_coords() you can run plot_local_correlation_map() to plot local correlations
#' within this UMAP space.
#' @export
add_umap_coords <- function()
{
  #run UMAP
  feature_cols <- colnames(data[ , grepl("^[GME][0-9]+$", names(data)) ])
  Xmat <- as.matrix(data[, feature_cols])
  set.seed(7)
  emb <- uwot::umap(
    Xmat,
    metric = "euclidean",
    scale = TRUE
  )
  umap_data <- data.frame(
    UMAP1 = emb[,1],
    UMAP2 = emb[,2],
    pseudotime = data$pseudotime,
    middle_alpha = data$middle_alpha,
    branch1_alpha = data$ALPHA_BRANCH_1,
    branch2_alpha = data$ALPHA_BRANCH_2,
    branch_id = factor(data$branch_id)
  ) 
  
  #add UMAP
}


#' Plot correlation for features by neighbourhoods
#' Therefore create a plotting space first (like UMAP coordinates) to then plot neighbourhoods into this space
#' @import ggplot2
#' @param correlations the CorrelationMatrix rreturned by run_loco: the matrix storing correlations for all neighbourhoods
#' @param space a data.frame the describes that space in which to plot the local correlations. This data.frame MUST have 
#'              following columns: N_ids, x, y. The first column 'N_ids' contains all neighbourhood IDs as in 'correlations', 
#'              and the following columns describe the x and a y coordinates for plotting. These coordinates can be for example UMAP, PCA
#'              coordinates. To generate UMAP coordinates just run add_umap_coords().
#' @export
plot_local_correlation_map <- function(correlations, space) 
{


  p +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::theme_minimal()
}

#' Plot correlation between featureA and featureB as a scatter plot using all cells contained in neighbourhoods that are within given space boundaries
#' [x_min ... x_max] and [y_min ... y_max].
#' @param featureA first feature (x-coordinate) of the plotted correlation
#' @param featureB second feature (y-coordinate) of the plotted correlation
#' @param space the data.frame that describes the space that should be used to select neighbourhoods. All cells within those neighbourhoods
#'              are then used for plotting. This data.frame must be in the same format as the data.frame created by add_umap_coords():
#'              N_ids = neighbourhood IDs, x = x coordinates, y = y coordinates
#' @param x_min The minimum x-coordinate to filter neighbourhoods that are used for plotting. Cells in neighbouhoods are inlcuded if the anchor cell
#'              of this neighbourhood has an x-coordinate >= x_min.
#' @param x_max The maximum x-coordinate to filter neighbourhoods that are used for plotting. Cells in neighbouhoods are inlcuded if the anchor cell
#'              of this neighbourhood has an x-coordinate <= x_max.
#' @param y_min The minimum y-coordinate to filter neighbourhoods that are used for plotting. Cells in neighbouhoods are inlcuded if the anchor cell
#'              of this neighbourhood has an y-coordinate >= y_min.
#' @param y_max The maximum y-coordinate to filter neighbourhoods that are used for plotting. Cells in neighbouhoods are inlcuded if the anchor cell
#'              of this neighbourhood has an y-coordinate <= y_max.
#' @export
plot_cell_level_correlation <- function(featureA, featureB, space, x_min. x_max, y_min, y_max)
{

}