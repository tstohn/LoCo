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
#' @param numberNeighbourhoods number of neighbourhoods to create. Per default (0), LoCo creates x = (totalCellNumber/ neighbourhoodSize) neighbourhoods.
#' @param neighbourhoodSize Size of neighbourhoods.
#' @param neighbourhoodKNN Number of nearest neighbours.
#' @param correlationCutoff Correlation threshold (minimum correlation between pairs to consider this correlation)
#' @param permutations Number of permutations for p-value calculation.
#' @param minSetSize Minimum set size.
#' @param corrSetAbundance Minimum abundance threshold (percentage of neighbourhoods that must contain a correlated pair to consider this pair).
#' @param correlationType spearman or pearson for the correlations to calculate. 
#' @return A named list containing LoCo results with four elements:
#'
#' \describe{
#'   \item{RawData}{
#'     A data.frame containing the input expression matrix in long form,
#'     with an additional column \code{cellID}. If the raw data did not contain named cellIDs the new
#'     cellIDs are enumerated in the form \code{C_<index>}, where the first index starts from 0.
#'     Rows correspond to cells and columns correspond to features.
#'   }
#'
#'   \item{LaplacianScores}{
#'     A data.frame summarising statistically significant feature pairs
#'     based on Laplacian scoring. Contains:
#'     \itemize{
#'       \item \code{FeaturePair}: names of feature pairs as \code{feature1_feature2}
#'       \item \code{LaplacianScore}: computed laplacian score
#'       \item \code{p_value}: permutation-based p-value
#'       \item \code{FeatureSet}: comma-separated list of features forming sets of co-correlated features
#'     }
#'   }
#'
#'   \item{Correlations}{
#'     A data.frame of correlation values per neighbourhood with 3 columns "CorrelationPair", "NeighbourhoodID", "Correlation".
#'     The first column \code{CorrelationPair} contains feature pairs (same name as in LaplacianScore),
#'     the second column \code{NeighbourhoodID} contains all neighbourhoodIDs (like \code{N_<number>}),
#'     and the third column \code{Correlation} contains the correlation of the pair in this neighbourhood.
#'   }
#'
#'   \item{Neighbourhoods}{
#'     A data.frame describing all constructed neighbourhoods:
#'     \itemize{
#'       \item \code{NeighborhoodID}: unique ID of each neighbourhood in the form \code{N_<index>} where the index starts from 0.
#'       \item \code{AnchorCellID}: central cell of the neighbourhood. LoCo starts by sampling random cells (anchor cells) and builds local neighbourhoods around them.
#'       \item \code{AllCellIDs}: comma-separated list of all cells in this neighbourhood. These were the closest cells to the anchor cell.
#'     }
#'   }
#' }
#'
#' @examples
#' # Load example dataset shipped with the package
#' infile <- system.file("example", "data_1.tsv", package = "loco")
#'
#' # Run LoCo on a small example dataset: find all local correlations with a correlation above 0.4
#' L1 <- run_loco(infile, correlationCutoff= 0.4, neighbourhoodSize = 25)
#' 
#' # Inspect the correlations with the lowest laplacian score / p-value - 
#' # these correlations change the most across the whole single-cell dataset but change only very little in their local neighbourhood
#' head(L1$LaplacianScores)
#' L2 <- add_umap_coords(L1)
#' plot <- plot_local_correlation_map(L2, L2$LaplacianScores$FeaturePair[1])
#' ggplot2::ggsave(filename = "local_correlation_map.png",
#' plot = plot,
#' width = 8,
#' height = 6,
#' dpi = 300)
#' 
#' @details
#' The result is a list of data.frame containing: 1.) the RawData, 
#' 2.) LaplacianScores: scored correlation pairs that vary across single-cell space, 
#' 3.) Correlations: all correlations in all neighbourhoods,
#' 4.) Neighbourhoods: the definition of neighbourhoods by their anchor cell and all contained cells.
#' @export
run_loco <- function(
  inFile,
  del = "\t",
  col = FALSE,
  row = FALSE,
  zscore = TRUE,
  thread = 1,
  correlatedSetMode = 1,
  numberCorrelations = 0,
  cellStateGeneFile = "",
  correlationStateGeneFile = "",
  numberNeighbourhoods = 0,
  neighbourhoodSize = 100,
  neighbourhoodKNN = 5,
  correlationCutoff = 0.5,
  permutations = 100,
  minSetSize = 2,
  corrSetAbundance = 0.01,
  correlationType = "spearman"
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
  if (!is.numeric(correlatedSetMode) || correlatedSetMode < 0) {
    stop("`correlatedSetMode` must be >= 0")
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
  if (!is.numeric(neighbourhoodSize) || neighbourhoodSize < 0) {
    stop("`neighbourhoodSize` must be >= 0")
  }
  if (!is.numeric(neighbourhoodKNN) || neighbourhoodKNN < 0) {
    stop("`neighbourhoodKNN` must be >= 0")
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
  if (!is.numeric(corrSetAbundance) || corrSetAbundance < 0 || corrSetAbundance > 1) {
    stop("`corrSetAbundance` must be >= 0 and <= 1")
  }

  if( (correlationType != "spearman") && (correlationType != "pearson") )
  {
        stop("`correlationType` must be 'spearman' or 'pearson'.")
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
    as.integer(neighbourhoodSize),
    as.integer(neighbourhoodKNN),
    correlationCutoff,
    as.integer(permutations),
    as.integer(minSetSize),
    corrSetAbundance,
    correlationType
  )

  return(res)
}


#' Add UMAP coordinates to LoCo result
#'
#' @description
#' Creates UMAP space for raw data and adds the UMAP-coordinates UMAP1/UMAP2 to RawData, Neighbourhoods and Correlations.
#' Run: \code{umapAddedLocoResult <- add_umap_coords(locoResult)} to add the UMAP coords to the loco result.
#' \code{umapAddedLocoResult} can be used to plot the correlations in neighbourhoods on a UMAP.
#' 
#' @details
#' To plot local correlations across neighbourhoods LoCo needs a space to plot those correlations in.
#' One possibility is to create the UMAP space of the raw data. Within this UMAP space you can inspect
#' original raw features and additionally assign UMAP coordinates to neighbourhoods by taking the central
#' cell of all neighbourhoods (the anchor cell) and assigning the UMAP coordinates of this cell to their
#' neighbourhood. After running \code{add_umap_coords()} you can run \code{plot_local_correlation_map()} to plot local correlations
#' within this UMAP space. You do not have to use these UMAP coordinates but this function makes it easy
#' to create the UMAP space for the raw data and the neighbourhoods. However, you can create your own
#' embedding (PCA, ...) and plot neighbourhoods (by assigning neighbourhood coordinates to the coordinates of their anchor cells).
#' \code{add_umap_coords()} adds the UMAP coordinates to the RawData data.frame, the Neighbourhoods data.frame and also to the Correlations data.frame.
#' 
#' @param locoResult the result of run_loco
#' @param n_pcs the number of components for PCA (default 0 = skip PCA embedding). If n_pcs > 0 the function firstly decomposes the input matrix into PCs and runs umap on the PCs. In this case scale is set to FALSE, as the PCA already scales the data.
#' @param scale z-score normalize feature counts before calculating UMAP coords (TRUE). This is ONLY considered when n_pcs is 0. If n_pcs > 0 the functions uses the first \code{n_pcs} PCs to run the UMAP embedding and does not z-score normalize in between!
#' @param metric distance-metric for UMAP (euclidean)
#' @param seed seed for UMAP (7)
#' @export
add_umap_coords <- function(locoResult,
                            n_pcs = 0,
                            scale = TRUE,
                            metric = "euclidean",
                            seed = 7)
{
  set.seed(seed)

  # -----------------------------
  # 1. sanity checks
  # -----------------------------
  if (!is.list(locoResult)) {
    stop("locoResult must be a list: the result of run_loco().")
  }

  if (is.null(locoResult$RawData)) {
    stop("locoResult has no RawData element. Please run run_loco and use add_umap_coords with this result.")
  }
  if (!is.data.frame(locoResult$RawData)) {
    stop("RawData (part of the result of run_loco) must be a data.frame or tibble. Please run run_loco and use add_umap_coords with this result.")
  }
  if (is.null(locoResult$Correlations)) {
    stop("locoResult has no Correlations element. Please run run_loco and use add_umap_coords with this result.")
  }
  if (!is.data.frame(locoResult$Correlations)) {
    stop("Correlations (part of the result of run_loco) must be a data.frame or tibble. Please run run_loco and use add_umap_coords with this result.")
  }

  rawData <- locoResult$RawData
  if (!"cellID" %in% colnames(rawData)) 
  {
    stop("RawData must contain a 'cellID' column. Please run run_loco and use add_umap_coords with this result.")
  }

  # -----------------------------
  # 2. split metadata vs matrix
  # -----------------------------

  cell_ids <- rawData$cellID
  feature_cols <- setdiff(colnames(rawData), "cellID")
  if (length(feature_cols) < 2) 
  {
    stop("Not enough feature columns for UMAP. There is only a single feature.")
  }
  X_mat <- as.matrix(rawData[, feature_cols])
  # ensure numeric
  mode(X_mat) <- "numeric"

  # -----------------------------
  # 3. PCA (optional)
  # -----------------------------
  if (n_pcs > 0) 
  {

    # center + scale is important
    pca <- prcomp(X_mat, center = TRUE, scale. = scale)

    n_use <- min(n_pcs, ncol(pca$x))

    emb_input <- pca$x[, 1:n_use, drop = FALSE]
    scale <- FALSE

  } else {
    emb_input <- X_mat
  }

  # -----------------------------
  # 3. run UMAP
  # -----------------------------

  # run UMAP on PCA data
  emb <- uwot::umap(
    emb_input,
    metric = metric,
    scale = scale
  )
  

  # -----------------------------
  # 4. add UMAP coords to rawData
  # -----------------------------
  rawData$UMAP1 <- emb[, 1]
  rawData$UMAP2 <- emb[, 2]
  locoResult$RawData <- rawData

  # -----------------------------
  # 4. add UMAP coords to Neighbourhoods
  # -----------------------------
  if ("Neighbourhoods" %in% names(locoResult)) 
  {

    neigh <- locoResult$Neighbourhoods

    if (!("AnchorCellID" %in% colnames(neigh))) 
    {
      stop("Neighbourhoods must contain 'AnchorCellID'")
    }

    # match AnchorCellID -> RawData$cellID
    idx <- match(neigh$AnchorCellID, rawData$cellID)

    # warn if something doesn't match
    if (any(is.na(idx))) 
    {
      warning("Some AnchorCellIDs not found in RawData.")
    }

    neigh$UMAP1 <- rawData$UMAP1[idx]
    neigh$UMAP2 <- rawData$UMAP2[idx]

    locoResult$Neighbourhoods <- neigh
  }

  # -----------------------------
  # 5. add UMAP coords to Correlations
  # -----------------------------
  if ("Correlations" %in% names(locoResult)) 
  {
    corr <- locoResult$Correlations

    if (!("NeighbourhoodID" %in% colnames(corr))) 
    {
      stop("Correlations must contain 'NeighbourhoodID'")
    }

    if (!("NeighbourhoodID" %in% colnames(neigh))) 
    {
      stop("Neighbourhoods must contain 'NeighbourhoodID'")
    }

    # match Correlations$NeighbourhoodID -> Neighbourhoods$NeighbourhoodID
    idx_corr <- match(corr$NeighbourhoodID, neigh$NeighbourhoodID)

    # warn if missing
    if (any(is.na(idx_corr))) 
    {
      warning("Some NeighbourhoodIDs in Correlations not found in Neighbourhoods.")
    }

    # assign coords (this automatically repeats correctly)
    corr$UMAP1 <- neigh$UMAP1[idx_corr]
    corr$UMAP2 <- neigh$UMAP2[idx_corr]

    locoResult$Correlations <- corr
  }

  return(locoResult)
}

#' Plot correlations for a correlationPair in all neighbourhoods
#' 
#' @description
#' Create a plotting space first (like UMAP coordinates) to then plot neighbourhoods into this space by, e.g., running 
#' \code{newLocoResult <- add_umap_coords(locoResult)} followed by
#' \code{plot_local_correlation_map(newLocoResult, "GENE1_GENE2")}
#' 
#' @import ggplot2
#' @param locoResult the result of run_loco after adding an embedding space to the data.frame Correlations in the loco-result. This can be generated by, e.g., running \code{locoResult2 <- add_umap_coords(locoResult1)} and then feeding \code{locoResult2} into this function.
#' @param correlationPair the name of the correlation that should be plotted. You can see all possible names in the column \code{FeaturePair} in \code{locoResult1$LaplacianScores}.
#' @param dim1 the x-axis dimension for the plot: per default the UMAP1 coordinate from add_umap_coords
#' @param dim2 the y-axis dimension for the plot: per default the UMAP2 coordinate from add_umap_coords
#' @export
plot_local_correlation_map <- function(locoResult, correlationPair, dim1 = "UMAP1", dim2 = "UMAP2") 
{
  # -----------------------------
  # checks
  # -----------------------------
  if (!is.list(locoResult)) {
    stop("locoResult must be a list")
  }

  if (is.null(locoResult$Correlations)) {
    stop("locoResult has no Correlations element")
  }

  # subset the Correlation pair
  correlation_col <- as.character(locoResult$Correlations$CorrelationPair)
  idx <- which(trimws(correlation_col) == trimws(correlationPair))
  if (length(idx) == 0) {
    stop("CorrelationPair not found: ", correlationPair)
  }
  corr <- locoResult$Correlations[idx, , drop = FALSE]

  # required columns: Correlation and the space dimensions
  req_cols <- c("Correlation", dim1, dim2)

  missing <- setdiff(req_cols, colnames(corr))
  if (length(missing) > 0) {
    stop("Missing required columns in Correlations: ",
         paste(missing, collapse = ", "))
  }

  # -----------------------------
  # plot
  # -----------------------------
  p <- ggplot2::ggplot(
    corr,
    ggplot2::aes(
      x = .data[[dim1]],
      y = .data[[dim2]],
      color = Correlation
    )
  ) +
    ggplot2::geom_point(
      size = 4,
      alpha = 0.75
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      x = dim1,
      y = dim2,
      color = "Correlation"
    ) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey90"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold")
    )+
    ggplot2::scale_color_gradient2(
      low = "#0000FF",
      mid = "white",
      high = "firebrick",
      midpoint = 0,
      limits = c(-1, 1),
      oob = scales::squish
    )

  return(p)
}


#' Plot cell-level correlation within a spatial window
#' 
#' @description
#' Plot correlation between featureA and featureB as a scatter plot using all cells contained in neighbourhoods that are within given space boundaries
#' [\code{x_min} ... \code{x_max}] and [\code{y_min} ... \code{y_max}].
#' 
#' @param locoResult The UMAP-annotated result of run_loco
#' @param featureA first feature (x-coordinate) of the plotted correlation
#' @param featureB second feature (y-coordinate) of the plotted correlation
#' @param x_min The minimum x-coordinate to filter neighbourhoods that are used for plotting. Cells in neighbourhoods are included if the anchor cell of this neighbourhood has an x-coordinate >= \code{x_min}.
#' @param x_max The maximum x-coordinate to filter neighbourhoods that are used for plotting. Cells in neighbourhoods are included if the anchor cell of this neighbourhood has an x-coordinate <= \code{x_max}.
#' @param y_min The minimum y-coordinate to filter neighbourhoods that are used for plotting. Cells in neighbourhoods are included if the anchor cell of this neighbourhood has a y-coordinate >= \code{y_min}.
#' @param y_max The maximum y-coordinate to filter neighbourhoods that are used for plotting. Cells in neighbourhoods are included if the anchor cell of this neighbourhood has a y-coordinate <= \code{y_max}.
#' @param dim1 the x-axis dimension for the plot: per default the UMAP1 coordinate from add_umap_coords
#' @param dim2 the y-axis dimension for the plot: per default the UMAP2 coordinate from add_umap_coords
#' @export
plot_cell_level_correlation <- function(locoResult, featureA, featureB, x_min, x_max, y_min, y_max, dim1 = "UMAP1", dim2 = "UMAP2")
{

  # -----------------------------
  # checks
  # -----------------------------
  if (!is.list(locoResult)) {
    stop("locoResult must be a list")
  }

  if (is.null(locoResult$RawData) || is.null(locoResult$Neighbourhoods)) {
    stop("locoResult must contain RawData and Neighbourhoods")
  }

  raw <- locoResult$RawData
  neigh <- locoResult$Neighbourhoods

  if (!all(c(featureA, featureB) %in% colnames(raw))) {
    stop("featureA / featureB not found in RawData")
  }

  if (!all(c(dim1, dim2) %in% colnames(neigh))) {
    stop("dim1/dim2 not found in Neighbourhoods")
  }

  # -----------------------------
  # 1. filter neighbourhoods by spatial window
  # -----------------------------
  neigh_filt <- neigh[
    neigh[[dim1]] >= x_min &
    neigh[[dim1]] <= x_max &
    neigh[[dim2]] >= y_min &
    neigh[[dim2]] <= y_max,
  , drop = FALSE]

  if (nrow(neigh_filt) == 0) {
    stop("No neighbourhoods in selected region")
  }

  # -----------------------------
  # 2. collect all cell IDs (comma-separated expansion)
  # -----------------------------
  cell_list <- unlist(strsplit(as.character(neigh_filt$AllCellIDs), ","))

  # clean whitespace
  cell_list <- trimws(cell_list)

  # unique cells
  cell_list <- unique(cell_list)

  # -----------------------------
  # 3. subset raw data
  # -----------------------------
  raw_filt <- raw[raw$cellID %in% cell_list, , drop = FALSE]

  if (nrow(raw_filt) == 0) {
    stop("No matching cells found in RawData")
  }

  # -----------------------------
  # 4. plot
  # -----------------------------
  p <- ggplot2::ggplot(
    raw_filt,
    ggplot2::aes(x = .data[[featureA]], y = .data[[featureB]])
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste0(featureA, " vs ", featureB),
      x = featureA,
      y = featureB
    )

  return(p)
}