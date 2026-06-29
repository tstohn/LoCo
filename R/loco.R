#' Run LoCo analysis
#'
#' @description
#' Runs the complete LoCo (Local Correlations) pipeline on a feature matrix.
#' LoCo samples random anchor cells, builds fixed-size neighbourhoods around
#' them, computes pairwise feature correlations within each neighbourhood, and
#' scores the resulting correlation pairs with a Laplacian smoothness criterion.
#' Pairs with a low Laplacian score vary systematically across the point cloud
#' (e.g. a gene pair that is tightly co-expressed only in one cell sub-population)
#' rather than appearing randomly.  A permutation test gives each pair a p-value.
#'
#' @param inFile Path to the input data file.  The file must be a delimited text
#'   matrix where rows are cells and columns are features (genes, proteins, ...),
#'   or vice versa -- see \code{col} and \code{row}.
#' @param del Single character used as the column delimiter in \code{inFile}.
#'   Default: \code{"\\t"} (tab-separated).
#' @param col Logical.  If \code{TRUE}, feature names are in the \emph{columns}
#'   (i.e. each row is a cell and each column is a feature -- the standard layout).
#'   If \code{FALSE} (default), LoCo treats the file as transposed (rows = features,
#'   columns = cells).
#' @param row Logical.  If \code{TRUE}, cell (row) names are present as the first
#'   column of the file and will be read as cell IDs.  If \code{FALSE} (default),
#'   LoCo assigns automatic IDs of the form \code{C_<index>} starting from 0.
#' @param zscore Logical.  If \code{TRUE} (default), each feature is z-score
#'   normalised (zero mean, unit variance) before correlation is computed.
#'   Recommended for data with very different dynamic ranges across features.
#' @param thread Integer.  Number of parallel threads to use.  Default: \code{1}.
#'   Increase for large datasets to speed up neighbourhood and correlation
#'   calculations.
#' @param correlatedSetMode Integer controlling which graph algorithm is used to
#'   group co-correlated features into \emph{feature sets}.  Only relevant when
#'   \code{calcFeatureSets = TRUE}.  Accepted values:
#'   \describe{
#'     \item{\code{0} -- Clique / fully-connected}{Finds maximal cliques using the
#'       Bron-Kerbosch algorithm.  Every feature in a returned set is directly
#'       correlated with every other feature in that set.  Produces the most
#'       stringent, compact sets but is computationally expensive for large graphs.}
#'     \item{\code{1} -- Connected components (default)}{Groups features into
#'       weakly connected subgraphs.  Two features end up in the same set if
#'       they are linked through any chain of pairwise correlations, even if
#'       not all pairs are directly correlated.  Fast and suitable for most use
#'       cases.}
#'     \item{\code{>= 2} -- Minimum-edge connectivity}{Each feature in a returned
#'       set must share at least \code{correlatedSetMode} direct correlations
#'       with other members of that set.  For example, \code{correlatedSetMode = 3}
#'       requires every feature to have at least 3 correlated partners within
#'       the set.  Higher values yield denser, more tightly interconnected sets.}
#'   }
#' @param numberCorrelations Integer.  Maximum number of feature pairs to
#'   evaluate.  \code{0} (default) computes all \eqn{n(n-1)/2} pairwise
#'   correlations.  Set to a positive integer to subsample pairs, which is
#'   useful for very large feature spaces where exhaustive computation would
#'   be prohibitive.
#' @param cellStateGeneFile Path to a plain-text file (one feature name per
#'   line, no header) listing the features to be used \emph{exclusively} for
#'   neighbourhood construction.  Only these features drive the KNN graph that
#'   defines which cells are neighbours.  Useful when a curated gene set (e.g.
#'   cell-type marker genes) should define cell proximity independently of the
#'   features being correlated.  Leave as \code{""} (default) to use all
#'   features for neighbourhood construction.
#' @param correlationStateGeneFile Path to a plain-text file (one feature name
#'   per line, no header) listing the features for which pairwise correlations
#'   are computed.  All other features are ignored during the correlation step.
#'   Leave as \code{""} (default) to compute correlations for all features.
#' @param numberNeighbourhoods Integer.  Total number of neighbourhoods to
#'   create.  \code{0} (default) lets LoCo automatically choose
#'   \eqn{n_{\mathrm{cells}} / \mathrm{neighbourhoodSize}} neighbourhoods so that
#'   the point cloud is covered roughly once.  Increase for denser coverage or
#'   set explicitly to cap compute time.
#' @param neighbourhoodSize Integer.  Number of cells in each neighbourhood
#'   (anchor cell + its nearest neighbours).  Default: \code{100}.  Larger
#'   neighbourhoods capture broader, more global co-expression patterns;
#'   smaller neighbourhoods capture finer local variation.  The neighbourhood
#'   size must be meaningfully smaller than the total number of cells.
#' @param neighbourhoodKNN Integer.  Number of nearest neighbours used to build
#'   the KNN graph from which neighbourhoods are grown.  Default: \code{5}.
#'   This controls the graph connectivity used when expanding from the anchor
#'   cell to fill the neighbourhood to \code{neighbourhoodSize} cells.
#' @param correlationCutoff Numeric in \eqn{[0, 1]}.  Minimum absolute
#'   correlation required for a feature pair to be recorded within a
#'   neighbourhood.  Default: \code{0.5}.  Pairs whose correlation falls below
#'   this threshold in a given neighbourhood are discarded for that
#'   neighbourhood.  Lower values retain weaker correlations; higher values
#'   restrict to strongly correlated pairs only.
#' @param permutations Integer.  Number of permutations used to estimate the
#'   null distribution for the Laplacian score p-value.  Default: \code{100}.
#'   Increase (e.g. to 1000) for more accurate p-values.
#' @param minSetSize Integer.  Minimum number of features required to report a
#'   feature set when \code{calcFeatureSets = TRUE}.  Sets smaller than this
#'   threshold are discarded.  Default: \code{2}.
#' @param corrSetAbundance Numeric in \eqn{[0, 1]}.  Minimum fraction of
#'   neighbourhoods in which a feature pair must show a correlation above
#'   \code{correlationCutoff} to be considered for Laplacian scoring.
#'   Default: \code{0.01} (1 \%).  Increase to focus on correlations that
#'   appear consistently across the dataset; decrease to also capture very
#'   rare, localised correlations.
#' @param correlationType Character string, either \code{"spearman"} (default)
#'   or \code{"pearson"}.  Type of correlation computed within each
#'   neighbourhood.  Spearman correlation is rank-based and more robust to
#'   outliers; Pearson correlation is faster but sensitive to non-linear
#'   relationships and extreme values.
#' @param calcFeatureSets Logical.  If \code{TRUE}, LoCo additionally groups
#'   co-correlated features into \emph{feature sets} (modules) using the
#'   algorithm selected by \code{correlatedSetMode}.  A feature-level graph is
#'   built across all neighbourhoods: an edge is drawn between two features if
#'   their pair passes the \code{corrSetAbundance} and \code{correlationCutoff}
#'   filters globally (not restricted to a single spatial region).  The
#'   resulting feature sets are reported in the \code{FeatureSet} column of
#'   \code{LaplacianScores}.  Default: \code{FALSE}.
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
#' \strong{Pipeline steps executed by \code{run_loco}:}
#' \enumerate{
#'   \item The input matrix is parsed and, if requested, z-score normalised.
#'   \item \code{numberNeighbourhoods} anchor cells are sampled at random.
#'         A KNN graph (k = \code{neighbourhoodKNN}) is built from the cell
#'         coordinates (using all features or only those in
#'         \code{cellStateGeneFile}), and each anchor cell collects its
#'         \code{neighbourhoodSize} nearest neighbours to form a neighbourhood.
#'   \item Within each neighbourhood, pairwise correlations
#'         (\code{correlationType}) are computed for all feature pairs (or the
#'         subset in \code{correlationStateGeneFile}).  Pairs below
#'         \code{correlationCutoff} are discarded for that neighbourhood.
#'   \item Pairs that survive the \code{corrSetAbundance} filter (present in
#'         at least that fraction of all neighbourhoods) are scored with the
#'         Laplacian smoothness criterion and given a permutation-based p-value.
#'   \item Optionally (\code{calcFeatureSets = TRUE}), features are grouped into
#'         co-correlated modules via the graph algorithm chosen by
#'         \code{correlatedSetMode}.
#' }
#'
#' \strong{Interpreting the Laplacian score:} A \emph{low} score means the
#' correlation is spatially smooth -- it varies gradually and systematically
#' across the point cloud, suggesting a biologically meaningful local
#' co-regulation.  A high score means the correlation is noisy or random.
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
  correlationType = "spearman",
  calcFeatureSets = FALSE
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
  if (!is.logical(calcFeatureSets)) {
    stop("`calcFeatureSets` must be TRUE or FALSE")
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
    cellStateGeneFile,
    correlationStateGeneFile,
    as.integer(numberNeighbourhoods),
    as.integer(neighbourhoodSize),
    as.integer(neighbourhoodKNN),
    correlationCutoff,
    as.integer(permutations),
    as.integer(minSetSize),
    corrSetAbundance,
    correlationType,
    calcFeatureSets
  )

  return(res)
}


#' Add UMAP coordinates to a LoCo result
#'
#' @description
#' Computes a UMAP embedding of the raw data returned by \code{\link{run_loco}}
#' and attaches the 2-D coordinates (\code{UMAP1}, \code{UMAP2}) to the
#' \code{RawData}, \code{Neighbourhoods}, and \code{Correlations} data frames
#' in the result list.  The enriched list can then be passed directly to
#' \code{\link{plot_local_correlation_map}} or
#' \code{\link{plot_cell_level_correlation}}.
#'
#' @details
#' To visualise where in the point cloud a correlation is strong or weak, LoCo
#' needs a 2-D layout.  This function provides one by running UMAP on the
#' feature matrix stored in \code{locoResult$RawData}.  Each neighbourhood is
#' then assigned the UMAP coordinates of its anchor cell, so that
#' neighbourhood-level statistics (e.g. per-neighbourhood correlation strength)
#' can be plotted in the same 2-D space as the individual cells.
#'
#' You are not required to use this function.  If you have an existing embedding
#' (PCA, t-SNE, spatial coordinates, ...), you can add columns \code{UMAP1} and
#' \code{UMAP2} (or any name of your choice) to \code{locoResult$Correlations}
#' and \code{locoResult$Neighbourhoods} manually and pass custom dimension names
#' via the \code{dim1}/\code{dim2} arguments of the plot functions.
#'
#' @param locoResult Named list returned by \code{\link{run_loco}}.
#' @param n_pcs Integer.  Number of principal components to compute before
#'   running UMAP.  Default \code{0} skips PCA and runs UMAP directly on the
#'   feature matrix.  Set to a positive integer (e.g. \code{30}) to first
#'   reduce dimensionality with PCA -- recommended for datasets with many
#'   features, as it speeds up UMAP and often improves the embedding quality.
#'   When \code{n_pcs > 0}, \code{scale} is automatically set to \code{FALSE}
#'   because PCA already centres and scales the data.
#' @param scale Logical.  If \code{TRUE} (default), z-score normalise the
#'   feature matrix before computing UMAP.  Ignored when \code{n_pcs > 0}.
#' @param metric Character string.  Distance metric passed to
#'   \code{\link[uwot]{umap}}.  Default: \code{"euclidean"}.  Other options
#'   supported by \pkg{uwot} (e.g. \code{"cosine"}) are also accepted.
#' @param seed Integer.  Random seed for reproducible UMAP layouts.
#'   Default: \code{7}.
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

#' Plot per-neighbourhood correlation strength across an embedding
#'
#' @description
#' Visualises how the correlation of a chosen feature pair varies spatially
#' across the point cloud.  Each neighbourhood is drawn as a point at the
#' position of its anchor cell in the 2-D embedding (e.g. UMAP), coloured by
#' the Spearman/Pearson correlation of the pair inside that neighbourhood.
#' Blue indicates negative correlation, white indicates no correlation, and
#' red indicates positive correlation.
#'
#' Run \code{\link{add_umap_coords}} first to attach UMAP coordinates, then
#' call this function:
#' \preformatted{
#' result2 <- add_umap_coords(result)
#' plot_local_correlation_map(result2, result2$LaplacianScores$FeaturePair[1])
#' }
#'
#' @import ggplot2
#' @param locoResult Named list returned by \code{\link{run_loco}} after
#'   embedding coordinates have been added (e.g. via
#'   \code{\link{add_umap_coords}}).  The \code{Correlations} data frame must
#'   contain the columns named by \code{dim1} and \code{dim2}.
#' @param correlationPair Character string.  Name of the feature pair to plot,
#'   formatted as \code{"featureA_featureB"}.  Valid names are listed in the
#'   \code{FeaturePair} column of \code{locoResult$LaplacianScores}.
#' @param dim1 Character string.  Name of the column in
#'   \code{locoResult$Correlations} to use as the x-axis.
#'   Default: \code{"UMAP1"} (set by \code{\link{add_umap_coords}}).
#' @param dim2 Character string.  Name of the column in
#'   \code{locoResult$Correlations} to use as the y-axis.
#'   Default: \code{"UMAP2"} (set by \code{\link{add_umap_coords}}).
#' @return A \code{ggplot2} object.  Save with \code{ggplot2::ggsave()}.
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
#' Produces a cell-level scatter plot of \code{featureA} vs \code{featureB}
#' using only the cells that belong to neighbourhoods whose anchor falls inside
#' the specified 2-D bounding box (\code{x_min}–\code{x_max},
#' \code{y_min}–\code{y_max}).
#'
#' This is the companion to \code{\link{plot_local_correlation_map}}: first use
#' the map to visually identify a region of the embedding where a pair shows
#' strong (or weak) correlation, then use this function to inspect the
#' underlying cell-level relationship in that region.
#'
#' @param locoResult Named list returned by \code{\link{run_loco}} after
#'   embedding coordinates have been added (e.g. via
#'   \code{\link{add_umap_coords}}).
#' @param featureA Character string.  Name of the first feature (plotted on the
#'   x-axis of the scatter).  Must be a column in \code{locoResult$RawData}.
#' @param featureB Character string.  Name of the second feature (plotted on
#'   the y-axis of the scatter).  Must be a column in \code{locoResult$RawData}.
#' @param x_min Numeric.  Lower bound of the x-axis filter applied to
#'   neighbourhood anchor positions.  Neighbourhoods whose anchor has a
#'   \code{dim1} coordinate \strong{below} this value are excluded.
#' @param x_max Numeric.  Upper bound of the x-axis filter.  Neighbourhoods
#'   whose anchor has a \code{dim1} coordinate \strong{above} this value are
#'   excluded.
#' @param y_min Numeric.  Lower bound of the y-axis filter applied to
#'   neighbourhood anchor positions.  Neighbourhoods whose anchor has a
#'   \code{dim2} coordinate \strong{below} this value are excluded.
#' @param y_max Numeric.  Upper bound of the y-axis filter.  Neighbourhoods
#'   whose anchor has a \code{dim2} coordinate \strong{above} this value are
#'   excluded.
#' @param dim1 Character string.  Column in \code{locoResult$Neighbourhoods}
#'   used as the x-axis coordinate for spatial filtering.
#'   Default: \code{"UMAP1"}.
#' @param dim2 Character string.  Column in \code{locoResult$Neighbourhoods}
#'   used as the y-axis coordinate for spatial filtering.
#'   Default: \code{"UMAP2"}.
#' @return A \code{ggplot2} scatter plot of \code{featureA} vs \code{featureB}
#'   for all cells in the selected spatial window.
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