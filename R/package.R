#' @useDynLib loco, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

#' loco: Local Correlation Detection in High-Dimensional Point-Cloud Data
#'
#' @description
#' \if{html}{\figure{loco_logo.png}{options: width=120 alt='LoCo logo'}}
#' \if{latex}{\figure{loco_logo.png}{options: width=3cm}}
#'
#' **LoCo** (Local Correlations) detects feature pairs whose correlation
#' changes systematically across a high-dimensional point cloud -- for example,
#' gene pairs that are co-expressed only in a specific sub-population of cells
#' in single-cell RNA-seq data.
#'
#' @section Algorithm overview:
#' LoCo works in three stages:
#'
#' \enumerate{
#'   \item \strong{Neighbourhood construction} -- Random anchor cells are
#'         sampled and a KNN graph is built to collect a fixed number of
#'         nearest neighbours around each anchor.  Each anchor plus its
#'         neighbours forms one \emph{neighbourhood}.
#'
#'   \item \strong{Local correlation} -- Within every neighbourhood, pairwise
#'         Spearman (or Pearson) correlations are computed for all feature
#'         pairs (or a user-supplied subset).  Only correlations above
#'         \code{correlationCutoff} that appear in at least
#'         \code{corrSetAbundance} fraction of neighbourhoods are retained.
#'
#'   \item \strong{Laplacian scoring} -- The retained correlations are scored
#'         with a Laplacian smoothness criterion.  A low Laplacian score
#'         indicates that the correlation is spatially coherent (i.e., it
#'         changes smoothly and systematically across the point cloud rather
#'         than appearing randomly).  A permutation test yields p-values for
#'         each scored pair.
#' }
#'
#' @section Typical workflow:
#' \preformatted{
#' library(loco)
#'
#' # 1. Run LoCo on your data matrix (cells x features TSV)
#' result <- run_loco("my_data.tsv", correlationCutoff = 0.4,
#'                    neighbourhoodSize = 50)
#'
#' # 2. Inspect the top-scoring (most spatially coherent) pairs
#' head(result$LaplacianScores)
#'
#' # 3. Embed into UMAP space for visualisation
#' result <- add_umap_coords(result)
#'
#' # 4. Plot a specific correlation pair across the UMAP
#' p <- plot_local_correlation_map(result, result$LaplacianScores$FeaturePair[1])
#'
#' # 5. Zoom into a spatial region to see cell-level scatter
#' p2 <- plot_cell_level_correlation(result,
#'         featureA = "GeneA", featureB = "GeneB",
#'         x_min = -3, x_max = 0, y_min = -2, y_max = 2)
#' }
#'
#' @section Functions:
#' \describe{
#'   \item{\code{\link{run_loco}}}{Main entry point -- runs the full LoCo
#'         pipeline and returns neighbourhoods, per-neighbourhood correlations,
#'         and Laplacian-scored feature pairs.}
#'   \item{\code{\link{add_umap_coords}}}{Computes a UMAP embedding of the
#'         raw data and attaches UMAP coordinates to the result list so that
#'         correlations can be plotted in 2-D space.}
#'   \item{\code{\link{plot_local_correlation_map}}}{Plots a chosen feature
#'         pair's per-neighbourhood correlation strength as a colour-coded
#'         scatter across any 2-D embedding.}
#'   \item{\code{\link{plot_cell_level_correlation}}}{Plots a cell-level
#'         scatter of two features, restricted to cells that fall inside a
#'         user-defined spatial window.}
#' }
#'
#' @docType package
#' @name loco-package
#' @aliases loco
"_PACKAGE"
