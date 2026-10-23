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
#' @param minNonNAFraction Numeric in \eqn{[0, 1]}.  Minimum fraction of
#'   neighbourhoods in which a feature pair must have a defined (non-NA)
#'   correlation.  Pairs with a non-NA fraction of at most this value are
#'   removed from \code{LaplacianScores} and \code{Correlations}.  A
#'   correlation is NA in a neighbourhood when it was not calculated there
#'   (e.g. not enough valid values).  For example \code{0.75} keeps only pairs
#'   with a correlation in more than 75\% of all neighbourhoods.  Default:
#'   \code{0} (no filtering).  This filter is applied right after the
#'   correlations and BEFORE the Laplacian score: removed pairs get no score
#'   and no p-value and do not count as tests for \code{p_adj}.  There is no
#'   other minimum: every remaining pair is scored.  Pairs whose Laplacian
#'   score is undefined (no two neighbouring neighbourhoods with a defined
#'   correlation, or no variation among them) are removed as well, so the
#'   result contains no pairs with an NA score.
#' @param numberCorrelations Integer.  Maximum number of SIGNIFICANT feature
#'   pairs to return.  \code{0} (default) disables this and returns all feature
#'   pairs.  If set to a positive integer, LoCo first calculates the Laplacian
#'   score, the permutation p-value (\code{p_value}) and the
#'   Benjamini-Hochberg adjusted p-value (\code{p_adj}) for ALL tested feature
#'   pairs (the filtering happens AFTER significance testing).  Pairs removed
#'   by \code{minNonNAFraction} were never tested and are not considered here.
#'   It then keeps only
#'   the significant pairs (\code{p_adj <= significanceCutoff}) and, of those,
#'   only the \code{numberCorrelations} pairs with the lowest Laplacian score
#'   (most local correlation).  All other pairs are removed from
#'   \code{LaplacianScores} and \code{Correlations}.  If fewer than
#'   \code{numberCorrelations} pairs are significant, all significant pairs are
#'   returned.
#'   Because the p-values and \code{p_adj} are calculated for all tested pairs
#'   before filtering, \code{p_adj} corrects for ALL tested pairs and is not
#'   affected by this filter.  The filter is applied BEFORE the feature sets
#'   are calculated, so the \code{FeatureSet} column (and the runtime of the
#'   feature set calculation) only uses the kept pairs.
#' @param significanceCutoff Numeric in \eqn{[0, 1]}.  Threshold on the
#'   adjusted p-value (\code{p_adj}) that defines a pair as significant.
#'   Only used when \code{numberCorrelations > 0}.  Default: \code{0.05}.
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
#'   the KNN graph of neighbourhoods.  Default: \code{5}.
#'   This controls the graph connectivity on which the Laplacian Score is calculated.
#'   Larger values result in correlations that are more smooth, smaller values
#'   can capture more local changes but are also more prone to capture noise.
#' @param correlationCutoff Numeric in \eqn{[0, 1]}.  Minimum absolute
#'   correlation required for a feature pair to be recorded within a
#'   neighbourhood.  Default: \code{0.5}.  Pairs whose correlation falls below
#'   this threshold in a given neighbourhood are discarded for that
#'   neighbourhood.  Lower values retain weaker correlations; higher values
#'   restrict to strongly correlated pairs only.
#' @param permutations Integer.  Number of permutations used to estimate the
#'   null distribution for the Laplacian score p-value.  Default: \code{500}.
#'   Increase (e.g. to 1000 or more) for more accurate p-values. The smallest possible
#'   p-value is \code{1 / (permutations + 1)}, so this also limits what survives
#'   multiple-testing correction.
#' @param minSetSize Integer.  Minimum number of features required to report a
#'   feature set when \code{calcFeatureSets = TRUE}.  Sets smaller than this
#'   threshold are discarded.  Default: \code{2}.
#' @param corrSetAbundance Numeric in \eqn{[0, 1]}.  Minimum fraction of
#'   neighbourhoods in which a feature pair must show a correlation above
#'   \code{correlationCutoff} to be considered for Laplacian scoring.
#'   Default: \code{0.01} (1 \%).  Increase to focus on correlations that
#'   appear consistently across the dataset; decrease to also capture very
#'   rare, localised correlations.
#' @param Nsampling Character string, either \code{"random"} (default) or
#'   \code{"kmeans"}.  Controls how neighbourhoods are chosen.
#'   \describe{
#'     \item{\code{"random"}}{Anchor cells are sampled uniformly at random;
#'       every neighbourhood is the anchor plus its \code{neighbourhoodSize}
#'       nearest cells (Manhattan distance).  Neighbourhoods can overlap.}
#'     \item{\code{"kmeans"}}{\emph{Non-overlapping} neighbourhoods by a
#'       size-constrained k-medians clustering (Manhattan distance) on the same
#'       cell-state features as the \code{"random"} mode (\code{cellStateGeneFile},
#'       or all features if no file is given).  The planned number of
#'       neighbourhoods is \code{numberNeighbourhoods} (default: number of cells /
#'       \code{neighbourhoodSize}).  Neighbourhood sizes may vary by
#'       \code{kmeansSlack} (default 5\%, i.e. 95 to 105 cells for
#'       \code{neighbourhoodSize = 100}); not every cell has to be in a
#'       neighbourhood.
#'       \enumerate{
#'         \item k-means++ seeding: the first centre is a random cell, every further
#'           centre a random cell drawn with probability proportional to its
#'           distance to the nearest centre so far (\code{seed}).
#'         \item k-medians iterations (at most 30, until at most 0.1\% of the cells
#'           change): all cells are assigned from scratch, cells with the strongest
#'           preference for one centre first, each to the nearest of its
#'           \code{kmeansR} (default 5) nearest centres with fewer than
#'           \code{neighbourhoodSize} cells; cells that do not fit then go to the
#'           nearest of these centres that is still below the maximum size.  A cell whose
#'           \code{kmeansR} nearest neighbourhoods are all full is \emph{not placed}
#'           in this round and considered again in the next one.  Neighbourhoods
#'           below the minimum size then take their nearest not-placed cells up to the
#'           minimum (so that no centre starves; cells placed far away this way are
#'           repaired in the swap phase).  From the 3rd iteration on, a neighbourhood
#'           that still stays below the minimum (per iteration only the smallest of
#'           neighbouring ones) is moved to the worst served cell (the cell furthest
#'           from its centre), so the other neighbourhoods of a small population take
#'           over its cells; in the last 5 iterations it is dissolved instead.  Every
#'           centre then moves to the median of its cells.
#'         \item Neighbourhoods with fewer than the minimum size are dissolved (they do
#'           not exist); their cells join one of their \code{kmeansR} nearest
#'           neighbourhoods with room if possible.
#'         \item Swap phase: every cell whose own neighbourhood is \emph{not among its
#'           \code{kmeansR} nearest neighbourhoods} is moved to the nearest of these
#'           with room, otherwise swapped with a cell of one of them if this reduces
#'           the total distance, otherwise left out of its neighbourhood (only if
#'           this keeps the minimum size).  Left-out cells join one of their
#'           \code{kmeansR} nearest neighbourhoods if one has room.
#'       }
#'       The anchor is the cell closest to the neighbourhood median.  A short summary
#'       is printed (number of neighbourhoods, \% of cells assigned to a
#'       neighbourhood, \% of these in their closest and in their 2nd to
#'       \code{kmeansR}-th closest neighbourhood); \code{NeighbourhoodStatus} in the
#'       result tells for every cell whether and why it is in no neighbourhood.}
#'   }
#' @param kmeansR Integer \eqn{\ge 1}.  Only for \code{Nsampling = "kmeans"}.  A
#'   cell is only placed into one of its \code{kmeansR} nearest neighbourhoods;
#'   cells outside them are repaired in the swap phase.  Default: \code{5}.
#' @param kmeansSlack Numeric in \eqn{[0, 1)}.  Only for
#'   \code{Nsampling = "kmeans"}.  Neighbourhoods may have between
#'   \eqn{\lceil s(1 - \code{kmeansSlack})\rceil} and
#'   \eqn{\lfloor s(1 + \code{kmeansSlack})\rfloor} cells
#'   (\eqn{s} = \code{neighbourhoodSize}).  Default: \code{0.05}.
#' @param seed Non-negative integer.  Seed of the random numbers of the
#'   neighbourhood construction (k-means++ seeding for \code{"kmeans"}, anchor
#'   cells for \code{"random"}).  \code{0}
#'   (default): a random seed.  The permutations of the significance test are
#'   not affected.
#' @param correlationType Character string, either \code{"spearman"} (default)
#'   or \code{"pearson"}.  Type of correlation computed within each
#'   neighbourhood.  Spearman correlation is rank-based and more robust to
#'   outliers; Pearson correlation is faster but sensitive to non-linear
#'   relationships and extreme values.
#' @param minVaryingFraction Numeric in \eqn{[0, 1]}.  Minimum fraction of the
#'   cells in a neighbourhood that must \emph{vary} for a feature, i.e. whose
#'   value differs from the most frequent value of this feature in the
#'   neighbourhood.  For sparse data such as scRNA-seq this is the fraction of
#'   cells expressing the gene (the most frequent value is 0).  In every
#'   neighbourhood where a feature is below this fraction, all correlations of
#'   this feature are set to NA (not calculated) in that neighbourhood.
#'   Default: \code{0.05}, so in a neighbourhood of 100 cells a gene must differ
#'   from its most frequent value in at least 5 cells.  \code{0} disables the
#'   filter (correlations are then only NA if a feature has no variation at all).
#'
#'   \strong{Why:} without this filter a feature that is, e.g., 0 in 99 cells
#'   and 1 in a single cell still has a defined correlation.  If the partner
#'   feature happens to be high in that one cell, the correlation is close to 1
#'   (Pearson and Spearman, since all zeroes share the same rank), although it
#'   is driven by a single cell.  Such correlations are spurious and add noise
#'   to the Laplacian score.
#'
#'   \strong{Effect on the other steps:} the filter is applied per feature and
#'   per neighbourhood WHILE the correlations are calculated, i.e. BEFORE the
#'   \code{correlationCutoff}/ \code{corrSetAbundance} filters, the Laplacian
#'   score and the permutation test.  The resulting NAs are treated like any
#'   other missing correlation (see \emph{Laplacian score and missing
#'   correlations} in Details): the neighbourhood is left out for that pair and
#'   NAs keep their position in the permutations.  Pairs of rarely expressed
#'   features therefore get more NAs and can be removed with
#'   \code{minNonNAFraction}.  Because the threshold is a fraction it adapts to
#'   the neighbourhood size (e.g. \code{Nsampling = "kmeans"}).  For continuous
#'   data without ties (e.g. protein intensities) almost every cell differs from
#'   the most frequent value, so the filter has no effect there.
#' @param considerOnlyNonZeroCells Logical.  If \code{TRUE}, the correlation of
#'   a pair in a neighbourhood is calculated only from the cells in which BOTH
#'   features are non-zero; cells with a zero in one of the two features are
#'   left out for this pair (for Spearman the ranks are calculated on these
#'   cells only).  The correlation is NA if fewer than 3 such cells remain.
#'   \strong{Only EXACT zeros are excluded}: a feature is zero in a cell only if
#'   its value in the input is exactly 0, so only real zero values (e.g. zero
#'   counts) are left out; small non-zero values are kept.  With
#'   \code{zscore = TRUE} there are no exact zeros left, so the filter has no
#'   effect -- use it with \code{zscore = FALSE} on counts or log-normalised data.  The correlation then measures
#'   how the LEVELS of two expressed features relate, and ignores whether they
#'   are switched on in the same cells (co-detection, e.g. mutually exclusive
#'   features).  Cells with only one feature expressed are not used, so the
#'   number of cells behind a correlation differs between pairs and
#'   neighbourhoods; combine it with \code{minFractionNonZeroCells} to require
#'   enough cells.  \code{minVaryingFraction} is still applied to all cells of
#'   the neighbourhood.  Default: \code{FALSE}.
#' @param minFractionNonZeroCells Numeric in \eqn{[0, 1]}.  Minimum fraction of
#'   the cells of a neighbourhood in which BOTH features of a pair are non-zero
#'   (\strong{only EXACT zeros} count as zero, i.e. values of exactly 0 in the input;
#'   no effect with \code{zscore = TRUE}, see \code{considerOnlyNonZeroCells}).  In neighbourhoods
#'   below it the correlation of this pair is set to NA.  This removes
#'   correlations that rest on very few co-expressing cells (e.g. 90\% of the
#'   cells zero in both features).  Applied per pair and neighbourhood WHILE the
#'   correlations are calculated, like \code{minVaryingFraction}; the NAs are
#'   treated like any other missing correlation.  Works with and without
#'   \code{considerOnlyNonZeroCells}.  Default: \code{0.1} (at least 10\% of the cells
#'   must have both features non-zero); \code{0} disables the filter.
#' @param calcFeatureSets Logical.  If \code{TRUE}, LoCo additionally groups
#'   co-correlated features into \emph{feature sets} (modules) using the
#'   algorithm selected by \code{correlatedSetMode}.  A feature-level graph is
#'   built across all neighbourhoods: an edge is drawn between two features if
#'   their pair passes the \code{corrSetAbundance} and \code{correlationCutoff}
#'   filters globally (not restricted to a single spatial region).  The
#'   resulting feature sets are reported in the \code{FeatureSet} column of
#'   \code{LaplacianScores}.  Default: \code{FALSE}.
#' @return A named list containing LoCo results with five elements:
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
#'       \item \code{LaplacianScore}: degree-weighted Laplacian score (0 = smooth, about 1 = no
#'         structure, maximum 2), calculated on the neighbourhoods with a defined correlation
#'       \item \code{p_value}: permutation-based p-value, calculated as
#'         \code{(k + 1) / (permutations + 1)} where \code{k} is the number of permutations
#'         with a Laplacian score at or below the observed one (so it is never 0). Permutations
#'         shuffle only the defined correlations, NAs keep their position.
#'       \item \code{p_adj}: Benjamini-Hochberg adjusted \code{p_value} (FDR) across all
#'         feature pairs with a defined p-value
#'       \item \code{FeatureSet}: comma-separated list of features forming sets of co-correlated features
#'     }
#'   }
#'
#'   \item{Correlations}{
#'     A numeric matrix of the correlation of every pair in every neighbourhood: one row per neighbourhood
#'     (row names = \code{NeighbourhoodID}, like \code{N_<number>}, in the same order as the rows of
#'     \code{Neighbourhoods}) and one column per feature pair (column names = pair names as in
#'     \code{LaplacianScores$FeaturePair}). \code{NA} where the correlation is not defined.
#'     \code{Correlations[, "A_B"]} are the correlations of pair \code{A_B} across all neighbourhoods.
#'     Use \code{\link{get_correlations}} for the long format (columns \code{CorrelationPair},
#'     \code{NeighbourhoodID}, \code{Correlation}).
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
#'
#'   \item{NeighbourhoodStatus}{
#'     A data.frame with one row per cell: \code{cellID} and \code{Status},
#'     \code{"placed"} if the cell is in a neighbourhood, otherwise why not:
#'     \code{"no_room"} (all its \code{kmeansR} nearest neighbourhoods were full),
#'     \code{"dissolved"} (its neighbourhood did not reach the minimum size),
#'     \code{"dropped"} (left out in the swap phase) for \code{"kmeans"}, and
#'     \code{"not_in_neighbourhood"} for \code{"random"}.
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
#'   \item Neighbourhoods are built in the cell-state space (all features or
#'         only those in \code{cellStateGeneFile}) as chosen by \code{Nsampling}:
#'         for \code{"random"}, \code{numberNeighbourhoods} anchor cells are
#'         sampled at random and each collects its \code{neighbourhoodSize}
#'         nearest cells; \code{"kmeans"} builds
#'         non-overlapping neighbourhoods.  A KNN graph (k = \code{neighbourhoodKNN})
#'         of the anchor cells connects the neighbourhoods.
#'   \item Within each neighbourhood, pairwise correlations
#'         (\code{correlationType}) are computed for all feature pairs (or the
#'         subset in \code{correlationStateGeneFile}).  Features that vary in
#'         fewer than \code{minVaryingFraction} of the cells of a neighbourhood
#'         get NA correlations in this neighbourhood.  Pairs below
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
#'
#' \strong{Laplacian score and missing correlations:} The score is the
#' Laplacian score of He, Cai and Niyogi (2005) of the correlations of a feature
#' pair across the neighbourhood graph:
#' \eqn{L = \sum_{i<j} w_{ij} (c_i - c_j)^2 / \sum_i d_i (c_i - \mu_w)^2}, where
#' \eqn{c_i} is the correlation in neighbourhood \eqn{i}, \eqn{w_{ij}} the edge
#' weight, \eqn{d_i = \sum_j w_{ij}} the degree and
#' \eqn{\mu_w = \sum_i d_i c_i / \sum_i d_i} the degree-weighted mean.  A score of
#' 0 is a perfectly smooth signal, about 1 means no relation to the graph, and
#' the maximum is 2.  A correlation is undefined (NA) in a neighbourhood if a
#' feature has no variation in it, or if it varies in fewer than
#' \code{minVaryingFraction} of the cells of this neighbourhood.  Such neighbourhoods are removed from the
#' graph for this pair (the remaining edges and weights are kept, nothing is
#' imputed), and the score, degrees and mean are calculated on the remaining
#' neighbourhoods only.  Pairs with a defined correlation in at most
#' \code{minNonNAFraction} of all neighbourhoods are removed before scoring, and
#' pairs whose score is undefined (no edge between two neighbourhoods with a
#' defined correlation, or no variation) are removed from the result, so no
#' returned pair has an NA score.  For the permutation test the
#' NAs stay at their position and only the defined correlations are shuffled
#' among the neighbourhoods with a defined correlation.
#' @export
run_loco <- function(
  inFile,
  del = "\t",
  col = FALSE,
  row = FALSE,
  zscore = TRUE,
  thread = 1,
  correlatedSetMode = 1,
  minNonNAFraction = 0,
  numberCorrelations = 0,
  significanceCutoff = 0.05,
  cellStateGeneFile = "",
  correlationStateGeneFile = "",
  numberNeighbourhoods = 0,
  neighbourhoodSize = 100,
  neighbourhoodKNN = 5,
  correlationCutoff = 0.5,
  permutations = 500,
  minSetSize = 2,
  corrSetAbundance = 0.01,
  correlationType = "spearman",
  minVaryingFraction = 0.05,
  considerOnlyNonZeroCells = FALSE,
  minFractionNonZeroCells = 0.1,
  calcFeatureSets = FALSE,
  Nsampling = "random",
  kmeansR = 5,
  kmeansSlack = 0.05,
  seed = 0
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
  if (!is.numeric(minNonNAFraction) || minNonNAFraction < 0 || minNonNAFraction > 1) {
    stop("`minNonNAFraction` must be between 0 and 1")
  }
  if (!is.numeric(numberCorrelations) || numberCorrelations < 0) {
    stop("`numberCorrelations` must be >= 0")
  }
  if (!is.numeric(significanceCutoff) || significanceCutoff < 0 || significanceCutoff > 1) {
    stop("`significanceCutoff` must be between 0 and 1")
  }
  if (!is.numeric(minVaryingFraction) || length(minVaryingFraction) != 1 ||
      minVaryingFraction < 0 || minVaryingFraction > 1) {
    stop("`minVaryingFraction` must be between 0 and 1")
  }
  if (!is.logical(considerOnlyNonZeroCells) || length(considerOnlyNonZeroCells) != 1 || is.na(considerOnlyNonZeroCells)) {
    stop("`considerOnlyNonZeroCells` must be TRUE or FALSE")
  }
  if (!is.numeric(minFractionNonZeroCells) || length(minFractionNonZeroCells) != 1 ||
      minFractionNonZeroCells < 0 || minFractionNonZeroCells > 1) {
    stop("`minFractionNonZeroCells` must be between 0 and 1")
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

  if (!is.character(Nsampling) || length(Nsampling) != 1 || !(Nsampling %in% c("random", "kmeans"))) {
    stop("`Nsampling` must be 'random' or 'kmeans'.")
  }
  if (!is.numeric(kmeansR) || length(kmeansR) != 1 || kmeansR < 1) {
    stop("`kmeansR` must be >= 1")
  }
  if (!is.numeric(kmeansSlack) || length(kmeansSlack) != 1 || kmeansSlack < 0 || kmeansSlack >= 1) {
    stop("`kmeansSlack` must be >= 0 and < 1")
  }
  if (!is.numeric(seed) || length(seed) != 1 || seed < 0 || seed != round(seed)) {
    stop("`seed` must be a non-negative integer (0 = random)")
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
    calcFeatureSets,
    Nsampling,
    minNonNAFraction,
    as.integer(numberCorrelations),
    significanceCutoff,
    minVaryingFraction,
    as.integer(kmeansR),
    kmeansSlack,
    as.integer(seed),
    considerOnlyNonZeroCells,
    minFractionNonZeroCells
  )
  # minNonNAFraction is applied in C++ BEFORE the Laplacian score (removed pairs are not scored/ tested);
  # p_adj (Benjamini-Hochberg over all tested pairs) and the numberCorrelations filter are applied after the
  # p-values and BEFORE the feature sets are calculated

  return(res)
}


#' Add UMAP coordinates to a LoCo result
#'
#' @description
#' Computes a UMAP embedding of the raw data returned by \code{\link{run_loco}}
#' and attaches the 2-D coordinates (\code{UMAP1}, \code{UMAP2}) to the
#' \code{RawData} and \code{Neighbourhoods} data frames in the result list
#' (\code{Correlations} is a matrix whose rows are the rows of \code{Neighbourhoods}).  The enriched list can then be passed directly to
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
#' \code{UMAP2} (or any name of your choice) to \code{locoResult$Neighbourhoods}
#' manually (e.g. the coordinates of the anchor cell of every neighbourhood) and pass
#' custom dimension names via the \code{dim1}/\code{dim2} arguments of the plot functions.
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

  return(locoResult)
}

#' Correlations of a LoCo result in long format
#'
#' @description
#' \code{run_loco} returns \code{Correlations} as a matrix (rows = neighbourhoods, columns = feature
#' pairs) because the long format needs one row per pair and neighbourhood (easily hundreds of millions of
#' rows). This function returns the long format for all or only some pairs.
#'
#' @param locoResult Named list returned by \code{\link{run_loco}}. Results of older versions of LoCo
#'   (\code{Correlations} already in long format) are accepted as well.
#' @param pairs Character vector of pair names (\code{"featureA_featureB"}, as in
#'   \code{LaplacianScores$FeaturePair}) or \code{NULL} (default) for all pairs. For large results
#'   select only the pairs you need: all pairs need about 3x the memory of the matrix.
#' @param neighbourhoodColumns Character vector of columns of \code{locoResult$Neighbourhoods} (e.g.
#'   \code{c("UMAP1", "UMAP2")} after \code{\link{add_umap_coords}}) that are added to every row by
#'   \code{NeighbourhoodID}, or \code{NULL} (default) for none.
#' @return A data.frame with the columns \code{CorrelationPair}, \code{NeighbourhoodID} and
#'   \code{Correlation} (plus \code{neighbourhoodColumns}): one row per pair and neighbourhood, ordered by
#'   neighbourhood (order of \code{Neighbourhoods}) and within a neighbourhood by pair (order of the matrix
#'   columns).
#' @examples
#' \dontrun{
#' corr <- get_correlations(result, pairs = result$LaplacianScores$FeaturePair[1:10])
#' result2 <- add_umap_coords(result)
#' corr_umap <- get_correlations(result2, pairs = "A_B", neighbourhoodColumns = c("UMAP1", "UMAP2"))
#' }
#' @export
get_correlations <- function(locoResult, pairs = NULL, neighbourhoodColumns = NULL)
{
  out <- correlations_long_format(locoResult, pairs)
  if (length(neighbourhoodColumns) > 0) {
    neigh <- locoResult$Neighbourhoods
    missing <- setdiff(neighbourhoodColumns, colnames(neigh))
    if (length(missing) > 0) stop("Missing columns in Neighbourhoods: ", paste(missing, collapse = ", "))
    idx <- match(out$NeighbourhoodID, neigh$NeighbourhoodID)
    if (anyNA(idx)) warning("Some NeighbourhoodIDs in Correlations not found in Neighbourhoods.")
    for (col in neighbourhoodColumns) out[[col]] <- neigh[[col]][idx]
  }
  out
}

# long format (CorrelationPair, NeighbourhoodID, Correlation) of the Correlations matrix for all or some pairs
correlations_long_format <- function(locoResult, pairs)
{
  if (!is.list(locoResult) || is.null(locoResult$Correlations)) {
    stop("locoResult must be the result of run_loco (with a Correlations element)")
  }
  corr <- locoResult$Correlations

  # result of an older LoCo version: already long format
  if (is.data.frame(corr)) {
    if (!is.null(pairs)) {
      missing <- setdiff(pairs, corr$CorrelationPair)
      if (length(missing) > 0) stop("CorrelationPair not found: ", paste(missing, collapse = ", "))
      corr <- corr[corr$CorrelationPair %in% pairs, , drop = FALSE]
    }
    return(corr)
  }

  if (!is.matrix(corr) || is.null(rownames(corr)) || is.null(colnames(corr))) {
    stop("Correlations must be a matrix with NeighbourhoodIDs as row names and pairs as column names")
  }
  if (!is.null(pairs)) {
    col_idx <- match(pairs, colnames(corr))
    if (anyNA(col_idx)) stop("CorrelationPair not found: ", paste(pairs[is.na(col_idx)], collapse = ", "))
    corr <- corr[, col_idx, drop = FALSE]
  }

  # a matrix without columns (no pairs) has no column names: keep the CorrelationPair column anyway
  pair_names <- if (is.null(colnames(corr))) character(0) else colnames(corr)
  data.frame(
    CorrelationPair = rep(pair_names, times = nrow(corr)),
    NeighbourhoodID = rep(rownames(corr), each = ncol(corr)),
    Correlation = as.vector(t(corr)),
    stringsAsFactors = FALSE
  )
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
#'   \code{\link{add_umap_coords}}).  The \code{Neighbourhoods} data frame must
#'   contain the columns named by \code{dim1} and \code{dim2}.
#' @param correlationPair Character string.  Name of the feature pair to plot,
#'   formatted as \code{"featureA_featureB"}.  Valid names are listed in the
#'   \code{FeaturePair} column of \code{locoResult$LaplacianScores}.
#' @param dim1 Character string.  Name of the column in
#'   \code{locoResult$Neighbourhoods} to use as the x-axis.
#'   Default: \code{"UMAP1"} (set by \code{\link{add_umap_coords}}).
#' @param dim2 Character string.  Name of the column in
#'   \code{locoResult$Neighbourhoods} to use as the y-axis.
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

  if (is.null(locoResult$Neighbourhoods)) {
    stop("locoResult has no Neighbourhoods element")
  }
  neigh <- locoResult$Neighbourhoods
  missing <- setdiff(c(dim1, dim2), colnames(neigh))
  if (length(missing) > 0) {
    stop("Missing required columns in Neighbourhoods: ", paste(missing, collapse = ", "),
         ". Run add_umap_coords() or add the coordinates of the anchor cells to Neighbourhoods.")
  }

  # correlations of this pair (one row per neighbourhood) + the coordinates of their neighbourhoods
  corr <- get_correlations(locoResult, pairs = trimws(correlationPair), neighbourhoodColumns = unique(c(dim1, dim2)))

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