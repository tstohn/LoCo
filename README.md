# LoCo
# Local Correlation Analysis

<img src=https://github.com/tstohn/LoCo/tree/master/docs/LoCo5.png width="200" />

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
