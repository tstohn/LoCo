#pragma once

#include <string>
#include <Rcpp.h>
using namespace Rcpp;

//redefine commands for cran
#include "loco_io.h"

// undefine PI as both R and nanoflann define it
#undef PI 

// include your existing headers
#include "Neighborhood.hpp"
#include "GraphData.hpp"
#include "SCParser.hpp"

//return R-object of LoCo results
// create a List of several dataframes
// 1.) raw data table
// 2.) laplacian scores
// 3.) correlations: numeric matrix, rows = neighbourhoods (NeighbourhoodID), cols = correlation pairs
// 4,) neighbourhoods: anchor cell + contained cells
Rcpp::List build_loco_object(const SingleCellData& rawData,
                             Neighborhood& neighborhood,
                            const bool calcFeatureSets) 
{

    //fill intermittend data that are used to write the R-result object
    std::vector<std::string> nIDs; //all neighborhoods IDs
    std::vector<std::string> nID_anchorCellID; //vector (for each nID) off all anchor cellIDs (same order as nIDs)
    std::vector<std::vector<std::string>> nID_allCellIDs;// vector (for each nID) of vector of all cell IDs 

    std::vector<std::string> correlation_pairs_laplacian; //all names of the correlation pairs for Alplacian
    std::vector<std::string> correlation_pairs_origional; //all names of the correlation pairs for amtrix


    std::vector<double> corrL;
    std::vector<double> pCorrL;
    std::vector<double> pAdjL;
    std::vector<std::string> featureSetStrings; // a string for each feature-pair listing all sets with comma and semicolon like A,B,C;A,B,D
    neighborhood.fill_result_data(
        nIDs, nID_anchorCellID,nID_allCellIDs,
        correlation_pairs_laplacian, correlation_pairs_origional,
        corrL,pCorrL,pAdjL,featureSetStrings, calcFeatureSets);

    // =========================
    // 1. safe raw data table
    // =========================

    int nrow = rawData.pointCloud.size();
    int ncol = rawData.pointCloud[0].size();

    // 1. Create a List with space for the cellID column + all gene columns
    Rcpp::List df_columns(ncol + 1);
    Rcpp::CharacterVector df_names(ncol + 1);

    // 2. Set up the first column as cellID
    df_columns[0] = Rcpp::wrap(rawData.cellIDs);
    df_names[0] = "cellID";

    // 3. Extract columns efficiently from your pointCloud matrix
    for (int j = 0; j < ncol; j++) 
    {
        Rcpp::NumericVector col_data(nrow);
        for (int i = 0; i < nrow; i++) 
        {
            col_data[i] = rawData.pointCloud[i][j];
        }
        df_columns[j + 1] = col_data;
        df_names[j + 1] = rawData.geneNames[j];
    }

    // 4. Finalize the list attributes to make R treat it exactly as a DataFrame
    df_columns.attr("names") = df_names;
    df_columns.attr("class") = "data.frame";
    df_columns.attr("row.names") = Rcpp::wrap(rawData.cellIDs); // Uses cellIDs as row names too

    Rcpp::DataFrame raw_df = df_columns;
    
    // =========================
    // 2.) Laplacian scores/ p-values/ cliques for all interesting pairs
    // =========================

    //check that the vector of feature_pairs and laplacian-scores is of same length
    int laplacian_size = correlation_pairs_laplacian.size();
    // Safety check (important for robustness)
    if (corrL.size() != laplacian_size || pCorrL.size() != laplacian_size || pAdjL.size() != laplacian_size) 
    {
        Rcpp::stop("Input vectors must have the same length");
    }
    if(calcFeatureSets)
    {
        if(featureSetStrings.size() != laplacian_size)
        {
            Rcpp::stop("Input vectors must have the same length: list of feature sets and length of laplacian scores does not match!");
        }
    }

    Rcpp::DataFrame laplacian_scores;
    // build DataFrame 
    if(calcFeatureSets)
    {
        laplacian_scores = Rcpp::DataFrame::create(
            Rcpp::Named("FeaturePair") = correlation_pairs_laplacian,
            Rcpp::Named("LaplacianScore") = corrL,
            Rcpp::Named("p_value") = pCorrL,
            Rcpp::Named("FeatureSet") = featureSetStrings,
            Rcpp::Named("p_adj") = pAdjL,
            Rcpp::_["stringsAsFactors"] = false
        );
    }
    else
    {
        laplacian_scores = Rcpp::DataFrame::create(
            Rcpp::Named("FeaturePair") = correlation_pairs_laplacian,
            Rcpp::Named("LaplacianScore") = corrL,
            Rcpp::Named("p_value") = pCorrL,
            Rcpp::Named("p_adj") = pAdjL,
            Rcpp::_["stringsAsFactors"] = false
        );
    }

    // =========================
    // 3. correlations: MATRIX rows = neighbourhoods (nIDs), cols = correlation pairs
    // =========================
    // C++ stores the correlations row-major with row = pair, col = neighbourhood: the values of one pair are
    // contiguous. An R matrix is column-major, so a (neighbourhoods x pairs) R matrix has exactly the same memory
    // layout -> a single copy, no transpose and no strings per value (get_correlations() in R creates the long format)
    const FlatMatrix& pairCorrs = neighborhood.get_pair_correlations();
    int nIDsize = nIDs.size();
    int nPairs  = correlation_pairs_origional.size();
    if (pairCorrs.rows != static_cast<size_t>(nPairs) || pairCorrs.cols != static_cast<size_t>(nIDsize)) {
        Rcpp::stop("Correlation matrix must have one row per correlation pair and one column per neighbourhood");
    }

    // R_xlen_t: neighbourhoods x pairs can exceed the int range (2^31 - 1), R supports long vectors
    R_xlen_t total = static_cast<R_xlen_t>(nIDsize) * static_cast<R_xlen_t>(nPairs);
    Rcpp::NumericVector corr_mat(total);
    std::copy(pairCorrs.data.begin(), pairCorrs.data.begin() + total, corr_mat.begin());
    corr_mat.attr("dim") = Rcpp::IntegerVector::create(nIDsize, nPairs);
    corr_mat.attr("dimnames") = Rcpp::List::create(Rcpp::wrap(nIDs), Rcpp::wrap(correlation_pairs_origional));

    // =========================
    // 4. safe neighbourhoodIDs - anchorCell - allCells
    //.  allCells is missing (empty col)
    // =========================

    // Size checks
    if (nID_anchorCellID.size() != nIDsize) {
        Rcpp::stop("We need one anchor cell per neighbourhood ID: anchor cells = " + std::to_string(nID_anchorCellID.size()) + ", neighbourhood IDs = " + std::to_string(nIDsize));
    }
    if ( nID_allCellIDs.size() != nIDsize) {
        Rcpp::stop("We need one list of cell per neighbourhood ID: cell lists= " + std::to_string(nID_allCellIDs.size()) + ", neighbourhood IDs = " + std::to_string(nIDsize));
    }

    //  Flatten AllCellIDs 
    std::vector<std::string> allCellsCollapsed(nIDsize);

    for (int i = 0; i < nIDsize; i++) {
        const auto& cells = nID_allCellIDs.at(i);

        std::string combined;
        for (size_t j = 0; j < cells.size(); j++) 
        {
            combined += cells.at(j);
            if (j < cells.size() - 1) 
            {
                combined += ",";
            }
        }

        allCellsCollapsed.at(i) = combined;
    }

    // build DataFrame 
    Rcpp::DataFrame n_df = Rcpp::DataFrame::create(
        Rcpp::Named("NeighbourhoodID") = nIDs,
        Rcpp::Named("AnchorCellID")   = nID_anchorCellID,
        Rcpp::Named("AllCellIDs")     = allCellsCollapsed,
        Rcpp::_["stringsAsFactors"] = false
    );

    // =========================
    // 5. per cell: in a neighbourhood ("placed") or why not
    // =========================
    const std::vector<std::string>& cellStatus = neighborhood.get_cell_neighbourhood_status();
    if (cellStatus.size() != rawData.cellIDs.size()) {
        Rcpp::stop("We need one neighbourhood status per cell");
    }
    Rcpp::DataFrame status_df = Rcpp::DataFrame::create(
        Rcpp::Named("cellID") = rawData.cellIDs,
        Rcpp::Named("Status") = cellStatus,
        Rcpp::_["stringsAsFactors"] = false
    );

    // =========================
    // CREATE LIST OF DATA TABLES
    // =========================
    Rcpp::List loco_result = Rcpp::List::create(
        Rcpp::Named("RawData") = raw_df,
        Rcpp::Named("LaplacianScores") = laplacian_scores,
        Rcpp::Named("Correlations") = corr_mat,
        Rcpp::Named("Neighbourhoods") = n_df,
        Rcpp::Named("NeighbourhoodStatus") = status_df
    );

    return(loco_result);
}

Rcpp::List run_correlation_propagation_across_graph(const SingleCellData& inFile, int thread,
                                              const unsigned int numberNeighbourhoods, const unsigned int neighborhoodSize, 
                                              const int neighborhoodKNN, const double& correlationCutoff,
                                               const std::vector<std::string>& cellStateGenes,
                                              const std::vector<std::string>& corrStateGenes, 
                                              const int permutations, const int minSetSize, const double corrSetAbundance, 
                                              const unsigned int correlatedSetMode, const std::string& correlationType,
                                            const bool calcFeatureSets, const std::string& nSampling,
                                            const double minNonNAFraction, const unsigned int numberCorrelations,
                                            const double significanceCutoff, const double minVaryingFraction,
                                            const double mdavAlpha, const double mdavRadiusFactor,
                                            const unsigned int kmeansR, const double kmeansSlack, const unsigned int seed)
{
    //we can store results for many Neighbourhood-size simultaneously
    Rcpp::List all_results;

    //generate cell-cell neighborhood graph
    std::vector<int> cellStateIdxs = get_indexlist_from_genenames(inFile, cellStateGenes);
    std::vector<int> corrIdxs = get_indexlist_from_genenames(inFile, corrStateGenes);

    // create graph of single-cell data
    // mdav decides the number of neighbourhoods itself (numberNeighbourhoods > 0 is only an upper limit there)
    unsigned int numberNeighbourhoodsCalculated = numberNeighbourhoods;
    if(numberNeighbourhoods == 0 && nSampling != "mdav")
    {
        numberNeighbourhoodsCalculated = inFile.pointCloud.size() / neighborhoodSize;
    }
    LOCO_OUT << "STEP 1\n";
    if(nSampling == "mdav")
    {
        LOCO_OUT << "\tCreating disjoint neighbourhoods with exactly " << neighborhoodSize << " cells (mdav)\n";
    }
    else
    {
        LOCO_OUT << "\tCreating " << numberNeighbourhoodsCalculated << " neighbourhoods with " << neighborhoodSize << " cells\n";
    }

    if(numberNeighbourhoodsCalculated == 1)
    {
        std::cout << "WARNING: YOU ARE RUNNING LoCo WITH A SINGLE NEIGHBOURHOOD!!!\n";
        std::cout << " LoCo still estimated correlations within this neihgbourhood but can not detect smooth changes\n";
    }

    bool printStatusUpdateCellDistCalc = true;
    unsigned int scGraphKnn = neighborhoodSize; //the KNN value is the number of cells in a neighborhood, we ONLY have to calcualte the knn closest neighbors, no need for more
    bool precalculateAllDistances = false;
    std::cout << "Calculate cell distances\n";
    std::shared_ptr<GraphData> scNormData = std::make_shared<GraphData>(inFile, cellStateIdxs, scGraphKnn, &GraphIni::cell_similarity_graph_manhattan_raw, thread, printStatusUpdateCellDistCalc, precalculateAllDistances);

    //create Neighborhoods
    Neighborhood neighborhood(scNormData, numberNeighbourhoodsCalculated, neighborhoodSize, neighborhoodKNN, 
                            inFile, cellStateIdxs, corrIdxs, permutations, corrSetAbundance,
                            correlatedSetMode, correlationType, nSampling, mdavAlpha, mdavRadiusFactor, thread,
                            kmeansR, kmeansSlack, seed);
    neighborhood.calculate_correlation_propagation(correlationCutoff, minSetSize, calcFeatureSets, thread,
                                                  minNonNAFraction, numberCorrelations, significanceCutoff,
                                                  minVaryingFraction);

    //return the RCPP data structure for loco
    Rcpp::List res = build_loco_object(inFile, neighborhood, calcFeatureSets);

    return(res);
}

// [[Rcpp::export]]
Rcpp::List run_loco_cpp(
    std::string inFile,
    char del,
    bool col,
    bool row,
    bool zscore,
    int thread,
    unsigned int correlatedSetMode,
    std::string cellStateGeneFile,
    std::string correlationStateGeneFile,
    unsigned int numberNeighbourhoods,
    unsigned int neighborhoodSize,
    int neighborhoodKNN,
    double correlationCutoff,
    int permutations,
    int minSetSize,
    double corrSetAbundance,
    std::string correlationType,
    bool calcFeatureSets,
    std::string nSampling,
    double minNonNAFraction,
    unsigned int numberCorrelations,
    double significanceCutoff,
    double minVaryingFraction,
    double mdavAlpha,
    double mdavRadiusFactor,
    unsigned int kmeansR,
    double kmeansSlack,
    unsigned int seed
){

    //READ IN DATA
    SCParser parser(inFile, del, col, row);
    SingleCellData inputDataRaw = parser.getData();

    if(zscore)
    {
        LOCO_OUT << "z-score normalize data (scale feature counts for each single-cell)\n";
        zscore_singleCelldata(inputDataRaw);
    }

    //Read in gene lists
    std::vector<std::string> cellStateGenes;
    std::vector<std::string> corrStateGenes;

    if (!cellStateGeneFile.empty())
        cellStateGenes = parse_list(cellStateGeneFile, ',');

    if (!correlationStateGeneFile.empty())
        corrStateGenes = parse_list(correlationStateGeneFile, ',');

    Rcpp::List result = run_correlation_propagation_across_graph(
        inputDataRaw,
        thread,
        numberNeighbourhoods,
        neighborhoodSize,
        neighborhoodKNN,
        correlationCutoff,
        cellStateGenes,
        corrStateGenes,
        permutations,
        minSetSize,
        corrSetAbundance,
        correlatedSetMode,
        correlationType,
        calcFeatureSets,
        nSampling,
        minNonNAFraction,
        numberCorrelations,
        significanceCutoff,
        minVaryingFraction,
        mdavAlpha,
        mdavRadiusFactor,
        kmeansR,
        kmeansSlack,
        seed
    );

    return(result);
}

