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
// 3.) correlations: tibble with 3 columns "CorrelationPair", "NeighbourhoodID", "Correlation"
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
    std::vector<std::string>& correlation_pairs_origional; //all names of the correlation pairs for amtrix

    std::vector<std::vector<double>> corrMat; //all correlations: rows = neighborhoods in nIDs, cols = correlationPairs in correlation_pairs

    std::vector<double> corrL;
    std::vector<double> pCorrL;
    std::vector<std::string> featureSetStrings; // a string for each feature-pair listing all sets with comma and semicolon like A,B,C;A,B,D
    neighborhood.fill_result_data(
        nIDs, nID_anchorCellID,nID_allCellIDs,
        correlation_pairs_laplacian, correlation_pairs_origional,corrMat,
        corrL,pCorrL,featureSetStrings, calcFeatureSets);

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
    int laplacian_size = correlation_pairs.size();
    // Safety check (important for robustness)
    if (corrL.size() != laplacian_size || pCorrL.size() != laplacian_size) 
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
            Rcpp::_["stringsAsFactors"] = false
        );
    }
    else
    {
        laplacian_scores = Rcpp::DataFrame::create(
            Rcpp::Named("FeaturePair") = correlation_pairs_laplacian,
            Rcpp::Named("LaplacianScore") = corrL,
            Rcpp::Named("p_value") = pCorrL,
            Rcpp::_["stringsAsFactors"] = false
        );
    }

    // =========================
    // 3. neighbourhoodIDs - correlations (LONG FORMAT)
    // =========================

    // Size checks
    int nIDsize = nIDs.size();
    int nPairs  = correlation_pairs_origional.size();
    int total   = nIDsize * nPairs; // Calculate exact final size

    if (corrMat.size() != nIDsize) {
        Rcpp::stop("corrMat must have same length as nIDs");
    }
    for (int i = 0; i < nIDsize; i++) {
        if (corrMat[i].size() != nPairs) {
            Rcpp::stop("Each corrMat row must match number of correlation pairs");
        }
    }

    // Pre-allocate the exact memory size UP FRONT
    Rcpp::CharacterVector out_pair(total);
    Rcpp::CharacterVector out_neigh(total);
    Rcpp::NumericVector    out_corr(total);

    int idx = 0; // Track the current flat index

    // Build LONG table efficiently using indices
    for (int i = 0; i < nIDsize; i++) {
        const std::string& neigh_id = nIDs[i];

        for (int p = 0; p < nPairs; p++) {
            out_pair[idx]  = correlation_pairs_origional[p];
            out_neigh[idx] = neigh_id;
            out_corr[idx]  = corrMat[i][p];
            idx++; // Advance flat index
        }
    }

    // Build DataFrame
    Rcpp::DataFrame corr_df = Rcpp::DataFrame::create(
        Rcpp::Named("CorrelationPair")   = out_pair,
        Rcpp::Named("NeighbourhoodID")   = out_neigh,
        Rcpp::Named("Correlation")       = out_corr,
        Rcpp::_["stringsAsFactors"]      = false
    );

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
    // CREATE LIST OF DATA TABLES
    // =========================
    Rcpp::List loco_result = Rcpp::List::create(
        Rcpp::Named("RawData") = raw_df,
        Rcpp::Named("LaplacianScores") = laplacian_scores,
        Rcpp::Named("Correlations") = corr_df,
        Rcpp::Named("Neighbourhoods") = n_df
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
                                            const bool calcFeatureSets)
{
    //we can store results for many Neighbourhood-size simultaneously
    Rcpp::List all_results;

    //generate cell-cell neighborhood graph
    std::vector<int> cellStateIdxs = get_indexlist_from_genenames(inFile, cellStateGenes);
    std::vector<int> corrIdxs = get_indexlist_from_genenames(inFile, corrStateGenes);

    // create graph of single-cell data
    unsigned int numberNeighbourhoodsCalculated = numberNeighbourhoods;
    if(numberNeighbourhoods == 0)
    {
        numberNeighbourhoodsCalculated = inFile.pointCloud.size() / neighborhoodSize;
    }
    LOCO_OUT << "STEP 1\n";
    LOCO_OUT << "\tCreating " << numberNeighbourhoodsCalculated << " neighbourhoods with " << neighborhoodSize << " cells\n";
    bool printStatusUpdateCellDistCalc = true;
    unsigned int scGraphKnn = neighborhoodSize; //the KNN value is the number of cells in a neighborhood, we ONLY have to calcualte the knn closest neighbors, no need for more
    bool precalculateAllDistances = false;
    std::cout << "Calculate cell distances\n";
    std::shared_ptr<GraphData> scNormData = std::make_shared<GraphData>(inFile, cellStateIdxs, scGraphKnn, &GraphIni::cell_similarity_graph_manhattan_raw, thread, printStatusUpdateCellDistCalc, precalculateAllDistances);

    //create Neighborhoods
    Neighborhood neighborhood(scNormData, numberNeighbourhoodsCalculated, neighborhoodSize, neighborhoodKNN, 
                            inFile, cellStateIdxs, corrIdxs, permutations, corrSetAbundance,
                            correlatedSetMode, correlationType);
    neighborhood.calculate_correlation_propagation(correlationCutoff, minSetSize, calcFeatureSets, thread);

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
    bool calcFeatureSets
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
        calcFeatureSets
    );

    return(result);
}

