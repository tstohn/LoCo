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
// 3.) correlations
// 4,) neighbourhoods: anchor cell + contained cells
Rcpp::List build_loco_object(const SingleCellData& rawData,
                             Neighborhood& neighborhood,
                             const int& numberCorrelations) 
{

    //fill intermittend data that are used to write the R-result object
    std::vector<std::string> nIDs; //all neighborhoods IDs
    std::vector<std::string> nID_anchorCellID; //vector (for each nID) off all anchor cellIDs (same order as nIDs)
    std::vector<std::vector<std::string>> nID_allCellIDs;// vector (for each nID) of vector of all cell IDs 

    std::vector<std::string> correlation_pairs; //all names of the correlation pairs
    std::vector<std::vector<double>> corrMat; //all correlations: rows = neighborhoods in nIDs, cols = correlationPairs in correlation_pairs


    std::vector<std::string> laplacian_correlation_pairs; //all names of the correlation pairs for laplacian
    std::vector<double> corrL;
    std::vector<double> pCorrL;
    std::vector<std::vector<std::string>> cliquesFlat;
    neighborhood.fill_result_data(
        numberCorrelations,
        nIDs, nID_anchorCellID,nID_allCellIDs,
        correlation_pairs,corrMat,
        laplacian_correlation_pairs,corrL,pCorrL,cliquesFlat);

    // =========================
    // 1. safe raw data table
    // =========================

    int nrow = rawData.pointCloud.size();
    int ncol = rawData.pointCloud[0].size();

    Rcpp::NumericMatrix mat(nrow, ncol);
    for (int i = 0; i < nrow; i++) 
    {
        for (int j = 0; j < ncol; j++) 
        {
            mat(i, j) = rawData.pointCloud[i][j];
        }
    }
    mat.attr("dimnames") = Rcpp::List::create(
    Rcpp::wrap(rawData.cellIDs),  // row names
    Rcpp::wrap(rawData.geneNames)   // column names
    );
    Rcpp::DataFrame raw_df = Rcpp::as<Rcpp::DataFrame>(mat);
    raw_df.push_front(Rcpp::wrap(rawData.cellIDs), "cellID");

    // =========================
    // 2.) Laplacian scores/ p-values/ cliques for all interesting pairs
    // =========================

    int laplacian_size = laplacian_correlation_pairs.size();
    // Safety check (important for robustness)
    if (corrL.size() != laplacian_size || pCorrL.size() != laplacian_size || cliquesFlat.size() != laplacian_size) 
    {
        Rcpp::stop("Input vectors must have the same length");
    }

    // Flatten cliquesFlat into comma-separated strings
    std::vector<std::string> cliquesCollapsed(laplacian_size);

    for (int i = 0; i < laplacian_size; i++) {
        const auto& cliqueVec = cliquesFlat[i];

        std::string combined;
        for (size_t j = 0; j < cliqueVec.size(); j++) 
        {
            combined += cliqueVec[j];
            if (j < cliqueVec.size() - 1) {
                combined += ";";
            }
        }

        cliquesCollapsed[i] = combined;
    }

    // build DataFrame 
    Rcpp::DataFrame laplacian_scores = Rcpp::DataFrame::create(
        Rcpp::Named("FeaturePair") = laplacian_correlation_pairs,
        Rcpp::Named("LaplacianScore") = corrL,
        Rcpp::Named("p_value") = pCorrL,
        Rcpp::Named("FeatureSet") = cliquesCollapsed,
        Rcpp::_["stringsAsFactors"] = false
    );

    // =========================
    // 3. neighbourhoodIDs - correlations
    // =========================

    // Size checks
    int nIDsize = nIDs.size();
    int nPairs = correlation_pairs.size();
    if (corrMat.size() != nIDsize) 
    {
        Rcpp::stop("corrMat must have same length as nIDs");
    }
    for (int i = 0; i < nIDsize; i++) 
    {
        if (corrMat[i].size() != nPairs) 
        {
            Rcpp::stop("Each corrMat row must match number of correlation pairs");
        }
    }

    //  Build correlation per neighbourhood list
    Rcpp::List df_corr;
    // First column: pair names
    df_corr["CorrelationPair"] = correlation_pairs;
    // Each column = one neighborhood
    for (int i = 0; i < nIDsize; i++) {
        std::vector<double> col(nPairs);
        for (int p = 0; p < nPairs; p++) 
        {
            col[p] = corrMat[i][p];  // row=i (neigh), col=p (pair)
        }
        df_corr[nIDs[i]] = col;
    }

    // Convert to DataFrame
    Rcpp::DataFrame corr_df(df_corr);
    corr_df.attr("stringsAsFactors") = false;


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
        Rcpp::Named("NeighborhoodID") = nIDs,
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
        Rcpp::Named("CorrelationMatrix") = corr_df,
        Rcpp::Named("Neighbourhoods") = n_df
    );

    return(loco_result);
}

Rcpp::List run_correlation_propagation_across_graph(const SingleCellData& inFile, int thread,
                                              const unsigned int& neighborhoodSize, 
                                              const int neighborhoodKNN, const double& correlationCutoff,
                                              int& numberCorrelations, const std::vector<std::string>& cellStateGenes,
                                              const std::vector<std::string>& corrStateGenes, 
                                              const int permutations, const int minSetSize, const double corrSetAbundance, 
                                              const unsigned int correlatedSetMode)
{
    //we can store results for many Neighbourhood-size simultaneously
    Rcpp::List all_results;

    //generate cell-cell neighborhood graph
    std::vector<int> cellStateIdxs = get_indexlist_from_genenames(inFile, cellStateGenes);
    std::vector<int> corrIdxs = get_indexlist_from_genenames(inFile, corrStateGenes);

    // create graph of single-cell data
    unsigned int numNeighborhoods = inFile.pointCloud.size() / neighborhoodSize;
    LOCO_OUT << "Creating " << numNeighborhoods << " neighbourhoods with " << neighborhoodSize << " cells\n";
    bool printStatusUpdateCellDistCalc = true;
    unsigned int scGraphKnn = neighborhoodSize; //the KNN value is the number of cells in a neighborhood, we ONLY have to calcualte the knn closest neighbors, no need for more
    bool precalculateAllDistances = false;
    std::shared_ptr<GraphData> scNormData = std::make_shared<GraphData>(inFile, cellStateIdxs, scGraphKnn, &GraphIni::cell_similarity_graph_manhattan_raw, thread, printStatusUpdateCellDistCalc, precalculateAllDistances);

    //create Neighborhoods
    Neighborhood neighborhood(scNormData, numNeighborhoods, neighborhoodSize, neighborhoodKNN, 
                            inFile, cellStateIdxs, corrIdxs, permutations, corrSetAbundance,
                            correlatedSetMode);
    neighborhood.calculate_correlation_propagation(correlationCutoff, minSetSize, thread);

    //return the RCPP data structure for loco
    Rcpp::List res = build_loco_object(inFile, neighborhood, numberCorrelations);

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
    int numberCorrelations,
    std::string cellStateGeneFile,
    std::string correlationStateGeneFile,
    unsigned int neighborhoodSize,
    int neighborhoodKNN,
    double correlationCutoff,
    int permutations,
    int minSetSize,
    double corrSetAbundance
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
        neighborhoodSize,
        neighborhoodKNN,
        correlationCutoff,
        numberCorrelations,
        cellStateGenes,
        corrStateGenes,
        permutations,
        minSetSize,
        corrSetAbundance,
        correlatedSetMode
    );

    return(result);
}

