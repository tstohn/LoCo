#include <Rcpp.h>
using namespace Rcpp;

//redefine commands for cran
#include "loco_io.h"

// undefine PI as both R and nanoflann define it
#undef PI 
#include "loco_R.h"

// include your existing headers
#include "Neighborhood.hpp"
#include "GraphData.hpp"
#include "SCParser.hpp"

//return R-object of LoCo results
// create a List of several dataframes
// 1.a) raw data table

// 2.a) neighbourhood - cells
// 2.b) laplacian scores

// 3a) neighbourhood - coords
// 3b) neighbourhood correlations data table
Rcpp::List build_loco_object(const SingleCellData& rawData,
                             Neighborhood& neighborhood,
                             const int& numberCorrelations) 
{

    //fill intermittend data that are used to write the R-result object
    std::vector<std::string> nIDs; //all neighborhoods IDs
    std::vector<std::string> nID_anchorCellID; //vector (for each nID) off all anchor cellIDs (same order as nIDs)
    std::vector<std::vector<std::string>> nID_allCellIDs;// vector (for each nID) of vector of all cell IDs 

    std::vector<std::string> correlation_pairs; //all names of the correlation pairs
    std::vector<std::vector<double>> corrMat; //all correlations

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

    std::cout << nrow << " " << ncol << " " << rawData.geneNames.size() << " " << rawData.cellIDs.size() << "\n";
    Rcpp::NumericMatrix mat(nrow, ncol);
    for (int i = 0; i < nrow; i++) 
    {
        for (int j = 0; j < ncol; j++) 
        {
            mat(i, j) = rawData.pointCloud[i][j];
        }
    }
    mat.attr("dimnames") = Rcpp::List::create(
    Rcpp::wrap(rawData.geneNames),  // row names
    Rcpp::wrap(rawData.cellIDs)     // column names
    );
    Rcpp::DataFrame raw_df = Rcpp::as<Rcpp::DataFrame>(mat);
    raw_df.push_front(Rcpp::wrap(rawData.geneNames), "gene");

    // =========================
    // 2. safe neighbourhoodIDs - cells
    // and laplacian scores
    // =========================


    // =========================
    // CREATE LIST OF DATA TABLES
    // =========================
    Rcpp::DataFrame loco_result = Rcpp::DataFrame::create(
        Rcpp::Named("RawData") = raw_df

    );

    return(loco_result);
}

Rcpp::List run_correlation_propagation_across_graph(const SingleCellData& inFile, const std::string& outFile, std::string& prefix, int thread,
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
    std::string outFile,
    std::string prefix,
    char del,
    bool col,
    bool row,
    bool zscore,
    int thread,
    unsigned int correlatedSetMode,
    int numberCorrelations,
    std::string cellStateGeneFile,
    std::string correlationStateGeneFile,
    unsigned int numNeighborhoods,
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
        outFile,
        prefix,
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

