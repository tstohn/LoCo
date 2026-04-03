#pragma once

#include <Rcpp.h>
using namespace Rcpp;

// undefine PI as both R and nanoflann define it
#undef PI 
#include "loco_R.h"

// include your existing headers
#include "Neighborhood.hpp"
#include "GraphData.hpp"
#include "SCParser.hpp"

std::vector<unsigned int> parseNeighborhoodSizes(const std::string& input)
{
    std::vector<unsigned int> result;

    std::string s = input;
    s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());

    if (s.empty())
        throw std::invalid_argument("Empty input for neighborhood size");

    if (s.front() == '[' && s.back() == ']')
    {
        s = s.substr(1, s.size() - 2);
        std::stringstream ss(s);
        std::string item;

        while (std::getline(ss, item, ','))
        {
            try {
                result.push_back(static_cast<unsigned int>(std::stoul(item)));
            } catch (const std::exception&) {
                throw std::invalid_argument("Invalid number in list: '" + item + "'");
            }
        }

        if (result.empty())
            throw std::invalid_argument("No numbers found inside brackets");
    }
    else
    {
        try {
            result.push_back(static_cast<unsigned int>(std::stoul(s)));
        } catch (const std::exception&) {
            throw std::invalid_argument("Invalid single number: '" + s + "'");
        }
    }

    return result;
}

//return R-object of LoCo results
Rcpp::List build_loco_object_fast(
    const std::vector<std::pair<int,int>>& pairs,
    const std::vector<std::string>& geneNames,
    const std::vector<int>& corrStateGenes,
    const std::vector<std::vector<double>>& corrMat,
    const std::vector<std::vector<double>>& slopeMat,
    const std::vector<std::vector<int>>& cellMat,
    const std::vector<std::vector<double>>& coords,
    const std::vector<std::string>& featureNames,
    const std::vector<std::string>& neighborhoodNames,
    const std::vector<double>& corrL,
    const std::vector<double>& pCorrL,
    const std::vector<double>& slopeL,
    const std::vector<double>& pSlopeL,
    const std::vector<std::string>& cliquesFlat,
    bool useCorrStateGenes
) {

    size_t P = pairs.size();
    size_t N = corrMat.size();
    size_t F = featureNames.size();

    // =========================
    // 1. LAPPLACIAN (preallocated vectors)
    // =========================
    std::vector<std::string> pairNames;
    std::vector<double> corrScore(P), pCorr(P), slopeScore(P), pSlope(P), ;
    std::vector<std::string> cliqueOut;

    pairNames.reserve(P);
    cliqueOut.reserve(P);

    for (size_t i = 0; i < P; i++) {

        int a = pairs[i].first;
        int b = pairs[i].second;

        if (useCorrStateGenes) {
            a = corrStateGenes[a];
            b = corrStateGenes[b];
        }

        pairNames.push_back(geneNames[a] + "_" + geneNames[b]);

        corrScore[i] = corrL[i];
        pCorr[i] = pCorrL[i];
        slopeScore[i] = slopeL[i];
        pSlope[i] = pSlopeL[i];
        cliqueOut.push_back(cliquesFlat[i]);
    }

    Rcpp::DataFrame laplacian = Rcpp::DataFrame::create(
        Named("ProteinPair") = pairNames,
        Named("CorrelationScore") = corrScore,
        Named("p_CorrelationScore") = pCorr,
        Named("SlopeScore") = slopeScore,
        Named("p_SlopeScore") = pSlope,
        Named("OrigionalCliques") = cliqueOut
    );

    // =========================
    // 2. CORRELATION MATRIX (NUMERIC MATRIX = FASTER)
    // =========================
    Rcpp::NumericMatrix corr(N, P);

    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < P; j++)
            corr(i, j) = corrMat[i][j];

    colnames(corr) = pairNames;

    // =========================
    // 3. SLOPE MATRIX
    // =========================
    Rcpp::NumericMatrix slope(N, P);

    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < P; j++)
            slope(i, j) = slopeMat[i][j];

    colnames(slope) = pairNames;

    // =========================
    // 4. COORDINATES (matrix not vector-of-vectors)
    // =========================
    Rcpp::NumericMatrix coord(coords.size(), F);

    for (size_t i = 0; i < coords.size(); i++)
        for (size_t j = 0; j < F; j++)
            coord(i, j) = coords[i][j];

    colnames(coord) = featureNames;
    rownames(coord) = neighborhoodNames;

    // =========================
    // 5. CELL MATRIX
    // =========================
    size_t C = cellMat[0].size();
    Rcpp::NumericMatrix cells(N, C);

    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < C; j++)
            cells(i, j) = cellMat[i][j];

    colnames(cells) = neighborhoodNames;

    // =========================
    // FINAL OBJECT
    // =========================
    return Rcpp::List::create(
        Named("laplacian") = laplacian,
        Named("correlations") = corr,
        Named("slopes") = slope,
        Named("neighborhoods") = coord,
        Named("cells") = cells
    );
}

void run_correlation_propagation_across_graph(const SingleCellData& inFile, const std::string& outFile, std::string& prefix, int thread,
                                              const std::vector<unsigned int>& neighborhoodSizes, 
                                              const int neighborhoodKNN, const double& correlationCutoff,
                                              int& numberCorrelations, const std::vector<std::string>& cellStateGenes,
                                              const std::vector<std::string>& corrStateGenes, 
                                              const int permutations, const int minSetSize, const double corrSetAbundance, 
                                              const unsigned int correlatedSetMode)
{
    //generate cell-cell neighborhood graph
    std::vector<int> cellStateIdxs = get_indexlist_from_genenames(inFile, cellStateGenes);
    std::vector<int> corrIdxs = get_indexlist_from_genenames(inFile, corrStateGenes);

    for(const unsigned int neighborhoodSize : neighborhoodSizes)
    {
        // create graph of single-cell data
        unsigned int numNeighborhoods = inFile.pointCloud.size() / neighborhoodSize;
        std::cout << "Creating " << numNeighborhoods << " neighbourhoods with " << neighborhoodSize << " cells\n";
        bool printStatusUpdateCellDistCalc = true;
        unsigned int scGraphKnn = neighborhoodSize; //the KNN value is the number of cells in a neighborhood, we ONLY have to calcualte the knn closest neighbors, no need for more
        bool precalculateAllDistances = false;
        std::shared_ptr<GraphData> scNormData = std::make_shared<GraphData>(inFile, cellStateIdxs, scGraphKnn, &GraphIni::cell_similarity_graph_manhattan_raw, thread, printStatusUpdateCellDistCalc, precalculateAllDistances);

        //create Neighborhoods
        Neighborhood neighborhood(scNormData, numNeighborhoods, neighborhoodSize, neighborhoodKNN, 
                                inFile, cellStateIdxs, corrIdxs, permutations, corrSetAbundance,
                                correlatedSetMode);
        neighborhood.calculate_correlation_propagation(correlationCutoff, minSetSize, thread);

        //write results to file:
        //neighborhood, coordinates, correlation, slope for every protein-pair
        //file for protein-pair to origional clique

        neighborhood.write_results_to_file(outFile, tmpPrefix, numberCorrelations);

        //replace this with RCPP data structes that loco should return
        Rcpp::List res = build_loco_object_fast(
            neighborhood.get_pairs(),
            inFile.geneNames,
            corrIdxs,
            neighborhood.get_corr_matrix(),
            neighborhood.get_slope_matrix(),
            neighborhood.get_cell_matrix(),
            neighborhood.get_coords(),
            neighborhood.get_feature_names(),
            neighborhood.get_neighborhood_names(),
            neighborhood.get_corr_laplacian(),
            neighborhood.get_p_corr(),
            neighborhood.get_slope_laplacian(),
            neighborhood.get_p_slope(),
            neighborhood.get_cliques_flat(),
            true
        );
    }
}

void run_loco(
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
    std::string neighborhoodSizeString,
    int neighborhoodKNN,
    double correlationCutoff,
    int permutations,
    int minSetSize,
    double corrSetAbundance
){

    std::vector<unsigned int> neighborhoodSizes = parseNeighborhoodSizes(neighborhoodSizeString);

    //READ IN DATA
    SCParser parser(inFile, del, col, row);
    SingleCellData inputDataRaw = parser.getData();

    if(zscore)
    {
        std::cout << "z-score normalize data (scale feature counts for each single-cell)\n";
        zscore_singleCelldata(inputDataRaw);
    }

    //Read in gene lists
    std::vector<std::string> cellStateGenes;
    std::vector<std::string> corrStateGenes;

    if (!cellStateGeneFile.empty())
        cellStateGenes = parse_list(cellStateGeneFile, ',');

    if (!correlationStateGeneFile.empty())
        corrStateGenes = parse_list(correlationStateGeneFile, ',');

    run_correlation_propagation_across_graph(
        inputDataRaw,
        outFile,
        prefix,
        thread,
        neighborhoodSizes,
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
}

