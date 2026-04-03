#pragma once

#include <Rcpp.h>
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
        std::string tmpPrefix = prefix;
        if(neighborhoodSizes.size() > 1)
        {
            tmpPrefix = prefix + "_Nsize_" + std::to_string(neighborhoodSize);
        }
        neighborhood.write_results_to_file(outFile, tmpPrefix, numberCorrelations);
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

