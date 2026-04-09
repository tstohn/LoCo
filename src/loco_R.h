#pragma once

#include <string>

// main entry exposed to R
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
);