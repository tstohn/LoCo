#pragma once
#include "loco_io.h"

#include "GraphHandler.hpp"
#include "generalUtils.hpp"
#include "threadPool.hpp"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <map>
#include <exception>
#include <numeric>
#include <vector>
#include <cmath>

struct VectorHash {
    size_t operator()(const std::vector<int>& v) const {
        size_t seed = 0;
        for (int i : v) {
            // "Golden Ratio" hashing to spread bits effectively
            seed ^= std::hash<int>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

// Hash combining function
template <typename T>
inline void hash_combine(std::size_t& seed, const T& value) 
{
    seed ^= std::hash<T>{}(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
struct pair_hash 
{
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const 
    {
        std::size_t hash = 0;
        hash_combine(hash, p.first);
        hash_combine(hash, p.second);
        return hash;
    }
};

// flat vector of vector for blazing fast access and copying
// default entries are -2 (as corrs are between -1 and 1)
struct FlatMatrix {
    size_t rows = 0;
    size_t cols = 0;
    std::vector<double> data;

    // Default Constructor (Starts empty for late-initialization)
    FlatMatrix() = default;

    // Size Constructor
    FlatMatrix(size_t r, size_t c) : rows(r), cols(c), data(r * c, -2.0) {}

    // Late Initialization Method
    void init(size_t r, size_t c) {
        rows = r;
        cols = c;
        data.assign(r * c, -2.0); // Allocates and sets everything to -2.0
    }

    // Overload () for ultra-fast writing/modifying data (Uses [] to avoid bound-check overhead)
    inline double& operator()(size_t row, size_t col) {
        return data[row * cols + col];
    }

    // Overload () for ultra-fast read-only access
    inline const double& operator()(size_t row, size_t col) const {
        return data[row * cols + col];
    }

    // C++17 Row Access: Returns a direct pointer to the start of the row (Zero copy!)
    inline const double* get_row(size_t row_idx) const {
        return &data[row_idx * cols];
    }
};

//HOW TO ITERATE THROUGH THE NIEGHBOURHHODS FOR A CORR-PAIR
/*
size_t num_subsets = 0;

// Call the function. It populates 'num_subsets' with the correct column count
const double* corrs = results.get_correlations_across_N_for_pair(4, num_subsets);

// The caller now uses 'num_subsets' to know when to stop!
for (size_t s = 0; s < num_subsets; ++s) {
    double c = corrs[s];
    // Process correlation...
}
*/

struct CorrelationResults
{
    // --- Data Storage ---
    std::vector<std::pair<int, int>> pairList;   // Global list of unique feature pairs
    std::vector<nodePtr> nodePtrList;           // Master ordered list of subsets (N)
    std::vector<std::string> featureNames;      // Maps a feature integer ID to its string name

    // Downstream Layout: Rows = Feature Pairs, Cols = Subsets (N)
    FlatMatrix pairToListOfAllNCorrs;

    // --- 1. Get Correlations Across N for a specific pair index ---
    // Returns a raw pointer to the start of the contiguous block of subset correlations.
    // Out parameter 'out_count' will tell the downstream loop how many elements exist.
    const double* get_correlations_across_N_for_pair(size_t pair_idx, size_t& out_count) const 
    {
        if (pair_idx >= pairList.size()) {
            throw std::out_of_range("Pair index out of bounds in CorrelationResults.");
        }
        
        out_count = pairToListOfAllNCorrs.cols; // The number of subsets (N)
        return pairToListOfAllNCorrs.get_row(pair_idx);
    }

    // --- 2. Feature Name Lookup Utilities ---
    // Gets the string names for a given feature pair
    std::pair<std::string, std::string> get_feature_names_for_pair(size_t pair_idx) const 
    {
        if (pair_idx >= pairList.size()) {
            throw std::out_of_range("Pair index out of bounds.");
        }
        
        int feat_a = pairList[pair_idx].first;
        int feat_b = pairList[pair_idx].second;
        
        return { featureNames.at(feat_a), featureNames.at(feat_b) };
    }

    // --- 3. Initialization & Blazing-Fast Multithreaded Transpose ---
    // Takes the raw, multi-threaded calculation matrix (Rows = Subsets, Cols = Pairs)
    // and transposes it beautifully into the internal downstream matrix layout.
    void init(const FlatMatrix& calculationMatrix, 
              const std::vector<std::pair<int, int>>& globalPairs,
              const std::vector<nodePtr>& globalNodes,
              const std::vector<std::string>& globalFeatureNames) 
    {
        // Copy lists and mappings
        pairList = globalPairs;
        nodePtrList = globalNodes;
        featureNames = globalFeatureNames;

        size_t num_subsets = nodePtrList.size();
        size_t num_pairs = pairList.size();

        // Safety check to ensure input matches dimensions
        if (calculationMatrix.rows != num_subsets || calculationMatrix.cols != num_pairs) {
            throw std::invalid_argument("Input calculation matrix dimensions do not match provided node/pair sizes.");
        }

        // Initialize our downstream matrix with transposed dimensions (Rows = Pairs, Cols = Subsets)
        pairToListOfAllNCorrs.init(num_pairs, num_subsets);

        // BLAZING FAST PARALLEL TRANSPOSE
        // Loop through pairs on the outer loop. If you use OpenMP, uncomment the line below:
        // #pragma omp parallel for schedule(static)
        for (size_t p = 0; p < num_pairs; ++p) 
        {
            for (size_t s = 0; s < num_subsets; ++s) 
            {
                // Read from calculation matrix (s, p) and write to downstream matrix (p, s)
                // This maximizes cache prefetching during writing!
                pairToListOfAllNCorrs(p, s) = calculationMatrix(s, p);
            }
        }
    }
};



//this is a result for every individual neighborhood
struct CorrelationPropagationResult
{
    //for all protein-protein pairs correlation/ slope in every neighborhood
    //pairs of int: int is the index for a protein. we can not use nodes bcs. they r specific to the neighborhood (different points in the neighborhoods)
    std::unordered_map<const std::pair<int, int>, const double, pair_hash> correlationResult; //pair -> correlation
    std::unordered_map<const std::pair<int, int>, const double, pair_hash> slopeResult; //pair -> slope
};

class Neighborhood
{

    public:
        //
        Neighborhood(const std::shared_ptr<const GraphData> data, unsigned int neighborhoodNumber, unsigned int neighborhoodSize, int neighborhoodKNN,
                     const SingleCellData& inputData,
                     const std::vector<int>& cellStateGenes, const std::vector<int>& corrStateGenes, int permutations,
                     const double& corrSetAbundance, const unsigned int correlatedSetMode,
                    const std::string& correlationType);

        //calculate how correlation cliques of proteins change smoothly along the
        //cell-cell neighborhood graph (from neighborhood to neighborhood)
        //fills the corrResult
        void calculate_correlation_propagation(double correlationStrengthCutoff, int minCliqueSize=2, int thread=5);
        void write_results_to_file(const std::string& output, const std::string& prefix, int& numberCorrelations);
        void write_shuffled_laplacians(const std::string& outFile, const std::string& prefix);
        void fill_result_data(
            const int& numberCorrelations,
            std::vector<std::string>& nIDs, // all neighborhoods IDs
            std::vector<std::string>& nID_anchorCellID, //achnor cell IDs for neighborhoods
            std::vector<std::vector<std::string>>& nID_allCellIDs, //vector off all cellIDs for all neighborhoods (same order as nIDs)
            std::vector<std::string>& correlation_pairs, //all names of the correlation pairs
            std::vector<std::vector<double>>& corrMat, //all correlations
            std::vector<std::string>& laplacian_correlation_pairs, //all names of the correlation pairs for laplacian
            std::vector<double>& corrL, 
            std::vector<double>& pCorrL, 
            std::vector<std::vector<std::string>>& cliquesFlat 
        );
    private:

        // new 4-step functions
        void step_2_calculate_correlation(const double& corrThreshold, const int threads);
        void calculate_correlations_for_N(nodePtr neighborhoodCenter, size_t neighborhood_idx, 
                                                const double& corrThreshold, FlatMatrix& tempCalculationMatrix,
                                                const std::vector<std::pair<int, int>>& tmpAllPairs , 
                                                std::vector<std::atomic<int>>& tmpCorrelationCountAboveThreshold, 
                                                std::atomic<int>& currentCount);

        // randomly select x elements from a vector
        std::vector<int> get_random_elements(int numbers, int maxNum);
        void create_neighborhood_graph(int knn);
        void bfs_enumerate_x_closest_neighborhoods(int x, int nodeID, std::vector<int> neighbors);
        std::vector<int> get_neighborhoodIds_by_distance(const int nHoodID);
        void extract_pairs_from_correlation_sets(std::unordered_map<nodePtr, std::shared_ptr<GraphData>>& neighborhoodCorrelations);
        void calculate_correlations_variance(unsigned int numberNodes, const std::pair<int, int>& correlationpair,
                                    std::unordered_map<const std::pair<int, int>, double, pair_hash>& corrVariance,
                                    int totalCount, double& currentCount);
        void calculate_slopes(unsigned int numberNodes, const std::pair<int, int>& correlationpair,
                            std::unordered_map<const std::pair<int, int>, double, pair_hash>& slopeVariance,
                            int totalCount, double& currentCount);
        void laplacian_score(const std::pair<int, int>& pair,
                             const std::unordered_map<const std::pair<int, int>, double, pair_hash>& corrVariance,
                             const std::unordered_map<const std::pair<int, int>, double, pair_hash>& slopeVariance,
                             int totalCount, double& currentCount);
        void laplacian_significance(const std::pair<int, int>& pair,
                                   const std::unordered_map<const std::pair<int, int>, double, pair_hash>& corrVariance,
                                   const std::unordered_map<const std::pair<int, int>, double, pair_hash>& slopeVariance,
                                   const std::vector<CorrelationPropagationResult>& vectorizedResults,
                                   int totalCount, double& currentCount);
        void calculate_laplacian_score(int threads);
        void filter_consistent_correlation_sets_sota(
                const std::unordered_map<nodePtr, std::vector<std::vector<int>>>& cliquesPerNeighborhood,
                int minCliqueSize);

        const std::vector<std::string> get_feaure_names() const
        {
            return(inputDataOrigional.geneNames);
        }
        void detect_cliques_in_neighborhood(nodePtr neighborhoodCenter, const double& correlationStrengthCutoff, int minCliqueSize,
                                            std::unordered_map<nodePtr, std::vector<std::vector<int>>>& cliquesPerNeighborhood, 
                                            std::unordered_map<nodePtr, std::shared_ptr<GraphData>>& neighborhoodCorrelations,
                                            int totalCount, double& currentCount);
        void filter_cliques_present(nodePtr neighborhoodCenter, std::unordered_map<nodePtr, std::vector<std::vector<int>>>& cliquesPerNeighborhood);
        std::vector<std::pair<int, int>> filter_best_pairs(size_t numberGenes);

        unsigned int neighborhoodSize;
        unsigned int neighbourhoodNum;
        /* A NEIGHBORHOOD is essentially a list of NODE IDs, which are the CENTER NODES that define each neighborhood,
        for all those NODE IDS we then aggregate surrounding neighbors to define their neoghborhood (we use GraphData which gets 
        <knn> nodes that r closest to this center point, to get an overlap of those neighborhoods we create another knn graph an get <knn'> surrounding
        neighborhoods that r connected)
        */
        //node IDs as in the origional data (important for later visualisations to map neighborhood center back to a cell)
        std::vector<nodePtr> centralNeighborhoodPtrs;
        //mapping of neighborhoods to all its internal nodes
        //this one has to stay index (int) and not nodePtr bcs we need to filter raw data to create correlation graph from raw pointClouds
        std::unordered_map< nodePtr, std::vector<int>> neighborhoods;
        //grpah of ONLY the nodes in neighborhoods: BE CAREFUL, the indexes of nodes in this graph ARE NOT ANYMORE the
        //same as in neighborhoodIDs and origional data
        std::shared_ptr<GraphHandler> neighborhoodGraph = nullptr;

        //temporary variable to store the inputFile format from which we generate graphs (TODO: make graph from nodePTr vector to make this step obsolete)
        SingleCellData inputDataOrigional;
        //genes used for cellState/ correlations
        std::vector<int> cellStateGenes;
        std::vector<int> corrStateGenes;

        //RESULTS: dict of neighborhood -> results (correlations etc.)
        std::unordered_map<nodePtr, CorrelationPropagationResult> corrResult;

        CorrelationResults neighbourhoodCorrs;

        //Laplacian scores for correlations/ slopes: 
        //it maps: <NodeName, NodeName> -> laplacian score. int is the index of a feature (protein) in the list of nodes of its Data/ Grpahhandler structure
        // (we have to use idx since every neighborhood has its OWN nodePtr due to different single cells in every neighborhood - however, the ID of the feature is the same)
        std::unordered_map< std::pair<int, int>, const double, pair_hash> correlationLaplacian;
        std::unordered_map< std::pair<int, int>, const double, pair_hash> slopeLaplacian;

        //significance values: simply vectors of laplacians for slope/ corr for every correlation/slope pair
        //this way we can calcualte a p-value for every found corr pair
        //the p-value is calcualted just right before writing
        int permutations;
        std::unordered_map< const std::pair<int, int>, const std::vector<double>, pair_hash> shuffledCorrLaplacians;
        std::unordered_map< const std::pair<int, int>, const std::vector<double>, pair_hash> shuffledSlopeLaplacians;

        std::unordered_map< const std::pair<int, int>, const double, pair_hash> p_corr;
        std::unordered_map< const std::pair<int, int>, const double, pair_hash> p_slope;

        //vector of interesting cliques
        std::vector<std::vector<int>> cliquesVector;

        //STORE ALL UNIQUE PAIRS OF FEATURE-CORRELATIONS (FROM ALL EXISTING CLIQUES) - vector of int-pairs (ordered)
        //int is the index of a feature as it was read from raw file, and its order in Neighborhood similarity graph
        //this can also be used to access the laplacians
        std::vector<std::pair<int, int>> pairs;
        //replace this maybe later, for now store in which cliues we observe correlations
        std::unordered_map<const std::pair<int, int>, std::vector<std::vector<std::string>>, pair_hash> pairToClique;
        double minimumCorrSetAbundance;
        unsigned int correlatedSetMode; //0 for finding a connected component (any path - even sparsly connected)
        //  / 1 for a whole clique (all connected)
        // 2 for a conected compoennt with min edges (threshold set to 2 at the moment)

        //spearman or pearson for the correlations to detect
        std::string correlationType;

        //thread variables
        std::mutex threadLock;
};




