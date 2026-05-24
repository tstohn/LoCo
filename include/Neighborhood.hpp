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
    // =========================================================================
    // DEFAULT ACCESS
    // =========================================================================

    // FAST MATRIX-CELL ACCESS
    // The 'const' version is automatically called by C++ when reading from a const matrix.
    // row = pair, col = N
    inline const double& operator()(size_t row, size_t col) const 
    {
        return data[row * cols + col];
    }

    // FAST ACCESS TO ALL CORRS ACROSS N FOR A PAIR: check the col entry to not go too fart
    // Returns a raw pointer to the start of the row. 
    // In C++, pointers can be indexed exactly like vectors (e.g., row_ptr[0], row_ptr[1]).
    inline const double* get_row(size_t row_idx) const 
    {
        if (row_idx >= rows) throw std::out_of_range("Row index out of bounds");
        return &data[row_idx * cols];
    }

    // =========================================================================
    // 2. Additional safe methods/ methods returning vectors
    // =========================================================================

    // FAST ACCESS (Unchecked)
    // Overloading () allows us to use matrix(row, col) syntax.
    // It maps 2D coordinates to a 1D index using: index = (row * total_columns) + col.
    // Returns a reference (&), meaning it acts as BOTH a Getter and a Setter!
    inline double& operator()(size_t row, size_t col) {
        return data[row * cols + col]; // Use [] to skip bounds-checking for blazing speed
    }

    // SAFE ACCESS (Checked - Use this for debugging!)
    // Will throw a loud error if you accidentally request an index out of bounds.
    inline double& at(size_t row, size_t col) {
        if (row >= rows || col >= cols) throw std::out_of_range("Matrix subscript out of bounds");
        return data[row * cols + col];
    }

    // OPTION B: DEEP COPY (Standard std::vector)
    // Extracts the entire row and copies it into a brand new standalone std::vector.
    // Use this if downstream functions explicitly require a 'std::vector<double>' object.
    inline std::vector<double> get_row_copy(size_t row_idx) const {
        if (row_idx >= rows) throw std::out_of_range("Row index out of bounds");
        
        auto row_start = data.begin() + (row_idx * cols);
        auto row_end   = row_start + cols;
        
        return std::vector<double>(row_start, row_end);
    }
};

//HOW TO ITERATE THROUGH THE NIEGHBOURHHODS FOR A CORR-PAIR
/*
const double* row_view = neighbourhoodCorrs.pairToListOfAllNCorrs.get_row_view(2);

// You can loop over this pointer EXACTLY like a vector!
for(size_t n = 0; n < myMat.cols; ++n) {
    double corr_val = row_view[n];
    // ... do math on corr_val ...
}
*/ 

struct CorrelationResults
{
    // --- Data Storage ---
    std::vector<std::pair<int, int>> pairList;   // list of unique feature pairs, the indices are alreayd inidces for the SUBSET of feauters uder for correlations
    std::vector<nodePtr> nodePtrList;           // Master ordered list of subsets (N)

    std::vector<std::string> featureNames; //list of single feature names (of the subset of features used for corr analysis)
    std::vector<std::string> pairNames;      // names of all pairs (e.g., []"A_B", "B_D"]), for quick printing later on

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

    //get the string for a pair
    std::string get_pair_name_string(size_t pair_idx) const
    {
        if (pair_idx >= pairList.size()) {
            throw std::out_of_range("Pair index out of bounds.");
        }

        int feat_a = pairList[pair_idx].first;
        int feat_b = pairList[pair_idx].second;

        return featureNames.at(feat_a) + "_" + featureNames.at(feat_b);
    }

    void init(const FlatMatrix& tempCalculationMatrix, 
              const std::vector<std::pair<int, int>>& allPairs,
              const std::vector<std::atomic<int>>& countsAboveThreshold,
              int minNeighborhoodsRequired,
              const std::vector<nodePtr>& centralNeighborhoodPtrs,
              const std::vector<std::string>& globalFeatureNames,
              const std::vector<int>& corrStateGenes) 
    {
        nodePtrList = centralNeighborhoodPtrs;
        size_t numNeighbourhoods = nodePtrList.size();
        size_t total_original_pairs = allPairs.size();

        // 1. MAP FEATURE NAMES
        size_t geneSize = corrStateGenes.empty() ? globalFeatureNames.size() : corrStateGenes.size();
        featureNames.reserve(geneSize);
        for (size_t i = 0; i < geneSize; ++i) 
        {
            int original_idx = corrStateGenes.empty() ? i : corrStateGenes[i];
            featureNames.push_back(globalFeatureNames[original_idx]);
        }

        // 2. FILTER SURVIVING PAIRS
        std::vector<size_t> originalPairIndices; 
        
        for (size_t p = 0; p < total_original_pairs; ++p) 
        {
            if (countsAboveThreshold[p].load() >= minNeighborhoodsRequired) 
            {
                pairList.push_back(allPairs[p]);      
                originalPairIndices.push_back(p);     
            }
        }

        size_t num_filtered_pairs = pairList.size();
        if (num_filtered_pairs == 0) 
        {
            LOCO_OUT << "Error: 0 feature-pairs passed the required threshold (minimum correlation in number of neighbourhoods)\n";
            pairToListOfAllNCorrs.init(0, 0); 
            LOCO_EXIT(EXIT_FAILURE);
        }

        // 3. PRE-BUILD THE PAIR NAMES VECTOR
        pairNames.reserve(num_filtered_pairs);
        for (const auto& pair : pairList)
        {
            pairNames.push_back(featureNames[pair.first] + "_" + featureNames[pair.second]);
        }

        // 4. TRANSPOSE MATRIX
        pairToListOfAllNCorrs.init(num_filtered_pairs, numNeighbourhoods);

        for (size_t new_p = 0; new_p < num_filtered_pairs; ++new_p) 
        {
            size_t old_p = originalPairIndices[new_p]; 
            
            for (size_t s = 0; s < numNeighbourhoods; ++s) 
            {
                pairToListOfAllNCorrs(new_p, s) = tempCalculationMatrix(s, old_p);
            }
        }
    }
};

struct LaplacianResults
{
    std::vector<std::string> pairNames;      // names of all pairs (e.g., []"A_B", "B_D"]), for quick printing later on

    //pre-calcualte variances, we need them many times when calculating significance
    std::vector<double> variances;
    std::vector<double> L;

    std::vector<double> p_values;
};

//small edge struct for quick shuffling in permutations for laplacian results
struct Edge { int i,j; double w; };


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
        void write_results_to_file(const std::string& output, const std::string& prefix, bool calcSets);
        void write_shuffled_laplacians(const std::string& outFile, const std::string& prefix);
        void fill_result_data(
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
        void laplacian_significance_for_pair(size_t pair_idx, const std::vector<Edge>& edges, std::atomic<int>& currentCount);
        void calculate_laplacian_score_for_pair(const int featurePairIdx);
        void calculate_pair_variance(size_t pair_idx);
        void step_2_calculate_correlation(const double& corrThreshold, const int threads);
        void step_3_calculate_laplacian_score(const int threads);

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
        LaplacianResults laplacianScores;

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




