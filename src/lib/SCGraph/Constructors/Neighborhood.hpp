#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <map>
#include <exception>
#include <numeric>
#include <vector>
#include <thread>
#include <pthread.h>
#include <cmath>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>

#include "SCGraph/Constructors/GraphHandler.hpp"
#include "Utils/generalUtils.hpp"

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
                     const double& corrSetAbundance, const unsigned int correlatedSetMode);

        //calculate how correlation cliques of proteins change smoothly along the
        //cell-cell neighborhood graph (from neighborhood to neighborhood)
        //fills the corrResult
        void calculate_correlation_propagation(double correlationStrengthCutoff, const bool printFoundCliquesPerNeighborhood, int minCliqueSize=2, int thread=5);
        void write_results_to_file(const std::string& output, const std::string& prefix, int& numberCorrelations);
        void write_shuffled_laplacians(const std::string& outFile);

    private:

        // randomly select x elements from a vector
        std::vector<int> get_random_elements(int numbers, int maxNum);
        void create_neighborhood_graph(int knn);
        void bfs_enumerate_x_closest_neighborhoods(int x, int nodeID, std::vector<int> neighbors);
        std::vector<int> get_neighborhoodIds_by_distance(const int nHoodID);
        void extract_pairs_from_correlation_sets(std::unordered_map<nodePtr, std::shared_ptr<GraphData>>& neighborhoodCorrelations);
        void calculate_correlations(unsigned int numberNodes, bool print, const std::pair<int, int>& correlationpair,
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
        void calculate_laplacian_score(bool print, int threads);
        void filter_consistent_correlation_sets_sota(
                const std::unordered_map<nodePtr, std::vector<std::vector<int>>>& cliquesPerNeighborhood,
                int minCliqueSize);

        const std::vector<std::string> get_feaure_names() const
        {
            return(inputDataOrigional.geneNames);
        }
        void detect_cliques_in_neighborhood(nodePtr neighborhoodCenter, const double& correlationStrengthCutoff, int minCliqueSize, const bool printFoundCliquesPerNeighborhood,
                                            std::unordered_map<nodePtr, std::vector<std::vector<int>>>& cliquesPerNeighborhood, 
                                            std::unordered_map<nodePtr, std::shared_ptr<GraphData>>& neighborhoodCorrelations,
                                            int totalCount, double& currentCount);
        void filter_cliques_present(nodePtr neighborhoodCenter, std::unordered_map<nodePtr, std::vector<std::vector<int>>>& cliquesPerNeighborhood);
        std::vector<std::pair<int, int>> filter_best_pairs(size_t numberGenes);

        unsigned int neighborhoodSize;
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

        //thread variables
        std::mutex threadLock;
};




