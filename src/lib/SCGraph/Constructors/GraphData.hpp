#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <map>
#include <exception>
#include <numeric>

#include "Node.hpp"
#include "Utils/kdTreeUtils.hpp"
#include "Parser/SCParser.hpp"
#include "Utils/correlationUtils.hpp"
#include "Utils/generalUtils.hpp"
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>

// vector of positions of the dimensions that contribute to a subspace
typedef std::vector<int> subSpace;

//functions to initialize the Graph
class GraphData;
namespace GraphIni
{
    // cell-cell similarity graph (maybe include euclidean/ manhattan/ cosine later)
    void cell_similarity_graph_euclidean(GraphData* graphData, const SingleCellData& inputData, const std::vector<int>& cellStateGenes, int threads, bool statusUpdate);
    void cell_similarity_graph_euclidean(GraphData* graphData, std::vector< nodePtr> inputNodes, const std::vector<int>& cellStateGenes, int threads, bool statusUpdate);
    void cell_similarity_graph_manhattan_raw(GraphData* graphData, const SingleCellData& inputData, const std::vector<int>& cellStateGenes, int threads, bool statusUpdate);
    void cell_similarity_graph_manhattan_nodes(GraphData* graphData, std::vector< nodePtr> inputNodes, const std::vector<int>& cellStateGenes, int threads, bool statusUpdate);

    // protein correlation graph
    void protein_correlation_graph(GraphData* graphData, const SingleCellData& inputData, const std::vector<int>& cellStateGenes, int threads, bool statusUpdate);    
}

// a map of neighbors to the actual distance, where neighbors are ordered according to insertion
//order is a vector containing the neighbors in descending order, map maps the neighbors to their actual distance
//iterating through it returns const std::pair<const nodePtr, const double> in descending distance order
template<typename K, typename V>
class OrderedHash {

    public:
        std::vector<K> order;
        std::unordered_map<K, V> map;

        //locates an elemnts, throws an error if not accessible
        const V& locate(K key) const 
        {
            auto iter = map.find(key);
            if (iter == map.end()) throw std::out_of_range("key not found");
            return iter->second;
        }

        //access distance to neighbor: return -1 if not present
        const V at(K key) const 
        {
            auto iter = map.find(key);
            if (iter == map.end()) return(-1);
            return iter->second;
        }

        //iterator to iterate through neighbors in ascending distance order
        typedef typename std::vector<K>::const_iterator VecConstIterator;
        class ConstIterator {
            public:
                ConstIterator(VecConstIterator it, const OrderedHash<K, V>* container)
                    : iterator(it), neighborDistances(container) {}

                ConstIterator& operator++() 
                {
                    ++iterator;
                    return *this;
                }

                ConstIterator operator++(int) 
                {
                    ConstIterator tmp = *this;
                    ++(*this);
                    return tmp;
                }

                bool operator==(const ConstIterator& other) const 
                {
                    return iterator == other.iterator;
                }

                bool operator!=(const ConstIterator& other) const 
                {
                    return iterator != other.iterator;
                }

                //content of iterator is: first=nodePtr, second=double
                const std::pair<const K, const V> operator*() const
                {
                    const K& key = *iterator;
                    const V& value = neighborDistances->map.at(key);
                    return std::make_pair(key, value);
                }

                K& neighbor() const
                {
                    return(*iterator);
                }

                const V& distance() const
                {
                    return(neighborDistances->map.at(*iterator));
                }

            private:
                typename std::vector<K>::const_iterator iterator;
                const OrderedHash<K, V>* neighborDistances;
            };

        ConstIterator begin() const
        {
            return ConstIterator(order.cbegin(), this);
        }

        ConstIterator end() const
        {
            return ConstIterator(order.cend(), this);
        }

        std::size_t size() const
        {
            return(order.size());
        }

        void set(K key, V value) 
        {
            auto iter = map.find(key);
            if (iter != map.end()) 
            {
                // no order change, just update value
                iter->second = value;
            }
            order.push_back(key);
            map.insert(std::make_pair(key, value));
        }

        void sort_values(bool descending = false)
        {
            if(!descending)
            {
                std::sort(order.begin(), order.end(), [this](const auto& a,const auto& b){ return map.at(a) < map.at(b); });
            }
            else
            {
                std::sort(order.begin(), order.end(), [this](const auto& a,const auto& b){ return map.at(a) > map.at(b); });
            }
        }

        void erase(K key) 
        {
            order.erase(locate(key).pos_iter);
            map.erase(key);
        }

        std::vector<K> nodes()
        {
            return(order);
        }

        V operator[](K key) const 
        {
            return locate(key);
        }  
};
typedef OrderedHash<nodePtr, double> OrderedNeighborDistanceHash;

/**
 * @brief Data of single cells: cell positions and their distances to each other WITHOUT already declaring connections
(this is the most basic class to store the single cell data from which graphs with connected cells are contructed)
this class is a specific base class for SCGraph, it stores already the nodes without connecting them, and enables quick
construction of SCGraphs with various granularities (connections)
  - vector of all nodes (single cells)
  - KDtree for quick look up of node distances

  optional:
  - adjacency matrix storing all distances between nodes
    Nodes are in same order as in Adjacency Matrix

If knn is >0 GrpahData create a KD tree, when getting adjacent nodes graphData performs a lookup in the KD tree for closest nodes
If knn is 0 graphData perform a bruteForce ALL AGAINST ALL computation of node-node distances and stores them in the adjacency matrix in DESCENDING ORDER!
This can be helpfull if we want to know the exact order of all neighbors... (only of number of nodes is small) 
**/
class GraphData
{
    public:
        GraphData(const SingleCellData& inputData, const std::vector<int>& cellStateGenes, unsigned int knn, void (*initialize_function)(GraphData*, const SingleCellData&, const std::vector<int>&,int , bool ) = &(GraphIni::cell_similarity_graph_manhattan_raw), int threads = 1, bool statusUpdate = false, bool precalcAllDistances = true);
        GraphData(std::vector< nodePtr> inputNodes, const std::vector<int>& cellStateGenes, unsigned int knn, void (*initialize_function)(GraphData*, std::vector< nodePtr>, const std::vector<int>&,int , bool ) = &(GraphIni::cell_similarity_graph_manhattan_nodes), int threads = 1, bool statusUpdate = false, bool precalcAllDistances = true);

        // NEW: function for knn initialization
        void build_knn_adjacency(int k);
        //KD TREE (not constructed by default)
        void create_kd_tree();

        // SETTER FUNCTIONS
        //###################
        void set_node(const std::vector<double>& point, const std::string name)
        {
            nodes.emplace_back(std::make_shared<node>(point, name));
            nodePtrToIdx.insert(std::pair<nodePtr, unsigned long>(nodes.back(), nodes.size() - 1));
            nodeNameToPtr.insert(std::pair<std::string, nodePtr>(name, nodes.back()));
        }
        void add_node(const nodePtr& node)
        {
            nodes.push_back(node);
            nodePtrToIdx.insert(std::pair<nodePtr, unsigned long>(nodes.back(), nodes.size() - 1));
            nodeNameToPtr.insert(std::pair<std::string, nodePtr>(node->get_name(), nodes.back()));
        }
        void add_distance(const int i, const int j, const double dist)
        {
            adjacencyList[nodes.at(i)].set(nodes.at(j), dist);
            adjacencyList[nodes.at(j)].set(nodes.at(i), dist);
        }
        void sort_adjacencyList(bool descending = false)
        {
            for(size_t i=0; i < nodes.size(); ++i)
            {
                adjacencyList.at(nodes.at(i)).sort_values(descending);    
            }
        }
        void set_neighbors_descending()
        {
            nodesOrderedDescending = true;
        }
        void set_feature_names(const std::vector<std::string>& oriFeatureNames)
        {
            featureNames = oriFeatureNames;
        }

        // GETTER FUNCTIONS
        //###################
        //get a vector of all nodes
        const nodePtrVector get_all_nodes() const
        {
            return(nodes);
        }
        //get a certain node
        const nodePtr get_node_at(const size_t i) const
        {
            return(nodes.at(i));
        }
        const nodePtr get_node_from_name(const std::string& name) const
        {
            return(nodeNameToPtr.at(name));
        }
        const std::vector<std::string> get_all_feature_names() const
        {
            return (featureNames);
        }
        unsigned long get_nodeIdx(const nodePtr n) const
        {
            return(nodePtrToIdx.at(n));
        }
        //get a vector of adjacent (node, distance) pairs for a ndoe
        const OrderedNeighborDistanceHash get_adjacent_nodes(nodePtr node) const
        {
            return(adjacencyList.at(node));
        }
        double get_distance_between_nodes(nodePtr nodeA, nodePtr nodeB) const
        {
            return(adjacencyList.at(nodeA).at(nodeB));
        }
        std::vector<int> get_adjacent_node_ids_knn(nodePtr node, size_t knn) const
        {
            if(!distancesPrecalculated)
            {
                std::cerr << "requesting to return distances from adjacency list, but they have not been precalcualted!\n";
                exit(EXIT_FAILURE);
            }

            std::vector<int> adjNodes;
            const OrderedNeighborDistanceHash neighbors = get_adjacent_nodes(node);
            if(knn > neighbors.size())
            {
                std::cerr << "Requesting more neighbors from node than it has!!!\n";
                std::cerr << knn << " requested, but only" << neighbors.size() << " neighbors exist\n";
            }
            size_t count = 0;
            for(const std::pair<nodePtr, const double> neighbor : neighbors)
            {
                if(count == knn) {break;}
                adjNodes.push_back(get_nodeIdx(neighbor.first));
                ++count;
            }

            return(adjNodes);
        }
        std::vector<nodePtr> get_adjacent_nodes_knn(nodePtr node, size_t knn) const
        {
            if(!distancesPrecalculated)
            {
                std::cerr << "requesting to return distances from adjacency list, but they have not been precalcualted!\n";
                exit(EXIT_FAILURE);
            }

            std::vector<nodePtr> adjNodes;
            const OrderedNeighborDistanceHash neighbors = get_adjacent_nodes(node);
            if(knn > neighbors.size())
            {
                std::cerr << "Requesting more neighbors from node than it has!!!\n";
            }
            size_t count = 0;
            for(const std::pair<const nodePtr, const double> neighbor : neighbors)
            {
                if(count == knn) {break;}
                adjNodes.push_back(neighbor.first);
                ++count;
            }

            return(adjNodes);
        }



        std::vector<int> get_adjacent_node_ids_knn_kdsearch(nodePtr node) const
        {
            if(distancesPrecalculated)
            {
                std::cout << "WARNING: requesting to return distances from KD-tree but distances were precalcualted!\n";
            }

            std::vector<int> adjNodes;

            const int dim   = get_node_at(0)->dimensions();
            std::vector<double> query_pt(dim);
            for (int d = 0; d < dim; ++d)
            {
                query_pt[d] = node->value_at(d);
            }

            std::vector<size_t> ret_indexes(graphKnn);
            std::vector<double> out_dists_sqr(graphKnn);

            nanoflann::KNNResultSet<double> resultSet(graphKnn);
            resultSet.init(ret_indexes.data(), out_dists_sqr.data());

            mat_indexPtr->index->findNeighbors(
                resultSet,
                query_pt.data(),
                nanoflann::SearchParams(10)
            );

            for (size_t r = 0; r < graphKnn; ++r)
            {
                size_t j = ret_indexes[r];

                nodePtr node_j = get_node_at(static_cast<int>(j));
                adjNodes.push_back(get_nodeIdx(node_j));
            }

            return(adjNodes);
        }
        std::vector<nodePtr> get_adjacent_nodes_knn_kdsearch(nodePtr node) const
        {
            if(distancesPrecalculated)
            {
                std::cout << "WARNING: requesting to return distances from KD-tree but distances were precalcualted!\n";
            }
            std::vector<nodePtr> adjNodes;
            
            const int dim   = get_node_at(0)->dimensions();
            std::vector<double> query_pt(dim);
            for (int d = 0; d < dim; ++d)
            {
                query_pt[d] = node->value_at(d);
            }

            std::vector<size_t> ret_indexes(graphKnn);
            std::vector<double> out_dists_sqr(graphKnn);

            nanoflann::KNNResultSet<double> resultSet(graphKnn);
            resultSet.init(ret_indexes.data(), out_dists_sqr.data());

            mat_indexPtr->index->findNeighbors(
                resultSet,
                query_pt.data(),
                nanoflann::SearchParams(10)
            );

            for (size_t r = 0; r < graphKnn; ++r)
            {
                size_t j = ret_indexes[r];

                nodePtr node_j = get_node_at(static_cast<int>(j));
                adjNodes.push_back(node_j);
            }

            return(adjNodes);
        }

        const std::unordered_map<nodePtr, OrderedNeighborDistanceHash>& return_adj_list() const
        {
            return adjacencyList;
        }
        size_t number_of_nodes() const
        {
            return(nodes.size());
        }
        unsigned int get_knn() const
        {
            return graphKnn;
        }
        bool distances_precalcualted() const
        {
            return distancesPrecalculated;
        }

        //writes edges between nodes to terminal
        void print_adjacency_by_name() const
        {            
            std::cout << "ADJ LIST:\n" << "_______________\n";
            for(nodePtr nodeTmp : nodes)
            {
                std::cout << "FROM: " << nodeTmp->get_name() << "\n";
                for(const std::pair<const nodePtr, const double> neighbor : get_adjacent_nodes(nodeTmp))
                {
                    std::cout << neighbor.first->get_name() << "(" << std::to_string(neighbor.second) << ")" << " - ";
                }
                std::cout << "\n";
            }
        }
        void print_adjacency_by_order() const 
        {
            for(size_t i = 0; i < nodes.size(); ++i)
            {
                const OrderedNeighborDistanceHash nodeNeighborList = adjacencyList.at(nodes.at(i));
                for(size_t j = 0; j < nodes.size(); ++j)
                {
                    double dist = nodeNeighborList.at(nodes.at(j));
                    std::cout << std::to_string(i) << " " << std::to_string(j) << " " << std::to_string(dist) << "\n"; 
                }
            }
        }

        //writes nodes and their values to terminal
        void print_data()
        {
            std::cout << "NODE VECTORS:\n" << "_______________\n";
            for(nodePtr nodeTmp : nodes)
            {
                std::cout << "#: " << nodeTmp->get_name() << "\n";
                for( const double& dim : nodeTmp->all_values())
                {
                    std::cout << std::to_string(dim) << "\n";
                }
            }
        }

        double get_average_radius()
        {
            double avgRad = 0;
            for(size_t i = 0; i < nodes.size(); ++i)
            {
                const OrderedNeighborDistanceHash nodeNeighborList = adjacencyList.at(nodes.at(i));
                for(size_t j = 0; j < nodes.size(); ++j)
                {
                    double dist = nodeNeighborList.at(nodes.at(j));
                    avgRad += dist;
                }
            }
            return(avgRad/(nodes.size()^2));
        }
        bool is_node_order_descending() const
        {
            return nodesOrderedDescending;
        }

        //search for close points in KD tree
        void search_kd_tree();
        std::vector<std::pair<unsigned long, double> >  get_points_within_radius(node node, double radius); // (node,distance) pairs
        void brute_force_get_points_within_radius(node node, double radius);

    private:

        void perform_pca(const std::vector<std::vector<double>>& pointCloud, 
                         const unsigned int& pcas);
        double calc_bandwidth();

        //kd tree variables
        typedef KDTreeVectorOfNodePtrAdaptor<> my_kd_tree_t;
        std::shared_ptr<my_kd_tree_t> mat_indexPtr;

        //bool beeing set of we use a function for graph creation that lists nodes in descending order
        bool nodesOrderedDescending = false;

        //ALL NODES
        nodePtrVector nodes;
        //mapping nodeName to nodePtr
        std::unordered_map<std::string, nodePtr> nodeNameToPtr; // sometimes you want to access a node by its name (e.g. when in NeighborhoodGraph we have not all nodes and the order does not exist anymore, but names do)
        //name of dimensions
        std::vector<std::string> featureNames;
        //cell state genes: nodes are created only with those
        std::vector<int> stateIndices;

        unsigned int graphKnn;

        //MAP OF NODEPTR TO ITS INDEX IN NODES (ABOVE)
        std::unordered_map<nodePtr, unsigned long> nodePtrToIdx; // we have a list of nodes as they r stored in the adjacency matrix (however in adj matrix nodes are stored in decreasing distance for quick look-ups)
        //ADJACENCY MATRIX: map of nodePtr to a vector of pairs that are <nodePtr, Distance>
        //the vector of neighbor-dist pairs is ordered from close to far neighbors
        std::unordered_map<nodePtr, OrderedNeighborDistanceHash> adjacencyList; 
        //mapping the names of genes/ proteins to the index of this feature in the value vector of every node 
        std::unordered_map<std::string, int> featureNameToIdx;

        // precalcualte all distances between points. 
        // If one uses KD-trees with few lookups on distances it makes sense to NOT precalculate and return distances on request
        bool distancesPrecalculated;
};
typedef std::shared_ptr<GraphData> GraphDataPtr;