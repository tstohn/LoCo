#pragma once

#include <iostream>
#include <memory>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <exception>
#include <numeric>
#include "loco_io.h"

#include "GraphData.hpp"

//graph structure for graphs in which we search for connected components
//this graph is not used for the neighborhoods graph, but only to create a quick/ fast accessable and traversable graph
// of proteins for all neighbourshoods
struct Graph {
    int num_nodes;
    //the starting positions for edges for every node
    std::vector<int> offsets;   // size = num_nodes + 1 because the last node needs and end (the offset of the next node which is just the end of the edges vector)
    // core element, a list of edges like [1,2,1,3], in this list all elements are nodes that the ofdfset node points to
    //these edges are ordered by nodes (for 1 edges, then 2 edges) and one can access the edges for a node with offset
    //so offset[0] until offset[1] are all the nodes that can be reached from node 0
    std::vector<int> edges;     
    std::vector<double> weights; // same as edges but containing the weights
};

/**
 * @brief creates an igraph_t from GraphData class. This is really JUST A BUILDER class. It is not intented to be used as
 *        a graph class itself.
* @param data SCData object that stores a single cell*feature matrix and adjacency matrix with all to all distances
* @param inputRadius radius that is used to include cells for neighborhood
* @param bandwidth bandwidth parameter to calculate an edge strength parameter between neighbors with a gaussian kernel
                   (0: estimate bandwidth by algorithm, -1: do not calcualte kernel but keep origional values)
*/
class GraphHandler
{

    public:
        // knn = 0 creates NO knn graph but a k-dist graph based on radius threshold, and scales the idstances with gaussian kerbel using bandwidth (if non negative)
        // bandwidth==0: estimate badnwidth, <0: skip gaussian kernel distance scaling, keepHighEdges is set if we keep only edges ABOVE a threshold
        GraphHandler(std::shared_ptr<const GraphData> data, int knn = 5, double inputRadius = 0, double bandwidth = 0);

        GraphHandler (const  GraphHandler &)  = delete;

        //Destructor
        ~GraphHandler() 
        {
            //free the weightmatrix that we used to initialize graph
            if(weightedAdjacencyMatrix == nullptr) return;
            // Cleanup the iGraph Constructor mess
            for (size_t i = 0; i < data->number_of_nodes(); i++) 
            {
                free((weightedAdjacencyMatrix)[i]);
            }
            free(weightedAdjacencyMatrix);
        }

        //returns the edge weight between ndoes if they exist, if they r not connected its zero
        double get_edge_weight_between_nodes(const nodePtr& nodeA, const nodePtr& nodeB) const;
        const std::vector<nodePtr> get_all_nodes() const
        {
            return(data->get_all_nodes());
        }
        const nodePtr get_node_at(const int i) const
        {
            return(data->get_node_at(i));
        }
        const std::vector<std::string> get_features() const
        {
            return(data->get_all_feature_names());
        }

        void create_graph();

        void dfs_search(const int start,
                        std::vector<char>& visited,
                        std::vector<int>& component,
                        std::vector<int>& stack,
                        const int minDegree);
        void bron_kerbosch(std::vector<std::vector<int>>& correlationSet,
                                std::vector<int>& R,
                                std::vector<int>& P,
                                std::vector<int>& X,
                                const int minSetSize);
        void calculate_fully_connected_components(std::vector<std::vector<int>>& cliqueVector, const int minCliqueSize);
        void calculate_connected_components(std::vector<std::vector<int>>& correlationSet, const size_t minSetSize);
        void calculate_min_edge_connected_components(std::vector<std::vector<int>>& cliqueVector,
                                           const size_t minCliqueSize, int kcoreThreshold);

        void fill_distance_matrix(double** weightMatrix);
        void fill_knn_matrix(double** weightMatrix);
        double calc_gaussian_kernel(const nodePtr& nodeA, const nodePtr& nodeB, const double& bandwidth) const;

        //getter functions
        const std::shared_ptr<const GraphData> return_data() const
        {
            return data;
        }
        std::vector<int> return_knn_neighbor_nodes(int id) const
        {
            //this only works for KNN graphs, not for K-DIST graphs
            if(knn == 0)
            {
                LOCO_ERR << "Error: Can not calculate neighbors when knn is set to zero";
                LOCO_EXIT(EXIT_FAILURE);
            }

            //get all neighbors from data (fastest way, it has all neighbors ordered)
            nodePtr node = data->get_node_at(id);
            return(data->get_adjacent_node_ids_knn(node, knn));
        }
        const std::shared_ptr<const GraphData> get_data() const
        {
            return data;
        }

    private:

        const int knn; //neighbors in case we construct a knn graph
        const double radius; //radius of nodes to use for graph construction
        std::shared_ptr<const GraphData> data;
        double bandwidth; //bandwidth to use for gaussian kernel on edges: 0 to calculate it, -1 to skip gaussian-kernel application

        //functions
        // calculate badnwidth
        double calc_bandwidth();
        std::vector<int> merge_similar_cliques(const std::vector<std::vector<int>>& cliques);
        double clique_similarity(const std::vector<int>& cliqueA, const std::vector<int>& cliqueB);
        void collapse_all_similar_cliques(std::vector<std::vector<int>>& cliques, double threshold);

        //INDEPENDANT PROCESSES
        std::vector<subSpace> subSpaces;

        //iGraph: nodes and adj matrix in the same order as ndoes in GraphData!!
        Graph graph;
        // adjacency matrix: symmetric 
        double** weightedAdjacencyMatrix = nullptr;
};
typedef std::shared_ptr<GraphHandler> GraphHandlerPtr;

/**
@brief: Class to handle a set of many SCGraphs. It can reconstruct many graphs with varying radii and store them.
@param: GraphData: the input data with distances between all points
@param: inoutRadii: vector of radii to construct graphs from GraphData
**/
class SCGraphSetHandler
{

    SCGraphSetHandler(std::shared_ptr<const GraphData> data, const std::vector<double>& inputRadii, double bandwidth = 0);

    private:
        std::vector<GraphHandlerPtr> graphHandlerVector;
};


