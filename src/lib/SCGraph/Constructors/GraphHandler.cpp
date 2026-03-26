#include <queue>

#include "GraphHandler.hpp"

double GraphHandler::get_edge_weight_between_nodes(const nodePtr& nodeA, const nodePtr& nodeB) const
{
    double edgeweight = 0;
    unsigned long i = data->get_nodeIdx(nodeA);
    unsigned long j = data->get_nodeIdx(nodeB);
    edgeweight = weightedAdjacencyMatrix[i][j];

    return(edgeweight);
}

double GraphHandler::clique_similarity(const std::vector<int>& cliqueA, const std::vector<int>& cliqueB)
{
    size_t cliqueAIdx = 0;
    size_t cliqueBIdx = 0;
    double score = 0;

    //count number of same elements
    while(cliqueBIdx < cliqueB.size())
    {
        //get next element of cliqueA if element in B is higher (if we reach end of A got to next element of cliqueB)
        while(  (cliqueAIdx < cliqueA.size()) && (cliqueB.at(cliqueBIdx) >= cliqueA.at(cliqueAIdx)) )
        {
            if(cliqueB.at(cliqueBIdx) == cliqueA.at(cliqueAIdx))
            {
                ++score;
                ++cliqueAIdx;
                break;
            }
            ++cliqueAIdx;
        }
        ++cliqueBIdx;
    }

    return((score/ std::min(cliqueA.size(), cliqueB.size())));
}

//combines cliques, if they overlap to x percent
std::vector<int> GraphHandler::merge_similar_cliques(const std::vector<std::vector<int>>& cliques)
{
   
   std::vector<int> mergedClique;
    //for both cliques: keep all proteins, and keep them only once
    for(std::vector<int> clique : cliques)
    {
        for(int protein : clique)
        {
            if(std::find(mergedClique.begin(), mergedClique.end(), protein) == mergedClique.end() )
            {
                mergedClique.push_back(protein);
            }
        }
    }
    std::sort(mergedClique.begin(), mergedClique.end());

    return(mergedClique);
}

void GraphHandler::collapse_all_similar_cliques(std::vector<std::vector<int>>& cliques, double threshold)
{

    //queue of cliques: we go through all of them and remove cliques when they r merged with another one
    std::queue<std::vector<int>> cliqueQueue;
    for(const std::vector<int>& clique : cliques)
    {
        cliqueQueue.push(clique);
    }
    //final results vector, it contins all cliques in the beginning, when cliques r merge, the cliues are removed and the merge is inserted
    std::vector<std::vector<int>> finalMergedCliques = cliques;

    //process all elements in this queue and combine them if similar
    //TODO: make queue a unordered_set, and remove all elements that were compared together
    while (!cliqueQueue.empty()) 
    {
        // Get the front element
        std::vector<int> cliqueInQuestion = cliqueQueue.front();
        std::vector<std::vector<int>> cliquesToCombine;
        for(std::vector<int> cliqueB : finalMergedCliques)
        {
            double similarity = clique_similarity(cliqueInQuestion, cliqueB);
            if(similarity >= threshold)
            {
                cliquesToCombine.push_back(cliqueB);
            }
        }
        if(cliquesToCombine.size() > 1) //it will for sure contain itself
        {
            //merge all similar cliques
            std::vector<int> mergedCliques = merge_similar_cliques(cliquesToCombine);
            //delete all cliques that were merged
            for (const std::vector<int>& element : cliquesToCombine) 
            {
                finalMergedCliques.erase(std::remove(finalMergedCliques.begin(), finalMergedCliques.end(), element), finalMergedCliques.end());
            }
            //insert the merged clique: it will be available for next comparison
            //it is no problem that this clique grows: we will comapre it to a smaller claique that can be contained to x-perc in this big one...
            finalMergedCliques.push_back(mergedCliques);
        }

        // Pop the front element for sure
        cliqueQueue.pop();
    }

    cliques = finalMergedCliques;
}

void GraphHandler::fill_distance_matrix(double** weightMatrix)
{
    bool descending = data->is_node_order_descending();
    if( (descending == true) && bandwidth != -1)
    {
        std::cerr << "When nodes are in descending order we can not apply gaussian kernel, we assumes this is use for protein \
        correclation graphs, where only un-scaled values make sense \n";
        std::exit(EXIT_FAILURE);
    }

    int numberNodes = data->number_of_nodes();
    //initialize matrix with zeroes
    for (int i = 0; i < numberNodes; i++) 
    {
        for(int  j = 0; j < numberNodes; j++)
        {
            weightMatrix[i][j] = 0;
        }
    }

    //calculate weight for every edge and insert it into new weighted adjacency matricx
    //we need to keep the initial order of nodes as in nodePtrVector...
    nodePtrVector nodes = data->get_all_nodes();
    for(size_t i = 0; i < nodes.size(); ++i)
    {
        //iterate through neighbors in ascending distance order (as soon as a distance exceeds threshood we can stop...)
        const OrderedNeighborDistanceHash neighbors = data->get_adjacent_nodes(nodes.at(i));
        for(const std::pair<const nodePtr, const double> neighbor : neighbors)
        {
            double weightedDist = 0;
            //as soon as we reach a nieghbor with a distance above threshold, all coming neighbors will be futher away as well
            //this is only done when there is a maximum edge strength threshold
            if( (descending == false) && (neighbor.second > radius)) {break;}
            else if( (descending == true) && (neighbor.second <= radius)) {break;}

            //only include distances that are below the radius
            //in case we keep only SMALLER edges (like cell-cell similiarty graphs)
            if(descending == false)
            {
                if( (neighbor.second >= 0.0) && (neighbor.second <= radius)) //dist might be -1 (default distance in GraphData for non existing edges like self-edges)
                {
                    //APPLY GAUSSIAN KERNEL
                    if(bandwidth > 0.0)
                    {
                        weightedDist = calc_gaussian_kernel(nodes.at(i), neighbor.first, bandwidth);
                    }
                    //do not further change weight
                    else
                    {
                        weightedDist = neighbor.second;
                    }
                }
            }
            else if(descending == true)
            {
                if(neighbor.second > radius)
                {
                    weightedDist = neighbor.second;
                }
            }
            //add to matrix
            unsigned long j = data->get_nodeIdx(neighbor.first);
            weightMatrix[i][j] = weightedDist;
            weightMatrix[j][i] = weightedDist;
        }
    }
}

void GraphHandler::fill_knn_matrix(double** weightMatrix)
{
    if(data->is_node_order_descending() == true)
    {
        std::cerr << "For KNN graphs we must order nodes in ascending order \n";
        std::exit(EXIT_FAILURE);
    }
    int numberNodes = data->number_of_nodes();
    //initialize matrix with zeroes
    for (int i = 0; i < numberNodes; i++) 
    {
        for(int  j = 0; j < numberNodes; j++)
        {
            weightMatrix[i][j] = 0;
        }
    }

    //go through closest neighbors, add them to adjacency matrix (be careful to set right index for nodePtrs)
    nodePtrVector nodes = data->get_all_nodes();
    for(size_t i = 0; i < nodes.size(); ++i)
    {
        const OrderedNeighborDistanceHash neighbors = data->get_adjacent_nodes(nodes.at(i));

        int edges = 0; //only include K nearest neighbors
        for(const std::pair<const nodePtr, const double> neighbor : neighbors)
        {
            if(edges == knn) break;
            double weightedDist = 0;            
            //APPLY GAUSSIAN KERNEL: bandwidth was already estimated in case it was initially zero
            if(bandwidth > 0.0)
            {
                weightedDist = calc_gaussian_kernel(nodes.at(i), neighbor.first, bandwidth);
            }
            //do not further change weight
            else
            {
                weightedDist = neighbor.second;
            }
            edges++;
            
            //add to matrix
            unsigned long j = data->get_nodeIdx(neighbor.first);
            weightMatrix[i][j] = weightedDist;
            weightMatrix[j][i] = weightedDist;
        }
    }
}


//estimating a bandwidth
// 1.) take 10% of closest neighbors
// 2.) calculate avgDist = average distance
// 3.) bandwidth = avgDist/ 2
double GraphHandler::calc_bandwidth()
{
    std::vector<double> distVector;
    unsigned long long closeNeighborSize = data->return_adj_list().size() / 10;

    if(closeNeighborSize == 0)
    {
        throw std::invalid_argument("Number of cells too small. In Bandwidth detection for edge-width estimation by gaussian kernel we selected 10\% of total number of cells which euquals zer0!");
    }
    // for every node
    for(const std::pair<const nodePtr, OrderedNeighborDistanceHash>& neighbor : data->return_adj_list())
    {
        //calculate average distance of ten perc. closest neighbors
        OrderedNeighborDistanceHash::ConstIterator it = neighbor.second.begin();
        double topTenDist = 0;
        for(size_t i =0; i < closeNeighborSize; ++i)
        {
            topTenDist += it.distance();
            ++it;
        }
        distVector.emplace_back(topTenDist/10);
    }

    //calculate average of distances for all nodes, divided by 5
    return( (accumulate(distVector.begin(), distVector.end(), 0))/2 );
}

double GraphHandler::calc_gaussian_kernel(const nodePtr& nodeA, const nodePtr& nodeB, const double& bandwidth) const
{
    double numerator = std::pow( nodeA->distance_to(nodeB), 2);
    double denominator = 2 * std::pow(bandwidth, 2);
    return( std::exp( -(numerator/denominator) ));
}

// TODO: we still do not calcualte new matrix after aplying gaussian kernel to it
// it should be a fully connected adj matrix: stpring also zeroes, se we can use it to initialize the iGraph
GraphHandler::GraphHandler(std::shared_ptr<const GraphData> data, int knn, double inputRadius, double bandwidth) :
            knn{knn}, radius{inputRadius}, data{data}, bandwidth{bandwidth}
{
    //create adjacency matrix with gaussian kernel
    //estimate a bandwidth
    if(knn == 0 && bandwidth == 0){bandwidth = calc_bandwidth();}

    //create the weightedAdjacencyMatrix used for iGraoph creation (ATTENTION: this data needs to be freed upon destruction)
    int numberNodes = data->number_of_nodes();
    weightedAdjacencyMatrix = (double**)malloc(numberNodes * sizeof(double*));
    for (int i = 0; i < numberNodes; i++) 
    {
        weightedAdjacencyMatrix[i] = (double*)malloc(numberNodes * sizeof(double));
    }

    if(knn > 0)
    {
        fill_knn_matrix(weightedAdjacencyMatrix);
    }
    else
    {
        fill_distance_matrix(weightedAdjacencyMatrix);
    }

}

void GraphHandler::print_adj_list()
{
    igraph_vector_int_t el;
    igraph_vector_int_init(&el, 0);
    igraph_get_edgelist(&igraphData.g, &el, 0);
    igraph_integer_t i, j, n;
    n = igraph_ecount(&igraphData.g);
    for (i = 0, j = 0; i < n; i++, j += 2) {
        printf("%" IGRAPH_PRId " --> %" IGRAPH_PRId ": %g\n", VECTOR(el)[j], VECTOR(el)[j + 1], VECTOR(igraphData.weights)[i]);
    }
}

void GraphHandler::create_graph()
{

    double** m = weightedAdjacencyMatrix;
    igraph_integer_t i, j;
    igraph_vector_int_t el;
    igraph_matrix_t mat;
    igraph_vector_int_init(&el, 0);
    igraph_vector_init(&igraphData.weights, 0);
    int numberNodes = data->number_of_nodes();
    igraph_matrix_init(&mat, numberNodes, numberNodes);

    for (i = 0; i < numberNodes; i++) 
    {
        for (j = 0; j < numberNodes; j++) 
        {
            MATRIX(mat, i, j) = m[i][j];
        }
    }

    //filling weights and graph
    igraph_weighted_adjacency(&igraphData.g, &mat, IGRAPH_ADJ_UNDIRECTED, &igraphData.weights, IGRAPH_LOOPS_ONCE);
    igraph_matrix_destroy(&mat);
    igraph_vector_int_destroy(&el);
}

void GraphHandler::create_clustering(double resolution)
{

    igraph_vector_int_t membership;
    igraph_integer_t nb_clusters;
    igraph_real_t quality;

    /* Set default seed to get reproducible results */
    igraph_rng_seed(igraph_rng_default(), 0);
    /*initialize cluster membership vector, which will be filled*/
    //igraph_vector_int_init(&membership, igraph_vcount(&igraphData.g));
    //igraph_vector_int_init(&membership, igraph_vcount(&igraphData.g));
    igraph_vector_int_init_range(&membership, 0, 10);

    if(resolution == 0) 
    {
        resolution = 1.0 / (2 * igraph_ecount(&igraphData.g));
    }
    //number iterations: -1 to run until it does no longer change
    int iterations = 10;
    igraph_community_leiden(&igraphData.g, &igraphData.weights, NULL, resolution, 0.01, 0, iterations, &membership, &nb_clusters, &quality);

    std::cout << "CLUSTER: " << std::to_string(nb_clusters) << " , QUALITY: " << std::to_string(quality) <<  "\n";

    igraph_vector_int_destroy(&membership);
}

void show_results(igraph_t *g, igraph_vector_int_t *membership, igraph_matrix_int_t *memberships, igraph_vector_t *modularity, FILE* f) {
    igraph_integer_t i, j, no_of_nodes = igraph_vcount(g);

    j = igraph_vector_which_max(modularity);
    for (i = 0; i < igraph_vector_int_size(membership); i++) {
        if (VECTOR(*membership)[i] != MATRIX(*memberships, j, i)) {
            fprintf(f, "WARNING: best membership vector element %" IGRAPH_PRId " does not match the best one in the membership matrix\n", i);
        }
    }

    fprintf(f, "Modularities:\n");
    igraph_vector_print(modularity);

    for (i = 0; i < igraph_matrix_int_nrow(memberships); i++) {
        for (j = 0; j < no_of_nodes; j++) {
            fprintf(f, "%" IGRAPH_PRId " ", MATRIX(*memberships, i, j));
        }
        fprintf(f, "||||\n");
    }

    fprintf(f, "\n");
}

// a method to create hierarchical clustering based on edge betweeness: quadric on numbr of edges!!!!
//number of nodes will be the first cluster ID of a cell-merge (it starts from zero with singletons)
void GraphHandler::create_edge_betweenness_clustering(std::vector<std::pair<int, int>>& mergesVector)
{
    igraph_vector_int_t edges;
    igraph_vector_t eb;
    igraph_matrix_int_t merges;
    igraph_vector_int_t membership;

    igraph_vector_init(&eb, 0);
    igraph_vector_int_init(&edges, 0);
    igraph_matrix_int_init(&merges, 0, 0);
    igraph_vector_int_init(&membership, 0);

    igraph_community_edge_betweenness(&igraphData.g, &edges, &eb, &merges,
                                      nullptr, nullptr,
                                      &membership,
                                      IGRAPH_UNDIRECTED,
                                      &igraphData.weights);

    for(long i = 0; i < merges.nrow; ++i)
    {
        mergesVector.emplace_back(std::pair<int, int>(MATRIX(merges, i, 0), MATRIX(merges, i, 1)));
    }

    igraph_matrix_int_print(&merges);

    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&eb);
    igraph_matrix_int_destroy(&merges);
    igraph_vector_int_destroy(&membership);
}

void GraphHandler::create_modularity_clustering()
{
    igraph_vector_t modularity;
    igraph_vector_int_t membership;
    igraph_matrix_int_t memberships;

    igraph_vector_init(&modularity, 0);
    igraph_vector_int_init(&membership, 0);
    igraph_matrix_int_init(&memberships, 0, 0);
    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_community_multilevel(&igraphData.g, 0, 1, &membership, &memberships, &modularity);
    show_results(&igraphData.g, &membership, &memberships, &modularity, stdout);

    igraph_vector_destroy(&modularity);
    igraph_vector_int_destroy(&membership);
    igraph_matrix_int_destroy(&memberships);
}

//find all sets of correlated features
// this is a subgraph of connected nodes (features) in the graph
void GraphHandler::find_correlation_sets(std::vector<std::vector<int>>& correlationSet, const size_t minSetSize)
{
    igraph_vector_int_t corrSet;
    igraph_vector_int_init(&corrSet, 0);
    igraph_integer_t no;
    igraph_connected_components(
        &igraphData.g,
        &corrSet,
        NULL,
        &no,
        IGRAPH_WEAK
    );

    std::vector<std::vector<int>> igraphCorrelationSets(no);
    for (int v = 0; v < igraph_vcount(&igraphData.g); ++v)
    {
        int comp = VECTOR(corrSet)[v];
        igraphCorrelationSets[comp].push_back(v);
    }

    for (auto& comp : igraphCorrelationSets)
    {
        if (comp.size() >= minSetSize)
        {
            std::sort(comp.begin(), comp.end());  // same as your clique sorting
            correlationSet.push_back(comp);
        }
    }

    //destroy igraph cliques vector
    igraph_vector_int_destroy(&corrSet);
}


//calculates cliques: SORTS every clique by its IDs in ascending order
void GraphHandler::calc_all_max_clique(std::vector<std::vector<int>>& cliqueVector, const int minCliqueSize, bool mergeSimilarCliques, bool print)
{

    igraph_vector_int_list_t cliques;
    igraph_integer_t no;

    igraph_vector_int_list_init(&cliques, 0);
    igraph_maximal_cliques(&igraphData.g, &cliques, /*min_size=*/ minCliqueSize,
                           /*max_size=*/ 0 /*no limit*/);
    igraph_maximal_cliques_count(&igraphData.g, &no, /*min_size=*/ minCliqueSize,
                                 /*max_size=*/ 0 /*no limit*/);

    //no clique found
    if (no != igraph_vector_int_list_size(&cliques)) 
    {
        throw std::logic_error("Error in Clique calculation!");
        exit(EXIT_FAILURE);
    }

    //print found cliques
    if(print)
    {
        std::cout << "_____\n";
        for (igraph_integer_t i = 0; i < igraph_vector_int_list_size(&cliques); i++) 
        {
            igraph_vector_int_t *v = igraph_vector_int_list_get_ptr(&cliques, i);
            igraph_vector_int_print(v);
        }
    }

    //safe the found cliques -igraph_vector- in a normal data type(vector of a vector of all proteins) (we ll afterwards delete the igraph object)
    igraph_integer_t igraphListSize = igraph_vector_int_list_size(&cliques);
    for (int i = 0; i < igraphListSize; ++i) 
    {
        igraph_vector_int_t* currentClique = igraph_vector_int_list_get_ptr(&cliques, i);
        // Convert the igraph_vector_int_t to std::vector<int>
        std::vector<int> currentCliqueVector(VECTOR(*currentClique), VECTOR(*currentClique) + igraph_vector_int_size(currentClique));
        //sort this vector for increasing number of proteinIDs
        std::sort(currentCliqueVector.begin(), currentCliqueVector.end());
        // Add the std::vector<int> to the cppVector
        cliqueVector.push_back(currentCliqueVector);
    }

    if(mergeSimilarCliques)
    {
        collapse_all_similar_cliques(cliqueVector, 0.8);
    }
    //destroy igraph cliques vector
    igraph_vector_int_list_destroy(&cliques);
}


void GraphHandler::calc_dense_groups_kcore(std::vector<std::vector<int>>& cliqueVector,
                                           const size_t minCliqueSize,
                                           int kcoreThreshold,
                                           bool mergeSimilarCliques,
                                           bool print)
{
    igraph_vector_int_t coreness;
    igraph_vector_int_init(&coreness, 0);

    // compute coreness
    igraph_coreness(&igraphData.g, &coreness, IGRAPH_ALL);

    igraph_integer_t nodeCount = igraph_vcount(&igraphData.g);

    // collect nodes that satisfy coreness threshold
    std::vector<int> validNodes;
    for (igraph_integer_t i = 0; i < nodeCount; i++)
    {
        if (VECTOR(coreness)[i] >= kcoreThreshold)
        {
            validNodes.push_back(i);
        }
    }

    igraph_vector_int_destroy(&coreness);

    if(validNodes.empty())
        return;

    // build subgraph with those nodes
    igraph_vs_t vs;
    igraph_vs_vector_small(&vs, validNodes.data(), validNodes.size());

    igraph_t subgraph;
    igraph_induced_subgraph(&igraphData.g, &subgraph, vs, IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH);

    igraph_vs_destroy(&vs);

    // find connected components in subgraph
    igraph_vector_int_t membership;
    igraph_vector_int_init(&membership, 0);

    igraph_vector_int_t csize;
    igraph_vector_int_init(&csize, 0);

    igraph_integer_t no;

    igraph_connected_components(&subgraph,
                                &membership,
                                &csize,
                                &no,
                                IGRAPH_WEAK);

    if(print)
        std::cout << "Found " << no << " dense components\n";

    // collect nodes per component
    std::vector<std::vector<int>> components(no);

    for(int i=0;i<igraph_vector_int_size(&membership);i++)
    {
        int comp = VECTOR(membership)[i];
        components[comp].push_back(validNodes[i]);
    }

    // filter by size
    for(auto& comp : components)
    {
        if(comp.size() >= minCliqueSize)
        {
            std::sort(comp.begin(), comp.end());
            cliqueVector.push_back(comp);

            if(print)
            {
                for(int v : comp)
                    std::cout << v << " ";
                std::cout << "\n";
            }
        }
    }

    if(mergeSimilarCliques)
    {
        collapse_all_similar_cliques(cliqueVector, 0.8);
    }

    igraph_vector_int_destroy(&membership);
    igraph_vector_int_destroy(&csize);
    igraph_destroy(&subgraph);
}
