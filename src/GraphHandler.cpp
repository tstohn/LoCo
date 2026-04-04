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
        LOCO_ERR << "When nodes are in descending order we can not apply gaussian kernel, we assumes this is use for protein \
        correclation graphs, where only un-scaled values make sense \n";
        LOCO_EXIT(EXIT_FAILURE);
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
        LOCO_ERR << "For KNN graphs we must order nodes in ascending order \n";
        LOCO_EXIT(EXIT_FAILURE);
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

void GraphHandler::create_graph()
{
    int n = data->number_of_nodes();
    double** m = weightedAdjacencyMatrix;

    graph.num_nodes = n;
    graph.offsets.resize(n + 1);

    // get number of edges for every nodes (exclude self edges)
    std::vector<int> degree(n, 0);

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) 
        {
            double w = m[i][j];

            if (w != 0.0) {
                degree[i]++;
                degree[j]++;
            }
        }
    }

    // set the offset of edges for nodes (we add a last element denoting the end of this element)
    graph.offsets[0] = 0;
    for (int i = 0; i < n; i++) {
        graph.offsets[i + 1] = graph.offsets[i] + degree[i];
    }

    // allocate size for edge/ weight vectors
    int total_edges = graph.offsets[n];
    graph.edges.resize(total_edges);
    graph.weights.resize(total_edges);

    // now store all edges (store them in BOTH direction)
    //cursor stores the current position in edges/ weights where the next node neighbour is stored into
    std::vector<int> cursor = graph.offsets;

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) 
        {
            double w = m[i][j];

            if (w != 0.0) {
                // i -> j
                int pos_i = cursor[i]++;
                graph.edges[pos_i] = j;
                graph.weights[pos_i] = w;

                // j -> i
                int pos_j = cursor[j]++;
                graph.edges[pos_j] = i;
                graph.weights[pos_j] = w;
            }
        }
    }
}

void GraphHandler::dfs_search(const int start,
                                std::vector<char>& visited,
                                std::vector<int>& component,
                                std::vector<int>& stack,
                                const int minDegree
                            )
{
    stack.clear();
    stack.push_back(start);
    visited[start] = 1;

    while (!stack.empty())
    {
        int u = stack.back();
        stack.pop_back();

        component.push_back(u);

        // iterate neighbors
        for (int k = graph.offsets[u]; k < graph.offsets[u + 1]; k++)
        {
            int v = graph.edges[k];

            if(minDegree > 0)
            {
                int deg_v = graph.offsets[v + 1] - graph.offsets[v];
                if(deg_v < minDegree){continue;}
            }

            if (!visited[v])
            {
                visited[v] = 1;
                stack.push_back(v);
            }
        }
    }
}

//find all sets of correlated features
// this is a subgraph of connected nodes (features) in the graph
void GraphHandler::calculate_connected_components(std::vector<std::vector<int>>& correlationSet, 
                                                    const size_t minSetSize)
{
    int n = graph.num_nodes;

    std::vector<char> visited(n, 0);
    std::vector<int> stack;
    stack.reserve(n);

    correlationSet.clear();
    correlationSet.reserve(n / 4); // heuristic

    for (int start = 0; start < n; start++)
    {
        if (visited[start]) continue;

        std::vector<int> component;
        component.reserve(32);

        // call reusable DFS
        dfs_search(start, visited, component, stack, 0);

        if (component.size() >= minSetSize)
        {
            std::sort(component.begin(), component.end());
            correlationSet.push_back(std::move(component));
        }
    }
}

void GraphHandler::bron_kerbosch(std::vector<std::vector<int>>& correlationSet,
                                 std::vector<int>& R,
                                 std::vector<int>& P,
                                 std::vector<int>& X,
                                 const int minSetSize)
{
    // Base case: maximal clique found
    if (P.empty() && X.empty())
    {
        if ((int)R.size() >= minSetSize)
        {
            std::vector<int> clique = R;
            std::sort(clique.begin(), clique.end());
            correlationSet.push_back(std::move(clique));
        }
        return;
    }

    // ---- Pivot Selection ----
    // To maximize pruning, we usually pick the pivot from (P U X) 
    // with the most neighbors in P. Let's stick to your simple choice for now.
    int pivot = -1;
    if (!P.empty()) pivot = P.back();
    else if (!X.empty()) pivot = X.back();

    // ---- Mark neighbors of pivot ----
    std::vector<char> isPivotNeighbor(graph.num_nodes, 0);
    if (pivot != -1)
    {
        for (int k = graph.offsets[pivot]; k < graph.offsets[pivot + 1]; k++)
        {
            isPivotNeighbor[graph.edges[k]] = 1;
        }
    }

    // ---- Pre-map P and X for O(1) intersections ----
    // This replaces the incredibly slow std::find
    std::vector<char> inP(graph.num_nodes, 0);
    for (int v : P) inP[v] = 1;

    std::vector<char> inX(graph.num_nodes, 0);
    for (int v : X) inX[v] = 1;

    // ---- Candidates are P \ N(pivot) ----
    std::vector<int> candidates;
    for (int v : P)
    {
        if (pivot == -1 || !isPivotNeighbor[v])
            candidates.push_back(v);
    }

    // Process candidates
    for (int v : candidates)
    {
        // P might have shrank in previous iterations of this loop!
        // We must make sure v is still in P before processing.
        if (!inP[v]) continue;

        R.push_back(v);

        std::vector<int> newP;
        std::vector<int> newX;

        // Intersect neighbors of v with P and X in O(deg(v)) instead of O(deg(v) * N)
        for (int k = graph.offsets[v]; k < graph.offsets[v + 1]; k++)
        {
            int u = graph.edges[k];
            if (inP[u]) newP.push_back(u);
            if (inX[u]) newX.push_back(u);
        }

        bron_kerbosch(correlationSet, R, newP, newX, minSetSize);

        R.pop_back();

        // Move v from P to X
        inP[v] = 0;
        inX[v] = 1;
        
        // Remove v from P physically
        P.erase(std::remove(P.begin(), P.end(), v), P.end());
        X.push_back(v);
    }
}

//calculates cliques: SORTS every clique by its IDs in ascending order
void GraphHandler::calculate_fully_connected_components(std::vector<std::vector<int>>& correlationSet, const int minSetSize)
{

    int n = graph.num_nodes;

    correlationSet.clear();

    std::vector<int> R;
    std::vector<int> P(n);
    std::vector<int> X;

    // P = all nodes
    for (int i = 0; i < n; i++)
        P[i] = i;

    bron_kerbosch(correlationSet, R, P, X, minSetSize);
}


void GraphHandler::calculate_min_edge_connected_components(std::vector<std::vector<int>>& correlationSet,
                                           const size_t minSetSize,
                                           int minDegree)
{
    int n = graph.num_nodes;

    std::vector<char> visited(n, 0);
    std::vector<int> stack;
    stack.reserve(n);

    correlationSet.clear();

    for (int start = 0; start < n; start++)
    {
        int deg_start = graph.offsets[start + 1] - graph.offsets[start];

        // skip invalid starting nodes
        if (visited[start] || deg_start < minDegree) continue;

        std::vector<int> component;
        component.reserve(32);

        dfs_search(start, visited, component, stack, minDegree);

        if (component.size() >= minSetSize)
        {
            std::sort(component.begin(), component.end());
            correlationSet.push_back(std::move(component));
        }
    }
}
