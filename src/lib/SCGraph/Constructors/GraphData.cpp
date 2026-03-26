#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "Constructors/GraphData.hpp"
// #include "opencv2/core.hpp" //for pca

using namespace nanoflann;

namespace InitializationHelper
{

    void calculate_distance(GraphData* graphData, 
                            const std::vector<int>& cellStateGenes,
                            const std::string& distance,
                            std::vector<std::pair<int, int>> nodeIdxPairs,
                            std::atomic<int>& count, 
                            double totalDist,
                            std::mutex& statusUpdateLock, bool statusUpdate)  
    {
        //fill this tuple with distances between node pairs
        std::vector<std::tuple<int, int, double>> distanceList;

        //iterate through pairs to calculate distance
        for(std::pair<int, int> nodePair : nodeIdxPairs)
        {
            int i = nodePair.first;
            int j = nodePair.second;
            double dist = 0;
            if(distance == "manhattan")
            {
                //if there r no specific state genes we can skip searching for relevant gene idxs
                if(cellStateGenes.empty()){dist = graphData->get_node_at(i)->distance_to(graphData->get_node_at(j), "manhattan");}
                else{dist = graphData->get_node_at(i)->distance_to(graphData->get_node_at(j), cellStateGenes, "manhattan");}
            }
            else if(distance == "euclidean")
            {
                //if there r no specific state genes we can skip searching for relevant gene idxs
                if(cellStateGenes.empty()){dist = graphData->get_node_at(i)->distance_to(graphData->get_node_at(j), "manhattan");}
                else{dist = graphData->get_node_at(i)->distance_to(graphData->get_node_at(j), cellStateGenes, "manhattan");}
            }
            else
            {
                std::cout << "Error: No distance metric for cell-similarity graph\n";
                exit(1);
            }

            //save the distance
            distanceList.push_back(std::make_tuple(i, j, dist));
            count++;
        }
        
        //Process the batch of distances
        //insert this distance into adjacency matrix of both nodes
        statusUpdateLock.lock();
        for (auto& values : distanceList) 
        {
            int i = std::get<0>(values);
            int j = std::get<1>(values);
            double dist = std::get<2>(values);
            graphData->add_distance(i, j, dist);
        }
        statusUpdateLock.unlock();

        if(statusUpdate)
        {
            statusUpdateLock.lock();
            double perc = count / totalDist;
            printProgress(perc);
            statusUpdateLock.unlock();
        }
    }  

//NEW
    void initialize_adjacency_list_knn(GraphData* graphData,
                               const std::vector<int>& cellStateGenes,
                               const std::string& distance,
                               int threads,
                               bool statusUpdate)
    {

        //precalcualte all distances between cells and sort them
        if(graphData->distances_precalcualted())
        {
            if (statusUpdate)
                std::cout << "\tSTEP[0a]:\tCalculate kNN cell-cell distances\n";

            graphData->build_knn_adjacency(graphData->get_knn());

            if (statusUpdate)
                std::cout << "\tSTEP[0b]:\tSorting cell-cell distances\n";

            graphData->sort_adjacencyList();
        }
        else
        {
            if (statusUpdate)
                std::cout << "\tSTEP[0a]:\tCreate KD-tree only\n";

            graphData->create_kd_tree();

            if (statusUpdate)
                std::cout << "Skipping distance calculations and do not fill adj. list\n";

        }



    }

    void initialize_adjacency_list(GraphData* graphData, 
                                    const std::vector<int>& cellStateGenes,
                                    const std::string& distance,
                                    int threads, bool statusUpdate)
    {

        if(statusUpdate)
        {
            std::cout << "\tSTEP[0a]:\tCalculate cell-cell distances\n";
            if(!graphData->distances_precalcualted())
            {
                std::cout << "WARNING: Distances get precalcualted despite seeting the flag to not precalcualte!!!\n";
                std::cout << "Distances can only be NOT calcualted when calcualting distances with the KD-tree\n";
            }
        }
        
        std::cout << "CALCULATING DISTANCES FROM \n";
        for(auto x : cellStateGenes){std::cout << x << " ";}
        std::cout << "\n";

        boost::asio::thread_pool pool_dist(threads);
        std::atomic<int> count = 0;
        double totalDist = (graphData->number_of_nodes() * graphData->number_of_nodes())/2 - graphData->number_of_nodes();
        std::mutex statusUpdateLock;
        //create AdjacencyMatrix
        std::vector<std::pair<int, int>> nodeIdxPairs; //current list of node-pair that we schedule to calcualte their distances
        int numberPairs = 0;
        for(int i=0; i < (graphData->number_of_nodes()-1); ++i)
        {
            for(int j=i+1; j < graphData->number_of_nodes(); ++j)
            {
                //make buckets of 10K pairs
                ++numberPairs;
                nodeIdxPairs.push_back(std::make_pair(i,j));

                //start a thread for this bucket, and empty the bucket
                //COPY THE PAIR BUCKET!!!
                if(numberPairs == 10000)
                {
                    boost::asio::post(pool_dist, std::bind(&calculate_distance, 
                            graphData, std::cref(cellStateGenes), std::cref(distance), nodeIdxPairs, std::ref(count), totalDist, std::ref(statusUpdateLock), statusUpdate));
                    numberPairs = 0;
                    nodeIdxPairs.clear();
                }
            }
        }
        if(!nodeIdxPairs.empty())
        {
            boost::asio::post(pool_dist, std::bind(&calculate_distance, 
                graphData, std::cref(cellStateGenes), std::cref(distance), nodeIdxPairs, std::ref(count), totalDist, std::ref(statusUpdateLock), statusUpdate));
             if(statusUpdate)
            {
                statusUpdateLock.lock();
                double perc = 1;
                printProgress(perc);
                statusUpdateLock.unlock();
            }
        }
        
        pool_dist.join();
        printProgress(1);
        std::cout << "\n";

        if(statusUpdate)
        {
            std::cout << "\tSTEP[0b]:\tSorting cell-cell distances\n";
        }
        //sort the values: before we simply inserted them all to the back
        graphData->sort_adjacencyList();
    }
    //CREATE DATA FROM RAW DATA OF TYPE SINGLECELLDATA
    void cell_similarity_graph(GraphData* graphData, const SingleCellData& inputData, 
                                const std::vector<int>& cellStateGenes, 
                                const std::string& distance,
                                int threads, bool statusUpdate)
    {
        //initialize nodes
        for(int i = 0; i < inputData.cellIDs.size(); ++i)
        {
            std::string nodeName = inputData.cellIDs.at(i);

            const std::vector<double>& point = inputData.pointCloud.at(i);
             //select only the state markers
            if(!cellStateGenes.empty())
            {
                std::vector<double> statePoint;
                statePoint.reserve(cellStateGenes.size());
                for (int idx : cellStateGenes)
                {
                    statePoint.push_back(point[idx]);
                }
                graphData->set_node(statePoint, nodeName);
            }
            else
            {
                graphData->set_node(point, nodeName);
            }
        }
        //create AdjacencyMatrix
        if(graphData->get_knn() > 0)
        {
            initialize_adjacency_list_knn(graphData, cellStateGenes, distance, threads, statusUpdate);
        }
        else
        {
            //TODO!!!
            initialize_adjacency_list(graphData, cellStateGenes, distance, threads, statusUpdate);
        }

        //safe feature names
        graphData->set_feature_names(inputData.geneNames);
    }

    //CREATE GRAPH FROM NODEPTR VECTOR WHEN A GRAPH HAS PREVIOUSLY BEEN CONSTRUCTED
    void cell_similarity_graph(GraphData* graphData, std::vector<nodePtr> inputNodes, 
                                const std::vector<int>& cellStateGenes,
                                const std::string& distance,
                                int threads, bool statusUpdate)
    {
        //initialize nodes
        for(const nodePtr& tmpNode : inputNodes)
        {
            graphData->add_node(tmpNode);
        }
        //create AdjacencyMatrix
        if(graphData->get_knn() > 0)
        {
            initialize_adjacency_list_knn(graphData, cellStateGenes, distance, threads, statusUpdate);
        }
        else
        {
            initialize_adjacency_list(graphData, cellStateGenes, distance, threads, statusUpdate);
        }
    }
}

// cell-cell similarity graph (maybe include euclidean/ manhattan/ cosine later)
void GraphIni::cell_similarity_graph_manhattan_raw(GraphData* graphData, const SingleCellData& inputData, const std::vector<int>& cellStateGenes,int threads, bool statusUpdate)
{
    InitializationHelper::cell_similarity_graph(graphData, inputData, cellStateGenes, "manhattan", threads,  statusUpdate);
}
void GraphIni::cell_similarity_graph_euclidean(GraphData* graphData, const SingleCellData& inputData, const std::vector<int>& cellStateGenes,int threads, bool statusUpdate)
{
    InitializationHelper::cell_similarity_graph(graphData, inputData, cellStateGenes, "euclidean", threads,  statusUpdate);
}
//Initializer funcitons with already created nodePtr instead of raw SingleCellData (e.g., when nodes were already created and subset for new graph)
void GraphIni::cell_similarity_graph_manhattan_nodes(GraphData* graphData, std::vector< nodePtr> inputNodes, const std::vector<int>& cellStateGenes,int threads, bool statusUpdate)
{
    InitializationHelper::cell_similarity_graph(graphData, inputNodes, cellStateGenes, "manhattan", threads,  statusUpdate);
}
void GraphIni::cell_similarity_graph_euclidean(GraphData* graphData, std::vector< nodePtr> inputNodes, const std::vector<int>& cellStateGenes,int threads, bool statusUpdate)
{
    InitializationHelper::cell_similarity_graph(graphData, inputNodes, cellStateGenes, "euclidean", threads,  statusUpdate);
}

// protein correlation graph: BE AWARE: EIGHBORS ARE ORDER HERE IN DESCENDING ORDER (big edges first)
//when we look in this graph for, e.g., cliques of connected nodes we want to have ALL THE HIGHLY correlated ones
//CAREFUL: inserted correlations are NOT ABSOLUTES: since we anywayswant only the HIGHLY correlated ones (not negatively correlated)
// UPDATE: We change correlations to ABSOLUTES: cliqyes can be of many pos and a neg vcorrealated gene, USER has to double check
//what this means!!!!
void GraphIni::protein_correlation_graph(GraphData* graphData, const SingleCellData& inputData, const std::vector<int>& corrStateGenes, int threads, bool statusUpdate)
{
    //set that the order of nodes is descending (default false for, e.g., cell-cell similarity graphs build with cell_similarity functions)
    graphData->set_neighbors_descending();

    //TRANSPOSE single-cell feature counts into protein single-cell counts
    std::vector<std::vector<double>> cellVector = inputData.pointCloud;
    int pointNumber = cellVector.size();

    int geneSize = inputData.geneNames.size();

    //was before: now we selesct the features in setNode
if(!corrStateGenes.empty()){geneSize = corrStateGenes.size();}
    std::vector<std::vector<double>> proteinVector(geneSize, std::vector<double>(pointNumber));

    //without subsetting we just keep all proteins and simply transpose the matrix
    if(corrStateGenes.empty())
    {
        for(int cellID=0; cellID < pointNumber; ++cellID)
        {
            for(int protID=0; protID < cellVector.at(cellID).size(); ++protID)
            {
                proteinVector.at(protID).at(cellID) = cellVector.at(cellID).at(protID);
            }
        }
        }
    else //otherwise we only keep certain genes
    {
        //the old protein ID and the ID in the new vector are no longer same, since we do not keep all proteins
        for(int cellID=0; cellID < pointNumber; ++cellID)
        {
            int newProteinId = 0;
            for(int protID=0; protID < cellVector.at(cellID).size(); ++protID)
            {
                if(std::find(corrStateGenes.begin(), corrStateGenes.end(), protID) == corrStateGenes.end() )
                {
                    continue;
                }
                proteinVector.at(newProteinId).at(cellID) = cellVector.at(cellID).at(protID);
                ++newProteinId;
            }
        }
    }

    //create node for every feature (protein, gene): values are values for every single cell
    int geneID = 0;
    for(const std::vector<double>& scCounts : proteinVector)
    {
        //for empty corrStateGenes we just iterate through the genes
        if(corrStateGenes.empty())
        {
            graphData->set_node(scCounts, inputData.geneNames.at(geneID));
        }
        else
        {
            //otherwise, we need to iterate through one corrStateGene after the other
            graphData->set_node(scCounts, inputData.geneNames.at(corrStateGenes.at(geneID)));
        }
        ++geneID;
    }

    //calcualte rank of all feature for spearman rank (do here only once, instead of everytime we calc a corr)
    // proteinVector is no longer in use, keep it to calcualte rankedList of cell-values per protein
    for(std::vector<double>& values : proteinVector)
    {
        rankify(values);
    }

    //calculate correlations between all features (across cells): adjacencyList
    //BE CAREFUL: proteinVector MUST STILL BE IN THE SAME ORDER AS NODES IN GRAPHDATA
    for(int i=0; i < (graphData->number_of_nodes()-1); ++i)
    {
        for(int j=i+1; j < graphData->number_of_nodes(); ++j)
        {
            std::vector<double> proteinA = proteinVector.at(i);
            std::vector<double> proteinB = proteinVector.at(j);
            double dist = std::abs(calcualte_correlation_coefficient(proteinA, proteinB));
            
            //insert this distance into adjacency matrix of both nodes
            graphData->add_distance(i, j, dist);
        }
    }

    //sort adjacencyList for quick look-up
    graphData->sort_adjacencyList(true);
}

/*void GraphData::perform_pca(const std::vector<std::vector<double>>& pointCloud, const unsigned int& pcas)
{
    //a pointCloud of the highest EigenVectors
    std::vector<std::vector<double>> newPointCloud;
    if(pcas == 0){return;}


    int sz = static_cast<int>(pointCloud.size());
    int dimensions = pointCloud.at(0).size();
    cv::Mat dataPts = cv::Mat(sz, dimensions, CV_64F);
    for(int i = 0; i < dataPts.rows; ++i)
    {
        for( int j = 0; j < dimensions; ++j)
        {
            dataPts.at<double>(i, j) = pointCloud.at(i).at(j);
        }
    }
    cv::PCA pca(dataPts, cv::Mat(), cv::PCA::DATA_AS_ROW, 3);

    for (int i = 0; i < dimensions; i++)
    {
        //eigenvalues ordered from highest to lowest
        double eigenVal = pca.eigenvalues.at<double>(i);
        std::cout << "NEXT EIGENVAL: " << eigenVal << "\n";
    }
}*/


//NEW; FUNCTION FOR DISTANCE INITIALIZATION WITH KNN
void GraphData::build_knn_adjacency(int k)
{
    // 1. Build kd-tree once
    create_kd_tree();

    const size_t N  = nodes.size();
    const int dim   = get_node_at(0)->dimensions();

    // 2. Reset adjacency map
    adjacencyList.clear();
    adjacencyList.reserve(N);  // unordered_map reserve

    std::cout << "creating cell-cell dist with k " << k << "\n";
    // 3. For each node, query its k nearest neighbors (the point itself is part of the results and will be added to cells for this N)
    for (size_t i = 0; i < N; ++i)
    {
        nodePtr node_i = get_node_at(static_cast<int>(i));

        std::vector<double> query_pt(dim);
        for (int d = 0; d < dim; ++d)
            query_pt[d] = node_i->value_at(d);

        std::vector<size_t> ret_indexes(k);
        std::vector<double> out_dists_sqr(k);

        nanoflann::KNNResultSet<double> resultSet(k);
        resultSet.init(ret_indexes.data(), out_dists_sqr.data());

        mat_indexPtr->index->findNeighbors(
            resultSet,
            query_pt.data(),
            nanoflann::SearchParams(10)
        );

        for (int r = 0; r < k; ++r)
        {
            size_t j = ret_indexes[r];

            nodePtr node_j = get_node_at(static_cast<int>(j));
            double dist = std::sqrt(out_dists_sqr[r]);

            // >>> INSERT NEIGHBOR ACCORDING TO HOW OrderedNeighborDistanceHash IS DEFINED <<<
            // Example A: if OrderedNeighborDistanceHash is std::vector<std::pair<nodePtr,double>>
            // neighbors.emplace_back(node_j, dist);

            // Example B: if it's std::map<double, nodePtr>
            // neighbors.emplace(dist, node_j);

            // Example C: if it's std::unordered_map<nodePtr, double>
            // neighbors.emplace(node_j, dist);

            adjacencyList[node_i].set(node_j, dist);
            adjacencyList[node_j].set(node_i, dist);

        }

        // 5. Store in adjacencyList
        //adjacencyList.emplace(node_i, std::move(neighbors));
    }
}

void GraphData::create_kd_tree()
{
    std::cout << "\t# BUILDING KD-TREE\n";
    const int dim = get_node_at(0)->dimensions();
       
    typedef KDTreeVectorOfNodePtrAdaptor<> my_kd_tree_t;
	//my_kd_tree_t   mat_index(dim /*dim*/, nodes, 10 /* max leaf */ );
    mat_indexPtr = std::make_shared<my_kd_tree_t>(dim /*dim*/, nodes, 10 /* max leaf */ );
	mat_indexPtr->index->buildIndex();
}

void GraphData::search_kd_tree()
{
    const int dim = get_node_at(0)->dimensions();

    // do a knn search for every node
    for(int i = 0; i < nodes.size(); ++i)
    {
        std::vector<double> query_pt(dim);
        int dimNum = get_node_at(0)->dimensions();
        for (int d = 0;d < dimNum; d++)
        {
            query_pt[d] = get_node_at(i)->value_at(d);
        }
        const int num_results = 10;
        std::vector<nanoflann::KNNResultSet<double>::IndexType> ret_indexes(num_results);
        std::vector<double> out_dists_sqr(num_results);

        nanoflann::KNNResultSet<double> resultSet(num_results);

        resultSet.init(&ret_indexes[0], &out_dists_sqr[0] );
        mat_indexPtr->index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

        if(i == nodes.size()-1)
        {
        	std::cout << "knnSearch(nn="<<num_results<<"): \n";
            for (int i = 0; i < num_results; i++)
                std::cout << "ret_index["<<i<<"]=" << ret_indexes[i] << " out_dist_sqr=" << out_dists_sqr[i] << std::endl;
        }
    }
}

//mathces is a vector of pairs: 1st index, 2nd is distance
std::vector<std::pair<unsigned long, double> > GraphData::get_points_within_radius(node node, double radius)
{

    std::vector<std::pair<unsigned long, double> > matches;
    nanoflann::SearchParams params;

    std::vector<double> searchPoint = node.all_values();

    mat_indexPtr->index->radiusSearch(&searchPoint[0], radius, matches, params); // returns const size_t variable that stores number of neighbors

    return(matches);
}

//for few comparisons it might be useful to not build a KD-tree and directly calcualte 
//distances
void GraphData::brute_force_get_points_within_radius(node node, double radius)
{
    std::vector<nodePtr> solutionNodes;

    for(nodePtr nodeB : nodes)
    {
        float manhattanDist = 0;
        for(int i = 0; i < nodeB->dimensions(); ++i)
        {
            manhattanDist += std::abs(nodeB->value_at(i) - node.value_at(i));
        }
        if(manhattanDist < radius)
        {
            solutionNodes.push_back(nodeB);
        }
    }

    //std::cout << "MATCHES FOUND: " << solutionNodes.size() << "\n";

}

GraphData::GraphData(const SingleCellData& inputData, const std::vector<int>& cellStateGenes, unsigned int knn, void (*initialize_function)(GraphData*, const SingleCellData&, const std::vector<int>&,int , bool ), int threads, bool statusUpdate, bool precalcAllDistances)
{
    graphKnn = knn;
    distancesPrecalculated = precalcAllDistances;
    stateIndices = cellStateGenes;
    //create nodes from the pointCloud: e.g. single-cell graph, protein correlation graph,...
    initialize_function(this, inputData, cellStateGenes, threads, statusUpdate);
}
GraphData::GraphData( std::vector< nodePtr> inputNodes, const std::vector<int>& cellStateGenes, unsigned int knn, void (*initialize_function)(GraphData*, std::vector< nodePtr>, const std::vector<int>&, int, bool),int threads, bool statusUpdate, bool precalcAllDistances)
{
    graphKnn = knn;
    distancesPrecalculated = precalcAllDistances;
    stateIndices = cellStateGenes;
    //create nodes from the pointCloud: e.g. single-cell graph, protein correlation graph,...
    initialize_function(this, inputNodes, cellStateGenes, threads, statusUpdate);
}