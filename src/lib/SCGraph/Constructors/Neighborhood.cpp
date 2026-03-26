#include <unordered_set>
#include <random>
#include <algorithm> 
#include<limits>
#include <set>

#include "Neighborhood.hpp"

namespace neighborhoodCalculations
{

    std::vector<int> get_random_elements(unsigned int numbers, unsigned int maxNum) 
    {
        if (numbers > maxNum)
            throw std::invalid_argument("numbers > maxNum");

        std::vector<int> pool(maxNum);
        std::iota(pool.begin(), pool.end(), 0);

        std::random_device rd;
        std::mt19937 gen(rd());

        std::shuffle(pool.begin(), pool.end(), gen);

        return std::vector<int>(pool.begin(), pool.begin() + numbers);
    }

    //return LOWEST values first: we want the values with LWOEST laplacian first
    bool sort_corr(const std::pair<int, int>& a, const std::pair<int, int>& b,
                   std::unordered_map< std::pair<int, int>, const double, pair_hash> correlationLaplacian)
    {
        return(correlationLaplacian.at(a) < correlationLaplacian.at(b));
    }

    std::vector<std::vector<int>> generate_subsets(
        const std::vector<int>& set,
        int minSize)
    {
        std::vector<std::vector<int>> subsets;

        int n = set.size();
        int maxMask = 1 << n;

        for(int mask = 0; mask < maxMask; ++mask)
        {
            int count = __builtin_popcount(mask);
            if(count < minSize) continue;

            std::vector<int> subset;
            subset.reserve(count);

            for(int i = 0; i < n; ++i)
            {
                if(mask & (1 << i))
                    subset.push_back(set[i]);
            }

            std::sort(subset.begin(), subset.end());
            subsets.push_back(subset);
        }

        return subsets;
    }


    bool contains_subset(const std::vector<int>& superset,
                        const std::vector<int>& subset)
    {
        return std::includes(
            superset.begin(), superset.end(),
            subset.begin(), subset.end());
    }
}

//it filters the x top pairs for CORRELATIONS and SLOPES
std::vector<std::pair<int, int>> Neighborhood::filter_best_pairs(size_t numberGenes)
{

    std::vector<std::pair<int, int>> correlations;
    for(auto& correlation : correlationLaplacian)
    {
        correlations.push_back(correlation.first);
    }
    //if the correlation filter is bigger than actual found correlation pairs reduce correlation filter to prevent segfault
    if(numberGenes > correlations.size() || numberGenes==0){numberGenes = correlations.size();}

    std::sort(correlations.begin(), correlations.end(),
    [&](const std::pair<int,int>& a, const std::pair<int,int>& b)
    {
        return neighborhoodCalculations::sort_corr(a, b, correlationLaplacian);
    });
    std::vector<std::pair<int, int>> filteredCorrelations(correlations.begin(), correlations.begin() + numberGenes);

    std::sort(correlations.begin(), correlations.end(),
    [&](const std::pair<int,int>& a, const std::pair<int,int>& b)
    {
        return neighborhoodCalculations::sort_corr(a, b, slopeLaplacian);
    });
    std::vector<std::pair<int, int>> filteredCorrelations_2(correlations.begin(), correlations.begin() + numberGenes);

    //insert all new correlatiopns from filteredCorrelations_2 into first vector
    for (const auto& correlation : filteredCorrelations_2) {
        // Check if element is not yet in vector b
        if (std::find(filteredCorrelations.begin(), filteredCorrelations.end(), correlation) == filteredCorrelations.end()) 
        {
            // Insert element into vector b
            filteredCorrelations.push_back(correlation);
        }
    }

    return(filteredCorrelations);
}

void Neighborhood::write_shuffled_laplacians(const std::string& outFile)
{
    std::ofstream outputFile;
    std::size_t found = outFile.find_last_of(".");

    // CORRELATION LAPLACIANS
    std::string outputCorr = outFile.substr(0,found) + "_shuffledCorr." + outFile.substr(found+1);
    std::remove(outputCorr.c_str());
    outputFile.open(outputCorr, std::ofstream::app);

    // collect all pairs
    std::vector<std::pair<int,int>> pairs;
    for(const auto& it : shuffledCorrLaplacians)
    {
        pairs.push_back(it.first);
    }

    if(pairs.empty())
    {
        outputFile.close();
        return;
    }

    int permutations = shuffledCorrLaplacians.at(pairs[0]).size();

    // HEADER
    outputFile << "Permutation";
    for(const std::pair<int,int>& pair : pairs)
    {
        int name_a_idx = pair.first;
        int name_b_idx = pair.second;
        if(!corrStateGenes.empty())
        {
            name_a_idx = corrStateGenes.at(name_a_idx);
            name_b_idx = corrStateGenes.at(name_b_idx);
        }
        outputFile << "\t"
                   << inputDataOrigional.geneNames.at(name_a_idx)
                   << "_"
                   << inputDataOrigional.geneNames.at(name_b_idx);
    }
    outputFile << "\n";

    // rows = permutations
    for(int p = 0; p < permutations; ++p)
    {
        outputFile << p;

        for(const std::pair<int,int>& pair : pairs)
        {
            outputFile << "\t" << shuffledCorrLaplacians.at(pair).at(p);
        }

        outputFile << "\n";
    }

    outputFile.close();


    // SLOPE LAPLACIANS
    std::string outputSlope = outFile.substr(0,found) + "_shuffledSlope." + outFile.substr(found+1);
    std::remove(outputSlope.c_str());
    outputFile.open(outputSlope, std::ofstream::app);

    pairs.clear();
    for(const auto& it : shuffledSlopeLaplacians)
    {
        pairs.push_back(it.first);
    }

    if(pairs.empty())
    {
        outputFile.close();
        return;
    }

    permutations = shuffledSlopeLaplacians.at(pairs[0]).size();

    // HEADER
    outputFile << "Permutation";
    for(const std::pair<int,int>& pair : pairs)
    {
        int name_a_idx = pair.first;
        int name_b_idx = pair.second;
        if(!corrStateGenes.empty())
        {
            name_a_idx = corrStateGenes.at(name_a_idx);
            name_b_idx = corrStateGenes.at(name_b_idx);
        }
        outputFile << "\t"
                   << inputDataOrigional.geneNames.at(name_a_idx)
                   << "_"
                   << inputDataOrigional.geneNames.at(name_b_idx);
    }
    outputFile << "\n";

    // rows = permutations
    for(int p = 0; p < permutations; ++p)
    {
        outputFile << p;

        for(const std::pair<int,int>& pair : pairs)
        {
            outputFile << "\t" << shuffledSlopeLaplacians.at(pair).at(p);
        }

        outputFile << "\n";
    }

    outputFile.close();
}

void Neighborhood::write_results_to_file(const std::string& outFile, const std::string& prefix, int& numberCorrelations)
{
    std::ofstream outputFile;   

    //get all the pairs that we want to write out
    std::vector<std::pair<int, int>> filteredPairs = pairs;
    //sort correlation pairs (and filter if the number of correlations to fitler is < total number of correlation pairs)
    filteredPairs = filter_best_pairs(numberCorrelations);

    //LAPLCACIAN SCORES FILE
    //1FILE with annotations: corrpair - LaplacianCorrscore - LaplacianSlopeScore - origionalCliqueItIsFrom
    std::string outputLaplace = outFile + "/" + prefix + "_laplacian.tsv";
    std::remove(outputLaplace.c_str());
    outputFile.open (outputLaplace, std::ofstream::app);
    //write HEADER
    outputFile << "ProteinPair" << "\t" << "CorrelationScore" << "\t" << "p_CorrelationScore" <<
    "\t" << "SlopeScore"  "\t" << "p_SlopeScore" <<
     "\t" << "OrigionalCliques" << "\n";
    //write LINES
    for(const std::pair<int, int>& pair : filteredPairs) 
    {
        //the name indices might have to be adjusted if we considered only a subset of the names
        int name_a_idx = pair.first;
        int name_b_idx = pair.second;
        if(!corrStateGenes.empty())
        {
            name_a_idx = corrStateGenes.at(name_a_idx);
            name_b_idx = corrStateGenes.at(name_b_idx);
        }

        outputFile << inputDataOrigional.geneNames.at(name_a_idx) << "_" << inputDataOrigional.geneNames.at(name_b_idx) << "\t" << 
        correlationLaplacian.at(pair) << "\t" << p_corr.at(pair) << "\t" <<
        slopeLaplacian.at(pair) << "\t" << p_slope.at(pair) << "\t";

        size_t count = 0;
        for(const auto& clique : pairToClique.at(pair))
        {
            for(size_t i = 0; i < clique.size(); ++i)
            {
                outputFile << clique.at(i);
                if(i != (clique.size()-1) )
                {
                    outputFile << ",";
                }
            }
            ++count;
            if(count != pairToClique.at(pair).size())
            {
                outputFile << ";";
            }
        }

        outputFile << "\n";
    }
    outputFile << "\n";
    outputFile.close();

    //2FILES with Matrix:  for NeighborhoodID * CorrpAir: value == corr/ slope
    //CORRELATIONS FILE
    std::string outputCorr = outFile + "/" + prefix + "_correlations.tsv";
    std::remove(outputCorr.c_str());
    outputFile.open (outputCorr, std::ofstream::app);
    //write HEADER
    outputFile << "Neighborhood";
    for(const std::pair<int, int>& pair : filteredPairs) 
    {
        int name_a_idx = pair.first;
        int name_b_idx = pair.second;
        if(!corrStateGenes.empty())
        {
            name_a_idx = corrStateGenes.at(name_a_idx);
            name_b_idx = corrStateGenes.at(name_b_idx);
        }
        outputFile << "\t" << inputDataOrigional.geneNames.at(name_a_idx) << "_" << inputDataOrigional.geneNames.at(name_b_idx);
    }

    outputFile << "\n";
    //write LINES (one line per neighborhood)
    for(const auto& neighborhoodResult : corrResult)
    {
        outputFile << neighborhoodResult.first->get_name();
        for(const std::pair<int, int>& pair : filteredPairs) 
        {
            outputFile << "\t" << neighborhoodResult.second.correlationResult.at(pair);
        }
        outputFile << "\n";
    }

    outputFile.close();

    //SLOPES FILE
    std::string outputSlope = outFile + "/" + prefix + "_slopes.tsv";
    std::remove(outputSlope.c_str());
    outputFile.open (outputSlope, std::ofstream::app);
    //write HEADER
    outputFile << "Neighborhood";
    for(const std::pair<int, int>& pair : filteredPairs) 
    {
        int name_a_idx = pair.first;
        int name_b_idx = pair.second;
        if(!corrStateGenes.empty())
        {
            name_a_idx = corrStateGenes.at(name_a_idx);
            name_b_idx = corrStateGenes.at(name_b_idx);
        }
        outputFile << "\t" << inputDataOrigional.geneNames.at(name_a_idx) << "_" << inputDataOrigional.geneNames.at(name_b_idx);
    }

    outputFile << "\n";
    //write LINES (one line per neighborhood)
    for(const auto& neighborhoodResult : corrResult)
    {
        outputFile << neighborhoodResult.first->get_name();
        for(const std::pair<int, int>& pair : filteredPairs) 
        {
            outputFile << "\t" << neighborhoodResult.second.slopeResult.at(pair);
        }
        outputFile << "\n";
    }
    outputFile.close();

    //1 File with annotations: NeighborhoodID - centralNodeDimensions[start ... end]
    std::string outputNeighborhoodCoordinates = outFile + "/" + prefix + "_coord.tsv";
    std::remove(outputNeighborhoodCoordinates.c_str());
    outputFile.open (outputNeighborhoodCoordinates, std::ofstream::app);
    //write HEADER
    outputFile << "Neighborhood";
    for(const std::string& feature : get_feaure_names()) 
    {
        outputFile << "\t" << feature;
    }
    outputFile << "\n";
    //write lines
    for(const nodePtr& node : neighborhoodGraph->get_all_nodes())
    {
        outputFile << node->get_name();
        for(const double& dim : node->all_values())
        {
            outputFile << "\t" << dim;
        }
        outputFile << "\n";
    }
    outputFile.close();

    //1 File with cellIDs in each neighborhood: NeighborhoodID in header - all cellIDs in rows for each neighborhood columns
    std::string outputNeighborhoodCells = outFile + "/" + prefix + "_cells.tsv";
    std::remove(outputNeighborhoodCells.c_str());
    outputFile.open (outputNeighborhoodCells, std::ofstream::app);
    //write HEADER (all neigborhood names)
    const std::vector<nodePtr> neighborHoods = neighborhoodGraph->get_all_nodes();
    for(size_t neiborhoodID = 0; neiborhoodID < neighborHoods.size(); ++neiborhoodID)
    {
        outputFile << neighborHoods.at(neiborhoodID)->get_name();
        if(neiborhoodID < (neighborHoods.size()-1)){outputFile << "\t";}
    }  

    outputFile << "\n";
    //write lines
    //iterate over cellIDs (0 - <number cells in neighborhood>), for every id write out the cell with this id for every neighborhood
    for(size_t cellID = 0; cellID < neighborhoodSize; ++cellID)
    {
        //for every neighborhood (in columns)
        for(size_t neiborhoodID = 0; neiborhoodID < neighborHoods.size(); ++neiborhoodID)
        {
            outputFile << neighborhoods.at(neighborHoods.at(neiborhoodID)).at(cellID);
            if(neiborhoodID < (neighborHoods.size()-1)){outputFile << "\t";}
        }
        outputFile << "\n";
    }
    outputFile.close();

    //WRITE PERMUTATION LAPLACIANS
    write_shuffled_laplacians(outFile);

}

void Neighborhood::calculate_correlations(unsigned int numberNodes, bool print, 
                                          const std::pair<int, int>& correlationpair,
                                          std::unordered_map<const std::pair<int, int>, double, pair_hash>& corrVariance,
                                          int totalCount, double& currentCount)
{
    unsigned int nonNanNodes = 0;
    //calcualte mean first
    double mean = 0;
    for(const nodePtr& node: neighborhoodGraph->get_all_nodes())
    {
        double correlationTmp = corrResult.at(node).correlationResult.at(correlationpair);
        if(!isnan(correlationTmp))
        {
            mean += correlationTmp;
            nonNanNodes++;
        }
    }

    double variance = 0;

    //if non-nan nodes are in less than 10% of the neighborhoods (>90% are nan-nodes), 
    // then assign the whole variance adn thereby laplacian to nan
    if(nonNanNodes < numberNodes/10)
    {
       variance = 0; //leave variance explicitely at zero. Laplacian score will therefore also be zero due to division by zero
    }
    else
    {
        mean /= nonNanNodes;
        //calcualte deviations from mean
        for(const nodePtr& node: neighborhoodGraph->get_all_nodes())
        {
            double correlationTmp = corrResult.at(node).correlationResult.at(correlationpair);
            if(!isnan(correlationTmp))
            {
                variance += pow((correlationTmp-mean), 2);
            }
        }

        variance /= nonNanNodes;
    }

    if(print)
    {
        std::cout << "VARIANCE: " << correlationpair.first << "-" << correlationpair.second << ": " << variance << "\n";
    }
    threadLock.lock();
    corrVariance.insert(std::pair<const std::pair<int, int>, double>(correlationpair, variance));
    ++currentCount; 
    double perc = currentCount / totalCount;
    printProgress(perc);
    threadLock.unlock();
}

void Neighborhood::calculate_slopes(unsigned int numberNodes,
                                    const std::pair<int, int>& correlationpair,
                                    std::unordered_map<const std::pair<int, int>, double, pair_hash>& slopeVariance,
                                    int totalCount, double& currentCount)
{
        //calcualte mean first
    double mean = 0;
    uint nonNanNodes = 0;
    for(const nodePtr& node: neighborhoodGraph->get_all_nodes())
    {
        double slopeTmp = corrResult.at(node).slopeResult.at(correlationpair);
        if(!isnan(slopeTmp))
        {
            mean += slopeTmp;
            nonNanNodes++;
        }
    }

    double variance = 0;
    if(nonNanNodes < numberNodes/10)
    {
       variance = 0; //enforce nan values further down by 0-division
    }
    else
    {
        mean /= nonNanNodes;
        //calcualte deviations from mean
        for(const nodePtr& node: neighborhoodGraph->get_all_nodes())
        {
            double slopeTmp = corrResult.at(node).slopeResult.at(correlationpair);
            if(!isnan(slopeTmp))
            {
                variance += pow((slopeTmp-mean), 2);
            }
        }
        variance /= nonNanNodes;
    }

    threadLock.lock();
    slopeVariance.insert(std::pair<const std::pair<int, int>, double>(correlationpair, variance));
    ++currentCount; 
    double perc = currentCount / totalCount;
    printProgress(perc);
    threadLock.unlock();
}

void Neighborhood::laplacian_score(const std::pair<int, int>& pair,
                                   const std::unordered_map<const std::pair<int, int>, double, pair_hash>& corrVariance,
                                   const std::unordered_map<const std::pair<int, int>, double, pair_hash>& slopeVariance,
                                   int totalCount, double& currentCount)
{

    //iterate over adjacency matrix (its a double** weightMatrix in GraphHandler)
    unsigned int nodeNumber = neighborhoodGraph->get_all_nodes().size();

    double corrWeightSum = 0;
    double slopeWeightSum = 0;

    for(unsigned int i = 0; i < (nodeNumber-1); ++i)
    {
        for(unsigned int j = i+1; j < nodeNumber; ++j)
        {
            nodePtr nodeA = neighborhoodGraph->get_node_at(i);
            nodePtr nodeB = neighborhoodGraph->get_node_at(j);

            //edge weight for closeness of values
            double weight = neighborhoodGraph->get_edge_weight_between_nodes(nodeA, nodeB);
            if(weight == 0){continue;}

            //actual feature values for nodes
            double featureNodeACorr = corrResult.at(nodeA).correlationResult.at(pair);
            double featureNodeBCorr = corrResult.at(nodeB).correlationResult.at(pair);
            double featureNodeASlope = corrResult.at(nodeA).slopeResult.at(pair);
            double featureNodeBSlope = corrResult.at(nodeB).slopeResult.at(pair);

            if( !isnan(featureNodeACorr) && !isnan(featureNodeBCorr))
            {
                corrWeightSum += weight * pow((featureNodeACorr - featureNodeBCorr), 2); //square the result
            }
            if( !isnan(featureNodeASlope) && !isnan(featureNodeBSlope))
            {
                slopeWeightSum += weight * pow((featureNodeASlope - featureNodeBSlope), 2); //square the result
            }
        }
    }

    //substract variance of the correlations
    corrWeightSum /= corrVariance.at(pair);
    slopeWeightSum /= slopeVariance.at(pair);

    //store the result, for every pair of nodes (here used pair in indices) - indices are the feature indices as in origional data/ and all node feature of Neighborhood Graph
    threadLock.lock();
    correlationLaplacian.insert(std::pair<const std::pair<int, int>, const double>(pair, corrWeightSum));
    slopeLaplacian.insert(std::pair<const std::pair<int, int>, const double>(pair, slopeWeightSum));
    ++currentCount; 
    double perc = currentCount / totalCount;
    printProgress(perc);
    threadLock.unlock();
}

//similar to laplacian score but with randomly shuffled correlations
//we shuffle the whole group of corelations for all nieghborhoods, that means calculated correlations per neighborhood stay together
//but get assigned to a new neighborhood
void Neighborhood::laplacian_significance(
    const std::pair<int,int>& pair,
    const std::unordered_map<const std::pair<int,int>,double,pair_hash>& corrVariance,
    const std::unordered_map<const std::pair<int,int>,double,pair_hash>& slopeVariance,
    const std::vector<CorrelationPropagationResult>& vectorizedResults,
    int totalCount,
    double& currentCount)
{
    thread_local std::mt19937 rng(std::random_device{}());
    const unsigned int nodeNumber = neighborhoodGraph->get_all_nodes().size();

    // ----- cache variances (1 map lookup only) -----
    const double corrVar  = corrVariance.at(pair);
    const double slopeVar = slopeVariance.at(pair);

    // ----- cache observed laplacians -----
    const double corrObserved  = correlationLaplacian.at(pair);
    const double slopeObserved = slopeLaplacian.at(pair);

    // ----- cache feature values (avoid map lookup in loops) -----
    std::vector<double> corrValues(nodeNumber);
    std::vector<double> slopeValues(nodeNumber);

    for(unsigned int i=0;i<nodeNumber;i++)
    {
        corrValues[i]  = vectorizedResults[i].correlationResult.at(pair);
        slopeValues[i] = vectorizedResults[i].slopeResult.at(pair);
    }

    // ----- build edge list once -----
    struct Edge { int i,j; double w; };
    std::vector<Edge> edges;
    edges.reserve(nodeNumber * 5); // guess due to KNN graph, will be more in reality since its not a directed graph

    for(unsigned int i=0;i<nodeNumber-1;i++)
    {
        nodePtr nodeA = neighborhoodGraph->get_node_at(i);

        for(unsigned int j=i+1;j<nodeNumber;j++)
        {
            nodePtr nodeB = neighborhoodGraph->get_node_at(j);

            double w = neighborhoodGraph->get_edge_weight_between_nodes(nodeA,nodeB);

            if(w!=0.0)
                edges.push_back({(int)i,(int)j,w});
        }
    }

    // ----- permutation vector -----
    std::vector<int> perm(nodeNumber);
    std::iota(perm.begin(),perm.end(),0);

    // ----- result storage -----
    std::vector<double> shuffledCorrLaplacianVector(permutations);
    std::vector<double> shuffledSlopeLaplacianVector(permutations);

    // ----- permutation loop -----
    for(int p=0;p<permutations;p++)
    {
        std::shuffle(perm.begin(),perm.end(),rng);

        double corrSum  = 0.0;
        double slopeSum = 0.0;

        for(const Edge& e : edges)
        {
            int pi = perm[e.i];
            int pj = perm[e.j];

            double corrA = corrValues[pi];
            double corrB = corrValues[pj];

            double slopeA = slopeValues[pi];
            double slopeB = slopeValues[pj];

            if(!std::isnan(corrA) && !std::isnan(corrB))
            {
                double d = corrA - corrB;
                corrSum += e.w * d * d;
            }

            if(!std::isnan(slopeA) && !std::isnan(slopeB))
            {
                double d = slopeA - slopeB;
                slopeSum += e.w * d * d;
            }
        }

        shuffledCorrLaplacianVector[p]  = corrSum  / corrVar;
        shuffledSlopeLaplacianVector[p] = slopeSum / slopeVar;
    }

    // ----- p-value calculation -----
    double p_corr_tmp  = 0.0;
    double p_slope_tmp = 0.0;

    for(int i=0;i<permutations;i++)
    {
        if(shuffledCorrLaplacianVector[i] <= corrObserved)
            p_corr_tmp++;

        if(shuffledSlopeLaplacianVector[i] <= slopeObserved)
            p_slope_tmp++;
    }

    p_corr_tmp  /= permutations;
    p_slope_tmp /= permutations;

    // ----- store results -----
    threadLock.lock();

    p_corr.insert({pair,p_corr_tmp});
    p_slope.insert({pair,p_slope_tmp});

    shuffledCorrLaplacians.insert({pair,shuffledCorrLaplacianVector});
    shuffledSlopeLaplacians.insert({pair,shuffledSlopeLaplacianVector});

    ++currentCount;
    double perc = currentCount / totalCount;
    printProgress(perc);

    threadLock.unlock();
}

//for now simply calculate sum of all edge weights
// -> smaller values mean tiny changes between nodes
//in origional paper is a matrix-multiplication formular
void Neighborhood::calculate_laplacian_score(bool print, int threads)
{   

    //calculate once the variance that we observe in all correlations/ slopes
    unsigned int numberNodes = neighborhoodGraph->get_all_nodes().size();

    //CORRELATION VARIANCE
    std::cout << "\tSTEP[6a]:\tCalculate variance of smooth Correlations\n";
    std::unordered_map<const std::pair<int, int>, double, pair_hash> corrVariance;
    boost::asio::thread_pool pool_corr(threads);
    double count = 0;
    for(const std::pair<int, int>& correlationpair : pairs)
    {
        boost::asio::post(pool_corr, std::bind(&Neighborhood::calculate_correlations, this, 
                          numberNodes, print, std::cref(correlationpair), std::ref(corrVariance), pairs.size(), std::ref(count)));
    }
    pool_corr.join();
    printProgress(1);
    std::cout << "\n";

    //SLOPE VARIANCE
    std::cout << "\tSTEP[6b]:\tCalculate variance of smooth Slopes\n";
    std::unordered_map<const std::pair<int, int>, double, pair_hash> slopeVariance;
    boost::asio::thread_pool pool_slope(threads);
    count = 0;
    for(const std::pair<int, int>& correlationpair : pairs)
    {
        boost::asio::post(pool_slope, std::bind(&Neighborhood::calculate_slopes, this, 
                          numberNodes, std::cref(correlationpair), std::ref(slopeVariance), pairs.size(), std::ref(count)));
    }
    pool_slope.join();
    printProgress(1);
    std::cout << "\n";

    //FULL LAPLACIAN SCORE
    std::cout << "\tSTEP[6c]:\tCalculate Laplacian Scores\n";
    boost::asio::thread_pool pool_lapl(threads);
    count = 0;
    for(const std::pair<int, int>& pair : pairs)
    {
        boost::asio::post(pool_lapl, std::bind(&Neighborhood::laplacian_score, this, 
                          std::cref(pair), std::cref(corrVariance), std::cref(slopeVariance),
                          pairs.size(), std::ref(count)));
    }
    pool_lapl.join();
    printProgress(1);
    std::cout << "\n";

    //calcualte significance: calcualte laplacian score for randomly shuffled features across N
    std::cout << "\tSTEP[6d]:\tCalculate Significance: Laplacian Scores after shuffling correlations\n";
    std::cout << "\t\t We keep the KNN graph and randomly reassign all correlations/slops of a neighborhood to a different neighborhood\n";
    count = 0;
    //firstly write the correlation/slope results to a vector that is easily shufflable
    std::vector<CorrelationPropagationResult> vectorizedResults;
    vectorizedResults.reserve(corrResult.size());
    for (const auto& [node, result] : corrResult)
    {
        vectorizedResults.push_back(result);
    }
    boost::asio::thread_pool pool_shuffle(threads);
    for(const std::pair<int, int>& pair : pairs)
    {
        boost::asio::post(pool_shuffle, std::bind(&Neighborhood::laplacian_significance, this, 
                          std::cref(pair), std::cref(corrVariance), std::cref(slopeVariance),
                          std::cref(vectorizedResults),
                          pairs.size(), std::ref(count)));
    }
    pool_shuffle.join();
    printProgress(1);
    std::cout << "\n";

}

//fills the reuslt per neighborhood:
//we take correlations from already calcualted neighborhood data
//slopes have to be calculated again

//only do it for interesting correlations within cliquesVector
void Neighborhood::extract_pairs_from_correlation_sets(std::unordered_map<nodePtr, std::shared_ptr<GraphData>>& neighborhoodCorrelations)
{

    //iterte through all cliques & insert all pairs in pair vector
    for(const std::vector<int>& cliques : cliquesVector)
    {
        //for all pairs of nodes from the clique
        for(size_t i = 0; i < (cliques.size()-1); ++i)
        {
            for(size_t j = (i+1); j < cliques.size(); ++j)
            {
                std::pair<int, int> tmpPair(cliques.at(i), cliques.at(j));
                //if the pair has a reported clique but needs to be updated
                if(pairToClique.find(tmpPair) != pairToClique.end())
                {
                    std::vector<std::string> features;
                    for(int cIdx : cliques)
                    {
                        std::string featureName = inputDataOrigional.geneNames.at(cIdx);
                        if(!corrStateGenes.empty())
                        {
                            int correctedIndex = corrStateGenes.at(cIdx);
                            featureName =  inputDataOrigional.geneNames.at(correctedIndex);
                        }
                        features.push_back(featureName);
                    }
                    pairToClique.at(tmpPair).push_back(features);
                }
                else 
                {
                    //if the pair has not been reported so far
                    std::vector<std::string> features;
                    for(int cIdx : cliques)
                    {
                        std::string featureName = inputDataOrigional.geneNames.at(cIdx);
                        if(!corrStateGenes.empty())
                        {
                            int correctedIndex = corrStateGenes.at(cIdx);
                            featureName =  inputDataOrigional.geneNames.at(correctedIndex);
                        }
                        features.push_back(featureName);
                    }
                    std::vector<std::vector<std::string>> cliquesVector;
                    cliquesVector.push_back(features);
                    pairToClique.insert(std::pair<const std::pair<int, int>, std::vector<std::vector<std::string>>>(tmpPair, cliquesVector));
                }
                if(std::find(pairs.begin(), pairs.end(), tmpPair) == pairs.end()) //if the pair was not stored yet
                {
                    pairs.push_back(tmpPair);
                }
            }
        }
    }

    //PROCESS ALL THOSE PAIRS
    for(const nodePtr& neighborhoodCenter : centralNeighborhoodPtrs)
    {
        //temporary result for this neighborhood
        CorrelationPropagationResult tmpResult;

        //all feature-pairs (correlations between two features)
        for(const std::pair<int, int>& pair : pairs)
        {
            //get current neighborhood Correlation Data
            std::shared_ptr<GraphData> tmpData = neighborhoodCorrelations.at(neighborhoodCenter);

            //store all CORRELATIONS for those pairs (get them from already calculated neighborhoodCorrelations)
            nodePtr featureNodeA = tmpData->get_node_at(pair.first);
            nodePtr featureNodeB = tmpData->get_node_at(pair.second);
            //re-calcualte correlations (before we had absolute values)
            //const double corr = tmpData->get_distance_between_nodes(featureNodeA, featureNodeB);
            const double corr = calcualte_correlation_coefficient(featureNodeA->all_values(), featureNodeB->all_values());
            tmpResult.correlationResult.insert(std::pair< std::pair<int, int>, double >(pair, corr));

            //store all SLOPES for those pairs (have to calculate new)
            const double slope = calculate_slope(featureNodeA->all_values(), featureNodeB->all_values());
            tmpResult.slopeResult.insert(std::pair< std::pair<int, int>, double >(pair, slope));
        }

        corrResult.insert(std::pair<nodePtr, CorrelationPropagationResult>(neighborhoodCenter, tmpResult));
    }
    
}

void Neighborhood::detect_cliques_in_neighborhood(nodePtr neighborhoodCenter, const double& correlationStrengthCutoff, int minCliqueSize, const bool printFoundCliquesPerNeighborhood,
                                                  std::unordered_map<nodePtr, std::vector<std::vector<int>>>& cliquesPerNeighborhood, 
                                                  std::unordered_map<nodePtr, std::shared_ptr<GraphData>>& neighborhoodCorrelations,
                                                  int totalCount, double& currentCount)
{
        //calculate all correlations/ slopes of them
        SingleCellData inputDataTmp = filter_singleCelldata(inputDataOrigional, neighborhoods.at(neighborhoodCenter));

        //for the protein graph calculate brute force all distances instead of reading form a KD-tree
        // this way we can filter correlations based on value
        unsigned int proteinGraphKnn = 0;
        std::shared_ptr<GraphData> correlationData = std::make_shared<GraphData>(inputDataTmp, corrStateGenes, proteinGraphKnn, &GraphIni::protein_correlation_graph);
        GraphHandler corrGraphBuilder = GraphHandler(correlationData, 0, correlationStrengthCutoff, -1);
        corrGraphBuilder.create_graph();

        //results are stored here afterwards
        std::vector<std::vector<int>> cliqueVectorRaw;
        if(correlatedSetMode == 0)
        {
            corrGraphBuilder.find_correlation_sets(cliqueVectorRaw, minCliqueSize);
        }
        else if(correlatedSetMode == 1)
        {
            corrGraphBuilder.calc_all_max_clique(cliqueVectorRaw, minCliqueSize, false, printFoundCliquesPerNeighborhood); //no clique merging, no printing
        }
        else if(correlatedSetMode == 2)
        {
            int minEdgeNumPerNode = 2;
            corrGraphBuilder.calc_dense_groups_kcore(cliqueVectorRaw, minCliqueSize, minEdgeNumPerNode, false, printFoundCliquesPerNeighborhood);
        }
        else
        {
            std::cout << "Invalid mode for detection of correlated sets: fallback to connected component -> mode 0\n";
            corrGraphBuilder.find_correlation_sets(cliqueVectorRaw, minCliqueSize);
        }

        //insert cliques in map<neighborhoodID -> cliques>
        threadLock.lock();
        cliquesPerNeighborhood.insert(std::make_pair(neighborhoodCenter, cliqueVectorRaw));
        neighborhoodCorrelations.insert(std::pair<nodePtr, std::shared_ptr<GraphData>>(neighborhoodCenter, correlationData));
        ++currentCount; 
        double perc = currentCount / totalCount;
        printProgress(perc);
        threadLock.unlock();
}

//filters cliques that do not 'randomly' occur in one neighborhood, but the ones that are shared also among the neighbors of this neighborhood
//(NOT IUSED AT THIS POINT: Since we do not look at subsets of cliques it almost never happens that the EXACT CLIQUE of several proteins is the same...)
void Neighborhood::filter_cliques_present(nodePtr neighborhoodCenter, std::unordered_map<nodePtr, std::vector<std::vector<int>>>& cliquesPerNeighborhood)
{
        //get closest neighborhoods (now only knn neighbors in graph)
        //TODO: connect via graph of whole data and BFS through it for neighbors (igraph_bfs_callback.c with own callback function)
        //for now just get 'numberNeighbors' closest neighbors for our main node (maybe later use graph and BFS)
        int numberNeighbors = 5;
        std::vector<nodePtr> neighborIDs = neighborhoodGraph->get_data()->get_adjacent_nodes_knn(neighborhoodCenter, numberNeighbors);

        //find common cliques in all neighborhoods
        //for every clique in actual nHood
        for(std::vector<int> clique : cliquesPerNeighborhood.at(neighborhoodCenter))
        {

            bool cliqueInNeighbors = true;
            //check if it is contained in all neighboring neighborhoods
            for(nodePtr neighboringNHood : neighborIDs)
            {
                std::vector<std::vector<int>> nCliques = cliquesPerNeighborhood.at(neighboringNHood);

                //the cliques are all sorted in ascending order
                if(std::find(nCliques.begin(), nCliques.end(), clique) == nCliques.end()) //if the clique is not in the neighboring cliques
                {
                    cliqueInNeighbors = false;
                    break;
                }
            }

            //if clique is not existing in neighbors reject it (filtering of cliques that are only in ONE neighborhood and probably outliers)
            if(!cliqueInNeighbors){continue;}

            threadLock.lock();
            //if we managed to run through for loop, this clique is in all neighboring-neighborHoods
            if(std::find(cliquesVector.begin(), cliquesVector.end(), clique) == cliquesVector.end()) //if the clique does not exist yet
            {
                cliquesVector.push_back(clique);
            }
            threadLock.unlock();
        }
}

void Neighborhood::filter_consistent_correlation_sets(
        const std::unordered_map<nodePtr, std::vector<std::vector<int>>>& cliquesPerNeighborhood,
        int minCliqueSize)
{
    std::set<std::vector<int>> candidateSets;

    // ----- generate all candidate subsets -----
    for(const auto& clique : cliquesVector)
    {
        std::vector<int> sortedClique = clique;
        std::sort(sortedClique.begin(), sortedClique.end());

        auto subsets = neighborhoodCalculations::generate_subsets(sortedClique, minCliqueSize);

        for(const auto& s : subsets)
            candidateSets.insert(s);
    }

    const int neighborhoodCount = cliquesPerNeighborhood.size();
    const int required = std::ceil(minimumCorrSetAbundance * neighborhoodCount);

    std::vector<std::vector<int>> result;

    // ----- check each candidate -----
    for(const auto& candidate : candidateSets)
    {
        int count = 0;

        for(const auto& [node, sets] : cliquesPerNeighborhood)
        {
            bool found = false;

            for(const auto& s : sets)
            {
                std::vector<int> sortedS = s;
                std::sort(sortedS.begin(), sortedS.end());

                if(neighborhoodCalculations::contains_subset(sortedS, candidate))
                {
                    found = true;
                    break;
                }
            }

            if(found)
                count++;

            if(count >= required)
                break;
        }

        if(count >= required)
            result.push_back(candidate);
    }

    cliquesVector = result;
}

void Neighborhood::calculate_correlation_propagation(double correlationStrengthCutoff, const bool printFoundCliquesPerNeighborhood, int minCliqueSize, int thread)
{

    //calculate all CLIQUES in all neighborhoods: TODO: maybe add storing slope/correlation
    std::cout << "STEP[1/6]:\tCalculate Cliques in all neighborhoods\n";
    //cliques in every neighborhood
    std::unordered_map<nodePtr, std::vector<std::vector<int>>> cliquesPerNeighborhood;
    //we calcualte all pairwise correlations already, store them here to not calc again later
    std::unordered_map<nodePtr, std::shared_ptr<GraphData>> neighborhoodCorrelations;
    //enqueue threads in pool for each thread to handle a separate neighborhood to calc all correlated cliques
    double count = 0;
    //create thread pool
    boost::asio::thread_pool pool_1(thread);
    for( nodePtr neighborhoodCenter : centralNeighborhoodPtrs)
    {
        boost::asio::post(pool_1, std::bind(&Neighborhood::detect_cliques_in_neighborhood, this, 
                          neighborhoodCenter, std::cref(correlationStrengthCutoff), minCliqueSize, printFoundCliquesPerNeighborhood,
                          std::ref(cliquesPerNeighborhood), std::ref(neighborhoodCorrelations),
                          centralNeighborhoodPtrs.size(), std::ref(count)));
    }
    pool_1.join();
    printProgress(1);
    std::cout << "\n";

    //cliquesVector is the final data structure storing all the correlated sets of features to process and from which to extract pairs
    std::cout << "STEP[2/6]:\tGet unique set of cliques\n";
    for(nodePtr neighborhoodCenter : centralNeighborhoodPtrs)
    {
        //for every clique in actual nHood
        for(std::vector<int> clique : cliquesPerNeighborhood.at(neighborhoodCenter))
        {
            if(std::find(cliquesVector.begin(), cliquesVector.end(), clique) == cliquesVector.end()) //if the clique does not exist yet
            {
                cliquesVector.push_back(clique);
            }
        }
    }

    //find consistent correlation sets
    std::cout << "STEP[4/6]:\tFilter consistent cliques that occur several neighborhoods\n";
    filter_consistent_correlation_sets(cliquesPerNeighborhood, minCliqueSize);

    //cliques might contain subset-cliques of other cliques - remove those
    std::cout << "STEP[3/6]:\tRemove duplicates/ subsets\n";
    if(cliquesVector.size() == 0)
    {
        std::cerr << "No Cliques detected with minimum clique size: " << minCliqueSize << "\n";
        std::cerr << "Exiting program: please reduce the number of a minimum correlation set or reduce minimum correlations\n";
        exit(EXIT_FAILURE);
    }
    remove_subsets(cliquesVector);

    std::cout << "STEP[5/6]:\tExtract pairwise correlations from correlated sets\n";
    //extract pairwise correlations from all found cliques.
    //pairwise correlations r extracted for cliquesVector data structure
    //results are correlations/ slopes for all pairs of features from cliques
    extract_pairs_from_correlation_sets(neighborhoodCorrelations);

    std::cout << "STEP[6/6]:\tCalculate Laplacian score\n";
    //make calculate_laplacian_score for all those slopes/ correlations
    calculate_laplacian_score(printFoundCliquesPerNeighborhood, thread);

    //PRINT ALL RESULT
    if(printFoundCliquesPerNeighborhood)
    {
        std::cout << "CORRELATIONS:\n";
        for (auto& it: correlationLaplacian) 
        {
            std::cout << inputDataOrigional.geneNames.at(it.first.first) << " " << inputDataOrigional.geneNames.at(it.first.second) << " " << it.second << "\n";
        }
        std::cout << "SLOPE:\n";
        for (auto& it: slopeLaplacian) 
        {
            std::cout << inputDataOrigional.geneNames.at(it.first.first) << " " << inputDataOrigional.geneNames.at(it.first.second) << " " << it.second << "\n";
        }
    }

}

void Neighborhood::create_neighborhood_graph(int knn)
{
    //create graph of neighborhoods: input is only the central nodes of neighborhoods
    std::shared_ptr<GraphData> scGraphData = std::make_shared<GraphData>(centralNeighborhoodPtrs, cellStateGenes, knn, &GraphIni::cell_similarity_graph_manhattan_nodes);
    //inout radius does not matter, but set bandwidth to neg. value to skip estimation
    neighborhoodGraph = std::make_shared<GraphHandler>(scGraphData, knn, 0, -1);
    neighborhoodGraph->create_graph();
}

Neighborhood::Neighborhood(const std::shared_ptr<const GraphData> scData, unsigned int neighborhoodNumber, 
                           unsigned int neighborhoodSize, int neighborhoodKNN,
                           const SingleCellData& inputData,
                           const std::vector<int>& cellStateGenes, const std::vector<int>& corrStateGenes, int permutations,
                           const double& corrSetAbundance, const uint correlatedSetMode) : 
                           neighborhoodSize(neighborhoodSize), inputDataOrigional(inputData),
                           cellStateGenes(cellStateGenes), corrStateGenes(corrStateGenes), permutations(permutations),
                           minimumCorrSetAbundance(corrSetAbundance), correlatedSetMode(correlatedSetMode)
{
    int cellIDRange = scData->number_of_nodes();
    //save all neighborhood IDs & node IDs making up neighborhoods
    std::vector<int> centralNodeIDs = neighborhoodCalculations::get_random_elements(neighborhoodNumber, cellIDRange);
    for(int centerNodeID : centralNodeIDs)
    {
        const nodePtr centerNode = scData->get_node_at(centerNodeID);
        centralNeighborhoodPtrs.push_back(centerNode);
        //std::vector<int> value = scData->get_adjacent_node_ids_knn(centerNode, neighborhoodSize);
        std::vector<int> value = scData->get_adjacent_node_ids_knn_kdsearch(centerNode);

        neighborhoods.insert(std::make_pair(centerNode, value));
    }

    //create the neighborhood graph (how to neighborhoods connect)
    // creates single-cell Graph and neighborhoods from that
    create_neighborhood_graph(neighborhoodKNN);
}