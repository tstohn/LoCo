#include <unordered_set>
#include <random>
#include <algorithm> 
#include <limits>
#include <set>
//#include <execution> // C++17 can improve runtime of std::sort with std::sort(std::execution::par, begin, end, ...), but we disabled it for R weaknesses to handle tbb
#include <atomic>

#include "Neighborhood.hpp"

namespace neighborhoodCalculations
{

    std::vector<int> get_random_elements(unsigned int numbers, unsigned int maxNum) 
    {
        if (numbers > maxNum)
            throw std::invalid_argument("Creating more neighbourhoods that cells is not allowed!");

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


// old function to write all the permuted laplacian scores: not logner done
void Neighborhood::write_shuffled_laplacians(const std::string& outFile, const std::string& prefix)
{
    std::ofstream outputFile;
    std::string folder_prefix = outFile + "/" + prefix;

    // CORRELATION LAPLACIANS
    std::string outputCorr = folder_prefix + "_shuffledCorr.tsv" ;
    std::remove(outputCorr.c_str());
    outputFile.open(outputCorr, std::ofstream::app);

    // collect all localPairs
    std::vector<std::pair<int,int>> localPairs;
    for(const auto& it : shuffledCorrLaplacians)
    {
        localPairs.push_back(it.first);
    }

    if(localPairs.empty())
    {
        outputFile.close();
        return;
    }

    int permutations = shuffledCorrLaplacians.at(localPairs[0]).size();

    // HEADER
    outputFile << "Permutation";
    for(const std::pair<int,int>& pair : localPairs)
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

        for(const std::pair<int,int>& pair : localPairs)
        {
            outputFile << "\t" << shuffledCorrLaplacians.at(pair).at(p);
        }

        outputFile << "\n";
    }

    outputFile.close();

}

void Neighborhood::write_results_to_file(const std::string& outFile, const std::string& prefix, bool calcSets)
{
    std::ofstream outputFile;   

    /////////////////////////////
    // 1.) LAPLCACIAN SCORES FILE
    /////////////////////////////
    // Create an array of indices [0, 1, 2, ..., N-1]
    std::vector<size_t> indices(laplacianScores.pairNames.size());
    std::iota(indices.begin(), indices.end(), 0);

    // Sort the indices based on the Laplacian values (L) in ascending order
    std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
        
        bool isNanA = std::isnan(laplacianScores.L[a]);
        bool isNanB = std::isnan(laplacianScores.L[b]);
        
        // Safety: Push NaNs to the very bottom of the file
        if (isNanA && isNanB) return false; 
        if (isNanA) return false;
        if (isNanB) return true;
        
        // Ascending order: smallest L goes first
        return laplacianScores.L[a] < laplacianScores.L[b];
    });
    //1FILE with annotations: corrpair - LaplacianCorrscore - LaplacianSlopeScore - origionalCliqueItIsFrom
    std::string outputLaplace = outFile + "/" + prefix + "_laplacian.tsv";
    std::remove(outputLaplace.c_str());
    outputFile.open (outputLaplace, std::ofstream::app);
    //write HEADER
    outputFile << "FeaturePair" << "\t" << "CorrelationScore" << "\t" << "p_CorrelationScore"; 
    if(calcSets){ outputFile << "\t" << "Feauture Set";} 
    outputFile << "\n";
    //write LINES
    for(size_t i = 0; i < indices.size(); ++i)
    {
        size_t idx = indices[i]; // Map loop iteration to the correct sorted row

        // write PAIR << Laplacian << p_value
        outputFile << laplacianScores.pairNames.at(idx) << "\t" << 
        std::to_string(laplacianScores.L[idx]) << "\t" << 
        std::to_string(laplacianScores.p_values[idx]);

        if(calcSets)
        {
            outputFile << "\t" << neighbourhoodCorrs.featureSetString[idx];
        }

        outputFile << "\n";
    }
    outputFile.close();

    /////////////////////////////
    // CORRELATIONS:  for NeighborhoodID * CorrpAir: value == corr
    /////////////////////////////
    //CORRELATIONS FILE
    std::string outputCorr = outFile + "/" + prefix + "_correlations.tsv";
    std::remove(outputCorr.c_str());
    outputFile.open (outputCorr, std::ofstream::app);
    //write HEADER: first col has all featurePair names and then every col-header is a neighbourhood (inverse to origional version)
    outputFile << "FeaturePair";
    for(const auto& neighbourhoddPtr : neighbourhoodCorrs.nodePtrList)
    {
        outputFile << "\t" << neighbourhoddPtr->get_name();
    }
    outputFile << "\n";

    for(int featurePairIdx =0; featurePairIdx < neighbourhoodCorrs.pairNames.size(); ++featurePairIdx) 
    {

        std::string featurePair = neighbourhoodCorrs.pairNames.at(featurePairIdx);
        //get all correlation values for all neighbourhood for this featurePair
        const double* pairCorrs = neighbourhoodCorrs.pairToListOfAllNCorrs.get_row(featurePairIdx);

        outputFile << featurePair;
        for(int Nidx = 0; Nidx < neighbourhoodNum; ++Nidx)
        {
            outputFile << "\t" << pairCorrs[Nidx];
        }
        outputFile << "\n";
    }
    outputFile << "\n";
    outputFile.close();
    
    /////////////////////////////
    // 3.)  File with annotations: NeighborhoodID - centralNodeDimensions[start ... end]
    /////////////////////////////
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

    /////////////////////////////
    // 4.) File with cellIDs in each neighborhood: NeighborhoodID in header - all cellIDs in rows for each neighborhood columns
    /////////////////////////////
    std::string outputNeighborhoodCells = outFile + "/" + prefix + "_cells.tsv";
    std::remove(outputNeighborhoodCells.c_str());
    outputFile.open (outputNeighborhoodCells, std::ofstream::app);
    //write HEADER (all neigborhood names)
    const std::vector<nodePtr> neighborHoodPtrVector = neighborhoodGraph->get_all_nodes();
    for(size_t neiborhoodID = 0; neiborhoodID < neighborHoodPtrVector.size(); ++neiborhoodID)
    {
        outputFile << neighborHoodPtrVector.at(neiborhoodID)->get_name();
        if(neiborhoodID < (neighborHoodPtrVector.size()-1)){outputFile << "\t";}
    }  
    outputFile << "\n";
    //write lines
    //iterate over cellIDs (0 - <number cells in neighborhood>), for every id write out the cell with this id for every neighborhood
    for(size_t cellID = 0; cellID < neighborhoodSize; ++cellID)
    {
        //for every neighborhood (in columns)
        for(size_t neiborhoodID = 0; neiborhoodID < neighborHoodPtrVector.size(); ++neiborhoodID)
        {
            outputFile << neighborhoods.at(neighborHoodPtrVector.at(neiborhoodID)).at(cellID);
            if(neiborhoodID < (neighborHoodPtrVector.size()-1)){outputFile << "\t";}
        }
        outputFile << "\n";
    }
    outputFile.close();
}

void Neighborhood::calculate_correlations_variance(unsigned int numberNodes, 
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
        if(!std::isnan(correlationTmp))
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
            if(!std::isnan(correlationTmp))
            {
                variance += pow((correlationTmp-mean), 2);
            }
        }

        variance /= nonNanNodes;
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
    unsigned int nonNanNodes = 0;
    for(const nodePtr& node: neighborhoodGraph->get_all_nodes())
    {
        double slopeTmp = corrResult.at(node).slopeResult.at(correlationpair);
        if(!std::isnan(slopeTmp))
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
            if(!std::isnan(slopeTmp))
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

            if( !std::isnan(featureNodeACorr) && !std::isnan(featureNodeBCorr))
            {
                corrWeightSum += weight * pow((featureNodeACorr - featureNodeBCorr), 2); //square the result
            }
            if( !std::isnan(featureNodeASlope) && !std::isnan(featureNodeBSlope))
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
void Neighborhood::calculate_laplacian_score(int threads)
{   

    //calculate once the variance that we observe in all correlations/ slopes
    unsigned int numberNodes = neighborhoodGraph->get_all_nodes().size();

    //CORRELATION VARIANCE
    LOCO_OUT << "\tSTEP[6a]:\tCalculate variance of smooth Correlations\n";
    std::unordered_map<const std::pair<int, int>, double, pair_hash> corrVariance;
    //boost::asio::thread_pool pool_corr(threads);
    ThreadPool pool_corr(threads);

    double count = 0;
    for(const std::pair<int, int>& correlationpair : pairs)
    {
        //boost::asio::post(pool_corr, std::bind(&Neighborhood::calculate_correlations_variance, this, 
        //                  numberNodes, print, std::cref(correlationpair), std::ref(corrVariance), pairs.size(), std::ref(count)));
        pool_corr.enqueue([
            this,                       
            numberNodes,                
            &correlationpair,           
            &corrVariance,              
            pairCount = pairs.size(),   
            &count                      
        ]() {
            calculate_correlations_variance(
                numberNodes,
                correlationpair,
                corrVariance,
                pairCount,
                count
            );
        });
    }
    pool_corr.wait_for_tasks();
    printProgress(1);
    LOCO_OUT << "\n";

    //SLOPE VARIANCE
    LOCO_OUT << "\tSTEP[6b]:\tCalculate variance of smooth Slopes\n";
    std::unordered_map<const std::pair<int, int>, double, pair_hash> slopeVariance;
    //boost::asio::thread_pool pool_slope(threads);
    ThreadPool pool_slope(threads);
    count = 0;
    for(const std::pair<int, int>& correlationpair : pairs)
    {
        //boost::asio::post(pool_slope, std::bind(&Neighborhood::calculate_slopes, this, 
        //                  numberNodes, std::cref(correlationpair), std::ref(slopeVariance), pairs.size(), std::ref(count)));
        pool_slope.enqueue([
            this,                       
            numberNodes,                
            &correlationpair,           
            &slopeVariance,              
            pairCount = pairs.size(),   
            &count                      
        ]() {
            calculate_slopes(
                numberNodes,
                correlationpair,
                slopeVariance,
                pairCount,
                count
            );
        });
    }
    pool_slope.wait_for_tasks();
    printProgress(1);
    LOCO_OUT << "\n";

    //FULL LAPLACIAN SCORE
    LOCO_OUT << "\tSTEP[6c]:\tCalculate Laplacian Scores\n";
    //boost::asio::thread_pool pool_lapl(threads);
    ThreadPool pool_lapl(threads);

    count = 0;
    for(const std::pair<int, int>& pair : pairs)
    {
        //boost::asio::post(pool_lapl, std::bind(&Neighborhood::laplacian_score, this, 
        //                  std::cref(pair), std::cref(corrVariance), std::cref(slopeVariance),
        //                  pairs.size(), std::ref(count)));
        pool_lapl.enqueue([
            this,                       
            &pair,                
            &corrVariance,           
            &slopeVariance,              
            pairCount = pairs.size(),   
            &count                      
        ]() {
            laplacian_score(
                pair,                
                corrVariance,           
                slopeVariance,              
                pairCount,   
                count  
            );
        });
    }
    pool_lapl.wait_for_tasks();
    printProgress(1);
    LOCO_OUT << "\n";

    //calcualte significance: calcualte laplacian score for randomly shuffled features across N
    LOCO_OUT << "\tSTEP[6d]:\tCalculate Significance: Laplacian Scores after shuffling correlations\n";
    LOCO_OUT << "\t\t We keep the KNN graph and randomly reassign all correlations/slops of a neighborhood to a different neighborhood\n";
    count = 0;
    //firstly write the correlation/slope results to a vector that is easily shufflable
    std::vector<CorrelationPropagationResult> vectorizedResults;
    vectorizedResults.reserve(corrResult.size());
    for (const auto& [node, result] : corrResult)
    {
        vectorizedResults.push_back(result);
    }
    ThreadPool pool_shuffle(threads);
    for(const std::pair<int, int>& pair : pairs)
    {
        //boost::asio::post(pool_shuffle, std::bind(&Neighborhood::laplacian_significance, this, 
        //                  std::cref(pair), std::cref(corrVariance), std::cref(slopeVariance),
        //                 std::cref(vectorizedResults),
        //                  pairs.size(), std::ref(count)));
        pool_shuffle.enqueue([
            this,                       
            &pair,                
            &corrVariance,           
            &slopeVariance,   
            &vectorizedResults,           
            pairCount = pairs.size(),   
            &count                      
        ]() {
            laplacian_significance(
                pair,                
                corrVariance,           
                slopeVariance, 
                vectorizedResults,             
                pairCount,   
                count  
            );
        });
    }
    pool_shuffle.wait_for_tasks();
    printProgress(1);
    LOCO_OUT << "\n";

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
                    std::vector<std::vector<std::string>> stringCliquesVector;
                    stringCliquesVector.push_back(features);
                    pairToClique.insert(std::pair<const std::pair<int, int>, std::vector<std::vector<std::string>>>(tmpPair, stringCliquesVector));
                }
                if(std::find(pairs.begin(), pairs.end(), tmpPair) == pairs.end()) //if the pair was not stored yet
                {
                    pairs.push_back(tmpPair);
                }
            }
        }
    }

    //pre-compute ranks for spearman correlation
    if(correlationType == "spearman")
    {
        for(auto const& [center, tmpData] : neighborhoodCorrelations) 
        {
            for(int i = 0; i < tmpData->number_of_nodes(); ++i) 
            {
                tmpData->get_node_at(i)->compute_ranks(); 
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

            double corr;
            if(correlationType == "spearman")
            {
                corr = calculate_correlation_coefficient(featureNodeA->all_ranked_values(), featureNodeB->all_ranked_values());
            }
            else if(correlationType == "pearson")
            {
                corr = calculate_correlation_coefficient(featureNodeA->all_values(), featureNodeB->all_values());
            }

            //if (std::isnan(corr)) 
            //{
            //    std::cout << "WARNING: Correlation value is NaN" << " between features: " << featureNodeA->get_name() << " and " << featureNodeB->get_name() << "\n";    
            //}
                
            tmpResult.correlationResult.insert(std::pair< std::pair<int, int>, double >(pair, corr));

            //store all SLOPES for those pairs (have to calculate new)
            const double slope = calculate_slope(featureNodeA->all_values(), featureNodeB->all_values());
            tmpResult.slopeResult.insert(std::pair< std::pair<int, int>, double >(pair, slope));
        }

        corrResult.insert(std::pair<nodePtr, CorrelationPropagationResult>(neighborhoodCenter, tmpResult));
    }
    
}

void Neighborhood::detect_cliques_in_neighborhood(nodePtr neighborhoodCenter, const double& correlationStrengthCutoff, int minCliqueSize,
                                                  std::unordered_map<nodePtr, std::vector<std::vector<int>>>& cliquesPerNeighborhood, 
                                                  std::unordered_map<nodePtr, std::shared_ptr<GraphData>>& neighborhoodCorrelations,
                                                  int totalCount, double& currentCount)
{
        //calculate all correlations/ slopes of them
        SingleCellData inputDataTmp = filter_singleCelldata(inputDataOrigional, neighborhoods.at(neighborhoodCenter));

        //for the protein graph calculate brute force all distances instead of reading form a KD-tree
        // this way we can filter correlations based on value
        unsigned int proteinGraphKnn = 0;
        std::shared_ptr<GraphData> correlationData;
        if(correlationType == "spearman")
        {
            correlationData = std::make_shared<GraphData>(inputDataTmp, corrStateGenes, proteinGraphKnn, &GraphIni::protein_correlation_graph_spearman);
        }
        else if(correlationType == "pearson")
        {
            correlationData = std::make_shared<GraphData>(inputDataTmp, corrStateGenes, proteinGraphKnn, &GraphIni::protein_correlation_graph_pearson);
        }
        
        GraphHandler corrGraphBuilder = GraphHandler(correlationData, 0, correlationStrengthCutoff, -1);
        corrGraphBuilder.create_graph();

        //results are stored here afterwards
        std::vector<std::vector<int>> cliqueVectorRaw;
        if(correlatedSetMode == 1)
        {
            corrGraphBuilder.calculate_connected_components(cliqueVectorRaw, minCliqueSize);
        }
        else if(correlatedSetMode == 0)
        {
            corrGraphBuilder.calculate_fully_connected_components(cliqueVectorRaw, minCliqueSize);
        }
        else if(correlatedSetMode >= 2)
        {
            corrGraphBuilder.calculate_min_edge_connected_components(cliqueVectorRaw, minCliqueSize, correlatedSetMode);
        }
        else
        {
            LOCO_OUT << "Invalid mode for detection of correlated sets: fallback to connected component -> mode 0\n";
            corrGraphBuilder.calculate_connected_components(cliqueVectorRaw, minCliqueSize);
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

// Helper to create a fast bitmask for pruning
uint64_t create_signature(const std::vector<int>& vec) {
    uint64_t mask = 0;
    for (int x : vec) mask |= (1ULL << (static_cast<uint32_t>(x) % 64));
    return mask;
}

void Neighborhood::filter_consistent_correlation_sets_sota(
        const std::unordered_map<nodePtr, std::vector<std::vector<int>>>& cliquesPerNeighborhood,
        int minCliqueSize)
{
    // 1. GENERATE CANDIDATES: all possible subsets of found sets of allowed size
    std::set<std::vector<int>> candidateSets;
    for(const auto& clique : cliquesVector) {
        std::vector<int> sortedClique = clique;
        std::sort(sortedClique.begin(), sortedClique.end());
        auto subsets = neighborhoodCalculations::generate_subsets(sortedClique, minCliqueSize);
        for(auto& s : subsets){ candidateSets.insert(std::move(s));}
    }

    // 2. PRE-PROCESS NEIGHBORHOODS (The Speed Boost)
    // We sort everything once and pre-calculate bitmasks
    struct IndexedNeighborhood {
        std::vector<std::vector<int>> sortedCliques;
        std::vector<uint64_t> masks;
    };
    
    std::vector<IndexedNeighborhood> processedData;
    processedData.reserve(cliquesPerNeighborhood.size());

    for (auto const& [node, cliques] : cliquesPerNeighborhood) {
        IndexedNeighborhood idx;
        for (auto c : cliques) {
            std::sort(c.begin(), c.end());
            idx.masks.push_back(create_signature(c));
            idx.sortedCliques.push_back(std::move(c));
        }
        processedData.push_back(std::move(idx));
    }

    const int neighborhoodCount = processedData.size();
    const int required = std::ceil(minimumCorrSetAbundance * neighborhoodCount);
    std::vector<std::vector<int>> result;

    // check all candidates
    for(const auto& candidate : candidateSets) {
        uint64_t candMask = create_signature(candidate);
        int count = 0;

        for(const auto& neighborhood : processedData) {
            bool found = false;
            for(size_t i = 0; i < neighborhood.sortedCliques.size(); ++i) {
                // PRUNING STEP: Bitmask check is nearly free
                if ((candMask & neighborhood.masks[i]) != candMask) continue;
                
                // Size check: Candidate can't be a subset if it's bigger
                if (candidate.size() > neighborhood.sortedCliques[i].size()) continue;

                // FINAL CHECK: std::includes is the fastest way to check subset on sorted vectors
                if (std::includes(neighborhood.sortedCliques[i].begin(), neighborhood.sortedCliques[i].end(),
                                  candidate.begin(), candidate.end())) {
                    found = true;
                    break;
                }
            }
            if(found) count++;
            if(count >= required) break; // Early exit for this candidate
        }

        if(count >= required) result.push_back(candidate);
    }

    cliquesVector = std::move(result);
}

//calcualte all apirwise corrs within a N and push it into shared tempCalculationMatrix
void Neighborhood::calculate_correlations_for_N(nodePtr neighborhoodCenter, size_t neighborhood_idx, 
                                                const double& corrThreshold, FlatMatrix& tempCalculationMatrix,
                                                const std::vector<std::pair<int, int>>& tmpAllPairs , 
                                                std::vector<std::atomic<int>>& tmpCorrelationCountAboveThreshold, 
                                                std::atomic<int>& currentCount)
{
    //1.) subset the data used for this neighbourhood

    // add filtering step to choose only corrGEnes
    std::vector<std::vector<double>> pointCloud = return_filtered_point_cloud(inputDataOrigional, neighborhoods.at(neighborhoodCenter), corrStateGenes);
    
    // make rank computation independent of nodes: do it on the transpose of raw data to get all features quickly
    //2.) compute ranks of values for spearman correlation
    std::vector<std::vector<double>> featureMajorData = transpose_to_feature_major(pointCloud);
    if(correlationType == "spearman")
    {
        for (size_t f = 0; f < featureMajorData.size(); ++f) 
        {
            rankify(featureMajorData[f]); 
        }
    }

    // 3.) Calculate correlation coefficients and update FlatMatrix safely!
    for(size_t p = 0; p < tmpAllPairs.size(); ++p)
    {
        int featA = tmpAllPairs[p].first;
        int featB = tmpAllPairs[p].second;

        // Calculate correlation (assuming your function takes two vectors)
        double corr = calculate_correlation_coefficient(featureMajorData[featA], featureMajorData[featB]);

        // Write directly to the matrix lock-free! (Using the neighborhood_idx as the row)
        tempCalculationMatrix(neighborhood_idx, p) = corr;

        // If above threshold, safely increment the global counter without a lock
        if (std::abs(corr) >= corrThreshold) 
        {
            // std::atomic safely handles multiple threads incrementing this at once
            tmpCorrelationCountAboveThreshold[p].fetch_add(1, std::memory_order_relaxed);
        }
    }

    //update all pairwise corrs for this N: add it to FlatMatrix
    currentCount++;
}

void Neighborhood::calculate_pair_variance(size_t pair_idx)
{
    // 1. Grab a zero-copy pointer to the entire row of correlations across all Ns
    size_t numNeighbourhoods = 0;
    const double* corrsAcrossN = neighbourhoodCorrs.get_correlations_across_N_for_pair(pair_idx, numNeighbourhoods);

    size_t nonNanNodes = 0;
    double mean = 0.0;

    // 2. FIRST PASS: Calculate the Mean and count valid entries
    for (size_t s = 0; s < numNeighbourhoods; ++s)
    {
        double correlationTmp = corrsAcrossN[s];
        
        // Ensure the value is valid (Not NaN and not the default -2.0 FlatMatrix value)
        if (!std::isnan(correlationTmp) && correlationTmp != -2.0)
        {
            mean += correlationTmp;
            nonNanNodes++;
        }
    }

    double variance = 0.0;

    // 3. Check for the 10% threshold
    if (nonNanNodes >= numNeighbourhoods / 10)
    {
        mean /= nonNanNodes;

        // SECOND PASS: Calculate deviations from the mean
        for (size_t s = 0; s < numNeighbourhoods; ++s)
        {
            double correlationTmp = corrsAcrossN[s];
            if (!std::isnan(correlationTmp) && correlationTmp != -2.0)
            {
                double deviation = correlationTmp - mean;
                variance += deviation * deviation; // Fast explicit multiplication instead of pow()
            }
        }

        variance /= nonNanNodes;
    }
    else 
    {
        variance = 0.0; // Left explicitly at zero per 10perc-rule
    }

    // 4. LOCK-FREE WRITE: Safe because every thread owns a unique pair_idx slot!
    laplacianScores.variances[pair_idx] = variance;
}

// neighbourhoods have order as in centralNeighborhoodPtrs
// corr-pairs as in tmpAllPairs which is then written to Flatmatrix
void Neighborhood::step_2_calculate_correlation(const double& corrThreshold, const int threads)
{

    //1.) INITIALIZE DATA STRUCTURES

    // a) Determine your feature size
    size_t geneSize = inputDataOrigional.geneNames.size();
    if (!corrStateGenes.empty()) 
    {
        geneSize = corrStateGenes.size();
    }
    size_t totalPairs = (geneSize * (geneSize - 1)) / 2;

    //b) temporary TRANSPOSE matrix to quickly iterae over N and calc corrs for every pair
    //(later transpose: neighbourhoodCorrs has outer vector for corr-pairs to quickly calc laplacian for corr-pairs across N)
    std::vector<std::pair<int, int>> tmpAllPairs; //later keep only the pairs with x% above corr-threshold
    tmpAllPairs.reserve(totalPairs);
    for (int i = 0; i < static_cast<int>(geneSize); ++i) 
    {
    for (int j = i + 1; j < static_cast<int>(geneSize); ++j) 
        {
            // emplace_back constructs the pair directly in place, 
            // which is slightly faster than push_back({i, j})
            tmpAllPairs.emplace_back(i, j); 
        }
    }
    FlatMatrix tempCalculationMatrix(neighbourhoodNum, totalPairs);

    //c) create structures
    // temporary structures to count correlation-pair with >threshold correlations
    //count in which feature-pair correlation is above threshold
    std::vector<std::atomic<int>> tmpCorrelationCountAboveThreshold(totalPairs); 

    //3.) calcualte all correlations in threads
    std::atomic<int> completedNCount{0};
    size_t neighborhood_idx = 0;
    {
        ThreadPool pool_correlations(threads);
        for( nodePtr neighborhoodCenter : centralNeighborhoodPtrs)
        {
            pool_correlations.enqueue([
                this,                              
                neighborhoodCenter,
                neighborhood_idx,                        
                &corrThreshold, 
                &tempCalculationMatrix,  
                &tmpAllPairs,                                                
                &tmpCorrelationCountAboveThreshold,    
                &completedNCount                        
            ]() {
                calculate_correlations_for_N(
                    neighborhoodCenter,      
                    neighborhood_idx,                  
                    corrThreshold,                              
                    tempCalculationMatrix,
                    tmpAllPairs,     
                    tmpCorrelationCountAboveThreshold,               
                    completedNCount  
                );
            });

            neighborhood_idx++; // Increment for the next loop iteration
        }

        // 2. The Main R Thread prints the progress while workers work!
        int last_val = -1;

        while (completedNCount < neighbourhoodNum) 
        {
            double percentage = (double)completedNCount.load() / neighbourhoodNum;
            int val = (int)(percentage * 100);

            // Only call the print function if the integer percentage actually changed
            if (val != last_val) 
            {
                printProgressMainThread(percentage);
                last_val = val;
            }

            // Sleep briefly so the main thread doesn't hog the CPU while waiting
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
        }

        // 3. Guarantee it prints 100% and a newline at the very end
        pool_correlations.wait_for_tasks();
    }
    printProgressMainThread(1.0);
    LOCO_OUT << "\n";

    //4.) Create result for correlations (per N)
    // get only corrs that are x% over threshood and safe them in neighbourhoodCorrs
    // do this inside the init function of the neighbourhoodCorrelations and instead of the matrix transpose do a 
    // in-place writing into the new amtrix (transposed) only for coors > threshold
    const int minRequired = std::ceil(minimumCorrSetAbundance * neighbourhoodNum);    
    neighbourhoodCorrs.init(tempCalculationMatrix, 
                      tmpAllPairs, 
                      tmpCorrelationCountAboveThreshold, 
                      minRequired, 
                      centralNeighborhoodPtrs, 
                      inputDataOrigional.geneNames,
                      corrStateGenes);

}

void Neighborhood::calculate_laplacian_score_for_pair(const int featurePairIdx)
{   
    double corrWeightSum = 0;
    const double* pairCorrs = neighbourhoodCorrs.pairToListOfAllNCorrs.get_row(featurePairIdx);

    //safely check variance: we set vairance of corrs to NAN if less than 10% of neighbourhoods had a valid corr defined
    double variance = laplacianScores.variances[featurePairIdx];
    if (variance == 0.0) 
    {
        laplacianScores.L[featurePairIdx] = std::numeric_limits<double>::quiet_NaN();
        return;
    }

    for(unsigned int i = 0; i < (neighbourhoodNum-1); ++i)
    {
        for(unsigned int j = i+1; j < neighbourhoodNum; ++j)
        {
            double weight = neighborhoodGraph->get_edge_weight_between_nodes(i, j);
            if(weight == 0){continue;}

            //actual feature values for nodes
            double featureNodeACorr = pairCorrs[i];
            double featureNodeBCorr = pairCorrs[j];

            if( !std::isnan(featureNodeACorr) && !std::isnan(featureNodeBCorr))
            {
                double diff = featureNodeACorr - featureNodeBCorr;
                corrWeightSum += weight * (diff * diff); //square the result
            }
        }
    }

    corrWeightSum /= laplacianScores.variances.at(featurePairIdx);
    laplacianScores.L[featurePairIdx] = corrWeightSum;
}

void Neighborhood::laplacian_significance_for_pair(size_t pair_idx, const std::vector<Edge>& edges, std::atomic<int>& currentCount)
{
    // Quick lookups
    double observedL = laplacianScores.L[pair_idx];
    double variance = laplacianScores.variances[pair_idx];

    // Safety Check: If L is NaN (due to 0 variance earlier), p-value is also NaN
    if (std::isnan(observedL) || variance == 0.0)
    {
        laplacianScores.p_values[pair_idx] = std::numeric_limits<double>::quiet_NaN();
        currentCount++;
        return;
    }

    // Zero-Copy access to feature values
    const double* pairCorrs = neighbourhoodCorrs.pairToListOfAllNCorrs.get_row(pair_idx);

    // create a vector of numbers from 0 to number of N
    thread_local std::mt19937 rng(std::random_device{}());
    std::vector<int> perm(neighbourhoodNum);
    std::iota(perm.begin(), perm.end(), 0);
    double p_count = 0.0;

    // Permutation Loop
    for(int p = 0; p < permutations; ++p)
    {
        //shuffle this vector of numbers: use these as indices for shuffled neighbourhood correlations
        std::shuffle(perm.begin(), perm.end(), rng);

        double corrSum = 0.0;

        // Iterate over the pre-built edge list (Massive speedup!)
        for(const Edge& e : edges)
        {
            // Map the graph's fixed indices to the shuffled correlation data
            double corrA = pairCorrs[perm[e.i]];
            double corrB = pairCorrs[perm[e.j]];

            if(!std::isnan(corrA) && !std::isnan(corrB))
            {
                double diff = corrA - corrB;
                corrSum += e.w * (diff * diff);
            }
        }

        double shuffledL = corrSum / variance;

        // Compare against observed
        if(shuffledL <= observedL)
        {
            p_count++;
        }
    }

    // Lock-free write to pre-allocated results
    laplacianScores.p_values[pair_idx] = p_count / permutations;
    currentCount++;
}

void Neighborhood::step_3_calculate_laplacian_score(const int threads)
{
    laplacianScores.pairNames = neighbourhoodCorrs.pairNames;

    // 1.) pre-calcualte variance for significance-calculations later
    // Get total number of filtered pairs from your global/member structure
    size_t totalPairs = laplacianScores.pairNames.size();
    // pre-allocate laplacian result vector for thread-safe writing
    laplacianScores.variances.resize(totalPairs, 0.0);
    // Initialize the atomic counter and the thread pool
    ThreadPool pool_variance(threads);
    // Enqueue tasks by index (Blazing fast, zero allocations inside the loop)
    for (size_t p = 0; p < totalPairs; ++p)
    {
        pool_variance.enqueue([this, p]() 
        {
            this->calculate_pair_variance(p);
        });
    }
    // Block main thread until all background workers finish their math
    pool_variance.wait_for_tasks();
    
    // 2.) calcualte laplacian score
    LOCO_OUT << "\t Calculate Laplacian Score for found correlations\n";
    // Pre-allocate sizes inside your global/member LaplacianResults struct
    laplacianScores.L.resize(totalPairs, 0.0);         // Holds final Laplacian scores
    ThreadPool pool_L(threads);
    {
        for (size_t p = 0; p < totalPairs; ++p)
        {
            pool_L.enqueue([this, p]() {
                this->calculate_laplacian_score_for_pair(p);
            });
        }
    }
    pool_L.wait_for_tasks();

    // 3.) calcualte significance values (only print status here)
    LOCO_OUT << "\t Calculate significance for Laplacian Score across " << std::to_string(permutations) << " permutations\n";
    // Build the global edge list exactly ONCE
    std::vector<Edge> globalEdges;
    globalEdges.reserve(neighbourhoodNum * 5); // guess due to KNN graph, will be more in reality since its not a directed graph
    for(unsigned int i = 0; i < neighbourhoodNum - 1; ++i)
    {
        for(unsigned int j = i + 1; j < neighbourhoodNum; ++j)
        {
            double w = neighborhoodGraph->get_edge_weight_between_nodes(i, j);
            if(w != 0.0)
            {
                globalEdges.push_back({(int)i, (int)j, w});
            }
        }
    }
    // Pre-allocate p-values
    laplacianScores.p_values.resize(totalPairs, 0.0);
    // Initialize the atomic counter
    std::atomic<int> completedPairsCount{0};
    {
        ThreadPool pool_signifiance(threads);
        // Enqueue all tasks
        for (size_t p = 0; p < totalPairs; ++p)
        {
            // Capture the atomic counter by reference
            pool_signifiance.enqueue([this, p, &globalEdges, &completedPairsCount]() 
            {
                this->laplacian_significance_for_pair(p, globalEdges, completedPairsCount);
            });
        }
        // MAIN THREAD PROGRESS LOOP
        // The main thread waits here and updates the UI instead of blocking entirely
        while (completedPairsCount < totalPairs)
        {
            // Calculate percentage
            double perc = static_cast<double>(completedPairsCount) / totalPairs;
            printProgressMainThread(perc);
            // Put the main thread to sleep for 100ms so it doesn't waste CPU cycles 
            // while checking the counter.
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        // 6. Ensure everything is cleanly wrapped up before destroying the pool
        pool_signifiance.wait_for_tasks(); 
    }
    // Print final 100% state and a newline
    printProgressMainThread(1.0);
    LOCO_OUT << "\n";

}


/*
OLD CODE AS EXMAPLE

      //calculate all correlations/ slopes of them
        SingleCellData inputDataTmp = filter_singleCelldata(inputDataOrigional, neighborhoods.at(neighborhoodCenter));

        //for the protein graph calculate brute force all distances instead of reading form a KD-tree
        // this way we can filter correlations based on value
        unsigned int proteinGraphKnn = 0;
        std::shared_ptr<GraphData> correlationData;
        if(correlationType == "spearman")
        {
            correlationData = std::make_shared<GraphData>(inputDataTmp, corrStateGenes, proteinGraphKnn, &GraphIni::protein_correlation_graph_spearman);
        }
        else if(correlationType == "pearson")
        {
            correlationData = std::make_shared<GraphData>(inputDataTmp, corrStateGenes, proteinGraphKnn, &GraphIni::protein_correlation_graph_pearson);
        }
        
        GraphHandler corrGraphBuilder = GraphHandler(correlationData, 0, correlationStrengthCutoff, -1);
        corrGraphBuilder.create_graph();

        //results are stored here afterwards
        std::vector<std::vector<int>> cliqueVectorRaw;
        if(correlatedSetMode == 1)
        {
            corrGraphBuilder.calculate_connected_components(cliqueVectorRaw, minCliqueSize);
        }
        else if(correlatedSetMode == 0)
        {
            corrGraphBuilder.calculate_fully_connected_components(cliqueVectorRaw, minCliqueSize);
        }
        else if(correlatedSetMode >= 2)
        {
            corrGraphBuilder.calculate_min_edge_connected_components(cliqueVectorRaw, minCliqueSize, correlatedSetMode);
        }
        else

*/

void Neighborhood::step_4_calculate_feature_sets(int minFeatureSetSize)
{
    // create graph with features as nodes
    // use the second GraphHandler initializer: it creates directly a graph from nodes/ edges
    // input are the feature names and the pairList
    GraphHandler consensusCorrelationGraph = GraphHandler(neighbourhoodCorrs.featureNames, neighbourhoodCorrs.pairList);
    std::vector<std::vector<int>> featureSetsRaw;

    //call subgraph algorithm on this data
    if(correlatedSetMode == 1)
    {
        consensusCorrelationGraph.calculate_connected_components(featureSetsRaw, minFeatureSetSize);
    }
    else if(correlatedSetMode == 0)
    {
        consensusCorrelationGraph.calculate_fully_connected_components(featureSetsRaw, minFeatureSetSize);
    }
    else if(correlatedSetMode >= 2)
    {
        consensusCorrelationGraph.calculate_min_edge_connected_components(featureSetsRaw, minFeatureSetSize, correlatedSetMode);
    }
    else
    {
        LOCO_OUT << "Invalid mode for detection of correlated sets: fallback to connected component -> mode 0\n";
        consensusCorrelationGraph.calculate_connected_components(featureSetsRaw, minFeatureSetSize);
    } 

    //iterate through the list of feature-pairs: for each pair check in which featureSetsRaw they are and create the string for it

    size_t numPairs = neighbourhoodCorrs.pairList.size();

    // 1. Setup fast lookup map (Pair -> Index in pairList)
    std::unordered_map<std::pair<int, int>, int, pair_hash> pairToIndexMap;
    pairToIndexMap.reserve(numPairs);    
    for (size_t i = 0; i < numPairs; ++i) {
        pairToIndexMap[neighbourhoodCorrs.pairList[i]] = i;
    }

    // 2. Storage for mapping pair index -> lists of feature sets that contain this pair
    std::vector<std::vector<std::vector<int>>> pairToSets(numPairs);

    // 3. The Smart Iteration: Loop sets, build pairs, map to index
    for (const auto& featureSet : featureSetsRaw) {
        
        if (featureSet.size() < 2) continue; // Skip sets too small to form a pair

        for (size_t i = 0; i < featureSet.size(); ++i) {
            for (size_t j = i + 1; j < featureSet.size(); ++j) {
                
                int u = featureSet[i];
                int v = featureSet[j];

                // Ensure the pair is always ordered [smaller, larger] to match pairList
                if (u > v) {
                    std::swap(u, v);
                }

                // Check if this generated pair exists in our main pairList
                auto it = pairToIndexMap.find({u, v});
                if (it != pairToIndexMap.end()) {
                    int pairIndex = it->second;
                    // Add this entire featureSet to the pair's storage
                    pairToSets[pairIndex].push_back(featureSet);
                }
            }
        }
    }

    // 4. Create the final vector of strings aligned exactly with pairList
    neighbourhoodCorrs.featureSetString.resize(numPairs);
    
    for (size_t i = 0; i < numPairs; ++i) {
        // Call your custom function, passing the collected sets for this specific pair
        neighbourhoodCorrs.featureSetString[i] = consensusCorrelationGraph.return_named_feature_set(pairToSets[i]);
    }

    // Now finalFeatureSetStrings[i] contains the correctly formatted string 
    // for the pair at neighbourhoodCorrs.pairList[i]
    
}

//centralneighbourhoodPtrs: all the neighbourhoodsPTrs, sotring all cells in enighbourhoods etc.
void Neighborhood::calculate_correlation_propagation(double correlationStrengthCutoff, int minFeatureSetSize, const bool calcSets, int threads)
{

    LOCO_OUT << "STEP 2";
    LOCO_OUT << "\tCalculate correlations in all neighbourhoods\n";
    //output: neighbourhoodCorrs now stores all correlations across neighbourhoods (outer vector is feature-pairs for quick L calculation)
    step_2_calculate_correlation(correlationStrengthCutoff, threads);

    LOCO_OUT << "STEP 3";
    LOCO_OUT << "\tCalculate Laplacian score\n";
    step_3_calculate_laplacian_score(threads);

    if(calcSets)
    {
        LOCO_OUT << "STEP 4";
        LOCO_OUT << "\tCalcualte correlated feature sets (biological programs)";
        std::vector<std::vector<int>> featureSetsRaw;
        step_4_calculate_feature_sets(minFeatureSetSize);
    }
}

void Neighborhood::create_neighborhood_graph(int knn)
{
    //create graph of neighborhoods: input is only the central nodes of neighborhoods
    std:;cout << "Create neighbourhood graph\n";
    std::shared_ptr<GraphData> scGraphData = std::make_shared<GraphData>(centralNeighborhoodPtrs, cellStateGenes, knn, &GraphIni::cell_similarity_graph_manhattan_nodes);
    //inout radius does not matter, but set bandwidth to neg. value to skip estimation
    neighborhoodGraph = std::make_shared<GraphHandler>(scGraphData, knn, 0, -1);
    neighborhoodGraph->create_graph();
}

Neighborhood::Neighborhood(const std::shared_ptr<const GraphData> scData, unsigned int neighborhoodNumber, 
                           unsigned int neighborhoodSize, int neighborhoodKNN,
                           const SingleCellData& inputData,
                           const std::vector<int>& cellStateGenes, const std::vector<int>& corrStateGenes, int permutations,
                           const double& corrSetAbundance, const unsigned int correlatedSetMode, const std::string& correlationType) : 
                           neighborhoodSize(neighborhoodSize), neighbourhoodNum(neighborhoodNumber), inputDataOrigional(inputData),
                           cellStateGenes(cellStateGenes), corrStateGenes(corrStateGenes), permutations(permutations),
                           minimumCorrSetAbundance(corrSetAbundance), correlatedSetMode(correlatedSetMode), correlationType(correlationType)
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

//the final R structure needs:
//return R-object of LoCo results
// create a List of several dataframes
// 1.a) raw data table
// 1.b) UMAP data table

// 2.a) neighbourhood - cells
// 2.b) laplacian scores
// 2c) shuffled results

// 3a) neighbourhood - coords
// 3b) neighbourhood correlations data table
void Neighborhood::fill_result_data(

    std::vector<std::string>& nIDs, // all neighborhoods IDs
    std::vector<std::string>& nID_anchorCellID, //achnor cell IDs for neighborhoods
    std::vector<std::vector<std::string>>& nID_allCellIDs, //vector off all cellIDs for all neighborhoods (same order as nIDs)

    std::vector<std::string>& correlation_pairs_sorted, //all names of the correlation pairs

    std::vector<std::string>& correlation_pairs_origional, //all names of the correlation pairs
    std::vector<std::vector<double>>& corrMat, //all correlations

    std::vector<double>& corrL, 
    std::vector<double>& pCorrL, 
    std::vector<std::string>& featureSetString,
    const bool calcFeatureSets
)
{

    // =========================
    // anchor cell IDs for neighbourhood IDs
    // =========================
    const std::vector<nodePtr> nodes = neighborhoodGraph->get_all_nodes();
    for(const auto& node : nodes)
    {
        std::string cellID = node->get_name();
        nIDs.push_back("N_" + cellID);
        nID_anchorCellID.push_back(cellID);
    }

    // =========================
    // LAPLACIAN RESULTS
    // =========================
    size_t numPairs = neighbourhoodCorrs.pairNames.size();
    std::vector<size_t> indices(numPairs);
    std::iota(indices.begin(), indices.end(), 0);
    // Sort the indices based on laplacianScores.L in increasing order
    std::sort(indices.begin(), indices.end(), 
        [this](size_t a, size_t b) {
            return laplacianScores.L[a] < laplacianScores.L[b];
        }
    );
    // Use the sorted indices to populate your lists
    for(size_t i = 0; i < numPairs; ++i)
    {
        size_t idx = indices[i]; // Grab the index that belongs in this sorted position

        //store correlations in order for Laplacian
        std::string featurePairTmp = neighbourhoodCorrs.pairNames[idx];
        correlation_pairs_sorted.push_back(featurePairTmp);

        // values
        corrL.push_back(laplacianScores.L[idx]);
        pCorrL.push_back(laplacianScores.p_values[idx]);

        // get feature sets
        if(calcFeatureSets)
        {
            featureSetString.push_back(neighbourhoodCorrs.featureSetString[idx]);
        }
    }

    // =========================
    // CORRELATION MATRIX
    // =========================
    size_t numNeighborhoods = nIDs.size();

    //store correlation in order of the correlation matrix
    correlation_pairs_origional = neighbourhoodCorrs.pairNames;

    //check order is as expected
    for (size_t i = 0; i < numNeighborhoods; ++i)
    {
        std::string graphName = "N_" + neighborhoodGraph->get_node_at(i)->get_name();
        
        if (nIDs[i] != graphName)
        {
            throw std::runtime_error(
                "Order assumption failed! Matrix column " + std::to_string(i) + 
                " is '" + graphName + "', but nIDs expects '" + nIDs[i] + "'."
            );
        }
    }

    // Pre-allocate target matrix to exact final size
    corrMat.resize(numNeighborhoods, std::vector<double>(numPairs, 0.0));

    // Do the direct memory transpose
    for (size_t p = 0; p < numPairs; ++p)
    {
        // Grab the flat row for this feature pair once
        const double* pairRow = neighbourhoodCorrs.pairToListOfAllNCorrs.get_row(p);

        for (size_t n = 0; n < numNeighborhoods; ++n)
        {
            // Because we verified the order above, we know that column 'n' 
            // perfectly matches the target row 'n'.
            corrMat[n][p] = pairRow[n]; 
        }
    }

    // =========================
    // CELL IDs PER NEIGHBORHOOD
    // =========================
    std::vector<std::string> nIDs_cell_order;
    for(const auto& neighbourhoodCellList : neighborhoods)
    {
        std::vector<std::string> cellIDList;
        //every cellID is an int from the vector of ints for cell-positions in raw data
        for(const auto& cellID : neighbourhoodCellList.second)
        {
            //get the cellID at position cell and push into cellIDList 
            cellIDList.push_back(inputDataOrigional.cellIDs.at(cellID));
        }
        nID_allCellIDs.push_back(cellIDList);
        nIDs_cell_order.push_back("N_" + neighbourhoodCellList.first->get_name());
    }

    //reorder according to nIDs:
    // Build lookup: current order -> index
    std::unordered_map<std::string, size_t> indexMap_allCellNids;
    for (size_t i = 0; i < nIDs_cell_order.size(); ++i)
    {
        indexMap_allCellNids[nIDs_cell_order[i]] = i;
    }
    // Prepare reordered containers
    std::vector<std::vector<std::string>> reorderedCellIDs;
    reorderedCellIDs.reserve(nIDs.size());
    // Reorder according to nIDs
    for (const auto& nID : nIDs)
    {
        auto it = indexMap_allCellNids.find(nID);
        if (it == indexMap_allCellNids.end())
        {
            throw std::runtime_error("nID not found in cellID data: " + nID);
        }
        size_t idx = it->second;
        reorderedCellIDs.push_back(nID_allCellIDs[idx]);
    }
    // Replace originals
    nID_allCellIDs = std::move(reorderedCellIDs);
}
