#include <iostream>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <unordered_set>
#include <random>

#include "Neighborhood.hpp"
#include "loco_io.h"

using namespace boost::program_options;

/*
output: path/file.tsv (program adds content specific string between file and .tsv)

- file of Corr/Slope Laplacians per Pair of features with their origional cliaues this pair was in
- 2 files of matrices for correlations/ slopes in every neighborhood 
- file of neighborhood coordinates

- corr/ cellstate gene files: files of comma seperated genes (one line). These genes are used for ether ONLY neighborhood
    building or ONLY the correlations between genes. 
- per default MAX top 200 corr and slpe genes (total 400) are written out

*/
bool parse_arguments(char** argv, int argc, std::string& inFile,  std::string& outFile, std::string& prefix,
                     int& threats, char& del, bool& col, bool& row,
                     unsigned int & numNeighborhoods, std::string& neighborhoodSizeStr, int& neighborhoodKNN, double& correlationCutoff,
                     std::string& cellStateGeneFile, std::string& correlationStateGeneFile,
                     bool& zscore, int& permutations, int&minSetSize, double& corrSetAbundance, unsigned int& correlatedSetMode,
                     std::string& correlationType, bool& calcSets, double& minVaryingFraction, std::string& nSampling,
                     bool& considerOnlyNonZeroCells, double& minFractionNonZeroCells)
{
    try
    {
        options_description desc("Options");
        desc.add_options()
            //INPUT/ OUTPUT
            ("input,i", value<std::string>(&inFile)->required(), "Input File in tsv format as a Cell * Gene matrix")
            ("output,o", value<std::string>(&outFile)->default_value("bin"), "Output Directory")
            ("prefix,p", value<std::string>(&prefix)->default_value("LoCo"), "Prefix for output files, e.g., names of the analysis.")

            //INPUT FORMAT
            ("delimiter,d", value<char>(&del)->default_value('\t'), "delimiter of input file")
            ("column,c", "explicitly parse column names (otherwise parsed only if string)")
            ("row,r", "explicitly parse row names (otherwise parsed only if string)")

            //NEIGHBORHOOD VARIABLES
            ("numNeighborhoods,n", value<unsigned int>(&numNeighborhoods)->default_value(0), "number of neighborhoods. By default this is the total number of cells divided by 50 \
                (the default number of cells per neighborhood). You can easily choose more neighborhoods/ cells per neighborhood but be aware that by doing so \
                neighborhoods will start to overlap, which artificially smoothes correlations between neighborhoods leading to an overestimation of p-values.")
            ("neighborhoodSize,s", value<std::string>(&neighborhoodSizeStr)->default_value("100"),
                "number of cells per neighborhood. Can be a single value (e.g. 100) or a range [size1, size2,...] to detect correlation patterns at different scales/ with different granularities.")
            ("correlationCutoff,x", value<double>(&correlationCutoff)->default_value(0.5), "threshold for minimum correlation strength between features")
            ("StateSpaceGenes,v", value<std::string>(&cellStateGeneFile)->default_value(""), "File with list of genes for state space (defines neighborhoods)")
            ("CorrSpaceGenes,w", value<std::string>(&correlationStateGeneFile)->default_value(""), "File with list of genes for correlations space (defines correlations that change through state space)")
            ("NeighborhoodKNN,y", value<int>(&neighborhoodKNN)->default_value(5),"KNN for the neighborhood graph: The number of nearest neighbors that get connected to every neighborhood")
            ("permutations,u", value<int>(&permutations)->default_value(500),"number of permutations to calculate p-values for laplacian scores. \
            This is the number of random assignments of correlations/ slopes to other neighborhoods")
            ("minimumCorrelationSetSize,m", value<int>(&minSetSize)->default_value(2), "The minimum size for sets of correlated features.\
            Loco reports the laplacian score for pairs of features. However, those features are first filtered by finding sets of correlated features.\
            These sets form a small subgraph of features that seem to interact with each other. If this parameter is set higher than 2 features are first filtered\
            by features that are somewhere in the single-cell space part of a bigger correlation set of interacting features, and only for those LoCo \
            calculates pairwise Laplacian scores.")
            ("correlationSetAbundance,a", value<double>(&corrSetAbundance)->default_value(0.01, "0.01"), "Abundance of correlations sets. This is a filter\
            for the set of correlated features to only retain correlated sets that are observed in more than one neighborhood. For all correlation pairs from correlated sets\
            Loco checks if this pair spans several connected neighborhoods and only retains pairs that do. This is the fraction of connected neighborhoods in which a correlation between pairs must be present. Imagine we find a correlation set of features A,B,C,D.\
            Loco then does a BFS from every node and tests which pairs of this set span a region in the Neighborhoodgraph that covers at least 5% of neighborhoods. \
            Imagine correlations between A,B,C are present in 10% but any correlation with D only in a single neighborhood: LoCo will use only the set A,B,C.")

            ("correlatedSetMode,q", value<unsigned int>(&correlatedSetMode)->default_value(1), "Mode to detect correlated sets. LoCo creates a graph of features\
            where features represent nodes and feature-nodes are connected if their correlation between them is >= correlationCutoff(-x flag). LoCo creates these \
            graphs for every neighbourhood, filters subgraphs based on different graph-connectivity and finally only retains features of these graphs for local \
            correlation detection. Graphs can be detected as follows: 0 = fully connected components (cliques), 1 = connected components (nodes with at least one edge), etc.")
            ("CorrelationType,b", value<std::string>(&correlationType)->default_value("spearman"), "spearman or pearson for the type of correlation to calculate")
            ("minVaryingFraction,g", value<double>(&minVaryingFraction)->default_value(0.05, "0.05"), "minimum fraction of cells (0-1) in a neighbourhood whose value \
            differs from the most frequent value of a feature (e.g. cells expressing a gene when most cells have 0 counts). In neighbourhoods where a feature is below it, \
            all correlations of this feature are NA (not calculated). This removes spurious correlations driven by single cells (e.g. all zeroes and a single 1). 0 disables the filter.")

            ("considerOnlyNonZeroCells,e", value<bool>(&considerOnlyNonZeroCells)->default_value(0), "1 = the correlation of a pair \
            in a neighbourhood uses only the cells in which BOTH features are non-zero. ONLY EXACT ZEROS (value == 0 in the input) are excluded, \
            so only real zero values are left out; with the z-score (-z 1, the default) there are no exact zeros: use -z 0. \
            NA if fewer than 3 such cells. Default 0 (all cells).")
            ("minFractionNonZeroCells,j", value<double>(&minFractionNonZeroCells)->default_value(0.1, "0.1"), "minimum fraction of cells (0-1) \
            in a neighbourhood in which BOTH features of a pair are non-zero (ONLY EXACT ZEROS count as zero; no effect with -z 1, use -z 0); \
            below it the correlation of the pair in this neighbourhood is NA. 0 disables the filter.")

            ("Nsampling,k", value<std::string>(&nSampling)->default_value("random"), "how neighbourhoods are built: random (default; \
            random anchor cells plus their nearest cells; neighbourhoods can overlap) or kmeans (non-overlapping neighbourhoods by a \
            size-constrained k-medians clustering)")

            //GENERAL
            ("FeatureSets,f", value<bool>(&calcSets)->default_value(0), "calculate feature sets: sets of features that co-correlate,default is false (0) (1=true, 0=false)")
            ("zScore,z", value<bool>(&zscore)->default_value(1), "z-score normalize data, default is true (1) (1=true, 0=false)")
            ("thread,t", value<int>(&threats)->default_value(5), "number of threads")
            ("help,h", "help message");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if(vm.count("help"))
        {
            std::cout << desc << "\n";
            std::cout << "EXAMPLE CALL:\n ./bin/loco -i <inFile>\n";
            return false;
        }

        notify(vm);

        //parse col names
        if (vm.count("column")){col = true;}
        //parse row names
        if (vm.count("row")){row = true;}
        //zscore: set by notify() from -z (0/1); vm.count("zScore") is always 1 because of the default value, so it
        //must not be used to switch the z-score on (that made -z 0 ineffective)
    }
    catch(std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return false;
    }
    return true;
}

std::vector<unsigned int> parseNeighborhoodSizes(const std::string& input)
{
    std::vector<unsigned int> result;

    std::string s = input;
    s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());

    if (s.empty())
        throw std::invalid_argument("Empty input for neighborhood size");

    if (s.front() == '[' && s.back() == ']')
    {
        s = s.substr(1, s.size() - 2);
        std::stringstream ss(s);
        std::string item;

        while (std::getline(ss, item, ','))
        {
            try {
                result.push_back(static_cast<unsigned int>(std::stoul(item)));
            } catch (const std::exception&) {
                throw std::invalid_argument("Invalid number in list: '" + item + "'");
            }
        }

        if (result.empty())
            throw std::invalid_argument("No numbers found inside brackets");
    }
    else
    {
        try {
            result.push_back(static_cast<unsigned int>(std::stoul(s)));
        } catch (const std::exception&) {
            throw std::invalid_argument("Invalid single number: '" + s + "'");
        }
    }

    return result;
}

void run_correlation_propagation_across_graph(const SingleCellData& inFile, const std::string& outFile, std::string& prefix, int thread,
                                              const unsigned int numNeighborhoods, const std::vector<unsigned int>& neighborhoodSizes, 
                                              const int neighborhoodKNN, const double& correlationCutoff,
                                              const std::vector<std::string>& cellStateGenes,
                                              const std::vector<std::string>& corrStateGenes, 
                                              const int permutations, const int minSetSize, const double corrSetAbundance, 
                                              const unsigned int correlatedSetMode, const std::string& correlationType, bool calcSets,
                                              const double minVaryingFraction, const std::string& nSampling,
                                              const bool considerOnlyNonZeroCells, const double minFractionNonZeroCells)
{
    //generate cell-cell neighborhood graph
    std::vector<int> cellStateIdxs = get_indexlist_from_genenames(inFile, cellStateGenes);
    std::vector<int> corrIdxs = get_indexlist_from_genenames(inFile, corrStateGenes);

    for(const unsigned int neighborhoodSize : neighborhoodSizes)
    {
        // create graph of single-cell data
        unsigned int numNeighborhoodsCalculated = numNeighborhoods;
        if(numNeighborhoods == 0)
        {
            numNeighborhoodsCalculated = inFile.pointCloud.size() / neighborhoodSize;
        }
        std::cout << "STEP 1\n";
        std::cout << "\tCreating " << numNeighborhoodsCalculated << " neighbourhoods with " << neighborhoodSize << " cells\n";
        
        if(numNeighborhoodsCalculated == 1)
        {
            std::cout << "WARNING: YOU ARE RUNNING LoCo WITH A SINGLE NEIGHBOURHOOD!!!\n";
            std::cout << " LoCo still estimates correlations within this neighbourhood but can not detect smooth changes\n";
        }
        bool printStatusUpdateCellDistCalc = true;
        unsigned int scGraphKnn = neighborhoodSize; //the KNN value is the number of cells in a neighborhood, we ONLY have to calcualte the knn closest neighbors, no need for more
        bool precalculateAllDistances = false;
        std::shared_ptr<GraphData> scNormData = std::make_shared<GraphData>(inFile, cellStateIdxs, scGraphKnn, &GraphIni::cell_similarity_graph_manhattan_raw, thread, printStatusUpdateCellDistCalc, precalculateAllDistances);

        //create Neighborhoods
        Neighborhood neighborhood(scNormData, numNeighborhoodsCalculated, neighborhoodSize, neighborhoodKNN, 
                                inFile, cellStateIdxs, corrIdxs, permutations, corrSetAbundance,
                                correlatedSetMode, correlationType, nSampling, thread);
        neighborhood.calculate_correlation_propagation(correlationCutoff, minSetSize, calcSets, thread,
                                                       0.0, 0, 0.05, minVaryingFraction,
                                                       considerOnlyNonZeroCells, minFractionNonZeroCells);

        //write results to file:
        //neighborhood, coordinates, correlation, slope for every protein-pair
        //file for protein-pair to origional clique
        std::string tmpPrefix = prefix;
        if(neighborhoodSizes.size() > 1)
        {
            tmpPrefix = prefix + "_Nsize_" + std::to_string(neighborhoodSize);
        }

        neighborhood.write_results_to_file(outFile, tmpPrefix, calcSets);
    }
}


/*
INPUT: cell-protein count matrix

OUTPUT: NEIGHBORHOOD - CORRELATION_VALUE matrix: 
    correlation value is a mixture of (correlation set_ actual correlations that change)
    algorithm finds correlations sets and then subsets within every correlations set that change differently across the cell-cell graph
    e.g.: for set A-B-C-D-E we might have the pair A-B-C that changes smoothly in all directions and the set D-E that changes slightly differently
    1.) correaltion:
    2.) slope:
additional annotation files:
    neighborhood to mean coordinates: mean coordinates for all genes
    correlation set ID to actual proteins that r correlated

*/


int main(int argc, char** argv)
{

    //PARSE ARGUMENTS
    std::string inFile;
    std::string outFile;
    std::string prefix;
    char del;
    bool col = false;
    bool row = false;
    bool zscore = false;
    bool calcSets = false;
    int thread;
    unsigned int correlatedSetMode;

    //gene lists for states/ correlations
    std::string cellStateGeneFile = "";
    std::string correlationStateGeneFile = "";

    //Neighborhood Variables
    unsigned int numNeighborhoods = 0; //number of neighborhoods
    std::vector<unsigned int> neighborhoodSizes; //number of cells in every neighborhood
    std::string neighborhoodSizeString;
    double correlationCutoff;
    int neighborhoodKNN; //number of nodes coonected to every node in the neighborhood graph (default 5)
    int permutations;
    int minSetSize;
    double corrSetAbundance;
    std::string correlationType;
    double minVaryingFraction;
    std::string nSampling;
    bool considerOnlyNonZeroCells = false;
    double minFractionNonZeroCells = 0.1;

    if(!parse_arguments(argv, argc, inFile, outFile, prefix,
                        thread, del, col, row, 
                        numNeighborhoods, neighborhoodSizeString, neighborhoodKNN, correlationCutoff,
                        cellStateGeneFile, correlationStateGeneFile, 
                        zscore, permutations, minSetSize, corrSetAbundance, correlatedSetMode, correlationType, calcSets, minVaryingFraction, nSampling,
                        considerOnlyNonZeroCells, minFractionNonZeroCells))
    {
        exit(EXIT_FAILURE);
    }

    if (!(correlationType == "spearman" || correlationType == "pearson")) 
    {
        std::cout << "Method '" << correlationType << "' is not a valid method for correlations (spearman|pearson)." << std::endl;
        exit(EXIT_FAILURE);
    } 

    if (minVaryingFraction < 0.0 || minVaryingFraction > 1.0)
    {
        std::cout << "minVaryingFraction must be between 0 and 1." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (minFractionNonZeroCells < 0.0 || minFractionNonZeroCells > 1.0)
    {
        std::cout << "minFractionNonZeroCells must be between 0 and 1." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (!(nSampling == "kmeans" || nSampling == "random"))
    {
        std::cout << "Nsampling must be 'kmeans' or 'random'." << std::endl;
        exit(EXIT_FAILURE);
    }

    neighborhoodSizes = parseNeighborhoodSizes(neighborhoodSizeString);

    //READ IN DATA
    SCParser parser(inFile, del, col, row);
    SingleCellData inputDataRaw = parser.getData();

    if(zscore)
    {
        std::cout << "z-score normalize data (per feature: (value - mean) / sd over all cells)\n";
        zscore_singleCelldata(inputDataRaw);
    }

    //Read in gene lists
    std::vector<std::string> cellStateGenes;
    std::vector<std::string> corrStateGenes;

    if(cellStateGeneFile != "")
    {
        cellStateGenes = parse_list(cellStateGeneFile, ',');
    }
    if(correlationStateGeneFile != "")
    {
        corrStateGenes = parse_list(correlationStateGeneFile, ',');
    }
    print_feature_roles(inputDataRaw, cellStateGenes, corrStateGenes);

    //CALL TOOL
    std::cout << "Running Correlation Analysis over single-cell graph\n";
    run_correlation_propagation_across_graph(inputDataRaw, outFile, prefix, thread,
                                              numNeighborhoods, neighborhoodSizes, 
                                              neighborhoodKNN, correlationCutoff,
                                              cellStateGenes, corrStateGenes, 
                                              permutations, minSetSize, corrSetAbundance, correlatedSetMode,
                                              correlationType, calcSets, minVaryingFraction, nSampling,
                                              considerOnlyNonZeroCells, minFractionNonZeroCells);

    return(EXIT_SUCCESS);
}