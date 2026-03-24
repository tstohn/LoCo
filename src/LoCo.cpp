#include <iostream>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <unordered_set>
#include <random>

#include "SCGraph/Constructors/Neighborhood.hpp"

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
                     int& threats, char& del, bool& trans, bool& col, bool& row,
                     unsigned int & numNeighborhoods, std::string& neighborhoodSizeStr, int& neighborhoodKNN, double& correlationCutoff,
                     int& numberCorrelations,
                     std::string& cellStateGeneFile, std::string& correlationStateGeneFile,
                     bool& zscore, int& permutations, int&minSetSize, double& corrSetAbundance, int& correlatedSetMode)
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
            ("column,c", "explicitely parse column names (otherwise parsed only if string)")
            ("row,r", "explicitely parse column names (otherwise parsed only if string)")
            ("transpose,trans", "matrix is transpose (default false)")
            ("filterCorrelations,f", value<int>(&numberCorrelations)->default_value(0), "filter the number of correlations to retain. Only write the gene-pairs of lowest laplacian for correlation/ slope. \
            Since we write the x lowest values for correlation & for slope the total number can be >x.")

            //NEIGHBORHOOD VARIABLES
            ("numNeighborhoods,n", value<unsigned int>(&numNeighborhoods)->default_value(0), "number of neighborhoods. By default this is the total number of cells divided by 50 \
                (the default number of cells per neighborhood). You can easily choose more neighborhoods/ cells per neighborhood but be aware that by doing so \
                neighborhoods will start to overlap, which artificially smoothes correlations between neighborhoods leading to an overestimation of p-values.")
            ("neighborhoodSize,s", value<std::string>(&neighborhoodSizeStr)->default_value("50"),
                "number of cells per neighborhood. Can be a single value (e.g. 50) or a range [size1, size2,...] to detect correlation patterns at different scales/ with different granularities.")
            ("correlationCutoff,x", value<double>(&correlationCutoff)->default_value(0.7), "threshold for minimum correlation strength between features")
            ("StateSpaceGenes,v", value<std::string>(&cellStateGeneFile)->default_value(""), "File with list of genes for state space (defines neighborhoods)")
            ("CorrSpaceGenes,w", value<std::string>(&correlationStateGeneFile)->default_value(""), "File with list of genes for correlations space (defines correlations that change through state space)")
            ("NeighborhoodKNN,y", value<int>(&neighborhoodKNN)->default_value(5),"KNN for the neighborhood graph: The number of nearest neighbors that get connected to every neighborhood")
            ("permutations,u", value<int>(&permutations)->default_value(100),"number of permutations to calcualte p-values for laplacian scores. \
            This is the number of random assignments of correlations/ slopes to other neighborhoods")
            ("minimumCorrelationSetSize,m", value<int>(&minSetSize)->default_value(2), "The minimum size for sets of correlated features.\
            Loco reports the laplacian score for pairs of features. However, those features are first filtered by finding sets of correlated features.\
            These sets form a small subgraph of features that seem to interact with each other. If this parameter is set higher than 2 features are first filtered\
            by features that are somewhere in the single-cell space part of a bigger correltion set of interacting features, and only for those Loco \
            calculates pairwise Laplcaians.")
            ("correlationSetAbundance,a", value<double>(&corrSetAbundance)->default_value(0.01), "Abundance of correlations sets. This is a filter\
            for the set of correlated features to only retain correlated sets that are observed in more than one neighborhood. For all correlation pairs from correlated sets\
            Loco checks if this pair spans several connected neighborhoods and only retains pairs that do. This is the fraction of connected neighborhoods in which a correlation between pairs must be present. Imagine we find a correlation set of features A,B,C,D.\
            Loco then does a BFS from every node and tests which pairs of this set span a region in the Neighborhoodgraph that covers at least 5% of neighborhoods. \
            Imagine correlations between A,B,C are present in 10% but any correlation with D only in a single neighborhood: Lco will use only the set A,B,C.")

            ("correlatedSetMode,q", value<int>(&correlatedSetMode)->default_value(2), "mode to detect correlated set: 0=connected components\
            1=fully connected components, 2=components with>=2edges per node")


            //GENERAL
            ("zScore,z", "z-score normalize the data. Set flag with -z if you want to ahve the data z-scored, (default isx  false)")
            ("thread,t", value<int>(&threats)->default_value(5), "number of threads")
            ("help,h", "help message");

        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);

        if(vm.count("help"))
        {
            std::cout << desc << "\n";
            std::cout << "EXAMPLE CALL:\n ./bin/scluster -i <inFile>\n";
            return false;
        }

        notify(vm);

        //transpose matrix
        if (vm.count("transpose")){trans = true;}
        //parse col names
        if (vm.count("column")){col = true;}
        //parse row names
        if (vm.count("row")){row = true;}
        //zscore data
        if (vm.count("zScore")){zscore = true;}
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
                                              const std::vector<unsigned int>& neighborhoodSizes, 
                                              const int neighborhoodKNN, const double& correlationCutoff,
                                              int& numberCorrelations, const std::vector<std::string>& cellStateGenes,
                                              const std::vector<std::string>& corrStateGenes, 
                                              const int permutations, const int minSetSize, const double corrSetAbundance, 
                                              const bool correlatedSetMode)
{
    //generate cell-cell neighborhood graph
    std::vector<int> cellStateIdxs = get_indexlist_from_genenames(inFile, cellStateGenes);
    std::vector<int> corrIdxs = get_indexlist_from_genenames(inFile, corrStateGenes);

    for(const unsigned int neighborhoodSize : neighborhoodSizes)
    {
        // create graph of single-cell data
        unsigned int numNeighborhoods = inFile.pointCloud.size() / neighborhoodSize;
        std::cout << "Creating " << numNeighborhoods << " neighbourhoods with " << neighborhoodSize << " cells\n";
        bool printStatusUpdateCellDistCalc = true;
        unsigned int scGraphKnn = neighborhoodSize; //the KNN value is the number of cells in a neighborhood, we ONLY have to calcualte the knn closest neighbors, no need for more
        bool precalculateAllDistances = false;
        std::shared_ptr<GraphData> scNormData = std::make_shared<GraphData>(inFile, cellStateIdxs, scGraphKnn, &GraphIni::cell_similarity_graph_manhattan_raw, thread, printStatusUpdateCellDistCalc, precalculateAllDistances);

        //create Neighborhoods
        Neighborhood neighborhood(scNormData, numNeighborhoods, neighborhoodSize, neighborhoodKNN, 
                                inFile, cellStateIdxs, corrIdxs, permutations, corrSetAbundance,
                                correlatedSetMode);
        bool printFoundCliquesPerNeighborhood = false;
        neighborhood.calculate_correlation_propagation(correlationCutoff, printFoundCliquesPerNeighborhood, minSetSize, thread);

        //write results to file:
        //neighborhood, coordinates, correlation, slope for every protein-pair
        //file for protein-pair to origional clique
        std::string tmpPrefix = prefix;
        if(neighborhoodSizes.size() > 1)
        {
            tmpPrefix = prefix + "_Nsize_" + std::to_string(neighborhoodSize);
        }
        neighborhood.write_results_to_file(outFile, tmpPrefix, numberCorrelations);
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
    bool trans = false;
    bool col = false;
    bool row = false;
    bool zscore = false;
    int thread;
    int correlatedSetMode;

    //gene lists for states/ correlations
    int numberCorrelations;
    std::string cellStateGeneFile = "";
    std::string correlationStateGeneFile = "";

    //Neighborhood Variables
    unsigned int numNeighborhoods; //number of neighborhoods
    std::vector<unsigned int> neighborhoodSizes; //number of cells in every neighborhood
    std::string neighborhoodSizeString;
    double correlationCutoff;
    int neighborhoodKNN; //number of nodes coonected to every node in the neighborhood graph (default 5)
    int permutations;
    int minSetSize;
    double corrSetAbundance;

    if(!parse_arguments(argv, argc, inFile, outFile, prefix,
                        thread, del, trans, col, row, 
                        numNeighborhoods, neighborhoodSizeString, neighborhoodKNN, correlationCutoff,
                        numberCorrelations, cellStateGeneFile, correlationStateGeneFile, 
                        zscore, permutations, minSetSize, corrSetAbundance, correlatedSetMode))
    {
        exit(EXIT_FAILURE);
    }

    neighborhoodSizes = parseNeighborhoodSizes(neighborhoodSizeString);

    //READ IN DATA
    SCParser parser(inFile, del, trans, col, row);
    SingleCellData inputDataRaw = parser.getData();

    if(zscore)
    {
        std::cout << "z-score normalize data (scale feature counts for each single-cell)\n";
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

    //CALL TOOL
    std::cout << "Running Correlation Analysis over single-cell graph\n";
    run_correlation_propagation_across_graph(inputDataRaw, outFile, prefix, thread,
                                              neighborhoodSizes, neighborhoodKNN, correlationCutoff,
                                              numberCorrelations, cellStateGenes, corrStateGenes, 
                                              permutations, minSetSize, corrSetAbundance, correlatedSetMode);

    return(EXIT_SUCCESS);
}