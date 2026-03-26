#include "SCParser.hpp"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

std::vector<int> get_indexlist_from_genenames(const SingleCellData& scData, const std::vector<std::string>& geneList)
{
    std::vector<int> geneIdxs;
    for(const std::string& gene : geneList)
    {
        int idx;
        try {
            idx = scData.geneNameToIdx.at(gene);
        }
        catch (const std::out_of_range& e) {
            std::cerr << "Gene name not found: " << gene << "\n";
            std::cerr << "Please double check your file for correlation/ state markers\n";
            exit(EXIT_FAILURE);
        }

        geneIdxs.push_back(idx);
    }
    return(geneIdxs);
}


std::vector<std::string> parse_list(const std::string& listFile, const char& sep)
{
    std::ifstream instream(listFile);
    std::vector<std::string> list;
    std::string line;

    if (instream.is_open()) 
    {
        while (std::getline(instream, line)) 
        { 
            std::stringstream ss(line);
            std::string substr;
            while(getline(ss, substr, sep))
            {
                list.push_back(substr);    
            }
        }
        instream.close();
    } 
    else 
    {
        std::cerr << "Unable to open file!" << std::endl;
        exit(EXIT_FAILURE);
    }

    return(list);
}

void zscore_singleCelldata(SingleCellData& data)
{
    if (data.pointCloud.empty())
        return;

    const size_t numCells = data.pointCloud.size();
    const size_t numFeatures = data.pointCloud[0].size();

    // ensure cells have same number of features
    for (const std::vector<double>& row : data.pointCloud)
    {
        if (row.size() != numFeatures)
            throw std::runtime_error("pointCloud is not rectangular. Some cells have more features than others.\n");
    }

    // For each feature (column)
    for (size_t j = 0; j < numFeatures; ++j)
    {
        double mean = 0.0;
        double sq_sum = 0.0;

        // Compute mean
        for (size_t i = 0; i < numCells; ++i)
            mean += data.pointCloud[i][j];

        mean /= static_cast<double>(numCells);

        // Compute variance
        for (size_t i = 0; i < numCells; ++i)
        {
            double diff = data.pointCloud[i][j] - mean;
            sq_sum += (diff * diff);
        }

        double variance = sq_sum / static_cast<double>(numCells);

        double stddev = std::sqrt(variance);

        // Avoid division by zero
        if (stddev == 0.0)
        {
            for (size_t i = 0; i < numCells; ++i)
                data.pointCloud[i][j] = 0.0;
        }
        else
        {
            for (size_t i = 0; i < numCells; ++i)
                data.pointCloud[i][j] = (data.pointCloud[i][j] - mean) / stddev;
        }
    }
}

//filters a singleCellData object with a vector of filter indices and returns the result
SingleCellData filter_singleCelldata(const SingleCellData& origionalScData, const std::vector<int>& indices) 
{
    SingleCellData resultData;
    int scDataEntries = origionalScData.cellIDs.size();

    for (int index : indices) 
    {
        // Check if the index is within bounds
        if (index < scDataEntries) 
        {
            // Add the corresponding element to the result vector
            resultData.cellIDs.push_back(origionalScData.cellIDs.at(index));
            resultData.pointCloud.push_back(origionalScData.pointCloud.at(index));
            //gene names stay the same
            resultData.geneNames = origionalScData.geneNames;
        } 
        else 
        {
            // Handle out-of-bounds index if needed
            std::cerr << "Warning: Index " << index << " is out of bounds." << std::endl;
        }
    }

    return(resultData);
}

namespace
{
inline bool endWith(std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    if((*it)=='-' || (*it)=='+'){++it;} //first sign can be + or -
    while (it != s.end() && (std::isdigit(*it) || (*it)=='.')) ++it;
    return !s.empty() && it == s.end();
}

} // end anonymous namespace

void SCParser::readGeneNames(const std::string& line, const char& del, 
                            bool& hasColumnNames, bool& hasRowNames, size_t& numDimensions)
{
    std::stringstream ss;
    ss.str(line);
    std::string substr;
    int count = 0;
    int nextGeneId = 0;
    //iterate over columns of the first row
    while(getline( ss, substr, del ))
    {
        //if we know that the data MUST have rownmes, we can skip the first entry as it must be emopty/ be sth. like CELL_ID_COLUMN
        if(hasRowNames == true)
        {
            continue;
        }

        //check only first Column if its a string even through no column names flag was given - in that case DO PARSE column names
        if(count == 0 && !is_number(substr))
        {
            //set col names to true, if there is NO NUMBER in first column elements
            hasColumnNames = true;
        }

        //if the first col is empty we assume that the data has rownames (cell IDs)
        if(count == 0 && substr=="")
        {
            hasRowNames = true;
            std::cout << "The first column entry in the first row is empty, we therefore assume that rownames (cell-names) are provided!\
            If this is not the case please remove the empty entry, give a proper gene name or remove the gene-name header entirely\n.";
            continue;
        }

        if(!is_number(substr) && count > 0 && !hasColumnNames)
        {
            std::cerr << "Encountered a string in a later column of the first row of input data - but loco expected NO column names.\n\
            Column names flag was not set and the first values in the row are numeric !!!!";
            exit(EXIT_FAILURE);
        }

        std::string colname = substr;
        if(hasColumnNames == false)
        {
            colname = std::to_string(nextGeneId);
        }
        std::cout << "adding gene: " << colname << "\n";
        data.geneNames.push_back( colname );
        ++nextGeneId;

        ++count;
    }

    numDimensions = data.geneNames.size();
}

void SCParser::readValueLine(const std::string& line, const char& del, bool& hasRowNames,
                            const bool& hasColumnNames, size_t& numDimensions, const size_t& lineNum)
{

    std::stringstream ss;
    ss.str(line);
    std::string substr;
    int columnNum = 0;
    std::vector<double> tmpDimensions;
    while(getline( ss, substr, del ))
    {
        //check for rownames is done in the second line if we have column names
        //if we have no column names its performed in the first line
        if( ((lineNum == 1 && hasColumnNames) || 
            (lineNum == 0 && !hasColumnNames)) && 
            columnNum == 0)
        {
            //check only in the first column of the first true row (only ONCE - to make sure hasRowNames is not suddenly overwritten later on)
            if(!is_number(substr))
            {
                //check for row names only in the first line (value == 0)
                // if the first value in the header (colnames row) was not empty we just realize now that we have rownames
                // in that case we need to also remove the first gene name
                if(lineNum == 0)
                {
                    data.geneNames.erase(data.geneNames.begin());
                    --numDimensions; //in case we have rownames we have to remove the dimension of the first column header
                    hasRowNames = true;
                }
                else if(hasRowNames == false)
                {
                    std::cerr << "We were expecting no row names (cell names). Previous lines did not seem to have cell names.\
                    Please provide input with names (strings) in every row of the first column or in none.\n";
                    exit(EXIT_FAILURE);
                }
            }
        }

        //setting empty column entries to zero
        if(substr == ""){substr = "0";}

        if(hasRowNames && columnNum==0) //if we have to parse cellIDs in first column of every row
        {
            data.cellIDs.push_back(substr);
        }
        else //otherwise parse the actual VALUE
        {
            tmpDimensions.emplace_back( std::stod(substr) );
        }
        ++columnNum;
    }
    if(!hasColumnNames && (lineNum == 0))
    {
        numDimensions = tmpDimensions.size();
    }
    if(tmpDimensions.size() != numDimensions)
    {
        std::string errorMessage = std::string("Dimensions of Header (") + std::to_string(numDimensions) +
                                std::string(") and line ") + std::to_string(lineNum) + 
                                std::string(" (") + std::to_string(tmpDimensions.size()) +
                                std::string(") do not match!\n This might be because if you provide column and row names (cell and feature names)\
                                    the very first column in the first row must contain some placeholder, an empty cell or a value like CELL_IDS.\n");
        throw std::runtime_error(errorMessage);
    }
    data.pointCloud.push_back(tmpDimensions);

    //if we had no rownames, we must add numbers for cellIDs
    if(!hasRowNames)
    {
        int currentCellId = lineNum;
        //in case of columnnames first line has names and true cellID is -1
        if(hasColumnNames)
        {
            currentCellId = lineNum - 1;
        }
        data.cellIDs.push_back(std::to_string(currentCellId));
    }
}

template<typename streamType>
unsigned int SCParser::get_number_of_lines(streamType& instream, std::string line)
{
    unsigned int numberOfLines = 0;
    while (std::getline(instream, line))
    {
        ++numberOfLines;
    }
    instream.clear();
    instream.seekg(0);
    return(numberOfLines);
}

//gets cell number
unsigned int SCParser::get_number_of_cells(const std::string& line, const char& del, bool& rowNameLine)
{
    std::stringstream ss;
    ss.str(line);
    std::string substr;
    int count = 0;
    while(getline( ss, substr, del ))
    {
        //if any column is not a number or is empty, assume that whole row only contains cellID-names
        if( !is_number(substr) || substr.empty() )
        {
            rowNameLine = true;
        }
        ++count;
    }
    return(count);
}

bool SCParser::has_gene_names(const std::string& line, const char& del)
{
    std::stringstream ss;
    ss.str(line);
    std::string substr;
    bool has_gene_names = false;
    getline( ss, substr, del );
    //in the first column check if we have a string or not
    if( !is_number(substr) )
    {
        has_gene_names = true;
    }
    return(has_gene_names);
}

//reads cell names (reads string in case we have cell names, or enumerates cellIDs from 0 to cellnum)
void SCParser::generate_cell_names(const std::string& firstLine, const char del, bool rowNames, bool colNames, int numCells)
{
    //parse actual cellIDs
    if(colNames)
    {
        std::stringstream ss;
        ss.str(firstLine);
        std::string substr;
        int count = 0;
        while(getline( ss, substr, del ))
        {
            if(rowNames && count >0)
            {
                data.cellIDs.push_back(substr);
            }
            else if(!rowNames)
            {
                data.cellIDs.push_back(substr);
            }
            ++count;
        }
    }
    else //simply enumerate
    {
        for(int i =0; i < numCells; ++i)
        {
            data.cellIDs.push_back(std::to_string(i));
        }
    }
}

template<typename streamType>
void SCParser::parseTsvFile(const char& del, streamType& instream,
                            const bool& col, const bool& row)
{
    std::string line;
    //reading lines
    
    size_t lineNum = 0;
    size_t numDimensions = 1;
    bool hasRowNames = row; //true: means PARSE, false: means do ONLY PARSE IF ITS A STRING
    bool hasColumnNames = col; //true: means PARSE, false: means do ONLY PARSE IF ITS A STRING
    while(std::getline(instream, line))
    {

        //if lineNum=0 read the header of the matrix if we have column names
        if(lineNum == 0)
        {
            readGeneNames(line, del, hasColumnNames,hasRowNames, numDimensions);
            if(!hasColumnNames)
            {
                //read value line
                readValueLine(line, del, hasRowNames, hasColumnNames, numDimensions, lineNum);
            }
            ++lineNum;
        }
        else
        {
            //read value line
            readValueLine(line, del, hasRowNames, hasColumnNames, numDimensions, lineNum);
            ++lineNum;
        }        
    }
    
}

SCParser::SCParser(std::string tsvFile, const char& del,
                   const bool& col, const bool& row)
{
    //read in tsvFile: create firstly only a pointCloud
    std::cout << "Reading Input file: delimiter=<" << del << ">, cell-ids=<" << row << ">, feature-ids=<" << col << ">\n";

    struct stat info;
    if( stat(tsvFile.c_str(), &info ) != 0 )
    {
        std::cout << "File " << tsvFile << " does not exist!\n";
        exit(EXIT_FAILURE);
    }

    if(endWith(tsvFile,".gz"))
    {
        std::ifstream instream(tsvFile, std::ios_base::in | std::ios_base::binary);
        boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
        inbuf.push(boost::iostreams::gzip_decompressor());
        inbuf.push(instream);

        std::istream istream(&inbuf);
        parseTsvFile<std::istream>(del, istream, col, row);

    }
    else
    {
        std::ifstream instream(tsvFile);
        parseTsvFile<std::ifstream>(del, instream, col, row);
    }

    //assert correct dimensionalities
    assert(data.cellIDs.size() == data.pointCloud.size());
    assert(data.geneNames.size() == data.pointCloud.at(0).size());

    //set the map from gene names -> idx
    int idx = 0;
    for(const std::string& name : data.geneNames)
    {
        data.geneNameToIdx.insert(std::pair<std::string, int>(name, idx));
        ++idx;
    }
}

template unsigned int SCParser::get_number_of_lines<std::ifstream>(std::ifstream& instream, std::string line);
template unsigned int SCParser::get_number_of_lines<std::istream>(std::istream& instream, std::string line);

template void SCParser::parseTsvFile<std::ifstream>(const char& del,
                                                   std::ifstream& instream, const bool& col, const bool& row);
template void SCParser::parseTsvFile<std::istream>(const char& del,
                                                  std::istream& instream, const bool& col, const bool& row);
