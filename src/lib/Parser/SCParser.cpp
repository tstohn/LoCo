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
            //cluster data might be empty
            if(!origionalScData.clusterIDs.empty())
            {
                resultData.clusterIDs.push_back(origionalScData.clusterIDs.at(index));
            }
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
                            bool& hasColumnNames, int& numDimensions)
{
    std::stringstream ss;
    ss.str(line);
    std::string substr;
    int count = 0;
    int nextGeneId = 0;
    bool clusterIDFound = false; //keep track if we parse cellIDs
    //iterate over columns of the first row
    while(getline( ss, substr, del ))
    {
        //check only first Column if its a string even through no column names flag was given - in that case DO PARSE column names
        if(count == 0 && !is_number(substr) && (substr != clusterID))
        {
            //set col names to true, if there is NO NUMBER in first column elements
            hasColumnNames = true;
        }
        if(!is_number(substr) && count > 0 && !hasColumnNames)
        {
            std::cerr << "Encountered a string in first row of input data - owever expecting no column names !!!!";
            exit(EXIT_FAILURE);
        }

        //check if column is the cellID column
        if(substr == clusterID)
        {
            //set the position in every row where to find cellID
            clusterIDPosition = count;
            clusterIDFound = true;
        }
        else
        {
            std::string colname = substr;
            if(hasColumnNames == false)
            {
                colname = std::to_string(nextGeneId);
            }
            std::cout << "adding gene: " << colname << "\n";
            data.geneNames.push_back( colname );
            ++nextGeneId;
        }

        ++count;
    }

    if(!clusterIDFound && clusterID != "")
    {
        std::cerr << "Did not find provided column for <clusterID> !!\n";
        exit(EXIT_FAILURE);
    }
    numDimensions = data.geneNames.size();
}

void SCParser::readValueLine(const std::string& line, const char& del, bool& hasRowNames,
                            const bool& hasColumnNames, int& numDimensions, const int& lineNum)
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
            if(!is_number(substr) && (columnNum != clusterIDPosition))
            {
                --numDimensions; //in case we have rownames we have to remove the dimension of the first column header
                hasRowNames = true;
            }
        }

        //setting empty column entries to zero
        if(substr == ""){substr = "0";}

        //if we are in the clusterID column
        if( (clusterIDPosition != -1) && (columnNum == clusterIDPosition))
        {
            data.clusterIDs.push_back(substr);
        }
        else if(hasRowNames && columnNum==0) //if we have to parse cellIDs in first column of every row
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
                                std::string(") do not match!");
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

void SCParser::readValueLineTranspose(const std::string& line, const char& del, 
                                     const int& row, const bool& rowNames)
{

    std::stringstream ss;
    ss.str(line);
    std::string substr;
    int columnNum = 0;
    while(getline( ss, substr, del ))
    {
        if(substr == ""){substr = "0";}
        if(rowNames && columnNum==0)
        {
            //if we have rownames store the current rownames (gene name)
            data.geneNames.push_back(substr);
            ++columnNum;
            continue;
        }
        else if(!rowNames && columnNum==0)
        {
            //insert an ID of current count for gene name
            data.geneNames.push_back(std::to_string(data.geneNames.size()));
        }
        double value = std::stod(substr);
        int vectorIndex = columnNum;
        if(rowNames){--vectorIndex;}
        data.pointCloud.at(vectorIndex).push_back(value);

        ++columnNum;
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
void SCParser::parseTsvFile(const std::string& file, const char& del, const bool& trans, streamType& instream,
                            const bool& col, const bool& row)
{
    std::string line;
    //reading lines
    if(!trans)
    {
        int lineNum = 0;
        int numDimensions = 1;
        bool hasRowNames = row; //true: means PARSE, false: means do ONLY PARSE IF ITS A STRING
        bool hasColumnNames = col; //true: means PARSE, false: means do ONLY PARSE IF ITS A STRING
        while(std::getline(instream, line))
        {

            //if lineNum=0 read the header of the matrix if we have column names
            if(lineNum == 0)
            {
                readGeneNames(line, del, hasColumnNames, numDimensions);
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
    else
    {
        //the transpose is slightly different since we need to count cells first to insert values at the right position
        int lineNum = 0;
        int row = 0;
        bool hasRowNames = row; //rownames are now the gene names
        bool hasColNames = col; //colnames are now cellIDs
        while(std::getline(instream, line))
        {
            if(lineNum == 0)
            {
                bool cellNameLine = false;
                //get number of cells, in stpore in cellNameLine if this row contains only strings (and therefore cell names)
                unsigned int numCells = get_number_of_cells(line, del, cellNameLine);
                std::string firstLine = line;
                //basically check if we do have row names (strings first col, even though the row name flag was not explicitely set)
                if(!cellNameLine && !hasRowNames) //if we were not in a column name line check this one
                {
                    //read value line
                    hasRowNames = has_gene_names(line, del);
                }
                else if(cellNameLine && !hasRowNames) //otherwise check next one
                {
                    //if it was a cellNameLine, get the next line to get the real number of cells
                    std::getline(instream, line);
                    hasRowNames = has_gene_names(line, del);
                    ++lineNum;
                }
                //substract one for cell number if first column was actually the row names
                if(hasRowNames){--numCells;}
                //generate the cellID (cell names) - this depends on if we have rowNames, colnames
                

                //resize temporary vector to right size
                data.pointCloud.resize(numCells, std::vector<double>(0));
                //current string in line contains first line with data, we still have to parse it
                readValueLineTranspose(line, del, row, hasRowNames);

                //in readLine: if hasCOl add all strings to CellIDs, otherwise add number
                //row and colanmes is still mixed: switch it
                //check in get_number_cell: check if we have colNames
                //then in readvalue if in first line add colnames otherf=wise if rownames add to rownames...

                ++lineNum;
                ++row;
            }
            else
            {
                readValueLineTranspose(line, del, row, hasRowNames);
                ++row;
            }
        }
    }
}

SCParser::SCParser(std::string tsvFile, const char& del, const bool& trans, 
                   const bool& col, const bool& row,
                   std::string clusterID) : clusterID{clusterID}
{
    //read in tsvFile: create firstly only a pointCloud
    std::cout << "Reading Input file: delimiter=<" << del << ">, transpose=<" << trans << ">\n";

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
        parseTsvFile<std::istream>(tsvFile, del, trans, istream, col, row);

    }
    else
    {
        std::ifstream instream(tsvFile);
        parseTsvFile<std::ifstream>(tsvFile, del, trans, instream, col, row);
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

template void SCParser::parseTsvFile<std::ifstream>(const std::string& file, const char& del, const bool& trans, 
                                                   std::ifstream& instream, const bool& col, const bool& row);
template void SCParser::parseTsvFile<std::istream>(const std::string& file, const char& del, const bool& trans, 
                                                  std::istream& instream, const bool& col, const bool& row);
