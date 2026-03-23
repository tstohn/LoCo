#ifndef SC_PARSER
#define SC_PARSER

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <exception>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <cmath>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

//struct storing the single cell data after parsing it from its raw file
struct SingleCellData
{
    //names of cells/ proteins: ONLY USED when COLUMN/ ROW names provided
    //checked by looking for non-numeric first COLUMNS, ROWS
    std::vector<std::string> geneNames;
    std::vector<std::string> cellIDs;
    std::vector<std::string> clusterIDs;

    //single-cells * protein counts
    std::vector<std::vector<double>> pointCloud;

    //map from gene name -> idx
    std::unordered_map<std::string, int> geneNameToIdx;
};
SingleCellData filter_singleCelldata(const SingleCellData& origionalScData, const std::vector<int>& indices);
void zscore_singleCelldata(SingleCellData& data);
std::vector<int> get_indexlist_from_genenames(const SingleCellData& scData, const std::vector<std::string>& geneList);

//generall parsing function (for gene lists)
std::vector<std::string> parse_list(const std::string& listFile, const char& sep);

/**
 * @brief Parser for single cell input data (matrix of CELL * VARIABLES[GENE/PROTEIN])
 * It can parse a cell*protein or the transpose matrix (when trans parameter set - transpose matrix can read ONLY gene*cell matrix)
 * IT CAN NOT READ ADDITIOANL ANNOTATIONS LIKE CLUSTERIDs
 * It can contain a header columns and a first row with row names...
 * Additionally is can contain another annotation variable in clusterID
 * 
 * When the first col, row is a string it assumes they r col, row names.
 * When they are not you can EXPLICITELY SET row, col names with -col/ -row,
 * otherwise the names are simply enumerated from 0
 * 
 * TODO: parse poitnCloudinto a shared PTR to not copy it when returning
 */
class SCParser
{
    public:
        SCParser(std::string tsvFile, const char& del, const bool& trans, 
                 const bool& col = false, const bool& row = false,
                 std::string clusterID = "");
        const std::vector<std::vector<double>> getPointVector() const
        {
            return data.pointCloud;
        }
        const std::vector<std::string> getGeneNames() const 
        {
            return data.geneNames;
        }
        const SingleCellData getData() const
        {
            return(data);
        }

    private:

        template<typename streamType>
        unsigned int get_number_of_lines(streamType& instream, std::string line);
        unsigned int get_number_of_cells(const std::string& line, const char& del, bool& rowNameLine);
        bool has_gene_names(const std::string& line, const char& del);
        void generate_cell_names(const std::string& firstLine, const char del, bool rowNames, bool colNames, int numCells);

        void readGeneNames(const std::string& line, const char& del, 
                           bool& hasColumnNames, int& numDimensions);
        void readValueLine(const std::string& line, const char& del, 
                           bool& hasRowNames, const bool& hasColumnNames,
                           int& numDimensions, const int& lineNum);
        void readValueLineTranspose(const std::string& line, const char& del, 
                                     const int& row, const bool& rowNames);

        //function to parse tsv file: if col, row are set to false, the col, row is still parsed if they are strings
        template<typename streamType>
        void parseTsvFile(const std::string& file, const char& del, const bool& trans, streamType& instream, const bool& col, const bool& row);
        int clusterIDPosition = -1;
        std::string clusterID;

        SingleCellData data;
};

#endif