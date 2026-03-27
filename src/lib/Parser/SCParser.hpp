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
#include <zlib.h>
#include <cassert>

class GzStreamBuf : public std::streambuf {
private:
    gzFile file = nullptr;
    std::vector<char> buffer;

public:
    explicit GzStreamBuf(const char* path)
        : buffer(32 * 1024) // 32KB buffer
    {
        file = gzopen(path, "rb");
        if (!file) {
            return; //  is_open() is checkd in Parser
        }

        setg(buffer.data(), buffer.data(), buffer.data());
    }

    ~GzStreamBuf() override {
        if (file) {
            gzclose(file);
        }
    }

    bool is_open() const {
        return file != nullptr;
    }

protected:
    int_type underflow() override {
        if (!file) return traits_type::eof();

        if (gptr() < egptr()) {
            return traits_type::to_int_type(*gptr());
        }

        int bytesRead = gzread(file, buffer.data(), buffer.size());

        if (bytesRead <= 0) {
            return traits_type::eof();
        }

        setg(buffer.data(),
             buffer.data(),
             buffer.data() + bytesRead);

        return traits_type::to_int_type(*gptr());
    }
};

//struct storing the single cell data after parsing it from its raw file
struct SingleCellData
{
    //names of cells/ proteins: ONLY USED when COLUMN/ ROW names provided
    //checked by looking for non-numeric first COLUMNS, ROWS
    std::vector<std::string> geneNames;
    std::vector<std::string> cellIDs;

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
 * It can contain a header columns and a first row with row names...
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
        SCParser(std::string tsvFile, const char& del, 
                 const bool& col = false, const bool& row = false);
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
                           bool& hasColumnNames, bool& hasRowNames, size_t& numDimensions);
        void readValueLine(const std::string& line, const char& del, 
                           bool& hasRowNames, const bool& hasColumnNames,
                           size_t& numDimensions, const size_t& lineNum);

        //function to parse tsv file: if col, row are set to false, the col, row is still parsed if they are strings
        template<typename streamType>
        void parseTsvFile(const char& del, streamType& instream, const bool& col, const bool& row);

        SingleCellData data;
};

#endif