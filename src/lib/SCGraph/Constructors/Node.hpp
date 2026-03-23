#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <exception>
#include <numeric>
#include <vector>
#include <cmath>

class node;
typedef std::shared_ptr<node> nodePtr;

class node
{
    public:
        node(std::vector<double> values, std::string name = "") : values(values), name(name){};

        //dimensions of node
        int dimensions()
        {
            return(values.size());
        }
        std::string get_name()
        {
            return(name);
        }
        //values of node in the dimension dim
        double value_at(int dim)
        {
            return(values.at(dim));
        }
        std::vector<double> all_values()
        {
            return(values);
        }
        double distance_to(const nodePtr& n, std::string method = "manhattan")
        {
            double dist = 0;
            if(method == "euclidean")
            {
                for(int i = 0; i < n->values.size(); ++i)
                {
                    dist += std::pow((values.at(i) - n->values.at(i)),2);
                }
                dist = std::sqrt(dist);
            }
            else if(method == "manhattan")
            {
                for(int i = 0; i < n->values.size(); ++i)
                {
                    dist += std::abs(values.at(i) - n->values.at(i));
                }
            }
            return(dist);
        }
        double distance_to(const nodePtr& n, const std::vector<int>& genes, std::string method = "manhattan")
        {
            double dist = 0;
            if(method == "euclidean")
            {
                for(int gene_idx = 0; gene_idx < genes.size(); ++gene_idx)
                {
                    dist += std::pow((values.at(genes.at(gene_idx)) - n->values.at(genes.at(gene_idx))),2);
                }
                dist = std::sqrt(dist);
            }
            else if(method == "manhattan")
            {
                for(int gene_idx = 0; gene_idx < genes.size(); ++gene_idx)
                {
                    dist += std::abs(values.at(genes.at(gene_idx)) - n->values.at(genes.at(gene_idx)));
                }
            }
            return(dist);
        }

    private:
        std::vector<double> values;
        std::string name;
};

typedef std::vector<nodePtr> nodePtrVector;
typedef std::pair<std::shared_ptr<node>, double> nodeDistPair;
