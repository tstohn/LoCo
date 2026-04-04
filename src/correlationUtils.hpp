#pragma once

#include <iostream>
#include <string>
#include <vector>
#include "loco_io.h"

namespace
{

inline std::vector<std::vector<double>> filterVector(const std::vector<std::vector<double>>& original, const std::vector<int>& indices) 
{
    std::vector<std::vector<double>> result;

    for (size_t index : indices) 
    {
        // Check if the index is within bounds
        if (index < original.size()) 
        {
            // Add the corresponding element to the result vector
            result.push_back(original.at(index));
        } 
        else 
        {
            // Handle out-of-bounds index if needed
            LOCO_ERR << "Warning: Index " << index << " is out of bounds." << std::endl;
        }
    }

    return result;
}

// Function returns the rank vector of the set of observations: used from geeksforgeeks: https://www.geeksforgeeks.org/program-spearmans-rank-correlation/
static void rankify(std::vector<double>& X) 
{
    int N = X.size();

    // Rank Vector
    std::vector<double> Rank_X(N);
    for(int i = 0; i < N; i++) 
    {
        int r = 1, s = 1;
        
        // Count no of smaller elements
        // in 0 to i-1
        for(int j = 0; j < i; j++) 
        {
            if (X[j] < X[i] ) r++;
            if (X[j] == X[i] ) s++;
        }
    
        // Count no of smaller elements
        // in i+1 to N-1
        for (int j = i+1; j < N; j++) 
        {
            if (X[j] < X[i] ) r++;
            if (X[j] == X[i] ) s++;
        }

        // Use Fractional Rank formula
        // fractional_rank = r + (n-1)/2
        Rank_X[i] = r + (s-1) * 0.5;        
    }
    
    // Return Rank Vector
    X = Rank_X;
}

}

inline void remove_subsets_sota(std::vector<std::vector<int>>& vectors) {
    if (vectors.empty()) return;

    // 1. Sort elements within each vector (Required for fast is_subset)
    for (auto& v : vectors) {
        std::sort(v.begin(), v.end());
    }

    // 2. Sort vectors by size DESCENDING (Large sets first)
    std::sort(vectors.begin(), vectors.end(), 
        [](const std::vector<int>& a, const std::vector<int>& b) {
            return a.size() > b.size(); 
        });

    std::vector<bool> toRemove(vectors.size(), false);
    
    // 3. Parallel Subset Check (O(N^2) but with massive pruning)
    // We use C++17 parallel policy. If not available, use a standard loop.
    for (size_t i = 1; i < vectors.size(); ++i) {
        // Only check if it's a subset of a LARGER set (indices 0 to i-1)
        for (size_t j = 0; j < i; ++j) {
            if (toRemove[j]) continue; // If the container is already gone, skip

            // should never be true bcs of ordered/ but keep for clarity
            if (vectors[i].size() > vectors[j].size()) continue;

            if (std::includes(vectors[j].begin(), vectors[j].end(),
                              vectors[i].begin(), vectors[i].end())) {
                toRemove[i] = true;
                break;
            }
        }
    }

    // 4. Efficient Erase
    size_t writeIdx = 0;
    for (size_t readIdx = 0; readIdx < vectors.size(); ++readIdx) {
        if (!toRemove[readIdx]) {
            if (readIdx != writeIdx) {
                vectors[writeIdx] = std::move(vectors[readIdx]);
            }
            writeIdx++;
        }
    }
    vectors.resize(writeIdx);
}

// function that returns Pearson correlation coefficient.
// for ranked correlations se raking to TRUE
inline double calcualte_correlation_coefficient(const std::vector<double>& A, const std::vector<double>& B, bool ranking = false)
{
    std::vector<double> X = A;
    std::vector<double> Y = B;

    if(ranking)
    {
        rankify(X);
        rankify(Y);
    }

    size_t n = X.size();
    double sum_X = 0, sum_Y = 0, sum_XY = 0;
    double squareSum_X = 0, squareSum_Y = 0;

    for (size_t i = 0; i < n; i++)
    {
        // sum of elements of array X.
        sum_X = sum_X + X[i];
        // sum of elements of array Y.
        sum_Y = sum_Y + Y[i];
        // sum of X[i] * Y[i].
        sum_XY = sum_XY + X[i] * Y[i];
        // sum of square of array elements.
        squareSum_X = squareSum_X + X[i] * X[i];
        squareSum_Y = squareSum_Y + Y[i] * Y[i];
    }

    // use formula for calculating correlation coefficient.
    double corr = (double)(n * sum_XY -  sum_X * sum_Y) / 
                    sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y - sum_Y * sum_Y));
    return corr;
}

//fit a linear function through points and return slope
inline double calculate_slope(const std::vector<double>& pointsA, const std::vector<double>& pointsB)
{
    if(pointsA.size() != pointsB.size())
    {
        LOCO_ERR << "In \'calculate_slope\' the size of the two vectors is not equal!!!\n";
        LOCO_EXIT(EXIT_FAILURE);
    }

    double sum_xy = 0;
    double sum_x = 0;
    double sum_y = 0;
    double sum_x_square = 0;
    double sum_y_square = 0;

    for(size_t i =0; i < pointsA.size(); ++i)
    {
        sum_xy += (pointsA.at(i) * pointsB.at(i));
        sum_x += pointsA.at(i);
        sum_y += pointsB.at(i);
        sum_x_square += (pointsA.at(i) * pointsA.at(i));
        sum_y_square += (pointsB.at(i) * pointsB.at(i));
    }

    double N = pointsA.size();
    double numerator = (N * sum_xy - sum_x * sum_y);
    double denominator = (N * sum_x_square - sum_x * sum_x);
    double coeff = numerator / denominator;

    return(coeff);
}