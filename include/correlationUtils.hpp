#pragma once
#include "loco_io.h"

#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <numeric>
#include <algorithm>

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
static void old_rankify(std::vector<double>& X) 
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

static void rankify(std::vector<double>& X) {
    size_t n = X.size();
    if (n <= 1) return;

    // Create a vector of indices [0, 1, 2, ... n-1]
    std::vector<size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);

    // Sort indices based on the values in X
    std::sort(idx.begin(), idx.end(), [&](size_t i, size_t j) {
        return X[i] < X[j];
    });

    std::vector<double> ranks(n);
    for (size_t i = 0; i < n; ) {
        size_t j = i + 1;
        // Find the range of identical values (ties)
        while (j < n && X[idx[j]] == X[idx[i]]) {
            j++;
        }

        // Calculate fractional rank: (start_rank + end_rank) / 2
        // We use i+1 and j because ranks are 1-based
        double fractional_rank = (i + 1 + j) / 2.0;
        
        for (size_t k = i; k < j; k++) {
            ranks[idx[k]] = fractional_rank;
        }
        i = j;
    }
    X = std::move(ranks);
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
inline double OLD_calcualte_correlation_coefficient(const std::vector<double>& A, const std::vector<double>& B, bool ranking = false)
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
    double denom_part1 = n * squareSum_X - sum_X * sum_X;
    double denom_part2 = n * squareSum_Y - sum_Y * sum_Y;

    if (denom_part1 <= 0 || denom_part2 <= 0)
        return 0.0;  // or NaN depending on your preference


    double corr = (double)(n * sum_XY -  sum_X * sum_Y) / 
                    sqrt(denom_part1 * denom_part2);
    return corr;
}

inline double calculate_correlation_coefficient(const std::vector<double>& X,
                                               const std::vector<double>& Y) {
    const size_t n = X.size();
    
    // Basic requirement for correlation
    if (n < 2 || n != Y.size()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // 1. Calculate Means
    double sum_x = 0.0, sum_y = 0.0;
    for (size_t i = 0; i < n; ++i) {
        sum_x += X[i];
        sum_y += Y[i];
    }
    double mu_x = sum_x / n;
    double mu_y = sum_y / n;

    // 2. Calculate Covariance and Variances
    double ss_xy = 0.0; // Sum of squares xy (numerator)
    double ss_xx = 0.0; // Sum of squares x
    double ss_yy = 0.0; // Sum of squares y

    for (size_t i = 0; i < n; ++i) {
        double dx = X[i] - mu_x;
        double dy = Y[i] - mu_y;
        ss_xy += dx * dy;
        ss_xx += dx * dx;
        ss_yy += dy * dy;
    }

    // 3. Robust Zero-Variance Guard
    // We check if the variance is effectively zero. 
    // Using a very small epsilon to allow for floating point noise in ranks.
    if (ss_xx < 1e-18 || ss_yy < 1e-18) 
    {
        return std::numeric_limits<double>::quiet_NaN();
        // return zero for correlations with no variance
        //return 0.0;
    }

    // 4. Final Calculation
    double denom = std::sqrt(ss_xx) * std::sqrt(ss_yy);
    double r = ss_xy / denom;

    // Numerical clamping to ensure results stay within [-1, 1]
    return std::max(-1.0, std::min(1.0, r));
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