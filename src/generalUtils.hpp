#pragma once

#include <iostream>
#include <string>
#include <vector>
#include "loco_io.h"

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

inline void printProgress(double percentage) 
{
    int val = (int) (percentage*100);
    int loadLength = (int) (percentage * PBWIDTH);
    int emptyLength = PBWIDTH - loadLength;
    LOCO_OUT << "\t\r[" << std::string(loadLength, '|') << std::string(emptyLength, ' ') << "] " << val << "%" << std::flush;
}

inline std::vector<int> smallerVector(const std::vector<int>& vec1, const std::vector<int>& vec2) 
{
    if (vec1.size() <= vec2.size()) 
    {
        return vec1;
    } else 
    {
        return vec2;
    }
}