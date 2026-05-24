#pragma once
#include "loco_io.h"

#include <iostream>
#include <string>
#include <vector>

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

inline void printProgress(double percentage) 
{

    //for RCPP do not print progress/ not thread safe
    #ifdef LOCO_RCPP
    return;
    #endif

    int val = (int) (percentage*100);
    //only print every 2 percent points
    if (val % 10 == 0) 
    {
        int loadLength = (int) (percentage * PBWIDTH);
        int emptyLength = PBWIDTH - loadLength;
        LOCO_OUT << "\t\r[" << std::string(loadLength, '|') << std::string(emptyLength, ' ') << "] " << val << "%" << std::flush;
    }
}

//cahnged the name to amke sure this is ONLY called in main thread, as this is not R-thread safe
inline void printProgressMainThread(double percentage) 
{
    int val = (int) (percentage * 100);
    
    // Print every 2 percent points
    if (val % 5 == 0) 
    {
        int loadLength = (int) (percentage * PBWIDTH);
        int emptyLength = PBWIDTH - loadLength;
        
        // LOCO_OUT automatically routes to Rcpp::Rcout or std::cout!
        LOCO_OUT << "\t\r[" << std::string(loadLength, '|') 
                 << std::string(emptyLength, ' ') << "] " << val << "%" << std::flush;
        
        // Force R/RStudio to draw the carriage return to the screen immediately
        #ifdef LOCO_RCPP
        R_FlushConsole(); 
        #endif
    }
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