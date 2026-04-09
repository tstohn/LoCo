#pragma once



#ifdef LOCO_RCPP

  #include <Rcpp.h>

  #ifdef PI
  #undef PI
  #endif

  #define LOCO_OUT Rcpp::Rcout
  #define LOCO_ERR Rcpp::Rcerr

  // IMPORTANT: preserve message
  #define LOCO_EXIT(msg) Rcpp::stop("Stopping LoCo")

#else

  #include <iostream>
  #include <cstdlib>

  #define LOCO_OUT std::cout
  #define LOCO_ERR std::cerr

  // CLI: print message + exit
  #define LOCO_EXIT(msg) do { \
    std::cerr << msg << std::endl; \
    std::exit(EXIT_FAILURE); \
  } while(0)

#endif