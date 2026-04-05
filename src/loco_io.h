#pragma once

#ifdef LOCO_RCPP

  #include <Rcpp.h>

  #define LOCO_OUT Rcpp::Rcout
  #define LOCO_ERR Rcpp::Rcerr

  #define LOCO_EXIT(msg) Rcpp::stop(msg)

#else

  #include <iostream>
  #include <cstdlib>

  #define LOCO_OUT std::cout
  #define LOCO_ERR std::cerr

  #define LOCO_EXIT(msg) { \
    std::cerr << msg << std::endl; \
    std::exit(1); \
  }

#endif