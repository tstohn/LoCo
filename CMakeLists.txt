cmake_minimum_required(VERSION 3.20)
project(LoCo LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ===============================
# Default build type
# ===============================
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# ===============================
# Flags (match Makefile)
# ===============================
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra)
endif()

# LTO
include(CheckIPOSupported)
check_ipo_supported(RESULT lto_supported)
if(lto_supported)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# ===============================
# Dependencies (via vcpkg)
# ===============================
find_package(Boost REQUIRED COMPONENTS program_options iostreams)
find_package(igraph CONFIG REQUIRED)
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

find_library(ARPACK_LIB arpack)

# ===============================
# Include dirs
# ===============================
include_directories(
    ${PROJECT_SOURCE_DIR}/src
)

file(GLOB_RECURSE ALL_INCLUDE_DIRS LIST_DIRECTORIES true ${PROJECT_SOURCE_DIR}/src/*)
foreach(dir ${ALL_INCLUDE_DIRS})
    if(IS_DIRECTORY ${dir})
        include_directories(${dir})
    endif()
endforeach()

# ===============================
# Sources
# ===============================
set(SOURCES
    src/LoCo.cpp
    src/lib/Parser/SCParser.cpp
    src/lib/SCGraph/Constructors/GraphData.cpp
    src/lib/SCGraph/Constructors/GraphHandler.cpp
    src/lib/SCGraph/Constructors/Neighborhood.cpp
)

add_executable(loco ${SOURCES})

# ===============================
# Linking
# ===============================
target_link_libraries(loco PRIVATE
    Boost::program_options
    Boost::iostreams
    igraph::igraph
    Threads::Threads
    ZLIB::ZLIB
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${ARPACK_LIB}
)

if(WIN32)
    target_link_libraries(loco PRIVATE ws2_32)
endif()

# ===============================
# nanoflann (header-only)
# ===============================
include(FetchContent)

FetchContent_Declare(
    nanoflann
    GIT_REPOSITORY https://github.com/jlblancoc/nanoflann.git
    GIT_TAG v1.3.2
)

FetchContent_MakeAvailable(nanoflann)
target_include_directories(loco PRIVATE ${nanoflann_SOURCE_DIR}/include)

# ===============================
# Output
# ===============================
set_target_properties(loco PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin
)