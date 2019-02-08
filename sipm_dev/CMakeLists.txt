# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(AnaSipm)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

SET(GCC_COVERAGE_COMPILE_FLAGS "-std=c++11 -O3 -Wall -Wextra")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )

ROOT_GENERATE_DICTIONARY(G__AnaSipm AnaSipm.h LINKDEF AnaSipmLinkDef.h)

#---Create a shared library with geneated dictionary
add_library(AnaSipm SHARED AnaSipm.cxx G__AnaSipm.cxx)
target_link_libraries(AnaSipm ${ROOT_LIBRARIES})

#---Create  a main program using the library
add_executable(decode_file decode_file.cxx)
add_executable(characterize_sipm characterize_sipm.cxx)
target_link_libraries(decode_file AnaSipm)
target_link_libraries(characterize_sipm AnaSipm)