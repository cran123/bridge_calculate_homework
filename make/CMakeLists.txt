cmake_minimum_required(VERSION 3.10)

PROJECT(stap++)

set (CMAKE_CXX_STANDARD 11)

OPTION(STAP++_DEBUG "Output detailed results for debugging." OFF)
OPTION(STAPPP_USE_EIGEN_SOLVER "Use Eigen SimplicialLDLT sparse solver." ON)
IF(STAP++_DEBUG)
   ADD_DEFINITIONS(-D_DEBUG_)
ENDIF()

set(STAPPP_SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/../src)

INCLUDE_DIRECTORIES(${STAPPP_SRC_DIR}/h)

if(STAPPP_USE_EIGEN_SOLVER)
   find_package(Eigen3 QUIET)
   if(NOT Eigen3_FOUND)
      set(STAPPP_BUNDLED_EIGEN_DIR ${CMAKE_CURRENT_LIST_DIR}/../third_party/eigen)
      if(EXISTS ${STAPPP_BUNDLED_EIGEN_DIR}/Eigen/Sparse)
         add_library(Eigen3::Eigen INTERFACE IMPORTED)
         set_target_properties(Eigen3::Eigen PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES ${STAPPP_BUNDLED_EIGEN_DIR})
      else()
         message(FATAL_ERROR
            "Eigen3 was not found. Install Eigen3 or add it as a submodule at third_party/eigen.")
      endif()
   endif()
   add_definitions(-DSTAPPP_USE_EIGEN)
   message(STATUS "Using Eigen3 sparse solver")
else()
   message(STATUS "Using legacy skyline LDLT solver")
endif()

AUX_SOURCE_DIRECTORY(${STAPPP_SRC_DIR}/cpp SRC)
list(FILTER SRC EXCLUDE REGEX "Hex8ConsistentLoad\\.cpp$")
FILE(GLOB_RECURSE HEAD ${STAPPP_SRC_DIR}/h/*.h)

source_group(SOURCE\ FILES FILES ${SRC})
source_group(HEADER\ FILES FILES ${HEAD})

ADD_EXECUTABLE(stap++ ${SRC} ${HEAD})
if(STAPPP_USE_EIGEN_SOLVER)
   target_link_libraries(stap++ Eigen3::Eigen)
endif()
ADD_EXECUTABLE(hex8_consistent_load ${STAPPP_SRC_DIR}/cpp/Hex8ConsistentLoad.cpp)
