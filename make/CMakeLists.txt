cmake_minimum_required(VERSION 3.10)

PROJECT(stap++)

set (CMAKE_CXX_STANDARD 11)

OPTION(STAP++_DEBUG "Output detailed results for debugging." OFF)
OPTION(STAPPP_USE_EIGEN_SOLVER "Use Eigen SimplicialLDLT sparse solver." ON)
OPTION(STAPPP_USE_PARDISO_SOLVER "Use Intel oneMKL PARDISO sparse solver when available." ON)
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
   if(STAPPP_USE_PARDISO_SOLVER)
      message(WARNING "PARDISO requires the Eigen sparse solver path; disabling PARDISO.")
      set(STAPPP_USE_PARDISO_SOLVER OFF)
   endif()
   message(STATUS "Using legacy skyline LDLT solver")
endif()

if(STAPPP_USE_PARDISO_SOLVER)
   find_package(MKL CONFIG QUIET HINTS "$ENV{MKLROOT}/lib/cmake/mkl")
   if(MKL_FOUND)
      add_definitions(-DSTAPPP_USE_PARDISO)
      message(STATUS "Using Intel oneMKL PARDISO solver by default")
   else()
      message(WARNING
         "Intel oneMKL was not found. PARDISO support is disabled; "
         "set MKL_DIR or source the oneMKL environment before configuring.")
      set(STAPPP_USE_PARDISO_SOLVER OFF)
   endif()
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
if(STAPPP_USE_PARDISO_SOLVER)
   target_link_libraries(stap++ MKL::MKL)
endif()
ADD_EXECUTABLE(hex8_consistent_load ${STAPPP_SRC_DIR}/cpp/Hex8ConsistentLoad.cpp)
