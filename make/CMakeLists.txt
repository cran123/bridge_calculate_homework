cmake_minimum_required(VERSION 3.10)

PROJECT(stap++)

set (CMAKE_CXX_STANDARD 11)

OPTION(STAP++_DEBUG "Output detailed results for debugging." OFF)
IF(STAP++_DEBUG)
   ADD_DEFINITIONS(-D_DEBUG_)
ENDIF()

set(STAPPP_SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/../src)

INCLUDE_DIRECTORIES(${STAPPP_SRC_DIR}/h)

AUX_SOURCE_DIRECTORY(${STAPPP_SRC_DIR}/cpp SRC)
list(FILTER SRC EXCLUDE REGEX "Hex8ConsistentLoad\\.cpp$")
FILE(GLOB_RECURSE HEAD ${STAPPP_SRC_DIR}/h/*.h)

source_group(SOURCE\ FILES FILES ${SRC})
source_group(HEADER\ FILES FILES ${HEAD})

ADD_EXECUTABLE(stap++ ${SRC} ${HEAD})
ADD_EXECUTABLE(hex8_consistent_load ${STAPPP_SRC_DIR}/cpp/Hex8ConsistentLoad.cpp)
