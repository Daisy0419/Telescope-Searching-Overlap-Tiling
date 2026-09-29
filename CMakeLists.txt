cmake_minimum_required(VERSION 3.16)
project(ts LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the build type" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

find_package(OpenMP REQUIRED)

# Common root for locally installed software.
set(SOFTWARE_ROOT "/home/daisy/Software" CACHE PATH
    "Root directory containing locally installed software")

# ---------- pybind11 ----------
set(PYBIND11_ROOT "${SOFTWARE_ROOT}/pybind11" CACHE PATH
    "Path to the cloned pybind11 repository")

if(NOT EXISTS "${PYBIND11_ROOT}/CMakeLists.txt")
    message(FATAL_ERROR
        "pybind11 was not found at ${PYBIND11_ROOT}.")
endif()

set(PYBIND11_FINDPYTHON ON)
add_subdirectory(
    "${PYBIND11_ROOT}"
    "${CMAKE_BINARY_DIR}/pybind11"
)

# ---------- Gurobi ----------
if(DEFINED ENV{GUROBI_HOME}
   AND NOT "$ENV{GUROBI_HOME}" STREQUAL "")
    set(GUROBI_DEFAULT_HOME "$ENV{GUROBI_HOME}")
else()
    set(GUROBI_DEFAULT_HOME
        "${SOFTWARE_ROOT}/gurobi1200/linux64")
endif()

set(GUROBI_HOME "${GUROBI_DEFAULT_HOME}" CACHE PATH
    "Path to the Gurobi linux64 directory")

set(GUROBI_INCLUDE_DIR "${GUROBI_HOME}/include")
set(GUROBI_LIBRARY_DIR "${GUROBI_HOME}/lib")
set(GUROBI_LIB "${GUROBI_LIBRARY_DIR}/libgurobi120.so")

# Retain the C++ library filename from your original configuration.
# Some Gurobi installations call it libgurobi_c++.a instead.
if(EXISTS "${GUROBI_LIBRARY_DIR}/libgurobi_g++8.5.a")
    set(GUROBI_CXX_LIB
        "${GUROBI_LIBRARY_DIR}/libgurobi_g++8.5.a")
elseif(EXISTS "${GUROBI_LIBRARY_DIR}/libgurobi_c++.a")
    set(GUROBI_CXX_LIB
        "${GUROBI_LIBRARY_DIR}/libgurobi_c++.a")
else()
    message(FATAL_ERROR
        "Gurobi C++ library was not found in "
        "${GUROBI_LIBRARY_DIR}")
endif()

if(NOT EXISTS "${GUROBI_LIB}")
    message(FATAL_ERROR
        "Gurobi library was not found: ${GUROBI_LIB}")
endif()

message(STATUS "Gurobi home: ${GUROBI_HOME}")
message(STATUS "Gurobi library: ${GUROBI_LIB}")
message(STATUS "Gurobi C++ library: ${GUROBI_CXX_LIB}")

# ---------- CPLEX ----------
set(CPLEX_STUDIO_ROOT
    "${SOFTWARE_ROOT}/ibm/ILOG/CPLEX_Studio2212"
    CACHE PATH "Path to CPLEX Optimization Studio")

set(CPLEX_ROOT "${CPLEX_STUDIO_ROOT}/cplex")
set(CONCERT_ROOT "${CPLEX_STUDIO_ROOT}/concert")

set(CPLEX_LIB_DIR
    "${CPLEX_ROOT}/lib/x86-64_linux/static_pic")
set(CONCERT_LIB_DIR
    "${CONCERT_ROOT}/lib/x86-64_linux/static_pic")

set(CPLEX_INCLUDE_DIR "${CPLEX_ROOT}/include")
set(CPLEX_CONCERT_INCLUDE_DIR "${CONCERT_ROOT}/include")

set(ILOCLEX_LIB "${CPLEX_LIB_DIR}/libilocplex.a")
set(CPLEX_LIB "${CPLEX_LIB_DIR}/libcplex.a")
set(CONCERT_LIB "${CONCERT_LIB_DIR}/libconcert.a")

foreach(CPLEX_REQUIRED_FILE
        "${ILOCLEX_LIB}"
        "${CPLEX_LIB}"
        "${CONCERT_LIB}")
    if(NOT EXISTS "${CPLEX_REQUIRED_FILE}")
        message(FATAL_ERROR
            "Required CPLEX library was not found: "
            "${CPLEX_REQUIRED_FILE}")
    endif()
endforeach()

message(STATUS "CPLEX Studio: ${CPLEX_STUDIO_ROOT}")
message(STATUS "CPLEX library: ${CPLEX_LIB}")
message(STATUS "Concert library: ${CONCERT_LIB}")

# ---------- LEMON ----------
set(LEMON_SOURCE_DIR
    "${SOFTWARE_ROOT}/lemon-1.3.1"
    CACHE PATH "Path to LEMON source directory")

set(LEMON_BUILD_DIR
    "${LEMON_SOURCE_DIR}/build"
    CACHE PATH "Path to LEMON build directory")

# # All ordinary implementation files are under src/.
# file(GLOB_RECURSE TS_CORE_SOURCES CONFIGURE_DEPENDS
#     "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
# )

# ---------- Sources ----------
# This is the only permitted main.cpp.
set(TS_MAIN_SOURCE
    "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp"
)

if(NOT EXISTS "${TS_MAIN_SOURCE}")
    message(FATAL_ERROR
        "Required entry point was not found: ${TS_MAIN_SOURCE}")
endif()

# Collect all C++ implementation files under src/.
file(GLOB_RECURSE TS_CORE_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

# main.cpp belongs only to the executable, not to the reusable core library.
list(REMOVE_ITEM TS_CORE_SOURCES "${TS_MAIN_SOURCE}")

set(TS_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/common
    ${CMAKE_CURRENT_SOURCE_DIR}/include/single_telescope/no_overlap
    ${CMAKE_CURRENT_SOURCE_DIR}/include/single_telescope/pixel_overlap
    ${CMAKE_CURRENT_SOURCE_DIR}/include/multi_telescope/min_telescopes
    ${CMAKE_CURRENT_SOURCE_DIR}/include/multi_telescope/max_probabilities
    ${CMAKE_CURRENT_SOURCE_DIR}/include/apps
    ${GUROBI_INCLUDE_DIR}
    ${CPLEX_INCLUDE_DIR}
    ${CPLEX_CONCERT_INCLUDE_DIR}
    ${LEMON_SOURCE_DIR}
    ${LEMON_BUILD_DIR}
    ${HIGHFIVE_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIRS}
)

set(TS_LINK_LIBRARIES
    OpenMP::OpenMP_CXX
    ${GUROBI_CXX_LIB}
    ${GUROBI_LIB}
    ${ILOCLEX_LIB}
    ${CONCERT_LIB}
    ${CPLEX_LIB}
    ${HDF5_LIBRARIES}
    pthread
    m
    dl
    stdc++fs
)

add_library(ts_core STATIC ${TS_CORE_SOURCES})
set_target_properties(ts_core PROPERTIES
    POSITION_INDEPENDENT_CODE ON
)

target_include_directories(ts_core PUBLIC ${TS_INCLUDE_DIRS})
target_link_libraries(ts_core PUBLIC ${TS_LINK_LIBRARIES})

# # Command-line executable.
# add_executable(ts "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")
# target_link_libraries(ts PRIVATE ts_core)
# Command-line executable: src/main.cpp is the only entry point.
add_executable(ts "${TS_MAIN_SOURCE}")
target_link_libraries(ts PRIVATE ts_core)

# Python extension module.
set(M4OPT_SOLVERS_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/m4opt_solvers")

pybind11_add_module(m4opt_solvers MODULE
    "${M4OPT_SOLVERS_DIR}/bindings.cpp"
    "${M4OPT_SOLVERS_DIR}/m4opt_solvers.cpp"
)

target_include_directories(m4opt_solvers PRIVATE
    "${M4OPT_SOLVERS_DIR}"
)

target_link_libraries(m4opt_solvers PRIVATE ts_core)

set_target_properties(m4opt_solvers PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
    BUILD_RPATH "${GUROBI_LIBRARY_DIR}"
    INSTALL_RPATH "${GUROBI_LIBRARY_DIR}"
)

set_target_properties(ts PROPERTIES
    BUILD_RPATH "${GUROBI_LIBRARY_DIR}"
    INSTALL_RPATH "${GUROBI_LIBRARY_DIR}"
)