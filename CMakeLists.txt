cmake_minimum_required(VERSION 3.16)
project(ts LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the build type" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

find_package(OpenMP REQUIRED)

# pybind11 is cloned once under the user's home directory. Because this is an
# out-of-tree source directory, add_subdirectory also needs a build directory.
set(PYBIND11_ROOT "/home/research/w.yanwang/pybind11" CACHE PATH
    "Path to the cloned pybind11 repository")

if(NOT EXISTS "${PYBIND11_ROOT}/CMakeLists.txt")
    message(FATAL_ERROR
        "pybind11 was not found at ${PYBIND11_ROOT}. Clone it there first.")
endif()

set(PYBIND11_FINDPYTHON ON)
add_subdirectory("${PYBIND11_ROOT}" "${CMAKE_BINARY_DIR}/pybind11")

# ---------- HDF5 ----------
if(DEFINED ENV{HDF5_ROOT})
    set(HDF5_ROOT "$ENV{HDF5_ROOT}")
endif()
find_package(HDF5 REQUIRED COMPONENTS C CXX HL)

message(STATUS "HDF5 include dirs: ${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5 libraries: ${HDF5_LIBRARIES}")

# ---------- HighFive (header-only) ----------
set(HIGHFIVE_ROOT "/home/research/w.yanwang/HighFive-2.10.1")
set(HIGHFIVE_INCLUDE_DIR "${HIGHFIVE_ROOT}/include")

# ---------- Gurobi ----------
set(GUROBI_HOME "$ENV{GUROBI_HOME}")
if(NOT GUROBI_HOME)
    message(FATAL_ERROR
        "GUROBI_HOME is not set. Please set it in your environment.")
endif()

set(GUROBI_INCLUDE_DIR "${GUROBI_HOME}/include")
set(GUROBI_LIBRARY_DIR "${GUROBI_HOME}/lib")
set(GUROBI_LIB "${GUROBI_LIBRARY_DIR}/libgurobi120.so")
set(GUROBI_CXX_LIB "${GUROBI_LIBRARY_DIR}/libgurobi_g++8.5.a")

message(STATUS "Gurobi Library: ${GUROBI_LIB}")
message(STATUS "Gurobi C++ Library: ${GUROBI_CXX_LIB}")

# ---------- CPLEX ----------
set(CPLEX_ROOT "/project/engineering/cplex22/cplex")
set(CPLEX_LIB_DIR "${CPLEX_ROOT}/lib/x86-64_linux/static_pic")
set(CONCERT_LIB_DIR
    "/project/engineering/cplex22/concert/lib/x86-64_linux/static_pic")
set(CPLEX_INCLUDE_DIR "${CPLEX_ROOT}/include")
set(CPLEX_CONCERT_INCLUDE_DIR
    "/project/engineering/cplex22/concert/include")
set(ILOCLEX_LIB "${CPLEX_LIB_DIR}/libilocplex.a")
set(CPLEX_LIB "${CPLEX_LIB_DIR}/libcplex.a")
set(CONCERT_LIB "${CONCERT_LIB_DIR}/libconcert.a")

# ---------- LEMON ----------
set(LEMON_SOURCE_DIR "/home/research/w.yanwang/lemon-1.3.1")
set(LEMON_BUILD_DIR "/home/research/w.yanwang/lemon-1.3.1/build")

# All ordinary implementation files are under src/. The root-level main.cpp
# is deliberately not included here.
file(GLOB_RECURSE TS_CORE_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

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

# Compile all reusable project code once. POSITION_INDEPENDENT_CODE is needed
# because this static library is also linked into a Python shared module.
add_library(ts_core STATIC ${TS_CORE_SOURCES})
set_target_properties(ts_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(ts_core PUBLIC ${TS_INCLUDE_DIRS})
target_link_libraries(ts_core PUBLIC ${TS_LINK_LIBRARIES})

if(HDF5_LIBRARY_DIRS)
    target_link_directories(ts_core PUBLIC ${HDF5_LIBRARY_DIRS})
endif()

# Existing command-line executable. main.cpp is at the project root.
add_executable(ts "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")
target_link_libraries(ts PRIVATE ts_core)

# Python extension module. These three files are outside src/, so they are not
# accidentally included in ts_core.
set(M4OPT_SOLVERS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/m4opt_solvers")

pybind11_add_module(m4opt_solvers MODULE
    "${M4OPT_SOLVERS_DIR}/bindings.cpp"
    "${M4OPT_SOLVERS_DIR}/m4opt_solvers.cpp"
)

target_include_directories(m4opt_solvers PRIVATE
    "${M4OPT_SOLVERS_DIR}"
)
target_link_libraries(m4opt_solvers PRIVATE ts_core)

# The executable is build/ts. The importable extension is written under
# build/python as m4opt_solvers.cpython-<version>-<platform>.so.
set_target_properties(m4opt_solvers PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
    BUILD_RPATH "${GUROBI_LIBRARY_DIR};${HDF5_LIBRARY_DIRS}"
    INSTALL_RPATH "${GUROBI_LIBRARY_DIR};${HDF5_LIBRARY_DIRS}"
)

set_target_properties(ts PROPERTIES
    BUILD_RPATH "${GUROBI_LIBRARY_DIR};${HDF5_LIBRARY_DIRS}"
    INSTALL_RPATH "${GUROBI_LIBRARY_DIR};${HDF5_LIBRARY_DIRS}"
)
