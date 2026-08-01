cmake_minimum_required(VERSION 3.11)
project(ts)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


find_package(OpenMP REQUIRED)

# ---------- HDF5 ----------
# set(HDF5_ROOT "/project/engineering/anaconda3")
# find_package(HDF5 REQUIRED COMPONENTS C CXX HL)
if(DEFINED ENV{HDF5_ROOT})
  set(HDF5_ROOT "$ENV{HDF5_ROOT}")
endif()
find_package(HDF5 REQUIRED COMPONENTS C CXX HL)

message(STATUS "HDF5 include dirs: ${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5 libraries: ${HDF5_LIBRARIES}")

# ---------- HighFive (header-only) ----------
set(HIGHFIVE_ROOT "/home/research/w.yanwang/HighFive-2.10.1")
set(HIGHFIVE_INCLUDE_DIR "${HIGHFIVE_ROOT}/include")

# -------------------------
# Gurobi config
# -------------------------
set(GUROBI_HOME "$ENV{GUROBI_HOME}")
if(NOT GUROBI_HOME)
  message(FATAL_ERROR "GUROBI_HOME is not set. Please set it in your environment.")
endif()
set(GUROBI_INCLUDE_DIR "${GUROBI_HOME}/include")
set(GUROBI_LIBRARY_DIR "${GUROBI_HOME}/lib")
set(GUROBI_LIB "${GUROBI_LIBRARY_DIR}/libgurobi120.so")
set(GUROBI_CXX_LIB "${GUROBI_LIBRARY_DIR}/libgurobi_g++8.5.a")

message(STATUS "Gurobi Library: ${GUROBI_LIB}")
message(STATUS "Gurobi C++ Library: ${GUROBI_CXX_LIB}")

set(CMAKE_INSTALL_RPATH ${GUROBI_LIBRARY_DIR})

# -------------------------
# CPLEX config
# -------------------------
set(CPLEX_ROOT "/project/engineering/cplex22/cplex")
set(CPLEX_LIB_DIR "${CPLEX_ROOT}/lib/x86-64_linux/static_pic")
set(CONCERT_LIB_DIR "/project/engineering/cplex22/concert/lib/x86-64_linux/static_pic")
set(CPLEX_INCLUDE_DIR "${CPLEX_ROOT}/include")
set(CPLEX_CONCERT_INCLUDE_DIR "/project/engineering/cplex22/concert/include")
set(ILOCLEX_LIB "${CPLEX_LIB_DIR}/libilocplex.a")
set(CPLEX_LIB "${CPLEX_LIB_DIR}/libcplex.a")
set(CONCERT_LIB "${CONCERT_LIB_DIR}/libconcert.a")

# -------------------------
# LEMON config
# -------------------------
set(LEMON_SOURCE_DIR "/home/research/w.yanwang/lemon-1.3.1")
set(LEMON_BUILD_DIR  "/home/research/w.yanwang/lemon-1.3.1/build")

# Collect all source files (works with your new src/ layout)
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

add_executable(ts ${SOURCES})

# target_include_directories(ts PRIVATE
#   ${CMAKE_CURRENT_SOURCE_DIR}/include
#   ${GUROBI_INCLUDE_DIR}
#   ${CPLEX_INCLUDE_DIR}
#   ${CPLEX_CONCERT_INCLUDE_DIR}
#   ${LEMON_SOURCE_DIR}
#   ${LEMON_BUILD_DIR}
#   ${HIGHFIVE_INCLUDE_DIR}
#   ${HDF5_INCLUDE_DIRS}
# )

target_include_directories(ts PRIVATE
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

target_link_libraries(ts PRIVATE
  OpenMP::OpenMP_CXX
  ${GUROBI_LIB}
  ${GUROBI_CXX_LIB}
  ${ILOCLEX_LIB}
  ${CONCERT_LIB}
  ${CPLEX_LIB}
  ${HDF5_LIBRARIES}
  pthread m dl
  stdc++fs
)

# If HDF5 needs extra link dirs on your system, this helps:
target_link_directories(ts PRIVATE ${HDF5_LIBRARY_DIRS})





# cmake_minimum_required(VERSION 3.11)
# project(ts)

# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

# if(NOT CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
# endif()
# message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# # option(ENABLE_ASAN "Enable AddressSanitizer and UBSan" OFF)

# find_package(OpenMP REQUIRED)

# # ---------- HDF5 ----------
# set(HDF5_ROOT "/project/engineering/anaconda3")
# find_package(HDF5 REQUIRED COMPONENTS C CXX HL)

# message(STATUS "HDF5 include dirs: ${HDF5_INCLUDE_DIRS}")
# message(STATUS "HDF5 libraries: ${HDF5_LIBRARIES}")

# # ---------- HighFive ----------
# set(HIGHFIVE_ROOT "/home/research/w.yanwang/HighFive-2.10.1")
# set(HIGHFIVE_INCLUDE_DIR "${HIGHFIVE_ROOT}/include")

# # ---------- Gurobi ----------
# set(GUROBI_HOME "$ENV{GUROBI_HOME}")
# if(NOT GUROBI_HOME)
#   message(FATAL_ERROR "GUROBI_HOME is not set. Please set it in your environment.")
# endif()

# set(GUROBI_INCLUDE_DIR "${GUROBI_HOME}/include")
# set(GUROBI_LIBRARY_DIR "${GUROBI_HOME}/lib")
# set(GUROBI_LIB "${GUROBI_LIBRARY_DIR}/libgurobi120.so")
# set(GUROBI_CXX_LIB "${GUROBI_LIBRARY_DIR}/libgurobi_g++8.5.a")

# message(STATUS "Gurobi Library: ${GUROBI_LIB}")
# message(STATUS "Gurobi C++ Library: ${GUROBI_CXX_LIB}")

# set(CMAKE_INSTALL_RPATH ${GUROBI_LIBRARY_DIR})

# # ---------- CPLEX ----------
# set(CPLEX_ROOT "/project/engineering/cplex22/cplex")
# set(CPLEX_LIB_DIR "${CPLEX_ROOT}/lib/x86-64_linux/static_pic")
# set(CONCERT_LIB_DIR "/project/engineering/cplex22/concert/lib/x86-64_linux/static_pic")
# set(CPLEX_INCLUDE_DIR "${CPLEX_ROOT}/include")
# set(CPLEX_CONCERT_INCLUDE_DIR "/project/engineering/cplex22/concert/include")

# set(ILOCLEX_LIB "${CPLEX_LIB_DIR}/libilocplex.a")
# set(CPLEX_LIB "${CPLEX_LIB_DIR}/libcplex.a")
# set(CONCERT_LIB "${CONCERT_LIB_DIR}/libconcert.a")

# # ---------- LEMON ----------
# set(LEMON_SOURCE_DIR "/home/research/w.yanwang/lemon-1.3.1")
# set(LEMON_BUILD_DIR  "/home/research/w.yanwang/lemon-1.3.1/build")

# # ---------- Sources ----------
# file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS
#   "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
# )

# add_executable(ts ${SOURCES})

# # # ---------- Sanitizers ----------
# # if(ENABLE_ASAN)
# #   if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
# #     message(STATUS "Enabling AddressSanitizer + UndefinedBehaviorSanitizer")

# #     target_compile_options(ts PRIVATE
# #       -g
# #       -O1
# #       -fno-omit-frame-pointer
# #       -fsanitize=address,undefined
# #     )

# #     target_link_libraries(ts PRIVATE
# #       -fsanitize=address,undefined
# #     )
# #   else()
# #     message(WARNING "ASan/UBSan requested, but compiler may not support it.")
# #   endif()
# # endif()

# target_include_directories(ts PRIVATE
#   ${CMAKE_CURRENT_SOURCE_DIR}/include
#   ${CMAKE_CURRENT_SOURCE_DIR}/include/common
#   ${CMAKE_CURRENT_SOURCE_DIR}/include/single_telescope
#   ${CMAKE_CURRENT_SOURCE_DIR}/include/multi_telescope/min_telescopes
#   ${CMAKE_CURRENT_SOURCE_DIR}/include/multi_telescope/max_probabilities
#   ${CMAKE_CURRENT_SOURCE_DIR}/include/apps

#   ${GUROBI_INCLUDE_DIR}
#   ${CPLEX_INCLUDE_DIR}
#   ${CPLEX_CONCERT_INCLUDE_DIR}
#   ${LEMON_SOURCE_DIR}
#   ${LEMON_BUILD_DIR}
#   ${HIGHFIVE_INCLUDE_DIR}
#   ${HDF5_INCLUDE_DIRS}
# )

# target_link_libraries(ts PRIVATE
#   OpenMP::OpenMP_CXX
#   ${GUROBI_LIB}
#   ${GUROBI_CXX_LIB}
#   ${ILOCLEX_LIB}
#   ${CONCERT_LIB}
#   ${CPLEX_LIB}
#   ${HDF5_LIBRARIES}
#   pthread
#   m
#   dl
#   stdc++fs
# )

# target_link_directories(ts PRIVATE ${HDF5_LIBRARY_DIRS})