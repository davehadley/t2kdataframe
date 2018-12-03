# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(t2kdataframe)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Debug)

set( T2KDATAFRAME_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src )
set( T2KDATAFRAME_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR} )

# Setup ROOT following instructions on https://root.cern.ch/how/integrate-root-my-project-cmake
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT 6.15 REQUIRED COMPONENTS RIO Net Tree ROOTVecOps ROOTDataFrame MathCore)
include_directories(${VDT_INCLUDE_DIRS})
#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

# Try to find a file to do make project with
set(MAKE_PROJECT_FILE_NAME "oaanalysisfile.root" CACHE FILEPATH "Path to oaAnalysis file used to by TFile::MakeProject.")
if(NOT EXISTS ${MAKE_PROJECT_FILE_NAME})
  message(FATAL_ERROR "ERROR: unable to find file to run TFile::MakeProject. Set this with -D MAKE_PROJECT_FILE_NAME=<filename.root>")
endif()
message(STATUS "using MAKE_PROJECT_FILE_NAME=${MAKE_PROJECT_FILE_NAME}")

# Run ROOT's TFile::MakeProject
set(MakeProjectLib ${T2KDATAFRAME_BINARY_DIR}/libT2KDataFrameMakeProject.so)
add_custom_command(
        OUTPUT ${MakeProjectLib}
        COMMAND $ENV{ROOTSYS}/bin/root -l -q ${MAKE_PROJECT_FILE_NAME} ${T2KDATAFRAME_SOURCE_DIR}/buildmakeproject.C && mv ${T2KDATAFRAME_BINARY_DIR}/libT2KDataFrameMakeProject/*.so ./
        DEPENDS ${MAKE_PROJECT_FILE_NAME}
)
add_custom_target(CustomMakeProjectTarget DEPENDS ${MakeProjectLib})

# Ensure all directories are accessible
include_directories( ${T2KDATAFRAME_BINARY_DIR})
include_directories( ${T2KDATAFRAME_BINARY_DIR}/libT2KDataFrameMakeProject)
include_directories("${T2KDATAFRAME_SOURCE_DIR}")

# Generate library to use NCFGD code in python
ROOT_GENERATE_DICTIONARY(T2KDataFrameDict src/T2KDataSource.hxx LINKDEF src/T2KDataSourceLinkDef.h DEPENDENCIES ${MakeProjectLib} CustomMakeProjectTarget)
add_library(T2KDataSource SHARED src/T2KDataSource.cxx ${MakeProjectLib} T2KDataFrameDict.cxx)
add_dependencies(T2KDataSource CustomMakeProjectTarget)
target_link_libraries(T2KDataSource ${MakeProjectLib} ${ROOT_LIBRARIES})

# Compile the main executable
add_executable(examplet2kdataframe src/example.cxx)
add_dependencies(examplet2kdataframe CustomMakeProjectTarget)
target_link_libraries(examplet2kdataframe T2KDataSource ${MakeProjectLib} ${ROOT_LIBRARIES} )

# Copy results to installation directory
install (FILES ${MakeProjectLib} DESTINATION lib)
install (TARGETS examplet2kdataframe DESTINATION bin)
install (FILES T2KDataSource.hxx DESTINATION include)
