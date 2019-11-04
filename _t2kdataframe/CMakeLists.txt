# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(t2kdataframecpp)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Debug)

set( T2KDATAFRAME_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
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
    # message(FATAL_ERROR "ERROR: unable to find file to run TFile::MakeProject. Set this with -D MAKE_PROJECT_FILE_NAME=<filename.root>")
    set(MAKE_PROJECT_FILE_NAME "/storage/epp2/t2k/data/nd280/production006/B/mcp/genie/2010-11-air/magnet/run4/anal/oa_gn_beam_91400000-0002_fdbfaupp6ljn_anal_000_prod006magnet201011airc-bsdv01_2.root")
endif()
message(STATUS "using MAKE_PROJECT_FILE_NAME=${MAKE_PROJECT_FILE_NAME}")

# Run ROOT's TFile::MakeProject
set(MakeProjectObj ${T2KDATAFRAME_BINARY_DIR}/libT2KDataFrameMakeProject/libT2KDataFrameMakeProjectProjectSource.o)
add_custom_command(
        OUTPUT ${MakeProjectObj}
        COMMAND $ENV{ROOTSYS}/bin/root -l -q ${MAKE_PROJECT_FILE_NAME} ${T2KDATAFRAME_SOURCE_DIR}/t2kdataframe/buildmakeproject.C
        DEPENDS ${MAKE_PROJECT_FILE_NAME}
)
add_custom_target(CustomMakeProjectTarget DEPENDS ${MakeProjectObj})

# Ensure all directories are accessible
include_directories( ${T2KDATAFRAME_BINARY_DIR})
include_directories( ${T2KDATAFRAME_BINARY_DIR}/libT2KDataFrameMakeProject)
include_directories("${T2KDATAFRAME_SOURCE_DIR}")

# Generate library to use NCFGD code in python
ROOT_GENERATE_DICTIONARY(T2KDataFrameDict t2kdataframe/T2KDataSource.hxx OPTIONS -inlineInputHeader LINKDEF t2kdataframe/T2KDataSourceLinkDef.h DEPENDENCIES CustomMakeProjectTarget)
add_dependencies(T2KDataFrameDict CustomMakeProjectTarget)
add_library(T2KDataSource SHARED t2kdataframe/T2KDataSource.cxx T2KDataFrameDict.cxx)
add_dependencies(T2KDataSource CustomMakeProjectTarget)
add_dependencies(T2KDataSource T2KDataFrameDict)
target_link_libraries(T2KDataSource ${MakeProjectObj} ${ROOT_LIBRARIES})

# Compile the main executable
add_executable(examplet2kdataframe t2kdataframe/example.cxx)
add_dependencies(examplet2kdataframe CustomMakeProjectTarget)
target_link_libraries(examplet2kdataframe T2KDataSource ${ROOT_LIBRARIES} )

set(EXAMPLE_SCRIPTS example.py)

foreach(_script ${EXAMPLE_SCRIPTS})
    configure_file(
            ${T2KDATAFRAME_SOURCE_DIR}/${_script}
            ${T2KDATAFRAME_BINARY_DIR}/${_script}
            COPYONLY
    )
endforeach()

# Copy results to installation directory
install (TARGETS T2KDataSource DESTINATION pyt2k)
install (TARGETS examplet2kdataframe DESTINATION bin)
install (FILES t2kdataframe/T2KDataSource.hxx DESTINATION include/t2kdataframe)