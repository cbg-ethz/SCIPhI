# This file is based on the CMakeLists.txt from the examples of dlib

cmake_minimum_required(VERSION 3.0.0)

project(SCIPhI)

option(BUILD_DIAGNOSTICS "Build diagnostic tools" OFF)

if (BUILD_DIAGNOSTICS)
	set(CMAKE_CXX_STANDARD 14)

	# Load the SeqAn module and fail if not found.
	find_package (SeqAn REQUIRED)
else()
	set(CMAKE_CXX_STANDARD 11)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for boost
find_package(Boost 1.47 COMPONENTS program_options graph filesystem iostreams unit_test_framework REQUIRED)

# Search for zlib as a dependency for SeqAn.
#find_package (ZLIB)

# Tell cmake we will need dlib.  This command will pull in dlib and compile it
# into your project.  Note that you don't need to compile or install dlib.  All
# cmake needs is the dlib source code folder and it will take care of everything.
add_subdirectory(${DLIB_INCLUDEDIR} dlib)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${Boost_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
message(${CMAKE_CURRENT_SOURCE_DIR})

# Add link directories
link_directories(${Boost_LIBRARY_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add executable and link against SeqAn dependencies.
add_executable(sciphi src/findBestTrees.cpp)
target_link_libraries (sciphi dlib::dlib ${Boost_LIBRARIES})

if (BUILD_DIAGNOSTICS)
	add_executable(convertTree2Map src/convertTreeIntoClusterMap.cpp)
	target_link_libraries (convertTree2Map dlib::dlib ${SEQAN_LIBRARIES} ${Boost_LIBRARIES})
	target_include_directories(convertTree2Map PRIVATE ${SEQAN_INCLUDE_DIRS})

	add_executable(createSingleCellTree src/createSingleCellSeqTree.cpp)
	target_link_libraries (createSingleCellTree dlib::dlib ${SEQAN_LIBRARIES} ${Boost_LIBRARIES})
	target_include_directories(createSingleCellTree PRIVATE ${SEQAN_INCLUDE_DIRS})

	add_executable(sciphi_test test/sciphi_test.cpp)
	target_link_libraries (sciphi_test dlib::dlib ${SEQAN_LIBRARIES} ${Boost_LIBRARIES})
	target_include_directories(sciphi_test PRIVATE ${SEQAN_INCLUDE_DIRS})
endif()
# To compile this program all you need to do is ask cmake.  You would type
# these commands from within the directory containing this CMakeLists.txt
# file:
#   mkdir build
#   cd build
#   cmake ..
#   cmake --build . --config Release
#
