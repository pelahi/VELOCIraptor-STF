# CMakeLists.txt for treefrog
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2018
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Rodrigo Tobar
#
# This file is part of VELOCIraptor.

cmake_minimum_required(VERSION 3.0)

# We have the version here
file(STRINGS ${CMAKE_SOURCE_DIR}/VERSION TF_VERSION)
string(REPLACE "." ";" VERSION_LIST ${TF_VERSION})
list(GET VERSION_LIST 0 TF_VERSION_MAJOR)
list(GET VERSION_LIST 1 TF_VERSION_MINOR)

# We have only C++ sources, but until cmake 3.4 the check_symbol_exists
# macro didn't work properly if only defining C++ as the project language
if( ${CMAKE_MAJOR_VERSION} EQUAL 3 AND ${CMAKE_MINOR_VERSION} LESS 4 )
	set(TF_LANGS C CXX)
else()
	set(TF_LANGS CXX)
endif()

project(TreeFrog VERSION ${TF_VERSION} LANGUAGES ${TF_LANGS})
set(CMAKE_CXX_STANDARD 11)

# Command-line treefrog-specific options
option(TF_BUILD_VELOCIRAPTOR "Whether VELOCiraptor should be built together with treefrog or not" OFF)

set(_treefrog_sources
    allvars.cxx
    crosscheck.cxx
    idroutines.cxx
    io.cxx
    main.cxx
    mpiroutines.cxx
    otherio.cxx
    stfio.cxx
    ui.cxx
    utilities.cxx)

# Include nbodylib, which is definitely needed
add_subdirectory(../../NBodylib NBodylib)
include_directories("${NBODYLIB_INCLUDE_DIRS}")

# Optionally include VELOCIraptor, which might be needed (and internally
# includes nbodylib)
if (TF_BUILD_VELOCIRAPTOR)
	add_subdirectory(../.. velociraptor)
	include_directories("${VR_INCLUDE_DIRS}")
endif()

add_executable(treefrog ${_treefrog_sources})
target_link_libraries(treefrog nbodylib)
