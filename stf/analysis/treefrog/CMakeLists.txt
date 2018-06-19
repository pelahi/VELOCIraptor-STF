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

set(TF_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}")
set(TF_DEFINES "")
set(TF_LIBS "")
set(TF_CXX_FLAGS "")
set(TF_LINK_FLAGS "")

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
list(APPEND TF_INCLUDE_DIRS ${NBODYLIB_INCLUDE_DIRS})
list(APPEND TF_DEFINES ${NBODYLIB_DEFINES})
list(APPEND TF_LIBS ${NBODYLIB_LIBS})
list(APPEND TF_CXX_FLAGS ${NBODYLIB_CXX_FLAGS})
list(APPEND TF_LINK_FLAGS ${NBODYLIB_LINK_FLAGS})
list(APPEND TF_DOC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/doc")

# Optionally include VELOCIraptor, which might be needed (and internally
# includes nbodylib)
if (TF_BUILD_VELOCIRAPTOR)
	add_subdirectory(../.. velociraptor)
	list(APPEND TF_INCLUDE_DIRS ${VR_INCLUDE_DIRS})
	list(APPEND TF_DEFINES ${VR_DEFINES})
	list(APPEND TF_CXX_FLAGS ${VR_CXX_FLAGS})
	list(APPEND TF_LINK_FLAGS ${VR_LINK_FLAGS})
	list(APPEND TF_LIBS ${VR_LIBS})
	list(INSERT TF_DOC_DIRS 0 ${VR_DOC_DIRS})
endif()

list(INSERT TF_DOC_DIRS 0 ${NBODYLIB_DOC_DIRS})

# cleanup duplicates and convert to single string
list(REMOVE_DUPLICATES TF_CXX_FLAGS)
string(REPLACE ";" " " TF_CXX_FLAGS "${TF_CXX_FLAGS}")
list(REMOVE_DUPLICATES TF_LINK_FLAGS)
string(REPLACE ";" " " TF_LINK_FLAGS "${TF_LINK_FLAGS}")

include_directories(${TF_INCLUDE_DIRS})
add_executable(treefrog ${_treefrog_sources})
target_link_libraries(treefrog nbodylib ${TF_LIBS})
target_compile_definitions(treefrog PRIVATE ${TF_DEFINES})
if (TF_CXX_FLAGS)
	set_target_properties(treefrog PROPERTIES COMPILE_FLAGS ${TF_CXX_FLAGS})
endif()
if (TF_LINK_FLAGS)
	set_target_properties(treefrog PROPERTIES LINK_FLAGS ${TF_LINK_FLAGS})
endif()

try_add_doc_target("${TF_DOC_DIRS}")