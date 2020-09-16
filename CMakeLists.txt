# CMakeLists.txt
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2018
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Rodrigo Tobar
#
# This file is part of VELOCIraptor.

cmake_minimum_required(VERSION 3.13)

# We have the version here
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION VR_VERSION)
string(REPLACE "." ";" VERSION_LIST ${VR_VERSION})
list(GET VERSION_LIST 0 VR_VERSION_MAJOR)
list(GET VERSION_LIST 1 VR_VERSION_MINOR)

# We have only C++ sources, but until cmake 3.4 the check_symbol_exists
# macro didn't work properly if only defining C++ as the project language
if( ${CMAKE_MAJOR_VERSION} EQUAL 3 AND ${CMAKE_MINOR_VERSION} LESS 4 )
	set(VR_LANGS C CXX)
else()
	set(VR_LANGS CXX C)
endif()

project(VELOCIraptor VERSION ${VR_VERSION} LANGUAGES ${VR_LANGS})
set(CMAKE_CXX_STANDARD 11)
# Set the code to always produce position independent code by default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#include cmake stuff
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
# include some general cmake utility functions
include (VRCMakeUtils)
# complation summary
include(VRCompilationMessages)

#
# Options users can give on the command line via -D
#
macro(vr_option optname optdesc status)
	option(VR_${optname} "${optdesc}" "${status}")
endmacro()

# Input formats
vr_option(HDF5 "Attempt to include HDF5 support in VELOCIraptor" ON)
vr_option(ALLOWCOMPRESSIONHDF5 "Attempt to include HDF5 compression support in VELOCIraptor" ON)
vr_option(ALLOWPARALLELHDF5 "Attempt to include parallel HDF5 support in VELOCIraptor" ON)
vr_option(ALLOWCOMPRESSIONPARALLELHDF5 "Attempt to include parallel HDF5 compression support in VELOCIraptor" OFF)
vr_option(XDR  "XDR input support (used by nchilada)" OFF)

# Precision options
vr_option(LONG_INT          "Use long ints to represent all integers. Needed if dealing with more than MAXINT number of particles" ON)

# OpenMP options
vr_option(OPENMP           "Attempt to include OpenMP support in VELOCIraptor" ON)

# MPI options
vr_option(MPI               "Attempt to include MPI support in VELOCIraptor" ON)
vr_option(MPI_REDUCE        "Reduce impact of MPI memory overhead at the cost of extra cpu cycles. Suggested this be turned on" ON)
vr_option(LARGE_MPI_DOMAIN  "Use large integers to represent MPI domains" OFF)

# Gadget handling details
vr_option(GADGET_LONGID  "Support for long IDSs in Gadget" OFF)
vr_option(GADGET_DPOS    "Support for double precision pos and vel in Gadget" OFF)
vr_option(GADGET_SMASS   "Support for single precision mass in Gadget" OFF)
vr_option(GADGET_HEAD2   "Support for header type 2 in Gadget" OFF)

# Gadget extra info
vr_option(GADGET_SPH_INFO   "Support for extra SPH information in Gadget" OFF)
vr_option(GADGET_STAR_INFO  "Support for extra star information in Gadget" OFF)
vr_option(GADGET_BH_INFO    "Support for extra black hole information in Gadget" OFF)

# Particle data options
vr_option(USE_STAR      "Use star particles" OFF)
vr_option(USE_GAS       "Use gas particles" OFF)
vr_option(USE_BH        "Use black holes" OFF)
vr_option(USE_EXTRA_DM_PROPERTIES "Store extra dark matter properties" OFF)
vr_option(USE_HYDRO     "Use all particle types, (gas, star, bh, etc)" OFF)
vr_option(NO_MASS       "Particles do not store mass (useful for pure N-body sims and reducing memory footprint)" OFF)
vr_option(USE_EXTRA_INPUT_INFO     "Store where particles are located in the input file" OFF)
vr_option(USE_EXTRA_FOF_INFO     "Store particles fof, group, SO ids" OFF)
vr_option(USE_LARGE_KDTREE      "Require large mem KDTree as there are more than > max 32-bit integer entries" OFF)

# VELOCiraptor behavior options
vr_option(STRUCTURE_DEN "Calculate local density dist. only for particles in field objects; uses all particles to estimate quantity" ON)
vr_option(HALO_DEN      "Like STRUCTURE_DEN, but se particles inside field objects" OFF)
vr_option(ZOOM_SIM      "Used to run against simulations with a high resolution region" OFF)
vr_option(USE_SWIFT_INTERFACE "Used to compile the code with hooks for being called from the SWIFT Hydro nbody code" OFF)


# Let's true to our word
if (VR_USE_SWIFT_INTERFACE)
	set(NBODY_USE_SWIFT_INTERFACE ON)
endif()
if (VR_USE_HYDRO)
	set(NBODY_USE_GAS ON)
	set(NBODY_USE_STARS ON)
	set(NBODY_USE_BH ON)
	set(VR_USE_BH ON)
	set(VR_USE_GAS ON)
	set(VR_USE_STAR ON)
endif()
if (VR_USE_STAR)
	set(NBODY_USE_STARS ON)
endif()
if (VR_USE_GAS)
	set(NBODY_USE_GAS ON)
endif()
if (VR_USE_BH)
	set(NBODY_USE_BH ON)
endif()
if (VR_USE_EXTRA_DM_PROPERTIES)
	set(NBODY_USE_EXTRA_DM_PROPERTIES ON)
endif()
if (VR_USE_EXTRA_INPUT_INFO)
	set(NBODY_USE_EXTRA_INPUT_INFO ON)
endif()
if (VR_USE_EXTRA_FOF_INFO)
	set(NBODY_USE_EXTRA_FOF_INFO ON)
endif()
if (VR_NO_MASS)
	set(NBODY_NO_MASS ON)
endif()
if (VR_LONG_INT)
	set(NBODY_LONG_INT ON)
endif()
if (VR_USE_LARGE_KDTREE)
	set(NBODY_USE_LARGE_KDTREE ON)
endif()
set(NBODY_OPENMP ${VR_OPENMP})


set(VR_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(VR_DEFINES "")
set(VR_LIBS "")
set(VR_CXX_FLAGS "")
set(VR_LINK_FLAGS "")
set(VR_DOC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/doc")

# Define macros depending on the user's choices
macro(vr_option_defines varname define)
	if (VR_${varname})
		list(APPEND VR_DEFINES ${define})
		set(VR_HAS_${varname} Yes)
	else()
		set(VR_HAS_${varname} No)
	endif()
endmacro()

vr_option_defines(XDR                       USEXDR)

vr_option_defines(LONG_INT                  LONGINT)

vr_option_defines(OPENMP                    USEOPENMP)

vr_option_defines(MPI_REDUCE                MPIREDUCEMEM)
vr_option_defines(LARGE_MPI_DOMAIN          HUGEMPI)

vr_option_defines(GADGET_LONGID             GADGETLONGID)
vr_option_defines(GADGET_DPOS               GADGETDOUBLEPRECISION)
vr_option_defines(GADGET_SMASS              GADGETSINGLEMASSPRECISION)
vr_option_defines(GADGET_HEAD2              GADGET2FORMAT)
vr_option_defines(GADGET_SPH_INFO           EXTRASPHINFO)
vr_option_defines(GADGET_STAR_INFO          EXTRASTARINFO)
vr_option_defines(GADGET_BH_INFO            EXTRABHINFO)

vr_option_defines(USE_GAS                   GASON)
vr_option_defines(USE_STAR                  STARON)
vr_option_defines(USE_BH                    BHON)
vr_option_defines(USE_EXTRA_DM_PROPERTIES   EXTRADMON)
vr_option_defines(NO_MASS                   NOMASS)
vr_option_defines(USE_EXTRA_INPUT_INFO      EXTRAINPUTINFO)
vr_option_defines(USE_EXTRA_FOF_INFO        EXTRAFOFINFO)

vr_option_defines(STRUCTURE_DEN             STRUCDEN)
vr_option_defines(HALO_DEN                  HALOONLYDEN)
vr_option_defines(ZOOM_SIM                  HIGHRES)
vr_option_defines(USE_LARGE_KDTREE          LARGETREE)

vr_option_defines(USE_SWIFT_INTERFACE       SWIFTINTERFACE)

if (VR_USE_SWIFT_INTERFACE)
        list(APPEND VR_INCLUDE_DIRS ${SWIFT_INCLUDE_DIR})
endif()


#
# Add git revision information in a generated .cxx file, which we then include
# as part of the library
#
include(cmake/GetGitRevisionDescription.cmake)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_local_changes(GIT_HAS_LOCAL_CHANGES)
set(git_revision_cxx "${CMAKE_CURRENT_BINARY_DIR}/git_revision.cxx")
configure_file("src/git_revision.cxx.in" "${git_revision_cxx}" @ONLY)


#
# Go and find our dependencies, depending on whether
# we actually need them or not
#
vr_ensure_git_submodules()
find_gsl()
vr_mpi()
vr_hdf5()
vr_config_errors()
vr_nbodylib()

# compilation state
vr_compilation_summary()

# Are we part of an outside build (i.e., of a treefrog build)?
# If so, we want to let the caller know what to include
set(_export OFF)
get_directory_property(_hasParent PARENT_DIRECTORY)
if (_hasParent)
	set(_export ON)
endif()

include_directories(${VR_INCLUDE_DIRS})

# This provides us with the velociraptor library and the stf binary
add_subdirectory(src)

# Export the include directories, if necessary
# If building on our own, add the "doc" target
if (_export)
	set(VR_VERSION "${VR_VERSION_MAJOR}.${VR_VERSION_MINOR}" PARENT_SCOPE)
	set(VR_VERSION_MAJOR "${VR_VERSION_MAJOR}" PARENT_SCOPE)
	set(VR_VERSION_MINOR "${VR_VERSION_MINOR}" PARENT_SCOPE)
	set(VR_INCLUDE_DIRS "${VR_INCLUDE_DIRS}" PARENT_SCOPE)
	set(VR_DEFINES "${VR_DEFINES}" PARENT_SCOPE)
	set(VR_CXX_FLAGS "${VR_CXX_FLAGS}" PARENT_SCOPE)
	set(VR_LINK_FLAGS "${VR_LINK_FLAGS}" PARENT_SCOPE)
	set(VR_LIBS "${VR_LIBS}" PARENT_SCOPE)
	set(VR_DOC_DIRS "${VR_DOC_DIRS}" PARENT_SCOPE)
else()
	try_add_doc_target("${VR_DOC_DIRS}")
endif()
