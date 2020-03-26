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

cmake_minimum_required(VERSION 3.0)

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
vr_option(SINGLE_PRECISION  "Use single point precision to store all properties and perform all calculations" OFF)
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
if (VR_SINGLE_PRECISION)
	set(NBODY_SINGLE_PRECISION ON)
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

vr_option_defines(SINGLE_PRECISION          SINGLEPRECISION)
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

#
# How we find GSL and set it up
#
macro(find_gsl)
	find_package(GSL REQUIRED)
	list(APPEND VR_INCLUDE_DIRS ${GSL_INCLUDE_DIRS})
	list(APPEND VR_LIBS ${GSL_LIBRARIES})
endmacro()

#
# How we find HDF5 and set it up
#

macro(find_hdf5)
	# FindHDF5 needs an environment variable, oddly, unlike
	# most other packages that use normal cmake variables
	if (HDF5_ROOT)
		set(ENV{HDF5_ROOT} ${HDF5_ROOT})
	endif()
	find_package(HDF5 COMPONENTS C)
	if (HDF5_FOUND)
		#
		list(APPEND VR_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
		list(APPEND VR_LIBS ${HDF5_LIBRARIES})
		list(APPEND VR_DEFINES USEHDF)
		set(VR_HAS_HDF5 Yes)
		#check if parallel hdf present
		if (HDF5_IS_PARALLEL AND VR_HAS_MPI AND VR_ALLOWPARALLELHDF5)
			set (ENV{HDF5_PREFER_PARALLEL} true)
			set(VR_HAS_PARALLEL_HDF5 Yes)
			list(APPEND VR_DEFINES USEPARALLELHDF)
			if (HDF5_VERSION VERSION_GREATER "1.10.0" AND VR_ALLOWCOMPRESSIONPARALLELHDF5)
				set(VR_HAS_COMPRESSED_HDF5 Yes)
				list(APPEND VR_DEFINES USEHDFCOMPRESSION)
				list(APPEND VR_DEFINES PARALLELCOMPRESSIONACTIVE)
			endif()
		else()
			if (VR_ALLOWCOMPRESSIONHDF5)
				set(VR_HAS_COMPRESSED_HDF5 Yes)
				list(APPEND VR_DEFINES USEHDFCOMPRESSION)
			endif()
		endif()
    endif()
endmacro()


#
# How we find MPI and set it up
#
macro(find_mpi)
	find_package(MPI)
	if (MPI_FOUND)
		list(APPEND VR_INCLUDE_DIRS ${MPI_CXX_INCLUDE_PATH})
		list(APPEND VR_LIBS ${MPI_CXX_LIBRARIES})
		list(APPEND VR_CXX_FLAGS ${MPI_CXX_FLAGS})
		list(APPEND VR_LINK_FLAGS ${MPI_CXX_FLAGS})
		list(APPEND VR_DEFINES USEMPI)
		set(VR_HAS_MPI Yes)
	endif()
endmacro()
if (VR_USE_SWIFT_INTERFACE)
        list(APPEND VR_INCLUDE_DIRS ${SWIFT_INCLUDE_DIR})
endif()


#
# Make sure we have the git submodules we need
#
macro(ensure_git_submodules)
	if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/NBodylib/CMakeLists.txt" OR NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tools/velociraptor_python_tools.py")
		find_package(Git QUIET)
		if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
			# Update submodules as needed
			message(STATUS "Updating NBodylib and tools submodule")
			execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init
			                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
			                RESULT_VARIABLE GIT_SUBMOD_RESULT)
			if(NOT GIT_SUBMOD_RESULT EQUAL "0")
			    message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
			endif()
		else()
			message(FATAL_ERROR "Cannot get NBodylib submodule or tools submodule automatically.
				  Make sure you get your submodules")
		endif()
	endif()
endmacro()


#
# Go and find our dependencies, depending on whether
# we actually need them or not
#
ensure_git_submodules()
find_gsl()

set(VR_HAS_MPI No)
if (VR_MPI)
	find_mpi()
endif()

set(VR_HAS_HDF5 No)
set(VR_HAS_COMPRESSED_HDF5 No)
set(VR_HAS_PARALLEL_HDF5 No)
if (VR_HDF5)
	find_hdf5()
endif()



# This provides us with the nbodylib library
# We need to add it unless it was already added by somebody else
if (NOT TARGET nbodylib)
	add_subdirectory(NBodylib)
	if (NBODYLIB_VERSION VERSION_LESS "1.28")
		message(FATAL_ERROR "NBodyLib version ${NBODYLIB_VERSION} unsupported,
		VELOCIraptor requires >= 1.28, try running git submodule update --recursive --remote")
	endif()
	list(INSERT VR_DOC_DIRS 0 ${NBODYLIB_DOC_DIRS})
endif()
list(APPEND VR_INCLUDE_DIRS ${NBODYLIB_INCLUDE_DIRS})
list(APPEND VR_DEFINES "${NBODYLIB_DEFINES}")
list(APPEND VR_CXX_FLAGS "${NBODYLIB_CXX_FLAGS}")
list(APPEND VR_LINK_FLAGS "${NBODYLIB_LINK_FLAGS}")
list(APPEND VR_LIBS "${NBODYLIB_LIBS}")


#
# Tell the world what what we are doing
#
macro(vr_report feature)

	# Output feature name and underscore it in the next line
	message("\n${feature}")
	string(REGEX REPLACE "." "-" _underscores ${feature})
	message("${_underscores}\n")

	set(_args "${ARGN}")
	list(LENGTH _args _nargs)
	math(EXPR _nargs "${_nargs} - 1")
	foreach(_idx RANGE 0 ${_nargs} 2)

		# Items in the list come with a message first, then the variable name
		list(GET _args ${_idx} _msg)
		math(EXPR _idx2 "${_idx} + 1")
		list(GET _args ${_idx2} _varname)

		# We try to keep things up to 80 cols
		string(LENGTH ${_msg} _len)
		math(EXPR _nspaces "75 - ${_len}")
		string(RANDOM LENGTH ${_nspaces} _spaces)
		string(REGEX REPLACE "." " " _spaces "${_spaces}")
		string(CONCAT _msg "${_msg}" ${_spaces})
		message(" ${_msg} ${VR_HAS_${_varname}}")
	endforeach()
endmacro()

message("\nVELOCIraptor successfully configured with the following settings:")
vr_report("File formats"
          "HDF5" HDF5
		  "Compressed HDF5" COMPRESSED_HDF5
		  "Parallel HDF5" PARALLEL_HDF5
          "nchilada" XDR)
if (VR_HAS_COMPRESSED_HDF5  AND VR_HAS_PARALLEL_HDF5)
	message("\n WARNING: Parallel Compression HDF5 active, use with caution as it is unstable!\n")
endif()
vr_report("Precision-specifics"
        "Long Integers" LONG_INT
        "Single precision floats" SINGLE_PRECISION)
vr_report("OpenMP-specifics"
        "OpenMP support" OPENMP)
vr_report("MPI-specifics"
          "MPI support" MPI
          "Reduce MPI memory overhead at the cost of extra CPU cycles" MPI_REDUCE
          "Use huge MPI domains" LARGE_MPI_DOMAIN)
vr_report("Gadget"
          "Use longs IDs" GADGET_LONGID "Use double precision pos and vel" GADGET_DPOS
          "Use single precision mass" GADGET_SMASS "Use header type 2" GADGET_HEAD2
          "Use extra SPH information" GADGET_SPH_INFO "Use extra star information" GADGET_STAR_INFO
          "Use extra black hole information" GADGET_BH_INFO)
vr_report("Particle-specifics"
          "Activate gas (& associated physics, properties calculated)" USE_GAS
		  "Activate stars (& associated physics, properties calculated)" USE_STAR
		  "Activate black holes (& associated physics, properties calculated)" USE_BH
		  "Activate extra dark matter properties (& associated properties)" USE_EXTRA_DM_PROPERTIES
		  "Mass not stored (for uniform N-Body sims, reduce mem footprint)" NO_MASS
		  "Large memory KDTree to handle > max 32-bit integer entries per tree" USE_LARGE_KDTREE
	  )
vr_report("Simulation-specifics"
          "Used to run against simulations with a high resolution region" ZOOM_SIM
	      "Build library for integration into SWIFT Sim code " USE_SWIFT_INTERFACE)
vr_report("Others"
          "Calculate local density dist. only for particles in field objects" STRUCTURE_DEN
	      "Like above, but use particles inside field objects only for calclation" HALO_DEN)

message("")
message("Compilation")
message("-----------")
message("")
message(" Include dirs: ${VR_INCLUDE_DIRS}")
message(" Macros defined: ${VR_DEFINES}")
message(" Libs: ${VR_LIBS}")
message(" C++ flags: ${VR_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
message(" Link flags: ${VR_LINK_FLAGS}")
message("")

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
