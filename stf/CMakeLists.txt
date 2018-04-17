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
file(STRINGS ${CMAKE_SOURCE_DIR}/VERSION VR_VERSION)
string(REPLACE "." ";" VERSION_LIST ${VR_VERSION})
list(GET VERSION_LIST 0 VELOCIRAPTOR_VERSION_MAJOR)
list(GET VERSION_LIST 1 VELOCIRAPTOR_VERSION_MINOR)

# We have only C++ sources, but until cmake 3.4 the check_symbol_exists
# macro didn't work properly if only defining C++ as the project language
if( ${CMAKE_MAJOR_VERSION} EQUAL 3 AND ${CMAKE_MINOR_VERSION} LESS 4 )
	set(VR_LANGS C CXX)
else()
	set(VR_LANGS CXX)
endif()

project(VELOCIraptor VERSION ${VR_VERSION} LANGUAGES ${VR_LANGS})
set(CMAKE_CXX_STANDARD 11)

#
# Options users can give on the command line via -D
#
macro(vr_option optname optdesc)
	option(VR_${optname} "${optdec}" ON)
endmacro()

# Optional dependencies
option(VR_NO_MPI    "Don't attempt to include MPI support in VELOCIraptor" OFF)
option(VR_NO_OPENMP "Don't attempt to include OpenMP support in VELOCIraptor" OFF)
option(VR_NO_HDF5   "Don't attempt to include HDF5 support in VELOCIraptor" OFF)
option(VR_NO_XDR    "Don't include XDF (nchilada) input support (which in turn requires HDF5 support)" OFF)

# MPI options
vr_option(MPI_REDUCE        "Reduce impact of MPI memory overhead at the cost of extra cpu cycles. Suggested this be turned on")
vr_option(LARGE_MPI_DOMAIN  "Use large integers to represent MPI domains")

# Gadget handling details
vr_option(GADGET_LONGID  "Support for long IDSs in Gadget")
vr_option(GADGET_DPOS    "Support for double precision pos and vel in Gadget")
vr_option(GADGET_SMASS   "Support for single precision mass in Gadget")
vr_option(GADGET_HEAD2   "Support for header type 2 in Gadget")

# Gadget extra info
vr_option(GADGET_SPH_INFO   "Support for extra SPH information in Gadget")
vr_option(GADGET_STAR_INFO  "Support for extra star information in Gadget")
vr_option(GADGET_BH_INFO    "Support for extra black hole information in Gadget")

# Others
vr_option(USE_BH        "Use black holes")
vr_option(USE_HYDRO     "NBODY_USE_GAS + NBODY_USE_STARS + VR_USE_BH")
vr_option(STRUCTURE_DEN "Calculate local density dist. only for particles in field objects; uses all particles to estimate quantity")
vr_option(HALO_DEN      "Like STRUCTURE_DEN, but se particles inside field objects")
vr_option(ZOOM_SIM      "Used to run against simulations with a high resolution region")

# Let's true to our word
if (VR_USE_HYDRO)
   set(NBODY_USE_GAS ON)
	set(NBODY_USE_STARS ON)
	set(VR_USE_BH ON)
endif()


set(VR_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(VR_DEFINES "")
set(VR_LIBS "")
set(VR_CXX_FLAGS "")
set(VR_LINK_FLAGS "")

# Define macros depending on the user's choices
macro(vr_option_defines varname define)
	if (VR_${varname})
		list(APPEND VR_DEFINES ${define})
		set(VR_HAS_${varname} Yes)
	else()
		set(VR_HAS_${varname} No)
	endif()
endmacro()

vr_option_defines(MPI_REDUCE                MPIREDUCEMEM)
vr_option_defines(LARGE_MPI_DOMAIN          HUGEMPI)

vr_option_defines(GADGET_LONGID             GADGETLONGID)
vr_option_defines(GADGET_DPOS               GADGETDOUBLEPRECISION)
vr_option_defines(GADGET_SMASS              GADGETSINGLEMASSPRECISION)
vr_option_defines(GADGET_HEAD2              GADGET2FORMAT)
vr_option_defines(GADGET_SPH_INFO           EXTRASPHINFO)
vr_option_defines(GADGET_STAR_INFO          EXTRASTARINFO)
vr_option_defines(GADGET_BH_INFO            EXTRABHINFO)

vr_option_defines(USE_BH                    BHON)
vr_option_defines(STRUCTURE_DEN             STRUCDEN)
vr_option_defines(HALO_DEN                  HALOONLYDEN)
vr_option_defines(ZOOM_SIM                  HIGHRES)

#
# How we find GSL and set it up
#
macro(find_gsl)
	# Easy...
	find_package(GSL REQUIRED)
	list(APPEND VR_INCLUDE_DIRS ${GSL_INCLUDE_DIRS})
	list(APPEND VR_LIBS ${GSL_LIBRARIES})
endmacro()


#
# How we find HDF5 and set it up
#
macro(find_hdf5)
	find_package(HDF5 COMPONENTS CXX)
	if (HDF5_FOUND)
		list(APPEND VR_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
		list(APPEND VR_LIBS ${HDF5_LIBRARIES})
		list(APPEND VR_DEFINES USEHDF)
		set(VR_HAS_HDF5 Yes)
		if (NOT VR_NO_XDR)
			list(APPEND VR_DEFINES USEXDR)
			set(VR_HAS_XDR Yes)
		endif()
	endif()
endmacro()


#
# How we find OpenMP and set it up
#
macro(find_openmp)
	find_package(OpenMP COMPONENTS CXX)
	if (OPENMP_FOUND)
		list(APPEND VR_CXX_FLAGS ${OpenMP_CXX_FLAGS})
		list(APPEND VR_LINK_FLAGS ${OpenMP_CXX_FLAGS})
		list(APPEND VR_DEFINES USEOPENMP)
		set(VR_HAS_OPENMP Yes)
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


#
# Go and find our dependencies, depending on whether
# we actually need them or not
#
find_gsl()

set(VR_HAS_HDF5 No)
set(VR_HAS_XDR No)
if (NOT VR_NO_HDF5)
	find_hdf5()
endif()

set(VR_HAS_OPENMP No)
if (NOT VR_NO_OPENMP)
	find_openmp()
endif()

set(VR_HAS_MPI No)
if (NOT VR_NO_MPI)
	find_mpi()
endif()


# This provides us with the nbodylib library
# We need to add it unless it was already added by somebody else
if (NOT TARGET nbodylib)
	add_subdirectory(NBodylib)
endif()
list(APPEND VR_INCLUDE_DIRS ${NBODYLIB_INCLUDE_DIRS})
list(APPEND VR_DEFINES "${NBODYLIB_DEFINES}")


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
vr_report(Dependencies
          "MPI" MPI "OpenMP" OPENMP "HDF5" HDF5 "nchilada" XDR)
vr_report("MPI-specifics (ignore if no MPI)"
          "Reduce MPI memory overhead at the cost of extra CPU cycles" MPI_REDUCE
          "Use huge MPI domains" LARGE_MPI_DOMAIN)
vr_report(Gadget
          "Use longs IDs" GADGET_LONGID "Use double precision pos and vel" GADGET_DPOS
          "Use single precision mass" GADGET_SMASS "Use header type 2" GADGET_HEAD2
          "Use extra SPH information" GADGET_SPH_INFO "Use extra star information" GADGET_STAR_INFO
          "Use extra black hole information" GADGET_BH_INFO)
vr_report(Others
          "Activate black hole physics" USE_BH
          "Calculate local density dist. only for particles in field objects" STRUCTURE_DEN
          "Like above, but use particles inside field objects only for calclation" HALO_DEN
          "Used to run against simulations with a high resolution region" ZOOM_SIM)
message("")
message("Compilation")
message("-----------")
message("")
message(" Include dirs: ${VR_INCLUDE_DIRS}")
message(" Macros defined: ${VR_DEFINES}")
message(" Libs: ${VR_LIBS}")
message(" C++ flags: ${VR_CXX_FLAGS}")
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
if (_export)
	set(VR_INCLUDE_DIRS "${VR_INCLUDE_DIRS}" PARENT_SCOPE)
	set(VR_DEFINES "${VR_DEFINES}" PARENT_SCOPE)
	set(VR_CXX_FLAGS "${VR_CXX_FLAGS}" PARENT_SCOPE)
	set(VR_LINK_FLAGS "${VR_LINK_FLAGS}" PARENT_SCOPE)
	set(VR_LIBS "${VR_LIBS}" PARENT_SCOPE)
endif()
