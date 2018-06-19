# CMakeLists.txt to build the NBodylib library
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

set(_subdirs Analysis Cosmology InitCond KDTree Math NBody)

# We have the version here
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION NBODY_VERSION)
string(REPLACE "." ";" VERSION_LIST ${NBODY_VERSION})
list(GET VERSION_LIST 0 NBODY_VERSION_MAJOR)
list(GET VERSION_LIST 1 NBODY_VERSION_MINOR)

# We have only C++ sources, but until cmake 3.4 the check_symbol_exists
# macro didn't work properly if only defining C++ as the project language
if( ${CMAKE_MAJOR_VERSION} EQUAL 3 AND ${CMAKE_MINOR_VERSION} LESS 4 )
	set(NBODY_LANGS C CXX)
else()
	set(NBODY_LANGS CXX)
endif()

project(NBodyLib VERSION ${NBODY_VERSION} LANGUAGES ${NBODY_LANGS})
set(CMAKE_CXX_STANDARD 11)


# Function to add a "doc" target, which will doxygen out the given
# list of directories
function(try_add_doc_target doc_dirs)
	find_program(DOXYGEN_FOUND doxygen)
	if (NOT DOXYGEN_FOUND)
		return()
	endif()

	# Create a target for each individual doc directory, then a final one
	# that depends on them all
	message("-- Adding doc target for directories: ${doc_dirs}")
	set(_dependencies "")
	set(x 1)
	foreach(_doc_dir IN ITEMS ${doc_dirs})
		add_custom_command(OUTPUT ${_doc_dir}/xml/index.xml
		                   COMMAND doxygen
		                   WORKING_DIRECTORY ${_doc_dir})
		add_custom_target(doc_${x}
		                  COMMAND doxygen
		                  WORKING_DIRECTORY ${_doc_dir})
		list(APPEND _dependencies doc_${x})
		math(EXPR x "${x} + 1")
	endforeach()
	add_custom_target(doc DEPENDS "${_dependencies}")
endfunction()

# Are we part of an outside build?
# If so, we want to let the caller know what to include, define, etc at the end
set(_export OFF)
get_directory_property(_hasParent PARENT_DIRECTORY)
if (_hasParent)
	set(_export ON)
endif()

# Options users can give on the command line via -D
macro(nbody_option optname optdesc status )
	option(NBODY_${optname} "${optdesc}" "${status}")
endmacro()

# Optional dependencies
nbody_option(OPENMP "Attempt to include OpenMP support in NBodyLib" ON)

# Precision options
nbody_option(SINGLE_PRECISION  "Use single point precision to store all properties and perform all calculations" OFF)
nbody_option(LONG_INT          "Use long ints to represent all integers. Needed if dealing with more than MAXINT number of particles" OFF)

# Particle class details options
nbody_option(NO_MASS                    "Do not store the mass as all particles are the same mass" OFF)
nbody_option(SINGLE_PARTICLE_PRECISION  "Do not store the mass as all particles are the same mass" OFF)
nbody_option(UNSIGNED_PARTICLE_PIDS     "Use unsigned particle PIDs" OFF)
nbody_option(UNSIGNED_PARTICLE_IDS      "Use unsigned particle IDs" OFF)
nbody_option(USE_GAS                    "Particle class has gas" OFF)
nbody_option(USE_STARS                  "Particle class has stars" OFF)
nbody_option(USE_BARYONS                "VR_USE_GAS + VR_USE_STARS" OFF)
nbody_option(USE_HYDRO                  "VR_USE_GAS + VR_USE_STARS+VR_USE_BH" OFF)

# Process high-level options
if (NBODY_USE_BARYONS)
	set(NBODY_USE_GAS ON)
	set(NBODY_USE_STARS ON)
endif()
if (NBODY_USE_HYDRO)
	set(NBODY_USE_GAS ON)
	set(NBODY_USE_STARS ON)
	set(NBODY_USE_BH ON)
endif()


# Flags used for compiling nbodylib.
# They are optionally exported at the end if there is an outer scope
set(NBODY_INCLUDE_DIRS "")
set(NBODY_DEFINES "")
set(NBODY_LIBS "")
set(NBODY_CXX_FLAGS "")
set(NBODY_LINK_FLAGS "")
set(NBODY_DOC_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/doc)

# Define macros depending on the user's choices
macro(nbody_option_defines varname define)
	if (NBODY_${varname})
		list(APPEND NBODY_DEFINES ${define})
		set(NBODY_HAS_${varname} Yes)
	else()
		set(NBODY_HAS_${varname} No)
	endif()
endmacro()

nbody_option_defines(SINGLE_PRECISION          SINGLEPRECISION)
nbody_option_defines(LONG_INT                  LONGINT)
nbody_option_defines(NO_MASS                   NOMASS)
nbody_option_defines(SINGLE_PARTICLE_PRECISION LOWPRECISIONPOS)
nbody_option_defines(UNSIGNED_PARTICLE_PIDS    PARTICLEUPIDS)
nbody_option_defines(UNSIGNED_PARTICLE_IDS     PARTICLEUIDS)
nbody_option_defines(LARGEKDTREE               LARGETREE)
nbody_option_defines(USE_GAS                   GASON)
nbody_option_defines(USE_STARS                 STARON)
nbody_option_defines(USE_BH                    BHON)

#
# How we find GSL and set it up
#
macro(find_gsl)
	find_package(GSL REQUIRED)
	list(APPEND NBODY_INCLUDE_DIRS ${GSL_INCLUDE_DIRS})
	list(APPEND NBODY_LIBS ${GSL_LIBRARIES})
endmacro()

# How we find OpenMP and set it up
#
macro(find_openmp)
	find_package(OpenMP COMPONENTS CXX)
	if (OPENMP_FOUND)
		list(APPEND NBODY_CXX_FLAGS ${OpenMP_CXX_FLAGS})
		list(APPEND NBODY_LINK_FLAGS ${OpenMP_CXX_FLAGS})
		list(APPEND NBODY_DEFINES USEOPENMP USEOMP)
		set(NBODY_HAS_OPENMP Yes)
	endif()
endmacro()

#
# Now go and find our dependencies, depending on whether
# we actually need them or not
#
find_gsl()

set(NBODY_HAS_OPENMP No)
if (NBODY_OPENMP)
	find_openmp()
endif()


# Add subdirectories to list of include directories
foreach(_subdir IN ITEMS ${_subdirs})
	list(APPEND NBODY_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src/${_subdir}")
endforeach()

#
# Tell the world what what we are doing
#
macro(nbody_report feature)

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
		message(" ${_msg} ${NBODY_HAS_${_varname}}")
	endforeach()
endmacro()

message("\nNBodyLib successfully configured with the following settings:")
nbody_report(Dependencies "OpenMP" OPENMP)
nbody_report(Types "All calculations/properties stored as float" SINGLE_PRECISION "All integeres are long int" LONG_INT)
nbody_report("Particle data" "Do not store mass, all particles are the same mass" NO_MASS
             "Use single precision to store positions, velocities, other props" SINGLE_PARTICLE_PRECISION
             "Use unsigned particle PIDs" UNSIGNED_PARTICLE_PIDS "Use unsigned particle IDs" UNSIGNED_PARTICLE_IDS
             "Activate gas" USE_GAS "Activate stars" USE_STARS)
message("")
message("Compilation")
message("-----------")
message("")
message(" Include directories: ${NBODY_INCLUDE_DIRS}")
message(" Macros defined: ${NBODY_DEFINES}")
message(" Libs: ${NBODY_LIBS}")
message(" C++ flags: ${NBODY_CXX_FLAGS}")
message(" Link flags: ${NBODY_LINK_FLAGS}")
message("")

#
# Macro used by subdirectories to build each a static library
#
macro(add_nbody_lib libname sources)
	add_library(${libname} OBJECT ${sources})
	target_compile_definitions(${libname} PRIVATE ${NBODY_DEFINES})
	target_include_directories(${libname} PRIVATE ${NBODY_INCLUDE_DIRS})
	if (NBODY_CXX_FLAGS)
		set_target_properties(${libname} PROPERTIES COMPILE_FLAGS ${NBODY_CXX_FLAGS})
	endif()
endmacro()

# Include all subdirectories,
# where the magic happens
foreach(_subdir IN ITEMS ${_subdirs})
	add_subdirectory(src/${_subdir})
endforeach()

# Build the final static libnbodylib.a library
add_library(nbodylib STATIC
    $<TARGET_OBJECTS:Analysis> $<TARGET_OBJECTS:Cosmology>
    $<TARGET_OBJECTS:InitCond> $<TARGET_OBJECTS:KDTree>
    $<TARGET_OBJECTS:Math> $<TARGET_OBJECTS:NBody>
)

# Export the include directories and definitions, if necessary
# If building on our own, add the "doc" target
if (_export)
	set(NBODYLIB_VERSION "${NBODY_VERSION_MAJOR}.${NBODY_VERSION_MINOR}" PARENT_SCOPE)
	set(NBODYLIB_VERSION_MAJOR "${NBODY_VERSION_MAJOR}" PARENT_SCOPE)
	set(NBODYLIB_VERSION_MINOR "${NBODY_VERSION_MINOR}" PARENT_SCOPE)
	set(NBODYLIB_INCLUDE_DIRS "${NBODY_INCLUDE_DIRS}" PARENT_SCOPE)
	set(NBODYLIB_DEFINES "${NBODY_DEFINES}" PARENT_SCOPE)
	set(NBODYLIB_CXX_FLAGS "${NBODY_CXX_FLAGS}" PARENT_SCOPE)
	set(NBODYLIB_LINK_FLAGS "${NBODY_LINK_FLAGS}" PARENT_SCOPE)
	set(NBODYLIB_LIBS "${NBODY_LIBS}" PARENT_SCOPE)
	set(NBODYLIB_DOC_DIRS "${NBODY_DOC_DIRS}" PARENT_SCOPE)
else()
	try_add_doc_target("${NBODY_DOC_DIRS}")
endif()
