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

# We have the version here
file(STRINGS ${CMAKE_SOURCE_DIR}/VERSION NBODY_VERSION)
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

# Options users can give on the command line via -D
macro(nbody_option optname optdesc)
	option(NBODY_${optname} "${optdec}" ON)
endmacro()

# Precision options
nbody_option(SINGLE_PRECISION  "Use single point precision to store all properties and perform all calculations")
nbody_option(LONG_INT          "Use long ints to represent all integers. Needed if dealing with more than MAXINT number of particles")

# Particle class details options
nbody_option(NO_MASS                    "Do not store the mass as all particles are the same mass")
nbody_option(SINGLE_PARTICLE_PRECISION  "Do not store the mass as all particles are the same mass")
nbody_option(UNSIGNED_PARTICLE_PIDS     "Use unsigned particle PIDs")
nbody_option(UNSIGNED_PARTICLE_IDS      "Use unsigned particle IDs")
nbody_option(USE_GAS                    "Particle class has gas")
nbody_option(USE_STARS                  "Particle class has stars")
nbody_option(USE_BARYONS                "VR_USE_GAS + VR_USE_STARS")

# Process high-level options
if (NBODY_USE_BARYONS)
	set(NBODY_USE_GAS ON)
	set(NBODY_USE_STARS ON)
endif()


# Flags used for compiling nbodylib.
# They are optionally exported at the end if there is an outer scope
set(NBODY_DEFINES "")
set(NBODY_LIBS "")
set(NBODY_CXX_FLAGS "")
set(NBODY_LINK_FLAGS "")

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
nbody_report(Types "All calculations/properties stored as float" SINGLE_PRECISION "All integeres are long int" LONG_INT)
nbody_report("Particle data" "Do not store mass, all particles are the same mass" NO_MASS
             "Use single precision to store positions, velocities, other props" SINGLE_PARTICLE_PRECISION
             "Use unsigned particle PIDs" UNSIGNED_PARTICLE_PIDS "Use unsigned particle IDs" UNSIGNED_PARTICLE_IDS
             "Activate gas" USE_GAS "Activate stars" USE_STARS)
message("")
message("Compilation")
message("-----------")
message("")
message(" Defines: ${NBODY_DEFINES}")
message(" Libs: ${NBODY_LIBS}")
message(" C++ flags: ${NBODY_CXX_FLAGS}")
message(" Link flags: ${NBODY_LINK_FLAGS}")
message("")


# Are we part of an outside build?
# If so, we want to let the caller know what to include, define, etc
set(_export OFF)
get_directory_property(_hasParent PARENT_DIRECTORY)
if (_hasParent)
	set(_export ON)
endif()

set(_subdirs Analysis Cosmology InitCond KDTree Math NBody)

# Calculate all subdirectories first
set(_inc_dirs "")
foreach(_subdir IN ITEMS ${_subdirs})
	list(APPEND _inc_dirs "${CMAKE_CURRENT_SOURCE_DIR}/src/${_subdir}")
endforeach()

#
# Macro used by subdirectories to build each a static library
#
macro(add_nbody_lib libname sources)
	add_library(${libname} OBJECT ${sources})
	target_compile_definitions(${libname} PRIVATE ${NBODY_DEFINES})
	target_include_directories(${libname} PRIVATE ${_inc_dirs})
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
if (_export)
	set(NBODYLIB_INCLUDE_DIRS "${_inc_dirs}" PARENT_SCOPE)
	set(NBODYLIB_DEFINES "${NBODY_DEFINES}" PARENT_SCOPE)
endif()
