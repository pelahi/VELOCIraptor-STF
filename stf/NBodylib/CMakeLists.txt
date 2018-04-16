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

# Are we part of an outside build?
# If so, we want to let the caller know what to include
set(_export_inc_dirs OFF)
get_directory_property(_hasParent PARENT_DIRECTORY)
if (_hasParent)
	set(_export_inc_dirs ON)
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
	target_compile_definitions(${libname} PRIVATE ${VR_NBODY_DEFINES})
	target_include_directories(${libname} PRIVATE ${_inc_dirs})
	if (VR_NBODY_CXX_FLAGS)
		set_target_properties(${libname} PROPERTIES COMPILE_FLAGS ${VR_NBODY_CXX_FLAGS})
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

# Export the include directories, if necessary
if (_export_inc_dirs)
	set(NBODYLIB_INCLUDE_DIRS "${_inc_dirs}" PARENT_SCOPE)
endif()
