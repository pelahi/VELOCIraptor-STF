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


#
# Macro used by subdirectories to build each a static library
#
macro(add_nbody_lib libname sources)
	add_library(${libname} STATIC ${sources})
	target_compile_definitions(${libname} PRIVATE ${VR_NBODY_DEFINES})
	if (VR_NBODY_CXX_FLAGS)
		set_target_properties(${libname} PROPERTIES COMPILE_FLAGS ${VR_NBODY_CXX_FLAGS})
	endif()
endmacro()

add_subdirectory(src/Analysis)
add_subdirectory(src/Cosmology)
add_subdirectory(src/InitCond)
add_subdirectory(src/KDTree)
add_subdirectory(src/Math)
add_subdirectory(src/NBody)


# Establish the dependencies between the different libraries
add_dependencies(Cosmology Math)
add_dependencies(NBody Math)
add_dependencies(InitCond Math NBody)
add_dependencies(KDTree Math NBody)
add_dependencies(Analysis Math NBody KDTree)
