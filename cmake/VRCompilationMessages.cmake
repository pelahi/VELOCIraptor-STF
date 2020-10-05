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

macro(vr_config_errors)
    if (VR_ZOOM_SIM AND VR_NO_MASS)
        message(FATAL_ERROR "VR compiled to not store mass and also compiled for zoom simulations.
        These options are incompatible. Use one or the other.
        Options are VR_ZOOM_SIM and VR_NO_MASS.")
    endif()
endmacro()

macro(vr_compilation_summary)
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
            "Long Integers" LONG_INT)
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

	string(TOUPPER CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE} CMAKE_BUILD_FLAGS)
	message("")
	message("Compilation")
	message("=========================")
	message("Compiler ID          : ${CMAKE_CXX_COMPILER_ID}" )
	message("Compiler             : ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ID}" )
	message("Build Type           : ${CMAKE_BUILD_TYPE}")
	message("Build Flags          : ${${CMAKE_BUILD_FLAGS}}")
	if(NOT ENABLE_SANITIZER STREQUAL "none")
		message("Sanitizer            : -fsanitize=${ENABLE_SANITIZER}")
	endif()
	message("----------------------")
	message("Include directories  : ${VR_INCLUDE_DIRS}")
	message("VR Macros         : ${VR_DEFINES}")
	message("VR Lib            : ${VR_LIBS}")
	message("VR C++ flags      : ${VR_CXX_FLAGS}")
	message("VR Link flags     : ${VR_LINK_FLAGS}")
	message("")

endmacro()
