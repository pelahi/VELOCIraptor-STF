#
# Some useful utils
#

#
# list values as bullet points
#
function(list_to_bulletpoints result)
    list(REMOVE_AT ARGV 0)
    set(temp "")
    foreach(item ${ARGV})
        set(temp "${temp}* ${item}")
    endforeach()
    set(${result} "${temp}" PARENT_SCOPE)
endfunction(list_to_bulletpoints)

#
# valid the option choosen based on allowed values
#
function(validate_option name values)
    string(TOLOWER ${${name}} needle_lower)
    string(TOUPPER ${${name}} needle_upper)
    list(FIND ${values} ${needle_lower} IDX_LOWER)
    list(FIND ${values} ${needle_upper} IDX_UPPER)
    if(${IDX_LOWER} LESS 0 AND ${IDX_UPPER} LESS 0)
        list_to_bulletpoints(POSSIBLE_VALUE_LIST ${${values}})
        message(FATAL_ERROR "\n########################################################################\n"
                            "Invalid value '${${name}}' for option ${name}\n"
                            "Possible values are : "
                            "${POSSIBLE_VALUE_LIST}"
                            "\n"
                            "########################################################################")
    endif()
endfunction(validate_option)

#
# Function to add a "doc" target, which will doxygen out the given
# list of directories
#
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

#
# - Prevent in-source builds.
# https://stackoverflow.com/questions/1208681/with-cmake-how-would-you-disable-in-source-builds/
#
macro(prevent_in_source_builds)
    # make sure the user doesn't play dirty with symlinks
    get_filename_component(srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
    get_filename_component(srcdir2 "${CMAKE_SOURCE_DIR}/.." REALPATH)
    get_filename_component(srcdir3 "${CMAKE_SOURCE_DIR}/../src" REALPATH)
    get_filename_component(bindir "${CMAKE_BINARY_DIR}" REALPATH)

    # disallow in-source builds
    if("${srcdir}" STREQUAL "${bindir}" OR "${srcdir2}" STREQUAL "${bindir}" OR "${srcdir3}" STREQUAL "${bindir}")
        message(FATAL_ERROR "\
            CMake must not to be run in the source directory. \
            Rather create a dedicated build directory and run CMake there. \
            To clean up after this aborted in-place compilation:
            rm -r CMakeCache.txt CMakeFiles
        ")
    endif()
endmacro()

#
# set the default build and also store the compilation flags
# as a string based on the currently choosen flags
#
macro(my_set_build_type)
	set(default_build_type "Release")
	if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
		message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
		set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
			STRING "Choose the type of build." FORCE)
		# Set the possible values of build type for cmake-gui
		set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
			"Debug" "Release" "MinSizeRel" "RelWithDebInfo")
	endif()
	#set(ACTIVE_COMPILE_OPTIONS )
endmacro()

macro(enable_santizer_option)
    set(ENABLE_SANITIZER "none" CACHE STRING "Select a code sanitizer option (none (default), address, leak, thread, undefined)")
    mark_as_advanced(ENABLE_SANITIZER)
    set(ENABLE_SANITIZER_VALUES none address leak thread undefined)
    set_property(CACHE ENABLE_SANITIZER PROPERTY STRINGS ${ENABLE_SANITIZER_VALUES})
    validate_option(ENABLE_SANITIZER ENABLE_SANITIZER_VALUES)
    string(TOLOWER ${ENABLE_SANITIZER} ENABLE_SANITIZER)
endmacro()

function(sanitizer_options mytarget)
    if(NOT ENABLE_SANITIZER STREQUAL "none")
        if((${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU") OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
            target_compile_options(${mytarget} PUBLIC -fsanitize=${ENABLE_SANITIZER})
            target_link_options(${mytarget} PUBLIC -fsanitize=${ENABLE_SANITIZER})
        else()
            message(WARNING "ENABLE_SANITIZER option not supported by ${CMAKE_CXX_COMPILER_ID} compilers. Ignoring.")
            set(ENABLE_SANITIZER "none")
        endif()
    endif()
endfunction()

#
# Make sure we have the git submodules we need
#
macro(vr_ensure_git_submodules)
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
# How we find GSL and set it up
#
macro(find_gsl)
	find_package(GSL REQUIRED)
	list(APPEND VR_INCLUDE_DIRS ${GSL_INCLUDE_DIRS})
	list(APPEND VR_LIBS ${GSL_LIBRARIES})
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

macro(vr_mpi)
    set(VR_HAS_MPI No)
    if (VR_MPI)
    	find_mpi()
    endif()
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

macro(vr_hdf5)
    set(VR_HAS_HDF5 No)
    set(VR_HAS_COMPRESSED_HDF5 No)
    set(VR_HAS_PARALLEL_HDF5 No)
    if (VR_HDF5)
    	find_hdf5()
    endif()
endmacro()

macro(vr_nbodylib)
    # This provides us with the nbodylib library
    # We need to add it unless it was already added by somebody else
    if (NOT TARGET nbodylib)
    	add_subdirectory(NBodylib)
    	if (NBODYLIB_VERSION VERSION_LESS "1.30")
    		message(FATAL_ERROR "NBodyLib version ${NBODYLIB_VERSION} unsupported,
    		VELOCIraptor requires >= 1.30, try running git submodule update --recursive --remote")
    	endif()
    	list(INSERT VR_DOC_DIRS 0 ${NBODYLIB_DOC_DIRS})
    endif()
    list(APPEND VR_INCLUDE_DIRS ${NBODYLIB_INCLUDE_DIRS})
    list(APPEND VR_DEFINES "${NBODYLIB_DEFINES}")
    list(APPEND VR_CXX_FLAGS "${NBODYLIB_CXX_FLAGS}")
    list(APPEND VR_LINK_FLAGS "${NBODYLIB_LINK_FLAGS}")
    list(APPEND VR_LIBS "${NBODYLIB_LIBS}")
endmacro()

#run some macros automatically
prevent_in_source_builds()
my_set_build_type()
enable_santizer_option()
