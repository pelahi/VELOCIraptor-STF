# CMakeLists.txt to build the VELOCIraptor library/binaries
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2018
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Rodrigo Tobar
#
# This file is part of VELOCIraptor.

set(VR_SOURCES
    allvars.cxx
    bgfield.cxx
    buildandsortarrays.cxx
    endianutils.cxx
    fofalgo.cxx
    gadgetio.cxx
    haloproperties.cxx
    hdfio.cxx
    io.cxx
    localbgcomp.cxx
    localfield.cxx
    mpigadgetio.cxx
    mpihdfio.cxx
    mpinchiladaio.cxx
    mpiramsesio.cxx
    mpiroutines.cxx
    mpitipsyio.cxx
    mpivar.cxx
    nchiladaio.cxx
    omproutines.cxx
    ramsesio.cxx
    search.cxx
    swiftinterface.cxx
    substructureproperties.cxx
    tipsyio.cxx
    ui.cxx
    unbind.cxx
    utilities.cxx
)

add_library(velociraptor STATIC ${VR_SOURCES})
target_compile_definitions(velociraptor PRIVATE ${VR_DEFINES})
if (VR_CXX_FLAGS)
	set_target_properties(velociraptor PROPERTIES COMPILE_FLAGS ${VR_CXX_FLAGS})
endif()
set_target_properties(velociraptor PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


# Bringing NBodyLib objects into the velociraptor library
target_link_libraries(velociraptor nbodylib_iface)

add_executable(stf main.cxx)
target_link_libraries(stf velociraptor ${VR_LIBS})

target_compile_definitions(stf PRIVATE ${VR_DEFINES})
if (VR_CXX_FLAGS)
	set_target_properties(stf PROPERTIES COMPILE_FLAGS ${VR_CXX_FLAGS})
endif()
if (VR_LINK_FLAGS)
	set_target_properties(stf PROPERTIES LINK_FLAGS ${VR_LINK_FLAGS})
endif()
set_target_properties(stf PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
