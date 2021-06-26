# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2018
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Rodrigo Tobar
#
# This file is part of VELOCIraptor.

set(tests
    test_h5_output_file
)

foreach(test ${tests})
  add_executable(${test} ${test}.cxx)
  target_link_libraries(${test} nbodylib_iface velociraptor ${VR_LIBS})
  if (VR_LINK_FLAGS)
    set_target_properties(${test} PROPERTIES LINK_FLAGS ${VR_LINK_FLAGS})
  endif()
endforeach()
