# CMakeLists.txt to build the Analysis library
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2018
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Rodrigo Tobar
#
# This file is part of VELOCIraptor.

set(ANALYSIS_SRC
    Density.cxx
    Morphology.cxx
    Power.cxx
)
add_nbody_lib(Analysis "${ANALYSIS_SRC}")
