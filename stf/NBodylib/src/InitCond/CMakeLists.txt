# CMakeLists.txt to build the InitCond library
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2018
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Rodrigo Tobar
#
# This file is part of VELOCIraptor.

set(INITCOND_SRC
    InitCond.cxx
)
add_nbody_lib(InitCond "${INITCOND_SRC}")
