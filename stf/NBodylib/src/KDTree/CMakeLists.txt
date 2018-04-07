# CMakeLists.txt to build the KDTree library
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2018
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Rodrigo Tobar
#
# This file is part of VELOCIraptor.

set(KDTREE_SRC
    KDCalcSmoothQuantities.cxx
    KDFindNearest.cxx
    KDFOF.cxx
    KDLeafNode.cxx
    KDSplitNode.cxx
    KDTree.cxx
)
add_nbody_lib(KDTree "${KDTREE_SRC}")
