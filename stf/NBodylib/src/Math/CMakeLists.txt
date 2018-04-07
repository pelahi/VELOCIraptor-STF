# CMakeLists.txt to build the Math library
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2018
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Rodrigo Tobar
#
# This file is part of VELOCIraptor.

set(MATH_SRC
    Coordinate2D.cxx
    Coordinate.cxx
    Fitting.cxx
    GMatrix.cxx
    Integrate.cxx
    Interpolate.cxx
    Matrix2D.cxx
    Matrix.cxx
    Random.cxx
    SpecialFunctions.cxx
    Statistics.cxx
)
add_nbody_lib(Math "${MATH_SRC}")
