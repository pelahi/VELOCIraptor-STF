#!/bin/bash
#
# Travis CI install script
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2018
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Pascal Elahi
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA 02111-1307  USA
#

fail() {
	echo $1 1>&2
	exit 1
}

# first we create a directory for the CMake binaries
DEPS_DIR="${TRAVIS_BUILD_DIR}/deps"
mkdir ${DEPS_DIR} && cd ${DEPS_DIR}
# we use wget to fetch the cmake binaries
travis_retry wget --no-check-certificate https://cmake.org/files/v3.13/cmake-3.13.4-Linux-x86_64.tar.gz
# extract the binaries; the output here is quite lengthy,
# so we swallow it to not clutter up the travis console
tar -xvf cmake-3.13.4-Linux-x86_64.tar.gz > /dev/null
mv cmake-3.13.4-Linux-x86_64 cmake-install
# add both the top-level directory and the bin directory from the archive
# to the system PATH. By adding it to the front of the path we hide the
# preinstalled CMake with our own.
PATH=${DEPS_DIR}/cmake-install:${DEPS_DIR}/cmake-install/bin:$PATH

cd ${TRAVIS_BUILD_DIR}
mkdir build
cd build

VR_CMAKE_OPTIONS="-DCMAKE_CXX_COMPILER=$COMPILER "
if [ "$USEGAS" = 1 ]; then
	VR_CMAKE_OPTIONS+=" -DVR_USE_GAS=ON "
fi
if [ "$USESTARS" = 1 ]; then
	VR_CMAKE_OPTIONS+=" -DVR_USE_STAR=ON "
fi
if [ "$USEBH" = 1 ]; then
	VR_CMAKE_OPTIONS+=" -DVR_USE_BH=ON "
fi
if [ "$USEEXTRADM" = 1 ]; then
	VR_CMAKE_OPTIONS+=" -DVR_USE_EXTRA_DM_PROPERTIES=ON "
fi
if [ "$USESWIFT" = 1 ]; then
	VR_CMAKE_OPTIONS+=" -DVR_USE_SWIFT_INTERFACE=ON -DCMAKE_CXX_FLAGS=-fPIC "
fi
if [ "$NOMASS" = 1 ]; then
	VR_CMAKE_OPTIONS+=" -DVR_NO_MASS=ON "
fi
if [ "$ZOOMSIM" = 1 ]; then
	VR_CMAKE_OPTIONS+=" -DVR_ZOOM_SIM=ON "
fi
if [ "$NOMPI" = 1 ]; then
	VR_CMAKE_OPTIONS+=" -DVR_MPI=OFF "
fi
if [ "$NOMPENMP" = 1 ]; then
	VR_CMAKE_OPTIONS+=" -DVR_OPENMP=OFF -DNBODY_OPENMP=OFF"
fi
if [ "$OPENACC" = 1 ]; then
	VR_CMAKE_OPTIONS+=" -DVR_OPENACC=ON -DNBODY_OPENACC=ON"
fi

# Go, go, go!
cmake .. ${VR_CMAKE_OPTIONS} || fail "cmake failed"
make all -j2 || fail "make failed"
cd ..
