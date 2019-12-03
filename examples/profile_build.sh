#!/bin/bash -l
# This script profiles stf using perf and hotspot to produce flame graphs
# It produces a build directory with a certian labels
# runs cmake with the desired options, runs the code on the desired input
# using the desired config and running the hotspot stuff
# it assumes the existence of perf that is in linux-tools-common
# also makes use of wget

#script that produces lots of qsub scripts to run velociraptor on simulation output
if [ $# -eq 0 ] || [ "$1" == "--help" ]
then
    echo "This script profiles a VR run."
    echo "The interface is as follows:"
    echo "buildlabel buildoptions VRargs VRconfig"
    echo "buildlabel: string, label for the build"
    echo "buildoptions: string, cmake options for build. Ex: \" -DVR_USE_GAS=ON \""
    echo "VRargs: string, options for VR. Ex: \"-i inputfile -I 2 -s 1 -o outputfile \""
    echo "VRconfig: string, path and file name of the config file"
    exit
fi


#initial and final snapshot numbers
buildlabel=$1
buildoptions=$2
VRargs=$3
VRconfig=$4

workingdir=`pwd`
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

#build exe
cd ${scriptdir}/../
mkdir build-${buildlabel}
cd build-${buildlabel}
rm -rf *
cmake ${buildoptions} -DCMAKE_BUILD_TYPE=RelWithDebugInfo ../
make -j

wget https://github.com/KDAB/hotspot/releases/download/v1.1.0/hotspot-v1.1.0-x86_64.AppImage
chmod +x hotspot-v1.1.0-x86_64.AppImage

#run code
#this will have produced a perf.data
perf record ./stf ${VRargs} -C ${VRconfig}

#run hotspot
./hotspot-v1.1.0-x86_64.AppImage perf.data
