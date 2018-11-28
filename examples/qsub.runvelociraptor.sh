#!/bin/bash -l

#script that produces lots of qsub scripts to run velociraptor on simulation output
if [ $# -eq 0 ] || [ "$1" == "--help" ]
then
    echo "This script submits lots of qsub jobs that run velociraptor."
    echo "The interface is as follows:"
    echo "inputdir outputdir snapshotbasename numfilespersnap initialsnap finalsnap stfconfigfile qsubbaseconfig nmpi nopenmp"
    exit 
fi

#pass the input and output directories
inputdir=$1
outputdir=$2
#also pass the number of files needed
snapbasename=$3
numfiles=$4
isnap=$5
fsnap=$6
#velociraptor config
velociraptorbaseconfig=$7
#basic qsub configuration not that given the diversity in PBS, qsubbaseconfig is not altered according to the input mpi and openmp 
#values passed here but they should agree to the ones present in the qsub script
qsubbaseconfig=$8
nmpi=$9
nopenmp=${10}
nprocs=$((${nmpi}*${nopenmp}))

#set the velociraptor executable
velociraptor=

mkdir ${outputdir}

for ((i=$isnap;i<=$fsnap;i++))
do
    ii=`printf %03d $i`
    file=${snapbasename}_${ii};
    velconfig=${outputdir}/velociraptor.${ii}.config
    cp ${velociraptorbaseconfig} ${velconfig}
    sed -i 's/SNVALUE/'"$i"'/g' ${velconfig}
    qsubconfig=${outputdir}/qsub.velociraptor.${ii}.sh
    cp ${qsubbaseconfig} ${qsubconfig}
    sed -i 's/JOBNAME/vel-'"${ii}"'/g' ${qsubconfig}
    echo "cd ${outputdir}" >> ${qsubconfig}
    #add the command to run to the qsub script
    echo "mpirun -np ${nmpi} ${velociraptor} -i ${inputdir}/${file} -s ${numfiles} -C ${velconfig} -o ${outputdir}/${file}.VELOCIraptor > ${outputdir}/${file}.VELOCIraptor.log" >> ${qsubconfig}
    qsub ${qsubconfig}
done

