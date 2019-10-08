#!/bin/bash
#this is a sample script to show you how you could run structure search on multiple snapshots and then build a structure (halo) tree.

#initial and final snapshot numbers
isnap=0
fsnap=100
nsnaps=`echo $isnap" "$fsnap|awk '{print $2-$1+1}'`
#number of input files
nfiles=1
#base config parameter file to be used
paramfile=stf.base.param
#cosmology or simulation name
simname=lcdm
#input directory
indir=./
#output dir
outdir=./
#code directory
codedir=./
#stf executable
stfexe=${vrdir}/bin/stf
#tree executable
treefrogexe=${treefrogdir}/bin/treefrog

echo $isnap,$fsnap,$nsnaps

for ((j=$isnap; j<=$fsnap; j++))
do
    jj=`printf "%03d" $j`
    cp $paramfile $outdir/$simname.sn$jj.param;
    sed -i .old 's/Output=OUTNAME/Output='"$outdir"'/'"$simname"'.c'"$i"'.sn'"$jj"'/g' $outdir/$simname.sn$jj.param;
    sed -i .old 's/Snapshot_value=SNVALUE/Snapshot_value='"$j"'/g' $outdir/$simname.sn$jj.param;
    ifile=`printf "%s/snapshot_%03d" $indir $j`
    $stfexe -i $ifile -s $nfiles -C $outdir/$simname.sn$jj.param > $outdir/$simname.sn$jj.log;
done

#treefrog commands

#largest particle ID value
Neff=1024
Nid=`echo $Neff | awk '{print $1^3.0}'`
#number of steps used when linking
numsteps=4
siglimit=0.1
#to make sure halo ids temporally unique, use this value times snapshot,
halotemporalidval=10000000000
#specify format, 0 ascii, 1 binary, 2 hdf5
ibinary=0
#specify no separate field and subhalo files
ifield=0
#number of input velociraptor files (set by number of mpi threads) per snapshot
numfiles=1

rm $outdir/halolist.txt
for ((j=$isnap; j<=$fsnap; j++))
do
    jj=`printf "%03d" $j`
    echo $outdir/$simname.sn$jj >> $outdir/halolist.txt
done
$treefrogexe -i $outdir/halolist.txt -s $nsnaps -N $numfiles -n $Nid -t $numsteps -h $halotemporalidval -B $ibinary -F $ifield -o $outdir/$simname.tree $outdir/tree.log
