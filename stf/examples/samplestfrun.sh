#!/bin/bash
#this is a sample script to show you how you could run structure search on multiple snapshots and then build a structure (halo) tree.

#initial and final snapshot numbers
isnap=0
fsnap=100
nsnaps=`echo $isnap" "$fsnap|awk '{print $2-$1+1}'`
nfiles=1
#base config parameter file to be used
paramfile=stf.base.param
#cosmology or simulation name
simname=lcdm
#input directory
indir=./
#output dir
outdir=./
#stf executable
stfexe=./bin/stf
#tree executable
hmt=./bin/halomergertree

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

#largest particle ID value
Neff=1024
Nid=`echo $Neff | awk '{print $1^3.0}'`
echo $Neff,$Nid
#rm $outdir/halolist.txt
for ((j=$isnap; j<=$fsnap; j++)) 
do
    jj=`printf "%03d" $j`
    echo $outdir/$simname.sn$jj >>$outdir/halolist.txt
done
$hmt -i $outdir/halolist.txt -s $nsnaps -n $Nid -o $outdir/$simname.tree $outdir/tree.log


