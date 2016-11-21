"""
This python library, Orb-Weaver, is designed to process VELOCIraptor Halo Finder + Merger Tree output and produce an orbit library and
even construct a new halo catalog (with modifications to CM positions, halo masses, identification of subhaloes). 
"""


import sys,os,os.path,string,time,re,struct
import math,operator
from pylab import *
import numpy as np
import h5py #import hdf5 interface
import tables as pytb #import pytables
import pandas as pd
from copy import deepcopy
from sklearn.neighbors import NearestNeighbors
import scipy.interpolate as scipyinterp
import scipy.spatial as spatial
import multiprocessing as mp
#append the path to the python tools directory
sys.path.append('../../tools/')
#import the python toosl in velociraptor_python_tools
from velociraptor_python_tools import *

"""
    Tools for parsing the halo catalog temporal information 
"""

def GetProgenLength(halodata,haloindex,halosnap,haloid,atime,HALOIDVAL,endreftime=-1):
    """
    Get the length of a halo's progenitors
    """
    proglen=1
    progid=halodata[halosnap]["Tail"][haloindex]
    progsnap=halodata[halosnap]["TailSnap"][haloindex]
    progindex=int(progid%HALOIDVAL-1)
    while (progid!=haloid):
        proglen+=1
        haloid=progid
        halosnap=progsnap
        haloindex=progindex
        progid=halodata[halosnap]["Tail"][haloindex]
        progsnap=halodata[halosnap]["TailSnap"][haloindex]
        progindex=int(progid%HALOIDVAL-1)
        if (atime[halosnap]<endreftime):break
    return proglen

"""
    Functions that add information to the halo catalog
"""

def IdentifyMergers(numsnaps,tree,numhalos,halodata,boxsize,hval,atime,MERGERMLIM=0.1,RADINFAC=1.2,RADOUTFAC=1.5,NPARTCUT=100, HALOIDVAL=1000000000000, iverbose=1,pos_tree=[]):
    """
    Using head/tail info in halodata dictionary identify mergers based on distance and mass ratios
    #todo still testing 

    """
    for j in range(numsnaps):
        #store id and snap and mass of last major merger and while we're at it, store number of major mergers
        halodata[j]["LastMerger"]=np.ones(numhalos[j],dtype=np.int64)*-1
        halodata[j]["LastMergerRatio"]=np.ones(numhalos[j],dtype=np.float64)*-1
        halodata[j]["LastMergerSnap"]=np.zeros(numhalos[j],dtype=np.uint32)
        halodata[j]["LastMergerDeltaSnap"]=np.zeros(numhalos[j],dtype=np.uint32)
        #halodata[j]["NumMergers"]=np.zeros(numhalos[j],dtype=np.uint32)
    #built KD tree to quickly search for near neighbours
    if (len(pos_tree)==0):
        pos=[[]for j in range(numsnaps)]
        pos_tree=[[]for j in range(numsnaps)]
        start=time.clock()
        if (iverbose): print "tree build"
        for j in range(numsnaps):
            if (numhalos[j]>0):
                boxval=boxsize*atime[j]/hval
                pos[j]=np.transpose(np.asarray([halodata[j]["Xc"],halodata[j]["Yc"],halodata[j]["Zc"]]))
                pos_tree[j]=spatial.cKDTree(pos[j],boxsize=boxval)
        if (iverbose): print "done ",time.clock()-start
    #else assume tree has been passed
    for j in range(numsnaps):
        if (numhalos[j]==0): continue
        #at snapshot look at all haloes that have not had a major merger set
        #note that only care about objects with certain number of particles
        partcutwdata=np.where(halodata[j]["npart"]>=NPARTCUT)
        mergercut=np.where(halodata[j]["LastMergerRatio"][partcutwdata]<0)
        hids=np.asarray(halodata[j]["ID"][partcutwdata][mergercut],dtype=np.uint64)
        start=time.clock()
        if (iverbose):print "Processing ", len(hids)
        if (len(hids)==0):continue

        for hidval in hids:
            #now for each object get the main progenitor
            haloid=np.uint64(hidval)
            haloindex=int(haloid%HALOIDVAL-1)
            halosnap=j
            originalhaloid=haloid
            progid=halodata[halosnap]["Tail"][haloindex]
            progsnap=halodata[halosnap]["TailSnap"][haloindex]
            progindex=int(progid%HALOIDVAL-1)
            numprog=tree[halosnap]["Num_progen"][haloindex]
            #if object has no progenitor set LastMergerRatio to 0 and LastMerger to 0
            if (numprog==0): 
                halodata[halosnap]["LastMerger"][haloindex]=0
                halodata[halosnap]["LastMergerRatio"][haloindex]=0
                continue
            #print "starting halos ",j, hidval
            #halo has main branch which we can wander on
            #while object is not its own progenitor move along tree to see how many major mergers it had across its history
            while (True):
                #now for each progenitor, lets find any nearby objects within a given mass/vmax interval
                posval=[halodata[progsnap]["Xc"][progindex],halodata[progsnap]["Yc"][progindex],halodata[progsnap]["Zc"][progindex]]
                radval=RADINFAC*halodata[progsnap]["R_200crit"][progindex]
                #get neighbour list within RADINFAC sorted by mass with most massive first
                NNlist=pos_tree[progsnap].query_ball_point(posval, radval)
                NNlist=[NNlist[ij] for ij in np.argsort(halodata[progsnap]["Mass_tot"][NNlist])[::-1]]
                #store boxval for periodic correction
                boxval=boxsize*atime[progsnap]/hval
                #now if list contains some objects, lets see if the velocity vectors are moving towards each other and mass/vmax ratios are okay
                if (len(NNlist)>0):
                    for NN in NNlist:
                        if (NN!=progindex):
                            mratio=halodata[progsnap]["Mass_tot"][NN]/halodata[progsnap]["Mass_tot"][progindex]
                            vratio=halodata[progsnap]["Vmax"][NN]/halodata[progsnap]["Vmax"][progindex]
                            #merger ratio is for object being larger of the two involved in merger
                            if (mratio>MERGERMLIM and mratio<1.0):
                                posvalrel=[halodata[progsnap]["Xc"][progindex]-halodata[progsnap]["Xc"][NN],halodata[progsnap]["Yc"][progindex]-halodata[progsnap]["Yc"][NN],halodata[progsnap]["Zc"][progindex]-halodata[progsnap]["Zc"][NN]]
                                for ij in range(3):
                                    if posvalrel[ij]<-0.5*boxval: posvalrel[ij]+=boxval
                                    elif posvalrel[ij]>0.5*boxval: posvalrel[ij]-=boxval
                                velvalrel=[halodata[progsnap]["VXc"][progindex]-halodata[progsnap]["VXc"][NN],halodata[progsnap]["VYc"][progindex]-halodata[progsnap]["VYc"][NN],halodata[progsnap]["VZc"][progindex]-halodata[progsnap]["VZc"][NN]]
                                radvelval=np.dot(posvalrel,velvalrel)/np.linalg.norm(posvalrel)
                                if (radvelval<0):
                                    #merger is happending
                                    #print "merger happening ", progsnap, NN

                                    #question of whether should move down the tree till merger no longer happening and define that as the start
                                    #this could also set the length of the merger
                                    #lets move along the tree of the infalling neighbour still it object is past the some factor of progenitor virial radius
                                    starthaloindex=progindex
                                    starthaloid=progid
                                    starthalosnap=progsnap
                                    startmergerindex=NN
                                    startmergerid=halodata[progsnap]["ID"][NN]
                                    startmergersnap=progsnap
                                    mergerstartindex=starthaloindex
                                    mergerstartid=starthaloid
                                    mergerstartsnap=starthalosnap
                                    while (tree[starthalosnap]["Num_progen"][starthaloindex]>0 and tree[startmergersnap]["Num_progen"][startmergerindex]>0):
                                        posvalrel=[halodata[starthalosnap]["Xc"][starthaloindex]-halodata[startmergersnap]["Xc"][startmergerindex],halodata[starthalosnap]["Yc"][starthaloindex]-halodata[startmergersnap]["Yc"][startmergerindex],halodata[starthalosnap]["Zc"][starthaloindex]-halodata[startmergersnap]["Zc"][startmergerindex]]
                                        boxval=boxsize*atime[starthalosnap]/hval
                                        for ij in range(3):
                                            if posvalrel[ij]<-0.5*boxval: posvalrel[ij]+=boxval
                                            elif posvalrel[ij]>0.5*boxval: posvalrel[ij]-=boxval
                                        radval=np.linalg.norm(posvalrel)/halodata[starthalosnap]["R_200crit"][starthaloindex]
                                        mratio=halodata[startmergersnap]["Mass_tot"][startmergerindex]/halodata[starthalosnap]["Mass_tot"][starthaloindex]

                                        #as moving back if halo now outside or too small, stop search and define this as start of merger
                                        if (radval>RADOUTFAC or mratio<MERGERMLIM):
                                            mergerstartindex=starthaloindex
                                            mergerstartid=starthaloid
                                            mergerstartsnap=starthalosnap
                                            break

                                        #move to next progenitors
                                        nextidval=halodata[starthalosnap]["Tail"][starthaloindex]
                                        nextsnapval=halodata[starthalosnap]["TailSnap"][starthaloindex]
                                        nextindexval=int(nextidval%HALOIDVAL-1)
                                        starthaloid=nextidval
                                        starthalosnap=nextsnapval
                                        starthaloindex=nextindexval

                                        nextidval=halodata[startmergersnap]["Tail"][startmergerindex]
                                        nextsnapval=halodata[startmergersnap]["TailSnap"][startmergerindex]
                                        nextindexval=int(nextidval%HALOIDVAL-1)
                                        startmergerid=nextidval
                                        startmergersnap=nextsnapval
                                        startmergerindex=nextindexval
                                    #store timescale of merger
                                    deltamergertime=(mergerstartsnap-progsnap)
                                    #set this as the merger for all halos from this point onwards till reach head or halo with non-zero merger
                                    merginghaloindex=mergerstartindex
                                    merginghaloid=mergerstartid
                                    merginghalosnap=mergerstartsnap
                                    oldmerginghaloid=merginghaloid
                                    #print "Merger found ",progsnap,mergerstartsnap, halodata[progsnap]["Mass_tot"][NN]/halodata[progsnap]["Mass_tot"][progindex], 
                                    #print halodata[startmergersnap]["Mass_tot"][startmergerindex]/halodata[starthalosnap]["Mass_tot"][starthaloindex]
                                    #now set merger time for all later haloes unless an new merger has happened
                                    while (oldmerginghaloid!=halodata[progsnap]["RootHead"][progindex] and halodata[merginghalosnap]["LastMergerRatio"][merginghaloindex]<0):
                                        halodata[merginghalosnap]["LastMerger"][merginghaloindex]=halodata[progsnap]["ID"][NN]
                                        halodata[merginghalosnap]["LastMergerRatio"][merginghaloindex]=halodata[progsnap]["Mass_tot"][NN]/halodata[progsnap]["Mass_tot"][progindex]
                                        halodata[merginghalosnap]["LastMergerSnap"][merginghaloindex]=progsnap
                                        halodata[merginghalosnap]["LastMergerDeltaSnap"][merginghaloindex]=deltamergertime

                                        oldmerginghaloid=merginghaloid
                                        mergingnextid=halodata[merginghalosnap]["Head"][merginghaloindex]
                                        mergingnextsnap=halodata[merginghalosnap]["HeadSnap"][merginghaloindex]
                                        mergingnextindex=int(mergingnextid%HALOIDVAL-1)
                                        merginghaloindex=mergingnextindex
                                        merginghaloid=mergingnextid
                                        merginghalosnap=mergingnextsnap
                #check if end of branch
                #move up and set last major merger to 0
                if (haloid==progid):
                    oldhaloid=haloid
                    currentsnap=halosnap
                    currentindex=haloindex
                    currentid=haloid
                    while (oldhaloid!=halodata[progsnap]["RootHead"][progindex] and halodata[currentsnap]["LastMergerRatio"][currentindex]<0):
                        halodata[currentsnap]["LastMerger"][currentindex]=0
                        halodata[currentsnap]["LastMergerRatio"][currentindex]=0
                        nextid=halodata[currentsnap]["Head"][currentindex]
                        nextsnap=halodata[currentsnap]["HeadSnap"][currentindex]
                        nextindex=int(nextid%HALOIDVAL-1)
                        oldhaloid=currentid
                        currentsnap=nextsnap
                        currentid=nextid
                        currentindex=nextindex
                    break
                #move to next step
                haloid=progid
                haloindex=progindex
                halosnap=progsnap
                progid=halodata[halosnap]["Tail"][haloindex]
                progsnap=halodata[halosnap]["TailSnap"][haloindex]
                progindex=int(progid%HALOIDVAL-1)
                numprog=tree[halosnap]["Num_progen"][haloindex]
        if (iverbose): print "Done snap",j,time.clock()-start
        return pos_tree

def IdentifyOrbits(numsnaps,tree,numhalos,halodata,boxsize,hval,atime,NPARTCUT=1000,HALOIDVAL=1000000000000, iverbose=1,pos_tree=[]):
    """
    Using head/tail info deteremine when subhaloes orbit their host, when they enter or leave etc. 
    #todo still testing 

    This adds the data blocks 
    OrbitID
    OrbitIDsnap
    NumOrbits
    OrbitPeriodDeltaA
    OrbitPeriodAstart 
    ClosestApproach
    MassAtAccretion
    VmaxAtAccretion

    and some crossing data

    NumInwardCrossing_R1.0
    NumInwardCrossing_R1.5
    NumInwardCrossing_R2.0
    NumInwardCrossing_R2.5
    NumInwardCrossing_R3.0

    NumOutwardCrossing_R1.0
    NumOutwardCrossing_R1.5
    NumOutwardCrossing_R2.0
    NumOutwardCrossing_R2.5
    NumOutwardCrossing_R3.0

    Note that the selection criteria is complex here. 
    A host halo is searched for all objects within some distance. It then tracks objects considered
    subhalos or halos that are within this radius to follow their paths. This ensures that if a 
    halo is falling into another and has its own subhaloes, these are not considered.

    The cut is made at 1.5 Rvir and traced backwards

    """

    #
    for j in range(numsnaps):
        #store the halo this measurments are in reference to
        halodata[j]["OrbitID"]=np.ones(numhalos[j],dtype=np.int64)*-1
        halodata[j]["OrbitIDSnap"]=np.zeros(numhalos[j],dtype=np.int32)
        #store the number of orbits an object went around this host 
        halodata[j]["NumOrbits"]=np.ones(numhalos[j],dtype=np.float64)*-1
        halodata[j]["OrbitPeriodDeltaA"]=np.ones(numhalos[j],dtype=np.float64)*-1
        halodata[j]["OrbitPeriodAstart"]=np.ones(numhalos[j],dtype=np.float64)*-1
        halodata[j]["ClosestApproach"]=np.ones(numhalos[j],dtype=np.float64)*-1
        #store accretion mass
        halodata[j]["MassAtAccretion"]=np.ones(numhalos[j],dtype=np.float64)*-1
        halodata[j]["VmaxAtAccretion"]=np.ones(numhalos[j],dtype=np.float64)*-1

        #store number of inward and outward crossing, makes it easy to identify backsplash galaxies
        halodata[j]["NumInwardCrossing_R1.0"]=np.zeros(numhalos[j],dtype=np.uint32)
        halodata[j]["NumInwardCrossing_R1.5"]=np.zeros(numhalos[j],dtype=np.uint32)
        halodata[j]["NumInwardCrossing_R2.0"]=np.zeros(numhalos[j],dtype=np.uint32)
        halodata[j]["NumInwardCrossing_R2.5"]=np.zeros(numhalos[j],dtype=np.uint32)
        halodata[j]["NumInwardCrossing_R3.0"]=np.zeros(numhalos[j],dtype=np.uint32)

        halodata[j]["NumOutwardCrossing_R1.0"]=np.zeros(numhalos[j],dtype=np.uint32)
        halodata[j]["NumOutwardCrossing_R1.5"]=np.zeros(numhalos[j],dtype=np.uint32)
        halodata[j]["NumOutwardCrossing_R2.0"]=np.zeros(numhalos[j],dtype=np.uint32)
        halodata[j]["NumOutwardCrossing_R2.5"]=np.zeros(numhalos[j],dtype=np.uint32)
        halodata[j]["NumOutwardCrossing_R3.0"]=np.zeros(numhalos[j],dtype=np.uint32)

    #built KD tree to quickly search for near neighbours
    if (len(pos_tree)==0):
        pos=[[]for j in range(numsnaps)]
        pos_tree=[[]for j in range(numsnaps)]
        start=time.clock()
        if (iverbose): print "tree build"
        for j in range(numsnaps):
            if (numhalos[j]>0):
                boxval=boxsize*atime[j]/hval
                pos[j]=np.transpose(np.asarray([halodata[j]["Xc"],halodata[j]["Yc"],halodata[j]["Zc"]]))
                pos_tree[j]=spatial.cKDTree(pos[j],boxsize=boxval)
        if (iverbose): print "done ",time.clock()-start
    #else assume tree has been passed
    for j in range(numsnaps):
        if (numhalos[j]==0): continue
        IdentifyOrbitsAtSnap(j, tree,numhalos,halodata,boxsize,hval,atime,NPARTCUT,HALOIDVAL, pos_tree, iverbose)

def IdentifyOrbitsAtSnap(snapval, tree,numhalos,halodata,boxsize,hval,atime,NPARTCUT,HALOIDVAL, pos_tree, iverbose):
    """
    Using head/tail info to calculate orbits of objects around halos at a given snapshot. See #ref IdentifyOrbits for data blocks
    produced. 
    """

    #at snapshot look host halos that are large enough to do some tracking of subhalos, secondary progenitors etc
    partcutwdata=np.where(halodata[snapval]["npart"]>=NPARTCUT)
    hosts=np.where(halodata[snapval]["hostHaloID"][partcutwdata]==-1)
    hids=np.asarray(halodata[snapval]["ID"][partcutwdata][hosts],dtype=np.uint64)
    start=time.clock()
    if (iverbose):print "Processing ", len(hids)
    if (len(hids)==0):
        if (iverbose): print "Done snap",snapval,time.clock()-start
        return 
    for hidval in hids:
        haloid=np.uint64(hidval)
        haloindex=int(haloid%HALOIDVAL-1)
        halosnap=snapval
        mainhaloid=haloid
        mainhalosnap=halosnap
        mainhaloindex=haloindex
        mainhalomass=halodata[halosnap]["Mass_tot"][haloindex]
        progid=halodata[halosnap]["Tail"][haloindex]
        progsnap=halodata[halosnap]["TailSnap"][haloindex]
        progindex=int(progid%HALOIDVAL-1)
        numprog=tree[halosnap]["Num_progen"][haloindex]
        mainprogid=progid
        mainprogsnap=progsnap
        mainprogindex=progindex

        #determine length of main branch
        proglength=GetProgenLength(halodata,haloindex,halosnap,haloid,atime,HALOIDVAL)
        #only keep going along halos that have existed along enough
        if (proglength<=10): continue

        haloid=mainhaloid
        halosnap=mainhalosnap
        haloindex=mainhaloindex
        progid=mainprogid
        progsnap=mainprogsnap
        progindex=mainprogindex
        #first store position of main halo as function of time
        endreftime=atime[halosnap]
        mainpos=np.zeros([6,proglength])
        mainatime=np.zeros(proglength)
        proglength=0
        while (haloid!=progid):
            afac=1.0/atime[halosnap]
            mainatime[proglength]=atime[halosnap]
            mainpos[0][proglength]=halodata[halosnap]["Xc"][haloindex]*afac
            mainpos[1][proglength]=halodata[halosnap]["Yc"][haloindex]*afac
            mainpos[2][proglength]=halodata[halosnap]["Zc"][haloindex]*afac
            mainpos[3][proglength]=halodata[halosnap]["VXc"][haloindex]*afac
            mainpos[4][proglength]=halodata[halosnap]["VYc"][haloindex]*afac
            mainpos[5][proglength]=halodata[halosnap]["VZc"][haloindex]*afac

            proglength+=1
            #store last time
            endreftime=atime[halosnap]

            haloid=progid
            halosnap=progsnap
            haloindex=progindex
            progid=halodata[halosnap]["Tail"][haloindex]
            progsnap=halodata[halosnap]["TailSnap"][haloindex]
            progindex=int(progid%HALOIDVAL-1)

        #interpolate position 
        posref=[scipyinterp.interp1d(mainatime,mainpos[ij]) for ij in range(6)]

        haloid=mainhaloid
        halosnap=mainhalosnap
        haloindex=mainhaloindex
        progid=mainprogid
        progsnap=mainprogsnap
        progindex=mainprogindex

        #now for each object get all surrounding subhalos within different radial bins
        posval=[halodata[halosnap]["Xc"][haloindex],halodata[halosnap]["Yc"][haloindex],halodata[halosnap]["Zc"][haloindex]]
        #now find all objects within this region
        radval=1.5*halodata[halosnap]["R_200crit"][haloindex]
        NNlist=pos_tree[halosnap].query_ball_point(posval, radval)
        radval=halodata[halosnap]["R_200crit"][haloindex]
        #proceed only if list is non zero
        if (len(NNlist)>0):
            #keep only smaller haloes and its own subhalos
            halolist=np.zeros(len(NNlist),dtype=np.uint64)
            nhalolist=0
            for NN in NNlist:
                hostid=halodata[halosnap]["hostHaloID"][NN]
                massval=halodata[halosnap]["Mass_tot"][NN]
                currentid=halodata[halosnap]["ID"][NN]
                if (hostid==-1 or hostid==haloid):
                    if (currentid!=mainhaloid and massval < mainhalomass):
                        halolist[nhalolist]=NN
                        nhalolist+=1
            #if there are halos that are subhalos or smaller halos in volume move on to next section
            if (nhalolist>0):
                #now have cleaned list of objects to examine
                for ihalo in range(nhalolist):
                    #then get the relative motion and fill orbit properties
                    GetHaloRelativeMotion(halolist[ihalo],mainhaloid,mainhalosnap,radval,halodata,boxsize,hval,atime,posref,endreftime,HALOIDVAL)

        #then look at secondary progenitors
        if (numprog>=2):
            wdata1=np.where(halodata[mainprogsnap]["Head"]==mainhaloid)
            progids=np.asarray(halodata[mainprogsnap]["ID"][wdata1],dtype=np.uint64)
            for progid in progids:
                if (progid!=mainprogid):
                    progindex=int(progid%HALOIDVAL-1)
                    GetHaloRelativeMotion(progindex,mainhaloid,mainprogsnap,radval,halodata,boxsize,hval,atime,posref,endreftime,HALOIDVAL)
    if (iverbose): print "Done snap",snapval,time.clock()-start

def GetHaloRelativeMotion(haloindexval,mainhaloid,mainhalosnap,mainhaloradval,halodata,boxsize,hval,atime,posref,endreftime, HALOIDVAL):
    """
    Examine the motion of an object relative to a reference position and some radial scale given by a mainhalo
    This is used to compute orbital parameters
    """
    #now for each halo trace its position relative to main host
    #idea is to trace back to point at which object 
    haloindex=haloindexval
    halosnap=mainhalosnap
    haloid=halodata[halosnap]["ID"][haloindex]

    #initialize number of orbits, orbit params etc
    halodata[halosnap]["NumOrbits"][haloindex]=0
    halodata[halosnap]["OrbitID"][haloindex]=mainhaloid
    halodata[halosnap]["OrbitIDSnap"][haloindex]=mainhalosnap

    progid=halodata[halosnap]["Tail"][haloindex]
    progsnap=halodata[halosnap]["TailSnap"][haloindex]
    progindex=int(progid%HALOIDVAL-1)

    #determine length of main branch
    proglength=GetProgenLength(halodata,haloindex,halosnap,haloid,atime,HALOIDVAL,endreftime)
    #if length is small do nothing
    if (proglength<=10): return

    poshalo=np.zeros([6,proglength])
    masshalo=np.zeros(proglength)
    vmaxhalo=np.zeros(proglength)
    atimehalo=np.zeros(proglength)
    jhalo=np.zeros([proglength,3])
    proglength=0

    haloindex=haloindexval
    halosnap=mainhalosnap
    haloid=halodata[halosnap]["ID"][haloindex]

    progid=halodata[halosnap]["Tail"][haloindex]
    progsnap=halodata[halosnap]["TailSnap"][haloindex]
    progindex=int(progid%HALOIDVAL-1)
    while (haloid!=progid):
        #move along till end
        if (atime[halosnap]<endreftime):break
        afac=1.0/atime[halosnap]
        boxval=boxsize*atime[halosnap]/hval
        refpos=[posref[ij](atime[halosnap]) for ij in range(6)]
        poshalo[0][proglength]=halodata[halosnap]["Xc"][haloindex]*afac-refpos[0]
        poshalo[1][proglength]=halodata[halosnap]["Yc"][haloindex]*afac-refpos[1]
        poshalo[2][proglength]=halodata[halosnap]["Zc"][haloindex]*afac-refpos[2]
        poshalo[3][proglength]=halodata[halosnap]["VXc"][haloindex]*afac-refpos[3]
        poshalo[4][proglength]=halodata[halosnap]["VYc"][haloindex]*afac-refpos[4]
        poshalo[5][proglength]=halodata[halosnap]["VZc"][haloindex]*afac-refpos[5]
        masshalo[proglength]=halodata[halosnap]["Mass_tot"][haloindex]
        vmaxhalo[proglength]=halodata[halosnap]["Vmax"][haloindex]
        atimehalo[proglength]=atime[halosnap]
        jhalo[proglength]=np.cross([poshalo[0][proglength],poshalo[1][proglength],poshalo[2][proglength]],[poshalo[3][proglength],poshalo[4][proglength],poshalo[5][proglength]])

        #correction positions for period
        if (poshalo[0][proglength]<-0.5*boxval): poshalo[0][proglength]+=boxval
        if (poshalo[1][proglength]<-0.5*boxval): poshalo[1][proglength]+=boxval
        if (poshalo[2][proglength]<-0.5*boxval): poshalo[2][proglength]+=boxval
        if (poshalo[0][proglength]>0.5*boxval): poshalo[0][proglength]-=boxval
        if (poshalo[1][proglength]>0.5*boxval): poshalo[1][proglength]-=boxval
        if (poshalo[2][proglength]>0.5*boxval): poshalo[2][proglength]-=boxval
        proglength+=1

        haloid=progid
        halosnap=progsnap
        haloindex=progindex
        progid=halodata[halosnap]["Tail"][haloindex]
        progsnap=halodata[halosnap]["TailSnap"][haloindex]
        progindex=int(progid%HALOIDVAL-1)

    radhalo=np.sqrt(poshalo[0]*poshalo[0]+poshalo[1]*poshalo[1]+poshalo[2]*poshalo[2])
    vrhalo=(poshalo[0]*poshalo[3]+poshalo[1]*poshalo[4]+poshalo[2]*poshalo[5])/radhalo

    #now have relative positions and velocities to construct orbits look at number of times object crosses particular 
    #boundaries
    haloindex=haloindexval
    halosnap=mainhalosnap
    haloid=halodata[halosnap]["ID"][haloindex]
    halodata[halosnap]["ClosestApproach"][haloindex]=min(radhalo)

    #look at sign of radial velocity and determine if it switches at any point
    #first look if all vr are negative (first infall)
    astart=aend=0
    if (len(np.where(vrhalo<0)[0])<proglength):
        #halo has at least had turn around so for all radii within 2.0 mainhaloradval lets look at number of orbits
        astart=atimehalo[0]
        for i in range(proglength-1):
            #stop when all radii are further away than 2.0 * radial scale
            if (np.sum(radhalo[i:]>=2.0*mainhaloradval)==proglength-i):
                if (aend==0): aend=atimehalo[i]
                break
            if (vrhalo[i]*vrhalo[i+1]<=0): 
                halodata[halosnap]["NumOrbits"][haloindex]+=0.5
            if (halodata[halosnap]["NumOrbits"][haloindex]==1):aend=atimehalo[i]
    halodata[halosnap]["OrbitPeriodAstart"][haloindex]=astart
    halodata[halosnap]["OrbitPeriodDeltaA"][haloindex]=astart-aend

    #now with orbits, the radial distance as a function of time could be fit with a sinusoid 
    #with a decaying amplitude and a trend towards smaller radii
    #however, can only be reasonably examined for objects that have completed two or more orbits
    #this would be a good way of selecting objects for further analysis

    haloindex=haloindexval
    for i in range(proglength-1):
        if (radhalo[i]<mainhaloradval and radhalo[i+1]>mainhaloradval):
            halodata[halosnap]["NumInwardCrossing_R1.0"][haloindex]+=1
        if (radhalo[i]>mainhaloradval and radhalo[i+1]<mainhaloradval):
            halodata[halosnap]["NumOutwardCrossing_R1.0"][haloindex]+=1

        if (radhalo[i]<=1.5*mainhaloradval and radhalo[i+1]>1.5*mainhaloradval):
            halodata[halosnap]["NumInwardCrossing_R1.5"][haloindex]+=1
        if (radhalo[i]>1.5*mainhaloradval and radhalo[i+1]<=1.5*mainhaloradval):
            halodata[halosnap]["NumOutwardCrossing_R1.5"][haloindex]+=1

        if (radhalo[i]<=2.0*mainhaloradval and radhalo[i+1]>2.0*mainhaloradval):
            halodata[halosnap]["NumInwardCrossing_R2.0"][haloindex]+=1
        if (radhalo[i]>2.0*mainhaloradval and radhalo[i+1]<=2.0*mainhaloradval):
            halodata[halosnap]["NumOutwardCrossing_R2.0"][haloindex]+=1

        if (radhalo[i]<=2.5*mainhaloradval and radhalo[i+1]>2.5*mainhaloradval):
            halodata[halosnap]["NumInwardCrossing_R2.5"][haloindex]+=1
        if (radhalo[i]>2.5*mainhaloradval and radhalo[i+1]<=2.5*mainhaloradval):
            halodata[halosnap]["NumOutwardCrossing_R2.5"][haloindex]+=1

        if (radhalo[i]<=3.0*mainhaloradval and radhalo[i+1]>3.0*mainhaloradval):
            halodata[halosnap]["NumInwardCrossing_R3.0"][haloindex]+=1
        if (radhalo[i]>3.0*mainhaloradval and radhalo[i+1]<=3.0*mainhaloradval):
            halodata[halosnap]["NumOutwardCrossing_R3.0"][haloindex]+=1

        #find point at which the halo is always further away than some factor times Rvir
        #that is defined as the accretion point
        if (min(radhalo[i:])>mainhaloradval and halodata[halosnap]["MassAtAccretion"][haloindex]==0):
            halodata[halosnap]["MassAtAccretion"][haloindex]=masshalo[i]
            halodata[halosnap]["VmaxAtAccretion"][haloindex]=vmaxhalo[i]
