#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""

    This python script reads a full tree and plots some tree statistics, like merger rates, merger trees
    
"""

import sys,os,string,time,re,struct
import math,operator
import random as rn
import matplotlib
matplotlib.use('Agg')
import numpy as np

#load python routines
#one of these paths should work
sys.path.append('/path/to/VELOCIraptor-STF/stf/tools/')
import velociraptor_python_tools as vpt

def dtda(a,Omegam,OmegaL):
    return 1.0/np.sqrt(Omegam*np.power(a,-3.0)+OmegaL)/a
def getCosmicTime(a1,a2,Omegam,OmegaL,hval):
    intval=scipyint.romberg(dtda,a1,a2,args=(Omegam,OmegaL,))
    return 9.77813/hval*intval

#load the data
#note that if you want you do not need to load everything into memory and can select a subset of fields/halo properites
#just set requestedfields=['ID','Mass_200crit'] for example to load only those properties
#otherwise, can leave requestedfields empty
requestedfields=[]
#not by default tis read routines places z=0 at index=0 but of course, you can reverse the order
atime,tree,numhalos,halodata,cosmodata,unitdata=vpt.ReadUnifiedTreeandHaloCatalog(treefilename,requestedfields)
numsnaps=len(atime)

""" 
    -----------------------------------------------
    Now for some sample analysis 
    -----------------------------------------------
"""
#this is the special value that allows the tree to be walked
TEMPORALHALOIDVAL=1000000000000

#here's an example of looking at the evolution, motion of a halo 
#lets select a sample of halos at z=0. 
selectedhalos=halodata[0]['ID'][np.where(halodata[0]['Mass_tot']*unitdata['Mass_to_solarmass']>1e12)]
len(numobjs)
for j in range(numobjs):
    #start off with some halo
    haloid=selectedhalos[j]
    #get the halo's index by looking at its id
    haloindex=int(haloid%TEMPORALHALOIDVAL-1)
    halosnap=0

    #store the values
    mainindex=haloindex
    mainsnap=halosnap
    mainid=haloid

    #move along main halo store properties
    progid=halodata[i][halosnap]["Tail"][haloindex]
    progindex=int(progid%TEMPORALHALOIDVAL-1)
    progsnap=int(np.floor(progid/TEMPORALHALOIDVAL))
    #now for a reverse ordered data set where you have the halos at z=0 located at [0] 
    progsnap=numsnaps-1-progsnap
    #otherwise, leave at is 
    
    #find out how long the object has lived for 
    proglength=1
    #we do this by checking if the progenitor id is not a halo's own id
    while (progid!=haloid):
        proglength+=1
        #move to the next halo in time
        haloid=progid
        halosnap=progsnap
        haloindex=progindex
        #get the next progenitor
        progid=halodata[i][halosnap]["Tail"][haloindex]
        progindex=int(progid%HALOIDVAL-1)
        progsnap=int(np.floor(progid/TEMPORALHALOIDVAL))
        progsnap=numsnaps-1-progsnap
    
    #now that we know how long an object has lived for,
    #lets do some stuff

    #reset 
    haloid=mainid
    haloindex=mainindex
    halosnap=mainsnap

    #allocate mem
    pos=np.zeros(proglength,3)
    vel=np.zeros(proglength,3)
    mass=np.zeros(proglength)

    #init the prog id
    progid=halodata[i][halosnap]["Tail"][haloindex]
    progindex=int(progid%HALOIDVAL-1)
    progsnap=int(np.floor(progid/TEMPORALHALOIDVAL))
    progsnap=numsnaps-1-progsnap

    isnapcount=0
    start=time.clock()
    while (haloid!=progid):
        pos[isnapcount]=np.array([halodata[i][halosnap]["Xc"][haloindex],halodata[i][halosnap]["Yc"][haloindex],halodata[i][halosnap]["Zc"][haloindex]])
        vel[isnapcount]=np.array([halodata[i][halosnap]["VXc"][haloindex],halodata[i][halosnap]["VYc"][haloindex],halodata[i][halosnap]["VZc"][haloindex]])
        mass[isnapcount]=[halodata[i][halosnap]["Mass_tot"][haloindex]
        isnapcount+=1
        haloid=progid
        halosnap=progsnap
        haloindex=progindex
        progid=halodata[i][halosnap]["Tail"][haloindex]
        progindex=int(progid%HALOIDVAL-1)
        progsnap=int(np.floor(progid/TEMPORALHALOIDVAL))
        progsnap=numsnaps-1-progsnap
    print "done",j,"with ",proglength,time.clock()-start
