# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 12:23:08 2018
This python script compares two (or more) catalogs using TreeFrog
and checks to see if there are consisten within some tolerance.

The general interface is to provide an input file that is a list of catalogs to compare
The code will then invoke a simple, single thread treefrog instance and do a
cross catalog comparison.

This is then read by python read tools and analysed

@author: pelahi
"""


import sys,os,string,time,re,struct
from subprocess import call
import math,operator
import random as rn
import numpy as np
from scipy.stats.mstats import mquantiles
from scipy.misc import comb
import scipy.interpolate as scipyinterp
import scipy.integrate as scipyint
import scipy.optimize as scipyopt
import scipy.special as scipysp
import itertools
from matplotlib import rc
from matplotlib.patches import Ellipse, Wedge, Arrow, FancyArrowPatch, Polygon, PathPatch, Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.path import Path
from mpl_toolkits.mplot3d import Axes3D
#import pp

#load python routines
#one of these paths should work
pathtovelociraptor=sys.argv[0].split('examples')[0]
sys.path.append(pathtovelociraptor+'/tools/')
import velociraptor_python_tools as vpt

#if __name__ == '__main__':

if (os.path.isfile(sys.argv[1])==False):
    print("Missing info file",sys.argv[1])
    exit(1)
if (os.path.isfile(sys.argv[2])==False):
    print("Missing tolerance file",sys.argv[2])
    exit(1)

#load the plot info file,
print("Reading info file", sys.argv[1])
infofile=open(sys.argv[1],"r")
#load the code and arguments that will be run
treefrogexe=infofile.readline().strip()
treefrogargs=infofile.readline().strip()

numsims=int(((infofile.readline()).split(" "))[0]) #number of simulations
labellist,filenamelist=[[] for i in range(numsims)],[[] for i in range(numsims)]
for i in range(numsims):
    data=(infofile.readline().strip()).split(' ')
    labellist[i]=data[0]
    filenamelist[i]=data[1]
infofile.close()

print("Reading tolerance file", sys.argv[2])
tolerancefile=open(sys.argv[2],"r")
tolerancefile.close()

tol=dict()
tol['numobjfrac']=0.005
tol['nomatchfrac']=0.005
tol['nomatchnpart']=30
tol['merit']=0.95
tol['nomatchnpartfracabove']=0.005

checksummary=dict()
checksummarykeylist=['numobjs','nomatch','nomatchsize','merits']

#run the treefrog code
#get all permutations of the list of inputs
for i in range(numsims-1):
    #permutelabelist=itertools.permutations(labellist[i+1:])
    #permutefilelist=itertools.permutations(filelist[i+1:])
    for j in range(i+1,numsims):
        catcompname="catcomplist."+labellist[i]+"."+labellist[j]+".txt"
        os.system("rm "+catcompname)
        os.system("echo "+filenamelist[j]+" >> "+catcompname)
        os.system("echo "+filenamelist[i]+" >> "+catcompname)
        outputname='catcomp.'+labellist[i]+'.'+labellist[j]
        if (os.path.exists(outputname)==False):
            cmd=treefrogexe+" "+treefrogargs+' -i '+catcompname+' -o '+outputname
            os.system(cmd)

#set the dictionary to store the catalogs
trees=dict()
for i in range(numsims):
    trees[labellist[i]]=dict()
    checksummary[labellist[i]]=dict()
    for j in range(i+1,numsims):
        trees[labellist[i]][labellist[j]]=[]
        checksummary[labellist[i]][labellist[j]]=dict()
for i in range(numsims):
    for j in range(i+1,numsims):
        outputname='catcomp.'+labellist[i]+'.'+labellist[j]
        trees[labellist[i]][labellist[j]]=vpt.ReadHaloMergerTree(outputname,0,0,True,True)

error = False

#now process and check to see if catcomps pass
for i in range(numsims):
    for j in range(i+1,numsims):
        n1=len(trees[labellist[i]][labellist[j]][0]['Num_progen'])
        n2=len(trees[labellist[i]][labellist[j]][1]['Num_progen'])

        print(labellist[i],labellist[j],'Missing matches:')
        wnomatch=np.where(trees[labellist[i]][labellist[j]][0]['Num_progen']==0)
        nnomatch=len(wnomatch[0])
        if (nnomatch>0):
            npartdata=np.array(trees[labellist[i]][labellist[j]][0]['Npart'][wnomatch])
            iddata=np.array(trees[labellist[i]][labellist[j]][0]['haloID'][wnomatch])
            npartstats=np.concatenate([np.array([max(npartdata),min(npartdata)], dtype=np.float32),np.array(np.percentile(npartdata,[16,50,84,2.5,97.5]),dtype=np.float32)])
        if (len(wnomatch[0])==0 and n1==n2):
            print('pass, all objects have match' )
        else:
            #check if number of object difference is issue
            if (n1!=n2):
                print('catalogs have different number of objects',n1,n2)
                if (np.abs(n1-n2)/float(n1+n2)>tol['numobjfrac']):
                    print('FAIL, too large a difference between catalogs')
                    error = True
                else:
                    print('PASS')
            #check if number of missing matches and issue
            if (nnomatch>0):
                #get fraction of missing matches above threshold
                fracabove=float(len(np.where(npartdata>=tol['nomatchnpart'])[0]))/float(nnomatch)
                print('catalog 1 -> 2 produces missing matches',nnomatch)
                if (nnomatch/float(n1)>tol['nomatchfrac']):
                    print('FAIL, too many missing matches')
                    error = True
                else:
                    print('PASS')
                #check if number of particles of objects with missing matches is an issue
                if (npartstats[1]>tol['nomatchnpart']):
                    print('FAIL, smallest missing matches too large',npartstats[1],tol['nomatchnpart'])
                    error = True
                if (npartstats[0]>tol['nomatchnpart']):
                    print('FAIL, largest missing matches too large',npartstats[0],tol['nomatchnpart'])
                    error = True
                if (npartstats[3]>tol['nomatchnpart']):
                    print('FAIL, average missing matches too large',npartstats[3],tol['nomatchnpart'])
                    error = True
                if (fracabove>tol['nomatchnpartfracabove']):
                    print('FAIL, fraction of groups above',tol['nomatchnpart'],' particles is too high',fracabove,tol['nomatchnpartfracabove'])
                    error = True
                print('Missing objects have :')
                print('Npart:',npartdata)
                print('ID:',iddata)

        #merit check
        wmatch=np.where(trees[labellist[i]][labellist[j]][0]['Num_progen']>0)
        nmatch=len(wmatch[0])
        print(labellist[i],labellist[j],'Merits of matches:')
        if (nmatch>0):
            meritdata=np.array([trees[labellist[i]][labellist[j]][0]['Merit'][w][0] for w in wmatch[0]])
            npartdata=np.array([trees[labellist[i]][labellist[j]][0]['Npart'][w] for w in wmatch[0]])
            iddata=np.array([trees[labellist[i]][labellist[j]][0]['haloID'][w] for w in wmatch[0]])
            meritstats=np.concatenate([np.array([max(meritdata),min(meritdata)], dtype=np.float32),np.array(np.percentile(meritdata,[16,50,84,2.5,97.5]),dtype=np.float32)])
            if (meritstats[1]<tol['merit']):
                print('FAIL, lowest merit too small',meritstats[1],tol['merit'])
                error = True
            if (meritstats[0]<tol['merit']):
                print('FAIL, largest merit too small',meritstats[0],tol['merit'])
                error = True
            if (meritstats[3]<tol['merit']):
                print('FAIL, average merit too small',meritstats[3],tol['merit'])
                error = True
            wdata=np.where(meritdata<tol['merit'])
            if (len(wdata[0])>0):
                print('Objects that fail tolerance are :')
                print('Merit:',meritdata[wdata])
                print('Npart:',npartdata[wdata])
                print('ID:',iddata[wdata])

# Return an overall PASS or FAIL
if error:
    print('\n*********************')
    print('* Comparison FAILED *')
    print('*********************\n')
    exit(1)
else:
    print('\n*********************')
    print('* Comparison PASSED *')
    print('*********************\n')
    exit(0)
