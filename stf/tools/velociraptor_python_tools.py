#Make backwards compatible with python 2, ignored in python 3
from __future__ import print_function

import sys,os,os.path,string,time,re,struct
import math,operator
from pylab import *
import numpy as np
import h5py #import hdf5 interface
import tables as pytb #import pytables
import pandas as pd
from copy import deepcopy
from collections import deque
import itertools
import scipy.interpolate as scipyinterp
import scipy.spatial as spatial
import multiprocessing as mp
from collections import deque

import cython
from cython.parallel import prange, parallel

#would be good to compile these routines with cython
#try to speed up search
#cimport numpy as np

"""

Routines for reading velociraptor output

"""

"""
	IO Routines
"""

def ReadPropertyFile(basefilename,ibinary=0,iseparatesubfiles=0,iverbose=0, desiredfields=[], isiminfo=True, iunitinfo=True):
	"""
	VELOCIraptor/STF files in various formats
	for example ascii format contains
	a header with
		filenumber number_of_files
		numhalos_in_file nnumhalos_in_total
	followed by a header listing the information contain. An example would be
		ID(1) ID_mbp(2) hostHaloID(3) numSubStruct(4) npart(5) Mvir(6) Xc(7) Yc(8) Zc(9) Xcmbp(10) Ycmbp(11) Zcmbp(12) VXc(13) VYc(14) VZc(15) VXcmbp(16) VYcmbp(17) VZcmbp(18) Mass_tot(19) Mass_FOF(20) Mass_200mean(21) Mass_200crit(22) Mass_BN97(23) Efrac(24) Rvir(25) R_size(26) R_200mean(27) R_200crit(28) R_BN97(29) R_HalfMass(30) Rmax(31) Vmax(32) sigV(33) veldisp_xx(34) veldisp_xy(35) veldisp_xz(36) veldisp_yx(37) veldisp_yy(38) veldisp_yz(39) veldisp_zx(40) veldisp_zy(41) veldisp_zz(42) lambda_B(43) Lx(44) Ly(45) Lz(46) q(47) s(48) eig_xx(49) eig_xy(50) eig_xz(51) eig_yx(52) eig_yy(53) eig_yz(54) eig_zx(55) eig_zy(56) eig_zz(57) cNFW(58) Krot(59) Ekin(60) Epot(61) n_gas(62) M_gas(63) Xc_gas(64) Yc_gas(65) Zc_gas(66) VXc_gas(67) VYc_gas(68) VZc_gas(69) Efrac_gas(70) R_HalfMass_gas(71) veldisp_xx_gas(72) veldisp_xy_gas(73) veldisp_xz_gas(74) veldisp_yx_gas(75) veldisp_yy_gas(76) veldisp_yz_gas(77) veldisp_zx_gas(78) veldisp_zy_gas(79) veldisp_zz_gas(80) Lx_gas(81) Ly_gas(82) Lz_gas(83) q_gas(84) s_gas(85) eig_xx_gas(86) eig_xy_gas(87) eig_xz_gas(88) eig_yx_gas(89) eig_yy_gas(90) eig_yz_gas(91) eig_zx_gas(92) eig_zy_gas(93) eig_zz_gas(94) Krot_gas(95) T_gas(96) Zmet_gas(97) SFR_gas(98) n_star(99) M_star(100) Xc_star(101) Yc_star(102) Zc_star(103) VXc_star(104) VYc_star(105) VZc_star(106) Efrac_star(107) R_HalfMass_star(108) veldisp_xx_star(109) veldisp_xy_star(110) veldisp_xz_star(111) veldisp_yx_star(112) veldisp_yy_star(113) veldisp_yz_star(114) veldisp_zx_star(115) veldisp_zy_star(116) veldisp_zz_star(117) Lx_star(118) Ly_star(119) Lz_star(120) q_star(121) s_star(122) eig_xx_star(123) eig_xy_star(124) eig_xz_star(125) eig_yx_star(126) eig_yy_star(127) eig_yz_star(128) eig_zx_star(129) eig_zy_star(130) eig_zz_star(131) Krot_star(132) tage_star(133) Zmet_star(134)

	then followed by data

	Note that a file will indicate how many files the total output has been split into

	Not all fields need be read in. If only want specific fields, can pass a string of desired fields like
	['ID', 'Mass_FOF', 'Krot']
	#todo still need checks to see if fields not present and if so, not to include them or handle the error
	"""
	#this variable is the size of the char array in binary formated data that stores the field names
	CHARSIZE=40

	start = time.clock()
	inompi=True
	if (iverbose): print("reading properties file",basefilename)
	filename=basefilename+".properties"
	#load header
	if (os.path.isfile(filename)==True):
		numfiles=0
	else:
		filename=basefilename+".properties"+".0"
		inompi=False
		if (os.path.isfile(filename)==False):
			print("file not found")
			return []
	byteoffset=0
	#used to store fields, their type, etc
	fieldnames=[]
	fieldtype=[]
	fieldindex=[]

	if (ibinary==0):
		#load ascii file
		halofile = open(filename, 'r')
		#read header information
		[filenum,numfiles]=halofile.readline().split()
		filenum=int(filenum);numfiles=int(numfiles)
		[numhalos, numtothalos]= halofile.readline().split()
		numhalos=np.uint64(numhalos);numtothalos=np.uint64(numtothalos)
		names = ((halofile.readline())).split()
		#remove the brackets in ascii file names
		fieldnames= [fieldname.split("(")[0] for fieldname in names]
		for i in np.arange(fieldnames.__len__()):
			fieldname=fieldnames[i]
			if fieldname in ["ID","hostHalo","numSubStruct","npart","n_gas","n_star"]:
				fieldtype.append(np.uint64)
			elif fieldname in ["ID_mbp"]:
				fieldtype.append(np.int64)
			else:
				fieldtype.append(np.float64)
		halofile.close()
		#if desiredfields is NULL load all fields
		#but if this is passed load only those fields
		if (len(desiredfields)>0):
			lend=len(desiredfields)
			fieldindex=np.zeros(lend,dtype=int)
			desiredfieldtype=[[] for i in range(lend)]
			for i in range(lend):
				fieldindex[i]=fieldnames.index(desiredfields[i])
				desiredfieldtype[i]=fieldtype[fieldindex[i]]
			fieldtype=desiredfieldtype
			fieldnames=desiredfields
		#to store the string containing data format
		fieldtypestring=''
		for i in np.arange(fieldnames.__len__()):
			if fieldtype[i]==np.uint64: fieldtypestring+='u8,'
			elif fieldtype[i]==np.int64: fieldtypestring+='i8,'
			elif fieldtype[i]==np.float64: fieldtypestring+='f8,'

	elif (ibinary==1):
		#load binary file
		halofile = open(filename, 'rb')
		[filenum,numfiles]=np.fromfile(halofile,dtype=np.int32,count=2)
		[numhalos,numtothalos]=np.fromfile(halofile,dtype=np.uint64,count=2)
		headersize=np.fromfile(halofile,dtype=np.int32,count=1)[0]
		byteoffset=np.dtype(np.int32).itemsize*3+np.dtype(np.uint64).itemsize*2+4*headersize
		for i in range(headersize):
			fieldnames.append(unpack('s', halofile.read(CHARSIZE)).strip())
		for i in np.arange(fieldnames.__len__()):
			fieldname=fieldnames[i]
			if fieldname in ["ID","hostHalo","numSubStruct","npart","n_gas","n_star"]:
				fieldtype.append(np.uint64)
			elif fieldname in ["ID_mbp"]:
				fieldtype.append(np.int64)
			else:
				fieldtype.append(np.float64)
		halofile.close()
		#if desiredfields is NULL load all fields
		#but if this is passed load only those fields
		if (len(desiredfields)>0):
			lend=len(desiredfields)
			fieldindex=np.zeros(lend,dtype=int)
			desiredfieldtype=[[] for i in range(lend)]
			for i in range(lend):
				fieldindex[i]=fieldnames.index(desiredfields[i])
				desiredfieldtype[i]=fieldtype[fieldindex[i]]
			fieldtype=desiredfieldtype
			fieldnames=desiredfields
		#to store the string containing data format
		fieldtypestring=''
		for i in np.arange(fieldnames.__len__()):
			if fieldtype[i]==np.uint64: fieldtypestring+='u8,'
			elif fieldtype[i]==np.int64: fieldtypestring+='i8,'
			elif fieldtype[i]==np.float64: fieldtypestring+='f8,'

	elif (ibinary==2):
		#load hdf file
		halofile = h5py.File(filename, 'r')
		filenum=int(halofile["File_id"][0])
		numfiles=int(halofile["Num_of_files"][0])
		numhalos=np.uint64(halofile["Num_of_groups"][0])
		numtothalos=np.uint64(halofile["Total_num_of_groups"][0])
		#atime=np.float(halofile.attrs["Time"])
		fieldnames=[str(n) for n in halofile.keys()]
		#clean of header info
		fieldnames.remove("File_id")
		fieldnames.remove("Num_of_files")
		fieldnames.remove("Num_of_groups")
		fieldnames.remove("Total_num_of_groups")
		fieldtype=[halofile[fieldname].dtype for fieldname in fieldnames]
		#if the desiredfields argument is passed only these fieds are loaded
		if (len(desiredfields)>0):
			if (iverbose):print("Loading subset of all fields in property file ", len(desiredfields), " instead of ", len(fieldnames))
			fieldnames=desiredfields
			fieldtype=[halofile[fieldname].dtype for fieldname in fieldnames]
		halofile.close()

	#allocate memory that will store the halo dictionary
	catalog={fieldnames[i]:np.zeros(numtothalos,dtype=fieldtype[i]) for i in range(len(fieldnames))}
	noffset=np.uint64(0)
	for ifile in range(numfiles):
		if (inompi==True): filename=basefilename+".properties"
		else: filename=basefilename+".properties"+"."+str(ifile)
		if (iverbose) : print("reading ",filename)
		if (ibinary==0):
			halofile = open(filename, 'r')
			halofile.readline()
			numhalos=np.uint64(halofile.readline().split()[0])
			halofile.close()
			if (numhalos>0):htemp = np.loadtxt(filename,skiprows=3, usecols=fieldindex, dtype=fieldtypestring, unpack=True, ndmin=1)
		elif(ibinary==1):
			halofile = open(filename, 'rb')
			np.fromfile(halofile,dtype=np.int32,count=2)
			numhalos=np.fromfile(halofile,dtype=np.uint64,count=2)[0]
			#halofile.seek(byteoffset);
			if (numhalos>0):htemp=np.fromfile(halofile, usecols=fieldindex, dtype=fieldtypestring, unpack=True)
			halofile.close()
		elif(ibinary==2):
			#here convert the hdf information into a numpy array
			halofile = h5py.File(filename, 'r')
			numhalos=np.uint64(halofile["Num_of_groups"][0])
			if (numhalos>0):htemp=[np.array(halofile[catvalue]) for catvalue in fieldnames]
			halofile.close()
		#numhalos=len(htemp[0])
		for i in range(len(fieldnames)):
			catvalue=fieldnames[i]
			if (numhalos>0): catalog[catvalue][noffset:noffset+numhalos]=htemp[i]
		noffset+=numhalos
	#if subhalos are written in separate files, then read them too
	if (iseparatesubfiles==1):
		for ifile in range(numfiles):
			if (inompi==True): filename=basefilename+".sublevels"+".properties"
			else: filename=basefilename+".sublevels"+".properties"+"."+str(ifile)
			if (iverbose) : print("reading ",filename)
			if (ibinary==0):
				halofile = open(filename, 'r')
				halofile.readline()
				numhalos=np.uint64(halofile.readline().split()[0])
				halofile.close()
				if (numhalos>0):htemp = np.loadtxt(filename,skiprows=3, usecols=fieldindex, dtype=fieldtypestring, unpack=True, ndmin=1)
			elif(ibinary==1):
				halofile = open(filename, 'rb')
				#halofile.seek(byteoffset);
				np.fromfile(halofile,dtype=np.int32,count=2)
				numhalos=np.fromfile(halofile,dtype=np.uint64,count=2)[0]
				if (numhalos>0):htemp=np.fromfile(halofile, usecols=fieldindex, dtype=fieldtypestring, unpack=True)
				halofile.close()
			elif(ibinary==2):
				halofile = h5py.File(filename, 'r')
				numhalos=np.uint64(halofile["Num_of_groups"][0])
				if (numhalos>0):htemp=[np.array(halofile[catvalue]) for catvalue in fieldnames]
				halofile.close()
			#numhalos=len(htemp[0])
			for i in range(len(fieldnames)):
				catvalue=fieldnames[i]
			if (numhalos>0): catalog[catvalue][noffset:noffset+numhalos]=htemp[i]
			noffset+=numhalos
	#load associated simulation info, time and units
	if (isiminfo):
		siminfoname=basefilename+".siminfo"
		siminfo=open(siminfoname,'r')
		catalog['SimulationInfo']=dict()
		for l in siminfo:
			d=l.strip().split(' : ')
			catalog['SimulationInfo'][d[0]]=float(d[1])
		siminfo.close()
	if (iunitinfo):
		unitinfoname=basefilename+".units"
		unitinfo=open(unitinfoname,'r')
		catalog['UnitInfo']=dict()
		for l in unitinfo:
			d=l.strip().split(' : ')
			catalog['UnitInfo'][d[0]]=float(d[1])
		unitinfo.close()

	if (iverbose): print("done reading properties file ",time.clock()-start)
	return catalog,numtothalos

def ReadPropertyFileMultiWrapper(basefilename,index,halodata,numhalos,atime,ibinary=0,iseparatesubfiles=0,iverbose=0,desiredfields=[]):
	"""
	Wrapper for multithreaded reading
	"""
	#call read routine and store the data
	halodata[index],numhalos[index],atime[index]=ReadPropertyFile(basefilename,ibinary,iseparatesubfiles,iverbose,desiredfields)

def ReadPropertyFileMultiWrapperNamespace(index,basefilename,ns,ibinary=0,iseparatesubfiles=0,iverbose=0,desiredfields=[]):
	#call read routine and store the data
	ns.hdata[index],ns.ndata[index],ns.adata[index]=ReadPropertyFile(basefilename,ibinary,iseparatesubfiles,iverbose,desiredfields)

def ReadHaloMergerTree(treefilename,ibinary=0,iverbose=0,imerit=False,inpart=False):
	"""
	VELOCIraptor/STF merger tree in ascii format contains
	a header with
		number_of_snapshots
		a description of how the tree was built
		total number of halos across all snapshots

	then followed by data
	for each snapshot
		snapshotvalue numhalos
		haloid_1 numprogen_1
		progenid_1
		progenid_2
		...
		progenid_numprogen_1
		haloid_2 numprogen_2
		.
		.
		.
	one can also have an output format that has an additional field for each progenitor, the meritvalue

	"""
	start = time.clock()
	tree=[]
	if (iverbose): print("reading Tree file",treefilename,os.path.isfile(treefilename))
	if (os.path.isfile(treefilename)==False):
		print("Error, file not found")
		return tree
	#if ascii format
	if (ibinary==0):
		treefile = open(treefilename, 'r')
		numsnap=int(treefile.readline())
		treefile.close()
	elif(ibinary==2):
		snaptreelist=open(treefilename,'r')
		numsnap = sum(1 for line in snaptreelist)
		snaptreelist.close()
	else:
		print("Unknown format, returning null")
		numsnap=0
		return tree

	tree=[{"haloID": [], "Num_progen": [], "Progen": []} for i in range(numsnap)]
	if (imerit):
		for i in range(numsnap):
			tree[i]['Merit']=[]
	if (inpart):
		for i in range(numsnap):
			tree[i]['Npart']=[]
			tree[i]['Npart_progen']=[]

	#if ascii format
	if (ibinary==0):
		treefile = open(treefilename, 'r')
		numsnap=int(treefile.readline())
		descrip=treefile.readline().strip()
		tothalos=int(treefile.readline())
		offset=0
		totalnumprogen=0
		for i in range(numsnap):
			[snapval,numhalos]=treefile.readline().strip().split('\t')
			snapval=int(snapval);numhalos=int(numhalos)
			#if really verbose
			if (iverbose==2): print(snapval,numhalos)
			tree[i]["haloID"]=np.zeros(numhalos, dtype=np.int64)
			tree[i]["Num_progen"]=np.zeros(numhalos, dtype=np.uint32)
			tree[i]["Progen"]=[[] for j in range(numhalos)]
			if (imerit): tree[i]["Merit"]=[[] for j in range(numhalos)]
			if (inpart):
				tree[i]["Npart"]=np.zeros(numhalos, dtype=np.uint32)
				tree[i]["Npart_progen"]=[[] for j in range(numhalos)]
			for j in range(numhalos):
				data=treefile.readline().strip().split('\t')
				hid=np.int64(data[0]);nprog=np.uint32(data[1])
				tree[i]["haloID"][j]=hid
				tree[i]["Num_progen"][j]=nprog
				if (inpart):tree[i]["Npart"][j]=np.uint32(data[2])
				totalnumprogen+=nprog
				if (nprog>0):
					tree[i]["Progen"][j]=np.zeros(nprog,dtype=np.int64)
					if (imerit): tree[i]["Merit"][j]=np.zeros(nprog,dtype=np.float32)
					if (inpart): tree[i]["Npart_progen"][j]=np.zeros(nprog,dtype=np.uint32)
					for k in range(nprog):
						data=treefile.readline().strip().split(' ')
						tree[i]["Progen"][j][k]=np.int64(data[0])
						if (imerit):tree[i]["Merit"][j][k]=np.float32(data[1])
						if (inpart):tree[i]["Npart_progen"][j][k]=np.uint32(data[2])

	elif(ibinary==2):

		snaptreelist=open(treefilename,'r')
		#read the first file, get number of snaps from hdf file
		snaptreename = snaptreelist.readline().strip()+".tree"
		treedata=h5py.File(snaptreename,"r")
		numsnaps=treedata.attrs['Number_of_snapshots']
		treedata.close()
		snaptreelist.close()

		snaptreelist=open(treefilename,'r')
		for snap in range(numsnaps):
			snaptreename = snaptreelist.readline().strip()+".tree"
			if (iverbose): print("Reading",snaptreename)
			treedata = h5py.File(snaptreename,"r")

			tree[snap]["haloID"] = np.asarray(treedata["ID"])
			tree[snap]["Num_progen"] = np.asarray(treedata["NumProgen"])
			if(inpart):tree[snap]["Npart"] = np.asarray(treedata["Npart"])

			#See if the dataset exits
			if("ProgenOffsets" in treedata.keys()):

				#Find the indices to split the array
				split = np.add(np.asarray(treedata["ProgenOffsets"]),tree[snap]["Num_progen"],dtype=np.uint64,casting="unsafe")

				#Read in the progenitors, splitting them as reading them in
				tree[snap]["Progen"] = np.split(treedata["Progenitors"][:],split[:-1])

				if(inpart): tree[snap]["Npart_progen"] = np.split(treedata["ProgenNpart"],split[:-1])
				if(imerit): tree[snap]["Merit"] =  np.split(treedata["Merits"],split[:-1])

		snaptreelist.close()
	if (iverbose): print("done reading tree file ",time.clock()-start)
	return tree

def ReadHaloMergerTreeDescendant(treefilename,ireverseorder=True,ibinary=0,iverbose=0,imerit=False,inpart=False):
	"""
	VELOCIraptor/STF descendant based merger tree in ascii format contains
	a header with
		number_of_snapshots
		a description of how the tree was built
		total number of halos across all snapshots

	then followed by data
	for each snapshot
		snapshotvalue numhalos
		haloid_1 numprogen_1
		progenid_1
		progenid_2
		...
		progenid_numprogen_1
		haloid_2 numprogen_2
		.
		.
		.
	one can also have an output format that has an additional field for each progenitor, the meritvalue

	"""
	start = time.clock()
	tree=[]
	if (iverbose): print("reading Tree file",treefilename,os.path.isfile(treefilename))
	if (os.path.isfile(treefilename)==False):
		print("Error, file not found")
		return tree
	#fine out how many snapshots there are
	#if ascii format
	if (ibinary==0):
		if (iverbose): print("Reading ascii input")
		treefile = open(treefilename, 'r')
		numsnap=int(treefile.readline())
		treefile.close()
	#hdf format, input file is a list of filenames
	elif(ibinary==2):
		if (iverbose): print("Reading HDF5 input")
		snaptreelist=open(treefilename,'r')
		numsnap = sum(1 for line in snaptreelist)
		snaptreelist.close()
	else:
		print("Unknown format, returning null")
		numsnap=0
		return tree

	tree=[{"haloID": [], "Num_descen": [], "Descen": [], "Rank": []} for i in range(numsnap)]
	if (imerit):
		for i in range(numsnap):
			tree[i]['Merit']=[]
	if (inpart):
		for i in range(numsnap):
			tree[i]['Npart']=[]
			tree[i]['Npart_descen']=[]

	if (ibinary==0):
		treefile = open(treefilename, 'r')
		numsnap=int(treefile.readline())
		descrip=treefile.readline().strip()
		tothalos=int(treefile.readline())
		offset=0
		totalnumdescen=0
		for i in range(numsnap):
			ii=i
			if (ireverseorder): ii=numsnap-1-i
			[snapval,numhalos]=treefile.readline().strip().split('\t')
			snapval=int(snapval);numhalos=int(numhalos)
			#if really verbose
			if (iverbose==2): print(snapval,numhalos)
			tree[ii]["haloID"]=np.zeros(numhalos, dtype=np.int64)
			tree[ii]["Num_descen"]=np.zeros(numhalos, dtype=np.uint32)
			tree[ii]["Descen"]=[[] for j in range(numhalos)]
			tree[ii]["Rank"]=[[] for j in range(numhalos)]
			if (imerit): tree[ii]["Merit"]=[[] for j in range(numhalos)]
			if (inpart):
				tree[i]["Npart"]=np.zeros(numhalos, dtype=np.uint32)
				tree[ii]["Npart_descen"]=[[] for j in range(numhalos)]
			for j in range(numhalos):
				data=treefile.readline().strip().split('\t')
				hid=np.int64(data[0]);ndescen=np.uint32(data[1])
				tree[ii]["haloID"][j]=hid
				tree[ii]["Num_descen"][j]=ndescen
				if (inpart):tree[ii]["Npart"][j]=np.uint32(data[2])
				totalnumdescen+=ndescen
				if (ndescen>0):
					tree[ii]["Descen"][j]=np.zeros(ndescen,dtype=np.int64)
					tree[ii]["Rank"][j]=np.zeros(ndescen,dtype=np.uint32)
					if (imerit): tree[ii]["Merit"][j]=np.zeros(ndescen,dtype=np.float32)
					if (inpart): tree[ii]["Npart_descen"][j]=np.zeros(ndescen,dtype=np.float32)
					for k in range(ndescen):
						data=treefile.readline().strip().split(' ')
						tree[ii]["Descen"][j][k]=np.int64(data[0])
						tree[ii]["Rank"][j][k]=np.uint32(data[1])
						if (imerit): tree[ii]["Merit"][j][k]=np.float32(data[2])
						if (inpart): tree[ii]["Npart_descen"][j][k]=np.uint32(data[3])

	#hdf format
	elif(ibinary==2):

		snaptreelist=open(treefilename,'r')
		#read the first file, get number of snaps from hdf file
		snaptreename = snaptreelist.readline().strip()+".tree"
		treedata=h5py.File(snaptreename,"r")
		numsnaps=treedata.attrs['Number_of_snapshots']
		treedata.close()
		snaptreelist.close()
		snaptreelist=open(treefilename,'r')
		for snap in range(numsnap):
			snaptreename = snaptreelist.readline().strip()+".tree"
			if (iverbose): print("Reading",snaptreename)
			treedata = h5py.File(snaptreename,"r")
			tree[snap]["haloID"] = np.array(treedata["ID"])
			tree[snap]["Num_descen"] = np.array(treedata["NumDesc"])
			if(inpart):tree[snap]["Npart"] = np.asarray(treedata["Npart"])

			#See if the dataset exits
			if("DescOffsets" in treedata.keys()):

				#Find the indices to split the array
				split = np.add(np.array(treedata["DescOffsets"]), tree[snap]["Num_descen"],dtype=np.uint64,casting="unsafe")

				# Read in the data splitting it up as reading it in
				tree[snap]["Rank"] = np.split(treedata["Ranks"][:],split[:-1])
				tree[snap]["Descen"] = np.split(treedata["Descendants"][:],split[:-1])

				if(inpart): tree[snap]["Npart_progen"] = np.split(treedata["ProgenNpart"][:],split[:-1])
				if(imerit): tree[snap]["Merit"] =  np.split(treedata["Merits"][:],split[:-1])

		snaptreelist.close()

	if (iverbose): print("done reading tree file ",time.clock()-start)
	return tree

def ReadHaloPropertiesAcrossSnapshots(numsnaps,snaplistfname,inputtype,iseperatefiles,iverbose=0,desiredfields=[]):
	"""
	read halo data from snapshots listed in file with snaplistfname file name
	"""
	halodata=[dict() for j in range(numsnaps)]
	ngtot=[0 for j in range(numsnaps)]
	atime=[0 for j in range(numsnaps)]
	start=time.clock()
	print("reading data")
	#if there are a large number of snapshots to read, read in parallel
	#only read in parallel if worthwhile, specifically if large number of snapshots and snapshots are ascii
	iparallel=(numsnaps>20 and inputtype==2)
	if (iparallel):
		#determine maximum number of threads
		nthreads=min(mp.cpu_count(),numsnaps)
		nchunks=int(np.ceil(numsnaps/float(nthreads)))
		print("Using", nthreads,"threads to parse ",numsnaps," snapshots in ",nchunks,"chunks")
		#load file names
		snapnamelist=open(snaplistfname,'r')
		catfilename=["" for j in range(numsnaps)]
		for j in range(numsnaps):
			catfilename[j]=snapnamelist.readline().strip()
		#allocate a manager
		manager = mp.Manager()
		#use manager to specify the dictionary and list that can be accessed by threads
		hdata=manager.list([manager.dict() for j in range(numsnaps)])
		ndata=manager.list([0 for j in range(numsnaps)])
		adata=manager.list([0 for j in range(numsnaps)])
		#now for each chunk run a set of proceses
		for j in range(nchunks):
			offset=j*nthreads
			#if last chunk then must adjust nthreads
			if (j==nchunks-1):
				nthreads=numsnaps-offset
			#when calling a process pass manager based proxies, which then are used to copy data back
			processes=[mp.Process(target=ReadPropertyFileMultiWrapper,args=(catfilename[offset+k],k+offset,hdata,ndata,adata,inputtype,iseperatefiles,iverbose,desiredfields)) for k in range(nthreads)]
			#start each process
			#store the state of each thread, alive or not, and whether it has finished
			activethreads=[[True,False] for k in range(nthreads)]
			count=0
			for p in processes:
				print("reading", catfilename[offset+count])
				p.start()
				#space threads apart (join's time out is 0.25 seconds
				p.join(0.2)
				count+=1
			totactivethreads=nthreads
			while(totactivethreads>0):
				count=0
				for p in processes:
					#join thread and see if still active
					p.join(0.5)
					if (p.is_alive()==False):
						#if thread nolonger active check if its been processed
						if (activethreads[count][1]==False):
							#make deep copy of manager constructed objects that store data
							#halodata[i][offset+count]=deepcopy(hdata[offset+count])
							#try instead init a dictionary
							halodata[offset+count]=dict(hdata[offset+count])
							ngtot[offset+count]=ndata[offset+count]
							atime[offset+count]=adata[offset+count]
							#effectively free the data in manager dictionary
							hdata[offset+count]=[]
							activethreads[count][0]=False
							activethreads[count][1]=True
							totactivethreads-=1
					count+=1
			#terminate threads
			for p in processes:
				p.terminate()

	else:
		snapnamelist=open(snaplistfname,'r')
		for j in range(0,numsnaps):
			catfilename=snapnamelist.readline().strip()
			print("reading ", catfilename)
			halodata[j],ngtot[j],atime[j] = ReadPropertyFile(catfilename,inputtype,iseperatefiles,iverbose,desiredfields)
	print("data read in ",time.clock()-start)
	return halodata,ngtot,atime

def ReadCrossCatalogList(fname,meritlim=0.1,iverbose=0):
	"""
	Reads a cross catalog produced by halomergertree,
	also allows trimming of cross catalog using a higher merit threshold than one used to produce catalog
	"""
	start = time.clock()
	if (iverbose): print("reading cross catalog")
	dfile=open(fname,"r")
	dfile.readline()
	dfile.readline()
	dataline=(dfile.readline().strip()).split('\t')
	ndata=np.int32(dataline[1])
	pdata=CrossCatalogList(ndata)
	for i in range(0,ndata):
		data=(dfile.readline().strip()).split('\t')
		nmatches=np.int32(data[1])
		for j in range(0,nmatches):
			data=(dfile.readline().strip()).split(' ')
			meritval=np.float32(data[1])
			nsharedval=np.float32(data[2])
			if(meritval>meritlim):
				nmatchid=np.int64(data[0])
				pdata.matches[i].append(nmatchid)
				pdata.matches[i].append(meritval)
				pdata.nsharedfrac[i].append(nsharedval)
				pdata.nmatches[i]+=1
	dfile.close()
	if (iverbose): print("done reading cross catalog ",time.clock()-start)
	return pdata

def ReadSimInfo(basefilename):
	"""
	Reads in the information in .siminfo and returns it as a dictionary
	"""

	filename = basefilename + ".siminfo"

	if (os.path.isfile(filename)==False):
		print("file not found")
		return []

	cosmodata = {}
	siminfofile = open(filename,"r")
	line = siminfofile.readline().strip().split(" : ")
	while(line[0]!=""):
		cosmodata[line[0]] = float(line[1])
		line = siminfofile.readline().strip().split(" : ")
	siminfofile.close()
	return cosmodata

def ReadUnitInfo(basefilename):
	"""
	Reads in the information in .units and returns it as a dictionary
	"""

	filename = basefilename + ".units"

	if (os.path.isfile(filename)==False):
		print("file not found")
		return []

	unitdata = {}
	unitsfile = open(filename,"r")
	line = unitsfile.readline().strip().split(" : ")
	while(line[0]!=""):
		unitdata[line[0]] = float(line[1])
		line = unitsfile.readline().strip().split(" : ")
	unitsfile.close()
	return unitdata

def ReadParticleDataFile(basefilename,ibinary=0,iseparatesubfiles=0,iparttypes=0,iverbose=0, binarydtype=np.int64):
	"""
	VELOCIraptor/STF catalog_group, catalog_particles and catalog_parttypes in various formats

	Note that a file will indicate how many files the total output has been split into

	"""
	inompi=True
	if (iverbose): print("reading particle data",basefilename)
	gfilename=basefilename+".catalog_groups"
	pfilename=basefilename+".catalog_particles"
	upfilename=pfilename+".unbound"
	tfilename=basefilename+".catalog_parttypes"
	utfilename=tfilename+".unbound"
	#check for file existence
	if (os.path.isfile(gfilename)==True):
		numfiles=0
	else:
		gfilename+=".0"
		pfilename+=".0"
		upfilename+=".0"
		tfilename+=".0"
		utfilename+=".0"
		inompi=False
		if (os.path.isfile(gfilename)==False):
			print("file not found")
			return []
	byteoffset=0

	#load header information from file to get total number of groups
	#ascii
	if (ibinary==0):
		gfile = open(gfilename, 'r')
		[filenum,numfiles]=gfile.readline().split()
		filenum=int(filenum);numfiles=int(numfiles)
		[numhalos, numtothalos]= gfile.readline().split()
		numhalos=np.uint64(numhalos);numtothalos=np.uint64(numtothalos)
	#binary
	elif (ibinary==1):
		gfile = open(gfilename, 'rb')
		[filenum,numfiles]=np.fromfile(gfile,dtype=np.int32,count=2)
		[numhalos,numtothalos]=np.fromfile(gfile,dtype=np.uint64,count=2)
	#hdf
	elif (ibinary==2):
		gfile = h5py.File(gfilename, 'r')
		filenum=int(gfile["File_id"][0])
		numfiles=int(gfile["Num_of_files"][0])
		numhalos=np.uint64(gfile["Num_of_groups"][0])
		numtothalos=np.uint64(gfile["Total_num_of_groups"][0])
	gfile.close()

	particledata=dict()
	particledata['Npart']=np.zeros(numtothalos,dtype=uint64)
	particledata['Npart_unbound']=np.zeros(numtothalos,dtype=uint64)
	particledata['Particle_IDs']=[[] for i in range(numtothalos)]
	if (iparttypes==1):
		particledata['Particle_Types']=[[] for i in range(numtothalos)]

	#now for all files
	counter=np.uint64(0)
	subfilenames=[""]
	if (iseparatesubfiles==1): subfilenames=["",".sublevels"]
	for ifile in range(numfiles):
		for subname in subfilenames:
			bfname=basefilename+subname
			gfilename=bfname+".catalog_groups"
			pfilename=bfname+".catalog_particles"
			upfilename=pfilename+".unbound"
			tfilename=bfname+".catalog_parttypes"
			utfilename=tfilename+".unbound"
			if (inompi==False):
				gfilename+="."+str(ifile)
				pfilename+="."+str(ifile)
				upfilename+="."+str(ifile)
				tfilename+="."+str(ifile)
				utfilename+="."+str(ifile)
			if (iverbose) : print("reading",bfname,ifile)

			#ascii
			if (ibinary==0):
				gfile = open(gfilename, 'r')
				#read header information
				gfile.readline()
				[numhalos,foo]= gfile.readline().split()
				numhalos=np.uint64(numhalos)
				gfile.close()
				#load data
				gdata=np.loadtxt(gfilename,skiprows=2,dtype=np.uint64)
				numingroup=gdata[:numhalos]
				offset=gdata[int(numhalos):int(2*numhalos)]
				uoffset=gdata[int(2*numhalos):int(3*numhalos)]
				#particle id data
				pfile=open(pfilename, 'r')
				pfile.readline()
				[npart,foo]= pfile.readline().split()
				npart=np.uint64(npart)
				pfile.close()
				piddata=np.loadtxt(pfilename,skiprows=2,dtype=np.int64)
				upfile= open(upfilename, 'r')
				upfile.readline()
				[unpart,foo]= upfile.readline().split()
				unpart=np.uint64(unpart)
				upfile.close()
				upiddata=np.loadtxt(upfilename,skiprows=2,dtype=np.int64)
				if (iparttypes==1):
					#particle id data
					tfile= open(tfilename, 'r')
					tfile.readline()
					[npart,foo]= tfile.readline().split()
					tfile.close()
					tdata=np.loadtxt(tfilename,skiprows=2,dtype=np.uint16)
					utfile= open(utfilename, 'r')
					utfile.readline()
					[unpart,foo]= utfile.readline().split()
					utfile.close()
					utdata=np.loadtxt(utfilename,skiprows=2,dtype=np.uint16)
			#binary
			elif (ibinary==1):
				gfile = open(gfilename, 'rb')
				np.fromfile(gfile,dtype=np.int32,count=2)
				[numhalos,foo]=np.fromfile(gfile,dtype=np.uint64,count=2)
				#need to generalise to
				numingroup=np.fromfile(gfile,dtype=binarydtype ,count=numhalos)
				offset=np.fromfile(gfile,dtype=binarydtype,count=numhalos)
				uoffset=np.fromfile(gfile,dtype=binarydtype,count=numhalos)
				gfile.close()
				pfile = open(pfilename, 'rb')
				np.fromfile(pfile,dtype=np.int32,count=2)
				[npart,foo]=np.fromfile(pfile,dtype=np.uint64,count=2)
				piddata=np.fromfile(pfile,dtype=binarydtype ,count=npart)
				pfile.close()
				upfile = open(upfilename, 'rb')
				np.fromfile(upfile,dtype=np.int32,count=2)
				[unpart,foo]=np.fromfile(upfile,dtype=np.uint64,count=2)
				upiddata=np.fromfile(upfile,dtype=binarydtype ,count=unpart)
				upfile.close()
				if (iparttypes==1):
					tfile = open(tfilename, 'rb')
					np.fromfile(tfile,dtype=np.int32,count=2)
					[npart,foo]=np.fromfile(tfile,dtype=np.uint16,count=2)
					tdata=np.fromfile(tfile,dtype=binarydtype ,count=npart)
					tfile.close()
					utfile = open(utfilename, 'rb')
					np.fromfile(utfile,dtype=np.int32,count=2)
					[unpart,foo]=np.fromfile(utfile,dtype=np.uint16,count=2)
					utdata=np.fromfile(utfile,dtype=binarydtype ,count=unpart)
					utfile.close()
			#hdf
			elif (ibinary==2):
				gfile = h5py.File(gfilename, 'r')
				numhalos=np.uint64(gfile["Num_of_groups"][0])
				numingroup=np.uint64(gfile["Group_Size"])
				offset=np.uint64(gfile["Offset"])
				uoffset=np.uint64(gfile["Offset_unbound"])
				gfile.close()
				pfile = h5py.File(pfilename, 'r')
				upfile = h5py.File(upfilename, 'r')
				piddata=np.int64(pfile["Particle_IDs"])
				upiddata=np.int64(upfile["Particle_IDs"])
				npart=len(piddata)
				unpart=len(upiddata)

				pfile.close()
				upfile.close()
				if (iparttypes==1):
					tfile = h5py.File(tfilename, 'r')
					utfile = h5py.File(utfilename, 'r')
					tdata=np.uint16(pfile["Particle_Types"])
					utdata=np.uint16(upfile["Particle_Types"])
					tfile.close()
					utfile.close()


			#now with data loaded, process it to produce data structure
			particledata['Npart'][counter:counter+numhalos]=numingroup
			unumingroup=np.zeros(numhalos,dtype=uint64)
			for i in range(int(numhalos-1)):
				unumingroup[i]=(uoffset[i+1]-uoffset[i]);
			unumingroup[-1]=(unpart-uoffset[-1])
			particledata['Npart_unbound'][counter:counter+numhalos]=unumingroup
			for i in range(numhalos):
				particledata['Particle_IDs'][int(i+counter)]=np.zeros(numingroup[i],dtype=np.int64)
				particledata['Particle_IDs'][int(i+counter)][:int(numingroup[i]-unumingroup[i])]=piddata[offset[i]:offset[i]+numingroup[i]-unumingroup[i]]
				particledata['Particle_IDs'][int(i+counter)][int(numingroup[i]-unumingroup[i]):numingroup[i]]=upiddata[uoffset[i]:uoffset[i]+unumingroup[i]]
				if (iparttypes==1):
					particledata['Particle_Types'][int(i+counter)]=np.zeros(numingroup[i],dtype=np.int64)
					particledata['Particle_Types'][int(i+counter)][:int(numingroup[i]-unumingroup[i])]=tdata[offset[i]:offset[i]+numingroup[i]-unumingroup[i]]
					particledata['Particle_Types'][int(i+counter)][int(numingroup[i]-unumingroup[i]):numingroup[i]]=utdata[uoffset[i]:uoffset[i]+unumingroup[i]]
			counter+=numhalos

	return particledata

"""
	Routines to build a hierarchy structure (both spatially and temporally)
"""

def BuildHierarchy(halodata,iverbose=0):
	"""
	the halo data stored in a velociraptor .properties file should store the id of its parent halo. Here
	this catalog is used to produce a hierarchy to quickly access the relevant subhaloes of a parent halo.
	#todo this should be deprecated as Hierarchy information is typically already contained in halo information
	"""
	halohierarchy=[]
	start=time.clock()
	if (iverbose): print("setting hierarchy")
	numhalos=len(halodata["npart"])
	subhaloindex=np.where(halodata["hostHaloID"]!=-1)
	lensub=len(subhaloindex[0])
	haloindex=np.where(halodata["hostHaloID"]==-1)
	lenhal=len(haloindex[0])
	halohierarchy=[[] for k in range(numhalos)]
	if (iverbose): print("prelims done ",time.clock()-start)
	for k in range(lenhal):
		halohierarchy[haloindex[0][k]]=np.where(halodata["hostHaloID"]==halodata["ID"][haloindex[0][k]])
	#NOTE: IMPORTANT this is only adding the subsub halos! I need to eventually parse the hierarchy
	#data first to deteremine the depth of the subhalo hierarchy and store how deep an object is in the hierarchy
	#then I can begin adding (sub)subhalos to parent subhalos from the bottom level up
	"""
	for k in range(0,len(halodata["npart"])):
		hid=np.int32(halodata["hostHaloID"][k])
		if (hid>-1 and halohierarchy[k]!=[]):
			halohierarchy[hid]=np.append(np.int32(halohierarchy[hid]),halohierarchy[k])
	"""
	if (iverbose): print("hierarchy set in read in ",time.clock()-start)
	return halohierarchy

def TraceMainProgen(istart,ihalo,numsnaps,numhalos,halodata,tree,TEMPORALHALOIDVAL):
	"""
	Follows a halo along tree to identify main progenitor
	"""
	#start at this snapshot
	k=istart
	#see if halo does not have a tail (descendant set).
	if (halodata[k]['Tail'][ihalo]==0):
		#if halo has not had a tail set the branch needs to be walked along the main branch
		haloid=halodata[k]['ID'][ihalo]
		#only set the head if it has not been set
		#otherwise it should have already been set and just need to store the root head
		if (halodata[k]['Head'][ihalo]==0):
			halodata[k]['Head'][ihalo]=haloid
			halodata[k]['HeadSnap'][ihalo]=k
			halodata[k]['RootHead'][ihalo]=haloid
			halodata[k]['RootHeadSnap'][ihalo]=k
			roothead,rootsnap,rootindex=haloid,k,ihalo
		else:
			roothead=halodata[k]['RootHead'][ihalo]
			rootsnap=halodata[k]['RootHeadSnap'][ihalo]
			rootindex=int(roothead%TEMPORALHALOIDVAL)-1
		#now move along tree first pass to store head and tails and root heads of main branch
		while (True):
			#instead of seraching array make use of the value of the id as it should be in id order
			#wdata=np.where(tree[k]['haloID']==haloid)
			#w2data=np.where(halodata[k]['ID']==haloid)[0][0]
			wdata=w2data=int(haloid%TEMPORALHALOIDVAL)-1
			halodata[k]['Num_progen'][wdata]=tree[k]['Num_progen'][wdata]
			#if no more progenitors, break from search
			#if (tree[k]['Num_progen'][wdata[0][0]]==0 or len(wdata[0])==0):
			if (tree[k]['Num_progen'][wdata]==0):
				#store for current halo its tail and root tail info (also store root tail for root head)
				halodata[k]['Tail'][w2data]=haloid
				halodata[k]['TailSnap'][w2data]=k
				halodata[k]['RootTail'][w2data]=haloid
				halodata[k]['RootTailSnap'][w2data]=k
				#only set the roots tail if it has not been set before (ie: along the main branch of root halo)
				#if it has been set then we are walking along a secondary branch of the root halo's tree
				if (halodata[rootsnap]['RootTail'][rootindex]==0):
					halodata[rootsnap]['RootTail'][rootindex]=haloid
					halodata[rootsnap]['RootTailSnap'][rootindex]=k
				break

			#store main progenitor
			#mainprog=tree[k]['Progen'][wdata[0][0]][0]
			mainprog=tree[k]['Progen'][wdata][0]
			#calculate stepsize based on the halo ids
			stepsize=int(((haloid-haloid%TEMPORALHALOIDVAL)-(mainprog-mainprog%TEMPORALHALOIDVAL))/TEMPORALHALOIDVAL)
			#store tail
			halodata[k]['Tail'][w2data]=mainprog
			halodata[k]['TailSnap'][w2data]=k+stepsize
			k+=stepsize

			#instead of searching array make use of the value of the id as it should be in id order
			#for progid in tree[k-stepsize]['Progen'][wdata[0][0]]:
			#    wdata3=np.where(halodata[k]['ID']==progid)[0][0]
			for progid in tree[k-stepsize]['Progen'][wdata]:
				wdata3=int(progid%TEMPORALHALOIDVAL)-1
				halodata[k]['Head'][wdata3]=haloid
				halodata[k]['HeadSnap'][wdata3]=k-stepsize
				halodata[k]['RootHead'][wdata3]=roothead
				halodata[k]['RootHeadSnap'][wdata3]=rootsnap

			#then store next progenitor
			haloid=mainprog

def TraceMainProgenParallelChunk(istart,ihalochunk,numsnaps,numhalos,halodata,tree,TEMPORALHALOIDVAL):
	"""
	Wrapper to allow for parallelisation
	"""
	for ihalo in ihalochunk:
		TraceMainProgen(istart,ihalo,numsnaps,numhalos,halodata,tree,TEMPORALHALOIDVAL)

def BuildTemporalHeadTail(numsnaps,tree,numhalos,halodata,TEMPORALHALOIDVAL=1000000000000, iverbose=1):
	"""
	Adds for each halo its Head and Tail and stores Roothead and RootTail to the halo
	properties file
	TEMPORALHALOIDVAL is used to parse the halo ids and determine the step size between descendant and progenitor
	"""
	print("Building Temporal catalog with head and tails")
	for k in range(numsnaps):
		halodata[k]['Head']=np.zeros(numhalos[k],dtype=np.int64)
		halodata[k]['Tail']=np.zeros(numhalos[k],dtype=np.int64)
		halodata[k]['HeadSnap']=np.zeros(numhalos[k],dtype=np.int32)
		halodata[k]['TailSnap']=np.zeros(numhalos[k],dtype=np.int32)
		halodata[k]['RootHead']=np.zeros(numhalos[k],dtype=np.int64)
		halodata[k]['RootTail']=np.zeros(numhalos[k],dtype=np.int64)
		halodata[k]['RootHeadSnap']=np.zeros(numhalos[k],dtype=np.int32)
		halodata[k]['RootTailSnap']=np.zeros(numhalos[k],dtype=np.int32)
		halodata[k]['Num_progen']=np.zeros(numhalos[k],dtype=np.uint32)
	#for each snapshot identify halos that have not had their tail set
	#for these halos, the main branch must be walked
	#allocate python manager to wrapper the tree and halo catalog so they can be altered in parallel
	manager=mp.Manager()
	chunksize=5000000 #have each thread handle this many halos at once
	#init to that at this point snapshots should be run in parallel
	if (numhalos[0]>2*chunksize): iparallel=1
	else: iparallel=-1 #no parallel at all
	iparallel=-1

	totstart=time.clock()

	if (iparallel==1):
		#need to copy halodata as this will be altered
		if (iverbose>0): print("copying halo")
		start=time.clock()
		mphalodata=manager.list([manager.dict(halodata[k]) for k in range(numsnaps)])
		if (iverbose>0): print("done",time.clock()-start)

	for istart in range(numsnaps):
		if (iverbose>0): print("Starting from halos at ",istart,"with",numhalos[istart])
		if (numhalos[istart]==0): continue
		#if the number of halos is large then run in parallel
		if (numhalos[istart]>2*chunksize and iparallel==1):
			#determine maximum number of threads
			nthreads=int(min(mp.cpu_count(),ceil(numhalos[istart]/float(chunksize))))
			nchunks=int(np.ceil(numhalos[istart]/float(chunksize)/float(nthreads)))
			if (iverbose>0): print("Using", nthreads,"threads to parse ",numhalos[istart]," halos in ",nchunks,"chunks, each of size", chunksize)
			#now for each chunk run a set of proceses
			for j in range(nchunks):
				start=time.clock()
				offset=j*nthreads*chunksize
				#if last chunk then must adjust nthreads
				if (j==nchunks-1):
					nthreads=int(ceil((numhalos[istart]-offset)/float(chunksize)))

				halochunk=[range(offset+k*chunksize,offset+(k+1)*chunksize) for k in range(nthreads)]
				#adjust last chunk
				if (j==nchunks-1):
					halochunk[-1]=range(offset+(nthreads-1)*chunksize,numhalos[istart])
				#when calling a process pass not just a work queue but the pointers to where data should be stored
				processes=[mp.Process(target=TraceMainProgenParallelChunk,args=(istart,halochunk[k],numsnaps,numhalos,mphalodata,tree,TEMPORALHALOIDVAL)) for k in range(nthreads)]
				count=0
				for p in processes:
					print(count+offset,k,min(halochunk[count]),max(halochunk[count]))
					p.start()
					count+=1
				for p in processes:
					#join thread and see if still active
					p.join()
				if (iverbose>1): print((offset+j*nthreads*chunksize)/float(numhalos[istart])," done in",time.clock()-start)
		#otherwise just single
		else :
			#if first time entering non parallel section copy data back from parallel manager based structure to original data structure
			#as parallel structures have been updated
			if (iparallel==1):
				#tree=[dict(mptree[k]) for k in range(numsnaps)]
				halodata=[dict(mphalodata[k]) for k in range(numsnaps)]
				#set the iparallel flag to 0 so that all subsequent snapshots (which should have fewer objects) not run in parallel
				#this is principly to minimize the amount of copying between manager based parallel structures and the halo/tree catalogs
				iparallel=0
			start=time.clock()
			chunksize=max(int(0.10*numhalos[istart]),10)
			for j in range(numhalos[istart]):
				#start at this snapshot
				#start=time.clock()
				TraceMainProgen(istart,j,numsnaps,numhalos,halodata,tree,TEMPORALHALOIDVAL)
				if (j%chunksize==0 and j>0):
					if (iverbose>1): print("done", j/float(numhalos[istart]), "in", time.clock()-start)
					start=time.clock()
	if (iverbose>0): print("done with first bit")
	#now have walked all the main branches and set the root head, head and tail values
	#and can set the root tail of all halos. Start at end of the tree and move in reverse setting the root tail
	#of a halo's head so long as that halo's tail is the current halo (main branch)
	for istart in range(numsnaps-1,-1,-1):
		for j in range(numhalos[istart]):
			#if a halo's root tail is itself then start moving up its along to its head (if its head is not itself as well
			k=istart
			#rootheadid,rootheadsnap=halodata[k]['RootHead'][j],halodata[k]['RootHeadSnap'][j]
			roottailid,roottailsnap=halodata[k]['RootTail'][j],halodata[k]['RootTailSnap'][j]
			headid,headsnap=halodata[k]['Head'][j],halodata[k]['HeadSnap'][j]
			if (roottailid==halodata[k]['ID'][j] and headid!=halodata[k]['ID'][j]):
				#headindex=np.where(halodata[headsnap]['ID']==headid)[0][0]
				headindex=int(headid%TEMPORALHALOIDVAL)-1
				headtailid,headtailsnap=halodata[headsnap]['Tail'][headindex],halodata[headsnap]['TailSnap'][headindex]
				haloid=halodata[k]['ID'][j]
				#only proceed in setting root tails of a head who's tail is the same as halo (main branch) till we reach a halo who is its own head
				while (headtailid==haloid and headid!=haloid):
					#set root tails
					halodata[headsnap]['RootTail'][headindex]=roottailid
					halodata[headsnap]['RootTailSnap'][headindex]=roottailsnap
					#move to next head
					haloid=halodata[headsnap]['ID'][headindex]
					#haloindex=np.where(halodata[headsnap]['ID']==haloid)[0][0]
					haloindex=int(haloid%TEMPORALHALOIDVAL)-1
					halosnap=headsnap
					headid,headsnap=halodata[halosnap]['Head'][haloindex],halodata[halosnap]['HeadSnap'][haloindex]
					headindex=int(headid%TEMPORALHALOIDVAL)-1
					#store the tail of the next head
					headtailid,headtailsnap=halodata[headsnap]['Tail'][headindex],halodata[headsnap]['TailSnap'][headindex]
	print("Done building", time.clock()-totstart)


def TraceMainDescendant(istart,ihalo,numsnaps,numhalos,halodata,tree,TEMPORALHALOIDVAL,ireverseorder=False):
	"""
	Follows a halo along descendant tree to root tails
	if reverse order than late times start at 0 and as one moves up in index
	one moves backwards in time
	"""

	#start at this snapshot
	halosnap=istart

	#see if halo does not have a Head set
	if (halodata[halosnap]['Head'][ihalo]==0):
		#if halo has not had a Head set the branch needs to be walked along the main branch
		haloid=halodata[halosnap]['ID'][ihalo]
		#only set the Root Tail if it has not been set. Here if halo has not had
		#tail set, then must be the the first progenitor
		#otherwise it should have already been set and just need to store the root tail
		if (halodata[halosnap]['Tail'][ihalo]==0):
			halodata[halosnap]['Tail'][ihalo]=haloid
			halodata[halosnap]['TailSnap'][ihalo]=halosnap
			halodata[halosnap]['RootTail'][ihalo]=haloid
			halodata[halosnap]['RootTailSnap'][ihalo]=halosnap
			roottail,rootsnap,rootindex=haloid,halosnap,ihalo
		else:
			roottail=halodata[halosnap]['RootTail'][ihalo]
			rootsnap=halodata[halosnap]['RootTailSnap'][ihalo]
			rootindex=int(roottail%TEMPORALHALOIDVAL)-1
		#now move along tree first pass to store head and tails and root tails of main branch
		while (True):
			#ids contain index information
			haloindex=int(haloid%TEMPORALHALOIDVAL)-1

			halodata[halosnap]['Num_descen'][haloindex]=tree[halosnap]['Num_descen'][haloindex]
			#if no more descendants, break from search
			if (halodata[halosnap]['Num_descen'][haloindex]==0):
				#store for current halo its tail and root tail info (also store root tail for root head)
				halodata[halosnap]['Head'][haloindex]=haloid
				halodata[halosnap]['HeadSnap'][haloindex]=halosnap
				halodata[halosnap]['RootHead'][haloindex]=haloid
				halodata[halosnap]['RootHeadSnap'][haloindex]=halosnap
				rootheadid,rootheadsnap,rootheadindex=haloid,halosnap,haloindex
				#only set the roots head of the root tail
				#if it has not been set before (ie: along the main branch of root halo)
				if (halodata[rootsnap]['RootHead'][rootindex]==0):
					halosnap,haloindex,haloid=rootsnap,rootindex,roottail
					#set the root head of the main branch
					while(True):
						halodata[halosnap]['RootHead'][haloindex]=rootheadid
						halodata[halosnap]['RootHeadSnap'][haloindex]=rootheadsnap
						descen=halodata[halosnap]['Head'][haloindex]
						descenindex=int(descen%TEMPORALHALOIDVAL)-1
						descensnap=int(((descen-descen%TEMPORALHALOIDVAL))/TEMPORALHALOIDVAL)
						if (ireverseorder):
							descensnap=numsnaps-1-descensnap
						if (haloid==descen):
							break
						halosnap,haloindex,haloid=descensnap,descenindex,descen
				break
			#now store the rank of the of the descandant.
			descenrank=tree[halosnap]['Rank'][haloindex][0]
			halodata[halosnap]['HeadRank'][haloindex]=descenrank
			#as we are only moving along main branches stop if object is rank is not 0
			if (descenrank>0):
				break
			#otherwise, get the descendant
			#store main progenitor
			maindescen=tree[halosnap]['Descen'][haloindex][0]
			maindescenindex=int(maindescen%TEMPORALHALOIDVAL)-1
			maindescensnap=int(((maindescen-maindescen%TEMPORALHALOIDVAL))/TEMPORALHALOIDVAL)
			#if reverse order, then higher snap values correspond to lower index
			if (ireverseorder):
				maindescensnap=numsnaps-1-maindescensnap
			#calculate stepsize in time based on the halo ids
			stepsize=maindescensnap-halosnap

			#store descendant
			halodata[halosnap]['Head'][haloindex]=maindescen
			halodata[halosnap]['HeadSnap'][haloindex]=maindescensnap

			#and update the root tails of the object
			halodata[maindescensnap]['Tail'][maindescenindex]=haloid
			halodata[maindescensnap]['TailSnap'][maindescenindex]=halosnap
			halodata[maindescensnap]['RootTail'][maindescenindex]=roottail
			halodata[maindescensnap]['RootTailSnap'][maindescenindex]=rootsnap
			halodata[maindescensnap]['Num_progen'][maindescenindex]+=1

			#then move to the next descendant
			haloid=maindescen
			halosnap=maindescensnap

def TraceMainDescendantParallelChunk(istart,ihalochunk,numsnaps,numhalos,halodata,tree,TEMPORALHALOIDVAL,ireverseorder):
	for ihalo in ihalochunk:
		TraceMainDescendant(istart,ihalo,numsnaps,numhalos,halodata,tree,TEMPORALHALOIDVAL,ireverseorder)

def BuildTemporalHeadTailDescendant(numsnaps,tree,numhalos,halodata,TEMPORALHALOIDVAL=1000000000000, ireverseorder=False, iverbose=1):
	"""
	Adds for each halo its Head and Tail and stores Roothead and RootTail to the halo
	properties file
	TEMPORALHALOIDVAL is used to parse the halo ids and determine the step size between descendant and progenitor
	"""
	print("Building Temporal catalog with head and tails using a descendant tree")
	for k in range(numsnaps):
		halodata[k]['Head']=np.zeros(numhalos[k],dtype=np.int64)
		halodata[k]['Tail']=np.zeros(numhalos[k],dtype=np.int64)
		halodata[k]['HeadSnap']=np.zeros(numhalos[k],dtype=np.int32)
		halodata[k]['TailSnap']=np.zeros(numhalos[k],dtype=np.int32)
		halodata[k]['RootHead']=np.zeros(numhalos[k],dtype=np.int64)
		halodata[k]['RootTail']=np.zeros(numhalos[k],dtype=np.int64)
		halodata[k]['RootHeadSnap']=np.zeros(numhalos[k],dtype=np.int32)
		halodata[k]['RootTailSnap']=np.zeros(numhalos[k],dtype=np.int32)
		halodata[k]['HeadRank']=np.zeros(numhalos[k],dtype=np.int64)
		halodata[k]['Num_descen']=np.zeros(numhalos[k],dtype=np.uint32)
		halodata[k]['Num_progen']=np.zeros(numhalos[k],dtype=np.uint32)
	#for each snapshot identify halos that have not had their tail set
	#for these halos, the main branch must be walked
	#allocate python manager to wrapper the tree and halo catalog so they can be altered in parallel
	manager=mp.Manager()
	chunksize=5000000 #have each thread handle this many halos at once
	#init to that at this point snapshots should be run in parallel
	if (numhalos[0]>2*chunksize): iparallel=1
	else: iparallel=-1 #no parallel at all
	iparallel=-1

	totstart=time.clock()

	if (ireverseorder):
		snaplist=range(numsnaps-1,-1,-1)
	else:
		snaplist=range(numsnaps)

	if (iparallel==1):
		#need to copy halodata as this will be altered
		if (iverbose>0): print("copying halo")
		start=time.clock()
		mphalodata=manager.list([manager.dict(halodata[k]) for k in range(numsnaps)])
		if (iverbose>0): print("done",time.clock()-start)

	for istart in snaplist:
		if (iverbose>0): print("Starting from halos at ",istart,"with",numhalos[istart])
		if (numhalos[istart]==0): continue
		#if the number of halos is large then run in parallel
		if (numhalos[istart]>2*chunksize and iparallel==1):
			#determine maximum number of threads
			nthreads=int(min(mp.cpu_count(),ceil(numhalos[istart]/float(chunksize))))
			nchunks=int(np.ceil(numhalos[istart]/float(chunksize)/float(nthreads)))
			if (iverbose>0): print("Using", nthreads,"threads to parse ",numhalos[istart]," halos in ",nchunks,"chunks, each of size", chunksize)
			#now for each chunk run a set of proceses
			for j in range(nchunks):
				start=time.clock()
				offset=j*nthreads*chunksize
				#if last chunk then must adjust nthreads
				if (j==nchunks-1):
					nthreads=int(ceil((numhalos[istart]-offset)/float(chunksize)))

				halochunk=[range(offset+k*chunksize,offset+(k+1)*chunksize) for k in range(nthreads)]
				#adjust last chunk
				if (j==nchunks-1):
					halochunk[-1]=range(offset+(nthreads-1)*chunksize,numhalos[istart])
				#when calling a process pass not just a work queue but the pointers to where data should be stored
				processes=[mp.Process(target=TraceMainDescendantParallelChunk,args=(istart,halochunk[k],numsnaps,numhalos,mphalodata,tree,TEMPORALHALOIDVAL,ireverseorder)) for k in range(nthreads)]
				count=0
				for p in processes:
					print(count+offset,k,min(halochunk[count]),max(halochunk[count]))
					p.start()
					count+=1
				for p in processes:
					#join thread and see if still active
					p.join()
				if (iverbose>1): print((offset+j*nthreads*chunksize)/float(numhalos[istart])," done in",time.clock()-start)
		#otherwise just single
		else :
			#if first time entering non parallel section copy data back from parallel manager based structure to original data structure
			#as parallel structures have been updated
			if (iparallel==1):
				#tree=[dict(mptree[k]) for k in range(numsnaps)]
				halodata=[dict(mphalodata[k]) for k in range(numsnaps)]
				#set the iparallel flag to 0 so that all subsequent snapshots (which should have fewer objects) not run in parallel
				#this is principly to minimize the amount of copying between manager based parallel structures and the halo/tree catalogs
				iparallel=0
			start=time.clock()
			chunksize=max(int(0.10*numhalos[istart]),10)
			for j in range(numhalos[istart]):
				#start at this snapshot
				#start=time.clock()
				TraceMainDescendant(istart,j,numsnaps,numhalos,halodata,tree,TEMPORALHALOIDVAL,ireverseorder)
				if (j%chunksize==0 and j>0):
					if (iverbose>1): print("done", j/float(numhalos[istart]), "in", time.clock()-start)
					start=time.clock()
	if (iverbose>0): print("done with first bit, setting the main branches walking forward in time")
	#now have walked all the main branches and set the root tail, head and tail values
	#in case halo data is with late times at beginning need to process items in reverse
	if (ireverseorder):
		snaplist=range(numsnaps)
	else:
		snaplist=range(numsnaps-1,-1,-1)

	for istart in snaplist:
		#identify all haloes which are not primary progenitors of their descendants, having a descendant rank >0
		wdata=np.where(halodata[istart]['HeadRank']>0)
		#sort this list based on descendant ranking
		sortedranking=np.argsort(halodata[istart]['HeadRank'][wdata])
		nrankedhalos=len(wdata[0])
		rankedhalos=halodata[istart]['ID'][wdata[0][sortedranking]]
		#for each of these haloes, set the head and use the root head information and root snap and set all the information
		#long its branch
		for ihalo in rankedhalos:
			haloid=ihalo
			haloindex=int(haloid%TEMPORALHALOIDVAL)-1
			halosnap=istart
			#now set the head of these objects
			maindescen=tree[halosnap]['Descen'][haloindex][0]
			maindescenindex=int(maindescen%TEMPORALHALOIDVAL)-1
			if (ireverseorder):
				maindescensnap=numsnaps-1-int((maindescen-maindescen%TEMPORALHALOIDVAL)/TEMPORALHALOIDVAL)
			else:
				maindescensnap=int((maindescen-maindescen%TEMPORALHALOIDVAL)/TEMPORALHALOIDVAL)
			#increase the number of progenitors of this descendant
			halodata[halosnap]['Head'][haloindex]=maindescen
			halodata[halosnap]['HeadSnap'][haloindex]=maindescensnap
			halodata[maindescensnap]['Num_progen'][maindescenindex]+=1
			#store the root head
			roothead=halodata[maindescensnap]['RootHead'][maindescenindex]
			rootsnap=halodata[maindescensnap]['RootHeadSnap'][maindescenindex]
			#now set the root head for all the progenitors of this object
			while (True):
				halodata[halosnap]['RootHead'][haloindex]=roothead
				halodata[halosnap]['RootHeadSnap'][haloindex]=rootsnap
				if (haloid==halodata[halosnap]['Tail'][haloindex]):
					break
				haloid=halodata[halosnap]['Tail'][haloindex]
				halosnap=halodata[halosnap]['TailSnap'][haloindex]
				haloindex=int(haloid%TEMPORALHALOIDVAL)-1
	print("Done building", time.clock()-totstart)

def GetProgenLength(halodata,haloindex,halosnap,haloid,atime,TEMPORALHALOIDVAL,endreftime=-1):
	"""
	Get the length of a halo's progenitors
	"""
	proglen=1
	progid=halodata[halosnap]["Tail"][haloindex]
	progsnap=halodata[halosnap]["TailSnap"][haloindex]
	progindex=int(progid%TEMPORALHALOIDVAL-1)
	while (progid!=haloid):
		if (atime[progsnap]<=endreftime):break
		proglen+=1
		haloid=progid
		halosnap=progsnap
		haloindex=progindex
		progid=halodata[halosnap]["Tail"][haloindex]
		progsnap=halodata[halosnap]["TailSnap"][haloindex]
		progindex=int(progid%TEMPORALHALOIDVAL-1)
	return proglen

def IdentifyMergers(numsnaps,tree,numhalos,halodata,boxsize,hval,atime,MERGERMLIM=0.1,RADINFAC=1.2,RADOUTFAC=1.5,NPARTCUT=100, TEMPORALHALOIDVAL=1000000000000, iverbose=1,pos_tree=[]):
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
		if (iverbose): print("tree build")
		for j in range(numsnaps):
			if (numhalos[j]>0):
				boxval=boxsize*atime[j]/hval
				pos[j]=np.transpose(np.asarray([halodata[j]["Xc"],halodata[j]["Yc"],halodata[j]["Zc"]]))
				pos_tree[j]=spatial.cKDTree(pos[j],boxsize=boxval)
		if (iverbose): print("done ",time.clock()-start)
	#else assume tree has been passed
	for j in range(numsnaps):
		if (numhalos[j]==0): continue
		#at snapshot look at all haloes that have not had a major merger set
		#note that only care about objects with certain number of particles
		partcutwdata=np.where(halodata[j]["npart"]>=NPARTCUT)
		mergercut=np.where(halodata[j]["LastMergerRatio"][partcutwdata]<0)
		hids=np.asarray(halodata[j]["ID"][partcutwdata][mergercut],dtype=np.uint64)
		start=time.clock()
		if (iverbose):print("Processing ", len(hids))
		if (len(hids)==0):continue

		for hidval in hids:
			#now for each object get the main progenitor
			haloid=np.uint64(hidval)
			haloindex=int(haloid%TEMPORALHALOIDVAL-1)
			halosnap=j
			originalhaloid=haloid
			progid=halodata[halosnap]["Tail"][haloindex]
			progsnap=halodata[halosnap]["TailSnap"][haloindex]
			progindex=int(progid%TEMPORALHALOIDVAL-1)
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
										nextindexval=int(nextidval%TEMPORALHALOIDVAL-1)
										starthaloid=nextidval
										starthalosnap=nextsnapval
										starthaloindex=nextindexval

										nextidval=halodata[startmergersnap]["Tail"][startmergerindex]
										nextsnapval=halodata[startmergersnap]["TailSnap"][startmergerindex]
										nextindexval=int(nextidval%TEMPORALHALOIDVAL-1)
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
										mergingnextindex=int(mergingnextid%TEMPORALHALOIDVAL-1)
										merginghaloindex=mergingnextindex
										merginghaloid=mergingnextid
										merginghalosnap=mergingnextsnap

				#move to next step
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
						nextindex=int(nextid%TEMPORALHALOIDVAL-1)
						oldhaloid=currentid
						currentsnap=nextsnap
						currentid=nextid
						currentindex=nextindex
					break
				haloid=progid
				haloindex=progindex
				halosnap=progsnap
				progid=halodata[halosnap]["Tail"][haloindex]
				progsnap=halodata[halosnap]["TailSnap"][haloindex]
				progindex=int(progid%TEMPORALHALOIDVAL-1)
				numprog=tree[halosnap]["Num_progen"][haloindex]
				#if at end of line then move up and set last major merger to 0
		if (iverbose): print("Done snap",j,time.clock()-start)

def generate_sublinks(numhalos,halodata):
	"""
	generate sublinks for specific time slice
	"""
	if (numhalos==0):
		return
	halos=np.where(halodata['hostHaloID']==-1)[0]
	for ihalo in halos:
		haloid=halodata['ID'][ihalo]
		halodata['PreviousSubhalo'][ihalo]=haloid
		w=np.where((halodata["hostHaloID"]==haloid))[0]
		if (len(w)>0):
			halodata['NextSubhalo'][ihalo]=halodata['ID'][w[0]]
			halodata['PreviousSubhalo'][w[0]]=halodata['ID'][w[0]]
			for isub in range(len(w)-1):
				subid=halodata['ID'][w[isub]]
				nextsubid=halodata['ID'][w[isub+1]]
				halodata['NextSubhalo'][w[isub]]=nextsubid
				halodata['PreviousSubhalo'][w[isub+1]]=subid
			halodata['NextSubhalo'][w[-1]]=halodata['ID'][w[-1]]
		else:
			halodata['NextSubhalo'][ihalo]=haloid

def GenerateSubhaloLinks(numsnaps,numhalos,halodata,TEMPORALHALOIDVAL=1000000000000, iverbose=0, iparallel=0):
	"""
	This code generates a quick way of moving across a halo's subhalo list

	The code is passed
	- the number of snapshots,
	- an array of the number of haloes per snapshot,
	- the halodata dictionary structure which must contain the halo merger tree based keys, Head, RootHead, etc, and mass, phase-space positions of haloes,
	and other desired properties
	"""
	for j in range(numsnaps):
		#store id and snap and mass of last major merger and while we're at it, store number of major mergers
		halodata[j]["NextSubhalo"]=np.zeros(numhalos[j],dtype=np.int64)
		halodata[j]["PreviousSubhalo"]=np.zeros(numhalos[j],dtype=np.int64)
	#iterate over all host halos and set their subhalo links
	start=time.clock()
	nthreads=1
	if (iparallel):
		manager=mp.Manager()
		nthreads=int(min(mp.cpu_count(),numsnaps))
		print("Number of threads is ",nthreads)
	for j in range(0,numsnaps,nthreads):
		start2=time.clock()
		if (iparallel):
			activenthreads=nthreads
			if (numsnaps-1-j<activenthreads): activenthreads=numsnaps-1-j
			processes=[mp.Process(target=generate_sublinks,args=(numhalos[j+k],halodata[j+k])) for k in range(activenthreads)]
			for p in processes:
				p.start()
			for p in processes:
				p.join()
			if (iverbose): print("Done snaps",j,"to",j+nthreads,time.clock()-start2)
		else:
			generate_sublinks(numhalos[j],halodata[j])
			if (iverbose): print("Done snap",j,time.clock()-start2)
	print("Done subhalolinks ",time.clock()-start)

def GenerateProgenitorLinks(numsnaps,numhalos,halodata,nsnapsearch=4,TEMPORALHALOIDVAL=1000000000000, iverbose=1):
	"""
	This code generates a quick way of moving across a halo's progenitor list storing a the next/previous progenitor

	The code is passed
	- the number of snapshots,
	- an array of the number of haloes per snapshot,
	- the halodata dictionary structure which must contain the halo merger tree based keys, Head, RootHead, etc, and mass, phase-space positions of haloes,
	and other desired properties
	"""
	for j in range(numsnaps):
		#store id and snap and mass of last major merger and while we're at it, store number of major mergers
		halodata[j]["NextProgenitor"]=np.ones(numhalos[j],dtype=np.int64)*-1
		halodata[j]["PreviousProgenitor"]=np.ones(numhalos[j],dtype=np.int64)*-1
	#move backward in time and identify all unique heads
	start=time.clock()
	for j in range(1,numsnaps):
		start2=time.clock()
		if (numhalos[j]==0): continue
		#find all unique heads
		heads=np.unique(np.array(np.int64(halodata[j]['Head'])))
		#for these heads identify all halos with this head
		for ihead in heads:
			currenttails=deque()
			for k in range(j,j+nsnapsearch):
				w=np.where(halodata[k]['Head']==ihead)
				if (len(w[0])>0):
					currenttails.extend(np.nditer(np.int64(halodata[k]["ID"][w])))
			if (len(currenttails)==0):
				continue
			haloid=currenttails[0]
			haloindex=int(haloid%TEMPORALHALOIDVAL-1)
			halosnap=numsnaps-1-(haloid-int(haloid%TEMPORALHALOIDVAL))/TEMPORALHALOIDVAL
			halodata[halosnap]['PreviousProgenitor'][haloindex]=haloid
			for itail in range(len(currenttails)-1):
				haloid=currenttails[itail]
				haloindex=int(haloid%TEMPORALHALOIDVAL-1)
				halosnap=numsnaps-1-(haloid-int(haloid%TEMPORALHALOIDVAL))/TEMPORALHALOIDVAL
				haloindex=int(currenttails[itail]%TEMPORALHALOIDVAL-1)
				nexthaloid=currenttails[itail+1]
				nexthaloindex=int(nexthaloid%TEMPORALHALOIDVAL-1)
				nexthalosnap=numsnaps-1-(nexthaloid-int(nexthaloid%TEMPORALHALOIDVAL))/TEMPORALHALOIDVAL
				halodata[halosnap]['NextProgenitor'][haloindex]=nexthaloid
				halodata[nexthalosnap]['PreviousProgenitor'][nexthaloindex]=haloid
			haloid=currenttails[-1]
			haloindex=int(haloid%TEMPORALHALOIDVAL-1)
			halosnap=numsnaps-1-(haloid-int(haloid%TEMPORALHALOIDVAL))/TEMPORALHALOIDVAL
			halodata[halosnap]['NextProgenitor'][haloindex]=haloid
		if (iverbose): print("Done snap",j,time.clock()-start2)
	print("Done progenitor links ",time.clock()-start)

def SetForestID(numsnaps,halodata,rootheadid,ForestID,AllRootHead,
	TEMPORALHALOIDVAL = 1000000000000,searchSnapLim = 5, ireversesnaporder=True):
	"""
	Sets the forest id of halos using a roothead as a start point.
	Given an initial root head and end snapshot,
	First append the roothead to the AllRootHead list.
	search all previous snapshots for any haloes that share the same roothead.
	Also at each snapshot, find all subhaloes of all haloes sharing the same
	root head
	if the roothead of a subhalo is not present in the AllRootHead list
	then recursively call SetForestID with this subhalo's root head as start point
	if a subhalo's current host is not within the tree defined by rootheadid
	then recursively call SetForestID with this host's root head as start point

	Parameters
	----------
	numsnaps : numpy.int32
		the number of snapshots
	halodata : dict
		the halodata dictionary structure which must contain the halo merger tree based keys (Head, RootHead), etc.
	rootheadid : numpy.int64
		the rootheadid of the tree that will be explored and have its forestID set
	AllRootHead : list
		a list that stores the current set of rootheadid values that have been searched

	Optional Parameters
	-------------------
	TEMPORALHALOIDVAL : numpy.int64
		Temporal ID value that makes Halo IDs temporally unique, adding a snapshot num* this value.
		Allows one to quickly parse a Halo ID to determine the snapshot it exists at and its index.
	searchSnapLim : numpy.int32
		Maximum number of snapshots to keep searching if no new halos are identified as beloning to
		a rootheadid's tree, moving backwards in time
	ireversesnaporder : bool
		Whether dictionary data has late times starting at 0 (True, default) or at end of dictionary (False)

	Returns
	-------
	AllRootHead : list
		Updated list
	halodata : dict
		Updated halo data

	"""


	if (ireversesnaporder): endSnap = numsnaps-int(rootheadid/TEMPORALHALOIDVAL)-1
	else : endSnap = int(rootheadid/TEMPORALHALOIDVAL)
	rootheadindex=int(rootheadid%TEMPORALHALOIDVAL-1)

	AllRootHead.append(rootheadid)

	#set the forest level of this searcheed
	#if this object is a host at final snap then set the forest level to 0
	#otherwise set the ForestLevel to 1
	ForestLevel=1*(halodata[endSnap]["hostHaloID"][rootheadindex]!=-1)

	#Indicator for amount of snapshots searcheed
	iSearchSnap = 0

	#set the direction of how the data will be processed
	if (ireversesnaporder): snaplist=np.arange(endSnap,numsnaps,dtype=np.int32)
	else : snaplist=np.arange(endsnap,-1,-1)
	for snap in snaplist:
		#Find which halos at this snapshot point to the RootDescedant
		sel = np.where(halodata[snap]["RootHead"]==rootheadid)[0]

		#keep track of how many snapshots there have been where there is nothing in the tree
		if(sel.size==0):
			iSearchSnap+=1
		if(iSearchSnap==searchSnapLim): break
		else: iSearchSnap = 0

		# Set all the halos within this tree within this snapshot to this forest ID
		halodata[snap]["ForestID"][sel] = ForestID
		halodata[snap]["ForestLevel"][sel] = ForestLevel

		#Lets find which halos are subhalos of the halos within the tree defined by
		#halos with the same rootheadid
		subHaloIndxs = np.where(np.in1d(halodata[snap]["hostHaloID"],halodata[snap]["ID"][sel]))[0]
		#Lets loop over all the subhalos within this selection, which contains
		#all subhalos of any host halos within the tree defined by rootheadid
		for subHaloIndx in subHaloIndxs:
			#See if this tree has already been set
			if(halodata[snap]["RootHead"][subHaloIndx] not in AllRootHead):
				#Lets walk the subhalo's tree setting the forest ID
				AllRootHead,halodata = SetForestID(numsnaps,halodata,halodata[snap]["RootHead"][subHaloIndx],ForestID,AllRootHead)

		#Extract the hosts of all subhalos in this selection that are not already in the tree defined by rootheadid
		treeSubhaloSel = (halodata[snap]["hostHaloID"][sel]!=-1) & (np.invert(np.in1d(halodata[snap]["hostHaloID"][sel],halodata[snap]["ID"][sel])))
		#Get the index of these hosts that lie outside the tree
		hostIndxs = np.unique(halodata[snap]["hostHaloID"][sel][treeSubhaloSel]%TEMPORALHALOIDVAL-1).astype(int)
		#Loop over all the index for the host halos
		for hostIndx in hostIndxs:
			#See if this tree has already been set
			if(halodata[snap]["RootHead"][hostIndx] not in AllRootHead):
				#Lets walk the hosts tree setting the forrest ID
				AllRootHead,halodata = SetForestID(numsnaps,halodata,halodata[snap]["RootHead"][hostIndx],ForestID,AllRootHead)

	return AllRootHead,halodata

def GenerateForest(numsnaps,numhalos,halodata,cosmo,atime,
	TEMPORALHALOIDVAL=1000000000000, iverbose=1, interactiontime=2, ispatialintflag=False, pos_tree=[]):
	"""
	This code traces all root heads back in time identifying all interacting haloes and bundles them together into the same forest id
	The idea is to have in the halodata dictionary an associated unique forest id for all related (sub)haloes. The code also allows
	for some cleaning of the forest, specifically if a (sub)halo is only interacting for some small fraction of time, then it is not
	assigned to the forest. This can limit the size of a forest, which could otherwise become the entire halo catalog.

	Parameters
	----------
	numsnaps : numpy.int32
		the number of snapshots
	numhalos : array
		array of the number of haloes per snapshot.
	halodata : dict
		the halodata dictionary structure which must contain the halo merger tree based keys (Head, RootHead), etc.
	cosmo : dict
		dictionary which has cosmological information such as box size, hval, Omega_m
	atime : array
		an array of scale factors

	Optional Parameters
	-------------------
	TEMPORALHALOIDVAL : numpy.int64
		Temporal ID value that makes Halo IDs temporally unique, adding a snapshot num* this value.
		Allows one to quickly parse a Halo ID to determine the snapshot it exists at and its index.
	iverbose : int
		verbosity of function (0, minimal, 1, verbose, 2 chatterbox)
	interactiontime : int
		Optional functionality not implemented yet. Allows forest to be split if connections do not span
		more than this number of snapshots
	ispatialintflag : bool
		Flag indicating whether spatial information should be used to join forests.
	pos_tree : scikit.spatial.cKDTree
		Optional functionality not implemented yet. Allows forests to be joined if haloes
		are spatially close.

	Returns
	-------
	ForestSize : numpy.array
		Update the halodata dictionary with ForestID information and also returns the size of
		the forests

	"""

	#initialize the dictionaries
	for j in range(numsnaps):
		#store id and snap and mass of last major merger and while we're at it, store number of major mergers
		halodata[j]["ForestID"]=np.ones(numhalos[j],dtype=np.int64)*-1
		halodata[j]["ForestLevel"]=np.ones(numhalos[j],dtype=np.int32)*-1
	#built KD tree to quickly search for near neighbours. only build if not passed.
	if (ispatialintflag):
		start=time.clock()
		boxsize=cosmo['BoxSize']
		hval=cosmo['Hubble_param']
		if (len(pos_tree)==0):
			pos=[[]for j in range(numsnaps)]
			pos_tree=[[]for j in range(numsnaps)]
			start=time.clock()
			if (iverbose): print("KD tree build")
			for j in range(numsnaps):
				if (numhalos[j]>0):
					boxval=boxsize*atime[j]/hval
					pos[j]=np.transpose(np.asarray([halodata[j]["Xc"],halodata[j]["Yc"],halodata[j]["Zc"]]))
					pos_tree[j]=spatial.cKDTree(pos[j],boxsize=boxval)
			if (iverbose): print("done ",time.clock()-start)

	#now start marching backwards in time from root heads
	#identifying all subhaloes that have every been subhaloes for long enough
	#and all progenitors and group them together into the same forest id
	forestidval=1
	start=time.clock()
	for j in range(numsnaps):
		start2=time.clock()
		if (numhalos[j]==0): continue
		#now with tree start at last snapshot and identify all root heads
		#only look at halos that are their own root head and are not subhalos
		rootheads=np.where((halodata[j]['ID']==halodata[j]['RootHead'])*(halodata[j]['hostHaloID']==-1)*(halodata[j]['ForestID']==-1))
		if (iverbose): print("At snapshot",j,len(rootheads[0]))
		for iroothead in rootheads[0]:
			#if a halo has been processed as part of a forest as a
			#result of walking the subhalo branches of a different root head
			#then move on to the next object
			if (halodata[j]['ForestID'][iroothead]!=-1): continue
			AllRootHead = []
			#begin recursively searching and setting the forest using the the roothead
			AllRootHead,halodata = SetForestID(numsnaps,halodata,halodata[j]["RootHead"][iroothead],forestidval,AllRootHead)
			#update forest id
			forestidval+=1
		if (iverbose): print("Done snap",j,time.clock()-start2)

	#get the size of each forest
	ForestSize=np.zeros(forestidval,dtype=int64)
	for j in range(numsnaps):
		if (numhalos[j]==0): continue
		uniqueforest,counts=np.unique(halodata[j]['ForestID'],return_counts=True)
		for icount in range(len(uniqueforest)):
			ForestSize[uniqueforest[icount]-1]+=counts[icount]
		if (iverbose): print("Finished processing forest size for snap",j)
	start2=time.clock()

	#first identify all subhalos and see if any have subhalo connections with different than their host
	for j in range(numsnaps):
		if (numhalos[j]==0): continue
		#now with tree start at last snapshot and identify all root heads
		#only look at halos that are their own root head and are not subhalos
		missingforest=np.where((halodata[j]['ForestID']==-1))
		rootheads=np.where((halodata[j]['ID']==halodata[j]['RootHead'])*(halodata[j]['ForestID']==-1))
		subrootheads=np.where((halodata[j]['ForestID']==-1)*(halodata[j]['hostHaloID']!=-1))
		if (iverbose): print("At snapshot",j," still have ",halodata[j]['ForestID'].size,len(missingforest[0]), " with no forest id ! Of which ",len(rootheads[0])," are root heads", len(subrootheads[0]),"are subhalos")
		#if (iverbose and len(missingforest[0])>0): print("At snapshot",j," still have ",len(missingforest[0]), " with no forest id ! Of which ",len(rootheads[0])," are root heads", len(subrootheads[0]),"are subhalos")
		if (len(subrootheads[0])>0):
			for isub in subrootheads[0]:
				hostid=halodata[j]['hostHaloID'][isub]
				hostindex=int(hostid%TEMPORALHALOIDVAL-1)
				halodata[j]['ForestID'][isub]=halodata[j]['ForestID'][hostindex]
				halodata[j]['ForestLevel'][isub]=halodata[j]['ForestLevel'][hostindex]+1
	#then return this
	print("Done generating forest",time.clock()-start)
	return ForestSize

"""
Adjust halo catalog for period, comoving coords, etc
"""
def AdjustforPeriod(numsnaps,numhalos,boxsize,hval,atime,halodata,icomove=0):
	"""
	Map halo positions from 0 to box size
	"""
	for i in range(numsnaps):
		if (icomove):
			boxval=boxsize/hval
		else:
			boxval=boxsize*atime[i]/hval
		wdata=np.where(halodata[i]["Xc"]<0)
		halodata[i]["Xc"][wdata]+=boxval
		wdata=np.where(halodata[i]["Yc"]<0)
		halodata[i]["Yc"][wdata]+=boxval
		wdata=np.where(halodata[i]["Zc"]<0)
		halodata[i]["Zc"][wdata]+=boxval

		wdata=np.where(halodata[i]["Xc"]>boxval)
		halodata[i]["Xc"][wdata]-=boxval
		wdata=np.where(halodata[i]["Yc"]>boxval)
		halodata[i]["Yc"][wdata]-=boxval
		wdata=np.where(halodata[i]["Zc"]>boxval)
		halodata[i]["Zc"][wdata]-=boxval

def AdjustComove(itocomovefromphysnumsnaps,numsnaps,numhalos,atime,halodata,igas=0,istar=0):
	"""
	Convert distances to/from physical from/to comoving
	"""
	for i in range(numsnaps):
		if (numhalos[i]==0): continue
		#converting from physical to comoving
		if (itocomovefromphysnumsnaps==1):
			fac=float(1.0/atime[i])
		#converting from comoving to physical
		else:
			fac=float(atime[i])
		if (fac==1): continue

		#convert physical distances
		halodata[i]["Xc"]*=fac
		halodata[i]["Yc"]*=fac
		halodata[i]["Zc"]*=fac
		halodata[i]["Xcmbp"]*=fac
		halodata[i]["Ycmbp"]*=fac
		halodata[i]["Zcmbp"]*=fac

		#sizes
		halodata[i]["Rvir"]*=fac
		halodata[i]["R_size"]*=fac
		halodata[i]["R_200mean"]*=fac
		halodata[i]["R_200crit"]*=fac
		halodata[i]["R_BN97"]*=fac
		halodata[i]["Rmax"]*=fac
		halodata[i]["R_HalfMass"]*=fac

		#if gas
		if (igas):
			halodata[i]["Xc_gas"]*=fac
			halodata[i]["Yc_gas"]*=fac
			halodata[i]["Zc_gas"]*=fac
			halodata[i]["R_HalfMass_gas"]*=fac

		#if stars
		if (istar):
			halodata[i]["Xc_star"]*=fac
			halodata[i]["Yc_star"]*=fac
			halodata[i]["Zc_star"]*=fac
			halodata[i]["R_HalfMass_star"]*=fac

"""
Code to use individual snapshot files and merge them together into a full unified hdf file containing information determined from the tree
"""

def ProduceUnifiedTreeandHaloCatalog(fname,numsnaps,tree,numhalos,halodata,atime,
	descripdata={'Title':'Tree and Halo catalog of sim', 'VELOCIraptor_version':1.15, 'Tree_version':1.1, 'Particle_num_threshold':20, 'Temporal_linking_length':1, 'Flag_gas':False, 'Flag_star':False, 'Flag_bh':False},
	cosmodata={'Omega_m':1.0, 'Omega_b':0., 'Omega_Lambda':0., 'Hubble_param':1.0,'BoxSize':1.0, 'Sigma8':1.0},
	unitdata={'UnitLength_in_Mpc':1.0, 'UnitVelocity_in_kms':1.0,'UnitMass_in_Msol':1.0, 'Flag_physical_comoving':True,'Flag_hubble_flow':False},
	partdata={'Flag_gas':False, 'Flag_star':False, 'Flag_bh':False},
	ibuildheadtail=0, icombinefile=1):

	"""

	produces a unifed HDF5 formatted file containing the full catalog plus information to walk the tree
	\ref BuildTemporalHeadTail must have been called before otherwise it is called.
	Code produces a file for each snapshot
	The keys are the same as that contained in the halo catalog dictionary with the addition of
	Num_of_snaps, and similar header info contain in the VELOCIraptor hdf files, ie Num_of_groups, Total_num_of_groups

	\todo don't know if I should use multiprocessing here to write files in parallel. IO might not be ideal

	"""
	if (ibuildheadtail==1):
		BuildTemporalHeadTail(numsnaps,tree,numhalos,halodata)
	totnumhalos=sum(numhalos)
	if (icombinefile==1):
		hdffile=h5py.File(fname+".snap.hdf.data",'w')
		headergrp=hdffile.create_group("Header")
		#store useful information such as number of snapshots, halos,
		#cosmology (Omega_m,Omega_b,Hubble_param,Omega_Lambda, Box size)
		#units (Physical [1/0] for physical/comoving flag, length in Mpc, km/s, solar masses, Gravity
		#and TEMPORALHALOIDVAL used to traverse tree information (converting halo ids to haloindex or snapshot), Reverse_order [1/0] for last snap listed first)
		#set the attributes of the header
		headergrp.attrs["NSnaps"]=numsnaps
		#overall description
		#simulation box size

		#cosmological params
		cosmogrp=headergrp.create_group("Cosmology")
		for key in cosmodata.keys():
			cosmogrp.attrs[key]=cosmodata[key]
		#unit params
		unitgrp=headergrp.create_group("Units")
		for key in unitdata.keys():
			unitgrp.attrs[key]=unitdata[key]
		#particle types
		partgrp=headergrp.create_group("Parttypes")
		partgrp.attrs["Flag_gas"]=descripdata["Flag_gas"]
		partgrp.attrs["Flag_star"]=descripdata["Flag_star"]
		partgrp.attrs["Flag_bh"]=descripdata["Flag_bh"]

		for i in range(numsnaps):
			snapgrp=hdffile.create_group("Snap_%03d"%(numsnaps-1-i))
			snapgrp.attrs["Snapnum"]=(numsnaps-1-i)
			snapgrp.attrs["NHalos"]=numhalos[i]
			snapgrp.attrs["scalefactor"]=atime[i]
			for key in halodata[i].keys():
				snapgrp.create_dataset(key,data=halodata[i][key])
		hdffile.close()
	else:
		for i in range(numsnaps):
			hdffile=h5py.File(fname+".snap_%03d.hdf.data"%(numsnaps-1-i),'w')
			hdffile.create_dataset("Snap_value",data=np.array([numsnaps-1-i],dtype=np.uint32))
			hdffile.create_dataset("NSnaps",data=np.array([numsnaps],dtype=np.uint32))
			hdffile.create_dataset("NHalos",data=np.array([numhalos[i]],dtype=np.uint64))
			hdffile.create_dataset("TotalNHalos",data=np.array([totnumhalos],dtype=np.uint64))
			hdffile.create_dataset("scalefactor",data=np.array([atime[i]],dtype=np.float64))
			for key in halodata[i].keys():
				hdffile.create_dataset(key,data=halodata[i][key])
			hdffile.close()

	hdffile=h5py.File(fname+".tree.hdf.data",'w')
	hdffile.create_dataset("NSnaps",data=np.array([numsnaps],dtype=np.uint32))
	hdffile.create_dataset("TotalNHalos",data=np.array([totnumhalos],dtype=np.uint64))
	hdffile.create_dataset("NHalos",data=np.array([numhalos],dtype=np.uint64))
	for i in range(numsnaps):
		snapgrp=hdffile.create_group("Snap_%03d"%(numsnaps-1-i))
		for key in tree[i].keys():
			"""
			#to be completed for progenitor list
			if (key=="Progen"):
				for j in range(numhalos[i]):
					halogrp=snapgrp.create_group("Halo"+str(j))
					halogrp.create_dataset(key,data=tree[i][key][j])
			else:
				snapgrp.create_dataset(key,data=tree[i][key])
			"""
			if ((key=="Progen") | (key=="Descen")): continue
			snapgrp.create_dataset(key,data=tree[i][key])
	hdffile.close()

def ProduceCombinedUnifiedTreeandHaloCatalog(fname,numsnaps,tree,numhalos,halodata,atime,
	descripdata={'Title':'Tree and Halo catalog of sim', 'VELOCIraptor_version':1.15, 'Tree_version':1.1, 'Particle_num_threshold':20, 'Temporal_linking_length':1, 'Flag_gas':False, 'Flag_star':False, 'Flag_bh':False},
	cosmodata={'Omega_m':1.0, 'Omega_b':0., 'Omega_Lambda':0., 'Hubble_param':1.0,'BoxSize':1.0, 'Sigma8':1.0},
	unitdata={'UnitLength_in_Mpc':1.0, 'UnitVelocity_in_kms':1.0,'UnitMass_in_Msol':1.0, 'Flag_physical_comoving':True,'Flag_hubble_flow':False},
	partdata={'Flag_gas':False, 'Flag_star':False, 'Flag_bh':False},
	ibuildheadtail=0,ibuildmajormergers=0, TEMPORALHALOIDVAL=1000000000000):

	"""
	produces a unifed HDF5 formatted file containing the full catalog plus information to walk the tree
	#ref BuildTemporalHeadTail must have been called before otherwise it is called.
	Code produces a file for each snapshot
	The keys are the same as that contained in the halo catalog dictionary with the addition of
	Num_of_snaps, and similar header info contain in the VELOCIraptor hdf files, ie Num_of_groups, Total_num_of_groups

	#todo don't know if I should use multiprocessing here to write files in parallel. IO might not be ideal

	Here the halodata is the dictionary contains the information

	"""

	if (ibuildheadtail==1):
		BuildTemporalHeadTail(numsnaps,tree,numhalos,halodata)
	if (ibuildmajormergers==1):
		IdentifyMergers(numsnaps,tree,numhalos,halodata,boxsize,hval,atime)
	hdffile=h5py.File(fname+".snap.hdf.data",'w')
	headergrp=hdffile.create_group("Header")
	#store useful information such as number of snapshots, halos,
	#cosmology (Omega_m,Omega_b,Hubble_param,Omega_Lambda, Box size)
	#units (Physical [1/0] for physical/comoving flag, length in Mpc, km/s, solar masses, Gravity
	#and TEMPORALHALOIDVAL used to traverse tree information (converting halo ids to haloindex or snapshot), Reverse_order [1/0] for last snap listed first)
	#set the attributes of the header
	headergrp.attrs["NSnaps"]=numsnaps
	#overall description
	headergrp.attrs["Title"]=descripdata["Title"]
	#simulation box size
	headergrp.attrs["BoxSize"]=cosmodata["BoxSize"]
	findergrp=headergrp.create_group("HaloFinder")
	findergrp.attrs["Name"]="VELOCIraptor"
	findergrp.attrs["Version"]=descripdata["VELOCIraptor_version"]
	findergrp.attrs["Particle_num_threshold"]=descripdata["Particle_num_threshold"]

	treebuildergrp=headergrp.create_group("TreeBuilder")
	treebuildergrp.attrs["Name"]="VELOCIraptor-Tree"
	treebuildergrp.attrs["Version"]=descripdata["Tree_version"]
	treebuildergrp.attrs["Temporal_linking_length"]=descripdata["Temporal_linking_length"]

	#cosmological params
	cosmogrp=headergrp.create_group("Cosmology")
	for key in cosmodata.keys():
		if (key!='BoxSize'): cosmogrp.attrs[key]=cosmodata[key]
	#unit params
	unitgrp=headergrp.create_group("Units")
	for key in unitdata.keys():
		unitgrp.attrs[key]=unitdata[key]
	#particle types
	partgrp=headergrp.create_group("Parttypes")
	partgrp.attrs["Flag_gas"]=descripdata["Flag_gas"]
	partgrp.attrs["Flag_star"]=descripdata["Flag_star"]
	partgrp.attrs["Flag_bh"]=descripdata["Flag_bh"]

	#now have finished with header

	#now need to create groups for halos and then a group containing tree information
	snapsgrp=hdffile.create_group("Snapshots")
	#internal tree keys
	treekeys=["RootHead", "RootHeadSnap", "Head", "HeadSnap", "Tail", "TailSnap", "RootTail", "RootTailSnap", "Num_progen"]

	for i in range(numsnaps):
		#note that I normally have information in reverse order so that might be something in the units
		snapgrp=snapsgrp.create_group("Snap_%03d"%(numsnaps-1-i))
		snapgrp.attrs["Snapnum"]=i
		snapgrp.attrs["NHalos"]=numhalos[i]
		snapgrp.attrs["scalefactor"]=atime[i]
	#now close file and use the pytables interface so as to write the table
	hdffile.close()
	#now write tables using pandas interface
	for i in range(numsnaps):
		#lets see if we can alter the code to write a table
		keys=halodata[i].keys()
		#remove tree keys
		for tkey in treekeys: keys.remove(tkey)
		#make temp dict
		dictval=dict()
		for key in keys:
			dictval[key]=halodata[i][key]
		#make a pandas DataFrame using halo dictionary
		df=pd.DataFrame.from_dict(dictval)
		df.to_hdf(fname+".snap.hdf.data","Snapshots/Snap_%03d/Halos"%(numsnaps-1-i), format='table', mode='a')

	#reopen with h5py interface
	hdffile=h5py.File(fname+".snap.hdf.data",'a')
	#then write tree information in separate group
	treegrp=hdffile.create_group("MergerTree")
	#Tree group should contain
	"""
		HaloSnapID
		HaloSnapNum
		HaloSnapIndex
		ProgenitorIndex
		ProgenitorSnapnum
		ProgenitorID
		DescendantIndex
		..
		..
		RootProgenitorIndex
		..
		..
		RootDescendantIndex
	"""
	#to save on memory, allocate each block separately
	#store halo information
	tothalos=sum(numhalos)
	tdata=np.zeros(tothalos,dtype=halodata[0]["ID"].dtype)
	count=0
	for i in range(numsnaps):
		tdata[count:int(numhalos[i])+count]=halodata[i]["ID"]
		count+=int(numhalos[i])
	treegrp.create_dataset("HaloSnapID",data=tdata)
	tdata=np.zeros(tothalos,dtype=np.uint32)
	count=0
	for i in range(numsnaps):
		tdata[count:int(numhalos[i])+count]=i
		count+=int(numhalos[i])
	treegrp.create_dataset("HaloSnapNum",data=tdata)
	tdata=np.zeros(tothalos,dtype=np.uint64)
	count=0
	for i in range(numsnaps):
		tdata[count:int(numhalos[i])+count]=range(int(numhalos[i]))
		count+=int(numhalos[i])
	treegrp.create_dataset("HaloSnapIndex",data=tdata)
	#store progenitors
	tdata=np.zeros(tothalos,dtype=halodata[0]["Tail"].dtype)
	count=0
	for i in range(numsnaps):
		tdata[count:int(numhalos[i])+count]=halodata[i]["Tail"]
		count+=int(numhalos[i])
	treegrp.create_dataset("ProgenitorID",data=tdata)
	tdata=np.zeros(tothalos,dtype=halodata[0]["TailSnap"].dtype)
	count=0
	for i in range(numsnaps):
		tdata[count:int(numhalos[i])+count]=halodata[i]["TailSnap"]
		count+=int(numhalos[i])
	treegrp.create_dataset("ProgenitorSnapnum",data=tdata)
	tdata=np.zeros(tothalos,dtype=np.uint64)
	count=0
	for i in range(numsnaps):
		tdata[count:int(numhalos[i])+count]=(halodata[i]["Tail"]%TEMPORALHALOIDVAL-1)
		count+=int(numhalos[i])
	treegrp.create_dataset("ProgenitorIndex",data=tdata)
	#store descendants
	tdata=np.zeros(tothalos,dtype=halodata[0]["Head"].dtype)
	count=0
	for i in range(numsnaps):
		tdata[count:int(numhalos[i])+count]=halodata[i]["Head"]
		count+=int(numhalos[i])
	treegrp.create_dataset("DescendantID",data=tdata)
	tdata=np.zeros(tothalos,dtype=halodata[0]["HeadSnap"].dtype)
	count=0
	for i in range(numsnaps):
		tdata[count:int(numhalos[i])+count]=halodata[i]["HeadSnap"]
		count+=int(numhalos[i])
	treegrp.create_dataset("DescendantSnapnum",data=tdata)
	tdata=np.zeros(tothalos,dtype=np.uint64)
	count=0
	for i in range(numsnaps):
		tdata[count:int(numhalos[i])+count]=(halodata[i]["Head"]%TEMPORALHALOIDVAL-1)
		count+=int(numhalos[i])
	treegrp.create_dataset("DescendantIndex",data=tdata)
	#store progenitors
	tdata=np.zeros(tothalos,dtype=halodata[0]["RootTail"].dtype)
	count=0
	for i in range(numsnaps):
		tdata[count:int(numhalos[i])+count]=halodata[i]["RootTail"]
		count+=int(numhalos[i])
	treegrp.create_dataset("RootProgenitorID",data=tdata)
	tdata=np.zeros(tothalos,dtype=halodata[0]["RootTailSnap"].dtype)
	count=0
	for i in range(numsnaps):
		tdata[count:int(numhalos[i])+count]=halodata[i]["RootTailSnap"]
		count+=int(numhalos[i])
	treegrp.create_dataset("RootProgenitorSnapnum",data=tdata)
	tdata=np.zeros(tothalos,dtype=np.uint64)
	count=0
	for i in range(numsnaps):
		tdata[count:int(numhalos[i])+count]=(halodata[i]["RootTail"]%TEMPORALHALOIDVAL-1)
		count+=int(numhalos[i])
	treegrp.create_dataset("RootProgenitorIndex",data=tdata)
	#store descendants
	tdata=np.zeros(tothalos,dtype=halodata[0]["RootHead"].dtype)
	count=0
	for i in range(numsnaps):
		tdata[count:int(numhalos[i])+count]=halodata[i]["RootHead"]
		count+=int(numhalos[i])
	treegrp.create_dataset("RootDescendantID",data=tdata)
	tdata=np.zeros(tothalos,dtype=halodata[0]["RootHeadSnap"].dtype)
	count=0
	for i in range(numsnaps):
		tdata[count:int(numhalos[i])+count]=halodata[i]["RootHeadSnap"]
		count+=int(numhalos[i])
	treegrp.create_dataset("RootDescendantSnapnum",data=tdata)
	tdata=np.zeros(tothalos,dtype=np.uint64)
	count=0
	for i in range(numsnaps):
		tdata[count:int(numhalos[i])+count]=(halodata[i]["RootHead"]%TEMPORALHALOIDVAL-1)
		count+=int(numhalos[i])
	treegrp.create_dataset("RootDescendantIndex",data=tdata)
	#store number of progenitors
	tdata=np.zeros(tothalos,dtype=np.uint32)
	count=0
	for i in range(numsnaps):
		tdata[count:int(numhalos[i])+count]=halodata[i]["Num_progen"]
		count+=int(numhalos[i])
	treegrp.create_dataset("NProgen",data=tdata)

	hdffile.close()

def ReadUnifiedTreeandHaloCatalog(fname, desiredfields=[], icombinedfile=1,iverbose=1):
	"""
	Read Unified Tree and halo catalog from HDF file with base filename fname.

	Parameters
	----------

	Returns
	-------
	"""
	if (icombinedfile):
		hdffile=h5py.File(fname+".snap.hdf.data",'r')

		#load data sets containing number of snaps
		headergrpname="Header/"
		numsnaps=hdffile[headergrpname].attrs["NSnaps"]

		#allocate memory
		halodata=[dict() for i in range(numsnaps)]
		numhalos=[0 for i in range(numsnaps)]
		atime=[0 for i in range(numsnaps)]
		tree=[[] for i in range(numsnaps)]
		cosmodata=dict()
		unitdata=dict()

		#load cosmology data
		cosmogrpname="Cosmology/"
		fieldnames=[str(n) for n in hdffile[headergrpname+cosmogrpname].attrs.keys()]
		for fieldname in fieldnames:
			cosmodata[fieldname]=hdffile[headergrpname+cosmogrpname].attrs[fieldname]

		#load unit data
		unitgrpname="Units/"
		fieldnames=[str(n) for n in hdffile[headergrpname+unitgrpname].attrs.keys()]
		for fieldname in fieldnames:
			unitdata[fieldname]=hdffile[headergrpname+unitgrpname].attrs[fieldname]

		#for each snap load the appropriate group
		start=time.clock()
		for i in range(numsnaps):
			snapgrpname="Snap_%03d/"%(numsnaps-1-i)
			if (iverbose==1):
				print("Reading ",snapgrpname)
			isnap=hdffile[snapgrpname].attrs["Snapnum"]
			atime[isnap]=hdffile[snapgrpname].attrs["scalefactor"]
			numhalos[isnap]=hdffile[snapgrpname].attrs["NHalos"]
			if (len(desiredfields)>0):
				fieldnames=desiredfields
			else:
				fieldnames=[str(n) for n in hdffile[snapgrpname].keys()]
			for catvalue in fieldnames:
				halodata[isnap][catvalue]=np.array(hdffile[snapgrpname+catvalue])
		hdffile.close()
		print("read halo data ",time.clock()-start)
	else :
		hdffile=h5py.File(fname+".snap_000.hdf.data",'r')
		numsnaps=int(hdffile["NSnaps"][0])
		#get field names
		fieldnames=[str(n) for n in hdffile.keys()]
		#clean of header info
		fieldnames.remove("Snapnum")
		fieldnames.remove("NSnaps")
		fieldnames.remove("NHalos")
		fieldnames.remove("TotalNHalos")
		fieldnames.remove("scalefactor")
		if (len(desiredfields)>0):
			fieldnames=desiredfields
		hdffile.close()
		halodata=[[] for i in range(numsnaps)]
		numhalos=[0 for i in range(numsnaps)]
		atime=[0 for i in range(numsnaps)]
		tree=[[] for i in range(numsnaps)]
		start=time.clock()
		for i in range(numsnaps):
			hdffile=h5py.File(fname+".snap_%03d.hdf.data"%(numsnaps-1-i),'r')
			atime[i]=(hdffile["scalefactor"])[0]
			numhalos[i]=(hdffile["NHalos"])[0]
			halodata[i]=dict()
			for catvalue in fieldnames:
				halodata[i][catvalue]=np.array(hdffile[catvalue])
			hdffile.close()
		print("read halo data ",time.clock()-start)
	#lets ignore the tree file for now
	for i in range(numsnaps):
		tree[i]=dict()
	return atime,tree,numhalos,halodata,cosmodata,unitdata
	if (icombinedfile==1):
		hdffile=h5py.File(fname+".tree.hdf.data",'r')
		treefields=["haloID", "Num_progen"]
		#do be completed for Progenitor list although information is contained in the halo catalog by searching for things with the same head
		#treefields=["haloID", "Num_progen", "Progen"]
		for i in range(numsnaps):
			snapgrpname="Snap_%03d/"%(numsnaps-1-i)
			tree[i]=dict()
			for catvalue in treefields:
				"""
				if (catvalue==treefields[-1]):
					tree[i][catvalue]=[[]for j in range(numhalos[i])]
					for j in range(numhalos[i]):
						halogrpname=snapgrpname+"/Halo"+str(j)
						tree[i][catvalue]=np.array(hdffile[halogrpname+catvalue])
				else:
					tree[i][catvalue]=np.array(hdffile[snapgrpname+catvalue])
				"""
				tree[i][catvalue]=np.array(hdffile[snapgrpname+catvalue])
		hdffile.close()
	return atime,tree,numhalos,halodata,cosmodata,unitdata

def ProduceHDFTree(fname,numsnaps,tree,numhalos,halodata,atime,
	descripdata={'Title':'Tree catalogue', 'VELOCIraptor_version':1.3, 'Tree_version':1.1, 'Particle_num_threshold':20, 'Temporal_linking_length':1, 'Flag_gas':False, 'Flag_star':False, 'Flag_bh':False}
	):
	"""
	Produces a HDF5 formatted file containing Reduced Tree information,
	ie; RootHead, Head, HeadSnap, Tail, RootTail, etc.

	Parameters
	----------
	fname : string
		filename of the hdf file to be written
	numsnaps : int
		the number of snapshots
	tree : dict
		the tree data
	numhalos : array
		array of number of halos per snapshot
	halodata : dict
		the halo data dictionary
	atime : array
		array of scalefactors/times of the snaphots
	discrptdata : dict
		stores a description of how the tree catalogue was produced

	Returns
	-------
	void :
		Only writes an hdf file. Nothing is returned.
	"""
	hdffile=h5py.File(fname,'w')
	headergrp=hdffile.create_group("Header")
	#store useful information such as number of snapshots, halos,
	#cosmology (Omega_m,Omega_b,Hubble_param,Omega_Lambda, Box size)
	#units (Physical [1/0] for physical/comoving flag, length in Mpc, km/s, solar masses, Gravity
	#and TEMPORALHALOIDVAL used to traverse tree information (converting halo ids to haloindex or snapshot), Reverse_order [1/0] for last snap listed first)
	#set the attributes of the header
	headergrp.attrs["NSnaps"]=numsnaps
	#overall description
	headergrp.attrs["Title"]=descripdata["Title"]
	findergrp=headergrp.create_group("HaloFinder")
	findergrp.attrs["Name"]="VELOCIraptor"
	findergrp.attrs["Version"]=descripdata["VELOCIraptor_version"]
	findergrp.attrs["Particle_num_threshold"]=descripdata["Particle_num_threshold"]

	treebuildergrp=headergrp.create_group("TreeBuilder")
	treebuildergrp.attrs["Name"]="VELOCIraptor-Tree"
	treebuildergrp.attrs["Version"]=descripdata["Tree_version"]
	treebuildergrp.attrs["Temporal_linking_length"]=descripdata["Temporal_linking_length"]

	#now need to create groups for halos and then a group containing tree information
	snapsgrp=hdffile.create_group("Snapshots")
	#internal tree keys
	halokeys=["RootHead", "RootHeadSnap", "Head", "HeadSnap", "Tail", "TailSnap", "RootTail", "RootTailSnap", "ID", "Num_progen"]

	for i in range(numsnaps):
		#note that I normally have information in reverse order so that might be something in the units
		snapgrp=snapsgrp.create_group("Snap_%03d"%i)
		snapgrp.attrs["Snapnum"]=i
		snapgrp.attrs["NHalos"]=numhalos[i]
		snapgrp.attrs["scalefactor"]=atime[i]
		for key in halokeys:
			snapgrp.create_dataset(key,data=halodata[i][key])
	hdffile.close()


"""
	Conversion Tools
"""


def ConvertASCIIPropertyFileToHDF(basefilename,iseparatesubfiles=0,iverbose=0):
	"""
	Reads an ASCII file and converts it to the HDF format for VELOCIraptor properties files

	"""
	inompi=True
	if (iverbose): print("reading properties file and converting to hdf",basefilename,os.path.isfile(basefilename))
	filename=basefilename+".properties"
	#load header
	if (os.path.isfile(basefilename)==True):
		numfiles=0
	else:
		filename=basefilename+".properties"+".0"
		inompi=False
		if (os.path.isfile(filename)==False):
			print("file not found")
			return []
	byteoffset=0
	#load ascii file
	halofile = open(filename, 'r')
	#read header information
	[filenum,numfiles]=halofile.readline().split()
	filenum=int(filenum);numfiles=int(numfiles)
	[numhalos, numtothalos]= halofile.readline().split()
	numhalos=np.uint64(numhalos);numtothalos=np.uint64(numtothalos)
	names = ((halofile.readline())).split()
	#remove the brackets in ascii file names
	fieldnames= [fieldname.split("(")[0] for fieldname in names]
	halofile.close()

	for ifile in range(numfiles):
		if (inompi==True):
			filename=basefilename+".properties"
			hdffilename=basefilename+".hdf.properties"
		else:
			filename=basefilename+".properties"+"."+str(ifile)
			hdffilename=basefilename+".hdf.properties"+"."+str(ifile)
		if (iverbose) : print("reading ",filename)
		halofile = open(filename, 'r')
		hdffile=h5py.File(hdffilename,'w')
		[filenum,numfiles]=halofile.readline().split()
		[numhalos, numtothalos]= halofile.readline().split()
		filenum=int(filenum);numfiles=int(numfiles)
		numhalos=np.uint64(numhalos);numtothalos=np.uint64(numtothalos)
		#write header info
		hdffile.create_dataset("File_id",data=np.array([filenum]))
		hdffile.create_dataset("Num_of_files",data=np.array([numfiles]))
		hdffile.create_dataset("Num_of_groups",data=np.array([numhalos]));
		hdffile.create_dataset("Total_num_of_groups",data=np.array([numtothalos]));
		halofile.close()
		if (numhalos>0): htemp = np.loadtxt(filename,skiprows=3).transpose()
		else: htemp=[[]for ikeys in range(len(fieldnames))]
		for ikeys in range(len(fieldnames)):
			if (fieldnames[ikeys]=="ID"):
				hdffile.create_dataset(fieldnames[ikeys],data=np.array(htemp[ikeys],dtype=np.uint64))
			elif (fieldnames[ikeys]=="ID_mbp"):
				hdffile.create_dataset(fieldnames[ikeys],data=np.array(htemp[ikeys],dtype=np.int64))
			elif (fieldnames[ikeys]=="hostHaloID"):
				hdffile.create_dataset(fieldnames[ikeys],data=np.array(htemp[ikeys],dtype=np.int64))
			elif fieldnames[ikeys] in ["numSubStruct","npart","n_gas","n_star"]:
				hdffile.create_dataset(fieldnames[ikeys],data=np.array(htemp[ikeys], dtype=np.uint64))
			else:
				hdffile.create_dataset(fieldnames[ikeys],data=np.array(htemp[ikeys], dtype=np.float64))

		hdffile.close()
	#if subhalos are written in separate files, then read them too
	if (iseparatesubfiles==1):
		for ifile in range(numfiles):
			if (inompi==True):
				filename=basefilename+".sublevels"+".properties"
				hdffilename=basefilename+".hdf"+".sublevels"+".properties"
			else:
				filename=basefilename+".sublevels"+".properties"+"."+str(ifile)
				hdffilename=basefilename+".hdf"+".sublevels"+".properties"+"."+str(ifile)
			if (iverbose) : print("reading ",filename)
			halofile = open(filename, 'r')
			hdffile=h5py.File(hdffilename,'w')
			[filenum,numfiles]=halofile.readline().split()
			[numhalos, numtothalos]= halofile.readline().split()
			filenum=int(filenum);numfiles=int(numfiles)
			numhalos=np.uint64(numhalos);numtothalos=np.uint64(numtothalos)
			#write header info
			hdffile.create_dataset("File_id",data=np.array([filenum]))
			hdffile.create_dataset("Num_of_files",data=np.array([numfiles]))
			hdffile.create_dataset("Num_of_groups",data=np.array([numhalos]));
			hdffile.create_dataset("Total_num_of_groups",data=np.array([numtothalos]));
			halofile.close()
			if (numhalos>0): htemp = np.loadtxt(filename,skiprows=3).transpose()
			else: htemp=[[]for ikeys in range(len(fieldnames))]
			for ikeys in range(len(fieldnames)):
				if (fieldnames[ikeys]=="ID"):
					hdffile.create_dataset(fieldnames[ikeys],data=np.array(htemp[ikeys],dtype=np.uint64))
				elif (fieldnames[ikeys]=="ID_mbp"):
					hdffile.create_dataset(fieldnames[ikeys],data=np.array(htemp[ikeys],dtype=np.int64))
				elif (fieldnames[ikeys]=="hostHaloID"):
					hdffile.create_dataset(fieldnames[ikeys],data=np.array(htemp[ikeys],dtype=np.int64))
				elif fieldnames[ikeys] in ["numSubStruct","npart","n_gas","n_star"]:
					hdffile.create_dataset(fieldnames[ikeys],data=np.array(htemp[ikeys], dtype=np.uint64))
				else:
					hdffile.create_dataset(fieldnames[ikeys],data=np.array(htemp[ikeys], dtype=np.float64))
			hdffile.close()

def ConvertASCIICatalogGroupsFileToHDF(basefilename,iseparatesubfiles=0,iverbose=0):
	"""
	Reads an ASCII file and converts it to the HDF format for VELOCIraptor files

	"""
	inompi=True
	if (iverbose): print("reading properties file and converting to hdf",basefilename,os.path.isfile(basefilename))
	filename=basefilename+".catalog_groups"
	#load header
	if (os.path.isfile(basefilename)==True):
		numfiles=0
	else:
		filename=basefilename+".catalog_groups"+".0"
		inompi=False
		if (os.path.isfile(filename)==False):
			print("file not found")
			return []
	byteoffset=0
	#load ascii file
	halofile = open(filename, 'r')
	#read header information
	[filenum,numfiles]=halofile.readline().split()
	filenum=int(filenum);numfiles=int(numfiles)
	[numhalos, numtothalos]= halofile.readline().split()
	numhalos=np.uint64(numhalos);numtothalos=np.uint64(numtothalos)
	halofile.close()

	fieldnames=["Group_Size","Offset","Offset_unbound","Number_of_substructures_in_halo","Parent_halo_ID"]
	fieldtype=[np.uint32,np.uint64,np.uint64,np.uint32,np.int64]

	for ifile in range(numfiles):
		if (inompi==True):
			filename=basefilename+".catalog_groups"
			hdffilename=basefilename+".hdf.catalog_groups"
		else:
			filename=basefilename+".catalog_groups"+"."+str(ifile)
			hdffilename=basefilename+".hdf.catalog_groups"+"."+str(ifile)
		if (iverbose) : print("reading ",filename)
		halofile = open(filename, 'r')
		hdffile=h5py.File(hdffilename,'w')
		[filenum,numfiles]=halofile.readline().split()
		[numhalos, numtothalos]= halofile.readline().split()
		filenum=int(filenum);numfiles=int(numfiles)
		numhalos=np.uint64(numhalos);numtothalos=np.uint64(numtothalos)
		#write header info
		hdffile.create_dataset("File_id",data=np.array([filenum]))
		hdffile.create_dataset("Num_of_files",data=np.array([numfiles]))
		hdffile.create_dataset("Num_of_groups",data=np.array([numhalos]));
		hdffile.create_dataset("Total_num_of_groups",data=np.array([numtothalos]));
		halofile.close()
		if (numhalos>0):
			#will look like one dimensional array of values split into
			#"Group_Size"
			#"Offset"
			#"Offset_unbound"
			#"Number_of_substructures_in_halo"
			#"Parent_halo_ID"
			#each of size numhalos
			cattemp = np.loadtxt(filename,skiprows=2).transpose()
			for ikeys in range(len(fieldnames)):
				hdffile.create_dataset(fieldnames[ikeys],data=np.array(cattemp[ikeys*numhalos:(ikeys+1)*numhalos],dtype=fieldtype[ikeys]))
		else:
			cattemp=[]
			for ikeys in range(len(fieldnames)):
				hdffile.create_dataset(fieldnames[ikeys],data=np.array([],dtype=fieldtype[ikeys]))
		hdffile.close()
	#if subhalos are written in separate files, then read them too
	if (iseparatesubfiles==1):
		for ifile in range(numfiles):
			if (inompi==True):
				filename=basefilename+".sublevels"+".catalog_groups"
				hdffilename=basefilename+".hdf"+".sublevels"+".catalog_groups"
			else:
				filename=basefilename+".sublevels"+".catalog_groups"+"."+str(ifile)
				hdffilename=basefilename+".hdf"+".sublevels"+".catalog_groups"+"."+str(ifile)
			if (iverbose) : print("reading ",filename)
			halofile = open(filename, 'r')
			hdffile=h5py.File(hdffilename,'w')
			[filenum,numfiles]=halofile.readline().split()
			[numhalos, numtothalos]= halofile.readline().split()
			filenum=int(filenum);numfiles=int(numfiles)
			numhalos=np.uint64(numhalos);numtothalos=np.uint64(numtothalos)
			#write header info
			hdffile.create_dataset("File_id",data=np.array([filenum]))
			hdffile.create_dataset("Num_of_files",data=np.array([numfiles]))
			hdffile.create_dataset("Num_of_groups",data=np.array([numhalos]));
			hdffile.create_dataset("Total_num_of_groups",data=np.array([numtothalos]));
			halofile.close()
			if (numhalos>0):
				cattemp = np.loadtxt(filename,skiprows=2).transpose()
				for ikeys in range(len(fieldnames)):
					hdffile.create_dataset(fieldnames[ikeys],data=np.array(cattemp[ikeys*numhalos:(ikeys+1)*numhalos],dtype=fieldtype[ikeys]))
			else:
				cattemp=[]
				for ikeys in range(len(fieldnames)):
					hdffile.create_dataset(fieldnames[ikeys],data=np.array([],dtype=fieldtype[ikeys]))
			hdffile.close()

def ConvertASCIICatalogParticleFileToHDF(basefilename,iunbound=0,iseparatesubfiles=0,iverbose=0):
	"""
	Reads an ASCII file and converts it to the HDF format for VELOCIraptor files
	"""
	inompi=True
	if (iverbose): print("reading properties file and converting to hdf",basefilename,os.path.isfile(basefilename))
	filename=basefilename+".catalog_particles"
	if (iunbound>0): filename+=".unbound"
	#load header
	if (os.path.isfile(basefilename)==True):
		numfiles=0
	else:
		filename=basefilename+".catalog_particles"
		if (iunbound>0): filename+=".unbound"
		filename+=".0"
		inompi=False
		if (os.path.isfile(filename)==False):
			print("file not found")
			return []
	byteoffset=0
	#load ascii file
	halofile = open(filename, 'r')
	#read header information
	[filenum,numfiles]=halofile.readline().split()
	filenum=int(filenum);numfiles=int(numfiles)
	[numhalos, numtothalos]= halofile.readline().split()
	numhalos=np.uint64(numhalos);numtothalos=np.uint64(numtothalos)
	halofile.close()

	for ifile in range(numfiles):
		if (inompi==True):
			filename=basefilename+".catalog_particles"
			hdffilename=basefilename+".hdf.catalog_particles"
			if (iunbound>0):
				filename+=".unbound"
				hdffilename+=".unbound"
		else:
			filename=basefilename+".catalog_particles"
			hdffilename=basefilename+".hdf.catalog_particles"
			if (iunbound>0):
				filename+=".unbound"
				hdffilename+=".unbound"
			filename+="."+str(ifile)
			hdffilename+="."+str(ifile)
		if (iverbose) : print("reading ",filename)
		halofile = open(filename, 'r')
		hdffile=h5py.File(hdffilename,'w')
		[filenum,numfiles]=halofile.readline().split()
		[numhalos, numtothalos]= halofile.readline().split()
		filenum=int(filenum);numfiles=int(numfiles)
		numhalos=np.uint64(numhalos);numtothalos=np.uint64(numtothalos)
		#write header info
		hdffile.create_dataset("File_id",data=np.array([filenum]))
		hdffile.create_dataset("Num_of_files",data=np.array([numfiles]))
		hdffile.create_dataset("Num_of_particles_in_groups",data=np.array([numhalos]));
		hdffile.create_dataset("Total_num_of_particles_in_all_groups",data=np.array([numtothalos]));
		halofile.close()
		if (numhalos>0): cattemp = np.loadtxt(filename,skiprows=2).transpose()
		else: cattemp=[]
		hdffile.create_dataset("Particle_IDs",data=np.array(cattemp,dtype=np.int64))
		hdffile.close()
	#if subhalos are written in separate files, then read them too
	if (iseparatesubfiles==1):
		for ifile in range(numfiles):
			if (inompi==True):
				filename=basefilename+".sublevels"+".catalog_particles"
				hdffilename=basefilename+".hdf"+".sublevels"+".catalog_particles"
			else:
				filename=basefilename+".sublevels"+".catalog_particles"+"."+str(ifile)
				hdffilename=basefilename+".hdf"+".sublevels"+".catalog_particles"+"."+str(ifile)
			if (iverbose) : print("reading ",filename)
			halofile = open(filename, 'r')
			hdffile=h5py.File(hdffilename,'w')
			[filenum,numfiles]=halofile.readline().split()
			[numhalos, numtothalos]= halofile.readline().split()
			filenum=int(filenum);numfiles=int(numfiles)
			numhalos=np.uint64(numhalos);numtothalos=np.uint64(numtothalos)
			#write header info
			hdffile.create_dataset("File_id",data=np.array([filenum]))
			hdffile.create_dataset("Num_of_files",data=np.array([numfiles]))
			hdffile.create_dataset("Num_of_particles_in_groups",data=np.array([numhalos]));
			hdffile.create_dataset("Total_num_of_particles_in_all_groups",data=np.array([numtothalos]));
			halofile.close()
			if (numhalos>0): cattemp = np.loadtxt(filename,skiprows=2).transpose()
			else: cattemp=[]
			hdffile.create_dataset("Particle_IDs",data=np.array(cattemp,dtype=np.int64))
			hdffile.close()

def ConvertASCIICatalogParticleTypeFileToHDF(basefilename,iunbound=0,iseparatesubfiles=0,iverbose=0):
	"""
	Reads an ASCII file and converts it to the HDF format for VELOCIraptor files
	"""
	inompi=True
	if (iverbose): print("reading properties file and converting to hdf",basefilename,os.path.isfile(basefilename))
	filename=basefilename+".catalog_parttypes"
	if (iunbound>0): filename+=".unbound"
	#load header
	if (os.path.isfile(basefilename)==True):
		numfiles=0
	else:
		filename=basefilename+".catalog_parttypes"
		if (iunbound>0): filename+=".unbound"
		filename+=".0"
		inompi=False
		if (os.path.isfile(filename)==False):
			print("file not found")
			return []
	byteoffset=0
	#load ascii file
	halofile = open(filename, 'r')
	#read header information
	[filenum,numfiles]=halofile.readline().split()
	filenum=int(filenum);numfiles=int(numfiles)
	[numhalos, numtothalos]= halofile.readline().split()
	numhalos=np.uint64(numhalos);numtothalos=np.uint64(numtothalos)
	halofile.close()

	for ifile in range(numfiles):
		if (inompi==True):
			filename=basefilename+".catalog_parttypes"
			hdffilename=basefilename+".hdf.catalog_parttypes"
			if (iunbound>0):
				filename+=".unbound"
				hdffilename+=".unbound"
		else:
			filename=basefilename+".catalog_parttypes"
			hdffilename=basefilename+".hdf.catalog_parttypes"
			if (iunbound>0):
				filename+=".unbound"
				hdffilename+=".unbound"
			filename+="."+str(ifile)
			hdffilename+="."+str(ifile)
		if (iverbose) : print("reading ",filename)
		halofile = open(filename, 'r')
		hdffile=h5py.File(hdffilename,'w')
		[filenum,numfiles]=halofile.readline().split()
		[numhalos, numtothalos]= halofile.readline().split()
		filenum=int(filenum);numfiles=int(numfiles)
		numhalos=np.uint64(numhalos);numtothalos=np.uint64(numtothalos)
		#write header info
		hdffile.create_dataset("File_id",data=np.array([filenum]))
		hdffile.create_dataset("Num_of_files",data=np.array([numfiles]))
		hdffile.create_dataset("Num_of_particles_in_groups",data=np.array([numhalos]));
		hdffile.create_dataset("Total_num_of_particles_in_all_groups",data=np.array([numtothalos]));
		halofile.close()
		if (numhalos>0): cattemp = np.loadtxt(filename,skiprows=2).transpose()
		else: cattemp=[]
		hdffile.create_dataset("Particle_Types",data=np.array(cattemp,dtype=np.int64))
		hdffile.close()
	#if subhalos are written in separate files, then read them too
	if (iseparatesubfiles==1):
		for ifile in range(numfiles):
			if (inompi==True):
				filename=basefilename+".sublevels"+".catalog_parttypes"
				hdffilename=basefilename+".hdf"+".sublevels"+".catalog_parttypes"
			else:
				filename=basefilename+".sublevels"+".catalog_parttypes"+"."+str(ifile)
				hdffilename=basefilename+".hdf"+".sublevels"+".catalog_parttypes"+"."+str(ifile)
			if (iverbose) : print("reading ",filename)
			halofile = open(filename, 'r')
			hdffile=h5py.File(hdffilename,'w')
			[filenum,numfiles]=halofile.readline().split()
			[numhalos, numtothalos]= halofile.readline().split()
			filenum=int(filenum);numfiles=int(numfiles)
			numhalos=np.uint64(numhalos);numtothalos=np.uint64(numtothalos)
			#write header info
			hdffile.create_dataset("File_id",data=np.array([filenum]))
			hdffile.create_dataset("Num_of_files",data=np.array([numfiles]))
			hdffile.create_dataset("Num_of_particles_in_groups",data=np.array([numhalos]));
			hdffile.create_dataset("Total_num_of_particles_in_all_groups",data=np.array([numtothalos]));
			halofile.close()
			if (numhalos>0): cattemp = np.loadtxt(filename,skiprows=2).transpose()
			else: cattemp=[]
			hdffile.create_dataset("Particle_Types",data=np.array(cattemp,dtype=np.int64))
			hdffile.close()

def ConvertASCIIToHDF(basefilename,iseparatesubfiles=0,itype=0,iverbose=0):
	ConvertASCIIPropertyFileToHDF(basefilename,iseparatesubfiles,iverbose)
	ConvertASCIICatalogGroupsFileToHDF(basefilename,iseparatesubfiles,iverbose)
	ConvertASCIICatalogParticleFileToHDF(basefilename,0,iseparatesubfiles,iverbose)
	ConvertASCIICatalogParticleFileToHDF(basefilename,1,iseparatesubfiles,iverbose)
	if (itype==1):
		ConvertASCIICatalogParticleTypeFileToHDF(basefilename,0,iseparatesubfiles,iverbose)
		ConvertASCIICatalogParticleTypeFileToHDF(basefilename,1,iseparatesubfiles,iverbose)
