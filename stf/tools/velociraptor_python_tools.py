import sys,os,os.path,string,time,re,struct
import math,operator
from pylab import *
import numpy as np
import h5py 

"""

Note that this code is compatible with python2. 
For python3 please search for all print statemetns and replace them with print() style calls

"""


def ReadPropertyFile(basefilename,ibinary=0,iseparatesubfiles=0,iverbose=0):
    """
    VELOCIraptor/STF files in ascii format contain 
    a header with 
        filenumber number_of_files
        nhalos_in_file nnhalos_in_total
    followed by a header listing the information contain. An example would be 
        ID(1) ID_mbp(2) hostHaloID(3) numSubStruct(4) npart(5) Mvir(6) Xc(7) Yc(8) Zc(9) Xcmbp(10) Ycmbp(11) Zcmbp(12) VXc(13) VYc(14) VZc(15) VXcmbp(16) VYcmbp(17) VZcmbp(18) Mass_tot(19) Mass_FOF(20) Mass_200mean(21) Mass_200crit(22) Mass_BN97(23) Efrac(24) Rvir(25) R_size(26) R_200mean(27) R_200crit(28) R_BN97(29) R_HalfMass(30) Rmax(31) Vmax(32) sigV(33) veldisp_xx(34) veldisp_xy(35) veldisp_xz(36) veldisp_yx(37) veldisp_yy(38) veldisp_yz(39) veldisp_zx(40) veldisp_zy(41) veldisp_zz(42) lambda_B(43) Lx(44) Ly(45) Lz(46) q(47) s(48) eig_xx(49) eig_xy(50) eig_xz(51) eig_yx(52) eig_yy(53) eig_yz(54) eig_zx(55) eig_zy(56) eig_zz(57) cNFW(58) Krot(59) Ekin(60) Epot(61) n_gas(62) M_gas(63) Xc_gas(64) Yc_gas(65) Zc_gas(66) VXc_gas(67) VYc_gas(68) VZc_gas(69) Efrac_gas(70) R_HalfMass_gas(71) veldisp_xx_gas(72) veldisp_xy_gas(73) veldisp_xz_gas(74) veldisp_yx_gas(75) veldisp_yy_gas(76) veldisp_yz_gas(77) veldisp_zx_gas(78) veldisp_zy_gas(79) veldisp_zz_gas(80) Lx_gas(81) Ly_gas(82) Lz_gas(83) q_gas(84) s_gas(85) eig_xx_gas(86) eig_xy_gas(87) eig_xz_gas(88) eig_yx_gas(89) eig_yy_gas(90) eig_yz_gas(91) eig_zx_gas(92) eig_zy_gas(93) eig_zz_gas(94) Krot_gas(95) T_gas(96) Zmet_gas(97) SFR_gas(98) n_star(99) M_star(100) Xc_star(101) Yc_star(102) Zc_star(103) VXc_star(104) VYc_star(105) VZc_star(106) Efrac_star(107) R_HalfMass_star(108) veldisp_xx_star(109) veldisp_xy_star(110) veldisp_xz_star(111) veldisp_yx_star(112) veldisp_yy_star(113) veldisp_yz_star(114) veldisp_zx_star(115) veldisp_zy_star(116) veldisp_zz_star(117) Lx_star(118) Ly_star(119) Lz_star(120) q_star(121) s_star(122) eig_xx_star(123) eig_xy_star(124) eig_xz_star(125) eig_yx_star(126) eig_yy_star(127) eig_yz_star(128) eig_zx_star(129) eig_zy_star(130) eig_zz_star(131) Krot_star(132) tage_star(133) Zmet_star(134) 

    then followed by data
    
    Note that a file will indicate how many files the total output has been split into

    """
    #this variable is the size of the char array in binary formated data that stores the field names
    CHARSIZE=40
    #define a data type that matches the binary data structure written by io.cxx of velociraptor
    dt = []

    start = time.clock()
    inompi=True
    if (iverbose): print "reading properties file",basefilename,os.path.isfile(basefilename)
    filename=basefilename+".properties"
    #load header
    if (os.path.isfile(basefilename)==True):
        numfiles=0
    else:
        filename=basefilename+".properties"+".0"
        inompi=False
        if (os.path.isfile(filename)==False):
            print "file not found"
            return []
    byteoffset=0
    #if ascii, binary or hdf5, first get field names
    if (ibinary==0):
        #load ascii file
        halofile = open(filename, 'r')
        #read header information
        [filenum,numfiles]=halofile.readline().split()
        filenum=int(filenum);numfiles=int(numfiles)
        [numhalos, numtothalos]= halofile.readline().split()
        numhalos=int(numhalos);numtothalos=int(numtothalos)
        names = ((halofile.readline())).split()
        #remove the brackets in ascii file names 
        fieldnames= [fieldname.split("(")[0] for fieldname in names]
        halofile.close()
    elif (ibinary==1):
        #load binary file
        halofile = open(filename, 'rb')
        [filenum,numfiles]=np.fromfile(halofile,dtype=np.int32,count=2)
        [numhalos,numtothalos]=np.fromfile(halofile,dtype=np.uint64,count=2)
        headersize=np.fromfile(halofile,dtype=np.int32,count=1)[0]
        fieldnames=[]
        byteoffset=np.dtype(np.int32).itemsize*3+np.dtype(np.uint64).itemsize*2+4*headersize
        for i in range(headersize):
            fieldnames.append(unpack('s', halofile.read(CHARSIZE)).strip())
        for i in np.arange(fieldnames.__len__()):
            if fieldname in ["ID","hostHalo","numSubStruct","npart","n_gas","n_star"]:
                dt.append((fieldname,np.uint64))
            else:
                dt.append((fieldname,np.float64))
        halofile.close()
    elif (ibinary==2):
        #load hdf file
        halofile = h5py.File(filename, 'r')
        filenum=int(halofile["File_id"][0])
        numfiles=int(halofile["Num_of_files"][0])
        numhalos=np.int64(halofile["Num_of_groups"][0])
        numtothalos=np.int64(halofile["Total_num_of_groups"][0])
        fieldnames=[str(n) for n in halofile.keys()]
        #clean of header info
        fieldnames.remove("File_id")
        fieldnames.remove("Num_of_files")
        fieldnames.remove("Num_of_groups")
        fieldnames.remove("Total_num_of_groups")
        
    #allocate memory that will store the halo dictionary
    catalog = dict()
    halos=[np.zeros(numtothalos) for catvalue in fieldnames]
    noffset=0
    for ifile in range(numfiles):
        if (inompi==True): filename=basefilename+".properties"
        else: filename=basefilename+".properties"+"."+str(ifile)
        if (iverbose) : print "reading ",filename
        if (ibinary==0): 
            htemp = np.loadtxt(filename,skiprows=3).transpose()
        elif(ibinary==1):
            halofile = open(filename, 'rb')
            halofile.seek(byteoffset);
            htemp=np.fromfile(halofile,dtype=dt).transpose()
        elif(ibinary==2):
            #here convert the hdf information into a numpy array
            htemp=[np.array(halofile[catvalue]) for catvalue in fieldnames]
        numhalos=len(htemp[0])
        for i in range(len(fieldnames)):
            catvalue=fieldnames[i]
            halos[i][noffset:noffset+numhalos]=htemp[i]
        noffset+=numhalos
    #if subhalos are written in separate files, then read them too
    if (iseparatesubfiles==1):
        for ifile in range(numfiles):
            if (inompi==True): filename=basefilename+".sublevels"+".properties"
            else: filename=basefilename+".sublevels"+".properties"+"."+str(ifile)
            if (iverbose) : print "reading ",filename
            if (ibinary==0): 
                htemp = np.loadtxt(filename,skiprows=3).transpose()
            elif(ibinary==1):
                halofile = open(filename, 'rb')
                halofile.seek(byteoffset);
                htemp=np.fromfile(halofile,dtype=dt).transpose()
            elif(ibinary==2):
                halofile = h5py.File(filename, 'r')
                htemp=[np.array(halofile[catvalue]) for catvalue in fieldnames]
            numhalos=len(htemp[0])
            for i in range(len(fieldnames)):
                catvalue=fieldnames[i]
                halos[i][noffset:noffset+numhalos]=htemp[i]
            noffset+=numhalos        

    #set dictionary
    for i in np.arange(fieldnames.__len__()):
        fieldname = fieldnames[i]
        catalog[fieldname] = halos[i]
        #only correct data format for ascii loading
        if (ibinary==0):
            if fieldname in ["ID","hostHalo","numSubStruct","npart","n_gas","n_star"]:
                catalog[fieldname] = int64(catalog[fieldname].round())
            
    if (iverbose): print "done reading properties file ",time.clock()-start
    return catalog,numtothalos

def ReadHaloMergerTree(treefilename,ibinary=0,iverbose=0):
    """
    VELOCIraptor/STF merger tree in ascii format contains 
    a header with 
        number_of_snapshots
        a description of how the tree was built
        total number of halos across all snapshots

    then followed by data
    for each snapshot 
        snapshotvalue nhalos
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
    if (iverbose): print "reading Tree file",treefilename,os.path.isfile(treefilename)
    if (os.path.isfile(treefilename)==False):
        print "Error, file not found"
        return []
    #if ascii format
    if (ibinary==0):
        treefile = open(treefilename, 'r')
        numsnap=int(treefile.readline())
        descrip=treefile.readline().strip()
        tothalos=int(treefile.readline())
        tree=[{"haloID" , "Num_progen", "Progen"} for i in range(numsnap)]
        offset=0
        for i in range(numsnap):
            [snapval,numhalos]=treefile.readline().strip().split(' ')
            snapval=int(snapval);numhalos=int(numhalos)
            tree[i]["haloID"]=np.zeros(numhalos, dtype=np.int64)
            tree[i]["Num_progen"]=np.zeros(numhalos, dtype=np.int32)
            tree[i]["Progen"]=[[] for j in range(numhalos)]
            for j in range(numhalos):
                [hid,nprog]=treefile.readline().strip().split(' ')
                hid=np.int64(hid);nprog=int(nprog)
                tree[i]["haloID"][j]=hid
                tree[i]["Num_progen"][j]=nprog
                if (nprog>0):
                    tree[i]["Progen"][j]=np.zeros(nprog,dtype=np.int64)
                    for k in range(nprog):
                        tree[i]["Progen"][j][k]=np.int64(treefile.readline())
    return tree            
            
        
        
    

def ReadCrossCatalogList(fname,meritlim=0.1,iverbose=0):
    """
    Reads a cross catalog produced by halomergertree, 
    also allows trimming of cross catalog using a higher merit threshold than one used to produce catalog
    """
    start = time.clock()
    if (iverbose): print "reading cross catalog"
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
    if (iverbose): print "done reading cross catalog ",time.clock()-start
    return pdata

def BuildHierarchy(halodata,iverbose=0):
    """
    the halo data stored in a velociraptor .properties file should store the id of its parent halo. Here 
    this catalog is used to produce a hierarchy to quickly access the relevant subhaloes of a parent halo.
    We could 
    """
    halohierarchy=[]
    start=time.clock()
    if (iverbose): print "setting hierarchy"
    nhalos=len(halodata["npart"])
    subhaloindex=np.where(halodata["hostHaloID"]!=-1)
    lensub=len(subhaloindex[0])
    haloindex=np.where(halodata["hostHaloID"]==-1)
    lenhal=len(haloindex[0])
    halohierarchy=[[] for k in range(nhalos)]
    if (iverbose): print "prelims done ",time.clock()-start
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
    if (iverbose): print "hierarchy set in read in ",time.clock()-start
    return halohierarchy


