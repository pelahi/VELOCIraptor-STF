/*! \file nchiladaio.cxx
 *  \brief this file contains routines for nchilada snapshot file io

 */


#include "stf.h"
#include "nchiladaitems.h"
#ifdef USEXDR

///\name XDR read routines, assist in reading nchilada data
//@{
int xdr_header(XDR *pxdrs,struct tipsy_dump *ph)
{
    int pad = 0;
    if (!xdr_double(pxdrs,&ph->time)) return 0;
    if (!xdr_int(pxdrs,&ph->nbodies)) return 0;
    if (!xdr_int(pxdrs,&ph->ndim)) return 0;
    if (!xdr_int(pxdrs,&ph->nsph)) return 0;
    if (!xdr_int(pxdrs,&ph->ndark)) return 0;
    if (!xdr_int(pxdrs,&ph->nstar)) return 0;
    if (!xdr_int(pxdrs,&pad)) return 0;
    return 1;
}

int xdr_gas(XDR *pxdrs,struct tipsy_gas_particle *ph)
{
    int i;
    if (!xdr_float(pxdrs,&ph->mass)) return 0;
    for(i=0; i<NCHILADAMAXDIM; i++) if (!xdr_float(pxdrs,&ph->pos[i])) return 0;
    for(i=0; i<NCHILADAMAXDIM; i++) if (!xdr_float(pxdrs,&ph->vel[i])) return 0;
    if (!xdr_float(pxdrs,&ph->rho)) return 0;
    if (!xdr_float(pxdrs,&ph->temp)) return 0;
    if (!xdr_float(pxdrs,&ph->hsmooth)) return 0;
    if (!xdr_float(pxdrs,&ph->metals)) return 0;
    if (!xdr_float(pxdrs,&ph->phi)) return 0;
    return 1;
}

int xdr_dark(XDR *pxdrs,struct tipsy_dark_particle *ph)
{
    int i;
    if (!xdr_float(pxdrs,&ph->mass)) return 0;
    for(i=0; i<NCHILADAMAXDIM; i++) if (!xdr_float(pxdrs,&ph->pos[i])) return 0;
    for(i=0; i<NCHILADAMAXDIM; i++) if (!xdr_float(pxdrs,&ph->vel[i])) return 0;
    if (!xdr_float(pxdrs,&ph->eps)) return 0;
    if (!xdr_float(pxdrs,&ph->phi)) return 0;
    return 1;
}

int xdr_star(XDR *pxdrs,struct tipsy_star_particle *ph)
{
    int i;
    if (!xdr_float(pxdrs,&ph->mass)) return 0;
    for(i=0; i<NCHILADAMAXDIM; i++) if (!xdr_float(pxdrs,&ph->pos[i])) return 0;
    for(i=0; i<NCHILADAMAXDIM; i++) if (!xdr_float(pxdrs,&ph->vel[i])) return 0;
    if (!xdr_float(pxdrs,&ph->metals)) return 0;
    if (!xdr_float(pxdrs,&ph->tform)) return 0;
    if (!xdr_float(pxdrs,&ph->eps)) return 0;
    if (!xdr_float(pxdrs,&ph->phi)) return 0;
    return 1;
}

int xdr_NCHeader(XDR *pxdrs,struct nchilada_dump *ph)
{
    if (!xdr_int(pxdrs,&ph->magic)) return 0;
    if (!xdr_double(pxdrs,&ph->time)) return 0;
    if (!xdr_int(pxdrs,&ph->iHighWord)) return 0;
    if (!xdr_int(pxdrs,&ph->nbodies)) return 0;
    if (!xdr_int(pxdrs,&ph->ndim)) return 0;
    if (!xdr_int(pxdrs,&ph->code)) return 0;
    return 1;
}

void xdr_NC_Open(XDR *xdrs, enum xdr_op op, FILE **fp, char *achInFile, long lstart)
{
    if(op == XDR_DECODE) *fp = fopen(achInFile,"r");
    else if(op == XDR_ENCODE && lstart) *fp = fopen(achInFile,"a");
    else if(op == XDR_ENCODE && lstart == 0) *fp = fopen(achInFile,"w");
    assert(*fp);
    fseek(*fp,lstart,0);
    xdrstdio_create(xdrs,*fp,op);
}
/*
int xdr_type(XDR *pxdrs, void *data, int type, int num)
{
    char *pachTmp;
    short *psTmp;
    int *piTmp;
    int iTmp;
    float *pfTmp;
    float fTmp;
    long *plTmp;
    long lTmp;
    double *pdTmp;
    double dTmp;
    int j;

    switch (type) {
        case int8:
            pachTmp = ((char *)data);
            for(j=0; j<num; j++) {
                iTmp = (int)pachTmp[j];
                if(!xdr_int(pxdrs, &iTmp)) return 0;
                }
            break;
        case uint8:
            pachTmp = ((char *)data);
            for(j=0; j<num; j++) {
                iTmp = (int)pachTmp[j];
                if(!xdr_int(pxdrs, &iTmp)) return 0;
                }
            break;
        case int16:
            psTmp = ((short *)data);
            for(j=0; j<num; j++) {
                iTmp = (int)psTmp[j];
                if(!xdr_int(pxdrs, &iTmp)) return 0;
                }
            break;
        case uint16:
            psTmp = ((short *)data);
            for(j=0; j<num; j++) {
                iTmp = (int)psTmp[j];
                if(!xdr_int(pxdrs, &iTmp)) return 0;
                }
            break;
        case int32:
            piTmp = ((int *)data);
            for(j=0; j<num; j++) if(!xdr_int(pxdrs, &piTmp[j])) return 0;
            break;
        case uint32:
            piTmp = ((int *)data);
            for(j=0; j<num; j++) if(!xdr_int(pxdrs, &piTmp[j])) return 0;
            break;
        case int64:
            plTmp = ((long *)data);
            for(j=0; j<num; j++) if(!xdr_long(pxdrs, &plTmp[j])) return 0;
            break;
        case uint64:
            plTmp = ((long *)data);
            for(j=0; j<num; j++) if(!xdr_long(pxdrs, &plTmp[j])) return 0;
            break;
        case float32:
            pfTmp = ((float *)data);
            for(j=0; j<num; j++) if(!xdr_float(pxdrs, &pfTmp[j])) return 0;
            break;
        case float64:
            pdTmp = ((double *)data);
            for(j=0; j<num; j++) if(!xdr_double(pxdrs, &pdTmp[j])) return 0;
            break;
        default:
            return 0;
        }
        return 1;
}
*/

///read a particular data field from a nchilada file
void *readFieldData(FILE *&infile, nchilada_dump &fh, unsigned int dim, u_int64_t numParticles, u_int64_t startParticle)
{
    XDR xdrs;
    xdrstdio_create(&xdrs, infile, XDR_DECODE);
    string message;
    if(!xdr_template(&xdrs, &fh)) {
        //throw XDRException("Couldn't read header from file!");
        message="Couldn't read header from file!";
        cout<<message<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,9);
#else
        exit(9);
#endif
    }
    /*
    if(fh.magic != FieldHeader::MagicNumber) {
        //throw XDRException("This file does not appear to be a field file (magic number doesn't match).");
        message="This file does not appear to be a field file (magic number doesn't match).";
        cout<<message<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,9);
#else
        exit(9);
#endif
    }
    */
    if(fh.ndim != dim) {
        //throw XDRException("Wrong dimension of positions.");
        message="Wrong dimension of positions.";
        cout<<message<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,9);
#else
        exit(9);
#endif
    }

    void* data = readField(fh, &xdrs, numParticles, startParticle);

    if(data == 0) {
        //throw XDRException("Had problems reading in the field");
        message="Had problems reading in the field";
        cout<<message<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,9);
#else
        exit(9);
#endif
    }
    xdr_destroy(&xdrs);
    return data;
}

///read a particular data field from a nchilada file
void *readFieldData(const string filename, nchilada_dump &fh, unsigned int dim, u_int64_t numParticles, u_int64_t startParticle)
{
    FILE* infile;
    string message;
    infile=fopen(filename.c_str(), "rb");
    if(!infile) {
        message="Couldn't open field file: ";
        message += filename;
        cout<<message<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,9);
#else
        exit(9);
#endif
        //throw XDRException(message);
    }
    void *data=readFieldData(infile,fh,dim,numParticles,startParticle);
    fclose(infile);
    return data;
}

/// Returns total number of particles in a given file
/// @param filename data file of interest. Criticially returns 0 if this file cannot be read
/// interpreted as no particles
Int_t ncGetCount(string filename)
{
    FILE* infile = fopen(filename.c_str(), "rb");
    if(!infile) {
        return 0;  // Assume there is none of this particle type
    }

    XDR xdrs;
    nchilada_dump fh;
    xdrstdio_create(&xdrs, infile, XDR_DECODE);

    if(!xdr_template(&xdrs, &fh)) {
        //throw XDRException("Couldn't read header from file!");
        cout<<"Couldn't read header from file!"<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,9);
#else
        exit(9);
#endif
    }
    if(fh.ndim != 3 && fh.ndim != 1) {
        //throw XDRException("Wrong dimension.");
        cout<<"Wrong dimension!"<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,9);
#else
        exit(9);
#endif
    }
    xdr_destroy(&xdrs);
    fclose(infile);
    return fh.nbodies;
}

//@}

///get the number of particles of the desired type
Int_t Nchilada_get_nbodies(char *fname, int ptype, Options &opt)
{
    Int_t j,nbodies,ndark,nsph,nstar;
    string filename(fname);
    /*
    int nusetypes,nbusetypes;
    int usetypes[NNCHILADATYPE];
    if (opt.partsearchtype==PSTALL) {
        //lets assume there are dm/stars/gas.
        nusetypes=3;
        usetypes[0]=NCHILADAGASTYPE;usetypes[1]=NCHILADADMTYPE;usetypes[2]=NCHILADASTARTYPE;
    }
    else if (opt.partsearchtype==PSTDARK) {
        nusetypes=1;usetypes[0]=NCHILADADMTYPE;
        if (opt.iBaryonSearch) {
            nbusetypes=2;usetypes[1]=NCHILADAGASTYPE;usetypes[2]=NCHILADASTARTYPE;
        }
    }
    else if (opt.partsearchtype==PSTGAS) {nusetypes=1;usetypes[0]=NCHILADAGASTYPE;}
    else if (opt.partsearchtype==PSTSTAR) {nusetypes=1;usetypes[0]=NCHILADASTARTYPE;}
    for(j=0, nbodies=0; j<nusetypes; j++) {
        k=usetypes[j];
        nbodies+=ncGetCount(filename + nchilada_part_name.part_names[usetypes[j]]+string("pos"));
    }
    */
    nsph = ncGetCount(filename + "/gas/pos");
    ndark = ncGetCount(filename + "/dark/pos");
    nstar = ncGetCount(filename + "/star/pos");
    if (opt.partsearchtype==PSTALL) nbodies=nsph+ndark+nstar;
    else if (opt.partsearchtype==PSTDARK) nbodies=ndark;
    else if (opt.partsearchtype==PSTGAS) nbodies=nsph;
    else if (opt.partsearchtype==PSTSTAR) nbodies=nstar;
    if (nbodies==0) {
        cout<<"Error. Zero particles of type "<<opt.partsearchtype<<" found. Either 0 nor can't find file"<<filename<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,9);
#else
        exit(9);
#endif
    }
    return nbodies;
}

///reads an nchilada formatted file.
void ReadNchilada(Options &opt, vector<Particle> &Part, const Int_t nbodies,Particle *&Pbaryons, Int_t nbaryons)
{
    Int_t i,j,k,n,nchunk;
    u_int64_t numParticles=nbodies,startParticle=0;
    //counters
    Int_t count,countsph,countstar,countbh,count2,bcount,bcount2,pc,pc_new, Ntot,indark,ingas,instar,inbh,Ntotfile;
    Int_t ntot_withmasses;
    char buf[2000];
    char DATA[5];
    string message;
    //store cosmology
    double z,aadjust,Hubble,Hubbleflow;
    Double_t mscale,lscale,lvscale;
    Double_t MP_DM=MAXVALUE,LN,N_DM,MP_B=0;
    int ifirstfile=0,*ireadfile,ireaderror=0;
    int *ireadtask,*readtaskID;
    int nusetypes,nbusetypes;
    int usetypes[NNCHILADATYPE];
    Nchilada_Part_Names nchilada_part_name;

    //common attributes
    /*
    fstream fpos[NNCHILADATYPE],fvel[NNCHILADATYPE],fmass[NNCHILADATYPE];
    fstream fgasu,fgasz,fgassfr;
    fstream fstartage,fstarz;
    */
    FILE *fpos[NNCHILADATYPE],*fvel[NNCHILADATYPE],*fmass[NNCHILADATYPE],*fid[NNCHILADATYPE];
    FILE *fgasu,*fgasz,*fgassfr;
    FILE *fstartage,*fstarz;
    nchilada_dump fhpos,fhvel,fhmass,fhid,fhgasu,fhgasz,fhgassfr,fhstartage,fhstarz;
    void *posdata,*veldata,*massdata,*iddata,*gasudata,*gaszdata,*gassfrdata,*startagedata,*starzdata;
    Double_t tageval;
    int *intbuff;
    long long *longbuff;
    unsigned int *uintbuff;
    unsigned long long *ulongbuff;
    float *posfloatbuff,*velfloatbuff,*massfloatbuff;
    float *gasufloatbuff,*gassfrfloatbuff,*gaszfloatbuff;
    float *starzfloatbuff,*startagefloatbuff;
    double *posdoublebuff,*veldoublebuff,*massdoublebuff;
    double *gasudoublebuff,*gassfrdoublebuff,*gaszdoublebuff;
    double *starzdoublebuff,*startagedoublebuff;

    if (opt.partsearchtype==PSTALL) {
        //lets assume there are dm/stars/gas.
        nusetypes=3;
        usetypes[0]=NCHILADAGASTYPE;usetypes[1]=NCHILADADMTYPE;usetypes[2]=NCHILADASTARTYPE;
    }
    else if (opt.partsearchtype==PSTDARK) {
        nusetypes=1;usetypes[0]=NCHILADADMTYPE;
        if (opt.iBaryonSearch) {
            nbusetypes=2;usetypes[1]=NCHILADAGASTYPE;usetypes[2]=NCHILADASTARTYPE;
        }
    }
    else if (opt.partsearchtype==PSTGAS) {nusetypes=1;usetypes[0]=NCHILADAGASTYPE;}
    else if (opt.partsearchtype==PSTSTAR) {nusetypes=1;usetypes[0]=NCHILADASTARTYPE;}


#ifndef USEMPI
    Int_t Ntotal;
    int ThisTask=0,NProcs=1;
    ireadfile=new int[opt.num_files];
    for (i=0;i<opt.num_files;i++) ireadfile[i]=1;
#endif

    //if verbose spit out the types of particles that are going to be searched for
    if (ThisTask==0 && opt.iverbose>1) {
        cout<<" --------------- "<<endl;
        cout<<"Expecting "<<nusetypes<<" types of particles to be read "<<endl;
        //for (i=0;i<nusetypes;i++) cout<<"Particle "<<usetypes[i]<<" with name "<<hdf_gnames.part_names[usetypes[i]]<<endl;
        if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
            cout<<"Additionally, as full separate baryon search , expecting "<<nbusetypes<<" baryon particles"<<endl;
            //for (i=1;i<=nbusetypes;i++) cout<<"Particle "<<usetypes[i]<<" with name "<<hdf_gnames.part_names[usetypes[i]]<<endl;
        }
    }

    //if MPI is used, read processors (all tasks with task numbers less than the number of snapshots) opens the file and loads the data into a particle buffer
    //this particle buffer is used to broadcast data to the appropriate processor
#ifdef USEMPI
    //since positions, velocities, masses are all at different points in the file,
    //to correctly assign particle to proccessor with correct velocities and mass must have several file pointers
    MPI_Status status;
    Particle *Pbuf;
    int mpi_ireaderror;

    //for parallel io
    Int_t BufSize=opt.mpiparticlebufsize;
    Int_t Nlocalbuf,ibuf=0,*Nbuf, *Nreadbuf,*nreadoffset;
    Int_t *Nlocalthreadbuf,Nlocaltotalbuf;
    int *irecv, sendTask,recvTask,irecvflag, *mpi_irecvflag;
    MPI_Request *mpi_request;
    Int_t *mpi_nsend_baryon;
    if (opt.iBaryonSearch) mpi_nsend_baryon=new Int_t[NProcs*NProcs];

    //extra blocks to store info
    /*
    float *velfloatbuff=new float[NCHILADACHUNKSIZE*3];
    double *veldoublebuff=new double[NCHILADACHUNKSIZE*3];
    float *massfloatbuff=new float[NCHILADACHUNKSIZE];
    double *massdoublebuff=new double[NCHILADACHUNKSIZE];
    float *ufloatbuff=new float[NCHILADACHUNKSIZE];
    double *udoublebuff=new double[NCHILADACHUNKSIZE];
    float *Zfloatbuff=new float[NCHILADACHUNKSIZE];
    double *Zdoublebuff=new double[NCHILADACHUNKSIZE];
    float *SFRfloatbuff=new float[NCHILADACHUNKSIZE];
    double *SFRdoublebuff=new double[NCHILADACHUNKSIZE];
    float *Tagefloatbuff=new float[NCHILADACHUNKSIZE];
    double *Tagedoublebuff=new double[NCHILADACHUNKSIZE];
    */
    Nbuf=new Int_t[NProcs];
    for (int j=0;j<NProcs;j++) Nbuf[j]=0;
    nreadoffset=new Int_t[opt.nsnapread];
    ireadtask=new int[NProcs];
    readtaskID=new int[opt.nsnapread];
    MPIDistributeReadTasks(opt,ireadtask,readtaskID);

    if (ThisTask==0) cout<<"There are "<<opt.nsnapread<<" threads reading "<<opt.num_files<<" files "<<endl;
    if (ireadtask[ThisTask]>=0)
    {
        //to temporarily store data from gadget file
        Pbuf=new Particle[BufSize*NProcs];
        Nreadbuf=new Int_t[opt.num_files];
        for (int j=0;j<opt.num_files;j++) Nreadbuf[j]=0;

        //to determine which files the thread should read
        ireadfile=new int[opt.num_files];
        ifirstfile=MPISetFilesRead(opt,ireadfile,ireadtask);
    }
    else {
        Nlocalthreadbuf=new Int_t[opt.nsnapread];
        irecv=new int[opt.nsnapread];
        mpi_irecvflag=new int[opt.nsnapread];
        for (i=0;i<opt.nsnapread;i++) irecv[i]=1;
        mpi_request=new MPI_Request[opt.nsnapread];
    }
    Nlocal=0;
    if (opt.iBaryonSearch) Nlocalbaryon[0]=0;

    //now begin reading
    if (ireadtask[ThisTask]>=0) {
#endif
    //now for each file lets open up the fstream file pointers
    /*
    for (j=0;j<nusetypes;j++) {
        fpos[j].open((string(opt.fname)+ nchilada_part_name.part_names[usetypes[j]]+string("pos")).c_str());
        fvel[j].open((string(opt.fname)+ nchilada_part_name.part_names[usetypes[j]]+string("vel")).c_str());
        fmass[j].open((string(opt.fname)+ nchilada_part_name.part_names[usetypes[j]]+string("mass")).c_str());
    }
    if (opt.partsearchtype==PSTALL || (opt.partsearchtype==PSTDARK && opt.iBaryonSearch>=1) || opt.partsearchtype==PSTGAS) {
        fgasu.open((string(opt.fname)+ nchilada_part_name.part_names[NCHILADAGASTYPE]+string("U")).c_str());
        fgassfr.open((string(opt.fname)+ nchilada_part_name.part_names[NCHILADAGASTYPE]+string("SFR")).c_str());
        fgasz.open((string(opt.fname)+ nchilada_part_name.part_names[NCHILADAGASTYPE]+string("Z")).c_str());
    }
    if (opt.partsearchtype==PSTALL || (opt.partsearchtype==PSTDARK && opt.iBaryonSearch>=1) || opt.partsearchtype==PSTGAS) {
        fstartage.open((string(opt.fname)+ nchilada_part_name.part_names[NCHILADASTARTYPE]+string("U")).c_str());
        fstarz.open((string(opt.fname)+ nchilada_part_name.part_names[NCHILADASTARTYPE]+string("Z")).c_str());
    }
    */
    for (j=0;j<nusetypes;j++) {
        fpos[j]=fopen((string(opt.fname)+ nchilada_part_name.part_names[usetypes[j]]+string("pos")).c_str(), "rb");
        fvel[j]=fopen((string(opt.fname)+ nchilada_part_name.part_names[usetypes[j]]+string("vel")).c_str(), "rb");
        fmass[j]=fopen((string(opt.fname)+ nchilada_part_name.part_names[usetypes[j]]+string("mass")).c_str(), "rb");
    }
#ifdef GASON
    if (opt.partsearchtype==PSTALL || (opt.partsearchtype==PSTDARK && opt.iBaryonSearch>=1) || opt.partsearchtype==PSTGAS) {
        fgasu=fopen((string(opt.fname)+ nchilada_part_name.part_names[NCHILADAGASTYPE]+string("U")).c_str(), "rb");
#ifdef STARON
        fgassfr=fopen((string(opt.fname)+ nchilada_part_name.part_names[NCHILADAGASTYPE]+string("SFR")).c_str(), "rb");
        fgasz=fopen((string(opt.fname)+ nchilada_part_name.part_names[NCHILADAGASTYPE]+string("Z")).c_str(), "rb");
#endif
    }
#endif
#if defined(STARON) || defined(BHON)
    if (opt.partsearchtype==PSTALL || (opt.partsearchtype==PSTDARK && opt.iBaryonSearch>=1) || opt.partsearchtype==PSTGAS) {
        fstartage=fopen((string(opt.fname)+ nchilada_part_name.part_names[NCHILADASTARTYPE]+string("U")).c_str(), "rb");
        fstarz=fopen((string(opt.fname)+ nchilada_part_name.part_names[NCHILADASTARTYPE]+string("Z")).c_str(), "rb");
    }
#endif
    //will need to think how best to read stuff in chunks

    for (j=0;j<nusetypes;j++) {
        k=usetypes[j];
        posdata=readFieldData(fpos[j],fhpos, 3, numParticles,startParticle);
        if (fhpos.code==float32) posfloatbuff=(float*)posdata;
        else posdoublebuff=(double*)posdata;
        veldata=readFieldData(fvel[j],fhvel, 3, numParticles,startParticle);
        if (fhvel.code==float32) velfloatbuff=(float*)veldata;
        else veldoublebuff=(double*)veldata;
        massdata=readFieldData(fmass[j],fhmass, 1, numParticles,startParticle);
        if (fhmass.code==float32) massfloatbuff=(float*)massdata;
        else massdoublebuff=(double*)massdata;
        iddata=readFieldData(fid[j],fhid, 1, numParticles,startParticle);
        if (fhid.code==int32) intbuff=(int*)iddata;
        else if(fhid.code==int64) longbuff=(long long*)iddata;
        else if(fhid.code==uint32) uintbuff=(unsigned int*)iddata;
        else if (fhid.code==uint64) ulongbuff=(unsigned long long*)iddata;
#ifdef GASON
        gasudata=readFieldData(fgasu,fhgasu, 1, numParticles,startParticle);
        if (fhgasu.code==float32) gasufloatbuff=(float*)gasudata;
        else gasudoublebuff=(double*)gasudata;
#ifdef STARON
        if (usetypes[j]==NCHILADAGASTYPE) {
            gaszdata=readFieldData(fgasz,fhgasz, 1, numParticles,startParticle);
            if (fhgasz.code==float32) gaszfloatbuff=(float*)gaszdata;
            else gaszdoublebuff=(double*)gaszdata;
            gassfrdata=readFieldData(fgassfr,fhgassfr, 1, numParticles,startParticle);
            if (fhgassfr.code==float32) gassfrfloatbuff=(float*)gassfrdata;
            else gassfrdoublebuff=(double*)gassfrdata;
        }
#endif
#endif
#ifdef STARON
        if (usetypes[j]==NCHILADASTARTYPE) {
            starzdata=readFieldData(fstarz,fhstarz, 1, numParticles,startParticle);
            if (fhstarz.code==float32) starzfloatbuff=(float*)starzdata;
            else starzdoublebuff=(double*)starzdata;
            startagedata=readFieldData(fstartage,fhstartage, 1, numParticles,startParticle);
            if (fhstartage.code==float32) startagefloatbuff=(float*)startagedata;
            else startagedoublebuff=(double*)startagedata;
        }
#endif
        for (i=0;i<fhpos.nbodies;i++) {
#ifdef USEMPI
            if (fhpos.code==float32) ibuf=MPIGetParticlesProcessor(posfloatbuff[i],posfloatbuff[i+nbodies],posfloatbuff[i+nbodies*2]);
            else ibuf=MPIGetParticlesProcessor(posdoublebuff[i],posdoublebuff[i+nbodies],posdoublebuff[i+nbodies*2]);
            if (fhpos.code==float32) {
                Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPosition(posfloatbuff[i],posfloatbuff[i+nbodies],posfloatbuff[i+nbodies*2]);
                Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetVelocity(velfloatbuff[i],velfloatbuff[i+nbodies],velfloatbuff[i+nbodies*2]);
                Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetMass(massfloatbuff[i]);
            }
            else {
                Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPosition(posdoublebuff[i],posdoublebuff[i+nbodies],posdoublebuff[i+nbodies*2]);
                Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetVelocity(veldoublebuff[i],veldoublebuff[i+nbodies],veldoublebuff[i+nbodies*2]);
                Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetMass(massdoublebuff[i]);
            }
            if (fhid.code==uint32) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPID(intbuff[i]);
            else Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPID(longbuff[i]);
            Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetID(i);
            if (k==NCHILADAGASTYPE) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetType(GASTYPE);
            else if (k==NCHILADADMTYPE) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetType(DARKTYPE);
            else if (k==NCHILADASTARTYPE) {
#ifdef BHON
                if (fhstartage.code==float32) tageval=startagefloatbuff[i];
                else tageval=startagedoublebuff[i];
                if (tageval>0) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetType(STARTYPE);
                else Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetType(BHTYPE);
#else
                Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetType(STARTYPE);
#endif
            }
#ifdef GASON
            if (k==NCHILADAGASTYPE) {
                if (fhgasu.code==float32) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetU(gasufloatbuff[i]);
                else Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetU(gasudoublebuff[i]);
#ifdef STARON
                if (fhgassfr.code==float32) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetSFR(gassfrfloatbuff[i]);
                else Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetSFR(gassfrdoublebuff[i]);
                if (fhgasz.code==float32) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetZmet(gaszfloatbuff[i]);
                else Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetZmet(gaszdoublebuff[i]);
#endif
            }
#endif
#ifdef STARON
            if (k==NCHILADASTARTYPE) {
                if (fhstartage.code==float32) tageval=startagefloatbuff[i];
                else tageval=startagedoublebuff[i];
#ifdef BHON
                if (tageval<0) tageval*=-1;
#endif
                if (fhstarz.code==float32) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetZmet(starzfloatbuff[i]);
                else Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetZmet(starzdoublebuff[i]);
                Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetTage(tageval);
            }
#endif
            if(ireadtask[ibuf]>=0&&ibuf!=ThisTask) Nreadbuf[ireadtask[ibuf]]++;
            Nbuf[ibuf]++;
            if (ibuf==ThisTask) {
                Nbuf[ibuf]--;
                Part[Nlocal++]=Pbuf[ibuf*BufSize+Nbuf[ibuf]];
            }
            else {
                //before a simple send was done because only Task zero was reading the data
                //but now if ibuf<opt.nsnapread, care must be taken.
                //blocking sends that are matched by non-blocking receives
                if(Nbuf[ibuf]==BufSize&&ireadtask[ibuf]<0) {
                    MPI_Send(&Nbuf[ibuf], 1, MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
                    MPI_Send(&Pbuf[ibuf*BufSize],sizeof(Particle)*Nbuf[ibuf],MPI_BYTE,ibuf,ibuf,MPI_COMM_WORLD);
                    Nbuf[ibuf]=0;
                }
                else if (Nbuf[ibuf]==BufSize&&ireadtask[ibuf]>=0) {
                    Nbuf[ibuf]=0;
                }
            }
//start of non mpi section
#else
            if (fhpos.code==float32) {
                Part[i].SetPosition(posfloatbuff[i],posfloatbuff[i+nbodies],posfloatbuff[i+nbodies*2]);
                Part[i].SetVelocity(velfloatbuff[i],velfloatbuff[i+nbodies],velfloatbuff[i+nbodies*2]);
                Part[i].SetMass(massfloatbuff[i]);
            }
            else {
                Part[i].SetPosition(posdoublebuff[i],posdoublebuff[i+nbodies],posdoublebuff[i+nbodies*2]);
                Part[i].SetVelocity(veldoublebuff[i],veldoublebuff[i+nbodies],veldoublebuff[i+nbodies*2]);
                Part[i].SetMass(massdoublebuff[i]);
            }
            if (fhid.code==uint32) Part[i].SetPID(intbuff[i]);
            else Part[i].SetPID(longbuff[i]);
            Part[i].SetID(i);
            if (k==NCHILADAGASTYPE) Part[i].SetType(GASTYPE);
            else if (k==NCHILADADMTYPE) Part[i].SetType(DARKTYPE);
            else if (k==NCHILADASTARTYPE) {
#ifdef BHON
                if (fhstartage.code==float32) tageval=startagefloatbuff[i];
                else tageval=startagedoublebuff[i];
                if (tageval>0) Part[i].SetType(STARTYPE);
                else Part[i].SetType(BHTYPE);
#else
                Part[i].SetType(STARTYPE);
#endif
            }
#ifdef GASON
            if (k==NCHILADAGASTYPE) {
                if (fhgasu.code==float32) Part[i].SetU(gasufloatbuff[i]);
                else Part[i].SetU(gasudoublebuff[i]);
#ifdef STARON
                if (fhgassfr.code==float32) Part[i].SetSFR(gassfrfloatbuff[i]);
                else Part[i].SetSFR(gassfrdoublebuff[i]);
                if (fhgasz.code==float32) Part[i].SetZmet(gaszfloatbuff[i]);
                else Part[i].SetZmet(gaszdoublebuff[i]);
#endif
            }
#endif
#ifdef STARON
            if (k==NCHILADASTARTYPE) {
                if (fhstartage.code==float32) tageval=startagefloatbuff[i];
                else tageval=startagedoublebuff[i];
#ifdef BHON
                if (tageval<0) tageval*=-1;
#endif
                if (fhstarz.code==float32) Part[i].SetZmet(starzfloatbuff[i]);
                else Part[i].SetZmet(starzdoublebuff[i]);
                Part[i].SetTage(tageval);
            }
#endif
#endif
//end of mpi ifdef
        }//end of loop over particles
    }//end of loop over particle types
    //close files
    for (j=0;j<nusetypes;j++) {
        fclose(fpos[j]);
        fclose(fvel[j]);
        fclose(fmass[j]);
    }
#ifdef GASON
    if (opt.partsearchtype==PSTALL || (opt.partsearchtype==PSTDARK && opt.iBaryonSearch>=1) || opt.partsearchtype==PSTGAS) {
        fclose(fgasu);
#ifdef STARON
        fclose(fgassfr);
        fclose(fgasz);
#endif
    }
#endif
#if defined(STARON) || defined(BHON)
    if (opt.partsearchtype==PSTALL || (opt.partsearchtype==PSTDARK && opt.iBaryonSearch>=1) || opt.partsearchtype==PSTGAS) {
        fclose(fstartage);
        fclose(fstarz);
    }
#endif

#ifdef USEMPI
    }//end of read task section
    else {
        MPIReceiveParticlesFromReadThreads(opt,Pbuf,Part.data(),readtaskID, irecv, mpi_irecvflag, Nlocalthreadbuf, mpi_request,Pbaryons);
    }

    for (i=0;i<NProcs;i++) Nbuf[i]=0;
    if (ireadtask[ThisTask]>=0 && opt.nsnapread>1) {
        delete[] Pbuf;
        Nlocalbuf=0;
        for (i=0;i<opt.nsnapread;i++) Nlocalbuf+=Nreadbuf[i];
        if (Nlocalbuf>0)
        {
            Pbuf=new Particle[Nlocalbuf];
            //determine offsets
            nreadoffset[0]=0;for (i=1;i<opt.nsnapread;i++)nreadoffset[i]=nreadoffset[i-1]+Nreadbuf[i-1];
            for(i=0;i<opt.num_files; i++)
            if (ireadfile[i])
            {
                for (j=0;j<nusetypes;j++) {
                }
            }
        }
        //gather all the items that must be sent.
        MPI_Allgather(Nbuf, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
        //if separate baryon search then sort the Pbuf array so that it is separated by type
        if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
            if (ThisTask<opt.nsnapread) {
            for(ibuf = 0; ibuf < opt.nsnapread; ibuf++) if (mpi_nsend[ThisTask * NProcs + ibuf] > 0)
            {
                Nbuf[ibuf]=0;
                for (i=0;i<mpi_nsend[ThisTask * NProcs + ibuf];i++) {
                    k=Pbuf[nreadoffset[ibuf]+i].GetType();
                    if (!(k==GASTYPE||k==STARTYPE||k==BHTYPE)) Pbuf[nreadoffset[ibuf]+i].SetID(0);
                    else {
                        if  (k==GASTYPE) Pbuf[nreadoffset[ibuf]+i].SetID(1);
                        else if  (k==STARTYPE) Pbuf[nreadoffset[ibuf]+i].SetID(2);
                        else if  (k==BHTYPE) Pbuf[nreadoffset[ibuf]+i].SetID(3);
                        Nbuf[ibuf]++;
                    }
                }
                qsort(&Pbuf[nreadoffset[ibuf]],mpi_nsend[ThisTask*NProcs+ibuf], sizeof(Particle), IDCompare);
            }
            }
            MPI_Allgather(Nbuf, NProcs, MPI_Int_t, mpi_nsend_baryon, NProcs, MPI_Int_t, MPI_COMM_WORLD);
            for (ibuf=0;ibuf<NProcs*NProcs;ibuf++) mpi_nsend[ibuf]-=mpi_nsend_baryon[ibuf];
        }
        //and then send all the data between the read threads
        MPISendParticlesBetweenReadThreads(opt, Pbuf, Part.data(), nreadoffset, ireadtask, readtaskID, Pbaryons, mpi_nsend_baryon);
        if (ireadtask[ThisTask]>=0) {
            delete[] Pbuf;
            if (opt.iBaryonSearch && opt.partsearchtype!=PSTALL) delete[] mpi_nsend_baryon;
            //set IDS
            for (i=0;i<Nlocal;i++) Part[i].SetID(i);
            if (opt.iBaryonSearch) for (i=0;i<Nlocalbaryon[0];i++) Pbaryons[i].SetID(i+Nlocal);
        }//end of read tasks
    }
#endif

    ///if gas found and Omega_b not set correctly (ie: ==0), assumes that
    ///lowest mass gas particle found corresponds to Omega_b
    ///Note that if there is mass evolution this WILL NOT WORK!
    if (opt.Omega_b==0 && MP_B==MAXVALUE){
        opt.Omega_b=MP_B/(MP_DM+MP_B)*opt.Omega_m;
        opt.Omega_cdm=opt.Omega_m-opt.Omega_b;
    }
    //adjust period
    if (opt.comove) opt.p*=opt.lengthinputconversion/opt.h;
    else opt.p*=opt.lengthinputconversion/opt.h*opt.a;
#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
    //update cosmological data and boundary in code units
    MPI_Bcast(&(opt.p),sizeof(opt.p),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.a),sizeof(opt.a),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.Omega_cdm),sizeof(opt.Omega_cdm),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.Omega_b),sizeof(opt.Omega_b),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.Omega_m),sizeof(opt.Omega_m),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.Omega_Lambda),sizeof(opt.Omega_Lambda),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.h),sizeof(opt.h),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.rhocrit),sizeof(opt.rhocrit),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.rhobg),sizeof(opt.rhobg),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.virlevel),sizeof(opt.virlevel),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.virBN98),sizeof(opt.virBN98),MPI_BYTE,0,MPI_COMM_WORLD);
#ifdef NOMASS
    MPI_Bcast(&(opt.MassValue),sizeof(opt.MassValue),MPI_BYTE,0,MPI_COMM_WORLD);
#endif
    MPI_Bcast(&(Ntotal),sizeof(Ntotal),MPI_BYTE,0,MPI_COMM_WORLD);

    MPI_Bcast(&(lscale),sizeof(lscale),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(lvscale),sizeof(lvscale),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(mscale),sizeof(mscale),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(Hubbleflow),sizeof(Hubbleflow),MPI_BYTE,0,MPI_COMM_WORLD);
#endif
    ///If compiled with HIGHRES, the code assumes that the gadget data is a multi-resolution simulation
    ///with the lowest mass dark matter particle corresponding to the highest resolution and
    ///thus the physical linking length is assumed to be in fraction of interparticle spacing
    ///and is adjusted to a physical distance. Note that if high res and neff is not passed
    ///code assumes that lowest mass gas particle can be used to determine Omega_b and thus can be used
    ///to calculate the mean interparticle spacing.
    ///but one can also pass opt.Neff to adjust what the code thinks is the average inter particle spacing
#ifdef HIGHRES
    if (opt.Neff==-1) {
        //Once smallest mass particle is found (which should correspond to highest resolution area,
        if (opt.Omega_b==0) MP_B=0;
        LN=pow(((MP_DM+MP_B)*opt.massinputconversion/opt.h)/(opt.Omega_m*3.0*opt.H*opt.h*opt.H*opt.h/(8.0*M_PI*opt.G)),1./3.)*opt.a;
    }
    else {
        LN=opt.p/(Double_t)opt.Neff;
    }
#endif
#ifdef USEMPI
    MPI_Bcast(&LN, 1, MPI_Real_t, 0, MPI_COMM_WORLD);
#endif
    opt.internalenergyinputconversion = opt.velocityinputconversion*opt.velocityinputconversion;
    ///if not an individual halo, assume cosmological and store scale of the highest resolution interparticle spacing to scale the physical FOF linking length
    if (opt.iSingleHalo==0)
    {
        opt.ellxscale=LN;
        opt.uinfo.eps*=LN;
    }

#ifdef USEMPI
    //finally adjust to appropriate units
    for (i=0;i<Nlocal;i++)
    {
        Part[i].SetMass(Part[i].GetMass()*mscale);
        for (int j=0;j<3;j++) Part[i].SetVelocity(j,Part[i].GetVelocity(j)*opt.velocityinputconversion*sqrt(opt.a)+Hubbleflow*Part[i].GetPosition(j));
        for (int j=0;j<3;j++) Part[i].SetPosition(j,Part[i].GetPosition(j)*lscale);
    }
    if (Pbaryons!=NULL && opt.iBaryonSearch==1) {
    for (i=0;i<Nlocalbaryon[0];i++)
    {
        Pbaryons[i].SetMass(Pbaryons[i].GetMass()*mscale);
        for (int j=0;j<3;j++) Pbaryons[i].SetVelocity(j,Pbaryons[i].GetVelocity(j)*opt.velocityinputconversion*sqrt(opt.a)+Hubbleflow*Pbaryons[i].GetPosition(j));
        for (int j=0;j<3;j++) Pbaryons[i].SetPosition(j,Pbaryons[i].GetPosition(j)*lscale);
    }
    }
#endif

}
#endif
