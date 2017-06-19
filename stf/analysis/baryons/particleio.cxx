/*! \file particleio.cxx
 *  \brief this file contains routines for reading in particle data
 */

#include "baryoniccontent.h"

#include "../../src/gadgetitems.h"
#include "../../src/tipsy_structs.h"
#include "../../src/endianutils.h"

///\name Read particle data files
//@{

///map ids to index, here for now coded such that index=i-1
//inline UInt_t MapPIDStoIndex(Options &opt, GADGETIDTYPE i){return i-1;}
inline UInt_t MapPIDStoIndex(Options &opt, GADGETIDTYPE i){return i;}

///Reads the header structure if its a tipsy file or get_nbodies if gadget
Int_t ReadHeader(Options &opt){
    InitEndian();
    if (opt.nfiles<1) {
        struct tipsy_dump tipsyheader;
        fstream Ftip(opt.fname, ios::in | ios::binary);
        if (!Ftip){cerr<<"ERROR: Unable to open " <<opt.fname<<endl;exit(8);}
        else cout<<"Reading tipsy format from "<<opt.fname<<endl;
        Ftip.read((char*)&tipsyheader,sizeof(tipsy_dump));
        return tipsyheader.nbodies;
    }
    else {
        return get_nbodies(opt.fname);
    }
}



///Reads particle data
///To add a new interface simply alter this to include the appropriate user written call
void ReadParticleData(Options &opt, const Int_t nbodies, Particle *Part, const Int_t ngroups, HaloParticleData *hp)
{
    InitEndian();
    if(opt.nfiles<1) ReadTipsy(opt,nbodies,Part,ngroups,hp);
    else ReadGadget(opt,nbodies,Part,ngroups,hp);
}

///reads a tipsy file, not implemented yet
void ReadTipsy(Options &opt, const Int_t nbodies, Particle *Part, const Int_t ngroups, HaloParticleData *hp)
{
}

/*! reads a gadget file.
    the approach is to load all the ids first, sort them then determine the index
    of the particle

    Note that gadget files contained
    header
    positions
    velocities
    ids
    masses (if necessary)
    for sph particles :
      internal energies u
      density
      if (cooling)
        mean molecular weight
        neutral hydrogen abundance
      smoothing lengths
      and further blocks of information, depending on included physics
*/
void ReadGadget(Options &opt, const Int_t nbodies, Particle *Part, const Int_t ngroups, HaloParticleData *hp)
{
    Int_t i,j,k,n,m,temp,count,count2,countsph,pc,pc_new,indark,ingas,instar,Ntot,Ntotfile;
    GADGETIDTYPE idval;
    Int_t indexval,groupval;
    Int_t ntot_withmasses;
    unsigned int dummy;
    double z,aadjust,Hubble,HubbleFlow;
    FLOAT ctemp[3],vtemp[3],sphtemp;
    REAL dtemp;
    char buf[200];
    char DATA[5];
    fstream *Fgad, *Fgadvel, *Fgadid, *Fgadmass, *Fgadsphdata;
    struct gadget_header *header;
    Double_t mscale,lscale,lvscale;
    Double_t MP_DM=MAXVALUE,LN;
    Int_t *pfof,*ncount;
    int ifirstfile=0;
#ifndef USEMPI
    //Int_t Nlocal;
    int ThisTask=0,NProcs=1;
#endif
    Int_t Nlocal,Ntotal;

    //if MPI is used, Processor zero opens the file and loads the data into a particle buffer
    //this particle buffer is used to broadcast data to the appropriate processor
    //since positions, velocities, masses are all at different points in the file,
    //to correctly assign particle to proccessor with correct velocities and mass must have several file pointers
    Particle *Pbuf;
    Int_t chunksize=GADGETCHUNKSIZE,nchunk;
    FLOAT *ctempchunk, *vtempchunk,*sphchunk;
    REAL *dtempchunk;
    GADGETIDTYPE *idvalchunk;
/*
#ifdef USEMPI
    MPI_Status status;
    if (ThisTask==0) Pbuf=new Particle[BufSize*NProcs];
    Int_t Nlocalbuf,ibuf=0,*Nbuf;
    if (ThisTask==0) {
        Nbuf=new Int_t[NProcs];
        for (int j=0;j<NProcs;j++) Nbuf[j]=0;
    }
    if (ThisTask==0) {
#endif
*/
    ctempchunk=new FLOAT[3*chunksize];
    vtempchunk=new FLOAT[3*chunksize];
    dtempchunk=new REAL[chunksize];
    sphchunk=new FLOAT[chunksize*opt.sphnio];//number of sph io data in the file
    idvalchunk=new GADGETIDTYPE[chunksize];
/*
#ifdef USEMPI
    }
#endif
*/
    Nlocal=0;
/*
#ifdef USEMPI
#ifndef MPIREDUCEMEM
    MPIDomainExtentGadget(opt);
    if (NProcs>1) {
    MPIDomainDecompositionGadget(opt);
    MPIInitialDomainDecomposition();
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    //if (ThisTask==0) {
#endif
*/
    //opening file

    Fgad=new fstream[opt.nfiles];
    Fgadvel=new fstream[opt.nfiles];
    Fgadid=new fstream[opt.nfiles];
    Fgadmass=new fstream[opt.nfiles];
    Fgadsphdata=new fstream[opt.nfiles*opt.sphnio];
    header=new gadget_header[opt.nfiles];
    //load the header and offset all the files to the appropriate place
    for(i=0; i<opt.nfiles; i++)
    {
        if(opt.nfiles>1) sprintf(buf,"%s.%d",opt.fname,i);
        else sprintf(buf,"%s",opt.fname);
        Fgad[i].open(buf,ios::in);
        Fgadvel[i].open(buf,ios::in);
        Fgadid[i].open(buf,ios::in);
        Fgadmass[i].open(buf,ios::in);
        for (int sphblocks=0;sphblocks<opt.sphnio;sphblocks++) Fgadsphdata[i+sphblocks].open(buf,ios::in);
        if(!Fgad[i]) {
            cout<<"can't open file "<<buf<<endl;
            exit(0);
        }
        else cout<<"reading "<<buf<<endl;
#ifdef GADGET2FORMAT
        Fgad[i].read((char*)&dummy, sizeof(dummy));
        Fgad[i].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
        Fgad[i].read((char*)&dummy, sizeof(dummy));
        Fgad[i].read((char*)&dummy, sizeof(dummy));
        fprintf(stderr,"reading... %s\n",DATA);
#endif
        Fgad[i].read((char*)&dummy, sizeof(dummy));
        Fgad[i].read((char*)&header[i], sizeof(gadget_header));
        Fgad[i].read((char*)&dummy, sizeof(dummy));
        //endian indep call
        header[i].Endian();
        count=0;for(k=0;k<NGTYPE;k++)count+=header[i].npart[k];
        countsph=header[i].npart[0];

#ifdef GADGET2FORMAT
        Fgadvel[i].read((char*)&dummy, sizeof(dummy));
        Fgadvel[i].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
        Fgadvel[i].read((char*)&dummy, sizeof(dummy));
        Fgadvel[i].read((char*)&dummy, sizeof(dummy));

        Fgadid[i].read((char*)&dummy, sizeof(dummy));
        Fgadid[i].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
        Fgadid[i].read((char*)&dummy, sizeof(dummy));
        Fgadid[i].read((char*)&dummy, sizeof(dummy));

        Fgadmass[i].read((char*)&dummy, sizeof(dummy));
        Fgadmass[i].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));

        for (int sphblocks=0;sphblocks<opt.sphnio;sphblocks++){
            Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
            Fgadsphdata[i+sphblocks].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
            Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
            Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
        }
#endif
        Fgadvel[i].read((char*)&dummy, sizeof(dummy));
        Fgadvel[i].read((char*)&header[i], sizeof(gadget_header));
        Fgadvel[i].read((char*)&dummy, sizeof(dummy));
        Fgadid[i].read((char*)&dummy, sizeof(dummy));
        Fgadid[i].read((char*)&header[i], sizeof(gadget_header));
        Fgadid[i].read((char*)&dummy, sizeof(dummy));
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));
        Fgadmass[i].read((char*)&header[i], sizeof(gadget_header));
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));
        for (int sphblocks=0;sphblocks<opt.sphnio;sphblocks++){
            Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
            Fgadsphdata[i+sphblocks].read((char*)&header[i], sizeof(gadget_header));
            Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
        }
        //endian indep call
        header[i].Endian();

        //move to velocities
#ifdef GADGET2FORMAT
        Fgadvel[i].read((char*)&dummy, sizeof(dummy));
        Fgadvel[i].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
        Fgadvel[i].read((char*)&dummy, sizeof(dummy));
        Fgadvel[i].read((char*)&dummy, sizeof(dummy));
#endif
        Fgadvel[i].read((char*)&dummy, sizeof(dummy));
        Fgadvel[i].seekp(count*3*sizeof(FLOAT),ios::cur);
        Fgadvel[i].read((char*)&dummy, sizeof(dummy));

        //move to ids
#ifdef GADGET2FORMAT
        Fgadid[i].read((char*)&dummy, sizeof(dummy));
        Fgadid[i].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
        Fgadid[i].read((char*)&dummy, sizeof(dummy));
        Fgadid[i].read((char*)&dummy, sizeof(dummy));
#endif
        Fgadid[i].read((char*)&dummy, sizeof(dummy));
        Fgadid[i].seekp(count*3*sizeof(FLOAT),ios::cur);
        Fgadid[i].read((char*)&dummy, sizeof(dummy));
#ifdef GADGET2FORMAT
        Fgadid[i].read((char*)&dummy, sizeof(dummy));
        Fgadid[i].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
        Fgadid[i].read((char*)&dummy, sizeof(dummy));
        Fgadid[i].read((char*)&dummy, sizeof(dummy));
#endif
        Fgadid[i].read((char*)&dummy, sizeof(dummy));
        Fgadid[i].seekp(count*3*sizeof(FLOAT),ios::cur);
        Fgadid[i].read((char*)&dummy, sizeof(dummy));

        //move to masses
#ifdef GADGET2FORMAT
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));
        Fgadmass[i].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));
#endif
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));
        Fgadmass[i].seekp(count*3*sizeof(FLOAT),ios::cur);
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));
#ifdef GADGET2FORMAT
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));
        Fgadmass[i].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));
#endif
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));
        Fgadmass[i].seekp(count*3*sizeof(FLOAT),ios::cur);
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));

#ifdef GADGET2FORMAT
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));
        Fgadmass[i].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));
#endif
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));
        Fgadmass[i].seekp(count*sizeof(idval),ios::cur);
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));

        //move to sph data
        for (int sphblocks=0;sphblocks<opt.sphnio;sphblocks++){
#ifdef GADGET2FORMAT
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
        Fgadsphdata[i+sphblocks].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
#endif
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
        Fgadsphdata[i+sphblocks].seekp(count*3*sizeof(FLOAT),ios::cur);
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
#ifdef GADGET2FORMAT
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
        Fgadsphdata[i+sphblocks].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
#endif
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
        Fgadsphdata[i+sphblocks].seekp(count*3*sizeof(FLOAT),ios::cur);
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));

#ifdef GADGET2FORMAT
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
        Fgadsphdata[i+sphblocks].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
#endif
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
        Fgadsphdata[i+sphblocks].seekp(count*sizeof(idval),ios::cur);
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));

#ifdef GADGET2FORMAT
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
        Fgadsphdata[i+sphblocks].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
#endif
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
        for (int k=0;k<NGTYPE;k++)
            if(header[i].mass[k]==0)
              Fgadsphdata[i+sphblocks].seekp(header[i].npart[k]*sizeof(REAL),ios::cur);
        Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
        }

#ifdef GADGET2FORMAT
        Fgad[i].read((char*)&dummy, sizeof(dummy));
        Fgad[i].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
        Fgad[i].read((char*)&dummy, sizeof(dummy));
        Fgad[i].read((char*)&dummy, sizeof(dummy));

        Fgadvel[i].read((char*)&dummy, sizeof(dummy));
        Fgadvel[i].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
        Fgadvel[i].read((char*)&dummy, sizeof(dummy));
        Fgadvel[i].read((char*)&dummy, sizeof(dummy));
        Fgadvel[i].read((char*)&dummy, sizeof(dummy));

        Fgadid[i].read((char*)&dummy, sizeof(dummy));
        Fgadid[i].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
        Fgadid[i].read((char*)&dummy, sizeof(dummy));
        Fgadid[i].read((char*)&dummy, sizeof(dummy));

        Fgadmass[i].read((char*)&dummy, sizeof(dummy));
        Fgadmass[i].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));
        Fgadmass[i].read((char*)&dummy, sizeof(dummy));

        for (int sphblocks=0;sphblocks<opt.sphnio;sphblocks++){
            Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
            Fgadsphdata[i+sphblocks].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
            Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
            Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
        }
#endif

        for (int sphblocks=1;sphblocks<opt.sphnio;sphblocks++){
            Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
            Fgadsphdata[i+sphblocks].seekp((sphblocks*countsph)*sizeof(FLOAT),ios::cur);
            Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
#ifdef GADGET2FORMAT
            Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
            Fgadsphdata[i+sphblocks].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
            Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
            Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
#endif
        }
    }
    //store useful cosmological quantities and header information
    opt.p=header[ifirstfile].BoxSize;
    z=header[ifirstfile].redshift;
    opt.a=1./(1.+z);
    opt.Omega_m=header[ifirstfile].Omega0;
    opt.Omega_Lambda=header[ifirstfile].OmegaLambda;
    opt.h=header[ifirstfile].HubbleParam;
    //Hubble flow
    if (opt.comove) aadjust=1.0;
    else aadjust=opt.a;
    Hubble=opt.h*opt.H*sqrt((1-opt.Omega_m-opt.Omega_Lambda)*pow(aadjust,-2.0)+opt.Omega_m*pow(aadjust,-3.0)+opt.Omega_Lambda);
    opt.rhobg=3.*Hubble*Hubble/(8.0*M_PI*opt.G)*opt.Omega_m;
    HubbleFlow=0.;
    //based on Bryan and Norman 1998 the virialization level is given by
    Double_t bnx=-((1-opt.Omega_m-opt.Omega_Lambda)*pow(aadjust,-2.0)+opt.Omega_Lambda)/((1-opt.Omega_m-opt.Omega_Lambda)*pow(aadjust,-2.0)+opt.Omega_m*pow(aadjust,-3.0)+opt.Omega_Lambda);
    opt.virlevel=(18.0*M_PI*M_PI+82.0*bnx-39*bnx*bnx)/opt.Omega_m;
    mscale=opt.M/opt.h;lscale=opt.L/opt.h*aadjust;lvscale=opt.L/opt.h*opt.a;
    opt.uinfo.eps*=opt.p*lscale;
    Ntot=0;
    for (j=0;j<6;j++) {
      //opt.numpart[j]=header[ifirstfile].npartTotal[j];
      Ntot+=header[ifirstfile].npartTotal[j];
      //adjust for HighWord part of gadget header
      //opt.numpart[j]+=((long long)header[ifirstfile].npartTotalHW[j]<<32);
      Ntot+=((long long)header[ifirstfile].npartTotalHW[j]<<32);
    }
    //if no mass is stored then assume mass equal and stored in first nonzero mass value in header.mass
#ifdef NOMASS
    for(k=0;k<6;k++) if(header[ifirstfile].mass[k]!=0) {opt.MassValue=header[ifirstfile].mass[k]*mscale;break;}
#endif
    cout<<"Snapshot contains "<<Ntot<<" particles at is at time "<<opt.a<<endl;
    cout<<"Particle system contains "<<nbodies<<" particles at is at time "<<opt.a<<endl;
    ///\todo will have to make all of this mpi compatible at some point.
    ///probably means that i need to carefully map the particle ids to local index
    ///when determining where in the local pfof array a particle should be placed
    ///and how the halos should be split.
#ifndef USEMPI
    Nlocal=Ntot;
#endif
    //build array to store haloids with mapping of ids to index
    if (opt.Neff==0) opt.Neff=Nlocal;
    //pfof=new Int_t[opt.Neff];
    pfof=new Int_t[opt.Neff];
    cout<<"Size allocated to store particle index data is "<<opt.Neff*sizeof(Int_t)/1024./1024./1024.<<endl;
    ncount=new Int_t[ngroups+1];
    ncount[0]=0;
    for (i=0;i<opt.Neff;i++) pfof[i]=0;

    cout<<"setting pfof for "<<ngroups<<endl;
    for (i=0;i<ngroups;i++) {
        ncount[i+1]=0;
        for (j=0;j<hp[i].NumberofParticles;j++) {
            pfof[MapPIDStoIndex(opt,hp[i].ParticleID[j])]=i+1;
        }
    }
    cout<<"done"<<endl;

    //now load the actual particle data
    for(i=0,count=0,pc=0;i<opt.nfiles; i++,pc=pc_new,count=count2)
    {
        for(k=0, ntot_withmasses=0,Ntotfile=0; k<6; k++) {if(header[i].mass[k]==0) ntot_withmasses+= header[i].npart[k];Ntotfile+=header[i].npart[k];}
        //now read positions, velocities, masses, etc
        Fgad[i].read((char*)&dummy, sizeof(dummy));
        if (dummy/Ntotfile/3!=sizeof(FLOAT)) {cout<<" mismatch in position type size, file has "<<dummy/Ntotfile/3<<" but using "<<sizeof(FLOAT)<<endl;exit(9);}
        Fgadvel[i].read((char*)&dummy, sizeof(dummy));
        Fgadid[i].read((char*)&dummy, sizeof(dummy));
        if (dummy/Ntotfile!=sizeof(GADGETIDTYPE)) {cout<<" mismatch in id type size, file has "<<dummy/Ntotfile<<" but using "<<sizeof(GADGETIDTYPE)<<endl;exit(9);}
        if (ntot_withmasses>0) Fgadmass[i].read((char*)&dummy, sizeof(dummy));
        if (ntot_withmasses>0) if (dummy/ntot_withmasses!=sizeof(REAL)) {cout<<" mismatch in mass type size, file has "<<dummy/ntot_withmasses<<" but using "<<sizeof(REAL)<<endl;exit(9);}
        for (int sphblocks=0;sphblocks<opt.sphnio;sphblocks++){
            Fgadsphdata[i+sphblocks].read((char*)&dummy, sizeof(dummy));
            if (dummy/header[i].npart[0]!=sizeof(FLOAT)) {cout<<" mismatch in sph quantity "<<sphblocks<<" type size, file has "<<(float)dummy/(float)header[i].npart[0]<<" but using "<<sizeof(FLOAT)<<endl;exit(9);}
        }
        int bcount=0,blahcount=0;
        int blahnbodies=0;
        for(k=1;k<4;k++)blahnbodies+=header[i].npart[k];
        cout<<blahnbodies<<endl;
        for(k=0,count2=count,pc_new=pc;k<6;k++)if (header[i].npart[k]>0)
        {
            //read gadget data in chunks
            if (header[i].npart[k]<chunksize)nchunk=header[i].npart[k];
            else nchunk=chunksize;
            for(n=0;n<header[i].npart[k];n+=nchunk)
            {

                if (header[i].npart[k]-n<chunksize&&header[i].npart[k]-n>0)nchunk=header[i].npart[k]-n;
                Fgad[i].read((char*)ctempchunk, sizeof(FLOAT)*3*nchunk);
                Fgadvel[i].read((char*)vtempchunk, sizeof(FLOAT)*3*nchunk);
                Fgadid[i].read((char*)idvalchunk, sizeof(GADGETIDTYPE)*nchunk);
#ifndef NOMASS
                if(header[i].mass[k]==0) Fgadmass[i].read((char*)dtempchunk, sizeof(REAL)*nchunk);
#endif
                ///sph data contains several quantities of interest
                if (k==0) {
                    for (int sphblocks=0;sphblocks<opt.sphnio;sphblocks++) {
                        Fgadsphdata[i+sphblocks].read((char*)&sphchunk[sphblocks*nchunk], sizeof(FLOAT)*nchunk);
                    }
                }
                for (int nn=0;nn<nchunk;nn++) {
                ctemp[0]=ctempchunk[0+3*nn];ctemp[1]=ctempchunk[1+3*nn];ctemp[2]=ctempchunk[2+3*nn];
                vtemp[0]=vtempchunk[0+3*nn];vtemp[1]=vtempchunk[1+3*nn];vtemp[2]=vtempchunk[2+3*nn];
                idval=idvalchunk[nn];
#ifndef NOMASS
                if(header[i].mass[k]==0) {
                    dtemp=dtempchunk[nn];
                }
                else dtemp=header[i].mass[k];
#endif
                for (m=0;m<3;m++) {ctemp[m]=LittleFLOAT(ctemp[m]);vtemp[m]=LittleFLOAT(vtemp[m]);}
#ifdef GADGETLONGID
                idval=LittleLongInt(idval);
#else
                idval=LittleInt(idval);
#endif
#ifndef NOMASS
                if(header[i].mass[k]==0) dtemp=LittleREAL(dtemp);
                if(dtemp<MP_DM&&dtemp>0) MP_DM=dtemp;
#else
                dtemp=1.0;
#endif
                indexval=MapPIDStoIndex(opt, idval);

                if (k==0 || k==4 || k==5) {indexval=bcount+blahnbodies;bcount++;}
                else {indexval=blahcount;blahcount++;}
                if (pfof[indexval]>0) {
                    groupval=pfof[indexval];
                    indexval=hp[groupval-1].noffset+ncount[groupval];
                    Part[indexval]=Particle(dtemp*mscale,
                    ctemp[0]*lscale,ctemp[1]*lscale,ctemp[2]*lscale,
                    vtemp[0]*opt.V*sqrt(opt.a)+HubbleFlow*ctemp[0]*lvscale,
                    vtemp[1]*opt.V*sqrt(opt.a)+HubbleFlow*ctemp[1]*lvscale,
                    vtemp[2]*opt.V*sqrt(opt.a)+HubbleFlow*ctemp[2]*lvscale);
                    if (k==GGASTYPE) {
                        opt.npart[GASTYPE]++;
                        Part[indexval].SetType(GASTYPE);
                    }
                    else if (k==GSTARTYPE) {
                        opt.npart[STARTYPE]++;
                        Part[indexval].SetType(STARTYPE);
                    }
                    else if (k==GBHTYPE) {
                        opt.npart[BHTYPE]++;
                        Part[indexval].SetType(BHTYPE);
                    }
                    else {
                        opt.npart[DMTYPE]++;
                        Part[indexval].SetType(DMTYPE);
                    }

                    Part[indexval].SetPID(idval);
                    Part[indexval].SetID(count2);
#ifdef GASON
                    if (k==GGASTYPE) {
                        //set self energy
                        sphtemp=sphchunk[nn+sphioblock_u*nchunk];
                        Part[indexval].SetU(sphtemp);
                        //set density
                        sphtemp=sphchunk[nn+sphioblock_d*nchunk];
                        Part[indexval].SetDensity(sphtemp);
                    }
#endif
                    ncount[groupval]++;
                    if (k==GGASTYPE) hp[groupval-1].NumofType[GASTYPE]++;
                    else if (k==GSTARTYPE) hp[groupval-1].NumofType[STARTYPE]++;
                    else if (k==GBHTYPE) hp[groupval-1].NumofType[BHTYPE]++;
                    else hp[groupval-1].NumofType[DMTYPE]++;

                }
                /*
                ibuf=MPIGetParticlesProcessor(ctemp[0],ctemp[1],ctemp[2]);
                Pbuf[ibuf*BufSize+Nbuf[ibuf]]=Particle(dtemp*mscale,
                    ctemp[0]*lscale,ctemp[1]*lscale,ctemp[2]*lscale,
                    vtemp[0]*opt.V*sqrt(opt.a)+Hubble*ctemp[0]*lvscale,
                    vtemp[1]*opt.V*sqrt(opt.a)+Hubble*ctemp[1]*lvscale,
                    vtemp[2]*opt.V*sqrt(opt.a)+Hubble*ctemp[2]*lvscale,
                    count2,k);
                //Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPID(count2);
                Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPID(idval);
                Nbuf[ibuf]++;
                if(ibuf==0){
                    Nbuf[ibuf]--;
                    Part[Nlocal++]=Pbuf[ibuf*BufSize+Nbuf[ibuf]];
                }
                else {
                    if(Nbuf[ibuf]==BufSize) {
                        MPI_Ssend(&Nbuf[ibuf], 1, MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
                        MPI_Ssend(&Pbuf[ibuf*BufSize],sizeof(Particle)*Nbuf[ibuf],MPI_BYTE,ibuf,ibuf,MPI_COMM_WORLD);
                        Nbuf[ibuf]=0;
                    }
                }*/
                count2++;
            }
        }
        }
        //more information contained in sph particles and if there is sf feed back but for the moment, ignore
        Fgad[i].close();
        Fgadvel[i].close();
        Fgadid[i].close();
        Fgadmass[i].close();
        for (int sphblocks=0;sphblocks<opt.sphnio;sphblocks++) Fgadsphdata[i+sphblocks].close();
    }
    for (i=0;i<ngroups;i++) {
        for (j=0;j<NPARTTYPES;j++)
             hp[i].AllNumofType[j]=hp[i].NumofType[j];
    }

    /*
#ifdef USEMPI
    //once finished reading the file if there are any particles left in the buffer broadcast them
    for(ibuf = 1; ibuf < NProcs; ibuf++)
    {
        while(Nbuf[ibuf])
        {
            MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
            MPI_Ssend(&Pbuf[ibuf*BufSize], sizeof(Particle)*Nbuf[ibuf], MPI_BYTE, ibuf, ibuf, MPI_COMM_WORLD);
            Nbuf[ibuf]=0;
        }
        //last broadcast with Nbuf[ibuf]=0 so that receiver knows no more particles are to be broadcast
        MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t,ibuf,ibuf+NProcs,MPI_COMM_WORLD);
    }
#endif
*/
    //}//end of reading task section
/*    
#ifdef USEMPI
    else {
        do {
            MPI_Recv(&Nlocalbuf, 1, MPI_Int_t, 0, ThisTask+NProcs, MPI_COMM_WORLD, &status);
            if(Nlocalbuf) {
                MPI_Recv(&Part[Nlocal],sizeof(Particle)*Nlocalbuf,MPI_BYTE,0,ThisTask, MPI_COMM_WORLD,&status);
                Nlocal+=Nlocalbuf;
            }
        } while(Nlocalbuf);
    }
#endif
    */
    //sort particles according to type
#ifdef USEOPENMP
#pragma omp parallel default(shared)  
{
    #pragma omp for
#endif
    for (i=0;i<ngroups;i++)
    {
        gsl_heapsort(&Part[hp[i].noffset],hp[i].NumberofParticles,sizeof(Particle),TypeCompare);
    }
#ifdef USEOPENMP
}
#endif

#ifdef USEMPI
    if (ThisTask==0) {
#endif
    //adjust period
    opt.p*=lscale;
#ifdef USEMPI
    }
#endif
#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
    //update cosmological data and boundary in code units
    MPI_Bcast(&(opt.p),sizeof(opt.p),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.a),sizeof(opt.a),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.Omega_m),sizeof(opt.Omega_m),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.Omega_Lambda),sizeof(opt.Omega_Lambda),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.h),sizeof(opt.h),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.rhobg),sizeof(opt.rhobg),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.virlevel),sizeof(opt.virlevel),MPI_BYTE,0,MPI_COMM_WORLD);
#ifdef NOMASS
    MPI_Bcast(&(opt.MassValue),sizeof(opt.MassValue),MPI_BYTE,0,MPI_COMM_WORLD);
#endif
    MPI_Bcast(&(Ntotal),sizeof(Ntotal),MPI_BYTE,0,MPI_COMM_WORLD);
#endif
    cout<<"Finished loading gadget"<<endl;
}

//@}

