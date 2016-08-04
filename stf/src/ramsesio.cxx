/*! \file ramsesio.cxx
 *  \brief this file contains routines for ramses snapshot file io
 * 
 * \todo need to check if amr file quantity ngrid_current is actually the number of cells in the file as 
 * an example fortran code I have been given seems to also use the ngridlevel array, which stores the number of cells 
 * at a given resolution level. 
 */

//-- RAMSES SPECIFIC IO

#include "stf.h"

#include "ramsesitems.h"
#include "endianutils.h"


int RAMSES_fortran_read(fstream &F, int &i){
    int dummy,byteoffset=0;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    F.read((char*)&i,sizeof(int)); byteoffset+=dummy;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    return byteoffset;
}
int RAMSES_fortran_read(fstream &F, RAMSESFLOAT &f){
    int dummy,byteoffset=0;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    F.read((char*)&f,sizeof(RAMSESFLOAT)); byteoffset+=dummy;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    return byteoffset;
}
int RAMSES_fortran_read(fstream &F, int *i){
    int dummy,byteoffset=0;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    F.read((char*)i,dummy); byteoffset+=dummy;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    return byteoffset;
}
int RAMSES_fortran_read(fstream &F, long long *i){
    int dummy,byteoffset=0;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    F.read((char*)i,dummy); byteoffset+=dummy;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    return byteoffset;
}
int RAMSES_fortran_read(fstream &F, RAMSESFLOAT *f){
    int dummy,byteoffset=0;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    F.read((char*)f,dummy); byteoffset+=dummy;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    return byteoffset;
}

int RAMSES_fortran_skip(fstream &F, int nskips){
    int dummy,byteoffset=0;
    for (int i=0;i<nskips;i++) {
        F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
        F.seekg(dummy,ios::cur); byteoffset+=dummy;
        F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    }
    return byteoffset;
}

Int_t RAMSES_get_nbodies(char *fname, int ptype, Options &opt)
{
    char buf[2000],buf1[2000],buf2[2000];
    string stringbuf;
    sprintf(buf1,"amr_%s.out00001",fname);
    sprintf(buf2,"amr_%s.out",fname);
    if (FileExists(buf1)) sprintf(buf,"%s",buf1);
    else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
    else {
        printf("Error. Can't find AMR data \nneither as `%s'\nnor as `%s'\n\n", buf1, buf2);
        exit(9);
    }
    //if gas searched in some fashion then load amr/hydro data
    if (opt.partsearchtype==PSTGAS||opt.partsearchtype==PSTALL||(opt.partsearchtype==PSTDARK&&opt.iBaryonSearch)) {
    sprintf(buf1,"hydro_%s.out00001",fname);
    sprintf(buf2,"hydro_%s.out",fname);
    if (FileExists(buf1)) sprintf(buf,"%s",buf1);
    else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
    else {
        printf("Error. Can't find Hydro data \nneither as `%s'\nnor as `%s'\n\n", buf1, buf2);
        exit(9);
    }
    }
    sprintf(buf1,"part_%s.out00001",fname);
    sprintf(buf2,"part_%s.out",fname);
    if (FileExists(buf1)) sprintf(buf,"%s",buf1);
    else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
    else {
        printf("Error. Can't find Particle data \nneither as `%s'\nnor as `%s'\n\n", buf1, buf2);
        exit(9);
    }

    fstream Framses;
    RAMSES_Header ramses_header_info;
    //buffers to load data
    int intbuff[NRAMSESTYPE];
    long long longbuff[NRAMSESTYPE];
    int i,j,k,ireaderror=0;
    Int_t nbodies=0;
    //IntType inttype;
    int dummy;

    int nusetypes,usetypes[NRAMSESTYPE];

    if (ptype==PSTALL) {nusetypes=4;usetypes[0]=RAMSESGASTYPE;usetypes[1]=RAMSESDMTYPE;usetypes[2]=RAMSESSTARTYPE;usetypes[3]=RAMSESBHTYPE;}
    else if (ptype==PSTDARK) {nusetypes=1;usetypes[0]=RAMSESDMTYPE;}
    else if (ptype==PSTGAS) {nusetypes=1;usetypes[0]=RAMSESGASTYPE;}
    else if (ptype==PSTSTAR) {nusetypes=1;usetypes[0]=RAMSESSTARTYPE;}
    else if (ptype==PSTBH) {nusetypes=1;usetypes[0]=RAMSESBHTYPE;}


    //Open the specified file and the specified dataset in the file.
    //first open amr data
    sprintf(buf1,"amr_%s.out00001",fname);
    sprintf(buf2,"amr_%s.out",fname);
    if (FileExists(buf1)) sprintf(buf,"%s",buf1);
    else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
    Framses.open(buf, ios::binary|ios::in);
    //read header info 
    //this is based on a ramsestotipsy fortran code that does not detail everything in the header block but at least its informative. The example has 
    /*
    read(10)ncpu
    read(10)ndim
    read(10)nx,ny,nz
    read(10)nlevelmax
    read(10)ngridmax
    read(10)nboundary
    read(10)ngrid_current
    read(10)boxlen
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)msph
    close(10)
    */
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.num_files, sizeof(int));
    Framses.read((char*)&dummy, sizeof(dummy));
    opt.num_files=ramses_header_info.num_files;

    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.ndim, sizeof(int));
    Framses.read((char*)&dummy, sizeof(dummy));

    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.nx, sizeof(int));
    Framses.read((char*)&ramses_header_info.ny, sizeof(int));
    Framses.read((char*)&ramses_header_info.nz, sizeof(int));
    Framses.read((char*)&dummy, sizeof(dummy));

    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.nlevelmax, sizeof(int));
    Framses.read((char*)&dummy, sizeof(dummy));

    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.ngridmax, sizeof(int));
    Framses.read((char*)&dummy, sizeof(dummy));

    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.nboundary, sizeof(int));
    Framses.read((char*)&dummy, sizeof(dummy));

    //this would be number of active grids but ignore for now
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.seekg(dummy,ios::cur);
    Framses.read((char*)&dummy, sizeof(dummy));
    
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.BoxSize, sizeof(RAMSESFLOAT));
    Framses.read((char*)&dummy, sizeof(dummy));

    //now skip 10 blocks
    for (i=0;i<10;i++) {
        Framses.read((char*)&dummy, sizeof(dummy));
        //skip block size given by dummy
        Framses.seekg(dummy,ios::cur);
        Framses.read((char*)&dummy, sizeof(dummy));
    }
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.mass[RAMSESGASTYPE], sizeof(RAMSESFLOAT));
    Framses.read((char*)&dummy, sizeof(dummy));
    
    Framses.close();

    //reopen to get number of amr cells might need to alter to read grid information and what cells have no so-called son cells
    if (opt.partsearchtype==PSTGAS||opt.partsearchtype==PSTALL||(opt.partsearchtype==PSTDARK&&opt.iBaryonSearch)) {
    for (i=0;i<ramses_header_info.num_files;i++) {
        sprintf(buf1,"amr_%s.out%05d",fname,i+1);
        sprintf(buf2,"amr_%s.out",fname);
        if (FileExists(buf1)) sprintf(buf,"%s",buf1);
        else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
        Framses.open(buf, ios::binary|ios::in);
        for (j=0;j<6;j++) {
            Framses.read((char*)&dummy, sizeof(dummy));
            Framses.seekg(dummy,ios::cur);
            Framses.read((char*)&dummy, sizeof(dummy));
        }
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.read((char*)&ramses_header_info.npart[RAMSESGASTYPE], sizeof(int));
        Framses.read((char*)&dummy, sizeof(dummy));
        ramses_header_info.npartTotal[RAMSESGASTYPE]+=ramses_header_info.npart[RAMSESGASTYPE];
    }

    //now hydro header data 
    sprintf(buf1,"hydro_%s.out00001",fname);
    sprintf(buf2,"hydro_%s.out",fname);
    if (FileExists(buf1)) sprintf(buf,"%s",buf1);
    else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
    Framses.open(buf, ios::binary|ios::in);

    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.seekg(dummy,ios::cur);
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.nvarh, sizeof(int));
    Framses.read((char*)&dummy, sizeof(dummy));
    for (i=0;i<3;i++) {
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.seekg(dummy,ios::cur);
        Framses.read((char*)&dummy, sizeof(dummy));
    }

    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.gamma_index, sizeof(RAMSESFLOAT));
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.close();
    }

    //now cosmological header data and unit info
    /*
    sprintf(buf1,"info_%s.txt",fname);
    Framses.open(buf1, ios::in);
    getline(Framses,stringbuf);
    getline(Framses,stringbuf);
    Framses>>stringbuf>>ramses_header_info.levelmin;
    getline(Framses,stringbuf);
    getline(Framses,stringbuf);
    getline(Framses,stringbuf);
    getline(Framses,stringbuf);
    getline(Framses,stringbuf);
    Framses>>stringbuf>>ramses_header_info.time;
    Framses>>stringbuf>>ramses_header_info.aexp;
    Framses>>stringbuf>>ramses_header_info.HubbleParam;
    Framses>>stringbuf>>ramses_header_info.Omegam;
    Framses>>stringbuf>>ramses_header_info.OmegaLambda;
    Framses>>stringbuf>>ramses_header_info.Omegak;
    Framses>>stringbuf>>ramses_header_info.Omegab;
    Framses>>stringbuf>>ramses_header_info.scale_l;
    Framses>>stringbuf>>ramses_header_info.scale_d;
    Framses>>stringbuf>>ramses_header_info.scale_t;
    */
    //fortran code has after this something to do with the order of something 
    /*
  read(10,*)
  read(10,'("ordering type=",A80)'),ordering
  read(10,*)
  allocate(cpu_list(1:ncpu))
  if(TRIM(ordering).eq.'hilbert')then
     allocate(bound_key(0:ncpu))
     allocate(cpu_read(1:ncpu))
     cpu_read=.false.
     do impi=1,ncpu
        read(10,'(I8,1X,E23.15,1X,E23.15)')i,bound_key(impi-1),bound_key(impi)
     end do
  endif
  */

    //now particle info
    for (i=0;i<ramses_header_info.num_files;i++) {
        sprintf(buf1,"part_%s.out%05d",fname,i+1);
        sprintf(buf2,"part_%s.out",fname);
        if (FileExists(buf1)) sprintf(buf,"%s",buf1);
        else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
        Framses.open(buf, ios::binary|ios::in);

        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.seekg(dummy,ios::cur);
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.seekg(dummy,ios::cur);
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.read((char*)&ramses_header_info.npart[RAMSESDMTYPE], sizeof(int));
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.seekg(dummy,ios::cur);
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.read((char*)&dummy, sizeof(dummy));
        //?? had total but probably not correct
        //Framses.read((char*)&ramses_header_info.npartTotal[RAMSESSTARTYPE], sizeof(int));
        Framses.read((char*)&ramses_header_info.npart[RAMSESSTARTYPE], sizeof(int));
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.seekg(dummy,ios::cur);
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.seekg(dummy,ios::cur);
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.read((char*)&ramses_header_info.npart[RAMSESSINKTYPE], sizeof(int));
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.close();
        ramses_header_info.npartTotal[RAMSESDMTYPE]+=ramses_header_info.npart[RAMSESDMTYPE];
        ramses_header_info.npartTotal[RAMSESSTARTYPE]+=ramses_header_info.npart[RAMSESSTARTYPE];
        ramses_header_info.npartTotal[RAMSESSINKTYPE]+=ramses_header_info.npart[RAMSESSINKTYPE];
        ramses_header_info.npartTotal[RAMSESDMTYPE]-=ramses_header_info.npart[RAMSESSTARTYPE];
        ramses_header_info.npartTotal[RAMSESDMTYPE]-=ramses_header_info.npart[RAMSESSINKTYPE];
    }
    ramses_header_info.npartTotal[RAMSESDMTYPE]-=ramses_header_info.npartTotal[RAMSESSTARTYPE];
    for(j=0, nbodies=0; j<nusetypes; j++) {
        k=usetypes[j];
        nbodies+=ramses_header_info.npartTotal[k];
        //nbodies+=((long long)(ramses_header_info.npartTotalHW[k]) << 32);
    }
    for (j=0;j<NPARTTYPES;j++) opt.numpart[j]=0;
    if (ptype==PSTALL || ptype==PSTDARK) opt.numpart[DARKTYPE]=ramses_header_info.npartTotal[RAMSESDMTYPE];
    if (ptype==PSTALL || ptype==PSTGAS) opt.numpart[GASTYPE]=ramses_header_info.npartTotal[RAMSESGASTYPE];
    if (ptype==PSTALL || ptype==PSTSTAR) opt.numpart[STARTYPE]=ramses_header_info.npartTotal[RAMSESSTARTYPE];
    if (ptype==PSTALL || ptype==PSTBH) opt.numpart[BHTYPE]=ramses_header_info.npartTotal[RAMSESSINKTYPE];
    return nbodies;

}

///reads a ramses file. If cosmological simulation uses cosmology (generally assuming LCDM or small deviations from this) to estimate the mean interparticle spacing
///and scales physical linking length passed by this distance. Also reads header and overrides passed cosmological parameters with ones stored in header.
///\todo need to implement the multiple reading threads send/receive to each other. An example of multiple read threads is seen in \ref gadgetio.cxx. It invovles
///having read threads read files over again and store all particles that need to be sent to other read threads, then sending this info
void ReadRamses(Options &opt, Particle *&Part, const Int_t nbodies,Particle *&Pbaryons, Int_t nbaryons)
{
    char buf[2000],buf1[2000],buf2[2000];
    string stringbuf,orderingstring;
    fstream Finfo;
    fstream *Famr;
    fstream *Fhydro;
    fstream *Fpart, *Fpartvel,*Fpartid,*Fpartmass, *Fpartlevel, *Fpartage, *Fpartmet;
    RAMSES_Header *header;
    int intbuff[NRAMSESTYPE];
    long long longbuff[NRAMSESTYPE];
    int i,j,k,n,idim,ivar,igrid,ireaderror=0;
    Int_t count2,bcount2;
    //IntType inttype;
    int dummy,byteoffset;
    Double_t MP_DM=MAXVALUE,LN,N_DM,MP_B=0;
    double z,aadjust,Hubble,Hubbleflow;
    Double_t mscale,lscale,lvscale,rhoscale;
    Double_t mtemp,utemp,rhotemp,Ztemp,Ttemp;
    Coordinate xpos,vpos;
    RAMSESFLOAT xtemp[3],vtemp[3];
    RAMSESIDTYPE idval;
    int typeval;
    RAMSESFLOAT ageval,metval;
    int *ngridlevel,*ngridbound,*ngridfile;
    int lmin=1000000,lmax=0;

    int ifirstfile=0,*ireadfile,ibuf=0;
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
    ireadfile=new int[opt.num_files];
    for (i=0;i<opt.num_files;i++) ireadfile[i]=1;
#endif
    ///\todo because of the stupid fortran format, easier if chunksize is BIG so that 
    ///number of particles local to a file are smaller
    Int_t chunksize=RAMSESCHUNKSIZE,nchunk;
    RAMSESFLOAT *xtempchunk, *vtempchunk, *mtempchunk, *sphtempchunk, *agetempchunk, *mettempchunk, *hydrotempchunk;
    RAMSESIDTYPE *idvalchunk;
    int *icellchunk;

    Famr=new fstream[opt.num_files];
    Fhydro=new fstream[opt.num_files];
    Fpart=new fstream[opt.num_files];
    Fpartvel=new fstream[opt.num_files];
    Fpartmass=new fstream[opt.num_files];
    Fpartid=new fstream[opt.num_files];
    Fpartlevel=new fstream[opt.num_files];
    Fpartage=new fstream[opt.num_files];
    Fpartmet=new fstream[opt.num_files];
    header=new RAMSES_Header[opt.num_files];

#ifdef USEMPI
    MPI_Status status;
    Particle *Pbuf;
    //for parallel io
    Int_t Nlocalbuf,*Nbuf, *Nreadbuf,*nreadoffset;
    Int_t *Nlocalthreadbuf,Nlocaltotalbuf;
    int *irecv, sendTask,recvTask,irecvflag, *mpi_irecvflag;
    MPI_Request *mpi_request;
    Int_t *mpi_nsend_baryon;
    if (opt.iBaryonSearch) mpi_nsend_baryon=new Int_t[NProcs*NProcs];

    Nbuf=new Int_t[NProcs];
    nreadoffset=new Int_t[opt.nsnapread];

    if (ThisTask<opt.nsnapread)
    {
        //to temporarily store data from gadget file
        Pbuf=new Particle[BufSize*NProcs];
        Nreadbuf=new Int_t[opt.num_files];
        for (int j=0;j<NProcs;j++) Nbuf[j]=0;
        for (int j=0;j<opt.num_files;j++) Nreadbuf[j]=0;

        //to determine which files the thread should read
        ireadfile=new int[opt.num_files];
        for (i=0;i<opt.num_files;i++) ireadfile[i]=0;
        int nread=opt.num_files/opt.nsnapread;
        int niread=ThisTask*nread,nfread=(ThisTask+1)*nread;
        if (ThisTask==opt.nsnapread-1) nfread=opt.num_files;
        for (i=niread;i<nfread;i++) ireadfile[i]=1;
        ifirstfile=niread;
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

#ifndef MPIREDUCEMEM
    MPIDomainExtentRAMSES(opt);
    if (NProcs>1) {
    MPIDomainDecompositionRAMSES(opt);
    MPIInitialDomainDecomposition();
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (ThisTask<opt.nsnapread) {
#endif

    //first read cosmological information
    sprintf(buf1,"info_%s.txt",opt.fname);
    Finfo.open(buf1, ios::in);
    getline(Finfo,stringbuf);
    getline(Finfo,stringbuf);
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].levelmin;
    getline(Finfo,stringbuf);
    getline(Finfo,stringbuf);
    getline(Finfo,stringbuf);
    getline(Finfo,stringbuf);
    getline(Finfo,stringbuf);
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].BoxSize;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].time;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].aexp;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].HubbleParam;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].Omegam;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].OmegaLambda;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].Omegak;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].Omegab;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].scale_l;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].scale_d;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].scale_t;
    //convert boxsize to comoving mpc
    header[ifirstfile].BoxSize*=header[ifirstfile].scale_l/3.08e24/header[ifirstfile].aexp;
    getline(Finfo,stringbuf);
    Finfo>>stringbuf>>orderingstring;
    getline(Finfo,stringbuf);
    Finfo.close();

    opt.p=header[ifirstfile].BoxSize;
    opt.a=exp(header[ifirstfile].aexp);
    opt.Omega_m=header[ifirstfile].Omegam;
    opt.Omega_Lambda=header[ifirstfile].OmegaLambda;
    opt.Omega_b=header[ifirstfile].Omegab;
    opt.h=header[ifirstfile].HubbleParam;
    opt.Omega_cdm=opt.Omega_m-opt.Omega_b;
    //Hubble flow
    if (opt.comove) aadjust=1.0;
    else aadjust=opt.a;
    Hubble=opt.h*opt.H*sqrt((1-opt.Omega_m-opt.Omega_Lambda)*pow(aadjust,-2.0)+opt.Omega_m*pow(aadjust,-3.0)+opt.Omega_Lambda);
    opt.rhobg=3.*Hubble*Hubble/(8.0*M_PI*opt.G)*opt.Omega_m;
    //if opt.virlevel<0, then use virial overdensity based on Bryan and Norman 1998 virialization level is given by
    if (opt.virlevel<0) 
    {
        Double_t bnx=-((1-opt.Omega_m-opt.Omega_Lambda)*pow(aadjust,-2.0)+opt.Omega_Lambda)/((1-opt.Omega_m-opt.Omega_Lambda)*pow(aadjust,-2.0)+opt.Omega_m*pow(aadjust,-3.0)+opt.Omega_Lambda);
        opt.virlevel=(18.0*M_PI*M_PI+82.0*bnx-39*bnx*bnx)/opt.Omega_m;
    }
    //adjust length scale so that convert from 0 to 1 (box units) to Mpc comoving
    opt.L*=header[ifirstfile].scale_l/3.08e24/header[ifirstfile].aexp;
    //adjust velocity scale to that ramses is converted to km/s from which you can convert again;
    opt.V*=header[ifirstfile].scale_l/header[ifirstfile].scale_t*1e-5;
    mscale=opt.M/opt.h;lscale=opt.L/opt.h*aadjust;lvscale=opt.L/opt.h*opt.a;rhoscale=mscale/(lscale*lscale*lscale);
    //ignore hubble flow
    Hubbleflow=0.;
    //for (int j=0;j<NPARTTYPES;j++) nbodies+=opt.numpart[j];
    cout<<"Particle system contains "<<nbodies<<" particles (of interest) at is at time "<<opt.a<<" in a box of size "<<opt.p<<endl;
    cout<<"Cosmology (h,Omega_m,Omega_cdm,Omega_b,Omega_L) = ("<< opt.h<<","<<opt.Omega_m<<","<<opt.Omega_cdm<<","<<opt.Omega_b<<","<<opt.Omega_Lambda<<")"<<endl;
    N_DM=opt.numpart[DARKTYPE];
    LN=(opt.p*lscale/pow(N_DM,1.0/3.0));

    //grab from the first particle file the dimensions of the arrays and also the number of cpus (should be number of files)
    sprintf(buf1,"part_%s.out00001",opt.fname);
    Fpart[ifirstfile].open(buf1, ios::binary|ios::in);
    RAMSES_fortran_read(Fpart[ifirstfile],header[ifirstfile].nfiles);
    RAMSES_fortran_read(Fpart[ifirstfile],header[ifirstfile].ndim);
    //adjust the number of files 
    opt.num_files=header[ifirstfile].nfiles;
    Fpart[ifirstfile].close();
#ifdef USEMPI
    //now read tasks prepped and can read files to send information
    }
#endif 
    

    //if not only gas being searched open particle data
    count2=bcount2=0;
    if (opt.partsearchtype!=PSTGAS) {
    if (ThisTask<opt.nsnapread) {
    //read particle files consists of positions,velocities, mass, id, and level (along with ages and met if some flags set)
    for (i=0;i<opt.num_files;i++) if (ireadfile[i]) {
        sprintf(buf1,"part_%s.out%05d",opt.fname,i+1);
        sprintf(buf2,"part_%s.out",opt.fname);
        if (FileExists(buf1)) sprintf(buf,"%s",buf1);
        else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
        Fpart[i].open(buf, ios::binary|ios::in);
        Fpartvel[i].open(buf, ios::binary|ios::in);
        Fpartmass[i].open(buf, ios::binary|ios::in);
        Fpartid[i].open(buf, ios::binary|ios::in);
        Fpartlevel[i].open(buf, ios::binary|ios::in);
        Fpartage[i].open(buf, ios::binary|ios::in);
        Fpartmet[i].open(buf, ios::binary|ios::in);

        //skip header information in each file save for number in the file
        //@{
        byteoffset=0;
        byteoffset+=RAMSES_fortran_skip(Fpart[i]);
        byteoffset+=RAMSES_fortran_skip(Fpart[i]);
        //store number of particles locally in file
        byteoffset+=RAMSES_fortran_read(Fpart[i],header[i].npartlocal);
        byteoffset+=RAMSES_fortran_skip(Fpart[i],5);

        //byteoffset now stores size of header offset for particles 
        Fpartvel[i].seekg(byteoffset,ios::cur);
        Fpartmass[i].seekg(byteoffset,ios::cur);
        Fpartid[i].seekg(byteoffset,ios::cur);
        Fpartlevel[i].seekg(byteoffset,ios::cur);
        Fpartage[i].seekg(byteoffset,ios::cur);
        Fpartmet[i].seekg(byteoffset,ios::cur);
        //skip positions
        for(idim=0;idim<header[ifirstfile].ndim;idim++)
        {
            RAMSES_fortran_skip(Fpartvel[i]);
            RAMSES_fortran_skip(Fpartmass[i]);
            RAMSES_fortran_skip(Fpartid[i]);
            RAMSES_fortran_skip(Fpartlevel[i]);
            RAMSES_fortran_skip(Fpartage[i]);
            RAMSES_fortran_skip(Fpartmet[i]);
        }
        //skip velocities
        for(idim=0;idim<header[ifirstfile].ndim;idim++)
        {
            RAMSES_fortran_skip(Fpartmass[i]);
            RAMSES_fortran_skip(Fpartid[i]);
            RAMSES_fortran_skip(Fpartlevel[i]);
            RAMSES_fortran_skip(Fpartage[i]);
            RAMSES_fortran_skip(Fpartmet[i]);
        }
        //skip mass
        RAMSES_fortran_skip(Fpartid[i]);
        RAMSES_fortran_skip(Fpartlevel[i]);
        RAMSES_fortran_skip(Fpartage[i]);
        RAMSES_fortran_skip(Fpartmet[i]);
        //skip ids;
        RAMSES_fortran_skip(Fpartlevel[i]);
        RAMSES_fortran_skip(Fpartage[i]);
        RAMSES_fortran_skip(Fpartmet[i]);
        //skip levels
        RAMSES_fortran_skip(Fpartage[i]);
        RAMSES_fortran_skip(Fpartmet[i]);
        //skip ages
        RAMSES_fortran_skip(Fpartmet[i]);
        //@}

        //data loaded into memory in chunks
        chunksize=nchunk=header[i].npartlocal;
        xtempchunk=new RAMSESFLOAT[3*chunksize];
        vtempchunk=new RAMSESFLOAT[3*chunksize];
        mtempchunk=new RAMSESFLOAT[chunksize];
        idvalchunk=new RAMSESIDTYPE[chunksize];
        agetempchunk=new RAMSESFLOAT[chunksize];
        mettempchunk=new RAMSESFLOAT[chunksize];
        for(idim=0;idim<header[ifirstfile].ndim;idim++)
        {
            RAMSES_fortran_read(Fpart[i],&xtempchunk[idim*nchunk]);
            RAMSES_fortran_read(Fpartvel[i],&vtempchunk[idim*nchunk]);
        }
        RAMSES_fortran_read(Fpartmass[i],mtempchunk);
        RAMSES_fortran_read(Fpartid[i],idvalchunk);
        for (int nn=0;nn<nchunk;nn++) {
            xtemp[0]=xtempchunk[nn];xtemp[1]=xtempchunk[nn+nchunk];xtemp[2]=xtempchunk[nn+2*nchunk];
            vtemp[0]=vtempchunk[nn];vtemp[1]=vtempchunk[nn+nchunk];vtemp[2]=vtempchunk[nn+2*nchunk];
            idval=idvalchunk[nn];
            for (int kk=0;kk<3;kk++) {xtemp[kk]=LittleFLOAT(xtemp[kk]);vtemp[kk]=LittleFLOAT(vtemp[kk]);}
#ifndef NOMASS
            mtemp=mtempchunk[nn];
#else
            mtemp=1.0;
#endif
            if (ageval==0 && idval>0) typeval=DARKTYPE;
            else if (idval>0) typeval=STARTYPE;
            else typeval=BHTYPE;
#ifdef USEMPI
            //determine processor this particle belongs on based on its spatial position
            ibuf=MPIGetParticlesProcessor(xtemp[0],xtemp[1],xtemp[2]);
#endif

            if (opt.partsearchtype==PSTALL) {
#ifdef USEMPI
                Pbuf[ibuf*BufSize+Nbuf[ibuf]]=Particle(mtemp*mscale,
                    xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
                    vtemp[0]*opt.V+Hubbleflow*xtemp[0],
                    vtemp[1]*opt.V+Hubbleflow*xtemp[1],
                    vtemp[2]*opt.V+Hubbleflow*xtemp[2],
                    count2,typeval);
                Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPID(idval);
                //assume that first sphblock is internal energy
                //ensure that store number of particles to be sent to the threads involved with reading snapshot files
                if(ibuf<opt.nsnapread&&ibuf!=ThisTask) Nreadbuf[ibuf]++;
                Nbuf[ibuf]++;
                if (ibuf==ThisTask) {
                    Nbuf[ibuf]--;
                    Part[Nlocal++]=Pbuf[ibuf*BufSize+Nbuf[ibuf]];
                }
                else {
                    //before a simple send was done because only Task zero was reading the data
                    //but now if ibuf<opt.nsnapread, care must be taken.
                    //blocking sends that are matched by non-blocking receives
                    if(Nbuf[ibuf]==BufSize&&ibuf>=opt.nsnapread) {
                        MPI_Send(&Nbuf[ibuf], 1, MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
                        MPI_Send(&Pbuf[ibuf*BufSize],sizeof(Particle)*Nbuf[ibuf],MPI_BYTE,ibuf,ibuf,MPI_COMM_WORLD);
                        Nbuf[ibuf]=0;
                    }
                    else if (Nbuf[ibuf]==BufSize&&ibuf<opt.nsnapread) {
                        Nbuf[ibuf]=0;
                    }
                }
#else
                Part[count2]=Particle(mtemp*mscale,
                    xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
                    vtemp[0]*opt.V+Hubbleflow*xtemp[0],
                    vtemp[1]*opt.V+Hubbleflow*xtemp[1],
                    vtemp[2]*opt.V+Hubbleflow*xtemp[2],
                    count2,typeval);
                Part[count2].SetPID(idval);
#endif
                count2++;
            }
            else if (opt.partsearchtype==PSTDARK) {
                if (!(typeval==STARTYPE||typeval==BHTYPE)) {
#ifdef USEMPI
                    Pbuf[ibuf*BufSize+Nbuf[ibuf]]=Particle(mtemp*mscale,
                        xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
                        vtemp[0]*opt.V+Hubbleflow*xtemp[0],
                        vtemp[1]*opt.V+Hubbleflow*xtemp[1],
                        vtemp[2]*opt.V+Hubbleflow*xtemp[2],
                        count2,DARKTYPE);
                    Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPID(idval);
                    //ensure that store number of particles to be sent to other reading threads
                    if(ibuf<opt.nsnapread&&ibuf!=ThisTask) Nreadbuf[ibuf]++;
                    Nbuf[ibuf]++;
                    //now determine what to do with the particle, local or must send
                    if (ibuf==ThisTask) {
                        Nbuf[ibuf]--;
                        Part[Nlocal++]=Pbuf[ibuf*BufSize+Nbuf[ibuf]];
                    }
                    else {
                        //if belongs to another mpi thread then see if buffer is full and send with appropriate flag
                        if(Nbuf[ibuf]==BufSize&&ibuf>=opt.nsnapread) {
                            MPI_Ssend(&Nbuf[ibuf], 1, MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
                            MPI_Ssend(&Pbuf[ibuf*BufSize],sizeof(Particle)*Nbuf[ibuf],MPI_BYTE,ibuf,ibuf,MPI_COMM_WORLD);
                            Nbuf[ibuf]=0;
                        }
                        else if (Nbuf[ibuf]==BufSize&&ibuf<opt.nsnapread) {
                            Nbuf[ibuf]=0;
                        }
                    }
#else
                    Part[count2]=Particle(mtemp*mscale,
                        xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
                        vtemp[0]*opt.V+Hubbleflow*xtemp[0],
                        vtemp[1]*opt.V+Hubbleflow*xtemp[1],
                        vtemp[2]*opt.V+Hubbleflow*xtemp[2],
                        count2,typeval);
                    Part[count2].SetPID(idval);
#endif
                    count2++;
                }
                else if (opt.iBaryonSearch) {
#ifdef USEMPI
                    Pbuf[ibuf*BufSize+Nbuf[ibuf]]=Particle(mtemp*mscale,
                        xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
                        vtemp[0]*opt.V+Hubbleflow*xtemp[0],
                        vtemp[1]*opt.V+Hubbleflow*xtemp[1],
                        vtemp[2]*opt.V+Hubbleflow*xtemp[2],
                        count2);
                    Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPID(idval);
                    if (typeval==STARTYPE) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetType(STARTYPE);
                    else if (typeval==BHTYPE) Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetType(BHTYPE);
                    //ensure that store number of particles to be sent to the reading threads
                    if(ibuf<opt.nsnapread&&ibuf!=ThisTask) Nreadbuf[ibuf]++;
                    Nbuf[ibuf]++;
                    if (ibuf==ThisTask) {
                        Nbuf[ibuf]--;
                        Pbaryons[Nlocalbaryon[0]++]=Pbuf[ibuf*BufSize+Nbuf[ibuf]];
                        if (k==RAMSESSTARTYPE) Nlocalbaryon[2]++;
                        else if (k==RAMSESSINKTYPE) Nlocalbaryon[3]++;
                    }
                    else {
                        if(Nbuf[ibuf]==BufSize&&ibuf>=opt.nsnapread) {
                            MPI_Ssend(&Nbuf[ibuf], 1, MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
                            MPI_Ssend(&Pbuf[ibuf*BufSize],sizeof(Particle)*Nbuf[ibuf],MPI_BYTE,ibuf,ibuf,MPI_COMM_WORLD);
                            Nbuf[ibuf]=0;
                        }
                        else if (Nbuf[ibuf]==BufSize&&ibuf<opt.nsnapread) {
                            Nbuf[ibuf]=0;
                        }
                    }
#else
                    Pbaryons[bcount2]=Particle(mtemp*mscale,
                        xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
                        vtemp[0]*opt.V+Hubbleflow*xtemp[0],
                        vtemp[1]*opt.V+Hubbleflow*xtemp[1],
                        vtemp[2]*opt.V+Hubbleflow*xtemp[2],
                        count2,typeval);
                    Pbaryons[bcount2].SetPID(idval);
#endif
                    bcount2++;
                }
            }
            else if (opt.partsearchtype==PSTSTAR) {
                if (typeval==STARTYPE) {
#ifdef USEMPI
                    //if using MPI, determine proccessor and place in ibuf, store particle in particle buffer and if buffer full, broadcast data
                    //unless ibuf is 0, then just store locally
                    Pbuf[ibuf*BufSize+Nbuf[ibuf]]=Particle(mtemp*mscale,
                        xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
                        vtemp[0]*opt.V+Hubbleflow*xtemp[0],
                        vtemp[1]*opt.V+Hubbleflow*xtemp[1],
                        vtemp[2]*opt.V+Hubbleflow*xtemp[2],
                        count2,STARTYPE);
                    //ensure that store number of particles to be sent to the reading threads
                    if(ibuf<opt.nsnapread&&ibuf!=ThisTask) Nreadbuf[ibuf]++;
                    Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPID(idval);
                    Nbuf[ibuf]++;
                    if (ibuf==ThisTask) {
                        Nbuf[ibuf]--;
                        Part[Nlocal++]=Pbuf[ibuf*BufSize+Nbuf[ibuf]];
                    }
                    else {
                        if(Nbuf[ibuf]==BufSize&&ibuf>=opt.nsnapread) {
                            MPI_Ssend(&Nbuf[ibuf], 1, MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
                            MPI_Ssend(&Pbuf[ibuf*BufSize],sizeof(Particle)*Nbuf[ibuf],MPI_BYTE,ibuf,ibuf,MPI_COMM_WORLD);
                            Nbuf[ibuf]=0;
                        }
                        else if (Nbuf[ibuf]==BufSize&&ibuf<opt.nsnapread) {
                            Nbuf[ibuf]=0;
                        }
                    }
#else
                Part[count2]=Particle(mtemp*mscale,
                    xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
                    vtemp[0]*opt.V+Hubbleflow*xtemp[0],
                    vtemp[1]*opt.V+Hubbleflow*xtemp[1],
                    vtemp[2]*opt.V+Hubbleflow*xtemp[2],
                    count2,typeval);
                Part[count2].SetPID(idval);
#endif
                    count2++;
                }
            }
        }
        delete[] xtempchunk;
        delete[] vtempchunk;
        delete[] mtempchunk;
        delete[] idvalchunk;
        delete[] agetempchunk;
        delete[] mettempchunk;
        Fpart[i].close();
        Fpartvel[i].close();
        Fpartmass[i].close();
        Fpartid[i].close();
        Fpartlevel[i].close();
        Fpartage[i].close();
        Fpartmet[i].close();
    }
#ifdef USEMPI
    //once finished reading the file if there are any particles left in the buffer broadcast them
    for(ibuf = opt.nsnapread; ibuf < NProcs; ibuf++)
    {
        MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
        if (Nbuf[ibuf]>0) {
            MPI_Ssend(&Pbuf[ibuf*BufSize], sizeof(Particle)*Nbuf[ibuf], MPI_BYTE, ibuf, ibuf, MPI_COMM_WORLD);
            Nbuf[ibuf]=0;
            //last broadcast with Nbuf[ibuf]=0 so that receiver knows no more particles are to be broadcast
            MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t,ibuf,ibuf+NProcs,MPI_COMM_WORLD);
        }
    }
#endif
    }
#ifdef USEMPI
    //if not reading information than waiting to receive information
    else {
        //for all threads not reading snapshots, simply receive particles as necessary from all threads involved with reading the data
        //first determine which threads are going to send information to this thread.
        for (i=0;i<opt.nsnapread;i++) if (irecv[i]) {
            mpi_irecvflag[i]=0;
            MPI_Irecv(&Nlocalthreadbuf[i], 1, MPI_Int_t, i, ThisTask+NProcs, MPI_COMM_WORLD, &mpi_request[i]);
        }
        Nlocaltotalbuf=0;
        //non-blocking receives for the number of particles one expects to receive
        do {
            irecvflag=0;
            for (i=0;i<opt.nsnapread;i++) if (irecv[i]) {
                if (mpi_irecvflag[i]==0) {
                    //test if a request has been sent for a Recv call by one of the read threads
                    MPI_Test(&mpi_request[i], &mpi_irecvflag[i], &status);
                    if (mpi_irecvflag[i]) {
                        if (Nlocalthreadbuf[i]>0) {
                            MPI_Recv(&Part[Nlocal],sizeof(Particle)*Nlocalthreadbuf[i],MPI_BYTE,i,ThisTask, MPI_COMM_WORLD,&status);
                            Nlocal+=Nlocalthreadbuf[i];
                            Nlocaltotalbuf+=Nlocalthreadbuf[i];
                            mpi_irecvflag[i]=0;
                            MPI_Irecv(&Nlocalthreadbuf[i], 1, MPI_Int_t, i, ThisTask+NProcs, MPI_COMM_WORLD, &mpi_request[i]);
                        }
                        else {
                            irecv[i]=0;
                        }
                    }
                }
            }
            for (i=0;i<opt.nsnapread;i++) irecvflag+=irecv[i];
        } while(irecvflag>0);
        //now that data is local, must adjust data iff a separate baryon search is required. 
        if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
            for (i=0;i<Nlocal;i++) {
                k=Part[i].GetType();
                if (!(k==GASTYPE||k==STARTYPE||k==BHTYPE)) Part[i].SetID(0);
                else {
                    Nlocalbaryon[0]++;
                    if  (k==GASTYPE) {Part[i].SetID(1);Nlocalbaryon[1]++;}
                    else if  (k==STARTYPE) {Part[i].SetID(2);Nlocalbaryon[2]++;}
                    else if  (k==BHTYPE) {Part[i].SetID(3);Nlocalbaryon[3]++;}
                }
            }
            //sorted so that dark matter particles first, baryons after
            qsort(Part,Nlocal, sizeof(Particle), IDCompare);
            Nlocal-=Nlocalbaryon[0];
            //index type separated
            for (i=0;i<Nlocal;i++) Part[i].SetID(i);
            for (i=0;i<Nlocalbaryon[0];i++) Part[i+Nlocal].SetID(i+Nlocal);
            //finally, need to move baryons forward by the Export Factor * Nlocal as need that extra buffer to copy data two and from mpi threads
//#ifndef MPIREDUCE
//            for (i=Nlocalbaryon[0]-1;i>=0;i--) Part[i+(Int_t)(Nlocal*MPIExportFac)]=Part[i+Nlocal];
//#endif
        }
    }
#endif
    }

    //if gas searched in some fashion then load amr/hydro data
    if (opt.partsearchtype==PSTGAS||opt.partsearchtype==PSTALL||(opt.partsearchtype==PSTDARK&&opt.iBaryonSearch)) {
    if (ThisTask<opt.nsnapread) {
    for (i=0;i<opt.num_files;i++) if (ireadfile[i]) {
        sprintf(buf1,"amr_.out%s%05d",opt.fname,i+1);
        sprintf(buf2,"amr_.out%s",opt.fname);
        if (FileExists(buf1)) sprintf(buf,"%s",buf1);
        else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
        Famr[i].open(buf, ios::binary|ios::in);
        sprintf(buf1,"hydro_%s.out%05d",opt.fname,i+1);
        sprintf(buf2,"hydro_%s.out",opt.fname);
        if (FileExists(buf1)) sprintf(buf,"%s",buf1);
        else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
        Fhydro[i].open(buf, ios::binary|ios::in);
        //read some of the amr header till get to number of cells in current file
        //@{
        byteoffset=0;
        byteoffset+=RAMSES_fortran_read(Famr[i],header[i].ndim);
        header[i].twotondim=pow(2,header[i].ndim);
        Famr[i].read((char*)&dummy, sizeof(dummy));
        Famr[i].read((char*)&header[i].nx, sizeof(int));
        Famr[i].read((char*)&header[i].ny, sizeof(int));
        Famr[i].read((char*)&header[i].nz, sizeof(int));
        Famr[i].read((char*)&dummy, sizeof(dummy));
        byteoffset+=RAMSES_fortran_read(Famr[i],header[i].nlevelmax);
        byteoffset+=RAMSES_fortran_read(Famr[i],header[i].ngridmax);
        byteoffset+=RAMSES_fortran_read(Famr[i],header[i].nboundary);
        byteoffset+=RAMSES_fortran_read(Famr[i],header[i].npart[RAMSESGASTYPE]);

        //then skip the rest
        for (j=0;j<14;j++) RAMSES_fortran_skip(Famr[i]);
        if (lmin>header[i].nlevelmax) lmin=header[i].nlevelmax;
        if (lmax<header[i].nlevelmax) lmax=header[i].nlevelmax;
        //@}
        //read header info from hydro files
        //@{
        RAMSES_fortran_skip(Fhydro[i]);
        RAMSES_fortran_read(Fhydro[i],header[i].nvarh);
        RAMSES_fortran_skip(Fhydro[i]);
        RAMSES_fortran_skip(Fhydro[i]);
        RAMSES_fortran_skip(Fhydro[i]);
        RAMSES_fortran_read(Fhydro[i],header[i].gamma_index);
        //@}
    }
    for (i=0;i<opt.num_files;i++) if (ireadfile[i]) {
        //then apparently read ngridlevels, which appears to be an array storing the number of grids at a given level
        ngridlevel=new int[header[i].nlevelmax];
        ngridfile=new int[(1+header[i].nboundary)*header[i].nlevelmax];
        RAMSES_fortran_read(Famr[i],ngridlevel);
        for (j=0;j<header[i].nlevelmax;j++) ngridfile[j]=ngridlevel[j];
        //skip some more 
        RAMSES_fortran_skip(Famr[i]);
        //if nboundary>0 then need two skip twice then read ngridbound
        if(header[i].nboundary>0) {
            ngridbound=new int[header[i].nboundary*header[i].nlevelmax];
            RAMSES_fortran_skip(Famr[i]);
            RAMSES_fortran_skip(Famr[i]);
            //ngridbound is an array of some sort but I don't see what it is used for
            RAMSES_fortran_read(Famr[i],ngridbound);
            for (j=0;j<header[i].nlevelmax;j++) ngridfile[header[i].nlevelmax+j]=ngridbound[j];
        }
        //skip some more 
        RAMSES_fortran_skip(Famr[i],2);
        //if odering list in info is bisection need to skip more 
        if (orderingstring==string("bisection")) RAMSES_fortran_skip(Famr[i],5);
        else RAMSES_fortran_skip(Famr[i],4);

        for (k=0;k<header[i].nboundary+1;k++) {
            for (j=0;j<header[i].nlevelmax;j++) {
                //first read amr for positions
                chunksize=nchunk=ngridfile[k*header[i].nlevelmax+j];
                if (chunksize>0) {
                    xtempchunk=new RAMSESFLOAT[3*chunksize];
                    //store son value in icell
                    icellchunk=new int[header[i].twotondim*chunksize];
                    //skip grid index, next index and prev index. 
                    RAMSES_fortran_skip(Famr[i],3);
                    //now read grid centre
                    for (idim=0;idim<header[i].ndim;idim++) {
                        RAMSES_fortran_read(Famr[i],&xtempchunk[idim*chunksize]);
                    }
                    //skip father index, then neighbours index
                    RAMSES_fortran_skip(Famr[i],1+2*header[i].ndim);
                    //read son index to determine if a cell in a specific grid is at the highest resolution and needs to be represented by a particle
                    for (idim=0;idim<header[i].twotondim;idim++) {
                        RAMSES_fortran_read(Famr[i],&icellchunk[idim*chunksize]);
                    }
                    //skip cpu map and refinement map (2^ndim*2)
                    RAMSES_fortran_skip(Famr[i],2*header[i].twotondim);
                }
                RAMSES_fortran_skip(Fhydro[i]);
                //then read hydro for other variables (first is density, then velocity, then pressure, then metallicity )
                if (chunksize>0) {
                    hydrotempchunk=new RAMSESFLOAT[chunksize*header[i].twotondim*header[i].nvarh];
                    //first read velocities (for 2 cells per number of dimensions (ie: cell corners?))
                    for (idim=0;idim<header[i].twotondim;idim++) {
                        for (ivar=0;ivar<header[i].nvarh;ivar++) {
                            RAMSES_fortran_read(Fhydro[i],&hydrotempchunk[idim*chunksize*header[i].nvarh+ivar*chunksize]);
                            for (igrid=0;igrid<chunksize;igrid++) {
                                //once we have looped over all the hydro data then can start actually storing it into the particle structures
                                if (ivar==header[i].nvarh-1) {
                                    //if cell has no internal cells or at maximum level produce a particle
                                    if (icellchunk[idim*chunksize+igrid]==0 || j==header[i].nlevelmax-1) {
                                        //first suggestion is to add some jitter to the particle positions
                                        double dx = pow(0.5, j);
                                        int ix, iy, iz;
                                        //below assumes three dimensions with 8 corners (? maybe cells) per grid
                                        iz = idim/4;
                                        iy = (idim - (4*iz))/2;
                                        ix = idim - (2*iy) - (4*iz);
                                        // Calculate absolute coordinates + jitter, and generate particle
                                        xpos[0] = ((((float)rand()/(float)RAND_MAX) * header[i].BoxSize * dx) +(header[i].BoxSize * (xtempchunk[igrid] + (double(ix)-0.5) * dx )) - (header[i].BoxSize*dx/2.0)) ;
                                        xpos[1] = ((((float)rand()/(float)RAND_MAX) * header[i].BoxSize * dx) +(header[i].BoxSize * (xtempchunk[igrid+1*chunksize] + (double(iy)-0.5) * dx )) - (header[i].BoxSize*dx/2.0)) ;
                                        xpos[2] = ((((float)rand()/(float)RAND_MAX) * header[i].BoxSize * dx) +(header[i].BoxSize * (xtempchunk[igrid+2*chunksize] + (double(iz)-0.5) * dx )) - (header[i].BoxSize*dx/2.0)) ;
#ifdef USEMPI
                                        //determine processor this particle belongs on based on its spatial position
                                        ibuf=MPIGetParticlesProcessor(xpos[0],xpos[1],xpos[2]);
#endif

                                        vpos[0]=hydrotempchunk[idim*chunksize*header[i].nvarh+1*chunksize+igrid];
                                        vpos[1]=hydrotempchunk[idim*chunksize*header[i].nvarh+2*chunksize+igrid];
                                        vpos[2]=hydrotempchunk[idim*chunksize*header[i].nvarh+3*chunksize+igrid];
                                        mtemp=dx*dx*dx*hydrotempchunk[idim*chunksize*header[i].nvarh+0*chunksize+igrid];
                                        //the self energy P/rho is given by 
                                        utemp=hydrotempchunk[idim*chunksize*header[i].nvarh+4*chunksize+igrid]/hydrotempchunk[idim*chunksize*header[i].nvarh+0*chunksize+igrid]/(header[i].gamma_index-1.0);
                                        rhotemp=hydrotempchunk[idim*chunksize*header[i].nvarh+0*chunksize+igrid]*rhoscale;
                                        Ztemp=hydrotempchunk[idim*chunksize*header[i].nvarh+5*chunksize+igrid];
                                        if (opt.partsearchtype==PSTALL) {
#ifdef USEMPI
                                            Pbuf[ibuf*BufSize+Nbuf[ibuf]]=Particle(mtemp*mscale,
                                                xpos[0]*lscale,xpos[1]*lscale,xpos[2]*lscale,
                                                vpos[0]*opt.V+Hubbleflow*xpos[0],
                                                vpos[1]*opt.V+Hubbleflow*xpos[1],
                                                vpos[2]*opt.V+Hubbleflow*xpos[2],
                                                count2,GASTYPE);
                                            Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPID(idval);
#ifdef GASON
                                            Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetU(utemp);
                                            Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetSPHDen(rhotemp);
#ifdef STARON
                                            Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetZmet(Ztemp);
#endif
#endif
                                            //ensure that store number of particles to be sent to the threads involved with reading snapshot files
                                            if(ibuf<opt.nsnapread&&ibuf!=ThisTask) Nreadbuf[ibuf]++;
                                            Nbuf[ibuf]++;
                                            if (ibuf==ThisTask) {
                                                Nbuf[ibuf]--;
                                                Part[Nlocal++]=Pbuf[ibuf*BufSize+Nbuf[ibuf]];
                                            }
                                            else {
                                                //before a simple send was done because only Task zero was reading the data
                                                //but now if ibuf<opt.nsnapread, care must be taken.
                                                //blocking sends that are matched by non-blocking receives
                                                if(Nbuf[ibuf]==BufSize&&ibuf>=opt.nsnapread) {
                                                    MPI_Send(&Nbuf[ibuf], 1, MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
                                                    MPI_Send(&Pbuf[ibuf*BufSize],sizeof(Particle)*Nbuf[ibuf],MPI_BYTE,ibuf,ibuf,MPI_COMM_WORLD);
                                                    Nbuf[ibuf]=0;
                                                }
                                                else if (Nbuf[ibuf]==BufSize&&ibuf<opt.nsnapread) {
                                                    Nbuf[ibuf]=0;
                                                }
                                            }
#else
                                            Part[count2]=Particle(mtemp*mscale,
                                                xpos[0]*lscale,xpos[1]*lscale,xpos[2]*lscale,
                                                vpos[0]*opt.V+Hubbleflow*xpos[0],
                                                vpos[1]*opt.V+Hubbleflow*xpos[1],
                                                vpos[2]*opt.V+Hubbleflow*xpos[2],
                                                count2,GASTYPE);
                                            Part[count2].SetPID(idval);
#ifdef GASON
                                            Part[count2].SetU(utemp);
                                            Part[count2].SetSPHDen(rhotemp);
#ifdef STARON
                                            Part[count2].SetZmet(Ztemp);
#endif
#endif

#endif
                                            count2++;
                                        }
                                        else if (opt.partsearchtype==PSTDARK&&opt.iBaryonSearch) {
#ifdef USEMPI
                                            Pbuf[ibuf*BufSize+Nbuf[ibuf]]=Particle(mtemp*mscale,
                                                xpos[0]*lscale,xpos[1]*lscale,xpos[2]*lscale,
                                                vpos[0]*opt.V+Hubbleflow*xpos[0],
                                                vpos[1]*opt.V+Hubbleflow*xpos[1],
                                                vpos[2]*opt.V+Hubbleflow*xpos[2],
                                                count2,GASTYPE);
                                            Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetPID(idval);
#ifdef GASON
                                            Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetU(utemp);
                                            Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetSPHDen(rhotemp);
#ifdef STARON
                                            Pbuf[ibuf*BufSize+Nbuf[ibuf]].SetZmet(Ztemp);
#endif
#endif
                                        //ensure that store number of particles to be sent to the reading threads
                                        if(ibuf<opt.nsnapread&&ibuf!=ThisTask) Nreadbuf[ibuf]++;
                                        Nbuf[ibuf]++;
                                        if (ibuf==ThisTask) {
                                            Nbuf[ibuf]--;
                                            Pbaryons[Nlocalbaryon[0]++]=Pbuf[ibuf*BufSize+Nbuf[ibuf]];
                                            Nlocalbaryon[1]++;
                                        }
                                        else {
                                            if(Nbuf[ibuf]==BufSize&&ibuf>=opt.nsnapread) {
                                                MPI_Ssend(&Nbuf[ibuf], 1, MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
                                                MPI_Ssend(&Pbuf[ibuf*BufSize],sizeof(Particle)*Nbuf[ibuf],MPI_BYTE,ibuf,ibuf,MPI_COMM_WORLD);
                                                Nbuf[ibuf]=0;
                                            }
                                            else if (Nbuf[ibuf]==BufSize&&ibuf<opt.nsnapread) {
                                                Nbuf[ibuf]=0;
                                            }
                                        }
#else
                                            Pbaryons[bcount2]=Particle(mtemp*mscale,
                                                xpos[0]*lscale,xpos[1]*lscale,xpos[2]*lscale,
                                                vpos[0]*opt.V+Hubbleflow*xpos[0],
                                                vpos[1]*opt.V+Hubbleflow*xpos[1],
                                                vpos[2]*opt.V+Hubbleflow*xpos[2],
                                                count2,GASTYPE);
                                            Pbaryons[bcount2].SetPID(idval);
#ifdef GASON
                                            Pbaryons[bcount2].SetU(utemp);
                                            Pbaryons[bcount2].SetSPHDen(rhotemp);
#ifdef STARON
                                            Pbaryons[bcount2].SetZmet(Ztemp);
#endif
#endif

#endif
                                        bcount2++;
                                    }
                                }
                                }
                            }
                        }
                    }
                }
                if (chunksize>0) {
                    delete[] xtempchunk;
                    delete[] hydrotempchunk;
                }
            }
        }
    }
#ifdef USEMPI
    //once finished reading the file if there are any particles left in the buffer broadcast them
    for(ibuf = opt.nsnapread; ibuf < NProcs; ibuf++)
    {
        MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
        if (Nbuf[ibuf]>0) {
            MPI_Ssend(&Pbuf[ibuf*BufSize], sizeof(Particle)*Nbuf[ibuf], MPI_BYTE, ibuf, ibuf, MPI_COMM_WORLD);
            Nbuf[ibuf]=0;
            //last broadcast with Nbuf[ibuf]=0 so that receiver knows no more particles are to be broadcast
            MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t,ibuf,ibuf+NProcs,MPI_COMM_WORLD);
        }
    }
#endif
    }//end of reading task
#ifdef USEMPI
    //if not reading information than waiting to receive information
    else {
        //for all threads not reading snapshots, simply receive particles as necessary from all threads involved with reading the data
        //first determine which threads are going to send information to this thread.
        for (i=0;i<opt.nsnapread;i++) if (irecv[i]) {
            mpi_irecvflag[i]=0;
            MPI_Irecv(&Nlocalthreadbuf[i], 1, MPI_Int_t, i, ThisTask+NProcs, MPI_COMM_WORLD, &mpi_request[i]);
        }
        Nlocaltotalbuf=0;
        //non-blocking receives for the number of particles one expects to receive
        do {
            irecvflag=0;
            for (i=0;i<opt.nsnapread;i++) if (irecv[i]) {
                if (mpi_irecvflag[i]==0) {
                    //test if a request has been sent for a Recv call by one of the read threads
                    MPI_Test(&mpi_request[i], &mpi_irecvflag[i], &status);
                    if (mpi_irecvflag[i]) {
                        if (Nlocalthreadbuf[i]>0) {
                            MPI_Recv(&Part[Nlocal],sizeof(Particle)*Nlocalthreadbuf[i],MPI_BYTE,i,ThisTask, MPI_COMM_WORLD,&status);
                            Nlocal+=Nlocalthreadbuf[i];
                            Nlocaltotalbuf+=Nlocalthreadbuf[i];
                            mpi_irecvflag[i]=0;
                            MPI_Irecv(&Nlocalthreadbuf[i], 1, MPI_Int_t, i, ThisTask+NProcs, MPI_COMM_WORLD, &mpi_request[i]);
                        }
                        else {
                            irecv[i]=0;
                        }
                    }
                }
            }
            for (i=0;i<opt.nsnapread;i++) irecvflag+=irecv[i];
        } while(irecvflag>0);
        //now that data is local, must adjust data iff a separate baryon search is required. 
        if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
            for (i=0;i<Nlocal;i++) {
                k=Part[i].GetType();
                if (!(k==GASTYPE||k==STARTYPE||k==BHTYPE)) Part[i].SetID(0);
                else {
                    Nlocalbaryon[0]++;
                    if  (k==GASTYPE) {Part[i].SetID(1);Nlocalbaryon[1]++;}
                    else if  (k==STARTYPE) {Part[i].SetID(2);Nlocalbaryon[2]++;}
                    else if  (k==BHTYPE) {Part[i].SetID(3);Nlocalbaryon[3]++;}
                }
            }
            //sorted so that dark matter particles first, baryons after
            qsort(Part,Nlocal, sizeof(Particle), IDCompare);
            Nlocal-=Nlocalbaryon[0];
            //index type separated
            for (i=0;i<Nlocal;i++) Part[i].SetID(i);
            for (i=0;i<Nlocalbaryon[0];i++) Part[i+Nlocal].SetID(i+Nlocal);
            //finally, need to move baryons forward by the Export Factor * Nlocal as need that extra buffer to copy data two and from mpi threads
//#ifndef MPIREDUCE
//            for (i=Nlocalbaryon[0]-1;i>=0;i--) Part[i+(Int_t)(Nlocal*MPIExportFac)]=Part[i+Nlocal];
//#endif
        }
    }
#endif
    }//end of check if gas loaded

//#ifdef USEMPI
//    }
//#endif
}
