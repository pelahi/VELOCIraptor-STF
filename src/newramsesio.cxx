/*! \file newramsesio.cxx
 *  \brief this file contains routines to read the NewCluster HDF5 variant of RAMSES snapshots
 *
 * Each snapshot is stored as a single part_%s.h5 file holding dark matter, star, sink, cloud and
 * tracer particles, rather than the ncpu-way split binary files used by the original ramses io
 * (see \ref ramsesio.cxx). Only dark matter and star particles are read here; gas cells, stored in a
 * companion cell_%s.h5 file, are not yet supported.
 *
 * Selected at runtime via the New_ramsesio cfg flag (opt.inewramsesio), see \ref GetParamFile in ui.cxx
 * and the dispatch in \ref io.cxx.
 */

#include "stf.h"

#include "endianutils.h"
#include "newramsesitems.h"

#ifdef USEHDF
#include "hdfitems.h"

///subset of the fields of the "data" compound dataset (under /dm or /star) that is actually read;
///both groups share this layout, additional star-only fields (birth_time, metallicity, ...) are
///present in the file but not yet consumed here
struct NewRamsesPartRecord {
    double position_x, position_y, position_z;
    float  velocity_x, velocity_y, velocity_z;
    float  mass;
    int    identity;
};

static hid_t NewRamsesPartH5Type()
{
    hid_t t = H5Tcreate(H5T_COMPOUND, sizeof(NewRamsesPartRecord));
    H5Tinsert(t, "position_x", HOFFSET(NewRamsesPartRecord, position_x), H5T_NATIVE_DOUBLE);
    H5Tinsert(t, "position_y", HOFFSET(NewRamsesPartRecord, position_y), H5T_NATIVE_DOUBLE);
    H5Tinsert(t, "position_z", HOFFSET(NewRamsesPartRecord, position_z), H5T_NATIVE_DOUBLE);
    H5Tinsert(t, "velocity_x", HOFFSET(NewRamsesPartRecord, velocity_x), H5T_NATIVE_FLOAT);
    H5Tinsert(t, "velocity_y", HOFFSET(NewRamsesPartRecord, velocity_y), H5T_NATIVE_FLOAT);
    H5Tinsert(t, "velocity_z", HOFFSET(NewRamsesPartRecord, velocity_z), H5T_NATIVE_FLOAT);
    H5Tinsert(t, "mass",       HOFFSET(NewRamsesPartRecord, mass),       H5T_NATIVE_FLOAT);
    H5Tinsert(t, "identity",   HOFFSET(NewRamsesPartRecord, identity),   H5T_NATIVE_INT);
    return t;
}

///builds the full part_*.h5 path from the New_ramsesio_filename cfg option (the directory holding the
///hdf5 files, independent of opt.fname which the original binary reader uses) and the snapshot number
///passed via -t (opt.ramsessnapname), exactly as the original reader combines opt.fname+opt.ramsessnapname.
///This way only -t needs to change between snapshots; the cfg file itself stays the same.
static void NewRamsesPartFilePath(Options &opt, char *buf)
{
    if (opt.newramsesfname == NULL) {
        cerr<<"Error: New_ramsesio requires New_ramsesio_filename to be set in the cfg file "
            <<"(path to the directory containing part_*.h5), terminating"<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,9);
#else
        exit(9);
#endif
    }
    sprintf(buf, "%s/part_%s.h5", opt.newramsesfname, opt.ramsessnapname);
}

#ifdef USEMPI
///position-only record, used by the MPI domain-decomposition pre-pass helpers where velocity/mass/id
///are not needed
struct NewRamsesPosRecord {
    double position_x, position_y, position_z;
};

static hid_t NewRamsesPosH5Type()
{
    hid_t t = H5Tcreate(H5T_COMPOUND, sizeof(NewRamsesPosRecord));
    H5Tinsert(t, "position_x", HOFFSET(NewRamsesPosRecord, position_x), H5T_NATIVE_DOUBLE);
    H5Tinsert(t, "position_y", HOFFSET(NewRamsesPosRecord, position_y), H5T_NATIVE_DOUBLE);
    H5Tinsert(t, "position_z", HOFFSET(NewRamsesPosRecord, position_z), H5T_NATIVE_DOUBLE);
    return t;
}
#endif

///reads [offset,offset+count) particles of one type out of group "dm" or "star" in chunks of
///NEWRAMSESCHUNKSIZE, converts to physical units and appends the particles to dest (Part or Pbaryons),
///routing them to the appropriate MPI domain when compiled with MPI support
///the finest-refinement-level DM particle mass, in the same raw (pre-mscale) code units as
///buf[k].mass, computed exactly as dmp_mass is in ramsesio.cxx (see RAMSES_get_nbodies there):
///for a uniform-resolution box of Neff^3 DM particles, each carries this fraction of OmegaM-OmegaB.
///Zoom-in simulations also contain coarser, more massive "buffer zone" DM particles (8x/64x/512x...
///this mass, one step per refinement level dropped) outside the high-res region; the original binary
///reader has always excluded them by this same mass-equality check (within 1e-5 relative tolerance)
///plus family==1, but this newer HDF5 reader had no equivalent filter, silently including them.
static inline Double_t NewRamsesDMPMass(const Options &opt, double omegam, double omegab)
{
    return 1.0/((Double_t)opt.Neff*opt.Neff*opt.Neff) * (omegam-omegab)/omegam;
}
static inline bool NewRamsesIsFineDM(Double_t rawmass, Double_t dmpmass)
{
    return fabs((rawmass-dmpmass)/dmpmass) < 1e-5;
}

///counts, without materializing any Particle objects, how many entries in the "dm" group's
///finest refinement level match dmpmass -- used by NewRAMSES_get_nbodies to get an exact
///a-priori count for allocation, mirroring the pre-scan RAMSES_get_nbodies does for the
///original binary reader.
static Int_t NewRamsesCountFineDM(hid_t Fhdf, Int_t ndmtotal, Double_t dmpmass)
{
    if (ndmtotal<=0) return 0;
    hid_t group     = H5Gopen2(Fhdf, "dm", H5P_DEFAULT);
    hid_t dataset   = H5Dopen2(group, "data", H5P_DEFAULT);
    hid_t filespace = H5Dget_space(dataset);
    //memory compound type containing only the "mass" field -- HDF5 only reads that field off disk
    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(float));
    H5Tinsert(memtype, "mass", 0, H5T_NATIVE_FLOAT);

    vector<float> buf;
    Int_t nread=0, nmatch=0;
    while (nread<ndmtotal) {
        Int_t thischunk = min((Int_t)NEWRAMSESCHUNKSIZE, ndmtotal-nread);
        hsize_t start = nread, hcount = thischunk;
        safe_hdf5<herr_t>(H5Sselect_hyperslab, filespace, H5S_SELECT_SET, &start, (const hsize_t*)NULL, &hcount, (const hsize_t*)NULL);
        hid_t memspace = H5Screate_simple(1, &hcount, NULL);
        buf.resize(thischunk);
        safe_hdf5<herr_t>(H5Dread, dataset, memtype, memspace, filespace, H5P_DEFAULT, buf.data());
        H5Sclose(memspace);
        for (Int_t k=0;k<thischunk;k++) if (NewRamsesIsFineDM((Double_t)buf[k], dmpmass)) nmatch++;
        nread+=thischunk;
    }
    H5Tclose(memtype);
    H5Sclose(filespace);
    H5Dclose(dataset);
    H5Gclose(group);
    return nmatch;
}

static void NewRamsesReadParticleGroup(
    Options &opt, hid_t Fhdf, const char *groupname, int ptype,
    Int_t offset, Int_t count,
    Double_t mscale, Double_t lscale, Double_t velscale, Double_t Hubbleflow,
    Particle *dest, Int_t &destcount,
    int *ireadtask, const Int_t BufSize, Int_t *Nbuf, Particle *Pbuf,
    Int_t *Nreadbuf, vector<Particle> *Preadbuf,
    Double_t dmpmass=-1)
{
    if (count<=0) return;
    hid_t group     = H5Gopen2(Fhdf, groupname, H5P_DEFAULT);
    hid_t dataset   = H5Dopen2(group, "data", H5P_DEFAULT);
    hid_t filespace = H5Dget_space(dataset);
    hid_t memtype   = NewRamsesPartH5Type();

    vector<NewRamsesPartRecord> buf;
    Int_t nread=0;
    while (nread<count) {
        Int_t thischunk = min((Int_t)NEWRAMSESCHUNKSIZE, count-nread);
        hsize_t start = offset+nread, hcount = thischunk;
        safe_hdf5<herr_t>(H5Sselect_hyperslab, filespace, H5S_SELECT_SET, &start, (const hsize_t*)NULL, &hcount, (const hsize_t*)NULL);
        hid_t memspace = H5Screate_simple(1, &hcount, NULL);
        buf.resize(thischunk);
        safe_hdf5<herr_t>(H5Dread, dataset, memtype, memspace, filespace, H5P_DEFAULT, buf.data());
        H5Sclose(memspace);

        for (Int_t k=0;k<thischunk;k++) {
            //exclude coarser zoom-in buffer-zone DM particles, matching the original binary reader
            if (dmpmass>=0 && !NewRamsesIsFineDM((Double_t)buf[k].mass, dmpmass)) continue;
            //raw code-unit positions (0 to 1); MPIGetParticlesProcessor expects these units (mpi_domain
            //boundaries are set up in code units by MPIDomainExtentRAMSES/MPIDomainDecompositionWithTree),
            //while the Particle itself is constructed with the usual lscale-scaled (kpc) positions
            Double_t rawx=buf[k].position_x, rawy=buf[k].position_y, rawz=buf[k].position_z;
            Double_t x=rawx*lscale, y=rawy*lscale, z=rawz*lscale;
            Double_t vx=buf[k].velocity_x*velscale+Hubbleflow*rawx;
            Double_t vy=buf[k].velocity_y*velscale+Hubbleflow*rawy;
            Double_t vz=buf[k].velocity_z*velscale+Hubbleflow*rawz;
            Double_t mass=buf[k].mass*mscale;
            Int_t idval=buf[k].identity;
#ifdef USEMPI
            int ibuf=MPIGetParticlesProcessor(opt,rawx,rawy,rawz);
            Int_t ibufindex=ibuf*BufSize+Nbuf[ibuf];
            Pbuf[ibufindex]=Particle(mass,x,y,z,vx,vy,vz,destcount,ptype);
            Pbuf[ibufindex].SetPID(idval);
            Nbuf[ibuf]++;
            MPIAddParticletoAppropriateBuffer(opt, ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, destcount, dest, Nreadbuf, Preadbuf);
#else
            dest[destcount]=Particle(mass,x,y,z,vx,vy,vz,destcount,ptype);
            dest[destcount].SetPID(idval);
            destcount++;
#endif
        }
        nread+=thischunk;
    }
    H5Tclose(memtype);
    H5Sclose(filespace);
    H5Dclose(dataset);
    H5Gclose(group);
}

Int_t NewRAMSES_get_nbodies(char *fname, int ptype, Options &opt)
{
    char buf[2000];
    NewRamsesPartFilePath(opt, buf);
    if (!FileExists(buf)) {
        printf("Error. Can't find new ramses particle data as `%s'\n\n", buf);
        exit(9);
    }
    if (ptype==PSTGAS || ptype==PSTBH) {
        cout<<"Warning: new ramses io does not yet support gas cells or sink particles, treating this particle count as 0"<<endl;
    }

    hid_t Fhdf = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
    long long ndm=0, nstar=0;
    if (ptype==PSTALL||ptype==PSTDARK) {
        ndm = read_attribute<long long>(Fhdf, "dm/size");
        double omegam = read_attribute<double>(Fhdf, "omega_m");
        double omegab = read_attribute<double>(Fhdf, "omega_b");
        Double_t dmpmass = NewRamsesDMPMass(opt, omegam, omegab);
        Int_t ndmfine = NewRamsesCountFineDM(Fhdf, (Int_t)ndm, dmpmass);
        if (ndmfine<ndm) cout<<"New ramses io: excluding "<<(ndm-ndmfine)<<" coarser zoom-buffer DM particles ("
            <<ndmfine<<" of "<<ndm<<" match the finest-level mass)"<<endl;
        ndm = ndmfine;
    }
    if (ptype==PSTALL||ptype==PSTSTAR) nstar = read_attribute<long long>(Fhdf, "star/size");
    H5Fclose(Fhdf);

    for (int j=0;j<NPARTTYPES;j++) opt.numpart[j]=0;
    if (ptype==PSTALL || ptype==PSTDARK) opt.numpart[DARKTYPE]=ndm;
    if (ptype==PSTALL || ptype==PSTSTAR) opt.numpart[STARTYPE]=nstar;

    Int_t nbodies=0;
    if (ptype==PSTALL) nbodies = ndm+nstar;
    else if (ptype==PSTDARK) nbodies = ndm;
    else if (ptype==PSTSTAR) nbodies = nstar;
    //PSTGAS, PSTBH not yet supported, return 0
    return nbodies;
}

/// Reads dark matter and star particles from the NewCluster HDF5 ramses format. See \ref ReadRamses for
/// the equivalent reader for the original binary format; this function mirrors its overall structure
/// (cosmology setup, particle-type dispatch, MPI read-thread buffering) but sources particle data from
/// a single HDF5 file read in contiguous row-range chunks instead of ncpu-many fortran binary files.
void ReadNewRamses(Options &opt, vector<Particle> &Part, const Int_t nbodies, Particle *&Pbaryons, Int_t nbaryons)
{
    char buf[2000];
    NewRamsesPartFilePath(opt, buf);
    Double_t mscale, lscale, velscale, Hubble, Hubbleflow=0.;
    double aadjust;
    Int_t count2=0, bcount2=0;

    if (!FileExists(buf)) {
        cerr<<"Error. Can't find new ramses particle data as `"<<buf<<"'"<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,9);
#else
        exit(9);
#endif
    }
    if (opt.partsearchtype==PSTGAS || opt.partsearchtype==PSTBH) {
        cerr<<"Error: new ramses io does not yet support gas cell or sink particles, terminating"<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,9);
#else
        exit(9);
#endif
    }

    //--- cosmology / units: every task reads these small root attributes independently ---
    hid_t Fhdf = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
    double boxlen  = read_attribute<double>(Fhdf, "boxlen");
    double aexp    = read_attribute<double>(Fhdf, "aexp");
    double H0      = read_attribute<double>(Fhdf, "H0");
    double omegam  = read_attribute<double>(Fhdf, "omega_m");
    double omegal  = read_attribute<double>(Fhdf, "omega_l");
    double omegab  = read_attribute<double>(Fhdf, "omega_b");
    double unit_l  = read_attribute<double>(Fhdf, "unit_l");
    double unit_d  = read_attribute<double>(Fhdf, "unit_d");
    double unit_t  = read_attribute<double>(Fhdf, "unit_t");
    Int_t ndmtotal   = (Int_t)read_attribute<long long>(Fhdf, "dm/size");
    Int_t nstartotal = (Int_t)read_attribute<long long>(Fhdf, "star/size");
    H5Fclose(Fhdf);

    opt.a            = aexp;
    opt.Omega_m      = omegam;
    opt.Omega_Lambda = omegal;
    opt.Omega_b      = omegab;
    opt.h            = H0/100.0;
    opt.Omega_cdm    = opt.Omega_m-opt.Omega_b;
    //set hubble unit to km/s/kpc
    opt.H = 0.1;
    //set Gravity to value for kpc (km/s)^2 / solar mass
    opt.G = 4.30211349e-6;
    //and for now fix the units
    opt.lengthtokpc=opt.velocitytokms=opt.masstosolarmass=1.0;

    //Hubble flow
    if (opt.comove) aadjust=1.0;
    else aadjust=opt.a;
    CalcOmegak(opt);
    Hubble=GetHubble(opt, aadjust);
    CalcCriticalDensity(opt, aadjust);
    CalcBackgroundDensity(opt, aadjust);
    CalcVirBN98(opt,aadjust);
    //if opt.virlevel<0, then use virial overdensity based on Bryan and Norman 1997
    if (opt.virlevel<0) opt.virlevel=opt.virBN98;
    PrintCosmology(opt);

    //box size in comoving kpc/h, matching the convention used by the original (binary) ramses reader
    opt.p = boxlen * (unit_l/3.086e21) / opt.a * (H0/100.0);
    opt.lengthinputconversion = opt.p;
    //convert ramses code velocities to km/s
    opt.velocityinputconversion = unit_l/unit_t*1e-5;

    //convert mass from code units to Solar Masses
    mscale = unit_d * (unit_l*unit_l*unit_l) / 1.988e33;
    //convert length from code units to kpc
    lscale = unit_l/3.086e21;
    velscale = opt.velocityinputconversion;

    //ignore hubble flow, matches ramsesio.cxx
    Hubbleflow=0.;

    //interparticle spacing (assuming a uniform resolution box)
    opt.ellxscale = lscale/(double)opt.Neff;
    //excludes coarser zoom-in buffer-zone DM particles, matching the original binary reader
    //(see NewRamsesDMPMass); ndmtotal itself is left as the raw HDF5 row count so the read
    //loops below still scan every row to test its mass, they just skip writing the ones that fail
    Double_t dmpmass = NewRamsesDMPMass(opt, omegam, omegab);

    cout<<"Particle system contains "<<nbodies<<" particles (of interest) at is at time "<<opt.a<<" in a box of size "<<opt.p<<endl;

    bool readdm       = (opt.partsearchtype==PSTALL || opt.partsearchtype==PSTDARK);
    bool readstar     = (opt.partsearchtype==PSTALL || opt.partsearchtype==PSTSTAR || (opt.partsearchtype==PSTDARK && opt.iBaryonSearch));
    bool startobaryon = (opt.partsearchtype==PSTDARK && opt.iBaryonSearch);

#ifndef USEMPI
    int ThisTask=0, NProcs=1;
    Fhdf = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (readdm) {
        NewRamsesReadParticleGroup(opt, Fhdf, "dm", DARKTYPE, 0, ndmtotal, mscale, lscale, velscale, Hubbleflow,
            Part.data(), count2, NULL, 0, NULL, NULL, NULL, NULL, dmpmass);
    }
    if (readstar) {
        Particle *stardest = startobaryon? Pbaryons : Part.data();
        Int_t &starcount = startobaryon? bcount2 : count2;
        NewRamsesReadParticleGroup(opt, Fhdf, "star", STARTYPE, 0, nstartotal, mscale, lscale, velscale, Hubbleflow,
            stardest, starcount, NULL, 0, NULL, NULL, NULL, NULL);
    }
    H5Fclose(Fhdf);
#else
    int ThisTask,NProcs;
    MPI_Comm_size(MPI_COMM_WORLD,&NProcs);
    MPI_Comm_rank(MPI_COMM_WORLD,&ThisTask);
    MPI_Comm mpi_comm_read;
    Int_t BufSize=opt.mpiparticlebufsize;
    Int_t *Nbuf, *Nreadbuf=NULL, *Nlocalthreadbuf=NULL;
    int *ireadtask, *readtaskID, *irecv=NULL, *mpi_irecvflag=NULL;
    MPI_Request *mpi_request=NULL;
    vector<Particle> *Preadbuf=NULL;
    Particle *Pbuf=NULL;
    Int_t *mpi_nsend_readthread=NULL, *mpi_nsend_readthread_baryon=NULL;

    //there is only a single physical file; (ab)use opt.num_files/opt.nsnapread so that
    //MPIDistributeReadTasks assigns opt.nsnapread tasks to independently read non-overlapping
    //row ranges of the dm/star datasets in parallel
    opt.num_files = opt.nsnapread;
    Nbuf=new Int_t[NProcs];
    ireadtask=new int[NProcs];
    readtaskID=new int[opt.nsnapread];
    MPIDistributeReadTasks(opt,ireadtask,readtaskID);
    MPI_Comm_split(MPI_COMM_WORLD, (ireadtask[ThisTask]>=0), ThisTask, &mpi_comm_read);

    Nlocal=0;
    if (opt.iBaryonSearch) Nlocalbaryon[0]=0;

    if (ireadtask[ThisTask]>=0) {
        Pbuf=new Particle[BufSize*NProcs];
        Nreadbuf=new Int_t[opt.nsnapread];
        for (int j=0;j<NProcs;j++) Nbuf[j]=0;
        for (int j=0;j<opt.nsnapread;j++) Nreadbuf[j]=0;
        if (opt.nsnapread>1) {
            Preadbuf=new vector<Particle>[opt.nsnapread];
            for (int j=0;j<opt.nsnapread;j++) Preadbuf[j].reserve(BufSize);
            mpi_nsend_readthread=new Int_t[opt.nsnapread*opt.nsnapread];
            if (opt.iBaryonSearch) mpi_nsend_readthread_baryon=new Int_t[opt.nsnapread*opt.nsnapread];
        }

        Int_t k = ireadtask[ThisTask];
        Fhdf = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (readdm) {
            Int_t nper=ndmtotal/opt.nsnapread, off=k*nper, cnt=(k==opt.nsnapread-1)?(ndmtotal-off):nper;
            NewRamsesReadParticleGroup(opt, Fhdf, "dm", DARKTYPE, off, cnt, mscale, lscale, velscale, Hubbleflow,
                Part.data(), Nlocal, ireadtask, BufSize, Nbuf, Pbuf, Nreadbuf, Preadbuf, dmpmass);
            if (opt.nsnapread>1) {
                MPI_Allgather(Nreadbuf, opt.nsnapread, MPI_Int_t, mpi_nsend_readthread, opt.nsnapread, MPI_Int_t, mpi_comm_read);
                MPISendParticlesBetweenReadThreads(opt, Preadbuf, Part.data(), ireadtask, readtaskID, Pbaryons, mpi_comm_read, mpi_nsend_readthread, mpi_nsend_readthread_baryon);
                for (int j=0;j<opt.nsnapread;j++) Nreadbuf[j]=0;
            }
        }
        if (readstar) {
            Int_t nper=nstartotal/opt.nsnapread, off=k*nper, cnt=(k==opt.nsnapread-1)?(nstartotal-off):nper;
            Particle *stardest = startobaryon? Pbaryons : Part.data();
            Int_t &starcount = startobaryon? Nlocalbaryon[0] : Nlocal;
            NewRamsesReadParticleGroup(opt, Fhdf, "star", STARTYPE, off, cnt, mscale, lscale, velscale, Hubbleflow,
                stardest, starcount, ireadtask, BufSize, Nbuf, Pbuf, Nreadbuf, Preadbuf);
            if (opt.nsnapread>1) {
                MPI_Allgather(Nreadbuf, opt.nsnapread, MPI_Int_t, mpi_nsend_readthread, opt.nsnapread, MPI_Int_t, mpi_comm_read);
                MPISendParticlesBetweenReadThreads(opt, Preadbuf, Part.data(), ireadtask, readtaskID, Pbaryons, mpi_comm_read, mpi_nsend_readthread, mpi_nsend_readthread_baryon);
                for (int j=0;j<opt.nsnapread;j++) Nreadbuf[j]=0;
            }
        }
        H5Fclose(Fhdf);

        //flush any particles still buffered for non-reading tasks
        for (int ib=0; ib<NProcs; ib++) if (ireadtask[ib]<0) {
            MPI_Ssend(&Nbuf[ib],1,MPI_Int_t, ib, ib+NProcs, MPI_COMM_WORLD);
            if (Nbuf[ib]>0) {
                MPI_Ssend(&Pbuf[ib*BufSize], sizeof(Particle)*Nbuf[ib], MPI_BYTE, ib, ib, MPI_COMM_WORLD);
                Nbuf[ib]=0;
                //last send with Nbuf[ib]=0 so that the receiver knows no more particles are coming
                MPI_Ssend(&Nbuf[ib],1,MPI_Int_t,ib,ib+NProcs,MPI_COMM_WORLD);
            }
        }
        if (opt.nsnapread>1) {
            MPI_Allgather(Nreadbuf, opt.nsnapread, MPI_Int_t, mpi_nsend_readthread, opt.nsnapread, MPI_Int_t, mpi_comm_read);
            MPISendParticlesBetweenReadThreads(opt, Preadbuf, Part.data(), ireadtask, readtaskID, Pbaryons, mpi_comm_read, mpi_nsend_readthread, mpi_nsend_readthread_baryon);
        }
    }
    else {
        Nlocalthreadbuf=new Int_t[opt.nsnapread];
        irecv=new int[opt.nsnapread];
        mpi_irecvflag=new int[opt.nsnapread];
        for (int i=0;i<opt.nsnapread;i++) irecv[i]=1;
        mpi_request=new MPI_Request[opt.nsnapread];
        MPIReceiveParticlesFromReadThreads(opt,Pbuf,Part.data(),readtaskID, irecv, mpi_irecvflag, Nlocalthreadbuf, mpi_request,Pbaryons);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    //update box size to match the convention expected by downstream code (see ramsesio.cxx)
    opt.p*=opt.a/opt.h;
    //store how to convert input internal energies to physical output internal energies;
    //not applicable here as no hydro quantities are read
    opt.internalenergyinputconversion = 1.0;

#ifdef USEMPI
    MPI_Bcast(&(opt.p),sizeof(opt.p),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Comm_free(&mpi_comm_read);
    if (opt.nsnapread>1) {
        if (mpi_nsend_readthread) delete[] mpi_nsend_readthread;
        if (mpi_nsend_readthread_baryon) delete[] mpi_nsend_readthread_baryon;
        if (ireadtask[ThisTask]>=0) delete[] Preadbuf;
    }
    delete[] Nbuf;
    if (ireadtask[ThisTask]>=0) {
        delete[] Nreadbuf;
        delete[] Pbuf;
    }
    else {
        delete[] Nlocalthreadbuf;
        delete[] irecv;
        delete[] mpi_irecvflag;
        delete[] mpi_request;
    }
    delete[] ireadtask;
    delete[] readtaskID;
#endif
}

#ifdef USEMPI
///counts, for a single group ("dm" or "star"), how many of its [offset,offset+count) particles fall in
///each MPI domain, used by \ref NewRAMSES_CountNumInDomain
static void NewRamsesCountParticleGroup(
    Options &opt, hid_t Fhdf, const char *groupname,
    Int_t offset, Int_t count,
    Int_t *Nbuf, Int_t *Nbaryonbuf, bool tobaryon)
{
    if (count<=0) return;
    hid_t group     = H5Gopen2(Fhdf, groupname, H5P_DEFAULT);
    hid_t dataset   = H5Dopen2(group, "data", H5P_DEFAULT);
    hid_t filespace = H5Dget_space(dataset);
    hid_t memtype   = NewRamsesPosH5Type();

    vector<NewRamsesPosRecord> buf;
    Int_t nread=0;
    while (nread<count) {
        Int_t thischunk = min((Int_t)NEWRAMSESCHUNKSIZE, count-nread);
        hsize_t start = offset+nread, hcount = thischunk;
        safe_hdf5<herr_t>(H5Sselect_hyperslab, filespace, H5S_SELECT_SET, &start, (const hsize_t*)NULL, &hcount, (const hsize_t*)NULL);
        hid_t memspace = H5Screate_simple(1, &hcount, NULL);
        buf.resize(thischunk);
        safe_hdf5<herr_t>(H5Dread, dataset, memtype, memspace, filespace, H5P_DEFAULT, buf.data());
        H5Sclose(memspace);

        for (Int_t k=0;k<thischunk;k++) {
            //raw code-unit positions (0 to 1); matches what MPIGetParticlesProcessor is fed with
            //everywhere else in the ramses MPI domain-decomposition code
            int ibuf=MPIGetParticlesProcessor(opt,buf[k].position_x,buf[k].position_y,buf[k].position_z);
            if (tobaryon) Nbaryonbuf[ibuf]++;
            else Nbuf[ibuf]++;
        }
        nread+=thischunk;
    }
    H5Tclose(memtype);
    H5Sclose(filespace);
    H5Dclose(dataset);
    H5Gclose(group);
}

///Only ThisTask==0 actually reads (the whole file, not a partition), leaving Nbuf/Nbaryonbuf at zero on
///every other task; the caller's existing MPI_Allreduce(...,MPI_SUM,...) then recovers task 0's exact
///counts on every task. This sidesteps needing any extra MPI synchronisation here.
void NewRAMSES_CountNumInDomain(Options &opt, Int_t *Nbuf, Int_t *Nbaryonbuf)
{
    if (ThisTask!=0) return;
    char buf[2000];
    NewRamsesPartFilePath(opt, buf);
    hid_t Fhdf = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
    Int_t ndmtotal   = (Int_t)read_attribute<long long>(Fhdf, "dm/size");
    Int_t nstartotal = (Int_t)read_attribute<long long>(Fhdf, "star/size");

    bool readdm       = (opt.partsearchtype==PSTALL || opt.partsearchtype==PSTDARK);
    bool readstar     = (opt.partsearchtype==PSTALL || opt.partsearchtype==PSTSTAR || (opt.partsearchtype==PSTDARK && opt.iBaryonSearch));
    bool startobaryon = (opt.partsearchtype==PSTDARK && opt.iBaryonSearch);

    if (readdm) NewRamsesCountParticleGroup(opt, Fhdf, "dm", 0, ndmtotal, Nbuf, Nbaryonbuf, false);
    if (readstar) NewRamsesCountParticleGroup(opt, Fhdf, "star", 0, nstartotal, Nbuf, Nbaryonbuf, startobaryon);
    H5Fclose(Fhdf);
}

///reads a single group ("dm" or "star") into Part_mpi/xtempall/famtempall starting at count_mpi, used by
///\ref NewRAMSES_ReadForDomainTree
static void NewRamsesReadForTreeGroup(
    Options &opt, hid_t Fhdf, const char *groupname, int ptype,
    Int_t offset, Int_t count, Double_t mscale, Double_t lscale, Int_t nbodies,
    vector<Particle> &Part_mpi, Int_t &count_mpi, RAMSESFLOAT *xtempall, int *famtempall)
{
    if (count<=0) return;
    hid_t group     = H5Gopen2(Fhdf, groupname, H5P_DEFAULT);
    hid_t dataset   = H5Dopen2(group, "data", H5P_DEFAULT);
    hid_t filespace = H5Dget_space(dataset);
    hid_t memtype   = NewRamsesPartH5Type();

    vector<NewRamsesPartRecord> buf;
    Int_t nread=0;
    while (nread<count) {
        Int_t thischunk = min((Int_t)NEWRAMSESCHUNKSIZE, count-nread);
        hsize_t start = offset+nread, hcount = thischunk;
        safe_hdf5<herr_t>(H5Sselect_hyperslab, filespace, H5S_SELECT_SET, &start, (const hsize_t*)NULL, &hcount, (const hsize_t*)NULL);
        hid_t memspace = H5Screate_simple(1, &hcount, NULL);
        buf.resize(thischunk);
        safe_hdf5<herr_t>(H5Dread, dataset, memtype, memspace, filespace, H5P_DEFAULT, buf.data());
        H5Sclose(memspace);

        for (Int_t k=0;k<thischunk;k++) {
            Double_t mass = buf[k].mass*mscale;
            //raw code-unit positions (0 to 1); matches the units MPIGetParticlesProcessor is fed with
            //in MPIDomainDecompositionWithTree once domain boundaries are converted back to code units
            RAMSESFLOAT x=buf[k].position_x, y=buf[k].position_y, z=buf[k].position_z;
            Part_mpi[count_mpi]=Particle(mass, x*lscale, y*lscale, z*lscale, 0.,0.,0., count_mpi, ptype);
            xtempall[count_mpi]           = x;
            xtempall[count_mpi+nbodies]   = y;
            xtempall[count_mpi+nbodies*2] = z;
            famtempall[count_mpi] = ptype;
            count_mpi++;
        }
        nread+=thischunk;
    }
    H5Tclose(memtype);
    H5Sclose(filespace);
    H5Dclose(dataset);
    H5Gclose(group);
}

///Every task independently reads the whole (small, single-file) dataset, so that every task ends up
///with an identical, complete Part_mpi/xtempall/famtempall to build its own (consistent) copy of the
///domain-decomposition tree from -- unlike the original binary reader, nothing here is partitioned
///across tasks, avoiding the need for any extra MPI synchronisation of the result.
void NewRAMSES_ReadForDomainTree(Options &opt, vector<Particle> &Part_mpi, Int_t nbodies, Int_t &count_mpi, RAMSESFLOAT *xtempall, int *famtempall, Double_t &lscale)
{
    char buf[2000];
    NewRamsesPartFilePath(opt, buf);
    hid_t Fhdf = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
    double unit_l = read_attribute<double>(Fhdf, "unit_l");
    double unit_d = read_attribute<double>(Fhdf, "unit_d");
    lscale = unit_l/3.086e21;
    Double_t mscale = unit_d*(unit_l*unit_l*unit_l)/1.988e33;
    Int_t ndmtotal   = (Int_t)read_attribute<long long>(Fhdf, "dm/size");
    Int_t nstartotal = (Int_t)read_attribute<long long>(Fhdf, "star/size");

    count_mpi=0;
    //mirrors MPIDomainDecompositionWithTree: PSTDARK with iBaryonSearch is not handled here
    if (opt.partsearchtype==PSTALL || opt.partsearchtype==PSTDARK) {
        NewRamsesReadForTreeGroup(opt, Fhdf, "dm", DARKTYPE, 0, ndmtotal, mscale, lscale, nbodies, Part_mpi, count_mpi, xtempall, famtempall);
    }
    if (opt.partsearchtype==PSTALL || opt.partsearchtype==PSTSTAR) {
        NewRamsesReadForTreeGroup(opt, Fhdf, "star", STARTYPE, 0, nstartotal, mscale, lscale, nbodies, Part_mpi, count_mpi, xtempall, famtempall);
    }
    H5Fclose(Fhdf);
}
#endif

#else

Int_t NewRAMSES_get_nbodies(char *fname, int ptype, Options &opt)
{
    cerr<<"Error: New_ramsesio requires VELOCIraptor to be compiled with HDF5 support (USEHDF), terminating"<<endl;
    exit(9);
    return 0;
}

void ReadNewRamses(Options &opt, vector<Particle> &Part, const Int_t nbodies, Particle *&Pbaryons, Int_t nbaryons)
{
    cerr<<"Error: New_ramsesio requires VELOCIraptor to be compiled with HDF5 support (USEHDF), terminating"<<endl;
    exit(9);
}

#ifdef USEMPI
void NewRAMSES_CountNumInDomain(Options &opt, Int_t *Nbuf, Int_t *Nbaryonbuf)
{
    cerr<<"Error: New_ramsesio requires VELOCIraptor to be compiled with HDF5 support (USEHDF), terminating"<<endl;
    exit(9);
}

void NewRAMSES_ReadForDomainTree(Options &opt, vector<Particle> &Part_mpi, Int_t nbodies, Int_t &count_mpi, RAMSESFLOAT *xtempall, int *famtempall, Double_t &lscale)
{
    cerr<<"Error: New_ramsesio requires VELOCIraptor to be compiled with HDF5 support (USEHDF), terminating"<<endl;
    exit(9);
}
#endif

#endif
