/*! \file swiftinterface.cxx
 *  \brief this file contains routines that allow the velociraptor library to interface with the swift N-body code from within swift.
 */

#include "ioutils.h"
#include "logging.h"
#include "swiftinterface.h"
#include "timer.h"

#ifdef SWIFTINTERFACE

Options libvelociraptorOpt;
Options libvelociraptorOptbackup;
Options libvelociraptorOptextra[10];
//KDTree *mpimeshtree;
//Particle *mpimeshinfo;

/// \defgroup SWIFTCONFIGERRORS Errors for swift
//@{
///
#define SWIFTCONFIGOPTMISSING 8
#define SWIFTCONFIGOPTERROR 9
#define SWIFTCONFIGOPTCONFLICT 10
//@}

///check configuration options with swift
inline int ConfigCheckSwift(Options &opt, Swift::siminfo &s)
{
    if (opt.iBaryonSearch && !(opt.partsearchtype==PSTALL || opt.partsearchtype==PSTDARK)) {
        LOG_RANK0(error) << "Conflict in config file: both gas/star/etc particle type search AND the separate baryonic (gas,star,etc) search flag are on. Check config";
        return SWIFTCONFIGOPTCONFLICT;
    }
    if (opt.iBoundHalos && opt.iKeepFOF) {
        LOG_RANK0(error) << "Conflict in config file: Asking for Bound Field objects but also asking to keep the 3DFOF/then run 6DFOF. This is incompatible. Check config";
        return SWIFTCONFIGOPTCONFLICT;
    }
    if (opt.HaloMinSize==-1) opt.HaloMinSize=opt.MinSize;

    if (s.idarkmatter == false && opt.partsearchtype==PSTDARK) {
        LOG_RANK0(error) << "Conflict in config file: simulation contains no dark matter but searching only dark matter. This is incompatible. Check config";
        return SWIFTCONFIGOPTCONFLICT;
    }
    if (s.idarkmatter == false && opt.partsearchtype==PSTALL && opt.iBaryonSearch) {
        LOG_RANK0(error) << "Conflict in config file: simulation contains no dark matter but using dark matter to define links when searching all particles . This is incompatible. Check config";
        return SWIFTCONFIGOPTCONFLICT;
    }
    if (s.igas == false && opt.partsearchtype==PSTGAS) {
        LOG_RANK0(error) << "Conflict in config file: simulation contains no gas but searching only gas. This is incompatible. Check config";
        return SWIFTCONFIGOPTCONFLICT;
    }
    if ((s.istar == false && opt.partsearchtype==PSTSTAR)) {
        LOG_RANK0(error) << "Conflict in config file: simulation contains no gas but searching only stars or using stars as basis for links. This is incompatible. Check config";
        return SWIFTCONFIGOPTCONFLICT;
    }

#ifndef USEHDF
    if (opt.ibinaryout==OUTHDF){
        LOG_RANK0(error) << "Code not compiled with HDF output enabled. Recompile with this enabled or change Binary_output";
        return SWIFTCONFIGOPTERROR;
    }
#endif

#ifndef USEADIOS
    if (opt.ibinaryout==OUTADIOS){
        LOG_RANK0(error) << "Code not compiled with ADIOS output enabled. Recompile with this enabled or change Binary_output";
        return SWIFTCONFIGOPTERROR;
    }
#endif

    //call the general ConfigCheck now
    ConfigCheck(opt);
    return 1;
}

inline bool CheckSwiftPartType(int type)
{
#ifdef HIGHRES
    return (type != DARKTYPE && type != DARK2TYPE && type != GASTYPE && type != STARTYPE && type != BHTYPE);
#else
    return (type != DARKTYPE && type != GASTYPE && type != STARTYPE && type != BHTYPE);
#endif
}


int InitVelociraptor(Options &opt, char* configname, unitinfo u, siminfo s, const int numthreads)
{
    // if mpi invokved, init the velociraptor tasks and openmp threads
#ifdef USEMPI
    //find out how big the SPMD world is
    MPI_Comm_size(MPI_COMM_WORLD,&NProcs);
    //mpi_domain=new MPI_Domain[NProcs];
    mpi_nlocal=new Int_t[NProcs];
    mpi_nsend=new Int_t[NProcs*NProcs];
    mpi_ngroups=new Int_t[NProcs];
    mpi_nhalos=new Int_t[NProcs];
    //and this processes' rank is
    MPI_Comm_rank(MPI_COMM_WORLD,&ThisTask);
    //store MinSize as when using mpi prior to stitching use min of 2;
    MinNumMPI=2;
#else
    int ThisTask = 0;
#endif
#ifdef USEOPENMP
    omp_set_num_threads(numthreads);
#endif

    // logging level is reset correctly later after reading the user configuration
    vr::init_logging(vr::LogLevel::trace);
    gsl_set_error_handler_off();

    int iconfigflag;
    ///read the parameter file
    opt.pname = configname;
    LOG_RANK0(info) << "Initialising VELOCIraptor git revision " << velociraptor::git_sha1();
    LOG_RANK0(info) << "Reading VELOCIraptor config file...";
    GetParamFile(opt);
    //on the fly finding and using swift's mesh mpi decomposition
    opt.iontheflyfinding = true;
    opt.impiusemesh = true;
    ///check configuration
    iconfigflag = ConfigCheckSwift(opt, s);
    if (iconfigflag != 1) return iconfigflag;

    LOG_RANK0(info) << "Setting cosmology, units, sim stuff";
    ///set units, here idea is to convert internal units so that have kpc, km/s, solar mass
    ///\todo switch this so run in reasonable swift units and store conversion
    opt.lengthtokpc=u.lengthtokpc;
    opt.velocitytokms=u.velocitytokms;
    opt.masstosolarmass=u.masstosolarmass;

    //run in swift internal units, don't convert units
    opt.lengthinputconversion=1.0;
    opt.massinputconversion=1.0;
    opt.velocityinputconversion=1.0;
    opt.energyinputconversion=1.0;
    opt.SFRinputconversion=1.0;
    opt.metallicityinputconversion=1.0;
    opt.istellaragescalefactor = 1;

    //set cosmological parameters that do not change
    ///these should be in units of kpc, km/s, and solar mass
    opt.G=u.gravity;
    opt.H=u.hubbleunit;

    //set if cosmological
    opt.icosmologicalin = s.icosmologicalsim;

    //store a general mass unit, useful if running uniform box with single mass
    //and saving memory by not storing mass per particle.
#ifdef NOMASS
    opt.MassValue = s.mass_uniform_box;
    NOMASSCheck(opt);
#endif

    //write velociraptor configuration info, appending .configuration to the input config file and writing every config option
    opt.outname = configname;

    //store list of names that
    WriteVELOCIraptorConfig(opt);

#ifdef USEMPI
    //initialize the mpi write communicator to comm world;
    MPIInitWriteComm();
#endif

    LOG_RANK0(info) << "Finished initialising VELOCIraptor";

    //return the configuration flag value
    return iconfigflag;

}

int InitVelociraptor(char* configname, unitinfo u, siminfo s, const int numthreads)
{
	return InitVelociraptor(libvelociraptorOpt, configname, u, s, numthreads);
}

int InitVelociraptorExtra(const int iextra, char* configname, unitinfo u, siminfo s, const int numthreads)
{
	return InitVelociraptor(libvelociraptorOptextra[iextra], configname, u, s, numthreads);
}

void SetVelociraptorCosmology(cosmoinfo c)
{
    ///set cosmology
    libvelociraptorOpt.h=c.littleh;
    libvelociraptorOpt.Omega_m=c.Omega_m;
    libvelociraptorOpt.Omega_b=c.Omega_b;
    libvelociraptorOpt.Omega_cdm=c.Omega_cdm;
    libvelociraptorOpt.Omega_Lambda=c.Omega_Lambda;
    libvelociraptorOpt.Omega_r=c.Omega_r;
    libvelociraptorOpt.Omega_nu=c.Omega_nu;
    libvelociraptorOpt.Omega_k=c.Omega_k;
    libvelociraptorOpt.Omega_de = 0.0;
    libvelociraptorOpt.w_de=c.w_de;
    if (libvelociraptorOpt.w_de != -1) {
        libvelociraptorOpt.Omega_de = libvelociraptorOpt.Omega_Lambda;
        libvelociraptorOpt.Omega_Lambda = 0;
    }
    PrintCosmology(libvelociraptorOpt);
}

//set simulation state
void SetVelociraptorSimulationState(cosmoinfo c, siminfo s)
{
    double Hubble;
    ///set cosmology
    if (libvelociraptorOpt.icosmologicalin) SetVelociraptorCosmology(c);

    //set current scalefactor
    libvelociraptorOpt.a=c.atime;

    //set some sim information
    libvelociraptorOpt.p=s.period;
    libvelociraptorOpt.zoomlowmassdm=s.zoomhigresolutionmass;
    libvelociraptorOpt.icosmologicalin=s.icosmologicalsim;
    libvelociraptorOpt.ellxscale=s.interparticlespacing;
    libvelociraptorOpt.uinfo.eps*=libvelociraptorOpt.ellxscale;

    libvelociraptorOpt.uinfo.icalculatepotential=true;

    // Set mesh information.
    libvelociraptorOpt.spacedimension[0] = s.spacedimension[0];
    libvelociraptorOpt.spacedimension[1] = s.spacedimension[1];
    libvelociraptorOpt.spacedimension[2] = s.spacedimension[2];
    libvelociraptorOpt.cellwidth[0] = s.cellwidth[0];
    libvelociraptorOpt.cellwidth[1] = s.cellwidth[1];
    libvelociraptorOpt.cellwidth[2] = s.cellwidth[2];
    libvelociraptorOpt.icellwidth[0] = s.icellwidth[0];
    libvelociraptorOpt.icellwidth[1] = s.icellwidth[1];
    libvelociraptorOpt.icellwidth[2] = s.icellwidth[2];
    libvelociraptorOpt.numcells = s.numcells;
    libvelociraptorOpt.numcellsperdim = cbrt(s.numcells);
    libvelociraptorOpt.cellloc = s.cellloc;

    //assume above is in comoving is cosmological simulation and then need to correct to physical
    if (libvelociraptorOpt.icosmologicalin) {
        libvelociraptorOpt.p*=libvelociraptorOpt.a;
        libvelociraptorOpt.ellxscale*=libvelociraptorOpt.a;
        libvelociraptorOpt.uinfo.eps*=libvelociraptorOpt.a;

        CalcOmegak(libvelociraptorOpt);
        Hubble=GetHubble(libvelociraptorOpt, libvelociraptorOpt.a);
        CalcCriticalDensity(libvelociraptorOpt, libvelociraptorOpt.a);
        CalcBackgroundDensity(libvelociraptorOpt, libvelociraptorOpt.a);
        CalcVirBN98(libvelociraptorOpt,libvelociraptorOpt.a);
        if (libvelociraptorOpt.virlevel<0) libvelociraptorOpt.virlevel=libvelociraptorOpt.virBN98;

        for (auto i=0;i<3;i++) {
            libvelociraptorOpt.spacedimension[i] *= libvelociraptorOpt.a;
            libvelociraptorOpt.cellwidth[i] *= libvelociraptorOpt.a;
            libvelociraptorOpt.icellwidth[i] /= libvelociraptorOpt.a;
        }
    }
    else {
        libvelociraptorOpt.rhocrit=1.0;
        libvelociraptorOpt.rhobg=1.0;
    }
    PrintSimulationState(libvelociraptorOpt);

#ifdef USEMPI
    //if single halo, use minsize to initialize the old minimum number
    //else use the halominsize since if mpi and not single halo, halos localized to mpi domain for substructure search
    if (libvelociraptorOpt.iSingleHalo) MinNumOld=libvelociraptorOpt.MinSize;
    else MinNumOld=libvelociraptorOpt.HaloMinSize;
    mpi_period = libvelociraptorOpt.p;
#endif
}

///\todo interface with swift is comoving positions, period, peculiar velocities and physical self-energy
extern "C" Swift::vr_return_data InvokeVelociraptor(const int snapnum, char* outputname,
    cosmoinfo c, siminfo s,
    const size_t num_gravity_parts, const size_t num_hydro_parts, const size_t num_star_parts,
    struct swift_vel_part *swift_parts, int *cell_node_ids,
    const int numthreads, const int ireturngroupinfoflag, const int ireturnmostbound
)
{
    const size_t num_bh_parts = 0;
    vr_return_data return_data = InvokeVelociraptorHydro(snapnum, outputname, c, s,
        num_gravity_parts, num_hydro_parts, num_star_parts, num_bh_parts,
        swift_parts, cell_node_ids,
        numthreads, ireturngroupinfoflag, ireturnmostbound);
    return return_data;
}
vr_return_data InvokeVelociraptorHydro(const int snapnum, char* outputname,
    cosmoinfo c, siminfo s,
    const size_t num_gravity_parts, const size_t num_hydro_parts, const size_t num_star_parts, const size_t num_bh_parts,
    struct swift_vel_part *swift_parts, int *cell_node_ids,
    const int numthreads,
    const int ireturngroupinfoflag, const int ireturnmostbound,
    struct swift_vel_gas_part *swift_gas_parts,
    struct swift_vel_star_part *swift_star_parts,
    struct swift_vel_bh_part *swift_bh_parts
)
{
#ifdef USEMPI
    LOG_RANK0(info) << "VELOCIraptor/STF running with MPI. Number of mpi threads: "<< NProcs;
#else
    int ThisTask=0;
    int NProcs=1;
    Int_t Nlocal, Ntotal, Nmemlocal;
    Int_t Nmemlocalbaryon, Nlocalbaryon[1];
#endif
    int nthreads;
#ifdef USEOPENMP
    omp_set_num_threads(numthreads);
#pragma omp parallel
{
    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
}
    LOG_RANK0(info) << "VELOCIraptor/STF running with OpenMP. Number of openmp threads: " << nthreads;
#else
    nthreads=1;
#endif

    // Construct default return value
    vr_return_data return_data;
    return_data.num_gparts_in_groups = 0;
    return_data.group_info = NULL;
    return_data.num_most_bound = 0;
    return_data.most_bound_index = NULL;

    libvelociraptorOpt.outname = outputname;
    libvelociraptorOpt.snapshotvalue = HALOIDSNVAL* snapnum;
    libvelociraptorOpt.memuse_peak = 0;
    libvelociraptorOpt.memuse_ave = 0;
    libvelociraptorOpt.memuse_nsamples = 0;

    //write associated units and simulation details (which contains scale factor/time information)
    SetVelociraptorSimulationState(c, s);
    WriteSimulationInfo(libvelociraptorOpt);
    WriteUnitInfo(libvelociraptorOpt);
    InitMemUsageLog(libvelociraptorOpt);

    vector<Particle> parts;
    #ifdef GASON
    HydroProperties hydro;
    #endif
    #ifdef STARON
    StarProperties star;
    #endif
    #ifdef BHON
    BHProperties bh;
    #endif
    Particle *pbaryons;
    Int_t *pfof = NULL, *pfofall = NULL, *pfofbaryons = NULL, *numingroup = NULL, **pglist = NULL;
    Int_t *nsub = NULL, *parentgid = NULL, *uparentgid =NULL, *stype = NULL;
    Int_t nbaryons, ndark, index;
    Int_t ngroup, nhalos;
    //KDTree *tree;
    //to store information about the group
    PropData *pdata = NULL,*pdatahalos = NULL;
    double time1;

    /// Set pointer to cell node IDs
    libvelociraptorOpt.cellnodeids = cell_node_ids;

    Nlocal=Nmemlocal=num_gravity_parts;
#ifdef USEMPI
    if(NProcs>1)Nmemlocal*=(1+libvelociraptorOpt.mpipartfac); /* JSW: Not set in parameter file. */
    MPI_Allreduce(&Nlocal, &Ntotal, 1, MPI_Int_t, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allgather(&Nlocal, 1, MPI_Int_t, mpi_nlocal, 1, MPI_Int_t, MPI_COMM_WORLD);
#else
    Ntotal=Nlocal;
#endif
    parts.resize(Nmemlocal);

    LOG(info) << "Copying particle data...";

    vr::Timer copy_timer;
    ndark = num_gravity_parts - num_hydro_parts - num_star_parts - num_bh_parts;
    nbaryons = num_hydro_parts + num_star_parts + num_bh_parts;
    Nlocalbaryon[0]=nbaryons;
    Nmemlocalbaryon=Nlocalbaryon[0];

#ifdef HIGHRES
    Int_t ninterloper = 0, n_dm = 0;
#endif

    /// If we are performing a baryon search, sort the particles so that the DM particles are at the start of the array followed by the gas particles.
    // note that we explicitly convert positions from comoving to physical as swift_vel_parts is in
    if (libvelociraptorOpt.iBaryonSearch>0 && libvelociraptorOpt.partsearchtype!=PSTALL) {
        size_t dmOffset = 0, baryonOffset = 0, gasOffset = 0, starOffset = 0, bhOffset = 0, otherparttype = 0;
        pbaryons=&(parts.data()[ndark]);
        LOG(info) << "There are " << nbaryons << " gas particles and " << ndark << " DM particles";
        for(auto i=0; i<Nlocal; i++)
        {
            for (auto j=0;j<3;j++) swift_parts[i].x[j]*=libvelociraptorOpt.a;
            if(swift_parts[i].type == DARKTYPE) {
                parts[dmOffset] = Particle(swift_parts[i]);
#ifdef HIGHRES
                if (dmOffset==0) libvelociraptorOpt.zoomlowmassdm = parts[dmOffset].GetMass();
#endif
                dmOffset++;
            }
#ifdef HIGHRES
            else if(swift_parts[i].type == DARK2TYPE) {
                parts[dmOffset] = Particle(swift_parts[i]);
                parts[dmOffset].SetType(DARK2TYPE);
                ninterloper++;
                dmOffset++;
            }
#endif
            else {
                if(swift_parts[i].type == GASTYPE) {
                    pbaryons[baryonOffset] = Particle(swift_parts[i]);
                    baryonOffset++;
                    gasOffset++;
                }
                else if(swift_parts[i].type == STARTYPE) {
                    pbaryons[baryonOffset] = Particle(swift_parts[i]);
                    baryonOffset++;
                    starOffset++;
                }
                else if(swift_parts[i].type == BHTYPE) {
                    pbaryons[baryonOffset] = Particle(swift_parts[i]);
                    baryonOffset++;
                    bhOffset++;
                }
                else {
                    LOG(warning) << "Unknown particle type found: index=" << i
                                 << " type=" << swift_parts[i].type
                                 << " while treating baryons differently. Exiting...";
                    return return_data;
                }
            }
        }
    }
    else {
        for(auto i=0; i<Nlocal; i++) {
            if (CheckSwiftPartType(swift_parts[i].type)){
                LOG(warning) << "Unknown particle type found: index=" << i
                             << " type=" << swift_parts[i].type
                             << " when loading particles. Exiting...";
                return return_data;
            }
            for (auto j=0;j<3;j++) swift_parts[i].x[j]*=libvelociraptorOpt.a;
            parts[i] = Particle(swift_parts[i]);
#ifdef HIGHRES
            if (swift_parts[i].type == DARKTYPE) libvelociraptorOpt.zoomlowmassdm = parts[i].GetMass();
            else if (swift_parts[i].type == DARK2TYPE) {
                parts[i].SetType(DARK2TYPE);
                ninterloper++;
            }
#endif
        }
    }
    //if extra information has been passed then store it
#ifdef GASON
    if (swift_gas_parts != NULL)
    {
        for (auto i=0; i<num_hydro_parts; i++)
        {
            index = swift_gas_parts[i].index;
            parts[index].SetHydroProperties(hydro);
        }
        free(swift_gas_parts);
    }
#endif
#ifdef STARON
    if (swift_star_parts != NULL)
    {
        for (auto i=0; i<num_star_parts; i++)
        {
            index = swift_star_parts[i].index;
            parts[index].SetStarProperties(star);
        }
        free(swift_star_parts);
    }
#endif
#ifdef BHON
    if (swift_bh_parts != NULL)
    {
        for (auto i=0; i<num_bh_parts; i++)
        {
            index = swift_bh_parts[i].index;
            parts[index].SetBHProperties(bh);
        }
        free(swift_bh_parts);
    }
#endif

    //lets free the memory of swift_parts
    free(swift_parts);

    LOG(info) << "Finished copying particle data";
#ifdef HIGHRES
    LOG(info) << "Zoom simulation where there are " << ninterloper << " low resolution interloper particles";
#endif
    LOG(info) << "It took " << copy_timer << " to copy " << Nlocal << " particles from SWIFT to a local format (out of " << Ntotal << ')';
    LOG(info) << "There are " << Nlocal << " particles and have allocated enough memory for "
              << Nmemlocal << " requiring "<< vr::memory_amount(Nmemlocal * sizeof(Particle));
    if (libvelociraptorOpt.iBaryonSearch > 0)
        LOG(info) << "There are " << Nlocalbaryon[0] << " baryon particles "
                  << "and have allocated enough memory for " << Nmemlocalbaryon
                  << " requiring " << vr::memory_amount(Nmemlocalbaryon*sizeof(Particle));
    LOG(info) << "Will also require additional memory for FOF algorithms and substructure search. "
              << "Largest mem needed for preliminary FOF search. "
              << "Rough estimate is " << vr::memory_amount(Nlocal * (sizeof(Int_tree_t) * 8));

    MEMORY_USAGE_REPORT(info, libvelociraptorOpt);

    //
    // Perform FOF search.
    //
    vr::Timer fof_timer;
    pfof=SearchFullSet(libvelociraptorOpt,Nlocal,parts,ngroup);
    LOG(info) << "FOF search with " << Nlocal << " particles and "
              << nthreads << " threads took " << fof_timer;
    nhalos=ngroup;
    //if caculating inclusive halo masses, then for simplicity, I assume halo id order NOT rearranged!
    //this is not necessarily true if baryons are searched for separately.
    if (libvelociraptorOpt.iInclusiveHalo>0 && libvelociraptorOpt.iInclusiveHalo<3) {
        pdatahalos=new PropData[nhalos+1];
        Int_t *numinhalos=BuildNumInGroup(Nlocal, nhalos, pfof);
        Int_t *sortvalhalos=new Int_t[Nlocal];
        Int_t *originalID=new Int_t[Nlocal];
        for (Int_t i=0;i<Nlocal;i++) {sortvalhalos[i]=pfof[i]*(pfof[i]>0)+Nlocal*(pfof[i]==0);originalID[i]=parts[i].GetID();parts[i].SetID(i);}
        Int_t *noffsethalos=BuildNoffset(Nlocal, parts.data(), nhalos, numinhalos, sortvalhalos);
        GetInclusiveMasses(libvelociraptorOpt, Nlocal, parts.data(), nhalos, pfof, numinhalos, pdatahalos, noffsethalos);
        qsort(parts.data(),Nlocal,sizeof(Particle),IDCompare);
        delete[] numinhalos;
        delete[] sortvalhalos;
        delete[] noffsethalos;
        for (Int_t i=0;i<Nlocal;i++) parts[i].SetID(originalID[i]);
        delete[] originalID;
    }

    //
    // Substructure search.
    //
    if (libvelociraptorOpt.iSubSearch) {
        vr::Timer timer;
        LOG(info) << "Searching subset";
        //if groups have been found (and localized to single MPI thread) then proceed to search for subsubstructures
        SearchSubSub(libvelociraptorOpt, Nlocal, parts, pfof,ngroup,nhalos,pdatahalos);
        LOG(info) << "Search for substructures " << Nlocal << " with " << nthreads
                  << " threads finished in " << timer;
    }
    pdata=new PropData[ngroup+1];
    //if inclusive halo mass required
    if (libvelociraptorOpt.iInclusiveHalo>0 && libvelociraptorOpt.iInclusiveHalo<3 && ngroup>0) {
        CopyMasses(libvelociraptorOpt,nhalos,pdatahalos,pdata);
        delete[] pdatahalos;
    }

    //
    // Search for baryons
    //

    //if only searching initially for dark matter groups, once found, search for associated baryonic structures if requried
    if (libvelociraptorOpt.iBaryonSearch>0) {
        vr::Timer timer;
        if (libvelociraptorOpt.partsearchtype==PSTDARK) {
            pfofall=SearchBaryons(libvelociraptorOpt, nbaryons, pbaryons, Nlocal, parts, pfof, ngroup,nhalos,libvelociraptorOpt.iseparatefiles,libvelociraptorOpt.iInclusiveHalo,pdata);
            pfofbaryons=&pfofall[Nlocal];
        }
        //if FOF search overall particle types then running sub search over just dm and need to associate baryons to just dm particles must determine number of baryons, sort list, run search, etc
        //but only need to run search if substructure has been searched
        else if (libvelociraptorOpt.iSubSearch==1){
            nbaryons=0;
            ndark=0;
            for (Int_t i=0;i<Nlocal;i++) {
                if (parts[i].GetType()==DARKTYPE)ndark++;
                else nbaryons++;
            }
            pbaryons=NULL;
            SearchBaryons(libvelociraptorOpt, nbaryons, pbaryons, ndark, parts, pfof, ngroup,nhalos,libvelociraptorOpt.iseparatefiles,libvelociraptorOpt.iInclusiveHalo,pdata);
        }
        LOG(info) << "Baryon search with " << nthreads << " threads finished in " << timer;
    }

    //get mpi local hierarchy
    nsub=new Int_t[ngroup+1];
    parentgid=new Int_t[ngroup+1];
    uparentgid=new Int_t[ngroup+1];
    stype=new Int_t[ngroup+1];
    Int_t nhierarchy=GetHierarchy(libvelociraptorOpt,ngroup,nsub,parentgid,uparentgid,stype);
    CopyHierarchy(libvelociraptorOpt,pdata,ngroup,nsub,parentgid,uparentgid,stype);

    //calculate data and output
    numingroup=BuildNumInGroup(Nlocal, ngroup, pfof);
    pglist=SortAccordingtoBindingEnergy(libvelociraptorOpt,Nlocal,parts.data(),ngroup,pfof,numingroup,pdata);//alters pglist so most bound particles first
    vr::Timer write_timer;
    WriteProperties(libvelociraptorOpt,ngroup,pdata);
    WriteGroupCatalog(libvelociraptorOpt, ngroup, numingroup, pglist, parts);
    //if baryons have been searched output related gas baryon catalogue
    if (libvelociraptorOpt.iBaryonSearch>0 || libvelociraptorOpt.partsearchtype==PSTALL){
      WriteGroupPartType(libvelociraptorOpt, ngroup, numingroup, pglist, parts);
    }
#ifdef EXTENDEDHALOOUTPUT
    if (opt.iExtendedOutput) WriteExtendedOutput (libvelociraptorOpt, ngroup, nbodies, pdata, parts, pfof);
#endif
    delete[] pfof;
    //if returning to swift as swift is writing a snapshot, then write for the groups where the particles are found in a file
    //assuming that the swift task and swift index can be used to determine where a particle will be written.
    if (ireturngroupinfoflag != 1 ) WriteSwiftExtendedOutput (libvelociraptorOpt, ngroup, numingroup, pglist, parts);
    LOG(info) << "Wrote all data in " << write_timer;

    // Find offset to first group on each MPI rank
    Int_t ngoffset=0,ngtot=0;
    ngtot=ngroup;
#ifdef USEMPI
    if (NProcs > 1) {
        ngtot = 0;
        for (auto j=0;j<ThisTask;j++)ngoffset+=mpi_ngroups[j];
        for (auto j=0;j<NProcs;j++)ngtot+=mpi_ngroups[j];
    }
#endif

    // Compute group membership for each particle and store in PID
    for (auto i=0; i<Nlocal; i++)parts[i].SetPID(0);
    for (auto i=1; i<=ngroup; i++) {
      for (auto j=0;j<numingroup[i];j++) {
        parts[pglist[i][j]].SetPID(i+ngoffset+libvelociraptorOpt.snapshotvalue);
      }
    }

    delete[] pdata;
    delete[] nsub;
    delete[] uparentgid;
    delete[] parentgid;
    delete[] stype;
    delete psldata;

    MEMORY_USAGE_REPORT(info, libvelociraptorOpt);

    // Make array of most bound particle indexes if we have groups and it was requested
    Int_t num_most_bound = 0;
    int *most_bound_index = NULL;
    if(ireturnmostbound && ngtot > 0) {
      
      // Make a vector containing just the most bound particle from each non-zero size group
      vector<Particle> most_bound_parts;
      most_bound_parts.reserve(ngroup);
      for (auto i=1; i<=ngroup; i++) {
        if(numingroup[i] > 0)most_bound_parts.push_back(parts[pglist[i][0]]);
      }
      num_most_bound = most_bound_parts.size();

      // Send each most bound particle to the MPI rank it was on in Swift
#ifdef USEMPI
      if (NProcs > 1) {
        for (auto i=0;i<num_most_bound; i++) most_bound_parts[i].SetID((most_bound_parts[i].GetSwiftTask()!=ThisTask));
        //now sort items according to whether local swift task
        qsort(most_bound_parts.data(), num_most_bound, sizeof(Particle), IDCompare);
        //communicate information
        MPISwiftExchange(most_bound_parts);
        num_most_bound = most_bound_parts.size();
      }
#endif
      // Make an array with the Swift indexes of the most bound particles.
      // This will be returned to Swift so it has to be allocated with malloc().
      most_bound_index = (int *) malloc(sizeof(int)*num_most_bound);
      for (auto i=0;i<num_most_bound; i++) {
        most_bound_index[i] = most_bound_parts[i].GetSwiftIndex();
      }
      most_bound_parts.clear();
    }
    // Add most bound particles to struct to return
    return_data.num_most_bound = num_most_bound;
    return_data.most_bound_index = most_bound_index;

    for (Int_t i=1;i<=ngroup;i++) delete[] pglist[i];
    delete[] pglist;
    delete[] numingroup;

    // Compute group_info array if we have groups and it was requested
    groupinfo *group_info = NULL;
    Int_t nig=0;
    if(ireturngroupinfoflag && ngtot > 0) {
      //first sort so all particles in groups first
      LOG(info) << "VELOCIraptor sorting info to return group ids to swift";
#ifdef USEMPI
      if (NProcs > 1) {
        for (auto i=0;i<Nlocal; i++) parts[i].SetID((parts[i].GetSwiftTask()!=ThisTask));
        //now sort items according to whether local swift task
        qsort(parts.data(), Nlocal, sizeof(Particle), IDCompare);
        //communicate information
        MPISwiftExchange(parts);
        Nlocal = parts.size(); 
      }
#endif
      qsort(parts.data(), Nlocal, sizeof(Particle), PIDCompare);
      nig=0;
      Int_t istart=0;
      for (auto i=0;i<Nlocal;i++) if (parts[i].GetPID()>0) {nig=Nlocal-i;istart=i;break;}
      if (nig>0)
        {
          group_info = (groupinfo *) malloc(nig*sizeof(struct groupinfo)); // Will be freed by Swift
          for (auto i=istart;i<Nlocal;i++) {
            group_info[i-istart].index=parts[i].GetSwiftIndex();
            group_info[i-istart].groupid=parts[i].GetPID();
          }
        }
    }
    // Add groupinfo to struct to return
    return_data.num_gparts_in_groups=nig;
    return_data.group_info = group_info;

    //free mem associate with mpi cell node ides
    libvelociraptorOpt.cellnodeids = NULL;
    libvelociraptorOpt.cellloc = NULL;
    free(s.cellloc);
    free(cell_node_ids);
    parts.clear();

    LOG(info) << "VELOCIraptor returning.";
    return return_data;
}


///\todo interface with swift is comoving positions, period, peculiar velocities and physical self-energy
vr_return_data InvokeVelociraptorExtra(const int iextra, const int snapnum, char* outputname,
    cosmoinfo c, siminfo s,
    const size_t num_gravity_parts, const size_t num_hydro_parts, const size_t num_star_parts,
    struct swift_vel_part *swift_parts, int *cell_node_ids,
    const int numthreads,
    const int ireturngroupinfoflag, const int ireturnmostbound)
{
    //store backup of the libvelociraptorOpt
    libvelociraptorOptbackup = libvelociraptorOpt;
    libvelociraptorOpt = libvelociraptorOptextra[iextra];
    vr_return_data return_data = InvokeVelociraptor(snapnum, outputname, c, s,
        num_gravity_parts, num_hydro_parts, num_star_parts,
        swift_parts, cell_node_ids,
        numthreads,
        ireturngroupinfoflag, ireturnmostbound);
    libvelociraptorOpt = libvelociraptorOptbackup;
    return return_data;
}

vr_return_data InvokeVelociraptorHydroExtra(const int iextra, const int snapnum, char* outputname,
    cosmoinfo c, siminfo s,
    const size_t num_gravity_parts, const size_t num_hydro_parts, const size_t num_star_parts, const size_t num_bh_parts,
    struct swift_vel_part *swift_parts, int *cell_node_ids,
    const int numthreads,
    const int ireturngroupinfoflag, const int ireturnmostbound,
    struct swift_vel_gas_part *swift_gas_parts,
    struct swift_vel_star_part *swift_star_parts,
    struct swift_vel_bh_part *swift_bh_parts
)
{
    //store backup of the libvelociraptorOpt
    libvelociraptorOptbackup = libvelociraptorOpt;
    libvelociraptorOpt = libvelociraptorOptextra[iextra];
    vr_return_data return_data = InvokeVelociraptorHydro(snapnum, outputname, c, s,
        num_gravity_parts, num_hydro_parts, num_star_parts, num_bh_parts,
        swift_parts, cell_node_ids, numthreads, ireturngroupinfoflag, ireturnmostbound,
        swift_gas_parts, swift_star_parts, swift_bh_parts);
    libvelociraptorOpt = libvelociraptorOptbackup;
    return return_data;
}

void CheckSwiftTasks(string message, const Int_t n, Particle *p){
    #ifndef USEMPI
    int ThisTask=0,NProcs=1;
    #endif
    int numbad = 0, numnonlocal=0;
    LOG(info) << "Checking swift tasks:" << message;
    for (auto i=0;i<n;i++) {
        int task=p[i].GetSwiftTask();
        int index = p[i].GetSwiftIndex();
        if (task<0 || task>NProcs) {
            LOG(error) << " has odd swift particle with nonsensical swift task (cur index, swift index, swift task)=("
                       << i << "," << index << "," << task << ")";
            numbad++;
        }
        if (task!=ThisTask) numnonlocal++;
    }
    LOG(info) << numnonlocal << " Swift particles found";
    if (numbad>0) {
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,9);
#else
        exit(9);
#endif
    }
}
#endif
