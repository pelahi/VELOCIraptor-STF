/*! \file main.cxx
 *  \brief Main program

    This code initialize MPI (if necessary) and reads parameters from the user input (see \ref ui.cxx)
    loads particle data from a N-body output (see \ref io.cxx) into an array of \ref NBody::Particle,
    analyzes the system, searches for structures then substructures (see \ref search.cxx),
    outputs the group ids (in two different fashions see \ref io.cxx) and can also analyze the structures
    (see \ref substructureproperties.cxx and \ref io.cxx).

    \todo remove array mpi_idlist and mpi_indexlist as these arrays are unnecessary
    \todo alter unbinding/sortbybindingenergy calls since seems a waste of cpu cycles
*/

#include "stf.h"

using namespace std;
using namespace Math;
using namespace NBody;
using namespace velociraptor;

int main(int argc,char **argv)
{
#ifdef SWIFTINTERFACE
    cout<<"Built with SWIFT interface enabled when running standalone VELOCIraptor. Should only be enabled when running VELOCIraptor as a library from SWIFT. Exiting..."<<endl;
    exit(0);
#endif
#ifdef USEMPI
    //start MPI
#ifdef USEOPENMP
    //if using hybrid then need to check that threads are available and use the correct initialization
    //Each thread will call MPI routines, but these calls will be coordinated to occur only one at a time within a process.
    int required=MPI_THREAD_FUNNELED;  // Required level of MPI threading support
    int provided; // Provided level of MPI threading support
    MPI_Init_thread(&argc, &argv, required, &provided);
#else
    MPI_Init(&argc,&argv);
#endif

    //find out how big the SPMD world is
    MPI_Comm_size(MPI_COMM_WORLD,&NProcs);
    //and this processes' rank is
    MPI_Comm_rank(MPI_COMM_WORLD,&ThisTask);

    if (ThisTask == 0) cout<<"Running VELOCIraptor "<<git_sha1()<<endl;
#ifdef USEOPENMP
    // Check the threading support level
    if (provided < required)
    {
        // Insufficient support, degrade to 1 thread and warn the user
        if (ThisTask == 0) cout << "Warning: This MPI implementation provides insufficient threading support. Required was " <<required<<" but provided was "<<provided<<endl;
        omp_set_num_threads(1);
        MPI_Finalize();
        exit(9);
    }
#endif

#else
    int ThisTask=0,NProcs=1;
    Int_t Nlocal,Ntotal;
    Int_t Ntotalbaryon[NBARYONTYPES], Nlocalbaryon[NBARYONTYPES];
#endif
    int nthreads;
#ifdef USEOPENMP
    nthreads = omp_get_max_threads();
#else
    nthreads=1;
#endif

#ifdef USEMPI
    if (ThisTask==0) cout<<"VELOCIraptor/STF running with MPI. Number of mpi threads: "<<NProcs<<endl;
#endif
#ifdef USEOPENMP
    if (ThisTask==0) cout<<"VELOCIraptor/STF running with OpenMP. Number of openmp threads: "<<nthreads<<endl;
    if (ThisTask==0) cout<<"VELOCIraptor/STF running with OpenMP version "<< _OPENMP << endl;
#endif

    gsl_set_error_handler_off();

    Options opt;
    //get arguments
    GetArgs(argc, argv, opt);
    cout.precision(10);

#ifdef USEMPI
    MPIInitWriteComm();
#ifdef USEADIOS
    //init adios
    adios_init_noxml(MPI_COMM_WORLD);
    //specify the buffer size (use the tot particle buffer size in bytes and convert to MB)
    adios_set_max_buffer_size(opt.mpiparticletotbufsize/1024/1024);
#endif
#endif

    InitMemUsageLog(opt);

    //variables
    //number of particles, (also number of baryons if use dm+baryon search)
    //to store (point to) particle data
    Int_t nbodies,nbaryons,ndark;
    vector<Particle> Part;
    Particle *Pbaryons;
    KDTree *tree;

    //number in subset, number of grids used if iSingleHalo==0;
    Int_t nsubset, ngrid;
    //number of groups, temp num groups, and number of halos
    Int_t ngroup, ng, nhalos;

    //to store group value (pfof), and also arrays to parse particles
    //vector<Int_t> pfof, pfofbaryons;
    Int_t *pfof, *pfofbaryons;;
    Int_t *numingroup,**pglist;
    Int_t *pfofall;
    //to store information about the group
    PropData *pdata=NULL,*pdatahalos=NULL;

    //to store time and output time taken
    auto tottime=MyGetTime();

    Coordinate cm,cmvel;
    Double_t Mtot;
    char fname1[1000],fname4[1000];

#ifdef USEMPI
    mpi_nlocal=new Int_t[NProcs];
    mpi_domain=new MPI_Domain[NProcs];
    mpi_nsend=new Int_t[NProcs*NProcs];
    mpi_ngroups=new Int_t[NProcs];
    mpi_nhalos=new Int_t[NProcs];
    //store MinSize as when using mpi prior to stitching use min of 2;
    MinNumMPI=2;
    //if single halo, use minsize to initialize the old minimum number
    //else use the halominsize since if mpi and not single halo, halos localized to mpi domain for substructure search
    if (opt.iSingleHalo) MinNumOld=opt.MinSize;
    else MinNumOld=opt.HaloMinSize;
#endif

    //read particle information and allocate memory
    auto time1=MyGetTime();
    //for MPI determine total number of particles AND the number of particles assigned to each processor
    if (ThisTask==0) {
        cout<<"Read header ... "<<endl;
        nbodies=ReadHeader(opt);
        if (opt.iBaryonSearch>0) {
            for (int i=0;i<NBARYONTYPES;i++) Ntotalbaryon[i]=Nlocalbaryon[i]=0;
            nbaryons=0;
            int pstemp=opt.partsearchtype;
            opt.partsearchtype=PSTGAS;
            nbaryons+=ReadHeader(opt);
            if (opt.iusestarparticles) {
                opt.partsearchtype=PSTSTAR;
                nbaryons+=ReadHeader(opt);
            }
            if (opt.iusesinkparticles) {
                opt.partsearchtype=PSTBH;
                nbaryons+=ReadHeader(opt);
            }
            opt.partsearchtype=pstemp;
        }
        else nbaryons=0;
    }

#ifdef USEMPI
    MPI_Bcast(&nbodies,1, MPI_Int_t,0,MPI_COMM_WORLD);
    //if MPI, possible desired particle types not present in file so update used particle types
    MPIUpdateUseParticleTypes(opt);
    if (opt.iBaryonSearch>0) MPI_Bcast(&nbaryons,1, MPI_Int_t,0,MPI_COMM_WORLD);
    //initial estimate need for memory allocation assuming that work balance is not greatly off
#endif
    if (ThisTask==0) {
        cout<<"There are "<<nbodies<<" particles in total that require "<<nbodies*sizeof(Particle)/1024./1024./1024.<<"GB of memory "<<endl;
        if (opt.iBaryonSearch>0) cout<<"There are "<<nbaryons<<" baryon particles in total that require "<<nbaryons*sizeof(Particle)/1024./1024./1024.<<"GB of memory "<<endl;
    }

    //note that for nonmpi particle array is a contiguous block of memory regardless of whether a separate baryon search is required
#ifndef USEMPI
    Nlocal=nbodies;
    if (opt.iBaryonSearch>0 && opt.partsearchtype!=PSTALL) {
        Part.resize(nbodies+nbaryons);
        Pbaryons=&(Part.data()[nbodies]);
        Nlocalbaryon[0]=nbaryons;
    }
    else {
        Part.resize(nbodies);
        Pbaryons=NULL;
        nbaryons=0;
    }
#else
    //for mpi however, it is not possible to have a simple contiguous block of memory IFF a separate baryon search is required.
    //for the simple reason that the local number of particles changes to ensure large fof groups are local to an mpi domain
    //however, when reading data, it is much simplier to have a contiguous block of memory, sort that memory (if necessary)
    //and then split afterwards the dm particles and the baryons
    if (NProcs==1) {Nlocal=Nmemlocal=nbodies;NExport=NImport=1;}
    else {
#ifdef MPIREDUCEMEM
        //if allocating reasonable amounts of memory, use MPIREDUCEMEM
        //this determines number of particles in the mpi domains
        MPINumInDomain(opt);
        cout<<ThisTask<<" There are "<<Nlocal<<" particles and have allocated enough memory for "<<Nmemlocal<<" requiring "<<Nmemlocal*sizeof(Particle)/1024./1024./1024.<<"GB of memory "<<endl;
        if (opt.iBaryonSearch>0) cout<<ThisTask<<"There are "<<Nlocalbaryon[0]<<" baryon particles and have allocated enough memory for "<<Nmemlocalbaryon<<" requiring "<<Nmemlocalbaryon*sizeof(Particle)/1024./1024./1024.<<"GB of memory "<<endl;
#else
        //otherwise just base on total number of particles * some factor and initialise the domains
        MPIDomainExtent(opt);
        MPIDomainDecomposition(opt);
        Nlocal=nbodies/NProcs*MPIProcFac;
        Nmemlocal=Nlocal;
        Nlocalbaryon[0]=nbaryons/NProcs*MPIProcFac;
        Nmemlocalbaryon=Nlocalbaryon[0];
        NExport=NImport=Nlocal*MPIExportFac;
        cout<<ThisTask<<" Have allocated enough memory for "<<Nmemlocal<<" requiring "<<Nmemlocal*sizeof(Particle)/1024./1024./1024.<<"GB of memory "<<endl;
        if (opt.iBaryonSearch>0) cout<<" Have allocated enough memory for "<<Nmemlocalbaryon<<" baryons particles requiring "<<Nmemlocalbaryon*sizeof(Particle)/1024./1024./1024.<<"GB of memory "<<endl;
#endif
    }
    cout<<ThisTask<<" will also require additional memory for FOF algorithms and substructure search. Largest mem needed for preliminary FOF search. Rough estimate is "<<Nlocal*(sizeof(Int_tree_t)*8)/1024./1024./1024.<<"GB of memory"<<endl;
    if (opt.iBaryonSearch>0 && opt.partsearchtype!=PSTALL) {
        Part.resize(Nmemlocal+Nmemlocalbaryon);
        Pbaryons=&(Part.data()[Nlocal]);
        nbaryons=Nlocalbaryon[0];
    }
    else {
        Part.resize(Nmemlocal);
        Pbaryons=NULL;
        nbaryons=0;
    }
#endif

    //now read particle data
    if (ThisTask==0) cout<<"Loading ... "<<endl;
    ReadData(opt, Part, nbodies, Pbaryons, nbaryons);
#ifdef USEMPI
    //if mpi and want separate baryon search then once particles are loaded into contigous block of memory and sorted according to type order,
    //allocate memory for baryons
    if (opt.iBaryonSearch>0 && opt.partsearchtype!=PSTALL) {
        Pbaryons=new Particle[Nmemlocalbaryon];
        nbaryons=Nlocalbaryon[0];

        for (Int_t i=0;i<Nlocalbaryon[0];i++) Pbaryons[i]=Part[i+Nlocal];
        Part.resize(Nlocal);
    }
#endif

#ifdef USEMPI
    Ntotal=nbodies;
    nbodies=Nlocal;
    NExport=NImport=Nlocal*MPIExportFac;
    mpi_period=opt.p;
    MPI_Allgather(&nbodies, 1, MPI_Int_t, mpi_nlocal, 1, MPI_Int_t, MPI_COMM_WORLD);
    MPI_Allreduce(&nbodies, &Ntotal, 1, MPI_Int_t, MPI_SUM, MPI_COMM_WORLD);
    cout<<"TIME::"<<ThisTask<<" took "<<MyElapsedTime(time1)<<" to load "<<Nlocal<<" of "<<Ntotal<<endl;
#else
    cout<<"TIME::"<<ThisTask<<" took "<<MyElapsedTime(time1)<<" to load "<<nbodies<<endl;
#endif

    //write out the configuration used by velociraptor having read in the data (as input data can contain cosmological information)
    WriteVELOCIraptorConfig(opt);
    WriteSimulationInfo(opt);
    WriteUnitInfo(opt);

    //set filenames if they have been passed
#ifdef USEMPI
    if (opt.smname!=NULL) sprintf(fname4,"%s.%d",opt.smname,ThisTask);
#else
    if (opt.smname!=NULL) sprintf(fname4,"%s",opt.smname);
#endif

    //read local velocity data or calculate it
    //(and if STRUCDEN flag or HALOONLYDEN is set then only calculate the velocity density function for objects within a structure
    //as found by SearchFullSet)
#if defined (STRUCDEN) || defined (HALOONLYDEN)
#else
    if (opt.iSubSearch==1) {
        time1=MyGetTime();
        if(FileExists(fname4)) ReadLocalVelocityDensity(opt, nbodies,Part);
        else  {
            GetVelocityDensity(opt, nbodies, Part.data());
            WriteLocalVelocityDensity(opt, nbodies,Part);
        }
        time1=MyGetTime()-time1;
        cout<<"TIME::"<<ThisTask<<" took "<<time1<<" to analyze/read local velocity density for "<<Nlocal<<" with "<<nthreads<<endl;
    }
#endif

    //here adjust Efrac to Omega_cdm/Omega_m from what it was before if baryonic search is separate
    if (opt.iBaryonSearch>0 && opt.partsearchtype!=PSTALL) opt.uinfo.Eratio*=opt.Omega_cdm/opt.Omega_m;

    //From here can either search entire particle array for "Halos" or if a single halo is loaded, then can just search for substructure
    if (!opt.iSingleHalo) {
#ifndef USEMPI
        time1=MyGetTime();
        pfof=SearchFullSet(opt,nbodies,Part,ngroup);
        nhalos=ngroup;
        cout<<"TIME:: took "<<MyElapsedTime(time1)<<" to search "<<nbodies<<" with "<<nthreads<<endl;
#else
        //nbodies=Ntotal;
        ///\todo Communication Buffer size determination and allocation. For example, eventually need something like FoFDataIn = (struct fofdata_in *) CommBuffer;
        ///At the moment just using NExport
        NExport=Nlocal*MPIExportFac;
        //Now when MPI invoked this returns pfof after local linking and linking across and also reorders groups
        //according to size and localizes the particles belong to the same group to the same mpi thread.
        //after this is called Nlocal is adjusted to the local subset where groups are localized to a given mpi thread.
        time1=MyGetTime();

        pfof=SearchFullSet(opt,Nlocal,Part,ngroup);
        cout<<"TIME::"<<ThisTask<<" took "<<MyElapsedTime(time1)<<" to search "<<Nlocal<<" with "<<nthreads<<endl;
        nbodies=Nlocal;
        nhalos=ngroup;
#endif
        //if compiled to determine inclusive halo masses, then for simplicity, I assume halo id order NOT rearranged!
        //this is not necessarily true if baryons are searched for separately.
        if (opt.iInclusiveHalo > 0 && opt.iInclusiveHalo < 3) {
            pdatahalos=new PropData[nhalos+1];
            Int_t *numinhalos=BuildNumInGroup(nbodies, nhalos, pfof);
            Int_t *sortvalhalos=new Int_t[nbodies];
            Int_t *originalID=new Int_t[nbodies];
            for (Int_t i=0;i<nbodies;i++) {sortvalhalos[i]=pfof[i]*(pfof[i]>0)+nbodies*(pfof[i]==0);originalID[i]=Part[i].GetID();Part[i].SetID(i);}
            Int_t *noffsethalos=BuildNoffset(nbodies, Part.data(), nhalos, numinhalos, sortvalhalos);
            ///here if inclusive halo flag is 3, then S0 masses are calculated after substructures are found for field objects
            ///and only calculate FOF masses. Otherwise calculate inclusive masses at this moment.
            GetInclusiveMasses(opt, nbodies, Part.data(), nhalos, pfof, numinhalos, pdatahalos, noffsethalos);
            sort(Part.begin(), Part.end(), IDCompareVec);
            delete[] numinhalos;
            delete[] sortvalhalos;
            delete[] noffsethalos;
            for (Int_t i=0;i<nbodies;i++) Part[i].SetID(originalID[i]);
            delete[] originalID;
        }
    }
    else {
        Coordinate *gvel;
        Matrix *gveldisp;
        GridCell *grid;
        ///\todo Scaling is still not MPI compatible
        if (opt.iScaleLengths) ScaleLinkingLengths(opt,nbodies,Part.data(),cm,cmvel,Mtot);
        opt.Ncell=opt.Ncellfac*nbodies;
        //build grid using leaf nodes of tree (which is guaranteed to be adaptive and have maximum number of particles in cell of tree bucket size)
        tree=InitializeTreeGrid(opt,nbodies,Part.data());
        ngrid=tree->GetNumLeafNodes();
        cout<<"Given "<<nbodies<<" particles, and max cell size of "<<opt.Ncell<<" there are "<<ngrid<<" leaf nodes or grid cells, with each node containing ~"<<nbodies/ngrid<<" particles"<<endl;
        grid=new GridCell[ngrid];
        //note that after this system is back in original order as tree has been deleted.
        FillTreeGrid(opt, nbodies, ngrid, tree, Part.data(), grid);
        //calculate cell quantities to get mean field
        gvel=GetCellVel(opt,nbodies,Part.data(),ngrid,grid);
        gveldisp=GetCellVelDisp(opt,nbodies,Part.data(),ngrid,grid,gvel);
        opt.HaloSigmaV=0;for (int j=0;j<ngrid;j++) opt.HaloSigmaV+=pow(gveldisp[j].Det(),1./3.);opt.HaloSigmaV/=(double)ngrid;

        //now that have the grid cell volume quantities and local volume density
        //can determine the logarithmic ratio between the particle velocity density and that predicted by the background velocity distribution
        GetDenVRatio(opt,nbodies,Part.data(),ngrid,grid,gvel,gveldisp);
        //WriteDenVRatio(opt,nbodies,Part);
        //and then determine how much of an outlier it is
        nsubset=GetOutliersValues(opt,nbodies,Part.data());
        //save the normalized denvratio and also determine how many particles lie above the threshold.
        //nsubset=WriteOutlierValues(opt, nbodies,Part);
        //Now check if any particles are above the threshold
        if (nsubset==0) {
            cout<<"no particles found above threshold of "<<opt.ellthreshold<<endl;
            cout<<"Exiting"<<endl;
            return 0;
        }
        else cout<<nsubset<< " above threshold of "<<opt.ellthreshold<<" to be searched"<<endl;
#ifndef USEMPI
        pfof=SearchSubset(opt,nbodies,nbodies,Part.data(),ngroup);
#else
        //nbodies=Ntotal;
        ///\todo Communication Buffer size determination and allocation. For example, eventually need something like FoFDataIn = (struct fofdata_in *) CommBuffer;
        ///At the moment just using NExport
        NExport=Nlocal*MPIExportFac;
        mpi_foftask=MPISetTaskID(nbodies);

        //Now when MPI invoked this returns pfof after local linking and linking across and also reorders groups
        //according to size and localizes the particles belong to the same group to the same mpi thread.
        //after this is called Nlocal is adjusted to the local subset where groups are localized to a given mpi thread.
        pfof=SearchSubset(opt,Nlocal,Nlocal,Part.data(),ngroup);
        nbodies=Nlocal;
#endif
    }
    if (opt.iSubSearch) {
        cout<<"Searching subset"<<endl;
        time1=MyGetTime();
        //if groups have been found (and localized to single MPI thread) then proceed to search for subsubstructures
        SearchSubSub(opt, nbodies, Part, pfof,ngroup,nhalos, pdatahalos);
        cout<<"TIME::"<<ThisTask<<" took "<<MyElapsedTime(time1)<<" to search for substructures "<<Nlocal<<" with "<<nthreads<<endl;
    }
    pdata=new PropData[ngroup+1];
    //if inclusive halo mass required
    if (opt.iInclusiveHalo > 0 && opt.iInclusiveHalo < 3 && ngroup>0) {
        CopyMasses(opt,nhalos,pdatahalos,pdata);
        delete[] pdatahalos;
    }

    //if only searching initially for dark matter groups, once found, search for associated baryonic structures if requried
    if (opt.iBaryonSearch>0) {
        time1=MyGetTime();
        if (opt.partsearchtype==PSTDARK) {
            pfofall=SearchBaryons(opt, nbaryons, Pbaryons, nbodies, Part, pfof, ngroup,nhalos,opt.iseparatefiles,opt.iInclusiveHalo,pdata);
            pfofbaryons=&pfofall[nbodies];
        }
        //if FOF search overall particle types then running sub search over just dm and need to associate baryons to just dm particles must determine number of baryons, sort list, run search, etc
        //but only need to run search if substructure has been searched
        else if (opt.iSubSearch==1){
            nbaryons=0;
            ndark=0;
            for (Int_t i=0;i<nbodies;i++) {
                if (Part[i].GetType()==DARKTYPE)ndark++;
                else nbaryons++;
            }
            Pbaryons=NULL;
            SearchBaryons(opt, nbaryons, Pbaryons, ndark, Part, pfof, ngroup,nhalos,opt.iseparatefiles,opt.iInclusiveHalo,pdata);
        }
        cout<<"TIME::"<<ThisTask<<" took "<<MyElapsedTime(time1)<<" to search baryons  with "<<nthreads<<endl;
    }

    //get mpi local hierarchy
    Int_t *nsub,*parentgid, *uparentgid,*stype;
    nsub=new Int_t[ngroup+1];
    parentgid=new Int_t[ngroup+1];
    uparentgid=new Int_t[ngroup+1];
    stype=new Int_t[ngroup+1];
    Int_t nhierarchy=GetHierarchy(opt,ngroup,nsub,parentgid,uparentgid,stype);
    CopyHierarchy(opt,pdata,ngroup,nsub,parentgid,uparentgid,stype);

    //if a separate baryon search has been run, now just place all particles together
    if (opt.iBaryonSearch>0 && opt.partsearchtype!=PSTALL) {
        delete[] pfof;
        pfof=&pfofall[0];
        nbodies+=nbaryons;
        Nlocal=nbodies;
    }

    //output results
    //if want to ignore any information regard particles themselves as particle PIDS are meaningless
    //which might be useful for runs where not interested in tracking just halo catalogues (save for
    //approximate methods like PICOLA. Here it writes desired output and exits
    if(opt.inoidoutput){
        numingroup=BuildNumInGroup(Nlocal, ngroup, pfof);
        CalculateHaloProperties(opt,Nlocal,Part.data(),ngroup,pfof,numingroup,pdata);
        WriteProperties(opt,ngroup,pdata);
        if (opt.iprofilecalc) WriteProfiles(opt, ngroup, pdata);
        delete[] numingroup;
        delete[] pdata;

        //get memory useage
        cout<<ThisTask<<" : finished running VR "<<endl;
        GetMemUsage(opt, __func__+string("--line--")+to_string(__LINE__), (opt.iverbose>=0));

#ifdef USEMPI
#ifdef USEADIOS
        adios_finalize(ThisTask);
#endif
        MPI_Finalize();
#endif
        return 0;
    }

    //if want a simple tipsy still array listing particles group ids in input order
    if(opt.iwritefof) {
#ifdef USEMPI
        if (ThisTask==0) {
            mpi_pfof=new Int_t[Ntotal];
            //since pfof is a local subset, not all pfof values have been set, thus initialize them to zero.
            for (Int_t i=0;i<Ntotal;i++) mpi_pfof[i]=0;
        }
        MPICollectFOF(Ntotal, pfof);
        if (ThisTask==0) WriteFOF(opt,Ntotal,mpi_pfof);
#else
        WriteFOF(opt,nbodies,pfof);
#endif
    }
    numingroup=BuildNumInGroup(Nlocal, ngroup, pfof);

    //if separate files explicitly save halos, associated baryons, and subhalos separately
    if (opt.iseparatefiles) {
        if (nhalos>0) {
            pglist=SortAccordingtoBindingEnergy(opt,Nlocal,Part.data(),nhalos,pfof,numingroup,pdata);//alters pglist so most bound particles first
            WriteProperties(opt,nhalos,pdata);
            WriteGroupCatalog(opt, nhalos, numingroup, pglist, Part,ngroup-nhalos);
            //if baryons have been searched output related gas baryon catalogue
            if (opt.partsearchtype==PSTALL){
                WriteGroupPartType(opt, nhalos, numingroup, pglist, Part);
            }
            WriteHierarchy(opt,ngroup,nhierarchy,psldata->nsinlevel,nsub,parentgid,stype);
            for (Int_t i=1;i<=nhalos;i++) delete[] pglist[i];
            delete[] pglist;
        }
        else {
#ifdef USEMPI
            //if calculating inclusive masses at end, must call SortAccordingtoBindingEnergy if
            //MPI as domain, despite having no groups might need to exchange particles
            if (opt.iInclusiveHalo==3) SortAccordingtoBindingEnergy(opt,Nlocal,Part.data(),nhalos,pfof,numingroup,pdata);
#endif
            WriteGroupCatalog(opt,nhalos,numingroup,NULL,Part);
            WriteHierarchy(opt,nhalos,nhierarchy,psldata->nsinlevel,nsub,parentgid,stype);
            if (opt.partsearchtype==PSTALL){
                WriteGroupPartType(opt, nhalos, numingroup, NULL, Part);
            }
        }
    }
    Int_t indexii=0;
    ng=ngroup;
    //if separate files, alter offsets
    if (opt.iseparatefiles) {
        sprintf(fname1,"%s.sublevels",opt.outname);
        sprintf(opt.outname,"%s",fname1);
        //alter index point to just output sublevels (assumes no reordering and assumes no change in nhalos as a result of unbinding in SubSubSearch)
        indexii=nhalos;
        ng=ngroup-nhalos;
    }

    if (ng>0) {
        pglist=SortAccordingtoBindingEnergy(opt,Nlocal,Part.data(),ng,pfof,&numingroup[indexii],&pdata[indexii],indexii);//alters pglist so most bound particles first
        WriteProperties(opt,ng,&pdata[indexii]);
        WriteGroupCatalog(opt, ng, &numingroup[indexii], pglist, Part);
        if (opt.iseparatefiles) WriteHierarchy(opt,ngroup,nhierarchy,psldata->nsinlevel,nsub,parentgid,stype,1);
        else WriteHierarchy(opt,ngroup,nhierarchy,psldata->nsinlevel,nsub,parentgid,stype,-1);
        if (opt.partsearchtype==PSTALL){
            WriteGroupPartType(opt, ng, &numingroup[indexii], pglist, Part);
        }
        for (Int_t i=1;i<=ng;i++) delete[] pglist[i];
        delete[] pglist;
    }
    else {
#ifdef USEMPI
        //if calculating inclusive masses at end, must call SortAccordingtoBindingEnergy if
        //MPI as domain, despite having no groups might need to exchange particles
        if (opt.iInclusiveHalo==3) SortAccordingtoBindingEnergy(opt,Nlocal,Part.data(),ng,pfof,numingroup,pdata);
#endif
        WriteProperties(opt,ng,NULL);
        WriteGroupCatalog(opt,ng,&numingroup[indexii],NULL,Part);
        if (opt.iseparatefiles) WriteHierarchy(opt,ngroup,nhierarchy,psldata->nsinlevel,nsub,parentgid,stype,1);
        else WriteHierarchy(opt,ngroup,nhierarchy,psldata->nsinlevel,nsub,parentgid,stype,-1);
        if (opt.partsearchtype==PSTALL){
            WriteGroupPartType(opt, ng, &numingroup[indexii], NULL, Part);
        }
    }

    if (opt.iprofilecalc) WriteProfiles(opt, ngroup, pdata);

#ifdef EXTENDEDHALOOUTPUT
    if (opt.iExtendedOutput) WriteExtendedOutput (opt, ngroup, Nlocal, pdata, Part, pfof);
#endif

    delete[] pfof;
    delete[] numingroup;
    delete[] pdata;
    delete psldata;


    delete[] nsub;
    delete[] parentgid;
    delete[] uparentgid;
    delete[] stype;

    cout<<"TIME::"<<ThisTask<<" took "<<MyElapsedTime(tottime)<<" in all"<<endl;

    //get memory useage
    cout<<ThisTask<<" : finished running VR "<<endl;
    GetMemUsage(opt, __func__+string("--line--")+to_string(__LINE__), (opt.iverbose>=0));

#ifdef USEMPI
#ifdef USEADIOS
    adios_finalize(ThisTask);
#endif
    MPI_Finalize();
#endif

    return 0;
}
