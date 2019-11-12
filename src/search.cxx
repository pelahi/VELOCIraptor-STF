/*! \file search.cxx
 *  \brief this file contains routines that search particle list using FOF routines
 */

//--  Suboutines that search particle list

#include "stf.h"

#include "swiftinterface.h"

/// \name Searches full system
//@{

/*!
    Search full system without finding outliers first
    Here search simulation for FOF haloes (either 3d or 6d). Note for 6D search, first 3d FOF halos are found, then velocity scale is set by largest 3DFOF
    halo, and 6DFOF search is done.
    If all particles searched but also want to treat baryons differently, then a FOF search is run linking DM particles and DM+baryons but not the other way
    around (no baryon+DM). That is DM particles are basis for generating links but NOT gas/star/bh particles
    Also start of implementation to keep the 3DFOF envelopes as separate structures.
    \todo 3DFOF envelop kept as separate structures is NOT fully tested nor truly implemented just yet.
    \todo OpenMP parallel finding likely has other opmisations that can be implemented to reduce compute time.
*/
Int_t* SearchFullSet(Options &opt, const Int_t nbodies, vector<Particle> &Part, Int_t &numgroups)
{
    Int_t i, *pfof,*pfoftemp, minsize;
    FOFcompfunc fofcmp;
    FOFcheckfunc fofcheck;
    fstream Fout;
    char fname[2000];
    Double_t param[20];
    Double_t *period=NULL;
    Double_t vscale2,mtotregion,vx,vy,vz;
    Double_t *vscale2array;
    Coordinate vmean(0,0,0);
    int maxnthreads,nthreads=1,tid;
    Int_tree_t *Len,*Head,*Next,*Tail;
    Int_t *storetype,*storeorgIndex;
    Int_t *ids, *numingroup=NULL, *noffset;
    Int_t *id_3dfof_of_6dfof;
    Int_t ng,npartingroups;
    Int_t totalgroups;
    Double_t time1,time2, time3;
    KDTree *tree;
    KDTree **tree3dfofomp;
    Int_t *p3dfofomp;
    int iorder = 1;
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
    Int_t Nlocal=nbodies;
#endif
#ifdef USEOPENMP
#pragma omp parallel
    {
    if (omp_get_thread_num()==0) maxnthreads=omp_get_num_threads();
    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
    OMP_Domain *ompdomain;
#endif
    if (opt.p>0) {
        period=new Double_t[3];
        for (int j=0;j<3;j++) period[j]=opt.p;
    }
    psldata=new StrucLevelData;

    minsize=opt.HaloMinSize;
#ifdef USEMPI
    //if using MPI, lower minimum number
    if (NProcs>1) minsize=MinNumMPI;
    iorder = 0;
#endif

    time1=MyGetTime();
    time2=MyGetTime();
    cout<<"Begin FOF search  of entire particle data set ... "<<endl;
    param[0]=tree->TPHYS;
    param[1]=(opt.ellxscale*opt.ellxscale)*(opt.ellphys*opt.ellphys)*(opt.ellhalophysfac*opt.ellhalophysfac);
    param[6]=param[1];
    cout<<"First build tree ... "<<endl;
#ifdef USEOPENMP
    //if using openmp produce tree with large buckets as a decomposition of the local mpi domain
    //to then run local fof searches on each domain before stitching
    int numompregions = ceil(nbodies/(float)opt.openmpfofsize);
    if (numompregions >= 4 && nthreads > 1 && opt.iopenmpfof == 1) {
        time3=MyGetTime();
        Double_t rdist = sqrt(param[1]);
        //determine the omp regions;
        tree = new KDTree(Part.data(),nbodies,opt.openmpfofsize,tree->TPHYS,tree->KEPAN,100);
        numompregions=tree->GetNumLeafNodes();
        ompdomain = OpenMPBuildDomains(opt, numompregions, tree, rdist);
        storeorgIndex = new Int_t[nbodies];
        for (i=0;i<nbodies;i++) storeorgIndex[i]=Part[i].GetID();
        //build local trees
        tree3dfofomp = OpenMPBuildLocalTrees(opt, numompregions, Part, ompdomain, period);
        if (opt.iverbose) cout<<ThisTask<<": finished building "<<numompregions<<" domains and trees "<<MyGetTime()-time3<<endl;
    }
    else {
        time3=MyGetTime();
        tree = new KDTree(Part.data(),nbodies,opt.Bsize,tree->TPHYS,tree->KEPAN,1000,0,0,0,period);
        if (opt.iverbose) cout<<ThisTask<<": finished building single tree with single OpenMP "<<MyGetTime()-time3<<endl;
    }

#else
    tree=new KDTree(Part.data(),nbodies,opt.Bsize,tree->TPHYS,tree->KEPAN,1000,0,0,0,period);
#endif
    cout<<"Done"<<endl;
    cout<<ThisTask<<" Search particles using 3DFOF in physical space"<<endl;
    cout<<ThisTask<<" Parameters used are : ellphys="<<sqrt(param[6])<<" Lunits (ell^2="<<param[6]<<" and likely "<<sqrt(param[6])/opt.ellxscale<<" in interparticle spacing"<<endl;
    if (opt.partsearchtype==PSTALL && opt.iBaryonSearch>1) {fofcmp=&FOF3dDM;param[7]=DARKTYPE;fofcheck=FOFchecktype;}
    else fofcmp=&FOF3d;


#ifdef USEMPI
    //if using mpi no need to locally sort just yet and might as well return the Head, Len, Next arrays
    Head=new Int_tree_t[nbodies];Next=new Int_tree_t[nbodies];
#else
    Head=NULL;Next=NULL;
#endif

#ifdef USEOPENMP
    //if enough regions then search each individually
    //then link across omp domains
    Int_t ompminsize = 2;
    if (numompregions>=4 && nthreads > 1 && opt.iopenmpfof == 1){
        time3=MyGetTime();
        Int_t orgIndex, omp_import_total;
        int omptask;
        Double_t rdist = sqrt(param[1]);
        pfof = new Int_t[nbodies];
        for (i=0;i<nbodies;i++)pfof[i]=0;
#ifndef USEMPI
        Head = new Int_tree_t[nbodies];
        Next = new Int_tree_t[nbodies];
#endif
        Int_t *omp_nrecv_total = new Int_t[numompregions];
        Int_t *omp_nrecv_offset = new Int_t[numompregions];
        OMP_ImportInfo *ompimport;

        //get fof in each region
        numgroups = OpenMPLocalSearch(opt,
            nbodies, Part, pfof, storeorgIndex,
            Head, Next,
            tree3dfofomp, param, rdist, ompminsize, fofcmp,
            numompregions, ompdomain);
        if (opt.iverbose) cout<<ThisTask<<": finished omp local search of "<<numompregions<<" containing total of "<<numgroups<<" groups "<<MyGetTime()-time3<<endl;
        if (numgroups > 0) {

        //then for each omp region determine the particles to "import" from other omp regions
        ompimport = OpenMPImportParticles(opt, nbodies, Part, pfof, storeorgIndex,
            numompregions, ompdomain, rdist,
            omp_nrecv_total, omp_nrecv_offset, omp_import_total);
        if (omp_import_total > 0) {
            OpenMPLinkAcross(opt, nbodies, Part, pfof, storeorgIndex, Head, Next,
                param, fofcheck, numompregions, ompdomain, tree3dfofomp,
                omp_nrecv_total, omp_nrecv_offset, ompimport);

            delete[] ompimport;
            }
        }
        //free memory
#ifndef USEMPI
        delete[] Head;
        delete[] Next;
#endif
        delete[] omp_nrecv_total;
        delete[] omp_nrecv_offset;

        #pragma omp parallel default(shared) \
        private(i)
        {
        #pragma omp for schedule(dynamic) nowait
        for (i=0;i<numompregions;i++) delete tree3dfofomp[i];
        }
        delete[] tree3dfofomp;
        delete[] ompdomain;

        //reset particle ids to before omp tree built
        for (i=0;i<nbodies;i++) Part[i].SetID(storeorgIndex[i]);
        delete[] storeorgIndex;

        //delete coarse omp tree and rebuild fine tree;
        delete tree;
        tree = NULL;
        //resort particles and group ids
        if (numgroups > 0) {
            numgroups = OpenMPResortParticleandGroups(nbodies, Part, pfof, minsize);
        }

        //and reallocate tree if required (that is only if not using MPI but searching for substructure
#if !defined(USEMPI) && defined(STRUCDEN)
        if (numgroups>0 && (opt.iSubSearch==1&&opt.foftype!=FOF6DCORE))
#endif
        tree = new KDTree(Part.data(),nbodies,opt.Bsize,tree->TPHYS,tree->KEPAN,1000,0,0,0,period);
        //if running MPI then need to pudate the head, next info
#ifdef USEMPI
        OpenMPHeadNextUpdate(nbodies, Part, numgroups, pfof, Head, Next);
#endif

    }
    else {
        //posible alteration for all particle search
        if (opt.partsearchtype==PSTALL && opt.iBaryonSearch>1) pfof=tree->FOFCriterionSetBasisForLinks(fofcmp,param,numgroups,minsize,iorder,0,FOFchecktype,Head,Next);
        else pfof=tree->FOF(sqrt(param[1]),numgroups,minsize,iorder,Head,Next);
    }
#else
    //posible alteration for all particle search
    if (opt.partsearchtype==PSTALL && opt.iBaryonSearch>1) pfof=tree->FOFCriterionSetBasisForLinks(fofcmp,param,numgroups,minsize,iorder,0,FOFchecktype,Head,Next);
    else pfof=tree->FOF(sqrt(param[1]),numgroups,minsize,iorder,Head,Next);
#endif

#ifndef USEMPI
    totalgroups=numgroups;
    //if this flag is set, calculate localfield value here for particles possibly resident in a field structure
#ifdef STRUCDEN
    if (numgroups>0 && (opt.iSubSearch==1&&opt.foftype!=FOF6DCORE)) {
    numingroup=BuildNumInGroup(nbodies, numgroups, pfof);
    storetype=new Int_t[nbodies];
    for (i=0;i<nbodies;i++) storetype[i]=Part[i].GetType();
    Int_t numinstrucs=0;
    //if not searching all particle then searching for baryons associated with substructures, then set type to group value
    //so that velocity density just calculated for particles in groups (type>0)
    if (!(opt.iBaryonSearch>=1 && opt.partsearchtype==PSTALL)) for (i=0;i<nbodies;i++) Part[i].SetType(numingroup[pfof[Part[i].GetID()]]>=MINSUBSIZE);
    //otherwise set type to group value for dark matter
    else {
        for (i=0;i<nbodies;i++) {
            if (Part[i].GetType()==DARKTYPE) Part[i].SetType(numingroup[pfof[Part[i].GetID()]]>=MINSUBSIZE);
            else Part[i].SetType(-1);
        }
    }
    for (i=0;i<nbodies;i++) if (Part[i].GetType()>0) numinstrucs++;
    if (opt.iverbose) cout<<"Number of particles in large subhalo searchable structures "<<numinstrucs<<endl;
    if (numinstrucs>0) GetVelocityDensity(opt, nbodies, Part.data(), tree);

    for (i=0;i<nbodies;i++) Part[i].SetType(storetype[i]);
    delete[] storetype;
    if (opt.fofbgtype>FOF6D) delete[] numingroup;
    }
#endif
    delete tree;
#endif

#ifdef USEMPI
    if (NProcs==1) {
        totalgroups=numgroups;
        if (tree != NULL) delete tree;
    }
    else {
    mpi_foftask=MPISetTaskID(Nlocal);

    Len=new Int_tree_t[nbodies];
    if (opt.iverbose>=2) {
        Int_t sum=0;
        for (i=0;i<nbodies;i++) sum+=(pfof[i]>0);
        cout<<ThisTask<<" has found locally "<<numgroups<<" with lower min size of "<<minsize<<", with  "<<sum<<" particles in all groups"<<endl;
    }
    if (numgroups) {
        numingroup=BuildNumInGroup(nbodies, numgroups, pfof);
        for (i=0;i<nbodies;i++) Len[i]=numingroup[pfof[Part[i].GetID()]];
        delete[] numingroup;
        numingroup=NULL;
    }
    time2=MyGetTime();

    //Also must ensure that group ids do not overlap between mpi threads so adjust group ids
    MPI_Allgather(&numgroups, 1, MPI_Int_t, mpi_ngroups, 1, MPI_Int_t, MPI_COMM_WORLD);
    MPIAdjustLocalGroupIDs(nbodies, pfof);
    //then determine export particles, declare arrays used to export data
#ifdef MPIREDUCEMEM
#ifdef SWIFTINTERFACE
    MPIGetExportNumUsingMesh(libvelociraptorOpt, nbodies, Part.data(), sqrt(param[1]));
#else
    MPIGetExportNum(nbodies, Part.data(), sqrt(param[1]));
#endif
#endif
    //allocate memory to store info
    cout<<ThisTask<<": Finished local search, nexport/nimport = "<<NExport<<" "<<NImport<<" in "<<MyGetTime()-time2<<endl;

    PartDataIn = new Particle[NExport];
    PartDataGet = new Particle[NImport];
    FoFDataIn = new fofdata_in[NExport];
    FoFDataGet = new fofdata_in[NImport];
    //if using MPI must determine which local particles need to be exported to other threads and used to search
    //that threads particles. This is done by seeing if the any particles have a search radius that overlaps with
    //the boundaries of another threads domain. Then once have exported particles must search local particles
    //relative to these exported particles. Iterate over search till no new links are found.

    //I have adjusted FOF data structure to have local group length and also seperated the export particles from export fof data
    //the reason is that will have to update fof data in iterative section but don't need to update particle information.
#ifdef SWIFTINTERFACE
    MPIBuildParticleExportListUsingMesh(libvelociraptorOpt, nbodies, Part.data(), pfof, Len, sqrt(param[1]));
#else
    MPIBuildParticleExportList(nbodies, Part.data(), pfof, Len, sqrt(param[1]));
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    //Now that have FoFDataGet (the exported particles) must search local volume using said particles
    //This is done by finding all particles in the search volume and then checking if those particles meet the FoF criterion
    //One must keep iterating till there are no new links.
    //Wonder if i don't need another loop and a final check
    Int_t links_across,links_across_total;
    cout<<ThisTask<<": Starting to linking across MPI domains"<<endl;
    do {
        if (opt.partsearchtype==PSTALL && opt.iBaryonSearch>1) {
            links_across=MPILinkAcross(nbodies, tree, Part.data(), pfof, Len, Head, Next, param[1], fofcheck, param);
        }
        else {
            links_across=MPILinkAcross(nbodies, tree, Part.data(), pfof, Len, Head, Next, param[1]);
        }
        if (opt.iverbose>=2) {
            cout<<ThisTask<<" has found "<<links_across<<" links to particles on other mpi domains "<<endl;
        }
        MPI_Allreduce(&links_across, &links_across_total, 1, MPI_Int_t, MPI_SUM, MPI_COMM_WORLD);
        MPIUpdateExportList(nbodies,Part.data(),pfof,Len);
    }while(links_across_total>0);
    if (ThisTask==0) cout<<ThisTask<<": finished linking across MPI domains in "<<MyGetTime()-time2<<endl;

    delete[] FoFDataIn;
    delete[] FoFDataGet;
    delete[] PartDataIn;
    delete[] PartDataGet;

    //reorder local particle array and delete memory associated with Head arrays, only need to keep Particles, pfof and some id and idexing information
    delete tree;
    delete[] Head;
    delete[] Next;
    delete[] Len;
    //Now redistribute groups so that they are local to a processor (also orders the group ids according to size
    opt.HaloMinSize=MinNumOld;//reset minimum size
    Int_t newnbodies=MPIGroupExchange(nbodies,Part.data(),pfof);
    //once groups are local, can free up memory. Might need to increase size
    //of vector
    if (Nmemlocal<Nlocal) {
        Part.resize(Nlocal);
        Nmemlocal=Nlocal;
    }
    delete[] mpi_foftask;
    delete[] pfof;
    pfof=new Int_t[newnbodies];
    //And compile the information and remove groups smaller than minsize
    numgroups=MPICompileGroups(newnbodies,Part.data(),pfof,opt.HaloMinSize);
    //and free up some memory if vector doesn't need to be as big
    if (Nmemlocal>Nlocal) {Part.resize(Nlocal);Nmemlocal=Nlocal;}
    cout<<"MPI thread "<<ThisTask<<" has found "<<numgroups<<endl;
    //free up memory now that only need to store pfof and global ids
    totalgroups=0;
    for (int j=0;j<NProcs;j++) totalgroups+=mpi_ngroups[j];
    Nlocal=newnbodies;
    }
#endif
    if (opt.iverbose>=2) {
        Int_t sum=0;
        for (i=0;i<Nlocal;i++) sum+=(pfof[i]>0);
        cout<<ThisTask<<" has found after full search "<<numgroups<<" with lower min size of "<<minsize<<", with  "<<sum<<" particles in all groups"<<endl;
    }
    if (ThisTask==0) cout<<"Total number of groups found is "<<totalgroups<<endl;
    if (ThisTask==0) cout<<ThisTask<<": finished FOF search in total time of "<<MyGetTime()-time1<<endl;

    //if calculating velocity density only of particles resident in field structures large enough for substructure search
#if defined(STRUCDEN) && defined(USEMPI)
    if (totalgroups>0&&(opt.iSubSearch==1&&opt.foftype!=FOF6DCORE))
    {
        storetype=new Int_t[Nlocal];
        Int_t numinstrucs=0,numlocalden=0;
        for (i=0;i<Nlocal;i++) storetype[i]=Part[i].GetType();
        if (!(opt.iBaryonSearch>=1 && opt.partsearchtype==PSTALL)) {
            numingroup=BuildNumInGroup(Nlocal, numgroups, pfof);
            for (i=0;i<Nlocal;i++) {
                Part[i].SetType((numingroup[pfof[i]]>=MINSUBSIZE));
                numlocalden += (Part[i].GetType()>0);
            }
        }
        //otherwise set type to group value for dark matter
        else {
            numingroup=BuildNumInGroupTyped(Nlocal,numgroups,pfof,Part.data(),DARKTYPE);
            for (i=0;i<Nlocal;i++) {
                if (Part[i].GetType()==DARKTYPE) Part[i].SetType(numingroup[pfof[Part[i].GetID()]]>=MINSUBSIZE);
                else Part[i].SetType(0);
                numlocalden += (Part[i].GetType()>0);
            }
        }
        for (i=0;i<Nlocal;i++) {numinstrucs+=(pfof[i]>0);}
        if (opt.iverbose) cout<<ThisTask<<" has "<<numinstrucs<<" particles for which density must be calculated"<<endl;
        cout<<ThisTask<<" Going to build tree "<<endl;
        tree=new KDTree(Part.data(),Nlocal,opt.Bsize,tree->TPHYS,tree->KEPAN,100,0,0,0,period);
        GetVelocityDensity(opt, Nlocal, Part.data(),tree);
        delete tree;
        for (i=0;i<Nlocal;i++) Part[i].SetType(storetype[i]);
        delete[] storetype;
    }
#endif

    //if search was periodic, alter particle positions in structures so substructure search no longer has to be
    //periodic
    if (opt.p>0&&numgroups>0) {
        if (numgroups>0)AdjustStructureForPeriod(opt,Nlocal,Part,numgroups,pfof);
        delete[] period;
    }

    //have now 3dfof groups local to a MPI thread and particles are back in index order that will be used from now on
    //note that from on, use Nlocal, which is altered in mpi but set to nbodies in non-mpi
    if (opt.fofbgtype<=FOF6D && totalgroups>0) {
    ///\todo In the 6DFOF, particles are sorted into group order but then resorted back into input index order
    ///this is not necessary if the tipsy still .grp array does not need to be constructed.
    ///we would removing storing the old ids alter how the local group lists are stiched together into the larger array
    ///and remove a slow sort back into original ID order
    time1=MyGetTime();
    time2=MyGetTime();

    //now if 6dfof search to overcome issues with 3DFOF by searching only particles that have been previously linked by the 3dfof
    //adjust physical linking length
    param[1]*=opt.ellhalo6dxfac*opt.ellhalo6dxfac;
    param[6]=param[1];
    //set min size and fof criterion
    minsize=opt.HaloMinSize;
    if (opt.fofbgtype!=FOFSTNOSUBSET) fofcmp=&FOF6d;
    else fofcmp=&FOFStream;

    cout<<ThisTask<<": Sorting particles for 6dfof/phase-space search "<<endl;
    npartingroups=0;
    //sort particles so that largest group is first, 2nd next, etc with untagged at end.
    Int_t iend=0;

    storetype=new Int_t[Nlocal];
    if (numingroup==NULL) numingroup=new Int_t[numgroups+1];
    noffset=new Int_t[numgroups+1];
    for (i=0;i<=numgroups;i++) numingroup[i]=noffset[i]=0;
    for (i=0;i<Nlocal;i++) {
        storetype[i]=Part[i].GetPID();
        Part[i].SetPID((pfof[i]==0)*Nlocal+(pfof[i]>0)*pfof[i]);
        npartingroups+=(Int_t)(pfof[i]>0);
        iend+=(pfof[i]==1);
        numingroup[pfof[i]]++;
    }
    for (i=2;i<=numgroups;i++) noffset[i]=noffset[i-1]+numingroup[i-1];
    qsort(Part.data(), Nlocal, sizeof(Particle), PIDCompare);
    //sort(Part.begin(),Part.end(),PIDCompareVec);
    for (i=0;i<Nlocal;i++) Part[i].SetPID(storetype[Part[i].GetID()]);
    delete[] storetype;
    //store index order
    ids=new Int_t[Nlocal];
    for (i=0;i<Nlocal;i++) ids[i]=Part[i].GetID();

    //if only using single velocity scale, use the largest "halo" to determine an appropriate velocity scale
    if (opt.fofbgtype==FOF6D && opt.iKeepFOF==0) {
        vscale2=mtotregion=vx=vy=vz=0;
        for (i=0;i<iend;i++) {
            vx+=Part[i].GetVelocity(0)*Part[i].GetMass();
            vy+=Part[i].GetVelocity(1)*Part[i].GetMass();
            vz+=Part[i].GetVelocity(2)*Part[i].GetMass();
        }
            mtotregion+=Part[i].GetMass();
        vmean[0]=vx/mtotregion;vmean[1]=vy/mtotregion;vmean[2]=vz/mtotregion;
        for (i=0;i<iend;i++) {
            for (int j=0;j<3;j++) vscale2+=pow(Part[i].GetVelocity(j)-vmean[j],2.0)*Part[i].GetMass();
        }
        if (mtotregion>0) vscale2/=mtotregion;

#ifdef USEMPI
        Double_t mpi_vscale2;
        MPI_Allreduce(&vscale2,&mpi_vscale2,1,MPI_Real_t,MPI_MAX,MPI_COMM_WORLD);
        vscale2=mpi_vscale2;
#endif
        //to account for fact that dispersion may not be high enough to link all outlying regions
        //scale by 6d velocity factor
        vscale2*=opt.ellhalo6dvfac*opt.ellhalo6dvfac;
        //set the velocity scale
        param[2]=(vscale2);
        param[7]=param[2];
        cout<<"Search "<<npartingroups<<" particles using 6DFOF with uniform velocity scale"<<endl;
        cout<<"Parameters used are : ellphys="<<sqrt(param[6])<<" Lunits, ellvel="<<sqrt(param[7])<<" Vunits."<<endl;
    }
    //otherwise each object has its own velocity scale
    else if(opt.fofbgtype==FOF6DADAPTIVE || opt.iKeepFOF){
        vscale2array=new Double_t[numgroups+1];
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,vscale2,mtotregion,vx,vy,vz,vmean)
{
#pragma omp for schedule(dynamic,1) nowait
#endif
        for (i=1;i<=numgroups;i++) {
            vscale2=mtotregion=vx=vy=vz=0;
            for (Int_t j=0;j<numingroup[i];j++) {
                vx+=Part[j+noffset[i]].GetVelocity(0)*Part[j+noffset[i]].GetMass();
                vy+=Part[j+noffset[i]].GetVelocity(1)*Part[j+noffset[i]].GetMass();
                vz+=Part[j+noffset[i]].GetVelocity(2)*Part[j+noffset[i]].GetMass();
                mtotregion+=Part[j+noffset[i]].GetMass();
            }
            vmean[0]=vx/mtotregion;vmean[1]=vy/mtotregion;vmean[2]=vz/mtotregion;
            for (Int_t j=0;j<numingroup[i];j++) {
                for (int k=0;k<3;k++) vscale2+=pow(Part[j+noffset[i]].GetVelocity(k)-vmean[k],2.0)*Part[j+noffset[i]].GetMass();
            }
            vscale2array[i]=vscale2/mtotregion*opt.ellhalo6dvfac*opt.ellhalo6dvfac;;
        }
#ifdef USEOPENMP
}
#endif
        cout<<"Search "<<npartingroups<<" particles using 6DFOF with adaptive velocity scale"<<endl;
        cout<<"Static parameters used are : ellphys="<<sqrt(param[6])<<" Lunits"<<endl;
    }
    //use the phase-space stream finding parameters
    else if (opt.fofbgtype==FOFSTNOSUBSET) {
        ///\todo not implemented yet so quit and spit error message
        if (ThisTask==0) cerr<<" THIS TYPE OF PHASE-SPACE HALO SEARCH not implemented yet, quiting. Please use halo search of 3,4, or 5 (adaptive 6d, 6d with single v scale, and 3d only)"<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,8);
#else
        exit(8);
#endif
    }
    //error, so write message and quit
    else {
        if (ThisTask==0) cerr<<" THIS TYPE OF PHASE-SPACE HALO SEARCH not implemented yet, quiting. Please use halo search of 3,4, or 5 (adaptive 6d, 6d with single v scale, and 3d only)"<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,8);
#else
        exit(8);
#endif
    }

    Head=new Int_tree_t[Nlocal];Next=new Int_tree_t[Nlocal];Tail=new Int_tree_t[Nlocal];Len=new Int_tree_t[Nlocal];
    KDTree *treeomp[nthreads];
    Double_t *paramomp=new Double_t[nthreads*20];
    Int_t **pfofomp;
    Int_t *ngomp;
    //now split search by group. If doing fixed velocity dispersion search can group together all small
    //3DFOF groups into a single search
    //otherwise don't group together all small groups
    iend=numgroups;
    if (opt.fofbgtype==FOF6D || opt.fofbgtype==FOFSTNOSUBSET) {
        //determine chunks to search, where groups are split down to objects that are small and can be grouped
        for (i=1;i<=numgroups;i++) if (numingroup[i]<=ompunbindnum) {iend=i;break;}
        numingroup[0]=0;
        for (i=iend;i<=numgroups;i++) numingroup[0]+=numingroup[i];
        numingroup[iend]=numingroup[0];
    }
    //copy the parameters appropriately to the paramsomp array which can be used by openmp
    for (i=0;i<nthreads;i++)
        for (int j=0;j<20;j++) paramomp[j+i*20]=param[j];

    ///\todo need to improve kdtree 6dfof construction to make use of scaling dimensions and running in 6d.
    ///before ran FOF criterion on physical tree but try scaling particles according to linking lengths, run
    ///6d phase tree and simple FOF ball search
    pfofomp=new Int_t*[iend+1];
    ngomp=new Int_t[iend+1];
    Double_t xscaling, vscaling;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,tid,xscaling,vscaling)
{
#pragma omp for schedule(dynamic,1) nowait
#endif
    for (i=1;i<=iend;i++) {
#ifdef USEOPENMP
        tid=omp_get_thread_num();
#else
        tid=0;
#endif
        //if adaptive 6dfof, set params
        if (opt.fofbgtype==FOF6DADAPTIVE) paramomp[2+tid*20]=paramomp[7+tid*20]=vscale2array[i];
        /*
        treeomp[tid]=new KDTree(&Part[noffset[i]],numingroup[i],opt.Bsize,treeomp[tid]->TPHYS,tree->KEPAN,100);
        pfofomp[i]=treeomp[tid]->FOFCriterion(fofcmp,&paramomp[tid*20],ngomp[i],minsize,1,0,Pnocheck,&Head[noffset[i]],&Next[noffset[i]],&Tail[noffset[i]],&Len[noffset[i]]);
        */
        //scale particle positions
        xscaling=1.0/sqrt(paramomp[1+tid*20]);vscaling=1.0/sqrt(paramomp[2+tid*20]);
        for (Int_t j=0;j<numingroup[i];j++) {
            Part[noffset[i]+j].ScalePhase(xscaling,vscaling);
        }
        xscaling=1.0/xscaling;vscaling=1.0/vscaling;
        treeomp[tid]=new KDTree(&(Part.data()[noffset[i]]),numingroup[i],opt.Bsize,treeomp[tid]->TPHS,tree->KEPAN,100);
        pfofomp[i]=treeomp[tid]->FOF(1.0,ngomp[i],minsize,1,&Head[noffset[i]],&Next[noffset[i]],&Tail[noffset[i]],&Len[noffset[i]]);
        delete treeomp[tid];
        for (Int_t j=0;j<numingroup[i];j++) {
            Part[noffset[i]+j].ScalePhase(xscaling,vscaling);
        }
    }
#ifdef USEOPENMP
}
#endif

    //now if keeping original 3DFOF structures (useful for stellar haloes search) then store original number of 3d fof haloes
    if (opt.iKeepFOF)
    {
        if (opt.iverbose>=2 && ThisTask==0) cout<<"Storing the 3D fof envelopes of the 6d fof structures found"<<endl;
        //store current number of 6dfof groups
        for (i=1;i<=iend;i++) ng+=ngomp[i];
        Int_t *pfof6dfof=new Int_t[Nlocal];
        for (i=0;i<Nlocal;i++) pfof6dfof[i]=0;
        ng=0;
        ///\todo restructure how data is stored in the pfof6dfof since will no longer necessarily keep stuff in index order
        for (i=1;i<=iend;i++) {
            for (int j=0;j<numingroup[i];j++) {
                pfof6dfof[ids[Part[noffset[i]+j].GetID()+noffset[i]]]=pfofomp[i][j]+(pfofomp[i][j]>0)*ng;
            }
            ng+=ngomp[i];
            delete[] pfofomp[i];
        }

        //now that all groups have been processed, need to adjust the number of 3DFOF groups as a 6dFOF group can span a entire 3DFOF group
        opt.num3dfof = numgroups;

        // store whether the 3dfof group remains
        int *i3dfofleft = new int[numgroups+1];
        //store 3dfof id of a 6dfof group
        Int_t *id3dfof = new Int_t[ng+1];

        // used to quickly access myigmid
        Int_t *n3dfofoffset = new Int_t[numgroups+1];
        n3dfofoffset[0] = n3dfofoffset[1] = 0;
        for (i=2;i<=numgroups;i++) n3dfofoffset[i] = n3dfofoffset[i-1] + ngomp[i-1];

        if (opt.iverbose) cout <<ThisTask<<" Before 6dfof search had "<< numgroups <<" 3dfof groups "<<endl;
        Int_t numin3dfof;
        for (i=1;i<=numgroups;i++)
        {
            numin3dfof=0;
            for (int j = 0; j < numingroup[i]; j++) numin3dfof += (pfof6dfof[ids[Part[noffset[i]+j].GetID()+noffset[i]]]==0);
            if (numin3dfof==0)
            {
                opt.num3dfof--;
                i3dfofleft[i]=0;
            }
            else i3dfofleft[i]=1;
        }
        if (opt.iverbose) cout<<ThisTask<<" now have" << opt.num3dfof << endl;

        // now have corrected the number of opt.num3dfof.
        // allocate the array storing the id of the 3dfof to which the 6dfof belongs
        id_3dfof_of_6dfof = new Int_t [opt.num3dfof + ng + 1];
        for (i = 1; i <= opt.num3dfof; i++)  id_3dfof_of_6dfof[i]=0;

        // now have corrected the number of opt.num3dfof. Adjust original pfof
        opt.num3dfof = 0;
        for (i=1;i<=numgroups;i++)
        {
            // if 3dfof remains
            if (i3dfofleft[i] == 1)
            {
                opt.num3dfof++;
                for (int j = 0; j < numingroup[i]; j++) pfof[ids[Part[noffset[i]+j].GetID()+noffset[i]]]=opt.num3dfof;

                for (int j=0;j<ngomp[i];j++) id3dfof[j+n3dfofoffset[i]]=opt.num3dfof;
            }
            // if no igm then all groups are not associated with an 3dfof
            else
                for (int j=0;j<ngomp[i];j++)  id3dfof[j+n3dfofoffset[i]]=0;
        }
        delete[] i3dfofleft;
        delete[] n3dfofoffset;
        for (i = 0; i < ng; i++) id_3dfof_of_6dfof[i + opt.num3dfof+1] = id3dfof[i];
        delete[] id3dfof;
        //now adjust 6dfof ids
        for (i = 0; i < Nlocal; i++) if (pfof6dfof[i] != 0) pfof[i] = pfof6dfof[i] + opt.num3dfof;
        //increase number of groups
        ng+=opt.num3dfof;
        delete[] pfof6dfof;
    }
    //if not keeping 3dfof just overwrite the pfof array
    else {
        ng=0;
        for (i=0;i<Nlocal;i++) pfof[i]=0;
        for (i=1;i<=iend;i++) {
            for (int j=0;j<numingroup[i];j++) {
                pfof[ids[Part[noffset[i]+j].GetID()+noffset[i]]]=pfofomp[i][j]+(pfofomp[i][j]>0)*ng;
            }
            ng+=ngomp[i];
            delete[] pfofomp[i];
        }
    }
    delete[] ngomp;
    delete[] paramomp;
    delete[] pfofomp;
    delete[] noffset;
    delete[] numingroup;
    delete[] Head;
    delete[] Tail;
    delete[] Next;
    delete[] Len;

    //reorder ids in descending group size order only if not keeping FOF
    if (ng>0 && opt.iKeepFOF==0) {
        if (opt.iverbose) cout<<" reordering "<<ng<<" groups "<<endl;
        numingroup=BuildNumInGroup(Nlocal, ng, pfof);
        Int_t **pglist=BuildPGList(Nlocal, ng, numingroup, pfof);
        ReorderGroupIDs(ng, ng, numingroup, pfof, pglist);
        for (i=1;i<=ng;i++) delete[] pglist[i];
        delete[] pglist;
        delete[] numingroup;
    }

    ///\todo only run this sort if necessary to keep id order
    for (i=0;i<npartingroups;i++) Part[i].SetID(ids[i]);
    gsl_heapsort(Part.data(), Nlocal, sizeof(Particle), IDCompare);
    //sort(Part.begin(), Part.end(), IDCompareVec);
    delete[] ids;
    numgroups=ng;

    }
    //end of extra 6dfof/phase-space search of 3D FOF envelops

    //now if not search for substructure but want bound halos need to check binding
    if (opt.iBoundHalos>=1) {
        CheckUnboundGroups(opt,Nlocal,Part.data(),numgroups,pfof);
#ifdef USEMPI
        if (ThisTask==0) cout<<ThisTask<<" After unnbinding halos"<<endl;
        //update number of groups if extra secondary search done
        if (opt.fofbgtype>FOF6D) {
            cout<<"MPI thread "<<ThisTask<<" has found "<<numgroups<<endl;
            MPI_Allgather(&numgroups, 1, MPI_Int_t, mpi_ngroups, 1, MPI_Int_t, MPI_COMM_WORLD);
            //free up memory now that only need to store pfof and global ids
            if (ThisTask==0) {
                int totalgroups=0;
                for (int j=0;j<NProcs;j++) totalgroups+=mpi_ngroups[j];
                cout<<"Total number of groups found is "<<totalgroups<<endl;
            }
        }
#endif
    }

#ifdef USEMPI
    //update number of groups if extra secondary search done
    if (opt.fofbgtype<=FOF6D) {
    cout<<"MPI thread "<<ThisTask<<" has found "<<numgroups<<endl;
    MPI_Allgather(&numgroups, 1, MPI_Int_t, mpi_ngroups, 1, MPI_Int_t, MPI_COMM_WORLD);
    //free up memory now that only need to store pfof and global ids
    if (ThisTask==0) {
        int totalgroups=0;
        for (int j=0;j<NProcs;j++) totalgroups+=mpi_ngroups[j];
        cout<<"Total number of groups found is "<<totalgroups<<endl;
    }
    if (ThisTask==0) cout<<ThisTask<<" finished 6d/phase-space fof search in "<<MyGetTime()-time2<<endl;
    }
    MPI_Allgather(&numgroups, 1, MPI_Int_t, mpi_nhalos, 1, MPI_Int_t, MPI_COMM_WORLD);
#endif

    //now that field structures have been identified, allocate enough memory for the psldata pointer,
    //allocate memory for lowest level in the substructure hierarchy, corresponding to field objects
    if (opt.iverbose) cout<<ThisTask<<" Now store hierarchy information "<<endl;
    //standard operation is to not keep the 3DFOF envelopes once 6dfof run
    if (!(opt.iKeepFOF && opt.fofbgtype<=FOF6D))
    {
        psldata->Allocate(numgroups);
        psldata->Initialize();
        for (i=0;i<Nlocal;i++) {
            if (pfof[i]>0) {
            if (psldata->gidhead[pfof[i]]==NULL) {
                //set the group id head pointer to the address within the pfof array
                psldata->gidhead[pfof[i]]=&pfof[i];
                //set particle pointer to particle address
                psldata->Phead[pfof[i]]=&Part[i];
                //set the parent pointers to appropriate addresss such that the parent and uber parent are the same as the groups head
                psldata->gidparenthead[pfof[i]]=&pfof[i];
                psldata->giduberparenthead[pfof[i]]=&pfof[i];
                //set structure type
                psldata->stypeinlevel[pfof[i]]=HALOSTYPE;
            }
            }
        }
        psldata->stype=HALOSTYPE;
    }
    else
    {
        // initialize the 3dfof level of the hierarchy
        psldata->Allocate(opt.num3dfof);
        psldata->Initialize();
        psldata->stype = FOF3DTYPE;

        for (i = 0; i < Nlocal; i++)
        {
            if ((pfof[i] <= opt.num3dfof) && (pfof[i] > 0)) {
            if (psldata->gidhead[pfof[i]] == NULL)
            {
                psldata->gidhead[pfof[i]] = &pfof[i];
                psldata->Phead[pfof[i]] = &Part[i];
                psldata->gidparenthead[pfof[i]] = &pfof[i];
                psldata->giduberparenthead[pfof[i]] = &pfof[i];
                psldata->stypeinlevel[pfof[i]] = FOF3DTYPE;
            }
            }
        }
        //adjust ng to 6dfof number
        ng=numgroups-opt.num3dfof;
        if (ng>0) {

            //reorder in descending 6dfof group size order only if keeping FOF
            numingroup=BuildNumInGroup(Nlocal, numgroups, pfof);
            Int_t **pglist=BuildPGList(Nlocal, numgroups, numingroup, pfof);
            Int_t *value6d3d=new Int_t[numgroups+1];
            //store size of largest 6dfof group to offset 3dfof when ordering 6d, leaving 3d order unchanged
            Int_t maxval=0;
            for (i=opt.num3dfof+1;i<=numgroups;i++) if (maxval<numingroup[i]) maxval=numingroup[i];
            for (i=1;i<=opt.num3dfof;i++) value6d3d[i]=(opt.num3dfof-i)+maxval+1;
            for (i=opt.num3dfof+1;i<=numgroups;i++) value6d3d[i]=numingroup[i];
            //need to adjust id_3dfof_of_6dfof mapping as well when reordering groups
            ReorderGroupIDsAndArraybyValue(numgroups, numgroups, numingroup, pfof, pglist,value6d3d,id_3dfof_of_6dfof);
            for (i=1;i<=numgroups;i++) delete[] pglist[i];
            delete[] pglist;
            delete[] numingroup;
            delete[] value6d3d;

            //initialize next level which stores halos
            psldata->nextlevel = new StrucLevelData();
            psldata->nextlevel->Phead = new Particle*[ng+1];
            psldata->nextlevel->gidhead = new Int_t*[ng+1];
            psldata->nextlevel->Pparenthead = new Particle*[ng+1];
            psldata->nextlevel->gidparenthead = new Int_t*[ng+1];
            psldata->nextlevel->giduberparenthead = new Int_t*[ng+1];
            psldata->nextlevel->stypeinlevel = new Int_t[ng+1];
            psldata->nextlevel->stype = HALOSTYPE;
            psldata->nextlevel->nsinlevel = ng;
            psldata->nextlevel->nextlevel = NULL;

            for (i = 0; i <= ng; i++)
            {
                psldata->nextlevel->Phead[i] = NULL;
                psldata->nextlevel->gidhead[i] = NULL;
                psldata->nextlevel->gidparenthead[i] = NULL;
                psldata->nextlevel->Pparenthead[i] = NULL;
                psldata->nextlevel->giduberparenthead[i] = NULL;
                psldata->nextlevel->stypeinlevel[i] = 0;
            }

            for (i = 0; i < Nlocal; i++)
            {
                if (pfof[i] > opt.num3dfof) {
                if (psldata->nextlevel->gidhead[pfof[i]-opt.num3dfof] == NULL)
                {
                    psldata->nextlevel->gidhead[pfof[i]-opt.num3dfof] = &pfof[i];
                    psldata->nextlevel->Phead[pfof[i]-opt.num3dfof] = &Part[i];

                    ///for 6dfof groups that are contained in a 3dfof group make sure to point to head
                    if (id_3dfof_of_6dfof[pfof[i]] > 0)
                    {
                        psldata->nextlevel->gidparenthead[pfof[i]-opt.num3dfof] = psldata->gidhead[id_3dfof_of_6dfof[pfof[i]]];
                        psldata->nextlevel->Pparenthead[pfof[i]-opt.num3dfof] = psldata->Phead[id_3dfof_of_6dfof[pfof[i]]];
                        psldata->nextlevel->giduberparenthead[pfof[i]-opt.num3dfof] = psldata->gidhead[id_3dfof_of_6dfof[pfof[i]]];
                    }
                    //if not in 3dfof then halo is its own uberparent (don't need to set particle head as NULL is fine)
                    else
                    {
                        psldata->nextlevel->gidparenthead[pfof[i]-opt.num3dfof] = &pfof[i];
                        psldata->nextlevel->giduberparenthead[pfof[i]-opt.num3dfof] = &pfof[i];
                        psldata->nextlevel->Pparenthead[pfof[i]-opt.num3dfof] = &Part[i];
                    }
                    psldata->nextlevel->stypeinlevel[pfof[i]-opt.num3dfof] = HALOSTYPE;
                }
                }
            }
        }//end of ng>0 check
    }

    if (opt.iverbose) cout<<ThisTask<<" Done storing halo substructre level data"<<endl;
    return pfof;
}

void AdjustStructureForPeriod(Options &opt, const Int_t nbodies, vector<Particle> &Part, Int_t numgroups, Int_t *pfof)
{
    Int_t i,j, gid;
    //Int_t *numingroup, **pglist;
    //Coordinate c;
    double diff;
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif
    vector<Coordinate> refpos(numgroups+1);
    vector<Int_t> irefpos(numgroups+1);
    if (opt.iverbose) cout<<ThisTask<<" Adjusting for period "<<opt.p<<endl;
    for (i=1;i<=numgroups;i++) irefpos[i]=-1;
    for (i=0;i<nbodies;i++) {
        if (pfof[i]==0) continue;
        if (irefpos[pfof[i]] != -1) continue;
        refpos[pfof[i]] = Coordinate(Part[i].GetPosition());
        irefpos[pfof[i]] = i;
    }
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,diff,gid)
{
#pragma omp for
#endif
    for (i=0;i<nbodies;i++)
    {
        gid = pfof[i];
        if (gid==0) continue;
        if (irefpos[gid] == i) continue;
        for (j=0;j<3;j++) {
            diff=refpos[gid][j]-Part[i].GetPosition(j);
            if (diff<-0.5*opt.p) Part[i].SetPosition(j, Part[i].GetPosition(j)-opt.p);
            else if (diff>0.5*opt.p) Part[i].SetPosition(j, Part[i].GetPosition(j)+opt.p);
        }
    }
#ifdef USEOPENMP
}
#endif
}

//@}

///\name Search using outliers from background velocity distribution.
//@{

///Search subset
/*!
    \todo there are issues with locality for iterative search, significance check, regarding MPI threads. Currently these processes assume all necessary particles
    are locally accessible to the threads domain. For iSingleHalo==0, that is a full search of all fof halos has been done, that is the case as halos are localized
    to a MPI domain before the subsearch is made. However, in the case a of a single halo which has been broken up into several mpi domains, I must think carefully
    how the search should be localized. It should definitely be localized prior to CheckSignificance and the search window across mpi domains should use the larger
    physical search window used by the iterative search if that has been called.
 */
Int_t* SearchSubset(Options &opt, const Int_t nbodies, const Int_t nsubset, Particle *Partsubset, Int_t &numgroups, Int_t sublevel, Int_t *pnumcores)
{
    KDTree *tree;
    Int_t *pfof, i, ii;
    FOFcompfunc fofcmp;
    fstream Fout;
    char fname[200];
    Double_t param[20];
    int nsearch=opt.Nvel;
    Int_t **nnID;
    Double_t **dist2;
    int nthreads=1,maxnthreads,tid;
    Int_t *numingroup, **pglist;
    Int_tree_t *GroupTail, *Head, *Next;
    Int_t bgoffset, *pfofbg, numgroupsbg=0;
    int maxhalocoresublevel;
    Int_t numsubs=0;
    //initialize
    numgroups=0;
    if (pnumcores!=NULL) *pnumcores=0;
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif

#ifdef USEMPI
    //if using MPI on a single halo, lower minimum number
    if (opt.iSingleHalo) opt.MinSize=MinNumMPI;
#endif

    //set the maximum sublevel for halo core search. Unless star particles really should only be doing this at first sublevel
    maxhalocoresublevel=opt.maxnlevelcoresearch;
    if (opt.partsearchtype==PSTSTAR) maxhalocoresublevel=100;

    int minsize=opt.MinSize;

    //for parallel environment store maximum number of threads
#ifdef USEOPENMP
#pragma omp parallel
    {
    if (omp_get_thread_num()==0) maxnthreads=omp_get_num_threads();
    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
#endif
    //set parameters
    param[1]=(opt.ellxscale*opt.ellxscale)*(opt.ellphys*opt.ellphys);
    param[2]=(opt.ellvscale*opt.ellvscale)*(opt.ellvel*opt.ellvel);
    param[6]=(opt.ellxscale*opt.ellxscale)*(opt.ellphys*opt.ellphys);
    param[7]=(opt.Vratio);
    if (opt.foftype==FOF6DSUBSET) {
        param[2] = opt.HaloSigmaV*(opt.halocorevfac * opt.halocorevfac);
        param[7] = param[2];
    }
    param[8]=cos(opt.thetaopen*M_PI);
    param[9]=opt.ellthreshold;
    //if iterating slightly increase constraints and decrease minimum number
    if (opt.iiterflag && opt.foftype==FOFSTPROB) {
        if (opt.iverbose>=2) cout<<"Increasing thresholds to search for initial list.\n";
        if (opt.foftype==FOF6DSUBSET) param[7]*=opt.vfac*opt.vfac;
        else param[7]*=opt.vfac;
        param[8]=cos(opt.thetaopen*M_PI*opt.thetafac);
        param[9]=opt.ellthreshold*opt.ellfac;
        minsize*=opt.nminfac;
    }

    //Set fof type
    //@{
    if (opt.foftype==FOFSTPROB) {
        if (opt.iverbose>=2) {
        cout<<"FOFSTPROB which uses: ellphys, vratio, thetaopen, and ellthreshold.\n";
        cout<<"Parameters used are : ellphys="<<sqrt(param[6])<<" Lunits, vratio="<<param[7]<<", cos(thetaopen)="<<param[8]<<", ";
        cout<<param[9]<<" is outlier threshold.\n";
        }
        fofcmp=&FOFStreamwithprob;
    }
    else if (opt.foftype==FOF6DSUBSET) {
        if (opt.iverbose>=2) {
        cout<<"FOF6D uses ellphys and ellvel.\n";
        cout<<"Parameters used are : ellphys="<<sqrt(param[6])<<" Lunits, ellvel="<<sqrt(param[7])<<" Vunits.\n";
        }
        fofcmp=&FOF6d;
    }
    else if (opt.foftype==FOFSTPROBNN) {
        if (opt.iverbose>=2) {
        cout<<"FOFSTPROBNN which uses: ellphys, vratio, thetaopen, and ellthreshold.\n";
        cout<<"Parameters used are : ellphys="<<sqrt(param[6])<<" Lunits, vratio="<<param[7]<<", cos(thetaopen)="<<param[8]<<", ";
        cout<<param[9]<<" is outlier threshold.\n";
        cout<<"Search is limited to nearest neighbours."<<endl;
        }
        fofcmp=&FOFStreamwithprobNN;
    }
    else if (opt.foftype==FOFSTPROBNNLX) {
        if (opt.iverbose>=2) {
        cout<<"FOFSTPROBNNLX which uses: ellphys, vratio, thetaopen, and ellthreshold.\n";
        cout<<"Parameters used are : ellphys="<<sqrt(param[6])<<" Lunits, vratio="<<param[7]<<", cos(thetaopen)="<<param[8]<<", ";
        cout<<param[9]<<" is outlier threshold.\n";
        cout<<"Search is limited to nearest neighbours."<<endl;
        }
        fofcmp=&FOFStreamwithprobNNLX;
    }
    else if (opt.foftype==FOFSTPROBNNNODIST) {
        if (opt.iverbose>=2) {
        cout<<"FOFSTPROBNN which uses: vratio, thetaopen, and ellthreshold.\n";
        cout<<"Parameters used are : vratio="<<param[7]<<", cos(thetaopen)="<<param[8]<<", ";
        cout<<param[9]<<" is outlier threshold.\n";
        cout<<"Search is limited to nearest neighbours."<<endl;
        }
        fofcmp=&FOFStreamwithprobNNNODIST;
    }
    else if (opt.foftype==FOF6DCORE) {
        if (opt.iverbose>=2) {
        cout<<"FOF6DCORE which identifies phase-space dense regions and assigns particles, ie core identification and growth\n";
        }
        //just build tree and initialize the pfof array
        tree=new KDTree(Partsubset,nsubset,opt.Bsize,tree->TPHYS);
        numgroups=0;
        pfof=new Int_t[nsubset];
        for (i=0;i<nsubset;i++) pfof[i]=0;
    }
    //@}
    //now actually search for dynamically distinct substructures
    //@{
    if (!(opt.foftype==FOFSTPROBNN||opt.foftype==FOFSTPROBNNLX||opt.foftype==FOFSTPROBNNNODIST||opt.foftype==FOF6DCORE)) {
        if (opt.iverbose>=2) cout<<"Building tree ... "<<endl;
        tree=new KDTree(Partsubset,nsubset,opt.Bsize,tree->TPHYS);
        param[0]=tree->GetTreeType();
        //if large enough for statistically significant structures to be found then search. This is a robust search
        if (nsubset>=MINSUBSIZE) {
            if (opt.iverbose>=2) cout<<"Now search ... "<<endl;
            pfof=tree->FOFCriterion(fofcmp,param,numgroups,minsize,1,1,FOFchecksub);
        }
        else {
            numgroups=0;
            pfof=new Int_t[nsubset];
            for (i=0;i<nsubset;i++) pfof[i]=0;
        }
        if (opt.iverbose>=2) cout<<"Done"<<endl;
    }
    else if (opt.foftype==FOFSTPROBNN||opt.foftype==FOFSTPROBNNLX||opt.foftype==FOFSTPROBNNNODIST) {
        //here idea is to use subset but only search NN neighbours in phase-space once that is built for each particle, look at first particle's NN and see if any meet stffof criteria for velocity
        //then examine first tagged particle that meets critera by examining its NN and so on till reach particle where all NN are either already tagged or do not meet criteria
        //delete tree;
        if (opt.iverbose>=2) cout<<"Building tree ... "<<endl;
        tree=new KDTree(Partsubset,nsubset,opt.Bsize,tree->TPHYS,tree->KEPAN,1000,1);
        if (opt.iverbose>=2) cout<<"Finding nearest neighbours"<<endl;
        nnID=new Int_t*[nsubset];
        for (i=0;i<nsubset;i++) nnID[i]=new Int_t[nsearch];
        dist2=new Double_t*[nthreads];
        for (int j=0;j<nthreads;j++)dist2[j]=new Double_t[nsearch];
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,tid)
{
#pragma omp for
#endif
        for (i=0;i<nsubset;i++) {
#ifdef USEOPENMP
            tid=omp_get_thread_num();
#else
            tid=0;
#endif
            tree->FindNearest(i,nnID[i],dist2[tid],nsearch);
        }
#ifdef USEOPENMP
}
#endif
        if (opt.iverbose>=2) cout<<"Done"<<endl;
        if (opt.iverbose>=2) cout<<"search nearest neighbours"<<endl;
        pfof=tree->FOFNNCriterion(fofcmp,param,nsearch,nnID,numgroups,minsize);
        for (i=0;i<nsubset;i++) delete[] nnID[i];
        delete[] nnID;
        for (i=0;i<nthreads;i++) delete[] dist2[i];
        delete[] dist2;
        if (opt.iverbose>=2) cout<<"Done"<<endl;
    }
    //@}

    //iteration to search region around streams using lower thresholds
    //determine number of groups
    if (opt.iiterflag&&numgroups>0 && opt.foftype!=FOF6DCORE) {
        Int_t ng=numgroups;
        int mergers;
        int *igflag,*ilflag;
        Int_t newlinks, intergrouplinks,*newlinksIndex,*intergrouplinksIndex;
        Int_t pid,ppid, ss,tail,startpoint;
        Int_t numgrouplinks,*numgrouplinksIndex,**newIndex,*oldnumingroup,**newintergroupIndex,**intergroupgidIndex;

        //to slowly expand search, must declare several arrays making it easier to move along a group.
        numingroup=BuildNumInGroup(nsubset, numgroups, pfof);
        //here the lists are index order, not id order
        pglist=BuildPGList(nsubset, numgroups, numingroup, pfof, Partsubset);
        Head=BuildHeadArray(nsubset,numgroups,numingroup,pglist);
        Next=BuildNextArray(nsubset,numgroups,numingroup,pglist);
        GroupTail=BuildGroupTailArray(nsubset,numgroups,numingroup,pglist);
        newlinksIndex=new Int_t[nsubset];
        intergrouplinksIndex=new Int_t[numgroups+1];
        igflag=new int[numgroups+1];
        ilflag=new int[numgroups+1];
        for (i=1;i<=numgroups;i++) igflag[i]=0;
        // NOTE for openmp, only significantly faster if nnID is 1 D array for reduction purposes
        // However, this means that in this MPI domain, only nthreads < 2^32 / nsubset can be used since larger numbers require 64 bit array addressing
        // This should not be an issue but a check is placed anyway
#ifdef USEOPENMP
#pragma omp parallel
        {
            if (omp_get_thread_num()==0) maxnthreads=omp_get_num_threads();
            if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
        }
        if (log((double)nthreads*nsubset)/log(2.0)>32.0) nthreads=(max(1,(int)(pow((double)2,32)/(double)nsubset)));
        omp_set_num_threads(nthreads);
#endif
        nnID=new Int_t*[1];
        nnID[0]=new Int_t[nsubset*nthreads];

        numgrouplinksIndex=new Int_t[numgroups+1];
        oldnumingroup=new Int_t[numgroups+1];
        newIndex=new Int_t*[numgroups+1];

        //first search groups near cell size as subhaloes near this scale may have central regions defining background and as such these particles will not appear to be
        //distinct outliers (will probably have ell~1-2) thus resulting in an underestimate in the mass and maximum circular velocity for compact massive (sub)halos
        //this linking uses a smaller physical linking length and searchs already tagged particles for any untagged particles included in the larger subset made with the lower ellthreshold
        //so long as one particle lies about the high ell threshold, the particles can be check to see if they lie in the same phase-space.
        param[1]=(opt.ellxscale*opt.ellxscale)*(opt.ellphys*opt.ellphys);
        param[2]=(opt.ellvscale*opt.ellvscale)*(opt.ellvel*opt.ellvel);
        param[6]=(opt.ellxscale*opt.ellxscale)*(opt.ellphys*opt.ellphys);
        param[7]=(opt.Vratio);
        param[8]=cos(opt.thetaopen*M_PI);
        param[9]=opt.ellthreshold;
        //if (opt.foftype==FOF6DSUBSET) param[7]/=opt.vfac*opt.vfac;

        fofcmp=&FOFStreamwithprobIterative;
        if (opt.iverbose>=2) {
        cout<<ThisTask<<" "<<"Begin iterative search"<<endl;
        cout<<ThisTask<<" "<<"Parameters used are : ellphys="<<sqrt(param[6])<<" Lunits, vratio="<<param[7]<<", cos(thetaopen)="<<param[8]<<", ";
        cout<<ThisTask<<" "<<param[9]<<" is outlier threshold.\n";
        }

        //first identify all particles to be searched
        newlinks=0;
        for (i=1;i<=numgroups;i++) if (numingroup[i]>opt.Ncell*0.1&&igflag[i]==0)
        {
            ss=Head[pglist[i][0]];
            do {newlinksIndex[newlinks++]=ss;} while((ss = Next[ss]) >= 0);
        }
        //initialize nnID list so that if particle tagged no need to check comparison function
        for (ii=0;ii<nsubset;ii++) nnID[0][ii]=0;
        for (ii=0;ii<newlinks;ii++) nnID[0][Partsubset[newlinksIndex[ii]].GetID()]=pfof[Partsubset[newlinksIndex[ii]].GetID()];
        //first search list and find any new links, not that here if a particle has a pfof value > that the pfofvalue of the reference particle, its nnID value is set to the reference
        //particle's group id, otherwise, its left alone (unless its untagged)
        //if number of searches is large, run search in parallel
        SearchForNewLinks(nsubset, tree, Partsubset, pfof, fofcmp, param, newlinks, newlinksIndex, nnID, nthreads);
        tid=0;
        //determine groups
        //newIndex=DetermineNewLinks(nsubset, Partsubset, pfof, numgroups, newlinks, newlinksIndex, numgrouplinksIndex, nnID[tid]);
        DetermineNewLinks(nsubset, Partsubset, pfof, numgroups, newlinks, newlinksIndex, numgrouplinksIndex, nnID[tid],&newIndex);
        //link all previously unlinked particles for large groups by searching for each group the number of links for that group and linking them (if unlinked)
        LinkUntagged(Partsubset, numgroups, pfof, numingroup, pglist, newlinks, numgrouplinksIndex, newIndex, Head, Next, GroupTail,nnID[tid]);
        for (Int_t j=1;j<=numgroups;j++) delete[] newIndex[j];

        for (i=0;i<nsubset;i++) if (Partsubset[i].GetPotential()<param[9]&&nnID[0][Partsubset[i].GetID()]==0) nnID[0][Partsubset[i].GetID()]=-1;

        fofcmp=&FOFStreamwithprob;
        param[1]=(opt.ellxscale*opt.ellxscale)*(opt.ellphys*opt.ellphys);// *opt.ellxfac*opt.ellxfac;//increase physical linking length slightly
        param[6]=(opt.ellxscale*opt.ellxscale)*(opt.ellphys*opt.ellphys);// *opt.ellxfac*opt.ellxfac;
        param[7]=(opt.Vratio)*opt.vfac;
        param[8]=cos(opt.thetaopen*M_PI*opt.thetafac);
        param[9]=opt.ellthreshold*opt.ellfac;
        if (opt.iverbose>=2) {
        cout<<ThisTask<<" "<<"Begin second expanded search with large linking length "<<endl;
        cout<<ThisTask<<" "<<"FOFSTPROB which uses: ellphys, vratio, thetaopen, and ellthreshold.\n";
        cout<<ThisTask<<" "<<"Parameters used are : ellphys="<<sqrt(param[6])<<" Lunits, vratio="<<param[7]<<", cos(thetaopen)="<<param[8]<<", ";
        cout<<ThisTask<<" "<<param[9]<<" is outlier threshold.\n";
        }
        //initialization, for merging, store old size prior to expanded search, set nnID
        newlinks=0;
        for (i=1;i<=numgroups;i++)
        {
            oldnumingroup[i]=numingroup[i];
            startpoint=Head[pglist[i][0]];
            ss=startpoint;
            do {if (Partsubset[ss].GetPotential()>=param[9]){newlinksIndex[newlinks++]=ss;}} while((ss = Next[ss]) >= 0);
        }
        //now continue to search all new links till there are no more links found
        tid=0;
        do {
            SearchForNewLinks(nsubset, tree, Partsubset, pfof, fofcmp, param, newlinks, newlinksIndex, nnID, nthreads);
            DetermineNewLinks(nsubset, Partsubset, pfof, numgroups, newlinks, newlinksIndex, numgrouplinksIndex, nnID[tid],&newIndex);
            LinkUntagged(Partsubset, numgroups, pfof, numingroup, pglist, newlinks, numgrouplinksIndex, newIndex, Head, Next, GroupTail,nnID[tid]);
            for (Int_t j=1;j<=numgroups;j++) delete[] newIndex[j];
        }while(newlinks);

        //then search for intergroup links. When merging, identify all intergroup links and define merging criterion based on mass ratios, number of links, old mass (prior to expansion), etc
        if (opt.iverbose>=2) cout<<ThisTask<<" "<<"Search for intergroup links"<<endl;
        newintergroupIndex=new Int_t*[numgroups+1];
        intergroupgidIndex=new Int_t*[numgroups+1];
        newlinks=0;
        tid=0;
        for (i=1;i<=numgroups;i++) if (igflag[i]==0)
        {
            startpoint=Head[pglist[i][0]];
            ss=startpoint;
            do {if (Partsubset[ss].GetDensity()>=param[9]){newlinksIndex[newlinks++]=ss;}} while((ss = Next[ss]) >= 0);
        }
        for (ii=0;ii<newlinks;ii++) nnID[0][Partsubset[newlinksIndex[ii]].GetID()]=pfof[Partsubset[newlinksIndex[ii]].GetID()];
        do {
            //initialize
            mergers=0;
            //first search list and find any new links, not that here if a particle has a pfof value > that the pfofvalue of the reference particle, its nnID value is set to the reference
            //particle's group id, otherwise, its left alone (unless its untagged)
            //if number of searches is large, run search in parallel
            SearchForNewLinks(nsubset, tree, Partsubset, pfof, fofcmp, param, newlinks, newlinksIndex, nnID, nthreads);
            DetermineGroupLinks(nsubset, Partsubset, pfof, numgroups, newlinks, newlinksIndex, numgrouplinksIndex, nnID[tid], &newIndex);
            DetermineGroupMergerConnections(Partsubset, numgroups, pfof, ilflag, numgrouplinksIndex, intergrouplinksIndex, nnID[tid], &newIndex, &newintergroupIndex, &intergroupgidIndex);
            newlinks=0;
            mergers+=MergeGroups(opt, Partsubset, numgroups, pfof, numingroup, oldnumingroup, pglist, numgrouplinksIndex, &intergroupgidIndex, &newintergroupIndex, intergrouplinksIndex, Head, Next, GroupTail, igflag, nnID[tid], newlinks, newlinksIndex);
        }while(mergers);

        param[1]=(opt.ellxscale*opt.ellxscale)*(opt.ellphys*opt.ellphys)*opt.ellxfac*opt.ellxfac;//increase physical linking length slightly again
        param[6]=(opt.ellxscale*opt.ellxscale)*(opt.ellphys*opt.ellphys)*opt.ellxfac*opt.ellxfac;
        param[7]=(opt.Vratio)*opt.vfac;
        param[8]=cos(opt.thetaopen*M_PI*opt.thetafac);
        param[9]=opt.ellthreshold*opt.ellfac;
        if (opt.iverbose>=2){
        cout<<ThisTask<<" "<<"Begin second expanded search with large linking length "<<endl;
        cout<<ThisTask<<" "<<"FOFSTPROB which uses: ellphys, vratio, thetaopen, and ellthreshold.\n";
        cout<<ThisTask<<" "<<"Parameters used are : ellphys="<<sqrt(param[6])<<" Lunits, vratio="<<param[7]<<", cos(thetaopen)="<<param[8]<<", ";
        cout<<ThisTask<<" "<<param[9]<<" is outlier threshold."<<endl;
        }
        //initialization, for merging, store old size prior to expanded search, set nnID
        newlinks=0;
        for (i=1;i<=numgroups;i++) if (igflag[i]==0)
        {
            startpoint=Head[pglist[i][0]];
            ss=startpoint;
            do {if (Partsubset[ss].GetPotential()>=param[9]){newlinksIndex[newlinks++]=ss;}} while((ss = Next[ss]) >= 0);
        }
        //now continue to search all new links till there are no more links found
        tid=0;
        do {
            SearchForNewLinks(nsubset, tree, Partsubset, pfof, fofcmp, param, newlinks, newlinksIndex, nnID, nthreads);
            DetermineNewLinks(nsubset, Partsubset, pfof, numgroups, newlinks, newlinksIndex, numgrouplinksIndex, nnID[tid],&newIndex);
            LinkUntagged(Partsubset, numgroups, pfof, numingroup, pglist, newlinks, numgrouplinksIndex, newIndex, Head, Next, GroupTail,nnID[tid]);
            for (Int_t j=1;j<=numgroups;j++) delete[] newIndex[j];
        }while(newlinks);

        //release memory
        for (i=1;i<=numgroups;i++)delete[] pglist[i];
        delete[] Head;
        delete[] Next;
        delete[] GroupTail;
        delete[] newlinksIndex;
        delete[] numgrouplinksIndex;
        delete[] newIndex;
        delete[] igflag;
        delete[] nnID[0];
        delete[] nnID;
        delete[] intergrouplinksIndex;
        delete[] intergroupgidIndex;
        delete[] newintergroupIndex;
        delete[] ilflag;

#ifdef USEOPENMP
#pragma omp parallel
    {
    omp_set_num_threads(maxnthreads);
    }
    nthreads=maxnthreads;
#endif
        //adjust groups
        for (i=1;i<=numgroups;i++) numingroup[i]=0;
        for (i=0;i<nsubset;i++) numingroup[pfof[i]]++;
        for (i=0;i<nsubset;i++) if (numingroup[pfof[i]]<opt.MinSize) pfof[i]=0;
        for (i=1;i<=numgroups;i++) numingroup[i]=0;
        for (i=0;i<nsubset;i++) numingroup[pfof[i]]++;
        //cout<<ThisTask<<" "<<"Now determine number of groups with non zero length"<<endl;
        for (i=numgroups;i>=1;i--) {
            if (numingroup[i]==0) ng--;
            else pglist[i]=new Int_t[numingroup[i]];
            numingroup[i]=0;
        }
        for (i=0;i<nsubset;i++) if (pfof[i]>0) pglist[pfof[i]][numingroup[pfof[i]]++]=i;
        if (ng) ReorderGroupIDs(numgroups, ng, numingroup, pfof, pglist);
        for (i=1;i<=numgroups;i++) if (numingroup[i]>0) delete[] pglist[i];
        delete[] pglist;
        delete[] numingroup;
        numgroups=ng;
        if (opt.iverbose>=2) cout<<ThisTask<<" "<<"After expanded search there are now "<< ng<<" groups"<<endl;
    }
    //once substructure groups are found, ensure all substructure groups are significant
    if (numgroups)
    {
        numingroup=BuildNumInGroup(nsubset, numgroups, pfof);
        //pglist is constructed without assuming particles are in index order
        pglist=BuildPGList(nsubset, numgroups, numingroup, pfof, Partsubset);
        CheckSignificance(opt,nsubset,Partsubset,numgroups,numingroup,pfof,pglist);
        for (i=1;i<=numgroups;i++) delete[] pglist[i];
        delete[] pglist;
        delete[] numingroup;
    }
    if (numgroups>0) if (opt.iverbose>=2) cout<<ThisTask<<": "<<numgroups<<" substructures found"<<endl;
    else {if (opt.iverbose>=2) cout<<ThisTask<<": "<<"NO SUBSTRUCTURES FOUND"<<endl;}

    //now search particle list for large compact substructures that are considered part of the background when using smaller grids
    if (nsubset>=MINSUBSIZE && opt.iLargerCellSearch && opt.foftype!=FOF6DCORE)
    {
        //first have to delete tree used in search so that particles are in original particle order
        //then construct a new grid with much larger cells so that new bg velocity dispersion can be estimated
        if (opt.iverbose>=2) cout<<" entering large cell search "<<opt.iLargerCellSearch<<endl;
        delete tree;
        Int_t ngrid;
        Coordinate *gvel;
        Matrix *gveldisp;
        GridCell *grid;
        Double_t nf, ncl=opt.Ncell;
        //adjust ncellfac locally
        nf=(opt.Ncellfac*8.0,MAXCELLFRACTION);
        opt.Ncell=nf*nsubset;

        //ONLY calculate grid quantities if substructures have been found
        if (numgroups>0) {
            tree=InitializeTreeGrid(opt,nsubset,Partsubset);
            ngrid=tree->GetNumLeafNodes();
            if (opt.iverbose>=2) cout<<ThisTask<<" "<<"bg search using "<<ngrid<<" grid cells, with each node containing ~"<<(opt.Ncell=nsubset/ngrid)<<" particles"<<endl;
            grid=new GridCell[ngrid];
            FillTreeGrid(opt, nsubset, ngrid, tree, Partsubset, grid);
            gvel=GetCellVel(opt,nsubset,Partsubset,ngrid,grid);
            gveldisp=GetCellVelDisp(opt,nsubset,Partsubset,ngrid,grid,gvel);
            GetDenVRatio(opt,nsubset,Partsubset,ngrid,grid,gvel,gveldisp);
            GetOutliersValues(opt,nsubset,Partsubset,-1);
        }
        ///produce tree to search for 6d phase space structures
        tree=new KDTree(Partsubset,nsubset,opt.Bsize,tree->TPHYS);

        //now begin fof6d search for large background objects that are missed using smaller grid cells ONLY IF substructures have been found
        //this search can identify merger excited radial shells so for the moment, disabled

        if (numgroups>0) {
            bgoffset=0;
            minsize=ncl*0.2;
            //use same linking length to find substructures
            param[1]=(opt.ellxscale*opt.ellxscale)*(opt.ellphys*opt.ellphys)*(opt.ellxfac*opt.ellxfac);
            param[6]=param[1];
            //velocity linking length from average sigmav from grid
            param[7]=opt.HaloLocalSigmaV;
            param[8]=cos(opt.thetaopen*M_PI);
            param[9]=opt.ellthreshold*opt.ellfac;
            if (opt.iverbose>=2) {
            cout<<ThisTask<<" "<<"FOF6D uses ellphys and ellvel."<<endl;
            cout<<ThisTask<<" "<<"Parameters used are : ellphys="<<sqrt(param[6])<<" Lunits, ellvel="<<sqrt(param[7])<<" Vunits."<<endl;
            cout<<ThisTask<<" "<<param[9]<<" is outlier threshold and contains more than "<<minsize<<endl;
            }
            //then reset to find large subhalo cores close to cell size limit
            fofcmp=&FOF6dbgup;
            //here this ensures that particles belong to a 6dfof substructure that is composed of particles
            //which are considered dynamical outliers using a very large grid
            pfofbg=tree->FOFCriterion(fofcmp,param,numgroupsbg,minsize,1,1,FOFchecksub);

            //now combine results such that if particle already belongs to a substructure ignore
            //otherwise leave tagged. Then check if these new background substructures share enough links with the
            //other substructures that they should be joined.
            if (numgroupsbg>=bgoffset+1) {
                int mergers;
                int *igflag,*ilflag;
                Int_t newlinks, oldlinks,intergrouplinks,*newlinksIndex,*oldlinksIndex,*intergrouplinksIndex;
                Int_t pid,ppid, ss,tail,startpoint;
                Int_t numgrouplinks,*numgrouplinksIndex,**newIndex,**newintergroupIndex,**intergroupgidIndex;

                Int_t ng=numgroups+(numgroupsbg-bgoffset),oldng=numgroups;
                for (i=0;i<nsubset;i++) if (pfof[i]==bgoffset&&pfofbg[i]>0) pfof[i]=oldng+(pfofbg[i]-bgoffset);
                numgroups+=(numgroupsbg-bgoffset);
                if (opt.iverbose>=2) cout<<ThisTask<<" "<<"Found "<<numgroups<<endl;

                //
                numingroup=BuildNumInGroup(nsubset, numgroups, pfof);
                pglist=BuildPGList(nsubset, numgroups, numingroup, pfof, Partsubset);
                Head=BuildHeadArray(nsubset,numgroups,numingroup,pglist);
                Next=BuildNextArray(nsubset,numgroups,numingroup,pglist);
                GroupTail=BuildGroupTailArray(nsubset,numgroups,numingroup,pglist);
                newlinksIndex=new Int_t[nsubset];
                oldlinksIndex=new Int_t[nsubset];
                intergrouplinksIndex=new Int_t[numgroups+1];
                igflag=new int[numgroups+1];
                ilflag=new int[numgroups+1];
                for (i=1;i<=numgroups;i++) igflag[i]=0;
#ifdef USEOPENMP
#pragma omp parallel
                {
                    if (omp_get_thread_num()==0) maxnthreads=omp_get_num_threads();
                    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
                }
                if (log((double)nthreads*nsubset)/log(2.0)>32.0) nthreads=(max(1,(int)(pow((double)2,32)/(double)nsubset)));
                omp_set_num_threads(nthreads);
#endif

                nnID=new Int_t*[1];
                nnID[0]=new Int_t[nsubset*nthreads];

                numgrouplinksIndex=new Int_t[numgroups+1];
                newIndex=new Int_t*[numgroups+1];

                //initialize nnID list so that if particle tagged no need to check comparison function
                for (ii=0;ii<nsubset;ii++) nnID[0][ii]=0;
                if (opt.iverbose>=2) cout<<ThisTask<<" "<<"Search for bg intergroup links"<<endl;
                newintergroupIndex=new Int_t*[numgroups+1];
                intergroupgidIndex=new Int_t*[numgroups+1];
                oldlinks=0;
                tid=0;
                for (i=1;i<=oldng;i++) {
                    startpoint=Head[pglist[i][0]];
                    ss=startpoint;
                    do {{oldlinksIndex[oldlinks++]=ss;}} while((ss = Next[ss]) >= 0);
                }
                for (ii=0;ii<oldlinks;ii++)
                    nnID[0][Partsubset[oldlinksIndex[ii]].GetID()]=pfof[Partsubset[oldlinksIndex[ii]].GetID()]+numgroups;
                newlinks=0;
                for (i=oldng+1;i<=numgroups;i++) if (numingroup[i]>0) {
                    startpoint=Head[pglist[i][0]];
                    ss=startpoint;
                    do {{newlinksIndex[newlinks++]=ss;}} while((ss = Next[ss]) >= 0);
                }
                else igflag[i]=1;
                for (ii=0;ii<newlinks;ii++)
                    nnID[0][Partsubset[newlinksIndex[ii]].GetID()]=pfof[Partsubset[newlinksIndex[ii]].GetID()];

                //now begin interative
                param[1]=(opt.ellxscale*opt.ellxscale)*(opt.ellphys*opt.ellphys);
                param[2]=(opt.ellvscale*opt.ellvscale)*(opt.ellvel*opt.ellvel);
                param[6]=(opt.ellxscale*opt.ellxscale)*(opt.ellphys*opt.ellphys);
                param[7]=(opt.Vratio);
                param[8]=cos(opt.thetaopen*M_PI);
                param[9]=opt.ellthreshold*opt.ellfac*0.8;
                fofcmp=&FOFStreamwithprobIterative;
                if (opt.iverbose>=2){
                cout<<ThisTask<<" "<<"Begin expanded search for groups near cell size"<<endl;
                cout<<ThisTask<<" "<<"Parameters used are : ellphys="<<sqrt(param[6])<<" Lunits, vratio="<<param[7]<<", cos(thetaopen)="<<param[8]<<", ";
                cout<<ThisTask<<" "<<param[9]<<" is outlier threshold.\n";
                }
                SearchForNewLinks(nsubset, tree, Partsubset, pfof, fofcmp, param, newlinks, newlinksIndex, nnID, nthreads);
                tid=0;
                DetermineNewLinks(nsubset, Partsubset, pfof, numgroups, newlinks, newlinksIndex, numgrouplinksIndex, nnID[tid],&newIndex);
                LinkUntagged(Partsubset, numgroups, pfof, numingroup, pglist, newlinks, numgrouplinksIndex, newIndex, Head, Next, GroupTail,nnID[tid]);
                for (Int_t j=1;j<=numgroups;j++) delete[] newIndex[j];
                newlinks=0;
                for (i=oldng+1;i<=numgroups;i++) if (numingroup[i]>0) {
                    startpoint=Head[pglist[i][0]];
                    ss=startpoint;
                    do {{newlinksIndex[newlinks++]=ss;}} while((ss = Next[ss]) >= 0);
                }
                else igflag[i]=1;

                fofcmp=&FOFStreamwithprob;
                param[1]=(opt.ellxscale*opt.ellxscale)*(opt.ellphys*opt.ellphys)*opt.ellxfac*opt.ellxfac;//increase physical linking length slightly
                param[6]=(opt.ellxscale*opt.ellxscale)*(opt.ellphys*opt.ellphys)*opt.ellxfac*opt.ellxfac;
                param[7]=(opt.Vratio)*opt.vfac;
                param[8]=cos(opt.thetaopen*M_PI*opt.thetafac);
                param[9]=-3.0;

                mergers=0;
                SearchForNewLinks(nsubset, tree, Partsubset, pfof, fofcmp, param, newlinks, newlinksIndex, nnID, nthreads);
                for (ii=0;ii<oldlinks;ii++)
                    if (nnID[0][Partsubset[oldlinksIndex[ii]].GetID()]==pfof[Partsubset[oldlinksIndex[ii]].GetID()]+numgroups)
                        nnID[0][Partsubset[oldlinksIndex[ii]].GetID()]-=numgroups;
                DetermineGroupLinks(nsubset, Partsubset, pfof, numgroups, oldlinks, oldlinksIndex, numgrouplinksIndex, nnID[tid], &newIndex);
                DetermineGroupMergerConnections(Partsubset, numgroups, pfof, ilflag, numgrouplinksIndex, intergrouplinksIndex, nnID[tid], &newIndex, &newintergroupIndex, &intergroupgidIndex);
                mergers+=MergeGroups(opt, Partsubset, numgroups, pfof, numingroup, numingroup, pglist, numgrouplinksIndex, &intergroupgidIndex, &newintergroupIndex, intergrouplinksIndex, Head, Next, GroupTail, igflag, nnID[tid], newlinks, newlinksIndex);

                //release memory
                for (i=1;i<=numgroups;i++)delete[] pglist[i];
                delete[] Head;
                delete[] Next;
                delete[] GroupTail;
                delete[] newlinksIndex;
                delete[] numgrouplinksIndex;
                delete[] newIndex;
                delete[] igflag;
                delete[] nnID[0];
                delete[] nnID;
                delete[] intergrouplinksIndex;
                delete[] intergroupgidIndex;
                delete[] newintergroupIndex;
                delete[] ilflag;

#ifdef USEOPENMP
#pragma omp parallel
    {
    omp_set_num_threads(maxnthreads);
    }
    nthreads=maxnthreads;
#endif

                //adjust sizes and reorder groups
                for (i=1;i<=numgroups;i++) numingroup[i]=0;
                for (i=0;i<nsubset;i++) numingroup[pfof[i]]++;
                for (i=0;i<nsubset;i++) if (numingroup[pfof[i]]<opt.MinSize) pfof[i]=0;
                for (i=1;i<=numgroups;i++) numingroup[i]=0;
                for (i=0;i<nsubset;i++) numingroup[pfof[i]]++;
                if (opt.iverbose>=2) cout<<ThisTask<<" "<<"Now determine number of groups with non zero length"<<endl;
                for (i=numgroups;i>=1;i--) {
                    if (numingroup[i]==0) ng--;
                    else pglist[i]=new Int_t[numingroup[i]];
                    numingroup[i]=0;
                }
                for (i=0;i<nsubset;i++) if (pfof[i]>0) pglist[pfof[i]][numingroup[pfof[i]]++]=i;
                if (ng) ReorderGroupIDs(numgroups, ng, numingroup, pfof, pglist);
                for (i=1;i<=numgroups;i++) if (numingroup[i]>0) delete[] pglist[i];
                delete[] pglist;
                delete[] numingroup;
                numgroups=ng;
                if (opt.iverbose>=2) cout<<ThisTask<<" "<<"After expanded search there are now "<< ng<<" groups"<<endl;
                numgroupsbg=0;
            }
            else if (opt.iverbose>=2) cout<<ThisTask<<" "<<"No large background substructure groups found"<<endl;
            delete[] pfofbg;
        }
        //output results of search
        if (numgroups>0) if (opt.iverbose>=2) cout<<numgroups<<" substructures found after large grid search"<<endl;
        else {if (opt.iverbose>=2) cout<<ThisTask<<": "<<"NO SUBSTRUCTURES FOUND"<<endl;}
    }
    numsubs=numgroups;

    //ONCE ALL substructures are found, search for cores of major mergers with minimum size set by cell size since grid is quite large after bg search
    //for missing large substructure cores
    if((opt.iHaloCoreSearch>0&&((!opt.iSingleHalo&&sublevel<=maxhalocoresublevel)||(opt.iSingleHalo&&sublevel==0)))||opt.foftype==FOF6DCORE)
    {
        if (opt.iverbose>=2) cout<<ThisTask<<" beginning 6dfof core search to find multiple cores"<<endl;
        bgoffset=1;
        //if adaptive core linking then need to calculate dispersion tensors in configuration and velocity space
        if (opt.iAdaptiveCoreLinking)
        {
            double xscale2;
            Double_t XX, YY, ZZ;
            Double_t VOL;
            Double_t nn;

            //calculate inertia tensor (assumes that system is alread in CM frame
            ///\todo should add checks to see if system is in CM frame and if it is not, move it to that frame
            Matrix myeigvec, myI;
            CalcPosSigmaTensor(nsubset, Partsubset, XX, YY, ZZ, myeigvec, myI, -1);
            xscale2=sqrt(XX+YY+ZZ);
            //volume of configuration ellipsoid
            VOL = 4/3.0 * acos(-1) * sqrt(XX*XX*XX)*pow(opt.halocoresigmafac,3.0);
            //number density
            nn = pow(VOL/nsubset,1./3.);
            //now adjust linking length
            //param[1] = xscale2 * nn * nn * opt.halocorexfac * opt.halocorexfac;
            param[1] = nn*nn*opt.halocorexfac * opt.halocorexfac;
            param[6] = param[1];
            //now the same for velocity space
            ///\todo still need to look at optimal scaling function for velocities
            CalcVelSigmaTensor(nsubset, Partsubset, XX, YY, ZZ, myeigvec, myI, -1);
            VOL = 4/3.0 * acos(-1) * sqrt(XX*YY*ZZ);
            nn = pow(VOL/nsubset,1./3.);
            param[2] = XX*(opt.halocorevfac * opt.halocorevfac);
            param[7] = param[2];
        }
        //not adaptive, using halo based linking lengths and halo dispersion
        else
        {
            if ((opt.iKeepFOF) && (opt.ellhalo6dxfac > 0))
                param[1] = (opt.ellxscale * opt.ellphys * opt.ellhalo6dxfac * opt.halocorexfac);
            else
                param[1] = (opt.ellxscale * opt.ellphys * opt.ellhalophysfac * opt.halocorexfac);
            //if at a sublevel then start at a higher density threshold
            param[1] *= pow(opt.halocorexfac,sublevel-1);
            param[1] = param[1] * param[1];
            param[6] = param[1];
            //velocity linking length from average sigmav from FINE SCALE grid
            param[2] = opt.HaloSigmaV * (opt.halocorevfac * opt.halocorevfac);
            param[7] = param[2];
        }

        //set minsize of the core.
        //for a galaxy search (ie: PSTSTAR, we set this to the min search size,
        //but for others, smaller (DM) substructures have already been cleaned.
        if (opt.partsearchtype!=PSTSTAR) {
            minsize=nsubset*opt.halocorenfac*pow(opt.halocorenumfaciter,sublevel-1);
            minsize=max(minsize,opt.MinSize);
        }
        else if (opt.foftype==FOF6DCORE || opt.partsearchtype==PSTSTAR){
            minsize=opt.MinSize;
        }
        if (opt.iverbose>=2) {
            cout<<ThisTask<<" "<<"Parameters used are : ellphys="<<sqrt(param[6])<<" Lunits, ellvel="<<sqrt(param[7])<<" Vunits"<<endl;
            cout<<"with minimum size of "<<minsize<<endl;
        }
        //start first 6d search
        fofcmp=&FOF6d;
        //here search for 6dfof groups, return ordered list. Also pass function check to ignore already tagged particles
        int iorder=1,icheck=1,numloops=0,numactiveloops=0;
        //for ignoring already tagged particles can just use the oultier potential and FOFcheckbg
        for (i=0;i<nsubset;i++) Partsubset[i].SetPotential(pfof[Partsubset[i].GetID()]);
        for (i=0;i<nsubset;i++) Partsubset[i].SetType(-1);
        param[9]=0.5;
        pfofbg=tree->FOFCriterion(fofcmp,param,numgroupsbg,minsize,iorder,icheck,FOFcheckbg);

        for (i=0;i<nsubset;i++) if (pfofbg[Partsubset[i].GetID()]<=1 && pfof[Partsubset[i].GetID()]==0) Partsubset[i].SetType(numactiveloops);

        //store the dispersion limit factor, which depends on level of the loop at which cores are found as the disperions criterion limits how large a dispersion a core can have when measured
        vector<Double_t> dispfac(numgroupsbg+1);
        //store the level a core is found at
        vector<int> corelevel(numgroupsbg+1);

        //now if searching for cores in fully adaptive fashion then process core (which will keep changing) till
        //no cores are found
        for (i=1;i<=numgroupsbg;i++) dispfac[i]=1.0;
        for (i=1;i<=numgroupsbg;i++) corelevel[i]=numactiveloops;
        if (opt.halocorenumloops>1)
        {
            //store the old velocity dispersion
            Double_t halocoreveldisp = param[7];
            //set the factor by which the minimum size is increased by
            //store the old number of groups and update the pfofbg list
            Int_t newnumgroupsbg=numgroupsbg;
            Int_t *pfofbgnew=new Int_t[nsubset];
            Int_t pid;
            //first copy pfofbg information
            newnumgroupsbg=numgroupsbg;
            for (i=0;i<nsubset;i++) pfofbgnew[i]=pfofbg[i];
            //for (i=1;i<=newnumgroupsbg;i++) dispfac[i]=param[7];

            //now keep doing this till zero new groups are found, copying over info of new groups above bgoffset as we go
            //store how much the dispersion will change given the allowed fof envelop;
            Double_t dispval=opt.halocorevfaciter*opt.halocorevfaciter*opt.halocorexfaciter*opt.halocorexfaciter;
            Double_t dispvaltot=1.0;
            do {
                numloops++;
                //free memory
                delete[] pfofbg;
                param[1]*=opt.halocorexfaciter*opt.halocorexfaciter;
                param[6]=param[1];
               	//adjust linking lengths in velocity space by the iterative factor, ~0.75
                param[2]*=opt.halocorevfaciter*opt.halocorevfaciter;
                param[7]=param[2];
               	//and alter the minsize
                minsize*=opt.halocorenumfaciter;
                if (minsize<opt.MinSize) minsize=opt.MinSize;
                dispvaltot*=dispval;
                //we adjust the particles potentials so as to ignore already tagged particles using FOFcheckbg
                //here since loop just iterates to search the largest core, we just set all previously tagged particles not belonging to main core as 1
                for (i=0;i<nsubset;i++) Partsubset[i].SetPotential((pfofbgnew[Partsubset[i].GetID()]!=1)+(pfof[Partsubset[i].GetID()]>0));
                pfofbg=tree->FOFCriterion(fofcmp,param,numgroupsbg,minsize,iorder,icheck,FOFcheckbg);
                //now if numgroupsbg is greater than one, need to update the pfofbgnew array
                if (numgroupsbg>1) {
                    numactiveloops++;
                    for (i=0;i<nsubset;i++) if (pfofbg[Partsubset[i].GetID()]==1 && pfofbgnew[Partsubset[i].GetID()]<=1 && pfof[Partsubset[i].GetID()]==0) Partsubset[i].SetType(numactiveloops);
                    dispfac[1]=dispvaltot;
                    corelevel[1]=numactiveloops;
                    for (i=2;i<=numgroupsbg;i++) dispfac.push_back(dispvaltot);
                    for (i=2;i<=numgroupsbg;i++) corelevel.push_back(numactiveloops);
                    for (i=0;i<nsubset;i++) {
                        pid=Partsubset[i].GetID();
                        //if the particle is untagged with new search, set it to be untagged
                        if (pfofbgnew[pid]==1 && pfofbg[pid]==0) pfofbgnew[pid]=0;
                        //if particle is tagged as new group > 1 then keep but offset id
                        else if (pfofbg[pid]>1) pfofbgnew[pid]=pfofbg[pid]-1+newnumgroupsbg;
                    }
                    newnumgroupsbg+=numgroupsbg-1;

                }
            }while (numgroupsbg > 0 && numloops<opt.halocorenumloops && minsize*opt.halocorenumfaciter<nsubset);
            //once the loop is finished, update info
            numgroupsbg=newnumgroupsbg;
            for (i=0;i<nsubset;i++) pfofbg[i]=pfofbgnew[i];
            delete[] pfofbgnew;
            param[7]=halocoreveldisp;
        }

        if (numgroupsbg>=bgoffset+1) {
            if (opt.iverbose>=2) cout<<"Number of cores: "<<numgroupsbg<<endl;
            //if cores are found, two options
            //a pnumcores pointer is passed, store the number of cores in this pointer assumes halo ids will not be reordered, making it easy to flag cores
            //if opt.iHaloCoreSearch==2 then start searching the halo to identify which particles belong to which core by linking particles to their closest
            //tagged neighbour in phase-space.
            if(pnumcores!=NULL) *pnumcores=numgroupsbg;
            if (opt.iHaloCoreSearch>=2) {
                HaloCoreGrowth(opt, nsubset, Partsubset, pfof, pfofbg, numgroupsbg, param,dispfac,numactiveloops,corelevel,nthreads);
                if(pnumcores!=NULL) *pnumcores=numgroupsbg;
                if (numgroupsbg>=bgoffset+1) {
                    for (i=0;i<nsubset;i++) if (pfofbg[i]>bgoffset) pfof[i]=numgroups+(pfofbg[i]-bgoffset);
                    numgroupsbg-=bgoffset;
                    //if (numgroups>0 && opt.coresubmergemindist>0) MergeSubstructuresCoresPhase(opt, nsubset, Partsubset, pfof, numgroups, numgroupsbg);
                    numgroups+=numgroupsbg;
                }
                if (opt.iverbose>=2) cout<<ThisTask<<": After 6dfof core search and assignment there are "<<numgroups<<" groups"<<endl;
            }
        }
        else {
            if (opt.iverbose>=2) cout<<ThisTask<<": has found no excess cores indicating mergers"<<endl;
        }

        delete[] pfofbg;
    }
    if (numgroups>0 && opt.coresubmergemindist>0 && nsubset>=MINSUBSIZE) MergeSubstructuresPhase(opt, nsubset, Partsubset, pfof, numgroups, numsubs, numgroupsbg);


#ifdef USEMPI
    //now if substructures are subsubstructures, then the region of interest has already been localized to a single MPI domain
    //so do not need to search other mpi domains
    if (sublevel==0) {

    //if using MPI must determine which local particles need to be exported to other threads and used to search
    //that threads particles. This is done by seeing if the any particles have a search radius that overlaps with
    //the boundaries of another threads domain. Then once have exported particles must search local particles
    //relative to these exported particles. Iterate over search till no new links are found.

    //First have barrier to ensure that all mpi tasks have finished the local search
    MPI_Barrier(MPI_COMM_WORLD);

    //To make it easy to search the local particle arrays, build several arrays, in particular, Head, Next, Len arrays
    Int_tree_t *Len;
    numingroup=BuildNumInGroup(nsubset, numgroups, pfof);
    pglist=BuildPGList(nsubset, numgroups, numingroup, pfof,Partsubset);
    Len=BuildLenArray(nsubset,numgroups,numingroup,pglist);
    Head=BuildHeadArray(nsubset,numgroups,numingroup,pglist);
    Next=BuildNextArray(nsubset,numgroups,numingroup,pglist);
    for (i=1;i<=numgroups;i++)delete[] pglist[i];
    delete[] pglist;
    //Also must ensure that group ids do not overlap between mpi threads so adjust group ids
    MPI_Allgather(&numgroups, 1, MPI_Int_t, mpi_ngroups, 1, MPI_Int_t, MPI_COMM_WORLD);
    MPIAdjustLocalGroupIDs(nsubset, pfof);

    //then determine export particles, declare arrays used to export data
#ifdef MPIREDUCEMEM
#ifdef SWIFTINTERFACE
    MPIGetExportNumUsingMesh(libvelociraptorOpt, nbodies, Partsubset, sqrt(param[1]));
#else
    MPIGetExportNum(nbodies, Partsubset, sqrt(param[1]));
#endif
#endif
    //then determine export particles, declare arrays used to export data
    PartDataIn = new Particle[NExport];
    PartDataGet = new Particle[NExport];
    FoFDataIn = new fofdata_in[NExport];
    FoFDataGet = new fofdata_in[NExport];
    //I have adjusted FOF data structure to have local group length and also seperated the export particles from export fof data
    //the reason is that will have to update fof data in iterative section but don't need to update particle information.
#ifdef SWIFTINTERFACE
    MPIBuildParticleExportListUsingMesh(libvelociraptorOpt, nsubset, Partsubset, pfof, Len, sqrt(param[1]));
#else
    MPIBuildParticleExportList(nsubset, Partsubset, pfof, Len, sqrt(param[1]));
#endif
    //Now that have FoFDataGet (the exported particles) must search local volume using said particles
    //This is done by finding all particles in the search volume and then checking if those particles meet the FoF criterion
    //One must keep iterating till there are no new links.
    //Wonder if i don't need another loop and a final check
    Int_t links_across,links_across_total;
    do {
        links_across=MPILinkAcross(nsubset, tree, Partsubset, pfof, Len, Head, Next, param[1], fofcmp, param);
        MPIUpdateExportList(nsubset, Partsubset,pfof,Len);
        MPI_Allreduce(&links_across, &links_across_total, 1, MPI_Int_t, MPI_SUM, MPI_COMM_WORLD);
    }while(links_across_total>0);

    //reorder local particle array and delete memory associated with Head arrays, only need to keep Particles, pfof and some id and idexing information
    delete tree;
    delete[] Head;
    delete[] Next;
    delete[] Len;
    delete[] numingroup;
    delete[] FoFDataIn;
    delete[] FoFDataGet;
    delete[] PartDataIn;
    delete[] PartDataGet;

    //Now redistribute groups so that they are local to a processor (also orders the group ids according to size
    if (opt.iSingleHalo) opt.MinSize=MinNumOld;//reset minimum size
    Int_t newnbodies=MPIGroupExchange(nsubset,Partsubset,pfof);
    ///\todo need to clean up this mpi section for single halo
/*
#ifndef MPIREDUCEMEM
    //once groups are local, can free up memory
    delete[] Partsubset;
    delete[] pfof;
    //delete[] mpi_idlist;//since particles have now moved, must generate new list
    //mpi_idlist=new Int_t[newnbodies];
    pfof=new Int_t[newnbodies];
    //store new particle data in mpi_Part1 as external variable ensures memory allocated is not deallocated when function returns
    mpi_Part1=new Particle[newnbodies];
    Partsubset=mpi_Part1;
#endif
*/
    ///\todo Before final compilation of data, should have unbind here but must adjust unbind so it
    ///does not call reordergroupids in it though it might be okay.
    //And compile the information and remove groups smaller than minsize
    numgroups=MPICompileGroups(newnbodies,Partsubset,pfof,opt.MinSize);
    MPI_Barrier(MPI_COMM_WORLD);
    cout<<"MPI thread "<<ThisTask<<" has found "<<numgroups<<endl;
    //free up memory now that only need to store pfof and global ids
    if (ThisTask==0) {
        int totalgroups=0;
        for (int j=0;j<NProcs;j++) totalgroups+=mpi_ngroups[j];
        cout<<"Total number of groups found is "<<totalgroups<<endl;
    }
    //free up memory now that only need to store pfof and global ids
    delete[] FoFGroupDataLocal;
    delete[] FoFGroupDataExport;
    Nlocal=newnbodies;
    }
    else{
#endif
        delete tree;
#ifdef USEMPI
    }
#endif
    if (opt.iverbose>=2) cout<<"Done search for substructure in this subset"<<endl;
    return pfof;
}

//search for unassigned background particles if cores have been found.
void HaloCoreGrowth(Options &opt, const Int_t nsubset, Particle *&Partsubset, Int_t *&pfof, Int_t *&pfofbg, Int_t &numgroupsbg, Double_t param[], vector<Double_t> &dispfac,
    int numactiveloops, vector<int> &corelevel,
    int nthreads){
    //for simplicity make a new particle array storing core particles
    Int_t nincore=0,nbucket=opt.Bsize,pid, pidcore;
    Particle *Pcore,*Pval;
    KDTree *tcore;
    Coordinate x1;
    Double_t D2,dval,mval;
    Double_t *mcore=new Double_t[numgroupsbg+1];
    Int_t *ncore=new Int_t[numgroupsbg+1];
    Int_t newnumgroupsbg=0,*newcore=new Int_t[numgroupsbg+1];
    int nsearch=opt.Nvel;
    int mincoresize;
    int tid,i;
    Int_t **nnID;
    Double_t **dist2;
    PriorityQueue *pq;
    Int_t nactivepart=nsubset;

    for (i=0;i<=numgroupsbg;i++)ncore[i]=mcore[i]=0;
    //determine the weights for the cores dispersions factors
    for (i=0;i<nsubset;i++) {
        if (pfofbg[i]>0) {
            nincore++;
            mcore[pfofbg[i]]++;
            ncore[pfofbg[i]]++;
        }
    }

    if (opt.iverbose>=2) {
        cout<<"Mass ratios of cores are "<<endl;
        for (i=1;i<=numgroupsbg;i++)cout<<i<<" "<<mcore[i]<<" "<<mcore[i]/mcore[1]<<" "<<ncore[i]<<" "<<dispfac[i]<<endl;
    }
    //if number of particles in core less than number in subset then start assigning particles
    if (nincore<nsubset) {
        Int_t *noffset=new Int_t[numgroupsbg+1];
        //if running fully adaptive core linking, then need to calculate phase-space dispersions for each core
        //about their centres and use this to determine distances
        if (opt.iPhaseCoreGrowth) {
            if (opt.iverbose>=2) cout<<"Searching untagged particles to assign to cores using full phase-space metrics"<<endl;
            vector<GMatrix> dist(nthreads,GMatrix(6,1));
            vector<GMatrix> cmphase(numgroupsbg+1,GMatrix(6,1));
            vector<GMatrix> invdisp(numgroupsbg+1,GMatrix(6,6));
            GMatrix coredist(6,1);
            Int_t nactive=0;

            //store particles
            Pcore=new Particle[nincore];
            nincore=0;
            for (i=0;i<nsubset;i++) if (pfofbg[Partsubset[i].GetID()]>0) {
                Pcore[nincore]=Partsubset[i];
                Pcore[nincore].SetType(pfofbg[Partsubset[i].GetID()]);
                nincore++;
            }
            qsort(Pcore,nincore,sizeof(Particle),TypeCompare);
            noffset[0]=noffset[1]=0;
            for (i=2;i<=numgroupsbg;i++) noffset[i]=noffset[i-1]+ncore[i-1];
            //now get centre of masses and dispersions
            for (i=1;i<=numgroupsbg;i++) {
                cmphase[i]=CalcPhaseCM(ncore[i], &Pcore[noffset[i]]);
                for (int j=0;j<ncore[i];j++) {
                    for (int k=0;k<6;k++) Pcore[noffset[i]+j].SetPhase(k,Pcore[noffset[i]+j].GetPhase(k)-cmphase[i](k,0));
                }
                CalcPhaseSigmaTensor(ncore[i], &Pcore[noffset[i]], invdisp[i]);
                ///\todo must be issue with either phase-space tensor or number of particles assigned as
                ///it is possible to get haloes of size 0
                invdisp[i]=invdisp[i].Inverse();
            }
            delete[] Pcore;

            //once phase-space centers and dispersions are calculated, check to see
            //if distance is significant. Here idea is get distance in dispersion of
            //candidate core and this must be by ND*halocoredistsig, where ND is number of dimensions, ie. 6
            //if core is not significant set its mcore to 0
            for (i=2;i<=numgroupsbg;i++) {
                for (int k=0;k<6;k++) coredist(k,0)=(cmphase[i](k,0)-cmphase[1](k,0));
                D2=(coredist.Transpose()*invdisp[i]*coredist)(0,0);
                if (D2<opt.halocorephasedistsig*opt.halocorephasedistsig*6.0) mcore[i]=0;
                else nactive++;
            }

            //if there are no active cores then return nothing
            if (nactive==0) {
                delete[] mcore;
                delete[] ncore;
                numgroupsbg=0;
                return;
            }
            //if doing simple phase core growth then not recalculating dispersion as one goes up loops
            if (opt.iPhaseCoreGrowth==1) numactiveloops=0;
            //otherwise recalculating dispersions at every level
            else if (opt.iPhaseCoreGrowth>=2) for (i=1;i<=numgroupsbg;i++) dispfac[i]=1.0;

            for (Int_t iloop=numactiveloops;iloop>=0;iloop--) {
#ifdef USEOPENMP
            //if particle number large enough to warrant parallel search
            if (nactivepart>ompperiodnum) {
            int nreduce=0;
#pragma omp parallel default(shared) \
private(i,tid,Pval,D2,dval,mval,pid)
{
#pragma omp for reduction(+:nreduce)
            for (i=0;i<nsubset;i++)
            {
                tid=omp_get_thread_num();
                Pval=&Partsubset[i];
                if (Pval->GetType()<iloop) continue;
                pid=Pval->GetID();
                if (pfofbg[pid]==0 && pfof[pid]==0) {
                    mval=mcore[1];
                    for (int k=0;k<6;k++) dist[tid](k,0)=Pval->GetPhase(k)-cmphase[1](k,0);
                    dval=(dist[tid].Transpose()*invdisp[1]*dist[tid])(0,0);
                    pfofbg[pid]=1;
                    for (int j=2;j<=numgroupsbg;j++) if (mcore[j]>0 && corelevel[j]>=iloop){
                        for (int k=0;k<6;k++) dist[tid](k,0)=Pval->GetPhase(k)-cmphase[j](k,0);
                        D2=(dist[tid].Transpose()*invdisp[j]*dist[tid])(0,0);
                        if (dval*dispfac[pfofbg[pid]]>D2*dispfac[j]) {dval=D2;mval=mcore[j];pfofbg[pid]=j;}
                    }
                    //if particle assigned to a core remove from search
                    Pval->SetType(-1);
                    nreduce++;
                }
            }
}
            nactivepart-=nreduce;
            }
            else {
#endif
            for (i=0;i<nsubset;i++)
            {
                tid=0;
                Pval=&Partsubset[i];
                if (Pval->GetType()<iloop) continue;
                pid=Pval->GetID();
                if (pfofbg[pid]==0 && pfof[pid]==0) {
                    mval=mcore[1];
                    for (int k=0;k<6;k++) dist[tid](k,0)=Pval->GetPhase(k)-cmphase[1](k,0);
                    dval=(dist[tid].Transpose()*invdisp[1]*dist[tid])(0,0);
                    pfofbg[pid]=1;
                    for (int j=2;j<=numgroupsbg;j++) if (mcore[j]>0 && corelevel[j]>=iloop){
                        for (int k=0;k<6;k++) dist[tid](k,0)=Pval->GetPhase(k)-cmphase[j](k,0);
                        D2=(dist[tid].Transpose()*invdisp[j]*dist[tid])(0,0);
                        if (dval*dispfac[pfofbg[pid]]>D2*dispfac[j]) {dval=D2;mval=mcore[j];pfofbg[pid]=j;}
                    }
                    Pval->SetType(-1);
                    nactivepart--;
                }
            }
#ifdef USEOPENMP
            }
#endif
            //otherwise, recalculate dispersions
            if (opt.iPhaseCoreGrowth>=2) {
                nincore=0;
                for (i=0;i<nsubset;i++)
                    if (pfofbg[i]>0)
                        nincore++;
                Pcore=new Particle[nincore];
                for (i=1;i<=numgroupsbg;i++) ncore[i]=0;
                nincore=0;
                for (i=0;i<nsubset;i++) if (pfofbg[Partsubset[i].GetID()]>0) {
                    Pcore[nincore]=Partsubset[i];
                    Pcore[nincore].SetType(pfofbg[Partsubset[i].GetID()]);
                    nincore++;
                    ncore[pfofbg[Partsubset[i].GetID()]]++;
                }
                qsort(Pcore,nincore,sizeof(Particle),TypeCompare);
                noffset[0]=noffset[1]=0;
                for (i=2;i<=numgroupsbg;i++) noffset[i]=noffset[i-1]+ncore[i-1];
                //now get centre of masses and dispersions
                for (i=1;i<=numgroupsbg;i++) if (corelevel[i]>=iloop) {
                    cmphase[i]=CalcPhaseCM(ncore[i], &Pcore[noffset[i]]);
                    for (int j=0;j<ncore[i];j++) {
                        for (int k=0;k<6;k++) Pcore[noffset[i]+j].SetPhase(k,Pcore[noffset[i]+j].GetPhase(k)-cmphase[i](k,0));
                    }
                    CalcPhaseSigmaTensor(ncore[i], &Pcore[noffset[i]], invdisp[i]);
                    invdisp[i]=invdisp[i].Inverse();
                }
                delete[] Pcore;
            }
        }//
        }//end of phase core growth
        //otherwise, use simplier calculation: find nearest particles belonging to cores, calculate distances to these particles and assign untagged
        //particle to the same group as the closest core particle
        else {
            if (opt.iverbose>=2) cout<<"Searching untagged particles to assign to cores using simple distance measure to nearest core particle"<<endl;
            if (nsearch>nincore) nsearch=nincore;
            if (nbucket>=nincore/8) nbucket=max(1,(int)nincore/8);
            Pcore=new Particle[nincore];
            nincore=0;
            for (i=0;i<nsubset;i++) if (pfofbg[Partsubset[i].GetID()]>0) {
                Pcore[nincore]=Partsubset[i];
                Pcore[nincore].SetType(pfofbg[Partsubset[i].GetID()]);
                nincore++;
            }
            tcore=new KDTree(Pcore,nincore,opt.Bsize,tcore->TPHYS);
            nnID=new Int_t*[nthreads];
            dist2=new Double_t*[nthreads];
            for (i=0;i<nthreads;i++) {
                nnID[i]=new Int_t[nsearch];
                dist2[i]=new Double_t[nsearch];
            }
            for (i=1;i<=numgroupsbg;i++) ncore[i]=0;
#ifdef USEOPENMP
            //if particle number large enough to warrant parallel search
            if (nsubset>ompperiodnum) {
#pragma omp parallel default(shared) \
private(i,tid,Pval,x1,D2,dval,mval,pid,pidcore)
{
#pragma omp for
            //for each particle in the subset if not assigned to any group (core or substructure) then assign particle
            //this is done using either a simple distance/sigmax+velocity distance/sigmav calculation to the nearest core particles
            //or if a more complex routine is required then ...
            for (i=0;i<nsubset;i++)
            {
                tid=omp_get_thread_num();
                Pval=&Partsubset[i];
                pid=Pval->GetID();
                if (pfofbg[pid]==0 && pfof[pid]==0) {
                    x1=Coordinate(Pval->GetPosition());
                    tcore->FindNearestPos(x1, nnID[tid], dist2[tid],nsearch);
                    dval=0;
                    pidcore=nnID[tid][0];
                    //calculat distance from current particle to core particle
                    for (int k=0;k<3;k++) {
                        dval+=(Pval->GetPosition(k)-Pcore[pidcore].GetPosition(k))*(Pval->GetPosition(k)-Pcore[pidcore].GetPosition(k))/param[6]+(Pval->GetVelocity(k)-Pcore[pidcore].GetVelocity(k))*(Pval->GetVelocity(k)-Pcore[pidcore].GetVelocity(k))/param[7];
                    }
                    //get the core particle mass ratio
                    mval=mcore[Pcore[pidcore].GetType()];
                    pfofbg[pid]=Pcore[pidcore].GetType();
                    //now initialized to first core particle, examine the rest to see if one is closer
                    for (int j=1;j<nsearch;j++) {
                        D2=0;
                        pidcore=nnID[tid][j];
                        for (int k=0;k<3;k++) {
                            D2+=(Pval->GetPosition(k)-Pcore[pidcore].GetPosition(k))*(Pval->GetPosition(k)-Pcore[pidcore].GetPosition(k))/param[6]+(Pval->GetVelocity(k)-Pcore[pidcore].GetVelocity(k))*(Pval->GetVelocity(k)-Pcore[pidcore].GetVelocity(k))/param[7];
                        }
                        //if distance * mass weight is smaller than current distance, reassign particle
                        if (dval>D2*mval/mcore[Pcore[pidcore].GetType()]) {dval=D2;mval=mcore[Pcore[pidcore].GetType()];pfofbg[pid]=Pcore[pidcore].GetType();}
                    }
                }
            }
}
            }
            else {
#endif
            for (i=0;i<nsubset;i++)
            {
                tid=0;
                Pval=&Partsubset[i];
                pid=Pval->GetID();
                if (pfofbg[pid]==0 && pfof[pid]==0) {
                    x1=Coordinate(Pval->GetPosition());
                    tcore->FindNearestPos(x1, nnID[tid], dist2[tid],nsearch);
                    dval=0;
                    pidcore=nnID[tid][0];
                    for (int k=0;k<3;k++) {
                        dval+=(Pval->GetPosition(k)-Pcore[pidcore].GetPosition(k))*(Pval->GetPosition(k)-Pcore[pidcore].GetPosition(k))/param[6]+(Pval->GetVelocity(k)-Pcore[pidcore].GetVelocity(k))*(Pval->GetVelocity(k)-Pcore[pidcore].GetVelocity(k))/param[7];
                    }
                    mval=mcore[Pcore[pidcore].GetType()];
                    pfofbg[pid]=Pcore[pidcore].GetType();
                    for (int j=1;j<nsearch;j++) {
                        D2=0;
                        pidcore=nnID[tid][j];
                        for (int k=0;k<3;k++) {
                            D2+=(Pval->GetPosition(k)-Pcore[pidcore].GetPosition(k))*(Pval->GetPosition(k)-Pcore[pidcore].GetPosition(k))/param[6]+(Pval->GetVelocity(k)-Pcore[pidcore].GetVelocity(k))*(Pval->GetVelocity(k)-Pcore[pidcore].GetVelocity(k))/param[7];
                        }
                        if (dval>D2*mval/mcore[Pcore[pidcore].GetType()]) {dval=D2;mval=mcore[Pcore[pidcore].GetType()];pfofbg[pid]=Pcore[pidcore].GetType();}
                    }
                }
            }
#ifdef USEOPENMP
            }
#endif
            //clean up memory
            delete tcore;
            delete[] Pcore;
            for (i=0;i<nthreads;i++) {
                delete [] nnID[i];
                delete [] dist2[i];
            }
            delete[] nnID;
            delete[] dist2;
        }
        //now that particles assigned to cores, remove if core too small
        if (opt.partsearchtype!=PSTSTAR&&opt.foftype!=FOF6DCORE) mincoresize=max((Int_t)(nsubset*opt.halocorenfac),(Int_t)opt.MinSize);//max((Int_t)(nsubset*MAXCELLFRACTION/2.0),(Int_t)opt.MinSize);
        else mincoresize=opt.MinSize;
        for (i=1;i<=numgroupsbg;i++) ncore[i]=0;
        for (i=0;i<nsubset;i++) {
            if (pfofbg[i]>0 && mcore[pfofbg[i]]==0) pfofbg[i]=0;
            ncore[pfofbg[i]]++;
        }
        pq=new PriorityQueue(numgroupsbg);
        newnumgroupsbg=0;
        for (i=1;i<=numgroupsbg;i++){
            if (ncore[i]<mincoresize) ncore[i]=0;
            else newnumgroupsbg++;
            pq->Push(i,(Double_t)ncore[i]);
        }
        //if no large enough cores are left, then halt
        if (newnumgroupsbg<=1) {
            numgroupsbg=0;
            delete pq;
            delete[] mcore;
            delete[] ncore;
            delete[] newcore;
            return;
        }
        newnumgroupsbg=0;
        for (i=1;i<=numgroupsbg;i++){
            if (pq->TopPriority()>=mincoresize) newcore[pq->TopQueue()]=++newnumgroupsbg;
            pq->Pop();
        }
        for (i=0;i<nsubset;i++) {
            pid=Partsubset[i].GetID();
            if (pfofbg[pid]>0 && ncore[pfofbg[pid]]>0) pfofbg[pid]=newcore[pfofbg[pid]];
            else if  (pfofbg[pid]>0 && ncore[pfofbg[pid]]==0) pfofbg[pid]=0;
        }
        numgroupsbg=newnumgroupsbg;
        delete pq;
        delete[] mcore;
        delete[] ncore;
        delete[] newcore;
    }
}


//Merge any groups that overlap in phase-space
void MergeSubstructuresCoresPhase(Options &opt, const Int_t nsubset, Particle *&Partsubset, Int_t *&pfof, Int_t &numsubs, Int_t &numcores)
{
    //get the phase centres of objects and see if they overlap
    Int_t pfofval, imerge, numlargesubs=0, newnumcores=numcores;
    Double_t disp, dist2, mindist2, fdist2=pow(opt.coresubmergemindist,2.0);
    Coordinate pos;
    vector<Int_t> numingroup, noffset, taggedsubs;
    vector<Particle> subs, cores;
    KDTree *tree;
    //vector<GMatrix> phasetensorsubs(numsubs,GMatrix(6,6)), phasetensorcores(numcores,GMatrix(6,6));
    vector<Double_t> sigXsubs(numsubs), sigVsubs(numsubs), sigXcores(numcores), sigVcores(numcores);
    struct indexfof {
        Int_t fofval;
        Int_t index;
    };
    vector<indexfof> indexing;
    subs.resize(numsubs);
    cores.resize(numcores);
    numingroup.resize(numsubs+numcores+1);
    noffset.resize(numsubs+numcores+1);
    indexing.resize(nsubset);
    for (auto &x:sigXsubs) x=0;
    for (auto &x:sigVsubs) x=0;
    for (auto &x:sigXcores) x=0;
    for (auto &x:sigVcores) x=0;
    //get center of mass in phase-space
    for (auto i=0;i<nsubset;i++) {
        pfofval = pfof[Partsubset[i].GetID()];
        indexing[i].fofval = pfofval;
        indexing[i].index = i;
        numingroup[pfofval]++;
        if (pfofval==0) continue;
        if (pfofval<=numsubs) {
            pfofval-=1;
            subs[pfofval].SetMass(subs[pfofval].GetMass()+Partsubset[i].GetMass());
            for (auto k=0;k<6;k++) subs[pfofval].SetPhase(k,subs[pfofval].GetPhase(k)+Partsubset[i].GetPhase(k)*Partsubset[i].GetMass());
        }
        else {
            pfofval-=numsubs+1;
            cores[pfofval].SetMass(cores[pfofval].GetMass()+Partsubset[i].GetMass());
            for (auto k=0;k<6;k++) cores[pfofval].SetPhase(k,cores[pfofval].GetPhase(k)+Partsubset[i].GetPhase(k)*Partsubset[i].GetMass());
        }
    }
    noffset[0]=0; for (auto i=1;i<=numsubs+numcores;i++) noffset[i]=numingroup[i-1]+noffset[i-1];
    pfofval=1;
    for (auto &x:subs) {
        x.SetPID(pfofval);
        x.SetID(pfofval);
        pfofval++;
        for (auto k=0;k<6;k++) x.SetPhase(k,x.GetPhase(k)/x.GetMass());
    }
    for (auto &x:cores) {
        x.SetPID(pfofval);
        x.SetID(pfofval);
        pfofval++;
        for (auto k=0;k<6;k++) x.SetPhase(k,x.GetPhase(k)/x.GetMass());
    }
    //sort indices by fof value
    sort(indexing.begin(), indexing.end(), [](indexfof &a, indexfof &b){
    return a.fofval < b.fofval;
    });
    //get the dispersions
    for (auto i=0;i<nsubset;i++) {
        pfofval = pfof[Partsubset[i].GetID()];
        if (pfofval==0) continue;
        if (pfofval<=numsubs) {
            pfofval-=1;
            disp=0; for (auto k=0;k<3;k++) disp+=pow(Partsubset[i].GetPosition(k)-subs[pfofval].GetPosition(k),2.0);
            sigXsubs[pfofval]+=disp*Partsubset[i].GetMass();
            disp=0; for (auto k=0;k<3;k++) disp+=pow(Partsubset[i].GetVelocity(k)-subs[pfofval].GetVelocity(k),2.0);
            sigVsubs[pfofval]+=disp*Partsubset[i].GetMass();
        }
        else {
            pfofval-=numsubs+1;
            disp=0; for (auto k=0;k<3;k++) disp+=pow(Partsubset[i].GetPosition(k)-cores[pfofval].GetPosition(k),2.0);
            sigXcores[pfofval]+=disp*Partsubset[i].GetMass();
            disp=0; for (auto k=0;k<3;k++) disp+=pow(Partsubset[i].GetVelocity(k)-cores[pfofval].GetVelocity(k),2.0);
            sigVcores[pfofval]+=disp*Partsubset[i].GetMass();
        }
    }
    for (auto i=0;i<numsubs;i++) {
        sigXsubs[i]*=1.0/subs[i].GetMass();
        sigVsubs[i]*=1.0/subs[i].GetMass();
    }
    for (auto i=0;i<numcores;i++) {
        sigXcores[i]*=1.0/cores[i].GetMass();
        sigVcores[i]*=1.0/cores[i].GetMass();
    }
    //now built tree on substructures
    tree = new KDTree(subs.data(),numsubs,1,tree->TPHYS,tree->KEPAN,100,0,0,0);
    //tree = new KDTree(subs.data(),numlargesubs,1,tree->TPHYS,tree->KEPAN,100,0,0,0);
    //check all cores to see if they overlap significantly with substructures
    newnumcores=0;
    for (auto i=0;i<numcores;i++) {
        for (auto k=0;k<3;k++) pos[k]=cores[i].GetPosition(k);
        taggedsubs = tree->SearchBallPosTagged(pos, sigXcores[i]*fdist2);
        if (taggedsubs.size()==0) {
            newnumcores++;
            cores[i].SetID(newnumcores+numsubs);
            continue;
        }
        //if objects are within search window of core, get min phase distance
        imerge=-1;
        mindist2=MAXVALUE;
        for (auto j=0;j<taggedsubs.size();j++) {
            disp = 0; for (auto k=0;k<3;k++) disp+=pow(subs[taggedsubs[j]].GetPosition(k)-cores[i].GetPosition(k),2.0);
            dist2 = disp/sigXcores[i];
            disp = 0; for (auto k=0;k<3;k++) disp+=pow(subs[taggedsubs[j]].GetVelocity(k)-cores[i].GetVelocity(k),2.0);
            dist2 += disp/sigVcores[i];
            if (dist2<fdist2 && dist2<mindist2){
                imerge=subs[taggedsubs[j]].GetPID();
                mindist2=dist2;
            }
        }
        //merging core with sub if one is found
        if (imerge!=-1) {
            pfofval=i+numsubs+1;
            for (auto j=noffset[pfofval];j<noffset[pfofval]+numingroup[pfofval];j++) {
                pfof[Partsubset[indexing[j].index].GetID()]=imerge;
            }
            cores[i].SetID(imerge);
            cores[i].SetPID(0);
        }
        else {
            newnumcores++;
            cores[i].SetID(newnumcores+numsubs);
        }
    }
    delete tree;
    //if object has been mergerged then must update the pfof values of cores
    if (newnumcores!=numcores && newnumcores>0) {
        for (auto i=0;i<numcores;i++) {
            pfofval=i+numsubs+1;
            if (cores[i].GetPID()==0) continue;
            if (cores[i].GetPID()==cores[i].GetID()) continue;
            for (auto j=noffset[pfofval];j<noffset[pfofval]+numingroup[pfofval];j++) {
                pfof[Partsubset[indexing[j].index].GetID()]=cores[i].GetID();
            }
        }
    }
    numcores=newnumcores;
}
///Testing a merge of all substructures
void MergeSubstructuresPhase(Options &opt, const Int_t nsubset, Particle *&Partsubset, Int_t *&pfof, Int_t &numgroups, Int_t &numsubs, Int_t &numcores)
{
    //get the phase centres of objects and see if they overlap
    Int_t pfofval, newpfofval, imerge, newnumgroups, newnumcores, nummerged=0, index1, index2;
    Double_t disp, dist2, dist2sub1,dist2sub2, mindist2, fdist2=pow(opt.coresubmergemindist,2.0);
    Coordinate pos;
    vector<Int_t> numingroup, noffset, taggedsubs;
    struct mergeinfo {
        Int_t pfofval;
        Int_t numingroup;
        int type;
        int nummerged;
        bool ismerged;
        vector<Int_t> mergedlist;
        mergeinfo(){
            nummerged=0;
            ismerged=false;
        };
        //Int_t numingroup;
    };
    vector<Particle> subs;
    vector<mergeinfo> minfo;
    KDTree *tree;
    //vector<GMatrix> phasetensorsubs(numsubs,GMatrix(6,6)), phasetensorcores(numcores,GMatrix(6,6));
    vector<Double_t> sigXsubs(numgroups), sigVsubs(numgroups), sigXcores(numgroups), sigVcores(numgroups);
    struct indexfof {
        Int_t fofval;
        Int_t index;
    };
    vector<indexfof> indexing;
//cout<<ThisTask<<" entering merging "<<numgroups<<" "<<numsubs<<" "<<numcores<<endl;
    subs.resize(numgroups);
    minfo.resize(numgroups);
    numingroup.resize(numgroups+1);
    noffset.resize(numgroups+1);
    indexing.resize(nsubset);
    for (auto &x:sigXsubs) x=0;
    for (auto &x:sigVsubs) x=0;
    //get center of mass in phase-space
    for (auto i=0;i<nsubset;i++) {
        pfofval = pfof[Partsubset[i].GetID()];
        indexing[i].fofval = pfofval;
        indexing[i].index = i;
        numingroup[pfofval]++;
        if (pfofval==0) continue;
        pfofval-=1;
        subs[pfofval].SetMass(subs[pfofval].GetMass()+Partsubset[i].GetMass());
        for (auto k=0;k<6;k++) subs[pfofval].SetPhase(k,subs[pfofval].GetPhase(k)+Partsubset[i].GetPhase(k)*Partsubset[i].GetMass());
    }
    noffset[0]=0; for (auto i=1;i<=numgroups;i++) noffset[i]=numingroup[i-1]+noffset[i-1];
    pfofval=1;
    for (auto &x:subs) {
        x.SetPID(pfofval);
        //x.SetID(pfofval);
        x.SetType((pfofval>numsubs));
        minfo[pfofval-1].pfofval=pfofval;
        minfo[pfofval-1].type = x.GetType();
        minfo[pfofval-1].numingroup=numingroup[pfofval];
        pfofval++;
        for (auto k=0;k<6;k++) x.SetPhase(k,x.GetPhase(k)/x.GetMass());
    }
    //sort indices by fof value
    sort(indexing.begin(), indexing.end(), [](indexfof &a, indexfof &b){
    return a.fofval < b.fofval;
    });
    //get the dispersions
    for (auto i=0;i<nsubset;i++) {
        pfofval = pfof[Partsubset[i].GetID()];
        if (pfofval==0) continue;
        pfofval-=1;
        disp=0; for (auto k=0;k<3;k++) disp+=pow(Partsubset[i].GetPosition(k)-subs[pfofval].GetPosition(k),2.0);
        sigXsubs[pfofval]+=disp*Partsubset[i].GetMass();
        disp=0; for (auto k=0;k<3;k++) disp+=pow(Partsubset[i].GetVelocity(k)-subs[pfofval].GetVelocity(k),2.0);
        sigVsubs[pfofval]+=disp*Partsubset[i].GetMass();
    }
    for (auto i=0;i<numgroups;i++) {
        sigXsubs[i]*=1.0/subs[i].GetMass();
        sigVsubs[i]*=1.0/subs[i].GetMass();
    }
    //now built tree on substructures
    tree = new KDTree(subs.data(),numgroups,1,tree->TPHYS,tree->KEPAN,100,0,0,0);
    //tree = new KDTree(subs.data(),numlargesubs,1,tree->TPHYS,tree->KEPAN,100,0,0,0);
    //check all cores to see if they overlap significantly with substructures
    for (auto i=0;i<numgroups;i++) {
        //ignore if already merged
        if (subs[i].GetPID()==0) continue;
        //don't merge substructures into cores
        if (subs[i].GetType()==1) continue;
        index1=subs[i].GetPID()-1;
        index1=subs[i].GetID();
        taggedsubs = tree->SearchBallPosTagged(i, sigXsubs[index1]*fdist2);
//cout<<i<<" sub tagged "<<subs[i].GetPID()<<" "<<numingroup[subs[i].GetPID()]<<" "<<taggedsubs.size()<<endl;
        if (taggedsubs.size()<=1) continue;
        //if objects are within search window of core, get min phase distance
        imerge=-1;
        mindist2=MAXVALUE;
        for (auto j=0;j<taggedsubs.size();j++) {
            if (i==taggedsubs[j]) continue;
            if (subs[taggedsubs[j]].GetPID()==0) continue;
            //index2=subs[taggedsubs[j]].GetPID()-1;
            index2=subs[taggedsubs[j]].GetID();
//double disp2;
            disp = 0; for (auto k=0;k<3;k++) disp+=pow(subs[taggedsubs[j]].GetPosition(k)-subs[i].GetPosition(k),2.0);
//disp2=disp;
            dist2sub1 = disp/sigXsubs[index1];
            dist2sub2 = disp/sigXsubs[index2];
            disp = 0; for (auto k=0;k<3;k++) disp+=pow(subs[taggedsubs[j]].GetVelocity(k)-subs[i].GetVelocity(k),2.0);
            dist2sub1 += disp/sigVsubs[index1];
            dist2sub2 = disp/sigVsubs[index2];
            dist2 = 0.5*(dist2sub1+dist2sub2);
//cout<<" dist "<<subs[i].GetPID()<<" "<<subs[taggedsubs[j]].GetPID()<<" "<<numingroup[subs[i].GetPID()]<<" "<<numingroup[subs[taggedsubs[j]].GetPID()]<<" "<<sqrt(dist2)<<" || "<<sqrt(disp2)<<" "<<sqrt(disp)<<" "<<sqrt(sigXsubs[index1])<<" "<<sqrt(sigVsubs[index1])<<" "<<sqrt(sigXsubs[index2])<<" "<<sqrt(sigVsubs[index2])<<endl;
            if (dist2sub1<fdist2 && dist2sub2<fdist2 && dist2<mindist2){
                imerge=taggedsubs[j];
                mindist2=dist2;
            }
        }
        //merging core with sub if one is found
        if (imerge!=-1) {
            nummerged++;
            pfofval=subs[imerge].GetPID();
            index2=subs[imerge].GetID();
            subs[imerge].SetPID(0);
            minfo[index2].ismerged=true;
            minfo[index1].nummerged = minfo[index2].nummerged+1;
            minfo[index1].numingroup+=minfo[index2].numingroup;
            minfo[index1].mergedlist.push_back(pfofval);
            for (auto j=0;j<minfo[index2].nummerged;j++) minfo[index1].mergedlist.push_back(minfo[index2].mergedlist[j]);
//cout<<"merging "<<subs[i].GetPID()<<" "<<pfofval<<" "<<minfo[index1].nummerged<<endl;
        }
    }
    delete tree;
    //if nothing has changed, do nothing
    if (nummerged==0) return;
//cout<<" merging subs has given  "<<nummerged<<" mergers "<<numgroups<<" "<<numcores<<endl;
    sort(minfo.begin(), minfo.end(), [](mergeinfo &a, mergeinfo &b){
        if (a.type<b.type) return true;
        else if (a.type==b.type) return (a.numingroup > b.numingroup);
        else return false;
    });
    newnumgroups=0;
    newnumcores=0;
    for (auto i=0;i<numgroups;i++) {
//cout<<i<<" "<<minfo[i].numingroup<<" "<<minfo[i].type<<" and before "<<numingroup[minfo[i].pfofval]<<endl;
        //if object has mergered do nothing
        if (minfo[i].ismerged==true) continue;
        newnumgroups++;
        if (minfo[i].type == 1) newnumcores++;
        //if object is still in same order and has not mergered with anything, do nothing
        if (minfo[i].pfofval == newnumgroups && minfo[i].nummerged==0) continue;
//cout<<"altering  "<<minfo[i].pfofval<<" to "<<newnumgroups<<endl;
        pfofval=minfo[i].pfofval;
        for (auto j=noffset[pfofval];j<noffset[pfofval]+numingroup[pfofval];j++) {
            pfof[Partsubset[indexing[j].index].GetID()]=newnumgroups;
        }
        if (minfo[i].nummerged==0) continue;
//cout<<"altering merged objects"<<minfo[i].nummerged<<" to "<<newnumgroups<<endl;
        for (auto &mergedgroup:minfo[i].mergedlist) {
//cout<<" halo that is merging "<<mergedgroup<<" "<<numingroup[mergedgroup]<<endl;
            for (auto j=noffset[mergedgroup];j<noffset[mergedgroup]+numingroup[mergedgroup];j++) {
                pfof[Partsubset[indexing[j].index].GetID()]=newnumgroups;
            }
        }
    }
//    for (auto i=0;i<nsubset;i++) if (pfof[Partsubset[i].GetID()]>newnumgroups) cout<<"WTF incorrect pfof value "<<i<<" "<<pfof[Partsubset[i].GetID()]<<" "<<newnumgroups<<endl;
//cout<<"new is "<<newnumgroups<<" "<<newnumcores<<endl;
    //exit(9);
    numcores=newnumcores;
    numgroups=newnumgroups;
    numsubs=numgroups-numcores;
//    cout<<ThisTask<<" after merging "<<numgroups<<" "<<numsubs<<" "<<numcores<<endl;
}

// ENCAPSULATED: ENCAPSULATION-01
void AdjustSubPartToPhaseCM(Int_t num, Particle *subPart, GMatrix &cmphase)
{
    int nthreads = 1;
#ifdef USEOPENMP
    cout << "OPTIMISATION-01: BEGIN PRINT" << endl;
    cout << "OPTIMISATION-01: num[" << num << "]" << endl;
    cout << "OPTIMISATION-01: omp_get_max_threads[" << omp_get_max_threads() << "]" << endl;
    cout << "OPTIMISATION-01: ompsearchnum[" << ompsearchnum << "]" << endl;
    nthreads = max(1, (int)(num/(float)ompsearchnum));
    cout << "OPTIMISATION-01: nthreads.1[" << nthreads << "]" << endl;
    nthreads = min(nthreads,omp_get_max_threads());
    cout << "OPTIMISATION-01: nthreads.2[" << nthreads << "]" << endl;
    cout << "OPTIMISATION-01: END PRINT" << endl;
#pragma omp parallel for \
default(shared)  \
num_threads(nthreads)
#endif
    for (auto j=0;j<num;j++)
    {
        for (int k=0;k<6;k++) subPart[j].SetPhase(k,subPart[j].GetPhase(k)-cmphase(k,0));
    }
}

// ENCAPSULATED: ENCAPSULATION-02
void PreCalcSearchSubSet(Int_t *subnumingroup, Int_t i, Options &opt, Particle *subPart, 
    Int_t sublevel, KDTree *tree, Int_t ngrid, int ThisTask, GridCell *grid, Coordinate *gvel, Matrix *gveldisp)
{
    if (subnumingroup[i]>=MINSUBSIZE&&opt.foftype!=FOF6DCORE) {
        //now if object is large enough for phase-space decomposition and search, compare local field to bg field
        opt.Ncell=opt.Ncellfac*subnumingroup[i];
        //if ncell is such that uncertainty would be greater than 0.5% based on Poisson noise, increase ncell till above unless cell would contain >25%
        while (opt.Ncell<MINCELLSIZE && subnumingroup[i]/4.0>opt.Ncell) opt.Ncell*=2;
        tree=InitializeTreeGrid(opt,subnumingroup[i],subPart);
        ngrid=tree->GetNumLeafNodes();
        if (opt.iverbose) cout<<ThisTask<<" Substructure "<<i<< " at sublevel "<<sublevel<<" with "<<subnumingroup[i]
            <<" particles split into are "<<ngrid<<" grid cells, with each node containing ~"<<subnumingroup[i]/ngrid<<" particles"<<endl;
        grid=new GridCell[ngrid];
        FillTreeGrid(opt, subnumingroup[i], ngrid, tree, subPart, grid);
        gvel=GetCellVel(opt,subnumingroup[i],subPart,ngrid,grid);
        gveldisp=GetCellVelDisp(opt,subnumingroup[i],subPart,ngrid,grid,gvel);
        opt.HaloLocalSigmaV=0;for (int j=0;j<ngrid;j++) opt.HaloLocalSigmaV+=pow(gveldisp[j].Det(),1./3.);opt.HaloLocalSigmaV/=(double)ngrid;

        Matrix eigvec(0.),I(0.);
        Double_t sigma2x,sigma2y,sigma2z;
        CalcVelSigmaTensor(subnumingroup[i], subPart, sigma2x, sigma2y, sigma2z, eigvec, I);
        opt.HaloSigmaV=pow(sigma2x*sigma2y*sigma2z,1.0/3.0);
        if (opt.HaloSigmaV>opt.HaloVelDispScale) opt.HaloVelDispScale=opt.HaloSigmaV;
#ifdef HALOONLYDEN
            GetVelocityDensity(opt,subnumingroup[i],subPart);
#endif
        GetDenVRatio(opt,subnumingroup[i],subPart,ngrid,grid,gvel,gveldisp);
        GetOutliersValues(opt,subnumingroup[i],subPart,sublevel);
        opt.idenvflag++;//largest field halo used to deteremine statistics of ratio
    }
    //otherwise only need to calculate a velocity scale for merger separation
    else {
        Matrix eigvec(0.),I(0.);
        Double_t sigma2x,sigma2y,sigma2z;
        CalcVelSigmaTensor(subnumingroup[i], subPart, sigma2x, sigma2y, sigma2z, eigvec, I);
        opt.HaloLocalSigmaV=opt.HaloSigmaV=pow(sigma2x*sigma2y*sigma2z,1.0/3.0);
    }
}

// ENCAPSULATED: ENCAPSULATION-03
void SearchSubStruct(Int_t *numcores, Int_t i, Options &opt, Int_t ng, bool iunbindflag, Int_t *subnumingroup,
    Particle *subPart, Int_t *subngroup, Int_t *subpfof, Int_t **subsubnumingroup, Int_t ***subsubpglist,
    Int_t *&pfof, Int_t **subpglist, Int_t &ngroup, Int_t ngroupidoffset) {
    Int_t *coreflag;
                ng=subngroup[i];
                subsubnumingroup[i]=BuildNumInGroup(subnumingroup[i], subngroup[i], subpfof);
                subsubpglist[i]=BuildPGList(subnumingroup[i], subngroup[i], subsubnumingroup[i], subpfof);
    if (opt.uinfo.unbindflag&&subngroup[i]>0) {
    //if also keeping track of cores then must allocate coreflag
    if (numcores[i]>0 && opt.iHaloCoreSearch>=1) {
        coreflag=new Int_t[ng+1];
        for (int icore=1;icore<=ng;icore++) coreflag[icore]=1+(icore>ng-numcores[i]);
    }
    else {coreflag=NULL;}
    iunbindflag=CheckUnboundGroups(opt,subnumingroup[i],subPart,subngroup[i],subpfof,subsubnumingroup[i],subsubpglist[i],1, coreflag);
    if (iunbindflag) {
        for (int j=1;j<=ng;j++) delete[] subsubpglist[i][j];
        delete[] subsubnumingroup[i];
        delete[] subsubpglist[i];
        if (subngroup[i]>0) {
            subsubnumingroup[i]=BuildNumInGroup(subnumingroup[i], subngroup[i], subpfof);
            subsubpglist[i]=BuildPGList(subnumingroup[i], subngroup[i], subsubnumingroup[i], subpfof);
        }
        //if need to update number of cores,
        if (numcores[i]>0 && opt.iHaloCoreSearch>=1) {
            numcores[i]=0;
            for (int icore=1;icore<=subngroup[i];icore++)numcores[i]+=(coreflag[icore]==2);
            delete[] coreflag;
        }
    }}
    for (Int_t j=0;j<subnumingroup[i];j++) if (subpfof[j]>0) pfof[subpglist[i][j]]=ngroup+ngroupidoffset+subpfof[j];

    //now alter subsubpglist so that index pointed is global subset index as global subset is used to get the particles to be searched for subsubstructure
    for (Int_t j=1;j<=subngroup[i];j++) for (Int_t k=0;k<subsubnumingroup[i][j];k++) subsubpglist[i][j][k]=subpglist[i][subsubpglist[i][j][k]];
}
/*!
    Given a initial ordered candidate list of substructures, find all substructures that are large enough to be searched.
    These substructures are used as a mean background velocity field and a new outlier list is found and searched.
    The subsubstructures must also be significant and if unbinding, must also be bound. If a particle in a subsubstructure is
    removed from the subsubstructure, its fof id is set to that of the substructure it is contained in.

    The minimum size is given by \ref MINSUBSIZE. This should be for objects composed of ~1000 particles or so as realistically,
    estimating the mean field needs at least 8 cells and ideally want about 100 particles per cell. Smaller objects will have
    cells used in mean field estimates that span too large a region in the phase-space of the object. The poorly resolved halos
    also have pretty noisy internal structures though they may have substructures composed of ~<10% of the particles (more
    realistially is ~<1%. So for 1000 that gives a maximum size of 100 particles.

    NOTE: if the code is altered and generalized to outliers in say the entropy distribution when searching for gas shocks,
    it might be possible to lower the cuts imposed.

    \todo ADACS optimisation request. Here the function could be altered to employ better parallelisation. Specifically, the loop over
    substructures at a given level could be parallelized (see for loop commented with ENCAPSULATE-01). Currently, at a given level in the substructure hierarchy
    each object is searched sequentially but this does not need to be the case. It would require restructureing the loop and some of calls within
    the loop so that the available pool of threads over which to run in parallel for the callled subroutines is adaptive. (Or it might be
    simply more useful to not have the functions called within this loop parallelised. This loop invokes a few routines that have OpenMP
    parallelisation: InitializeTreeGrid, GetCellVel, GetCellVelDisp, CalcVelSigmaTensor, etc. 
*/
void SearchSubSub(Options &opt, const Int_t nsubset, vector<Particle> &Partsubset, Int_t *&pfof, Int_t &ngroup, Int_t &nhalos, PropData *pdata)
{
    //now build a sublist of groups to search for substructure
    Int_t nsubsearch, oldnsubsearch,sublevel,maxsublevel,ngroupidoffset,ngroupidoffsetold,ngrid;
    bool iflag,iunbindflag;
    Particle *subPart;
    Int_t firstgroup,firstgroupoffset;
    Int_t ng,*numingroup,**pglist;
    Int_t *subpfof,*subngroup;
    Int_t *subnumingroup,**subpglist;
    Int_t **subsubnumingroup, ***subsubpglist;
    Int_t *numcores,*coreflag;
    Int_t *subpfofold;
    Coordinate *gvel;
    Matrix *gveldisp;
    KDTree *tree;
    GridCell *grid;
    Coordinate cm,cmvel;
    //variables to keep track of structure level, pfof values (ie group ids) and their parent structure
    //use to point to current level
    StrucLevelData *pcsld;
    //use to store total number in sublevel;
    Int_t ns;
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif
    cout<<ThisTask<<" Beginning substructure search "<<endl;
    if (ngroup>0) {
    //point to current structure level
    pcsld=psldata;
    //ns=0;

    //if a single halo was not passed then store fof halos for unbinding purposes
    if (!opt.iSingleHalo) nhalos=ngroup;

    nsubsearch=ngroup;sublevel=1;ngroupidoffset=ngroupidoffsetold=0;
    if (opt.iBaryonSearch>=1 && opt.partsearchtype==PSTALL) numingroup=BuildNumInGroupTyped(nsubset, ngroup, pfof, Partsubset.data(), DARKTYPE);
    else numingroup=BuildNumInGroup(nsubset, ngroup, pfof);
    //since initially groups in order find index of smallest group that can be searched for substructure
    //for (Int_t i=1;i<=ngroup;i++) if (numingroup[i]<MINSUBSIZE) {nsubsearch=i-1;break;}
    //for (Int_t i=1;i<=ngroup;i++) if (numingroup[i]<MINCELLSIZE) {nsubsearch=i-1;break;}
    firstgroup=1;
    firstgroupoffset=0;
    //if keeping fof as a leve, then don't start substructure for all 3dfof halos and move the structure level pointer upwards
    ///\todo check offsets of firstgroup and old
    if (opt.iKeepFOF) {
        firstgroup=opt.num3dfof+1;
        ngroupidoffsetold=opt.num3dfof;
        pcsld=psldata->nextlevel;
        nsubsearch=ngroup-opt.num3dfof;
    }
    for (Int_t i=firstgroup;i<=ngroup;i++) if (numingroup[i]<opt.MinSize*2) {nsubsearch=i-firstgroup;break;}
    iflag=(nsubsearch>0);

    if (iflag) {
    if (opt.iBaryonSearch>=1 && opt.partsearchtype==PSTALL) pglist=BuildPGListTyped(nsubset, ngroup, numingroup, pfof,Partsubset.data(),DARKTYPE);
    else pglist=BuildPGList(nsubset, ngroup, numingroup, pfof);
    //now store group ids of (sub)structures that will be searched for (sub)substructure.
    //since at level zero, the particle group list that is going to be used to calculate the background, outliers and searched through is simple pglist here
    //also the group size is simple numingroup
    subnumingroup=new Int_t[nsubsearch+1];
    subpglist=new Int_t*[nsubsearch+1];
    for (Int_t i=1;i<=nsubsearch;i++) {
        subnumingroup[i]=numingroup[i+firstgroupoffset];
        subpglist[i]=new Int_t[subnumingroup[i]];
        for (Int_t j=0;j<subnumingroup[i];j++) subpglist[i][j]=pglist[i+firstgroupoffset][j];
    }
    for (Int_t i=1;i<=ngroup;i++) delete[] pglist[i];
    delete[] pglist;
    delete[] numingroup;
    //now start searching while there are still sublevels to be searched

        while (iflag) {
        if (opt.iverbose) cout<<ThisTask<<" There are "<<nsubsearch<<" substructures large enough to search for other substructures at sub level "<<sublevel<<endl;
        oldnsubsearch=nsubsearch;
        subsubnumingroup=new Int_t*[nsubsearch+1];
        subsubpglist=new Int_t**[nsubsearch+1];
        subngroup=new Int_t[nsubsearch+1];
        numcores=new Int_t[nsubsearch+1];
        subpfofold=new Int_t[nsubsearch+1];
        ns=0;
        // START: ENCAPSULATE
        //here loop over all sublevel groups that need to be searched for substructure
        // ADACS: this loop proceses halos (that are independent) sequentially. This is unecessary.
        // #pragma omp parallel for...
        // parallelise loop, collection, and/or calculations per collection
        // TODO: give a collection a pool of threads thread_pool
        for (Int_t i=1;i<=oldnsubsearch;i++) {
            subpfofold[i]=pfof[subpglist[i][0]];
            subPart=new Particle[subnumingroup[i]];
            for (Int_t j=0;j<subnumingroup[i];j++) subPart[j]=Partsubset[subpglist[i][j]];
            //now if low statistics, then possible that very central regions of subhalo will be higher due to cell size used and Nv search
            //so first determine centre of subregion
            // ADACS: here is an example of unecessary parallelisation in most cases
            // ADACS: (save for very high res zooms of individual objects containing billions of particles 
            Double_t cmx=0.,cmy=0.,cmz=0.,cmvelx=0.,cmvely=0.,cmvelz=0.;
            Double_t mtotregion=0.0;
            Int_t j;
            //suggested encapsulation by PJE
            if (opt.icmrefadjust) {
                //this routine is in substructureproperties.cxx. Has internal parallelisation
                GMatrix cmphase = CalcPhaseCM(subnumingroup[i], subPart);
                //this routine is within this file, also has internal parallelisation
                AdjustSubPartToPhaseCM(subnumingroup[i], subPart, cmphase);
            }
 
            //ADACS: for large objects, extra processing steps are requried
            //ADACS: Some of these subroutines make use of OpenMP. For this to continue
		    //ADACS: the pool of threads would have to be changed
            // PreCalcSearchSubSet()
            // subnumingroup, opt, subpart, sublevel, tree, ngrid, thistask, grid, gvel, gveldisp
            PreCalcSearchSubSet(subnumingroup, i, opt, subPart, sublevel, tree, ngrid, ThisTask, grid, gvel, gveldisp);
            //ADACS: Here the object is searched. Not much of this uses OpenMP but there are 
		    //  one or two subroutines called within SearchSubset that do make use of OpenMP. 

            subpfof=SearchSubset(opt,subnumingroup[i],subnumingroup[i],subPart,subngroup[i],sublevel,&numcores[i]);
            //now if subngroup>0 change the pfof ids of these particles in question and see if there are any substrucures that can be searched again.
            //the group ids must be stored along with the number of groups in this substructure that will be searched at next level.
            //now check if self bound and if not, id doesn't change from original subhalo,ie: subpfof[j]=0
            // SearchSubStruct(subngroup, i, ng, subsubnumingroup, subnumingroup, subpfof, subsubpglist, opt, 
            //     numcores, coreflag, iunbindflag, subPart, pfof, subpglist, ngroup, ngroupidoffset);
            if (subngroup[i]) {
                SearchSubStruct(numcores, i, opt, ng, iunbindflag, subnumingroup,
                    subPart, subngroup, subpfof, subsubnumingroup, subsubpglist,
                    pfof, subpglist, ngroup, ngroupidoffset);
                
                ngroupidoffset+=subngroup[i];
            }
            
            delete[] subpfof; // freeze
            delete[] subPart; // freeze
            //increase tot num of objects at sublevel
            //ADACS: this would need a reduction at the end.
            ns+=subngroup[i]; // reduction
        }
        // END: ENCAPSULATE

        //if objects have been found adjust the StrucLevelData
        //this stores the address of the parent particle and pfof along with child substructure particle and pfof
        if (ns>0) {
            pcsld->nextlevel=new StrucLevelData();
            //pcsld->nextlevel->Initialize(ns);
            pcsld->nextlevel->Phead=new Particle*[ns+1];
            pcsld->nextlevel->gidhead=new Int_t*[ns+1];
            pcsld->nextlevel->Pparenthead=new Particle*[ns+1];
            pcsld->nextlevel->gidparenthead=new Int_t*[ns+1];
            pcsld->nextlevel->giduberparenthead=new Int_t*[ns+1];
            pcsld->nextlevel->stypeinlevel=new Int_t[ns+1];
            pcsld->nextlevel->stype=HALOSTYPE+10*sublevel;
            pcsld->nextlevel->nsinlevel=ns;
            Int_t nscount=1;
            for (Int_t i=1;i<=oldnsubsearch;i++) {
                Int_t ii=0,iindex;
                Int_t *gidparentheadval,*giduberparentheadval;
                Particle *Pparentheadval;
                //here adjust head particle of parent structure if necessary. Search for first instance where
                //the pfof value of the particles originally associated with the parent structure have a value
                //less than the expected values for substructures
                while (pfof[subpglist[i][ii]]>ngroup+ngroupidoffset-ns && ii<subnumingroup[i]) ii++;
                //if a (sub)structure has been fully decomposed into (sub)substructures then possible no particles
                //remaining with a halo id, this must be handled
                if (pfof[subpglist[i][ii]]>ngroup+ngroupidoffset-ns) {
                    //in the case that no particles remain part of the parent (sub)structure
                    //then remove the (sub)structure from the current structure level, add the new structures to the next structure level
                    //first find index of (sub)structure
                    Particle *Pval=&Partsubset[subpglist[i][0]];
                    iindex=0;
                    while (pcsld->Phead[iindex++]!=Pval);
                    //store the parent/uber parent info if not a field halo
                    if(sublevel==1&&opt.iKeepFOF==0) {
                        Pparentheadval=NULL;
                        gidparentheadval=NULL;
                        giduberparentheadval=NULL;
                    }
                    else {
                        Pparentheadval=pcsld->Phead[iindex];
                        gidparentheadval=pcsld->gidhead[iindex];
                        giduberparentheadval=pcsld->giduberparenthead[iindex];
                    }
                    //now copy information
                    for (Int_t jj=iindex;jj<pcsld->nsinlevel;jj++) {
                        pcsld->Phead[jj]=pcsld->Phead[jj+1];
                        pcsld->gidhead[jj]=pcsld->gidhead[jj+1];
                        pcsld->Pparenthead[jj]=pcsld->Pparenthead[jj+1];
                        pcsld->gidparenthead[jj]=pcsld->gidparenthead[jj+1];
                        pcsld->giduberparenthead[jj]=pcsld->giduberparenthead[jj+1];
                        pcsld->stypeinlevel[jj]=pcsld->stypeinlevel[jj+1];
                    }
                    //decrease number of structures in level, adjust ids, offsets, etc
                    for (Int_t jj=0;jj<nsubset;jj++) if (pfof[jj]>subpfofold[i]) pfof[jj]--;
                    pcsld->nsinlevel--;
                    ngroupidoffsetold--;
                }
                else {
                    //otherwise, just a matter of updating some pointers
                    //store index in the structure list to access the parent (sub)structure
                    //iindex=pfof[subpglist[i][ii]]-ngroupidoffsetold-firstgroupoffset;
                    //don't need to offset group as already taken care of.
                    iindex=pfof[subpglist[i][ii]]-ngroupidoffsetold;
                    pcsld->gidhead[iindex]=&pfof[subpglist[i][ii]];
                    pcsld->Phead[iindex]=&Partsubset[subpglist[i][ii]];
                    //only for field haloes does the gidparenthead and giduberparenthead need to be adjusted
                    //but only if 3DFOFs are not kept as uber parents
                    if(sublevel==1&&opt.iKeepFOF==0) {
                        pcsld->gidparenthead[iindex]=&pfof[subpglist[i][ii]];
                        pcsld->giduberparenthead[iindex]=&pfof[subpglist[i][ii]];
                    }
                    Pparentheadval=&Partsubset[subpglist[i][ii]];
                    gidparentheadval=pcsld->gidhead[iindex];
                    giduberparentheadval=pcsld->giduberparenthead[iindex];
                }
                for (Int_t j=1;j<=subngroup[i];j++) {
                    //need to restructure pointers so that they point to what the parent halo
                    //points to and point to the head of the structure
                    //this is after viable structures have been identfied
                    pcsld->nextlevel->Phead[nscount]=&Partsubset[subsubpglist[i][j][0]];
                    pcsld->nextlevel->gidhead[nscount]=&pfof[subsubpglist[i][j][0]];
                    pcsld->nextlevel->Pparenthead[nscount]=Pparentheadval;
                    pcsld->nextlevel->gidparenthead[nscount]=gidparentheadval;
                    pcsld->nextlevel->giduberparenthead[nscount]=giduberparentheadval;
                    //if multiple core region then set the appropriate structure level
                    if (j>=(subngroup[i]-numcores[i])+1) pcsld->nextlevel->stypeinlevel[nscount]=HALOSTYPE+10*(sublevel-1)+HALOCORESTYPE;
                    else pcsld->nextlevel->stypeinlevel[nscount]=HALOSTYPE+10*sublevel;
                    nscount++;
                }
            }
            ngroupidoffsetold+=pcsld->nsinlevel;
            pcsld=pcsld->nextlevel;
        }
        if (opt.iverbose) cout<<ThisTask<<"Finished searching substructures to sublevel "<<sublevel<<endl;
        sublevel++;
        for (Int_t i=1;i<=oldnsubsearch;i++) delete[] subpglist[i];
        delete[] subpglist;
        delete[] subnumingroup;
        nsubsearch=0;
        //after looping over all level sublevel substructures adjust nsubsearch, set subpglist subnumingroup, so that can move to next level.
        for (Int_t i=1;i<=oldnsubsearch;i++)
            for (Int_t j=1;j<=subngroup[i];j++) if (subsubnumingroup[i][j]>MINSUBSIZE)
                nsubsearch++;
        if (nsubsearch>0) {
            subnumingroup=new Int_t[nsubsearch+1];
            subpglist=new Int_t*[nsubsearch+1];
            nsubsearch=1;
            for (Int_t i=1;i<=oldnsubsearch;i++) {
                for (Int_t j=1;j<=subngroup[i];j++)
                    if (subsubnumingroup[i][j]>MINSUBSIZE) {
                        subnumingroup[nsubsearch]=subsubnumingroup[i][j];
                        subpglist[nsubsearch]=new Int_t[subnumingroup[nsubsearch]];
                        for (Int_t k=0;k<subnumingroup[nsubsearch];k++) subpglist[nsubsearch][k]=subsubpglist[i][j][k];
                        nsubsearch++;
                    }
            }
            nsubsearch--;
        }
        else iflag=false;
        //free memory
        for (Int_t i=1;i<=oldnsubsearch;i++) {
            if (subngroup[i]>0) {
                for (Int_t j=1;j<=subngroup[i];j++) delete[] subsubpglist[i][j];
                delete[] subsubnumingroup[i];
                delete[] subsubpglist[i];
            }
        }
        delete[] subsubnumingroup;
        delete[] subsubpglist;
        delete[] subngroup;
        delete[] numcores;
        delete[] subpfofold;
        if (opt.iverbose) cout<<ThisTask<<"Finished storing next level of substructures to be searched for subsubstructure"<<endl;
    }

    ngroup+=ngroupidoffset;
    cout<<ThisTask<<"Done searching substructure to "<<sublevel-1<<" sublevels "<<endl;
    }
    else delete[] numingroup;
    //if not an idividual halo and want bound haloes after substructure search (and not searching for baryons afterwards)
    //unbind halo population
    if (!opt.iSingleHalo&&opt.iBoundHalos>1&&!opt.iBaryonSearch) {
        //begin by storing information of the current hierarchy
        Int_t nhaloidoffset=0,nhierarchy,gidval;
        Int_t *nsub,*parentgid,*uparentgid,*stype;
        StrucLevelData *ppsldata,**papsldata;

        ng=nhalos;
        numingroup=BuildNumInGroup(nsubset, ngroup, pfof);
        pglist=BuildPGList(nsubset, ngroup, numingroup, pfof);

        //store the parent and uber parent info
        nsub=new Int_t[ngroup+1];
        parentgid=new Int_t[ngroup+1];
        uparentgid=new Int_t[ngroup+1];
        stype=new Int_t[ngroup+1];
        GetHierarchy(opt,ngroup,nsub,parentgid,uparentgid,stype);
        //store pointer to the various substructure levels
        ppsldata=psldata;
        nhierarchy=0;
        //determine the number of levels in the hierarchy
        while (ppsldata->nextlevel!=NULL){nhierarchy++;ppsldata=ppsldata->nextlevel;}
        ppsldata=psldata;
        papsldata=new StrucLevelData*[nhierarchy];
        nhierarchy=0;
        while (ppsldata!=NULL) {papsldata[nhierarchy++]=ppsldata;ppsldata=ppsldata->nextlevel;}

        if(CheckUnboundGroups(opt,nsubset,Partsubset.data(),nhalos,pfof,numingroup,pglist,0)) {
            //if haloes adjusted then need to update the StrucLevelData
            //first update just halos (here ng=old nhalos)
            //by setting NULL values in structure level and moving all the unbound halos the end of array
            for (Int_t i=1;i<=ng;i++) {
                if (numingroup[i]==0) {
                    psldata->Phead[i]=NULL;
                    psldata->gidhead[i]=NULL;
                }
                else {
                    psldata->Phead[i]=&Partsubset[pglist[i][0]];
                    psldata->gidhead[i]=&pfof[pglist[i][0]];
                    //set the parent pointers to appropriate addresss such that the parent and uber parent are the same as the groups head
                    psldata->gidparenthead[i]=psldata->gidhead[i];
                    psldata->giduberparenthead[i]=psldata->gidhead[i];
                }
            }
            for (Int_t i=1,j=0;i<=nhalos;i++) {
                if (psldata->Phead[i]==NULL) {
                    j=i+1;while(psldata->Phead[j]==NULL) j++;
                    //copy the first non-NULL pointer to current NULL pointers,
                    //ie: for a non-existing structure which has been unbound, remove it from the structure list
                    //by copying the pointers of the next still viable structure to that address and setting
                    //the pointers at the new position, j, to NULL
                    psldata->Phead[i]=psldata->Phead[j];
                    psldata->gidhead[i]=psldata->gidhead[j];
                    psldata->gidparenthead[i]=psldata->gidparenthead[j];
                    psldata->giduberparenthead[i]=psldata->giduberparenthead[j];
                    psldata->stypeinlevel[i]=psldata->stypeinlevel[j];
                    psldata->Phead[j]=NULL;
                    psldata->gidhead[j]=NULL;
                    psldata->gidparenthead[j]=NULL;
                    psldata->giduberparenthead[j]=NULL;
                    psldata->stypeinlevel[j]=BGTYPE;
                }
            }
            psldata->nsinlevel=nhalos;
            //then adjust sublevels so that they point to the appropriate parent and uberparent values
            for (Int_t i=nhierarchy-1;i>0;i--){
                for (int j=1;j<=papsldata[i]->nsinlevel;j++) {
                    gidval=(*papsldata[i]->gidhead[j]);
                    if (parentgid[gidval]!=GROUPNOPARENT) {
                        if (numingroup[parentgid[gidval]]>0) papsldata[i]->gidparenthead[j]=&pfof[pglist[parentgid[gidval]][0]];
                        else papsldata[i]->gidparenthead[j]=NULL;
                    }
                    if (uparentgid[gidval]!=GROUPNOPARENT) {
                        if (numingroup[uparentgid[gidval]]>0) papsldata[i]->giduberparenthead[j]=&pfof[pglist[uparentgid[gidval]][0]];
                        else papsldata[i]->giduberparenthead[j]=NULL;
                    }
                }
            }
            //if inclusive halo masses calculated then also need to adjust the ordering of these properties
            if (opt.iInclusiveHalo) ReorderInclusiveMasses(ng,nhalos,numingroup,pdata);
            //adjust halo ids
            ReorderGroupIDs(ng,nhalos, numingroup, pfof,pglist);
            nhaloidoffset=ng-nhalos;
            for (Int_t i=0;i<nsubset;i++) if (pfof[i]>ng) pfof[i]-=nhaloidoffset;
        }
        for (Int_t i=1;i<=ngroup;i++) delete[] pglist[i];
        delete[] pglist;
        delete[] numingroup;
        delete[] nsub;
        delete[] parentgid;
        delete[] uparentgid;
        delete[] papsldata;
        ngroup-=nhaloidoffset;
    }
    }
    //update the number of local groups found
#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(&ngroup, 1, MPI_Int_t, mpi_ngroups, 1, MPI_Int_t, MPI_COMM_WORLD);
    cout<<ThisTask<<" has found a total of "<<ngroup<<endl;
#endif
}

/*!
    Check significance of group using significance parameter.  A group is considered significant if (ave/expected ave) (and possibly (max-min)/varexpected] is significant relative to Poisson noise.
    If a group is not start removing particle with lowest ell value.
*/
int CheckSignificance(Options &opt, const Int_t nsubset, Particle *Partsubset, Int_t &numgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist)
{
    Int_t i;
    Double_t ellaveexp,ellvarexp, ellvallim;
    Double_t *aveell,*varell,*maxell,*minell, *ptypefrac, *betaave;
    Int_t iflag=0,ng=numgroups;
    aveell=new Double_t[numgroups+1];
    maxell=new Double_t[numgroups+1];
    minell=new Double_t[numgroups+1];
    //varell=new Double_t[numgroups+1];
    betaave=new Double_t[numgroups+1];

    /*
    if (opt.iiterflag) ellvallim=opt.ellthreshold*opt.ellfac;
    else ellvallim=opt.ellthreshold;
    */
    ellvallim=opt.ellthreshold;
    ellaveexp=sqrt(2.0/M_PI)*exp(-ellvallim*ellvallim)*exp(0.5*ellvallim*ellvallim)/(1.0-gsl_sf_erf(ellvallim/sqrt(2.0)));

    //first adjust system and store pList for each group so one can access pfof appropriately
    if (opt.iverbose) {
        cout<<"Checking that groups have a significance level of "<<opt.siglevel<<" and contain more than "<<opt.MinSize<<" members"<<endl;
    }
    for (i=0;i<=numgroups;i++) {aveell[i]=0.;maxell[i]=-MAXVALUE;minell[i]=MAXVALUE;//numingroup[i]=0;
    }
    for (i=0;i<nsubset;i++)
    {
        Int_t gid=pfof[Partsubset[i].GetID()];
        Double_t ellvalue=Partsubset[i].GetPotential();
        aveell[gid]+=ellvalue;
        if(maxell[gid]<ellvalue)maxell[gid]=ellvalue;
        if(minell[gid]>ellvalue)minell[gid]=ellvalue;
    }
    for (i=1;i<=numgroups;i++) {
        aveell[i]/=(Double_t)numingroup[i];
        betaave[i]=(aveell[i]/ellaveexp-1.0)*sqrt((Double_t)numingroup[i]);
        //flag indicating that group ids need to be adjusted
        if(betaave[i]<opt.siglevel) iflag=1;
    }
    if (opt.iverbose) cout<<"Done"<<endl;
    if (iflag){
        if (opt.iverbose) cout<<"Remove groups below significance level"<<endl;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i)
{
#pragma omp for
#endif
        for (i=1;i<=numgroups;i++) {
            if(betaave[i]<opt.siglevel) {
                do {
                    if ((numingroup[i])<opt.MinSize) {
                        for (Int_t j=0;j<numingroup[i];j++) pfof[Partsubset[pglist[i][j]].GetID()]=0;
                        numingroup[i]=-1;
                        break;
                    }
                    Int_t iminell;
                    Double_t vminell=minell[i];
                    aveell[i]=(aveell[i]*(Double_t)numingroup[i]-minell[i])/(Double_t)(numingroup[i]-1.0);
                    minell[i]=maxell[i];
                    for (Int_t j=0;j<numingroup[i];j++)
                        if(Partsubset[pglist[i][j]].GetPotential()>vminell&&Partsubset[pglist[i][j]].GetPotential()<minell[i])minell[i]=Partsubset[pglist[i][j]].GetPotential();
                        else if (Partsubset[pglist[i][j]].GetPotential()==vminell) iminell=j;
                    pfof[Partsubset[pglist[i][iminell]].GetID()]=0;
                    if (iminell!=numingroup[i]-1)pglist[i][iminell]=pglist[i][--numingroup[i]];
                    betaave[i]=(aveell[i]/ellaveexp-1.0)*sqrt((Double_t)numingroup[i]);
                } while(betaave[i]<opt.siglevel);
            }
            if ((numingroup[i])<opt.MinSize) {
                for (Int_t j=0;j<numingroup[i];j++) pfof[Partsubset[pglist[i][j]].GetID()]=0;
                numingroup[i]=-1;
            }
        }
#ifdef USEOPENMP
}
#endif
        if (opt.iverbose) cout<<"Done"<<endl;
        for (i=1;i<=numgroups;i++) if (numingroup[i]==-1) ng--;
        if (ng) ReorderGroupIDs(numgroups, ng, numingroup, pfof, pglist, Partsubset);
        else if (opt.iverbose) printf("No groups of significance found\n");
    }

    //free memory
    delete[] aveell;
    delete[] maxell;
    delete[] minell;
    delete[] betaave;

    numgroups=ng;
    return iflag;
}

//@}

/// \name Routines searches baryonic or other components separately based on initial dark matter (or other) search
//@{
/*!
 * Searches star and gas particles separately to see if they are associated with any dark matter particles belonging to a substructure
 *
 * Assumes baryonic particles that could be associated with haloes located in MPI's thread are also on that thread if all particles were
 * searched initial for field objects.
 * Otherwise, sends dark matter particles that overlap other mpi domains across to see if DM particle with which baryon particle associated with
 * lies on another mpi domain.
 *
 * \todo might use full phase-space tensor association.
*/
Int_t* SearchBaryons(Options &opt, Int_t &nbaryons, Particle *&Pbaryons, const Int_t ndark, vector<Particle> &Part, Int_t *&pfofdark, Int_t &ngroupdark, Int_t &nhalos, int ihaloflag, int iinclusive, PropData *pdata)
{
    KDTree *tree;
    Double_t *period;
    Int_t *pfofbaryons, *pfofall, *pfofold;
    Int_t i,pindex,npartingroups,ng,nghalos,nhalosold=nhalos, baryonfofold;
    Int_t *ids, *storeval,*storeval2;
    Double_t D2,dval,rval;
    Coordinate x1;
    Particle p1;
    int icheck;
    FOFcompfunc fofcmp;
    Double_t param[20];
    int nsearch=opt.Nvel;
    Int_t *nnID=NULL,*numingroup;
    Double_t *dist2=NULL, *localdist;
    int nthreads=1,maxnthreads,tid;
    int minsize;
    Int_t nparts=ndark+nbaryons;
    Int_t nhierarchy=1,gidval;
    StrucLevelData *ppsldata,**papsldata;
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif

    cout<<ThisTask<<" search baryons "<<nparts<<" "<<ndark<<endl;
    //if searched all particles in FOF, reorder particles and also the pfof group id array
    if (opt.partsearchtype==PSTALL) {
        cout<<" of only baryons in FOF structures as baryons have already been grouped in FOF search "<<endl;
        //store original pfof value
        pfofall=pfofdark;
#ifndef USEMPI
        //if no objects are found then stop (or if no substructures found as baryons won't change association)
        if (ngroupdark==0 || ngroupdark==nhalos) return pfofall;
#endif
        pfofbaryons=&pfofall[ndark];
        storeval=new Int_t[nparts];
        storeval2=new Int_t[nparts];
        for (i=0;i<nparts;i++) {
            if (Part[i].GetType()==DARKTYPE) Part[i].SetType(-1);
            storeval2[i]=Part[i].GetPID();
            Part[i].SetPID(pfofdark[i]);
        }
        qsort(Part.data(),nparts,sizeof(Particle),TypeCompare);
        //sort(Part.begin(),Part.end(),TypeCompareVec);
        Pbaryons=&Part[ndark];
        for (i=0;i<nparts;i++) {
            //store id order after type sort
            storeval[i]=Part[i].GetID();
            Part[i].SetID(i);
            pfofdark[i]=Part[i].GetPID();
        }
    }
#ifdef USEOPENMP
#pragma omp parallel
    {
    if (omp_get_thread_num()==0) maxnthreads=omp_get_num_threads();
    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
#endif
    //if serial can allocate pfofall at this point as the number of particles in local memory will not change
#ifndef USEMPI
    if (opt.partsearchtype!=PSTALL) {
        pfofall=new Int_t[nbaryons+ndark];
        for (i=0;i<ndark;i++) pfofall[i]=pfofdark[i];
        pfofbaryons=&pfofall[ndark];
        for (i=0;i<nbaryons;i++) pfofbaryons[i]=0;
        if (ngroupdark==0) {
            delete[] storeval;
            delete[] storeval2;
            return pfofall;
        }
    }
#else
    //if using mpi but all particle FOF search, then everything is localized
    if (opt.partsearchtype!=PSTALL) {
        pfofbaryons=new Int_t[nbaryons];
        for (i=0;i<nbaryons;i++) pfofbaryons[i]=0;
        localdist=new Double_t[nbaryons];
        for (i=0;i<nbaryons;i++) localdist[i]=MAXVALUE;
    }
#endif

    //build a particle list containing only dark matter particles
    numingroup=BuildNumInGroup(ndark, ngroupdark, pfofdark);

    //determine total number of particles in groups and set search appropriately
    npartingroups=0;
    for (i=1;i<=ngroupdark;i++) npartingroups+=numingroup[i];
    if (npartingroups<=2*opt.MinSize) nsearch=npartingroups-1;
    else nsearch=2*opt.MinSize;
    if (opt.p>0) {
        period=new Double_t[3];
        for (int j=0;j<3;j++) period[j]=opt.p;
    }

    //sort dark matter particles so that particles belonging to a group are first, then all other dm particles
    if (opt.iverbose) cout<<"sort particles so that tree only uses particles in groups "<<npartingroups<<endl;

    //search all dm particles in structures
    for (i=0;i<ndark;i++) Part[i].SetPotential(2*(pfofdark[i]==0)+(pfofdark[i]>1));
    qsort(Part.data(), ndark, sizeof(Particle), PotCompare);
    //sort(Part.begin(),Part.begin()+ndark,PotCompareVec);
    ids=new Int_t[ndark+1];
    //store the original order of the dark matter particles
    for (i=0;i<ndark;i++) ids[i]=Part[i].GetID();

    if (npartingroups>0) {
    //set parameters, first to some fraction of halo linking length
    param[1]=(opt.ellxscale*opt.ellxscale)*(opt.ellphys*opt.ellphys)*(opt.ellhalophysfac*opt.ellhalophysfac);
    param[6]=param[1];
    //also check to see if velocity scale still zero find dispersion of largest halo
    //otherwise search uses the largest average "local" velocity dispersion of halo identified in the mpi domain.
    if (opt.HaloVelDispScale==0) {
        Double_t mtotregion,vx,vy,vz;
        Coordinate vmean;
        mtotregion=vx=vy=vz=0;
        for (i=0;i<numingroup[1];i++) {
            vx+=Part[i].GetVelocity(0)*Part[i].GetMass();
            vy+=Part[i].GetVelocity(1)*Part[i].GetMass();
            vz+=Part[i].GetVelocity(2)*Part[i].GetMass();
            mtotregion+=Part[i].GetMass();
        }
        vmean[0]=vx/mtotregion;vmean[1]=vy/mtotregion;vmean[2]=vz/mtotregion;
        for (i=0;i<numingroup[1];i++) {
            for (int j=0;j<3;j++) opt.HaloVelDispScale+=pow(Part[i].GetVelocity(j)-vmean[j],2.0)*Part[i].GetMass();
        }
        opt.HaloVelDispScale/=mtotregion;
        param[2]=opt.HaloVelDispScale;
    }
    else param[2]=opt.HaloVelDispScale*16.0;//here use factor of 4 in local dispersion //could remove entirely and just use global dispersion but this will over compensate.
    param[7]=param[2];

    //Set fof type
    fofcmp=&FOF6d;
    if (opt.iverbose) {
        cout<<"Baryon search "<<nbaryons<<endl;
        cout<<"FOF6D uses ellphys and ellvel.\n";
        cout<<"Parameters used are : ellphys="<<sqrt(param[6])<<" Lunits, ellvel="<<sqrt(param[7])<<" Vunits.\n";
        cout<<"Building tree to search dm containing "<<npartingroups<<endl;
    }
    //build tree of baryon particles (in groups if a full particle search was done, otherwise npartingroups=nbaryons
    tree=new KDTree(Part.data(),npartingroups,nsearch/2,tree->TPHYS,tree->KEPAN,100,0,0,0,period);
    //allocate memory for search
    //find the closest dm particle that belongs to the largest dm group and associate the baryon with that group (including phase-space window)
    if (opt.iverbose) cout<<"Searching ..."<<endl;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,tid,p1,pindex,x1,D2,dval,rval,icheck,nnID,dist2,baryonfofold)
{
    nnID=new Int_t[nsearch];
    dist2=new Double_t[nsearch];
#pragma omp for
#else
    nnID=new Int_t[nsearch];
    dist2=new Double_t[nsearch];
#endif
    for (i=0;i<nbaryons;i++)
    {
#ifdef USEOPENMP
        tid=omp_get_thread_num();
#else
        tid=0;
#endif

        //if all particles have been searched for field objects then ignore baryons not associated with a group
        if (opt.partsearchtype==PSTALL && pfofbaryons[i]==0) continue;
        p1=Pbaryons[i];
        x1=Coordinate(p1.GetPosition());
        rval=dval=MAXVALUE;
        baryonfofold=pfofbaryons[i];
        tree->FindNearestPos(x1, nnID, dist2,nsearch);
        if (dist2[0]<param[6]) {
        for (int j=0;j<nsearch;j++) {
            D2=0;
            pindex=ids[Part[nnID[j]].GetID()];
            //determine if baryonic particle needs to be searched. Note that if all particles have been searched during FOF
            //then particle is checked regardless. If that is not the case, particle is searched only if its current group
            //is smaller than the group of the dm particle being examined.
            //But do not allow baryons to switch between fof structures
            if (opt.partsearchtype==PSTALL) icheck=((pfofdark[pindex]>nhalos)||(pfofdark[pindex]==baryonfofold));
            else icheck=(numingroup[pfofbaryons[i]]<numingroup[pfofdark[pindex]]);
            if (icheck) {
                if (fofcmp(p1,Part[nnID[j]],param)) {
                    for (int k=0;k<3;k++) {
                        D2+=(p1.GetPosition(k)-Part[nnID[j]].GetPosition(k))*(p1.GetPosition(k)-Part[nnID[j]].GetPosition(k))/param[6]+(p1.GetVelocity(k)-Part[nnID[j]].GetVelocity(k))*(p1.GetVelocity(k)-Part[nnID[j]].GetVelocity(k))/param[7];
                    }
                    //if gas thermal properties stored then also add self-energy to distance measure
#ifdef GASON
                    D2+=p1.GetU()/param[7];
#endif
                    //check to see if phase-space distance is small
                    if (dval>D2) {
                        dval=D2;pfofbaryons[i]=pfofdark[pindex];
                        rval=dist2[j];
#ifdef USEMPI
                        if (opt.partsearchtype!=PSTALL) localdist[i]=dval;
#endif
                    }
                }
            }
        }
        }
    }
    delete[] nnID;
    delete[] dist2;
#ifdef USEOPENMP
}
#endif
    }

#ifdef USEMPI
    //if mpi then baryons are not necessarily local if opt.partsearchtype!=PSTALL
    //in that case must search other mpi domains.
    //if all particles are searched then just need to reset the particle order
    ///\todo need to update this for mpi vector
    if (opt.partsearchtype!=PSTALL) {
        if (opt.iverbose) cout<<ThisTask<<" finished local search"<<endl;
        MPI_Barrier(MPI_COMM_WORLD);
        //determine all tagged dark matter particles that have search areas that overlap another mpi domain
#ifdef SWIFTINTERFACE
        MPIGetExportNumUsingMesh(libvelociraptorOpt, npartingroups, Part.data(), sqrt(param[1]));
#else
        MPIGetExportNum(npartingroups, Part.data(), sqrt(param[1]));
#endif
        //to store local mpi task
        mpi_foftask=MPISetTaskID(nbaryons);
        //then determine export particles, declare arrays used to export data
        PartDataIn = new Particle[NExport+1];
        PartDataGet = new Particle[NImport+1];
        FoFDataIn = new fofdata_in[NExport+1];
        FoFDataGet = new fofdata_in[NImport+1];
        //exchange particles

        MPIBuildParticleExportBaryonSearchList(npartingroups, Part.data(), pfofdark, ids, numingroup, sqrt(param[1]));

        //now dark matter particles associated with a group existing on another mpi domain are local and can be searched.
        NExport=MPISearchBaryons(nbaryons, Pbaryons, pfofbaryons, numingroup, localdist, nsearch, param, period);

        //reset order
        delete tree;
        for (i=0;i<ndark;i++) Part[i].SetID(ids[i]);
        qsort(Part.data(), ndark, sizeof(Particle), IDCompare);
        //sort(Part.begin(), Part.begin()+ndark, IDCompareVec);
        delete[] ids;

        //reorder local particle array and delete memory associated with Head arrays, only need to keep Particles, pfof and some id and idexing information
        delete[] FoFDataIn;
        delete[] FoFDataGet;
        delete[] PartDataIn;
        delete[] PartDataGet;

        Int_t newnbaryons=MPIBaryonGroupExchange(nbaryons,Pbaryons,pfofbaryons);
        //once baryons are correctly associated to the appropriate mpi domain and are local (either in Pbaryons or in the \ref fofid_in structure, specifically FOFGroupData arrays) must then copy info correctly.
//#ifdef MPIREDUCEMEM
        if (Nmemlocalbaryon<newnbaryons)
//#endif
        {
        //note that if mpireduce is not set then all info is copied into the FOFGroupData structure and must deallocated and reallocate Pbaryon array
            delete[] Pbaryons;
            Pbaryons=new Particle[newnbaryons];
            delete[] pfofbaryons;
            pfofbaryons=new Int_t[newnbaryons];
        }
        //then compile groups and if inclusive halo masses not calculated, reorder group ids
        MPIBaryonCompileGroups(newnbaryons,Pbaryons,pfofbaryons,opt.MinSize,(opt.iInclusiveHalo==0));
        delete[] mpi_foftask;
        if (opt.iverbose) cout<<ThisTask<<" finished search across domains"<<endl;
        //now allocate pfofall and store info
        pfofall=new Int_t[newnbaryons+ndark];
        for (i=0;i<ndark;i++) pfofall[i]=pfofdark[i];
        for (i=0;i<newnbaryons;i++) pfofall[i+ndark]=pfofbaryons[i];
        delete[] pfofbaryons;

        //update nbaryons
        nbaryons=newnbaryons;
        Nlocalbaryon[0]=newnbaryons;

        Part.resize(nparts);
        for (i=0;i<nbaryons;i++)Part[i+ndark]=Pbaryons[i];
        delete[] Pbaryons;
        Pbaryons=&Part.data()[ndark];
        for (i=0;i<nbaryons;i++) Pbaryons[i].SetID(i+ndark);
        Nlocal=nparts;
    } // end of if preliminary search is NOT all particles
    else {
        //reset order
        if (npartingroups>0) delete tree;
        for (i=0;i<ndark;i++) Part[i].SetID(ids[i]);
        qsort(Part.data(), ndark, sizeof(Particle), IDCompare);
        //sort(Part.begin(), Part.end(), IDCompareVec);
        delete[] ids;
        for (i=0;i<nparts;i++) {Part[i].SetPID(pfofall[Part[i].GetID()]);Part[i].SetID(storeval[i]);}
        qsort(Part.data(), nparts, sizeof(Particle), IDCompare);
        //sort(Part.begin(), Part.end(), IDCompareVec);
        for (i=0;i<nparts;i++) {
            pfofall[i]=Part[i].GetPID();Part[i].SetPID(storeval2[i]);
            if (Part[i].GetType()==-1)Part[i].SetType(DARKTYPE);
        }
        delete[] storeval;
        delete[] storeval2;
    }
//if NOT mpi
#else
    //now that search has finished, reset dark matter particles back to prior to search order
    //extra sorts are needed to reset the particles back to input order if all particles were searched initially
    if (opt.partsearchtype==PSTALL) {
        //reset order
        if (npartingroups>0) delete tree;
        for (i=0;i<ndark;i++) Part[i].SetID(ids[i]);
        qsort(Part.data(), ndark, sizeof(Particle), IDCompare);
        //sort(Part.begin(), Part.begin()+ndark, IDCompareVec);
        delete[] ids;
        for (i=0;i<nparts;i++) {Part[i].SetPID(pfofall[Part[i].GetID()]);Part[i].SetID(storeval[i]);}
        qsort(Part.data(), nparts, sizeof(Particle), IDCompare);
        //sort(Part.begin(), Part.end(), IDCompareVec);
        for (i=0;i<nparts;i++) {
            pfofall[i]=Part[i].GetPID();Part[i].SetPID(storeval2[i]);
            if (Part[i].GetType()==-1)Part[i].SetType(DARKTYPE);
        }
        delete[] storeval;
        delete[] storeval2;
    }
    else {
        delete tree;
        for (i=0;i<ndark;i++) Part[i].SetID(ids[i]);
        qsort(Part.data(), ndark, sizeof(Particle), IDCompare);
        //sort(Part.begin(), Part.begin()+ndark, IDCompareVec);
        delete[] ids;
        for (i=0;i<nbaryons;i++) Pbaryons[i].SetID(i+ndark);
    }
#endif
    //and free up memory
    delete[] numingroup;

    cout<<"Done"<<endl;

    //if unbinding go back and redo full unbinding after including baryons
    //note that here if all particles are searched, the particle array has been reorderd from the input order
    //to dark matter first followed by baryons as has the pfofall array
    if (opt.uinfo.unbindflag&&ngroupdark>0) {
        if (opt.iverbose) cout<<ThisTask<<" starting unbind of dm+baryons"<<endl;
        //build new arrays based on pfofall.
        Int_t *ningall,**pglistall;
        ningall=BuildNumInGroup(nparts, ngroupdark, pfofall);
        pglistall=BuildPGList(nparts, ngroupdark, ningall, pfofall);
        //first repoint everything in the structure pointer to pfofall. See GetHierarchy for some guidance on how to parse the
        //the structure data
        nhierarchy=1;
        ppsldata=psldata;
        //determine the number of levels in the hierarchy
        while (ppsldata->nextlevel!=NULL){nhierarchy++;ppsldata=ppsldata->nextlevel;}
        ppsldata=psldata;
        papsldata=new StrucLevelData*[nhierarchy];
        nhierarchy=0;
        while (ppsldata!=NULL) {papsldata[nhierarchy++]=ppsldata;ppsldata=ppsldata->nextlevel;}
        //now parse hierarchy and repoint stuff to the pfofall pointer instead of pfof
        if (opt.partsearchtype!=PSTALL) {
            for (i=nhierarchy-1;i>=0;i--){
                for (int j=1;j<=papsldata[i]->nsinlevel;j++) {
                    gidval=(*papsldata[i]->gidhead[j]);
                    papsldata[i]->gidhead[j]=&pfofall[pglistall[gidval][0]];
                    if (papsldata[i]->gidparenthead[j]!=NULL) {
                        gidval=(*papsldata[i]->gidparenthead[j]);
                        if (numingroup[gidval]>0) papsldata[i]->gidparenthead[j]=&pfofall[pglistall[gidval][0]];
                        else papsldata[i]->gidparenthead[j]=NULL;
                    }
                    if (papsldata[i]->giduberparenthead[j]!=NULL) {
                        gidval=(*papsldata[i]->giduberparenthead[j]);
                        if (numingroup[gidval]>0) papsldata[i]->giduberparenthead[j]=&pfofall[pglistall[gidval][0]];
                    }
                    else papsldata[i]->giduberparenthead[j]=NULL;
                }
            }
        }
        //store old number of groups
        ng=ngroupdark;
        //if any structures have been changed need to update the structure pointers.
        //For now, do NOT reorder group ids based on number of particles that a group is composed of
        //as need to adjust the hierarchy data structure.
        //before unbinding store the parent and uberparent of an object
        Int_t *nsub,*parentgid,*uparentgid,*stype;
        nsub=new Int_t[ngroupdark+1];
        parentgid=new Int_t[ngroupdark+1];
        uparentgid=new Int_t[ngroupdark+1];
        stype=new Int_t[ngroupdark+1];
        GetHierarchy(opt,ngroupdark,nsub,parentgid,uparentgid,stype);
        //store the old pfof values of all so that in case particle unbound from a
        //substructure they are reassigned to the uber parent halo
        pfofold=new Int_t[nparts];
        for (i=0;i<nparts;i++) pfofold[i]=pfofall[i];
        //Only unbind subhalos if desired.
        if (opt.iBoundHalos==0) {
            //this can be done by setting the number in the group to the negative value
            //as values below zero will ignored
            for (i=1;i<=nhalos;i++) {
                ningall[i]=-ningall[i];
            }
        }
        if (CheckUnboundGroups(opt,nparts, Part.data(), ngroupdark, pfofall, ningall,pglistall,0)) {
            if (opt.iBoundHalos==0) {
                for (i=1;i<=nhalos;i++) {
                    ningall[i]=-ningall[i];
                }
            }
            //now if pfofall is zero but was a substructure reassign back to uber parent
            //so long as that uber parent still exists.
            for (i=0;i<nparts;i++)
            {
                if (pfofall[Part[i].GetID()]==0 && pfofold[Part[i].GetID()]>nhalos) {
                    if (ningall[uparentgid[pfofold[Part[i].GetID()]]]>0) {
                        pfofall[Part[i].GetID()]=uparentgid[pfofold[Part[i].GetID()]];
                        ningall[uparentgid[pfofold[Part[i].GetID()]]]++;
                    }
                }
            }
            //must rebuild pglistall
            for (i=1;i<=ng;i++) delete[] pglistall[i];
            delete[] pglistall;
            pglistall=BuildPGList(nparts, ng, ningall, pfofall);
            //now adjust the structure pointers after unbinding where groups are NOT reordered
            //first find groups that have been removed, tag their head as NULL
            //otherwise update the head, parent and uber parent
            int nlevel=0,ninleveloff=0,ii;
            for (i=1;i<=ng;i++) {
                if (i-ninleveloff>papsldata[nlevel]->nsinlevel) {
                    ninleveloff+=papsldata[nlevel]->nsinlevel;
                    nlevel++;
                }
                ii=i-ninleveloff;
                if (ningall[i]==0) {
                    papsldata[nlevel]->Phead[ii]=NULL;
                    papsldata[nlevel]->gidhead[ii]=NULL;
                    papsldata[nlevel]->gidparenthead[ii]=NULL;
                    papsldata[nlevel]->giduberparenthead[ii]=NULL;
                }
                else {
                    papsldata[nlevel]->Phead[ii]=&Part[pglistall[i][0]];
                    papsldata[nlevel]->gidhead[ii]=&pfofall[pglistall[i][0]];
                    if (parentgid[i]!=GROUPNOPARENT) {
                        if (ningall[parentgid[i]]>0) papsldata[nlevel]->gidparenthead[ii]=&pfofall[pglistall[parentgid[i]][0]];
                        else papsldata[nlevel]->gidparenthead[ii]=NULL;
                    }
                    if (uparentgid[i]!=GROUPNOPARENT) {
                        if (ningall[uparentgid[i]]>0) papsldata[nlevel]->giduberparenthead[ii]=&pfofall[pglistall[uparentgid[i]][0]];
                        else papsldata[nlevel]->giduberparenthead[ii]=NULL;
                    }
                }
            }
            //then clean up the structure list so that removed groups no longer in structure list
            for (i=nhierarchy-1;i>=0;i--) {
                Int_t ninlevel=0, iend=papsldata[i]->nsinlevel;
                for (Int_t k=1;k<=papsldata[i]->nsinlevel;k++) if (papsldata[i]->Phead[k]!=NULL) ninlevel++;
                for (Int_t k=1,j=0;k<=ninlevel;k++) {
                    if (papsldata[i]->Phead[k]==NULL) {
                        //j=k+1;while(papsldata[i]->Phead[j]==NULL) j++;
                        //find last entry that is not NULL
                        j=iend;while(papsldata[i]->Phead[j]==NULL) j--;
                        iend--;
                        //copy the first non-NULL pointer to current NULL pointers,
                        //ie: for a non-existing structure which has been unbound, remove it from the structure list
                        //by copying the pointers of the next still viable structure to that address and setting
                        //the pointers at the new position, j, to NULL
                        papsldata[i]->Phead[k]=papsldata[i]->Phead[j];
                        papsldata[i]->gidhead[k]=papsldata[i]->gidhead[j];
                        papsldata[i]->gidparenthead[k]=papsldata[i]->gidparenthead[j];
                        papsldata[i]->giduberparenthead[k]=papsldata[i]->giduberparenthead[j];
                        papsldata[i]->stypeinlevel[k]=papsldata[i]->stypeinlevel[j];
                        papsldata[i]->Phead[j]=NULL;
                        papsldata[i]->gidhead[j]=NULL;
                        papsldata[i]->gidparenthead[j]=NULL;
                        papsldata[i]->giduberparenthead[j]=NULL;
                        papsldata[i]->stypeinlevel[j]=BGTYPE;
                    }
                }
                papsldata[i]->nsinlevel=ninlevel;
            }
            if (opt.iverbose) cout<<ThisTask<<" Reorder after finding baryons and unbinding, previously had "<<ng<<" groups and now have "<<ngroupdark<<endl;
            //reorder group ids, keeping the ordering unchanged.
            map<Int_t, Int_t> remap;
            Int_t newng=0, oldpid, newpid;
            remap[0]=0;
            for (i=1;i<=ng;i++) {
                if (ningall[i]>0) {
                    newng++;
                    remap[i]=newng;
                }
                else  remap[i]=0;
            }
            for (i=0;i<nparts;i++)
            {
                oldpid = pfofall[Part[i].GetID()];
                if (oldpid==0) continue;
                newpid = remap[oldpid];
                if (oldpid==newpid) continue;
                pfofall[Part[i].GetID()] = newpid;
            }
            if (opt.iverbose) cout<<ThisTask<<" Done"<<endl;
            delete[] ningall;
            for (i=1;i<=ng;i++) delete[] pglistall[i];
            delete[] pglistall;
            for (i=nhierarchy-1;i>=0;i--) papsldata[i]=NULL;
            delete[] papsldata;
        }
        else {
            delete[] ningall;
            for (i=1;i<=ng;i++) delete[] pglistall[i];
            delete[] pglistall;
        }
        delete[] pfofold;
        delete[] nsub;
        delete[] parentgid;
        delete[] uparentgid;
        delete[] stype;
    }

#ifdef USEMPI
    //if number of groups has changed then update
    if (opt.uinfo.unbindflag) {
    cout<<"MPI thread "<<ThisTask<<" has found "<<ngroupdark<<endl;
    MPI_Allgather(&ngroupdark, 1, MPI_Int_t, mpi_ngroups, 1, MPI_Int_t, MPI_COMM_WORLD);
    //free up memory now that only need to store pfof and global ids
    if (ThisTask==0) {
        int totalgroups=0;
        for (int j=0;j<NProcs;j++) totalgroups+=mpi_ngroups[j];
        cout<<"Total number of groups found is "<<totalgroups<<endl;
    }
    }
#endif

    return pfofall;
}
//@}

/// \name Routines used to determine substructure hierarchy
//@{
Int_t GetHierarchy(Options &opt,Int_t ngroups, Int_t *nsub, Int_t *parentgid, Int_t *uparentgid, Int_t* stype)
{
    if (opt.iverbose) cout<<"Getting Hierarchy "<<ngroups<<endl;
    Int_t ng=0,nhierarchy=1,noffset=0;
    StrucLevelData *ppsldata,**papsldata;
    ppsldata=psldata;
    while (ppsldata->nextlevel!=NULL){nhierarchy++;ppsldata=ppsldata->nextlevel;}
    for (Int_t i=1;i<=ngroups;i++) nsub[i]=0;
    for (Int_t i=1;i<=ngroups;i++) parentgid[i]=GROUPNOPARENT;
    for (Int_t i=1;i<=ngroups;i++) uparentgid[i]=GROUPNOPARENT;
    ppsldata=psldata;
    papsldata=new StrucLevelData*[nhierarchy];
    nhierarchy=0;
    while (ppsldata!=NULL) {papsldata[nhierarchy++]=ppsldata;ppsldata=ppsldata->nextlevel;}
    for (int i=nhierarchy-1;i>=1;i--){
        //store number of substructures
        for (int j=1;j<=papsldata[i]->nsinlevel;j++) {
	           if (papsldata[i]->gidparenthead[j]!=NULL&&papsldata[i]->gidparenthead[j]!=papsldata[i]->gidhead[j]) nsub[*(papsldata[i]->gidparenthead[j])]++;
        }
        //then add these to parent substructure
        for (int j=1;j<=papsldata[i]->nsinlevel;j++) {
    	    if (papsldata[i]->gidparenthead[j]!=NULL&&papsldata[i]->gidparenthead[j]!=papsldata[i]->gidhead[j]){
                nsub[*(papsldata[i]->gidparenthead[j])]+=nsub[*(papsldata[i]->gidhead[j])];
                parentgid[*(papsldata[i]->gidhead[j])]=*(papsldata[i]->gidparenthead[j]);
            }
            if (papsldata[i]->giduberparenthead[j]!=NULL &&papsldata[i]->gidparenthead[j]!=papsldata[i]->gidhead[j]) uparentgid[*(papsldata[i]->gidhead[j])]=*(papsldata[i]->giduberparenthead[j]);
            stype[*(papsldata[i]->gidhead[j])]=(papsldata[i]->stypeinlevel[j]);
        }
    }
    //store field structures (top treel level) types
    for (int j=1;j<=papsldata[0]->nsinlevel;j++) {
        stype[*(papsldata[0]->gidhead[j])]=(papsldata[0]->stypeinlevel[j]);
    }
    for (int i=0;i<nhierarchy;i++)papsldata[i]=NULL;
    delete[] papsldata;
    if(opt.iverbose) cout<<"Done"<<endl;
    return nhierarchy;
}

void CopyHierarchy(Options &opt,PropData *pdata, Int_t ngroups, Int_t *nsub, Int_t *parentgid, Int_t *uparentgid, Int_t* stype)
{
    Int_t i,haloidoffset=0;
#ifdef USEMPI
    for (int j=0;j<ThisTask;j++)haloidoffset+=mpi_ngroups[j];
#endif
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i)
{
    #pragma omp for nowait
#endif
    for (i=1;i<=ngroups;i++) {
        pdata[i].haloid=opt.snapshotvalue+i;
        //if using mpi than ids must be offset
#ifdef USEMPI
        pdata[i].haloid+=haloidoffset;
#endif
        if (parentgid[i]!=GROUPNOPARENT) {
            pdata[i].hostid=opt.snapshotvalue+uparentgid[i];
            pdata[i].directhostid=opt.snapshotvalue+parentgid[i];
            if(opt.iKeepFOF && uparentgid[i]<=opt.num3dfof) pdata[i].hostfofid=opt.snapshotvalue+uparentgid[i];
            if(opt.iKeepFOF && uparentgid[i]>opt.num3dfof) pdata[i].hostfofid=GROUPNOPARENT;
#ifdef USEMPI
            pdata[i].hostid+=haloidoffset;
            pdata[i].directhostid+=haloidoffset;
            if(opt.iKeepFOF && uparentgid[i]<=opt.num3dfof) pdata[i].hostfofid+=haloidoffset;
#endif
        }
        else {
            pdata[i].hostid=GROUPNOPARENT;
            pdata[i].directhostid=GROUPNOPARENT;
            pdata[i].hostfofid=GROUPNOPARENT;
        }
        pdata[i].stype=stype[i];
        pdata[i].numsubs=nsub[i];
    }
#ifdef USEOPENMP
}
#endif
}
//@}

/// \name Routines used for interative search
//@{
///for halo mergers check, get mean velocity dispersion
inline Double_t GetVelDisp(Particle *Part, Int_t numgroups, Int_t *numingroup, Int_t **pglist){
    Int_t i;
    Double_t mt,sigmav=0;
    Coordinate cmval;
    Matrix dispmatrix;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,mt,cmval,dispmatrix)
{
#pragma omp for reduction(+:sigmav)
#endif
    for (i=1;i<=numgroups;i++) {
        for (int k=0;k<3;k++) cmval[k]=0.;
        mt=0.;
        for (Int_t j=0;j<numingroup[i];j++) {
            for (int k=0;k<3;k++) cmval[k]+=Part[pglist[i][j]].GetVelocity(k)*Part[pglist[i][j]].GetMass();
            mt+=Part[pglist[i][j]].GetMass();
        }
        mt=1.0/mt;
        for (int k=0;k<3;k++) cmval[k]*=mt;
        for (int k=0;k<3;k++) for (int l=0;l<3;l++) dispmatrix(k,l)=0.;
        for (Int_t j=0;j<numingroup[i];j++) {
            for (int k=0;k<3;k++) for (int l=k;l<3;l++) dispmatrix(k,l)+=Part[pglist[i][j]].GetMass()*(Part[pglist[i][j]].GetVelocity(k)-cmval[k])*(Part[pglist[i][j]].GetVelocity(l)-cmval[l]);
        }
        for (int k=0;k<3;k++) for (int l=k;l<3;l++) dispmatrix(k,l)*=mt;
        for (int k=1;k<3;k++) for (int l=0;l<k;l++) dispmatrix(k,l)=dispmatrix(l,k);
        sigmav+=pow(dispmatrix.Det(),1./3.);
    }
#ifdef USEOPENMP
}
#endif
    sigmav/=(double)numgroups;
    return sigmav;
}

///Create array for each group listing all previously unlinked particles that should now be linked. Note that the triple pointer newIndex is a double pointer, extra pointer layer is to ensure that
///it is passed by reference so that the memory allocated in the function is associated with the newIndex pointer once the function as returned.
inline void DetermineNewLinks(Int_t nsubset, Particle *Partsubset, Int_t *pfof, Int_t numgroups, Int_t &newlinks, Int_t *newlinksIndex, Int_t *numgrouplinksIndex, Int_t *nnID, Int_t ***newIndex) {
    Int_t pid,ppid;
    newlinks=0;
    for (Int_t j=1;j<=numgroups;j++) numgrouplinksIndex[j]=0;
    for (Int_t j=0;j<nsubset;j++) {
        ppid=Partsubset[j].GetID();
        if (nnID[ppid]!=pfof[ppid]) {
            if (nnID[ppid]>0&&pfof[ppid]==0) {newlinksIndex[newlinks++]=j;numgrouplinksIndex[nnID[ppid]]++;}
        }
    }
    for (Int_t j=1;j<=numgroups;j++) {
        (*newIndex)[j]=new Int_t[numgrouplinksIndex[j]+1];
        numgrouplinksIndex[j]=0;
    }
    for (Int_t j=0;j<nsubset;j++) {
        ppid=Partsubset[j].GetID();
        if (nnID[ppid]>0&&pfof[ppid]==0) {(*newIndex)[nnID[ppid]][numgrouplinksIndex[nnID[ppid]]++]=j;}
    }
}

///Create array for each group listing all possibly intergroup links
inline void DetermineGroupLinks(Int_t nsubset, Particle *Partsubset, Int_t *pfof, Int_t numgroups, Int_t &newlinks, Int_t *newlinksIndex, Int_t *numgrouplinksIndex, Int_t *nnID, Int_t ***newIndex) {
    Int_t i,pid,ppid;
    //store number of links between groups in newIndex
    for (i=1;i<=numgroups;i++) numgrouplinksIndex[i]=0;
    //first determine number of links a group has with particles identified as belonging to another group
    for (i=0;i<newlinks;i++) {
        ppid=Partsubset[newlinksIndex[i]].GetID();
        if (nnID[ppid]!=pfof[ppid]&&nnID[ppid]>-1) {numgrouplinksIndex[nnID[ppid]]++;}
    }
    //then go through list and store all these particles
    for (i=1;i<=numgroups;i++) {
        (*newIndex)[i]=new Int_t[numgrouplinksIndex[i]+1];
        numgrouplinksIndex[i]=0;
    }
    for (i=0;i<newlinks;i++) {
        ppid=Partsubset[newlinksIndex[i]].GetID();
        if (nnID[ppid]!=pfof[ppid]&&nnID[ppid]>-1) {(*newIndex)[nnID[ppid]][numgrouplinksIndex[nnID[ppid]]++]=newlinksIndex[i];}
    }
}
///once Group links are found, reduce the data so that one can determine for each group, number of other groups linked, their gids and the number of that group linked
///these are stored in numgrouplinksIndex, intergroupgidIndex, & newintergroupIndex
inline void DetermineGroupMergerConnections(Particle *Partsubset, Int_t numgroups, Int_t *pfof, int *ilflag, Int_t *numgrouplinksIndex, Int_t *intergrouplinksIndex, Int_t *nnID, Int_t ***newIndex, Int_t ***newintergroupIndex, Int_t ***intergroupgidIndex) {
    Int_t i,ii,pid,ppid, intergrouplinks;
    //ilflag ensures that group only assocaited once to current group being searched for intergroup links
    for (i=1;i<=numgroups;i++) ilflag[i]=0;
    for (i=1;i<=numgroups;i++) {
        intergrouplinks=0;
        for (ii=0;ii<numgrouplinksIndex[i];ii++) {
            ppid=Partsubset[(*newIndex)[i][ii]].GetID();
            if (ilflag[pfof[ppid]]==0) {intergrouplinksIndex[intergrouplinks++]=pfof[ppid];ilflag[pfof[ppid]]=1;}
        }
        if (intergrouplinks>0) {
            (*newintergroupIndex)[i]=new Int_t[intergrouplinks];
            (*intergroupgidIndex)[i]=new Int_t[intergrouplinks];
            for (ii=0;ii<intergrouplinks;ii++) {
                (*intergroupgidIndex)[i][ii]=intergrouplinksIndex[ii];
                (*newintergroupIndex)[i][ii]=0;
                ilflag[intergrouplinksIndex[ii]]=0;
            }
            for (ii=0;ii<numgrouplinksIndex[i];ii++) {
                ppid=Partsubset[(*newIndex)[i][ii]].GetID();
                for (Int_t j=0;j<intergrouplinks;j++) if (pfof[ppid]==(*intergroupgidIndex)[i][j]) {(*newintergroupIndex)[i][j]++;}
            }
        }
        numgrouplinksIndex[i]=intergrouplinks;
        delete[] (*newIndex)[i];
    }
}

///Routine uses a tree to mark all particles meeting the search criteria given by the FOF comparison function and the paramaters in param. The particles are marked in nnID in ID order
///Only the set of particles in newlinksIndex are searched for other particles meeting this criteria. Note that the tree search is set up such that if the marker value (here given by pfof) greater than the particle's
///current id marker or the particle's current marker is 0, the particle's nnID value is set to the marker value.
inline void SearchForNewLinks(Int_t nsubset, KDTree *tree, Particle *Partsubset, Int_t *pfof, FOFcompfunc &fofcmp, Double_t *param, Int_t newlinks, Int_t *newlinksIndex, Int_t **nnID, Double_t **dist2, int nthreads) {
    Int_t i,ii;
    int tid;
#ifdef USEOPENMP
        if (newlinks>ompsearchnum) {
        for (ii=0;ii<nsubset;ii++) for (int j=1;j<nthreads;j++) nnID[0][ii+j*nsubset]=nnID[0][ii];
#pragma omp parallel default(shared) \
private(ii,tid)
{
#pragma omp for schedule(dynamic,1) nowait
        for (ii=0;ii<newlinks;ii++) {
            tid=omp_get_thread_num();
            tree->SearchCriterion(newlinksIndex[ii], fofcmp,param, pfof[Partsubset[newlinksIndex[ii]].GetID()], &nnID[0][tid*nsubset], &dist2[0][tid*nsubset]);
        }
}
#pragma omp parallel default(shared) \
private(ii,tid)
{
#pragma omp for
        for (ii=0;ii<nsubset;ii++) {
            for (int j=1;j<nthreads;j++) {
                if (nnID[0][ii]==0) nnID[0][ii]=nnID[0][ii+j*nsubset];
                else if (nnID[0][ii+j*nsubset]>0) {
                    if (nnID[0][ii+j*nsubset]<nnID[0][ii]) nnID[0][ii]=nnID[0][ii+j*nsubset];
                }
            }
        }
}
        }
        else
#endif
        {
            tid=0;
            for (ii=0;ii<newlinks;ii++) {
                tree->SearchCriterion(newlinksIndex[ii], fofcmp,param, pfof[Partsubset[newlinksIndex[ii]].GetID()], nnID[tid], dist2[tid]);
            }
        }
}
inline void SearchForNewLinks(Int_t nsubset, KDTree *tree, Particle *Partsubset, Int_t *pfof, FOFcompfunc &fofcmp, Double_t *param, Int_t newlinks, Int_t *newlinksIndex, Int_t **nnID, int nthreads) {
    Int_t i,ii;
    int tid;
#ifdef USEOPENMP
        if (newlinks>ompsearchnum) {
        for (ii=0;ii<nsubset;ii++) for (int j=1;j<nthreads;j++) nnID[0][ii+j*nsubset]=nnID[0][ii];
#pragma omp parallel default(shared) \
private(ii,tid)
{
#pragma omp for schedule(dynamic,1) nowait
        for (ii=0;ii<newlinks;ii++) {
            tid=omp_get_thread_num();
            tree->SearchCriterion(newlinksIndex[ii], fofcmp,param, pfof[Partsubset[newlinksIndex[ii]].GetID()], &nnID[0][tid*nsubset]);
        }
}
#pragma omp parallel default(shared) \
private(ii,tid)
{
#pragma omp for
        for (ii=0;ii<nsubset;ii++) {
            for (int j=1;j<nthreads;j++) {
                if (nnID[0][ii]==0) nnID[0][ii]=nnID[0][ii+j*nsubset];
                else if (nnID[0][ii+j*nsubset]>0) {
                    if (nnID[0][ii+j*nsubset]<nnID[0][ii]) nnID[0][ii]=nnID[0][ii+j*nsubset];
                }
            }
        }
}
        }
        else
#endif
        {
            tid=0;
            for (ii=0;ii<newlinks;ii++) {
                tree->SearchCriterion(newlinksIndex[ii], fofcmp,param, pfof[Partsubset[newlinksIndex[ii]].GetID()], nnID[tid]);
            }
        }
}

///This routine links particle's stored in newIndex to a group. Routine assumes that each group's list of new particle's is independent.
inline void LinkUntagged(Particle *Partsubset, Int_t numgroups, Int_t *pfof, Int_t *numingroup, Int_t **pglist, Int_t newlinks, Int_t *numgrouplinksIndex, Int_t **newIndex, Int_tree_t *Head, Int_tree_t *Next, Int_tree_t *GroupTail, Int_t *nnID) {
    Int_t i,pindex,tail,pid,ppid;
    if (newlinks>0) {
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,pindex,tail,pid,ppid)
{
#pragma omp for schedule(dynamic,1) nowait
#endif
        for (i=1;i<=numgroups;i++) if (numgrouplinksIndex[i]>0) {
            pindex=pglist[i][0];
            tail=GroupTail[i];
            for (Int_t k=0;k<numgrouplinksIndex[i];k++) {
                pid=newIndex[i][k];
                ppid=Partsubset[pid].GetID();
                Head[pid]=Head[pindex];
                Next[tail]=pid;
                tail=pid;
                pfof[ppid]=i;
                nnID[ppid]=i;
            }
            numingroup[i]+=numgrouplinksIndex[i];
            GroupTail[i]=newIndex[i][numgrouplinksIndex[i]-1];
        }
#ifdef USEOPENMP
}
#endif
    }
}

///This routine merges candidate substructre groups together so long as the groups share enough links compared to the groups original size given in oldnumingroup
///(which is generally the group size prior to the expanded search).
inline Int_t MergeGroups(Options &opt, Particle *Partsubset, Int_t numgroups, Int_t *pfof, Int_t *numingroup, Int_t *oldnumingroup, Int_t **pglist, Int_t *numgrouplinksIndex, Int_t ***intergroupgidIndex, Int_t ***newintergroupIndex, Int_t *intergrouplinksIndex, Int_tree_t *Head, Int_tree_t *Next, Int_tree_t *GroupTail, int *igflag, Int_t *nnID, Int_t &newlinks, Int_t *newlinksIndex) {
    Int_t i,gid,gidjoin,pindex,lindex,sindex,tail,ss,intergrouplinks,mergers=0;
    for (i=1;i<=numgroups;i++) if (igflag[i]==0&&numgrouplinksIndex[i]>0) {
        intergrouplinks=0;
        for (Int_t j=0;j<numgrouplinksIndex[i];j++) {
            gidjoin=(*intergroupgidIndex)[i][j];
            if (igflag[gidjoin]==0&&i!=gidjoin) {
            //merge if enough links relative to original size of group
            int merge=((*newintergroupIndex)[i][j]>opt.fmerge*oldnumingroup[gidjoin]);
            if (merge) {
                mergers++;
                intergrouplinksIndex[gidjoin]=0;numingroup[gidjoin]=0;
                igflag[gidjoin]=1;
                pindex=pglist[i][0];
                sindex=pglist[gidjoin][0];lindex=pindex;
                gid=pfof[Partsubset[Head[lindex]].GetID()];
                //first adjust length of the larger group
                tail=GroupTail[i];
                //then adjust the next of the tail of the larger group
                Next[tail]=Head[sindex];
                ss=Head[sindex];
                //then adjust the group id, head and group length of the smaller group
                do {
                    pfof[Partsubset[ss].GetID()]=gid;
                    nnID[Partsubset[ss].GetID()]=gid;
                    Head[ss]=Head[lindex];
                    intergrouplinks++;
                    newlinksIndex[newlinks++]=ss;
                } while((ss = Next[ss]) >= 0);
                GroupTail[i]=GroupTail[gidjoin];
                if (numgrouplinksIndex[gidjoin]>0) {delete[] (*newintergroupIndex)[gidjoin];delete[] (*intergroupgidIndex)[gidjoin];}
            }
            }
        }
        numingroup[i]+=intergrouplinks;
        delete[] (*newintergroupIndex)[i];delete[] (*intergroupgidIndex)[i];
    }
    return mergers;
}

///This routine merges candidate halo groups together during the merger check so long as the groups share enough links compared to the groups original size given in oldnumingroup
///or one group is significantly larger than the other (in which case the secondary object is probably a substructure of the other)
inline Int_t MergeHaloGroups(Options &opt, Particle *Partsubset, Int_t numgroups, Int_t *pfof, Int_t *numingroup, Int_t *oldnumingroup, Int_t **pglist, Int_t *numgrouplinksIndex, Int_t ***intergroupgidIndex, Int_t ***newintergroupIndex, Int_t *intergrouplinksIndex, Int_t *Head, Int_t *Next, Int_t *GroupTail, int *igflag, Int_t *nnID, Int_t &newlinks, Int_t *newlinksIndex) {
    Int_t i,gid,gidjoin,pindex,lindex,sindex,tail,ss,intergrouplinks,mergers=0;
    for (i=1;i<=numgroups;i++) if (igflag[i]==0&&numgrouplinksIndex[i]>0) {
        intergrouplinks=0;
        for (Int_t j=0;j<numgrouplinksIndex[i];j++) {
            gidjoin=(*intergroupgidIndex)[i][j];
            if (igflag[gidjoin]==0&&i!=gidjoin) {
            //merge if enough links relative to original size of group
            int merge=(((*newintergroupIndex)[i][j]>opt.fmergebg*oldnumingroup[gidjoin])||((Double_t)oldnumingroup[gidjoin]/(Double_t)oldnumingroup[i]<opt.HaloMergerRatio*opt.fmergebg));
            if (merge) {
                mergers++;
                intergrouplinksIndex[gidjoin]=0;numingroup[gidjoin]=0;
                igflag[gidjoin]=1;
                pindex=pglist[i][0];
                sindex=pglist[gidjoin][0];lindex=pindex;
                gid=pfof[Partsubset[Head[lindex]].GetID()];
                //first adjust length of the larger group
                tail=GroupTail[i];
                //then adjust the next of the tail of the larger group
                Next[tail]=Head[sindex];
                ss=Head[sindex];
                //then adjust the group id, head and group length of the smaller group
                do {
                    pfof[Partsubset[ss].GetID()]=gid;
                    nnID[Partsubset[ss].GetID()]=gid;
                    Head[ss]=Head[lindex];
                    intergrouplinks++;
                    newlinksIndex[newlinks++]=ss;
                } while((ss = Next[ss]) >= 0);
                GroupTail[i]=GroupTail[gidjoin];
                if (numgrouplinksIndex[gidjoin]>0) {delete[] (*newintergroupIndex)[gidjoin];delete[] (*intergroupgidIndex)[gidjoin];}
            }
            }
        }
        numingroup[i]+=intergrouplinks;
        delete[] (*newintergroupIndex)[i];delete[] (*intergroupgidIndex)[i];
    }
    return mergers;
}

//@}
