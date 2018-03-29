/*! \file swiftinterface.cxx
 *  \brief this file contains routines that allow the velociraptor library to interface with the swift N-body code from within swift.
 */


#include "swiftinterface.h"

#ifdef SWIFTINTERFACE

Options libvelociraptorOpt;
//KDTree *mpimeshtree;
//Particle *mpimeshinfo;

void InitVelociraptor(char* configname, char* outputname, cosmoinfo c, unitinfo u, siminfo s)
{

#ifdef USEMPI
    //find out how big the SPMD world is
    MPI_Comm_size(MPI_COMM_WORLD,&NProcs);
    //mpi_domain=new MPI_Domain[NProcs];
    mpi_nlocal=new Int_t[NProcs];
    mpi_nsend=new Int_t[NProcs*NProcs];
    mpi_ngroups=new Int_t[NProcs];
    //and this processes' rank is
    MPI_Comm_rank(MPI_COMM_WORLD,&ThisTask);
    //store MinSize as when using mpi prior to stitching use min of 2;
    MinNumMPI=2;

#endif
    cout<<"Initialising VELOCIraptor..."<< endl;

    libvelociraptorOpt.pname = configname;
    libvelociraptorOpt.outname = outputname;

    cout<<"Reading VELOCIraptor config file..."<< endl;
    GetParamFile(libvelociraptorOpt);
    cout<<"Setting cosmology, units, sim stuff "<<endl;
    ///set units, here idea is to convert internal units so that have kpc, km/s, solar mass
    libvelociraptorOpt.lengthtokpc=1.0;
    libvelociraptorOpt.velocitytokms=1.0;
    libvelociraptorOpt.masstosolarmass=1.0;
    libvelociraptorOpt.L=u.lengthtokpc;
    libvelociraptorOpt.M=u.masstosolarmass;
    libvelociraptorOpt.V=u.velocitytokms;
    ///these should be in units of kpc, km/s, and solar mass
    libvelociraptorOpt.G=u.gravity;
    libvelociraptorOpt.H=u.hubbleunit;
    ///set cosmology
    libvelociraptorOpt.a=c.atime;
    libvelociraptorOpt.h=c.littleh;
    libvelociraptorOpt.Omega_m=c.Omega_m;
    libvelociraptorOpt.Omega_b=c.Omega_b;
    libvelociraptorOpt.Omega_cdm=c.Omega_cdm;
    libvelociraptorOpt.Omega_Lambda=c.Omega_Lambda;
    libvelociraptorOpt.w_de=c.w_de;

    //if libvelociraptorOpt.virlevel<0, then use virial overdensity based on Bryan and Norman 1998 virialization level is given by
    if (libvelociraptorOpt.virlevel<0)
    {
        Double_t bnx=-((1-libvelociraptorOpt.Omega_m-libvelociraptorOpt.Omega_Lambda)*pow(libvelociraptorOpt.a,-2.0)+libvelociraptorOpt.Omega_Lambda)/((1-libvelociraptorOpt.Omega_m-libvelociraptorOpt.Omega_Lambda)*pow(libvelociraptorOpt.a,-2.0)+libvelociraptorOpt.Omega_m*pow(libvelociraptorOpt.a,-3.0)+libvelociraptorOpt.Omega_Lambda);
        libvelociraptorOpt.virlevel=(18.0*M_PI*M_PI+82.0*bnx-39*bnx*bnx)/libvelociraptorOpt.Omega_m;
    }
    //set some sim information
    libvelociraptorOpt.p=s.period;
    libvelociraptorOpt.zoomlowmassdm=s.zoomhigresolutionmass;
    libvelociraptorOpt.icosmologicalin=s.icosmologicalsim;
    libvelociraptorOpt.ellxscale=s.interparticlespacing;
    libvelociraptorOpt.uinfo.eps*=libvelociraptorOpt.ellxscale;
    if (libvelociraptorOpt.icosmologicalin) {
        libvelociraptorOpt.ellxscale*=libvelociraptorOpt.a;
        libvelociraptorOpt.uinfo.eps*=libvelociraptorOpt.a;
        double Hubble=libvelociraptorOpt.h*libvelociraptorOpt.H*sqrt((1-libvelociraptorOpt.Omega_m-libvelociraptorOpt.Omega_Lambda)*pow(libvelociraptorOpt.a,-2.0)+libvelociraptorOpt.Omega_m*pow(libvelociraptorOpt.a,-3.0)+libvelociraptorOpt.Omega_Lambda);
        libvelociraptorOpt.rhobg=3.*Hubble*Hubble/(8.0*M_PI*libvelociraptorOpt.G)*libvelociraptorOpt.Omega_m;
    }
    else {
        libvelociraptorOpt.rhobg=1.0;
    }
    //assume above is in comoving is cosmological simulation and then need to correct to physical
    if (libvelociraptorOpt.icosmologicalin) {
        libvelociraptorOpt.p*=libvelociraptorOpt.a;
        libvelociraptorOpt.ellxscale*=libvelociraptorOpt.a;
        libvelociraptorOpt.uinfo.eps*=libvelociraptorOpt.a;
    }
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

    cout<<"Finished initialising VELOCIraptor"<<endl;
    if (libvelociraptorOpt.HaloMinSize==-1) libvelociraptorOpt.HaloMinSize=libvelociraptorOpt.MinSize;

#ifdef USEMPI
    //if single halo, use minsize to initialize the old minimum number
    //else use the halominsize since if mpi and not single halo, halos localized to mpi domain for substructure search
    if (libvelociraptorOpt.iSingleHalo) MinNumOld=libvelociraptorOpt.MinSize;
    else MinNumOld=libvelociraptorOpt.HaloMinSize;
    mpi_period = libvelociraptorOpt.p;
#endif
    //allocate memory to store metis mesh info
    //mpimeshinfo=new Particle[libvelociraptorOpt.numcells];

    //write velociraptor info
    WriteVELOCIraptorConfig(libvelociraptorOpt);
    WriteSimulationInfo(libvelociraptorOpt);
    WriteUnitInfo(libvelociraptorOpt);

}

void InvokeVelociraptor(const int num_gravity_parts, struct gpart *gravity_parts, const int *cell_node_ids) {
#ifndef USEMPI
    int ThisTask=0;
    int NProcs=1;
    Int_t Nlocal, Ntotal, Nmemlocal;
#endif
    int nthreads;
#ifdef USEOPENMP
#pragma omp parallel
{
    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
}
#else
    nthreads=1;
#endif

    Particle *parts;
    Int_t *pfof,*numingroup,**pglist;
    Int_t ngroup, nhalos;
    //Int_t nbodies,nbaryons,ndark;
    //KDTree *tree;
    //to store information about the group
    PropData *pdata=NULL,*pdatahalos=NULL;
    double time1;

    /// Set pointer to cell node IDs
    libvelociraptorOpt.cellnodeids = cell_node_ids;

    Nlocal=Nmemlocal=num_gravity_parts;
    Nmemlocal*=(1+libvelociraptorOpt.mpipartfac); /* JSW: Not set in parameter file. */
#ifdef USEMPI
    MPI_Allreduce(&Nlocal, &Ntotal, 1, MPI_Int_t, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allgather(&Nlocal, 1, MPI_Int_t, mpi_nlocal, 1, MPI_Int_t, MPI_COMM_WORLD);
#else
    Ntotal=Nlocal;
#endif
    parts=new Particle[Nmemlocal];
    cout<<"Copying particle data..."<< endl;
    time1=MyGetTime();
    for(auto i=0; i<Nlocal; i++) {
        parts[i] = Particle(gravity_parts[i], libvelociraptorOpt.L, libvelociraptorOpt.V, libvelociraptorOpt.M, libvelociraptorOpt.icosmologicalin,libvelociraptorOpt.a,libvelociraptorOpt.h);
        parts[i].SetType(DARKTYPE);
    }
    time1=MyGetTime()-time1;
    cout<<"Finished copying particle data."<< endl;
    Nlocal=num_gravity_parts; /* JSW: Already set previously to same value. */
    cout<<"TIME::"<<ThisTask<<" took "<<time1<<" to copy "<<Nlocal<<" particles from SWIFT to a local format. Out of "<<Ntotal<<endl;
    cout<<ThisTask<<" There are "<<Nlocal<<" particles and have allocated enough memory for "<<Nmemlocal<<" requiring "<<Nmemlocal*sizeof(Particle)/1024./1024./1024.<<"GB of memory "<<endl;
    //if (libvelociraptorOpt.iBaryonSearch>0) cout<<ThisTask<<"There are "<<Nlocalbaryon[0]<<" baryon particles and have allocated enough memory for "<<Nmemlocalbaryon<<" requiring "<<Nmemlocalbaryon*sizeof(Particle)/1024./1024./1024.<<"GB of memory "<<endl;
    cout<<ThisTask<<" will also require additional memory for FOF algorithms and substructure search. Largest mem needed for preliminary FOF search. Rough estimate is "<<Nlocal*(sizeof(Int_tree_t)*8)/1024./1024./1024.<<"GB of memory"<<endl;

    //build binary tree to store cells for quick search
    /*
    Double_t *period=NULL;
    if (libvelociraptorOpt.p>0) {
        period=new Double_t[3];
        for (int j=0;j<3;j++) period[j]=libvelociraptorOpt.p;
    }
    for (auto i=0;i<libvelociraptorOpt.numcells; i++) for (auto j=0;j<3;j++) mpimeshinfo[i].SetPosition(j,libvelociraptorOpt.cellloc[i].loc[j]);
    mpimeshtree=new KDTree(mpimeshinfo,libvelociraptorOpt.numcells,1,mpimeshtree->TPHYS,mpimeshtree->KEPAN,100,0,0,0,period);
    */

    //
    // Perform FOF search.
    //
    time1=MyGetTime();
    pfof=SearchFullSet(libvelociraptorOpt,Nlocal,parts,ngroup);
    time1=MyGetTime()-time1;
    cout<<"TIME::"<<ThisTask<<" took "<<time1<<" to search "<<Nlocal<<" with "<<nthreads<<endl;
    nhalos=ngroup;
    //if caculating inclusive halo masses, then for simplicity, I assume halo id order NOT rearranged!
    //this is not necessarily true if baryons are searched for separately.
    if (libvelociraptorOpt.iInclusiveHalo) {
        pdatahalos=new PropData[nhalos+1];
        Int_t *numinhalos=BuildNumInGroup(Nlocal, nhalos, pfof);
        Int_t *sortvalhalos=new Int_t[Nlocal];
        Int_t *originalID=new Int_t[Nlocal];
        for (Int_t i=0;i<Nlocal;i++) {sortvalhalos[i]=pfof[i]*(pfof[i]>0)+Nlocal*(pfof[i]==0);originalID[i]=parts[i].GetID();parts[i].SetID(i);}
        Int_t *noffsethalos=BuildNoffset(Nlocal, parts, nhalos, numinhalos, sortvalhalos);
        GetInclusiveMasses(libvelociraptorOpt, Nlocal, parts, nhalos, pfof, numinhalos, pdatahalos, noffsethalos);
        qsort(parts,Nlocal,sizeof(Particle),IDCompare);
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
        cout<<"Searching subset"<<endl;
        time1=MyGetTime();
        //if groups have been found (and localized to single MPI thread) then proceed to search for subsubstructures
        SearchSubSub(libvelociraptorOpt, Nlocal, parts, pfof,ngroup,nhalos,pdatahalos);
        time1=MyGetTime()-time1;
        cout<<"TIME::"<<ThisTask<<" took "<<time1<<" to search for substructures "<<Nlocal<<" with "<<nthreads<<endl;
    }
    pdata=new PropData[ngroup+1];
    //if inclusive halo mass required
    if (libvelociraptorOpt.iInclusiveHalo && ngroup>0) {
        CopyMasses(nhalos,pdatahalos,pdata);
        delete[] pdatahalos;
    }
    //need to add baryon interface

    //get mpi local hierarchy
    Int_t *nsub,*parentgid, *uparentgid,*stype;
    nsub=new Int_t[ngroup+1];
    parentgid=new Int_t[ngroup+1];
    uparentgid=new Int_t[ngroup+1];
    stype=new Int_t[ngroup+1];
    Int_t nhierarchy=GetHierarchy(libvelociraptorOpt,ngroup,nsub,parentgid,uparentgid,stype);
    CopyHierarchy(libvelociraptorOpt,pdata,ngroup,nsub,parentgid,uparentgid,stype);

    //calculate data and output
    numingroup=BuildNumInGroup(Nlocal, ngroup, pfof);
    pglist=SortAccordingtoBindingEnergy(libvelociraptorOpt,Nlocal,parts,ngroup,pfof,numingroup,pdata);//alters pglist so most bound particles first
    WriteProperties(libvelociraptorOpt,ngroup,pdata);
    WriteGroupCatalog(libvelociraptorOpt, ngroup, numingroup, pglist, parts);
    for (Int_t i=1;i<=ngroup;i++) delete[] pglist[i];
    delete[] pglist;
    delete[] parts;

    //delete tree used to search mpi mesh
    //delete mpimeshtree;

    ///\todo need to return fof and substructure information back to swift

    cout<<"VELOCIraptor returning."<< endl;

}

#endif
