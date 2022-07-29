/*! \file generate.cxx
 *  \brief this file contains routines to generate a test input file
 */


#include "logging.h"
#include "stf.h"
#include "timer.h"

#include <random>
#ifdef USEHDF
#include "hdfitems.h"
#endif

/// generate input file
void GenerateInput(Options &opt, vector<Particle> &Part) {
#ifndef USEMPI
    int ThisTask = 0, NProcs = 1;
    Int_t Nlocal;
#endif
    if (opt.iGenerateInput == false) return;
    vr::Timer total_timer;
    vector<GaussianDistrib> Gaus;
    LOG_RANK0(info)<<" Generating Input ... ";
    opt.Ngenerate = pow(opt.Ngenerate, 3.0);
    //set some cosmology
    opt.h = 1;
    opt.p = 100.0;
    opt.Omega_m = 1.0;
    opt.Omega_Lambda = 0.0;
    opt.Omega_b = 0.0;
    opt.Omega_cdm = opt.Omega_m;
    opt.Omega_k = 0;
    opt.Omega_r = 0.0;
    opt.Omega_nu = 0.0;
    opt.Omega_de = 0.0;
    opt.lengthtokpc=1.0;
    opt.velocitytokms=1.0;
    opt.masstosolarmass=1.0;
    if (opt.fbackground > 1) {
        opt.fbackground = 1;
        opt.Ngeneratehalos = 0;
    }
    if (opt.fbackground < 0) {
        opt.fbackground = 0;
        opt.Ngeneratehalos = opt.Ngenerate/1000;
    }

    opt.G = CalcGravitationalConstant(opt);
    opt.H = CalcHubbleUnit(opt);
    CalcOmegak(opt);
    CalcCriticalDensity(opt, opt.a);
    CalcBackgroundDensity(opt, opt.a);
    CalcVirBN98(opt, opt.a);
    opt.mpgenerate = opt.rhobg * pow(opt.p, 3.0)/opt.Ngenerate;

    Nlocal = opt.Ngenerate/NProcs;
    opt.Ngeneratehalos = opt.Ngeneratehalos/static_cast<double>(NProcs);
    if (NProcs > 1) {
        if (ThisTask == NProcs - 1)
        {
            Nlocal = opt.Ngenerate-Nlocal*(NProcs-1);
        }
    }
    opt.Nbackground = opt.fbackground * opt.Ngenerate;
    LOG_RANK0(info)<<" producing positions for "<<opt.Ngenerate<<" particles, echo of size "<<sizeof(Particle)<<" requiring "<<opt.Ngenerate*sizeof(Particle)/1024./1024./1024.<<" GB of memory";
    LOG_RANK0(info)<<" Number of iterations "<<opt.generate_Niter<<" and produce output "<<opt.igenerate_output;
    for (auto iter=0;iter<opt.generate_Niter;iter++) {
    LOG_RANK0(info)<<" Iteration "<<iter;
    Part.resize(Nlocal);
    for (auto p:Part) p.SetMass(opt.mpgenerate);
    LOG(info) << GetMemUsage(__func__+string("--line--")+to_string(__LINE__));

    Gaus = ProduceGaussians(opt, Part);
    Int_t npoints = 0;
    for (auto &x:Gaus) npoints += x.npoints;
    //update the background number of points
    opt.Nbackground = Nlocal - npoints;
    LOG(info)<<" producing positions for "<<opt.Ngeneratehalos<<" halos containing total of "<<npoints<<" particles and a background of "<<opt.Nbackground;
    opt.fbackground = opt.Nbackground/static_cast<double>(Nlocal);
    PopulateGaussians(opt, Part, Gaus);
    ProduceBackground(opt, Part, npoints);
    #ifdef USEMPI
    opt.numcellsperdim=320;
    for (auto i=0;i<3;i++) {
        mpi_xlim[i][1] = opt.p;  
        mpi_xlim[i][0] = 0;
    }
    MPIInitialDomainDecompositionWithMesh(opt);
    PlaceGeneratedInMesh(opt, Part);
    #endif
    if (iter != opt.generate_Niter -1) {
        Part.clear();
        Part.shrink_to_fit();
        Gaus.clear();
        Gaus.shrink_to_fit();
    }
    }
    if (opt.igenerate_output) WriteGeneratedInput(opt, Part, Gaus);
    LOG(info)<<" Generating input: Done " << total_timer;
#ifdef USEMPI
    MPI_Finalize();
#endif
    exit(0);
}

/// produce a background uniform distribution of particles within the volume
void ProduceBackground(Options &opt, vector<Particle> &Part, Int_t noffset)
{
#ifndef USEMPI
    int ThisTask = 0, NProcs = 1;
    unsigned long long Nlocal = Part.size();
#endif
    vr::Timer local_timer;
    LOG(info)<<" Produce Background ... ";

    if (opt.fbackground == 0) return;
    Int_t nbg = Nlocal - noffset;
#if defined(USEOPENMP)
#pragma omp parallel default(shared)
{
#endif
    unsigned seed = 4320;
    #ifdef USEMPI 
    seed *= (ThisTask+1);
    #endif
    default_random_engine generator(seed);
    uniform_real_distribution<double> pos(0.0,opt.p);
    uniform_real_distribution<double> vel(-opt.H,opt.H);
    Coordinate x,v;
#if defined(USEOPENMP)
    #pragma omp for schedule(static)
#endif
    for (auto i = 0; i < nbg; i++) {
        for (auto j = 0; j < 3; j++) {
            x[j] = pos(generator);
            v[j] = vel(generator);
        }
        Part[i+noffset].SetPosition(x.GetCoord());
        Part[i+noffset].SetVelocity(v.GetCoord());
    }
#if defined(USEOPENMP)
}
#endif
    LOG(info)<<" Produce Background: Done " << local_timer;
    LOG(info)<<GetMemUsage(__func__+string("--line--")+to_string(__LINE__));
}

inline Matrix NormalizeDisp(double sX[3], double ssX[3], double norm, uniform_real_distribution<double> &disp, default_random_engine& generator){
    for (auto j = 0; j < 3; j++) {
        do {
            sX[j] = disp(generator);
        }while (sX[j]==0);
        ssX[j] = disp(generator) - 0.5;
    }
    Matrix X;
    for (auto j = 0; j < 3; j++) {
        for (auto k = j+1; k < 3; k++) {
            X(j,k) = X(k,j) = ssX[j]*sX[j]*sX[k];
        }
    }
    for (auto j = 0; j < 3; j++) X(j,j) = sX[j]*sX[j];
    norm /= pow(X.Det(), 1.0/6.0);
    X = X*norm*norm;
    return X;

}

/// produce Gaussian
vector<GaussianDistrib> ProduceGaussians(Options &opt, vector<Particle> &Part)
{
#ifndef USEMPI
    int ThisTask = 0, NProcs = 1;
    unsigned long long Nlocal = Part.size();
#endif
    vr::Timer local_timer;
    LOG(info)<<" Produce Gaussians ... ";

    vector<GaussianDistrib> Gaus;
    if (opt.fbackground == 1 || opt.Ngeneratehalos == 0) return Gaus;
    Gaus.resize(opt.Ngeneratehalos);
    auto Ninhalos = Nlocal*(1-opt.fbackground);
    Int_t maxnpart = max(static_cast<double>(Ninhalos)/(100*static_cast<double>(opt.Ngeneratehalos)),100.0);
    Int_t npoints = 0;
    unsigned seed=4320;
    #ifdef USEMPI 
    seed *= (ThisTask+1);
    #endif
    default_random_engine generator(seed);
    double mmin = log(20.0), mmax = log(maxnpart);
    uniform_real_distribution<double> npart(mmin, mmax);
    for (auto i = 0; i < opt.Ngeneratehalos; i++) 
    {
        if (npoints + exp(mmax)/2 > Ninhalos) {
            opt.Ngeneratehalos = i;
            Gaus.resize(opt.Ngeneratehalos);
            break;
        }
        do {
            Gaus[i].npoints = exp(npart(generator));
        } while (Gaus[i].npoints + npoints > Ninhalos);
        npoints += Gaus[i].npoints;
    }

#if defined(USEOPENMP)
#pragma omp parallel default(shared) firstprivate(generator)
{
#endif
    uniform_real_distribution<double> pos(0.0, opt.p);
    uniform_real_distribution<double> vel(-opt.H, opt.H);
    uniform_real_distribution<double> scaling(0.1, 0.5);
    uniform_real_distribution<double> disp(0, 1.0);
    Coordinate x,v;
    double mass, sigX, sigV, sX[3], ssX[3], sV[3], ssV[3];
    Matrix X, V;
#if defined(USEOPENMP)
    #pragma omp for schedule(static) reduction(+:npoints)
#endif
    for (auto i = 0; i < opt.Ngeneratehalos; i++) {
        for (auto j = 0; j < 3; j++) {
            Gaus[i].mean[j] = pos(generator);
            Gaus[i].mean[j+3] = vel(generator);
        }
        for (auto j = 0; j < 36; j++) Gaus[i].covar[j] = 0;
        mass = Gaus[i].npoints*opt.mpgenerate;
        sigX = pow(mass/(200.0*4.0*M_PI/3.0*opt.rhocrit),1.0/3.0)*scaling(generator);
        sigV = sqrt(opt.G*mass/sigX);
        X = NormalizeDisp(sX, ssX, sigX, disp, generator);
        V = NormalizeDisp(sV, ssV, sigV, disp, generator);
        Gaus[i].mass = mass;
        Gaus[i].sigX = sigX;
        Gaus[i].sigV = sigV;
        for (auto j = 0; j < 3; j++) {
            for (auto k = 0; k <3; k++) {
                Gaus[i].covar[j*6+k] = X(j,k);
                Gaus[i].covar[(j+3)*6+k+3] = V(j,k);
            }
        }
    }
#if defined(USEOPENMP)
}
#endif
    LOG(info)<<" Produce Gaussians: Done "<<local_timer;
    LOG(info)<<GetMemUsage(__func__+string("--line--")+to_string(__LINE__));
    return Gaus;

}

/// populate gauassians
void PopulateGaussians(Options &opt, vector<Particle> &Part, vector<GaussianDistrib> &Gaus)
{
    if (Gaus.size() == 0) return;
#ifndef USEMPI
    int ThisTask = 0, NProcs = 1;
#endif
    int nthreads = 1;
    vr::Timer local_timer;
    LOG(info)<<" Populate Gaussians ... ";
    if (opt.Ngeneratehalos == 0) return;

    // get random numbers first then transform them as necessary
    Int_t npoints = 0;
    // vector<Int_t> noffset(Gaus.size(),0);
    Int_t *noffset = new Int_t[Gaus.size()];
    for (auto i = 0; i < opt.Ngeneratehalos; i++)
    {
        noffset[i] = npoints;
        npoints += Gaus[i].npoints;
    }
    LOG(info)<<" Gaussians contain "<<npoints;
    float *rn = new float[npoints*6];
    LOG(info)<<GetMemUsage(__func__+string("--line--")+to_string(__LINE__));
#if defined(USEOPENMP)
#pragma omp parallel default(shared)
{
#endif
    unsigned seed = 4320;
    #ifdef USEMPI 
    seed *= (ThisTask+1);
    #endif
    default_random_engine generator(seed);
    normal_distribution<double> g(0.0, 1.0);
#if defined(USEOPENMP)
    #pragma omp for schedule(static)
#endif
    for (auto i = 0; i < npoints*6; i++) {
        rn[i] = g(generator);
    }
#if defined(USEOPENMP)
}
#endif

    // now transform points based on mean and variance
    //
#if defined(USEOPENACC)
{
//if calculating with GPUS, then use float arrays to calc positions
unsigned long long N = 1000000000;
float *xx = new float[N];
float *yy = new float[N];
unsigned long long ii;
for (ii=0;ii<N;ii++) xx[ii] = ii;
#pragma acc kernels copyin(xx[0:N]), copyout(yy[0:N])
for (ii = 0; ii < N; ii++)
{
    yy[ii] = xx[ii]*xx[ii];
}
delete[] xx;
delete[] yy;
}
#endif

//
#if defined(USEOPENMP) && defined(USEOPENMPGPU)
{
//if calculating with GPUS, then use float arrays to calc positions
unsigned long long N = 1000000000;
float *xx = new float[N];
float *yy = new float[N];
unsigned long long ii;
for (ii=0;ii<N;ii++) xx[ii] = ii;
#pragma omp target parallel for map(to: xx) map(from: yy)
for (ii = 0; ii < N; ii++)
{
    yy[ii] = xx[ii]*xx[ii];
}
delete[] xx;
delete[] yy;
}
#endif

// #if defined(USEOPENACC)
// #else

#if defined(USEOPENMP) && defined(USEOPENMPGPU)
    //with GPUS, need to store eigenvectors and eigenvalues of all the halos
    float *eigenvalarray, *eigenvecarray, *newphase;
    Int_t *eigoffset;
    eigenvalarray = new float[opt.Ngeneratehalos * 6];
    eigenvecarray = new float[opt.Ngeneratehalos * 6 * 6];
    newphase = new float[npoints * 6];
    eigoffset = new Int_t[npoints];
#endif
#ifdef USEOPENMP
    int chunksize = 10;
#pragma omp parallel if (opt.Ngeneratehalos>ompgeneratehalos && npoints > ompgeneratenum) default(shared)
{
#endif
    GMatrix y(1,6), x(6,1), s(6,6), eigenval(6,1), eigenvec(6,6);
#ifdef USEOPENMP
    #pragma omp for schedule(dynamic,chunksize) nowait
#endif
    for (auto i = 0; i < opt.Ngeneratehalos; i++)
    {
        //get the eigenvalues and eigenvectors of the covariance matrix
        //then take normal random numbers and scale by sqrt of eigen values
        //transform with eigenvectors and then add mean
        for (auto j = 0; j<6; j++)
        {
            for (auto k = j; k<6; k++)
            {
                s(k,j) = s(j, k) = Gaus[i].covar[j*6+k];
            }
        }
        s.Eigenvalvec(eigenval, eigenvec);
        for (auto j = 0; j < 6; j++) eigenval(j,0) = sqrt(eigenval(j,0));
        eigenvec = eigenvec.Inverse();
#if defined(USEOPENMP) && defined(USEOPENMPGPU)
        //if calculating using GPU, need to copy data to arrays
        for (auto j = 0; j < 6; j++) {
            eigenvalarray[6*noffset[i]+j] = eigenval(j,0);
            for (auto k = 0; k < 6; k++) {
                eigenvecarray[36*noffset[i]+6*j+k] = eigenvec(j,k);
            }
        }
        for (auto j = 0; j < Gaus[i].npoints; j++)
        {
            eigoffset[noffset[i]+j] = i;
            for (auto k = 0; k < 6; k++) newphase[(noffset[i]+j)*6+k] = Gaus[i].mean[k];
        }
#else
        //otherwise simply compute the position of the particles
        for (auto j = 0; j < Gaus[i].npoints; j++)
        {
            for (auto k = 0; k < 6; k++) y(0,k) = rn[(noffset[i]+j)*6+k];
            x = eigenvec * eigenval * y;
            for (auto k = 0; k < 6; k++) x(k,0) += Gaus[i].mean[k];
            for (auto k = 0; k < 3; k++) {
                if (x(k,0)<0) x(k,0) += opt.p;
                else if (x(k,0)>opt.p) x(k,0) -= opt.p;
            }
            for (auto k = 0; k < 6; k++) Part[noffset[i]+j].SetPhase(k, x(k,0));
        }
#endif
    }
#ifdef USEOPENMP
}
#endif


    //if calculating with GPUS, then use float arrays to calc positions
#if defined(USEOPENMP) && defined(USEOPENMPGPU)
#pragma omp target if (npoints>ompgpunum) map(to: opt.p, rn, eigoffset, eigenvalarray, eigenvecarray) map(newphase)
{
#pragma omp parallel for
    for (auto i = 0; i < npoints; i++)
    {
        Int_t index, eigindex1, eigindex2;
        index = 6*i;
        eigindex1 = eigoffset[i]*6;
        eigindex2 = eigoffset[i]*36;
        float val[6];
        for (auto j = 0; j < 6; j++)
        {
            val[j] = rn[index+j];
            for (auto k = 0; k<6; k++) {
                newphase[index+j] += eigenvecarray[eigindex2+6*j+k]*eigenvalarray[eigindex1+j]*val[j];
            }
            if (j<3) {
                if (newphase[index+j]>opt.p) newphase[index+j] -= opt.p;
                else if (newphase[index+j]<0) newphase[index+j] += opt.p;
            }
        }
    }
}
    delete[] eigoffset;
    delete[] eigenvalarray;
    delete[] eigenvecarray;
#pragma omp parallel for schedule(static)
    for (auto i = 0; i < npoints; i++)
    {
        Int_t index = 6*i;
        for (auto k = 0; k < 6; k++) Part[i].SetPhase(k, newphase[index+k]);
    }
    delete[] newphase;
#endif // end of if openmpgpu
// #endif //end of if openacc else
    delete[] rn;
    delete[] noffset;
    LOG(info)<<" Populate Gaussians: Done "<< local_timer;
    LOG(info)<<GetMemUsage(__func__+string("--line--")+to_string(__LINE__));
}

/// check distribution of particles in a mesh;
void PlaceGeneratedInMesh(Options &opt, vector<Particle> &Part)
{
    for (auto &p:Part) {
        MPIGetParticlesProcessor(opt, p.X(), p.Y(), p.Z());
    }
    unsigned long long *buff = new unsigned long long[opt.numcells];
    for (auto i=0;i<opt.numcells;i++) buff[i]=0;
    MPI_Allreduce(opt.cellnodenumparts.data(), buff, opt.numcells, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    for (auto i=0;i<opt.numcells;i++) opt.cellnodenumparts[i]=buff[i];
    delete[] buff;
    auto stats = statsofdata(opt.cellnodenumparts);
    LOG_RANK0(info)<< "The cell stats are "<<stats[1]<<" +/- "<<stats[2]<<" with cells with zero "<<stats[0];
    LOG_RANK0(info)<< "Quantiles (min,2.5,16,50,86,97.5,max) = "
    <<stats[3]<<" "
    <<stats[4]<<" "
    <<stats[5]<<" "
    <<stats[6]<<" "
    <<stats[7]<<" "
    <<stats[8]<<" "
    <<stats[9];

    vector<unsigned long long> mpinumparts(NProcs, 0);
    for (auto i=0;i<opt.numcells;i++)
    {
        auto itask = opt.cellnodeids[i];
        mpinumparts[itask] += opt.cellnodenumparts[i];
    }
    stats = statsofdata(mpinumparts);
    LOG_RANK0(info)<< "The MPI particle stats are "<<stats[1]<<" +/- "<<stats[2]<<" with number of MPI with zero "<<stats[0];
    LOG_RANK0(info)<< "Quantiles (min,2.5,16,50,86,97.5,max) = "
    <<stats[3]<<" "
    <<stats[4]<<" "
    <<stats[5]<<" "
    <<stats[6]<<" "
    <<stats[7]<<" "
    <<stats[8]<<" "
    <<stats[9];
}

//write a file
void WriteGeneratedInput(Options &opt, vector<Particle> &Part, vector<GaussianDistrib> &Gaus)
{
#ifdef USEHDF
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
    unsigned long long Nlocal = Part.size();
#endif
    vr::Timer local_timer;

    LOG(info)<<" Writing input ... ";

    fstream Fout;
    string fname;
    ostringstream os;
    char buf[40];
    long long unsigned nwritecommtot=0;
    long long unsigned noffset = 0;
    vector<unsigned long long> npart(6,0);
    vector<double> mass(6,0);
    void *data;
    vector<hsize_t> dims;
    hid_t datatype;

#ifdef USEMPI
    MPIBuildWriteComm(opt);
    MPI_Allreduce(&Nlocal, &nwritecommtot, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi_comm_write);
    if (NProcs > 1) {
        std::vector<Int_t> mpi_nparts(NProcsWrite);
        MPI_Allgather(&Nlocal, 1, MPI_Int_t, mpi_nparts.data(), 1, MPI_Int_t, mpi_comm_write);
        for (int j=0;j<ThisWriteTask;j++)noffset+=mpi_nparts[j];
    }
#endif

    H5OutputFile Fhdf;
    Fhdf.set_verbose(true);
    int ival;

    os << opt.outname <<".generatedinput";
#ifdef USEMPI
#ifdef USEPARALLELHDF
        os<<"."<<ThisWriteComm;
#else
        os<<"."<<ThisTask;
#endif
#else
#endif
    os <<".hdf5";
    fname = os.str();
#ifdef USEMPI
#endif
    LOG(info) << ThisWriteTask <<" and "<<ThisWriteComm<<" writing "<<fname<< "data is "<<Nlocal<<" with offset "<<noffset;
    LOG(info) << GetMemUsage(__func__+string("--line--")+to_string(__LINE__));
    LOG(info) <<"Expected memory needed to write data is max "<<3*sizeof(float)*Nlocal/1024./1024./1024.<<" GB";
    // create groups
    Fhdf.create(string(fname));
    Fhdf.close_group(Fhdf.create_group("Header"));
    Fhdf.close_group(Fhdf.create_group("Cosmology"));
    Fhdf.close_group(Fhdf.create_group("PartType1"));
    Fhdf.close_group(Fhdf.create_group("PartType1/Gaussians"));
    Fhdf.close();

#ifdef USEPARALLELHDF
    if(opt.mpinprocswritesize>1)
    {
        //if parallel then open file in serial so task 0 writes attributes
        Fhdf.append(string(fname),H5F_ACC_TRUNC, 0, false);
        npart[1] = nwritecommtot;
        //if parallel HDF then only
        if (ThisWriteTask==0) 
        {
            Fhdf.write_attribute(string("Header"), "BoxSize", opt.p);
            Fhdf.write_attribute(string("Header"), "NumFilesPerSnapshot", NWriteComms);
            npart[1] = nwritecommtot;
            Fhdf.write_attribute(string("Header"), "NumPart_ThisFile", npart);
            npart[1] = opt.Ngenerate;
            Fhdf.write_attribute(string("Header"), "NumPart_Total", npart);
            npart[1] = 0;
            Fhdf.write_attribute(string("Header"), "NumPart_Total_HighWord", npart);
            mass[1] = opt.mpgenerate;
            Fhdf.write_attribute(string("Header"), "MassTable", mass);
            Fhdf.write_attribute(string("Header"), "Redshift", opt.a);
            Fhdf.write_attribute(string("Header"), "Time", opt.a);
            Fhdf.write_attribute(string("Cosmology"), "Cosmological run", 1);
            Fhdf.write_attribute(string("Cosmology"), "Omega_m", opt.Omega_m);
            Fhdf.write_attribute(string("Cosmology"), "Omega_lambda", opt.Omega_Lambda);
            Fhdf.write_attribute(string("Cosmology"), "Omega_b", opt.Omega_b);
            Fhdf.write_attribute(string("Cosmology"), "Omega_r", opt.Omega_r);
            Fhdf.write_attribute(string("Cosmology"), "Omega_k", opt.Omega_k);
            Fhdf.write_attribute(string("Cosmology"), "h", opt.h);
            Fhdf.write_attribute(string("Header"), "Code", string("SWIFT"));
        }
        Fhdf.close();
        MPI_Barrier(mpi_comm_write);
        Fhdf.append(string(fname));
    }
    else{
        Fhdf.append(string(fname));
        Fhdf.write_attribute(string("/Header"), "BoxSize", opt.p);
        Fhdf.write_attribute(string("/Header"), "NumFilesPerSnapshot", NProcs);
        npart[1] = Nlocal;
        Fhdf.write_attribute(string("/Header"), "NumPart_ThisFile", npart);
        npart[1] = opt.Ngenerate;
        Fhdf.write_attribute(string("/Header"), "NumPart_Total", npart);
        npart[1] = 0;
        Fhdf.write_attribute(string("/Header"), "NumPart_Total_HighWord", npart);
        mass[1] = opt.mpgenerate;
        Fhdf.write_attribute(string("/Header"), "MassTable", mass);
        Fhdf.write_attribute(string("/Header"), "Redshift", opt.a);
        Fhdf.write_attribute(string("/Header"), "Time", opt.a);
        Fhdf.write_attribute(string("/Cosmology"), "Cosmological run", 1);
        Fhdf.write_attribute(string("/Cosmology"), "Omega_m", opt.Omega_m);
        Fhdf.write_attribute(string("/Cosmology"), "Omega_lambda", opt.Omega_Lambda);
        Fhdf.write_attribute(string("/Cosmology"), "Omega_b", opt.Omega_b);
        Fhdf.write_attribute(string("/Cosmology"), "Omega_r", opt.Omega_r);
        Fhdf.write_attribute(string("/Cosmology"), "Omega_k", opt.Omega_k);
        Fhdf.write_attribute(string("/Cosmology"), "h", opt.h);
        Fhdf.write_attribute(string("/Header"), "Code", string("SWIFT"));
    }
#else
    Fhdf.append(string(fname));
    Fhdf.write_attribute(string("/Header"), "BoxSize", opt.p);
    Fhdf.write_attribute(string("/Header"), "NumFilesPerSnapshot", NProcs);
    npart[1] = Nlocal;
    Fhdf.write_attribute(string("/Header"), "NumPart_ThisFile", npart);
    npart[1] = opt.Ngenerate;
    Fhdf.write_attribute(string("/Header"), "NumPart_Total", npart);
    npart[1] = 0;
    Fhdf.write_attribute(string("/Header"), "NumPart_Total_HighWord", npart);
    mass[1] = opt.mpgenerate;
    Fhdf.write_attribute(string("/Header"), "MassTable", mass);
    Fhdf.write_attribute(string("/Header"), "Redshift", opt.a);
    Fhdf.write_attribute(string("/Header"), "Time", opt.a);
    Fhdf.write_attribute(string("/Cosmology"), "Cosmological run", 1);
    Fhdf.write_attribute(string("/Cosmology"), "Omega_m", opt.Omega_m);
    Fhdf.write_attribute(string("/Cosmology"), "Omega_lambda", opt.Omega_Lambda);
    Fhdf.write_attribute(string("/Cosmology"), "Omega_b", opt.Omega_b);
    Fhdf.write_attribute(string("/Cosmology"), "Omega_r", opt.Omega_r);
    Fhdf.write_attribute(string("/Cosmology"), "Omega_k", opt.Omega_k);
    Fhdf.write_attribute(string("/Cosmology"), "h", opt.h);
    Fhdf.write_attribute(string("/Header"), "Code", string("SWIFT"));
#endif

    datatype = H5T_NATIVE_ULLONG;
    npart.clear();
    npart.resize(Nlocal);
    for (auto i=0; i<Nlocal; i++) npart[i] = noffset+i+1;
    Fhdf.write_dataset(opt, "PartType1/ParticleIDs", Nlocal, npart.data(), datatype);
    npart.clear();
    npart.shrink_to_fit();
    datatype = H5T_NATIVE_DOUBLE;
    data= ::operator new(sizeof(double)*(3*Nlocal));
    dims.resize(2);dims[0]=Nlocal;dims[1]=3;
    for (auto i=0;i<Nlocal;i++) {
        for (auto j=0;j<3;j++) {
            ((double*)data)[i*3+j]=Part[i].GetPosition(j);
        }
    }
    Fhdf.write_dataset_nd(opt, "PartType1/Coordinates", 2, dims.data(), data, datatype);
    for (auto i=0;i<Nlocal;i++) {
        for (auto j=0;j<3;j++) {
            ((double*)data)[i*3l+j]=Part[i].GetVelocity(j);
        }
    }
    Fhdf.write_dataset_nd(opt, "PartType1/Velocities", 2, dims.data(), data, datatype);
    ::operator delete(data);
    // mass.clear();
    // mass.resize(Nlocal,opt.mpgenerate);
    // Fhdf.write_dataset(opt, "/PartType1/Masses", Nlocal, mass.data(), datatype);
    /// write the halo information 
    datatype = H5T_NATIVE_FLOAT;
    data= ::operator new(sizeof(float)*(3*opt.Ngeneratehalos));
    dims.resize(2);dims[0]=opt.Ngeneratehalos;dims[1]=3;
    for (auto i=0;i<opt.Ngeneratehalos;i++) {
        for (auto j=0;j<3;j++) {
            ((float*)data)[i*3+j]=Gaus[i].mean[j];
        }
    }
    Fhdf.write_dataset_nd(opt, "PartType1/Gaussians/Coordinates", 2, dims.data(), data, datatype);
    for (auto i=0;i<opt.Ngeneratehalos;i++) {
        for (auto j=0;j<3;j++) {
            ((float*)data)[i*3+j]=Gaus[i].mean[j+3];
        }
    }
    Fhdf.write_dataset_nd(opt, "PartType1/Gaussians/Velocities", 2, dims.data(), data, datatype);
    ::operator delete(data);

    datatype = H5T_NATIVE_ULLONG;
    data= ::operator new(sizeof(unsigned long long)*(opt.Ngeneratehalos));
    dims.resize(1);dims[0]=opt.Ngeneratehalos;
    for (auto i=0;i<opt.Ngeneratehalos;i++) {
        ((unsigned long long*)data)[i]=Gaus[i].npoints;
    }
    Fhdf.write_dataset(opt, "PartType1/Gaussians/N", opt.Ngeneratehalos, data, datatype);
    ::operator delete(data);

    datatype = H5T_NATIVE_FLOAT;
    data= ::operator new(sizeof(float)*(opt.Ngeneratehalos));
    dims.resize(1);dims[0]=opt.Ngeneratehalos;
    for (auto i=0;i<opt.Ngeneratehalos;i++) {
        ((float*)data)[i]=Gaus[i].sigX;
    }
    Fhdf.write_dataset(opt, "PartType1/Gaussians/sigX", opt.Ngeneratehalos, data, datatype);
    for (auto i=0;i<opt.Ngeneratehalos;i++) {
        ((float*)data)[i]=Gaus[i].sigV;
    }
    Fhdf.write_dataset(opt, "PartType1/Gaussians/sigV", opt.Ngeneratehalos, data, datatype);
    ::operator delete(data);

    Fhdf.close();
    LOG(info)<<" Done writing "<<local_timer;
#endif
}

void VectorizationTest(Options &opt) 
{
    if (!opt.test_vectorization) return;
    vr::Timer total_timer;

    // silly compute to test performance 
    LOG(info)<<" Running silly calcs to test vectorization ";

    int nthreads = 1;
#ifdef USEOPENMP
    nthreads = omp_get_max_threads();
#endif

    vector<int> x_int, y_int;
    vector<float> x_float, y_float;
    vector<double> x_double, y_double;
    unsigned int Nentries = 20.0*1024.0*1024.0*1024.0/8.0/6.0;
    LOG(info)<<" Vectorization tests running with "<<Nentries;
    vr::Timer timer_mem;
    x_int.resize(Nentries);
    y_int.resize(Nentries);
    x_float.resize(Nentries);
    y_float.resize(Nentries);
    x_double.resize(Nentries);
    y_double.resize(Nentries);
    LOG(info)<<" mem alloc "<<timer_mem;
    LOG(info) << GetMemUsage(__func__+string("--line--")+to_string(__LINE__));
    vr::Timer timer_vec;

#pragma omp parallel for \
default(none) shared(x_int, y_int, x_float, y_float, x_double, y_double, Nentries) \
schedule(static) if (nthreads > 1)
   for (auto i=0u; i<Nentries; i++) {
      x_int[i] = i;
      auto temp = x_int[i];
      y_int[i] = temp+temp*pow(temp,2) + temp/(temp+1);
     //y_int[i] = x_int[i]+x_int[i]*pow(x_int[i],2) + x_int[i]/(x_int[i]+1);
//   }

//#pragma omp parallel for \
//default(none) shared(x_float, y_float, Nentries) \
//schedule(static) if (nthreads > 1) 
 //  for (auto i=0u; i<Nentries; i++) {
      x_float[i] = i;
      auto tempf = x_float[i];
      y_float[i] = tempf+tempf*pow(tempf,2) + tempf/(tempf+1);
//   }

//#pragma omp parallel for \
//default(none) shared(x_double, y_double, Nentries) \
//schedule(static) if (nthreads > 1)
//   for (auto i=0u; i<Nentries; i++) {
      x_double[i] = i;
      auto tempd = x_double[i];
      y_double[i] = tempd+tempd*pow(tempd,2) + tempd/(tempd+1);
   }
    LOG(info)<<" Finished "<< timer_vec;

    x_int.clear();
    x_float.clear();
    x_double.clear();
    y_int.clear();
    y_float.clear();
    y_double.clear();
    x_int.shrink_to_fit();
    x_float.shrink_to_fit();
    x_double.shrink_to_fit();
    y_int.shrink_to_fit();
    y_float.shrink_to_fit();
    y_double.shrink_to_fit();
 
    LOG(info)<<" Finished and mem cleared "<<timer_vec;
}
