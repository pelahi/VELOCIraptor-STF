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
    LOG(info)<<" Generating Input ... ";
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
    opt.Ngeneratehalos = opt.Ngeneratehalos/(float)NProcs;
    if (NProcs > 1) {
        if (ThisTask == NProcs - 1)
        {
            Nlocal = opt.Ngenerate*Nlocal*(NProcs-1);
        }
    }
    opt.Nbackground = opt.fbackground * opt.Ngenerate;
    Part.resize(Nlocal);
    for (auto p:Part) p.SetMass(opt.mpgenerate);
    LOG(info) << GetMemUsage(__func__+string("--line--")+to_string(__LINE__));

    vector<GaussianDistrib> Gaus = ProduceGaussians(opt, Part);
    Int_t npoints = 0;
    for (auto &x:Gaus) npoints += x.npoints;
    //update the background number of points
    opt.Nbackground = Nlocal - npoints;
    LOG(info)<<" producing positions for "<<opt.Ngeneratehalos<<" locally containing "<<npoints<<" particles and a background of "<<opt.Nbackground;
#ifdef USEMPI
#endif
    opt.fbackground = opt.Nbackground/(double)opt.Ngenerate;
    PopulateGaussians(opt, Part, Gaus);
    ProduceBackground(opt, Part, npoints);
    //WriteGeneratedInput(opt, Part, Gaus);
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
    default_random_engine generator;
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
    Int_t maxnpart = max((opt.Ngenerate-opt.Nbackground)/(100.0*opt.Ngeneratehalos),100.0);
    Int_t npoints = 0;
    default_random_engine generator;
    double mmin = log(20.0), mmax = log(maxnpart);
    uniform_real_distribution<double> npart(mmin, mmax);
    for (auto i = 0; i < opt.Ngeneratehalos; i++) {
        if (npoints + exp(mmax)/2 > Nlocal) {
            opt.Ngeneratehalos = i;
            Gaus.resize(opt.Ngeneratehalos);
            break;
        }
        do {
            Gaus[i].npoints = exp(npart(generator));
        } while (Gaus[i].npoints + npoints > Nlocal);
        npoints += Gaus[i].npoints;
    }

#if defined(USEOPENMP)
#pragma omp parallel default(shared) private(generator)
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
    float *rn = new float[npoints*6];
#if defined(USEOPENMP)
#pragma omp parallel default(shared)
{
#endif
    default_random_engine generator;
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
    vector<unsigned long long> npart(6,0);
    vector<double> mass(6,0);
    void *data;
    vector<hsize_t> dims;
    hid_t datatype;

#ifdef USEMPI
    MPIBuildWriteComm(opt);
#endif

    H5OutputFile Fhdf;
    int itemp=0;
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

#ifdef USEPARALLELHDF
        if(opt.mpinprocswritesize>1){
            //if parallel then open file in serial so task 0 writes header
            Fhdf.create(string(fname),H5F_ACC_TRUNC, 0, false);
        }
        else{
             Fhdf.create(string(fname),H5F_ACC_TRUNC, ThisWriteComm, false);
        }
#else
        Fhdf.create(string(fname));
#endif
        itemp=0;

#ifdef USEPARALLELHDF
        if(opt.mpinprocswritesize>1){
            MPI_Allreduce(&Nlocal, &nwritecommtot, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi_comm_write);
            npart[1] = nwritecommtot;
            //if parallel HDF then only
            if (ThisWriteTask==0) {
                Fhdf.create_group("Header");
                Fhdf.create_group("Cosmology");
                Fhdf.create_group("PartType1");
                Fhdf.create_group("PartType1/Gaussians");
                Fhdf.write_attribute(string("/Header"), "BoxSize", opt.p);
                Fhdf.write_attribute(string("/Header"), "NumFilesPerSnapshot", NWriteComms);
                npart[1] = nwritecommtot;
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
            Fhdf.close();
            MPI_Barrier(MPI_COMM_WORLD);
            //reopen for parallel write
            Fhdf.append(string(fname));
        }
        else{
            Fhdf.create_group("Header");
            Fhdf.create_group("Cosmology");
            Fhdf.create_group("PartType1");
            Fhdf.create_group("PartType1/Gaussians");
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
        Fhdf.create_group("Header");
        Fhdf.create_group("Cosmology");
        Fhdf.create_group("PartType1");
        Fhdf.create_group("PartType1/Gaussians");
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
        datatype = H5T_NATIVE_FLOAT;
        data= ::operator new(sizeof(float)*(3*Nlocal));
        dims.resize(2);dims[0]=Nlocal;dims[1]=3;
        for (auto i=0;i<Nlocal;i++) {
            for (auto j=0;j<3;j++) {
                ((float*)data)[i*3+j]=Part[i].GetPosition(j);
            }
        }
        Fhdf.write_dataset_nd(opt, "/PartType1/Coordinates", 2, dims.data(), data, datatype);
        for (auto i=0;i<Nlocal;i++) {
            for (auto j=0;j<3;j++) {
                ((float*)data)[i*3l+j]=Part[i].GetVelocity(j);
            }
        }
        Fhdf.write_dataset_nd(opt, "/PartType1/Velocities", 2, dims.data(), data, datatype);
        ::operator delete(data);
        mass.clear();
        mass.resize(Nlocal,opt.mpgenerate);
        Fhdf.write_dataset(opt, "/PartType1/Masses", Nlocal, mass.data(), datatype);
        datatype = H5T_NATIVE_ULLONG;
        npart.clear();
        npart.resize(Nlocal);
        for (auto i=0; i<Nlocal; i++) npart[i] = i+1;
        Fhdf.write_dataset(opt, "/PartType1/ParticleIDs", Nlocal, npart.data(), datatype);

        datatype = H5T_NATIVE_FLOAT;
        data= ::operator new(sizeof(float)*(3*opt.Ngeneratehalos));
        dims.resize(2);dims[0]=opt.Ngeneratehalos;dims[1]=3;
        for (auto i=0;i<opt.Ngeneratehalos;i++) {
            for (auto j=0;j<3;j++) {
                ((float*)data)[i*3+j]=Gaus[i].mean[j];
            }
        }
        Fhdf.write_dataset_nd(opt, "/PartType1/Gaussians/Coordinates", 2, dims.data(), data, datatype);
        for (auto i=0;i<opt.Ngeneratehalos;i++) {
            for (auto j=0;j<3;j++) {
                ((float*)data)[i*3+j]=Gaus[i].mean[j+3];
            }
        }
        Fhdf.write_dataset_nd(opt, "/PartType1/Gaussians/Velocities", 2, dims.data(), data, datatype);
        ::operator delete(data);

        datatype = H5T_NATIVE_ULLONG;
        data= ::operator new(sizeof(unsigned long long)*(opt.Ngeneratehalos));
        dims.resize(1);dims[0]=opt.Ngeneratehalos;
        for (auto i=0;i<opt.Ngeneratehalos;i++) {
            ((unsigned long long*)data)[i]=Gaus[i].npoints;
        }
        Fhdf.write_dataset(opt, "/PartType1/Gaussians/N", opt.Ngeneratehalos, data, datatype);
        ::operator delete(data);

        datatype = H5T_NATIVE_FLOAT;
        data= ::operator new(sizeof(float)*(opt.Ngeneratehalos));
        dims.resize(1);dims[0]=opt.Ngeneratehalos;
        for (auto i=0;i<opt.Ngeneratehalos;i++) {
            ((float*)data)[i]=Gaus[i].sigX;
        }
        Fhdf.write_dataset(opt, "/PartType1/Gaussians/sigX", opt.Ngeneratehalos, data, datatype);
        for (auto i=0;i<opt.Ngeneratehalos;i++) {
            ((float*)data)[i]=Gaus[i].sigV;
        }
        Fhdf.write_dataset(opt, "/PartType1/Gaussians/sigV", opt.Ngeneratehalos, data, datatype);
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

    vector<int> x_int, y_int;
    vector<float> x_float, y_float;
    vector<double> x_double, y_double;
    unsigned int Nentries = 120.0*1024.0*1024.0*1024.0/8.0/6.0;
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
    exit(9);
}