/*! \file generate.cxx
 *  \brief this file contains routines to generate a test input file
 */


#include "stf.h"
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
    if (NProcs > 1) {
        if (ThisTask == NProcs - 1)
        {
            Nlocal = opt.Ngenerate*Nlocal*(NProcs-1);
        }
    }
    opt.Nbackground = opt.fbackground * opt.Ngenerate;
    Part.resize(Nlocal);
    for (auto p:Part) p.SetMass(opt.mpgenerate);

    cout<<ThisTask<<" "<<opt.Ngenerate<<" "<<Nlocal<<endl;
    GetMemUsage(opt, __func__+string("--line--")+to_string(__LINE__), (opt.iverbose>=0));

    vector<GaussianDistrib> Gaus = ProduceGaussians(opt);
    PopulateGaussians(opt, Part, Gaus);
    ProduceBackground(opt, Part, opt.Nbackground);
    WriteGeneratedInput(opt, Part, Gaus);
}

/// produce a background uniform distribution of particles within the volume
void ProduceBackground(Options &opt, vector<Particle> &Part, Int_t noffset)
{
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
    GetMemUsage(opt, __func__+string("--line--")+to_string(__LINE__), (opt.iverbose>=0));
}

inline Matrix NormalizeDisp(double sX[3], double ssX[3], double norm, uniform_real_distribution<double> &disp, default_random_engine& generator){
    for (auto j = 0; j < 3; j++) {
        do {
            sX[j] = disp(generator);
        }while (sX[j]==0);
        ssX[j] = pow(disp(generator),2.0);
    }
    Matrix X;
    for (auto j = 0; j < 3; j++) {
        for (auto k = j+1; k < 3; k++) {
            X(j,k) = X(k,j) = ssX[j]*sX[j]*sX[k];
        }
    }
    for (auto j = 0; j < 3; j++) X(j,j) = sX[j]*sX[j];
    norm /= pow(X.Det(), 1.0/6.0);
    X = X*norm;
    return X;

}

/// produce Gaussian
vector<GaussianDistrib> ProduceGaussians(Options &opt)
{
    vector<GaussianDistrib> Gaus;
    if (opt.fbackground == 1 || opt.Ngeneratehalos == 0) return Gaus;
    Gaus.resize(opt.Ngeneratehalos);
    Int_t maxnpart = max((opt.Ngenerate-opt.Nbackground)/(100.0*opt.Ngeneratehalos),100.0);

#if defined(USEOPENMP)
#pragma omp parallel default(shared)
{
#endif
    default_random_engine generator;
    uniform_real_distribution<double> pos(0.0, opt.p);
    uniform_real_distribution<double> vel(-opt.H, opt.H);
    uniform_real_distribution<double> scaling(0.1, 0.5);
    uniform_real_distribution<double> disp(-1.0, 1.0);
    double mmin = log(20.0), mmax = log(maxnpart);
    uniform_real_distribution<double> npart(mmin, mmax);
    Coordinate x,v;
    double mass, sigX, sigV, sX[3], ssX[3], sV[3], ssV[3];
    Matrix X, V;
#if defined(USEOPENMP)
    #pragma omp for schedule(static)
#endif
    for (auto i = 0; i < opt.Ngeneratehalos; i++) {
        for (auto j = 0; j < 3; j++) {
            Gaus[i].mean[j] = pos(generator);
            Gaus[i].mean[j+3] = vel(generator);
        }
        for (auto j = 0; j < 36; j++) Gaus[i].covar[j] = 0;
        Gaus[i].npoints = exp(npart(generator));
        mass = Gaus[i].npoints*opt.mpgenerate;
        sigX = pow(mass/(200.0*4.0*M_PI/3.0*opt.rhocrit),1.0/3.0)*scaling(generator);
        sigV = sqrt(opt.G*mass/sigX);
        X = NormalizeDisp(sX, ssX, sigX, disp, generator);
        V = NormalizeDisp(sV, ssV, sigV, disp, generator);
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
    GetMemUsage(opt, __func__+string("--line--")+to_string(__LINE__), (opt.iverbose>=0));
    return Gaus;

}

/// populate gauassians
void PopulateGaussians(Options &opt, vector<Particle> &Part, vector<GaussianDistrib> &Gaus)
{
    if (Gaus.size() == 0) return;
    // get random numbers first then transform them as necessary
    Int_t npoints = 0;
    vector<Int_t> noffset(Gaus.size(),0);
    for (auto i = 0; i < opt.Ngeneratehalos; i++)
    {
        noffset[i] = npoints;
        npoints += Gaus[i].npoints;
    }
    vector<double> rn(npoints*6);
#if defined(USEOPENMP)
#pragma omp parallel default(shared)
{
#endif
    default_random_engine generator;
    normal_distribution<double> g(0.0, 1.0);
#if defined(USEOPENMP)
    #pragma omp for schedule(static)
#endif
    for (auto i = 0; i < rn.size(); i++) {
        rn[i] = g(generator);
    }
#if defined(USEOPENMP)
}
#endif
    // now transform points based on mean and variance
    //
// #if defined(USEOPENACC)
// #else
#if defined(USEOPENMP)
#pragma omp parallel default(shared)
{
#endif
    GMatrix y(1,6), x(6,1), s(6,6), eigenval(6,1), eigenvec(6,6);
#if defined(USEOPENMP)
    #pragma omp for schedule(static)
#endif
    for (auto i = 0; i < opt.Ngeneratehalos; i++) {
        //get the eigenvalues and eigenvectors of the covariance matrix
        //then take normal random numbers and scale by sqrt of eigen values
        //transform with eigenvectors and then add mean
        for (auto j = 0; j<6; j++)
        {
            for (auto k = 0; k<6; k++)
            {
                s(j, k) = Gaus[i].covar[j*6+k];
            }
        }
        s.Eigenvalvec(eigenval, eigenvec);
        for (auto j = 0; j < 6; j++) eigenval(j,0) = sqrt(eigenval(j,0));
        for (auto j = 0; j < Gaus[i].npoints; j++)
        {
            for (auto k = 0; k < 6; k++) y(0,k) = rn[(noffset[i]+j)*6+k];
            x = eigenvec * eigenval * y;
            for (auto k = 0; k < 6; k++)
            {
                Part[noffset[i]+j].SetPhase(k, x(k,0));
            }
        }
    }
#ifdef USEOPENMP
}
#endif
// #endif
    opt.Nbackground = opt.Ngenerate - npoints;
    opt.fbackground = opt.Nbackground/(double)opt.Ngenerate;
    GetMemUsage(opt, __func__+string("--line--")+to_string(__LINE__), (opt.iverbose>=0));
}

//write a file
void WriteGeneratedInput(Options &opt, vector<Particle> &Part, vector<GaussianDistrib> &Gaus)
{
#ifdef USEHDF
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
    unsigned long long Nlocal = Part.size();
#endif

    fstream Fout;
    string fname;
    ostringstream os;
    char buf[40];
    long long unsigned nwritecommtot=0;
    vector<unsigned long long> npart(6,0);
    vector<double> mass(6,0);

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
        // names[itemp++]=string("Header/BoxSize");
        // names[itemp++]=string("Header/MassTable");
        // names[itemp++]=string("Header/NumPart_ThisFile");
        // names[itemp++]=string("Header/NumPart_Total");
        // names[itemp++]=string("Header/NumPart_Total_HighWord");
        // names[itemp++]=string("Cosmology/Omega_m");
        // names[itemp++]=string("Cosmology/Omega_lambda");
        // names[itemp++]=string("Header/Redshift");
        // names[itemp++]=string("Header/Time");
        // names[itemp++]=string("Header/NumFilesPerSnapshot");
        // names[itemp++]=string("Cosmology/h");
        // names[itemp++]=string("Cosmology/Cosmological run");
        // names[itemp++]=string("Cosmology/Omega_b");
        // names[itemp++]=string("Cosmology/Omega_r");
        // names[itemp++]=string("Cosmology/Omega_k");
        // names[itemp++]=string("Cosmology/w");
        // names[itemp++]=string("Cosmology/w_0");
        // names[itemp++]=string("Cosmology/w_a");


#ifdef USEPARALLELHDF
        if(opt.mpinprocswritesize>1){
            MPI_Allreduce(&Nlocal, &nwritecommtot, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi_comm_write);
            npart[1] = nwritecommtot;
            //if parallel HDF then only
            if (ThisWriteTask==0) {
                Fhdf.write_attribute(string("/Header/"), "BoxSize", opt.p);
                Fhdf.write_attribute(string("/Header/"), "NumFilesPerSnapshot", NWriteComms);
                npart[1] = nwritecommtot;
                Fhdf.write_attribute(string("/Header/"), "NumPart_ThisFile", npart);
                npart[1] = opt.Ngenerate;
                Fhdf.write_attribute(string("/Header/"), "NumPart_Total", npart);
                npart[1] = 0;
                Fhdf.write_attribute(string("/Header/"), "NumPart_Total_HighWord", npart);
            }
            Fhdf.close();
            MPI_Barrier(MPI_COMM_WORLD);
            //reopen for parallel write
            Fhdf.append(string(fname));
        }
        else{
            Fhdf.write_attribute(string("/Header/"), "BoxSize", opt.p);
            Fhdf.write_attribute(string("/Header/"), "NumFilesPerSnapshot", NProcs);
            npart[1] = Nlocal;
            Fhdf.write_attribute(string("/Header/"), "NumPart_ThisFile", npart);
            npart[1] = opt.Ngenerate;
            Fhdf.write_attribute(string("/Header/"), "NumPart_Total", npart);
            npart[1] = 0;
            Fhdf.write_attribute(string("/Header/"), "NumPart_Total_HighWord", npart);
        }
#else
        Fhdf.write_attribute(string("/Header/"), "BoxSize", opt.p);
        Fhdf.write_attribute(string("/Header/"), "NumFilesPerSnapshot", NProcs);
        npart[1] = Nlocal;
        Fhdf.write_attribute(string("/Header/"), "NumPart_ThisFile", npart);
        npart[1] = opt.Ngenerate;
        Fhdf.write_attribute(string("/Header/"), "NumPart_Total", npart);
        npart[1] = 0;
        Fhdf.write_attribute(string("/Header/"), "NumPart_Total_HighWord", npart);
        // Fhdf.write_dataset(opt, datagroupnames.prop[itemp++], 1, &ThisTask);
        // Fhdf.write_dataset(opt, datagroupnames.prop[itemp++], 1, &NProcs);
        // Fhdf.write_dataset(opt, datagroupnames.prop[itemp++], 1, &ng);
        // Fhdf.write_dataset(opt, datagroupnames.prop[itemp++], 1, &ngtot);
        // Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.icosmologicalin);
        // Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.icomoveunit);
        // Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.p);
        // Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.a);
        // Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.lengthtokpc);
        // Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.velocitytokms);
        // Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.masstosolarmass);
        Fhdf.close();
#endif

#endif
}
