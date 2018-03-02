/*! \file swiftinterface.cxx
 *  \brief this file contains routines that allow the velociraptor library to interface with the swift N-body code from within swift.
 */


#include "swiftinterface.h"

Options libvelociraptorOpt;
void InitVelociraptor(char* configname, char* outputname, cosmoinfo c, unitinfo u, siminfo s)
{

    cout<<"Initialising VELOCIraptor..."<< endl;

    libvelociraptorOpt.pname = configname;
    libvelociraptorOpt.outname = outputname;

    cout<<"Reading VELOCIraptor config file..."<< endl;
    GetParamFile(libvelociraptorOpt);
    cout<<"Setting cosmology, units, sim stuff "<<endl;
    ///set units
    libvelociraptorOpt.lengthtokpc=u.lengthtokpc;
    libvelociraptorOpt.velocitytokms=u.velocitytokms;
    libvelociraptorOpt.masstosolarmass=u.masstosolarmass;
    libvelociraptorOpt.L=1.0;
    libvelociraptorOpt.M=1.0;
    libvelociraptorOpt.V=1.0;
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
    if (libvelociraptorOpt.icosmologicalin) {
        libvelociraptorOpt.ellxscale=s.interparticlespacing;
        libvelociraptorOpt.uinfo.eps*=libvelociraptorOpt.ellxscale;
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
    cout<<"Finished initialising VELOCIraptor"<<endl;
}

void InvokeVelociraptor(const int num_gravity_parts, struct gpart *gravity_parts) {

    Particle *parts;
    parts=new Particle[num_gravity_parts];

    cout<<"Copying particle data..."<< endl;
    for(int i=0; i<num_gravity_parts; i++) {
        parts[i] = Particle(&gravity_parts[i]);

        //cout<<"Particle " << i << " ID: " << parts[i].GetPID() << ", gpart ID: " << gravity_parts[i].id_or_neg_offset << endl;

    }
    cout<<"Finished copying particle data."<< endl;

#ifdef USEMPI
    //find out how big the SPMD world is
    MPI_Comm_size(MPI_COMM_WORLD,&NProcs);
    mpi_domain=new MPI_Domain[NProcs];
#endif

    Int_t *pfof, ngroup;

    pfof=SearchFullSet(libvelociraptorOpt,num_gravity_parts,parts,ngroup);

    cout<<"Finished searching for haloes, no. of groups: " << ngroup << ". Freeing particle data."<< endl;

    delete[] parts;

    cout<<"VELOCIraptor returning."<< endl;
}
