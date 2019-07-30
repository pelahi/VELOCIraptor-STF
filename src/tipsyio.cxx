/*! \file tipsyio.cxx
 *  \brief this file contains routines for tipsy io
 */

//-- TIPSY SPECIFIC IO

#include "stf.h"

#include "tipsy_structs.h"
#include "endianutils.h"

///reads a tipsy file
void ReadTipsy(Options &opt, vector<Particle> &Part, const Int_t nbodies,Particle *&Pbaryons, Int_t nbaryons)
{
    struct tipsy_dump tipsyheader;
    struct tipsy_gas_particle gas;
    struct tipsy_dark_particle dark;
    struct tipsy_star_particle star;
    Int_t  count,oldcount,ngas,nstar,ndark, Ntot;
    double time,aadjust,z,Hubble,Hubbleflow,mtotold;
    double MP_DM=MAXVALUE, MP_B=MAXVALUE;
    int temp;
    Double_t mscale,lscale,lvscale,LN=1.0;
    Double_t posfirst[3];
    fstream Ftip;
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif

    //if MPI is used, Processor zero opens the file and loads the data into a particle buffer
    //this particle buffer is used to broadcast data to the appropriate processor
#ifdef USEMPI
    MPI_Status status;
    Particle *Pbuf;
    Int_t BufSize=opt.mpiparticlebufsize;
    Pbuf=new Particle[BufSize*NProcs];
    if (ThisTask==0) Pbuf=new Particle[BufSize*NProcs];
    Int_t Nlocalbuf,ibuf=0,*Nbuf;
    if (ThisTask==0) {
        Nbuf=new Int_t[NProcs];
        for (int j=0;j<NProcs;j++) Nbuf[j]=0;
    }
    Nlocal=0;
    MPIDomainExtentTipsy(opt);
    MPIInitialDomainDecomposition();
    if (ThisTask==0) {
#endif

    Ftip.open(opt.fname, ios::in | ios::binary);
    if (!Ftip){cerr<<"ERROR: Unable to open " <<opt.fname<<endl;exit(8);}
    else cout<<"Reading tipsy format from "<<opt.fname<<endl;

    InitEndian();

    //read tipsy header.
    Ftip.read((char*)&tipsyheader,sizeof(tipsy_dump));
    tipsyheader.SwitchtoBigEndian();
    Ftip.close();
    //offset stream by a double (time),  an integer (nbodies) ,integer (ndim), an integer (ngas)
    //read an integer (ndark), skip an integer (nstar), then data begins.
    time=tipsyheader.time;
    if ((opt.a-time)/opt.a>1e-2)
    {
        cout<<"Note that atime provided != to time in tipsy file (a,t): "<<opt.a<<","<<time<<endl;
        cout<<"Setting atime to that in file "<<endl;
        opt.a=time;
    }
    if (opt.comove) aadjust=1.0;
    else aadjust=opt.a;
    Ntot=tipsyheader.nbodies;
    ngas=tipsyheader.nsph;
    nstar=tipsyheader.nstar;
    ndark=tipsyheader.ndark;
    opt.numpart[GASTYPE]=ngas;
    opt.numpart[DARKTYPE]=ndark;
    opt.numpart[STARTYPE]=nstar;

    //Hubble flow and scale units
    z=1./opt.a-1.;
    Hubble=opt.h*opt.H*sqrt((1.0-opt.Omega_m-opt.Omega_Lambda)*pow(1.0+z,2.0)+opt.Omega_m*pow(1.0+z,3.0)+opt.Omega_Lambda);
    opt.rhobg=3.*Hubble*Hubble/8.0/M_PI/opt.G*opt.Omega_m;
    Double_t bnx=-((1-opt.Omega_m-opt.Omega_Lambda)*pow(aadjust,-2.0)+opt.Omega_Lambda)/((1-opt.Omega_m-opt.Omega_Lambda)*pow(aadjust,-2.0)+opt.Omega_m*pow(aadjust,-3.0)+opt.Omega_Lambda);
    opt.virBN98=(18.0*M_PI*M_PI+82.0*bnx-39*bnx*bnx)/opt.Omega_m;
    //if opt.virlevel<0, then use virial overdensity based on Bryan and Norman 1997 virialization level is given by
    if (opt.virlevel<0) opt.virlevel=opt.virBN98;

    mscale=opt.massinputconversion;lscale=opt.lengthinputconversion*aadjust;lvscale=opt.lengthinputconversion*opt.a;
    //normally Hubbleflow=lvscale*Hubble but we only care about peculiar velocities
    //ignore hubble flow
    Hubbleflow=0.;

    cout<<"File contains "<<Ntot<<" particles at is at time "<<opt.a<<endl;
    cout<<"There "<<ngas<<" gas, "<<ndark<<" dark, "<<nstar<<" stars."<<endl;
    cout<<"System to be searched contains "<<nbodies<<" particles of type "<<opt.partsearchtype<<" at time "<<opt.a<<endl;

    count=0;mtotold=0;

    //determine first particle about which use period
    if (opt.p>0) {
        count=0;
        Ftip.open(opt.fname, ios::in | ios::binary);
        Ftip.read((char*)&tipsyheader,sizeof(tipsy_dump));
        tipsyheader.SwitchtoBigEndian();
        for (Int_t i=0;i<ngas;i++)
        {
            Ftip.read((char*)&gas,sizeof(tipsy_gas_particle));
            gas.SwitchtoBigEndian();
            if ((opt.partsearchtype==PSTALL||opt.partsearchtype==PSTGAS)&&count==0) {
                posfirst[0]=gas.pos[0];posfirst[1]=gas.pos[1];posfirst[2]=gas.pos[2];
                count++;
                break;
            }
        }
        if (count==0) {
        for (Int_t i=0;i<ndark;i++)
        {
            Ftip.read((char*)&dark,sizeof(tipsy_dark_particle));
            dark.SwitchtoBigEndian();
            if ((opt.partsearchtype==PSTALL||opt.partsearchtype==PSTDARK)&&count==0) {
                posfirst[0]=dark.pos[0];posfirst[1]=dark.pos[1];posfirst[2]=dark.pos[2];
                count++;
                break;
            }
        }
        }
        if (count==0) {
        for (Int_t i=0;i<nstar;i++)
        {
            Ftip.read((char*)&star,sizeof(tipsy_star_particle));
            star.SwitchtoBigEndian();
            if ((opt.partsearchtype==PSTALL||opt.partsearchtype==PSTSTAR)&&count==0) {
                posfirst[0]=star.pos[0];posfirst[1]=star.pos[1];posfirst[2]=star.pos[2];
                count++;
                break;
            }
        }
        }
        Ftip.close();
    }

    oldcount=count=0;
    Ftip.open(opt.fname, ios::in | ios::binary);
    Ftip.read((char*)&tipsyheader,sizeof(tipsy_dump));
    tipsyheader.SwitchtoBigEndian();
    for (Int_t i=0;i<ngas;i++)
    {
        Ftip.read((char*)&gas,sizeof(tipsy_gas_particle));
        gas.SwitchtoBigEndian();
        //if particle is closer do to periodicity then alter position
        if (opt.p>0.0)
        {
            for (int j=0;j<3;j++) {
                if (gas.pos[j]-posfirst[j]>opt.p/2.0) gas.pos[j]-=opt.p;
                else if (gas.pos[j]-posfirst[j]<-opt.p/2.0) gas.pos[j]+=opt.p;
            }
        }
        if (opt.partsearchtype==PSTALL||opt.partsearchtype==PSTGAS) {
#ifndef USEMPI
            Part[count]=Particle(gas.mass*mscale,
                gas.pos[0]*lscale,gas.pos[1]*lscale,gas.pos[2]*lscale,
                gas.vel[0]*opt.velocityinputconversion+Hubbleflow*gas.pos[0],
                gas.vel[1]*opt.velocityinputconversion+Hubbleflow*gas.pos[1],
                gas.vel[2]*opt.velocityinputconversion+Hubbleflow*gas.pos[2],
                count,GASTYPE);
#else
            //if using MPI, determine ibuf, store particle in particle buffer and if buffer full, broadcast data
            //unless ibuf is 0, then just store locally
            ibuf=MPIGetParticlesProcessor(gas.pos[0],gas.pos[1],gas.pos[2]);
            Pbuf[ibuf*BufSize+Nbuf[ibuf]]=Particle(gas.mass*mscale,
                gas.pos[0]*lscale,gas.pos[1]*lscale,gas.pos[2]*lscale,
                gas.vel[0]*opt.velocityinputconversion+Hubbleflow*gas.pos[0],
                gas.vel[1]*opt.velocityinputconversion+Hubbleflow*gas.pos[1],
                gas.vel[2]*opt.velocityinputconversion+Hubbleflow*gas.pos[2],
                count,GASTYPE);
            Nbuf[ibuf]++;
            if(ibuf==0){
                Nbuf[ibuf]--;
                Part[Nlocal++]=Pbuf[ibuf*BufSize+Nbuf[ibuf]];
            }
            else {
                if(Nbuf[ibuf]==BufSize) {
                    MPI_Ssend(&Nbuf[ibuf], 1, MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
                    MPI_Ssend(&Pbuf[ibuf*BufSize],sizeof(Particle)*Nbuf[ibuf],MPI_BYTE,ibuf,ibuf,MPI_COMM_WORLD);
                    Nbuf[ibuf] = 0;
                }
            }
#endif
        count++;
        }
    }
    cout<<"Finished storing "<<count-oldcount<<" gas particles"<<endl;
    oldcount=count;
    for (Int_t i=0;i<ndark;i++)
    {
        Ftip.read((char*)&dark,sizeof(tipsy_dark_particle));
        dark.SwitchtoBigEndian();
        if (MP_DM>dark.mass) MP_DM=dark.mass;
        //if particle is closer do to periodicity then alter position
        if (opt.partsearchtype==PSTALL||opt.partsearchtype==PSTDARK) {
#ifndef USEMPI
        Part[count]=Particle(dark.mass*mscale,
            dark.pos[0]*lscale,dark.pos[1]*lscale,dark.pos[2]*lscale,
            dark.vel[0]*opt.velocityinputconversion+Hubbleflow*dark.pos[0],
            dark.vel[1]*opt.velocityinputconversion+Hubbleflow*dark.pos[1],
            dark.vel[2]*opt.velocityinputconversion+Hubbleflow*dark.pos[2],
            count,DARKTYPE);
#else
            ibuf=MPIGetParticlesProcessor(dark.pos[0],dark.pos[1],dark.pos[2]);
            Pbuf[ibuf*BufSize+Nbuf[ibuf]]=Particle(dark.mass*mscale,
                dark.pos[0]*lscale,dark.pos[1]*lscale,dark.pos[2]*lscale,
                dark.vel[0]*opt.velocityinputconversion+Hubbleflow*dark.pos[0],
                dark.vel[1]*opt.velocityinputconversion+Hubbleflow*dark.pos[1],
                dark.vel[2]*opt.velocityinputconversion+Hubbleflow*dark.pos[2],
                count,DARKTYPE);
            Nbuf[ibuf]++;
            if(ibuf==0){
                Nbuf[ibuf]--;
                Part[Nlocal++]=Pbuf[ibuf*BufSize+Nbuf[ibuf]];
            }
            else {
                if(Nbuf[ibuf]==BufSize) {
                    MPI_Ssend(&Nbuf[ibuf], 1, MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
                    MPI_Ssend(&Pbuf[ibuf*BufSize],sizeof(Particle)*Nbuf[ibuf],MPI_BYTE,ibuf,ibuf,MPI_COMM_WORLD);
                    Nbuf[ibuf] = 0;
                }
            }
#endif
        count++;
        }
    }
    cout<<"Finished storing "<<count-oldcount<<" dark particles"<<endl;
    oldcount=count;
    for (Int_t i=0;i<nstar;i++)
    {
        Ftip.read((char*)&star,sizeof(tipsy_star_particle));
        star.SwitchtoBigEndian();
        //if particle is closer do to periodicity then alter position
        if (opt.p>0.0)
        {
            for (int j=0;j<3;j++) {
                if (star.pos[j]-posfirst[j]>opt.p/2.0) star.pos[j]-=opt.p;
                else if (star.pos[j]-posfirst[j]<-opt.p/2.0) star.pos[j]+=opt.p;
            }
        }
        if (opt.partsearchtype==PSTALL||opt.partsearchtype==PSTSTAR) {
#ifndef USEMPI
        Part[count]=Particle(star.mass*mscale,
            star.pos[0]*lscale,star.pos[1]*lscale,star.pos[2]*lscale,
            star.vel[0]*opt.velocityinputconversion+Hubbleflow*star.pos[0],
            star.vel[1]*opt.velocityinputconversion+Hubbleflow*star.pos[1],
            star.vel[2]*opt.velocityinputconversion+Hubbleflow*star.pos[2],
            count,STARTYPE);
#else
            ibuf=MPIGetParticlesProcessor(star.pos[0],star.pos[1],star.pos[2]);
            Pbuf[ibuf*BufSize+Nbuf[ibuf]]=Particle(star.mass*mscale,
                star.pos[0]*lscale,star.pos[1]*lscale,star.pos[2]*lscale,
                star.vel[0]*opt.velocityinputconversion+Hubbleflow*star.pos[0],
                star.vel[1]*opt.velocityinputconversion+Hubbleflow*star.pos[1],
                star.vel[2]*opt.velocityinputconversion+Hubbleflow*star.pos[2],
                count,STARTYPE);
            Nbuf[ibuf]++;
            if(ibuf==0){
                Nbuf[ibuf]--;
                Part[Nlocal++]=Pbuf[ibuf*BufSize+Nbuf[ibuf]];
            }
            else {
                if(Nbuf[ibuf]==BufSize) {
                    MPI_Ssend(&Nbuf[ibuf], 1, MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
                    MPI_Ssend(&Pbuf[ibuf*BufSize],sizeof(Particle)*Nbuf[ibuf],MPI_BYTE,ibuf,ibuf,MPI_COMM_WORLD);
                    Nbuf[ibuf]=0;
                }
            }
#endif
        count++;
        }
    }

    cout<<"Finished storing "<<count-oldcount<<" star particles"<<endl;
    //once finished reading the file if there are any particles left in the buffer broadcast them
#ifdef USEMPI
    for(ibuf = 1; ibuf < NProcs; ibuf++)
    {
        while(Nbuf[ibuf])
        {
            MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
            MPI_Ssend(&Pbuf[ibuf*BufSize], sizeof(Particle)*Nbuf[ibuf], MPI_BYTE, ibuf, ibuf, MPI_COMM_WORLD);
            Nbuf[ibuf]=0;
        }
        //last broadcast with Nbuf[ibuf]=0 so that receiver knows no more particles are to be broadcast
        MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t,ibuf,ibuf+NProcs,MPI_COMM_WORLD);
    }
#endif
    Ftip.close();

#ifdef USEMPI
    }
    else {
        do {
            MPI_Recv(&Nlocalbuf, 1, MPI_Int_t, 0, ThisTask+NProcs, MPI_COMM_WORLD, &status);
            if(Nlocalbuf) {
                MPI_Recv(&Part[Nlocal],sizeof(Particle)*Nlocalbuf,MPI_BYTE,0,ThisTask, MPI_COMM_WORLD,&status);
                Nlocal+=Nlocalbuf;
            }
        } while(Nlocalbuf);
    }
    //store global ids and index locally
    //for (Int_t i=0;i<Nlocal;i++) mpi_indexlist[i]=Part[i].GetID();
    //for (Int_t i=0;i<Nlocal;i++) mpi_idlist[i]=Part[i].GetPID();
    //if (ThisTask==0) {
        //MPI_Bcast(&cm,sizeof(Coordinate),MPI_BYTE,0,MPI_COMM_WORLD);
        //MPI_Bcast(&cmvel,sizeof(Coordinate),MPI_BYTE,0,MPI_COMM_WORLD);
    //}
#endif

    //calculate the interparticle spacing
#ifdef HIGHRES
    if (opt.Neff==-1) {
        //Once smallest mass particle is found (which should correspond to highest resolution area,
        LN=pow(((MP_DM)*opt.massinputconversion)/(opt.Omega_cdm*3.0*opt.H*opt.h*opt.H*opt.h/(8.0*M_PI*opt.G)),1./3.)*opt.a;
    }
    else {
        LN=opt.p/(Double_t)opt.Neff;
    }
#endif
    opt.internalenergyinputconversion = opt.velocityinputconversion*opt.velocityinputconversion;

    //adjust physical scales by the inferred interparticle spacing
    opt.ellxscale=LN;
    opt.uinfo.eps*=LN;

}
