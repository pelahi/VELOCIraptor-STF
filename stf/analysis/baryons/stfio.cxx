/*! \file stfio.cxx
 *  \brief this file contains routines for io for reading in STF/VELOCIraptor output
 */

#include "baryoniccontent.h"
/// \name read the output produced by STF
//@{

///Read halo data from an idividual snapshot;
///\todo May want to adjust stio output so that one also outputs the particle ids making indexing easy but this is not trivial
///as particles get moved around all the time. 
HaloParticleData *ReadHaloGroupCatalogData(char* infile, Int_t &numhalos, Int_t &nbodies,int mpi_ninput, int ibinary, int ifieldhalos)
{
    cout<<"Read halo group and particle catalogs"<<endl;
    HaloParticleData *Halo;
    Int_t itemp, noffset,noldhalopartoffset;
    long unsigned nparts,haloid;
    long unsigned TotalNumberofHalos;
    char fname1[1000],fname2[1000],fname3[1000];
    char fname4[1000],fname5[1000],fname6[1000];
    fstream Fgroup,Fpart,Fupart; //field objects
    fstream Fsgroup,Fspart,Fsupart; //sublevels

    Int_t nids,nuids,nsids,nsuids,nglocal,nsglocal,nidsoff,nn;
    Int_t nidstot,nuidstot,nsidstot,nsuidstot;
    Int_t *numingroup,*numingroupbound,*offset,*uoffset;
    int j=0,nmpicount,itask,nprocs;
    Int_t *idval;
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif
    //initialize number of particles in structures
    nbodies=0;

    ///\todo need to fully implement MPI routines. At the moment, only ThisTask==0 does anything.
    //if (ThisTask==0) {

    ///adjust files read based on whether mpi output or not
    if (mpi_ninput==0) nmpicount=1;
    else nmpicount=mpi_ninput;
    //to offset the halo structure pointer given that data maybe split among multiple files
    noffset=0;
    for (j=0;j<nmpicount;j++) {
    nglocal=nsglocal=0;
    //load the group catalogues and the related particle files
    if (mpi_ninput==0) {
        sprintf(fname1,"%s.catalog_groups",infile);
        sprintf(fname2,"%s.catalog_particles",infile);
        sprintf(fname3,"%s.catalog_particles.unbound",infile);
        sprintf(fname4,"%s.sublevels.catalog_groups",infile);
        sprintf(fname5,"%s.sublevels.catalog_particles",infile);
        sprintf(fname6,"%s.sublevels.catalog_particles.unbound",infile);
    }
    else {
        sprintf(fname1,"%s.catalog_groups.%d",infile,j);
        sprintf(fname2,"%s.catalog_particles.%d",infile,j);
        sprintf(fname3,"%s.catalog_particles.unbound.%d",infile,j);
        sprintf(fname4,"%s.sublevels.catalog_groups.%d",infile,j);
        sprintf(fname5,"%s.sublevels.catalog_particles.%d",infile,j);
        sprintf(fname6,"%s.sublevels.catalog_particles.unbound.%d",infile,j);
    }

    if(ibinary) {
    Fgroup.open(fname1,ios::in|ios::binary);
    Fpart.open(fname2,ios::in|ios::binary);
    Fupart.open(fname3,ios::in|ios::binary);
    if (ifieldhalos) {
    Fsgroup.open(fname4,ios::in|ios::binary);
    Fspart.open(fname5,ios::in|ios::binary);
    Fsupart.open(fname6,ios::in|ios::binary);
    }
    }
    else {
    Fgroup.open(fname1,ios::in);
    Fpart.open(fname2,ios::in);
    Fupart.open(fname3,ios::in);
    if(ifieldhalos) {
    Fsgroup.open(fname4,ios::in);
    Fspart.open(fname5,ios::in);
    Fsupart.open(fname6,ios::in);
    }
    }
    //check header make sure that number of mpi values considered agree
    if (ibinary) {
    Fgroup.read((char*)&itask,sizeof(int));
    Fgroup.read((char*)&nprocs,sizeof(int));
    if(ifieldhalos){
    Fsgroup.read((char*)&itask,sizeof(int));
    Fsgroup.read((char*)&nprocs,sizeof(int));
    }
    }
    else {
    Fgroup>>itask>>nprocs;
    if(ifieldhalos) {
    Fsgroup>>itask>>nprocs;
    }
    }
    if (nprocs!=mpi_ninput&&mpi_ninput>0) {
        cout<<"Error, number of mpi outputs was set to "<<mpi_ninput<<" but file indicates there are "<<nprocs<<endl;
        /*
        cout<<"Terminating"<<endl;
#ifdef USEMPI
        MPI_Finalize();
#endif
        exit(9);
        */
        nmpicount=mpi_ninput=nprocs;
    }
    else if (!(mpi_ninput==0&&nprocs==1)) {
        cout<<"Error, number of mpi outputs was set to zero but file indicates there are more than one mpi output"<<endl;
        cout<<"Terminating"<<endl;
        cout<<mpi_ninput<<" "<<nprocs<<endl;
#ifdef USEMPI
        MPI_Finalize();
#endif
        exit(9);
    }

    if (ibinary) {
    Fgroup.read((char*)&itemp,sizeof(Int_t));
    nglocal=itemp;
    Fgroup.read((char*)&itemp,sizeof(Int_t));
    TotalNumberofHalos=itemp;
    if(ifieldhalos) {
    Fsgroup.read((char*)&itemp,sizeof(Int_t));
    nsglocal=itemp;
    Fsgroup.read((char*)&itemp,sizeof(Int_t));
    }
    }
    else {
    Fgroup>>nglocal;
    Fgroup>>TotalNumberofHalos;
    if(ifieldhalos) {
    Fsgroup>>nsglocal;
    Fsgroup>>itemp;
    }
    }

    cout<<infile<<" has "<<nglocal<<endl;
    //allocate memory for halos
    if (j==0) {
        cout<<" and the total number of halos in all files is "<<TotalNumberofHalos<<endl;
        Halo=new HaloParticleData[TotalNumberofHalos];
        numhalos=TotalNumberofHalos;
        Halo[0].noffset=noldhalopartoffset=0;
    }

    //read the header info before reading lengths of the halos
    if(ibinary) {
    Fpart.read((char*)&itask,sizeof(int));
    Fpart.read((char*)&nprocs,sizeof(int));
    Fupart.read((char*)&itask,sizeof(int));
    Fupart.read((char*)&nprocs,sizeof(int));

    Fpart.read((char*)&nids,sizeof(Int_t));
    Fpart.read((char*)&nidstot,sizeof(Int_t));//actually total number across all files
    Fupart.read((char*)&nuids,sizeof(Int_t));
    Fupart.read((char*)&nuidstot,sizeof(Int_t));//total across all files
    if(ifieldhalos) {
    Fspart.read((char*)&itask,sizeof(int));
    Fspart.read((char*)&nprocs,sizeof(int));
    Fsupart.read((char*)&itask,sizeof(int));
    Fsupart.read((char*)&nprocs,sizeof(int));

    Fspart.read((char*)&nsids,sizeof(Int_t));
    Fspart.read((char*)&nsidstot,sizeof(Int_t));//total across all files
    Fsupart.read((char*)&nsuids,sizeof(Int_t));
    Fsupart.read((char*)&nsuidstot,sizeof(Int_t));//total across all files
    }
    }
    else {
    Fpart>>itask>>nprocs;
    Fupart>>itask>>nprocs;
    Fpart>>nids>>nidstot;
    Fupart>>nuids>>nuidstot;
    if (ifieldhalos) {
    Fspart>>itask>>nprocs;
    Fsupart>>itask>>nprocs;
    Fspart>>nsids>>nsidstot;
    Fsupart>>nsuids>>nsuidstot;
    }
    }
    //cout<<infile<<" has "<<nglocal<<" "<<nids<<" "<<nuids<<" | "<<nsglocal<<" "<<nsids<<" "<<nsuids<<" | "<<nbglocal<<" "<<nbids<<" "<<nbuids<<" | "<<nbsglocal<<" "<<nbsids<<" "<<nbsuids<<endl;
    cout<<infile<<" has "<<nglocal<<" "<<nids<<" "<<nuids<<" | "<<nsglocal<<" "<<nsids<<" "<<nsuids<<endl;
    if (nglocal>0) {
        numingroup=new Int_t[nglocal+1];
        offset=new Int_t[nglocal+1];
        uoffset=new Int_t[nglocal+1];
        numingroupbound=new Int_t[nglocal+1];
        idval=new Int_t[nids+nuids+1];

        if(ibinary) {
        Fgroup.read((char*)numingroup,sizeof(Int_t)*nglocal);
        for (Int_t i=0;i<nglocal;i++)
            Halo[i+noffset].Alloc(numingroup[i]);
        Fgroup.read((char*)offset,sizeof(Int_t)*nglocal);
        Fgroup.read((char*)uoffset,sizeof(Int_t)*nglocal);
        }
        else {
        for (Int_t i=0;i<nglocal;i++) {
            Fgroup>>numingroup[i];
            Halo[i+noffset].Alloc(numingroup[i]);
        }
        for (Int_t i=0;i<nglocal;i++) Fgroup>>offset[i];
        for (Int_t i=0;i<nglocal;i++) Fgroup>>uoffset[i];
        }
        for (Int_t i=0;i<nglocal;i++) nbodies+=numingroup[i];
        //set offsets
        for (Int_t i=0;i<nglocal;i++) {
                Halo[i+noffset].noffset=noldhalopartoffset;
                noldhalopartoffset+=numingroup[i];
        }
        ///store offsets for halo for access to particle array
        //now read bound/unbound particle list
        if(ibinary){
        ///\todo must eventually make ids format a standard size like using GADGETIDTYPE
        nidsoff=0;
        Fpart.read((char*)&idval[nidsoff],sizeof(Int_t)*nids);
        nidsoff+=nids;
        Fupart.read((char*)&idval[nidsoff],sizeof(Int_t)*nuids);
        }
        else {
        nidsoff=0;
        for (Int_t i=0;i<nids;i++) Fpart>>idval[i];
        nidsoff+=nids;
        for (Int_t i=0;i<nuids;i++) Fupart>>idval[i+nidsoff];
        }

        for (Int_t i=0;i<nglocal-1;i++) numingroupbound[i]=numingroup[i]-(uoffset[i+1]-uoffset[i]);
        numingroupbound[nglocal-1]=numingroup[nglocal-1]-(nuids-uoffset[nglocal-1]);
        for (Int_t i=0;i<nglocal;i++) {
            nidsoff=offset[i];
            nn=numingroupbound[i];
            for (Int_t k=0;k<nn;k++) Halo[i+noffset].ParticleID[k]=idval[nidsoff+k];
            nn=numingroup[i]-numingroupbound[i];
            nidsoff=nids+uoffset[i];
            for (Int_t k=0;k<nn;k++) Halo[i+noffset].ParticleID[k+numingroupbound[i]]=idval[nidsoff+k];
        }
        delete[] idval;
        delete[] offset;
        delete[] uoffset;
        delete[] numingroupbound;
        delete[] numingroup;
    }
    for (Int_t i=0;i<nglocal+nsglocal;i++) Halo[i+noffset].haloID=i+1+noffset;

    if (ifieldhalos) {
    if (nsglocal>0) {
        //now read substructure data
        numingroup=new Int_t[nsglocal];
        offset=new Int_t[nsglocal];
        uoffset=new Int_t[nsglocal];
        numingroupbound=new Int_t[nsglocal];
        idval=new Int_t[nsids+nsuids];
        if(ibinary) {
        Fsgroup.read((char*)numingroup,sizeof(Int_t)*nsglocal);
        for (Int_t i=0;i<nsglocal;i++)
            Halo[i+nglocal+noffset].Alloc(numingroup[i]);
        Fsgroup.read((char*)offset,sizeof(Int_t)*nsglocal);
        Fsgroup.read((char*)uoffset,sizeof(Int_t)*nsglocal);
        }
        else {
        for (Int_t i=0;i<nsglocal;i++) {
            Fsgroup>>numingroup[i];
            Halo[i+nglocal+noffset].Alloc(numingroup[i]);
        }
        for (Int_t i=0;i<nsglocal;i++) Fsgroup>>offset[i];
        for (Int_t i=0;i<nsglocal;i++) Fsgroup>>uoffset[i];
        }

        for (Int_t i=0;i<nsglocal;i++) nbodies+=numingroup[i];
        ///store offsets for halo for access to particle array
        for (Int_t i=0;i<nsglocal;i++) {
            Halo[i+noffset+nglocal].noffset=noldhalopartoffset;
            noldhalopartoffset+=numingroup[i];
        }

        if(ibinary) {
        nidsoff=0;
        Fspart.read((char*)&idval[nidsoff],sizeof(Int_t)*nsids);
        nidsoff+=nsids;
        Fsupart.read((char*)&idval[nidsoff],sizeof(Int_t)*nsuids);
        }
        else {
        nidsoff=0;
        for (Int_t i=0;i<nsids;i++) Fspart>>idval[nidsoff+i];
        nidsoff+=nsids;
        for (Int_t i=0;i<nsuids;i++) Fsupart>>idval[nidsoff+i];
        }

        for (Int_t i=0;i<nsglocal-1;i++) numingroupbound[i]=numingroup[i]-(uoffset[i+1]-uoffset[i]);
        numingroupbound[nsglocal-1]=numingroup[nsglocal-1]-(nsuids-uoffset[nsglocal-1]);
        for (Int_t i=0;i<nsglocal;i++) {
            nidsoff=offset[i];
            nn=numingroupbound[i];
            for (Int_t k=0;k<nn;k++) Halo[i+nglocal+noffset].ParticleID[k]=idval[nidsoff+k];
            nn=numingroup[i]-numingroupbound[i];
            nidsoff=nsids+uoffset[i];
            for (Int_t k=0;k<nn;k++) Halo[i+nglocal+noffset].ParticleID[k+numingroupbound[i]]=idval[nidsoff+k];
        }
        delete[] idval;
        delete[] offset;
        delete[] uoffset;
        delete[] numingroupbound;
        delete[] numingroup;
    }
    }

    noffset+=nglocal+nsglocal;

    Fgroup.close();
    Fpart.close();
    Fupart.close();
    if (ifieldhalos) {
        Fsgroup.close();
        Fspart.close();
        Fsupart.close();
    }

    }

    //}//end of ThisTask==0 
    cout<<"Done"<<endl;
    return Halo;
}

///Read halo property data
///\todo update the properties data file structure, particularly the ascii format
PropData *ReadHaloPropertiesData(char* infile, const Int_t numhalos, int mpi_ninput, int ibinary, int ifieldhalos)
{
    cout<<"Read halo properties"<<endl;
    PropData *haloprop;
    Int_t noffset,ngroupoffset,ngroups,nsgroups,ngtot,nfiles,ninfile[2];
    char fname1[1000],fname2[1000];
    fstream Fprop,Fsubprop,*Fp[2];
    int j=0,nmpicount,itask,nprocs;
    //to read correctly sized values
    long long idbound;
    float value,ctemp[3],mtemp[9];
    double dvalue;
    int ivalue;
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif
    ///\todo need to fully implement MPI routines. At the moment, only ThisTask==0 does anything.
    //if (ThisTask==0) {
    haloprop=new PropData[numhalos];
    ///adjust files read based on whether mpi output or not
    if (mpi_ninput==0) nmpicount=1;
    else nmpicount=mpi_ninput;
    //to offset the halo structure pointer given that data maybe split among multiple files
    noffset=0;
    for (j=0;j<nmpicount;j++) {
    ngroups=nsgroups=0;
#ifdef USEMPI
    sprintf(fname1,"%s.properties.%d",infile,j);
    sprintf(fname2,"%s.sublevels.properties.%d",infile,j);
#else
    sprintf(fname1,"%s.properties",infile);
    sprintf(fname2,"%s.sublevels.properties",infile);
#endif
    //write header
    if(ibinary) {
    Fprop.open(fname1,ios::in|ios::binary);
    Fprop.read((char*)&itask,sizeof(int));
    Fprop.read((char*)&nmpicount,sizeof(int));
    Fprop.read((char*)&ngtot,sizeof(Int_t));
    Fprop.read((char*)&ngroups,sizeof(Int_t));

    if(ifieldhalos) {
    Fsubprop.open(fname2,ios::in|ios::binary);
    Fsubprop.read((char*)&itask,sizeof(int));
    Fsubprop.read((char*)&nmpicount,sizeof(int));
    Fsubprop.read((char*)&ngtot,sizeof(Int_t));
    Fsubprop.read((char*)&nsgroups,sizeof(Int_t));
    }
    }
    else {
    Fprop.open(fname1,ios::in);
    Fprop>>itask>>nmpicount;
    Fprop>>ngtot>>ngroups;
    ///read header string from line that dictactes the columns printed 
    Fprop.ignore(10000,'\n');
    if(ifieldhalos) {
    Fsubprop.open(fname2,ios::in);
    Fsubprop>>itask>>nmpicount;
    Fsubprop>>ngtot>>nsgroups;
    Fsubprop.ignore(10000,'\n');
    }
    }

    cout<<"Starting to read halo properties "<<fname1<<" "<<ngroups<<" "<<nsgroups<<endl;
    nfiles=1;Fp[0]=&Fprop;ninfile[0]=ngroups;
    if (ifieldhalos) {nfiles=2;Fp[1]=&Fsubprop;ninfile[1]=nsgroups;}
    ngroupoffset=0;
    for (Int_t ii=0;ii<nfiles;ii++) {
    for (Int_t i=0;i<ninfile[ii];i++) {
        //set halo properties to dm type
        haloprop[i+noffset].ptype=DMTYPE;
        //and read data
        if (ibinary) {
        (*Fp[ii]).read((char*)&idbound,sizeof(long long));
        (*Fp[ii]).read((char*)&dvalue,sizeof(dvalue));
        haloprop[i+noffset+ngroupoffset].gMvir=dvalue;
        (*Fp[ii]).read((char*)&ivalue,sizeof(ivalue));
        haloprop[i+noffset+ngroupoffset].num=ivalue;
        (*Fp[ii]).read((char*)ctemp,sizeof(float)*3);
        for (int k=0;k<3;k++) haloprop[i+noffset+ngroupoffset].gcm[k]=ctemp[k];
        (*Fp[ii]).read((char*)ctemp,sizeof(float)*3);
        for (int k=0;k<3;k++) haloprop[i+noffset+ngroupoffset].gpos[k]=ctemp[k];
        (*Fp[ii]).read((char*)ctemp,sizeof(float)*3);
        for (int k=0;k<3;k++) haloprop[i+noffset+ngroupoffset].gcmvel[k]=ctemp[k];
        (*Fp[ii]).read((char*)ctemp,sizeof(float)*3);
        for (int k=0;k<3;k++) haloprop[i+noffset+ngroupoffset].gvel[k]=ctemp[k];

        (*Fp[ii]).read((char*)&value,sizeof(value));
        haloprop[i+noffset+ngroupoffset].gRvir=value;
        (*Fp[ii]).read((char*)&value,sizeof(value));
        haloprop[i+noffset+ngroupoffset].gRmbp=value;
        (*Fp[ii]).read((char*)&value,sizeof(value));
        haloprop[i+noffset+ngroupoffset].gRmaxvel=value;
        (*Fp[ii]).read((char*)&value,sizeof(value));
        haloprop[i+noffset+ngroupoffset].gmaxvel=value;
        (*Fp[ii]).read((char*)&value,sizeof(value));
        haloprop[i+noffset+ngroupoffset].gsigma_v=value;

        (*Fp[ii]).read((char*)ctemp,sizeof(float)*3);
        for (int k=0;k<3;k++) haloprop[i+noffset+ngroupoffset].gJ[k]=ctemp[k];
        (*Fp[ii]).read((char*)&value,sizeof(value));
        haloprop[i+noffset+ngroupoffset].gq=value;
        (*Fp[ii]).read((char*)&value,sizeof(value));
        haloprop[i+noffset+ngroupoffset].gs=value;
        (*Fp[ii]).read((char*)mtemp,sizeof(float)*9);
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) haloprop[i+noffset+ngroupoffset].geigvec(k,n)=mtemp[k*3+n];
        //afterwards would be padding for 8 more floats (chars), add extra stuff like total mass, mass enclosed in Rmax, T, Pot, Efrac

        (*Fp[ii]).read((char*)&dvalue,sizeof(dvalue));
        haloprop[i+noffset+ngroupoffset].gmass=dvalue;
        (*Fp[ii]).read((char*)&dvalue,sizeof(dvalue));
        haloprop[i+noffset+ngroupoffset].gMmaxvel=dvalue;
        (*Fp[ii]).read((char*)&value,sizeof(value));
        haloprop[i+noffset+ngroupoffset].T=value;
        (*Fp[ii]).read((char*)&value,sizeof(value));
        haloprop[i+noffset+ngroupoffset].Pot=value;
        (*Fp[ii]).read((char*)&value,sizeof(value));
        haloprop[i+noffset+ngroupoffset].Efrac=value;
        for (int k=0;k<1;k++) (*Fp[ii]).read((char*)&value,sizeof(value));//for fact that 2*2*4+3*4
        }
        else {
            ///\todo needs adjusting!!
        (*Fp[ii])>>idbound;
        (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].gMvir;
        (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].num;
        for (int k=0;k<3;k++) (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].gcm[k];
        for (int k=0;k<3;k++) (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].gpos[k];
        for (int k=0;k<3;k++) (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].gcmvel[k];
        for (int k=0;k<3;k++) (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].gvel[k];
        (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].gRvir;
        (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].gRmbp;
        (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].gRmaxvel;
        (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].gmaxvel;
        (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].gsigma_v;
        for (int k=0;k<3;k++) (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].gJ[k];
        (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].gq;
        (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].gs;
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].geigvec(k,n);
        (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].gmass;
        (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].gMmaxvel;
        (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].T;
        (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].Pot;
        (*Fp[ii])>>haloprop[i+noffset+ngroupoffset].Efrac;
        }
    }
    ngroupoffset=ngroups;
    }
        /*
    for (Int_t i=0;i<ngroups;i++) {
        //set halo properties to dm type
        haloprop[i+noffset].ptype=DMTYPE;
        //and read data
        if (ibinary) {
        Fprop.read((char*)&idbound,sizeof(long long));
        Fprop.read((char*)&dvalue,sizeof(dvalue));
        haloprop[i+noffset].gMvir=dvalue;
        Fprop.read((char*)&ivalue,sizeof(ivalue));
        haloprop[i+noffset].num=ivalue;
        Fprop.read((char*)ctemp,sizeof(float)*3);
        for (int k=0;k<3;k++) haloprop[i+noffset].gcm[k]=ctemp[k];
        Fprop.read((char*)ctemp,sizeof(float)*3);
        for (int k=0;k<3;k++) haloprop[i+noffset].gpos[k]=ctemp[k];
        Fprop.read((char*)ctemp,sizeof(float)*3);
        for (int k=0;k<3;k++) haloprop[i+noffset].gcmvel[k]=ctemp[k];
        Fprop.read((char*)ctemp,sizeof(float)*3);
        for (int k=0;k<3;k++) haloprop[i+noffset].gvel[k]=ctemp[k];

        Fprop.read((char*)&value,sizeof(value));
        haloprop[i+noffset].gRvir=value;
        Fprop.read((char*)&value,sizeof(value));
        haloprop[i+noffset].gRmbp=value;
        Fprop.read((char*)&value,sizeof(value));
        haloprop[i+noffset].gRmaxvel=value;
        Fprop.read((char*)&value,sizeof(value));
        haloprop[i+noffset].gmaxvel=value;
        Fprop.read((char*)&value,sizeof(value));
        haloprop[i+noffset].gsigma_v=value;

        Fprop.read((char*)ctemp,sizeof(float)*3);
        for (int k=0;k<3;k++) haloprop[i+noffset].gJ[k]=ctemp[k];
        Fprop.read((char*)&value,sizeof(value));
        haloprop[i+noffset].gq=value;
        Fprop.read((char*)&value,sizeof(value));
        haloprop[i+noffset].gs=value;
        Fprop.read((char*)mtemp,sizeof(float)*9);
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) haloprop[i+noffset].geigvec(k,n)=mtemp[k*3+n];
        //afterwards would be padding for 8 more floats (chars), add extra stuff like total mass, mass enclosed in Rmax, T, Pot, Efrac

        Fprop.read((char*)&dvalue,sizeof(dvalue));
        haloprop[i+noffset].gmass=dvalue;
        Fprop.read((char*)&dvalue,sizeof(dvalue));
        haloprop[i+noffset].gMmaxvel=dvalue;
        Fprop.read((char*)&value,sizeof(value));
        haloprop[i+noffset].T=value;
        Fprop.read((char*)&value,sizeof(value));
        haloprop[i+noffset].Pot=value;
        Fprop.read((char*)&value,sizeof(value));
        haloprop[i+noffset].Efrac=value;
        for (int k=0;k<1;k++) Fprop.read((char*)&value,sizeof(value));//for fact that 2*2*4+3*4
        }
        else {
            ///\todo needs adjusting!!
        Fprop>>idbound;
        Fprop>>haloprop[i+noffset].gMvir;
        Fprop>>haloprop[i+noffset].num;
        for (int k=0;k<3;k++) Fprop>>haloprop[i+noffset].gcm[k];
        for (int k=0;k<3;k++) Fprop>>haloprop[i+noffset].gpos[k];
        for (int k=0;k<3;k++) Fprop>>haloprop[i+noffset].gcmvel[k];
        for (int k=0;k<3;k++) Fprop>>haloprop[i+noffset].gvel[k];
        Fprop>>haloprop[i+noffset].gRvir;
        Fprop>>haloprop[i+noffset].gRmbp;
        Fprop>>haloprop[i+noffset].gRmaxvel;
        Fprop>>haloprop[i+noffset].gmaxvel;
        Fprop>>haloprop[i+noffset].gsigma_v;
        for (int k=0;k<3;k++) Fprop>>haloprop[i+noffset].gJ[k];
        Fprop>>haloprop[i+noffset].gq;
        Fprop>>haloprop[i+noffset].gs;
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fprop>>haloprop[i+noffset].geigvec(k,n);
        Fprop>>haloprop[i+noffset].gmass;
        Fprop>>haloprop[i+noffset].gMmaxvel;
        Fprop>>haloprop[i+noffset].T;
        Fprop>>haloprop[i+noffset].Pot;
        Fprop>>haloprop[i+noffset].Efrac;
        }

        ///\todo need to make error check for reference frame more robust
        //if gcmvel nan it is because poor interation in stf. for quick fix, set gcmvel to gvel if nan
        if (isnan(haloprop[i+noffset].gcmvel[0])) for (int k=0;k<3;k++) haloprop[i+noffset].gcmvel[k]=haloprop[i+noffset].gvel[k];
    }
    Fprop.close();
    cout<<"Properties read for halos"<<endl;
    if (ifieldhalos) {
    for (Int_t i=0;i<nsgroups;i++) {
        if (ibinary) {
        Fsubprop.read((char*)&idbound,sizeof(long long));
        Fsubprop.read((char*)&dvalue,sizeof(dvalue));
        haloprop[i+ngroups+noffset].gMvir=dvalue;
        Fsubprop.read((char*)&ivalue,sizeof(ivalue));
        haloprop[i+ngroups+noffset].num=ivalue;
        Fsubprop.read((char*)ctemp,sizeof(float)*3);
        for (int k=0;k<3;k++) haloprop[i+ngroups+noffset].gcm[k]=ctemp[k];
        Fsubprop.read((char*)ctemp,sizeof(float)*3);
        for (int k=0;k<3;k++) haloprop[i+ngroups+noffset].gpos[k]=ctemp[k];
        Fsubprop.read((char*)ctemp,sizeof(float)*3);
        for (int k=0;k<3;k++) haloprop[i+ngroups+noffset].gcmvel[k]=ctemp[k];
        Fsubprop.read((char*)ctemp,sizeof(float)*3);
        for (int k=0;k<3;k++) haloprop[i+ngroups+noffset].gvel[k]=ctemp[k];

        Fsubprop.read((char*)&value,sizeof(value));
        haloprop[i+ngroups+noffset].gRvir=value;
        Fsubprop.read((char*)&value,sizeof(value));
        haloprop[i+ngroups+noffset].gRmbp=value;
        Fsubprop.read((char*)&value,sizeof(value));
        haloprop[i+ngroups+noffset].gRmaxvel=value;
        Fsubprop.read((char*)&value,sizeof(value));
        haloprop[i+ngroups+noffset].gmaxvel=value;
        Fsubprop.read((char*)&value,sizeof(value));
        haloprop[i+ngroups+noffset].gsigma_v=value;

        Fsubprop.read((char*)ctemp,sizeof(float)*3);
        for (int k=0;k<3;k++) haloprop[i+ngroups+noffset].gJ[k]=ctemp[k];
        Fsubprop.read((char*)&value,sizeof(value));
        haloprop[i+ngroups+noffset].gq=value;
        Fsubprop.read((char*)&value,sizeof(value));
        haloprop[i+ngroups+noffset].gs=value;
        Fsubprop.read((char*)mtemp,sizeof(float)*9);
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) haloprop[i+ngroups+noffset].geigvec(k,n)=mtemp[k*3+n];
        //afterwards would be padding for 8 more floats (chars), add extra stuff like total mass, mass enclosed in Rmax, T, Pot, Efrac

        Fsubprop.read((char*)&dvalue,sizeof(dvalue));
        haloprop[i+ngroups+noffset].gmass=dvalue;
        Fsubprop.read((char*)&dvalue,sizeof(dvalue));
        haloprop[i+ngroups+noffset].gMmaxvel=dvalue;

        Fsubprop.read((char*)&value,sizeof(value));
        haloprop[i+ngroups+noffset].T=value;
        Fsubprop.read((char*)&value,sizeof(value));
        haloprop[i+ngroups+noffset].Pot=value;
        Fsubprop.read((char*)&value,sizeof(value));
        haloprop[i+ngroups+noffset].Efrac=value;
        for (int k=0;k<1;k++) Fsubprop.read((char*)&value,sizeof(value));//for fact that 2*2*4+3*4
        }
        else {
        Fsubprop>>idbound;
        Fsubprop>>haloprop[i+ngroups+noffset].gMvir;
        Fsubprop>>haloprop[i+ngroups+noffset].num;
        for (int k=0;k<3;k++) Fsubprop>>haloprop[i+ngroups+noffset].gcm[k];
        for (int k=0;k<3;k++) Fsubprop>>haloprop[i+ngroups+noffset].gpos[k];
        for (int k=0;k<3;k++) Fsubprop>>haloprop[i+ngroups+noffset].gcmvel[k];
        for (int k=0;k<3;k++) Fsubprop>>haloprop[i+ngroups+noffset].gvel[k];
        Fsubprop>>haloprop[i+ngroups+noffset].gRvir;
        Fsubprop>>haloprop[i+ngroups+noffset].gRmbp;
        Fsubprop>>haloprop[i+ngroups+noffset].gRmaxvel;
        Fsubprop>>haloprop[i+ngroups+noffset].gmaxvel;
        Fsubprop>>haloprop[i+ngroups+noffset].gsigma_v;
        for (int k=0;k<3;k++) Fsubprop>>haloprop[i+ngroups+noffset].gJ[k];
        Fsubprop>>haloprop[i+ngroups+noffset].gq;
        Fsubprop>>haloprop[i+ngroups+noffset].gs;
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fsubprop>>haloprop[i+ngroups+noffset].geigvec(k,n);
        Fsubprop>>haloprop[i+ngroups+noffset].gmass;
        Fsubprop>>haloprop[i+ngroups+noffset].gMmaxvel;
        Fsubprop>>haloprop[i+ngroups+noffset].T;
        Fsubprop>>haloprop[i+ngroups+noffset].Pot;
        Fsubprop>>haloprop[i+ngroups+noffset].Efrac;
        }
        ///\todo need to make error check for reference frame more robust
        //if gcmvel nan it is because poor interation in stf. for quick fix, set gcmvel to gvel if nan
        if (isnan(haloprop[i+ngroups+noffset].gcmvel[0])) for (int k=0;k<3;k++) haloprop[i+ngroups+noffset].gcmvel[k]=haloprop[i+ngroups+noffset].gvel[k];
    }
    Fsubprop.close();
    cout<<"Properties read for subhalos"<<endl;
    }
    */

    noffset+=ngroups+nsgroups;
    }
    cout<<"Done"<<endl;
    //}//ThisTask==0
    return haloprop;
}

//@}
