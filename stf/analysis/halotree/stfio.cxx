/*! \file stfio.cxx
 *  \brief this file contains routines for reading stf style group catalog data
 */

#include "halomergertree.h"
int CheckType(unsigned int t, int tmatch){
    if (tmatch>ALLTYPEMATCH) return (t==tmatch);
    else if (tmatch==DMGASTYPEMATCH) return (t==GASTYPEMATCH || t==DMTYPEMATCH);
    else return 0;
}

#ifdef USEMPI
///Read halo data from an idividual snapshot;
Int_t MPIReadHaloGroupCatalogDataNum(char* infile, int mpi_ninput, int ibinary, int ifieldhalos, int itypematch)
{
    Int_t itemp;
    unsigned long ltemp;
    Int_t noffset,numfiletypes;
    long unsigned nparts,haloid;
    long unsigned TotalNumberofHalos;
    char fname1[1000],fname2[1000],fname3[1000];
    char fname4[1000],fname5[1000],fname6[1000];
    char fname7[1000],fname8[1000];
    char fname9[1000],fname10[1000];
    fstream Fgroup,Fpart,Fupart; //field objects
    fstream Fsgroup,Fspart,Fsupart; //sublevels
    fstream Fparttype,Fuparttype; //field objects
    fstream Fsparttype,Fsuparttype; //sublevels
    char *fnamearray[10];
    fstream *Farray[10];
    unsigned long nids,nuids,nsids,nsuids,nglocal,nsglocal;
    unsigned long nidstot,nuidstot,nsidstot,nsuidstot;
    Int_t *numingroup,*numingroupbound,*offset,*uoffset;
    Int_t counter,nn;
    int nmpicount,itask,nprocs;
    Int_t *idval;
    UInt_t *typeval;
    numfiletypes=3;
    if (ifieldhalos) numfiletypes+=3;
    if (itypematch!=ALLTYPEMATCH) {
        numfiletypes+=2;
        if (ifieldhalos) numfiletypes+=2;
    }

    ///adjust files read based on whether mpi output or not
    if (mpi_ninput==0) nmpicount=1;
    else nmpicount=mpi_ninput;
    //to offset the halo structure pointer given that data maybe split among multiple files
    noffset=0;

    nglocal=nsglocal=0;
    //load the group catalogues and the related particle files
    if (mpi_ninput==0) {
        sprintf(fname1,"%s.catalog_groups",infile);
        sprintf(fname2,"%s.catalog_particles",infile);
        sprintf(fname3,"%s.catalog_particles.unbound",infile);
        sprintf(fname4,"%s.sublevels.catalog_groups",infile);
        sprintf(fname5,"%s.sublevels.catalog_particles",infile);
        sprintf(fname6,"%s.sublevels.catalog_particles.unbound",infile);
        sprintf(fname7,"%s.catalog_parttypes",infile);
        sprintf(fname8,"%s.catalog_parttypes.unbound",infile);
        sprintf(fname9,"%s.sublevels.catalog_parttypes",infile);
        sprintf(fname10,"%s.sublevels.catalog_parttypes.unbound",infile);
    }
    else {
        sprintf(fname1,"%s.catalog_groups.%d",infile,0);
        sprintf(fname2,"%s.catalog_particles.%d",infile,0);
        sprintf(fname3,"%s.catalog_particles.unbound.%d",infile,0);
        sprintf(fname4,"%s.sublevels.catalog_groups.%d",infile,0);
        sprintf(fname5,"%s.sublevels.catalog_particles.%d",infile,0);
        sprintf(fname6,"%s.sublevels.catalog_particles.unbound.%d",infile,0);
        sprintf(fname7,"%s.catalog_parttypes.%d",infile,0);
        sprintf(fname8,"%s.catalog_parttypes.unbound.%d",infile,0);
        sprintf(fname9,"%s.sublevels.catalog_parttypes.%d",infile,0);
        sprintf(fname10,"%s.sublevels.catalog_parttypes.unbound.%d",infile,0);
    }
    itemp=0;
    fnamearray[itemp++]=fname1;fnamearray[itemp++]=fname2;fnamearray[itemp++]=fname3;
    if (ifieldhalos) {
        fnamearray[itemp++]=fname4;fnamearray[itemp++]=fname5;fnamearray[itemp++]=fname6;
    }
    if (itypematch!=ALLTYPEMATCH) { 
        fnamearray[itemp++]=fname7;fnamearray[itemp++]=fname8;
        if (ifieldhalos) {
            fnamearray[itemp++]=fname9;fnamearray[itemp++]=fname10;
        }
    }
    if (ibinary) {
        Fgroup.open(fname1,ios::in|ios::binary);
        Fpart.open(fname2,ios::in|ios::binary);
        Fupart.open(fname3,ios::in|ios::binary);
        if (ifieldhalos) {
            Fsgroup.open(fname4,ios::in|ios::binary);
            Fspart.open(fname5,ios::in|ios::binary);
            Fsupart.open(fname6,ios::in|ios::binary);
        }
        if (itypematch!=ALLTYPEMATCH) {
            Fparttype.open(fname7,ios::in|ios::binary);
            Fuparttype.open(fname8,ios::in|ios::binary);
            if (ifieldhalos) {
                Fsparttype.open(fname9,ios::in|ios::binary);
                Fsuparttype.open(fname10,ios::in|ios::binary);
            }
        }
    }
    else {
        Fgroup.open(fname1,ios::in);
        Fpart.open(fname2,ios::in);
        Fupart.open(fname3,ios::in);
        if (ifieldhalos) {
            Fsgroup.open(fname4,ios::in);
            Fspart.open(fname5,ios::in);
            Fsupart.open(fname6,ios::in);
        }
        if (itypematch!=ALLTYPEMATCH) {
            Fparttype.open(fname7,ios::in);
            Fuparttype.open(fname8,ios::in);
            if (ifieldhalos) {
                Fsparttype.open(fname9,ios::in);
                Fsuparttype.open(fname10,ios::in);
            }
        }
    }
    itemp=0;
    Farray[itemp++]=&Fgroup;Farray[itemp++]=&Fpart;Farray[itemp++]=&Fupart;
    if (ifieldhalos) {
        Farray[itemp++]=&Fsgroup;Farray[itemp++]=&Fspart;Farray[itemp++]=&Fsupart;
    }
    if (itypematch!=ALLTYPEMATCH) {
        Farray[itemp++]=&Fparttype;Farray[itemp++]=&Fuparttype;
        if (ifieldhalos) {
            Farray[itemp++]=&Fsparttype;Farray[itemp++]=&Fsuparttype;
        }
    }
    for (int i=0;i<numfiletypes;i++) {
        if(!Farray[i]->is_open()){
            cerr<<"can't open "<<fnamearray[i]<<endl;
            MPI_Abort(MPI_COMM_WORLD,9);
        }
    }
    //check header make sure that number of mpi values considered agree
    if (ibinary) {
        Fgroup.read((char*)&itask,sizeof(int));
        Fgroup.read((char*)&nprocs,sizeof(int));
        if (ifieldhalos) {
            Fsgroup.read((char*)&itask,sizeof(int));
            Fsgroup.read((char*)&nprocs,sizeof(int));
        }
    }
    else {
        Fgroup>>itask>>nprocs;
        if (ifieldhalos) {
            Fsgroup>>itask>>nprocs;
        }
    }
    if (nprocs!=mpi_ninput&&mpi_ninput>0) {
        cout<<"Error, number of mpi outputs was set to "<<mpi_ninput<<" but file indicates there are "<<nprocs<<endl;
        cout<<"Correcting to this number and proceeding"<<endl;
        nmpicount=mpi_ninput=nprocs;
    }
    else if ((mpi_ninput==0&&nprocs!=1)) {
        cout<<"Error, number of mpi outputs was set to zero but file indicates there are more than one mpi output"<<endl;
        cout<<"Terminating"<<endl;
        exit(9);
    }

    if (ibinary) {
        Fgroup.read((char*)&ltemp,sizeof(unsigned long));
        nglocal=ltemp;
        Fgroup.read((char*)&ltemp,sizeof(unsigned long));
        TotalNumberofHalos=ltemp;
        if (ifieldhalos) {
            Fsgroup.read((char*)&ltemp,sizeof((unsigned long));
            nsglocal=ltemp;
            Fsgroup.read((char*)&ltemp,sizeof((unsigned long));
        }
    }
    else {
        Fgroup>>nglocal;
        Fgroup>>TotalNumberofHalos;
        if (ifieldhalos) {
            Fsgroup>>nsglocal;
            Fsgroup>>itemp;
        }
    }
    for (int i=0;i<numfiletypes;i++) Farray[i]->close();
    
    return TotalNumberofHalos;
}

///Read the number of particles in halos from an idividual snapshot, useful for splitting snaphots across threads in load balanced fashion
Int_t MPIReadHaloGroupCatalogDataParticleNum(char* infile, int mpi_ninput, int ibinary, int ifieldhalos, int itypematch)
{
    
    HaloData *Halo;
    Int_t itemp;
    unsigned long ltemp;
    Int_t noffset,numfiletypes;
    long unsigned nparts,haloid;
    long unsigned TotalNumberofParticles=0;
    char fname1[1000],fname2[1000],fname3[1000];
    char fname4[1000],fname5[1000],fname6[1000];
    char fname7[1000],fname8[1000];
    char fname9[1000],fname10[1000];
    fstream Fgroup,Fpart,Fupart; //field objects
    fstream Fsgroup,Fspart,Fsupart; //sublevels
    fstream Fparttype,Fuparttype; //field objects
    fstream Fsparttype,Fsuparttype; //sublevels
    char *fnamearray[10];
    fstream *Farray[10];
    unsigned long nids,nuids,nsids,nsuids,nglocal,nsglocal;
    unsigned long nidstot,nuidstot,nsidstot,nsuidstot;
    Int_t *numingroup,*numingroupbound,*offset,*uoffset;
    Int_t counter,nn;
    int nmpicount,itask,nprocs;
    Int_t *idval;
    UInt_t *typeval;
    numfiletypes=3;
    if (ifieldhalos) numfiletypes+=3;
    if (itypematch!=ALLTYPEMATCH) {
        numfiletypes+=2;
        if (ifieldhalos) numfiletypes+=2;
    }

    ///adjust files read based on whether mpi output or not
    if (mpi_ninput==0) nmpicount=1;
    else nmpicount=mpi_ninput;
    //to offset the halo structure pointer given that data maybe split among multiple files
    noffset=0;
    for (int k=0;k<nmpicount;k++) {
        nglocal=nsglocal=0;
        //load the group catalogues and the related particle files
        if (mpi_ninput==0) {
            sprintf(fname1,"%s.catalog_groups",infile);
            sprintf(fname2,"%s.catalog_particles",infile);
            sprintf(fname3,"%s.catalog_particles.unbound",infile);
            sprintf(fname4,"%s.sublevels.catalog_groups",infile);
            sprintf(fname5,"%s.sublevels.catalog_particles",infile);
            sprintf(fname6,"%s.sublevels.catalog_particles.unbound",infile);
            sprintf(fname7,"%s.catalog_parttypes",infile);
            sprintf(fname8,"%s.catalog_parttypes.unbound",infile);
            sprintf(fname9,"%s.sublevels.catalog_parttypes",infile);
            sprintf(fname10,"%s.sublevels.catalog_parttypes.unbound",infile);
        }
        else {
            sprintf(fname1,"%s.catalog_groups.%d",infile,k);
            sprintf(fname2,"%s.catalog_particles.%d",infile,k);
            sprintf(fname3,"%s.catalog_particles.unbound.%d",infile,k);
            sprintf(fname4,"%s.sublevels.catalog_groups.%d",infile,k);
            sprintf(fname5,"%s.sublevels.catalog_particles.%d",infile,k);
            sprintf(fname6,"%s.sublevels.catalog_particles.unbound.%d",infile,k);
            sprintf(fname7,"%s.catalog_parttypes.%d",infile,k);
            sprintf(fname8,"%s.catalog_parttypes.unbound.%d",infile,k);
            sprintf(fname9,"%s.sublevels.catalog_parttypes.%d",infile,k);
            sprintf(fname10,"%s.sublevels.catalog_parttypes.unbound.%d",infile,k);
        }
        itemp=0;
        fnamearray[itemp++]=fname1;fnamearray[itemp++]=fname2;fnamearray[itemp++]=fname3;
        if (ifieldhalos) {
            fnamearray[itemp++]=fname4;fnamearray[itemp++]=fname5;fnamearray[itemp++]=fname6;
        }
        if (itypematch!=ALLTYPEMATCH) { 
            fnamearray[itemp++]=fname7;fnamearray[itemp++]=fname8;
            if (ifieldhalos) {
                fnamearray[itemp++]=fname9;fnamearray[itemp++]=fname10;
            }
        }
        if (ibinary) {
            Fgroup.open(fname1,ios::in|ios::binary);
            Fpart.open(fname2,ios::in|ios::binary);
            Fupart.open(fname3,ios::in|ios::binary);
            if (ifieldhalos) {
                Fsgroup.open(fname4,ios::in|ios::binary);
                Fspart.open(fname5,ios::in|ios::binary);
                Fsupart.open(fname6,ios::in|ios::binary);
            }
            if (itypematch!=ALLTYPEMATCH) {
                Fparttype.open(fname7,ios::in|ios::binary);
                Fuparttype.open(fname8,ios::in|ios::binary);
                if (ifieldhalos) {
                    Fsparttype.open(fname9,ios::in|ios::binary);
                    Fsuparttype.open(fname10,ios::in|ios::binary);
                }
            }
        }
        else {
            Fgroup.open(fname1,ios::in);
            Fpart.open(fname2,ios::in);
            Fupart.open(fname3,ios::in);
            if (ifieldhalos) {
                Fsgroup.open(fname4,ios::in);
                Fspart.open(fname5,ios::in);
                Fsupart.open(fname6,ios::in);
            }
            if (itypematch!=ALLTYPEMATCH) {
                Fparttype.open(fname7,ios::in);
                Fuparttype.open(fname8,ios::in);
                if (ifieldhalos) {
                    Fsparttype.open(fname9,ios::in);
                    Fsuparttype.open(fname10,ios::in);
                }
            }
        }
        itemp=0;
        Farray[itemp++]=&Fgroup;Farray[itemp++]=&Fpart;Farray[itemp++]=&Fupart;
        if (ifieldhalos) {
            Farray[itemp++]=&Fsgroup;Farray[itemp++]=&Fspart;Farray[itemp++]=&Fsupart;
        }
        if (itypematch!=ALLTYPEMATCH) {
            Farray[itemp++]=&Fparttype;Farray[itemp++]=&Fuparttype;
            if (ifieldhalos) {
                Farray[itemp++]=&Fsparttype;Farray[itemp++]=&Fsuparttype;
            }
        }
        for (int i=0;i<numfiletypes;i++) {
            if(!Farray[i]->is_open()){
                cerr<<"can't open "<<fnamearray[i]<<endl;
                MPI_Abort(MPI_COMM_WORLD,9);
            }
        }
        //check header make sure that number of mpi values considered agree
        if (ibinary) {
            Fgroup.read((char*)&itask,sizeof(int));
            Fgroup.read((char*)&nprocs,sizeof(int));
            if (ifieldhalos) {
                Fsgroup.read((char*)&itask,sizeof(int));
                Fsgroup.read((char*)&nprocs,sizeof(int));
            }
        }
        else {
            Fgroup>>itask>>nprocs;
            if (ifieldhalos) {
                Fsgroup>>itask>>nprocs;
            }
        }

        if (ibinary) {
            Fgroup.read((char*)&ltemp,sizeof(unsigned long));
            nglocal=ltemp;
            Fgroup.read((char*)&ltemp,sizeof(unsigned long));
            if (ifieldhalos) {
                Fsgroup.read((char*)&ltemp,sizeof(unsigned long));
                nsglocal=ltemp;
                Fsgroup.read((char*)&ltemp,sizeof(unsigned long));
            }
        }
        else {
            Fgroup>>nglocal;
            Fgroup>>itemp;
            if (ifieldhalos) {
                Fsgroup>>nsglocal;
                Fsgroup>>itemp;
            }
        }
        if (nglocal>0) {
            numingroup=new Int_t[nglocal+1];

            if (ibinary) {
                Fgroup.read((char*)numingroup,sizeof(Int_t)*nglocal);
                for (Int_t i=0;i<nglocal;i++) TotalNumberofParticles+=numingroup[i];
            }
            else {
                for (Int_t i=0;i<nglocal;i++) {Fgroup>>numingroup[i];TotalNumberofParticles+=numingroup[i];}
            }
            delete[] numingroup;
        }

        if (ifieldhalos) {
            if (nsglocal>0) {
                //now read substructure data
                numingroup=new Int_t[nsglocal];
                if (ibinary) {
                    Fsgroup.read((char*)numingroup,sizeof(Int_t)*nsglocal);
                    for (Int_t i=0;i<nsglocal;i++) TotalNumberofParticles+=numingroup[i];
                }
                else {
                    for (Int_t i=0;i<nsglocal;i++) {Fsgroup>>numingroup[i];TotalNumberofParticles+=numingroup[i];}
                }
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
        if (itypematch!=ALLTYPEMATCH) {
            Fparttype.close();
            Fuparttype.close();
            if (ifieldhalos) {
                Fsparttype.close();
                Fsuparttype.close();
            }
        }

    }

    return TotalNumberofParticles;
}

///Allocate memory to store halo information
HaloData *MPIReadHaloGroupCatalogDataAllocation(char* infile, Int_t &numhalos, int mpi_ninput, int ibinary, int ifieldhalos, int itypematch)
{
    HaloData *Halo;
    Int_t itemp;
    unsigned long ltemp;
    Int_t noffset,numfiletypes;
    long unsigned nparts,haloid;
    long unsigned TotalNumberofHalos;
    char fname1[1000],fname2[1000],fname3[1000];
    char fname4[1000],fname5[1000],fname6[1000];
    char fname7[1000],fname8[1000];
    char fname9[1000],fname10[1000];
    fstream Fgroup,Fpart,Fupart; //field objects
    fstream Fsgroup,Fspart,Fsupart; //sublevels
    fstream Fparttype,Fuparttype; //field objects
    fstream Fsparttype,Fsuparttype; //sublevels
    char *fnamearray[10];
    fstream *Farray[10];
    unsigned long nids,nuids,nsids,nsuids,nglocal,nsglocal;
    unsigned long nidstot,nuidstot,nsidstot,nsuidstot;
    Int_t *numingroup,*numingroupbound,*offset,*uoffset;
    Int_t counter,nn;
    int nmpicount,itask,nprocs;
    Int_t *idval;
    UInt_t *typeval;
    numfiletypes=3;
    if (ifieldhalos) numfiletypes+=3;
    if (itypematch!=ALLTYPEMATCH) {
        numfiletypes+=2;
        if (ifieldhalos) numfiletypes+=2;
    }

    ///adjust files read based on whether mpi output or not
    if (mpi_ninput==0) nmpicount=1;
    else nmpicount=mpi_ninput;
    //to offset the halo structure pointer given that data maybe split among multiple files
    noffset=0;
    for (int k=0;k<nmpicount;k++) {
    nglocal=nsglocal=0;
    //load the group catalogues and the related particle files
    if (mpi_ninput==0) {
        sprintf(fname1,"%s.catalog_groups",infile);
        sprintf(fname2,"%s.catalog_particles",infile);
        sprintf(fname3,"%s.catalog_particles.unbound",infile);
        sprintf(fname4,"%s.sublevels.catalog_groups",infile);
        sprintf(fname5,"%s.sublevels.catalog_particles",infile);
        sprintf(fname6,"%s.sublevels.catalog_particles.unbound",infile);
        sprintf(fname7,"%s.catalog_parttypes",infile);
        sprintf(fname8,"%s.catalog_parttypes.unbound",infile);
        sprintf(fname9,"%s.sublevels.catalog_parttypes",infile);
        sprintf(fname10,"%s.sublevels.catalog_parttypes.unbound",infile);
    }
    else {
        sprintf(fname1,"%s.catalog_groups.%d",infile,k);
        sprintf(fname2,"%s.catalog_particles.%d",infile,k);
        sprintf(fname3,"%s.catalog_particles.unbound.%d",infile,k);
        sprintf(fname4,"%s.sublevels.catalog_groups.%d",infile,k);
        sprintf(fname5,"%s.sublevels.catalog_particles.%d",infile,k);
        sprintf(fname6,"%s.sublevels.catalog_particles.unbound.%d",infile,k);
        sprintf(fname7,"%s.catalog_parttypes.%d",infile,k);
        sprintf(fname8,"%s.catalog_parttypes.unbound.%d",infile,k);
        sprintf(fname9,"%s.sublevels.catalog_parttypes.%d",infile,k);
        sprintf(fname10,"%s.sublevels.catalog_parttypes.unbound.%d",infile,k);
    }
    itemp=0;
    fnamearray[itemp++]=fname1;fnamearray[itemp++]=fname2;fnamearray[itemp++]=fname3;
    if (ifieldhalos) {
        fnamearray[itemp++]=fname4;fnamearray[itemp++]=fname5;fnamearray[itemp++]=fname6;
    }
    if (itypematch!=ALLTYPEMATCH) { 
        fnamearray[itemp++]=fname7;fnamearray[itemp++]=fname8;
        if (ifieldhalos) {
            fnamearray[itemp++]=fname9;fnamearray[itemp++]=fname10;
        }
    }
    if (ibinary) {
        Fgroup.open(fname1,ios::in|ios::binary);
        Fpart.open(fname2,ios::in|ios::binary);
        Fupart.open(fname3,ios::in|ios::binary);
        if (ifieldhalos) {
        Fsgroup.open(fname4,ios::in|ios::binary);
        Fspart.open(fname5,ios::in|ios::binary);
        Fsupart.open(fname6,ios::in|ios::binary);
        }
        if (itypematch!=ALLTYPEMATCH) {
            Fparttype.open(fname7,ios::in|ios::binary);
            Fuparttype.open(fname8,ios::in|ios::binary);
            if (ifieldhalos) {
                Fsparttype.open(fname9,ios::in|ios::binary);
                Fsuparttype.open(fname10,ios::in|ios::binary);
            }
        }
    }
    else {
        Fgroup.open(fname1,ios::in);
        Fpart.open(fname2,ios::in);
        Fupart.open(fname3,ios::in);
        if (ifieldhalos) {
        Fsgroup.open(fname4,ios::in);
        Fspart.open(fname5,ios::in);
        Fsupart.open(fname6,ios::in);
        }
        if (itypematch!=ALLTYPEMATCH) {
            Fparttype.open(fname7,ios::in);
            Fuparttype.open(fname8,ios::in);
            if (ifieldhalos) {
                Fsparttype.open(fname9,ios::in);
                Fsuparttype.open(fname10,ios::in);
            }
        }
    }
    itemp=0;
    Farray[itemp++]=&Fgroup;Farray[itemp++]=&Fpart;Farray[itemp++]=&Fupart;
    if (ifieldhalos) {
        Farray[itemp++]=&Fsgroup;Farray[itemp++]=&Fspart;Farray[itemp++]=&Fsupart;
    }
    if (itypematch!=ALLTYPEMATCH) {
        Farray[itemp++]=&Fparttype;Farray[itemp++]=&Fuparttype;
        if (ifieldhalos) {
            Farray[itemp++]=&Fsparttype;Farray[itemp++]=&Fsuparttype;
        }
    }
    for (int i=0;i<numfiletypes;i++) {
        if(!Farray[i]->is_open()){
            cerr<<"can't open "<<fnamearray[i]<<endl;
            MPI_Abort(MPI_COMM_WORLD,9);
        }
    }
    //check header make sure that number of mpi values considered agree
    if (ibinary) {
        Fgroup.read((char*)&itask,sizeof(int));
        Fgroup.read((char*)&nprocs,sizeof(int));
        if (ifieldhalos) {
            Fsgroup.read((char*)&itask,sizeof(int));
            Fsgroup.read((char*)&nprocs,sizeof(int));
        }
    }
        else {
        Fgroup>>itask>>nprocs;
        if (ifieldhalos) {
        Fsgroup>>itask>>nprocs;
        }
    }

    if (ibinary) {
        Fgroup.read((char*)&ltemp,sizeof(unsigned long));
        nglocal=ltemp;
        Fgroup.read((char*)&ltemp,sizeof(unsigned long));
        TotalNumberofHalos=ltemp;
        if (ifieldhalos) {
            Fsgroup.read((char*)&ltemp,sizeof(unsigned long));
            nsglocal=ltemp;
            Fsgroup.read((char*)&ltemp,sizeof(unsigned long));
        }
    }
    else {
        Fgroup>>nglocal;
        Fgroup>>TotalNumberofHalos;
        if (ifieldhalos) {
            Fsgroup>>nsglocal;
            Fsgroup>>itemp;
        }
    }

    //allocate memory for halos
    if (k==0) {
        Halo=new HaloData[TotalNumberofHalos];
        numhalos=TotalNumberofHalos;
    }
    if (nglocal>0) {
        numingroup=new Int_t[nglocal+1];

        if (ibinary) {
            Fgroup.read((char*)numingroup,sizeof(Int_t)*nglocal);
            for (Int_t i=0;i<nglocal;i++) {
                Halo[i+noffset].Alloc(numingroup[i]);
            }
        }
        else {
            for (Int_t i=0;i<nglocal;i++) {
                Fgroup>>numingroup[i];
                Halo[i+noffset].Alloc(numingroup[i]);
            }
        }

        delete[] numingroup;
    }

    if (ifieldhalos) {
        if (nsglocal>0) {
            //now read substructure data
            numingroup=new Int_t[nsglocal];
            if (ibinary) {
                Fsgroup.read((char*)numingroup,sizeof(Int_t)*nsglocal);
                for (Int_t i=0;i<nsglocal;i++) {
                    Halo[i+nglocal+noffset].Alloc(numingroup[i]);
                }
            }
            else {
                for (Int_t i=0;i<nsglocal;i++) {
                    Fsgroup>>numingroup[i];
                    Halo[i+nglocal+noffset].Alloc(numingroup[i]);
                }
            }
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
    if (itypematch!=ALLTYPEMATCH) {
        Fparttype.close();
        Fuparttype.close();
        if (ifieldhalos) {
            Fsparttype.close();
            Fsuparttype.close();
        }
    }

    }

    return Halo;
}

///Read halo data from an individual snapshot;
void MPIReadHaloGroupCatalogData(char* infile, Int_t &numhalos, HaloData *&Halo, int mpi_ninput, int ibinary, int ifieldhalos, int itypematch, int iverbose)
{
    Int_t itemp;
    unsigned long ltemp;
    Int_t noffset,numfiletypes;
    long unsigned nparts,haloid;
    long unsigned TotalNumberofHalos;
    char fname1[1000],fname2[1000],fname3[1000];
    char fname4[1000],fname5[1000],fname6[1000];
    char fname7[1000],fname8[1000];
    char fname9[1000],fname10[1000];
    fstream Fgroup,Fpart,Fupart; //field objects
    fstream Fsgroup,Fspart,Fsupart; //sublevels
    fstream Fparttype,Fuparttype; //field objects
    fstream Fsparttype,Fsuparttype; //sublevels
    char *fnamearray[10];
    fstream *Farray[10];
    unsigned long nids,nuids,nsids,nsuids,nglocal,nsglocal;
    unsigned long nidstot,nuidstot,nsidstot,nsuidstot;
    Int_t *numingroup,*numingroupbound,*offset,*uoffset;
    Int_t counter,nn;
    int nmpicount,itask,nprocs;
    Int_t *idval;
    UInt_t *typeval;
    numfiletypes=3;
    if (ifieldhalos) numfiletypes+=3;
    if (itypematch!=ALLTYPEMATCH) {
        numfiletypes+=2;
        if (ifieldhalos) numfiletypes+=2;
    }

    ///adjust files read based on whether mpi output or not
    if (mpi_ninput==0) nmpicount=1;
    else nmpicount=mpi_ninput;
    //to offset the halo structure pointer given that data maybe split among multiple files
    noffset=0;
    for (int k=0;k<nmpicount;k++) {
    nglocal=nsglocal=0;
    //load the group catalogues and the related particle files
    if (mpi_ninput==0) {
        sprintf(fname1,"%s.catalog_groups",infile);
        sprintf(fname2,"%s.catalog_particles",infile);
        sprintf(fname3,"%s.catalog_particles.unbound",infile);
        sprintf(fname4,"%s.sublevels.catalog_groups",infile);
        sprintf(fname5,"%s.sublevels.catalog_particles",infile);
        sprintf(fname6,"%s.sublevels.catalog_particles.unbound",infile);
        sprintf(fname7,"%s.catalog_parttypes",infile);
        sprintf(fname8,"%s.catalog_parttypes.unbound",infile);
        sprintf(fname9,"%s.sublevels.catalog_parttypes",infile);
        sprintf(fname10,"%s.sublevels.catalog_parttypes.unbound",infile);
    }
    else {
        sprintf(fname1,"%s.catalog_groups.%d",infile,k);
        sprintf(fname2,"%s.catalog_particles.%d",infile,k);
        sprintf(fname3,"%s.catalog_particles.unbound.%d",infile,k);
        sprintf(fname4,"%s.sublevels.catalog_groups.%d",infile,k);
        sprintf(fname5,"%s.sublevels.catalog_particles.%d",infile,k);
        sprintf(fname6,"%s.sublevels.catalog_particles.unbound.%d",infile,k);
        sprintf(fname7,"%s.catalog_parttypes.%d",infile,k);
        sprintf(fname8,"%s.catalog_parttypes.unbound.%d",infile,k);
        sprintf(fname9,"%s.sublevels.catalog_parttypes.%d",infile,k);
        sprintf(fname10,"%s.sublevels.catalog_parttypes.unbound.%d",infile,k);
    }
    itemp=0;
    fnamearray[itemp++]=fname1;fnamearray[itemp++]=fname2;fnamearray[itemp++]=fname3;
    if (ifieldhalos) {
        fnamearray[itemp++]=fname4;fnamearray[itemp++]=fname5;fnamearray[itemp++]=fname6;
    }
    if (itypematch!=ALLTYPEMATCH) { 
        fnamearray[itemp++]=fname7;fnamearray[itemp++]=fname8;
        if (ifieldhalos) {
            fnamearray[itemp++]=fname9;fnamearray[itemp++]=fname10;
        }
    }
    if (ibinary) {
        Fgroup.open(fname1,ios::in|ios::binary);
        Fpart.open(fname2,ios::in|ios::binary);
        Fupart.open(fname3,ios::in|ios::binary);
        if (ifieldhalos) {
            Fsgroup.open(fname4,ios::in|ios::binary);
            Fspart.open(fname5,ios::in|ios::binary);
            Fsupart.open(fname6,ios::in|ios::binary);
        }
        if (itypematch!=ALLTYPEMATCH) {
            Fparttype.open(fname7,ios::in|ios::binary);
            Fuparttype.open(fname8,ios::in|ios::binary);
            if (ifieldhalos) {
                Fsparttype.open(fname9,ios::in|ios::binary);
                Fsuparttype.open(fname10,ios::in|ios::binary);
            }
        }
    }
    else {
        Fgroup.open(fname1,ios::in);
        Fpart.open(fname2,ios::in);
        Fupart.open(fname3,ios::in);
        if (ifieldhalos) {
            Fsgroup.open(fname4,ios::in);
            Fspart.open(fname5,ios::in);
            Fsupart.open(fname6,ios::in);
        }
        if (itypematch!=ALLTYPEMATCH) {
            Fparttype.open(fname7,ios::in);
            Fuparttype.open(fname8,ios::in);
            if (ifieldhalos) {
                Fsparttype.open(fname9,ios::in);
                Fsuparttype.open(fname10,ios::in);
            }
        }
    }
    itemp=0;
    Farray[itemp++]=&Fgroup;Farray[itemp++]=&Fpart;Farray[itemp++]=&Fupart;
    if (ifieldhalos) {
        Farray[itemp++]=&Fsgroup;Farray[itemp++]=&Fspart;Farray[itemp++]=&Fsupart;
    }
    if (itypematch!=ALLTYPEMATCH) {
        Farray[itemp++]=&Fparttype;Farray[itemp++]=&Fuparttype;
        if (ifieldhalos) {
            Farray[itemp++]=&Fsparttype;Farray[itemp++]=&Fsuparttype;
        }
    }
    for (int i=0;i<numfiletypes;i++) {
        if(!Farray[i]->is_open()){
            cerr<<"can't open "<<fnamearray[i]<<endl;
            MPI_Abort(MPI_COMM_WORLD,9);
        }
        else if (iverbose) cout<<"open "<<fnamearray[i]<<endl;
    }
    //check header make sure that number of mpi values considered agree
    if (ibinary) {
        Fgroup.read((char*)&itask,sizeof(int));
        Fgroup.read((char*)&nprocs,sizeof(int));
        if (ifieldhalos) {
            Fsgroup.read((char*)&itask,sizeof(int));
            Fsgroup.read((char*)&nprocs,sizeof(int));
        }
    }
    else {
        Fgroup>>itask>>nprocs;
        if (ifieldhalos) {
            Fsgroup>>itask>>nprocs;
        }
    }
    if (nprocs!=mpi_ninput&&mpi_ninput>0) {
        cout<<"Error, number of mpi outputs was set to "<<mpi_ninput<<" but file indicates there are "<<nprocs<<endl;
        cout<<"Correcting to this number and proceeding"<<endl;
        nmpicount=mpi_ninput=nprocs;
    }
    else if ((mpi_ninput==0&&nprocs!=1)) {
        cout<<"Error, number of mpi outputs was set to zero but file indicates there are more than one mpi output"<<endl;
        cout<<"Terminating"<<endl;
        exit(9);
    }

    if (ibinary) {
        Fgroup.read((char*)&ltemp,sizeof(unsigned long));
        nglocal=ltemp;
        Fgroup.read((char*)&ltemp,sizeof(unsigned long));
        TotalNumberofHalos=ltemp;
        if (ifieldhalos) {
            Fsgroup.read((char*)&ltemp,sizeof(unsigned long));
            nsglocal=ltemp;
            Fsgroup.read((char*)&ltemp,sizeof(unsigned long));
        }
    }
    else {
        Fgroup>>nglocal;
        Fgroup>>TotalNumberofHalos;
        if (ifieldhalos) {
            Fsgroup>>nsglocal;
            Fsgroup>>itemp;
        }
    }

    if (iverbose) cout<<infile<<" has "<<nglocal<<endl;
    //allocate memory for halos
    if (k==0) {
        if (iverbose) cout<<" and the total number of halos in all files is "<<TotalNumberofHalos<<endl;
        numhalos=TotalNumberofHalos;
    }

    //read the header info before reading lengths of the halos
    if (ibinary) {
        Fpart.read((char*)&itask,sizeof(int));
        Fpart.read((char*)&nprocs,sizeof(int));
        Fupart.read((char*)&itask,sizeof(int));
        Fupart.read((char*)&nprocs,sizeof(int));

        Fpart.read((char*)&nids,sizeof(unsigned long));
        Fpart.read((char*)&nidstot,sizeof(unsigned long));//total ids
        Fupart.read((char*)&nuids,sizeof(unsigned long));
        Fupart.read((char*)&nuidstot,sizeof(unsigned long));//total ids

        if (ifieldhalos) {
            Fspart.read((char*)&itask,sizeof(int));
            Fspart.read((char*)&nprocs,sizeof(int));
            Fsupart.read((char*)&itask,sizeof(int));
            Fsupart.read((char*)&nprocs,sizeof(int));

            Fspart.read((char*)&nsids,sizeof(unsigned long));
            Fspart.read((char*)&nsidstot,sizeof(unsigned long));//total ids
            Fsupart.read((char*)&nsuids,sizeof(unsigned long));
            Fsupart.read((char*)&nsuidstot,sizeof(unsigned long));//total ids
        }

        if (itypematch!=ALLTYPEMATCH) {
            Fparttype.read((char*)&itask,sizeof(int));
            Fparttype.read((char*)&nprocs,sizeof(int));
            Fuparttype.read((char*)&itask,sizeof(int));
            Fuparttype.read((char*)&nprocs,sizeof(int));

            Fparttype.read((char*)&nids,sizeof(unsigned long));
            Fparttype.read((char*)&nidstot,sizeof(unsigned long));//total ids
            Fuparttype.read((char*)&nuids,sizeof(unsigned long));
            Fuparttype.read((char*)&nuidstot,sizeof(unsigned long));//total ids

            if (ifieldhalos) {
                Fsparttype.read((char*)&itask,sizeof(int));
                Fsparttype.read((char*)&nprocs,sizeof(int));
                Fsuparttype.read((char*)&itask,sizeof(int));
                Fsuparttype.read((char*)&nprocs,sizeof(int));

                Fsparttype.read((char*)&nsids,sizeof(unsigned long));
                Fsparttype.read((char*)&nsidstot,sizeof(unsigned long));//total ids
                Fsuparttype.read((char*)&nsuids,sizeof(unsigned long));
                Fsuparttype.read((char*)&nsuidstot,sizeof(unsigned long));//total ids
            }
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
        if (itypematch!=ALLTYPEMATCH) {
            Fparttype>>itask>>nprocs;
            Fuparttype>>itask>>nprocs;
            Fparttype>>nids>>nidstot;
            Fuparttype>>nuids>>nuidstot;
            if (ifieldhalos) {
                Fsparttype>>itask>>nprocs;
                Fsuparttype>>itask>>nprocs;
                Fsparttype>>nsids>>nsidstot;
                Fsuparttype>>nsuids>>nsuidstot;
            }
        }
    }

    //cout<<infile<<" has "<<nglocal<<" "<<nids<<" "<<nuids<<endl;
    if (nglocal>0) {
        numingroup=new Int_t[nglocal+1];
        offset=new Int_t[nglocal+1];
        uoffset=new Int_t[nglocal+1];
        numingroupbound=new Int_t[nglocal+1];
        idval=new Int_t[nids+nuids+1];
        if (itypematch!=ALLTYPEMATCH) typeval=new UInt_t[nids+nuids+1];

        if (ibinary) {
        Fgroup.read((char*)numingroup,sizeof(Int_t)*nglocal);
        Fgroup.read((char*)offset,sizeof(Int_t)*nglocal);
        Fgroup.read((char*)uoffset,sizeof(Int_t)*nglocal);
        }
        else {
        for (Int_t i=0;i<nglocal;i++) Fgroup>>numingroup[i];
        for (Int_t i=0;i<nglocal;i++) Fgroup>>offset[i];
        for (Int_t i=0;i<nglocal;i++) Fgroup>>uoffset[i];
        }

        //now read bound particle list
        if (ibinary) {
            Fpart.read((char*)idval,sizeof(Int_t)*nids);
            Fupart.read((char*)&idval[nids],sizeof(Int_t)*nuids);
            if (itypematch!=ALLTYPEMATCH) {
                Fparttype.read((char*)typeval,sizeof(UInt_t)*nids);
                Fuparttype.read((char*)&typeval[nids],sizeof(UInt_t)*nuids);
            }
        }
        else {
            for (Int_t i=0;i<nids;i++) Fpart>>idval[i];
            for (Int_t i=0;i<nuids;i++) Fupart>>idval[i+nids];
            if (itypematch!=ALLTYPEMATCH) {
                for (Int_t i=0;i<nids;i++) Fparttype>>typeval[i];
                for (Int_t i=0;i<nuids;i++) Fuparttype>>typeval[i+nids];
            }
        }

        for (Int_t i=0;i<nglocal-1;i++) {
            numingroupbound[i]=Halo[i+noffset].NumberofParticles-(uoffset[i+1]-uoffset[i]);
        }
        numingroupbound[nglocal-1]=Halo[nglocal-1+noffset].NumberofParticles-(nuids-uoffset[nglocal-1]);
        //if don't care about particle type, simply store info
        if (itypematch==ALLTYPEMATCH) {
            for (Int_t i=0;i<nglocal;i++) {
                for (Int_t j=0;j<numingroupbound[i];j++) Halo[i+noffset].ParticleID[j]=idval[offset[i]+j];
                nn=numingroup[i]-numingroupbound[i];
                for (Int_t j=0;j<nn;j++) Halo[i+noffset].ParticleID[j+numingroupbound[i]]=idval[uoffset[i]+j+nids];
            }
        }
        else {
            for (Int_t i=0;i<nglocal;i++) {
                counter=0;
                for (Int_t j=0;j<numingroupbound[i];j++) if (CheckType(typeval[offset[i]+j],itypematch)) Halo[i+noffset].ParticleID[counter++]=idval[offset[i]+j];
                nn=numingroup[i]-numingroupbound[i];
                for (Int_t j=0;j<nn;j++) if (CheckType(typeval[uoffset[i]+j+nids],itypematch)) Halo[i+noffset].ParticleID[counter++]=idval[uoffset[i]+j+nids];
                Halo[i+noffset].NumberofParticles=counter;
            }
        }
        delete[] idval;
        delete[] offset;
        delete[] uoffset;
        delete[] numingroupbound;
        delete[] numingroup;
        if (itypematch!=ALLTYPEMATCH) delete[] typeval;
    }

    for (Int_t i=0;i<nglocal+nsglocal;i++) Halo[i+noffset].haloID=i+1+noffset;
    if (ifieldhalos) {
        //cout<<infile<<" has sublevels with "<<nsglocal<<" "<<nsids<<" "<<nsuids<<endl;
        if (nsglocal>0) {
            //now read substructure data
            numingroup=new Int_t[nsglocal];
            offset=new Int_t[nsglocal];
            uoffset=new Int_t[nsglocal];
            numingroupbound=new Int_t[nsglocal];
            idval=new Int_t[nsids+nsuids];
            if (itypematch!=ALLTYPEMATCH) typeval=new UInt_t[nids+nuids+1];
            if (ibinary) {
            Fsgroup.read((char*)numingroup,sizeof(Int_t)*nsglocal);
            Fsgroup.read((char*)offset,sizeof(Int_t)*nsglocal);
            Fsgroup.read((char*)uoffset,sizeof(Int_t)*nsglocal);
            }
            else {
            for (Int_t i=0;i<nsglocal;i++) Fsgroup>>numingroup[i];
            for (Int_t i=0;i<nsglocal;i++) Fsgroup>>offset[i];
            for (Int_t i=0;i<nsglocal;i++) Fsgroup>>uoffset[i];
            }
            //now read bound particle list
            if (ibinary) {
            Fspart.read((char*)idval,sizeof(Int_t)*nsids);
            Fsupart.read((char*)&idval[nsids],sizeof(Int_t)*nsuids);
            if (itypematch!=ALLTYPEMATCH) {
                Fsparttype.read((char*)typeval,sizeof(UInt_t)*nsids);
                Fsuparttype.read((char*)&typeval[nsids],sizeof(UInt_t)*nsuids);
            }
            }
            else {
            for (Int_t i=0;i<nsids;i++) Fspart>>idval[i];
            for (Int_t i=0;i<nsuids;i++) Fsupart>>idval[i+nsids];
            if (itypematch!=ALLTYPEMATCH) {
                for (Int_t i=0;i<nsids;i++) Fsparttype>>typeval[i];
                for (Int_t i=0;i<nsuids;i++) Fsuparttype>>typeval[i+nsids];
            }
            }

            for (Int_t i=0;i<nsglocal-1;i++) {
                numingroupbound[i]=Halo[i+nglocal+noffset].NumberofParticles-(uoffset[i+1]-uoffset[i]);
            }
            numingroupbound[nsglocal-1]=Halo[nglocal+nsglocal-1+noffset].NumberofParticles-(nsuids-uoffset[nsglocal-1]);
            if (itypematch==ALLTYPEMATCH) {
                for (Int_t i=0;i<nsglocal;i++) {
                    for (Int_t j=0;j<numingroupbound[i];j++) Halo[i+nglocal+noffset].ParticleID[j]=idval[offset[i]+j];
                    nn=numingroup[i]-numingroupbound[i];
                    for (Int_t j=0;j<nn;j++) Halo[i+nglocal+noffset].ParticleID[j+numingroupbound[i]]=idval[uoffset[i]+j+nsids];
                }
            }
            else {
                for (Int_t i=0;i<nsglocal;i++) {
                    counter=0;
                    for (Int_t j=0;j<numingroupbound[i];j++) if (CheckType(typeval[offset[i]+j],itypematch)) Halo[i+nglocal+noffset].ParticleID[counter++]=idval[offset[i]+j];
                    nn=numingroup[i]-numingroupbound[i];
                    for (Int_t j=0;j<nn;j++) if (CheckType(typeval[uoffset[i]+j+nids],itypematch)) Halo[i+nglocal+noffset].ParticleID[counter++]=idval[uoffset[i]+j+nids];
                    Halo[i+nglocal+noffset].NumberofParticles=counter;
                }
            }
            delete[] idval;
            delete[] offset;
            delete[] uoffset;
            delete[] numingroupbound;
            delete[] numingroup;
            if (itypematch!=ALLTYPEMATCH) delete[] typeval;
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
    if (itypematch!=ALLTYPEMATCH) {
        Fparttype.close();
        Fuparttype.close();
        if (ifieldhalos) {
            Fsparttype.close();
            Fsuparttype.close();
        }
    }

    }
}
#endif


///Read halo catalog data from an individual snapshot;
HaloData *ReadHaloGroupCatalogData(char* infile, Int_t &numhalos, int mpi_ninput, int ibinary, int ifieldhalos, int itypematch, int iverbose)
{
    HaloData *Halo;
    Int_t itemp;
    unsigned long ltemp;
    Int_t noffset,numfiletypes;
    long unsigned nparts,haloid;
    long unsigned TotalNumberofHalos;
    char fname1[1000],fname2[1000],fname3[1000];
    char fname4[1000],fname5[1000],fname6[1000];
    char fname7[1000],fname8[1000];
    char fname9[1000],fname10[1000];
    fstream Fgroup,Fpart,Fupart; //field objects
    fstream Fsgroup,Fspart,Fsupart; //sublevels
    fstream Fparttype,Fuparttype; //field objects
    fstream Fsparttype,Fsuparttype; //sublevels
    char *fnamearray[10];
    fstream *Farray[10];
    unsigned long nids,nuids,nsids,nsuids,nglocal,nsglocal;
    unsigned long nidstot,nuidstot,nsidstot,nsuidstot;
    Int_t *numingroup,*numingroupbound,*offset,*uoffset;
    Int_t counter,nn;
    int nmpicount,itask,nprocs;
    Int_t *idval;
    UInt_t *typeval;
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif
    numfiletypes=3;
    if (ifieldhalos) numfiletypes+=3;
    if (itypematch!=ALLTYPEMATCH) {
        numfiletypes+=2;
        if (ifieldhalos) numfiletypes+=2;
    }

    ///adjust files read based on whether mpi output or not
    if (mpi_ninput==0) nmpicount=1;
    else nmpicount=mpi_ninput;
    //to offset the halo structure pointer given that data maybe split among multiple files
    noffset=0;
    for (int k=0;k<nmpicount;k++) {
        nglocal=nsglocal=0;
        //load the group catalogues and the related particle files
        if (mpi_ninput==0) {
            sprintf(fname1,"%s.catalog_groups",infile);
            sprintf(fname2,"%s.catalog_particles",infile);
            sprintf(fname3,"%s.catalog_particles.unbound",infile);
            sprintf(fname4,"%s.sublevels.catalog_groups",infile);
            sprintf(fname5,"%s.sublevels.catalog_particles",infile);
            sprintf(fname6,"%s.sublevels.catalog_particles.unbound",infile);
            sprintf(fname7,"%s.catalog_parttypes",infile);
            sprintf(fname8,"%s.catalog_parttypes.unbound",infile);
            sprintf(fname9,"%s.sublevels.catalog_parttypes",infile);
            sprintf(fname10,"%s.sublevels.catalog_parttypes.unbound",infile);
        }
        else {
            sprintf(fname1,"%s.catalog_groups.%d",infile,k);
            sprintf(fname2,"%s.catalog_particles.%d",infile,k);
            sprintf(fname3,"%s.catalog_particles.unbound.%d",infile,k);
            sprintf(fname4,"%s.sublevels.catalog_groups.%d",infile,k);
            sprintf(fname5,"%s.sublevels.catalog_particles.%d",infile,k);
            sprintf(fname6,"%s.sublevels.catalog_particles.unbound.%d",infile,k);
            sprintf(fname7,"%s.catalog_parttypes.%d",infile,k);
            sprintf(fname8,"%s.catalog_parttypes.unbound.%d",infile,k);
            sprintf(fname9,"%s.sublevels.catalog_parttypes.%d",infile,k);
            sprintf(fname10,"%s.sublevels.catalog_parttypes.unbound.%d",infile,k);
        }
        itemp=0;
        fnamearray[itemp++]=fname1;fnamearray[itemp++]=fname2;fnamearray[itemp++]=fname3;
        if (ifieldhalos) {
            fnamearray[itemp++]=fname4;fnamearray[itemp++]=fname5;fnamearray[itemp++]=fname6;
        }
        if (itypematch!=ALLTYPEMATCH) { 
            fnamearray[itemp++]=fname7;fnamearray[itemp++]=fname8;
            if (ifieldhalos) {
                fnamearray[itemp++]=fname9;fnamearray[itemp++]=fname10;
            }
        }
        if (ibinary) {
            Fgroup.open(fname1,ios::in|ios::binary);
            Fpart.open(fname2,ios::in|ios::binary);
            Fupart.open(fname3,ios::in|ios::binary);
            if (ifieldhalos) {
                Fsgroup.open(fname4,ios::in|ios::binary);
                Fspart.open(fname5,ios::in|ios::binary);
                Fsupart.open(fname6,ios::in|ios::binary);
            }
            if (itypematch!=ALLTYPEMATCH) {
                Fparttype.open(fname7,ios::in|ios::binary);
                Fuparttype.open(fname8,ios::in|ios::binary);
                if (ifieldhalos) {
                    Fsparttype.open(fname9,ios::in|ios::binary);
                    Fsuparttype.open(fname10,ios::in|ios::binary);
                }
            }
        }
        else {
            Fgroup.open(fname1,ios::in);
            Fpart.open(fname2,ios::in);
            Fupart.open(fname3,ios::in);
            if (ifieldhalos) {
                Fsgroup.open(fname4,ios::in);
                Fspart.open(fname5,ios::in);
                Fsupart.open(fname6,ios::in);
            }
            if (itypematch!=ALLTYPEMATCH) {
                Fparttype.open(fname7,ios::in);
                Fuparttype.open(fname8,ios::in);
                if (ifieldhalos) {
                    Fsparttype.open(fname9,ios::in);
                    Fsuparttype.open(fname10,ios::in);
                }
            }
        }
        itemp=0;
        Farray[itemp++]=&Fgroup;Farray[itemp++]=&Fpart;Farray[itemp++]=&Fupart;
        if (ifieldhalos) {
            Farray[itemp++]=&Fsgroup;Farray[itemp++]=&Fspart;Farray[itemp++]=&Fsupart;
        }
        if (itypematch!=ALLTYPEMATCH) {
            Farray[itemp++]=&Fparttype;Farray[itemp++]=&Fuparttype;
            if (ifieldhalos) {
                Farray[itemp++]=&Fsparttype;Farray[itemp++]=&Fsuparttype;
            }
        }
        for (int i=0;i<numfiletypes;i++) {
            if(!Farray[i]->is_open()){
                cerr<<"can't open "<<fnamearray[i]<<endl;
    #ifdef USEMPI
                MPI_Abort(MPI_COMM_WORLD,9);
    #else
                exit(9);
    #endif
            }
            else if (iverbose) cout<<"open "<<fnamearray[i]<<endl;
        }
        //check header make sure that number of mpi values considered agree
        if (ibinary) {
            Fgroup.read((char*)&itask,sizeof(int));
            Fgroup.read((char*)&nprocs,sizeof(int));
            if (ifieldhalos) {
                Fsgroup.read((char*)&itask,sizeof(int));
                Fsgroup.read((char*)&nprocs,sizeof(int));
            }
        }
        else {
            Fgroup>>itask>>nprocs;
            if (ifieldhalos) {
                Fsgroup>>itask>>nprocs;
            }
        }
        if (nprocs!=mpi_ninput&&mpi_ninput>0) {
            cout<<"Error, number of mpi outputs was set to "<<mpi_ninput<<" but file indicates there are "<<nprocs<<endl;
            cout<<"Correcting to this number and proceeding"<<endl;
            nmpicount=mpi_ninput=nprocs;
        }
        else if ((mpi_ninput==0&&nprocs!=1)) {
            cout<<"Error, number of mpi outputs was set to zero but file indicates there are more than one mpi output"<<endl;
            cout<<"Terminating"<<endl;
            exit(9);
        }

        if (ibinary) {
            Fgroup.read((char*)&ltemp,sizeof(unsigned long));
            nglocal=ltemp;
            Fgroup.read((char*)&ltemp,sizeof(unsigned long));
            TotalNumberofHalos=ltemp;
            if (ifieldhalos) {
                Fsgroup.read((char*)&ltemp,sizeof(unsigned long));
                nsglocal=ltemp;
                Fsgroup.read((char*)&ltemp,sizeof(unsigned long));
            }
        }
        else {
            Fgroup>>nglocal;
            Fgroup>>TotalNumberofHalos;
            if (ifieldhalos) {
                Fsgroup>>nsglocal;
                Fsgroup>>itemp;
            }
        }

        if (iverbose) cout<<infile<<" has "<<nglocal<<endl;
        //allocate memory for halos
        if (k==0) {
            if (iverbose) cout<<" and the total number of halos in all files is "<<TotalNumberofHalos<<endl;
            Halo=new HaloData[TotalNumberofHalos];
            numhalos=TotalNumberofHalos;
        }

        //read the header info before reading lengths of the halos
        if (ibinary) {
            Fpart.read((char*)&itask,sizeof(int));
            Fpart.read((char*)&nprocs,sizeof(int));
            Fupart.read((char*)&itask,sizeof(int));
            Fupart.read((char*)&nprocs,sizeof(int));

            Fpart.read((char*)&nids,sizeof(unsigned long));
            Fpart.read((char*)&nidstot,sizeof(unsigned long));//total ids
            Fupart.read((char*)&nuids,sizeof(unsigned long));
            Fupart.read((char*)&nuidstot,sizeof(unsigned long));//total ids

            if (ifieldhalos) {
                Fspart.read((char*)&itask,sizeof(int));
                Fspart.read((char*)&nprocs,sizeof(int));
                Fsupart.read((char*)&itask,sizeof(int));
                Fsupart.read((char*)&nprocs,sizeof(int));

                Fspart.read((char*)&nsids,sizeof(unsigned long));
                Fspart.read((char*)&nsidstot,sizeof(unsigned long));//total ids
                Fsupart.read((char*)&nsuids,sizeof(unsigned long));
                Fsupart.read((char*)&nsuidstot,sizeof(unsigned long));//total ids
            }

            if (itypematch!=ALLTYPEMATCH) {
                Fparttype.read((char*)&itask,sizeof(int));
                Fparttype.read((char*)&nprocs,sizeof(int));
                Fuparttype.read((char*)&itask,sizeof(int));
                Fuparttype.read((char*)&nprocs,sizeof(int));

                Fparttype.read((char*)&nids,sizeof(unsigned long));
                Fparttype.read((char*)&nidstot,sizeof(unsigned long));//total ids
                Fuparttype.read((char*)&nuids,sizeof(unsigned long));
                Fuparttype.read((char*)&nuidstot,sizeof(unsigned long));//total ids

                if (ifieldhalos) {
                    Fsparttype.read((char*)&itask,sizeof(int));
                    Fsparttype.read((char*)&nprocs,sizeof(int));
                    Fsuparttype.read((char*)&itask,sizeof(int));
                    Fsuparttype.read((char*)&nprocs,sizeof(int));

                    Fsparttype.read((char*)&nsids,sizeof(unsigned long));
                    Fsparttype.read((char*)&nsidstot,sizeof(unsigned long));//total ids
                    Fsuparttype.read((char*)&nsuids,sizeof(unsigned long));
                    Fsuparttype.read((char*)&nsuidstot,sizeof(unsigned long));//total ids
                }
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
            if (itypematch!=ALLTYPEMATCH) {
                Fparttype>>itask>>nprocs;
                Fuparttype>>itask>>nprocs;
                Fparttype>>nids>>nidstot;
                Fuparttype>>nuids>>nuidstot;
                if (ifieldhalos) {
                    Fsparttype>>itask>>nprocs;
                    Fsuparttype>>itask>>nprocs;
                    Fsparttype>>nsids>>nsidstot;
                    Fsuparttype>>nsuids>>nsuidstot;
                }
            }
        }

        //cout<<infile<<" has "<<nglocal<<" "<<nids<<" "<<nuids<<endl;
        if (nglocal>0) {
            numingroup=new Int_t[nglocal+1];
            offset=new Int_t[nglocal+1];
            uoffset=new Int_t[nglocal+1];
            numingroupbound=new Int_t[nglocal+1];
            idval=new Int_t[nids+nuids+1];
            if (itypematch!=ALLTYPEMATCH) typeval=new UInt_t[nids+nuids+1];

            if (ibinary) {
                Fgroup.read((char*)numingroup,sizeof(Int_t)*nglocal);
                for (Int_t i=0;i<nglocal;i++) {
                    Halo[i+noffset].Alloc(numingroup[i]);
                }
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

            //now read bound particle list
            if (ibinary) {
                Fpart.read((char*)idval,sizeof(Int_t)*nids);
                Fupart.read((char*)&idval[nids],sizeof(Int_t)*nuids);
                if (itypematch!=ALLTYPEMATCH) {
                    Fparttype.read((char*)typeval,sizeof(UInt_t)*nids);
                    Fuparttype.read((char*)&typeval[nids],sizeof(UInt_t)*nuids);
                }
            }
            else {
                for (Int_t i=0;i<nids;i++) Fpart>>idval[i];
                for (Int_t i=0;i<nuids;i++) Fupart>>idval[i+nids];
                if (itypematch!=ALLTYPEMATCH) {
                    for (Int_t i=0;i<nids;i++) Fparttype>>typeval[i];
                    for (Int_t i=0;i<nuids;i++) Fuparttype>>typeval[i+nids];
                }
            }

            for (Int_t i=0;i<nglocal-1;i++) {
                numingroupbound[i]=Halo[i+noffset].NumberofParticles-(uoffset[i+1]-uoffset[i]);
            }
            numingroupbound[nglocal-1]=Halo[nglocal-1+noffset].NumberofParticles-(nuids-uoffset[nglocal-1]);
            //if don't care about particle type, simply store info
            if (itypematch==ALLTYPEMATCH) {
                for (Int_t i=0;i<nglocal;i++) {
                    for (Int_t j=0;j<numingroupbound[i];j++) Halo[i+noffset].ParticleID[j]=idval[offset[i]+j];
                    nn=numingroup[i]-numingroupbound[i];
                    for (Int_t j=0;j<nn;j++) Halo[i+noffset].ParticleID[j+numingroupbound[i]]=idval[uoffset[i]+j+nids];
                }
            }
            else {
                for (Int_t i=0;i<nglocal;i++) {
                    counter=0;
                    for (Int_t j=0;j<numingroupbound[i];j++) if (CheckType(typeval[offset[i]+j],itypematch)) Halo[i+noffset].ParticleID[counter++]=idval[offset[i]+j];
                    nn=numingroup[i]-numingroupbound[i];
                    for (Int_t j=0;j<nn;j++) if (CheckType(typeval[uoffset[i]+j+nids],itypematch)) Halo[i+noffset].ParticleID[counter++]=idval[uoffset[i]+j+nids];
                    Halo[i+noffset].NumberofParticles=counter;
                }
            }
            delete[] idval;
            delete[] offset;
            delete[] uoffset;
            delete[] numingroupbound;
            delete[] numingroup;
            if (itypematch!=ALLTYPEMATCH) delete[] typeval;
        }

        for (Int_t i=0;i<nglocal+nsglocal;i++) Halo[i+noffset].haloID=i+1+noffset;
        if (ifieldhalos) {
        //cout<<infile<<" has sublevels with "<<nsglocal<<" "<<nsids<<" "<<nsuids<<endl;
        if (nsglocal>0) {
            //now read substructure data
            numingroup=new Int_t[nsglocal];
            offset=new Int_t[nsglocal];
            uoffset=new Int_t[nsglocal];
            numingroupbound=new Int_t[nsglocal];
            idval=new Int_t[nsids+nsuids];
            if (itypematch!=ALLTYPEMATCH) typeval=new UInt_t[nids+nuids+1];
            if (ibinary) {
                Fsgroup.read((char*)numingroup,sizeof(Int_t)*nsglocal);
                for (Int_t i=0;i<nsglocal;i++) {
                    Halo[i+nglocal+noffset].Alloc(numingroup[i]);
                }
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
            //now read bound particle list
            if (ibinary) {
                Fspart.read((char*)idval,sizeof(Int_t)*nsids);
                Fsupart.read((char*)&idval[nsids],sizeof(Int_t)*nsuids);
                if (itypematch!=ALLTYPEMATCH) {
                    Fsparttype.read((char*)typeval,sizeof(UInt_t)*nsids);
                    Fsuparttype.read((char*)&typeval[nsids],sizeof(UInt_t)*nsuids);
                }
            }
            else {
                for (Int_t i=0;i<nsids;i++) Fspart>>idval[i];
                for (Int_t i=0;i<nsuids;i++) Fsupart>>idval[i+nsids];
                if (itypematch!=ALLTYPEMATCH) {
                    for (Int_t i=0;i<nsids;i++) Fsparttype>>typeval[i];
                    for (Int_t i=0;i<nsuids;i++) Fsuparttype>>typeval[i+nsids];
                }
            }

            for (Int_t i=0;i<nsglocal-1;i++) {
                numingroupbound[i]=Halo[i+nglocal+noffset].NumberofParticles-(uoffset[i+1]-uoffset[i]);
            }
            numingroupbound[nsglocal-1]=Halo[nglocal+nsglocal-1+noffset].NumberofParticles-(nsuids-uoffset[nsglocal-1]);
            if (itypematch==ALLTYPEMATCH) {
                for (Int_t i=0;i<nsglocal;i++) {
                    for (Int_t j=0;j<numingroupbound[i];j++) Halo[i+nglocal+noffset].ParticleID[j]=idval[offset[i]+j];
                    nn=numingroup[i]-numingroupbound[i];
                    for (Int_t j=0;j<nn;j++) Halo[i+nglocal+noffset].ParticleID[j+numingroupbound[i]]=idval[uoffset[i]+j+nsids];
                }
            }
            else {
                for (Int_t i=0;i<nsglocal;i++) {
                    counter=0;
                    for (Int_t j=0;j<numingroupbound[i];j++) if (CheckType(typeval[offset[i]+j],itypematch)) Halo[i+nglocal+noffset].ParticleID[counter++]=idval[offset[i]+j];
                    nn=numingroup[i]-numingroupbound[i];
                    for (Int_t j=0;j<nn;j++) if (CheckType(typeval[uoffset[i]+j+nids],itypematch)) Halo[i+nglocal+noffset].ParticleID[counter++]=idval[uoffset[i]+j+nids];
                    Halo[i+nglocal+noffset].NumberofParticles=counter;
                }
            }
            delete[] idval;
            delete[] offset;
            delete[] uoffset;
            delete[] numingroupbound;
            delete[] numingroup;
            if (itypematch!=ALLTYPEMATCH) delete[] typeval;
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
        if (itypematch!=ALLTYPEMATCH) {
            Fparttype.close();
            Fuparttype.close();
            if (ifieldhalos) {
                Fsparttype.close();
                Fsuparttype.close();
            }
        }

    }

    return Halo;
}

