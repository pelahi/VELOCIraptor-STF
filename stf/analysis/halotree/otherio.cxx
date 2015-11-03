/*! \file stfio.cxx
 *  \brief this file contains routines for other types of particle lists
 */

#include "halomergertree.h"

///Read halo data from an idividual snapshot;
HaloData *ReadHaloData(char* infile, Int_t &numhalos)
{
    HaloData *Halo;
    long unsigned i, j,nparts,haloid;
    long unsigned TotalNumberofHalos;

    fprintf(stderr,"- reading %s ... ",infile);
    FILE *f = fopen(infile,"r");
    if(f == NULL) {
      fprintf(stderr,"could not open %s\nABORTING\n",infile);
      exit(1);
    }
    fscanf(f, "%ld", &TotalNumberofHalos);
    Halo = new HaloData[TotalNumberofHalos];
    for(i=0; i<TotalNumberofHalos; i++){
      fscanf(f, "%ld %ld",&(nparts),&(haloid));
      Halo[i].Alloc(nparts);
      Halo[i].haloID=haloid;
      double e,x,y,z,vx,vy,vz;
      for(j=0; j<Halo[i].NumberofParticles; j++){
        fscanf(f, "%ld %f %f %f %f %f %f %f",
               &(Halo[i].ParticleID[j]),
               &e,&x,&y,&z,&vx,&vy,&vz
             );
      }
    }
    numhalos=TotalNumberofHalos;
    fclose(f);
    fprintf(stderr,"done\n");
    return Halo;
}

///Read halo data from an idividual snapshot;
///as new AHF just outputs particle id and particle type AND particle type
///is important as only DM particles have continuous ids, adjust read so that 
///though total is allocated only use particles of type 1 for zoom simulations
HaloData *ReadNIFTYData(char* infile, Int_t &numhalos, int idcorrectflag)
{
    HaloData *Halo;
    long unsigned i, j,nparts,haloid;
    long unsigned TotalNumberofHalos;
    long long idvaloffset=0;
    int ncount=0,type;
    long long idval;

    fprintf(stderr,"- reading %s ... \n",infile);
    FILE *f = fopen(infile,"r");
    if(f == NULL) {
      fprintf(stderr,"could not open %s\nABORTING\n",infile);
      exit(1);
    }
    if (idcorrectflag==1) {
        if (strcmp(infile,"somefilename")==0) {
            idvaloffset+=1;//offset some value
            }
    }
    if (idvaloffset!=0) cout<<infile<<" offseting ids by "<<idvaloffset<<endl;
    fscanf(f, "%ld", &TotalNumberofHalos);
    Halo = new HaloData[TotalNumberofHalos];
    for(i=0; i<TotalNumberofHalos; i++){
      fscanf(f, "%ld %ld",&(nparts),&(haloid));
      Halo[i].Alloc(nparts);
      Halo[i].haloID=haloid;
      ncount=0;
      for(j=0; j<Halo[i].NumberofParticles; j++){
        fscanf(f, "%ld %d",
               &idval,&type
               );
           if (type==NIFTYDMTYPE) {
            Halo[i].ParticleID[ncount++]=idval+idvaloffset;
           }
      }
      Halo[i].NumberofParticles=ncount;
    }
    numhalos=TotalNumberofHalos;
    fclose(f);
    fprintf(stderr,"done\n");
    return Halo;
}
