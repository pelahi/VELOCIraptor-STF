/*! \file stfio.cxx
 *  \brief this file contains routines for other types of particle lists
 */

#include "TreeFrog.h"

///Read halo data from an idividual snapshot;
HaloData *ReadHaloData(string &infile, Int_t &numhalos)
{
    HaloData *Halo;
    long unsigned i, j,nparts,haloid;
    long unsigned TotalNumberofHalos;

    fprintf(stderr,"- reading %s ... ",infile.c_str());
    FILE *f = fopen(infile.c_str(),"r");
    if(f == NULL) {
      fprintf(stderr,"could not open %s\nABORTING\n",infile.c_str());
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
HaloData *ReadNIFTYData(string &infile, Int_t &numhalos, int idcorrectflag, int hidoffset)
{
    HaloData *Halo;
    long unsigned i, j,nparts,haloid;
    long unsigned TotalNumberofHalos;
    long long idvaloffset=0;
    int ncount=0,type;
    long long idval;

    fprintf(stderr,"- reading %s ... \n",infile.c_str());
    FILE *f = fopen(infile.c_str(),"r");
    if(f == NULL) {
      fprintf(stderr,"could not open %s\nABORTING\n",infile.c_str());
      exit(1);
    }
    if (idcorrectflag==1) {
        if (strcmp(infile.c_str(),"somefilename")==0) {
            idvaloffset+=1;//offset some value
            }
    }
    if (idvaloffset!=0) cout<<infile<<" offseting ids by "<<idvaloffset<<endl;
    fscanf(f, "%ld", &TotalNumberofHalos);
    Halo = new HaloData[TotalNumberofHalos];
    for(i=0; i<TotalNumberofHalos; i++){
      fscanf(f, "%ld %ld",&(nparts),&(haloid));
      Halo[i].Alloc(nparts);
      Halo[i].haloID=haloid+hidoffset;
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


///Read void data from an idividual snapshot;
///This format is CellID, CellType, VoidID, ParticleID
///note that a void consists of void AND sheet particles
HaloData *ReadVoidData(string &infile, Int_t &numhalos, int idcorrectflag, int hidoffset)
{
    HaloData *Halo;
    long unsigned i, j,nparts,haloid;
    long unsigned TotalNumberofHalos;
    long long idvaloffset=0;
    int ncount=0,type;
    long long idval;
    fstream Fin;

    Fin.open(infile.c_str(),ios::in);
    if (!Fin) {
      cerr<< "could not open "<<infile<<"\nABORTING\n";
      exit(1);
    }
    else cout<<"reading "<<infile<<endl;
    if (idvaloffset!=0) cout<<infile<<" offseting ids by "<<idvaloffset<<endl;
    ///read the entire file to determine the number of lines (number of particles) and 
    ///the maximum void id
    string str;
    nparts=0;
    if (Fin.good())
    {
        while(getline(Fin,str)) 
        {
            istringstream ss(str);
            //now have string parse it 
            ss >> type >> type>> haloid>> idval;
            if (haloid>TotalNumberofHalos) TotalNumberofHalos=haloid;
            nparts++;
        }
    }
    Fin.clear();
    Fin.seekg(0, ios::beg);
    Halo = new HaloData[TotalNumberofHalos];
    int *numingroup=new int[TotalNumberofHalos];
    for(i=0; i<TotalNumberofHalos; i++) numingroup[i]=0;
    if (Fin.good())
    {
        while(getline(Fin,str)) 
        {
            istringstream ss(str);
            //now have string parse it 
            ss >> type >> type>> haloid>> idval;
            if (type==VOIDSTYPE) numingroup[haloid-1]++;
        }
    }
    Fin.clear();
    Fin.seekg(0, ios::beg);
    for(i=0; i<TotalNumberofHalos; i++) {Halo[i].Alloc(numingroup[i]);numingroup[i]=0;}
    if (Fin.good())
    {
        while(getline(Fin,str)) 
        {
            istringstream ss(str);
            //now have string parse it 
            ss >> type >> type>> haloid>> idval;
            if (type==VOIDSTYPE) Halo[haloid-1].ParticleID[numingroup[haloid-1]++]=idval+idvaloffset;
        }
    }
    Fin.close();
    numhalos=TotalNumberofHalos;
    return Halo;
}
