#ifndef GADGETITEMS_H
#define GADGETITEMS_H

//for endian independance
#include "endianutils.h"

///for gadget coords
#ifdef GADGETDOUBLEPRECISION
#define FLOAT double
#else
#define FLOAT float
#endif

///for gadget mass
#ifdef GADGETSINGLEMASSPRECISION
#define REAL float
#else
#define REAL double
#endif

///for gadget IDS
#ifdef GADGETLONGID
#define GADGETIDTYPE long long
#else
#define GADGETIDTYPE unsigned int
#endif

///number of particle types
#define NGTYPE 6
///define gadget particle types
//@{
#define GGASTYPE 0
#define GDMTYPE 1
#define GDM2TYPE 2
#define GDM3TYPE 3
#define GSTARTYPE 4
#define GBHTYPE 5
//@}

//for gadget integers
//#define INTEGER int
#define INTEGER unsigned int
//#define INTEGER long long

//gadget stuff

//how many chunks of a gadget array to read in one go
#define GADGETCHUNKSIZE 200000

///for waves data, u, rho, Ne, Nh, HSML contiguous block
#define NUMGADGETSPHBLOCKS 5 
///for waves data, star only has stellar age
#define NUMGADGETSTARBLOCKS 1 
///for waves data, bh not present
#define NUMGADGETBHBLOCKS 1

struct gadget_header
{
    INTEGER     npart[NGTYPE];
    double      mass[NGTYPE];
    double      time;
    double      redshift;
    int         flag_sfr;
    int         flag_feedback;
    INTEGER     npartTotal[NGTYPE];
    int         flag_cooling;
    int         num_files;
    double      BoxSize;
    double      Omega0;
    double      OmegaLambda;
    double      HubbleParam;
    int         flag_stellar_age;
    int         flag_metals;
    INTEGER     npartTotalHW[NGTYPE];
    int         flag_entropy_ICs;
    char        fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8-2*4-6*4-1*4];  /* fills to 256 Bytes, here calculation is based on 6 gadget particle types with default header stuff*/
    //make endian independent when read
    void Endian()
    {
        for(int i = 0; i < 6; i++)
        {
        npart[i]= LittleInt( npart[i] );
        mass[i] = LittleDouble( mass[i] );
        npartTotal[i]=LittleInt(npartTotal[i]);
        }
        flag_sfr = LittleInt(flag_sfr);
        flag_feedback = LittleInt(flag_feedback);
        flag_cooling = LittleInt(flag_cooling);
        num_files = LittleInt(num_files);
        time = LittleDouble(time);
        redshift = LittleDouble(redshift);
        BoxSize = LittleDouble(BoxSize);
        Omega0 = LittleDouble(Omega0);
        OmegaLambda = LittleDouble(OmegaLambda);
        HubbleParam = LittleDouble(HubbleParam);
    }
};

struct gadget_particle_data 
{
  FLOAT  Pos[3];
  FLOAT  Vel[3];
  REAL  Mass;
  int    Type;

  FLOAT  Rho, U, Temp, Ne;

  void Endian()
  {
    for(int i = 0; i < 3; ++i) //our compiler will unroll this
    {
      Pos[i]= LittleDouble(Pos[i]);
      Vel[i]= LittleDouble(Vel[i]);
    }
    //type is read in from header so does not
    //need to be altered.
    //Type = LittleInt(Type);
    //because mass can also be read from header
    //should not be automatically altered.
    //Mass = LittleFloat(Mass);
    Rho = LittleDouble(Rho);
    U = LittleDouble(U);
    Temp = LittleDouble(Temp);
    Ne = LittleDouble(Ne);
  }
};

inline int find_files(char *fname)
{
  FILE *fd;
  char buf[200], buf1[200];
  int dummy;
  struct gadget_header header1;


  sprintf(buf, "%s.%d", fname, 0);
  sprintf(buf1, "%s", fname);

  if((fd = fopen(buf, "r")))
    {
      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);
      fclose(fd);
      header1.Endian(); //change endian

      return header1.num_files;
    }

  if((fd = fopen(buf1, "r")))
    {
      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);
      fclose(fd);
      header1.Endian(); //change endian
      return 1;
    }

  printf("Error. Can't find snapshot!\nneither as `%s'\nnor as `%s'\n\n", buf, buf1);

  exit(1);
  return 0;
}

//get number of particles (-1 is all, -2 is all dark, otherwise specify specific gadget type
inline Int_t get_nbodies(char *fname, int ptype=-1)
{
    FILE *fd;
    char buf[2000], buf1[2000];
    int j, dummy;
    Int_t nbodies=0;
    struct gadget_header header1;
    int ihwflag=0;
    char DATA[5];

    sprintf(buf, "%s.%d", fname, 0);
    sprintf(buf1, "%s", fname);
    //try first file name with .0 and then single
    if(!(fd = fopen(buf, "r"))) {
        if(!(fd = fopen(buf1, "r"))) {
            printf("Error. Can't find snapshot!\nneither as `%s'\nnor as `%s'\n\n", buf, buf1);
            exit(1);
            return 0;
        }
    }

#ifdef GADGET2FORMAT
    fread(&dummy, sizeof(dummy), 1, fd);
    fread(&DATA[0], sizeof(char)*4, 1, fd);DATA[4] = '\0';
    fread(&dummy, sizeof(dummy), 1, fd);
    fread(&dummy, sizeof(dummy), 1, fd);
#endif
    fread(&dummy, sizeof(dummy), 1, fd);
    fread(&header1, sizeof(header1), 1, fd);
    fread(&dummy, sizeof(dummy), 1, fd);
    fclose(fd);
    header1.Endian(); //change endian
    //first check to see if highword ntotal is zero
    for(j=0, nbodies=0; j<NGTYPE; j++) if (header1.npartTotalHW[j]!=0) ihwflag=1;
    if (ptype==-1) {
        for(j=0, nbodies=0; j<NGTYPE; j++) {
            nbodies+=header1.npartTotal[j];
            if(ihwflag) nbodies+=((long long)(header1.npartTotalHW[j]) << 32);
        }
        return nbodies;
    }
    //assumes that by default dark matter is between 1 and type 4
    else if (ptype==-2) {
        for(j=GGASTYPE+1, nbodies=0; j<GSTARTYPE; j++) {
            nbodies+=header1.npartTotal[j];
            if(ihwflag) nbodies+=((long long)(header1.npartTotalHW[j]) << 32);
        }
        return nbodies;
    }
    else {
        nbodies=header1.npartTotal[ptype];
        if(ihwflag) nbodies+=((long long)(header1.npartTotalHW[ptype]) << 32);
        return nbodies;
    }
}

#endif 
