/*! \file ramsesitems.h
 *  \brief this file contains definitions and routines for reading RAMSES files
 *  
 * RAMSES files come in 4 separate files, 
 * amr_ file containing information about the cells (number, positions) etc
 * hydro_ file containing information about the hydrodynamical properties of the cells (density, velocity, pressure/gamma, then passive scalars like metallicity)
 * part_ file containing information about number of dark matter/star particles/sink (BH?) and their position/velocity info
 * info_ file containing cosmological information
 */

#ifndef RAMSESITEMS_H
#define RAMSESITEMS_H

#ifdef RAMSESSINGLEPRECISION
#define RAMSESFLOAT float
#else
#define RAMSESFLOAT double
#endif
#ifdef RAMSESLONGIDS
#define RAMSESIDTYPE int
#else
#define RAMSESIDTYPE  unsigned int
#endif
#define RAMSESINTTYPE int 

///number of particle types
#define NRAMSESTYPE 5
///define gadget particle types
//@{
#define RAMSESGASTYPE 0
#define RAMSESDMTYPE 1
#define RAMSESSTARTYPE 2
#define RAMSESBHTYPE 3
#define RAMSESSINKTYPE 4
//@}

///how many particle properties are read from file in one go
#define RAMSESCHUNKSIZE 100000

#define NUMRAMSESSPHBLOCKS 1 

///\name Structures for the RAMSES interface
//@{

///data stored in the header group structure in the RAMSES format
struct RAMSES_Header {

    ///\name amr information
    //@{
    ///number of files ramses output was written with (is the number of cpus used)
    int         nfiles;
    ///number of dimensions of the cells
    int         ndim;
    ///2 to the power of number of dimensions
    int         twotondim;
    ///base mesh resolution vector
    int         nx,ny,nz;
    ///max level of refinements
    int         nlevelmax;
    ///max number of grids in file (per domain?)
    int         ngridmax;
    ///max number of boundary (ghost) cells in file (per domain?)
    int         nboundary;
    ///number of active grids 
    int         ngridactive;
    //@}
    ///\name hydro info 
    //@{
    ///number of hydro variables
    int         nvarh;
    int         levelmin, levelmax;
    RAMSESFLOAT gamma_index;
    //@}

    RAMSESFLOAT BoxSize;
    int         npartlocal;
    int         npart[NRAMSESTYPE];
    int         npartTotal[NRAMSESTYPE];
    int         npartTotalHW[NRAMSESTYPE];
    int         nstarTotal;
    RAMSESFLOAT mass[NRAMSESTYPE];
    double      Omegam, OmegaLambda, Omegak, Omegab, HubbleParam;
    double      aexp, time;
    double      scale_l,scale_d,scale_t;
    int         num_files;

    ///constructor
    RAMSES_Header() {
        for (int k=0;k<NRAMSESTYPE;k++) npartTotal[k]=npartTotalHW[k]=0;
    }
};
//@}

int RAMSES_fortran_read(fstream &, int &);
int RAMSES_fortran_read(fstream &, RAMSESFLOAT &);
int RAMSES_fortran_read(fstream &, int*);
int RAMSES_fortran_read(fstream &, RAMSESFLOAT *);
int RAMSES_fortran_read(fstream &, RAMSESIDTYPE *);
int RAMSES_fortran_skip(fstream &, int nskips=1);

/// \name Get the number of particles in the ramses files
//@{
Int_t RAMSES_get_nbodies(char *fname, int ptype, Options &opt);
//@}

#endif 
