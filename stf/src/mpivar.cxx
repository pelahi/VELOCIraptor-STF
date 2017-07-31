/*! \file mpivar.cxx
 *  \brief this file contains external variables listed in mpivar
 */

#ifdef USEMPI

//-- For MPI

#include "stf.h"

#include "gadgetitems.h"
#include "tipsy_structs.h"
#include "endianutils.h"

///\name Variabiles for mpi tasks
/// defined in \ref mpivar.h
//@{

int ThisTask, NProcs;
Int_t Ntotal, NExport, NImport, Nlocal, Nlocaltot,Ngridtotal,Ngridlocal;
Int_t Ntotalbaryon[NBARYONTYPES], Nlocalbaryon[NBARYONTYPES];
Int_t Nmemlocal,Noldlocal,Nmemlocalbaryon;

Double_t mpi_xlim[3][2],mpi_dxsplit[3];
int mpi_nxsplit[3],mpi_ideltax[3];
Double_t mpi_period;
MPI_Domain *mpi_domain;
Int_t *mpi_nlocal,*mpi_nsend,*mpi_idlist;
short_mpi_t *mpi_foftask;
Int_t *mpi_ngroups, *mpi_pfof, *mpi_indexlist;
int *mpi_part_send_domain;

Int_t mpi_maxgid,mpi_gidoffset;
fofdata_in *FoFDataIn, *FoFDataGet;
fofid_in *FoFGroupDataLocal, *FoFGroupDataExport;
nndata_in *NNDataIn, *NNDataGet;
Particle *PartDataIn, *PartDataGet;
Particle *mpi_Part1=NULL, *mpi_Part2=NULL;

Int_t MinNumMPI,MinNumOld;

GridCell *mpi_grid;
Coordinate *mpi_gvel;
Matrix *mpi_gveldisp;

//@}



#endif
