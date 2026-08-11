/*! \file newramsesitems.h
 *  \brief definitions and routines for reading RAMSES snapshots stored in the NewCluster HDF5 format
 *
 * Unlike the original multi-file binary RAMSES format (see \ref ramsesitems.h and \ref ramsesio.cxx),
 * each snapshot in this format is stored as a single part_%s.h5 file below opt.fname (dark matter,
 * star, sink, cloud and tracer particles), named using the same snapshot string passed via
 * opt.ramsessnapname (the -t option). Gas cells live in a companion cell_%s.h5 file and are not yet
 * read here.
 *
 * Each particle type sits under its own group ("/dm", "/star", ...) with a "data" dataset holding one
 * HDF5 compound (struct-of-arrays-free, i.e. array-of-structs) record per particle, plus
 * "hilbert_boundary"/"chunk_boundary"/"level_boundary" datasets describing an internal spatial chunking
 * of the data that is not used by this reader (the whole per-type dataset is instead read out in
 * contiguous row-range chunks).
 */

#ifndef NEWRAMSESITEMS_H
#define NEWRAMSESITEMS_H

///pulls in the RAMSESFLOAT typedef used by the MPI domain-decomposition helpers below
#include "ramsesitems.h"

///number of particle records read from a compound dataset in one hyperslab call
#define NEWRAMSESCHUNKSIZE 5000000

/// \name Get the number of particles of interest stored in the new (HDF5) ramses format
//@{
Int_t NewRAMSES_get_nbodies(char *fname, int ptype, Options &opt);
//@}

///Reads dark matter and star particles from the new (HDF5) ramses format. Mirrors the interface of
///\ref ReadRamses so it can be used as a drop-in alternative selected via opt.inewramsesio.
void ReadNewRamses(Options &opt, vector<Particle> &Part, const Int_t nbodies, Particle *&Pbaryons, Int_t nbaryons=0);

#ifdef USEMPI
/// \name MPI domain-decomposition pre-pass helpers, used by mpiramsesio.cxx (see MPINumInDomainRAMSES and
/// MPIDomainDecompositionWithTree) as the HDF5-format equivalent of their direct binary-file reads.
//@{
///counts, into Nbuf[ibuf] (and Nbaryonbuf[ibuf] when opt.partsearchtype==PSTDARK && opt.iBaryonSearch),
///the number of dm/star particles of interest that fall into each MPI domain. Caller MPI_Allreduce's
///Nbuf/Nbaryonbuf afterwards.
void NewRAMSES_CountNumInDomain(Options &opt, Int_t *Nbuf, Int_t *Nbaryonbuf);

///reads dm/star particles (position, mass, type; velocity left at 0) into Part_mpi starting at
///count_mpi, and records the same particles' raw code-unit positions/type into xtempall (sized
///3*nbodies) / famtempall (sized nbodies) at the same index, for use by the tree-based MPI domain
///decomposition in MPIDomainDecompositionWithTree. Also returns the length scale (code units to kpc)
///used, needed by the caller to convert the tree's node boundaries back to code units. PSTDARK with
///iBaryonSearch is not handled (mirrors the unimplemented case already present in
///MPIDomainDecompositionWithTree).
void NewRAMSES_ReadForDomainTree(Options &opt, vector<Particle> &Part_mpi, Int_t nbodies, Int_t &count_mpi, RAMSESFLOAT *xtempall, int *famtempall, Double_t &lscale);
//@}
#endif

#endif
