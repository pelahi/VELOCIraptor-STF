/*! \file ui.cxx
 *  \brief user interface
 */

#include "halomergertree.h"

///routine to get arguments from command line
/// \todo alter interface as now need to be able to specify only smdata file (no grid, res, normalized res) and functionality to specify eps, background fof search, etc
void GetArgs(int argc, char *argv[], Options &opt)
{
    int option;
    int NumArgs = 0;
    while ((option = getopt(argc, argv, ":i:s:n:o:C:c:S:I:N:B:F:M:H:D:O:T:")) != EOF)
    {
        switch(option)
        {
            case 'i': 
                opt.fname = optarg;
                NumArgs += 2;
                break;
            case 's': 
                opt.numsnapshots = atoi(optarg);
                NumArgs += 2;
                break;
            case 'n': 
                opt.NumPart = atol(optarg);
                NumArgs += 2;
                break;
            case 'o': 
                opt.outname = optarg;
                NumArgs += 2;
                break;
            case 'C': 
                opt.matchtype = atoi(optarg);
                NumArgs += 2;
                break;
            case 'c': 
                opt.icatalog = atoi(optarg);
                NumArgs += 2;
                break;
            case 'D': 
                opt.idcorrectflag = atoi(optarg);
                NumArgs += 2;
                break;
            case 'S': 
                opt.mlsig = atof(optarg);
                NumArgs += 2;
                break;
            case 'I': 
                opt.ioformat = atoi(optarg);
                NumArgs += 2;
                break;
            case 'O': 
                opt.outputformat = atoi(optarg);
                NumArgs += 2;
                break;
            case 'N': 
                opt.nmpifiles = atoi(optarg);
                NumArgs += 2;
                break;
            case 'B': 
                opt.ibinary = atoi(optarg);
                NumArgs += 2;
                break;
            case 'F': 
                opt.ifield = atoi(optarg);
                NumArgs += 2;
                break;
            case 'M': 
                opt.imapping = atoi(optarg);
                NumArgs += 2;
                break;
            case 'H': 
                opt.snapshotvaloffset = atoi(optarg);
                NumArgs += 2;
                break;
            case 'T': 
                opt.itypematch = atoi(optarg);
                NumArgs += 2;
                break;
            case '?':
                usage();
        }
    }
    if (opt.fname==NULL||opt.outname==NULL){
#ifdef USEMPI
    if (ThisTask==0)
#endif
        cerr<<"Must provide input and output file names\n";
#ifdef USEMPI
        MPI_Finalize();
#endif
        exit(8);
    }

    if (opt.imapping==DSIMPLEMAP) opt.mappingfunc=simplemap;
    //else if (opt.imapping==???) opt.mappingfunc=???;
}

///Outputs the usage to stdout 
void usage(void)
{
    Options opt;
#ifdef USEMPI
    if (ThisTask==0) {
#endif
    cerr<<"USAGE:\n";
    cerr<<"\n";
    cerr<<"-i <file containing filelist>\n";
    cerr<<"-s <number of files/snapshots>\n";
    cerr<<"-n <number of particles>\n";
    cerr<<"-o <output filename>\n";
    cerr<<"-O <output format>\n";
    cerr<<"-C <cross correlation function type to identify main progenitor ("<<opt.matchtype<<" ["<<NsharedN1N2<<" "<<Nshared<<"])\n";
    cerr<<"-c <produce cross catalog match (0/1) default ("<<opt.icatalog<<")\n";
    cerr<<"-D <adjust IDs for nIFTY cross catalogs across simulations ("<<opt.idcorrectflag<<")\n";
    cerr<<"-T <type of particles to cross correlate ("<<opt.itypematch<<" ["<<ALLTYPEMATCH<<" is all particle types, ";
    cerr<<DMTYPEMATCH<<" is DM particle types, ";
    cerr<<GASTYPEMATCH<<" is GAS particle types, ";
    cerr<<STARTYPEMATCH<<" is STAR particle types, ";
    cerr<<DMGASTYPEMATCH<<" is both DM and GAS particle types, ";
    cerr<<",])\n";
    cerr<<"-S <significance of cross match relative to Poisson noise ("<<opt.mlsig<<")\n";
    cerr<<"-I <Input format ("<<opt.ioformat<<" [Sussing "<<DSUSSING<<", normal cat "<<DCATALOG<<", nIFTY "<<DNIFTY<<" ])\n";
    cerr<<"-N <if output is split between multiple files due to mpi, number of files written ("<<opt.nmpifiles<<")>\n";
    cerr<<"-B <binary or ascii format ("<<opt.ibinary<<")>\n";
    cerr<<"-F <field objects in separate file ("<<opt.ifield<<")>\n";
    cerr<<"-M <Mapping of particle ids to index ("<<opt.imapping<<" [ no maping "<<DNOMAP<<", simple mapping "<<DSIMPLEMAP<<"])\n";
#ifdef USEMPI
    }
    MPI_Finalize();
#endif
    exit(1);
}

