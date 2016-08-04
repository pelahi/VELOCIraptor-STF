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
    while ((option = getopt(argc, argv, ":i:s:t:n:o:C:c:S:I:N:B:F:M:H:h:D:O:T:v:m:d:")) != EOF)
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
            case 't': 
                opt.numsteps = atoi(optarg);
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
            case 'h': 
                opt.haloidval= atol(optarg);
                NumArgs += 2;
                break;
            case 'd': 
                opt.haloidoffset= atoi(optarg);
                NumArgs += 2;
                break;
            case 'T': 
                opt.itypematch = atoi(optarg);
                NumArgs += 2;
                break;
            case 'v': 
                opt.iverbose = atoi(optarg);
                NumArgs += 2;
                break;
#ifdef USEMPI
            case 'm': 
                opt.numpermpi = atol(optarg);
                NumArgs += 2;
                break;
#endif
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
    if (opt.numsnapshots<2){
#ifdef USEMPI
    if (ThisTask==0)
#endif
        cerr<<"Number of snapshots must be >=2\n";
#ifdef USEMPI
        MPI_Finalize();
#endif
        exit(8);
    }
    if (opt.numsteps<1){
#ifdef USEMPI
    if (ThisTask==0)
#endif
        cerr<<"Number of steps over which to produce links must be >=1\n";
#ifdef USEMPI
        MPI_Finalize();
#endif
        exit(8);
    }
    if (opt.icatalog==DCROSSCAT) {
        if (opt.numsnapshots>2) {cerr<<"Cross catalog, yet more than two snapshots compared, reseting and only comparing two"<<endl;opt.numsnapshots=2;}
        if (opt.numsteps>1) {cerr<<"Cross catalog, yet asking to use more than a single step when linking, reseting and only linking across one "<<endl;opt.numsteps=1;}
    }
    //else if (opt.imapping==???) opt.mappingfunc=???;
    opt.description=(char*)"VELOCIraptor halo merger tree constructed by identifying the main progenitor with the highest value of ";
    if(opt.matchtype==NsharedN1N2)      opt.description+=(char*)"Nshared^2/Nh/Np |";
    else if(opt.matchtype==NsharedN1)   opt.description+=(char*)"Nshared/Nh | ";
    else if(opt.matchtype==Nshared)     opt.description+=(char*)"Nshared |";
    else if (opt.matchtype==Nsharedcombo) opt.description=(char*)"Nshared/Nh+(Nshared^2/Nh/Np) so as to weight progenitors that contribute similar amounts by how much of their mass contributes to the new object | ";
    opt.description+=(char*)"Tree built using ";
    opt.description+=static_cast<ostringstream*>( &(ostringstream() << opt.numsteps) )->str();
    opt.description+=(char*)" temporal steps | ";
    opt.description+=(char*)"Particle types for matching limited to ";
    if (opt.itypematch==ALLTYPEMATCH) opt.description+=(char*)" all |";
    else {opt.description+=(char*)" part type ";opt.description+=static_cast<ostringstream*>( &(ostringstream() << opt.itypematch) )->str();}
    opt.description+=(char*)" | ";
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
    cerr<<"-t <number of steps integrated over to find links ("<<opt.numsteps<<")\n";
    cerr<<"-n <number of particles>\n";
    cerr<<"-o <output filename>\n";
    cerr<<"-O <output format>\n";
    cerr<<"-C <cross correlation function type to identify main progenitor ("<<opt.matchtype<<" ["<<NsharedN1N2<<" "<<Nshared<<"])\n";
    cerr<<"-c <produce cross catalog match (0 halo tree ,1 cross catalog ,2 full graph) default ("<<opt.icatalog<<")\n";
    cerr<<"-T <type of particles to cross correlate ("<<opt.itypematch<<" ["<<ALLTYPEMATCH<<" is all particle types, ";
    cerr<<DMTYPEMATCH<<" is DM particle types, ";
    cerr<<GASTYPEMATCH<<" is GAS particle types, ";
    cerr<<STARTYPEMATCH<<" is STAR particle types, ";
    cerr<<DMGASTYPEMATCH<<" is both DM and GAS particle types, ";
    cerr<<",])\n";
    cerr<<"-S <significance of cross match relative to Poisson noise ("<<opt.mlsig<<")\n";
    cerr<<"-I <Input format ("<<opt.ioformat<<" [Sussing "<<DSUSSING<<", normal cat "<<DCATALOG<<", nIFTY "<<DNIFTY<<", Void "<<DVOID<<" ])\n";
    cerr<<"-N <if output is split between multiple files due to mpi, number of files written ("<<opt.nmpifiles<<")>\n";
    cerr<<"-B <binary or ascii format ("<<opt.ibinary<<")>\n";
    cerr<<"-F <field objects in separate file ("<<opt.ifield<<")>\n";
    cerr<<"-H <offset snapshot values by this number ("<<opt.snapshotvaloffset<<")>\n";
    cerr<<"-h <adjust Halo ID values stored in group catalog, useful for matching these values to those stored in .properties files produced by the halo finder. output is ID+(snap+snapshotvaloffset)*haloIDval ("<<opt.haloidval<<")\n";
    cerr<<"-d <adjust Halo ID by this offset value ("<<opt.haloidoffset<<")\n";

    cerr<<"-D <adjust particle IDs for nIFTY cross catalogs across simulations ("<<opt.idcorrectflag<<")\n";
    cerr<<"-M <Mapping of particle ids to index ("<<opt.imapping<<" [ no maping "<<DNOMAP<<", simple mapping "<<DSIMPLEMAP<<"])\n";
    cerr<<"-v <verbose flag 1/0 ("<<opt.iverbose<<")\n";
#ifdef USEMPI
    cerr<<"-m <number of items per mpi thead, use for load balacing. If 0, based on input ("<<opt.numpermpi<<")\n";
    
#endif
#ifdef USEMPI
    }
    MPI_Finalize();
#endif
    exit(1);
}

