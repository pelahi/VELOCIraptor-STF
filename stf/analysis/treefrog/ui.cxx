/*! \file ui.cxx
 *  \brief user interface
 */

#include "TreeFrog.h"

///routine to get arguments from command line
void GetArgs(int argc, char *argv[], Options &opt)
{
    int option;
    int NumArgs = 0;
    while ((option = getopt(argc, argv, ":i:s:t:n:f:p:o:C:c:S:I:N:B:F:M:H:h:D:O:T:v:m:d:z:Z:a:x:u:")) != EOF)
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
                opt.MaxIDValue = atol(optarg);
                NumArgs += 2;
                break;
            case 'f':
                opt.particle_frac = atof(optarg);
                NumArgs += 2;
                break;
            case 'p':
                opt.min_numpart = atoi(optarg);
                NumArgs += 2;
                break;
            case 'o':
                opt.outname = optarg;
                NumArgs += 2;
                break;
            case 'C':
                opt.imerittype = atoi(optarg);
                NumArgs += 2;
                break;
            case 'x':
                opt.imultsteplinkcrit= atoi(optarg);
                NumArgs += 2;
                break;
            case 'u':
                opt.iopttemporalmerittype= atoi(optarg);
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
            case 'a':
                opt.meritlimit = atof(optarg);
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
            case 'z':
                opt.ndesiredmpithreads = atoi(optarg);
                NumArgs += 2;
                break;
            case 'Z':
                opt.iwriteparallel = atoi(optarg);
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
            MPI_Abort(MPI_COMM_WORLD,8);
#else
            exit(8);
#endif
    }

    if (opt.imapping==DSIMPLEMAP) opt.mappingfunc=simplemap;
    if (opt.numsnapshots<2){
#ifdef USEMPI
    if (ThisTask==0)
#endif
        cerr<<"Number of snapshots must be >=2\n";
#ifdef USEMPI
            MPI_Abort(MPI_COMM_WORLD,8);
#else
            exit(8);
#endif
    }
    if (opt.numsteps<1){
#ifdef USEMPI
    if (ThisTask==0)
#endif
        cerr<<"Number of steps over which to produce links must be >=1\n";
#ifdef USEMPI
            MPI_Abort(MPI_COMM_WORLD,8);
#else
            exit(8);
#endif
    }
    if (opt.icatalog==DCROSSCAT) {
        if (opt.numsnapshots>2) {cerr<<"Cross catalog, yet more than two snapshots compared, reseting and only comparing two"<<endl;opt.numsnapshots=2;}
        if (opt.numsteps>1) {cerr<<"Cross catalog, yet asking to use more than a single step when linking, reseting and only linking across one "<<endl;opt.numsteps=1;}
    }
    //else if (opt.imapping==???) opt.mappingfunc=???;
    opt.description=(char*)"VELOCIraptor halo merger tree constructed by identifying the main progenitor with the highest value of ";
    if(opt.imerittype==NsharedN1N2)      opt.description+=(char*)"Nshared^2/Nh/Np |";
    else if(opt.imerittype==NsharedN1)   opt.description+=(char*)"Nshared/Nh | ";
    else if(opt.imerittype==Nshared)     opt.description+=(char*)"Nshared |";
    else if (opt.imerittype==Nsharedcombo) opt.description=(char*)"Nshared/Nh+(Nshared^2/Nh/Np) so as to weight progenitors that contribute similar amounts by how much of their mass contributes to the new object | ";
    opt.description=(char*)"Optimal temporal merits are set by  ";
    if(opt.iopttemporalmerittype==GENERALIZEDMERITTIME)  opt.description+=(char*)"a generalized temporal merit taking into account average time evolution |";
    else if(opt.iopttemporalmerittype==GENERALIZEDMERITTIMEPROGEN)  opt.description+=(char*)"a generalized temporal merit taking into account average time evolution and maximise the ranking of the progenitor so that links always point to primary progen/descen |";
    opt.description+=(char*)"Tree built using ";
    opt.description+=static_cast<ostringstream*>( &(ostringstream() << opt.numsteps) )->str();
    opt.description+=(char*)" temporal steps | ";
    opt.description+=(char*)"Particle types for matching limited to ";
    if (opt.itypematch==ALLTYPEMATCH) opt.description+=(char*)" all |";
    else {opt.description+=(char*)" part type ";opt.description+=static_cast<ostringstream*>( &(ostringstream() << opt.itypematch) )->str();}
    opt.description+=(char*)" | ";
    if (opt.particle_frac<1 && opt.particle_frac>0) {
    opt.description+=(char*)" Weighted merit with ";opt.description+=static_cast<ostringstream*>( &(ostringstream() << opt.particle_frac) )->str();
    opt.description+=(char*)" fraction of most bound particles with min particle num of  ";opt.description+=static_cast<ostringstream*>( &(ostringstream() << opt.min_numpart) )->str();
    opt.description+=(char*)" | ";
    }
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
    cerr<<" ========================= "<<endl;
    cerr<<"-i <file containing filelist>\n";
    cerr<<"-I <Input format ("<<opt.ioformat<<" [Sussing "<<DSUSSING<<", normal velociraptor catalog "<<DCATALOG<<", nIFTY "<<DNIFTY<<", Void "<<DVOID<<" ])\n";
    cerr<<"-s <number of files/snapshots>\n";
    cerr<<"-N <if output is split between multiple files due to mpi, number of files written ("<<opt.nmpifiles<<")>\n";
    cerr<<"-c <produce cross catalog match (0 halo tree ,1 cross catalog ,2 full graph) default ("<<opt.icatalog<<")\n";
    cerr<<"-o <output filename>\n";
    cerr<<"-O <output format>\n";
    cerr<<" ========================= "<<endl<<endl;

    cerr<<" ========================= "<<endl;
    cerr<<" For normal catalog produced by velociraptor "<<endl;
    cerr<<" ========================= "<<endl;
    cerr<<"-B <input format for (binary 1, hdf 2, or ascii 0) ("<<opt.ibinary<<")>\n";
    cerr<<"-F <field objects in separate file ("<<opt.ifield<<")>\n";
    cerr<<" ========================= "<<endl<<endl;

    cerr<<" ========================= "<<endl;
    cerr<<" Cross matching options "<<endl;
    cerr<<" ========================= "<<endl;
    cerr<<"-n <Max ID value of particles [Must specify] ("<<opt.MaxIDValue<<")>\n";
    cerr<<"-t <number of steps integrated over to find links ("<<opt.numsteps<<")\n";
    cerr<<"-x <criteria for when to keep searching for new links if multiple steps invoked ("<<opt.imultsteplinkcrit<<"). Possibilities are :";
    cerr<<MSLCMISSING<<" Only missing ,";
    cerr<<MSLCMERIT<<" Missing & low merit given by merit limit, ";
    cerr<<endl;
    cerr<<"-C <cross correlation function type to identify main progenitor ("<<opt.imerittype<<" ["<<NsharedN1N2<<" "<<Nshared<<"])\n";
    cerr<<"-T <type of particles to cross correlate ("<<opt.itypematch<<" ["<<ALLTYPEMATCH<<" is all particle types, ";
    cerr<<DMTYPEMATCH<<" is DM particle types, ";
    cerr<<GASTYPEMATCH<<" is GAS particle types, ";
    cerr<<STARTYPEMATCH<<" is STAR particle types, ";
    cerr<<DMGASTYPEMATCH<<" is both DM and GAS particle types, ";
    cerr<<",])\n";
    cerr<<"-S <significance of cross match relative to Poisson noise ("<<opt.mlsig<<")\n";
    cerr<<"-a <merit limit to search for new links ("<<opt.meritlimit<<")\n";
    cerr<<"-f <fraction of particles to use to calculate weigted merit>\n";
    cerr<<"-p <minimum number of particles used in weighted merit>\n";
    cerr<<" ========================= "<<endl<<endl;

    cerr<<" ========================= "<<endl;
    cerr<<" ID to index mapping if required"<<endl;
    cerr<<" ========================= "<<endl;
    cerr<<"-D <adjust particle IDs for nIFTY cross catalogs across simulations ("<<opt.idcorrectflag<<")\n";
    cerr<<"-M <Mapping of particle ids to index ("<<opt.imapping<<" [ no maping "<<DNOMAP<<", simple mapping "<<DSIMPLEMAP<<", computational expensive but memory efficient adaptive map "<<DMEMEFFICIENTMAP<<"])\n";
    cerr<<" ========================= "<<endl<<endl;

    cerr<<" ========================= "<<endl;
    cerr<<" Ouptut options for halo ID values"<<endl;
    cerr<<" ========================= "<<endl;
    cerr<<"-h <adjust Halo ID values stored in group catalog, useful for matching these values to those stored in .properties files produced by the halo finder. output is ID+(snap+snapshotvaloffset)*haloIDval ("<<opt.haloidval<<")\n";
    cerr<<"-d <adjust Halo ID by this offset value ("<<opt.haloidoffset<<")\n";
    cerr<<"-H <offset snapshot values by this number ("<<opt.snapshotvaloffset<<")>\n";
    cerr<<" ========================= "<<endl<<endl;

    cerr<<" ========================= "<<endl;
    cerr<<" Other options"<<endl;
    cerr<<" ========================= "<<endl;
    cerr<<"-v <verbose flag 1/0 ("<<opt.iverbose<<")\n";
    cerr<<" ========================= "<<endl<<endl;

#ifdef USEMPI
    cerr<<" ========================= "<<endl;
    cerr<<" For mpi load balancing "<<endl;
    cerr<<" ========================= "<<endl;
    cerr<<"-m <number of items per mpi thead, use for load balacing. If 0, based on input ("<<opt.numpermpi<<")\n";
    cerr<<"-z <number of mpi theads used to calculate load balacing. If >0 this used with one actual mpi thread but determines load balancing based on desired number of mpi threads. Write load balancing file and terminates. If 0 (default) normal operation \n";
    cerr<<"-Z <whether to write output in parallel (0/1). \n";
    cerr<<" ========================= "<<endl<<endl;
#endif
#ifdef USEMPI
    }
    MPI_Finalize();
#endif
    exit(1);
}

/*!
    \page args Command line options
    The following commands are accepted as command line arguments (more info can be found in \ref Options struct and \ref ui.cxx for user interface or the sample configuration file
in the examples directory). \n \n

    \section ioconfig Input/Output related.
    See \ref io.cxx, \ref stfio.cxx, and \ref otherio.cxx for implementation of the code, \ref allvars.h for definitions
    \arg \b \e -i < file containing filelist >
    \arg \b \e -s < number of files/snapshots to be processed >
    \arg \b \e -I < Input format, can be \ref DSUSSING (ascii) \ref DCATALOG (VELOCIraptor output which can be ascii, binary or hdf see -B),  \ref DNIFTY (ascii), \ref DVOID (ascii) >
    \arg \b \e -N < 0/1 if output is split between multiple files due to mpi, number of files written, useful for reading VELOCIraptor output >
    \arg \b \e -B < for VELOCIraptor catalogs can be \ref INASCII (ascii) \ref INBINARY (binary) or \ref INHDF (hdf5) >
    \arg \b \e -F < 0/1 field objects in VELOCIraptor output in separate file >
    \arg \b \e -o < output filename >
    \arg \b \e -O < output type, whether merger tree, full graph, whether to output merit of each match >
    \arg \b \e -c < produce cross catalog match (0 halo tree ,1 cross catalog where code acts like particle correlator, matching analogues between simulations >

    \section linkingconfig Linking options
    \arg \b \e -t < number of steps over to find links >
    \arg \b \e -n < maximum id of particles, used to allocate an array of this size so that ids are mapped to an index, allows quick check of what groups a particle belongs to>
    \arg \b \e -C < cross correlation function type to identify main progenitor (see \ref Options.imerittype, which can be  \ref NsharedN1N2 standard merit function \f[ \mathcal{M}=N_{\rm sh}^2/N_1/N_2 \f], \ref Nshared >
    \arg \b \e -T < type of particles to cross correlate (\ref opt.itypematch, which can be \ref ALLTYPEMATCH is all particle types, \ref DMTYPEMATCH is DM particle types
    \ref GASTYPEMATCH is GAS particle types \ref STARTYPEMATCH is STAR particle types, \ref DMGASTYPEMATCH is both DM and GAS particle types >
    \arg \b \e -S < significance of cross match relative to Poisson noise >
    \arg \b \e -a < merit limit below which will search for links at earlier snapshot (as if no match is found)>

    \section outputconfig Output options
    \arg \b \e -H < offset snapshot values by this number >
    \arg \b \e -h < adjust Halo ID values stored in group catalog, useful for matching these values to those stored in .properties files produced by the halo finder. output is ID+(snap+snapshotvaloffset)*haloIDval >
    \arg \b \e -d < adjust Halo ID by this offset value >

    \section otherconfig Other options
    \arg \b \e -D < adjust particle IDs for nIFTY cross catalogs across simulations >
    \arg \b \e -M < Mapping of particle ids to index  (\ref opt.imapping with \ref DNOMAP no mapping of ids to indices, all others must be implemented >
    \arg \b \e -v < verbose flag 1/0 >
    \section mpiconfig MPI options
    \arg \b \e -m < number of items per mpi thead, use for load balacing. Use 0 if no mpi used when building halo catalog >
    \arg \b \e -z < number of mpi theads used to calculate load balacing. If >0 this used with one actual mpi thread but determines load balancing based on desired number of mpi threads. Write load balancing file and terminates. If 0 (default) normal operation >
    \arg \b \e -Z < whether to write output in parallel (0/1) >

*/
