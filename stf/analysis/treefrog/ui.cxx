/*! \file ui.cxx
 *  \brief user interface
 */

#include "TreeFrog.h"

///Read parameters from a parameter file. For list of currently implemented options see \ref configopt
///\todo still more parameters that can be adjusted
void GetParamFile(Options &opt)
{
    string line,sep="=";
    string tag,val;
    char buff[1024],*pbuff,tbuff[1024],vbuff[1024],fname[1024];
    fstream paramfile,cfgfile;
    if (!FileExists(opt.configname)){
            cerr<<"Config file does not exist or can't be read, terminating"<<endl;
#ifdef USEMPI
            MPI_Abort(MPI_COMM_WORLD,9);
#else
            exit(9);
#endif
    }
    paramfile.open(opt.configname, ios::in);
    unsigned j,k;
    if (paramfile.is_open())
    {
        while (paramfile.good()){
            getline(paramfile,line);
            //if line is not commented out or empty
            if (line[0]!='#'&&line.length()!=0) {
                if (j=line.find(sep)){
                    //clean up string
                    tag=line.substr(0,j);
                    strcpy(buff, tag.c_str());
                    pbuff=strtok(buff," ");
                    strcpy(tbuff, pbuff);
                    val=line.substr(j+1,line.length()-(j+1));
                    strcpy(buff, val.c_str());
                    pbuff=strtok(buff," ");
                    if (pbuff==NULL) continue;
                    strcpy(vbuff, pbuff);
                    //config search type
                    if (strcmp(tbuff, "Tree_direction")==0) {
                        opt.isearchdirection = atoi(vbuff);
                    }
                    else if (strcmp(tbuff, "Particle_type_criterion")==0) {
                        opt.itypematch = atoi(vbuff);
                    }
                    else if (strcmp(tbuff, "Merit_type")==0) {
                        opt.imerittype = atoi(vbuff);
                        opt.idefaultvalues = 0;
                    }
                    else if (strcmp(tbuff, "Multistep_criterion")==0) {
                        opt.imultsteplinkcrit = atoi(vbuff);
                        opt.idefaultvalues = 0;
                    }
                    else if (strcmp(tbuff, "Optimal_temporal_merit_criterion")==0) {
                        opt.iopttemporalmerittype = atoi(vbuff);
                        opt.idefaultvalues = 0;
                    }

                    //config search parameters
                    else if (strcmp(tbuff, "Number_of_linking_steps")==0)
                        opt.numsteps = atoi(vbuff);
                    else if (strcmp(tbuff, "Shared_particle_signal_to_noise_limit")==0) {
                        opt.mlsig = atof(vbuff);
                        opt.idefaultvalues = 0;
                    }
                    else if (strcmp(tbuff, "Merit_limit_continuing_search")==0) {
                        opt.meritlimit = atof(vbuff);
                        opt.idefaultvalues = 0;
                    }
                    else if (strcmp(tbuff, "Particle_core_fraction")==0) {
                        opt.particle_frac = atof(vbuff);
                        opt.idefaultvalues = 0;
                    }
                    else if (strcmp(tbuff, "Particle_core_min_num")==0) {
                        opt.min_numpart = atoi(vbuff);
                        opt.idefaultvalues = 0;
                    }
                    else if (strcmp(tbuff, "Particle_core_max_num")==0){
                        opt.max_numpart = atoi(vbuff);
                        opt.idefaultvalues = 0;
                    }
                }
            }
        }
        paramfile.close();
    }
}

///routine to get arguments from command line
void GetArgs(int argc, char *argv[], Options &opt)
{
    int option;
    int NumArgs = 0;
    while ((option = getopt(argc, argv, ":i:s:I:N:B:F:o:O:d:T:D:M:X:E:U:C:l:m:n:t:f:p:b:a:h:H:g:v:y:z:Z:")) != EOF)
    {
        switch(option)
        {
            ///input related options such as input name, number of snaps to process, input format, outputformat
            case 'i':
                opt.fname = optarg;
                NumArgs += 2;
                break;
            case 's':
                opt.numsnapshots = atoi(optarg);
                NumArgs += 2;
                break;
            case 'I':
                opt.ioformat = atoi(optarg);
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
            case 'o':
                opt.outname = optarg;
                NumArgs += 2;
                break;
            case 'O':
                opt.outputformat = atoi(optarg);
                NumArgs += 2;
                break;
            case 'd':
                opt.outdataformat = atoi(optarg);
                NumArgs += 2;
                break;

            //general operations
            case 'T':
                opt.itypematch = atoi(optarg);
                NumArgs += 2;
                break;
            case 'D':
                NumArgs += 2;
                opt.isearchdirection = atoi(optarg);
                break;
            case 'C':
                opt.icatalog = atoi(optarg);
                opt.idefaultvalues = 0;
                NumArgs += 2;
                break;

            //cross matching options
            case 'M':
                opt.imerittype = atoi(optarg);
                opt.idefaultvalues = 0;
                NumArgs += 2;
                break;
            case 'X':
                opt.imultsteplinkcrit= atoi(optarg);
                opt.idefaultvalues = 0;
                NumArgs += 2;
                break;
            case 'E':
                opt.icorematchtype= atoi(optarg);
                opt.idefaultvalues = 0;
                NumArgs += 2;
                break;
            case 'U':
                opt.iopttemporalmerittype= atoi(optarg);
                opt.idefaultvalues = 0;
                NumArgs += 2;
                break;

            //id index mapping
            case 'l':
                opt.idcorrectflag = atoi(optarg);
                NumArgs += 2;
                break;
            case 'm':
                opt.imapping = atoi(optarg);
                NumArgs += 2;
                break;
            case 'n':
                opt.MaxIDValue = atol(optarg);
                NumArgs += 2;
                break;

            //configuration of operation of tree
            case 't':
                opt.numsteps = atoi(optarg);
                NumArgs += 2;
                break;
            case 'f':
                opt.particle_frac = atof(optarg);
                opt.idefaultvalues = 0;
                NumArgs += 2;
                break;
            case 'p':
                opt.min_numpart = atoi(optarg);
                opt.idefaultvalues = 0;
                NumArgs += 2;
                break;
            case 'b':
                opt.mlsig = atof(optarg);
                opt.idefaultvalues = 0;
                NumArgs += 2;
                break;
            case 'a':
                opt.meritlimit = atof(optarg);
                opt.idefaultvalues = 0;
                NumArgs += 2;
                break;

            //output offsets
            case 'H':
                opt.snapshotvaloffset = atoi(optarg);
                NumArgs += 2;
                break;
            case 'h':
                opt.haloidval= atol(optarg);
                NumArgs += 2;
                break;
            case 'g':
                opt.haloidoffset= atoi(optarg);
                NumArgs += 2;
                break;
            //other options
            case 'v':
                opt.iverbose = atoi(optarg);
                NumArgs += 2;
                break;
#ifdef USEMPI
            //mpi related
            case 'y':
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
    /*
    if(configflag){
        cout<<"Reading config file"<<endl;
        GetParamFile(opt);
    }
    else {
        cout<<"NO CONFIG FILE PASSED! Using default values"<<endl;
    }
    */
#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    ConfigCheck(opt);
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
    cerr<<" Basic input/output commands "<<endl;
    cerr<<" ========================= "<<endl;
    cerr<<"-i <file containing filelist>\n";
    cerr<<"-I <Input format ("<<opt.ioformat<<" [Sussing "<<DSUSSING<<", normal velociraptor catalog "<<DCATALOG<<", nIFTY "<<DNIFTY<<", Void "<<DVOID<<" ])\n";
    cerr<<"-s <number of files/snapshots>\n";
    cerr<<"-c <produce cross catalog match (0 halo tree ,1 cross catalog ,2 full graph) default ("<<opt.icatalog<<")\n";
    cerr<<"-o <output filename>\n";
    cerr<<"-O <output format, ASCII, HDF5 ("<<OUTASCII<<","<<" "<<OUTHDF<<"), with default "<<opt.outputformat<<">\n";
    cerr<<"-d <output data, 0 for minimal, 1 for outputing merits as well, with default "<<opt.outdataformat<<">\n";
    cerr<<" ========================= "<<endl<<endl;

    cerr<<" ========================= "<<endl;
    cerr<<" input modifiers for normal catalog produced by velociraptor "<<endl;
    cerr<<" ========================= "<<endl;
    cerr<<"-B <input format for (binary 1, hdf 2, or ascii 0) ("<<opt.ibinary<<")>\n";
    cerr<<"-F <field objects in separate file ("<<opt.ifield<<")>\n";
    cerr<<"-N <if output is split between multiple files due to mpi, number of files written ("<<opt.nmpifiles<<")>\n";
    cerr<<" ========================= "<<endl<<endl;

    cerr<<" ========================= "<<endl;
    cerr<<" General tree construction options "<<endl;
    cerr<<" ========================= "<<endl;
    cerr<<"-D <Direction of tree. Default (" <<opt.isearchdirection<<") Possibilities are : \n";
    cerr<<'\t'<<SEARCHPROGEN<<" progenitor, \n";
    cerr<<'\t'<<SEARCHDESCEN<<" descendant, \n";
    cerr<<'\t'<<SEARCHALL<<" both directions \n";
    cerr<<"-T <type of particles to cross correlate. Default ("<<opt.itypematch<<"). Possibilities are :\n";
    cerr<<'\t'<<ALLTYPEMATCH<<" is all particle types, \n";
    cerr<<'\t'<<DMTYPEMATCH<<" is DM particle types, \n";
    cerr<<'\t'<<GASTYPEMATCH<<" is GAS particle types, \n";
    cerr<<'\t'<<STARTYPEMATCH<<" is STAR particle types, \n";
    cerr<<'\t'<<DMGASTYPEMATCH<<" is both DM and GAS particle types, \n";
    cerr<<"-M <cross correlation function type to identify main progenitor/descendant/link. Default ("<<opt.imerittype<<"). Possibilities are :\n";
    cerr<<'\t'<<NsharedN1N2<<" standard merit of Nshared^2/N1/N2, \n";
    cerr<<'\t'<<NsharedN1<<" fraction merit of Nshared/N1, \n";
    cerr<<"-X <criteria for when to keep searching for new links if multiple steps invoked ("<<opt.imultsteplinkcrit<<"). Possibilities are :\n";
    cerr<<'\t'<<MSLCMISSING<<" Only missing ,\n";
    cerr<<'\t'<<MSLCMERIT<<" Missing & low merit given by merit limit, \n";
    cerr<<'\t'<<MSLCPRIMARYPROGEN<<" Missing & or low ranking descendant when constructing descendant tree, \n";
    cerr<<'\t'<<MSLCMERITPRIMARYPROGEN<<" Missing & low merit & or low ranking descendant when constructing descendant tree, \n";
    cerr<<endl;
    cerr<<"-E <core matching criteria when searching  ("<<opt.icorematchtype<<"). Possibilities are :\n";
    cerr<<'\t'<<PARTLISTNOCORE<<" no core matching refinement ,\n";
    cerr<<'\t'<<PARTLISTCORE<<" for descendant tree, identify initial links using core of current halo and all particles with biased Nsh^2/N1/N2 merit of other haloes followed by core-core with unbiased Nsh^2/N1/N2 merit.";
    cerr<<" For progenitor tree, does core to all with unbiased Nsh^2/N1/N2 merit. \n";
    cerr<<'\t'<<PARTLISTCORECORE<<" for descendant and progenitor trees, does normal all to all particle match followed by a core-to-core match with unbiased merits.\n";
    cerr<<endl;
    cerr<<"-U <generalized temporal merit ("<<opt.iopttemporalmerittype<<").";
    cerr<<endl;
    cerr<<" ========================= "<<endl<<endl;

    cerr<<" ========================= "<<endl;
    cerr<<" Cross matching options "<<endl;
    cerr<<" ========================= "<<endl;
    cerr<<"-t <number of steps integrated over to find links ("<<opt.numsteps<<")\n";
    cerr<<"-b <significance of cross match relative to Poisson noise ("<<opt.mlsig<<")\n";
    cerr<<"-a <merit limit to search for new links ("<<opt.meritlimit<<")\n";
    cerr<<"-f <fraction of particles to use to calculate weigted merit>\n";
    cerr<<"-p <minimum number of particles used in weighted merit>\n";
    cerr<<" ========================= "<<endl<<endl;

    cerr<<" ========================= "<<endl;
    cerr<<" ID related options "<<endl;
    cerr<<" ========================= "<<endl;
    cerr<<"-n <Max ID value of particles [Must specify if not mapping ids to index] ("<<opt.MaxIDValue<<")>\n";
    cerr<<"-D <adjust particle IDs for nIFTY cross catalogs across simulations ("<<opt.idcorrectflag<<")\n";
    cerr<<"-m <Mapping of particle ids to index ("<<opt.imapping<<" [ no maping "<<DNOMAP<<", simple mapping "<<DSIMPLEMAP<<", computational expensive but memory efficient adaptive map "<<DMEMEFFICIENTMAP<<"])\n";
    cerr<<" ========================= "<<endl<<endl;

    cerr<<" ========================= "<<endl;
    cerr<<" Ouptut options for halo ID values"<<endl;
    cerr<<" ========================= "<<endl;
    cerr<<"-h <adjust Halo ID values stored in group catalog, useful for matching these values to those stored in .properties files produced by the halo finder. output is ID+(snap+snapshotvaloffset)*haloIDval ("<<opt.haloidval<<")\n";
    cerr<<"-g <adjust Halo ID by this offset value ("<<opt.haloidoffset<<")\n";
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
    cerr<<"-y <number of items per mpi thead, use for load balacing. If 0, based on input ("<<opt.numpermpi<<")\n";
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

//check to see if config options are acceptable.
inline void ConfigCheck(Options &opt)
{
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
    if (opt.idefaultvalues) {
        cout<<"Parameters specifying special matching criteria not passed by user, using default values"<<endl;
        if (opt.icatalog==DCROSSCAT){
            if(opt.isearchdirection==SEARCHALL) {
                opt.icorematchtype=PARTLISTNOCORE;
                opt.min_numpart=20;
                opt.particle_frac=0;
                opt.meritlimit=0;
                opt.imerittype=NsharedN1N2;
            }
        }
        else {
            if(opt.isearchdirection==SEARCHDESCEN) {
                opt.icorematchtype=PARTLISTCORE;
                opt.min_numpart=20;
                opt.particle_frac=0.1;
                opt.meritlimit=0.1;
                opt.imerittype=NsharedN1N2;
                opt.imultsteplinkcrit=MSLCPRIMARYPROGEN;
                opt.iopttemporalmerittype=GENERALIZEDMERITTIMEPROGEN;
            }
            else if (opt.isearchdirection==SEARCHPROGEN) {
                opt.icorematchtype=PARTLISTCORECORE;
                opt.min_numpart=20;
                opt.particle_frac=0.1;
                opt.meritlimit=0.1;
                opt.imerittype=NsharedN1N2;
                opt.imultsteplinkcrit=MSLCMERIT;
                opt.iopttemporalmerittype=GENERALIZEDMERITTIME;
            }
        }
    }
    if (opt.icatalog==DCROSSCAT) {
        if (opt.numsnapshots>2) {cerr<<"Cross catalog, yet more than two snapshots compared, reseting and only comparing two"<<endl;opt.numsnapshots=2;}
        if (opt.numsteps>1) {cerr<<"Cross catalog, yet asking to use more than a single step when linking, reseting and only linking across one "<<endl;opt.numsteps=1;}
    }
    if (opt.outputformat<OUTASCII && opt.outputformat>OUTHDF){
        cerr<<"Output format not valid, defaulting to ascii"<<endl; opt.outputformat=OUTASCII;
    }
    if (opt.outdataformat<0){
        cerr<<"Output data requested not valid, defaulting to minimal output"<<endl; opt.outdataformat=0;
    }

    //now set description
    opt.description=(char*)"Produce tree in direction of  ";
    if(opt.isearchdirection==SEARCHPROGEN)      opt.description+=(char*)"progenitors |";
    else if(opt.isearchdirection==SEARCHDESCEN)      opt.description+=(char*)"descendants |";
    else if(opt.isearchdirection==SEARCHALL)      opt.description+=(char*)"both forward (descendant) and backward (progenitor) |";

    opt.description=(char*)"TreeFrog Tree constructed by identifying the link with the highest value of ";
    if(opt.imerittype==NsharedN1N2)      opt.description+=(char*)"Nshared^2/Nh/Np |";
    else if(opt.imerittype==NsharedN1)   opt.description+=(char*)"Nshared/Nh | ";
    else if(opt.imerittype==Nshared)     opt.description+=(char*)"Nshared |";
    else if (opt.imerittype==Nsharedcombo) opt.description=(char*)"Nshared/Nh+(Nshared^2/Nh/Np) so as to weight progenitors that contribute similar amounts by how much of their mass contributes to the new object | ";

    opt.description=(char*)"Optimal temporal merits are set by ";
    if(opt.iopttemporalmerittype==GENERALIZEDMERITTIME)  opt.description+=(char*)"a generalized temporal merit taking into account average time evolution |";
    else if(opt.iopttemporalmerittype==GENERALIZEDMERITTIMEPROGEN)  opt.description+=(char*)"a generalized temporal merit taking into account average time evolution and maximise the ranking of the progenitor so that links always point to primary progen/descen |";

    opt.description=(char*)"Multistep criterion to keep searching ";
    if(opt.imultsteplinkcrit==MSLCMISSING)  opt.description+=(char*)" only if missing a link |";
    else if(opt.imultsteplinkcrit==MSLCMERIT)  opt.description+=(char*)" if missing a link or only link low merit |";
    else if(opt.imultsteplinkcrit==MSLCPRIMARYPROGEN)  opt.description+=(char*)" if missing a link or if link is secondary progenitor when constructing descendant tree |";
    else if(opt.imultsteplinkcrit==MSLCMERITPRIMARYPROGEN)  opt.description+=(char*)" if missing a link, low merit or if link is secondary progenitor when constructing descendant tree |";

    opt.description+=(char*)"Tree built using ";
    opt.description+=static_cast<ostringstream*>( &(ostringstream() << opt.numsteps) )->str();
    opt.description+=(char*)" temporal steps | ";

    opt.description+=(char*)"Particle types for matching limited to ";
    if (opt.itypematch==ALLTYPEMATCH) opt.description+=(char*)" all |";
    else {opt.description+=(char*)" part type ";opt.description+=static_cast<ostringstream*>( &(ostringstream() << opt.itypematch) )->str();}
    opt.description+=(char*)" | ";

    if (opt.particle_frac<1 && opt.particle_frac>0) {
        opt.description+=(char*)" Fractions of paritcles from which merit calculated with ";opt.description+=static_cast<ostringstream*>( &(ostringstream() << opt.particle_frac) )->str();
        opt.description+=(char*)" with lower particle number limit of ";opt.description+=static_cast<ostringstream*>( &(ostringstream() << opt.min_numpart) )->str();
        opt.description+=(char*)" | ";
    }

    opt.description+=(char*)"Merit threshold is  ";
    opt.description+=static_cast<ostringstream*>( &(ostringstream() << opt.meritlimit) )->str();
    opt.description+=(char*)" | ";

    cout<<"TreeFrog running with "<<endl<<"---------"<<endl<<opt.description<<endl<<"---------"<<endl<<endl;
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
    \arg \b \e -B < for VELOCIraptor catalogs can be \ref INASCII (ascii) \ref INBINARY (binary) or \ref INHDF (hdf5) or \ref INADIOS (adios) >
    \arg \b \e -F < 0/1 field objects in VELOCIraptor output in separate file >
    \arg \b \e -o < output filename >
    \arg \b \e -O < output type, whether ASCII, HDF5, etc >
    \arg \b \e -d < output data format, what is included. Minimal data is 0, if merits desired, 1>  >
    \arg \b \e -c < produce cross catalog match (0 halo tree ,1 cross catalog where code acts like particle correlator, matching analogues between simulations >

    \section linkingconfig Linking options
    \arg \b \e -t < number of steps over to find links >
    \arg \b \e -n < maximum id of particles, used to allocate an array of this size so that ids are mapped to an index, allows quick check of what groups a particle belongs to>
    \arg \b \e -C < cross correlation function type to identify main progenitor (see \ref Options.imerittype, which can be  \ref NsharedN1N2 standard merit function \f[ \mathcal{M}=N_{\rm sh}^2/N_1/N_2 \f], \ref Nshared >
    \arg \b \e -T < type of particles to cross correlate (\ref opt.itypematch, which can be \ref ALLTYPEMATCH is all particle types, \ref DMTYPEMATCH is DM particle types
    \ref GASTYPEMATCH is GAS particle types \ref STARTYPEMATCH is STAR particle types, \ref DMGASTYPEMATCH is both DM and GAS particle types >
    \arg \b \e -S < significance of cross match relative to Poisson noise >
    \arg \b \e -a < merit limit below which will search for links at earlier snapshot (as if no match is found)>
    \arg \b \e -f < fraction of particles to use to calculate merit, assumes that particles are sorted in some meaningful order in the input files. For example, VELOCIraptor sorts particles by binding energy
    so particle 0 is more bound than the last particle associated with the object. These particles would be represented of the dynamically coldest most dense core of the object.
    \arg \b \e -p <minimum number of particles used when fractions of particles are invoked>;

    \section outputconfig Output options
    \arg \b \e -H < offset snapshot values by this number >
    \arg \b \e -h < adjust Halo ID values stored in group catalog, useful for matching these values to those stored in .properties files produced by the halo finder. output is ID+(snap+snapshotvaloffset)*haloIDval >
    \arg \b \e -d < adjust Halo ID by this offset value >

    \section otherconfig Other options
    \arg \b \e -D < adjust particle IDs for nIFTY cross catalogs across simulations >
    \arg \b \e -M < Mapping of particle ids to index  (\ref opt.imapping with \ref DNOMAP no mapping of ids to indices, \ref DMEMEFFICIENTMAP produces a unique ID to index mapping, needs more memory while building the map>
    \arg \b \e -v < verbose flag 1/0 >

    \section mpiconfig MPI options
    \arg \b \e -m < number of items per mpi thead, use for load balacing. Use 0 if no mpi used when building halo catalog >
    \arg \b \e -z < number of mpi theads used to calculate load balacing. If >0 this used with one actual mpi thread but determines load balancing based on desired number of mpi threads. Write load balancing file and terminates. If 0 (default) normal operation >
    \arg \b \e -Z < whether to write output in parallel (0/1) >

*/
