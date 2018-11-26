/*! \file ui.cxx
 *  \brief user interface
 */

#include "stf.h"

///routine to get arguments from command line
/// \todo alter interface as now need to be able to specify only smdata file (no grid, res, normalized res) and functionality to specify eps, background fof search, etc
void GetArgs(int argc, char *argv[], Options &opt)
{
    int option;
    int NumArgs = 0;
    int configflag=0;
    while ((option = getopt(argc, argv, ":C:I:i:s:Z:o:G:S:B:t:")) != EOF)
    {
        switch(option)
        {
            case 'C':
                opt.pname = optarg;
                configflag=1;
                NumArgs += 2;
                break;
            case 'I':
                opt.inputtype = atoi(optarg);
                NumArgs += 2;
                break;
            case 'i':
                opt.fname = optarg;
                NumArgs += 2;
                break;
            case 's':
                opt.num_files = atoi(optarg);
                NumArgs += 2;
                break;
            case 'Z':
                opt.nsnapread = atoi(optarg);
                NumArgs += 2;
                break;
            case 'o':
                opt.outname = optarg;
                NumArgs += 2;
                break;
            case 'G':
                opt.gnsphblocks = atoi(optarg);
                NumArgs += 2;
                break;
            case 'S':
                opt.gnstarblocks = atoi(optarg);
                NumArgs += 2;
                break;
            case 'B':
                opt.gnbhblocks = atoi(optarg);
                NumArgs += 2;
                break;
            case 't':
                opt.ramsessnapname = optarg;
                NumArgs += 2;
                break;
            case '?':
                usage();
        }
    }
    if(configflag){
        cout<<"Reading config file"<<endl;
        GetParamFile(opt);
    }
    else {
        cout<<"NO CONFIG FILE PASSED! Using default values"<<endl;
    }
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
    cerr<<"-C <configuration file (overrides other options)> "<<endl;
    cerr<<"-I <input format [Gadget (Default) "<<IOGADGET<<", HDF (if implemented)"<<IOHDF<<", TIPSY "<<IOTIPSY<<", RAMSES "<<IORAMSES<<", HDF "<<IOHDF<<", NCHILADA "<<IONCHILADA<<">"<<endl;
    cerr<<"-i <input file> "<<endl;
    cerr<<"-s <number of files per output for gadget input 1 [default]>"<<endl;
    cerr<<"-Z <number of threads used in parallel read ("<<opt.nsnapread<<")>"<<endl;
    cerr<<"-o <output filename>"<<endl;
    cerr<<" ===== EXTRA OPTIONS FOR GADGET INPUT ====== "<<endl;
    cerr<<"-g <number of extra sph/gas blocks for gadget>"<<endl;
    cerr<<"-s <number of extra star blocks for gadget>"<<endl;
    cerr<<"-b <number of extra bh blocks for gadget>"<<endl;
    cerr<<" ===== EXTRA OPTIONS REQUIRED FOR RAMSES INPUT ====== "<<endl;
    cerr<<"-t <ramses snapnumber>"<<endl;

#ifdef USEMPI
    }
    MPI_Finalize();
#endif
    exit(1);
}

/*!
    \page configopt  Configuration options

    An example configuration file can be found the examples directory (see for instance <a href="../../examples/sample.cfg">sample</a>). within the repository.
    Only the keywords listed here will be used, all other words/characters are ignored.
    Note that if misspell a keyword it will not be used. One can check the options used by examining <b><em>foo</em>.configuration</b>,
    where <em>foo</em> is your base output filename \ref Options.outname. This file lists all the options. Since this file is
    always written <b><em>DO NOT</em></b> name your input configuration file <b><em>foo</em>.configuration</b>.

    There are numerous key words that can be passed to change the values stored in \ref Options. Here we list them
    \section outconfig Output related.
    See \ref outputs for types of output (and \ref io.cxx for implementation of the code)

    \arg <b> \e Output </b> Output base name. Overrides the name passed with the command line argument <b> \e -o </b>. Only implemented for completeness. \ref Options.outname \n
    \arg <b> \e Write_group_array_file </b> When producing output also produce a file which lists for every particle the group they belong to. Can be used with \b tipsy format or to tag every particle. \ref Options.iwritefof
    \arg <b> \e Output_den </b> A filename for storing the intermediate step of calculating local densities. This is particularly useful if the code is not compiled with \b STRUCDEN & \b HALOONLYDEN (see \ref STF-makeflags). \ref Options.smname \n

    \section searchconfig Parameters related to search type.
    See \ref io.cxx (and related ios like \ref gadgetio.cxx), \ref search.cxx, \ref fofalgo.h for extra details
    \arg <b> \e Particle_search_type </b> A integer describing what types of particles are searched \ref Options.partsearchtype. (See \ref SEARCHTYPES) > \n
        - \b 1 \e all particles are searched
        - \b 2 \e DarkMatter particles (which are typically defined as type 1,2,3 for gadget) are searched
        - \b 3 \e Star particles (which are typically defined as type 4 for gadget) are searched
        - \b 4 \e Gas particles (which are typically defined as type 0 for gadget) are searched
    \n
    \arg <b> \e Baryon_searchflag </b> 0/1/2 indicating gas/stellar search done separately from DM search \ref Options.iBaryonSearch. \n
        - \b 1 field search run as normal and then substructure search for baryons run using baryons identified in field search.
        - \b 2 field search also altered to treat baryons differently, allowing only DM particles to be used as head links (ie link dm-dm, dm-baryon, but not baryon-baryon nor baryon-dm).
        - \b 0 is do nothing special for baryon particles.
    \n
    \arg <b> \e Singlehalo_search </b> 0/1 flag indicates that no field search is going to be run and the entire volume will be treated as a background region. Useful if searching for substructures in non-cosmological simulations. But can also be co-opted for other searches using different outlier criteria and FOF algorithms. \ref Options.iSingleHalo \n

    \section localdensityconfig Parameters related to local density estimator.
    See \ref localfield.cxx, \ref bgfield.cxx & \ref localbgcomp.cxx for more details

    \arg <b> \e Nsearch_velocity </b> number of velocity neighbours used to calculate velocity density, adjust \ref Options.Nvel (suggested value is 32) \n
    \arg <b> \e Nsearch_physical </b> number of physical neighbours searched for Nv to calculate velocity density  \ref Options.Nsearch (suggested value is 256) \n
    \arg <b> \e Cell_fraction </b> fraction of a halo contained in a subvolume used to characterize the background  \ref Options.Ncellfac \n
    \arg <b> \e Grid_type </b> integer describing type of grid used to decompose volume for substructure search  \ref Options.gridtype (see \ref GRIDTYPES) \n
        - \b 1 \e standard physical shannon entropy, balanced KD tree volume decomposition into cells
        - \b 2 \e phase phase-space shannon entropy, balanced KD tree volume decomposition into cells
        - \b 3 \e simple simple physical balanced KD tree decomposition of volume into cells

    \section fofconfig Parameters related to FOF search
    See \ref search.cxx \ref fofalgo.h for more details

    \subsection fofsubconfig Configuration for substructure search
    \arg <b> \e Search_for_substructure </b> By default field objects are searched for internal substructures but can disable this by setting this to 0 \n
    \arg <b> \e Keep_FOF </b> if field 6DFOF search is done, allows to keep structures found in 3DFOF (can be interpreted as the inter halo stellar mass when only stellar search is used).\n
    \arg <b> \e FoF_search_type </b> There are several substructure FOF criteria implemented (see \ref FOFTYPES for more types and \ref fofalgo.h for implementation) \n
        - \b 1 \e standard phase-space based, well tested VELOCIraptor criterion.

    \arg <b> \e Outlier_threshold </b> threshold of sigma level of outliers to be searched (default is 2.5) \ref Options.ellthreshold \n
    \arg <b> \e Significance_level </b> minimum significance level of group (default is 1) \ref Options.siglevel \n
    \arg <b> \e Velocity_ratio </b> speed ratio used in linking particles \ref Options.Vratio \n
    \arg <b> \e Velocity_opening_angle </b> angle between velocities when linking (in units of \f$ \pi \f$) \ref Options.thetaopen \n
    \arg <b> \e Physical_linking_length </b> physical linking length used in fof (if cosmological gadget file then assumed to be in units of inter particle spacing \ref gadgetio.cxx, if loading in a single halo then based on average interparticle spacing calculated in \ref haloproperties.cxx) \ref Options.ellphys \n
    \arg <b> \e Minimum_size </b> Minimum group (substructure) size \ref Options.MinSize \n

    \arg <b> \e Iterative_searchflag </b> 1/0 use interactive substructure search which is designed to first identify spatially compact candidate outlier regions and then relaxes the criteria to find the more diffuse (in phase-space) regions associate with these candidate structures \ref Options.iiterflag \n
    \arg <b> \e Iterative_threshold_factor </b> factor multiplied with \ref Options.ellthreshold when using iterative method and identifying outlier regions associated with the initial candidate list of spatially compact outlier groups. Typical values are \f$ \sim 1 \f$ \ref Options.ellfac \n
    \arg <b> \e Iterative_linking_length_factor </b> factor multiplied with \ref Options.ellphys when using iterative method and identifying outlier regions associated with the initial candidate list of spatially compact outlier groups. Typical values are \f$ \sim 2 \f$ \ref Options.ellxfac \n
    \arg <b> \e Iterative_Vratio_length_factor </b> factor multiplied with \ref Options.Vratio when using iterative method and identifying outlier regions associated with the initial candidate list of spatially compact outlier groups. Typical values are \f$ \sim 1 \f$ \ref Options.vfac \n
    \arg <b> \e Iterative_ThetaOp_length_factor </b> factor multiplied with \ref Options.thetaopen when using iterative method and identifying outlier regions associated with the initial candidate list of spatially compact outlier groups. Typical values are \f$ \sim 1 \f$ \ref Options.thetafac \n
    \arg <b> \e CMrefadjustsubsearch_flag </b> 1/0 flag indicating whether particles are moved to the rough CM velocity frame of the background before substructures are searched for. \ref Options.icmrefadjust \n

    \subsection foffieldconfig Configuration for field search
    \arg <b> \e FoF_Field_search_type </b> There are several FOF criteria implemented to search for so-called field objects (see \ref FOFTYPES for more types and \ref fofalgo.h for implementation) \ref Options.fofbgtype \n
        - \b 5 \e standard 3D FOF based algorithm
        - \b 4 \e standard 3D FOF based algorithm <b> FOLLOWED </b> by 6D FOF search using the velocity scale defined by the largest halo on <b> ONLY </b> particles in 3DFOF groups
        - \b 3 \e standard 3D FOF based algorithm <b> FOLLOWED </b> by 6D FOF search using the velocity scale for each 3DFOF group
    \arg <b> \e Minimum_halo_size </b> Allows field objects (or so-called halos) to require a different minimum size (typically would be <= \ref Options.MinSize. Default is -1 which sets it to \ref Options.MinSize) \ref Options.HaloMinSize \n
    \arg <b> \e Halo_linking_length_factor </b> allows one to use different physical linking lengths between field objects and substructures.  (Typically for 3DFOF searches of dark matter haloes, set to value such that this times \ref Options.ellphys = 0.2 the interparticle spacing when examining cosmological simulations ) \ref Options.ellhalophysfac \n
    \arg <b> \e Halo_velocity_linking_length_factor </b> allows one to use different velocity linking lengths between field objects and substructures when using 6D FOF searches.  (Since in such cases the general idea is to use the local velocity dispersion to define a scale, \f$ \geq5 \f$ times this value seems to correctly scale searches) \ref Options.ellhalovelfac \n
    \arg <b> \e Halo_6D_linking_length_factor </b> allows one to use different linking lengths between 3DFOF and 6DFOF field search. Typically values are \f$ \sim 1 \f$ \n
    \arg <b> \e Halo_6D_vel_linking_length_factor </b> scaling applied to dispersions used in 6DFOF field search. Typical values are \f$ \geq 1.25 \f$ \n

    \arg <b> \e Halo_core_search </b> 0/1/2 flag allows one to explicitly search for large 6D FOF cores that are indicative of a recent major merger. Since substructure is defined on the scale of the maximum cell size and major mergers typically result two or more phase-space dense regions that are \e larger than the cell size used in reasonable substructure searches, one can identify them using this search. The overall goal is to treat these objects differently than a substructure. However, if 2 is set, then smaller core is treated as substruture and all particles within the FOF envelop are assigned to the cores based on their phase-space distance to core particles \ref Options.iHaloCoreSearch \n
    \arg <b> \e Use_adaptive_core_search </b> 0/1 flag allows one to run complex adaptive phase-space search for large 6D FOF cores and then use these linking lengths to separate mergers. 0 is simple high density dispersively cold cores with v scale adaptive, 1 is adaptive with simple dispersions in both x and v.
    \arg <b> \e Use_phase_tensor_core_growth </b> 0/1 flag allows one to run complex phase-space growth of merger remnants (6D FOF cores found). 0 is assignment with simple x and v dispersion to nearest core particle, 1 is phase-space tensor distance assignemnt to CM of core.
    \arg <b> \e Halo_core_ellx_fac </b> scaling applied to linking length when identifying merger remnants. Typically values are \f$ \sim0.5 \f$  \ref Options.halocorexfac
    \arg <b> \e Halo_core_ellv_fac </b> scaling applied to local dispersion to define the velocity scale used to identify merger remnants. Typically values are \f$ \sim1 \f$  \ref Options.halocorevfac
    \arg <b> \e Halo_core_ncellfac </b> used to determine the minimum number of particles a merger remnants is composed of, specifically \f$ N_{\rm min}= f_{\rm ncell}* N_{\rm S} \f$. Typically values are 0.005
    \arg <b> \e Halo_core_adaptive_sigma_fac </b> used when running fully adaptive core search, specifies the width of the physical linking length in configuration space dispersion (think of this as how many sigma to include). Typically values are \f$ \sim2 \f$. This has been tested on hydrodynamnical simulations to separate galaxy mergers. \ref Options.halocoresigmafac
    \arg <b> \e Halo_core_num_loops </b> allows the core search to iterate, shrinking the velocity linking length to used till the number of cores identified decreases or this limit is reached. Allows apative search with larger linking length to be robust.  Typically values are \f$ \sim5 \f$ with loops only running often twice. \ref Options.halocorenumloops.
    \arg <b> \e Halo_core_loop_ellx_fac </b> Factor by which configuration linking length is decreased when running loops for core search.  Typically values are \f$ \sim0.75 \f$. \ref Options.halocorexfaciter
    \arg <b> \e Halo_core_loop_ellv_fac </b> Factor by which velocity linking length is decreased when running loops for core search.  Typically values are \f$ \sim0.75 \f$. \ref Options.halocorevfaciter
    \arg <b> \e Halo_core_loop_elln_fac </b> Factor by which min group size is changed when running loops for core search.  Typically values are \f$ \sim1.25 \f$. \ref Options.halocorenumfaciter
    \arg <b> \e Halo_core_phase_significance </b> Significance a core must be in terms of phase-space distance scaled by dispersions (sigma).

    \section unbindconfig Unbinding Parameters
    See \ref unbinding and \ref unbind.cxx for more details

    \arg <b> \e Unbind_flag  </b> whether or not substructures should be passed through an unbinding routine. \ref Options.uinfo & \ref UnbindInfo.unbindflag \n
    \arg <b> \e Allowed_kinetic_potential_ratio  </b> ratio of kinetic to potential energy at which a particle is still considered bound, ie: particle is still bound
    if \f$ \alpha T+W<0 \f$, so \f$ \alpha=1 \f$ would be standard unbinding and \f$ \alpha<1 \f$ allows one to identify unbound tidal debris.
    Given that <b> VELOCIraptor </b> was designed to identify tidal streams, it makes little sense to have this set to 1 unless explicitly required.
    Note that the code still separates particles into bound and unbound. Typical values of \f$ \alpha\geq 0.2 \f$ seems to minimize the number of false positives
    in tidal debris while still identifying completely unbound tidal debris. \ref Options.uinfo & \ref UnbindInfo.Eratio \n
    \arg <b> \e Min_bound_mass_frac </b> Designed to demand a substructure still have a minimum amount of self-bound mass. \ref Options.uinfo & \ref UnbindInfo.minEfrac \n
    \arg <b> \e Bound_halos </b> 0/1/2 flag to make sure field objects such as haloes are self bound before (use 1) and also after (use 2) substructures have been identified
    and extracted from the halo. Demanding boundness after substructure search can have interesting consequences as it is possible that a multiple merger will appear as
    a single FOF halo, however all with all the cores removed, the FOF halo is actually an unbound structure. \ref Options.iBoundHalos \n
    \arg <b> \e Keep_background_potential </b> 1/0 flag When determining whether a structure is self-bound, the approach taken is to treat the candidate structure in isolation. Then determine the velocity reference frame to determine the kinetic energy of each particle and remove them. However, it is possible one wishes to keep the background particles when determining the potential, that is once one starts unbinding, don't treat the candidate structure in isolation but in a background sea. When finding tidal debris, it is useful to keep the background. \ref Options.uinfo & \ref UnbindInfo.bgpot \n
    \arg <b> \e Kinetic_reference_frame_type </b> specify kinetic frame when determining whether particle is bound.
    Default is to use the centre-of-mass velocity frame (0) but can also use region around minimum of the potential (1). \ref Options.uinfo & \ref UnbindInfo.cmvelreftype \n
    \arg <b> \e Min_npot_ref </b> Set the minimum number of particles used to calculate the velocity of the minimum of the potential (10). \ref Options.uinfo & \ref UnbindInfo.Npotref \n
    \arg <b> \e Frac_pot_ref </b> Set the fraction of particles used to calculate the velocity of the minimum of the potential (0.1). \ref Options.uinfo & \ref UnbindInfo.fracpotref \n
    \arg <b> \e Unbinding_type </b> Set the unbinding criteria, either just remove particles deemeed "unbound", that is those with \f$ \alpha T+W>0\f$, choosing \ref UPART. Or with \ref USYSANDPART
    removes "unbound" particles till system also has a true bound fraction > \ref UnbindInfo.minEfrac.
    \arg <b> \e Softening_length </b> Set the (simple plummer) gravitational softening length. \ref UnbindInfo.eps 

    \section cosmoconfig Units & Cosmology
    \subsection unitconfig Units
    \arg <b> \e Length_unit </b> change the input length unit \ref Options.L \n
    \arg <b> \e Velocity_unit </b> change the input velocity unit \ref Options.V \n
    \arg <b> \e Mass_unit </b> change the input mass unit \ref Options.M \n
    \arg <b> \e Gravity </b> change the gravity unit. should be set such that \f$ v^2=GM/r \f$.   \ref Options.G \n
    \arg <b> \e Hubble_unit </b> change the value of Hubble expansion (from normal 100 km/s/Mpc to units used, that is velocity unit/length unit) \ref Options.G \n
    \arg <b> \e Mass_value </b> if not mass is stored using the option \b NOMASS (see \ref STF-makeflags) then this is the mass of the particles \ref Options.MassValue \n
    \arg <b> \e Length_unit_to_kpc </b> Specify the conversion factor from the output unit to kpc \ref Options.lengthtokpc \n
    \arg <b> \e Velocity_unit_to_kms </b> Specify the conversion factor from the output unit to kms \ref Options.velocitytokms \n
    \arg <b> \e Mass_unit_to_solarmass </b> Specify the conversion factor from the output unit to solar masses \ref Options.masstosolarmass \n

    \subsection cosmologyconfig Cosmology
    \arg <b> \e Period </b> if not defined in the input data one can pass the period in the input units of the data. This is not always necessary, for instance the gadget snapshot format has this in the header.  \ref Options.p \n
    \arg <b> \e Scale_factor </b> scale factor time if not defined in the input data. This is not always necessary, for instance the gadget snapshot format has this in the header.  \ref Options.a \n
    \arg <b> \e h_val </b> the "little h" value often used in cosmological simulation when defining expansion rates, distance scales, etc. This is not always necessary, for instance the gadget snapshot format has this in the header.  \ref Options.h \n
    \arg <b> \e Omega_m </b> matter density in units of the critical density used in cosmological simulations. This is not always necessary, for instance the gadget snapshot format has this in the header.  \ref Options.Omega_m \n
    \arg <b> \e Omega_Lambda </b> energy density of the cosmological constant (or can technically be used for a dark energy fluid) in units of the critical density used in cosmological simulations. This is not always necessary, for instance the gadget snapshot format has this in the header.  \ref Options.Omega_Lambda \n
    \arg <b> \e Omega_cdm </b> dark matter density in units of the critical density used in cosmological simulations. This is not always necessary since for simple pure DM cosmological simulations this is easily determined by the total mass in the simulation. <em> However, for multiple resolution or non-standard DM models, this should be provided. </em> \ref Options.Omega_cdm \n
    \arg <b> \e Omega_b </b> baryon density in units of the critical density used in cosmological simulations. This is not always necessary since for simple pure adiabatic cosmological simulations this is easily determined by the total gas mass in the simulation.  \ref Options.Omega_b \n
    \arg <b> \e w_of_DE </b> equation of state of the dark energy fluid, \f$ w=\frac{p}{\rho} \f$. This is not necessary unless one is using a cosmological simulation with a \f$ w\neq -1 \f$ \ref Options.w_de \n
    \arg <b> \e Virial_density </b> virial overdensity in units of the background matter density used in cosmological simulations. This is not absolutely necessary for gadget format as for this input, the Bryan & Norman 1998 virial density is calculated based on a LCDM cosmology as \ref Options.virlevel is set to -1 and if left, a virial density is calculated. If set and a gadget input is used, this overrides the Bryan & Norman calculation. \n
    \arg <b> \e Critical_density </b> critical density in input units used in cosmological simulations. This is not absolutely necessary for gadget format as it can be calculated from the quantities in the header for typical GR cosmologies. \ref Options.rhobg \n

    \section otherconfigs Other configuration options
    \arg <b> \e Effective_Resolution </b> If running a multiple resolution cosmological zoom simulation, simple method of scaling the linking length by using the period, ie: \f$ p/N_{\rm eff} \f$ \ref Options.Neff \n
    \arg <b> \e Snapshot_value </b> If halo ids need to be offset to some starting value based on the snapshot of the output, which is useful for some halo merger tree codes, one can specific a snapshot number, and all halo ids will be listed as internal haloid + \f$ sn\times10^{12}\f$. \ref Options.snapshotvalue \n
    \arg <b> \e Verbose </b> 2/1/0 flag indicating how talkative the code is (2 very verbose, 1 verbose, 0 quiet). \ref Options.iverbose \n
    \arg <b> \e Inclusive_halo_mass </b> 1/0 flag indicating whether inclusive masses are calculated for field objects. \ref Options.iInclusiveHalo \n

    \section ioconfigs I/O options
    \arg <b> \e Cosmological_input </b> 1/0 indicating that input simulation is cosmological or not. With cosmological input, a variety of length/velocity scales are set to determine such things as the virial overdensity, linking length. \ref Options.icosmologicalin \n
    \arg <b> \e Input_chunk_size </b> Amount of information to read from input file in one go (100000). \ref Options.inputbufsize \n
    \arg <b> \e Write_group_array_file </b> 0/1 flag indicating whether write a single large tipsy style group assignment file is written. \ref Options.iwritefof \n
    \arg <b> \e Separate_output_files </b> 1/0 flag indicating whether separate files are written for field and subhalo groups. \ref Options.iseparatefiles \n
    \arg <b> \e Binary_output </b> 3/2/1/0 flag indicating whether output is hdf, binary or ascii. \ref Options.ibinaryout, \ref OUTADIOS, \ref OUTHDF, \ref OUTBINARY, \ref OUTASCII \n
    \arg <b> \e Extensive_halo_properties_output </b> 1/0 flag indicating whether to calculate/output even more halo properties. \ref Options.iextrahalooutput \n
    \arg <b> \e Extended_output </b> 1/0 flag indicating whether produce extended output for quick particle extraction from input catalog of particles in structures \ref Options.iextendedoutput \n
    \arg <b> \e Comoving_units </b> 1/0 flag indicating whether the properties output is in physical or comoving little h units. \ref Options.icomoveunit \n

    \section inputflags input flags related to varies input formats
    \arg <b> \e NSPH_extra_blocks </b> If gadget snapshot is loaded one can specific the number of extra <b> SPH </b> blocks are read/in the file. \ref Options.gnsphblocks \n
    \arg <b> \e NStar_extra_blocks </b> If gadget snapshot is loaded one can specific the number of extra <b> Star </b> blocks are read/in the file. \ref Options.gnstarblocks \n
    \arg <b> \e NBH_extra_blocks </b> If gadget snapshot is loaded one can specific the number of extra <b> Black hole </b> blocks are read/in the file. \ref Options.gnbhblocks \n

    \arg <b> \e HDF_name_convention </b> HDF dataset naming convection. See \ref hdfitems.h for what naming conventions are available and what names exist. Currently have \ref HDFNUMNAMETYPES. \ref Options.ihdfnameconvention \n
    \arg <b> \e Input_includes_star_particle </b> If star particle specific information is in the input file. \ref Options.iusestarparticles \n
    \arg <b> \e Input_includes_bh_particle </b> If bh/sink particle specific information is in the input file. \ref Options.iusesinkparticles \n
    \arg <b> \e Input_includes_wind_particle </b> If wind particle specific information is in the input file. \ref Options.iusewindparticles \n
    \arg <b> \e Input_includes_tracer_particle </b> If tracer particle specific information is in the input file. \ref Options.iusetracerparticles \n
    \arg <b> \e Input_includes_star_particle </b> If star particle specific information is in the input file. \ref Options.iusestarparticles \n

    \section mpiconfigs MPI specific options
    \arg <b> \e MPI_part_allocation_fac </b> Memory allocated in mpi mode to store particles is (1+factor)* the memory need for the initial mpi decomposition.
    This factor should be >0 and is mean to allow a little room for particles to be exchanged between mpi threads withouth having to require new memory allocations and copying
    of data. \ref Options.mpipartfac \n
    \arg <b> \e MPI_particle_total_buf_size </b> Total memory size in bytes used to store particles in temporary buffer such that
    particles are sent to non-reading mpi processes in one communication round in chunks of size buffer_size/NProcs/sizeof(Particle). \ref Options.mpiparticlebufsize \n



    */

///Read parameters from a parameter file. For list of currently implemented options see \ref configopt
///\todo still more parameters that can be adjusted
void GetParamFile(Options &opt)
{
    string line,sep="=";
    string tag,val;
    char buff[1024],*pbuff,tbuff[1024],vbuff[1024],fname[1024];
    fstream paramfile,cfgfile;
    if (!FileExists(opt.pname)){
            cerr<<"Config file: "<<opt.pname <<" does not exist or can't be read, terminating"<<endl;
#ifdef USEMPI
            MPI_Abort(MPI_COMM_WORLD,9);
#else
            exit(9);
#endif
    }
    paramfile.open(opt.pname, ios::in);
    unsigned j,k;
    //first find output name, determine the number of valid entries in
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
                    if (strcmp(tbuff, "Output")==0){
                        opt.outname=new char[1024];
                        strcpy(opt.outname,vbuff);
                        break;
                    }
                }
            }
        }
#ifndef SWIFTINTERFACE
        if (opt.outname==NULL) {
            cerr<<"No output name given, terminating"<<endl;
#ifdef USEMPI
            MPI_Finalize();
#endif
            exit(9);
        }
#endif
    }
    //sprintf(fname,"%s.cfg",opt.outname);
    //cfgfile.open(fname, ios::out);
    paramfile.clear();
    paramfile.seekg(0, ios::beg);
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
                    //cfgfile<<tbuff<<"="<<vbuff<<endl;
                    //store/read local density distribution funciton values
                    if (strcmp(tbuff, "Output_den")==0){
                        opt.smname=new char[1024];
                        sprintf(opt.smname,"%s.localden",opt.outname);
                    }
                    //config search type
                    else if (strcmp(tbuff, "Particle_search_type")==0)
                        opt.partsearchtype = atoi(vbuff);
                    else if (strcmp(tbuff, "FoF_search_type")==0)
                        opt.foftype = atoi(vbuff);
                    else if (strcmp(tbuff, "FoF_Field_search_type")==0)
                        opt.fofbgtype = atoi(vbuff);
                    else if (strcmp(tbuff, "Search_for_substructure")==0)
                        opt.iSubSearch = atoi(vbuff);
                    else if (strcmp(tbuff, "Keep_FOF")==0)
                        opt.iKeepFOF = atoi(vbuff);
                    else if (strcmp(tbuff, "Iterative_searchflag")==0)
                        opt.iiterflag = atoi(vbuff);
                    else if (strcmp(tbuff, "Unbind_flag")==0)
                        opt.uinfo.unbindflag = atoi(vbuff);
                    else if (strcmp(tbuff, "Baryon_searchflag")==0)
                        opt.iBaryonSearch = atoi(vbuff);
                    else if (strcmp(tbuff, "CMrefadjustsubsearch_flag")==0)
                        opt.icmrefadjust = atoi(vbuff);
                    else if (strcmp(tbuff, "Halo_core_search")==0)
                        opt.iHaloCoreSearch = atoi(vbuff);
                    else if (strcmp(tbuff, "Use_adaptive_core_search")==0)
                        opt.iAdaptiveCoreLinking = atof(vbuff);
                    else if (strcmp(tbuff, "Use_phase_tensor_core_growth")==0)
                        opt.iPhaseCoreGrowth = atof(vbuff);

                    //bg and fof parameters
                    else if (strcmp(tbuff, "Cell_fraction")==0)
                        opt.Ncellfac = atof(vbuff);
                    else if (strcmp(tbuff, "Grid_type")==0)
                        opt.gridtype = atoi(vbuff);
                    else if (strcmp(tbuff, "Nsearch_velocity")==0)
                        opt.Nvel = atoi(vbuff);
                    else if (strcmp(tbuff, "Nsearch_physical")==0)
                        opt.Nsearch = atoi(vbuff);
                    else if (strcmp(tbuff, "Outlier_threshold")==0)
                        opt.ellthreshold = atof(vbuff);
                    else if (strcmp(tbuff, "Significance_level")==0)
                        opt.siglevel = atof(vbuff);
                    else if (strcmp(tbuff, "Velocity_ratio")==0)
                        opt.Vratio = atof(vbuff);
                    else if (strcmp(tbuff, "Velocity_opening_angle")==0)
                        opt.thetaopen = atof(vbuff);
                    else if (strcmp(tbuff, "Physical_linking_length")==0)
                        opt.ellphys = atof(vbuff);
                    else if (strcmp(tbuff, "Velocity_linking_length")==0)
                        opt.ellvel = atof(vbuff);
                    else if (strcmp(tbuff, "Minimum_size")==0)
                        opt.MinSize = atoi(vbuff);
                    else if (strcmp(tbuff, "Minimum_halo_size")==0)
                        opt.HaloMinSize = atoi(vbuff);
                    else if (strcmp(tbuff, "Halo_linking_length_factor")==0)
                        opt.ellhalophysfac = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_velocity_linking_length_factor")==0)
                        opt.ellhalovelfac = atof(vbuff);
                    //specific to 6DFOF field search
                    else if (strcmp(tbuff, "Halo_6D_linking_length_factor")==0)
                        opt.ellhalo6dxfac = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_6D_vel_linking_length_factor")==0)
                        opt.ellhalo6dvfac = atof(vbuff);
                    //specific search for 6d fof core searches
                    else if (strcmp(tbuff, "Halo_core_ellx_fac")==0)
                        opt.halocorexfac = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_core_ellv_fac")==0)
                        opt.halocorevfac = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_core_ncellfac")==0)
                        opt.halocorenfac = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_core_adaptive_sigma_fac")==0)
                        opt.halocoresigmafac = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_core_num_loops")==0)
                        opt.halocorenumloops = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_core_loop_ellx_fac")==0)
                        opt.halocorexfaciter = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_core_loop_ellv_fac")==0)
                        opt.halocorevfaciter = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_core_loop_elln_fac")==0)
                        opt.halocorenumfaciter = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_core_phase_significance")==0)
                        opt.halocorephasedistsig = atof(vbuff);

                    //for changing factors used in iterative search
                    else if (strcmp(tbuff, "Iterative_threshold_factor")==0)
                        opt.ellfac = atof(vbuff);
                    else if (strcmp(tbuff, "Iterative_linking_length_factor")==0)
                        opt.ellxfac = atof(vbuff);
                    else if (strcmp(tbuff, "Iterative_Vratio_factor")==0)
                        opt.vfac = atof(vbuff);
                    else if (strcmp(tbuff, "Iterative_ThetaOp_factor")==0)
                        opt.thetafac = atof(vbuff);
                    //for changing effective resolution when rescaling linking lengh
                    else if (strcmp(tbuff, "Effective_resolution")==0)
                        opt.Neff = atof(vbuff);
                    //for changing effective resolution when rescaling linking lengh
                    else if (strcmp(tbuff, "Singlehalo_search")==0)
                        opt.iSingleHalo = atoi(vbuff);

                    //units, cosmology
                    else if (strcmp(tbuff, "Length_unit")==0)
                        opt.L = atof(vbuff);
                    else if (strcmp(tbuff, "Velocity_unit")==0)
                        opt.V = atof(vbuff);
                    else if (strcmp(tbuff, "Mass_unit")==0)
                        opt.M = atof(vbuff);
                    else if (strcmp(tbuff, "Hubble_unit")==0)
                        opt.H = atof(vbuff);
                    else if (strcmp(tbuff, "Gravity")==0)
                        opt.G = atof(vbuff);
                    else if (strcmp(tbuff, "Mass_value")==0)
                        opt.MassValue = atof(vbuff);
                    else if (strcmp(tbuff, "Period")==0)
                        opt.p = atof(vbuff);
                    else if (strcmp(tbuff, "Scale_factor")==0)
                        opt.a = atof(vbuff);
                    else if (strcmp(tbuff, "h_val")==0)
                        opt.h = atof(vbuff);
                    else if (strcmp(tbuff, "Omega_m")==0)
                        opt.Omega_m = atof(vbuff);
                    else if (strcmp(tbuff, "Omega_Lambda")==0)
                        opt.Omega_Lambda = atof(vbuff);
                    else if (strcmp(tbuff, "Critical_density")==0)
                        opt.rhobg = atof(vbuff);
                    else if (strcmp(tbuff, "Virial_density")==0)
                        opt.virlevel = atof(vbuff);
                    else if (strcmp(tbuff, "Omega_cdm")==0)
                        opt.Omega_cdm= atof(vbuff);
                    else if (strcmp(tbuff, "Omega_b")==0)
                        opt.Omega_b= atof(vbuff);
                    else if (strcmp(tbuff, "w_of_DE")==0)
                        opt.w_de= atof(vbuff);
                    //so units can be specified to convert to kpc, km/s, solar mass
                    //not necessarily the code units data is reported but the
                    //translation of these code units to these units
                    else if (strcmp(tbuff, "Length_unit_to_kpc")==0)
                        opt.lengthtokpc = atof(vbuff);
                    else if (strcmp(tbuff, "Velocity_to_kms")==0)
                        opt.velocitytokms = atof(vbuff);
                    else if (strcmp(tbuff, "Mass_to_solarmass")==0)
                        opt.masstosolarmass = atof(vbuff);
                    //unbinding
                    else if (strcmp(tbuff, "Softening_length")==0)
                        opt.uinfo.eps = atof(vbuff);
                    else if (strcmp(tbuff, "Allowed_kinetic_potential_ratio")==0)
                        opt.uinfo.Eratio = atof(vbuff);
                    else if (strcmp(tbuff, "Min_bound_mass_frac")==0)
                        opt.uinfo.minEfrac = atof(vbuff);
                    else if (strcmp(tbuff, "Bound_halos")==0)
                        opt.iBoundHalos = atoi(vbuff);
                    else if (strcmp(tbuff, "Keep_background_potential")==0)
                        opt.uinfo.bgpot = atoi(vbuff);
                    else if (strcmp(tbuff, "Kinetic_reference_frame_type")==0)
                        opt.uinfo.cmvelreftype = atoi(vbuff);
                    else if (strcmp(tbuff, "Min_npot_ref")==0)
                        opt.uinfo.Npotref = atoi(vbuff);
                    else if (strcmp(tbuff, "Frac_pot_ref")==0)
                        opt.uinfo.fracpotref = atof(vbuff);
                    else if (strcmp(tbuff, "Unbinding_type")==0)
                        opt.uinfo.unbindtype = atoi(vbuff);

                    //other options
                    else if (strcmp(tbuff, "Verbose")==0)
                        opt.iverbose = atoi(vbuff);
                    else if (strcmp(tbuff, "Write_group_array_file")==0)
                        opt.iwritefof = atoi(vbuff);
                    else if (strcmp(tbuff, "Snapshot_value")==0)
                        opt.snapshotvalue = HALOIDSNVAL*atoi(vbuff);
                    else if (strcmp(tbuff, "Inclusive_halo_masses")==0)
                        opt.iInclusiveHalo = atoi(vbuff);

                    //input related
                    else if (strcmp(tbuff, "Cosmological_input")==0)
                        opt.icosmologicalin = atoi(vbuff);
                    //input read related
                    else if (strcmp(tbuff, "Input_chunk_size")==0)
                        opt.inputbufsize = atol(vbuff);
                    else if (strcmp(tbuff, "MPI_particle_total_buf_size")==0)
                        opt.mpiparticletotbufsize = atol(vbuff);
                    //mpi memory related
                    else if (strcmp(tbuff, "MPI_part_allocation_fac")==0)
                        opt.mpipartfac = atof(vbuff);

                    //output related
                    else if (strcmp(tbuff, "Separate_output_files")==0)
                        opt.iseparatefiles = atoi(vbuff);
                    else if (strcmp(tbuff, "Binary_output")==0)
                        opt.ibinaryout = atoi(vbuff);
                    else if (strcmp(tbuff, "Comoving_units")==0)
                        opt.icomoveunit = atoi(vbuff);
                    else if (strcmp(tbuff, "Extensive_halo_properties_output")==0)
                        opt.iextrahalooutput = atoi(vbuff);
                    else if (strcmp(tbuff, "Extensive_gas_properties_output")==0)
                        opt.iextragasoutput = atoi(vbuff);
                    else if (strcmp(tbuff, "Extended_output")==0)
                        opt.iextendedoutput = atoi(vbuff);
                    else if (strcmp(tbuff, "Spherical_overdensity_halo_particle_list_output")==0)
                        opt.iSphericalOverdensityPartList = atoi(vbuff);

                    //gadget io related to extra info for sph, stars, bhs,
                    else if (strcmp(tbuff, "NSPH_extra_blocks")==0)
                        opt.gnsphblocks = atoi(vbuff);
                    else if (strcmp(tbuff, "NStar_extra_blocks")==0)
                        opt.gnstarblocks = atoi(vbuff);
                    else if (strcmp(tbuff, "NBH_extra_blocks")==0)
                        opt.gnbhblocks = atoi(vbuff);


                    //input related to info for stars, bhs, winds/tracers, etc
                    else if (strcmp(tbuff, "HDF_name_convention")==0)
                        opt.ihdfnameconvention = atoi(vbuff);
                    else if (strcmp(tbuff, "Input_includes_star_particle")==0)
                        opt.iusestarparticles = atoi(vbuff);
                    else if (strcmp(tbuff, "Input_includes_bh_particle")==0)
                        opt.iusesinkparticles = atoi(vbuff);
                    else if (strcmp(tbuff, "Input_includes_wind_particle")==0)
                        opt.iusewindparticles = atoi(vbuff);
                    else if (strcmp(tbuff, "Input_includes_tracer_particle")==0)
                        opt.iusetracerparticles = atoi(vbuff);
                    else if (strcmp(tbuff, "Input_includes_extradm_particle")==0)
                        opt.iuseextradarkparticles = atoi(vbuff);

                }
            }
        }
        paramfile.close();
        //cfgfile.close();
    }
}

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
    if (opt.iBaryonSearch && !(opt.partsearchtype==PSTALL || opt.partsearchtype==PSTDARK)) {
#ifdef USEMPI
    if (ThisTask==0)
#endif
        cerr<<"Conflict in config file: both gas/star/etc particle type search AND the separate baryonic (gas,star,etc) search flag are on. Check config\n";
#ifdef USEMPI
            MPI_Abort(MPI_COMM_WORLD,8);
#else
            exit(8);
#endif
    }
    if (opt.iBoundHalos && opt.iKeepFOF) {
#ifdef USEMPI
    if (ThisTask==0)
#endif
        cerr<<"Conflict in config file: Asking for Bound Field objects but also asking to keep the 3DFOF/then run 6DFOF. This is incompatible. Check config\n";
#ifdef USEMPI
            MPI_Abort(MPI_COMM_WORLD,8);
#else
            exit(8);
#endif
    }
    if (opt.HaloMinSize==-1) opt.HaloMinSize=opt.MinSize;

    if (opt.num_files<1){
#ifdef USEMPI
    if (ThisTask==0)
#endif
        cerr<<"Invalid number of input files (<1) \n";
#ifdef USEMPI
            MPI_Abort(MPI_COMM_WORLD,8);
#else
            exit(8);
#endif
    }

    if (opt.inputbufsize<1){
#ifdef USEMPI
    if (ThisTask==0)
#endif
        cerr<<"Invalid read buf size (<1)\n";
#ifdef USEMPI
            MPI_Abort(MPI_COMM_WORLD,8);
#else
            exit(8);
#endif
    }

    if (opt.lengthtokpc<=0){
#ifdef USEMPI
    if (ThisTask==0)
#endif
        cerr<<"Invalid unit conversion, length unit to kpc is <=0 or was not set. Update config file\n";
#ifdef USEMPI
            MPI_Abort(MPI_COMM_WORLD,8);
#else
            exit(8);
#endif
    }
    if (opt.velocitytokms<=0){
#ifdef USEMPI
    if (ThisTask==0)
#endif
        cerr<<"Invalid unit conversion, velocity unit to km/s is <=0 or was not set. Update config file\n";
#ifdef USEMPI
            MPI_Abort(MPI_COMM_WORLD,8);
#else
            exit(8);
#endif
    }
    if (opt.masstosolarmass<=0){
#ifdef USEMPI
    if (ThisTask==0)
#endif
        cerr<<"Invalid unit conversion, mass unit to solar mass is <=0 or was not set. Update config file\n";
#ifdef USEMPI
            MPI_Abort(MPI_COMM_WORLD,8);
#else
            exit(8);
#endif
    }

#ifdef USEMPI
    if (opt.mpiparticletotbufsize<(long int)(sizeof(Particle)*NProcs) && opt.mpiparticletotbufsize!=-1){
        if (ThisTask==0) {
            cerr<<"Invalid input particle buffer send size, mininmum input buffer size given paritcle byte size ";
            cerr<<sizeof(Particle)<<" and have "<<NProcs<<" mpi processes is "<<sizeof(Particle)*NProcs<<endl;
        }
        MPI_Abort(MPI_COMM_WORLD,8);
    }
    //if total buffer size is -1 then calculate individual buffer size based on default mpi size
    else if (opt.mpiparticletotbufsize==-1) {
        opt.mpiparticletotbufsize=MPIPartBufSize*NProcs*sizeof(Particle);
        opt.mpiparticlebufsize=MPIPartBufSize;
    }
    else {
        opt.mpiparticlebufsize=opt.mpiparticletotbufsize/NProcs/sizeof(Particle);
    }
    if (opt.mpipartfac<0){
        if (ThisTask==0) {
            cerr<<"Invalid MPI particle allocation factor, must be >0."<<endl;
        }
        MPI_Abort(MPI_COMM_WORLD,8);
    }
    else if (opt.mpipartfac>1){
        if (ThisTask==0) {
            cerr<<"WARNING: MPI Particle allocation factor is high (>1)."<<endl;
        }
    }
#endif

#ifndef USEHDF
    if (opt.ibinaryout==OUTHDF){
#ifdef USEMPI
    if (ThisTask==0)
#endif
        cerr<<"Code not compiled with HDF output enabled. Recompile with this enabled or change Binary_output.\n";
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,8);
#else
        exit(8);
#endif
}
#endif

#ifndef USEADIOS
    if (opt.ibinaryout==OUTADIOS){
#ifdef USEMPI
    if (ThisTask==0)
#endif
        cerr<<"Code not compiled with ADIOS output enabled. Recompile with this enabled or change Binary_output.\n";
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,8);
#else
        exit(8);
#endif
}
#endif

#ifdef USEMPI
    if (ThisTask==0) {
#endif
    cout<<"CONFIG INFO SUMMARY -------------------------- "<<endl;
    switch(opt.partsearchtype)
    {
            case PSTALL:
                cout<<"Searching all particles regardless of type"<<endl;
                break;
            case PSTDARK:
                cout<<"Searching dark matter particles "<<endl;
                break;
            case PSTGAS:
                cout<<"Searching gas particles "<<endl;
                break;
            case PSTSTAR:
                cout<<"Searching star particles"<<endl;
                break;
            case PSTBH:
                cout<<"Searching BH particles?! Really? Are there enough?"<<endl;
                break;
    }
    if (opt.fofbgtype==FOF6D) cout<<"Field objects found with 3d FOF"<<endl;
    else if (opt.fofbgtype==FOF6D) cout<<"Field objects found with initial 3d FOF then use single dispersion to find 6d FOFs"<<endl;
    else if (opt.fofbgtype==FOF6DADAPTIVE) cout<<"Field objects found with initial 3d FOF then use adaptive dispersion to find 6d FOFs"<<endl;
    else if (opt.fofbgtype==FOFSTNOSUBSET) cout<<"Field objects found with initial 3d FOF then use phase-space stream finding algorithm"<<endl;
    if (opt.fofbgtype<=FOF6D && opt.iKeepFOF) cout<<"Field objects found with initial 3d FOF then use dispersion to find 6d FOFs and 3dFOF objects kept as base objects (consider them inter halo stellar mass or ICL for stellar only searches) "<<endl;
    if (opt.iBaryonSearch==1) cout<<"Baryons treated separately for substructure search"<<endl;
    else if (opt.iBaryonSearch==2) cout<<"Baryons treated separately for substructure search and also FOF field search"<<endl;
    if (opt.iSingleHalo) cout<<"Field objects NOT searched for, assuming single Halo and subsearch using mean field first step"<<endl;
    cout<<"Allowed potential to kinetic ratio when unbinding particles "<<opt.uinfo.Eratio<<endl;
    if (opt.HaloMinSize!=opt.MinSize) cout<<"Field objects (aka Halos) have different minimum required size than substructures: "<<opt.HaloMinSize<<" vs "<<opt.MinSize<<endl;
    cout<<"Units: L="<<opt.L<<", M="<<opt.M<<", V="<<opt.V<<", G="<<opt.G<<endl;
    if (opt.ibinaryout) cout<<"Binary output"<<endl;
    if (opt.iseparatefiles) cout<<"Separate files output"<<endl;
    if (opt.iextendedoutput) cout<<"Extended output for particle extraction from input files"<<endl;
    if (opt.iHaloCoreSearch) cout<<"Searching for 6dfof cores so as to disentangle mergers"<<endl;
    if (opt.iHaloCoreSearch && opt.iAdaptiveCoreLinking) cout<<"With adaptive linking lengths"<<endl;
    if (opt.iHaloCoreSearch && opt.iPhaseCoreGrowth) cout<<"With with phase-space tensor core assignment"<<endl;
    if (opt.inputtype==IOGADGET) cout<<"Gadget file particle input "<<endl;
    else if (opt.inputtype==IOTIPSY) cout<<"Tipsy file particle input "<<endl;
    else if (opt.inputtype==IORAMSES) cout<<"RAMSES file particle input "<<endl;
#ifdef USEHDF
    else if (opt.inputtype==IOHDF) cout<<"HDF file particle input "<<endl;
#endif
#ifdef USEXDR
    else if (opt.inputtype==IONCHILADA) cout<<"NCHILADA file particle input "<<endl;
#endif
#ifdef USEMPI
    cout<<"MPI input communication requires allocation of temporary particle buffer of size "<<opt.mpiparticletotbufsize/1024./1024./1024.<<" GB."<<endl;
#endif
    cout<<" -------------------------- "<<endl;
#ifdef USEMPI
    }
#endif
#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}
