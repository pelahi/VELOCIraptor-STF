.. _usage:

Using **VELOCIraptor**
######################

An example configuration file can be found the examples directory within the repository
(see for instance :download:`sample <../examples/sample.cfg>`). This sample file lists
all the options. *Only the keywords listed here will be used, all other words/characters
are ignored*. One can check the options used by examining **foo.configuration**, where **foo** is
your base output filename.

.. warning:: Note that if misspell a keyword it will not be used.
.. warning:: Since this file is always written **DO NOT** name your input configuration file **foo.configuration**.


There are numerous key words that can be passed. Here we list them

.. topic:: Output related

    * Output base name. Overrides the name passed with the command line argument **-o**. Only implemented for completeness.
        ``Output = filename``
    * A filename for storing the intermediate step of calculating local densities. This is particularly useful if the code is not compiled with **STRUCDEN** & **HALOONLYDEN** (see :ref:`getting`:`Compilation Options`_).
        ``Output_den = filename``
    * Flag indicating whether separate files are written for field and subhalo groups.
        ``Separate_output_files = 1/0``
    * Flag indicating whether to producing a file which lists for every particle the group they belong to. Can be used with **tipsy** format or to tag every particle.
        ``Write_group_array_file = 1/0``
    * Integer indicating whether output is hdf (2), binary (1), ascii (0) or adios (3).
        ``Binary_output = 3/2/1/0``
    * Flag indicating whether to calculate/output even more halo properties.
        ``Extensive_halo_properties_output = 1/0``
    * Flag indicating whether produce extended output for quick particle extraction from input catalog of particles in structures
        ``Extended_output = 1/0``
    * Flag indicating whether the properties output is in physical or comoving little h units.
        ``Comoving_units = 1/0``

.. topic:: Input related

    * Flag indicating that input simulation is cosmological or not. With cosmological input, a variety of length/velocity scales are set to determine such things as the virial overdensity, linking length.
        ``Cosmological_input = 1/0``
    * Amount of information to read from input file in one go (100000).
        ``Input_chunk_size = 100000``
        \arg <b> \e HDF_name_convention </b> HDF dataset naming convection. See \ref hdfitems.h for what naming conventions are available and what names exist. Currently have \ref HDFNUMNAMETYPES. \ref Options.ihdfnameconvention \n
        \arg <b> \e Input_includes_star_particle </b> If star particle specific information is in the input file. \ref Options.iusestarparticles \n
        \arg <b> \e Input_includes_bh_particle </b> If bh/sink particle specific information is in the input file. \ref Options.iusesinkparticles \n
        \arg <b> \e Input_includes_wind_particle </b> If wind particle specific information is in the input file. \ref Options.iusewindparticles \n
        \arg <b> \e Input_includes_tracer_particle </b> If tracer particle specific information is in the input file. \ref Options.iusetracerparticles \n
        \arg <b> \e Input_includes_star_particle </b> If star particle specific information is in the input file. \ref Options.iusestarparticles \n

.. topic:: Parameters related to search type.

    * An integer describing what types of particles are searched. Options are:
        - **1** *All* particles are searched
        - **2** *DarkMatter* particles (which are typically defined as type 1,2,3 for gadget) are searched
        - **3** Star particles (which are typically defined as type 4 for gadget) are searched
        - **4** Gas particles (which are typically defined as type 0 for gadget) are searched
        ``Particle_search_type = 1/2/3/4``
    * An integer indicating gas/stellar search done separately from DM search.
        - **1** field search run as normal and then substructure search for baryons run using baryons identified in field search.
        - **2** field search also altered to treat baryons differently, allowing only DM particles to be used as head links (ie link dm-dm, dm-baryon, but not baryon-baryon nor baryon-dm). Then DM substructure search with baryons associated to closest DM particle in phase-space.
        - **0** do nothing special for baryon particles.
        ``Baryon_searchflag = 0/1/2``
    * Flag indicating whether field objects are searched for internal substructures. Default is 1 (on)
        ``Search_for_substructure = 1/0``
    * Flag indicates that no field search is going to be run and the entire volume will be treated as a background region (halo). Useful if searching for substructures in non-cosmological simulations. But can also be co-opted for other searches using different outlier criteria and FOF algorithms
        ``Singlehalo_search_search = 0/1``

.. topic:: Parameters related to field (halo) search.

    * An integer indicating what type of field search is run. There are several
        - **5** standard 3D FOF based algorithm
        - **4** standard 3D FOF based algorithm :strong:`FOLLOWED` by 6D FOF search using the velocity scale defined by the largest halo on particles in 3DFOF groups
        - **3** standard 3D FOF based algorithm :strong:`FOLLOWED` by 6D FOF search using :emphasis:`adaptive` velocity scale for each 3DFOF group on particles in these groups.
        ``FoF_Field_search_type = 5/4/3``
    * Multiplicative factor of order unity that allows one to use different physical linking lengths between field objects and substructures. :strong:`Note`: substructure search defines the base linking length via ```Physical_linking_length``. Typically for standard 3DFOF searches of dark matter haloes, set to 2.0, as typical base linking length is 0.1 times the interparticle spacing when examining cosmological simulations.
        ``Halo_linking_length_factor = 2.0``
    * Multiplicative factor of order unity for the dispersions used in 6D searches. Typical values are order unity as velocity dispersions are used to define the velocity linking length scale.
        ``Halo_velocity_linking_length_factor =``
    * Multiplicative factor of order unity that allows one to use different configuration space linking lengths between 3DFOF and 6DFOF field search. Typically this is 1.0
        ``Halo_6D_linking_length_factor = 1.0``
    * Multiplicative factor of order unity scaling applied to dispersions used in 6DFOF field search. Typical values are 1.25.
        ``Halo_6D_vel_linking_length_factor = 1.25``
    * Flag that keeps the 3DFOF if field 6DFOF search is done. This is typically invoked when searching for galaxies as the 3DFOF can be interpreted as the inter halo stellar mass and 6DFOF galaxies.
        ``Keep_FOF = 0/1``
    * Integer that allows field objects (or so-called halos) to require a different minimum size than all other substructures. Ignored if not passed or <0, with halos defaulting to ``Minimum_size`` value.
        ``Minimum_halo_size =-1``


.. topic:: Parameters related to local density estimator used to identify particles in substructures

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

\arg <b> \e Halo_core_search </b> 0/1/2 flag allows one to explicitly search for large 6D FOF cores that are indicative of a recent major merger. Since substructure is defined on the scale of the maximum cell size and major mergers typically result two or more phase-space dense regions that are \e larger than the cell size used in reasonable substructure searches, one can identify them using this search. The overall goal is to treat these objects differently than a substructure. However, if 2 is set, then smaller core is treated as substruture and all particles within the FOF envelop are assigned to the cores based on their phase-space distance to core particles \ref Options.iHaloCoreSearch \n
\arg <b> \e Use_adaptive_core_search </b> 0/1 flag allows one to run complex adaptive phase-space search for large 6D FOF cores and then use these linking lengths to separate mergers. 0 is simple high density dispersively cold cores with v scale adaptive, 1 is adaptive with simple dispersions in both x and v.
\arg <b> \e Use_phase_tensor_core_growth </b> 0/1 flag allows one to run complex phase-space growth of merger remnants (6D FOF cores found). 0 is assignment with simple x and v dispersion to nearest core particle, 1 is phase-space tensor distance assignemnt to CM of core.
\arg <b> \e Halo_core_ellx_fac </b> scaling applied to linking length when identifying merger remnants. Typically values are \f$ \sim0.5 \f$  \ref Options.halocorexfac
\arg <b> \e Halo_core_ellv_fac </b> scaling applied to local dispersion to define the velocity scale used to identify merger remnants. Typically values are \f$ \sim1 \f$  \ref Options.halocorevfac
\arg <b> \e Halo_core_ncellfac </b> used to determine the minimum number of particles a merger remnants is composed of, specifically \f$ N_{\rm min}= f_{\rm ncell}* N_{\rm S} \f$. Typically values are     \arg <b> \e Halo_core_adaptive_sigma_fac </b> used when running fully adaptive core search, specifies the width of the physical linking length in configuration space dispersion (think of this as how many sigma to include). Typically values are \f$ \sim2 \f$. This has been tested on hydrodynamnical simulations to separate galaxy mergers. \ref Options.halocoresigmafac
\arg <b> \e Halo_core_num_loops </b> allows the core search to iterate, shrinking the velocity linking length to used till the number of cores identified decreases or this limit is reached. Allows apative search with larger linking length to be robust.  Typically values are \f$ \sim5 \f$ with loops only running often twice. \ref Options.halocorenumloops.
\arg <b> \e Halo_core_loop_ellx_fac </b> Factor by which configuration linking length is decreased when running loops for core search.  Typically values are \f$ \sim0.75 \f$. \ref Options.halocorexfaciter
\arg <b> \e Halo_core_loop_ellv_fac </b> Factor by which velocity linking length is decreased when running loops for core search.  Typically values are \f$ \sim0.75 \f$. \ref Options.halocorevfaciter
\arg <b> \e Halo_core_loop_elln_fac </b> Factor by which min group size is changed when running loops for core search.  Typically values are \f$ \sim0.25 \f$. \ref Options.halocorenumfaciter

\subsection fofotherconfig Other modifiers for search
\arg <b> \e CM_refadjustflag </b> 0/1 flag indicates that one moves to the CM frame of the structure to search for substructures.

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

\section cosmoconfig Units & Cosmology
\subsection unitconfig Units
\arg <b> \e Length_unit </b> change the input length unit \ref Options.L \n
\arg <b> \e Velocity_unit </b> change the input velocity unit \ref Options.V \n
\arg <b> \e Mass_unit </b> change the input mass unit \ref Options.M \n
\arg <b> \e Gravity </b> change the gravity unit. should be set such that \f$ v^2=GM/r \f$.   \ref Options.G \n
\arg <b> \e Hubble_unit </b> change the value of Hubble expansion (from normal 100 km/s/Mpc to units used, that is velocity unit/length unit) \ref Options.G \n
\arg <b> \e Mass_value </b> if not mass is stored using the option \b NOMASS (see \ref STF-makeflags) then this is the mass of the particles \ref Options.MassValue \n
\arg <b> \e Length_unit_to_kpc </b> specify the desired return unit in kpc \ref Options.lengthtokpc \n
\arg <b> \e Velocity_unit_to_kms </b> specify the desired return unit in kms \ref Options.velocitytokms \n
\arg <b> \e Mass_unit_to_solarmass </b> specify the desired return unit in solar masses \ref Options.masstosolarmass \n

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

\section inputflags input flags related to varies input formats
\arg <b> \e NSPH_extra_blocks </b> If gadget snapshot is loaded one can specific the number of extra <b> SPH </b> blocks are read/in the file. \ref Options.gnsphblocks \n
\arg <b> \e NStar_extra_blocks </b> If gadget snapshot is loaded one can specific the number of extra <b> Star </b> blocks are read/in the file. \ref Options.gnstarblocks \n
\arg <b> \e NBH_extra_blocks </b> If gadget snapshot is loaded one can specific the number of extra <b> Black hole </b> blocks are read/in the file. \ref Options.gnbhblocks \n


\section mpiconfigs MPI specific options
\arg <b> \e MPI_part_allocation_fac </b> Memory allocated in mpi mode to store particles is (1+factor)* the memory need for the initial mpi decomposition.
This factor should be >0 and is mean to allow a little room for particles to be exchanged between mpi threads withouth having to require new memory allocations and copying
of data. \ref Options.mpipartfac \n
\arg <b> \e MPI_particle_total_buf_size </b> Total memory size in bytes used to store particles in temporary buffer such that
particles are sent to non-reading mpi processes in one communication round in chunks of size buffer_size/NProcs/sizeof(Particle). \ref Options.mpiparticlebufsize \n


Using **TreeFrog**
##################
