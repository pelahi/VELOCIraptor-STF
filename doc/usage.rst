.. _usage:

Using **VELOCIraptor**
######################

.. _running:

Running the code
================

**Velociraptor** is a stand alone executable (with the executable named stf (or STructure Finder for historical reasons)).
It can be run in serial, with OpenMP, or MPI APIs. A typical command to start the code looks like:
::

 ./stf < args >

When compiled with OpenMP, setting the environment variable ``OMP_NUM_THREADS`` will set the number of threads in the openmp sections.

With MPI using 8 MPI threads:
::

 mpirun -np 8 ./stf < args >

where here we assume that the parallel
environment uses the ``mpirun`` command to start MPI
applications. Depending on the operating system, other commands may be
required for this task, e.g. ``srun`` on some Cray machines. Note that
the code can in principle be started using an arbitrary number of
mpi threads, but the mpi decomposition is most efficient for powers of 2.

The output produced by VELOCIraptor will typically consist of several files containing:
bulk properties of structures found; particles belonging to these structures; and several
additional files containing configuration information.

When running in MPI, currently each mpi thread writes its own output unless
the code has been compiled with a parallel HDF5 library and HDF5 output is requested.
In that case, a single file is written containing data from all threads for each type
of output requested.

.. note:: At the moment, mpirun assumes that a single structure can fit onto the memory local to the mpi thread. If larger field objects (haloes) are to beanalyzed such that they are unlikely to fit into local memory, it is suggested another machine be used. Revision is in the works to use the Singlehalo_search option after field halos have been identified.

.. _cmdargs:

Arguments
---------

The code has several command line arguements. To list the arguments, type
::

    ./stf -?

.. topic:: The arguments that can be passed are:

    **-i** ``< file name of input file >``

    **-s** ``< number of files over which input is split >``

    **-I** ``< input format [1 Gadget, 2 HDF5, 3 Tipsy, 4 RAMSES, 5 NCHILADA] >``

    **-Z** ``< number of files to read in parallel (when mpi is invoked) >``

    **-o** ``< output base name (this can be overwritten by a configuration option in the config file. Suggestion would be to not use this option in the config file, use explicit command>``

    **-C** ``< configuration file name (see`` :ref:`configoptions` ``) >``

Of these arguements, only an input file and an output name must be provided.
In such a case, it is assumed that there is only 1 input file, 1 read mpi thread,
and default values for all other confirguration options. We suggest you do NOT run the code in this fashion.
Instead we suggest the code be run with at least a configuration file passed.
::

    ./stf -i input -o output -C configfile.txt

This configuration file is an ascii file that lists keywords and values.
A list of keywords, along with a description is presented below in :ref:`configoptions`.
A more typical command for a large cosmological simulation might be something like
::

    export OMP_NUM_THREADS=4
    mpirun -np 64 ./stf -i somehdfbasename -s 128 -I 2 -Z 64 -o output -C configfile.txt > stf.log

.. _swiftintegration:

Running within swiftsim
-----------------------

|vr| is also able to be called from within an N-body/Hydrodynamnical code as a library.
Currently the code has been integrated in to **swifsim**. Details can be found
in the **swiftsim** `documentation <http://icc.dur.ac.uk/swift/docs/index.html>`_.
The key is that the **swiftsim** code's configuration file  lists the |vr| configuration file
used to run |vr|.

.. _briefoutput:

Output
------

Here we provide a *brief* description of the standard data products provided by |vr|.
For a more detailed discussion and some sample analysis using these data products see :ref:`output`.

When operating in a typical configuration with typical compile time options,
the executable (or each mpi thread) will produce several files (with the mpi
threads appending their rank to the end of the file name, unless parallel HDF5 output is requested).
The files typically produced are :

.. topic:: Output files

    * ``.properties``: a file containing the bulk properties of all structures identified.
    * ``.catalog_groups``: a file containing the size of the structures (in number of particles associated) & information need to read particle information produced by velociraptor
    * ``.catalog_particles``: a file containing a list of particle IDs of those in structures. Information contained in ``.catalog_groups`` is used to parse this data.
    * ``.catalog_particles.unbound``: similar to ``catalog_particles`` but lists particles in structures but are formally unbound. Information contained in ``.catalog_groups`` is used to parse this data.

.. topic:: Extra output files

    * ``.profiles``: a file containing the radial mass profiles of (sub)halos
    * ``.catalog_parttypes``: a file similar to ``.catalog_particles`` but stores particle type instead of paricle id.
    * ``.catalog_parttypes.unbound``: a file similar to ``.catalog_parttypes`` but for unbound particles.
    * ``.extendedinfo``: a file containing extra information on where particles are located in the input file for quick extraction from said input file of particles within groups. Still in alpha
    * ``.catalog_SOlist``: a file containing particle IDs within the spherical overdensity region of halos.

.. _configoptions:

Configuration File
------------------

An example configuration file can be found the examples directory within the repository
(see for instance :download:`sample <../examples/sample_dmcosmological_run.cfg>`).
This sample file lists all the options. *Only the keywords listed here will be used, all other words/characters are ignored*. One can check the options used by examining **foo.configuration**, where **foo** is your base output filename.

We suggest the following files as a basis:
    * :download:`N-body simulations configuration <../examples/sample_hydrocosmological_run.cfg>`
    * :download:`Hydro simulations configuration <../examples/sample_hydrocosmological_run.cfg>`
    * :download:`SWIFT N-body simulation configuration <../examples/sample_swiftdm_3dfof_subhalo.cfg>`
    * :download:`SWIFT Hydro simulation configuration <../examples/sample_swifthydro_3dfof_subhalo_extra_properties.cfg>`

.. warning:: Note that if misspell a keyword it will not be used.
.. warning:: Since this file is always written **DO NOT** name your input configuration file **foo.configuration**.

There are numerous key words that can be passed. Here we list them, grouped into several categories:

    - :ref:`IO <config_io>`

        - :ref:`Inputs <config_input>`
        - :ref:`Outputs <config_output>`

    - :ref:`Parameters related to type of search <config_search>`

        - :ref:`Field search <config_field_search>`
        - :ref:`Substructure search <config_sub_search>`
        - :ref:`Local Velocity Density <config_local_vden>`
        - :ref:`Core search <config_core_search>`

    - :ref:`Unbinding <config_unbinding>`
    - :ref:`Properties <config_properties>`
    - :ref:`Units/Cosmology <config_siminfo>`

        - :ref:`Units <config_units>`
        - :ref:`Cosmology <config_cosmology>`

    - :ref:`Parallel <config_parallel>`

        - :ref:`MPI <config_mpi>`
        - :ref:`OpenMP <config_openmp>`

    - :ref:`Miscellaneous <config_misc>`


.. _config_io:

I/O
^^^

Input and output related options

.. _config_input:

..topic:: Input related

    ``Cosmological_input = 1/0``
        * Flag indicating that input simulation is cosmological or not. With cosmological input, a variety of length/velocity scales are set to determine such things as the virial overdensity, linking length.
    ``Input_chunk_size = 100000``
        * Amount of information to read from input file in one go (100000).
    ``HDF_name_convention =``
        * Integer describing HDF dataset naming convection. Currently implemented values can be found in :ref:`subsection_hdfnames`.
    ``Input_includes_dm_particle = 1/0``
        * Flag indicating whether file contains dark matter/N-body particles in input file.
    ``Input_includes_gas_particle = 1/0``
        * Flag indicating whether file contains gas particles in input file.
    ``Input_includes_star_particle = 1/0``
        * Flag indicating whether file contains star particles in input file.
    ``Input_includes_bh_particle = 1/0``
        * Flag indicating whether file contains black hole particles in input file.
    ``Input_includes_wind_particle = 1/0``
        * Flag indicating whether file contains wind particles in input file.
    ``Input_includes_tracer_particle = 1/0``
        * Flag indicating whether file contains tracer particles in input file.
    ``Input_includes_extradm_particle = 1/0``
        * Flag indicating whether file contains extra (low resolution) N-body particles in input file from a zoom simulation.
    Gas related input
        ``Gas_internal_property_names = ,``
            * Comma separated list of strings listing extra gas properties to be read from HDF file for which bulk mean/total properties are calculated for objects. Useful way of passing properties like molecular H2 fraction, etc.
        ``Gas_chemistry_names = ,``
            * Comma separated list of strings listing extra chemical properties to be read from HDF file for which bulk mean/total properties are calculated for objects. Useful way of passing properties like molecular H2 fraction, etc.
        ``Gas_chemistry_production_names = ,``
            * Comma separated list of strings listing extra production channels for metals to be read from HDF file for which bulk mean/total properties are calculated for objects. Useful way of passing properties like molecular H2 fraction, etc.
    Star related input
        ``Star_internal_property_names = ,``
            * Comma separated list of strings listing extra star properties to be read from HDF file for which bulk mean/total properties are calculated for objects. Useful way of passing properties like molecular H2 fraction, etc.
        ``Star_chemistry_names = ,``
            * Comma separated list of strings listing extra chemical properties to be read from HDF file for which bulk mean/total properties are calculated for objects. Useful way of passing properties like molecular H2 fraction, etc.
        ``Star_chemistry_production_names = ,``
            * Comma separated list of strings listing extra production channels for metals to be read from HDF file for which bulk mean/total properties are calculated for objects. Useful way of passing properties like molecular H2 fraction, etc.
    Black hole related input
        ``BH_internal_property_names = ,``
            * Comma separated list of strings listing extra black properties to be read from HDF file for which bulk mean/total properties are calculated for objects. Useful way of passing properties like molecular H2 fraction, etc.
        ``BH_chemistry_names = ,``
            * Comma separated list of strings listing extra chemical properties to be read from HDF file for which bulk mean/total properties are calculated for objects. Useful way of passing properties like molecular H2 fraction, etc.
        ``BH_chemistry_production_names = ,``
            * Comma separated list of strings listing extra production channels for metals to be read from HDF file for which bulk mean/total properties are calculated for objects. Useful way of passing properties like molecular H2 fraction, etc.
    Extra DM related input
        ``Extra_dm_internal_property_names = ,``
            * Comma separated list of strings listing extra dm properties to be read from HDF file for which bulk mean/total properties are calculated for objects. Useful for modified dark matter simulations, such as annihilating and self-interactive dark matter.
    Gadget related input
        ``NSPH_extra_blocks =``
            * Integer inticading  the number of extra **SPH** blocks are read in the file if gadget input.
        ``NStar_extra_blocks =``
            * Integer inticading  the number of extra **star** blocks are read in the file if gadget input.
        ``NBH_extra_blocks =``
            * Integer inticading  the number of extra **BH** blocks are read in the file if gadget input.

.. _config_output:

.. topic:: Output related

    ``Output = filename``
        * Output base name. Overrides the name passed with the command line argument **-o**. Only implemented for completeness.
    ``Output_den = filename``
        * A filename for storing the intermediate step of calculating local densities. This is particularly useful if the code is not compiled with **STRUCDEN** & **HALOONLYDEN** (see :ref:`compileoptions`).
    ``Separate_output_files = 1/0``
        * Flag indicating whether separate files are written for field and subhalo groups.
    ``Write_group_array_file = 1/0``
        * Flag indicating whether to producing a file which lists for every particle the group they belong to. Can be used with **tipsy** format or to tag every particle.
    ``Binary_output = 2/1/0``
        * Integer flag indicating type of output.
            - **2** self-describing binar format of HDF5. **Recommended**.
            - **1** raw binary.
            - **0** ASCII.
    ``Extended_output = 1/0``
        * Flag indicating whether produce extended output for quick particle extraction from input catalog of particles in structures
    ``Spherical_overdensity_halo_particle_list_output = 1/0``
        * Flag indicating whether particle IDs identified within the spherical overdensity of field halos is written (to a .catalog_SOlist). Useful if looking at evolution of particles within spherical overdensities.
    ``Sort_by_binding_energy = 1/0``
        * Flag indicating whether particle IDs written in .catalog_particles are sorted by binding energy (1) or potential energy (0).
    ``No_particle_ID_list_output = 1/0``
        * Flag indicating whether particle IDs written (i.e., write the .catalog_\* files). Default is 1. Particle ID files are necessary for constructing merger trees but if just properties of (sub)halos, then turn off.

.. _config_search:

Searching for Structures
^^^^^^^^^^^^^^^^^^^^^^^^

Options related to searching for (sub)halos. General search parameters set particles to be search and the overall type of search.

    ``Particle_search_type = 1/2/3/4``
        * An integer describing what types of particles are searched. A full list of options is in :ref:`subsection_searchtypes`. Typical options are:
            - **1** *All* particles are searched
            - **2** *DarkMatter* particles (which are typically defined as type 1,2,3 for gadget) are searched
            - **3** Star particles (which are typically defined as type 4 for gadget) are searched
            - **4** Gas particles (which are typically defined as type 0 for gadget) are searched
    ``Baryon_searchflag = 0/1/2``
        * An integer indicating gas/stellar search done separately from DM search.
            - **2** field search also altered to treat baryons differently, allowing only DM particles to be used as head links (ie link dm-dm, dm-baryon, but not baryon-baryon nor baryon-dm). Then DM substructure search with baryons associated to closest DM particle in phase-space. **Recommended**.
            - **1** field search run as normal and then substructure search for baryons run using baryons identified in field search.
            - **0** do nothing special for baryon particles.
    ``Search_for_substructure = 1/0``
        * Flag indicating whether field objects are searched for internal substructures. Default is 1 (on)
    ``Singlehalo_search_search = 0/1``
        * Flag indicates that no field search is going to be run and the entire volume will be treated as a background region (halo). Useful if searching for substructures in non-cosmological simulations. But can also be co-opted for other searches using different outlier criteria and FOF algorithms

.. _config_field_search:

.. topic:: Parameters related to field (halo) search

    ``FoF_Field_search_type = 5/4/3``
        * An integer indicating what type of field search is run. There are several
            - **5** standard 3D FOF based algorithm
            - **4** standard 3D FOF based algorithm :strong:`FOLLOWED` by 6D FOF search using the velocity scale defined by the largest halo on particles in 3DFOF groups
            - **3** standard 3D FOF based algorithm :strong:`FOLLOWED` by 6D FOF search using :emphasis:`adaptive` velocity scale for each 3DFOF group on particles in these groups.
    ``Halo_3D_linking_length = 0.2``
        * Linking length used to find configuration space 3D FOF halos. If cosmological file then assumed to be in units of inter particle spacing, if loading in a single halo then can be based on average interparticle spacing calculated, otherwise in input units. Default is 0.2 in interpaticle spacing units.
    ``Halo_velocity_linking_length_factor = 1.0``
        * Multiplicative factor of order unity for the dispersions used in 6D searches. Typical values are order unity as velocity dispersions are used to define the velocity linking length scale.
    ``Halo_6D_linking_length_factor = 1.0``
        * Multiplicative factor of order unity that allows one to use different configuration space linking lengths between 3DFOF and 6DFOF field search. Typically this is 1.0
    ``Halo_6D_vel_linking_length_factor = 1.25``
        * Multiplicative factor of order unity scaling applied to dispersions used in 6DFOF field search. Typical values are 1.25.
    ``Keep_FOF = 0/1``
        * Flag that keeps the 3DFOF if field 6DFOF search is done. This is typically invoked when searching for galaxies as the 3DFOF can be interpreted as the inter halo stellar mass and 6DFOF galaxies.
    ``Minimum_halo_size =-1``
        * Integer that allows field objects (or so-called halos) to require a different minimum size than all other substructures. Ignored if not passed or <0, with halos defaulting to ``Minimum_size`` value.

.. _config_sub_search:

.. topic:: Parameters related to substructure search

    **Note**: default values are fine and typically do not need to be set in the configuration file. Exception would be Minimum_size

    ``FoF_search_type = 1``
        * Integer indicating what type of FOF algorithm to use. Several substructure FOF criteria are implemented (see :ref:`subsection_foftypes` for complete list). Suggested value is **1**, the standard phase-space based, well tested VELOCIraptor criterion.
    ``Outlier_threshold = 2.5``
        * Threshold of sigma level of outliers to be searched which should be order unity but > 1 (default is 2.5)
    ``Significance_level = 1.0``
        * Minimum significance level of a substructure which should be order unity (default is 1)
    ``Velocity_ratio = 2.0``
        * Speed ratio used in linking particles which should be order unity and > 1 (default is 2)
    ``Velocity_opening_angle = 0.10``
        * Angle between velocities when linking (in units of :math:`\pi`) (default is 0.10)
    ``Substructure_physical_linking_length = 0.1``
        * Physical linking length used in phase-space substructure FOF. If cosmological file then assumed to be in units of inter particle spacing, if loading in a single halo then can be based on average interparticle spacing calculated, otherwise in input units. Default is 0.1 in interpaticle spacing units.
    ``CMrefadjustsubsearch_flag = 1/0``
        * Flag indicating whether particles are moved to the rough CM velocity frame of the background before substructures are searched for (default is on)
    ``Iterative_searchflag = 1/0``
        * Flag to use interactive substructure search which is designed to first identify spatially compact candidate outlier regions and then relaxes the criteria to find the more diffuse (in phase-space) regions associate with these candidate structures (default is on)
    ``Iterative_linking_length_factor = 2.0``
        * Factor multiplied with linking length when using iterative method and identifying outlier regions associated with the initial candidate list of spatially compact outlier groups. Typical values are order ``Halo_linking_length_factor`` (2.0)
    ``Iterative_threshold_factor = 1.0``
        * Factor multiplied with threshold when using iterative method and identifying outlier regions associated with the initial candidate list of spatially compact outlier groups. Typical values are order unity.
    ``Iterative_Vratio_length_factor = 1.0``
        * Factor multiplied with speed ratio when using iterative method and identifying outlier regions associated with the initial candidate list of spatially compact outlier groups. Typical values are order unity.
    ``Iterative_ThetaOp_length_factor = 1.0``
        * Factor multiplied with opening angle when using iterative method and identifying outlier regions associated with the initial candidate list of spatially compact outlier groups. Typical values are order unity.
    ``Minimum_size = 20``
        * Minimum number of particles in a (sub)structure (default is 20).

.. _config_local_vden:

.. topic:: Configuration for local density calculation used to identify substructures

    **Note**: default values are fine and typically do not need to be set in the configuration file.

    ``Local_velocity_density_approximate_calculation = 2/1/0``
        * Flag indicating how to calculate computationally expensive local velocity densities.
            - **2** approximative search limited to particles in halos (requires no mpi communication). **Recommended**.
            - **1** approximative search, group particles in leaf nodes of tree
            - **0** full search per particle.
    ``Nsearch_velocity = 32``
        * Number of velocity neighbours used to calculate velocity density (suggested value is 32)
    ``Nsearch_physical = 32``
        * Number of physical neighbours searched to calculate velocity density (suggested value is 256)
    ``Cell_fraction = 0.1``
        * Fraction of a halo contained in a subvolume used to characterize the background (suggested value is 0.01)
    ``Grid_type = 1``
        * Integer describing type of grid used to decompose volume for substructure search (suggested value is 1)
            - **1** standard physical shannon entropy, balanced KD tree volume decomposition into cells. **Recommended**
            - **2** phase phase-space shannon entropy, balanced KD tree volume decomposition into cells
            - **3** simple simple physical balanced KD tree decomposition of volume into cells

.. _config_core_search:

.. topic:: Configuration for core search and growth.

    This either identifies major mergers in DM simulations or used to find galaxies when searching for stars.

    ``Halo_core_search = 0/1/2``
        * Integer allows one to explicitly search for large 6D FOF cores that are indicative of a recent major merger. Since substructure is defined on the scale of the maximum cell size and major mergers typically result two or more phase-space dense regions that are *larger* than the cell size used in reasonable substructure searches, one can identify them using this search. The overall goal is to treat these objects differently than a substructure. However, if 2 is set, then smaller core is treated as substructure and all particles within the FOF envelop are assigned to the cores based on their phase-space distance to core particles.
            - **2** search for cores and growth them. **Recommended**.
            - **1**
            - **0** do not search cores.
    ``Use_adaptive_core_search = 0/1``
        * Flag allows one to run complex adaptive phase-space search for large 6D FOF cores and then use these linking lengths to separate mergers. 0 is simple high density dispersively cold cores with velocity scale adaptive, 1 is adaptive in both configuration & velocity.
    ``Use_phase_tensor_core_growth = 0/1``
        * Flag allows one to run complex phase-space growth of merger remnants (6D FOF cores found). 0 is assignment with simple x and v dispersion to nearest core particle, 1 is phase-space tensor distance assignemnt to CM of core.
    ``Halo_core_ellx_fac =``
        * Factor applied to linking length when identifying merger remnants. Typically values are 0.5
    ``Halo_core_ellv_fac =``
        * Factor applied to local dispersion to define the velocity scale used to identify merger remnants. Typically values are order unity
    ``Halo_core_ncellfac = 0.005``
        * Factor used to determine the minimum number of particles a merger remnants is composed of using number of particles in the halo times this factor. For DM typically values are 0.005.
    ``Halo_core_adaptive_sigma_fac = 2.0``
        * Factor used when running fully adaptive core search, specifies the width of the physical linking length in configuration space dispersion (think of this as how many sigma to include). Typically values are 2. This has been tested on hydrodynamnical simulations to separate galaxy mergers.
    ``Halo_core_num_loops = 10``
        * Allows the core search to iterate, shrinking linking lengths used till the number of cores identified reaches zero or this limit is reached. Allows apative search with larger linking length to be robust.  Typically values are 10, though typically loops run twice.
    ``Halo_core_loop_ellx_fac = 0.75``
        * Factor by which configuration linking length is decreased when running loops for core search.  Typically values are 0.75
    ``Halo_core_loop_ellv_fac = 1.0``
        * Factor by which velocity linking length is decreased when running loops for core search.  Typically values are 1.
    ``Halo_core_loop_elln_fac = 1.2``
        * Factor by which min group size is changed when running loops for core search.  Typically values are order unity & > 1.
    ``Halo_core_phase_significance = 2.0``
        * Significance a core must be in terms of phase-space distance scaled by dispersions (sigma). Typical values are order unity & > 1.

.. topic:: Configuration for cleaning up substructuers that overlap in phase-space.

    Substructures can be merged together if they overlap in phase space.

    ``Structure_phase_merge_dist = 0.25``
        * Phase-distance normalised by dispersions below which structures are merged together. Typical valuse are < 1.
    ``Apply_phase_merge_to_host = 1``
        * Flag whether to also check substructures can be merged with the host background. 1 is on.


.. _config_unbinding:

Unbinding
^^^^^^^^^

Particles in strutures can be checked to see if they are bound relative to a kinetic reference frame (CM of the structure).
This cleans the (sub)structures of spurious objects and particles.

    ``Unbind_flag = 1/0``
        * Flag indciating whether substructures passed through an unbinding routine.
    ``Unbinding_type = 1/0``
        * Integer setting the unbinding criteria used. Either just remove particles deemeed "unbound" (**1**), that is those with :math:`\alpha T+W>0` given by ``Allowed_kinetic_potential_ratio``, or (**0**) additionally removes "unbound" and least bound particles till system also has a true bound fraction > ``Min_bound_mass_frac``.
    ``Allowed_kinetic_potential_ratio =``
        * Ratio of kinetic to potential energy at which a particle is still considered bound, ie: particle is still bound if :math:`\alpha T+W<0`, so :math:`\alpha=1` would be standard unbinding and :math:`\alpha<1` allows one to identify unbound tidal debris. Given that **VELOCIraptor** was designed to identify tidal streams, it makes little sense to have this set to 1 unless explicitly required. Note that the code still separates particles into bound and unbound. Values of :math:`\alpha\geq 0.2` seems to minimize the number of false positives in tidal debris while still identifying completely unbound tidal debris.
    ``Min_bound_mass_frac =``
        * Minimum fraction of particles that must be self-bound. If interested in identifying tidal debris, ues values of 0.2, for self-bound substructures, use :math:`\gtrsim 0.5`
    ``Bound_halos = 0/1/2``
        * Integer that ignores the boundness of field structures (haloes) (**0**), checks if they are self bound only before (**1**) or also after (**2**) substructures have been identified and extracted from the halo. Demanding boundness after substructure search can have interesting consequences as it is possible that a multiple merger will appear as a single FOF halo, however all with all the cores removed, the FOF halo is actually an unbound structure.
    ``Keep_background_potential = 1/0``
        * Flag indicating whether while checking if a structure is bound, to treat the candidate structure in isolation, updating the potential continuously, or leave the background potential.  background sea. When finding tidal debris, it is useful to keep the background. \ref Options.uinfo & \ref UnbindInfo.bgpot \n
    ``Kinetic_reference_frame_type = 0/1``
        * Integer that sets the kinetic frame when determining whether particle is bound. Default is to use the centre-of-mass velocity frame (0) but can also use region around minimum of the potential (1).
    ``Min_npot_ref = 10``
        * The minimum number of particles used to calculate the velocity of the minimum of the potential (default is 10).
    ``Frac_pot_ref = 0.1``
        * Fraction of particles used to calculate the velocity of the minimum of the potential (0.1). If smaller than ``Min_npot_ref``, that is used.
    ``Unbinding_max_unbound_removal_fraction_per_iteration = 0.5``
        * Maximum fraction of unbound particles removed per iteration in unbinding process.
    ``Unbinding_max_unbound_fraction = 0.95``
        * Maximum fraction of particles that can be considered unbound before group removed entirely and is not processed iteratively.
    ``Unbinding_max_unbound_fraction_allowed = 0.005``
        * Maximum fraction of unbound particles allowed after unbinding. If set to zero, all unbound particles removed.


.. _config_properties:

Properties
^^^^^^^^^^

Configuration options related to the bulk properties calculated.

    ``Inclusive_halo_mass = 3/2/1/0``
        * Flag indicating whether inclusive masses are calculated for field objects.
            - **3** indicates inclusive SO masses are calculated after substructure is found.
            - **2** indicates inclusive SO masses are calculated before substructure is found.
            - **1** indicates inclusive SO masses are calculated before substructure is found but limited to particles in the halo.
            - **0** indicates masses exclusive.
    ``Iterate_cm_flag = 0``
        * Flag indicating whether to iteratively find the centre-of-mass of an object (1) or simply deterine bulk centre of mass and centre of mass velocity (0). Calculation is based on all particles exclusively belonging to the object.
    ``Reference_frame_for_properties = 2``
        * Flag indicating what reference position to use when calculating radially dependent properties.
            - **2** use the position of the particle with the minimum potential.
            - **1** use the position of the most bound particle.
            - **0** use the centre-of-mass.
    ``Extensive_halo_properties_output = 1``
        * Flag indicating that one should calculate more properties for objects, such as angular momentum in spherical overdensity apertures.
    ``Extensive_gas_properties_output = 1``
        * Flag indicating that in addition to calculating extra halo properties also calculate gas content in spherical overdensity apertures as well as their angular momentum. Must be used in conjunction with ``Extensive_halo_properties_output = 1``.
    ``Extensive_star_properties_output = 1``
        * Flag indicating that in addition to calculating extra halo properties also calculate stellar content in spherical overdensity apertures as well as their angular momentum. Must be used in conjunction with ``Extensive_halo_properties_output = 1``.
    Aperture related config options
        ``Calculate_aperture_quantities = 1``
            * Flag on whether to calculate aperture related masses, dispersions, metallicities
        ``Number_of_apertures = 6``
            * Number of spherical apertures
        ``Aperture_values_in_kpc = 3,5,10,30,50,100,``
            * Comma separated list of values in kpc
        ``Number_of_projected_apertures = 3``
            * Number of projected apertures. Code calculates 3 projections per aperture: x, y, z.
        ``Projected_aperture_values_in_kpc=10,50,100,``
            * Comma separated list of values in kpc
    Spherical overdensity related config options
        ``Number_of_overdensities = 5``
            * Number of spherical overdensities
        ``Overdensity_values_in_critical_density=25,100,500,1000,2500,``
            * Comma separated list of spherical overdensity thresholds in units of the critical density in cosmological simulations
    Radial profile related config options
        ``Calculate_radial_profiles = 1``
            * Flag on whether to calculate radial profiles of masses
        ``Radial_profile_norm = 0``
            * Flag setting the radial normalisation and scaling. Default is log rad bins, in proper kpc
        ``Number_of_radial_profile_bin_edges = 9``
            * Number of bin edges listed. Assumes lowest bin edge is r=0.
        ``Radial_profile_bin_edges = -2.,-1.50,-1.00,-0.50,0.00,0.50,1.00,1.50,2.00``
            * Comma separated list of (log) r bin edges. Here example is for log r in proper kpc binning so values are log(r).

.. topic:: Configuration for Extra Properties

    These are configuration options related to the bulk properties calculated based on extra
    properties of the particles. For instance, if hydro particles have a field called ``Turbulence``
    that contains some quantity of the internal turbulent energy and one wanted to calculate the
    average of this value for an object, one would use these options to load data from an HDF5 file
    (other inputs are not so easily parsed, making this not an option). One needs to provide what calculation
    to do (in the form of an integer flag specifying the calculation) and a string indicating the units.
    If the input is in the form of a 2D array from which a particular column is to be used, one can also
    set an index. The result is sorted in an output field that contains the name of the input field,
    the index (if >0), and a simple string describing the function along with the units and ending with particle type,
    ie: ``Turbulence_average_km/s^2_gas``
    These config options are combinations of particle type, categories and entry types.
    A full entry must be provided in a comma separated list and terminate in a comma.

    * Currently implemented are options for
        * ``Gas_``
        * ``Stars_``
        * ``BH_``
        *  ``Extra_DM_``
    * The currently catagories of properties are (except for ``Extra_DM`` which only has the first listed). An input field can be specified a number of times with different desired calculations to be run.
        * ``_internal_property``
        * ``_chemistry``
        * ``_chemistry_production``
    * The entries are
        * ``_names``
        * ``_index_in_file``
        * ``_calculation_type``
        * ``_input_output_unit_conversion_factors``
        * ``_output_units``

    Calculations allowed are as follows. You can add **massweighted** to any
    entry to calculate the mass weighted quantity. Note at entries should be lower case.
        * average
        * total
        * std (standard deviation)
        * min
        * max
        * logaverage (average(log(x)))
        * logstd (std(log(x)))

    One can also calculate total or average in aperutures provided aperture
    quantitites are being calculated.
        * aperture_total
        * aperture_average

    Example extra hydro Properties related config options
        ``Gas_internal_property_names = ,``
            * Names of fields to be read from an input HDF5 file that relate to hydro quantities, for which calculations can be done
        ``Gas_internal_property_index_in_file = ,``
            * Index in 2d array to be read from an input HDF5, useful for fields like metallicity where it is common to have an entry for each element
        ``Gas_internal_property_calculation_type = ,``
            * Integer flag indicating what calculation is to be done.
        ``Gas_internal_property_input_output_unit_conversion_factors = ,``
            * Float storing the conversion factor (if not 1.0) to take input units to output units.
        ``Gas_internal_property_output_units = ,``
            * String storing the units of the output.

.. _config_siminfo:

Simulation Info
^^^^^^^^^^^^^^^

Options related to the input and output units and cosmology.

.. _config_units:

.. topic:: Units

    Set internal (and output) units and conversion factors to well known units

    ``Length_unit =``
        * Factor by which input length unit is scaled, setting the internal code and output unit
    ``Velocity_unit =``
        * Factor by which input velocity unit is scaled, setting the internal code and output unit
    ``Mass_unit =``
        * Factor by which input mass unit is scaled, setting the internal code and output unit
    ``Gravity =``
        * Gravity in the internal output units, that is should be set such that :math:`v^2=Gm/r`, where v,m,r are the internal velocity, mass and length units.
    ``Hubble_unit =``
        * Unit of Hubble expansion in internal output units (from normal km/s/Mpc use 100). This is ignored if non-cosmological input
    ``Mass_value =``
        * If code is compiled not to store mass using the option **NOMASS** (see :ref:`compileoptions`) then set this value.
    ``Length_unit_to_kpc =``
        * Specify the conversion factor from the output unit to kpc
    ``Velocity_unit_to_kms =``
        * Specify the conversion factor from the output unit to km/s
    ``Mass_unit_to_solarmass =``
        * Specify the conversion factor from the output unit to solar masses
    ``Comoving_units = 1/0``
        * Flag indicating whether the properties output is in physical or comoving little h units.

.. _config_cosmology:

.. topic:: Cosmology

    If input is cosmological, then for some input formats (gadget, HDF), these quantites can be read from the input file. Tipsy formats require that these be set in the configuration file.

    ``Period = 0``
        * Period of the box in input units.
    ``Scale_factor = 1.0``
        * Scale factor time
    ``h_val = 1.0``
        * The "little h" value often used in cosmological simulations.
    ``Omega_m = 1.0``
        * Matter density in units of the critical density at z=0 used in cosmological simulations.
    ``Omega_Lambda = 0.0``
        * Energy density of the cosmological constant (or dark energy ) in units of the critical density at z=0 used in cosmological simulations.
    ``Omega_cdm = 1.0``
        * Dark matter density in units of the critical density at z=0 used in cosmological simulations. For non-standard DM models (annihilating, decaying, coupled), may be useful to provide the current DM density.
    ``Omega_b = 0.0``
        * Baryon density in units of the critical density at z=0 used in cosmological simulations.
    ``Omega_r = 0.0``
        * Radiation density in units of the critical density at z=0 used in cosmological simulations. Typically 0 (negligible).
    ``Omega_nu = 0.0``
        * Neutrino density in units of the critical density at z=0 used in cosmological simulations. Typically 0 (negligible).
    ``Omega_k = 0.0``
        * Curvature density in units of the critical density at z=0 used in cosmological simulations. Typically 0 (flat).
    ``Omega_DE = 0.0``
        * Dark Energy density in units of the critical density at z=0 used in cosmological simulations. This is addition to (or replacing) the energy density of the cosmological constant and has an associated equation of state, :math:`w_{DE}`.
    ``w_of_DE = -1.0``
        * Equation of state of the dark energy fluid, :math:`w=\frac{p}{\rho}`. This is not necessary unless one is using a cosmological simulation with :math:`w\neq -1`. Currently not fully implemented.
    ``Virial_density = 200.0``
        * Virial overdensity in units of the background matter density used in cosmological simulations. If -1, then the Bryan & Norman 1998 virial density is calculated based on a LCDM cosmology, otherwise overrides the Bryan & Norman calculation.
    ``Critical_density = 1.0``
        * Critical density in input units used in cosmological simulations.

.. _config_parallel:

Parallel
^^^^^^^^

Options related to MPI/OpenMP/Pthread parallelisation.

.. _config_mpi:

.. topic:: MPI

    MPI specific options

    ``MPI_part_allocation_fac = 0.1``
        * Factor used in memory allocated in mpi mode to store particles is (1+factor)* the memory need for the initial mpi decomposition. This factor should be >0 and is mean to allow a little room for particles to be exchanged between mpi threads withouth having to require new memory allocations and copying of data.
    ``MPI_particle_total_buf_size =``
        * Total memory size in bytes used to store particles in temporary buffer such that particles are sent to non-reading mpi processes in chunks of size buffer_size/NProcs/sizeof(Particle).
    ``MPI_number_of_tasks_per_write =``
        * Number of mpi tasks that are grouped for collective HDF5 writes is parallel HDF5 is enabled. Net result is that the total number of files written is ceiling(Number of MPI tasks)/(Number of tasks per write)

.. _config_openmp:

.. topic:: OpenMP specific parallelisation options

        ``OMP_run_fof = 1``
            * Flag indicating whether to run FOF searches with OpenMP threads.
        ``OMP_fof_region_size = 100000000``
            * Number of particles per OpenMP region.

.. _config_misc:

Miscellaneous
^^^^^^^^^^^^^

    Other configuration options

    ``Snapshot_value =``
        * If halo ids need to be offset to some starting value based on the snapshot of the output (say to make temporally unique halo ids that are useful for some halo merger tree codes), one can specific a snapshot number. All halo ids will be listed as internal haloid + snapnum * :math:`10^{12}` (or if using 32 bit integers and 64 bit integers, then ids offset by :math:`10^{6}`).
    ``Effective_Resolution =``
        * If running a multiple resolution zoom simulation, simple method of scaling the linking length by using the period and this effective resolution, ie: :math:`p/N_{\rm eff}`
    ``Verbose = 0/1/2``
        * Integer indicating how talkative the code is (2 very verbose, 1 verbose, 0 quiet).


.. _subsection_searchtypes:

Particle search types
^^^^^^^^^^^^^^^^^^^^^

List of particle types (and combinations) that can be searched are

.. doxygengroup:: SEARCHTYPES

.. _subsection_foftypes:

FOF search types
^^^^^^^^^^^^^^^^

List of fof aglorithms implemented are

.. doxygengroup:: FOFTYPES


.. _subsection_hdfnames:

HDF Input
^^^^^^^^^

List of naming conventions are

.. doxygengroup:: HDFNAMES

For complete discussion of implementation see :download:`../src/hdfitems.h`
