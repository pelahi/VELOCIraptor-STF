.. _output:

Understanding and Analysing |vr| Output
###################################################

|vr| produces several different types of output files.

(with the mpi threads appending their rank to the end of the file name unless not compiled with MPI or if Parallel HDF5 is used.):

.. topic:: Standard files

    * ``.properties``: a file containing the bulk properties of all structures identified.
    * ``.catalog_groups``: a file containing the size of the structures (in number of particles associated) & information need to read particle information produced by velociraptor
    * ``.catalog_particles``: a file containing a list of particle IDs of those in structures. Information contained in ``.catalog_groups`` is used to parse this data.
    * ``.catalog_particles.unbound``: similar to ``catalog_particles`` but lists particles in structures but are formally unbound. Information contained in ``.catalog_groups`` is used to parse this data.

.. topic:: Extra files

    * ``.catalog_parttypes``: a file similar to ``.catalog_particles`` but containing a list of particle types of those in structures. Information contained in ``.catalog_groups`` is used to parse this data. Produced if multiple particle types are processed by |vr|.
    * ``.catalog_parttypes.unbound``: similar to ``catalog_parttypes`` but lists particles in structures but are formally unbound.
    * ``.profiles``  : a file containing the radial profiles of groups. Produced if radial profiles are requested.
    * ``.catalog_SOlist`` : a file containing the a list of particle IDs of particles found within a large Spherical region around Field halos. Produced if a list of paritcles wihtin so regions is requested.

Properties
==========

There are a variety of properties calculated for each object found. Some are
typical of all halo finders such as the mass of an object (which can be a halo,
subhalo, tidal debris), along with more complex properties such as the
eigenvectors and eigenvalues of the mass distribution defined by the reduced
inertia tensor. The number of properties also varies with the type of run. For
hydrodynamic simulations where |vr| has been compiled to use gas, star and black hole
properties, such as masses, temperatures, etc are also calculated. The code
will also calculate properties based on loading specific extra fields associated
with particle types (this interface requires HDF5 input or on the fly invocation
and outputs properties with the same name as the loaded property, see :ref:`usage`).

We give an almost complete list of properties and the keyword associate with the property (in ASCII and HDF5).
For clarity we list properties in several tables corresponding to

    - :ref:`standardhaloprops`,
    - :ref:`gasprops`,
    - :ref:`starprops`,
    - :ref:`bhprops`,
    - :ref:`interloperprops`,
    - :ref:`extradmprops`,

.. _standardhaloprops:

Standard Properties
-------------------

This is a list of standard properties that are alwas calculated unless indicated
otherwise (some require an extra configuration option).

.. _standardhalopropstable:

+-------------------+-------------------------------------------------------------------------------------------------------+
| Name              | Comments                                                                                              |
+===================+=======================================================================================================+
| **ID and Type information**                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| ID                | Halo ID. ID = index of halo + 1 + TEMPORALHALOIDVAL * Snapshot_value,                                 |
|                   | giving a temporally unique halo id that can be quickly parsed for an                                  |
|                   | index and a snapshot number.                                                                          |
+-------------------+-------------------------------------------------------------------------------------------------------+
| ID_mbp            | Particle ID of the most bound particle in the group.                                                  |
+-------------------+-------------------------------------------------------------------------------------------------------+
| hostHaloID        | ID of the host field halo. If an object is a field halo, this is -1.                                  |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Structuretype     | Structure types contain information on how the object was found and at                                |
|                   | what level in the subhalo hierarchy. Field halos are 10. Substructures                                |
|                   | identified using the local velocity field are type 10+10=20,                                          |
|                   | substructures identified using cores are type 10+5=15. For structures                                 |
|                   | found at level 2 (ie: subhalos within subhalos), the type offset is 20,                               |
|                   | and so on.                                                                                            |
+-------------------+-------------------------------------------------------------------------------------------------------+
| numSubStruct      | Number of substructures. Subhalos can have subsubhalos.                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| **Mass and radius properties**: `All properties are in output units.`                                                     |
+-------------------+-------------------------------------------------------------------------------------------------------+
| npart             | Number of particles belonging exclusively to the object.                                              |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Mass_tot          | Total mass of particles belonging exclusively to the object,                                          |
+-------------------+-------------------------------------------------------------------------------------------------------+
|                   |:math:`M_{\rm tot}`.                                                                                   |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Mass_FOF          | Total mass of particles in the FOF, :math:`M_{\rm FOF}`. Is zero for                                  |
|                   | substructure.                                                                                         |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Mass_200mean      | Overdensity mass defined by mean matter density, :math:`M_{200\rho_m}`.                               |
|                   | For field halos, if inclusive masses are desired, this is based on the                                |
|                   | particles in the FOF. If full spherical overdensity masses are desired,                               |
|                   | then includes all particles (whether they belong to the object, the                                   |
|                   | background or another object) within a spherical region. For subhalos,                                |
|                   | this is based on particles belonging exclusively to the object.                                       |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Mass_200crit      | Overdensity mass defined by critical density, :math:`M_{200\rho_c}`.                                  |
|                   | Behaviour like Mass_200mean.                                                                          |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Mass_BN98         | Overdensity mass defined by mean matter density and :math:`\Delta(z)`                                 |
|                   | given by Bryan & Norman (1998), :math:`M_{\Delta(z)\rho_c}`.                                          |
|                   | Behaviour like Mass_200mean.                                                                          |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Mvir              | User defined virial mass, :math:`M_{\rm vir}`. Behaviour like                                         |
|                   | Mass_200mean.                                                                                         |
+-------------------+-------------------------------------------------------------------------------------------------------+
| R_size            | Maximum distance of particles belonging exclusively to the object and                                 |
|                   | the centre-of-mass.                                                                                   |
+-------------------+-------------------------------------------------------------------------------------------------------+
| R_200mean         | Radius related to overdensity mass Mass_200mean.                                                      |
+-------------------+-------------------------------------------------------------------------------------------------------+
| R_200crit         | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| R_BN98            | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Rvir              | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| R_HalfMass        | Half mass radius based on the Mass_tot.                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| **Angular Momentum in Spherical Overdensity**: `Calculate if extra halo properties are requested`                         |
| `by setting the config option ` **Extensive_halo_properties_output=1**                                                    |
| `Also calculates inclusive spherical overdensity and also exclusive to halo as _exclusive.`                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Lx_200c           | :math:`x` component of the total angular momentum all the mass within :math:`R_{200\rho_c}`.          |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Ly_200c           | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Lz_200c           | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Lx_200m           | :math:`x` component of the total angular momentum all the mass within :math:`R_{200\rho_m}`.          |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Ly_200m           | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Lz_200m           | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Lx_BN98           | :math:`x` component of the total angular momentum all the mass within :math:`R_{BN98}`.               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Ly_BN98           | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Lz_BN98           | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| **Position and Velocity**: `All properties are in output units.`                                                          |
| `Objects have positions periodically wrapped.`                                                                            |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Xc                | :math:`x` coordinate of centre-of-mass.                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Yc                | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Zc                | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Xcmbp             | :math:`x` coordinate of most bound particle.                                                          |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Ycmbp             | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Zcmbp             | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Xcminpot          | :math:`x` coordinate of the minimum potential.                                                        |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Ycminpot          | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Zcminpot          | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| VXc               | :math:`v_x` velocity of centre-of-mass.                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| VYc               | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| VZc               | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| VXcmbp            | :math:`v_x` velocity of most bound particle.                                                          |
+-------------------+-------------------------------------------------------------------------------------------------------+
| VYcmbp            | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| VZcmbp            | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| VXcminpot         | :math:`v_x` velocity of the particle with the minimum potential.                                      |
+-------------------+-------------------------------------------------------------------------------------------------------+
| VYcminpot         | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| VZcminpot         | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| **Velocity and Angular Momentum**: `All properties are in output units.`                                                  |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Vmax              | Maximum circular velocity based on particles belonging exclusively to                                 |
|                   | the object, where circular velocities are defined by                                                  |
|                   | :math:`V_{\rm circ}^2=GM/R`.                                                                          |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Rmax              | Radius of maximum circular velocity.                                                                  |
+-------------------+-------------------------------------------------------------------------------------------------------+
| sigV              | Velocity dispersion based on the velocity dispersion tensor                                           |
|                   | :math:`\sigma_v=|\Sigma|^{1/6}`, where :math:`\Sigma` is the velocity                                 |
|                   | dispersion tensor.                                                                                    |
+-------------------+-------------------------------------------------------------------------------------------------------+
| veldisp_xx        | The :math:`x,x` component of the velocity dispersion tensor.                                          |
+-------------------+-------------------------------------------------------------------------------------------------------+
| veldisp_xy        | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| veldisp_xz        | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| veldisp_yx        | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| veldisp_yy        | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| veldisp_yz        | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| veldisp_zx        | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| veldisp_zy        | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| veldisp_zz        | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Lx                | :math:`x` component of the total angular momentum about the                                           |
|                   | centre-of-mass using particles belonging exclusively to the object.                                   |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Ly                | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Lz                | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| lambda_B          | Bullock et al (2001) like spin parameter :math:`\lambda_B` using total                                |
|                   | angular momentum and the spherical overdensity mass,                                                  |
|                   | :math:`\lambda_B=\frac{J}{\sqrt{2}MVR}`.                                                              |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Krot              | Measure of rotational support about the angular momentum axis                                         |
|                   | :math:`\kappa_{\rm rot}=\frac{\sum_i 1/2 m_i j_{z,i}r_i}{\sum_i T_i}`,                                |
|                   | where the first sum is over the motion of particles along the angular                                 |
|                   | momentum axis and the second sum is over kinetic energies                                             |
|                   | (see Sales et al (2010)).                                                                             |
+-------------------+-------------------------------------------------------------------------------------------------------+
| **Morphology**: `All properties are in output units.`                                                                     |
+-------------------+-------------------------------------------------------------------------------------------------------+
|                   | following Prada et al, (2012a) where we solve                                                         |
| cNFW              | Calculated assuming an NFW profile (Navarro, Frenk, & White 1997)                                     |
|                   | :math:`\frac{V_{\rm max}^2}{GM_\Delta/R_\Delta}-\frac{0.216c}{\ln(1+c)-c/(1+c)}=0.`                   |
+-------------------+-------------------------------------------------------------------------------------------------------+
| q                 | We calculate the shape using the reduced inertia tensor (Dubinski et al, 1991; Allgood et al, 2006),  |
|                   | :math:`\tilde{I}_{j,k}=\sum\limits_n \frac{m_n x^\prime_{j,n} x^\prime_{k,n}}{(r^\prime_{n})^2}`      |
|                   | where the sum is over particles exclusively belonging to the object                                   |
|                   | and, :math:`(r^\prime_n)^2=(x^\prime_n)^2+(y^\prime_n/q)^2+(z^\prime_n/s)^2`                          |
|                   | is the ellipsoidal distance between the halo's centre-of-mass and the                                 |
|                   | :math:`n_{\rm th}` particle, primed coordinates are in the eigenvector                                |
|                   | frame of the reduced inertia tensor and :math:`q` & :math:`s` are the                                 |
|                   | semi-major and minor axis ratios respectively. Thus :math:`q` is the                                  |
|                   | semi-major axis ratio. In eigenvector frame, x axis is major, y is semi-major, and z minor.           |
+-------------------+-------------------------------------------------------------------------------------------------------+
| s                 | Minor axis ratio.                                                                                     |
+-------------------+-------------------------------------------------------------------------------------------------------+
| eig_xx            | Eigenvectors of morphology.                                                                           |
+-------------------+-------------------------------------------------------------------------------------------------------+
| eig_xy            | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| eig_xz            | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| eig_yx            | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| eig_yy            | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| eig_yz            | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| eig_zx            | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| eig_zy            | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| eig_zz            | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| **Energy**: `All properties are in output units.`                                                                         |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Ekin              | The total kinetic energy, :math:`\sum T_i`.                                                           |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Epot              | The total gravitational potential energy :math:`1/2\sum W_i`, where  1/2 comes from double counting.  |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Efrac             | The fraction of particles that are formally bound (i.e., have :math:`W_i+T_i<0`).                     |
+-------------------+-------------------------------------------------------------------------------------------------------+
| **Quantities within** :math:`R(V_{\rm max})`: Properties based on particles within :math:`r\leq R(V_{\rm max})`.          |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_sigV        | Dispersion, like sigV for :math:`r\leq R(V_{\rm max})`.                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_veldisp_xx  | Dispersion tensor, like veldisp_xx for :math:`r\leq R(V_{\rm max})`.                                  |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_veldisp_xy  | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_veldisp_xz  | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_veldisp_yx  | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_veldisp_yy  | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_veldisp_yz  | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_veldisp_zx  | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_veldisp_zy  | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_veldisp_zz  | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_lambda_B    | Spin parameter, like lambda_B for :math:`r\leq R(V_{\rm max})`.                                       |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_Lx          | Total angular momentum, like Lx for :math:`r\leq R(V_{\rm max})`.                                     |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_Ly          | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_Lz          | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_q           | Semi-major axis ratio, like q for :math:`r\leq R(V_{\rm max})`.                                       |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_s           | Minor axisratio, like s for :math:`r\leq R(V_{\rm max})`.                                             |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_eig_xx      | Eigenvectors of morphology, like eig_xx for :math:`r\leq R(V_{\rm max})`.                             |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_eig_xy      | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_eig_xz      | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_eig_yx      | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_eig_yy      | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_eig_yz      | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_eig_zx      | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_eig_zy      | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| RVmax_eig_zz      | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
+-------------------+-------------------------------------------------------------------------------------------------------+
| **Additional Spherical Overdensity Mass/radius**: `If extra spherical overdensity values are requested via`               |
| Overdensity_values_in_critical_density `config option, code calculates masses/radii/angular momentum following`           |
| `a naming convention of` SO_property_rhocrivalue_rhocrit `where rhocritvalue is the overdensity value in units of the`    |
| `critical density, e.g.,` SO_mass_100_rhocrit.                                                                            |
| `The code will also calculate quantities based on particle type: gas, star, interloper, following`                        |
| SO_property_parttype_rhocrivalue_rhocrit                                                                                  |
+-------------------+-------------------------------------------------------------------------------------------------------+
| mass              | Mass enclosing a average density of the associated SO value.                                          |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Lx                | Angular momentum of enclosed mass in x-direction                                                      |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Ly                | |ditto| in y-direction                                                                                |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Lz                | |ditto| in z-direction                                                                                |
+-------------------+-------------------------------------------------------------------------------------------------------+
| **Aperture quantities**: `If aperture quantities are requested code calculates a variety of properties`                   |
| `within spherical aperture in pkpc.`                                                                                      |
| `Naming convention is` Aperture_quantity_radiusvalue_kpc, `or for a specific` `particle type`                             |
| Aperture_quantity_parttype_radiusvalue_kpc, `e.g.` Aperture_mass_10_kpc.                                                  |
| `Particle types where individual quantities are calculated: gas, star, bh, interloper.`                                   |
| `We list the property names here.`                                                                                        |
+-------------------+-------------------------------------------------------------------------------------------------------+
| mass              | Total mass in aperture.                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| npart             | Total number of particles.                                                                            |
+-------------------+-------------------------------------------------------------------------------------------------------+
| rhalfmass         | Radius enclosing half the mass within the aperture.                                                   |
+-------------------+-------------------------------------------------------------------------------------------------------+
| veldisp           | Velocity disperion                                                                                    |
+-------------------+-------------------------------------------------------------------------------------------------------+
| **Projected aperture quantities**: `Similar to aperture quantitites but for 3 different projections based on particles`   |
| `within a projected radius in pkpc. Naming convention is` Projected_aperture_i_quantity_radiusvalue_kpc, `where`          |
| `i is from 0, 1, 2 for a x,y,z projection.`                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| mass              | Total mass in aperture.                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| rhalfmass         | Radius enclosing half the mass within the aperture.                                                   |
+-------------------+-------------------------------------------------------------------------------------------------------+

.. _gasprops:

Gas Properties
--------------

This is a list of gas properties that are calculated if code is compiled with
**USE_GAS**. Some require an extra configuration option. Also, Spherical overdensity
masses + angular momentum, aperture properties, projected aperture properties are calculated
for gas particles along along with some extra gas only properties.

.. _gaspropstable:

+-------------------+-------------------------------------------------------------------------------------------------------+
| Name              | Comments                                                                                              |
+===================+=======================================================================================================+
| **Gas quantities**: `Bulk properties of gas particles/tracers when compiled to process gas properties. Properties unique` |
| `to gas are T_gas and SFR_gas.`                                                                                           |
+-------------------+-------------------------------------------------------------------------------------------------------+
| n_gas             | Number of gas particles.                                                                              |
+-------------------+-------------------------------------------------------------------------------------------------------+
| M_gas             | Total gas mass :math:`M_{\rm gas}`.                                                                   |
+-------------------+-------------------------------------------------------------------------------------------------------+
| M_gas_Rvmax       | Gas mass within :math:`R(V_{\rm max})`.                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| M_gas_30kpc       | Gas mass within 30 pkpc.                                                                              |
+-------------------+-------------------------------------------------------------------------------------------------------+
| M_gas_500c        | Gas mass within a spherical overdensity of :math:`500\rho_c`.                                         |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Xc_gas            | :math:`x` coordinate of centre-of-mass of gas particles relative to Xc.                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Yc_gas            | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Zc_gas            | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| VXc_gas           | :math:`x` coordinate of centre-of-mass velocity of gas particles relative to VXc.                     |
+-------------------+-------------------------------------------------------------------------------------------------------+
| VYc_gas           | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| VZc_gas           | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Efrac_gas         | Like Efrac but for gas particles only.                                                                |
+-------------------+-------------------------------------------------------------------------------------------------------+
| R_HalfMass_gas    | Like R_HalfMass but for gas particles only.                                                           |
+-------------------+-------------------------------------------------------------------------------------------------------+
| veldisp_xx_gas    | Like veldisp_xx but for gas particles only and relative to the centre-of-mass.                        |
+-------------------+-------------------------------------------------------------------------------------------------------+
| veldisp_xy_gas    | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| veldisp_xz_gas    | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| veldisp_yx_gas    | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| veldisp_yy_gas    | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| veldisp_yz_gas    | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| veldisp_zx_gas    | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| veldisp_zy_gas    | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| veldisp_zz_gas    | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Lx_gas            | Like Lx but for gas particles only and relative to the centre-of-mass.                                |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Ly_gas            | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Lz_gas            | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| q_gas             | Like q but for gas particles only and relative to the centre-of-mass.                                 |
+-------------------+-------------------------------------------------------------------------------------------------------+
| s_gas             | Like s but for gas particles only and relative to the centre-of-mass.                                 |
+-------------------+-------------------------------------------------------------------------------------------------------+
| eig_xx_gas        | Like eig_xx but for gas particles only and relative to the centre-of-mass.                            |
+-------------------+-------------------------------------------------------------------------------------------------------+
| eig_xy_gas        | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| eig_xz_gas        | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| eig_yx_gas        | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| eig_yy_gas        | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| eig_yz_gas        | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| eig_zx_gas        | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| eig_zy_gas        | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| eig_zz_gas        | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Krot_gas          | Like Krot but for gas particles only and relative to the centre-of-mass.                              |
+-------------------+-------------------------------------------------------------------------------------------------------+
| T_gas             | Average temperature of gas.                                                                           |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Zmet_gas          | Average metallicity of gas.                                                                           |
+-------------------+-------------------------------------------------------------------------------------------------------+
| SFR_gas           | Total star formation rate of gas.                                                                     |
+-------------------+-------------------------------------------------------------------------------------------------------+
| **Star Forming (sf)/Non Star Forming (nsf) Gas quantities**: `Similar to gas properties but split by sf/nsf gas.`         |
| `For brevity, we list only quantities unique to sf, as the nsf gas is similar but with _nsf naming convention.`           |
| `Only calculated if` **USE_GAS** `and` **USE_STAR** `flags on.`                                                           |
+-------------------+-------------------------------------------------------------------------------------------------------+
| M_gas_sf          | Total gas mass :math:`M_{\rm gas}`.                                                                   |
+-------------------+-------------------------------------------------------------------------------------------------------+
| R_HalfMass_gas_sf | Half mass radii.                                                                                      |
+-------------------+-------------------------------------------------------------------------------------------------------+
| sigV_gas_sf       | Velocity dispersion of the gas.                                                                       |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Lx_gas_sf         | Like Lx_gas but for star forming gas.                                                                 |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Ly_gas_sf         | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Lz_gas_sf         | |ditto|                                                                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Krot_gas_sf       | Like Krot_gas but for star forming gas                                                                |
+-------------------+-------------------------------------------------------------------------------------------------------+
| T_gas_sf          | Average temperature of star forming gas.                                                              |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Zmet_gas_sf       | Average metallicity of star forming gas.                                                              |
+-------------------+-------------------------------------------------------------------------------------------------------+
| **Aperture quantities**: `If aperture quantities are requested code calculates a variety of properties`                   |
| `within spherical aperture in pkpc.`                                                                                      |
| `Naming convention is` Aperture_quantity_gas_radiusvalue_kpc.                                                             |
| `We list the additional properties calculated for gas here (which are in addition to mass,rhalfmass, etc).`               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Zmet              | Average gas metallicity in aperture.                                                                  |
+-------------------+-------------------------------------------------------------------------------------------------------+
| SFR               | Total star formation rate of gas in aperture.                                                         |
+-------------------+-------------------------------------------------------------------------------------------------------+
| **Projected aperture quantities**: `Similar to aperture quantitites but for 3 different projections based on particles`   |
| `within a projected radius in pkpc. Naming convention is` Projected_aperture_i_quantity_gas_radiusvalue_kpc, `where`      |
| `i is from 0, 1, 2 for a x,y,z projection.`                                                                               |
| `We list the additional properties calculated for gas here (which are in addition to mass,rhalfmass, etc).`               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Zmet              | Average gas metallicity in projected aperture.                                                        |
+-------------------+-------------------------------------------------------------------------------------------------------+
| SFR               | Total star formation rate of gas in projected aperture.                                               |
+-------------------+-------------------------------------------------------------------------------------------------------+

+--------------------------------------------------------+------------------------------------------------------------------+
| Name                                                   | Comments                                                         |
+========================================================+==================================================================+
| **Extra Gas Properties**: `If extra gas fields are loaded by listing them using` Gas_internal_property_names              |
| Gas_chemistry_names `and/or` Gas_chemistry_production_names. `The are associated input options related to the input index |
| calclation type done and output units. The output will have the following naming convention:`                             |
| nameoffield_index_#_calculation_units_gas `e.g.``, AlphaElements_index_0_average_unitless_gas.                            |
| `Also requires that code is compiled with the` **USE_GAS** `flag`                                                         |
| `As an example we show the fields if`                                                                                     |
| Gas_internal_property_names=Pressure,MetalMassFractionFromSNIa,                                                           |
| Gas_internal_property_index=0,1,                                                                                          |
| Gas_internal_property_output_units=kPa,unitless,                                                                          |
| Gas_internal_property_calculation_type=max,average,                                                                       |
+--------------------------------------------------------+------------------------------------------------------------------+
| Pressure_index_0_max_kPa_gas                           | maximum pressure of gas in object.                               |
+--------------------------------------------------------+------------------------------------------------------------------+
| MetalMassFractionFromSNIa_index_1_average_unitless_gas | average of this field.                                           |
+--------------------------------------------------------+------------------------------------------------------------------+
| `One can also specify` aperture_total `and` aperture_average `as functions if aperture quantities are calcualed.          |
| The output will have a simlar naming convention to above but with` Aperture_ `at the start and ending with the aperture   |
| aperture itself` #_kpc` for each aperture listed.                                                                         |
+--------------------------------------------------------+------------------------------------------------------------------+

.. _starprops:

Star Properties
---------------

This is a list of stellar properties that are calculated if code is compiled with
**USE_STAR**. Some require an extra configuration option.

.. _starpropstable:

+-------------------+-------------------------------------------------------------------------------------------------------+
| Name              | Comments                                                                                              |
+===================+=======================================================================================================+
| **Star quantities**: `Bulk stellar properties when compiled to process star properties. Similar to gas properties`        |
| `but has _star instead of _ gas. For brevity, we list only quantities unique to star particles.`                          |
+-------------------+-------------------------------------------------------------------------------------------------------+
| tage_star          | Average stellar age.                                                                                 |
+-------------------+-------------------------------------------------------------------------------------------------------+
| **Aperture quantities**: `If aperture quantities are requested code calculates a variety of properties`                   |
| `within spherical aperture in pkpc.`                                                                                      |
| `Naming convention is` Aperture_quantity_star_radiusvalue_kpc.                                                            |
| `We list the additional properties calculated for star here (which are in addition to mass,rhalfmass, etc).`              |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Zmet              | Average stellar metallicity in aperture.                                                              |
+-------------------+-------------------------------------------------------------------------------------------------------+
| **Projected aperture quantities**: `Similar to aperture quantitites but for 3 different projections based on particles`   |
| `within a projected radius in pkpc. Naming convention is` Projected_aperture_i_quantity_star_radiusvalue_kpc, `where`     |
| `i is from 0, 1, 2 for a x,y,z projection.`                                                                               |
| `We list the additional properties calculated for gas here (which are in addition to mass,rhalfmass, etc).`               |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Zmet              | Average stellar metallicity in projected aperture.                                                    |
+-------------------+-------------------------------------------------------------------------------------------------------+
| **Extra Star Properties**: `Like the extra gas properties but calculated if ` Star_internal_property_names                |
| Star_chemistry_names `and/or` Star_chemistry_production_names.                                                            |
| `Naming convention is the same but ends with _star`                                                                       |
| `Also requires that code is compiled with the` **USE_STAR** `flag`                                                        |
+-------------------+-------------------------------------------------------------------------------------------------------+

.. _bhprops:

Black Hole Properties
---------------------

This is a list of black hole properties that are calculated if code is compiled with
**USE_BH**. Some require an extra configuration option.

.. _bhpropstable:

+-------------------+-------------------------------------------------------------------------------------------------------+
| Name              | Comments                                                                                              |
+===================+=======================================================================================================+
| **Black hole quantities**: `Bulk properties of black hole particles when compiled to process black hole properties.`      |
+-------------------+-------------------------------------------------------------------------------------------------------+
| n_bh              | Number of black hole particles.                                                                       |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Mass_bh           | Total mass of black hole particles.                                                                   |
+-------------------+-------------------------------------------------------------------------------------------------------+
| **Extra Black hole Properties**: `Like the extra gas properties but calculated if ` BH_internal_property_names            |
| BH_chemistry_names `and/or` BH_chemistry_production_names.                                                                |
| `Naming convention is simialr save ends with _bh`                                                                         |
| `Also requires that code is compiled with the` **USE_BH** `flag`                                                          |
+-------------------+-------------------------------------------------------------------------------------------------------+

.. _interloperprops:

Interloper Properties
---------------------

This is a list of interloper DM properties that are calculated if code is compiled with
**ZOOM_SIM**. These properties are based on low resolution particles and can be
used to gauge the level of contamination

.. _interloperpropstable:

+-------------------+-------------------------------------------------------------------------------------------------------+
| Name              | Comments                                                                                              |
+===================+=======================================================================================================+
| **Interloper particles**: `If analysing multi-resolution simulations, low resolution particles are often treated as`      |
| `contaminants. These are bulk properties of low resolution contaminant particles.`                                        |
+-------------------+-------------------------------------------------------------------------------------------------------+
| n_interloper      | Number of low resolution, interloper particles.                                                       |
+-------------------+-------------------------------------------------------------------------------------------------------+
| Mass_interloper   | Total mass of low resolution, interloper particles.                                                   |
+-------------------+-------------------------------------------------------------------------------------------------------+


.. _extradmprops:

Extra DM Properties
-------------------

This is a list of Extra DM properties that are calculated if code is compiled with
**USE_EXTRADM**. These properties are useful if running on standard dark matter.

.. _extradmpropstable:

+-------------------+-------------------------------------------------------------------------------------------------------+
| Name              | Comments                                                                                              |
+===================+=======================================================================================================+
| **Extra DM Properties**: `Like the extra gas properties but calculated if ` Extra_DM_internal_property_names              |
| `has fields specified. Useful for nonstandard dark matter runs, such as annihilating or interacting dark matter.`         |
| `Naming convention is nameoffield_extra_dm`                                                                               |
| `Also requires that code is compiled with the` **USE_EXTRADM** `flag`                                                     |
+-------------------+-------------------------------------------------------------------------------------------------------+

.. |ditto| unicode:: U+03003 .. ditto mark
