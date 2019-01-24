.. _output:

Understanding and Analysing **VELOCIraptor** Output
###################################################

**VELOCIraptor** produces produce several different types of output files.

(with the mpi threads appending their rank to the end of the file name):

.. topic:: Standard files

    * ``.properties``: a file containing the bulk properties of all structures identified.
    * ``.catalog_groups``: a file containing the size of the structures (in number of particles associated) & information need to read particle information produced by velociraptor
    * ``.catalog_particles``: a file containing a list of particle IDs of those in structures. Information contained in ``.catalog_groups`` is used to parse this data.
    * ``.catalog_particles.unbound``: similar to ``catalog_particles`` but lists particles in structures but are formally unbound. Information contained in ``.catalog_groups`` is used to parse this data.


Properties
==========

There are a variety of properties calculated for each object found. Some are typical of all halo finders
such as the mass of an object (which can be a halo, subhalo, tidal debris), along with more complex properties
such as the eigenvectors and eigenvalues of the mass distribution defined by the reduced inertia tensor.
The number of properties also varies with the type of run. For hydrodynamic simulations where **VELOCIraptor**
has been compiled to use gas properties and star properties, gas masses, temperatures, etc are also calculated.

We give an almost complete list of properties and the keyword associate with the property (in ASCII, HDF5 and ADIOS outputs)
.. table:: properties_table

+===================+========================================================================+
|Name               |Comments                                                                |
:------------------:|:-----------------------------------------------------------------------:
|Center             |Left                                                                    |
+===================+========================================================================+
|ID and Type information                                                                     |
+===================+========================================================================+
| ID                | Halo ID. ID = index of halo + 1 + TEMPORALHALOIDVAL * Snapshot_value,  |
|                   | giving a temporally unique halo id that can be quickly parsed for an   |
|                   | index and a snapshot number.                                           |
+-------------------+------------------------------------------------------------------------+
| ID_mbp             |Particle ID of the most bound particle in the group.                   |
+-------------------+------------------------------------------------------------------------+
| hostHal oID       | ID of the host field halo. If an object is a field halo, this is -1.   |
| Structuretype     | Structure types contain information on how the object was found and at |
|                   | what level in the subhalo hierarchy. Field halos are 10. Substructures |
|                   | identified using the local velocity field are type 10+10=20,           |
|                   | substructures identified using cores are type 10+5=15. For structures  |
|                   | found at level 2 (ie: subhalos within subhalos), the type offset is 20,|
|                   | and so on.                                                             |
+-------------------+------------------------------------------------------------------------+
| numSubStruct      | Number of substructures. Subhalos can have subsubhalos.                |
+===================+========================================================================+
|Mass and radius properties. All properties are in output units.                             |
+===================+========================================================================+
| npart             | Number of particles belonging exclusively to the object.               |
+-------------------+------------------------------------------------------------------------+
| Mass_tot          | Total mass of particles belonging exclusively to the object,           |
+-------------------+------------------------------------------------------------------------+
|                   |:math:`M_{\rm tot}`.                                                    |
+-------------------+------------------------------------------------------------------------+
| Mass_FOF          | Total mass of particles in the FOF, :math:`M_{\rm FOF}`. Is zero for   |
|                   | substructure.                                                          |
+-------------------+------------------------------------------------------------------------+
| Mass_200mean      | Overdensity mass defined by mean matter density, :math:`M_{200\rho_m}`.|
|                   | For field halos, if inclusive masses are desired, this is based on the |
|                   | particles in the FOF. If full spherical overdensity masses are desired,|
|                   | then includes all particles (whether they belong to the object, the    |
|                   | background or another object) within a spherical region. For subhalos, |
|                   | this is based on particles belonging exclusively to the object.        |
+-------------------+------------------------------------------------------------------------+
| Mass_200crit      | Overdensity mass defined by critical density, :math:`M_{200\rho_c}`.   |
|                   | Behaviour like Mass_200mean.                                           |
+-------------------+------------------------------------------------------------------------+
| Mass_BN98         | Overdensity mass defined by mean matter density and :math:`\Delta(z)`  |
|                   | given by Bryan & Norman (1998), :math:`M_{\Delta(z)\rho_c}`.           |
|                   | Behaviour like Mass_200mean.                                           |
+-------------------+------------------------------------------------------------------------+
| Mvir              | User defined virial mass, :math:`M_{\rm vir}`. Behaviour like Mass_200mean.|
+-------------------+------------------------------------------------------------------------+
| R_size            | Maximum distance of particles belonging exclusively to the object and  |
|                   | the centre-of-mass.                                                    |
+-------------------+------------------------------------------------------------------------+
| R_200mean         | Radius related to overdensity mass Mass_200mean.                       |
+-------------------+------------------------------------------------------------------------+
| R_200crit         | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| R_BN98            | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| Rvir              | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| R_HalfMass        | Half mass radius based on the Mass_tot.                                |
+===================+========================================================================+
|Position and Velocity : All properties are in output units.                                 |
|Objects need not have positions periodically wrapped.                                       |
+===================+========================================================================+
| Xc                | :math:`x` coordinate of centre-of-mass.                                |
+-------------------+------------------------------------------------------------------------+
| Yc                | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| Zc                | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| Xcmbp             | $x$ coordinate of most bound particle.                                 |
+-------------------+------------------------------------------------------------------------+
| Ycmbp             | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| Zcmbp             | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| VXc               | $v_x$ velocity of centre-of-mass.                                      |
+-------------------+------------------------------------------------------------------------+
| VYc               | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| VZc               | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| VXcmbp            | $v_x$ velocity of most bound particle.                                 |
+-------------------+------------------------------------------------------------------------+
| VYcmbp            | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| VZcmbp            | :math:`\ditto`                                                         |
+===================+========================================================================+
Velocity and Angular Momentum : All properties are in output units.                          |
+===================+========================================================================+
| Vmax              | Maximum circular velocity based on particles belonging exclusively to  |
|                   | the object, where circular velocities are defined by                   |
|                   | :math:`V_{\rm circ}^2=GM/R`.                                           |
+-------------------+------------------------------------------------------------------------+
| Rmax              | Radius of maximum circular velocity.                                   |
+-------------------+------------------------------------------------------------------------+
| sigV              | Velocity dispersion based on the velocity dispersion tensor            |
|                   | :math:`\sigma_v=|\Sigma|^{1/6}`, where :math:`\Sigma` is the velocity  |
|                   | dispersion tensor.                                                     |
+-------------------+------------------------------------------------------------------------+
| veldisp_xx        | The :math:`x,x` component of the velocity dispersion tensor.           |
+-------------------+------------------------------------------------------------------------+
| veldisp_xy        | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| veldisp_xz        | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| veldisp_yx        | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| veldisp_yy        | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| veldisp_yz        | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| veldisp_zx        | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| veldisp_zy        | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| veldisp_zz        | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| Lx                | :math:`x` component of the total angular momentum about the            |
|                   | centre-of-mass using particles belonging exclusively to the object.    |
+-------------------+------------------------------------------------------------------------+
| Ly                | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| Lz                | :math:`\ditto`                                                         |
+-------------------+------------------------------------------------------------------------+
| lambda_B          | Bullock et al (2001) like spin parameter :math:`\lambda_B` using total |
|                   | angular momentum and the spherical overdensity mass,                   |
|                   | :math:`\lambda_B=\frac{J}{\sqrt{2}MVR}`.                               |
+-------------------+------------------------------------------------------------------------+
| Krot              | Measure of rotational support about the angular momentum axis          |
|                   | :math:`\kappa_{\rm rot}=\frac{\sum_i 1/2 m_i j_{z,i}r_i}\sum_i T_i`,   |
|                   | where the first sum is over the motion of particles along the angular  |
|                   | momentum axis and the second sum is over kinetic energies              |
|                   | (see Sales et al (2010)).                                              |
+===================+========================================================================+
| Morphology : All properties are in output units.
+===================+========================================================================+
