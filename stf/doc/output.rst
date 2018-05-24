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
