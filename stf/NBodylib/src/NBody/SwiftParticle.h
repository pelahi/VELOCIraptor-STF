/*! \file SwiftParticle.h
 *  \brief header file for the SWIFT particle type.

*/

#ifndef SWIFT_PARTICLE_H
#define SWIFT_PARTICLE_H

#ifdef SWIFTINTERFACE

namespace Swift
{
  
  /* SWIFT enum of part types. Should match VELOCIraptor type values. */
  enum part_type {
    swift_type_gas = 0,
    swift_type_dark_matter = 1,
    swift_type_star = 4,
    swift_type_black_hole = 5,
    swift_type_count
  } __attribute__((packed));

  /* SWIFT/VELOCIraptor particle. */
    struct swift_vel_part {

      /*! Particle ID. */
      long long id;

      /*! Particle position. */
      double x[3];

      /*! Particle velocity. */
      float v[3];

      /*! Particle mass. */
      float mass;

      /*! Gravitational potential */
      float potential;

      /*! Internal energy of gas particle */
      float u;

      /*! Type of the #gpart (DM, gas, star, ...) */
      enum part_type type;

    };

}

#endif

#endif
