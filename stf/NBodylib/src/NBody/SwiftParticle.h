/*! \file SwiftParticle.h
 *  \brief header file for the SWIFT particle type.

*/

#ifndef SWIFT_PARTICLE_H
#define SWIFT_PARTICLE_H

#ifdef SWIFTINTERFACE
  /**
 * @brief The default struct alignment in SWIFT.
 */
#define SWIFT_STRUCT_ALIGNMENT 32

namespace Swift
{
    /* Include some struct definitions from SWIFT. */
    extern "C" {
        #include "align.h"
        #include "timeline.h"
        #include "part_type.h"
        #include "gravity_part.h"
        #include "hydro_part.h"
    }

    /* SWIFT/VELOCIraptor particle. */
    struct swift_vel_part {

      /*! Particle ID. If negative, it is the negative offset of the #part with
        which this gpart is linked. */
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

    } SWIFT_STRUCT_ALIGN;

}

#endif

#endif
