/*! \file SwiftParticle.h
 *  \brief header file for the SWIFT particle type.

*/

#ifndef SWIFT_PARTICLE_H
#define SWIFT_PARTICLE_H

#define SWIFT_STRUCT_ALIGNMENT 32 
#define SWIFT_STRUCT_ALIGN __attribute__((aligned(SWIFT_STRUCT_ALIGNMENT)))

namespace Swift
{

    /* The different types of particles a #gpart can link to.
       Note we use the historical values from Gadget for these fields. */
    enum part_type {
        swift_type_gas = 0,
        swift_type_dark_matter = 1,
        swift_type_star = 4,
        swift_type_black_hole = 5,
        swift_type_count
    };

    /* Gravity particle. */
    struct gpart {

        /* Particle ID. If negative, it is the negative offset of the #part with
           which this gpart is linked. */
        long long id_or_neg_offset;

        /* Particle position. */
        double x[3];

        /* Offset between current position and position at last tree rebuild. */
        float x_diff[3];

        /* Particle velocity. */
        float v_full[3];

        /* Particle acceleration. */
        float a_grav[3];

        /* Particle mass. */
        float mass;

        /* Gravitational potential */
        float potential;

        /* Softening length */
        float epsilon;

        /* Time-step length */
        //timebin_t time_bin;

        /* Type of the #gpart (DM, gas, star, ...) */
        enum part_type type;

        //#ifdef SWIFT_DEBUG_CHECKS
        //
        //    /* Numer of gparts this gpart interacted with */
        //    long long num_interacted;
        //
        //    /* Time of the last drift */
        //    integertime_t ti_drift;
        //
        //    /* Time of the last kick */
        //    integertime_t ti_kick;
        //
        //#endif
        //
        //#ifdef SWIFT_GRAVITY_FORCE_CHECKS
        //
        //    /* Brute-force particle acceleration. */
        //    double a_grav_exact[3];
        //
        //#endif

    } SWIFT_STRUCT_ALIGN;

}

#endif
