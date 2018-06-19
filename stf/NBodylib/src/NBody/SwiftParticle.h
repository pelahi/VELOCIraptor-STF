/*! \file SwiftParticle.h
 *  \brief header file for the SWIFT particle type.

*/

#ifndef SWIFT_PARTICLE_H
#define SWIFT_PARTICLE_H

#ifdef SWIFTINTERFACE

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
}

#endif

#endif
