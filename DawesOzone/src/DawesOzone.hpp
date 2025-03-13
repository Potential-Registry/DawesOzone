#ifndef RYNLIB_DAWESOZONE_POT_HPP

#include "Python.h"
#include "RynTypes.hpp"


        Real_t DawesOzone_calc_potential_(
            const FlatCoordinates,
            const Names,
            const ExtraBools,
            const ExtraInts,
            const ExtraFloats
            );
        extern "C" { void calc_potential_(RawWalkerBuffer,int*,double*,double*);} 

#define RYNLIB_DAWESOZONE_POT_HPP

#endif //RYNLIB_DAWESOZONE_POT_HPP