#ifndef CONST_H
#define CONST_H (1)

#include "sciTypes.h"

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

namespace streamUns {
  const real na  = 6.02252e23 ;           // Avogadro's number
  const real Rh  = 8314.3     ;           // Universal Gas Constant 
  const real pi  = M_PI ;
  
  const real EPSILON = 1e-30 ;
}

#endif
