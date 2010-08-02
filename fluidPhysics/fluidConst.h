#ifndef FLUID_CONST_H
#define FLUID_CONST_H (1)

#include <Tools/tools.h>

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

namespace fluidPhysics {
  const double na  = 6.02252e23 ;           // Avogadro's number
  // US Standard Atmosphere Constant
  const double Rh  = 8314.3     ;           // Universal Gas Constant
  // Generally accepted value
  //  const double Rh  = 8314.472 ; // correct values
  const double EPSILON = 1e-30 ;
}

#endif
