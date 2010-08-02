#ifndef RBF_H
#define RBF_H

#include <sciTypes.h>
#include <cmath>

typedef double real ;

namespace streamUns {
	
// RBF: 1-8 compact, 9-14 global
real radialBasisFunction(const real x, const real r, const real a, const int rbfNr) ;

}
#endif