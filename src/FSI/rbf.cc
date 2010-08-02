#include <sciTypes.h>
#include <cmath>
#include "rbf.h"

namespace streamUns {
	
// RBF: 1-8 compact, 9-14 global
real radialBasisFunction(const real x, const real r, const real a, const int rbfNr) {
	real zeta = x / r ;
	real phi = 0. ;
	real zeta2, zeta3, zeta4, zeta5 ;
	zeta2 = zeta*zeta ;

	if ((zeta >= 1.0) && ( rbfNr > 0 && rbfNr < 9)) {
		return 0. ;
	}

	switch (rbfNr) {
		case 1 : // CP C0
			phi = pow((1. - zeta), 2.) ;
			break ;
		case 2 : // CP C2
			phi = pow((1. - zeta), 4.) * (4. * zeta + 1) ;
			break ;
		case 3 : // CP C4
			phi = pow((1. - zeta), 6.) * (35./3. * zeta2 + 6. * zeta + 1.) ;
			break ;
		case 4 : // CP C6
			zeta3 = zeta2 * zeta ;
			phi = pow((1. - zeta), 8.) * (32. * zeta3 + 25.*zeta2 + 8. * zeta + 1.) ;
			break ;
		case 5 : // CTPS C0
			phi = pow((1. - zeta), 5.) ;
			break ;
		case 6 : // CTPS C1
			zeta3 = zeta2 * zeta ;
			zeta4 = zeta3 * zeta ;
			zeta5 = zeta4 * zeta ;
			phi = 1. + 80./3. * zeta2 - 40. * zeta3 + 15.*zeta4 - 8./3. * zeta5 + 20. * zeta2 * log(zeta) ;
			break ;
		case 7 : // CTPS C2a
			zeta3 = zeta2 * zeta ;
			zeta4 = zeta3 * zeta ;
			zeta5 = zeta4 * zeta ;
			phi = 1. - 30. * zeta2 - 10. * zeta3 + 45.*zeta4 - 6. * zeta5 - 60. * zeta3 * log(zeta) ;
			break ;
		case 8 : // CTPS C2b
			zeta3 = zeta2 * zeta ;
			zeta4 = zeta3 * zeta ;
			zeta5 = zeta4 * zeta ;
			phi = 1. - 20. * zeta2 + 80. * zeta3 - 45.*zeta4 - 16. * zeta5 + 60. * zeta4 * log(zeta) ;
			break ;
/*		case 9 : // TPS
			if (x > 0.) {
				phi = x*x * log(x) ;
			} else {
				phi = 0.;
			}
			break ;
		case 10 : // MQB
			phi = sqrt(a*a + x*x) ;
			break ;
		case 11 : // IMQB
			phi = sqrt(1. / (a*a + x*x)) ;
			break ;
		case 12 : // QB
			phi = 1. + x*x ;
			break ;
		case 13 : // IQB
			phi = 1. / (1. + x*x) ;
			break ;
		case 14 : // Gauss
			phi = exp(-x*x) ;
			break ;
*/
		case 9 : // TPS
			if (zeta > 0.) {
				phi = zeta2 * log(zeta) ;
			} else {
				phi = 0.;
			}
			break ;
		case 10 : // MQB
			phi = sqrt(a*a + zeta2) ;
			break ;
		case 11 : // IMQB
			phi = sqrt(1. / (a*a + zeta2)) ;
			break ;
		case 12 : // QB
			phi = 1. + zeta2 ;
			break ;
		case 13 : // IQB
			phi = 1. / (1. + zeta2) ;
			break ;
		case 14 : // Gauss
			phi = exp(-zeta2) ;
			break ;
/*		case 9 : // TPS
			if (x > 0.) {
				phi = x*x * log(x) ;
			} else {
				phi = 0.;
			}
			break ;
		case 10 : // MQB
			phi = sqrt(a*a + x*x) ;
			break ;
		case 11 : // IMQB
			phi = sqrt(1. / (a*a + x*x)) ;
			break ;
		case 12 : // QB
			phi = 1. + x*x ;
			break ;
		case 13 : // IQB
			phi = 1. / (1. + x*x) ;
			break ;
		case 14 : // Gauss
			phi = exp(-x*x) ;
			break ;
			*/
	}
//	if (phi < 1.0e-50) {
//		cout << "RBF phi too small! -> 0" << endl ;
//		phi = 0.0 ;
//	}
	return phi ; 
} ;

}
