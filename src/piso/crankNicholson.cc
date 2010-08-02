//-----------------------------------------------------------------------------
// Description: This file contains rules for the fully-implicit
//   Crank-Nicholson scheme for all the governing equations.
//-----------------------------------------------------------------------------

// Standard library includes.
#include <vector>
using std::vector ;

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// CHEM includes.
#include "eos.h"

// StreamUns includes.
#include "const.h"
#include "referenceFrame.h"
#include "residual.h"
#include "sciTypes.h"
#include "varsFileInputs.h"

namespace streamUns {

//-----------------------------------------------------------------------------
// General rules.

  // Creates the time integrator constraints.
  class ThetaParameter : public singleton_rule {
    private:
      const_param<string> timeIntegrator ;
      param<real> thetaParameter ;
    public:

      // Define input and output.
      ThetaParameter() {
        name_store("timeIntegrator",timeIntegrator) ;
        name_store("thetaParameter",thetaParameter) ;
        input("timeIntegrator") ;
        output("thetaParameter") ;
      }

      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        *thetaParameter=1.0 ; if(*timeIntegrator=="CN") *thetaParameter=0.5 ;
      }
  } ;

  register_rule<ThetaParameter> registerThetaParameter ;

}
