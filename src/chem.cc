// This file contains miscellaneous rules required for maintaining
// compatability between CHEM and STREAM.

// Standard library includes.
#include <string>
using std::string ;

// Loci includes.
#include <Loci.h>

// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;

namespace streamUns {

  // Sets the default background pressure.
  class DefaultP0 : public default_rule {
    private:
      param<double> p0 ;
    public:

      // Define input and output.
      DefaultP0() {
        name_store("p0",p0) ;
        output("p0") ;
        comments("The background pressure used in gage pressure ") ;
        comments("computations.  Pa = Pg+p0.") ;
      }

      // Set the default value.
      virtual void compute(const sequence& seq) { *p0=0 ; }
  } ;

  register_rule<DefaultP0> registerDefaultP0 ;

  // Gets the number of species.
  class GetNumSpecies : public singleton_rule {
    private:
      const_param<EOS> eos ;
      param<int> numSpecies ;
    public:

      GetNumSpecies() {
        name_store("eos",eos) ;
        name_store("numSpecies",numSpecies) ;
        input("eos") ;
        output("numSpecies") ;
      }

      virtual void compute(const sequence &seq) {
        *numSpecies=eos->numSpecies() ;
      }
  } ;

  register_rule<GetNumSpecies> registerGetNumSpecies ;
}
