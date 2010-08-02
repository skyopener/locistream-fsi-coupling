// Loci includes.
#include <Loci.h>
                                                                                
// StreamUns includes.
#include "const.h"
#include "flameletInitialCondition.h"
#include "sciTypes.h"
#include "flameletVarsFileInputs.h"
                                                                                
namespace streamUns {

  // Options for the flamelet equations.
  class OptionalFlameletEquations : public optional_rule {
    private:
      param<FlameletEquationOptions> flameletEquationOptions ;
    public:

      // Define input and output.
      OptionalFlameletEquations() {
        name_store("flameletEquationOptions",flameletEquationOptions) ;
        output("flameletEquationOptions") ;
      }

      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;

  register_rule<OptionalFlameletEquations> registerOptionalFlameletEquations ;

  
  class OptionalFlameletInitialCondition : public optional_rule {
    private:
      param<FlameletInitialCondition> flameletInitialCondition ;
    public:

      // Define input and output.
      OptionalFlameletInitialCondition() {
        name_store("flameletInitialCondition",flameletInitialCondition) ;
        output("flameletInitialCondition") ;
      }

      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;

  register_rule<OptionalFlameletInitialCondition> registerOptionalFlameletInitialCondition ;

  
  // Sets the default flamelet inviscid flux.
  class DefaultFlameletInviscidFlux : public default_rule {
    private:
      param<string> flameletInviscidFlux ;
    public:

      // Define input and output.
      DefaultFlameletInviscidFlux() {
        name_store("flameletInviscidFlux",flameletInviscidFlux) ;
        output("flameletInviscidFlux") ;           
        comments("Sets the default flameletInviscidFlux to 'FOU'.") ;
        comments("Valid options are 'FOU' and 'SOU'.") ;
      }    

      // Set the default value.  
      virtual void compute(const sequence& seq) {
        *flameletInviscidFlux="FOU" ;               
      }
  } ;

  register_rule<DefaultFlameletInviscidFlux>
    registerDefaultFlameletInviscidFlux ;
  
}
