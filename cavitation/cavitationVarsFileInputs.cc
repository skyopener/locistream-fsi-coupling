// Loci includes.
#include <Loci.h>
                                                                                
// StreamUns includes.
#include "const.h"
#include "cavitationInitialCondition.h"
#include "sciTypes.h"
#include "cavitationVarsFileInputs.h"
                                                                                
namespace streamUns {

  // Options for the cavitation equations.
  class OptionalCavitationEquations : public optional_rule {
    private:
      param<CavitationEquationOptions> cavitationEquationOptions ;
    public:

      // Define input and output.
      OptionalCavitationEquations() {
        name_store("cavitationEquationOptions",cavitationEquationOptions) ;
        output("cavitationEquationOptions") ;
      }

      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;

  register_rule<OptionalCavitationEquations> registerOptionalCavitationEquations ;

  
  class OptionalCavitationInitialCondition : public optional_rule {
    private:
      param<CavitationInitialCondition> cavitationInitialCondition ;
    public:

      // Define input and output.
      OptionalCavitationInitialCondition() {
        name_store("cavitationInitialCondition",cavitationInitialCondition) ;
        output("cavitationInitialCondition") ;
      }

      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;

  register_rule<OptionalCavitationInitialCondition> registerOptionalCavitationInitialCondition ;

  
  // Sets the default cavitation inviscid flux.
  class DefaultCavitationInviscidFlux : public default_rule {
    private:
      param<string> cavitationInviscidFlux ;
    public:

      // Define input and output.
      DefaultCavitationInviscidFlux() {
        name_store("cavitationInviscidFlux",cavitationInviscidFlux) ;
        output("cavitationInviscidFlux") ;           
        comments("Sets the default cavitationInviscidFlux to 'FOU'.") ;
        comments("Valid options are 'FOU' and 'SOU'.") ;
      }    

      // Set the default value.  
      virtual void compute(const sequence& seq) {
        *cavitationInviscidFlux="FOU" ;               
      }
  } ;

  register_rule<DefaultCavitationInviscidFlux>
    registerDefaultCavitationInviscidFlux ;
  
}
