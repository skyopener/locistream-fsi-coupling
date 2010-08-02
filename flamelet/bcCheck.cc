// Loci includes.
#include <Loci.h>
                                                                                
// StreamUns includes.
#include "gridReader/readGrid.h"
#include "sciTypes.h"
                                                                                
namespace streamUns {
	
  // Checks incompressible inlet boundaries.
  class FlameletCheckIncompressibleInlet : public BC_Check {
    private:
      string errorMessage ;
    public:
      FlameletCheckIncompressibleInlet() : errorMessage("") {}
    public:
      string BoundaryConditions() { return "incompressibleInlet" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        param<string> flowRegime=facts.get_fact("flowRegime") ;
        if(*flowRegime!="turbulent"){
          errorMessage="flowRegime must be turbulent for flamelet model." ;
          return false ;
        }
        
        if(!bc_options.optionExists("Z")) {
          errorMessage="Must specify 'Z'." ;
          return false ;
        }
        
        if(!bc_options.optionExists("Zvar")) {
          errorMessage="Must specify 'Zvar'." ;
          return false ;
        }
        
        return true ;
      }
      ostream &ErrorMessage(std::ostream &s) {
        s << errorMessage << endl ; return s ;
      }
      string VariablesChecked(fact_db &facts) {
        string s="Z,Zvar" ;
        return s ;
      }
  } ;

  register_BC<FlameletCheckIncompressibleInlet> registerFlameletCheckIncompressibleInlet ;
  

}