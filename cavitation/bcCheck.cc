// Loci includes.
#include <Loci.h>
                                                                                
// StreamUns includes.
#include "gridReader/readGrid.h"
#include "sciTypes.h"
                                                                                
namespace streamUns {
	
  // Checks incompressible inlet boundaries.
  class CavitationCheckIncompressibleInlet : public BC_Check {
    private:
      string errorMessage ;
    public:
      CavitationCheckIncompressibleInlet() : errorMessage("") {}
    public:
      string BoundaryConditions() { return "incompressibleInlet" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        param<string> flowRegime=facts.get_fact("flowRegime") ;
        if(*flowRegime!="turbulent"){
          errorMessage="flowRegime must be turbulent for cavitation model." ;
          return false ;
        }
        
        if(!bc_options.optionExists("Alpha")) {
          errorMessage="Must specify 'Alpha'." ;
          return false ;
        }
        
        return true ;
      }
      ostream &ErrorMessage(std::ostream &s) {
        s << errorMessage << endl ; return s ;
      }
      string VariablesChecked(fact_db &facts) {
        string s="Alpha" ;
        return s ;
      }
  } ;

  register_BC<CavitationCheckIncompressibleInlet> registerCavitationCheckIncompressibleInlet ;
  

}