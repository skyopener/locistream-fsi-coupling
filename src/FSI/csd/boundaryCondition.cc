// Standard library includes.
#include <vector>
using std::vector ;

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// StreamUns includes.
#include "gridReader/readGrid.h"
#include "sciTypes.h"

namespace streamUns {

  // Adds 'top' and 'bottom' as valid tags to the no-slip boundary
  // condition. These tags are use to mark the top and bottom of a
  // wing.
  class CheckNoslipCSD : public BC_Check {
    public:
      string BoundaryConditions() { return "noslip" ; }
      string VariablesChecked(fact_db &facts) { return "top,bottom" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        return true ;
      }
      ostream &ErrorMessage(ostream &s) { return s ; }
  } ;

  register_BC<CheckNoslipCSD> registerCheckNoslipCSD ;

}
