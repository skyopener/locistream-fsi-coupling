// This file contains all rules that are specific to interfacing with the NLAMS
// structural dynamics package.

// Standard library includes.
#include <vector>
using std::vector ;

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// StreamUns includes.
#include "sciTypes.h"

// Global data for each process that will be exchanged between Loci-Stream
// and the structural dynamics solver.
int numNodeTop,numNodeBottom,numCellTop,numCellBottom ;
float *xCTop,*yCTop,*zCTop,*pCTop ;
float *xCBottom,*yCBottom,*zCBottom,*pCBottom ;
float *xTop,*yTop,*zTop,*dXTop,*dYTop,*dZTop ;
float *xBottom,*yBottom,*zBottom,*dXBottom,*dYBottom,*dZBottom ;

// Forward declarations using C linkage to avoid mangling the names of
// the Fortran functions that will be called.
extern "C" {
 void allocatemvbmoduledata_(int*,int*,int*,int*) ;
}

namespace streamUns {

//-----------------------------------------------------------------------------
// Rules for creating a list of the 'top' node numbers. Note that we use a
// temporary variable to assemble the list and then have a permanent variable
// which is the unique and sorted list.

  class TopNodesTempUnit : public unit_rule {
    private:
      param<vector<size_t> > topNodeTemp ;
    public:

      // Define input and output.
      TopNodesTempUnit() {
        name_store("topNodeTemp",topNodeTemp) ;
        output("topNodeTemp") ;
        constraint("UNIVERSE") ;
      }

      // Initialize the vector.
      void compute(const sequence &seq) { topNodeTemp->clear() ; }
  } ;

  register_rule<TopNodesTempUnit> registerTopNodesTempUnit ;

}
