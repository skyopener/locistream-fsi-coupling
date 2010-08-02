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

  // Allocates memory for variables in the mvb module.
  void AllocateMVBModuleData() {

    // Call the Fortran 90 function.
    allocatemvbmoduledata_(&numNodeTop,&numNodeBottom,&numCellTop,
      &numCellBottom) ;
  }

  // Passes the top and bottom node coordinates to the CSD solver via
  // the mvb module.
  void NodeCoordinatesToCSD() {

    // Call the Fortran 90 function.
  }

}
