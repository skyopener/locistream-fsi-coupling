#ifndef INITIALCONDITIONREGIONS_HEADER
#define INITIALCONDITIONREGIONS_HEADER
#include <Loci.h>
#include "sciTypes.h"
#include "eos.h"
#include "fluidState.h"
#include "Vector.h"
#include <vector>

namespace streamUns {

  class geomTest : public Loci::CPTR_type {
  public:
    virtual bool inGeomPt(vect3d pt) const = 0 ;
  } ;
  Loci::CPTR<geomTest> geomTestFactory(string name, const options_list ol) ;

  struct ICstate_info {
    Loci::CPTR<geomTest> geomTestFunc ;
    string name ;
    fluidState regionState ;
    std::vector<double> q,qp ;
  } ;

  struct ICparsedInitRegion {
    fluidState defaultState ;
    std::vector<ICstate_info> fluidRegions ;
    std::vector<double> default_q ;
    std::vector<double> default_qp ;
  } ;
  
}

#endif
