#ifndef CAVITATIONVARSFILEINPUTS_H
#define CAVITATIONVARSFILEINPUTS_H
                                                                                
// Standard library includes.
#include <iostream>
#include <map>
using std::cerr ;
using std::endl ;

// Loci includes.
#include <Loci.h>
                                                                                
// StreamUns includes.
#include "sciTypes.h"
                                                                                
namespace streamUns {

  class CavitationEquationOptions : public options_list {
    public:
      CavitationEquationOptions() : options_list
        ("model:Cdest:Cprod:rhol:rhov:pv:vInfinite:tInfinite:linearSolver:relaxationFactor:maxIterations") {}
  } ;

}

namespace Loci {

  template<> struct data_schema_traits<streamUns::CavitationEquationOptions> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<streamUns::CavitationEquationOptions>
      Converter_Type ;
  } ;
  
}

#endif
