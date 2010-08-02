#ifndef FLAMELETVARSFILEINPUTS_H
#define FLAMELETVARSFILEINPUTS_H
                                                                                
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

  class FlameletEquationOptions : public options_list {
    public:
      FlameletEquationOptions() : options_list
        ("table:linearSolver:relaxationFactor:maxIterations") {}
  } ;

}

namespace Loci {

  template<> struct data_schema_traits<streamUns::FlameletEquationOptions> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<streamUns::FlameletEquationOptions>
      Converter_Type ;
  } ;
  
}

#endif
