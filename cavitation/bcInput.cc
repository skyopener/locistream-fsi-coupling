// Loci includes.
#include <Loci.h>
using Loci::Area ;

// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;

// StreamUns includes.
#include "sciTypes.h"

namespace streamUns {

  class BC_Alpha_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<real> Alpha_BC ;
    public:

      BC_Alpha_compute() {
        name_store("BC_options",BC_options) ;
        name_store("Alpha_BC",Alpha_BC) ;
        input("BC_options") ;
        output("Alpha_BC") ;
        constraint("Alpha_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOption("Alpha",Alpha_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_Alpha_compute> register_BC_Alpha_compute ;

}
