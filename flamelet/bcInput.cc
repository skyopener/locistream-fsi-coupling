// Loci includes.
#include <Loci.h>
using Loci::Area ;

// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;

// StreamUns includes.
#include "sciTypes.h"

namespace streamUns {

  class BC_Z_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<real> Z_BC ;
    public:

      BC_Z_compute() {
        name_store("BC_options",BC_options) ;
        name_store("Z_BC",Z_BC) ;
        input("BC_options") ;
        output("Z_BC") ;
        constraint("Z_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOption("Z",Z_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_Z_compute> register_BC_Z_compute ;

  class BC_Zvar_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<real> Zvar_BC ;
    public:

      BC_Zvar_compute() {
        name_store("BC_options",BC_options) ;
        name_store("Zvar_BC",Zvar_BC) ;
        input("BC_options") ;
        output("Zvar_BC") ;
        constraint("Zvar_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOption("Zvar",Zvar_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_Zvar_compute> register_BC_Zvar_compute ;

}
