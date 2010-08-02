// Standard library includes.
#include <string>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <math.h>

using std::string ;
using std::vector ;

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;

// StreamUns includes.
#include "sciTypes.h"

namespace streamUns {

// Rules for the computing the Chi

  class ChiInterior : public pointwise_rule {
    private:
      const_store<real> omega,Zvar ;
      store<real> Chi ;
    public:

      // Define input and output.
      ChiInterior() {
        name_store("omega",omega) ;
        name_store("Zvar",Zvar) ;
        name_store("Chi",Chi) ;
        input("omega,Zvar") ;
        output("Chi") ;
        constraint("geom_cells,flameletModel") ;
      }

      // Compute for a single cell.
      void calculate(Entity cell) {
      	Chi[cell]=log10(2.0*0.09*omega[cell]*Zvar[cell]+1.e-30); 
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ChiInterior> registerChiInterior ;

  class ChiBoundary : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> Chi ;
      store<real> Chi_f ;
    public:

      // Define input and output.
      ChiBoundary() {
        name_store("ci",ci) ;
        name_store("Chi",Chi) ;
        name_store("Chi_f",Chi_f) ;
        input("ci->Chi") ;
        output("Chi_f") ;
        constraint("boundaryFaces,flameletModel") ;
      }

      // Compute for a single face. Simple extrapolation for now.
      void calculate(Entity face) {
        Chi_f[face]=Chi[ci[face]] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ChiBoundary> registerChiBoundary ;

}
