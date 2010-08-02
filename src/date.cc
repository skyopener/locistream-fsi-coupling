//-----------------------------------------------------------------------------
// Description: This file contains rules for implementing Date's method for
//   the pressure-correction equation.
//
// Authors: Siddharth Thakur and Jeff Wright
//          Streamline Numerics Inc.
//          Gainesville, Florida 32609
//-----------------------------------------------------------------------------

// Standard library includes.
#include <vector>
using std::vector ;

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// StreamUns includes.
#include "residual.h"
#include "sciTypes.h"

namespace streamUns {

  // Rule to compute the mass flux on interior faces for the FOU and
  // SOU convection schemes.
  class StageOneMassFluxInteriorDate : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<vect3d> vStar ;
      const_store<real> faceDensity ;
      const_store<Area> area ;
      store<real> stageOneMassFlux ;
    public:

      // Define input and output.
      StageOneMassFluxInteriorDate() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("vStar",vStar) ;
        name_store("faceDensity",faceDensity) ;
        name_store("area",area) ;
        name_store("stageOneMassFlux",stageOneMassFlux) ;
        input("(cl,cr)->vStar,faceDensity,area") ;
        output("stageOneMassFlux") ;
        constraint("internalFaces,dateMassFlux") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        stageOneMassFlux[face]=0.5*faceDensity[face]*dot(vStar[cl[face]]+
          vStar[cr[face]],area[face].n)*area[face].sada ;
      }
                                                                                
      // Calculate mass flux for all faces in a sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;
                                                                                
  register_rule<StageOneMassFluxInteriorDate>
    registerStageOneMassFluxInteriorDate ;

  // Rule for corrected mass flux on interior faces.
  class CorrectMassFluxInteriorDate : public pointwise_rule {
    private:
      const_store<real> stageOneMassFlux ;
      store<real> massFluxCorrected ;
    public:
                                                                                
      // Define input and output.
      CorrectMassFluxInteriorDate() {
        name_store("stageOneMassFlux",stageOneMassFlux) ;
        name_store("massFluxCorrected",massFluxCorrected) ;
        input("stageOneMassFlux") ;
        output("massFluxCorrected") ;
        constraint("internalFaces,dateMassFlux") ;
      }
                                                                                
      // No correction.
      void calculate(Entity face) {
        massFluxCorrected[face]=stageOneMassFlux[face] ;
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<CorrectMassFluxInteriorDate>
    registerCorrectMassFluxInteriorDate ;
}
