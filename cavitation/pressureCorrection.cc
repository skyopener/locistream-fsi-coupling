// Standard library includes.
#include <fstream>

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;
    
// StreamUns includes.
#include "residual.h"
#include "sciTypes.h"
#include "varsFileInputs.h"
       
namespace streamUns {

  class PressureCorrectionDiagonalCavitationInternal : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> rho ;
      const_store<real> Alpha ;
      const_store<real> massFlux ;
      store<real> D ;
    public:
 
       PressureCorrectionDiagonalCavitationInternal() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("rho",rho) ;
        name_store("Alpha",Alpha) ;
        name_store("pPrime_D",D) ;
        name_store("stageOneMassFlux",massFlux) ;
        input("(cl,cr)->(rho,Alpha)") ;
        input("stageOneMassFlux") ;
        output("(cl,cr)->pPrime_D") ;
        constraint("internalFaces,incompressibleFlow,cavitationModel") ;
      }

      void calculate(Entity face) {
        if(massFlux[face]>0.0){
          real cRho=4.0*(1.0-Alpha[cl[face]]) ;
          real compressibleCoefficient=cRho*massFlux[face]/rho[cl[face]] ;
          D[cl[face]]+=compressibleCoefficient ;
        }else{
          real cRho=4.0*(1.0-Alpha[cr[face]]) ;
          real compressibleCoefficient=-cRho*massFlux[face]/rho[cr[face]] ;
          D[cr[face]]+=compressibleCoefficient ;
        }
      }

      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
 
  register_rule<PressureCorrectionDiagonalCavitationInternal>
    registerPressureCorrectionDiagonalCavitationInternal ;

  class PressureCorrectionDiagonalCavitationBoundary : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> rho ;
      const_store<real> Alpha ;
      const_store<real> massFlux ;
      store<real> D ;
    public:
 
       PressureCorrectionDiagonalCavitationBoundary() {
        name_store("ci",ci) ;
        name_store("rho",rho) ;
        name_store("Alpha",Alpha) ;
        name_store("stageOneMassFlux",massFlux) ;
        name_store("pPrime_D",D) ;
        input("ci->(rho,Alpha),stageOneMassFlux") ;
        output("ci->pPrime_D") ;
        constraint("nonSpecifiedMassFlux_BC,cavitationModel") ;
      }

      void calculate(Entity face) {
        if(massFlux[face]>0.0){
          real cRho=4.0*(1.0-Alpha[ci[face]]) ;
          D[ci[face]]+=cRho*massFlux[face]/rho[ci[face]] ;
        }
      }
 
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
 
  register_rule<PressureCorrectionDiagonalCavitationBoundary>
    registerPressureCorrectionDiagonalCavitationBoundary ;
 
  // Pseudo-compressible part of pressure correction equation.
  class ComputePressureCorrectionLowerUpperCavitation : public
  pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<real> rho ;
      const_store<real> Alpha ;
      const_store<real> pPrimeCoefficient ;
      const_store<real> massFlux ;
      store<real> L,U ;
    public:
 
       // Define input and output.
       ComputePressureCorrectionLowerUpperCavitation() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("rho",rho) ;
        name_store("Alpha",Alpha) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        name_store("stageOneMassFlux",massFlux) ;
        name_store("cavitation::pPrime_L",L) ;
        name_store("cavitation::pPrime_U",U) ;
        input("(cl,cr)->(rho,Alpha)") ;
        input("pPrimeCoefficient,stageOneMassFlux") ;
        output("cavitation::pPrime_L,cavitation::pPrime_U") ;
        constraint("internalFaces,incompressibleFlow,cavitationModel") ;
      }

      void calculate(Entity face) {
        if(massFlux[face]>0.0){
          real cRho=4.0*(1.0-Alpha[cl[face]]) ;
	  real compressibleCoefficient=cRho*massFlux[face]/rho[cl[face]] ;
          L[face]=-pPrimeCoefficient[face]-compressibleCoefficient ;
          U[face]=-pPrimeCoefficient[face] ;
        }else{
          real cRho=4.0*(1.0-Alpha[cr[face]]) ;
	  real compressibleCoefficient=-cRho*massFlux[face]/rho[cr[face]] ;
          L[face]=-pPrimeCoefficient[face] ;
          U[face]=-pPrimeCoefficient[face]-compressibleCoefficient ;
        }
      }

      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
 
  register_rule<ComputePressureCorrectionLowerUpperCavitation>
    registerComputePressureCorrectionLowerUpperCavitation ;
 
}
