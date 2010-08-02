//-----------------------------------------------------------------------------
// Description: This file contains extra source terms required for using BDF2
//   with deforming meshes.
//-----------------------------------------------------------------------------

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// StreamUns includes.
#include "const.h"
#include "sciTypes.h"

namespace streamUns {

  // Rule to add extra source to replace vol{n,it} with vol{n-1} in temporal
  // source term for the BDF2 scheme. This term is currently activated by the
  // presence of a gridMover that is specified in the .vars file. If there is
  // no grid mover specified, then the code assumes that no mesh deformation
  // will be occurring.
  class TemporalToVelocitySourceTermDeformingMeshBDF2 : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_param<string> gridMover ;
      const_param<real> timeStep ;
      const_param<real> timeIntegratorFactor1 ;
      const_store<real> timeStepFactor ;
      const_store<real> rhoOld ;
      const_store<vect3d> vOld,v ;
      const_store<real> volOld,vol ;
      const_store<real> cellRadius ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      TemporalToVelocitySourceTermDeformingMeshBDF2() {
        name_store("gridMover{n,it}",gridMover) ;
        name_store("timeStep{n,it}",timeStep) ;
        name_store("timeIntegratorFactor1{n}",timeIntegratorFactor1) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n-1}",rhoOld) ;
        name_store("v{n-1}",vOld) ;
        name_store("v{n,it}",v) ;
        name_store("vol{n-1}",volOld) ;
        name_store("vol{n,it}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("vSourceTerm{n,it}",vSourceTerm) ;
        input("gridMover{n,it}") ;
        input("rho{n-1},v{n-1},v{n,it},vol{n-1},vol{n,it},timeStep{n,it}") ;
        input("cellRadius{n,it},timeStepFactor{n},timeIntegratorFactor1{n}") ;
        output("vSourceTerm{n,it}") ;
        constraint("geom_cells,BDF2Integrator") ;
      }

      // Add temporal component to source term for a single cell.
      void calculate(Entity cell) {
        vSourceTerm[cell]+=(0.5*(*timeIntegratorFactor1)*rhoOld[cell]*
         (volOld[cell]-vol[cell])*cellRadius[cell]/((*timeStep)*
         timeStepFactor[cell]))*(v[cell]-vOld[cell]) ;
      }

      // Add temporal component to source term for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToVelocitySourceTermDeformingMeshBDF2>
    registerTemporalToVelocitySourceTermDeformingMeshBDF2 ;

  // Rule to add extra source to replace vol{n,it} with vol{n-1} in temporal
  // source term for the BDF2 scheme. This term is currently activated by the
  // presence of a gridMover that is specified in the .vars file. If there is
  // no grid mover specified, then the code assumes that no mesh deformation
  // will be occurring.
  class TemporalToKSourceTermDeformingMeshBDF2 : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_param<string> gridMover ;
      const_param<real> timeStep ;
      const_param<real> timeIntegratorFactor1 ;
      const_store<real> timeStepFactor ;
      const_store<real> rhoOld ;
      const_store<real> kOld,k ;
      const_store<real> volOld,vol ;
      const_store<real> cellRadius ;
      store<real> kSourceTerm ;
    public:

      // Define input and output.
      TemporalToKSourceTermDeformingMeshBDF2() {
        name_store("gridMover{n,it}",gridMover) ;
        name_store("timeStep{n,it}",timeStep) ;
        name_store("timeIntegratorFactor1{n}",timeIntegratorFactor1) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n-1}",rhoOld) ;
        name_store("k{n-1}",kOld) ;
        name_store("k{n,it}",k) ;
        name_store("vol{n-1}",volOld) ;
        name_store("vol{n,it}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("kSourceTerm{n,it}",kSourceTerm) ;
        input("gridMover{n,it}") ;
        input("rho{n-1},k{n-1},k{n,it},vol{n-1},vol{n,it},timeStep{n,it}") ;
        input("cellRadius{n,it},timeStepFactor{n},timeIntegratorFactor1{n}") ;
        output("kSourceTerm{n,it}") ;
        constraint("geom_cells,BDF2Integrator") ;
      }

      // Add temporal component to source term for a single cell.
      void calculate(Entity cell) {
        kSourceTerm[cell]+=(0.5*(*timeIntegratorFactor1)*rhoOld[cell]*
          (volOld[cell]-vol[cell])*cellRadius[cell]/((*timeStep)*
          timeStepFactor[cell]))*(k[cell]-kOld[cell]) ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToKSourceTermDeformingMeshBDF2>
    registerTemporalToKSourceTermDeformingMeshBDF2 ;

  class TemporalToOmegaSourceTermDeformingMeshBDF2 : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_param<string> gridMover ;
      const_param<real> timeStep ;
      const_param<real> timeIntegratorFactor1 ;
      const_store<real> timeStepFactor ;
      const_store<real> rhoOld ;
      const_store<real> omegaOld,omega ;
      const_store<real> volOld,vol ;
      const_store<real> cellRadius ;
      store<real> omegaSourceTerm ;
    public:

      // Define input and output.
      TemporalToOmegaSourceTermDeformingMeshBDF2() {
        name_store("gridMover{n,it}",gridMover) ;
        name_store("timeStep{n,it}",timeStep) ;
        name_store("timeIntegratorFactor1{n}",timeIntegratorFactor1) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n-1}",rhoOld) ;
        name_store("omega{n-1}",omegaOld) ;
        name_store("omega{n,it}",omega) ;
        name_store("vol{n-1}",volOld) ;
        name_store("vol{n,it}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("omegaSourceTerm{n,it}",omegaSourceTerm) ;
        input("gridMover{n,it},rho{n-1},omega{n-1},omega{n,it}") ;
        input("vol{n-1},vol{n,it},timeStep{n,it}") ;
        input("cellRadius{n,it},timeStepFactor{n},timeIntegratorFactor1{n}") ;
        output("omegaSourceTerm{n,it}") ;
        constraint("geom_cells,BDF2Integrator") ;
      }

      // Add temporal component to source term for a single cell.
      void calculate(Entity cell) {
        omegaSourceTerm[cell]+=(0.5*(*timeIntegratorFactor1)*rhoOld[cell]*
          (volOld[cell]-vol[cell])*cellRadius[cell]/((*timeStep)*
          timeStepFactor[cell]))*(omega[cell]-omegaOld[cell]) ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToOmegaSourceTermDeformingMeshBDF2>
    registerTemporalToOmegaSourceTermDeformingMeshBDF2 ;

}
