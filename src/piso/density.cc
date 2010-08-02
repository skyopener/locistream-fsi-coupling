//-----------------------------------------------------------------------------
// Description: This file contains rules for density.
//
// Authors: Siddharth Thakur and Jeff Wright
//          Streamline Numerics Inc.
//          Gainesville, Florida 32609
//-----------------------------------------------------------------------------

// Standard library includes.
#include<vector>
using std::vector ;

// Loci includes.
#include <Loci.h>

// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;

// StreamUns includes.
#include "sciTypes.h"

namespace streamUns {

//-----------------------------------------------------------------------------
// Density boundary condition rules.

  // Rule for boundary faces with specified density. Assigns density value
  // to all boundary faces that have the property rho_BC.
  class BoundaryDensitySpecification : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<real> rho_BC ;
      store<real> rho_f ;
    public:

      // Define input and output.
      BoundaryDensitySpecification() {
        name_store("ref",ref) ;
        name_store("rho_BC",rho_BC) ;
        name_store("rho_f",rho_f) ;
        input("ref->rho_BC") ;
        output("rho_f") ;
      }

      // Calculate density for a single face.
      void calculate(Entity face) { rho_f[face]=rho_BC[ref[face]] ; }

      // Calculate density for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryDensitySpecification>
    registerBoundaryDensitySpecification ;

  // Rule for extrapolating density to boundary faces.
  class BoundaryDensityExtrapolation : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> rho ;
      store<real> rho_f ;
    public:

      // Define input and output.
      BoundaryDensityExtrapolation() {
        name_store("ci",ci) ;
        name_store("rho",rho) ;
        name_store("extrapolated::rho_f",rho_f) ;
        input("ci->rho") ;
        output("extrapolated::rho_f") ;
        constraint("extrapolatedDensity_BC") ;
      }

      // Calculate density for a single face.
      void calculate(Entity face) { rho_f[face]=rho[ci[face]] ; }

      // Calculate density for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryDensityExtrapolation>
    registerBoundaryDensityExtrapolation ;

  // Computes density for compressible flow on boundaries where it has not
  // been otherwise specified. New on 020705.
  class BoundaryDensityComputation: public pointwise_rule {
    private:
      const_store<EOS::State> eos_state_f ;
      store<real> rho_f ;
    public:

      // Define input and output.
      BoundaryDensityComputation() {
        name_store("eos_state_f",eos_state_f) ;
        name_store("rho_f",rho_f) ;
        input("eos_state_f") ;
        output("rho_f") ;
        constraint("eos_state_f,compressibleFlow") ;
      }

      // Calculate density for a single face.
      void calculate(Entity face) { rho_f[face]=eos_state_f[face].density() ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryDensityComputation> registerBoundaryDensityComputation ;

  // Assigns boundary density from interpolated interface value.
  class BoundaryDensityInterpolation : public pointwise_rule {
    private:
      const_store<real> rhoI ;
      store<real> rho_f ;
    public:

      // Define input and output.
      BoundaryDensityInterpolation() {
        name_store("interpolateFace(rho)",rhoI) ;
        name_store("interpolated::rho_f",rho_f) ;
        input("interpolateFace(rho)") ;
        output("interpolated::rho_f") ;
        constraint("interface_BC") ;
      }

      // Calculate density for a single face.
      void calculate(Entity face) { rho_f[face]=rhoI[face] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryDensityInterpolation>
    registerBoundaryDensityInterpolation ;

//-----------------------------------------------------------------------------
// Rules for marching the density.

  // Time build rule for density when using BDF2 time integrator. Although
  // this rule sets rho{n=-1} from rho_ic, these values are not really
  // used since BDF is used on the first timestep for non-restarts.
  class TimeBuildDensityBDF2 : public pointwise_rule {
    private:
      const_store<real> rho_ic ;
      store<real> rho ;
    public:

      // Define input and output.
      TimeBuildDensityBDF2() {
        name_store("rho_ic",rho_ic) ;
        name_store("rho{n=-1}",rho) ;
        input("rho_ic") ;
        output("rho{n=-1}") ;
        constraint("geom_cells") ;
      }

      // Assign density at time zero for a single cell.
      inline void calculate(Entity cell) { rho[cell]=rho_ic[cell] ; }

      // Assign density at time zero for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildDensityBDF2> registerTimeBuildDensityBDF2 ;

  // Time build rule for density.
  class TimeBuildDensity : public pointwise_rule {
    private:
      const_store<real> rho_ic ;
      store<real> rhoTimeStepZero ;
    public:

      // Define input and output.
      TimeBuildDensity() {
        name_store("rho_ic",rho_ic) ;
        name_store("rho{n=0}",rhoTimeStepZero) ;
        input("rho_ic") ;
        output("rho{n=0}") ;
        constraint("geom_cells") ;
      }

      // Assign density at time zero for a single cell.
      inline void calculate(Entity cell) {
        rhoTimeStepZero[cell]=rho_ic[cell] ;
      }

      // Assign density at time zero for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildDensity> registerTimeBuildDensity ;

  // Convenience rule for incompressible flow so that we don't have
  // to have different rules later for compressible and incompressible
  // flow.
  class UpdatedDensityPredictorIncompressible : public pointwise_rule {
    private:
      const_store<real> rho ;
      store<real> rhoStar ;
    public:

      // Define input and output.
      UpdatedDensityPredictorIncompressible() {
        name_store("rho{n}",rho) ;
        name_store("rhoStar{n}",rhoStar) ;
        input("rho{n}") ;
        output("rhoStar{n}") ;
        constraint("incompressibleFlow{n},geom_cells{n}") ;
      }

      // Set for a single cell.
      void calculate(Entity cell) { rhoStar[cell]=rho[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<UpdatedDensityPredictorIncompressible>
    registerUpdatedDensityPredictorIncompressible ;

  // Iteration build rule for density.
  class IterationBuildDensity : public pointwise_rule {
    private:
      const_store<real> rhoStar ;
      store<real> rho ;
    public:

      // Define input and output.
      IterationBuildDensity() {
        name_store("rhoStar{n}",rhoStar) ;
        name_store("rho{n,it=0}",rho) ;
        input("rhoStar{n}") ;
        output("rho{n,it=0}") ;
        constraint("geom_cells{n}") ;
      }

      // Set for a single cell.
      void calculate(Entity cell) { rho[cell]=rhoStar[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationBuildDensity> registerIterationBuildDensity ;

  // Convenience rule for incompressible flow so that we don't have
  // to have different rules later for compressible and incompressible
  // flow.
  class UpdatedDensityCorrectorIncompressible : public pointwise_rule {
    private:
      const_store<real> rho ;
      store<real> rhoCorrected ;
    public:

      // Define input and output.
      UpdatedDensityCorrectorIncompressible() {
        name_store("rho{n,it}",rho) ;
        name_store("rhoCorrected{n,it}",rhoCorrected) ;
        input("rho{n,it}") ;
        output("rhoCorrected{n,it}") ;
        constraint("incompressibleFlow{n,it},geom_cells{n,it}") ;
      }

      // Set for a single cell.
      void calculate(Entity cell) { rhoCorrected[cell]=rho[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<UpdatedDensityCorrectorIncompressible>
    registerUpdatedDensityCorrectorIncompressible ;

  // Iteration advance rule for density for compressible flow.
  class IterationAdvanceDensity : public pointwise_rule {
    private:
      const_store<real> rhoCorrected ;
      store<real> rho ;
    public:

      // Define input and output.
      IterationAdvanceDensity() {
        name_store("rhoCorrected{n,it}",rhoCorrected) ;
        name_store("rho{n,it+1}",rho) ;
        input("rhoCorrected{n,it}") ;
        output("rho{n,it+1}") ;
        constraint("geom_cells{n,it}") ;
      }

      // Assign density at end of iteration for a single cell.
      void calculate(Entity cell) { rho[cell]=rhoCorrected[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationAdvanceDensity> registerIterationAdvanceDensity ;

  // Iteration collapse rule for density.
  class IterationCollapseDensity : public pointwise_rule {
    private:
      store<real> rho ;
    public:

      // Define input and output.
      IterationCollapseDensity() {
        name_store("rho{n,it}",rho) ;
        input("rho{n,it}") ;
        output("rho{n+1}=rho{n,it}") ;
        conditional("iterationFinished{n,it-1}") ;
        constraint("geom_cells{n,it}") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<IterationCollapseDensity> registerIterationCollapseDensity ;
}
