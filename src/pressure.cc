//-----------------------------------------------------------------------------
// Description: This file contains rules for pressure.
//
// Authors: Siddharth Thakur and Jeff Wright
//          Streamline Numerics Inc.
//          Gainesville, Florida 32609
//-----------------------------------------------------------------------------

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// StreamUns includes.
#include "sciTypes.h"

namespace streamUns {

//-----------------------------------------------------------------------------
// Pressure boundary condition rules.

  // Assigns boundary pressure from interpolated interface value. This is only
  // done for compressible flow. For incompressible flow, pressure is
  // extrapolated to the boundary.
  class BoundaryPressureInterpolation : public pointwise_rule {
    private:
      const_store<real> pI ;
      store<real> p_f ;
    public:

      // Define input and output.
      BoundaryPressureInterpolation() {
        name_store("interpolateFace(p)",pI) ;
        name_store("p_f",p_f) ;
        input("interpolateFace(p)") ;
        output("p_f") ;
        constraint("interface_BC") ;
//      constraint("interface_BC,compressibleFlow") ; ORIGINAL
      }

      // Calculate pressure for a single face.
      void calculate(Entity face) { p_f[face]=pI[face] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryPressureInterpolation>
    registerBoundaryPressureInterpolation ;

  // Rule for boundary faces with specified pressure. Assigns pressure value
  // to all boundary faces that have the property p_BC.
  class BoundaryPressureSpecification : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<real> p_BC ;
      store<real> p_f ;
    public:

      // Define input and output.
      BoundaryPressureSpecification() {
        name_store("ref",ref) ;
        name_store("p_BC",p_BC) ;
        name_store("p_f",p_f) ;
        input("ref->p_BC") ;
        output("p_f") ;
      }

      // Calculate pressure for a single face.
      void calculate(Entity face) { p_f[face]=p_BC[ref[face]] ; }

      // Calculate pressure for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryPressureSpecification>
    registerBoundaryPressureSpecification ;

  // Rule for boundary faces with specified pressure. Assigns pressure value
  // to all boundary faces that have the property p_BC. This priority version
  // converts to a total pressure inlet when backflow occurs.
  class BoundaryPressureSpecificationEntrainment : public pointwise_rule {
    private:
      const_Map ci,ref ;
      const_store<real> p_BC ;
      const_store<real> rho ;
      const_store<vect3d> v ;
      const_store<Area> area ;
      store<real> p_f ;
    public:

      // Define input and output.
      BoundaryPressureSpecificationEntrainment() {
        name_store("ci",ci) ;
        name_store("ref",ref) ;
        name_store("p_BC",p_BC) ;
        name_store("rho",rho) ;
        name_store("v",v) ;
        name_store("area",area) ;
        name_store("entrainment::p_f",p_f) ;
        input("ref->p_BC,ci->(rho,v),area") ;
        output("entrainment::p_f") ;
        constraint("ref->entrainment_BCoption") ;
      }

      // Calculate pressure for a single face.
      void calculate(Entity face) {
        p_f[face]=(dot(v[ci[face]],area[face].n)>=0.0)? p_BC[ref[face]]:
          p_BC[ref[face]]-0.5*rho[ci[face]]*dot(v[ci[face]],v[ci[face]]) ;
      }

      // Calculate pressure for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryPressureSpecificationEntrainment>
    registerBoundaryPressureSpecificationEntrainment ;

  // Unit rule to sum extrapolated pressure times area.
  class BoundaryPressureAreaUnit : public unit_rule {
    private:
      param<real> boundaryPressureAreaSum ;
    public:

      // Define input and output.
      BoundaryPressureAreaUnit() {
        name_store("boundaryPressureAreaSum(X)",boundaryPressureAreaSum) ;
        output("boundaryPressureAreaSum(X)") ;
        constraint("X") ;
      }

      // Initialize the value.
      void compute(const sequence &seq) { *boundaryPressureAreaSum=0.0 ; }
  } ;

  register_rule<BoundaryPressureAreaUnit> registerBoundaryPressureAreaUnit ;

  // Unit rule to sum extrapolated pressure times area.
  class BoundaryPressureAreaApply : public apply_rule<param<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> p ;
      const_store<Area> area ;
      param<real> boundaryPressureAreaSum ;
    public:

      // Define input and output.
      BoundaryPressureAreaApply() {
        name_store("ci",ci) ;
        name_store("p",p) ;
        name_store("area",area) ;
        name_store("boundaryPressureAreaSum(X)",boundaryPressureAreaSum) ;
        input("ci->p,area") ;
        output("boundaryPressureAreaSum(X)") ;
        constraint("X") ;
      }

      // Add the face value to the total.
      void calculate(Entity face) {
        join(*boundaryPressureAreaSum,p[ci[face]]*area[face].sada) ;
      }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryPressureAreaApply> registerBoundaryPressureAreaApply ;

  // Rule for boundary faces with specified mean pressure. Assigns pressure
  // value to all boundary faces that have the property pMean_BC.
  class BoundaryMeanPressureSpecification : public pointwise_rule {
    private:
      const_Map ref,ci ;
      const_param<real> boundaryPressureAreaSum ;
      const_param<real> totalArea ;
      const_store<real> p ;
      const_store<real> pMean_BC ;
      store<real> p_f ;
    private:
      real pC ;
    public:

      // Define input and output.
      BoundaryMeanPressureSpecification() {
        name_store("ref",ref) ;
        name_store("ci",ci) ;
        name_store("boundaryPressureAreaSum(fixedPressureOutlet_BC)",
          boundaryPressureAreaSum) ;
        name_store("totalArea(fixedPressureOutlet_BC)",totalArea) ;
        name_store("p",p) ;
        name_store("pMean_BC",pMean_BC) ;
        name_store("p_f",p_f) ;
        input("boundaryPressureAreaSum(fixedPressureOutlet_BC)") ;
        input("totalArea(fixedPressureOutlet_BC)") ;
        input("ref->pMean_BC,ci->p") ;
        output("p_f") ;
      }

      // Extrapolate and add correction for a face. Compute the pressure
      // correction which is added to all faces so that the mean pressure is
      // pMean. IMPORTANT: There is an implicit assumption here that pMean is
      // the same for all fixedPressureOutlet boundaries. Tried to put this in
      // compute(), but was getting seg. fault when using seq.begin() as the
      // first face in the sequence for multi-process runs.
      void calculate(Entity face) {
        pC=pMean_BC[ref[face]]-(*boundaryPressureAreaSum)/(*totalArea) ;
        p_f[face]=p[ci[face]]+pC ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryMeanPressureSpecification>
    registerBoundaryMeanPressureSpecification ;

  // Rule for extrapolating pressure to boundary faces. This occurs for
  // incompressible and subsonic inlets, extrapolated pressure outlets,
  // and no-slip and slip boundaries.
  class BoundaryPressureExtrapolation : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> p ;
      store<real> p_f ;
    public:

      // Define input and output.
      BoundaryPressureExtrapolation() {
        name_store("ci",ci) ;
        name_store("p",p) ;
        name_store("p_f",p_f) ;
        input("ci->p") ;
        output("p_f") ;
        constraint("extrapolatedPressure_BC") ;
      }

      // Calculate pressure for a single face.
      void calculate(Entity face) { p_f[face]=p[ci[face]] ; }

      // Calculate pressure for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryPressureExtrapolation>
    registerBoundaryPressureExtrapolation ;

  // Priority rule for computing pressure from total-pressure for
  // total-pressure inlets.
  class BoundaryPressureTotalPressureInletCompressible : public pointwise_rule {
    private:
      const_Map ci ;
      const_Map ref ;
      const_store<real> cp,cv ;
      const_store<real> machNumber ;
      const_store<real> p0_BC ;
      store<real> p_f ;
    public:

      // Define input and output.
      BoundaryPressureTotalPressureInletCompressible() {
        name_store("ci",ci) ;
        name_store("ref",ref) ;
        name_store("cp",cp) ;
        name_store("cv",cv) ;
        name_store("machNumber(v)",machNumber) ;
        name_store("p0_BC",p0_BC) ;
        name_store("totalPressureInlet::p_f",p_f) ;
        input("ci->(cp,cv,machNumber(v)),ref->p0_BC") ;
        output("totalPressureInlet::p_f") ;
        constraint("compressibleFlow,totalPressureInlet_BC") ;
      }

      // Calculate pressure for a single face.
      void calculate(Entity face) {
        real gamma=cp[ci[face]]/cv[ci[face]] ;
        p_f[face]=p0_BC[ref[face]]/pow(1.0+0.5*(gamma-1.0)*machNumber[ci[face]]*
          machNumber[ci[face]],gamma/(gamma-1.0)) ;
      }

      // Calculate pressure for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  //register_rule<BoundaryPressureTotalPressureInletCompressible>
  //  registerBoundaryPressureTotalPressureInletCompressible ;

//-----------------------------------------------------------------------------
// Rules for marching the pressure.

  // Time build rule for pressure when using BDF2 time integrator. Although
  // this rule sets p{n=-1} from p_ic, these values are not really used since
  // BDF is used on the first timestep for non-restarts.
  class TimeBuildPressureBDF2: public pointwise_rule {
    private:
      const_store<real> p_ic ;
      store<real> p ;
    public:

      // Define input and output.
      TimeBuildPressureBDF2() {
        name_store("p_ic",p_ic) ;
        name_store("p{n=-1}",p) ;
        input("p_ic") ;
        output("p{n=-1}") ;
        constraint("geom_cells") ;
      }

      // Assign p for a single cell.
      void calculate(Entity cell) { p[cell]=p_ic[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildPressureBDF2> registerTimeBuildPressureBDF2 ;

// Time build rule for pressure.
  class TimeBuildPressure : public pointwise_rule {
    private:
      const_store<real> p_ic ;
      store<real> pTimeStepZero ;
    public:

      // Define input and output.
      TimeBuildPressure() {
        name_store("p_ic",p_ic) ;
        name_store("p{n=0}",pTimeStepZero) ;
        input("p_ic") ;
        output("p{n=0}") ;
        constraint("geom_cells") ;
      }

      // Assign pressure at time zero for a single cell.
      void calculate(Entity cell) { pTimeStepZero[cell]=p_ic[cell] ; }

      // Assign pressure at time zero for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildPressure> registerTimeBuildPressure ;

// Iteration build rule for pressure.
  class IterationBuildPressure: public pointwise_rule {
    private:
      const_store<real> pTimeStepN ;
      store<real> pIterationZero ;
    public:

      // Define input and output.
      IterationBuildPressure() {
        name_store("p{n}",pTimeStepN) ;
        name_store("p{n,it=0}",pIterationZero) ;
        input("p{n}") ;
        output("p{n,it=0}") ;
        constraint("geom_cells{n}") ;
      }

      // Assign pressure at iteration zero for a single cell.
      void calculate(Entity cell) { pIterationZero[cell]=pTimeStepN[cell] ; }

      // Assign pressure at iteration zero for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationBuildPressure> registerIterationBuildPressure ;

  // Iteration advance rule for pressure.
  class IterationAdvancePressure : public pointwise_rule {
    private:
      const_store<real> pCorrected ;
      store<real> pIterationPlusOne ;
    public:

      // Define input and output.
      IterationAdvancePressure() {
        name_store("pCorrected{n,it}",pCorrected) ;
        name_store("p{n,it+1}",pIterationPlusOne) ;
        input("pCorrected{n,it}") ;
        output("p{n,it+1}") ;
        constraint("geom_cells{n,it}") ;
      }

      // Assign pressure at end of iteration for a single cell.
      void calculate(Entity cell) { pIterationPlusOne[cell]=pCorrected[cell] ; }

      // Assign pressure at end of iteration for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationAdvancePressure> registerIterationAdvancePressure ;

  // Iteration collapse rule for pressure. Note that we have changed
  // this rule from "p{n+1}=p{n,it}" to the current form. With the old form
  // we could not do one iteration per time step, because p{n,it} would not
  // get updated and thus the residuals would never change.
  class IterationCollapsePressure : public pointwise_rule {
    private:
      store<real> p ;
    public:

      // Define input and output.
      IterationCollapsePressure() {
        name_store("p{n,it}",p) ;
        input("p{n,it}") ;
        output("p{n+1}=p{n,it}") ;
        conditional("iterationFinished{n,it-1}") ;
        constraint("geom_cells{n,it}") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<IterationCollapsePressure> registerIterationCollapsePressure ;
}
