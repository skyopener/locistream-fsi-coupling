//-----------------------------------------------------------------------------
// Description: This file contains rules for the inviscid and viscous fluxes.
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

  // Creates the PISO constraint, which will activate certain features in the
  // gridMotion module.
  class PISOConstraint : public constraint_rule {
    private:
      Constraint PISO ;
    public:

      // Define input and output.
      PISOConstraint() {
        name_store("PISO",PISO) ;
        output("PISO") ;
        constraint("UNIVERSE") ;
      }

      // Set up the constraint.
      virtual void compute(const sequence& seq) { PISO=~EMPTY ; }
  } ;

  register_rule<PISOConstraint> registerPISOConstraint ;

  // Default rule for grid mass flux. We must have a separate variable for this
  // so we can be GCL on the first iteration. The variable massFlux will now
  // only contain the flux due to the fluid velocity.
  class DefaultGridMassFlux : public pointwise_rule {
    private:
      store<real> gridMassFlux ;
    public:

      // Define input and output.
      DefaultGridMassFlux() {
        name_store("gridMassFlux",gridMassFlux) ;
        output("gridMassFlux") ;
        constraint("faces") ;
      }

      // Set to zero.
      void calculate(Entity face) { gridMassFlux[face]=0.0 ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DefaultGridMassFlux> registerDefaultGridMassFlux ;

  // Rule to create variable internalFaces which is used solely as a constraint.
  class InternalFaces : public pointwise_rule {
    private:
      store<bool> internalFaces ;
    public:

      // Define input and output.
      InternalFaces() {
        name_store("internalFaces",internalFaces) ;
        output("internalFaces") ;
        constraint("(cl,cr)->cellcenter") ;
      }

      // Assign flag for a single internal face.
      void calculate(Entity face) {
        internalFaces[face]=true ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InternalFaces> registerInternalFaces ;

  // Rule to create variable boundaryFaces which is used solely as a constraint.
  class BoundaryFaces : public pointwise_rule {
    private:
      store<bool> boundaryFaces ;
    public:

      // Define input and output.
      BoundaryFaces() {
        name_store("boundaryFaces",boundaryFaces) ;
        output("boundaryFaces") ;
        constraint("ci->geom_cells") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { boundaryFaces[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryFaces> registerBoundaryFaces ;

  // Rule to create variable boundaryMassFluxCorrected. This variable, which is
  // solely used as a constraint identifies boundaries over which the mass flux
  // is corrected.
  class BoundaryMassFluxCorrectedSpecifiedPressure : public pointwise_rule {
    private:
      store<bool> boundaryMassFluxCorrected ;
    public:

      // Define input and output.
      BoundaryMassFluxCorrectedSpecifiedPressure() {
        name_store("boundaryMassFluxCorrected",boundaryMassFluxCorrected) ;
        output("boundaryMassFluxCorrected") ;
        constraint("nonSpecifiedMassFlux_BC,specifiedPressure_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<BoundaryMassFluxCorrectedSpecifiedPressure>
    registerBoundaryMassFluxCorrectedSpecifiedPressure ;

  // This and subsequent rules create the variable boundaryMassFluxNotCorrected.
  // Boundaries with specified velocity do not have mass flux corrected.
  class BoundaryMassFluxNotCorrectedSpecifiedVelocity : public pointwise_rule {
    private:
      store<bool> boundaryMassFluxNotCorrected ;
    public:

      // Define input and output.
      BoundaryMassFluxNotCorrectedSpecifiedVelocity() {
        name_store("boundaryMassFluxNotCorrected",boundaryMassFluxNotCorrected) ;
        output("boundaryMassFluxNotCorrected") ;
        constraint("ref->v_BCoption") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<BoundaryMassFluxNotCorrectedSpecifiedVelocity>
    registerBoundaryMassFluxNotCorrectedSpecifiedVelocity ;

  // Boundaries with specified mass flux do not have mass flux corrected.
  class BoundaryMassFluxNotCorrectedSpecifiedMassFlux : public pointwise_rule {
    private:
      store<bool> boundaryMassFluxNotCorrected ;
    public:

      // Define input and output.
      BoundaryMassFluxNotCorrectedSpecifiedMassFlux() {
        name_store("boundaryMassFluxNotCorrected",boundaryMassFluxNotCorrected) ;
        output("boundaryMassFluxNotCorrected") ;
        constraint("ref->massFlux_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<BoundaryMassFluxNotCorrectedSpecifiedMassFlux>
    registerBoundaryMassFluxNotCorrectedSpecifiedMassFlux ;

  // Noslip boundaries do not have mass flux corrected.
  class BoundaryMassFluxNotCorrectedNoslip : public pointwise_rule {
    private:
      store<bool> boundaryMassFluxNotCorrected ;
    public:

      // Define input and output.
      BoundaryMassFluxNotCorrectedNoslip() {
        name_store("boundaryMassFluxNotCorrected",boundaryMassFluxNotCorrected) ;
        output("boundaryMassFluxNotCorrected") ;
        constraint("noslip_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<BoundaryMassFluxNotCorrectedNoslip>
    registerBoundaryMassFluxNotCorrectedNoslip ;

  // Slip boundaries do not have mass flux corrected.
  class BoundaryMassFluxNotCorrectedSlip : public pointwise_rule {
    private:
      store<bool> boundaryMassFluxNotCorrected ;
    public:

      // Define input and output.
      BoundaryMassFluxNotCorrectedSlip() {
        name_store("boundaryMassFluxNotCorrected",boundaryMassFluxNotCorrected) ;
        output("boundaryMassFluxNotCorrected") ;
        constraint("slip_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<BoundaryMassFluxNotCorrectedSlip>
    registerBoundaryMassFluxNotCorrectedSlip ;

  // Symmetry boundaries do not have mass flux corrected.
  class BoundaryMassFluxNotCorrectedSymmetry : public pointwise_rule {
    private:
      store<bool> boundaryMassFluxNotCorrected ;
    public:

      // Define input and output.
      BoundaryMassFluxNotCorrectedSymmetry() {
        name_store("boundaryMassFluxNotCorrected",boundaryMassFluxNotCorrected) ;
        output("boundaryMassFluxNotCorrected") ;
        constraint("symmetry_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<BoundaryMassFluxNotCorrectedSymmetry>
    registerBoundaryMassFluxNotCorrectedSymmetry ;

  // Extrapolated pressure outlet boundaries do not have mass flux corrected.
  class BoundaryMassFluxNotCorrectedExtrapolatedPressureOutlet : public pointwise_rule {
    private:
      store<bool> boundaryMassFluxNotCorrected ;
    public:

      // Define input and output.
      BoundaryMassFluxNotCorrectedExtrapolatedPressureOutlet() {
        name_store("boundaryMassFluxNotCorrected",boundaryMassFluxNotCorrected) ;
        output("boundaryMassFluxNotCorrected") ;
        constraint("extrapolatedPressureOutlet_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<BoundaryMassFluxNotCorrectedExtrapolatedPressureOutlet>
    registerBoundaryMassFluxNotCorrectedExtrapolatedPressureOutlet ;

  // Constraint for all faces where mass flux is corrected. This adds interior
  // faces.
  class MassFluxCorrectedInterior : public pointwise_rule {
    private:
      store<bool> massFluxCorrected ;
    public:

      // Define input and output.
      MassFluxCorrectedInterior() {
        name_store("massFluxCorrected",massFluxCorrected) ;
        output("massFluxCorrected") ;
        constraint("internalFaces") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { massFluxCorrected[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<MassFluxCorrectedInterior>
    registerMassFluxCorrectedInterior ;

  // Constraint for all faces where mass flux is corrected. This adds
  // boundary faces.
  class MassFluxCorrectedBoundary : public pointwise_rule {
    private:
      store<bool> massFluxCorrected ;
    public:

      // Define input and output.
      MassFluxCorrectedBoundary() {
        name_store("massFluxCorrected",massFluxCorrected) ;
        output("massFluxCorrected") ;
        constraint("boundaryMassFluxCorrected") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { massFluxCorrected[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<MassFluxCorrectedBoundary>
    registerMassFluxCorrectedBoundary ;

  // Rule to compute the diffusion product for interior faces. Checked.
  class DiffusionProductInterior : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<vect3d> cellCenter ;
      const_store<Area> area ;
      store<real> diffusionProduct ;
    public:

      // Define input and output.
      DiffusionProductInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("cellcenter",cellCenter) ;
        name_store("area",area) ;
        name_store("diffusionProduct",diffusionProduct) ;
        input("(cl,cr)->cellcenter,area") ;
        output("diffusionProduct") ;
        constraint("internalFaces") ;
      }

      // Compute the diffusion product for a single internal face.
      void calculate(Entity face) {
        diffusionProduct[face]=area[face].sada/dot(area[face].n,
          cellCenter[cr[face]]-cellCenter[cl[face]]) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusionProductInterior> registerDiffusionProductInterior ;

  // Rule to compute the diffusion product for boundary faces. Checked.
  class DiffusionProductBoundary : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<Area> area ;
      store<real> diffusionProduct ;
    public:

      // Define input and output.
      DiffusionProductBoundary() {
        name_store("ci",ci) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("area",area) ;
        name_store("diffusionProduct",diffusionProduct) ;
        input("ci->cellcenter,facecenter,area") ;
        output("diffusionProduct") ;
        constraint("boundaryFaces") ;
      }

      // Compute the diffusion product for a single internal face.
      void calculate(Entity face) {
        diffusionProduct[face]=area[face].sada/dot(area[face].n,
          faceCenter[face]-cellCenter[ci[face]]) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusionProductBoundary> registerDiffusionProductBoundary ;

  // Time build rule for mass flux on interior faces.
  class TimeBuildInteriorMassFlux : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<real> rho_ic ;
      const_store<vect3d> v_ic ;
      const_store<Area> area_ic ;
      const_store<real> faceRadius ;
      store<real> massFlux ;
    public:

      // Define input and output.
      TimeBuildInteriorMassFlux() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("rho_ic",rho_ic) ;
        name_store("v_ic",v_ic) ;
        name_store("area_ic",area_ic) ;
        name_store("faceRadius",faceRadius) ;
        name_store("massFlux{n=0}",massFlux) ;
        input("(cl,cr)->(rho_ic,v_ic),area_ic,faceRadius") ;
        output("massFlux{n=0}") ;
        constraint("internalFaces,noRestart") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        massFlux[face]=0.25*(rho_ic[cl[face]]+rho_ic[cr[face]])*
          dot(v_ic[cl[face]]+v_ic[cr[face]],area_ic[face].n)*
          area_ic[face].sada*faceRadius[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<TimeBuildInteriorMassFlux> registerTimeBuildInteriorMassFlux ;

  // Time build rule for mass flux on boundary faces where mass flux is
  // iterating. This currently includes boundaries with a) specified pressure
  // where mass flux is not specified and b) total pressure inlets.
  class TimeBuildBoundaryMassFlux : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> rho_ic ;
      const_store<vect3d> v_ic ;
      const_store<Area> area_ic ;
      const_store<real> faceRadius ;
      store<real> massFlux ;
    public:

      // Define input and output.
      TimeBuildBoundaryMassFlux() {
        name_store("ci",ci) ;
        name_store("rho_ic",rho_ic) ;
        name_store("v_ic",v_ic) ;
        name_store("area_ic",area_ic) ;
        name_store("faceRadius",faceRadius) ;
        name_store("massFlux{n=0}",massFlux) ;
        input("ci->(rho_ic,v_ic),area_ic,faceRadius") ;
        output("massFlux{n=0}") ;
        constraint("boundaryMassFluxCorrected,noRestart") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        massFlux[face]=rho_ic[ci[face]]*dot(v_ic[ci[face]],
          area_ic[face].n)*area_ic[face].sada*faceRadius[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<TimeBuildBoundaryMassFlux> registerTimeBuildBoundaryMassFlux ;

  // Mass flux for boundaries where velocity is specified.
  class BoundaryMassFluxSpecifiedVelocity : public pointwise_rule {
    private:
      const_store<real> rho_f ;
      const_store<vect3d> v_f ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<real> massFlux ;
    public:

      // Define input and output.
      BoundaryMassFluxSpecifiedVelocity() {
        name_store("rho_f",rho_f) ;
        name_store("v_f",v_f) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("massFlux",massFlux) ;
        input("rho_f,v_f,area,faceRadius") ;
        output("massFlux") ;
        constraint("ref->v_BCoption") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        massFlux[face]=rho_f[face]*dot(v_f[face],area[face].n)*
          area[face].sada*faceRadius[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryMassFluxSpecifiedVelocity>
    registerBoundaryMassFluxSpecifiedVelocity ;

  // Mass flux for boundaries where mass flux is specified.
  class BoundaryMassFluxSpecifiedMassFlux : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<real> massFlux_BC ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<real> massFlux ;
    public:

      // Define input and output.
      BoundaryMassFluxSpecifiedMassFlux() {
        name_store("ref",ref) ;
        name_store("massFlux_BC",massFlux_BC) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("massFlux",massFlux) ;
        input("ref->massFlux_BC,area,faceRadius") ;
        output("massFlux") ;
        constraint("ref->massFlux_BC") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        massFlux[face]=massFlux_BC[ref[face]]*area[face].sada*
          faceRadius[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryMassFluxSpecifiedMassFlux>
    registerBoundaryMassFluxSpecifiedMassFlux ;

  // Mass flux for extrapolated pressure outlets.
  class BoundaryMassFluxExtrapolatedPressureOutlet : public pointwise_rule {
    private:
      const_store<real> rho_f ;
      const_store<vect3d> v_f ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<real> massFlux ;
    public:

      // Define input and output.
      BoundaryMassFluxExtrapolatedPressureOutlet() {
        name_store("rho_f",rho_f) ;
        name_store("v_f",v_f) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("massFlux",massFlux) ;
        input("rho_f,v_f,area,faceRadius") ;
        output("massFlux") ;
        constraint("extrapolatedPressureOutlet_BC") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        massFlux[face]=rho_f[face]*dot(v_f[face],area[face].n)*
          area[face].sada*faceRadius[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryMassFluxExtrapolatedPressureOutlet>
    registerBoundaryMassFluxExtrapolatedPressureOutlet ;

  // Mass flux for no-slip boundaries. Setting to the value of the grid mass
  // flux gives identical cancellation.
  class BoundaryMassFluxNoSlip : public pointwise_rule {
    private:
      const_store<real> gridMassFlux ;
      store<real> massFlux ;
    public:

      // Define input and output.
      BoundaryMassFluxNoSlip() {
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("massFlux",massFlux) ;
        input("gridMassFlux") ;
        output("massFlux") ;
        constraint("noslip_BC") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) { massFlux[face]=gridMassFlux[face] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryMassFluxNoSlip> registerBoundaryMassFluxNoSlip ;

  // Mass flux for no-slip boundaries. Setting to the value of the grid mass
  // flux gives identical cancellation.
  class BoundaryMassFluxSlip : public pointwise_rule {
    private:
      const_store<real> gridMassFlux ;
      store<real> massFlux ;
    public:

      // Define input and output.
      BoundaryMassFluxSlip() {
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("massFlux",massFlux) ;
        input("gridMassFlux") ;
        output("massFlux") ;
        constraint("slip_BC") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) { massFlux[face]=gridMassFlux[face] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryMassFluxSlip> registerBoundaryMassFluxSlip ;

  // Mass flux for symmetry boundaries. Setting to the value of the grid mass
  // flux gives identical cancellation.
  class BoundaryMassFluxSymmetry : public pointwise_rule {
    private:
      const_store<real> gridMassFlux ;
      store<real> massFlux ;
    public:

      // Define input and output.
      BoundaryMassFluxSymmetry() {
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("massFlux",massFlux) ;
        input("gridMassFlux") ;
        output("massFlux") ;
        constraint("symmetry_BC") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) { massFlux[face]=gridMassFlux[face] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryMassFluxSymmetry> registerBoundaryMassFluxSymmetry ;

  // TEMPORARY.
  class TempMassFluxCorrectedPredictor : public pointwise_rule {
    private:
      const_store<real> massFluxCorrected_p_temp ;
      store<real> massFluxCorrected_p ;
    public:

      // Define input and output.
      TempMassFluxCorrectedPredictor() {
        name_store("massFluxCorrected_p_temp",massFluxCorrected_p_temp) ;
        name_store("massFluxCorrected_p",massFluxCorrected_p) ;
        input("massFluxCorrected_p_temp") ;
        output("massFluxCorrected_p") ;
        constraint("massFluxCorrected") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        massFluxCorrected_p[face]=massFluxCorrected_p_temp[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<TempMassFluxCorrectedPredictor>
    registerTempMassFluxCorrectedPredictor ;

  // Passes massFlux through to massFluxCorrected_p so we can use one variable
  // in subsequent rules.
  class BoundaryMassFluxCorrectedPredictor : public pointwise_rule {
    private:
      const_store<real> massFlux ;
      store<real> massFluxCorrected_p ;
    public:

      // Define input and output.
      BoundaryMassFluxCorrectedPredictor() {
        name_store("massFlux",massFlux) ;
        name_store("massFluxCorrected_p",massFluxCorrected_p) ;
        input("massFlux") ;
        output("massFluxCorrected_p") ;
        constraint("boundaryMassFluxNotCorrected") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) { massFluxCorrected_p[face]=massFlux[face] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; } 

  } ;

  register_rule<BoundaryMassFluxCorrectedPredictor>
    registerBoundaryMassFluxCorrectedPredictor ;

  // TEMPORARY.
  class TempMassFluxCorrectedCorrector : public pointwise_rule {
    private:
      const_store<real> massFluxCorrected_c_temp ;
      store<real> massFluxCorrected_c ;
    public:

      // Define input and output.
      TempMassFluxCorrectedCorrector() {
        name_store("massFluxCorrected_c_temp",massFluxCorrected_c_temp) ;
        name_store("massFluxCorrected_c",massFluxCorrected_c) ;
        input("massFluxCorrected_c_temp") ;
        output("massFluxCorrected_c") ;
        constraint("massFluxCorrected") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        massFluxCorrected_c[face]=massFluxCorrected_c_temp[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<TempMassFluxCorrectedCorrector>
    registerTempMassFluxCorrectedCorrector ;

  // Passes massFlux through to massFluxCorrected_c so we can use one variable
  // in subsequent rules.
  class BoundaryMassFluxCorrectedCorrector : public pointwise_rule {
    private:
      const_store<real> massFlux ;
      store<real> massFluxCorrected_c ;
    public:

      // Define input and output.
      BoundaryMassFluxCorrectedCorrector() {
        name_store("massFlux",massFlux) ;
        name_store("massFluxCorrected_c",massFluxCorrected_c) ;
        input("massFlux") ;
        output("massFluxCorrected_c") ;
        constraint("boundaryMassFluxNotCorrected") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) { massFluxCorrected_c[face]=massFlux[face] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryMassFluxCorrectedCorrector>
    registerBoundaryMassFluxCorrectedCorrector ;

  // Iteration build rule for mass flux.
  class IterationBuildMassFlux : public pointwise_rule {
    private:
      const_store<real> massFluxN ;
      store<real> massFlux ;
    public:

      // Define input and output.
      IterationBuildMassFlux() {
        name_store("massFluxCorrected_p{n}",massFluxN) ;
        name_store("massFlux{n,it=0}",massFlux) ;
        input("massFluxCorrected_p{n}") ;
        output("massFlux{n,it=0}") ;
        constraint("massFluxCorrected{n}") ;
      }

      // Assign mass flux at iteration zero for a single face.
      void calculate(Entity face) { massFlux[face]=massFluxN[face] ; }

      // Loop over faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationBuildMassFlux> registerIterationBuildMassFlux ;

  // Rule to compute the density on interior faces for incompressible flow.
  // Note that we already have the values for the boundary faces (rho_f) from
  // boundary condition rules for density. Checked.
  class InteriorFaceDensityIncompressible : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<real> rho ;
      store<real> faceDensity ;
    public:

      // Define input and output.
      InteriorFaceDensityIncompressible() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("rho",rho) ;
        name_store("faceDensity",faceDensity) ;
        input("(cl,cr)->rho") ;
        output("faceDensity") ;
        constraint("internalFaces,incompressibleFlow") ;
      }

      // Calculate density for a single face.
      void calculate(Entity face) {
        faceDensity[face]=0.5*(rho[cl[face]]+rho[cr[face]]) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<InteriorFaceDensityIncompressible>
    registerInteriorFaceDensityIncompressible ;

  // Rule to compute the face mach number. Note that we have now changed this
  // to be the Mach number based on the face normal velocity. Thus we do not
  // need to change this for turbomachinery flows. JW 7/12/2007
  class FaceMachNumber : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<real> soundSpeed ;
      const_store<vect3d> v ;
      const_store<Area> area ;
      store<real> faceMachNumber ;
    public:

      // Define input and output.
      FaceMachNumber() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("v",v) ;
        name_store("area",area) ;
        name_store("soundSpeed",soundSpeed) ;
        name_store("faceMachNumber",faceMachNumber) ;
        input("(cl,cr)->(v,soundSpeed),area") ;
        output("faceMachNumber") ;
        constraint("internalFaces,compressibleFlow") ;
      } ;

      // Calculate mach number for a single face.
      void calculate(Entity face) {
        faceMachNumber[face]=abs(dot(v[cl[face]]+v[cr[face]],area[face].n))/
          (soundSpeed[cl[face]]+soundSpeed[cr[face]]) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<FaceMachNumber> registerFaceMachNumber ;

  // Rule to compute the density on interior faces for compressible flow and
  // first-order upwinding. Note that we already have the values for rho_f for
  // the boundary faces from boundary condition rules for density. No change
  // required for multiple reference frames since we are dotting the face
  // velocity with the face normal, which is assumed always to be paralled to
  // the axis of rotation.
  class InteriorFaceDensityCompressibleFOU : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<real> rho ;
      const_store<vect3d> v ;
      const_store<real> faceMachNumber ;
      const_store<Area> area ;
      store<real> faceDensity ;
    public:

      // Define input and output.
      InteriorFaceDensityCompressibleFOU() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("rho",rho) ;
        name_store("v",v) ;
        name_store("faceMachNumber",faceMachNumber) ;
        name_store("area",area) ;
        name_store("faceDensity",faceDensity) ;
        input("(cl,cr)->(rho,v)") ;
        input("faceMachNumber,area") ;
        output("faceDensity") ;
        constraint("internalFaces,compressibleFlow,fouInviscidFlux") ;
      }

      // Calculate face density for a single face.
      void calculate(Entity face) {
        if(faceMachNumber[face]<0.3){
          faceDensity[face]=0.5*(rho[cl[face]]+rho[cr[face]]) ;
        }else{
          vect3d faceVelocity=0.5*(v[cl[face]]+v[cr[face]]) ;
          faceDensity[face]=(dot(faceVelocity,area[face].n)>0.0)? rho[cl[face]]:
            rho[cr[face]] ;
        }
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<InteriorFaceDensityCompressibleFOU>
    registerInteriorFaceDensityCompressibleFOU ;

  // Rule to compute the density on interior faces for compressible flow and
  // either second-order upwind or the Roe scheme. Note that we already have
  // the values for rho_f for the boundary faces from boundary condition rules
  // for density. No change required for turbomachinery.
  class InteriorFaceDensityCompressibleSOUROE : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<real> rho ;
      const_store<vect3d> rhoGradient ;
      const_store<real> rhoLimiter ;
      const_store<vect3d> v ;
      const_store<vect3d> cellCenter ;
      const_store<Area> area ;
      const_store<vect3d> faceCenter ;
      const_store<real> faceMachNumber ;
      store<real> faceDensity ;
    public:

      // Define input and output.
      InteriorFaceDensityCompressibleSOUROE() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("rho",rho) ;
        name_store("grads(rho)",rhoGradient) ;
        name_store("limiters(rho)",rhoLimiter) ;
        name_store("v",v) ;
        name_store("cellcenter",cellCenter) ;
        name_store("area",area) ;
        name_store("facecenter",faceCenter) ;
        name_store("faceMachNumber",faceMachNumber) ;
        name_store("faceDensity",faceDensity) ;
        input("(cl,cr)->(rho,grads(rho),limiters(rho),v,cellcenter)") ;
        input("area,facecenter,faceMachNumber") ;
        output("faceDensity") ;
        constraint("internalFaces,compressibleFlow,souOrRoeInviscidFlux") ;
      }

      // Calculate face density for a single face.
      void calculate(Entity face) {
        if(faceMachNumber[face]<0.3){
          faceDensity[face]=0.5*(rho[cl[face]]+rho[cr[face]]) ;
        }else{
          vect3d faceVelocity=0.5*(v[cl[face]]+v[cr[face]]) ;
          faceDensity[face]=(dot(faceVelocity,area[face].n)>0.0)? rho[cl[face]]+
            rhoLimiter[cl[face]]*dot(rhoGradient[cl[face]],faceCenter[face]-
            cellCenter[cl[face]]):rho[cr[face]]+rhoLimiter[cr[face]]*
            dot(rhoGradient[cr[face]],faceCenter[face]-cellCenter[cr[face]]) ;
        }
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<InteriorFaceDensityCompressibleSOUROE>
    registerInteriorFaceDensityCompressibleSOUROE ;

  // Rule to compute the Mach number.
  class ComputeMachNumber : public pointwise_rule {
    private:
      const_store<vect3d> v ;
      const_store<real> soundSpeed ;
      store<real> machNumber ;
    public:

      // Define input and output.
      ComputeMachNumber() {
        name_store("X",v) ;
        name_store("soundSpeed",soundSpeed) ;
        name_store("machNumber(X)",machNumber) ;
        input("X,soundSpeed") ;
        output("machNumber(X)") ;
      }

      // Compute Mach number for a single entity.
      void calculate(Entity e) { machNumber[e]=norm(v[e])/soundSpeed[e] ; }

      // Loop over entities.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeMachNumber> registerComputeMachNumber ;

  // Rule to compute the new mass flux on interior faces for the FOU and
  // SOU convection schemes. New form of momentum interpolation added by
  // ST to prevent decoupling observed with BDF2 on cylinder case.
  class MassFluxStarInteriorFOUSOU : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<vect3d> vStar ;
      const_store<real> p ;
      const_store<vect3d> pGradient ;
      const_store<vect3d> cellCenter ;
      const_store<real> vol ;
      const_store<real> faceDensity ;
      const_store<Area> area ;
      const_store<real> diffusionProduct ;
      const_store<real> pPrimeCoefficient ;
      const_store<real> faceRadius ;
      store<real> massFluxStar ;
    public:

      // Define input and output.
      MassFluxStarInteriorFOUSOU() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("vStar",vStar) ;
        name_store("p",p) ;
        name_store("grads(p)",pGradient) ;
        name_store("cellcenter",cellCenter) ;
        name_store("vol",vol) ;
        name_store("faceDensity",faceDensity) ;
        name_store("area",area) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        name_store("faceRadius",faceRadius) ;
        name_store("massFluxStar",massFluxStar) ;
        name_store("diffusionProduct",diffusionProduct) ;
        input("(cl,cr)->(vStar,p,grads(p),cellcenter,vol)") ;
        input("faceDensity,diffusionProduct,pPrimeCoefficient,faceRadius,area") ;
        output("massFluxStar") ;
        constraint("internalFaces,fouOrSouInviscidFlux") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        vect3d averagePressureGradient=(pGradient[cl[face]]*vol[cr[face]]+
          pGradient[cr[face]]*vol[cl[face]])/(vol[cl[face]]+vol[cr[face]]) ;
//      massFluxStar[face]=0.5*faceDensity[face]*dot(vStar[cl[face]]+
//        vStar[cr[face]],area[face].n)*area[face].sada*faceRadius[face]-
//        pPrimeCoefficient[face]*(p[cr[face]]-p[cl[face]]-
//        dot(averagePressureGradient,area[face].n)*area[face].sada/
//        diffusionProduct[face]) ;
        massFluxStar[face]=0.5*faceDensity[face]*dot(vStar[cl[face]]+
          vStar[cr[face]],area[face].n)*area[face].sada*faceRadius[face]-
          pPrimeCoefficient[face]*(p[cr[face]]-p[cl[face]]-
          dot(averagePressureGradient,cellCenter[cr[face]]-
          cellCenter[cl[face]])) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<MassFluxStarInteriorFOUSOU> registerMassFluxStarInteriorFOUSOU ;

  // Rule to compute the new mass flux on boundary faces. Only need this for
  // boundary faces where the mass flux is going to be corrected.
  class MassFluxStarBoundary : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<vect3d> vStar ;
      const_store<real> rho_f ;
      const_store<Area> area ;
      store<real> massFluxStar ;
    public:

      // Define input and output.
      MassFluxStarBoundary() {
        name_store("ci",ci) ;
        name_store("vStar",vStar) ;
        name_store("rho_f",rho_f) ;
        name_store("area",area) ;
        name_store("massFluxStar",massFluxStar) ;
        input("ci->vStar,rho_f,area") ;
        output("massFluxStar") ;
        constraint("boundaryMassFluxCorrected") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        massFluxStar[face]=rho_f[face]*dot(vStar[ci[face]],
          area[face].n)*area[face].sada ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<MassFluxStarBoundary> registerMassFluxStarBoundary ;

  // Rule to compute the new mass flux in the corrector stage on interior
  // faces for the FOU and SOU convection schemes.
  class MassFluxStarHatInteriorFOUSOU : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<vect3d> vHat;
      const_store<real> vMainCoefficient_it,vMainCoefficient ;
      const_store<real> vol ;
      const_store<real> faceDensity,massFlux ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<real> massFluxStarHat ;
    public:

      // Define input and output.
      MassFluxStarHatInteriorFOUSOU() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("vHat",vHat) ;
        name_store("vMainCoefficient_it",vMainCoefficient_it) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        name_store("vol",vol) ;
        name_store("faceDensity",faceDensity) ;
        name_store("massFlux",massFlux) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("massFluxStarHat",massFluxStarHat) ;
        input("(cl,cr)->(vHat,vMainCoefficient_it,vMainCoefficient,vol)") ;
        input("faceDensity,massFlux,faceRadius,area") ;
        output("massFluxStarHat") ;
        constraint("internalFaces,fouOrSouInviscidFlux") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        real vL=vol[cl[face]],vR=vol[cr[face]] ;
        real vMainCoefficientRatio=(vR*(vMainCoefficient_it[cl[face]]/
          vMainCoefficient[cl[face]])+vL*(vMainCoefficient_it[cr[face]]/
          vMainCoefficient[cr[face]]))/(vL+vR) ;
        vect3d vHatAverage=(vR*vHat[cl[face]]+vL*vHat[cr[face]])/(vL+vR) ;
        massFluxStarHat[face]=vMainCoefficientRatio*massFlux[face]+
          faceDensity[face]*dot(vHatAverage,area[face].n)*area[face].sada ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<MassFluxStarHatInteriorFOUSOU>
    registerMassFluxStarHatInteriorFOUSOU ;

  // Rule to compute the new mass flux on boundary faces. Only need this for
  // boundary faces where the mass flux is going to be corrected.
  class MassFluxStarHatBoundary : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<vect3d> vHat;
      const_store<real> vMainCoefficient_it,vMainCoefficient ;
      const_store<real> rho_f,massFlux ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<real> massFluxStarHat ;
    public:

      // Define input and output.
      MassFluxStarHatBoundary() {
        name_store("ci",ci) ;
        name_store("vHat",vHat) ;
        name_store("vMainCoefficient_it",vMainCoefficient_it) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        name_store("rho_f",rho_f) ;
        name_store("massFlux",massFlux) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("massFluxStarHat",massFluxStarHat) ;
        input("ci->(vHat,vMainCoefficient_it,vMainCoefficient)") ;
        input("rho_f,massFlux,faceRadius,area") ;
        output("massFluxStarHat") ;
        constraint("boundaryMassFluxCorrected") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        massFluxStarHat[face]=vMainCoefficient_it[ci[face]]/
          vMainCoefficient[ci[face]]*massFlux[face]+rho_f[face]*
          dot(vHat[ci[face]],area[face].n)*area[face].sada ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<MassFluxStarHatBoundary>
    registerMassFluxStarHatBoundary ;

  // Iteration advance rule for mass flux.
  class IterationAdvanceMassFlux : public pointwise_rule {
    private:
      const_store<real> massFluxCorrected ;
      store<real> massFlux ;
    public:

      // Define input and output.
      IterationAdvanceMassFlux() {
        name_store("massFluxCorrected_c{n,it}",massFluxCorrected) ;
        name_store("massFlux{n,it+1}",massFlux) ;
        input("massFluxCorrected_c{n,it}") ;
        output("massFlux{n,it+1}") ;
        constraint("massFluxCorrected{n,it}") ;
      }

      // Assign mass flux at end of iteration for a single face.
      void calculate(Entity face) {
        massFlux[face]=massFluxCorrected[face] ;
      }

      // Loop over faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationAdvanceMassFlux> registerIterationAdvanceMassFlux ;

  // Iteration collapse rule for mass flux.
  class IterationCollapseMassFlux : public pointwise_rule {
    private:
      store<real> massFlux ;
    public:

      // Define input and output.
      IterationCollapseMassFlux() {
        name_store("massFlux{n,it}",massFlux) ;
        input("massFlux{n,it}") ;
        output("massFlux{n+1}=massFlux{n,it}") ;
        conditional("iterationFinished{n,it-1}") ;
        constraint("massFluxCorrected{n,it}") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<IterationCollapseMassFlux> registerIterationCollapseMassFlux ;

  // Time collapse rule for face data. Old rule.
  class TimeCollapseFace : public pointwise_rule {
    private:
      const_store<real> massFlux ;
      store<int> solution ;
    public:

      // Define input and output.
      TimeCollapseFace() {
        name_store("massFlux{n}",massFlux) ;
        name_store("solution",solution) ;
        input("massFlux{n}") ;
        output("solution") ;
        conditional("timeSteppingFinished{n}") ;
        constraint("faces{n}") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}
  } ;

//register_rule<TimeCollapseFace> registerTimeCollapseFace ;

  // Unit rule to initialize net inlet mass flux.
  class InletMassFluxUnit : public unit_rule {
    private:
      param<real> inletMassFlux ;
    public:

      // Define input and output.
      InletMassFluxUnit() {
        name_store("inletMassFlux",inletMassFlux) ;
        output("inletMassFlux") ;
        constraint("inlet_BC") ;
      }

      // Initialize.
      virtual void compute(const sequence &seq) { *inletMassFlux=0.0 ; }
  } ;

  register_rule<InletMassFluxUnit> registerInletMassFluxUnit ;

  // Apply rule to compute the net mass flux through all inlets.
  class InletMassFluxApply : public apply_rule<param<real>,Loci::Summation
  <real> > {
    private:
      const_store<real> massFlux ;
      param<real> inletMassFlux ;
    public:

      // Define input and output.
      InletMassFluxApply() {
        name_store("massFlux",massFlux) ;
        name_store("inletMassFlux",inletMassFlux) ;
        input("massFlux") ;
        output("inletMassFlux") ;
        constraint("inlet_BC") ;
      }

      // Add the face increment.
      void calculate(Entity face) { *inletMassFlux+=massFlux[face] ; }

      // Sum over all inlet faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InletMassFluxApply> registerInletMassFluxApply ;

  // Write the inlet mass flux.
  class WriteInletMassFlux : public pointwise_rule {
    private:
      const_param<real> inletMassFlux ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      WriteInletMassFlux() {
        name_store("inletMassFlux{n,it}",inletMassFlux) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("inletMassFlux{n,it}") ;
        output("OUTPUT{n,it}") ;
        conditional("do_print{n,it}") ;
      }

      // Initialize.
      virtual void compute(const sequence &seq) {
        if(Loci::MPI_rank==0){
          cout.setf(ios::scientific,ios::floatfield) ; cout.precision(4) ;
          cout << "Inlet Mass Transfer = " << *inletMassFlux << " kg/s"
            << endl ;
        }
      }
  } ;

  //register_rule<WriteInletMassFlux> registerWriteInletMass ;

  // Unit rule to initialize net outlet mass flux.
  class OutletMassFluxUnit : public unit_rule {
    private:
      param<real> outletMassFlux ;
    public:

      // Define input and output.
      OutletMassFluxUnit() {
        name_store("outletMassFlux",outletMassFlux) ;
        output("outletMassFlux") ;
        constraint("outlet_BC") ;
      }

      // Initialize.
      virtual void compute(const sequence &seq) { *outletMassFlux=0.0 ; }
  } ;

  register_rule<OutletMassFluxUnit> registerOutletMassFluxUnit ;

  // Apply rule to compute the net mass flux through all outlets.
  class OutletMassFluxApply : public apply_rule<param<real>,Loci::Summation
  <real> > {
    private:
      const_store<real> massFlux ;
      param<real> outletMassFlux ;
    public:

      // Define input and output.
      OutletMassFluxApply() {
        name_store("massFlux",massFlux) ;
        name_store("outletMassFlux",outletMassFlux) ;
        input("massFlux") ;
        output("outletMassFlux") ;
        constraint("outlet_BC") ;
      }

      // Add the face increment.
      void calculate(Entity face) { *outletMassFlux+=massFlux[face] ; }

      // Sum over all outlet faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<OutletMassFluxApply> registerOutletMassFluxApply ;

  // Write the outlet mass flux.
  class WriteOutletMassFlux : public pointwise_rule {
    private:
      const_param<real> outletMassFlux ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      WriteOutletMassFlux() {
        name_store("outletMassFlux{n,it}",outletMassFlux) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("outletMassFlux{n,it}") ;
        output("OUTPUT{n,it}") ;
        conditional("do_print{n,it}") ;
      }

      // Initialize.
      virtual void compute(const sequence &seq) {
        if(Loci::MPI_rank==0){
          cout.setf(ios::scientific,ios::floatfield) ; cout.precision(4) ;
          cout << "Outlet Mass Transfer = " << *outletMassFlux << " kg/s"
            << endl ;
        }
      }
  } ;

  //register_rule<WriteOutletMassFlux> registerWriteOutletMass ;
}
