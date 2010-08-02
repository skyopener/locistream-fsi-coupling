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

  // Creates the SIMPLE constraint, which will activate certain features in the
  // gridMotion module.
  class SIMPLEConstraint : public constraint_rule {
    private:
      Constraint SIMPLE ;
    public:

      // Define input and output.
      SIMPLEConstraint() {
        name_store("SIMPLE",SIMPLE) ;
        output("SIMPLE") ;
        constraint("UNIVERSE") ;
      }

      // Set up the constraint.
      virtual void compute(const sequence& seq) { SIMPLE=~EMPTY ; }
  } ;

  register_rule<SIMPLEConstraint> registerSIMPLEConstraint ;

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

      // Compute the diffusion product for a sequence of faces.
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

      // Compute the diffusion product for a sequence of faces.
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
      store<real> massFluxTimeStepZero ;
    public:

      // Define input and output.
      TimeBuildInteriorMassFlux() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("rho_ic",rho_ic) ;
        name_store("v_ic",v_ic) ;
        name_store("area_ic",area_ic) ;
        name_store("faceRadius",faceRadius) ;
        name_store("massFlux{n=0}",massFluxTimeStepZero) ;
        input("(cl,cr)->(rho_ic,v_ic),area_ic,faceRadius") ;
        output("massFlux{n=0}") ;
        constraint("internalFaces,noRestart") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        massFluxTimeStepZero[face]=0.25*(rho_ic[cl[face]]+rho_ic[cr[face]])*
          dot(v_ic[cl[face]]+v_ic[cr[face]],area_ic[face].n)*
          area_ic[face].sada*faceRadius[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<TimeBuildInteriorMassFlux> registerTimeBuildInteriorMassFlux ;

  // Time build rule for mass flux on boundary faces. This rule provides values
  // for faces where boundary density and velocity values have not been provided
  // in the .vars file.
  class TimeBuildBoundaryMassFlux : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> rho_ic ;
      const_store<vect3d> v_ic ;
      const_store<Area> area_ic ;
      const_store<real> faceRadius ;
      store<real> massFluxTimeStepZero ;
    public:

      // Define input and output.
      TimeBuildBoundaryMassFlux() {
        name_store("ci",ci) ;
        name_store("rho_ic",rho_ic) ;
        name_store("v_ic",v_ic) ;
        name_store("area_ic",area_ic) ;
        name_store("faceRadius",faceRadius) ;
        name_store("massFlux{n=0}",massFluxTimeStepZero) ;
        input("ci->(rho_ic,v_ic),area_ic,faceRadius") ;
        output("massFlux{n=0}") ;
        constraint("boundaryFaces,noRestart") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        massFluxTimeStepZero[face]=rho_ic[ci[face]]*dot(v_ic[ci[face]],
          area_ic[face].n)*area_ic[face].sada*faceRadius[face] ;
      }

      // Calculate mass flux for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<TimeBuildBoundaryMassFlux> registerTimeBuildBoundaryMassFlux ;

  // Priority build rule for mass flux on incompressible inlets with
  // time-independent velocity boundary conditions. The constraint
  // "ref->v_BCoption" prevents conflict with another mass flux build
  // rule when we have specified "mdot" for the incompressible inlet.
  class TimeBuildBoundaryMassFluxIncompressibleInletSteady : public
  pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> rho_ic ;
      const_store<vect3d> v_f ;
      const_store<Area> area_ic ;
      const_store<real> faceRadius ;
      store<real> massFluxTimeStepZero ;
    public:

      // Define input and output.
      TimeBuildBoundaryMassFluxIncompressibleInletSteady() {
        name_store("ci",ci) ;
        name_store("rho_ic",rho_ic) ;
        name_store("v_f",v_f) ;
        name_store("area_ic",area_ic) ;
        name_store("faceRadius",faceRadius) ;
        name_store("incompressibleInlet::massFlux{n=0}",massFluxTimeStepZero) ;
        input("ci->rho_ic,v_f,area_ic,faceRadius") ;
        output("incompressibleInlet::massFlux{n=0}") ;
        constraint("incompressibleInlet_BC,noRestart") ;
        constraint("ref->v_BCoption,ref->timeIndependent_BC") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        massFluxTimeStepZero[face]=rho_ic[ci[face]]*dot(v_f[face],area_ic[face].n)*
          area_ic[face].sada*faceRadius[face] ;
      }

      // Calculate mass flux for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<TimeBuildBoundaryMassFluxIncompressibleInletSteady>
    registerTimeBuildBoundaryMassFluxIncompressibleInletSteady ;

  // Similar rule for time-dependent velocity boundary conditions. Here
  // we use a differnet variable v_f_ic since at stationary time we
  // have no "stime" to use in the rules for v_f.
  class TimeBuildBoundaryMassFluxIncompressibleInletUnsteady : public
  pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> rho_ic ;
      const_store<vect3d> v_f ;
      const_store<Area> area_ic ;
      const_store<real> faceRadius ;
      store<real> massFluxTimeStepZero ;
    public:

      // Define input and output.
      TimeBuildBoundaryMassFluxIncompressibleInletUnsteady() {
        name_store("ci",ci) ;
        name_store("rho_ic",rho_ic) ;
        name_store("v_f_ic",v_f) ;
        name_store("area_ic",area_ic) ;
        name_store("faceRadius",faceRadius) ;
        name_store("incompressibleInlet::massFlux{n=0}",massFluxTimeStepZero) ;
        input("ci->rho_ic,v_f_ic,area_ic,faceRadius") ;
        output("incompressibleInlet::massFlux{n=0}") ;
        constraint("incompressibleInlet_BC,noRestart") ;
        constraint("ref->v_BCoption,ref->timeDependent_BC") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        massFluxTimeStepZero[face]=rho_ic[ci[face]]*dot(v_f[face],area_ic[face].n)*
          area_ic[face].sada*faceRadius[face] ;
      }

      // Calculate mass flux for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<TimeBuildBoundaryMassFluxIncompressibleInletUnsteady>
    registerTimeBuildBoundaryMassFluxIncompressibleInletUnsteady ;

  // Time build rule for mass flux on boundary faces. This prioriry rule
  // overrides the standard rule for faces on boundaries where the mass flux
  // has been explicitly specified in the .vars file. Note: Uncommented
  // constraint on 01/31/05 when adding restart capability.
  class TimeBuildBoundaryMassFluxSpecified : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<real> massFlux_BC ;
      const_store<Area> area_ic ;
      const_store<real> faceRadius ;
      store<real> massFluxTimeStepZero ;
    public:

      // Define input and output.
      TimeBuildBoundaryMassFluxSpecified() {
        name_store("ref",ref) ;
        name_store("massFlux_BC",massFlux_BC) ;
        name_store("area_ic",area_ic) ;
        name_store("faceRadius",faceRadius) ;
        name_store("specified::massFlux{n=0}",massFluxTimeStepZero) ;
        input("ref->massFlux_BC,area_ic,faceRadius") ;
        output("specified::massFlux{n=0}") ;
        constraint("ref->massFlux_BC,noRestart") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        massFluxTimeStepZero[face]=massFlux_BC[ref[face]]*area_ic[face].sada*
          faceRadius[face] ;
      }

      // Calculate mass flux for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<TimeBuildBoundaryMassFluxSpecified>
    registerTimeBuildBoundaryMassFluxSpecified ;

  // Time build rule for mass flux on boundary faces. This prioriry rule
  // overrides the standard rule for no-slip faces.
  class TimeBuildBoundaryMassFluxNoSlip : public pointwise_rule {
    private:
      store<real> massFluxTimeStepZero ;
    public:

      // Define input and output.
      TimeBuildBoundaryMassFluxNoSlip() {
        name_store("noslip::massFlux{n=0}",massFluxTimeStepZero) ;
        output("noslip::massFlux{n=0}") ;
        constraint("noslip_BC,noRestart") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) { massFluxTimeStepZero[face]=0.0 ; }

      // Calculate mass flux for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<TimeBuildBoundaryMassFluxNoSlip>
    registerTimeBuildBoundaryMassFluxNoSlip ;

  // Time build rule for mass flux on boundary faces. This prioriry rule
  // overrides the standard rule for slip faces.
  class TimeBuildBoundaryMassFluxSlip : public pointwise_rule {
    private:
      store<real> massFluxTimeStepZero ;
    public:

      // Define input and output.
      TimeBuildBoundaryMassFluxSlip() {
        name_store("slip::massFlux{n=0}",massFluxTimeStepZero) ;
        output("slip::massFlux{n=0}") ;
        constraint("slip_BC,noRestart") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) { massFluxTimeStepZero[face]=0.0 ; }

      // Calculate mass flux for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<TimeBuildBoundaryMassFluxSlip>
    registerTimeBuildBoundaryMassFluxSlip ;

  // Time build rule for mass flux on boundary faces. This prioriry rule
  // overrides the standard rule for symmetry faces.
  class TimeBuildBoundaryMassFluxSymmetry : public pointwise_rule {
    private:
      store<real> massFluxTimeStepZero ;
    public:

      // Define input and output.
      TimeBuildBoundaryMassFluxSymmetry() {
        name_store("symmetry::massFlux{n=0}",massFluxTimeStepZero) ;
        output("symmetry::massFlux{n=0}") ;
        constraint("symmetry_BC,noRestart") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) { massFluxTimeStepZero[face]=0.0 ; }

      // Calculate mass flux for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<TimeBuildBoundaryMassFluxSymmetry>
    registerTimeBuildBoundaryMassFluxSymmetry ;

  // Iteration build rule for mass flux.
  class IterationBuildMassFlux : public pointwise_rule {
    private:
      const_store<real> massFluxTimeStepN ;
      store<real> massFluxIterationZero ;
    public:

      // Define input and output.
      IterationBuildMassFlux() {
        name_store("massFlux{n}",massFluxTimeStepN) ;
        name_store("massFlux{n,it=0}",massFluxIterationZero) ;
        input("massFlux{n}") ;
        output("massFlux{n,it=0}") ;
        constraint("massFlux{n}") ;
      }

      // Assign mass flux at iteration zero for a single face.
      inline void calculate(Entity face) {
        massFluxIterationZero[face]=massFluxTimeStepN[face] ;
      }

      // Assign mass flux at iteration zero for a sequence of faces.
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

      // Calculate face density for all faces in sequence.
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

      // Compute density for a sequence of entities.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeMachNumber> registerComputeMachNumber ;

  // Rule to compute the stage 0 mass flux on interior faces for the FOU and
  // SOU convection schemes. No change required for turbomachinery.
  class StageZeroMassFluxInteriorFOUSOU : public pointwise_rule {
    private:
      const_param<real> vRelaxationFactor ;
      const_Map cl,cr ;
      const_store<vect3d> v ;
      const_store<real> faceDensity ;
      const_store<real> massFlux ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<real> stageZeroMassFlux ;
    public:

      // Define input and output.
      StageZeroMassFluxInteriorFOUSOU() {
        name_store("vRelaxationFactor",vRelaxationFactor) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("v",v) ;
        name_store("faceDensity",faceDensity) ;
        name_store("massFlux",massFlux) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("stageZeroMassFlux",stageZeroMassFlux) ;
        input("vRelaxationFactor") ;
        input("(cl,cr)->v") ;
        input("faceDensity,massFlux,area,faceRadius") ;
        output("stageZeroMassFlux") ;
        constraint("internalFaces") ;
      }

      // Calculate the mass flux for a single face.
      void calculate(Entity face) {
        stageZeroMassFlux[face]=(1.0-(*vRelaxationFactor))*(massFlux[face]-
          0.5*faceDensity[face]*dot(v[cl[face]]+v[cr[face]],area[face].n)*
          area[face].sada*faceRadius[face]) ;
      }

      // Calculate the mass flux for all faces in a sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<StageZeroMassFluxInteriorFOUSOU>
    registerStageZeroMassFluxInteriorFOUSOU ;

  // Rule to compute the stage 1 mass flux on interior faces for the FOU and
  // SOU convection schemes. No change required for turbomachninery.
  class StageOneMassFluxInteriorFOUSOU : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<vect3d> vStar ;
      const_store<real> p ;
      const_store<vect3d> pGradient ;
      const_store<vect3d> cellCenter ;
      const_store<real> vol ;
      const_store<real> faceDensity ;
      const_store<real> stageZeroMassFlux ;
      const_store<Area> area ;
      const_store<real> pPrimeCoefficient ;
      const_store<real> faceRadius ;
      store<real> stageOneMassFlux ;
    public:

      // Define input and output.
      StageOneMassFluxInteriorFOUSOU() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("vStar",vStar) ;
        name_store("p",p) ;
        name_store("grads(p)",pGradient) ;
        name_store("cellcenter",cellCenter) ;
        name_store("vol",vol) ;
        name_store("faceDensity",faceDensity) ;
        name_store("stageZeroMassFlux",stageZeroMassFlux) ;
        name_store("area",area) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        name_store("faceRadius",faceRadius) ;
        name_store("stageOneMassFlux",stageOneMassFlux) ;
        input("(cl,cr)->(vStar,p,grads(p),cellcenter,vol)") ;
        input("faceDensity,stageZeroMassFlux,pPrimeCoefficient") ;
        input("faceRadius,area") ;
        output("stageOneMassFlux") ;
        constraint("internalFaces,fouOrSouInviscidFlux") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {

        // Compute the averaged pressure gradient.
        vect3d averagePressureGradient=(pGradient[cl[face]]*vol[cr[face]]+
          pGradient[cr[face]]*vol[cl[face]])/(vol[cl[face]]+vol[cr[face]]) ;

        // Compute the stage one mass flux.
        stageOneMassFlux[face]=0.5*faceDensity[face]*dot(vStar[cl[face]]+
          vStar[cr[face]],area[face].n)*area[face].sada*faceRadius[face]-
          pPrimeCoefficient[face]*(p[cr[face]]-p[cl[face]]-
          dot(averagePressureGradient,cellCenter[cr[face]]-
          cellCenter[cl[face]]))+stageZeroMassFlux[face] ;
      }

      // Calculate mass flux for all faces in a sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<StageOneMassFluxInteriorFOUSOU>
    registerStageOneMassFluxInteriorFOUSOU ;

  // Rule to compute the stage 1 mass flux on interior faces for the ROE scheme.
  // No change required for turbomachinery.
  class StageOneMassFluxInteriorRoe : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<vect3d> vStar ;
      const_store<real> p ;
      const_store<vect3d> pGradient ;
      const_store<vect3d> cellCenter ;
      const_store<real> vol ;
      const_store<real> faceDensity ;
      const_store<real> faceMachNumber ;
      const_store<real> stageZeroMassFlux ;
      const_store<real> roeMassFlux ;
      const_store<Area> area ;
      const_store<real> pPrimeCoefficient ;
      store<real> stageOneMassFlux ;
    public:

      // Define input and output.
      StageOneMassFluxInteriorRoe() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("vStar",vStar) ;
        name_store("p",p) ;
        name_store("grads(p)",pGradient) ;
        name_store("cellcenter",cellCenter) ;
        name_store("vol",vol) ;
        name_store("faceDensity",faceDensity) ;
        name_store("faceMachNumber",faceMachNumber) ;
        name_store("stageZeroMassFlux",stageZeroMassFlux) ;
        name_store("roeMassFlux",roeMassFlux) ;
        name_store("area",area) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        name_store("stageOneMassFlux",stageOneMassFlux) ;
        input("(cl,cr)->(vStar,p,grads(p),cellcenter,vol)") ;
        input("faceDensity,faceMachNumber,stageZeroMassFlux,roeMassFlux") ;
        input("pPrimeCoefficient,area") ;
        output("stageOneMassFlux") ;
        constraint("internalFaces,roeInviscidFlux") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        if(faceMachNumber[face]>0.3){
          stageOneMassFlux[face]=roeMassFlux[face] ;
        }else{
          vect3d averagePressureGradient=(pGradient[cl[face]]*vol[cr[face]]+
            pGradient[cr[face]]*vol[cl[face]])/(vol[cl[face]]+vol[cr[face]]) ;
          stageOneMassFlux[face]=0.5*faceDensity[face]*dot(vStar[cl[face]]+
            vStar[cr[face]],area[face].n)*area[face].sada-
            pPrimeCoefficient[face]*(p[cr[face]]-p[cl[face]]-
            dot(averagePressureGradient,cellCenter[cr[face]]-
            cellCenter[cl[face]]))+stageZeroMassFlux[face] ;
        }
      }

      // Calculate mass flux for all faces in a sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<StageOneMassFluxInteriorRoe>
    registerStageOneMassFluxInteriorRoe ;

  // Rule to compute the stage 1 mass flux on boundary faces for the FOU and
  // SOU convection schemes.
  class StageOneMassFluxBoundaryFOUSOU : public pointwise_rule {
    private:
      const_store<real> rho_f ;
      const_store<vect3d> v_f ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<real> stageOneMassFlux ;
    public:

      // Define input and output.
      StageOneMassFluxBoundaryFOUSOU() {
        name_store("rho_f",rho_f) ;
        name_store("v_f",v_f) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("stageOneMassFlux",stageOneMassFlux) ;
        input("rho_f,v_f,area,faceRadius") ;
        output("stageOneMassFlux") ;
        constraint("boundaryFaces") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        stageOneMassFlux[face]=rho_f[face]*dot(v_f[face],area[face].n)*
          area[face].sada*faceRadius[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<StageOneMassFluxBoundaryFOUSOU>
    registerStageOneMassFluxBoundaryFOUSOU ;

  // Rule to compute the stage 1 mass flux for noslip faces. Set to
  // the grid mass flux so there is identical cancellation.
  class StageOneMassFluxBoundaryNoslip : public pointwise_rule {
    private:
      const_store<real> gridMassFlux ;
      store<real> stageOneMassFlux ;
    public:

      // Define input and output.
      StageOneMassFluxBoundaryNoslip() {
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("noslip::stageOneMassFlux",stageOneMassFlux) ;
        input("gridMassFlux") ;
        output("noslip::stageOneMassFlux") ;
        constraint("noslip_BC") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        stageOneMassFlux[face]=gridMassFlux[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<StageOneMassFluxBoundaryNoslip>
    registerStageOneMassFluxBoundaryNoslip ;

  // Rule to compute the stage 1 mass flux for slip faces. Set to
  // the grid mass flux so there is identical cancellation.
  class StageOneMassFluxBoundarySlip : public pointwise_rule {
    private:
      const_store<real> gridMassFlux ;
      store<real> stageOneMassFlux ;
    public:

      // Define input and output.
      StageOneMassFluxBoundarySlip() {
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("slip::stageOneMassFlux",stageOneMassFlux) ;
        input("gridMassFlux") ;
        output("slip::stageOneMassFlux") ;
        constraint("slip_BC") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        stageOneMassFlux[face]=gridMassFlux[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<StageOneMassFluxBoundarySlip>
    registerStageOneMassFluxBoundarySlip ;

  // Iteration advance rule for mass flux.

  // Rule to compute the stage 1 mass flux for symmetry faces.  Set to
  // the grid mass flux so there is identical cancellation.
  class StageOneMassFluxBoundarySymmetry : public pointwise_rule {
    private:
      const_store<real> gridMassFlux ;
      store<real> stageOneMassFlux ;
    public:

      // Define input and output.
      StageOneMassFluxBoundarySymmetry() {
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("symmetry::stageOneMassFlux",stageOneMassFlux) ;
        input("gridMassFlux") ;
        output("symmetry::stageOneMassFlux") ;
        constraint("symmetry_BC") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) { stageOneMassFlux[face]=gridMassFlux[face] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<StageOneMassFluxBoundarySymmetry>
    registerStageOneMassFluxBoundarySymmetry ;

  // Iteration advance rule for mass flux.
  class IterationAdvanceMassFlux : public pointwise_rule {
    private:
      const_store<real> massFluxCorrected ;
      store<real> massFluxIterationPlusOne ;
    public:

      // Define input and output.
      IterationAdvanceMassFlux() {
        name_store("massFluxCorrected{n,it}",massFluxCorrected) ;
        name_store("massFlux{n,it+1}",massFluxIterationPlusOne) ;
        input("massFluxCorrected{n,it}") ;
        output("massFlux{n,it+1}") ;
        constraint("massFluxCorrected{n,it}") ;
      }

      // Assign mass flux at end of iteration for a single face.
      void calculate(Entity face) {
        massFluxIterationPlusOne[face]=massFluxCorrected[face] ;
      }

      // Assign mass flux at end of iteration for a sequence of faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationAdvanceMassFlux> registerIterationAdvanceMassFlux ;

  // Iteration collapse rule for mass flux. NOTE: This rule worked for
  // incompressible flow, but for compressible, the sequence was the empty set.
  // Turns out we should not use iterationFinished variable as input, only
  // as a constraint. Why did it work for incompressible flow though. Strange.
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
        constraint("massFlux{n,it}") ;
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

  register_rule<TimeCollapseFace> registerTimeCollapseFace ;

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
