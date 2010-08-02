//-----------------------------------------------------------------------------
// Description: This file contains rules for the fully-implicit
//   Crank-Nicholson scheme for all the governing equations.
//-----------------------------------------------------------------------------

// Standard library includes.
#include <vector>
using std::vector ;

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// CHEM includes.
#include "eos.h"

// StreamUns includes.
#include "const.h"
#include "referenceFrame.h"
#include "residual.h"
#include "sciTypes.h"
#include "varsFileInputs.h"

namespace streamUns {

//-----------------------------------------------------------------------------
// General rules.

  // Creates the time integrator constraints.
  class ThetaParameter : public singleton_rule {
    private:
      const_param<string> timeIntegrator ;
      param<real> thetaParameter ;
    public:

      // Define input and output.
      ThetaParameter() {
        name_store("timeIntegrator",timeIntegrator) ;
        name_store("thetaParameter",thetaParameter) ;
        input("timeIntegrator") ;
        output("thetaParameter") ;
      }

      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        *thetaParameter=1.0 ; if(*timeIntegrator=="CN") *thetaParameter=0.5 ;
      }
  } ;

  register_rule<ThetaParameter> registerThetaParameter ;

//-----------------------------------------------------------------------------
// Rules for the momentum equation.

  // Rule to initialize the inviscid contribution from the previous time
  // level. This is defined as the inviscid flux out of cell 'cl' .
  class InitializeOldInviscidFlux : public unit_rule {
    private:
      store<vect3d> oldInviscidFlux ;
    public:
                                                                                
      // Define input and output.
      InitializeOldInviscidFlux() {
        name_store("oldInviscidFlux",oldInviscidFlux) ;
        output("oldInviscidFlux") ;
        constraint("UNIVERSE") ;
      }
                                                                                
      // Initialize.
      void calculate(Entity face) {
        oldInviscidFlux[face]=vector3d<real>(0.0,0.0,0.0) ;
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<InitializeOldInviscidFlux> registerInitializeOldInviscidFlux ;

  // Rule to add the first-order contribution from the previous
  // time level to the inviscid flux for interior faces.
  class OldFOUInviscidFluxInterior : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
      const_store<vect3d> v ;
      const_store<real> massFlux ;
      store<vect3d> oldInviscidFlux ;
    public:

      // Define input and output.
      OldFOUInviscidFluxInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("v",v) ;
        name_store("massFlux",massFlux) ;
        name_store("oldInviscidFlux",oldInviscidFlux) ;
        input("(cl,cr)->v,massFlux") ;
        output("oldInviscidFlux") ;
        constraint("internalFaces") ;
      }

      // Increment the flux for a single face.
      void calculate(Entity face) {
        oldInviscidFlux[face]+=(massFlux[face]>0.0)? massFlux[face]*v[cl[face]]:
          massFlux[face]*v[cr[face]] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<OldFOUInviscidFluxInterior> registerOldFOUInviscidFluxInterior ;

  // Rule to add the second-order contribution from the previous
  // time level to the inviscid flux for interior faces.
  class OldSOUInviscidFluxInterior : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
      const_store<vect3d> v ;
      const_store<tens3d> vGradient ;
      const_store<vect3d> vLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      store<vect3d> oldInviscidFlux ;
    public:

      // Define input and output.
      OldSOUInviscidFluxInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("v",v) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("limiterv3d(v)",vLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("massFlux",massFlux) ;
        name_store("oldInviscidFlux",oldInviscidFlux) ;
        input("(cl,cr)->(v,gradv3d(v),limiterv3d(v),cellcenter)") ;
        input("facecenter,massFlux") ;
        output("oldInviscidFlux") ;
        constraint("souOrRoeInviscidFlux,internalFaces") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        oldInviscidFlux[face]+=(massFlux[face]>0.0)? massFlux[face]*
          ComponentProduct(vLimiter[cl[face]],dotTemp(vGradient[cl[face]],
          faceCenter[face]-cellCenter[cl[face]])):massFlux[face]*
          ComponentProduct(vLimiter[cr[face]],dotTemp(vGradient[cr[face]],
          faceCenter[face]-cellCenter[cr[face]])) ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<OldSOUInviscidFluxInterior> registerOldSOUInviscidFluxInterior ;

  // Rule to add the first-order inviscid contribution from the previous time
  // level to the source term for boundary faces.
  class OldFOUInviscidFluxBoundary : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_Map ci ;
      const_store<vect3d> v ;
      const_store<vect3d> v_f ;
      const_store<real> massFlux ;
      store<vect3d> oldInviscidFlux ;
    public:

      // Define input and output.
      OldFOUInviscidFluxBoundary() {
        name_store("ci",ci) ;
        name_store("v",v) ;
        name_store("v_f",v_f) ;
        name_store("massFlux",massFlux) ;
        name_store("oldInviscidFlux",oldInviscidFlux) ;
        input("ci->v,v_f,massFlux") ;
        output("oldInviscidFlux") ;
        constraint("boundaryFaces") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        oldInviscidFlux[face]+=(massFlux[face]>0.0)? massFlux[face]*v[ci[face]]:
          massFlux[face]*v_f[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<OldFOUInviscidFluxBoundary> registerOldFOUInviscidFluxBoundary ;

  // Rule to add the second-order contribution from the previous
  // time level to the inviscid flux for boundary faces.
  class OldSOUInviscidFluxBoundary : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_Map ci ;
      const_store<vect3d> v ;
      const_store<tens3d> vGradient ;
      const_store<vect3d> vLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      store<vect3d> oldInviscidFlux ;
    public:

      // Define input and output.
      OldSOUInviscidFluxBoundary() {
        name_store("ci",ci) ;
        name_store("v",v) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("limiterv3d(v)",vLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("massFlux",massFlux) ;
        name_store("oldInviscidFlux",oldInviscidFlux) ;
        input("ci->(v,gradv3d(v),limiterv3d(v),cellcenter)") ;
        input("facecenter,massFlux") ;
        output("oldInviscidFlux") ;
        constraint("souOrRoeInviscidFlux,boundaryFaces") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        if(massFlux[face]>0.0) oldInviscidFlux[face]+=massFlux[face]*
          ComponentProduct(vLimiter[ci[face]],dotTemp(vGradient[ci[face]],
          faceCenter[face]-cellCenter[ci[face]])) ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<OldSOUInviscidFluxBoundary> registerOldSOUInviscidFluxBoundary ;

  // Rule to compute the pressure for interior faces.
  class FacePressureInterior : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<real> p ;
      const_store<vect3d> pGradient ;
      const_store<real> pLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      store<real> facePressure ;
    public:
                                                                                
      // Define input and output.
      FacePressureInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("p",p) ;
        name_store("grads(p)",pGradient) ;
        name_store("limiters(p)",pLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("facePressure",facePressure) ;
        input("(cl,cr)->(p,grads(p),limiters(p),cellcenter),facecenter") ;
        output("facePressure") ;
        constraint("internalFaces") ;
      }
                                                                                
      // Compute the pressure at the face using extrapolation with limiters.
      void calculate(Entity face) {
        facePressure[face]=0.5*((p[cl[face]]+pLimiter[cl[face]]*
          dot(pGradient[cl[face]],(faceCenter[face]-cellCenter[cl[face]])))+
          (p[cr[face]]+pLimiter[cr[face]]*dot(pGradient[cr[face]],
          (faceCenter[face]-cellCenter[cr[face]])))) ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<FacePressureInterior> registerFacePressureInterior ;

  // Rule to compute the pressure for boundary faces.
  class FacePressureBoundary : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> p ;
      const_store<vect3d> pGradient ;
      const_store<real> pLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      store<real> facePressure ;
    public:
                                                                                
      // Define input and output.
      FacePressureBoundary() {
        name_store("ci",ci) ;
        name_store("p",p) ;
        name_store("grads(p)",pGradient) ;
        name_store("limiters(p)",pLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("facePressure",facePressure) ;
        input("ci->(p,grads(p),limiters(p),cellcenter),facecenter") ;
        output("facePressure") ;
        constraint("boundaryFaces") ;
      }
                                                                                
      // Compute the pressure at the face using extrapolation with limiters.
      void calculate(Entity face) {
        facePressure[face]=(p[ci[face]]+pLimiter[ci[face]]*
          dot(pGradient[ci[face]],(faceCenter[face]-cellCenter[ci[face]]))) ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<FacePressureBoundary> registerFacePressureBoundary ;

  // Rule to compute the total diffusive flux through each interior face. The
  // flux is defined as the flux into cell 'cl'.
  class DiffusiveFluxInterior : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<vect3d> cellCenter ;
      const_store<Area> area ;
      const_store<real> diffusionProduct ;
      const_store<real> laminarViscosity ;
      const_store<real> eddyViscosity ;
      const_store<vect3d> v ;
      const_store<tens3d> vGradient ;
      store<vect3d> diffusiveFlux ;
    public:
                                                                                
      // Define input and output.
      DiffusiveFluxInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("cellcenter",cellCenter) ;
        name_store("area",area) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("v",v) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("diffusiveFlux",diffusiveFlux) ;
        input("area,diffusionProduct") ;
        input("(cl,cr)->(cellcenter,laminarViscosity,eddyViscosity)") ;
        input("(cl,cr)->(v,gradv3d(v))") ;
        output("diffusiveFlux") ;
        constraint("internalFaces,viscousFlow") ;
      }
                                                                                
      // Compute the diffusive flux for a single face.
      void calculate(Entity face) {
        real temp=0.5*(laminarViscosity[cl[face]]+eddyViscosity[cl[face]]+
          laminarViscosity[cr[face]]+eddyViscosity[cr[face]])*
          diffusionProduct[face] ;
// Change required for Loci-3-1. Seems implementing tensor3d in terms of
// vector3d causes some unexpected compiler issues.
//      vect3d secondarySourceTerm=0.5*(laminarViscosity[cl[face]]+eddyViscosity
//        [cl[face]]+laminarViscosity[cr[face]]+eddyViscosity[cr[face]])*
//        dotTemp(0.5*(vGradient[cl[face]]+vGradient[cr[face]]),(area[face].n*
//        area[face].sada-diffusionProduct[face]*(cellCenter[cr[face]]-
//        cellCenter[cl[face]]))) ;
        tens3d tempTensor=(vGradient[cl[face]]+vGradient[cr[face]]) ;
        vect3d secondarySourceTerm=0.25*(laminarViscosity[cl[face]]+
          eddyViscosity[cl[face]]+laminarViscosity[cr[face]]+
          eddyViscosity[cr[face]])*dotTemp(tempTensor,(area[face].n*
          area[face].sada-diffusionProduct[face]*(cellCenter[cr[face]]-
          cellCenter[cl[face]]))) ;
        diffusiveFlux[face]=temp*(v[cr[face]]-v[cl[face]])+secondarySourceTerm ;
      }
                                                                                
      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<DiffusiveFluxInterior> registerDiffusiveFluxInterior ;

  // Rule to compute the total diffusive flux through each boundary face that
  // is not using wall functions. The flux is defined as the flux into 'ci'.
  class DiffusiveFluxBoundary : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> v ;
      const_store<tens3d> vGradient ;
      const_store<vect3d> faceCenter ;
      const_store<real> diffusionProduct ;
      const_store<vect3d> v_f ;
      const_store<real> laminarViscosity_f ;
      const_store<real> eddyViscosity_f ;
      const_store<Area> area ;
      store<vect3d> diffusiveFlux ;
    public:
                                                                                
      // Define input and output.
      DiffusiveFluxBoundary() {
        name_store("ci",ci) ;
        name_store("cellcenter",cellCenter) ;
        name_store("v",v) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("facecenter",faceCenter) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("v_f",v_f) ;
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        name_store("eddyViscosity_f",eddyViscosity_f) ;
        name_store("area",area) ;
        name_store("diffusiveFlux",diffusiveFlux) ;
        input("ci->(cellcenter,v,gradv3d(v))") ;
        input("facecenter,diffusionProduct,v_f") ;
        input("laminarViscosity_f,eddyViscosity_f,area") ;
        output("diffusiveFlux") ;
        constraint("boundaryVelocityDiffusion,viscousFlow,noWallFunction_BC") ;
      }

      // Compute the diffusive flux for a single face.
      void calculate(Entity face) {
        diffusiveFlux[face]=(laminarViscosity_f[face]+eddyViscosity_f[face])*
          ((v_f[face]-v[ci[face]])*diffusionProduct[face]+
          dotTemp(vGradient[ci[face]],(area[face].n*area[face].sada-
          diffusionProduct[face]*(faceCenter[face]-cellCenter[ci[face]])))) ;
      }
                                                                                
      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<DiffusiveFluxBoundary> registerDiffusiveFluxBoundary ;

  // Rule to compute the total diffusive flux through each boundary face that
  // is using wall functions. The flux is defined as the flux into 'ci'.
  class DiffusiveFluxBoundaryWallFunction : public pointwise_rule {
    private:
      const_store<vect3d> tauWall ;
      const_store<Area> area ;
      store<vect3d> diffusiveFlux ;
    public:
                                                                                
      // Define input and output.
      DiffusiveFluxBoundaryWallFunction() {
        name_store("tauWall",tauWall) ;
        name_store("area",area) ;
        name_store("diffusiveFlux",diffusiveFlux) ;
        input("tauWall,area") ;
        output("diffusiveFlux") ;
        constraint("viscousFlow,ref->wallFunction_BCoption") ;
      }
                                                                                
      // Compute the diffusive flux for a single face.
      void calculate(Entity face) {
        diffusiveFlux[face]=tauWall[face]*area[face].sada ;
      }
                                                                                
      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<DiffusiveFluxBoundaryWallFunction>
    registerDiffusiveFluxBoundaryWallFunction ;

  // Rule to initialize the source term from the previous time level.
  class InitializeOldVelocitySourceTerm : public unit_rule {
    private:
      store<vect3d> oldVSourceTerm ;
    public:
                                                                                
      // Define input and output.
      InitializeOldVelocitySourceTerm() {
        name_store("oldVSourceTerm",oldVSourceTerm) ;
        output("oldVSourceTerm") ;
        constraint("UNIVERSE") ;
      }
                                                                                
      // Initialize.
      void calculate(Entity cell) {
        oldVSourceTerm[cell]=vector3d<real>(0.0,0.0,0.0) ;
      }
                                                                                
      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<InitializeOldVelocitySourceTerm>
    registerInitializeOldVelocitySourceTerm ;

  // Rule to add the inviscid contribution to the source term from the previous
  // time level for interior faces.
  class InviscidFluxToOldVelocitySourceTermInterior : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
      const_param<real> thetaParameter ;
      const_store<real> facePressure ;
      const_store<vect3d> oldInviscidFlux ;
      const_store<Area> area ;
      store<vect3d> oldVSourceTerm ;
    public:
                                                                                
      // Define input and output.
      InviscidFluxToOldVelocitySourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("facePressure",facePressure) ;
        name_store("oldInviscidFlux",oldInviscidFlux) ;
        name_store("area",area) ;
        name_store("oldVSourceTerm",oldVSourceTerm) ;
        input("thetaParameter,facePressure,oldInviscidFlux,area") ;
        output("(cl,cr)->oldVSourceTerm") ;
        constraint("internalFaces") ;
      }
                                                                                
      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        vect3d temp=(facePressure[face]*area[face].n*area[face].sada+
          oldInviscidFlux[face])*(*thetaParameter) ;
        oldVSourceTerm[cl[face]]-=temp ; oldVSourceTerm[cr[face]]+=temp ;
      }

                                                                                
      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<InviscidFluxToOldVelocitySourceTermInterior>
    registerInviscidFluxToOldVelocitySourceTermInterior ;

  // Rule to add the inviscid contribution to the source term from the previous
  // time level for boundary faces.
  class InviscidFluxToOldVelocitySourceTermBoundary : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map ci ;
      const_param<real> thetaParameter ;
      const_store<real> facePressure ;
      const_store<vect3d> oldInviscidFlux ;
      const_store<Area> area ;
      store<vect3d> oldVSourceTerm ;
    public:
                                                                                
      // Define input and output.
      InviscidFluxToOldVelocitySourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("facePressure",facePressure) ;
        name_store("oldInviscidFlux",oldInviscidFlux) ;
        name_store("area",area) ;
        name_store("oldVSourceTerm",oldVSourceTerm) ;
        input("thetaParameter,facePressure,oldInviscidFlux,area") ;
        output("ci->oldVSourceTerm") ;
        constraint("boundaryFaces") ;
      }
                                                                                
      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        vect3d temp=(facePressure[face]*area[face].n*area[face].sada+
          oldInviscidFlux[face])*(*thetaParameter) ;
        oldVSourceTerm[ci[face]]-=temp ;
      }

                                                                                
      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<InviscidFluxToOldVelocitySourceTermBoundary>
    registerInviscidFluxToOldVelocitySourceTermBoundary ;

  // Rule to add the diffusive contribution to the source term from the previous
  // time level for interior faces.
  class DiffusiveFluxToOldVelocitySourceTermInterior : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
      const_param<real> thetaParameter ;
      const_store<vect3d> diffusiveFlux ;
      store<vect3d> oldVSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToOldVelocitySourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("diffusiveFlux",diffusiveFlux) ;
        name_store("oldVSourceTerm",oldVSourceTerm) ;
        input("thetaParameter,diffusiveFlux") ;
        output("(cl,cr)->oldVSourceTerm") ;
        constraint("internalFaces") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        oldVSourceTerm[cl[face]]+=(*thetaParameter)*diffusiveFlux[face] ;
        oldVSourceTerm[cr[face]]-=(*thetaParameter)*diffusiveFlux[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToOldVelocitySourceTermInterior>
    registerDiffusiveFluxToOldVelocitySourceTermInterior ;

  // Rule to add the diffusive contribution to the source term from the previous
  // time level for boundary faces.
  class DiffusiveFluxToOldVelocitySourceTermBoundary : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map ci ;
      const_param<real> thetaParameter ;
      const_store<vect3d> diffusiveFlux ;
      store<vect3d> oldVSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToOldVelocitySourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("diffusiveFlux",diffusiveFlux) ;
        name_store("oldVSourceTerm",oldVSourceTerm) ;
        input("thetaParameter,diffusiveFlux") ;
        output("ci->oldVSourceTerm") ;
        constraint("boundaryVelocityDiffusion") ;
      }

      // Increment the source term for the cell attach to a single face.
      void calculate(Entity face) {
        oldVSourceTerm[ci[face]]+=(*thetaParameter)*diffusiveFlux[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToOldVelocitySourceTermBoundary>
    registerDiffusiveFluxToOldVelocitySourceTermBoundary ;

  // Rule to add the source term from the previous time level to the current
  // source term.
  class OldSourceTermToVelocitySourceTerm : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_store<vect3d> oldVSourceTerm ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      OldSourceTermToVelocitySourceTerm() {
        name_store("oldVSourceTerm{n}",oldVSourceTerm) ;
        name_store("vSourceTerm{n,it}",vSourceTerm) ;
        input("oldVSourceTerm{n}") ;
        output("vSourceTerm{n,it}") ;
        constraint("geom_cells{n,it},crankNicholsonIntegrator{n,it}") ;
      }

      // Increment the source term for a cell.
      void calculate(Entity cell) { vSourceTerm[cell]+=oldVSourceTerm[cell] ; }

      // Call calculate for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<OldSourceTermToVelocitySourceTerm>
    registerOldSourceTermToVelocitySourceTerm ;

  // Rule to add an additional term from subtracting the continuity equation.
  class ContinuityToVelocitySourceTerm : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_param<real> thetaParameter ;
      const_store<real> netMassFlux ;
      const_store<vect3d> v ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      ContinuityToVelocitySourceTerm() {
        name_store("thetaParameter{n,it}",thetaParameter) ;
        name_store("netMassFlux{n}",netMassFlux) ;
        name_store("v{n,it}",v) ;
        name_store("vSourceTerm{n,it}",vSourceTerm) ;
        input("thetaParameter{n,it},netMassFlux{n},v{n,it}") ;
        output("vSourceTerm{n,it}") ;
        constraint("geom_cells{n,it},crankNicholsonIntegrator{n,it}") ;
      }

      // Increment the source term for a cell.
      void calculate(Entity cell) {
        vSourceTerm[cell]+=(*thetaParameter)*netMassFlux[cell]*v[cell] ;
      }

      // Call calculate for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ContinuityToVelocitySourceTerm>
    registerContinuityToVelocitySourceTerm ;

//-----------------------------------------------------------------------------
// Rules for the pressure-correction equation.

  // Rule to add the old mass flux to the rhs.
  class PressureCorrectionOldMassFluxRHSInternal : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_param<real> thetaParameter ;
      const_store<real> massFlux ;
      store<real> B ;
    public:
                                                                                
      // Define input and output.
      PressureCorrectionOldMassFluxRHSInternal() {
        name_store("cl{n,it}",cl) ;
        name_store("cr{n,it}",cr) ;
        name_store("thetaParameter{n,it}",thetaParameter) ;
        name_store("massFlux{n}",massFlux) ;
        name_store("pPrime_B{n,it}",B) ;
        input("thetaParameter{n,it},massFlux{n}") ;
        output("(cl{n,it},cr{n,it})->pPrime_B{n,it}") ;
        constraint("internalFaces{n,it},crankNicholsonIntegrator{n,it}") ;
      }
                                                                                
      // Distribute contribution for a single face.
      void calculate(Entity face) {
        B[cl[face]]-=massFlux[face]*(*thetaParameter) ;
        B[cr[face]]+=massFlux[face]*(*thetaParameter) ;
      }
                                                                                
      // Distribute contribution for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<PressureCorrectionOldMassFluxRHSInternal>
    registerPressureCorrectionOldMassFluxRHSInternal ;

  // Rule to assemble the rhs term for the linear system for both incompressible  // and compressible flow. Assembling over boundary faces. Checked.
  class PressureCorrectionOldMassFluxRHSBoundary : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_param<real> thetaParameter ;
      const_store<real> massFlux ;
      store<real> B ;
    public:
                                                                                
      // Define input and output.
      PressureCorrectionOldMassFluxRHSBoundary() {
        name_store("ci{n,it}",ci) ;
        name_store("thetaParameter{n,it}",thetaParameter) ;
        name_store("massFlux{n}",massFlux) ;
        name_store("pPrime_B{n,it}",B) ;
        input("thetaParameter{n,it},massFlux{n}") ;
        output("ci{n,it}->pPrime_B{n,it}") ;
        constraint("boundaryFaces{n,it},crankNicholsonIntegrator{n,it}") ;
      }
                                                                                
      // Distribute contribution for a single face.
      void calculate(Entity face) {
        B[ci[face]]-=massFlux[face]*(*thetaParameter) ;
      }
                                                                                
      // Distribute contribution for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<PressureCorrectionOldMassFluxRHSBoundary>
    registerPressureCorrectionOldMassFluxRHSBoundary ;
}
