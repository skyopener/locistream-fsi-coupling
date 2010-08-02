//-----------------------------------------------------------------------------
// Description: This file contains rules for the momentum equation to allow
//   for reference frame changes for cells on either side of an iterior face.
//-----------------------------------------------------------------------------

// Standard library includes.
#include <vector>
using std::vector ;
                                                                                
// Loci includes.
#include <Loci.h>
using Loci::Area ;
                                                                                
// StreamUns includes.
#include "referenceFrame.h"
#include "rotorStator.h"
#include "sciTypes.h"
                                                                                
namespace streamUns {

  // Rule to compute the absolute velocity for all cells.
  class AbsoluteVelocity : public pointwise_rule {
    private:
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<unsigned int> cellReferenceFrame ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> v ;
      store<vect3d> vAbs ;
    public:
                                                                                
      // Define input and output.
      AbsoluteVelocity() {
        name_store("referenceFrame",referenceFrame) ;
        name_store("cellReferenceFrame",cellReferenceFrame) ;
        name_store("cellcenter",cellCenter) ;
        name_store("v",v) ;
        name_store("vAbs",vAbs) ;
        input("referenceFrame,cellReferenceFrame,cellcenter,v") ;
        output("vAbs") ;
        constraint("geom_cells") ;
      } ;
                                                                                
      // Calculate the absolute velocity for a cell.
      void calculate(Entity cell) {
        unsigned int n=cellReferenceFrame[cell] ;
        vect3d delta=(*referenceFrame)[n].axisEnd-(*referenceFrame)[n].
          axisStart,omega=((*referenceFrame)[n].omega/norm(delta))*delta,
          r=cellCenter[cell]-(*referenceFrame)[n].axisStart ;
        vAbs[cell]=v[cell]+cross(omega,r) ;
      }
                                                                                
      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;
                                                                                
  register_rule<AbsoluteVelocity> registerAbsoluteVelocity ;

  // Rule to compute the absolute velocity for all boundary faces.
  class AbsoluteVelocityBoundary : public pointwise_rule {
    private:
      const_Map ci ;
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<unsigned int> cellReferenceFrame ;
      const_store<vect3d> faceCenter ;
      const_store<vect3d> v_f ;
      store<vect3d> vAbs_f ;
    public:
                                                                                
      // Define input and output.
      AbsoluteVelocityBoundary() {
        name_store("ci",ci) ;
        name_store("referenceFrame",referenceFrame) ;
        name_store("cellReferenceFrame",cellReferenceFrame) ;
        name_store("facecenter",faceCenter) ;
        name_store("v_f",v_f) ;
        name_store("vAbs_f",vAbs_f) ;
        input("referenceFrame,ci->cellReferenceFrame,facecenter,v_f") ;
        output("vAbs_f") ;
        constraint("boundaryFaces") ;
      } ;
                                                                                
      // Calculate the absolute velocity for a face.
      void calculate(Entity face) {
        unsigned int n=cellReferenceFrame[ci[face]] ;
        vect3d delta=(*referenceFrame)[n].axisEnd-(*referenceFrame)[n].
          axisStart,omega=((*referenceFrame)[n].omega/norm(delta))*delta,
          r=faceCenter[face]-(*referenceFrame)[n].axisStart ;
        vAbs_f[face]=v_f[face]+cross(omega,r) ;
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;
                                                                                
  register_rule<AbsoluteVelocityBoundary> registerAbsoluteVelocityBoundary ;

  // Rule to modify the convection flux for faces across which a change of
  // reference frame occurs.
  class ConvectiveFluxCorrectionToVelocitySourceTerm : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
      const_param<real> thetaParameter ;
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<unsigned int> cellReferenceFrame ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      store<vect3d> vSourceTerm ;
    public:
                                                                                
      // Define input and output.
      ConvectiveFluxCorrectionToVelocitySourceTerm() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("referenceFrame",referenceFrame) ;
        name_store("cellReferenceFrame",cellReferenceFrame) ;
        name_store("facecenter",faceCenter) ;
        name_store("massFlux",massFlux) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("thetaParameter,referenceFrame,(cl,cr)->cellReferenceFrame") ;
        input("facecenter,massFlux") ;
        output("(cl,cr)->vSourceTerm") ;
        constraint("internalFaces,turbomachinery") ;
      }
                                                                                
      // Need to modify the convection term only on faces where a change of
      // reference frame occurs. For positive mass flux, only need to modify
      // flux into cr and vice-versa for cl.
      void calculate(Entity face) {
        if(cellReferenceFrame[cl[face]]!=cellReferenceFrame[cr[face]]){
          unsigned int nL=cellReferenceFrame[cl[face]],nR=cellReferenceFrame
            [cr[face]] ;
          vect3d deltaL=(*referenceFrame)[nL].axisEnd-(*referenceFrame)[nL].
            axisStart,omegaL=((*referenceFrame)[nL].omega/norm(deltaL))*deltaL,
            rL=faceCenter[face]-(*referenceFrame)[nL].axisStart ;
          vect3d deltaR=(*referenceFrame)[nR].axisEnd-(*referenceFrame)[nR].
            axisStart,omegaR=((*referenceFrame)[nR].omega/norm(deltaR))*deltaR,
            rR=faceCenter[face]-(*referenceFrame)[nR].axisStart ;
          if(massFlux[face]>0.0){
            vSourceTerm[cr[face]]+=massFlux[face]*(cross(omegaL,rL)-
              cross(omegaR,rR))*(*thetaParameter) ;
          }else{
            vSourceTerm[cl[face]]-=massFlux[face]*(cross(omegaR,rR)-
              cross(omegaL,rL))*(*thetaParameter) ;
          }
        }
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<ConvectiveFluxCorrectionToVelocitySourceTerm>
    registerConvectiveFluxCorrectionToVelocitySourceTerm ;

  // Rule to modify the implicit part of the diffusive flux for faces across
  //  which a change of reference frame occurs.
  class ImplicitDiffusiveFluxCorrectionToVelocitySourceTerm : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
      const_param<real> thetaParameter ;
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<unsigned int> cellReferenceFrame ;
      const_store<vect3d> cellCenter ;
      const_store<real> faceRadius ;
      const_store<real> diffusionProduct ;
      const_store<real> laminarViscosity ;
      const_store<real> eddyViscosity ;
      store<vect3d> vSourceTerm ;
    public:
                                                                                
      // Define input and output.
      ImplicitDiffusiveFluxCorrectionToVelocitySourceTerm() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("referenceFrame",referenceFrame) ;
        name_store("cellReferenceFrame",cellReferenceFrame) ;
        name_store("cellcenter",cellCenter) ;
        name_store("faceRadius",faceRadius) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("muu(temperature,p,y)",laminarViscosity) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("thetaParameter,faceRadius,diffusionProduct") ;
        input("referenceFrame,(cl,cr)->(cellReferenceFrame,cellcenter)") ;
        input("(cl,cr)->(muu(temperature,p,y),eddyViscosity)") ;
        output("(cl,cr)->vSourceTerm") ;
        constraint("internalFaces,viscousFlow,turbomachinery") ;
      }
                                                                                
      // Need to modify the diffusion term only on faces where a change of
      // reference frame occurs.
      void calculate(Entity face) {
        if(cellReferenceFrame[cl[face]]!=cellReferenceFrame[cr[face]]){
          unsigned int nL=cellReferenceFrame[cl[face]],nR=cellReferenceFrame
            [cr[face]] ;
          vect3d deltaL=(*referenceFrame)[nL].axisEnd-(*referenceFrame)[nL].
            axisStart,omegaL=((*referenceFrame)[nL].omega/norm(deltaL))*deltaL,
            rLInL=cellCenter[cl[face]]-(*referenceFrame)[nL].axisStart,
            rLInR=cellCenter[cl[face]]-(*referenceFrame)[nR].axisStart ;
          vect3d deltaR=(*referenceFrame)[nR].axisEnd-(*referenceFrame)[nR].
            axisStart,omegaR=((*referenceFrame)[nR].omega/norm(deltaR))*deltaR,
            rRInR=cellCenter[cr[face]]-(*referenceFrame)[nR].axisStart,
            rRInL=cellCenter[cr[face]]-(*referenceFrame)[nL].axisStart ;
          real temp=0.5*(laminarViscosity[cl[face]]+eddyViscosity[cl[face]]+
            laminarViscosity[cr[face]]+eddyViscosity[cr[face]])*
            diffusionProduct[face]*(*thetaParameter)*faceRadius[face] ;
          vSourceTerm[cl[face]]+=temp*(cross(omegaR,rRInR)-cross(omegaL,
            rRInL)) ;
          vSourceTerm[cr[face]]+=temp*(cross(omegaL,rLInL)-cross(omegaR,
            rLInR)) ;
        }
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<ImplicitDiffusiveFluxCorrectionToVelocitySourceTerm>
    registerImplicitDiffusiveFluxCorrectionToVelocitySourceTerm ;

  // Rule to REMOVE the diffusive flux contribution to the velocity source term
  // for interior faces that separate two different reference frames. Use of
  // this rule amounts to neglecting the secondary diffusion term on these
  // faces which should result in an extremely small unnoticable error. This
  // rule simply removes the contribution that was added in the /src directory
  // rule, but only at rotor/stator interfaces.
  class RemoveDiffusiveFluxFromVelocitySourceTermInterior : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
      const_param<real> thetaParameter ;
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<unsigned int> cellReferenceFrame ;
      const_store<vect3d> cellCenter ;
      const_store<real> vol ;
      const_store<tens3d> vGradient ;
      const_store<real> laminarViscosity ;
      const_store<real> eddyViscosity ;
      const_store<real> diffusionProduct ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<vect3d> vSourceTerm ;
    public:
                                                                                
      // Define input and output.
      RemoveDiffusiveFluxFromVelocitySourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("referenceFrame",referenceFrame) ;
        name_store("cellReferenceFrame",cellReferenceFrame) ;
        name_store("cellcenter",cellCenter) ;
        name_store("vol",vol) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("muu(temperature,p,y)",laminarViscosity) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("thetaParameter,referenceFrame") ;
        input("(cl,cr)->(cellReferenceFrame,cellcenter,vol,gradv3d(v))") ;
        input("(cl,cr)->muu(temperature,p,y),(cl,cr)->eddyViscosity") ;
        input("area,diffusionProduct,faceRadius") ;
        output("(cl,cr)->vSourceTerm") ;
        constraint("internalFaces,viscousFlow,turbomachinery") ;
      }
      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        if(cellReferenceFrame[cl[face]]!=cellReferenceFrame[cr[face]]){
          tens3d tempTensor=product(vol[cr[face]],vGradient[cl[face]])+
            product(vol[cl[face]],vGradient[cr[face]]) ;
          vect3d secondarySourceTerm=0.5*(laminarViscosity[cl[face]]+
            eddyViscosity[cl[face]]+laminarViscosity[cr[face]]+
            eddyViscosity[cr[face]])*dotTemp(tempTensor,(area[face].n*
            area[face].sada-diffusionProduct[face]*(cellCenter[cr[face]]-
            cellCenter[cl[face]])))*(*thetaParameter)*faceRadius[face]/
            (vol[cl[face]]+vol[cr[face]]) ;
          vSourceTerm[cl[face]]-=secondarySourceTerm ;
          vSourceTerm[cr[face]]+=secondarySourceTerm ;
        }
      }
                                                                                
      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<RemoveDiffusiveFluxFromVelocitySourceTermInterior>
    registerRemoveDiffusiveFluxFromVelocitySourceTermInterior ;

  // Rule to remove the components of viscous flux that are non-zero only for
  // compressible flows at the rotor/stator interfaces.
  class RemoveDiffusiveFluxFromVelocitySourceTermCompressibleInterior : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
      const_param<real> thetaParameter ;
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<unsigned int> cellReferenceFrame ;
      const_store<real> vol ;
      const_store<tens3d> vGradient ;
      const_store<real> laminarViscosity ;
      const_store<real> eddyViscosity ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<vect3d> vSourceTerm ;
    public:
                                                                                
      // Define input and output.
      RemoveDiffusiveFluxFromVelocitySourceTermCompressibleInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("referenceFrame",referenceFrame) ;
        name_store("cellReferenceFrame",cellReferenceFrame) ;
        name_store("vol",vol) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("muu(temperature,p,y)",laminarViscosity) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("thetaParameter,referenceFrame,faceRadius") ;
        input("(cl,cr)->(cellReferenceFrame,vol,gradv3d(v))") ;
        input("(cl,cr)->muu(temperature,p,y),(cl,cr)->eddyViscosity,area") ;
        output("(cl,cr)->vSourceTerm") ;
        constraint("internalFaces,viscousFlow,compressibleFlow") ;
        constraint("turbomachinery") ;
      }
                                                                                
      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        if(cellReferenceFrame[cl[face]]!=cellReferenceFrame[cr[face]]){
          real vR=vol[cr[face]]/(vol[cl[face]]+vol[cr[face]]) ;
          real vL=vol[cl[face]]/(vol[cl[face]]+vol[cr[face]]) ;
                                                                                
// Change for Loci-3-1.
//        tens3d tempTensor=Transpose(vGradient[cl[face]]*vR+vGradient
//          [cr[face]]*vL) ;
          tens3d tempTensor=Transpose(product(vR,vGradient[cl[face]])+
            product(vL,vGradient[cr[face]])) ;
                                                                                
          real temp=2.0*Trace(tempTensor)/3.0 ;
          tempTensor.x.x-=temp ; tempTensor.y.y-=temp ; tempTensor.z.z-=temp ;
          vect3d secondarySourceTerm=0.5*(laminarViscosity[cl[face]]+
            eddyViscosity[cl[face]]+laminarViscosity[cr[face]]+
            eddyViscosity[cr[face]])*dotTemp(tempTensor,area[face].n*
            area[face].sada)*(*thetaParameter)*faceRadius[face] ;
          vSourceTerm[cl[face]]-=secondarySourceTerm ;
          vSourceTerm[cr[face]]+=secondarySourceTerm ;
        }
      }
                                                                                
      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<RemoveDiffusiveFluxFromVelocitySourceTermCompressibleInterior>
    registerRemoveDiffusiveFluxFromVelocitySourceTermCompressibleInterior ;

  // Rule to add centipedal acceleration term to source term.
  class CoriolisCentripetalToVelocitySourceTerm : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<unsigned int> cellReferenceFrame ;
      const_store<vect3d> cellCenter ;
      const_store<real> rho ;
      const_store<real> vol ;
      const_store<vect3d> v ;
      store<vect3d> vSourceTerm ;
    public:
                                                                                
      // Define input and output.
      CoriolisCentripetalToVelocitySourceTerm() {
        name_store("referenceFrame",referenceFrame) ;
        name_store("cellReferenceFrame",cellReferenceFrame) ;
        name_store("cellcenter",cellCenter) ;
        name_store("rho",rho) ;
        name_store("vol",vol) ;
        name_store("v",v) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("referenceFrame,cellReferenceFrame") ;
        input("cellcenter,rho,vol,v") ;
        output("vSourceTerm") ;
        constraint("geom_cells,turbomachinery") ;
      } ;
                                                                                
      // Increment the source term for a cell.
      void calculate(Entity cell) {
        unsigned int n=cellReferenceFrame[cell] ;
        vect3d delta=(*referenceFrame)[n].axisEnd-(*referenceFrame)[n].
          axisStart,omega=((*referenceFrame)[n].omega/norm(delta))*delta,
          r=cellCenter[cell]-(*referenceFrame)[n].axisStart ;
        vSourceTerm[cell]-=(2.0*cross(omega,v[cell])+cross(omega,
          cross(omega,r)))*rho[cell]*vol[cell] ;
      }
                                                                                
      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;
                                                                                
  register_rule<CoriolisCentripetalToVelocitySourceTerm>
    registerCoriolisCentripetalToVelocitySourceTerm ;

  // Rule to determine the reference frame for each cell.
  class ImposedPressureGradientToVelocitySourceTerm : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_param<RotorData> rotorData ;
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<real> vol ;
      const_store<vect3d> cellCenter ;
      const_store<unsigned int> cellReferenceFrame ;
      store<vect3d> vSourceTerm ;
    public:
                                                                                
      // Define input and output.
      ImposedPressureGradientToVelocitySourceTerm() {
        name_store("rotorData",rotorData) ;
        name_store("referenceFrame",referenceFrame) ;
        name_store("vol",vol) ;
        name_store("cellcenter",cellCenter) ;
        name_store("cellReferenceFrame",cellReferenceFrame) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("rotorData,referenceFrame,vol,cellcenter") ;
        input("cellReferenceFrame") ;
        output("vSourceTerm") ;
        constraint("geom_cells") ;
      }
                                                                                
      // Impose only specified non-zero pressure gradient.
      void calculate(Entity cell) {
        const real dpdzA=rotorData->DPDZA(cellCenter[cell],*referenceFrame),
          dpdzB=rotorData->DPDZB(cellCenter[cell],*referenceFrame) ;
        if(dpdzA!=0.0 || dpdzB!=0.0){

          // Compute the axis unit vector.
          int referenceFrameNumber=cellReferenceFrame[cell] ;
          vect3d axis=(*referenceFrame)[referenceFrameNumber].axisEnd-
            (*referenceFrame)[referenceFrameNumber].axisStart ;
          axis/=norm(axis) ;

          // Compute dpdz at the cell radius.
          vect3d r=cellCenter[cell]-(*referenceFrame)[referenceFrameNumber].
            axisStart ;
          real radius=norm(r-dot(r,axis)*axis) ;
          real dpdz=dpdzA+radius*dpdzB ;

          // Remove old pressure gradient and add new.
          vSourceTerm[cell]-=dpdz*axis*vol[cell] ;
        }
      }
                                                                                
      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<ImposedPressureGradientToVelocitySourceTerm>
    registerImposedPressureGradientToVelocitySourceTerm ;

}
