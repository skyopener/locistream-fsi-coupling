//-----------------------------------------------------------------------------
// Description: This file contains rules for the energy equation to allow
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
#include "sciTypes.h"
                                                                                
namespace streamUns {

// NOTE: For form "energy3" of the energy equation, no corrections to the
// source term are necessary since in this form the diffusion term is
// artificially obtained by adding zero to both sides, the l.h.s being
// implicit and the r.h.s. being explicit. At convergence these terms will
// cancel out. Until we see that there are any intermediate problems caused
// by this approach, we will go with this.

  // Rename rule to allow us to override the total enthalpy gradient for
  // turbomachinery flows. This rule assigns grads(h) while a priority rule
  // in the turbomachinery module assigns gradsTurbo(h).
  class HGradientTurbo : public pointwise_rule {
    private:
      store<vect3d> hGradient ;
    public:
                                                                                
      // Define input and output.
      HGradientTurbo() {
        name_store("gradsTurbo(h)",hGradient) ;
        input("gradsTurbo(h)") ;
        output("turbo::hGradient=gradsTurbo(h)") ;
        constraint("geom_cells,turbomachinery") ;
      }
                                                                                
      // Empty compute method.
      virtual void compute(const sequence &seq) {}
                                                                                
  } ;
                                                                                
  register_rule<HGradientTurbo> registerHGradientTurbo ;

  // Rule to modify the convection flux for faces across which a change of
  // reference frame occurs.
  class ConvectiveFluxCorrectionToTotalEnthalpySourceTerm : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_param<real> thetaParameter ;
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<unsigned int> cellReferenceFrame ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      store<real> hSourceTerm ;
    public:

      // Define input and output.
      ConvectiveFluxCorrectionToTotalEnthalpySourceTerm() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("referenceFrame",referenceFrame) ;
        name_store("cellReferenceFrame",cellReferenceFrame) ;
        name_store("facecenter",faceCenter) ;
        name_store("massFlux",massFlux) ;
        name_store("hSourceTerm",hSourceTerm) ;
        input("thetaParameter,referenceFrame,(cl,cr)->cellReferenceFrame") ;
        input("facecenter,massFlux") ;
        output("(cl,cr)->hSourceTerm") ;
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
          vect3d tempL=cross(omegaL,rL),tempR=cross(omegaR,rR) ;
          if(massFlux[face]>0.0){
            hSourceTerm[cr[face]]+=massFlux[face]*(dot(tempR,tempR)-
              dot(tempL,tempL))*(*thetaParameter) ;
          }else{
            hSourceTerm[cl[face]]-=massFlux[face]*(dot(tempL,tempL)-
              dot(tempR,tempR))*(*thetaParameter) ;
          }
        }
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ConvectiveFluxCorrectionToTotalEnthalpySourceTerm>
    registerConvectiveFluxCorrectionToTotalEnthalpySourceTerm ;

  // Rule to remove the viscous stress contribution from the source term for
  // interior faces separating two different frames of reference.
  class RemoveViscousStressFromTotalEnthalpySourceTermInterior : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<unsigned int> cellReferenceFrame ;
      const_store<real> vol ;
      const_store<vect3d> v ;
      const_store<tens3d> vGradient ;
      const_store<real> laminarViscosity,eddyViscosity ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<real> hSourceTerm ;
    public:

      // Define input and output.
      RemoveViscousStressFromTotalEnthalpySourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("referenceFrame",referenceFrame) ;
        name_store("cellReferenceFrame",cellReferenceFrame) ;
        name_store("vol",vol) ;
        name_store("v",v) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("muu(temperature,p,y)",laminarViscosity) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("hSourceTerm",hSourceTerm) ;
        input("referenceFrame,(cl,cr)->(cellReferenceFrame,vol,v,gradv3d(v))") ;
        input("(cl,cr)->muu(temperature,p,y),(cl,cr)->eddyViscosity") ;
        input("area,faceRadius") ;
        output("(cl,cr)->hSourceTerm") ;
        constraint("internalFaces,viscousFlow,energy34,turbomachinery") ;
      }

      // Remove what was previosly added for the cells attach to a single face.
      void calculate(Entity face) {
        if(cellReferenceFrame[cl[face]]!=cellReferenceFrame[cr[face]]){
          real vR=vol[cr[face]]/(vol[cl[face]]+vol[cr[face]]) ;
          real vL=vol[cl[face]]/(vol[cl[face]]+vol[cr[face]]) ;

// Change for Loci-3-1.
//        tens3d tempVGradient=vGradient[cl[face]]*vR+vGradient[cr[face]]*vL ;
          tens3d tempVGradient=product(vR,vGradient[cl[face]])+product(vL,
            vGradient[cr[face]]) ;

          tens3d tempStress=tempVGradient+Transpose(tempVGradient) ;
          real temp=2.0*Trace(tempVGradient)/3.0 ;
          tempStress.x.x-=temp ; tempStress.y.y-=temp ; tempStress.z.z-=temp ;
          real sourceTerm=0.25*area[face].sada*(laminarViscosity[cl[face]]+
            eddyViscosity[cl[face]]+laminarViscosity[cr[face]]+eddyViscosity
            [cr[face]])*dot(dotTemp(tempStress,v[cl[face]]+v[cr[face]]),
            area[face].n)*faceRadius[face] ;
          hSourceTerm[cl[face]]-=sourceTerm ;
          hSourceTerm[cr[face]]+=sourceTerm ;
        }
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<RemoveViscousStressFromTotalEnthalpySourceTermInterior>
    registerRemoveViscousStressFromTotalEnthalpySourceTermInterior ;

}
