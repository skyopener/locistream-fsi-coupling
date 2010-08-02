//-----------------------------------------------------------------------------
// Description: This file contains rules for the momentum, pressure, energy
// and species equations when using the Roe scheme for the inviscid flux.
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

//-----------------------------------------------------------------------------
// Rules for computing Roe-averaged face values.

  // Rule for computing the face density used for Roe-averaging.
  class RoeFaceDensity : public pointwise_rule {
    private:
      const_store<real> rhoLeft ;
      const_store<real> rhoRight ;
      store<real> rhoTilde ;
    public:
                                                                                
      // Define input and output.
      RoeFaceDensity() {
        name_store("lefts(rho)",rhoLeft) ;
        name_store("rights(rho)",rhoRight) ;
        name_store("rhoTilde",rhoTilde) ;
        input("lefts(rho),rights(rho)") ;
        output("rhoTilde") ;
        constraint("internalFaces") ;
      }
                                                                                
      // Limited extrapolation from the cell to the face.
      void calculate(Entity face) {
        rhoTilde[face]=sqrt(rhoLeft[face]*rhoRight[face]) ;
      }
                                                                                
      // Call calculate for a sequence of faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<RoeFaceDensity> registerRoeFaceDensity ;

  // Rule for computing Roe-averaged scalar values at the faces.
  class RoeAveragedScalar : public pointwise_rule {
    private:
      const_store<real> rhoTilde ;
      const_store<real> xLeft ;
      const_store<real> xRight ;
      store<real> xTilde ;
    public:

      // Define input and output.
      RoeAveragedScalar() {
        name_store("rhoTilde",rhoTilde) ;
        name_store("lefts(X)",xLeft) ;
        name_store("rights(X)",xRight) ;
        name_store("roeAveragedScalar(X)",xTilde) ;
        input("rhoTilde,lefts(X),rights(X)") ;
        output("roeAveragedScalar(X)") ;
        constraint("internalFaces") ;
      }

      // Limited extrapolation from the cell to the face.
      void calculate(Entity face) {
        xTilde[face]=(xLeft[face]+xRight[face]*rhoTilde[face])/(1.0+
          rhoTilde[face]) ;
      }

      // Call calculate for a sequence of faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<RoeAveragedScalar> registerRoeAveragedScalar ;

  // Rule for computing Roe-averaged vector values at the faces.
  class RoeAveragedVector : public pointwise_rule {
    private:
      const_store<real> rhoTilde ;
      const_store<vect3d> xLeft ;
      const_store<vect3d> xRight ;
      store<vect3d> xTilde ;
    public:

      // Define input and output.
      RoeAveragedVector() {
        name_store("rhoTilde",rhoTilde) ;
        name_store("leftv3d(X)",xLeft) ;
        name_store("rightv3d(X)",xRight) ;
        name_store("roeAveragedVector(X)",xTilde) ;
        input("rhoTilde,leftv3d(X),rightv3d(X)") ;
        output("roeAveragedVector(X)") ;
        constraint("internalFaces") ;
      }

      // Limited extrapolation from the cell to the face.
      void calculate(Entity face) {
        xTilde[face]=(xLeft[face]+xRight[face]*rhoTilde[face])/(1.0+
          rhoTilde[face]) ;
      }

      // Call calculate for a sequence of faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<RoeAveragedVector> registerRoeAveragedVector ;

//-----------------------------------------------------------------------------
// Rules for computing the various variables for the Roe scheme.

  // Rule for computing a number of variables. 
  class ComputeRoeVariables : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<real> soundSpeed ;
      const_store<real> rhoTilde ;
      const_store<real> pLeft,pRight ;
      const_store<vect3d> vLeft,vRight,vTilde ;
      const_store<Area> area ;
      store<real> vNormalTilde,deltaVNormal,delVNormal ;
      store<real> deltaP,delP,cStar ;
    public:

      // Define input and output.
      ComputeRoeVariables() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("soundSpeed",soundSpeed) ;
        name_store("rhoTilde",rhoTilde) ;
        name_store("lefts(p)",pLeft) ;
        name_store("rights(p)",pRight) ;
        name_store("leftv3d(v)",vLeft) ;
        name_store("rightv3d(v)",vRight) ;
        name_store("roeAveragedVector(v)",vTilde) ;
        name_store("area",area) ;
        name_store("vNormalTilde",vNormalTilde) ;
        name_store("deltaVNormal",deltaVNormal) ;
        name_store("delVNormal",delVNormal) ;
        name_store("deltaP",deltaP) ;
        name_store("delP",delP) ;
        name_store("cStar",cStar) ;
        input("(cl,cr)->soundSpeed") ;
        input("rhoTilde,lefts(p),rights(p),leftv3d(v),rightv3d(v)") ;
        input("roeAveragedVector(v),area") ;
        output("vNormalTilde,deltaVNormal,delVNormal") ;
        output("delP,deltaP,cStar") ;
        constraint("internalFaces") ;
      }

      // Compute all the Roe variables for a face.
      void calculate(Entity face) {
        real cFace=0.5*(soundSpeed[cl[face]]+soundSpeed[cr[face]]) ;
        vNormalTilde[face]=dot(vTilde[face],area[face].n) ;
        deltaVNormal[face]=dot(vRight[face]-vLeft[face],area[face].n) ;
        cStar[face]=0.5*(abs(vNormalTilde[face]+cFace)+abs(vNormalTilde[face]-
          cFace)) ;
        real machNumberStar=0.5*(abs(vNormalTilde[face]+cFace)-
          abs(vNormalTilde[face]-cFace))/cFace ;
        deltaP[face]=pRight[face]-pLeft[face] ;
        delVNormal[face]=machNumberStar*deltaVNormal[face]+(cStar[face]-
          abs(vNormalTilde[face]))*deltaP[face]/(rhoTilde[face]*cFace*cFace) ;
        delP[face]=machNumberStar*deltaP[face]+(cStar[face]-abs(vNormalTilde
          [face]))*rhoTilde[face]*delVNormal[face] ;
      }

      // Call calculate for a sequence of faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeRoeVariables> registerComputeRoeVariables ;

  // Rule for computing the Roe mass flux. 
  class ComputeRoeMassFlux : public pointwise_rule {
    private:
      const_store<real> faceDensity ;
      const_store<real> rhoTilde ;
      const_store<vect3d> vLeft,vRight ;
      const_store<real> delVNormal ;
      const_store<Area> area ;
      store<real> roeMassFlux ;
    public:

      // Define input and output.
      ComputeRoeMassFlux() {
        name_store("faceDensity",faceDensity) ;
        name_store("rhoTilde",rhoTilde) ;
        name_store("leftv3d(v)",vLeft) ;
        name_store("rightv3d(v)",vRight) ;
        name_store("delVNormal",delVNormal) ;
        name_store("area",area) ;
        name_store("roeMassFlux",roeMassFlux) ;
        input("faceDensity,rhoTilde,leftv3d(v),rightv3d(v)") ;
        input("delVNormal,area") ;
        output("roeMassFlux") ;
        constraint("internalFaces,roeInviscidFlux") ;
      }

      // Compute all the Roe variables for a face.
      void calculate(Entity face) {
        roeMassFlux[face]=0.5*area[face].sada*(faceDensity[face]*
          dot(vLeft[face]+vRight[face],area[face].n)-rhoTilde[face]*
          delVNormal[face]) ;
      }

      // Call calculate for a sequence of faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeRoeMassFlux> registerComputeRoeMassFlux ;

}
