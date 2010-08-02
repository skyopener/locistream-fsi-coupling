//-----------------------------------------------------------------------------
// Description: This file contains rules for outputting wall values to a file
//   which can be post-processed with Ed's Fieldview extractor.
//
// Authors: Original code by Ed Luke. Modified by Jeff Wright
//-----------------------------------------------------------------------------
                                                                                
// Loci includes.
#include <Loci.h>
using Loci::Area ;
                                                                                
// StreamUns includes.
#include "sciTypes.h"
                                                                                
namespace streamUns {

  // Rule to compute the wall shear stress for cases in which wall functions
  // are not used.
  class TauNoWallFunction : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<vect3d> v ;
      const_store<real> laminarViscosity ;
      const_store<vect3d> faceCenter,cellCenter ;
      const_store<Area> area ;
      store<vect3d> tauWall ;
    public:

      TauNoWallFunction() {
        name_store("ci",ci) ;
        name_store("v",v) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("facecenter",faceCenter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("area",area) ;
        name_store("tauWall",tauWall) ;
        input("ci->(v,laminarViscosity,cellcenter),facecenter,area") ;
        output("tauWall") ;
      }

      // Compute the wall shear stress for a single face.
      void calculate (Entity face) {
        const real vNormal= dot(v[ci[face]],area[face].n) ;
        const vect3d vTangential=v[ci[face]]-vNormal*area[face].n ;
        const real yNormal=dot(faceCenter[face]-cellCenter[ci[face]],
          area[face].n) ;
        tauWall[face]=vTangential*(-laminarViscosity[ci[face]]/yNormal) ;
      }

      // Compute the wall shear stress for all faces
      virtual void compute (const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TauNoWallFunction> registerTauNoWallFunction ;

  // Priority rule to compute tau on wall-function boundaries.
  class TauWallFunction : public pointwise_rule {
    private:
      const_store<vect3d> tauWallTemp ;
      store<vect3d> tauWall ;
    public:

      TauWallFunction() {
        name_store("tauWallTemp",tauWallTemp) ;
        name_store("wallFunction::tauWall",tauWall) ;
        input("tauWallTemp") ;
        output("wallFunction::tauWall") ;
        constraint("noslip_BC,ref->wallFunction_BCoption") ;
      }

      // Copy value.
      void calculate (Entity face) { tauWall[face]=tauWallTemp[face] ; }

      // Loop over faces.
      virtual void compute (const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TauWallFunction> registerTauWallFunction ;

  // Rule to compute y-plus at the wall for cases in which wall functions
  // are not used. Note that we are computing tauWall here as well, since
  // using the rule above to provide it as input creates a cycle somehow.
  // Point this out to Ed, since I don't think this should be the case.
  class YPlusNoWallFunction : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> rho ;
      const_store<real> laminarViscosity ;
      const_store<vect3d> faceCenter,cellCenter ;
      const_store<Area> area ;
      const_store<vect3d> tauWall ;
      store<real> yPlusWall ;
    public:

      YPlusNoWallFunction() {
        name_store("ci",ci) ;
        name_store("rho",rho) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("facecenter",faceCenter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("area",area) ;
        name_store("tauWall",tauWall) ;
        name_store("yPlusWall",yPlusWall) ;
        input("ci->(rho,laminarViscosity,cellcenter)") ;
        input("facecenter,tauWall,area") ;
        output("yPlusWall") ;
        constraint("noslip_BC") ;
      }

      // Compute y-plus for a single face.
      void calculate (Entity face) {
        const real tauWallNorm=norm(tauWall[face]) ;
        const real yNormal=dot(faceCenter[face]-cellCenter[ci[face]],
          area[face].n) ;
        yPlusWall[face]=sqrt(rho[ci[face]]*tauWallNorm)*yNormal/
          laminarViscosity[ci[face]] ;
      }

      // Loop over faces.
      virtual void compute (const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<YPlusNoWallFunction> registerYPlusNoWallFunction ;

  // Priority rule to compute yPlus on wall-function boundaries.
  class YPlusWallFunction : public pointwise_rule {
    private:
      const_store<real> yPlusWallTemp ;
      store<real> yPlusWall ;
    public:

      YPlusWallFunction() {
        name_store("yPlusWallTemp",yPlusWallTemp) ;
        name_store("wallFunction::yPlusWall",yPlusWall) ;
        input("yPlusWallTemp") ;
        output("wallFunction::yPlusWall") ;
        constraint("noslip_BC,ref->wallFunction_BCoption") ;
      }

      // Copy value.
      void calculate (Entity face) {
        yPlusWall[face]=yPlusWallTemp[face] ;
      }

      // Loop over faces.
      virtual void compute (const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<YPlusWallFunction> registerYPlusWallFunction ;

}


