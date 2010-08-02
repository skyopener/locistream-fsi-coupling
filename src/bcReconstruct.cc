// Loci includes.
#include <Loci.h>
#include <Tools/tools.h>
using Loci::Area ;

// StreamUns includes.
#include "const.h"
#include "sciTypes.h"

namespace streamUns {

  class boundary_reconstruct : public pointwise_rule {
    const_Map ci,cl,cr ;
    const_multiMap upper,lower ;
    const_store<real> X ;
    const_store<vect3d> cellcenter,facecenter ;
    const_store<Area> area ;
    store<real> boundaryX ;
  public:
    boundary_reconstruct() {
      name_store("ci",ci) ;
      name_store("cl",cl) ;
      name_store("cr",cr) ;
      name_store("upper",upper) ;
      name_store("lower",lower) ;
      name_store("X",X) ;
      name_store("cellcenter",cellcenter) ;
      name_store("facecenter",facecenter) ;
      name_store("area",area) ;
      name_store("boundaryValue(X)",boundaryX) ;
      input("ci->upper->cr->(X,cellcenter)") ;
      input("ci->lower->cl->(X,cellcenter)") ;
      input("ci->(X,cellcenter)") ;
      input("area,facecenter") ;
      output("boundaryValue(X)") ;
    }

    void calculate(Entity fc) {
      Entity c = ci[fc] ;
      const vect3d n = area[fc].n ;
      real sx2i = 0 ;
      real suxi = 0 ;
      double maxX= X[c] ;
      double minX= X[c] ;
      // Allow a 20% change 
      minX = min(minX,1.2*X[c]) ;
      minX = min(minX,.8*X[c]) ;
      maxX = max(maxX,1.2*X[c]) ;
      maxX = max(maxX,.8*X[c]) ;
      for(const Entity *lp = lower.begin(c);lp!=lower.end(c);++lp) {
        const Entity cc = cl[*lp] ;
        const real x = -dot((cellcenter[cc]-cellcenter[c]),n) ;
        minX =min(minX,X[cc]) ;
        maxX =max(maxX,X[cc]) ;
        if(x > 0) {
          real u = X[cc]-X[c] ;
          sx2i += x*x ;
          suxi += u*x ;
        }
      }
      for(const Entity *up = upper.begin(c);up!=upper.end(c);++up) {
        const Entity cc = cr[*up] ;
        real x = -dot((cellcenter[cc]-cellcenter[c]),n) ;
        if(x > 0) {
          minX = min(minX,X[cc]) ;
          maxX = max(maxX,X[cc]) ;
          real u = X[cc]-X[c] ;
          sx2i += x*x ;
          suxi += u*x ;
        }
      }
      real g = suxi/(sx2i+EPSILON) ;
      real xb = -dot((facecenter[fc]-cellcenter[c]),n) ;
      boundaryX[fc] = X[c] + g*xb ;

      /* Limit to first order if change too much
      if(boundaryX[fc] < minX)
        boundaryX[fc] = X[c] ;
      if(boundaryX[fc] > maxX)
        boundaryX[fc] = X[c] ;*/
    }
    void compute(const sequence &seq) {
      do_loop(seq,this) ;
    }
  } ;

  register_rule<boundary_reconstruct> register_boundary_reconstruct ;

  class boundary_reconstruct_v3d : public pointwise_rule {
    const_Map ci,cl,cr ;
    const_multiMap upper,lower ;
    const_store<vect3d> X ;
    const_store<vect3d> cellcenter,facecenter ;
    const_store<Area> area ;
    store<vect3d> boundaryX ;
  public:
    boundary_reconstruct_v3d() {
      name_store("ci",ci) ;
      name_store("cl",cl) ;
      name_store("cr",cr) ;
      name_store("upper",upper) ;
      name_store("lower",lower) ;
      name_store("X",X) ;
      name_store("cellcenter",cellcenter) ;
      name_store("facecenter",facecenter) ;
      name_store("area",area) ;
      name_store("boundaryValueV3D(X)",boundaryX) ;
      input("ci->upper->cr->(X,cellcenter)") ;
      input("ci->lower->cl->(X,cellcenter)") ;
      input("ci->(X,cellcenter)") ;
      input("area,facecenter") ;
      output("boundaryValueV3D(X)") ;
    }

    void calculate(Entity fc) {
      const Entity c = ci[fc] ;
      const vect3d n = area[fc].n ;
      real sx2i = 0 ;
      vect3d suxi = vect3d(0.0,0.0,0.0) ;

      double maxU = norm(X[c]) ;
      double minU = maxU ;

      maxU = maxU*1.2 ;
      minU = minU*.8 ;
      
      for(const Entity *lp = lower.begin(c);lp!=lower.end(c);++lp) {
        const Entity cc = cl[*lp] ;
        const real x = -dot((cellcenter[cc]-cellcenter[c]),n) ;
        if(x > 0) {
          const vect3d u = X[cc]-X[c] ;
          sx2i += x*x ;
          suxi += u*x ;
        }
      }
      for(const Entity *up = upper.begin(c);up!=upper.end(c);++up) {
        const Entity cc = cr[*up] ;
        const real x = -dot((cellcenter[cc]-cellcenter[c]),n) ;
        if(x > 0) {
          const vect3d u = X[cc]-X[c] ;
          sx2i += x*x ;
          suxi += u*x ;
        }
      }
      const vect3d g = suxi/(sx2i+EPSILON) ;
      const real xb = -dot((facecenter[fc]-cellcenter[c]),n) ;
      boundaryX[fc] = X[c] + g*xb ;
      /*double nb = norm(boundaryX[fc]) ;
      If things change too mcuh, switch to first order
      if(nb > maxU || nb < minU)
        boundaryX[fc] = X[c] ;*/
    }
    void compute(const sequence &seq) {
      do_loop(seq,this) ;
    }
  } ;

  register_rule<boundary_reconstruct_v3d> register_boundary_reconstruct_v3d ;

}
