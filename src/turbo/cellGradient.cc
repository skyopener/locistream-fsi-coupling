//-----------------------------------------------------------------------------
// Description: This file contains a modified version of gradv3d which takes
//   into account reference frame changes of neighbor cells.
//-----------------------------------------------------------------------------

// Standard library includes.
#include <vector>
using std::vector ;

// Loci includes.
#include <Loci.h>
using Loci::real_t ;

// StreamUns includes.
#include "referenceFrame.h"
#include "sciTypes.h"
                                                                                
namespace streamUns {

  // Note that all boundary values used to compute a cell gradient are
  // assumed to be in the same reference frame as the cell. Reference frame
  // change should occur across only internal faces. This function should
  // only be used for total enthalpy.
  class GradScalarTurbo : public pointwise_rule {
    private:
      const_multiMap lower, upper,boundary_map ;
      const_Map cl,cr ;
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<unsigned int> cellReferenceFrame ;
      const_store<real> X,X_f ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> Wf_l, Wf_r ;
      store<vect3d> grads ;
      int vs ;

    public:

      // Define input and output.
      GradScalarTurbo() {
        name_store("lower",lower) ;
        name_store("upper",upper) ;
        name_store("boundary_map",boundary_map) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("referenceFrame",referenceFrame) ;
        name_store("cellReferenceFrame",cellReferenceFrame) ;
        name_store("X",X) ;
        name_store("X_f",X_f) ;
        name_store("cellcenter",cellCenter) ;
        name_store("gradsTurbo(X)",grads) ;
        name_store("Wf_l",Wf_l) ;
        name_store("Wf_r",Wf_r) ;
        input("lower->Wf_r") ;
        input("upper->Wf_l") ;
        input("boundary_map->Wf_l") ;
        input("lower->cl->(X,cellReferenceFrame,cellcenter)") ;
        input("upper->cr->(X,cellReferenceFrame,cellcenter)") ;
        input("boundary_map->X_f") ;
        input("X,cellReferenceFrame,cellcenter") ;
        output("gradsTurbo(X)") ;
        constraint("geom_cells,stencilStandard,turbomachinery") ;
      }

      // Calculate gradient for a cell.
      void calculate(Entity cell) {

        vect3d Qt_b = vect3d(0,0,0) ;
        real X_center = X[cell] ;
        for(const int *li=lower.begin(cell);li!=lower.end(cell);++li) {
          unsigned int n=cellReferenceFrame[cell] ;
          unsigned int nNbr=cellReferenceFrame[cl[*li]] ;
          if(n==nNbr){
            real df = X[cl[*li]] - X_center ;
            Qt_b += df*Wf_r[*li] ;
          }else{
            vect3d delta=(*referenceFrame)[n].axisEnd-(*referenceFrame)[n].
              axisStart,omega=((*referenceFrame)[n].omega/norm(delta))*delta,
              r=cellCenter[cl[*li]]-(*referenceFrame)[n].axisStart ;
            vect3d deltaNbr=(*referenceFrame)[nNbr].axisEnd-(*referenceFrame)
              [nNbr].axisStart,omegaNbr=((*referenceFrame)[nNbr].omega/
              norm(delta))*delta,rNbr=cellCenter[cl[*li]]-(*referenceFrame)
              [nNbr].axisStart ;
            const vector3d<real_t> a=cross(omegaNbr,rNbr),b=cross(omega,r) ;
            const real df=X[cl[*li]]+0.5*dot(a,a)-0.5*dot(b,b)-X_center ;
            Qt_b += df*Wf_r[*li] ;
          }
        }
        for(const int *ui=upper.begin(cell);ui!=upper.end(cell);++ui) {
          unsigned int n=cellReferenceFrame[cell] ;
          unsigned int nNbr=cellReferenceFrame[cr[*ui]] ;
          if(n==nNbr){
            real df = X[cr[*ui]] - X_center ;
            Qt_b += df*Wf_l[*ui] ;
          }else{
            vect3d delta=(*referenceFrame)[n].axisEnd-(*referenceFrame)[n].
              axisStart,omega=((*referenceFrame)[n].omega/norm(delta))*delta,
              r=cellCenter[cr[*ui]]-(*referenceFrame)[n].axisStart ;
            vect3d deltaNbr=(*referenceFrame)[nNbr].axisEnd-(*referenceFrame)
              [nNbr].axisStart,omegaNbr=((*referenceFrame)[nNbr].omega/
              norm(delta))*delta,rNbr=cellCenter[cr[*ui]]-(*referenceFrame)
              [nNbr].axisStart ;
            const vector3d<real_t> a=cross(omegaNbr,rNbr),b=cross(omega,r) ;
            const real df=X[cr[*ui]]+0.5*dot(a,a)-0.5*dot(b,b)-X_center ;
            Qt_b += df*Wf_l[*ui] ;
          }
        }
        for(const int *bi=boundary_map.begin(cell);bi!=boundary_map.end(cell);
        ++bi) {
          real df = X_f[*bi] - X_center ;
          Qt_b += df*Wf_l[*bi] ;
        }

        grads[cell] = Qt_b ;

      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GradScalarTurbo> registerGradScalarTurbo ;

  // Note that all boundary values used to compute a cell gradient are
  // assumed to be in the same reference frame as the cell. Reference frame
  // change should occur across only internal faces.
  class GradVec3dTurbo : public pointwise_rule {
    private:
      const_multiMap lower,upper,boundary_map ;
      const_Map cl,cr ;
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<unsigned int> cellReferenceFrame ;
      const_store<vector3d<real_t> > X,X_f,cellCenter ;
      const_store<vector3d<real_t> > Wf_l,Wf_r ;
      store<tensor3d<real_t> > grads ;
    public:

      // Define input and output.
      GradVec3dTurbo() {
        name_store("lower",lower) ;
        name_store("upper",upper) ;
        name_store("boundary_map",boundary_map) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("referenceFrame",referenceFrame) ;
        name_store("cellReferenceFrame",cellReferenceFrame) ;
        name_store("X",X) ;
        name_store("X_f",X_f) ;
        name_store("cellcenter",cellCenter) ;
        name_store("turbo::gradv3d(X)",grads) ;
        name_store("Wf_l",Wf_l) ;
        name_store("Wf_r",Wf_r) ;
        input("lower->Wf_r") ;
        input("upper->Wf_l") ;
        input("boundary_map->Wf_l") ;
        input("lower->cl->(X,cellReferenceFrame,cellcenter)") ;
        input("upper->cr->(X,cellReferenceFrame,cellcenter)") ;
        input("boundary_map->X_f") ;
        input("X,cellReferenceFrame,cellcenter") ;
        output("turbo::gradv3d(X)") ;
        constraint("geom_cells,stencilStandard,turbomachinery") ;
      }

      // Calculate gradient for a cell.
      void calculate(Entity cell) {
        tensor3d<real_t> Qt_b=tensor3d<real_t>(vector3d<real_t>(0,0,0),
          vector3d<real_t>(0,0,0),vector3d<real_t>(0,0,0)) ;
        const vector3d<real_t> X_center=X[cell] ;
        for(const int *li=lower.begin(cell);li!=lower.end(cell);++li) {
          unsigned int n=cellReferenceFrame[cell] ;
          unsigned int nNbr=cellReferenceFrame[cl[*li]] ;
          if(n==nNbr){
            const vector3d<real_t> df=X[cl[*li]] - X_center ;
            Qt_b.x += df.x*Wf_r[*li] ; Qt_b.y += df.y*Wf_r[*li] ;
            Qt_b.z += df.z*Wf_r[*li] ;
          }else{
            vect3d delta=(*referenceFrame)[n].axisEnd-(*referenceFrame)[n].
              axisStart,omega=((*referenceFrame)[n].omega/norm(delta))*delta,
              r=cellCenter[cl[*li]]-(*referenceFrame)[n].axisStart ;
            vect3d deltaNbr=(*referenceFrame)[nNbr].axisEnd-(*referenceFrame)
              [nNbr].axisStart,omegaNbr=((*referenceFrame)[nNbr].omega/
              norm(delta))*delta,rNbr=cellCenter[cl[*li]]-(*referenceFrame)
              [nNbr].axisStart ;
            const vector3d<real_t> df=X[cl[*li]]+cross(omegaNbr,rNbr)-
              cross(omega,r)-X_center ;
            Qt_b.x+=df.x*Wf_r[*li] ; Qt_b.y+=df.y*Wf_r[*li] ;
            Qt_b.z+=df.z*Wf_r[*li] ;
          }
        }
        for(const int *ui=upper.begin(cell);ui!=upper.end(cell);++ui) {
          unsigned int n=cellReferenceFrame[cell] ;
          unsigned int nNbr=cellReferenceFrame[cr[*ui]] ;
          if(n==nNbr){
            const vector3d<real_t> df=X[cr[*ui]]-X_center ;
            Qt_b.x+=df.x*Wf_l[*ui] ; Qt_b.y+=df.y*Wf_l[*ui] ;
            Qt_b.z+=df.z*Wf_l[*ui] ;
          }else{
            vect3d delta=(*referenceFrame)[n].axisEnd-(*referenceFrame)[n].
              axisStart,omega=((*referenceFrame)[n].omega/norm(delta))*delta,
              r=cellCenter[cr[*ui]]-(*referenceFrame)[n].axisStart ;
            vect3d deltaNbr=(*referenceFrame)[nNbr].axisEnd-(*referenceFrame)
              [nNbr].axisStart,omegaNbr=((*referenceFrame)[nNbr].omega/
              norm(delta))*delta,rNbr=cellCenter[cr[*ui]]-(*referenceFrame)
              [nNbr].axisStart ;
            const vector3d<real_t> df=X[cr[*ui]]+cross(omegaNbr,rNbr)-
              cross(omega,r)-X_center ;
            Qt_b.x+=df.x*Wf_l[*ui] ; Qt_b.y+=df.y*Wf_l[*ui] ;
            Qt_b.z+=df.z*Wf_l[*ui] ;
          }
        }
        for(const int *bi = boundary_map.begin(cell);bi!=boundary_map.
        end(cell);++bi) {
          const vector3d<real_t> df=X_f[*bi]-X_center ;
          Qt_b.x+=df.x*Wf_l[*bi] ; Qt_b.y+=df.y*Wf_l[*bi] ;
          Qt_b.z+=df.z*Wf_l[*bi] ;
        }
        grads[cell]=Qt_b ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GradVec3dTurbo> registerGradVec3dTurbo ;

}
