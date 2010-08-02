// Standard library includes.
#include <vector>
using std::vector ;
                                                                                
// Loci includes.
#include <Loci.h>
using Loci::Area ;
                                                                                
// StreamUns includes.
#include "gridReader/readGrid.h"
#include "referenceFrame.h"
#include "sciTypes.h"
                                                                                
namespace streamUns {

  // Adds turbomachinery options to noslip boundaries.
  class CheckNoslipTurbo : public BC_Check {
    private:
      string errorMessage ;
    public:
      CheckNoslipTurbo() : errorMessage("") {}
    public:
      string BoundaryConditions() { return "noslip" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        return true ;
      }
      ostream &ErrorMessage(std::ostream &s) {
        s << errorMessage << endl ; return s;
      }
      string VariablesChecked(fact_db &facts) {
        string s="stationary" ; return s ;
      }
  } ;
                                                                                
  register_BC<CheckNoslipTurbo> registerCheckNoslipTurbo ;

  // Rule for boundary faces with specified velocity. This priority rule
  // automatically adds on the -cross(omega,r) component so that the use need
  // only specify the velocity from the absolute frame of reference.
  class BoundaryVelocitySpecificationTurbo : public pointwise_rule {
    private:
      const_Map ci,ref ;
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<unsigned int> cellReferenceFrame ;
      const_store<vect3d> faceCenter ;
      const_store<vect3d> constantV_BC ;
      store<vect3d> v_f ;
    public:
                                                                                
      // Define input and output.
      BoundaryVelocitySpecificationTurbo() {
        name_store("ci",ci) ;
        name_store("ref",ref) ;
        name_store("referenceFrame",referenceFrame) ;
        name_store("cellReferenceFrame",cellReferenceFrame) ;
        name_store("facecenter",faceCenter) ;
        name_store("constantV_BC",constantV_BC) ;
        name_store("turbo::v_f",v_f) ;
        input("ref->constantV_BC,referenceFrame,ci->cellReferenceFrame") ;
        input("facecenter") ;
        output("turbo::v_f") ;
      }
                                                                                
      // Calculate face velocity for a single face.
      void calculate(Entity face) {
        unsigned int n=cellReferenceFrame[ci[face]] ;
        vect3d delta=(*referenceFrame)[n].axisEnd-(*referenceFrame)[n].
          axisStart,omega=((*referenceFrame)[n].omega/norm(delta))*delta,
          r=faceCenter[face]-(*referenceFrame)[n].axisStart ;
        v_f[face]=constantV_BC[ref[face]]-cross(omega,r) ;
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;
                                                                                
  register_rule<BoundaryVelocitySpecificationTurbo>
    registerBoundaryVelocitySpecificationTurbo ;

  // Priority rule for assigning velocity on no-slip boundary faces that are
  // stationary. Since the equations for turbomachinery applications are always
  // solved in the relative frame, surfaces rotating with the frame of
  // reference have zero velocity and stationary no-slip boundaries have a
  // velocity of -cross(omega,r) . The reference frame of a boundary face is
  // assumed to be the same as that of its boundary cell.
  class BoundaryVelocityNoSlipStationary : public pointwise_rule {
    private:
      const_Map ci ;
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<unsigned int> cellReferenceFrame ;
      const_store<vect3d> faceCenter ;
      store<vect3d> v_f ;
    public:
                                                                                
      // Define input and output.
      BoundaryVelocityNoSlipStationary() {
        name_store("ci",ci) ;
        name_store("referenceFrame",referenceFrame) ;
        name_store("cellReferenceFrame",cellReferenceFrame) ;
        name_store("facecenter",faceCenter) ;
        name_store("stationary::v_f",v_f) ;
        input("referenceFrame,ci->cellReferenceFrame,facecenter") ;
        output("stationary::v_f") ;
        constraint("noslip_BC,ref->stationary_BCoption,turbomachinery") ;
      }
                                                                                
      // Calculate face velocity for a single face.
      void calculate(Entity face) {
        unsigned int n=cellReferenceFrame[ci[face]] ;
        vect3d delta=(*referenceFrame)[n].axisEnd-(*referenceFrame)[n].
          axisStart,omega=((*referenceFrame)[n].omega/norm(delta))*delta,
          r=faceCenter[face]-(*referenceFrame)[n].axisStart ;
        v_f[face]=-1.0*cross(omega,r) ;
      }
                                                                                
      // Calculate face velocity for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<BoundaryVelocityNoSlipStationary>
    registerBoundaryVelocityNoSlipStationary ;

  // Rule for for assigning velocity on inlet faces with specified mass flux.
  class BoundaryVelocityFromMassFluxNormalTurbo : public pointwise_rule {
    private:
      const_Map ci,ref ;
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<unsigned int> cellReferenceFrame ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux_BC ;
      const_store<real> rho_f ;
      const_store<Area> area ;
      store<vect3d> v_f ;
    public:

      // Define input and output.
      BoundaryVelocityFromMassFluxNormalTurbo() {
        name_store("ci",ci) ;
        name_store("ref",ref) ;
        name_store("referenceFrame",referenceFrame) ;
        name_store("cellReferenceFrame",cellReferenceFrame) ;
        name_store("facecenter",faceCenter) ;
        name_store("massFlux_BC",massFlux_BC) ;
        name_store("rho_f",rho_f) ;
        name_store("area",area) ;
        name_store("turbo::v_f",v_f) ;
        input("ref->massFlux_BC,rho_f,area") ;
        input("referenceFrame,ci->cellReferenceFrame,facecenter") ;
        output("turbo::v_f") ;
      }

      // Calculate face velocity for a single face.
      void calculate(Entity face) {
        unsigned int n=cellReferenceFrame[ci[face]] ;
        vect3d delta=(*referenceFrame)[n].axisEnd-(*referenceFrame)[n].
          axisStart,omega=((*referenceFrame)[n].omega/norm(delta))*delta,
          r=faceCenter[face]-(*referenceFrame)[n].axisStart ;
        v_f[face]=(massFlux_BC[ref[face]]/rho_f[face])*area[face].n-
          cross(omega,r) ;
      }

      // Calculate face velocity for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryVelocityFromMassFluxNormalTurbo>
    registerBoundaryVelocityFromMassFluxNormalTurbo ;

  // Rule for for assigning velocity on inlet faces with specified mass flux
  // and flow direction.
  class BoundaryVelocityFromMassFluxTurbo : public pointwise_rule {
    private:
      const_Map ci,ref ;
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<unsigned int> cellReferenceFrame ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux_BC ;
      const_store<vect3d> flowDirection_BC ;
      const_store<real> rho_f ;
      store<vect3d> v_f ;
    public:

      // Define input and output.
      BoundaryVelocityFromMassFluxTurbo() {
        name_store("ci",ci) ;
        name_store("ref",ref) ;
        name_store("referenceFrame",referenceFrame) ;
        name_store("cellReferenceFrame",cellReferenceFrame) ;
        name_store("facecenter",faceCenter) ;
        name_store("massFlux_BC",massFlux_BC) ;
        name_store("flowDirection_BC",flowDirection_BC) ;
        name_store("rho_f",rho_f) ;
        name_store("turbo::specifiedFlowDirection::v_f",v_f) ;
        input("ref->(massFlux_BC,flowDirection_BC),rho_f") ;
        input("referenceFrame,ci->cellReferenceFrame,facecenter") ;
        output("turbo::specifiedFlowDirection::v_f") ;
      }

      // Calculate face velocity for a single face.
      void calculate(Entity face) {
        unsigned int n=cellReferenceFrame[ci[face]] ;
        vect3d delta=(*referenceFrame)[n].axisEnd-(*referenceFrame)[n].
          axisStart,omega=((*referenceFrame)[n].omega/norm(delta))*delta,
          r=faceCenter[face]-(*referenceFrame)[n].axisStart ;
        v_f[face]=(-massFlux_BC[ref[face]]/rho_f[face])*
          flowDirection_BC[ref[face]]-cross(omega,r) ;
      }

      // Calculate face velocity for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryVelocityFromMassFluxTurbo>
    registerBoundaryVelocityFromMassFluxTurbo ;

}
