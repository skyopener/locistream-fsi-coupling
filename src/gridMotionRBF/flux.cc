//-----------------------------------------------------------------------------
// Description: This file contains rules for the GCL mass fluxes required for
//   grid movement.
//-----------------------------------------------------------------------------

// Standard library includes.
#include <vector>
using std::vector ;

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// StreamUns includes.
#include "sciTypes.h"

namespace streamUns {

  // Add inlets to inletOutlet_BC constraint. This contraint is used in
  // setting the grid mass flux at inlet and outlet boundaries. All other
  // boundaries have zero grid mass flux.
  class AddInletsToInletOutletBC : public pointwise_rule {
    private:
      store<bool> inletOutlet_BC ;
    public:

      // Define input and output.
      AddInletsToInletOutletBC() {
        name_store("inletOutlet_BC",inletOutlet_BC) ;
        output("inletOutlet_BC") ;
        constraint("inlet_BC") ;
      }

      // Set the flag to true.
      void calculate(Entity face) { inletOutlet_BC[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<AddInletsToInletOutletBC> registerAddInletsToInletOutletBC ;

  // Add outlets to inletOutlet_BC constraint. This contraint is used in
  // setting the grid mass flux at inlet and outlet boundaries. All other
  // boundaries have zero grid mass flux.
  class AddOutletsToInletOutletBC : public pointwise_rule {
    private:
      store<bool> inletOutlet_BC ;
    public:

      // Define input and output.
      AddOutletsToInletOutletBC() {
        name_store("inletOutlet_BC",inletOutlet_BC) ;
        output("inletOutlet_BC") ;
        constraint("outlet_BC") ;
      }

      // Set the flag to true.
      void calculate(Entity face) { inletOutlet_BC[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<AddOutletsToInletOutletBC> registerAddOutletsToInletOutletBC ;

  // Rule to compute the facet normal sum and sum of the dot product of facet
  // velocity with the facet normal for each face.
  class GCLOldFacetSums : public pointwise_rule {
    private:
      const_multiMap face2node ;
      const_param<real> timeStep ;
      const_store<vect3d> node_vOld,posOld,pos ;
      const_store<vect3d> faceCenter,oldFaceCenter ;
      store<vect3d> oldFacetNormalSum ;
      store<real> oldFacetVDotNormalSum ;
    public:

      // Define input and output.
      GCLOldFacetSums() {
        name_store("face2node",face2node) ;
        name_store("timeStep{n}",timeStep) ;
        name_store("node_vOld{n}",node_vOld) ;
        name_store("pos{n-1}",posOld) ;
        name_store("pos{n}",pos) ;
        name_store("facecenter{n}",faceCenter) ;
        name_store("facecenter{n-1}",oldFaceCenter) ;
        name_store("oldFacetNormalSum{n}",oldFacetNormalSum) ;
        name_store("oldFacetVDotNormalSum{n}",oldFacetVDotNormalSum) ;
        input("timeStep{n},face2node->(node_vOld{n},pos{n-1},pos{n})") ;
        input("facecenter{n},facecenter{n-1}") ;
        output("oldFacetNormalSum{n},oldFacetVDotNormalSum{n}") ;
      }

      // Calculate the mass flux for a single face.
      void calculate(Entity face) {

        // Save face centers.
        const vect3d center=faceCenter[face],oldCenter=oldFaceCenter[face] ;

        // Compute face center velocity. Fixed the BUG where the face center
        // velocity was computed incorrectly by averaging the face nodal
        // velocities. Must use the actual old and new face centers since the
        // coordinates of the face center are obtained elsewhere, NOT by using
        // the average of the face node coordinates.
        const vect3d face_vOld=(1.0/(*timeStep))*(center-oldCenter) ;

        // Face must be broken into its component facets to be GCL compliant.
        // Add the physical velocity and grid velocity components to the mass
        // flux on a facet-by-facet basis. Note that we need to use the exact
        // velocity of each facet center, but can still just use the old face
        // fluid velocity averaged from the centers.
        int n0=*(face2node.begin(face)),n1=*(face2node.end(face)-1) ;
        vect3d facetNormal=(cross(posOld[n1]-oldCenter,posOld[n0]-oldCenter)+
          cross(pos[n1]-center,pos[n0]-center)+0.5*(cross(posOld[n1]-oldCenter,
          pos[n0]-center)+cross(pos[n1]-center,posOld[n0]-oldCenter)))/6.0 ;
        vect3d facetVelocity=(node_vOld[n0]+node_vOld[n1]+face_vOld)/3.0 ;
        oldFacetNormalSum[face]=facetNormal ;
        oldFacetVDotNormalSum[face]=dot(facetVelocity,facetNormal) ;
        for(const int *ni=face2node.begin(face)+1;ni!=(face2node.end(face));
        ni++) {
          n0=*(ni),n1=*(ni-1) ;
          facetNormal=(cross(posOld[n1]-oldCenter,posOld[n0]-oldCenter)+
            cross(pos[n1]-center,pos[n0]-center)+0.5*(cross(posOld[n1]-
            oldCenter,pos[n0]-center)+cross(pos[n1]-center,posOld[n0]-
            oldCenter)))/6.0 ;
          facetVelocity=(node_vOld[n0]+node_vOld[n1]+face_vOld)/3.0 ;
          oldFacetNormalSum[face]+=facetNormal ;
          oldFacetVDotNormalSum[face]+=dot(facetVelocity,facetNormal) ;
        }
      }

      // Calculate the mass flux for all faces in a sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<GCLOldFacetSums> registerGCLOldFacetSums ;

  // Rule to compute the facet normal sum and sum of the dot product of facet
  // velocity with the facet normal for each face.
  class GCLFacetSums : public pointwise_rule {
    private:
      const_multiMap face2node ;
      const_param<real> timeStep ;
      const_store<vect3d> node_v,posOld,pos ;
      const_store<vect3d> faceCenter,oldFaceCenter ;
      store<vect3d> facetNormalSum ;
      store<real> facetVDotNormalSum ;
    public:

      // Define input and output.
      GCLFacetSums() {
        name_store("face2node",face2node) ;
        name_store("timeStep{n}",timeStep) ;
        name_store("node_v{n,it}",node_v) ;
        name_store("pos{n}",posOld) ;
        name_store("pos{n,it}",pos) ;
        name_store("facecenter{n,it}",faceCenter) ;
        name_store("facecenter{n}",oldFaceCenter) ;
        name_store("facetNormalSum{n,it}",facetNormalSum) ;
        name_store("facetVDotNormalSum{n,it}",facetVDotNormalSum) ;
        input("timeStep{n}") ;
        input("face2node->(node_v{n,it},pos{n},pos{n,it})") ;
        input("facecenter{n,it},facecenter{n}") ;
        output("facetNormalSum{n,it}") ;
        output("facetVDotNormalSum{n,it}") ;
      }

      // Calculate the mass flux for a single face.
      void calculate(Entity face) {

        // Save face centers.
        const vect3d center=faceCenter[face],oldCenter=oldFaceCenter[face] ;

        // Compute face center velocity. Fixed the BUG where the face center
        // velocity was computed incorrectly by averaging the face nodal
        // velocities. Must use the actual old and new face centers since the
        // coordinates of the face center are obtained elsewhere, NOT by using
        // the average of the face node coordinates.
        const vect3d face_v=(1.0/(*timeStep))*(center-oldCenter) ;

        // Face must be broken into its component facets to be GCL compliant.
        // Add the physical velocity and grid velocity components to the mass
        // flux on a facet-by-facet basis. Note that we need to use the exact
        // velocity of each facet center, but can still just use the old face
        // fluid velocity averaged from the centers.
        int n0=*(face2node.begin(face)),n1=*(face2node.end(face)-1) ;
        vect3d facetNormal=(cross(posOld[n1]-oldCenter,posOld[n0]-oldCenter)+
          cross(pos[n1]-center,pos[n0]-center)+0.5*(cross(posOld[n1]-oldCenter,
          pos[n0]-center)+cross(pos[n1]-center,posOld[n0]-oldCenter)))/6.0 ;
        vect3d facetVelocity=(node_v[n0]+node_v[n1]+face_v)/3.0 ;
        facetNormalSum[face]=facetNormal ;
        facetVDotNormalSum[face]=dot(facetVelocity,facetNormal) ;
        for(const int *ni=face2node.begin(face)+1;ni!=(face2node.end(face));
        ni++) {
          n0=*(ni),n1=*(ni-1) ;
          facetNormal=(cross(posOld[n1]-oldCenter,posOld[n0]-oldCenter)+
            cross(pos[n1]-center,pos[n0]-center)+0.5*(cross(posOld[n1]-
            oldCenter,pos[n0]-center)+cross(pos[n1]-center,posOld[n0]-
            oldCenter)))/6.0 ;
          facetVelocity=(node_v[n0]+node_v[n1]+face_v)/3.0 ;
          facetNormalSum[face]+=facetNormal ;
          facetVDotNormalSum[face]+=dot(facetVelocity,facetNormal) ;
        }
      }

      // Calculate the mass flux for all faces in a sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;

  register_rule<GCLFacetSums> registerGCLFacetSums ;

  // Priority rule for grid mass flux over interior faces.
  class PriorityGridMassFluxInteriorBDF : public pointwise_rule {
    private:
      const_store<real> facetVDotNormalSum ;
      const_store<real> faceDensity ;
      store<real> gridMassFlux ;
    public:

      // Define input and output.
      PriorityGridMassFluxInteriorBDF() {
        name_store("facetVDotNormalSum{n,it}",facetVDotNormalSum) ;
        name_store("faceDensity{n,it}",faceDensity) ;
        name_store("priority::gridMassFlux{n,it}",gridMassFlux) ;
        input("facetVDotNormalSum{n,it},faceDensity{n,it}") ;
        output("priority::gridMassFlux{n,it}") ;
        constraint("BDFIntegrator{n,it},internalFaces{n,it}") ;
      }

      // Calculate the mass flux for a single face.
      void calculate(Entity face) {
        gridMassFlux[face]=faceDensity[face]*facetVDotNormalSum[face] ;
      }

      // Calculate the mass flux for all faces in a sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;

  register_rule<PriorityGridMassFluxInteriorBDF>
    registerPriorityGridMassFluxInteriorBDF ;

  // Priority rule for grid mass flux over inlet and outlet faces.
  class PriorityGridMassFluxBoundaryBDF : public pointwise_rule {
    private:
      const_store<real> facetVDotNormalSum ;
      const_store<real> rho_f ;
      store<real> gridMassFlux ;
    public:

      // Define input and output.
      PriorityGridMassFluxBoundaryBDF() {
        name_store("facetVDotNormalSum{n,it}",facetVDotNormalSum) ;
        name_store("rho_f{n,it}",rho_f) ;
        name_store("priority::gridMassFlux{n,it}",gridMassFlux) ;
        input("facetVDotNormalSum{n,it},rho_f{n,it}") ;
        output("priority::gridMassFlux{n,it}") ;
        constraint("BDFIntegrator{n,it},inletOutlet_BC{n,it}") ;
      }

      // Calculate the mass flux for a single face.
      void calculate(Entity face) {
        gridMassFlux[face]=rho_f[face]*facetVDotNormalSum[face] ;
      }

      // Calculate the mass flux for all faces in a sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;

  register_rule<PriorityGridMassFluxBoundaryBDF>
    registerPriorityGridMassFluxBoundaryBDF ;

  // Priority rule for grid mass flux over interior faces.
  class PriorityGridMassFluxInteriorBDF2 : public pointwise_rule {
    private:
      const_param<int> n ;
      const_store<real> oldFacetVDotNormalSum ;
      const_store<real> facetVDotNormalSum ;
      const_store<real> faceDensity ;
      store<real> gridMassFlux ;
    public:

      // Define input and output.
      PriorityGridMassFluxInteriorBDF2() {
        name_store("$n{n}",n) ;
        name_store("oldFacetVDotNormalSum{n}",oldFacetVDotNormalSum) ;
        name_store("facetVDotNormalSum{n,it}",facetVDotNormalSum) ;
        name_store("faceDensity{n,it}",faceDensity) ;
        name_store("priority::gridMassFlux{n,it}",gridMassFlux) ;
        input("oldFacetVDotNormalSum{n},facetVDotNormalSum{n,it}") ;
        input("faceDensity{n,it}") ;
        output("priority::gridMassFlux{n,it}") ;
        constraint("BDF2Integrator{n,it},internalFaces{n,it}") ;
      }

      // Calculate the mass flux for a single face.
      void calculate(Entity face) {
        if((*n)!=0){
          gridMassFlux[face]=faceDensity[face]*(1.5*facetVDotNormalSum[face]-
            0.5*oldFacetVDotNormalSum[face]) ;
        }else{
          gridMassFlux[face]=faceDensity[face]*facetVDotNormalSum[face] ;
        }
      }

      // Calculate the mass flux for all faces in a sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;

  register_rule<PriorityGridMassFluxInteriorBDF2>
    registerPriorityGridMassFluxInteriorBDF2 ;

  // Priority rule for grid mass flux over inlet and outlet faces.
  class PriorityGridMassFluxBoundaryBDF2 : public pointwise_rule {
    private:
      const_param<int> n ;
      const_store<real> oldFacetVDotNormalSum ;
      const_store<real> facetVDotNormalSum ;
      const_store<real> rho_f ;
      store<real> gridMassFlux ;
    public:

      // Define input and output.
      PriorityGridMassFluxBoundaryBDF2() {
        name_store("$n{n}",n) ;
        name_store("oldFacetVDotNormalSum{n}",oldFacetVDotNormalSum) ;
        name_store("facetVDotNormalSum{n,it}",facetVDotNormalSum) ;
        name_store("rho_f{n,it}",rho_f) ;
        name_store("priority::gridMassFlux{n,it}",gridMassFlux) ;
        input("oldFacetVDotNormalSum{n},facetVDotNormalSum{n,it},rho_f{n,it}") ;
        output("priority::gridMassFlux{n,it}") ;
        constraint("BDF2Integrator{n,it},inletOutlet_BC{n,it}") ;
      }

      // Calculate the mass flux for a single face.
      void calculate(Entity face) {
        if((*n)!=0){
          gridMassFlux[face]=rho_f[face]*(1.5*facetVDotNormalSum[face]-
            0.5*oldFacetVDotNormalSum[face]) ;
        }else{
          gridMassFlux[face]=rho_f[face]*facetVDotNormalSum[face] ;
        }
      }

      // Calculate the mass flux for all faces in a sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;

  register_rule<PriorityGridMassFluxBoundaryBDF2>
    registerPriorityGridMassFluxBoundaryBDF2 ;

  // GCL mass flux for stage 0. We are not going to allow axisymmetric flow
  // with moving meshes, so we do not incorporate face radius here. No change
  // required for turbomachinery.
  class GCLStageZeroMassFluxInteriorFOUSOU : public pointwise_rule {
    private:
      const_param<int> it ;
      const_param<real> vRelaxationFactor ;
      const_Map cl,cr ;
      const_store<vect3d> v ;
      const_store<vect3d> facetNormalSum ;
      const_store<real> faceDensity ;
      const_store<real> massFlux ;
      store<real> stageZeroMassFlux ;
    public:

      // Define input and output.
      GCLStageZeroMassFluxInteriorFOUSOU() {
        name_store("$it{n,it}",it) ;
        name_store("vRelaxationFactor{n,it}",vRelaxationFactor) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("v{n,it}",v) ;
        name_store("facetNormalSum{n,it}",facetNormalSum) ;
        name_store("faceDensity{n,it}",faceDensity) ;
        name_store("massFlux{n,it}",massFlux) ;
        name_store("priority::stageZeroMassFlux{n,it}",stageZeroMassFlux) ;
        input("$it{n,it},vRelaxationFactor{n,it},(cl,cr)->v{n,it}") ;
        input("facetNormalSum{n,it},faceDensity{n,it},massFlux{n,it}") ;
        output("priority::stageZeroMassFlux{n,it}") ;
        constraint("BDFIntegrator{n,it},internalFaces{n,it}") ;
      }

      // Calculate the mass flux for a single face.
      void calculate(Entity face) {
        if(*it==0){
          stageZeroMassFlux[face]=0.0 ;
        }else{
          stageZeroMassFlux[face]=(1.0-(*vRelaxationFactor))*(massFlux[face]-
            0.5*faceDensity[face]*dot(v[cl[face]]+v[cr[face]],
            facetNormalSum[face])) ;
        }
      }

      // Calculate the mass flux for all faces in a sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;

  register_rule<GCLStageZeroMassFluxInteriorFOUSOU>
    registerGCLStageZeroMassFluxInteriorFOUSOU ;

  // GCL version of the above rule. No change required for turbomachinery.
  class GCLStageOneMassFluxInteriorFOUSOU : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<vect3d> facetNormalSum ;
      const_store<vect3d> vStar ;
      const_store<real> p ;
      const_store<vect3d> pGradient ;
      const_store<vect3d> cellCenter ;
      const_store<real> vol ;
      const_store<real> faceDensity ;
      const_store<real> stageZeroMassFlux ;
      const_store<real> pPrimeCoefficient ;
      store<real> stageOneMassFlux ;
    public:

      // Define input and output.
      GCLStageOneMassFluxInteriorFOUSOU() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("facetNormalSum{n,it}",facetNormalSum) ;
        name_store("vStar{n,it}",vStar) ;
        name_store("p{n,it}",p) ;
        name_store("grads(p){n,it}",pGradient) ;
        name_store("cellcenter{n,it}",cellCenter) ;
        name_store("vol{n,it}",vol) ;
        name_store("faceDensity{n,it}",faceDensity) ;
        name_store("stageZeroMassFlux{n,it}",stageZeroMassFlux) ;
        name_store("pPrimeCoefficient{n,it}",pPrimeCoefficient) ;
        name_store("priority::stageOneMassFlux{n,it}",stageOneMassFlux) ;
        input("facetNormalSum{n,it}") ;
        input("(cl,cr)->(vStar{n,it},p{n,it},grads(p){n,it})") ;
        input("(cl,cr)->(cellcenter{n,it},vol{n,it})") ;
        input("faceDensity{n,it},stageZeroMassFlux{n,it}") ;
        input("pPrimeCoefficient{n,it}") ;
        output("priority::stageOneMassFlux{n,it}") ;
        constraint("BDFIntegrator,internalFaces,fouOrSouInviscidFlux") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        vect3d averagePressureGradient=(pGradient[cl[face]]*vol[cr[face]]+
          pGradient[cr[face]]*vol[cl[face]])/(vol[cl[face]]+vol[cr[face]]) ;
        stageOneMassFlux[face]=stageZeroMassFlux[face]-pPrimeCoefficient[face]*
          (p[cr[face]]-p[cl[face]]-dot(averagePressureGradient,
          cellCenter[cr[face]]-cellCenter[cl[face]]))+0.5*faceDensity[face]*
          dot(vStar[cl[face]]+vStar[cr[face]],facetNormalSum[face]) ;
      }

      // Calculate mass flux for all faces in a sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;

  register_rule<GCLStageOneMassFluxInteriorFOUSOU>
    registerGCLStageOneMassFluxInteriorFOUSOU ;

  // GCL mass flux for boundary faces.
  class GCLStageOneMassFluxBoundaryFOUSOU : public pointwise_rule {
    private:
      const_multiMap face2node ;
      const_store<vect3d> facetNormalSum ;
      const_store<real> rho_f ;
      const_store<vect3d> v_f ;
      store<real> stageOneMassFlux ;
                                                                                
    public:

      // Define input and output.
      GCLStageOneMassFluxBoundaryFOUSOU() {
        name_store("facetNormalSum{n,it}",facetNormalSum) ;
        name_store("rho_f{n,it}",rho_f) ;
        name_store("v_f{n,it}",v_f) ;
        name_store("priority::stageOneMassFlux{n,it}",stageOneMassFlux) ;
        input("facetNormalSum{n,it},rho_f{n,it},v_f{n,it}") ;
        output("priority::stageOneMassFlux{n,it}") ;
        constraint("BDFIntegrator,inletOutlet_BC") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        stageOneMassFlux[face]=rho_f[face]*dot(v_f[face],facetNormalSum[face]) ;
      }

      // Calculate mass flux for all faces in a sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;

  register_rule<GCLStageOneMassFluxBoundaryFOUSOU>
    registerGCLStageOneMassFluxBoundaryFOUSOU ;

  // GCL mass flux for stage 0. We are not going to allow axisymmetric flow
  // with moving meshes, so we do not incorporate face radius here. No change
  // required for turbomachinery.
  class GCLStageZeroMassFluxInteriorFOUSOUBDF2 : public pointwise_rule {
    private:
      const_param<int> it ;
      const_param<int> n ;
      const_param<real> vRelaxationFactor ;
      const_Map cl,cr ;
      const_store<vect3d> v ;
      const_store<vect3d> oldFacetNormalSum,facetNormalSum ;
      const_store<real> faceDensity ;
      const_store<real> massFlux ;
      store<real> stageZeroMassFlux ;
    public:

      // Define input and output.
      GCLStageZeroMassFluxInteriorFOUSOUBDF2() {
        name_store("$it{n,it}",it) ;
        name_store("$n{n}",n) ;
        name_store("vRelaxationFactor{n,it}",vRelaxationFactor) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("v{n,it}",v) ;
        name_store("oldFacetNormalSum{n}",oldFacetNormalSum) ;
        name_store("facetNormalSum{n,it}",facetNormalSum) ;
        name_store("faceDensity{n,it}",faceDensity) ;
        name_store("massFlux{n,it}",massFlux) ;
        name_store("priority::stageZeroMassFlux{n,it}",stageZeroMassFlux) ;
        input("$it{n,it},$n{n},vRelaxationFactor{n,it},(cl,cr)->v{n,it}") ;
        input("oldFacetNormalSum{n},facetNormalSum{n,it}") ;
        input("faceDensity{n,it},massFlux{n,it}") ;
        output("priority::stageZeroMassFlux{n,it}") ;
        constraint("BDF2Integrator{n,it},internalFaces{n,it}") ;
      }

      // Calculate the mass flux for a single face. Degenerates to BDF on the
      // first timestep.
      void calculate(Entity face) {

        // We had this in the eqivalent BDF rule. Forgot why we are doing this.
        // Try to refresh you memory sometime and see if this is necessary.
        if(*it==0){ stageZeroMassFlux[face]=0.0 ; return ; }

        if((*n)!=0){
          stageZeroMassFlux[face]=(1.0-(*vRelaxationFactor))*(massFlux[face]-
            faceDensity[face]*0.5*dot(v[cl[face]]+v[cr[face]],
            1.5*facetNormalSum[face]-0.5*oldFacetNormalSum[face])) ;
        }else{
          stageZeroMassFlux[face]=(1.0-(*vRelaxationFactor))*(massFlux[face]-
            faceDensity[face]*0.5*dot(v[cl[face]]+v[cr[face]],
            facetNormalSum[face])) ;
        }
      }

      // Calculate the mass flux for all faces in a sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;

  register_rule<GCLStageZeroMassFluxInteriorFOUSOUBDF2>
    registerGCLStageZeroMassFluxInteriorFOUSOUBDF2 ;

  // GCL version of the above rule. No change required for turbomachinery.
  class GCLStageOneMassFluxInteriorFOUSOUBDF2 : public pointwise_rule {
    private:
      const_param<int> n ;
      const_Map cl,cr ;
      const_store<vect3d> oldFacetNormalSum,facetNormalSum ;
      const_store<vect3d> vStar ;
      const_store<real> p ;
      const_store<vect3d> pGradient ;
      const_store<vect3d> cellCenter ;
      const_store<real> vol ;
      const_store<real> faceDensity ;
      const_store<real> stageZeroMassFlux ;
      const_store<real> pPrimeCoefficient ;
      store<real> stageOneMassFlux ;
    public:

      // Define input and output.
      GCLStageOneMassFluxInteriorFOUSOUBDF2() {
        name_store("$n{n}",n) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("oldFacetNormalSum{n}",oldFacetNormalSum) ;
        name_store("facetNormalSum{n,it}",facetNormalSum) ;
        name_store("vStar{n,it}",vStar) ;
        name_store("p{n,it}",p) ;
        name_store("grads(p){n,it}",pGradient) ;
        name_store("cellcenter{n,it}",cellCenter) ;
        name_store("vol{n,it}",vol) ;
        name_store("faceDensity{n,it}",faceDensity) ;
        name_store("stageZeroMassFlux{n,it}",stageZeroMassFlux) ;
        name_store("pPrimeCoefficient{n,it}",pPrimeCoefficient) ;
        name_store("priority::stageOneMassFlux{n,it}",stageOneMassFlux) ;
        input("$n{n},oldFacetNormalSum{n},facetNormalSum{n,it}") ;
        input("(cl,cr)->(vStar{n,it},p{n,it},grads(p){n,it})") ;
        input("(cl,cr)->(cellcenter{n,it},vol{n,it})") ;
        input("faceDensity{n,it},stageZeroMassFlux{n,it}") ;
        input("pPrimeCoefficient{n,it}") ;
        output("priority::stageOneMassFlux{n,it}") ;
        constraint("BDF2Integrator,internalFaces,fouOrSouInviscidFlux") ;
      }

      // Calculate mass flux for a single face. Degenerates to BDF on the
      // first timestep.
      void calculate(Entity face) {
        vect3d averagePressureGradient=(pGradient[cl[face]]*vol[cr[face]]+
          pGradient[cr[face]]*vol[cl[face]])/(vol[cl[face]]+vol[cr[face]]) ;
        if((*n)!=0){
          stageOneMassFlux[face]=stageZeroMassFlux[face]-
            pPrimeCoefficient[face]*(p[cr[face]]-p[cl[face]]-
            dot(averagePressureGradient,cellCenter[cr[face]]-
            cellCenter[cl[face]]))+0.5*faceDensity[face]*dot(vStar[cl[face]]+
            vStar[cr[face]],1.5*facetNormalSum[face]-
            0.5*oldFacetNormalSum[face]) ;
        }else{
          stageOneMassFlux[face]=stageZeroMassFlux[face]-
            pPrimeCoefficient[face]*(p[cr[face]]-p[cl[face]]-
            dot(averagePressureGradient,cellCenter[cr[face]]-
            cellCenter[cl[face]]))+0.5*faceDensity[face]*dot(vStar[cl[face]]+
            vStar[cr[face]],facetNormalSum[face]) ;
        }
      }

      // Calculate mass flux for all faces in a sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;

  register_rule<GCLStageOneMassFluxInteriorFOUSOUBDF2>
    registerGCLStageOneMassFluxInteriorFOUSOUBDF2 ;

  // GCL mass flux for boundary faces.
  class GCLStageOneMassFluxBoundaryFOUSOUBDF2 : public pointwise_rule {
    private:
      const_param<int> n ;
      const_store<vect3d> oldFacetNormalSum,facetNormalSum ;
      const_store<real> rho_f ;
      const_store<vect3d> v_f ;
      store<real> stageOneMassFlux ;
                                                                                
    public:

      // Define input and output.
      GCLStageOneMassFluxBoundaryFOUSOUBDF2() {
        name_store("$n{n}",n) ;
        name_store("oldFacetNormalSum{n}",oldFacetNormalSum) ;
        name_store("facetNormalSum{n,it}",facetNormalSum) ;
        name_store("rho_f{n,it}",rho_f) ;
        name_store("v_f{n,it}",v_f) ;
        name_store("priority::stageOneMassFlux{n,it}",stageOneMassFlux) ;
        input("$n{n},oldFacetNormalSum{n},facetNormalSum{n,it}") ;
        input("rho_f{n,it},v_f{n,it}") ;
        output("priority::stageOneMassFlux{n,it}") ;
        constraint("BDF2Integrator,inletOutlet_BC") ;
      }

      // Calculate mass flux for a single face. Degenerates to BDF on the
      // first timestep.
      void calculate(Entity face) {
        if((*n)!=0){
          stageOneMassFlux[face]=rho_f[face]*dot(v_f[face],
            1.5*facetNormalSum[face]-0.5*oldFacetNormalSum[face]) ;
        }else{
          stageOneMassFlux[face]=rho_f[face]*dot(v_f[face],
            facetNormalSum[face]) ;
        }
      }

      // Calculate mass flux for all faces in a sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;

  register_rule<GCLStageOneMassFluxBoundaryFOUSOUBDF2>
    registerGCLStageOneMassFluxBoundaryFOUSOUBDF2 ;
}
