//-----------------------------------------------------------------------------
// Description: This file contains rules required for the PISO algorithm.
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

  // Computes the new face center based on the updated nodal positions.
  class NewFaceCenter : public pointwise_rule {
    private:
      const_multiMap face2node ;
      const_store<vect3d> posStar ;
      store<vect3d> newFaceCenter ;
    public :

      // Define input and output.
     NewFaceCenter() {
        name_store("face2node",face2node) ;
        name_store("posStar",posStar) ;
        name_store("newFaceCenter",newFaceCenter) ;
        input("face2node->posStar") ;
        output("newFaceCenter") ;
     }

     // Compute center for a single face.
     void calculate(Entity fc) {
       vect3d nodesum(0.0,0.0,0.0) ; real lensum=0.0 ;
       warn(face2node.begin(fc)==face2node.end(fc) || face2node.begin(fc)+1==
         face2node.end(fc)) ;
       for(const int *id=face2node.begin(fc);id+1!=face2node.end(fc);++id) {
         vect3d pos0=posStar[*id],pos1=posStar[*(id+1)] ;
         vect3d edge_loc=0.5*(pos0+pos1),edge_vec=pos0-pos1 ;
         real len=sqrt(dot(edge_vec,edge_vec)) ;
         nodesum+=len*edge_loc ; lensum += len ;
       }
       const int *id=face2node.begin(fc),*idend=face2node.end(fc)-1 ;
       vect3d pos0=posStar[*id],pos1=posStar[*idend] ;
       vect3d edge_loc=0.5*(pos0+pos1),edge_vec=pos0-pos1 ;
       real len=sqrt(dot(edge_vec,edge_vec)) ;
       nodesum+=len*edge_loc ; lensum += len ;
       newFaceCenter[fc]=nodesum/lensum ;
     }

     // Loop through faces.
     virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<NewFaceCenter> registerNewFaceCenter ;

  // Compute new node velocities for solution-dependent grid motion.
  class NodeVelocityPISO : public pointwise_rule {
    private:
      const_param<real> timeStep ;
      const_store<vect3d> pos,posStar ;
      store<vect3d> node_v ;
    public:

      // Define input and output.
      NodeVelocityPISO() {
        name_store("timeStep{n}",timeStep) ;
        name_store("pos{n}",pos) ;
        name_store("posStar{n}",posStar) ;
        name_store("node_v{n}",node_v) ;
        input("timeStep{n},pos{n},posStar{n}") ;
        output("node_v{n}") ;
        constraint("PISO,pos") ;
      }

      // Position for single node.
      void calculate(Entity n) { node_v[n]=(posStar[n]-pos[n])/(*timeStep) ; }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<NodeVelocityPISO> registerNodeVelocityPISO ;

  // Rule to compute the facet normal sum and sum of the dot product of facet
  // velocity with the facet normal for each face.
  class GCLFacetSumsPISO : public pointwise_rule {
    private:
      const_multiMap face2node ;
      const_param<real> timeStep ;
      const_store<vect3d> node_v,pos,posStar ;
      const_store<vect3d> oldFaceCenter,newFaceCenter ;
      store<vect3d> facetNormalSum ;
      store<real> facetVDotNormalSum ;
    public:

      // Define input and output.
      GCLFacetSumsPISO() {
        name_store("face2node",face2node) ;
        name_store("timeStep{n}",timeStep) ;
        name_store("node_v{n}",node_v) ;
        name_store("pos{n}",pos) ;
        name_store("posStar{n}",posStar) ;
        name_store("facecenter{n}",oldFaceCenter) ;
        name_store("newFaceCenter{n}",newFaceCenter) ;
        name_store("facetNormalSum{n}",facetNormalSum) ;
        name_store("facetVDotNormalSum{n}",facetVDotNormalSum) ;
        input("timeStep{n}") ;
        input("face2node->(node_v{n},pos{n},posStar{n})") ;
        input("facecenter{n},newFaceCenter{n}") ;
        output("facetNormalSum{n},facetVDotNormalSum{n}") ;
        constraint("PISO{n},faces{n}") ;
      }

      // Calculate the mass flux for a single face.
      void calculate(Entity face) {

        // Save face centers.
        const vect3d center=newFaceCenter[face],oldCenter=oldFaceCenter[face] ;

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
        vect3d facetNormal=(cross(pos[n1]-oldCenter,pos[n0]-oldCenter)+
          cross(posStar[n1]-center,posStar[n0]-center)+0.5*(cross(pos[n1]-oldCenter,
          posStar[n0]-center)+cross(posStar[n1]-center,pos[n0]-oldCenter)))/6.0 ;
        vect3d facetVelocity=(node_v[n0]+node_v[n1]+face_v)/3.0 ;
        facetNormalSum[face]=facetNormal ;
        facetVDotNormalSum[face]=dot(facetVelocity,facetNormal) ;
        for(const int *ni=face2node.begin(face)+1;ni!=(face2node.end(face));
        ni++) {
          n0=*(ni),n1=*(ni-1) ;
          facetNormal=(cross(pos[n1]-oldCenter,pos[n0]-oldCenter)+
            cross(posStar[n1]-center,posStar[n0]-center)+0.5*(cross(pos[n1]-
            oldCenter,posStar[n0]-center)+cross(posStar[n1]-center,pos[n0]-
            oldCenter)))/6.0 ;
          facetVelocity=(node_v[n0]+node_v[n1]+face_v)/3.0 ;
          facetNormalSum[face]+=facetNormal ;
          facetVDotNormalSum[face]+=dot(facetVelocity,facetNormal) ;
        }
      }
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<GCLFacetSumsPISO> registerGCLFacetSumsPISO ;

  // Rule to compute the new mass flux on interior faces for the FOU and
  // SOU convection schemes. New form of momentum interpolation added by
  // ST to prevent decoupling observed with BDF2 on cylinder case.
  class GCLMassFluxStarInteriorFOUSOU : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<vect3d> facetNormalSum ;
      const_store<vect3d> vStar ;
      const_store<real> p ;
      const_store<vect3d> pGradient ;
      const_store<vect3d> cellCenter ;
      const_store<real> vol ;
      const_store<real> faceDensity ;
      const_store<Area> area ;
      const_store<real> diffusionProduct ;
      const_store<real> pPrimeCoefficient ;
      store<real> massFluxStar ;
    public:

      // Define input and output.
      GCLMassFluxStarInteriorFOUSOU() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("facetNormalSum",facetNormalSum) ;
        name_store("vStar",vStar) ;
        name_store("p",p) ;
        name_store("grads(p)",pGradient) ;
        name_store("cellcenter",cellCenter) ;
        name_store("vol",vol) ;
        name_store("faceDensity",faceDensity) ;
        name_store("area",area) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        name_store("GCL::massFluxStar",massFluxStar) ;
        name_store("diffusionProduct",diffusionProduct) ;
        input("facetNormalSum,(cl,cr)->(vStar,p,grads(p),cellcenter,vol)") ;
        input("faceDensity,diffusionProduct,pPrimeCoefficient,area") ;
        output("GCL::massFluxStar") ;
        constraint("PISO,BDFIntegrator") ;
        constraint("internalFaces,fouOrSouInviscidFlux") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        vect3d averagePressureGradient=(pGradient[cl[face]]*vol[cr[face]]+
          pGradient[cr[face]]*vol[cl[face]])/(vol[cl[face]]+vol[cr[face]]) ;
        massFluxStar[face]=0.5*faceDensity[face]*dot(vStar[cl[face]]+
          vStar[cr[face]],facetNormalSum[face])-pPrimeCoefficient[face]*
          (p[cr[face]]-p[cl[face]]-dot(averagePressureGradient,area[face].n)*
          area[face].sada/diffusionProduct[face]) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<GCLMassFluxStarInteriorFOUSOU>
    registerGCLMassFluxStarInteriorFOUSOU ;

  // Rule to compute the new mass flux in the corrector stage on interior
  // faces for the FOU and SOU convection schemes.
  class GCLMassFluxStarHatInteriorFOUSOU : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<vect3d> facetNormalSum ;
      const_store<vect3d> vHat;
      const_store<real> vMainCoefficient_it,vMainCoefficient ;
      const_store<real> vol ;
      const_store<real> faceDensity,massFlux ;
      const_store<real> faceRadius ;
      store<real> massFluxStarHat ;
    public:

      // Define input and output.
      GCLMassFluxStarHatInteriorFOUSOU() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("facetNormalSum",facetNormalSum) ;
        name_store("vHat",vHat) ;
        name_store("vMainCoefficient_it",vMainCoefficient_it) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        name_store("vol",vol) ;
        name_store("faceDensity",faceDensity) ;
        name_store("massFlux",massFlux) ;
        name_store("faceRadius",faceRadius) ;
        name_store("GCL::massFluxStarHat",massFluxStarHat) ;
        input("(cl,cr)->(vHat,vMainCoefficient_it,vMainCoefficient,vol)") ;
        input("faceDensity,massFlux,faceRadius,facetNormalSum") ;
        output("GCL::massFluxStarHat") ;
        constraint("PISO,BDFIntegrator") ;
        constraint("internalFaces,fouOrSouInviscidFlux") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        real vL=vol[cl[face]],vR=vol[cr[face]] ;
        real vMainCoefficientRatio=(vR*(vMainCoefficient_it[cl[face]]/
          vMainCoefficient[cl[face]])+vL*(vMainCoefficient_it[cr[face]]/
          vMainCoefficient[cr[face]]))/(vL+vR) ;
        vect3d vHatAverage=(vR*vHat[cl[face]]+vL*vHat[cr[face]])/(vL+vR) ;
        massFluxStarHat[face]=vMainCoefficientRatio*massFlux[face]+
          faceDensity[face]*dot(vHatAverage,facetNormalSum[face]) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<GCLMassFluxStarHatInteriorFOUSOU>
    registerGCLMassFluxStarHatInteriorFOUSOU ;

  // Mass flux for boundaries where velocity is specified.
  class GCLBoundaryMassFluxSpecifiedVelocityBDF : public pointwise_rule {
    private:
      const_store<real> rho_f ;
      const_store<vect3d> v_f ;
      const_store<vect3d> facetNormalSum ;
      const_store<real> faceRadius ;
      store<real> massFlux ;
    public:

      // Define input and output.
      GCLBoundaryMassFluxSpecifiedVelocityBDF() {
        name_store("rho_f",rho_f) ;
        name_store("v_f",v_f) ;
        name_store("facetNormalSum",facetNormalSum) ;
        name_store("faceRadius",faceRadius) ;
        name_store("GCL::massFlux",massFlux) ;
        input("rho_f,v_f,facetNormalSum,faceRadius") ;
        output("GCL::massFlux") ;
        constraint("PISO,BDFIntegrator,ref->v_BCoption") ;
      }

      // Calculate mass flux for a single face.
      void calculate(Entity face) {
        massFlux[face]=rho_f[face]*dot(v_f[face],facetNormalSum[face]) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<GCLBoundaryMassFluxSpecifiedVelocityBDF>
    registerGCLBoundaryMassFluxSpecifiedVelocityBDF ;
}
