//-----------------------------------------------------------------------------
// Description: This file contains rules for assembling and solving the linear
//   elasticity equations for the movement of interior nodes given the
//   boundary node displacements.
//-----------------------------------------------------------------------------
                                                                                
// Standard library includes.
#include <map>
using std::map ;
#include <vector>
using std::vector ;
                                                                                
// Loci includes.
#include <Loci.h>
                                                                                
// StreamUns includes.
#include "const.h"
#include "FSI_move.h"
#include "residual.h"
#include "sciTypes.h"
                                                                                
namespace streamUns {

  // Rule for creating constraint for all the nodes.
  class NodeConstraint : public pointwise_rule {
    private:
      const_multiMap face2node ;
      store<bool> nodes ;
    public:

      // Define input and output.
      NodeConstraint() {
        name_store("face2node",face2node) ;
        name_store("nodes",nodes) ;
        input("face2node") ;
        output("face2node->nodes") ;
        constraint("faces") ;
      }

      // Do nothing for a face since this rule is only used to establish
      // a constraint.
      void calculate(Entity face) {}

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<NodeConstraint> registerNodeConstraint ;

  // Rule for creating a boundary node constraint.
  class BoundaryNodeConstraint : public pointwise_rule {
    private:
      const_multiMap face2node ;
      store<bool> boundaryNodes ;
    public:

      // Define input and output.
      BoundaryNodeConstraint() {
        name_store("face2node",face2node) ;
        name_store("boundaryNodes",boundaryNodes) ;
        input("face2node") ;
        output("face2node->boundaryNodes") ;
        constraint("boundaryFaces") ;
      }

      // Do nothing for a face since this rule is only used to establish
      // a constraint.
      void calculate(Entity face) {}

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryNodeConstraint> registerBoundaryNodeConstraint ;

  // Rule for printing out node positions. Only used as a check.
  class NodePosition : public pointwise_rule {
    private:
      const_store<vect3d> pos ;
      store<vect3d> nodePosition ;
      unsigned int numNode ;
    public:

      // Define input and output.
      NodePosition() {
        name_store("pos",pos) ;
        name_store("nodePosition",nodePosition) ;
        output("nodePosition") ;
        constraint("nodes") ;
      }

      // Do nothing for a node.
      void calculate(Entity node) { ++numNode ; }

      // Loop over nodes.
      virtual void compute(const sequence &seq) {
        numNode=0 ; do_loop(seq,this) ;
        cout << "numNode: " << numNode << endl ;
      }
  } ;

  register_rule<NodePosition> registerNodePosition ;

  // Rule for printing out boundary node positions. Only used as a check.
  class BoundaryNodePosition : public pointwise_rule {
    private:
      const_store<vect3d> pos ;
      store<vect3d> boundaryNodePosition ;
      unsigned int numNode ;
    public:

      // Define input and output.
      BoundaryNodePosition() {
        name_store("pos",pos) ;
        name_store("boundaryNodePosition",boundaryNodePosition) ;
        output("boundaryNodePosition") ;
        constraint("boundaryNodes") ;
      }

      // Do nothing for a node.
      void calculate(Entity node) { ++numNode ; }

      // Loop over nodes.
      virtual void compute(const sequence &seq) {
        numNode=0 ; do_loop(seq,this) ;
        cout << "numBoundaryNode: " << numNode << endl ;
      }
  } ;

  register_rule<BoundaryNodePosition> registerBoundaryNodePosition ;

//-----------------------------------------------------------------------------
// Rules for assembling the linear elasticity equation for node displacement.

  // Factor commonly used.
  class GridMoverFactor0 : public pointwise_rule {
    private:
      const_store<real> chi ;
      const_store<real> diffusionVolume ;
      store<real> gridMoverFactor0 ;
    public:

      // Define input and output.
      GridMoverFactor0() {
        name_store("chi",chi) ;
        name_store("diffusionVolume",diffusionVolume) ;
        name_store("gridMoverFactor0",gridMoverFactor0) ;
        input("chi,diffusionVolume") ;
        output("gridMoverFactor0") ;
        constraint("edges") ;
      }

      // Set the edge value.
      void calculate(Entity edge) {
	//cout << " diffusion volume = " << diffusionVolume[edge] ;
        gridMoverFactor0[edge]=pow(diffusionVolume[edge],-chi[edge]) ;
      }

      // Loop over edges.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GridMoverFactor0> registerGridMoverFactor0 ;

  // Another combination commonly used. Checked thoroughly on 11/19/2008.
  class ScaledDualArea : public pointwise_rule {
    private:
      const_param<real> nu ;
      const_store<vect3d> dualArea ;
      const_store<real> gridMoverFactor0 ;
      store<vect3d> scaledDualArea ;
    public:

      // Define input and output.
      ScaledDualArea() {
        name_store("gridMoverNu",nu) ;
        name_store("dualArea",dualArea) ;
        name_store("gridMoverFactor0",gridMoverFactor0) ;
        name_store("scaledDualArea",scaledDualArea) ;
        input("gridMoverNu,dualArea,gridMoverFactor0") ;
        output("scaledDualArea") ;
        constraint("edges") ;
      }

      // Set the edge value.
      void calculate(Entity edge) {
        scaledDualArea[edge]=(gridMoverFactor0[edge]/(1.0-2.0*(*nu)))*
          dualArea[edge] ;
      }

      // Loop over edges.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ScaledDualArea> registerScaledDualArea ;

  // Facet area sums for internal faces.
  class InternalFacetAreaSum : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_multiMap face2edge ;
      const_MapVec<2> edge2node ;
      const_store<vect3d> pos ;
      const_store<vect3d> faceCenter ;
      const_store<vect3d> cellCenter ;
      const_store<real> diffusionVolume ;
      store<vector<vect3d> > dASumLeft,dASumRight ;
    public:

      // Define input and output.
      InternalFacetAreaSum() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("face2edge",face2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("pos",pos) ;
        name_store("facecenter",faceCenter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("diffusionVolume",diffusionVolume) ;
        name_store("dASumLeft",dASumLeft) ;
        name_store("dASumRight",dASumRight) ;
        input("face2edge->edge2node->pos,(cl,cr)->cellcenter,facecenter") ;
        input("face2edge->diffusionVolume") ;
        output("dASumLeft,dASumRight") ;
        constraint("internalFaces") ;
      }

      // Contributions to nodes for a given face.
      void calculate(Entity face) {
        dASumLeft[face].clear() ; dASumRight[face].clear() ;
        unsigned int numEdge=face2edge.num_elems(face) ;
        for(unsigned int i=0;i<numEdge;++i){
                                                                                
          // Get node positions and edge center.
          unsigned int edgeNum=face2edge[face][i] ;
          unsigned int n0=edge2node[edgeNum][0],n1=edge2node[edgeNum][1] ;
          vect3d pos0=pos[n0],pos1=pos[n1] ;
          vect3d edgeCenter=0.5*(pos0+pos1),edgeVector=pos1-pos0 ;
                                                                                
          // Note that we are ensuring here that both the front and backward
          // area vectors are outward to the diffusion sub-control volume.
          vect3d dAFrontLeft=cross(cellCenter[cl[face]]-pos1,faceCenter[face]-
            pos1),dABackLeft=cross(cellCenter[cl[face]]-pos0,faceCenter[face]-
            pos0),dAFrontRight=cross(cellCenter[cr[face]]-pos1,faceCenter[face]-
            pos1),dABackRight=cross(cellCenter[cr[face]]-pos0,faceCenter[face]-
            pos0) ;
          if(dot(dAFrontLeft,edgeVector)<0.0) dAFrontLeft=-1.0*dAFrontLeft ;
          if(dot(dABackLeft,edgeVector)>0.0) dABackLeft=-1.0*dABackLeft ;
          if(dot(dAFrontRight,edgeVector)<0.0) dAFrontRight=-1.0*dAFrontRight ;
          if(dot(dABackRight,edgeVector)>0.0) dABackRight=-1.0*dABackRight ;

          // Add the front the the back..
          dASumLeft[face].push_back((dAFrontLeft+dABackLeft)/(6.0*
            diffusionVolume[edgeNum])) ;
          dASumRight[face].push_back((dAFrontRight+dABackRight)/(6.0*
            diffusionVolume[edgeNum])) ;
        }
      }

      // Loop over edges.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InternalFacetAreaSum> registerInternalFacetAreaSum ;

  // Facet area sums for boundary faces.
  class BoundaryFacetAreaSum : public pointwise_rule {
    private:
      const_Map ci ;
      const_multiMap face2edge ;
      const_MapVec<2> edge2node ;
      const_store<vect3d> pos ;
      const_store<vect3d> faceCenter ;
      const_store<vect3d> cellCenter ;
      const_store<real> diffusionVolume ;
      store<vector<vect3d> > dASum,dABoundary ;
    public:

      // Define input and output.
      BoundaryFacetAreaSum() {
        name_store("ci",ci) ;
        name_store("face2edge",face2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("pos",pos) ;
        name_store("facecenter",faceCenter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("diffusionVolume",diffusionVolume) ;
        name_store("dASum",dASum) ;
        name_store("dABoundary",dABoundary) ;
        input("face2edge->edge2node->pos,ci->cellcenter,facecenter") ;
        input("face2edge->diffusionVolume") ;
        output("dASum,dABoundary") ;
        constraint("boundaryFaces") ;
      }

      // Contributions to nodes for a given face.
      void calculate(Entity face) {
        dASum[face].clear() ; dABoundary[face].clear() ;
        unsigned int numEdge=face2edge.num_elems(face) ;
        for(unsigned int i=0;i<numEdge;++i){
                                                                                
          // Get node positions and edge center.
          unsigned int edgeNum=face2edge[face][i] ;
          unsigned int n0=edge2node[edgeNum][0],n1=edge2node[edgeNum][1] ;
          vect3d pos0=pos[n0],pos1=pos[n1] ;
          vect3d edgeCenter=0.5*(pos0+pos1),edgeVector=pos1-pos0 ;

          // Note that we are ensuring here that both the front and backward
          // area vectors are outward to the diffusion sub-control volume.
          vect3d dAFront=cross(cellCenter[ci[face]]-pos1,faceCenter[face]-
            pos1),dABack=cross(cellCenter[ci[face]]-pos0,faceCenter[face]-
            pos0) ;
          if(dot(dAFront,edgeVector)<0.0) dAFront=-1.0*dAFront ;
          if(dot(dABack,edgeVector)>0.0) dABack=-1.0*dABack ;

          // Add the front the the back..
          dASum[face].push_back((dAFront+dABack)/
            (6.0*diffusionVolume[edgeNum])) ;

          // Compute the boundary area vector. Again, ensure that it is outward
          // to the diffusion sub-control volume.
          vect3d temp=cross(pos1-faceCenter[face],pos0-faceCenter[face]) ;
          if(dot(temp,faceCenter[face]-cellCenter[ci[face]])<0.0)
            temp=-1.0*temp ;
          dABoundary[face].push_back(temp/(6.0*diffusionVolume[edgeNum])) ;
        }
      }

      // Loop over edges.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryFacetAreaSum> registerBoundaryFacetAreaSum ;

  // Rule to initialize the main coefficient.
  class InitializeDisplacementMainCoefficient : public unit_rule {
    private:
      store<real> sMainCoefficient ;
    public:
                                                                                
      // Define input and output.
      InitializeDisplacementMainCoefficient() {
        name_store("sMainCoefficient{n}",sMainCoefficient) ;
        output("sMainCoefficient{n}") ;
        constraint("nodes{n}") ;
      }
                                                                                
      // Set the main coefficient to zero for a single node.
      void calculate(Entity node) { sMainCoefficient[node]=0.0 ; }
                                                                                
      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<InitializeDisplacementMainCoefficient>
    registerInitializeDisplacementMainCoefficient ;

  // Rule to add the implicit contribution from the diffusion term. While this
  // contribution will be distributed to the main coefficient of all nodes,
  // ultimately only the contributions to interior nodes are used. Checked
  // thoroughly on 11/19/2008.
  class DiffusionToDisplacementMainCoefficient : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_MapVec<2> edge2node ;
      const_store<real> gridMoverFactor0 ;
      const_store<vect3d> dualArea ;
      const_store<real> diffusionVolume ;
      store<real> sMainCoefficient ;
    public:

      // Define input and output.
      DiffusionToDisplacementMainCoefficient() {
        name_store("edge2node{n}",edge2node) ;
        name_store("gridMoverFactor0{n}",gridMoverFactor0) ;
        name_store("dualArea{n}",dualArea) ;
        name_store("diffusionVolume{n}",diffusionVolume) ;
        name_store("sMainCoefficient{n}",sMainCoefficient) ;
        input("gridMoverFactor0{n},dualArea{n},diffusionVolume{n}") ;
        output("edge2node{n}->sMainCoefficient{n}") ;
        constraint("edges{n}") ;
      }

      // Increment the main coefficient for the nodes attached to the edge.
      void calculate(Entity edge) {
if(diffusionVolume[edge]<0.0)
cout << "edge,dv: " << edge << " " << diffusionVolume[edge] << endl ;
        real temp=dot(dualArea[edge],dualArea[edge])*gridMoverFactor0[edge]/
          (3.0*diffusionVolume[edge]) ;
        sMainCoefficient[edge2node[edge][0]]+=temp ;
        sMainCoefficient[edge2node[edge][1]]+=temp ;
      }

      // Loop over edges.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusionToDisplacementMainCoefficient>
    registerDiffusionToDisplacementMainCoefficient ;

  // Rule to initialize the source term.
  class InitializeDisplacementSourceTerm : public unit_rule {
    private:
      store<vect3d> sSourceTerm ;
    public:
                                                                                
      // Define input and output.
      InitializeDisplacementSourceTerm() {
        name_store("sSourceTerm{n,itg}",sSourceTerm) ;
        output("sSourceTerm{n,itg}") ;
        constraint("nodes{n,itg}") ;
      }
                                                                                
      // Set the source term to zero for a single node.
      void calculate(Entity node) {
        sSourceTerm[node]=vector3d<real>(0.0,0.0,0.0) ;
      }
                                                                                
      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<InitializeDisplacementSourceTerm>
    registerInitializeDisplacementSourceTerm ;

  // Apply rule for adding the main dilatation edge contribution to the
  // source term. Checked thoroughly on 11/19/2008.
  class DilatationToDisplacementSourceTerm : public apply_rule
  <store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_param<real> nu ;
      const_MapVec<2> edge2node ;
      const_store<real> gridMoverFactor0 ;
      const_store<vect3d> dualArea,scaledDualArea ;
      const_store<real> diffusionVolume ;
      const_store<vect3d> node_s ;
      store<vect3d> sSourceTerm ;
    public:

      // Define input and output.
      DilatationToDisplacementSourceTerm() {
        name_store("gridMoverNu{n,itg}",nu) ;
        name_store("edge2node{n,itg}",edge2node) ;
        name_store("gridMoverFactor0{n}",gridMoverFactor0) ;
        name_store("dualArea{n}",dualArea) ;
        name_store("scaledDualArea{n}",scaledDualArea) ;
        name_store("diffusionVolume{n}",diffusionVolume) ;
        name_store("node_sStar{n,itg}",node_s) ;
        name_store("sSourceTerm{n,itg}",sSourceTerm) ;
        input("gridMoverNu{n,itg},gridMoverFactor0{n},dualArea{n}") ;
        input("diffusionVolume{n},edge2node{n,itg}->node_sStar{n,itg}") ;
        output("edge2node{n,itg}->sSourceTerm{n,itg}") ;
        constraint("edges{n,itg}") ;
      }

      // Increment the main coefficient for the nodes attached to the edge.
      void calculate(Entity edge) {
        Entity n0=edge2node[edge][0],n1=edge2node[edge][1] ;
        vect3d temp=(dot(node_s[n1]-node_s[n0],dualArea[edge])/(3.0*
          diffusionVolume[edge]))*scaledDualArea[edge] ;
        sSourceTerm[n0]+=temp ; sSourceTerm[n1]-=temp ;
      }

      // Loop over edges.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DilatationToDisplacementSourceTerm>
    registerDilatationToDisplacementSourceTerm ;

  // Add the interior face cross-diffusion and cross-dilatation terms to the
  // nodal source terms. Checked thoroughly on 11/19/2008.
  class CrossDiffusionDilatationApplyInterior : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
      const_multiMap face2edge ;
      const_MapVec<2> edge2node ;
      const_param<real> nu ;
      const_store<real> gridMoverFactor0 ;
      const_store<vect3d> dualArea,scaledDualArea ;
      const_store<real> diffusionVolume ;
      const_store<vect3d> face_s ;
      const_store<vect3d> cell_s ;
      const_store<vector<vect3d> > dASumLeft,dASumRight ;
      store<vect3d> sSourceTerm ;
    public:

      // Define input and output.
      CrossDiffusionDilatationApplyInterior() {
        name_store("cl{n,itg}",cl) ;
        name_store("cr{n,itg}",cr) ;
        name_store("face2edge{n,itg}",face2edge) ;
        name_store("edge2node{n,itg}",edge2node) ;
        name_store("gridMoverNu{n,itg}",nu) ;
        name_store("gridMoverFactor0{n}",gridMoverFactor0) ;
        name_store("dualArea{n}",dualArea) ;
        name_store("scaledDualArea{n}",scaledDualArea) ;
        name_store("diffusionVolume{n}",diffusionVolume) ;
        name_store("face(node_sStar){n,itg}",face_s) ;
        name_store("cell(node_sStar){n,itg}",cell_s) ;
        name_store("dASumLeft{n}",dASumLeft) ;
        name_store("dASumRight{n}",dASumRight) ;
        name_store("sSourceTerm{n,itg}",sSourceTerm) ;
        input("gridMoverNu{n,itg},face2edge{n,itg}->gridMoverFactor0{n}") ;
        input("face2edge{n,itg}->(dualArea{n},diffusionVolume{n})") ;
        input("face2edge{n,itg}->scaledDualArea{n}") ;
        input("face(node_sStar){n,itg}") ;
        input("(cl{n,itg},cr{n,itg})->cell(node_sStar){n,itg}") ;
        input("dASumLeft{n},dASumRight{n}") ;
        output("face2edge{n,itg}->edge2node{n,itg}->sSourceTerm{n,itg}") ;
        constraint("internalFaces{n,itg}") ;
      }

      // Contributions to nodes for a given face.
      void calculate(Entity face) {
        unsigned int numEdge=face2edge.num_elems(face) ;
        for(unsigned int i=0;i<numEdge;++i){

          // Get node positions and edge center.
          unsigned int edgeNum=face2edge[face][i] ;
          unsigned int n0=edge2node[edgeNum][0],n1=edge2node[edgeNum][1] ;

          // Compute the cross-diffusion coefficients associated with the left
          // and right cells of this face-edge pair.
          real crossDiffusionCoefficientLeft=dot(dASumLeft[face][i],
            dualArea[edgeNum]) ;
          real crossDiffusionCoefficientRight=dot(dASumRight[face][i],
            dualArea[edgeNum]) ;

          // Assemble cross-diffusion term to the source.
          vect3d tempLeft(cell_s[cl[face]]+face_s[face]) ;
          vect3d tempRight(cell_s[cr[face]]+face_s[face]) ;
          vect3d crossDiffusionTerm=gridMoverFactor0[edgeNum]*
            (crossDiffusionCoefficientLeft*tempLeft+
            crossDiffusionCoefficientRight*tempRight) ;

          // Assemble cross-dilatation term to the source.
          vect3d crossDilatationTerm=(dot(tempLeft,dASumLeft[face][i])+
            dot(tempRight,dASumRight[face][i]))*scaledDualArea[edgeNum] ;

          sSourceTerm[n0]+=(crossDiffusionTerm+crossDilatationTerm) ;
          sSourceTerm[n1]-=(crossDiffusionTerm+crossDilatationTerm) ;

        }
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CrossDiffusionDilatationApplyInterior>
    registerCrossDiffusionDilatationApplyInterior ;

  // Add the boundary face cross-diffusion and cross-dilatation terms to the
  // nodal source terms. This includes two contributions, one from the interior
  // facet and one from the boundary facet. Checked throughly on 11/19/2008.
  class CrossDiffusionDilatationApplyBoundary : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_Map ci ;
      const_multiMap face2edge ;
      const_MapVec<2> edge2node ;
      const_param<real> nu ;
      const_store<real> gridMoverFactor0 ;
      const_store<vect3d> dualArea,scaledDualArea ;
      const_store<real> diffusionVolume ;
      const_store<vect3d> node_s ;
      const_store<vect3d> face_s ;
      const_store<vect3d> cell_s ;
      const_store<vector<vect3d> > dASum,dABoundary ;
      store<vect3d> sSourceTerm ;
    public:

      // Define input and output.
      CrossDiffusionDilatationApplyBoundary() {
        name_store("ci{n,itg}",ci) ;
        name_store("face2edge{n,itg}",face2edge) ;
        name_store("edge2node{n,itg}",edge2node) ;
        name_store("gridMoverNu{n,itg}",nu) ;
        name_store("gridMoverFactor0{n}",gridMoverFactor0) ;
        name_store("dualArea{n}",dualArea) ;
        name_store("scaledDualArea{n}",scaledDualArea) ;
        name_store("diffusionVolume{n}",diffusionVolume) ;
        name_store("node_sStar{n,itg}",node_s) ;
        name_store("face(node_sStar){n,itg}",face_s) ;
        name_store("cell(node_sStar){n,itg}",cell_s) ;
        name_store("dASum{n}",dASum) ;
        name_store("dABoundary{n}",dABoundary) ;
        name_store("sSourceTerm{n,itg}",sSourceTerm) ;
        input("gridMoverNu{n,itg}") ;
        input("face2edge{n,itg}->edge2node{n,itg}->node_sStar{n,itg}");
        input("face2edge{n,itg}->gridMoverFactor0{n}") ;
        input("face2edge{n,itg}->(dualArea{n},diffusionVolume{n})") ;
        input("face2edge{n,itg}->scaledDualArea{n}") ;
        input("face(node_sStar){n,itg},ci{n,itg}->cell(node_sStar){n,itg}") ;
        input("dASum{n},dABoundary{n}") ;
        output("face2edge{n,itg}->edge2node{n,itg}->sSourceTerm{n,itg}") ;
        constraint("boundaryFaces{n,itg}") ;
      }

      // Contributions to nodes for a given face.
      void calculate(Entity face) {
        unsigned int numEdge=face2edge.num_elems(face) ;
        for(unsigned int i=0;i<numEdge;++i){

          // Get node positions and edge center.
          unsigned int edgeNum=face2edge[face][i] ;
          unsigned int n0=edge2node[edgeNum][0],n1=edge2node[edgeNum][1] ;

          // Compute the cross-diffusion coefficient associated with the
          // boundary cell of this face-edge pair.
          real crossDiffusionCoefficient=dot(dASum[face][i],dualArea[edgeNum]) ;

          // Assemble cross-diffusion term.
          vect3d temp(cell_s[ci[face]]+face_s[face]) ;
          vect3d crossDiffusionTerm=gridMoverFactor0[edgeNum]*
            crossDiffusionCoefficient*temp ;

          // Assemble cross-dilatation term.
          vect3d crossDilatationTerm=dot(temp,dASum[face][i])*
            scaledDualArea[edgeNum] ;

          // Assemble the boundary part of the cross diffusion contribution.
          vect3d temp2(0.5*(node_s[n0]+node_s[n1])+face_s[face]) ;
          vect3d boundaryCrossDiffusionTerm=temp2*dot(dABoundary[face][i],
            dualArea[edgeNum])*gridMoverFactor0[edgeNum] ;

          // Assemble the boundary part of the cross dilation contribution.
          vect3d boundaryCrossDilatationTerm=dot(temp2,dABoundary[face][i])*
            scaledDualArea[edgeNum] ;

          // Add contributions to source term.
          sSourceTerm[n0]+=(crossDiffusionTerm+crossDilatationTerm+
            boundaryCrossDiffusionTerm+boundaryCrossDilatationTerm) ;
          sSourceTerm[n1]-=(crossDiffusionTerm+crossDilatationTerm+
            boundaryCrossDiffusionTerm+boundaryCrossDilatationTerm) ;

          // Diffusion flux through the boundary set to zero and makes no
          // contribution to nodes 0 or 1. This is equivalent to setting all
          // components of the displacement gradient tensor to zero at the
          // boundary nodes. Consequently, nodes with symmetry displacement BCs
          // are computed under stress free boundary condition in the plane
          // of symmetry and zero displacement normal to the plane of symmetry.
          // Dilation is also set to zero at the boundary faces and makes no
          // contribution to nodes 0 or 1. Original comment by MPH 1/14/05.

        }
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CrossDiffusionDilatationApplyBoundary>
    registerCrossDiffusionDilatationApplyBoundary ;

  // Apply rule for adding the rigid-body rotation diffusion term to the
  // source term. This term helps ensure absolution rigid-body rotation of
  // cells near the surface of moving bodies.
  class RigidBodyRotationDiffusionToDisplacementSourceTerm : public apply_rule
  <store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_param<real> gridMoverKappa ;
      const_MapVec<2> edge2node ;
      const_store<real> gridMoverFactor0 ;
      const_store<vect3d> dualArea ;
      const_store<real> diffusionVolume ;
      const_store<tens3d> node_sGrad ;
      store<vect3d> sSourceTerm ;
    public:

      // Define input and output.
      RigidBodyRotationDiffusionToDisplacementSourceTerm() {
        name_store("gridMoverKappa{n,itg}",gridMoverKappa) ;
        name_store("gridMoverFactor0{n}",gridMoverFactor0) ;
        name_store("edge2node{n,itg}",edge2node) ;
        name_store("dualArea{n}",dualArea) ;
        name_store("diffusionVolume{n}",diffusionVolume) ;
        name_store("nodeGrad(node_sStar){n,itg}",node_sGrad) ;
        name_store("sSourceTerm{n,itg}",sSourceTerm) ;
        input("gridMoverKappa{n,itg},gridMoverFactor0{n},dualArea{n}") ;
        input("diffusionVolume{n}") ;
        input("edge2node{n,itg}->nodeGrad(node_sStar){n,itg}") ;
        output("edge2node{n,itg}->sSourceTerm{n,itg}") ;
        constraint("edges{n,itg}") ;
      }

      // Increment the main coefficient for the nodes attached to the edge.
      void calculate(Entity edge) {
        Entity n0=edge2node[edge][0],n1=edge2node[edge][1] ;
        tens3d sGrad0=Transpose(node_sGrad[n0]) ;
        sGrad0.x.x=sGrad0.y.y=sGrad0.z.z=0.0 ;
        tens3d sGrad1=Transpose(node_sGrad[n1]) ;
        sGrad1.x.x=sGrad1.y.y=sGrad1.z.z=0.0 ;
        vect3d temp=-0.5*dotTemp(sGrad0+sGrad1,dualArea[edge])*
          gridMoverFactor0[edge]*(*gridMoverKappa) ;
        sSourceTerm[n0]+=temp ; sSourceTerm[n1]-=temp ;
      }

      // Loop over edges.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<RigidBodyRotationDiffusionToDisplacementSourceTerm>
    registerRigidBodyRotationDiffusionToDisplacementSourceTerm ;

  // Rule to compute the diagonal term for the linear system.
  class DisplacementMatrixDiagonal : public pointwise_rule {
    private:
      const_param<real> gridMoverRelaxation ;
      const_store<real> sMainCoefficient ;
      store<real> D ;
    public:

      // Define input and output.
      DisplacementMatrixDiagonal() {
        name_store("gridMoverRelaxation{n}",gridMoverRelaxation) ;
        name_store("sMainCoefficient{n}",sMainCoefficient) ;
        name_store("sStar_D{n}",D) ;
        input("gridMoverRelaxation{n},sMainCoefficient{n}") ;
        output("sStar_D{n}") ;
        constraint("nodes{n}") ;
      }

      // Set coefficient for a single node.
      void calculate(Entity node) {
        D[node]=sMainCoefficient[node]/(*gridMoverRelaxation) ;
      }

      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DisplacementMatrixDiagonal> registerDisplacementMatrixDiagonal ;

  // Note that we do not need a separate rule to handle sXZero, sYZero and
  // sZZero nodes because the rhs is zero and any non-zero main diagonal
  // coefficient will do. We are guaranteed by the formulation to have a
  // non-zero main diagonal for all nodes on symmetry boundaries.

  // Priority rule to compute the diagonal term for the linear system for
  // boundary nodes where the displacement is specified.
  class DisplacementMatrixDiagonalSpecified : public pointwise_rule {
    private:
      store<real> D ;
    public:

      // Define input and output.
      DisplacementMatrixDiagonalSpecified() {
        name_store("priority::sStar_D{n}",D) ;
        output("priority::sStar_D{n}") ;
        constraint("boundaryDisplacement{n}") ;
      }

      // Set coefficient for a single node.
      void calculate(Entity node) { D[node]=1.0 ; }

      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DisplacementMatrixDiagonalSpecified>
    registerDisplacementMatrixDiagonalSpecified ;

  // Default edge coefficient multiplier.
  class EdgeCoefficientMultiplierDefault : public pointwise_rule {
    private:
      store<real> edgeCoefficientMultiplier ;
    public:

      // Define input and output.
      EdgeCoefficientMultiplierDefault() {
        name_store("edgeCoefficientMultiplier",edgeCoefficientMultiplier) ;
        output("edgeCoefficientMultiplier") ;
        constraint("nodes") ;
      }

      // Set the value for a node.
      void calculate(Entity node) { edgeCoefficientMultiplier[node]=1.0 ; }

      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<EdgeCoefficientMultiplierDefault>
    registerEdgeCoefficientMultiplierDefault ;

  // Priority edge coefficient multiplier.
  class EdgeCoefficientMultiplierPriority : public pointwise_rule {
    private:
      store<real> edgeCoefficientMultiplier ;
    public:

      // Define input and output.
      EdgeCoefficientMultiplierPriority() {
        name_store("boundaryDisplacement::edgeCoefficientMultiplier",
          edgeCoefficientMultiplier) ;
        output("boundaryDisplacement::edgeCoefficientMultiplier") ;
        constraint("boundaryDisplacement") ;
      }

      // Set the value for a node.
      void calculate(Entity node) { edgeCoefficientMultiplier[node]=0.0 ; }

      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<EdgeCoefficientMultiplierPriority>
    registerEdgeCoefficientMultiplierPriority ;

  // Rule to compute the edge coefficient that relates the two nodes. Checked
  // thoroughly on 11/19/2008.
  class DisplacementMatrixEdgeCoefficient : public pointwise_rule {
    private:
      const_MapVec<2> edge2node ;
      const_store<real> gridMoverFactor0 ;
      const_store<vect3d> dualArea ;
      const_store<real> diffusionVolume ;
      const_store<real> edgeCoefficientMultiplier ;
      storeVec<real> sStar_E ;
    public:

      // Define input and output.
      DisplacementMatrixEdgeCoefficient() {
        name_store("edge2node{n}",edge2node) ;
        name_store("gridMoverFactor0{n}",gridMoverFactor0) ;
        name_store("dualArea{n}",dualArea) ;
        name_store("diffusionVolume{n}",diffusionVolume) ;
        name_store("edgeCoefficientMultiplier{n}",edgeCoefficientMultiplier) ;
        name_store("sStar_E{n}",sStar_E) ;
        input("gridMoverFactor0{n},dualArea{n},diffusionVolume{n}") ;
        input("edge2node{n}->edgeCoefficientMultiplier{n}") ;
        output("sStar_E{n}") ;
        constraint("edges{n}") ;
      }

      // Set the coefficients for the nodes attached to the edge.
      void calculate(Entity edge) {
        unsigned int n0=edge2node[edge][0],n1=edge2node[edge][1] ;
        real temp=dot(dualArea[edge],dualArea[edge])*gridMoverFactor0[edge]/
          (3.0*diffusionVolume[edge]) ;
        sStar_E[edge][0]=-temp*edgeCoefficientMultiplier[n0] ;
        sStar_E[edge][1]=-temp*edgeCoefficientMultiplier[n1] ;
      }

      // Loop over edges.
      virtual void compute(const sequence &seq) {
        sStar_E.setVecSize(2) ; do_loop(seq,this) ;
      }
  } ;

  register_rule<DisplacementMatrixEdgeCoefficient>
    registerDisplacementMatrixEdgeCoefficient ;

  // Priority rule to zero out the edge coefficients for boundary edges that
  // do not contain any symmetry nodes. Note that all edges common to a
  // symmetry boundary and another boundary type are marked non-symmetry since
  // the nodes attached to these edges will have the displacement specified.
  class DisplacementMatrixEdgeCoefficientNoSymmetry : public pointwise_rule {
    private:
      storeVec<real> sStar_E ;
    public:

      // Define input and output.
      DisplacementMatrixEdgeCoefficientNoSymmetry() {
        name_store("boundaryNoSymmetryEdges::sStar_E{n}",sStar_E) ;
        output("boundaryNoSymmetryEdges::sStar_E{n}") ;
        constraint("boundaryNoSymmetryEdges{n}") ;
      }

      // Set the coefficients to zero.
      void calculate(Entity edge) { sStar_E[edge][0]=sStar_E[edge][1]=0.0 ; }

      // Loop over edges.
      virtual void compute(const sequence &seq) {
        sStar_E.setVecSize(2) ; do_loop(seq,this) ;
      }
  } ;

//register_rule<DisplacementMatrixEdgeCoefficientNoSymmetry>
//  registerDisplacementMatrixEdgeCoefficientNoSymmetry ;

  // Rule to compute the right-hand side for the linear system.
  class DisplacementRHS : public pointwise_rule {
    private:
      const_param<real> twoDimensionFactor ;
      const_param<real> gridMoverRelaxation ;
      const_store<vect3d> node_sStar ;
      const_store<real> sMainCoefficient ;
      const_store<vect3d> sSourceTerm ;
      store<vect3d> B ;
    public:
                                                                                
      // Define input and output.
      DisplacementRHS() {
        name_store("twoDimensionFactor{n,itg}",twoDimensionFactor) ;
        name_store("gridMoverRelaxation{n,itg}",gridMoverRelaxation) ;
        name_store("node_sStar{n,itg}",node_sStar) ;
        name_store("sMainCoefficient{n,itg}",sMainCoefficient) ;
        name_store("sSourceTerm{n,itg}",sSourceTerm) ;
        name_store("sStar_B{n,itg}",B) ;
        input("gridMoverRelaxation{n,itg},node_sStar{n,itg}") ;
        input("sMainCoefficient{n,itg}") ;
        input("sSourceTerm{n,itg},twoDimensionFactor{n,itg}") ;
        output("sStar_B{n,itg}") ;
        constraint("nodes{n,itg}") ;
      }
                                                                                
      // Add relaxation for a single node.
      void calculate(Entity node) {
        B[node]=sSourceTerm[node]+(1.0-(*gridMoverRelaxation))*
          sMainCoefficient[node]*node_sStar[node]/(*gridMoverRelaxation) ;
        B[node].z*=(*twoDimensionFactor) ;
      }
                                                                                
      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<DisplacementRHS> registerDisplacementRHS ;

  // Priority rule to constrain the x-displacement of boundary nodes with
  // the "sXZero" option.
  class DisplacementRHSSXZero : public pointwise_rule {
    private:
      const_param<real> twoDimensionFactor ;
      const_param<real> gridMoverRelaxation ;
      const_store<vect3d> node_sStar ;
      const_store<real> sMainCoefficient ;
      const_store<vect3d> sSourceTerm ;
      store<vect3d> B ;
    public:
                                                                                
      // Define input and output.
      DisplacementRHSSXZero() {
        name_store("twoDimensionFactor{n,itg}",twoDimensionFactor) ;
        name_store("gridMoverRelaxation{n,itg}",gridMoverRelaxation) ;
        name_store("node_sStar{n,itg}",node_sStar) ;
        name_store("sMainCoefficient{n,itg}",sMainCoefficient) ;
        name_store("sSourceTerm{n,itg}",sSourceTerm) ;
        name_store("priority::sStar_B{n,itg}",B) ;
        input("gridMoverRelaxation{n,itg},node_sStar{n,itg}") ;
        input("sMainCoefficient{n,itg},sSourceTerm{n,itg}") ;
        input("twoDimensionFactor{n,itg}") ;
        output("priority::sStar_B{n,itg}") ;
        constraint("sXZeroNodes{n,itg}") ;
      }
                                                                                
      // Add relaxation for a single node.
      void calculate(Entity node) {
        B[node]=sSourceTerm[node]+(1.0-(*gridMoverRelaxation))*
          sMainCoefficient[node]*node_sStar[node]/(*gridMoverRelaxation) ;
        B[node].x=0.0 ; B[node].z*=(*twoDimensionFactor) ;
      }
                                                                                
      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  //register_rule<DisplacementRHSSXZero> registerDisplacementRHSSXZero ;

  // Priority rule to constrain the y-displacement of boundary nodes with
  // the "sYZero" option.
  class DisplacementRHSSYZero : public pointwise_rule {
    private:
      const_param<real> twoDimensionFactor ;
      const_param<real> gridMoverRelaxation ;
      const_store<vect3d> node_sStar ;
      const_store<real> sMainCoefficient ;
      const_store<vect3d> sSourceTerm ;
      store<vect3d> B ;
    public:
                                                                                
      // Define input and output.
      DisplacementRHSSYZero() {
        name_store("twoDimensionFactor{n,itg}",twoDimensionFactor) ;
        name_store("gridMoverRelaxation{n,itg}",gridMoverRelaxation) ;
        name_store("node_sStar{n,itg}",node_sStar) ;
        name_store("sMainCoefficient{n,itg}",sMainCoefficient) ;
        name_store("sSourceTerm{n,itg}",sSourceTerm) ;
        name_store("priority::sStar_B{n,itg}",B) ;
        input("gridMoverRelaxation{n,itg},node_sStar{n,itg}") ;
        input("sMainCoefficient{n,itg},sSourceTerm{n,itg}") ;
        input("twoDimensionFactor{n,itg}") ;
        output("priority::sStar_B{n,itg}") ;
        constraint("sYZeroNodes{n,itg}") ;
      }
                                                                                
      // Add relaxation for a single node.
      void calculate(Entity node) {
        B[node]=sSourceTerm[node]+(1.0-(*gridMoverRelaxation))*
          sMainCoefficient[node]*node_sStar[node]/(*gridMoverRelaxation) ;
        B[node].y=0.0 ; B[node].z*=(*twoDimensionFactor) ;
      }
                                                                                
      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  //register_rule<DisplacementRHSSYZero> registerDisplacementRHSSYZero ;

  // Priority rule to constrain the z-displacement of boundary nodes with
  // the "sZZero" option.
  class DisplacementRHSSZZero : public pointwise_rule {
    private:
      const_param<real> gridMoverRelaxation ;
      const_store<vect3d> node_sStar ;
      const_store<real> sMainCoefficient ;
      const_store<vect3d> sSourceTerm ;
      store<vect3d> B ;
    public:
                                                                                
      // Define input and output.
      DisplacementRHSSZZero() {
        name_store("gridMoverRelaxation{n,itg}",gridMoverRelaxation) ;
        name_store("node_sStar{n,itg}",node_sStar) ;
        name_store("sMainCoefficient{n,itg}",sMainCoefficient) ;
        name_store("sSourceTerm{n,itg}",sSourceTerm) ;
        name_store("priority::sStar_B{n,itg}",B) ;
        input("gridMoverRelaxation{n,itg},node_sStar{n,itg}") ;
        input("sMainCoefficient{n,itg},sSourceTerm{n,itg}") ;
        output("priority::sStar_B{n,itg}") ;
        constraint("sZZeroNodes{n,itg}") ;
      }
                                                                                
      // Add relaxation for a single node.
      void calculate(Entity node) {
        B[node]=sSourceTerm[node]+(1.0-(*gridMoverRelaxation))*
          sMainCoefficient[node]*node_sStar[node]/(*gridMoverRelaxation) ;
        B[node].z=0.0 ;
      }
                                                                                
      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
    register_rule<DisplacementRHSSZZero> registerDisplacementRHSSZZero ;

  // Priority rule to compute the rhs term for the linear system for
  // boundary nodes where the displacement is completely specified. This rule
  // overrides the rules above for nodes that they have in common.
  class DisplacementRHSSpecified : public pointwise_rule {
    private:
      const_store<vect3d> node_s_b ;
      store<vect3d> B ;
    public:

      // Define input and output.
      DisplacementRHSSpecified() {
        name_store("node_s_b{n,itg}",node_s_b) ;
        name_store("priority::priority::sStar_B{n,itg}",B) ;
        input("node_s_b{n,itg}") ;
        output("priority::priority::sStar_B{n,itg}") ;
        constraint("boundaryDisplacement{n,itg}") ;
      }

      // Set source term for a single node.
      void calculate(Entity node) { B[node]=node_s_b[node] ; }

      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DisplacementRHSSpecified> registerDisplacementRHSSpecified ;

//-----------------------------------------------------------------------------
// Rules for computing the residual of the displacement equation.

  // Rule to initialize the residual.
  class InitializeDisplacementResidual : public unit_rule {
    private:
      store<vect3d> sResidualTemp ;
    public:

      // Define input and output.
      InitializeDisplacementResidual() {
        name_store("sResidualTemp",sResidualTemp) ;
        output("sResidualTemp") ;
        constraint("nodes") ;
      }

      // Initialize the residual for a single node.
      void calculate(Entity node) { sResidualTemp[node]=vect3d(0.0,0.0,0.0) ; }

      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeDisplacementResidual>
    registerInitializeDisplacementResidual ;

  // Adds source and diagonal terms to residual.
  class ComputeDisplacementResidualOne : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_store<real> D ;
      const_store<vect3d> node_s ;
      const_store<vect3d> B ;
      store<vect3d> sResidualTemp ;
    public:

      // Define input and output.
      ComputeDisplacementResidualOne() {
        name_store("sStar_D",D) ;
        name_store("node_sStar",node_s) ;
        name_store("sStar_B",B) ;
        name_store("sResidualTemp",sResidualTemp) ;
        input("sStar_D,node_sStar,sStar_B") ;
        output("sResidualTemp") ;
        constraint("nodes") ;
      }

      // Add the source and diagonal terms.
      void calculate(Entity node) {
        sResidualTemp[node]+=B[node]-D[node]*node_s[node] ;
      }

      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeDisplacementResidualOne>
    registerComputeDisplacementResidualOne ;

  // Adds edge contributions to residual.
  class ComputeDisplacementResidualTwo : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_MapVec<2> edge2node ;
      const_store<vect3d> node_s ;
      const_storeVec<real> sStar_E ;
      store<vect3d> sResidualTemp ;
    public:

      // Define input and output.
      ComputeDisplacementResidualTwo() {
        name_store("edge2node",edge2node) ;
        name_store("node_sStar",node_s) ;
        name_store("sStar_E",sStar_E) ;
        name_store("sResidualTemp",sResidualTemp) ;
        input("edge2node->node_sStar,sStar_E") ;
        output("edge2node->sResidualTemp") ;
        constraint("edges") ;
      }

      // Add the edge neighbor contribution to the residual for each of the two
      // nodes on either side of the edge.
      void calculate(Entity edge) {
        Entity n0=edge2node[edge][0],n1=edge2node[edge][1] ;
        sResidualTemp[n0]-=sStar_E[edge][0]*node_s[n1] ;
        sResidualTemp[n1]-=sStar_E[edge][1]*node_s[n0] ;
      }

      // Loop over edges.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeDisplacementResidualTwo>
    registerComputeDisplacementResidualTwo ;

  // Default rule for the residual.
  class ComputeFinalDisplacementResidualDefault : public pointwise_rule {
    private:
      const_store<vect3d> sResidualTemp ;
      store<vect3d> sResidual ;
    public:

      // Define input and output.
      ComputeFinalDisplacementResidualDefault() {
        name_store("sResidualTemp",sResidualTemp) ;
        name_store("sResidual",sResidual) ;
        input("sResidualTemp") ;
        output("sResidual") ;
        constraint("nodes") ;
      }

      // Copy the value.
      void calculate(Entity node) { sResidual[node]=sResidualTemp[node] ; }

      // Loop over edges.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeFinalDisplacementResidualDefault>
    registerComputeFinalDisplacementResidualDefault ;

  // Rule for nodes constrained in x-direction.
  class ComputeFinalDisplacementResidualSXZero : public pointwise_rule {
    private:
      const_store<vect3d> sResidualTemp ;
      store<vect3d> sResidual ;
    public:

      // Define input and output.
      ComputeFinalDisplacementResidualSXZero() {
        name_store("sResidualTemp",sResidualTemp) ;
        name_store("priority::sResidual",sResidual) ;
        input("sResidualTemp") ;
        output("priority::sResidual") ;
        constraint("sXZeroNodes") ;
      }

      // Copy the value.
      void calculate(Entity node) {
        sResidual[node]=sResidualTemp[node] ; sResidual[node].x=0.0 ;
      }

      // Loop over edges.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeFinalDisplacementResidualSXZero>
    registerComputeFinalDisplacementResidualSXZero ;

  // Rule for nodes constrained in y-direction.
  class ComputeFinalDisplacementResidualSYZero : public pointwise_rule {
    private:
      const_store<vect3d> sResidualTemp ;
      store<vect3d> sResidual ;
    public:

      // Define input and output.
      ComputeFinalDisplacementResidualSYZero() {
        name_store("sResidualTemp",sResidualTemp) ;
        name_store("priority::sResidual",sResidual) ;
        input("sResidualTemp") ;
        output("priority::sResidual") ;
        constraint("sYZeroNodes") ;
      }

      // Copy the value.
      void calculate(Entity node) {
        sResidual[node]=sResidualTemp[node] ; sResidual[node].y=0.0 ;
      }

      // Loop over edges.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeFinalDisplacementResidualSYZero>
    registerComputeFinalDisplacementResidualSYZero ;

  // Rule for nodes constrained in z-direction.
  class ComputeFinalDisplacementResidualSZZero : public pointwise_rule {
    private:
      const_store<vect3d> sResidualTemp ;
      store<vect3d> sResidual ;
    public:

      // Define input and output.
      ComputeFinalDisplacementResidualSZZero() {
        name_store("sResidualTemp",sResidualTemp) ;
        name_store("priority::sResidual",sResidual) ;
        input("sResidualTemp") ;
        output("priority::sResidual") ;
        constraint("sZZeroNodes") ;
      }

      // Copy the value.
      void calculate(Entity node) {
        sResidual[node]=sResidualTemp[node] ; sResidual[node].z=0.0 ;
      }

      // Loop over edges.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeFinalDisplacementResidualSZZero>
    registerComputeFinalDisplacementResidualSZZero ;

  // Rule for specified displacment.
  class ComputeFinalDisplacementResidualBoundaryDisplacement : public
  pointwise_rule {
    private:
      store<vect3d> sResidual ;
    public:

      // Define input and output.
      ComputeFinalDisplacementResidualBoundaryDisplacement() {
        name_store("priority::priority::sResidual",sResidual) ;
        output("priority::priority::sResidual") ;
        constraint("boundaryDisplacement") ;
      }

      // Copy the value.
      void calculate(Entity node) { sResidual[node]=vect3d(0.0,0.0,0.0) ; }

      // Loop over edges.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeFinalDisplacementResidualBoundaryDisplacement>
    registerComputeFinalDisplacementResidualBoundaryDisplacement ;

  // Rule to initialize the total displacement residual.
  class InitializeTotalDisplacementResidual : public unit_rule {
    private:
      param<VectorResidual> sResidualData ;
    public:
                                                                                
      // Define input and output.
      InitializeTotalDisplacementResidual() {
        name_store("sResidualData",sResidualData) ;
        output("sResidualData") ;
        constraint("nodes") ;
      }
                                                                                
      // Initialize the residual.
      virtual void compute(const sequence &seq) {
        *sResidualData=VectorResidual() ;
      }
  } ;
                                                                                
  register_rule<InitializeTotalDisplacementResidual>
    registerInitializeTotalDisplacementResidual ;

  // Rule to compute the total displacement residual.
  class ComputeTotalDisplacementResidual : public apply_rule
  <param<VectorResidual>,VectorResidualJoin> {
    private:
      const_store<vect3d> sResidual ;
      const_store<vect3d> pos ;
      param<VectorResidual> sResidualData ;
    public:
                                                                                
      // Define input and output.
      ComputeTotalDisplacementResidual() {
        name_store("sResidual",sResidual) ;
        name_store("pos",pos) ;
        name_store("sResidualData",sResidualData) ;
        input("sResidual,pos") ;
        output("sResidualData") ;
        constraint("nodes") ;
      }
                                                                                
      // Add the contribution to the residual for a single node.
      void calculate(Entity node) {
        VectorResidual temp ; temp.maxResidual=sResidual[node] ;
        temp.totalResidual=vect3d(abs(sResidual[node].x),
          abs(sResidual[node].y),abs(sResidual[node].z)) ;
        temp.maxResidualLocation=pos[node] ;
        join(*sResidualData,temp) ;
      }
                                                                                
      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<ComputeTotalDisplacementResidual>
    registerComputeTotalDisplacementResidual ;

//-----------------------------------------------------------------------------
// Rules to determine when iteration is finished for node displacement.

  // Build rule.
  class NodeIterationFinishedBuild : public singleton_rule {
    private:
      param<bool> nodeIterationFinished ;
    public:

      // Define input and output.
      NodeIterationFinishedBuild() {
        name_store("nodeIterationFinished{n,itg=-1}",nodeIterationFinished) ;
        output("nodeIterationFinished{n,itg=-1}") ;
//      constraint("UNIVERSE,gridMover,gridMotionTimeDependent") ;
        constraint("UNIVERSE{n},gridMover{n},gridMotionTimeDependent{n}") ;
      }

      // Set to false.
      void compute(const sequence &seq) { *nodeIterationFinished=false ; }
  } ;

  register_rule<NodeIterationFinishedBuild> registerNodeIterationFinishedBuild ;

  // Class to determine if iteration is finished for node displacement.
  // NOTE: Loci refuses to promote this rule to nn loop so we have had to
  // explicitly include the n.
  class CheckNodeIterationFinished : public singleton_rule {
    private:
      const_param<int> itg ;
      const_param<int> gridMoverMaxIterationsPerTimeStep ;
      const_param<real> convergenceTolerance ;
      const_param<VectorResidual> sResidualData ;
      param<bool> nodeIterationFinished ;
    public:

      // Define input and output.
      CheckNodeIterationFinished() {
        name_store("$itg{n,itg}",itg) ;
        name_store("gridMoverMaxIterationsPerTimeStep{n,itg}",
          gridMoverMaxIterationsPerTimeStep) ;
        name_store("convergenceTolerance{n,itg}",convergenceTolerance) ;
        name_store("sResidualData{n,itg}",sResidualData) ;
        name_store("nodeIterationFinished{n,itg}",nodeIterationFinished) ;
        input("$itg{n,itg},gridMoverMaxIterationsPerTimeStep{n,itg}") ;
        input("convergenceTolerance{n,itg},sResidualData{n,itg}") ;
        output("nodeIterationFinished{n,itg}") ;
        constraint("sResidualData{n,itg},gridMover{n,itg}") ;
        constraint("gridMotionTimeDependent{n,itg}") ;
      }

      // Check if iteration is finished.
      void compute(const sequence &seq) {

        // Set the output format and write out the node displacement residuals.
        if(Loci::MPI_rank==0){
          cout.setf(ios::scientific,ios::floatfield) ; cout.precision(6) ;
          cout <<"RS: " << *itg << " " << sResidualData->totalResidual
            << endl ;
        }

        // See if we are converged.
        *nodeIterationFinished=(*itg==*gridMoverMaxIterationsPerTimeStep-1 ||
          (sResidualData->totalResidual.x<*convergenceTolerance &&
          sResidualData->totalResidual.y<*convergenceTolerance &&
          sResidualData->totalResidual.z<*convergenceTolerance)) ;

        // TEMPORARY.
        if(sResidualData->totalResidual.x==0.0 && sResidualData->
        totalResidual.y==0.0 && sResidualData->totalResidual.z==0.0)
          *nodeIterationFinished=false ;
      }
  } ;

  register_rule<CheckNodeIterationFinished> registerCheckNodeIterationFinished ;

}





