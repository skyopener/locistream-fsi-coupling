// Contains rules for computing that gradient of a vect3d at the nodes given
// the nodal values of the variable. Green's theorem is used.

// Standard library includes.
#include <algorithm>
using std::find ;
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

  // Creates the cell-to-node map.
  class CellToNodeMap : public pointwise_rule {
    private:
      const_multiMap upper,lower,boundary_map,face2node ;
      store<vector<int> > cellToNode ;
    public:

      // Define input and output.
      CellToNodeMap() {
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("boundary_map",boundary_map) ;
        name_store("face2node",face2node) ;
        name_store("cellToNode",cellToNode) ;
        input("(upper,lower,boundary_map)->face2node") ;
        output("cellToNode") ;
      }

      // Find the nodes for each cell using face-to-node data.
      void calculate(Entity cell) {

        // Add nodes from upper faces.
        for(const Entity *face=upper[cell].begin();face!=upper[cell].end();
        ++face){
          for(const int *ni=face2node.begin(*face);ni!=face2node.end(*face);
            ++ni) if(find(cellToNode[cell].begin(),cellToNode[cell].end(),
            *ni)==cellToNode[cell].end()) cellToNode[cell].push_back(*ni) ;
        }

        // Add nodes from lower faces. Those previously added by upper will
        // be overwritten.
        for(const Entity *face=lower[cell].begin();face!=lower[cell].end();
        ++face){
          for(const int *ni=face2node.begin(*face);ni!=face2node.end(*face);
            ++ni) if(find(cellToNode[cell].begin(),cellToNode[cell].end(),
            *ni)==cellToNode[cell].end()) cellToNode[cell].push_back(*ni) ;
        }

        // Add nodes from boundary_map faces. Those previously added by upper
        // and lower  will be overwritten.
        for(const Entity *face=boundary_map[cell].begin();face!=
        boundary_map[cell].end();++face){
          for(const int *ni=face2node.begin(*face);ni!=face2node.end(*face);
            ++ni) if(find(cellToNode[cell].begin(),cellToNode[cell].end(),
            *ni)==cellToNode[cell].end()) cellToNode[cell].push_back(*ni) ; 
        }
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CellToNodeMap> registerCellToNodeMap ;

  // Compute a cell value given the values at the nodes. Note that we are
  // keeping the cell to face maps so Loci will give us the range of nodes
  // that we need.
  class CellValueFromNodeValueVector : public pointwise_rule {
    private:
      const_multiMap upper,lower,boundary_map,face2node ;
      const_store<vector<int> > cellToNode ;
      const_store<vect3d> X ;
      store<vect3d> cellX ;
    public:

      // Define input and output.
      CellValueFromNodeValueVector() {
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("boundary_map",boundary_map) ;
        name_store("face2node",face2node) ;
        name_store("cellToNode",cellToNode) ;
        name_store("X",X) ;
        name_store("cell(X)",cellX) ;
        input("cellToNode,(upper,lower,boundary_map)->face2node->X") ;
        output("cell(X)") ;
      }

      // Average the node values for a single cell. The algorithm used here is
      // very inefficient, but it is the best we can do in light of not having
      // cell-to-node pointer data.
      void calculate(Entity cell) {
        vect3d sum=vect3d(0.0,0.0,0.0) ;
        for(vector<int>::const_iterator p=cellToNode[cell].begin();p!=
        cellToNode[cell].end();++p) sum+=X[*p] ;
        cellX[cell]=(1.0/(real)(cellToNode[cell].size()))*sum ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CellValueFromNodeValueVector>
    registerCellValueFromNodeValueVector ;

  // Compute a face value given the values at the nodes.
  class FaceFromNodeValueVector : public pointwise_rule {
    private:
      const_multiMap face2node ;
      const_store<vect3d> X ;
      store<vect3d> faceX ;
    public:
                                                                                
      // Define input and output.
      FaceFromNodeValueVector() {
        name_store("face2node",face2node) ;
        name_store("X",X) ;
        name_store("face(X)",faceX) ;
        input("face2node->X") ;
        output("face(X)") ;
        constraint("faces") ;
      }
                                                                                
      // Compute the displacement of the face center.
      void calculate(Entity face) {
        vect3d sum(0.0,0.0,0.0) ;
        for(const int *ni=face2node.begin(face);ni!=face2node.end(face);++ni)
          sum+=X[*ni] ;
        faceX[face]=(1.0/real(face2node.num_elems(face)))*sum ;
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<FaceFromNodeValueVector> registerFaceFromNodeValueVector ;

  // Rule to initialize the main coefficient.
  class InitializeNodeGradientVector : public unit_rule {
    private:
      store<tens3d> gradX ;
    public:

      // Define input and output.
      InitializeNodeGradientVector() {
        name_store("nodeGrad(X)",gradX) ;
        output("nodeGrad(X)") ;
        constraint("nodes") ;
      }

      // Initialize to zero.
      void calculate(Entity node) {
        gradX[node]=tens3d(vect3d(0.0,0.0,0.0),vect3d(0.0,0.0,0.0),
          vect3d(0.0,0.0,0.0)) ;
      }

      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeNodeGradientVector>
    registerInitializeNodeGradientVector ;

  // Add contributions from interior faces to node gradient.
  class InternalNodeGradientVector : public apply_rule<store<tens3d>,
  Loci::Summation<tens3d> > {
    private:
      const_Map cl,cr ;
      const_multiMap face2node,face2edge ;
      const_MapVec<2> edge2node ;
      const_store<vect3d> pos ;
      const_store<real> dualVolume ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> cellX ;
      const_store<vect3d> X ;
      const_store<vect3d> faceX ;
      const_store<vect3d> faceCenter ;
      store<tens3d> gradX ;
    public:

      // Define input and output.
      InternalNodeGradientVector() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("face2node",face2node) ;
        name_store("face2edge",face2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("pos",pos) ;
        name_store("dualVolume",dualVolume) ;
        name_store("cellcenter",cellCenter) ;
        name_store("cell(X)",cellX) ;
        name_store("X",X) ;
        name_store("face(X)",faceX) ;
        name_store("facecenter",faceCenter) ;
        name_store("nodeGrad(X)",gradX) ;
        input("(cl,cr)->(cellcenter,cell(X))") ;
        input("face2node,facecenter,face(X)") ;
        input("face2edge,face2edge->edge2node->(pos,X,dualVolume)") ;
        output("face2edge->edge2node->nodeGrad(X)") ;
        constraint("internalFaces") ;
      }

      // Distribute contributions from a face.
      void calculate(Entity face) {

        /* Compute the value of X at the face center.
        vect3d sum(0.0,0.0,0.0) ;
        for(const int *ni=face2node.begin(face);ni!=face2node.end(face);++ni)
          sum+=X[*ni] ;
        vect3d faceX=(1.0/real(face2node.num_elems(face)))*sum ;*/

        // Loop over face edges. Each triangular facet composed of the edge
        // center, face center and cell center contributes to the integral.
        unsigned int numEdge=face2edge.num_elems(face) ;
        for(unsigned int i=0;i<numEdge;++i){
          unsigned int edgeNum=face2edge[face][i] ;
          int n0=edge2node[edgeNum][0],n1=edge2node[edgeNum][1] ;
          vect3d pos0=pos[n0],pos1=pos[n1],edgeCenter=0.5*(pos0+pos1) ;
          vect3d dA0=cross(cellCenter[cl[face]]-edgeCenter,faceCenter[face]-
            edgeCenter) ;
          vect3d dA1=cross(cellCenter[cr[face]]-edgeCenter,faceCenter[face]-
            edgeCenter) ;
          vect3d tempX=0.5*(X[n0]+X[n1])+faceX[face] ;
          tens3d gradIncrement=Dyad(tempX+cellX[cl[face]],dA0)-Dyad(tempX+
            cellX[cr[face]],dA1) ;
          if(dot(dA0,pos1-pos0)>0.0){
            gradX[n0]+=gradIncrement/(6.0*dualVolume[n0]) ;
            gradX[n1]-=gradIncrement/(6.0*dualVolume[n1]) ;
          }else{
            gradX[n0]-=gradIncrement/(6.0*dualVolume[n0]) ;
            gradX[n1]+=gradIncrement/(6.0*dualVolume[n1]) ;
          }
        }
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InternalNodeGradientVector> registerInternalNodeGradientVector ;

  // Add contributions from boundary faces to node gradient. This includes a
  // contribution from an interior facet and a facet on the boundary.
  class BoundaryNodeGradientVector : public apply_rule<store<tens3d>,
  Loci::Summation<tens3d> > {
    private:
      const_Map ci ;
      const_multiMap face2node,face2edge ;
      const_MapVec<2> edge2node ;
      const_store<vect3d> pos ;
      const_store<real> dualVolume ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> cellX ;
      const_store<vect3d> X ;
      const_store<vect3d> faceX ;
      const_store<vect3d> faceCenter ;
      store<tens3d> gradX ;
    public:

      // Define input and output.
      BoundaryNodeGradientVector() {
        name_store("ci",ci) ;
        name_store("face2node",face2node) ;
        name_store("face2edge",face2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("pos",pos) ;
        name_store("dualVolume",dualVolume) ;
        name_store("cellcenter",cellCenter) ;
        name_store("cell(X)",cellX) ;
        name_store("X",X) ;
        name_store("face(X)",faceX) ;
        name_store("facecenter",faceCenter) ;
        name_store("nodeGrad(X)",gradX) ;
        input("ci->(cellcenter,cell(X)),face2node,facecenter,face(X)") ;
        input("face2edge,face2edge->edge2node->(pos,X,dualVolume)") ;
        output("face2edge->edge2node->nodeGrad(X)") ;
        constraint("boundaryFaces") ;
      }

      // Distribute contributions from a face.
      void calculate(Entity face) {

        /* Compute the value of X at the face center.
        vect3d sum(0.0,0.0,0.0) ;
        for(const int *ni=face2node.begin(face);ni!=face2node.end(face);++ni)
          sum+=X[*ni] ;
        vect3d faceX=(1.0/real(face2node.num_elems(face)))*sum ;*/

        // Loop over face edges. Each triangular facet composed of the edge
        // center, face center and cell center contributes to the integral.
        unsigned int numEdge=face2edge.num_elems(face) ;
        for(unsigned int i=0;i<numEdge;++i){
          unsigned int edgeNum=face2edge[face][i] ;
          int n0=edge2node[edgeNum][0],n1=edge2node[edgeNum][1] ;
          vect3d pos0=pos[n0],pos1=pos[n1],edgeCenter=0.5*(pos0+pos1) ;

          // Interior facet contribution.
          vect3d dA=cross(cellCenter[ci[face]]-edgeCenter,faceCenter[face]-
            edgeCenter) ;
          vect3d tempX=0.5*(X[n0]+X[n1])+faceX[face] ;
          tens3d gradIncrement=Dyad(tempX+cellX[ci[face]],dA) ;
          if(dot(dA,pos1-pos0)>0.0){
            gradX[n0]+=gradIncrement/(6.0*dualVolume[n0]) ;
            gradX[n1]-=gradIncrement/(6.0*dualVolume[n1]) ;
          }else{
            gradX[n0]-=gradIncrement/(6.0*dualVolume[n0]) ;
            gradX[n1]+=gradIncrement/(6.0*dualVolume[n1]) ;
          }

          // Boundary facet contribution. Different for the two edge nodes.
          vect3d dA0=cross(edgeCenter-pos0,faceCenter[face]-pos0) ;
          vect3d dA1=cross(edgeCenter-pos1,faceCenter[face]-pos1) ;
          vect3d X0=tempX+X[n0],X1=tempX+X[n1] ;
          if(dot(dA0,cellCenter[ci[face]]-edgeCenter)>0.0){
            gradX[n0]-=Dyad(X0,dA0)/(6.0*dualVolume[n0]) ;
            gradX[n1]+=Dyad(X1,dA1)/(6.0*dualVolume[n1]) ;
          }else{
            gradX[n0]+=Dyad(X0,dA0)/(6.0*dualVolume[n0]) ;
            gradX[n1]-=Dyad(X1,dA1)/(6.0*dualVolume[n1]) ;
          }
        }
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryNodeGradientVector> registerBoundaryNodeGradientVector ;

}
