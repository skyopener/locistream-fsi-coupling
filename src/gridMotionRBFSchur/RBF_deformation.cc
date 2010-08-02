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
#include <mpi.h>
                                                                                
// StreamUns includes.
#include "const.h"
#include "move.h"
#include "residual.h"
#include "sciTypes.h"
                                                                                
namespace streamUns {

//---------------------------------------------------------------------
// Define some useful constraints
//
  // Rule for creating constraint for all the nodes.
//  class NodeConstraint : public pointwise_rule {
//    private:
//      const_multiMap face2node ;
//      store<bool> nodes ;
//    public:
//
//      // Define input and output.
//      NodeConstraint() {
//        name_store("face2node",face2node) ;
//        name_store("nodes",nodes) ;
//        input("face2node") ;
//        output("face2node->nodes") ;
//        constraint("faces") ;
//      }
//
//      // Do nothing for a face since this rule is only used to establish
//      // a constraint.
//      void calculate(Entity face) {}
//
//      // Loop over faces.
//      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
//  } ;
//
//  register_rule<NodeConstraint> registerNodeConstraint ;

//---------------------------------------------------------------------
// Rules for printing out node and boundary node positions
//
  // Rule for printing out node positions. Only used as a check.
//  class NodePosition : public pointwise_rule {
//    private:
//      const_store<vect3d> pos ;
//      store<vect3d> nodePosition ;
//      unsigned int numNode ;
//    public:
//
//      // Define input and output.
//      NodePosition() {
//        name_store("pos",pos) ;
//        name_store("nodePosition",nodePosition) ;
//        output("nodePosition") ;
//        constraint("nodes") ;
//      }
//
//      // Do nothing for a node.
//      void calculate(Entity node) { ++numNode ; }
//
//      // Loop over nodes.
//      virtual void compute(const sequence &seq) {
//        numNode=0 ; do_loop(seq,this) ;
//        cout << "numNode: " << numNode << endl ;
//      }
//  } ;
//
//  register_rule<NodePosition> registerNodePosition ;
//
//  // Rule for printing out boundary node positions. Only used as a check.
//  class BoundaryNodePosition : public pointwise_rule {
//    private:
//      const_store<vect3d> pos ;
//      store<vect3d> boundaryNodePosition ;
//      unsigned int numNode ;
//    public:
//
//      // Define input and output.
//      BoundaryNodePosition() {
//        name_store("pos",pos) ;
//        name_store("boundaryNodePosition",boundaryNodePosition) ;
//        output("boundaryNodePosition") ;
//        constraint("boundaryDisplacement") ;
//      }
//
//      // Do nothing for a node.
//      void calculate(Entity node) { ++numNode ; }
//
//      // Loop over nodes.
//      virtual void compute(const sequence &seq) {
//        numNode=0 ; do_loop(seq,this) ;
//        cout << "numBoundaryNode: " << numNode << endl ;
//      }
//  } ;
//
//  register_rule<BoundaryNodePosition> registerBoundaryNodePosition ;

//-----------------------------------------------------------------------------

  // Rule to assemble the RBF RHS matrix in the x-direction 
  // with the "sXZero" option and the specified boundary displacement nodes
  class AssembleRHSRBF : public pointwise_rule {
    private:
      const_store<vect3d> node_s_b ;
      store<vect3d> B ;
    public:
                                                                                
      // Define input and output.
      AssembleRHSRBF() {
        name_store("sStar_B{n}",B) ;
        name_store("node_s_b{n}", node_s_b) ;
        input("node_s_b{n}") ;
        output("sStar_B{n}") ;
        constraint("boundaryDisplacement{n}") ;
      }
                                                                                
      // Add relaxation for a single node.
      void calculate(Entity node) {
		  B[node] = node_s_b[node] ;
      }
                                                                                
      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<AssembleRHSRBF> registerAssembleRHSRBF ;

//-----------------------------------------------------------------------------
// Rules for setting up the PETSC data.

  // Gets the number of boundary nodes assigned on all processes. This rule is
  // definitely not in the Loci spirit since we are collecting data from
  // other processes.
  class RbfGetLocalBoundaryNodeDataUnit : public unit_rule {
    private:
      param<vector<int> > rbfNumBoundaryNode ;
    public:

      // Define input and output.
      RbfGetLocalBoundaryNodeDataUnit() {
        name_store("rbfNumBoundaryNode",rbfNumBoundaryNode) ;
        output("rbfNumBoundaryNode") ;
        constraint("pos") ;
        disable_threading() ;
      }

      // Get the number of nodes for each process.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<RbfGetLocalBoundaryNodeDataUnit> registerRbfGetLocalBoundaryNodeDataUnit ;

  // Gets the number of boundary nodes assigned on all processes. This rule is
  // definitely not in the Loci spirit since we are collecting data from
  // other processes.
  class RbfGetLocalBoundaryNodeDataApply : public apply_rule<param<vector<int> >,Loci::NullOp<vector<int> > > {
    private:
      param<vector<int> > rbfNumBoundaryNode ;
    public:

      // Define input and output.
      RbfGetLocalBoundaryNodeDataApply() {
        name_store("rbfNumBoundaryNode",rbfNumBoundaryNode) ;
        output("rbfNumBoundaryNode") ;
        constraint("boundaryDisplacement") ;
//        constraint("RBFxConstraints") ;
        disable_threading() ;
      }

      // Get the number of nodes for each process.
      virtual void compute(const sequence &seq) {

		const int rank=Loci::MPI_rank ;
		const int p=Loci::MPI_processes ;

		entitySet rbfLocalBoundaryNode ;
        // Get the collection of entities assigned to this processor
        Loci::storeRepP myEntities=Loci::exec_current_fact_db->get_variable
          ("my_entities") ;
        entitySet localEntities=~EMPTY ;
        if(myEntities!=0) localEntities=(*myEntities).domain() ;

        // Get the local nodes.
        rbfLocalBoundaryNode=entitySet(seq) ;

        // Distribute the number of nodes to all processes.
        *rbfNumBoundaryNode=Loci::all_collect_sizes((rbfLocalBoundaryNode).size()) ;
	/*	
		(*rbfLocalNumNodeP4).resize(p);
		(*rbfLocalStart).resize(p);

		(*rbfGlobalNumNode)=0 ;
        for(unsigned int i=0;i<(*rbfNumBoundaryNode).size();++i) (*rbfGlobalNumNode)+=(*rbfNumBoundaryNode)[i] ;
		(*rbfGlobalNumNodeP4) = (*rbfGlobalNumNode) + 4 ;
//		cout << "0:rank, localNumNode, globalNumNode: " << Loci::MPI_rank << ", " << ", " << localNumNode << ", " << globalNumNode << endl;
		
		// include the polynomial coeff. computation
		for(int i=0;i<p;++i) {
			(*rbfLocalNumNodeP4)[i]=((i==p-1)?((*rbfNumBoundaryNode)[i]+4):(*rbfNumBoundaryNode)[i]) ;
			for(int j=0;j<i;++j) (*rbfLocalStart)[i]+=(*rbfNumBoundaryNode)[j] ;
		}
	*/				
//		cout << "4:rank, localStart, localNumNode, globalNumNode: " << Loci::MPI_rank << ", " << localStart << ", " << localNumNode << ", " << globalNumNode << endl;
      }
  } ;

  register_rule<RbfGetLocalBoundaryNodeDataApply> registerRbfGetLocalBoundaryNodeDataApply ;

  // Creates the node-to-row map. Each boundary node now has a unique row index
  class RBFBoundaryNodeToRow : public pointwise_rule {
    private:
      const_param<vector<int> > rbfNumBoundaryNode ;
      store<int> rbfBoundaryNodeToRow ;
    public:

      // Define input and output.
      RBFBoundaryNodeToRow() {
        name_store("rbfNumBoundaryNode",rbfNumBoundaryNode) ;
        name_store("rbfBoundaryNodeToRow",rbfBoundaryNodeToRow) ;
        input("rbfNumBoundaryNode") ;
        output("rbfBoundaryNodeToRow") ;
        constraint("boundaryDisplacement") ;
//        constraint("RBFxConstraints") ;
        disable_threading() ;
      }

      // Set the Petsc row for each node.
      virtual void compute(const sequence & seq) {

        // Compute the row offset for this process.
//        int offset=(*rbfLocalStart)[Loci::MPI_rank] ;
        int offset=0 ;
        for(int i=0;i<Loci::MPI_rank;++i){ offset+=(*rbfNumBoundaryNode)[i] ; }

        // Assign row number.
        sequence::const_iterator nodePtr=seq.begin() ;
        for(int i=0;i<(*rbfNumBoundaryNode)[Loci::MPI_rank];++nodePtr,++i){
          rbfBoundaryNodeToRow[*nodePtr]=offset+i ;
        }
      }
  } ;

  register_rule<RBFBoundaryNodeToRow> registerRBFBoundaryNodeToRow ;
  
//-----------------------------Q---------------------------------------------------------
//

  // Sets up the Q vector.
  class RBFBoundaryNodeAssembleQUnit : public unit_rule { // see interpolateFile.cc
    private:
      const_param<vector<int> > rbfNumBoundaryNode ;
	  const_store<vect3d> pos ;
	  param<vector<real> > rbfQ ;
    public:

      // Define input and output.
      RBFBoundaryNodeAssembleQUnit() {
		name_store("rbfBoundaryNodeQ",rbfQ) ;
        output("rbfBoundaryNodeQ") ;
        constraint("sStar_PETSCLinearSolver") ;
		constraint("pos"); // creates Q for all nodes 
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {}
  } ;
 
  register_rule<RBFBoundaryNodeAssembleQUnit> registerRBFBoundaryNodeAssembleQUnit ;

  // Assemble Q -> unit rule to output blackbox
  class RBFBoundaryNodeAssembleQApply : public apply_rule<param<vector<real> >,Loci::NullOp<vector<real> > > {
    private:
      const_param<vector<int> > rbfNumBoundaryNode ;
	  const_store<vect3d> pos ;
	  param<vector<real> > rbfQ ;
    public:

      // Define input and output.
      RBFBoundaryNodeAssembleQApply() {
        name_store("rbfNumBoundaryNode",rbfNumBoundaryNode) ;
		name_store("pos",pos) ;
        name_store("rbfBoundaryNodeQ", rbfQ) ;
        input("rbfNumBoundaryNode") ;
        input("pos") ;
        output("rbfBoundaryNodeQ") ;
        constraint("sStar_PETSCLinearSolver") ;
		constraint("boundaryDisplacement");
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {

        // Get the number of local and global nodes.
		const int p = Loci::MPI_processes ;
		const int rank = Loci::MPI_rank ;
        int globalNumNode=0 ;
        for(unsigned int i=0;i<(*rbfNumBoundaryNode).size();++i) globalNumNode+=(*rbfNumBoundaryNode)[i] ;
		
		// Compute the row offset for this process.
        int localStart=0 ; 
        for(int i=0;i<rank;++i) localStart+=(*rbfNumBoundaryNode)[i]  ;  // 3 coordinates
		int localNumNode=(*rbfNumBoundaryNode)[rank]; 
					
		vector<int> recvcounts(p,0) ;
		vector<int> displs(p,0) ;
		for(int i=0;i<p;++i) {
			for(int j=0;j<i;++j) displs[i]+=(3*(*rbfNumBoundaryNode)[j]);
			recvcounts[i] = 3*(*rbfNumBoundaryNode)[i] ;
		}
				
		// Allocate rbfQ
		(*rbfQ).resize(3*globalNumNode) ; // number of nodes * 3 coordinates
//		(*rbfQ).resize(3*(*rbfGlobalNumNode)) ; // number of nodes * 3 coordinates

		// Fill in rbfQ
        sequence::const_iterator nodePtr=seq.begin() ;
//		for(int i=localStart; i<localStart+localNum;++nodePtr, ++i) {
//		int counter=(*rbfLocalStart)[rank];
		int counter=localStart;
		for(nodePtr=seq.begin(); nodePtr!=seq.end();++nodePtr, ++counter) { // Q[row,col]=Q[3*row+col], Q=[0,1,2;3,4,5;...]
			(*rbfQ)[3*counter+0] = pos[*nodePtr].x ;
			(*rbfQ)[3*counter+1] = pos[*nodePtr].y ;
			(*rbfQ)[3*counter+2] = pos[*nodePtr].z ;
		}

		// Allgatherv
//		for(int i=0; i<(*rbfQ).size(); ++i) {
//			cout << "before gather: p,r,rbfQ[" << i << "]: " << p << "," << r << ", " << (*rbfQ)[i] << endl ;
//		}
		MPI_Allgatherv(&(*rbfQ)[3*(localStart)], 3*localNumNode, MPI_DOUBLE, &(*rbfQ)[0], &recvcounts[0], &displs[0], MPI_DOUBLE, MPI_COMM_WORLD);
      }
  } ;

  register_rule<RBFBoundaryNodeAssembleQApply> registerRBFBoundaryNodeAssembleQApply ;

}





