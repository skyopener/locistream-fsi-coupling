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
#include "move.h"
#include "residual.h"
#include "sciTypes.h"
                                                                                
namespace streamUns {

//---------------------------------------------------------------------
// Define some useful constraints
//
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

  // Rule for creating a internal node constraint.
  class InternalNodeConstraint : public pointwise_rule {
    private:
      store<bool> internalNodes ;
    public:

      // Define input and output.
      InternalNodeConstraint() {
        name_store("internalNodes",boundaryNodes) ;
        output("internalNodes") ;
        constraint("nodes,NOT(boundaryNodes)") ;
      }

      // Do nothing for a face since this rule is only used to establish
      // a constraint.
      void calculate(Entity node) {}

      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InternalNodeConstraint> registerInternalNodeConstraint ;

  // Rule to create constraint for RBF-x-direction: specified displacement
  class RbfXConstraintSpec : public pointwise_rule {
    private:
      store<bool> rbfXconstraints ;
    public:

      // Define input and output.
      RbfXConstraintSpec() {
        name_store("rbfXconstraints",rbfXconstraints) ;
        output("rbfXconstraints") ;
        constraint("boundaryDisplacement") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity node) { rbfXconstraints[node]=true ; }

      // Assign flag for a sequence of nodes with specified boundary displacements.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<RbfXConstraintSpec> registerRbfXConstraintSpec ;

  // Rule to create constraint for RBF-x-direction: SXZero
  class RbfXConstraintSym : public pointwise_rule {
    private:
      store<bool> rbfXconstraints ;
    public:

      // Define input and output.
      RbfXConstraintSym() {
        name_store("rbfXconstraints",rbfXconstraints) ;
        output("rbfXconstraints") ;
        constraint("sXZeroNodes") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity node) { rbfXconstraints[node]=true ; }

      // Assign flag for a sequence of nodes with sym.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<RbfXConstraintSym> registerRbfXConstraintSym ;


  // Rule to create constraint for RBF-y-direction: specified displacement
  class RbfYConstraintSpec : public pointwise_rule {
    private:
      store<bool> rbfYconstraints ;
    public:

      // Define input and output.
      RbfXConstraintSpec() {
        name_store("rbfYconstraints",rbfYconstraints) ;
        output("rbfYconstraints") ;
        constraint("boundaryDisplacement") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity node) { rbfYconstraints[node]=true ; }

      // Assign flag for a sequence of nodes with specified boundary displacements.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<RbfYConstraintSpec> registerRbfYConstraintSpec ;

  // Rule to create constraint for RBF-y-direction: SYZero
  class RbfYConstraintSym : public pointwise_rule {
    private:
      store<bool> rbfYconstraints ;
    public:

      // Define input and output.
      RbfYConstraintSym() {
        name_store("rbfYconstraints",rbfYconstraints) ;
        output("rbfYconstraints") ;
        constraint("sYZeroNodes") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity node) { rbfYconstraints[node]=true ; }

      // Assign flag for a sequence of nodes with sym.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<RbfYConstraintSym> registerRbfYConstraintSym ;

  // Rule to create constraint for RBF-z-direction: specified displacement
  class RbfZConstraintSpec : public pointwise_rule {
    private:
      store<bool> rbfZconstraints ;
    public:

      // Define input and output.
      RbfZConstraintSpec() {
        name_store("rbfZconstraints",rbfZconstraints) ;
        output("rbfZconstraints") ;
        constraint("boundaryDisplacement") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity node) { rbfZconstraints[node]=true ; }

      // Assign flag for a sequence of nodes with specified boundary displacements.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<RbfZConstraintSpec> registerRbfZConstraintSpec ;

  // Rule to create constraint for RBF-z-direction: SZZero
  class RbfZConstraintSym : public pointwise_rule {
    private:
      store<bool> rbfZconstraints ;
    public:

      // Define input and output.
      RbfZConstraintSym() {
        name_store("rbfZconstraints",rbfZconstraints) ;
        output("rbfZconstraints") ;
        constraint("sZZeroNodes") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity node) { rbfZconstraints[node]=true ; }

      // Assign flag for a sequence of nodes with sym.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<RbfZConstraintSym> registerRbfZConstraintSym ;

//---------------------------------------------------------------------
// Rules for printing out node and boundary node positions
//
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

  // Rule to assemble the RBF RHS matrix in the x-direction 
  // with the "sXZero" option and the specified boundary displacement nodes
  class AssembleRHSRBFX : public pointwise_rule {
    private:
      store<real> Bx ;
      const_store<vect3d> node_s_b ;
    public:
                                                                                
      // Define input and output.
      AssembleRHSRBFX() {
        name_store("Bx{n}",Bx) ;
        input("node_s_b{n}") ;
        output("Bx{n}") ;
        constraint("rbfXconstraints{n}") ;
      }
                                                                                
      // Add relaxation for a single node.
      void calculate(Entity node) {
		  Bx[node] = node_s_b.x;
      }
                                                                                
      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<AssembleRHSRBFX> registerAssembleRHSRBFX ;
                                                                                

  // Rule to assemble the RBF RHS matrix in the y-direction 
  // with the "sYZero" option and the specified boundary displacement nodes
  class AssembleRHSRBFY : public pointwise_rule {
    private:
      store<real> By ;
      const_store<vect3d> node_s_b ;
    public:
                                                                                
      // Define input and output.
      AssembleRHSRBFY() {
        name_store("By{n}",By) ;
        input("node_s_b{n}") ;
        output("By{n}") ;
        constraint("rbfYconstraints{n}") ;
      }
                                                                                
      // Add relaxation for a single node.
      void calculate(Entity node) {
		  Bz[node] = node_s_b.z;
      }
                                                                                
      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<AssembleRHSRBFY> registerAssembleRHSRBFY ;
                                                                                

  // Rule to assemble the RBF RHS matrix in the z-direction 
  // with the "sZZero" option and the specified boundary displacement nodes
  class AssembleRHSRBFZ : public pointwise_rule {
    private:
      store<real> Bz ;
      const_store<vect3d> node_s_b ;
    public:
                                                                                
      // Define input and output.
      AssembleRHSRBFZ() {
        name_store("Bz{n}",Bz) ;
        input("node_s_b{n}") ;
        output("Bz{n}") ;
        constraint("rbfZconstraints{n}") ;
      }
                                                                                
      // Add relaxation for a single node.
      void calculate(Entity node) {
		  Bz[node] = node_s_b.z;
		  
      }
                                                                                
      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<AssembleRHSRBFZ> registerAssembleRHSRBFZ ;

}





