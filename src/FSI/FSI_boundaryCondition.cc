//-----------------------------------------------------------------------------
// Description: This file contains some modifications to boundary conditions
//   due to grid movement.
//-----------------------------------------------------------------------------
                                                                                
// Standard library includes.
#include <vector>
using std::vector ;
                                                                                
// Loci includes.
#include <Loci.h>
using Loci::Area ;

// StreamUns includes.
#include "const.h"
#include "gridReader/readGrid.h"
#include "FSI_move.h"
#include "residual.h"
#include "sciTypes.h"
#include "varsFileInputs.h"
                                                                                
namespace streamUns {

//-----------------------------------------------------------------------------
// Boundary condition check rule for the new "general" boundary condition type.

	
  // Checks CSD  boundaries.
  class CheckNoslipCSD : public BC_Check {
    public:
      string BoundaryConditions() { return "noslip" ; }
      string VariablesChecked(fact_db &facts) { return "top,bottom,tip,s" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        return true ;
      }
      ostream &ErrorMessage(ostream &s) { return s ; }
  } ;

  register_BC<CheckNoslipCSD>
    registerCheckNoslipCSD ;

  // Checks fixed-pressure outlet boundaries.
  class CheckFixedPressureOutletGridMotion : public BC_Check {
    public:
      string BoundaryConditions() { return "fixedPressureOutlet" ; }
      string VariablesChecked(fact_db &facts) { return "s" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        return true ;
      }
      ostream &ErrorMessage(ostream &s) { return s ; }
  } ;

  register_BC<CheckFixedPressureOutletGridMotion>
    registerCheckFixedPressureOutletGridMotion ;

  // Checks general boundaries.
  class CheckGeneral : public BC_Check {
    private:
      string errorMessage ;
    public:
      CheckGeneral() : errorMessage("") {}
    public:
      string BoundaryConditions() { return "general" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        return true ;
      }
      ostream &ErrorMessage(std::ostream &s) {
        s << errorMessage << endl ; return s;
      }
      string VariablesChecked(fact_db &facts) {
        string s="s" ; return s ;
      }
  } ;
                                                                                
  register_BC<CheckGeneral> registerCheckGeneral ;

  // Checks incompressible inlet boundaries.
  class CheckIncompressibleInletGridMotion : public BC_Check {
    public:
      string BoundaryConditions() { return "incompressibleInlet" ; }
      string VariablesChecked(fact_db &facts) { return "s" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        return true ;
      }
      ostream &ErrorMessage(ostream &s) { return s ; }
  } ;

  register_BC<CheckIncompressibleInletGridMotion>
    registerCheckIncompressibleInletGridMotion ;
/*
  // Checks noslip boundaries.
  class CheckNoslipGridMotion : public BC_Check {
    public:
      string BoundaryConditions() { return "noslip" ; }
      string VariablesChecked(fact_db &facts) { return "s" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        return true ;
      }
      ostream &ErrorMessage(ostream &s) { return s ; }
  } ;

  register_BC<CheckNoslipGridMotion> registerCheckNoslipGridMotion ;
*/
  // Checks symmetry and slip boundaries.
  class CheckSymmetryGridMotion : public BC_Check {
    public:
      string BoundaryConditions() { return "symmetry" ; }
      string VariablesChecked(fact_db &facts) {
        return "sXZero,sYZero,sZZero,s" ;
      }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        return true ;
      }
      ostream &ErrorMessage(ostream &s) { return s ; }
  } ;
                                                                                
  register_BC<CheckSymmetryGridMotion> registerCheckSymmetryGridMotion ;

  
//---------------------------------------------------------------------------
// Rules to create constraints which are used to constrain the top and bottom 
// no-slip surfaces for communication with the structural solver
//
//
  // Constraint used to set top surface
  class TopCSDConstraint : public pointwise_rule {
    private:
 //     const_multiMap face2node ;
      store<bool> topCSDfaces ;
    public:

      // Define input and output.
      TopCSDConstraint() {
  //      name_store("face2node",face2node) ;
        name_store("topCSDfaces",topCSDfaces) ;
        output("topCSDfaces") ; // topCSDnodes -> nodes
        constraint("ref->top_BCoption") ; // topCSD_BCoption -> faces
      }

	  void calculate(Entity face) {topCSDfaces[face]=true ;}

      // Don't really have to do anything here since we are just setting
      // up a constraint.
      virtual void compute(const sequence &seq) {
      	do_loop(seq,this) ;
      	}

  } ;

  register_rule<TopCSDConstraint> registerTopCSDConstraint ;

  // Constraint used to set bottom surface
  class BottomCSDConstraint : public pointwise_rule {
    private:
 //     const_multiMap face2node ;
      store<bool> bottomCSDfaces ;
    public:

      // Define input and output.
      BottomCSDConstraint() {
  //      name_store("face2node",face2node) ;
        name_store("bottomCSDfaces",bottomCSDfaces) ;
        output("bottomCSDfaces") ;
        constraint("ref->bottom_BCoption") ;
      }

	  void calculate(Entity face) {bottomCSDfaces[face]=true ;}

      // Don't really have to do anything here since we are just setting
      // up a constraint.
      virtual void compute(const sequence &seq) {
      	do_loop(seq,this) ;
      	}

  } ;

  register_rule<BottomCSDConstraint> registerBottomCSDConstraint ;

// Constraint used to set bottom surface
  class TipCSDConstraint : public pointwise_rule {
    private:
 //     const_multiMap face2node ;
      store<bool> tipCSDfaces ;
    public:

      // Define input and output.
      TipCSDConstraint() {
  //      name_store("face2node",face2node) ;
        name_store("tipCSDfaces",tipCSDfaces) ;
        output("tipCSDfaces") ;
        constraint("ref->tip_BCoption") ;
      }

	  void calculate(Entity face) {tipCSDfaces[face]=true ;}

      // Don't really have to do anything here since we are just setting
      // up a constraint.
      virtual void compute(const sequence &seq) {
      	do_loop(seq,this) ;
      	}

  } ;

  register_rule<TipCSDConstraint> registerTipCSDConstraint ;
 

//-----------------------------------------------------------------------------
// Rules to create constraints which are used to constrain the displacement in
// a certain Cartesian direction for symmetry boundaries.

  // Constraint used to set top noslip
  class SXZeroConstraint : public pointwise_rule {
    private:
      const_multiMap face2node ;
      store<bool> sXZeroNodes ;
    public:

      // Define input and output.
      SXZeroConstraint() {
        name_store("face2node",face2node) ;
        name_store("sXZeroNodes",sXZeroNodes) ;
        output("face2node->sXZeroNodes") ;
        constraint("ref->sXZero_BCoption") ;
      }

	  void calculate(Entity node) {sXZeroNodes[node]=true ;}

      // Don't really have to do anything here since we are just setting
      // up a constraint.
      virtual void compute(const sequence &seq) {
      	do_loop(seq,this) ;
      	}

  } ;

  register_rule<SXZeroConstraint> registerSXZeroConstraint ;

  // Constraint used to set y-component of displacement to zero.
  class SYZeroConstraint : public pointwise_rule {
    private:
      const_multiMap face2node ;
      store<bool> sYZeroNodes ;
    public:

      // Define input and output.
      SYZeroConstraint() {
        name_store("face2node",face2node) ;
        name_store("sYZeroNodes",sYZeroNodes) ;
        output("face2node->sYZeroNodes") ;
        constraint("ref->sYZero_BCoption") ;
      }

	  void calculate(Entity node) {sYZeroNodes[node]=true ;}

      // Don't really have to do anything here since we are just setting
      // up a constraint.
      virtual void compute(const sequence &seq) {
      	do_loop(seq,this) ;
      }

  } ;

  register_rule<SYZeroConstraint> registerSYZeroConstraint ;

  // Constraint used to set z-component of displacement to zero.
  class SZZeroConstraint : public pointwise_rule {
    private:
      const_multiMap face2node ;
      store<bool> sZZeroNodes ;
    public:

      // Define input and output.
      SZZeroConstraint() {
        name_store("face2node",face2node) ;
        name_store("sZZeroNodes",sZZeroNodes) ;
        output("face2node->sZZeroNodes") ;
        constraint("ref->sZZero_BCoption") ;
      }

	  void calculate(Entity node) {sZZeroNodes[node]=true ;}

      // Don't really have to do anything here since we are just setting
      // up a constraint.
      virtual void compute(const sequence &seq) {
      	do_loop(seq,this) ;
      }

  } ;

  register_rule<SZZeroConstraint> registerSZZeroConstraint ;

//-----------------------------------------------------------------------------
// Rules to create a certain boundary constraints.

  // Create a constraint for node with specified displacement. This constraint
  // is used during final matrix assembly.
  class RigidBoundaryDisplacementConstraint : public pointwise_rule {
    private:
      const_multiMap face2node ;
      store<bool> rigidBoundaryDisplacement ;
    public:

      // Define input and output.
      RigidBoundaryDisplacementConstraint() {
        name_store("face2node",face2node) ;
        name_store("rigidBoundaryDisplacement",rigidBoundaryDisplacement) ;
        output("face2node->rigidBoundaryDisplacement") ;
        constraint("ref->s_BCoption") ; // face -> boundary surface specified with s=..
      }

			
//      // Assign flag for a single boundary face.
//      void calculate(Entity face) { rigidBoundaryDisplacement[node]=true ; }

      // Don't really have to do anything here since we are just setting
      // up a constraint.
      virtual void compute(const sequence &seq) {  
      	cout << "Rank: " << Loci::MPI_rank << ".. Setting the rigidBoundaryDisplacement rule" << endl ;
      	}

  } ;

  register_rule<RigidBoundaryDisplacementConstraint>
    registerRigidBoundaryDisplacementConstraint ;

 // All inlet faces have non-zero diffusion flux.
  class BoundaryFlexibleBoundaryDisplacementFaceTop : public pointwise_rule {
    private:
      store<bool> FlexibleBoundaryDisplacementFace ;
    public:

      // Define input and output.
      BoundaryFlexibleBoundaryDisplacementFaceTop() {
        name_store("FlexibleBoundaryDisplacementFace",FlexibleBoundaryDisplacementFace) ;
        output("FlexibleBoundaryDisplacementFace") ;
        constraint("ref->top_BCoption") ;
      }

      // Assign flag for a single boundary face.
   //   void calculate(Entity face) { FlexibleBoundaryDisplacementFace[face]=true ; }

      // Assign flag for a sequence of boundary faces.
      virtual void compute(const sequence &seq) { 
      	// do_loop(seq,this) ; 
      }
  } ;

  register_rule<BoundaryFlexibleBoundaryDisplacementFaceTop>
    registerBoundaryFlexibleBoundaryDisplacementFaceTop ;

  // All outlet faces have non-zero diffusion flux.
  class BoundaryFlexibleBoundaryDisplacementFaceBottom : public pointwise_rule {
    private:
      store<bool> FlexibleBoundaryDisplacementFace ;
    public:

      // Define input and output.
      BoundaryFlexibleBoundaryDisplacementFaceBottom() {
        name_store("FlexibleBoundaryDisplacementFace",FlexibleBoundaryDisplacementFace) ;
        output("FlexibleBoundaryDisplacementFace") ;
        constraint("ref->bottom_BCoption") ;
      }

      // Assign flag for a single boundary face.
  //    void calculate(Entity face) { FlexibleBoundaryDisplacementFace[face]=true ; }

      // Assign flag for a sequence of boundary faces.
      virtual void compute(const sequence &seq) { 
      	// do_loop(seq,this) ; 
      }
  } ;

  register_rule<BoundaryFlexibleBoundaryDisplacementFaceBottom>
    registerBoundaryFlexibleBoundaryDisplacementFaceBottom ;
    
  // All outlet faces have non-zero diffusion flux.
  class BoundaryFlexibleBoundaryDisplacementFaceTip : public pointwise_rule {
    private:
      store<bool> FlexibleBoundaryDisplacementFace ;
    public:

      // Define input and output.
      BoundaryFlexibleBoundaryDisplacementFaceTip() {
        name_store("FlexibleBoundaryDisplacementFace",FlexibleBoundaryDisplacementFace) ;
        output("FlexibleBoundaryDisplacementFace") ;
        constraint("ref->tip_BCoption") ;
      }

      // Assign flag for a single boundary face.
  //    void calculate(Entity face) { FlexibleBoundaryDisplacementFace[face]=true ; }

      // Assign flag for a sequence of boundary faces.
      virtual void compute(const sequence &seq) { 
      	// do_loop(seq,this) ; 
      }
  } ;

  register_rule<BoundaryFlexibleBoundaryDisplacementFaceTip>
    registerBoundaryFlexibleBoundaryDisplacementFaceTip ;  

  class FlexibleBoundaryDisplacementConstraint : public pointwise_rule {
    private:
      const_multiMap face2node ;
      store<bool> FlexibleBoundaryDisplacement ;
    public:

      // Define input and output.
      FlexibleBoundaryDisplacementConstraint() {
        name_store("face2node",face2node) ;
        name_store("FlexibleBoundaryDisplacement",FlexibleBoundaryDisplacement) ;
        output("face2node->FlexibleBoundaryDisplacement") ;
        constraint("FlexibleBoundaryDisplacementFace") ; // face -> boundary surface
      }

      // Don't really have to do anything here since we are just setting
      // up a constraint.
      virtual void compute(const sequence &seq) {  
      	cout << "Rank: " << Loci::MPI_rank << ".. Setting the FlexibleBoundaryDisplacement rule" << endl ;
      	}

  } ;

  register_rule<FlexibleBoundaryDisplacementConstraint>
    registerFlexibleBoundaryDisplacementConstraint ;
    
  class BoundaryDisplacementConstraintFlexibleFace : public pointwise_rule {
    private:
      store<bool> boundaryDisplacementFace ;
    public:

      // Define input and output.
      BoundaryDisplacementConstraintFlexibleFace() {
        name_store("boundaryDisplacementFace",boundaryDisplacementFace) ;
        output("boundaryDisplacementFace") ;
        constraint("FlexibleBoundaryDisplacementFace") ; // face -> boundary surface
      }

      // Don't really have to do anything here since we are just setting
      // up a constraint.
      virtual void compute(const sequence &seq) {  
      	cout << "Rank: " << Loci::MPI_rank << ".. Setting the boundaryDisplacement rule (flexible)" << endl ;
      	}

  } ;

  register_rule<BoundaryDisplacementConstraintFlexibleFace>
    registerBoundaryDisplacementConstraintFlexibleFace ;
    
  class BoundaryDisplacementConstraintRigidFace : public pointwise_rule {
    private:
      store<bool> boundaryDisplacementFace ;
    public:

      // Define input and output.
      BoundaryDisplacementConstraintRigidFace() {
        name_store("boundaryDisplacementFace",boundaryDisplacementFace) ;
        output("boundaryDisplacementFace") ;
        constraint("ref->s_BCoption") ; // face -> boundary surface
      }

      // Don't really have to do anything here since we are just setting
      // up a constraint.
      virtual void compute(const sequence &seq) {  
      	cout << "Rank: " << Loci::MPI_rank << ".. Setting the boundaryDisplacement rule (rigid)" << endl ;
      	}
  } ;

  register_rule<BoundaryDisplacementConstraintRigidFace>
    registerBoundaryDisplacementConstraintRigidFace ;

 class BoundaryDisplacementConstraint : public pointwise_rule {
    private:
      const_multiMap face2node ;
      store<bool> boundaryDisplacement ;
    public:

      // Define input and output.
      BoundaryDisplacementConstraint() {
        name_store("face2node",face2node) ;
        name_store("boundaryDisplacement",boundaryDisplacement) ;
        output("face2node->boundaryDisplacement") ;
        constraint("boundaryDisplacementFace") ; // face -> boundary surface
      }

      // Don't really have to do anything here since we are just setting
      // up a constraint.
      virtual void compute(const sequence &seq) {  
      	cout << "Rank: " << Loci::MPI_rank << ".. Setting the FlexibleBoundaryDisplacement rule" << endl ;
      	}

  } ;    

  register_rule<BoundaryDisplacementConstraint>
    registerBoundaryDisplacementConstraint ;
    
//-----------------------------------------------------------------------------
// Rules to create a constraint for all boundary edges belonging to non-symmetry
// boundaries. Note that this constraint also includes edges which are common
// between a symmetry and a non-symmetry boundary. This constraint will be
// determined based on displacement specification options rather than flow
// solver boundary condition types. All non-symmetry boundaries must have the
// displacement explicitly specified in some manner.

  // Rule to add nodes on boundaries with the "s=" option to the displacement
  // boundary condition constraint.
  class BoundaryNoSymmetryEdgeConstraint : public pointwise_rule {
    private:
      const_multiMap face2edge ;
      store<bool> boundaryNoSymmetryEdges ;
    public:

      // Define input and output.
      BoundaryNoSymmetryEdgeConstraint() {
        name_store("face2edge",face2edge) ;
        name_store("boundaryNoSymmetryEdges",boundaryNoSymmetryEdges) ;
        output("face2edge->boundaryNoSymmetryEdges") ;
        constraint("ref->s_BCoption") ;
      }

      // Don't really have to do anything here since we are just setting
      // up a constraint.
      virtual void compute(const sequence &seq) {}

  } ;

  register_rule<BoundaryNoSymmetryEdgeConstraint>
    registerBoundaryNoSymmetryEdgeConstraint ;

//-----------------------------------------------------------------------------
// Displacement boundary conditions. Note that we are assembling and solving
// the PDE for nodes attached to faces on "symmetry" boundaries.

  // Rule to get the displacement boundary condition value for the "s=" option.
  class BC_s_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<BoundaryDisplacement> s_BC ;
    public:

      // Define input and output.
      BC_s_compute() {
        name_store("BC_options",BC_options) ;
        name_store("s_BC",s_BC) ;
        input("BC_options") ;
        output("s_BC") ;
        constraint("s_BCoption") ;
      }

      // Get value for a single entity.
      void calculate(Entity e) {
        Loci::option_value_type optionValueType=BC_options[e].getOptionValueType("s") ;
        if(optionValueType==Loci::REAL){
          real temp ; BC_options[e].getOptionUnits("s","m",temp) ;
          s_BC[e].SetConstantValue(vect3d(temp,0.0,0.0)) ; s_BC[e].SetType(0) ;
          cout << "constantly (non)-moving boundary conditions [real] found" << endl ;
        }else if(optionValueType==Loci::LIST){
          vect3d temp ; get_vect3dOption(BC_options[e],"s","m",temp) ;
          s_BC[e].SetConstantValue(temp) ; s_BC[e].SetType(0) ;
          cout << "constantly (non)-moving boundary conditions [list] found" << endl ;
        }else if(optionValueType==Loci::FUNCTION){
          Loci::options_list::arg_list options ; string name ;
          BC_options[e].getOption("s",name,options) ;
          if(name=="rigidBody"){
          cout << "rigid-body displacement specification" << endl ;

            // Create an options list from the list of arguments so we can
            // use Ed's parsing stuff to retreive values.
            string validOptions="sConst:tDir:tMag:tFreq:tPhi:rCenter:rAxis:" ;
            validOptions+="rAlphaBar:rMag:rFreq:rPhi:func" ;
            Loci::options_list rigidBodyOptions(validOptions) ;
            rigidBodyOptions.Input(options) ;

            // Now extract the options and set the value.
            char *optionName[]={"sConst","tDir","tMag","tFreq","tPhi",
              "rCenter","rAxis","rAlphaBar","rMag","rFreq","rPhi","func"} ;
            RigidBodyDisplacement rigidBodyDisplacement ;
            for(unsigned int i=0;i<12;++i){
              if(rigidBodyOptions.optionExists(optionName[i])){

                // Scalars.
                if(i==0 || i==2 || i==3 || i==4 || i==7 || i==8 || i==9 || i==10){
                  real value ; rigidBodyOptions.getOption(optionName[i],value) ;
                  if(i==0){ rigidBodyDisplacement.sConst=value ;
                  }else if(i==2){ rigidBodyDisplacement.tMag=value ;
                  }else if(i==3){ rigidBodyDisplacement.tFreq=value ;
                  }else if(i==4){ rigidBodyDisplacement.tPhi=value ;
                  }else if(i==7){ rigidBodyDisplacement.rAlphaBar=value ;
                  }else if(i==8){ rigidBodyDisplacement.rMag=value ;
                  }else if(i==9){ rigidBodyDisplacement.rFreq=value ;
                  }else if(i==10){ rigidBodyDisplacement.rPhi=value ;
                  }
                }
                // Vectors.
                if(i==1 || i==5 || i==6){
                  vect3d value=get_vect3d(rigidBodyOptions,optionName[i],"m") ;
                  if(i==1){ rigidBodyDisplacement.tDir=value ;
                  }else if(i==5){ rigidBodyDisplacement.rCenter=value ;
                  }else if(i==6){ rigidBodyDisplacement.rAxis=value ;
                  }
                }
                // Integers
		if(i==11) {
		  real value; rigidBodyOptions.getOption(optionName[i],value) ;
		  if (i==11) { rigidBodyDisplacement.func=(value) ;
		  }
		}
              }
            }
            s_BC[e].SetRigidBodyValue(rigidBodyDisplacement) ;
            s_BC[e].SetType(1) ; // rigidBodyMotion
//          }else if(name=="flexible"){
//          	cout << "Flexible boundary conditions found" << endl ;
//          	s_BC[e].SetType(2) ; 
	  }else{
            cerr << "Bad displacement boundary condition in .vars file.: " << name  << endl ; Loci::Abort() ;
          }
        }else{
          cerr << "Bad displacement option type in .vars file." << endl ;
          Loci::Abort() ;
        }
      }

      // Loop over entities.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_s_compute> register_BC_s_compute ;

  // Rule for boundary faces with specified displacement. Assigns displacement
  // value to all nodes connected to boundary faces that have the property s_BC.
  // Note that we have no way of determining in what order the faces will be
  // processed, so the user is responsible for ensuring that all boundary
  // patches that have a common line of nodes have the same displacement
  // boundary condition.
//  class BoundaryDisplacementSpecification : public pointwise_rule {
  class BoundaryDisplacementSpecification : public pointwise_rule {	
    private:
      const_Map ref ;
      const_multiMap face2node ;
      const_param<real> solutionTime ;
      const_param<real> timeStep ;
      const_store<BoundaryDisplacement> s_BC ;
      const_store<vect3d> pos ;
      store<vect3d> node_s_b ;
//      store<bool> boundaryNodes ;
    public:

      // Define input and output.
      BoundaryDisplacementSpecification() {
        name_store("ref",ref) ;
        name_store("face2node",face2node) ;
        name_store("stime{n}",solutionTime) ;
        name_store("timeStep{n}",timeStep) ;
        name_store("s_BC{n}",s_BC) ;
        name_store("pos",pos) ;
        name_store("node_s_b{n}",node_s_b) ;
 //       name_store("boundaryNodes",boundaryNodes) ;
        input("stime{n},timeStep{n},ref->s_BC{n},face2node->pos") ;
        output("face2node->node_s_b{n}") ;
//        output("face2node->boundaryNodes") ;
	constraint("gridMotionTimeDependent{n}") ;
	constraint("ref->s_BCoption") ;
      }

      // Set boundary displacement for nodes attached to this face. Note that
      // this rule is actually setting each node's displacement value multiple
      // times (the same value each time). We should fix this, but for now
      // we'll just leave it.
      void calculate(Entity face) {
		for(const int *ni=face2node.begin(face);ni!=face2node.end(face);++ni){
		  node_s_b[*ni]=s_BC[ref[face]].Value(*solutionTime+*timeStep,pos[*ni]) ;
//cout << "n,s: " << *ni << " " << node_s_b[*ni] << endl ;
//			boundaryNodes[*ni]=true ;
		}
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;

  register_rule<BoundaryDisplacementSpecification>
    registerBoundaryDisplacementSpecification ;
    
 // class BoundaryDisplacementSpecificationSolutionDependentApply : public apply_rule<store<vect3d>,Loci::NullOp<vect3d> > {
  class BoundaryDisplacementSpecificationSolutionDependentApply : public pointwise_rule {
    private:
      const_Map ref ;
      const_multiMap face2node ;
      const_param<real> solutionTime ;
      const_param<real> timeStep ;
      const_param<int> CSDstartingTimeStep ;
      const_store<BoundaryDisplacement> s_BC ;
      const_store<vect3d> pos ;
      store<vect3d> node_s_b_bc ;
//      store<bool> boundaryNodes ;
      real laggedTime ;
    public:

      // Define input and output.
      BoundaryDisplacementSpecificationSolutionDependentApply() {
        name_store("ref",ref) ;
        name_store("face2node",face2node) ;
        name_store("stime{n}",solutionTime) ;
        name_store("timeStep{n}",timeStep) ;
        name_store("s_BC{n}",s_BC) ;
        name_store("pos",pos) ;
        name_store("node_s_b_bc{n}",node_s_b_bc) ;
	name_store("CSDstartingTimeStep", CSDstartingTimeStep) ;
	input("CSDstartingTimeStep") ;
 //       name_store("boundaryNodes",boundaryNodes) ;
        input("stime{n},timeStep{n},ref->s_BC{n},face2node->pos") ;
        output("face2node->node_s_b_bc{n}") ;
	constraint("gridMotionSolutionDependent{n}") ;
	//constraint("rigidBoundaryDisplacement{n}") ;				
	constraint("ref->s_BCoption") ; // top and bottom nodes done in FSI_CSD2CFD.loci, here just 0.0
//				disable_threading() ;
//        output("face2node->boundaryNodes") ;
      }

      // Set boundary displacement for nodes attached to this face. Note that
      // this rule is actually setting each node's displacement value multiple
      // times (the same value each time). We should fix this, but for now
      // we'll just leave it.
      void calculate(Entity face) {
	//real laggedTime = (*solutionTime) - (*timeStep) * (*CSDstartingTimeStep) ;
	for(const int *ni=face2node.begin(face);ni!=face2node.end(face);++ni){
	  if (laggedTime >= 0.0) {
//	    node_s_b_bc[*ni]=s_BC[ref[face]].Value(*solutionTime+*timeStep,pos[*ni]) ;
	    node_s_b_bc[*ni]=s_BC[ref[face]].Value(laggedTime,pos[*ni]) ;
	  } else {
	    node_s_b_bc[*ni]=vect3d(0.0,0.0,0.0);
	  }
//cout << "n,s: " << *ni << " " << node_s_b[*ni] << endl ;
//			boundaryNodes[*ni]=true ;
	}
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { 
	laggedTime = (*solutionTime) - (*timeStep) * ((*CSDstartingTimeStep)+1) + (*timeStep);
	if (Loci::MPI_rank==0) cout << "SolutionTime, TimeStep, CSDstartingTimeStep, laggedTime = " << (*solutionTime) << ", " << (*timeStep) << ", " << (*CSDstartingTimeStep) << ", " << (laggedTime) << endl ;
	do_loop(seq,this) ; 
      }
                                                                                
  } ;

  register_rule<BoundaryDisplacementSpecificationSolutionDependentApply>
    registerBoundaryDisplacementSpecificationSolutionDependentApply ;

//-----------------------------------------------------------------------------
// Velocity boundary conditions.

  // Rule for assigning velocity on no-slip boundary faces.
//  class BoundaryVelocityNoSlipGridMotion : public pointwise_rule {
//    private:
//      const_store<vect3d> face_v ;
//      store<vect3d> v_f ;
//    public:
//                                                                                
//      // Define input and output.
//      BoundaryVelocityNoSlipGridMotion() {
//        name_store("face_v",face_v) ;
//        name_store("priority::gridMotion::v_f",v_f) ;
//        input("face_v") ;
//        output("priority::gridMotion::v_f") ;
//        constraint("noslip_BC") ;
//      }
//                                                                                
//      // Calculate face velocity for a single face.
//      void calculate(Entity face) { v_f[face]=face_v[face] ; 
//      	//if (Loci::MPI_rank==0) cout << "face_v=v_f=" << v_f[face] << endl ;
//      	}
//                                                                                
//      // Calculate face velocity for all faces in sequence.
//      virtual void compute(const sequence &seq) { 
//      	do_loop(seq,this) ; 
//      	if (Loci::MPI_rank==0) cout << "gridMotion::v_f " << endl ;
//      	}
//  } ;
//                                                                                
//  register_rule<BoundaryVelocityNoSlipGridMotion>
//    registerBoundaryVelocityNoSlipGridMotion ;

  // Rule for assigning velocity on no-slip boundary faces.
  class BoundaryVelocityboundaryDisplacementFaceGridMotion : public pointwise_rule {
    private:
      const_store<vect3d> face_v ;
      store<vect3d> v_f ;
    public:
                                                                                
      // Define input and output.
      BoundaryVelocityboundaryDisplacementFaceGridMotion() {
        name_store("face_v",face_v) ;
        name_store("gridMotion::v_f",v_f) ;
        input("face_v") ;
        output("gridMotion::v_f") ;
     //   constraint("boundaryDisplacementFace") ;
        constraint("faces") ;
        constraint("noslip_BC") ;
      }
                                                                                
      // Calculate face velocity for a single face.
      void calculate(Entity face) { v_f[face]=face_v[face] ; 
      	//if (Loci::MPI_rank==0) 
	  cout << "face_v=v_f=" << v_f[face] << endl ;
      	}
                                                                                
      // Calculate face velocity for all faces in sequence.
      virtual void compute(const sequence &seq) { 
      	do_loop(seq,this) ; 
      	if (Loci::MPI_rank==0) cout << "boundaryDisplacementFace->gridMotion::v_f " << endl ;
      	}
  } ;
                                                                                
  register_rule<BoundaryVelocityboundaryDisplacementFaceGridMotion>
    registerBoundaryVelocityboundaryDisplacementFaceGridMotion ;

  // Rule for for assigning velocity on slip boundary faces. We have replaced
  // the old slip bc that we never used for stability reasons with this new
  // one since we need to distinguish between slip and symmetry for moving
  // boundary problems (i.e. the normal component of the fluid velocity on a
  // slip surface is equal to the face velocity, but zero for a symmetry
  // surface).
  class BoundaryVelocitySlipGridMotion : public pointwise_rule {
    private:
      const_store<Area> area ;
      const_store<vect3d> face_v ;
      const_Map ci ;
      const_store<vect3d> v ;
      store<vect3d> v_f ;
    public:
                                                                                
      // Define input and output.
      BoundaryVelocitySlipGridMotion() {
        name_store("area",area) ;
        name_store("face_v",face_v) ;
        name_store("ci",ci) ;
        name_store("v",v) ;
        name_store("priority::v_f",v_f) ;
        input("area,face_v,ci->v") ;
        output("priority::v_f") ;
        constraint("slip_BC") ;
      }
                                                                                
      // Calculate face velocity for a single face.
      void calculate(Entity face) {
        const vect3d vCell=v[ci[face]] ;
        v_f[face]=vCell-dot(vCell,area[face].n)*area[face].n+face_v[face] ;
      }
                                                                                
      // Calculate face velocity for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<BoundaryVelocitySlipGridMotion>
    registerBoundaryVelocitySlipGridMotion ;


}