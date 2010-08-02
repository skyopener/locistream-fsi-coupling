//-----------------------------------------------------------------------------
// Description: This file contains some of the basic rules common to all
//   grid movement schemes.
//-----------------------------------------------------------------------------

// Standard library includes.
#include <vector>
using std::vector ;

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// PETSC includes
//#include "petsc.h"
//#include <petscerror.h>
//#include <petscksp.h>

// StreamUns includes.
#include "const.h"
#include "move.h"
#include "residual.h"
#include "sciTypes.h"
#include "varsFileInputs.h"

// Include PETSC wrappers
#include "FSI_petsc.h"

// Forward declarations using C linkage to avoid mangling the names of
// the Fortran functions that will be called.
extern "C" {
 void getNumberOfNodesFromCSD_(int*) ;
 void getNodalCoordinatesFromCSD_(double*) ;
 void getNodalDisplacementsFromCSD_(double*) ;
}

namespace streamUns {

  //---------------------------------------------------------------------------
  // Implementation for FlexibleBodyDisplacement.

  // Assemble matrix Q which is the matrix with the nodal coordinates
  void FlexibleBodyDisplacement::AssembleQ() {

    int N;
	getNumberOfNodesFromCSD_(&N) ;
	setNb(N) ;
	matrixQ.resize(Nb, d) ;
	matrixB.resize(Nb, d) ;
	getNodalDisplacementsFromCSD_(&matrixB[0]) ;
	getNodalCoordinatesFromCSD_(&matrixQ[0]) ;
  }

  // Assemble the system matrix A, RHS b, solve, for the Weight coefficients
  void getFSIrbfWeight() {
	
	PetscScalar valuex, valuey, valuez ;
	real distance=0. ;
	int columnIndex ;

	// Assemble the RHS vectors
	PETSCMultiVector FSI_RHS_vector, FSI_phi ; // declare RHS vector
	FSI_RHS_vector.Create(Nbd1, Nbd1) ; // create RHS vector
	for(int i=0;i<Nb;++i) {
		valuex=matrixB(i,0);
		valuey=matrixB(i,1);
		valuez=matrixB(i,2);
		FSI_RHS_vector.SetValue(&i, &valuex, &valuey, &valuez) ;
	}


	// Determine the number of non-zero row-values
	vector<int> numberOfNonZeros(Nbd1,0)
	for(int i=0;i<Nb;++i) {
		for(int j=0;j<Nb;++j) {
			if (rbfNr > 0 && rbfNr < 9) {
				distance = pow(matrixQ(i,0) - matrixQ(j,0), 2.) ;
				distance += pow(matrixQ(i,1) - matrixQ(j,1), 2.) ;
				distance += pow(matrixQ(i,2) - matrixQ(j,2), 2.) ;
				distance /= r * r ;
				if (distance <= 1.) {
					++numberOfNonZeros(i);
				}
			} else {
				++numberOfNonZeros(i);
			}
		}
	}
	numberOfNonZeros(Nbd1-4)=Nb;
	numberOfNonZeros(Nbd1-3)=Nb;
	numberOfNonZeros(Nbd1-2)=Nb;
	numberOfNonZeros(Nbd1-1)=Nb;


	// Setup the system matrix A
    PETSCMatrix FSI_A_matrix ;
	int temp=0;
	FSI_A_matrix.Create(Nbd1, Nbd1, &numberOfNonZeros[0], &temp) ;
	numberOfNonZeros.erase(numberOfNonZeros.begin(), numberOfNonZeros.end() ); // delete all elements of numberOfNonZeros

	// Fill in the matrix A
	for(int i=0;i<Nb;++i) {
		for(int j=0;j<Nb;++j) {
			distance = pow(matrixQ(i,0) - matrixQ(j,0), 2.) ;
			distance += pow(matrixQ(i,1) - matrixQ(j,1), 2.) ;
			distance += pow(matrixQ(i,2) - matrixQ(j,2), 2.) ;
			distance = sqrt(distance) ;
			if (!(rbfNr > 0 && rbfNr < 9) && (distance/r > 1.)) {
				valuex = radialBasisFunction(distance, r, a, rbfNr) ;
				columnIndex=j; //columnValue=valuex ;
				FSI_A_matrix.SetRowValues(i,1,&columnIndex,&valuex);
			}
		}
		for(int j=0; j<4; ++j) {
			columnIndex=Nb+j; valuex=((j==0)?1.:matrixQ(i,j-1)); 
			FSI_A_matrix.SetRowValues(i,1,&columnIndex,&valuex) ;
			columnIndex=i;
			FSI_A_matrix.SetRowValues(Nb+j,1,&columnIndex,&valuex) ;
		}
	}
	FSI_A_matrix.GetInfo() ;
	FSI_A_matrix.AssemblyBegin() ; FSI_A_matrix.AssemblyEnd() ;

	// Setup the KSP solver
	PETSCKsp FSI_ksp ;
	FSI_ksp.Create() ;
	FSI_ksp.SetTolerances(linearSolverTolerance,1.0e-30,linearSolverMaxIterations) ;
	PC pc ; FSI_ksp.GetPC(&pc) ; PCSetType(pc,PCJACOBI) ;
	FSI_ksp.SetInitialGuessNonzero() ;
	FSI_ksp.SetFromOptions() ;

	// Solve the system
	FSI_ksp.SetOperators(FSI_A_matrix) ; FSI_ksp.






			}
		}
	}




	

	

  // Returns the serialized buffer size.
  int FlexibleBodyDisplacement::BufferSize() const {
    int bufferSize=5 ; return bufferSize ;
  }

  // Packs the data into a buffer.
  void FlexibleBodyDisplacement::PackBuffer(real *buffer,int size) {
    int i=0 ;
    buffer[i++]=r ; buffer[i++]=a ; buffer[i++]=rbfNr ;
    buffer[i++]=linearSolverTolerance ; buffer[i++]=linearSolverMaxIterations ;
  }

  // Unpacks the data from a buffer.
  void FlexibleBodyDisplacement::UnpackBuffer(real *buffer,int size) {
    int i=0 ;
    r=buffer[i++] ; a=buffer[i++] ; rbfNr=buffer[i++] ;
    linearSolverTolerance=buffer[i++] ; linearSolverMaxIterations=buffer[i++] ;
  }

  vect3d FlexibleBodyDisplacement::Value(vect3d position) const {

    // Compute the rotation origin, which is the point on the axis of
    // rotation that when connected to this node is perpendicular to the
    // axis of rotation.
    vect3d r=position-rCenter ; real s0=dot(r,rAxis)/dot(rAxis,rAxis) ;
    vect3d rOrigin=rCenter+s0*rAxis ;

    // Compute normalized position vector from rotation center. Compute
    // vector perpendicular to postion vector and rotation axis.
    vect3d rr=position-rOrigin ; real R=norm(rr) ; vect3d rHat=(1.0/R)*rr ;
    vect3d tHat=cross(rHat,rAxis) ;

    // Compute displacement.
    real theta=rAlphaBar+rMag*sin(rFreq*time+rPhi) ;
    return (sConst+tMag*sin(tFreq*time+tPhi))*tDir+R*((cos(theta)-1.0)*rHat-
      sin(theta)*tHat) ;
  }

  //---------------------------------------------------------------------------
  // Implementation for RigidBodyDisplacement.

  // Returns the serialized buffer size.
  int RigidBodyDisplacement::BufferSize() const {
    int bufferSize=17 ; return bufferSize ;
  }

  // Packs the data into a buffer.
  void RigidBodyDisplacement::PackBuffer(real *buffer,int size) {
    int i=0 ;
    buffer[i++]=sConst ; buffer[i++]=tMag ; buffer[i++]=tFreq ;
    buffer[i++]=tPhi ; buffer[i++]=rAlphaBar ; buffer[i++]=rMag ;
    buffer[i++]=rFreq ; buffer[i++]=rPhi ;
    buffer[i++]=tDir.x ; buffer[i++]=tDir.y ; buffer[i++]=tDir.z ;
    buffer[i++]=rCenter.x ; buffer[i++]=rCenter.y ; buffer[i++]=rCenter.z ;
    buffer[i++]=rAxis.x ; buffer[i++]=rAxis.y ; buffer[i++]=rAxis.z ;
  }

  // Unpacks the data from a buffer.
  void RigidBodyDisplacement::UnpackBuffer(real *buffer,int size) {
    int i=0 ;
    sConst=buffer[i++] ; tMag=buffer[i++] ; tFreq=buffer[i++] ;
    tPhi=buffer[i++] ; rAlphaBar=buffer[i++] ; rMag=buffer[i++] ;
    rFreq=buffer[i++] ; rPhi=buffer[i++] ;
    tDir.x=buffer[i++] ; tDir.y=buffer[i++] ; tDir.z=buffer[i++] ;
    rCenter.x=buffer[i++] ; rCenter.y=buffer[i++] ; rCenter.z=buffer[i++] ;
    rAxis.x=buffer[i++] ; rAxis.y=buffer[i++] ; rAxis.z=buffer[i++] ;
  }

  vect3d RigidBodyDisplacement::Value(real time,vect3d position) const {

    // Compute the rotation origin, which is the point on the axis of
    // rotation that when connected to this node is perpendicular to the
    // axis of rotation.
    vect3d r=position-rCenter ; real s0=dot(r,rAxis)/dot(rAxis,rAxis) ;
    vect3d rOrigin=rCenter+s0*rAxis ;

    // Compute normalized position vector from rotation center. Compute
    // vector perpendicular to postion vector and rotation axis.
    vect3d rr=position-rOrigin ; real R=norm(rr) ; vect3d rHat=(1.0/R)*rr ;
    vect3d tHat=cross(rHat,rAxis) ;

    // Compute displacement.
    real theta=rAlphaBar+rMag*sin(rFreq*time+rPhi) ;
    return (sConst+tMag*sin(tFreq*time+tPhi))*tDir+R*((cos(theta)-1.0)*rHat-
      sin(theta)*tHat) ;
  }

  //---------------------------------------------------------------------------
  // Implementation for BoundaryDisplacement.

  // Returns the serialized buffer size.
  int BoundaryDisplacement::BufferSize() const {
    int bufferSize=4+rigidBody.BufferSize() ; return bufferSize ;
  }

  // Overridden virtual method called by fact database when a fact of type
  // BoundaryDisplacment is read. Will never be used.
  istream& BoundaryDisplacement::Input(istream &in) {
    return in ;
  }

  // Packs the data into a buffer.
  void BoundaryDisplacement::PackBuffer(real *buffer,int size) {
    int i=0 ;
    buffer[i++]=real(type) ;
    buffer[i++]=constant.x ; buffer[i++]=constant.y ; buffer[i++]=constant.z ;
    rigidBody.PackBuffer(buffer+4,size) ;
  }

  // Prints the values. Will never be used.
  ostream& BoundaryDisplacement::Print(ostream &out) const {
    return out ;
  }

  // Unpacks the data from a buffer.
  void BoundaryDisplacement::UnpackBuffer(real *buffer,int size) {
    int i=0 ;
    type=int(buffer[i++]) ;
    constant.x=buffer[i++] ; constant.y=buffer[i++] ; constant.z=buffer[i++] ;
    rigidBody.UnpackBuffer(buffer+4,size) ;
  }

  //---------------------------------------------------------------------------
  // Setup Rules.

  // The following unit/apply sequence of rules identifies if the mesh
  // deformation is time-dependent or solution dependent. The value is 1
  // for time-dependent and 0 for solution-dependent. If there are any
  // boundary conditions that are solution-dependent, then the logial AND
  // operation for the apply rule will cause this value to be false.
  class MeshDeformationTimeDependentUnit : public unit_rule {
    private:
      param<string> gridMoverRbf ;
      param<int> meshDeformationTimeDependent ;
    public:
                                                                                
      // Define input and output.
      MeshDeformationTimeDependentUnit() {
        name_store("gridMoverRbf",gridMoverRbf) ;
        name_store("meshDeformationTimeDependent",meshDeformationTimeDependent) ;
        input("gridMoverRbf") ;
        output("meshDeformationTimeDependent") ;
      }
                                                                                
      // Set the default value.
      virtual void compute(const sequence &seq) {
        *meshDeformationTimeDependent=1 ; // RBF time dependent
      }
  } ;
                                                                                
  register_rule<MeshDeformationTimeDependentUnit>
    registerMeshDeformationTimeDependentUnit ;

  // The apply rule that checks for solution-dependent deforming boundaries.
  class MeshDeformationTimeDependentApply : public apply_rule<param<int>,Loci::Product<int> > {
    private:
      const_store<BoundaryDisplacement> s_BC ;
      param<int> meshDeformationTimeDependent ;
    public:
                                                                                
      // Define input and output.
      MeshDeformationTimeDependentApply() {
        name_store("s_BC",s_BC) ;
        name_store("meshDeformationTimeDependent",meshDeformationTimeDependent) ;
        input("s_BC") ;
        output("meshDeformationTimeDependent") ;
        constraint("s_BCoption") ;
      }

      // Look at a single boundary.
      void calculate(Entity e) {
        int flag=(s_BC[e].IsSolutionDependent())? 0:1 ;
        join(*meshDeformationTimeDependent,flag) ;
      }
                                                                                
      // Loop over the boundaries.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<MeshDeformationTimeDependentApply>
    registerMeshDeformationTimeDependentApply ;

  // Set the grid motion type for deforming meshes.
  class DeformingMeshGridMotionType : public constraint_rule {
    private:
      const_param<int> meshDeformationTimeDependent ;
      Constraint gridMotionTimeDependent,gridMotionSolutionDependent ;
    public:

      // Define input and output.
      DeformingMeshGridMotionType() {
        name_store("meshDeformationTimeDependent",meshDeformationTimeDependent);
        name_store("gridMotionTimeDependent",gridMotionTimeDependent) ;
        name_store("gridMotionSolutionDependent",gridMotionSolutionDependent) ;
        input("meshDeformationTimeDependent") ;
        output("gridMotionTimeDependent,gridMotionSolutionDependent") ;
      }

      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        gridMotionTimeDependent=(*meshDeformationTimeDependent==1)?
          ~EMPTY:EMPTY ;
        gridMotionSolutionDependent=(*meshDeformationTimeDependent==0)?
          ~EMPTY:EMPTY ;
      }
  } ;

  register_rule<DeformingMeshGridMotionType> registerDeformingMeshGridMotionType ;

  // Rule to declare "r" as valid .vars file variable. This variable is used
  // to normalize the distances in the rbfs
  class DefaultRBFr : public default_rule {
    private:
      param<real> rbfR ;
    public:
                                                                                
      // Define input and output.
      DefaultRBFr() {
        name_store("rbfR",rbfR) ;
        output("rbfR") ;
      }
                                                                                
      // Set the default value.
      virtual void compute(const sequence &seq) { *rbfR=1.0 ; }
  } ;
                                                                                
  register_rule<DefaultRBFr> registerDefaultRBFr ;

  // Rule to declare "a" as valid .vars file variable. This variable is used
  // as a parameters in the global RBFs
  class DefaultRBFa : public default_rule {
    private:
      param<real> rbfA ;
    public:
                                                                                
      // Define input and output.
      DefaultRBFa() {
        name_store("rbfA",rbfA) ;
        output("rbfA") ;
      }
                                                                                
      // Set the default value.
      virtual void compute(const sequence &seq) { *rbfA=1.0 ; }
  } ;
                                                                                
  register_rule<DefaultRBFa> registerDefaultRBFa ;

  // Rule to declare "rbfLinearSolver" as valid .vars file variable.
  class DefaultRbfLinearSolver : public default_rule {
    private:
      param<string> rbfLinearSolver ;
    public:
                                                                                
      // Define input and output.
      DefaultRbfLinearSolver() {
        name_store("rbfLinearSolver",rbfLinearSolver) ;
        output("rbfLinearSolver") ;
      }
                                                                                
      // Set the default value.
      virtual void compute(const sequence &seq) {
        *rbfLinearSolver="PETSC" ;
			if (Loci::MPI_rank==0) cout << "RBFsolver = " << *rbfLinearSolver << endl ;
      }
  } ;
                                                                                
  register_rule<DefaultRbfLinearSolver>
    registerDefaultRbfLinearSolver ;

  // Rule to declare "RBFMaxLinearSolverIterations" as valid .vars file
  // variable. This variable is passed to the linear solver.
  class DefaultRbfMaxLinearSolverIterations : public default_rule {
    private:
      param<int> rbfMaxLinearSolverIterations ;
    public:
                                                                                
      // Define input and output.
      DefaultRbfMaxLinearSolverIterations() {
        name_store("rbfMaxLinearSolverIterations",
          rbfMaxLinearSolverIterations) ;
        output("rbfMaxLinearSolverIterations") ;
      }
                                                                                
      // Set the default value.
      virtual void compute(const sequence &seq) {
        *rbfMaxLinearSolverIterations=10 ;
      }
  } ;
                                                                                
  register_rule<DefaultRbfMaxLinearSolverIterations>
    registerDefaultRbfMaxLinearSolverIterations ;

  // Rule to declare "RBFnumber" as valid .vars file variable. This variable
  // represents the RBF choice
  class DefaultRbfNumber : public default_rule {
    private:
      param<int> rbfNumber ;
    public:

      // Define input and output.
      DefaultRbfNumber() {
        name_store("rbfNumber",rbfNumber) ;
        output("rbfNumber") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) { *rbfNumber=1 ; }
  } ;

  register_rule<DefaultRbfNumber> registerDefaultRbfNumber ;

  // Rule to declare "rbfTolerance" as valid .vars file variable.
  // This variable is passed to the linear solver.
  class DefaultRbfTolerance : public default_rule {
    private:
      param<real> rbfTolerance ;
    public:
                                                                                
      // Define input and output.
      DefaultRbfTolerance() {
        name_store("rbfTolerance",rbfTolerance) ;
        output("rbfTolerance") ;
      }
                                                                                
      // Set the default value.
      virtual void compute(const sequence &seq) {
        *rbfTolerance=1.0e-03 ;
      }
  } ;
                                                                                
  register_rule<DefaultRbfTolerance> registerDefaultRbfTolerance ;
                                                                                
  // Creates the linear solver constraint.
  class rbfLinearSolverConstraint : public constraint_rule {
    private:
      const_param<string> rbfLinearSolver ;
      Constraint s_JacobiLinearSolver,s_PETSCLinearSolver,s_HYPRELinearSolver ;
    public:
                                                                                
      // Define input and output.
      rbfLinearSolverConstraint() {
        name_store("rbfLinearSolver",rbfLinearSolver) ;
        name_store("sStar_JacobiLinearSolver",s_JacobiLinearSolver) ;
        name_store("sStar_PETSCLinearSolver",s_PETSCLinearSolver) ;
        name_store("sStar_HYPRELinearSolver",s_HYPRELinearSolver) ;
        input("rbfLinearSolver") ;
        output("sStar_JacobiLinearSolver,sStar_PETSCLinearSolver") ;
        output("sStar_HYPRELinearSolver") ;
      }
                                                                                
      // Set up the constraint.
      virtual void compute(const sequence& seq) {
			if (Loci::MPI_rank==0) cout << "constraint: RBFsolver = " << *rbfLinearSolver << endl ;
        if(*rbfLinearSolver=="Jacobi"){
          s_JacobiLinearSolver=~EMPTY ; s_PETSCLinearSolver=EMPTY ;
          s_HYPRELinearSolver=EMPTY ;
        }else if(*rbfLinearSolver=="PETSC"){
			if (Loci::MPI_rank==0) cout << "RBFsolver == PETSC " << endl ;
          s_JacobiLinearSolver=EMPTY ; s_PETSCLinearSolver=~EMPTY ;
          s_HYPRELinearSolver=EMPTY ;
        }else if(*rbfLinearSolver=="HYPRE"){
          s_JacobiLinearSolver=EMPTY ; s_PETSCLinearSolver=EMPTY ;
          s_HYPRELinearSolver=~EMPTY ;
        }else{
          cerr << "ERROR: Bad linear solver type for grid mover." << endl ;
          Loci::Abort() ;
        }
      }
  } ;
                                                                                
  register_rule<rbfLinearSolverConstraint>
    registerrbfLinearSolverConstraint ;

  // Creates the maximum iterations variable used by the linear solvers.
  class LinearSolverMaxIterations : public singleton_rule {
    private:
      const_param<int> rbfMaxLinearSolverIterations ;
      param<int> s_maxLinearSolverIterations ;
    public:
                                                                                
      // Define input and output.
      LinearSolverMaxIterations() {
        name_store("rbfMaxLinearSolverIterations",
          rbfMaxLinearSolverIterations) ;
        name_store("sStar_maxLinearSolverIterations",
          s_maxLinearSolverIterations) ;
        input("rbfMaxLinearSolverIterations") ;
        output("sStar_maxLinearSolverIterations") ;
      }
                                                                                
      // Set the value.
      virtual void compute(const sequence &seq) {
        *s_maxLinearSolverIterations=*rbfMaxLinearSolverIterations ;
      }
  } ;
                                                                                
  register_rule<LinearSolverMaxIterations> registerLinearSolverMaxIterations ;

  // Creates the tolerance variable used by the linear solvers.
  class LinearSolverTolerance : public singleton_rule {
    private:
      const_param<real> rbfTolerance ;
      param<real> s_linearSolverTolerance ;
    public:
                                                                                
      // Define input and output.
      LinearSolverTolerance() {
        name_store("rbfTolerance",rbfTolerance) ;
        name_store("sStar_linearSolverTolerance",s_linearSolverTolerance) ;
        input("rbfTolerance") ;
        output("sStar_linearSolverTolerance") ;
      }
                                                                                
      // Set the value.
      virtual void compute(const sequence &seq) {
        *s_linearSolverTolerance=*rbfTolerance ;
      }
  } ;
                                                                                
  register_rule<LinearSolverTolerance> registerLinearSolverTolerance ;

  // Rule to declare "gridMover" as valid .vars file variable.
  class OptionalGridMover : public optional_rule {
    private:
      param<string> gridMoverRbf ;
    public:
                                                                                
      // Define input and output.
      OptionalGridMover() {
        name_store("gridMoverRbf",gridMoverRbf) ;
        output("gridMoverRbf") ;
      }
                                                                                
      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;
                                                                                
  register_rule<OptionalGridMover> registerOptionalGridMover ;

//-----------------------------------------------------------------------------
// Priority rule to compute initial face areas when there is grid movement.

  class InitialFaceCenterGridMotion : public pointwise_rule {
    private:
      const_multiMap face2node ;
      const_store<vect3d> pos ;
      store<vect3d> facecenter ;
    public :

      // Define input and output.
     InitialFaceCenterGridMotion() {
        name_store("face2node",face2node) ;
        name_store("pos_ic",pos) ;
        name_store("facecenter_ic",facecenter) ;
        input("face2node->pos_ic") ;
        output("facecenter_ic") ;
     }

     // Compute center for a single face.
     void calculate(Entity face) {
       vect3d nodesum(0.0,0.0,0.0) ; real lensum=0.0 ;
       for(const int *id=face2node.begin(face);id+1!=face2node.end(face);++id) {
         vect3d pos0=pos[*id],pos1=pos[*(id+1)] ;
         vect3d edge_loc=0.5*(pos0+pos1),edge_vec=pos0-pos1 ;
         real len=sqrt(dot(edge_vec,edge_vec)) ;
         nodesum+=len*edge_loc ; lensum += len ;
       }
       const int *id=face2node.begin(face),*idend=face2node.end(face)-1 ;
       vect3d pos0=pos[*id],pos1=pos[*idend] ;
       vect3d edge_loc=0.5*(pos0+pos1),edge_vec=pos0-pos1 ;
       real len=sqrt(dot(edge_vec,edge_vec)) ;
       nodesum+=len*edge_loc ; lensum += len ;
       facecenter[face]=nodesum/lensum ;
     }

     // Loop through faces.
     virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitialFaceCenterGridMotion>
    registerInitialFaceCenterGridMotion ;

  class InitialAreaGridMotion : public pointwise_rule {
    private:
      const_store<vect3d> pos_ic ;
      const_store<vect3d> facecenter_ic ;
      const_multiMap face2node ;
      store<Area> area_ic ;
    public:
                                                                                
      // Define input and output.
      InitialAreaGridMotion() {
        name_store("pos_ic",pos_ic) ;
        name_store("face2node",face2node) ;
        name_store("facecenter_ic",facecenter_ic) ;
        name_store("gridMotion::area_ic",area_ic) ;
        input("facecenter_ic,face2node->pos_ic") ;
        output("gridMotion::area_ic") ;
      }
                                                                                
      // Area for a single face.
      void calculate(Entity fc) {
        const vect3d center=facecenter_ic[fc] ;
        const int first=*(face2node.begin(fc)) ;
        const int last=*(face2node.end(fc)-1) ;
        vect3d sum(cross(pos_ic[last]-center,pos_ic[first]-center)) ;
        for(const int* ni=face2node.begin(fc)+1;ni!=(face2node.end(fc));ni++) {
          sum+=cross(pos_ic[*(ni-1)]-center,pos_ic[*(ni)]-center) ;
        }
        area_ic[fc].n=sum ;
        const real sada=sqrt(dot(area_ic[fc].n,area_ic[fc].n)) ;
        area_ic[fc].n*=1./(sada+EPSILON) ; area_ic[fc].sada=0.5*sada ;
      }
                                                                                
      // Loop over faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<InitialAreaGridMotion> registerInitialAreaGridMotion ;


//-----------------------------------------------------------------------------
// Initial conditions for node displacement.

  // Default node displacement initial condition.
  class NodeDisplacementInitialCondition : public pointwise_rule {
    private:
      store<vect3d> s_ic ;
    public:

      // Define input and output.
      NodeDisplacementInitialCondition() {
        name_store("node_s_ic",s_ic) ;
        output("node_s_ic") ;
        constraint("pos,noRestart") ;
      }

      // Set displacement for single node.
      void calculate(Entity n) { s_ic[n]=vect3d(0.0,0.0,0.0) ; }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<NodeDisplacementInitialCondition>
    registerNodeDisplacementInitialCondition ;

  // Node displacement initial condition for restart.
  class NodeDisplacementInitialConditionRestart : public pointwise_rule {
    private:
      const_param<string> restartNum ;
      store<vect3d> s_ic ;
    public:

      // Define input and output.
      NodeDisplacementInitialConditionRestart() {
        name_store("restartNum",restartNum) ;
        name_store("node_s_ic",s_ic) ;
        input("restartNum") ;
        output("node_s_ic") ;
        constraint("restart,pos") ;
      }
                                                                                
      // Read the node displacement values from file.
      void compute(const sequence &seq) {
        ostringstream oss ; oss << "restart/node_s_hdf5." << *restartNum ;
        string fileName=oss.str() ;
        if(Loci::MPI_rank==0) cout << "Reading node displacements "
          << "from restart file: " << fileName << endl ;
        entitySet dom=entitySet(seq) ;
        hid_t fileID=Loci::hdf5OpenFile(fileName.c_str(),H5F_ACC_RDONLY,
          H5P_DEFAULT);
        Loci::readContainer(fileID,"node_s",s_ic.Rep(),dom) ;
        Loci::hdf5CloseFile(fileID) ;
      }
  } ;
                                                                                
  register_rule<NodeDisplacementInitialConditionRestart>
    registerNodeDisplacementInitialConditionRestart ;

//-----------------------------------------------------------------------------
// Rules that define the iteration history of the node position.

  // Default node position initial condition.
  class NodePositionInitialCondition : public pointwise_rule {
    private:
      const_store<vect3d> pos ;
      const_store<vect3d> s_ic ;
      store<vect3d> pos_ic ;
    public:

      // Define input and output.
      NodePositionInitialCondition() {
        name_store("pos",pos) ;
        name_store("node_s_ic",s_ic) ;
        name_store("pos_ic",pos_ic) ;
        input("pos,node_s_ic") ;
        output("pos_ic") ;
      }

      // Set position for single node.
      void calculate(Entity n) { pos_ic[n]=pos[n]+s_ic[n] ; }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<NodePositionInitialCondition>
    registerNodePositionInitialCondition ;

  // Time build rule for node position.
  class TimeBuildNodePositionBDF2 : public pointwise_rule {
    private:
      const_store<vect3d> pos_ic ;
      store<vect3d> pos ;
    public:

      // Define input and output.
      TimeBuildNodePositionBDF2() {
        name_store("pos_ic",pos_ic) ;
        name_store("pos{n=-1}",pos) ;
        input("pos_ic") ;
        output("pos{n=-1}") ;
      }

      // Set position for single node.
      void calculate(Entity n) { pos[n]=pos_ic[n] ; }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildNodePositionBDF2> registerTimeBuildNodePositionBDF2 ;

  // Time build rule for node position.
  class TimeBuildNodePosition : public pointwise_rule {
    private:
      const_store<vect3d> pos_ic ;
      store<vect3d> pos ;
    public:

      // Define input and output.
      TimeBuildNodePosition() {
        name_store("pos_ic",pos_ic) ;
        name_store("pos{n=0}",pos) ;
        input("pos_ic") ;
        output("pos{n=0}") ;
      }

      // Set position for single node.
      void calculate(Entity n) { pos[n]=pos_ic[n] ; }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildNodePosition> registerTimeBuildNodePosition ;

  // /* CK 4/9/2010
  // Iteration build rule for node position for time-dependent grid motion,
  // where the node positions have already been updated before going
  // into the "it" loop.
  class IterationBuildNodePositionGridMotionTimeDependent : public pointwise_rule {
    private:
      const_store<vect3d> posStar ;
      store<vect3d> pos ;
    public:

      // Define input and output.
      IterationBuildNodePositionGridMotionTimeDependent() {
        name_store("posStar{n}",posStar) ;
        name_store("pos{n+1}",pos) ;
        input("posStar{n}") ;
        output("pos{n+1}") ;
        constraint("pos,gridMotionTimeDependent") ;
      }

      // Assign position for a single node.
      void calculate(Entity n) { pos[n]=posStar[n] ; }

      // Loop over nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationBuildNodePositionGridMotionTimeDependent>
    registerIterationBuildNodePositionGridMotionTimeDependent ;
  //*/

//-----------------------------------------------------------------------------
// Rules that define the iteration history of the node displacement.

  // Time build rule for node displacement.
  class TimeBuildNodeDisplacement : public pointwise_rule {
    private:
      const_store<vect3d> s_ic ;
      store<vect3d> s ;
    public:

      // Define input and output.
      TimeBuildNodeDisplacement() {
        name_store("node_s_ic",s_ic) ;
        name_store("node_s{n=0}",s) ;
        input("node_s_ic") ;
        output("node_s{n=0}") ;
      }

      // Set displacement for single node.
      void calculate(Entity n) { s[n]=s_ic[n] ; }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildNodeDisplacement> registerTimeBuildNodeDisplacement ;


  // Updated node positions for time-dependent grid motion.
  class UpdateNodePositionGridMotionTimeDependent : public pointwise_rule {
    private:
      const_store<vect3d> posRef ;
      const_store<vect3d> node_sStar ;
      store<vect3d> posStar ;
    public:

      // Define input and output.
      UpdateNodePositionGridMotionTimeDependent() {
        name_store("pos",posRef) ;
        name_store("node_sStar{n}",node_sStar) ;
        name_store("posStar{n}",posStar) ;
        input("pos,node_sStar{n}") ;
        output("posStar{n}") ;
        constraint("pos,gridMotionTimeDependent") ;
      }

      // Position for single node.
      void calculate(Entity n) { posStar[n]=posRef[n]+node_sStar[n] ; }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<UpdateNodePositionGridMotionTimeDependent>
    registerUpdateNodePositionGridMotionTimeDependent ;

  // Time advance rule for node displacement for time-dependent grid motion.
  class TimeAdvanceNodeDisplacementGridMotionTimeDependent : public
  pointwise_rule {
    private:
      store<vect3d> s ;
    public:
                                                                                
      // Define input and output.
      TimeAdvanceNodeDisplacementGridMotionTimeDependent() {
        name_store("node_sStar{n}",s) ;
        input("node_sStar{n}") ;
        output("node_s{n+1}=node_sStar{n}") ;
      }

      // Empty compute.
      virtual void compute(const sequence &seq) {}
  } ;
                                                                                
  register_rule<TimeAdvanceNodeDisplacementGridMotionTimeDependent>
    registerTimeAdvanceNodeDisplacementGridMotionTimeDependent ;

//-----------------------------------------------------------------------------
// Generic rules independent of the grid motion type.

  // Compute old node velocities. This rule assumes that we have a constant
  // timestep.
  class OldNodeVelocity : public pointwise_rule {
    private:
      const_param<real> timeStep ;
      const_store<vect3d> posOld,pos ;
      store<vect3d> node_vOld ;
    public:

      // Define input and output.
      OldNodeVelocity() {
        name_store("timeStep{n}",timeStep) ;
        name_store("pos{n-1}",posOld) ;
        name_store("pos{n}",pos) ;
        name_store("node_vOld{n}",node_vOld) ;
        input("timeStep{n},pos{n-1},pos{n}") ;
        output("node_vOld{n}") ;
      }

      // Position for single node.
      void calculate(Entity n) {
        node_vOld[n]=(pos[n]-posOld[n])/(*timeStep) ;
      }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<OldNodeVelocity> registerOldNodeVelocity ;

  // Compute new node velocities for solution-dependent grid motion. Note
  // that pos{n,it} are the updated positions.
  class NodeVelocity : public pointwise_rule {
    private:
      const_param<real> timeStep ;
      const_store<vect3d> posOld ;
      const_store<vect3d> pos ;
      store<vect3d> node_v ;
    public:

      // Define input and output.
      NodeVelocity() {
        name_store("timeStep{n}",timeStep) ;
        name_store("pos{n}",posOld) ;
        name_store("posStar{n}",pos) ;
        name_store("node_v{n}",node_v) ;
        input("timeStep{n},pos{n},posStar{n}") ;
        output("node_v{n}") ;
        constraint("pos") ;
      }

      // Position for single node.
      void calculate(Entity n) { node_v[n]=(pos[n]-posOld[n])/(*timeStep) ; }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<NodeVelocity> registerNodeVelocity ;

  // Default rule for face velocities.
  class GridFaceVelocityDefault : public pointwise_rule {
    private:
      store<vect3d> face_v ;
    public:

      // Define input and output.
      GridFaceVelocityDefault() {
        name_store("face_v",face_v) ;
        output("face_v") ;
        constraint("faces") ;
      }

      // Velocity for a single face.
      void calculate(Entity face) { face_v[face]=vect3d(0.0,0.0,0.0) ; }

      // Loop through faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GridFaceVelocityDefault> registerGridFaceVelocityDefault ;

  // Compute grid face velocity.
  class GridFaceVelocityPriority : public pointwise_rule {
    private:
      const_multiMap face2node ;
      const_store<vect3d> node_v ;
      store<vect3d> face_v ;
    public:

      // Define input and output.
      GridFaceVelocityPriority() {
        name_store("face2node",face2node) ;
        name_store("node_v",node_v) ;
        name_store("priority::face_v",face_v) ;
        input("face2node->node_v") ;
        output("priority::face_v") ;
      }

      // Velocity for a single face.
      void calculate(Entity face) {
        vect3d sum(0.0,0.0,0.0) ;
        for(const int *ni=face2node.begin(face);ni!=face2node.end(face);++ni)
          sum+=node_v[*ni] ;
        face_v[face]=(1.0/real(face2node.num_elems(face)))*sum ;
      }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GridFaceVelocityPriority> registerGridFaceVelocityPriority ;

//-----------------------------------------------------------------------------
// Rules for producing output data.

  // Magnitude of node displacement.
  class NodeDisplacementMagnitude : public pointwise_rule {
    private:
      const_store<vect3d> node_s ;
      store<float> sMag ;
    public:

      // Define input and output.
      NodeDisplacementMagnitude() {
        name_store("node_s",node_s) ;
        name_store("sMag",sMag) ;
        input("node_s") ;
        output("sMag") ;
        constraint("pos") ;
      }

      // Set displacement magnitude for single node.
      void calculate(Entity n) { sMag[n]=norm(node_s[n]) ; }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<NodeDisplacementMagnitude> registerNodeDisplacementMagnitude ;
}
