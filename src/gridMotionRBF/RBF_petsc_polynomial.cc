//-----------------------------------------------------------------------------
// Description: This file contains rules for implementing a RBF linear solver
// using PETSC.
//-----------------------------------------------------------------------------

// mpi
#include<mpi.h>

// Standard library includes.
#include<vector>
#include<string>
using std::vector ;

// Loci includes.
#include<Loci>
using Loci::singleton_rule ;
using Loci::pointwise_rule ;
using Loci::unit_rule ;
using Loci::apply_rule ;
using Loci::sequence ;
using Loci::entitySet ;
using Loci::EMPTY ;
using Loci::Entity ;
using Loci::register_rule ;
using Loci::blackbox ;
using Loci::const_blackbox ;
using Loci::param ;
using Loci::const_param ;
using Loci::store ;
using Loci::const_store ;
using Loci::const_Map ;
using Loci::const_multiMap ;
using Loci::const_MapVec ;
using Loci::const_storeVec ;

// PETSC includes.
#include <petsc.h>
#include <petscerror.h>
#include <petscksp.h>
#include <petscpc.h>

// StreamUns includes.
#include "sciTypes.h"

// RBF includes.
#include "rbf.h"

// Create a namespace for the PETSC objects so there is no confusion. Note that
// there is already a type "Mat" defined by Loci.
namespace Petsc {
  typedef struct _p_Vec* Vec ;
  typedef struct _p_Mat* Mat ;
  typedef struct _p_KSP* KSP ;
} ;

namespace streamUns {

	// delimiter
	string delim=", " ;
//-----------------------------------------------------------------------------
// Wrapper classes for PETSC objects. These have been created primarily to
// handle memory management.
  
  // Wrapper for Vec. Note that we cannot check the success of the destroy in
  // the destructor since the destructor cannot return a value.
  class PETSCVector {
    private:
      bool created ;
      mutable Petsc::Vec v ;
    public:
      PETSCVector() : created(false) {}
      ~PETSCVector() { if(created) VecDestroy(v) ; }
    public:
      int AssemblyBegin() const {
        int ierr=VecAssemblyBegin(v) ; CHKERRQ(ierr) ; return 0 ;
      }
      int AssemblyEnd() const {
        int ierr=VecAssemblyEnd(v) ; CHKERRQ(ierr) ; return 0 ;
      }
      int Create(int localSize,int globalSize) {
        if(created) return 0 ;
        int ierr=VecCreate(PETSC_COMM_WORLD,&v) ; CHKERRQ(ierr) ;
        ierr=VecSetSizes(v,localSize,globalSize) ; CHKERRQ(ierr) ;
        ierr=VecSetFromOptions(v) ; CHKERRQ(ierr) ;
        created=true ; return 0 ;
      }
      const Petsc::Vec& Data() const { return v ; }
      int DuplicateFrom(const PETSCVector &a) {
        if(created) return 0 ; int ierr=VecDuplicate(a.v,&v) ; CHKERRQ(ierr) ;
        created=true ; return 0 ;
      }
      int GetArray(PetscScalar **array) const {
         int ierr=VecGetArray(v,array) ; CHKERRQ(ierr) ; return 0 ;
      }
      int GetOwnershipRange(int *low,int *high) const {
        int ierr=VecGetOwnershipRange(v,low,high) ; CHKERRQ(ierr) ; return 0 ;
      }
      int RestoreArray(PetscScalar **array) const {
        int ierr=VecRestoreArray(v,array) ; CHKERRQ(ierr) ; return 0 ;
      }
      int SetValue(const int *index,const PetscScalar *value) const {
        int ierr=VecSetValues(v,1,index,value,INSERT_VALUES) ; CHKERRQ(ierr) ;
        return 0 ;
      }
  } ;

  // Wrapper for Vec. In this version we have three vectors which are used in
  // the vector solver rules. Note that we cannot check the success of the
  // destroy in the destructor since the destructor cannot return a value.
  class PETSCMultiVector {
    private:
      bool created ;
      mutable Petsc::Vec vX,vY,vZ ;
    public:
      PETSCMultiVector() : created(false) {}
      ~PETSCMultiVector() {
        if(created){ VecDestroy(vX) ; VecDestroy(vY) ; VecDestroy(vZ) ; }
      }
    public:
      int AssemblyBegin() const {
        int ierr=VecAssemblyBegin(vX) ; CHKERRQ(ierr) ;
        ierr=VecAssemblyBegin(vY) ; CHKERRQ(ierr) ;
        ierr=VecAssemblyBegin(vZ) ; CHKERRQ(ierr) ; return 0 ;
      }
      int AssemblyEnd() const {
        int ierr=VecAssemblyEnd(vX) ; CHKERRQ(ierr) ;
        ierr=VecAssemblyEnd(vY) ; CHKERRQ(ierr) ;
        ierr=VecAssemblyEnd(vZ) ; CHKERRQ(ierr) ; return 0 ;
      }
      int Create(int localSize,int globalSize) {
        if(created) return 0 ;
        int ierr=VecCreate(PETSC_COMM_WORLD,&vX) ; CHKERRQ(ierr) ;
        ierr=VecSetSizes(vX,localSize,globalSize) ; CHKERRQ(ierr) ;
        ierr=VecSetFromOptions(vX) ; CHKERRQ(ierr) ;
        ierr=VecCreate(PETSC_COMM_WORLD,&vY) ; CHKERRQ(ierr) ;
        ierr=VecSetSizes(vY,localSize,globalSize) ; CHKERRQ(ierr) ;
        ierr=VecSetFromOptions(vY) ; CHKERRQ(ierr) ;
        ierr=VecCreate(PETSC_COMM_WORLD,&vZ) ; CHKERRQ(ierr) ;
        ierr=VecSetSizes(vZ,localSize,globalSize) ; CHKERRQ(ierr) ;
        ierr=VecSetFromOptions(vZ) ; CHKERRQ(ierr) ;
        created=true ; return 0 ;
      }
      const Petsc::Vec& DataX() const { return vX ; }
      const Petsc::Vec& DataY() const { return vY ; }
      const Petsc::Vec& DataZ() const { return vZ ; }
      int DuplicateFrom(const PETSCMultiVector &a) {
        if(created) return 0 ;
        int ierr=VecDuplicate(a.vX,&vX) ; CHKERRQ(ierr) ;
        ierr=VecDuplicate(a.vY,&vY) ; CHKERRQ(ierr) ;
        ierr=VecDuplicate(a.vZ,&vZ) ; CHKERRQ(ierr) ;
        created=true ; return 0 ;
      }
      int GetArray(PetscScalar **arrayX,PetscScalar **arrayY,PetscScalar
      **arrayZ) const {
         int ierr=VecGetArray(vX,arrayX) ; CHKERRQ(ierr) ;
         ierr=VecGetArray(vY,arrayY) ; CHKERRQ(ierr) ;
         ierr=VecGetArray(vZ,arrayZ) ; CHKERRQ(ierr) ; return 0 ;
      }
      int GetOwnershipRange(int *low,int *high) const {
        int ierr=VecGetOwnershipRange(vX,low,high) ; CHKERRQ(ierr) ; return 0 ;
      }
      int RestoreArray(PetscScalar **arrayX,PetscScalar **arrayY,PetscScalar
      **arrayZ) const {
        int ierr=VecRestoreArray(vX,arrayX) ; CHKERRQ(ierr) ;
        ierr=VecRestoreArray(vY,arrayY) ; CHKERRQ(ierr) ;
        ierr=VecRestoreArray(vZ,arrayZ) ; CHKERRQ(ierr) ; return 0 ;
      }
      int SetValue(const int *index,const PetscScalar *valueX,const PetscScalar
      *valueY,const PetscScalar *valueZ) const {
        int ierr=VecSetValues(vX,1,index,valueX,INSERT_VALUES) ; CHKERRQ(ierr) ;
        ierr=VecSetValues(vY,1,index,valueY,INSERT_VALUES) ; CHKERRQ(ierr) ;
        ierr=VecSetValues(vZ,1,index,valueZ,INSERT_VALUES) ; CHKERRQ(ierr) ;
        return 0 ;
      }
  } ;

  // Wrapper for Mat. This class assumes a square matrix. Note that we cannot
  // check the success of the destroy in the destructor since the destructor
  // cannot return a value.
  class PETSCMatrix {
    private:
      bool created ;
//      Petsc::Mat m ;
      mutable Petsc::Mat m ;
    public:
      PETSCMatrix() : created(false) {}
      ~PETSCMatrix() { if(created) MatDestroy(m) ; }
    public:
      int AssemblyBegin() const {
//		cout << "i'm here assembly begin,r: " << Loci::MPI_rank << endl ;
        int ierr=MatAssemblyBegin(m,MAT_FINAL_ASSEMBLY) ;
//		cout << "i'm here assembly begin after,r: " << Loci::MPI_rank << endl ;
        CHKERRQ(ierr) ; return ierr ;
      }
      int AssemblyEnd() const {
//		cout << "i'm here assembly end,r: " << Loci::MPI_rank << endl ;
       int ierr=MatAssemblyEnd(m,MAT_FINAL_ASSEMBLY) ;
//		cout << "i'm here assembly end after,r: " << Loci::MPI_rank << endl ;
       CHKERRQ(ierr) ; return ierr ;
      }
      int Create(int numLocalRow,int numGlobalRow,int *numDiagonalNonZero,int
      *numOffDiagonalNonZero) {
        if(created) return 0 ;
		int ierr ;
        if(Loci::MPI_processes>1){
          ierr=MatCreateMPIAIJ(PETSC_COMM_WORLD,numLocalRow,numLocalRow,
            numGlobalRow,numGlobalRow,0,numDiagonalNonZero,0,
            numOffDiagonalNonZero,&m) ; CHKERRQ(ierr) ;
        }else{
          ierr=MatCreateSeqAIJ(PETSC_COMM_WORLD,numGlobalRow,numGlobalRow,
            0,numDiagonalNonZero,&m) ; CHKERRQ(ierr) ;
        }
        created=true ; return ierr ;
      }
      const Petsc::Mat& Data() const { return m ; }
	  int GetInfo() {
		  MatInfo info;
		  double mal, nz_alloc, nz_used, nz_un ;
		  int ierr = MatGetInfo(m,MAT_LOCAL,&info); CHKERRQ(ierr) ;
//		cout << "i'm here getinfo, r: " << Loci::MPI_rank << endl ;
		  mal = info.mallocs; nz_alloc=info.nz_allocated; nz_used=info.nz_used; nz_un=info.nz_unneeded ;
		  cout << "PETSC info: rank= " << Loci::MPI_rank << " Number of mallocs during MatSetValues() --> " << mal << delim << "Number of nonzeros (allocated, used, unneeded): " << nz_alloc << delim << nz_used << delim << nz_un << endl ;
	return ierr ;
	  }
      int GetOwnershipRange(int *firstRow,int *lastRow) {
        int ierr=MatGetOwnershipRange(m,firstRow,lastRow) ; CHKERRQ(ierr) ;
        return 0 ;
      }
      int SetRowValues(int rowNum,int numColumn,int *columnIndices,PetscScalar
      *columnValues) {
        int ierr=MatSetValues(m,1,&rowNum,numColumn,columnIndices,columnValues,
          INSERT_VALUES) ; CHKERRQ(ierr) ; return ierr ;
      }
  } ;
  
  // Wrapper for KSP.
  class PETSCKsp {
    private:
      bool created ;
      mutable Petsc::KSP ksp ;
    public:
      PETSCKsp() : created(false) {}
      ~PETSCKsp() { if(created) KSPDestroy(ksp) ; }
    public:
      int Create() {
        if(created) return 0 ;
        int ierr=KSPCreate(PETSC_COMM_WORLD,&ksp) ; CHKERRQ(ierr) ;
        created=true ; return 0 ;
      }
	  int GetConvergedReason() const {
		  KSPConvergedReason reason ;
		  int ierr = KSPGetConvergedReason(ksp, &reason) ;
		  if (Loci::MPI_rank==0) {
			  if (reason==KSP_CONVERGED_RTOL) {
				  cout << "PETSC solver: GetConvergedReason = " << reason << " KSP_CONVERGED_RTOL" << endl ;
			  } else if (reason==KSP_CONVERGED_ATOL) {
				  cout << "PETSC solver: GetConvergedReason = " << reason << " KSP_CONVERGED_ATOL" << endl ;
			  } else if (reason==KSP_CONVERGED_ITS) {
				  cout << "PETSC solver: GetConvergedReason = " << reason << " KSP_CONVERGED_ITS" << endl ;
			  } else if (reason==KSP_DIVERGED_NULL) {
				  cout << "PETSC solver: GetConvergedReason = " << reason << " KSP_DIVERGED_NULL" << endl ;
			  } else if (reason==KSP_DIVERGED_ITS) {
				  cout << "PETSC solver: GetConvergedReason = " << reason << " KSP_DIVERGED_ITS" << endl ;
			  } else if (reason==KSP_DIVERGED_DTOL) {
				  cout << "PETSC solver: GetConvergedReason = " << reason << " KSP_DIVERGED_DTOL" << endl ;
			  } else {
				  cout << "PETSC solver: GetConvergedReason = " << reason << endl ;
			  }
		  }
		  return 0;
	  }
      int GetIterationNumber() const {
        int numIterations ; KSPGetIterationNumber(ksp,&numIterations) ;
        return numIterations ;
      }
      int GetPC(PC *pc) const {
        int ierr=KSPGetPC(ksp,pc) ; CHKERRQ(ierr) ; return 0 ;
      }
	  int SetInitialGuessNonzero() const {
        int ierr=KSPSetInitialGuessNonzero(ksp, PETSC_TRUE) ; CHKERRQ(ierr) ; return 0 ;
      }
      int SetFromOptions() {
        int ierr=KSPSetFromOptions(ksp) ; CHKERRQ(ierr) ; return 0 ;
      }
	  int SetRestart(int Restart) {
		  int ierr=KSPGMRESSetRestart(ksp, Restart) ;
		  CHKERRQ(ierr) ; return ierr ;
	  }
      int SetOperators(const PETSCMatrix &m) const {
        int ierr=KSPSetOperators(ksp,m.Data(),m.Data(),DIFFERENT_NONZERO_PATTERN) ; // SAME_NONZERO_PATTERN
        CHKERRQ(ierr) ; return 0 ;
      }
      void SetTolerances(real relativeTolerance,real absoluteTolerance,int maxIterations) const {
        KSPSetTolerances(ksp,relativeTolerance,absoluteTolerance,PETSC_DEFAULT,maxIterations) ;
      }
	  void SetMINRES() {
		  KSPSetType(ksp, KSPMINRES) ;
	  }
      int Solve(const PETSCVector &rhs, PETSCVector &sol) const {
//      int numIteration ;
        int ierr=KSPSolve(ksp,rhs.Data(),sol.Data()) ; CHKERRQ(ierr) ;
//KSPGetIterationNumber(ksp,&numIteration) ;
//cout << "petsc iterations: " << numIteration << endl ;
        return ierr ;
      }
      int SolveVector(const PETSCMultiVector &rhs, PETSCMultiVector &sol)
      const {
//      int numIteration ;
        int ierr=KSPSolve(ksp,rhs.DataX(),sol.DataX()) ; CHKERRQ(ierr) ;
        ierr=KSPSolve(ksp,rhs.DataY(),sol.DataY()) ; CHKERRQ(ierr) ;
        ierr=KSPSolve(ksp,rhs.DataZ(),sol.DataZ()) ; CHKERRQ(ierr) ;
//KSPGetIterationNumber(ksp,&numIteration) ;
//cout << "petsc iterations: " << numIteration << endl ;
        return ierr ;
      }
  } ;

  
//-----------------------------Nonzero entries---------------------------------------------------------
//
  // Determines the number of non-zero entries in the local portion of the
  // Petsc matrix for each node. The local portion is defined as the square
  // row/column sub-block of the PETSC matrix whose rows map to nodes local
  // to the process.
  class rbfBoundaryNodeNumDiagonalNonZeroEntries : public pointwise_rule {
    private:
	  const_param<real> r ;
	  const_param<int> rbfNumber ;	  
	  const_param<vector<real> > Q ;
      const_param<vector<int> > rbfNumBoundaryNode ;
      store<int> rbfBoundaryNodeNumDiagonalNonZero ;
	  store<int> rbfBoundaryNodeNumOffDiagonalNonZero ;
    public:

      // Define input and output.
      rbfBoundaryNodeNumDiagonalNonZeroEntries() {
        name_store("rbfNumBoundaryNode",rbfNumBoundaryNode) ;
        name_store("rbfBoundaryNodeNumDiagonalNonZero",rbfBoundaryNodeNumDiagonalNonZero) ; 
		name_store("rbfBoundaryNodeNumOffDiagonalNonZero",rbfBoundaryNodeNumOffDiagonalNonZero) ;
		name_store("rbfR",r) ;
		name_store("rbfNumber",rbfNumber) ;
		name_store("rbfBoundaryNodeQ",Q) ;
		input("rbfR, rbfNumber") ;
		input("rbfBoundaryNodeQ") ;
        output("rbfBoundaryNodeNumDiagonalNonZero") ;
		output("rbfBoundaryNodeNumOffDiagonalNonZero") ;
        constraint("boundaryDisplacement") ;
        constraint("sStar_PETSCLinearSolver") ;
		disable_threading() ;
      }

      // Loop over nodes.
      void compute(const sequence & seq) { 
		
		int rank=Loci::MPI_rank ;
		
		// Get the number of local and global nodes.
        int globalNumNode=0, localStart=0  ;
        for(unsigned int i=0;i<(*rbfNumBoundaryNode).size();++i) globalNumNode+=(*rbfNumBoundaryNode)[i] ;
        for(int i=0;i<rank;++i) localStart+=(*rbfNumBoundaryNode)[i]  ;  // 3 coordinates
		int localNum=(*rbfNumBoundaryNode)[rank]; 
		
		double distance = 0.;
		int localEnd = localStart+localNum ;
		sequence::const_iterator nodePtr;
//		for(int i=localStart;i<localEnd;++i){	// Q[row,col]=Q[3*row+col], Q=[0,1,2;3,4,5;...]
		int counter=localStart ;
		int I,I1,I2,J,J1,J2;
		int counterD=0; int counterOD=0;
		for(nodePtr=seq.begin();nodePtr!=seq.end();++nodePtr,++counter){	// Q[row,col]=Q[3*row+col], Q=[0,1,2;3,4,5;...]
			I=3*counter; I1=3*counter+1; I2=3*counter+2 ;
			for(int j=0;j<globalNumNode;++j){
				J=3*j; J1=3*j+1; J2=3*j+2 ;
				rbfBoundaryNodeNumDiagonalNonZero[*nodePtr]=0;
				rbfBoundaryNodeNumOffDiagonalNonZero[*nodePtr]=0;
				if( (j>=localStart) && (j<localEnd) ) { // nodePtr2 is in the same partition as nodePtr1 -> diagonal
					if (*rbfNumber > 0 && *rbfNumber < 9) { // compact
						distance = pow((*Q)[I] - (*Q)[J], 2.) ; // ((*Q)[I] - (*Q)[J])*((*Q)[I] - (*Q)[J]) ;
						distance += pow((*Q)[I1] - (*Q)[J1], 2.) ; // ((*Q)[I1] - (*Q)[J1])*((*Q)[I1] - (*Q)[J1]) ;
						distance += pow((*Q)[I2] - (*Q)[J2], 2.) ; // ((*Q)[I2] - (*Q)[J2])*((*Q)[I2] - (*Q)[J2]) ;
						distance = sqrt(distance) ;
						distance /= (*r) ;
						if (distance <=1.) {
//							cout << "before plus?" << counterD << endl ;
//							++(rbfBoundaryNodeNumDiagonalNonZero[*nodePtr]);
//							(rbfBoundaryNodeNumDiagonalNonZero[*nodePtr])=(rbfBoundaryNodeNumDiagonalNonZero[*nodePtr])+1;
							++counterD;
//							cout << "after plus? " << counterD << endl ;
//							cout << "after plus? " << rbfBoundaryNodeNumDiagonalNonZero[*nodePtr] << endl ;
						} 
//						cout << "counter,j,distance,diagNonZero: " << counter << ", " << j << ", " << distance << ", " << rbfBoundaryNodeNumDiagonalNonZero[*nodePtr] << endl ;
					} else { // global
//						++(rbfBoundaryNodeNumDiagonalNonZero[*nodePtr]);
						++counterD;
					}
				} else { // nodePtr2 is not in the same partition as nodePtr2 -> off diagonal
					if (*rbfNumber > 0 && *rbfNumber < 9) { // compact
						distance = pow((*Q)[I] - (*Q)[J], 2.) ; // ((*Q)[I] - (*Q)[J])*((*Q)[I] - (*Q)[J]) ;
						distance += pow((*Q)[I1] - (*Q)[J1], 2.) ; // ((*Q)[I1] - (*Q)[J1])*((*Q)[I1] - (*Q)[J1]) ;
						distance += pow((*Q)[I2] - (*Q)[J2], 2.) ; // ((*Q)[I2] - (*Q)[J2])*((*Q)[I2] - (*Q)[J2]) ;
						distance /= ((*r)*(*r)) ;
						if (distance <=1) {
//							++(rbfBoundaryNodeNumOffDiagonalNonZero[*nodePtr]);
							++counterOD;
						}
					} else { // global
//						++(rbfBoundaryNodeNumOffDiagonalNonZero[*nodePtr]);
						++counterOD;
					}
				}
			}			
			rbfBoundaryNodeNumDiagonalNonZero[*nodePtr]=counterD;
			rbfBoundaryNodeNumOffDiagonalNonZero[*nodePtr]=counterOD;
//			cout << "counter,counterD,rbfNumDiagNonZEro: " << counter << ", " << counterD << ", " << rbfBoundaryNodeNumDiagonalNonZero[*nodePtr] << endl ;
			counterD=0; counterOD=0; // reset
		}
	  }
  } ;

  register_rule<rbfBoundaryNodeNumDiagonalNonZeroEntries>
    registerrbfBoundaryNodeNumDiagonalNonZeroEntries ;

//-----------------------------Setup linear solver---------------------------------------------------------
//
  // Sets up the PETSC linear solver.
  class PETSCBoundaryNodeSetupSolver : public singleton_rule {
    private:
      const_param<int> maxLinearSolverIterations ;
      const_param<real> linearSolverTolerance ;
      blackbox<PETSCKsp> ksp ;
    public:

      // Define input and output.
      PETSCBoundaryNodeSetupSolver() {
        name_store("sStar_maxLinearSolverIterations",maxLinearSolverIterations) ;
        name_store("sStar_linearSolverTolerance",linearSolverTolerance) ;
        name_store("petscBoundaryNodeKSP",ksp) ;
        input("sStar_maxLinearSolverIterations") ;
        input("sStar_linearSolverTolerance") ;
        output("petscBoundaryNodeKSP") ;
        constraint("sStar_PETSCLinearSolver") ;
		constraint("boundaryDisplacement") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {
        (*ksp).Create() ;
//		(*ksp).SetMINRES() ; // use minres
//		  (*ksp).SetInitialGuessNonzero() ; // Set initial guess using the values
        (*ksp).SetTolerances(*linearSolverTolerance,1.0e-10,*maxLinearSolverIterations) ;
	  (*ksp).SetRestart(30) ;
        (*ksp).SetFromOptions() ;
      // Valve case.
//      PC pc ; (*ksp).GetPC(&pc) ; PCSetType(pc,PCILU) ; // PCFactorSetShiftType(pc,MAT_SHIFT_POSITIVE_DEFINITE) ;//PCILUSetLevels(pc,3) ;
      // Injector case.
//      PC pc ; (*ksp).GetPC(&pc) ; PCSetType(pc,KSPCG); // PCSetType(pc,PCJACOBI) ;
      PC pc ; (*ksp).GetPC(&pc) ; PCSetType(pc,PCJACOBI); // PCSetType(pc,PCJACOBI) ;
//      PC pc ; (*ksp).GetPC(&pc) ; PCSetType(pc,PCSOR); // PCSetType(pc,PCJACOBI) ;
      // PCILUSetLevels(pc,0) ;
      }
  } ;

  register_rule<PETSCBoundaryNodeSetupSolver> registerPETSCBoundaryNodeSetupSolver ;

//---------------------------------------------b-------------------------------------------------
//
//
  // Sets up the PETSC right-hand-side vector.
  class PETSCBoundaryNodeUnitRHS : public unit_rule {
    private:
      const_param<vector<int> > rbfNumBoundaryNode ;
      const_store<vect3d> B ;
      blackbox<PETSCMultiVector> b;
    public:

      // Define input and output.
      PETSCBoundaryNodeUnitRHS() {
        name_store("rbfNumBoundaryNode",rbfNumBoundaryNode) ;
        name_store("petscBoundaryNodeB",b) ;
        name_store("sStar_B",B);
        input("rbfNumBoundaryNode") ;
		input("sStar_B");
        output("petscBoundaryNodeB") ;
        constraint("sStar_PETSCLinearSolver") ;
		constraint("boundaryDisplacement") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {

        // Get the number of local and global nodes.
        int localNumNode=(*rbfNumBoundaryNode)[Loci::MPI_rank],globalNumNode=0 ;
        for(unsigned int i=0;i<(*rbfNumBoundaryNode).size();++i) globalNumNode+=(*rbfNumBoundaryNode)[i] ;
		
		// include the polynomial coeff. computation
		if (Loci::MPI_rank == Loci::MPI_processes-1) { // last rank takes care of beta
			localNumNode += 4; // d+1 = 4
		}
					
		globalNumNode+=4; // for beta (d+1)
		
        // Allocate the unknown and rhs vectors.
		(*b).Create(localNumNode,globalNumNode) ;

		// Compute the row offset for this process.
        int localStart=0 ;
        for(int i=0;i<Loci::MPI_rank;++i){ localStart+=(*rbfNumBoundaryNode)[i] ; }
		//
					
		sequence::const_iterator nodePtr;
		PetscScalar valuex, valuey, valuez;

//		cout << "rank, localStart, localNumNode, globalNumNode: " << Loci::MPI_rank << ", " << localStart << ", " << localNumNode << ", " << globalNumNode << endl;
		
		int counter=localStart ;
		for(nodePtr=seq.begin();nodePtr!=seq.end();++nodePtr,++counter){			
				// RHS vector
				valuex=B[*nodePtr].x;
				valuey=B[*nodePtr].y;
				valuez=B[*nodePtr].z;
//				cout << "rank, i,B.x,valuex,valuey " << Loci::MPI_rank << ", " << counter << ", " << B[*nodePtr].x << ", " << valuex << delim << B[*nodePtr].y << delim << valuey << endl;
				(*b).SetValue(&counter,&valuex, &valuey, &valuez) ;
		}		
      }
  } ;
 
  register_rule<PETSCBoundaryNodeUnitRHS> registerPETSCBoundaryNodeUnitRHS ;



  // Sets up the PETSC right-hand-side vector.
  class PETSCBoundaryNodeApplyRHS : public apply_rule<blackbox<PETSCVector>,Loci::NullOp<PETSCVector> > {
    private:
      const_param<vector<int> > rbfNumBoundaryNode ;
      const_store<vect3d> B ;
      blackbox<PETSCMultiVector> b;
    public:

      // Define input and output.
      PETSCBoundaryNodeApplyRHS() {
        name_store("rbfNumBoundaryNode",rbfNumBoundaryNode) ;
        name_store("petscBoundaryNodeB",b) ;
        name_store("sStar_B",B);
        input("rbfNumBoundaryNode") ;
		input("sStar_B");
        output("petscBoundaryNodeB") ;
        constraint("sStar_PETSCLinearSolver") ;
		constraint("boundaryDisplacement") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) { }
  } ;
 
  register_rule<PETSCBoundaryNodeApplyRHS> registerPETSCBoundaryNodeApplyRHS ;


//---------------------------------------------A-------------------------------------------------
//
  
  // Sets up the PETSC matrix. Note that this is a unit_rule, since this is the
  // only way we can have stores as input to a rule outputting blackboxes.
  class PETSCBoundaryNodeSetupMatrixUnit : public unit_rule {
    private:
	  const_param<real> r,a ;
	  const_param<int> rbfNr ;	  
	  const_store<vect3d> pos ;
//      const_store<int> rbfBoundaryNodeToRow ;
      const_store<int> rbfBoundaryNodeNumDiagonalNonZero ;
      const_store<int> rbfBoundaryNodeNumOffDiagonalNonZero ;
      const_param<vector<int> > rbfNumBoundaryNode ;
      const_param<vector<real> > Q ;
      blackbox<PETSCMatrix> A ;
	private:
      int columnIndex ;
      PetscScalar columnValue ;
    public:

      // Define input and output.
      PETSCBoundaryNodeSetupMatrixUnit() {
        name_store("rbfBoundaryNodeNumDiagonalNonZero",rbfBoundaryNodeNumDiagonalNonZero) ;
        name_store("rbfBoundaryNodeNumOffDiagonalNonZero",rbfBoundaryNodeNumOffDiagonalNonZero) ;
        name_store("rbfNumBoundaryNode",rbfNumBoundaryNode) ;
        name_store("petscBoundaryNodeA",A) ;		
		name_store("rbfBoundaryNodeQ",Q) ;
//		name_store("rbfBoundaryNodeToRow",rbfBoundaryNodeToRow) ;
		name_store("rbfR",r) ;
		name_store("rbfA",a) ;
		name_store("rbfNumber",rbfNr) ;
		name_store("pos",pos) ;
		input("rbfNumber, rbfR, rbfA");
//        input("rbfBoundaryNodeToRow") ;
        input("rbfBoundaryNodeQ") ;
		input("pos");
        input("rbfBoundaryNodeNumDiagonalNonZero,rbfBoundaryNodeNumOffDiagonalNonZero") ;
        input("rbfNumBoundaryNode") ;
//		input("rbfBoundaryNodeToRow") ;
        output("petscBoundaryNodeA") ;
        constraint("sStar_PETSCLinearSolver");
		constraint("boundaryDisplacement") ;
//		constraint("RBFxConstraints") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {

		int rank=Loci::MPI_rank ;
		int p=Loci::MPI_processes ;

        // Get the number of local and global nodes.
        int localNumNode=(*rbfNumBoundaryNode)[rank],globalNumNode=0 ;
        for(unsigned int i=0;i<(*rbfNumBoundaryNode).size();++i) globalNumNode+=(*rbfNumBoundaryNode)[i] ;
		  		
		// include the polynomial coeff. computation
		if (rank == p-1) { // last rank takes care of beta
			localNumNode += 4; // d+1 = 4
		}		
		globalNumNode+=4; // for beta (d+1)

		// Compute the row offset for this process.
        int localStart=0 ; 
        for(int i=0;i<rank;++i) localStart+=(*rbfNumBoundaryNode)[i]  ;  // 3 coordinates

//		cout << "Start Allocation A" << endl ;
        // Make temporary copy of matrix allocation data.
        int count=0,*numDiagonalNonZero=new int[localNumNode], *numOffDiagonalNonZero=new int[localNumNode] ;
		sequence::const_iterator nodePtr ;
        for(nodePtr=seq.begin();nodePtr!=seq.end();++nodePtr,++count){
          numDiagonalNonZero[count]=rbfBoundaryNodeNumDiagonalNonZero[*nodePtr] + ((rank==p-1)?4:0) ; // + Q, only r==p1 gets Q in diagNonZero, otherwise in offDiagNonZero
          numOffDiagonalNonZero[count]=rbfBoundaryNodeNumOffDiagonalNonZero[*nodePtr] + ((rank==p-1)?0:4) ;
//		  cout << "r,count,rbfDNZ,numDiag = " << rank << ", " << count << ", " << rbfBoundaryNodeNumDiagonalNonZero[*nodePtr] << ", " << numDiagonalNonZero[count] <<  endl ;
//          numDiagonalNonZero[rbfBoundaryNodeToRow[*nodePtr]]=rbfBoundaryNodeNumDiagonalNonZero[*nodePtr] ;
//          numOffDiagonalNonZero[rbfBoundaryNodeToRow[*nodePtr]]=rbfBoundaryNodeNumOffDiagonalNonZero[*nodePtr] ;
        }
		
		// Account for the last four rows: Q
		if(rank==p-1) {
			for(int i=0;i<4;++i) {
//				cout << "p,r,i,count,localNumNode = " << p << delim << r << ", " << i << ", " << count << delim << localNumNode << endl ;
				numOffDiagonalNonZero[count+i]=localStart; // globalNumNode - 4 (zeros) = # boundary Nodes
				numDiagonalNonZero[count+i]=globalNumNode-localStart-4;
			}
		}
		for(int i=0; i<localNumNode; ++i) {
//			cout << "before gather: p,r,diagNonZero,offDiagNonZero[" << i << "]: " << p << "," << rank << ", " << (numDiagonalNonZero)[i] << ", " << numOffDiagonalNonZero[i]  << endl ;
		}
		int sumDiag=0, sumOffDiag=0;
		for(int i=0;i<localNumNode;++i) {sumDiag+=numDiagonalNonZero[i]; sumOffDiag+=numOffDiagonalNonZero[i]; }
//		cout << "p,r,Diag,OffDiag" << p << delim << rank << delim << sumDiag << delim << sumOffDiag << endl ;
//		cout << "Start Allocation A:creating" << endl ;
        // Allocate the matrix.
        int ierr = (*A).Create(localNumNode,globalNumNode,numDiagonalNonZero,numOffDiagonalNonZero) ;
		if (!ierr) {
//			cout << "A matrix succesfully created for rank: " << rank << endl ;
		} else {
			cout << "A matrix creation not successful, r,ierr: " << rank << delim << ierr << endl ;
		}
//		cout << "p,r,localNumNode,globalNumNode,diagNonZero,offDiagNonZero" << p << delim << rank << delim << localNumNode << delim << globalNumNode << delim << numDiagonalNonZero << delim << numOffDiagonalNonZero << endl ;

//		cout << "Start Allocation A:created" << endl ;
        // Deallocate temporary copy of matrix allocation data.
        delete [] numDiagonalNonZero ; delete [] numOffDiagonalNonZero ;

/*		
		double distance = 0.;
		int localEnd = localStart+localNumNode ;
		double valueQ[4] ;
		PetscScalar value=0.0 ;
//		cout << "localStart, localEnd, globalNumNode, sizeQ: " << localStart << ", " << localEnd << ", " << globalNumNode << ", " << (*Q).size() << endl ;
//		cout << "here ok?" << endl ;
		// Q[row,col]=Q[3*row+col], Q=[0,1,2;3,4,5;...]
		int counter=localStart ;
		int I,I1,I2,J,J1,J2;
//		cout << "Start Assembling A,r=" << rank <<  endl ;
		for(nodePtr=seq.begin();nodePtr!=seq.end();++nodePtr,++counter){	// Q[row,col]=Q[3*row+col], Q=[0,1,2;3,4,5;...]
			I=3*counter; I1=3*counter+1; I2=3*counter+2 ;
//		for(int i=localStart;i<localEnd;++i){			
			for(int j=0;j<globalNumNode;++j){
				J=3*j; J1=3*j+1; J2=3*j+2 ;
//				cout << "I,j " << I << ", " << j << ": " ;
				distance = pow((*Q)[I] - (*Q)[J], 2.) ; //)*((*Q)[I] - (*Q)[J]) ;
				distance += pow((*Q)[I1] - (*Q)[J1], 2.) ; //*((*Q)[I1] - (*Q)[J1]) ;
				distance += pow((*Q)[I2] - (*Q)[J2], 2.) ; //)*((*Q)[I2] - (*Q)[J2]) ;
if (distance < 0.0) cout << "distance < 0" << endl ;
				distance = sqrt(distance) ;
				if (!( (*rbfNr > 0 && *rbfNr < 9) && (distance/(*r) > 1.))) { // compact
					value = radialBasisFunction(distance, *r, *a, *rbfNr) ;				
//					cout << "distance,value: " << distance << ", " << value << endl ;
					columnIndex=j; columnValue=value;
//					cout << "p,r,columnIndex,columnValue: " << p << delim << rank << delim << columnIndex << ", " << columnValue << endl ;
					(*A).SetRowValues(counter,1,&columnIndex,&columnValue) ; // i:row, j:column
				}
			}			
//			cout << "after inserting phi_bb" << endl ;
			valueQ[0] = 1.0; valueQ[1]=(*Q)[I]; valueQ[2]=(*Q)[I1]; valueQ[3]=(*Q)[I2];
			for(int j=0; j<4; ++j) {
				columnIndex=globalNumNode-4+j; columnValue=valueQ[j];
				(*A).SetRowValues(counter,1,&columnIndex,&columnValue) ; // i:row, j:column
				columnIndex=counter; 
				(*A).SetRowValues(globalNumNode-4+j,1,&columnIndex, &columnValue); // Q^T
			}			
		}
//		cout << "Start Assembling A: petsc assemble start, p,r: " << p << delim << rank << endl ;
		ierr = (*A).GetInfo() ; cout << "ierr = " << rank << delim << ierr << endl ;
		cout << "Start Assembling A,r=" << rank <<  endl ;
//		MPI_Barrier(MPI_COMM_WORLD) ;
		ierr = (*A).AssemblyBegin() ; cout << "ierr = " << rank << delim << ierr << endl ;
		ierr = (*A).AssemblyEnd() ; cout << "ierr = " << rank << delim << ierr << endl ;
		cout << "After Assembling A: petsc assemble end, p,r: " << p << delim << rank << endl ;
		*/
      }
  } ;

  register_rule<PETSCBoundaryNodeSetupMatrixUnit> registerPETSCBoundaryNodeSetupMatrixUnit ;


  // Assemble PETSC-A- system.
  class PETSCBoundaryNodeSetupMatrixApply : public apply_rule<blackbox<PETSCMatrix>,Loci::NullOp<PETSCMatrix> > {
    private:
	  const_param<real> r,a ;
	  const_param<int> rbfNr ;	  
	  const_store<vect3d> pos ;
      const_store<int> rbfBoundaryNodeNumDiagonalNonZero ;
      const_store<int> rbfBoundaryNodeNumOffDiagonalNonZero ;
      const_param<vector<int> > rbfNumBoundaryNode ;
      const_param<vector<real> > Q ;
      blackbox<PETSCMatrix> A ;
	private:
      int columnIndex,rowIndex ;
      PetscScalar columnValue ;
    public:

      // Define input and output.
      PETSCBoundaryNodeSetupMatrixApply() {
        name_store("rbfBoundaryNodeNumDiagonalNonZero",rbfBoundaryNodeNumDiagonalNonZero) ;
        name_store("rbfBoundaryNodeNumOffDiagonalNonZero",rbfBoundaryNodeNumOffDiagonalNonZero) ;
        name_store("rbfNumBoundaryNode",rbfNumBoundaryNode) ;
        name_store("petscBoundaryNodeA",A) ;		
		name_store("rbfBoundaryNodeQ",Q) ;
		name_store("rbfR",r) ;
		name_store("rbfA",a) ;
		name_store("rbfNumber",rbfNr) ;
		name_store("pos",pos) ;
		input("rbfNumber, rbfR, rbfA");
//        input("rbfBoundaryNodeToRow") ;
        input("rbfBoundaryNodeQ") ;
		input("pos");
        input("rbfBoundaryNodeNumDiagonalNonZero,rbfBoundaryNodeNumOffDiagonalNonZero") ;
        input("rbfNumBoundaryNode") ;
        output("petscBoundaryNodeA") ;
        constraint("sStar_PETSCLinearSolver");
		constraint("boundaryDisplacement") ;
        disable_threading() ;
      }

      // Assemble and solve.
	  virtual void compute(const sequence &seq) { 
	  
		int rank=Loci::MPI_rank ;
		int p=Loci::MPI_processes ;

        // Get the number of local and global nodes.
//        int localNumNode=(*rbfNumBoundaryNode)[rank],globalNumNode=0 ;
        int globalNumNode = 0;
        for(unsigned int i=0;i<(*rbfNumBoundaryNode).size();++i) globalNumNode+=(*rbfNumBoundaryNode)[i] ;
		  		
		// include the polynomial coeff. computation
//		if (rank == p-1) { // last rank takes care of beta
//			localNumNode += 4; // d+1 = 4
//		}		
		//globalNumNode+=4; // for beta (d+1)

		// Compute the row offset for this process.
        int localStart=0 ; 
        for(int i=0;i<rank;++i) localStart+=(*rbfNumBoundaryNode)[i]  ;  // 3 coordinates

		double distance = 0.;
//		int localEnd = localStart+localNumNode ;
		double valueQ[4] ;
		int ierr=0;
		PetscScalar value=0.0 ;
//		cout << "localStart, localEnd, globalNumNode, sizeQ: " << localStart << ", " << localEnd << ", " << globalNumNode << ", " << (*Q).size() << endl ;
//		cout << "here ok?" << endl ;
		// Q[row,col]=Q[3*row+col], Q=[0,1,2;3,4,5;...]
		int counter=localStart ;
		int I,I1,I2,J,J1,J2;
//		cout << "Start Assembling A inside the rule,r=" << rank << "RBF: r,a,rnfNr = " << *r << delim << *a << delim << *rbfNr << " localStart,localEnd,Qsize: " << localStart << delim << localEnd << delim << (*Q).size() << endl ;
		sequence::const_iterator nodePtr ;
		for(nodePtr=seq.begin();nodePtr!=seq.end();++nodePtr,++counter){	// Q[row,col]=Q[3*row+col], Q=[0,1,2;3,4,5;...]
			I=3*counter; I1=3*counter+1; I2=3*counter+2 ;
//		for(int i=localStart;i<localEnd;++i){			
			for(int j=0;j<globalNumNode;++j){
				J=3*j; J1=3*j+1; J2=3*j+2 ;
//				cout << "r, I,J " << rank << delim << I << ", " << J << ": " ;
				distance = pow((*Q)[I] - (*Q)[J], 2.) ; //)*((*Q)[I] - (*Q)[J]) ;
				distance += pow((*Q)[I1] - (*Q)[J1], 2.) ; //*((*Q)[I1] - (*Q)[J1]) ;
				distance += pow((*Q)[I2] - (*Q)[J2], 2.) ; //)*((*Q)[I2] - (*Q)[J2]) ;
				if (distance < 0.0) cout << "distance < 0" << endl ;
				distance = sqrt(distance) ;
				if (!( ((*rbfNr) > 0 && (*rbfNr) < 9) && (distance/(*r) > 1.))) { // compact
					value = radialBasisFunction(distance, *r, *a, *rbfNr) ;				
//					cout << "r, distance,value: " << rank << delim <<  distance << ", " << value << endl ;
					rowIndex=counter; columnIndex=j; columnValue=value;
					ierr=(*A).SetRowValues(rowIndex,1,&columnIndex,&columnValue) ; // i:row, j:column
					if (ierr) cout << "r,ierr,rowIndex,columnIndex,columnValue,sizeQ: " << rank << delim << ierr << delim << rowIndex << delim << columnIndex << delim << columnValue << delim << (*Q).size() << endl ;
				}
			}			
//			cout << "after inserting phi_bb, r, counter" << rank << delim << counter << endl ;
			valueQ[0] = 1.0; valueQ[1]=(*Q)[I]; valueQ[2]=(*Q)[I1]; valueQ[3]=(*Q)[I2];
			for(int j=0; j<4; ++j) {
//				rowIndex=counter; columnIndex=globalNumNode-4+j; columnValue=valueQ[j] ;
				rowIndex=counter; columnIndex=globalNumNode+j; columnValue=valueQ[j] ;
				ierr = (*A).SetRowValues(rowIndex,1,&columnIndex,&columnValue) ; // i:row, j:column
//				rowIndex=globalNumNode-4+j; columnIndex=counter; 
				rowIndex=globalNumNode+j; columnIndex=counter; 
				(*A).SetRowValues(rowIndex,1,&columnIndex, &columnValue); // Q^T
			}			
		}
//		cout << "Start Assembling A: petsc assemble start, p,r: " << p << delim << rank << endl ;
//		MPI_Barrier(MPI_COMM_WORLD) ; cout << "After MPI_Barrier, r: " << rank << delim << ierr << endl ;
		ierr = (*A).GetInfo() ; //cout << "ierr = " << rank << delim << ierr << endl ;
//		cout << "Start Assembling A,r=" << rank <<  endl ;
		ierr = (*A).AssemblyBegin() ; //cout << "ierr = " << rank << delim << ierr << endl ;
		ierr = (*A).AssemblyEnd() ; //cout << "ierr = " << rank << delim << ierr << endl ;
//		cout << "After Assembling A: petsc assemble end, p,r: " << p << delim << rank << endl ;
	  
	  }
  } ;

  register_rule<PETSCBoundaryNodeSetupMatrixApply> registerPETSCBoundaryNodeSetupMatrixApply ;

//---------------------------------------------solution-------------------------------------------------
//
  
  // Sets up the PETSC solution vector.
  class PETSCBoundaryNodeSetupPhi : public unit_rule {
    private:
      const_param<vector<int> > rbfNumBoundaryNode ;
      const_blackbox<PETSCMultiVector> b ;
      const_blackbox<PETSCMatrix> A ;
      const_blackbox<PETSCKsp> ksp ;
      blackbox<PETSCMultiVector> phi;
    public:

      // Define input and output.
      PETSCBoundaryNodeSetupPhi() {
        name_store("rbfNumBoundaryNode",rbfNumBoundaryNode) ;
        name_store("petscBoundaryNodePhi",phi) ;
        name_store("petscBoundaryNodeA",A) ;
        name_store("petscBoundaryNodeB",b) ;
        name_store("petscBoundaryNodeKSP",ksp) ;
        input("petscBoundaryNodeA,petscBoundaryNodeB,petscBoundaryNodeKSP") ;
        input("rbfNumBoundaryNode") ;
        output("petscBoundaryNodePhi") ;
        constraint("sStar_PETSCLinearSolver");
        constraint("boundaryDisplacement");
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {

        // Get the number of local and global nodes.
        int localNumNode=(*rbfNumBoundaryNode)[Loci::MPI_rank],globalNumNode=0 ;
        for(unsigned int i=0;i<(*rbfNumBoundaryNode).size();++i) globalNumNode+=(*rbfNumBoundaryNode)[i] ;
		
		// include the polynomial coeff. computation
		if (Loci::MPI_rank == Loci::MPI_processes-1) { // last rank takes care of beta
			localNumNode += 4; // d+1 = 4
		}
					
		globalNumNode+=4; // for beta (d+1)
		
        // Allocate the unknown and rhs vectors.
		(*phi).Create(localNumNode,globalNumNode) ;
		
		double starttime, endtime;
		// Solve the linear system
		cout << "p,r: Start Solving" << Loci::MPI_processes << delim << Loci::MPI_rank << endl ;
		starttime = MPI_Wtime();
		(*ksp).SetOperators(*A) ; (*ksp).SolveVector(*b,*phi) ;
		endtime = MPI_Wtime();
		cout << "p,r: " << Loci::MPI_processes << delim << Loci::MPI_rank << delim << " solution time: " << endtime - starttime << " sec"  << endl ;
		(*ksp).GetConvergedReason() ;
		int its = (*ksp).GetIterationNumber() ;
		if (Loci::MPI_rank==0) cout << ", number of iterations: " << its << endl ;
      }
  } ;
 
  register_rule<PETSCBoundaryNodeSetupPhi> registerPETSCBoundaryNodeSetupPhi ;


  // Assemble and solve the PETSC system.
  class PETSCBoundaryNodeSolveApply : public apply_rule<blackbox<PETSCVector>,Loci::NullOp<PETSCVector> > {
    private:
      const_param<vector<int> > rbfNumBoundaryNode ;
      const_blackbox<PETSCMultiVector> b ;
      const_blackbox<PETSCMatrix> A ;
      const_blackbox<PETSCKsp> ksp ;
      blackbox<PETSCMultiVector> phi;
    public:

      // Define input and output.
      PETSCBoundaryNodeSolveApply() {
        name_store("rbfNumBoundaryNode",rbfNumBoundaryNode) ;
        name_store("petscBoundaryNodePhi",phi) ;
        name_store("petscBoundaryNodeA",A) ;
        name_store("petscBoundaryNodeB",b) ;
        name_store("petscBoundaryNodeKSP",ksp) ;
        input("petscBoundaryNodeA,petscBoundaryNodeB,petscBoundaryNodeKSP") ;
        input("rbfNumBoundaryNode") ;
        output("petscBoundaryNodePhi") ;
        constraint("sStar_PETSCLinearSolver");
		constraint("boundaryDisplacement") ;
        disable_threading() ;
      }

      // Assemble and solve.
	  virtual void compute(const sequence &seq) { }
  } ;

  register_rule<PETSCBoundaryNodeSolveApply> registerPETSCBoundaryNodeSolveApply ;

//-----------------------------Weight---------------------------------------------------------
//

  // Sets up the rbfWeight vector.
  class RBFBoundaryNodeAssembleWeightUnit : public unit_rule { // see interpolateFile.cc
    private:
	  param<vector<real> > rbfWeight ;
    public:

      // Define input and output.
      RBFBoundaryNodeAssembleWeightUnit() {
		name_store("rbfBoundaryNodeWeight",rbfWeight) ;
        output("rbfBoundaryNodeWeight") ;
        constraint("sStar_PETSCLinearSolver") ;
        constraint("pos") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {}
  } ;
 
  register_rule<RBFBoundaryNodeAssembleWeightUnit> registerRBFBoundaryNodeAssembleWeightUnit ;

  // Assemble Q -> unit rule to output blackbox
  class RBFBoundaryNodeAssembleWeightApply : public apply_rule<param<vector<real> >,Loci::NullOp<vector<real> > > {
    private:
      const_param<vector<int> > rbfNumBoundaryNode ;
	  const_store<vect3d> pos ;
      const_blackbox<PETSCMultiVector> phi;
	  param<vector<real> > rbfWeight ;
    public:
      // Define input and output.
      RBFBoundaryNodeAssembleWeightApply() {
        name_store("rbfNumBoundaryNode",rbfNumBoundaryNode) ;
		name_store("pos",pos) ;
        name_store("rbfBoundaryNodeWeight", rbfWeight) ;
        name_store("petscBoundaryNodePhi", phi) ;
        input("rbfNumBoundaryNode") ;
        input("pos") ;
        input("petscBoundaryNodePhi") ;
        output("rbfBoundaryNodeWeight") ;
        constraint("sStar_PETSCLinearSolver") ;
		constraint("boundaryDisplacement");
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {

        // Get the number of local and global nodes.
		const int p = Loci::MPI_processes ;
		const int r = Loci::MPI_rank ;
        int localNumNode=(*rbfNumBoundaryNode)[r],globalNumNode=0 ;
        for(unsigned int i=0;i<(*rbfNumBoundaryNode).size();++i) globalNumNode+=(*rbfNumBoundaryNode)[i] ;
		
		// Compute the row offset for this process.
        int localStart=0 ; 
        for(int i=0;i<r;++i) localStart+=(*rbfNumBoundaryNode)[i]  ;  // 3 coordinates

		// include the polynomial coeff. computation
		if (r == p-1) { // last rank takes care of beta
			localNumNode += 4; // d+1 = 4
		}
					
		globalNumNode+=4; // for beta (d+1)
					
		vector<int> recvcounts(p,0) ;
		vector<int> displs(p,0) ;
		for(int i=0;i<p;++i) {
			for(int j=0;j<i;++j) displs[i]+=(3*(*rbfNumBoundaryNode)[j]);
			recvcounts[i] = 3*((i==p-1)?((*rbfNumBoundaryNode)[i]+4):(*rbfNumBoundaryNode)[i]) ;
//			cout << "working first? " << (*rbfNumBoundaryNode)[i] << endl ;
//			cout << "working? " << ((i==p-1)?(*rbfNumBoundaryNode)[i]+4:(*rbfNumBoundaryNode)[i]) << endl ;
//			cout << "r,p,displs,recvcounts[" << i <<"] = " << r << ", " << p << ", " << displs[i] << ", " << recvcounts[i] << endl ;
		}
				
		// Allocate rbfQ
		(*rbfWeight).resize(3*globalNumNode) ; // number of nodes * 3 coordinates

//        int minRowNum,maxRowNum ;
		int counter=0;
//        (*phi).GetOwnershipRange(&minRowNum,&maxRowNum) ;
        PetscScalar *phixCopy,*phiyCopy,*phizCopy ;
        (*phi).GetArray(&phixCopy,&phiyCopy,&phizCopy) ;
		// Fill in rbfQ
//		cout << "r,localStart,LocalEnd,minRowNum,maxRowNum" << r << ", " << localStart << ", " << localStart+localNumNode << ", " << minRowNum << ", " << maxRowNum << endl ;
//		cout << "r,rbfWeightSize" << r << ", " << (*rbfWeight).size() << endl ;
		for(int i=localStart; i<localStart+localNumNode; ++i,++counter) {
			(*rbfWeight)[3*i+0] = phixCopy[counter] ;
			(*rbfWeight)[3*i+1] = phiyCopy[counter] ;
			(*rbfWeight)[3*i+2] = phizCopy[counter] ;
		}
        (*phi).RestoreArray(&phixCopy,&phiyCopy,&phizCopy) ; // according to manual inexpensive

		// Allgatherv
//		for(int i=0; i<(*rbfWeight).size(); ++i) {
//			cout << "before gather: p,r,rbfWeight[" << i << "]: " << p << "," << r << ", " << (*rbfWeight)[i] << endl ;
//		}
		MPI_Allgatherv(&(*rbfWeight)[3*localStart], 3*localNumNode, MPI_DOUBLE, &(*rbfWeight)[0], &recvcounts[0], &displs[0], MPI_DOUBLE, MPI_COMM_WORLD);
//		for(int i=0; i<(*rbfWeight).size(); ++i) {
//			cout << "after gather: p,r,rbfWeight[" << i << "]: " << p << "," << r << ", " << (*rbfWeight)[i] << endl ;
//		}
      }
  } ;

  register_rule<RBFBoundaryNodeAssembleWeightApply> registerRBFBoundaryNodeAssembleWeightApply ;
//---------------------------------------------apply-------------------------------------------------
//
//
  // Rule to apply interpolation to the internal nodes
  class PETSCInternalNodeApply : public unit_rule {
    private:
	  const_param<real> r, a ;
	  const_param<int> rbfNr ;
	  const_store<vect3d> pos ;
	  const_param<vector<real> > Q ;
	  const_param<vector<real> > rbfWeight ;
	  const_param<vector<int> > rbfNumBoundaryNode ;
	  store<vect3d> node_sStar ;
    private:
	  int globalNumNode;
    public:

      // Define input and output.
      PETSCInternalNodeApply() {
		name_store("rbfBoundaryNodeQ",Q) ;
		name_store("rbfBoundaryNodeWeight",rbfWeight) ;
		name_store("rbfR",r) ;
		name_store("rbfA",a) ;
		name_store("rbfNumber",rbfNr) ;
		name_store("pos",pos) ;		
        name_store("rbfNumBoundaryNode",rbfNumBoundaryNode) ;
        name_store("node_sStar",node_sStar) ;
		input("rbfNumber, rbfR, rbfA");
		input("pos");
		input("rbfNumBoundaryNode");
		input("rbfBoundaryNodeQ") ;
		input("rbfBoundaryNodeWeight") ;
        output("node_sStar") ;
        constraint("sStar_PETSCLinearSolver") ; //internal nodes
        constraint("pos") ; //internal nodes
        disable_threading() ;
      }
	  
	  // Loop over internal nodes and apply RBF interpolation
      void calculate(Entity node) {
        	
		// globalNumNode now without pol
		// polynomial
//		cout << "globalNumNode = " << globalNumNode << endl ;
		globalNumNode = 0 ;
        for(unsigned int i=0;i<(*rbfNumBoundaryNode).size();++i) globalNumNode+=(*rbfNumBoundaryNode)[i] ;
//		cout << "after globalNumNode = " << globalNumNode << endl ;
		double distance = 0. ;
		double *value = new double[globalNumNode];
		int I,I1,I2;
		for(int i=0;i<globalNumNode;++i){
			I=3*i; I1=3*i+1; I2=3*i+2;
			distance = pow(pos[node].x - (*Q)[I], 2.) ;
			distance += pow(pos[node].y -(*Q)[I1], 2.) ;
			distance += pow(pos[node].z - (*Q)[I2], 2.) ;
			distance = sqrt(distance);
			value[i] = radialBasisFunction(distance, *r, *a, *rbfNr) ;
		}
		for(int d=0;d<3;++d) {
			double polynom = (*rbfWeight)[3*(globalNumNode+0)+d];
			polynom+=(*rbfWeight)[3*(globalNumNode+1)+d]*pos[node].x ;
			polynom+=(*rbfWeight)[3*(globalNumNode+2)+d]*pos[node].y ;
			polynom+=(*rbfWeight)[3*(globalNumNode+3)+d]*pos[node].z ;
			if (d==0) {
				node_sStar[node].x = polynom ;
			} else if (d==1) {
				node_sStar[node].y = polynom ;
			} else if (d==2) {
				node_sStar[node].z = polynom ;
			} else {
				cout << "We shouldn't be here " << endl ;
			}
		}
			
			// RBF contribution
			for(int i=0;i<globalNumNode;++i){		
		//		I=3*i; I1=3*i+1; I2=3*i+2;
		//		distance = pow(pos[node].x - (*Q)[I], 2.) ;
		//		distance += pow(pos[node].y -(*Q)[I1], 2.) ;
		//		distance += pow(pos[node].z - (*Q)[I2], 2.) ;
		//		distance = sqrt(distance);
//				value = radialBasisFunction(distance, *r, *a, *rbfNr) ;
				node_sStar[node].x += (*rbfWeight)[3*i+0] * value[i] ;
				node_sStar[node].y += (*rbfWeight)[3*i+1] * value[i] ;
				node_sStar[node].z += (*rbfWeight)[3*i+2] * value[i] ;
			}
			delete [] value ; 
      }

      // Assemble and solve.
	  virtual void compute(const sequence &seq) {
			
        // Get the number of local and global nodes.
		const int p = Loci::MPI_processes ;
		const int rank = Loci::MPI_rank ;
        int localNumNode=(*rbfNumBoundaryNode)[rank],globalNumNode=0 ;
        for(unsigned int i=0;i<(*rbfNumBoundaryNode).size();++i) globalNumNode+=(*rbfNumBoundaryNode)[i] ;
//		cout << "p,r,localNumNode,globalNumNode: " << p << delim << rank << delim << localNumNode << delim << globalNumNode << endl ;

		double starttime, endtime;
		// Solve the linear system
		cout << "p,r: Start Applying" << Loci::MPI_processes << delim << Loci::MPI_rank << endl ;
		starttime = MPI_Wtime();
		do_loop(seq,this) ;				
		endtime = MPI_Wtime();
		cout << "p,r: " << p << delim << rank << delim << "Applied in " << endtime - starttime << " sec" << endl ;
	  }	
  } ;

  register_rule<PETSCInternalNodeApply> registerPETSCInternalNodeApply ;

  // Rule to apply interpolation to the internal nodes
  class PETSCInternalNodeApplyDummy : public apply_rule<store<vect3d>,Loci::NullOp<vect3d> > {
    private:
	  const_param<real> r, a ;
	  const_param<int> rbfNr ;
	  const_store<vect3d> pos ;
	  const_param<vector<real> > Q ;
	  const_param<vector<real> > rbfWeight ;
	  const_param<vector<int> > rbfNumBoundaryNode ;
	  store<vect3d> node_sStar ;
    private:
	  int globalNumNode;
    public:

      // Define input and output.
      PETSCInternalNodeApplyDummy() {
		name_store("rbfBoundaryNodeQ",Q) ;
		name_store("rbfBoundaryNodeWeight",rbfWeight) ;
		name_store("rbfR",r) ;
		name_store("rbfA",a) ;
		name_store("rbfNumber",rbfNr) ;
		name_store("pos",pos) ;		
        name_store("rbfNumBoundaryNode",rbfNumBoundaryNode) ;
        name_store("node_sStar",node_sStar) ;
		input("rbfNumber, rbfR, rbfA");
		input("pos");
		input("rbfNumBoundaryNode");
		input("rbfBoundaryNodeQ") ;
		input("rbfBoundaryNodeWeight") ;
        output("node_sStar") ;
        constraint("sStar_PETSCLinearSolver") ; //internal nodes
        constraint("pos") ; //internal nodes
        disable_threading() ;
      }
	  
	  // Loop over internal nodes and apply RBF interpolation
      //void calculate(Entity node) {
      //}

      // Assemble and solve.
	  virtual void compute(const sequence &seq) {
	  }	
  } ;

  register_rule<PETSCInternalNodeApplyDummy> registerPETSCInternalNodeApplyDummy ;

  // Rule to update zero to RBF-x-direction-displacement: SXZero
  class RbfXSetZeroSym : public apply_rule<store<vect3d>,Loci::NullOp<vect3d> > {
    private:
      store<vect3d> node_sStar ;
    public:

      // Define input and output.
      RbfXSetZeroSym() {
        name_store("node_sStar",node_sStar) ;
        output("node_sStar") ;
        constraint("pos") ; //internal nodes
        constraint("sXZeroNodes") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity node) { 
	  	  node_sStar[node].x=0.0;
	  }

      // Assign flag for a sequence of nodes with sym.
      virtual void compute(const sequence &seq) { 
//		  cout << "inside RbfXSetZeroSym" << endl;
		  do_loop(seq,this) ; 
	  }
  } ;

  register_rule<RbfXSetZeroSym> registerRbfXSetZeroSym ;

  // Rule to update zero to RBF-y-direction-displacement: SZZero
  class RbfYSetZeroSym : public apply_rule<store<vect3d>,Loci::NullOp<vect3d> > {
    private:
      store<vect3d> node_sStar ;
    public:

      // Define input and output.
      RbfYSetZeroSym() {
        name_store("node_sStar",node_sStar) ;
        output("node_sStar") ;
        constraint("pos") ; //internal nodes
        constraint("sYZeroNodes") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity node) { 
	  	  node_sStar[node].y=0.0;
	  }

      // Assign flag for a sequence of nodes with sym.
      virtual void compute(const sequence &seq) { 
//		  cout << "inside RbfYSetZeroSym" << endl;
		  do_loop(seq,this) ; 
	  }
  } ;

  register_rule<RbfYSetZeroSym> registerRbfYSetZeroSym ;

  // Rule to update zero to RBF-z-direction-displacement: SZZero
  class RbfZSetZeroSym : public apply_rule<store<vect3d>,Loci::NullOp<vect3d> > {
    private:
      store<vect3d> node_sStar ;
    public:

      // Define input and output.
      RbfZSetZeroSym() {
        name_store("node_sStar",node_sStar) ;
        output("node_sStar") ;
        constraint("pos") ; //internal nodes
        constraint("sZZeroNodes") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity node) { 
	  	  node_sStar[node].z=0.0;
	  }

      // Assign flag for a sequence of nodes with sym.
      virtual void compute(const sequence &seq) { 
//		  cout << "inside RbfZSetZeroSym" << endl;
		  do_loop(seq,this) ; 
	  }
  } ;

  register_rule<RbfZSetZeroSym> registerRbfZSetZeroSym ;
}


