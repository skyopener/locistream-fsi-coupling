#line 1 "FSI_CFD2CSDInterpolation.loci"
//-----------------------------------------------------------------------------
// Description: This file contains rules for implementing a FSI linear solver
//
// using FSI.
//-----------------------------------------------------------------------------

// mpi
#include<mpi.h>

// Standard library includes.
#include<vector>
#include<string>
using std::vector ;
#include<cmath>

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


// StreamUns includes.
#include "sciTypes.h"

// RBF includes.
#include "rbf.h"
        
// boost::ublas includes
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas ;
	
// Include FSI wrappers

// FSI includes.
#include <petsc.h>
#include <petscerror.h>
#include <petscksp.h>
#include <petscvec.h>

// StreamUns includes.
#include "sciTypes.h"

// delim
#include "delim.h"


// Create a namespace for the FSI objects so there is no confusion. Note that
// there is already a type "Mat" defined by Loci.
namespace Petsc {
  typedef struct _p_Vec* Vec ;
  typedef struct _p_Mat* Mat ;
  typedef struct _p_KSP* KSP ;
} ;



namespace streamUns {

	// delimiter
	//delim = ", " ;
//-----------------------------------------------------------------------------
// Wrapper classes for FSI objects. These have been created primarily to
// handle memory management.
  
  // Wrapper for Vec. Note that we cannot check the success of the destroy in
  // the destructor since the destructor cannot return a value.
  class PETSCCFD2CSDVector {
    private:
      bool created ;
      mutable Petsc::Vec v ;
    public:
      PETSCCFD2CSDVector() : created(false) {}
      ~PETSCCFD2CSDVector() { if(created) VecDestroy(v) ; }
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
      int DuplicateFrom(const PETSCCFD2CSDVector &a) {
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
  class PETSCCFD2CSDMultiVector {
    private:
      bool created ;
      mutable Petsc::Vec vX,vY,vZ ;
    public:
      PETSCCFD2CSDMultiVector() : created(false) {}
      ~PETSCCFD2CSDMultiVector() {
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
      int DuplicateFrom(const PETSCCFD2CSDMultiVector &a) {
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
      
//      int VecView() const {
//      	int ierr=VecView(vX, PETSC_VIEWER_STDOUT_WORLD) ; CHKERRQ(ierr) ;
//      	ierr=VecView(vY, PETSC_VIEWER_STDOUT_WORLD) ; CHKERRQ(ierr) ;
//      	ierr=VecView(vZ, PETSC_VIEWER_STDOUT_WORLD) ; CHKERRQ(ierr) ;
//      	return ierr ;
//      }
  } ;

  // Wrapper for Mat. This class assumes a square matrix. Note that we cannot
  // check the success of the destroy in the destructor since the destructor
  // cannot return a value.
  class PETSCCFD2CSDMatrix {
    private:
      bool created ;
      mutable Petsc::Mat m ;
    public:
      PETSCCFD2CSDMatrix() : created(false) {}
      ~PETSCCFD2CSDMatrix() { if(created) MatDestroy(m) ; }
    public:
      int AssemblyBegin() const {
        int ierr=MatAssemblyBegin(m,MAT_FINAL_ASSEMBLY) ;
        CHKERRQ(ierr) ; return 0 ;
      }
      int AssemblyEnd() const {
       int ierr=MatAssemblyEnd(m,MAT_FINAL_ASSEMBLY) ;
       CHKERRQ(ierr) ; return 0 ;
      }
      int Create(int numLocalRow,int numGlobalRow,int *numDiagonalNonZero,int
      *numOffDiagonalNonZero) {
        if(created) return 0 ;
        if(Loci::MPI_processes>1){
          int ierr=MatCreateMPIAIJ(PETSC_COMM_WORLD,numLocalRow,numLocalRow,
            numGlobalRow,numGlobalRow,0,numDiagonalNonZero,0,
            numOffDiagonalNonZero,&m) ; CHKERRQ(ierr) ;
        }else{
          int ierr=MatCreateSeqAIJ(PETSC_COMM_WORLD,numGlobalRow,numGlobalRow,
            0,numDiagonalNonZero,&m) ; CHKERRQ(ierr) ;
        }
        created=true ; return 0 ;
      }
      const Petsc::Mat& Data() const { return m ; }
      int Destroy() const {
      	if(created) MatDestroy(m) ;
      }
	  int GetInfo() {
		  MatInfo info;
		  double mal, nz_alloc, nz_used, nz_un ;
		  MatGetInfo(m,MAT_LOCAL,&info);
		  mal = info.mallocs; nz_alloc=info.nz_allocated; nz_used=info.nz_used; nz_un=info.nz_unneeded ;
		  cout << "Force Interpolation CFD2CSD: rank = " << Loci::MPI_rank << ", Number of mallocs during MatSetValues(): " << mal << delim << "Number of nonzeros (allocated, used, unneeded): " << nz_alloc << delim << nz_used << delim << nz_un << endl ;
		  return 0 ;
	  }
      int GetOwnershipRange(int *firstRow,int *lastRow) {
        int ierr=MatGetOwnershipRange(m,firstRow,lastRow) ; CHKERRQ(ierr) ;
        return 0 ;
      }
      int SetRowValues(int rowNum,int numColumn,int *columnIndices,PetscScalar
      *columnValues) const {
        int ierr=MatSetValues(m,1,&rowNum,numColumn,columnIndices,columnValues,
          INSERT_VALUES) ; CHKERRQ(ierr) ; return 0 ;
      }
  } ;
  
  // Wrapper for KSP.
  class PETSCCFD2CSDKsp {
    private:
      bool created ;
      mutable Petsc::KSP ksp ;
    public:
      PETSCCFD2CSDKsp() : created(false) {}
      ~PETSCCFD2CSDKsp() { if(created) KSPDestroy(ksp) ; }
    public:
      int Create() {
        if(created) return 0 ;
        int ierr=KSPCreate(PETSC_COMM_WORLD,&ksp) ; CHKERRQ(ierr) ;
        created=true ; return 0 ;
      }
      int Destroy() {
      	if(created) KSPDestroy(ksp) ;
      }
      	  int GetConvergedReason() const {
		  KSPConvergedReason reason ;
		  KSPGetConvergedReason(ksp, &reason) ;
		  if (Loci::MPI_rank==0) {
			  if (reason==KSP_CONVERGED_RTOL) {
				  cout << "CFD2CSD: PETSC solver: GetConvergedReason = " << reason << " KSP_CONVERGED_RTOL" << endl ;
			  } else if (reason==KSP_CONVERGED_ATOL) {
				  cout << "CFD2CSD: PETSC solver: GetConvergedReason = " << reason << " KSP_CONVERGED_ATOL" << endl ;
			  } else if (reason==KSP_CONVERGED_ITS) {
				  cout << "CFD2CSD: PETSC solver: GetConvergedReason = " << reason << " KSP_CONVERGED_ITS" << endl ;
			  } else if (reason==KSP_DIVERGED_NULL) {
				  cout << "CFD2CSD: PETSC solver: GetConvergedReason = " << reason << " KSP_DIVERGED_NULL" << endl ;
			  } else if (reason==KSP_DIVERGED_ITS) {
				  cout << "CFD2CSD: PETSC solver: GetConvergedReason = " << reason << " KSP_DIVERGED_ITS" << endl ;
			  } else if (reason==KSP_DIVERGED_DTOL) {
				  cout << "CFD2CSD: PETSC solver: GetConvergedReason = " << reason << " KSP_DIVERGED_DTOL" << endl ;
			  } else {
				  cout << "CFD2CSD: PETSC solver: GetConvergedReason = " << reason << endl ;
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
      int SetOperators(const PETSCCFD2CSDMatrix &m) const {
        int ierr=KSPSetOperators(ksp,m.Data(),m.Data(),DIFFERENT_NONZERO_PATTERN) ; // SAME_NONZERO_PATTERN
        CHKERRQ(ierr) ; return 0 ;
      }
      void SetTolerances(real relativeTolerance,real absoluteTolerance,int maxIterations) const {
        KSPSetTolerances(ksp,relativeTolerance,absoluteTolerance,PETSC_DEFAULT,maxIterations) ;
      }
      int SetRestart(int Restart) {
		  int ierr=KSPGMRESSetRestart(ksp, Restart) ;
		  CHKERRQ(ierr) ; return ierr ;
	  }
      int Solve(const PETSCCFD2CSDVector &rhs, PETSCCFD2CSDVector &sol) const {
//      int numIteration ;
        int ierr=KSPSolve(ksp,rhs.Data(),sol.Data()) ; CHKERRQ(ierr) ;
//KSPGetIterationNumber(ksp,&numIteration) ;
//cout << "fsi iterations: " << numIteration << endl ;
        return ierr ;
      }
      int SolveVector(const PETSCCFD2CSDMultiVector &rhs, PETSCCFD2CSDMultiVector &sol)
      const {
//      int numIteration ;
		//		cout << "Inside solve vector 1 " << Loci::MPI_processes << delim << Loci::MPI_rank << endl ;
        int ierr=KSPSolve(ksp,rhs.DataX(),sol.DataX()) ; CHKERRQ(ierr) ;
  //      cout << "Inside solve vector 2 " << Loci::MPI_processes << delim << Loci::MPI_rank << endl ;
        ierr=KSPSolve(ksp,rhs.DataY(),sol.DataY()) ; CHKERRQ(ierr) ;
   //      cout << "Inside solve vector 3 " << Loci::MPI_processes << delim << Loci::MPI_rank << endl ;
        ierr=KSPSolve(ksp,rhs.DataZ(),sol.DataZ()) ; CHKERRQ(ierr) ;
        //KSPGetIterationNumber(ksp,&numIteration) ;
//cout << "fsi iterations: " << numIteration << endl ;
        return ierr ;
      }
  } ;

 


//---------------------------------------------b-------------------------------------------------
//
//
  // Sets up the FSI right-hand-side vector.
  //class FSIBoundaryFaceUnitRHS : public unit_rule {
  class FSIBoundaryFaceUnitRHS : public singleton_rule {  	
    private:
      const_param<vector<int> > fsiNumBoundaryFace ;
      const_store<vect3d> B ;
	  	const_store<vect3d> faceCenter ;
      const_param<bool> CFDIterationFinished ;
      blackbox<PETSCCFD2CSDMultiVector> b;
    public:

      // Define input and output.
      FSIBoundaryFaceUnitRHS() {
        name_store("fsiNumBoundaryFace(X){n,it}",fsiNumBoundaryFace) ;
        name_store("fsiBoundaryFaceB(X){n,it}",b) ;
        name_store("FSI_B(X){n,it}",B);
        name_store("CFDIterationFinished{n,it-1}",CFDIterationFinished) ;
   //     name_store("facecenter{n,it-1}",faceCenter) ;
  //      input("facecenter{n,it-1}") ;
    		input("CFDIterationFinished{n,it-1}") ;
        input("fsiNumBoundaryFace(X){n,it}") ;
	input("FSI_B(X){n,it}");
        output("fsiBoundaryFaceB(X){n,it}") ;
 //       constraint("fsi_FSILinearSolver") ;
	constraint("X{n,it}") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {
      	
	
	
      	if (*CFDIterationFinished) { 
      		
      	int rank=	Loci::MPI_rank;
      	int p=Loci::MPI_processes;

        // Get the number of local and global nodes.
        int localNumFace=(*fsiNumBoundaryFace)[rank],globalNumFace=0 ;
        for(unsigned int i=0;i<(*fsiNumBoundaryFace).size();++i) globalNumFace+=(*fsiNumBoundaryFace)[i] ;
		
		// include the polynomial coeff. computation
				if (rank == p-1) { // last rank takes care of beta
					localNumFace += 4; // d+1 = 4
				}
					
				globalNumFace+=4; // for beta (d+1)
		
        // Allocate the unknown and rhs vectors.
				(*b).Create(localNumFace,globalNumFace) ;

		// Compute the row offset for this process.
        int localStart=0 ;
        for(int i=0;i<rank;++i){ localStart+=(*fsiNumBoundaryFace)[i] ; }
		//
					
				sequence::const_iterator facePtr;
				PetscScalar valuex, valuey, valuez;

	//	cout << "rank, localStart, localNumFace, globalNumFace: " << Loci::MPI_rank << ", " << localStart << ", " << localNumFace << ", " << globalNumFace << endl;
		
				int counter=localStart ;
				for(facePtr=seq.begin();facePtr!=seq.end();++facePtr,++counter){			
						// RHS vector
						valuex=B[*facePtr].x;
						valuey=B[*facePtr].y;
						valuez=B[*facePtr].z;
	//					cout << "rank, i,B.x,valuex,B.y,valuey,B.z,valuez" << Loci::MPI_rank << ", " << counter << ", " << valuex << delim << valuey << delim << valuez << endl;
						(*b).SetValue(&counter,&valuex, &valuey, &valuez) ;
				}		

		// assemble
				int rowIndex=0;
				if (rank==p-1) { // the last four rows of B are zero
					for(int i=0; i<4; ++i) {
							rowIndex=globalNumFace-4+i; valuex=0.0; valuey=0.0; valuez=0.0;
						(*b).SetValue(&rowIndex,&valuex, &valuey, &valuez) ;
	//					cout << "rank = " << rank << ", localNumFace, globalNumFace, rowIndex, valuex, valuey, valuez: " << localNumFace << delim << globalNumFace << delim << rowIndex << delim << valuex << delim << valuey << delim << valuez << endl ;
					}
				}
				(*b).AssemblyBegin(); (*b).AssemblyEnd() ;
				
	//			if (rank==0) cout << "petsc vecview" << endl ;
				//(*b).VecView() ;
	//			VecView((*b).DataX(), PETSC_VIEWER_STDOUT_WORLD);
	//			VecView((*b).DataY(), PETSC_VIEWER_STDOUT_WORLD);
	//			VecView((*b).DataZ(), PETSC_VIEWER_STDOUT_WORLD);
		 } // CFDIterationFinished
		}
  } ;
 
  register_rule<FSIBoundaryFaceUnitRHS> registerFSIBoundaryFaceUnitRHS ;


//---------------------------------------------A-------------------------------------------------
//

//$type fsiBoundaryFaceA(X) blackbox<PETSCCFD2CSDMatrix>  ;
// $type fsiBoundaryFacePhi(X0) blackbox<PETSCCFD2CSDMultiVector> 
// $type fsiBoundaryFaceB(X0) blackbox<PETSCCFD2CSDMultiVector> 
//$type fsiBoundaryFaceNumDiagonalNonZero(X) store<int> ; 
//$type fsiBoundaryFaceNumOffDiagonalNonZero(X) store<int> ;
// $type fsiNumBoundaryFace(X0) param<vector<int> > 
// $type fsiBoundaryFaceQ(X0) blackbox<vector<real> > 
// $type FSIRBFr param<real> 
// $type FSIRBFa param<real> 
// $type FSIRBFnr param<int> 
// $type CFDIterationFinished param<bool> 
// $type FSIRBFMaxLinearSolverIterations param<int> 
// $type FSIIterationTolerance param<real> 
// $type fsiBoundaryFaceWeight(X0) blackbox<vector<real> > 

namespace {class file_FSI_CFD2CSDInterpolation000_1280694170m614 : public Loci::blackbox_rule {
#line 463 "FSI_CFD2CSDInterpolation.loci"
    Loci::const_param<real>  L_FSIRBFr_ ; 
#line 463 "FSI_CFD2CSDInterpolation.loci"
    Loci::const_param<real>  L_FSIRBFa_ ; 
#line 463 "FSI_CFD2CSDInterpolation.loci"
    Loci::const_param<int>  L_FSIRBFnr_ ; 
#line 463 "FSI_CFD2CSDInterpolation.loci"
    Loci::const_param<int>  L_FSIRBFMaxLinearSolverIterations_ ; 
#line 463 "FSI_CFD2CSDInterpolation.loci"
    Loci::const_param<real>  L_FSIIterationTolerance_ ; 
#line 463 "FSI_CFD2CSDInterpolation.loci"
    Loci::const_param<vector<int> >  L_fsiNumBoundaryFaceX_nit__ ; 
#line 463 "FSI_CFD2CSDInterpolation.loci"
    Loci::const_blackbox<vector<real> >  L_fsiBoundaryFaceQX_nit__ ; 
#line 463 "FSI_CFD2CSDInterpolation.loci"
    Loci::const_blackbox<PETSCCFD2CSDMultiVector>  L_fsiBoundaryFaceBX_nit__ ; 
#line 463 "FSI_CFD2CSDInterpolation.loci"
    Loci::const_param<bool>  L_CFDIterationFinished_nit_M_1__ ; 
#line 463 "FSI_CFD2CSDInterpolation.loci"
    Loci::blackbox<vector<real> >  L_fsiBoundaryFaceWeightX_nit__ ; 
#line 463 "FSI_CFD2CSDInterpolation.loci"
public:
#line 463 "FSI_CFD2CSDInterpolation.loci"
    file_FSI_CFD2CSDInterpolation000_1280694170m614() {
#line 463 "FSI_CFD2CSDInterpolation.loci"
       name_store("FSIRBFr",L_FSIRBFr_) ;
#line 463 "FSI_CFD2CSDInterpolation.loci"
       name_store("FSIRBFa",L_FSIRBFa_) ;
#line 463 "FSI_CFD2CSDInterpolation.loci"
       name_store("FSIRBFnr",L_FSIRBFnr_) ;
#line 463 "FSI_CFD2CSDInterpolation.loci"
       name_store("FSIRBFMaxLinearSolverIterations",L_FSIRBFMaxLinearSolverIterations_) ;
#line 463 "FSI_CFD2CSDInterpolation.loci"
       name_store("FSIIterationTolerance",L_FSIIterationTolerance_) ;
#line 463 "FSI_CFD2CSDInterpolation.loci"
       name_store("fsiNumBoundaryFace(X){n,it}",L_fsiNumBoundaryFaceX_nit__) ;
#line 463 "FSI_CFD2CSDInterpolation.loci"
       name_store("fsiBoundaryFaceQ(X){n,it}",L_fsiBoundaryFaceQX_nit__) ;
#line 463 "FSI_CFD2CSDInterpolation.loci"
       name_store("fsiBoundaryFaceB(X){n,it}",L_fsiBoundaryFaceBX_nit__) ;
#line 463 "FSI_CFD2CSDInterpolation.loci"
       name_store("CFDIterationFinished{n,it-1}",L_CFDIterationFinished_nit_M_1__) ;
#line 463 "FSI_CFD2CSDInterpolation.loci"
       name_store("fsiBoundaryFaceWeight(X){n,it}",L_fsiBoundaryFaceWeightX_nit__) ;
#line 463 "FSI_CFD2CSDInterpolation.loci"
       input("fsiNumBoundaryFace(X){n,it},fsiBoundaryFaceQ(X){n,it},fsiBoundaryFaceB(X){n,it},CFDIterationFinished{n,it-1},FSIRBFr,FSIRBFa,FSIRBFnr,FSIRBFMaxLinearSolverIterations,FSIIterationTolerance") ;
#line 463 "FSI_CFD2CSDInterpolation.loci"
       output("fsiBoundaryFaceWeight(X){n,it}") ;
#line 463 "FSI_CFD2CSDInterpolation.loci"
       constraint("X{n,it}") ;
#line 463 "FSI_CFD2CSDInterpolation.loci"
       disable_threading() ;
#line 463 "FSI_CFD2CSDInterpolation.loci"
    }
#line 463 "FSI_CFD2CSDInterpolation.loci"
    void prelude(const Loci::sequence &seq) { 
				
  
  int rank=Loci::MPI_rank ;
  int p=Loci::MPI_processes ;
  real r=(*L_FSIRBFr_), a=(*L_FSIRBFa_) ;
  int fsiNr = (*L_FSIRBFnr_) ;	
  int I,I1,I2,J,J1,J2;
  real distance ;
  real linearSolverTolerance = (*L_FSIIterationTolerance_);
  int linearSolverMaxIteration = (*L_FSIRBFMaxLinearSolverIterations_) ;

  if (Loci::MPI_rank==0) cout << "fsiBoundaryFaceWeight<-fsiNumBoundaryFace CFDIterationFinished" << (*L_CFDIterationFinished_nit_M_1__) << endl ;
  
  if (*L_CFDIterationFinished_nit_M_1__) {	
		//		cout << "p, rank, r, a, fsiNr: "<< p << delim << rank << delim << r  << delim << a << delim << fsiNr << endl ;

    // Get the number of local and global nodes.
        int localNumFace=(*L_fsiNumBoundaryFaceX_nit__)[rank],globalNumFace=0 ;
        for(unsigned int i=0;i<(*L_fsiNumBoundaryFaceX_nit__).size();++i) globalNumFace+=(*L_fsiNumBoundaryFaceX_nit__)[i] ;
		  		  		
    //    cout << "p, r, localNumface, globalnumFace: "<< p << delim << rank << delim << localNumFace  << delim << globalNumFace << endl ;


		// Compute the row offset for this process.
        int localStart=0 ; 
        for(int i=0;i<rank;++i) localStart+=(*L_fsiNumBoundaryFaceX_nit__)[i]  ;  // 3 coordinates

    //    cout << "p, r, localstart: "<< p << delim << rank << delim << localStart  << endl ;

				
				vector<int> numberOfDiagNonZeros(localNumFace+((rank==p-1)?4:0),0) ;
				vector<int> numberOfOffDiagNonZeros(localNumFace+((rank==p-1)?4:0),0) ;
//				vector<int> numberOfOffDiagNonZeros(localNumFace+((rank==p-1)?4:0),0) ;
					
				int localEnd=localStart+localNumFace ;
				//if (rank == p-1) localEnd-=4 ;
				
		//		cout << "Start Allocation A: size of numberOfDiagNonZeros ->, rank = " << rank << ", " << (numberOfDiagNonZeros).size() << endl ;
		//		cout << "Start Allocation A: size of numberOfOffDiagNonZeros ->, rank = " << rank << ", " << (numberOfOffDiagNonZeros).size() << endl ;
				
				
		// Determine the number of nonzero local elements		
		
		    cout << "Start Allocation A: rank = " << rank << endl ; //", " << (*$fsiBoundaryFaceQ(X){n,it}{n}).size() << endl ;
		  //  for(int i=localStart;i<localEnd;++i) {
		 //   	I=3*i; I1=3*i+1; I2=3*i+2 ;
		  //  	cout << "rank: " << rank << ", i: " << i << "Qx, Qy, Qz: " << (*$fsiBoundaryFaceQ(X){n,it}{n})[I] << delim << (*$fsiBoundaryFaceQ(X){n,it}{n})[I1] << delim << (*$fsiBoundaryFaceQ(X){n,it}{n})[I2] << endl ;
		  //  }
		    
		
				for(int i=localStart;i<localEnd;++i) { // Looping over Q
					I=3*i; I1=3*i+1; I2=3*i+2 ;
					int counterD=0; int counterOD=0;
					for(int j=0;j<globalNumFace;++j){		
						J=3*j; J1=3*j+1; J2=3*j+2 ;				
						if( (j>=localStart) && (j<localEnd) ) {		// DIAGONAL					
							if (fsiNr > 0 && fsiNr < 9) { // compact								
								distance = pow((*L_fsiBoundaryFaceQX_nit__)[I] - (*L_fsiBoundaryFaceQX_nit__)[J], 2.) ; //)*((*Q)[I] - (*Q)[J]) ;
								distance += pow((*L_fsiBoundaryFaceQX_nit__)[I1] - (*L_fsiBoundaryFaceQX_nit__)[J1], 2.) ; //*((*Q)[I1] - (*Q)[J1]) ;
								distance += pow((*L_fsiBoundaryFaceQX_nit__)[I2] - (*L_fsiBoundaryFaceQX_nit__)[J2], 2.) ; //)*((*Q)[I2] - (*Q)[J2]) ;
								distance = sqrt(distance) / r ;
								if (distance <=1.) {
									++counterD;
									//numberOfDiagNonZeros[i]+=1;
								} 
							} else { // global
								++counterD;
								//numberOfDiagNonZeros[i]+=1;
							}
						} else {																	// OFFDIAGONAL		
							if (fsiNr > 0 && fsiNr < 9) { // compact
								distance = pow((*L_fsiBoundaryFaceQX_nit__)[I] - (*L_fsiBoundaryFaceQX_nit__)[J], 2.) ; //)*((*Q)[I] - (*Q)[J]) ;
								distance += pow((*L_fsiBoundaryFaceQX_nit__)[I1] - (*L_fsiBoundaryFaceQX_nit__)[J1], 2.) ; //*((*Q)[I1] - (*Q)[J1]) ;
								distance += pow((*L_fsiBoundaryFaceQX_nit__)[I2] - (*L_fsiBoundaryFaceQX_nit__)[J2], 2.) ; //)*((*Q)[I2] - (*Q)[J2]) ;
								distance = sqrt(distance) / r ;
								if (distance <=1.) {
									++counterOD ;
									//numberOfOffDiagNonZeros[i]+=1;
								}
							} else { // global
								++counterOD ;
								//numberOfOffDiagNonZeros[i]+=1;
							}
						}
					}
					numberOfDiagNonZeros[i-localStart]=counterD;
					numberOfOffDiagNonZeros[i-localStart]=counterOD;
			//		cout << "rank, i, counterD, counterOD = " << rank << delim << i << delim << counterD << delim << counterOD << endl ;
				}
				
//				for(int i=0;i<numberOfDiagNonZeros.size();++i) {
//						//sumDiag+=numberOfDiagNonZeros[i]; sumOffDiag+=numberOfOffDiagNonZeros[i]; 
//						cout << "r, i, numOfDiagNonZero, numOfOFfdiag = " << rank << delim << i << delim << numberOfDiagNonZeros[i] << delim << numberOfOffDiagNonZeros[i] << endl ;						
//				}
				
				int count = 0 ; int tempInt ;
				for(count=0; count<localNumFace;++count){
					tempInt = numberOfDiagNonZeros[count] ;
          numberOfDiagNonZeros[count] = tempInt + ((rank==p-1)?4:0) ; // + Q, only r==p1 gets Q in diagNonZero, otherwise in offDiagNonZero
          tempInt = numberOfOffDiagNonZeros[count] ;	
          numberOfOffDiagNonZeros[count] = tempInt + ((rank==p-1)?0:4) ;
//          numberOfOffDiagNonZeros[count] += ((rank==p-1)?4:0);
        }					
							
		//	    cout << "p, r, count is " << p << delim << rank << delim << count << endl;
			    		
				// Account for the last four rows: Q
				if(rank==p-1) {
					for(int i=0;i<4;++i) {
						numberOfOffDiagNonZeros[localNumFace+i]=localStart; // globalNumFace - 4 (zeros) = # boundary Faces
						numberOfDiagNonZeros[localNumFace+i]=globalNumFace+4-localStart;
//						numberOfDiagNonZeros[count+i]=globalNumFace-localStart-4;
					}
				}	
				
				// include the polynomial coeff. computation
				if (rank == p-1) { // last rank takes care of beta
					localNumFace += 4; // d+1 = 4
				}		
				globalNumFace+=4; // for beta (d+1)
				//localEnd=localStart+localNumFace ;// for beta (d+1), localNumFace updated for rank==p-1
				
				
    		int sumDiag=0, sumOffDiag=0;
				for(int i=0;i<localNumFace;++i) {
						sumDiag+=numberOfDiagNonZeros[i]; sumOffDiag+=numberOfOffDiagNonZeros[i]; 
		//				cout << "r, i, numOfDiagNonZero, numOfOFfdiag = " << rank << delim << i << delim << numberOfDiagNonZeros[i] << delim << numberOfOffDiagNonZeros[i] << endl ;						
				}
		//		cout << "Force Interpolation CFD2CSD: p,r,Diag,OffDiag,localNumFace,globalNumFace " << p << delim << rank << delim << sumDiag << delim << sumOffDiag << delim << localNumFace << delim << globalNumFace << endl ;
				
				PETSCCFD2CSDMatrix CFD2CSD_A_matrix ;
        // Allocate the matrix.
        int ierr = (CFD2CSD_A_matrix).Create(localNumFace,globalNumFace,&numberOfDiagNonZeros[0],&numberOfOffDiagNonZeros[0]) ;
				if (!ierr) {
					cout << "Force Interpolation: A matrix succesfully created for rank, localNumNode, globalNumNode, sumDiag, sumOffDiag " << rank << delim << localNumFace << delim << globalNumFace << delim << sumDiag << delim << sumOffDiag << endl ;
				} else {
					cout << "A matrix creation not successful, r,ierr: " << rank << delim << ierr << endl ;
				}
        // Deallocate temporary copy of matrix allocation data.
        //delete [] numDiagonalNonZero ; delete [] numOffDiagonalNonZero ;
		
				cout << "Start Assembling A, rank" << rank << endl ;
				int columnIndex, rowIndex ;
        PetscScalar columnValue ;
				//double distance = 0.;
				//int localEnd = localStart+localNumFace ;
				double valueQ[4] ;
				PetscScalar value=0.0 ;
				//int counter=localStart ;
				
				
				
				for(int counter=localStart; counter<localEnd;++counter){
				//for(nodePtr=seq.begin();nodePtr!=seq.end();++nodePtr,++counter){	// Q[row,col]=Q[3*row+col], Q=[0,1,2;3,4,5;...]
					I=3*counter; I1=3*counter+1; I2=3*counter+2 ;
					for(int j=0;j<globalNumFace-4;++j){
						J=3*j; J1=3*j+1; J2=3*j+2 ;
						distance = pow((*L_fsiBoundaryFaceQX_nit__)[I] - (*L_fsiBoundaryFaceQX_nit__)[J], 2.) ; //)*((*Q)[I] - (*Q)[J]) ;
						distance += pow((*L_fsiBoundaryFaceQX_nit__)[I1] - (*L_fsiBoundaryFaceQX_nit__)[J1], 2.) ; //*((*Q)[I1] - (*Q)[J1]) ;
						distance += pow((*L_fsiBoundaryFaceQX_nit__)[I2] - (*L_fsiBoundaryFaceQX_nit__)[J2], 2.) ; //)*((*Q)[I2] - (*Q)[J2]) ;
						distance = sqrt(distance) ;
						if (!( (fsiNr > 0 && fsiNr < 9) && (distance/(r) > 1.))) { // compact
							value = radialBasisFunction(distance, r, a, fsiNr) ;				
							rowIndex=counter; columnIndex=j; columnValue=value;
							(CFD2CSD_A_matrix).SetRowValues(rowIndex,1,&columnIndex,&columnValue) ; // i:row, j:column
				//			cout << "columnvalueM = " << columnValue << endl ;
						}
					}			
					valueQ[0] = 1.0; valueQ[1]=(*L_fsiBoundaryFaceQX_nit__)[I]; valueQ[2]=(*L_fsiBoundaryFaceQX_nit__)[I1]; valueQ[3]=(*L_fsiBoundaryFaceQX_nit__)[I2];
					for(int j=0; j<4; ++j) {
						rowIndex=counter; columnIndex=globalNumFace-4+j; columnValue=valueQ[j];
						(CFD2CSD_A_matrix).SetRowValues(rowIndex,1,&columnIndex,&columnValue) ; // i:row, j:column
						rowIndex=globalNumFace-4+j; columnIndex=counter; 
						(CFD2CSD_A_matrix).SetRowValues(rowIndex,1,&columnIndex, &columnValue); // Q^T
		//				cout << "columnvalueQ = " << columnValue << endl ;
					}			
				}
				
				
			(CFD2CSD_A_matrix).GetInfo() ;
			(CFD2CSD_A_matrix).AssemblyBegin() ; (CFD2CSD_A_matrix).AssemblyEnd() ;
			cout << "Force Interpolation CFD2CSD: Matrix Assembled, r = " << rank << endl ;
			
			//cout << "Inside CFD2CSD phi: p,r,localNumFace,GlobalNumFace" << Loci::MPI_processes << delim << Loci::MPI_rank << delim << localNumFace << delim << globalNumFace << endl ;
				
		    // Allocate the unknown and rhs vectors.
		//     cout << "Force Interpolation CFD2CSD: p,r: Creating phi" << Loci::MPI_processes << delim << Loci::MPI_rank << endl ;   
		    
		    PETSCCFD2CSDMultiVector phi ;
		    (phi).Create(localNumFace,globalNumFace) ;
				//(*$fsiBoundaryFacePhi(X){n}).Create(localNumFace,globalNumFace) ;
	//			cout << "Force Interpolation CFD2CSD: p,r: phi created" << Loci::MPI_processes << delim << Loci::MPI_rank << endl ;
				
				PETSCCFD2CSDKsp CFD2CSDKsp_ksp ;
				CFD2CSDKsp_ksp.Create() ;
				CFD2CSDKsp_ksp.SetTolerances(linearSolverTolerance,1.0e-30,linearSolverMaxIteration) ;
				PC pc ; CFD2CSDKsp_ksp.GetPC(&pc) ; PCSetType(pc,PCJACOBI) ; // PCSetFromOptions(pc) ;
				//PC pc ; FSI_ksp.GetPC(&pc) ; PCSetType(pc,PCHYPRE) ; PCHYPRESetType(pc,"boomeramg") ; PCSetFromOptions(pc);
				CFD2CSDKsp_ksp.SetRestart(1000000) ;
				CFD2CSDKsp_ksp.SetInitialGuessNonzero() ;
				CFD2CSDKsp_ksp.SetFromOptions() ;
				
				double starttime, endtime;
				// Solve the linear system
				cout << "Force Interpolation CFD2CSD: p,r: Start Solving" << Loci::MPI_processes << delim << Loci::MPI_rank << endl ;
				starttime = MPI_Wtime();
		//		cout << "Force Interpolation CFD2CSD: p,r: Before operator" << Loci::MPI_processes << delim << Loci::MPI_rank << endl ;
				CFD2CSDKsp_ksp.SetOperators(CFD2CSD_A_matrix) ; 
		//		cout << "Force Interpolation CFD2CSD: p,r: after operator" << Loci::MPI_processes << delim << Loci::MPI_rank << endl ;
				CFD2CSDKsp_ksp.SolveVector((*L_fsiBoundaryFaceBX_nit__),(phi)) ;
				endtime = MPI_Wtime();
				cout << "Force Interpolation CFD2CSD: p,r: " << Loci::MPI_processes << delim << Loci::MPI_rank << delim << " solved in " << endtime - starttime << " sec" << endl ;
				CFD2CSDKsp_ksp.GetConvergedReason() ;
				int its = CFD2CSDKsp_ksp.GetIterationNumber() ;
				if (Loci::MPI_rank==0) cout << "Force Interpolation CFD2CSD: Petsc solver converged in " << its << " iterations" << endl ;
		//		CFD2CSDKsp_ksp.Destroy() ;
		//		CFD2CSD_A_matrix.Destroy() ;
				
				
				vector<int> recvcounts(p,0) ;
				vector<int> displs(p,0) ;
				for(int i=0;i<p;++i) {
					for(int j=0;j<i;++j) displs[i]+=(3 * ((*L_fsiNumBoundaryFaceX_nit__)[j]) ) ;
					recvcounts[i] = 3*((i==p-1)?((*L_fsiNumBoundaryFaceX_nit__)[i]+4):(*L_fsiNumBoundaryFaceX_nit__)[i]) ;
				}
				
				// Allocate fsiQ
				//(*fsiWeight).resize(3*globalNumFace) ; // number of nodes * 3 coordinates
				(*L_fsiBoundaryFaceWeightX_nit__).resize(3*globalNumFace) ; // number of nodes * 3 coordinates
		
				int counter=0;
		        PetscScalar *phixCopy,*phiyCopy,*phizCopy ;
		        (phi).GetArray(&phixCopy,&phiyCopy,&phizCopy) ;
				// Fill in fsiQ
				for(int i=localStart; i<localStart+localNumFace; ++i,++counter) {
					(*L_fsiBoundaryFaceWeightX_nit__)[3*i+0] = phixCopy[counter] ; //cout << "fsiWeight[3*" << i << "+0]=" << (*fsiWeight)[3*i+0] << ", " ;
					(*L_fsiBoundaryFaceWeightX_nit__)[3*i+1] = phiyCopy[counter] ; //cout << "fsiWeight[3*" << i << "+1]=" << (*fsiWeight)[3*i+1] << ", " ;
					(*L_fsiBoundaryFaceWeightX_nit__)[3*i+2] = phizCopy[counter] ; //cout << "fsiWeight[3*" << i << "+2]=" << (*fsiWeight)[3*i+2] << endl ;
				}
		        (phi).RestoreArray(&phixCopy,&phiyCopy,&phizCopy) ; // according to manual inexpensive
		
				// Allgatherv
				MPI_Allgatherv(&(*L_fsiBoundaryFaceWeightX_nit__)[3*localStart], 3*localNumFace, MPI_DOUBLE, &(*L_fsiBoundaryFaceWeightX_nit__)[0], &recvcounts[0], &displs[0], MPI_DOUBLE, MPI_COMM_WORLD);
//				PetscScalar *phixCopy,*phiyCopy,*phizCopy ;
//        (*$fsiBoundaryFacePhi(X){n}).GetArray(&phixCopy,&phiyCopy,&phizCopy) ;
//		// Fill in fsiQ
//			int counter=0;
//			for(int i=localStart; i<localStart+localNumFace; ++i,++counter) {
//				cout << "rank = " << rank << " phixCopy[" << i << "]=" << phixCopy[counter] << ", " ;
//				cout << "rank = " << rank << " phiyCopy[" << i << "]=" << phiyCopy[counter] << ", " ;
//				cout << "rank = " << rank << " phizCopy[" << i << "]=" << phizCopy[counter] << endl ;
//			}
//        (*$fsiBoundaryFacePhi(X){n}).RestoreArray(&phixCopy,&phiyCopy,&phizCopy) ; // according to manual inexpensive
		} // CFDIterationFinished
  }    void compute(const Loci::sequence &seq) { 
#line 807 "FSI_CFD2CSDInterpolation.loci"
      prelude(seq) ;
#line 807 "FSI_CFD2CSDInterpolation.loci"
    }
#line 807 "FSI_CFD2CSDInterpolation.loci"
} ;
#line 807 "FSI_CFD2CSDInterpolation.loci"
Loci::register_rule<file_FSI_CFD2CSDInterpolation000_1280694170m614> register_file_FSI_CFD2CSDInterpolation000_1280694170m614 ;
#line 807 "FSI_CFD2CSDInterpolation.loci"
}
#line 807 "FSI_CFD2CSDInterpolation.loci"
class FSIInternalFaceApply : public unit_rule {
    private:
	  const_param<real> r, a ;
	  const_param<int> fsiNr ;
	  const_blackbox<ublas::matrix<real,ublas::column_major> > CSDnodes_it ;
	  const_blackbox<ublas::matrix<int,ublas::column_major> > CSDConnectivity ;
	  const_blackbox<vector<real> > Qtop ;
	  const_blackbox<vector<real> > Qbottom ;
	  const_blackbox<vector<real> > fsiWeightTop ;
	  const_blackbox<vector<real> > fsiWeightBottom ;
	  const_param<vector<int> > fsiNumBoundaryFaceTop ;
	  const_param<vector<int> > fsiNumBoundaryFaceBottom ;
	  const_param<int> CSDdimension ;
	  const_param<real> CSD2dSpanCenter ;
	  const_param<bool> CFDIterationFinished ;
	  blackbox<ublas::matrix<real,ublas::column_major> > CSDForce ;
    private:
	  int globalNumFace;
    public:

      // Define input and output.
      FSIInternalFaceApply() {
	name_store("fsiBoundaryFaceQ(topCSDfaces){n,it}",Qtop) ;
	name_store("fsiBoundaryFaceQ(bottomCSDfaces){n,it}",Qbottom) ;
	name_store("fsiBoundaryFaceWeight(topCSDfaces){n,it}",fsiWeightTop) ;
	name_store("fsiBoundaryFaceWeight(bottomCSDfaces){n,it}",fsiWeightBottom) ;
	name_store("FSIRBFr",r) ;
	name_store("FSIRBFa",a) ;
	name_store("FSIRBFnr",fsiNr) ;
	name_store("fsiNumBoundaryFace(topCSDfaces)",fsiNumBoundaryFaceTop) ;
	name_store("fsiNumBoundaryFace(bottomCSDfaces)",fsiNumBoundaryFaceBottom) ;
	name_store("CSDnodes{n,it}",CSDnodes_it) ;
	name_store("CSDConnectivity",CSDConnectivity) ;
	name_store("CSDForce{n,it}",CSDForce) ;
	name_store("CSDdimension",CSDdimension) ;
	name_store("CSD2dSpanCenter",CSD2dSpanCenter) ;    
	name_store("CFDIterationFinished{n,it-1}",CFDIterationFinished) ;
	input("CFDIterationFinished{n,it-1}") ;
	input("CSDConnectivity") ;
	input("CSD2dSpanCenter") ;
	input("CSDdimension") ;
	input("FSIRBFnr, FSIRBFr, FSIRBFa");		
	input("fsiNumBoundaryFace(topCSDfaces)") ;
	input("fsiNumBoundaryFace(bottomCSDfaces)") ;
	input("fsiBoundaryFaceQ(topCSDfaces){n,it}") ;
	input("fsiBoundaryFaceQ(bottomCSDfaces){n,it}") ;
	input("fsiBoundaryFaceWeight(topCSDfaces){n,it}") ;
	input("fsiBoundaryFaceWeight(bottomCSDfaces){n,it}") ;
	input("CSDnodes{n,it}") ;
	output("CSDForce{n,it}") ;
	constraint("FSINLAMS") ;
	constraint("FSICoupling") ; //internal nodes
	constraint("topCSDfaces{n,it},bottomCSDfaces{n,it}") ; //internal nodes
  //  disable_threading() ;
      }
	  
	  // Loop over internal nodes and apply FSI interpolation
      void calculate(Entity node) {
        	

      }

      // Assemble and solve.
	  virtual void compute(const sequence &seq) {
			
			
			
        // Get the number of local and global nodes.
			const int p = Loci::MPI_processes ;
			const int rank = Loci::MPI_rank ;
	    int localNumFaceTop=(*fsiNumBoundaryFaceTop)[rank],globalNumFaceTop=0 ;
	    for(unsigned int i=0;i<(*fsiNumBoundaryFaceTop).size();++i) globalNumFaceTop+=(*fsiNumBoundaryFaceTop)[i] ;
	    int localNumFaceBottom=(*fsiNumBoundaryFaceBottom)[rank],globalNumFaceBottom=0 ;
	    for(unsigned int i=0;i<(*fsiNumBoundaryFaceBottom).size();++i) globalNumFaceBottom+=(*fsiNumBoundaryFaceBottom)[i] ;
//		cout << "p,r,localNumFace,globalNumFace: " << p << delim << rank << delim << localNumFace << delim << globalNumFace << endl ;


		

			double starttime, endtime;
			// Solve the linear system
	//		cout << "Force Interpolation CFD2CSD -> p,r: Start FSI CFD2CSD Applying" << Loci::MPI_processes << delim << Loci::MPI_rank << endl ;
			starttime = MPI_Wtime();			
			double distance = 0. ;
			vector<real> valueTop ; valueTop.resize(globalNumFaceTop); std::fill( valueTop.begin(),valueTop.end(), 0. );
			vector<real> valueBottom ; valueBottom.resize(globalNumFaceBottom); std::fill( valueBottom.begin(),valueBottom.end(), 0. );
			int CSDnn = (*CSDnodes_it).size1() ; // if (rank==0) cout << "Force Interpolation CFD2CSD -> CSDnodes.nb = " << CSDnn << endl ;
			(*CSDForce).resize(6*CSDnn,1) ; std::fill( (*CSDForce).data().begin(), (*CSDForce).data().end(), 0. );  
			
//			if (Loci::MPI_rank==0) {
//					for(int i=0; i<CSDnn; ++i) {
//						cout << "CSDForce: " << (*CSDForce)(6*i+0,0) << "\t" << (*CSDForce)(6*i+1,0) << "\t" << (*CSDForce)(6*i+2,0) << "\t" << (*CSDForce)(6*i+3,0) << "\t" << (*CSDForce)(6*i+4,0) << "\t" << (*CSDForce)(6*i+5,0) <<endl ;
//					}
//				}
			if (Loci::MPI_rank==0) cout << "FSIInternalFaceApply CFDIterationFinished" << (*CFDIterationFinished) << endl ;	
			if (*CFDIterationFinished) { 
				
		//			if (Loci::MPI_rank==0) cout << "globalNumFaceTop, globalNumFaceBottom: " << globalNumFaceTop << ", " << globalNumFaceBottom << endl ;
	//		if (Loci::MPI_rank==0) {
		//		for(int i=0; i<globalNumFaceTop; ++i) {
		//			cout << "Qtop:\t " << (*Qtop)[3*i+0] << "\t" << (*Qtop)[3*i+1] << "\t" << (*Qtop)[3*i+2] << "\tphiCFD2CSDtop:\t " << (*fsiWeightTop)[3*i+0] << "\t" << (*fsiWeightTop)[3*i+1] << "\t" << (*fsiWeightTop)[3*i+2] << endl ;
			//	}
		//	}
//			if (Loci::MPI_rank==0) {
//				for(int i=0; i<globalNumFaceBottom; ++i) {
//					cout << "Qbottom:\t " << (*Qbottom)[3*i+0] << "\t" << (*Qbottom)[3*i+1] << "\t" << (*Qbottom)[3*i+2] << "\tphiCFD2CSDbottom: " << (*fsiWeightBottom)[3*i+0] << "\t" << (*fsiWeightBottom)[3*i+1] << "\t" << (*fsiWeightBottom)[3*i+2] << endl ;
//				}
//			}

			// Compute CSD mesh element centers
			int CSDelemN = (*CSDConnectivity).size1() ;
			ublas::matrix<real,ublas::column_major> CSDelementcenters ;
			ublas::matrix<real,ublas::column_major> CSDelementnormalvectors ;
			ublas::matrix<real,ublas::column_major> CSDelementforce ;
			CSDelementcenters.resize(CSDelemN, 3) ; std::fill( CSDelementcenters.data().begin(), CSDelementcenters.data().end(), 0.0);
			CSDelementnormalvectors.resize(CSDelemN, 3) ; std::fill( CSDelementnormalvectors.data().begin(), CSDelementnormalvectors.data().end(), 0.0);
			CSDelementforce.resize(CSDelemN, 3) ; std::fill( CSDelementforce.data().begin(), CSDelementforce.data().end(), 0.0);
			//real CSDelementforce(CSDelemN) ;
			int node1, node2, node3;
			real node1v[3]; real node2v[3]; real node3v[3]; real temp1v[3]; real temp2v[3] ;
			real tempMag ;

	//		if (rank==0) cout << "[I] Force Interpolation starting, CSDelemN=" << CSDelemN << endl ;

			for(int i=0; i<CSDelemN; ++i) {
				node1=(*CSDConnectivity)(i,0)-1 ; node2=(*CSDConnectivity)(i,1)-1 ; node3=(*CSDConnectivity)(i,2)-1 ;
			//	cout << "Force Int i = " << i << endl ;cout << node1 << ", " << node2 << ", " << node3 << endl ;
				
				for(int d=0; d<3; ++d) {	
					node1v[d] = (*CSDnodes_it)( node1, d); node2v[d] = (*CSDnodes_it)( node2, d); node3v[d] = (*CSDnodes_it)( node3, d); 
					CSDelementcenters(i, d) = 1.0 / 3.0 * ( node1v[d] + node2v[d] +node3v[d] ); 
					temp1v[d] = node2v[d] - node1v[d];
					temp2v[d] = node3v[d] - node1v[d];
				}
			//	cout << "Force Int i = " << i << endl ; for (int j=0; j<3; ++j) cout << "node1v, node2v, node3v, temp1v, temp2v: " << node1v[j] << ", " << node2v[j] << ", " << node3v[j] << ", " << temp1v[j] << ", " << temp2v[j] << endl ;
				CSDelementnormalvectors(i, 0) = temp1v[1]*temp2v[2] - temp1v[2]*temp2v[1] ;
				CSDelementnormalvectors(i, 1) = temp1v[2]*temp2v[0] - temp1v[0]*temp2v[2] ;
				CSDelementnormalvectors(i, 2) = temp1v[0]*temp2v[1] - temp1v[1]*temp2v[0] ;
				tempMag = sqrt( pow(CSDelementnormalvectors(i, 0), 2.0) + pow(CSDelementnormalvectors(i, 1), 2.0) + pow(CSDelementnormalvectors(i, 2), 2.0) ) ;
				for(int d=0; d<3; ++d) CSDelementnormalvectors(i,d) /= tempMag ;
			//	cout << "Force Int i = " << i << endl ; for (int j=0; j<3; ++j) cout << "normalvector: " << CSDelementnormalvectors(i,j) << endl ;
			}
			
		//	cout << "Force Interpolation CFD2CSD -> p,r: CSD element centers and normal vectors assembled for: " << Loci::MPI_processes << delim << Loci::MPI_rank << endl ;
				
							
				for(int n=0;n<CSDelemN;++n) {
					int I,I1,I2;
					for(int i=0;i<globalNumFaceTop;++i){
						I=3*i; I1=3*i+1; I2=3*i+2;
						distance = pow((CSDelementcenters)(n,0) - (*Qtop)[I], 2.) ; // CSDnodes.y = CFD.x
						distance += pow((CSDelementcenters)(n,1) -(*Qtop)[I1], 2.) ; // no y component
						if (*CSDdimension == 3) {
							distance += pow((CSDelementcenters)(n,2) - (*Qtop)[I2], 2.) ; // CSDnodes.x = CFD.z
						} else {
							distance += pow((*CSD2dSpanCenter) - (*Qtop)[I2], 2.) ; // CSDnodes.x = CFD.z
						}
						distance = sqrt(distance); 
						valueTop[i] = radialBasisFunction(distance, *r, *a, *fsiNr) ;  
					}
					for(int i=0;i<globalNumFaceBottom;++i){
						I=3*i; I1=3*i+1; I2=3*i+2;
						distance = pow((CSDelementcenters)(n,0) - (*Qbottom)[I], 2.) ; // CSDnodes.y = CFD.x
						distance += pow((CSDelementcenters)(n,1) -(*Qbottom)[I1], 2.) ; // no y component
						if (*CSDdimension == 3) {
							distance += pow((CSDelementcenters)(n,2) - (*Qbottom)[I2], 2.) ; // CSDnodes.x = CFD.z
						} else {
							distance += pow((*CSD2dSpanCenter) - (*Qbottom)[I2], 2.) ; // CSDnodes.x = CFD.z
						}
			//			if (Loci::MPI_rank==0) cout << "Q[I2], centerspan: " << (*Qbottom)[I2] << ", " << *CSD2dSpanCenter << endl ;
						distance = sqrt(distance); 
						valueBottom[i] = radialBasisFunction(distance, *r, *a, *fsiNr) ;
					}
					for(int d=0;d<3;++d) {
						double polynomTop = (*fsiWeightTop)[3*(globalNumFaceTop+0)+d];
						polynomTop+=(*fsiWeightTop)[3*(globalNumFaceTop+1)+d]*(CSDelementcenters)(n,0) ; // CSDnodes.y = CFD.x
						polynomTop+=(*fsiWeightTop)[3*(globalNumFaceTop+2)+d]*(CSDelementcenters)(n,1) ; // CSDnodes.z = CFD.y
						if (*CSDdimension == 3) {
							polynomTop+=(*fsiWeightTop)[3*(globalNumFaceTop+3)+d]*(CSDelementcenters)(n,2) ; // CSDnodes.x = CFD.z
						} else {
							polynomTop+=(*fsiWeightTop)[3*(globalNumFaceTop+3)+d]*(*CSD2dSpanCenter) ; // CSDnodes.x = CFD.z
						}
			//			cout << "polytop = " << polynomTop << endl ;
						double polynomBottom = (*fsiWeightBottom)[3*(globalNumFaceBottom+0)+d]; //cout << "polynomBottom0 = " << polynomBottom << "size, globalNumFacebottom, d" << (*fsiWeightBottom).size() << ", " << globalNumFaceBottom << ", " << d << endl ;
						polynomBottom+=(*fsiWeightBottom)[3*(globalNumFaceBottom+1)+d]*(CSDelementcenters)(n,0) ; //cout << "polynomBottom1 = " << polynomBottom << endl ; // CSDnodes.y = CFD.x
						polynomBottom+=(*fsiWeightBottom)[3*(globalNumFaceBottom+2)+d]*(CSDelementcenters)(n,1) ; //cout << "polynomBottom2 = " << polynomBottom << endl ;// CSDnodes.z = CFD.y
						if (*CSDdimension == 3) {
							polynomBottom+=(*fsiWeightBottom)[3*(globalNumFaceBottom+3)+d]*(CSDelementcenters)(n,2) ;// cout << "polynomBottom3a = " << polynomBottom << endl ;// CSDnodes.x = CFD.z
						} else {
							polynomBottom+=(*fsiWeightBottom)[3*(globalNumFaceBottom+3)+d]*(*CSD2dSpanCenter) ; //cout << "polynomBottom3b = " << polynomBottom << endl ;// CSDnodes.x = CFD.z
						}
				//		cout << "polynomBottom = " << polynomBottom << endl ;
						(CSDelementforce)(n,d) = - (polynomTop - polynomBottom) ;   // i = 6 * n + d
					}
					
				// FSI contribution
				real temp = 0.0 ;
					for(int i=0;i<globalNumFaceTop;++i){		
						for(int d=0;d<3;++d) {
							//temp = (*fsiWeightTop)[3*i+d] * valueTop[i] ; //cout << "temp = " << temp << "CSDForce = " << (*CSDForce)(6*n+d,0) ;; 							
							(CSDelementforce)(n,d) -= (*fsiWeightTop)[3*i+d] * valueTop[i] ;
							//cout << "CSDForce(6*" << n << "+" << d << ") = " << (*CSDForce)(6*n+d,0) << " i: " << i << " fsiWeigtTop = " << (*fsiWeightTop)[3*i+d] << " valueTop = " << valueTop[i] <<  endl;
						}
					}
					for(int i=0;i<globalNumFaceBottom;++i){		
						for(int d=0;d<3;++d) {
							(CSDelementforce)(n,d) += (*fsiWeightBottom)[3*i+d] * valueBottom[i] ;
							//cout << "CSDForce(6*" << n << "+" << d << ") = " << (*CSDForce)(6*n+d,0) << " i: " << i << "fsiWeightBottom = " << (*fsiWeightBottom)[3*i+d] << " valueBottom = " << valueBottom[i] <<  endl;
						}
					}
				}			
				endtime = MPI_Wtime();

				if (rank==0) cout << "[I] Force Interpolated onto the CSD element centers" << endl ;

				int nodeNR;

				for(int e=0; e<CSDelemN; ++e) {
					for(int n=0; n<3; ++n) {
						nodeNR = (*CSDConnectivity)(e,n)-1 ;
						//if (rank==0) cout << "nodeNR, CSDforce_index, CSDnn: " << nodeNR << ", " << 6*nodeNR << ", " << 6*nodeNR+1 << ", " << 6*nodeNR+2 << ", " << CSDnn << ", " << endl ;
						//if (rank==0) cout  << "csdelementforce: size, value" << (CSDelementforce).size1() << ", " << (CSDelementforce)(nodeNR,1) << endl ;
						//if (rank==0) cout  << "CSDelementnormalvectors: " << CSDelementnormalvectors(nodeNR,0) << ", " <<  CSDelementnormalvectors(nodeNR,1) << ", "  << CSDelementnormalvectors(nodeNR,2) << endl ;
						for(int d=0; d<3; ++d) {
							(*CSDForce)(6*nodeNR+d,0) += 1.0 / 3.0 * (CSDelementforce)(e,1) * CSDelementnormalvectors(e,d) ; // CSDForce already initialized, only y-component important
						}
						//if (rank==0) cout << "nodeNR, CSDforce: " << nodeNR << ", " << (*CSDForce)(6*nodeNR+0,0) << ", " << (*CSDForce)(6*nodeNR+1,0) << ", " << (*CSDForce)(6*nodeNR+2,0) << endl ;
					}
				}			
				
				if (rank==0) cout << "[I] Force Interpolated onto the CSD nodes" << endl ;
				
				if (*CSDdimension == 2) { // make explicitly zero
					for(int i=0; i<CSDnn; ++i) {
						(*CSDForce)(6*i+2,0) = 0.0 ;
					}
				}
//				if (Loci::MPI_rank==0) {
//					for(int i=0; i<CSDnn; ++i) {
//						cout << "valuetop: " << valueTop[i] << "\t" << "valuebottom: " << valueBottom[i] << "CSDnodes_it" << (*CSDnodes_it)(i,0) << "CSDnodes_it" << (*CSDnodes_it)(i,1) << "CSDnodes_it" << (*CSDnodes_it)(i,2) << endl ;
//					}
//				}
				if (Loci::MPI_rank==0) {
					for(int i=0; i<CSDnn; ++i) {
		//				cout << "CSDnodes: " << (*CSDnodes_it)(i,0) << "\t" << (*CSDnodes_it)(i,1) << "\t" << (*CSDnodes_it)(i,2) << "\t" << "CSDForce: " << (*CSDForce)(6*i+0,0) << "\t" << (*CSDForce)(6*i+1,0) << "\t" << (*CSDForce)(6*i+2,0) << "\t" << (*CSDForce)(6*i+3,0) << "\t" << (*CSDForce)(6*i+4,0) << "\t" << (*CSDForce)(6*i+5,0) <<endl ;
					}
				}
				//cout << "Force Interpolation CFD2CSD -> p,r: End FSI CFD2CSD" << p << delim << rank << delim << "Applied in " << endtime - starttime << " sec" << endl ;
				if (rank==0) cout << "[I] Force Interpolation completed" << endl ;
				//if (rank==0) cout << "Force Interpolation, CSDnodes: " << (*CSDnodes_it) << endl; 
				//if (rank==0) cout << "Force Interpolation, CSDForce: " << (*CSDForce) << endl; 
		} // CFDIterationFinished
			
		}	
  } ;

  register_rule<FSIInternalFaceApply> registerFSIInternalFaceApply ;

 // Rule to apply interpolation to the internal nodes
  class FSIInternalFaceApplyDummy : public apply_rule<blackbox<ublas::matrix<real,ublas::column_major> >,Loci::NullOp<ublas::matrix<real,ublas::column_major> > >  {
    private:
	  const_param<real> r, a ;
	  const_param<int> fsiNr ;
	  const_blackbox<ublas::matrix<real,ublas::column_major> > CSDnodes_it ;
	  const_blackbox<vector<real> > Qtop ;
	  const_blackbox<vector<real> > Qbottom ;
	  const_blackbox<vector<real> > fsiWeightTop ;
	  const_blackbox<vector<real> > fsiWeightBottom ;
	  const_param<vector<int> > fsiNumBoundaryFaceTop ;
	  const_param<vector<int> > fsiNumBoundaryFaceBottom ;
	  const_param<int> CSDdimension ;
	  const_param<real> CSD2dSpanCenter ;
//	  const_param<bool> CFDIterationFinished ;
	  blackbox<ublas::matrix<real,ublas::column_major> > CSDForce ;
    private:
	  	int globalNumFace;
    public:

      // Define input and output.
      FSIInternalFaceApplyDummy() {
		name_store("fsiBoundaryFaceQ(topCSDfaces){n,it}",Qtop) ;
		name_store("fsiBoundaryFaceQ(bottomCSDfaces){n,it}",Qbottom) ;
		name_store("fsiBoundaryFaceWeight(topCSDfaces){n,it}",fsiWeightTop) ;
		name_store("fsiBoundaryFaceWeight(bottomCSDfaces){n,it}",fsiWeightBottom) ;
		name_store("FSIRBFr",r) ;
		name_store("FSIRBFa",a) ;
		name_store("FSIRBFnr",fsiNr) ;
    name_store("fsiNumBoundaryFace(topCSDfaces)",fsiNumBoundaryFaceTop) ;
    name_store("fsiNumBoundaryFace(bottomCSDfaces)",fsiNumBoundaryFaceBottom) ;
    name_store("CSDnodes{n,it}",CSDnodes_it) ;
    name_store("CSDForce{n,it}",CSDForce) ;
    name_store("CSDdimension",CSDdimension) ;
    name_store("CSD2dSpanCenter",CSD2dSpanCenter) ;    
//		name_store("CFDIterationFinished{n,it-1}",CFDIterationFinished) ;
 //   input("CFDIterationFinished{n,it-1}") ;
    input("CSD2dSpanCenter") ;
    input("CSDdimension") ;
		input("FSIRBFnr, FSIRBFr, FSIRBFa");		
		input("fsiNumBoundaryFace(topCSDfaces)") ;
    input("fsiNumBoundaryFace(bottomCSDfaces)") ;
		input("fsiBoundaryFaceQ(topCSDfaces){n,it}") ;
		input("fsiBoundaryFaceQ(bottomCSDfaces){n,it}") ;
		input("fsiBoundaryFaceWeight(topCSDfaces){n,it}") ;
		input("fsiBoundaryFaceWeight(bottomCSDfaces){n,it}") ;
		input("CSDnodes{n,it}") ;
    output("CSDForce{n,it}") ;
    constraint("FSICoupling") ; //internal nodes
	constraint("FSINLAMS") ;
    constraint("topCSDfaces{n,it},bottomCSDfaces{n,it}") ; //internal nodes
	//	    disable_threading() ;
      }
	  
	  // Loop over internal nodes and apply FSI interpolation
      void calculate(Entity node) {
        	

      }

      // Assemble and solve.
	  virtual void compute(const sequence &seq) {}	
  } ;

  register_rule<FSIInternalFaceApplyDummy> registerFSIInternalFaceApplyDummy ;

class FSIInternalFaceApplyFSIEULERBEAM : public unit_rule {
    private:
	  const_param<real> r, a ;
	  const_param<int> fsiNr ;
	  const_blackbox<ublas::matrix<real,ublas::column_major> > CSDnodes_it ;
	  const_blackbox<ublas::matrix<int,ublas::column_major> > CSDConnectivity ;
	  const_blackbox<vector<real> > Qtop ;
	  const_blackbox<vector<real> > Qbottom ;
	  const_blackbox<vector<real> > fsiWeightTop ;
	  const_blackbox<vector<real> > fsiWeightBottom ;
	  const_param<vector<int> > fsiNumBoundaryFaceTop ;
	  const_param<vector<int> > fsiNumBoundaryFaceBottom ;
	  const_param<int> CSDdimension ;
	  const_param<real> CSD2dSpanCenter ;
	  const_param<bool> CFDIterationFinished ;
	  blackbox<ublas::matrix<real,ublas::column_major> > CSDForce ;
    private:
	  int globalNumFace;
    public:

      // Define input and output.
      FSIInternalFaceApplyFSIEULERBEAM() {
		name_store("fsiBoundaryFaceQ(topCSDfaces){n,it}",Qtop) ;
		name_store("fsiBoundaryFaceQ(bottomCSDfaces){n,it}",Qbottom) ;
		name_store("fsiBoundaryFaceWeight(topCSDfaces){n,it}",fsiWeightTop) ;
		name_store("fsiBoundaryFaceWeight(bottomCSDfaces){n,it}",fsiWeightBottom) ;
		name_store("FSIRBFr",r) ;
		name_store("FSIRBFa",a) ;
		name_store("FSIRBFnr",fsiNr) ;
    name_store("fsiNumBoundaryFace(topCSDfaces)",fsiNumBoundaryFaceTop) ;
    name_store("fsiNumBoundaryFace(bottomCSDfaces)",fsiNumBoundaryFaceBottom) ;
    name_store("CSDnodes{n,it}",CSDnodes_it) ;
//	name_store("CSDConnectivity{n}",CSDConnectivity) ;
    name_store("CSDForce{n,it}",CSDForce) ;
    name_store("CSDdimension",CSDdimension) ;
    name_store("CSD2dSpanCenter",CSD2dSpanCenter) ;    
    name_store("CFDIterationFinished{n,it-1}",CFDIterationFinished) ;
    input("CFDIterationFinished{n,it-1}") ;
//	input("CSDConnectivity{n}") ;
    input("CSD2dSpanCenter") ;
    input("CSDdimension") ;
    input("FSIRBFnr, FSIRBFr, FSIRBFa");		
    input("fsiNumBoundaryFace(topCSDfaces)") ;
    input("fsiNumBoundaryFace(bottomCSDfaces)") ;
    input("fsiBoundaryFaceQ(topCSDfaces){n,it}") ;
    input("fsiBoundaryFaceQ(bottomCSDfaces){n,it}") ;
    input("fsiBoundaryFaceWeight(topCSDfaces){n,it}") ;
    input("fsiBoundaryFaceWeight(bottomCSDfaces){n,it}") ;
    input("CSDnodes{n,it}") ;
    output("CSDForce{n,it}") ;
    constraint("FSICoupling") ; //internal nodes	
    constraint("FSIEULERBEAM") ;
    constraint("topCSDfaces{n,it},bottomCSDfaces{n,it}") ; //internal nodes
  //  disable_threading() ;
      }
	  
	  // Loop over internal nodes and apply FSI interpolation
      void calculate(Entity node) {
        	

      }

      // Assemble and solve.
	  virtual void compute(const sequence &seq) {
			
			
			
        // Get the number of local and global nodes.
			const int p = Loci::MPI_processes ;
			const int rank = Loci::MPI_rank ;
	    int localNumFaceTop=(*fsiNumBoundaryFaceTop)[rank],globalNumFaceTop=0 ;
	    for(unsigned int i=0;i<(*fsiNumBoundaryFaceTop).size();++i) globalNumFaceTop+=(*fsiNumBoundaryFaceTop)[i] ;
	    int localNumFaceBottom=(*fsiNumBoundaryFaceBottom)[rank],globalNumFaceBottom=0 ;
	    for(unsigned int i=0;i<(*fsiNumBoundaryFaceBottom).size();++i) globalNumFaceBottom+=(*fsiNumBoundaryFaceBottom)[i] ;
//		cout << "p,r,localNumFace,globalNumFace: " << p << delim << rank << delim << localNumFace << delim << globalNumFace << endl ;


		

			double starttime, endtime;
			// Solve the linear system
	//		cout << "Force Interpolation CFD2CSD -> p,r: Start FSI CFD2CSD Applying" << Loci::MPI_processes << delim << Loci::MPI_rank << endl ;
			starttime = MPI_Wtime();			
			double distance = 0. ;
			vector<real> valueTop ; valueTop.resize(globalNumFaceTop); std::fill( valueTop.begin(),valueTop.end(), 0. );
			vector<real> valueBottom ; valueBottom.resize(globalNumFaceBottom); std::fill( valueBottom.begin(),valueBottom.end(), 0. );
			int CSDnn = (*CSDnodes_it).size1() ; // if (rank==0) cout << "Force Interpolation CFD2CSD -> CSDnodes.nb = " << CSDnn << endl ;
			(*CSDForce).resize(6*CSDnn,1) ; std::fill( (*CSDForce).data().begin(), (*CSDForce).data().end(), 0. );  
			
//			if (Loci::MPI_rank==0) {
//					for(int i=0; i<CSDnn; ++i) {
//						cout << "CSDForce: " << (*CSDForce)(6*i+0,0) << "\t" << (*CSDForce)(6*i+1,0) << "\t" << (*CSDForce)(6*i+2,0) << "\t" << (*CSDForce)(6*i+3,0) << "\t" << (*CSDForce)(6*i+4,0) << "\t" << (*CSDForce)(6*i+5,0) <<endl ;
//					}
//				}
			if (Loci::MPI_rank==0) cout << "FSIInternalFaceApplyFSIEULERBEAM CFDIterationFinished" << (*CFDIterationFinished) << endl ;		
			if (*CFDIterationFinished) { 
				
		//			if (Loci::MPI_rank==0) cout << "globalNumFaceTop, globalNumFaceBottom: " << globalNumFaceTop << ", " << globalNumFaceBottom << endl ;
	//		if (Loci::MPI_rank==0) {
		//		for(int i=0; i<globalNumFaceTop; ++i) {
		//			cout << "Qtop:\t " << (*Qtop)[3*i+0] << "\t" << (*Qtop)[3*i+1] << "\t" << (*Qtop)[3*i+2] << "\tphiCFD2CSDtop:\t " << (*fsiWeightTop)[3*i+0] << "\t" << (*fsiWeightTop)[3*i+1] << "\t" << (*fsiWeightTop)[3*i+2] << endl ;
			//	}
		//	}
//			if (Loci::MPI_rank==0) {
//				for(int i=0; i<globalNumFaceBottom; ++i) {
//					cout << "Qbottom:\t " << (*Qbottom)[3*i+0] << "\t" << (*Qbottom)[3*i+1] << "\t" << (*Qbottom)[3*i+2] << "\tphiCFD2CSDbottom: " << (*fsiWeightBottom)[3*i+0] << "\t" << (*fsiWeightBottom)[3*i+1] << "\t" << (*fsiWeightBottom)[3*i+2] << endl ;
//				}
//			}

			// Compute CSD mesh element centers
//			int CSDelemN = (*CSDConnectivity).size1() ;
//			ublas::matrix<real,ublas::column_major> CSDelementcenters ;
//			ublas::matrix<real,ublas::column_major> CSDelementnormalvectors ;
//			ublas::matrix<real,ublas::column_major> CSDelementforce ;
//			CSDelementcenters.resize(CSDelemN, 3) ; std::fill( CSDelementcenters.data().begin(), CSDelementcenters.data().end(), 0.0);
//			CSDelementnormalvectors.resize(CSDelemN, 3) ; std::fill( CSDelementnormalvectors.data().begin(), CSDelementnormalvectors.data().end(), 0.0);
//			CSDelementforce.resize(CSDelemN, 3) ; std::fill( CSDelementforce.data().begin(), CSDelementforce.data().end(), 0.0);
//			//real CSDelementforce(CSDelemN) ;
//			int node1, node2, node3;
//			real node1v[3]; real node2v[3]; real node3v[3]; real temp1v[3]; real temp2v[3] ;
//			real tempMag ;
//
//	//		if (rank==0) cout << "[I] Force Interpolation starting, CSDelemN=" << CSDelemN << endl ;
//
//			for(int i=0; i<CSDelemN; ++i) {
//				node1=(*CSDConnectivity)(i,0)-1 ; node2=(*CSDConnectivity)(i,1)-1 ; node3=(*CSDConnectivity)(i,2)-1 ;
//			//	cout << "Force Int i = " << i << endl ;cout << node1 << ", " << node2 << ", " << node3 << endl ;
//				
//				for(int d=0; d<3; ++d) {	
//					node1v[d] = (*CSDnodes_it)( node1, d); node2v[d] = (*CSDnodes_it)( node2, d); node3v[d] = (*CSDnodes_it)( node3, d); 
//					CSDelementcenters(i, d) = 1.0 / 3.0 * ( node1v[d] + node2v[d] +node3v[d] ); 
//					temp1v[d] = node2v[d] - node1v[d];
//					temp2v[d] = node3v[d] - node1v[d];
//				}
//			//	cout << "Force Int i = " << i << endl ; for (int j=0; j<3; ++j) cout << "node1v, node2v, node3v, temp1v, temp2v: " << node1v[j] << ", " << node2v[j] << ", " << node3v[j] << ", " << temp1v[j] << ", " << temp2v[j] << endl ;
//				CSDelementnormalvectors(i, 0) = temp1v[1]*temp2v[2] - temp1v[2]*temp2v[1] ;
//				CSDelementnormalvectors(i, 1) = temp1v[2]*temp2v[0] - temp1v[0]*temp2v[2] ;
//				CSDelementnormalvectors(i, 2) = temp1v[0]*temp2v[1] - temp1v[1]*temp2v[0] ;
//				tempMag = sqrt( pow(CSDelementnormalvectors(i, 0), 2.0) + pow(CSDelementnormalvectors(i, 1), 2.0) + pow(CSDelementnormalvectors(i, 2), 2.0) ) ;
//				for(int d=0; d<3; ++d) CSDelementnormalvectors(i,d) /= tempMag ;
//			//	cout << "Force Int i = " << i << endl ; for (int j=0; j<3; ++j) cout << "normalvector: " << CSDelementnormalvectors(i,j) << endl ;
//			}
//			
		//	cout << "Force Interpolation CFD2CSD -> p,r: CSD element centers and normal vectors assembled for: " << Loci::MPI_processes << delim << Loci::MPI_rank << endl ;
				
							
				for(int n=0;n<CSDnn;++n) {
					int I,I1,I2;
					for(int i=0;i<globalNumFaceTop;++i){
						I=3*i; I1=3*i+1; I2=3*i+2;
						distance = pow((*CSDnodes_it)(n,0) - (*Qtop)[I], 2.) ; // CSDnodes.y = CFD.x
						distance += pow((*CSDnodes_it)(n,1) -(*Qtop)[I1], 2.) ; // no y component
						if (*CSDdimension == 3) {
							distance += pow((*CSDnodes_it)(n,2) - (*Qtop)[I2], 2.) ; // CSDnodes.x = CFD.z
						} else {
							distance += pow((*CSD2dSpanCenter) - (*Qtop)[I2], 2.) ; // CSDnodes.x = CFD.z
						}
						distance = sqrt(distance); 
						valueTop[i] = radialBasisFunction(distance, *r, *a, *fsiNr) ;  
					}
					for(int i=0;i<globalNumFaceBottom;++i){
						I=3*i; I1=3*i+1; I2=3*i+2;
						distance = pow((*CSDnodes_it)(n,0) - (*Qbottom)[I], 2.) ; // CSDnodes.y = CFD.x
						distance += pow((*CSDnodes_it)(n,1) -(*Qbottom)[I1], 2.) ; // no y component
						if (*CSDdimension == 3) {
							distance += pow((*CSDnodes_it)(n,2) - (*Qbottom)[I2], 2.) ; // CSDnodes.x = CFD.z
						} else {
							distance += pow((*CSD2dSpanCenter) - (*Qbottom)[I2], 2.) ; // CSDnodes.x = CFD.z
						}
			//			if (Loci::MPI_rank==0) cout << "Q[I2], centerspan: " << (*Qbottom)[I2] << ", " << *CSD2dSpanCenter << endl ;
						distance = sqrt(distance); 
						valueBottom[i] = radialBasisFunction(distance, *r, *a, *fsiNr) ;
					}
					for(int d=0;d<3;++d) {
						double polynomTop = (*fsiWeightTop)[3*(globalNumFaceTop+0)+d];
						polynomTop+=(*fsiWeightTop)[3*(globalNumFaceTop+1)+d]*(*CSDnodes_it)(n,0) ; // CSDnodes.y = CFD.x
						polynomTop+=(*fsiWeightTop)[3*(globalNumFaceTop+2)+d]*(*CSDnodes_it)(n,1) ; // CSDnodes.z = CFD.y
						if (*CSDdimension == 3) {
							polynomTop+=(*fsiWeightTop)[3*(globalNumFaceTop+3)+d]*(*CSDnodes_it)(n,2) ; // CSDnodes.x = CFD.z
						} else {
							polynomTop+=(*fsiWeightTop)[3*(globalNumFaceTop+3)+d]*(*CSD2dSpanCenter) ; // CSDnodes.x = CFD.z
						}
			//			cout << "polytop = " << polynomTop << endl ;
						double polynomBottom = (*fsiWeightBottom)[3*(globalNumFaceBottom+0)+d]; //cout << "polynomBottom0 = " << polynomBottom << "size, globalNumFacebottom, d" << (*fsiWeightBottom).size() << ", " << globalNumFaceBottom << ", " << d << endl ;
						polynomBottom+=(*fsiWeightBottom)[3*(globalNumFaceBottom+1)+d]*(*CSDnodes_it)(n,0) ; //cout << "polynomBottom1 = " << polynomBottom << endl ; // CSDnodes.y = CFD.x
						polynomBottom+=(*fsiWeightBottom)[3*(globalNumFaceBottom+2)+d]*(*CSDnodes_it)(n,1) ; //cout << "polynomBottom2 = " << polynomBottom << endl ;// CSDnodes.z = CFD.y
						if (*CSDdimension == 3) {
							polynomBottom+=(*fsiWeightBottom)[3*(globalNumFaceBottom+3)+d]*(*CSDnodes_it)(n,2) ;// cout << "polynomBottom3a = " << polynomBottom << endl ;// CSDnodes.x = CFD.z
						} else {
							polynomBottom+=(*fsiWeightBottom)[3*(globalNumFaceBottom+3)+d]*(*CSD2dSpanCenter) ; //cout << "polynomBottom3b = " << polynomBottom << endl ;// CSDnodes.x = CFD.z
						}
				//		cout << "polynomBottom = " << polynomBottom << endl ;
						(*CSDForce)(6*n+d,0) = - (polynomTop - polynomBottom) ;   // i = 6 * n + d
					}
					
				// FSI contribution
				real temp = 0.0 ;
					for(int i=0;i<globalNumFaceTop;++i){		
						for(int d=0;d<3;++d) {
							//temp = (*fsiWeightTop)[3*i+d] * valueTop[i] ; //cout << "temp = " << temp << "CSDForce = " << (*CSDForce)(6*n+d,0) ;; 							
							(*CSDForce)(6*n+d,0) -= (*fsiWeightTop)[3*i+d] * valueTop[i] ;
							//cout << "CSDForce(6*" << n << "+" << d << ") = " << (*CSDForce)(6*n+d,0) << " i: " << i << " fsiWeigtTop = " << (*fsiWeightTop)[3*i+d] << " valueTop = " << valueTop[i] <<  endl;
						}
					}
					for(int i=0;i<globalNumFaceBottom;++i){		
						for(int d=0;d<3;++d) {
							(*CSDForce)(6*n+d,0) += (*fsiWeightBottom)[3*i+d] * valueBottom[i] ;
							//cout << "CSDForce(6*" << n << "+" << d << ") = " << (*CSDForce)(6*n+d,0) << " i: " << i << "fsiWeightBottom = " << (*fsiWeightBottom)[3*i+d] << " valueBottom = " << valueBottom[i] <<  endl;
						}
					}
				}			
				endtime = MPI_Wtime();

				if (rank==0) cout << "[I] Force Interpolated onto the CSD element centers" << endl ;

//				int nodeNR;
//
//				for(int e=0; e<CSDelemN; ++e) {
//					for(int n=0; n<3; ++n) {
//						nodeNR = (*CSDConnectivity)(e,n)-1 ;
//						//if (rank==0) cout << "nodeNR, CSDforce_index, CSDnn: " << nodeNR << ", " << 6*nodeNR << ", " << 6*nodeNR+1 << ", " << 6*nodeNR+2 << ", " << CSDnn << ", " << endl ;
//						//if (rank==0) cout  << "csdelementforce: size, value" << (CSDelementforce).size1() << ", " << (CSDelementforce)(nodeNR,1) << endl ;
//						//if (rank==0) cout  << "CSDelementnormalvectors: " << CSDelementnormalvectors(nodeNR,0) << ", " <<  CSDelementnormalvectors(nodeNR,1) << ", "  << CSDelementnormalvectors(nodeNR,2) << endl ;
//						for(int d=0; d<3; ++d) {
//							(*CSDForce)(6*nodeNR+d,0) += 1.0 / 3.0 * (CSDelementforce)(e,1) * CSDelementnormalvectors(e,d) ; // CSDForce already initialized, only y-component important
//						}
//						//if (rank==0) cout << "nodeNR, CSDforce: " << nodeNR << ", " << (*CSDForce)(6*nodeNR+0,0) << ", " << (*CSDForce)(6*nodeNR+1,0) << ", " << (*CSDForce)(6*nodeNR+2,0) << endl ;
//					}
//				}			
//				
//				if (rank==0) cout << "[I] Force Interpolated onto the CSD nodes" << endl ;
				
				if (*CSDdimension == 2) { // make explicitly zero
					for(int i=0; i<CSDnn; ++i) {
						(*CSDForce)(6*i+2,0) = 0.0 ;
					}
				}
//				if (Loci::MPI_rank==0) {
//					for(int i=0; i<CSDnn; ++i) {
//						cout << "valuetop: " << valueTop[i] << "\t" << "valuebottom: " << valueBottom[i] << "CSDnodes_it" << (*CSDnodes_it)(i,0) << "CSDnodes_it" << (*CSDnodes_it)(i,1) << "CSDnodes_it" << (*CSDnodes_it)(i,2) << endl ;
//					}
//				}
				if (Loci::MPI_rank==0) {
					for(int i=0; i<CSDnn; ++i) {
		//				cout << "CSDnodes: " << (*CSDnodes_it)(i,0) << "\t" << (*CSDnodes_it)(i,1) << "\t" << (*CSDnodes_it)(i,2) << "\t" << "CSDForce: " << (*CSDForce)(6*i+0,0) << "\t" << (*CSDForce)(6*i+1,0) << "\t" << (*CSDForce)(6*i+2,0) << "\t" << (*CSDForce)(6*i+3,0) << "\t" << (*CSDForce)(6*i+4,0) << "\t" << (*CSDForce)(6*i+5,0) <<endl ;
					}
				}
				//cout << "Force Interpolation CFD2CSD -> p,r: End FSI CFD2CSD" << p << delim << rank << delim << "Applied in " << endtime - starttime << " sec" << endl ;
				if (rank==0) cout << "[I] Force Interpolation completed" << endl ;
				//if (rank==0) cout << "Force Interpolation, CSDnodes: " << (*CSDnodes_it) << endl; 
				//if (rank==0) cout << "Force Interpolation, CSDForce: " << (*CSDForce) << endl; 
		} // CFDIterationFinished 
			
		}	
  } ;

  register_rule<FSIInternalFaceApplyFSIEULERBEAM> registerFSIInternalFaceApplyFSIEULERBEAM ;

 // Rule to apply interpolation to the internal nodes
  class FSIInternalFaceApplyDummyFSIEULERBEAM : public apply_rule<blackbox<ublas::matrix<real,ublas::column_major> >,Loci::NullOp<ublas::matrix<real,ublas::column_major> > >  {
    private:
	  const_param<real> r, a ;
	  const_param<int> fsiNr ;
	  const_blackbox<ublas::matrix<real,ublas::column_major> > CSDnodes_it ;
	  const_blackbox<vector<real> > Qtop ;
	  const_blackbox<vector<real> > Qbottom ;
	  const_blackbox<vector<real> > fsiWeightTop ;
	  const_blackbox<vector<real> > fsiWeightBottom ;
	  const_param<vector<int> > fsiNumBoundaryFaceTop ;
	  const_param<vector<int> > fsiNumBoundaryFaceBottom ;
	  const_param<int> CSDdimension ;
	  const_param<real> CSD2dSpanCenter ;
//	  const_param<bool> CFDIterationFinished ;
	  blackbox<ublas::matrix<real,ublas::column_major> > CSDForce ;
    private:
	  	int globalNumFace;
    public:

      // Define input and output.
      FSIInternalFaceApplyDummyFSIEULERBEAM() {
		name_store("fsiBoundaryFaceQ(topCSDfaces){n,it}",Qtop) ;
		name_store("fsiBoundaryFaceQ(bottomCSDfaces){n,it}",Qbottom) ;
		name_store("fsiBoundaryFaceWeight(topCSDfaces){n,it}",fsiWeightTop) ;
		name_store("fsiBoundaryFaceWeight(bottomCSDfaces){n,it}",fsiWeightBottom) ;
		name_store("FSIRBFr",r) ;
		name_store("FSIRBFa",a) ;
		name_store("FSIRBFnr",fsiNr) ;
    name_store("fsiNumBoundaryFace(topCSDfaces)",fsiNumBoundaryFaceTop) ;
    name_store("fsiNumBoundaryFace(bottomCSDfaces)",fsiNumBoundaryFaceBottom) ;
    name_store("CSDnodes{n,it}",CSDnodes_it) ;
    name_store("CSDForce{n,it}",CSDForce) ;
    name_store("CSDdimension",CSDdimension) ;
    name_store("CSD2dSpanCenter",CSD2dSpanCenter) ;    
//		name_store("CFDIterationFinished{n,it-1}",CFDIterationFinished) ;
 //   input("CFDIterationFinished{n,it-1}") ;
    input("CSD2dSpanCenter") ;
    input("CSDdimension") ;
		input("FSIRBFnr, FSIRBFr, FSIRBFa");		
		input("fsiNumBoundaryFace(topCSDfaces)") ;
    input("fsiNumBoundaryFace(bottomCSDfaces)") ;
		input("fsiBoundaryFaceQ(topCSDfaces){n,it}") ;
		input("fsiBoundaryFaceQ(bottomCSDfaces){n,it}") ;
		input("fsiBoundaryFaceWeight(topCSDfaces){n,it}") ;
		input("fsiBoundaryFaceWeight(bottomCSDfaces){n,it}") ;
		input("CSDnodes{n,it}") ;
    output("CSDForce{n,it}") ;
    constraint("FSICoupling") ; //internal nodes
	 constraint("FSIEULERBEAM") ;
    constraint("topCSDfaces{n,it},bottomCSDfaces{n,it}") ; //internal nodes
	//	    disable_threading() ;
      }
	  
	  // Loop over internal nodes and apply FSI interpolation
      void calculate(Entity node) {
        	

      }

      // Assemble and solve.
	  virtual void compute(const sequence &seq) {}	
  } ;

  register_rule<FSIInternalFaceApplyDummyFSIEULERBEAM> registerFSIInternalFaceApplyDummyFSIEULERBEAM ;

}


