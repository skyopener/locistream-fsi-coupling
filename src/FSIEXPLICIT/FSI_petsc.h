//-----------------------------------------------------------------------------
// Description: This file contains rules for implementing a RBF linear solver
// using PETSC.
//-----------------------------------------------------------------------------

#ifndef FSI_PETSC_H
#define FSI_PETSC_H

// mpi
#include<mpi.h>

// Standard library includes.
#include<vector>
#include<string>
using std::vector ;

// PETSC includes.
#include <petsc.h>
#include <petscerror.h>
#include <petscksp.h>

// StreamUns includes.
#include "sciTypes.h"

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
      mutable Petsc::Mat m ;
    public:
      PETSCMatrix() : created(false) {}
      ~PETSCMatrix() { if(created) MatDestroy(m) ; }
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
	  int GetInfo() {
		  MatInfo info;
		  double mal, nz_alloc, nz_used, nz_un ;
		  MatGetInfo(m,MAT_LOCAL,&info);
		  mal = info.mallocs; nz_alloc=info.nz_allocated; nz_used=info.nz_used; nz_un=info.nz_unneeded ;
		  cout << "PETSC info: Number of mallocs during MatSetValues(): " << mal << delim << "Number of nonzeros (allocated, used, unneeded): " << nz_alloc << delim << nz_used << delim << nz_un << endl ;
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
      int SetOperators(const PETSCMatrix &m) const {
        int ierr=KSPSetOperators(ksp,m.Data(),m.Data(),DIFFERENT_NONZERO_PATTERN) ; // SAME_NONZERO_PATTERN
        CHKERRQ(ierr) ; return 0 ;
      }
      void SetTolerances(real relativeTolerance,real absoluteTolerance,int maxIterations) const {
        KSPSetTolerances(ksp,relativeTolerance,absoluteTolerance,PETSC_DEFAULT,maxIterations) ;
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
}


#endif
