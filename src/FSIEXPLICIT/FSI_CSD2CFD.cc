#line 1 "FSI_CSD2CFD.loci"
//-----------------------------------------------------------------------------
// Description: This file contains rules for assembling and solving the linear
//   elasticity equations for the movement of interior nodes given the
//   boundary node displacements.
//-----------------------------------------------------------------------------
                                                                                
// Standard library includes.
#include <iostream>
#include <iomanip>
#include <fstream>
//#include <cstdlib> // for exit function
#include <string>
#include <sstream>
//#include <cmath>
#include <vector>
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

//#include "hypre.h"

// Create a namespace for the PETSC objects so there is no confusion. Note that
// there is already a type "Mat" defined by Loci.
namespace Petsc {
  typedef struct _p_Vec* Vec ;
  typedef struct _p_Mat* Mat ;
  typedef struct _p_KSP* KSP ;
} ;

                                      
// boost::ublas includes
#include <boost/numeric/ublas/matrix.hpp>

namespace ublas = boost::numeric::ublas ;                        
                        
#include <mpi.h>
                                                                                
// StreamUns includes.
#include "sciTypes.h"                        
                        
// CSDvariables
#line 1 "FSI_CSDvariables.lh"
// $type CSDE1 param<real>  // Young's modulus 1
// $type CSDE2 param<real>  // Young's modulus 2
// $type CSDnu12 param<real>  // Poisson ratio 12
// $type CSDnu21 param<real>  // Poisson ratio 21
// $type CSDG12 param<real>  // Shear modulus 12
// $type CSDrhoStructure param<real>  // structure density
// $type CSDthicknessStructure param<real>  // structure thickness
// $type CSDintegrationScheme param<int>  // integrationScheme for the CSD-> 1:newmark, 2:generalized alpha
// $type CSDdelta param<real>  // CSD convergence criteria for the inner iteration
// $type CSDswitchStiffening param<int>  // 1: on, 2: off -> stiffening tangent stiffness matrix
//$type CSDexcitationType param<int> ; // 0: plunging, 1: flapping <- not used
// $type CSDflappingType param<int>  // 0: switch off, 1: sin, 2: cos, 3: 1-cos, 4: delay ( 1-exp(-4000*t^2))*sin
// $type CSDplungingType param<int>  // 0: switch off, 1: sin, 2: cos, 3: 1-cos
// $type CSDfrequency param<real>  // excitation frequency: ex. sin(2*pi*f*t) -> f (NOT omega!)
// $type CSDgenAlphaCoeff param<real>  // alpha for generalized alpha time integration scheme
// $type CSDnewmarkGammaCoeff param<real>  // Newmark beta coefficient1
// $type CSDnewmarkBetaCoeff param<real>  // Newmark beta coefficient2
// $type CSDdampingCoeff1 param<real>  // damping coefficient for the Rayleigh damping 1
// $type CSDdampingCoeff2 param<real>  // damping coefficient for the Rayleigh damping 2
// $type CSDplungeAmplitudeX param<real>  // plunge amplitude X
// $type CSDplungeAmplitudeY param<real>  // plunge amplitude Y
// $type CSDplungeAmplitudeZ param<real>  // plunge amplitude Z
// $type CSDflappingAmplitudeX param<real>  // flapping amplitude X <- pitch, different coordinate system in CSD, x <- spanwise
// $type CSDflappingAmplitudeY param<real>  // flapping amplitude Y <- flap, different coordinate system in CSD, y <- streamwise
// $type CSDflappingAmplitudeZ param<real>  // flapping amplitude Z <- lag, , different coordinate system in CSD, z <- -gravity 
// $type CSDMeshFilename param<string>  // CSD mesh filename
// $type CSDConnectivityFilename param<string>  // CSD connectivity filename
// $type CSDBCFilename param<string>  // CSD boundary conditions filename
// $type CSDstartingTimeStep param<int>  // 
// $type CSDtipNode param<int>  // Tip Node number inside CSD mesh
// $type CSDdimension param<int>  // CSD solver dimension constraint 2 -> 2D, 3 -> 3D
// $type CSD2dSpanCenter param<real>  // CSD solver dimension == 2D --> set the coordinate for the midspan point
// $type CSDYorigin param<real>  // CSD mesh assumed to be a plate with its Y (height) coordinate being this value. Then when interpolating to&from CFD the difference between the top and bottom CFD surfaces can be eliminated

// $type CFDMaxTotalInnerIterations param<int>  // Maximum allowed CFD iterations per time step

// $type FSICouplingMethod param<string>  // parameter to turn on the FSI coupling
// $type FSIRBFr param<real>  // FSI RBF interpolation r
// $type FSIRBFa param<real>  // FSI RBF interpolation a
// $type FSIRBFnr param<int>  // FSI RBF interpolation function nr
// $type FSIRBFMaxLinearSolverIterations param<int>  // FSI RBF interpolation number of maximum iterations for the linear solver
// $type FSIRBFTolerance param<real>  // FSI RBF interpolation tolerance for the linear solver
// $type FSICSD2CFDRBFr param<real>  // FSI CSD2CFD RBF interpolation r
// $type FSICSD2CFDRBFa param<real>  // FSI CSD2CFD RBF interpolation a
// $type FSICSD2CFDRBFnr param<int>  // FSI CSD2CFD RBF interpolation function nr
// $type FSICSD2CFDRBFMaxLinearSolverIterations param<int>  // FSI CSD2CFD RBF interpolation number of maximum iterations for the linear solver
// $type FSICSD2CFDRBFTolerance param<real>  // FSI RBF interpolation tolerance for the linear solver
// $type FSIIterationTolerance param<real>  // FSI inner iteration tolerance
// $type FSIIterationMinimum param<int>  // FSI minimum number of inner iterations 
// $type itfsi param<int>  // FSI inner iteration counter
// $type itfsi_ic param<int>  // FSI inner iteration counter initializer
// $type FSICoupling Constraint
// $type FSINLAMS Constraint
// $type FSIEULERBEAM Constraint
// $type FSI3DCONT Constraint
// $type FSIPRESCRIBED Constraint

  
// $type CSDNumberOfNodes param<int> 
// $type CSDNumberOfElements param<int> 
// $type CSDNumberOfBCs param<int> 
// $type CSDnodes_ic blackbox<ublas::matrix<real,ublas::column_major> >  // CSD initial 
// $type CSDnodes blackbox<ublas::matrix<real,ublas::column_major> >  // CSD current x
// $type CSDnodes_it blackbox<ublas::matrix<real,ublas::column_major> >  // CSD current x
// $type CSDnodesDisp_ic blackbox<ublas::matrix<real,ublas::column_major> >  // CSD initial generalized displacements
// $type CSDnodesDisp blackbox<ublas::matrix<real,ublas::column_major> >  // CSD current generalized displacements	
// $type CSDnodesVel_ic blackbox<ublas::matrix<real,ublas::column_major> >  // CSD initial generalized dot(displacements)
// $type CSDnodesVel blackbox<ublas::matrix<real,ublas::column_major> >  // CSD current generalized dot(displacements)
// $type CSDnodesAcc_ic blackbox<ublas::matrix<real,ublas::column_major> >  // CSD initial generalized ddot(displacements) 
// $type CSDnodesAcc blackbox<ublas::matrix<real,ublas::column_major> >  // CSD current generalized ddot(displacements)
// $type CSDConnectivity blackbox<ublas::matrix<int,ublas::column_major> >  // CSD connectivity matrix
// $type CSDBCdof param<vector<int> > 
// $type CSDBCZeroConstraint param<vector<real> > 

// $type CSDdisplacementsStar blackbox<ublas::matrix<real,ublas::column_major> >  // CSD displacements from NLAMS
// $type CSDdisplacementsOld blackbox<ublas::matrix<real,ublas::column_major> >  // CSD displacements from NLAMS, previous iteration
// $type CSDdisplacementsOldStar blackbox<ublas::matrix<real,ublas::column_major> >  // CSD displacements from NLAMS, previous iteration
// $type CSDdisplacementsOld_ic blackbox<ublas::matrix<real,ublas::column_major> >  // CSD displacements from NLAMS, previous iteration
// $type CSDnodesDispStar blackbox<ublas::matrix<real,ublas::column_major> >  // CSD generalized displacements from NLAMS
// $type CSDnodesVelStar blackbox<ublas::matrix<real,ublas::column_major> >  // CSD generalized dot(displacements) from NLAMS
// $type CSDnodesAccStar blackbox<ublas::matrix<real,ublas::column_major> >  // CSD generalized ddot(displacements) from NLAMS
// $type CSDForce blackbox<ublas::matrix<real,ublas::column_major> >  // Forces from CFD to CSD
// $type CSDForcePreStar blackbox<ublas::matrix<real,ublas::column_major> >  // CFD2CSD Force at the previous fsi iteration
// $type CSDForcePre blackbox<ublas::matrix<real,ublas::column_major> >  // CFD2CSD Force at the previous fsi iteration
// $type CSDForcePre_ic blackbox<ublas::matrix<real,ublas::column_major> >  // CFD2CSD Force at the previous fsi iteration initialization
// $type CSDnodesSysStar blackbox<boost::multi_array<real,4> >  // Nodal coordinates of each elements: from NLAMS
// $type CSDnodesSys blackbox<boost::multi_array<real,4> >  // Nodal coordinates of each elements
// $type CSDnodesSys_ic blackbox<boost::multi_array<real,4> >  // Nodal coordinates of each elements initialization
// $type CSDRBFweights blackbox<ublas::matrix<real,ublas::row_major> > 
	
// EULERBEAM
// $type CSDEulerXstart param<real> 
// $type CSDEulerXend param<real> 
// $type CSDEulerChord param<real> 
// $type CSDEulerXnum param<int> 
// $type CSDEulerAxis param<int>  // beam direction: 0->x, 1->y, 2->z, displacement probably always in y direction
// $type CSDEulerBeamDirection param<int> 
// $type CSDEulerSpanDirection param<int> 
// $type CSDxStar blackbox<ublas::vector<real> > 
// $type CSDxdotStar blackbox<ublas::vector<real> > 
// $type CSDxddotStar blackbox<ublas::vector<real> > 
// $type CSDx_ic blackbox<ublas::vector<real> > 
// $type CSDxdot_ic blackbox<ublas::vector<real> > 
// $type CSDxddot_ic blackbox<ublas::vector<real> > 
// $type CSDx blackbox<ublas::vector<real> > 
// $type CSDxdot blackbox<ublas::vector<real> > 
// $type CSDxddot blackbox<ublas::vector<real> > 
//$type CSDForcePre_ic blackbox<ublas::vector<real> > ;
//$type CSDForcePre blackbox<ublas::vector<real> > ;
//$type CSDForcePreStar blackbox<ublas::vector<real> > ;
#line 70 "FSI_CSD2CFD.loci"


// RBF
#include "rbf.h"

// delim
#include "delim.h"

namespace streamUns {


	// delimiter
	//delim = ", " ;
	
//-----------------------------------------------------------------------------
// Wrapper classes for PETSC objects. These have been created primarily to
// handle memory management.
  
  // Wrapper for Vec. Note that we cannot check the success of the destroy in
  // the destructor since the destructor cannot return a value.
  class PETSCCSD2CFDVector { // sequential vector
    private:
      bool created ;
      mutable Petsc::Vec v ;
    public:
      PETSCCSD2CFDVector() : created(false) {}
      ~PETSCCSD2CFDVector() { if(created) VecDestroy(v) ; }
    public:
      int AssemblyBegin() const {
        int ierr=VecAssemblyBegin(v) ; CHKERRQ(ierr) ; return 0 ;
      }
      int AssemblyEnd() const {
        int ierr=VecAssemblyEnd(v) ; CHKERRQ(ierr) ; return 0 ;
      }
      int Create(int localSize,int globalSize) {
        if(created) return 0 ;
        int ierr=VecCreateSeq(PETSC_COMM_SELF, localSize, &v) ; CHKERRQ(ierr) ;  
        //int ierr=VecCreate(PETSC_COMM_WORLD,&v) ; CHKERRQ(ierr) ;        
        //ierr=VecSetSizes(v,localSize,globalSize) ; CHKERRQ(ierr) ;
        //ierr=VecSetFromOptions(v) ; CHKERRQ(ierr) ;
        created=true ; return 0 ;
      }
	  int Destroy() const {
		  int ierr=VecDestroy(v) ; CHKERRQ(ierr) ;		  
		  return ierr ;
	  }
      const Petsc::Vec& Data() const { return v ; }
      int DuplicateFrom(const PETSCCSD2CFDVector &a) {
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
  class PETSCCSD2CFDMultiVector {
    private:
      bool created ;
      mutable Petsc::Vec vX,vY,vZ ;
    public:
      PETSCCSD2CFDMultiVector() : created(false) {}
      ~PETSCCSD2CFDMultiVector() {
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
        int ierr ;
        ierr=VecCreateSeq(PETSC_COMM_SELF, localSize, &vX) ; CHKERRQ(ierr) ; 	
        ierr=VecCreateSeq(PETSC_COMM_SELF, localSize, &vY) ; CHKERRQ(ierr) ; 	
        ierr=VecCreateSeq(PETSC_COMM_SELF, localSize, &vZ) ; CHKERRQ(ierr) ; 	
        //int ierr=VecCreate(PETSC_COMM_WORLD,&vX) ; CHKERRQ(ierr) ;
        //ierr=VecSetSizes(vX,localSize,globalSize) ; CHKERRQ(ierr) ;
        //ierr=VecSetFromOptions(vX) ; CHKERRQ(ierr) ;
        //ierr=VecCreate(PETSC_COMM_WORLD,&vY) ; CHKERRQ(ierr) ;
        //ierr=VecSetSizes(vY,localSize,globalSize) ; CHKERRQ(ierr) ;
        //ierr=VecSetFromOptions(vY) ; CHKERRQ(ierr) ;
        //ierr=VecCreate(PETSC_COMM_WORLD,&vZ) ; CHKERRQ(ierr) ;
        //ierr=VecSetSizes(vZ,localSize,globalSize) ; CHKERRQ(ierr) ;
       // ierr=VecSetFromOptions(vZ) ; CHKERRQ(ierr) ;
        created=true ; return 0 ;
      }
	  int Destroy() const {
		  int ierr=VecDestroy(vX) ; CHKERRQ(ierr) ;
		  ierr=VecDestroy(vY) ; CHKERRQ(ierr) ;
		  ierr=VecDestroy(vZ) ; CHKERRQ(ierr) ;
		  return ierr ;
	  }
      const Petsc::Vec& DataX() const { return vX ; }
      const Petsc::Vec& DataY() const { return vY ; }
      const Petsc::Vec& DataZ() const { return vZ ; }
      int DuplicateFrom(const PETSCCSD2CFDMultiVector &a) {
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
  class PETSCCSD2CFDMatrix {
    private:
      bool created ;
      mutable Petsc::Mat m ;
    public:
      PETSCCSD2CFDMatrix() : created(false) {}
      ~PETSCCSD2CFDMatrix() { if(created) MatDestroy(m) ; }
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
        //if(Loci::MPI_processes>1){
        //  int ierr=MatCreateMPIAIJ(PETSC_COMM_WORLD,numLocalRow,numLocalRow,
        //    numGlobalRow,numGlobalRow,0,numDiagonalNonZero,0,
        //    numOffDiagonalNonZero,&m) ; CHKERRQ(ierr) ;
        //}else{
          int ierr=MatCreateSeqAIJ(PETSC_COMM_SELF,numGlobalRow,numGlobalRow,0,numDiagonalNonZero,&m) ; CHKERRQ(ierr) ;
        //}
        created=true ; return 0 ;
      }
	  int Destroy() const {
		  int ierr=MatDestroy(m) ; CHKERRQ(ierr) ;
		  return ierr ;
	  }
      const Petsc::Mat& Data() const { return m ; }
	  int GetInfo() {
		  MatInfo info;
		  double mal, nz_alloc, nz_used, nz_un ;
		  MatGetInfo(m,MAT_LOCAL,&info);
		  mal = info.mallocs; nz_alloc=info.nz_allocated; nz_used=info.nz_used; nz_un=info.nz_unneeded ;
		  cout << "Displacement Interpolation CSD -> CFD: PETSC info: Number of mallocs during MatSetValues(): " << mal << delim << "Number of nonzeros (allocated, used, unneeded): " << nz_alloc << delim << nz_used << delim << nz_un << endl ;
		  return 0;
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
  class PETSCCSD2CFDKsp {
    private:
      bool created ;
      mutable Petsc::KSP ksp ;
    public:
      PETSCCSD2CFDKsp() : created(false) {}
      ~PETSCCSD2CFDKsp() { if(created) KSPDestroy(ksp) ; }
    public:
      int Create() {
        if(created) return 0 ;
        //int ierr=KSPCreate(PETSC_COMM_WORLD,&ksp) ; CHKERRQ(ierr) ;
        int ierr=KSPCreate(PETSC_COMM_SELF,&ksp) ; CHKERRQ(ierr) ;
        created=true ; return 0 ;
      }
	  int Destroy() const {
		  int ierr=KSPDestroy(ksp) ; CHKERRQ(ierr) ;
		  return ierr ;
	  }
	  int GetConvergedReason() const {
		  KSPConvergedReason reason ;
		  KSPGetConvergedReason(ksp, &reason) ;
		  if (Loci::MPI_rank==0) {
			  if (reason==KSP_CONVERGED_RTOL) {
				  cout << "CSD2CFD: PETSC solver: GetConvergedReason = " << reason << " KSP_CONVERGED_RTOL" << endl ;
			  } else if (reason==KSP_CONVERGED_ATOL) {
				  cout << "CSD2CFD: PETSC solver: GetConvergedReason = " << reason << " KSP_CONVERGED_ATOL" << endl ;
			  } else if (reason==KSP_CONVERGED_ITS) {
				  cout << "CSD2CFD: PETSC solver: GetConvergedReason = " << reason << " KSP_CONVERGED_ITS" << endl ;
			  } else if (reason==KSP_DIVERGED_NULL) {
				  cout << "CSD2CFD: PETSC solver: GetConvergedReason = " << reason << " KSP_DIVERGED_NULL" << endl ;
			  } else if (reason==KSP_DIVERGED_ITS) {
				  cout << "CSD2CFD: PETSC solver: GetConvergedReason = " << reason << " KSP_DIVERGED_ITS" << endl ;
			  } else if (reason==KSP_DIVERGED_DTOL) {
				  cout << "CSD2CFD: PETSC solver: GetConvergedReason = " << reason << " KSP_DIVERGED_DTOL" << endl ;
			  } else {
				  cout << "CSD2CFD: PETSC solver: GetConvergedReason = " << reason << endl ;
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
      int SetOperators(const PETSCCSD2CFDMatrix &m) const {
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
      int Solve(const PETSCCSD2CFDVector &rhs, PETSCCSD2CFDVector &sol) const {
//      int numIteration ;
        int ierr=KSPSolve(ksp,rhs.Data(),sol.Data()) ; CHKERRQ(ierr) ;
//KSPGetIterationNumber(ksp,&numIteration) ;
//cout << "petsc iterations: " << numIteration << endl ;
        return ierr ;
      }
      int SolveVector(const PETSCCSD2CFDMultiVector &rhs, PETSCCSD2CFDMultiVector &sol)
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


//$type FSICoupling Constraint ;
//$type bottomCSDnodes store<bool> ;
//$type topCSDnodes store<bool> ;

// $type CFDIterationFinished param<bool> 
// $type ncycle param<int> 

namespace {class file_FSI_CSD2CFD000_1280811401m0 : public Loci::blackbox_rule {
#line 366 "FSI_CSD2CFD.loci"
    Loci::const_param<real>  L_FSICSD2CFDRBFr_ ; 
#line 366 "FSI_CSD2CFD.loci"
    Loci::const_param<real>  L_FSICSD2CFDRBFa_ ; 
#line 366 "FSI_CSD2CFD.loci"
    Loci::const_param<int>  L_FSICSD2CFDRBFnr_ ; 
#line 366 "FSI_CSD2CFD.loci"
    Loci::const_param<int>  L_FSICSD2CFDRBFMaxLinearSolverIterations_ ; 
#line 366 "FSI_CSD2CFD.loci"
    Loci::const_param<real>  L_FSICSD2CFDRBFTolerance_ ; 
#line 366 "FSI_CSD2CFD.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 366 "FSI_CSD2CFD.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDdisplacementsStar_n__ ; 
#line 366 "FSI_CSD2CFD.loci"
    Loci::const_param<int>  L_ncycle_n__ ; 
#line 366 "FSI_CSD2CFD.loci"
    Loci::blackbox<ublas::matrix<real,ublas::row_major> >  L_CSDRBFweights_n__ ; 
#line 366 "FSI_CSD2CFD.loci"
public:
#line 366 "FSI_CSD2CFD.loci"
    file_FSI_CSD2CFD000_1280811401m0() {
#line 366 "FSI_CSD2CFD.loci"
       name_store("FSICSD2CFDRBFr",L_FSICSD2CFDRBFr_) ;
#line 366 "FSI_CSD2CFD.loci"
       name_store("FSICSD2CFDRBFa",L_FSICSD2CFDRBFa_) ;
#line 366 "FSI_CSD2CFD.loci"
       name_store("FSICSD2CFDRBFnr",L_FSICSD2CFDRBFnr_) ;
#line 366 "FSI_CSD2CFD.loci"
       name_store("FSICSD2CFDRBFMaxLinearSolverIterations",L_FSICSD2CFDRBFMaxLinearSolverIterations_) ;
#line 366 "FSI_CSD2CFD.loci"
       name_store("FSICSD2CFDRBFTolerance",L_FSICSD2CFDRBFTolerance_) ;
#line 366 "FSI_CSD2CFD.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 366 "FSI_CSD2CFD.loci"
       name_store("CSDdisplacementsStar{n}",L_CSDdisplacementsStar_n__) ;
#line 366 "FSI_CSD2CFD.loci"
       name_store("ncycle{n}",L_ncycle_n__) ;
#line 366 "FSI_CSD2CFD.loci"
       name_store("CSDRBFweights{n}",L_CSDRBFweights_n__) ;
#line 366 "FSI_CSD2CFD.loci"
       input("CSDnodes_ic,CSDdisplacementsStar{n},ncycle{n},FSICSD2CFDRBFr,FSICSD2CFDRBFa,FSICSD2CFDRBFnr,FSICSD2CFDRBFMaxLinearSolverIterations,FSICSD2CFDRBFTolerance") ;
#line 366 "FSI_CSD2CFD.loci"
       output("CSDRBFweights{n}") ;
#line 366 "FSI_CSD2CFD.loci"
       constraint("gridMotionTimeDependent{n}") ;
#line 366 "FSI_CSD2CFD.loci"
       disable_threading() ;
#line 366 "FSI_CSD2CFD.loci"
    }
#line 366 "FSI_CSD2CFD.loci"
    void prelude(const Loci::sequence &seq) { 
	
		if (Loci::MPI_rank==0) cout << "Displacement Interpolation CSD -> CFD: Starting" << endl ;
		
		if (Loci::MPI_rank==0) cout << "Displacement Interpolation CSD -> CFD: Starting, CFDIterationFinished" << endl ;
		int p = Loci::MPI_processes; int rank = Loci::MPI_rank ;
		PetscScalar valuex, valuey, valuez ;
		real distance=0. ;
		int columnIndex ;
		real r=(*L_FSICSD2CFDRBFr_), a=(*L_FSICSD2CFDRBFa_) ;
		int rbfNumber = (*L_FSICSD2CFDRBFnr_) ; // we get input as reals, convert to int
		int linearSolverMaxIteration = (*L_FSICSD2CFDRBFMaxLinearSolverIterations_) ; 
		real linearSolverTolerance = *L_FSICSD2CFDRBFTolerance_;
		int dimension = 3;
		int Nb = (*L_CSDnodes_ic_).size1() ;
		int Nbd1 = Nb + dimension +1 ;
		
	//	cout << "rank = " << rank << ", " << "Displacement Interpolation CSD -> CFD: Nb = " << Nb << ", Nbd1 = " << Nbd1 << endl ;
	//	cout << "rank = " << rank << ", " << "Displacement Interpolation CSD -> CFD: Assembling RHS" << endl ;
		// Assemble the RHS vectors
		PETSCCSD2CFDMultiVector FSI_RHS_vector, FSI_phi ; // declare RHS vector
		FSI_RHS_vector.Create(Nbd1, Nbd1) ; // create RHS vector
		for(int i=0;i<Nb;++i) {
			valuex=(*L_CSDdisplacementsStar_n__)(i,0);
			valuey=(*L_CSDdisplacementsStar_n__)(i,1);
			valuez=(*L_CSDdisplacementsStar_n__)(i,2);
			FSI_RHS_vector.SetValue(&i, &valuex, &valuey, &valuez) ;
		//	cout << "rank = " << rank << ", " << "Displacement Interpolation CSD -> CFD: RHS, i = " << i << ", valuex = " << valuex << ", valuey = " << valuey << ", valuez = " << valuez << ", CSDnodes.x = " << (*$CSDnodes_ic)(i,0) << ", CSDnodes.y = " << (*$CSDnodes_ic)(i,1) << ", CSDnodes.z = " << (*$CSDnodes_ic)(i,2) << endl ;
		}

		int tempTimeStepNumber = (*L_ncycle_n__) ;
	
	//	if (rank==0) {
//			std::stringstream filename ;	
//			filename << "CSD2CFDRHS" << std::setfill('0') << std::setw(5) << tempTimeStepNumber << ".dat" ;
//			ofstream CSD2CFDRHS ;
//			CSD2CFDRHS.open(filename.str().c_str(), ios::out) ;
////			CSD2CFDRHS << "CSDnodes.x" << ", " << "CSDnodes.y" << ", " << "CSDnodes.z" << ", " << "CSDdisp.x" << ", " << "CSDdisp.y" << ", " << "CSDdisp.z" << endl ;
//			for(int i=0; i<Nb;++i) { // in CFD coordinates
//				CSD2CFDRHS << (*$CSDnodes_ic)(i,0)  << ", " << (*$CSDnodes_ic)(i,1) << ", " << (*$CSDnodes_ic)(i,2) << ", " << (*$CSDdisplacementsStar{n})(i,0) << ", "  << (*$CSDdisplacementsStar{n})(i,1) << ", " << (*$CSDdisplacementsStar{n})(i,2)  << endl ;
//			}
//			CSD2CFDRHS.close();
	//	}
	
		// Determine the number of non-zero row-values
	//	cout << "rank = " << rank << ", " << "Displacement Interpolation CSD -> CFD: Determining non-zero row values" << endl ;
		vector<int> numberOfNonZeros(Nbd1,0) ;
		for(int i=0;i<Nb;++i) {
			for(int j=0;j<Nb;++j) {
				if (rbfNumber > 0 && rbfNumber < 9) {
					distance = pow((*L_CSDnodes_ic_)(i,0) -(*L_CSDnodes_ic_)(j,0), 2.) ;
					distance += pow((*L_CSDnodes_ic_)(i,1) - (*L_CSDnodes_ic_)(j,1), 2.) ;
					distance += pow((*L_CSDnodes_ic_)(i,2) - (*L_CSDnodes_ic_)(j,2), 2.) ;
					distance /= r * r ;
					if (distance <= 1.) {
						numberOfNonZeros[i]+=1;
					}
				} else {
					numberOfNonZeros[i]+=1;
				}
			}
		}
		numberOfNonZeros[Nbd1-4]=Nb;
		numberOfNonZeros[Nbd1-3]=Nb;
		numberOfNonZeros[Nbd1-2]=Nb;
		numberOfNonZeros[Nbd1-1]=Nb;
	
	
		// Setup the system matrix A
	//	cout << "rank = " << rank << ", " << "Displacement Interpolation CSD -> CFD: Setting up the system matrix" << endl ;
	  PETSCCSD2CFDMatrix FSI_A_matrix ;
		int temp=0;
		FSI_A_matrix.Create(Nbd1, Nbd1, &numberOfNonZeros[0], &temp) ;
		numberOfNonZeros.erase(numberOfNonZeros.begin(), numberOfNonZeros.end() ); // delete all elements of numberOfNonZeros
	//	cout << "rank = " << rank << ", " << "Displacement Interpolation CSD -> CFD: System matrix created.. Filling in the elements" << endl ;
			
		// Fill in the matrix A
		for(int i=0;i<Nb;++i) {
			for(int j=0;j<Nb;++j) {
				distance = pow((*L_CSDnodes_ic_)(i,0) - (*L_CSDnodes_ic_)(j,0), 2.) ;
				distance += pow((*L_CSDnodes_ic_)(i,1) - (*L_CSDnodes_ic_)(j,1), 2.) ;
				distance += pow((*L_CSDnodes_ic_)(i,2) - (*L_CSDnodes_ic_)(j,2), 2.) ;
				distance = sqrt(distance) ;
				if ( !( (rbfNumber > 0 && rbfNumber < 9) && (distance/r > 1.)) ) {
					valuex = radialBasisFunction(distance, r, a, rbfNumber) ;
					columnIndex=j; //columnValue=valuex ;
					FSI_A_matrix.SetRowValues(i,1,&columnIndex,&valuex);
				}
			}
			for(int j=0; j<4; ++j) {
				columnIndex=Nb+j; valuex=((j==0)?1.:(*L_CSDnodes_ic_)(i,j-1)); 
				FSI_A_matrix.SetRowValues(i,1,&columnIndex,&valuex) ;
				columnIndex=i;
				FSI_A_matrix.SetRowValues(Nb+j,1,&columnIndex,&valuex) ;
			}
		}
		if (rank==0) FSI_A_matrix.GetInfo() ;
		FSI_A_matrix.AssemblyBegin() ; FSI_A_matrix.AssemblyEnd() ;
	//	cout << "rank = " << rank << ", " << "Displacement Interpolation CSD -> CFD: Matrix Assembly done, solving" << endl ;
		MatView((FSI_A_matrix).Data(), PETSC_VIEWER_STDOUT_WORLD);
		
		// Setup the KSP solver
		PETSCCSD2CFDKsp FSI_ksp ;
		FSI_ksp.Create() ;
		FSI_ksp.SetTolerances(linearSolverTolerance,1.0e-30,linearSolverMaxIteration) ;
		PC pc ; FSI_ksp.GetPC(&pc) ; PCSetType(pc,PCJACOBI) ;  PCSetFromOptions(pc) ;
		//PC pc ; FSI_ksp.GetPC(&pc) ; PCSetType(pc,PCHYPRE) ; PCHYPRESetType(pc,"boomeramg") ; PCSetFromOptions(pc);
		FSI_ksp.SetRestart(1000000) ;
		FSI_ksp.SetInitialGuessNonzero() ;
		FSI_ksp.SetFromOptions() ;
	
	//	cout << "rank = " << rank << ", " << "Displacement Interpolation CSD -> CFD: Solving the system" << endl ;
		// Solve the system
		FSI_phi.Create(Nbd1, Nbd1) ;
		FSI_ksp.SetOperators(FSI_A_matrix) ; FSI_ksp.SolveVector(FSI_RHS_vector, FSI_phi) ;
		
	//	cout << "rank = " << rank << ", " << "Displacement Interpolation CSD -> CFD: system solved" << endl ;
	
		// Copy the solution to the local matrix
		(*L_CSDRBFweights_n__).resize(Nbd1, 3) ;
		PetscScalar *phixCopy, *phiyCopy, *phizCopy ;
		FSI_phi.GetArray(&phixCopy,&phiyCopy,&phizCopy) ;
		for(int i=0;i<Nbd1;++i) {
			(*L_CSDRBFweights_n__)(i,0) = phixCopy[i] ;
			(*L_CSDRBFweights_n__)(i,1) = phiyCopy[i] ;
			(*L_CSDRBFweights_n__)(i,2) = phizCopy[i] ;
		}
		FSI_ksp.GetConvergedReason() ;
		int its = FSI_ksp.GetIterationNumber() ;
		cout << "rank = " << rank << ", " << "CSD2CFD: solve -> number of iterations: " << its << endl ;
		cout << "rank = " << rank << ", " << "Displacement Interpolation CSD -> CFD: Complete" << endl ;
		//cout << "Displacement Interpolation CSD -> CFD: Complete, rank=" << rank << endl ;
//		FSI_ksp.Destroy();
//		FSI_phi.Destroy();
//		FSI_A_matrix.Destroy();
//		FSI_RHS_vector.Destroy() ;	
}    void compute(const Loci::sequence &seq) { 
#line 504 "FSI_CSD2CFD.loci"
      prelude(seq) ;
#line 504 "FSI_CSD2CFD.loci"
    }
#line 504 "FSI_CSD2CFD.loci"
} ;
#line 504 "FSI_CSD2CFD.loci"
Loci::register_rule<file_FSI_CSD2CFD000_1280811401m0> register_file_FSI_CSD2CFD000_1280811401m0 ;
#line 504 "FSI_CSD2CFD.loci"
}
#line 504 "FSI_CSD2CFD.loci"
// $type pos store<vect3d> 
// $type node_s_b_flex store<vect3d> 

namespace {class file_FSI_CSD2CFD001_1280811401m4 : public Loci::pointwise_rule {
#line 509 "FSI_CSD2CFD.loci"
    Loci::const_param<int>  L_CSDdimension_ ; 
#line 509 "FSI_CSD2CFD.loci"
    Loci::const_param<real>  L_CSD2dSpanCenter_ ; 
#line 509 "FSI_CSD2CFD.loci"
    Loci::const_param<real>  L_CSDYorigin_ ; 
#line 509 "FSI_CSD2CFD.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 509 "FSI_CSD2CFD.loci"
    Loci::const_param<int>  L_ncycle_n__ ; 
#line 509 "FSI_CSD2CFD.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::row_major> >  L_CSDRBFweights_n__ ; 
#line 509 "FSI_CSD2CFD.loci"
    Loci::const_store<vect3d>  L_pos_ ; 
#line 509 "FSI_CSD2CFD.loci"
    Loci::const_param<real>  L_FSICSD2CFDRBFr_n__ ; 
#line 509 "FSI_CSD2CFD.loci"
    Loci::const_param<real>  L_FSICSD2CFDRBFa_n__ ; 
#line 509 "FSI_CSD2CFD.loci"
    Loci::const_param<int>  L_FSICSD2CFDRBFnr_n__ ; 
#line 509 "FSI_CSD2CFD.loci"
    Loci::store<vect3d>  L_node_s_b_flex_n__ ; 
#line 509 "FSI_CSD2CFD.loci"
public:
#line 509 "FSI_CSD2CFD.loci"
    file_FSI_CSD2CFD001_1280811401m4() {
#line 509 "FSI_CSD2CFD.loci"
       name_store("CSDdimension",L_CSDdimension_) ;
#line 509 "FSI_CSD2CFD.loci"
       name_store("CSD2dSpanCenter",L_CSD2dSpanCenter_) ;
#line 509 "FSI_CSD2CFD.loci"
       name_store("CSDYorigin",L_CSDYorigin_) ;
#line 509 "FSI_CSD2CFD.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 509 "FSI_CSD2CFD.loci"
       name_store("ncycle{n}",L_ncycle_n__) ;
#line 509 "FSI_CSD2CFD.loci"
       name_store("CSDRBFweights{n}",L_CSDRBFweights_n__) ;
#line 509 "FSI_CSD2CFD.loci"
       name_store("pos",L_pos_) ;
#line 509 "FSI_CSD2CFD.loci"
       name_store("FSICSD2CFDRBFr{n}",L_FSICSD2CFDRBFr_n__) ;
#line 509 "FSI_CSD2CFD.loci"
       name_store("FSICSD2CFDRBFa{n}",L_FSICSD2CFDRBFa_n__) ;
#line 509 "FSI_CSD2CFD.loci"
       name_store("FSICSD2CFDRBFnr{n}",L_FSICSD2CFDRBFnr_n__) ;
#line 509 "FSI_CSD2CFD.loci"
       name_store("node_s_b_flex{n}",L_node_s_b_flex_n__) ;
#line 509 "FSI_CSD2CFD.loci"
       input("pos,CSDnodes_ic,CSDYorigin,ncycle{n},CSD2dSpanCenter,CSDdimension,CSDRBFweights{n},FSICSD2CFDRBFr{n},FSICSD2CFDRBFa{n},FSICSD2CFDRBFnr{n}") ;
#line 509 "FSI_CSD2CFD.loci"
       output("node_s_b_flex{n}") ;
#line 509 "FSI_CSD2CFD.loci"
       constraint("gridMotionTimeDependent{n},FlexibleBoundaryDisplacement{n}") ;
#line 509 "FSI_CSD2CFD.loci"
       disable_threading() ;
#line 509 "FSI_CSD2CFD.loci"
    }
#line 509 "FSI_CSD2CFD.loci"
    void calculate(Entity _e_) { 
#line 510 "FSI_CSD2CFD.loci"

	//if ($CFDIterationFinished{n,it-1}) {

		//int p = Loci::MPI_processes; 
		int rank = Loci ::MPI_rank ;
		real distance =0. ;
		real r =(L_FSICSD2CFDRBFr_n__[_e_]), a =(L_FSICSD2CFDRBFa_n__[_e_]) ;
		int rbfNumber = (L_FSICSD2CFDRBFnr_n__[_e_]) ; // we get input as reals, convert to int
		int d = 3;
		int Nb = (L_CSDnodes_ic_[_e_]).size1 () ;
		//int Nbd1 = Nb + d + 1;	
		
		real posy = (L_CSDYorigin_[_e_]) ;
		real posx = (L_pos_[_e_]).x ;
		//real posy = ($pos).y  ;
		real posz = ((L_CSDdimension_[_e_]==2)? (L_CSD2dSpanCenter_[_e_]) : (L_pos_[_e_]).z )  ;
				
		//cout << "Displacement Interpolation CSD -> CFD: Start applying, r=" << rank << endl ;
		vector <real > value (Nb ,0.);  
		  for (int i =0;i <Nb ;++i ) {
			  distance = pow (posx -  (L_CSDnodes_ic_[_e_])(i ,0), 2.) ; // 
			  distance += pow (posy - (L_CSDnodes_ic_[_e_])(i ,1), 2.) ; // removing the y-component between the CFD and the CSD 
			  distance += pow (posz - (L_CSDnodes_ic_[_e_])(i ,2), 2.) ; // 
			  distance = sqrt (distance ) ;
			  value [i ] = radialBasisFunction (distance , r , a , rbfNumber ) ;
		  }
			for (int di =0;di <d ;++di ) {
				double polynom = (L_CSDRBFweights_n__[_e_])(Nb , di ) ;
				polynom += (L_CSDRBFweights_n__[_e_])(Nb +1,di ) * posx ;
				polynom += (L_CSDRBFweights_n__[_e_])(Nb +2,di ) * posy ;				
				polynom += (L_CSDRBFweights_n__[_e_])(Nb +3,di ) * posz ;
				if (di ==0) {
					(L_node_s_b_flex_n__[_e_]).x = polynom ;
				} else if (di ==1) {
					(L_node_s_b_flex_n__[_e_]).y = polynom ;
				} else if (di ==2) {
					(L_node_s_b_flex_n__[_e_]).z = polynom ;
				} else {
					cout << "We shouldn't be here " << endl ;
				}
			}
			
			// RBF contribution
			for (int i =0;i <Nb ;++i ){		
				(L_node_s_b_flex_n__[_e_]).x += (L_CSDRBFweights_n__[_e_])(i ,0) * value [i ] ;
				(L_node_s_b_flex_n__[_e_]).y += (L_CSDRBFweights_n__[_e_])(i ,1) * value [i ] ;
				(L_node_s_b_flex_n__[_e_]).z += (L_CSDRBFweights_n__[_e_])(i ,2) * value [i ] ;
			}
			value .erase (value .begin (), value .end ()) ; 

			int tempTimeStepNumber = (L_ncycle_n__[_e_]) ;
			std ::stringstream filename ;	
//			filename << "CSD2CFDSOL" << std::setfill('0') << std::setw(5) << tempTimeStepNumber << ".dat" ;
//			ofstream CSD2CFDRHS ;
//			CSD2CFDRHS.open(filename.str().c_str(), ios::app) ;
////			CSD2CFDRHS << "CSDnodes.x" << ", " << "CSDnodes.y" << ", " << "CSDnodes.z" << ", " << "CSDdisp.x" << ", " << "CSDdisp.y" << ", " << "CSDdisp.z" << endl ;
////			for(int i=0; i<Nb;++i) { // in CFD coordinates
//				CSD2CFDRHS << posx  << ", " << posy << ", " << posz << ", " << ($node_s_b_flex{n}).x << ", "  << ($node_s_b_flex{n}).y << ", " << ($node_s_b_flex{n}).z  << endl ;
////			}
//			CSD2CFDRHS.close();

	//		cout << "rank = " << rank << ", " << "Displacement Interpolation CSD -> CFD: x,y,z=" << $pos.x << ", " << $pos.y << ", " << $pos.z << ", node_s_b_flex = " << ($node_s_b_flex{n}).x << ", " << ($node_s_b_flex{n}).y << ", " << ($node_s_b_flex{n}).z << ", " << endl ;
	//		cout << "CSD2CFD: Displacement Interpolation CSD -> CFD: r,x,y,z=" << rank << ", " << $pos.x << ", " << $pos.y << ", " << $pos.z << ", node_s_b_flex = " << ($node_s_b_flex{n}).x << ", " << ($node_s_b_flex{n}).y << ", " << ($node_s_b_flex{n}).z << ", " << endl ;
	//		cout << "rank = " << rank << ", " << "Displacement Interpolation CSD -> CFD: End applying" << endl ;
	//	} else {
	//		($node_s_b_flex{n}) = vect3d(0.0, 0.0, 0.0) ;
	//	}
	}    void compute(const Loci::sequence &seq) { 
#line 577 "FSI_CSD2CFD.loci"
      do_loop(seq,this) ;
#line 577 "FSI_CSD2CFD.loci"
    }
#line 577 "FSI_CSD2CFD.loci"
} ;
#line 577 "FSI_CSD2CFD.loci"
Loci::register_rule<file_FSI_CSD2CFD001_1280811401m4> register_file_FSI_CSD2CFD001_1280811401m4 ;
#line 577 "FSI_CSD2CFD.loci"
}
#line 577 "FSI_CSD2CFD.loci"
 
}





