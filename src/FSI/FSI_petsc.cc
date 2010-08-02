//-----------------------------------------------------------------------------
// Description: This file contains rules for implementing a FSI linear solver
//
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

// StreamUns includes.
#include "sciTypes.h"

// FSI includes.
#include "fsi.h"

// Include PETSC wrappers
#include "FSI_petsc.h"
  
//-----------------------------Nonzero entries---------------------------------------------------------
//
  // Determines the number of non-zero entries in the local portion of the
  // Petsc matrix for each node. The local portion is defined as the square
  // row/column sub-block of the PETSC matrix whose rows map to nodes local
  // to the process.
  class fsiBoundaryFaceNumDiagonalNonZeroEntries : public pointwise_rule {
    private:
	  const_param<real> r ;
	  const_param<int> fsiNumber ;	  
	  const_param<vector<real> > Q ;
      const_param<vector<int> > fsiNumBoundaryFace ;
      store<int> fsiBoundaryFaceNumDiagonalNonZero ;
	  store<int> fsiBoundaryFaceNumOffDiagonalNonZero ;
    public:

      // Define input and output.
      fsiBoundaryFaceNumDiagonalNonZeroEntries() {
        name_store("fsiNumBoundaryFace(X)",fsiNumBoundaryFace) ;
        name_store("fsiBoundaryFaceNumDiagonalNonZero(X)",fsiBoundaryFaceNumDiagonalNonZero) ; 
		name_store("fsiBoundaryFaceNumOffDiagonalNonZero(X)",fsiBoundaryFaceNumOffDiagonalNonZero) ;
		name_store("fsiR",r) ;
		name_store("fsiNumber",fsiNumber) ;
		name_store("fsiBoundaryFaceQ(X)",Q) ;
		input("fsiR, fsiNumber") ;
		input("fsiBoundaryFaceQ(X)") ;
        output("fsiBoundaryFaceNumDiagonalNonZero(X)") ;
		output("fsiBoundaryFaceNumOffDiagonalNonZero(X)") ;
        constraint("X") ;
		disable_threading() ;
      }

      // Loop over nodes.
      void compute(const sequence & seq) { 
		
		int rank=Loci::MPI_rank ;
		
		// Get the number of local and global nodes.
        int globalNumFace=0, localStart=0  ;
        for(unsigned int i=0;i<(*fsiNumBoundaryFace).size();++i) globalNumFace+=(*fsiNumBoundaryFace)[i] ;
        for(int i=0;i<rank;++i) localStart+=(*fsiNumBoundaryFace)[i]  ;  // 3 coordinates
		int localNum=(*fsiNumBoundaryFace)[rank]; 
		
		double distance = 0.;
		int localEnd = localStart+localNum ;
		sequence::const_iterator nodePtr;
//		for(int i=localStart;i<localEnd;++i){	// Q[row,col]=Q[3*row+col], Q=[0,1,2;3,4,5;...]
		int counter=localStart ;
		int I,I1,I2,J,J1,J2;
		int counterD=0; int counterOD=0;
		for(nodePtr=seq.begin();nodePtr!=seq.end();++nodePtr,++counter){	// Q[row,col]=Q[3*row+col], Q=[0,1,2;3,4,5;...]
			I=3*counter; I1=3*counter+1; I2=3*counter+2 ;
			for(int j=0;j<globalNumFace;++j){
				J=3*j; J1=3*j+1; J2=3*j+2 ;
				fsiBoundaryFaceNumDiagonalNonZero[*nodePtr]=0;
				fsiBoundaryFaceNumOffDiagonalNonZero[*nodePtr]=0;
				if( (j>=localStart) && (j<localEnd) ) { // nodePtr2 is in the same partition as nodePtr1 -> diagonal
					if (*fsiNumber > 0 && *fsiNumber < 9) { // compact
						distance = pow((*Q)[I] - (*Q)[J], 2.) ; // ((*Q)[I] - (*Q)[J])*((*Q)[I] - (*Q)[J]) ;
						distance += pow((*Q)[I1] - (*Q)[J1], 2.) ; // ((*Q)[I1] - (*Q)[J1])*((*Q)[I1] - (*Q)[J1]) ;
						distance += pow((*Q)[I2] - (*Q)[J2], 2.) ; // ((*Q)[I2] - (*Q)[J2])*((*Q)[I2] - (*Q)[J2]) ;
						distance /= ((*r)*(*r)) ;
						if (distance <=1.) {
							++counterD;
						} 
					} else { // global
						++counterD;
					}
				} else { // nodePtr2 is not in the same partition as nodePtr2 -> off diagonal
					if (*fsiNumber > 0 && *fsiNumber < 9) { // compact
						distance = pow((*Q)[I] - (*Q)[J], 2.) ; // ((*Q)[I] - (*Q)[J])*((*Q)[I] - (*Q)[J]) ;
						distance += pow((*Q)[I1] - (*Q)[J1], 2.) ; // ((*Q)[I1] - (*Q)[J1])*((*Q)[I1] - (*Q)[J1]) ;
						distance += pow((*Q)[I2] - (*Q)[J2], 2.) ; // ((*Q)[I2] - (*Q)[J2])*((*Q)[I2] - (*Q)[J2]) ;
						distance /= ((*r)*(*r)) ;
						if (distance <=1) {
							++counterOD;
						}
					} else { // global
						++counterOD;
					}
				}
			}			
			fsiBoundaryFaceNumDiagonalNonZero[*nodePtr]=counterD;
			fsiBoundaryFaceNumOffDiagonalNonZero[*nodePtr]=counterOD;
			counterD=0; counterOD=0; // reset
		}
	  }
  } ;

  register_rule<fsiBoundaryFaceNumDiagonalNonZeroEntries>
    registerfsiBoundaryFaceNumDiagonalNonZeroEntries ;

//-----------------------------Setup linear solver---------------------------------------------------------
//
  // Sets up the PETSC linear solver.
  class PETSCBoundaryFaceSetupSolver : public singleton_rule {
    private:
      const_param<int> maxLinearSolverIterations ;
      const_param<real> linearSolverTolerance ;
      blackbox<PETSCKsp> ksp ;
    public:

      // Define input and output.
      PETSCBoundaryFaceSetupSolver() {
        name_store("fsi_maxLinearSolverIterations",maxLinearSolverIterations) ;
        name_store("fsi_linearSolverTolerance",linearSolverTolerance) ;
        name_store("petscBoundaryFaceKSP(X)",ksp) ;
        input("fsi_maxLinearSolverIterations") ;
        input("fsi_linearSolverTolerance") ;
        output("petscBoundaryFaceKSP(X)") ;
        constraint("fsi_PETSCLinearSolver") ;
		constraint("X") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {
        (*ksp).Create() ;
        (*ksp).SetTolerances(*linearSolverTolerance,1.0e-30,*maxLinearSolverIterations) ;
      // Valve case.
      //PC pc ; (*ksp).GetPC(&pc) ; PCSetType(pc,PCILU) ; PCILUSetLevels(pc,3) ;
      // Injector case.
//      PC pc ; (*ksp).GetPC(&pc) ; PCSetType(pc,KSPCG); // PCSetType(pc,PCJACOBI) ;
      PC pc ; (*ksp).GetPC(&pc) ; PCSetType(pc,PCJACOBI); // PCSetType(pc,PCJACOBI) ;
	  (*ksp).SetInitialGuessNonzero() ; // Set initial guess using the values
      // PCILUSetLevels(pc,0) ;
        (*ksp).SetFromOptions() ;
      }
  } ;

  register_rule<PETSCBoundaryFaceSetupSolver> registerPETSCBoundaryFaceSetupSolver ;

//---------------------------------------------b-------------------------------------------------
//
//
  // Sets up the PETSC right-hand-side vector.
  class PETSCBoundaryFaceUnitRHS : public unit_rule {
    private:
      const_param<vector<int> > fsiNumBoundaryFace ;
      const_store<vect3d> B ;
      blackbox<PETSCMultiVector> b;
    public:

      // Define input and output.
      PETSCBoundaryFaceUnitRHS() {
        name_store("fsiNumBoundaryFace(X)",fsiNumBoundaryFace) ;
        name_store("petscBoundaryFaceB(X)",b) ;
        name_store("FSI_B",B);
        input("fsiNumBoundaryFace(X)") ;
				input("FSI_B");
        output("petscBoundaryFaceB") ;
        constraint("fsi_PETSCLinearSolver") ;
				constraint("X") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {

        // Get the number of local and global nodes.
        int localNumFace=(*fsiNumBoundaryFace)[Loci::MPI_rank],globalNumFace=0 ;
        for(unsigned int i=0;i<(*fsiNumBoundaryFace).size();++i) globalNumFace+=(*fsiNumBoundaryFace)[i] ;
		
		// include the polynomial coeff. computation
		if (Loci::MPI_rank == Loci::MPI_processes-1) { // last rank takes care of beta
			localNumFace += 4; // d+1 = 4
		}
					
		globalNumFace+=4; // for beta (d+1)
		
        // Allocate the unknown and rhs vectors.
		(*b).Create(localNumFace,globalNumFace) ;

		// Compute the row offset for this process.
        int localStart=0 ;
        for(int i=0;i<Loci::MPI_rank;++i){ localStart+=(*fsiNumBoundaryFace)[i] ; }
		//
					
		sequence::const_iterator nodePtr;
		PetscScalar valuex, valuey, valuez;

//		cout << "rank, localStart, localNumFace, globalNumFace: " << Loci::MPI_rank << ", " << localStart << ", " << localNumFace << ", " << globalNumFace << endl;
		
		int counter=localStart ;
		for(nodePtr=seq.begin();nodePtr!=seq.end();++nodePtr,++counter){			
				// RHS vector
				valuex=B[*nodePtr].x;
				valuey=B[*nodePtr].y;
				valuez=B[*nodePtr].z;
//				cout << "rank, i,B.x,valuex,valuey " << Loci::MPI_rank << ", " << counter << ", " << B[*nodePtr].x << ", " << valuex << delim << B[*nodePtr].y << delim << valuey << endl;
				(*b).SetValue(&counter,&valuex, &valuey, &valuez) ;
		}		

		// assemble
		(*b).AssemblyBegin(); (*b).AssemblyEnd() ;
      }
  } ;
 
  register_rule<PETSCBoundaryFaceUnitRHS> registerPETSCBoundaryFaceUnitRHS ;



  // Sets up the PETSC right-hand-side vector.
  class PETSCBoundaryFaceApplyRHS : public apply_rule<blackbox<PETSCVector>,Loci::NullOp<PETSCVector> > {
    private:
      const_param<vector<int> > fsiNumBoundaryFace ;
      const_store<vect3d> B ;
      blackbox<PETSCMultiVector> b;
    public:

      // Define input and output.
      PETSCBoundaryFaceApplyRHS() {
        name_store("fsiNumBoundaryFace(X)",fsiNumBoundaryFace) ;
        name_store("petscBoundaryFaceB(X)",b) ;
        name_store("FSI_B",B);
        input("fsiNumBoundaryFace(X)") ;
		input("FSI_B(X)");
        output("petscBoundaryFaceB(X)") ;
        constraint("fsi_PETSCLinearSolver") ;
		constraint("X") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) { }
  } ;
 
  register_rule<PETSCBoundaryFaceApplyRHS> registerPETSCBoundaryFaceApplyRHS ;


//---------------------------------------------A-------------------------------------------------
//
  
  // Sets up the PETSC matrix. Note that this is a unit_rule, since this is the
  // only way we can have stores as input to a rule outputting blackboxes.
  class PETSCBoundaryFaceSetupMatrixUnit : public unit_rule {
    private:
	  const_param<real> r,a ;
	  const_param<int> fsiNr ;	  
//      const_store<int> fsiBoundaryFaceToRow ;
      const_store<int> fsiBoundaryFaceNumDiagonalNonZero ;
      const_store<int> fsiBoundaryFaceNumOffDiagonalNonZero ;
      const_param<vector<int> > fsiNumBoundaryFace ;
      const_param<vector<real> > Q ;
      blackbox<PETSCMatrix> A ;
	private:
      int columnIndex ;
      PetscScalar columnValue ;
    public:

      // Define input and output.
      PETSCBoundaryFaceSetupMatrixUnit() {
        name_store("fsiBoundaryFaceNumDiagonalNonZero(X)",fsiBoundaryFaceNumDiagonalNonZero) ;
        name_store("fsiBoundaryFaceNumOffDiagonalNonZero(X)",fsiBoundaryFaceNumOffDiagonalNonZero) ;
        name_store("fsiNumBoundaryFace(X)",fsiNumBoundaryFace) ;
        name_store("petscBoundaryFaceA(X)",A) ;		
		name_store("fsiBoundaryFaceQ(X)",Q) ;
		name_store("fsiR",r) ;
		name_store("fsiA",a) ;
		name_store("fsiNumber",fsiNr) ;
		input("fsiNumber, fsiR, fsiA");
        input("fsiBoundaryFaceQ(X)") ;
        input("fsiBoundaryFaceNumDiagonalNonZero(X),fsiBoundaryFaceNumOffDiagonalNonZero(X)") ;
        input("fsiNumBoundaryFace(X)") ;
        output("petscBoundaryFaceA(X)") ;
		constraint("X") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {

		int rank=Loci::MPI_rank ;
		int p=Loci::MPI_processes ;

        // Get the number of local and global nodes.
        int localNumFace=(*fsiNumBoundaryFace)[rank],globalNumFace=0 ;
        for(unsigned int i=0;i<(*fsiNumBoundaryFace).size();++i) globalNumFace+=(*fsiNumBoundaryFace)[i] ;
		  		
		// include the polynomial coeff. computation
		if (rank == p-1) { // last rank takes care of beta
			localNumFace += 4; // d+1 = 4
		}		
		globalNumFace+=4; // for beta (d+1)

		// Compute the row offset for this process.
        int localStart=0 ; 
        for(int i=0;i<rank;++i) localStart+=(*fsiNumBoundaryFace)[i]  ;  // 3 coordinates

//		cout << "Start Allocation A" << endl ;
        // Make temporary copy of matrix allocation data.
        int count=0,*numDiagonalNonZero=new int[localNumFace], *numOffDiagonalNonZero=new int[localNumFace] ;
		sequence::const_iterator nodePtr ;
        for(nodePtr=seq.begin();nodePtr!=seq.end();++nodePtr,++count){
          numDiagonalNonZero[count]=fsiBoundaryFaceNumDiagonalNonZero[*nodePtr] + ((rank==p-1)?4:0) ; // + Q, only r==p1 gets Q in diagNonZero, otherwise in offDiagNonZero
          numOffDiagonalNonZero[count]=fsiBoundaryFaceNumOffDiagonalNonZero[*nodePtr] + ((rank==p-1)?0:4) ;
        }
		
		// Account for the last four rows: Q
		if(rank==p-1) {
			for(int i=0;i<4;++i) {
				numOffDiagonalNonZero[count+i]=localStart; // globalNumFace - 4 (zeros) = # boundary Faces
				numDiagonalNonZero[count+i]=globalNumFace-localStart-4;
			}
		}
		for(int i=0; i<localNumFace; ++i) {
		}
		int sumDiag=0, sumOffDiag=0;
		for(int i=0;i<localNumFace;++i) {sumDiag+=numDiagonalNonZero[i]; sumOffDiag+=numOffDiagonalNonZero[i]; }
		cout << "p,r,Diag,OffDiag" << p << delim << rank << delim << sumDiag << delim << sumOffDiag << endl ;
        // Allocate the matrix.
        (*A).Create(localNumFace,globalNumFace,numDiagonalNonZero,numOffDiagonalNonZero) ;

        // Deallocate temporary copy of matrix allocation data.
        delete [] numDiagonalNonZero ; delete [] numOffDiagonalNonZero ;

		
		double distance = 0.;
		int localEnd = localStart+localNumFace ;
		double valueQ[4] ;
		PetscScalar value=0.0 ;
		int counter=localStart ;
		int I,I1,I2,J,J1,J2;
		for(nodePtr=seq.begin();nodePtr!=seq.end();++nodePtr,++counter){	// Q[row,col]=Q[3*row+col], Q=[0,1,2;3,4,5;...]
			I=3*counter; I1=3*counter+1; I2=3*counter+2 ;
			for(int j=0;j<globalNumFace-4;++j){
				J=3*j; J1=3*j+1; J2=3*j+2 ;
				distance = pow((*Q)[I] - (*Q)[J], 2.) ; //)*((*Q)[I] - (*Q)[J]) ;
				distance += pow((*Q)[I1] - (*Q)[J1], 2.) ; //*((*Q)[I1] - (*Q)[J1]) ;
				distance += pow((*Q)[I2] - (*Q)[J2], 2.) ; //)*((*Q)[I2] - (*Q)[J2]) ;
				distance = sqrt(distance) ;
				if (!( (*fsiNr > 0 && *fsiNr < 9) && (distance/(*r) > 1.))) { // compact
					value = radialBasisFunction(distance, *r, *a, *fsiNr) ;				
					columnIndex=j; columnValue=value;
					(*A).SetRowValues(counter,1,&columnIndex,&columnValue) ; // i:row, j:column
				}
			}			
			valueQ[0] = 1.0; valueQ[1]=(*Q)[I]; valueQ[2]=(*Q)[I1]; valueQ[3]=(*Q)[I2];
			for(int j=0; j<4; ++j) {
				columnIndex=globalNumFace-4+j; columnValue=valueQ[j];
				(*A).SetRowValues(counter,1,&columnIndex,&columnValue) ; // i:row, j:column
				columnIndex=counter; 
				(*A).SetRowValues(globalNumFace-4+j,1,&columnIndex, &columnValue); // Q^T
			}			
		}
		(*A).GetInfo() ;
		(*A).AssemblyBegin() ; (*A).AssemblyEnd() ;
      }
  } ;

  register_rule<PETSCBoundaryFaceSetupMatrixUnit> registerPETSCBoundaryFaceSetupMatrixUnit ;


  // Assemble PETSC-A- system.
  class PETSCBoundaryFaceSetupMatrixApply : public apply_rule<blackbox<PETSCMatrix>,Loci::NullOp<PETSCMatrix> > {
    private:
	  const_param<real> r,a ;
	  const_param<int> fsiNr ;	  
      const_store<int> fsiBoundaryFaceNumDiagonalNonZero ;
      const_store<int> fsiBoundaryFaceNumOffDiagonalNonZero ;
      const_param<vector<int> > fsiNumBoundaryFace ;
      const_param<vector<real> > Q ;
      blackbox<PETSCMatrix> A ;
	private:
//      int *columnIndex ;
//      PetscScalar *columnValue ;
    public:

      // Define input and output.
      PETSCBoundaryFaceSetupMatrixApply() {
        name_store("fsiBoundaryFaceNumDiagonalNonZero(X)",fsiBoundaryFaceNumDiagonalNonZero) ;
        name_store("fsiBoundaryFaceNumOffDiagonalNonZero(X)",fsiBoundaryFaceNumOffDiagonalNonZero) ;
        name_store("fsiNumBoundaryFace(X)",fsiNumBoundaryFace) ;
        name_store("petscBoundaryFaceA(X)",A) ;		
		name_store("fsiBoundaryFaceQ(X)",Q) ;
		name_store("fsiR",r) ;
		name_store("fsiA",a) ;
		name_store("fsiNumber",fsiNr) ;
		input("fsiNumber, fsiR, fsiA");
//        input("fsiBoundaryFaceToRow") ;
        input("fsiBoundaryFaceQ(X)") ;
        input("fsiBoundaryFaceNumDiagonalNonZero(X),fsiBoundaryFaceNumOffDiagonalNonZero(X)") ;
        input("fsiNumBoundaryFace(X)") ;
        output("petscBoundaryFaceA(X)") ;
		constraint("X") ;
        disable_threading() ;
      }

      // Assemble and solve.
	  virtual void compute(const sequence &seq) { }
  } ;

  register_rule<PETSCBoundaryFaceSetupMatrixApply> registerPETSCBoundaryFaceSetupMatrixApply ;

//---------------------------------------------solution-------------------------------------------------
//
  
  // Sets up the PETSC solution vector.
  class PETSCBoundaryFaceSetupPhi : public unit_rule {
    private:
      const_param<vector<int> > fsiNumBoundaryFace ;
      const_blackbox<PETSCMultiVector> b ;
      const_blackbox<PETSCMatrix> A ;
      const_blackbox<PETSCKsp> ksp ;
      blackbox<PETSCMultiVector> phi;
    public:

      // Define input and output.
      PETSCBoundaryFaceSetupPhi() {
        name_store("fsiNumBoundaryFace(X)",fsiNumBoundaryFace) ;
        name_store("petscBoundaryFacePhi(X)",phi) ;
        name_store("petscBoundaryFaceA(X)",A) ;
        name_store("petscBoundaryFaceB(X)",b) ;
        name_store("petscBoundaryFaceKSP(X)",ksp) ;
        input("petscBoundaryFaceA(X),petscBoundaryFaceB(X),petscBoundaryFaceKSP(X)") ;
        input("fsiNumBoundaryFace(X)") ;
        output("petscBoundaryFacePhi(X)") ;
        constraint("X");
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {

        // Get the number of local and global nodes.
        int localNumFace=(*fsiNumBoundaryFace)[Loci::MPI_rank],globalNumFace=0 ;
        for(unsigned int i=0;i<(*fsiNumBoundaryFace).size();++i) globalNumFace+=(*fsiNumBoundaryFace)[i] ;
		
		// include the polynomial coeff. computation
		if (Loci::MPI_rank == Loci::MPI_processes-1) { // last rank takes care of beta
			localNumFace += 4; // d+1 = 4
		}
					
		globalNumFace+=4; // for beta (d+1)
		
        // Allocate the unknown and rhs vectors.
		(*phi).Create(localNumFace,globalNumFace) ;
		
		double starttime, endtime;
		// Solve the linear system
		cout << "p,r: Start Solving" << Loci::MPI_processes << delim << Loci::MPI_rank << endl ;
		starttime = MPI_Wtime();
		(*ksp).SetOperators(*A) ; (*ksp).SolveVector(*b,*phi) ;
		endtime = MPI_Wtime();
		cout << "p,r: " << Loci::MPI_processes << delim << Loci::MPI_rank << delim << " solved in " << endtime - starttime << " sec" << endl ;
		int its = (*ksp).GetIterationNumber() ;
		cout << "Petsc solver converged in " << its << " iterations" << endl ;
      }
  } ;
 
  register_rule<PETSCBoundaryFaceSetupPhi> registerPETSCBoundaryFaceSetupPhi ;


  // Assemble and solve the PETSC system.
  class PETSCBoundaryFaceSolveApply : public apply_rule<blackbox<PETSCVector>,Loci::NullOp<PETSCVector> > {
    private:
      const_param<vector<int> > fsiNumBoundaryFace ;
      const_blackbox<PETSCMultiVector> b ;
      const_blackbox<PETSCMatrix> A ;
      const_blackbox<PETSCKsp> ksp ;
      blackbox<PETSCMultiVector> phi;
    public:

      // Define input and output.
      PETSCBoundaryFaceSolveApply() {
        name_store("fsiNumBoundaryFace(X)",fsiNumBoundaryFace) ;
        name_store("petscBoundaryFacePhi(X)",phi) ;
        name_store("petscBoundaryFaceA(X)",A) ;
        name_store("petscBoundaryFaceB(X)",b) ;
        name_store("petscBoundaryFaceKSP(X)",ksp) ;
        input("petscBoundaryFaceA(X),petscBoundaryFaceB(X),petscBoundaryFaceKSP(X)") ;
        input("fsiNumBoundaryFace(X)") ;
        output("petscBoundaryFacePhi(X)") ;
		constraint("X") ;
        disable_threading() ;
      }

      // Assemble and solve.
	  virtual void compute(const sequence &seq) { }
  } ;

  register_rule<PETSCBoundaryFaceSolveApply> registerPETSCBoundaryFaceSolveApply ;

//-----------------------------Weight---------------------------------------------------------
//

  // Sets up the fsiWeight vector.
  class FSIBoundaryFaceAssembleWeightUnit : public unit_rule { // see interpolateFile.cc
    private:
	  param<vector<real> > fsiWeight ;
    public:

      // Define input and output.
      FSIBoundaryFaceAssembleWeightUnit() {
		name_store("fsiBoundaryFaceWeight(X)",fsiWeight) ;
        output("fsiBoundaryFaceWeight(X)") ;
        constraint("faces") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {}
  } ;
 
  register_rule<FSIBoundaryFaceAssembleWeightUnit> registerFSIBoundaryFaceAssembleWeightUnit ;

  // Assemble Q -> unit rule to output blackbox
  class FSIBoundaryFaceAssembleWeightApply : public apply_rule<param<vector<real> >,Loci::NullOp<vector<real> > > {
    private:
      const_param<vector<int> > fsiNumBoundaryFace ;
      const_blackbox<PETSCMultiVector> phi;
	  param<vector<real> > fsiWeight ;
    public:
      // Define input and output.
      FSIBoundaryFaceAssembleWeightApply() {
        name_store("fsiNumBoundaryFace(X)",fsiNumBoundaryFace) ;
        name_store("fsiBoundaryFaceWeight(X)", fsiWeight) ;
        name_store("petscBoundaryFacePhi(X)", phi) ;
        input("fsiNumBoundaryFace(X)") ;
        input("petscBoundaryFacePhi(X)") ;
        output("fsiBoundaryFaceWeight(X)") ;
				constraint("X");
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {

        // Get the number of local and global nodes.
		const int p = Loci::MPI_processes ;
		const int r = Loci::MPI_rank ;
        int localNumFace=(*fsiNumBoundaryFace)[r],globalNumFace=0 ;
        for(unsigned int i=0;i<(*fsiNumBoundaryFace).size();++i) globalNumFace+=(*fsiNumBoundaryFace)[i] ;
		
		// Compute the row offset for this process.
        int localStart=0 ; 
        for(int i=0;i<r;++i) localStart+=(*fsiNumBoundaryFace)[i]  ;  // 3 coordinates

		// include the polynomial coeff. computation
		if (r == p-1) { // last rank takes care of beta
			localNumFace += 4; // d+1 = 4
		}
					
		globalNumFace+=4; // for beta (d+1)
					
		vector<int> recvcounts(p,0) ;
		vector<int> displs(p,0) ;
		for(int i=0;i<p;++i) {
			for(int j=0;j<i;++j) displs[i]+=(3*(*fsiNumBoundaryFace)[j]);
			recvcounts[i] = 3*((i==p-1)?((*fsiNumBoundaryFace)[i]+4):(*fsiNumBoundaryFace)[i]) ;
		}
				
		// Allocate fsiQ
		(*fsiWeight).resize(3*globalNumFace) ; // number of nodes * 3 coordinates

		int counter=0;
        PetscScalar *phixCopy,*phiyCopy,*phizCopy ;
        (*phi).GetArray(&phixCopy,&phiyCopy,&phizCopy) ;
		// Fill in fsiQ
		for(int i=localStart; i<localStart+localNumFace; ++i,++counter) {
			(*fsiWeight)[3*i+0] = phixCopy[counter] ;
			(*fsiWeight)[3*i+1] = phiyCopy[counter] ;
			(*fsiWeight)[3*i+2] = phizCopy[counter] ;
		}
        (*phi).RestoreArray(&phixCopy,&phiyCopy,&phizCopy) ; // according to manual inexpensive

		// Allgatherv
		MPI_Allgatherv(&(*fsiWeight)[3*localStart], 3*localNumFace, MPI_DOUBLE, &(*fsiWeight)[0], &recvcounts[0], &displs[0], MPI_DOUBLE, MPI_COMM_WORLD);
      }
  } ;

  register_rule<FSIBoundaryFaceAssembleWeightApply> registerFSIBoundaryFaceAssembleWeightApply ;
//---------------------------------------------apply-------------------------------------------------
//
  // Rule to apply interpolation to the internal nodes
  class PETSCInternalFaceApply : public blackbox_rule {
    private:
	  const_param<real> r, a ;
	  const_param<int> fsiNr ;
	  const_blackbox<ublas::matrix<real,ublas::column_major> > CSDnodes ;
	  const_param<vector<real> > Qtop ;
	  const_param<vector<real> > Qbottom ;
	  const_param<vector<real> > fsiWeightTop ;
	  const_param<vector<real> > fsiWeightBottom ;
	  const_param<vector<int> > fsiNumBoundaryFace ;
	  store_blackbox<ublas::matrix<real,ublas::column_major> > CSDForce ;
    private:
	  int globalNumFace;
    public:

      // Define input and output.
      PETSCInternalFaceApply() {
		name_store("fsiBoundaryFaceQ(topCSDNodes)",Qtop) ;
		name_store("fsiBoundaryFaceQ(bottomCSDNodes)",Qbottom) ;
		name_store("fsiBoundaryFaceWeight(topCSDNodes)",fsiWeightTop) ;
		name_store("fsiBoundaryFaceWeight(bottomCSDNodes)",fsiWeightBottom) ;
		name_store("fsiR",r) ;
		name_store("fsiA",a) ;
		name_store("fsiNumber",fsiNr) ;
		name_store("facecenter",faceCenter) ;
		name_store("pos",pos) ;		
    name_store("fsiNumBoundaryFace",fsiNumBoundaryFace) ;
    name_store("CSDnodes",CSDNodes) ;
    name_store("CSDForce",CSDForce) ;
		input("fsiNumber, fsiR, fsiA");		
		input("facecenter") ;
		input("fsiNumBoundaryFace");
		input("fsiBoundaryFaceQ(topCSDNodes)") ;
		input("fsiBoundaryFaceQ(bottomCSDNodes)") ;
		input("fsiBoundaryFaceWeight(topCSDNodes)") ;
		input("fsiBoundaryFaceWeight(bottomCSDNodes)") ;
		input("CSDNodes") ;
    output("CSDForce") ;
    constraint("fsiCoupling") ; //internal nodes
    constraint("topCSDNodes,bottomCSDNodes") ; //internal nodes
    disable_threading() ;
      }
	  
	  // Loop over internal nodes and apply FSI interpolation
      void calculate(Entity node) {
        	

      }

      // Assemble and solve.
	  virtual void compute(const sequence &seq) {
			
        // Get the number of local and global nodes.
		const int p = Loci::MPI_processes ;
		const int rank = Loci::MPI_rank ;
    int localNumFace=(*fsiNumBoundaryFace)[rank],globalNumFace=0 ;
    for(unsigned int i=0;i<(*fsiNumBoundaryFace).size();++i) globalNumFace+=(*fsiNumBoundaryFace)[i] ;
//		cout << "p,r,localNumFace,globalNumFace: " << p << delim << rank << delim << localNumFace << delim << globalNumFace << endl ;

		double starttime, endtime;
		// Solve the linear system
		cout << "p,r: Start FSI CFD2CSD Applying" << Loci::MPI_processes << delim << Loci::MPI_rank << endl ;
		starttime = MPI_Wtime();			
		double distance = 0. ;
		ublas::matrix<real,ublas::row_major> value ; value.resize(CSDNodes.size1(), 2);
		for(int n=0;n<value.size1();++n) {
			int I,I1,I2;
			for(int i=0;i<globalNumFace;++i){
				I=3*i; I1=3*i+1; I2=3*i+2;
				distance = pow(CSDNodes(n,0) - (*QTop)[I], 2.) ;
				distance += pow(CSDNodes(n,1) -(*QTop)[I1], 2.) ;
				distance += pow(CSDNodes(n,2) - (*QTop)[I2], 2.) ;
				distance = sqrt(distance);
				value(i,0) = radialBasisFunction(distance, *r, *a, *fsiNr) ;
				distance = pow(CSDNodes(n,0) - (*QBottom)[I], 2.) ;
				distance += pow(CSDNodes(n,1) -(*QBottom)[I1], 2.) ;
				distance += pow(CSDNodes(n,2) - (*QBottom)[I2], 2.) ;
				distance = sqrt(distance);
				value(i,1) = radialBasisFunction(distance, *r, *a, *fsiNr) ;
			}
			for(int d=0;d<3;++d) {
				double polynomTop = (*fsiWeightTop)[3*(globalNumFace+0)+d];
				polynomTop+=(*fsiWeightTop)[3*(globalNumFace+1)+d]*CSDNodes(n,0) ;
				polynomTop+=(*fsiWeightTop)[3*(globalNumFace+2)+d]*CSDNodes(n,1) ;
				polynomTop+=(*fsiWeightTop)[3*(globalNumFace+3)+d]*CSDNodes(n,2) ;
				double polynomBottom = (*fsiWeightBottom)[3*(globalNumFace+0)+d];
				polynomBottom+=(*fsiWeightBottom)[3*(globalNumFace+1)+d]*CSDNodes(n,0) ;
				polynomBottom+=(*fsiWeightBottom)[3*(globalNumFace+2)+d]*CSDNodes(n,1) ;
				polynomBottom+=(*fsiWeightBottom)[3*(globalNumFace+3)+d]*CSDNodes(n,2) ;
				CSDForce(n,d) = polynomTop-polynomBottom ;
			}
				
			// FSI contribution
			for(int i=0;i<globalNumFace;++i){		
				for(int d=0;d<3;++d) {
					CSDForce(n,d) += (*fsiWeightTop)[3*i+d] * value[i,0] + (*fsiWeightBottom)[3*i+d] * value[i,1] ;
				}
			}
		}			
		endtime = MPI_Wtime();
		cout << "p,r: End FSI CFD2CSD" << p << delim << rank << delim << "Applied in " << endtime - starttime << " sec" << endl ;
	  }	
  } ;

  register_rule<PETSCInternalFaceApply> registerPETSCInternalFaceApply ;

}


