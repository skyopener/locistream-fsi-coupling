//-----------------------------------------------------------------------------
// Description: This file contains rules for the PETSC linear solver. We have
//   temporarily reinstated this file due to a possible bug in Ed's version
//   which comes with the finite-volume module.
//
// Authors: Siddharth Thakur and Jeff Wright
//          Streamline Numerics Inc.
//          Gainesville, Florida 32609
//-----------------------------------------------------------------------------

// Standard library includes.
#include<vector>
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
      int SetFromOptions() {
        int ierr=KSPSetFromOptions(ksp) ; CHKERRQ(ierr) ; return 0 ;
      }
      int SetOperators(const PETSCMatrix &m) const {
        int ierr=KSPSetOperators(ksp,m.Data(),m.Data(),SAME_NONZERO_PATTERN) ;
        CHKERRQ(ierr) ; return 0 ;
      }
      void SetTolerances(real relativeTolerance,real absoluteTolerance,int
      maxIterations) const {
        KSPSetTolerances(ksp,relativeTolerance,absoluteTolerance,PETSC_DEFAULT,
          maxIterations) ;
      }
      int Solve(const PETSCVector &rhs, PETSCVector &sol) const {
//      int numIteration ;
        int ierr=KSPSolve(ksp,rhs.Data(),sol.Data()) ; CHKERRQ(ierr) ;
//KSPGetIterationNumber(ksp,&numIteration) ;
//cout << "petsc iterations: " << numIteration << endl ;
        return ierr ;
      }
  } ;

//-----------------------------------------------------------------------------
// Rules for setting up the PETSC data.

  // Gets the number of cells assigned on all processes. This rule is
  // definitely not in the Loci spirit since we are collecting data from
  // other processes.
  class PETSCGetLocalCellData : public singleton_rule {
    private:
      param<vector<int> > petscNumCell ;
      param<entitySet> petscLocalCell ;
    public:

      // Define input and output.
      PETSCGetLocalCellData() {
        name_store("petscNumCellTemp",petscNumCell) ;
        name_store("petscLocalCellTemp",petscLocalCell) ;
        output("petscNumCellTemp,petscLocalCellTemp") ;
        constraint("geom_cells") ;
        disable_threading() ;
      }

      // Get the number of cells for each process.
      virtual void compute(const sequence &seq) {

        // Get the collection of entities assigned to this processor
        Loci::storeRepP myEntities=Loci::exec_current_fact_db->get_variable
          ("my_entities") ;
        entitySet localEntities=~EMPTY ;
        if(myEntities!=0) localEntities=(*myEntities).domain() ;

        // Get the local cells, not including the ghost cells.
        entitySet localCellWithGhost=entitySet(seq) ;
        *petscLocalCell=(localEntities & localCellWithGhost) ;

        // Distribute the number of cells to all processes.
        *petscNumCell=Loci::all_collect_sizes((*petscLocalCell).size()) ;
      }
  } ;

  register_rule<PETSCGetLocalCellData> registerPETSCGetLocalCellData ;

  // Creates the cell-to-row map.
  class PETSCCellToRow : public pointwise_rule {
    private:
      const_param<vector<int> > petscNumCell ;
      store<int> petscCellToRow ;
    public:

      // Define input and output.
      PETSCCellToRow() {
        name_store("petscNumCellTemp",petscNumCell) ;
        name_store("petscCellToRowTemp",petscCellToRow) ;
        input("petscNumCellTemp") ;
        output("petscCellToRowTemp") ;
        constraint("geom_cells") ;
        disable_threading() ;
      }

      // Set the Petsc row for each cell.
      virtual void compute(const sequence & seq) {

        // Compute the row offset for this process.
        int offset=0 ;
        for(int i=0;i<Loci::MPI_rank;++i){ offset+=(*petscNumCell)[i] ; }

        // Assign row number.
        sequence::const_iterator cellPtr=seq.begin() ;
        for(int i=0;i<(*petscNumCell)[Loci::MPI_rank];++cellPtr,++i){
          petscCellToRow[*cellPtr]=offset+i ;
        }
      }
  } ;

  register_rule<PETSCCellToRow> registerPETSCCellToRow ;

  // Rule that copies petscCellToRow for periodic boundaries
  class PETSCCellToRowPeriodic : public pointwise_rule {
    private:
      const_Map cl,cr,pmap ;
      store<int> petscCellToRow ;
    public:

      // Define input and output.
      PETSCCellToRowPeriodic() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pmap",pmap) ;
        name_store("petscCellToRowTemp",petscCellToRow) ;
        input("pmap->cl->petscCellToRowTemp") ;
        output("cr->petscCellToRowTemp") ;
      }

      // Copy the values to ghost cells.
      void calculate(Entity face) {
        petscCellToRow[cr[face]]=petscCellToRow[cl[pmap[face]]] ;
      }

      // Call calculate for all periodic faces. 
      void compute(const sequence & seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PETSCCellToRowPeriodic> registerPETSCCellToRowPeriodic ;

  // Determines the number of non-zero entries in the local portion of the
  // Petsc matrix for each cell. The local portion is defined as the square
  // row/column sub-block of the PETSC matrix whose rows map to cells local
  // to the process.
  class PETSCNumDiagonalNonZeroEntries : public pointwise_rule {
    private:
      const_param<entitySet> petscLocalCell ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      store<int> petscNumDiagonalNonZero ;
    public:

      // Define input and output.
      PETSCNumDiagonalNonZeroEntries() {
        name_store("petscLocalCellTemp",petscLocalCell) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("petscNumDiagonalNonZeroTemp",petscNumDiagonalNonZero) ;
        input("petscLocalCellTemp,upper->cr,lower->cl") ;
        output("petscNumDiagonalNonZeroTemp") ;
        constraint("geom_cells") ;
      }

      // Initialize count with the diagonal entry.
      void addDiagonal(Entity cell) { petscNumDiagonalNonZero[cell]=1 ; }

      // Add the lower neighbor entry to the count if the neighbor is local.
      // A neighbor is local if it appears in the set of local entities.
      void addLower(Entity cell) {
        int numNeighbor=lower.num_elems(cell) ;
        for(int i=0;i<numNeighbor;++i) {
          if((*petscLocalCell).inSet(cl[lower[cell][i]]))
            ++petscNumDiagonalNonZero[cell] ;
        }
      }

      // Add the upper neighbor entry to the count if the neighbor is local.
      // A neighbor is local if it appears in the set of local entities.
      void addUpper(Entity cell) {
        int numNeighbor=upper.num_elems(cell) ;
        for(int i=0;i<numNeighbor;++i){
          if((*petscLocalCell).inSet(cr[upper[cell][i]]))
            ++petscNumDiagonalNonZero[cell] ;
        }
      }

      // Compute the number of local non-zero entries.
      void compute(const sequence & seq) {
        do_loop(seq,this,&PETSCNumDiagonalNonZeroEntries::addDiagonal) ;
        do_loop(seq,this,&PETSCNumDiagonalNonZeroEntries::addLower) ;
        do_loop(seq,this,&PETSCNumDiagonalNonZeroEntries::addUpper) ;
      }
  } ;

  register_rule<PETSCNumDiagonalNonZeroEntries>
    registerPETSCNumDiagonalNonZeroEntries ;

  // Determines the number of non-zero entries in the non-local portion of the
  // Petsc matrix for each cell. For each row, this is just the total number
  // of row entries minus the previously computed number of local entries.
  class PETSCNumOffDiagonalNonZeroEntries : public pointwise_rule {
    private:
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<int> petscNumDiagonalNonZero ;
      store<int> petscNumOffDiagonalNonZero ;
    public:

      // Define input and output.
      PETSCNumOffDiagonalNonZeroEntries() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("petscNumDiagonalNonZeroTemp",petscNumDiagonalNonZero) ;
        name_store("petscNumOffDiagonalNonZeroTemp",
          petscNumOffDiagonalNonZero) ;
        input("upper->cr,lower->cl,petscNumDiagonalNonZeroTemp") ;
        output("petscNumOffDiagonalNonZeroTemp") ;
        constraint("geom_cells") ;
      }

      // Compute the number of off-diagonal entries for each cell.
      void calculate(Entity cell) {
        petscNumOffDiagonalNonZero[cell]=upper.num_elems(cell)+lower.
          num_elems(cell)+1-petscNumDiagonalNonZero[cell] ;
      }

      // Call calculate for all cells.
      void compute(const sequence & seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PETSCNumOffDiagonalNonZeroEntries>
    registerPETSCNumOffDiagonalNonZeroEntries ;

  // Sets up the PETSC right-hand-side vector.
  class PETSCSetupRHS : public singleton_rule {
    private:
      const_param<vector<int> > petscNumCell ;
      blackbox<PETSCVector> b ;
    public:

      // Define input and output.
      PETSCSetupRHS() {
        name_store("petscNumCellTemp",petscNumCell) ;
        name_store("petscBTemp(X)",b) ;
        input("petscNumCellTemp") ;
        output("petscBTemp(X)") ;
        constraint("X_PETSCLinearSolver,geom_cells") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {

        // Get the number of local and global cells.
        int localNumCell=(*petscNumCell)[Loci::MPI_rank],globalNumCell=0 ;
        for(unsigned int i=0;i<(*petscNumCell).size();++i) globalNumCell+=
          (*petscNumCell)[i] ;

        // Allocate the unknown and rhs vectors.
        (*b).Create(localNumCell,globalNumCell) ;
      }
  } ;

  register_rule<PETSCSetupRHS> registerPETSCSetupRHS ;

  // Sets up the PETSC matrix. Note that this is a unit_rule, since this is the
  // only way we can have stores as input to a rule outputting blackboxes.
  class PETSCSetupMatrixUnit : public unit_rule {
    private:
      const_store<int> petscNumDiagonalNonZero ;
      const_store<int> petscNumOffDiagonalNonZero ;
      const_param<vector<int> > petscNumCell ;
      blackbox<PETSCMatrix> A ;
    public:

      // Define input and output.
      PETSCSetupMatrixUnit() {
        name_store("petscNumDiagonalNonZeroTemp",petscNumDiagonalNonZero) ;
        name_store("petscNumOffDiagonalNonZeroTemp",
          petscNumOffDiagonalNonZero) ;
        name_store("petscNumCellTemp",petscNumCell) ;
        name_store("petscATemp(X)",A) ;
        input("petscNumDiagonalNonZeroTemp,petscNumOffDiagonalNonZeroTemp") ;
        input("petscNumCellTemp") ;
        output("petscATemp(X)") ;
        constraint("X_PETSCLinearSolver,geom_cells") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {

        // Get the number of local and global cells.
        int localNumCell=(*petscNumCell)[Loci::MPI_rank],globalNumCell=0 ;
        for(unsigned int i=0;i<(*petscNumCell).size();++i) globalNumCell+=
          (*petscNumCell)[i] ;

        // Make temporary copy of matrix allocation data.
        int count=0,*numDiagonalNonZero=new int[localNumCell],
          *numOffDiagonalNonZero=new int[localNumCell] ;
        for(sequence::const_iterator cellPtr=seq.begin();cellPtr!=seq.end();
        ++cellPtr,++count){
          numDiagonalNonZero[count]=petscNumDiagonalNonZero[*cellPtr] ;
          numOffDiagonalNonZero[count]=petscNumOffDiagonalNonZero[*cellPtr] ;
        }

        // Allocate the matrix.
        (*A).Create(localNumCell,globalNumCell,numDiagonalNonZero,
          numOffDiagonalNonZero) ;

        // Deallocate temporary copy of matrix allocation data.
        delete [] numDiagonalNonZero ; delete [] numOffDiagonalNonZero ;
      }
  } ;

  register_rule<PETSCSetupMatrixUnit> registerPETSCSetupMatrixUnit ;

  // Empty apply rule required by Loci. The data type and operator do not
  // matter since nothing is done by this rule. Keep the same inputs and
  // outputs as the unit rule, even though we don't have to.
  class PETSCSetupMatrixApply : public apply_rule<blackbox<PETSCMatrix>,
  Loci::NullOp<PETSCMatrix> > {
    private:
      const_store<int> petscNumDiagonalNonZero ;
      const_store<int> petscNumOffDiagonalNonZero ;
      const_param<vector<int> > petscNumCell ;
      blackbox<PETSCMatrix> A ;
    public:

      // Define input and output.
      PETSCSetupMatrixApply() {
        name_store("petscNumDiagonalNonZeroTemp",petscNumDiagonalNonZero) ;
        name_store("petscNumOffDiagonalNonZeroTemp",
          petscNumOffDiagonalNonZero) ;
        name_store("petscNumCellTemp",petscNumCell) ;
        name_store("petscATemp(X)",A) ;
        input("petscNumDiagonalNonZeroTemp,petscNumOffDiagonalNonZeroTemp") ;
        input("petscNumCellTemp") ;
        output("petscATemp(X)") ;
        constraint("X_PETSCLinearSolver,geom_cells") ;
        disable_threading() ;
      }

      // Do nothing.
      void compute(const sequence & seq) {}
  } ;

  register_rule<PETSCSetupMatrixApply> registerPETSCSetupMatrixApply ;

  // Sets up the PETSC linear solver.
  class PETSCSetupSolver : public singleton_rule {
    private:
      const_param<int> maxLinearSolverIterations ;
      const_param<real> petscConvergenceTolerance ;
      const_blackbox<PETSCMatrix> A ;
      blackbox<PETSCKsp> ksp ;
    public:

      // Define input and output.
      PETSCSetupSolver() {
        name_store("X_maxLinearSolverIterations",maxLinearSolverIterations) ;
        name_store("petscConvergenceTolerance",petscConvergenceTolerance) ;
        name_store("petscATemp(X)",A) ;
        name_store("petscKSPTemp(X)",ksp) ;
        input("petscATemp(X),X_maxLinearSolverIterations") ;
        input("petscConvergenceTolerance") ;
        output("petscKSPTemp(X)") ;
        constraint("X_PETSCLinearSolver,geom_cells") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {
        (*ksp).Create() ;
        (*ksp).SetTolerances(*petscConvergenceTolerance,1.0e-30,
          *maxLinearSolverIterations) ;
        PC pc ; (*ksp).GetPC(&pc) ; PCSetType(pc,PCJACOBI) ;
        (*ksp).SetFromOptions() ;
      }
  } ;

  register_rule<PETSCSetupSolver> registerPETSCSetupSolver ;

  // Assemble and solve the PETSC system.
  class PETSCSolveUnit : public unit_rule {
    private:
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> xU,xL,xB,xD ;
      const_store<int> petscCellToRow ;
      const_blackbox<PETSCVector> b ;
      const_blackbox<PETSCMatrix> A ;
      const_blackbox<PETSCKsp> ksp ;
      blackbox<PETSCVector> phi ;
    private:
      int *columnIndex ;
      PetscScalar *columnValue ;
    public:

      // Define input and output.
      PETSCSolveUnit() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("X_U",xU) ;
        name_store("X_L",xL) ;
        name_store("X_B",xB) ;
        name_store("X_D",xD) ;
        name_store("petscCellToRowTemp",petscCellToRow) ;
        name_store("petscBTemp(X)",b) ;
        name_store("petscATemp(X)",A) ;
        name_store("petscKSPTemp(X)",ksp) ;
        name_store("petscPhiTemp(X)",phi) ;
        input("upper->cr->petscCellToRowTemp") ;
        input("lower->cl->petscCellToRowTemp") ;
        input("upper->X_U,lower->X_L,X_B,X_D") ;
        input("petscATemp(X),petscBTemp(X),petscKSPTemp(X)") ;
        output("petscPhiTemp(X)") ;
        constraint("X_PETSCLinearSolver,geom_cells") ;
        disable_threading() ;
      }

      // Assemble and solve.
      virtual void compute(const sequence &seq) {

        // Create the solution vector. Duplication only occurs the first time
        // this function is called. Later calls just return.
        (*phi).DuplicateFrom(*b) ;

        // Allocate memory for column index and values arrays. There should
        // never be more than 100 entries.
        columnIndex=new int[100] ; columnValue=new PetscScalar[100] ;

        // Assemble the matrix.
        do_loop(seq,this,&PETSCSolveUnit::AssembleMatrix) ;
        (*A).AssemblyBegin() ; (*A).AssemblyEnd() ;
        
        // Dellocate memory.
        delete [] columnIndex ; delete [] columnValue ;

        // Assemble the rhs vector.
        do_loop(seq,this,&PETSCSolveUnit::AssembleRHS) ;
        (*b).AssemblyBegin() ; (*b).AssemblyEnd() ;

        // Solve the linear system
        (*ksp).SetOperators(*A) ; (*ksp).Solve(*b,*phi) ;
      }

      // Sets matrix values for a cell's row.
      void AssembleMatrix(Entity cell) {

        // Get the row number for this cell.
        int rowIndex=petscCellToRow[cell] ;

        // Set the diagonal entry for the row.
        columnIndex[0]=rowIndex ; columnValue[0]=xD[cell] ;

        // Set the column entries for "upper" neighbors. Note the new
        // if check for periodic boundaries which avoids possible
        // self-reference.
        const int numUpperEntry=upper.num_elems(cell) ; int count=1 ;
        for(int i=0;i<numUpperEntry;++i,++count){
          if(petscCellToRow[cell]!=petscCellToRow[cr[upper[cell][i]]]){
            columnIndex[count]=petscCellToRow[cr[upper[cell][i]]] ;
            columnValue[count]=xU[upper[cell][i]] ;
          }else{
            cout << "Self-reference in AssembleMatrix" << endl ;
          }
        }

        // Set the column entries for "lower" neighbors.
        const int numLowerEntry=lower.num_elems(cell) ;
        for(int i=0;i<numLowerEntry;++i,++count){
          columnIndex[count]=petscCellToRow[cl[lower[cell][i]]] ;
          columnValue[count]=xL[lower[cell][i]] ;
        }

        // Insert row entries into the Petsc matrix.
        (*A).SetRowValues(rowIndex,numUpperEntry+numLowerEntry+1,columnIndex,
          columnValue) ;
      }

      // Sets matrix values for a cell's row.
      void AssembleRHS(Entity cell) {
        const PetscScalar value=xB[cell] ;
        (*b).SetValue(&(petscCellToRow[cell]),&value) ;
      }
  } ;

  register_rule<PETSCSolveUnit> registerPETSCSolveUnit ;

  // Empty apply rule required by Loci. The data type and operator do not
  // matter since nothing is done by this rule. Keep the same inputs and
  // outputs as the unit rule, even though we don't have to.
  class PETSCSolveApply : public apply_rule<blackbox<PETSCVector>,
  Loci::NullOp<PETSCVector> > {
    private:
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> xU,xL,xB,xD ;
      const_store<int> petscCellToRow ;
      const_blackbox<PETSCVector> b ;
      const_blackbox<PETSCMatrix> A ;
      const_blackbox<PETSCKsp> ksp ;
      blackbox<PETSCVector> phi ;
    public:

      // Define input and output.
      PETSCSolveApply() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("X_U",xU) ;
        name_store("X_L",xL) ;
        name_store("X_B",xB) ;
        name_store("X_D",xD) ;
        name_store("petscCellToRowTemp",petscCellToRow) ;
        name_store("petscBTemp(X)",b) ;
        name_store("petscATemp(X)",A) ;
        name_store("petscKSPTemp(X)",ksp) ;
        name_store("petscPhiTemp(X)",phi) ;
        input("upper->cr->petscCellToRowTemp") ;
        input("lower->cl->petscCellToRowTemp") ;
        input("upper->X_U,lower->X_L,X_B,X_D") ;
        input("petscATemp(X),petscBTemp(X),petscKSPTemp(X)") ;
        output("petscPhiTemp(X)") ;
        constraint("X_PETSCLinearSolver,geom_cells") ;
        disable_threading() ;
      }

      // Do nothing.
      void compute(const sequence & seq) {}
  } ;

  register_rule<PETSCSolveApply> registerPETSCSolveApply ;

  // Non-parametric rule for extracting pPrime from PETSC.
  class PETSCCopyPPrime : public pointwise_rule {
    private:
      const_blackbox<PETSCVector> petscPPrime ;
      store<real> pPrime ;
    public:

      // Define input and output.
      PETSCCopyPPrime() {
        name_store("petscPhiTemp(pPrime)",petscPPrime) ;
        name_store("pPrime",pPrime) ;
        input("petscPhiTemp(pPrime)") ;
        output("pPrime") ;
        constraint("pPrime_PETSCLinearSolver,geom_cells") ;
        disable_threading() ;
      }

      // Copy the solution back from PETSC.
      void compute(const sequence & seq) {
        int minRowNum,maxRowNum ;
        (*petscPPrime).GetOwnershipRange(&minRowNum,&maxRowNum) ;
        PetscScalar *pPrimeCopy ; (*petscPPrime).GetArray(&pPrimeCopy) ;
        sequence::const_iterator cellPtr=seq.begin() ;
        for(int row=minRowNum;row<maxRowNum;++row,++cellPtr) pPrime[*cellPtr]=
          pPrimeCopy[row-minRowNum] ;
        (*petscPPrime).RestoreArray(&pPrimeCopy) ;
      }
  } ;

  register_rule<PETSCCopyPPrime> registerPETSCCopyPPrime ;

  // Non-parametric rule for extracting pPrimeStar from PETSC.
  class PETSCCopyPPrimeStar : public pointwise_rule {
    private:
      const_blackbox<PETSCVector> petscPPrimeStar ;
      store<real> pPrimeStar ;
    public:

      // Define input and output.
      PETSCCopyPPrimeStar() {
        name_store("petscPhiTemp(pPrimeStar)",petscPPrimeStar) ;
        name_store("pPrimeStar",pPrimeStar) ;
        input("petscPhiTemp(pPrimeStar)") ;
        output("pPrimeStar") ;
        constraint("pPrimeStar_PETSCLinearSolver,geom_cells") ;
        disable_threading() ;
      }

      // Copy the solution back from PETSC.
      void compute(const sequence & seq) {
        int minRowNum,maxRowNum ;
        (*petscPPrimeStar).GetOwnershipRange(&minRowNum,&maxRowNum) ;
        PetscScalar *pPrimeStarCopy ;
        (*petscPPrimeStar).GetArray(&pPrimeStarCopy) ;
        sequence::const_iterator cellPtr=seq.begin() ;
        for(int row=minRowNum;row<maxRowNum;++row,++cellPtr)
          pPrimeStar[*cellPtr]=pPrimeStarCopy[row-minRowNum] ;
        (*petscPPrimeStar).RestoreArray(&pPrimeStarCopy) ;
      }
  } ;

  register_rule<PETSCCopyPPrimeStar> registerPETSCCopyPPrimeStar ;
}







