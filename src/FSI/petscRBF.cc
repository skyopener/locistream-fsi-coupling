//-----------------------------------------------------------------------------
// Description: This file contains rules for implementing a RBF linear solver
// using PETSC.
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
using Loci::const_MapVec ;
using Loci::const_storeVec ;

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

//-----------------------------------------------------------------------------
// Rules for setting up the PETSC data.

  // Gets the number of boundary nodes assigned on all processes. This rule is
  // definitely not in the Loci spirit since we are collecting data from
  // other processes.
  class PETSCGetLocalBoundaryNodeDataX : public singleton_rule {
    private:
      param<vector<int> > petscNumBoundaryNodeX ;
      param<entitySet> petscLocalBoundaryNodeX ;
    public:

      // Define input and output.
      PETSCGetLocalBoundaryNodeDataX() {
        name_store("petscNumNode",petscNumNode) ;
        name_store("petscLocalNode",petscLocalNode) ;
        output("petscNumNode,petscLocalNode") ;
        constraint("pos") ;
        disable_threading() ;
      }

      // Get the number of nodes for each process.
      virtual void compute(const sequence &seq) {

        // Get the collection of entities assigned to this processor
        Loci::storeRepP myEntities=Loci::exec_current_fact_db->get_variable
          ("my_entities") ;
        entitySet localEntities=~EMPTY ;
        if(myEntities!=0) localEntities=(*myEntities).domain() ;

        // Get the local nodes.
        *petscLocalNode=entitySet(seq) ;

        // Distribute the number of nodes to all processes.
        *petscNumNode=Loci::all_collect_sizes((*petscLocalNode).size()) ;
      }
  } ;

  register_rule<PETSCGetLocalNodeData> registerPETSCGetLocalNodeData ;

  // Creates the node-to-row map.
  class PETSCNodeToRow : public pointwise_rule {
    private:
      const_param<vector<int> > petscNumNode ;
      store<int> petscNodeToRow ;
    public:

      // Define input and output.
      PETSCNodeToRow() {
        name_store("petscNumNode",petscNumNode) ;
        name_store("petscNodeToRow",petscNodeToRow) ;
        input("petscNumNode") ;
        output("petscNodeToRow") ;
        constraint("pos") ;
        disable_threading() ;
      }

      // Set the Petsc row for each node.
      virtual void compute(const sequence & seq) {

        // Compute the row offset for this process.
        int offset=0 ;
        for(int i=0;i<Loci::MPI_rank;++i){ offset+=(*petscNumNode)[i] ; }

        // Assign row number.
        sequence::const_iterator nodePtr=seq.begin() ;
        for(int i=0;i<(*petscNumNode)[Loci::MPI_rank];++nodePtr,++i){
          petscNodeToRow[*nodePtr]=offset+i ;
        }
      }
  } ;

  register_rule<PETSCNodeToRow> registerPETSCNodeToRow ;

  // Determines the number of non-zero entries in the local portion of the
  // Petsc matrix for each node. The local portion is defined as the square
  // row/column sub-block of the PETSC matrix whose rows map to nodes local
  // to the process.
  class PETSCNodeNumDiagonalNonZeroEntries : public pointwise_rule {
    private:
      const_param<entitySet> petscLocalNode ;
      const_multiMap node2edge ;
      const_MapVec<2> edge2node ;
      store<int> petscNodeNumDiagonalNonZero ;
    public:

      // Define input and output.
      PETSCNodeNumDiagonalNonZeroEntries() {
        name_store("petscLocalNode",petscLocalNode) ;
        name_store("node2edge",node2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("petscNodeNumDiagonalNonZero",petscNodeNumDiagonalNonZero) ;
        input("petscLocalNode,node2edge->edge2node") ;
        output("petscNodeNumDiagonalNonZero") ;
        constraint("pos") ;
      }

      // Count local nodes for this node.
      void calculate(Entity node) {
                                                                                
        // Initialize count with the diagonal entry.
        petscNodeNumDiagonalNonZero[node]=1 ;

        // Add neighbor node for each edge if it is local.
        int numEdge=node2edge.num_elems(node) ;
        for(int i=0;i<numEdge;++i) {
          if((*petscLocalNode).inSet((node==edge2node[node2edge[node][i]][0])?
          edge2node[node2edge[node][i]][1]:edge2node[node2edge[node][i]][0]))
            ++petscNodeNumDiagonalNonZero[node] ;
        }
      }

      // Loop over nodes.
      void compute(const sequence & seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PETSCNodeNumDiagonalNonZeroEntries>
    registerPETSCNodeNumDiagonalNonZeroEntries ;

  // Determines the number of non-zero entries in the non-local portion of the
  // Petsc matrix for each node.
  class PETSCNodeNumOffDiagonalNonZeroEntries : public pointwise_rule {
    private:
      const_param<entitySet> petscLocalNode ;
      const_multiMap node2edge ;
      const_MapVec<2> edge2node ;
      store<int> petscNodeNumOffDiagonalNonZero ;
    public:

      // Define input and output.
      PETSCNodeNumOffDiagonalNonZeroEntries() {
        name_store("petscLocalNode",petscLocalNode) ;
        name_store("node2edge",node2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("petscNodeNumOffDiagonalNonZero",
          petscNodeNumOffDiagonalNonZero) ;
        input("petscLocalNode,node2edge->edge2node") ;
        output("petscNodeNumOffDiagonalNonZero") ;
        constraint("pos") ;
      }

      // Count non-local nodes for this node.
      void calculate(Entity node) {
                                                                                
        // Initialize count.
        petscNodeNumOffDiagonalNonZero[node]=0 ;

        // Add neighbor node for each edge if it is not local.
        int numEdge=node2edge.num_elems(node) ;
        for(int i=0;i<numEdge;++i) {
          if(!(*petscLocalNode).inSet((node==edge2node[node2edge[node][i]][0])?
          edge2node[node2edge[node][i]][1]:edge2node[node2edge[node][i]][0]))
            ++petscNodeNumOffDiagonalNonZero[node] ;
        }
      }

      // Loop over nodes.
      void compute(const sequence & seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PETSCNodeNumOffDiagonalNonZeroEntries>
    registerPETSCNodeNumOffDiagonalNonZeroEntries ;

  // Sets up the PETSC right-hand-side vector.
  class PETSCNodeSetupRHS : public singleton_rule {
    private:
      const_param<vector<int> > petscNumNode ;
      blackbox<PETSCVector> b ;
    public:

      // Define input and output.
      PETSCNodeSetupRHS() {
        name_store("petscNumNode",petscNumNode) ;
        name_store("petscNodeB(X)",b) ;
        input("petscNumNode") ;
        output("petscNodeB(X)") ;
        constraint("X_PETSCLinearSolver,pos") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {

        // Get the number of local and global nodes.
        int localNumNode=(*petscNumNode)[Loci::MPI_rank],globalNumNode=0 ;
        for(unsigned int i=0;i<(*petscNumNode).size();++i) globalNumNode+=
          (*petscNumNode)[i] ;

        // Allocate the unknown and rhs vectors.
        (*b).Create(localNumNode,globalNumNode) ;
      }
  } ;

  register_rule<PETSCNodeSetupRHS> registerPETSCNodeSetupRHS ;

  // Sets up the PETSC right-hand-side vector for a vector equation.
  class PETSCNodeSetupVectorRHS : public singleton_rule {
    private:
      const_param<vector<int> > petscNumNode ;
      blackbox<PETSCMultiVector> b ;
    public:

      // Define input and output.
      PETSCNodeSetupVectorRHS() {
        name_store("petscNumNode",petscNumNode) ;
        name_store("petscNodeVectorB(X)",b) ;
        input("petscNumNode") ;
        output("petscNodeVectorB(X)") ;
        constraint("X_PETSCLinearSolver,pos") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {

        // Get the number of local and global nodes.
        int localNumNode=(*petscNumNode)[Loci::MPI_rank],globalNumNode=0 ;
        for(unsigned int i=0;i<(*petscNumNode).size();++i) globalNumNode+=
          (*petscNumNode)[i] ;

        // Allocate the unknown and rhs vectors.
        (*b).Create(localNumNode,globalNumNode) ;
      }
  } ;

  register_rule<PETSCNodeSetupVectorRHS> registerPETSCNodeSetupVectorRHS ;

  // Sets up the PETSC matrix. Note that this is a unit_rule, since this is the
  // only way we can have stores as input to a rule outputting blackboxes.
  class PETSCNodeSetupMatrixUnit : public unit_rule {
    private:
      const_store<int> petscNodeNumDiagonalNonZero ;
      const_store<int> petscNodeNumOffDiagonalNonZero ;
      const_param<vector<int> > petscNumNode ;
      blackbox<PETSCMatrix> A ;
    public:

      // Define input and output.
      PETSCNodeSetupMatrixUnit() {
        name_store("petscNodeNumDiagonalNonZero",petscNodeNumDiagonalNonZero) ;
        name_store("petscNodeNumOffDiagonalNonZero",
          petscNodeNumOffDiagonalNonZero) ;
        name_store("petscNumNode",petscNumNode) ;
        name_store("petscNodeA(X)",A) ;
        input("petscNodeNumDiagonalNonZero,petscNodeNumOffDiagonalNonZero") ;
        input("petscNumNode") ;
        output("petscNodeA(X)") ;
        constraint("X_PETSCLinearSolver,pos") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {

        // Get the number of local and global nodes.
        int localNumNode=(*petscNumNode)[Loci::MPI_rank],globalNumNode=0 ;
        for(unsigned int i=0;i<(*petscNumNode).size();++i) globalNumNode+=
          (*petscNumNode)[i] ;

        // Make temporary copy of matrix allocation data.
        int count=0,*numDiagonalNonZero=new int[localNumNode],
          *numOffDiagonalNonZero=new int[localNumNode] ;
        for(sequence::const_iterator nodePtr=seq.begin();nodePtr!=seq.end();
        ++nodePtr,++count){
          numDiagonalNonZero[count]=petscNodeNumDiagonalNonZero[*nodePtr] ;
          numOffDiagonalNonZero[count]=petscNodeNumOffDiagonalNonZero
            [*nodePtr] ;
        }

        // Allocate the matrix.
        (*A).Create(localNumNode,globalNumNode,numDiagonalNonZero,
          numOffDiagonalNonZero) ;

        // Deallocate temporary copy of matrix allocation data.
        delete [] numDiagonalNonZero ; delete [] numOffDiagonalNonZero ;
      }
  } ;

  register_rule<PETSCNodeSetupMatrixUnit> registerPETSCNodeSetupMatrixUnit ;

  // Empty apply rule required by Loci. The data type and operator do not
  // matter since nothing is done by this rule. Keep the same inputs and
  // outputs as the unit rule, even though we don't have to.
  class PETSCNodeSetupMatrixApply : public apply_rule<blackbox<PETSCMatrix>,
  Loci::NullOp<PETSCMatrix> > {
    private:
      const_store<int> petscNodeNumDiagonalNonZero ;
      const_store<int> petscNodeNumOffDiagonalNonZero ;
      const_param<vector<int> > petscNumNode ;
      blackbox<PETSCMatrix> A ;
    public:

      // Define input and output.
      PETSCNodeSetupMatrixApply() {
        name_store("petscNodeNumDiagonalNonZero",petscNodeNumDiagonalNonZero) ;
        name_store("petscNodeNumOffDiagonalNonZero",
          petscNodeNumOffDiagonalNonZero) ;
        name_store("petscNumNode",petscNumNode) ;
        name_store("petscNodeA(X)",A) ;
        input("petscNodeNumDiagonalNonZero,petscNodeNumOffDiagonalNonZero") ;
        input("petscNumNode") ;
        output("petscNodeA(X)") ;
        constraint("X_PETSCLinearSolver,pos") ;
        disable_threading() ;
      }

      // Do nothing.
      void compute(const sequence & seq) {}
  } ;

  register_rule<PETSCNodeSetupMatrixApply> registerPETSCNodeSetupMatrixApply ;

  // Sets up the PETSC linear solver.
  class PETSCNodeSetupSolver : public singleton_rule {
    private:
      const_param<int> maxLinearSolverIterations ;
      const_param<real> linearSolverTolerance ;
      const_blackbox<PETSCMatrix> A ;
      blackbox<PETSCKsp> ksp ;
    public:

      // Define input and output.
      PETSCNodeSetupSolver() {
        name_store("X_maxLinearSolverIterations",maxLinearSolverIterations) ;
        name_store("X_linearSolverTolerance",linearSolverTolerance) ;
        name_store("petscNodeA(X)",A) ;
        name_store("petscNodeKSP(X)",ksp) ;
        input("petscNodeA(X),X_maxLinearSolverIterations") ;
        input("X_linearSolverTolerance") ;
        output("petscNodeKSP(X)") ;
        constraint("X_PETSCLinearSolver,pos") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {
        (*ksp).Create() ;
        (*ksp).SetTolerances(*linearSolverTolerance,1.0e-30,
          *maxLinearSolverIterations) ;
      // Valve case.
      //PC pc ; (*ksp).GetPC(&pc) ; PCSetType(pc,PCILU) ; PCILUSetLevels(pc,3) ;
      // Injector case.
      PC pc ; (*ksp).GetPC(&pc) ; PCSetType(pc,PCJACOBI) ;
      // PCILUSetLevels(pc,0) ;
        (*ksp).SetFromOptions() ;
      }
  } ;

  register_rule<PETSCNodeSetupSolver> registerPETSCNodeSetupSolver ;

  // Assemble and solve the PETSC system.
  class PETSCNodeSolveUnit : public unit_rule {
    private:
      const_multiMap node2edge ;
      const_MapVec<2> edge2node ;
      const_storeVec<real> xE ;
      const_store<real> xB,xD ;
      const_store<int> petscNodeToRow ;
      const_blackbox<PETSCVector> b ;
      const_blackbox<PETSCMatrix> A ;
      const_blackbox<PETSCKsp> ksp ;
      blackbox<PETSCVector> phi ;
    private:
      int *columnIndex ;
      PetscScalar *columnValue ;
    public:

      // Define input and output.
      PETSCNodeSolveUnit() {
        name_store("node2edge",node2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("X_E",xE) ;
        name_store("X_B",xB) ;
        name_store("X_D",xD) ;
        name_store("petscNodeToRow",petscNodeToRow) ;
        name_store("petscNodeB(X)",b) ;
        name_store("petscNodeA(X)",A) ;
        name_store("petscNodeKSP(X)",ksp) ;
        name_store("petscNodePhi(X)",phi) ;
        input("node2edge->edge2node,node2edge->X_E,X_B,X_D") ;
        input("node2edge->edge2node->petscNodeToRow") ;
        input("petscNodeA(X),petscNodeB(X),petscNodeKSP(X)") ;
        output("petscNodePhi(X)") ;
        constraint("X_PETSCLinearSolver,pos") ;
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
        do_loop(seq,this,&PETSCNodeSolveUnit::AssembleMatrix) ;
        (*A).AssemblyBegin() ; (*A).AssemblyEnd() ;
        
        // Dellocate memory.
        delete [] columnIndex ; delete [] columnValue ;

        // Assemble the rhs vector.
        do_loop(seq,this,&PETSCNodeSolveUnit::AssembleRHS) ;
        (*b).AssemblyBegin() ; (*b).AssemblyEnd() ;

        // Solve the linear system
        (*ksp).SetOperators(*A) ; (*ksp).Solve(*b,*phi) ;
      }

      // Sets matrix values for a node's row.
      void AssembleMatrix(Entity node) {

        // Get the row number for this node.
        int rowIndex=petscNodeToRow[node] ;

        // Set the diagonal entry for the row.
        columnIndex[0]=rowIndex ; columnValue[0]=xD[node] ;

        // Set the column entries for edge neighbors.
        const int numEdge=node2edge.num_elems(node) ; int count=1 ;
        for(int i=0;i<numEdge;++i){
          bool iAmNodeZero=(node==edge2node[node2edge[node][i]][0]) ;
          if(iAmNodeZero){
            columnIndex[count]=petscNodeToRow[edge2node[node2edge[node][i]]
              [1]] ;
            columnValue[count]=xE[node2edge[node][i]][0] ; ++count ;
          }else{
            columnIndex[count]=petscNodeToRow[edge2node[node2edge[node][i]]
              [0]] ;
            columnValue[count]=xE[node2edge[node][i]][1] ; ++count ;
          }
        }

        // Insert row entries into the Petsc matrix.
        (*A).SetRowValues(rowIndex,numEdge+1,columnIndex,columnValue) ;
      }

      // Sets matrix values for a node's row.
      void AssembleRHS(Entity node) {
        const PetscScalar value=xB[node] ;
        (*b).SetValue(&(petscNodeToRow[node]),&value) ;
      }
  } ;

  register_rule<PETSCNodeSolveUnit> registerPETSCNodeSolveUnit ;

  // Empty apply rule required by Loci. The data type and operator do not
  // matter since nothing is done by this rule. Keep the same inputs and
  // outputs as the unit rule, even though we don't have to.
  class PETSCNodeSolveApply : public apply_rule<blackbox<PETSCVector>,
  Loci::NullOp<PETSCVector> > {
    private:
      const_multiMap node2edge ;
      const_MapVec<2> edge2node ;
      const_storeVec<real> xE ;
      const_store<real> xB,xD ;
      const_store<int> petscNodeToRow ;
      const_blackbox<PETSCVector> b ;
      const_blackbox<PETSCMatrix> A ;
      const_blackbox<PETSCKsp> ksp ;
      blackbox<PETSCVector> phi ;
    public:

      // Define input and output.
      PETSCNodeSolveApply() {
        name_store("node2edge",node2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("X_E",xE) ;
        name_store("X_B",xB) ;
        name_store("X_D",xD) ;
        name_store("petscNodeToRow",petscNodeToRow) ;
        name_store("petscNodeB(X)",b) ;
        name_store("petscNodeA(X)",A) ;
        name_store("petscNodeKSP(X)",ksp) ;
        name_store("petscNodePhi(X)",phi) ;
        input("node2edge->edge2node,node2edge->X_E,X_B,X_D") ;
        input("node2edge->edge2node->petscNodeToRow") ;
        input("petscNodeA(X),petscNodeB(X),petscNodeKSP(X)") ;
        output("petscNodePhi(X)") ;
        constraint("X_PETSCLinearSolver,pos") ;
        disable_threading() ;
      }

      // Do nothing.
      void compute(const sequence & seq) {}
  } ;

  register_rule<PETSCNodeSolveApply> registerPETSCNodeSolveApply ;

  // Assemble and solve the PETSC system for a vector rhs.
  class PETSCNodeSolveVectorUnit : public unit_rule {
    private:
      const_multiMap node2edge ;
      const_MapVec<2> edge2node ;
      const_storeVec<real> xE ;
      const_store<real> xD ;
      const_store<vect3d> xB ;
      const_store<int> petscNodeToRow ;
      const_blackbox<PETSCMultiVector> b ;
      const_blackbox<PETSCMatrix> A ;
      const_blackbox<PETSCKsp> ksp ;
      blackbox<PETSCMultiVector> phi ;
    private:
      int *columnIndex ;
      PetscScalar *columnValue ;
    public:

      // Define input and output.
      PETSCNodeSolveVectorUnit() {
        name_store("node2edge",node2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("X_E",xE) ;
        name_store("X_B",xB) ;
        name_store("X_D",xD) ;
        name_store("petscNodeToRow",petscNodeToRow) ;
        name_store("petscNodeVectorB(X)",b) ;
        name_store("petscNodeA(X)",A) ;
        name_store("petscNodeKSP(X)",ksp) ;
        name_store("petscNodeVectorPhi(X)",phi) ;
        input("node2edge->edge2node,node2edge->X_E,X_B,X_D") ;
        input("node2edge->edge2node->petscNodeToRow") ;
        input("petscNodeA(X),petscNodeVectorB(X),petscNodeKSP(X)") ;
        output("petscNodeVectorPhi(X)") ;
        constraint("X_PETSCLinearSolver,pos") ;
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
        do_loop(seq,this,&PETSCNodeSolveVectorUnit::AssembleMatrix) ;
        (*A).AssemblyBegin() ; (*A).AssemblyEnd() ;
        
        // Dellocate memory.
        delete [] columnIndex ; delete [] columnValue ;

        // Assemble the rhs vector.
        do_loop(seq,this,&PETSCNodeSolveVectorUnit::AssembleRHS) ;
        (*b).AssemblyBegin() ; (*b).AssemblyEnd() ;

        // Solve the linear system
        (*ksp).SetOperators(*A) ; (*ksp).SolveVector(*b,*phi) ;
      }

      // Sets matrix values for a node's row.
      void AssembleMatrix(Entity node) {

        // Get the row number for this node.
        int rowIndex=petscNodeToRow[node] ;

        // Set the diagonal entry for the row.
        columnIndex[0]=rowIndex ; columnValue[0]=xD[node] ;

        // Set the column entries for edge neighbors.
        const int numEdge=node2edge.num_elems(node) ; int count=1 ;
        for(int i=0;i<numEdge;++i){
          bool iAmNodeZero=(node==edge2node[node2edge[node][i]][0]) ;
          if(iAmNodeZero){
            columnIndex[count]=petscNodeToRow[edge2node[node2edge[node][i]]
              [1]] ;
            columnValue[count]=xE[node2edge[node][i]][0] ; ++count ;
          }else{
            columnIndex[count]=petscNodeToRow[edge2node[node2edge[node][i]]
              [0]] ;
            columnValue[count]=xE[node2edge[node][i]][1] ; ++count ;
          }
        }

        // Insert row entries into the Petsc matrix.
        (*A).SetRowValues(rowIndex,numEdge+1,columnIndex,columnValue) ;
      }

      // Sets matrix values for a node's row.
      void AssembleRHS(Entity node) {
        const PetscScalar valueX=xB[node].x,valueY=xB[node].y,valueZ=
          xB[node].z ;
        (*b).SetValue(&(petscNodeToRow[node]),&valueX,&valueY,&valueZ) ;
      }
  } ;

  register_rule<PETSCNodeSolveVectorUnit> registerPETSCNodeSolveVectorUnit ;

  // Empty apply rule required by Loci. The data type and operator do not
  // matter since nothing is done by this rule. Keep the same inputs and
  // outputs as the unit rule, even though we don't have to.
  class PETSCNodeSolveVectorApply : public apply_rule<blackbox<PETSCVector>,
  Loci::NullOp<PETSCVector> > {
    private:
      const_multiMap node2edge ;
      const_MapVec<2> edge2node ;
      const_storeVec<real> xE ;
      const_store<real> xD ;
      const_store<vect3d> xB ;
      const_store<int> petscNodeToRow ;
      const_blackbox<PETSCMultiVector> b ;
      const_blackbox<PETSCMatrix> A ;
      const_blackbox<PETSCKsp> ksp ;
      blackbox<PETSCMultiVector> phi ;
    public:

      // Define input and output.
      PETSCNodeSolveVectorApply() {
        name_store("node2edge",node2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("X_E",xE) ;
        name_store("X_B",xB) ;
        name_store("X_D",xD) ;
        name_store("petscNodeToRow",petscNodeToRow) ;
        name_store("petscNodeVectorB(X)",b) ;
        name_store("petscNodeA(X)",A) ;
        name_store("petscNodeKSP(X)",ksp) ;
        name_store("petscNodeVectorPhi(X)",phi) ;
        input("node2edge->edge2node,node2edge->X_E,X_B,X_D") ;
        input("node2edge->edge2node->petscNodeToRow") ;
        input("petscNodeA(X),petscNodeVectorB(X),petscNodeKSP(X)") ;
        output("petscNodeVectorPhi(X)") ;
        constraint("X_PETSCLinearSolver,pos") ;
        disable_threading() ;
      }

      // Do nothing.
      void compute(const sequence & seq) {}
  } ;

  register_rule<PETSCNodeSolveVectorApply> registerPETSCNodeSolveVectorApply ;

  // Non-parametric rule for extracting sStar from PETSC.
  class PETSCNodeCopyVector : public pointwise_rule {
    private:
      const_blackbox<PETSCMultiVector> petscSStar ;
      store<vect3d> sStar ;
    public:

      // Define input and output.
      PETSCNodeCopyVector() {
        name_store("petscNodeVectorPhi(sStar)",petscSStar) ;
        name_store("sStar",sStar) ;
        input("petscNodeVectorPhi(sStar)") ;
        output("sStar") ;
        constraint("sStar_PETSCLinearSolver,pos,gridMover") ;
        disable_threading() ;
      }

      // Copy the solution back from PETSC.
      void compute(const sequence & seq) {
        int minRowNum,maxRowNum ;
        (*petscSStar).GetOwnershipRange(&minRowNum,&maxRowNum) ;
        PetscScalar *sStarXCopy,*sStarYCopy,*sStarZCopy ;
        (*petscSStar).GetArray(&sStarXCopy,&sStarYCopy,&sStarZCopy) ;
        sequence::const_iterator nodePtr=seq.begin() ;
        for(int row=minRowNum;row<maxRowNum;++row,++nodePtr) sStar[*nodePtr]=
          vect3d(sStarXCopy[row-minRowNum],sStarYCopy[row-minRowNum],
          sStarZCopy[row-minRowNum]) ;
        (*petscSStar).RestoreArray(&sStarXCopy,&sStarYCopy,&sStarZCopy) ;
      }
  } ;

  register_rule<PETSCNodeCopyVector> registerPETSCNodeCopyVector ;

}


