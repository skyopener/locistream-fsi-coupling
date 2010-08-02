//-----------------------------------------------------------------------------
// Description: This file contains rules for implementing a nodal linear solver
// using HYPRE.
//-----------------------------------------------------------------------------
                                                                                
// StreamUns includes.
#include "hypre.h"

namespace streamUns {

  // Gets the number of nodes assigned on all processes. This rule is
  // definitely not in the Loci spirit since we are collecting data from
  // other processes.
  class HypreGetLocalNodeData : public singleton_rule {
    private:
      param<vector<int> > hypreNumNode ;
      param<vector<int> > hypreNodeRowNumStart,hypreNodeRowNumEnd ;
      param<entitySet> hypreLocalNode ;
    public:

      // Define input and output.
      HypreGetLocalNodeData() {
        name_store("hypreNumNode",hypreNumNode) ;
        name_store("hypreNodeRowNumStart",hypreNodeRowNumStart) ;
        name_store("hypreNodeRowNumEnd",hypreNodeRowNumEnd) ;
        name_store("hypreLocalNode",hypreLocalNode) ;
        output("hypreNumNode,hypreNodeRowNumStart,hypreNodeRowNumEnd") ;
        output("hypreLocalNode") ;
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
        *hypreLocalNode=entitySet(seq) ;

        // Distribute the number of nodes to all processes.
        *hypreNumNode=Loci::all_collect_sizes((*hypreLocalNode).size()) ;

        // Save the starting and ending row numbers.
        *hypreNodeRowNumStart=vector<int>((*hypreNumNode).size()) ;
        *hypreNodeRowNumEnd=vector<int>((*hypreNumNode).size()) ;
        (*hypreNodeRowNumStart)[Loci::MPI_rank]=0 ;
        for(int i=0;i<Loci::MPI_rank;++i)
          (*hypreNodeRowNumStart)[Loci::MPI_rank]+=(*hypreNumNode)[i] ;
        (*hypreNodeRowNumEnd)[Loci::MPI_rank]=(*hypreNodeRowNumStart)
          [Loci::MPI_rank]+(*hypreNumNode)[Loci::MPI_rank]-1 ;
      }
  } ;

  register_rule<HypreGetLocalNodeData> registerHypreGetLocalNodeData ;

  // Creates the node-to-row map.
  class HypreNodeToRow : public pointwise_rule {
    private:
      const_param<vector<int> > hypreNumNode ;
      store<int> hypreNodeToRow ;
    public:
                                                                                
      // Define input and output.
      HypreNodeToRow() {
        name_store("hypreNumNode",hypreNumNode) ;
        name_store("hypreNodeToRow",hypreNodeToRow) ;
        input("hypreNumNode") ;
        output("hypreNodeToRow") ;
        constraint("pos") ;
        disable_threading() ;
      }
                                                                                
      // Set the Hypre row for each node.
      virtual void compute(const sequence & seq) {
                                                                                
        // Compute the row offset for this process.
        int offset=0 ;
        for(int i=0;i<Loci::MPI_rank;++i){ offset+=(*hypreNumNode)[i] ; }
                                                                                
        // Assign row number.
        sequence::const_iterator nodePtr=seq.begin() ;
        for(int i=0;i<(*hypreNumNode)[Loci::MPI_rank];++nodePtr,++i){
          hypreNodeToRow[*nodePtr]=offset+i ;
        }
      }
  } ;
                                                                                
  register_rule<HypreNodeToRow> registerHypreNodeToRow ;

  // Determines the number of non-zero entries in the local portion of the
  // HYPRE matrix for each node. The local portion is defined as the square
  // row/column sub-block of the HYPRE matrix whose rows map to nodes local
  // to the process.
  class HypreNodeNumDiagonalNonZeroEntries : public pointwise_rule {
    private:
      const_param<entitySet> hypreLocalNode ;
      const_multiMap node2edge ;
      const_multiMap edge2node ;
      store<int> hypreNodeNumDiagonalNonZero ;
    public:

      // Define input and output.
      HypreNodeNumDiagonalNonZeroEntries() {
        name_store("hypreLocalNode",hypreLocalNode) ;
        name_store("node2edge",node2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("hypreNodeNumDiagonalNonZero",hypreNodeNumDiagonalNonZero) ;
        input("hypreLocalNode,node2edge->edge2node") ;
        output("hypreNodeNumDiagonalNonZero") ;
        constraint("pos") ;
      }

      // Count local nodes for this node.
      void calculate(Entity node) {

        // Initialize count with the diagonal entry.
        hypreNodeNumDiagonalNonZero[node]=1 ;

        // Add neighbor node for each edge if it is local.
        int numEdge=node2edge.num_elems(node) ;
        for(int i=0;i<numEdge;++i) {
          if((*hypreLocalNode).inSet((node==edge2node[node2edge[node][i]][0])?
          edge2node[node2edge[node][i]][1]:edge2node[node2edge[node][i]][0]))
            ++hypreNodeNumDiagonalNonZero[node] ;
        }
      }

      // Loop over nodes.
      void compute(const sequence & seq) { do_loop(seq,this) ; }

  } ;

  register_rule<HypreNodeNumDiagonalNonZeroEntries>
    registerHypreNodeNumDiagonalNonZeroEntries ;

  // Determines the number of non-zero entries in the non-local portion of the
  // HYPRE matrix for each node.
  class HypreNodeNumOffDiagonalNonZeroEntries : public pointwise_rule {
    private:
      const_param<entitySet> hypreLocalNode ;
      const_MapVec<2> edge2node ;
      const_multiMap node2edge ;
      const_store<int> hypreNodeNumDiagonalNonZero ;
      store<int> hypreNodeNumOffDiagonalNonZero ;
    public:

      // Define input and output.
      HypreNodeNumOffDiagonalNonZeroEntries() {
        name_store("hypreLocalNode",hypreLocalNode) ;
        name_store("edge2node",edge2node) ;
        name_store("node2edge",node2edge) ;
        name_store("hypreNodeNumOffDiagonalNonZero",
          hypreNodeNumOffDiagonalNonZero) ;
        input("hypreLocalNode,node2edge") ;
        output("hypreNodeNumOffDiagonalNonZero") ;
        constraint("pos") ;
      }

      // Count non-local nodes for this node.
      void calculate(Entity node) {

        // Initialize count.
        hypreNodeNumOffDiagonalNonZero[node]=0 ;

        // Add neighbor node for each edge if it is not local.
        int numEdge=node2edge.num_elems(node) ;
        for(int i=0;i<numEdge;++i) {
          if(!(*hypreLocalNode).inSet((node==edge2node[node2edge[node][i]][0])?
            edge2node[node2edge[node][i]][1]:edge2node[node2edge[node][i]][0]))
            ++hypreNodeNumOffDiagonalNonZero[node] ;
        }
      }

      // Loop over nodes.
      void compute(const sequence & seq) { do_loop(seq,this) ; }
  } ;

  register_rule<HypreNodeNumOffDiagonalNonZeroEntries>
    registerHypreNodeNumOffDiagonalNonZeroEntries ;

  // Sets up the HYPRE right-hand-side vector.
  class HypreNodeSetupRHS : public singleton_rule {
    private:
      const_param<vector<int> > hypreNodeRowNumStart,hypreNodeRowNumEnd ;
      blackbox<HypreVector> b ;
    public:

      // Define input and output.
      HypreNodeSetupRHS() {
        name_store("hypreNodeRowNumStart",hypreNodeRowNumStart) ;
        name_store("hypreNodeRowNumEnd",hypreNodeRowNumEnd) ;
        name_store("hypreNodeB(X)",b) ;
        input("hypreNodeRowNumStart,hypreNodeRowNumEnd") ;
        output("hypreNodeB(X)") ;
        constraint("X_HYPRELinearSolver,pos") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {
        (*b).Create((*hypreNodeRowNumStart)[Loci::MPI_rank],
          (*hypreNodeRowNumEnd)[Loci::MPI_rank]) ;
      }
  } ;

  register_rule<HypreNodeSetupRHS> registerHypreNodeSetupRHS ;

  // Sets up the HYPRE matrix. Note that this is a unit_rule, since this is the
  // only way we can have stores as input to a rule outputting blackboxes.
  class HypreNodeSetupMatrixUnit : public unit_rule {
    private:
      const_store<int> hypreNodeToRow ;
      const_param<vector<int> > hypreNodeRowNumStart,hypreNodeRowNumEnd ;
      const_store<int> hypreNodeNumDiagonalNonZero ;
      const_store<int> hypreNodeNumOffDiagonalNonZero ;
      const_param<vector<int> > hypreNumNode ;
      blackbox<HypreMatrix> A ;
    public:

      // Define input and output.
      HypreNodeSetupMatrixUnit() {
        name_store("hypreNodeToRow",hypreNodeToRow) ;
        name_store("hypreNodeRowNumStart",hypreNodeRowNumStart) ;
        name_store("hypreNodeRowNumEnd",hypreNodeRowNumEnd) ;
        name_store("hypreNodeNumDiagonalNonZero",hypreNodeNumDiagonalNonZero) ;
        name_store("hypreNodeNumOffDiagonalNonZero",
          hypreNodeNumOffDiagonalNonZero) ;
        name_store("hypreNumNode",hypreNumNode) ;
        name_store("hypreNodeA(X)",A) ;
        input("hypreNodeToRow,hypreNodeRowNumStart,hypreNodeRowNumEnd") ;
        input("hypreNodeNumDiagonalNonZero,hypreNodeNumOffDiagonalNonZero") ;
        input("hypreNumNode") ;
        output("hypreNodeA(X)") ;
        constraint("X_HYPRELinearSolver,pos") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {

        // Make temporary copy of matrix allocation data.
        int count=0,*numDiagonalNonZero=hypre_CTAlloc(int,seq.size()),
          *numOffDiagonalNonZero=hypre_CTAlloc(int,seq.size()) ;
        for(sequence::const_iterator nodePtr=seq.begin();nodePtr!=seq.end();
        ++nodePtr,++count){
          numDiagonalNonZero[count]=hypreNodeNumDiagonalNonZero[*nodePtr] ;
          numOffDiagonalNonZero[count]=hypreNodeNumOffDiagonalNonZero
            [*nodePtr] ;
        }
                                                                                
        // Allocate the matrix.
        (*A).Create((*hypreNodeRowNumStart)[Loci::MPI_rank],
          (*hypreNodeRowNumEnd)[Loci::MPI_rank],numDiagonalNonZero,
          numOffDiagonalNonZero) ;
                                                                                
        // Deallocate temporary copy of matrix allocation data.
        hypre_TFree(numDiagonalNonZero) ; hypre_TFree(numOffDiagonalNonZero) ;
      }
  } ;
                                                                                
  register_rule<HypreNodeSetupMatrixUnit> registerHypreNodeSetupMatrixUnit ;

  // Empty apply rule required by Loci. The data type and operator do not
  // matter since nothing is done by this rule. Keep the same inputs and
  // outputs as the unit rule, even though we don't have to.
  class HypreNodeSetupMatrixApply : public apply_rule<blackbox<HypreMatrix>,
  Loci::NullOp<HypreMatrix> > {
    private:
      const_store<int> hypreNodeToRow ;
      const_store<int> hypreNodeNumDiagonalNonZero ;
      const_store<int> hypreNodeNumOffDiagonalNonZero ;
      const_param<vector<int> > hypreNumNode ;
      blackbox<HypreMatrix> A ;
    public:
                                                                                
      // Define input and output.
      HypreNodeSetupMatrixApply() {
        name_store("hypreNodeToRow",hypreNodeToRow) ;
        name_store("hypreNodeNumDiagonalNonZero",hypreNodeNumDiagonalNonZero) ;
        name_store("hypreNodeNumOffDiagonalNonZero",
          hypreNodeNumOffDiagonalNonZero) ;
        name_store("hypreNumNode",hypreNumNode) ;
        name_store("hypreNodeA(X)",A) ;
        input("hypreNodeToRow,hypreNodeNumDiagonalNonZero") ;
        input("hypreNodeNumOffDiagonalNonZero,hypreNumNode") ;
        output("hypreNodeA(X)") ;
        constraint("X_HYPRELinearSolver,pos") ;
        disable_threading() ;
      }
                                                                                
      // Do nothing.
      void compute(const sequence & seq) {}
  } ;
                                                                                
  register_rule<HypreNodeSetupMatrixApply> registerHypreNodeSetupMatrixApply ;

  // Sets up the HYPRE linear solver.
  class HypreNodeSetupSolver : public singleton_rule {
    private:
      const_param<int> maxLinearSolverIterations ;
      const_param<real> linearSolverTolerance ;
      blackbox<HypreSolver> hypreSolver ;
    public:

      // Define input and output.
      HypreNodeSetupSolver() {
        name_store("X_maxLinearSolverIterations",maxLinearSolverIterations) ;
        name_store("X_linearSolverTolerance",linearSolverTolerance) ;
        name_store("hypreNodeSolver(X)",hypreSolver) ;
        input("X_maxLinearSolverIterations,X_linearSolverTolerance") ;
        output("hypreNodeSolver(X)") ;
        constraint("X_HYPRELinearSolver,pos") ;
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {
        (*hypreSolver).Create() ;
        (*hypreSolver).SetCoarsenType(6) ;
        (*hypreSolver).SetCycleNumSweeps(2,2,2,1) ;
        (*hypreSolver).SetCycleRelaxType(3,3,3,3) ;
        (*hypreSolver).SetMaxIterations(*maxLinearSolverIterations) ;
        (*hypreSolver).SetMaxLevels(16) ;
        (*hypreSolver).SetStrongThreshold(0.5) ;
        (*hypreSolver).SetTolerance(*linearSolverTolerance) ;
      }
  } ;
                                                                                
  register_rule<HypreNodeSetupSolver> registerHypreNodeSetupSolver ;

  // Assemble and solve the HYPRE system.
  class HypreNodeSolveUnit : public unit_rule {
    private:
      const_param<vector<int> > hypreNodeRowNumStart,hypreNodeRowNumEnd ;
      const_multiMap node2edge ;
      const_MapVec<2> edge2node ;
      const_storeVec<real> xE ;
      const_store<real> xB,xD ;
      const_store<int> hypreNodeToRow ;
      const_blackbox<HypreVector> b ;
      const_blackbox<HypreMatrix> A ;
      const_blackbox<HypreSolver> hypreSolver ;
      blackbox<HypreVector> phi ;
    public:

      // Define input and output.
      HypreNodeSolveUnit() {
        name_store("hypreNodeRowNumStart",hypreNodeRowNumStart) ;
        name_store("hypreNodeRowNumEnd",hypreNodeRowNumEnd) ;
        name_store("node2edge",node2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("X_E",xE) ;
        name_store("X_B",xB) ;
        name_store("X_D",xD) ;
        name_store("hypreNodeToRow",hypreNodeToRow) ;
        name_store("hypreNodeB(X)",b) ;
        name_store("hypreNodeA(X)",A) ;
        name_store("hypreNodeSolver(X)",hypreSolver) ;
        name_store("hypreNodePhi(X)",phi) ;
        input("hypreNodeRowNumStart,hypreNodeRowNumEnd") ;
        input("node2edge->edge2node,node2edge->X_E,X_B,X_D") ;
        input("node2edge->edge2node->hypreNodeToRow") ;
        input("hypreNodeA(X),hypreNodeB(X),hypreNodeSolver(X)") ;
        output("hypreNodePhi(X)") ;
        constraint("X_HYPRELinearSolver,pos") ;
        disable_threading() ;
      }

      // Assemble and solve.
      virtual void compute(const sequence &seq) {

        // Create and initialize the solution vector.
        (*phi).Create((*hypreNodeRowNumStart)[Loci::MPI_rank],
          (*hypreNodeRowNumEnd)[Loci::MPI_rank]) ;
        InitializeSolutionVector(seq) ;

        // Assemble the matrix.
        AssembleMatrix(seq) ;
//cout << "AssembleMatrix" << endl ;
        // Assemble the rhs vector.
        AssembleRHS(seq) ;
//cout << "AssembleRHS" << endl ;
        // Solve the linear system
        (*hypreSolver).Solve(*A,*b,*phi) ;
//cout << "Solve" << endl ;
      }

      // Sets matrix values.
      void AssembleMatrix(const sequence &seq) {

        // Initialize the matrix.
        (*A).Initialize() ;

        // Loop over nodes. Note the fixed dimension of 100. There should
        // always be way less than 100 neighbors for any cell. This may not
        // always be true as in the case of a pipe grid where the centerline
        // is degenerate. Then the nodes on the centerline have as many
        // neighbors as in the theta direction.
        int numRow=1 ;
        int *colNum=hypre_CTAlloc(int,100) ;
        double *colValue=hypre_CTAlloc(double,100) ;
        sequence::const_iterator nodePtr=seq.begin() ;
        for(int j=0;j<seq.size();++j,++nodePtr){

          // Get the row number for this node.
          int rowNum[1]={hypreNodeToRow[*nodePtr]} ;

          // Set the diagonal entry for the row.
          colNum[0]=rowNum[0] ; colValue[0]=xD[*nodePtr] ;

          // Set the column entries for edge neighbors.
          const int numEdge=node2edge.num_elems(*nodePtr) ; int count=1 ;
          for(int i=0;i<numEdge;++i){
            bool iAmNodeZero=(*nodePtr==edge2node[node2edge[*nodePtr][i]][0]) ;
            if(iAmNodeZero){
              colNum[count]=hypreNodeToRow[edge2node[node2edge[*nodePtr][i]]
                [1]] ;
              colValue[count]=xE[node2edge[*nodePtr][i]][0] ; ++count ;
            }else{
              colNum[count]=hypreNodeToRow[edge2node[node2edge[*nodePtr][i]]
                [0]] ;
              colValue[count]=xE[node2edge[*nodePtr][i]][1] ; ++count ;
            }
          }

          // Insert row entries into the HYPRE matrix.
          int numCol[1]={numEdge+1} ;
          (*A).SetValues(numRow,numCol,rowNum,colNum,colValue) ;
        }
        hypre_TFree(colNum) ; hypre_TFree(colValue) ;

        // Assemble the matrix.
        (*A).Assemble() ;
      }

      // Sets right-hand side values.
      void AssembleRHS(const sequence &seq) {
        (*b).Initialize() ; double *v=hypre_CTAlloc(double,seq.size()) ;
        sequence::const_iterator nodePtr=seq.begin() ;
        for(int i=0;i<seq.size();++i,++nodePtr) v[i]=xB[*nodePtr] ;
        (*b).SetValues(v) ; (*b).Assemble() ; hypre_TFree(v) ;
      }

      // Initializes the unknown vector.
      void InitializeSolutionVector(const sequence &seq) {
        (*phi).Initialize() ; double *v=hypre_CTAlloc(double,seq.size()) ;
        for(int i=0;i<seq.size();++i) v[i]=0.0 ;
        (*phi).SetValues(v) ; (*phi).Assemble() ; hypre_TFree(v) ;
      }
  } ;

  register_rule<HypreNodeSolveUnit> registerHypreNodeSolveUnit ;

  // Empty apply rule required by Loci. The data type and operator do not
  // matter since nothing is done by this rule. Keep the same inputs and
  // outputs as the unit rule, even though we don't have to.
  class HypreNodeSolveApply : public apply_rule<blackbox<HypreVector>,
  Loci::NullOp<HypreVector> > {
    private:
      const_param<vector<int> > hypreNodeRowNumStart,hypreNodeRowNumEnd ;
      const_multiMap node2edge ;
      const_MapVec<2> edge2node,xE ;
      const_store<real> xB,xD ;
      const_store<int> hypreNodeToRow ;
      const_blackbox<HypreVector> b ;
      const_blackbox<HypreMatrix> A ;
      const_blackbox<HypreSolver> hypreSolver ;
      blackbox<HypreVector> phi ;
    public:

      // Define input and output.
      HypreNodeSolveApply() {
        name_store("hypreNodeRowNumStart",hypreNodeRowNumStart) ;
        name_store("hypreNodeRowNumEnd",hypreNodeRowNumEnd) ;
        name_store("node2edge",node2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("X_E",xE) ;
        name_store("X_B",xB) ;
        name_store("X_D",xD) ;
        name_store("hypreNodeToRow",hypreNodeToRow) ;
        name_store("hypreNodeB(X)",b) ;
        name_store("hypreNodeA(X)",A) ;
        name_store("hypreNodeSolver(X)",hypreSolver) ;
        name_store("hypreNodePhi(X)",phi) ;
        input("hypreNodeRowNumStart,hypreNodeRowNumEnd") ;
        input("node2edge->edge2node,node2edge->X_E,X_B,X_D") ;
        input("node2edge->edge2node->hypreNodeToRow") ;
        input("hypreNodeA(X),hypreNodeB(X),hypreNodeSolver(X)") ;
        output("hypreNodePhi(X)") ;
        constraint("X_HYPRELinearSolver,pos") ;
        disable_threading() ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;
                                                                                
  register_rule<HypreNodeSolveApply> registerHypreNodeSolveApply ;

  // Non-parametric rule for extracting sStar from HYPRE.
  class HypreNodeCopy : public pointwise_rule {
    private:
      const_store<int> hypreNodeToRow ;
      const_blackbox<HypreVector> hypreSStar ;
      store<real> sStar ;
    public:
                                                                                
      // Define input and output.
      HypreNodeCopy() {
        name_store("hypreNodeToRow",hypreNodeToRow) ;
        name_store("hypreNodePhi(sStar)",hypreSStar) ;
        name_store("sStar",sStar) ;
        input("hypreNodeToRow,hypreNodePhi(sStar)") ;
        output("sStar") ;
        constraint("s_HYPRELinearSolver,pos,gridMover") ;
        disable_threading() ;
      }
                                                                                
      // Copy the solution back from HYPRE.
      void compute(const sequence &seq) {
        double *v=hypre_CTAlloc(double,seq.size()) ;
        (*hypreSStar).GetValues(v) ;
        sequence::const_iterator nodePtr=seq.begin() ;
        for(int i=0;i<seq.size();++i,++nodePtr) sStar[*nodePtr]=v[i] ;
        hypre_TFree(v) ;
      }
  } ;
                                                                                
  register_rule<HypreNodeCopy> registerHypreNodeCopy ;
}
