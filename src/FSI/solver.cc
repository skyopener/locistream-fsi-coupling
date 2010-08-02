//-----------------------------------------------------------------------------
// Description: This file contains rule for a Jacobi solver.
//-----------------------------------------------------------------------------
                                                                                
// Standard library includes.
#include <vector>
using std::vector ;
                                                                                
// Loci includes.
#include <Loci.h>
using Loci::Area ;
                                                                                
// StreamUns includes.
#include "sciTypes.h"

namespace streamUns {

  // Build rule for the Jacobi iteration process.
  class BuildJacobiDisplacement : public pointwise_rule {
    private:
      const_store<vect3d> node_sStar ;
      store<vect3d> sStar ;
    public:

      // Define input and output.
      BuildJacobiDisplacement() {
        name_store("node_sStar",node_sStar) ;
        name_store("sStar{igss=0}",sStar) ;
        input("node_sStar") ;
        output("sStar{igss=0}") ;
        constraint("nodes,gridMover,sStar_JacobiLinearSolver") ;
      }

      // Initialize value for a single node.
      void calculate(Entity node) { sStar[node]=node_sStar[node] ; }

      // Loop over nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BuildJacobiDisplacement> registerBuildJacobiDisplacement ;

  // Rule to initialize the sum.
  class JacobiSumUnit : public unit_rule {
    private:
      store<vect3d> sum ;
    public:

      // Define input and output.
      JacobiSumUnit() {
        name_store("nodeJacobiSum",sum) ;
        output("nodeJacobiSum") ;
        constraint("UNIVERSE,gridMover,sStar_JacobiLinearSolver") ;
      }

      // Initialize.
      void calculate(Entity node) { sum[node]=vect3d(0.0,0.0,0.0) ; }

      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<JacobiSumUnit> registerJacobiSumUnit ;

  // Apply rule to add rhs.
  class JacobiSumNodeApply : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_store<vect3d> sStar_B ;
      store<vect3d> sum ;
    public:

      // Define input and output.
      JacobiSumNodeApply() {
        name_store("sStar_B",sStar_B) ;
        name_store("nodeJacobiSum",sum) ;
        input("sStar_B") ;
        output("nodeJacobiSum") ;
        constraint("nodes,gridMover,sStar_JacobiLinearSolver") ;
      }

      // Add rhs to sum.
      void calculate(Entity node) { sum[node]+=sStar_B[node] ; }

      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<JacobiSumNodeApply> registerJacobiSumNodeApply ;

  // Apply rule to add neighbor contribution.
  class JacobiSumEdgeApply : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_MapVec<2> edge2node ;
      const_store<vect3d> sStar ;
      const_storeVec<real> sStar_E ;
      store<vect3d> sum ;
    public:

      // Define input and output.
      JacobiSumEdgeApply() {
        name_store("edge2node",edge2node) ;
        name_store("sStar",sStar) ;
        name_store("sStar_E",sStar_E) ;
        name_store("nodeJacobiSum",sum) ;
        input("sStar_E,edge2node->sStar") ;
        output("edge2node->nodeJacobiSum") ;
        constraint("edges,gridMover,sStar_JacobiLinearSolver") ;
      }

      // Add rhs to sum.
      void calculate(Entity edge) {
        int n0=edge2node[edge][0],n1=edge2node[edge][1] ;
        sum[n0]-=sStar_E[edge][0]*sStar[n1] ;
        sum[n1]-=sStar_E[edge][1]*sStar[n0] ;
      }

      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<JacobiSumEdgeApply> registerJacobiSumEdgeApply ;

  // Rule for advancing the Jacobi iteration.
  class AdvanceJacobiDisplacement : public pointwise_rule {
    private:
      const_store<real> sStar_D ;
      const_store<vect3d> sum ;
      store<vect3d> sStar ;
    public:

      // Define input and output.
      AdvanceJacobiDisplacement() {
        name_store("sStar_D{igss}",sStar_D) ;
        name_store("nodeJacobiSum{igss}",sum) ;
        name_store("sStar{igss+1}",sStar) ;
        input("sStar_D{igss},nodeJacobiSum{igss}") ;
        output("sStar{igss+1}") ;
        constraint("nodes,gridMover,sStar_JacobiLinearSolver") ;
      }

      // Divide sum by main coefficient for a single node.
      void calculate(Entity node) {
        sStar[node]=(1.0/sStar_D[node])*sum[node] ;
      }

      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<AdvanceJacobiDisplacement> registerAdvanceJacobiDisplacement ;

  // Rule for advancing the Jacobi iteration for nodes constrained not to
  // move in the x-direction.
  class AdvanceJacobiDisplacementSXZero : public pointwise_rule {
    private:
      const_store<real> sStar_D ;
      const_store<vect3d> sum ;
      store<vect3d> sStar ;
    public:

      // Define input and output.
      AdvanceJacobiDisplacementSXZero() {
        name_store("sStar_D{igss}",sStar_D) ;
        name_store("nodeJacobiSum{igss}",sum) ;
        name_store("sXZero::sStar{igss+1}",sStar) ;
        input("sStar_D{igss},nodeJacobiSum{igss}") ;
        output("sXZero::sStar{igss+1}") ;
        constraint("nodes,gridMover,sStar_JacobiLinearSolver") ;
        constraint("sXZeroNodes") ;
      }

      // Divide sum by main coefficient for a single node.
      void calculate(Entity node) {
        sStar[node]=(1.0/sStar_D[node])*sum[node] ;
        sStar[node].x=0.0 ;
      }

      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<AdvanceJacobiDisplacementSXZero>
    registerAdvanceJacobiDisplacementSXZero ;

  // Rule for advancing the Jacobi iteration for nodes constrained not to
  // move in the y-direction.
  class AdvanceJacobiDisplacementSYZero : public pointwise_rule {
    private:
      const_store<real> sStar_D ;
      const_store<vect3d> sum ;
      store<vect3d> sStar ;
    public:

      // Define input and output.
      AdvanceJacobiDisplacementSYZero() {
        name_store("sStar_D{igss}",sStar_D) ;
        name_store("nodeJacobiSum{igss}",sum) ;
        name_store("sYZero::sStar{igss+1}",sStar) ;
        input("sStar_D{igss},nodeJacobiSum{igss}") ;
        output("sYZero::sStar{igss+1}") ;
        constraint("nodes,gridMover,sStar_JacobiLinearSolver") ;
        constraint("sYZeroNodes") ;
      }

      // Divide sum by main coefficient for a single node.
      void calculate(Entity node) {
        sStar[node]=(1.0/sStar_D[node])*sum[node] ;
        sStar[node].y=0.0 ;
      }

      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<AdvanceJacobiDisplacementSYZero>
    registerAdvanceJacobiDisplacementSYZero ;

  // Rule for advancing the Jacobi iteration for nodes constrained not to
  // move in the z-direction.
  class AdvanceJacobiDisplacementSZZero : public pointwise_rule {
    private:
      const_store<real> sStar_D ;
      const_store<vect3d> sum ;
      store<vect3d> sStar ;
    public:

      // Define input and output.
      AdvanceJacobiDisplacementSZZero() {
        name_store("sStar_D{igss}",sStar_D) ;
        name_store("nodeJacobiSum{igss}",sum) ;
        name_store("sZZero::sStar{igss+1}",sStar) ;
        input("sStar_D{igss},nodeJacobiSum{igss}") ;
        output("sZZero::sStar{igss+1}") ;
        constraint("nodes,gridMover,sStar_JacobiLinearSolver") ;
        constraint("sZZeroNodes") ;
      }

      // Divide sum by main coefficient for a single node.
      void calculate(Entity node) {
        sStar[node]=(1.0/sStar_D[node])*sum[node] ;
        sStar[node].z=0.0 ;
      }

      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<AdvanceJacobiDisplacementSZZero>
    registerAdvanceJacobiDisplacementSZZero ;

  // Unit rule to intialize the iteration termination flag.
  class CheckJacobiDisplacementFinishedUnit : public unit_rule {
    private:
      param<bool> jacobiDisplacementFinished ;
    public :

      CheckJacobiDisplacementFinishedUnit() {
        name_store("jacobiDisplacementFinished",jacobiDisplacementFinished) ;
        constraint("UNIVERSE,gridMover,sStar_JacobiLinearSolver") ;
        output("jacobiDisplacementFinished") ;
      }

      // Initialize the iteration termination flag.
      void compute(const sequence &seq) { *jacobiDisplacementFinished=false ; }
  } ;

  register_rule<CheckJacobiDisplacementFinishedUnit>
    registerCheckJacobiDisplacementFinishedUnit ;

  // Join operator for the following apply rule.
  struct chk_join {
    void operator() (bool &r,const bool &s) { r=r||s ; }
  } ;

  // Apply rule to set termination flag when the maximum number of iterations
  // has been reached.
  class CheckJacobiDisplacementFinishedMaxIterations : public
  apply_rule<param<bool>,chk_join> {
    private:
      param<bool> jacobiDisplacementFinished ;
      const_param<int> igss,gridMoverMaxLinearSolverIterations ;
    public :
                                                                                
      // Define input and output.
      CheckJacobiDisplacementFinishedMaxIterations() {
        name_store("jacobiDisplacementFinished",jacobiDisplacementFinished) ;
        name_store("gridMoverMaxLinearSolverIterations",
          gridMoverMaxLinearSolverIterations) ;
        name_store("$igss",igss) ;
        input("$igss,gridMoverMaxLinearSolverIterations") ;
        output("jacobiDisplacementFinished") ;
        constraint("$igss,gridMoverMaxLinearSolverIterations,gridMover") ;
        constraint("sStar_JacobiLinearSolver") ;
      }
                                                                                
      // Set the iteration termination flag.
      void compute(const sequence &seq) {
        *jacobiDisplacementFinished=(*igss==
          *gridMoverMaxLinearSolverIterations) ;
      }
  } ;

  register_rule<CheckJacobiDisplacementFinishedMaxIterations>
    registerCheckJacobiDisplacementFinishedMaxIterations ;

  // Collapse rule for the Jacobi iteration.
  class CollapseJacobiDisplacement : public pointwise_rule {
    private:
      store<vect3d> sStar ;
    public:

      // Define input and output.
      CollapseJacobiDisplacement() {
        name_store("sStar{igss}",sStar) ;
        input("sStar{igss}") ;
        output("sStar=sStar{igss}") ;
        constraint("nodes,gridMover,sStar_JacobiLinearSolver") ;
        conditional("jacobiDisplacementFinished{igss}") ;
      }

      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<CollapseJacobiDisplacement> registerCollapseJacobiDisplacement ;


}
