// Loci includes.
#include <Loci.h>

// StreamUns includes.
#include "sciTypes.h"

// Gauss-Seidel rules for Alpha.
namespace streamUns {

  // Build rule for the Gauss-Seidel iteration process.
  class BuildGaussSeidelAlpha : public pointwise_rule {
    private:
      const_store<real> Alpha ;
      store<real> AlphaStar ;
    public:

      // Define input and output.
      BuildGaussSeidelAlpha() {
        name_store("Alpha{n,it}",Alpha) ;
        name_store("AlphaStar{n,it,igsz=0}",AlphaStar) ;
        input("Alpha{n,it}") ;
        output("AlphaStar{n,it,igsz=0}") ;
        constraint("geom_cells{n,it}") ;
      }

      // Initialize value for a single cell.
      void calculate(Entity cell) { AlphaStar[cell]=Alpha[cell] ; }

      // Initialize values for a sequence of cells.
      void compute(const sequence &seq) {
        do_loop(seq,this,&BuildGaussSeidelAlpha::calculate) ;
      }
  } ;
  
  // Rule for forward pass of Gauss-Seidel iteration process.
  class GaussSeidelAlphaForward : public pointwise_rule {
    private:
      const_store<real> L,U ;
      const_store<real> D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> AlphaStarForward ;
      store<real> AlphaStar ;
    public:

      // Define input and output.
      GaussSeidelAlphaForward() {
        set_relaxed_recursion() ;
        name_store("AlphaStar_L",L) ;
        name_store("AlphaStar_U",U) ;
        name_store("AlphaStar_D",D) ;
        name_store("AlphaStar_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("AlphaStar",AlphaStar) ;
        name_store("AlphaStarForward",AlphaStarForward) ;
        input("AlphaStar_D,AlphaStar_B") ;
        input("upper->AlphaStar_U,lower->AlphaStar_L") ;
        input("upper->cr->AlphaStar") ;
        input("lower->cl->AlphaStarForward") ;
        output("AlphaStarForward=AlphaStar") ;
      }

      // Forward pass for a single cell.
      void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*AlphaStar[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*AlphaStarForward[cl[*li]] ;

        // Compute AlphaStar. This gets copied into AlphaStarForward later.
        AlphaStar[cell]=temp*(1.0/D[cell]) ;
      }

      // Forward pass for a sequence of cells.
      void compute(const sequence &seq) {
        do_loop(seq,this,&GaussSeidelAlphaForward::calculate) ;
      }
  } ;

  // Rule for backward pass of Gauss-Seidel iteration process.
  class GaussSeidelAlphaBackward : public pointwise_rule {
    private:
      const_store<real> L,U,D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> AlphaStarBackward ;
      store<real> AlphaStarForward ;
    public:

      // Define input and output.
      GaussSeidelAlphaBackward() {
        set_relaxed_recursion() ;
        name_store("AlphaStar_L",L) ;
        name_store("AlphaStar_U",U) ;
        name_store("AlphaStar_D",D) ;
        name_store("AlphaStar_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("AlphaStarForward",AlphaStarForward) ;
        name_store("AlphaStarBackward",AlphaStarBackward) ;
        input("AlphaStar_B,AlphaStar_D") ;
        input("upper->AlphaStar_U,lower->AlphaStar_L") ;
        input("upper->cr->AlphaStarBackward") ;
        input("lower->cl->AlphaStarForward") ;
        output("AlphaStarBackward=AlphaStarForward") ;
      }

      // Backward pass for a single cell.
      void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*AlphaStarBackward[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*AlphaStarForward[cl[*li]] ;

        // Compute AlphaStarForward. This gets copied into AlphaStarBackward
        // later.
        AlphaStarForward[cell]=temp*(1.0/D[cell]) ;
      }

      // Backward pass for a sequence of cells.
      void compute(const sequence &seq) {
        do_loop(seq,this,&GaussSeidelAlphaBackward::calculate) ;
      }
  } ;

  // Rule for advancing the Gauss-Seidel iteration.
  class AdvanceGaussSeidelAlpha : public pointwise_rule {
    private:
      store<real> AlphaStar ;
    public:

      // Define input and output.
      AdvanceGaussSeidelAlpha() {
        name_store("AlphaStarBackward{n,it,igsz}",AlphaStar) ;
        input("AlphaStarBackward{n,it,igsz}") ;
        output("AlphaStar{n,it,igsz+1}=AlphaStarBackward{n,it,igsz}") ;
        constraint("geom_cells{n,it}") ;
      }

      virtual void compute(const sequence &seq) {}
  } ;

  // Unit rule to intialize the iteration termination flag.
  class CheckGaussSeidelAlphaUnit : public unit_rule {
    private:
      param<bool> gaussSeidelAlphaFinished ;
    public :
 
      // Define input and output.
      CheckGaussSeidelAlphaUnit() {
        name_store("gaussSeidelAlphaFinished",gaussSeidelAlphaFinished) ;
        constraint("UNIVERSE") ;
        output("gaussSeidelAlphaFinished") ;
      }

      // Initialize the iteration termination flag.
      void compute(const sequence &seq) { *gaussSeidelAlphaFinished=false ; }
  } ;

  // Join operator for the following apply rule.    
  struct chk_join {          
    void operator() (bool &r,const bool &s) { r=r||s ; }
  } ; 
  
  // Apply rule to set termination flag when the maximum number of iterations
  // has been reached.
  class CheckGaussSeidelAlphaMaxIterations : public apply_rule<param<bool>,chk_join> {
    private:
      param<bool> gaussSeidelAlphaFinished ;
      const_param<int> igsz,AlphaMaxIterations ;
    public :
 
      // Define input and output.
      CheckGaussSeidelAlphaMaxIterations() {
        name_store("gaussSeidelAlphaFinished",gaussSeidelAlphaFinished) ;
        name_store("$igsz",igsz) ;
        name_store("cavitationMaxIterations",AlphaMaxIterations) ;
//        input("$igsz,cavitationMaxIterations") ;
//        input("gaussSeidelAlphaFinished") ;
        input("$igsz,cavitationMaxIterations") ;
        output("gaussSeidelAlphaFinished") ;
        constraint("$igsz,cavitationMaxIterations") ;
      }

      // Set the iteration termination flag.
      void compute(const sequence &seq) {
        *gaussSeidelAlphaFinished=(*igsz==*AlphaMaxIterations) ;
      }
  } ;

  // Collapse rule for the Gauss-Seidel iteration.
  class CollapseGaussSeidelAlpha : public pointwise_rule {
    private:
      store<real> AlphaStar ;
    public:

      // Define input and output.
      CollapseGaussSeidelAlpha() {
        name_store("AlphaStar{n,it,igsz}",AlphaStar) ;
        input("AlphaStar{n,it,igsz}") ;
        output("AlphaStar{n,it}=AlphaStar{n,it,igsz}") ;
        constraint("geom_cells{n,it}") ;
        conditional("gaussSeidelAlphaFinished{n,it,igsz}") ;
      }

      virtual void compute(const sequence &seq) {}
  } ;

  // Register the rules.
  register_rule<BuildGaussSeidelAlpha> registerBuildGaussSeidelAlpha ;
  register_rule<GaussSeidelAlphaForward> registerGaussSeidelAlphaForward ;
  register_rule<GaussSeidelAlphaBackward> registerGaussSeidelAlphaBackward ;
  register_rule<AdvanceGaussSeidelAlpha> registerAdvanceGaussSeidelAlpha ;
  register_rule<CheckGaussSeidelAlphaUnit> registerCheckGaussSeidelAlphaUnit ;
  register_rule<CheckGaussSeidelAlphaMaxIterations>
  registerCheckGaussSeidelAlphaMaxIterations ;
  register_rule<CollapseGaussSeidelAlpha> registerCollapseGaussSeidelAlpha ;
  
}
