// Loci includes.
#include <Loci.h>

// StreamUns includes.
#include "sciTypes.h"

// Gauss-Seidel rules for Z.
namespace streamUns {

  // Build rule for the Gauss-Seidel iteration process.
  class BuildGaussSeidelZ : public pointwise_rule {
    private:
      const_store<real> Z ;
      store<real> ZStar ;
    public:

      // Define input and output.
      BuildGaussSeidelZ() {
        name_store("Z{n,it}",Z) ;
        name_store("ZStar{n,it,igsz=0}",ZStar) ;
        input("Z{n,it}") ;
        output("ZStar{n,it,igsz=0}") ;
        constraint("geom_cells{n,it}") ;
      }

      // Initialize value for a single cell.
      void calculate(Entity cell) { ZStar[cell]=Z[cell] ; }

      // Initialize values for a sequence of cells.
      void compute(const sequence &seq) {
        do_loop(seq,this,&BuildGaussSeidelZ::calculate) ;
      }
  } ;

  // Rule for forward pass of Gauss-Seidel iteration process.
  class GaussSeidelZForward : public pointwise_rule {
    private:
      const_store<real> L,U ;
      const_store<real> D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> ZStarForward ;
      store<real> ZStar ;
    public:

      // Define input and output.
      GaussSeidelZForward() {
        set_relaxed_recursion() ;
        name_store("ZStar_L",L) ;
        name_store("ZStar_U",U) ;
        name_store("ZStar_D",D) ;
        name_store("ZStar_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("ZStar",ZStar) ;
        name_store("ZStarForward",ZStarForward) ;
        input("ZStar_D,ZStar_B") ;
        input("upper->ZStar_U,lower->ZStar_L") ;
        input("upper->cr->ZStar") ;
        input("lower->cl->ZStarForward") ;
        output("ZStarForward=ZStar") ;
      }

      // Forward pass for a single cell.
      void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*ZStar[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*ZStarForward[cl[*li]] ;

        // Compute ZStar. This gets copied into ZStarForward later.
        ZStar[cell]=temp*(1.0/D[cell]) ;
      }

      // Forward pass for a sequence of cells.
      void compute(const sequence &seq) {
        do_loop(seq,this,&GaussSeidelZForward::calculate) ;
      }
  } ;

  // Rule for backward pass of Gauss-Seidel iteration process.
  class GaussSeidelZBackward : public pointwise_rule {
    private:
      const_store<real> L,U,D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> ZStarBackward ;
      store<real> ZStarForward ;
    public:

      // Define input and output.
      GaussSeidelZBackward() {
        set_relaxed_recursion() ;
        name_store("ZStar_L",L) ;
        name_store("ZStar_U",U) ;
        name_store("ZStar_D",D) ;
        name_store("ZStar_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("ZStarForward",ZStarForward) ;
        name_store("ZStarBackward",ZStarBackward) ;
        input("ZStar_B,ZStar_D") ;
        input("upper->ZStar_U,lower->ZStar_L") ;
        input("upper->cr->ZStarBackward") ;
        input("lower->cl->ZStarForward") ;
        output("ZStarBackward=ZStarForward") ;
      }

      // Backward pass for a single cell.
      void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*ZStarBackward[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*ZStarForward[cl[*li]] ;

        // Compute ZStarForward. This gets copied into ZStarBackward
        // later.
        ZStarForward[cell]=temp*(1.0/D[cell]) ;
      }

      // Backward pass for a sequence of cells.
      void compute(const sequence &seq) {
        do_loop(seq,this,&GaussSeidelZBackward::calculate) ;
      }
  } ;

  // Rule for advancing the Gauss-Seidel iteration.
  class AdvanceGaussSeidelZ : public pointwise_rule {
    private:
      store<real> ZStar ;
    public:

      // Define input and output.
      AdvanceGaussSeidelZ() {
        name_store("ZStarBackward{n,it,igsz}",ZStar) ;
        input("ZStarBackward{n,it,igsz}") ;
        output("ZStar{n,it,igsz+1}=ZStarBackward{n,it,igsz}") ;
        constraint("geom_cells{n,it}") ;
      }

      virtual void compute(const sequence &seq) {}
  } ;

  // Unit rule to intialize the iteration termination flag.
  class CheckGaussSeidelZUnit : public unit_rule {
    private:
      param<bool> gaussSeidelZFinished ;
    public :
 
      // Define input and output.
      CheckGaussSeidelZUnit() {
        name_store("gaussSeidelZFinished",gaussSeidelZFinished) ;
        constraint("UNIVERSE") ;
        output("gaussSeidelZFinished") ;
      }

      // Initialize the iteration termination flag.
      void compute(const sequence &seq) { *gaussSeidelZFinished=false ; }
  } ;

  // Join operator for the following apply rule.    
  struct chk_join {          
    void operator() (bool &r,const bool &s) { r=r||s ; }
  } ; 
  
  // Apply rule to set termination flag when the maximum number of iterations
  // has been reached.
  class CheckGaussSeidelZMaxIterations : public apply_rule<param<bool>,chk_join> {
    private:
      param<bool> gaussSeidelZFinished ;
      const_param<int> igsz,ZMaxIterations ;
    public :
 
      // Define input and output.
      CheckGaussSeidelZMaxIterations() {
        name_store("gaussSeidelZFinished",gaussSeidelZFinished) ;
        name_store("$igsz",igsz) ;
        name_store("flameletMaxIterations",ZMaxIterations) ;
        input("$igsz,flameletMaxIterations") ;
        input("gaussSeidelZFinished") ;
        output("gaussSeidelZFinished") ;
        constraint("$igsz,flameletMaxIterations") ;
      }

      // Set the iteration termination flag.
      void compute(const sequence &seq) {
        *gaussSeidelZFinished=(*igsz==*ZMaxIterations) ;
      }
  } ;

  // Collapse rule for the Gauss-Seidel iteration.
  class CollapseGaussSeidelZ : public pointwise_rule {
    private:
      store<real> ZStar ;
    public:

      // Define input and output.
      CollapseGaussSeidelZ() {
        name_store("ZStar{n,it,igsz}",ZStar) ;
        input("ZStar{n,it,igsz}") ;
        output("ZStar{n,it}=ZStar{n,it,igsz}") ;
        constraint("geom_cells{n,it}") ;
        conditional("gaussSeidelZFinished{n,it,igsz}") ;
      }

      virtual void compute(const sequence &seq) {}
  } ;

  // Register the rules.
  register_rule<BuildGaussSeidelZ> registerBuildGaussSeidelZ ;
  register_rule<GaussSeidelZForward> registerGaussSeidelZForward ;
  register_rule<GaussSeidelZBackward> registerGaussSeidelZBackward ;
  register_rule<AdvanceGaussSeidelZ> registerAdvanceGaussSeidelZ ;
  register_rule<CheckGaussSeidelZUnit> registerCheckGaussSeidelZUnit ;
  register_rule<CheckGaussSeidelZMaxIterations>
  registerCheckGaussSeidelZMaxIterations ;
  register_rule<CollapseGaussSeidelZ> registerCollapseGaussSeidelZ ;
}



// Gauss-Seidel rules for Zvar.
namespace streamUns {

  // Build rule for the Gauss-Seidel iteration process.
  class BuildGaussSeidelZvar : public pointwise_rule {
    private:
      const_store<real> Zvar ;
      store<real> ZvarStar ;
    public:

      // Define input and output.
      BuildGaussSeidelZvar() {
        name_store("Zvar{n,it}",Zvar) ;
        name_store("ZvarStar{n,it,igszv=0}",ZvarStar) ;
        input("Zvar{n,it}") ;
        output("ZvarStar{n,it,igszv=0}") ;
        constraint("geom_cells{n,it}") ;
      }

      // Initialize value for a single cell.
      void calculate(Entity cell) { ZvarStar[cell]=Zvar[cell] ; }

      // Initialize values for a sequence of cells.
      void compute(const sequence &seq) {
        do_loop(seq,this,&BuildGaussSeidelZvar::calculate) ;
      }
  } ;

  // Rule for forward pass of Gauss-Seidel iteration process.
  class GaussSeidelZvarForward : public pointwise_rule {
    private:
      const_store<real> L,U ;
      const_store<real> D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> ZvarStarForward ;
      store<real> ZvarStar ;
    public:

      // Define input and output.
      GaussSeidelZvarForward() {
        set_relaxed_recursion() ;
        name_store("ZvarStar_L",L) ;
        name_store("ZvarStar_U",U) ;
        name_store("ZvarStar_D",D) ;
        name_store("ZvarStar_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("ZvarStar",ZvarStar) ;
        name_store("ZvarStarForward",ZvarStarForward) ;
        input("ZvarStar_D,ZvarStar_B") ;
        input("upper->ZvarStar_U,lower->ZvarStar_L") ;
        input("upper->cr->ZvarStar") ;
        input("lower->cl->ZvarStarForward") ;
        output("ZvarStarForward=ZvarStar") ;
      }

      // Forward pass for a single cell.
      void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*ZvarStar[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*ZvarStarForward[cl[*li]] ;

        // Compute ZvarStar. This gets copied into ZvarStarForward later.
        ZvarStar[cell]=temp*(1.0/D[cell]) ;
      }

      // Forward pass for a sequence of cells.
      void compute(const sequence &seq) {
        do_loop(seq,this,&GaussSeidelZvarForward::calculate) ;
      }
  } ;

  // Rule for backward pass of Gauss-Seidel iteration process.
  class GaussSeidelZvarBackward : public pointwise_rule {
    private:
      const_store<real> L,U,D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> ZvarStarBackward ;
      store<real> ZvarStarForward ;
    public:

      // Define input and output.
      GaussSeidelZvarBackward() {
        set_relaxed_recursion() ;
        name_store("ZvarStar_L",L) ;
        name_store("ZvarStar_U",U) ;
        name_store("ZvarStar_D",D) ;
        name_store("ZvarStar_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("ZvarStarForward",ZvarStarForward) ;
        name_store("ZvarStarBackward",ZvarStarBackward) ;
        input("ZvarStar_B,ZvarStar_D") ;
        input("upper->ZvarStar_U,lower->ZvarStar_L") ;
        input("upper->cr->ZvarStarBackward") ;
        input("lower->cl->ZvarStarForward") ;
        output("ZvarStarBackward=ZvarStarForward") ;
      }

      // Backward pass for a single cell.
      void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*ZvarStarBackward[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*ZvarStarForward[cl[*li]] ;

        // Compute ZvarStarForward. This gets copied into ZvarStarBackward
        // later.
        ZvarStarForward[cell]=temp*(1.0/D[cell]) ;
      }

      // Backward pass for a sequence of cells.
      void compute(const sequence &seq) {
        do_loop(seq,this,&GaussSeidelZvarBackward::calculate) ;
      }
  } ;

  // Rule for advancing the Gauss-Seidel iteration.
  class AdvanceGaussSeidelZvar : public pointwise_rule {
    private:
      store<real> ZvarStar ;
    public:

      // Define input and output.
      AdvanceGaussSeidelZvar() {
        name_store("ZvarStarBackward{n,it,igszv}",ZvarStar) ;
        input("ZvarStarBackward{n,it,igszv}") ;
        output("ZvarStar{n,it,igszv+1}=ZvarStarBackward{n,it,igszv}") ;
        constraint("geom_cells{n,it}") ;
      }

      virtual void compute(const sequence &seq) {}
  } ;

  // Unit rule to intialize the iteration termination flag.
  class CheckGaussSeidelZvarUnit : public unit_rule {
    private:
      param<bool> gaussSeidelZvarFinished ;
    public :
 
      // Define input and output.
      CheckGaussSeidelZvarUnit() {
        name_store("gaussSeidelZvarFinished",gaussSeidelZvarFinished) ;
        constraint("UNIVERSE") ;
        output("gaussSeidelZvarFinished") ;
      }

      // Initialize the iteration termination flag.
      void compute(const sequence &seq) { *gaussSeidelZvarFinished=false ; }
  } ;
  
  // Apply rule to set termination flag when the maximum number of iterations
  // has been reached.
  class CheckGaussSeidelZvarMaxIterations : public apply_rule<param<bool>,chk_join> {
    private:
      param<bool> gaussSeidelZvarFinished ;
      const_param<int> igszv,ZvarMaxIterations ;
    public :
 
      // Define input and output.
      CheckGaussSeidelZvarMaxIterations() {
        name_store("gaussSeidelZvarFinished",gaussSeidelZvarFinished) ;
        name_store("$igszv",igszv) ;
        name_store("flameletMaxIterations",ZvarMaxIterations) ;
        input("$igszv,flameletMaxIterations") ;
        input("gaussSeidelZvarFinished") ;
        output("gaussSeidelZvarFinished") ;
        constraint("$igszv,flameletMaxIterations") ;
      }

      // Set the iteration termination flag.
      void compute(const sequence &seq) {
        *gaussSeidelZvarFinished=(*igszv==*ZvarMaxIterations) ;
      }
  } ;

  // Collapse rule for the Gauss-Seidel iteration.
  class CollapseGaussSeidelZvar : public pointwise_rule {
    private:
      store<real> ZvarStar ;
    public:

      // Define input and output.
      CollapseGaussSeidelZvar() {
        name_store("ZvarStar{n,it,igszv}",ZvarStar) ;
        input("ZvarStar{n,it,igszv}") ;
        output("ZvarStar{n,it}=ZvarStar{n,it,igszv}") ;
        constraint("geom_cells{n,it}") ;
        conditional("gaussSeidelZvarFinished{n,it,igszv}") ;
      }

      virtual void compute(const sequence &seq) {}
  } ;

  // Register the rules.
  register_rule<BuildGaussSeidelZvar> registerBuildGaussSeidelZvar ;
  register_rule<GaussSeidelZvarForward> registerGaussSeidelZvarForward ;
  register_rule<GaussSeidelZvarBackward> registerGaussSeidelZvarBackward ;
  register_rule<AdvanceGaussSeidelZvar> registerAdvanceGaussSeidelZvar ;
  register_rule<CheckGaussSeidelZvarUnit> registerCheckGaussSeidelZvarUnit ;
  register_rule<CheckGaussSeidelZvarMaxIterations>
  registerCheckGaussSeidelZvarMaxIterations ;
  register_rule<CollapseGaussSeidelZvar> registerCollapseGaussSeidelZvar ;
}
