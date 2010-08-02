// Loci includes.
#include <Loci.h>

// StreamUns includes.
#include "sciTypes.h"

// Gauss-Seidel rules for vStar.
namespace streamUns {

  // Build rule for the Gauss-Seidel iteration process.
  class BuildGaussSeidelVector : public pointwise_rule {
    private:
      const_store<vect3d> v ;
      store<vect3d> vStar ;
    public:

      // Define input and output.
      BuildGaussSeidelVector() {
        name_store("v",v) ;
        name_store("vStar{igsv=0}",vStar) ;
        input("v") ;
        output("vStar{igsv=0}") ;
        constraint("geom_cells") ;
      }

      // Initialize value for a single cell.
      inline void calculate(Entity cell) {
        vStar[cell]=v[cell] ;
      }

      // Initialize values for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BuildGaussSeidelVector> registerBuildGaussSeidelVector ;

  // Rule for forward pass of Gauss-Seidel iteration process. Note that the
  // rule to handle periodic boundaries is located in periodic.loci, since
  // we also need vStar in the ghost cells for other rules.
  class GaussSeidelVectorForward : public pointwise_rule {
    private:
      const_store<real> L,U ;
      const_store<real> D ;
      const_store<vect3d> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<vect3d> vStarForward ;
      store<vect3d> vStar ;
    public:

      // Define input and output.
      GaussSeidelVectorForward() {
        name_store("vStar_L",L) ;
        name_store("vStar_U",U) ;
        name_store("vStar_D",D) ;
        name_store("vStar_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("vStar",vStar) ;
        name_store("vStarForward",vStarForward) ;
        input("vStar_D,vStar_B") ;
        input("upper->vStar_U,lower->vStar_L") ;
        input("upper->cr->vStar") ;
        input("lower->cl->vStarForward") ;
        output("vStarForward=vStar") ;
        set_relaxed_recursion() ;
      }

      // Forward pass for a sequence of cells.
      void compute(const sequence &seq) {

        for(sequence::const_iterator si=seq.begin();si!=seq.end();++si) {

          // Save the right-hand side.
          const Entity cell=*si ; vect3d temp=B[cell] ;

// Debug
//cout << "cell,D,-(L+U)" << cell << " " << D[cell] << " " ;
//real sum=0.0 ;
//for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
//sum-=U[*ui] ;
//for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
//sum-=L[*li] ;
//cout << sum ;
//if(D[cell]<sum) cout << " ERROR" ; cout << endl ;
// End debug

          // Add the contributions from the lower-numbered cells.
          for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
            temp-=U[*ui]*vStar[cr[*ui]] ;

          // Add the contributions from the higher-numbered cells.
          for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
            temp-=L[*li]*vStarForward[cl[*li]] ;

          // Compute vStar. This gets copied into vStarForward later.
          vStar[cell]=temp*(1.0/D[cell]) ;
        }
      }
  } ;

  register_rule<GaussSeidelVectorForward> registerGaussSeidelVectorForward ;

  // Rule for backward pass of Gauss-Seidel iteration process.
  class GaussSeidelVectorBackward : public pointwise_rule {
    private:
      const_store<real> L,U,D ;
      const_store<vect3d> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<vect3d> vStarBackward ;
      store<vect3d> vStarForward ;
    public:

      // Define input and output.
      GaussSeidelVectorBackward() {
        name_store("vStar_L",L) ;
        name_store("vStar_U",U) ;
        name_store("vStar_D",D) ;
        name_store("vStar_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("vStarForward",vStarForward) ;
        name_store("vStarBackward",vStarBackward) ;
        input("vStar_B,vStar_D") ;
        input("upper->vStar_U,lower->vStar_L") ;
        input("upper->cr->vStarBackward") ;
        input("lower->cl->vStarForward") ;
        output("vStarBackward=vStarForward") ;
        set_relaxed_recursion() ;
      }

      // Backward pass for a sequence of cells.
      void compute(const sequence &seq) {

        for(sequence::const_iterator si=seq.begin();si!=seq.end();++si) {

          // Save the right-hand side.
          const Entity cell=*si ; vect3d temp=B[cell] ;

          // Add the contributions from the lower-numbered cells.
          for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
            temp-=U[*ui]*vStarBackward[cr[*ui]] ;

          // Add the contributions from the higher-numbered cells.
          for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
            temp-=L[*li]*vStarForward[cl[*li]] ;

          // Compute vStarForward. This gets copied into vStarBackward later.
          vStarForward[cell]=temp*(1.0/D[cell]) ;
        }
      }
  } ;

  register_rule<GaussSeidelVectorBackward> registerGaussSeidelVectorBackward ;

  // Rule for backward pass of Gauss-Seidel iteration process for periodic
  // boundaries.
  class GaussSeidelVectorBackwardPeriodic : public pointwise_rule {
    private:
      const_Map cl,cr,pmap ;
      store<vect3d> vStarForward ;
    public:

      // Define input and output.
      GaussSeidelVectorBackwardPeriodic() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pmap",pmap) ;
        name_store("vStarForward",vStarForward) ;
        input("pmap->cl->vStarForward") ;
        output("cr->(vStarBackward=vStarForward)") ;
        constraint("periodicFaces") ;
      }

      // For a single face.
      void calculate(Entity face) {
        vStarForward[cr[face]]=vStarForward[cl[pmap[face]]] ; }

      // Loop over faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelVectorBackwardPeriodic>
    registerGaussSeidelVectorBackwardPeriodic ;

  // Rule for advancing the Gauss-Seidel iteration.
  class AdvanceGaussSeidelVector : public pointwise_rule {
    private:
      store<vect3d> vStar ;
    public:

      // Define input and output.
      AdvanceGaussSeidelVector() {
        name_store("vStarBackward{igsv}",vStar) ;
        input("vStarBackward{igsv}") ;
        output("vStar{igsv+1}=vStarBackward{igsv}") ;
        constraint("geom_cells") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<AdvanceGaussSeidelVector> registerAdvanceGaussSeidelVector ;

  // Unit rule to intialize the iteration termination flag.
  class CheckGaussSeidelVectorUnit : public unit_rule {
    private:
      param<bool> gaussSeidelVectorFinished ;
    public :
 
      // Define input and output. Replaced the UNIVERSE constraint with one for
      // vMaxIterations so this rule doesn't get stranded for cases where we are
      // not solving for velocity.
      CheckGaussSeidelVectorUnit() {
        name_store("gaussSeidelVectorFinished",gaussSeidelVectorFinished) ;
        constraint("vMaxIterations") ;
        output("gaussSeidelVectorFinished") ;
      }

      // Initialize the iteration termination flag.
      void compute(const sequence &seq) { *gaussSeidelVectorFinished=false ; }
  } ;

  register_rule<CheckGaussSeidelVectorUnit> registerCheckGaussSeidelVectorUnit ;

  // Join operator for the following apply rule.
  struct chk_join {
    void operator() (bool &r,const bool &s) { r=r||s ; }
  } ;

  // Apply rule to set termination flag when the maximum number of iterations
  // has been reached.
  class CheckGaussSeidelVectorMaxIterations : public apply_rule<param<bool>,
  chk_join> {
    private:
      param<bool> gaussSeidelVectorFinished ;
      const_param<int> igsv,vMaxIterations ;
    public :
 
      // Define input and output.
      CheckGaussSeidelVectorMaxIterations() {
        name_store("gaussSeidelVectorFinished",gaussSeidelVectorFinished) ;
        name_store("$igsv",igsv) ;
        name_store("vMaxIterations",vMaxIterations) ;
        input("$igsv,vMaxIterations") ;
        input("gaussSeidelVectorFinished") ;
        output("gaussSeidelVectorFinished") ;
        constraint("$igsv,vMaxIterations") ;
      }

      // Set the iteration termination flag.
      void compute(const sequence &seq) {
        *gaussSeidelVectorFinished=(*igsv==*vMaxIterations) ;
      }
  } ;

  register_rule<CheckGaussSeidelVectorMaxIterations>
    registerCheckGaussSeidelVectorMaxIterations ;

  // Collapse rule for the Gauss-Seidel iteration.
  class CollapseGaussSeidelVector : public pointwise_rule {
    private:
      store<vect3d> vStar ;
    public:

      // Define input and output.
      CollapseGaussSeidelVector() {
        name_store("vStar{igsv}",vStar) ;
        input("vStar{igsv}") ;
        output("vStar=vStar{igsv}") ;
        constraint("geom_cells") ;
        conditional("gaussSeidelVectorFinished{igsv}") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<CollapseGaussSeidelVector> registerCollapseGaussSeidelVector ;
}

// Gauss-Seidel rules for pPrime.
namespace streamUns {

  // Build rule for the Gauss-Seidel iteration process. Note that pressure is
  // passed, but is not used. This is only required in order to get the build
  // rule to schedule. This should not be necessary. Ed is going to fix this
  // in a future version of Loci.
  class BuildGaussSeidelPPrime : public pointwise_rule {
    private:
      const_store<real> p ;
      store<real> pPrime ;
    public:

      // Define input and output.
      BuildGaussSeidelPPrime() {
        name_store("p",p) ;
        name_store("pPrime{igsp=0}",pPrime) ;
        input("p") ;
        output("pPrime{igsp=0}") ;
        constraint("pPrime_SGSLinearSolver,geom_cells") ;
      }

      // Initialize value for a single cell.
      void calculate(Entity cell) { pPrime[cell]=0.0 ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BuildGaussSeidelPPrime>
    registerBuildGaussSeidelPPrime ;

  // Rule for forward pass of Gauss-Seidel iteration process.
  class GaussSeidelPPrimeForward : public pointwise_rule {
    private:
      const_store<real> L,U ;
      const_store<real> D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> pPrimeForward ;
      store<real> pPrime ;
    public:

      // Define input and output.
      GaussSeidelPPrimeForward() {
        name_store("pPrime_L",L) ;
        name_store("pPrime_U",U) ;
        name_store("pPrime_D",D) ;
        name_store("pPrime_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pPrime",pPrime) ;
        name_store("pPrimeForward",pPrimeForward) ;
        input("pPrime_D,pPrime_B") ;
        input("upper->pPrime_U,lower->pPrime_L") ;
        input("upper->cr->pPrime") ;
        input("lower->cl->pPrimeForward") ;
        output("pPrimeForward=pPrime") ;
        set_relaxed_recursion() ;
      }

      // Forward pass for a single cell.
      void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*pPrime[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*pPrimeForward[cl[*li]] ;

        // Compute pPrime. This gets copied into pPrimeForward later.
        pPrime[cell]=temp*(1.0/D[cell]) ;
      }

      // Forward pass for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelPPrimeForward>
    registerGaussSeidelPPrimeForward ;

  // Rule for backward pass of Gauss-Seidel iteration process.
  class GaussSeidelPPrimeBackward : public pointwise_rule {
    private:
      const_store<real> L,U,D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> pPrimeBackward ;
      store<real> pPrimeForward ;
    public:

      // Define input and output.
      GaussSeidelPPrimeBackward() {
        name_store("pPrime_L",L) ;
        name_store("pPrime_U",U) ;
        name_store("pPrime_D",D) ;
        name_store("pPrime_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pPrimeForward",pPrimeForward) ;
        name_store("pPrimeBackward",pPrimeBackward) ;
        input("pPrime_B,pPrime_D") ;
        input("upper->pPrime_U,lower->pPrime_L") ;
        input("upper->cr->pPrimeBackward") ;
        input("lower->cl->pPrimeForward") ;
        output("pPrimeBackward=pPrimeForward") ;
        set_relaxed_recursion() ;
      }

      // Backward pass for a single cell.
       void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*pPrimeBackward[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*pPrimeForward[cl[*li]] ;

        // Compute pPrimeForward. This gets copied into pPrimeBackward later.
        pPrimeForward[cell]=temp*(1.0/D[cell]) ;
      }

      // Backward pass for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelPPrimeBackward>
    registerGaussSeidelPPrimeBackward ;

  // Rule for backward pass of Gauss-Seidel iteration process for periodic
  // boundaries.
  class GaussSeidelPPrimeBackwardPeriodic : public pointwise_rule {
    private:
      const_Map cl,cr,pmap ;
      store<real> pPrimeForward ;
    public:

      // Define input and output.
      GaussSeidelPPrimeBackwardPeriodic() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pmap",pmap) ;
        name_store("pPrimeForward",pPrimeForward) ;
        input("pmap->cl->pPrimeForward") ;
        output("cr->(pPrimeBackward=pPrimeForward)") ;
        constraint("periodicFaces") ;
      }

      // For a single face.
      void calculate(Entity face) {
        pPrimeForward[cr[face]]=pPrimeForward[cl[pmap[face]]] ; }

      // Loop over faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelPPrimeBackwardPeriodic>
    registerGaussSeidelPPrimeBackwardPeriodic ;

  // Rule for advancing the Gauss-Seidel iteration.
  class AdvanceGaussSeidelPPrime : public pointwise_rule {
    private:
      store<real> pPrime ;
    public:

      // Define input and output.
      AdvanceGaussSeidelPPrime() {
        name_store("pPrimeBackward{igsp}",pPrime) ;
        input("pPrimeBackward{igsp}") ;
        output("pPrime{igsp+1}=pPrimeBackward{igsp}") ;
        constraint("pPrime_SGSLinearSolver,geom_cells") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<AdvanceGaussSeidelPPrime>
    registerAdvanceGaussSeidelPPrime ;

  // Unit rule to intialize the iteration termination flag.
  class CheckGaussSeidelPPrimeUnit : public unit_rule {
    private:
      param<bool> gaussSeidelPPrimeFinished ;
    public :
 
      // Define input and output.
      CheckGaussSeidelPPrimeUnit() {
        name_store("gaussSeidelPPrimeFinished",gaussSeidelPPrimeFinished) ;
        constraint("UNIVERSE") ;
        output("gaussSeidelPPrimeFinished") ;
      }

      // Initialize the iteration termination flag.
      void compute(const sequence &seq) { *gaussSeidelPPrimeFinished=false ; }
  } ;

  register_rule<CheckGaussSeidelPPrimeUnit>
    registerCheckGaussSeidelPPrimeUnit ;

  // Apply rule to set termination flag when the maximum number of iterations
  // has been reached.
  class CheckGaussSeidelPPrimeMaxIterations : public apply_rule<param<bool>,
  chk_join> {
    private:
      param<bool> gaussSeidelPPrimeFinished ;
      const_param<int> igsp,pMaxIterations ;
    public :
 
      // Define input and output.
      CheckGaussSeidelPPrimeMaxIterations() {
        name_store("gaussSeidelPPrimeFinished",gaussSeidelPPrimeFinished) ;
        name_store("$igsp",igsp) ;
        name_store("pPrime_maxLinearSolverIterations",pMaxIterations) ;
        input("$igsp,pPrime_maxLinearSolverIterations") ;
        input("gaussSeidelPPrimeFinished") ;
        output("gaussSeidelPPrimeFinished") ;
        constraint("$igsp,pPrime_maxLinearSolverIterations") ;
      }

      // Set the iteration termination flag.
      void compute(const sequence &seq) {
        *gaussSeidelPPrimeFinished=(*igsp==*pMaxIterations) ;
      }
  } ;

  register_rule<CheckGaussSeidelPPrimeMaxIterations>
    registerCheckGaussSeidelPPrimeMaxIterations ;

  // Collapse rule for the Gauss-Seidel iteration.
  class CollapseGaussSeidelPPrime : public pointwise_rule {
    private:
      store<real> pPrime ;
    public:

      // Define input and output.
      CollapseGaussSeidelPPrime() {
        name_store("pPrime{igsp}",pPrime) ;
        input("pPrime{igsp}") ;
        output("pPrime=pPrime{igsp}") ;
        constraint("pPrime_SGSLinearSolver,geom_cells") ;
        conditional("gaussSeidelPPrimeFinished{igsp}") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<CollapseGaussSeidelPPrime>
    registerCollapseGaussSeidelPPrime ;
}

// Gauss-Seidel rules for pPrimeStar.
namespace streamUns {

  // Build rule for the Gauss-Seidel iteration process. Note that pressure is
  // passed, but is not used. This is only required in order to get the build
  // rule to schedule. This should not be necessary. Ed is going to fix this
  // in a future version of Loci.
  class BuildGaussSeidelPPrimeStar: public pointwise_rule {
    private:
      const_store<real> p ;
      store<real> pPrimeStar ;
    public:

      // Define input and output.
      BuildGaussSeidelPPrimeStar() {
        name_store("p",p) ;
        name_store("pPrimeStar{igspp=0}",pPrimeStar) ;
        input("p") ;
        output("pPrimeStar{igspp=0}") ;
        constraint("pPrimeStar_SGSLinearSolver,geom_cells") ;
      }

      // Initialize value for a single cell.
      void calculate(Entity cell) { pPrimeStar[cell]=0.0 ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BuildGaussSeidelPPrimeStar>
    registerBuildGaussSeidelPPrimeStar ;

  // Rule for forward pass of Gauss-Seidel iteration process.
  class GaussSeidelPPrimeStarForward : public pointwise_rule {
    private:
      const_store<real> L,U ;
      const_store<real> D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> pPrimeStarForward ;
      store<real> pPrimeStar ;
    public:

      // Define input and output.
      GaussSeidelPPrimeStarForward() {
        name_store("pPrimeStar_L",L) ;
        name_store("pPrimeStar_U",U) ;
        name_store("pPrimeStar_D",D) ;
        name_store("pPrimeStar_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pPrimeStar",pPrimeStar) ;
        name_store("pPrimeStarForward",pPrimeStarForward) ;
        input("pPrimeStar_D,pPrimeStar_B") ;
        input("upper->pPrimeStar_U,lower->pPrimeStar_L") ;
        input("upper->cr->pPrimeStar") ;
        input("lower->cl->pPrimeStarForward") ;
        output("pPrimeStarForward=pPrimeStar") ;
        set_relaxed_recursion() ;
      }

      // Forward pass for a single cell.
      void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*pPrimeStar[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*pPrimeStarForward[cl[*li]] ;

        // Compute pPrimeStar. This gets copied into pPrimeStarForward later.
        pPrimeStar[cell]=temp*(1.0/D[cell]) ;
      }

      // Forward pass for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelPPrimeStarForward>
    registerGaussSeidelPPrimeStarForward ;

  // Rule for backward pass of Gauss-Seidel iteration process.
  class GaussSeidelPPrimeStarBackward : public pointwise_rule {
    private:
      const_store<real> L,U,D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> pPrimeStarBackward ;
      store<real> pPrimeStarForward ;
    public:

      // Define input and output.
      GaussSeidelPPrimeStarBackward() {
        name_store("pPrimeStar_L",L) ;
        name_store("pPrimeStar_U",U) ;
        name_store("pPrimeStar_D",D) ;
        name_store("pPrimeStar_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pPrimeStarForward",pPrimeStarForward) ;
        name_store("pPrimeStarBackward",pPrimeStarBackward) ;
        input("pPrimeStar_B,pPrimeStar_D") ;
        input("upper->pPrimeStar_U,lower->pPrimeStar_L") ;
        input("upper->cr->pPrimeStarBackward") ;
        input("lower->cl->pPrimeStarForward") ;
        output("pPrimeStarBackward=pPrimeStarForward") ;
        set_relaxed_recursion() ;
      }

      // Backward pass for a single cell.
       void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*pPrimeStarBackward[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*pPrimeStarForward[cl[*li]] ;

        // Compute pPrimeForward. This gets copied into pPrimeBackward later.
        pPrimeStarForward[cell]=temp*(1.0/D[cell]) ;
      }

      // Backward pass for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelPPrimeStarBackward>
    registerGaussSeidelPPrimeStarBackward ;

  // Rule for backward pass of Gauss-Seidel iteration process for periodic
  // boundaries.
  class GaussSeidelPPrimeStarBackwardPeriodic : public pointwise_rule {
    private:
      const_Map cl,cr,pmap ;
      store<real> pPrimeStarForward ;
    public:

      // Define input and output.
      GaussSeidelPPrimeStarBackwardPeriodic() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pmap",pmap) ;
        name_store("pPrimeStarForward",pPrimeStarForward) ;
        input("pmap->cl->pPrimeStarForward") ;
        output("cr->(pPrimeStarBackward=pPrimeStarForward)") ;
        constraint("periodicFaces") ;
      }

      // For a single face.
      void calculate(Entity face) {
        pPrimeStarForward[cr[face]]=pPrimeStarForward[cl[pmap[face]]] ; }

      // Loop over faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelPPrimeStarBackwardPeriodic>
    registerGaussSeidelPPrimeStarBackwardPeriodic ;

  // Rule for advancing the Gauss-Seidel iteration.
  class AdvanceGaussSeidelPPrimeStar : public pointwise_rule {
    private:
      store<real> pPrimeStar ;
    public:

      // Define input and output.
      AdvanceGaussSeidelPPrimeStar() {
        name_store("pPrimeStarBackward{igspp}",pPrimeStar) ;
        input("pPrimeStarBackward{igspp}") ;
        output("pPrimeStar{igspp+1}=pPrimeStarBackward{igspp}") ;
        constraint("pPrimeStar_SGSLinearSolver,geom_cells") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<AdvanceGaussSeidelPPrimeStar>
    registerAdvanceGaussSeidelPPrimeStar ;

  // Unit rule to intialize the iteration termination flag.
  class CheckGaussSeidelPPrimeStarUnit : public unit_rule {
    private:
      param<bool> gaussSeidelPPrimeStarFinished ;
    public :
 
      // Define input and output.
      CheckGaussSeidelPPrimeStarUnit() {
        name_store("gaussSeidelPPrimeStarFinished",
          gaussSeidelPPrimeStarFinished) ;
        constraint("UNIVERSE") ;
        output("gaussSeidelPPrimeStarFinished") ;
      }

      // Initialize the iteration termination flag.
      void compute(const sequence &seq) {
        *gaussSeidelPPrimeStarFinished=false ;
      }
  } ;

  register_rule<CheckGaussSeidelPPrimeStarUnit>
    registerCheckGaussSeidelPPrimeStarUnit ;

  // Apply rule to set termination flag when the maximum number of iterations
  // has been reached.
  class CheckGaussSeidelPPrimeStarMaxIterations : public apply_rule<param<bool>,
  chk_join> {
    private:
      param<bool> gaussSeidelPPrimeStarFinished ;
      const_param<int> igspp,pMaxIterations ;
    public :
 
      // Define input and output.
      CheckGaussSeidelPPrimeStarMaxIterations() {
        name_store("gaussSeidelPPrimeStarFinished",
          gaussSeidelPPrimeStarFinished) ;
        name_store("$igspp",igspp) ;
        name_store("pPrimeStar_maxLinearSolverIterations",pMaxIterations) ;
        input("$igspp,pPrimeStar_maxLinearSolverIterations") ;
        input("gaussSeidelPPrimeStarFinished") ;
        output("gaussSeidelPPrimeStarFinished") ;
        constraint("$igspp,pPrimeStar_maxLinearSolverIterations") ;
      }

      // Set the iteration termination flag.
      void compute(const sequence &seq) {
        *gaussSeidelPPrimeStarFinished=(*igspp==*pMaxIterations) ;
      }
  } ;

  register_rule<CheckGaussSeidelPPrimeStarMaxIterations>
    registerCheckGaussSeidelPPrimeStarMaxIterations ;

  // Collapse rule for the Gauss-Seidel iteration.
  class CollapseGaussSeidelPPrimeStar : public pointwise_rule {
    private:
      store<real> pPrimeStar ;
    public:

      // Define input and output.
      CollapseGaussSeidelPPrimeStar() {
        name_store("pPrimeStar{igspp}",pPrimeStar) ;
        input("pPrimeStar{igspp}") ;
        output("pPrimeStar=pPrimeStar{igspp}") ;
        constraint("pPrimeStar_SGSLinearSolver,geom_cells") ;
        conditional("gaussSeidelPPrimeStarFinished{igspp}") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<CollapseGaussSeidelPPrimeStar>
    registerCollapseGaussSeidelPPrimeStar ;
}

// Gauss-Seidel rules for total enthalpy.
namespace streamUns {

  // Build rule for the Gauss-Seidel iteration process.
  class BuildGaussSeidelTotalEnthalpy : public pointwise_rule {
    private:
      const_store<real> h ;
      store<real> hStar ;
    public:

      // Define input and output.
      BuildGaussSeidelTotalEnthalpy() {
        name_store("h",h) ;
        name_store("hStar{igsh=0}",hStar) ;
        input("h") ;
        output("hStar{igsh=0}") ;
        constraint("compressibleFlow,geom_cells") ;
        constraint("hStar_SGSLinearSolver") ;
      }

      // Initialize value for a single cell.
      void calculate(Entity cell) { hStar[cell]=h[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BuildGaussSeidelTotalEnthalpy>
    registerBuildGaussSeidelTotalEnthalpy ;

  // Rule for forward pass of Gauss-Seidel iteration process.
  class GaussSeidelTotalEnthalpyForward : public pointwise_rule {
    private:
      const_store<real> L,U ;
      const_store<real> D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> hStarForward ;
      store<real> hStar ;
    public:

      // Define input and output.
      GaussSeidelTotalEnthalpyForward() {
        name_store("hStar_L(massFluxCorrected_p)",L) ;
        name_store("hStar_U(massFluxCorrected_p)",U) ;
        name_store("hStar_D",D) ;
        name_store("hStar_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("hStar",hStar) ;
        name_store("hStarForward",hStarForward) ;
        input("hStar_D,hStar_B") ;
        input("upper->hStar_U(massFluxCorrected_p)") ;
        input("lower->hStar_L(massFluxCorrected_p)") ;
        input("upper->cr->hStar") ;
        input("lower->cl->hStarForward") ;
        output("hStarForward=hStar") ;
        set_relaxed_recursion() ;
      }

      // Forward pass for a single cell.
      void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*hStar[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*hStarForward[cl[*li]] ;

        // Compute hStar. This gets copied into hStarForward later.
        hStar[cell]=temp*(1.0/D[cell]) ;
      }

      // Forward pass for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelTotalEnthalpyForward>
    registerGaussSeidelTotalEnthalpyForward ;

  // Rule for backward pass of Gauss-Seidel iteration process.
  class GaussSeidelTotalEnthalpyBackward : public pointwise_rule {
    private:
      const_store<real> L,U,D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> hStarBackward ;
      store<real> hStarForward ;
    public:

      // Define input and output.
      GaussSeidelTotalEnthalpyBackward() {
        name_store("hStar_L(massFluxCorrected_p)",L) ;
        name_store("hStar_U(massFluxCorrected_p)",U) ;
        name_store("hStar_D",D) ;
        name_store("hStar_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("hStarForward",hStarForward) ;
        name_store("hStarBackward",hStarBackward) ;
        input("hStar_B,hStar_D") ;
        input("upper->hStar_U(massFluxCorrected_p)") ;
        input("lower->hStar_L(massFluxCorrected_p)") ;
        input("upper->cr->hStarBackward") ;
        input("lower->cl->hStarForward") ;
        output("hStarBackward=hStarForward") ;
        set_relaxed_recursion() ;
      }

      // Backward pass for a single cell.
      void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*hStarBackward[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*hStarForward[cl[*li]] ;

        // Compute hStarForward. This gets copied into hStarBackward later.
        hStarForward[cell]=temp*(1.0/D[cell]) ;
      }

      // Backward pass for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelTotalEnthalpyBackward>
    registerGaussSeidelTotalEnthalpyBackward ;

  // Rule for backward pass of Gauss-Seidel iteration process for periodic
  // boundaries.
  class GaussSeidelTotalEnthalpyBackwardPeriodic : public pointwise_rule {
    private:
      const_Map cl,cr,pmap ;
      store<real> hStarForward ;
    public:

      // Define input and output.
      GaussSeidelTotalEnthalpyBackwardPeriodic() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pmap",pmap) ;
        name_store("hStarForward",hStarForward) ;
        input("pmap->cl->hStarForward") ;
        output("cr->(hStarBackward=hStarForward)") ;
        constraint("periodicFaces") ;
      }

      // For a single face.
      void calculate(Entity face) {
        hStarForward[cr[face]]=hStarForward[cl[pmap[face]]] ; }

      // Loop over faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelTotalEnthalpyBackwardPeriodic>
    registerGaussSeidelTotalEnthalpyBackwardPeriodic ;

  // Rule for advancing the Gauss-Seidel iteration.
  class AdvanceGaussSeidelTotalEnthalpy : public pointwise_rule {
    private:
      store<real> hStar ;
    public:

      // Define input and output.
      AdvanceGaussSeidelTotalEnthalpy() {
        name_store("hStarBackward{igsh}",hStar) ;
        input("hStarBackward{igsh}") ;
        output("hStar{igsh+1}=hStarBackward{igsh}") ;
        constraint("hStar_SGSLinearSolver,geom_cells") ;
      }

      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<AdvanceGaussSeidelTotalEnthalpy>
    registerAdvanceGaussSeidelTotalEnthalpy ;

  // Unit rule to intialize the iteration termination flag.
  class CheckGaussSeidelTotalEnthalpyUnit : public unit_rule {
    private:
      param<bool> gaussSeidelTotalEnthalpyFinished ;
    public :
 
      // Define input and output.
      CheckGaussSeidelTotalEnthalpyUnit() {
        name_store("gaussSeidelTotalEnthalpyFinished",
          gaussSeidelTotalEnthalpyFinished) ;
        constraint("UNIVERSE") ;
        output("gaussSeidelTotalEnthalpyFinished") ;
      }

      // Initialize the iteration termination flag.
      void compute(const sequence &seq) {
        *gaussSeidelTotalEnthalpyFinished=false ;
      }
  } ;

  register_rule<CheckGaussSeidelTotalEnthalpyUnit>
    registerCheckGaussSeidelTotalEnthalpyUnit ;

  // Apply rule to set termination flag when the maximum number of iterations
  // has been reached.
  class CheckGaussSeidelTotalEnthalpyMaxIterations : public
  apply_rule<param<bool>,chk_join> {
    private:
      param<bool> gaussSeidelTotalEnthalpyFinished ;
      const_param<int> igsh,hMaxIterations ;
    public :
 
      // Define input and output.
      CheckGaussSeidelTotalEnthalpyMaxIterations()
      {
        name_store("gaussSeidelTotalEnthalpyFinished",
          gaussSeidelTotalEnthalpyFinished) ;
        name_store("$igsh",igsh) ;
        name_store("hStar_maxLinearSolverIterations",hMaxIterations) ;
        input("$igsh,hStar_maxLinearSolverIterations") ;
        input("gaussSeidelTotalEnthalpyFinished") ;
        output("gaussSeidelTotalEnthalpyFinished") ;
        constraint("$igsh,hStar_maxLinearSolverIterations") ;
      }

      // Set the iteration termination flag.
      void compute(const sequence &seq) {
        *gaussSeidelTotalEnthalpyFinished=(*igsh==*hMaxIterations) ;
      }
  } ;

  register_rule<CheckGaussSeidelTotalEnthalpyMaxIterations>
    registerCheckGaussSeidelTotalEnthalpyMaxIterations ;

  // Collapse rule for the Gauss-Seidel iteration.
  class CollapseGaussSeidelTotalEnthalpy : public pointwise_rule {
    private:
      store<real> hStar ;
    public:

      // Define input and output.
      CollapseGaussSeidelTotalEnthalpy() {
        name_store("hStar{igsh}",hStar) ;
        input("hStar{igsh}") ;
        output("hStar=hStar{igsh}") ;
        constraint("hStar_SGSLinearSolver,geom_cells") ;
        conditional("gaussSeidelTotalEnthalpyFinished{igsh}") ;
      }

      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<CollapseGaussSeidelTotalEnthalpy>
    registerCollapseGaussSeidelTotalEnthalpy ;
}

/* Gauss-Seidel rules for total enthalpy in corrector stage.
namespace streamUns {

  // Build rule for the Gauss-Seidel iteration process.
  class BuildGaussSeidelTotalEnthalpyCorrector : public pointwise_rule {
    private:
      const_store<real> h ;
      store<real> hCorrected ;
    public:

      // Define input and output.
      BuildGaussSeidelTotalEnthalpyCorrector() {
        name_store("h",h) ;
        name_store("hCorrected{igshc=0}",hCorrected) ;
        input("h") ;
        output("hCorrected{igshc=0}") ;
        constraint("compressibleFlow,geom_cells") ;
        constraint("hStar_SGSLinearSolver") ;
      }

      // Initialize value for a single cell.
      void calculate(Entity cell) { hCorrected[cell]=h[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BuildGaussSeidelTotalEnthalpyCorrector>
    registerBuildGaussSeidelTotalEnthalpyCorrector ;

  // Rule for forward pass of Gauss-Seidel iteration process.
  class GaussSeidelTotalEnthalpyForwardCorrector : public pointwise_rule {
    private:
      const_store<real> L,U ;
      const_store<real> D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> hCorrectedForward ;
      store<real> hCorrected ;
    public:

      // Define input and output.
      GaussSeidelTotalEnthalpyForwardCorrector() {
        name_store("hStar_L(massFluxCorrected_c)",L) ;
        name_store("hStar_U(massFluxCorrected_c)",U) ;
        name_store("hCorrected_D",D) ;
        name_store("hCorrected_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("hCorrected",hCorrected) ;
        name_store("hCorrectedForward",hCorrectedForward) ;
        input("hCorrected_D,hCorrected_B") ;
        input("upper->hStar_U(massFluxCorrected_c)") ;
        input("lower->hStar_L(massFluxCorrected_c)") ;
        input("upper->cr->hCorrected") ;
        input("lower->cl->hCorrectedForward") ;
        output("hCorrectedForward=hCorrected") ;
        set_relaxed_recursion() ;
      }

      // Forward pass for a single cell.
      void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*hCorrected[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*hCorrectedForward[cl[*li]] ;

        // Compute hCorrected. This gets copied into hCorrectedForward later.
        hCorrected[cell]=temp*(1.0/D[cell]) ;
      }

      // Forward pass for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelTotalEnthalpyForwardCorrector>
    registerGaussSeidelTotalEnthalpyForwardCorrector ;

  // Rule for backward pass of Gauss-Seidel iteration process.
  class GaussSeidelTotalEnthalpyBackwardCorrector : public pointwise_rule {
    private:
      const_store<real> L,U,D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> hCorrectedBackward ;
      store<real> hCorrectedForward ;
    public:

      // Define input and output.
      GaussSeidelTotalEnthalpyBackwardCorrector() {
        name_store("hStar_L(massFluxCorrected_c)",L) ;
        name_store("hStar_U(massFluxCorrected_c)",U) ;
        name_store("hCorrected_D",D) ;
        name_store("hCorrected_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("hCorrectedForward",hCorrectedForward) ;
        name_store("hCorrectedBackward",hCorrectedBackward) ;
        input("hCorrected_B,hCorrected_D") ;
        input("upper->hStar_U(massFluxCorrected_c)") ;
        input("lower->hStar_L(massFluxCorrected_c)") ;
        input("upper->cr->hCorrectedBackward") ;
        input("lower->cl->hCorrectedForward") ;
        output("hCorrectedBackward=hCorrectedForward") ;
        set_relaxed_recursion() ;
      }

      // Backward pass for a single cell.
      void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*hCorrectedBackward[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*hCorrectedForward[cl[*li]] ;

        // Compute hStarForward. This gets copied into hCorrectedBackward later.
        hCorrectedForward[cell]=temp*(1.0/D[cell]) ;
      }

      // Backward pass for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelTotalEnthalpyBackwardCorrector>
    registerGaussSeidelTotalEnthalpyBackwardCorrector ;

  // Rule for backward pass of Gauss-Seidel iteration process for periodic
  // boundaries.
  class GaussSeidelTotalEnthalpyBackwardPeriodicCorrector : public pointwise_rule {
    private:
      const_Map cl,cr,pmap ;
      store<real> hCorrectedForward ;
    public:

      // Define input and output.
      GaussSeidelTotalEnthalpyBackwardPeriodicCorrector() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pmap",pmap) ;
        name_store("hCorrectedForward",hCorrectedForward) ;
        input("pmap->cl->hCorrectedForward") ;
        output("cr->(hCorrectedBackward=hCorrectedForward)") ;
        constraint("periodicFaces") ;
      }

      // For a single face.
      void calculate(Entity face) {
        hCorrectedForward[cr[face]]=hCorrectedForward[cl[pmap[face]]] ; }

      // Loop over faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelTotalEnthalpyBackwardPeriodicCorrector>
    registerGaussSeidelTotalEnthalpyBackwardPeriodicCorrector ;

  // Rule for advancing the Gauss-Seidel iteration.
  class AdvanceGaussSeidelTotalEnthalpyCorrector : public pointwise_rule {
    private:
      store<real> hCorrected ;
    public:

      // Define input and output.
      AdvanceGaussSeidelTotalEnthalpyCorrector() {
        name_store("hCorrectedBackward{igshc}",hCorrected) ;
        input("hCorrectedBackward{igshc}") ;
        output("hCorrected{igshc+1}=hCorrectedBackward{igshc}") ;
        constraint("hStar_SGSLinearSolver,geom_cells") ;
      }

      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<AdvanceGaussSeidelTotalEnthalpyCorrector>
    registerAdvanceGaussSeidelTotalEnthalpyCorrector ;

  // Unit rule to intialize the iteration termination flag.
  class CheckGaussSeidelTotalEnthalpyUnitCorrector : public unit_rule {
    private:
      param<bool> gaussSeidelTotalEnthalpyFinished ;
    public :
 
      // Define input and output.
      CheckGaussSeidelTotalEnthalpyUnitCorrector() {
        name_store("gaussSeidelTotalEnthalpyFinished_c",
          gaussSeidelTotalEnthalpyFinished) ;
        constraint("UNIVERSE") ;
        output("gaussSeidelTotalEnthalpyFinished_c") ;
      }

      // Initialize the iteration termination flag.
      void compute(const sequence &seq) {
        *gaussSeidelTotalEnthalpyFinished=false ;
      }
  } ;

  register_rule<CheckGaussSeidelTotalEnthalpyUnitCorrector>
    registerCheckGaussSeidelTotalEnthalpyUnitCorrector ;

  // Apply rule to set termination flag when the maximum number of iterations
  // has been reached.
  class CheckGaussSeidelTotalEnthalpyMaxIterationsCorrector : public
  apply_rule<param<bool>,chk_join> {
    private:
      param<bool> gaussSeidelTotalEnthalpyFinished ;
      const_param<int> igshc,hMaxIterations ;
    public :
 
      // Define input and output.
      CheckGaussSeidelTotalEnthalpyMaxIterationsCorrector()
      {
        name_store("gaussSeidelTotalEnthalpyFinished_c",
          gaussSeidelTotalEnthalpyFinished) ;
        name_store("$igshc",igshc) ;
        name_store("hStar_maxLinearSolverIterations",hMaxIterations) ;
        input("$igshc,hStar_maxLinearSolverIterations") ;
        input("gaussSeidelTotalEnthalpyFinished_c") ;
        output("gaussSeidelTotalEnthalpyFinished_c") ;
        constraint("$igshc,hStar_maxLinearSolverIterations") ;
      }

      // Set the iteration termination flag.
      void compute(const sequence &seq) {
        *gaussSeidelTotalEnthalpyFinished=(*igshc==*hMaxIterations) ;
      }
  } ;

  register_rule<CheckGaussSeidelTotalEnthalpyMaxIterationsCorrector>
    registerCheckGaussSeidelTotalEnthalpyMaxIterationsCorrector ;

  // Collapse rule for the Gauss-Seidel iteration.
  class CollapseGaussSeidelTotalEnthalpyCorrector : public pointwise_rule {
    private:
      store<real> hCorrected ;
    public:

      // Define input and output.
      CollapseGaussSeidelTotalEnthalpyCorrector() {
        name_store("hCorrected{igshc}",hCorrected) ;
        input("hCorrected{igshc}") ;
        output("hCorrected=hCorrected{igshc}") ;
        constraint("hStar_SGSLinearSolver,geom_cells") ;
        conditional("gaussSeidelTotalEnthalpyFinished_c{igshc}") ;
      }

      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<CollapseGaussSeidelTotalEnthalpyCorrector>
    registerCollapseGaussSeidelTotalEnthalpyCorrector ;
}*/

// Gauss-Seidel rules for k.
namespace streamUns {

  // Build rule for the Gauss-Seidel iteration process.
  class BuildGaussSeidelK : public pointwise_rule {
    private:
      const_store<real> k ;
      store<real> kStar ;
    public:

      // Define input and output.
      BuildGaussSeidelK() {
        name_store("k",k) ;
        name_store("kStar{igsk=0}",kStar) ;
        input("k") ;
        output("kStar{igsk=0}") ;
        constraint("geom_cells") ;
      }

      // Initialize value for a single cell.
      void calculate(Entity cell) { kStar[cell]=k[cell] ; }

      // Initialize values for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BuildGaussSeidelK> registerBuildGaussSeidelK ;

  // Rule for forward pass of Gauss-Seidel iteration process.
  class GaussSeidelKForward : public pointwise_rule {
    private:
      const_store<real> L,U ;
      const_store<real> D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> kStarForward ;
      store<real> kStar ;
    public:

      // Define input and output.
      GaussSeidelKForward() {
        name_store("kStar_L(massFluxCorrected_p)",L) ;
        name_store("kStar_U(massFluxCorrected_p)",U) ;
        name_store("kStar_D",D) ;
        name_store("kStar_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("kStar",kStar) ;
        name_store("kStarForward",kStarForward) ;
        input("kStar_D,kStar_B") ;
        input("upper->kStar_U(massFluxCorrected_p)") ;
        input("lower->kStar_L(massFluxCorrected_p)") ;
        input("upper->cr->kStar") ;
        input("lower->cl->kStarForward") ;
        output("kStarForward=kStar") ;
        set_relaxed_recursion() ;
      }

      // Forward pass for a single cell.
      void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*kStar[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*kStarForward[cl[*li]] ;

        // Compute kStar. This gets copied into kStarForward later.
        kStar[cell]=temp*(1.0/D[cell]) ;
      }

      // Forward pass for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelKForward> registerGaussSeidelKForward ;

  // Rule for backward pass of Gauss-Seidel iteration process.
  class GaussSeidelKBackward : public pointwise_rule {
    private:
      const_store<real> L,U,D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> kStarBackward ;
      store<real> kStarForward ;
    public:

      // Define input and output.
      GaussSeidelKBackward() {
        name_store("kStar_L(massFluxCorrected_p)",L) ;
        name_store("kStar_U(massFluxCorrected_p)",U) ;
        name_store("kStar_D",D) ;
        name_store("kStar_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("kStarForward",kStarForward) ;
        name_store("kStarBackward",kStarBackward) ;
        input("kStar_B,kStar_D") ;
        input("upper->kStar_U(massFluxCorrected_p)") ;
        input("lower->kStar_L(massFluxCorrected_p)") ;
        input("upper->cr->kStarBackward") ;
        input("lower->cl->kStarForward") ;
        output("kStarBackward=kStarForward") ;
        set_relaxed_recursion() ;
      }

      // Backward pass for a single cell.
      void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*kStarBackward[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*kStarForward[cl[*li]] ;

        // Compute kStarForward. This gets copied into kStarBackward later.
        kStarForward[cell]=temp*(1.0/D[cell]) ;
      }

      // Backward pass for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelKBackward> registerGaussSeidelKBackward ;

  // Rule for backward pass of Gauss-Seidel iteration process for periodic
  // boundaries.
  class GaussSeidelKBackwardPeriodic : public pointwise_rule {
    private:
      const_Map cl,cr,pmap ;
      store<real> kStarForward ;
    public:

      // Define input and output.
      GaussSeidelKBackwardPeriodic() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pmap",pmap) ;
        name_store("kStarForward",kStarForward) ;
        input("pmap->cl->kStarForward") ;
        output("cr->(kStarBackward=kStarForward)") ;
        constraint("periodicFaces") ;
      }

      // For a single face.
      void calculate(Entity face) {
        kStarForward[cr[face]]=kStarForward[cl[pmap[face]]] ; }

      // Loop over faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelKBackwardPeriodic>
    registerGaussSeidelKBackwardPeriodic ;

  // Rule for advancing the Gauss-Seidel iteration.
  class AdvanceGaussSeidelK : public pointwise_rule {
    private:
      store<real> kStar ;
    public:

      // Define input and output.
      AdvanceGaussSeidelK() {
        name_store("kStarBackward{igsk}",kStar) ;
        input("kStarBackward{igsk}") ;
        output("kStar{igsk+1}=kStarBackward{igsk}") ;
        constraint("geom_cells") ;
      }

      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<AdvanceGaussSeidelK> registerAdvanceGaussSeidelK ;

  // Unit rule to intialize the iteration termination flag.
  class CheckGaussSeidelKUnit : public unit_rule {
    private:
      param<bool> gaussSeidelKFinished ;
    public :
 
      // Define input and output.
      CheckGaussSeidelKUnit() {
        name_store("gaussSeidelKFinished",gaussSeidelKFinished) ;
        constraint("UNIVERSE") ;
        output("gaussSeidelKFinished") ;
      }

      // Initialize the iteration termination flag.
      void compute(const sequence &seq) { *gaussSeidelKFinished=false ; }
  } ;

  register_rule<CheckGaussSeidelKUnit> registerCheckGaussSeidelKUnit ;

  // Apply rule to set termination flag when the maximum number of iterations
  // has been reached.
  class CheckGaussSeidelKMaxIterations : public apply_rule<param<bool>,
  chk_join> {
    private:
      param<bool> gaussSeidelKFinished ;
      const_param<int> igsk,kMaxIterations ;
    public :
 
      // Define input and output.
      CheckGaussSeidelKMaxIterations() {
        name_store("gaussSeidelKFinished",gaussSeidelKFinished) ;
        name_store("$igsk",igsk) ;
        name_store("turbulenceMaxIterations",kMaxIterations) ;
        input("$igsk,turbulenceMaxIterations") ;
        input("gaussSeidelKFinished") ;
        output("gaussSeidelKFinished") ;
        constraint("$igsk,turbulenceMaxIterations") ;
      }

      // Set the iteration termination flag.
      void compute(const sequence &seq) {
        *gaussSeidelKFinished=(*igsk==*kMaxIterations) ;
      }
  } ;

  register_rule<CheckGaussSeidelKMaxIterations>
    registerCheckGaussSeidelKMaxIterations ;

  // Collapse rule for the Gauss-Seidel iteration.
  class CollapseGaussSeidelK : public pointwise_rule {
    private:
      store<real> kStar ;
    public:

      // Define input and output.
      CollapseGaussSeidelK() {
        name_store("kStar{igsk}",kStar) ;
        input("kStar{igsk}") ;
        output("kStar=kStar{igsk}") ;
        constraint("geom_cells") ;
        conditional("gaussSeidelKFinished{igsk}") ;
      }

      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<CollapseGaussSeidelK> registerCollapseGaussSeidelK ;
}

// Gauss-Seidel rules for omega.
namespace streamUns {

  // Build rule for the Gauss-Seidel iteration process.
  class BuildGaussSeidelOmega : public pointwise_rule {
    private:
      const_store<real> omega ;
      store<real> omegaStar ;
    public:

      // Define input and output.
      BuildGaussSeidelOmega() {
        name_store("omega",omega) ;
        name_store("omegaStar{igso=0}",omegaStar) ;
        input("omega") ;
        output("omegaStar{igso=0}") ;
        constraint("geom_cells") ;
      }

      // Initialize value for a single cell.
      void calculate(Entity cell) { omegaStar[cell]=omega[cell] ; }

      // Initialize values for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BuildGaussSeidelOmega> registerBuildGaussSeidelOmega ;

  // Rule for forward pass of Gauss-Seidel iteration process.
  class GaussSeidelOmegaForward : public pointwise_rule {
    private:
      const_store<real> L,U ;
      const_store<real> D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> omegaStarForward ;
      store<real> omegaStar ;
    public:

      // Define input and output.
      GaussSeidelOmegaForward() {
        name_store("omegaStar_L(massFluxCorrected_p)",L) ;
        name_store("omegaStar_U(massFluxCorrected_p)",U) ;
        name_store("omegaStar_D",D) ;
        name_store("omegaStar_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("omegaStar",omegaStar) ;
        name_store("omegaStarForward",omegaStarForward) ;
        input("omegaStar_D,omegaStar_B") ;
        input("upper->omegaStar_U(massFluxCorrected_p)") ;
        input("lower->omegaStar_L(massFluxCorrected_p)") ;
        input("upper->cr->omegaStar") ;
        input("lower->cl->omegaStarForward") ;
        output("omegaStarForward=omegaStar") ;
        set_relaxed_recursion() ;
      }

      // Forward pass for a single cell.
      void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*omegaStar[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*omegaStarForward[cl[*li]] ;

        // Compute omegaStar. This gets copied into omegaStarForward later.
        omegaStar[cell]=temp*(1.0/D[cell]) ;
      }

      // Forward pass for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelOmegaForward> registerGaussSeidelOmegaForward ;

  // Rule for backward pass of Gauss-Seidel iteration process.
  class GaussSeidelOmegaBackward : public pointwise_rule {
    private:
      const_store<real> L,U,D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> omegaStarBackward ;
      store<real> omegaStarForward ;
    public:

      // Define input and output.
      GaussSeidelOmegaBackward() {
        name_store("omegaStar_L(massFluxCorrected_p)",L) ;
        name_store("omegaStar_U(massFluxCorrected_p)",U) ;
        name_store("omegaStar_D",D) ;
        name_store("omegaStar_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("omegaStarForward",omegaStarForward) ;
        name_store("omegaStarBackward",omegaStarBackward) ;
        input("omegaStar_B,omegaStar_D") ;
        input("upper->omegaStar_U(massFluxCorrected_p)") ;
        input("lower->omegaStar_L(massFluxCorrected_p)") ;
        input("upper->cr->omegaStarBackward") ;
        input("lower->cl->omegaStarForward") ;
        output("omegaStarBackward=omegaStarForward") ;
        set_relaxed_recursion() ;
      }

      // Backward pass for a single cell.
      void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*omegaStarBackward[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*omegaStarForward[cl[*li]] ;

        // Compute omegaStarForward. This gets copied into omegaStarBackward
        // later.
        omegaStarForward[cell]=temp*(1.0/D[cell]) ;
      }

      // Backward pass for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelOmegaBackward> registerGaussSeidelOmegaBackward ;

  // Rule for backward pass of Gauss-Seidel iteration process for periodic
  // boundaries.
  class GaussSeidelOmegaBackwardPeriodic : public pointwise_rule {
    private:
      const_Map cl,cr,pmap ;
      store<real> omegaStarForward ;
    public:

      // Define input and output.
      GaussSeidelOmegaBackwardPeriodic() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pmap",pmap) ;
        name_store("omegaStarForward",omegaStarForward) ;
        input("pmap->cl->omegaStarForward") ;
        output("cr->(omegaStarBackward=omegaStarForward)") ;
        constraint("periodicFaces") ;
      }

      // For a single face.
      void calculate(Entity face) {
        omegaStarForward[cr[face]]=omegaStarForward[cl[pmap[face]]] ; }

      // Loop over faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelOmegaBackwardPeriodic>
    registerGaussSeidelOmegaBackwardPeriodic ;

  // Rule for advancing the Gauss-Seidel iteration.
  class AdvanceGaussSeidelOmega : public pointwise_rule {
    private:
      store<real> omegaStar ;
    public:

      // Define input and output.
      AdvanceGaussSeidelOmega() {
        name_store("omegaStarBackward{igso}",omegaStar) ;
        input("omegaStarBackward{igso}") ;
        output("omegaStar{igso+1}=omegaStarBackward{igso}") ;
        constraint("geom_cells") ;
      }

      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<AdvanceGaussSeidelOmega> registerAdvanceGaussSeidelOmega ;

  // Unit rule to intialize the iteration termination flag.
  class CheckGaussSeidelOmegaUnit : public unit_rule {
    private:
      param<bool> gaussSeidelOmegaFinished ;
    public :
 
      // Define input and output.
      CheckGaussSeidelOmegaUnit() {
        name_store("gaussSeidelOmegaFinished",gaussSeidelOmegaFinished) ;
        constraint("UNIVERSE") ;
        output("gaussSeidelOmegaFinished") ;
      }

      // Initialize the iteration termination flag.
      void compute(const sequence &seq) { *gaussSeidelOmegaFinished=false ; }
  } ;

  register_rule<CheckGaussSeidelOmegaUnit> registerCheckGaussSeidelOmegaUnit ;

  // Apply rule to set termination flag when the maximum number of iterations
  // has been reached.
  class CheckGaussSeidelOmegaMaxIterations : public apply_rule<param<bool>,
  chk_join> {
    private:
      param<bool> gaussSeidelOmegaFinished ;
      const_param<int> igso,omegaMaxIterations ;
    public :
 
      // Define input and output.
      CheckGaussSeidelOmegaMaxIterations() {
        name_store("gaussSeidelOmegaFinished",gaussSeidelOmegaFinished) ;
        name_store("$igso",igso) ;
        name_store("turbulenceMaxIterations",omegaMaxIterations) ;
        input("$igso,turbulenceMaxIterations") ;
        input("gaussSeidelOmegaFinished") ;
        output("gaussSeidelOmegaFinished") ;
        constraint("$igso,turbulenceMaxIterations") ;
      }

      // Set the iteration termination flag.
      void compute(const sequence &seq) {
        *gaussSeidelOmegaFinished=(*igso==*omegaMaxIterations) ;
      }
  } ;

  register_rule<CheckGaussSeidelOmegaMaxIterations>
    registerCheckGaussSeidelOmegaMaxIterations ;

  // Collapse rule for the Gauss-Seidel iteration.
  class CollapseGaussSeidelOmega : public pointwise_rule {
    private:
      store<real> omegaStar ;
    public:

      // Define input and output.
      CollapseGaussSeidelOmega() {
        name_store("omegaStar{igso}",omegaStar) ;
        input("omegaStar{igso}") ;
        output("omegaStar=omegaStar{igso}") ;
        constraint("geom_cells") ;
        conditional("gaussSeidelOmegaFinished{igso}") ;
      }

      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<CollapseGaussSeidelOmega> registerCollapseGaussSeidelOmega ;
}

// Gauss-Seidel rules for the species.
namespace streamUns {

  // Build rule for the Gauss-Seidel iteration process.
  class BuildGaussSeidelSpecies : public pointwise_rule {
    private:
      const_store<real> yCurr ;
      store<real> yCurrStar ;
    public:

      // Define input and output.
      BuildGaussSeidelSpecies() {
        name_store("yCurr",yCurr) ;
        name_store("yCurrStar{igsy=0}",yCurrStar) ;
        input("yCurr") ;
        output("yCurrStar{igsy=0}") ;
        constraint("geom_cells") ;
      }

      // Initialize value for a single cell.
      void calculate(Entity cell) { yCurrStar[cell]=yCurr[cell] ; }

      // Initialize values for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BuildGaussSeidelSpecies> registerBuildGaussSeidelSpecies ;

  // Rule for forward pass of Gauss-Seidel iteration process.
  class GaussSeidelSpeciesForward : public pointwise_rule {
    private:
      const_store<real> L,U ;
      const_store<real> D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> yCurrStarForward ;
      store<real> yCurrStar ;
    public:

      // Define input and output.
      GaussSeidelSpeciesForward() {
        name_store("yCurrStar_L(massFluxCorrected_p)",L) ;
        name_store("yCurrStar_U(massFluxCorrected_p)",U) ;
        name_store("yCurrStar_D",D) ;
        name_store("yCurrStar_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("yCurrStar",yCurrStar) ;
        name_store("yCurrStarForward",yCurrStarForward) ;
        input("yCurrStar_D,yCurrStar_B") ;
        input("upper->yCurrStar_U(massFluxCorrected_p)") ;
        input("lower->yCurrStar_L(massFluxCorrected_p)") ;
        input("upper->cr->yCurrStar") ;
        input("lower->cl->yCurrStarForward") ;
        output("yCurrStarForward=yCurrStar") ;
        set_relaxed_recursion() ;
      }

      // Forward pass for a single cell.
      void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*yCurrStar[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*yCurrStarForward[cl[*li]] ;

        // Compute yCurrStar.
        yCurrStar[cell]=temp*(1.0/D[cell]) ;
      }

      // Forward pass for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelSpeciesForward> registerGaussSeidelSpeciesForward ;

  // Rule for backward pass of Gauss-Seidel iteration process.
  class GaussSeidelSpeciesBackward : public pointwise_rule {
    private:
      const_store<real> L,U,D ;
      const_store<real> B ;
      const_multiMap upper,lower ;
      const_Map cl,cr ;
      const_store<real> yCurrStarBackward ;
      store<real> yCurrStarForward ;
    public:

      // Define input and output.
      GaussSeidelSpeciesBackward() {
        name_store("yCurrStar_L(massFluxCorrected_p)",L) ;
        name_store("yCurrStar_U(massFluxCorrected_p)",U) ;
        name_store("yCurrStar_D",D) ;
        name_store("yCurrStar_B",B) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("yCurrStarForward",yCurrStarForward) ;
        name_store("yCurrStarBackward",yCurrStarBackward) ;
        input("yCurrStar_B,yCurrStar_D") ;
        input("upper->yCurrStar_U(massFluxCorrected_p)") ;
        input("lower->yCurrStar_L(massFluxCorrected_p)") ;
        input("upper->cr->yCurrStarBackward") ;
        input("lower->cl->yCurrStarForward") ;
        output("yCurrStarBackward=yCurrStarForward") ;
        set_relaxed_recursion() ;
      }

      // Backward pass for a single cell.
      void calculate(Entity cell) {

        // Save the right-hand side.
        real temp=B[cell] ;

        // Add the contributions from the lower-numbered cells.
        for(const int *ui = upper.begin(cell);ui!=upper.end(cell);++ui)
          temp-=U[*ui]*yCurrStarBackward[cr[*ui]] ;

        // Add the contributions from the higher-numbered cells.
        for(const int *li = lower.begin(cell);li!=lower.end(cell);++li)
          temp-=L[*li]*yCurrStarForward[cl[*li]] ;

        // Compute yCurrStarForward.
        yCurrStarForward[cell]=temp*(1.0/D[cell]) ;
      }

      // Backward pass for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelSpeciesBackward> registerGaussSeidelSpeciesBackward ;

  // Rule for backward pass of Gauss-Seidel iteration process for periodic
  // boundaries.
  class GaussSeidelSpeciesBackwardPeriodic : public pointwise_rule {
    private:
      const_Map cl,cr,pmap ;
      store<real> yCurrStarForward ;
    public:

      // Define input and output.
      GaussSeidelSpeciesBackwardPeriodic() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pmap",pmap) ;
        name_store("yCurrStarForward",yCurrStarForward) ;
        input("pmap->cl->yCurrStarForward") ;
        output("cr->(yCurrStarBackward=yCurrStarForward)") ;
        constraint("periodicFaces") ;
      }

      // For a single face.
      void calculate(Entity face) {
        yCurrStarForward[cr[face]]=yCurrStarForward[cl[pmap[face]]] ; }

      // Loop over faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GaussSeidelSpeciesBackwardPeriodic>
    registerGaussSeidelSpeciesBackwardPeriodic ;

  // Rule for advancing the Gauss-Seidel iteration.
  class AdvanceGaussSeidelSpecies : public pointwise_rule {
    private:
      store<real> yCurrStar ;
    public:

      // Define input and output.
      AdvanceGaussSeidelSpecies() {
        name_store("yCurrStarBackward{igsy}",yCurrStar) ;
        input("yCurrStarBackward{igsy}") ;
        output("yCurrStar{igsy+1}=yCurrStarBackward{igsy}") ;
        constraint("geom_cells") ;
      }

      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<AdvanceGaussSeidelSpecies> registerAdvanceGaussSeidelSpecies ;

  // Unit rule to intialize the iteration termination flag.
  class CheckGaussSeidelSpeciesUnit : public unit_rule {
    private:
      param<bool> gaussSeidelSpeciesFinished ;
    public :
 
      // Define input and output.
      CheckGaussSeidelSpeciesUnit() {
        name_store("gaussSeidelSpeciesFinished",gaussSeidelSpeciesFinished) ;
        constraint("UNIVERSE") ;
        output("gaussSeidelSpeciesFinished") ;
      }

      // Initialize the iteration termination flag.
      void compute(const sequence &seq) { *gaussSeidelSpeciesFinished=false ; }
  } ;

  register_rule<CheckGaussSeidelSpeciesUnit>
    registerCheckGaussSeidelSpeciesUnit ;

  // Apply rule to set termination flag when the maximum number of iterations
  // has been reached.
  class CheckGaussSeidelSpeciesMaxIterations : public apply_rule<param<bool>,
  chk_join> {
    private:
      param<bool> gaussSeidelSpeciesFinished ;
      const_param<int> igsy,speciesMaxIterations ;
    public :
 
      // Define input and output.
      CheckGaussSeidelSpeciesMaxIterations() {
        name_store("gaussSeidelSpeciesFinished",gaussSeidelSpeciesFinished) ;
        name_store("$igsy",igsy) ;
        name_store("speciesMaxIterations",speciesMaxIterations) ;
        input("$igsy,speciesMaxIterations") ;
        input("gaussSeidelSpeciesFinished") ;
        output("gaussSeidelSpeciesFinished") ;
        constraint("$igsy,speciesMaxIterations") ;
      }

      // Set the iteration termination flag.
      void compute(const sequence &seq) {
        *gaussSeidelSpeciesFinished=(*igsy==*speciesMaxIterations) ;
      }
  } ;

  register_rule<CheckGaussSeidelSpeciesMaxIterations>
    registerCheckGaussSeidelSpeciesMaxIterations ;

  // Collapse rule for the Gauss-Seidel iteration.
  class CollapseGaussSeidelSpecies : public pointwise_rule {
    private:
      store<real> yCurrStar ;
    public:

      // Define input and output.
      CollapseGaussSeidelSpecies() {
        name_store("yCurrStar{igsy}",yCurrStar) ;
        input("yCurrStar{igsy}") ;
        output("yCurrStar=yCurrStar{igsy}") ;
        constraint("geom_cells") ;
        conditional("gaussSeidelSpeciesFinished{igsy}") ;
      }

      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<CollapseGaussSeidelSpecies> registerCollapseGaussSeidelSpecies ;
}

