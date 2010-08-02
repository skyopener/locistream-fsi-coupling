//-----------------------------------------------------------------------------
// Description: This file contains rules for the species equations.
//
// Authors: Siddharth Thakur and Jeff Wright
//          Streamline Numerics Inc.
//          Gainesville, Florida 32609
//-----------------------------------------------------------------------------

// Standard library includes.
#include <vector>
using std::vector ;

// Loci includes.
#include <Loci.h>

// Fluid physics library includes.
#include "eos.h" 
#include "reaction.h"
using fluidPhysics::EOS ;
using fluidPhysics::reaction ;

// StreamUns includes.
#include "sciTypes.h"
#include "varsFileInputs.h"

// CVODE includes.
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <cvode/cvode_spgmr.h>
#include <cvode/cvode_spbcgs.h>
#include <cvode/cvode_sptfqmr.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

namespace streamUns {

  using Loci::Area ;

//-----------------------------------------------------------------------------
// Classes and functions used to communicate with CVODE.

  // This class is used to pass data to CVODE.
  class CVODEReactionData {
    public:
      int numSpecies ;
      const EOS &eos ;
      EOS::State eosState ;
      real rho,P,T,h ;
      const reaction &reactor ;
      reaction::rates *reactionRate ;
      tmp_array<real> rhoTemp,yTemp,w,dwdr ;
      float *hint ;
      bool invertEOS ;
    public:
      CVODEReactionData(int numSpecies,const EOS &eos,real rho,real P,real T,
      real h,const reaction &reactor,bool invertEOS) : numSpecies(numSpecies),
      eos(eos),rho(rho),P(P),T(T),h(h),reactor(reactor),rhoTemp(numSpecies),
      yTemp(numSpecies),w(numSpecies),dwdr(numSpecies),invertEOS(invertEOS) {
        reactionRate=new reaction::rates[reactor.num_rates()] ;
        hint=new float[numSpecies+2] ;
      }

      ~CVODEReactionData() {
        delete [] reactionRate ; delete [] hint ;
      }
  } ;

  // The right-hand side function called by the CVODE integrator.
  static int cvodeRHS(realtype t,N_Vector y,N_Vector ydot,void *f_data) {
                                                                                
    // Recast the reaction data.
    CVODEReactionData *reactionData=(CVODEReactionData*)f_data ;
    int n=reactionData->numSpecies ;

    // The values in the argument "y" of this function are species masses,
    // but we must convert to mass fraction to pass to the EOS inversion
    // function.
    for(int i=0;i<n;++i){
      if(reactionData->invertEOS) reactionData->rhoTemp[i]=NV_Ith_S(y,i) ;
      reactionData->yTemp[i]=NV_Ith_S(y,i)/reactionData->rho ;
    }

    if(reactionData->invertEOS){

      // Set the hint values from the current temperature and pressure.
      reactionData->hint[0]=reactionData->T ;
      reactionData->hint[1]=reactionData->P ;

      // Solve the EOS given rho, h and y to get T and P.
      reactionData->eosState=reactionData->eos.State_from_rho_h(reactionData->
        rhoTemp,reactionData->h,reactionData->hint) ;
      reactionData->P=reactionData->eosState.pressure() ;
      reactionData->T=reactionData->eosState.temperature() ;

      // Compute the reaction rates based on the new temperature and pressure.
      reactionData->reactor.extract_rates(reactionData->reactionRate,
        reactionData->eos,reactionData->T) ;
    }

    // Compute the production rates based on the new mass fractions that
    // have been passed into this function.
    reactionData->reactor.compute_w(reactionData->w,reactionData->
      reactionRate,reactionData->yTemp,reactionData->eosState,
      reactionData->eos) ;

    // Copy computed production rate to CVODE data structure.
    for(int i=0;i<n;++i) NV_Ith_S(ydot,i)=reactionData->w[i] ;

    // Return success.
    return 0 ;
  }

  /* The right-hand side function called by the CVODE integrator.
  static int cvodeJacobian(long int N,DenseMat J,realtype t,N_Vector y,N_Vector
  fy,void *jac_data,N_Vector tmp1,N_Vector tmp2,N_Vector tmp3) {
                                                                                
    // Recast the reaction data.
    CVODEReactionData *reactionData=(CVODEReactionData*)jac_data ;
    int n=reactionData->numSpecies ;

    // The values in the argument "y" of this function are species masses,
    // but we must convert to mass fraction to pass to the EOS inversion
    // function.
    for(int i=0;i<n;++i){
      reactionData->rhoTemp[i]=NV_Ith_S(y,i) ;
      reactionData->yTemp[i]=NV_Ith_S(y,i)/reactionData->rho ;
    }

    // Solve the EOS given rho, h and y to get T and P.
    reactionData->hint[0]=reactionData->T ;
    reactionData->hint[1]=reactionData->P ;
    EOS::State eosState=reactionData->eos.State_from_rho_h(reactionData->
      rhoTemp,reactionData->h,reactionData->hint) ;
    reactionData->P=eosState.pressure() ;
    reactionData->T=eosState.temperature() ;

    // Copute the reaction rates based on the new temperature and pressure.
    reactionData->reactor.extract_rates(reactionData->reactionRate,
      reactionData->eos,reactionData->T) ;

    // Compute the Jacobian entries for each species based on the new mass
    // fractions that have been passed into this function.
    for(int i=0;i<n;++i){
      reactionData->reactor.compute_dwdrs(reactionData->dwdr,i,reactionData->
        reactionRate,0.0,NULL,reactionData->yTemp,eosState,reactionData->eos) ;
      for(int j=0;j<n;++j) DENSE_ELEM(J,i,j)=reactionData->dwdr[j] ;
    }

    // Return success.
    return 0 ;
  }*/

//-----------------------------------------------------------------------------
// Rule to set operator splitting constraints.
                                                                                
  class OperatorSplittingConstraints : public constraint_rule {
    private:
      const_param<string> operatorSplitting ;
      Constraint sourceSplitting,strangSplitting ;
    public:
                                                                                
      // Define input and output.
      OperatorSplittingConstraints() {
        name_store("operatorSplitting",operatorSplitting) ;
        name_store("sourceSplitting",sourceSplitting) ;
        name_store("strangSplitting",strangSplitting) ;
        input("operatorSplitting") ;
        output("sourceSplitting,strangSplitting") ;
      }
                                                                                
      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        if(*operatorSplitting=="source"){
          sourceSplitting=~EMPTY ; strangSplitting=EMPTY ;
        }else if(*operatorSplitting=="strang"){
          sourceSplitting=EMPTY ; strangSplitting=~EMPTY ;
        }else{
          cerr << "Bad operatorSplitting type in .vars file." << endl ;
          Loci::Abort() ;
        }
      }
                                                                                
  } ;
                                                                                
  register_rule<OperatorSplittingConstraints>
    registerOperatorSplittingConstraints ;

//-----------------------------------------------------------------------------
// Rule to set species production constraints.
                                                                                
  class SpeciesProductionConstraints : public constraint_rule {
    private:
      const_param<string> speciesProduction ;
      Constraint cvode,instantaneous ;
    public:
                                                                                
      // Define input and output.
      SpeciesProductionConstraints() {
        name_store("speciesProduction",speciesProduction) ;
        name_store("speciesProductionCVODE",cvode) ;
        name_store("speciesProductionInstantaneous",instantaneous) ;
        input("speciesProduction") ;
        output("speciesProductionCVODE,speciesProductionInstantaneous") ;
      }
                                                                                
      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        if(*speciesProduction=="CVODE"){
          cvode=~EMPTY ; instantaneous=EMPTY ;
        }else if(*speciesProduction=="instantaneous"){
          cvode=EMPTY ; instantaneous=~EMPTY ;
        }else{
          cerr << "Bad speciesProduction type in .vars file." << endl ;
          Loci::Abort() ;
        }
      }
                                                                                
  } ;
                                                                                
  register_rule<SpeciesProductionConstraints>
    registerSpeciesProductionConstraints ;

//-----------------------------------------------------------------------------
// Rules to process species equation options from the .vars file.
                                                                                
  // Creates the species equation solver constraints.
  class SpeciesEquationSolverConstraints : public constraint_rule {
    private:
      const_param<SpeciesEquationOptions> speciesEquationOptions ;
      Constraint yCurrStar_SGSLinearSolver,yCurrStar_PETSCLinearSolver ;
      Constraint yCurrStar_HYPRELinearSolver ;
    public:
                                                                                
      // Define input and output.
      SpeciesEquationSolverConstraints() {
        name_store("speciesEquationOptions",speciesEquationOptions) ;
        name_store("yCurrStar_SGSLinearSolver",yCurrStar_SGSLinearSolver) ;
        name_store("yCurrStar_PETSCLinearSolver",yCurrStar_PETSCLinearSolver) ;
        name_store("yCurrStar_HYPRELinearSolver",yCurrStar_HYPRELinearSolver) ;
        input("speciesEquationOptions") ;
        output("yCurrStar_SGSLinearSolver,yCurrStar_PETSCLinearSolver") ;
        output("yCurrStar_HYPRELinearSolver") ;
      }
                                                                                
      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        if((*speciesEquationOptions).optionExists("linearSolver")){
          Loci::option_value_type optionValueType=speciesEquationOptions->
            getOptionValueType("linearSolver") ;
          switch(optionValueType){
            case Loci::NAME:
              {
                Loci::option_values optionValues=speciesEquationOptions->
                  getOption("linearSolver") ;
                string name ; optionValues.get_value(name) ;
                if(name=="SGS"){
                  yCurrStar_SGSLinearSolver=~EMPTY ;
                  yCurrStar_PETSCLinearSolver=EMPTY ;
                  yCurrStar_HYPRELinearSolver=EMPTY ;
                }else if(name=="PETSC"){
                  yCurrStar_SGSLinearSolver=EMPTY ;
                  yCurrStar_PETSCLinearSolver=~EMPTY ;
                  yCurrStar_HYPRELinearSolver=EMPTY ;
                }else if(name=="HYPRE"){
                  yCurrStar_SGSLinearSolver=EMPTY ;
                  yCurrStar_PETSCLinearSolver=EMPTY ;
                  yCurrStar_HYPRELinearSolver=~EMPTY ;
                }else{
                  cerr << "Bad linearSolver for speciesEquation." << endl ;
                  Loci::Abort() ;
                }
              }
              break ;
            default:
              cerr << "Bad type for linearSolver in speciesEquation." << endl ;
              Loci::Abort() ;
          }
        }else{
          yCurrStar_SGSLinearSolver=~EMPTY ;
        }
      }
                                                                                
  } ;
                                                                                
  register_rule<SpeciesEquationSolverConstraints>
    registerSpeciesEquationSolverConstraints ;

  // Creates the species equation solver parameters.
  class SpeciesEquationSolverParameters : public singleton_rule {
    private:
      const_param<SpeciesEquationOptions> speciesEquationOptions ;
      param<int> yCurrStar_maxLinearSolverIterations ;
      param<real> speciesRelaxationFactor ;
    public:
                                                                                
      // Define input and output.
      SpeciesEquationSolverParameters() {
        name_store("speciesEquationOptions",speciesEquationOptions) ;
        name_store("yCurrStar_maxLinearSolverIterations",
          yCurrStar_maxLinearSolverIterations) ;
        name_store("speciesRelaxationFactor",speciesRelaxationFactor) ;
        input("speciesEquationOptions") ;
        output("yCurrStar_maxLinearSolverIterations,speciesRelaxationFactor") ;
      }
                                                                                
      // Set up the parameters.
      virtual void compute(const sequence& seq) {
                                                                                
        // Relaxation factor.
        if((*speciesEquationOptions).optionExists("relaxationFactor")){
          Loci::option_value_type optionValueType=speciesEquationOptions->
            getOptionValueType("relaxationFactor") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                speciesEquationOptions->getOption("relaxationFactor",temp) ;
                if(temp<=0.0 || temp>1.0){
                  cerr << "Bad relaxationFactor for speciesEquation." << endl ;
                  Loci::Abort() ;
                }
                *speciesRelaxationFactor=temp ;
              }
              break ;
            default:
              cerr << "Bad type for relaxationFactor in speciesEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *speciesRelaxationFactor=0.5 ;
        }
                                                                                
        // Maximum number of iterations.
        if((*speciesEquationOptions).optionExists("maxIterations")){
          Loci::option_value_type optionValueType=speciesEquationOptions->
            getOptionValueType("maxIterations") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                speciesEquationOptions->getOption("maxIterations",temp) ;
                if(int(temp)<0){
                  cerr << "Bad maxIterations value for speciesEquation."
                    << endl ; Loci::Abort() ;
                }
                *yCurrStar_maxLinearSolverIterations=int(temp) ;
              }
              break ;
            default:
              cerr << "Bad type for maxIterations in speciesEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *yCurrStar_maxLinearSolverIterations=5 ;
        }
      }
  } ;
                                                                                
  register_rule<SpeciesEquationSolverParameters>
    registerSpeciesEquationSolverParamters ;

//-----------------------------------------------------------------------------
// Rules to provide species for flows with no species transport. Simplifies the
// coding in several other places.

  // Initial species.
  class InitialSpeciesInterior : public pointwise_rule {
    private:
      storeVec<real> y_ic ;
    public:

      // Define input and output.
      InitialSpeciesInterior() {
        name_store("y_ic",y_ic) ;
        output("y_ic") ;
        constraint("noSpeciesTransport,geom_cells") ;
      }

      // Now, for a pure material we have a since species mass fraction
      // with a value of one.
      void calculate(Entity cell) { y_ic[cell][0]=1.0 ; }

      // Assign species mass fractions for all faces in sequence.
      virtual void compute(const sequence &seq) {
        y_ic.setVecSize(1) ; do_loop(seq,this) ;
      }
  } ;

  register_rule<InitialSpeciesInterior> registerInitialSpeciesInterior ;

  // Species for cells.
  class SpeciesInterior : public pointwise_rule {
    private:
      storeVec<real> y ;
    public:

      // Define input and output.
      SpeciesInterior() {
        name_store("y",y) ;
        output("y") ;
        constraint("noSpeciesTransport,geom_cells") ;
      }

      // Now, for a pure material we have a since species mass fraction
      // with a value of one.
      void calculate(Entity cell) { y[cell][0]=1.0 ; }

      // Assign species mass fractions for all faces in sequence.
      virtual void compute(const sequence &seq) {
        y.setVecSize(1) ; do_loop(seq,this) ;
      }
  } ;

  register_rule<SpeciesInterior> registerSpeciesInterior ;

  // Species for cells.
  class SpeciesBoundary : public pointwise_rule {
    private:
      storeVec<real> y_f ;
    public:

      // Define input and output.
      SpeciesBoundary() {
        name_store("y_f",y_f) ;
        output("y_f") ;
        constraint("noSpeciesTransport,boundaryFaces") ;
      }

      // Now, for a pure material we have a since species mass fraction
      // with a value of one.
      void calculate(Entity face) { y_f[face][0]=1.0 ; }

      // Assign species mass fractions for all faces in sequence.
      virtual void compute(const sequence &seq) {
        y_f.setVecSize(1) ; do_loop(seq,this) ;
      }
  } ;

  register_rule<SpeciesBoundary> registerSpeciesBoundary ;

  // Starred species for use in temperature inversion from enthalpy.
  class StarredSpeciesInterior : public pointwise_rule {
    private:
      storeVec<real> y ;
    public:

      // Define input and output.
      StarredSpeciesInterior() {
        name_store("yStar",y) ;
        output("yStar") ;
        constraint("noSpeciesTransport,geom_cells") ;
      }

      // Now, for a pure material we have a since species mass fraction
      // with a value of one.
      void calculate(Entity cell) { y[cell][0]=1.0 ; }

      // Assign species mass fractions for all faces in sequence.
      virtual void compute(const sequence &seq) {
        y.setVecSize(1) ; do_loop(seq,this) ;
      }
  } ;

  register_rule<StarredSpeciesInterior> registerStarredSpeciesInterior ;

//-----------------------------------------------------------------------------
// Species boundary condition rules.

  // Rule for boundary faces with specified species mass fractions. Assigns
  // species mass fraction values to all boundary faces that have the property
  // y_BC. Checked.
  class BoundarySpeciesMassFractionSpecification : public pointwise_rule {
    private:
      const_param<int> numSpecies ;
      const_Map ref ;
      const_storeVec<real> y_BC ;
      storeVec<real> y_f ;
    public:

      // Define input and output.
      BoundarySpeciesMassFractionSpecification() {
        name_store("numSpecies",numSpecies) ;
        name_store("ref",ref) ;
        name_store("y_BC",y_BC) ;
        name_store("y_f",y_f) ;
        input("numSpecies,ref->y_BC") ;
        output("y_f") ;
      }

      // Assign species mass fractions for a single face.
      void calculate(Entity face) { y_f[face]=y_BC[ref[face]] ; }

      // Assign species mass fractions for all faces in sequence.
      virtual void compute(const sequence &seq) {
        y_f.setVecSize(*numSpecies) ; do_loop(seq,this) ;
      }

  } ;

  register_rule<BoundarySpeciesMassFractionSpecification>
    registerBoundarySpeciesMassFractionSpecification ;

  // Rule for extrapolating species mass fractions to boundary faces. This
  // occurs for all outlets, no-slip, slip and symmetry boundaries. Right now
  // we are using the low-order method. Checked.
  class BoundarySpeciesMassFractionExtrapolation : public pointwise_rule {
    private:
      const_param<int> numSpecies ;
      const_Map ci ;
      const_storeVec<real> y ;
      storeVec<real> y_f ;
    public:

      // Define input and output.
      BoundarySpeciesMassFractionExtrapolation() {
        name_store("numSpecies",numSpecies) ;
        name_store("ci",ci) ;
        name_store("y",y) ;
        name_store("y_f",y_f) ;
        input("numSpecies,ci->y") ;
        output("y_f") ;
        constraint("speciesTransport,extrapolatedSpeciesMassFraction_BC") ;
      }

      // Calculate species mass fractions for a single face.
      void calculate(Entity face) { y_f[face]=y[ci[face]] ; }

      // Calculate species mass fractions for all faces in sequence.
      virtual void compute(const sequence &seq) {
        y_f.setVecSize(*numSpecies) ; do_loop(seq,this) ;
      }
  } ;

  register_rule<BoundarySpeciesMassFractionExtrapolation>
    registerBoundarySpeciesMassFractionExtrapolation ;

//-----------------------------------------------------------------------------
// Time-scheme independent rules for the species loop. Note that each time-
// stepping scheme will have a different build rule for the species loop.

  // Rule to extract the interior values for the current species. Checked.
  class ExtractCurrentSpeciesInterior : public pointwise_rule {
    private:
      const_param<int> is ;
      const_storeVec<real> yStar ;
      store<real> yCurr ;
    public:

      // Define input and output.
      ExtractCurrentSpeciesInterior() {
        name_store("$is{n,it,is}",is) ;
        name_store("yStar{n,it,is}",yStar) ;
        name_store("yCurr{n,it,is}",yCurr) ;
        input("$is{n,it,is},yStar{n,it,is}") ;
        output("yCurr{n,it,is}") ;
        constraint("speciesTransport{n,it,is},geom_cells{n,it,is}") ;
      }

      // Extract current species for a single cell.
      void calculate(Entity cell) { yCurr[cell]=yStar[cell][*is] ; }

      // Extract current species for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ExtractCurrentSpeciesInterior>
    registerExtractCurrentSpeciesInterior ;

  // Rule to extract the boundary values for the current species. Checked.
  class ExtractCurrentSpeciesBoundary : public pointwise_rule {
    private:
      const_param<int> is ;
      const_storeVec<real> y_f ;
      store<real> yCurr_f ;
    public:

      // Define input and output.
      ExtractCurrentSpeciesBoundary() {
        name_store("$is{n,it,is}",is) ;
        name_store("y_f{n,it,is}",y_f) ;
        name_store("yCurr_f{n,it,is}",yCurr_f) ;
        input("$is{n,it,is},y_f{n,it,is}") ;
        output("yCurr_f{n,it,is}") ;
        constraint("speciesTransport{n,it,is},boundaryFaces{n,it,is}") ;
      }

      // Extract current species for a single face.
      void calculate(Entity face) { yCurr_f[face]=y_f[face][*is] ; }

      // Extract current species for a sequence of faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ExtractCurrentSpeciesBoundary>
    registerExtractCurrentSpeciesBoundary ;

  // Advance rule for species iteration loop. Checked.
  class SpeciesIterationLoopAdvance : public pointwise_rule {
    private:
      const_param<int> is ;
      const_param<int> numSpecies ;
      const_store<real> yCurrStar ;
      storeVec<real> yStar ;
    public:

      // Define input and output.
      SpeciesIterationLoopAdvance() {
        name_store("$is{n,it,is}",is) ;
        name_store("numSpecies{n,it,is}",numSpecies) ;
        name_store("yCurrStar{n,it,is}",yCurrStar) ;
        name_store("yStar{n,it,is}",yStar) ;
        input("$is{n,it,is},numSpecies{n,it,is},yCurrStar{n,it,is}") ;
        output("yStar{n,it,is+1}=yStar{n,it,is}") ;
        constraint("speciesTransport{n,it,is},geom_cells{n,it,is}") ;
      }

      // Advance species mass fractions for a single cell.
      void calculate(Entity cell) { yStar[cell][*is]=yCurrStar[cell] ; }

      // Advance species mass fractions for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SpeciesIterationLoopAdvance>
    registerSpeciesIterationLoopAdvance ;

  // Collapse rule for species iteration loop. Checked.
  class SpeciesIterationLoopCollapse : public pointwise_rule {
    private:
      const_param<int> numSpecies ; // Added 07/07/2004.
      const_param<bool> iterationFinished ;
      storeVec<real> yStar ;
    public:

      // Define input and output.
      SpeciesIterationLoopCollapse() {
        name_store("numSpecies{n,it,is}",numSpecies) ; // Added 07/07/2004.
        name_store("speciesIterationFinished{n,it,is}",iterationFinished) ;
        name_store("yStar{n,it,is}",yStar) ;
        input("numSpecies{n,it,is}") ;
        input("speciesIterationFinished{n,it,is}") ;
        input("yStar{n,it,is}") ;
        output("yStar{n,it}=yStar{n,it,is}") ;
        conditional("speciesIterationFinished{n,it,is}") ;
        constraint("speciesTransport{n,it,is},geom_cells{n,it,is}") ;
      }

      // Re-normalize species mass fractions for a cell. In order to prevent
      // floating-point exceptions in the divide (which occurs on the Pentiums
      // when yStar[cell][i]<4.94066e-324 is divided by ySum>1.0) we must
      // limit the species to a reasonable range. Added 07/07/2004.
      void calculate(Entity cell) {
        real ySum=0.0 ;
        for(int i=0;i<*numSpecies;++i){
          if(yStar[cell][i]<1.0e-30) yStar[cell][i]=0.0 ;
          else if(yStar[cell][i]>1.0) yStar[cell][i]=1.0 ;
          ySum+=yStar[cell][i] ;
        }
        for(int i=0;i<*numSpecies;++i) yStar[cell][i]/=ySum ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SpeciesIterationLoopCollapse>
    registerSpeciesIterationLoopCollapse ;

  // Class to determine it the species iteration is finished.
  class CheckSpeciesIterationFinished : public singleton_rule {
    private:
      const_param<int> numSpecies,is ;
      param<bool> iterationFinished ;
    public:

      // Define input and output.
      CheckSpeciesIterationFinished() {
        name_store("numSpecies",numSpecies) ;
        name_store("$is",is) ;
        name_store("speciesIterationFinished",iterationFinished) ;
        input("numSpecies,$is") ;
        output("speciesIterationFinished") ;
      }

      // Check if iteration is finished.
      void compute(const sequence &seq) {
        *iterationFinished=(*is==*numSpecies) ;
      }
  } ;

  register_rule<CheckSpeciesIterationFinished>
    registerCheckSpeciesIterationFinished ;

//-----------------------------------------------------------------------------
// Rules to create a constraint for boundary faces with non-zero species 
// diffusion flux.

  // All inlet faces have non-zero diffusion flux.
  class BoundarySpeciesDiffusionInlet : public pointwise_rule {
    private:
      store<bool> boundarySpeciesDiffusion ;
    public:

      // Define input and output.
      BoundarySpeciesDiffusionInlet() {
        name_store("boundarySpeciesDiffusion",boundarySpeciesDiffusion) ;
        output("boundarySpeciesDiffusion") ;
        constraint("inlet_BC") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { boundarySpeciesDiffusion[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundarySpeciesDiffusionInlet>
    registerBoundarySpeciesDiffusionInlet ;

  // All outlet faces have non-zero diffusion flux.
  class BoundarySpeciesDiffusionOutlet : public pointwise_rule {
    private:
      store<bool> boundarySpeciesDiffusion ;
    public:

      // Define input and output.
      BoundarySpeciesDiffusionOutlet() {
        name_store("boundarySpeciesDiffusion",boundarySpeciesDiffusion) ;
        output("boundarySpeciesDiffusion") ;
        constraint("outlet_BC") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { boundarySpeciesDiffusion[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundarySpeciesDiffusionOutlet>
    registerBoundarySpeciesDiffusionOutlet ;

  // Interface faces have non-zero diffusion flux.
  class BoundarySpeciesDiffusionInterface : public pointwise_rule {
    private:
      store<bool> boundarySpeciesDiffusion ;
    public:

      // Define input and output.
      BoundarySpeciesDiffusionInterface() {
        name_store("boundarySpeciesDiffusion",boundarySpeciesDiffusion) ;
        output("boundarySpeciesDiffusion") ;
        constraint("interface_BC") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { boundarySpeciesDiffusion[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundarySpeciesDiffusionInterface>
    registerBoundarySpeciesDiffusionInterface ;

//-----------------------------------------------------------------------------
// Scheme independent rules for assembling the species equations.

  // Species viscosity for the interior faces.
  class YViscosityInterior : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_param<real> laminarSchmidtNumber,turbulentSchmidtNumber ;
      const_store<real> laminarViscosity,eddyViscosity ;
      store<real> yViscosity ;
    public:

      // Define input and output.
      YViscosityInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("laminarSchmidtNumber",laminarSchmidtNumber) ;
        name_store("turbulentSchmidtNumber",turbulentSchmidtNumber) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("yViscosity",yViscosity) ;
        input("laminarSchmidtNumber,turbulentSchmidtNumber") ;
        input("(cl,cr)->(laminarViscosity,eddyViscosity)") ;
        output("yViscosity") ;
        constraint("internalFaces") ;
      }

      // Compute for each face.
      void calculate(Entity face) {
        yViscosity[face]=0.5*((laminarViscosity[cl[face]]+
          laminarViscosity[cr[face]])/(*laminarSchmidtNumber)+(eddyViscosity
          [cl[face]]+eddyViscosity[cr[face]])/(*turbulentSchmidtNumber)) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<YViscosityInterior> registerYViscosityInterior ;

  // Species viscosity for the boundary faces.
  class YViscosityBoundary : public pointwise_rule {
    private:
      const_param<real> laminarSchmidtNumber,turbulentSchmidtNumber ;
      const_store<real> laminarViscosity_f,eddyViscosity_f ;
      store<real> yViscosity ;
    public:

      // Define input and output.
      YViscosityBoundary() {
        name_store("laminarSchmidtNumber",laminarSchmidtNumber) ;
        name_store("turbulentSchmidtNumber",turbulentSchmidtNumber) ;
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        name_store("eddyViscosity_f",eddyViscosity_f) ;
        name_store("yViscosity",yViscosity) ;
        input("laminarSchmidtNumber,turbulentSchmidtNumber") ;
        input("laminarViscosity_f,eddyViscosity_f") ;
        output("yViscosity") ;
        constraint("boundaryFaces") ;
      }

      // Compute for each face.
      void calculate(Entity face) {
        yViscosity[face]=laminarViscosity_f[face]/(*laminarSchmidtNumber)+
          eddyViscosity_f[face]/(*turbulentSchmidtNumber) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<YViscosityBoundary> registerYViscosityBoundary ;

  // Rule to initialize the main coefficient. Checked.
  class InitializeSpeciesMainCoefficient : public unit_rule {
    private:
      store<real> yCurrMainCoefficient ;
    public:

      // Define input and output.
      InitializeSpeciesMainCoefficient() {
        name_store("yCurrMainCoefficient",yCurrMainCoefficient) ;
        output("yCurrMainCoefficient") ;
        constraint("vol") ;
      }

      // Set the main coefficient to zero for a single cell.
      void calculate(Entity cell) { yCurrMainCoefficient[cell]=0.0 ; }

      // Set the main coefficient to zero for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeSpeciesMainCoefficient>
    registerInitializeSpeciesMainCoefficient ;

  // Rule to add the first-order inviscid flux contribution to the main
  // coefficient for interior faces. Checked.
  class FOUInviscidFluxToSpeciesMainCoefficientInterior : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> massFlux ;
      store<real> yCurrMainCoefficient ;
    public:

      // Define input and output.
      FOUInviscidFluxToSpeciesMainCoefficientInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("massFlux",massFlux) ;
        name_store("yCurrMainCoefficient",yCurrMainCoefficient) ;
        input("massFlux") ;
        output("(cl,cr)->yCurrMainCoefficient") ;
        constraint("internalFaces") ;
      }

      // Increment the main coefficient for cells attached to a single face.
      void calculate(Entity face) {
        if(massFlux[face]>0.0){
          yCurrMainCoefficient[cr[face]]+=massFlux[face] ;
        }else{
          yCurrMainCoefficient[cl[face]]-=massFlux[face] ;
        }
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToSpeciesMainCoefficientInterior>
    registerFOUInviscidFluxToSpeciesMainCoefficientInterior ;

  // Rule to add the first-order inviscid flux contribution to the main
  // coefficient for boundary faces. Checked.
  class FOUInviscidFluxToSpeciesMainCoefficientBoundary : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> massFlux ;
      store<real> yCurrMainCoefficient ;
    public:

      // Define input and output.
      FOUInviscidFluxToSpeciesMainCoefficientBoundary() {
        name_store("ci",ci) ;
        name_store("massFlux",massFlux) ;
        name_store("yCurrMainCoefficient",yCurrMainCoefficient) ;
        input("massFlux") ;
        output("ci->yCurrMainCoefficient") ;
        constraint("boundaryFaces") ;
      }

      // Increment the main coefficient for the cell attach to a single face.
      void calculate(Entity face) {
        if(massFlux[face]<0.0) yCurrMainCoefficient[ci[face]]-=massFlux[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToSpeciesMainCoefficientBoundary>
    registerFOUInviscidFluxToSpeciesMainCoefficientBoundary ;

  // Rule to add the diffusive flux contribution to the main coefficient for
  // interior faces. Checked 07/10/04.
  class DiffusiveFluxToSpeciesMainCoefficientInterior : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> faceRadius ;
      const_store<real> diffusionProduct ;
      const_store<real> yViscosity ;
      store<real> yCurrMainCoefficient ;
    public:

      // Define input and output.
      DiffusiveFluxToSpeciesMainCoefficientInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("faceRadius",faceRadius) ;
        name_store("yViscosity",yViscosity) ;
        name_store("yCurrMainCoefficient",yCurrMainCoefficient) ;
        input("diffusionProduct,faceRadius,yViscosity") ;
        output("(cl,cr)->yCurrMainCoefficient") ;
        constraint("internalFaces,viscousFlow") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real temp=yViscosity[face]*diffusionProduct[face]*
          faceRadius[face] ;
        yCurrMainCoefficient[cl[face]]+=temp ;
        yCurrMainCoefficient[cr[face]]+=temp ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToSpeciesMainCoefficientInterior>
    registerDiffusiveFluxToSpeciesMainCoefficientInterior ;

  // Rule to add the diffusive flux contribution to the main coefficient for
  // boundary faces. Checked 07/10/04.
  class DiffusiveFluxToSpeciesMainCoefficientBoundary : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> diffusionProduct ;
      const_store<real> faceRadius ;
      const_store<real> yViscosity ;
      store<real> yCurrMainCoefficient ;
    public:

      // Define input and output.
      DiffusiveFluxToSpeciesMainCoefficientBoundary() {
        name_store("ci",ci) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("faceRadius",faceRadius) ;
        name_store("yViscosity",yViscosity) ;
        name_store("yCurrMainCoefficient",yCurrMainCoefficient) ;
        input("diffusionProduct,faceRadius,yViscosity") ;
        output("ci->yCurrMainCoefficient") ;
        constraint("boundarySpeciesDiffusion,viscousFlow") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        yCurrMainCoefficient[ci[face]]+=yViscosity[face]*
          diffusionProduct[face]*faceRadius[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToSpeciesMainCoefficientBoundary>
    registerDiffusiveFluxToSpeciesMainCoefficientBoundary ;

  // Rule to ensure main coefficient is not zero for inviscid flows.
  class NetMassFluxToSpeciesMainCoefficient : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_store<real> netMassFlux ;
      store<real> yCurrMainCoefficient ;
    public:
                                                                                
      // Define input and output.
      NetMassFluxToSpeciesMainCoefficient() {
        name_store("netMassFlux",netMassFlux) ;
        name_store("yCurrMainCoefficient",yCurrMainCoefficient) ;
        input("netMassFlux") ;
        output("yCurrMainCoefficient") ;
        constraint("inviscidFlow,geom_cells") ;
      }
                                                                                
      // Add net mass flux for a single cell.
      void calculate(Entity cell) {
        if(netMassFlux[cell]<0.0) yCurrMainCoefficient[cell]-=
          netMassFlux[cell] ;
      }
                                                                                
      // Add net mass flux for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<NetMassFluxToSpeciesMainCoefficient>
    registerNetMassFluxToSpeciesMainCoefficient ;

  // Rule to initialize the source term. Checked.
  class InitializeSpeciesSourceTerm : public unit_rule {
    private:
      store<real> yCurrSourceTerm ;
    public:

      // Define input and output.
      InitializeSpeciesSourceTerm() {
        name_store("yCurrSourceTerm",yCurrSourceTerm) ;
        output("yCurrSourceTerm") ;
        constraint("vol") ;
      }

      // Set the source term to zero for a single cell.
      void calculate(Entity cell) { yCurrSourceTerm[cell]=0.0 ; }

      // Set the source term to zero for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeSpeciesSourceTerm>
    registerInitializeSpeciesSourceTerm ;

  // Rule to add the first-order inviscid flux contribution to the source term
  // for boundary faces. Checked.
  class FOUInviscidFluxToSpeciesSourceTermBoundary : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> yCurr ; // 5/13/05
      const_store<real> massFlux ;
      const_store<real> yCurr_f ;
      store<real> yCurrSourceTerm ;
    public:

      // Define input and output.
      FOUInviscidFluxToSpeciesSourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("yCurr",yCurr) ;
        name_store("massFlux",massFlux) ;
        name_store("yCurr_f",yCurr_f) ;
        name_store("yCurrSourceTerm",yCurrSourceTerm) ;
        input("ci->yCurr") ; // 5/13/05
        input("massFlux,yCurr_f") ;
        output("ci->yCurrSourceTerm") ;
        constraint("boundaryFaces,fouInviscidFlux") ;
      }

      // Increment the source term for the cell attach to a single face.
      void calculate(Entity face) {
        if(massFlux[face]<=0.0) yCurrSourceTerm[ci[face]]-=massFlux[face]*
          yCurr_f[face] ;
        else yCurrSourceTerm[ci[face]]-=massFlux[face]*(yCurr_f[face]-
          yCurr[ci[face]]) ;// 5/13/05
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToSpeciesSourceTermBoundary>
    registerFOUInviscidFluxToSpeciesSourceTermBoundary ;

  // Rule to add the second-order convection contribution to the source term
  // for interior faces. Checked.
  class SOUInviscidFluxToSpeciesSourceTermInterior : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<vect3d> yCurrGradient ;
      const_store<real> yCurrLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      store<real> yCurrSourceTerm ;
    public:

      // Define input and output.
      SOUInviscidFluxToSpeciesSourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("grads(yCurr)",yCurrGradient) ;
        name_store("limiters(yCurr)",yCurrLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("massFlux",massFlux) ;
        name_store("yCurrSourceTerm",yCurrSourceTerm) ;
        input("(cl,cr)->(grads(yCurr),limiters(yCurr),cellcenter)") ;
        input("facecenter,massFlux") ;
        output("(cl,cr)->yCurrSourceTerm") ;
        constraint("internalFaces,souOrRoeInviscidFlux") ;
      }

      // Increment the source term for the cells attached to a single face.
      void calculate(Entity face) {
        real secondOrderSource=(massFlux[face]>0.0)? massFlux[face]*
          yCurrLimiter[cl[face]]*dot(yCurrGradient[cl[face]],faceCenter[face]-
          cellCenter[cl[face]]):massFlux[face]*yCurrLimiter[cr[face]]*
          dot(yCurrGradient[cr[face]],faceCenter[face]-cellCenter[cr[face]]) ;
        yCurrSourceTerm[cl[face]]-=secondOrderSource ;
        yCurrSourceTerm[cr[face]]+=secondOrderSource ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SOUInviscidFluxToSpeciesSourceTermInterior>
    registerSOUInviscidFluxToSpeciesSourceTermInterior ;

  // Rule to add the second-order convection contribution to the source term
  // for boundary faces. Checked.
  class SOUInviscidFluxToSpeciesSourceTermBoundary : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<vect3d> yCurrGradient ;
      const_store<real> yCurrLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      const_store<real> yCurr_f ;
      store<real> yCurrSourceTerm ;
    public:

      // Define input and output.
      SOUInviscidFluxToSpeciesSourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("grads(yCurr)",yCurrGradient) ;
        name_store("limiters(yCurr)",yCurrLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("massFlux",massFlux) ;
        name_store("yCurr_f",yCurr_f) ;
        name_store("yCurrSourceTerm",yCurrSourceTerm) ;
        input("ci->(grads(yCurr),limiters(yCurr),cellcenter)") ;
        input("facecenter,massFlux,yCurr_f") ;
        output("ci->yCurrSourceTerm") ;
        constraint("boundaryFaces,souOrRoeInviscidFlux") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        yCurrSourceTerm[ci[face]]-=(massFlux[face]>0.0)? massFlux[face]*
          yCurrLimiter[ci[face]]*dot(yCurrGradient[ci[face]],faceCenter[face]-
          cellCenter[ci[face]]):massFlux[face]*yCurr_f[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SOUInviscidFluxToSpeciesSourceTermBoundary>
    registerSOUInviscidFluxToSpeciesSourceTermBoundary ;

  // Rule to add the diffusive flux contribution to the source term for
  // interior faces. Checked 07/10/04.
  class DiffusiveFluxToSpeciesSourceTermInterior : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> vol ;
      const_store<vect3d> yCurrGradient ;
      const_store<real> diffusionProduct ;
      const_store<vect3d> geometryFactor0 ;
      const_store<real> faceRadius ;
      const_store<real> yViscosity ;
      store<real> yCurrSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToSpeciesSourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("vol",vol) ;
        name_store("grads(yCurr)",yCurrGradient) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("geometryFactor0",geometryFactor0) ;
        name_store("faceRadius",faceRadius) ;
        name_store("yViscosity",yViscosity) ;
        name_store("yCurrSourceTerm",yCurrSourceTerm) ;
        input("(cl,cr)->(vol,grads(yCurr))") ;
        input("diffusionProduct,geometryFactor0,faceRadius,yViscosity") ;
        output("(cl,cr)->yCurrSourceTerm") ;
        constraint("internalFaces,viscousFlow") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real secondarySourceTerm=yViscosity[face]*
          dot((yCurrGradient[cl[face]]*vol[cr[face]]+yCurrGradient[cr[face]]*
          vol[cl[face]])/(vol[cl[face]]+vol[cr[face]]),geometryFactor0[face])*
          faceRadius[face] ;
        yCurrSourceTerm[cl[face]]+=secondarySourceTerm ;
        yCurrSourceTerm[cr[face]]-=secondarySourceTerm ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToSpeciesSourceTermInterior>
    registerDiffusiveFluxToSpeciesSourceTermInterior ;

  // Rule to add the diffusive flux contribution to the source term for
  // boundary faces. Checked 07/10/04.
  class DiffusiveFluxToSpeciesSourceTermBoundary : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> yCurrGradient ;
      const_store<vect3d> faceCenter ;
      const_store<real> diffusionProduct ;
      const_store<real> yCurr_f ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      const_store<real> yViscosity ;
      store<real> yCurrSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToSpeciesSourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("cellcenter",cellCenter) ;
        name_store("grads(yCurr)",yCurrGradient) ;
        name_store("facecenter",faceCenter) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("yCurr_f",yCurr_f) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("yViscosity",yViscosity) ;
        name_store("yCurrSourceTerm",yCurrSourceTerm) ;
        input("ci->(cellcenter,grads(yCurr)),facecenter,diffusionProduct") ;
        input("yCurr_f,area,faceRadius,yViscosity") ;
        output("ci->yCurrSourceTerm") ;
        constraint("boundarySpeciesDiffusion,viscousFlow") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real sourceTerm=yViscosity[face]*(yCurr_f[face]*
          diffusionProduct[face]+dot(yCurrGradient[ci[face]],(area[face].n*
          area[face].sada-diffusionProduct[face]*(faceCenter[face]-
          cellCenter[ci[face]]))))*faceRadius[face] ;
        yCurrSourceTerm[ci[face]]+=sourceTerm ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToSpeciesSourceTermBoundary>
    registerDiffusiveFluxToSpeciesSourceTermBoundary ;

  // Rule to add the instantaneous species production rate to the source term.
  // Unlike the CVODE production rate which is computes every timestep, this
  // production rate is computed every iteration.
  class InstantaneousProductionRateToSpeciesSourceTerm : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_param<int> is ;
      const_store<real> cellRadius ;
      const_storeVec<real> yProductionRate ;
      store<real> yCurrSourceTerm ;
    public:

      // Define input and output.
      InstantaneousProductionRateToSpeciesSourceTerm() {
        name_store("$is{n,it,is}",is) ;
        name_store("cellRadius{n}",cellRadius) ;
        name_store("yProductionRate{n,it}",yProductionRate) ;
        name_store("yCurrSourceTerm{n,it,is}",yCurrSourceTerm) ;
        input("$is{n,it,is},cellRadius{n},yProductionRate{n,it}") ;
        output("yCurrSourceTerm{n,it,is}") ;
        constraint("geom_cells,speciesProductionInstantaneous") ;
      }

      // Add production rate to source term for a single cell.
      void calculate(Entity cell) {
        yCurrSourceTerm[cell]+=cellRadius[cell]*yProductionRate[cell][*is] ;
      }

      // Add production rate to source term for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InstantaneousProductionRateToSpeciesSourceTerm>
    registerInstantaneousProductionRateToSpeciesSourceTerm ;

  // Rule to compute the diagonal term for the linear system. Checked.
  class ComputeSpeciesMatrixDiagonal : public pointwise_rule {
    private:
      const_param<real> speciesRelaxationFactor ;
      const_store<real> yCurrMainCoefficient ;
      store<real> D ;
    public:

      // Define input and output.
      ComputeSpeciesMatrixDiagonal() {
        name_store("speciesRelaxationFactor",speciesRelaxationFactor) ;
        name_store("yCurrMainCoefficient",yCurrMainCoefficient) ;
        name_store("yCurrStar_D",D) ;
        input("speciesRelaxationFactor,yCurrMainCoefficient") ;
        output("yCurrStar_D") ;
        constraint("geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) {
        D[cell]=yCurrMainCoefficient[cell]/(*speciesRelaxationFactor) ;
      }

      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeSpeciesMatrixDiagonal>
    registerComputeSpeciesMatrixDiagonal ;

  // Rule to copy yCurrStar_D for periodic faces.
  class ComputeSpeciesMatrixDiagonalPeriodic : public pointwise_rule {
    private:
      const_Map cl,cr,pmap ;
      store<real> D ;
    public:

      // Define input and output.
      ComputeSpeciesMatrixDiagonalPeriodic() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pmap",pmap) ;
        name_store("yCurrStar_D",D) ;
        input("pmap->cl->yCurrStar_D") ;
        output("cr->yCurrStar_D") ;
      }

      // For a face.
      void calculate(Entity face) { D[cr[face]]=D[cl[pmap[face]]] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; } } ;

  register_rule<ComputeSpeciesMatrixDiagonalPeriodic>
    registerComputeSpeciesMatrixDiagonalPeriodic ;

  // Rule to initialize the lower terms for the linear system. Checked.
  class InitializeSpeciesMatrixLower : public unit_rule {
    private:
      store<real> L ;
    public:

      // Define input and output.
      InitializeSpeciesMatrixLower() {
        name_store("yCurrStar_L",L) ;
        output("yCurrStar_L") ;
        constraint("internalFaces") ;
      }

      // Initialize for a single face.
      void calculate(Entity face) { L[face]=0.0 ; }

      // Initialize for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeSpeciesMatrixLower>
    registerInitializeSpeciesMatrixLower ;

  // Rule to add the first-order inviscid flux contribution to the lower terms
  // for the linear system. Checked.
  class FOUInviscidFluxToSpeciesMatrixLower : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> massFlux ;
      store<real> L ;
    public:

      // Define input and output.
      FOUInviscidFluxToSpeciesMatrixLower() {
        name_store("massFlux",massFlux) ;
        name_store("yCurrStar_L",L) ;
        input("massFlux") ;
        output("yCurrStar_L") ;
        constraint("internalFaces") ;
      }

      // Increment the lower term for a single face. Note that the increment
      // is the negative of the one in streamUns since this coefficient is
      // for a term on the lhs of the equation. In streamUns, the coefficient
      // is associated with a term on the rhs.
      void calculate(Entity face) {
        if(massFlux[face]>0.0) L[face]-=massFlux[face] ;
      }

      // Increment the lower term for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToSpeciesMatrixLower>
    registerFOUInviscidFluxToSpeciesMatrixLower ;

  // Rule to add the diffusive flux contribution to the lower terms for the
  // linear system. Checked 07/10/04.
  class DiffusiveFluxToSpeciesMatrixLower : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> diffusionProduct ;
      const_store<real> faceRadius ;
      const_store<real> yViscosity ;
      store<real> L ;
    public:

      // Define input and output.
      DiffusiveFluxToSpeciesMatrixLower() {
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("faceRadius",faceRadius) ;
        name_store("yViscosity",yViscosity) ;
        name_store("yCurrStar_L",L) ;
        input("diffusionProduct,faceRadius,yViscosity") ;
        output("yCurrStar_L") ;
        constraint("internalFaces,viscousFlow") ;
      }

      // Increment the lower term for a single face.
      void calculate(Entity face) {
        L[face]-=yViscosity[face]*diffusionProduct[face]*faceRadius[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToSpeciesMatrixLower>
    registerDiffusiveFluxToSpeciesMatrixLower ;

  // Rule to initialize the upper terms for the linear system. Checked.
  class InitializeSpeciesMatrixUpper : public unit_rule {
    private:
      store<real> U ;
    public:

      // Define input and output.
      InitializeSpeciesMatrixUpper() {
        name_store("yCurrStar_U",U) ;
        output("yCurrStar_U") ;
        constraint("internalFaces") ;
      }

      // Initialize for a single face.
      void calculate(Entity face) { U[face]=0.0 ; }

      // Initialize for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeSpeciesMatrixUpper>
    registerInitializeSpeciesMatrixUpper ;

  // Rule to add the first-order inviscid flux contribution to the upper terms
  // for the linear system. Checked.
  class FOUInviscidFluxToSpeciesMatrixUpper : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> massFlux ;
      store<real> U ;
    public:

      // Define input and output.
      FOUInviscidFluxToSpeciesMatrixUpper() {
        name_store("massFlux",massFlux) ;
        name_store("yCurrStar_U",U) ;
        input("massFlux") ;
        output("yCurrStar_U") ;
        constraint("internalFaces") ;
      }

      // Increment the upper term for a single face. Contribution is opposite
      // to that in streamUns, as noted above for the lower terms.
      void calculate(Entity face) {
        if(massFlux[face]<=0.0) U[face]+=massFlux[face] ;
      }

      // Increment the upper term for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToSpeciesMatrixUpper>
    registerFOUInviscidFluxToSpeciesMatrixUpper ;

  // Rule to add the diffusive flux contribution to the upper terms for the
  // linear system. Checked 07/10/04.
  class DiffusiveFluxToSpeciesMatrixUpper : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> diffusionProduct ;
      const_store<real> faceRadius ;
      const_store<real> yViscosity ;
      store<real> U ;
    public:

      // Define input and output.
      DiffusiveFluxToSpeciesMatrixUpper() {
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("faceRadius",faceRadius) ;
        name_store("yViscosity",yViscosity) ;
        name_store("yCurrStar_U",U) ;
        input("diffusionProduct,faceRadius,yViscosity") ;
        output("yCurrStar_U") ;
        constraint("internalFaces,viscousFlow") ;
      }

      // Increment the upper term for a single face.
      void calculate(Entity face) {
        U[face]-=yViscosity[face]*diffusionProduct[face]*faceRadius[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToSpeciesMatrixUpper>
    registerDiffusiveFluxToSpeciesMatrixUpper ;

  // Rule to compute the right-hand side for the linear system. Checked.
  class ComputeSpeciesRHS : public pointwise_rule {
    private:
      const_param<real> speciesRelaxationFactor ;
      const_store<real> yCurr ;
      const_store<real> yCurrMainCoefficient ;
      const_store<real> yCurrSourceTerm ;
      store<real> B ;
    public:

      // Define input and output.
      ComputeSpeciesRHS() {
        name_store("speciesRelaxationFactor",speciesRelaxationFactor) ;
        name_store("yCurr",yCurr) ;
        name_store("yCurrMainCoefficient",yCurrMainCoefficient) ;
        name_store("yCurrSourceTerm",yCurrSourceTerm) ;
        name_store("yCurrStar_B",B) ;
        input("speciesRelaxationFactor,yCurr") ;
        input("yCurrMainCoefficient,yCurrSourceTerm") ;
        output("yCurrStar_B") ;
        constraint("geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) {
        B[cell]=yCurrSourceTerm[cell]+(1.0-(*speciesRelaxationFactor))*
          yCurrMainCoefficient[cell]*yCurr[cell]/(*speciesRelaxationFactor) ;
      }

      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeSpeciesRHS> registerComputeSpeciesRHS ;

  // Priority rule for inviscid flow to compute the right-hand side for the
  // linear system. Includes the net mass flux term which cancels the term
  // added to the main coefficient. Checked.
  class ComputeSpeciesRHSInviscid : public pointwise_rule {
    private:
      const_param<real> speciesRelaxationFactor ;
      const_store<real> yCurr ;
      const_store<real> yCurrMainCoefficient ;
      const_store<real> yCurrSourceTerm ;
      const_store<real> netMassFlux ;
      store<real> B ;
    public:

      // Define input and output.
      ComputeSpeciesRHSInviscid() {
        name_store("speciesRelaxationFactor",speciesRelaxationFactor) ;
        name_store("yCurr",yCurr) ;
        name_store("yCurrMainCoefficient",yCurrMainCoefficient) ;
        name_store("yCurrSourceTerm",yCurrSourceTerm) ;
        name_store("netMassFlux",netMassFlux) ;
        name_store("inviscidFlow::yCurrStar_B",B) ;
        input("speciesRelaxationFactor,yCurr") ;
        input("yCurrMainCoefficient,yCurrSourceTerm,netMassFlux") ;
        output("inviscidFlow::yCurrStar_B") ;
        constraint("inviscidFlow,geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) {
        B[cell]=yCurrSourceTerm[cell]+(1.0-(*speciesRelaxationFactor))*
          yCurrMainCoefficient[cell]*yCurr[cell]/(*speciesRelaxationFactor) ;
        if(netMassFlux[cell]<0.0) B[cell]-=netMassFlux[cell]*yCurr[cell] ;
      }

      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeSpeciesRHSInviscid> registerComputeSpeciesRHSInviscid ;

  // Rule to copy yCurrStar_B for periodic faces.
  class ComputeSpeciesRHSPeriodic : public pointwise_rule {
    private:
      const_Map cl,cr,pmap ;
      store<real> B ;
    public:

      // Define input and output.
      ComputeSpeciesRHSPeriodic() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pmap",pmap) ;
        name_store("yCurrStar_B",B) ;
        input("pmap->cl->yCurrStar_B") ;
        output("cr->yCurrStar_B") ;
      }

      // For a face.
      void calculate(Entity face) { B[cr[face]]=B[cl[pmap[face]]] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeSpeciesRHSPeriodic> registerComputeSpeciesRHSPeriodic ;

//-----------------------------------------------------------------------------
// Rules for computing the residual of the current species equation.

  // Rule to initialize the species residual. Checked.
  class InitializeSpeciesResidual : public unit_rule {
    private:
      store<real> speciesResidual ;
    public:

      // Define input and output.
      InitializeSpeciesResidual() {
        name_store("speciesResidual",speciesResidual) ;
        output("speciesResidual") ;
        constraint("speciesTransport,vol") ;
      }

      // Initialize the residual for a single cell.
      void calculate(Entity cell) { speciesResidual[cell]=0.0 ; }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeSpeciesResidual> registerInitializeSpeciesResidual ;

  // Rule to add cell contributions to the current species residual. Checked.
  class ComputeSpeciesResidualOne : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> rhoScale,vScale,lScale ;
      const_store<real> D,yCurr,B ;
      store<real> speciesResidual ;
    private:
      real yFactor ;
    public:

      // Define input and output.
      ComputeSpeciesResidualOne() {
        name_store("rhoScale",rhoScale) ;
        name_store("vScale",vScale) ;
        name_store("lScale",lScale) ;
        name_store("yCurrStar_D",D) ;
        name_store("yCurr",yCurr) ;
        name_store("yCurrStar_B",B) ;
        name_store("speciesResidual",speciesResidual) ;
        input("rhoScale,vScale,lScale,yCurrStar_D,yCurr,yCurrStar_B") ;
        output("speciesResidual") ;
        constraint("vol") ;
      }

      // Initialize the residual for a single cell.
      void calculate(Entity cell) {
//cout << "cell,B,D,yCurr: " << cell << " " << B[cell] << " " << D[cell]
//<< " " << yCurr[cell] << endl ;
        speciesResidual[cell]+=(B[cell]-D[cell]*yCurr[cell])/yFactor ;
      }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) {
        yFactor=(*rhoScale)*(*vScale)*(*lScale)*(*lScale) ; do_loop(seq,this) ;
      }
  } ;

  register_rule<ComputeSpeciesResidualOne> registerComputeSpeciesResidualOne ;

  // Rule to add face contributions to the current species residual. Checked.
  class ComputeSpeciesResidualTwo : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_param<real> rhoScale,vScale,lScale ;
      const_store<real> yCurr,L,U ;
      store<real> speciesResidual ;
    private:
      real yFactor ;
    public:

      // Define input and output.
      ComputeSpeciesResidualTwo() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("rhoScale",rhoScale) ;
        name_store("vScale",vScale) ;
        name_store("lScale",lScale) ;
        name_store("yCurr",yCurr) ;
        name_store("yCurrStar_L",L) ;
        name_store("yCurrStar_U",U) ;
        name_store("speciesResidual",speciesResidual) ;
        input("rhoScale,vScale,lScale,(cl,cr)->yCurr,yCurrStar_L,yCurrStar_U") ;
        output("(cl,cr)->speciesResidual") ;
        constraint("internalFaces") ;
      }

      // Add the neighbor contribution to the residual for each of the two
      // cells on either side of the face.
      void calculate(Entity face) {
        speciesResidual[cl[face]]-=U[face]*yCurr[cr[face]]/yFactor ;
        speciesResidual[cr[face]]-=L[face]*yCurr[cl[face]]/yFactor ;
      }

      // Add the neighbor contributions for a sequence of faces.
      virtual void compute(const sequence &seq) {
        yFactor=(*rhoScale)*(*vScale)*(*lScale)*(*lScale) ; do_loop(seq,this) ;
      }
  } ;

  register_rule<ComputeSpeciesResidualTwo> registerComputeSpeciesResidualTwo ;

  // Rule to initialize the total residual for the current species. Checked.
  class InitializeTotalCurrentSpeciesResidual : public unit_rule {
    private:
      param<real> totalCurrentSpeciesResidual ;
    public:

      // Define input and output.
      InitializeTotalCurrentSpeciesResidual() {
        name_store("totalCurrentSpeciesResidual",totalCurrentSpeciesResidual) ;
        output("totalCurrentSpeciesResidual") ;
        constraint("speciesTransport,geom_cells") ;
      }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) {
        *totalCurrentSpeciesResidual=0.0 ;
      }
  } ;

  register_rule<InitializeTotalCurrentSpeciesResidual>
    registerInitializeTotalCurrentSpeciesResidual ;

  // Rule to compute the total current species residual. Checked.
  class ComputeTotalCurrentSpeciesResidual : public apply_rule<param<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> speciesResidual ;
      param<real> totalCurrentSpeciesResidual ;
    public:

      // Define input and output.
      ComputeTotalCurrentSpeciesResidual() {
        name_store("speciesResidual",speciesResidual) ;
        name_store("totalCurrentSpeciesResidual",totalCurrentSpeciesResidual) ;
        input("speciesResidual") ;
        output("totalCurrentSpeciesResidual") ;
        constraint("geom_cells") ;
      }

      // Add the cell contribution to the residual for a single cell.
      void calculate(Entity cell) {
        *totalCurrentSpeciesResidual+=(speciesResidual[cell]<0.0)?
          -speciesResidual[cell]:speciesResidual[cell] ;
      }

      // Add the cell contribution to the residual for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeTotalCurrentSpeciesResidual>
    registerComputeTotalCurrentSpeciesResidual ;

//-----------------------------------------------------------------------------
// Rule to set the species residuals to zero for flows with no species
// transport. Checked.

  class TotalSpeciesResidualNoSpeciesTransport : public singleton_rule {
    private:
      const_param<int> numSpecies ;
      param<vector<real> > totalSpeciesResidual ;
    public:

      // Define input and output.
      TotalSpeciesResidualNoSpeciesTransport() {
        name_store("numSpecies",numSpecies) ;
        name_store("totalSpeciesResidual",totalSpeciesResidual) ;
        input("numSpecies") ;
        output("totalSpeciesResidual") ;
        constraint("noSpeciesTransport,geom_cells") ;
      }

      // Set the residual.
      virtual void compute(const sequence &seq) {
        *totalSpeciesResidual=vector<real>(*numSpecies,0.0) ;
      }
  } ;

  register_rule<TotalSpeciesResidualNoSpeciesTransport>
    registerTotalSpeciesResidualNoSpeciesTransport ;

//-----------------------------------------------------------------------------
// Rules for computing the species residuals for flows with species transport.

  // Build rule.
  class TotalSpeciesResidualBuild : public singleton_rule {
    private:
      const_param<int> numSpecies ;
      param<vector<real> > totalSpeciesResidual ;
    public:

      // Define input and output.
      TotalSpeciesResidualBuild() {
        name_store("numSpecies{n,it}",numSpecies) ;
        name_store("totalSpeciesResidual{n,it,is=0}",totalSpeciesResidual) ;
        input("numSpecies{n,it}") ;
        output("totalSpeciesResidual{n,it,is=0}") ;
        constraint("speciesTransport{n,it},geom_cells{n,it}") ;
      }

      // Initialize the residual.
      virtual void compute(const sequence &seq) {
        *totalSpeciesResidual=vector<real>(*numSpecies,0.0) ;
      }
  } ;

  register_rule<TotalSpeciesResidualBuild> registerTotalSpeciesResidualBuild ;

  // Advance rule. Note that since we cannot use an update-in-place rule for
  // param variables (bug), we must copy the old values.
  class TotalSpeciesResidualAdvance : public singleton_rule {
    private:
      const_param<int> is ;
      const_param<real> totalCurrentSpeciesResidual ;
      const_param<vector<real> > totalSpeciesResidualOld ;
      param<vector<real> > totalSpeciesResidualNew ;
    public:

      // Define input and output.
      TotalSpeciesResidualAdvance() {
        name_store("$is{n,it,is}",is) ;
        name_store("totalCurrentSpeciesResidual{n,it,is}",
          totalCurrentSpeciesResidual) ;
        name_store("totalSpeciesResidual{n,it,is}",totalSpeciesResidualOld) ;
        name_store("totalSpeciesResidual{n,it,is+1}",totalSpeciesResidualNew) ;
        input("$is{n,it,is},totalCurrentSpeciesResidual{n,it,is}") ;
        input("totalSpeciesResidual{n,it,is}") ;
        output("totalSpeciesResidual{n,it,is+1}");
        constraint("geom_cells{n,it,is},speciesTransport{n,it,is}") ;
      }

      // Save current species total residual.
      void compute(const sequence &seq) {
         *totalSpeciesResidualNew=*totalSpeciesResidualOld ;
        (*totalSpeciesResidualNew)[*is]=*totalCurrentSpeciesResidual ;
      }
  } ;

  register_rule<TotalSpeciesResidualAdvance>
    registerTotalSpeciesResidualAdvance ;

  // Collapse rule.
  class TotalSpeciesResidualCollapse : public singleton_rule {
    private:
      const_param<bool> iterationFinished ;
      const_param<vector<real> > totalSpeciesResidual ;
      param<vector<real> > totalSpeciesResidualFinal ;
    public:

      // Define input and output.
      TotalSpeciesResidualCollapse() {
        name_store("speciesIterationFinished{n,it,is}",iterationFinished) ;
        name_store("totalSpeciesResidual{n,it,is}",totalSpeciesResidual) ;
        name_store("totalSpeciesResidual{n,it}",totalSpeciesResidualFinal) ;
        input("speciesIterationFinished{n,it,is}") ;
        input("totalSpeciesResidual{n,it,is}") ; // Trouble is here.
        output("totalSpeciesResidual{n,it}") ;
        conditional("speciesIterationFinished{n,it,is}") ;
        constraint("geom_cells{n,it,is},speciesTransport{n,it,is}") ;
      }

      // Copy residual.
      void compute(const sequence &seq) {
        *totalSpeciesResidualFinal=*totalSpeciesResidual ;
      }
  } ;

  register_rule<TotalSpeciesResidualCollapse>
    registerTotalSpeciesResidualCollapse ;

//-----------------------------------------------------------------------------
// Rules for computing the energy equation source term due to species
// diffusion.

  // Rule to extract current species enthalpy from EOS for the cells. This rule
  // was added to prevent having to define eos_state and eos_mixture_state over
  // "vol" for cases with periodic boundaries, which causes trouble with
  // generating a schedule due to the heavy use of rename rules with these
  // variables used to minimize memory usage.
  class CurrentSpeciesEnthalpyInterior : public pointwise_rule {
    private:
      const_param<int> is ;
      const_param<EOS> eos ;
      const_store<EOS::State> eosState ;
      const_storeVec<real> eosMixtureState ;
      store<real> yCurrEnthalpyCell ;
    public:

      // Define input and output.
      CurrentSpeciesEnthalpyInterior() {
        name_store("$is",is) ;
        name_store("eos",eos) ;
        name_store("eos_state",eosState) ;
        name_store("eos_mixture_state",eosMixtureState) ;
        name_store("yCurrEnthalpyCell",yCurrEnthalpyCell) ;
        input("$is,eos,eos_state,eos_mixture_state") ;
        output("yCurrEnthalpyCell") ;
        constraint("geom_cells") ;
      }

      // Compute the species enthalpy for the cell.
      void calculate(Entity cell) {
        yCurrEnthalpyCell[cell]=eos->speciesEnthalpy(*is,eosState[cell],
          eosMixtureState[cell]) ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CurrentSpeciesEnthalpyInterior>
    registerCurrentSpeciesEnthalpyInterior ;

  // Build rule.
  class SpeciesEnthalpyDiffusionBuild : public pointwise_rule {
    private:
      store<real> hDiffusion ;
    public:

      // Define input and output.
      SpeciesEnthalpyDiffusionBuild() {
        name_store("hDiffusion{n,it,is=0}",hDiffusion) ;
        output("hDiffusion{n,it,is=0}") ;
        constraint("geom_cells{n,it}") ;
      }

      // Initialize for a single cell.
      void calculate(Entity cell) { hDiffusion[cell]=0.0 ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SpeciesEnthalpyDiffusionBuild>
    registerSpeciesEnthalpyDiffusionBuild ;

  // Rule to compute the enthalpy for a species for interior faces.
  class SpeciesEnthalpyInterior : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<real> yCurrEnthalpyCell ;
      store<real> yCurrEnthalpy ;
    public:

      // Define input and output.
      SpeciesEnthalpyInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("yCurrEnthalpyCell",yCurrEnthalpyCell) ;
        name_store("yCurrEnthalpy",yCurrEnthalpy) ;
        input("(cl,cr)->yCurrEnthalpyCell") ;
        output("yCurrEnthalpy") ;
        constraint("internalFaces") ;
      }

      // Compute the species enthalpy for the face.
      void calculate(Entity face) {
        yCurrEnthalpy[face]=0.5*(yCurrEnthalpyCell[cl[face]]+
          yCurrEnthalpyCell[cr[face]]) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SpeciesEnthalpyInterior> registerSpeciesEnthalpyInterior ;

  // Rule to compute the enthalpy for a species for symmetry faces.
  class SpeciesEnthalpySymmetry : public pointwise_rule {
    private:
      const_Map ci ;
      const_param<int> is ;
      const_param<EOS> eos ;
      const_store<EOS::State> eosState ;
      const_storeVec<real> eosMixtureState ;
      store<real> yCurrEnthalpy_f ;
    public:

      // Define input and output.
      SpeciesEnthalpySymmetry() {
        name_store("ci",ci) ;
        name_store("$is",is) ;
        name_store("eos",eos) ;
        name_store("eos_state",eosState) ;
        name_store("eos_mixture_state",eosMixtureState) ;
        name_store("symmetry::yCurrEnthalpy_f",yCurrEnthalpy_f) ;
        input("$is,eos,ci->(eos_state,eos_mixture_state)") ;
        output("symmetry::yCurrEnthalpy_f") ;
        constraint("symmetry_BC") ;
      }

      // Compute the species enthalpy for the face.
      void calculate(Entity face) {
        yCurrEnthalpy_f[face]=eos->speciesEnthalpy(*is,eosState[ci[face]],
          eosMixtureState[ci[face]]) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SpeciesEnthalpySymmetry> registerSpeciesEnthalpySymmetry ;

  // Rule to compute the enthalpy for a species for boundary faces.
  class SpeciesEnthalpyBoundary : public pointwise_rule {
    private:
      const_param<int> is ;
      const_param<EOS> eos ;
      const_store<EOS::State> eosState_f ;
      const_storeVec<real> eosMixtureState_f ;
      store<real> yCurrEnthalpy_f ;
    public:

      // Define input and output.
      SpeciesEnthalpyBoundary() {
        name_store("$is",is) ;
        name_store("eos",eos) ;
        name_store("eos_state_f",eosState_f) ;
        name_store("eos_mixture_state_f",eosMixtureState_f) ;
        name_store("yCurrEnthalpy_f",yCurrEnthalpy_f) ;
        input("$is,eos,eos_state_f,eos_mixture_state_f") ;
        output("yCurrEnthalpy_f") ;
//      constraint("boundaryFaces") ; 1/28/2009
      }

      // Compute the species enthalpy for the face.
      void calculate(Entity face) {
        yCurrEnthalpy_f[face]=eos->speciesEnthalpy(*is,eosState_f[face],
          eosMixtureState_f[face]) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SpeciesEnthalpyBoundary> registerSpeciesEnthalpyBoundary ;

  // Unit rule to compute the diffusion of enthalpy due to the diffusion of the
  // current species.
  class CurrentSpeciesEnthalpyDiffusionUnit : public unit_rule {
    private:
      store<real> yCurrHDiffusion ;
    public:

      // Define input and output.
      CurrentSpeciesEnthalpyDiffusionUnit() {
        name_store("yCurrHDiffusion",yCurrHDiffusion) ;
        output("yCurrHDiffusion") ;
        constraint("viscousFlow,vol") ;
      }

      // Compute the enthalpy diffusion into each cell
      void calculate(Entity cell) { yCurrHDiffusion[cell]=0.0 ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CurrentSpeciesEnthalpyDiffusionUnit>
    registerCurrentSpeciesEnthalpyDiffusionUnit ;

  // Apply rule to compute the diffusion of enthalpy due to the diffusion of
  // the current species for interion faces.
  class CurrentSpeciesEnthalpyDiffusionInterior : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> yCurr ;
      const_store<real> vol ;
      const_store<vect3d> yCurrGradient ;
      const_store<real> yCurrEnthalpy ;
      const_store<real> diffusionProduct ;
      const_store<vect3d> geometryFactor0 ;
      const_store<real> faceRadius ;
      const_store<real> yViscosity ;
      store<real> yCurrHDiffusion ;
    public:

      // Define input and output.
      CurrentSpeciesEnthalpyDiffusionInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("yCurr",yCurr) ;
        name_store("vol",vol) ;
        name_store("grads(yCurr)",yCurrGradient) ;
        name_store("yCurrEnthalpy",yCurrEnthalpy) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("geometryFactor0",geometryFactor0) ;
        name_store("faceRadius",faceRadius) ;
        name_store("yViscosity",yViscosity) ;
        name_store("yCurrHDiffusion",yCurrHDiffusion) ;
        input("(cl,cr)->(grads(yCurr),yCurr,vol),yViscosity") ;
        input("diffusionProduct,faceRadius,yCurrEnthalpy,geometryFactor0") ;
        output("(cl,cr)->yCurrHDiffusion") ;
        constraint("internalFaces") ;
      }

      // Compute the enthalpy diffusion into the cells. Note that we are now
      // volume averaging the cell gradients to get the face value.
      void calculate(Entity face) {
        real sourceTerm=yViscosity[face]*((yCurr[cr[face]]-
          yCurr[cl[face]])*diffusionProduct[face]+dot((vol[cr[face]]*
          yCurrGradient[cl[face]]+vol[cl[face]]*yCurrGradient[cr[face]])/
          (vol[cl[face]]+vol[cr[face]]),geometryFactor0[face]))*
          yCurrEnthalpy[face]*faceRadius[face] ;
        yCurrHDiffusion[cl[face]]+=sourceTerm ;
        yCurrHDiffusion[cr[face]]-=sourceTerm ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CurrentSpeciesEnthalpyDiffusionInterior>
    registerCurrentSpeciesEnthalpyDiffusionInterior ;

  // Apply rule to compute the diffusion of enthalpy due to the diffusion of
  // the current species for boundary faces. Enthalpy diffusion due to species
  // diffusion only occurs on inlets and outlets.
  class CurrentSpeciesEnthalpyDiffusionBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> yCurr ;
      const_store<real> yCurr_f ;
      const_store<real> diffusionProduct ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<vect3d> yCurrGradient ;
      const_store<real> yCurrEnthalpy_f ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      const_store<real> yViscosity ;
      store<real> yCurrHDiffusion ;
    public:

      // Define input and output.
      CurrentSpeciesEnthalpyDiffusionBoundary() {
        name_store("ci",ci) ;
        name_store("yCurr",yCurr) ;
        name_store("yCurr_f",yCurr_f) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("grads(yCurr)",yCurrGradient) ;
        name_store("yCurrEnthalpy_f",yCurrEnthalpy_f) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("yViscosity",yViscosity) ;
        name_store("yCurrHDiffusion",yCurrHDiffusion) ;
        input("ci->(grads(yCurr),yCurr,cellcenter)") ;
        input("yCurr_f,diffusionProduct,facecenter") ;
        input("area,faceRadius,yCurrEnthalpy_f,yViscosity") ;
        output("ci->yCurrHDiffusion") ;
        constraint("boundarySpeciesDiffusion") ;
      }

      // Compute the enthalpy diffusion into the cells. 
      void calculate(Entity face) {
        real sourceTerm=yViscosity[face]*((yCurr_f[face]-yCurr[ci[face]])*
          diffusionProduct[face]+dot(yCurrGradient[ci[face]],(area[face].n*
          area[face].sada-diffusionProduct[face]*(faceCenter[face]-
          cellCenter[ci[face]]))))*yCurrEnthalpy_f[face]*faceRadius[face] ;
        yCurrHDiffusion[ci[face]]+=sourceTerm ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CurrentSpeciesEnthalpyDiffusionBoundary>
    registerCurrentSpeciesEnthalpyDiffusionBoundary ;

  // Advance rule.
  class SpeciesEnthalpyDiffusionAdvance : public pointwise_rule {
    private:
      const_param<int> is ;
      const_store<real> yCurrHDiffusion ;
      store<real> hDiffusion ;
    public:

      // Define input and output.
      SpeciesEnthalpyDiffusionAdvance() {
        name_store("$is{n,it,is}",is) ;
        name_store("yCurrHDiffusion{n,it,is}",yCurrHDiffusion) ;
        name_store("hDiffusion{n,it,is}",hDiffusion) ;
        input("$is{n,it,is},yCurrHDiffusion{n,it,is},hDiffusion{n,it,is}") ;
        output("hDiffusion{n,it,is+1}=hDiffusion{n,it,is}");
        constraint("geom_cells{n,it,is}") ;
      }

      // Add the enthalpy diffusion for the current species for a single cell.
      void calculate(Entity cell) { hDiffusion[cell]+=yCurrHDiffusion[cell] ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SpeciesEnthalpyDiffusionAdvance>
    registerSpeciesEnthalpyDiffusionAdvance ;

  // Collapse rule.
  class SpeciesEnthalpyDiffusionCollapse : public pointwise_rule {
    private:
      const_param<bool> iterationFinished ;
      store<real> hDiffusion ;
    public:

      // Define input and output.
      SpeciesEnthalpyDiffusionCollapse() {
        name_store("speciesIterationFinished{n,it,is}",iterationFinished) ;
        name_store("hDiffusion{n,it,is}",hDiffusion) ;
        input("speciesIterationFinished{n,it,is}") ;
        input("hDiffusion{n,it,is}") ;
        output("hDiffusion{n,it}=hDiffusion{n,it,is}") ;
        conditional("speciesIterationFinished{n,it,is}") ;
        constraint("geom_cells{n,it,is}") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<SpeciesEnthalpyDiffusionCollapse>
    registerSpeciesEnthalpyDiffusionCollapse ;

//-----------------------------------------------------------------------------
// Rules for marching the species.

  // Time build rule for species mass fractions when using BDF2 time integrator.
  class TimeBuildSpeciesBDF2 : public pointwise_rule {
    private:
      const_param<int> numSpecies ;
      const_storeVec<real> y_ic ;
      storeVec<real> y ;
    public:

      // Define input and output.
      TimeBuildSpeciesBDF2() {
        name_store("numSpecies",numSpecies) ;
        name_store("y_ic",y_ic) ;
        name_store("y{n=-1}",y) ;
        input("numSpecies,y_ic") ;
        output("y{n=-1}") ;
        constraint("speciesTransport,geom_cells") ;
      }

      // Assign species mass fractions at time zero for a single cell.
      void calculate(Entity cell) { y[cell]=y_ic[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) {
        y.setVecSize(*numSpecies) ; do_loop(seq,this) ;
      }
  } ;

  register_rule<TimeBuildSpeciesBDF2> registerTimeBuildSpeciesBDF2 ;

  // Time build rule for species mass fractions. Checked.
  class TimeBuildSpecies : public pointwise_rule {
    private:
      const_param<int> numSpecies ;
      const_storeVec<real> y_ic ;
      storeVec<real> yTimeStepZero ;
    public:

      // Define input and output.
      TimeBuildSpecies() {
        name_store("numSpecies",numSpecies) ;
        name_store("y_ic",y_ic) ;
        name_store("y{n=0}",yTimeStepZero) ;
        input("numSpecies,y_ic") ;
        output("y{n=0}") ;
        constraint("speciesTransport,geom_cells") ;
      }

      // Assign species mass fractions at time zero for a single cell.
      void calculate(Entity cell) { yTimeStepZero[cell]=y_ic[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) {
        yTimeStepZero.setVecSize(*numSpecies) ; do_loop(seq,this) ;
      }
  } ;

  register_rule<TimeBuildSpecies> registerTimeBuildSpecies ;

  // Iteration build rule for species mass fractions. Checked.
  class IterationBuildSpecies : public pointwise_rule {
    private:
      const_param<int> numSpecies ;
      const_storeVec<real> yTimeStepN ;
      storeVec<real> yIterationZero ;
    public:

      // Define input and output.
      IterationBuildSpecies() {
        name_store("numSpecies{n}",numSpecies) ;
        name_store("y{n}",yTimeStepN) ;
        name_store("y{n,it=0}",yIterationZero) ;
        input("numSpecies{n},y{n}") ;
        output("y{n,it=0}") ;
        constraint("speciesTransport{n},geom_cells{n}") ;
      }

      // Assign species mass fractions for a single cell.
      void calculate(Entity cell) { yIterationZero[cell]=yTimeStepN[cell] ; }

      // Assign species mass fractions for a sequence of cells.
      void compute(const sequence &seq) {
        yIterationZero.setVecSize(*numSpecies) ; do_loop(seq,this) ;
      }
  } ;

  register_rule<IterationBuildSpecies> registerIterationBuildSpecies ;

  // Iteration build rule for species mass fractions when using first-order
  // Strang splitting.
  class IterationBuildSpeciesStrang : public pointwise_rule {
    private:
      const_param<int> numSpecies ;
      const_storeVec<real> yODE ;
      storeVec<real> yIterationZero ;
    public:

      // Define input and output.
      IterationBuildSpeciesStrang() {
        name_store("numSpecies{n}",numSpecies) ;
        name_store("yODE{n}",yODE) ;
        name_store("strang::y{n,it=0}",yIterationZero) ;
        input("numSpecies{n},yODE{n}") ;
        output("strang::y{n,it=0}") ;
        constraint("speciesTransport{n},geom_cells{n}") ;
        constraint("strangSplitting{n}") ;
      }

      // Assign species mass fractions for a single cell.
      void calculate(Entity cell) { yIterationZero[cell]=yODE[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) {
        yIterationZero.setVecSize(*numSpecies) ; do_loop(seq,this) ;
      }
  } ;

  register_rule<IterationBuildSpeciesStrang>
    registerIterationBuildSpeciesStrang ;

  // Species iteration build rule.
  class SpeciesIterationLoopBuild : public pointwise_rule {
    private:
      const_param<int> numSpecies ;
      const_storeVec<real> yTimeStepNIterationIt ;
      storeVec<real> ySpeciesIterationZero ;
    public:

      // Define input and output.
      SpeciesIterationLoopBuild() {
        name_store("numSpecies{n,it}",numSpecies) ;
        name_store("y{n,it}",yTimeStepNIterationIt) ;
        name_store("yStar{n,it,is=0}",ySpeciesIterationZero) ;
        input("numSpecies{n,it},y{n,it}") ;
        output("yStar{n,it,is=0}") ;
        constraint("speciesTransport{n,it},geom_cells{n,it}") ;
      }

      // Assign species mass fractions for a single cell.
      void calculate(Entity cell) {
        ySpeciesIterationZero[cell]=yTimeStepNIterationIt[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) {
        ySpeciesIterationZero.setVecSize(*numSpecies) ; do_loop(seq,this) ;
      }
  } ;

  register_rule<SpeciesIterationLoopBuild>
    registerSpeciesIterationLoopBuild ;

  // Rule to add temporal component of species equation to the main coefficient.
  class TemporalToSpeciesMainCoefficient: public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> timeIntegratorFactor0 ;
      const_store<real> timeStepFactor ;
      const_store<real> rho ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> yCurrMainCoefficient ;
    public:

      // Define input and output.
      TemporalToSpeciesMainCoefficient() {
        name_store("timeIntegratorFactor0{n}",timeIntegratorFactor0) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n}",rho) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("yCurrMainCoefficient{n,it}",yCurrMainCoefficient) ;
        input("rho{n},vol{n},cellRadius{n,it}") ;
        input("timeStepFactor{n},timeIntegratorFactor0{n}") ;
        output("yCurrMainCoefficient{n,it}") ;
        constraint("geom_cells{n,it}") ;
      }

      // Add temporal component for a single cell.
      void calculate(Entity cell) {
        yCurrMainCoefficient[cell]+=rho[cell]*vol[cell]*cellRadius[cell]*
          (*timeIntegratorFactor0)/timeStepFactor[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToSpeciesMainCoefficient>
    registerTemporalToSpeciesMainCoefficient ;

  // Rule to add temporal component of species equation to the source term.
  class TemporalToSpeciesSourceTerm : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<int> is ;
      const_param<real> timeIntegratorFactor0 ;
      const_store<real> timeStepFactor ;
      const_store<real> rho ;
      const_storeVec<real> y ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> yCurrSourceTerm ;
    public:

      // Define input and output.
      TemporalToSpeciesSourceTerm() {
        name_store("$is{n,it,is}",is) ;
        name_store("timeIntegratorFactor0{n}",timeIntegratorFactor0) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n}",rho) ;
        name_store("y{n}",y) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("yCurrSourceTerm{n,it,is}",yCurrSourceTerm) ;
        input("$is{n,it,is},rho{n},y{n},vol{n},cellRadius{n,it}") ;
        input("timeStepFactor{n},timeIntegratorFactor0{n}") ;
        output("yCurrSourceTerm{n,it,is}") ;
        constraint("geom_cells") ;
      }

      // Add temporal component to source term for a single cell.
      void calculate(Entity cell) {
        yCurrSourceTerm[cell]+=(rho[cell]*vol[cell]*cellRadius[cell]*
          (*timeIntegratorFactor0)/timeStepFactor[cell])*y[cell][*is] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToSpeciesSourceTerm>
    registerTemporalToSpeciesSourceTerm ;

  // Rule to add temporal component of species equation to the source term.
  class TemporalToSpeciesSourceTermBDF2 : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<int> is ;
      const_param<real> timeIntegratorFactor1 ;
      const_store<real> timeStepFactor ;
      const_store<real> rhoOld ;
      const_storeVec<real> yOld ;
      const_store<real> yCurr ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> yCurrSourceTerm ;
    public:

      // Define input and output.
      TemporalToSpeciesSourceTermBDF2() {
        name_store("$is{n,it,is}",is) ;
        name_store("timeIntegratorFactor1{n}",timeIntegratorFactor1) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n-1}",rhoOld) ;
        name_store("y{n-1}",yOld) ;
        name_store("yCurr{n,it,is}",yCurr) ;
        name_store("vol{n,it}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("yCurrSourceTerm{n,it,is}",yCurrSourceTerm) ;
        input("$is{n,it,is},rho{n-1},y{n-1},yCurr{n,it,is},vol{n,it}") ;
        input("cellRadius{n,it},timeStepFactor{n}") ;
        input("timeIntegratorFactor1{n}") ;
        output("yCurrSourceTerm{n,it,is}") ;
        constraint("geom_cells,BDF2Integrator") ;
      }

      // Add temporal component to source term for a single cell.
      void calculate(Entity cell) {
        yCurrSourceTerm[cell]+=((*timeIntegratorFactor1)*rhoOld[cell]*
          vol[cell]*cellRadius[cell]/timeStepFactor[cell])*
          (yCurr[cell]-yOld[cell][*is]) ;
      }

      // Add temporal component to source term for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToSpeciesSourceTermBDF2>
    registerTemporalToSpeciesSourceTermBDF2 ;

  // For Strang splitting, add a correction so that we get y{n+1}-yODE{n}
  // rather than y{n+1}-y{n} for BDF. For BDF2, the correction will give
  // us 3y{n+1}-4yODE{n}+y{n-1} instead of 3y{n+1}-4y{n}+y{n-1} . See
  // notes of March. 12, 2009.
  class TemporalToSpeciesSourceTermStrang : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<int> is ;
      const_store<real> timeStepFactor ;
      const_param<real> timeIntegratorFactor0 ;
      const_store<real> rho ;
      const_storeVec<real> y,yODE ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> yCurrSourceTerm ;
    public:

      // Define input and output.
      TemporalToSpeciesSourceTermStrang() {
        name_store("$is{n,it,is}",is) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("timeIntegratorFactor0{n}",timeIntegratorFactor0) ;
        name_store("rho{n}",rho) ;
        name_store("y{n}",y) ;
        name_store("yODE{n}",yODE) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("yCurrSourceTerm{n,it,is}",yCurrSourceTerm) ;
        input("$is{n,it,is},rho{n},y{n},yODE{n},vol{n},cellRadius{n,it}") ;
        input("timeStepFactor{n},timeIntegratorFactor0{n}") ;
        output("yCurrSourceTerm{n,it,is}") ;
        constraint("strangSplitting,geom_cells") ;
      }

      // Add temporal component to source term for a single cell.
      void calculate(Entity cell) {
        yCurrSourceTerm[cell]+=(rho[cell]*vol[cell]*cellRadius[cell]/
          timeStepFactor[cell])*(*timeIntegratorFactor0)*
          (yODE[cell][*is]-y[cell][*is]) ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToSpeciesSourceTermStrang>
    registerTemporalToSpeciesSourceTermStrang ;

  // Rule to add the species production rate to the source term. This is put,
  // here since we will have different rules for different time-integration
  // schemes.
  class ProductionRateToSpeciesSourceTerm : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<int> is ;
      const_store<real> cellRadius ;
      const_storeVec<real> yProductionRate ;
      store<real> yCurrSourceTerm ;
    public:

      // Define input and output.
      ProductionRateToSpeciesSourceTerm() {
        name_store("$is{n,it,is}",is) ;
        name_store("cellRadius{n}",cellRadius) ;
        name_store("yProductionRate{n}",yProductionRate) ;
        name_store("yCurrSourceTerm{n,it,is}",yCurrSourceTerm) ;
        input("$is{n,it,is},cellRadius{n},yProductionRate{n}") ;
        output("yCurrSourceTerm{n,it,is}") ;
        constraint("geom_cells,speciesProductionCVODE,sourceSplitting") ;
      }

      // Add production rate to source term for a single cell.
      void calculate(Entity cell) {
        yCurrSourceTerm[cell]+=cellRadius[cell]*yProductionRate[cell][*is] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ProductionRateToSpeciesSourceTerm>
    registerProductionRateToSpeciesSourceTerm ;

  // Iteration advance rule for species mass fractions.
  class IterationAdvanceSpecies : public pointwise_rule {
    private:
      const_param<int> numSpecies ;
      const_storeVec<real> yStar ;
      storeVec<real> yIterationPlusOne ;
    public:

      // Define input and output.
      IterationAdvanceSpecies() {
        name_store("numSpecies{n,it}",numSpecies) ;
        name_store("yStar{n,it}",yStar) ;
        name_store("y{n,it+1}",yIterationPlusOne) ;
        input("numSpecies{n,it},yStar{n,it}") ;
        output("y{n,it+1}") ;
        constraint("geom_cells{n,it}") ;
        constraint("speciesTransport{n,it}") ;
      }

      // Assign species mass fractions for a single cell.
      void calculate(Entity cell) { yIterationPlusOne[cell]=yStar[cell] ; }

      // Assign species mass fractions for a sequence of cells.
      void compute(const sequence &seq) {
        yIterationPlusOne.setVecSize(*numSpecies) ; do_loop(seq,this) ;
      }
  } ;

  register_rule<IterationAdvanceSpecies> registerIterationAdvanceSpecies ;

  // Iteration collapse rule for species mass fractions. Note that we
  // have changed this rule from "y{n+1}=y{n,it}" to the current form. With
  // the old form we could not do one iteration per time step, because y{n,it}
  // would not get updated and thus the residuals would never change. NOTE:
  // The above comment is now invalid, JW 12/7/2006.
  class IterationCollapseSpecies : public pointwise_rule {
    private:
      const_param<bool> iterationFinished ;
      storeVec<real> y ;
    public:

      // Define input and output.
      IterationCollapseSpecies() {
        name_store("iterationFinished{n,it-1}",iterationFinished) ;
        name_store("y{n,it}",y) ;
        input("iterationFinished{n,it-1}") ;
        input("y{n,it}") ;
        output("y{n+1}=y{n,it}") ;
        conditional("iterationFinished{n,it-1}") ;
        constraint("speciesTransport{n,it},geom_cells{n,it}") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<IterationCollapseSpecies> registerIterationCollapseSpecies ;

//-----------------------------------------------------------------------------
// Rules for computing the finite-rate chemistry species production rates.

  // Compute the species production rates for use in the "source splitting"
  // method.
  class SpeciesProductionRates : public pointwise_rule {
    private:
      const_param<int> operatorSplittingEOS ;
      const_param<real> dt ;
      const_param<real> maxTemperatureChange ;
      const_param<real> cvodeTolerance ;
      const_param<double> Tf ;
      const_param<EOS> eos ;
      const_param<reaction> reactor ;
      const_param<int> numSpecies ;
      const_store<real> vol ;
      const_store<real> rho,p,temperature,h ;
      const_store<vect3d> v ;
      const_storeVec<real> y ;
      const_storeVec<float> hint ;
      storeVec<real> yProductionRate ;
      store<real> timeStepFactor ;
#ifdef ODE_PERFORMANCE
      store<float> odeIntegrationTime ;
#endif
    private:
      int flag ;
      realtype relativeTolerance ;
      N_Vector absoluteTolerance,rhoNew,rhoOld ;
      realtype time ;
      void *cvodeMemory ;
    public:

      // Define input and output.
      SpeciesProductionRates() {
        name_store("operatorSplittingEOS",operatorSplittingEOS) ;
        name_store("dt",dt) ;
        name_store("maxTemperatureChange",maxTemperatureChange) ;
        name_store("cvodeTolerance",cvodeTolerance) ;
        name_store("Tf",Tf) ;
        name_store("eos",eos) ;
        name_store("reactor",reactor) ;
        name_store("numSpecies",numSpecies) ;
        name_store("vol",vol) ;
        name_store("rho",rho) ;
        name_store("p",p) ;
        name_store("temperature",temperature) ;
        name_store("h",h) ;
        name_store("v",v) ;
        name_store("y",y) ;
        name_store("hint_n",hint) ;
        name_store("yProductionRate",yProductionRate) ;
        name_store("source::timeStepFactor",timeStepFactor) ;
#ifdef ODE_PERFORMANCE
        name_store("odeIntegrationTime",odeIntegrationTime) ;
#endif
        input("operatorSplittingEOS") ;
        input("dt,cvodeTolerance,Tf,eos,reactor,numSpecies,vol") ;
        input("rho,p,temperature,h,v,y,hint_n,maxTemperatureChange") ;
        output("yProductionRate,source::timeStepFactor") ;
#ifdef ODE_PERFORMANCE
        output("odeIntegrationTime") ;
#endif
        constraint("combustion,sourceSplitting,geom_cells") ;
        constraint("speciesProductionCVODE") ;
        enable_dynamic_scheduling() ;
      }

      // Compute the species production rates for a single cell.
      void calculate(Entity cell) {

        // Get the starting time.
#ifdef ODE_PERFORMANCE
        double startTime=MPI_Wtime() ;
#endif

        // HARDCODE
        real localTimeStep=*dt ;
//      vect3d loc=cellCenter[cell] ;
//      if(loc.x>0.0 && loc.x<1.0 && loc.y>0.0 && loc.y<1.0)
//        localTimeStep=1.0e-03*(*dt) ;

        // New logic to improve robustness. A number of tries will be made,
        // each time reducing the relative tolerance if the previow try
        // failed. In the CVODE documentation, they recommend not using a
        // relative tolerance higher than 1.0e-03 .
        realtype localRelativeTolerance=relativeTolerance ;
        int numTries=0,maxNumTries=5 ;
        while(numTries<maxNumTries){

          // Initialize the old and new species densities. Set the absolute
          // error tolerance based on the old species mass fractions.
          for(int i=0;i<*numSpecies;++i){
            NV_Ith_S(rhoOld,i)=rho[cell]*y[cell][i] ;
            NV_Ith_S(rhoNew,i)=NV_Ith_S(rhoOld,i) ;
            NV_Ith_S(absoluteTolerance,i)=(y[cell][i]==0.0)? 1.0e-20:0.0 ;
          }

          // Turn reactions off below a cutoff temperature. Below the freezing
          // temperature we set the production rates to zero. We must also set
          // the time-step factor here.
          if(temperature[cell]<*Tf){
            for(int i=0;i<*numSpecies;++i) yProductionRate[cell][i]=0.0 ;
            timeStepFactor[cell]=1.0 ;
#ifdef ODE_PERFORMANCE
            odeIntegrationTime[cell]=0.0 ;
#endif
            return ;
          }

          // Allocate and assign the reaction data to be used by CVODE.
          real staticEnthalpy=h[cell]-0.5*dot(v[cell],v[cell]) ;
          bool invertEOS=(*operatorSplittingEOS==3)? true:false ;
          CVODEReactionData f_data(*numSpecies,*eos,rho[cell],p[cell],
            temperature[cell],staticEnthalpy,*reactor,invertEOS) ;

          // Set up CVODE.
          if(cvodeMemory==NULL){

            // Create solver memory and request BDF with Newton iteration.
            cvodeMemory=CVodeCreate(CV_BDF,CV_NEWTON) ;
            if(cvodeMemory==NULL){
              cerr << "ERROR: CVodeCreate() failed!" << endl ; Loci::Abort() ;
            }

            // Allocate CVODE memory.
            flag=CVodeMalloc(cvodeMemory,cvodeRHS,0.0,rhoOld,CV_SV,
              localRelativeTolerance,absoluteTolerance) ;
            if(flag!=CV_SUCCESS){
              cerr << "ERROR: CVodeMalloc() failed!" << endl ; Loci::Abort() ;
            }
          }else{

            // Reinitialize CVODE solver.
            flag=CVodeReInit(cvodeMemory,cvodeRHS,0.0,rhoOld,CV_SV,
              localRelativeTolerance,absoluteTolerance) ;
            if(flag!=CV_SUCCESS){
              cerr << "ERROR: CVodeReInit() failed!" << endl ; Loci::Abort() ;
            }
          }

          // Set data to be used by RHS function
          CVodeSetFdata(cvodeMemory,&f_data) ;

          // Choose the linear solver.
          flag=CVDense(cvodeMemory,*numSpecies) ;
          if(flag!=CVDENSE_SUCCESS){
            cerr << "ERROR: CVDense() failed!" << endl ; Loci::Abort() ;
          }

          // No EOS inversion.
          if(*operatorSplittingEOS==1){

            // Put the density, temperature and pressure into eosState, which
            // is used by Ed's compute_w method in cvodeRHS.
            f_data.eosState.rho=rho[cell] ; f_data.eosState.P=p[cell] ;
            f_data.eosState.T=temperature[cell] ;
                                                                                
            // Compute the reaction rates.
            f_data.reactor.extract_rates(f_data.reactionRate,f_data.eos,
              f_data.T) ;
          }

          // Call CVODE to integrate the system of ODEs to get the new species
          // densites (rhoM*y[i]). Regarding the maximum temperature change,
          // we only limit increases in temperature.
          time=0.0 ;
          if(*operatorSplittingEOS==1){
            flag=CVode(cvodeMemory,localTimeStep,rhoNew,&time,CV_NORMAL) ;
          }else{
            bool maxTemperatureChangeReached=false ;
            while(time<localTimeStep && !maxTemperatureChangeReached){

              // Invert EOS once every internal CVODE timestep.
              if(*operatorSplittingEOS==2){

                // Set the hint values from current temperature and pressure.
                f_data.hint[0]=f_data.T ; f_data.hint[1]=f_data.P ;

                // Copy species densities so we can pass to EOS inversion
                // function.
                tmp_array<real> rhoTemp(*numSpecies) ;
                for(int i=0;i<*numSpecies;++i) rhoTemp[i]=NV_Ith_S(rhoNew,i) ;

                // Solve the EOS given rho, h and y to get T and P.
                f_data.eosState=f_data.eos.State_from_rho_h(rhoTemp,
//                staticEnthalpy,f_data.hint) ;
                  staticEnthalpy,hint[cell]) ;

                // Set the new pressure and temperature.
                f_data.P=f_data.eosState.pressure() ;
                f_data.T=f_data.eosState.temperature() ;

                // Compute the reaction rates based on the new temperature and
                // pressure.
                f_data.reactor.extract_rates(f_data.reactionRate,f_data.eos,
                  f_data.T) ;
              }

              flag=CVode(cvodeMemory,localTimeStep,rhoNew,&time,CV_ONE_STEP) ;
              if(flag<0) break ;
              if((*maxTemperatureChange)>0.0){
                maxTemperatureChangeReached=((f_data.T-temperature[cell])/
                  temperature[cell]>(*maxTemperatureChange)/100.0)? true:false ;
              }else{
                maxTemperatureChangeReached=false ;
              }
//          cout << "P,T: " << f_data.P << " " << f_data.T << endl ;
            }
          }

          if(numTries==maxNumTries-1 && flag!=CV_SUCCESS){
            if(flag==CV_TOO_MUCH_WORK){
              cerr << "ERROR: CVODE could not reach end of itegration time!"
                << endl ;
              cerr << "integrationTime: " << localTimeStep << endl ;
            }else if(flag==CV_TOO_MUCH_WORK){
              cerr << "ERROR: CVODE could not satisfy itegration accuracy!"
                << endl ;
            }else if(flag==CV_ERR_FAILURE){
              cerr << "ERROR: Too many error test failures in CVODE!" << endl ;
            }else if(flag==CV_CONV_FAILURE){
              cerr << "ERROR: Too many convergence test failures in CVODE!"
                << endl ;
            }

            // Get output information from the solver.
            long int numSolverSteps,numFCall,numNewtonIterations ;
            realtype lastStepSize,currentTime,toleranceScaleFactor ;
            CVodeGetNumSteps(cvodeMemory,&numSolverSteps) ;
            CVodeGetNumRhsEvals(cvodeMemory,&numFCall) ;
            CVodeGetNumNonlinSolvIters(cvodeMemory,&numNewtonIterations) ;
            CVodeGetLastStep(cvodeMemory,&lastStepSize) ;
            CVodeGetCurrentTime(cvodeMemory,&currentTime) ;
            CVodeGetTolScaleFactor(cvodeMemory,&toleranceScaleFactor) ;
            real integrationTime=(time<localTimeStep)? time:localTimeStep ;
            cerr << "integrationTime: " << integrationTime << endl ;
            cerr << "Number of internal solver steps: " << numSolverSteps
              << endl ;
            cerr << "Number of call to f: " << numFCall << endl ;
            cerr << "Number of Newton iterations: " << numNewtonIterations
              << endl ;
            cerr << "Last internal step size: " << lastStepSize << endl ;
            cerr << "Current time reached: " << currentTime << endl ;
            cerr << "Suggested relative tolerance scaling factor: "
              << toleranceScaleFactor << endl ;
            cerr << "cell,rho,p,t,species: " << cell << " " << rho[cell] << " "
              << p[cell] << " " << temperature[cell] << " " ;
            for(int i=0;i<*numSpecies;++i) cerr << y[cell][i] << " " ;
            cerr << endl ; Loci::Abort() ;
          }else if(flag!=CV_SUCCESS){
            localRelativeTolerance/=10.0 ;
          }else{
            real integrationTime=(time<localTimeStep)? time:localTimeStep ;
            real factor=vol[cell]/integrationTime ;
            for(int i=0;i<*numSpecies;++i) yProductionRate[cell][i]=factor*
              (NV_Ith_S(rhoNew,i)-NV_Ith_S(rhoOld,i)) ;

            // If we are using the temperature cut-off, set the time-step
            // factor as well for the fluids integration.
            timeStepFactor[cell]=(time<localTimeStep)? time/(*dt):
              localTimeStep/(*dt) ;

#ifdef ODE_PERFORMANCE
            // Get the ending time and assign total computation time in seconds.
            double endTime=MPI_Wtime() ;
            odeIntegrationTime[cell]=endTime-startTime ;
#endif

            return ;
          }
          ++numTries ;
        }

#ifdef ODE_PERFORMANCE
        // Get the ending time and assign total computation time in seconds.
        double endTime=MPI_Wtime() ;
        odeIntegrationTime[cell]=endTime-startTime ;
#endif
      }

      // Assign species mass fractions for a sequence of cells.
      void compute(const sequence &seq) {

        // Set the relative tolerance.
        relativeTolerance=(*cvodeTolerance) ;

        // Must set the size of the production rate array.
        yProductionRate.setVecSize(*numSpecies) ;

        // Initialize the variables required by CVODE.
        absoluteTolerance=rhoNew=rhoOld=NULL ; cvodeMemory=NULL ;
        absoluteTolerance=N_VNew_Serial(*numSpecies) ;
        if(absoluteTolerance==NULL){
          cerr << "ERROR: Inside SpeciesProductionRates()" << endl ;
          cerr << "Problem allocating memory for absoluteTolerance" << endl ;
          Loci::Abort() ;
        }
        rhoNew=N_VNew_Serial(*numSpecies) ;
        if(rhoNew==NULL){
          cerr << "ERROR: Inside SpeciesProductionRates()" << endl ;
          cerr << "Problem allocating memory for rhoNew" << endl ;
          Loci::Abort() ;
        }
        rhoOld=N_VNew_Serial(*numSpecies) ;
        if(rhoOld==NULL){
          cerr << "ERROR: Inside SpeciesProductionRates()" << endl ;
          cerr << "Problem allocating memory for rhoOld" << endl ;
          Loci::Abort() ;
        }

        // Compute the production rates.
        do_loop(seq,this) ;

        // Deallocate memory associated with CVODE.
        N_VDestroy_Serial(absoluteTolerance) ; N_VDestroy_Serial(rhoNew) ;
        N_VDestroy_Serial(rhoOld) ;
        if(cvodeMemory) CVodeFree(&cvodeMemory) ;
      }
  } ;

  register_rule<SpeciesProductionRates> registerSpeciesProductionRates ;

  // Call the mechanism library to compute the instantaneous species production
  // rates.
  class InstantaneousSpeciesProductionRates : public pointwise_rule {
    private:
      const_param<int> numSpecies ;
      const_param<EOS> eos ;
      const_param<reaction> reactor ;
      const_param<double> Tf ;
      const_store<real> vol ;
      const_store<real> rho,p,T ;
      const_storeVec<real> y ;
      storeVec<real> yProductionRate ;
    private:
      reaction::rates *reactionRate ;
    public:

      // Define input and output.
      InstantaneousSpeciesProductionRates() {
        name_store("numSpecies",numSpecies) ;
        name_store("eos",eos) ;
        name_store("reactor",reactor) ;
        name_store("Tf",Tf) ;
        name_store("vol",vol) ;
        name_store("rho",rho) ;
        name_store("p",p) ;
        name_store("temperature",T) ;
        name_store("y",y) ;
        name_store("yProductionRate",yProductionRate) ;
        input("numSpecies,eos,reactor,Tf,vol,rho,p,temperature,y") ;
        output("yProductionRate") ;
        constraint("geom_cells,speciesProductionInstantaneous") ;
      }

      // Compute the species production rates for a single cell.
      void calculate(Entity cell) {

        // Copy mass fractions.
        tmp_array<real> yNew(*numSpecies) ;
        for(int i=0;i<*numSpecies;++i) yNew[i]=y[cell][i] ;

        // Pre-compute the reaction rates. Below the freezing temperature we
        // set the species production rates to zero and return.
        if(T[cell]>*Tf){
          reactor->extract_rates(reactionRate,*eos,T[cell]) ;
        }else{
          for(int i=0;i<*numSpecies;++i) yProductionRate[cell][i]=0.0 ; return ;
        }

        // Create an eos_state so we can just use Ed's original compute_w()
        // which only uses this to get rho, P and T.
        EOS::State eosState ;
        eosState.rho=rho[cell] ; eosState.P=p[cell] ; eosState.T=T[cell] ;

        // Call the combustion library to get the species production rates.
        tmp_array<real> w(*numSpecies) ;
//      reactor->compute_w(w,reactionRate,yNew,rho[cell],p[cell],T[cell],*eos) ;
        reactor->compute_w(w,reactionRate,yNew,eosState,*eos) ;
        for(int i=0;i<*numSpecies;++i){
          yProductionRate[cell][i]=w[i]*vol[cell] ;
        }
      }

      // Assign species mass fractions for a sequence of cells.
      void compute(const sequence &seq) {

        // Must set the size of the production rate array.
        yProductionRate.setVecSize(*numSpecies) ;
                                                                                
        // Set up memory for pre-computed reaction rates.
        reactionRate=new reaction::rates[reactor->num_rates()] ;

        // Compute the production rates.
        do_loop(seq,this) ;

        // Delete other memory.
        delete [] reactionRate ;
      }
  } ;

  register_rule<InstantaneousSpeciesProductionRates>
    registerInstantaneousSpeciesProductionRates ;

  // Performs the ODE step for the Strang splitting algorithm.
  class ODEIntegration : public pointwise_rule {
    private:
      const_param<int> operatorSplittingEOS ;
      const_param<real> maxTemperatureChange ;
      const_param<real> dt ;
      const_param<real> cvodeTolerance ;
      const_param<double> Tf ;
      const_param<EOS> eos ;
      const_param<reaction> reactor ;
      const_param<int> numSpecies ;
      const_store<real> vol ;
      const_store<real> rho,p,temperature,h ;
      const_store<vect3d> v ;
      const_storeVec<real> y ;
      const_storeVec<float> hint ;
      storeVec<real> yODE ;
      store<real> timeStepFactor ;
#ifdef ODE_PERFORMANCE
      store<float> odeIntegrationTime ;
#endif
    private:
      int flag ;
      realtype relativeTolerance ;
      N_Vector absoluteTolerance,rhoNew,rhoOld ;
      realtype time ;
      void *cvodeMemory ;
      reaction::rates *reactionRate ;
    public:

      // Define input and output.
      ODEIntegration() {
        name_store("operatorSplittingEOS",operatorSplittingEOS) ;
        name_store("maxTemperatureChange",maxTemperatureChange) ;
        name_store("dt",dt) ;
        name_store("cvodeTolerance",cvodeTolerance) ;
        name_store("Tf",Tf) ;
        name_store("eos",eos) ;
        name_store("reactor",reactor) ;
        name_store("numSpecies",numSpecies) ;
        name_store("vol",vol) ;
        name_store("rho",rho) ;
        name_store("p",p) ;
        name_store("temperature",temperature) ;
        name_store("h",h) ;
        name_store("v",v) ;
        name_store("y",y) ;
        name_store("hint_n",hint) ;
        name_store("yODE",yODE) ;
        name_store("strang::timeStepFactor",timeStepFactor) ;
#ifdef ODE_PERFORMANCE
        name_store("odeIntegrationTime",odeIntegrationTime) ;
#endif
        input("operatorSplittingEOS") ;
        input("maxTemperatureChange,dt,cvodeTolerance,Tf,eos,reactor") ;
        input("numSpecies,vol,rho,p,temperature,h,v,y,hint_n") ;
        output("yODE,strang::timeStepFactor") ;
#ifdef ODE_PERFORMANCE
        output("odeIntegrationTime") ;
#endif
        constraint("combustion,strangSplitting,geom_cells") ;
        enable_dynamic_scheduling() ;
      }

      // Integrate the ODE for a single cell.
      void calculate(Entity cell) {

        // Get the starting time.
#ifdef ODE_PERFORMANCE
        double startTime=MPI_Wtime() ;
#endif
        // Below the freezing temperature we set the outgoing species mass
        // fraction to what came in. Must also time-step factor here.
        if(temperature[cell]<*Tf){
          for(int i=0;i<*numSpecies;++i) yODE[cell][i]=y[cell][i] ;
          timeStepFactor[cell]=1.0 ;
#ifdef ODE_PERFORMANCE
          odeIntegrationTime[cell]=0.0 ;
#endif
          return ;
        }

        // New logic to improve robustness. A number of tries will be made,
        // each time reducing the relative tolerance if the previos try
        // failed. In the CVODE documentation, they recommend not using a
        // relative tolerance higher than 1.0e-03 .
        realtype localRelativeTolerance=relativeTolerance ;
        int numTries=0,maxNumTries=5 ;
        while(numTries<maxNumTries){

          // Initialize the old and new species densities. Set the absolute
          // error tolerance based on the old species mass fractions.
          for(int i=0;i<*numSpecies;++i){
            NV_Ith_S(rhoOld,i)=rho[cell]*y[cell][i] ;
            NV_Ith_S(rhoNew,i)=NV_Ith_S(rhoOld,i) ;
            NV_Ith_S(absoluteTolerance,i)=(y[cell][i]==0.0)? 1.0e-20:0.0 ;
          }

          // Allocate and assign the reaction data to be used by CVODE.
          real staticEnthalpy=h[cell]-0.5*dot(v[cell],v[cell]) ;
          //real staticEnergy=staticEnthalpy-p[cell]/rho[cell] ;
          bool invertEOS=(*operatorSplittingEOS==3)? true:false ;
          CVODEReactionData f_data(*numSpecies,*eos,rho[cell],p[cell],
            temperature[cell],staticEnthalpy,*reactor,invertEOS) ;

          // Set up CVODE.
          if(cvodeMemory==NULL){

            // Create solver memory and request BDF with Newton iteration.
            cvodeMemory=CVodeCreate(CV_BDF,CV_NEWTON) ;
            if(cvodeMemory==NULL){
              cerr << "ERROR: CVodeCreate() failed!" << endl ; Loci::Abort() ;
            }

            // Allocate CVODE memory.
            flag=CVodeMalloc(cvodeMemory,cvodeRHS,0.0,rhoOld,CV_SV,
              localRelativeTolerance,absoluteTolerance) ;
            if(flag!=CV_SUCCESS){
              cerr << "ERROR: CVodeMalloc() failed!" << endl ; Loci::Abort() ;
            }
          }else{

            // Reinitialize CVODE solver.
            flag=CVodeReInit(cvodeMemory,cvodeRHS,0.0,rhoOld,CV_SV,
              localRelativeTolerance,absoluteTolerance) ;
            if(flag!=CV_SUCCESS){
              cerr << "ERROR: CVodeReInit() failed!" << endl ; Loci::Abort() ;
            }
          }

          // Set the maximum number of steps.
          CVodeSetMaxNumSteps(cvodeMemory,100000) ;

          // Set data to be used by RHS function
          CVodeSetFdata(cvodeMemory,&f_data) ;

          // Choose the linear solver.
          flag=CVDense(cvodeMemory,*numSpecies) ;
          if(flag!=CVDENSE_SUCCESS){
            cerr << "ERROR: CVDense() failed!" << endl ; Loci::Abort() ;
          }

          // No EOS inversion during solution of ODE.
          if(*operatorSplittingEOS==1){

            // Put the density, temperature and pressure into eosState, which
            // is used by Ed's compute_w method in cvodeRHS.
            f_data.eosState.rho=rho[cell] ; f_data.eosState.P=p[cell] ;
            f_data.eosState.T=temperature[cell] ;
                                                                                
            // Compute the reaction rates.
            f_data.reactor.extract_rates(f_data.reactionRate,f_data.eos,
              f_data.T) ;
          }

          // Call CVODE to integrate the system of ODEs to get the new species
          // densites (rhoM*y[i]). Regarding the maximum temperature change,
          // we only limit increases in temperature.
          time=0.0 ;
          if(*operatorSplittingEOS==1){
            flag=CVode(cvodeMemory,*dt,rhoNew,&time,CV_NORMAL) ;
          }else{
            bool maxTemperatureChangeReached=false ;
            while(time<(*dt) && !maxTemperatureChangeReached){

              // Invert EOS once every internal CVODE timestep.
              if(*operatorSplittingEOS==2){

                // Set the hint values from current temperature and pressure.
                f_data.hint[0]=f_data.T ; f_data.hint[1]=f_data.P ;

                // Copy species densities so we can pass to EOS inversion
                // function.
                tmp_array<real> rhoTemp(*numSpecies) ;
                for(int i=0;i<*numSpecies;++i) rhoTemp[i]=NV_Ith_S(rhoNew,i) ;
                                                                                
                // Solve the EOS given rho, h and y to get T and P.
                f_data.eosState=f_data.eos.State_from_rho_h(rhoTemp,
                  staticEnthalpy,hint[cell]) ;
//              f_data.eosState=f_data.eos.State_from_rho_e(rhoTemp,
//                staticEnergy,hint[cell]) ;

                // Set the new pressure and temperature.
                f_data.P=f_data.eosState.pressure() ;
                f_data.T=f_data.eosState.temperature() ;
                                                                                
                // Compute the reaction rates based on the new temperature and
                // pressure.
                f_data.reactor.extract_rates(f_data.reactionRate,f_data.eos,
                  f_data.T) ;
              }

              flag=CVode(cvodeMemory,*dt,rhoNew,&time,CV_ONE_STEP) ;
              if(flag<0) break ;

              // Limit species densities that come out of ODE solver. Since the
              // ODE is integrated at constant density, there should be no change
              // in total mixture density during the ODE solve.
              tmp_array<real> yTemp(*numSpecies) ;
              for(int i=0;i<*numSpecies;++i) yTemp[i]=NV_Ith_S(rhoNew,i)/rho[cell] ;
              real ySum=0.0 ;
              for(int i=0;i<*numSpecies;++i){
                if(yTemp[i]<1.0e-30) yTemp[i]=0.0 ;
                else if(yTemp[i]>1.0) yTemp[i]=1.0 ;
                ySum+=yTemp[i] ;
              }
              for(int i=0;i<*numSpecies;++i){
                yTemp[i]/=ySum ; NV_Ith_S(rhoNew,i)=rho[cell]*yTemp[i] ;
              }

              if((*maxTemperatureChange)>0.0){
                maxTemperatureChangeReached=((f_data.T-temperature[cell])/
                  temperature[cell]>(*maxTemperatureChange)/100.0)? true:false ;
              }else{
                maxTemperatureChangeReached=false ;
              }
//            cout << "P,T: " << f_data.P << " " << f_data.T << endl ;
            }
          }

          if(numTries==maxNumTries-1 && flag!=CV_SUCCESS){
            if(flag==CV_TOO_MUCH_WORK){
              cerr << "ERROR: CVODE could not reach end of itegration time!"
                << endl ;
              cerr << "integrationTime: " << *dt << endl ;
            }else if(flag==CV_TOO_MUCH_WORK){
              cerr << "ERROR: CVODE could not satisfy itegration accuracy!"
                << endl ;
            }else if(flag==CV_ERR_FAILURE){
              cerr << "ERROR: Too many error test failures in CVODE!" << endl ;
            }else if(flag==CV_CONV_FAILURE){
              cerr << "ERROR: Too many convergence test failures in CVODE!"
                << endl ;
            }

            // Get output information from the solver.
            long int numSolverSteps,numFCall,numNewtonIterations ;
            realtype lastStepSize,currentTime,toleranceScaleFactor ;
            CVodeGetNumSteps(cvodeMemory,&numSolverSteps) ;
            CVodeGetNumRhsEvals(cvodeMemory,&numFCall) ;
            CVodeGetNumNonlinSolvIters(cvodeMemory,&numNewtonIterations) ;
            CVodeGetLastStep(cvodeMemory,&lastStepSize) ;
            CVodeGetCurrentTime(cvodeMemory,&currentTime) ;
            CVodeGetTolScaleFactor(cvodeMemory,&toleranceScaleFactor) ;
//          cerr << "integrationTime: " << integrationTime[cell] << endl ;
            cerr << "Number of internal solver steps: " << numSolverSteps
              << endl ;
            cerr << "Number of call to f: " << numFCall << endl ;
            cerr << "Number of Newton iterations: " << numNewtonIterations
              << endl ;
            cerr << "Last internal step size: " << lastStepSize << endl ;
            cerr << "Current time reached: " << currentTime << endl ;
            cerr << "Suggested relative tolerance scaling factor: "
              << toleranceScaleFactor << endl ;
            cerr << "cell,rho,p,t,species: " << cell << " " << rho[cell] << " "
              << p[cell] << " " << temperature[cell] << " " ;
            for(int i=0;i<*numSpecies;++i) cerr << y[cell][i] << " " ;
            cerr << endl ; Loci::Abort() ;
          }else if(flag!=CV_SUCCESS){
            localRelativeTolerance/=10.0 ;
          }else{

            // Compute new species mass fractions. The mixture density should
            // not change from the ODE integration, so use the incoming
            // density to convert from species densities to mass fractions.
            for(int i=0;i<*numSpecies;++i) yODE[cell][i]=NV_Ith_S(rhoNew,i)/
              rho[cell] ;

            // If we are using the temperature cut-off, set the time-step
            // factor as well for the fluids integration.
            timeStepFactor[cell]=(time<(*dt))? time/(*dt):1.0 ;

#ifdef ODE_PERFORMANCE
            // Get the ending time and assign total computation time in seconds.
            double endTime=MPI_Wtime() ;
            odeIntegrationTime[cell]=endTime-startTime ;
#endif

            return ;
          }
          ++numTries ;
        }

#ifdef ODE_PERFORMANCE
        // Get the ending time and assign total computation time in seconds.
        double endTime=MPI_Wtime() ;
        odeIntegrationTime[cell]=endTime-startTime ;
#endif
      }

      // Assign species mass fractions for a sequence of cells.
      void compute(const sequence &seq) {

        // Set the relative tolerance.
        relativeTolerance=(*cvodeTolerance) ;

        // Must set the size of the production rate array.
        yODE.setVecSize(*numSpecies) ;

        // Initialize the variables required by CVODE.
        absoluteTolerance=rhoNew=rhoOld=NULL ; cvodeMemory=NULL ;
        absoluteTolerance=N_VNew_Serial(*numSpecies) ;
        if(absoluteTolerance==NULL){
          cerr << "ERROR: Inside ODEIntegration()" << endl ;
          cerr << "Problem allocating memory for absoluteTolerance" << endl ;
          Loci::Abort() ;
        }
        rhoNew=N_VNew_Serial(*numSpecies) ;
        if(rhoNew==NULL){
          cerr << "ERROR: Inside ODEIntegration()" << endl ;
          cerr << "Problem allocating memory for rhoNew" << endl ;
          Loci::Abort() ;
        }
        rhoOld=N_VNew_Serial(*numSpecies) ;
        if(rhoOld==NULL){
          cerr << "ERROR: Inside ODEIntegration()" << endl ;
          cerr << "Problem allocating memory for rhoOld" << endl ;
          Loci::Abort() ;
        }


        // Compute the production rates.
        do_loop(seq,this) ;

        // Deallocate memory associated with CVODE.
        N_VDestroy_Serial(absoluteTolerance) ; N_VDestroy_Serial(rhoNew) ;
        N_VDestroy_Serial(rhoOld) ;
        if(cvodeMemory) CVodeFree(&cvodeMemory) ;
      }
  } ;

  register_rule<ODEIntegration> registerODEIntegration ;

#ifdef ODE_PERFORMANCE
  // Unit/apply sequence to compute total cpu time for ODE integration
  // for each process.
  class TotalODEIntegrationTimeUnit : public unit_rule {
    private:
      param<vector<real> > totalODEIntegrationTime ;
    public:

      // Define input and output.
      TotalODEIntegrationTimeUnit() {
        name_store("totalODEIntegrationTime",totalODEIntegrationTime) ;
        output("totalODEIntegrationTime") ;
        constraint("UNIVERSE") ;
      }

      // Assign species mass fractions for all faces in sequence.
      virtual void compute(const sequence &seq) {
        *totalODEIntegrationTime=vector<real>(Loci::MPI_processes,0.0) ;
      }
  } ;

  register_rule<TotalODEIntegrationTimeUnit> registerTotalODEIntegrationTimeUnit ;

  // Summation join operator for standard library vector<real>.
  class StandardVectorSummation {
    public:
      void operator()(vector<real> &a,const vector<real> &b){
        vector<real>::const_iterator bPtr=b.begin() ;
        for(vector<real>::iterator aPtr=a.begin();aPtr!=a.end();++aPtr,++bPtr)
          *aPtr+=*bPtr ;
      }
  } ;

  // Apply rule for net inlet species mass flow rate.
  class TotalODEIntegrationTimeApply : public
  apply_rule<param<vector<real> >,StandardVectorSummation> {
    private:
      const_store<float> odeIntegrationTime ;
      param<vector<real> > totalODEIntegrationTime ;
    private:
      vector<real> cellIntegrationTime ;
    public:

      // Define input and output.
      TotalODEIntegrationTimeApply() {
        name_store("odeIntegrationTime",odeIntegrationTime) ;
        name_store("totalODEIntegrationTime",totalODEIntegrationTime) ;
        input("odeIntegrationTime") ;
        output("totalODEIntegrationTime") ;
        constraint("geom_cells") ;
      }

      // Add time for a single cell.
      void calculate(Entity cell) {
        cellIntegrationTime[Loci::MPI_rank]=odeIntegrationTime[cell] ;
        join(*totalODEIntegrationTime,cellIntegrationTime) ;
      }

      // Loop over cells.
      void compute(const sequence &seq) {
        cellIntegrationTime=vector<real>(Loci::MPI_processes,0.0) ;
        do_loop(seq,this) ;
      }
  } ;

  register_rule<TotalODEIntegrationTimeApply>
    registerTotalODEIntegrationTimeApply ;

  class ODEIntegratorData : public singleton_rule {
    private:
      const_param<int> n ;
      const_param<bool> debug ;
      const_param<vector<real> > totalODEIntegrationTime ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      ODEIntegratorData() {
        name_store("$n{n}",n) ;
        name_store("debug{n}",debug) ;
        name_store("totalODEIntegrationTime{n}",totalODEIntegrationTime) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        input("$n{n},debug{n},totalODEIntegrationTime{n}") ;
        output("OUTPUT{n}") ;
      }

      // Write the data.
      virtual void compute(const sequence &seq) {
        if(Loci::MPI_rank==0){
          cout << "ODE Integration Data for timestep: " << *n << endl ;
          for(size_t i=0;i<(*totalODEIntegrationTime).size();++i)
            cout << "  Proc " << i << ": " << (*totalODEIntegrationTime)[i] << endl ;
        }
      }
  } ;

  register_rule<ODEIntegratorData> registerODEIntegratorData ;
#endif
}

