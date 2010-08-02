//-----------------------------------------------------------------------------
// Description: This file contains rules for cavitation chemistry model
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
using Loci::Area ;

// StreamUns includes.
#include "const.h"
#include "referenceFrame.h"
#include "residual.h"
#include "sciTypes.h"
#include "varsFileInputs.h"
#include "cavitationVarsFileInputs.h"

//// Include the cavitation table class and create one
//#include "cavitation_table.h"
//cavitation_Table cavitation_table;

namespace streamUns {

//-----------------------------------------------------------------------------
// Rules to process cavitation equation options from the .vars file.

  // Creates the constraints for the cavitation equations.
  class CavitationModelConstraints : public constraint_rule {
    private:
      const_param<CavitationEquationOptions> cavitationEquationOptions ;
      Constraint cavitationModel ;
    public:
                                                                                
      // Define input and output.
      CavitationModelConstraints() {
        name_store("cavitationEquationOptions",cavitationEquationOptions) ;
        name_store("cavitationModel",cavitationModel) ;
        input("cavitationEquationOptions") ;
        output("cavitationModel") ;
      }
                                                                                
      // Set up the constraints.
      virtual void compute(const sequence& seq) {
	    if((*cavitationEquationOptions).optionExists("model")) {
		cavitationModel=~EMPTY ;
		}
      }
  } ;
                                                                                
  register_rule<CavitationModelConstraints>
    registerCavitationModelConstraints ;
	
  // Creates the cavitation equation solver constraints.
  class CavitationEquationSolverConstraints : public constraint_rule {
    private:
      const_param<CavitationEquationOptions> cavitationEquationOptions ;
      Constraint cavitationSGSLinearSolver ;
    public:
                                                                                
      // Define input and output.
      CavitationEquationSolverConstraints() {
        name_store("cavitationEquationOptions",cavitationEquationOptions) ;
        name_store("cavitationSGSLinearSolver",cavitationSGSLinearSolver) ;
        input("cavitationEquationOptions") ;
        output("cavitationSGSLinearSolver") ;
      }
                                                                                
      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        if((*cavitationEquationOptions).optionExists("linearSolver")){
          Loci::option_value_type optionValueType=cavitationEquationOptions->
            getOptionValueType("linearSolver") ;
          switch(optionValueType){
            case Loci::NAME:
              {
                Loci::option_values optionValues=cavitationEquationOptions->
                  getOption("linearSolver") ;
                string name ; optionValues.get_value(name) ;
                if(name=="SGS"){
                  cavitationSGSLinearSolver=~EMPTY ;
                }else{
                  cerr << "Bad linearSolver for cavitationEquation." << endl ;
                  Loci::Abort() ;
                }
              }
              break ;
            default:
              cerr << "Bad type for linearSolver in cavitationEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          cavitationSGSLinearSolver=~EMPTY ;
        }
      }
  } ;
																				
  register_rule<CavitationEquationSolverConstraints>
    registerCavitationEquationSolverConstraints ;

  // Creates the cavitation inviscid flux constraints.
  class CavitationInviscidFluxConstraints : public constraint_rule {
    private:
      const_param<string> cavitationInviscidFlux ;
      Constraint fouCavitation,souCavitation ;
    public:

      // Define input and output.
      CavitationInviscidFluxConstraints() {
        name_store("cavitationInviscidFlux",cavitationInviscidFlux) ;
        name_store("fouCavitation",fouCavitation) ;
        name_store("souCavitation",souCavitation) ;
        input("cavitationInviscidFlux") ;
        output("fouCavitation,souCavitation") ;
      }

      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        if(*cavitationInviscidFlux!="FOU" && *cavitationInviscidFlux!="SOU"){
          cerr << "Bad cavitationInviscidFlux in .vars file." << endl ;
          Loci::Abort() ;
        }
        fouCavitation=(*cavitationInviscidFlux=="FOU")? ~EMPTY:EMPTY ;
        souCavitation=(*cavitationInviscidFlux=="SOU")? ~EMPTY:EMPTY ;
      }
  } ;       

  register_rule<CavitationInviscidFluxConstraints>
    registerCavitationInviscidFluxConstraints ;
	
  // Creates the cavitation equation solver parameters.
  class CavitationEquationSolverParameters : public singleton_rule {
    private:
      const_param<CavitationEquationOptions> cavitationEquationOptions ;
      param<int> cavitationMaxIterations ;
      param<real> cavitationRelaxationFactor ;
      param<real> Cdest,Cprod,rhol,rhov,pv,vInfinite,tInfinite ;
    public:
                                                                                
      // Define input and output.
      CavitationEquationSolverParameters() {
        name_store("cavitationEquationOptions",cavitationEquationOptions) ;
        name_store("Cdest",Cdest) ;
        name_store("Cprod",Cprod) ;
        name_store("rhol",rhol) ;
        name_store("rhov",rhov) ;
        name_store("pv",pv) ;
        name_store("vInfinite",vInfinite) ;
        name_store("tInfinite",tInfinite) ;
        name_store("cavitationMaxIterations",cavitationMaxIterations) ;
        name_store("cavitationRelaxationFactor",cavitationRelaxationFactor) ;
        input("cavitationEquationOptions") ;
        output("Cdest,Cprod,rhol,rhov,pv,vInfinite,tInfinite,cavitationMaxIterations,cavitationRelaxationFactor") ;
      }
                                                                                
      // Set up the parameters.
      virtual void compute(const sequence& seq) {
                
//	// Get table name and initialize the table
//	Loci::option_values optionValues=cavitationEquationOptions->getOption("table") ;
//        string tableFileName;
//        optionValues.get_value(tableFileName) ;
//        cavitation_table.read(tableFileName);
                                                                
        // Emperical constant Cdest for cavitation model.
        if((*cavitationEquationOptions).optionExists("Cdest")){
          Loci::option_value_type optionValueType=cavitationEquationOptions->
            getOptionValueType("Cdest") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                cavitationEquationOptions->getOption("Cdest",temp) ;
                if(int(temp)<0){
                  cerr << "Bad Cdest value for cavitationEquation."
                    << endl ; Loci::Abort() ;
                }
                *Cdest=temp ;
              }
              break ;
            default:
              cerr << "Bad type for Cdest in cavitationEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *Cdest=0.0 ;
        }

        // Emperical constant Cprod for cavitation model.
        if((*cavitationEquationOptions).optionExists("Cprod")){
          Loci::option_value_type optionValueType=cavitationEquationOptions->
            getOptionValueType("Cprod") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                cavitationEquationOptions->getOption("Cprod",temp) ;
                if(int(temp)<0){
                  cerr << "Bad Cprod value for cavitationEquation."
                    << endl ; Loci::Abort() ;
                }
                *Cprod=temp ;
              }
              break ;
            default:
              cerr << "Bad type for Cprod in cavitationEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *Cprod=0.0 ;
        }

        // Liquid density at given temperature.
        if((*cavitationEquationOptions).optionExists("rhol")){
          Loci::option_value_type optionValueType=cavitationEquationOptions->
            getOptionValueType("rhol") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                cavitationEquationOptions->getOption("rhol",temp) ;
                if(int(temp)<0){
                  cerr << "Bad rhol value for cavitationEquation."
                    << endl ; Loci::Abort() ;
                }
                *rhol=temp ;
              }
              break ;
            default:
              cerr << "Bad type for rhol in cavitationEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *rhol=1.0 ;
        }

        // Vapor density at given temperature.
        if((*cavitationEquationOptions).optionExists("rhov")){
          Loci::option_value_type optionValueType=cavitationEquationOptions->
            getOptionValueType("rhov") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                cavitationEquationOptions->getOption("rhov",temp) ;
                if(int(temp)<0){
                  cerr << "Bad rhov value for cavitationEquation."
                    << endl ; Loci::Abort() ;
                }
                *rhov=temp ;
              }
              break ;
            default:
              cerr << "Bad type for rhov in cavitationEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *rhov=0.001 ;
        }

        // Vapor pressure at given temperature.
        if((*cavitationEquationOptions).optionExists("pv")){
          Loci::option_value_type optionValueType=cavitationEquationOptions->
            getOptionValueType("pv") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                cavitationEquationOptions->getOption("pv",temp) ;
                //for dimensionless form
                if(int(temp)>0){
                  cerr << "Bad pv value for cavitationEquation."
                    << endl ; Loci::Abort() ;
                }
                *pv=temp ;
              }
              break ;
            default:
              cerr << "Bad type for pv in cavitationEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *pv=-0.2 ;
        }

        // Referenced velocity scale.
        if((*cavitationEquationOptions).optionExists("vInfinite")){
          Loci::option_value_type optionValueType=cavitationEquationOptions->
            getOptionValueType("vInfinite") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                cavitationEquationOptions->getOption("vInfinite",temp) ;
                if(int(temp)<0){
                  cerr << "Bad vInfinite value for cavitationEquation."
                    << endl ; Loci::Abort() ;
                }
                *vInfinite=temp ;
              }
              break ;
            default:
              cerr << "Bad type for vInfinite in cavitationEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *vInfinite=1.0 ;
        }

        // Referenced time scale.
        if((*cavitationEquationOptions).optionExists("tInfinite")){
          Loci::option_value_type optionValueType=cavitationEquationOptions->
            getOptionValueType("tInfinite") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                cavitationEquationOptions->getOption("tInfinite",temp) ;
                if(int(temp)<0){
                  cerr << "Bad tInfinite value for cavitationEquation."
                    << endl ; Loci::Abort() ;
                }
                *tInfinite=temp ;
              }
              break ;
            default:
              cerr << "Bad type for tInfinite in cavitationEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *tInfinite=1.0 ;
        }
		
        // Relaxation factor.
        if((*cavitationEquationOptions).optionExists("relaxationFactor")){
          Loci::option_value_type optionValueType=cavitationEquationOptions->
            getOptionValueType("relaxationFactor") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                cavitationEquationOptions->getOption("relaxationFactor",temp) ;
                if(temp<=0.0 || temp>1.0){
                  cerr << "Bad relaxationFactor for cavitationEquation."
                    << endl ; Loci::Abort() ;
                }
                *cavitationRelaxationFactor=temp ;
              }
              break ;
            default:
              cerr << "Bad type for relaxationFactor in cavitationEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *cavitationRelaxationFactor=0.5 ;
        }
                                                                                
        // Maximum number of iterations.
        if((*cavitationEquationOptions).optionExists("maxIterations")){
          Loci::option_value_type optionValueType=cavitationEquationOptions->
            getOptionValueType("maxIterations") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                cavitationEquationOptions->getOption("maxIterations",temp) ;
                if(int(temp)<0){
                  cerr << "Bad maxIterations value for cavitationEquation."
                    << endl ; Loci::Abort() ;
                }
                *cavitationMaxIterations=int(temp) ;
              }
              break ;
            default:
              cerr << "Bad type for maxIterations in cavitationEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *cavitationMaxIterations=5 ;
        }
      }
  } ;
                                                                                
  register_rule<CavitationEquationSolverParameters>
    registerCavitationEquationSolverParamters ;

//-----------------------------------------------------------------------------
// Boundary condition rules for Alpha.

  // Rule for boundary faces with specified Alpha. Assigns value to all boundary
  // faces that have the property Alpha_BC.
  class BoundaryAlphaSpecification : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<real> Alpha_BC ;
      store<real> Alpha_f ;
    public:

      // Define input and output.
      BoundaryAlphaSpecification() {
        name_store("ref",ref) ;
        name_store("Alpha_BC",Alpha_BC) ;
        name_store("Alpha_f",Alpha_f) ;
        input("ref->Alpha_BC") ;
        output("Alpha_f") ;
      }

      // Calculate Alpha for a single face.
      void calculate(Entity face) { Alpha_f[face]=Alpha_BC[ref[face]] ; }

      // Calculate Alpha for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryAlphaSpecification>
    registerBoundaryAlphaSpecification ;

  // Rule for extrapolating Alpha to boundary faces. This occurs for all outlets,
  // noslip, slip and symmetry boundaries. Right now we are using the low-order method.
  class BoundaryAlphaExtrapolation : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> Alpha ;
      store<real> Alpha_f ;
    public:

      // Define input and output.
      BoundaryAlphaExtrapolation() {
        name_store("ci",ci) ;
        name_store("Alpha",Alpha) ;
        name_store("Alpha_f",Alpha_f) ;
        input("ci->Alpha") ;
        output("Alpha_f") ;
        constraint("extrapolatedAlpha_BC") ;
      }

      // Calculate Alpha for a single face.
      void calculate(Entity face) {
        Alpha_f[face]=Alpha[ci[face]] ;
      }

      // Calculate Alpha for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryAlphaExtrapolation> registerBoundaryAlphaExtrapolation ;

//-----------------------------------------------------------------------------
// Rules for assembling the Alpha (liquid volume fraction) equation.

  // Rule to initialize the main coefficient.
  class initializeAlphaMainCoefficient : public unit_rule {
    private:
      store<real> AlphaMainCoefficient ;
    public:

      // Define input and output.
      initializeAlphaMainCoefficient() {
        name_store("AlphaMainCoefficient",AlphaMainCoefficient) ;
        output("AlphaMainCoefficient") ;
//        constraint("vol") ;
        constraint("geom_cells") ;
      }

      // Set the main coefficient to Zero for a single cell.
      void calculate(Entity cell) { AlphaMainCoefficient[cell]=0.0 ; }

      // Set the main coefficient to Zero for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<initializeAlphaMainCoefficient> registerinitializeAlphaMainCoefficient ;

  // Rule to add the first-order inviscid flux contribution to the main
  // coefficient for interior faces.
  class FOUInviscidFluxToAlphaMainCoefficientInterior : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> massFlux ;
      store<real> AlphaMainCoefficient ;
    public:

      // Define input and output.
      FOUInviscidFluxToAlphaMainCoefficientInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("massFlux",massFlux) ;
        name_store("AlphaMainCoefficient",AlphaMainCoefficient) ;
        input("massFlux") ;
        output("cl->AlphaMainCoefficient,cr->AlphaMainCoefficient") ;
        constraint("internalFaces") ;
      }

      // Increment the main coefficient for cells attached to a single face.
      void calculate(Entity face) {
        if(massFlux[face]>0.0){
          AlphaMainCoefficient[cr[face]]+=massFlux[face] ;
        }else{
          AlphaMainCoefficient[cl[face]]-=massFlux[face] ;
        }
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToAlphaMainCoefficientInterior>
    registerFOUInviscidFluxToAlphaMainCoefficientInterior ;

  // Rule to add the first-order inviscid flux contribution to the main
  // coefficient for boundary faces.
  class FOUInviscidFluxToAlphaMainCoefficientBoundary : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> massFlux ;
      store<real> AlphaMainCoefficient ;
    public:

      // Define input and output.
      FOUInviscidFluxToAlphaMainCoefficientBoundary() {
        name_store("ci",ci) ;
        name_store("massFlux",massFlux) ;
        name_store("AlphaMainCoefficient",AlphaMainCoefficient) ;
        input("massFlux") ;
        output("ci->AlphaMainCoefficient") ;
        constraint("boundaryFaces") ;
      }

      // Increment the main coefficient for the cell attach to a single face.
      void calculate(Entity face) {
        if(massFlux[face]<0.0) AlphaMainCoefficient[ci[face]]-=massFlux[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToAlphaMainCoefficientBoundary>
    registerFOUInviscidFluxToAlphaMainCoefficientBoundary ;

	
  // Rule to initialize the source term.
  class initializeAlphaSourceTerm : public unit_rule {
	  private:
		store<real> AlphaSourceTerm ;
	  public:

      // Define input and output.
	initializeAlphaSourceTerm() {
		name_store("AlphaSourceTerm",AlphaSourceTerm) ;
		output("AlphaSourceTerm") ;
//		constraint("vol") ;
        constraint("geom_cells") ;
	}

      // Set the source term to Zero for a single cell.
	void calculate(Entity cell) { AlphaSourceTerm[cell]=0.0 ; }

      // Set the source term to Zero for a sequence of cells.
	virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<initializeAlphaSourceTerm> registerinitializeAlphaSourceTerm ;
  
  
  // Rule to add the first-order inviscid flux contribution to the source term
  // for boundary faces.
  class FOUInviscidFluxToAlphaSourceTermBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> Alpha ;
      const_store<real> massFlux ;
      const_store<real> Alpha_f ;
      store<real> AlphaSourceTerm ;
    public:

      // Define input and output.
      FOUInviscidFluxToAlphaSourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("Alpha",Alpha) ;
        name_store("massFlux",massFlux) ;
        name_store("Alpha_f",Alpha_f) ;
        name_store("AlphaSourceTerm",AlphaSourceTerm) ;
        input("ci->Alpha,massFlux,Alpha_f") ;
        output("ci->AlphaSourceTerm") ;
        constraint("boundaryFaces,fouCavitation") ;
      }

      // Increment the source term for the cell attach to a single face.
      void calculate(Entity face) {
        if(massFlux[face]<=0.0) AlphaSourceTerm[ci[face]]-=massFlux[face]*Alpha_f[face] ;
        else AlphaSourceTerm[ci[face]]-=massFlux[face]*(Alpha_f[face]-Alpha[ci[face]]) ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToAlphaSourceTermBoundary>
    registerFOUInviscidFluxToAlphaSourceTermBoundary ;

  // Rule to add the second-order convection contribution to the source term
  // for interior faces.
  class SOUInviscidFluxToAlphaSourceTermInterior : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> Alpha ;
      const_store<vect3d> AlphaGradient ;
      const_store<real> AlphaLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      store<real> AlphaSourceTerm ;
    public:

      // Define input and output.
      SOUInviscidFluxToAlphaSourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("Alpha",Alpha) ;
        name_store("grads(Alpha)",AlphaGradient) ;
        name_store("limiters(Alpha)",AlphaLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("massFlux",massFlux) ;
        name_store("AlphaSourceTerm",AlphaSourceTerm) ;
        input("(cl,cr)->(Alpha,grads(Alpha),limiters(Alpha),cellcenter)") ;
        input("facecenter,massFlux") ;
        output("(cl,cr)->AlphaSourceTerm") ;
        constraint("internalFaces,souCavitation") ;
      }

      // Increment the source term for the cells attach to a single face. Note
      // that here we are insisting that the reconstructed value of Alpha at the
      // face be positive.
      void calculate(Entity face) {
        real secondOrderSource=(massFlux[face]>0.0)? massFlux[face]*
          max(AlphaLimiter[cl[face]]*dot(AlphaGradient[cl[face]],faceCenter[face]-
          cellCenter[cl[face]]),-Alpha[cl[face]]):massFlux[face]*
          max(AlphaLimiter[cr[face]]*dot(AlphaGradient[cr[face]],faceCenter[face]-
          cellCenter[cr[face]]),-Alpha[cr[face]]) ;
        AlphaSourceTerm[cl[face]]-=secondOrderSource ;
        AlphaSourceTerm[cr[face]]+=secondOrderSource ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SOUInviscidFluxToAlphaSourceTermInterior>
    registerSOUInviscidFluxToAlphaSourceTermInterior ;

	
  // Rule to add the second-order convection contribution to the source term
  // for boundary faces.
  class SOUInviscidFluxToAlphaSourceTermBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> Alpha ;
      const_store<vect3d> AlphaGradient ;
      const_store<real> AlphaLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      const_store<real> Alpha_f ;
      store<real> AlphaSourceTerm ;
    public:

      // Define input and output.
      SOUInviscidFluxToAlphaSourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("Alpha",Alpha) ;
        name_store("grads(Alpha)",AlphaGradient) ;
        name_store("limiters(Alpha)",AlphaLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("massFlux",massFlux) ;
        name_store("Alpha_f",Alpha_f) ;
        name_store("AlphaSourceTerm",AlphaSourceTerm) ;
        input("ci->(Alpha,grads(Alpha),limiters(Alpha),cellcenter)") ;
        input("facecenter,massFlux,Alpha_f") ;
        output("ci->AlphaSourceTerm") ;
        constraint("boundaryFaces,souCavitation") ;
      }

      // Increment the source term for the cells attach to a single face. Note
      // that here we are insisting that the reconstructed value of Alpha at the
      // face be positive.
      void calculate(Entity face) {
        AlphaSourceTerm[ci[face]]-=(massFlux[face]>0.0)? massFlux[face]*
          max(AlphaLimiter[ci[face]]*dot(AlphaGradient[ci[face]],faceCenter[face]-
          cellCenter[ci[face]]),-Alpha[ci[face]]):massFlux[face]*Alpha_f[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SOUInviscidFluxToAlphaSourceTermBoundary>
    registerSOUInviscidFluxToAlphaSourceTermBoundary ;

  // Rule to add the evaporation and condensation term to the main coefficient.
  class EvapCondToAlphaMainCoefficient : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_store<real> rho,p,vol,cellRadius ;
      const_param<real> Cdest,Cprod,pv,tInfinite,vInfinite,rhol,rhov ; 
      store<real> AlphaMainCoefficient ;
    public:

      // Define input and output.
      EvapCondToAlphaMainCoefficient() {
        name_store("rhol",rhol) ;
        name_store("rhov",rhov) ;
        name_store("rho",rho) ;
        name_store("p",p) ;
        name_store("Cdest",Cdest) ;
        name_store("Cprod",Cprod) ;
        name_store("pv",pv) ;
        name_store("tInfinite",tInfinite) ;
        name_store("vInfinite",vInfinite) ;
        name_store("vol",vol) ;
        name_store("cellRadius",cellRadius) ;
        name_store("AlphaMainCoefficient",AlphaMainCoefficient) ;
        input("rho,rhol,rhov,p,Cdest,Cprod,pv,tInfinite,vInfinite,vol,cellRadius") ;
        output("AlphaMainCoefficient") ;
        constraint("geom_cells,cavitationModel") ;
      }

      // Increment the main coefficient for a single cell.
      void calculate(Entity cell) {
        if(p[cell]<=pv[cell]){
          AlphaMainCoefficient[cell]-=rho[cell]*(*Cdest)*(p[cell]-(*pv))*vol[cell]*cellRadius[cell]/(0.5*(*rhov)*(*vInfinite)*(*vInfinite)*(*tInfinite)) ;}
        else{
          AlphaMainCoefficient[cell]+=rho[cell]*(*Cprod)*(p[cell]-(*pv))*vol[cell]*cellRadius[cell]/(0.5*(*rhol)*(*vInfinite)*(*vInfinite)*(*tInfinite)) ;}
      }

      // Increment the main coefficient for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<EvapCondToAlphaMainCoefficient> registerEvapCondToAlphaMainCoefficient ;
	
  // Rule to add the evaporation and condensation term to the source term.
  class EvapCondToAlphaSourceTerm : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> rho,p,vol,cellRadius ;
      const_param<real> Cdest,Cprod,pv,tInfinite,vInfinite,rhol,rhov ; 
      store<real> AlphaSourceTerm ;
    public:

      // Define input and output.
      EvapCondToAlphaSourceTerm() {
        name_store("rhol",rhol) ;
        name_store("rhov",rhov) ;
        name_store("rho",rho) ;
        name_store("p",p) ;
        name_store("Cdest",Cdest) ;
        name_store("Cprod",Cprod) ;
        name_store("pv",pv) ;
        name_store("tInfinite",tInfinite) ;
        name_store("vInfinite",vInfinite) ;
        name_store("vol",vol) ;
        name_store("cellRadius",cellRadius) ;
        name_store("AlphaSourceTerm",AlphaSourceTerm) ;
        input("rho,rhol,rhov,p,Cdest,Cprod,pv,tInfinite,vInfinite,vol,cellRadius") ;
        output("AlphaSourceTerm") ;
        constraint("geom_cells,cavitationModel") ;
      }

      // Increment the source term for a single cell.
      void calculate(Entity cell) {
        if(p[cell]<=pv[cell]){
          AlphaSourceTerm[cell]+=0.0 ;}
        else{
          AlphaSourceTerm[cell]+=rho[cell]*(*Cprod)*(p[cell]-(*pv))*vol[cell]*cellRadius[cell]/(0.5*(*rhol)*(*vInfinite)*(*vInfinite)*(*tInfinite)) ;}
      }

      // Increment the source term for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<EvapCondToAlphaSourceTerm> registerEvapCondToAlphaSourceTerm ;
	
  // Rule to compute the diagonal term for the linear system.
  class ComputeAlphaMatrixDiagonal : public pointwise_rule {
    private:
      const_param<real> AlphaRelaxationFactor ;
      const_store<real> AlphaMainCoefficient ;
      store<real> D ;
    public:

      // Define input and output.
      ComputeAlphaMatrixDiagonal() {
        name_store("cavitationRelaxationFactor",AlphaRelaxationFactor) ;
        name_store("AlphaMainCoefficient",AlphaMainCoefficient) ;
        name_store("AlphaStar_D",D) ;
        input("cavitationRelaxationFactor,AlphaMainCoefficient") ;
        output("AlphaStar_D") ;
        constraint("geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) {
        D[cell]=AlphaMainCoefficient[cell]/(*AlphaRelaxationFactor) ;
      }

      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeAlphaMatrixDiagonal> registerComputeAlphaMatrixDiagonal ;

  
  // Rule to initialize the lower terms for the linear system.
  class initializeAlphaMatrixLower : public unit_rule {
    private:
      store<real> L ;
    public:

      // Define input and output.
      initializeAlphaMatrixLower() {
        name_store("AlphaStar_L",L) ;
        output("AlphaStar_L") ;
        constraint("internalFaces") ;
      }

      // initialize for a single face.
      void calculate(Entity face) { L[face]=0.0 ; }

      // initialize for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<initializeAlphaMatrixLower> registerinitializeAlphaMatrixLower ;

  // Rule to add the first-order inviscid flux contribution to the lower terms
  // for the linear system.
  class FOUInviscidFluxToAlphaMatrixLower : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> massFlux ;
      store<real> L ;
    public:

      // Define input and output.
      FOUInviscidFluxToAlphaMatrixLower() {
        name_store("massFlux",massFlux) ;
        name_store("AlphaStar_L",L) ;
        input("massFlux") ;
        output("AlphaStar_L") ;
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

  register_rule<FOUInviscidFluxToAlphaMatrixLower>
    registerFOUInviscidFluxToAlphaMatrixLower ;


  // Rule to initialize the upper terms for the linear system.
  class initializeAlphaMatrixUpper : public unit_rule {
    private:
      store<real> U ;
    public:

      // Define input and output.
      initializeAlphaMatrixUpper() {
        name_store("AlphaStar_U",U) ;
        output("AlphaStar_U") ;
        constraint("internalFaces") ;
      }

      // initialize for a single face.
      void calculate(Entity face) { U[face]=0.0 ; }

      // initialize for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<initializeAlphaMatrixUpper> registerinitializeAlphaMatrixUpper ;

  // Rule to add the first-order inviscid flux contribution to the upper terms
  // for the linear system.
  class FOUInviscidFluxToAlphaMatrixUpper : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> massFlux ;
      store<real> U ;
    public:

      // Define input and output.
      FOUInviscidFluxToAlphaMatrixUpper() {
        name_store("massFlux",massFlux) ;
        name_store("AlphaStar_U",U) ;
        input("massFlux") ;
        output("AlphaStar_U") ;
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

  register_rule<FOUInviscidFluxToAlphaMatrixUpper>
    registerFOUInviscidFluxToAlphaMatrixUpper ;


  // Rule to compute the right-hand side for the linear system.
  class ComputeAlphaRHS : public pointwise_rule {
    private:
      const_param<real> AlphaRelaxationFactor ;
      const_store<real> Alpha ;
      const_store<real> AlphaMainCoefficient ;
      const_store<real> AlphaSourceTerm ;
      store<real> B ;
    public:

      // Define input and output.
      ComputeAlphaRHS() {
        name_store("cavitationRelaxationFactor",AlphaRelaxationFactor) ;
        name_store("Alpha",Alpha) ;
        name_store("AlphaMainCoefficient",AlphaMainCoefficient) ;
        name_store("AlphaSourceTerm",AlphaSourceTerm) ;
        name_store("AlphaStar_B",B) ;
        input("cavitationRelaxationFactor,Alpha,AlphaMainCoefficient,AlphaSourceTerm") ;
        output("AlphaStar_B") ;
        constraint("geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) {
        B[cell]=AlphaSourceTerm[cell]+(1.0-(*AlphaRelaxationFactor))*
          AlphaMainCoefficient[cell]*Alpha[cell]/(*AlphaRelaxationFactor) ;
      }

      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeAlphaRHS> registerComputeAlphaRHS ;

//-----------------------------------------------------------------------------
// Rules for computing the residual of the Alpha equation.

  // Rule to initialize the Alpha residual.
  class initializeAlphaResidual : public unit_rule {
    private:
      store<real> AlphaResidual ;
    public:

      // Define input and output.
      initializeAlphaResidual() {
        name_store("AlphaResidual",AlphaResidual) ;
        output("AlphaResidual") ;
//        constraint("vol") ;
        constraint("geom_cells") ;
      }

      // initialize the residual for a single cell.
      void calculate(Entity cell) { AlphaResidual[cell]=0.0 ; }

      // initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<initializeAlphaResidual> registerinitializeAlphaResidual ;

  // Rule to compute the Alpha residual for each cell.
  class ComputeAlphaResidualOne : public apply_rule<store<real>,Loci::Summation
  <real> > {
    private:
      const_store<real> D ;
      const_store<real> Alpha ;
      const_store<real> B ;
      store<real> AlphaResidual ;
    public:

      // Define input and output.
      ComputeAlphaResidualOne() {
        name_store("AlphaStar_D",D) ;
        name_store("Alpha",Alpha) ;
        name_store("AlphaStar_B",B) ;
        name_store("AlphaResidual",AlphaResidual) ;
        input("AlphaStar_D,Alpha,AlphaStar_B") ;
        output("AlphaResidual") ;
//        constraint("vol") ;
        constraint("geom_cells") ;
      }

      // initialize the residual for a single cell.
      void calculate(Entity cell) {
        AlphaResidual[cell]+=B[cell]-D[cell]*Alpha[cell] ;
      }

      // initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ;}
  } ;

  register_rule<ComputeAlphaResidualOne> registerComputeAlphaResidualOne ;

  // Rule to compute the Alpha residual for each cell.
  class ComputeAlphaResidualTwo : public apply_rule<store<real>,Loci::Summation
  <real> > {
    private:
      const_Map cl,cr ;
      const_store<real> Alpha ;
      const_store<real> L,U ;
      store<real> AlphaResidual ;
    public:

      // Define input and output.
      ComputeAlphaResidualTwo() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("Alpha",Alpha) ;
        name_store("AlphaStar_L",L) ;
        name_store("AlphaStar_U",U) ;
        name_store("AlphaResidual",AlphaResidual) ;
        input("(cl,cr)->Alpha,AlphaStar_L,AlphaStar_U") ;
        output("(cl,cr)->AlphaResidual") ;
        constraint("internalFaces") ;
      }

      // Add the neighbor contribution to the residual for each of the two
      // cells on either side of the face.
      void calculate(Entity face) {
        AlphaResidual[cl[face]]-= U[face]*Alpha[cr[face]] ;
        AlphaResidual[cr[face]]-= L[face]*Alpha[cl[face]] ;
      }

      // Add the neighbor contributions for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ;}
  } ;

  register_rule<ComputeAlphaResidualTwo> registerComputeAlphaResidualTwo ;


  // Rule to initialize the total Alpha residual.
  class InitializeTotalAlphaResidual : public unit_rule {
    private:
      param<ScalarResidual> AlphaResidualData ;
    public:

      // Define input and output.
      InitializeTotalAlphaResidual() {
        name_store("AlphaResidualData",AlphaResidualData) ;
        output("AlphaResidualData") ;
        constraint("cavitationModel,geom_cells") ;
      }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) {
        AlphaResidualData=ScalarResidual() ;
      }
  } ;

  register_rule<InitializeTotalAlphaResidual> registerInitializeTotalAlphaResidual ;

    // Rule to compute the total Alpha residual.
  class ComputeTotalAlphaResidual : public apply_rule<param<ScalarResidual>,
  ScalarResidualJoin> {
    private:
      const_store<real> AlphaResidual ;
      const_store<real> D ;
      const_store<vect3d> cellCenter ;
      param<ScalarResidual> AlphaResidualData ;
    public:

      // Define input and output.
      ComputeTotalAlphaResidual() {
        name_store("AlphaResidual",AlphaResidual) ;
        name_store("cellcenter",cellCenter) ;
        name_store("AlphaResidualData",AlphaResidualData) ;
        input("AlphaResidual,cellcenter") ;
        output("AlphaResidualData") ;
        constraint("cavitationModel,geom_cells") ;
      }

      // Add the cell contribution to the residual for a single cell.
      void calculate(Entity cell) {
        ScalarResidual temp ;
        temp.maxResidual=AlphaResidual[cell] ;
        temp.totalResidual=abs(AlphaResidual[cell]) ;
        temp.maxResidualLocation=cellCenter[cell] ;
        join(*AlphaResidualData,temp) ;
      }

      // Add the cell contribution to the residual for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeTotalAlphaResidual> registerComputeTotalAlphaResidual ;

//-----------------------------------------------------------------------------
// Rules to limit Alpha.

  // Limit Alpha to positive values.
  class LimitAlpha: public pointwise_rule {
    private:
      store<real> AlphaStar ;
    public:

      // Define input and output.
      LimitAlpha() {
        name_store("AlphaStar",AlphaStar) ;
        input("AlphaStar") ;
        output("AlphaStarLimited=AlphaStar") ;
        constraint("geom_cells") ;
      }

      // Limit Alpha for a single cell.
      void calculate(Entity cell) {
        if(AlphaStar[cell]<0.0) AlphaStar[cell]=0. ;
        if(AlphaStar[cell]>1.0) AlphaStar[cell]=1. ;
      }

      // Limit Alpha for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<LimitAlpha> registerLimitAlpha ;

//-----------------------------------------------------------------------------
// Rules for marching the Alpha equations. 

  // Time build rule for Alpha when using BDF2 time integrator. Although
  // this rule sets Alpha{n=-1} from Alpha_ic, these values are not really
  // used since BDF is used on the first timestep for non-restarts. 
  // The only purpose for this rule is to let Loci know that there are
  // two previous time-levels that need to be maintained for BDF2.
  class TimeBuildAlphaBDF2: public pointwise_rule {
    private:
      const_store<real> Alpha_ic ;
      store<real> Alpha;
    public:

      // Define input and output.
      TimeBuildAlphaBDF2() {
        name_store("Alpha_ic",Alpha_ic) ;
        name_store("Alpha{n=-1}",Alpha) ;
        input("Alpha_ic") ;
        output("Alpha{n=-1}") ;
        constraint("geom_cells,cavitationModel") ;
      }

      // Assign Alpha for a single cell.
      void calculate(Entity cell) {
        Alpha[cell]=Alpha_ic[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildAlphaBDF2> registerTimeBuildAlphaBDF2 ;

  // Time build rule for Alpha.
  class TimeBuildAlpha: public pointwise_rule {
    private:
      const_store<real> Alpha_ic ;
      store<real> AlphaTimeStepZero ;
    public:

      // Define input and output.
      TimeBuildAlpha() {
        name_store("Alpha_ic",Alpha_ic) ;
        name_store("Alpha{n=0}",AlphaTimeStepZero) ;
        input("Alpha_ic") ;
        output("Alpha{n=0}") ;
        constraint("geom_cells,cavitationModel") ;
      }

      // Assign Alpha at time Zero for a single cell.
      void calculate(Entity cell) {
        AlphaTimeStepZero[cell]=Alpha_ic[cell] ;
      }

      // Assign Alpha at time Zero for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildAlpha> registerTimeBuildAlpha ;

  // Iteration build rule for Alpha.
  class IterationBuildAlpha : public pointwise_rule {
    private:
      const_store<real> AlphaTimeStepN ;
      store<real> AlphaIterationZero ;
    public:

      // Define input and output.
      IterationBuildAlpha() {
        name_store("Alpha{n}",AlphaTimeStepN) ;
        name_store("Alpha{n,it=0}",AlphaIterationZero) ;
        input("Alpha{n}") ;
        output("Alpha{n,it=0}") ;
        constraint("geom_cells{n},cavitationModel{n}") ;
      }

      // Assign Alpha at iteration Zero for a single cell.
      void calculate(Entity cell) {
        AlphaIterationZero[cell]=AlphaTimeStepN[cell] ;
      }

      // Assign Alpha at iteration Zero for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationBuildAlpha> registerIterationBuildAlpha ;

  // Rule to add temporal component of the Alpha equation to the main coefficient.
  class TemporalToAlphaMainCoefficient: public apply_rule<store<real>,Loci::
  Summation<real> > {
    private:
      const_param<real> timeIntegratorFactor0 ;
      const_store<real> timeStepFactor ;
      const_store<real> vol ;
      const_store<real> rho ;
      const_store<real> cellRadius ;
      store<real> AlphaMainCoefficient ;
    public:

      // Define input and output.
      TemporalToAlphaMainCoefficient() {
        name_store("timeIntegratorFactor0{n}",timeIntegratorFactor0) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n}",rho) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("AlphaMainCoefficient{n,it}",AlphaMainCoefficient) ;
        input("rho{n},vol{n},cellRadius{n,it}") ;
        input("timeStepFactor{n},timeIntegratorFactor0{n}") ;
        output("AlphaMainCoefficient{n,it}") ;
        constraint("geom_cells") ;
      }

      // Add temporal component for a single cell.
      void calculate(Entity cell) {
        AlphaMainCoefficient[cell]+=rho[cell]*vol[cell]*cellRadius[cell]*
          (*timeIntegratorFactor0)/timeStepFactor[cell] ;
      }

      // Add temporal component for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToAlphaMainCoefficient> registerTemporalToAlphaMainCoefficient ;

  // Rule to add temporal component of the Alpha equation to the source term.
  class TemporalToAlphaSourceTerm : public apply_rule<store<real>,Loci::
  Summation<real> > {
    private:
      const_param<real> timeIntegratorFactor0 ;
      const_store<real> timeStepFactor ;
      const_store<real> rho ;
      const_store<real> Alpha ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> AlphaSourceTerm ;
    public:

      // Define input and output.
      TemporalToAlphaSourceTerm() {
        name_store("timeIntegratorFactor0{n}",timeIntegratorFactor0) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n}",rho) ;
        name_store("Alpha{n}",Alpha) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("AlphaSourceTerm{n,it}",AlphaSourceTerm) ;
        input("rho{n},Alpha{n},vol{n},cellRadius{n,it}") ;
        input("timeStepFactor{n},timeIntegratorFactor0{n}") ;
        output("AlphaSourceTerm{n,it}") ;
        constraint("geom_cells") ;
      }

      // Add temporal component to source term for a single cell.
      void calculate(Entity cell) {
        AlphaSourceTerm[cell]+=(rho[cell]*vol[cell]*cellRadius[cell]*
          (*timeIntegratorFactor0)/timeStepFactor[cell])*Alpha[cell] ;
      }

      // Add temporal component to source term for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToAlphaSourceTerm> registerTemporalToAlphaSourceTerm ;

  // Rule to add temporal component of momentum equation to the source term
  // for the BDF2 scheme.
  class TemporalToAlphaSourceTermBDF2 : public apply_rule<store<real>,Loci::
  Summation<real> > {
    private:
      const_param<real> timeIntegratorFactor1 ;
      const_store<real> timeStepFactor ;
      const_store<real> rhoOld ;
      const_store<real> AlphaOld,Alpha ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> AlphaSourceTerm ;
    public:

      // Define input and output.
      TemporalToAlphaSourceTermBDF2() {
        name_store("timeIntegratorFactor1{n}",timeIntegratorFactor1) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n-1}",rhoOld) ;
        name_store("Alpha{n-1}",AlphaOld) ;
        name_store("Alpha{n,it}",Alpha) ;
        name_store("vol{n,it}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("AlphaSourceTerm{n,it}",AlphaSourceTerm) ;
        input("rho{n-1},Alpha{n-1},Alpha{n,it},vol{n,it}") ;
        input("cellRadius{n,it},timeStepFactor{n},timeIntegratorFactor1{n}") ;
        output("AlphaSourceTerm{n,it}") ;
        constraint("geom_cells,BDF2Integrator") ;
      }

      // Add temporal component to source term for a single cell.
      void calculate(Entity cell) {
        AlphaSourceTerm[cell]+=((*timeIntegratorFactor1)*rhoOld[cell]*
          vol[cell]*cellRadius[cell]/timeStepFactor[cell])*
          (Alpha[cell]-AlphaOld[cell]) ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToAlphaSourceTermBDF2> registerTemporalToAlphaSourceTermBDF2 ;

  // Iteration advance rule for Alpha. 
  class IterationAdvanceAlpha : public pointwise_rule {
    private:
      const_store<real> AlphaStar ;
      store<real> AlphaIterationPlusOne ;
    public:

      // Define input and output.
      IterationAdvanceAlpha() {
        name_store("AlphaStarLimited{n,it}",AlphaStar) ;
        name_store("Alpha{n,it+1}",AlphaIterationPlusOne) ;
        input("AlphaStarLimited{n,it}") ;
        output("Alpha{n,it+1}") ;
        constraint("geom_cells{n,it},cavitationModel{n,it}") ;
      }

      // Assign Alpha at end of iteration for a single cell.
      void calculate(Entity cell) {
        AlphaIterationPlusOne[cell]=AlphaStar[cell] ;
      }

      // Assign Alpha at end of iteration for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationAdvanceAlpha> registerIterationAdvanceAlpha ;

  // Iteration collapse rule for Alpha.
  class IterationCollapseAlpha : public pointwise_rule {
    private:
      const_param<bool> iterationFinished ;
      store<real> Alpha ;
    public:

      // Define input and output.
      IterationCollapseAlpha() {
        name_store("iterationFinished{n,it-1}",iterationFinished) ;
        name_store("Alpha{n,it}",Alpha) ;
        input("iterationFinished{n,it-1}") ;
        input("Alpha{n,it}") ;
        output("Alpha{n+1}=Alpha{n,it}") ;
        conditional("iterationFinished{n,it-1}") ;
        constraint("geom_cells{n,it},cavitationModel{n,it}") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<IterationCollapseAlpha> registerIterationCollapseAlpha ;

//-----------------------------------------------------------------------------

  // Iteration advance rule for density
  class IterationAdvanceDensityCavitation : public pointwise_rule {
    private:
      const_param<real> rhol,rhov ;
      store<real> rho ;
      store<real> Alpha ;
    public:

     // Define input and output.
     IterationAdvanceDensityCavitation() {
        name_store("rhol",rhol) ;
        name_store("rhov",rhov) ;
        name_store("Alpha{n,it}",Alpha) ;
        name_store("cavitation::rho{n,it+1}",rho) ;
        input("rhol,rhov,Alpha{n,it}") ;
        output("cavitation::rho{n,it+1}") ;
        constraint("geom_cells{n,it},cavitationModel{n,it}") ;
      }

      // Assign density at end of iteration for a single cell.
      void calculate(Entity cell) { rho[cell]=(*rhol)*Alpha[cell]+(*rhov)*(1.0-Alpha[cell]) ; }

      // Assign density at end of iteration for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationAdvanceDensityCavitation>
    registerIterationAdvanceDensityCavitation ;
	
}
