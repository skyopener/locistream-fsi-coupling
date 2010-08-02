//-----------------------------------------------------------------------------
// Description: This file contains rules for flamelet chemistry model
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
#include "flameletVarsFileInputs.h"

// Include the flamelet table class and create one
#include "flamelet_table.h"
Flamelet_Table flamelet_table;

namespace streamUns {

//-----------------------------------------------------------------------------
// Rules to process flamelet equation options from the .vars file.

  // Creates the constraints for the flamelet equations.
  class FlameletModelConstraints : public constraint_rule {
    private:
      const_param<FlameletEquationOptions> flameletEquationOptions ;
      Constraint flameletModel ;
    public:
                                                                                
      // Define input and output.
      FlameletModelConstraints() {
        name_store("flameletEquationOptions",flameletEquationOptions) ;
        name_store("flameletModel",flameletModel) ;
        input("flameletEquationOptions") ;
        output("flameletModel") ;
      }
                                                                                
      // Set up the constraints.
      virtual void compute(const sequence& seq) {
         if((*flameletEquationOptions).optionExists("table")){
		 flameletModel=~EMPTY;
         }
       }
  } ;
                                                                                
  register_rule<FlameletModelConstraints>
    registerFlameletModelConstraints ;


  // Creates the flamelet equation solver constraints.
  class FlameletEquationSolverConstraints : public constraint_rule {
    private:
      const_param<FlameletEquationOptions> flameletEquationOptions ;
      Constraint flameletSGSLinearSolver ;
    public:
                                                                                
      // Define input and output.
      FlameletEquationSolverConstraints() {
        name_store("flameletEquationOptions",flameletEquationOptions) ;
        name_store("flameletSGSLinearSolver",flameletSGSLinearSolver) ;
        input("flameletEquationOptions") ;
        output("flameletSGSLinearSolver") ;
      }
                                                                                
      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        if((*flameletEquationOptions).optionExists("linearSolver")){
          Loci::option_value_type optionValueType=flameletEquationOptions->
            getOptionValueType("linearSolver") ;
          switch(optionValueType){
            case Loci::NAME:
              {
                Loci::option_values optionValues=flameletEquationOptions->
                  getOption("linearSolver") ;
                string name ; optionValues.get_value(name) ;
                if(name=="SGS"){
                  flameletSGSLinearSolver=~EMPTY ;
                }else{
                  cerr << "Bad linearSolver for flameletEquation." << endl ;
                  Loci::Abort() ;
                }
              }
              break ;
            default:
              cerr << "Bad type for linearSolver in flameletEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          flameletSGSLinearSolver=~EMPTY ;
        }
      }
                                                                                
  } ;
                                                                                
  register_rule<FlameletEquationSolverConstraints>
    registerFlameletEquationSolverConstraints ;

  // Creates the flamelet inviscid flux constraints.
  class FlameletInviscidFluxConstraints : public constraint_rule {
    private:
      const_param<string> flameletInviscidFlux ;
      Constraint fouFlamelet,souFlamelet ;
    public:

      // Define input and output.
      FlameletInviscidFluxConstraints() {
        name_store("flameletInviscidFlux",flameletInviscidFlux) ;
        name_store("fouFlamelet",fouFlamelet) ;
        name_store("souFlamelet",souFlamelet) ;
        input("flameletInviscidFlux") ;
        output("fouFlamelet,souFlamelet") ;
      }

      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        if(*flameletInviscidFlux!="FOU" && *flameletInviscidFlux!="SOU"){
          cerr << "Bad flameletInviscidFlux in .vars file." << endl ;
          Loci::Abort() ;
        }
        fouFlamelet=(*flameletInviscidFlux=="FOU")? ~EMPTY:EMPTY ;
        souFlamelet=(*flameletInviscidFlux=="SOU")? ~EMPTY:EMPTY ;
      }
  } ;       

  register_rule<FlameletInviscidFluxConstraints>
    registerFlameletInviscidFluxConstraints ;
	
  // Creates the flamelet equation solver parameters.
  class FlameletEquationSolverParameters : public singleton_rule {
    private:
      const_param<FlameletEquationOptions> flameletEquationOptions ;
      param<int> flameletMaxIterations ;
      param<real> flameletRelaxationFactor ;
    public:
                                                                                
      // Define input and output.
      FlameletEquationSolverParameters() {
        name_store("flameletEquationOptions",flameletEquationOptions) ;
        name_store("flameletMaxIterations",flameletMaxIterations) ;
        name_store("flameletRelaxationFactor",flameletRelaxationFactor) ;
        input("flameletEquationOptions") ;
        output("flameletMaxIterations,flameletRelaxationFactor") ;
      }
                                                                                
      // Set up the parameters.
      virtual void compute(const sequence& seq) {
                
	// Get table name and initialize the table
	Loci::option_values optionValues=flameletEquationOptions->getOption("table") ;
        string tableFileName;
        optionValues.get_value(tableFileName) ;
        flamelet_table.read(tableFileName);

                                                                
        // Relaxation factor.
        if((*flameletEquationOptions).optionExists("relaxationFactor")){
          Loci::option_value_type optionValueType=flameletEquationOptions->
            getOptionValueType("relaxationFactor") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                flameletEquationOptions->getOption("relaxationFactor",temp) ;
                if(temp<=0.0 || temp>1.0){
                  cerr << "Bad relaxationFactor for flameletEquation."
                    << endl ; Loci::Abort() ;
                }
                *flameletRelaxationFactor=temp ;
              }
              break ;
            default:
              cerr << "Bad type for relaxationFactor in flameletEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *flameletRelaxationFactor=0.5 ;
        }
                                                                                
        // Maximum number of iterations.
        if((*flameletEquationOptions).optionExists("maxIterations")){
          Loci::option_value_type optionValueType=flameletEquationOptions->
            getOptionValueType("maxIterations") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                flameletEquationOptions->getOption("maxIterations",temp) ;
                if(int(temp)<0){
                  cerr << "Bad maxIterations value for flameletEquation."
                    << endl ; Loci::Abort() ;
                }
                *flameletMaxIterations=int(temp) ;
              }
              break ;
            default:
              cerr << "Bad type for maxIterations in flameletEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *flameletMaxIterations=5 ;
        }
		
      }
  } ;
                                                                                
  register_rule<FlameletEquationSolverParameters>
    registerFlameletEquationSolverParamters ;

//-----------------------------------------------------------------------------
// Rules for the constants for the flamelet model.

  class FlameletModelConstants : public singleton_rule {
    private:
      param<real> Cg ; // Constant for production of Zvar
      param<real> Cd ; // Constant for destruction of Zvar
    public:

      // Define input and output.
      FlameletModelConstants() {
        name_store("Cg",Cg) ;
        name_store("Cd",Cd) ;
        output("Cg,Cd") ;
        constraint("flameletModel") ;
      }

      // Compute the model constants.
      virtual void compute(const sequence &seq) {
	 *Cg=2.86;
	 *Cd=2.0;
      }
  } ;

  register_rule<FlameletModelConstants> registerFlameletModelConstants ;
  

//-----------------------------------------------------------------------------
// Viscosity of Z and Zvar equations

  // Species viscosity for the interior faces.
  class ZViscosityInterior : public pointwise_rule {
	  private:
		  const_Map cl,cr ;
		  const_param<real> turbulentSchmidtNumber ;
		  const_store<real> diffusivity,eddyViscosity ;
		  store<real> ZViscosity ;
	  public:

      // Define input and output.
	ZViscosityInterior() {
		name_store("cl",cl) ;
		name_store("cr",cr) ;
		name_store("turbulentSchmidtNumber",turbulentSchmidtNumber) ;
		name_store("eddyViscosity",eddyViscosity) ;
		name_store("diffusivity",diffusivity) ;
		name_store("ZViscosity",ZViscosity) ;
		input("turbulentSchmidtNumber") ;
		input("(cl,cr)->(diffusivity,eddyViscosity)") ;
		output("ZViscosity") ;
		constraint("internalFaces") ;
	}

      // Compute for each face.
	void calculate(Entity face) {
		ZViscosity[face]=0.5*(diffusivity[cl[face]]+
				diffusivity[cr[face]]+(eddyViscosity
				[cl[face]]+eddyViscosity[cr[face]])/(*turbulentSchmidtNumber)) ;
	}

      // Loop over faces.
	virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ZViscosityInterior> registerZViscosityInterior ;

   // Species viscosity for the boundary faces.
  class ZViscosityBoundary : public pointwise_rule {
	  private:
		  const_param<real> turbulentSchmidtNumber ;
		  const_store<real> diffusivity_f,eddyViscosity_f ;
		  store<real> ZViscosity ;
	  public:

         // Define input and output.
	ZViscosityBoundary() {
		name_store("turbulentSchmidtNumber",turbulentSchmidtNumber) ;
		name_store("eddyViscosity_f",eddyViscosity_f) ;
		name_store("diffusivity_f",diffusivity_f) ;
		name_store("ZViscosity",ZViscosity) ;
		input("turbulentSchmidtNumber") ;
		input("diffusivity_f,eddyViscosity_f") ;
		output("ZViscosity") ;
		constraint("boundaryFaces") ;
	}

      // Compute for each face.
	void calculate(Entity face) {
		ZViscosity[face]=diffusivity_f[face]+
				eddyViscosity_f[face]/(*turbulentSchmidtNumber) ;
	}

      // Loop over faces.
	virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ZViscosityBoundary> registerZViscosityBoundary ;

//-----------------------------------------------------------------------------
// Boundary condition rules for Z.

  // Rule for boundary faces with specified Z. Assigns value to all boundary
  // faces that have the property Z_BC.
  class BoundaryZSpecification : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<real> Z_BC ;
      store<real> Z_f ;
    public:

      // Define input and output.
      BoundaryZSpecification() {
        name_store("ref",ref) ;
        name_store("Z_BC",Z_BC) ;
        name_store("Z_f",Z_f) ;
        input("ref->Z_BC") ;
        output("Z_f") ;
      }

      // Calculate Z for a single face.
      void calculate(Entity face) { Z_f[face]=Z_BC[ref[face]] ; }

      // Calculate Z for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryZSpecification>
    registerBoundaryZSpecification ;

  // Rule for extrapolating Z to boundary faces. This occurs for all outlets,
  // noslip, slip and symmetry boundaries. Right now we are using the low-order method.
  class BoundaryZExtrapolation : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> Z ;
      store<real> Z_f ;
    public:

      // Define input and output.
      BoundaryZExtrapolation() {
        name_store("ci",ci) ;
        name_store("Z",Z) ;
        name_store("Z_f",Z_f) ;
        input("ci->Z") ;
        output("Z_f") ;
        constraint("extrapolatedZ_BC") ;
      }

      // Calculate Z for a single face.
      void calculate(Entity face) {
        Z_f[face]=Z[ci[face]] ;
      }

      // Calculate Z for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryZExtrapolation> registerBoundaryZExtrapolation ;


//-----------------------------------------------------------------------------
// Boundary condition rules for Zvar.

  // Rule for boundary faces with specified Zvar. Assigns Zvar value to all
  // boundary faces that have the property Zvar_BC.
  class BoundaryZvarSpecification : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<real> Zvar_BC ;
      store<real> Zvar_f ;
    public:

      // Define input and output.
      BoundaryZvarSpecification() {
        name_store("ref",ref) ;
        name_store("Zvar_BC",Zvar_BC) ;
        name_store("Zvar_f",Zvar_f) ;
        input("ref->Zvar_BC") ;
        output("Zvar_f") ;
      }

      // Calculate Zvar for a single face.
      void calculate(Entity face) { Zvar_f[face]=Zvar_BC[ref[face]] ; }

      // Calculate Zvar for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryZvarSpecification>
    registerBoundaryZvarSpecification ;

  // Rule for extrapolating Zvar to boundary faces. This occurs for all
  // outlets, nolip, slip and symmetry boundaries. Right now we are using the
  // low-order method.
  class BoundaryZvarExtrapolation : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> Zvar ;
      store<real> Zvar_f ;
    public:

      // Define input and output.
    BoundaryZvarExtrapolation() {
        name_store("ci",ci) ;
        name_store("Zvar",Zvar) ;
		name_store("Zvar_f",Zvar_f) ;
		input("ci->Zvar") ;
		output("Zvar_f") ;
		constraint("extrapolatedZvar_BC") ;
      }

      // Calculate Zvar for a single face.
      void calculate(Entity face) {
	      Zvar_f[face]=Zvar[ci[face]] ;
      }

      // Calculate Zvar for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryZvarExtrapolation> registerBoundaryZvarExtrapolation ;

//-----------------------------------------------------------------------------
// Rules to create a constraint for boundary faces with non-zero diffusion
// flux for flamelet quantities.

  // All inlet faces have non-zero diffusion flux.
  class BoundaryFlameletDiffusionInlet : public pointwise_rule {
    private:
      store<bool> boundaryFlameletDiffusion ;
    public:

      // Define input and output.
      BoundaryFlameletDiffusionInlet() {
        name_store("boundaryFlameletDiffusion",boundaryFlameletDiffusion) ;
        output("boundaryFlameletDiffusion") ;
        constraint("inlet_BC") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { boundaryFlameletDiffusion[face]=true ; }

      // Assign flag for a sequence of boundary faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryFlameletDiffusionInlet>
    registerBoundaryFlameletDiffusionInlet ;

  // All outlet faces have non-zero diffusion flux.
  class BoundaryFlameletDiffusionOutlet : public pointwise_rule {
    private:
      store<bool> boundaryFlameletDiffusion ;
    public:

      // Define input and output.
      BoundaryFlameletDiffusionOutlet() {
        name_store("boundaryFlameletDiffusion",boundaryFlameletDiffusion) ;
        output("boundaryFlameletDiffusion") ;
        constraint("outlet_BC") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { boundaryFlameletDiffusion[face]=true ; }

      // Assign flag for a sequence of boundary faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryFlameletDiffusionOutlet>
    registerBoundaryFlameletDiffusionOutlet ;

//-----------------------------------------------------------------------------
// Rules for assembling the Z equation.

  // Rule to initialize the main coefficient.
  class InitializeZMainCoefficient : public unit_rule {
    private:
      store<real> ZMainCoefficient ;
    public:

      // Define input and output.
      InitializeZMainCoefficient() {
        name_store("ZMainCoefficient",ZMainCoefficient) ;
        output("ZMainCoefficient") ;
        constraint("geom_cells") ;
      }

      // Set the main coefficient to zero for a single cell.
      void calculate(Entity cell) { ZMainCoefficient[cell]=0.0 ; }

      // Set the main coefficient to zero for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeZMainCoefficient> registerInitializeZMainCoefficient ;

  // Rule to add the first-order inviscid flux contribution to the main
  // coefficient for interior faces.
  class FOUInviscidFluxToZMainCoefficientInterior : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> massFlux ;
      store<real> ZMainCoefficient ;
    public:

      // Define input and output.
      FOUInviscidFluxToZMainCoefficientInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("massFlux",massFlux) ;
        name_store("ZMainCoefficient",ZMainCoefficient) ;
        input("massFlux") ;
        output("cl->ZMainCoefficient,cr->ZMainCoefficient") ;
        constraint("internalFaces") ;
      }

      // Increment the main coefficient for cells attached to a single face.
      void calculate(Entity face) {
        if(massFlux[face]>0.0){
          ZMainCoefficient[cr[face]]+=massFlux[face] ;
        }else{
          ZMainCoefficient[cl[face]]-=massFlux[face] ;
        }
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToZMainCoefficientInterior>
    registerFOUInviscidFluxToZMainCoefficientInterior ;

  // Rule to add the first-order inviscid flux contribution to the main
  // coefficient for boundary faces.
  class FOUInviscidFluxToZMainCoefficientBoundary : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> massFlux ;
      store<real> ZMainCoefficient ;
    public:

      // Define input and output.
      FOUInviscidFluxToZMainCoefficientBoundary() {
        name_store("ci",ci) ;
        name_store("massFlux",massFlux) ;
        name_store("ZMainCoefficient",ZMainCoefficient) ;
        input("massFlux") ;
        output("ci->ZMainCoefficient") ;
        constraint("boundaryFaces") ;
      }

      // Increment the main coefficient for the cell attach to a single face.
      void calculate(Entity face) {
        if(massFlux[face]<0.0) ZMainCoefficient[ci[face]]-=massFlux[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToZMainCoefficientBoundary>
    registerFOUInviscidFluxToZMainCoefficientBoundary ;

  // Rule to add the diffusive flux contribution to the main coefficient for
  // interior faces.
  class DiffusiveFluxToZMainCoefficientInterior : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> faceRadius ;
      const_store<real> diffusionProduct ;
      const_store<real> ZViscosity ;
      store<real> ZMainCoefficient ;
    public:

      // Define input and output.
      DiffusiveFluxToZMainCoefficientInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("faceRadius",faceRadius) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("ZViscosity",ZViscosity) ;
        name_store("ZMainCoefficient",ZMainCoefficient) ;
        input("diffusionProduct,faceRadius,ZViscosity") ;
        output("(cl,cr)->ZMainCoefficient") ;
        constraint("internalFaces") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real temp=ZViscosity[face]*diffusionProduct[face]*faceRadius[face] ;
        ZMainCoefficient[cl[face]]+=temp ; ZMainCoefficient[cr[face]]+=temp ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToZMainCoefficientInterior>
    registerDiffusiveFluxToZMainCoefficientInterior ;

  // Rule to add the diffusive flux contribution to the main coefficient for
  // boundary faces.
  class DiffusiveFluxToZMainCoefficientBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
	    const_Map ci ;
	    const_store<real> faceRadius ;
	    const_store<real> diffusionProduct ;
	    const_store<real> ZViscosity ;
	    store<real> ZMainCoefficient ;
    public:

      // Define input and output.
      DiffusiveFluxToZMainCoefficientBoundary() {	
	name_store("ci",ci) ;
	name_store("faceRadius",faceRadius) ;
	name_store("diffusionProduct",diffusionProduct) ;
	name_store("ZViscosity",ZViscosity) ;
	name_store("ZMainCoefficient",ZMainCoefficient) ;
	input("diffusionProduct,faceRadius,ZViscosity") ;
	output("ci->ZMainCoefficient") ;
        constraint("boundaryFlameletDiffusion") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
	ZMainCoefficient[ci[face]]+=ZViscosity[face]*diffusionProduct[face]*faceRadius[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToZMainCoefficientBoundary>
    registerDiffusiveFluxToZMainCoefficientBoundary ;

  // Rule to initialize the source term.
  class InitializeZSourceTerm : public unit_rule {
	  private:
		store<real> ZSourceTerm ;
	  public:

      // Define input and output.
	InitializeZSourceTerm() {
		name_store("ZSourceTerm",ZSourceTerm) ;
		output("ZSourceTerm") ;
		constraint("geom_cells") ;
	}

      // Set the source term to zero for a single cell.
	void calculate(Entity cell) { ZSourceTerm[cell]=0.0 ; }

      // Set the source term to zero for a sequence of cells.
	virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeZSourceTerm> registerInitializeZSourceTerm ;
  
  
  // Rule to add the first-order inviscid flux contribution to the source term
  // for boundary faces.
  class FOUInviscidFluxToZSourceTermBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> Z ;
      const_store<real> massFlux ;
      const_store<real> Z_f ;
      store<real> ZSourceTerm ;
    public:

      // Define input and output.
      FOUInviscidFluxToZSourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("Z",Z) ;
        name_store("massFlux",massFlux) ;
        name_store("Z_f",Z_f) ;
        name_store("ZSourceTerm",ZSourceTerm) ;
        input("ci->Z,massFlux,Z_f") ;
        output("ci->ZSourceTerm") ;
        constraint("boundaryFaces,fouFlamelet") ;
      }

      // Increment the source term for the cell attach to a single face.
      void calculate(Entity face) {
        if(massFlux[face]<=0.0) ZSourceTerm[ci[face]]-=massFlux[face]*Z_f[face] ;
        else ZSourceTerm[ci[face]]-=massFlux[face]*(Z_f[face]-Z[ci[face]]) ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToZSourceTermBoundary>
    registerFOUInviscidFluxToZSourceTermBoundary ;

  // Rule to add the second-order convection contribution to the source term
  // for interior faces.
  class SOUInviscidFluxToZSourceTermInterior : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> Z ;
      const_store<vect3d> ZGradient ;
      const_store<real> ZLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      store<real> ZSourceTerm ;
    public:

      // Define input and output.
      SOUInviscidFluxToZSourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("Z",Z) ;
        name_store("grads(Z)",ZGradient) ;
        name_store("limiters(Z)",ZLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("massFlux",massFlux) ;
        name_store("ZSourceTerm",ZSourceTerm) ;
        input("(cl,cr)->(Z,grads(Z),limiters(Z),cellcenter)") ;
        input("facecenter,massFlux") ;
        output("(cl,cr)->ZSourceTerm") ;
        constraint("internalFaces,souFlamelet") ;
      }

      // Increment the source term for the cells attach to a single face. Note
      // that here we are insisting that the reconstructed value of Z at the
      // face be positive.
      void calculate(Entity face) {
        real secondOrderSource=(massFlux[face]>0.0)? massFlux[face]*
          max(ZLimiter[cl[face]]*dot(ZGradient[cl[face]],faceCenter[face]-
          cellCenter[cl[face]]),-Z[cl[face]]):massFlux[face]*
          max(ZLimiter[cr[face]]*dot(ZGradient[cr[face]],faceCenter[face]-
          cellCenter[cr[face]]),-Z[cr[face]]) ;
        ZSourceTerm[cl[face]]-=secondOrderSource ;
        ZSourceTerm[cr[face]]+=secondOrderSource ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SOUInviscidFluxToZSourceTermInterior>
    registerSOUInviscidFluxToZSourceTermInterior ;

  // Rule to add the second-order convection contribution to the source term
  // for boundary faces.
  class SOUInviscidFluxToZSourceTermBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> Z ;
      const_store<vect3d> ZGradient ;
      const_store<real> ZLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      const_store<real> Z_f ;
      store<real> ZSourceTerm ;
    public:

      // Define input and output.
      SOUInviscidFluxToZSourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("Z",Z) ;
        name_store("grads(Z)",ZGradient) ;
        name_store("limiters(Z)",ZLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("massFlux",massFlux) ;
        name_store("Z_f",Z_f) ;
        name_store("ZSourceTerm",ZSourceTerm) ;
        input("ci->(Z,grads(Z),limiters(Z),cellcenter)") ;
        input("facecenter,massFlux,Z_f") ;
        output("ci->ZSourceTerm") ;
        constraint("boundaryFaces,souFlamelet") ;
      }

      // Increment the source term for the cells attach to a single face. Note
      // that here we are insisting that the reconstructed value of Z at the
      // face be positive.
      void calculate(Entity face) {
        ZSourceTerm[ci[face]]-=(massFlux[face]>0.0)? massFlux[face]*
          max(ZLimiter[ci[face]]*dot(ZGradient[ci[face]],faceCenter[face]-
          cellCenter[ci[face]]),-Z[ci[face]]):massFlux[face]*Z_f[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SOUInviscidFluxToZSourceTermBoundary>
    registerSOUInviscidFluxToZSourceTermBoundary ;

  
  // Rule to add the diffusive flux contribution to the source term for
  // interior faces.
  class DiffusiveFluxToZSourceTermInterior : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> vol ;
      const_store<vect3d> ZGradient ;
      const_store<vect3d> geometryFactor0 ;
      const_store<real> faceRadius ;
      const_store<real> ZViscosity ;
      store<real> ZSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToZSourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("vol",vol) ;
        name_store("grads(Z)",ZGradient) ;
        name_store("geometryFactor0",geometryFactor0) ;
        name_store("faceRadius",faceRadius) ;
        name_store("ZViscosity",ZViscosity) ;
        name_store("ZSourceTerm",ZSourceTerm) ;
        input("(cl,cr)->(vol,grads(Z))") ;
        input("geometryFactor0,faceRadius,ZViscosity") ;
        output("(cl,cr)->ZSourceTerm") ;
        constraint("internalFaces") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real secondarySourceTerm=ZViscosity[face]*dot((ZGradient[cl[face]]*
          vol[cr[face]]+ZGradient[cr[face]]*vol[cl[face]]),
          geometryFactor0[face])*faceRadius[face]/(vol[cl[face]]+vol[cr[face]]);
        ZSourceTerm[cl[face]]+=secondarySourceTerm ;
        ZSourceTerm[cr[face]]-=secondarySourceTerm ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToZSourceTermInterior>
    registerDiffusiveFluxToZSourceTermInterior ;

  // Rule to add the diffusive flux contribution to the source term for
  // boundary faces.
  class DiffusiveFluxToZSourceTermBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> ZGradient ;
      const_store<vect3d> faceCenter ;
      const_store<real> diffusionProduct ;
      const_store<real> Z_f ;
      const_store<real> ZViscosity ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<real> ZSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToZSourceTermBoundary() {
        name_store("ci",ci) ;
		name_store("ZViscosity",ZViscosity) ;
        name_store("cellcenter",cellCenter) ;
        name_store("grads(Z)",ZGradient) ;
        name_store("facecenter",faceCenter) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("Z_f",Z_f) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("ZSourceTerm",ZSourceTerm) ;
        input("ci->cellcenter,ci->grads(Z)") ;
        input("facecenter,diffusionProduct,Z_f") ;
        input("ZViscosity,area,faceRadius") ;
        output("ci->ZSourceTerm") ;
        constraint("boundaryFlameletDiffusion") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real sourceTerm=ZViscosity[face]*
          (Z_f[face]*diffusionProduct[face]+dot(ZGradient[ci[face]],
          (area[face].n*area[face].sada-diffusionProduct[face]*
          (faceCenter[face]-cellCenter[ci[face]]))))*faceRadius[face] ;
        ZSourceTerm[ci[face]]+=sourceTerm ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToZSourceTermBoundary>
    registerDiffusiveFluxToZSourceTermBoundary ;

  // Rule to compute the diagonal term for the linear system.
  class ComputeZMatrixDiagonal : public pointwise_rule {
    private:
      const_param<real> ZRelaxationFactor ;
      const_store<real> ZMainCoefficient ;
      store<real> D ;
    public:

      // Define input and output.
      ComputeZMatrixDiagonal() {
        name_store("flameletRelaxationFactor",ZRelaxationFactor) ;
        name_store("ZMainCoefficient",ZMainCoefficient) ;
        name_store("ZStar_D",D) ;
        input("flameletRelaxationFactor,ZMainCoefficient") ;
        output("ZStar_D") ;
        constraint("geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) {
        D[cell]=ZMainCoefficient[cell]/(*ZRelaxationFactor) ;
      }

      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeZMatrixDiagonal> registerComputeZMatrixDiagonal ;

  
  // Rule to initialize the lower terms for the linear system.
  class InitializeZMatrixLower : public unit_rule {
    private:
      store<real> L ;
    public:

      // Define input and output.
      InitializeZMatrixLower() {
        name_store("ZStar_L",L) ;
        output("ZStar_L") ;
        constraint("internalFaces") ;
      }

      // Initialize for a single face.
      void calculate(Entity face) { L[face]=0.0 ; }

      // Initialize for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeZMatrixLower> registerInitializeZMatrixLower ;

  // Rule to add the first-order inviscid flux contribution to the lower terms
  // for the linear system.
  class FOUInviscidFluxToZMatrixLower : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> massFlux ;
      store<real> L ;
    public:

      // Define input and output.
      FOUInviscidFluxToZMatrixLower() {
        name_store("massFlux",massFlux) ;
        name_store("ZStar_L",L) ;
        input("massFlux") ;
        output("ZStar_L") ;
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

  register_rule<FOUInviscidFluxToZMatrixLower>
    registerFOUInviscidFluxToZMatrixLower ;

  // Rule to add the diffusive flux contribution to the lower terms for the
  // linear system.
  class DiffusiveFluxToZMatrixLower : public apply_rule<store<real>,Loci::
  Summation<real> > {
    private:
      const_store<real> diffusionProduct ;
      const_store<real> faceRadius ;
      const_store<real> ZViscosity ;
      store<real> L ;
    public:

      // Define input and output.
      DiffusiveFluxToZMatrixLower() {
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("faceRadius",faceRadius) ;
        name_store("ZViscosity",ZViscosity) ;
        name_store("ZStar_L",L) ;
        input("diffusionProduct,faceRadius,ZViscosity") ;
        output("ZStar_L") ;
        constraint("internalFaces,flameletModel") ;
      }

      // Increment the lower term for a single face.
      void calculate(Entity face) {
        L[face]-=ZViscosity[face]*diffusionProduct[face]*faceRadius[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToZMatrixLower>
    registerDiffusiveFluxToZMatrixLower ;

  // Rule to initialize the upper terms for the linear system.
  class InitializeZMatrixUpper : public unit_rule {
    private:
      store<real> U ;
    public:

      // Define input and output.
      InitializeZMatrixUpper() {
        name_store("ZStar_U",U) ;
        output("ZStar_U") ;
        constraint("internalFaces") ;
      }

      // Initialize for a single face.
      void calculate(Entity face) { U[face]=0.0 ; }

      // Initialize for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeZMatrixUpper> registerInitializeZMatrixUpper ;

  // Rule to add the first-order inviscid flux contribution to the upper terms
  // for the linear system.
  class FOUInviscidFluxToZMatrixUpper : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> massFlux ;
      store<real> U ;
    public:

      // Define input and output.
      FOUInviscidFluxToZMatrixUpper() {
        name_store("massFlux",massFlux) ;
        name_store("ZStar_U",U) ;
        input("massFlux") ;
        output("ZStar_U") ;
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

  register_rule<FOUInviscidFluxToZMatrixUpper>
    registerFOUInviscidFluxToZMatrixUpper ;

  // Rule to add the diffusive flux contribution to the upper terms for the
  // linear system.
  class DiffusiveFluxToZMatrixUpper : public apply_rule<store<real>,Loci::
  Summation<real> > {
    private:
      const_store<real> diffusionProduct ;
      const_store<real> faceRadius ;
      const_store<real> ZViscosity ;
      store<real> U ;
    public:

      // Define input and output.
      DiffusiveFluxToZMatrixUpper() {
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("faceRadius",faceRadius) ;
        name_store("ZViscosity",ZViscosity) ;
        name_store("ZStar_U",U) ;
        input("diffusionProduct,faceRadius,ZViscosity") ;
        output("ZStar_U") ;
        constraint("internalFaces,flameletModel") ;
      }

      // Increment the lower term for a single face.
      void calculate(Entity face) {
        U[face]-=ZViscosity[face]*diffusionProduct[face]*faceRadius[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToZMatrixUpper>
    registerDiffusiveFluxToZMatrixUpper ;

  // Rule to compute the right-hand side for the linear system.
  class ComputeZRHS : public pointwise_rule {
    private:
      const_param<real> ZRelaxationFactor ;
      const_store<real> Z ;
      const_store<real> ZMainCoefficient ;
      const_store<real> ZSourceTerm ;
      store<real> B ;
    public:

      // Define input and output.
      ComputeZRHS() {
        name_store("flameletRelaxationFactor",ZRelaxationFactor) ;
        name_store("Z",Z) ;
        name_store("ZMainCoefficient",ZMainCoefficient) ;
        name_store("ZSourceTerm",ZSourceTerm) ;
        name_store("ZStar_B",B) ;
        input("flameletRelaxationFactor,Z,ZMainCoefficient,ZSourceTerm") ;
        output("ZStar_B") ;
        constraint("geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) {
        B[cell]=ZSourceTerm[cell]+(1.0-(*ZRelaxationFactor))*
          ZMainCoefficient[cell]*Z[cell]/(*ZRelaxationFactor) ;
      }

      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeZRHS> registerComputeZRHS ;

//-----------------------------------------------------------------------------
// Rules for computing the residual of the Z equation.

  // Rule to initialize the Z residual.
  class InitializeZResidual : public unit_rule {
    private:
      store<real> ZResidual ;
    public:

      // Define input and output.
      InitializeZResidual() {
        name_store("ZResidual",ZResidual) ;
        output("ZResidual") ;
        constraint("geom_cells") ;
      }

      // Initialize the residual for a single cell.
      void calculate(Entity cell) { ZResidual[cell]=0.0 ; }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeZResidual> registerInitializeZResidual ;

  // Rule to compute the Z residual for each cell.
  class ComputeZResidualOne : public apply_rule<store<real>,Loci::Summation
  <real> > {
    private:
      const_store<real> D ;
      const_store<real> Z ;
      const_store<real> B ;
      store<real> ZResidual ;
    public:

      // Define input and output.
      ComputeZResidualOne() {
        name_store("ZStar_D",D) ;
        name_store("Z",Z) ;
        name_store("ZStar_B",B) ;
        name_store("ZResidual",ZResidual) ;
        input("ZStar_D,Z,ZStar_B") ;
        output("ZResidual") ;
        constraint("geom_cells") ;
      }

      // Initialize the residual for a single cell.
      void calculate(Entity cell) {
        ZResidual[cell]+=B[cell]-D[cell]*Z[cell] ;
      }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ;}
  } ;

  register_rule<ComputeZResidualOne> registerComputeZResidualOne ;

  // Rule to compute the Z residual for each cell.
  class ComputeZResidualTwo : public apply_rule<store<real>,Loci::Summation
  <real> > {
    private:
      const_Map cl,cr ;
      const_store<real> Z ;
      const_store<real> L,U ;
      store<real> ZResidual ;
    public:

      // Define input and output.
      ComputeZResidualTwo() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("Z",Z) ;
        name_store("ZStar_L",L) ;
        name_store("ZStar_U",U) ;
        name_store("ZResidual",ZResidual) ;
        input("(cl,cr)->Z,ZStar_L,ZStar_U") ;
        output("(cl,cr)->ZResidual") ;
        constraint("internalFaces") ;
      }

      // Add the neighbor contribution to the residual for each of the two
      // cells on either side of the face.
      void calculate(Entity face) {
        ZResidual[cl[face]]-= U[face]*Z[cr[face]] ;
        ZResidual[cr[face]]-= L[face]*Z[cl[face]] ;
      }

      // Add the neighbor contributions for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ;}
  } ;

  register_rule<ComputeZResidualTwo> registerComputeZResidualTwo ;


  // Rule to initialize the total Z residual.
  class InitializeTotalZResidual : public unit_rule {
    private:
      param<ScalarResidual> ZResidualData ;
    public:

      // Define input and output.
      InitializeTotalZResidual() {
        name_store("ZResidualData",ZResidualData) ;
        output("ZResidualData") ;
        constraint("flameletModel,geom_cells") ;
      }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) {
        ZResidualData=ScalarResidual() ;
      }
  } ;

  register_rule<InitializeTotalZResidual> registerInitializeTotalZResidual ;

  // Rule to compute the total Z residual.
  class ComputeTotalZResidual : public apply_rule<param<ScalarResidual>,
  ScalarResidualJoin> {
    private:
      const_store<real> ZResidual ;
      const_store<real> D ;
      const_store<vect3d> cellCenter ;
      param<ScalarResidual> ZResidualData ;
    public:

      // Define input and output.
      ComputeTotalZResidual() {
        name_store("ZResidual",ZResidual) ;
        name_store("cellcenter",cellCenter) ;
        name_store("ZResidualData",ZResidualData) ;
        input("ZResidual,cellcenter") ;
        output("ZResidualData") ;
        constraint("flameletModel,geom_cells") ;
      }

      // Add the cell contribution to the residual for a single cell.
      void calculate(Entity cell) {
        ScalarResidual temp ;
        temp.maxResidual=ZResidual[cell] ;
        temp.totalResidual=abs(ZResidual[cell]) ;
        temp.maxResidualLocation=cellCenter[cell] ;
        join(*ZResidualData,temp) ;
      }

      // Add the cell contribution to the residual for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeTotalZResidual> registerComputeTotalZResidual ;

//-----------------------------------------------------------------------------
// Scheme independent rules for assembling the Zvar equation.

  // Rule to initialize the main coefficient.
  class InitializeZvarMainCoefficient : public unit_rule {
    private:
      store<real> ZvarMainCoefficient ;
    public:

      // Define input and output.
     InitializeZvarMainCoefficient() {
        name_store("ZvarMainCoefficient",ZvarMainCoefficient) ;
	output("ZvarMainCoefficient") ;
        constraint("geom_cells") ;
      }

      // Set the main coefficient to zero for a single cell.
      void calculate(Entity cell) { ZvarMainCoefficient[cell]=0.0 ; }

      // Set the main coefficient to zero for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeZvarMainCoefficient>
    registerInitializeZvarMainCoefficient ;

  // Rule to add the first-order inviscid flux contribution to the main
  // coefficient for interior faces.
  class FOUInviscidFluxToZvarMainCoefficientInterior : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> massFlux ;
      store<real> ZvarMainCoefficient ;
    public:

      // Define input and output.
      FOUInviscidFluxToZvarMainCoefficientInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("massFlux",massFlux) ;
        name_store("ZvarMainCoefficient",ZvarMainCoefficient) ;
        input("massFlux") ;
        output("cl->ZvarMainCoefficient,cr->ZvarMainCoefficient") ;
        constraint("internalFaces") ;
      }

      // Increment the main coefficient for cells attached to a single face.
      void calculate(Entity face) {
        if(massFlux[face]>0.0){
          ZvarMainCoefficient[cr[face]]+=massFlux[face] ;
        }else{
          ZvarMainCoefficient[cl[face]]-=massFlux[face] ;
        }
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToZvarMainCoefficientInterior>
    registerFOUInviscidFluxToZvarMainCoefficientInterior ;

  // Rule to add the first-order inviscid flux contribution to the main
  // coefficient for boundary faces.
  class FOUInviscidFluxToZvarMainCoefficientBoundary : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> massFlux ;
      store<real> ZvarMainCoefficient ;
    public:

      // Define input and output.
      FOUInviscidFluxToZvarMainCoefficientBoundary() {
        name_store("ci",ci) ;
        name_store("massFlux",massFlux) ;
        name_store("ZvarMainCoefficient",ZvarMainCoefficient) ;
        input("massFlux") ;
        output("ci->ZvarMainCoefficient") ;
        constraint("boundaryFaces") ;
      }

      // Increment the main coefficient for the cell attach to a single face.
      void calculate(Entity face) {
        if(massFlux[face]<0.0) ZvarMainCoefficient[ci[face]]-=massFlux[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToZvarMainCoefficientBoundary>
    registerFOUInviscidFluxToZvarMainCoefficientBoundary ;

  // Rule to add the destruction term to the main coefficient.
  class DestructionToZvarMainCoefficient : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> Cd;
      const_param<real> betaStar ;
      const_store<real> rho ;
      const_store<real> omega ;
      const_store<real> Zvar ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> ZvarMainCoefficient ;
    public:

      // Define input and output.
      DestructionToZvarMainCoefficient() {
        name_store("Cd",Cd) ;
        name_store("betaStar",betaStar) ;
        name_store("rho",rho) ;
        name_store("omega",omega) ;
        name_store("Zvar",Zvar) ;
        name_store("vol",vol) ;
        name_store("cellRadius",cellRadius) ;
        name_store("ZvarMainCoefficient",ZvarMainCoefficient) ;
        input("Cd,betaStar,rho,omega,Zvar,vol,cellRadius") ;
        output("ZvarMainCoefficient") ;
        constraint("geom_cells,flameletModel") ;
      }

      // Increment the main coefficient for a single cell.
      void calculate(Entity cell) {
        ZvarMainCoefficient[cell]+=(*Cd)*rho[cell]*(*betaStar)*omega[cell]*vol[cell]*
          cellRadius[cell] ;
      }

      // Increment the main coefficient for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DestructionToZvarMainCoefficient>
    registerDestructionToZvarMainCoefficient ;

  // Rule to add the diffusive flux contribution to the main coefficient for
  // interior faces.
  class DiffusiveFluxToZvarMainCoefficientInterior : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> faceRadius ;
      const_store<real> diffusionProduct ;
      const_store<real> ZViscosity ;
      store<real> ZvarMainCoefficient ;
    public:

      // Define input and output.
      DiffusiveFluxToZvarMainCoefficientInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("faceRadius",faceRadius) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("ZViscosity",ZViscosity) ;
        name_store("ZvarMainCoefficient",ZvarMainCoefficient) ;
        input("faceRadius,diffusionProduct,ZViscosity") ;
        output("(cl,cr)->ZvarMainCoefficient") ;
        constraint("internalFaces") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real temp=ZViscosity[face]*diffusionProduct[face]*faceRadius[face] ;
        ZvarMainCoefficient[cl[face]]+=temp ;
        ZvarMainCoefficient[cr[face]]+=temp ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToZvarMainCoefficientInterior>
    registerDiffusiveFluxToZvarMainCoefficientInterior ;

  // Rule to add the diffusive flux contribution to the main coefficient for
  // boundary faces.
  class DiffusiveFluxToZvarMainCoefficientBoundary : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> faceRadius ;
      const_store<real> diffusionProduct ;
      const_store<real> ZViscosity ;
      store<real> ZvarMainCoefficient ;
    public:

      // Define input and output.
      DiffusiveFluxToZvarMainCoefficientBoundary() {
        name_store("ci",ci) ;
        name_store("faceRadius",faceRadius) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("ZViscosity",ZViscosity) ;
        name_store("ZvarMainCoefficient",ZvarMainCoefficient) ;
        input("ci->faceRadius") ;
        input("diffusionProduct,ZViscosity") ;
        output("ci->ZvarMainCoefficient") ;
        constraint("boundaryFlameletDiffusion") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        ZvarMainCoefficient[ci[face]]+=ZViscosity[face]*diffusionProduct[face]*faceRadius[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToZvarMainCoefficientBoundary>
    registerDiffusiveFluxToZvarMainCoefficientBoundary ;

  // Rule to initialize the source term.
  class InitializeZvarSourceTerm : public unit_rule {
    private:
      store<real> ZvarSourceTerm ;
    public:

      // Define input and output.
      InitializeZvarSourceTerm() {
        name_store("ZvarSourceTerm",ZvarSourceTerm) ;
        output("ZvarSourceTerm") ;
        constraint("geom_cells") ;
      }

      // Set the source term to zero for a single cell.
      void calculate(Entity cell) { ZvarSourceTerm[cell]=0.0 ; }

      // Set the source term to zero for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeZvarSourceTerm> registerInitializeZvarSourceTerm ;

  // Rule to add the first-order inviscid flux contribution to the source term
  // for boundary faces.
  class FOUInviscidFluxToZvarSourceTermBoundary : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> Zvar ;
      const_store<real> massFlux ;
      const_store<real> Zvar_f ;
      store<real> ZvarSourceTerm ;
    public:

      // Define input and output.
      FOUInviscidFluxToZvarSourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("Zvar",Zvar) ;
        name_store("massFlux",massFlux) ;
        name_store("Zvar_f",Zvar_f) ;
        name_store("ZvarSourceTerm",ZvarSourceTerm) ;
        input("ci->Zvar,massFlux,Zvar_f") ;
        output("ci->ZvarSourceTerm") ;
        constraint("boundaryFaces,fouFlamelet") ;
      }

      // Increment the source term for the cell attach to a single face.
      void calculate(Entity face) {
        if(massFlux[face]<=0.0) ZvarSourceTerm[ci[face]]-=massFlux[face]*Zvar_f[face] ;
        else ZvarSourceTerm[ci[face]]-=massFlux[face]*(Zvar_f[face]-Zvar[ci[face]]) ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToZvarSourceTermBoundary>
    registerFOUInviscidFluxToZvarSourceTermBoundary ;

  // Rule to add the second-order convection contribution to the source term
  // for interior faces.
  class SOUInviscidFluxToZvarSourceTermInterior : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> Zvar ;
      const_store<vect3d> ZvarGradient ;
      const_store<real> ZvarLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      store<real> ZvarSourceTerm ;
    public:

      // Define input and output.
      SOUInviscidFluxToZvarSourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("Zvar",Zvar) ;
        name_store("grads(Zvar)",ZvarGradient) ;
        name_store("limiters(Zvar)",ZvarLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("massFlux",massFlux) ;
        name_store("ZvarSourceTerm",ZvarSourceTerm) ;
        input("(cl,cr)->(Zvar,grads(Zvar),limiters(Zvar),cellcenter)") ;
        input("facecenter,massFlux") ;
        output("(cl,cr)->ZvarSourceTerm") ;
        constraint("internalFaces,souFlamelet") ;
      }

      // Increment the source term for the cells attach to a single face. Note
      // that here we are insisting that the reconstructed value of Zvar at the
      // face be positive.
      void calculate(Entity face) {
        real secondOrderSource=(massFlux[face]>0.0)? massFlux[face]*
          max(ZvarLimiter[cl[face]]*dot(ZvarGradient[cl[face]],
          faceCenter[face]-cellCenter[cl[face]]),-Zvar[cl[face]]):
          massFlux[face]*max(ZvarLimiter[cr[face]]*dot(ZvarGradient[cr[face]],
          faceCenter[face]-cellCenter[cr[face]]),-Zvar[cr[face]]) ;
        ZvarSourceTerm[cl[face]]-=secondOrderSource ;
        ZvarSourceTerm[cr[face]]+=secondOrderSource ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SOUInviscidFluxToZvarSourceTermInterior>
    registerSOUInviscidFluxToZvarSourceTermInterior ;

  // Rule to add the second-order convection contribution to the source term
  // for boundary faces.
  class SOUInviscidFluxToZvarSourceTermBoundary : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> Zvar ;
      const_store<vect3d> ZvarGradient ;
      const_store<real> ZvarLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      const_store<real> Zvar_f ;
      store<real> ZvarSourceTerm ;
    public:

      // Define input and output.
      SOUInviscidFluxToZvarSourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("Zvar",Zvar) ;
        name_store("grads(Zvar)",ZvarGradient) ;
        name_store("limiters(Zvar)",ZvarLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("massFlux",massFlux) ;
        name_store("Zvar_f",Zvar_f) ;
        name_store("ZvarSourceTerm",ZvarSourceTerm) ;
        input("ci->(Zvar,grads(Zvar),limiters(Zvar),cellcenter)") ;
        input("facecenter,massFlux,Zvar_f") ;
        output("ci->ZvarSourceTerm") ;
        constraint("boundaryFaces,souFlamelet") ;
      }

      // Increment the source term for the cells attach to a single face. Note
      // that here we are insisting that the reconstructed value of Zvar at the
      // face be positive.
      void calculate(Entity face) {
        ZvarSourceTerm[ci[face]]-=(massFlux[face]>0.0)? massFlux[face]*
          max(ZvarLimiter[ci[face]]*dot(ZvarGradient[ci[face]],
          faceCenter[face]-cellCenter[ci[face]]),-Zvar[ci[face]]):
          massFlux[face]*Zvar_f[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SOUInviscidFluxToZvarSourceTermBoundary>
    registerSOUInviscidFluxToZvarSourceTermBoundary ;

  // Rule to add the production term to the source term.
  class ProductionToZvarSourceTerm : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_param<real> Cg ;
      const_store<real> eddyViscosity ;
      const_store<vect3d> ZGradient ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> ZvarSourceTerm ;
    public:

      // Define input and output.
      ProductionToZvarSourceTerm() {
        name_store("Cg",Cg) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("grads(Z)",ZGradient) ;
        name_store("vol",vol) ;
        name_store("cellRadius",cellRadius) ;
        name_store("ZvarSourceTerm",ZvarSourceTerm) ;
        input("Cg") ;
        input("eddyViscosity") ;
        input("grads(Z),vol,cellRadius") ;
        output("ZvarSourceTerm") ;
        constraint("geom_cells,flameletModel") ;
      }

      // Increment the source term for a single cell. Note that kProduction
      // already has been multiplied by the volume.
      void calculate(Entity cell) {
        ZvarSourceTerm[cell]+=(*Cg)*eddyViscosity[cell]*dot(ZGradient[cell],ZGradient[cell])*vol[cell]*cellRadius[cell];
      }

      // Increment the source term for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ProductionToZvarSourceTerm>
    registerProductionToZvarSourceTerm ;

  // Rule to add the diffusive flux contribution to the source term for
  // interior faces.
  class DiffusiveFluxToZvarSourceTermInterior : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> vol ;
      const_store<vect3d> ZvarGradient ;
      const_store<vect3d> geometryFactor0 ;
      const_store<real> faceRadius ;
      const_store<real> ZViscosity ;
      store<real> ZvarSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToZvarSourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("vol",vol) ;
        name_store("grads(Zvar)",ZvarGradient) ;
        name_store("geometryFactor0",geometryFactor0) ;
        name_store("faceRadius",faceRadius) ;
        name_store("ZViscosity",ZViscosity) ;
        name_store("ZvarSourceTerm",ZvarSourceTerm) ;
        input("(cl,cr)->(vol,grads(Zvar))");
        input("geometryFactor0,faceRadius,ZViscosity") ;
        output("(cl,cr)->ZvarSourceTerm") ;
        constraint("internalFaces") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real secondarySourceTerm=ZViscosity[face]*
          dot((ZvarGradient[cl[face]]*vol[cr[face]]+ZvarGradient[cr[face]]*
          vol[cl[face]]),geometryFactor0[face])*faceRadius[face]/
          (vol[cl[face]]+vol[cr[face]]) ;
        ZvarSourceTerm[cl[face]]+=secondarySourceTerm ;
        ZvarSourceTerm[cr[face]]-=secondarySourceTerm ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToZvarSourceTermInterior>
    registerDiffusiveFluxToZvarSourceTermInterior ;

  // Rule to add the diffusive flux contribution to the source term for
  // boundary faces.
  class DiffusiveFluxToZvarSourceTermBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> ZvarGradient ;
      const_store<vect3d> faceCenter ;
      const_store<real> diffusionProduct ;
      const_store<real> Zvar_f ;
      const_store<real> ZViscosity ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<real> ZvarSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToZvarSourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("cellcenter",cellCenter) ;
        name_store("grads(Zvar)",ZvarGradient) ;
        name_store("facecenter",faceCenter) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("Zvar_f",Zvar_f) ;
        name_store("ZViscosity",ZViscosity) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("ZvarSourceTerm",ZvarSourceTerm) ;
        input("ci->cellcenter,ci->grads(Zvar)") ;
        input("facecenter,diffusionProduct,Zvar_f") ;
        input("ZViscosity,area,faceRadius") ;
        output("ci->ZvarSourceTerm") ;
        constraint("boundaryFlameletDiffusion") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real sourceTerm=ZViscosity[face]*(Zvar_f[face]*diffusionProduct[face]+dot(ZvarGradient
          [ci[face]],(area[face].n*area[face].sada-diffusionProduct[face]*
          (faceCenter[face]-cellCenter[ci[face]]))))*faceRadius[face] ;
        ZvarSourceTerm[ci[face]]+=sourceTerm ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToZvarSourceTermBoundary>
    registerDiffusiveFluxToZvarSourceTermBoundary ;

  // Rule to compute the diagonal term for the linear system.
  class ComputeZvarMatrixDiagonal : public pointwise_rule {
    private:
      const_param<real> ZvarRelaxationFactor ;
      const_store<real> ZvarMainCoefficient ;
      store<real> D ;
    public:

      // Define input and output.
      ComputeZvarMatrixDiagonal() {
        name_store("flameletRelaxationFactor",ZvarRelaxationFactor) ;
        name_store("ZvarMainCoefficient",ZvarMainCoefficient) ;
        name_store("ZvarStar_D",D) ;
        input("flameletRelaxationFactor,ZvarMainCoefficient") ;
        output("ZvarStar_D") ;
        constraint("geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) {
        D[cell]=ZvarMainCoefficient[cell]/(*ZvarRelaxationFactor) ;
      }

      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeZvarMatrixDiagonal> registerComputeZvarMatrixDiagonal ;

  // Rule to initialize the lower terms for the linear system.
  class InitializeZvarMatrixLower : public unit_rule {
    private:
      store<real> L ;
    public:

      // Define input and output.
      InitializeZvarMatrixLower() {
        name_store("ZvarStar_L",L) ;
        output("ZvarStar_L") ;
        constraint("internalFaces") ;
      }

      // Initialize for a single face.
      void calculate(Entity face) { L[face]=0.0 ; }

      // Initialize for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeZvarMatrixLower> registerInitializeZvarMatrixLower ;

  // Rule to add the first-order inviscid flux contribution to the lower terms
  // for the linear system.
  class FOUInviscidFluxToZvarMatrixLower : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> massFlux ;
      store<real> L ;
    public:

      // Define input and output.
      FOUInviscidFluxToZvarMatrixLower() {
        name_store("massFlux",massFlux) ;
        name_store("ZvarStar_L",L) ;
        input("massFlux") ;
        output("ZvarStar_L") ;
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

  register_rule<FOUInviscidFluxToZvarMatrixLower>
    registerFOUInviscidFluxToZvarMatrixLower ;

  // Rule to add the diffusive flux contribution to the lower terms for the
  // linear system.
  class DiffusiveFluxToZvarMatrixLower : public apply_rule<store<real>,Loci::
  Summation<real> > {
    private:
      const_store<real> diffusionProduct ;
      const_store<real> faceRadius ;
      const_store<real> ZViscosity ;
      store<real> L ;
    public:

      // Define input and output.
      DiffusiveFluxToZvarMatrixLower() {
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("faceRadius",faceRadius) ;
        name_store("ZViscosity",ZViscosity) ;
        name_store("ZvarStar_L",L) ;
        input("diffusionProduct,faceRadius,ZViscosity") ;
        output("ZvarStar_L") ;
        constraint("internalFaces,flameletModel") ;
      }

      // Increment the lower term for a single face.
      void calculate(Entity face) {
        L[face]-=ZViscosity[face]*diffusionProduct[face]*faceRadius[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToZvarMatrixLower>
    registerDiffusiveFluxToZvarMatrixLower ;

  // Rule to initialize the upper terms for the linear system.
  class InitializeZvarMatrixUpper : public unit_rule {
    private:
      store<real> U ;
    public:

      // Define input and output.
      InitializeZvarMatrixUpper() {
        name_store("ZvarStar_U",U) ;
        output("ZvarStar_U") ;
        constraint("internalFaces") ;
      }

      // Initialize for a single face.
      void calculate(Entity face) { U[face]=0.0 ; }

      // Initialize for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeZvarMatrixUpper> registerInitializeZvarMatrixUpper ;

  // Rule to add the first-order inviscid flux contribution to the upper terms
  // for the linear system.
  class FOUInviscidFluxToZvarMatrixUpper : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> massFlux ;
      store<real> U ;
    public:

      // Define input and output.
      FOUInviscidFluxToZvarMatrixUpper() {
        name_store("massFlux",massFlux) ;
        name_store("ZvarStar_U",U) ;
        input("massFlux") ;
        output("ZvarStar_U") ;
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

  register_rule<FOUInviscidFluxToZvarMatrixUpper>
    registerFOUInviscidFluxToZvarMatrixUpper ;

  // Rule to add the diffusive flux contribution to the lower terms for the
  // linear system.
  class DiffusiveFluxToZvarMatrixUpper : public apply_rule<store<real>,Loci::
  Summation<real> > {
    private:
      const_store<real> diffusionProduct ;
      const_store<real> faceRadius ;
      const_store<real> ZViscosity ;
      store<real> U ;
    public:

      // Define input and output.
      DiffusiveFluxToZvarMatrixUpper() {
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("faceRadius",faceRadius) ;
        name_store("ZViscosity",ZViscosity) ;
        name_store("ZvarStar_U",U) ;
        input("diffusionProduct,faceRadius,ZViscosity") ;
        output("ZvarStar_U") ;
        constraint("internalFaces,flameletModel") ;
      }

      // Increment the lower term for a single face.
      void calculate(Entity face) {
        U[face]-=ZViscosity[face]*diffusionProduct[face]*faceRadius[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToZvarMatrixUpper>
    registerDiffusiveFluxToZvarMatrixUpper ;

  // Rule to compute the right-hand side for the linear system.
  class ComputeZvarRHS : public pointwise_rule {
    private:
      const_param<real> ZvarRelaxationFactor ;
      const_store<real> Zvar ;
      const_store<real> ZvarMainCoefficient ;
      const_store<real> ZvarSourceTerm ;
      store<real> B ;
    public:

      // Define input and output.
      ComputeZvarRHS() {
        name_store("flameletRelaxationFactor",ZvarRelaxationFactor) ;
        name_store("Zvar",Zvar) ;
        name_store("ZvarMainCoefficient",ZvarMainCoefficient) ;
        name_store("ZvarSourceTerm",ZvarSourceTerm) ;
        name_store("ZvarStar_B",B) ;
        input("flameletRelaxationFactor,Zvar,ZvarMainCoefficient") ;
        input("ZvarSourceTerm") ;
        output("ZvarStar_B") ;
        constraint("geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) {
        B[cell]=ZvarSourceTerm[cell]+(1.0-(*ZvarRelaxationFactor))*
          ZvarMainCoefficient[cell]*Zvar[cell]/(*ZvarRelaxationFactor) ;
      }

      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeZvarRHS> registerComputeZvarRHS ;

//-----------------------------------------------------------------------------
// Rules for computing the residual of the Zvar equation.

  // Rule to initialize the Zvar residual.
  class InitializeZvarResidual : public unit_rule {
    private:
      store<real> ZvarResidual ;
    public:

      // Define input and output.
      InitializeZvarResidual() {
        name_store("ZvarResidual",ZvarResidual) ;
        output("ZvarResidual") ;
        constraint("flameletModel,geom_cells") ;
      }

      // Initialize the residual for a single cell.
      void calculate(Entity cell) { ZvarResidual[cell]=0.0 ; }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeZvarResidual> registerInitializeZvarResidual ;

  // Rule to compute the Zvar residual for each cell.
  class ComputeZvarResidualOne : public apply_rule<store<real>,Loci::Summation
  <real> > {
    private:
      const_store<real> D ;
      const_store<real> Zvar ;
      const_store<real> B ;
      store<real> ZvarResidual ;
    public:

      // Define input and output.
      ComputeZvarResidualOne() {
        name_store("ZvarStar_D",D) ;
        name_store("Zvar",Zvar) ;
        name_store("ZvarStar_B",B) ;
        name_store("ZvarResidual",ZvarResidual) ;
        input("ZvarStar_D,Zvar,ZvarStar_B") ;
        output("ZvarResidual") ;
        constraint("geom_cells") ;
      }

      // Initialize the residual for a single cell.
      void calculate(Entity cell) {
        ZvarResidual[cell]+=B[cell]-D[cell]*Zvar[cell];
      }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) {
        do_loop(seq,this) ;
      }
  } ;

  register_rule<ComputeZvarResidualOne> registerComputeZvarResidualOne ;

  // Rule to compute the Zvar residual for each cell.
  class ComputeZvarResidualTwo : public apply_rule<store<real>,Loci::Summation
  <real> > {
    private:
      const_Map cl,cr ;
      const_store<real> Zvar ;
      const_store<real> L,U ;
      store<real> ZvarResidual ;
    public:

      // Define input and output.
      ComputeZvarResidualTwo() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("Zvar",Zvar) ;
        name_store("ZvarStar_L",L) ;
        name_store("ZvarStar_U",U) ;
        name_store("ZvarResidual",ZvarResidual) ;
        input("(cl,cr)->Zvar,ZvarStar_L,ZvarStar_U") ;
        output("(cl,cr)->ZvarResidual") ;
        constraint("internalFaces") ;
      }

      // Add the neighbor contribution to the residual for each of the two
      // cells on either side of the face.
      void calculate(Entity face) {
        ZvarResidual[cl[face]]-= U[face]*Zvar[cr[face]];
        ZvarResidual[cr[face]]-= L[face]*Zvar[cl[face]];
      }

      // Add the neighbor contributions for a sequence of faces.
      virtual void compute(const sequence &seq) {
        do_loop(seq,this) ;
      }
  } ;

  register_rule<ComputeZvarResidualTwo> registerComputeZvarResidualTwo ;

  // Rule to initialize the total Zvar residual.
  class InitializeTotalZvarResidual : public unit_rule {
    private:
      param<ScalarResidual> ZvarResidualData ;
    public:

      // Define input and output.
      InitializeTotalZvarResidual() {
        name_store("ZvarResidualData",ZvarResidualData) ;
        output("ZvarResidualData") ;
        constraint("flameletModel,geom_cells") ;
      }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) {
        ZvarResidualData=ScalarResidual() ;
      }
  } ;

  register_rule<InitializeTotalZvarResidual>
    registerInitializeTotalZvarResidual ;

  // Rule to compute the total Zvar residual. Note the use of ZvarStar_D to
  // determine the wall-function cells to exclude from the sum.
  class ComputeTotalZvarResidual : public apply_rule<param<ScalarResidual>,
  ScalarResidualJoin> {
    private:
      const_store<real> ZvarResidual ;
      const_store<real> D ;
      const_store<vect3d> cellCenter ;
      param<ScalarResidual> ZvarResidualData ;
    public:

      // Define input and output.
      ComputeTotalZvarResidual() {
        name_store("ZvarResidual",ZvarResidual) ;
        name_store("ZvarStar_D",D) ;
        name_store("cellcenter",cellCenter) ;
        name_store("ZvarResidualData",ZvarResidualData) ;
        input("ZvarResidual,ZvarStar_D,cellcenter") ;
        output("ZvarResidualData") ;
        constraint("geom_cells") ;
      }

      // Add the cell contribution to the residual for a single cell.
      void calculate(Entity cell) {
        ScalarResidual temp ;
        temp.maxResidual=ZvarResidual[cell] ;
        temp.totalResidual=abs(ZvarResidual[cell]) ;
        temp.maxResidualLocation=cellCenter[cell] ;
        join(*ZvarResidualData,temp) ;
      }

      // Add the cell contribution to the residual for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeTotalZvarResidual> registerComputeTotalZvarResidual ;

//-----------------------------------------------------------------------------
// Rules to limit Z and Zvar.

  // Limit Z and Zvar to positive values.
  class LimitZZvar: public pointwise_rule {
    private:
      store<real> ZStar,ZvarStar ;
    public:

      // Define input and output.
      LimitZZvar() {
        name_store("ZStar",ZStar) ;
        name_store("ZvarStar",ZvarStar) ;
        input("ZStar,ZvarStar") ;
        output("ZStarLimited=ZStar,ZvarStarLimited=ZvarStar") ;
        constraint("geom_cells") ;
      }

      // Limit Z and Zvar for a single cell.
      void calculate(Entity cell) {
        if(ZStar[cell]<0.0) ZStar[cell]=0. ;
        if(ZStar[cell]>1.0) ZStar[cell]=1. ;
        if(ZvarStar[cell]<0.0) ZvarStar[cell]=0. ;
        if(ZvarStar[cell]>0.25) ZvarStar[cell]=0.25 ;
      }

      // Limit Z and Zvar for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<LimitZZvar> registerLimitZZvar ;

//-----------------------------------------------------------------------------
// Rules for marching the Z and Zvar equations. 

  // Time build rule for Z and Zvar when using BDF2 time integrator. Although
  // this rule sets k{n=-1} from Z_ic and Zvar{n=-1} from Zvar_ic, these
  // values are not really used since BDF is used on the first timestep for
  // non-restarts. The only purpose for this rule is to let Loci know that
  // there are two previous time-levels that need to be maintained for BDF2.
  class TimeBuildZZvarBDF2: public pointwise_rule {
    private:
      const_store<real> Z_ic,Zvar_ic ;
      store<real> Z,Zvar;
    public:

      // Define input and output.
      TimeBuildZZvarBDF2() {
        name_store("Z_ic",Z_ic) ;
        name_store("Zvar_ic",Zvar_ic) ;
        name_store("Z{n=-1}",Z) ;
        name_store("Zvar{n=-1}",Zvar) ;
        input("Z_ic,Zvar_ic") ;
        output("Z{n=-1},Zvar{n=-1}") ;
        constraint("geom_cells,flameletModel") ;
      }

      // Assign Z and Zvar for a single cell.
      void calculate(Entity cell) {
        Z[cell]=Z_ic[cell] ; Zvar[cell]=Zvar_ic[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildZZvarBDF2> registerTimeBuildZZvarBDF2 ;

  // Time build rule for Z and Zvar.
  class TimeBuildZZvar: public pointwise_rule {
    private:
      const_store<real> Z_ic,Zvar_ic ;
      store<real> ZTimeStepZero,ZvarTimeStepZero ;
    public:

      // Define input and output.
      TimeBuildZZvar() {
        name_store("Z_ic",Z_ic) ;
        name_store("Zvar_ic",Zvar_ic) ;
        name_store("Z{n=0}",ZTimeStepZero) ;
        name_store("Zvar{n=0}",ZvarTimeStepZero) ;
        input("Z_ic,Zvar_ic") ;
        output("Z{n=0},Zvar{n=0}") ;
        constraint("geom_cells,flameletModel") ;
      }

      // Assign Z and Zvar at time zero for a single cell.
      void calculate(Entity cell) {
        ZTimeStepZero[cell]=Z_ic[cell] ;
        ZvarTimeStepZero[cell]=Zvar_ic[cell] ;
      }

      // Assign Z and Zvar at time zero for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildZZvar> registerTimeBuildZZvar ;

  // Iteration build rule for Z and Zvar.
  class IterationBuildZZvar : public pointwise_rule {
    private:
      const_store<real> ZTimeStepN,ZvarTimeStepN ;
      store<real> ZIterationZero,ZvarIterationZero ;
    public:

      // Define input and output.
      IterationBuildZZvar() {
        name_store("Z{n}",ZTimeStepN) ;
        name_store("Zvar{n}",ZvarTimeStepN) ;
        name_store("Z{n,it=0}",ZIterationZero) ;
        name_store("Zvar{n,it=0}",ZvarIterationZero) ;
        input("Z{n},Zvar{n}") ;
        output("Z{n,it=0},Zvar{n,it=0}") ;
        constraint("geom_cells{n},flameletModel{n}") ;
      }

      // Assign Z and Zvar at iteration zero for a single cell.
      void calculate(Entity cell) {
        ZIterationZero[cell]=ZTimeStepN[cell] ;
        ZvarIterationZero[cell]=ZvarTimeStepN[cell] ;
      }

      // Assign Z and Zvar at iteration zero for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationBuildZZvar> registerIterationBuildZZvar ;

  // Rule to add temporal component of the Z equation to the main coefficient.
  class TemporalToZMainCoefficient: public apply_rule<store<real>,Loci::
  Summation<real> > {
    private:
      const_param<real> timeIntegratorFactor0 ;
      const_store<real> timeStepFactor ;
      const_store<real> rho ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> ZMainCoefficient ;
    public:

      // Define input and output.
      TemporalToZMainCoefficient() {
        name_store("timeIntegratorFactor0{n}",timeIntegratorFactor0) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n}",rho) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("ZMainCoefficient{n,it}",ZMainCoefficient) ;
        input("rho{n},vol{n},cellRadius{n,it}") ;
        input("timeStepFactor{n},timeIntegratorFactor0{n}") ;
        output("ZMainCoefficient{n,it}") ;
        constraint("geom_cells") ;
      }

      // Add temporal component for a single cell.
      void calculate(Entity cell) {
        ZMainCoefficient[cell]+=rho[cell]*vol[cell]*cellRadius[cell]*
          (*timeIntegratorFactor0)/timeStepFactor[cell] ;
      }

      // Add temporal component for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToZMainCoefficient> registerTemporalToZMainCoefficient ;

  // Rule to add temporal component of the Zvar equation to the main
  // coefficient.
  class TemporalToZvarMainCoefficient: public apply_rule<store<real>,Loci::
  Summation<real> > {
    private:
      const_param<real> timeIntegratorFactor0 ;
      const_store<real> timeStepFactor ;
      const_store<real> rho ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> ZvarMainCoefficient ;
    public:

      // Define input and output.
      TemporalToZvarMainCoefficient() {
        name_store("timeIntegratorFactor0{n}",timeIntegratorFactor0) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n}",rho) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("ZvarMainCoefficient{n,it}",ZvarMainCoefficient) ;
        input("rho{n},vol{n},cellRadius{n,it}") ;
        input("timeStepFactor{n},timeIntegratorFactor0{n}") ;
        output("ZvarMainCoefficient{n,it}") ;
        constraint("geom_cells") ;
      }

      // Add temporal component for a single cell.
      void calculate(Entity cell) {
        ZvarMainCoefficient[cell]+=rho[cell]*vol[cell]*cellRadius[cell]*
          (*timeIntegratorFactor0)/timeStepFactor[cell] ;
      }

      // Add temporal component for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToZvarMainCoefficient>
    registerTemporalToZvarMainCoefficient ;

  // Rule to add temporal component of the k equation to the source term.
  class TemporalToZSourceTerm : public apply_rule<store<real>,Loci::
  Summation<real> > {
    private:
      const_param<real> timeIntegratorFactor0 ;
      const_store<real> timeStepFactor ;
      const_store<real> rho ;
      const_store<real> Z ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> ZSourceTerm ;
    public:

      // Define input and output.
      TemporalToZSourceTerm() {
        name_store("timeIntegratorFactor0{n}",timeIntegratorFactor0) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n}",rho) ;
        name_store("Z{n}",Z) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("ZSourceTerm{n,it}",ZSourceTerm) ;
        input("rho{n},Z{n},vol{n},cellRadius{n,it}") ;
        input("timeStepFactor{n},timeIntegratorFactor0{n}") ;
        output("ZSourceTerm{n,it}") ;
        constraint("geom_cells") ;
      }

      // Add temporal component to source term for a single cell.
      void calculate(Entity cell) {
        ZSourceTerm[cell]+=(rho[cell]*vol[cell]*cellRadius[cell]*
          (*timeIntegratorFactor0)/timeStepFactor[cell])*Z[cell] ;
      }

      // Add temporal component to source term for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToZSourceTerm> registerTemporalToZSourceTerm ;

  // Rule to add temporal component of momentum equation to the source term
  // for the BDF2 scheme.
  class TemporalToZSourceTermBDF2 : public apply_rule<store<real>,Loci::
  Summation<real> > {
    private:
      const_param<real> timeIntegratorFactor1 ;
      const_store<real> timeStepFactor ;
      const_store<real> rhoOld ;
      const_store<real> ZOld,Z ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> ZSourceTerm ;
    public:

      // Define input and output.
      TemporalToZSourceTermBDF2() {
        name_store("timeIntegratorFactor1{n}",timeIntegratorFactor1) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n-1}",rhoOld) ;
        name_store("Z{n-1}",ZOld) ;
        name_store("Z{n,it}",Z) ;
        name_store("vol{n,it}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("ZSourceTerm{n,it}",ZSourceTerm) ;
        input("rho{n-1},Z{n-1},Z{n,it},vol{n,it}") ;
        input("cellRadius{n,it},timeStepFactor{n},timeIntegratorFactor1{n}") ;
        output("ZSourceTerm{n,it}") ;
        constraint("geom_cells,BDF2Integrator") ;
      }

      // Add temporal component to source term for a single cell.
      void calculate(Entity cell) {
        ZSourceTerm[cell]+=((*timeIntegratorFactor1)*rhoOld[cell]*
          vol[cell]*cellRadius[cell]/timeStepFactor[cell])*
          (Z[cell]-ZOld[cell]) ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToZSourceTermBDF2> registerTemporalToZSourceTermBDF2 ;

  // Rule to add temporal component of the Zvar equation to the source term.
  class TemporalToZvarSourceTerm : public apply_rule<store<real>,Loci::
  Summation<real> > {
    private:
      const_param<real> timeIntegratorFactor0 ;
      const_store<real> timeStepFactor ;
      const_store<real> rho ;
      const_store<real> Zvar ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> ZvarSourceTerm ;
    public:

      // Define input and output.
      TemporalToZvarSourceTerm() {
        name_store("timeIntegratorFactor0{n}",timeIntegratorFactor0) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n}",rho) ;
        name_store("Zvar{n}",Zvar) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("ZvarSourceTerm{n,it}",ZvarSourceTerm) ;
        input("rho{n},Zvar{n},vol{n},cellRadius{n,it}") ;
        input("timeStepFactor{n},timeIntegratorFactor0{n}") ;
        output("ZvarSourceTerm{n,it}") ;
        constraint("geom_cells") ;
      }

      // Add temporal component to source term for a single cell.
      void calculate(Entity cell) {
        ZvarSourceTerm[cell]+=(rho[cell]*vol[cell]*cellRadius[cell]*
          (*timeIntegratorFactor0)/timeStepFactor[cell])*Zvar[cell] ;
      }

      // Add temporal component to source term for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToZvarSourceTerm> registerTemporalToZvarSourceTerm ;

  // Rule to add temporal component of momentum equation to the source term
  // for the BDF2 scheme.
  class TemporalToZvarSourceTermBDF2 : public apply_rule<store<real>,Loci::
  Summation<real> > {
    private:
      const_param<real> timeIntegratorFactor1 ;
      const_store<real> timeStepFactor ;
      const_store<real> rhoOld ;
      const_store<real> ZvarOld,Zvar ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> ZvarSourceTerm ;
    public:

      // Define input and output.
      TemporalToZvarSourceTermBDF2() {
        name_store("timeIntegratorFactor1{n}",timeIntegratorFactor1) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n-1}",rhoOld) ;
        name_store("Zvar{n-1}",ZvarOld) ;
        name_store("Zvar{n,it}",Zvar) ;
        name_store("vol{n,it}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("ZvarSourceTerm{n,it}",ZvarSourceTerm) ;
        input("rho{n-1},Zvar{n-1},Zvar{n,it},vol{n,it}") ;
        input("cellRadius{n,it},timeStepFactor{n},timeIntegratorFactor1{n}") ;
        output("ZvarSourceTerm{n,it}") ;
        constraint("geom_cells,BDF2Integrator") ;
      }

      // Add temporal component to source term for a single cell.
      void calculate(Entity cell) {
        ZvarSourceTerm[cell]+=((*timeIntegratorFactor1)*rhoOld[cell]*
          vol[cell]*cellRadius[cell]/timeStepFactor[cell])*
          (Zvar[cell]-ZvarOld[cell]) ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToZvarSourceTermBDF2>
    registerTemporalToZvarSourceTermBDF2 ;

  // Iteration advance rule for Z and Zvar. 
  class IterationAdvanceZZvar : public pointwise_rule {
    private:
      const_store<real> ZStar,ZvarStar ;
      store<real> ZIterationPlusOne,ZvarIterationPlusOne ;
    public:

      // Define input and output.
      IterationAdvanceZZvar() {
        name_store("ZStarLimited{n,it}",ZStar) ;
        name_store("ZvarStarLimited{n,it}",ZvarStar) ;
        name_store("Z{n,it+1}",ZIterationPlusOne) ;
        name_store("Zvar{n,it+1}",ZvarIterationPlusOne) ;
        input("ZStarLimited{n,it},ZvarStarLimited{n,it}") ;
        output("Z{n,it+1},Zvar{n,it+1}") ;
        constraint("geom_cells{n,it},flameletModel{n,it}") ;
      }

      // Assign Z and Zvar at end of iteration for a single cell.
      void calculate(Entity cell) {
        ZIterationPlusOne[cell]=ZStar[cell] ;
        ZvarIterationPlusOne[cell]=ZvarStar[cell] ;
      }

      // Assign Z and Zvar at end of iteration for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationAdvanceZZvar> registerIterationAdvanceZZvar ;

  // Iteration collapse rule for Z and Zvar.
  class IterationCollapseZZvar : public pointwise_rule {
    private:
      const_param<bool> iterationFinished ;
      store<real> Z,Zvar ;
    public:

      // Define input and output.
      IterationCollapseZZvar() {
        name_store("iterationFinished{n,it-1}",iterationFinished) ;
        name_store("Z{n,it}",Z) ;
        name_store("Zvar{n,it}",Zvar) ;
        input("iterationFinished{n,it-1}") ;
        input("Z{n,it},Zvar{n,it}") ;
        output("Z{n+1}=Z{n,it},Zvar{n+1}=Zvar{n,it}") ;
        conditional("iterationFinished{n,it-1}") ;
        constraint("flameletModel{n,it},geom_cells{n,it}") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<IterationCollapseZZvar> registerIterationCollapseZZvar ;

   // Iteration advance rule for density
  class IterationAdvanceDensityFlamelet : public pointwise_rule {
    private:
      const_store<real> rhoStar ;
      store<real> rho ;
    public:

      // Define input and output.
      IterationAdvanceDensityFlamelet() {
        name_store("rhoStar{n,it}",rhoStar) ;
        name_store("flamelet::rho{n,it+1}",rho) ;
        input("rhoStar{n,it}") ;
        output("flamelet::rho{n,it+1}") ;
        constraint("geom_cells{n,it}") ;
        constraint("flameletModel{n,it}") ;
      }

      // Assign density at end of iteration for a single cell.
      void calculate(Entity cell) { rho[cell]=rhoStar[cell] ; }

      // Assign density at end of iteration for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationAdvanceDensityFlamelet>
    registerIterationAdvanceDensityFlamelet ;
  
}

