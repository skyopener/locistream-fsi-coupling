//-----------------------------------------------------------------------------
// Description: This file contains rules for assembling and solving the
//   pressure correction equation, as well as rules for correcting the
//   velocity, pressure and mass flux.
//
// Authors: Siddharth Thakur and Jeff Wright
//          Streamline Numerics Inc.
//          Gainesville, Florida 32609
//-----------------------------------------------------------------------------

// Standard library includes.
#include <fstream>

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;

// StreamUns includes.
#include "residual.h"
#include "sciTypes.h"
#include "varsFileInputs.h"

namespace streamUns {

//-----------------------------------------------------------------------------
// Rule to extract pPrime from PETSC. This is here now since we are using Ed's
// version of the PETSC solver in the FVM module.

  // Had a disable_threading() in this rule, but that was wrong. This still has
  // not fixed the bad behaviour after we have switched to Ed's PETSC solver.
  class PETSCCopyPPrime : public pointwise_rule {
    private:
      const_store<real> petscPPrime ;
      store<real> pPrime ;
    public:
                                                                                
      // Define input and output.
      PETSCCopyPPrime() {
        name_store("petscScalarSolve(pPrime)",petscPPrime) ;
        name_store("pPrime",pPrime) ;
        input("petscScalarSolve(pPrime)") ;
        output("pPrime") ;
        constraint("geom_cells,pPrime_PETSCLinearSolver") ;
      }
                                                                                
      // Copy the solution back from PETSC.
      void compute(const sequence & seq) {}
  } ;
                                                                                
//  register_rule<PETSCCopyPPrime> registerPETSCCopyPPrime ;

//-----------------------------------------------------------------------------
// Rules to process pressure correction equation options from the .vars file.

  // Creates the momentum equation solver constraints.
  class PressureEquationSolverConstraints : public constraint_rule {
    private:
      const_param<PressureEquationOptions> pressureEquationOptions ;
      Constraint pSGSLinearSolver,pPETSCLinearSolver,pHYPRELinearSolver ;
      Constraint compressiblePPrime ;
    public:
                                                                                
      // Define input and output.
      PressureEquationSolverConstraints() {
        name_store("pressureEquationOptions",pressureEquationOptions) ;
        name_store("pPrime_SGSLinearSolver",pSGSLinearSolver) ;
        name_store("pPrime_PETSCLinearSolver",pPETSCLinearSolver) ;
        name_store("pPrime_HYPRELinearSolver",pHYPRELinearSolver) ;
        name_store("compressiblePPrime",compressiblePPrime) ;
        input("pressureEquationOptions") ;
        output("pPrime_SGSLinearSolver,pPrime_PETSCLinearSolver") ;
        output("pPrime_HYPRELinearSolver,compressiblePPrime") ;
      }
                                                                                
      // Set up the constraints.
      virtual void compute(const sequence& seq) {

        // Linear solver.
        if((*pressureEquationOptions).optionExists("linearSolver")){
          Loci::option_value_type optionValueType=pressureEquationOptions->
            getOptionValueType("linearSolver") ;
          switch(optionValueType){
            case Loci::NAME:
              {
                Loci::option_values optionValues=pressureEquationOptions->
                  getOption("linearSolver") ;
                string name ; optionValues.get_value(name) ;
                if(name=="SGS"){
                  pSGSLinearSolver=~EMPTY ; pPETSCLinearSolver=EMPTY ;
                  pHYPRELinearSolver=EMPTY ;
                }else if(name=="PETSC"){
                  pSGSLinearSolver=EMPTY ; pPETSCLinearSolver=~EMPTY ;
                  pHYPRELinearSolver=EMPTY ;
                }else if(name=="HYPRE"){
                  pSGSLinearSolver=EMPTY ; pPETSCLinearSolver=EMPTY ;
                  pHYPRELinearSolver=~EMPTY ;
                }else{
                  cerr << "Bad linearSolver for pressureEquation." << endl ;
                  Loci::Abort() ;
                }
              }
              break ;
            default:
              cerr << "Bad type for linearSolver in pressureEquation." << endl ;
              Loci::Abort() ;
          }
        }else{
          pSGSLinearSolver=~EMPTY ; pPETSCLinearSolver=EMPTY ;
          pHYPRELinearSolver=EMPTY ;
        }

        // Compressible terms.
        if((*pressureEquationOptions).optionExists("incompressibleForm")){
          compressiblePPrime=EMPTY ;
        }else{
          compressiblePPrime=~EMPTY ;
        }
      }
  } ;
                                                                                
  register_rule<PressureEquationSolverConstraints>
    registerPressureEquationSolverConstraints ;

  // Creates the pressure equation solver maximum iterations parameter.
  class PressureEquationSolverParameters : public singleton_rule {
    private:
      const_param<PressureEquationOptions> pressureEquationOptions ;
      param<int> pMaxIterations ;
      param<real> pRelaxationFactor ;
    public:

      // Define input and output.
      PressureEquationSolverParameters() {
        name_store("pressureEquationOptions",pressureEquationOptions) ;
        name_store("pPrime_maxLinearSolverIterations",pMaxIterations) ;
        name_store("pRelaxationFactor",pRelaxationFactor) ;
        input("pressureEquationOptions") ;
        output("pPrime_maxLinearSolverIterations,pRelaxationFactor") ;
      }

      // Set up the parameter.
      virtual void compute(const sequence& seq) {

        // Maximum number of iterations.
        if((*pressureEquationOptions).optionExists("maxIterations")){
          Loci::option_value_type optionValueType=pressureEquationOptions->
            getOptionValueType("maxIterations") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                pressureEquationOptions->getOption("maxIterations",temp) ;
                if(int(temp)<0){
                  cerr << "Bad maxIterations value for pressureEquation."
                    << endl ; Loci::Abort() ;
                }
                *pMaxIterations=int(temp) ;
              }
              break ;
            default:
              cerr << "Bad type for maxIterations in pressureEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *pMaxIterations=5 ;
        }

        // Relaxation factor.
        if((*pressureEquationOptions).optionExists("relaxationFactor")){
          Loci::option_value_type optionValueType=pressureEquationOptions->
            getOptionValueType("relaxationFactor") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                pressureEquationOptions->getOption("relaxationFactor",temp) ;
                if(temp<=0.0 || temp>1.0){
                  cerr << "Bad relaxationFactor for pressureEquation." << endl ;
                  Loci::Abort() ;
                }
                *pRelaxationFactor=temp ;
              }
              break ;
            default:
              cerr << "Bad type for relaxationFactor in pressureEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *pRelaxationFactor=0.5 ;
        }
      }
  } ;
                                                                                
  register_rule<PressureEquationSolverParameters>
    registerPressureEquationSolverParamters ;

//-----------------------------------------------------------------------------
// Scheme independent rules for assembling the pressure-correction equation.

  // Rule to compute the pressure-correction coefficient on interior faces for
  // the FOU and SOU convection schemes.
  class PPrimeCoefficientInteriorFOUSOU : public pointwise_rule {
    private:
      const_param<real> vRelaxationFactor ;
      const_param<real> thetaParameter ;
      const_Map cl,cr ;
      const_store<real> vMainCoefficient ;
      const_store<real> vol ;
      const_store<real> faceDensity ;
      const_store<real> diffusionProduct ;
      const_store<real> faceRadius ;
      store<real> pPrimeCoefficient ;
    public:

      // Define input and output.
      PPrimeCoefficientInteriorFOUSOU() {
        name_store("vRelaxationFactor",vRelaxationFactor) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        name_store("vol",vol) ;
        name_store("faceDensity",faceDensity) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("faceRadius",faceRadius) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        input("thetaParameter,vRelaxationFactor,faceRadius") ;
        input("(cl,cr)->(vMainCoefficient,vol),faceDensity,diffusionProduct") ;
        output("pPrimeCoefficient") ;
        constraint("internalFaces,fouOrSouInviscidFlux") ;
      }

      // Calculate the pressure correction coefficient for a single face.
      void calculate(Entity face) {
        pPrimeCoefficient[face]=0.5*faceDensity[face]*(*vRelaxationFactor)*
          (vol[cl[face]]/vMainCoefficient[cl[face]]+vol[cr[face]]/
          vMainCoefficient[cr[face]])*diffusionProduct[face]*(*thetaParameter)*
          faceRadius[face]*faceRadius[face] ;
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<PPrimeCoefficientInteriorFOUSOU>
    registerPPrimeCoefficientInteriorFOUSOU ;

  // Rule to compute the pressure-correction coefficient on interior faces for
  // the Roe scheme.
  class PPrimeCoefficientInteriorRoe : public pointwise_rule {
    private:
      const_param<real> vRelaxationFactor ;
      const_Map cl,cr ;
      const_store<real> vMainCoefficient ;
      const_store<real> vol ;
      const_store<real> soundSpeed ;
      const_store<real> faceDensity ;
      const_store<real> faceMachNumber ;
      const_store<real> cStar ;
      const_store<real> vNormalTilde ;
      const_store<real> diffusionProduct ;
      const_store<Area> area ;
      store<real> pPrimeCoefficient ;
    public:

      // Define input and output.
      PPrimeCoefficientInteriorRoe() {
        name_store("vRelaxationFactor",vRelaxationFactor) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        name_store("vol",vol) ;
        name_store("soundSpeed",soundSpeed) ;
        name_store("faceDensity",faceDensity) ;
        name_store("faceMachNumber",faceMachNumber) ;
        name_store("cStar",cStar) ;
        name_store("vNormalTilde",vNormalTilde) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("area",area) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        input("vRelaxationFactor") ;
        input("(cl,cr)->(vMainCoefficient,vol,soundSpeed)") ;
        input("faceDensity,faceMachNumber,cStar,vNormalTilde") ;
        input("diffusionProduct,area") ;
        output("pPrimeCoefficient") ;
        constraint("internalFaces,roeInviscidFlux") ;
      }

      // Calculate the pressure correction coefficient for a single face.
      void calculate(Entity face) {
        if(faceMachNumber[face]>0.3){
          real cFace=0.5*(soundSpeed[cl[face]]+soundSpeed[cr[face]]) ;
          pPrimeCoefficient[face]=0.5*area[face].sada*(cStar[face]-
            abs(vNormalTilde[face]))/(cFace*cFace) ;
        }else{
          pPrimeCoefficient[face]=0.5*faceDensity[face]*(*vRelaxationFactor)*
          (vol[cl[face]]/vMainCoefficient[cl[face]]+vol[cr[face]]/
          vMainCoefficient[cr[face]])*diffusionProduct[face] ;
        }
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<PPrimeCoefficientInteriorRoe>
    registerPPrimeCoefficientInteriorRoe ;

  // Rule to compute the pressure-correction coefficient on boundary faces for
  // the FOU and SOU convection schemes.
  class PPrimeCoefficientBoundary : public pointwise_rule {
    private:
      const_param<real> vRelaxationFactor ;
      const_param<real> thetaParameter ;
      const_Map ci ;
      const_store<real> vMainCoefficient ;
      const_store<real> vol ;
      const_store<real> rho_f ;
      const_store<real> diffusionProduct ;
      const_store<real> faceRadius ;
      store<real> pPrimeCoefficient ;
    public:
                                                                                
      // Define input and output.
      PPrimeCoefficientBoundary() {
        name_store("vRelaxationFactor",vRelaxationFactor) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("ci",ci) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        name_store("vol",vol) ;
        name_store("rho_f",rho_f) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("faceRadius",faceRadius) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        input("vRelaxationFactor,thetaParameter,ci->(vMainCoefficient,vol)") ;
        input("rho_f,diffusionProduct,faceRadius") ;
        output("pPrimeCoefficient") ;
        constraint("boundaryFaces") ;
      }
                                                                                
      // Calculate the pressure-correction coefficient for a single face.
      void calculate(Entity face) {
        pPrimeCoefficient[face]=rho_f[face]*(*vRelaxationFactor)*
          vol[ci[face]]/vMainCoefficient[ci[face]]*diffusionProduct[face]*
          (*thetaParameter)*faceRadius[face]*faceRadius[face] ;
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;
                                                                                
  register_rule<PPrimeCoefficientBoundary> registerPPrimeCoefficientBoundary ;

  // Rule to initialize the diagonal term for the linear system. Checked.
  class InitializePressureCorrectionDiagonal : public unit_rule {
    private:
      store<real> D ;
    public:

      // Define input and output.
      InitializePressureCorrectionDiagonal() {
        name_store("pPrime_D",D) ;
        output("pPrime_D") ;
        constraint("vol") ;
      }

      // Initialize for a single cell.
      void calculate(Entity cell) { D[cell]=0.0 ; }

      // Initialize for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializePressureCorrectionDiagonal>
    registerInitializePressureCorrectionDiagonal ;

  // Rule to assemble the diagonal term for the linear system for
  // incompressible flow. Assembling over internal faces. Checked.
  class PressureCorrectionDiagonalIncompressibleInternal : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> pPrimeCoefficient ;
      store<real> D ;
    public:

      // Define input and output.
      PressureCorrectionDiagonalIncompressibleInternal() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pPrime_D",D) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        input("pPrimeCoefficient") ;
        output("(cl,cr)->pPrime_D") ;
        constraint("internalFaces") ;
      }

      // Distribute contribution for a single face.
      void calculate(Entity face) {
        D[cl[face]]+=pPrimeCoefficient[face] ;
        D[cr[face]]+=pPrimeCoefficient[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PressureCorrectionDiagonalIncompressibleInternal>
    registerPressureCorrectionDiagonalIncompressibleInternal ;

  // Rule to assemble the diagonal term for the linear system. Assembling over
  // boundary faces. Note that this rule now applies for both incompressible
  // and compressible flows. NOTE: Removed the nonSpecifiedMassFlux_BC constraint
  // so that this rule will be used for interpolated interface boundaries as
  // well where the pressure is specified via interpolation.
  class PressureCorrectionDiagonalBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> pPrimeCoefficient ;
      store<real> D ;
    public:

      // Define input and output.
      PressureCorrectionDiagonalBoundary() {
        name_store("ci",ci) ;
        name_store("pPrime_D",D) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        input("pPrimeCoefficient") ;
        output("ci->pPrime_D") ;
        constraint("specifiedPressure_BC") ;
      }

      // Distribute contribution for a single face.
      void calculate(Entity face) { D[ci[face]]+=pPrimeCoefficient[face] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PressureCorrectionDiagonalBoundary>
    registerPressureCorrectionDiagonalBoundary ;

  // Rule to assemble the diagonal term for the linear system for
  // compressible flow. Assembling over internal faces. Replaced "massFlux"
  // with "stageOneMassFlux" on 01/07/05.
  class PressureCorrectionDiagonalCompressibleInternal : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> rho ;
      const_store<real> soundSpeed ;
      const_store<real> massFlux ;
      store<real> D ;
    public:

      // Define input and output.
      PressureCorrectionDiagonalCompressibleInternal() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("rho",rho) ;
        name_store("soundSpeed",soundSpeed) ;
        name_store("pPrime_D",D) ;
        name_store("stageOneMassFlux",massFlux) ;
        input("(cl,cr)->rho,(cl,cr)->soundSpeed") ;
        input("stageOneMassFlux") ;
        output("(cl,cr)->pPrime_D") ;
        constraint("internalFaces,compressibleFlow,compressiblePPrime") ;
      }

      // Distribute contribution for a single face.
      void calculate(Entity face) {
        if(massFlux[face]>0.0){
          real compressibleCoefficient=massFlux[face]/(rho[cl[face]]*
            soundSpeed[cl[face]]*soundSpeed[cl[face]]) ;
          D[cl[face]]+=compressibleCoefficient ;
        }else{
          real compressibleCoefficient=-massFlux[face]/(rho[cr[face]]*
            soundSpeed[cr[face]]*soundSpeed[cr[face]]) ;
          D[cr[face]]+=compressibleCoefficient ;
        }
      }

      // Distribute contribution for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PressureCorrectionDiagonalCompressibleInternal>
    registerPressureCorrectionDiagonalCompressibleInternal ;

  // Rule to assemble the diagonal term for the linear system for
  // compressible flow. Assembling over boundary faces. Replaced "massFlux"
  // with "stageOneMassFlux" on 01/07/05.
  class PressureCorrectionDiagonalCompressibleBoundary : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> rho ;
      const_store<real> soundSpeed ;
      const_store<real> massFlux ;
      store<real> D ;
    public:

      // Define input and output.
      PressureCorrectionDiagonalCompressibleBoundary() {
        name_store("ci",ci) ;
        name_store("rho",rho) ;
        name_store("soundSpeed",soundSpeed) ;
        name_store("stageOneMassFlux",massFlux) ;
        name_store("pPrime_D",D) ;
        input("ci->rho,ci->soundSpeed,stageOneMassFlux") ;
        output("ci->pPrime_D") ;
        constraint("nonSpecifiedMassFlux_BC,compressibleFlow") ;
        constraint("compressiblePPrime") ;
      }

      // Distribute contribution for a single face.
      void calculate(Entity face) {
        if(massFlux[face]>0.0){
          D[ci[face]]+=massFlux[face]/(rho[ci[face]]*soundSpeed[ci[face]]*
            soundSpeed[ci[face]]) ;
        }
      }

      // Distribute contribution for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PressureCorrectionDiagonalCompressibleBoundary>
    registerPressureCorrectionDiagonalCompressibleBoundary ;

  // Rule to add the unsteady component to the diagonal term for
  // compressible flow. Assembling over cells. NOTE: Leaving the
  // time-step factor out for now as it seems to cause problems.
  class PressureCorrectionDiagonalUnsteadyCompressible : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_param<real> dt ;
      const_store<real> soundSpeed ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> D ;
    public:

      // Define input and output.
      PressureCorrectionDiagonalUnsteadyCompressible() {
        name_store("dt{n}",dt) ;
        name_store("soundSpeed{n,it}",soundSpeed) ;
        name_store("vol{n,it}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("pPrime_D{n,it}",D) ;
        input("dt{n},soundSpeed{n,it}") ;
        input("vol{n,it},cellRadius{n,it}") ;
        output("pPrime_D{n,it}") ;
        constraint("geom_cells,compressibleFlow") ;
      }

      // Distribute contribution for a single face.
      void calculate(Entity cell) {
        D[cell]+=vol[cell]/((*dt)*soundSpeed[cell]*soundSpeed[cell])*
          cellRadius[cell] ;
      }

      // Distribute contribution for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PressureCorrectionDiagonalUnsteadyCompressible>
    registerPressureCorrectionDiagonalUnsteadyCompressible ;

  // Rule to add the contribution to the diagonal near total-pressure inlets.
  // Need to UPDATE for real fluids.
  class PPrimeCoefficientTotalPressureInlet: public pointwise_rule {
    private:
      const_Map ci ;
      const_Map ref ;
      const_store<real> p0_BC ;
      const_store<EOS::State> eosState ;
      const_store<vect3d> v ;
      const_store<vect3d> v_f ;
      const_store<real> temperature_f ;
      const_store<Area> area ;
      store<real> c ;
    public:

      // Define input and output.
      PPrimeCoefficientTotalPressureInlet() {
        name_store("ci",ci) ;
        name_store("ref",ref) ;
        name_store("p0_BC",p0_BC) ;
        name_store("eos_state",eosState) ;
        name_store("v_f",v_f) ;
        name_store("v",v) ;
        name_store("temperature_f",temperature_f) ;
        name_store("area",area) ;
        name_store("pPrimeCoefficientTotalPressureInlet",c) ;
        input("ref->p0_BC,v_f,temperature_f,area,ci->(v,eos_state)") ;
        output("pPrimeCoefficientTotalPressureInlet") ;
        constraint("totalPressureInlet_BC") ;
      }

      // Distribute contribution for a single face.
      void calculate(Entity face) {

        // Compute the velocity component normal to the inlet.
        const real vNormal=-dot(v_f[face],area[face].n) ;
        const real vSquared=dot(v_f[face],v_f[face]) ;

        // Compute the coefficient that approximates the change in normal
        // velocity with respect to pressure.
        real gamma=eosState[ci[face]].Gamma(),gasConstant=eosState[ci[face]].
          gasConstant() ;
        c[face]=gasConstant*temperature_f[face]/(p0_BC[ref[face]]*
          vNormal*pow(1.0+(gamma-1.0)*vSquared/(2.0*gamma*gasConstant*
          temperature_f[face]),(1.0-2.0*gamma)/(gamma-1.0))) ;
      }

      // Distribute contribution for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeCoefficientTotalPressureInlet>
    registerPPrimeCoefficientTotalPressureInlet ;

  // Rule to add the contribution to the diagonal near total-pressure inlets.
  class PressureCorrectionDiagonalTotalPressureInlet: public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> rho_f ;
      const_store<Area> area ;
      const_store<real> c ;
      store<real> D ;
    public:

      // Define input and output.
      PressureCorrectionDiagonalTotalPressureInlet() {
        name_store("ci",ci) ;
        name_store("rho_f",rho_f) ;
        name_store("area",area) ;
        name_store("pPrimeCoefficientTotalPressureInlet",c) ;
        name_store("pPrime_D",D) ;
        input("rho_f,area,pPrimeCoefficientTotalPressureInlet") ;
        output("ci->pPrime_D") ;
        constraint("totalPressureInlet_BC") ;
      }

      // Distribute contribution for a single face.
      void calculate(Entity face) {
        D[ci[face]]+=rho_f[face]*c[face]*area[face].sada ;
      }

      // Distribute contribution for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PressureCorrectionDiagonalTotalPressureInlet>
    registerPressureCorrectionDiagonalTotalPressureInlet ;

  // Rule to compute the lower and upper terms for the linear system for
  // incompressible flow. Checked.
  class ComputePressureCorrectionLowerUpperIncompressible : public
  pointwise_rule {
    private:
      const_store<real> pPrimeCoefficient ;
      store<real> L,U ;
    public:

      // Define input and output.
      ComputePressureCorrectionLowerUpperIncompressible() {
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        name_store("pPrime_L",L) ;
        name_store("pPrime_U",U) ;
        input("pPrimeCoefficient") ;
        output("pPrime_L,pPrime_U") ;
        constraint("internalFaces") ;
      }

      // Initialize for a single face.
      void calculate(Entity face) {
        L[face]=-pPrimeCoefficient[face] ; U[face]=-pPrimeCoefficient[face] ;
      }

      // Initialize for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputePressureCorrectionLowerUpperIncompressible>
    registerComputePressureCorrectionLowerUpperIncompressible ;

  // Rule to compute the lower and upper terms for the linear system for
  // compressible flow. Replaced "massFlux" with "stageOneMassFlux" on
  // 01/07/2005.
  class ComputePressureCorrectionLowerUpperCompressible : public
  pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<real> rho ;
      const_store<real> soundSpeed ;
      const_store<real> pPrimeCoefficient ;
      const_store<real> massFlux ;
      store<real> L,U ;
    public:

      // Define input and output.
      ComputePressureCorrectionLowerUpperCompressible() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("rho",rho) ;
        name_store("soundSpeed",soundSpeed) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        name_store("stageOneMassFlux",massFlux) ;
        name_store("compressible::pPrime_L",L) ;
        name_store("compressible::pPrime_U",U) ;
        input("(cl,cr)->rho,(cl,cr)->soundSpeed") ;
        input("pPrimeCoefficient,stageOneMassFlux") ;
        output("compressible::pPrime_L,compressible::pPrime_U") ;
        constraint("internalFaces,compressibleFlow,compressiblePPrime") ;
      }

      // Initialize for a single face.
      void calculate(Entity face) {
        if(massFlux[face]>0.0){
          real compressibleCoefficient=massFlux[face]/(rho[cl[face]]*
            soundSpeed[cl[face]]*soundSpeed[cl[face]]) ;
          L[face]=-pPrimeCoefficient[face]-compressibleCoefficient ;
          U[face]=-pPrimeCoefficient[face] ;
        }else{
          real compressibleCoefficient=-massFlux[face]/(rho[cr[face]]*
            soundSpeed[cr[face]]*soundSpeed[cr[face]]) ;
          L[face]=-pPrimeCoefficient[face] ;
          U[face]=-pPrimeCoefficient[face]-compressibleCoefficient ;
        }
      }

      // Initialize for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputePressureCorrectionLowerUpperCompressible>
    registerComputePressureCorrectionLowerUpperCompressible ;

  // Rule to initialize the rhs term for the linear system. Checked.
  class InitializePressureCorrectionRHS : public unit_rule {
    private:
      store<real> B ;
    public:

      // Define input and output.
      InitializePressureCorrectionRHS() {
        name_store("pPrime_B",B) ;
        output("pPrime_B") ;
        constraint("vol") ;
      }

      // Initialize for a single cell.
      void calculate(Entity cell) { B[cell]=0.0 ; }

      // Initialize for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializePressureCorrectionRHS>
    registerInitializePressureCorrectionRHS ;

  // Rule to assemble the rhs term for the linear system. Right now
  // this rule is only good for incompressible flows with moving grids.
  class PressureCorrectionRHSInternal : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_param<real> thetaParameter ;
      const_store<real> stageOneMassFlux ;
      store<real> B ;
    public:

      // Define input and output.
      PressureCorrectionRHSInternal() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("stageOneMassFlux",stageOneMassFlux) ;
        name_store("pPrime_B",B) ;
        input("thetaParameter,stageOneMassFlux") ;
        output("(cl,cr)->pPrime_B") ;
        constraint("internalFaces") ;
      }

      // Distribute contribution for a single face.
      void calculate(Entity face) {
        B[cl[face]]-=stageOneMassFlux[face]*(*thetaParameter) ;
        B[cr[face]]+=stageOneMassFlux[face]*(*thetaParameter) ;
      }

      // Distribute contribution for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PressureCorrectionRHSInternal>
    registerPressureCorrectionRHSInternal ;

  // Rule to assemble the rhs term for the linear system. Right now
  // this rule is only good for incompressible flows with moving grids.
  class PressureCorrectionRHSBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_param<real> thetaParameter ;
      const_store<real> stageOneMassFlux ;
      store<real> B ;
    public:

      // Define input and output.
      PressureCorrectionRHSBoundary() {
        name_store("ci",ci) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("stageOneMassFlux",stageOneMassFlux) ;
        name_store("pPrime_B",B) ;
        input("thetaParameter,stageOneMassFlux") ;
        output("ci->pPrime_B") ;
        constraint("boundaryFaces") ;
      }

      // Distribute contribution for a single face.
      void calculate(Entity face) {
        B[ci[face]]-=stageOneMassFlux[face]*(*thetaParameter) ;
      }

      // Distribute contribution for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PressureCorrectionRHSBoundary>
    registerPressureCorrectionRHSBoundary ;

  // Rule to add the unsteady component to the right-hand side term.
  class PressureCorrectionRHSUnsteadyBDF : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> dt ;
      const_store<real> rhoOld,rhoNew ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> B ;
    public:

      // Define input and output.
      PressureCorrectionRHSUnsteadyBDF() {
        name_store("dt{n}",dt) ;
        name_store("rho{n}",rhoOld) ;
        name_store("rho{n,it}",rhoNew) ;
        name_store("vol{n,it}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("pPrime_B{n,it}",B) ;
        input("dt{n},rho{n},rho{n,it},vol{n,it},cellRadius{n,it}") ;
        output("pPrime_B{n,it}") ;
        constraint("compressibleFlow,BDFIntegrator,vol") ;
      }

      // Compute contribution for a single cell.
      void calculate(Entity cell) {
        B[cell]+=vol[cell]*(rhoOld[cell]-rhoNew[cell])/(*dt)*
          cellRadius[cell] ;
      }

      // Compute contribution for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PressureCorrectionRHSUnsteadyBDF>
    registerPressureCorrectionRHSUnsteadyBDF ;

  // Rule to add the unsteady component to the right-hand side term.
  class PressureCorrectionRHSUnsteadyBDF2 : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> dtOld,dt ;
      const_param<real> timeIntegratorFactor2 ;
      const_store<real> rhoOlder,rhoOld,rhoNew ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> B ;
    private:
      real f0,f1,f2 ;
    public:

      // Define input and output.
      PressureCorrectionRHSUnsteadyBDF2() {
        name_store("dt{n-1}",dtOld) ;
        name_store("dt{n}",dt) ;
        name_store("timeIntegratorFactor2{n}",timeIntegratorFactor2) ;
        name_store("rho{n-1}",rhoOlder) ;
        name_store("rho{n}",rhoOld) ;
        name_store("rho{n,it}",rhoNew) ;
        name_store("vol{n,it}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("pPrime_B{n,it}",B) ;
        input("dt{n-1},dt{n},timeIntegratorFactor2{n}") ;
        input("rho{n-1},rho{n},rho{n,it},vol{n,it},cellRadius{n,it}") ;
        output("pPrime_B{n,it}") ;
        constraint("compressibleFlow,BDF2Integrator,vol") ;
      }

      // Compute contribution for a single cell.
      void calculate(Entity cell) {
        B[cell]+=vol[cell]*cellRadius[cell]*((rhoOld[cell]-rhoNew[cell])/(*dt)-
          (*timeIntegratorFactor2)*(f2*rhoNew[cell]-f1*rhoOld[cell]+f0*
          rhoOlder[cell])) ;
      }

      // Compute contribution for a sequence of cells.
      virtual void compute(const sequence &seq) {
        real P=(*dt)+(*dtOld),Q=P*P/((*dt)*(*dt)) ;
        f0=1.0/((*dt)*Q-P),f1=Q*f0-(1.0/(*dt)),f2=f1-f0 ;
        do_loop(seq,this) ;
      }
  } ;

  register_rule<PressureCorrectionRHSUnsteadyBDF2>
    registerPressureCorrectionRHSUnsteadyBDF2 ;

//-----------------------------------------------------------------------------
// Rules for computing the residual of the pressure correction equation.

  // Create a variable for the residual of the pressure correction equation
  // for use with output rules.
  class PPrimeResidual : public pointwise_rule {
    private:
      const_param<real> rhoScale,vScale,lScale ;
      const_store<real> B ;
      store<real> pResidual ;
    private:
      real pPrimeFactor ;
    public:

      // Define input and output.
      PPrimeResidual() {
        name_store("rhoScale",rhoScale) ;
        name_store("vScale",vScale) ;
        name_store("lScale",lScale) ;
        name_store("pPrime_B",B) ;
        name_store("pResidual",pResidual) ;
        input("rhoScale,vScale,lScale,pPrime_B") ;
        output("pResidual") ;
        constraint("geom_cells") ;
      }

      // Assign cell residual value.
      void calculate(Entity cell) { pResidual[cell]=B[cell]/pPrimeFactor ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) {
        pPrimeFactor=(*rhoScale)*(*vScale)*(*lScale)*(*lScale) ;
        do_loop(seq,this) ;
      }
  } ;

  register_rule<PPrimeResidual> registerPPrimeResidual ;

  // Rule to initialize the total pressure correction residual. Checked.
  class InitializeTotalPressureCorrectionResidual : public unit_rule {
    private:
      param<ScalarResidual> pPrimeResidualData ;
    public:

      // Define input and output.
      InitializeTotalPressureCorrectionResidual() {
        name_store("pPrimeResidualData",pPrimeResidualData) ;
        output("pPrimeResidualData") ;
        constraint("geom_cells") ;
      }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) {
        *pPrimeResidualData=ScalarResidual() ;
      }
  } ;

  register_rule<InitializeTotalPressureCorrectionResidual>
    registerInitializeTotalPressureCorrectionResidual ;

  // Rule to compute the total pressure correction residual. Checked.
  class ComputeTotalPressureCorrectionResidual : public apply_rule
  <param<ScalarResidual>,ScalarResidualJoin> {
    private:
      const_store<real> B ;
      const_store<vect3d> cellCenter ;
      param<ScalarResidual> pPrimeResidualData ;
    public:

      // Define input and output.
      ComputeTotalPressureCorrectionResidual() {
        name_store("pPrime_B",B) ;
        name_store("cellcenter",cellCenter) ;
        name_store("pPrimeResidualData",pPrimeResidualData) ;
        input("pPrime_B,cellcenter") ;
        output("pPrimeResidualData") ;
        constraint("geom_cells") ;
      }

      // Add the cell contribution to the residual for a single cell.
      void calculate(Entity cell) {
        ScalarResidual temp ; temp.maxResidual=B[cell] ;
        temp.totalResidual=abs(B[cell]) ;
        temp.maxResidualLocation=cellCenter[cell] ;
        join(*pPrimeResidualData,temp) ;
      }

      // Add the cell contribution to the residual for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeTotalPressureCorrectionResidual>
    registerComputeTotalPressureCorrectionResidual ;

//-----------------------------------------------------------------------------
// Rules for correcting velocity, pressure and mass flux.

  // Default rule to set face pressure correction for boundary faces to zero.
  // Checked.
  class FacePressureCorrectionBoundaryDefault : public pointwise_rule {
    private:
      store<real> pPrimeFace ;
    public:

      // Define input and output.
      FacePressureCorrectionBoundaryDefault() {
        name_store("pPrimeFace",pPrimeFace) ;
        output("pPrimeFace") ;
        constraint("boundaryFaces") ;
      }

      // Set face pressure correction for a single face.
      void calculate(Entity face) { pPrimeFace[face]=0.0 ; }

      // Set face pressure correction for a sequence of faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FacePressureCorrectionBoundaryDefault>
    registerFacePressureCorrectionBoundaryDefault ;

  // Priority rule to set face pressure correction for boundary faces where
  // the pressure has not been specified. Checked.
  class FacePressureCorrectionBoundaryPriority : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> pPrime ;
      store<real> pPrimeFace ;
    public:

      // Define input and output.
      FacePressureCorrectionBoundaryPriority() {
        name_store("ci",ci) ;
        name_store("pPrime",pPrime) ;
        name_store("extrapolatedPressure_BC::pPrimeFace",pPrimeFace) ;
        input("ci->pPrime") ;
        output("extrapolatedPressure_BC::pPrimeFace") ;
        constraint("extrapolatedPressure_BC") ;
      }

      // Set face pressure correction for a single face.
      void calculate(Entity face) {
        pPrimeFace[face]=pPrime[ci[face]] ;
      }

      // Set face pressure correction for a sequence of faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FacePressureCorrectionBoundaryPriority>
    registerFacePressureCorrectionBoundaryPriority ;

  // Rule to initialize the velocity correction.
  class InitializeCorrectedVelocity : public unit_rule {
    private:
      store<vect3d> vCorrected ;
    public:

      // Define input and output. Added vMaxIterations constraint so this rule
      // does not get stranded for problems where we are not solving for the
      // velocity.
      InitializeCorrectedVelocity() {
        name_store("vCorrected",vCorrected) ;
        output("vCorrected") ;
        constraint("vol,vMaxIterations") ;
      }

      // Initialize velocity correction for a single cell.
      void calculate(Entity cell) { vCorrected[cell]=vect3d(0.0,0.0,0.0) ; }

      // Initialize velocity correction for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeCorrectedVelocity>
    registerInitializeCorrectedVelocity ;

  // Rule to initialize the velocity correction to vStar. This must be done as
  // an apply rule since the unit rule above must initialize to zero.
  class CorrectedVelocityFromVStar : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_store<vect3d> vStar ;
      store<vect3d> vCorrected ;
    public:

      // Define input and output.
      CorrectedVelocityFromVStar() {
        name_store("vStar",vStar) ;
        name_store("vCorrected",vCorrected) ;
        input("vStar") ;
        output("vCorrected") ;
        constraint("vol") ;
      }

      // Initialize velocity correction for a single cell.
      void calculate(Entity cell) { vCorrected[cell]+=vStar[cell] ; }

      // Initialize velocity correction for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CorrectedVelocityFromVStar>
    registerCorrectedVelocityFromVStar ;

  // Rule to add cell velocity corrections from interior faces. Checked.
  class CorrectVelocityInterior : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_param<real> twoDimensionFactor ;
      const_param<real> vRelaxationFactor ;
      const_Map cl,cr ;
      const_store<real> pPrime ;
      const_store<real> vMainCoefficient ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<vect3d> vCorrected ;
    public:

      // Define input and output.
      CorrectVelocityInterior() {
        name_store("twoDimensionFactor",twoDimensionFactor) ;
        name_store("vRelaxationFactor",vRelaxationFactor) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pPrime",pPrime) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("vCorrected",vCorrected) ;
        input("twoDimensionFactor,vRelaxationFactor,area,faceRadius") ;
        input("(cl,cr)->(pPrime,vMainCoefficient)") ;
        output("(cl,cr)->vCorrected") ;
        constraint("internalFaces") ;
      }

      // Correct the velocity for the cells attach to a single face.
      void calculate(Entity face) {
        vect3d temp=0.5*(pPrime[cl[face]]+pPrime[cr[face]])*
          (*vRelaxationFactor)*area[face].n*area[face].sada*faceRadius[face] ;
        temp.z*=(*twoDimensionFactor) ;
        vCorrected[cl[face]]-=temp/vMainCoefficient[cl[face]] ;
        vCorrected[cr[face]]+=temp/vMainCoefficient[cr[face]] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CorrectVelocityInterior>
    registerCorrectVelocityInterior ;

  // Rule to add cell velocity corrections from boundary faces. Checked.
  class CorrectVelocityBoundary : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_param<real> twoDimensionFactor ;
      const_param<real> vRelaxationFactor ;
      const_Map ci ;
      const_store<real> vMainCoefficient ;
      const_store<real> pPrimeFace ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<vect3d> vCorrected ;
    public:

      // Define input and output.
      CorrectVelocityBoundary() {
        name_store("twoDimensionFactor",twoDimensionFactor) ;
        name_store("vRelaxationFactor",vRelaxationFactor) ;
        name_store("ci",ci) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        name_store("pPrimeFace",pPrimeFace) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("vCorrected",vCorrected) ;
        input("twoDimensionFactor,vRelaxationFactor,pPrimeFace") ;
        input("area,faceRadius,ci->vMainCoefficient") ;
        output("ci->vCorrected") ;
        constraint("ci->geom_cells") ;
      }

      // Correct the velocity for the cells attach to a single face.
      void calculate(Entity face) {
        vect3d temp=(pPrimeFace[face]*(*vRelaxationFactor)/
          vMainCoefficient[ci[face]])*area[face].n*area[face].sada*
          faceRadius[face] ;
        temp.z*=(*twoDimensionFactor) ;
        vCorrected[ci[face]]-=temp ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CorrectVelocityBoundary> registerCorrectVelocityBoundary ;

  // Rule for correcting the cell pressure.
  class CorrectPressure : public pointwise_rule {
    private:
      const_param<real> pRelaxationFactor ;
      const_store<real> p ;
      const_store<real> pPrime ;
      store<real> pCorrected ;
    public:

      // Define input and output.
      CorrectPressure() {
        name_store("pRelaxationFactor",pRelaxationFactor) ;
        name_store("p",p) ;
        name_store("pPrime",pPrime) ;
        name_store("pCorrected",pCorrected) ;
        input("pRelaxationFactor,p,pPrime") ;
        output("pCorrected") ;
        constraint("geom_cells") ;
      }

      // Correct the pressure for a single cell.
      void calculate(Entity cell) {
        pCorrected[cell]=p[cell]+(*pRelaxationFactor)*pPrime[cell] ;
      }

      // Correct the pressure for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CorrectPressure> registerCorrectPressure ;

  // Rule for correcting the cell pressure next to total pressure inlets.
  class CorrectPressureTotalPressureInlet : public pointwise_rule {
    private:
      const_Map ci ;
      const_Map ref ;
      const_param<real> pRelaxationFactor ;
      const_store<real> p0_BC ;
      const_store<real> p ;
      const_store<real> pPrime ;
      store<real> pCorrected ;
    public:

      // Define input and output.
      CorrectPressureTotalPressureInlet() {
        name_store("ci",ci) ;
        name_store("ref",ref) ;
        name_store("pRelaxationFactor",pRelaxationFactor) ;
        name_store("p0_BC",p0_BC) ;
        name_store("p",p) ;
        name_store("pPrime",pPrime) ;
        name_store("totalPressureInlet::pCorrected",pCorrected) ;
        input("pRelaxationFactor,ci->p,ci->pPrime,ref->p0_BC") ;
        output("ci->(totalPressureInlet::pCorrected)") ;
        constraint("totalPressureInlet_BC") ;
      }

      // Correct the pressure for a single cell.
      void calculate(Entity face) {
        real newPressure=p[ci[face]]+(*pRelaxationFactor)*pPrime[ci[face]] ;
        pCorrected[ci[face]]=(newPressure<p0_BC[ref[face]])? newPressure:
          0.99*p0_BC[ref[face]] ;
      }

      // Correct the pressure for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CorrectPressureTotalPressureInlet>
    registerCorrectPressureTotalPressureInlet ;

  // Rule for corrected mass flux on interior faces. Checked.
  class CorrectMassFluxInterior : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<real> pPrime ;
      const_store<real> pPrimeCoefficient ;
      const_store<real> stageOneMassFlux ;
      store<real> massFluxCorrected ;
    public:

      // Define input and output.
      CorrectMassFluxInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pPrime",pPrime) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        name_store("stageOneMassFlux",stageOneMassFlux) ;
        name_store("massFluxCorrected",massFluxCorrected) ;
        input("(cl,cr)->pPrime,pPrimeCoefficient,stageOneMassFlux") ;
        output("massFluxCorrected") ;
        constraint("internalFaces") ;
      }

      // Compute corrected mass flux for a single face.
      void calculate(Entity face) {
        massFluxCorrected[face]=stageOneMassFlux[face]+pPrimeCoefficient[face]*
          (pPrime[cl[face]]-pPrime[cr[face]]) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CorrectMassFluxInterior> registerCorrectMassFluxInterior ;

  // Rule for corrected mass flux on boundary faces. Checked.
  class CorrectMassFluxBoundaryDefault : public pointwise_rule {
    private:
      const_store<real> stageOneMassFlux ;
      store<real> massFluxCorrected ;
    public:

      // Define input and output.
      CorrectMassFluxBoundaryDefault() {
        name_store("stageOneMassFlux",stageOneMassFlux) ;
        name_store("massFluxCorrected",massFluxCorrected) ;
        input("stageOneMassFlux") ;
        output("massFluxCorrected") ;
        constraint("boundaryFaces") ;
      }

      // Compute corrected mass flux for a single face.
      void calculate(Entity face) {
        massFluxCorrected[face]=stageOneMassFlux[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CorrectMassFluxBoundaryDefault>
    registerCorrectMassFluxBoundaryDefault ;

  // Rule for corrected mass flux on boundary faces where mass flux is not
  // specified and pressure is specified. Checked.
  class CorrectMassFluxBoundaryPriority : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> pPrime ;
      const_store<real> pPrimeCoefficient ;
      const_store<real> stageOneMassFlux ;
      store<real> massFluxCorrected ;
    public:

      // Define input and output.
      CorrectMassFluxBoundaryPriority() {
        name_store("ci",ci) ;
        name_store("pPrime",pPrime) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        name_store("stageOneMassFlux",stageOneMassFlux) ;
        name_store("priority::massFluxCorrected",massFluxCorrected) ;
        input("ci->pPrime,pPrimeCoefficient,stageOneMassFlux") ;
        output("priority::massFluxCorrected") ;
        constraint("nonSpecifiedMassFlux_BC,specifiedPressure_BC") ;
      }

      // Compute corrected mass flux for a single face.
      void calculate(Entity face) {
        massFluxCorrected[face]=stageOneMassFlux[face]+pPrimeCoefficient[face]*
          pPrime[ci[face]] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CorrectMassFluxBoundaryPriority>
    registerCorrectMassFluxBoundaryPriority ;

  // Rule for corrected mass flux on total pressure inlets.
  class CorrectMassFluxTotalPressureInlet : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> pPrime ;
      const_store<real> pPrimeCoefficient ;
      const_store<real> rho_f ;
      const_store<real> stageOneMassFlux ;
      const_store<Area> area ;
      store<real> massFluxCorrected ;
    public:

      // Define input and output.
      CorrectMassFluxTotalPressureInlet() {
        name_store("ci",ci) ;
        name_store("pPrime",pPrime) ;
        name_store("pPrimeCoefficientTotalPressureInlet",pPrimeCoefficient) ;
        name_store("rho_f",rho_f) ;
        name_store("stageOneMassFlux",stageOneMassFlux) ;
        name_store("area",area) ;
        name_store("totalPressureInlet::massFluxCorrected",massFluxCorrected) ;
        input("ci->pPrime,pPrimeCoefficientTotalPressureInlet") ;
        input("rho_f,area,stageOneMassFlux") ;
        output("totalPressureInlet::massFluxCorrected") ;
        constraint("totalPressureInlet_BC") ;
      }

      // Compute corrected mass flux for a single face.
      void calculate(Entity face) {
        massFluxCorrected[face]=stageOneMassFlux[face]+rho_f[face]*
          area[face].sada*pPrimeCoefficient[face]*pPrime[ci[face]] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

//register_rule<CorrectMassFluxTotalPressureInlet>
//  registerCorrectMassFluxTotalPressureInlet ;

}
