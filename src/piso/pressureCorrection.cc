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

// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;

// StreamUns includes.
#include "residual.h"
#include "sciTypes.h"
#include "varsFileInputs.h"

namespace streamUns {

  using Loci::Area ;

//-----------------------------------------------------------------------------
// Rules to process pressure correction equation options from the .vars file.

  // Creates the momentum equation solver constraints.
  class PressureEquationSolverConstraints : public constraint_rule {
    private:
      const_param<PressureEquationOptions> pressureEquationOptions ;
      Constraint pSGSLinearSolver,pPETSCLinearSolver,pHYPRELinearSolver ;
      Constraint pStarSGSLinearSolver,pStarPETSCLinearSolver ;
      Constraint pStarHYPRELinearSolver,compressiblePPrime ;
    public:
                                                                                
      // Define input and output.
      PressureEquationSolverConstraints() {
        name_store("pressureEquationOptions",pressureEquationOptions) ;
        name_store("pPrime_SGSLinearSolver",pSGSLinearSolver) ;
        name_store("pPrime_PETSCLinearSolver",pPETSCLinearSolver) ;
        name_store("pPrime_HYPRELinearSolver",pHYPRELinearSolver) ;
        name_store("pPrimeStar_SGSLinearSolver",pStarSGSLinearSolver) ;
        name_store("pPrimeStar_PETSCLinearSolver",pStarPETSCLinearSolver) ;
        name_store("pPrimeStar_HYPRELinearSolver",pStarHYPRELinearSolver) ;
        name_store("compressiblePPrime",compressiblePPrime) ;
        input("pressureEquationOptions") ;
        output("pPrime_SGSLinearSolver,pPrime_PETSCLinearSolver") ;
        output("pPrime_HYPRELinearSolver,compressiblePPrime") ;
        output("pPrimeStar_SGSLinearSolver,pPrimeStar_PETSCLinearSolver") ;
        output("pPrimeStar_HYPRELinearSolver") ;
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
                  pStarSGSLinearSolver=~EMPTY ; pStarPETSCLinearSolver=EMPTY ;
                  pStarHYPRELinearSolver=EMPTY ;
                }else if(name=="PETSC"){
                  pSGSLinearSolver=EMPTY ; pPETSCLinearSolver=~EMPTY ;
                  pHYPRELinearSolver=EMPTY ;
                  pStarSGSLinearSolver=EMPTY ; pStarPETSCLinearSolver=~EMPTY ;
                  pStarHYPRELinearSolver=EMPTY ;
                }else if(name=="HYPRE"){
                  pSGSLinearSolver=EMPTY ; pPETSCLinearSolver=EMPTY ;
                  pHYPRELinearSolver=~EMPTY ;
                  pStarSGSLinearSolver=EMPTY ; pStarPETSCLinearSolver=EMPTY ;
                  pStarHYPRELinearSolver=~EMPTY ;
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
          pStarSGSLinearSolver=~EMPTY ; pStarPETSCLinearSolver=EMPTY ;
          pStarHYPRELinearSolver=EMPTY ;
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

  // Creates the pressure equation solver parameters.
  class PressureEquationSolverParameters : public singleton_rule {
    private:
      const_param<PressureEquationOptions> pressureEquationOptions ;
      param<int> pMaxIterations,pStarMaxIterations,pPrime_numStages ;
    public:

      // Define input and output.
      PressureEquationSolverParameters() {
        name_store("pressureEquationOptions",pressureEquationOptions) ;
        name_store("pPrime_maxLinearSolverIterations",pMaxIterations) ;
        name_store("pPrimeStar_maxLinearSolverIterations",pStarMaxIterations) ;
        name_store("pPrime_numStages",pPrime_numStages) ;
        input("pressureEquationOptions") ;
        output("pPrime_maxLinearSolverIterations") ;
        output("pPrimeStar_maxLinearSolverIterations,pPrime_numStages") ;
      }

      // Set up the parameters.
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
                *pStarMaxIterations=int(temp) ;
              }
              break ;
            default:
              cerr << "Bad type for maxIterations in pressureEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *pMaxIterations=5 ;
          *pStarMaxIterations=5 ;
        }

        // Number of stages. More than one indicates that deferred-correction is
        // happening.
        if((*pressureEquationOptions).optionExists("numStages")){
          Loci::option_value_type optionValueType=pressureEquationOptions->
            getOptionValueType("numStages") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                pressureEquationOptions->getOption("numStages",temp) ;
                if(int(temp)<0){
                  cerr << "Bad numStages value for pressureEquation."
                    << endl ; Loci::Abort() ;
                }
                *pPrime_numStages=int(temp) ;
              }
              break ;
            default:
              cerr << "Bad type for numStages in pressureEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *pPrime_numStages=1 ;
        }
      }
  } ;
                                                                                
  register_rule<PressureEquationSolverParameters>
    registerPressureEquationSolverParamters ;

//-----------------------------------------------------------------------------
// Boundary conditions for pPrime and pPrimeStar.

  // Default rule for boundaries where pressure is specified.
  class BoundaryPPrimeSpecified : public pointwise_rule {
    private:
      store<real> pPrimeOld_f ;
    public:

      // Define input and output.
      BoundaryPPrimeSpecified() {
        name_store("pPrimeOld_f",pPrimeOld_f) ;
        output("pPrimeOld_f") ;
        constraint("boundaryFaces") ;
      }

      // Calculate for a single face.
      void calculate(Entity face) { pPrimeOld_f[face]=0.0 ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryPPrimeSpecified>
    registerBoundaryPPrimeSpecified ;

  // Priority rule for boundaries where pressure is extrapolated.
  class BoundaryPPrimeExtrapolation : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> pPrimeOld ;
      store<real> pPrimeOld_f ;
    public:

      // Define input and output.
      BoundaryPPrimeExtrapolation() {
        name_store("ci",ci) ;
        name_store("pPrimeOld",pPrimeOld) ;
        name_store("extrapolated::pPrimeOld_f",pPrimeOld_f) ;
        input("ci->pPrimeOld") ;
        output("extrapolated::pPrimeOld_f") ;
        constraint("extrapolatedPressure_BC") ;
      }

      // Calculate for a single face.
      void calculate(Entity face) {
        pPrimeOld_f[face]=pPrimeOld[ci[face]] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryPPrimeExtrapolation>
    registerBoundaryPPrimeExtrapolation ;

  // Default rule for boundaries where pressure is specified.
  class BoundaryPPrimeStarSpecified : public pointwise_rule {
    private:
      store<real> pPrimeStarOld_f ;
    public:

      // Define input and output.
      BoundaryPPrimeStarSpecified() {
        name_store("pPrimeStarOld_f",pPrimeStarOld_f) ;
        output("pPrimeStarOld_f") ;
        constraint("boundaryFaces") ;
      }

      // Calculate for a single face.
      void calculate(Entity face) { pPrimeStarOld_f[face]=0.0 ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryPPrimeStarSpecified>
    registerBoundaryPPrimeStarSpecified ;

  // Priority rule for boundaries where pressure is extrapolated.
  class BoundaryPPrimeStarExtrapolation : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> pPrimeStarOld ;
      store<real> pPrimeStarOld_f ;
    public:

      // Define input and output.
      BoundaryPPrimeStarExtrapolation() {
        name_store("ci",ci) ;
        name_store("pPrimeStarOld",pPrimeStarOld) ;
        name_store("extrapolated::pPrimeStarOld_f",pPrimeStarOld_f) ;
        input("ci->pPrimeStarOld") ;
        output("extrapolated::pPrimeStarOld_f") ;
        constraint("extrapolatedPressure_BC") ;
      }

      // Calculate for a single face.
      void calculate(Entity face) {
        pPrimeStarOld_f[face]=pPrimeStarOld[ci[face]] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryPPrimeStarExtrapolation>
    registerBoundaryPPrimeStarExtrapolation ;

//-----------------------------------------------------------------------------
// Rules common to the solution of pPrime in the predictor and pPrimeStar in
// the corrector.

  // Rule to compute the pressure-correction coefficient on interior faces for
  // the FOU and SOU convection schemes.
  class PPrimeCoefficientInteriorFOUSOU : public pointwise_rule {
    private:
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
        name_store("thetaParameter",thetaParameter) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        name_store("vol",vol) ;
        name_store("faceDensity",faceDensity) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("faceRadius",faceRadius) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        input("thetaParameter,faceRadius") ;
        input("(cl,cr)->(vMainCoefficient,vol),faceDensity,diffusionProduct") ;
        output("pPrimeCoefficient") ;
        constraint("internalFaces,fouOrSouInviscidFlux") ;
      }

      // Calculate the pressure correction coefficient for a single face.
      void calculate(Entity face) {
        pPrimeCoefficient[face]=0.5*faceDensity[face]*(vol[cl[face]]/
          vMainCoefficient[cl[face]]+vol[cr[face]]/vMainCoefficient[cr[face]])*
          diffusionProduct[face]*(*thetaParameter)*faceRadius[face]*
          faceRadius[face] ;
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
        name_store("thetaParameter",thetaParameter) ;
        name_store("ci",ci) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        name_store("vol",vol) ;
        name_store("rho_f",rho_f) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("faceRadius",faceRadius) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        input("thetaParameter,ci->(vMainCoefficient,vol)") ;
        input("rho_f,diffusionProduct,faceRadius") ;
        output("pPrimeCoefficient") ;
        constraint("boundaryFaces") ;
      }
                                                                                
      // Calculate the pressure-correction coefficient for a single face.
      void calculate(Entity face) {
        pPrimeCoefficient[face]=rho_f[face]*vol[ci[face]]/
          vMainCoefficient[ci[face]]*diffusionProduct[face]*(*thetaParameter)*
          faceRadius[face]*faceRadius[face] ;
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;
                                                                                
  register_rule<PPrimeCoefficientBoundary> registerPPrimeCoefficientBoundary ;

//-----------------------------------------------------------------------------
// Rules to establish the deferred-correction loops for pressure-correction in
// the predictor stage.

  // Build rule for pressure-correction iteration flag.
  class PPrimeIterationFinishedBuild : public singleton_rule {
    private:
      param<bool> pPrimeIterationFinished ;
    public:

      // Define input and output.
      PPrimeIterationFinishedBuild() {
        name_store("pPrimeIterationFinished{n,ip=-1}",pPrimeIterationFinished) ;
        output("pPrimeIterationFinished{n,ip=-1}") ;
        constraint("UNIVERSE{n}") ;
      }

      // Set to false.
      void compute(const sequence &seq) { *pPrimeIterationFinished=false ; }
  } ;

  register_rule<PPrimeIterationFinishedBuild>
    registerPPrimeIterationFinishedBuild ;

  // Class to determine if pressure-correction iteration is finished.
  class CheckPPrimeIterationFinished : public singleton_rule {
    private:
      const_param<int> ip ;
      const_param<int> pPrime_numStages ;
      param<bool> pPrimeIterationFinished ;
    public:

      // Define input and output.
      CheckPPrimeIterationFinished() {
        name_store("$ip{n,ip}",ip) ;
        name_store("pPrime_numStages{n,ip}",pPrime_numStages) ;
        name_store("pPrimeIterationFinished{n,ip}",pPrimeIterationFinished) ;
        input("$ip{n,ip},pPrime_numStages{n,ip}") ;
        output("pPrimeIterationFinished{n,ip}") ;
      }

      // Check if iteration is finished.
      void compute(const sequence &seq) {
        *pPrimeIterationFinished=(*ip==*pPrime_numStages-1) ;
      }
  } ;

  register_rule<CheckPPrimeIterationFinished> registerCheckPPrimeIterationFinished ;

  // Build rule for corrected density.
  class PredictorBuildDensityCorrected : public pointwise_rule {
    private:
      const_store<real> rho ;
      store<real> rhoCorrected ;
    public:

      // Define input and output.
      PredictorBuildDensityCorrected() {
        name_store("rho{n}",rho) ;
        name_store("rhoCorrected_p{n,ip=0}",rhoCorrected) ;
        input("rho{n}") ;
        output("rhoCorrected_p{n,ip=0}") ;
        constraint("geom_cells{n}") ;
      }

      // Set for a single cell.
      void calculate(Entity cell) { rhoCorrected[cell]=rho[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PredictorBuildDensityCorrected> registerPredictorBuildDensityCorrected ;

  // Advance rule for corrected density.
  class PredictorAdvanceDensityCorrected : public pointwise_rule {
    private:
      const_store<real> rhoPrime ;
      const_store<real> rhoCorrectedOld ;
      store<real> rhoCorrected ;
    public:

      // Define input and output.
      PredictorAdvanceDensityCorrected() {
        name_store("rhoPrime(pPrime){n,ip}",rhoPrime) ;
        name_store("rhoCorrected_p{n,ip}",rhoCorrectedOld) ;
        name_store("rhoCorrected_p{n,ip+1}",rhoCorrected) ;
        input("rhoPrime(pPrime){n,ip},rhoCorrected_p{n,ip}") ;
//      output("rhoCorrected_p{n,ip+1}=rhoCorrected_p{n,ip}") ;
        output("rhoCorrected_p{n,ip+1}") ;
        constraint("geom_cells{n,ip}") ;
      }

      // Set for a single cell.
//    void calculate(Entity cell) { rhoCorrected[cell]+=rhoPrime[cell] ; }
      void calculate(Entity cell) {
        rhoCorrected[cell]=rhoCorrectedOld[cell]+rhoPrime[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PredictorAdvanceDensityCorrected>
    registerPredictorAdvanceDensityCorrected ;

  // Build rule for corrected velocity.
  class PredictorBuildVCorrected : public pointwise_rule {
    private:
      const_store<vect3d> vStar ;
      store<vect3d> vCorrected ;
    public:

      // Define input and output.
      PredictorBuildVCorrected() {
        name_store("vStar{n}",vStar) ;
        name_store("vCorrected_p{n,ip=0}",vCorrected) ;
        input("vStar{n}") ;
        output("vCorrected_p{n,ip=0}") ;
        constraint("geom_cells{n}") ;
      }

      // Set for a single cell.
      void calculate(Entity cell) { vCorrected[cell]=vStar[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PredictorBuildVCorrected> registerPredictorBuildVCorrected ;

  // Advance rule for corrected velocity.
  class PredictorAdvanceVCorrected : public pointwise_rule {
    private:
      const_store<vect3d> vPrime ;
      const_store<vect3d> vCorrectedOld ;
      store<vect3d> vCorrected ;
    public:

      // Define input and output.
      PredictorAdvanceVCorrected() {
        name_store("vPrime(pPrime){n,ip}",vPrime) ;
        name_store("vCorrected_p{n,ip}",vCorrectedOld) ;
        name_store("vCorrected_p{n,ip+1}",vCorrected) ;
        input("vPrime(pPrime){n,ip},vCorrected_p{n,ip}") ;
//      output("vCorrected_p{n,ip+1}=vCorrected_p{n,ip}") ;
        output("vCorrected_p{n,ip+1}") ;
        constraint("geom_cells{n,ip}") ;
      }

      // Set for a single cell.
//    void calculate(Entity cell) { vCorrected[cell]+=vPrime[cell] ; }
      void calculate(Entity cell) {
        vCorrected[cell]=vCorrectedOld[cell]+vPrime[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PredictorAdvanceVCorrected> registerPredictorAdvanceVCorrected ;

  // Collapse rule for corrected velocity.
  class PredictorCollapseVCorrected : public pointwise_rule {
    private:
      store<vect3d> vCorrected ;
    public:

      // Define input and output.
      PredictorCollapseVCorrected() {
        name_store("vCorrected_p{n,ip}",vCorrected) ;
        input("vCorrected_p{n,ip}") ;
        output("vCorrected_p{n}=vCorrected_p{n,ip}") ;
        conditional("pPrimeIterationFinished{n,ip-1}") ;
        constraint("geom_cells{n,ip}") ;
      }

      // Empty compute method.
      void compute(const sequence &seq) {}
  } ;

  register_rule<PredictorCollapseVCorrected>
    registerPredictorCollapseVCorrected ;

  // Build rule for corrected pressure.
  class PredictorBuildPCorrected : public pointwise_rule {
    private:
      const_store<real> p ;
      store<real> pCorrected ;
    public:

      // Define input and output.
      PredictorBuildPCorrected() {
        name_store("p{n}",p) ;
        name_store("pCorrected_p{n,ip=0}",pCorrected) ;
        input("p{n}") ;
        output("pCorrected_p{n,ip=0}") ;
        constraint("geom_cells{n}") ;
      }

      // Set for a single cell.
      void calculate(Entity cell) { pCorrected[cell]=p[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PredictorBuildPCorrected> registerPredictorBuildPCorrected ;

  // Advance rule for corrected pressure.
  class PredictorAdvancePCorrected : public pointwise_rule {
    private:
      const_store<real> pPrime ;
      const_store<real> pCorrectedOld ;
      store<real> pCorrected ;
    public:

      // Define input and output.
      PredictorAdvancePCorrected() {
        name_store("pPrime{n,ip}",pPrime) ;
        name_store("pCorrected_p{n,ip}",pCorrectedOld) ;
        name_store("pCorrected_p{n,ip+1}",pCorrected) ;
        input("pPrime{n,ip},pCorrected_p{n,ip}") ;
//      output("pCorrected_p{n,ip+1}=pCorrected_p{n,ip}") ;
        output("pCorrected_p{n,ip+1}") ;
        constraint("geom_cells{n,ip}") ;
      }

      // Set for a single cell.
//    void calculate(Entity cell) { pCorrected[cell]+=pPrime[cell] ; }
      void calculate(Entity cell) {
        pCorrected[cell]=pCorrectedOld[cell]+pPrime[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PredictorAdvancePCorrected> registerPredictorAdvancePCorrected ;

  // Collapse rule for corrected pressure.
  class PredictorCollapsePCorrected : public pointwise_rule {
    private:
      store<real> pCorrected ;
    public:

      // Define input and output.
      PredictorCollapsePCorrected() {
        name_store("pCorrected_p{n,ip}",pCorrected) ;
        input("pCorrected_p{n,ip}") ;
        output("pCorrected_p{n}=pCorrected_p{n,ip}") ;
        conditional("pPrimeIterationFinished{n,ip-1}") ;
        constraint("geom_cells{n,ip}") ;
      }

      // Empty compute method.
      void compute(const sequence &seq) {}
  } ;

  register_rule<PredictorCollapsePCorrected>
    registerPredictorCollapsePCorrected ;

  // Build rule for corrected mass flux.
  class PredictorBuildMassFluxCorrected : public pointwise_rule {
    private:
      const_store<real> massFluxStar ;
      store<real> massFluxCorrected ;
    public:

      // Define input and output.
      PredictorBuildMassFluxCorrected() {
        name_store("massFluxStar{n}",massFluxStar) ;
        name_store("massFluxCorrected_p_temp{n,ip=0}",massFluxCorrected) ;
        input("massFluxStar{n}") ;
        output("massFluxCorrected_p_temp{n,ip=0}") ;
        constraint("massFluxCorrected{n}") ;
      }

      // Set for a single cell.
      void calculate(Entity cell) {
        massFluxCorrected[cell]=massFluxStar[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PredictorBuildMassFluxCorrected>
    registerPredictorBuildMassFluxCorrected ;

  // Advance rule for corrected mass flux.
  class PredictorAdvanceMassFluxCorrected : public pointwise_rule {
    private:
      const_store<real> massFluxPrime ;
      const_store<real> massFluxCorrectedOld ;
      store<real> massFluxCorrected ;
    public:

      // Define input and output.
      PredictorAdvanceMassFluxCorrected() {
        name_store("massFluxPrime(pPrime,pPrimeOld){n,ip}",massFluxPrime) ;
        name_store("massFluxCorrected_p_temp{n,ip}",massFluxCorrectedOld) ;
        name_store("massFluxCorrected_p_temp{n,ip+1}",massFluxCorrected) ;
        input("massFluxPrime(pPrime,pPrimeOld){n,ip},massFluxCorrected_p_temp{n,ip}") ;
//      output("massFluxCorrected_p{n,ip+1}=massFluxCorrected_p{n,ip}") ;
        output("massFluxCorrected_p_temp{n,ip+1}") ;
        constraint("massFluxCorrected{n,ip}") ;
      }

      // Set for a single cell.
//    void calculate(Entity cell) {
//      massFluxCorrected[cell]+=massFluxPrime[cell] ;
//    }
      void calculate(Entity cell) {
        massFluxCorrected[cell]=massFluxCorrectedOld[cell]+massFluxPrime[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PredictorAdvanceMassFluxCorrected>
    registerPredictorAdvanceMassFluxCorrected ;

  // Collapse rule for corrected mass flux.
  class PredictorCollapseMassFluxCorrected : public pointwise_rule {
    private:
      store<real> massFluxCorrected ;
    public:

      // Define input and output.
      PredictorCollapseMassFluxCorrected() {
        name_store("massFluxCorrected_p_temp{n,ip}",massFluxCorrected) ;
        input("massFluxCorrected_p_temp{n,ip}") ;
        output("massFluxCorrected_p_temp{n}=massFluxCorrected_p_temp{n,ip}") ;
        conditional("pPrimeIterationFinished{n,ip-1}") ;
        constraint("massFluxCorrected{n,ip}") ;
      }

      // Empty compute method.
      void compute(const sequence &seq) {}
  } ;

  register_rule<PredictorCollapseMassFluxCorrected>
    registerPredictorCollapseMassFluxCorrected ;

  // Build rule for lagged pressure correction.
  class PredictorBuildPPrimeOld : public pointwise_rule {
    private:
      store<real> pPrimeOld ;
    public:

      // Define input and output.
      PredictorBuildPPrimeOld() {
        name_store("pPrimeOld{n,ip=0}",pPrimeOld) ;
        output("pPrimeOld{n,ip=0}") ;
        constraint("geom_cells{n}") ;
      }

      // Set for a single cell.
      void calculate(Entity cell) { pPrimeOld[cell]=0.0 ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PredictorBuildPPrimeOld> registerPredictorBuildPPrimeOld ;

  // Advance rule for lagged pressure correction.
  class PredictorAdvancePPrimeOld : public pointwise_rule {
    private:
      const_store<real> pPrime ;
      store<real> pPrimeOld ;
    public:

      // Define input and output.
      PredictorAdvancePPrimeOld() {
        name_store("pPrime{n,ip}",pPrime) ;
        name_store("pPrimeOld{n,ip+1}",pPrimeOld) ;
        input("pPrime{n,ip}") ;
        output("pPrimeOld{n,ip+1}") ;
        constraint("geom_cells{n,ip}") ;
      }

      // Set for a single cell.
      void calculate(Entity cell) { pPrimeOld[cell]=pPrime[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PredictorAdvancePPrimeOld> registerPredictorAdvancePPrimeOld ;

//-----------------------------------------------------------------------------
// Rules to establish the deferred-correction loops for pressure-correction in
// the corrector stage.

  // Build rule for pressure-correction iteration flag.
  class PPrimeStarIterationFinishedBuild : public singleton_rule {
    private:
      param<bool> pPrimeStarIterationFinished ;
    public:

      // Define input and output.
      PPrimeStarIterationFinishedBuild() {
        name_store("pPrimeStarIterationFinished{n,it,ips=-1}",
          pPrimeStarIterationFinished) ;
        output("pPrimeStarIterationFinished{n,it,ips=-1}") ;
        constraint("UNIVERSE{n,it}") ;
      }

      // Set to false.
      void compute(const sequence &seq) { *pPrimeStarIterationFinished=false ; }
  } ;

  register_rule<PPrimeStarIterationFinishedBuild>
    registerPPrimeStarIterationFinishedBuild ;

  // Class to determine if pressure-correction iteration is finished.
  class CheckPPrimeStarIterationFinished : public singleton_rule {
    private:
      const_param<int> ips ;
      const_param<int> pPrime_numStages ;
      param<bool> pPrimeStarIterationFinished ;
    public:

      // Define input and output.
      CheckPPrimeStarIterationFinished() {
        name_store("$ips{n,it,ips}",ips) ;
        name_store("pPrime_numStages{n,it,ips}",pPrime_numStages) ;
        name_store("pPrimeStarIterationFinished{n,it,ips}",pPrimeStarIterationFinished) ;
        input("$ips{n,it,ips},pPrime_numStages{n,it,ips}") ;
        output("pPrimeStarIterationFinished{n,it,ips}") ;
      }

      // Check if iteration is finished.
      void compute(const sequence &seq) {
        *pPrimeStarIterationFinished=(*ips==*pPrime_numStages-1) ;
      }
  } ;

  register_rule<CheckPPrimeStarIterationFinished> registerCheckPPrimeStarIterationFinished ;

  // Build rule for the residual data.
  class CorrectorBuildPPrimeStarResidualData : public singleton_rule {
    private:
      param<ScalarResidual> pPrimeStarResidualData ;
    public:

      // Define input and output.
      CorrectorBuildPPrimeStarResidualData() {
        name_store("pPrimeStarResidualData{n,it,ips=0}",pPrimeStarResidualData) ;
        output("pPrimeStarResidualData{n,it,ips=0}") ;
        constraint("UNIVERSE{n,it}") ;
      }

      // Initialize the value.
      void compute(const sequence &seq) { *pPrimeStarResidualData=ScalarResidual() ; }
  } ;

  register_rule<CorrectorBuildPPrimeStarResidualData>
    registerCorrectorBuildPPrimeStarResidualData ;

  // Advance rule for the residual data.
  class CorrectorAdvancePPrimeStarResidualData : public singleton_rule {
    private:
      const_param<int> ips ;
      const_param<ScalarResidual> pPrimeStarResidualDataOld ;
      const_param<ScalarResidual> pPrimeStarResidualDataNew ;
      param<ScalarResidual> pPrimeStarResidualData ;
    public:

      // Define input and output.
      CorrectorAdvancePPrimeStarResidualData() {
        name_store("$ips{n,it,ips}",ips) ;
        name_store("pPrimeStarResidualData{n,it,ips}",pPrimeStarResidualDataOld) ;
        name_store("pPrimeStarResidualDataNew{n,it,ips}",pPrimeStarResidualDataNew) ;
        name_store("pPrimeStarResidualData{n,it,ips+1}",pPrimeStarResidualData) ;
        input("$ips{n,it,ips},pPrimeStarResidualData{n,it,ips}") ;
        input("pPrimeStarResidualDataNew{n,it,ips}") ;
        output("pPrimeStarResidualData{n,it,ips+1}") ;
      }

      // Only save the value for the first iteration.
      void compute(const sequence &seq) {
        if(*ips==0) *pPrimeStarResidualData=*pPrimeStarResidualDataNew ;
        else *pPrimeStarResidualData=*pPrimeStarResidualDataOld ;
      }
  } ;

  register_rule<CorrectorAdvancePPrimeStarResidualData>
    registerCorrectorAdvancePPrimeStarResidualData ;

  // Collapse rule for the residual data.
  class CorrectorCollapsePPrimeStarResidualData : public singleton_rule {
    private:
      const_param<ScalarResidual> pPrimeStarResidualDataOld ;
      param<ScalarResidual> pPrimeStarResidualData ;
    public:

      // Define input and output.
      CorrectorCollapsePPrimeStarResidualData() {
        name_store("pPrimeStarResidualData{n,it,ips}",pPrimeStarResidualDataOld) ;
        name_store("pPrimeStarResidualData{n,it}",pPrimeStarResidualData) ;
        input("pPrimeStarResidualData{n,it,ips}") ;
        output("pPrimeStarResidualData{n,it}") ;
        conditional("pPrimeStarIterationFinished{n,it,ips-1}") ;
      }

      // Do nothing.
      void compute(const sequence &seq) {
        *pPrimeStarResidualData=*pPrimeStarResidualDataOld ;
      }
  } ;

  register_rule<CorrectorCollapsePPrimeStarResidualData>
    registerCorrectorCollapsePPrimeStarResidualData ;

  // Build rule for corrected density.
  class CorrectorBuildDensityCorrected : public pointwise_rule {
    private:
      const_store<real> rho ;
      store<real> rhoCorrected ;
    public:

      // Define input and output.
      CorrectorBuildDensityCorrected() {
        name_store("rho{n,it}",rho) ;
        name_store("rhoCorrected_c{n,it,ips=0}",rhoCorrected) ;
        input("rho{n,it}") ;
        output("rhoCorrected_c{n,it,ips=0}") ;
        constraint("geom_cells{n,it}") ;
      }

      // Set for a single cell.
      void calculate(Entity cell) { rhoCorrected[cell]=rho[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CorrectorBuildDensityCorrected>
    registerCorrectorBuildDensityCorrected ;

  // Advance rule for corrected density.
  class CorrectorAdvanceDensityCorrected : public pointwise_rule {
    private:
      const_store<real> rhoPrime ;
      const_store<real> rhoCorrectedOld ;
      store<real> rhoCorrected ;
    public:

      // Define input and output.
      CorrectorAdvanceDensityCorrected() {
        name_store("rhoPrime(pPrimeStar){n,it,ips}",rhoPrime) ;
        name_store("rhoCorrected_c{n,it,ips}",rhoCorrectedOld) ;
        name_store("rhoCorrected_c{n,it,ips+1}",rhoCorrected) ;
        input("rhoPrime(pPrimeStar){n,it,ips},rhoCorrected_c{n,it,ips}") ;
//      output("rhoCorrected_c{n,it,ips+1}=rhoCorrected_c{n,it,ips}") ;
        output("rhoCorrected_c{n,it,ips+1}") ;
        constraint("geom_cells{n,it,ips}") ;
      }

      // Set for a single cell.
//    void calculate(Entity cell) { rhoCorrected[cell]+=rhoPrime[cell] ; }
      void calculate(Entity cell) {
        rhoCorrected[cell]=rhoCorrectedOld[cell]+rhoPrime[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CorrectorAdvanceDensityCorrected>
    registerCorrectorAdvanceDensityCorrected ;

  // Build rule for corrected velocity.
  class CorrectorBuildVCorrected : public pointwise_rule {
    private:
      const_store<vect3d> vStarHat ;
      store<vect3d> vCorrected ;
    public:

      // Define input and output.
      CorrectorBuildVCorrected() {
        name_store("vStarHat{n,it}",vStarHat) ;
        name_store("vCorrected_c{n,it,ips=0}",vCorrected) ;
        input("vStarHat{n,it}") ;
        output("vCorrected_c{n,it,ips=0}") ;
        constraint("geom_cells{n,it}") ;
      }

      // Set for a single cell.
      void calculate(Entity cell) { vCorrected[cell]=vStarHat[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CorrectorBuildVCorrected> registerCorrectorBuildVCorrected ;

  // Advance rule for corrected velocity.
  class CorrectorAdvanceVCorrected : public pointwise_rule {
    private:
      const_store<vect3d> vPrime ;
      const_store<vect3d> vCorrectedOld ;
      store<vect3d> vCorrected ;
    public:

      // Define input and output.
      CorrectorAdvanceVCorrected() {
        name_store("vPrime(pPrimeStar){n,it,ips}",vPrime) ;
        name_store("vCorrected_c{n,it,ips}",vCorrectedOld) ;
        name_store("vCorrected_c{n,it,ips+1}",vCorrected) ;
        input("vPrime(pPrimeStar){n,it,ips},vCorrected_c{n,it,ips}") ;
//      output("vCorrected_c{n,it,ips+1}=vCorrected_c{n,it,ips}") ;
        output("vCorrected_c{n,it,ips+1}") ;
        constraint("geom_cells{n,it,ips}") ;
      }

      // Set for a single cell.
//    void calculate(Entity cell) { vCorrected[cell]+=vPrime[cell] ; }
      void calculate(Entity cell) {
        vCorrected[cell]=vCorrectedOld[cell]+vPrime[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CorrectorAdvanceVCorrected> registerCorrectorAdvanceVCorrected ;

  // Collapse rule for corrected velocity.
  class CorrectorCollapseVCorrected : public pointwise_rule {
    private:
      store<vect3d> vCorrected ;
    public:

      // Define input and output.
      CorrectorCollapseVCorrected() {
        name_store("vCorrected_c{n,it,ips}",vCorrected) ;
        input("vCorrected_c{n,it,ips}") ;
        output("vCorrected_c{n,it}=vCorrected_c{n,it,ips}") ;
        conditional("pPrimeStarIterationFinished{n,it,ips-1}") ;
        constraint("geom_cells{n,it,ips}") ;
      }

      // Empty compute method.
      void compute(const sequence &seq) {}
  } ;

  register_rule<CorrectorCollapseVCorrected>
    registerCorrectorCollapseVCorrected ;

  // Build rule for corrected pressure.
  class CorrectorBuildPCorrected : public pointwise_rule {
    private:
      const_store<real> p ;
      store<real> pCorrected ;
    public:

      // Define input and output.
      CorrectorBuildPCorrected() {
        name_store("p{n,it}",p) ;
        name_store("pCorrected_c{n,it,ips=0}",pCorrected) ;
        input("p{n,it}") ;
        output("pCorrected_c{n,it,ips=0}") ;
        constraint("geom_cells{n,it}") ;
      }

      // Set for a single cell.
      void calculate(Entity cell) { pCorrected[cell]=p[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CorrectorBuildPCorrected> registerCorrectorBuildPCorrected ;

  // Advance rule for corrected pressure.
  class CorrectorAdvancePCorrected : public pointwise_rule {
    private:
      const_store<real> pPrimeStar ;
      const_store<real> pCorrectedOld ;
      store<real> pCorrected ;
    public:

      // Define input and output.
      CorrectorAdvancePCorrected() {
        name_store("pPrimeStar{n,it,ips}",pPrimeStar) ;
        name_store("pCorrected_c{n,it,ips}",pCorrectedOld) ;
        name_store("pCorrected_c{n,it,ips+1}",pCorrected) ;
        input("pPrimeStar{n,it,ips},pCorrected_c{n,it,ips}") ;
//      output("pCorrected_c{n,it,ips+1}=pCorrected_c{n,it,ips}") ;
        output("pCorrected_c{n,it,ips+1}") ;
        constraint("geom_cells{n,it,ips}") ;
      }

      // Set for a single cell.
//    void calculate(Entity cell) { pCorrected[cell]+=pPrimeStar[cell] ; }
      void calculate(Entity cell) {
        pCorrected[cell]=pCorrectedOld[cell]+pPrimeStar[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CorrectorAdvancePCorrected> registerCorrectorAdvancePCorrected ;

  // Collapse rule for corrected pressure.
  class CorrectorCollapsePCorrected : public pointwise_rule {
    private:
      store<real> pCorrected ;
    public:

      // Define input and output.
      CorrectorCollapsePCorrected() {
        name_store("pCorrected_c{n,it,ips}",pCorrected) ;
        input("pCorrected_c{n,it,ips}") ;
        output("pCorrected_c{n,it}=pCorrected_c{n,it,ips}") ;
        conditional("pPrimeStarIterationFinished{n,it,ips-1}") ;
        constraint("geom_cells{n,it,ips}") ;
      }

      // Empty compute method.
      void compute(const sequence &seq) {}
  } ;

  register_rule<CorrectorCollapsePCorrected>
    registerCorrectorCollapsePCorrected ;

  // Build rule for corrected mass flux.
  class CorrectorBuildMassFluxCorrected : public pointwise_rule {
    private:
      const_store<real> massFluxStarHat ;
      store<real> massFluxCorrected ;
    public:

      // Define input and output.
      CorrectorBuildMassFluxCorrected() {
        name_store("massFluxStarHat{n,it}",massFluxStarHat) ;
        name_store("massFluxCorrected_c_temp{n,it,ips=0}",massFluxCorrected) ;
        input("massFluxStarHat{n,it}") ;
        output("massFluxCorrected_c_temp{n,it,ips=0}") ;
        constraint("massFluxCorrected{n,it}") ;
      }

      // Set for a single cell.
      void calculate(Entity cell) {
        massFluxCorrected[cell]=massFluxStarHat[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CorrectorBuildMassFluxCorrected>
    registerCorrectorBuildMassFluxCorrected ;

  // Advance rule for corrected mass flux.
  class CorrectorAdvanceMassFluxCorrected : public pointwise_rule {
    private:
      const_store<real> massFluxPrime ;
      const_store<real> massFluxCorrectedOld ;
      store<real> massFluxCorrected ;
    public:

      // Define input and output.
      CorrectorAdvanceMassFluxCorrected() {
        name_store("massFluxPrime(pPrimeStar,pPrimeStarOld){n,it,ips}",massFluxPrime) ;
        name_store("massFluxCorrected_c_temp{n,it,ips}",massFluxCorrectedOld) ;
        name_store("massFluxCorrected_c_temp{n,it,ips+1}",massFluxCorrected) ;
        input("massFluxPrime(pPrimeStar,pPrimeStarOld){n,it,ips},massFluxCorrected_c_temp{n,it,ips}") ;
//      output("massFluxCorrected_c{n,it,ips+1}=massFluxCorrected_c{n,it,ips}") ;
        output("massFluxCorrected_c_temp{n,it,ips+1}") ;
        constraint("massFluxCorrected{n,it,ips}") ;
      }

      // Set for a single cell.
//    void calculate(Entity cell) {
//      massFluxCorrected[cell]+=massFluxPrime[cell] ;
//    }
      void calculate(Entity cell) {
        massFluxCorrected[cell]=massFluxCorrectedOld[cell]+massFluxPrime[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CorrectorAdvanceMassFluxCorrected>
    registerCorrectorAdvanceMassFluxCorrected ;

  // Collapse rule for corrected mass flux.
  class CorrectorCollapseMassFluxCorrected : public pointwise_rule {
    private:
      store<real> massFluxCorrected ;
    public:

      // Define input and output.
      CorrectorCollapseMassFluxCorrected() {
        name_store("massFluxCorrected_c_temp{n,it,ips}",massFluxCorrected) ;
        input("massFluxCorrected_c_temp{n,it,ips}") ;
        output("massFluxCorrected_c_temp{n,it}=massFluxCorrected_c_temp{n,it,ips}") ;
        conditional("pPrimeStarIterationFinished{n,it,ips-1}") ;
        constraint("massFluxCorrected{n,it,ips}") ;
      }

      // Empty compute method.
      void compute(const sequence &seq) {}
  } ;

  register_rule<CorrectorCollapseMassFluxCorrected>
    registerCorrectorCollapseMassFluxCorrected ;

  // Build rule for lagged pressure correction.
  class CorrectorBuildPPrimeStarOld : public pointwise_rule {
    private:
      store<real> pPrimeStarOld ;
    public:

      // Define input and output.
      CorrectorBuildPPrimeStarOld() {
        name_store("pPrimeStarOld{n,it,ips=0}",pPrimeStarOld) ;
        output("pPrimeStarOld{n,it,ips=0}") ;
        constraint("geom_cells{n,it}") ;
      }

      // Set for a single cell.
      void calculate(Entity cell) { pPrimeStarOld[cell]=0.0 ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CorrectorBuildPPrimeStarOld>
    registerCorrectorBuildPPrimeStarOld ;

  // Advance rule for lagged pressure correction.
  class CorrectorAdvancePPrimeStarOld : public pointwise_rule {
    private:
      const_store<real> pPrimeStar ;
      store<real> pPrimeStarOld ;
    public:

      // Define input and output.
      CorrectorAdvancePPrimeStarOld() {
        name_store("pPrimeStar{n,it,ips}",pPrimeStar) ;
        name_store("pPrimeStarOld{n,it,ips+1}",pPrimeStarOld) ;
        input("pPrimeStar{n,it,ips}") ;
        output("pPrimeStarOld{n,it,ips+1}") ;
        constraint("geom_cells{n,it,ips}") ;
      }

      // Set for a single cell.
      void calculate(Entity cell) { pPrimeStarOld[cell]=pPrimeStar[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CorrectorAdvancePPrimeStarOld>
    registerCorrectorAdvancePPrimeStarOld ;

//-----------------------------------------------------------------------------
// Rules for assembling coefficients for pPrime equation in the predictor.

  // Rule to initialize the diagonal term for the linear system.
  class InitializePPrimeDiagonal : public unit_rule {
    private:
      store<real> D ;
    public:

      // Define input and output.
      InitializePPrimeDiagonal() {
        name_store("pPrime_D",D) ;
        output("pPrime_D") ;
        constraint("vol") ;
      }

      // Initialize for a single cell.
      void calculate(Entity cell) { D[cell]=0.0 ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializePPrimeDiagonal>
    registerInitializePPrimeDiagonal ;

  // Rule to assemble the diagonal term for the linear system for
  // incompressible flow.
  class PPrimeDiagonalIncompressibleInternal : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> pPrimeCoefficient ;
      store<real> D ;
    public:

      // Define input and output.
      PPrimeDiagonalIncompressibleInternal() {
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

  register_rule<PPrimeDiagonalIncompressibleInternal>
    registerPPrimeDiagonalIncompressibleInternal ;

  // Rule to assemble the diagonal term for the linear system. Assembling over
  // boundary faces. Note that this rule now applies for both incompressible
  // and compressible flows.  NOTE: Removed the nonSpecifiedMassFlux_BC constraint
  // so that this rule will be used for interpolated interface boundaries as
  // well where the pressure is specified via interpolation.
  class PPrimeDiagonalBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> pPrimeCoefficient ;
      store<real> D ;
    public:

      // Define input and output.
      PPrimeDiagonalBoundary() {
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

  register_rule<PPrimeDiagonalBoundary> registerPPrimeDiagonalBoundary ;

  // Rule to assemble the diagonal term for the linear system for
  // compressible flow.
  class PPrimeDiagonalCompressibleInternal : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> rho ;
      const_store<real> soundSpeed ;
      const_store<real> massFluxStar ;
      store<real> D ;
    public:

      // Define input and output.
      PPrimeDiagonalCompressibleInternal() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("rho",rho) ;
        name_store("soundSpeed",soundSpeed) ;
        name_store("pPrime_D",D) ;
        name_store("massFluxStar",massFluxStar) ;
        input("(cl,cr)->(rho,soundSpeed)") ;
        input("massFluxStar") ;
        output("(cl,cr)->pPrime_D") ;
        constraint("internalFaces,compressibleFlow,compressiblePPrime") ;
      }

      // Distribute contribution for a single face.
      void calculate(Entity face) {
        if(massFluxStar[face]>0.0){
          real compressibleCoefficient=massFluxStar[face]/(rho[cl[face]]*
            soundSpeed[cl[face]]*soundSpeed[cl[face]]) ;
          D[cl[face]]+=compressibleCoefficient ;
        }else{
          real compressibleCoefficient=-massFluxStar[face]/(rho[cr[face]]*
            soundSpeed[cr[face]]*soundSpeed[cr[face]]) ;
          D[cr[face]]+=compressibleCoefficient ;
        }
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeDiagonalCompressibleInternal>
    registerPPrimeDiagonalCompressibleInternal ;

  // Rule to assemble the diagonal term for the linear system for
  // compressible flow.
  class PPrimeDiagonalCompressibleBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> rho ;
      const_store<real> soundSpeed ;
      const_store<real> massFlux ;
      store<real> D ;
    public:

      // Define input and output.
      PPrimeDiagonalCompressibleBoundary() {
        name_store("ci",ci) ;
        name_store("rho",rho) ;
        name_store("soundSpeed",soundSpeed) ;
        name_store("massFlux",massFlux) ;
        name_store("pPrime_D",D) ;
        input("ci->(rho,soundSpeed),massFlux") ;
        output("ci->pPrime_D") ;
        constraint("nonSpecifiedMassFlux_BC,compressibleFlow") ;
        constraint("compressiblePPrime") ;
      }

      // Distribute contribution for a single face.
      void calculate(Entity face) {
        if(massFlux[face]>0.0) D[ci[face]]+=massFlux[face]/(rho[ci[face]]*
          soundSpeed[ci[face]]*soundSpeed[ci[face]]) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeDiagonalCompressibleBoundary>
    registerPPrimeDiagonalCompressibleBoundary ;

  // Rule to add the unsteady component to the diagonal term for
  // compressible flow. NOTE: Leaving the time-step factor out for
  // now as it seems to cause problems.
  class PPrimeDiagonalUnsteadyCompressible : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> timeStep ;
      const_param<real> timeIntegratorFactor0 ;
      const_store<real> soundSpeed ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> D ;
    public:

      // Define input and output.
      PPrimeDiagonalUnsteadyCompressible() {
        name_store("timeStep",timeStep) ;
        name_store("timeIntegratorFactor0",timeIntegratorFactor0) ;
        name_store("soundSpeed",soundSpeed) ;
        name_store("vol",vol) ;
        name_store("cellRadius",cellRadius) ;
        name_store("pPrime_D",D) ;
        input("timeStep,timeIntegratorFactor0,soundSpeed,vol,cellRadius") ;
        output("pPrime_D") ;
        constraint("geom_cells,compressibleFlow") ;
      }

      // Distribute contribution for a single face.
      void calculate(Entity cell) {
        D[cell]+=0.5*((*timeIntegratorFactor0)+1.0)*vol[cell]/((*timeStep)*
          soundSpeed[cell]*soundSpeed[cell])*cellRadius[cell] ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeDiagonalUnsteadyCompressible>
    registerPPrimeDiagonalUnsteadyCompressible ;

  // Rule to compute the lower and upper terms for the linear system for
  // incompressible flow.
  class PPrimeLowerUpperIncompressible : public pointwise_rule {
    private:
      const_store<real> pPrimeCoefficient ;
      store<real> L,U ;
    public:

      // Define input and output.
      PPrimeLowerUpperIncompressible() {
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

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeLowerUpperIncompressible>
    registerPPrimeLowerUpperIncompressible ;

  // Rule to compute the lower and upper terms for the linear system for
  // compressible flow.
  class PPrimeLowerUpperCompressible : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<real> rho ;
      const_store<real> soundSpeed ;
      const_store<real> pPrimeCoefficient ;
      const_store<real> massFluxStar ;
      store<real> L,U ;
    public:

      // Define input and output.
      PPrimeLowerUpperCompressible() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("rho",rho) ;
        name_store("soundSpeed",soundSpeed) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
//      name_store("massFluxStar",massFluxStar) ;
        name_store("massFluxCorrected_p_temp",massFluxStar) ;
        name_store("compressible::pPrime_L",L) ;
        name_store("compressible::pPrime_U",U) ;
        input("(cl,cr)->(rho,soundSpeed)") ;
//      input("pPrimeCoefficient,massFluxStar") ;
        input("pPrimeCoefficient,massFluxCorrected_p_temp") ;
        output("compressible::pPrime_L,compressible::pPrime_U") ;
        constraint("internalFaces,compressibleFlow,compressiblePPrime") ;
      }

      // Compute for a single face.
      void calculate(Entity face) {
        if(massFluxStar[face]>0.0){
          real compressibleCoefficient=massFluxStar[face]/(rho[cl[face]]*
            soundSpeed[cl[face]]*soundSpeed[cl[face]]) ;
          L[face]=-pPrimeCoefficient[face]-compressibleCoefficient ;
          U[face]=-pPrimeCoefficient[face] ;
        }else{
          real compressibleCoefficient=-massFluxStar[face]/(rho[cr[face]]*
            soundSpeed[cr[face]]*soundSpeed[cr[face]]) ;
          L[face]=-pPrimeCoefficient[face] ;
          U[face]=-pPrimeCoefficient[face]-compressibleCoefficient ;
        }
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeLowerUpperCompressible>
    registerPPrimeLowerUpperCompressible ;

  // Rule to initialize the rhs term for the linear system.
  class InitializePPrimeRHS : public unit_rule {
    private:
      store<real> B ;
    public:

      // Define input and output.
      InitializePPrimeRHS() {
        name_store("pPrime_B",B) ;
        output("pPrime_B") ;
        constraint("vol") ;
      }

      // Initialize for a single cell.
      void calculate(Entity cell) { B[cell]=0.0 ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializePPrimeRHS> registerInitializePPrimeRHS ;

  // Rule to assemble the rhs term for the linear system.
  class PPrimeRHSInternal : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<vect3d> cellCenter,pPrimeGradient,faceCenter ;
      const_store<Area> area ;
      const_store<real> pPrimeCoefficient,massFluxStar ;
      store<real> B ;
    public:

      // Define input and output.
      PPrimeRHSInternal() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("cellcenter",cellCenter) ;
        name_store("grads(pPrimeOld)",pPrimeGradient) ;
        name_store("facecenter",faceCenter) ;
        name_store("area",area) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        name_store("massFluxCorrected_p_temp",massFluxStar) ;
        name_store("pPrime_B",B) ;
        input("(cl,cr)->(cellcenter,grads(pPrimeOld))") ;
        input("facecenter,area,pPrimeCoefficient,massFluxCorrected_p_temp") ;
        output("(cl,cr)->pPrime_B") ;
        constraint("internalFaces") ;
      }

      // Distribute contribution for a single face.
      void calculate(Entity face) {
        vect3d rL=faceCenter[face]-dot(faceCenter[face]-cellCenter[cl[face]],
          area[face].n)*area[face].n ;
        vect3d rR=faceCenter[face]-dot(faceCenter[face]-cellCenter[cr[face]],
          area[face].n)*area[face].n ;
        real laggedSourceTerm=pPrimeCoefficient[face]*(dot(pPrimeGradient[cl[face]],
          rL-cellCenter[cl[face]])-dot(pPrimeGradient[cr[face]],rR-cellCenter[cr[face]])) ;
        B[cl[face]]-=massFluxStar[face]+laggedSourceTerm ;
        B[cr[face]]+=massFluxStar[face]+laggedSourceTerm ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeRHSInternal> registerPPrimeRHSInternal ;

  // Rule to assemble the rhs term for the linear system.
  class PPrimeRHSBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_param<real> thetaParameter ;
      const_store<real> massFlux ;
      store<real> B ;
    public:

      // Define input and output.
      PPrimeRHSBoundary() {
        name_store("ci",ci) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("massFlux",massFlux) ;
        name_store("pPrime_B",B) ;
        input("thetaParameter,massFlux") ;
        output("ci->pPrime_B") ;
        constraint("boundaryFaces") ;
      }

      // Compute for a single face.
      void calculate(Entity face) {
        B[ci[face]]-=massFlux[face]*(*thetaParameter) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeRHSBoundary> registerPPrimeRHSBoundary ;

  // Rule to add the unsteady component to the right-hand side term.
  // This rule will have to be updated for deforming meshes.
  class PPrimeRHSUnsteadyBDF : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> timeStep ;
      const_store<real> rho,rhoOld ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> B ;
    public:

      // Define input and output. Note the {n,it,ips} on the compressible
      // flow constraint. Without this, Loci complains that this should be
      // a build rule.
      PPrimeRHSUnsteadyBDF() {
        name_store("timeStep{n}",timeStep) ;
        name_store("rho{n}",rhoOld) ;
        name_store("rhoCorrected_p{n,ip}",rho) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n}",cellRadius) ;
        name_store("pPrime_B{n,ip}",B) ;
        input("timeStep{n},rho{n},rhoCorrected_p{n,ip},vol{n},cellRadius{n}") ;
        output("pPrime_B{n,ip}") ;
        constraint("compressibleFlow{n,ip},BDFIntegrator,vol") ;
      }

      // Compute contribution for a single cell.
      void calculate(Entity cell) {
        B[cell]+=vol[cell]*cellRadius[cell]*(rhoOld[cell]-rho[cell])/
          (*timeStep) ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeRHSUnsteadyBDF> registerPPrimeRHSUnsteadyBDF ;

  // Rule to add the unsteady component to the right-hand side term.
  // This rule will have to be updated for deforming meshes.
  class PPrimeRHSUnsteady : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<int> n ;
      const_param<real> timeStep ;
      const_store<real> rho,rhoOld ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> B ;
    public:

      // Define input and output.
      PPrimeRHSUnsteady() {
        name_store("$n{n}",n) ;
        name_store("timeStep{n}",timeStep) ;
        name_store("rho{n-1}",rhoOld) ;
        name_store("rho{n}",rho) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n}",cellRadius) ;
        name_store("pPrime_B{n}",B) ;
        input("$n{n},timeStep{n},rho{n-1},rho{n},vol{n},cellRadius{n}") ;
        output("pPrime_B{n}") ;
        constraint("compressibleFlow,BDF2Integrator,vol") ;
      }

      // Compute contribution for a single cell. Use BDF on first timestep,
      // which amounts to adding nothing for n=0.
      void calculate(Entity cell) {
        if((*n)!=0.0) B[cell]+=0.5*vol[cell]*cellRadius[cell]*(rho[cell]-
          rhoOld[cell])/(*timeStep) ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeRHSUnsteady> registerPPrimeRHSUnsteady ;

//-----------------------------------------------------------------------------
// Rules for assembling coefficients for pPrimeStar equation in the corrector.

  // Rule to initialize the diagonal term for the linear system.
  class InitializePPrimeStarDiagonal : public unit_rule {
    private:
      store<real> D ;
    public:

      // Define input and output.
      InitializePPrimeStarDiagonal() {
        name_store("pPrimeStar_D",D) ;
        output("pPrimeStar_D") ;
        constraint("vol") ;
      }

      // Initialize for a single cell.
      void calculate(Entity cell) { D[cell]=0.0 ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializePPrimeStarDiagonal>
    registerInitializePPrimeStarDiagonal ;

  // Rule to assemble the diagonal term for the linear system for
  // incompressible flow.
  class PPrimeStarDiagonalIncompressibleInternal : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> pPrimeCoefficient ;
      store<real> D ;
    public:

      // Define input and output.
      PPrimeStarDiagonalIncompressibleInternal() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pPrimeStar_D",D) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        input("pPrimeCoefficient") ;
        output("(cl,cr)->pPrimeStar_D") ;
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

  register_rule<PPrimeStarDiagonalIncompressibleInternal>
    registerPPrimeStarDiagonalIncompressibleInternal ;

  // Rule to assemble the diagonal term for the linear system. Assembling over
  // boundary faces. Note that this rule now applies for both incompressible
  // and compressible flows.  NOTE: Removed the nonSpecifiedMassFlux_BC constraint
  // so that this rule will be used for interpolated interface boundaries as
  // well where the pressure is specified via interpolation.
  class PPrimeStarDiagonalBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> pPrimeCoefficient ;
      store<real> D ;
    public:

      // Define input and output.
      PPrimeStarDiagonalBoundary() {
        name_store("ci",ci) ;
        name_store("pPrimeStar_D",D) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        input("pPrimeCoefficient") ;
        output("ci->pPrimeStar_D") ;
        constraint("specifiedPressure_BC") ;
      }

      // Distribute contribution for a single face.
      void calculate(Entity face) { D[ci[face]]+=pPrimeCoefficient[face] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeStarDiagonalBoundary> registerPPrimeStarDiagonalBoundary ;

  // Rule to assemble the diagonal term for the linear system for
  // compressible flow.
  class PPrimeStarDiagonalCompressibleInternal : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> rho ;
      const_store<real> soundSpeed ;
      const_store<real> massFluxStarHat ;
      store<real> D ;
    public:

      // Define input and output.
      PPrimeStarDiagonalCompressibleInternal() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("rho",rho) ;
        name_store("soundSpeed",soundSpeed) ;
        name_store("pPrimeStar_D",D) ;
        name_store("massFluxStarHat",massFluxStarHat) ;
        input("(cl,cr)->(rho,soundSpeed)") ;
        input("massFluxStarHat") ;
        output("(cl,cr)->pPrimeStar_D") ;
        constraint("internalFaces,compressibleFlow,compressiblePPrime") ;
      }

      // Distribute contribution for a single face.
      void calculate(Entity face) {
        if(massFluxStarHat[face]>0.0){
          real compressibleCoefficient=massFluxStarHat[face]/(rho[cl[face]]*
            soundSpeed[cl[face]]*soundSpeed[cl[face]]) ;
          D[cl[face]]+=compressibleCoefficient ;
        }else{
          real compressibleCoefficient=-massFluxStarHat[face]/(rho[cr[face]]*
            soundSpeed[cr[face]]*soundSpeed[cr[face]]) ;
          D[cr[face]]+=compressibleCoefficient ;
        }
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeStarDiagonalCompressibleInternal>
    registerPPrimeStarDiagonalCompressibleInternal ;

  // Rule to assemble the diagonal term for the linear system for
  // compressible flow.
  class PPrimeStarDiagonalCompressibleBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> rho ;
      const_store<real> soundSpeed ;
      const_store<real> massFlux ;
      store<real> D ;
    public:

      // Define input and output.
      PPrimeStarDiagonalCompressibleBoundary() {
        name_store("ci",ci) ;
        name_store("rho",rho) ;
        name_store("soundSpeed",soundSpeed) ;
        name_store("massFlux",massFlux) ;
        name_store("pPrimeStar_D",D) ;
        input("ci->(rho,soundSpeed),massFlux") ;
        output("ci->pPrimeStar_D") ;
        constraint("nonSpecifiedMassFlux_BC,compressibleFlow") ;
        constraint("compressiblePPrime") ;
      }

      // Distribute contribution for a single face.
      void calculate(Entity face) {
        if(massFlux[face]>0.0) D[ci[face]]+=massFlux[face]/(rho[ci[face]]*
            soundSpeed[ci[face]]*soundSpeed[ci[face]]) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeStarDiagonalCompressibleBoundary>
    registerPPrimeStarDiagonalCompressibleBoundary ;

  // Rule to add the unsteady component to the diagonal term for
  // compressible flow. NOTE: Leaving the time-step factor out for
  // now as it seems to cause problems.
  class PPrimeStarDiagonalUnsteadyCompressible : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> timeStep ;
      const_param<real> timeIntegratorFactor0 ;
      const_store<real> soundSpeed ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> D ;
    public:

      // Define input and output.
      PPrimeStarDiagonalUnsteadyCompressible() {
        name_store("timeStep",timeStep) ;
        name_store("timeIntegratorFactor0",timeIntegratorFactor0) ;
        name_store("soundSpeed",soundSpeed) ;
        name_store("vol",vol) ;
        name_store("cellRadius",cellRadius) ;
        name_store("pPrimeStar_D",D) ;
        input("timeStep,timeIntegratorFactor0,soundSpeed,vol,cellRadius") ;
        output("pPrimeStar_D") ;
        constraint("geom_cells,compressibleFlow") ;
      }

      // Distribute contribution for a single face.
      void calculate(Entity cell) {
        D[cell]+=0.5*((*timeIntegratorFactor0)+1.0)*vol[cell]/((*timeStep)*
          soundSpeed[cell]*soundSpeed[cell])*cellRadius[cell] ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeStarDiagonalUnsteadyCompressible>
    registerPPrimeStarDiagonalUnsteadyCompressible ;

  // Rule to compute the lower and upper terms for the linear system for
  // incompressible flow.
  class PPrimeStarLowerUpperIncompressible : public pointwise_rule {
    private:
      const_store<real> pPrimeCoefficient ;
      store<real> L,U ;
    public:

      // Define input and output.
      PPrimeStarLowerUpperIncompressible() {
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        name_store("pPrimeStar_L",L) ;
        name_store("pPrimeStar_U",U) ;
        input("pPrimeCoefficient") ;
        output("pPrimeStar_L,pPrimeStar_U") ;
        constraint("internalFaces") ;
      }

      // Initialize for a single face.
      void calculate(Entity face) {
        L[face]=-pPrimeCoefficient[face] ; U[face]=-pPrimeCoefficient[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeStarLowerUpperIncompressible>
    registerPPrimeStarLowerUpperIncompressible ;

  // Rule to compute the lower and upper terms for the linear system for
  // compressible flow.
  class PPrimeStarLowerUpperCompressible : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<real> rho ;
      const_store<real> soundSpeed ;
      const_store<real> pPrimeCoefficient ;
      const_store<real> massFluxStarHat ;
      store<real> L,U ;
    public:

      // Define input and output.
      PPrimeStarLowerUpperCompressible() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("rho",rho) ;
        name_store("soundSpeed",soundSpeed) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
//      name_store("massFluxStarHat",massFluxStarHat) ;
        name_store("massFluxCorrected_c_temp",massFluxStarHat) ;
        name_store("compressible::pPrimeStar_L",L) ;
        name_store("compressible::pPrimeStar_U",U) ;
        input("(cl,cr)->(rho,soundSpeed)") ;
//      input("pPrimeCoefficient,massFluxStarHat") ;
        input("pPrimeCoefficient,massFluxCorrected_c_temp") ;
        output("compressible::pPrimeStar_L,compressible::pPrimeStar_U") ;
        constraint("internalFaces,compressibleFlow,compressiblePPrime") ;
      }

      // Compute for a single face.
      void calculate(Entity face) {
        if(massFluxStarHat[face]>0.0){
          real compressibleCoefficient=massFluxStarHat[face]/(rho[cl[face]]*
            soundSpeed[cl[face]]*soundSpeed[cl[face]]) ;
          L[face]=-pPrimeCoefficient[face]-compressibleCoefficient ;
          U[face]=-pPrimeCoefficient[face] ;
        }else{
          real compressibleCoefficient=-massFluxStarHat[face]/(rho[cr[face]]*
            soundSpeed[cr[face]]*soundSpeed[cr[face]]) ;
          L[face]=-pPrimeCoefficient[face] ;
          U[face]=-pPrimeCoefficient[face]-compressibleCoefficient ;
        }
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeStarLowerUpperCompressible>
    registerPPrimeStarLowerUpperCompressible ;

  // Rule to initialize the rhs term for the linear system.
  class InitializePPrimeStarRHS : public unit_rule {
    private:
      store<real> B ;
    public:

      // Define input and output.
      InitializePPrimeStarRHS() {
        name_store("pPrimeStar_B",B) ;
        output("pPrimeStar_B") ;
        constraint("vol") ;
      }

      // Initialize for a single cell.
      void calculate(Entity cell) { B[cell]=0.0 ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializePPrimeStarRHS> registerInitializePPrimeStarRHS ;

  // Rule to assemble the rhs term for the linear system.
  class PPrimeStarRHSInternal : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> massFlux ;
      store<real> B ;
    public:

      // Define input and output.
      PPrimeStarRHSInternal() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("massFluxCorrected_c_temp",massFlux) ;
        name_store("pPrimeStar_B",B) ;
        input("massFluxCorrected_c_temp") ;
        output("(cl,cr)->pPrimeStar_B") ;
        constraint("internalFaces") ;
      }

      // Distribute contribution for a single face.
      void calculate(Entity face) {
        B[cl[face]]-=massFlux[face] ; B[cr[face]]+=massFlux[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeStarRHSInternal> registerPPrimeStarRHSInternal ;

  // Rule to added the lagged pressure-correction source term.
  class PPrimeStarRHSLaggedInternal : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<vect3d> cellCenter,pPrimeGradient,faceCenter ;
      const_store<Area> area ;
      const_store<real> pPrimeCoefficient ;
      store<real> B ;
    public:

      // Define input and output.
      PPrimeStarRHSLaggedInternal() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("cellcenter",cellCenter) ;
        name_store("grads(pPrimeStarOld)",pPrimeGradient) ;
        name_store("facecenter",faceCenter) ;
        name_store("area",area) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        name_store("pPrimeStar_B",B) ;
        input("(cl,cr)->(cellcenter,grads(pPrimeStarOld))") ;
        input("facecenter,area,pPrimeCoefficient") ;
        output("(cl,cr)->pPrimeStar_B") ;
        constraint("internalFaces") ;
      }

      // Distribute contribution for a single face.
      void calculate(Entity face) {
        vect3d rL=faceCenter[face]-dot(faceCenter[face]-cellCenter[cl[face]],
          area[face].n)*area[face].n ;
        vect3d rR=faceCenter[face]-dot(faceCenter[face]-cellCenter[cr[face]],
          area[face].n)*area[face].n ;
        real laggedSourceTerm=pPrimeCoefficient[face]*(dot(pPrimeGradient[cl[face]],
          rL-cellCenter[cl[face]])-dot(pPrimeGradient[cr[face]],rR-cellCenter[cr[face]])) ;
        B[cl[face]]-=laggedSourceTerm ; B[cr[face]]+=laggedSourceTerm ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeStarRHSLaggedInternal> registerPPrimeStarRHSLaggedInternal ;

  // Rule to assemble the rhs term for the linear system.
  class PPrimeStarRHSBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> massFlux ;
      store<real> B ;
    public:

      // Define input and output.
      PPrimeStarRHSBoundary() {
        name_store("ci",ci) ;
        name_store("massFlux",massFlux) ;
        name_store("pPrimeStar_B",B) ;
        input("massFlux") ;
        output("ci->pPrimeStar_B") ;
        constraint("boundaryFaces") ;
      }

      // Compute for a single face.
      void calculate(Entity face) { B[ci[face]]-=massFlux[face] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeStarRHSBoundary> registerPPrimeStarRHSBoundary ;

  // Rule to add the unsteady component to the right-hand side term.
  // This rule will have to be updated for deforming meshes.
  class PPrimeStarRHSUnsteadyBDF : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> timeStep ;
      const_store<real> rho,rhoOld ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> B ;
    public:

      // Define input and output. Note the {n,it,ips} on the compressible
      // flow constraint. Without this, Loci complains that this should be
      // a build rule.
      PPrimeStarRHSUnsteadyBDF() {
        name_store("timeStep{n}",timeStep) ;
        name_store("rho{n}",rhoOld) ;
//      name_store("rho{n,it}",rho) ;
        name_store("rhoCorrected_c{n,it,ips}",rho) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n}",cellRadius) ;
        name_store("pPrimeStar_B{n,it,ips}",B) ;
//      input("timeStep{n},rho{n},rho{n,it},vol{n},cellRadius{n}") ;
        input("timeStep{n},rho{n},rhoCorrected_c{n,it,ips},vol{n},cellRadius{n}") ;
        output("pPrimeStar_B{n,it,ips}") ;
        constraint("compressibleFlow{n,it,ips},BDFIntegrator,vol") ;
      }

      // Compute contribution for a single cell.
      void calculate(Entity cell) {
        B[cell]+=vol[cell]*cellRadius[cell]*(rhoOld[cell]-rho[cell])/
          (*timeStep) ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeStarRHSUnsteadyBDF> registerPPrimeStarRHSUnsteadyBDF ;

  // Rule to add the unsteady component to the right-hand side term.
  // This rule will have to be updated for deforming meshes.
  class PPrimeStarRHSUnsteady : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<int> n ;
      const_param<real> timeStep ;
      const_store<real> rho,rhoOld,rhoOlder ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> B ;
    public:

      // Define input and output. Note the {n,it,ips} on the compressible
      // flow constraint. Without this, Loci complains that this should be
      // a build rule.
      PPrimeStarRHSUnsteady() {
        name_store("$n{n}",n) ;
        name_store("timeStep{n}",timeStep) ;
        name_store("rho{n-1}",rhoOlder) ;
        name_store("rho{n}",rhoOld) ;
        name_store("rho{n,it}",rho) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n}",cellRadius) ;
        name_store("pPrimeStar_B{n,it,ips}",B) ;
        input("$n{n},timeStep{n},rho{n-1},rho{n},rho{n,it}") ;
        input("vol{n},cellRadius{n}") ;
        output("pPrimeStar_B{n,it,ips}") ;
        constraint("compressibleFlow{n,it,ips},BDF2Integrator,vol") ;
      }

      // Compute contribution for a single cell. Use BDF on first timestep.
      void calculate(Entity cell) {
        if((*n)!=0){
          B[cell]+=0.5*vol[cell]*cellRadius[cell]*(4.0*rhoOld[cell]-3.0*
            rho[cell]-rhoOlder[cell])/(*timeStep) ;
        }else{
          B[cell]+=vol[cell]*cellRadius[cell]*(rhoOld[cell]-rho[cell])/
            (*timeStep) ;
        }
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PPrimeStarRHSUnsteady> registerPPrimeStarRHSUnsteady ;

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

//-----------------------------------------------------------------------------
// Rules for computing the residual of the pPrimeStar equation which is solved
// in the corrector stage.

  // Create a variable for the residual of the pPrimeStar equation for use with
  // output rules.
  class PPrimeStarResidual : public pointwise_rule {
    private:
      const_param<real> rhoScale,vScale,lScale ;
      const_store<real> B ;
      store<real> pResidual ;
    private:
      real pPrimeFactor ;
    public:

      // Define input and output.
      PPrimeStarResidual() {
        name_store("rhoScale",rhoScale) ;
        name_store("vScale",vScale) ;
        name_store("lScale",lScale) ;
        name_store("pPrimeStar_B",B) ;
        name_store("pPrimeStarResidual",pResidual) ;
        input("rhoScale,vScale,lScale,pPrimeStar_B") ;
        output("pPrimeStarResidual") ;
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

  register_rule<PPrimeStarResidual> registerPPrimeStarResidual ;

  // Rule to initialize the total pPrimeStar residual.
  class InitializeTotalPPrimeStarResidual : public unit_rule {
    private:
      param<ScalarResidual> pPrimeStarResidualData ;
    public:

      // Define input and output.
      InitializeTotalPPrimeStarResidual() {
        name_store("pPrimeStarResidualDataNew",pPrimeStarResidualData) ;
        output("pPrimeStarResidualDataNew") ;
        constraint("geom_cells") ;
      }

      // Initialize the residual.
      virtual void compute(const sequence &seq) {
        *pPrimeStarResidualData=ScalarResidual() ;
      }
  } ;

  register_rule<InitializeTotalPPrimeStarResidual>
    registerInitializeTotalPPrimeStarResidual ;

  // Rule to compute the total pPrimeStar residual.
  class ComputeTotalPPrimeStarResidual : public apply_rule
  <param<ScalarResidual>,ScalarResidualJoin> {
    private:
      const_store<real> B ;
      const_store<vect3d> cellCenter ;
      param<ScalarResidual> pPrimeStarResidualData ;
    public:

      // Define input and output.
      ComputeTotalPPrimeStarResidual() {
        name_store("pPrimeStar_B",B) ;
        name_store("cellcenter",cellCenter) ;
        name_store("pPrimeStarResidualDataNew",pPrimeStarResidualData) ;
        input("pPrimeStar_B,cellcenter") ;
        output("pPrimeStarResidualDataNew") ;
        constraint("geom_cells") ;
      }

      // Add the cell contribution to the residual for a single cell.
      void calculate(Entity cell) {
        ScalarResidual temp ; temp.maxResidual=B[cell] ;
        temp.totalResidual=abs(B[cell]) ;
        temp.maxResidualLocation=cellCenter[cell] ;
        join(*pPrimeStarResidualData,temp) ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeTotalPPrimeStarResidual>
    registerComputeTotalPPrimeStarResidual ;

//-----------------------------------------------------------------------------
// Rules for computing velocity and mass flux corrections.

  // Rule to compute the density correction which is required for the deferred-
  // correction loop for pressure correction.
  class DensityCorrection : public pointwise_rule {
    private:
      const_store<real> pPrime,soundSpeed ;
      store<real> rhoPrime ;
    public:

      // Define input and output.
      DensityCorrection() {
        name_store("X",pPrime) ;
        name_store("soundSpeed",soundSpeed) ;
        name_store("rhoPrime(X)",rhoPrime) ;
        input("X,soundSpeed") ;
        output("rhoPrime(X)") ;
        constraint("vol") ;
      }

      // Set value for a single cell.
      void calculate(Entity cell) {
        rhoPrime[cell]=pPrime[cell]/(soundSpeed[cell]*soundSpeed[cell]) ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DensityCorrection> registerDensityCorrection ;

  // Rule to initialize the velocity correction. This is now a parametric
  // variable with X representing the pressure correction.
  class InitializeVelocityCorrection : public unit_rule {
    private:
      store<vect3d> vPrime ;
    public:

      // Define input and output.
      InitializeVelocityCorrection() {
        name_store("vPrime(X)",vPrime) ;
        output("vPrime(X)") ;
        constraint("vol,X") ;
      }

      // Initialize velocity correction for a single cell.
      void calculate(Entity cell) { vPrime[cell]=vect3d(0.0,0.0,0.0) ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeVelocityCorrection>
    registerInitializeVelocityCorrection ;

  // Rule to add cell velocity corrections from interior faces.
  class VelocityCorrectionInterior : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_param<real> twoDimensionFactor ;
      const_Map cl,cr ;
      const_store<real> pPrime ;
      const_store<real> vMainCoefficient ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<vect3d> vPrime ;
    public:

      // Define input and output.
      VelocityCorrectionInterior() {
        name_store("twoDimensionFactor",twoDimensionFactor) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("X",pPrime) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("vPrime(X)",vPrime) ;
        input("twoDimensionFactor,area,faceRadius") ;
        input("(cl,cr)->(X,vMainCoefficient)") ;
        output("(cl,cr)->vPrime(X)") ;
        constraint("internalFaces") ;
      }

      // Correct the velocity for the cells attach to a single face.
      void calculate(Entity face) {
        vect3d temp=0.5*(pPrime[cl[face]]+pPrime[cr[face]])*area[face].n*
          area[face].sada*faceRadius[face] ;
        temp.z*=(*twoDimensionFactor) ;
        vPrime[cl[face]]-=temp/vMainCoefficient[cl[face]] ;
        vPrime[cr[face]]+=temp/vMainCoefficient[cr[face]] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<VelocityCorrectionInterior>
    registerVelocityCorrectionInterior ;

  // Rule to add the higher-order components of the pressure correction
  // to the velocity correction.
  class VelocityCorrectionHighOrderInterior : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
      const_param<real> twoDimensionFactor ;
      const_param<real> thetaParameter ;
      const_store<real> faceRadius ;
      const_store<vect3d> pPrimeGradient ;
      const_store<real> pPrimeLimiter ;
      const_store<real> vMainCoefficient ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<Area> area ;
      store<vect3d> vPrime ;
    public:

      // Define input and output.
      VelocityCorrectionHighOrderInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("twoDimensionFactor",twoDimensionFactor) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("faceRadius",faceRadius) ;
        name_store("grads(X)",pPrimeGradient) ;
        name_store("limiters(X)",pPrimeLimiter) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("area",area) ;
        name_store("vPrime(X)",vPrime) ;
        input("twoDimensionFactor,thetaParameter,faceRadius") ;
        input("(cl,cr)->(grads(X),limiters(X),vMainCoefficient,cellcenter)") ;
        input("facecenter,area") ;
        output("(cl,cr)->vPrime(X)") ;
        constraint("internalFaces") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        vect3d temp=0.5*(pPrimeLimiter[cl[face]]*dot(pPrimeGradient[cl[face]],
          faceCenter[face]-cellCenter[cl[face]])+pPrimeLimiter[cr[face]]*
          dot(pPrimeGradient[cr[face]],faceCenter[face]-cellCenter[cr[face]]))*
          area[face].n*area[face].sada*faceRadius[face] ;
        temp.z*=(*twoDimensionFactor) ;
        vPrime[cl[face]]-=temp/vMainCoefficient[cl[face]] ;
        vPrime[cr[face]]+=temp/vMainCoefficient[cr[face]] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

//register_rule<VelocityCorrectionHighOrderInterior>
//  registerVelocityCorrectionHighOrderInterior ;

  // Rule to add cell velocity corrections from boundary faces.
  class VelocityCorrectionBoundary : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_param<real> twoDimensionFactor ;
      const_Map ci ;
      const_store<real> vMainCoefficient ;
      const_store<real> pPrime ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<vect3d> vPrime ;
    public:

      // Define input and output.
      VelocityCorrectionBoundary() {
        name_store("twoDimensionFactor",twoDimensionFactor) ;
        name_store("ci",ci) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        name_store("X",pPrime) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("vPrime(X)",vPrime) ;
        input("twoDimensionFactor,ci->(X,vMainCoefficient),area,faceRadius") ;
        output("ci->vPrime(X)") ;
      }

      // Correct the velocity for the cells attach to a single face.
      void calculate(Entity face) {
        vect3d temp=(pPrime[ci[face]]/vMainCoefficient[ci[face]])*
          area[face].n*area[face].sada*faceRadius[face] ;
        temp.z*=(*twoDimensionFactor) ;
        vPrime[ci[face]]-=temp ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<VelocityCorrectionBoundary> registerVelocityCorrectionBoundary ;

  // Rule to add cell velocity corrections from boundary faces.
  class VelocityCorrectionHighOrderBoundary : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_param<real> twoDimensionFactor ;
      const_Map ci ;
      const_store<vect3d> pPrimeGradient ;
      const_store<real> pPrimeLimiter ;
      const_store<real> vMainCoefficient ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<vect3d> vPrime ;
    public:

      // Define input and output.
      VelocityCorrectionHighOrderBoundary() {
        name_store("twoDimensionFactor",twoDimensionFactor) ;
        name_store("ci",ci) ;
        name_store("grads(X)",pPrimeGradient) ;
        name_store("limiters(X)",pPrimeLimiter) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("vPrime(X)",vPrime) ;
        input("ci->(grads(X),limiters(X),vMainCoefficient,cellcenter)") ;
        input("facecenter,area,faceRadius,twoDimensionFactor") ;
        output("ci->vPrime(X)") ;
      }

      // Correct the velocity for the cells attach to a single face.
      void calculate(Entity face) {
        vect3d temp=(pPrimeLimiter[ci[face]]*dot(pPrimeGradient[ci[face]],
          faceCenter[face]-cellCenter[ci[face]])/vMainCoefficient[ci[face]]*
          area[face].sada*faceRadius[face])*area[face].n ;
        temp.z*=(*twoDimensionFactor) ;
        vPrime[ci[face]]-=temp ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

//register_rule<VelocityCorrectionHighOrderBoundary>
//  registerVelocityCorrectionHighOrderBoundary ;

  // Parametric rule for corrected mass flux on interior faces.
  class MassFluxCorrectionInterior : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<vect3d> cellCenter ;
      const_store<real> pPrime ;
      const_store<vect3d> pPrimeGradient ;
      const_store<vect3d> faceCenter ;
      const_store<Area> area ;
      const_store<real> pPrimeCoefficient ;
      store<real> massFluxPrime ;
    public:

      // Define input and output.
      MassFluxCorrectionInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("cellcenter",cellCenter) ;
        name_store("X",pPrime) ;
        name_store("grads(Y)",pPrimeGradient) ;
        name_store("facecenter",faceCenter) ;
        name_store("area",area) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        name_store("massFluxPrime(X,Y)",massFluxPrime) ;
        input("(cl,cr)->(cellcenter,X,grads(Y))") ;
        input("facecenter,area,pPrimeCoefficient") ;
        output("massFluxPrime(X,Y)") ;
        constraint("internalFaces") ;
      }

      // Compute corrected mass flux for a single face.
      void calculate(Entity face) {
        vect3d rL=faceCenter[face]-dot(faceCenter[face]-cellCenter[cl[face]],
          area[face].n)*area[face].n ;
        vect3d rR=faceCenter[face]-dot(faceCenter[face]-cellCenter[cr[face]],
          area[face].n)*area[face].n ;
        massFluxPrime[face]=pPrimeCoefficient[face]*(pPrime[cl[face]]-
          pPrime[cr[face]]+dot(pPrimeGradient[cl[face]],rL-cellCenter[cl[face]])-
          dot(pPrimeGradient[cr[face]],rR-cellCenter[cr[face]])) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<MassFluxCorrectionInterior> registerMassFluxCorrectionInterior ;

  // Rule for computing mass flux correction on boundary faces where mass
  // flux is iterating.
  class MassFluxCorrectionBoundary : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<vect3d> cellCenter ;
      const_store<real> pPrime ;
      const_store<vect3d> pPrimeGradient ;
      const_store<vect3d> faceCenter ;
      const_store<Area> area ;
      const_store<real> pPrimeCoefficient ;
      store<real> massFluxPrime ;
    public:

      // Define input and output.
      MassFluxCorrectionBoundary() {
        name_store("ci",ci) ;
        name_store("cellcenter",cellCenter) ;
        name_store("X",pPrime) ;
        name_store("grads(Y)",pPrimeGradient) ;
        name_store("facecenter",faceCenter) ;
        name_store("area",area) ;
        name_store("pPrimeCoefficient",pPrimeCoefficient) ;
        name_store("massFluxPrime(X,Y)",massFluxPrime) ;
        input("ci->(cellcenter,X,grads(Y))") ;
        input("facecenter,area,pPrimeCoefficient") ;
        output("massFluxPrime(X,Y)") ;
        constraint("boundaryMassFluxCorrected") ;
      }

      // Compute mass flux correction for a single face.
      void calculate(Entity face) {
        vect3d r=faceCenter[face]-dot(faceCenter[face]-cellCenter[ci[face]],
          area[face].n)*area[face].n ;
        massFluxPrime[face]=pPrimeCoefficient[face]*(pPrime[ci[face]]+
          dot(pPrimeGradient[ci[face]],r-cellCenter[ci[face]])) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<MassFluxCorrectionBoundary>
    registerCorrectMassFluxBoundary ;
}
