//-----------------------------------------------------------------------------
// Description: This file contains rules for Menter's family of k-omega
//   turbulence models.
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

namespace streamUns {

//-----------------------------------------------------------------------------
// Pre-computed quantities.

  // Multiplier to allow for DES turbulence models. Here set to one as default.
  class FDESDefault : public pointwise_rule {
    private:
      store<real> fDES ;
    public:

      // Define input and output.
      FDESDefault() {
        name_store("fDES",fDES) ;
        output("fDES") ;
        constraint("geom_cells") ;
      }

      // Compute for each cell.
      void calculate(Entity cell) { fDES[cell]=1.0 ; }
                                                                                
      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FDESDefault> registerFDESDefault ;

  // TKE viscosity for the interior faces.
  class KViscosityInterior : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_param<real> sigmaK1,sigmaK2 ;
      const_store<real> laminarViscosity,eddyViscosity ;
      const_store<real> f1 ;
      store<real> kViscosity ;
    public:

      // Define input and output.
      KViscosityInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("sigmaK1",sigmaK1) ;
        name_store("sigmaK2",sigmaK2) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("f1",f1) ;
        name_store("kViscosity",kViscosity) ;
        input("sigmaK1,sigmaK2") ;
        input("(cl,cr)->(laminarViscosity,eddyViscosity,f1)") ;
        output("kViscosity") ;
        constraint("internalFaces") ;
      }

      // Compute for each face.
      void calculate(Entity face) {
        real sigmaKLeft=f1[cl[face]]*(*sigmaK1)+(1.0-f1[cl[face]])*(*sigmaK2) ;
        real sigmaKRight=f1[cr[face]]*(*sigmaK1)+(1.0-f1[cr[face]])*(*sigmaK2) ;
        kViscosity[face]=0.5*(laminarViscosity[cl[face]]+sigmaKLeft*
          eddyViscosity[cl[face]]+laminarViscosity[cr[face]]+sigmaKRight*
          eddyViscosity[cr[face]]) ;
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<KViscosityInterior> registerKViscosityInterior ;

  // Omega viscosity for the interior faces.
  class OmegaViscosityInterior : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_param<real> sigmaOmega1,sigmaOmega2 ;
      const_store<real> laminarViscosity,eddyViscosity ;
      const_store<real> f1 ;
      store<real> omegaViscosity ;
    public:

      // Define input and output.
      OmegaViscosityInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("sigmaOmega1",sigmaOmega1) ;
        name_store("sigmaOmega2",sigmaOmega2) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("omegaViscosity",omegaViscosity) ;
        name_store("f1",f1) ;
        input("sigmaOmega1,sigmaOmega2") ;
        input("(cl,cr)->(laminarViscosity,eddyViscosity,f1)") ;
        output("omegaViscosity") ;
        constraint("internalFaces") ;
      }

      // Compute for each face.
      void calculate(Entity face) {
        real sigmaOmegaLeft=f1[cl[face]]*(*sigmaOmega1)+(1.0-f1[cl[face]])*
          (*sigmaOmega2) ;
        real sigmaOmegaRight=f1[cr[face]]*(*sigmaOmega1)+(1.0-f1[cr[face]])*
          (*sigmaOmega2) ;
        omegaViscosity[face]=0.5*(laminarViscosity[cl[face]]+sigmaOmegaLeft*
          eddyViscosity[cl[face]]+laminarViscosity[cr[face]]+sigmaOmegaRight*
          eddyViscosity[cr[face]]) ;
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<OmegaViscosityInterior> registerOmegaViscosityInterior ;

//-----------------------------------------------------------------------------
// Rules to process turbulence equation options from the .vars file.

  // Creates the turbulence equation solver constraints.
  class TurbulenceEquationSolverConstraints : public constraint_rule {
    private:
      const_param<TurbulenceEquationOptions> turbulenceEquationOptions ;
      Constraint turbulenceSGSLinearSolver ;
    public:
                                                                                
      // Define input and output.
      TurbulenceEquationSolverConstraints() {
        name_store("turbulenceEquationOptions",turbulenceEquationOptions) ;
        name_store("turbulenceSGSLinearSolver",turbulenceSGSLinearSolver) ;
        input("turbulenceEquationOptions") ;
        output("turbulenceSGSLinearSolver") ;
      }
                                                                                
      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        if((*turbulenceEquationOptions).optionExists("linearSolver")){
          Loci::option_value_type optionValueType=turbulenceEquationOptions->
            getOptionValueType("linearSolver") ;
          switch(optionValueType){
            case Loci::NAME:
              {
                Loci::option_values optionValues=turbulenceEquationOptions->
                  getOption("linearSolver") ;
                string name ; optionValues.get_value(name) ;
                if(name=="SGS"){
                  turbulenceSGSLinearSolver=~EMPTY ;
                }else{
                  cerr << "Bad linearSolver for turbulenceEquation." << endl ;
                  Loci::Abort() ;
                }
              }
              break ;
            default:
              cerr << "Bad type for linearSolver in turbulenceEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          turbulenceSGSLinearSolver=~EMPTY ;
        }
      }
                                                                                
  } ;
                                                                                
  register_rule<TurbulenceEquationSolverConstraints>
    registerTurbulenceEquationSolverConstraints ;
                                                                                
  // Creates the constraints for the turbulence equations.
  class TurbulenceModelConstraints : public constraint_rule {
    private:
      const_param<TurbulenceEquationOptions> turbulenceEquationOptions ;
      Constraint kOmegaTurbulenceModel,menterTurbulenceModel ;
      Constraint menterBSLTurbulenceModel,menterSSTTurbulenceModel ;
      Constraint menterBSLSSTTurbulenceModel,menterSST2003TurbulenceModel ;
    public:
                                                                                
      // Define input and output.
      TurbulenceModelConstraints() {
        name_store("turbulenceEquationOptions",turbulenceEquationOptions) ;
        name_store("kOmegaTurbulenceModel",kOmegaTurbulenceModel) ;
        name_store("menterTurbulenceModel",menterTurbulenceModel) ;
        name_store("menterBSLTurbulenceModel",menterBSLTurbulenceModel) ;
        name_store("menterSSTTurbulenceModel",menterSSTTurbulenceModel) ;
        name_store("menterBSLSSTTurbulenceModel",menterBSLSSTTurbulenceModel) ;
        name_store("menterSST2003TurbulenceModel",
          menterSST2003TurbulenceModel) ;
        input("turbulenceEquationOptions") ;
        output("kOmegaTurbulenceModel,menterTurbulenceModel") ;
        output("menterBSLTurbulenceModel,menterSSTTurbulenceModel") ;
        output("menterBSLSSTTurbulenceModel,menterSST2003TurbulenceModel") ;
      }
                                                                                
      // Set up the constraints.
      virtual void compute(const sequence& seq) {

        if((*turbulenceEquationOptions).optionExists("model")){
          Loci::option_value_type optionValueType=turbulenceEquationOptions->
            getOptionValueType("model") ;
          switch(optionValueType){
            case Loci::NAME:
              {
                Loci::option_values optionValues=turbulenceEquationOptions->
                   getOption("model") ;
                string modelName ; optionValues.get_value(modelName) ;
                if(modelName=="menterBSL" || modelName=="menterSST" ||
                modelName=="menterSST2003"){
                  kOmegaTurbulenceModel=~EMPTY ; menterTurbulenceModel=~EMPTY ;
                  if(modelName=="menterBSL"){
                    menterBSLTurbulenceModel=~EMPTY ;
                    menterSSTTurbulenceModel=EMPTY ;
                    menterBSLSSTTurbulenceModel=~EMPTY ;
                    menterSST2003TurbulenceModel=EMPTY ;
                  }else if(modelName=="menterSST"){
                    menterBSLTurbulenceModel=EMPTY ;
                    menterSSTTurbulenceModel=~EMPTY ;
                    menterBSLSSTTurbulenceModel=~EMPTY ;
                    menterSST2003TurbulenceModel=EMPTY ;
                  }else{
                    menterBSLTurbulenceModel=EMPTY ;
                    menterSSTTurbulenceModel=EMPTY ;
                    menterBSLSSTTurbulenceModel=EMPTY ;
                    menterSST2003TurbulenceModel=~EMPTY ;
                  }
                }else{
                  cerr << "ERROR: Bad model for turbulenceEquation." << endl ;
                  Loci::Abort() ;
                }
              }
              break ;
            default:
              cerr << "Bad type for model in turbulenceEquation." << endl ;
              Loci::Abort() ;
          }
        }else{
          cerr << "ERROR: Must supply a model in turbulenceEquation." << endl ;
          Loci::Abort() ;
        }
      }
  } ;
                                                                                
  register_rule<TurbulenceModelConstraints>
    registerTurbulenceModelConstraints ;

  // Creates the turbulence equation solver parameters.
  class TurbulenceEquationSolverParameters : public singleton_rule {
    private:
      const_param<TurbulenceEquationOptions> turbulenceEquationOptions ;
      param<int> turbulenceMaxIterations ;
    public:
                                                                                
      // Define input and output.
      TurbulenceEquationSolverParameters() {
        name_store("turbulenceEquationOptions",turbulenceEquationOptions) ;
        name_store("turbulenceMaxIterations",turbulenceMaxIterations) ;
        input("turbulenceEquationOptions") ;
        output("turbulenceMaxIterations") ;
      }
                                                                                
      // Set up the parameters.
      virtual void compute(const sequence& seq) {

        // Maximum number of iterations.
        if((*turbulenceEquationOptions).optionExists("maxIterations")){
          Loci::option_value_type optionValueType=turbulenceEquationOptions->
            getOptionValueType("maxIterations") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                turbulenceEquationOptions->getOption("maxIterations",temp) ;
                if(int(temp)<0){
                  cerr << "Bad maxIterations value for turbulenceEquation."
                    << endl ; Loci::Abort() ;
                }
                *turbulenceMaxIterations=int(temp) ;
              }
              break ;
            default:
              cerr << "Bad type for maxIterations in turbulenceEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *turbulenceMaxIterations=5 ;
        }
      }
  } ;

  register_rule<TurbulenceEquationSolverParameters>
    registerTurbulenceEquationSolverParamters ;

//-----------------------------------------------------------------------------
// New rules for a noWallFunction multiplication factor. This factor has been
// created so that we can get rid of the noWallFunction constraint, which is
// hard to implement since we can no longer have code in the boundary condition
// setup function.

  class NoWallFunctionDefault : public pointwise_rule {
    private:
      store<real> noWallFunction ;
    public:

      // Define input and output.
      NoWallFunctionDefault() {
        name_store("noWallFunction",noWallFunction) ;
        output("noWallFunction") ;
        constraint("boundaryFaces") ;
      }

      // Set default value.
      void calculate (Entity face) { noWallFunction[face]=1.0 ; }

      // Loop over faces.
      void compute (const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<NoWallFunctionDefault> registerNoWallFunctionDefault ;

  // Priority rule sets multiplier to zero for wallFunction faces.
  class NoWallFunctionPriority : public pointwise_rule {
    private:
      store<real> noWallFunction ;
    public:

      // Define input and output.
      NoWallFunctionPriority() {
        name_store("priority::noWallFunction",noWallFunction) ;
        output("priority::noWallFunction") ;
        constraint("ref->wallFunction_BCoption") ;
      }

      // Set default value.
      void calculate (Entity face) { noWallFunction[face]=0.0 ; }

      // Loop over faces.
      void compute (const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<NoWallFunctionPriority> registerNoWallFunctionPriority ;

//-----------------------------------------------------------------------------
// Ed's rule to compute the distance from the cell to the closest no-slip face.

  class DistanceToNoSlip : public pointwise_rule {
    private:
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<Area> area ;
      const_Map minCellToNoSlip ;
      store<real> distanceToNoSlip ;
    public:

      // Define input and output.
      DistanceToNoSlip() {
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("min_cell2noslip",minCellToNoSlip) ;
        name_store("priority::distanceToNoSlip",distanceToNoSlip) ;
        name_store("area",area) ;
        input("min_cell2noslip->area") ;
        input("cellcenter,min_cell2noslip->facecenter") ;
        output("priority::distanceToNoSlip") ;
      }

      // Compute the distance for a single cell.
      void calculate (Entity cell) {
        distanceToNoSlip[cell]=abs(dot(faceCenter[minCellToNoSlip[cell]]-
          cellCenter[cell],area[minCellToNoSlip[cell]].n)) ;
      }

      // Compute the distance for a sequence of cells.
      void compute (const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DistanceToNoSlip> registerDistanceToNoSlip ;

  // Default rule to provide a large value when there are no noslip surfaces.
  // Did not used to have this rule, and when one solved a freestream
  // problem, there was no destruction term in omega equation.
  class DefaultDistanceToNoSlip : public pointwise_rule {
    private:
      store<real> distanceToNoSlip ;
    public:

      // Define input and output.
      DefaultDistanceToNoSlip() {
        name_store("distanceToNoSlip",distanceToNoSlip) ;
        output("distanceToNoSlip") ;
        constraint("UNIVERSE") ;
      }

      // Set distance to large value.
      void calculate (Entity cell) { distanceToNoSlip[cell]=1.0e+30 ; }

      // Compute the distance for a sequence of cells.
      void compute (const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DefaultDistanceToNoSlip> registerDefaultDistanceToNoSlip ;

//-----------------------------------------------------------------------------
// Rules for the constants for the Baseline and SST models. Checked.

  class BaselineModelConstants : public singleton_rule {
    private:
      param<real> betaStar,beta1,beta2 ;
      param<real> sigmaK1,sigmaK2,sigmaOmega1,sigmaOmega2 ;
      param<real> gamma1,gamma2 ;
    public:

      // Define input and output.
      BaselineModelConstants() {
        name_store("betaStar",betaStar) ;
        name_store("beta1",beta1) ;
        name_store("beta2",beta2) ;
        name_store("sigmaK1",sigmaK1) ;
        name_store("sigmaK2",sigmaK2) ;
        name_store("sigmaOmega1",sigmaOmega1) ;
        name_store("sigmaOmega2",sigmaOmega2) ;
        name_store("gamma1",gamma1) ;
        name_store("gamma2",gamma2) ;
        output("betaStar,beta1,beta2,sigmaK1,sigmaK2,sigmaOmega1,sigmaOmega2") ;
        output("gamma1,gamma2") ;
        constraint("menterBSLTurbulenceModel") ;
      }

      // Compute the model constants.
      virtual void compute(const sequence &seq) {
        *betaStar=0.09 ; *beta1=0.0750 ; *beta2=0.0828 ;
        *sigmaK1=0.5 ; *sigmaK2=1.0 ; *sigmaOmega1=0.5 ; *sigmaOmega2=0.856 ;
        *gamma1=(*beta1)/(*betaStar)-(*sigmaOmega1)*0.41*0.41/sqrt(*betaStar) ;
        *gamma2=(*beta2)/(*betaStar)-(*sigmaOmega2)*0.41*0.41/sqrt(*betaStar) ;
      }
  } ;

  register_rule<BaselineModelConstants> registerBaselineModelConstants ;

  class SSTModelConstants : public singleton_rule {
    private:
      param<real> betaStar,beta1,beta2 ;
      param<real> sigmaK1,sigmaK2,sigmaOmega1,sigmaOmega2 ;
      param<real> gamma1,gamma2 ;
    public:

      // Define input and output.
      SSTModelConstants() {
        name_store("betaStar",betaStar) ;
        name_store("beta1",beta1) ;
        name_store("beta2",beta2) ;
        name_store("sigmaK1",sigmaK1) ;
        name_store("sigmaK2",sigmaK2) ;
        name_store("sigmaOmega1",sigmaOmega1) ;
        name_store("sigmaOmega2",sigmaOmega2) ;
        name_store("gamma1",gamma1) ;
        name_store("gamma2",gamma2) ;
        output("betaStar,beta1,beta2,sigmaK1,sigmaK2,sigmaOmega1,sigmaOmega2") ;
        output("gamma1,gamma2") ;
        constraint("menterSSTTurbulenceModel") ;
      }

      // Compute the model constants.
      virtual void compute(const sequence &seq) {
        *betaStar=0.09 ; *beta1=0.0750 ; *beta2=0.0828 ;
        *sigmaK1=0.85 ; *sigmaK2=1.0 ; *sigmaOmega1=0.5 ; *sigmaOmega2=0.856 ;
        *gamma1=(*beta1)/(*betaStar)-(*sigmaOmega1)*0.41*0.41/sqrt(*betaStar) ;
        *gamma2=(*beta2)/(*betaStar)-(*sigmaOmega2)*0.41*0.41/sqrt(*betaStar) ;
      }
  } ;

  register_rule<SSTModelConstants> registerSSTModelConstants ;

  class SST2003ModelConstants : public singleton_rule {
    private:
      param<real> betaStar,beta1,beta2 ;
      param<real> sigmaK1,sigmaK2,sigmaOmega1,sigmaOmega2 ;
      param<real> gamma1,gamma2 ;
    public:

      // Define input and output.
      SST2003ModelConstants() {
        name_store("betaStar",betaStar) ;
        name_store("beta1",beta1) ;
        name_store("beta2",beta2) ;
        name_store("sigmaK1",sigmaK1) ;
        name_store("sigmaK2",sigmaK2) ;
        name_store("sigmaOmega1",sigmaOmega1) ;
        name_store("sigmaOmega2",sigmaOmega2) ;
        name_store("gamma1",gamma1) ;
        name_store("gamma2",gamma2) ;
        output("betaStar,beta1,beta2,sigmaK1,sigmaK2,sigmaOmega1,sigmaOmega2") ;
        output("gamma1,gamma2") ;
        constraint("menterSST2003TurbulenceModel") ;
      }

      // Compute the model constants.
      virtual void compute(const sequence &seq) {
        *betaStar=0.09 ; *beta1=0.0750 ; *beta2=0.0828 ;
        *sigmaK1=0.85 ; *sigmaK2=1.0 ; *sigmaOmega1=0.5 ; *sigmaOmega2=0.856 ;
        *gamma1=5.0/9.0 ; *gamma2=0.44 ;
      }
  } ;

  register_rule<SST2003ModelConstants> registerSST2003ModelConstants ;

//-----------------------------------------------------------------------------
// Rules for the computing the eddy viscosity for the Baseline and SST models.

  // Rule for the baseline model to assign eddy viscosity for cells. Checked.
  class EddyViscosityBaselineInterior : public pointwise_rule {
    private:
      const_param<real> eddyViscosityLimit ;
      const_store<real> rho ;
      const_store<real> k ;
      const_store<real> omega ;
      const_store<real> laminarViscosity ;
      store<real> eddyViscosity ;
    public:

      // Define input and output.
      EddyViscosityBaselineInterior() {
        name_store("eddyViscosityLimit",eddyViscosityLimit) ;
        name_store("rho",rho) ;
        name_store("k",k) ;
        name_store("omega",omega) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("eddyViscosity",eddyViscosity) ;
        input("eddyViscosityLimit,rho,k,omega,laminarViscosity") ;
        output("eddyViscosity") ;
        constraint("vol,menterBSLTurbulenceModel") ;
      }

      // Set eddy viscosity to zero for a single cell.
      void calculate(Entity cell) {
        eddyViscosity[cell]=min(rho[cell]*k[cell]/omega[cell],
          *eddyViscosityLimit*laminarViscosity[cell]) ;
      }

      // Set eddy viscosity to zero for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<EddyViscosityBaselineInterior>
    registerEddyViscosityBaselineInterior ;

  // Rule for the baseline model to assign eddy viscosity for boundary faces.
  // Checked.
  class EddyViscosityBaselineBoundary : public pointwise_rule {
    private:
      const_param<real> eddyViscosityLimit ;
      const_store<real> rho_f ;
      const_store<real> k_f ;
      const_store<real> omega_f ;
      const_store<real> laminarViscosity_f ;
      store<real> eddyViscosity_f ;
    public:

      // Define input and output.
      EddyViscosityBaselineBoundary() {
        name_store("eddyViscosityLimit",eddyViscosityLimit) ;
        name_store("rho_f",rho_f) ;
        name_store("k_f",k_f) ;
        name_store("omega_f",omega_f) ;
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        name_store("eddyViscosity_f",eddyViscosity_f) ;
        input("eddyViscosityLimit,rho_f,k_f,omega_f,laminarViscosity_f") ;
        output("eddyViscosity_f") ;
        constraint("boundaryFaces,menterBSLTurbulenceModel") ;
      }

      // Compute the eddy viscosity for a single face.
      void calculate(Entity face) {
        eddyViscosity_f[face]=min(rho_f[face]*k_f[face]/omega_f[face],
          *eddyViscosityLimit*laminarViscosity_f[face]) ;
      }

      // Compute the eddy viscosity for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<EddyViscosityBaselineBoundary>
    registerEddyViscosityBaselineBoundary ;

  // Rule for the SST model to assign eddy viscosity for cells.
  class EddyViscositySSTInterior : public pointwise_rule {
    private:
      const_param<real> eddyViscosityLimit ;
      const_store<real> rho ;
      const_store<real> k ;
      const_store<real> omega ;
      const_store<tens3d> vGradient ;
      const_store<real> f2 ;
      const_store<real> laminarViscosity ;
      store<real> eddyViscosity ;
    public:

      // Define input and output.
      EddyViscositySSTInterior() {
        name_store("eddyViscosityLimit",eddyViscosityLimit) ;
        name_store("rho",rho) ;
        name_store("k",k) ;
        name_store("omega",omega) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("f2",f2) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("eddyViscosity",eddyViscosity) ;
        input("eddyViscosityLimit,rho,k,omega,gradv3d(v),f2") ;
        input("laminarViscosity") ;
        output("eddyViscosity") ;
        constraint("vol,menterSSTTurbulenceModel") ;
      }

      // Set eddy viscosity to zero for a single cell.
      void calculate(Entity cell) {
        real vorticityMagnitude=norm(vect3d(vGradient[cell].z.y-
          vGradient[cell].y.z,vGradient[cell].x.z-vGradient[cell].z.x,
          vGradient[cell].y.x-vGradient[cell].x.y)) ;
        real argA=0.31*omega[cell],argB=vorticityMagnitude*f2[cell] ;
        real kinematicEddyViscosity=(argA>argB)? k[cell]/omega[cell]:0.31*
          k[cell]/argB ;
        eddyViscosity[cell]=min(rho[cell]*kinematicEddyViscosity,
          *eddyViscosityLimit*laminarViscosity[cell]) ;
      }

      // Set eddy viscosity to zero for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<EddyViscositySSTInterior> registerEddyViscositySSTInterior ;

  // Rule for the SST model to assign eddy viscosity for boundary faces.
  class EddyViscositySSTBoundary : public pointwise_rule {
    private:
      const_Map ci ;
      const_param<real> eddyViscosityLimit ;
      const_store<tens3d> vGradient ;
      const_store<real> f2 ;
      const_store<real> rho_f ;
      const_store<real> k_f ;
      const_store<real> omega_f ;
      const_store<real> laminarViscosity_f ;
      store<real> eddyViscosity_f ;
    public:

      // Define input and output.
      EddyViscositySSTBoundary() {
        name_store("ci",ci) ;
        name_store("eddyViscosityLimit",eddyViscosityLimit) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("f2",f2) ;
        name_store("rho_f",rho_f) ;
        name_store("k_f",k_f) ;
        name_store("omega_f",omega_f) ;
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        name_store("eddyViscosity_f",eddyViscosity_f) ;
        input("eddyViscosityLimit") ;
        input("ci->gradv3d(v),ci->f2,rho_f,k_f,omega_f") ;
        input("laminarViscosity_f") ;
        output("eddyViscosity_f") ;
        constraint("boundaryFaces,menterSSTTurbulenceModel") ;
      }

      // Compute the eddy viscosity for a single face.
      void calculate(Entity face) {
        int e=ci[face] ;
        real vorticityMagnitude=norm(vect3d(vGradient[e].z.y-vGradient[e].y.z,
          vGradient[e].x.z-vGradient[e].z.x,vGradient[e].y.x-
          vGradient[e].x.y)) ;
        real argA=0.31*omega_f[face],argB=vorticityMagnitude*f2[e] ;
        real kinematicEddyViscosity=(argA>argB)? k_f[face]/omega_f[face]:0.31*
          k_f[face]/argB ;
        eddyViscosity_f[face]=min(rho_f[face]*kinematicEddyViscosity,
          *eddyViscosityLimit*laminarViscosity_f[face]) ;
      }

      // Compute the eddy viscosity for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<EddyViscositySSTBoundary> registerEddyViscositySSTBoundary ;
 
  // Rule for the 2003 SST model to assign eddy viscosity for cells.
  class EddyViscositySST2003Interior : public pointwise_rule {
    private:
      const_param<real> eddyViscosityLimit ;
      const_store<real> rho ;
      const_store<real> k ;
      const_store<real> omega ;
      const_store<tens3d> vGradient ;
      const_store<real> f2 ;
      const_store<real> laminarViscosity ;
      store<real> eddyViscosity ;
    public:

      // Define input and output.
      EddyViscositySST2003Interior() {
        name_store("eddyViscosityLimit",eddyViscosityLimit) ;
        name_store("rho",rho) ;
        name_store("k",k) ;
        name_store("omega",omega) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("f2",f2) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("eddyViscosity",eddyViscosity) ;
        input("eddyViscosityLimit,rho,k,omega,gradv3d(v),f2") ;
        input("laminarViscosity") ;
        output("eddyViscosity") ;
        constraint("vol,menterSST2003TurbulenceModel") ;
      }

      // Set eddy viscosity to zero for a single cell.
      void calculate(Entity cell) {
        tens3d strainRate=vGradient[cell]+Transpose(vGradient[cell]) ;
        real S=sqrt(ScalarProduct(strainRate,strainRate)/2.0) ;
        real argA=0.31*omega[cell],argB=S*f2[cell] ;
        real kinematicEddyViscosity=(argA>argB)? k[cell]/omega[cell]:0.31*
          k[cell]/argB ;
        eddyViscosity[cell]=min(rho[cell]*kinematicEddyViscosity,
          *eddyViscosityLimit*laminarViscosity[cell]) ;
      }

      // Set eddy viscosity to zero for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<EddyViscositySST2003Interior>
    registerEddyViscositySST2003Interior ;
 
  // Rule for the 2003 SST model to assign eddy viscosity for boundary faces.
  class EddyViscositySST2003Boundary : public pointwise_rule {
    private:
      const_Map ci ;
      const_param<real> eddyViscosityLimit ;
      const_store<tens3d> vGradient ;
      const_store<real> f2 ;
      const_store<real> rho_f ;
      const_store<real> k_f ;
      const_store<real> omega_f ;
      const_store<real> laminarViscosity_f ;
      store<real> eddyViscosity_f ;
    public:

      // Define input and output.
      EddyViscositySST2003Boundary() {
        name_store("ci",ci) ;
        name_store("eddyViscosityLimit",eddyViscosityLimit) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("f2",f2) ;
        name_store("rho_f",rho_f) ;
        name_store("k_f",k_f) ;
        name_store("omega_f",omega_f) ;
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        name_store("eddyViscosity_f",eddyViscosity_f) ;
        input("eddyViscosityLimit") ;
        input("ci->gradv3d(v),ci->f2,rho_f,k_f,omega_f") ;
        input("laminarViscosity_f") ;
        output("eddyViscosity_f") ;
        constraint("boundaryFaces,menterSST2003TurbulenceModel") ;
      }

      // Compute the eddy viscosity for a single face.
      void calculate(Entity face) {
        int e=ci[face] ;
        tens3d strainRate=vGradient[e]+Transpose(vGradient[e]) ;
        real S=sqrt(ScalarProduct(strainRate,strainRate)/2.0) ;
        real argA=0.31*omega_f[face],argB=S*f2[e] ;
        real kinematicEddyViscosity=(argA>argB)? k_f[face]/omega_f[face]:0.31*
          k_f[face]/argB ;
        eddyViscosity_f[face]=min(rho_f[face]*kinematicEddyViscosity,
          *eddyViscosityLimit*laminarViscosity_f[face]) ;
      }

      // Compute the eddy viscosity for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<EddyViscositySST2003Boundary>
    registerEddyViscositySST2003Boundary ;

  // Temporary hack to get around possible bug in using maps to set priority
  // variable values.
  class EddyViscosityCellHack : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> eddyViscosityWall ;
      store<real> eddyViscosityCell ;
    public:
                                                                                
      // Define input and output.
      EddyViscosityCellHack() {
        name_store("ci",ci) ;
        name_store("eddyViscosityWall",eddyViscosityWall) ;
        name_store("eddyViscosityCell",eddyViscosityCell) ;
        input("eddyViscosityWall") ;
        output("ci->eddyViscosityCell") ;
        constraint("ref->wallFunction_BCoption") ;
      }
                                                                                
      // Compute the eddy viscosity for a single face.
      void calculate(Entity face) {
        eddyViscosityCell[ci[face]]=eddyViscosityWall[face] ;
      }
                                                                                
      // Compute the eddy viscosity for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<EddyViscosityCellHack> registerEddyViscosityCellHack ;

  // Priority rule to assign eddy viscosity for cells next to wall-function
  // boundaries.
  class EddyViscosityWallFunction : public pointwise_rule {
    private:
      const_param<real> eddyViscosityLimit ;
      const_store<real> eddyViscosityCell ;
      store<real> eddyViscosity ;
    public:

      // Define input and output.
      EddyViscosityWallFunction() {
        name_store("eddyViscosityLimit",eddyViscosityLimit) ;
        name_store("eddyViscosityCell",eddyViscosityCell) ;
        name_store("wallFunction::eddyViscosity",eddyViscosity) ;
        input("eddyViscosityLimit,eddyViscosityCell") ;
        output("wallFunction::eddyViscosity") ;
        constraint("turbulentFlow,wallFunctionCells") ;
      }

      // Compute the eddy viscosity for a single face.
      void calculate(Entity cell) {
        eddyViscosity[cell]=min(eddyViscosityCell[cell],*eddyViscosityLimit) ;
      }

      // Compute the eddy viscosity for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<EddyViscosityWallFunction> registerEddyViscosityWallFunction ;

//-----------------------------------------------------------------------------
// Rules for the computing the laminar viscosity for post-processing analysis.
// Must rename in order to use output rules in scalarOutput.cc .

  // Laminar viscosity for cells.
  class LaminarViscosityInterior : public pointwise_rule {
    private:
      const_store<real> muu ;
      store<real> laminarViscosity ;
    public:

      // Define input and output.
      LaminarViscosityInterior() {
        name_store("muu(temperature,p,y)",muu) ;
        name_store("laminarViscosity",laminarViscosity) ;
        input("muu(temperature,p,y)") ;
        output("laminarViscosity") ;
        constraint("geom_cells") ;
      }

      // Compute laminar viscosity for a single cell.
      void calculate(Entity cell) { laminarViscosity[cell]=muu[cell] ; }

      // Compute laminar viscosity for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<LaminarViscosityInterior> registerLaminarViscosityInterior ;

  // Laminar viscosity for boundary faces.
  class LaminarViscosityBoundary : public pointwise_rule {
    private:
      const_store<real> muu_f ;
      store<real> laminarViscosity_f ;
    public:

      // Define input and output.
      LaminarViscosityBoundary() {
        name_store("muu(temperature_f,p_f,y_f)",muu_f) ;
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        input("muu(temperature_f,p_f,y_f)") ;
        output("laminarViscosity_f") ;
        constraint("boundaryFaces") ;
      }

      // Compute laminar viscosity for a single face.
      void calculate(Entity face) { laminarViscosity_f[face]=muu_f[face] ; }

      // Compute laminar viscosity for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<LaminarViscosityBoundary> registerLaminarViscosityBoundary ;

//-----------------------------------------------------------------------------
// Rules for the computing the ratio of eddy viscosity to laminar viscosity
// for post-processing analysis.

  // Viscosity ratio for cells.
  class ViscosityRatioInterior : public pointwise_rule {
    private:
      const_store<real> laminarViscosity ;
      const_store<real> eddyViscosity ;
      store<real> viscosityRatio ;
    public:

      // Define input and output.
      ViscosityRatioInterior() {
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("viscosityRatio",viscosityRatio) ;
        input("laminarViscosity,eddyViscosity") ;
        output("viscosityRatio") ;
        constraint("geom_cells") ;
      }

      // Compute viscosity ratio for a single cell.
      void calculate(Entity cell) {
        viscosityRatio[cell]=eddyViscosity[cell]/laminarViscosity[cell] ;
      }

      // Compute viscosity ratio for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ViscosityRatioInterior> registerViscosityRatioInterior ;

  // Viscosity ratio for boundary faces.
  class ViscosityRatioBoundary : public pointwise_rule {
    private:
      const_store<real> laminarViscosity_f ;
      const_store<real> eddyViscosity_f ;
      store<real> viscosityRatio_f ;
    public:

      // Define input and output.
      ViscosityRatioBoundary() {
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        name_store("eddyViscosity_f",eddyViscosity_f) ;
        name_store("viscosityRatio_f",viscosityRatio_f) ;
        input("laminarViscosity_f,eddyViscosity_f") ;
        output("viscosityRatio_f") ;
        constraint("boundaryFaces") ;
      }

      // Compute viscosity ratio for a single face.
      void calculate(Entity face) {
        viscosityRatio_f[face]=eddyViscosity_f[face]/laminarViscosity_f[face] ;
      }

      // Compute viscosity ratio for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ViscosityRatioBoundary> registerViscosityRatioBoundary ;

//-----------------------------------------------------------------------------
// Boundary condition rules for k.

  // Rule for boundary faces with specified k. Assigns value to all boundary
  // faces that have the property k_BC.
  class BoundaryKSpecification : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<real> k_BC ;
      store<real> k_f ;
    public:

      // Define input and output.
      BoundaryKSpecification() {
        name_store("ref",ref) ;
        name_store("k_BC",k_BC) ;
        name_store("k_f",k_f) ;
        input("ref->k_BC") ;
        output("k_f") ;
      }

      // Calculate k for a single face.
      void calculate(Entity face) { k_f[face]=k_BC[ref[face]] ; }

      // Calculate k for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryKSpecification>
    registerBoundaryKSpecification ;

  // Rule for specifying k with a profile. A single Cartesian coordinate is
  // used for the interpolation.
  class BoundaryKProfileCartesian : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<vect3d> faceCenter ;
      const_store<string> cartesianK_BC ;
      store<real> k_f ;
    public:

      // Define input and output.
      BoundaryKProfileCartesian() {
        name_store("ref",ref) ;
        name_store("facecenter",faceCenter) ;
        name_store("cartesianK_BC",cartesianK_BC) ;
        name_store("k_f",k_f) ;
        input("facecenter,ref->cartesianK_BC") ;
        output("k_f") ;
        disable_threading() ;
      }

      // Calculate k for all faces in sequence.
      virtual void compute(const sequence &seq) {

        if(Loci::GLOBAL_OR((seq != EMPTY))){

          // Create a map to organize ref values for the faces.
          std::map<int, Loci::entitySet> bcmap ;
          for(sequence::const_iterator si=seq.begin();si!=seq.end();++si){
            bcmap[ref[*si]]+=*si ;
          }

          // Loop through the map. Each map entry has differnet input file.
          std::map<int,Loci::entitySet>::iterator bci ;
          for(bci=bcmap.begin();bci!=bcmap.end();++bci) {

            // Open the file containing the profile.
            string fileName=cartesianK_BC[bci->first] ;
            ifstream in(fileName.c_str(),ios::in) ;
            if(in.fail()) {
              cerr << "Open failed on " << fileName.c_str() << endl ;
              Loci::Abort() ;
            }

            // Skip spaces.
            while(!in.eof() && isspace(in.peek())) in.get() ;

            // Read in the number of points on the boundary and the variable
            // flag which indicates the coordinate direction to use in the
            // interpolation.
            int np,coordFlag ; in >> np >> coordFlag ;
            if(np<2){
              cerr << "Bad number of data points in bc_k.dat." << endl ;
              Loci::Abort() ;
            }
            if(coordFlag<0 || coordFlag>2){
              cerr << "Bad coordinate flag in bc_k.dat." << endl ;
              Loci::Abort() ;
            }

            // Read in the coordinates and k values.
            vect3d *center=new vect3d[np] ; real *k =new real[np] ;
            for(int i=0;i<np;++i){
              in >> center[i] ; in >> k[i] ;
              if(k[i]<0.0){
                cerr << "Negative k for data point " << i << " in bc_k.dat."
                  << endl ; Loci::Abort() ;
              }
            }

            vect3d first = center[0] ; sequence s=sequence(bci->second) ;
            for(sequence::const_iterator si=s.begin();si!=s.end();++si) {
              int face=*si ;
              const vect3d &fcenter=faceCenter[face] ;
              int id1=0 ; real d1=0.0 ;
              switch(coordFlag){
                case 0: d1=fabs(first.x-fcenter.x) ; break ;
                case 1: d1=fabs(first.y-fcenter.y) ; break ;
                case 2: d1=fabs(first.z-fcenter.z) ; break ;
              }

              // Find the nearest point.
              for(int i=1;i<np;++i) {
                real dis=0.0 ;
                switch(coordFlag){
                  case 0: dis=fabs(fcenter.x-center[i].x) ; break ;
                  case 1: dis=fabs(fcenter.y-center[i].y) ; break ;
                  case 2: dis=fabs(fcenter.z-center[i].z) ; break ;
                }
                if(d1>dis){ d1=dis ; id1=i ; }
              }

              // Find the second nearest point
              int id2=0 ; real dd=0.0 ;
              switch(coordFlag){
                case 0: dd=center[id1].x-fcenter.x ; break ;
                case 1: dd=center[id1].y-fcenter.y ; break ;
                case 2: dd=center[id1].z-fcenter.z ; break ;
              }
              if(fabs(dd)>EPSILON)
                if(dd<0.0)
                  if(id1<np-1){
                    id2=id1+1 ;
                  }else{
                    k_f[face]=k[np-1] ; continue ;
                  }
                else
                  if(id1>0){
                    id2=id1-1 ;
                  }else{
                    k_f[face]=k[0] ; continue ;
                  }
              else {
                k_f[face]=k[id1] ; continue ;
              }
	      real d2=0.0 ;
              switch(coordFlag){
                case 0: d2=fabs(center[id2].x-fcenter.x) ; break ;
                case 1: d2=fabs(center[id2].y-fcenter.y) ; break ;
                case 2: d2=fabs(center[id2].z-fcenter.z) ; break ;
              }
              k_f[face]=(d2*k[id1]+d1*k[id2])/(d1+d2) ;
            }
            delete [] center ; delete [] k ;
          }
        }
      }
  } ;
                                                                                
  register_rule<BoundaryKProfileCartesian> registerBoundaryKProfileCartesian ;

  // Rule for specifying an axisymmetric k from a radial profile.
  class BoundaryKProfileAxisymmetric : public pointwise_rule {
    private:
      const_Map ref ;
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<string> axisymmetricK_BC ;
      const_store<unsigned int> referenceFrame_BC ;
      const_store<vect3d> faceCenter ;
      store<real> k_f ;
    public:

      // Define input and output.
      BoundaryKProfileAxisymmetric() {
        name_store("ref",ref) ;
        name_store("referenceFrame",referenceFrame) ;
        name_store("axisymmetricK_BC",axisymmetricK_BC) ;
        name_store("referenceFrame_BC",referenceFrame_BC) ;
        name_store("facecenter",faceCenter) ;
        name_store("k_f",k_f) ;
        input("referenceFrame,ref->(referenceFrame_BC,axisymmetricK_BC)") ;
        input("facecenter") ;
        output("k_f") ;
        disable_threading() ;
      }

      // Calculate k for all faces in sequence.
      virtual void compute(const sequence &seq) {

        if(Loci::GLOBAL_OR((seq != EMPTY))){

          // Create a map to organize ref values for the faces.
          std::map<int, Loci::entitySet> bcmap ;
          for(sequence::const_iterator si=seq.begin();si!=seq.end();++si){
            bcmap[ref[*si]]+=*si ;
          }
                                                                                
          // Loop through the map. Each map entry has differnet input file.
          std::map<int,Loci::entitySet>::iterator bci ;
          for(bci=bcmap.begin();bci!=bcmap.end();++bci) {
                                                                                
            // Open the file containing the profile.
            string fileName=axisymmetricK_BC[bci->first] ;
            ifstream in(fileName.c_str(),ios::in) ;
            if(in.fail()) {
              cerr << "Open failed on " << fileName.c_str() << endl ;
              Loci::Abort() ;
            }

            // Skip spaces.
            while(!in.eof() && isspace(in.peek())) in.get() ;

            // Read in the number of points in the profile.
            int np ; in >> np ;
            if(np<2){
              cerr << "Bad number of data points in bc_k.dat." << endl ;
              Loci::Abort() ;
            }

            // Read in the radius and k values.
            real *r=new real[np] ; real *k =new real[np] ;
            for(int i=0;i<np;++i){ in >> r[i] ; in >> k[i] ; }

            real first=r[0] ; sequence s=sequence(bci->second) ;
            for(sequence::const_iterator si=s.begin();si!=s.end();++si) {

              // Compute the axial, radial and theta unit vectors.
              int face=*si ; const vect3d &fcenter=faceCenter[face] ;
              vect3d delta=(*referenceFrame)[referenceFrame_BC[bci->first]].
                axisEnd-(*referenceFrame)[referenceFrame_BC[bci->first]].
                axisStart ;
              vect3d pos=fcenter-(*referenceFrame)[referenceFrame_BC
                [bci->first]].axisStart ;
              vect3d iHatZ=(1.0/norm(delta))*delta,iR=pos-dot(pos,iHatZ)*iHatZ ;

              // Find the nearest point.
              int id1=0 ; real d1=fabs(first-norm(iR)) ;
              for(int i=1;i<np;++i) {
                real dis=fabs(norm(iR)-r[i]) ; if(d1>dis){ d1=dis ; id1=i ; }
              }

              // Find the second nearest point
              int id2 ; real dd=r[id1]-norm(iR) ;
              if(fabs(dd)>EPSILON)
                if(dd<0.0)
                  if(id1<np-1){
                    id2=id1+1 ;
                  }else{
                    k_f[face]=k[np-1] ; continue ;
                  }
                else
                  if(id1>0){
                    id2=id1-1 ;
                  }else{
                    k_f[face]=k[0] ; continue ;
                  }
              else {
                k_f[face]=k[id1] ; continue ;
              }
              real d2=fabs(r[id2]-norm(iR)) ;
              k_f[face]=(d2*k[id1]+d1*k[id2])/(d1+d2) ;
            }
            delete [] r ; delete [] k ;
          }
        }
      }
  } ;

  register_rule<BoundaryKProfileAxisymmetric>
    registerBoundaryKProfileAxisymmetric ;

  // Rule for extrapolating k to boundary faces. This occurs for all outlets,
  // slip and symmetry boundaries. Right now we are using the low-order method.
  class BoundaryKExtrapolation : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> k ;
      store<real> k_f ;
    public:

      // Define input and output.
      BoundaryKExtrapolation() {
        name_store("ci",ci) ;
        name_store("k",k) ;
        name_store("k_f",k_f) ;
        input("ci->k") ;
        output("k_f") ;
        constraint("extrapolatedK_BC") ;
      }

      // Calculate k for a single face.
      void calculate(Entity face) { k_f[face]=k[ci[face]] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryKExtrapolation> registerBoundaryKExtrapolation ;

  // Rule for assigning k on no-slip boundary faces.
  class BoundaryKNoSlip : public pointwise_rule {
    private:
      store<real> k_f ;
    public:

      // Define input and output.
      BoundaryKNoSlip() {
        name_store("k_f",k_f) ;
        output("k_f") ;
        constraint("noslip_BC") ;
      }

      // Calculate k for a single face.
      void calculate(Entity face) { k_f[face]=0.0 ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryKNoSlip> registerBoundaryKNoSlip ;

//-----------------------------------------------------------------------------
// Boundary condition rules for omega.

  // Rule for boundary faces with specified omega. Assigns omega value to all
  // boundary faces that have the property omega_BC.
  class BoundaryOmegaSpecification : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<real> omega_BC ;
      store<real> omega_f ;
    public:

      // Define input and output.
      BoundaryOmegaSpecification() {
        name_store("ref",ref) ;
        name_store("omega_BC",omega_BC) ;
        name_store("omega_f",omega_f) ;
        input("ref->omega_BC") ;
        output("omega_f") ;
      }

      // Calculate omega for a single face.
      void calculate(Entity face) { omega_f[face]=omega_BC[ref[face]] ; }

      // Calculate omega for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryOmegaSpecification>
    registerBoundaryOmegaSpecification ;

  // Rule for specifying omega with a profile. A single Cartesian coordinate is
  // used for the interpolation.
  class BoundaryOmegaProfileCartesian : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<vect3d> faceCenter ;
      const_store<string> cartesianOmega_BC ;
      store<real> omega_f ;
    public:

      // Define input and output.
      BoundaryOmegaProfileCartesian() {
        name_store("ref",ref) ;
        name_store("facecenter",faceCenter) ;
        name_store("cartesianOmega_BC",cartesianOmega_BC) ;
        name_store("omega_f",omega_f) ;
        input("facecenter,ref->cartesianOmega_BC") ;
        output("omega_f") ;
        disable_threading() ;
      }

      // Calculate omega for all faces in sequence.
      virtual void compute(const sequence &seq) {

        if(Loci::GLOBAL_OR((seq != EMPTY))){

          // Create a map to organize ref values for the faces.
          std::map<int, Loci::entitySet> bcmap ;
          for(sequence::const_iterator si=seq.begin();si!=seq.end();++si){
            bcmap[ref[*si]]+=*si ;
          }

          // Loop through the map. Each map entry has differnet input file.
          std::map<int,Loci::entitySet>::iterator bci ;
          for(bci=bcmap.begin();bci!=bcmap.end();++bci) {

            // Open the file containing the profile.
            string fileName=cartesianOmega_BC[bci->first] ;
            ifstream in(fileName.c_str(),ios::in) ;
            if(in.fail()) {
              cerr << "Open failed on " << fileName.c_str() << endl ;
              Loci::Abort() ;
            }

            // Skip spaces.
            while(!in.eof() && isspace(in.peek())) in.get() ;

            // Read in the number of points on the boundary and the variable
            // flag which indicates the coordinate direction to use in the
            // interpolation.
            int np,coordFlag ; in >> np >> coordFlag ;
            if(np<2){
              cerr << "Bad number of data points in bc_omega.dat." << endl ;
              Loci::Abort() ;
            }
            if(coordFlag<0 || coordFlag>2){
              cerr << "Bad coordinate flag in bc_omega.dat." << endl ;
              Loci::Abort() ;
            }

            // Read in the coordinates and k values.
            vect3d *center=new vect3d[np] ; real *omega =new real[np] ;
            for(int i=0;i<np;++i){
              in >> center[i] ; in >> omega[i] ;
              if(omega[i]<0.0){
                cerr << "Negative omega for data point " << i
                  << " in bc_omega.dat." << endl ; Loci::Abort() ;
              }
            }

            vect3d first = center[0] ; sequence s=sequence(bci->second) ;
            for(sequence::const_iterator si=s.begin();si!=s.end();++si) {
              int face=*si ;
              const vect3d &fcenter=faceCenter[face] ;
              int id1=0 ; real d1=0.0 ;
              switch(coordFlag){
                case 0: d1=fabs(first.x-fcenter.x) ; break ;
                case 1: d1=fabs(first.y-fcenter.y) ; break ;
                case 2: d1=fabs(first.z-fcenter.z) ; break ;
              }

              // Find the nearest point.
              for(int i=1;i<np;++i) {
                real dis=0.0 ;
                switch(coordFlag){
                  case 0: dis=fabs(fcenter.x-center[i].x) ; break ;
                  case 1: dis=fabs(fcenter.y-center[i].y) ; break ;
                  case 2: dis=fabs(fcenter.z-center[i].z) ; break ;
                }
                if(d1>dis){ d1=dis ; id1=i ; }
              }

              // Find the second nearest point
              int id2=0 ; real dd=0.0 ;
              switch(coordFlag){
                case 0: dd=center[id1].x-fcenter.x ; break ;
                case 1: dd=center[id1].y-fcenter.y ; break ;
                case 2: dd=center[id1].z-fcenter.z ; break ;
              }
              if(fabs(dd)>EPSILON)
                if(dd<0.0)
                  if(id1<np-1){
                    id2=id1+1 ;
                  }else{
                    omega_f[face]=omega[np-1] ; continue ;
                  }
                else
                  if(id1>0){
                    id2=id1-1 ;
                  }else{
                    omega_f[face]=omega[0] ; continue ;
                  }
              else {
                omega_f[face]=omega[id1] ; continue ;
              }
              real d2=0.0 ;
              switch(coordFlag){
                case 0: d2=fabs(center[id2].x-fcenter.x) ; break ;
                case 1: d2=fabs(center[id2].y-fcenter.y) ; break ;
                case 2: d2=fabs(center[id2].z-fcenter.z) ; break ;
              }
              omega_f[face]=(d2*omega[id1]+d1*omega[id2])/(d1+d2) ;
            }
            delete [] center ; delete [] omega ;
          }
        }
      }
  } ;
                                                                                
  register_rule<BoundaryOmegaProfileCartesian>
    registerBoundaryOmegaProfileCartesian ;

  // Rule for specifying an axisymmetric omega from a radial profile.
  class BoundaryOmegaProfileAxisymmetric : public pointwise_rule {
    private:
      const_Map ref ;
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<string> axisymmetricOmega_BC ;
      const_store<unsigned int> referenceFrame_BC ;
      const_store<vect3d> faceCenter ;
      store<real> omega_f ;
    public:

      // Define input and output.
      BoundaryOmegaProfileAxisymmetric() {
        name_store("ref",ref) ;
        name_store("referenceFrame",referenceFrame) ;
        name_store("axisymmetricOmega_BC",axisymmetricOmega_BC) ;
        name_store("referenceFrame_BC",referenceFrame_BC) ;
        name_store("facecenter",faceCenter) ;
        name_store("omega_f",omega_f) ;
        input("referenceFrame,ref->(axisymmetricOmega_BC,referenceFrame_BC)") ;
        input("facecenter") ;
        output("omega_f") ;
        disable_threading() ;
      }

      // Calculate omega for all faces in sequence.
      virtual void compute(const sequence &seq) {

        if(Loci::GLOBAL_OR((seq != EMPTY))){
                                                                                
          // Create a map to organize ref values for the faces.
          std::map<int, Loci::entitySet> bcmap ;
          for(sequence::const_iterator si=seq.begin();si!=seq.end();++si){
            bcmap[ref[*si]]+=*si ;
          }
                                                                                
          // Loop through the map. Each map entry has differnet input file.
          std::map<int,Loci::entitySet>::iterator bci ;
          for(bci=bcmap.begin();bci!=bcmap.end();++bci) {
                                                                                
            // Open the file containing the profile.
            string fileName=axisymmetricOmega_BC[bci->first] ;
            ifstream in(fileName.c_str(),ios::in) ;
            if(in.fail()) {
              cerr << "Open failed on " << fileName.c_str() << endl ;
              Loci::Abort() ;
            }

            // Skip spaces.
            while(!in.eof() && isspace(in.peek())) in.get() ;

            // Read in the number of points in the profile.
            int np ; in >> np ;
            if(np<2){
              cerr << "Bad number of data points in bc_omega.dat." << endl ;
              Loci::Abort() ;
            }

            // Read in the radius and k values.
            real *r=new real[np] ; real *omega =new real[np] ;
            for(int i=0;i<np;++i){ in >> r[i] ; in >> omega[i] ; }

            real first=r[0] ; sequence s=sequence(bci->second) ;
            for(sequence::const_iterator si=s.begin();si!=s.end();++si) {

              // Compute the axial, radial and theta unit vectors.
              int face=*si ; const vect3d &fcenter=faceCenter[face] ;
              vect3d delta=(*referenceFrame)[referenceFrame_BC[bci->first]].
                axisEnd-(*referenceFrame)[referenceFrame_BC[bci->first]].
                axisStart ;
              vect3d pos=fcenter-(*referenceFrame)[referenceFrame_BC
                [bci->first]].axisStart ;
              vect3d iHatZ=(1.0/norm(delta))*delta,iR=pos-dot(pos,iHatZ)*iHatZ ;

              // Find the nearest point.
              int id1=0 ; real d1=fabs(first-norm(iR)) ;
              for(int i=1;i<np;++i) {
                real dis=fabs(norm(iR)-r[i]) ; if(d1>dis){ d1=dis ; id1=i ; }
              }

              // Find the second nearest point
              int id2 ; real dd=r[id1]-norm(iR) ;
              if(fabs(dd)>EPSILON)
                if(dd<0.0)
                  if(id1<np-1){
                    id2=id1+1 ;
                  }else{
                    omega_f[face]=omega[np-1] ; continue ;
                  }
                else
                  if(id1>0){
                    id2=id1-1 ;
                  }else{
                    omega_f[face]=omega[0] ; continue ;
                  }
              else {
                omega_f[face]=omega[id1] ; continue ;
              }
              real d2=fabs(r[id2]-norm(iR)) ;
              omega_f[face]=(d2*omega[id1]+d1*omega[id2])/(d1+d2) ;
            }
            delete [] r ; delete [] omega ;
          }
        }
      }
  } ;

  register_rule<BoundaryOmegaProfileAxisymmetric>
    registerBoundaryOmegaProfileAxisymmetric ;

  // Rule for extrapolating omega to boundary faces. This occurs for all
  // outlets, slip and symmetry boundaries. Right now we are using the
  // low-order method.
  class BoundaryOmegaExtrapolation : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> omega ;
      store<real> omega_f ;
    public:

      // Define input and output.
      BoundaryOmegaExtrapolation() {
        name_store("ci",ci) ;
        name_store("omega",omega) ;
        name_store("omega_f",omega_f) ;
        input("ci->omega") ;
        output("omega_f") ;
        constraint("extrapolatedOmega_BC") ;
      }

      // Calculate omega for a single face.
      void calculate(Entity face) { omega_f[face]=omega[ci[face]] ; }

      // Calculate omega for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryOmegaExtrapolation> registerBoundaryOmegaExtrapolation ;

  // Rule for assigning omega on no-slip boundary faces for Menter's family of
  // k-omega turbulence models.
  class BoundaryOmegaNoSlipMenter : public pointwise_rule {
    private:
      const_Map ci ;
      const_param<real> beta1 ;
      const_store<real> distanceToNoSlip ;
      const_store<real> rho_f ;
      const_store<real> laminarViscosity_f ;
      store<real> omega_f ;
    public:

      // Define input and output.
      BoundaryOmegaNoSlipMenter() {
        name_store("ci",ci) ;
        name_store("beta1",beta1) ;
        name_store("distanceToNoSlip",distanceToNoSlip) ;
        name_store("rho_f",rho_f) ;
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        name_store("omega_f",omega_f) ;
        input("beta1,ci->distanceToNoSlip,rho_f,laminarViscosity_f") ;
        output("omega_f") ;
        constraint("noslip_BC,menterTurbulenceModel") ;
      }

      // Calculate omega for a single face.
      void calculate(Entity face) {
        omega_f[face]=60.0*laminarViscosity_f[face]/(*beta1*distanceToNoSlip
          [ci[face]]*distanceToNoSlip[ci[face]]*rho_f[face]) ;
      }

      // Calculate omega for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryOmegaNoSlipMenter> registerBoundaryOmegaNoSlipMenter ;

//-----------------------------------------------------------------------------
// Rules to create a constraint for boundary faces with non-zero diffusion
// flux for turbulence quantities.

  // All inlet faces have non-zero diffusion flux.
  class BoundaryTurbulenceDiffusionInlet : public pointwise_rule {
    private:
      store<bool> boundaryTurbulenceDiffusion ;
    public:

      // Define input and output.
      BoundaryTurbulenceDiffusionInlet() {
        name_store("boundaryTurbulenceDiffusion",boundaryTurbulenceDiffusion) ;
        output("boundaryTurbulenceDiffusion") ;
        constraint("inlet_BC") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { boundaryTurbulenceDiffusion[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryTurbulenceDiffusionInlet>
    registerBoundaryTurbulenceDiffusionInlet ;

  // All outlet faces have non-zero diffusion flux.
  class BoundaryTurbulenceDiffusionOutlet : public pointwise_rule {
    private:
      store<bool> boundaryTurbulenceDiffusion ;
    public:

      // Define input and output.
      BoundaryTurbulenceDiffusionOutlet() {
        name_store("boundaryTurbulenceDiffusion",boundaryTurbulenceDiffusion) ;
        output("boundaryTurbulenceDiffusion") ;
        constraint("outlet_BC") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { boundaryTurbulenceDiffusion[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryTurbulenceDiffusionOutlet>
    registerBoundaryTurbulenceDiffusionOutlet ;

  // No-slip faces have non-zero diffusion flux.
  class BoundaryTurbulenceDiffusionNoSlip : public pointwise_rule {
    private:
      store<bool> boundaryTurbulenceDiffusion ;
    public:

      // Define input and output.
      BoundaryTurbulenceDiffusionNoSlip() {
        name_store("boundaryTurbulenceDiffusion",boundaryTurbulenceDiffusion) ;
        output("boundaryTurbulenceDiffusion") ;
        constraint("noslip_BC") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { boundaryTurbulenceDiffusion[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryTurbulenceDiffusionNoSlip>
    registerBoundaryTurbulenceDiffusionNoSlip ;

  // Interface faces have non-zero diffusion flux.
  class BoundaryTurbulenceDiffusionInterface : public pointwise_rule {
    private:
      store<bool> boundaryTurbulenceDiffusion ;
    public:

      // Define input and output.
      BoundaryTurbulenceDiffusionInterface() {
        name_store("boundaryTurbulenceDiffusion",boundaryTurbulenceDiffusion) ;
        output("boundaryTurbulenceDiffusion") ;
        constraint("interface_BC") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { boundaryTurbulenceDiffusion[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryTurbulenceDiffusionInterface>
    registerBoundaryTurbulenceDiffusionInterface ;

//-----------------------------------------------------------------------------
// Scheme independent rules for assembling the k equation.

  // Rule to initialize the main coefficient. This is now a parametric variable,
  // with X=massFlux and Y=density.
  class InitializeKMainCoefficient : public unit_rule {
    private:
      store<real> kMainCoefficient ;
    public:

      // Define input and output.
      InitializeKMainCoefficient() {
        name_store("kMainCoefficient(X,Y)",kMainCoefficient) ;
        output("kMainCoefficient(X,Y)") ;
        constraint("vol") ;
      }

      // Set the main coefficient to zero for a single cell.
      void calculate(Entity cell) { kMainCoefficient[cell]=0.0 ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeKMainCoefficient> registerInitializeKMainCoefficient ;

  // Rule to add the first-order inviscid flux contribution to the main
  // coefficient for interior faces.
  class FOUInviscidFluxToKMainCoefficientInterior : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      store<real> kMainCoefficient ;
    public:

      // Define input and output.
      FOUInviscidFluxToKMainCoefficientInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("X",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("kMainCoefficient(X,Y)",kMainCoefficient) ;
        input("X,gridMassFlux") ;
        output("(cl,cr)->kMainCoefficient(X,Y)") ;
        constraint("internalFaces") ;
      }

      // Increment the main coefficient for cells attached to a single face.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        if(netMassFlux>0.0){
          kMainCoefficient[cl[face]]+=netMassFlux ;
        }else{
          kMainCoefficient[cr[face]]-=netMassFlux ;
        }
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToKMainCoefficientInterior>
    registerFOUInviscidFluxToKMainCoefficientInterior ;

  // Rule to add the first-order inviscid flux contribution to the main
  // coefficient for boundary faces.
  class FOUInviscidFluxToKMainCoefficientBoundary : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      store<real> kMainCoefficient ;
    public:

      // Define input and output.
      FOUInviscidFluxToKMainCoefficientBoundary() {
        name_store("ci",ci) ;
        name_store("X",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("kMainCoefficient(X,Y)",kMainCoefficient) ;
        input("X,gridMassFlux") ;
        output("ci->kMainCoefficient(X,Y)") ;
        constraint("boundaryFaces") ;
      }

      // Increment the main coefficient for the cell attach to a single face.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        if(netMassFlux>0.0) kMainCoefficient[ci[face]]+=netMassFlux ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToKMainCoefficientBoundary>
    registerFOUInviscidFluxToKMainCoefficientBoundary ;

  // Rule to add the destruction term to the main coefficient. Note that
  // we have now added a multiplier here which allows for the implementation
  // of DES models. 
  class DestructionToKMainCoefficient : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> betaStar ;
      const_store<real> rho ;
      const_store<real> omega ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      const_store<real> fDES ;
      store<real> kMainCoefficient ;
    public:

      // Define input and output.
      DestructionToKMainCoefficient() {
        name_store("betaStar",betaStar) ;
        name_store("rho",rho) ;
        name_store("omega",omega) ;
        name_store("vol",vol) ;
        name_store("cellRadius",cellRadius) ;
        name_store("fDES",fDES) ;
        name_store("kMainCoefficient(X,Y)",kMainCoefficient) ;
        input("betaStar,rho,omega,vol,cellRadius,fDES") ;
        output("kMainCoefficient(X,Y)") ;
        constraint("geom_cells,menterTurbulenceModel") ;
      }

      // Increment the main coefficient for a single cell.
      void calculate(Entity cell) {
        kMainCoefficient[cell]+=(*betaStar)*rho[cell]*omega[cell]*vol[cell]*
          cellRadius[cell]*fDES[cell] ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DestructionToKMainCoefficient>
    registerDestructionToKMainCoefficient ;

  // Rule to add the diffusive flux contribution to the main coefficient for
  // interior faces.
  class DiffusiveFluxToKMainCoefficientInterior : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> faceRadius ;
      const_store<real> diffusionProduct ;
      const_store<real> kViscosity ;
      store<real> kMainCoefficient ;
    public:

      // Define input and output.
      DiffusiveFluxToKMainCoefficientInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("faceRadius",faceRadius) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("kViscosity",kViscosity) ;
        name_store("kMainCoefficient(X,Y)",kMainCoefficient) ;
        input("diffusionProduct,faceRadius,kViscosity") ;
        output("(cl,cr)->kMainCoefficient(X,Y)") ;
        constraint("internalFaces") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real temp=kViscosity[face]*diffusionProduct[face]*faceRadius[face] ;
        kMainCoefficient[cl[face]]+=temp ; kMainCoefficient[cr[face]]+=temp ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToKMainCoefficientInterior>
    registerDiffusiveFluxToKMainCoefficientInterior ;

  // Rule to add the diffusive flux contribution to the main coefficient for
  // boundary faces.
  class DiffusiveFluxToKMainCoefficientBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_param<real> sigmaK1,sigmaK2 ;
      const_store<real> faceRadius ;
      const_store<real> diffusionProduct ;
      const_store<real> laminarViscosity_f ;
      const_store<real> eddyViscosity_f ;
      const_store<real> f1 ;
      store<real> kMainCoefficient ;
    public:

      // Define input and output.
      DiffusiveFluxToKMainCoefficientBoundary() {
        name_store("ci",ci) ;
        name_store("sigmaK1",sigmaK1) ;
        name_store("sigmaK2",sigmaK2) ;
        name_store("faceRadius",faceRadius) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        name_store("eddyViscosity_f",eddyViscosity_f) ;
        name_store("f1",f1) ;
        name_store("kMainCoefficient(X,Y)",kMainCoefficient) ;
        input("sigmaK1,sigmaK2,ci->f1,faceRadius") ;
        input("diffusionProduct,laminarViscosity_f,eddyViscosity_f") ;
        output("ci->kMainCoefficient(X,Y)") ;
        constraint("boundaryTurbulenceDiffusion") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real sigmaK=f1[ci[face]]*(*sigmaK1)+(1.0-f1[ci[face]])*(*sigmaK2) ;
        kMainCoefficient[ci[face]]+=(laminarViscosity_f[face]+sigmaK*
          eddyViscosity_f[face])*diffusionProduct[face]*faceRadius[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToKMainCoefficientBoundary>
    registerDiffusiveFluxToKMainCoefficientBoundary ;

  // Rule to add temporal component of momentum equation to the main
  // coefficient. This is now an iteration independent rule so that it can
  // be used at both {n} and {n,it} where it is needed.
  class TemporalToKMainCoefficient: public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> timeStep ;
      const_param<real> timeIntegratorFactor0 ;
      const_store<real> timeStepFactor ;
      const_store<real> rho ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> kMainCoefficient ;
    public:

      // Define input and output.
      TemporalToKMainCoefficient() {
        name_store("timeStep",timeStep) ;
        name_store("timeIntegratorFactor0",timeIntegratorFactor0) ;
        name_store("timeStepFactor",timeStepFactor) ;
        name_store("Y",rho) ;
        name_store("vol",vol) ;
        name_store("cellRadius",cellRadius) ;
        name_store("kMainCoefficient(X,Y)",kMainCoefficient) ;
        input("timeStep,timeStepFactor,timeIntegratorFactor0") ;
        input("Y,vol,cellRadius") ;
        output("kMainCoefficient(X,Y)") ;
        constraint("geom_cells") ;
      }

      // Add temporal component for a single cell.
      void calculate(Entity cell) {
        kMainCoefficient[cell]+=0.5*rho[cell]*vol[cell]*cellRadius[cell]*
          ((*timeIntegratorFactor0)+1.0)/((*timeStep)*timeStepFactor[cell]) ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToKMainCoefficient> registerTemporalToKMainCoefficient ;

  // Rule to initialize the source term. This is now a parametric variable
  // with X=velocity and Y=massFlux.
  class InitializeKSourceTerm : public unit_rule {
    private:
      store<real> kSourceTerm ;
    public:

      // Define input and output.
      InitializeKSourceTerm() {
        name_store("kSourceTerm(X,Y)",kSourceTerm) ;
        output("kSourceTerm(X,Y)") ;
        constraint("vol") ;
      }

      // Set the source term to zero for a single cell.
      void calculate(Entity cell) { kSourceTerm[cell]=0.0 ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeKSourceTerm> registerInitializeKSourceTerm ;

  // Rule to add the first-order inviscid flux contribution to the source term
  // for boundary faces.
  class FOUInviscidFluxToKSourceTermBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      const_store<real> k_f ;
      store<real> kSourceTerm ;
    public:

      // Define input and output.
      FOUInviscidFluxToKSourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("Y",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("k_f",k_f) ;
        name_store("kSourceTerm(X,Y)",kSourceTerm) ;
        input("Y,gridMassFlux,k_f") ;
        output("ci->kSourceTerm(X,Y)") ;
        constraint("boundaryFaces") ;
      }

      // Increment the source term for the cell attach to a single face.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        if(netMassFlux<0.0) kSourceTerm[ci[face]]-=netMassFlux*k_f[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToKSourceTermBoundary>
    registerFOUInviscidFluxToKSourceTermBoundary ;

  // Rule to add the second-order convection contribution to the source term
  // for interior faces.
  class SOUInviscidFluxToKSourceTermInterior : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> k ;
      const_store<vect3d> kGradient ;
      const_store<real> kLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      store<real> kSourceTerm ;
    public:

      // Define input and output.
      SOUInviscidFluxToKSourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("k",k) ;
        name_store("grads(k)",kGradient) ;
        name_store("limiters(k)",kLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("Y",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("kSourceTerm(X,Y)",kSourceTerm) ;
        input("(cl,cr)->(k,grads(k),limiters(k),cellcenter)") ;
        input("facecenter,Y,gridMassFlux") ;
        output("(cl,cr)->kSourceTerm(X,Y)") ;
        constraint("internalFaces,souTurbulence") ;
      }

      // Increment the source term for the cells attach to a single face. Note
      // that here we are insisting that the reconstructed value of k at the
      // face be positive.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        real secondOrderSource=(netMassFlux>0.0)? netMassFlux*
          max(kLimiter[cl[face]]*dot(kGradient[cl[face]],faceCenter[face]-
          cellCenter[cl[face]]),-k[cl[face]]):netMassFlux*
          max(kLimiter[cr[face]]*dot(kGradient[cr[face]],faceCenter[face]-
          cellCenter[cr[face]]),-k[cr[face]]) ;
        kSourceTerm[cl[face]]-=secondOrderSource ;
        kSourceTerm[cr[face]]+=secondOrderSource ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SOUInviscidFluxToKSourceTermInterior>
    registerSOUInviscidFluxToKSourceTermInterior ;

  // Rule to add the second-order convection contribution to the source term
  // for boundary faces.
  class SOUInviscidFluxToKSourceTermBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> k ;
      const_store<vect3d> kGradient ;
      const_store<real> kLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      const_store<real> k_f ;
      store<real> kSourceTerm ;
    public:

      // Define input and output.
      SOUInviscidFluxToKSourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("k",k) ;
        name_store("grads(k)",kGradient) ;
        name_store("limiters(k)",kLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("Y",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("k_f",k_f) ;
        name_store("kSourceTerm(X,Y)",kSourceTerm) ;
        input("ci->(k,grads(k),limiters(k),cellcenter)") ;
        input("facecenter,Y,gridMassFlux,k_f") ;
        output("ci->kSourceTerm(X,Y)") ;
        constraint("boundaryFaces,souTurbulence") ;
      }

      // Increment the source term for the cells attach to a single face. Note
      // that here we are insisting that the reconstructed value of k at the
      // face be positive.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        if(netMassFlux>0.0) kSourceTerm[ci[face]]-=netMassFlux*
          max(kLimiter[ci[face]]*dot(kGradient[ci[face]],faceCenter[face]-
          cellCenter[ci[face]]),-k[ci[face]]) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SOUInviscidFluxToKSourceTermBoundary>
    registerSOUInviscidFluxToKSourceTermBoundary ;

  // Rule to add the production term to the source term. Checked.
  class ProductionToKSourceTerm : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> cellRadius ;
      const_store<real> kProduction ;
      store<real> kSourceTerm ;
    public:

      // Define input and output.
      ProductionToKSourceTerm() {
        name_store("cellRadius",cellRadius) ;
        name_store("kProduction(X)",kProduction) ;
        name_store("kSourceTerm(X,Y)",kSourceTerm) ;
        input("cellRadius,kProduction(X)") ;
        output("kSourceTerm(X,Y)") ;
        constraint("geom_cells") ;
      }

      // Increment the source term for a single cell.
      void calculate(Entity cell) {
        kSourceTerm[cell]+=cellRadius[cell]*kProduction[cell] ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ProductionToKSourceTerm> registerProductionToKSourceTerm ;

  // Rule to add the diffusive flux contribution to the source term for
  // interior faces.
  class DiffusiveFluxToKSourceTermInterior : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> vol ;
      const_store<vect3d> kGradient ;
      const_store<vect3d> geometryFactor0 ;
      const_store<real> faceRadius ;
      const_store<real> kViscosity ;
      store<real> kSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToKSourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("vol",vol) ;
        name_store("grads(k)",kGradient) ;
        name_store("geometryFactor0",geometryFactor0) ;
        name_store("faceRadius",faceRadius) ;
        name_store("kViscosity",kViscosity) ;
        name_store("kSourceTerm(X,Y)",kSourceTerm) ;
        input("(cl,cr)->(vol,grads(k))") ;
        input("geometryFactor0,faceRadius,kViscosity") ;
        output("(cl,cr)->kSourceTerm(X,Y)") ;
        constraint("internalFaces") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real secondarySourceTerm=kViscosity[face]*dot((kGradient[cl[face]]*
          vol[cr[face]]+kGradient[cr[face]]*vol[cl[face]]),
          geometryFactor0[face])*faceRadius[face]/(vol[cl[face]]+vol[cr[face]]);
        kSourceTerm[cl[face]]+=secondarySourceTerm ;
        kSourceTerm[cr[face]]-=secondarySourceTerm ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToKSourceTermInterior>
    registerDiffusiveFluxToKSourceTermInterior ;

  // Rule to add the diffusive flux contribution to the source term for
  // boundary faces.
  class DiffusiveFluxToKSourceTermBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_param<real> sigmaK1,sigmaK2 ;
      const_store<real> f1 ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> kGradient ;
      const_store<vect3d> faceCenter ;
      const_store<real> diffusionProduct ;
      const_store<real> k_f ;
      const_store<real> laminarViscosity_f ;
      const_store<real> eddyViscosity_f ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<real> kSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToKSourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("sigmaK1",sigmaK1) ;
        name_store("sigmaK2",sigmaK2) ;
        name_store("f1",f1) ;
        name_store("cellcenter",cellCenter) ;
        name_store("grads(k)",kGradient) ;
        name_store("facecenter",faceCenter) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("k_f",k_f) ;
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        name_store("eddyViscosity_f",eddyViscosity_f) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("kSourceTerm(X,Y)",kSourceTerm) ;
        input("sigmaK1,sigmaK2,ci->f1,ci->cellcenter,ci->grads(k)") ;
        input("facecenter,diffusionProduct,k_f") ;
        input("laminarViscosity_f,eddyViscosity_f,area,faceRadius") ;
        output("ci->kSourceTerm(X,Y)") ;
        constraint("boundaryTurbulenceDiffusion") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real sigmaK=f1[ci[face]]*(*sigmaK1)+(1.0-f1[ci[face]])*(*sigmaK2) ;
        real sourceTerm=(laminarViscosity_f[face]+sigmaK*eddyViscosity_f[face])*
          (k_f[face]*diffusionProduct[face]+dot(kGradient[ci[face]],
          (area[face].n*area[face].sada-diffusionProduct[face]*
          (faceCenter[face]-cellCenter[ci[face]]))))*faceRadius[face] ;
        kSourceTerm[ci[face]]+=sourceTerm ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToKSourceTermBoundary>
    registerDiffusiveFluxToKSourceTermBoundary ;

  // Rule to compute the temporal source term for BDF. Note that this
  // term only includes the 'old' components from the {n} level. This
  // rule will compute and be used at {n}. In addition, the value will
  // promote to {n,it} where it is also needed.
  class KTemporalSourceTermBDF : public pointwise_rule {
    private:
      const_param<real> timeStep ;
      const_store<real> timeStepFactor ;
      const_store<real> rho,vol,cellRadius ;
      const_store<real> k ;
      store<real> kTemporalSourceTerm ;
    public:

      // Define input and output.
      KTemporalSourceTermBDF() {
        name_store("timeStep{n}",timeStep) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n}",rho) ;
        name_store("k{n}",k) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n}",cellRadius) ;
        name_store("kTemporalSourceTerm{n}",kTemporalSourceTerm) ;
        input("rho{n},k{n},vol{n},cellRadius{n}") ;
        input("timeStep{n},timeStepFactor{n}") ;
        output("kTemporalSourceTerm{n}") ;
        constraint("BDFIntegrator,geom_cells") ;
      }

      // Compute for a single cell.
      void calculate(Entity cell) {
        kTemporalSourceTerm[cell]=(rho[cell]*vol[cell]*cellRadius[cell]/
          ((*timeStep)*timeStepFactor[cell]))*k[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<KTemporalSourceTermBDF> registerKTemporalSourceTermBDF ;

  // Rule to compute the temporal source term for BDF2. Note that this
  // term only includes the 'old' components from the {n-1} and {n} levels.
  // This rule will compute and be used at {n}. In addition, the value will
  // promote to {n,it} where it is also needed.
  class KTemporalSourceTermBDF2 : public pointwise_rule {
    private:
      const_param<int> n ;
      const_param<real> timeStep ;
      const_store<real> timeStepFactor ;
      const_store<real> rhoOld,rho ;
      const_store<real> kOld,k ;
      const_store<real> volOld,vol ;
      const_store<real> cellRadius ;
      store<real> kTemporalSourceTerm ;
    public:

      // Define input and output.
      KTemporalSourceTermBDF2() {
        name_store("$n{n}",n) ;
        name_store("timeStep{n}",timeStep) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n-1}",rhoOld) ;
        name_store("rho{n}",rho) ;
        name_store("k{n-1}",kOld) ;
        name_store("k{n}",k) ;
//      name_store("vol{n-1}",volOld) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n}",cellRadius) ;
        name_store("kTemporalSourceTerm{n}",kTemporalSourceTerm) ;
//      input("$n{n},rho{n-1},rho{n},k{n-1},k{n},vol{n-1},vol{n}") ;
        input("$n{n},rho{n-1},rho{n},k{n-1},k{n},vol{n}") ;
        input("cellRadius{n},timeStep{n},timeStepFactor{n}") ;
        output("kTemporalSourceTerm{n}") ;
        constraint("BDF2Integrator,geom_cells") ;
      }

      // Compute for a single cell. NOTE: Not using old volume until issue
      // with vol{n=1} for static meshes is resolved.
      void calculate(Entity cell) {
        if((*n)!=0){
          kTemporalSourceTerm[cell]=((2.0*rho[cell]*vol[cell])*k[cell]-(0.5*
            rhoOld[cell]*vol[cell])*kOld[cell])*cellRadius[cell]/
           ((*timeStep)*timeStepFactor[cell]) ;
        }else{
          kTemporalSourceTerm[cell]=(rho[cell]*vol[cell]*cellRadius[cell]/
            ((*timeStep)*timeStepFactor[cell]))*k[cell] ;
        }
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<KTemporalSourceTermBDF2>
    registerKTemporalSourceTermBDF2 ;

  // Rule to add temporal component of momentum equation to the source term.
  // This is now an iteration independent rule so that it can be used at both
  // {n} and {n,it} where it is needed.
  class TemporalToKSourceTerm : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> kTemporalSourceTerm ;
      store<real> kSourceTerm ;
    public:

      // Define input and output.
      TemporalToKSourceTerm() {
        name_store("kTemporalSourceTerm",kTemporalSourceTerm) ;
        name_store("kSourceTerm(X,Y)",kSourceTerm) ;
        input("kTemporalSourceTerm") ;
        output("kSourceTerm(X,Y)") ;
        constraint("geom_cells") ;
      }

      // Add for a single cell.
      void calculate(Entity cell) {
        kSourceTerm[cell]+=kTemporalSourceTerm[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToKSourceTerm> registerTemporalToKSourceTerm ;

  // Rule to compute the diagonal term for the linear system.
  class ComputeKMatrixDiagonal : public pointwise_rule {
    private:
      const_store<real> kMainCoefficient ;
      store<real> D ;
    public:

      // Define input and output.
      ComputeKMatrixDiagonal() {
        name_store("kMainCoefficient(massFluxCorrected_p,rhoStar)",kMainCoefficient) ;
        name_store("kStar_D",D) ;
        input("kMainCoefficient(massFluxCorrected_p,rhoStar)") ;
        output("kStar_D") ;
        constraint("geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) { D[cell]=kMainCoefficient[cell] ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeKMatrixDiagonal> registerComputeKMatrixDiagonal ;

  // Rule to compute the diagonal term for the linear system for cell next
  // to wall function boundaries.
  class ComputeKMatrixDiagonalWallFunction : public pointwise_rule {
    private:
      store<real> D ;
    public:

      // Define input and output.
      ComputeKMatrixDiagonalWallFunction() {
        name_store("wallFunction::kStar_D",D) ;
        output("wallFunction::kStar_D") ;
        constraint("wallFunctionCells") ;
      }

      // Set main coefficient to a large number.
      void calculate(Entity cell) { D[cell]=1.0e30 ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeKMatrixDiagonalWallFunction>
    registerComputeKMatrixDiagonalWallFunction ;

  // Rule to copy kStar_D for periodic faces.
  class ComputeKMatrixDiagonalPeriodic : public pointwise_rule {
    private:
      const_Map cl,cr,pmap ;
      store<real> D ;
    public:

      // Define input and output.
      ComputeKMatrixDiagonalPeriodic() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pmap",pmap) ;
        name_store("kStar_D",D) ;
        input("pmap->cl->kStar_D") ;
        output("cr->kStar_D") ;
      }

      // For a face.
      void calculate(Entity face) { D[cr[face]]=D[cl[pmap[face]]] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; } } ;

  register_rule<ComputeKMatrixDiagonalPeriodic>
    registerComputeKMatrixDiagonalPeriodic ;

  // Rule to initialize the lower terms for the linear system. This is now
  // a parametric variable with X=massFlux.
  class InitializeKMatrixLower : public unit_rule {
    private:
      store<real> L ;
    public:

      // Define input and output.
      InitializeKMatrixLower() {
        name_store("kStar_L(X)",L) ;
        output("kStar_L(X)") ;
        constraint("internalFaces") ;
      }

      // Initialize for a single face.
      void calculate(Entity face) { L[face]=0.0 ; }

      // Initialize for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeKMatrixLower> registerInitializeKMatrixLower ;

  // Rule to add the first-order inviscid flux contribution to the lower terms
  // for the linear system.
  class FOUInviscidFluxToKMatrixLower : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      store<real> L ;
    public:

      // Define input and output.
      FOUInviscidFluxToKMatrixLower() {
        name_store("X",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("kStar_L(X)",L) ;
        input("X,gridMassFlux") ;
        output("kStar_L(X)") ;
        constraint("internalFaces") ;
      }

      // Increment the lower term for a single face. Note that the increment
      // is the negative of the one in streamUns since this coefficient is
      // for a term on the lhs of the equation. In streamUns, the coefficient
      // is associated with a term on the rhs.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        if(netMassFlux>0.0) L[face]-=netMassFlux ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToKMatrixLower>
    registerFOUInviscidFluxToKMatrixLower ;

  // Rule to add the diffusive flux contribution to the lower terms for the
  // linear system.
  class DiffusiveFluxToKMatrixLower : public apply_rule<store<real>,Loci::
  Summation<real> > {
    private:
      const_store<real> diffusionProduct ;
      const_store<real> faceRadius ;
      const_store<real> kViscosity ;
      store<real> L ;
    public:

      // Define input and output.
      DiffusiveFluxToKMatrixLower() {
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("faceRadius",faceRadius) ;
        name_store("kViscosity",kViscosity) ;
        name_store("kStar_L(X)",L) ;
        input("diffusionProduct,faceRadius,kViscosity") ;
        output("kStar_L(X)") ;
        constraint("internalFaces,menterTurbulenceModel") ;
      }

      // Increment the lower term for a single face.
      void calculate(Entity face) {
        L[face]-=kViscosity[face]*diffusionProduct[face]*faceRadius[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToKMatrixLower>
    registerDiffusiveFluxToKMatrixLower ;

  // Rule to initialize the upper terms for the linear system. This is now
  // a parametric variable with X=massFlux.
  class InitializeKMatrixUpper : public unit_rule {
    private:
      store<real> U ;
    public:

      // Define input and output.
      InitializeKMatrixUpper() {
        name_store("kStar_U(X)",U) ;
        output("kStar_U(X)") ;
        constraint("internalFaces") ;
      }

      // Initialize for a single face.
      void calculate(Entity face) { U[face]=0.0 ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeKMatrixUpper> registerInitializeKMatrixUpper ;

  // Rule to add the first-order inviscid flux contribution to the upper terms
  // for the linear system.
  class FOUInviscidFluxToKMatrixUpper : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      store<real> U ;
    public:

      // Define input and output.
      FOUInviscidFluxToKMatrixUpper() {
        name_store("X",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("kStar_U(X)",U) ;
        input("X,gridMassFlux") ;
        output("kStar_U(X)") ;
        constraint("internalFaces") ;
      }

      // Increment the upper term for a single face. Contribution is opposite
      // to that in streamUns, as noted above for the lower terms.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        if(netMassFlux<0.0) U[face]+=netMassFlux ;
      }

      // Increment the upper term for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToKMatrixUpper>
    registerFOUInviscidFluxToKMatrixUpper ;

  // Rule to add the diffusive flux contribution to the upper terms for the
  // linear system.
  class DiffusiveFluxToKMatrixUpper : public apply_rule<store<real>,Loci::
  Summation<real> > {
    private:
      const_store<real> diffusionProduct ;
      const_store<real> faceRadius ;
      const_store<real> kViscosity ;
      store<real> U ;
    public:

      // Define input and output.
      DiffusiveFluxToKMatrixUpper() {
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("faceRadius",faceRadius) ;
        name_store("kViscosity",kViscosity) ;
        name_store("kStar_U(X)",U) ;
        input("diffusionProduct,faceRadius,kViscosity") ;
        output("kStar_U(X)") ;
        constraint("internalFaces,menterTurbulenceModel") ;
      }

      // Increment the lower term for a single face.
      void calculate(Entity face) {
        U[face]-=kViscosity[face]*diffusionProduct[face]*faceRadius[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToKMatrixUpper>
    registerDiffusiveFluxToKMatrixUpper ;

  // Rule to compute the right-hand side for the linear system.
  class ComputeKRHS : public pointwise_rule {
    private:
      const_store<real> kSourceTerm ;
      store<real> B ;
    public:

      // Define input and output.
      ComputeKRHS() {
        name_store("kSourceTerm(vCorrected_p,massFluxCorrected_p)",kSourceTerm) ;
        name_store("kStar_B",B) ;
        input("kSourceTerm(vCorrected_p,massFluxCorrected_p)") ;
        output("kStar_B") ;
        constraint("geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) { B[cell]=kSourceTerm[cell] ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeKRHS> registerComputeKRHS ;

  // Rule to copy kStar_B for periodic faces.
  class ComputeKRHSPeriodic : public pointwise_rule {
    private:
      const_Map cl,cr,pmap ;
      store<real> B ;
    public:

      // Define input and output.
      ComputeKRHSPeriodic() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pmap",pmap) ;
        name_store("kStar_B",B) ;
        input("pmap->cl->kStar_B") ;
        output("cr->kStar_B") ;
      }

      // For a face.
      void calculate(Entity face) { B[cr[face]]=B[cl[pmap[face]]] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeKRHSPeriodic> registerComputeKRHSPeriodic ;

  // Temporary hack to get around possible bug in using maps to set priority
  // variable values.
  class KCellHack : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> kWall ;
      store<real> kCell ;
    public:

      // Define input and output.
      KCellHack() {
        name_store("ci",ci) ;
        name_store("kWall",kWall) ;
        name_store("kCell",kCell) ;
        input("kWall") ;
        output("ci->kCell") ;
        constraint("ref->wallFunction_BCoption") ;
      }

      // Compute k for a single face.
      void calculate(Entity face) { kCell[ci[face]]=kWall[face] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<KCellHack> registerKCellHack ;

  // Rule to compute the right-hand side for the linear system for cells next
  // to wall function boundaries.
  class ComputeKRHSWallFunction : public pointwise_rule {
    private:
      const_store<real> kCell ;
      store<real> B ;
    public:

      // Define input and output.
      ComputeKRHSWallFunction() {
        name_store("kCell",kCell) ;
        name_store("wallFunction::kStar_B",B) ;
        input("kCell") ;
        output("wallFunction::kStar_B") ;
        constraint("wallFunctionCells") ;
      }

      // Set right-hand side to large number times kCell.
      void calculate(Entity cell) { B[cell]=1.0e30*kCell[cell] ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeKRHSWallFunction> registerComputeKRHSWallFunction ;

  // Rule to compute the production term for the k equation for incompressible
  // flow. This is now a parametric variable, with X=velocity.
  class ComputeKProductionIncompressible : public pointwise_rule {
    private:
      const_store<real> rho,k,eddyViscosity ;
      const_store<tens3d> vGradient ;
      const_store<real> vol ;
      store<real> kProduction ;
    public:

      // Define input and output.
      ComputeKProductionIncompressible() {
        name_store("rho",rho) ;
        name_store("k",k) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("gradv3d(X)",vGradient) ;
        name_store("vol",vol) ;
        name_store("kProduction(X)",kProduction) ;
        input("rho,k,eddyViscosity,gradv3d(X),vol") ;
        output("kProduction(X)") ;
        constraint("incompressibleFlow,geom_cells") ;
        constraint("menterBSLSSTTurbulenceModel") ;
      }

      // Compute for a single cell.
      void calculate(Entity cell) {
        tens3d stress=product(eddyViscosity[cell],vGradient[cell]+
          Transpose(vGradient[cell])) ;
        real temp=2.0*rho[cell]*k[cell]/3.0 ;
        stress.x.x-=temp ; stress.y.y-=temp ; stress.z.z-=temp ;
        kProduction[cell]=ScalarProduct(stress,vGradient[cell])*vol[cell] ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeKProductionIncompressible>
    registerComputeKProductionIncompressible ;
 
  // Rule to compute the production term for the k equation for incompressible
  // flow when using Menter's 2003 version of SST.
  class ComputeKProductionSST2003Incompressible : public pointwise_rule {
    private:
      const_param<real> betaStar ;
      const_store<real> rho,k,omega,eddyViscosity ;
      const_store<tens3d> vGradient ;
      const_store<real> vol ;
      store<real> kProduction ;
    public:

      // Define input and output.
      ComputeKProductionSST2003Incompressible() {
        name_store("betaStar",betaStar) ;
        name_store("rho",rho) ;
        name_store("k",k) ;
        name_store("omega",omega) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("gradv3d(X)",vGradient) ;
        name_store("vol",vol) ;
        name_store("kProduction(X)",kProduction) ;
        input("betaStar,rho,k,omega,eddyViscosity,gradv3d(X),vol") ;
        output("kProduction(X)") ;
        constraint("incompressibleFlow,geom_cells") ;
        constraint("menterSST2003TurbulenceModel") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) {
        tens3d stress=product(eddyViscosity[cell],vGradient[cell]+
          Transpose(vGradient[cell])) ;
        real temp=2.0*rho[cell]*k[cell]/3.0 ;
        stress.x.x-=temp ; stress.y.y-=temp ; stress.z.z-=temp ;
        kProduction[cell]=min(ScalarProduct(stress,vGradient[cell]),
          10.0*(*betaStar)*rho[cell]*omega[cell]*k[cell])*vol[cell] ;
      }

      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeKProductionSST2003Incompressible>
    registerComputeKProductionSST2003Incompressible ;

  // Rule to compute the production term for the k equation.
  class ComputeKProductionCompressible : public pointwise_rule {
    private:
      const_store<real> rho,k,eddyViscosity ;
      const_store<tens3d> vGradient ;
      const_store<real> vol ;
      store<real> kProduction ;
    public:

      // Define input and output.
      ComputeKProductionCompressible() {
        name_store("rho",rho) ;
        name_store("k",k) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("gradv3d(X)",vGradient) ;
        name_store("vol",vol) ;
        name_store("kProduction(X)",kProduction) ;
        input("rho,k,eddyViscosity,gradv3d(X),vol") ;
        output("kProduction(X)") ;
        constraint("compressibleFlow,geom_cells") ;
        constraint("menterBSLSSTTurbulenceModel") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) {
        tens3d stress=product(eddyViscosity[cell],vGradient[cell]+
          Transpose(vGradient[cell])) ;
        real temp=2.0*(rho[cell]*k[cell]+eddyViscosity[cell]*
          Trace(vGradient[cell]))/3.0 ; ;
        stress.x.x-=temp ; stress.y.y-=temp ; stress.z.z-=temp ;
        kProduction[cell]=ScalarProduct(stress,vGradient[cell])*vol[cell] ;
      }

      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeKProductionCompressible>
    registerComputeKProductionCompressible ;
 
  // Rule to compute the production term for the k equation when using
  // Menter's 2003 version of SST
  class ComputeKProductionSST2003Compressible : public pointwise_rule {
    private:
      const_param<real> betaStar ;
      const_store<real> rho,k,omega,eddyViscosity ;
      const_store<tens3d> vGradient ;
      const_store<real> vol ;
      store<real> kProduction ;
    public:

      // Define input and output.
      ComputeKProductionSST2003Compressible() {
        name_store("betaStar",betaStar) ;
        name_store("rho",rho) ;
        name_store("k",k) ;
        name_store("omega",omega) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("gradv3d(X)",vGradient) ;
        name_store("vol",vol) ;
        name_store("kProduction(X)",kProduction) ;
        input("betaStar,rho,k,omega,eddyViscosity,gradv3d(X),vol") ;
        output("kProduction(X)") ;
        constraint("compressibleFlow,geom_cells") ;
        constraint("menterSST2003TurbulenceModel") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) {
        tens3d stress=product(eddyViscosity[cell],vGradient[cell]+
          Transpose(vGradient[cell])) ;
        real temp=2.0*(rho[cell]*k[cell]+eddyViscosity[cell]*
          Trace(vGradient[cell]))/3.0 ;
        stress.x.x-=temp ; stress.y.y-=temp ; stress.z.z-=temp ;
        kProduction[cell]=min(ScalarProduct(stress,vGradient[cell]),
          10.0*(*betaStar)*rho[cell]*k[cell]*omega[cell])*vol[cell] ;
      }

      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeKProductionSST2003Compressible>
    registerComputeKProductionSST2003Compressible ;

  // Rule to compute the blending function F1. Checked.
  class ComputeF1 : public pointwise_rule {
    private:
      const_param<real> sigmaOmega2 ;
      const_store<real> rho,k,omega,laminarViscosity ;
      const_store<vect3d> kGradient,omegaGradient ;
      const_store<real> distanceToNoSlip ;
      store<real> f1 ;
    public:

      // Define input and output.
      ComputeF1() {
        name_store("sigmaOmega2",sigmaOmega2) ;
        name_store("rho",rho) ;
        name_store("k",k) ;
        name_store("omega",omega) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("grads(k)",kGradient) ;
        name_store("grads(omega)",omegaGradient) ;
        name_store("distanceToNoSlip",distanceToNoSlip) ;
        name_store("f1",f1) ;
        input("sigmaOmega2,rho,k,omega,laminarViscosity") ;
        input("grads(k),grads(omega),distanceToNoSlip") ;
        output("f1") ;
        constraint("geom_cells") ;
      }

      // Compute F1 for a single cell.
      void calculate(Entity cell) {

        // Compute the cross-diffusion term.
        real crossDiffusion ;
        {
          real argA=2.0*rho[cell]*(*sigmaOmega2)*dot(kGradient[cell],
            omegaGradient[cell])/omega[cell] ;
          crossDiffusion=(argA>1.0e-20)? argA:1.0e-20 ;
        }

        // Finish the computation of F1.
        real argA=sqrt(k[cell])/(0.09*omega[cell]*distanceToNoSlip[cell]) ;
        real argB=500.0*laminarViscosity[cell]/(rho[cell]*distanceToNoSlip
          [cell]*distanceToNoSlip[cell]*omega[cell]) ;
        real argC=4.0*rho[cell]*(*sigmaOmega2)*k[cell]/(crossDiffusion*
          distanceToNoSlip[cell]*distanceToNoSlip[cell]) ;
        f1[cell]=(argA>argB)? argA:argB ;
        f1[cell]=(f1[cell]<argC)? tanh(f1[cell]*f1[cell]*f1[cell]*f1[cell]):
          tanh(argC*argC*argC*argC) ;
      }

      // Compute F1 for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeF1> registerComputeF1 ;
 
  // Rule to compute the blending function F1 for Menter's SST model of 2003.
  class ComputeF1SST2003 : public pointwise_rule {
    private:
      const_param<real> sigmaOmega2 ;
      const_store<real> rho,k,omega,laminarViscosity ;
      const_store<vect3d> kGradient,omegaGradient ;
      const_store<real> distanceToNoSlip ;
      store<real> f1 ;
    public:

      // Define input and output.
      ComputeF1SST2003() {
        name_store("sigmaOmega2",sigmaOmega2) ;
        name_store("rho",rho) ;
        name_store("k",k) ;
        name_store("omega",omega) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("grads(k)",kGradient) ;
        name_store("grads(omega)",omegaGradient) ;
        name_store("distanceToNoSlip",distanceToNoSlip) ;
        name_store("menterSST2003::f1",f1) ;
        input("sigmaOmega2,rho,k,omega,laminarViscosity") ;
        input("grads(k),grads(omega),distanceToNoSlip") ;
        output("menterSST2003::f1") ;
        constraint("geom_cells,menterSST2003TurbulenceModel") ;
      }

      // Compute F1 for a single cell.
      void calculate(Entity cell) {

        // Compute the cross-diffusion term.
        real crossDiffusion ;
        {
          real argA=2.0*rho[cell]*(*sigmaOmega2)*dot(kGradient[cell],
            omegaGradient[cell])/omega[cell] ;
          crossDiffusion=(argA>1.0e-10)? argA:1.0e-10 ;
        }

        // Finish the computation of F1.
        real argA=sqrt(k[cell])/(0.09*omega[cell]*distanceToNoSlip[cell]) ;
        real argB=500.0*laminarViscosity[cell]/(rho[cell]*distanceToNoSlip
          [cell]*distanceToNoSlip[cell]*omega[cell]) ;
        real argC=4.0*rho[cell]*(*sigmaOmega2)*k[cell]/(crossDiffusion*
          distanceToNoSlip[cell]*distanceToNoSlip[cell]) ;
        f1[cell]=(argA>argB)? argA:argB ;
        f1[cell]=(f1[cell]<argC)? tanh(f1[cell]*f1[cell]*f1[cell]*f1[cell]):
          tanh(argC*argC*argC*argC) ;
      }

      // Compute F1 for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeF1SST2003> registerComputeF1SST2003 ;

  // Rule to compute the function F2 used in the modified definition of the
  // kinematic eddy viscosity for the SST model.
  class ComputeF2 : public pointwise_rule {
    private:
      const_store<real> rho,k,omega,laminarViscosity ;
      const_store<real> distanceToNoSlip ;
      store<real> f2 ;
    public:

      // Define input and output.
      ComputeF2() {
        name_store("rho",rho) ;
        name_store("k",k) ;
        name_store("omega",omega) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("distanceToNoSlip",distanceToNoSlip) ;
        name_store("f2",f2) ;
        input("rho,k,omega,laminarViscosity") ;
        input("distanceToNoSlip") ;
        output("f2") ;
        constraint("geom_cells") ;
      }

      // Compute F2 for a single cell.
      void calculate(Entity cell) {
        real argA=2.0*sqrt(k[cell])/(0.09*omega[cell]*distanceToNoSlip[cell]) ;
        real argB=500.0*laminarViscosity[cell]/(rho[cell]*distanceToNoSlip
          [cell]*distanceToNoSlip[cell]*omega[cell]) ;
        f2[cell]=(argA>argB)? tanh(argA*argA):tanh(argB*argB) ;
      }

      // Compute F2 for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeF2> registerComputeF2 ;

//-----------------------------------------------------------------------------
// Rules for computing the residual of the k equation.

  // Rule to initialize the k residual.
  class InitializeKResidual : public unit_rule {
    private:
      store<real> kResidual ;
    public:

      // Define input and output.
      InitializeKResidual() {
        name_store("kResidual",kResidual) ;
        output("kResidual") ;
        constraint("vol") ;
      }

      // Initialize the residual for a single cell.
      void calculate(Entity cell) { kResidual[cell]=0.0 ; }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

//register_rule<InitializeKResidual> registerInitializeKResidual ;

  // Rule to compute the k residual for each cell.
  class ComputeKResidualOne : public apply_rule<store<real>,Loci::Summation
  <real> > {
    private:
      const_param<real> rhoScale,vScale,lScale,kScale ;
      const_store<real> D ;
      const_store<real> k ;
      const_store<real> B ;
      store<real> kResidual ;
    private:
      real kFactor ;
    public:

      // Define input and output.
      ComputeKResidualOne() {
        name_store("rhoScale",rhoScale) ;
        name_store("vScale",vScale) ;
        name_store("lScale",lScale) ;
        name_store("kScale",kScale) ;
        name_store("kStar_D",D) ;
        name_store("k",k) ;
        name_store("kStar_B",B) ;
        name_store("kResidual",kResidual) ;
        input("rhoScale,vScale,lScale,kScale,kStar_D,k,kStar_B") ;
        output("kResidual") ;
        constraint("vol") ;
      }

      // Initialize the residual for a single cell.
      void calculate(Entity cell) {
        kResidual[cell]+=(D[cell]!=1.0e30)? (B[cell]-D[cell]*k[cell])/
          kFactor:0.0 ;
      }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) {
        kFactor=(*rhoScale)*(*vScale)*(*kScale)*(*lScale)*(*lScale) ;
        do_loop(seq,this) ;
      }
  } ;

//register_rule<ComputeKResidualOne> registerComputeKResidualOne ;

  // Rule to compute the k residual for each cell.
  class ComputeKResidualTwo : public apply_rule<store<real>,Loci::Summation
  <real> > {
    private:
      const_Map cl,cr ;
      const_param<real> rhoScale,vScale,lScale,kScale ;
      const_store<real> k ;
      const_store<real> D,L,U ;
      store<real> kResidual ;
    private:
      real kFactor ;
    public:

      // Define input and output.
      ComputeKResidualTwo() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("rhoScale",rhoScale) ;
        name_store("vScale",vScale) ;
        name_store("lScale",lScale) ;
        name_store("kScale",kScale) ;
        name_store("k",k) ;
        name_store("kStar_D",D) ;
        name_store("kStar_L",L) ;
        name_store("kStar_U",U) ;
        name_store("kResidual",kResidual) ;
        input("rhoScale,vScale,lScale,kScale") ;
        input("(cl,cr)->(k,kStar_D),kStar_L,kStar_U") ;
        output("(cl,cr)->kResidual") ;
        constraint("internalFaces") ;
      }

      // Add the neighbor contribution to the residual for each of the two
      // cells on either side of the face.
      void calculate(Entity face) {
        kResidual[cl[face]]-=(D[cl[face]]!=1.0e30)? U[face]*k[cr[face]]/
          kFactor:0.0 ;
        kResidual[cr[face]]-=(D[cr[face]]!=1.0e30)? L[face]*k[cl[face]]/
          kFactor:0.0 ;
      }

      // Add the neighbor contributions for a sequence of faces.
      virtual void compute(const sequence &seq) {
        kFactor=(*rhoScale)*(*vScale)*(*kScale)*(*lScale)*(*lScale) ;
        do_loop(seq,this) ;
      }
  } ;

//register_rule<ComputeKResidualTwo> registerComputeKResidualTwo ;

  // Rule to set total k residual to zero for laminar flow.
  class TotalKResidualLaminar : public singleton_rule {
    private:
      param<ScalarResidual> kResidualData ;
    public:

      // Define input and output.
      TotalKResidualLaminar() {
        name_store("kResidualData",kResidualData) ;
        output("kResidualData") ;
        constraint("laminarFlow,geom_cells") ;
      }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) {
        kResidualData=ScalarResidual() ;
      }
  } ;

//register_rule<TotalKResidualLaminar> registerTotalKResidualLaminar ;

  // Rule to set total k residual to zero for inviscid flow.
  class TotalKResidualInviscid : public singleton_rule {
    private:
      param<ScalarResidual> kResidualData ;
    public:

      // Define input and output.
      TotalKResidualInviscid() {
        name_store("kResidualData",kResidualData) ;
        output("kResidualData") ;
        constraint("inviscidFlow,geom_cells") ;
      }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) {
        kResidualData=ScalarResidual() ;
      }
  } ;

//register_rule<TotalKResidualInviscid> registerTotalKResidualInviscid ;

  // Rule to initialize the total k residual.
  class InitializeTotalKResidual : public unit_rule {
    private:
      param<ScalarResidual> kResidualData ;
    public:

      // Define input and output.
      InitializeTotalKResidual() {
        name_store("kResidualData",kResidualData) ;
        output("kResidualData") ;
        constraint("turbulentFlow,geom_cells") ;
      }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) {
        kResidualData=ScalarResidual() ;
      }
  } ;

//register_rule<InitializeTotalKResidual> registerInitializeTotalKResidual ;

  // Rule to compute the total k residual. Note the use of kStar_D to determine
  // the wall-function cells to exclude from the sum.
  class ComputeTotalKResidual : public apply_rule<param<ScalarResidual>,
  ScalarResidualJoin> {
    private:
      const_store<real> kResidual ;
      const_store<real> D ;
      const_store<vect3d> cellCenter ;
      param<ScalarResidual> kResidualData ;
    public:

      // Define input and output.
      ComputeTotalKResidual() {
        name_store("kResidual",kResidual) ;
        name_store("kStar_D",D) ;
        name_store("cellcenter",cellCenter) ;
        name_store("kResidualData",kResidualData) ;
        input("kResidual,kStar_D,cellcenter") ;
        output("kResidualData") ;
        constraint("turbulentFlow,geom_cells") ;
      }

      // Add the cell contribution to the residual for a single cell.
      void calculate(Entity cell) {
        ScalarResidual temp ;
        temp.maxResidual=kResidual[cell] ;
        temp.totalResidual=abs(kResidual[cell]) ;
        temp.maxResidualLocation=cellCenter[cell] ;
        join(*kResidualData,temp) ;
      }

      // Add the cell contribution to the residual for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

//register_rule<ComputeTotalKResidual> registerComputeTotalKResidual ;

//-----------------------------------------------------------------------------
// Scheme independent rules for assembling the omega equation.

  // Rule to initialize the main coefficient.  This is now a parametric
  // variable, with X=massFlux and Y=density.
  class InitializeOmegaMainCoefficient : public unit_rule {
    private:
      store<real> omegaMainCoefficient ;
    public:

      // Define input and output.
      InitializeOmegaMainCoefficient() {
        name_store("omegaMainCoefficient(X,Y)",omegaMainCoefficient) ;
        output("omegaMainCoefficient(X,Y)") ;
        constraint("vol") ;
      }

      // Set the main coefficient to zero for a single cell.
      void calculate(Entity cell) { omegaMainCoefficient[cell]=0.0 ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeOmegaMainCoefficient>
    registerInitializeOmegaMainCoefficient ;

  // Rule to add the first-order inviscid flux contribution to the main
  // coefficient for interior faces.
  class FOUInviscidFluxToOmegaMainCoefficientInterior : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      store<real> omegaMainCoefficient ;
    public:

      // Define input and output.
      FOUInviscidFluxToOmegaMainCoefficientInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("X",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("omegaMainCoefficient(X,Y)",omegaMainCoefficient) ;
        input("X,gridMassFlux") ;
        output("(cl,cr)->omegaMainCoefficient(X,Y)") ;
        constraint("internalFaces") ;
      }

      // Increment the main coefficient for cells attached to a single face.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        if(netMassFlux>0.0){
          omegaMainCoefficient[cl[face]]+=netMassFlux ;
        }else{
          omegaMainCoefficient[cr[face]]-=netMassFlux ;
        }
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToOmegaMainCoefficientInterior>
    registerFOUInviscidFluxToOmegaMainCoefficientInterior ;

  // Rule to add the first-order inviscid flux contribution to the main
  // coefficient for boundary faces.
  class FOUInviscidFluxToOmegaMainCoefficientBoundary : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      store<real> omegaMainCoefficient ;
    public:

      // Define input and output.
      FOUInviscidFluxToOmegaMainCoefficientBoundary() {
        name_store("ci",ci) ;
        name_store("X",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("omegaMainCoefficient(X,Y)",omegaMainCoefficient) ;
        input("X,gridMassFlux") ;
        output("ci->omegaMainCoefficient(X,Y)") ;
        constraint("boundaryFaces") ;
      }

      // Increment the main coefficient for the cell attach to a single face.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        if(netMassFlux>0.0) omegaMainCoefficient[ci[face]]+=netMassFlux ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToOmegaMainCoefficientBoundary>
    registerFOUInviscidFluxToOmegaMainCoefficientBoundary ;

  // Rule to add the destruction term to the main coefficient.
  class DestructionToOmegaMainCoefficient : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> beta1,beta2 ;
      const_store<real> f1 ;
      const_store<real> rho ;
      const_store<real> omega ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> omegaMainCoefficient ;
    public:

      // Define input and output.
      DestructionToOmegaMainCoefficient() {
        name_store("beta1",beta1) ;
        name_store("beta2",beta2) ;
        name_store("f1",f1) ;
        name_store("rho",rho) ;
        name_store("omega",omega) ;
        name_store("vol",vol) ;
        name_store("cellRadius",cellRadius) ;
        name_store("omegaMainCoefficient(X,Y)",omegaMainCoefficient) ;
        input("beta1,beta2,f1,rho,omega,vol,cellRadius") ;
        output("omegaMainCoefficient(X,Y)") ;
        constraint("geom_cells,menterTurbulenceModel") ;
      }

      // Increment the main coefficient for a single cell.
      void calculate(Entity cell) {
        real beta=f1[cell]*(*beta1)+(1.0-f1[cell])*(*beta2) ;
        omegaMainCoefficient[cell]+=beta*rho[cell]*omega[cell]*vol[cell]*
          cellRadius[cell] ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DestructionToOmegaMainCoefficient>
    registerDestructionToOmegaMainCoefficient ;

  // Rule to add the diffusive flux contribution to the main coefficient for
  // interior faces.
  class DiffusiveFluxToOmegaMainCoefficientInterior : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> faceRadius ;
      const_store<real> diffusionProduct ;
      const_store<real> omegaViscosity ;
      store<real> omegaMainCoefficient ;
    public:

      // Define input and output.
      DiffusiveFluxToOmegaMainCoefficientInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("faceRadius",faceRadius) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("omegaViscosity",omegaViscosity) ;
        name_store("omegaMainCoefficient(X,Y)",omegaMainCoefficient) ;
        input("faceRadius,diffusionProduct,omegaViscosity") ;
        output("(cl,cr)->omegaMainCoefficient(X,Y)") ;
        constraint("internalFaces") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real temp=omegaViscosity[face]*diffusionProduct[face]*faceRadius[face] ;
        omegaMainCoefficient[cl[face]]+=temp ;
        omegaMainCoefficient[cr[face]]+=temp ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToOmegaMainCoefficientInterior>
    registerDiffusiveFluxToOmegaMainCoefficientInterior ;

  // Rule to add the diffusive flux contribution to the main coefficient for
  // boundary faces.
  class DiffusiveFluxToOmegaMainCoefficientBoundary : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_param<real> sigmaOmega1,sigmaOmega2 ;
      const_store<real> faceRadius ;
      const_store<real> diffusionProduct ;
      const_store<real> laminarViscosity_f ;
      const_store<real> eddyViscosity_f ;
      const_store<real> f1 ;
      store<real> omegaMainCoefficient ;
    public:

      // Define input and output.
      DiffusiveFluxToOmegaMainCoefficientBoundary() {
        name_store("ci",ci) ;
        name_store("sigmaOmega1",sigmaOmega1) ;
        name_store("sigmaOmega2",sigmaOmega2) ;
        name_store("faceRadius",faceRadius) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        name_store("eddyViscosity_f",eddyViscosity_f) ;
        name_store("f1",f1) ;
        name_store("omegaMainCoefficient(X,Y)",omegaMainCoefficient) ;
        input("sigmaOmega1,sigmaOmega2,ci->f1,faceRadius") ;
        input("diffusionProduct,laminarViscosity_f,eddyViscosity_f") ;
        output("ci->omegaMainCoefficient(X,Y)") ;
        constraint("boundaryTurbulenceDiffusion") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real sigmaOmega=f1[ci[face]]*(*sigmaOmega1)+(1.0-f1[ci[face]])*
          (*sigmaOmega2) ;
        omegaMainCoefficient[ci[face]]+=(laminarViscosity_f[face]+sigmaOmega*
          eddyViscosity_f[face])*diffusionProduct[face]*faceRadius[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToOmegaMainCoefficientBoundary>
    registerDiffusiveFluxToOmegaMainCoefficientBoundary ;

  // Rule to add temporal component of momentum equation to the main
  // coefficient. This is now an iteration independent rule so that it can
  // be used at both {n} and {n,it} where it is needed.
  class TemporalToOmegaMainCoefficient: public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> timeStep ;
      const_param<real> timeIntegratorFactor0 ;
      const_store<real> timeStepFactor ;
      const_store<real> rho ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> omegaMainCoefficient ;
    public:

      // Define input and output.
      TemporalToOmegaMainCoefficient() {
        name_store("timeStep",timeStep) ;
        name_store("timeIntegratorFactor0",timeIntegratorFactor0) ;
        name_store("timeStepFactor",timeStepFactor) ;
        name_store("Y",rho) ;
        name_store("vol",vol) ;
        name_store("cellRadius",cellRadius) ;
        name_store("omegaMainCoefficient(X,Y)",omegaMainCoefficient) ;
        input("timeStep,timeStepFactor,timeIntegratorFactor0") ;
        input("Y,vol,cellRadius") ;
        output("omegaMainCoefficient(X,Y)") ;
        constraint("geom_cells") ;
      }

      // Add temporal component for a single cell.
      void calculate(Entity cell) {
        omegaMainCoefficient[cell]+=0.5*rho[cell]*vol[cell]*cellRadius[cell]*
          ((*timeIntegratorFactor0)+1.0)/((*timeStep)*timeStepFactor[cell]) ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToOmegaMainCoefficient>
    registerTemporalToOmegaMainCoefficient ;

  // Rule to initialize the source term.  This is now a parametric variable
  // with X=velocity and Y=massFlux.
  class InitializeOmegaSourceTerm : public unit_rule {
    private:
      store<real> omegaSourceTerm ;
    public:

      // Define input and output.
      InitializeOmegaSourceTerm() {
        name_store("omegaSourceTerm(X,Y)",omegaSourceTerm) ;
        output("omegaSourceTerm(X,Y)") ;
        constraint("vol") ;
      }

      // Set the source term to zero for a single cell.
      void calculate(Entity cell) { omegaSourceTerm[cell]=0.0 ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeOmegaSourceTerm> registerInitializeOmegaSourceTerm ;

  // Rule to add the first-order inviscid flux contribution to the source term
  // for boundary faces.
  class FOUInviscidFluxToOmegaSourceTermBoundary : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> omega ;
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      const_store<real> omega_f ;
      store<real> omegaSourceTerm ;
    public:

      // Define input and output.
      FOUInviscidFluxToOmegaSourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("omega",omega) ;
        name_store("Y",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("omega_f",omega_f) ;
        name_store("omegaSourceTerm(X,Y)",omegaSourceTerm) ;
        input("ci->omega,Y,gridMassFlux,omega_f") ;
        output("ci->omegaSourceTerm(X,Y)") ;
        constraint("boundaryFaces") ;
      }

      // Increment the source term for the cell attach to a single face.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        if(netMassFlux<0.0) omegaSourceTerm[ci[face]]-=netMassFlux*
          omega_f[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToOmegaSourceTermBoundary>
    registerFOUInviscidFluxToOmegaSourceTermBoundary ;

  // Rule to add the second-order convection contribution to the source term
  // for interior faces.
  class SOUInviscidFluxToOmegaSourceTermInterior : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> omega ;
      const_store<vect3d> omegaGradient ;
      const_store<real> omegaLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      store<real> omegaSourceTerm ;
    public:

      // Define input and output.
      SOUInviscidFluxToOmegaSourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("omega",omega) ;
        name_store("grads(omega)",omegaGradient) ;
        name_store("limiters(omega)",omegaLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("Y",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("omegaSourceTerm(X,Y)",omegaSourceTerm) ;
        input("(cl,cr)->(omega,grads(omega),limiters(omega),cellcenter)") ;
        input("facecenter,Y,gridMassFlux") ;
        output("(cl,cr)->omegaSourceTerm(X,Y)") ;
        constraint("internalFaces,souTurbulence") ;
      }

      // Increment the source term for the cells attach to a single face. Note
      // that here we are insisting that the reconstructed value of omega at the
      // face be positive.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        real secondOrderSource=(netMassFlux>0.0)? netMassFlux*
          max(omegaLimiter[cl[face]]*dot(omegaGradient[cl[face]],
          faceCenter[face]-cellCenter[cl[face]]),-omega[cl[face]]):
          netMassFlux*max(omegaLimiter[cr[face]]*dot(omegaGradient[cr[face]],
          faceCenter[face]-cellCenter[cr[face]]),-omega[cr[face]]) ;
        omegaSourceTerm[cl[face]]-=secondOrderSource ;
        omegaSourceTerm[cr[face]]+=secondOrderSource ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SOUInviscidFluxToOmegaSourceTermInterior>
    registerSOUInviscidFluxToOmegaSourceTermInterior ;

  // Rule to add the second-order convection contribution to the source term
  // for boundary faces.
  class SOUInviscidFluxToOmegaSourceTermBoundary : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> omega ;
      const_store<vect3d> omegaGradient ;
      const_store<real> omegaLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      const_store<real> omega_f ;
      store<real> omegaSourceTerm ;
    public:

      // Define input and output.
      SOUInviscidFluxToOmegaSourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("omega",omega) ;
        name_store("grads(omega)",omegaGradient) ;
        name_store("limiters(omega)",omegaLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("Y",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("omega_f",omega_f) ;
        name_store("omegaSourceTerm(X,Y)",omegaSourceTerm) ;
        input("ci->(omega,grads(omega),limiters(omega),cellcenter)") ;
        input("facecenter,Y,gridMassFlux,omega_f") ;
        output("ci->omegaSourceTerm(X,Y)") ;
        constraint("boundaryFaces,souTurbulence") ;
      }

      // Increment the source term for the cells attach to a single face. Note
      // that here we are insisting that the reconstructed value of omega at the
      // face be positive.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        if(netMassFlux>0.0) omegaSourceTerm[ci[face]]-=netMassFlux*
          max(omegaLimiter[ci[face]]*dot(omegaGradient[ci[face]],
          faceCenter[face]-cellCenter[ci[face]]),-omega[ci[face]]) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SOUInviscidFluxToOmegaSourceTermBoundary>
    registerSOUInviscidFluxToOmegaSourceTermBoundary ;

  // Rule to add the production and cross-diffusion terms to the source term.
  // Checked.
  class ProductionAndCrossDiffusionToOmegaSourceTerm : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_param<real> beta1,beta2,gamma1,gamma2,sigmaOmega2 ;
      const_store<real> rho ;
      const_store<real> omega ;
      const_store<real> eddyViscosity ;
      const_store<real> kProduction ;
      const_store<vect3d> kGradient ;
      const_store<vect3d> omegaGradient ;
      const_store<real> f1 ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> omegaSourceTerm ;
    public:

      // Define input and output.
      ProductionAndCrossDiffusionToOmegaSourceTerm() {
        name_store("beta1",beta1) ;
        name_store("beta2",beta2) ;
        name_store("gamma1",gamma1) ;
        name_store("gamma2",gamma2) ;
        name_store("sigmaOmega2",sigmaOmega2) ;
        name_store("rho",rho) ;
        name_store("omega",omega) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("kProduction(X)",kProduction) ;
        name_store("grads(k)",kGradient) ;
        name_store("grads(omega)",omegaGradient) ;
        name_store("f1",f1) ;
        name_store("vol",vol) ;
        name_store("cellRadius",cellRadius) ;
        name_store("omegaSourceTerm(X,Y)",omegaSourceTerm) ;
        input("beta1,beta2,gamma1,gamma2,sigmaOmega2") ;
        input("rho,omega,eddyViscosity,kProduction(X)") ;
        input("grads(k),grads(omega),f1,vol,cellRadius") ;
        output("omegaSourceTerm(X,Y)") ;
        constraint("geom_cells,menterTurbulenceModel") ;
      }

      // Increment the source term for a single cell. Note that kProduction
      // already has been multiplied by the volume. Note the new stunt to
      // prevent dividing by a zero eddy viscosity. Added 4/2/2009 .
      void calculate(Entity cell) {
        if(eddyViscosity[cell]>1.0e-10){
        real gamma=f1[cell]*(*gamma1)+(1.0-f1[cell])*(*gamma2) ;
        omegaSourceTerm[cell]+=(gamma*kProduction[cell]*rho[cell]/
          eddyViscosity[cell]+2.0*(1.0-f1[cell])*rho[cell]*(*sigmaOmega2)*
          dot(kGradient[cell],omegaGradient[cell])*vol[cell]/omega[cell])*
          cellRadius[cell] ;
        }
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ProductionAndCrossDiffusionToOmegaSourceTerm>
    registerProductionAndCrossDiffusionToOmegaSourceTerm ;

  // Rule to add the diffusive flux contribution to the source term for
  // interior faces.
  class DiffusiveFluxToOmegaSourceTermInterior : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> vol ;
      const_store<vect3d> omegaGradient ;
      const_store<vect3d> geometryFactor0 ;
      const_store<real> faceRadius ;
      const_store<real> omegaViscosity ;
      store<real> omegaSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToOmegaSourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("vol",vol) ;
        name_store("grads(omega)",omegaGradient) ;
        name_store("geometryFactor0",geometryFactor0) ;
        name_store("faceRadius",faceRadius) ;
        name_store("omegaViscosity",omegaViscosity) ;
        name_store("omegaSourceTerm(X,Y)",omegaSourceTerm) ;
        input("(cl,cr)->(vol,grads(omega))");
        input("geometryFactor0,faceRadius,omegaViscosity") ;
        output("(cl,cr)->omegaSourceTerm(X,Y)") ;
        constraint("internalFaces") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real secondarySourceTerm=omegaViscosity[face]*
          dot((omegaGradient[cl[face]]*vol[cr[face]]+omegaGradient[cr[face]]*
          vol[cl[face]]),geometryFactor0[face])*faceRadius[face]/
          (vol[cl[face]]+vol[cr[face]]) ;
        omegaSourceTerm[cl[face]]+=secondarySourceTerm ;
        omegaSourceTerm[cr[face]]-=secondarySourceTerm ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToOmegaSourceTermInterior>
    registerDiffusiveFluxToOmegaSourceTermInterior ;

  // Rule to add the diffusive flux contribution to the source term for
  // boundary faces.
  class DiffusiveFluxToOmegaSourceTermBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_param<real> sigmaOmega1,sigmaOmega2 ;
      const_store<real> f1 ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> omegaGradient ;
      const_store<vect3d> faceCenter ;
      const_store<real> diffusionProduct ;
      const_store<real> omega_f ;
      const_store<real> laminarViscosity_f ;
      const_store<real> eddyViscosity_f ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<real> omegaSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToOmegaSourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("sigmaOmega1",sigmaOmega1) ;
        name_store("sigmaOmega2",sigmaOmega2) ;
        name_store("f1",f1) ;
        name_store("cellcenter",cellCenter) ;
        name_store("grads(omega)",omegaGradient) ;
        name_store("facecenter",faceCenter) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("omega_f",omega_f) ;
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        name_store("eddyViscosity_f",eddyViscosity_f) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("omegaSourceTerm(X,Y)",omegaSourceTerm) ;
        input("sigmaOmega1,sigmaOmega2") ;
        input("ci->f1,ci->cellcenter,ci->grads(omega)") ;
        input("facecenter,diffusionProduct,omega_f") ;
        input("laminarViscosity_f,eddyViscosity_f,area,faceRadius") ;
        output("ci->omegaSourceTerm(X,Y)") ;
        constraint("boundaryTurbulenceDiffusion") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real sigmaOmega=f1[ci[face]]*(*sigmaOmega1)+(1.0-f1[ci[face]])*
          (*sigmaOmega2) ;
        real sourceTerm=(laminarViscosity_f[face]+sigmaOmega*eddyViscosity_f
          [face])*(omega_f[face]*diffusionProduct[face]+dot(omegaGradient
          [ci[face]],(area[face].n*area[face].sada-diffusionProduct[face]*
          (faceCenter[face]-cellCenter[ci[face]]))))*faceRadius[face] ;
        omegaSourceTerm[ci[face]]+=sourceTerm ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToOmegaSourceTermBoundary>
    registerDiffusiveFluxToOmegaSourceTermBoundary ;

  // Rule to compute the temporal source term for BDF. Note that this
  // term only includes the 'old' components from the {n} level. This
  // rule will compute and be used at {n}. In addition, the value will
  // promote to {n,it} where it is also needed.
  class OmegaTemporalSourceTermBDF : public pointwise_rule {
    private:
      const_param<real> timeStep ;
      const_store<real> timeStepFactor ;
      const_store<real> rho,vol,cellRadius ;
      const_store<real> omega ;
      store<real> omegaTemporalSourceTerm ;
    public:

      // Define input and output.
      OmegaTemporalSourceTermBDF() {
        name_store("timeStep{n}",timeStep) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n}",rho) ;
        name_store("omega{n}",omega) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n}",cellRadius) ;
        name_store("omegaTemporalSourceTerm{n}",omegaTemporalSourceTerm) ;
        input("rho{n},omega{n},vol{n},cellRadius{n}") ;
        input("timeStep{n},timeStepFactor{n}") ;
        output("omegaTemporalSourceTerm{n}") ;
        constraint("BDFIntegrator,geom_cells") ;
      }

      // Compute for a single cell.
      void calculate(Entity cell) {
        omegaTemporalSourceTerm[cell]=(rho[cell]*vol[cell]*cellRadius[cell]/
          ((*timeStep)*timeStepFactor[cell]))*omega[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<OmegaTemporalSourceTermBDF> registerOmegaTemporalSourceTermBDF ;

  // Rule to compute the temporal source term for BDF2. Note that this
  // term only includes the 'old' components from the {n-1} and {n} levels.
  // This rule will compute and be used at {n}. In addition, the value will
  // promote to {n,it} where it is also needed.
  class OmegaTemporalSourceTermBDF2 : public pointwise_rule {
    private:
      const_param<int> n ;
      const_param<real> timeStep ;
      const_store<real> timeStepFactor ;
      const_store<real> rhoOld,rho ;
      const_store<real> omegaOld,omega ;
      const_store<real> volOld,vol ;
      const_store<real> cellRadius ;
      store<real> omegaTemporalSourceTerm ;
    public:

      // Define input and output.
      OmegaTemporalSourceTermBDF2() {
        name_store("$n{n}",n) ;
        name_store("timeStep{n}",timeStep) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n-1}",rhoOld) ;
        name_store("rho{n}",rho) ;
        name_store("omega{n-1}",omegaOld) ;
        name_store("omega{n}",omega) ;
//      name_store("vol{n-1}",volOld) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n}",cellRadius) ;
        name_store("omegaTemporalSourceTerm{n}",omegaTemporalSourceTerm) ;
//      input("$n{n},rho{n-1},rho{n},omega{n-1},omega{n},vol{n-1},vol{n}") ;
        input("$n{n},rho{n-1},rho{n},omega{n-1},omega{n},vol{n}") ;
        input("cellRadius{n},timeStep{n},timeStepFactor{n}") ;
        output("omegaTemporalSourceTerm{n}") ;
        constraint("BDF2Integrator,geom_cells") ;
      }

      // Compute for a single cell. Use BDF on first timestep.
      void calculate(Entity cell) {
        if((*n)!=0){
          omegaTemporalSourceTerm[cell]=((2.0*rho[cell]*vol[cell])*omega[cell]-
           (0.5*rhoOld[cell]*vol[cell])*omegaOld[cell])*cellRadius[cell]/
           ((*timeStep)*timeStepFactor[cell]) ;
        }else{
          omegaTemporalSourceTerm[cell]=(rho[cell]*vol[cell]*cellRadius[cell]/
            ((*timeStep)*timeStepFactor[cell]))*omega[cell] ;
        }
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<OmegaTemporalSourceTermBDF2>
    registerOmegaTemporalSourceTermBDF2 ;

  // Rule to add temporal component of momentum equation to the source term.
  // This is now an iteration independent rule so that it can be used at both
  // {n} and {n,it} where it is needed.
  class TemporalToOmegaSourceTerm : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> omegaTemporalSourceTerm ;
      store<real> omegaSourceTerm ;
    public:

      // Define input and output.
      TemporalToOmegaSourceTerm() {
        name_store("omegaTemporalSourceTerm",omegaTemporalSourceTerm) ;
        name_store("omegaSourceTerm(X,Y)",omegaSourceTerm) ;
        input("omegaTemporalSourceTerm") ;
        output("omegaSourceTerm(X,Y)") ;
        constraint("geom_cells") ;
      }

      // Add for a single cell.
      void calculate(Entity cell) {
        omegaSourceTerm[cell]+=omegaTemporalSourceTerm[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToOmegaSourceTerm> registerTemporalToOmegaSourceTerm ;

  // Rule to compute the diagonal term for the linear system.
  class ComputeOmegaMatrixDiagonal : public pointwise_rule {
    private:
      const_store<real> omegaMainCoefficient ;
      store<real> D ;
    public:

      // Define input and output.
      ComputeOmegaMatrixDiagonal() {
        name_store("omegaMainCoefficient(massFluxCorrected_p,rhoStar)",omegaMainCoefficient) ;
        name_store("omegaStar_D",D) ;
        input("omegaMainCoefficient(massFluxCorrected_p,rhoStar)") ;
        output("omegaStar_D") ;
        constraint("geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) { D[cell]=omegaMainCoefficient[cell] ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeOmegaMatrixDiagonal> registerComputeOmegaMatrixDiagonal ;

  // Rule to compute the diagonal term for the linear system for cell next
  // to wall function boundaries.
  class ComputeOmegaMatrixDiagonalWallFunction : public pointwise_rule {
    private:
      store<real> D ;
    public:

      // Define input and output.
      ComputeOmegaMatrixDiagonalWallFunction() {
        name_store("wallFunction::omegaStar_D",D) ;
        output("wallFunction::omegaStar_D") ;
        constraint("wallFunctionCells") ;
      }

      // Set main coefficient to a large number.
      void calculate(Entity cell) { D[cell]=1.0e30 ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeOmegaMatrixDiagonalWallFunction>
    registerComputeOmegaMatrixDiagonalWallFunction ;

  // Rule to copy omegaStar_D for periodic faces.
  class ComputeOmegaMatrixDiagonalPeriodic : public pointwise_rule {
    private:
      const_Map cl,cr,pmap ;
      store<real> D ;
    public:

      // Define input and output.
      ComputeOmegaMatrixDiagonalPeriodic() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pmap",pmap) ;
        name_store("omegaStar_D",D) ;
        input("pmap->cl->omegaStar_D") ;
        output("cr->omegaStar_D") ;
      }

      // For a face.
      void calculate(Entity face) { D[cr[face]]=D[cl[pmap[face]]] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; } } ;

  register_rule<ComputeOmegaMatrixDiagonalPeriodic>
    registerComputeOmegaMatrixDiagonalPeriodic ;

  // Rule to initialize the lower terms for the linear system. This is now
  // a parametric variable with X=massFlux.
  class InitializeOmegaMatrixLower : public unit_rule {
    private:
      store<real> L ;
    public:

      // Define input and output.
      InitializeOmegaMatrixLower() {
        name_store("omegaStar_L(X)",L) ;
        output("omegaStar_L(X)") ;
        constraint("internalFaces") ;
      }

      // Initialize for a single face.
      void calculate(Entity face) { L[face]=0.0 ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeOmegaMatrixLower> registerInitializeOmegaMatrixLower ;

  // Rule to add the first-order inviscid flux contribution to the lower terms
  // for the linear system.
  class FOUInviscidFluxToOmegaMatrixLower : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      store<real> L ;
    public:

      // Define input and output.
      FOUInviscidFluxToOmegaMatrixLower() {
        name_store("X",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("omegaStar_L(X)",L) ;
        input("X,gridMassFlux") ;
        output("omegaStar_L(X)") ;
        constraint("internalFaces") ;
      }

      // Increment the lower term for a single face. Note that the increment
      // is the negative of the one in streamUns since this coefficient is
      // for a term on the lhs of the equation. In streamUns, the coefficient
      // is associated with a term on the rhs.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        if(netMassFlux>0.0) L[face]-=netMassFlux ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToOmegaMatrixLower>
    registerFOUInviscidFluxToOmegaMatrixLower ;

  // Rule to add the diffusive flux contribution to the lower terms for the
  // linear system.
  class DiffusiveFluxToOmegaMatrixLower : public apply_rule<store<real>,Loci::
  Summation<real> > {
    private:
      const_store<real> diffusionProduct ;
      const_store<real> faceRadius ;
      const_store<real> omegaViscosity ;
      store<real> L ;
    public:

      // Define input and output.
      DiffusiveFluxToOmegaMatrixLower() {
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("faceRadius",faceRadius) ;
        name_store("omegaViscosity",omegaViscosity) ;
        name_store("omegaStar_L(X)",L) ;
        input("diffusionProduct,faceRadius,omegaViscosity") ;
        output("omegaStar_L(X)") ;
        constraint("internalFaces,menterTurbulenceModel") ;
      }

      // Increment the lower term for a single face.
      void calculate(Entity face) {
        L[face]-=omegaViscosity[face]*diffusionProduct[face]*faceRadius[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToOmegaMatrixLower>
    registerDiffusiveFluxToOmegaMatrixLower ;

  // Rule to initialize the upper terms for the linear system. This is now
  // a parametric variable with X=massFlux.
  class InitializeOmegaMatrixUpper : public unit_rule {
    private:
      store<real> U ;
    public:

      // Define input and output.
      InitializeOmegaMatrixUpper() {
        name_store("omegaStar_U(X)",U) ;
        output("omegaStar_U(X)") ;
        constraint("internalFaces") ;
      }

      // Initialize for a single face.
      void calculate(Entity face) { U[face]=0.0 ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeOmegaMatrixUpper> registerInitializeOmegaMatrixUpper ;

  // Rule to add the first-order inviscid flux contribution to the upper terms
  // for the linear system.
  class FOUInviscidFluxToOmegaMatrixUpper : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      store<real> U ;
    public:

      // Define input and output.
      FOUInviscidFluxToOmegaMatrixUpper() {
        name_store("X",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("omegaStar_U(X)",U) ;
        input("X,gridMassFlux") ;
        output("omegaStar_U(X)") ;
        constraint("internalFaces") ;
      }

      // Increment the upper term for a single face. Contribution is opposite
      // to that in streamUns, as noted above for the lower terms.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        if(netMassFlux<0.0) U[face]+=netMassFlux ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToOmegaMatrixUpper>
    registerFOUInviscidFluxToOmegaMatrixUpper ;

  // Rule to add the diffusive flux contribution to the lower terms for the
  // linear system.
  class DiffusiveFluxToOmegaMatrixUpper : public apply_rule<store<real>,Loci::
  Summation<real> > {
    private:
      const_store<real> diffusionProduct ;
      const_store<real> faceRadius ;
      const_store<real> omegaViscosity ;
      store<real> U ;
    public:

      // Define input and output.
      DiffusiveFluxToOmegaMatrixUpper() {
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("faceRadius",faceRadius) ;
        name_store("omegaViscosity",omegaViscosity) ;
        name_store("omegaStar_U(X)",U) ;
        input("diffusionProduct,faceRadius,omegaViscosity") ;
        output("omegaStar_U(X)") ;
        constraint("internalFaces,menterTurbulenceModel") ;
      }

      // Increment the lower term for a single face.
      void calculate(Entity face) {
        U[face]-=omegaViscosity[face]*diffusionProduct[face]*faceRadius[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToOmegaMatrixUpper>
    registerDiffusiveFluxToOmegaMatrixUpper ;

  // Rule to compute the right-hand side for the linear system.
  class ComputeOmegaRHS : public pointwise_rule {
    private:
      const_store<real> omegaSourceTerm ;
      store<real> B ;
    public:

      // Define input and output.
      ComputeOmegaRHS() {
        name_store("omegaSourceTerm(vCorrected_p,massFluxCorrected_p)",omegaSourceTerm) ;
        name_store("omegaStar_B",B) ;
        input("omegaSourceTerm(vCorrected_p,massFluxCorrected_p)") ;
        output("omegaStar_B") ;
        constraint("geom_cells") ;
      }

      // Set source term for a single cell.
      void calculate(Entity cell) {
        B[cell]=omegaSourceTerm[cell] ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeOmegaRHS> registerComputeOmegaRHS ;

  // Rule to copy omegaStar_B for periodic faces.
  class ComputeOmegaRHSPeriodic : public pointwise_rule {
    private:
      const_Map cl,cr,pmap ;
      store<real> B ;
    public:

      // Define input and output.
      ComputeOmegaRHSPeriodic() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pmap",pmap) ;
        name_store("omegaStar_B",B) ;
        input("pmap->cl->omegaStar_B") ;
        output("cr->omegaStar_B") ;
      }

      // For a face.
      void calculate(Entity face) { B[cr[face]]=B[cl[pmap[face]]] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeOmegaRHSPeriodic> registerComputeOmegaRHSPeriodic ;

  // Temporary hack to get around possible bug in using maps to set priority
  // variable values.
  class OmegaCellHack : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> omegaWall ;
      store<real> omegaCell ;
    public:
                                                                                
      // Define input and output.
      OmegaCellHack() {
        name_store("ci",ci) ;
        name_store("omegaWall",omegaWall) ;
        name_store("omegaCell",omegaCell) ;
        input("omegaWall") ;
        output("ci->omegaCell") ;
        constraint("ref->wallFunction_BCoption") ;
      }
                                                                                
      // Compute omega for a single face.
      void calculate(Entity face) { omegaCell[ci[face]]=omegaWall[face] ; }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<OmegaCellHack> registerOmegaCellHack ;

  // Rule to compute the right-hand side for the linear system for cells next
  // to wall function boundaries.
  class ComputeOmegaRHSWallFunction : public pointwise_rule {
    private:
      const_store<real> omegaCell ;
      store<real> B ;
    public:

      // Define input and output.
      ComputeOmegaRHSWallFunction() {
        name_store("omegaCell",omegaCell) ;
        name_store("wallFunction::omegaStar_B",B) ;
        input("omegaCell") ;
        output("wallFunction::omegaStar_B") ;
        constraint("wallFunctionCells") ;
      }

      // Set right-hand side to large number times omegaCell.
      void calculate(Entity cell) { B[cell]=1.0e30*omegaCell[cell] ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeOmegaRHSWallFunction>
    registerComputeOmegaRHSWallFunction ;

//-----------------------------------------------------------------------------
// Rules for k and omega correctors.

  // Rule to compute the corrected value of k.
  class ComputeKCorrected : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_multiMap upper,lower ;
      const_store<real> k ;
      const_store<real> kMainCoefficient ;
      const_store<real> kStar_L,kStar_U ;
      const_store<real> kSourceTerm ;
      store<real> kCorrected ;
    public:

      // Define input and output.
      ComputeKCorrected() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("k",k) ;
        name_store("kMainCoefficient(massFluxCorrected_c,rhoCorrected)",kMainCoefficient) ;
        name_store("kStar_L(massFluxCorrected_c)",kStar_L) ;
        name_store("kStar_U(massFluxCorrected_c)",kStar_U) ;
        name_store("kSourceTerm(vCorrected_c,massFluxCorrected_c)",kSourceTerm) ;
        name_store("kCorrected",kCorrected) ;
        input("upper->cr->k,lower->cl->k") ;
        input("kMainCoefficient(massFluxCorrected_c,rhoCorrected)") ;
        input("lower->kStar_L(massFluxCorrected_c)") ;
        input("upper->kStar_U(massFluxCorrected_c)") ;
        input("kSourceTerm(vCorrected_c,massFluxCorrected_c)") ;
        output("kCorrected") ;
        constraint("geom_cells") ;
      }

      // Add explicit and implicit contributions for the cell, then divide by the
      // main coefficient.
      void calculate(Entity cell) {
        kCorrected[cell]=kSourceTerm[cell] ;
        for(const int *ui=upper.begin(cell);ui!=upper.end(cell);++ui)
          kCorrected[cell]-=kStar_U[*ui]*k[cr[*ui]] ;
        for(const int *li=lower.begin(cell);li!=lower.end(cell);++li)
          kCorrected[cell]-=kStar_L[*li]*k[cl[*li]] ;
        kCorrected[cell]/=kMainCoefficient[cell] ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeKCorrected> registerComputeKCorrected ;

  // Priority rule for kCorrected for wall-function cells.
  class ComputeKCorrectedWallFunction : public pointwise_rule {
    private:
      const_store<real> kCell ;
      store<real> kCorrected ;
    public:

      // Define input and output.
      ComputeKCorrectedWallFunction() {
        name_store("kCell",kCell) ;
        name_store("wallFunction::kCorrected",kCorrected) ;
        input("kCell") ;
        output("wallFunction::kCorrected") ;
        constraint("wallFunctionCells") ;
      }

      // Set value for a cell.
      void calculate(Entity cell) { kCorrected[cell]=kCell[cell] ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeKCorrectedWallFunction>
    registerComputeKCorrectedWallFunction ;

  // Rule to compute the corrected value of omega.
  class ComputeOmegaCorrected : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_multiMap upper,lower ;
      const_store<real> omega ;
      const_store<real> omegaMainCoefficient ;
      const_store<real> omegaStar_L,omegaStar_U ;
      const_store<real> omegaSourceTerm ;
      store<real> omegaCorrected ;
    public:

      // Define input and output.
      ComputeOmegaCorrected() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("omega",omega) ;
        name_store("omegaMainCoefficient(massFluxCorrected_c,rhoCorrected)",omegaMainCoefficient) ;
        name_store("omegaStar_L(massFluxCorrected_c)",omegaStar_L) ;
        name_store("omegaStar_U(massFluxCorrected_c)",omegaStar_U) ;
        name_store("omegaSourceTerm(vCorrected_c,massFluxCorrected_c)",omegaSourceTerm) ;
        name_store("omegaCorrected",omegaCorrected) ;
        input("upper->cr->omega,lower->cl->omega") ;
        input("omegaMainCoefficient(massFluxCorrected_c,rhoCorrected)") ;
        input("lower->omegaStar_L(massFluxCorrected_c)") ;
        input("upper->omegaStar_U(massFluxCorrected_c)") ;
        input("omegaSourceTerm(vCorrected_c,massFluxCorrected_c)") ;
        output("omegaCorrected") ;
        constraint("geom_cells") ;
      }

      // Add explicit and implicit contributions for the cell, then divide by the
      // main coefficient.
      void calculate(Entity cell) {
        omegaCorrected[cell]=omegaSourceTerm[cell] ;
        for(const int *ui=upper.begin(cell);ui!=upper.end(cell);++ui)
          omegaCorrected[cell]-=omegaStar_U[*ui]*omega[cr[*ui]] ;
        for(const int *li=lower.begin(cell);li!=lower.end(cell);++li)
          omegaCorrected[cell]-=omegaStar_L[*li]*omega[cl[*li]] ;
        omegaCorrected[cell]/=omegaMainCoefficient[cell] ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeOmegaCorrected> registerComputeOmegaCorrected ;

  // Priority rule for omegaCorrected for wall-function cells.
  class ComputeOmegaCorrectedWallFunction : public pointwise_rule {
    private:
      const_store<real> omegaCell ;
      store<real> omegaCorrected ;
    public:

      // Define input and output.
      ComputeOmegaCorrectedWallFunction() {
        name_store("omegaCell",omegaCell) ;
        name_store("wallFunction::omegaCorrected",omegaCorrected) ;
        input("omegaCell") ;
        output("wallFunction::omegaCorrected") ;
        constraint("wallFunctionCells") ;
      }

      // Set value for a cell.
      void calculate(Entity cell) { omegaCorrected[cell]=omegaCell[cell] ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeOmegaCorrectedWallFunction>
    registerComputeOmegaCorrectedWallFunction ;

//-----------------------------------------------------------------------------
// Rules for computing the residual of the omega equation.

  // Rule to initialize the omega residual.
  class InitializeOmegaResidual : public unit_rule {
    private:
      store<real> omegaResidual ;
    public:

      // Define input and output.
      InitializeOmegaResidual() {
        name_store("omegaResidual",omegaResidual) ;
        output("omegaResidual") ;
        constraint("kOmegaTurbulenceModel,vol") ;
      }

      // Initialize the residual for a single cell.
      void calculate(Entity cell) { omegaResidual[cell]=0.0 ; }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

//register_rule<InitializeOmegaResidual> registerInitializeOmegaResidual ;

  // Rule to compute the omega residual for each cell.
  class ComputeOmegaResidualOne : public apply_rule<store<real>,Loci::Summation
  <real> > {
    private:
      const_param<real> rhoScale,vScale,lScale,omegaScale ;
      const_store<real> D ;
      const_store<real> omega ;
      const_store<real> B ;
      store<real> omegaResidual ;
    private:
      real omegaFactor ;
    public:

      // Define input and output.
      ComputeOmegaResidualOne() {
        name_store("rhoScale",rhoScale) ;
        name_store("vScale",vScale) ;
        name_store("lScale",lScale) ;
        name_store("omegaScale",omegaScale) ;
        name_store("omegaStar_D",D) ;
        name_store("omega",omega) ;
        name_store("omegaStar_B",B) ;
        name_store("omegaResidual",omegaResidual) ;
        input("rhoScale,vScale,lScale,omegaScale") ;
        input("omegaStar_D,omega,omegaStar_B") ;
        output("omegaResidual") ;
        constraint("vol") ;
      }

      // Initialize the residual for a single cell.
      void calculate(Entity cell) {
        omegaResidual[cell]+=(D[cell]!=1.0e30)? (B[cell]-D[cell]*omega[cell])/
          omegaFactor:0.0 ;
      }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) {
        omegaFactor=(*rhoScale)*(*vScale)*(*omegaScale)*(*lScale)*(*lScale) ;
        do_loop(seq,this) ;
      }
  } ;

//register_rule<ComputeOmegaResidualOne> registerComputeOmegaResidualOne ;

  // Rule to compute the omega residual for each cell.
  class ComputeOmegaResidualTwo : public apply_rule<store<real>,Loci::Summation
  <real> > {
    private:
      const_Map cl,cr ;
      const_param<real> rhoScale,vScale,lScale,omegaScale ;
      const_store<real> omega ;
      const_store<real> D,L,U ;
      store<real> omegaResidual ;
    private:
      real omegaFactor ;
    public:

      // Define input and output.
      ComputeOmegaResidualTwo() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("rhoScale",rhoScale) ;
        name_store("vScale",vScale) ;
        name_store("lScale",lScale) ;
        name_store("omegaScale",omegaScale) ;
        name_store("omega",omega) ;
        name_store("omegaStar_D",D) ;
        name_store("omegaStar_L",L) ;
        name_store("omegaStar_U",U) ;
        name_store("omegaResidual",omegaResidual) ;
        input("rhoScale,vScale,lScale,omegaScale") ;
        input("(cl,cr)->(omega,omegaStar_D),omegaStar_L,omegaStar_U") ;
        output("(cl,cr)->omegaResidual") ;
        constraint("internalFaces") ;
      }

      // Add the neighbor contribution to the residual for each of the two
      // cells on either side of the face.
      void calculate(Entity face) {
        omegaResidual[cl[face]]-=(D[cl[face]]!=1.0e30)? U[face]*
          omega[cr[face]]/omegaFactor:0.0 ;
        omegaResidual[cr[face]]-=(D[cr[face]]!=1.0e30)? L[face]*
          omega[cl[face]]/omegaFactor:0.0 ;
      }

      // Add the neighbor contributions for a sequence of faces.
      virtual void compute(const sequence &seq) {
        omegaFactor=(*rhoScale)*(*vScale)*(*omegaScale)*(*lScale)*(*lScale) ;
        do_loop(seq,this) ;
      }
  } ;

//register_rule<ComputeOmegaResidualTwo> registerComputeOmegaResidualTwo ;

  // Rule to set total omega residual to zero for laminar flow.
  class TotalOmegaResidualLaminar : public singleton_rule {
    private:
      param<ScalarResidual> omegaResidualData ;
    public:

      // Define input and output.
      TotalOmegaResidualLaminar() {
        name_store("omegaResidualData",omegaResidualData) ;
        output("omegaResidualData") ;
        constraint("laminarFlow,geom_cells") ;
      }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) {
        omegaResidualData=ScalarResidual() ;
      }
  } ;

//register_rule<TotalOmegaResidualLaminar> registerTotalOmegaResidualLaminar ;

  // Rule to set total omega residual to zero for inviscid flow.
  class TotalOmegaResidualInviscid : public singleton_rule {
    private:
      param<ScalarResidual> omegaResidualData ;
    public:

      // Define input and output.
      TotalOmegaResidualInviscid() {
        name_store("omegaResidualData",omegaResidualData) ;
        output("omegaResidualData") ;
        constraint("inviscidFlow,geom_cells") ;
      }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) {
        omegaResidualData=ScalarResidual() ;
      }
  } ;

//register_rule<TotalOmegaResidualInviscid> registerTotalOmegaResidualInviscid ;

  // Rule to initialize the total omega residual.
  class InitializeTotalOmegaResidual : public unit_rule {
    private:
      param<ScalarResidual> omegaResidualData ;
    public:

      // Define input and output.
      InitializeTotalOmegaResidual() {
        name_store("omegaResidualData",omegaResidualData) ;
        output("omegaResidualData") ;
        constraint("turbulentFlow,geom_cells") ;
      }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) {
        omegaResidualData=ScalarResidual() ;
      }
  } ;

//register_rule<InitializeTotalOmegaResidual>
//  registerInitializeTotalOmegaResidual ;

  // Rule to compute the total omega residual. Note the use of omegaStar_D to
  // determine the wall-function cells to exclude from the sum.
  class ComputeTotalOmegaResidual : public apply_rule<param<ScalarResidual>,
  ScalarResidualJoin> {
    private:
      const_store<real> omegaResidual ;
      const_store<real> D ;
      const_store<vect3d> cellCenter ;
      param<ScalarResidual> omegaResidualData ;
    public:

      // Define input and output.
      ComputeTotalOmegaResidual() {
        name_store("omegaResidual",omegaResidual) ;
        name_store("omegaStar_D",D) ;
        name_store("cellcenter",cellCenter) ;
        name_store("omegaResidualData",omegaResidualData) ;
        input("omegaResidual,omegaStar_D,cellcenter") ;
        output("omegaResidualData") ;
        constraint("turbulentFlow,geom_cells") ;
      }

      // Add the cell contribution to the residual for a single cell.
      void calculate(Entity cell) {
        ScalarResidual temp ;
        temp.maxResidual=omegaResidual[cell] ;
        temp.totalResidual=abs(omegaResidual[cell]) ;
        temp.maxResidualLocation=cellCenter[cell] ;
        join(*omegaResidualData,temp) ;
      }

      // Add the cell contribution to the residual for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

//register_rule<ComputeTotalOmegaResidual> registerComputeTotalOmegaResidual ;

//-----------------------------------------------------------------------------
// Rules to limit k and omega.

  // This and the follwing rules compute the maximum value of omega assigned
  // on the noslip walls.
  class MaxOmegaWallUnit : public unit_rule {
    private:
      param<real> maxOmegaWall ;
    public:

      // Define input and output.
      MaxOmegaWallUnit() {
        name_store("maxOmegaWall",maxOmegaWall) ;
        output("maxOmegaWall") ;
        constraint("UNIVERSE") ;
      }

      // Initialize.
      virtual void compute(const sequence &seq) { *maxOmegaWall=0.0 ; }
  } ;

  register_rule<MaxOmegaWallUnit> registerMaxOmegaWallUnit ;

  // Find the maximum value. If we are using only wall functions, then the
  // maximum value will be zero which means no limiting will occur.
  class MaxOmegaWallApply : public apply_rule<param<real>,
  Loci::Maximum<real> > {
    private:
      const_store<real> noWallFunction ;
      const_store<real> omega_f ;
      param<real> maxOmegaWall ;
    public:

      // Define input and output.
      MaxOmegaWallApply() {
        name_store("noWallFunction",noWallFunction) ;
        name_store("omega_f",omega_f) ;
        name_store("maxOmegaWall",maxOmegaWall) ;
        input("noWallFunction,omega_f") ;
        output("maxOmegaWall") ;
        constraint("noslip_BC") ;
      }

      // Compare.
      void calculate(Entity face) {
        join(*maxOmegaWall,noWallFunction[face]*omega_f[face]) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<MaxOmegaWallApply> registerMaxOmegaWallApply ;

  // Limit k and omega to positive values using freestream turbulence
  // intensity. If either k or omega is less than zero, we set a
  // consistent freestream state using the current local density and
  // velocity and assuming that the eddy viscosity is 100 times the
  // current local laminar viscosity.
  class LimitK : public pointwise_rule {
    private:
      const_param<real> turbulenceIntensityFreestream ;
      const_store<real> rho,k,omega ;
      const_store<vect3d> v ;
      const_store<real> laminarViscosity ;
      store<real> kLimited ;
    public:

      // Define input and output.
      LimitK() {
        name_store("turbulenceIntensityFreestream",
          turbulenceIntensityFreestream) ;
        name_store("rho",rho) ;
        name_store("X",v) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("Y",k) ;
        name_store("Z",omega) ;
        name_store("kLimited(X,Y,Z)",kLimited) ;
        input("turbulenceIntensityFreestream") ;
        input("rho,X,laminarViscosity,Y,Z") ;
        output("kLimited(X,Y,Z)") ;
        constraint("geom_cells") ;
      }

      // Limit k for a single cell.
      void calculate(Entity cell) {
        if(k[cell]<0.0 || omega[cell]<=0.0){
          real vMag=norm(v[cell]) ;
          kLimited[cell]=1.5*pow((*turbulenceIntensityFreestream)*vMag,2) ;
        }else{
          kLimited[cell]=k[cell] ;
        }
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<LimitK> registerLimitK ;

  // Have to make this a separate rule now since we cannot have a single
  // rule outputting multiple parametric variables.
  class LimitOmega : public pointwise_rule {
    private:
      const_param<bool> debug ;
      const_store<vect3d> cellCenter ;
      const_store<real> rho,k,omega ;
      const_store<vect3d> v ;
      const_store<real> laminarViscosity ;
      const_store<real> kLimited ;
      store<real> omegaLimited ;
    public:

      // Define input and output.
      LimitOmega() {
        name_store("debug",debug) ;
        name_store("cellcenter",cellCenter) ;
        name_store("rho",rho) ;
        name_store("X",v) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("Y",k) ;
        name_store("Z",omega) ;
        name_store("kLimited(X,Y,Z)",kLimited) ;
        name_store("omegaLimited(X,Y,Z)",omegaLimited) ;
        input("debug") ;
        input("cellcenter,rho,X,laminarViscosity,Y,Z,kLimited(X,Y,Z)") ;
        output("omegaLimited(X,Y,Z)") ;
        constraint("geom_cells") ;
      }

      // Limit omega for a single cell. We have put a fix in for the
      // case where the cell velocity is zero, which was causing omega to be
      // assigned to zero, which later caused a divide by zero. JW 6/25/2009
      void calculate(Entity cell) {
        if(k[cell]<0.0 || omega[cell]<=0.0){
          if(*debug){
            cout << "WARNING: k or omega < 0 for cell " << cell << endl ;
            cout << "  cellcenter: " << cellCenter[cell] << endl ;
            cout << "  rho: " << rho[cell] << endl ;
            cout << "  v: " << v[cell] << endl ;
            cout << "  k: " << k[cell] << endl ;
            cout << "  omega: " << omega[cell] << endl ;
            cout << "  laminarViscosity: " << laminarViscosity[cell] << endl ;
          }
          real vMag=norm(v[cell]) ;
          if(vMag>1.0e-10){
            omegaLimited[cell]=rho[cell]*kLimited[cell]/(100.0*
              laminarViscosity[cell]) ;
            if(omegaLimited[cell]<1.0e-15) omegaLimited[cell]=1.0e-15 ;
          }else{
            omegaLimited[cell]=1.0 ;
          }
        }else if(omega[cell]>1.0e+15){
          if(*debug){
            cout << "WARNING: Clipping omega for cell " << cell << endl ;
            cout << "  cellcenter: " << cellCenter[cell] << endl ;
            cout << "  rho: " << rho[cell] << endl ;
            cout << "  v: " << v[cell] << endl ;
            cout << "  k: " << k[cell] << endl ;
            cout << "  omega: " << omega[cell] << endl ;
            cout << "  laminarViscosity: " << laminarViscosity[cell] << endl ;
          }
          omegaLimited[cell]=1.0e+15 ;
        }else{
          omegaLimited[cell]=omega[cell] ;
        }
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<LimitOmega> registerLimitOmega ;

//-----------------------------------------------------------------------------
// Rules for marching the k and omega equations. Note that these rules also
// are used by the other families of k-omega models.

  // Time build rule for k and omega when using BDF2 time integrator. Although
  // this rule sets k{n=-1} from k_ic and omega{n=-1} from omega_ic, these
  // values are not really used since BDF is used on the first timestep for
  // non-restarts.
  class TimeBuildKOmegaBDF2: public pointwise_rule {
    private:
      const_store<real> k_ic,omega_ic ;
      store<real> k,omega;
    public:

      // Define input and output.
      TimeBuildKOmegaBDF2() {
        name_store("k_ic",k_ic) ;
        name_store("omega_ic",omega_ic) ;
        name_store("k{n=-1}",k) ;
        name_store("omega{n=-1}",omega) ;
        input("k_ic,omega_ic") ;
        output("k{n=-1},omega{n=-1}") ;
        constraint("geom_cells,kOmegaTurbulenceModel") ;
      }

      // Assign k and omega for a single cell.
      void calculate(Entity cell) {
        k[cell]=k_ic[cell] ; omega[cell]=omega_ic[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildKOmegaBDF2> registerTimeBuildKOmegaBDF2 ;

  // Time build rule for k and omega.
  class TimeBuildKOmega: public pointwise_rule {
    private:
      const_store<real> k_ic,omega_ic ;
      store<real> kTimeStepZero,omegaTimeStepZero ;
    public:

      // Define input and output.
      TimeBuildKOmega() {
        name_store("k_ic",k_ic) ;
        name_store("omega_ic",omega_ic) ;
        name_store("k{n=0}",kTimeStepZero) ;
        name_store("omega{n=0}",omegaTimeStepZero) ;
        input("k_ic,omega_ic") ;
        output("k{n=0},omega{n=0}") ;
        constraint("geom_cells,kOmegaTurbulenceModel") ;
      }

      // Assign k and omega at time zero for a single cell.
      void calculate(Entity cell) {
        kTimeStepZero[cell]=k_ic[cell] ;
        omegaTimeStepZero[cell]=omega_ic[cell] ;
      }

      // Assign k and omega at time zero for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildKOmega> registerTimeBuildKOmega ;

  // Iteration build rule for k and omega.
  class IterationBuildKOmega : public pointwise_rule {
    private:
      const_store<real> kLimited,omegaLimited ;
      store<real> k,omega ;
    public:

      // Define input and output.
      IterationBuildKOmega() {
        name_store("kLimited(vCorrected_p,kStar,omegaStar){n}",kLimited) ;
        name_store("omegaLimited(vCorrected_p,kStar,omegaStar){n}",omegaLimited) ;
        name_store("k{n,it=0}",k) ;
        name_store("omega{n,it=0}",omega) ;
        input("kLimited(vCorrected_p,kStar,omegaStar){n}") ;
        input("omegaLimited(vCorrected_p,kStar,omegaStar){n}") ;
        output("k{n,it=0},omega{n,it=0}") ;
        constraint("geom_cells{n},kOmegaTurbulenceModel{n}") ;
      }

      // Assign k and omega for a single cell.
      void calculate(Entity cell) {
        k[cell]=kLimited[cell] ; omega[cell]=omegaLimited[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationBuildKOmega> registerIterationBuildKOmega ;

  // Iteration advance rule for k and omega. Added limited variables on
  // 07/09/2004. Used to have just kStar and omegaStar.
  class IterationAdvanceKOmega : public pointwise_rule {
    private:
      const_store<real> kLimited,omegaLimited ;
      store<real> k,omega ;
    public:

      // Define input and output.
      IterationAdvanceKOmega() {
        name_store("kLimited(vCorrected_c,kCorrected,omegaCorrected){n,it}",kLimited) ;
        name_store("omegaLimited(vCorrected_c,kCorrected,omegaCorrected){n,it}",omegaLimited) ;
        name_store("k{n,it+1}",k) ;
        name_store("omega{n,it+1}",omega) ;
        input("kLimited(vCorrected_c,kCorrected,omegaCorrected){n,it}") ;
        input("omegaLimited(vCorrected_c,kCorrected,omegaCorrected){n,it}") ;
        output("k{n,it+1},omega{n,it+1}") ;
        constraint("geom_cells{n,it},kOmegaTurbulenceModel{n,it}") ;
      }

      // Assign k and omega for a single cell.
      void calculate(Entity cell) {
        k[cell]=kLimited[cell] ; omega[cell]=omegaLimited[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationAdvanceKOmega> registerIterationAdvanceKOmega ;

  // Iteration collapse rule for k and omega. Note that we have changed
  // this rule from "k{n+1}=k{n,it}, ..." to the current form. With the old form
  // we could not do one iteration per time step, because k{n,it}, ... would
  // not get updated and thus the residuals would never change.
  class IterationCollapseKOmega : public pointwise_rule {
    private:
      const_param<bool> iterationFinished ;
      store<real> k,omega ;
    public:

      // Define input and output.
      IterationCollapseKOmega() {
        name_store("iterationFinished{n,it-1}",iterationFinished) ;
        name_store("k{n,it}",k) ;
        name_store("omega{n,it}",omega) ;
        input("iterationFinished{n,it-1}") ;
        input("k{n,it},omega{n,it}") ;
        output("k{n+1}=k{n,it},omega{n+1}=omega{n,it}") ;
        conditional("iterationFinished{n,it-1}") ;
        constraint("kOmegaTurbulenceModel{n,it},geom_cells{n,it}") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<IterationCollapseKOmega> registerIterationCollapseKOmega ;

}
