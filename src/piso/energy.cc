//-----------------------------------------------------------------------------
// Description: This file contains rules for the energy equation, including
//    rules for temperature and total enthalpy boundary conditions.
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

// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;

// StreamUns includes.
#include "const.h"
#include "residual.h"
#include "sciTypes.h"
#include "varsFileInputs.h"

namespace streamUns {

//-----------------------------------------------------------------------------
// Rules to process energy equation options from the .vars file.

  // Creates the energy equation solver constraints.
  class EnergyEquationSolverConstraints : public constraint_rule {
    private:
      const_param<EnergyEquationOptions> energyEquationOptions ;
      Constraint hSGSLinearSolver,hPETSCLinearSolver,hHYPRELinearSolver ;
    public:
                                                                                
      // Define input and output.
      EnergyEquationSolverConstraints() {
        name_store("energyEquationOptions",energyEquationOptions) ;
        name_store("hStar_SGSLinearSolver",hSGSLinearSolver) ;
        name_store("hStar_PETSCLinearSolver",hPETSCLinearSolver) ;
        name_store("hStar_HYPRELinearSolver",hHYPRELinearSolver) ;
        input("energyEquationOptions") ;
        output("hStar_SGSLinearSolver,hStar_PETSCLinearSolver") ;
        output("hStar_HYPRELinearSolver") ;
      }
                                                                                
      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        if((*energyEquationOptions).optionExists("linearSolver")){
          Loci::option_value_type optionValueType=energyEquationOptions->
            getOptionValueType("linearSolver") ;
          switch(optionValueType){
            case Loci::NAME:
              {
                Loci::option_values optionValues=energyEquationOptions->
                  getOption("linearSolver") ;
                string name ; optionValues.get_value(name) ;
                if(name=="SGS"){
                  hSGSLinearSolver=~EMPTY ; hPETSCLinearSolver=EMPTY ;
                  hHYPRELinearSolver=EMPTY ;
                }else if(name=="PETSC"){
                  hSGSLinearSolver=EMPTY ; hPETSCLinearSolver=~EMPTY ;
                  hHYPRELinearSolver=EMPTY ;
                }else if(name=="HYPRE"){
                  hSGSLinearSolver=EMPTY ; hPETSCLinearSolver=EMPTY ;
                  hHYPRELinearSolver=~EMPTY ;
                }else{
                  cerr << "Bad linearSolver for energyEquation." << endl ;
                  Loci::Abort() ;
                }
              }
              break ;
            default:
              cerr << "Bad type for linearSolver in energyEquation." << endl ;
              Loci::Abort() ;
          }
        }else{
          hSGSLinearSolver=~EMPTY ;
        }
      }

  } ;
                                                                                
  register_rule<EnergyEquationSolverConstraints>
    registerEnergyEquationSolverConstraints ;

  // Creates the energy equation solver parameters.
  class EnergyEquationSolverParameters : public singleton_rule {
    private:
      const_param<EnergyEquationOptions> energyEquationOptions ;
      param<int> hMaxIterations ;
    public:
                                                                                
      // Define input and output.
      EnergyEquationSolverParameters() {
        name_store("energyEquationOptions",energyEquationOptions) ;
        name_store("hStar_maxLinearSolverIterations",hMaxIterations) ;
        input("energyEquationOptions") ;
        output("hStar_maxLinearSolverIterations") ;
      }
                                                                                
      // Set up the parameters.
      virtual void compute(const sequence& seq) {
                                                                                
        // Maximum number of iterations.
        if((*energyEquationOptions).optionExists("maxIterations")){
          Loci::option_value_type optionValueType=energyEquationOptions->
            getOptionValueType("maxIterations") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                energyEquationOptions->getOption("maxIterations",temp) ;
                if(int(temp)<0){
                  cerr << "Bad maxIterations value for energyEquation."
                    << endl ; Loci::Abort() ;
                }
                *hMaxIterations=int(temp) ;
              }
              break ;
            default:
              cerr << "Bad type for maxIterations in energyEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *hMaxIterations=5 ;
        }
      }
  } ;
                                                                                
  register_rule<EnergyEquationSolverParameters>
    registerEnergyEquationSolverParamters ;

//-----------------------------------------------------------------------------
// Rules to create a constraint which indicates that the temperature is
// specified via some means.

  // Temperature specified with T=value.
  class BoundaryTemperatureConstantSpecification : public pointwise_rule {
    private:
      store<bool> specifiedTemperature_BC ;
    public:

      // Define input and output.
      BoundaryTemperatureConstantSpecification() {
        name_store("specifiedTemperature_BC",specifiedTemperature_BC) ;
        output("specifiedTemperature_BC") ;
        constraint("ref->T_BCoption") ;
      }

      // Set the flag to true.
      void calculate(Entity face) { specifiedTemperature_BC[face]=true ; }

      // Loop over all faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryTemperatureConstantSpecification>
    registerBoundaryTemperatureConstantSpecification ;

  // Temperature specified with cartesianT="bc_temp.dat".
  class BoundaryTemperatureCartesianSpecification : public pointwise_rule {
    private:
      store<bool> specifiedTemperature_BC ;
    public:

      // Define input and output.
      BoundaryTemperatureCartesianSpecification() {
        name_store("specifiedTemperature_BC",specifiedTemperature_BC) ;
        output("specifiedTemperature_BC") ;
        constraint("ref->cartesianT_BCoption") ;
      }

      // Set the flag to true.
      void calculate(Entity face) { specifiedTemperature_BC[face]=true ; }

      // Loop over all faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryTemperatureCartesianSpecification>
    registerBoundaryTemperatureCartesianSpecification ;

//-----------------------------------------------------------------------------
// Temperature boundary condition rules.

  // Rule for boundary faces with specified temperature. Assigns temperature
  // value to all boundary faces that have the property T_BC. Checked 04/20/04.
  class BoundaryTemperatureSpecification : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<real> T_BC ;
      store<real> temperature_f ;
    public:

      // Define input and output.
      BoundaryTemperatureSpecification() {
        name_store("ref",ref) ;
        name_store("T_BC",T_BC) ;
        name_store("temperature_f",temperature_f) ;
        input("ref->T_BC") ;
        output("temperature_f") ;
      }

      // Calculate temperature for a single face.
      void calculate(Entity face) {
        temperature_f[face]=T_BC[ref[face]] ;
      }

      // Calculate temperature for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryTemperatureSpecification>
    registerBoundaryTemperatureSpecification ;

  // Rule for extrapolating temperature to boundary faces. This occurs for
  // all outlets, slip and symmetry boundaries. Right now we are
  // using the low-order method. Checked 04/20/04.
  class BoundaryTemperatureExtrapolation : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> temperature ;
      store<real> temperature_f ;
    public:

      // Define input and output.
      BoundaryTemperatureExtrapolation() {
        name_store("ci",ci) ;
        name_store("temperature",temperature) ;
        name_store("temperature_f",temperature_f) ;
        input("ci->temperature") ;
        output("temperature_f") ;
        constraint("compressibleFlow,extrapolatedTemperature_BC") ;
      }

      // Calculate temperature for a single face.
      void calculate(Entity face) {
        temperature_f[face]=temperature[ci[face]] ;
      }

      // Calculate temperature for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryTemperatureExtrapolation>
    registerBoundaryTemperatureExtrapolation ;

  // Rule for adiabatic faces. Assigns temperature value to all boundary faces
  // that have the property adiabatic_BC.
  class BoundaryTemperatureAdiabatic : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> temperature ;
      store<real> temperature_f ;
    public:

      // Define input and output.
      BoundaryTemperatureAdiabatic() {
        name_store("ci",ci) ;
        name_store("temperature",temperature) ;
        name_store("temperature_f",temperature_f) ;
        input("ci->temperature") ;
        output("temperature_f") ;
        constraint("compressibleFlow,ref->adiabatic_BCoption") ;
      }

      // Calculate temperature for a single face.
      void calculate(Entity face) {
        temperature_f[face]=temperature[ci[face]] ;
      }

      // Calculate temperature for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryTemperatureAdiabatic>
    registerBoundaryTemperatureAdiabatic ;

  // Rule for specifying temperature with a profile. A single Cartesian
  // coordinate is used for the interpolation.
  class BoundaryTemperatureProfileCartesian : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<vect3d> faceCenter ;
      const_store<string> cartesianT_BC ;
      store<real> temperature_f ;
    public:

      // Define input and output.
      BoundaryTemperatureProfileCartesian() {
        name_store("ref",ref) ;
        name_store("facecenter",faceCenter) ;
        name_store("cartesianT_BC",cartesianT_BC) ;
        name_store("temperature_f",temperature_f) ;
        input("facecenter,ref->cartesianT_BC") ;
        output("temperature_f") ;
        disable_threading() ;
      }

      // Calculate temperature for all faces in sequence.
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
            string fileName=cartesianT_BC[bci->first] ;
            ifstream in(fileName.c_str(),ios::in) ;
            if(in.fail()) {
              cerr << "Open failed on " << fileName.c_str() << endl ;
              Loci::Abort() ;
            }

            // Skip spaces.
            while(!in.eof() && isspace(in.peek())) in.get() ;

            // Read in the number of points on the boundary.
            int np,coordFlag ; in >> np >> coordFlag ;
            if(np<2){
              cerr << "Bad number of data points in bc_temp.dat." << endl ;
              Loci::Abort() ;
            }
            if(coordFlag<0 || coordFlag>2){
              cerr << "Bad coordinate flag in bc_temp.dat." << endl ;
              Loci::Abort() ;
            }

            // Read in the coordinates and temperatures.
            vect3d *center=new vect3d[np] ; real *T =new real[np] ;
            for(int i=0;i<np;++i){
              in >> center[i] ; in >> T[i] ;
              if(T[i]<=0.0){
                cerr << "Negative or zero temperature for data point " << i
                  << " in bc_temp.dat." << endl ; Loci::Abort() ;
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
                    temperature_f[face]=T[np-1] ; continue ;
                  }
                else
                  if(id1>0){
                    id2=id1-1 ;
                  }else{
                    temperature_f[face]=T[0] ; continue ;
                  }
              else {
                temperature_f[face]=T[id1] ; continue ;
              }
              real d2=0.0 ;
              switch(coordFlag){
                case 0: d2=fabs(center[id2].x-fcenter.x) ; break ;
                case 1: d2=fabs(center[id2].y-fcenter.y) ; break ;
                case 2: d2=fabs(center[id2].z-fcenter.z) ; break ;
              }
              temperature_f[face]=(d2*T[id1]+d1*T[id2])/(d1+d2) ;
            }
            delete [] center ; delete [] T ;
          }
        }
      }
  } ;

  register_rule<BoundaryTemperatureProfileCartesian>
    registerBoundaryTemperatureProfileCartesian ;

//-----------------------------------------------------------------------------
// Total enthalpy boundary condition rules.

  // Rule for extrapolating total enthalpy to symmetry faces.
  class BoundaryTotalEnthalpySymmetry : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> h ;
      store<real> h_f ;
    public:

      // Define input and output.
      BoundaryTotalEnthalpySymmetry() {
        name_store("ci",ci) ;
        name_store("h",h) ;
        name_store("symmetry::h_f",h_f) ;
        input("ci->h") ;
        output("symmetry::h_f") ;
        constraint("compressibleFlow,symmetry_BC") ;
      }

      // Calculate total enthalpy for a single face.
      void calculate(Entity face) { h_f[face]=h[ci[face]] ; }

      // Calculate total enthalpy for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryTotalEnthalpySymmetry>
    registerBoundaryTotalEnthalpySymmetry ;

  // Rule for computing total enthalpy on boundary faces.
  class BoundaryTotalEnthalpyComputation : public pointwise_rule {
    private:
      const_store<vect3d> v_f ;
      const_store<EOS::State> eosState_f ;
      store<real> h_f ;
    public:

      // Define input and output.
      BoundaryTotalEnthalpyComputation() {
        name_store("v_f",v_f) ;
        name_store("eos_state_f",eosState_f) ;
        name_store("h_f",h_f) ;
        input("v_f,eos_state_f") ;
        output("h_f") ;
      }

      // Calculate total enthalpy for a single face.
      void calculate(Entity face) {
        h_f[face]=eosState_f[face].enthalpy()+0.5*dot(v_f[face],v_f[face]) ;
      }

      // Calculate total enthalpy for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryTotalEnthalpyComputation>
    registerBoundaryTotalEnthalpyComputation ;

//-----------------------------------------------------------------------------
// Rules to provide a zero temperature for incompressible flows. Simplifies the
// coding in several other places at the expense of a single variable.

  // Temperature for cells.
  class IncompressibleTemperatureInterior : public pointwise_rule {
    private:
      store<real> temperature ;
    public:

      // Define input and output.
      IncompressibleTemperatureInterior() {
        name_store("temperature",temperature) ;
        output("temperature") ;
        constraint("incompressibleFlow,geom_cells") ;
      }

      // Set the temperature to zero for a single cell.
      void calculate(Entity cell) { temperature[cell]=0.0 ; }

      // Set the temperature to zero for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IncompressibleTemperatureInterior>
    registerIncompressibleTemperatureInterior ;

  // Temperature for boundary faces.
  class IncompressibleTemperatureBoundary : public pointwise_rule {
    private:
      store<real> temperature_f ;
    public:

      // Define input and output.
      IncompressibleTemperatureBoundary() {
        name_store("temperature_f",temperature_f) ;
        output("temperature_f") ;
        constraint("incompressibleFlow,boundaryFaces") ;
      }

      // Set the temperature to zero for a single face.
      void calculate(Entity face) { temperature_f[face]=0.0 ; }

      // Set the temperature to zero for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IncompressibleTemperatureBoundary>
    registerIncompressibleTemperatureBoundary ;

//-----------------------------------------------------------------------------
// Rules for computing the diffusion coefficient for temperature.

  // Temperature diffusion coefficient. This has been made parametric so we
  // can use it for both cells and boundary faces.
  class TemperatureDiffusionCoefficient : public pointwise_rule {
    private:
      const_param<real> turbulentPrandtlNumber ;
      const_store<real> T,p,eddyViscosity ;
      const_storeVec<real> y ;
      const_store<real> cp,thermalConductivity ;
      store<real> tDiffusionCoeff ;
    public:

      // Define input and output.
      TemperatureDiffusionCoefficient() {
        name_store("turbulentPrandtlNumber",turbulentPrandtlNumber) ;
        name_store("T",T) ;
        name_store("P",p) ;
        name_store("Y",y) ;
        name_store("V",eddyViscosity) ;
        name_store("cp",cp) ;
        name_store("kconduct(T,P,Y)",thermalConductivity) ;
        name_store("tDiffusionCoeff(T,P,Y,V)",tDiffusionCoeff) ;
        input("turbulentPrandtlNumber,T,P,Y,V,cp,kconduct(T,P,Y)") ;
        output("tDiffusionCoeff(T,P,Y,V)") ;
      }

      // Compute for a single entity.
      void calculate(Entity e) {
        tDiffusionCoeff[e]=thermalConductivity[e]+eddyViscosity[e]*cp[e]/
          (*turbulentPrandtlNumber) ;
      }

      // Compute for a sequence of entities.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemperatureDiffusionCoefficient>
    registerTemperatureDiffusionCoefficient ;

//-----------------------------------------------------------------------------
// Rules for computing the kinetic energy which is use in the diffusion term
// for total enthalpy. Cannot make this parametric to handle both cells and
// boundary faces since we must have separate name for computing the gradient.

  // Kinetic energy per unit mass for cells.
  class ComputeKineticEnergyInterior : public pointwise_rule {
    private:
      const_store<vect3d> v ;
      store<real> kineticEnergy ;
    public:

      // Define input and output.
      ComputeKineticEnergyInterior() {
        name_store("v",v) ;
        name_store("kineticEnergy",kineticEnergy) ;
        input("v") ;
        output("kineticEnergy") ;
        constraint("geom_cells") ;
      }

      // Compute for a single cell.
      void calculate(Entity cell) {
        kineticEnergy[cell]=0.5*dot(v[cell],v[cell]) ;
      }

      // Compute for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeKineticEnergyInterior>
    registerComputeKineticEnergyInterior ;

  // Kinetic energy per unit mass for boundary faces.
  class ComputeKineticEnergyBoundary : public pointwise_rule {
    private:
      const_store<vect3d> v_f ;
      store<real> kineticEnergy_f ;
    public:

      // Define input and output.
      ComputeKineticEnergyBoundary() {
        name_store("v_f",v_f) ;
        name_store("kineticEnergy_f",kineticEnergy_f) ;
        input("v_f") ;
        output("kineticEnergy_f") ;
        constraint("boundaryFaces") ;
      }

      // Compute for a single face.
      void calculate(Entity face) {
        kineticEnergy_f[face]=0.5*dot(v_f[face],v_f[face]) ;
      }

      // Compute for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeKineticEnergyBoundary>
    registerComputeKineticEnergyBoundary ;

//-----------------------------------------------------------------------------
// Rules to create a constraint for boundary faces with non-zero energy
// diffusion flux.

  // All inlet faces have non-zero diffusion flux.
  class BoundaryEnergyDiffusionInlet : public pointwise_rule {
    private:
      store<bool> boundaryEnergyDiffusion ;
    public:

      // Define input and output.
      BoundaryEnergyDiffusionInlet() {
        name_store("boundaryEnergyDiffusion",boundaryEnergyDiffusion) ;
        output("boundaryEnergyDiffusion") ;
        constraint("inlet_BC") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { boundaryEnergyDiffusion[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryEnergyDiffusionInlet>
    registerBoundaryEnergyDiffusionInlet ;

  // All outlet faces have non-zero diffusion flux.
  class BoundaryEnergyDiffusionOutlet : public pointwise_rule {
    private:
      store<bool> boundaryEnergyDiffusion ;
    public:

      // Define input and output.
      BoundaryEnergyDiffusionOutlet() {
        name_store("boundaryEnergyDiffusion",boundaryEnergyDiffusion) ;
        output("boundaryEnergyDiffusion") ;
        constraint("outlet_BC") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { boundaryEnergyDiffusion[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryEnergyDiffusionOutlet>
    registerBoundaryEnergyDiffusionOutlet ;

  // No-slip faces with specified temperature have non-zero diffusion flux.
  class BoundaryEnergyDiffusionNoSlip : public pointwise_rule {
    private:
      store<bool> boundaryEnergyDiffusion ;
    public:

      // Define input and output.
      BoundaryEnergyDiffusionNoSlip() {
        name_store("boundaryEnergyDiffusion",boundaryEnergyDiffusion) ;
        output("boundaryEnergyDiffusion") ;
        constraint("noslip_BC,ref->T_BCoption") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { boundaryEnergyDiffusion[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryEnergyDiffusionNoSlip>
    registerBoundaryEnergyDiffusionNoSlip ;

  // No-slip faces with prescibed temperature profile have non-zero diffusion
  // flux.
  class BoundaryEnergyDiffusionNoSlipPrescribed : public pointwise_rule {
    private:
      store<bool> boundaryEnergyDiffusion ;
    public:

      // Define input and output.
      BoundaryEnergyDiffusionNoSlipPrescribed() {
        name_store("boundaryEnergyDiffusion",boundaryEnergyDiffusion) ;
        output("boundaryEnergyDiffusion") ;
        constraint("noslip_BC,ref->cartesianT_BCoption") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { boundaryEnergyDiffusion[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryEnergyDiffusionNoSlipPrescribed>
    registerBoundaryEnergyDiffusionNoSlipPrescribed ;

  // Interface faces have non-zero diffusion flux.
  class BoundaryEnergyDiffusionInterface : public pointwise_rule {
    private:
      store<bool> boundaryEnergyDiffusion ;
    public:

      // Define input and output.
      BoundaryEnergyDiffusionInterface() {
        name_store("boundaryEnergyDiffusion",boundaryEnergyDiffusion) ;
        output("boundaryEnergyDiffusion") ;
        constraint("interface_BC") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { boundaryEnergyDiffusion[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryEnergyDiffusionInterface>
    registerBoundaryEnergyDiffusionInterface ;

//-----------------------------------------------------------------------------
// Rules for computing the heat flux on walls.

  // Default rule computes heat flux assuming no wall function. This rule is
  // ultimately meant to provide values for noslip walls where there is no
  // wall function.
  class QWallDefault : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> temperature ;
      const_store<real> thermalConductivity ;
      const_store<real> temperature_f ;
      const_store<vect3d> faceCenter,cellCenter ;
      const_store<Area> area ;
      store<real> qWall ;
    public:
                                                                                
      // Define input and output.
      QWallDefault() {
        name_store("ci",ci) ;
        name_store("temperature",temperature) ;
        name_store("kconduct(temperature_f,p_f,y_f)",
          thermalConductivity) ;
        name_store("temperature_f",temperature_f) ;
        name_store("facecenter",faceCenter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("area",area) ;
        name_store("qWall",qWall) ;
        input("ci->(temperature,cellcenter)") ;
        input("kconduct(temperature_f,p_f,y_f)") ;
        input("temperature_f,facecenter,area") ;
        output("qWall") ;
        constraint("boundaryFaces") ;
      }
                                                                                
      // Heat flux for the face.
      void calculate (Entity face) {
        const real yNormal=dot(faceCenter[face]-cellCenter[ci[face]],
          area[face].n) ;
        qWall[face]=thermalConductivity[face]*(temperature[ci[face]]-
          temperature_f[face])/yNormal ;
      }
                                                                                
      // Loop over faces.
      virtual void compute (const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<QWallDefault> registerQWallDefault ;

  // Rule to compute the heat flux at inlet boundaries.
  class QWallInlet : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<vect3d> temperatureGradient ;
      const_store<real> tDiffusionCoeff_f ;
      const_store<Area> area ;
      store<real> qWall ;
    public:
                                                                                
      // Define input and output.
      QWallInlet() {
        name_store("ci",ci) ;
        name_store("grads(temperature)",temperatureGradient) ;
        name_store("tDiffusionCoeff(temperature_f,p_f,y_f,eddyViscosity_f)",
          tDiffusionCoeff_f) ;
        name_store("area",area) ;
        name_store("inlet::qWall",qWall) ;
        input("ci->grads(temperature),area") ;
        input("tDiffusionCoeff(temperature_f,p_f,y_f,eddyViscosity_f)") ;
        output("inlet::qWall") ;
        constraint("viscousFlow,inlet_BC") ;
      }
                                                                                
      // Compute the heat flux for a face.
      void calculate (Entity face) {
        qWall[face]=-tDiffusionCoeff_f[face]*
          dot(temperatureGradient[ci[face]],area[face].n) ;
      }
                                                                                
      // Loop over faces
      virtual void compute (const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<QWallInlet> registerQWallInlet ;

  // Rule to compute the heat flux at outlet boundaries.
  class QWallOutlet : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<vect3d> temperatureGradient ;
      const_store<real> tDiffusionCoeff_f ;
      const_store<Area> area ;
      store<real> qWall ;
    public:
                                                                                
      // Define input and output.
      QWallOutlet() {
        name_store("ci",ci) ;
        name_store("grads(temperature)",temperatureGradient) ;
        name_store("tDiffusionCoeff(temperature_f,p_f,y_f,eddyViscosity_f)",
          tDiffusionCoeff_f) ;
        name_store("area",area) ;
        name_store("outlet::qWall",qWall) ;
        input("ci->grads(temperature),area") ;
        input("tDiffusionCoeff(temperature_f,p_f,y_f,eddyViscosity_f)") ;
        output("outlet::qWall") ;
        constraint("viscousFlow,outlet_BC") ;
      }
                                                                                
      // Compute the heat flux for a face.
      void calculate (Entity face) {
        qWall[face]=-tDiffusionCoeff_f[face]*
          dot(temperatureGradient[ci[face]],area[face].n) ;
      }
                                                                                
      // Loop over faces
      virtual void compute (const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<QWallOutlet> registerQWallOutlet ;

  // Priority rule for heat flux at a symmetry boundary.
  class QWallSymmetry : public pointwise_rule {
    private:
      store<real> qWall ;
    public:
                                                                                
      // Define input and output.
      QWallSymmetry() {
        name_store("symmetry::qWall",qWall) ;
        output("symmetry::qWall") ;
        constraint("viscousFlow,symmetry_BC") ;
      }
                                                                                
      // No heat flux.
      void calculate (Entity face) { qWall[face]=0.0 ; }
                                                                                
      // Loop over faces.
      virtual void compute (const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<QWallSymmetry> registerQWallSymmetry ;

  // Priority rule for heat flux at a slip boundary.
  class QWallSlip : public pointwise_rule {
    private:
      store<real> qWall ;
    public:
                                                                                
      // Define input and output.
      QWallSlip() {
        name_store("slip::qWall",qWall) ;
        output("slip::qWall") ;
        constraint("viscousFlow,slip_BC") ;
      }
                                                                                
      // No heat flux.
      void calculate (Entity face) { qWall[face]=0.0 ; }
                                                                                
      // Loop over faces.
      virtual void compute (const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<QWallSlip> registerQWallSlip ;

  // Priority rule for adiabatic boundaries.
  class QWallAdiabatic : public pointwise_rule {
    private:
      store<real> qWall ;
    public:
                                                                                
      // Define input and output.
      QWallAdiabatic() {
        name_store("adiabatic::qWall",qWall) ;
        output("adiabatic::qWall") ;
        constraint("viscousFlow,ref->adiabatic_BCoption") ;
      }
                                                                                
      // No heat flux.
      void calculate (Entity face) { qWall[face]=0.0 ; }
                                                                                
      // Loop over faces.
      virtual void compute (const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<QWallAdiabatic> registerQWallAdiabatic ;

  // Priority rule for wall function boundaries.
  class QWallWallFunction : public pointwise_rule {
    private:
      const_store<real> qWallTemp ;
      store<real> qWall ;
    public:
                                                                                
      // Define input and output.
      QWallWallFunction() {
        name_store("qWallTemp",qWallTemp) ;
        name_store("wallFunction::qWall",qWall) ;
        input("qWallTemp") ;
        output("wallFunction::qWall") ;
        constraint("viscousFlow,ref->wallFunction_BCoption") ;
        constraint("specifiedTemperature_BC") ;
      }
                                                                                
      // Copy temp value from wall function rule.
      void calculate (Entity face) { qWall[face]=qWallTemp[face] ; }
                                                                                
      // Loop over faces.
      virtual void compute (const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<QWallWallFunction> registerQWallWallFunction ;

//-----------------------------------------------------------------------------
// Scheme independent rules for assembling the energy equation.

  // Rule to initialize the main coefficient. This is now a parametric variable
  // with X representing the mass flux and Y representing density.
  class InitializeTotalEnthalpyMainCoefficient : public unit_rule {
    private:
      store<real> hMainCoefficient ;
    public:

      // Define input and output.
      InitializeTotalEnthalpyMainCoefficient() {
        name_store("hMainCoefficient(X,Y)",hMainCoefficient) ;
        output("hMainCoefficient(X,Y)") ;
        constraint("vol") ;
      }

      // Set the main coefficient to zero for a single cell.
      void calculate(Entity cell) { hMainCoefficient[cell]=0.0 ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeTotalEnthalpyMainCoefficient>
    registerInitializeTotalEnthalpyMainCoefficient ;

  // Rule to add the first-order inviscid flux contribution to the main
  // coefficient for interior faces.
  class FOUInviscidFluxToTotalEnthalpyMainCoefficientInterior : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> massFlux ;
      store<real> hMainCoefficient ;
    public:

      // Define input and output.
      FOUInviscidFluxToTotalEnthalpyMainCoefficientInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("X",massFlux) ;
        name_store("hMainCoefficient(X,Y)",hMainCoefficient) ;
        input("X") ;
        output("(cl,cr)->hMainCoefficient(X,Y)") ;
        constraint("internalFaces") ;
      }

      // Increment the main coefficient for cells attached to a single face.
      void calculate(Entity face) {
        if(massFlux[face]>0.0){
          hMainCoefficient[cl[face]]+=massFlux[face] ;
        }else{
          hMainCoefficient[cr[face]]-=massFlux[face] ;
        }
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToTotalEnthalpyMainCoefficientInterior>
    registerFOUInviscidFluxToTotalEnthalpyMainCoefficientInterior ;

  // Rule to convert matrix form from 'natural' to 'diagonalDominance'.
  class DiagonalDominanceToTotalEnthalpyMainCoefficientInterior : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> massFlux ;
      store<real> hMainCoefficient ;
    public:

      // Define input and output.
      DiagonalDominanceToTotalEnthalpyMainCoefficientInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("X",massFlux) ;
        name_store("hMainCoefficient(X,Y)",hMainCoefficient) ;
        input("X") ;
        output("(cl,cr)->hMainCoefficient(X,Y)") ;
        constraint("internalFaces,diagonalDominance") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        hMainCoefficient[cl[face]]-=massFlux[face] ;
        hMainCoefficient[cr[face]]+=massFlux[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiagonalDominanceToTotalEnthalpyMainCoefficientInterior>
    registerDiagonalDominanceToTotalEnthalpyMainCoefficientInterior ;

  // Rule to add the first-order inviscid flux contribution to the main
  // coefficient for boundary faces. NOTE: Using variable massFlux results
  // in using old mass flux for specified pressure boundaries instead of
  // newly corrected boundary mass flux. Will figure how to get around
  // this later.
  class FOUInviscidFluxToTotalEnthalpyMainCoefficientBoundary : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> massFlux ;
      store<real> hMainCoefficient ;
    public:

      // Define input and output.
      FOUInviscidFluxToTotalEnthalpyMainCoefficientBoundary() {
        name_store("ci",ci) ;
        name_store("X",massFlux) ;
        name_store("hMainCoefficient(X,Y)",hMainCoefficient) ;
        input("X") ;
        output("ci->hMainCoefficient(X,Y)") ;
        constraint("boundaryFaces") ;
      }

      // Increment the main coefficient for the cell attach to a single face.
      void calculate(Entity face) {
        if(massFlux[face]>0.0) hMainCoefficient[ci[face]]+=massFlux[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToTotalEnthalpyMainCoefficientBoundary>
    registerFOUInviscidFluxToTotalEnthalpyMainCoefficientBoundary ;

  // Rule to add temporal component of the total enthalpy equation to the main
  // coefficient.
  class TemporalToTotalEnthalpyMainCoefficient: public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> timeStep ;
      const_param<real> timeIntegratorFactor0 ;
      const_store<real> timeStepFactor ;
      const_store<real> rho ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> hMainCoefficient ;
    public:

      // Define input and output.
      TemporalToTotalEnthalpyMainCoefficient() {
        name_store("timeStep",timeStep) ;
        name_store("timeIntegratorFactor0",timeIntegratorFactor0) ;
        name_store("timeStepFactor",timeStepFactor) ;
        name_store("Y",rho) ;
        name_store("vol",vol) ;
        name_store("cellRadius",cellRadius) ;
        name_store("hMainCoefficient(X,Y)",hMainCoefficient) ;
        input("Y,vol,cellRadius,timeStep") ;
        input("timeStepFactor,timeIntegratorFactor0") ;
        output("hMainCoefficient(X,Y)") ;
        constraint("geom_cells") ;
      }

      // Add temporal component for a single cell.
      void calculate(Entity cell) {
        hMainCoefficient[cell]+=0.5*rho[cell]*vol[cell]*cellRadius[cell]*
          ((*timeIntegratorFactor0)+1.0)/((*timeStep)*timeStepFactor[cell]) ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToTotalEnthalpyMainCoefficient>
    registerTemporalToTotalEnthalpyMainCoefficient ;

  // Rule to initialize the source term. This is now a parametric variable
  // with X=velocity, Y=pressure and Z=massFlux.
  class InitializeTotalEnthalpySourceTerm : public unit_rule {
    private:
      store<real> hSourceTerm ;
    public:

      // Define input and output.
      InitializeTotalEnthalpySourceTerm() {
        name_store("hSourceTerm(X,Y,Z)",hSourceTerm) ;
        output("hSourceTerm(X,Y,Z)") ;
        constraint("vol") ;
      }

      // Set the source term to zero for a single cell.
      void calculate(Entity cell) { hSourceTerm[cell]=0.0 ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeTotalEnthalpySourceTerm>
    registerInitializeTotalEnthalpySourceTerm ;

  // Rule to add the net mass flux times the cell value to offset the
  // contribution to the lhs that was added to get a diagonally
  // dominant matrix.
  class DiagonalDominanceToTotalEnthalpySourceTerm : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> h ;
      const_store<real> netMassFlux ;
      store<real> hSourceTerm ;
    public:

      // Define input and output.
      DiagonalDominanceToTotalEnthalpySourceTerm() {
        name_store("h",h) ;
        name_store("netMassFlux(Z)",netMassFlux) ;
        name_store("hSourceTerm(X,Y,Z)",hSourceTerm) ;
        input("h,netMassFlux(Z)") ;
        output("hSourceTerm(X,Y,Z)") ;
        constraint("geom_cells,diagonalDominance") ;
      }

      // Increment the source term for the cell.
      void calculate(Entity cell) {
        hSourceTerm[cell]-=netMassFlux[cell]*h[cell] ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiagonalDominanceToTotalEnthalpySourceTerm>
    registerDiagonalDominanceToTotalEnthalpySourceTerm ;

  // Rule to add the first-order inviscid flux contribution to the source term
  // for boundary faces.
  class FOUInviscidFluxToTotalEnthalpySourceTermBoundary : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> massFlux ;
      const_store<real> h_f ;
      store<real> hSourceTerm ;
    public:

      // Define input and output.
      FOUInviscidFluxToTotalEnthalpySourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("Z",massFlux) ;
        name_store("h_f",h_f) ;
        name_store("hSourceTerm(X,Y,Z)",hSourceTerm) ;
        input("Z,h_f") ;
        output("ci->hSourceTerm(X,Y,Z)") ;
        constraint("boundaryFaces") ;
      }

      // Increment the source term for the cell attach to a single face.
      void calculate(Entity face) {
        if(massFlux[face]<0.0) hSourceTerm[ci[face]]-=massFlux[face]*
          h_f[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToTotalEnthalpySourceTermBoundary>
    registerFOUInviscidFluxToTotalEnthalpySourceTermBoundary ;

  // Rename rule to allow us to override the total enthalpy gradient for
  // turbomachinery flows. This rule assigns grads(h) while a priority rule
  // in the turbomachinery module assigns gradsTurbo(h).
  class HGradient : public pointwise_rule {
    private:
      store<vect3d> hGradient ;
    public:

      // Define input and output.
      HGradient() {
        name_store("grads(h)",hGradient) ;
        input("grads(h)") ;
        output("hGradient=grads(h)") ;
        constraint("vol") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}

  } ;

  register_rule<HGradient> registerHGradient ;

  // Rule to add the second-order convection contribution to the source term
  // for interior faces. Checked.
  class SOUInviscidFluxToTotalEnthalpySourceTermInterior : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<vect3d> hGradient ;
      const_store<real> hLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      store<real> hSourceTerm ;
    public:

      // Define input and output.
      SOUInviscidFluxToTotalEnthalpySourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("hGradient",hGradient) ;
        name_store("limiters(h)",hLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("Z",massFlux) ;
        name_store("hSourceTerm(X,Y,Z)",hSourceTerm) ;
        input("(cl,cr)->(hGradient,limiters(h),cellcenter)") ;
        input("facecenter,Z") ;
        output("(cl,cr)->hSourceTerm(X,Y,Z)") ;
        constraint("internalFaces,souOrRoeInviscidFlux") ;
      }

      // Increment the source term for the cells attached to a single face.
      void calculate(Entity face) {
        real secondOrderSource=(massFlux[face]>0.0)? massFlux[face]*
          hLimiter[cl[face]]*dot(hGradient[cl[face]],faceCenter[face]-
          cellCenter[cl[face]]):massFlux[face]*hLimiter[cr[face]]*
          dot(hGradient[cr[face]],faceCenter[face]-cellCenter[cr[face]]) ;
        hSourceTerm[cl[face]]-=secondOrderSource ;
        hSourceTerm[cr[face]]+=secondOrderSource ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SOUInviscidFluxToTotalEnthalpySourceTermInterior>
    registerSOUInviscidFluxToTotalEnthalpySourceTermInterior ;

  // Rule to add the second-order convection contribution to the source term
  // for boundary faces. Checked.
  class SOUInviscidFluxToTotalEnthalpySourceTermBoundary : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<vect3d> hGradient ;
      const_store<real> hLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      const_store<real> h_f ;
      store<real> hSourceTerm ;
    public:

      // Define input and output.
      SOUInviscidFluxToTotalEnthalpySourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("hGradient",hGradient) ;
        name_store("limiters(h)",hLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("Z",massFlux) ;
        name_store("h_f",h_f) ;
        name_store("hSourceTerm(X,Y,Z)",hSourceTerm) ;
        input("ci->(hGradient,limiters(h),cellcenter)") ;
        input("facecenter,Z,h_f") ;
        output("ci->hSourceTerm(X,Y,Z)") ;
        constraint("boundaryFaces,souOrRoeInviscidFlux") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        if(massFlux[face]>0.0) hSourceTerm[ci[face]]-=massFlux[face]*
          hLimiter[ci[face]]*dot(hGradient[ci[face]],faceCenter[face]-
          cellCenter[ci[face]]) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SOUInviscidFluxToTotalEnthalpySourceTermBoundary>
    registerSOUInviscidFluxToTotalEnthalpySourceTermBoundary ;

  // Rule to add the heat flux contribution to the source term for interior
  // faces.
  class HeatFluxToTotalEnthalpySourceTermInterior : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> vol ;
      const_store<vect3d> temperatureGradient ;
      const_store<real> tDiffusionCoeff ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<real> hSourceTerm ;
    public:

      // Define input and output.
      HeatFluxToTotalEnthalpySourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("vol",vol) ;
        name_store("grads(temperature)",temperatureGradient) ;
        name_store("tDiffusionCoeff(temperature,p,y,eddyViscosity)",
          tDiffusionCoeff) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("hSourceTerm(X,Y,Z)",hSourceTerm) ;
        input("(cl,cr)->(vol,grads(temperature))") ;
        input("(cl,cr)->tDiffusionCoeff(temperature,p,y,eddyViscosity)") ;
        input("area,faceRadius") ;
        output("(cl,cr)->hSourceTerm(X,Y,Z)") ;
        constraint("internalFaces,viscousFlow") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        real sourceTerm=0.5*area[face].sada*(tDiffusionCoeff[cl[face]]+
          tDiffusionCoeff[cr[face]])*dot(temperatureGradient[cl[face]]*
          vol[cr[face]]+temperatureGradient[cr[face]]*vol[cl[face]],
          area[face].n)*faceRadius[face]/(vol[cl[face]]+vol[cr[face]]) ;
        hSourceTerm[cl[face]]+=sourceTerm ; hSourceTerm[cr[face]]-=sourceTerm ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<HeatFluxToTotalEnthalpySourceTermInterior>
    registerHeatFluxToTotalEnthalpySourceTermInterior ;

  // Rule to add the heat flux contribution to the source term for boundary
  // faces.
  class HeatFluxToTotalEnthalpySourceTermBoundary : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> qWall ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<real> hSourceTerm ;
    public:

      // Define input and output.
      HeatFluxToTotalEnthalpySourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("qWall",qWall) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("hSourceTerm(X,Y,Z)",hSourceTerm) ;
        input("qWall,area,faceRadius") ;
        output("ci->hSourceTerm(X,Y,Z)") ;
        constraint("viscousFlow,boundaryFaces") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        hSourceTerm[ci[face]]-=qWall[face]*area[face].sada*faceRadius[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<HeatFluxToTotalEnthalpySourceTermBoundary>
    registerHeatFluxToTotalEnthalpySourceTermBoundary ;

  // Rule to add the viscous stress contribution to the source term for
  // interior faces.
  class ViscousStressToTotalEnthalpySourceTermInterior : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> vol ;
      const_store<vect3d> v ;
      const_store<tens3d> vGradient ;
      const_store<real> viscosity ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<real> hSourceTerm ;
    public:

      // Define input and output.
      ViscousStressToTotalEnthalpySourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("vol",vol) ;
        name_store("X",v) ;
        name_store("gradv3d(X)",vGradient) ;
        name_store("viscosity",viscosity) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("hSourceTerm(X,Y,Z)",hSourceTerm) ;
        input("(cl,cr)->(vol,X,gradv3d(X)),viscosity") ;
        input("area,faceRadius") ;
        output("(cl,cr)->hSourceTerm(X,Y,Z)") ;
        constraint("internalFaces,viscousFlow") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real vR=vol[cr[face]]/(vol[cl[face]]+vol[cr[face]]) ;
        real vL=vol[cl[face]]/(vol[cl[face]]+vol[cr[face]]) ;
        tens3d tempVGradient=product(vR,vGradient[cl[face]])+product(vL,
          vGradient[cr[face]]) ;
        tens3d tempStress=tempVGradient+Transpose(tempVGradient) ;
        real temp=2.0*Trace(tempVGradient)/3.0 ;
        tempStress.x.x-=temp ; tempStress.y.y-=temp ; tempStress.z.z-=temp ;
        real sourceTerm=0.5*area[face].sada*viscosity[face]*
          dot(dotTemp(tempStress,v[cl[face]]+v[cr[face]]),area[face].n)*
          faceRadius[face] ;
        hSourceTerm[cl[face]]+=sourceTerm ; hSourceTerm[cr[face]]-=sourceTerm ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ViscousStressToTotalEnthalpySourceTermInterior>
    registerViscousStressToTotalEnthalpySourceTermInterior ;

  // Rule to add the viscous stress contribution to the source term for
  // boundary faces. Note that we are always using the non-wall function form
  // without the wall stress coming via a wall function. Since the face
  // velocity is always zero for no-slip faces, there is no contribution
  // anyway.
  class ViscousStressToTotalEnthalpySourceTermBoundary : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<tens3d> vGradient ;
      const_store<vect3d> v_f ;
      const_store<real> viscosity ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<real> hSourceTerm ;
    public:

      // Define input and output.
      ViscousStressToTotalEnthalpySourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("gradv3d(X)",vGradient) ;
        name_store("X_f",v_f) ;
        name_store("viscosity",viscosity) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("hSourceTerm(X,Y,Z)",hSourceTerm) ;
        input("ci->gradv3d(X),X_f,viscosity") ;
        input("area,faceRadius") ;
        output("ci->hSourceTerm(X,Y,Z)") ;
        constraint("boundaryEnergyDiffusion,viscousFlow") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        tens3d tempStress=vGradient[ci[face]]+Transpose(vGradient[ci[face]]) ;
        real temp=2.0*Trace(vGradient[ci[face]])/3.0 ; ;
        tempStress.x.x-=temp ; tempStress.y.y-=temp ; tempStress.z.z-=temp ;
        real sourceTerm=area[face].sada*viscosity[face]*
          dot(dotTemp(tempStress,v_f[face]),area[face].n)*faceRadius[face] ;
        hSourceTerm[ci[face]]+=sourceTerm ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ViscousStressToTotalEnthalpySourceTermBoundary>
    registerViscousStressToTotalEnthalpySourceTermBoundary ;

  // Rule to add the Roe dissipation contribution to the source term for
  // interior faces.
  class RoeDissipationToTotalEnthalpySourceTermInterior : public apply_rule
  <store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
      const_store<real> delP,deltaP,vNormalTilde ;
      const_store<Area> area ;
      store<real> hSourceTerm ;
    public:
                                                                                
      // Define input and output.
      RoeDissipationToTotalEnthalpySourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("delP",delP) ;
        name_store("deltaP",deltaP) ;
        name_store("vNormalTilde",vNormalTilde) ;
        name_store("area",area) ;
        name_store("hSourceTerm(X,Y,Z)",hSourceTerm) ;
        input("delP,deltaP,vNormalTilde,area") ;
        output("(cl,cr)->hSourceTerm(X,Y,Z)") ;
        constraint("internalFaces,roeInviscidFlux") ;
      }
                                                                                
      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        real source=0.5*(vNormalTilde[face]*delP[face]-abs(vNormalTilde[face])*
          deltaP[face])*area[face].sada ;
        hSourceTerm[cl[face]]+=source ; hSourceTerm[cr[face]]-=source ;
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<RoeDissipationToTotalEnthalpySourceTermInterior>
    registerRoeDissipationToTotalEnthalpySourceTermInterior ;

  // Rule to add the buoyancy contribution to the source term.
  class BuoyancyToTotalEnthalpySourceTerm : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<vect3d> gravity ;
      const_param<real> rhoReference ;
      const_store<real> rho ;
      const_store<vect3d> v ;
      const_store<real> vol ;
      store<real> hSourceTerm ;
    public:
                                                                                
      // Define input and output.
      BuoyancyToTotalEnthalpySourceTerm() {
        name_store("gravityAcceleration",gravity) ;
        name_store("gravityRhoRef",rhoReference) ;
        name_store("rho",rho) ;
        name_store("X",v) ;
        name_store("vol",vol) ;
        name_store("hSourceTerm(X,Y,Z)",hSourceTerm) ;
        input("gravityAcceleration,gravityRhoRef,rho,X,vol") ;
        output("hSourceTerm(X,Y,Z)") ;
        constraint("geom_cells") ;
      }
                                                                                
      // Add buoyance force for single cell.
      void calculate(Entity cell) {
        hSourceTerm[cell]+=(rho[cell]-*rhoReference)*dot(*gravity,v[cell])*
          vol[cell] ;
      }
                                                                                
      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<BuoyancyToTotalEnthalpySourceTerm>
    registerBuoyancyToTotalEnthalpySourceTerm ;

  // Rule to add the species enthalpy diffusion to the source term.
  class SpeciesDiffusionToTotalEnthalpySourceTerm : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_store<real> hDiffusion ;
      store<real> hSourceTerm ;
    public:

      // Define input and output.
      SpeciesDiffusionToTotalEnthalpySourceTerm() {
        name_store("hDiffusion",hDiffusion) ;
        name_store("hSourceTerm(X,Y,Z)",hSourceTerm) ;
        input("hDiffusion") ;
        output("hSourceTerm(X,Y,Z)") ;
        constraint("speciesTransport,viscousFlow,geom_cells") ;
      }

      // Add species enthalpy diffusion for single cell.
      void calculate(Entity cell) { hSourceTerm[cell]+=hDiffusion[cell] ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SpeciesDiffusionToTotalEnthalpySourceTerm>
    registerSpeciesDiffusionToTotalEnthalpySourceTerm ;

  // Rule to compute the temporal source term for BDF. Note that this
  // term only includes the 'old' components from the {n} level. This
  // rule will compute and be used at {n}. In addition, the value will
  // promote to {n,it} where it is also needed.
  class TotalEnthalpyTemporalSourceTermBDF : public pointwise_rule {
    private:
      const_param<real> timeStep ;
      const_store<real> timeStepFactor ;
      const_store<real> rho,p,h,vol,cellRadius ;
      store<real> hTemporalSourceTerm ;
    public:

      // Define input and output.
      TotalEnthalpyTemporalSourceTermBDF() {
        name_store("timeStep{n}",timeStep) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n}",rho) ;
        name_store("p{n}",p) ;
        name_store("h{n}",h) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n}",cellRadius) ;
        name_store("hTemporalSourceTerm{n}",hTemporalSourceTerm) ;
        input("rho{n},p{n},h{n},vol{n},cellRadius{n}") ;
        input("timeStep{n},timeStepFactor{n}") ;
        output("hTemporalSourceTerm{n}") ;
        constraint("BDFIntegrator,geom_cells") ;
      }

      // Compute for a single cell.
      void calculate(Entity cell) {
        hTemporalSourceTerm[cell]=(rho[cell]*h[cell]-p[cell])*vol[cell]*
          cellRadius[cell]/((*timeStep)*timeStepFactor[cell]) ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TotalEnthalpyTemporalSourceTermBDF>
    registerTotalEnthalpyTemporalSourceTermBDF ;

  // Rule to compute the temporal source term for BDF2. Note that this
  // term only includes the 'old' components from the {n-1} and {n} levels.
  // This rule will compute and be used at {n}. In addition, the value will
  // promote to {n,it} where it is also needed.
  class TotalEnthalpyTemporalSourceTermBDF2 : public pointwise_rule {
    private:
      const_param<int> n ;
      const_param<real> timeStep ;
      const_store<real> timeStepFactor ;
      const_store<real> rhoOld,rho ;
      const_store<real> pOld,p ;
      const_store<real> hOld,h ;
      const_store<real> volOld,vol ;
      const_store<real> cellRadius ;
      store<real> hTemporalSourceTerm ;
    public:

      // Define input and output.
      TotalEnthalpyTemporalSourceTermBDF2() {
        name_store("$n{n}",n) ;
        name_store("timeStep{n}",timeStep) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n-1}",rhoOld) ;
        name_store("rho{n}",rho) ;
        name_store("p{n-1}",pOld) ;
        name_store("p{n}",p) ;
        name_store("h{n-1}",hOld) ;
        name_store("h{n}",h) ;
//      name_store("vol{n-1}",volOld) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n}",cellRadius) ;
        name_store("hTemporalSourceTerm{n}",hTemporalSourceTerm) ;
//      input("$n{n},rho{n-1},rho{n},p{n-1},p{n},h{n-1},h{n},vol{n-1},vol{n}") ;
        input("$n{n},rho{n-1},rho{n},p{n-1},p{n},h{n-1},h{n},vol{n}") ;
        input("cellRadius{n},timeStep{n},timeStepFactor{n}") ;
        output("hTemporalSourceTerm{n}") ;
        constraint("BDF2Integrator,geom_cells") ;
      }

      // Compute for a single cell.  Use BDF on first timestep.
      void calculate(Entity cell) {
        if((*n)!=0){
          hTemporalSourceTerm[cell]=(2.0*(rho[cell]*h[cell]+p[cell])*vol[cell]-
            0.5*(rhoOld[cell]*hOld[cell]-pOld[cell])*vol[cell])*
            cellRadius[cell]/((*timeStep)*timeStepFactor[cell]) ;
        }else{
          hTemporalSourceTerm[cell]=(rho[cell]*h[cell]-p[cell])*vol[cell]*
            cellRadius[cell]/((*timeStep)*timeStepFactor[cell]) ;
        }
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TotalEnthalpyTemporalSourceTermBDF2>
    registerTotalEnthalpyTemporalSourceTermBDF2 ;

  // Rule to add temporal component of the total enthalpy equation to the
  // source term. This adds in the terms from the old time level(s) plus the
  // additional term for the pressure time derivative based on the new
  // pressure.
  class TemporalToTotalEnthalpySourceTerm : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> timeStep ;
      const_param<real> timeIntegratorFactor0 ;
      const_store<real> timeStepFactor ;
      const_store<real> p ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      const_store<real> hTemporalSourceTerm ;
      store<real> hSourceTerm ;
    public:

      // Define input and output.
      TemporalToTotalEnthalpySourceTerm() {
        name_store("timeStep",timeStep) ;
        name_store("timeIntegratorFactor0",timeIntegratorFactor0) ;
        name_store("timeStepFactor",timeStepFactor) ;
        name_store("Y",p) ;
        name_store("vol",vol) ;
        name_store("cellRadius",cellRadius) ;
        name_store("hTemporalSourceTerm",hTemporalSourceTerm) ;
        name_store("hSourceTerm(X,Y,Z)",hSourceTerm) ;
        input("timeStep,timeStepFactor,timeIntegratorFactor0") ;
        input("Y,vol,cellRadius,hTemporalSourceTerm") ;
        output("hSourceTerm(X,Y,Z)") ;
        constraint("geom_cells") ;
      }

      // Add temporal component to source term for a single cell.
      void calculate(Entity cell) {
        hSourceTerm[cell]+=hTemporalSourceTerm[cell]+0.5*p[cell]*vol[cell]*
         ((*timeIntegratorFactor0)+1.0)*cellRadius[cell]/((*timeStep)*
         timeStepFactor[cell]) ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToTotalEnthalpySourceTerm>
    registerTemporalToTotalEnthalpySourceTerm ;

  // Rule to compute the diagonal term for the linear system.
  class ComputeTotalEnthalpyMatrixDiagonal : public pointwise_rule {
    private:
      const_store<real> hMainCoefficient ;
      store<real> D ;
    public:

      // Define input and output.
      ComputeTotalEnthalpyMatrixDiagonal() {
        name_store("hMainCoefficient(massFluxCorrected_p,rhoStar)",hMainCoefficient) ;
        name_store("hStar_D",D) ;
        input("hMainCoefficient(massFluxCorrected_p,rhoStar)") ;
        output("hStar_D") ;
        constraint("geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) { D[cell]=hMainCoefficient[cell] ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeTotalEnthalpyMatrixDiagonal>
    registerComputeTotalEnthalpyMatrixDiagonal ;

  // Rule to copy hStar_D for periodic faces.
  class ComputeTotalEnthalpyMatrixDiagonalPeriodic : public pointwise_rule {
    private:
      const_Map cl,cr,pmap ;
      store<real> D ;
    public:

      // Define input and output.
      ComputeTotalEnthalpyMatrixDiagonalPeriodic() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pmap",pmap) ;
        name_store("hStar_D",D) ;
        input("pmap->cl->hStar_D") ;
        output("cr->hStar_D") ;
      }

      // For a face.
      void calculate(Entity face) { D[cr[face]]=D[cl[pmap[face]]] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; } } ;

  register_rule<ComputeTotalEnthalpyMatrixDiagonalPeriodic>
    registerComputeTotalEnthalpyMatrixDiagonalPeriodic ;

  // Rule to compute the diagonal term for the linear system.
  class ComputeTotalEnthalpyMatrixDiagonalCorrector : public pointwise_rule {
    private:
      const_store<real> hMainCoefficient ;
      store<real> D ;
    public:

      // Define input and output.
      ComputeTotalEnthalpyMatrixDiagonalCorrector() {
        name_store("hMainCoefficient(massFluxCorrected_c,rhoCorrected)",hMainCoefficient) ;
        name_store("hCorrected_D",D) ;
        input("hMainCoefficient(massFluxCorrected_c,rhoCorrected)") ;
        output("hCorrected_D") ;
        constraint("geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) { D[cell]=hMainCoefficient[cell] ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeTotalEnthalpyMatrixDiagonalCorrector>
    registerComputeTotalEnthalpyMatrixDiagonalCorrector ;

  // Rule to initialize the lower terms for the linear system.  This is now
  // a parametric variable with X=massFlux.
  class InitializeTotalEnthalpyMatrixLower : public unit_rule {
    private:
      store<real> L ;
    public:

      // Define input and output.
      InitializeTotalEnthalpyMatrixLower() {
        name_store("hStar_L(X)",L) ;
        output("hStar_L(X)") ;
        constraint("internalFaces") ;
      }

      // Initialize for a single face.
      void calculate(Entity face) { L[face]=0.0 ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeTotalEnthalpyMatrixLower>
    registerInitializeTotalEnthalpyMatrixLower ;

  // Rule to add the first-order inviscid flux contribution to the lower terms
  // for the linear system. Checked.
  class FOUInviscidFluxToTotalEnthalpyMatrixLower : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_store<real> massFlux ;
      store<real> L ;
    public:

      // Define input and output.
      FOUInviscidFluxToTotalEnthalpyMatrixLower() {
        name_store("X",massFlux) ;
        name_store("hStar_L(X)",L) ;
        input("X") ;
        output("hStar_L(X)") ;
        constraint("internalFaces") ;
      }

      // Increment the lower term for a single face. Note that the increment
      // is the negative of the one in streamUns since this coefficient is
      // for a term on the lhs of the equation. In streamUns, the coefficient
      // is associated with a term on the rhs.
      void calculate(Entity face) {
        if(massFlux[face]>0.0) L[face]-=massFlux[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToTotalEnthalpyMatrixLower>
    registerFOUInviscidFluxToTotalEnthalpyMatrixLower ;

  // Rule to initialize the upper terms for the linear system.  This is now
  // a parametric variable with X=massFlux.
  class InitializeTotalEnthalpyMatrixUpper : public unit_rule {
    private:
      store<real> U ;
    public:

      // Define input and output.
      InitializeTotalEnthalpyMatrixUpper() {
        name_store("hStar_U(X)",U) ;
        output("hStar_U(X)") ;
        constraint("internalFaces") ;
      }

      // Initialize for a single face.
      void calculate(Entity face) { U[face]=0.0 ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeTotalEnthalpyMatrixUpper>
    registerInitializeTotalEnthalpyMatrixUpper ;

  // Rule to add the first-order inviscid flux contribution to the upper terms
  // for the linear system. Checked.
  class FOUInviscidFluxToTotalEnthalpyMatrixUpper : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_store<real> massFlux ;
      store<real> U ;
    public:

      // Define input and output.
      FOUInviscidFluxToTotalEnthalpyMatrixUpper() {
        name_store("X",massFlux) ;
        name_store("hStar_U(X)",U) ;
        input("X") ;
        output("hStar_U(X)") ;
        constraint("internalFaces") ;
      }

      // Increment the upper term for a single face.
      void calculate(Entity face) {
        if(massFlux[face]<0.0) U[face]+=massFlux[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToTotalEnthalpyMatrixUpper>
    registerFOUInviscidFluxToTotalEnthalpyMatrixUpper ;

  // Rule to compute the right-hand side for the linear system.
  class ComputeTotalEnthalpyRHS : public pointwise_rule {
    private:
      const_store<real> hSourceTerm ;
      store<real> B ;
    public:

      // Define input and output.
      ComputeTotalEnthalpyRHS() {
        name_store("hSourceTerm(vCorrected_p,pCorrected_p,massFluxCorrected_p)",
          hSourceTerm) ;
        name_store("hStar_B",B) ;
        input("hSourceTerm(vCorrected_p,pCorrected_p,massFluxCorrected_p)") ;
        output("hStar_B") ;
        constraint("geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) { B[cell]=hSourceTerm[cell] ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeTotalEnthalpyRHS> registerComputeTotalEnthalpyRHS ;

  // Rule to copy hStar_B for periodic faces.
  class ComputeTotalEnthalpyRHSPeriodic : public pointwise_rule {
    private:
      const_Map cl,cr,pmap ;
      store<real> B ;
    public:

      // Define input and output.
      ComputeTotalEnthalpyRHSPeriodic() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pmap",pmap) ;
        name_store("hStar_B",B) ;
        input("pmap->cl->hStar_B") ;
        output("cr->hStar_B") ;
      }

      // For a face.
      void calculate(Entity face) { B[cr[face]]=B[cl[pmap[face]]] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeTotalEnthalpyRHSPeriodic>
    registerComputeTotalEnthalpyRHSPeriodic ;

  // Rule to compute the right-hand side for the linear system.
  class ComputeTotalEnthalpyRHSCorrector : public pointwise_rule {
    private:
      const_store<real> hSourceTerm ;
      store<real> B ;
    public:

      // Define input and output.
      ComputeTotalEnthalpyRHSCorrector() {
        name_store("hSourceTerm(vCorrected_c,pCorrected_c,massFluxCorrected_c)",
          hSourceTerm) ;
        name_store("hCorrected_B",B) ;
        input("hSourceTerm(vCorrected_c,pCorrected_c,massFluxCorrected_c)") ;
        output("hCorrected_B") ;
        constraint("geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) { B[cell]=hSourceTerm[cell] ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeTotalEnthalpyRHSCorrector>
    registerComputeTotalEnthalpyRHSCorrector ;

  // Total enthalpy corrector.
  class ComputeTotalEnthalpyCorrected : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_multiMap upper,lower ;
      const_store<real> h ;
      const_store<real> hMainCoefficient ;
      const_store<real> hStar_L,hStar_U ;
      const_store<real> hSourceTerm ;
      store<real> hCorrected ;
    public:

      // Define input and output.
      ComputeTotalEnthalpyCorrected() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("h",h) ;
        name_store("hMainCoefficient(massFluxCorrected_c,rhoCorrected)",hMainCoefficient) ;
        name_store("hStar_L(massFluxCorrected_c",hStar_L) ;
        name_store("hStar_U(massFluxCorrected_c)",hStar_U) ;
        name_store("hSourceTerm(vCorrected_c,pCorrected_c,massFluxCorrected_c)",hSourceTerm) ;
        name_store("hCorrected",hCorrected) ;
        input("upper->cr->h,lower->cl->h") ;
        input("hMainCoefficient(massFluxCorrected_c,rhoCorrected)") ;
        input("lower->hStar_L(massFluxCorrected_c)") ;
        input("upper->hStar_U(massFluxCorrected_c)") ;
        input("hSourceTerm(vCorrected_c,pCorrected_c,massFluxCorrected_c)") ;
        output("hCorrected") ;
        constraint("geom_cells") ;
      }

      // Add explicit and implicit contributions for the cell, then divide by the
      // main coefficient.
      void calculate(Entity cell) {
        hCorrected[cell]=hSourceTerm[cell] ;
        for(const int *ui=upper.begin(cell);ui!=upper.end(cell);++ui)
          hCorrected[cell]-=hStar_U[*ui]*h[cr[*ui]] ;
        for(const int *li=lower.begin(cell);li!=lower.end(cell);++li)
          hCorrected[cell]-=hStar_L[*li]*h[cl[*li]] ;
        hCorrected[cell]/=hMainCoefficient[cell] ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeTotalEnthalpyCorrected>
    registerComputeTotalEnthalpyCorrected ;

//-----------------------------------------------------------------------------
// Rules for computing the residual of the energy equation.

  // Rule to initialize the energy residual. Checked.
  class InitializeEnergyResidual : public unit_rule {
    private:
      store<real> hResidual ;
    public:

      // Define input and output.
      InitializeEnergyResidual() {
        name_store("hResidual",hResidual) ;
        output("hResidual") ;
        constraint("compressibleFlow,vol") ;
      }

      // Initialize the residual for a single cell.
      void calculate(Entity cell) { hResidual[cell]=0.0 ; }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

//register_rule<InitializeEnergyResidual> registerInitializeEnergyResidual ;

  // Rule to compute the energy residual for each cell. Checked.
  class ComputeEnergyResidualOne : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> rhoScale,vScale,lScale,hScale ;
      const_store<real> D ;
      const_store<real> h ;
      const_store<real> B ;
      store<real> hResidual ;
    private:
      real hFactor ;
    public:

      // Define input and output.
      ComputeEnergyResidualOne() {
        name_store("rhoScale",rhoScale) ;
        name_store("vScale",vScale) ;
        name_store("lScale",lScale) ;
        name_store("hScale",hScale) ;
        name_store("hStar_D",D) ;
        name_store("h",h) ;
        name_store("hStar_B",B) ;
        name_store("hResidual",hResidual) ;
        input("rhoScale,vScale,lScale,hScale,hStar_D,h,hStar_B") ;
        output("hResidual") ;
        constraint("vol") ;
      }

      // Initialize the residual for a single cell.
      void calculate(Entity cell) {
        hResidual[cell]+=(B[cell]-D[cell]*h[cell])/hFactor ;
      }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) {
        hFactor=(*rhoScale)*(*vScale)*(*hScale)*(*lScale)*(*lScale) ;
        do_loop(seq,this) ;
      }
  } ;

//register_rule<ComputeEnergyResidualOne> registerComputeEnergyResidualOne ;

  // Rule to compute the energy residual for each cell. Checked.
  class ComputeEnergyResidualTwo : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_param<real> rhoScale,vScale,lScale,hScale ;
      const_store<real> h ;
      const_store<real> L,U ;
      store<real> hResidual ;
    private:
      real hFactor ;
    public:

      // Define input and output.
      ComputeEnergyResidualTwo() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("rhoScale",rhoScale) ;
        name_store("vScale",vScale) ;
        name_store("lScale",lScale) ;
        name_store("hScale",hScale) ;
        name_store("h",h) ;
        name_store("hStar_L",L) ;
        name_store("hStar_U",U) ;
        name_store("hResidual",hResidual) ;
        input("rhoScale,vScale,lScale,hScale,(cl,cr)->h,hStar_L,hStar_U") ;
        output("(cl,cr)->hResidual") ;
        constraint("internalFaces") ;
      }

      // Add the neighbor contribution to the residual for each of the two
      // cells on either side of the face.
      void calculate(Entity face) {
        hResidual[cl[face]]-=U[face]*h[cr[face]]/hFactor ;
        hResidual[cr[face]]-=L[face]*h[cl[face]]/hFactor ;
      }

      // Add the neighbor contributions for a sequence of faces.
      virtual void compute(const sequence &seq) {
        hFactor=(*rhoScale)*(*vScale)*(*hScale)*(*lScale)*(*lScale) ;
        do_loop(seq,this) ;
      }
  } ;

//register_rule<ComputeEnergyResidualTwo> registerComputeEnergyResidualTwo ;

  // Rule to set the total energy residual to zero for incompressible flow.
  class TotalEnergyResidualIncompressible : public singleton_rule {
    private:
      param<ScalarResidual> hResidualData ;
    public:

      // Define input and output.
      TotalEnergyResidualIncompressible() {
        name_store("hResidualData",hResidualData) ;
        output("hResidualData") ;
        constraint("incompressibleFlow,geom_cells") ;
      }

      // Set the residual.
      virtual void compute(const sequence &seq) {
        hResidualData=ScalarResidual() ;
      }
  } ;

  register_rule<TotalEnergyResidualIncompressible>
    registerTotalEnergyResidualIncompressible ;

  // Rule to initialize the total energy residual. Checked.
  class InitializeTotalEnergyResidual : public unit_rule {
    private:
      param<ScalarResidual> hResidualData ;
    public:

      // Define input and output.
      InitializeTotalEnergyResidual() {
        name_store("hResidualData",hResidualData) ;
        output("hResidualData") ;
        constraint("compressibleFlow,geom_cells") ;
      }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) {
        hResidualData=ScalarResidual() ;
      }
  } ;

//register_rule<InitializeTotalEnergyResidual>
//  registerInitializeTotalEnergyResidual ;

  // Rule to compute the total energy residual. Checked.
  class ComputeTotalEnergyResidual : public apply_rule<param<ScalarResidual>,
  ScalarResidualJoin> {
    private:
      const_store<real> hResidual ;
      const_store<vect3d> cellCenter ;
      param<ScalarResidual> hResidualData ;
    public:

      // Define input and output.
      ComputeTotalEnergyResidual() {
        name_store("hResidual",hResidual) ;
        name_store("cellcenter",cellCenter) ;
        name_store("hResidualData",hResidualData) ;
        input("hResidual,cellcenter") ;
        output("hResidualData") ;
        constraint("compressibleFlow,geom_cells") ;
      }

      // Add the cell contribution to the residual for a single cell.
      void calculate(Entity cell) {
        ScalarResidual temp ;
        temp.maxResidual=hResidual[cell] ;
        temp.totalResidual=abs(hResidual[cell]) ;
        temp.maxResidualLocation=cellCenter[cell] ;
        join(*hResidualData,temp) ;
      }

      // Add the cell contribution to the residual for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

//register_rule<ComputeTotalEnergyResidual>
//  registerComputeTotalEnergyResidual ;

//-----------------------------------------------------------------------------
// Rules for marching the temperature.

  // Time build rule for temperature.
  class TimeBuildTemperature : public pointwise_rule {
    private:
      const_store<real> T_ic ;
      store<real> temperatureTimeStepZero ;
    public:

      // Define input and output.
      TimeBuildTemperature() {
        name_store("T_ic",T_ic) ;
        name_store("temperature{n=0}",temperatureTimeStepZero) ;
        input("T_ic") ;
        output("temperature{n=0}") ;
        constraint("compressibleFlow,geom_cells") ;
      }

      // Assign temperature at time zero for a single cell.
      void calculate(Entity cell) {
        temperatureTimeStepZero[cell]=T_ic[cell] ;
      }

      // Assign temperature at time zero for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildTemperature> registerTimeBuildTemperature ;

  // Iteration build rule for temperature.
  class IterationBuildTemperature : public pointwise_rule {
    private:
      const_store<real> temperatureStar ;
      store<real> temperature ;
    public:

      // Define input and output.
      IterationBuildTemperature() {
        name_store("temperatureStar{n}",temperatureStar) ;
        name_store("temperature{n,it=0}",temperature) ;
        input("temperatureStar{n}") ;
        output("temperature{n,it=0}") ;
        constraint("compressibleFlow{n},geom_cells{n}") ;
      }

      // Assign temperature at iteration zero for a single cell.
      void calculate(Entity cell) {
        temperature[cell]=temperatureStar[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationBuildTemperature> registerIterationBuildTemperature ;

  // Iteration advance rule for temperature.
  class IterationAdvanceTemperature : public pointwise_rule {
    private:
      const_store<real> temperatureCorrected ;
      store<real> temperature ;
    public:

      // Define input and output.
      IterationAdvanceTemperature() {
        name_store("temperatureCorrected{n,it}",temperatureCorrected) ;
        name_store("temperature{n,it+1}",temperature) ;
        input("temperatureCorrected{n,it}") ;
        output("temperature{n,it+1}") ;
        constraint("geom_cells{n,it},compressibleFlow{n,it}") ;
      }

      // Assign temperature at end of iteration for a single cell.
      void calculate(Entity cell) {
        temperature[cell]=temperatureCorrected[cell] ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationAdvanceTemperature>
    registerIterationAdvanceTemperature ;

  // Iteration collapse rule for temperature.
  class IterationCollapseTemperature : public pointwise_rule {
    private:
      store<real> temperature ;
    public:

      // Define input and output.
      IterationCollapseTemperature() {
        name_store("temperature{n,it}",temperature) ;
        input("temperature{n,it}") ;
        output("temperature{n+1}=temperature{n,it}") ;
        conditional("iterationFinished{n,it-1}") ;
        constraint("compressibleFlow{n,it},geom_cells{n,it}") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<IterationCollapseTemperature>
    registerIterationCollapseTemperature ;

//-----------------------------------------------------------------------------
// Rules for marching the total enthalpy.

  // Time build rule for total enthalpy when using BDF2 time integrator.
  class TimeBuildTotalEnthalpyBDF2 : public pointwise_rule {
    private:
      const_store<real> h_ic ;
      store<real> h ;
    public:

      // Define input and output.
      TimeBuildTotalEnthalpyBDF2() {
        name_store("h_ic",h_ic) ;
        name_store("h{n=-1}",h) ;
        input("h_ic") ;
        output("h{n=-1}") ;
        constraint("compressibleFlow,geom_cells") ;
      }

      // Assign total enthalpy at time zero for a single cell.
      void calculate(Entity cell) { h[cell]=h_ic[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildTotalEnthalpyBDF2>
    registerTimeBuildTotalEnthalpyBDF2 ;

  // Time build rule for total enthalpy.
  class TimeBuildTotalEnthalpy : public pointwise_rule {
    private:
      const_store<real> h_ic ;
      store<real> hTimeStepZero ;
    public:

      // Define input and output.
      TimeBuildTotalEnthalpy() {
        name_store("h_ic",h_ic) ;
        name_store("h{n=0}",hTimeStepZero) ;
        input("h_ic") ;
        output("h{n=0}") ;
        constraint("compressibleFlow,geom_cells") ;
      }

      // Assign total enthalpy at time zero for a single cell.
      void calculate(Entity cell) {
        hTimeStepZero[cell]=h_ic[cell] ;
      }

      // Assign total enthalpy at time zero for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildTotalEnthalpy> registerTimeBuildTotalEnthalpy ;

  // Iteration build rule for total enthalpy.
  class IterationBuildTotalEnthalpy : public pointwise_rule {
    private:
      const_store<real> hStar ;
      store<real> h ;
    public:

      // Define input and output.
      IterationBuildTotalEnthalpy() {
        name_store("hStar{n}",hStar) ;
        name_store("h{n,it=0}",h) ;
        input("hStar{n}") ;
        output("h{n,it=0}") ;
        constraint("geom_cells{n},compressibleFlow{n}") ;
      }

      // Set for a single cell.
      void calculate(Entity cell) { h[cell]=hStar[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationBuildTotalEnthalpy>
    registerIterationBuildTotalEnthalpy ;

  // Iteration advance rule for total enthalpy.
  class IterationAdvanceTotalEnthalpy : public pointwise_rule {
    private:
      const_store<real> hCorrected ;
      store<real> h ;
    public:

      // Define input and output.
      IterationAdvanceTotalEnthalpy() {
        name_store("hCorrected{n,it}",hCorrected) ;
        name_store("h{n,it+1}",h) ;
        input("hCorrected{n,it}") ;
        output("h{n,it+1}") ;
        constraint("geom_cells{n,it},compressibleFlow{n,it}") ;
      }

      // Assign total enthalpy at end of iteration for a single cell.
      void calculate(Entity cell) { h[cell]=hCorrected[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationAdvanceTotalEnthalpy>
    registerIterationAdvanceTotalEnthalpy ;

  // Iteration collapse rule for total enthalpy.
  class IterationCollapseTotalEnthalpy : public pointwise_rule {
    private:
      store<real> h ;
    public:

      // Define input and output.
      IterationCollapseTotalEnthalpy() {
        name_store("h{n,it}",h) ;
        input("h{n,it}") ;
        output("h{n+1}=h{n,it}") ;
        conditional("iterationFinished{n,it-1}") ;
        constraint("compressibleFlow{n,it},geom_cells{n,it}") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<IterationCollapseTotalEnthalpy>
    registerIterationCollapseTotalEnthalpy ;
}

