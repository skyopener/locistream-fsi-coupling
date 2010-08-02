//-----------------------------------------------------------------------------
// Description: This file contains rules for the energy equation, including
//    rules for temperature and total enthalpy boundary conditions.
//
// Authors: Siddharth Thakur and Jeff Wright
//          Streamline Numerics Inc.
//          Gainesville, Florida 32609
//-----------------------------------------------------------------------------

// Standard library includes.
#include <map>
#include <vector>
using std::map ;
using std::vector ;

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;

// StreamUns includes.
#include "bcInput.h"
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
      param<real> hRelaxationFactor ;
    public:
                                                                                
      // Define input and output.
      EnergyEquationSolverParameters() {
        name_store("energyEquationOptions",energyEquationOptions) ;
        name_store("hStar_maxLinearSolverIterations",hMaxIterations) ;
        name_store("hStar_RelaxationFactor",hRelaxationFactor) ;
        input("energyEquationOptions") ;
        output("hStar_maxLinearSolverIterations,hStar_RelaxationFactor") ;
      }
                                                                                
      // Set up the parameters.
      virtual void compute(const sequence& seq) {
                                                                                
        // Relaxation factor.
        if((*energyEquationOptions).optionExists("relaxationFactor")){
          Loci::option_value_type optionValueType=energyEquationOptions->
            getOptionValueType("relaxationFactor") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                energyEquationOptions->getOption("relaxationFactor",temp) ;
                if(temp<=0.0 || temp>1.0){
                  cerr << "Bad relaxationFactor for energyEquation." << endl ;
                  Loci::Abort() ;
                }
                *hRelaxationFactor=temp ;
              }
              break ;
            default:
              cerr << "Bad type for relaxationFactor in energyEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *hRelaxationFactor=0.5 ;
        }
                                                                                
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
// Reservoir temperature boundary condition rules.

  // Rule for specifying constant Treservoir. Assigns value to all boundary
  // faces that have the property constantTReservoir_BC.
  class BoundaryTreservoirConstant : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<real> constantTreservoir_BC ;
      store<real> Treservoir_f ;
    public:

      // Define input and output.
      BoundaryTreservoirConstant() {
        name_store("ref",ref) ;
        name_store("constantTreservoir_BC",constantTreservoir_BC) ;
        name_store("Treservoir_f",Treservoir_f) ;
        input("ref->constantTreservoir_BC") ;
        output("Treservoir_f") ;
      }

      // Assign reservoir temperature for a single face.
      void calculate(Entity face) {
        Treservoir_f[face]=constantTreservoir_BC[ref[face]] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryTreservoirConstant> registerBoundaryTreservoirConstant ;

  // Rule for specifying Treservoir as a function of x, y and z.
  class BoundaryTreservoirFunctionSteady : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<Loci::options_list> BC_options ;
      const_store<vect3d> faceCenter ;
      store<real> Treservoir_f ;
    public:

      // Define input and output.
      BoundaryTreservoirFunctionSteady() {
        name_store("ref",ref) ;
        name_store("facecenter",faceCenter) ;
        name_store("BC_options",BC_options) ;
        name_store("Treservoir_f",Treservoir_f) ;
        input("facecenter,ref->BC_options") ;
        output("Treservoir_f") ;
        constraint("ref->Treservoir_functionSteady_BC") ;
      }

      // Calculate reservoir temperature for a single face.
      void calculate(Entity face) {
        map<string,real> varMap ; varMap["x"]=faceCenter[face].x ;
        varMap["y"]=faceCenter[face].y ; varMap["z"]=faceCenter[face].z ;
        string function ;
        BC_options[ref[face]].getOption("Treservoir",function) ;
        Loci::exprP p=Loci::expression::create(function) ;
        Treservoir_f[face]=p->evaluate(varMap) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryTreservoirFunctionSteady>
    registerBoundaryTreservoirFunctionSteady ;

//-----------------------------------------------------------------------------
// Rwall boundary condition rules.

  // Rule for specifying constant Rwall. Assigns value to all boundary
  // faces that have the property constantRwall_BC.
  class BoundaryRwallConstant : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<real> constantRwall_BC ;
      store<real> Rwall_f ;
    public:

      // Define input and output.
      BoundaryRwallConstant() {
        name_store("ref",ref) ;
        name_store("constantRwall_BC",constantRwall_BC) ;
        name_store("Rwall_f",Rwall_f) ;
        input("ref->constantRwall_BC") ;
        output("Rwall_f") ;
      }

      // Assign Rwall for a single face.
      void calculate(Entity face) {
        Rwall_f[face]=constantRwall_BC[ref[face]] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryRwallConstant> registerBoundaryRwallConstant ;

//-----------------------------------------------------------------------------
// Rules to create a constraint which indicates that the temperature is
// specified via some means.

  // Temperature specified with T=value or Twall=value.
  class BoundaryTemperatureConstantSpecification : public pointwise_rule {
    private:
      store<bool> specifiedTemperature_BC ;
    public:

      // Define input and output.
      BoundaryTemperatureConstantSpecification() {
        name_store("specifiedTemperature_BC",specifiedTemperature_BC) ;
        output("specifiedTemperature_BC") ;
        constraint("ref->T_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}

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

      // Do nothing.
      virtual void compute(const sequence &seq) {}

  } ;

  register_rule<BoundaryTemperatureCartesianSpecification>
    registerBoundaryTemperatureCartesianSpecification ;

  // Temperature specified indirectly via an energy balance between the
  // fluid and solid wall.
  class BoundaryTemperatureReservoirSpecification : public pointwise_rule {
    private:
      store<bool> specifiedTemperature_BC ;
    public:

      // Define input and output.
      BoundaryTemperatureReservoirSpecification() {
        name_store("specifiedTemperature_BC",specifiedTemperature_BC) ;
        output("specifiedTemperature_BC") ;
        constraint("ref->Treservoir_BCoption,ref->Rwall_BCoption") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}

  } ;

  register_rule<BoundaryTemperatureReservoirSpecification>
    registerBoundaryTemperatureReservoirSpecification ;

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

  // Compute temperature at boundary based on a specified wall heat flux.
  class BoundaryTemperatureFromHeatFlux : public pointwise_rule {
    private:
      const_Map ref,ci ;
      const_store<real> qwall_BC ;
      const_store<Area> area ;
      const_store<vect3d> cellCenter,faceCenter ;
      const_store<real> thermalConductivity,temperature ;
      store<real> temperature_f ;
    public:

      // Define input and output.
      BoundaryTemperatureFromHeatFlux() {
        name_store("ref",ref) ;
        name_store("ci",ci) ;
        name_store("area",area) ;
        name_store("qwall_BC",qwall_BC) ;
        name_store("cellcenter",cellCenter) ;
        name_store("temperature",temperature) ;
        name_store("kconduct(temperature,p,y)",thermalConductivity) ;
        name_store("facecenter",faceCenter) ;
        name_store("temperature_f",temperature_f) ;
        input("ci->(cellcenter,temperature,kconduct(temperature,p,y))") ;
        input("facecenter,area,ref->qwall_BC") ;
        output("temperature_f") ;
      }

      // Calculate temperature for a single face. Since we use the exact value
      // for qwall that was specified by the user in assembly of the energy
      // equation, this resulting wall value is only used for computation of the
      // cell temperature gradient as well as post-processing.
      void calculate(Entity face) {
        const real dx=dot((faceCenter[face]-cellCenter[ci[face]]),area[face].n) ;
        temperature_f[face]=-qwall_BC[ref[face]]*dx/thermalConductivity[ci[face]]+
          temperature[ci[face]] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryTemperatureFromHeatFlux>
    registerBoundaryTemperatureFromHeatFlux ;

  // Compute temperature at boundary from the wall energy balance between
  // the fluid and the solid. The wall heat flux for the fluid is later backed out
  // from this assigned wall temperature value. For consistency, it is imperative
  // that the formula used here and in the later rule for the heat flux from the
  // fluid side be the same, and it is. This rule is for constant Rwall, in which
  // case we do not need to do a Newton iteration.
  class BoundaryTemperatureFromWallEnergyBalance : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> Treservoir_f,Rwall_f ;
      const_store<Area> area ;
      const_store<vect3d> cellCenter,faceCenter ;
      const_store<real> thermalConductivity,temperature ;
      store<real> temperature_f ;
    public:

      // Define input and output.
      BoundaryTemperatureFromWallEnergyBalance() {
        name_store("ci",ci) ;
        name_store("area",area) ;
        name_store("Treservoir_f",Treservoir_f) ;
        name_store("Rwall_f",Rwall_f) ;
        name_store("cellcenter",cellCenter) ;
        name_store("temperature",temperature) ;
        name_store("kconduct(temperature,p,y)",thermalConductivity) ;
        name_store("facecenter",faceCenter) ;
        name_store("temperature_f",temperature_f) ;
        input("ci->(cellcenter,temperature,kconduct(temperature,p,y))") ;
        input("facecenter,area,Treservoir_f,Rwall_f") ;
        output("temperature_f") ;
        constraint("ref->Rwall_const_BC") ;
      }

      // Calculate temperature for a single face. Since we use the exact value
      // for qwall that was specified by the user in assembly of the energy
      // equation, this resulting wall value is only used for computation of the
      // cell temperature gradient as well as post-processing.
      void calculate(Entity face) {
        const real dy=dot((faceCenter[face]-cellCenter[ci[face]]),area[face].n),
          Rfluid=dy/thermalConductivity[ci[face]] ;
        temperature_f[face]=(Rwall_f[face]*temperature[ci[face]]+Rfluid*
          Treservoir_f[face])/(Rfluid+Rwall_f[face]) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryTemperatureFromWallEnergyBalance>
    registerBoundaryTemperatureFromWallEnergyBalance ;

  // Rule for specifying temperature with a profile. A single Cartesian
  // coordinate is used for the interpolation.
  class BoundaryTemperatureProfileCartesian : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<vect3d> faceCenter ;
      const_store<CartesianBoundaryCondition> cartesianT_BC ;
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
      }

      // Calculate temperature for all faces in sequence.
      virtual void compute(const sequence &seq) {

        // Create a map to organize ref values for the faces.
        std::map<int, Loci::entitySet> bcmap ;
        for(sequence::const_iterator si=seq.begin();si!=seq.end();++si){
          bcmap[ref[*si]]+=*si ;
        }

        // Loop through the map. Each map entry has differnet input file.
        std::map<int,Loci::entitySet>::iterator bci ;
        for(bci=bcmap.begin();bci!=bcmap.end();++bci) {

          // Get the boundary condition.
          const CartesianBoundaryCondition &c=cartesianT_BC[bci->first] ;

          // Search for location in profile and interpolate value.
          const vect3d first=c.Coord(0) ; sequence s=sequence(bci->second) ;
          for(sequence::const_iterator si=s.begin();si!=s.end();++si) {
            int face=*si ;
            const vect3d &fcenter=faceCenter[face] ;
            size_t id1=0 ; real d1=0.0 ;
            switch(c.CoordFlag()){
              case 0: d1=fabs(first.x-fcenter.x) ; break ;
              case 1: d1=fabs(first.y-fcenter.y) ; break ;
              case 2: d1=fabs(first.z-fcenter.z) ; break ;
            }

            // Find the nearest point.
            for(size_t i=1;i<c.NumPoints();++i) {
              real dis=0.0 ;
              switch(c.CoordFlag()){
                case 0: dis=fabs(fcenter.x-c.Coord(i).x) ; break ;
                case 1: dis=fabs(fcenter.y-c.Coord(i).y) ; break ;
                case 2: dis=fabs(fcenter.z-c.Coord(i).z) ; break ;
              }
              if(d1>dis){ d1=dis ; id1=i ; }
            }

            // Find the second nearest point
            size_t id2=0 ; real dd=0.0 ;
            switch(c.CoordFlag()){
              case 0: dd=c.Coord(id1).x-fcenter.x ; break ;
              case 1: dd=c.Coord(id1).y-fcenter.y ; break ;
              case 2: dd=c.Coord(id1).z-fcenter.z ; break ;
            }
            if(fabs(dd)>EPSILON)
              if(dd<0.0)
                if(id1<c.NumPoints()-1){
                  id2=id1+1 ;
                }else{
                  temperature_f[face]=c.Value(c.NumPoints()-1) ; continue ;
                }
              else
                if(id1>0){
                  id2=id1-1 ;
                }else{
                  temperature_f[face]=c.Value(0) ; continue ;
                }
            else {
              temperature_f[face]=c.Value(id1) ; continue ;
            }
            real d2=0.0 ;
            switch(c.CoordFlag()){
              case 0: d2=fabs(c.Coord(id2).x-fcenter.x) ; break ;
              case 1: d2=fabs(c.Coord(id2).y-fcenter.y) ; break ;
              case 2: d2=fabs(c.Coord(id2).z-fcenter.z) ; break ;
            }
            temperature_f[face]=(d2*c.Value(id1)+d1*c.Value(id2))/(d1+d2) ;
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
//      constraint("boundaryFaces,compressibleFlow") ; 1/28/2009
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
// Rules for computing the diffusion coefficient for total enthalpy.

  // Total enthalpy diffusion coefficient. This has been made parametric so we
  // can use it for both cells and boundary faces.
  class TotalEnthalpyDiffusionCoefficient : public pointwise_rule {
    private:
      const_param<real> turbulentPrandtlNumber ;
      const_store<real> T,p,eddyViscosity ;
      const_storeVec<real> y ;
      const_store<real> cp,thermalConductivity ;
      store<real> hDiffusionCoeff ;
    public:

      // Define input and output.
      TotalEnthalpyDiffusionCoefficient() {
        name_store("turbulentPrandtlNumber",turbulentPrandtlNumber) ;
        name_store("T",T) ;
        name_store("P",p) ;
        name_store("Y",y) ;
        name_store("V",eddyViscosity) ;
        name_store("cp",cp) ;
        name_store("kconduct(T,P,Y)",thermalConductivity) ;
        name_store("hDiffusionCoeff(T,P,Y,V)",hDiffusionCoeff) ;
        input("turbulentPrandtlNumber,T,P,Y,V,cp") ;
        input("kconduct(T,P,Y)") ;
        output("hDiffusionCoeff(T,P,Y,V)") ;
      }

      // Compute for a single entity.
      void calculate(Entity e) {
        hDiffusionCoeff[e]=thermalConductivity[e]/cp[e]+eddyViscosity[e]/
          (*turbulentPrandtlNumber) ;
      }

      // Compute for a sequence of entities.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TotalEnthalpyDiffusionCoefficient>
    registerTotalEnthalpyDiffusionCoefficient ;

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
        constraint("noslip_BC,ref->T_BC") ;
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
        name_store("qWall",qWall) ;
        input("ci->grads(temperature),area") ;
        input("tDiffusionCoeff(temperature_f,p_f,y_f,eddyViscosity_f)") ;
        output("qWall") ;
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
        name_store("qWall",qWall) ;
        input("ci->grads(temperature),area") ;
        input("tDiffusionCoeff(temperature_f,p_f,y_f,eddyViscosity_f)") ;
        output("qWall") ;
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
        name_store("qWall",qWall) ;
        output("qWall") ;
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
        name_store("qWall",qWall) ;
        output("qWall") ;
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
        name_store("qWall",qWall) ;
        output("qWall") ;
        constraint("viscousFlow,ref->adiabatic_BCoption") ;
      }
                                                                                
      // No heat flux.
      void calculate (Entity face) { qWall[face]=0.0 ; }
                                                                                
      // Loop over faces.
      virtual void compute (const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<QWallAdiabatic> registerQWallAdiabatic ;

  // Rule computes heat flux assuming no wall function. This rule provides
  // values for noslip walls where there is no wall function.
  class QWallNoSlip : public pointwise_rule {
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
      QWallNoSlip() {
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
        constraint("noslip_BC,specifiedTemperature_BC") ;
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

  register_rule<QWallNoSlip> registerQWallNoSlip ;

  // Priority rule for wall function boundaries.
  class QWallNoSlipWallFunction : public pointwise_rule {
    private:
      const_store<real> qWallTemp ;
      store<real> qWall ;
    public:
                                                                                
      // Define input and output.
      QWallNoSlipWallFunction() {
        name_store("qWallTemp",qWallTemp) ;
        name_store("wallFunction::qWall",qWall) ;
        input("qWallTemp") ;
        output("wallFunction::qWall") ;
        constraint("viscousFlow,ref->wallFunction_BCoption") ;
        constraint("noslip_BC,specifiedTemperature_BC") ;
      }
                                                                                
      // Copy temp value from wall function rule.
      void calculate (Entity face) { qWall[face]=qWallTemp[face] ; }
                                                                                
      // Loop over faces.
      virtual void compute (const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<QWallNoSlipWallFunction> registerQWallNoSlipWallFunction ;

  // Priority rule for specified heat flux at a boundary. The wall heat
  // flux specified by the user represents the flux coming to the wall from
  // the fluid domain. This corresponds to a flux out of the fluid domain.
  class QWallSpecified : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<real> qwall_BC ;
      store<real> qWall ;
    public:

      // Define input and output.
      QWallSpecified() {
        name_store("ref",ref) ;
        name_store("qwall_BC",qwall_BC) ;
        name_store("qWall",qWall) ;
        input("ref->qwall_BC") ;
        output("qWall") ;
      }

      // Set heat flux.
      void calculate (Entity face) { qWall[face]=qwall_BC[ref[face]] ; }

      // Loop over faces.
      virtual void compute (const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<QWallSpecified> registerQWallSpecified ;

//-----------------------------------------------------------------------------
// Scheme independent rules for assembling the energy equation.

  // Rule to initialize the main coefficient. Checked.
  class InitializeTotalEnthalpyMainCoefficient : public unit_rule {
    private:
      store<real> hMainCoefficient ;
    public:

      // Define input and output.
      InitializeTotalEnthalpyMainCoefficient() {
        name_store("hMainCoefficient",hMainCoefficient) ;
        output("hMainCoefficient") ;
        constraint("vol") ;
      }

      // Set the main coefficient to zero for a single cell.
      void calculate(Entity cell) { hMainCoefficient[cell]=0.0 ; }

      // Set the main coefficient to zero for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeTotalEnthalpyMainCoefficient>
    registerInitializeTotalEnthalpyMainCoefficient ;

  // Rule to add the first-order inviscid flux contribution to the main
  // coefficient for interior faces. Checked.
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
        name_store("massFlux",massFlux) ;
        name_store("hMainCoefficient",hMainCoefficient) ;
        input("massFlux") ;
        output("cl->hMainCoefficient,cr->hMainCoefficient") ;
        constraint("internalFaces") ;
      }

      // Increment the main coefficient for cells attached to a single face.
      void calculate(Entity face) {
        if(massFlux[face]>0.0){
          hMainCoefficient[cr[face]]+=massFlux[face] ;
        }else{
          hMainCoefficient[cl[face]]-=massFlux[face] ;
        }
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToTotalEnthalpyMainCoefficientInterior>
    registerFOUInviscidFluxToTotalEnthalpyMainCoefficientInterior ;

  // Rule to add the first-order inviscid flux contribution to the main
  // coefficient for boundary faces. Checked.
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
        name_store("massFlux",massFlux) ;
        name_store("hMainCoefficient",hMainCoefficient) ;
        input("massFlux") ;
        output("ci->hMainCoefficient") ;
        constraint("boundaryFaces") ;
      }

      // Increment the main coefficient for the cell attach to a single face.
      void calculate(Entity face) {
        if(massFlux[face]<0.0) hMainCoefficient[ci[face]]-=massFlux[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToTotalEnthalpyMainCoefficientBoundary>
    registerFOUInviscidFluxToTotalEnthalpyMainCoefficientBoundary ;

  // Rule to add the diffusive flux contribution to the main coefficient for
  // interior faces.
  class DiffusiveFluxToTotalEnthalpyMainCoefficientInterior : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> faceRadius ;
      const_store<real> diffusionProduct ;
      const_store<real> hDiffusionCoeff ;
      store<real> hMainCoefficient ;
    public:

      // Define input and output.
      DiffusiveFluxToTotalEnthalpyMainCoefficientInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("faceRadius",faceRadius) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("hDiffusionCoeff(temperature,p,y,eddyViscosity)",
          hDiffusionCoeff) ;
        name_store("hMainCoefficient",hMainCoefficient) ;
        input("faceRadius,diffusionProduct") ;
        input("(cl,cr)->hDiffusionCoeff(temperature,p,y,eddyViscosity)") ;
        output("(cl,cr)->hMainCoefficient") ;
        constraint("internalFaces,viscousFlow") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real temp=0.5*(hDiffusionCoeff[cl[face]]+hDiffusionCoeff[cr[face]])*
          diffusionProduct[face]*faceRadius[face] ;
        hMainCoefficient[cl[face]]+=temp ; hMainCoefficient[cr[face]]+=temp ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToTotalEnthalpyMainCoefficientInterior>
    registerDiffusiveFluxToTotalEnthalpyMainCoefficientInterior ;

  // Rule to ensure main coefficient is not zero for inviscid flows.
  class NetMassFluxToTotalEnthalpyMainCoefficient : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_store<real> netMassFlux ;
      store<real> hMainCoefficient ;
    public:
                                                                                
      // Define input and output.
      NetMassFluxToTotalEnthalpyMainCoefficient() {
        name_store("netMassFlux",netMassFlux) ;
        name_store("hMainCoefficient",hMainCoefficient) ;
        input("netMassFlux") ;
        output("hMainCoefficient") ;
        constraint("inviscidFlow,geom_cells") ;
      }
                                                                                
      // Add net mass flux for a single cell.
      void calculate(Entity cell) {
        if(netMassFlux[cell]<0.0) hMainCoefficient[cell]-=netMassFlux[cell] ;
      }
                                                                                
      // Add net mass flux for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<NetMassFluxToTotalEnthalpyMainCoefficient>
    registerNetMassFluxToTotalEnthalpyMainCoefficient ;

  // Rule to initialize the source term. Checked.
  class InitializeTotalEnthalpySourceTerm : public unit_rule {
    private:
      store<real> hSourceTerm ;
    public:

      // Define input and output.
      InitializeTotalEnthalpySourceTerm() {
        name_store("hSourceTerm",hSourceTerm) ;
        output("hSourceTerm") ;
        constraint("vol") ;
      }

      // Set the source term to zero for a single cell.
      void calculate(Entity cell) { hSourceTerm[cell]=0.0 ; }

      // Set the source term to zero for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeTotalEnthalpySourceTerm>
    registerInitializeTotalEnthalpySourceTerm ;

  // Rule to add the first-order inviscid flux contribution to the source term
  // for boundary faces. Checked.
  class FOUInviscidFluxToTotalEnthalpySourceTermBoundary : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> h ;
      const_store<real> massFlux ;
      const_store<real> h_f ;
      store<real> hSourceTerm ;
    public:

      // Define input and output.
      FOUInviscidFluxToTotalEnthalpySourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("h",h) ;
        name_store("massFlux",massFlux) ;
        name_store("h_f",h_f) ;
        name_store("hSourceTerm",hSourceTerm) ;
        input("ci->h") ;
        input("massFlux,h_f") ;
        output("ci->hSourceTerm") ;
        constraint("boundaryFaces,fouInviscidFlux") ;
      }

      // Increment the source term for the cell attach to a single face.
      void calculate(Entity face) {
        if(massFlux[face]<=0.0) hSourceTerm[ci[face]]-=massFlux[face]*
          h_f[face] ;
        else hSourceTerm[ci[face]]-=massFlux[face]*(h_f[face]-h[ci[face]]) ;
      }

      // Call calculate for a sequence of faces.
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
        name_store("massFlux",massFlux) ;
        name_store("hSourceTerm",hSourceTerm) ;
        input("(cl,cr)->(hGradient,limiters(h),cellcenter)") ;
        input("facecenter,massFlux") ;
        output("(cl,cr)->hSourceTerm") ;
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

      // Call calculate for a sequence of faces.
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
        name_store("massFlux",massFlux) ;
        name_store("h_f",h_f) ;
        name_store("hSourceTerm",hSourceTerm) ;
        input("ci->(hGradient,limiters(h),cellcenter)") ;
        input("facecenter,massFlux,h_f") ;
        output("ci->hSourceTerm") ;
        constraint("boundaryFaces,souOrRoeInviscidFlux") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        hSourceTerm[ci[face]]-=(massFlux[face]>0.0)? massFlux[face]*
          hLimiter[ci[face]]*dot(hGradient[ci[face]],faceCenter[face]-
          cellCenter[ci[face]]):massFlux[face]*h_f[face] ;
      }

      // Call calculate for a sequence of faces.
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
        name_store("hSourceTerm",hSourceTerm) ;
        input("(cl,cr)->(vol,grads(temperature))") ;
        input("(cl,cr)->tDiffusionCoeff(temperature,p,y,eddyViscosity)") ;
        input("area,faceRadius") ;
        output("(cl,cr)->hSourceTerm") ;
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

      // Call calculate for a sequence of faces.
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
        name_store("hSourceTerm",hSourceTerm) ;
        input("qWall,area,faceRadius") ;
        output("ci->hSourceTerm") ;
        constraint("viscousFlow,boundaryFaces") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        hSourceTerm[ci[face]]-=qWall[face]*area[face].sada*faceRadius[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<HeatFluxToTotalEnthalpySourceTermBoundary>
    registerHeatFluxToTotalEnthalpySourceTermBoundary ;

  // Rule to subtract the implicit diffusion term from the source term.
  class DiffusiveFluxToTotalEnthalpySourceTermInteriorEnergy3 : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> faceRadius ;
      const_store<real> diffusionProduct ;
      const_store<real> hDiffusionCoeff ;
      const_store<real> h ;
      store<real> hSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToTotalEnthalpySourceTermInteriorEnergy3() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("faceRadius",faceRadius) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("hDiffusionCoeff(temperature,p,y,eddyViscosity)",
          hDiffusionCoeff) ;
        name_store("h",h) ;
        name_store("hSourceTerm",hSourceTerm) ;
        input("faceRadius,diffusionProduct") ;
        input("(cl,cr)->hDiffusionCoeff(temperature,p,y,eddyViscosity)") ;
        input("(cl,cr)->h") ;
        output("(cl,cr)->hSourceTerm") ;
        constraint("internalFaces,viscousFlow") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        real temp=0.5*(hDiffusionCoeff[cl[face]]+hDiffusionCoeff[cr[face]])*
          diffusionProduct[face]*faceRadius[face] ;
        hSourceTerm[cl[face]]+=temp*(h[cl[face]]-h[cr[face]]) ;
        hSourceTerm[cr[face]]+=temp*(h[cr[face]]-h[cl[face]]) ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToTotalEnthalpySourceTermInteriorEnergy3>
    registerDiffusiveFluxToTotalEnthalpySourceTermInteriorEnergy3 ;

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
        name_store("v",v) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("viscosity",viscosity) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("hSourceTerm",hSourceTerm) ;
        input("(cl,cr)->(vol,v,gradv3d(v)),viscosity") ;
        input("area,faceRadius") ;
        output("(cl,cr)->hSourceTerm") ;
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

      // Call calculate for a sequence of faces.
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
        name_store("gradv3d(v)",vGradient) ;
        name_store("v_f",v_f) ;
        name_store("viscosity",viscosity) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("hSourceTerm",hSourceTerm) ;
        input("ci->gradv3d(v),v_f,viscosity,area,faceRadius") ;
        output("ci->hSourceTerm") ;
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

      // Call calculate for a sequence of faces.
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
        name_store("hSourceTerm",hSourceTerm) ;
        input("delP,deltaP,vNormalTilde,area") ;
        output("(cl,cr)->hSourceTerm") ;
        constraint("internalFaces,roeInviscidFlux") ;
      }
                                                                                
      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        real source=0.5*(vNormalTilde[face]*delP[face]-abs(vNormalTilde[face])*
          deltaP[face])*area[face].sada ;
        hSourceTerm[cl[face]]+=source ; hSourceTerm[cr[face]]-=source ;
      }
                                                                                
      // Call calculate for a sequence of faces.
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
        name_store("v",v) ;
        name_store("vol",vol) ;
        name_store("hSourceTerm",hSourceTerm) ;
        input("gravityAcceleration,gravityRhoRef,rho,v,vol") ;
        output("hSourceTerm") ;
        constraint("geom_cells") ;
      }
                                                                                
      // Add buoyance force for single cell.
      void calculate(Entity cell) {
        hSourceTerm[cell]+=(rho[cell]-*rhoReference)*dot(*gravity,v[cell])*
          vol[cell] ;
      }
                                                                                
      // Call calculate for a sequence of cells.
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
        name_store("hSourceTerm",hSourceTerm) ;
        input("hDiffusion") ;
        output("hSourceTerm") ;
        constraint("speciesTransport,viscousFlow,geom_cells") ;
      }
                                                                                
      // Add species enthalpy diffusion for single cell.
      void calculate(Entity cell) {
        hSourceTerm[cell]+=hDiffusion[cell] ;
      }
                                                                                
      // Call calculate for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<SpeciesDiffusionToTotalEnthalpySourceTerm>
    registerSpeciesDiffusionToTotalEnthalpySourceTerm ;

  // Rule to compute the diagonal term for the linear system. Checked.
  class ComputeTotalEnthalpyMatrixDiagonal : public pointwise_rule {
    private:
      const_param<real> hRelaxationFactor ;
      const_store<real> hMainCoefficient ;
      store<real> D ;
    public:

      // Define input and output.
      ComputeTotalEnthalpyMatrixDiagonal() {
        name_store("hStar_RelaxationFactor",hRelaxationFactor) ;
        name_store("hMainCoefficient",hMainCoefficient) ;
        name_store("hStar_D",D) ;
        input("hStar_RelaxationFactor,hMainCoefficient") ;
        output("hStar_D") ;
        constraint("geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) {
        D[cell]=hMainCoefficient[cell]/(*hRelaxationFactor) ;
      }

      // Add relaxation for a sequence of cells.
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

  // Rule to initialize the lower terms for the linear system. Checked.
  class InitializeTotalEnthalpyMatrixLower : public unit_rule {
    private:
      store<real> L ;
    public:

      // Define input and output.
      InitializeTotalEnthalpyMatrixLower() {
        name_store("hStar_L",L) ;
        output("hStar_L") ;
        constraint("internalFaces") ;
      }

      // Initialize for a single face.
      void calculate(Entity face) { L[face]=0.0 ; }

      // Initialize for a sequence of faces.
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
        name_store("massFlux",massFlux) ;
        name_store("hStar_L",L) ;
        input("massFlux") ;
        output("hStar_L") ;
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

  register_rule<FOUInviscidFluxToTotalEnthalpyMatrixLower>
    registerFOUInviscidFluxToTotalEnthalpyMatrixLower ;

  // Rule to add the diffusive flux contribution to the lower terms for the
  // linear system.
  class DiffusiveFluxToTotalEnthalpyMatrixLower : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> hDiffusionCoeff ;
      const_store<real> diffusionProduct ;
      const_store<real> faceRadius ;
      store<real> L ;
    public:

      // Define input and output.
      DiffusiveFluxToTotalEnthalpyMatrixLower() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("hDiffusionCoeff(temperature,p,y,eddyViscosity)",
          hDiffusionCoeff) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("faceRadius",faceRadius) ;
        name_store("hStar_L",L) ;
        input("(cl,cr)->hDiffusionCoeff(temperature,p,y,eddyViscosity)") ;
        input("diffusionProduct,faceRadius") ;
        output("hStar_L") ;
        constraint("internalFaces,viscousFlow") ;
      }

      // Increment the lower term for a single face.
      void calculate(Entity face) {
        L[face]-=0.5*(hDiffusionCoeff[cl[face]]+hDiffusionCoeff[cr[face]])*
          diffusionProduct[face]*faceRadius[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToTotalEnthalpyMatrixLower>
    registerDiffusiveFluxToTotalEnthalpyMatrixLower ;

  // Rule to initialize the upper terms for the linear system. Checked.
  class InitializeTotalEnthalpyMatrixUpper : public unit_rule {
    private:
      store<real> U ;
    public:

      // Define input and output.
      InitializeTotalEnthalpyMatrixUpper() {
        name_store("hStar_U",U) ;
        output("hStar_U") ;
        constraint("internalFaces") ;
      }

      // Initialize for a single face.
      void calculate(Entity face) { U[face]=0.0 ; }

      // Initialize for a sequence of faces.
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
        name_store("massFlux",massFlux) ;
        name_store("hStar_U",U) ;
        input("massFlux") ;
        output("hStar_U") ;
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

  register_rule<FOUInviscidFluxToTotalEnthalpyMatrixUpper>
    registerFOUInviscidFluxToTotalEnthalpyMatrixUpper ;

  // Rule to add the diffusive flux contribution to the upper terms for the
  // linear system.
  class DiffusiveFluxToTotalEnthalpyMatrixUpper : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> hDiffusionCoeff ;
      const_store<real> diffusionProduct ;
      const_store<real> faceRadius ;
      store<real> U ;
    public:

      // Define input and output.
      DiffusiveFluxToTotalEnthalpyMatrixUpper() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("hDiffusionCoeff(temperature,p,y,eddyViscosity)",
          hDiffusionCoeff) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("faceRadius",faceRadius) ;
        name_store("hStar_U",U) ;
        input("(cl,cr)->hDiffusionCoeff(temperature,p,y,eddyViscosity)") ;
        input("diffusionProduct,faceRadius") ;
        output("hStar_U") ;
        constraint("internalFaces,viscousFlow") ;
      }

      // Increment the lower term for a single face.
      void calculate(Entity face) {
        U[face]-=0.5*(hDiffusionCoeff[cl[face]]+hDiffusionCoeff[cr[face]])*
          diffusionProduct[face]*faceRadius[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToTotalEnthalpyMatrixUpper>
    registerDiffusiveFluxToTotalEnthalpyMatrixUpper ;

  // Rule to compute the right-hand side for the linear system. Checked.
  class ComputeTotalEnthalpyRHS : public pointwise_rule {
    private:
      const_param<real> hRelaxationFactor ;
      const_store<real> h ;
      const_store<real> hMainCoefficient ;
      const_store<real> hSourceTerm ;
      store<real> B ;
    public:

      // Define input and output.
      ComputeTotalEnthalpyRHS() {
        name_store("hStar_RelaxationFactor",hRelaxationFactor) ;
        name_store("h",h) ;
        name_store("hMainCoefficient",hMainCoefficient) ;
        name_store("hSourceTerm",hSourceTerm) ;
        name_store("hStar_B",B) ;
        input("hStar_RelaxationFactor,h,hMainCoefficient,hSourceTerm") ;
        output("hStar_B") ;
        constraint("geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) {
        B[cell]=hSourceTerm[cell]+(1.0-(*hRelaxationFactor))*
          hMainCoefficient[cell]*h[cell]/(*hRelaxationFactor) ;
      }

      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeTotalEnthalpyRHS> registerComputeTotalEnthalpyRHS ;

  // Priority rule for inviscid flow to compute the right-hand side for the
  // linear system. Includes the net mass flux term which cancels the term
  // added to the main coefficient. Checked.
  class ComputeTotalEnthalpyRHSInviscid : public pointwise_rule {
    private:
      const_param<real> hRelaxationFactor ;
      const_store<real> h ;
      const_store<real> hMainCoefficient ;
      const_store<real> hSourceTerm ;
      const_store<real> netMassFlux ;
      store<real> B ;
    public:
                                                                                
      // Define input and output.
      ComputeTotalEnthalpyRHSInviscid() {
        name_store("hStar_RelaxationFactor",hRelaxationFactor) ;
        name_store("h",h) ;
        name_store("hMainCoefficient",hMainCoefficient) ;
        name_store("hSourceTerm",hSourceTerm) ;
        name_store("netMassFlux",netMassFlux) ;
        name_store("inviscidFlow::hStar_B",B) ;
        input("hStar_RelaxationFactor,h,hMainCoefficient") ;
        input("hSourceTerm,netMassFlux") ;
        output("inviscidFlow::hStar_B") ;
        constraint("inviscidFlow,geom_cells") ;
      }
                                                                                
      // Add relaxation for a single cell.
      void calculate(Entity cell) {
        B[cell]=hSourceTerm[cell]+(1.0-(*hRelaxationFactor))*
          hMainCoefficient[cell]*h[cell]/(*hRelaxationFactor) ;
        if(netMassFlux[cell]<0.0) B[cell]-=netMassFlux[cell]*h[cell] ;
      }
                                                                                
      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<ComputeTotalEnthalpyRHSInviscid>
    registerComputeTotalEnthalpyRHSInviscid ;

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

  register_rule<InitializeEnergyResidual> registerInitializeEnergyResidual ;

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

  register_rule<ComputeEnergyResidualOne> registerComputeEnergyResidualOne ;

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

  register_rule<ComputeEnergyResidualTwo> registerComputeEnergyResidualTwo ;

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

  register_rule<InitializeTotalEnergyResidual>
    registerInitializeTotalEnergyResidual ;

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

  register_rule<ComputeTotalEnergyResidual>
    registerComputeTotalEnergyResidual ;

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
      const_store<real> temperatureTimeStepN ;
      store<real> temperatureIterationZero ;
    public:

      // Define input and output.
      IterationBuildTemperature() {
        name_store("temperature{n}",temperatureTimeStepN) ;
        name_store("temperature{n,it=0}",temperatureIterationZero) ;
        input("temperature{n}") ;
        output("temperature{n,it=0}") ;
        constraint("compressibleFlow{n},geom_cells{n}") ;
      }

      // Assign temperature at iteration zero for a single cell.
      void calculate(Entity cell) {
        temperatureIterationZero[cell]=temperatureTimeStepN[cell] ;
      }

      // Assign temperature at iteration zero for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationBuildTemperature> registerIterationBuildTemperature ;

  // Iteration advance rule for temperature. TEMPORARY.
  class IterationAdvanceTemperature : public pointwise_rule {
    private:
      const_store<real> temperatureStar ;
      store<real> temperatureIterationPlusOne ;
    public:

      // Define input and output.
      IterationAdvanceTemperature() {
        name_store("temperatureStar{n,it}",temperatureStar) ;
        name_store("temperature{n,it+1}",temperatureIterationPlusOne) ;
        input("temperatureStar{n,it}") ;
        output("temperature{n,it+1}") ;
        constraint("geom_cells{n,it},compressibleFlow{n,it}") ;
      }

      // Assign temperature at end of iteration for a single cell.
      void calculate(Entity cell) {
        temperatureIterationPlusOne[cell]=temperatureStar[cell] ;
      }

      // Assign temperature at end of iteration for a sequence of cells.
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
      const_store<real> hTimeStepN ;
      store<real> hIterationZero ;
    public:

      // Define input and output.
      IterationBuildTotalEnthalpy() {
        name_store("h{n}",hTimeStepN) ;
        name_store("h{n,it=0}",hIterationZero) ;
        input("h{n}") ;
        output("h{n,it=0}") ;
        constraint("geom_cells{n},compressibleFlow{n}") ;
      }

      // Assign total enthalpy at iteration zero for a single cell.
      void calculate(Entity cell) {
        hIterationZero[cell]=hTimeStepN[cell] ;
      }

      // Assign total enthalpy at iteration zero for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationBuildTotalEnthalpy>
    registerIterationBuildTotalEnthalpy ;

  // Rule to add temporal component of the total enthalpy equation to the main
  // coefficient.
  class TemporalToTotalEnthalpyMainCoefficient: public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> timeIntegratorFactor0 ;
      const_store<real> timeStepFactor ;
      const_store<real> rho ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> hMainCoefficient ;
    public:

      // Define input and output.
      TemporalToTotalEnthalpyMainCoefficient() {
        name_store("timeIntegratorFactor0{n}",timeIntegratorFactor0) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n}",rho) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("hMainCoefficient{n,it}",hMainCoefficient) ;
        input("rho{n},vol{n},cellRadius{n,it}") ;
        input("timeStepFactor{n},timeIntegratorFactor0{n}") ;
        output("hMainCoefficient{n,it}") ;
        constraint("geom_cells") ;
      }

      // Add temporal component for a single cell.
      void calculate(Entity cell) {
        hMainCoefficient[cell]+=rho[cell]*vol[cell]*cellRadius[cell]*
          (*timeIntegratorFactor0)/timeStepFactor[cell] ;
      }

      // Add temporal component for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToTotalEnthalpyMainCoefficient>
    registerTemporalToTotalEnthalpyMainCoefficient ;

  // Rule to add temporal component of the total enthalpy equation to the
  // source term. Note that we still need dt{n} in here since we have the
  // temporal pressure derivative on the rhs.
  class TemporalToTotalEnthalpySourceTerm : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> dt ;
      const_param<real> timeIntegratorFactor0 ;
      const_store<real> timeStepFactor ;
      const_store<real> rho ;
      const_store<real> pNew,pOld ;
      const_store<real> h ;
      const_store<real> oldVol,vol ;
      const_store<real> cellRadius ;
      store<real> hSourceTerm ;
    public:

      // Define input and output.
      TemporalToTotalEnthalpySourceTerm() {
        name_store("dt{n}",dt) ;
        name_store("timeIntegratorFactor0{n}",timeIntegratorFactor0) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n}",rho) ;
        name_store("p{n}",pOld) ;
        name_store("p{n,it}",pNew) ;
        name_store("h{n}",h) ;
        name_store("vol{n}",oldVol) ;
        name_store("vol{n,it}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("hSourceTerm{n,it}",hSourceTerm) ;
        input("rho{n},p{n},p{n,it},h{n},vol{n},vol{n,it},cellRadius{n,it}") ;
        input("dt{n},timeStepFactor{n},timeIntegratorFactor0{n}") ;
        output("hSourceTerm{n,it}") ;
        constraint("geom_cells") ;
      }

      // Add temporal component to source term for a single cell.
      void calculate(Entity cell) {
        hSourceTerm[cell]+=(rho[cell]*h[cell]*oldVol[cell]*
          (*timeIntegratorFactor0)+(pNew[cell]-pOld[cell])*vol[cell]/
          (*dt))*cellRadius[cell]/timeStepFactor[cell] ;
      }

      // Add temporal component to source term for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToTotalEnthalpySourceTerm>
    registerTemporalToTotalEnthalpySourceTerm ;

  // Rule to add temporal component of total enthalpy equation to the
  // source term for the BDF2 scheme.
  class TemporalToTotalEnthalpySourceTermBDF2 : public apply_rule
  <store<real>,Loci::Summation<real> > {
    private:
      const_param<real> dtOld,dt ;
      const_param<real> timeIntegratorFactor1 ;
      const_param<real> timeIntegratorFactor2 ;
      const_store<real> timeStepFactor ;
      const_store<real> rhoOld ;
      const_store<real> pNew,pOld,pOlder ;
      const_store<real> hOld,h ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> hSourceTerm ;
    private:
      real f0,f1,f2 ;
    public:

      // Define input and output.
      TemporalToTotalEnthalpySourceTermBDF2() {
        name_store("dt{n-1}",dtOld) ;
        name_store("dt{n}",dt) ;
        name_store("timeIntegratorFactor1{n}",timeIntegratorFactor1) ;
        name_store("timeIntegratorFactor2{n}",timeIntegratorFactor2) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n-1}",rhoOld) ;
        name_store("p{n-1}",pOlder) ;
        name_store("p{n}",pOld) ;
        name_store("p{n,it}",pNew) ;
        name_store("h{n-1}",hOld) ;
        name_store("h{n,it}",h) ;
        name_store("vol{n,it}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("hSourceTerm{n,it}",hSourceTerm) ;
        input("rho{n-1},p{n-1},p{n},p{n,it},h{n-1},h{n,it}") ;
        input("vol{n,it},dt{n-1},dt{n},cellRadius{n,it}") ;
        input("timeStepFactor{n},timeIntegratorFactor1{n}") ;
        input("timeIntegratorFactor2{n}") ;
        output("hSourceTerm{n,it}") ;
        constraint("geom_cells,BDF2Integrator") ;
      }

      // Add temporal component to source term for a single cell. Checked
      // correctness of pressure term assuming constant dt. JW 3/3/2010
      void calculate(Entity cell) {
        hSourceTerm[cell]+=((*timeIntegratorFactor1)*rhoOld[cell]*
          vol[cell]*cellRadius[cell]/timeStepFactor[cell])*
          (h[cell]-hOld[cell]) ;
        hSourceTerm[cell]+=(f2*pNew[cell]-f1*pOld[cell]+f0*pOlder[cell])*
          vol[cell]*cellRadius[cell]/timeStepFactor[cell]*
          (*timeIntegratorFactor2) ;
      }

      // Loop over cells.
      void compute(const sequence &seq) {
        real P=(*dt)+(*dtOld),Q=P*P/((*dt)*(*dt)) ;
        f0=1.0/((*dt)*Q-P),f1=Q*f0-(1.0/(*dt)),f2=f1-f0 ;
        do_loop(seq,this) ;
      }
  } ;

  register_rule<TemporalToTotalEnthalpySourceTermBDF2>
    registerTemporalToTotalEnthalpySourceTermBDF2 ;

  // Iteration advance rule for total enthalpy.
  class IterationAdvanceTotalEnthalpy : public pointwise_rule {
    private:
      const_store<real> hStar ;
      store<real> hIterationPlusOne ;
    public:

      // Define input and output.
      IterationAdvanceTotalEnthalpy() {
        name_store("hStar{n,it}",hStar) ;
        name_store("h{n,it+1}",hIterationPlusOne) ;
        input("hStar{n,it}") ;
        output("h{n,it+1}") ;
        constraint("geom_cells{n,it},compressibleFlow{n,it}") ;
      }

      // Assign total enthalpy at end of iteration for a single cell.
      void calculate(Entity cell) {
        hIterationPlusOne[cell]=hStar[cell] ;
      }

      // Assign total enthalpy at end of iteration for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationAdvanceTotalEnthalpy>
    registerIterationAdvanceTotalEnthalpy ;

  // Iteration collapse rule for total enthalpy. Note that we have changed
  // this rule from "h{n+1}=h{n,it}" to the current form. With the old form
  // we could not do one iteration per time step, because h{n,it} would not
  // get updated and thus the residuals would never change. NOTE: The above
  // comment is now invalid, JW 12/7/2006.
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

