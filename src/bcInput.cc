// Standard library includes.
#include <set>
#include <vector>
using std::set ;
using std::vector ;

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;

// StreamUns includes.
#include "bcInput.h"
#include "sciTypes.h"

namespace streamUns {

  class BC_rho_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<real> rho_BC ;
    public:

      BC_rho_compute() {
        name_store("BC_options",BC_options) ;
        name_store("rho_BC",rho_BC) ;
        input("BC_options") ;
        output("rho_BC") ;
        constraint("rho_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOptionUnits("rho","kg/m/m/m",rho_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_rho_compute> register_BC_rho_compute ;

  // Creates the new velocity boundary condition constraints.
  class VelocityBoundaryConditionConstraints : public constraint_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      Constraint v_axisymmetric_BC,v_cartesian_BC,v_const_BC ;
      Constraint v_function_BC,v_functionSteady_BC,v_functionUnsteady_BC ;
      Constraint timeDependent_BC,timeIndependent_BC ;
    private:
      entitySet vAxisymmetric,vCartesian,vConstant ;
      entitySet vFunction,vFunctionSteady,vFunctionUnsteady ;
    public:

      VelocityBoundaryConditionConstraints() {
        name_store("BC_options",BC_options) ;
        name_store("v_axisymmetric_BC",v_axisymmetric_BC) ;
        name_store("v_cartesian_BC",v_cartesian_BC) ;
        name_store("v_const_BC",v_const_BC) ;
        name_store("v_function_BC",v_function_BC) ;
        name_store("v_functionSteady_BC",v_functionSteady_BC) ;
        name_store("v_functionUnsteady_BC",v_functionUnsteady_BC) ;
        name_store("timeDependent_BC",timeDependent_BC) ;
        name_store("timeIndependent_BC",timeIndependent_BC) ;
        input("BC_options") ;
        output("v_axisymmetric_BC,v_cartesian_BC,v_const_BC") ;
        output("v_function_BC,v_functionSteady_BC,v_functionUnsteady_BC") ;
        output("timeDependent_BC,timeIndependent_BC") ;
        constraint("v_BCoption") ;
      }

      // Note that we need to include UNIT_VALUE in the if check below
      // since including units on a scalar value for velocity causes the
      // option value type to be set to UNIT_VALUE instead of REAL.
      void calculate(Entity e) {
        Loci::option_value_type optionValueType=BC_options[e].
          getOptionValueType("v") ;
        if(optionValueType==Loci::REAL || optionValueType==Loci::UNIT_VALUE ||
        optionValueType==Loci::LIST){
          vConstant+=e ;
        }else if(optionValueType==Loci::FUNCTION){
          Loci::options_list::arg_list options ; string name ;
          BC_options[e].getOption("v",name,options) ;
          if(name=="axisymmetric"){
            vAxisymmetric+=e ;
          }else if(name=="cartesian"){
            vCartesian+=e ;
          }else if(name=="function"){
            bool isUnsteady=false ; string validOptions="vX:vY:vZ" ;
            Loci::options_list functionOptions(validOptions) ;
            functionOptions.Input(options) ;
            vector<string> optionName(3) ;
            optionName[0]="vX",optionName[1]="vY",optionName[2]="vZ" ;
            for(unsigned int i=0;i<3;++i){
              if(functionOptions.optionExists(optionName[i])){
                Loci::option_value_type optionValueType=functionOptions.
                  getOptionValueType(optionName[i]) ;
                if(optionValueType==Loci::STRING){
                  string value ;
                  functionOptions.getOption(optionName[i],value) ;
                  Loci::exprP p=Loci::expression::create(value) ;
                  set<string> namelist ; Loci::getVarNames(p,namelist) ;
                  if(namelist.find(string("t")) != namelist.end())
                    isUnsteady=true ;
                }else{
                  cerr << "ERROR: Bad option type for function value "
                    << optionName[i] << " !" ; Loci::Abort() ;
                }
              }
            }
            if(isUnsteady) vFunctionUnsteady+=e ; else vFunctionSteady+=e ;
            vFunction+=e ;
          }else{
            cerr << "ERROR: Velocity boundary condition /'" << name
              << "/' unknown!" << endl ; Loci::Abort() ;
          }
        }
      }

      void compute(const sequence &seq) {
        do_loop(seq,this) ;
        v_axisymmetric_BC=vAxisymmetric ; v_cartesian_BC=vCartesian ;
        v_const_BC=vConstant ; v_function_BC=vFunction ;
        v_functionSteady_BC=vFunctionSteady ;
        v_functionUnsteady_BC=vFunctionUnsteady ;
        timeDependent_BC=vFunctionUnsteady ;
        timeIndependent_BC=vAxisymmetric+vCartesian+vConstant+vFunctionSteady ;
      }
  } ;

  register_rule<VelocityBoundaryConditionConstraints>
    registerVelocityBoundaryConditionConstraints ;

  class ConstantVBCCompute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<vect3d> constantV_BC ;
    public:

      ConstantVBCCompute() {
        name_store("BC_options",BC_options) ;
        name_store("constantV_BC",constantV_BC) ;
        input("BC_options") ;
        output("constantV_BC") ;
        constraint("v_const_BC") ;
      }

      void calculate(Entity e) {
        get_vect3dOption(BC_options[e],"v","m/s",constantV_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ConstantVBCCompute> registerConstantVBCCompute ;

  class CartesianVBCCompute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<string> cartesianV_BC ;
    public:

      CartesianVBCCompute() {
        name_store("BC_options",BC_options) ;
        name_store("cartesianV_BC",cartesianV_BC) ;
        input("BC_options") ;
        output("cartesianV_BC") ;
        constraint("v_cartesian_BC") ;
      }

      void calculate(Entity e) {
        Loci::options_list::arg_list options ; string name ;
        BC_options[e].getOption("v",name,options) ;
        string validOptions="file" ;
        Loci::options_list cartesianOptions(validOptions) ;
        cartesianOptions.Input(options) ;
        Loci::option_value_type optionValueType=cartesianOptions.
          getOptionValueType("file") ;
        if(optionValueType==Loci::STRING){
          cartesianOptions.getOption("file",cartesianV_BC[e]) ;
        }else{
          cerr << "ERROR: Bad option type for Cartesian velocity "
            << "option 'file' in boundary conditions." << endl ;
        }
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CartesianVBCCompute> registerCartesianVBCCompute ;

  class AxisymmetricVBCCompute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<string> axisymmetricV_BC ;
    public:

      AxisymmetricVBCCompute() {
        name_store("BC_options",BC_options) ;
        name_store("axisymmetricV_BC",axisymmetricV_BC) ;
        input("BC_options") ;
        output("axisymmetricV_BC") ;
        constraint("v_axisymmetric_BC") ;
      }

      void calculate(Entity e) {
        Loci::options_list::arg_list options ; string name ;
        BC_options[e].getOption("v",name,options) ;
        string validOptions="file" ;
        Loci::options_list axisymmetricOptions(validOptions) ;
        axisymmetricOptions.Input(options) ;
        Loci::option_value_type optionValueType=axisymmetricOptions.
          getOptionValueType("file") ;
        if(optionValueType==Loci::STRING){
          axisymmetricOptions.getOption("file",axisymmetricV_BC[e]) ;
        }else{
          cerr << "ERROR: Bad option type for axisymmetric velocity "
            << "option 'file' in boundary conditions." << endl ;
        }
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<AxisymmetricVBCCompute> registerAxisymmetricVBCCompute ;

  // Stores the function strings for the velocity components for later use.
  // If the user does not specify a component, set it to the zero string.
  class FunctionVBCCompute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<string> functionVX_BC,functionVY_BC,functionVZ_BC ;
    public:

      FunctionVBCCompute() {
        name_store("BC_options",BC_options) ;
        name_store("functionVX_BC",functionVX_BC) ;
        name_store("functionVY_BC",functionVY_BC) ;
        name_store("functionVZ_BC",functionVZ_BC) ;
        input("BC_options") ;
        output("functionVX_BC,functionVY_BC,functionVZ_BC") ;
        constraint("v_function_BC") ;
      }

      void calculate(Entity e) {
        Loci::options_list::arg_list options ; string name ;
        BC_options[e].getOption("v",name,options) ;
        string validOptions="vX:vY:vZ" ;
        Loci::options_list functionOptions(validOptions) ;
        functionOptions.Input(options) ;
        if(functionOptions.optionExists("vX"))
          functionOptions.getOption("vX",functionVX_BC[e]) ;
        else
          functionVX_BC[e]="0.0" ;
        if(functionOptions.optionExists("vY"))
          functionOptions.getOption("vY",functionVY_BC[e]) ;
        else
          functionVY_BC[e]="0.0" ;
        if(functionOptions.optionExists("vZ"))
          functionOptions.getOption("vZ",functionVZ_BC[e]) ;
        else
          functionVZ_BC[e]="0.0" ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FunctionVBCCompute> registerFunctionVBCCompute ;

  class BC_p_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<real> p_BC ;
    public:

      BC_p_compute() {
        name_store("BC_options",BC_options) ;
        name_store("p_BC",p_BC) ;
        input("BC_options") ;
        output("p_BC") ;
        constraint("p_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOptionUnits("p","Pa",p_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_p_compute> register_BC_p_compute ;

  class BC_pMean_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<real> pMean_BC ;
    public:

      BC_pMean_compute() {
        name_store("BC_options",BC_options) ;
        name_store("pMean_BC",pMean_BC) ;
        input("BC_options") ;
        output("pMean_BC") ;
        constraint("pMean_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOptionUnits("pMean","Pa",pMean_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_pMean_compute> register_BC_pMean_compute ;

  class BC_p0_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<real> p0_BC ;
    public:

      BC_p0_compute() {
        name_store("BC_options",BC_options) ;
        name_store("p0_BC",p0_BC) ;
        input("BC_options") ;
        output("p0_BC") ;
        constraint("p0_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOptionUnits("p0","Pa",p0_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_p0_compute> register_BC_p0_compute ;

  class BC_T_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<real> T_BC ;
    public:

      BC_T_compute() {
        name_store("BC_options",BC_options) ;
        name_store("T_BC",T_BC) ;
        input("BC_options") ;
        output("T_BC") ;
        constraint("T_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOptionUnits("T","kelvin",T_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_T_compute> register_BC_T_compute ;

  // Support for Twall for compatibility with CHEM.
  class BC_Twall_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<real> T_BC ;
    public:

      BC_Twall_compute() {
        name_store("BC_options",BC_options) ;
        name_store("T_BC",T_BC) ;
        input("BC_options") ;
        output("T_BC") ;
        constraint("Twall_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOptionUnits("Twall","kelvin",T_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_Twall_compute> register_BC_Twall_compute ;

  // Support for qwall for compatibility with CHEM.
  class BC_qwall_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<real> qwall_BC ;
    public:

      BC_qwall_compute() {
        name_store("BC_options",BC_options) ;
        name_store("qwall_BC",qwall_BC) ;
        input("BC_options") ;
        output("qwall_BC") ;
        constraint("qwall_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOptionUnits("qwall","watt/m/m",qwall_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_qwall_compute> register_BC_qwall_compute ;

  class BC_cartesianT_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<CartesianBoundaryCondition> cartesianT_BC ;
    public:

      BC_cartesianT_compute() {
        name_store("BC_options",BC_options) ;
        name_store("cartesianT_BC",cartesianT_BC) ;
        input("BC_options") ;
        output("cartesianT_BC") ;
        constraint("cartesianT_BCoption") ;
        disable_threading() ;
      }

      void calculate(Entity e) {
        string fileName ; BC_options[e].getOption("cartesianT",fileName) ;
        cartesianT_BC[e].ReadFile(fileName) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_cartesianT_compute> register_BC_cartesianT_compute ;

  class BC_axisymmetricT_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<string> axisymmetricT_BC ;
    public:

      BC_axisymmetricT_compute() {
        name_store("BC_options",BC_options) ;
        name_store("axisymmetricT_BC",axisymmetricT_BC) ;
        input("BC_options") ;
        output("axisymmetricT_BC") ;
        constraint("axisymmetricT_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOption("axisymmetricT",axisymmetricT_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_axisymmetricT_compute> register_BC_axisymmetricT_compute ;

  class BC_T0_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<real> T0_BC ;
    public:

      BC_T0_compute() {
        name_store("BC_options",BC_options) ;
        name_store("T0_BC",T0_BC) ;
        input("BC_options") ;
        output("T0_BC") ;
        constraint("T0_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOptionUnits("T0","kelvin",T0_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_T0_compute> register_BC_T0_compute ;

//-----------------------------------------------------------------------------
// Options for Treservoir, which is the temperature of a reservoir next to a
// a finite thickness wall separating the fluid domain from the reservoir. This
// option can be specified in any of the following ways:
//   1) Treservoir=constant
//   2) Treservoir="f(x,y,z)"

  // Creates the reservoir temperature boundary condition constraints.
  class TreservoirBoundaryConditionConstraints : public constraint_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      Constraint Treservoir_const_BC ;
      Constraint Treservoir_function_BC ;
      Constraint Treservoir_functionSteady_BC,Treservoir_functionUnsteady_BC ;
    private:
      entitySet TreservoirConstant,TreservoirFunction ;
      entitySet TreservoirFunctionSteady,TreservoirFunctionUnsteady ;
    public:

      TreservoirBoundaryConditionConstraints() {
        name_store("BC_options",BC_options) ;
        name_store("Treservoir_const_BC",Treservoir_const_BC) ;
        name_store("Treservoir_function_BC",Treservoir_function_BC) ;
        name_store("Treservoir_functionSteady_BC",Treservoir_functionSteady_BC) ;
        name_store("Treservoir_functionUnsteady_BC",Treservoir_functionUnsteady_BC) ;
        input("BC_options") ;
        output("Treservoir_const_BC,Treservoir_function_BC") ;
        output("Treservoir_functionSteady_BC,Treservoir_functionUnsteady_BC") ;
        constraint("Treservoir_BCoption") ;
      }

      // Note that we need to include UNIT_VALUE in the if check below
      // since including units on a scalar value for TReservoir causes the
      // option value type to be set to UNIT_VALUE instead of REAL.
      void calculate(Entity e) {
        Loci::option_value_type optionValueType=BC_options[e].
          getOptionValueType("Treservoir") ;
        if(optionValueType==Loci::REAL || optionValueType==Loci::UNIT_VALUE){
          TreservoirConstant+=e ;
        }else if(optionValueType==Loci::STRING){
          string function ;
          BC_options[e].getOption("Treservoir",function) ;
          bool isUnsteady=false ;
          Loci::exprP p=Loci::expression::create(function) ;
          set<string> namelist ; Loci::getVarNames(p,namelist) ;
          if(namelist.find(string("t")) != namelist.end()) isUnsteady=true ;
          if(isUnsteady) TreservoirFunctionUnsteady+=e ;
          else TreservoirFunctionSteady+=e ;
          TreservoirFunction+=e ;
        }else{
          cerr << "ERROR: Bad type for Treservoir boundary condition!" << endl ;
          Loci::Abort() ;
        }
      }

      void compute(const sequence &seq) {
        do_loop(seq,this) ; Treservoir_const_BC=TreservoirConstant ;
        Treservoir_function_BC=TreservoirFunction ;
        Treservoir_functionSteady_BC=TreservoirFunctionSteady ;
        Treservoir_functionUnsteady_BC=TreservoirFunctionUnsteady ;
      }
  } ;

  register_rule<TreservoirBoundaryConditionConstraints>
    registerTreservoirBoundaryConditionConstraints ;

  class ConstantTreservoirBCCompute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<real> constantTreservoir_BC ;
    public:

      ConstantTreservoirBCCompute() {
        name_store("BC_options",BC_options) ;
        name_store("constantTreservoir_BC",constantTreservoir_BC) ;
        input("BC_options") ;
        output("constantTreservoir_BC") ;
        constraint("Treservoir_const_BC") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOptionUnits("Treservoir","K",constantTreservoir_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ConstantTreservoirBCCompute> registerConstantTreservoirBCCompute ;

//-----------------------------------------------------------------------------
// Options for Rwall which is the thermal resistance associated with a finite
// thickness wall. This option can be specified in any of the following ways:
//   1) Rwall=constant
//   2) Rwall="f(T)", where T is temperature.

  // Creates the reservoir temperature boundary condition constraints.
  class RwallBoundaryConditionConstraints : public constraint_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      Constraint Rwall_const_BC ;
      Constraint Rwall_function_BC ;
    private:
      entitySet RwallConstant,RwallFunction ;
    public:

      RwallBoundaryConditionConstraints() {
        name_store("BC_options",BC_options) ;
        name_store("Rwall_const_BC",Rwall_const_BC) ;
        name_store("Rwall_function_BC",Rwall_function_BC) ;
        input("BC_options") ;
        output("Rwall_const_BC,Rwall_function_BC") ;
        constraint("Rwall_BCoption") ;
      }

      // Note that we need to include UNIT_VALUE in the if check below
      // since including units on a scalar value for TReservoir causes the
      // option value type to be set to UNIT_VALUE instead of REAL.
      void calculate(Entity e) {
        Loci::option_value_type optionValueType=BC_options[e].
          getOptionValueType("Rwall") ;
        if(optionValueType==Loci::REAL || optionValueType==Loci::UNIT_VALUE){
          RwallConstant+=e ;
        }else if(optionValueType==Loci::STRING){
          RwallFunction+=e ;
        }else{
          cerr << "ERROR: Bad type for Rwall boundary condition!" << endl ;
          Loci::Abort() ;
        }
      }

      void compute(const sequence &seq) {
        do_loop(seq,this) ; Rwall_const_BC=RwallConstant ;
        Rwall_function_BC=RwallFunction ;
      }
  } ;

  register_rule<RwallBoundaryConditionConstraints>
    registerRwallBoundaryConditionConstraints ;

  class ConstantRwallBCCompute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<real> constantRwall_BC ;
    public:

      ConstantRwallBCCompute() {
        name_store("BC_options",BC_options) ;
        name_store("constantRwall_BC",constantRwall_BC) ;
        input("BC_options") ;
        output("constantRwall_BC") ;
        constraint("Rwall_const_BC") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOptionUnits("Rwall","m*m*K/watt",constantRwall_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ConstantRwallBCCompute> registerConstantRwallBCCompute ;

//-----------------------------------------------------------------------------
// Options for turbulence quantities.

  class BC_k_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<real> k_BC ;
    public:

      BC_k_compute() {
        name_store("BC_options",BC_options) ;
        name_store("k_BC",k_BC) ;
        input("BC_options") ;
        output("k_BC") ;
        constraint("k_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOption("k",k_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_k_compute> register_BC_k_compute ;

  class BC_cartesianK_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<string> cartesianK_BC ;
    public:

      BC_cartesianK_compute() {
        name_store("BC_options",BC_options) ;
        name_store("cartesianK_BC",cartesianK_BC) ;
        input("BC_options") ;
        output("cartesianK_BC") ;
        constraint("cartesianK_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOption("cartesianK",cartesianK_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_cartesianK_compute> register_BC_cartesianK_compute ;

  class BC_axisymmetricK_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<string> axisymmetricK_BC ;
    public:

      BC_axisymmetricK_compute() {
        name_store("BC_options",BC_options) ;
        name_store("axisymmetricK_BC",axisymmetricK_BC) ;
        input("BC_options") ;
        output("axisymmetricK_BC") ;
        constraint("axisymmetricK_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOption("axisymmetricK",axisymmetricK_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_axisymmetricK_compute> register_BC_axisymmetricK_compute ;

  class BC_omega_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<real> omega_BC ;
    public:

      BC_omega_compute() {
        name_store("BC_options",BC_options) ;
        name_store("omega_BC",omega_BC) ;
        input("BC_options") ;
        output("omega_BC") ;
        constraint("omega_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOption("omega",omega_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_omega_compute> register_BC_omega_compute ;

  class BC_cartesianOmega_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<string> cartesianOmega_BC ;
    public:

      BC_cartesianOmega_compute() {
        name_store("BC_options",BC_options) ;
        name_store("cartesianOmega_BC",cartesianOmega_BC) ;
        input("BC_options") ;
        output("cartesianOmega_BC") ;
        constraint("cartesianOmega_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOption("cartesianOmega",cartesianOmega_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_cartesianOmega_compute> register_BC_cartesianOmega_compute ;

  class BC_axisymmetricOmega_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<string> axisymmetricOmega_BC ;
    public:

      BC_axisymmetricOmega_compute() {
        name_store("BC_options",BC_options) ;
        name_store("axisymmetricOmega_BC",axisymmetricOmega_BC) ;
        input("BC_options") ;
        output("axisymmetricOmega_BC") ;
        constraint("axisymmetricOmega_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOption("axisymmetricOmega",axisymmetricOmega_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_axisymmetricOmega_compute>
    register_BC_axisymmetricOmega_compute ;

  // In the original StreamUns, the massFlux boundary condition is mass per
  // unit time. In the LOCI version, it is mass per unit area per unit time.
  // This allows us to get around having to store the total boundary area for
  // each face associated with the fixed-mass boundary. Bug fix on 02/16/2005.
  // Added negative to give proper mass flux value for incoming flow.
  class BC_massFlux_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<real> massFlux_BC ;
    public:

      BC_massFlux_compute() {
        name_store("BC_options",BC_options) ;
        name_store("massFlux_BC",massFlux_BC) ;
        input("BC_options") ;
        output("massFlux_BC") ;
        constraint("massFlux_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOptionUnits("massFlux","kg/s/m/m",massFlux_BC[e]) ;
        massFlux_BC[e]*=-1.0 ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_massFlux_compute> register_BC_massFlux_compute ;

  // Temporary hack here since vect3d options don't seem to be supported.
  // Here we just read it as a real and assume (1.0,0.0,0.0).
  class FlowDirectionBoundaryCondition : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<vect3d> flowDirection_BC ;
    public:
      FlowDirectionBoundaryCondition() {
        name_store("BC_options",BC_options) ;
        name_store("flowDirection_BC",flowDirection_BC) ;
        input("BC_options") ;
        output("flowDirection_BC") ;
        constraint("flowDirection_BCoption") ;
      }

      // Get the flow direction. Normalize it for the user.
      void calculate(Entity e) {
        get_vect3dOption(BC_options[e],"flowDirection","",flowDirection_BC[e]) ;
        real flowDirectionNorm=norm(flowDirection_BC[e]) ;
        if(flowDirectionNorm!=0.0){
          flowDirection_BC[e]/=norm(flowDirection_BC[e]) ;
        }else{
          cerr << "ERROR: Flow direction has zero magnitude!" ; exit(1) ;
        }
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FlowDirectionBoundaryCondition>
    registerFlowDirectionBoundaryCondition ;

  class BC_mdot_compute : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<real> mdot_BC ;
    public:

      BC_mdot_compute() {
        name_store("BC_options",BC_options) ;
        name_store("mdot_BC",mdot_BC) ;
        input("BC_options") ;
        output("mdot_BC") ;
        constraint("mdot_BCoption") ;
      }

      void calculate(Entity e) {
        BC_options[e].getOptionUnits("mdot","kg/s",mdot_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BC_mdot_compute> register_BC_mdot_compute ;

  // Rule to initialize the boundary area.
  class BoundaryAreaUnit: public unit_rule {
    private:
      store<real> boundary_area ;
    public:

      // Define input and output.
      BoundaryAreaUnit() {
        name_store("boundary_area",boundary_area) ;
        output("boundary_area") ;
        constraint("UNIVERSE") ;
      }

      // Initialize for a single entity.
      void calculate(Entity e) { boundary_area[e]=0.0 ; }

      // Initialize for a sequence of entities.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryAreaUnit> registerBoundaryAreaUnit ;

  // Rule to finish computation of the boundary area.
  class BoundaryAreaApply : public apply_rule<store<real>,Loci::
  Summation<real> > {
    private:
      store<real> boundary_area ;
      const_Map ref ;
      const_store<Area> area ;
    public:

      // Define input and output.
      BoundaryAreaApply() {
        name_store("boundary_area",boundary_area) ;
        name_store("ref",ref) ;
        name_store("area",area) ;
        input("area") ;
        input("ref->boundary_area") ;
        output("ref->boundary_area") ;
        constraint("ref,ci") ;
      }

      // Compute for a single entity.
      void calculate(Entity face) {
        join(boundary_area[ref[face]],area[face].sada) ;
      }

      // Compute for a sequence of entities.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryAreaApply> registerBoundaryAreaApply ;

  // Rule to compute the massFlux boundary condition given the mass flow rate
  // boundary condition. Bug fix on 02/16/2005. Added negative to give proper
  // mass flux value for incoming flow.
  class MassFluxFromMassFlowRate: public pointwise_rule {
    private:
      const_store<real> mdot,boundary_area ;
      store<real> massFlux ;
    public:

      // Define input and output.
      MassFluxFromMassFlowRate() {
        name_store("mdot_BC",mdot) ;
        name_store("massFlux_BC",massFlux) ;
        name_store("boundary_area",boundary_area) ;
        input("mdot_BC,boundary_area") ;
        output("massFlux_BC") ;
      }

      // Calculate the mass flux for a single entity.
      void calculate(Entity e) { massFlux[e]=-mdot[e]/boundary_area[e] ; }

      // Calculate the mass flux for a sequence of entities.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<MassFluxFromMassFlowRate> registerMassFluxFromMassFlowRate ;

  class BC_y_compute : public pointwise_rule {
    private:
      const_param<EOS> eos ;
      const_store<Loci::options_list> BC_options ;
      storeVec<real> y_BC ;
    public:

      // Define input and output.
      BC_y_compute() {
        name_store("eos",eos) ;
        name_store("BC_options",BC_options) ;
        name_store("y_BC",y_BC) ;
        input("eos,BC_options") ;
        output("y_BC") ;
        constraint("mixture_BCoption") ;
      }

      // Assign species mass fractions for a single entity.
      void calculate(Entity e) {
        const Loci::options_list &optionsList=BC_options[e] ;
        if(optionsList.optionExists("mixture")){
          if(optionsList.getOptionValueType("mixture")!=Loci::LIST){
            cerr << "Mixture must be specified as a species list." << endl ;
            Loci::Abort() ;
          }else{
            Loci::options_list::arg_list speciesList ;
            optionsList.getOption("mixture",speciesList) ;
            Loci::options_list::arg_list::iterator yPtr ;
            for(int i=0;i<eos->numSpecies();++i) y_BC[e][i]=0.0 ;
            for(Loci::options_list::arg_list::iterator yPtr=speciesList.begin();
            yPtr!=speciesList.end();++yPtr) {
              Loci::option_values::value_list_type speciesArg ;
              yPtr->get_value(speciesArg) ;
              if(yPtr->type_of()!=Loci::NAME_ASSIGN || speciesArg.size()!=1 ||
              speciesArg.front().type_of()!=Loci::REAL){
                cerr << "Error in mixture assignment." << endl ; Loci::Abort() ;
              }else{
                string speciesName ; int speciesIndex ; double speciesValue ;
                yPtr->get_value(speciesName) ;
                speciesArg.front().get_value(speciesValue) ;
                speciesIndex=eos->speciesIndex(speciesName) ;
                if(speciesIndex==-1){
                  cerr << "ERROR: Species " << speciesName << " does not exist."
                    << endl ; Loci::Abort() ;
                }
                y_BC[e][speciesIndex]=speciesValue ;
              }
            }
          }
        }
      }

      // Assign species mass fractions for a sequence of entities.
      void compute(const sequence &seq) {
        y_BC.setVecSize(eos->numSpecies()) ; do_loop(seq,this) ;
      }
  } ;

  register_rule<BC_y_compute> register_BC_y_compute ;

  class ReferenceFrameBC : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<unsigned int> referenceFrame_BC ;
    public:

      ReferenceFrameBC() {
        name_store("BC_options",BC_options) ;
        name_store("referenceFrame_BC",referenceFrame_BC) ;
        input("BC_options") ;
        output("referenceFrame_BC") ;
        constraint("referenceFrame_BCoption") ;
      }

      void calculate(Entity e) {
        real temp ; BC_options[e].getOption("referenceFrame",temp) ;
        if(temp<0.0){
          cerr << "ERROR: Bad reference frame number in .vars file." << endl ;
          Loci::Abort() ;
        }
        referenceFrame_BC[e]=(unsigned int)(temp+0.01) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ReferenceFrameBC> registerReferenceFrameBC ;

  class MomentCenterBC : public pointwise_rule {
    private:
      const_store<Loci::options_list> BC_options ;
      store<vect3d> momentCenter_BC ;
    public:

      MomentCenterBC() {
        name_store("BC_options",BC_options) ;
        name_store("momentCenter_BC",momentCenter_BC) ;
        input("BC_options") ;
        output("momentCenter_BC") ;
        constraint("momentCenter_BCoption") ;
      }

      void calculate(Entity e) {
        get_vect3dOption(BC_options[e],"momentCenter","m",momentCenter_BC[e]) ;
      }

      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<MomentCenterBC> registerMomentCenterBC ;
}
