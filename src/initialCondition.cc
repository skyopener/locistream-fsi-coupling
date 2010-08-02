// System includes.
#include <rpc/rpc.h>
#include <rpc/xdr.h>

// Standard library includes.
#include <vector>
using std::vector ;

// Loci includes.
#include <fact_db.h>
#include <Tools/parse.h>
using namespace Loci::parse ;

// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;

// StreamUns includes.
#include "initialCondition.h"
#include "name_var.h"
#include "sciTypes.h"

namespace streamUns {

  //---------------------------------------------------------------------------
  // Class InitialCondition.

  // Returns the serialized buffer size.
  int InitialCondition::BufferSize() const {
    int bufferSize=16+2*speciesMassFraction.size() ;
    for(map<string,real>::const_iterator m=speciesMassFraction.begin();m!=
      speciesMassFraction.end();++m) bufferSize+=m->first.size() ;
    return bufferSize ;
  }

  // Overridden virtual method called by fact database when a fact of type
  // InitialCondition is read.
  istream& InitialCondition::Input(istream &in) {
    InitialConditionOptions optionsList ;
    in >> optionsList ; Input(optionsList) ; return in ;
  }

  // Gets the values for the initial conditions from the specified options.
  void InitialCondition::Input(const options_list &optionsList) {

    // Density.
    if(optionsList.optionExists("rho")){
      if(optionsList.getOptionValueType("rho")==Loci::REAL) {
        optionsList.getOption("rho",density) ; densityDefined=true ;
      }else if(optionsList.getOptionValueType("rho")==Loci::UNIT_VALUE){
        Loci::UNIT_type densityUnitType ;
        optionsList.getOption("rho",densityUnitType) ;
        if(!densityUnitType.is_compatible("kg/m/m/m")){
          cerr << "Wrong type of unit for density: " << densityUnitType
            << endl ;
        }else{
          density=densityUnitType.get_value_in("kg/m/m/m") ;
          densityDefined=true ;
        }
      }else{
        cerr << "Incorrect type for 'rho' in InitialCondition." << endl ;
      }
    }

    // Velocity.
    if(optionsList.optionExists("v")){
      velocity=get_vect3d(optionsList,"v","m/s") ; velocityDefined=true ;
    }

    // Temperature.
    if(optionsList.optionExists("T")){
      if(optionsList.getOptionValueType("T")==Loci::REAL) {
        optionsList.getOption("T",temperature) ; temperatureDefined=true ;
      }else if(optionsList.getOptionValueType("T")==Loci::UNIT_VALUE){
        Loci::UNIT_type temperatureUnitType ;
        optionsList.getOption("T",temperatureUnitType) ;
        if(!temperatureUnitType.is_compatible("kelvin")){
          cerr << "Wrong type of unit for temperature: " << temperatureUnitType 
            << endl ;
        }else{
          temperature=temperatureUnitType.get_value_in("kelvin") ;
          temperatureDefined=true ;
        }
      }else{
        cerr << "Incorrect type for 'T' in InitialCondition." << endl ;
      }
    }

    // Pressure.
    if(optionsList.optionExists("p")){
      if(optionsList.getOptionValueType("p")==Loci::REAL) {
        optionsList.getOption("p",pressure) ; pressureDefined=true ;
      }else if(optionsList.getOptionValueType("p")==Loci::UNIT_VALUE){
        Loci::UNIT_type pressureUnitType ;
        optionsList.getOption("p",pressureUnitType) ;
        if(!pressureUnitType.is_compatible("Pa")){
          cerr << "Wrong type of unit for pressure: " << pressureUnitType 
            << endl ;
        }else{
          pressure=pressureUnitType.get_value_in("Pa") ;
          pressureDefined=true ;
        }
      }else{
        cerr << "Incorrect type for 'p' in InitialCondition." << endl ;
      }
    }

    // K.
    if(optionsList.optionExists("k")){
      if(optionsList.getOptionValueType("k")==Loci::REAL) {
        optionsList.getOption("k",k) ; kDefined=true ;
      }else if(optionsList.getOptionValueType("k")==Loci::UNIT_VALUE){
        Loci::UNIT_type kUnitType ;
        optionsList.getOption("k",kUnitType) ;
        if(!kUnitType.is_compatible("m*m/s/s")){
          cerr << "Wrong type of unit for k: " << kUnitType << endl ;
        }else{
          k=kUnitType.get_value_in("m*m/s/s") ; kDefined=true ;
        }
      }else{
        cerr << "Incorrect type for 'k' in InitialCondition." << endl ;
      }
    }

    // Omega.
    if(optionsList.optionExists("omega")){
      if(optionsList.getOptionValueType("omega")==Loci::REAL) {
        optionsList.getOption("omega",omega) ; omegaDefined=true ;
      }else if(optionsList.getOptionValueType("omega")==Loci::UNIT_VALUE){
        Loci::UNIT_type omegaUnitType ;
        optionsList.getOption("omega",omegaUnitType) ;
        if(!omegaUnitType.is_compatible("s^-1")){
          cerr << "Wrong type of unit for omega: " << omegaUnitType << endl ;
        }else{
          omega=omegaUnitType.get_value_in("s^-1") ; omegaDefined=true ;
        }
      }else{
        cerr << "Incorrect type for 'omega' in InitialCondition." << endl ;
      }
    }

    // Species mass fractions.
    if(optionsList.optionExists("mixture")){
      if(optionsList.getOptionValueType("mixture")!=Loci::LIST){
        cerr << "Mixture must be specified as a species list." << endl ;
        Loci::Abort() ;
      }else{
        speciesMassFractionDefined=true ;
        Loci::options_list::arg_list speciesList ;
        optionsList.getOption("mixture",speciesList) ;
        Loci::options_list::arg_list::iterator yPtr ;
        for(Loci::options_list::arg_list::iterator yPtr=speciesList.begin();
        yPtr!=speciesList.end();++yPtr) {
          Loci::option_values::value_list_type speciesArg ;
          yPtr->get_value(speciesArg) ;
          if(yPtr->type_of()!=Loci::NAME_ASSIGN || speciesArg.size()!=1 ||
          speciesArg.front().type_of()!=Loci::REAL){
            cerr << "Error in mixture assignment." << endl ; Loci::Abort() ;
          }else{
            string speciesName ; double speciesValue ;
            yPtr->get_value(speciesName) ;
            speciesArg.front().get_value(speciesValue) ;
            speciesMassFraction[speciesName]=speciesValue ;
          }
        }
      }
    }
  }

  // Packs the data into a buffer.
  void InitialCondition::PackBuffer(real *buffer,int size) {
    int i=0 ;
    buffer[i++]=density ; buffer[i++]=temperature ; buffer[i++]=pressure ;
    buffer[i++]=k ; buffer[i++]=omega ;
    buffer[i++]=velocity.x ; buffer[i++]=velocity.y ; buffer[i++]=velocity.z ;
    buffer[i++]=speciesMassFraction.size() ;
    for(map<string,real>::const_iterator m=speciesMassFraction.begin();m!=
    speciesMassFraction.end();++m){
      string name=m->first ; int size=name.size() ; buffer[i++]=real(size) ;
      for(int j=0;j<size;++j) buffer[i++]=real(name[j]) ;
      buffer[i++]=m->second ;
    }
    buffer[i++]=real(densityDefined) ; buffer[i++]=real(temperatureDefined) ;
    buffer[i++]=real(pressureDefined) ; buffer[i++]=real(kDefined) ;
    buffer[i++]=real(omegaDefined) ; buffer[i++]=real(velocityDefined) ;
    buffer[i++]=real(speciesMassFractionDefined) ;
  }

  // Prints the initial condition.
  ostream& InitialCondition::Print(ostream &out) const {
    int outputCount=0 ;
    out << "<" ;
    if(densityDefined){ out << "rho=" << density ; ++outputCount ; }
    if(velocityDefined){
      if(outputCount!=0) out << "," ;
      out << "v=[" << velocity.x << "," << velocity.y << "," << velocity.z
        << "]" ;
      ++outputCount ;
    }
    if(pressureDefined){
      if(outputCount!=0) out << "," ;
      out << "p=" << pressure ; ++outputCount ;
    }
    if(temperatureDefined){
      if(outputCount!=0) out << "," ;
      out << "T=" << temperature ; ++outputCount ;
    }
    if(kDefined){
      if(outputCount!=0) out << "," ;
      out << "k=" << k ; ++outputCount ;
    }
    if(omegaDefined){
      if(outputCount!=0) out << "," ;
      out << "omega=" << omega ; ++outputCount ;
    }
    out << ">" ; return out ;
  }

  // Unpacks the data from a buffer.
  void InitialCondition::UnpackBuffer(real *buffer,int size) {
    int i=0 ;
    density=buffer[i++] ; temperature=buffer[i++] ; pressure=buffer[i++] ;
    k=buffer[i++] ; omega=buffer[i++] ;
    velocity.x=buffer[i++] ; velocity.y=buffer[i++] ; velocity.z=buffer[i++] ;
    int numSpecies=int(buffer[i++]) ; speciesMassFraction.clear() ;
    for(int j=0;j<numSpecies;++j){
      int nameSize=int(buffer[i++]) ; string name ;
      for(int k=0;k<nameSize;++k){ char c=char(buffer[i++]) ; name+=c ; }
      speciesMassFraction[name]=buffer[i++] ;
    }
    densityDefined=bool(buffer[i++]) ; temperatureDefined=bool(buffer[i++]) ;
    pressureDefined=bool(buffer[i++]) ; kDefined=bool(buffer[i++]) ;
    omegaDefined=bool(buffer[i++]) ; velocityDefined=bool(buffer[i++]) ;
    speciesMassFractionDefined=bool(buffer[i++]) ;
  }

  //---------------------------------------------------------------------------
  // Rules.

  // Rule to assign the initial condition for density for incompressible flows.
  class DensityInitialConditionIncompressible : public pointwise_rule {
    private:
      const_param<InitialCondition> initialCondition ;
      store<real> rho_ic ;
    public:

      // Define input and output.
      DensityInitialConditionIncompressible() {
        name_store("initialCondition",initialCondition) ;
        name_store("rho_ic",rho_ic) ;
        input("initialCondition") ;
        output("rho_ic") ;
        constraint("geom_cells,incompressibleFlow") ;
      }

      // Assign density for a single cell.
      void calculate(Entity cell) {
        rho_ic[cell]=initialCondition->Density() ;
      }

      // Assign density for a sequence of cells.
      void compute(const sequence &seq) {

        // Check that density was provided by the user.
        if(!initialCondition->IsDensityDefined())
          cerr << "ERROR: Initial condition for density required!" << endl ;

        // Loop through the sequence of cells.
        do_loop(seq,this) ;
      }
  } ;

  register_rule<DensityInitialConditionIncompressible>
    registerDensityInitialConditionIncompressible ;

  // Rule to assign the initial condition for density from ql and qr for
  // incompressible flows.
  class DensityInitialConditionQLQR : public pointwise_rule {
    private:
      const_param<InitialCondition> ql,qr ;
      const_param<real> xMid ;
      const_store<vect3d> cellCenter ;
      store<real> rho_ic ;
    public:

      // Define input and output.
      DensityInitialConditionQLQR() {
        name_store("ql",ql) ;
        name_store("qr",qr) ;
        name_store("xmid",xMid) ;
        name_store("cellcenter",cellCenter) ;
        name_store("rho_ic",rho_ic) ;
        input("ql,qr,xmid,cellcenter") ;
        output("rho_ic") ;
        constraint("geom_cells,incompressibleFlow") ;
      }

      // Assign density for a single cell.
      void calculate(Entity cell) {
        if(cellCenter[cell].x > *xMid) rho_ic[cell]=qr->Density() ;
        else rho_ic[cell]=ql->Density() ;
      }

      // Assign density for a sequence of cells.
      void compute(const sequence &seq) {

        // Check that density was provided by the user.
        if(!ql->IsDensityDefined())
          cerr << "ERROR: Missing density in ql !" << endl ;
        if(!qr->IsDensityDefined())
          cerr << "ERROR: Missing density in qr !" << endl ;

        // Loop through the sequence of cells.
        do_loop(seq,this) ;
      }
  } ;

  register_rule<DensityInitialConditionQLQR>
    registerDensityInitialConditionQLQR ;

  // Rule to assign the initial condition for density for compressible flows.
  class DensityInitialConditionCompressible : public pointwise_rule {
    private:
      const_store<EOS::State> eos_state_ic ;
      store<real> rho_ic ;
    public:
                                                                                                    
      // Define input and output.
      DensityInitialConditionCompressible() {
        name_store("eos_state_ic",eos_state_ic) ;
        name_store("rho_ic",rho_ic) ;
        input("eos_state_ic") ;
        output("rho_ic") ;
        constraint("geom_cells,compressibleFlow") ;
      }
                                                                                                    
      // Assign density for a single cell.
      void calculate(Entity cell) {
        rho_ic[cell]=eos_state_ic[cell].density() ;
      }
                                                                                                    
      // Assign density for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DensityInitialConditionCompressible>
    registerDensityInitialConditionCompressible ;

  // Rule to assign the initial condition for velocity.
  class VelocityInitialCondition : public pointwise_rule {
    private:
      const_param<InitialCondition> initialCondition ;
      store<vect3d> v_ic ;
    public:

      // Define input and output.
      VelocityInitialCondition() {
        name_store("initialCondition",initialCondition) ;
        name_store("v_ic",v_ic) ;
        input("initialCondition") ;
        output("v_ic") ;
        constraint("geom_cells") ;
      }

      // Assign velocity for a single cell.
      void calculate(Entity cell) {
        v_ic[cell]=initialCondition->Velocity() ;
      }

      // Assign velocity for a sequence of cells.
      void compute(const sequence &seq) {

        // Check that velocity was provided by the user.
        if(!initialCondition->IsVelocityDefined())
          cerr << "ERROR: Initial condition for velocity required!" << endl ;

        // Loop through the sequence of cells.
        do_loop(seq,this) ;
      }
  } ;

  register_rule<VelocityInitialCondition> registerVelocityInitialCondition ;

  // Rule to assign the initial condition for velocity from ql and qr.
  class VelocityInitialConditionQLQR : public pointwise_rule {
    private:
      const_param<InitialCondition> ql,qr ;
      const_param<real> xMid ;
      const_store<vect3d> cellCenter ;
      store<vect3d> v_ic ;
    public:

      // Define input and output.
      VelocityInitialConditionQLQR() {
        name_store("ql",ql) ;
        name_store("qr",qr) ;
        name_store("xmid",xMid) ;
        name_store("cellcenter",cellCenter) ;
        name_store("v_ic",v_ic) ;
        input("ql,qr,xmid,cellcenter") ;
        output("v_ic") ;
        constraint("geom_cells") ;
      }

      // Assign velocity for a single cell.
      void calculate(Entity cell) {
        if(cellCenter[cell].x > *xMid) v_ic[cell]=qr->Velocity() ;
        else v_ic[cell]=ql->Velocity() ;
      }

      // Assign velocity for a sequence of cells.
      void compute(const sequence &seq) {

        // Check that velocity was provided by the user.
        if(!ql->IsVelocityDefined())
          cerr << "ERROR: Missing velocity in ql !" << endl ;
        if(!qr->IsVelocityDefined())
          cerr << "ERROR: Missing velocity in qr !" << endl ;

        // Loop through the sequence of cells.
        do_loop(seq,this) ;
      }
  } ;

  register_rule<VelocityInitialConditionQLQR>
    registerVelocityInitialConditionQLQR ;

  // Had to do this, but maybe will get rid of later.
  class BoundaryVelocityInitialCondition : public pointwise_rule {
    private:
      const_Map ref ;
      const_param<real> timeStep ;
      const_store<vect3d> faceCenter ;
      const_store<string> functionVX_BC,functionVY_BC,functionVZ_BC ;
      store<vect3d> v_f_ic ;
    public:

      // Define input and output.
      BoundaryVelocityInitialCondition() {
        name_store("ref",ref) ;
        name_store("timeStep",timeStep) ;
        name_store("facecenter",faceCenter) ;
        name_store("functionVX_BC",functionVX_BC) ;
        name_store("functionVY_BC",functionVY_BC) ;
        name_store("functionVZ_BC",functionVZ_BC) ;
        name_store("v_f_ic",v_f_ic) ;
        input("timeStep,facecenter") ;
        input("ref->(functionVX_BC,functionVY_BC,functionVZ_BC)") ;
        output("v_f_ic") ;
        constraint("ref->v_functionUnsteady_BC") ;
      }

      // Calculate face velocity for a single face.
      void calculate(Entity face) {
        map<string,real> varMap ; varMap["x"]=faceCenter[face].x ;
        varMap["y"]=faceCenter[face].y ; varMap["z"]=faceCenter[face].z ;
        varMap["t"]=(*timeStep) ;
        Loci::exprP pX=Loci::expression::create(functionVX_BC[ref[face]]) ;
        Loci::exprP pY=Loci::expression::create(functionVY_BC[ref[face]]) ;
        Loci::exprP pZ=Loci::expression::create(functionVZ_BC[ref[face]]) ;
        real vX=pX->evaluate(varMap),vY=pY->evaluate(varMap),
          vZ=pZ->evaluate(varMap) ;
        v_f_ic[face]=vect3d(vX,vY,vZ) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryVelocityInitialCondition>
    registerBoundaryVelocityInitialCondition ;

  // Rule to assign the initial condition for pressure.
  class PressureInitialCondition : public pointwise_rule {
    private:
      const_param<InitialCondition> initialCondition ;
      store<real> p_ic ;
    public:

      // Define input and output.
      PressureInitialCondition() {
        name_store("initialCondition",initialCondition) ;
        name_store("p_ic",p_ic) ;
        input("initialCondition") ;
        output("p_ic") ;
        constraint("geom_cells") ;
      }

      // Assign pressure for a single cell.
      void calculate(Entity cell) {
        p_ic[cell]=initialCondition->Pressure() ;
      }

      // Assign pressure for a sequence of cells.
      void compute(const sequence &seq) {

        // Check that pressure was provided by the user.
        if(!initialCondition->IsPressureDefined())
          cerr << "ERROR: Initial condition for pressure required!" << endl ;

        // Loop through the sequence of cells.
        do_loop(seq,this) ;
      }
  } ;

  register_rule<PressureInitialCondition> registerPressureInitialCondition ;

  // Rule to assign the initial condition for pressure from ql and qr.
  class PressureInitialConditionQLQR : public pointwise_rule {
    private:
      const_param<InitialCondition> ql,qr ;
      const_param<real> xMid ;
      const_store<vect3d> cellCenter ;
      store<real> p_ic ;
    public:

      // Define input and output.
      PressureInitialConditionQLQR() {
        name_store("ql",ql) ;
        name_store("qr",qr) ;
        name_store("xmid",xMid) ;
        name_store("cellcenter",cellCenter) ;
        name_store("p_ic",p_ic) ;
        input("ql,qr,xmid,cellcenter") ;
        output("p_ic") ;
        constraint("geom_cells") ;
      }

      // Assign pressure for a single cell.
      void calculate(Entity cell) {
        if(cellCenter[cell].x > *xMid) p_ic[cell]=qr->Pressure() ;
        else p_ic[cell]=ql->Pressure() ;
      }

      // Assign pressure for a sequence of cells.
      void compute(const sequence &seq) {

        // Check that pressure was provided by the user.
        if(!ql->IsPressureDefined())
          cerr << "ERROR: Missing pressure in ql !" << endl ;
        if(!qr->IsPressureDefined())
          cerr << "ERROR: Missing pressure in qr !" << endl ;

        // Loop through the sequence of cells.
        do_loop(seq,this) ;
      }
  } ;

  register_rule<PressureInitialConditionQLQR>
    registerPressureInitialConditionQLQR ;

  // Rule to assign the initial condition for temperature.
  class TemperatureInitialCondition : public pointwise_rule {
    private:
      const_param<InitialCondition> initialCondition ;
      store<real> T_ic ;
    public:

      // Define input and output.
      TemperatureInitialCondition() {
        name_store("initialCondition",initialCondition) ;
        name_store("T_ic",T_ic) ;
        input("initialCondition") ;
        output("T_ic") ;
        constraint("compressibleFlow,geom_cells") ;
      }

      // Assign temperature for a single cell.
      void calculate(Entity cell) {
        T_ic[cell]=initialCondition->Temperature() ;
      }

      // Assign temperature for a sequence of cells.
      void compute(const sequence &seq) {

        // Check that temperature was provided by the user.
        if(!initialCondition->IsTemperatureDefined())
          cerr << "ERROR: Initial condition for temperature required!" << endl ;

        // Loop through the sequence of cells.
        do_loop(seq,this) ;
      }
  } ;

  register_rule<TemperatureInitialCondition>
    registerTemperatureInitialCondition ;

  // Rule to assign the initial condition for temperature from ql and qr.
  class TemperatureInitialConditionQLQR : public pointwise_rule {
    private:
      const_param<InitialCondition> ql,qr ;
      const_param<real> xMid ;
      const_store<vect3d> cellCenter ;
      store<real> T_ic ;
    public:

      // Define input and output.
      TemperatureInitialConditionQLQR() {
        name_store("ql",ql) ;
        name_store("qr",qr) ;
        name_store("xmid",xMid) ;
        name_store("cellcenter",cellCenter) ;
        name_store("T_ic",T_ic) ;
        input("ql,qr,xmid,cellcenter") ;
        output("T_ic") ;
        constraint("compressibleFlow,geom_cells") ;
      }

      // Assign temperature for a single cell.
      void calculate(Entity cell) {
        if(cellCenter[cell].x > *xMid) T_ic[cell]=qr->Temperature() ;
        else T_ic[cell]=ql->Temperature() ;
      }

      // Assign temperature for a sequence of cells.
      void compute(const sequence &seq) {

        // Check that temperature was provided by the user.
        if(!ql->IsTemperatureDefined())
          cerr << "ERROR: Missing temperature in ql !" << endl ;
        if(!qr->IsTemperatureDefined())
          cerr << "ERROR: Missing temperature in qr !" << endl ;

        // Loop through the sequence of cells.
        do_loop(seq,this) ;
      }
  } ;

  register_rule<TemperatureInitialConditionQLQR>
    registerTemperatureInitialConditionQLQR ;

  // Rule to assign the initial condition for total enthalpy.
  class TotalEnthalpyInitialCondition : public pointwise_rule {
    private:
      const_store<EOS::State> eos_state_ic ;
      const_store<vect3d> v_ic ;
      store<real> h_ic ;
    public:

      // Define input and output.
      TotalEnthalpyInitialCondition() {
        name_store("eos_state_ic",eos_state_ic) ;
        name_store("v_ic",v_ic) ;
        name_store("h_ic",h_ic) ;
        input("eos_state_ic,v_ic") ;
        output("h_ic") ;
        constraint("compressibleFlow,geom_cells") ;
      }

      // Assign total enthalpy for a single cell.
      void calculate(Entity cell) {
        h_ic[cell]=eos_state_ic[cell].enthalpy()+0.5*dot(v_ic[cell],
          v_ic[cell]) ;
      }

      // Assign total enthalpy for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TotalEnthalpyInitialCondition>
    registerTotalEnthalpyInitialCondition ;

  // Rule to assign the initial condition for internal energy.
  class InternalEnergyInitialCondition : public pointwise_rule {
    private:
      const_store<EOS::State> eos_state_ic ;
      store<real> h_ic ;
    public:

      // Define input and output.
      InternalEnergyInitialCondition() {
        name_store("eos_state_ic",eos_state_ic) ;
        name_store("h_ic",h_ic) ;
        input("eos_state_ic") ;
        output("h_ic") ;
        constraint("compressibleFlow,geom_cells,energy5") ;
      }

      // Assign internal energy for a single cell.
      void calculate(Entity cell) { h_ic[cell]=eos_state_ic[cell].energy() ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InternalEnergyInitialCondition>
    registerInternalEnergyInitialCondition ;

  // Rule to assign the initial condition for k. Checked.
  class KInitialCondition : public pointwise_rule {
    private:
      const_param<InitialCondition> initialCondition ;
      store<real> k_ic ;
    public:

      // Define input and output.
      KInitialCondition() {
        name_store("initialCondition",initialCondition) ;
        name_store("k_ic",k_ic) ;
        input("initialCondition") ;
        output("k_ic") ;
        constraint("geom_cells") ;
      }

      // Assign k for a single cell.
      void calculate(Entity cell) { k_ic[cell]=initialCondition->K() ; }

      // Assign k for a sequence of cells.
      void compute(const sequence &seq) {

        // Check that k was provided by the user.
        if(!initialCondition->IsKDefined())
          cerr << "ERROR: Initial condition for k required!" << endl ;

        // Loop through the sequence of cells.
        do_loop(seq,this) ;
      }
  } ;

  register_rule<KInitialCondition> registerKInitialCondition ;

  // Rule to assign the initial condition for k from ql and qr.
  class KInitialConditionQLQR : public pointwise_rule {
    private:
      const_param<InitialCondition> ql,qr ;
      const_param<real> xMid ;
      const_store<vect3d> cellCenter ;
      store<real> k_ic ;
    public:

      // Define input and output.
      KInitialConditionQLQR() {
        name_store("ql",ql) ;
        name_store("qr",qr) ;
        name_store("xmid",xMid) ;
        name_store("cellcenter",cellCenter) ;
        name_store("k_ic",k_ic) ;
        input("ql,qr,xmid,cellcenter") ;
        output("k_ic") ;
        constraint("geom_cells") ;
      }

      // Assign k for a single cell.
      void calculate(Entity cell) {
        if(cellCenter[cell].x > *xMid) k_ic[cell]=qr->K() ;
        else k_ic[cell]=ql->K() ;
      }

      // Assign k for a sequence of cells.
      void compute(const sequence &seq) {

        // Check that k was provided by the user.
        if(!ql->IsKDefined())
          cerr << "ERROR: Missing k in ql !" << endl ;
        if(!qr->IsKDefined())
          cerr << "ERROR: Missing k in qr !" << endl ;

        // Loop through the sequence of cells.
        do_loop(seq,this) ;
      }
  } ;

  register_rule<KInitialConditionQLQR> registerKInitialConditionQLQR ;

  // Rule to assign the initial condition for omega. Checked.
  class OmegaInitialCondition : public pointwise_rule {
    private:
      const_param<InitialCondition> initialCondition ;
      store<real> omega_ic ;
    public:

      // Define input and output.
      OmegaInitialCondition() {
        name_store("initialCondition",initialCondition) ;
        name_store("omega_ic",omega_ic) ;
        input("initialCondition") ;
        output("omega_ic") ;
        constraint("geom_cells") ;
      }

      // Assign omega for a single cell.
      void calculate(Entity cell) { omega_ic[cell]=initialCondition->Omega() ; }

      // Assign omega for a sequence of cells.
      void compute(const sequence &seq) {

        // Check that omega was provided by the user.
        if(!initialCondition->IsOmegaDefined())
          cerr << "ERROR: Initial condition for omega required!" << endl ;

        // Loop through the sequence of cells.
        do_loop(seq,this) ;
      }
  } ;

  register_rule<OmegaInitialCondition> registerOmegaInitialCondition ;

  // Rule to assign the initial condition for omega from ql and qr.
  class OmegaInitialConditionQLQR : public pointwise_rule {
    private:
      const_param<InitialCondition> ql,qr ;
      const_param<real> xMid ;
      const_store<vect3d> cellCenter ;
      store<real> omega_ic ;
    public:

      // Define input and output.
      OmegaInitialConditionQLQR() {
        name_store("ql",ql) ;
        name_store("qr",qr) ;
        name_store("xmid",xMid) ;
        name_store("cellcenter",cellCenter) ;
        name_store("omega_ic",omega_ic) ;
        input("ql,qr,xmid,cellcenter") ;
        output("omega_ic") ;
        constraint("geom_cells") ;
      }

      // Assign omega for a single cell.
      void calculate(Entity cell) {
        if(cellCenter[cell].x > *xMid) omega_ic[cell]=qr->Omega() ;
        else omega_ic[cell]=ql->Omega() ;
      }

      // Assign omega for a sequence of cells.
      void compute(const sequence &seq) {

        // Check that omega was provided by the user.
        if(!ql->IsOmegaDefined())
          cerr << "ERROR: Missing omega in ql !" << endl ;
        if(!qr->IsOmegaDefined())
          cerr << "ERROR: Missing omega in qr !" << endl ;

        // Loop through the sequence of cells.
        do_loop(seq,this) ;
      }
  } ;

  register_rule<OmegaInitialConditionQLQR> registerOmegaInitialConditionQLQR ;

  // Rule to assign the initial condition for the species.
  class SpeciesInitialCondition : public pointwise_rule {
    private:
      const_param<EOS> eos ;
      const_param<InitialCondition> initialCondition ;
      storeVec<real> y_ic ;
    private:
      vector<real> speciesValue ;
    public:

      // Define input and output.
      SpeciesInitialCondition() {
        name_store("eos",eos) ;
        name_store("initialCondition",initialCondition) ;
        name_store("y_ic",y_ic) ;
        input("eos,initialCondition") ;
        output("y_ic") ;
        constraint("speciesTransport,geom_cells") ;
      }

      // Assign the species mass fractions for a single cell.
      void calculate(Entity cell) {
        for(int i=0;i<eos->numSpecies();++i) y_ic[cell][i]=speciesValue[i] ;
      }

      // Assign the species mass fractions for a sequence of cells.
      void compute(const sequence &seq) {

        // Check that the species were provided by the user.
        if(!initialCondition->IsSpeciesMassFractionDefined()){
          cerr << "ERROR: Initial condition for species required!" << endl ;
          Loci::Abort() ;
        }

        // Set the species values. The mixture material definition sets the
        // order of the species. Note that any unspecified species mass
        // fractions are initialized to zero by default.
        speciesValue=vector<real>(eos->numSpecies(),0.0) ;
        const map<string,real> &speciesMassFraction=initialCondition->
          SpeciesMassFraction() ;
        for(map<string,real>::const_iterator m=speciesMassFraction.begin();
        m!=speciesMassFraction.end();++m){
          int speciesIndex=eos->speciesIndex(m->first) ;
          if(speciesIndex==-1){
            cerr << "ERROR: Species " << m->first << " does not exist."
              << endl ; Loci::Abort() ;
          }
          speciesValue[speciesIndex]=m->second ;
        }

        // Set the number of species.
        y_ic.setVecSize(eos->numSpecies()) ;

        // Loop through the sequence of cells.
        do_loop(seq,this) ;
      }
  } ;

  register_rule<SpeciesInitialCondition> registerSpeciesInitialCondition ;

  // Rule to assign the initial condition for the species for ql and qr.
  class SpeciesInitialConditionQLQR : public pointwise_rule {
    private:
      const_param<EOS> eos ;
      const_param<InitialCondition> ql,qr ;
      const_param<real> xMid ;
      const_store<vect3d> cellCenter ;
      storeVec<real> y_ic ;
    private:
      vector<real> speciesValueQL,speciesValueQR ;
    public:

      // Define input and output.
      SpeciesInitialConditionQLQR() {
        name_store("eos",eos) ;
        name_store("ql",ql) ;
        name_store("qr",qr) ;
        name_store("xmid",xMid) ;
        name_store("cellcenter",cellCenter) ;
        name_store("y_ic",y_ic) ;
        input("eos,ql,qr,xmid,cellcenter") ;
        output("y_ic") ;
        constraint("speciesTransport,geom_cells") ;
      }

      // Assign the species mass fractions for a single cell.
      void calculate(Entity cell) {
        if(cellCenter[cell].x > *xMid){
          for(int i=0;i<eos->numSpecies();++i) y_ic[cell][i]=speciesValueQR[i] ;
        }else{
          for(int i=0;i<eos->numSpecies();++i) y_ic[cell][i]=speciesValueQL[i] ;
        }
      }

      // Assign the species mass fractions for a sequence of cells.
      void compute(const sequence &seq) {

        // Check that the species were provided by the user.
        if(!ql->IsSpeciesMassFractionDefined()){
          cerr << "ERROR: Missing species mass fractions in ql !" << endl ;
          Loci::Abort() ;
        }
        if(!qr->IsSpeciesMassFractionDefined()){
          cerr << "ERROR: Missing species mass fractions in qr !" << endl ;
          Loci::Abort() ;
        }

        // Set the species values. The mixture material definition sets the
        // order of the species. Note that any unspecified species mass
        // fractions are initialized to zero by default.
        speciesValueQL=vector<real>(eos->numSpecies(),0.0) ;
        speciesValueQR=vector<real>(eos->numSpecies(),0.0) ;
        const map<string,real> &speciesMassFractionQL=ql->
          SpeciesMassFraction() ;
        const map<string,real> &speciesMassFractionQR=qr->
          SpeciesMassFraction() ;
        for(map<string,real>::const_iterator m=speciesMassFractionQL.begin();
        m!=speciesMassFractionQL.end();++m){
          int speciesIndex=eos->speciesIndex(m->first) ;
          if(speciesIndex==-1){
            cerr << "ERROR: Species " << m->first << " does not exist."
              << endl ; Loci::Abort() ;
          }
          speciesValueQL[speciesIndex]=m->second ;
        }
        for(map<string,real>::const_iterator m=speciesMassFractionQR.begin();
        m!=speciesMassFractionQR.end();++m){
          int speciesIndex=eos->speciesIndex(m->first) ;
          if(speciesIndex==-1){
            cerr << "ERROR: Species " << m->first << " does not exist."
              << endl ; Loci::Abort() ;
          }
          speciesValueQR[speciesIndex]=m->second ;
        }

        // Set the number of species.
        y_ic.setVecSize(eos->numSpecies()) ;

        // Loop through the sequence of cells.
        do_loop(seq,this) ;
      }
  } ;

  register_rule<SpeciesInitialConditionQLQR>
    registerSpeciesInitialConditionQLQR ;

  //---------------------------------------------------------------------------
  // EOS initial condition.

  // Rule to assign the initial condition for density for compressible flows.
  class EOSInitialConditionFromYPT : public pointwise_rule {
    private:
      const_param<EOS> eos ;
      const_store<real> p_ic ;
      const_store<real> T_ic ;
      const_storeVec<real> y_ic ;
      store<EOS::State> eos_state_ic ;
      storeVec<real> eos_mixture_state_ic ;
      storeVec<float> hint_n_ic ;
    public:

      // Define input and output.
      EOSInitialConditionFromYPT() {
        name_store("eos",eos) ;
        name_store("p_ic",p_ic) ;
        name_store("T_ic",T_ic) ;
        name_store("y_ic",y_ic) ;
        name_store("eos_state_ic",eos_state_ic) ;
        name_store("eos_mixture_state_ic",eos_mixture_state_ic) ;
        name_store("hint_n_ic",hint_n_ic) ;
        input("eos,p_ic,T_ic,y_ic") ;
        output("eos_state_ic,eos_mixture_state_ic,hint_n_ic") ;
        constraint("geom_cells,compressibleFlow") ;
      }

      // Get the EOS state for a single cell.
      void calculate(Entity cell) {
        eos_state_ic[cell]=eos->State_from_mixture_p_T(y_ic[cell],p_ic[cell],
          T_ic[cell],eos_mixture_state_ic[cell]) ;
        eos[cell].getHint(hint_n_ic[cell],eos_state_ic[cell],
          eos_mixture_state_ic[cell]) ;
      }

      // Get the EOS state for a sequence of cells.
      void compute(const sequence &seq) {
        eos_mixture_state_ic.setVecSize(eos->mixtureStateSize()) ;
        hint_n_ic.setVecSize(eos->hintSize()) ; do_loop(seq,this) ;
      }
  } ;

  register_rule<EOSInitialConditionFromYPT> registerEOSInitialConditionFromYPT ;
}




