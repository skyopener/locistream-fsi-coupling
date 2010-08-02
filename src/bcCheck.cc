// Loci includes.
#include <Loci.h>
                                                                                
// StreamUns includes.
#include <readGrid.h>
#include "sciTypes.h"
                                                                                
namespace streamUns {

  // Checks symmetry and slip boundaries.
  class CheckSlipSymmetry : public BC_Check {
    public:
      string BoundaryConditions() { return "slip,symmetry" ; }
      string VariablesChecked(fact_db &facts) { return "" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        return true ;
      }
      ostream &ErrorMessage(ostream &s) { return s ; }
  } ;

  register_BC<CheckSlipSymmetry> registerCheckSlipSymmetry ;

  // Checks incompressible inlet boundaries.
  class CheckIncompressibleInlet : public BC_Check {
    private:
      string errorMessage ;
    public:
      CheckIncompressibleInlet() : errorMessage("") {}
    public:
      string BoundaryConditions() { return "incompressibleInlet" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        param<string> flowRegime=facts.get_fact("flowRegime") ;
        if(!bc_options.optionExists("mixture") && !bc_options.
        optionExists("rho")) {
          errorMessage="Must specify either 'rho' or 'mixture'." ;
          return false ;
        }
        if(!bc_options.optionExists("v") && !bc_options.optionExists("mdot") &&
        !bc_options.optionExists("cartesianV") && !bc_options.optionExists
        ("axisymmetricV") && !bc_options.optionExists("mms")) {
          errorMessage="Must specify velocity or mass flow rate." ;
          return false ;
        }
        if(*flowRegime=="turbulent" && !bc_options.optionExists("k") &&
        !bc_options.optionExists("cartesianK") && !bc_options.optionExists
        ("axisymmetricK") && !bc_options.optionExists("mms")){
          errorMessage="Must specify turbulent kinetic energy (k)." ;
          return false ;
        }
        if(*flowRegime=="turbulent" && !bc_options.optionExists("omega") &&
        !bc_options.optionExists("cartesianOmega") && !bc_options.optionExists
        ("axisymmetricOmega") && !bc_options.optionExists("mms")){
          errorMessage="Must specify turbulent dissipation rate (omega)." ;
          return false ;
        }
        if((bc_options.optionExists("axisymmetricV") || bc_options.
        optionExists("axisymmetricK") || bc_options.optionExists
        ("axisymmetricOmega")) && !bc_options.optionExists("referenceFrame")){
          errorMessage="Must specify reference frame when using " ;
          errorMessage+="axisymmetric profiles." ; return false ;
        }
        return true ;
      }
      ostream &ErrorMessage(std::ostream &s) {
        s << errorMessage << endl ; return s ;
      }
      string VariablesChecked(fact_db &facts) {
        string s="rho,v,cartesianV,axisymmetricV,mdot,massFlux,flowDirection" ;
        s+=",k,cartesianK,axisymmetricK" ;
        s+=",omega,cartesianOmega,axisymmetricOmega,referenceFrame" ;
        return s ;
      }
  } ;

  register_BC<CheckIncompressibleInlet> registerCheckIncompressibleInlet ;

  // Checks subsonic inlet boundaries.
  class CheckSubsonicInlet : public BC_Check {
    private:
      string errorMessage ;
    public:
      CheckSubsonicInlet() : errorMessage("") {}
    public:
      string BoundaryConditions() { return "subsonicInlet" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        param<string> flowRegime=facts.get_fact("flowRegime") ;
        if(!bc_options.optionExists("massFlux") && !bc_options.
        optionExists("mdot") && !bc_options.optionExists("v") &&
        !bc_options.optionExists("cartesianV") && !bc_options.
        optionExists("axisymmetricV")) {
          errorMessage="Must specify either 'massFlux','mdot' or velocity." ;
          return false ;
        }
        if(!bc_options.optionExists("T") && !bc_options.
        optionExists("cartesianT") && !bc_options.optionExists
        ("axisymmetricT")) {
          errorMessage="Must specify temperature." ; return false ;
        }
        if(*flowRegime=="turbulent" && !bc_options.optionExists("k") &&
        !bc_options.optionExists("cartesianK") && !bc_options.optionExists
        ("axisymmetricK")){
          errorMessage="Must specify turbulent kinetic energy (k)." ;
          return false ;
        }
        if(*flowRegime=="turbulent" && !bc_options.optionExists("omega") &&
        !bc_options.optionExists("cartesianOmega") && !bc_options.optionExists
        ("axisymmetricOmega")){
          errorMessage="Must specify turbulent dissipation rate (omega)." ;
          return false ;
        }
        if((bc_options.optionExists("axisymmetricV") || bc_options.
        optionExists("axisymmetricK") || bc_options.optionExists
        ("axisymmetricOmega") || bc_options.optionExists("axisymmetricT")) &&
        !bc_options.optionExists("referenceFrame")){
          errorMessage="Must specify reference frame when using " ;
          errorMessage+="axisymmetric options." ; return false ;
        }
        return true ;
      }
      ostream &ErrorMessage(std::ostream &s) {
        s << errorMessage << endl ; return s;
      }
      string VariablesChecked(fact_db &facts) {
        string s="massFlux,mdot,flowDirection,v,cartesianV,axisymmetricV," ;
        s+="T,cartesianT,axisymmetricT,k,cartesianK,axisymmetricK,omega," ;
        s+="cartesianOmega,axisymmetricOmega,mixture,referenceFrame" ;
        return s ;
      }
  } ;

  register_BC<CheckSubsonicInlet> registerCheckSubsonicInlet ;

  // Checks supersonic inlet boundaries.
  class CheckSupersonicInlet : public BC_Check {
    private:
      string errorMessage ;
    public:
      CheckSupersonicInlet() : errorMessage("") {}
    public:
      string BoundaryConditions() { return "supersonicInlet" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        param<string> flowRegime=facts.get_fact("flowRegime") ;
        if(!bc_options.optionExists("v") && !bc_options.optionExists
        ("cartesianV") && !bc_options.optionExists("axisymmetricV")) {
          errorMessage="Must specify velocity." ; return false ;
        }
        if(bc_options.optionExists("axisymmetricV") && !bc_options.
        optionExists("referenceFrame")){
          errorMessage="Must specify reference frame when using " ;
          errorMessage+="axisymmetricV." ; return false ;
        }
        if(!bc_options.optionExists("p")) {
          errorMessage="Must specify pressure." ; return false ;
        }
        if(!bc_options.optionExists("T") && !bc_options.
        optionExists("cartesianT") && !bc_options.optionExists
        ("axisymmetricT")) {
          errorMessage="Must specify temperature." ; return false ;
        }
        if(*flowRegime=="turbulent" && !bc_options.optionExists("k") &&
        !bc_options.optionExists("cartesianK") && !bc_options.optionExists
        ("axisymmetricK")){
          errorMessage="Must specify turbulent kinetic energy (k)." ;
          return false ;
        }
        if(*flowRegime=="turbulent" && !bc_options.optionExists("omega") &&
        !bc_options.optionExists("cartesianOmega") && !bc_options.optionExists
        ("axisymmetricOmega")){
          errorMessage="Must specify turbulent dissipation rate (omega)." ;
          return false ;
        }
        if((bc_options.optionExists("axisymmetricV") || bc_options.
        optionExists("axisymmetricK") || bc_options.optionExists
        ("axisymmetricOmega") || bc_options.optionExists("axisymmetricT")) &&
        !bc_options.optionExists("referenceFrame")){
          errorMessage="Must specify reference frame when using " ;
          errorMessage+="axisymmetric options." ; return false ;
        }
        return true ;
      }
      ostream &ErrorMessage(std::ostream &s) {
        s << errorMessage << endl ; return s;
      }
      string VariablesChecked(fact_db &facts) {
        string s="v,cartesianV,axisymmetricV,p,T,cartesianT,axisymmetricT," ;
        s+="k,cartesianK,axisymmetricK,omega,cartesianOmega," ;
        s+="axisymmetricOmega,mixture,referenceFrame" ;
        return s ;
      }
  } ;

  register_BC<CheckSupersonicInlet> registerCheckSupersonicInlet ;

  // Checks total pressure inlet boundaries.
  class CheckTotalPressureInlet : public BC_Check {
    private:
      string errorMessage ;
    public:
      CheckTotalPressureInlet() : errorMessage("") {}
    public:
      string BoundaryConditions() { return "totalPressureInlet" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        param<string> flowRegime=facts.get_fact("flowRegime"),
          flowCompressibility=facts.get_fact("flowCompressibility") ;
        if(!bc_options.optionExists("p0")) {
          errorMessage="Must specify total pressure." ; return false ;
        }
        if(!bc_options.optionExists("flowDirection")) {
          errorMessage="Must specify flow direction." ; return false ;
        }
        if(*flowCompressibility=="compressible" && !bc_options.optionExists
        ("T") && !bc_options.optionExists("cartesianT") && !bc_options.
        optionExists("axisymmetricT") && !bc_options.optionExists("T0")) {
          errorMessage="Must specify temperature or total temperature." ;
          return false ;
        }
        if(*flowRegime=="turbulent" && !bc_options.optionExists("k") &&
        !bc_options.optionExists("cartesianK") && !bc_options.optionExists
        ("axisymmetricK")){
          errorMessage="Must specify turbulent kinetic energy (k)." ;
          return false ;
        }
        if(*flowRegime=="turbulent" && !bc_options.optionExists("omega") &&
        !bc_options.optionExists("cartesianOmega") && !bc_options.optionExists
        ("axisymmetricOmega")){
          errorMessage="Must specify turbulent dissipation rate (omega)." ;
          return false ;
        }
        if((bc_options.optionExists("axisymmetricK") || bc_options.
        optionExists("axisymmetricOmega") || bc_options.optionExists
        ("axisymmetricT")) && !bc_options.optionExists("referenceFrame")){
          errorMessage="Must specify reference frame when using " ;
          errorMessage+="axisymmetric options." ; return false ;
        }
        return true ;
      }
      ostream &ErrorMessage(std::ostream &s) {
        s << errorMessage << endl ; return s;
      }
      string VariablesChecked(fact_db &facts) {
        string s="p0,flowDirection,T,T0,k,cartesianK,axisymmetricK," ;
        s+="omega,cartesianOmega,axisymmetricOmega,mixture" ; return s ;
      }
  } ;

  register_BC<CheckTotalPressureInlet> registerCheckTotalPressureInlet ;

  // Checks extrapolated pressure outlet boundaries.
  class CheckExtrapolatedPressureOutlet : public BC_Check {
    private:
      string errorMessage ;
    public:
      CheckExtrapolatedPressureOutlet() : errorMessage("") {}
    public:
      string BoundaryConditions() { return "extrapolatedPressureOutlet" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        return true ;
      }
      ostream &ErrorMessage(std::ostream &s) { return s; }
      string VariablesChecked(fact_db &facts) { return "entrainment" ; }
  } ;

  register_BC<CheckExtrapolatedPressureOutlet>
    registerCheckExtrapolatedPressureOutlet ;

  // Allows for extrapolated pressure outlet boundaries with no arguments.
  class CheckExtrapolatedPressureOutletNoArguments : public BC_Check {
    private:
      string errorMessage ;
    public:
      CheckExtrapolatedPressureOutletNoArguments() : errorMessage("") {}
    public:
      string BoundaryConditions() { return "extrapolatedPressureOutlet" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        return true ;
      }
      ostream &ErrorMessage(std::ostream &s) { return s; }
      string VariablesChecked(fact_db &facts) { return "" ; }
  } ;

  register_BC<CheckExtrapolatedPressureOutletNoArguments>
    registerCheckExtrapolatedPressureOutletNoArguments ;

  // Checks fixed pressure outlet boundaries.
  class CheckFixedPressureOutlet : public BC_Check {
    private:
      string errorMessage ;
    public:
      CheckFixedPressureOutlet() : errorMessage("") {}
    public:
      string BoundaryConditions() { return "fixedPressureOutlet" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        if(!bc_options.optionExists("p") && !bc_options.optionExists("pMean")) {
          errorMessage="Must specify pressure or mean pressure." ;
          return false ;
        }
        return true ;
      }
      ostream &ErrorMessage(std::ostream &s) {
        s << errorMessage << endl ; return s;
      }
      string VariablesChecked(fact_db &facts) {
        string s="p,pMean,entrainment,ID1,ID2,ID3,ID4,ID5" ;
        return s ;
      }
  } ;

  register_BC<CheckFixedPressureOutlet> registerCheckFixedPressureOutlet ;

  // Checks noslip boundaries. Removed the code that forced the user to
  // specifiy either "adiabatic" or temperature since this is now inconsistent
  // with the new "AP^n" boundary condition.
  class CheckNoslip : public BC_Check {
    private:
      string errorMessage ;
    public:
      CheckNoslip() : errorMessage("") {}
    public:
      string BoundaryConditions() { return "noslip" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        param<string> flowCompressibility=facts.get_fact
          ("flowCompressibility") ;
        if(*flowCompressibility=="compressible"){
          int count=0 ;
          if(bc_options.optionExists("adiabatic")) count++ ;
          if(bc_options.optionExists("T")) count++ ;
          if(bc_options.optionExists("cartesianT")) count++ ;
          if(bc_options.optionExists("Twall")) count++ ;
          if(bc_options.optionExists("qwall")) count++ ;
          if(bc_options.optionExists("Treservoir") && bc_options.optionExists("Rwall")) count++ ;
          if(count!=1) {
            errorMessage="Must specify only one of T, cartesianT, Twall, adiabatic, qwall" ;
            errorMessage+=" or the pair (Treservoir, Rwall)." ;
            return false ;
          }
          if(bc_options.optionExists("T")) {
            if(!check_scalar_units(bc_options,"T","K")) {
              errorMessage = "Wrong units for 'T'." ;
              return false ;
            }
          }
          if(bc_options.optionExists("Twall")) {
            if(!check_scalar_units(bc_options,"Twall","K")) {
              errorMessage = "Wrong units for 'Twall'." ;
              return false ;
            }
          }
          if(bc_options.optionExists("qwall")) {
            if(!check_scalar_units(bc_options,"qwall","watt/m/m")) {
              errorMessage = "Wrong units for 'qwall'." ;
              return false ;
            }
          }
        }
        return true ;
      }
      ostream &ErrorMessage(std::ostream &s) {
        s << errorMessage << endl ; return s;
      }
      string VariablesChecked(fact_db &facts) {
        param<string> flowCompressibility=facts.get_fact
          ("flowCompressibility") ;
        string s="wallFunction" ; s+=",referenceFrame" ; s+=",momentCenter" ;
        if(*flowCompressibility=="compressible")
          s+=",adiabatic,T,cartesianT,Twall,qwall,Treservoir,Rwall" ;
        return s ;
      }
  } ;

  register_BC<CheckNoslip> registerCheckNoslip ;

  // Allows for a noslip boundary with no arguments.
  class CheckNoslipNoArguments : public BC_Check {
    private:
      string errorMessage ;
    public:
      CheckNoslipNoArguments() : errorMessage("") {}
    public:
      string BoundaryConditions() { return "noslip" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        return true ;
      }
      ostream &ErrorMessage(std::ostream &s) {
        s << errorMessage << endl ; return s;
      }
      string VariablesChecked(fact_db &facts) { return "" ; }
  } ;

  register_BC<CheckNoslipNoArguments> registerCheckNoslipNoArguments ;

  // New interpolated interface boundary condition.
  class CheckInterface : public BC_Check {
    private:
      string errorMessage ;
    public:
      CheckInterface() : errorMessage("") {}
    public:
      string BoundaryConditions() { return "interface" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        return true ;
      }
      ostream &ErrorMessage(std::ostream &s) {
        s << errorMessage << endl ; return s;
      }
      string VariablesChecked(fact_db &facts) { return "" ; }
  } ;

  register_BC<CheckInterface> registerCheckInterface ;

}
