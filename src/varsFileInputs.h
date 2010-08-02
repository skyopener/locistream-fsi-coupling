#ifndef VARSFILEINPUTS_H
#define VARSFILEINPUTS_H
                                                                                
// Standard library includes.
#include <iostream>
#include <map>
using std::cerr ;
using std::endl ;

// Loci includes.
#include <Loci.h>
                                                                                
// StreamUns includes.
#include "sciTypes.h"
                                                                                
namespace streamUns {

  // Class for momentum equation options.
  class MomentumEquationOptions : public options_list {
    public:
      MomentumEquationOptions() :
        options_list("linearSolver:relaxationFactor:maxIterations") {}
  } ;

  class PressureEquationOptions : public options_list {
    public:
      PressureEquationOptions() : options_list
        ("linearSolver:relaxationFactor:maxIterations:incompressibleForm:numStages") {}
  } ;

  class EnergyEquationOptions : public options_list {
    public:
      EnergyEquationOptions() : options_list
        ("linearSolver:relaxationFactor:maxIterations") {}
  } ;

  class TurbulenceEquationOptions : public options_list {
    public:
      TurbulenceEquationOptions() : options_list
        ("model:des:linearSolver:relaxationFactor:maxIterations") {}
  } ;

  class SpeciesEquationOptions : public options_list {
    public:
      SpeciesEquationOptions() : options_list
        ("linearSolver:relaxationFactor:maxIterations") {}
  } ;

  class probe_options : public options_list {
    public:
      probe_options() : options_list ("probe0:probe1:probe2:probe3:probe4:probe5:probe6:probe7:probe8:probe9") {} ;
  } ;

  class GravityOptions : public options_list {
    public:
      GravityOptions() : options_list("rhoref:g") {} ;
  } ;

  class IgnitionOptions : public options_list {
    public:
      IgnitionOptions() :
        options_list("heat:position:delta:nstart:nstop:Tmax:axis:dz:rmin:rmax") {} ;
  } ;

  class PerturbationOptions : public options_list {
    public:
      PerturbationOptions() : options_list("position:delta:nstart:nstop:dv")
        {} ;
  } ;

  class ReferenceValue : public options_list {
    public:
      ReferenceValue() : options_list("L:rho:v:h:k:omega") {} ;
  } ;
}

namespace Loci {

  template<> struct data_schema_traits<streamUns::MomentumEquationOptions> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<streamUns::MomentumEquationOptions>
      Converter_Type ;
  } ;

  template<> struct data_schema_traits<streamUns::PressureEquationOptions> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<streamUns::PressureEquationOptions>
      Converter_Type ;
  } ;

  template<> struct data_schema_traits<streamUns::EnergyEquationOptions> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<streamUns::EnergyEquationOptions>
      Converter_Type ;
  } ;

  template<> struct data_schema_traits<streamUns::TurbulenceEquationOptions> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<streamUns::TurbulenceEquationOptions>
      Converter_Type ;
  } ;

  template<> struct data_schema_traits<streamUns::SpeciesEquationOptions> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<streamUns::SpeciesEquationOptions>
      Converter_Type ;
  } ;

  template<> struct data_schema_traits<streamUns::probe_options> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<streamUns::probe_options> Converter_Type ;
  } ;

  template<> struct data_schema_traits<streamUns::GravityOptions> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<streamUns::GravityOptions> Converter_Type ;
  } ;

  template<> struct data_schema_traits<streamUns::IgnitionOptions> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<streamUns::IgnitionOptions> Converter_Type ;
  } ;

  template<> struct data_schema_traits<streamUns::PerturbationOptions> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<streamUns::PerturbationOptions>
      Converter_Type ;
  } ;

  template<> struct data_schema_traits<streamUns::ReferenceValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<streamUns::ReferenceValue>
      Converter_Type ;
  } ;
}

#endif
