#ifndef INITIALCONDITION_H
#define INITIALCONDITION_H

// Standard library includes.
#include <iostream>
#include <map>
#include <vector>
using std::cerr ;
using std::endl ;
using std::istream ;
using std::ostream ;
using std::map ;
using std::vector ;

// Loci includes.
#include <Loci.h>

// StreamUns includes.
#include "sciTypes.h"

namespace streamUns {

  // Class for initial condition options.
  class InitialConditionOptions : public options_list {
    public:
      InitialConditionOptions() : options_list("rho:v:T:p:k:omega:mixture") {}
  } ;

  // Class for the initial conditions.
  class InitialCondition {
    private:
      real density,temperature,pressure,k,omega ;
      vect3d velocity ;
      map<string,real> speciesMassFraction ;
      bool densityDefined,temperatureDefined,pressureDefined,velocityDefined ;
      bool speciesMassFractionDefined,kDefined,omegaDefined ;
    public:
      InitialCondition() : densityDefined(false),temperatureDefined(false),
        pressureDefined(false),velocityDefined(false),
        speciesMassFractionDefined(false),kDefined(false),omegaDefined(false) {}
      virtual ~InitialCondition() {}
    public:
      virtual void Input(const options_list &optionsList) ;
      virtual istream& Input(istream &in) ;
      virtual ostream& Print(ostream &out) const ;
    public:
      virtual int BufferSize() const ;
      virtual void PackBuffer(real *buffer,int size) ;
      virtual void UnpackBuffer(real *buffer,int size) ;
    public:
      real Density() const { return density ; }
      real K() const { return k ; }
      real Omega() const { return omega ; }
      real Pressure() const { return pressure ; }
      const map<string,real>& SpeciesMassFraction() const {
        return speciesMassFraction ;
      }
      real Temperature() const { return temperature ; }
      const vect3d& Velocity() const { return velocity ; }
    public:
      bool IsDensityDefined() const { return densityDefined ; }
      bool IsKDefined() const { return kDefined ; }
      bool IsOmegaDefined() const { return omegaDefined ; }
      bool IsPressureDefined() const { return pressureDefined ; }
      bool IsSpeciesMassFractionDefined() const {
        return speciesMassFractionDefined ;
      }
      bool IsTemperatureDefined() const { return temperatureDefined ; }
      bool IsVelocityDefined() const { return velocityDefined ; }
  } ;

  // Output operator.
  inline ostream& operator<<(ostream &out, const InitialCondition
  &initialCondition) {
    return initialCondition.Print(out) ;
  }

  // Input operator.
  inline istream& operator>>(istream &in,InitialCondition &initialCondition) {
    return initialCondition.Input(in) ;
  }

}

namespace Loci {

  // Class to serialize InitialCondition.
  class InitialConditionConverter {
    private:
      streamUns::InitialCondition &ref ;
    public:
      explicit InitialConditionConverter(streamUns::InitialCondition &ref) :
        ref(ref) {}
    public:
      int getSize() const {
        return ref.BufferSize() ;
      }
      void getState(streamUns::real *buf,int &size) {
        size=getSize() ; ref.PackBuffer(buf,size) ;
      }
      void setState(streamUns::real *buf,int size) {
        ref.UnpackBuffer(buf,size) ;
      }
  } ;

  // Schema converter for InitialCondition.
  template<> struct data_schema_traits<streamUns::InitialCondition> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef streamUns::real Converter_Base_Type ;
    typedef InitialConditionConverter Converter_Type ;
  } ;

}

#endif
