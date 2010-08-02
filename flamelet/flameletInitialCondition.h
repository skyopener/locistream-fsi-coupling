#ifndef FLAMELETINITIALCONDITION_H
#define FLAMELETINITIALCONDITION_H

// Standard library includes.
#include <iostream>
#include <vector>
using std::cerr ;
using std::endl ;
using std::istream ;
using std::ostream ;
using std::vector ;

// Loci includes.
#include <Loci.h>

// StreamUns includes.
#include "sciTypes.h"

namespace streamUns {

  // Class for initial condition options.
  class FlameletInitialConditionOptions : public options_list {
    public:
      FlameletInitialConditionOptions() : options_list("Z:Zvar") {}
  } ;

  // Class for the initial conditions.
  class FlameletInitialCondition {
    private:
      real Z,Zvar ;
      bool ZDefined,ZvarDefined ;
    public:
      FlameletInitialCondition() : ZDefined(false),ZvarDefined(false) {}
      virtual ~FlameletInitialCondition() {}
    public:
      virtual void Input(const options_list &optionsList) ;
      virtual istream& Input(istream &in) ;
      virtual ostream& Print(ostream &out) const ;
    public:
      virtual int BufferSize() const ;
      virtual void PackBuffer(real *buffer,int size) ;
      virtual void UnpackBuffer(real *buffer,int size) ;
    public:
      real z() const { return Z ; }
      real zvar() const { return Zvar ; }
    public:
      bool IsZDefined() const { return ZDefined ; }
      bool IsZvarDefined() const { return ZvarDefined ; }
  } ;

  // Output operator.
  inline ostream& operator<<(ostream &out, const FlameletInitialCondition
  &flameletInitialCondition) {
    return flameletInitialCondition.Print(out) ;
  }

  // Input operator.
  inline istream& operator>>(istream &in,FlameletInitialCondition &flameletInitialCondition) {
    return flameletInitialCondition.Input(in) ;
  }

}

namespace Loci {

  // Class to serialize InitialCondition.
  class FlameletInitialConditionConverter {
    private:
      streamUns::FlameletInitialCondition &ref ;
    public:
      explicit FlameletInitialConditionConverter(streamUns::FlameletInitialCondition &ref) :
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
  template<> struct data_schema_traits<streamUns::FlameletInitialCondition> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef streamUns::real Converter_Base_Type ;
    typedef FlameletInitialConditionConverter Converter_Type ;
  } ;

}

#endif
