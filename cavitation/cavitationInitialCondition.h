#ifndef CAVITATIONINITIALCONDITION_H
#define CAVITATIONINITIALCONDITION_H

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
  class CavitationInitialConditionOptions : public options_list {
    public:
      CavitationInitialConditionOptions() : options_list("Alpha") {}
  } ;

  // Class for the initial conditions.
  class CavitationInitialCondition {
    private:
      real Alpha ;
      bool AlphaDefined ;
    public:
      CavitationInitialCondition() : AlphaDefined(false) {}
      virtual ~CavitationInitialCondition() {}
    public:
      virtual void Input(const options_list &optionsList) ;
      virtual istream& Input(istream &in) ;
      virtual ostream& Print(ostream &out) const ;
    public:
      virtual int BufferSize() const ;
      virtual void PackBuffer(real *buffer,int size) ;
      virtual void UnpackBuffer(real *buffer,int size) ;
    public:
      real alpha() const { return Alpha ; }
    public:
      bool IsAlphaDefined() const { return AlphaDefined ; }
  } ;

  // Output operator.
  inline ostream& operator<<(ostream &out, const CavitationInitialCondition
  &cavitationInitialCondition) {
    return cavitationInitialCondition.Print(out) ;
  }

  // Input operator.
  inline istream& operator>>(istream &in,CavitationInitialCondition &cavitationInitialCondition) {
    return cavitationInitialCondition.Input(in) ;
  }

}

namespace Loci {

  // Class to serialize InitialCondition.
  class CavitationInitialConditionConverter {
    private:
      streamUns::CavitationInitialCondition &ref ;
    public:
      explicit CavitationInitialConditionConverter(streamUns::CavitationInitialCondition &ref) :
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
  template<> struct data_schema_traits<streamUns::CavitationInitialCondition> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef streamUns::real Converter_Base_Type ;
    typedef CavitationInitialConditionConverter Converter_Type ;
  } ;

}

#endif
