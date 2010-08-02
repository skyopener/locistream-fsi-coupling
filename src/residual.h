#ifndef RESIDUAL_H
#define RESIDUAL_H
#include "sciTypes.h"

namespace streamUns {

  class ScalarResidual {
    public:
      real maxResidual,totalResidual ;
      vect3d maxResidualLocation ;
    public:
      ScalarResidual() : maxResidual(0.0),totalResidual(0.0),
        maxResidualLocation(0.0,0.0,0.0) {}
  } ;

  class ScalarResidualJoin {
    public:
      void operator()(ScalarResidual &a,const ScalarResidual &b){
        if(abs(b.maxResidual)>abs(a.maxResidual)){
          a.maxResidual=b.maxResidual ;
          a.maxResidualLocation=b.maxResidualLocation ;
        }
        a.totalResidual+=abs(b.totalResidual) ;
      }
  } ;

  class VectorResidual {
    public:
      vect3d maxResidual,maxResidualLocation,totalResidual ;
    public:
      VectorResidual() : maxResidual(0.0,0.0,0.0),
        maxResidualLocation(0.0,0.0,0.0),totalResidual(0.0,0.0,0.0) {}
  } ;

  class VectorResidualJoin {
    public:
      void operator()(VectorResidual &a,const VectorResidual &b){
        if(norm(b.maxResidual)>norm(a.maxResidual)){
          a.maxResidual=b.maxResidual ;
          a.maxResidualLocation=b.maxResidualLocation ;
        }
        a.totalResidual+=b.totalResidual ;
      }
  } ;
}

namespace Loci {
  template<> struct data_schema_traits<streamUns::ScalarResidual> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct=CompoundFactory(streamUns::ScalarResidual()) ;
      LOCI_INSERT_TYPE(ct,streamUns::ScalarResidual,maxResidual) ;
      LOCI_INSERT_TYPE(ct,streamUns::ScalarResidual,maxResidualLocation) ;
      LOCI_INSERT_TYPE(ct,streamUns::ScalarResidual,totalResidual) ;
      return DatatypeP(ct) ;
    }
  } ;

  template<> struct data_schema_traits<streamUns::VectorResidual> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct=CompoundFactory(streamUns::VectorResidual()) ;
      LOCI_INSERT_TYPE(ct,streamUns::VectorResidual,maxResidual) ;
      LOCI_INSERT_TYPE(ct,streamUns::VectorResidual,maxResidualLocation) ;
      LOCI_INSERT_TYPE(ct,streamUns::VectorResidual,totalResidual) ;
      return DatatypeP(ct) ;
    }
  } ;
}

#endif

