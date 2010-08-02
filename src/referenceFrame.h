#ifndef REFERENCE_FRAME_H
#define REFERENCE_FRAME_H
#include "sciTypes.h"

namespace streamUns {

  class ReferenceFrame {
    public:
      real omega ;
      vect3d axisStart,axisEnd ;
    public:
      ReferenceFrame() : omega(0.0),axisStart(0.0,0.0,0.0),
        axisEnd(0.0,0.0,0.0) {}
  } ;

  inline ostream& operator<<(ostream &s,const ReferenceFrame &r) {
    s << r.omega << ' ' << r.axisStart << ' ' << r.axisEnd ; return s ;
  }

  inline istream& operator>>(istream &s,ReferenceFrame &r) {
    s >> r.omega >> r.axisStart >> r.axisEnd ; return s ;
  }
}

namespace Loci {
  template<> struct data_schema_traits<streamUns::ReferenceFrame> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(streamUns::ReferenceFrame()) ;
      LOCI_INSERT_TYPE(ct,streamUns::ReferenceFrame,omega) ;
      LOCI_INSERT_TYPE(ct,streamUns::ReferenceFrame,axisStart) ;
      LOCI_INSERT_TYPE(ct,streamUns::ReferenceFrame,axisEnd) ;
      return DatatypeP(ct) ;
    }
  } ;
}

#endif

