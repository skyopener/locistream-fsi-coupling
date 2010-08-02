#ifndef INTEGRATE_H
#define INTEGRATE_H

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

  // Class used to hold species masses.
  class StandardVector {
    public:
      vector<real> data ;
    public:
      StandardVector(unsigned int size=0) : data(size,0.0) {}
      real Data(unsigned int i) const { return data[i] ; }
      unsigned int Size() const { return data.size() ; }
    public:

      // Returns the serialized buffer size.
      int BufferSize() const {
        int bufferSize=1+data.size() ; return bufferSize ;
      }

     // Packs the data into a buffer.
     void PackBuffer(real *buffer,int size) {
       int i=0 ;
       buffer[i++]=data.size() ;
       for(unsigned int j=0;j<data.size();++j) buffer[i++]=data[j] ;
     }

     // Prints the contents of the vector.
     ostream& Print(ostream &out) const { return out ; }

     // Unpacks the data from a buffer.
     void UnpackBuffer(real *buffer,int size) {
       int i=0 ; int n=int(buffer[i++]) ; data=vector<real>(n) ;
       for(int j=0;j<n;++j) data[j]=buffer[i++] ;
     }
  } ;

  // Output operator.
  inline ostream& operator<<(ostream &out, const StandardVector &v) {
    return v.Print(out) ;
  }
                                                                                
  // Input operator.
  inline istream& operator>>(istream &in,StandardVector &v) {
    return in ;
  }

  // Class used to join species masses.
  class StandardVectorJoin {
    public:
      void operator()(StandardVector &a,const StandardVector &b){
        for(unsigned int i=0;i<a.data.size();++i) a.data[i]+=b.data[i] ;
      }
  } ;
}

namespace Loci {
                                                                                
  // Class to serialize StandardVector.
  class StandardVectorConverter {
    private:
      streamUns::StandardVector &ref ;
    public:
      explicit StandardVectorConverter(streamUns::StandardVector &ref) :
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
                                                                                
  // Schema converter for StandardVector.
  template<> struct data_schema_traits<streamUns::StandardVector> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
                                                                                
    typedef streamUns::real Converter_Base_Type ;
    typedef StandardVectorConverter Converter_Type ;
  } ;
                                                                                
}

#endif
