#ifndef ROTORDEFINITION_H
#define ROTORDEFINITION_H 
                                                                                
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
#include "referenceFrame.h"
#include "sciTypes.h"
                                                                                
namespace streamUns {

  // Class for rotor options.
  class RotorOptions : public options_list {
    public:
      RotorOptions() : options_list ("refFrameNum:center:radius:length:dpdz") {}
  } ;

  // Class which defines a rotor.
  class Rotor {
    private:
      real Lref ;
      unsigned int referenceFrameNumber ;
      vect3d center ;
      real radius,length,dpdzA,dpdzB ;
    private:
      bool Inside(const vect3d &coord) const {
        return false ;
      }
    public:
      Rotor(real Lref=0.0) : Lref(Lref) {}
    public:
      real DPDZA() const { return dpdzA ; }
      real DPDZB() const { return dpdzB ; }
      bool InReferenceFrame(const vect3d& cellCenter,const vector
        <ReferenceFrame>& referenceFrame) const ;
      unsigned int ReferenceFrameNumber() const {
        return referenceFrameNumber ;
      }
    public:
      void Input(const options_list &optionsList) ;
      istream& Input(istream &in) ;
      ostream& Print(ostream &out) const ;
    public:
      int BufferSize() const ;
      void PackBuffer(real *buffer,int size) ;
      void UnpackBuffer(real *buffer,int size) ;
  } ;
                                                                                
  // Output operator.
  inline ostream& operator<<(ostream &out, const Rotor &rotor) {
    return rotor.Print(out) ;
  }
                                                                                
  // Input operator.
  inline istream& operator>>(istream &in,Rotor &rotor) {
    return rotor.Input(in) ;
  }

  // Class for defining the rotors from a file .
  class RotorData {
    private:
      vector<Rotor> rotor ;
    public:
      istream& Input(istream &in) ;
      ostream& Print(ostream &out) const ;
    public:
      int BufferSize() const ;
      void PackBuffer(real *buffer,int size) ;
      void UnpackBuffer(real *buffer,int size) ;
    public:
      real DPDZA(const vect3d &cellCenter,const vector<ReferenceFrame>
        &referenceFrame) const ;
      real DPDZB(const vect3d &cellCenter,const vector<ReferenceFrame>
        &referenceFrame) const ;
      unsigned int ReferenceFrameNumber(const vect3d &cellCenter,const
        vector<ReferenceFrame> &referenceFrame) const ;
    public:
      istream& Read(istream &in,const real &Lref) ;
  } ;
                                                                                
  // Output operator.
  inline ostream& operator<<(ostream &out, const RotorData &rotorData) {
    return rotorData.Print(out) ;
  }
                                                                                
  // Input operator.
  inline istream& operator>>(istream &in,RotorData &rotorData) {
    return rotorData.Input(in) ;
  }

}

namespace Loci {

  // Class to serialize RotorData.
  class RotorDataConverter {
    private:
      streamUns::RotorData &ref ;
    public:
      explicit RotorDataConverter(streamUns::RotorData &ref) : ref(ref) {}
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
                                                                                
  // Schema converter for RotorData.
  template<> struct data_schema_traits<streamUns::RotorData> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
                                                                                
    typedef streamUns::real Converter_Base_Type ;
    typedef RotorDataConverter Converter_Type ;
  } ;

}

#endif
