#ifndef PROPERTY_H
#define PROPERTY_H

// Standard library includes.
#include <iostream>
using std::istream ;
using std::ostream ;
#include <map>
using std::map ;
#include <vector>
using std::vector ;


namespace fluidPhysics {

//-----------------------------------------------------------------------------
// Polynomial segment class.
//-----------------------------------------------------------------------------

  class PolynomialSegment {
    private:
      double tMin,tMax,a[4] ;
    public:
      bool InsideRange(double t) const {
        return (t>=tMin && t<=tMax)? true:false ;
      }
      double Value(double t) const { return a[0]+a[1]*t+a[2]*t*t+a[3]*t*t*t ; }
      double TMax() const { return tMax ; }
      double TMin() const { return tMin ; }
    public:
      PolynomialSegment() {}
      PolynomialSegment(double tMin,double tMax,const double *b) : tMin(tMin),
        tMax(tMax) { a[0]=b[0] ; a[1]=b[1] ; a[2]=b[2] ; a[3]=b[3] ; }
    public:
      int BufferSize() const ;
      void PackBuffer(double *buffer,int size) const ;
      void UnpackBuffer(double *buffer,int size) ;
  } ;

//-----------------------------------------------------------------------------
// Piecewise polynomial property class.
//-----------------------------------------------------------------------------

  class PiecewisePolynomialProperty {
    private:
      vector<PolynomialSegment> polynomialSegment ;
    public:
      double Value(double t) const ;
    public:
      void AddSegment(PolynomialSegment s) { polynomialSegment.push_back(s) ; }
      int NumSegment() const { return polynomialSegment.size() ; }
    public:
      int BufferSize() const ;
      void PackBuffer(double *buffer,int size) const ;
      void UnpackBuffer(double *buffer,int size) ;
  } ;

//-----------------------------------------------------------------------------
// Mapped piecewise polynomial property class.
//-----------------------------------------------------------------------------

  class MappedPiecewisePolynomialProperty {
    private:
      map<double,PiecewisePolynomialProperty> value ;
    public:
      void AddProperty(double p,PiecewisePolynomialProperty newProperty) {
        value[p]=newProperty ;
      }
    public:
      int BufferSize() const ;
      void PackBuffer(double *buffer,int size) ;
      void UnpackBuffer(double *buffer,int size) ;
    public:
      double Value(double p,double t) const ;
  } ;

//-----------------------------------------------------------------------------
// Transport property options for the species.
//-----------------------------------------------------------------------------

  class SpeciesTransportPropertyOptions : public options_list {
    public:
      SpeciesTransportPropertyOptions() : options_list("mu:k") {}
  } ;

//-----------------------------------------------------------------------------
// Transport property database class.
//-----------------------------------------------------------------------------

  class TransportPropertyDatabase {
    private:
      vector<double> molecularMass ;
      double *phi ;
      vector<MappedPiecewisePolynomialProperty> viscosity ;
      vector<MappedPiecewisePolynomialProperty> conductivity ;
    public:
      TransportPropertyDatabase() { phi=0 ; }
      virtual ~TransportPropertyDatabase() { if(phi) delete [] phi ; }
    public:
      virtual void Input(const SpeciesTransportPropertyOptions &options) ;
      virtual istream& Input(istream &in) ;
      virtual ostream& Print(ostream &out) const ;
    public:
      int BufferSize() const ;
      void PackBuffer(double *buffer,int size) ;
      void UnpackBuffer(double *buffer,int size) ;
    public:
      void AddConductivityEntry(MappedPiecewisePolynomialProperty c) {
        conductivity.push_back(c) ;
      }
      void AddMolecularMass(double m) { molecularMass.push_back(m) ; }
      void AddViscosityEntry(MappedPiecewisePolynomialProperty v) {
        viscosity.push_back(v) ;
      }
      void Setup() ;
    public:
      double Conductivity(double p,double t,const double *y) const ;
      double Viscosity(double p,double t,const double *y) const  ;
  } ;

  // Output operator.
  inline ostream &operator<<(ostream &s,const TransportPropertyDatabase
  &transportPropertyDatabase) {
    return transportPropertyDatabase.Print(s) ;
  }

  // Input operator.
  inline istream &operator>>(istream &s,TransportPropertyDatabase
  &transportPropertyDatabase) {
    return transportPropertyDatabase.Input(s) ;
  }

}

namespace Loci {

  // Class to serialize InitialCondition.
  class TransportPropertyDatabaseConverter {
    private:
    fluidPhysics::TransportPropertyDatabase &ref ;
    public:
      explicit TransportPropertyDatabaseConverter(fluidPhysics::
      TransportPropertyDatabase &ref) : ref(ref) {}
    public:
      int getSize() const {
        return ref.BufferSize() ;
      }
      void getState(double *buf,int &size) {
        size=getSize() ; ref.PackBuffer(buf,size) ;
      }
      void setState(double *buf,int size) {
        ref.UnpackBuffer(buf,size) ;
      }
  } ;
                                                                                
  // Schema converter for TransportPropertyDatabase.
  template<> struct data_schema_traits<fluidPhysics::TransportPropertyDatabase> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
                                                                                
    typedef double Converter_Base_Type ;
    typedef TransportPropertyDatabaseConverter Converter_Type ;
  } ;

}

#endif
