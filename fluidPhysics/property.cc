// Loci includes.
#include <Loci.h>

#include "eos.h"
#include "property.h"

#include <string>
using std::string ;
#include <iostream>
using std::cerr ;
using std::endl ;
using std::cout ;
#include <fstream>
using std::ifstream ;

namespace fluidPhysics {

  
//-----------------------------------------------------------------------------
// Class PolynomialSegment.
//-----------------------------------------------------------------------------

  // Returns the size of the buffer.
  int PolynomialSegment::BufferSize() const {
    return 6 ;
  }

  // Packs the data into a buffer.
  void PolynomialSegment::PackBuffer(double *buffer,int size) const {
    int i=0 ; buffer[i++]=tMin ; buffer[i++]=tMax ;
    for(int j=0;j<4;++j) buffer[i++]=a[j] ;
  }

  // Unpacks the data from a buffer.
  void PolynomialSegment::UnpackBuffer(double *buffer,int size) {
    int i=0 ; tMin=buffer[i++] ; tMax=buffer[i++] ;
    for(int j=0;j<4;++j) a[j]=buffer[i++] ;
  }

//-----------------------------------------------------------------------------
// Class PiecewisePolynomialProperty.
//-----------------------------------------------------------------------------

  // Returns the size of the buffer.
  int PiecewisePolynomialProperty::BufferSize() const {
    int size=1 ;
    for(vector<PolynomialSegment>::const_iterator p=polynomialSegment.begin();
      p!=polynomialSegment.end();++p) size+=p->BufferSize() ;
    return size ;
  }

  // Packs the database into a buffer.
  void PiecewisePolynomialProperty::PackBuffer(double *buffer,int size) const {
    int i=0 ; buffer[i++]=polynomialSegment.size() ;
    for(vector<PolynomialSegment>::const_iterator p=polynomialSegment.begin();
      p!=polynomialSegment.end();++p) p->PackBuffer(buffer,size) ;
  }

  // Unpacks the database from a buffer.
  void PiecewisePolynomialProperty::UnpackBuffer(double *buffer,int size) {
    int i=0 ; polynomialSegment=vector<PolynomialSegment>((int)buffer[i++]) ;
    for(vector<PolynomialSegment>::iterator p=polynomialSegment.begin();
      p!=polynomialSegment.end();++p) p->UnpackBuffer(buffer,size) ;
  }

  // Returns the value for a single value of the independent variable. The
 // return at the end is just to make some compilers happy. Will never get
 // there.
  inline double PiecewisePolynomialProperty::Value(double t) const {
    if(t<polynomialSegment.front().TMin())
      return polynomialSegment.front().Value(polynomialSegment.front().TMin()) ;
    if(t>polynomialSegment.back().TMax())
      return polynomialSegment.back().Value(polynomialSegment.back().TMax()) ;
    for(unsigned int i=0;i<polynomialSegment.size();++i)
      if(polynomialSegment[i].InsideRange(t))
        return polynomialSegment[i].Value(t) ;
    return 0.0 ;
  }

//-----------------------------------------------------------------------------
// Class MappedPiecewisePolynomialProperty.
//-----------------------------------------------------------------------------

  // Returns the size of the buffer.
  int MappedPiecewisePolynomialProperty::BufferSize() const {
    int size=1 ;
    for(map<double,PiecewisePolynomialProperty>::const_iterator m=value.
    begin();m!=value.end();++m){
      size+=1+m->second.BufferSize() ;
    }
    return size ;
  }

  // Packs the database into a buffer.
  void MappedPiecewisePolynomialProperty::PackBuffer(double *buffer,int size) {
    int i=0 ; buffer[i++]=value.size() ;
    for(map<double,PiecewisePolynomialProperty>::const_iterator m=value.
    begin();m!=value.end();++m){
      buffer[i++]=m->first ; m->second.PackBuffer(buffer,size) ;
    }
  }

  // Unpacks the database from a buffer.
  void MappedPiecewisePolynomialProperty::UnpackBuffer(double *buffer,int size) {
    int i=0,mapSize=int(buffer[i++]) ;
    for(int j=0;j<mapSize;++j){
      double temp=buffer[i++] ;
      PiecewisePolynomialProperty p ; p.UnpackBuffer(buffer,size) ;
      value[temp]=p ;
    }
  }

  // Returns the value for a single value of the independent variable.
  inline double MappedPiecewisePolynomialProperty::Value(double p,double t) const {

    // Find the first map entry with key greater than "p".
    map<double,PiecewisePolynomialProperty>::const_iterator m=value.
      upper_bound(p) ;

    // Get the value. Returns lowest pressure range value for pressure below
    // map range. Similar for higher.
    if(m==value.begin()){
      return m->second.Value(t) ;
    }else if(m==value.end()){
      --m ; return m->second.Value(t) ;
    }else{
      double upperP=m->first,upperValue=m->second.Value(t) ;
      --m ; double lowerP=m->first,lowerValue=m->second.Value(t) ;
      return (1.0/(upperP-lowerP))*((p-lowerP)*upperValue+
        (upperP-p)*lowerValue) ;
    }
  }

//-----------------------------------------------------------------------------
// Class TransportPropertyDatabase.
//-----------------------------------------------------------------------------

  // Returns the size of the buffer.
  int TransportPropertyDatabase::BufferSize() const {
    int numSpecies=molecularMass.size() ;
    int size=1+numSpecies+numSpecies*numSpecies ;
    for(int i=0;i<numSpecies;++i){
      size+=viscosity[i].BufferSize() ; size+=conductivity[i].BufferSize() ;
    }
    return size ;
  }
  
  // Returns the mixture conductivity given pressure, temperature and
  // species mass fractions.
  double TransportPropertyDatabase::Conductivity(double p,double t,const double *y)
  const {

    // Simple calculation for a single species.
    int numSpecies=molecularMass.size() ;
    if(numSpecies==1) return conductivity[0].Value(p,t) ;

    // Convert mass fractions to mole fractions.
    double sum=0.0 ; vector<double> x(numSpecies) ;
    for(int i=0;i<numSpecies;++i){
      x[i]=y[i]/molecularMass[i] ; sum+=x[i] ;
    }
    double rsum=1.0/sum ; for(int i=0;i<numSpecies;++i) x[i]*=rsum ;

    // Compute the mixture viscosity using Wilke's rule.
    double mixtureConductivity=0.0 ;
    for(int i=0;i<numSpecies;++i){
      double sum=x[0]*phi[i] ;
      for(int j=1;j<numSpecies;++j) sum+=x[j]*phi[j*numSpecies+i] ;
      mixtureConductivity+=x[i]*conductivity[i].Value(p,t)/sum ;
    }
    return mixtureConductivity ;
  }

  // Reads data from an options list.
  void TransportPropertyDatabase::Input(const SpeciesTransportPropertyOptions
  &options) {
    string err("TRANSPORT PROPERTY DATABASE ERROR: ") ;

    // Entry for viscosity.
    MappedPiecewisePolynomialProperty tempMu ;

    // Read the viscosity.
    if(!options.optionExists("mu")){
      cerr << err << "mu not specified for species" << endl ; Loci::Abort() ;
    }
    options_list::arg_list argListMu ; options.getOption("mu",argListMu) ;
    unsigned int count=0 ;
    while(count<argListMu.size()){

      // New piecewise polynomial property.
      PiecewisePolynomialProperty property ;

      // Read a pressure value.
      double p = 0;
      if(argListMu[count].type_of()==Loci::FUNCTION){
        string name ; argListMu[count].get_value(name) ;
        if(name=="p"){
          options_list::arg_list argListP ;
          argListMu[count].get_value(argListP) ;
          if(argListP.size()!=1){
            cerr << err << "Too many pressure arguments in mu" << endl ;
            Loci::Abort() ;
          }else{
            if(argListP[0].type_of()!=Loci::REAL){
              cerr << err << "Bad pressure argument in mu" << endl ;
              Loci::Abort() ;
            }else{
              argListP[0].get_value(p) ;
            }
          }
        }else{
          cerr << err << "Bad pressure specification in mu" << endl ;
          Loci::Abort() ;
        }
      }else{
        cerr << err << "Bad pressure specification in mu" << endl ;
        Loci::Abort() ;
      }
      ++count ;

      // Read the first temperature value.
      double tMin=0 ;
      if(argListMu[count].type_of()==Loci::REAL){
        argListMu[count].get_value(tMin) ;
      }else{
        cerr << err << "Bad temperature specification in mu" << endl ;
        Loci::Abort() ;
      }
      ++count ;

      // Read polynomial and temperature pairs until we get to the
      // end for this pressure value.
      bool done=false ;
      while(!done){

        // See if we have exhausted the arguments for mu.
        if(count>=argListMu.size()) done=true ;

        // See if we have read all polynomials.
        if(!done && argListMu[count].type_of()==Loci::FUNCTION){
          string name ; argListMu[count].get_value(name) ;
          if(name=="p") done=true ;
        }

        // Read the next polynomial/temperature pair.
        if(!done){

          // Polynomial curve fit. No need to test for function since this
          // was just done above.
          double coeff[4] ;
          options_list::arg_list argListPoly ;
          argListMu[count].get_value(argListPoly) ;
          if(argListPoly.size()!=4){
            cerr << err << "Bad number of polynomial coefficients in mu"
              << endl ; Loci::Abort() ;
          }else{
            for(int i=0;i<4;++i){
              if(argListPoly[i].type_of()!=Loci::REAL){
                cerr << err << "Bad polynomial coefficient in mu" << endl ;
                Loci::Abort() ;
              }else{
                argListPoly[i].get_value(coeff[i]) ;
              }
            }
          }
          ++count ;

          // Temperature.
          double tMax=0 ;
          if(argListMu[count].type_of()==Loci::REAL){
            argListMu[count].get_value(tMax) ;
          }else{
            cerr << err << "Bad tMax specification in mu" << endl ;
            Loci::Abort() ;
          }
          ++count ;

          // New polynomial segment.
          property.AddSegment(PolynomialSegment(tMin,tMax,coeff)) ;
        }
      }

      // Save the piecewise polynomial segment.
      tempMu.AddProperty(p,property) ; 
    }

    // Save the property for viscosity.
    AddViscosityEntry(tempMu) ;

    // Entry for conductivity.
    MappedPiecewisePolynomialProperty tempK ;

    // Read the thermal conductivity.
    if(!options.optionExists("k")){
      cerr << err << "k not specified for species" << endl ; Loci::Abort() ;
    }
    options_list::arg_list argListK ; options.getOption("k",argListK) ;
    count=0 ;
    while(count<argListK.size()){

      // New piecewise polynomial property.
      PiecewisePolynomialProperty property ;

      // Read a pressure value.
      double p = 0 ;
      if(argListK[count].type_of()==Loci::FUNCTION){
        string name ; argListK[count].get_value(name) ;
        if(name=="p"){
          options_list::arg_list argListP ;
          argListK[count].get_value(argListP) ;
          if(argListP.size()!=1){
            cerr << err << "Too many pressure arguments in k" << endl ;
            Loci::Abort() ;
          }else{
            if(argListP[0].type_of()!=Loci::REAL){
              cerr << err << "Bad pressure argument in k" << endl ;
              Loci::Abort() ;
            }else{
              argListP[0].get_value(p) ;
            }
          }
        }else{
          cerr << err << "Bad pressure specification in k" << endl ;
          Loci::Abort() ;
        }
      }else{
        cerr << err << "Bad pressure specification in k" << endl ;
        Loci::Abort() ;
      }
      ++count ;

      // Read the first temperature value.
      double tMin = 0;
      if(argListK[count].type_of()==Loci::REAL){
        argListK[count].get_value(tMin) ;
      }else{
        cerr << err << "Bad temperature specification in k" << endl ;
        Loci::Abort() ;
      }
      ++count ;

      // Read polynomial and temperature pairs until we get to the
      // end for this pressure value.
      bool done=false ;
      while(!done){

        // See if we have exhausted the arguments for k.
        if(count>=argListK.size()) done=true ;

        // See if we have read all polynomials.
        if(!done && argListK[count].type_of()==Loci::FUNCTION){
          string name ; argListK[count].get_value(name) ;
          if(name=="p") done=true ;
        }

        // Read the next polynomial/temperature pair.
        if(!done){

          // Polynomial curve fit. No need to test for function since this
          // was just done above.
          double coeff[4] ;
          options_list::arg_list argListPoly ;
          argListK[count].get_value(argListPoly) ;
          if(argListPoly.size()!=4){
            cerr << err << "Bad number of polynomial coefficients in k"
              << endl ; Loci::Abort() ;
          }else{
            for(int i=0;i<4;++i){
              if(argListPoly[i].type_of()!=Loci::REAL){
                cerr << err << "Bad polynomial coefficient in k" << endl ;
                Loci::Abort() ;
              }else{
                argListPoly[i].get_value(coeff[i]) ;
              }
            }
          }
          ++count ;

          // Temperature.
          double tMax = 0;
          if(argListK[count].type_of()==Loci::REAL){
            argListK[count].get_value(tMax) ;
          }else{
            cerr << err << "Bad tMax specification in mu" << endl ;
            Loci::Abort() ;
          }
          ++count ;

          // New polynomial segment.
          property.AddSegment(PolynomialSegment(tMin,tMax,coeff)) ;
        }
      }

      // Save the piecewise polynomial segment.
      tempK.AddProperty(p,property) ; 
    }

    // Save the property for conductivity.
    AddConductivityEntry(tempK) ;
  }

  // Overridden virtual method called by fact database when a fact of type
  // InitialCondition is read.
  istream& TransportPropertyDatabase::Input(istream &in) {
    SpeciesTransportPropertyOptions options ;
    in >> options ; Input(options) ; return in ;
  }

  // Packs the database into a buffer.
  void TransportPropertyDatabase::PackBuffer(double *buffer,int size) {
    int i=0,numSpecies=molecularMass.size() ;
    buffer[i++]=double(numSpecies) ;
    for(int j=0;j<numSpecies;++j) buffer[i++]=molecularMass[j] ;
    for(int j=0;j<numSpecies*numSpecies;++j) buffer[i++]=phi[j] ;
    for(int j=0;j<numSpecies;++j){
      viscosity[j].PackBuffer(buffer,size) ;
      conductivity[j].PackBuffer(buffer,size) ;
    }
  }

  // Prints the data in the database.
  ostream& TransportPropertyDatabase::Print(ostream &out) const {
    return out ;
  }

  // Sets up additional data.
  void TransportPropertyDatabase::Setup() {

    // Compute the interaction parameter for the species pairs.
    int numSpecies=molecularMass.size() ;
    phi=new double[numSpecies*numSpecies] ;
    for(int j=0;j<numSpecies;++j) for(int i=0;i<numSpecies;++i)
      phi[j*numSpecies+i]=sqrt(molecularMass[j]/molecularMass[i]) ;
  }

  // Unpacks the database from a buffer.
  void TransportPropertyDatabase::UnpackBuffer(double *buffer,int size) {
    int i=0,numSpecies=int(buffer[i++]) ;
    molecularMass.clear() ; phi=new double[numSpecies*numSpecies] ;
    for(int j=0;j<numSpecies;++j) molecularMass.push_back(buffer[i++]) ;
    for(int j=0;j<numSpecies*numSpecies;++j) phi[j]=buffer[i++] ;
    for(int j=0;j<numSpecies;++j){
      viscosity[j].UnpackBuffer(buffer,size) ;
      conductivity[j].UnpackBuffer(buffer,size) ;
    }
  }

  // Returns the mixture viscosity given pressure, temperature and
  // species mass fractions.
  double TransportPropertyDatabase::Viscosity(double p,double t,const double *y) const {

    // Simple calculation for a single species.
    int numSpecies=molecularMass.size() ;
    if(numSpecies==1) return viscosity[0].Value(p,t) ;

    // Convert mass fractions to mole fractions.
    double sum=0.0 ; vector<double> x(numSpecies) ;
    for(int i=0;i<numSpecies;++i){
      x[i]=y[i]/molecularMass[i] ; sum+=x[i] ;
    }
    double rsum=1.0/sum ; for(int i=0;i<numSpecies;++i) x[i]*=rsum ;

    // Compute the mixture viscosity using Wilke's rule.
    double mixtureViscosity=0.0 ;
    for(int i=0;i<numSpecies;++i){
      double sum=x[0]*phi[i] ;
      for(int j=1;j<numSpecies;++j) sum+=x[j]*phi[j*numSpecies+i] ;
      mixtureViscosity+=x[i]*viscosity[i].Value(p,t)/sum ;
//cout << "i,mu: " << i << " " << viscosity[i].Value(p,t) << endl ;
    }
    return mixtureViscosity ;
  }

}

