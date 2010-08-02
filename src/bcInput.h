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

//-----------------------------------------------------------------------------
// Data structure to hold a Cartesian BC specification.
//-----------------------------------------------------------------------------

  class CartesianBoundaryCondition {
    private:
      size_t numPoints,coordFlag ;
      vector<vect3d> coord ;
      vector<real> value ;
    public:
      const vect3d& Coord(size_t i) const { return coord[i] ; }
      size_t CoordFlag() const { return coordFlag ; }
      size_t NumPoints() const { return numPoints ; }
      real Value(size_t i) const { return value[i] ; }
    public:
      void ReadFile(string fileName) {
        ifstream in(fileName.c_str(),ios::in) ;
        if(in.fail()) {
          cerr << "Open failed on " << fileName.c_str() << endl ;
          Loci::Abort() ;
        }

        // Skip spaces.
        while(!in.eof() && isspace(in.peek())) in.get() ;

        // Read in the number of points on the boundary and the coordinate
        // flag which indicates which coordinate value is changing.
        in >> numPoints >> coordFlag ;
        if(numPoints<2){
          cerr << "Bad number of data points in " << fileName << endl ;
          Loci::Abort() ;
        }
        if(coordFlag<0 || coordFlag>2){
          cerr << "Bad coordinate flag in " << fileName << endl ;
          Loci::Abort() ;
        }

        // Read in the coordinates and values.
        coord=vector<vect3d>(numPoints),value=vector<real>(numPoints) ;
        for(size_t i=0;i<numPoints;++i){
          in >> coord[i] ; in >> value[i] ;
        }
      }
    public:
      virtual istream& Input(istream &in) { return in ; }
      virtual ostream& Print(ostream &out) const { return out ; }
    public:
      // Returns the size of the buffer.
      int BufferSize() const {
        return 2+4*numPoints ;
      }
      // Packs the data into a buffer.
      void PackBuffer(double *buffer,int size) const {
        int i=0 ; buffer[i++]=numPoints ; buffer[i++]=coordFlag ;
        for(size_t j=0;j<numPoints;++j){
          buffer[i++]=coord[j].x ; buffer[i++]=coord[j].y ; buffer[i++]=coord[j].z ;
          buffer[i++]=value[j] ;
        }
      }
      // Unpacks the data from a buffer.
      void UnpackBuffer(double *buffer,int size) {
        int i=0 ;
        numPoints=int(buffer[i++]+0.001) ; coordFlag=int(buffer[i++]+0.001) ;
        coord=vector<vect3d>(numPoints) ; value=vector<real>(numPoints) ;
        for(size_t j=0;j<numPoints;++j){
          coord[j].x=buffer[i++] ; coord[j].y=buffer[i++] ; coord[j].z=buffer[i++] ;
          value[j]=buffer[i++] ;
        }
      }
  } ;

  // Output operator.
  inline ostream &operator<<(ostream &s,const CartesianBoundaryCondition &c) {
    return c.Print(s) ;
  }

  // Input operator.
  inline istream &operator>>(istream &s,CartesianBoundaryCondition &c) {
    return c.Input(s) ;
  }
}

namespace Loci {

  // Class to serialize CartesianBoundaryCondition.
  class CartesianBoundaryConditionConverter {
    private:
      streamUns::CartesianBoundaryCondition &ref ;
    public:
      explicit CartesianBoundaryConditionConverter(streamUns::CartesianBoundaryCondition
        &ref) : ref(ref) {}
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

  // Schema converter for CartesianBoundaryCondition.
  template<> struct data_schema_traits<streamUns::CartesianBoundaryCondition> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef double Converter_Base_Type ;
    typedef CartesianBoundaryConditionConverter Converter_Type ;
  } ;

}
#endif
