#ifndef FSIMOVE_H
#define FSIMOVE_H
                                                                                
// Standard library includes.
#include <iostream>
using std::cerr ;
using std::endl ;
                                                                                
// Loci includes.
#include <Loci.h>
                                                                                
// StreamUns includes.
#include "sciTypes.h"

namespace streamUns {

  // Options for rigid body grid motion.
  class RigidBodyGridMotionOptions : public options_list {
    public:
      RigidBodyGridMotionOptions() : options_list("sMag:sDir:tDir:tMag:tFreq:tPhi:rCenter:rAxis:rAlphaBar:rMag:rFreq:rPhi") {} ;
  } ;

  // Class to specifiy a time-dependent rigid-body displacement.
  class RigidBodyDisplacement {
    public:
      real func ; // actually int
      real sConst,tMag,tFreq,tPhi,rAlphaBar,rMag,rFreq,rPhi ;
      vect3d tDir,rCenter,rAxis ;
    public:
      RigidBodyDisplacement() : sConst(0.0),tMag(0.0),tFreq(0.0),tPhi(0.0),
        rAlphaBar(0.0),rMag(0.0),rFreq(0.0),rPhi(0.0),tDir(0.0,1.0,0.0),
        rCenter(0.0,0.0,0.0),rAxis(0.0,0.0,-1.0),func(1) {}
    public:
      vect3d Value(real time,vect3d position) const ;
    public:
      virtual int BufferSize() const ;
      virtual void PackBuffer(real *buffer,int size) ;
      virtual void UnpackBuffer(real *buffer,int size) ;
  } ;
/*
  // Class to specifiy a flexible-body displacement.
  class FlexibleBodyDisplacement {
	private:
	  int Nb, Nbd1, d ;
	  real rbfNr, linearSolverMaxIterations ; // input as reals, but use them as ints inside relevant functions
	  real r, a, linearSolverTolerance ;
    public:
	  ublas::matrix<real,ublas::column_major> matrixQ ;
	  ublas::matrix<real,ublas::column_major> matrixB ;
	  ublas::matrix<real,ublas::column_major> matrixForceTop ;
	  ublas::matrix<real,ublas::column_major> matrixForceBottom ;
	  ublas::matrix<real,ublas::column_major> matrixForce ;
	  ublas::matrix<real,ublas::row_major> FSIrbfWeight ;
	public::
	  void setNb(int N) { Nb=N; d=3; Nbd1=Nb+1+d} ; } ;
	  void AssembleQ();
	  void getFSIrbfWeight() ;
      vect3d Value(vect3d position) const ;
    public:
      virtual int BufferSize() const ;
      virtual void PackBuffer(real *buffer,int size) ;
      virtual void UnpackBuffer(real *buffer,int size) ;
  } ;
*/
  // Class for boundary displacement specification options. The types are as
  // follows: 0-constant,1-rigid-body,2-FSI(TOP),3=FSI(BOTTOM). Normally I would do all this in terms of
  // base and derived classes, but since we would then need to have a base
  // class pointer as the data type to the Loci rule, this gets messy according
  // to Ed.
  class BoundaryDisplacement {
    private:
      int type ;
      vect3d constant ;
      RigidBodyDisplacement rigidBody ;
    public:
      BoundaryDisplacement() : type(0),constant(0.0,0.0,0.0) {}
    public:
      bool IsSolutionDependent() const { return (type>1) ; } // type=2 -> solution dependent
      void SetConstantValue(vect3d v) { constant=v ; }
      void SetRigidBodyValue(const RigidBodyDisplacement v) { rigidBody=v ; }
      void SetType(int t) { type=t ; }
	  	int GetType() { return type ; }
      vect3d Value(real time,vect3d position) const {
        if(type==0) { 
        	return constant  ;
		} else if(type==1) {
	        return rigidBody.Value(time,position) ;
		} else if(type==2) {
			return vect3d(0.0, 0.0, 0.0) ;
		}
				//  return flexibleBody.Value(position) ;
      	//}
      }
    public:
      virtual istream& Input(istream &in) ;
      virtual ostream& Print(ostream &out) const ;
    public:
      virtual int BufferSize() const ;
      virtual void PackBuffer(real *buffer,int size) ;
      virtual void UnpackBuffer(real *buffer,int size) ;

  } ;

  // Output operator.
  inline ostream& operator<<(ostream &out, const BoundaryDisplacement &boundaryDisplacement) {
    return boundaryDisplacement.Print(out) ;
  }

  // Input operator.
  inline istream& operator>>(istream &in,BoundaryDisplacement &boundaryDisplacement) {
    return boundaryDisplacement.Input(in) ;
  }

}


namespace Loci {

  // Class to serialize BoundaryDisplacement.
  class BoundaryDisplacementConverter {
    private:
      streamUns::BoundaryDisplacement &ref ;
    public:
      explicit BoundaryDisplacementConverter(streamUns::BoundaryDisplacement
        &ref) : ref(ref) {}
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
                                                                                
  // Schema converter for BoundaryDisplacement.
  template<> struct data_schema_traits<streamUns::BoundaryDisplacement> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
                                                                                
    typedef streamUns::real Converter_Base_Type ;
    typedef BoundaryDisplacementConverter Converter_Type ;
  } ;

  template<> struct data_schema_traits<streamUns::RigidBodyGridMotionOptions> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<streamUns::RigidBodyGridMotionOptions>
      Converter_Type ;
  } ;
}

#endif
