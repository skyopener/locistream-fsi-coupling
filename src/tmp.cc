// Loci includes.
#include <Loci.h>
using Loci::Area ;

// StreamUns includes.
#include "sciTypes.h"

namespace streamUns {

  // Unit rule to compute the total area of a boundary. The parameters specify
  // the boundary constraints.
  class TotalBoundaryAreaNewUnit : public unit_rule {
    private:
      param<real> totalArea ;
    public:

      // Define input and output.
      TotalBoundaryAreaNewUnit() {
        name_store("totalArea(X,Y)",totalArea) ;
        constraint("X,ref->Y_BCoption") ;
        output("totalArea(X,Y)") ;
      }

      // Initialize the value.
      void compute(const sequence &seq) { *totalArea=0.0 ; }
  } ;

  register_rule<TotalBoundaryAreaNewUnit> registerTotalBoundaryAreaNewUnit ;

  // Apply rule to compute the total area of a boundary. The parameters
  // specify the boundary constraints.
  class TotalBoundaryAreaNewApply : public apply_rule<param<real>,
  Loci::Summation<real> > {
    private:
      const_store<Area> area ;
      param<real> totalArea ;
    public:

      // Define input and output.
      TotalBoundaryAreaNewApply() {
        name_store("area",area) ;
        name_store("totalArea(X,Y)",totalArea) ;
        input("area") ;
        output("totalArea(X,Y)") ;
        constraint("X,ref->Y_BCoption") ;
      }

      // Add the face value to the total.
      void calculate(Entity face) { join(*totalArea,area[face].sada) ; }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TotalBoundaryAreaNewApply> registerTotalBoundaryAreaNewApply ;

  // Unit rule to sum extrapolated pressure times area.
  class BoundaryPressureAreaWithIDUnit : public unit_rule {
    private:
      param<real> boundaryPressureAreaSum ;
    public:

      // Define input and output.
      BoundaryPressureAreaWithIDUnit() {
        name_store("boundaryPressureAreaSum(X,Y)",boundaryPressureAreaSum) ;
        output("boundaryPressureAreaSum(X,Y)") ;
        constraint("X,ref->Y_BCoption") ;
      }

      // Initialize the value.
      void compute(const sequence &seq) { *boundaryPressureAreaSum=0.0 ; }
  } ;

  register_rule<BoundaryPressureAreaWithIDUnit>
    registerBoundaryPressureAreaWithIDUnit ;

  // Unit rule to sum extrapolated pressure times area.
  class BoundaryPressureAreaWithIDApply : public apply_rule<param<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> p ;
      const_store<Area> area ;
      param<real> boundaryPressureAreaSum ;
    public:

      // Define input and output.
      BoundaryPressureAreaWithIDApply() {
        name_store("ci",ci) ;
        name_store("p",p) ;
        name_store("area",area) ;
        name_store("boundaryPressureAreaSum(X,Y)",boundaryPressureAreaSum) ;
        input("ci->p,area") ;
        output("boundaryPressureAreaSum(X,Y)") ;
        constraint("X,ref->Y_BCoption") ;
      }

      // Add the face value to the total.
      void calculate(Entity face) {
        join(*boundaryPressureAreaSum,p[ci[face]]*area[face].sada) ;
      }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryPressureAreaWithIDApply>
    registerBoundaryPressureAreaWithIDApply ;

  // Rule for boundary faces with specified mean pressure. Assigns pressure
  // value to all boundary faces that have the property pMean_BC.
  class BoundaryMeanPressureWithIDSpecification : public pointwise_rule {
    private:
      const_Map ref,ci ;
      const_param<real> boundaryPressureAreaSum ;
      const_param<real> totalArea ;
      const_store<real> p ;
      const_store<real> pMean_BC ;
      store<real> p_f ;
    private:
      real pC ;
    public:

      // Define input and output.
      BoundaryMeanPressureWithIDSpecification(const char *id) {
        string bPAS=string("boundaryPressureAreaSum(fixedPressureOutlet_BC,")+
          id+")" ;
        string tA=string("totalArea(fixedPressureOutlet_BC,")+id+")" ;
        name_store("ref",ref) ;
        name_store("ci",ci) ;
        name_store(bPAS,boundaryPressureAreaSum) ;
        name_store(tA,totalArea) ;
        name_store("p",p) ;
        name_store("pMean_BC",pMean_BC) ;
        name_store("id::p_f",p_f) ;
        input(bPAS) ;
        input(tA) ;
        input("ref->pMean_BC,ci->p") ;
        output("id::p_f") ;
      }

      // Extrapolate and add correction for a face. Compute the pressure
      // correction which is added to all faces so that the mean pressure is
      // pMean. IMPORTANT: There is an implicit assumption here that pMean is
      // the same for all fixedPressureOutlet boundaries. Tried to put this in
      // compute(), but was getting seg. fault when using seq.begin() as the
      // first face in the sequence for multi-process runs.
      void calculate(Entity face) {
        pC=pMean_BC[ref[face]]-(*boundaryPressureAreaSum)/(*totalArea) ;
        p_f[face]=p[ci[face]]+pC ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  // Use a macro to make it easy to create new rules for each ID.
  #define BOUNDARY_MEAN_PRESSURE_ID(X) class BOUNDARY_MEAN_PRESSURE_##X :\
     public BoundaryMeanPressureWithIDSpecification {\
      public:\
      BOUNDARY_MEAN_PRESSURE_##X() :\
      BoundaryMeanPressureWithIDSpecification(#X) {}\
    };\
    register_rule<BOUNDARY_MEAN_PRESSURE_##X>\
    register_BOUNDARY_MEAN_PRESSURE_##X;

  // One rule for each ID
  BOUNDARY_MEAN_PRESSURE_ID(ID1) ;
  BOUNDARY_MEAN_PRESSURE_ID(ID2) ;
  BOUNDARY_MEAN_PRESSURE_ID(ID3) ;
  BOUNDARY_MEAN_PRESSURE_ID(ID4) ;
  BOUNDARY_MEAN_PRESSURE_ID(ID5) ;

}
