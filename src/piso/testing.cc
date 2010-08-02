//-----------------------------------------------------------------------------
// Description: This file contains rules for various test problems.
//
// Authors: Siddharth Thakur and Jeff Wright
//          Streamline Numerics Inc.
//          Gainesville, Florida 32609
//-----------------------------------------------------------------------------

// Standard library includes.
#include <string>
#include <sstream>
#include <stdlib.h>
#include <vector>
using std::string ;
using std::vector ;

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// StreamUns includes.
#include "gridReader/readGrid.h"
#include "sciTypes.h"

namespace streamUns {

//-----------------------------------------------------------------------------
// Unsteady cylinder problem.

  // Rule to write out the recirculation length behind a cylinder at Re=40.
  class CylinderData : public pointwise_rule {
    private:
      const_param<real> sTime ;
      const_Map cl,cr ;
      const_store<vect3d> v ;
      const_store<vect3d> cellCenter ;
      param<bool> OUTPUT ;
    private:
      real xMin,xMax,vXMin,vXMax ;
    public:

      // Define input and output.
      CylinderData() {
        name_store("stime{n}",sTime) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("v{n}",v) ;
        name_store("cellcenter{n}",cellCenter) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        input("stime{n},(cl,cr)->v{n},(cl,cr)->cellcenter{n}") ;
        output("OUTPUT{n}") ;
        conditional("do_print{n}") ;
        constraint("internalFaces{n},cylinder{n}") ;
      }

      // Write data for the cell if it is at the desired location.
      void calculate(Entity face) {
        real epsilon=1.0e-06 ;
        real xFace=0.5*(cellCenter[cl[face]].x+cellCenter[cr[face]].x) ;
        real yFace=0.5*(cellCenter[cl[face]].y+cellCenter[cr[face]].y) ;
        if(abs(yFace)<epsilon){
          if(xFace>0.5){
            real vX=0.5*(v[cl[face]].x+v[cr[face]].x) ;
            if(vX<0.0 && xFace>xMin){ xMin=xFace ; vXMin=vX ; }
            if(vX>0.0 && xFace<xMax){ xMax=xFace ; vXMax=vX ; }
          }
        }
      }

      // Write the data.
      virtual void compute(const sequence &seq) {
        xMin=0.5 ; xMax=20.0 ; do_loop(seq,this) ;
        real x=(-vXMin*xMax+vXMax*xMin)/(vXMax-vXMin) ;
        ofstream out("recirculationLength.txt",ios::app) ;
        out.setf(ios::scientific,ios::floatfield) ; out.precision(6) ;
        out << *sTime << " " << x-0.5 << endl ;
      }
  } ;

  register_rule<CylinderData> registerCylinderData ;

//-----------------------------------------------------------------------------
// Taylor decaying vortex problem.

  // This rule creates an "exact" boundary condition type with the allowable
  // options "taylor" and "rho".
  class CheckExactTaylor : public BC_Check {
    private:
      string errorMessage ;
    public:
      CheckExactTaylor() : errorMessage("") {}
    public:
      string BoundaryConditions() { return "exact" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        return true ;
      }
      ostream &ErrorMessage(std::ostream &s) {
        s << errorMessage << endl ; return s;
      }
      string VariablesChecked(fact_db &facts) {
        string s="taylor,rho" ; return s ;
      }
  } ;

  register_BC<CheckExactTaylor> registerCheckExactTaylor ;

  // Add additional constraints for the exact boundary condition.
  class ExactConstraints : public pointwise_rule {
    private:
      store<bool> boundaryMassFluxNotCorrected,extrapolatedPressure_BC ;
    public:

      // Define input and output.
      ExactConstraints() {
        name_store("boundaryMassFluxNotCorrected",boundaryMassFluxNotCorrected) ;
        name_store("extrapolatedPressure_BC",extrapolatedPressure_BC) ;
        output("boundaryMassFluxNotCorrected,extrapolatedPressure_BC") ;
        constraint("exact_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<ExactConstraints> registerExactConstraints ;

  // Sets the default Reynolds number.
  class DefaultReynoldsNumber : public default_rule {
    private:
      param<real> Re ;
    public:

      // Define input and output.
      DefaultReynoldsNumber() {
        name_store("Re",Re) ;
        output("Re") ;
        comments("Sets the default Reynolds number to 100.0.") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) { *Re=100.0 ; }
  } ;

  register_rule<DefaultReynoldsNumber> registerDefaultReynoldsNumber ;

  // Priority rule to assign exact initial condition.
  class InitialConditionTaylor : public pointwise_rule {
    private:
      const_param<real> timeStep ;
      const_param<real> Re ;
      const_store<vect3d> cellCenter ;
      store<vect3d> v_ic ;
      store<real> rho_ic,p_ic ;
    public:

      // Define input and output.
      InitialConditionTaylor() {
        name_store("timeStep",timeStep) ;
        name_store("Re",Re) ;
        name_store("cellcenter",cellCenter) ;
        name_store("taylor::rho_ic",rho_ic) ;
        name_store("taylor::v_ic",v_ic) ;
        name_store("taylor::p_ic",p_ic) ;
        input("timeStep,Re,cellcenter") ;
        constraint("taylor,geom_cells") ;
        output("taylor::rho_ic,taylor::v_ic,taylor::p_ic") ;
      }

      // Set velocity and pressure for a single cell.
      void calculate(Entity cell) {
        real x=cellCenter[cell].x,y=cellCenter[cell].y ;
        rho_ic[cell]=1.0 ;
        v_ic[cell]=vect3d(-cos(x)*sin(y),sin(x)*cos(y),0.0) ;
        p_ic[cell]=-0.25*(cos(2.0*x)+cos(2.0*y)) ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<InitialConditionTaylor> registerInitialConditionTaylor ;

  // Priority rule to assign exact solution for velocity at the boundary.
  class BoundarySpecificationVelocityTaylor : public pointwise_rule {
    private:
      const_param<real> oldSolutionTime ;
      const_param<real> timeStep ;
      const_param<real> Re ;
      const_store<vect3d> faceCenter ;
      const_store<Area> area ;
      store<vect3d> v_f ;
    public:

      // Define input and output.
      BoundarySpecificationVelocityTaylor() {
        name_store("Re",Re) ;
        name_store("stime",oldSolutionTime) ;
        name_store("timeStep",timeStep) ;
        name_store("facecenter",faceCenter) ;
        name_store("area",area) ;
        name_store("taylor::v_f",v_f) ;
        input("Re,stime,timeStep,facecenter,area") ;
        constraint("ref->taylor_BCoption") ;
        output("taylor::v_f") ;
      }

      // Set velocity, pressure and mass flux for a single face. Assumes
      // rho=1.0.
      void calculate(Entity face) {
        real x=faceCenter[face].x,y=faceCenter[face].y ;
        real time=*oldSolutionTime+*timeStep ;
        v_f[face]=vect3d(-cos(x)*sin(y)*exp(-2.0*time/(*Re)),sin(x)*cos(y)*
          exp(-2.0*time/(*Re)),0.0) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundarySpecificationVelocityTaylor>
    registerBoundarySpecificationVelocityTaylor ;

  // Priority rule to assign exact solution for mass flux at the boundary.
  // IMPORTANT NOTE: The constraint "ci->v" is required in order to ensure
  // that the rule for massFlux{n,it} is actually placed inside {n,it},
  // instead of being brought out in front of {n,it}, which Loci can do
  // since the inputs are not a function of {n,it}. Not having the rule
  // execute inside {n,it} causes scheduling error by Loci which results in
  // FPE's. Ed is going to fix the Loci scheduler so this will not be
  // required in the future.
  class BoundarySpecificationMassFluxTaylor : public pointwise_rule {
    private:
      const_param<real> oldSolutionTime ;
      const_param<real> timeStep ;
      const_param<real> Re ;
      const_store<vect3d> faceCenter ;
      const_store<Area> area ;
      store<real> massFlux ;
    public:

      // Define input and output.
      BoundarySpecificationMassFluxTaylor() {
        name_store("Re",Re) ;
        name_store("stime",oldSolutionTime) ;
        name_store("timeStep",timeStep) ;
        name_store("facecenter",faceCenter) ;
        name_store("area",area) ;
        name_store("massFlux",massFlux) ;
        input("Re,stime,timeStep,facecenter,area") ;
        constraint("ref->taylor_BCoption,ci->v") ;
        output("massFlux") ;
      }

      // Set velocity, pressure and mass flux for a single face. Assumes
      // rho=1.0.
      void calculate(Entity face) {
        real x=faceCenter[face].x,y=faceCenter[face].y ;
        real time=*oldSolutionTime+*timeStep ;
        vect3d v_f=vect3d(-cos(x)*sin(y)*exp(-2.0*time/(*Re)),sin(x)*cos(y)*
          exp(-2.0*time/(*Re)),0.0) ;
        massFlux[face]=dot(v_f,area[face].n)*area[face].sada ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundarySpecificationMassFluxTaylor>
    registerBoundarySpecificationMassFluxTaylor ;

  // Rule to write out the total kinetic energy in the domain.
  class TaylorKineticEnergy : public pointwise_rule {
    private:
      const_param<real> timeStep ;
      const_param<real> solutionTime ;
      const_param<real> Re ;
      const_store<vect3d> v ;
      const_store<real> vol ;
      const_store<vect3d> cellCenter ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      TaylorKineticEnergy() {
        name_store("timeStep{n}",timeStep) ;
        name_store("stime{n}",solutionTime) ;
        name_store("Re{n}",Re) ;
        name_store("v{n}",v) ;
        name_store("vol{n}",vol) ;
        name_store("cellcenter{n}",cellCenter) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        input("timeStep{n},stime{n},Re{n}") ;
        input("v{n},vol{n},cellcenter{n}") ;
        output("OUTPUT{n}") ;
        conditional("do_print{n}") ;
        constraint("geom_cells{n},taylor{n}") ;
      }

      // Write the data.
      virtual void compute(const sequence &seq) {

        // Write out numerical solution.
        ostringstream ossKE,ossKEExact ;
        ossKE << "output/KE.dat" ; ossKEExact << "output/KEExact.dat" ;
        string kEFile=ossKE.str(),kEExactFile=ossKEExact.str() ;
        ofstream outKE(kEFile.c_str(),ios::app) ;
        ofstream outKEExact(kEExactFile.c_str(),ios::app) ;
        outKE.setf(ios::scientific,ios::floatfield) ; outKE.precision(12) ;
        outKEExact.setf(ios::scientific,ios::floatfield) ;
        outKEExact.precision(12) ;
        cout << "Writing taylor kinetic energy." << endl ;
        real kineticEnergy=0.0 ;
        real t=*solutionTime ; //+*timeStep ;
        real kineticEnergyExact=acos(-1.0)*acos(-1.0)*exp(-4.0*t/(*Re)) ;
        for(sequence::const_iterator si=seq.begin();si!=seq.end();++si) {
          Entity cell=*si ;
//cout << "cell,vol,v,dKE: " << cell << " " << vol[cell] << " " << v[cell]
//  << " " << 0.5*vol[cell]*dot(v[cell],v[cell]) << endl ;
          kineticEnergy+=0.5*vol[cell]*dot(v[cell],v[cell]) ;
        }
        outKE << t << " " << kineticEnergy << endl ;
        outKEExact << t << " " << kineticEnergyExact << endl ;
      }
  } ;

  register_rule<TaylorKineticEnergy> registerTaylorKineticEnergy ;

  // Rule to write out the total kinetic energy in the domain.
  class WriteVelocityField : public pointwise_rule {
    private:
      const_store<vector3d<float> > v ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      WriteVelocityField() {
        name_store("cell2node_v3d(v){n}",v) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        input("cell2node_v3d(v){n}") ;
        output("OUTPUT{n}") ;
        conditional("do_print{n}") ;
        constraint("pos{n},taylor{n}") ;
      }

      // Write the data.
      virtual void compute(const sequence &seq) {
        for(sequence::const_iterator si=seq.begin();si!=seq.end();++si) {
          Entity n=*si ;
//cout << "node,v: " << n << " " << v[n] << endl ;
        }
      }
  } ;

//  register_rule<WriteVelocityField> registerWriteVelocityField ;

}
