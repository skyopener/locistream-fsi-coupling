//-----------------------------------------------------------------------------
// Description: This file contains rules for implementing the axisymmetric
//   form of the governing equations.
//-----------------------------------------------------------------------------
                                                                                
// Standard library includes.
#include <vector>
using std::vector ;

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// CHEM includes.
#include "eos.h"

// StreamUns includes.
#include "const.h"
#include "referenceFrame.h"
#include "residual.h"
#include "sciTypes.h"
#include "varsFileInputs.h"
                                                                                
namespace streamUns {

//-----------------------------------------------------------------------------
// Rules for the momentum equation.

  // Rule to add components of viscous flux that are non-zero only for
  // compressible flows. This rule tacks on the axisymmetric contribution to
  // the divergence term.
  class DiffusiveFluxToVelocitySourceTermAxisymmetricInterior : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
      const_param<real> thetaParameter ;
      const_store<real> vol ;
      const_store<vect3d> v ;
      const_store<real> laminarViscosity ;
      const_store<real> eddyViscosity ;
      const_store<Area> area ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToVelocitySourceTermAxisymmetricInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("vol",vol) ;
        name_store("v",v) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("area",area) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("thetaParameter,(cl,cr)->(vol,v)") ;
        input("(cl,cr)->laminarViscosity,(cl,cr)->eddyViscosity,area") ;
        output("(cl,cr)->vSourceTerm") ;
        constraint("internalFaces,viscousFlow,compressibleFlow") ;
        constraint("axisymmetric") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real faceV=(v[cl[face]].y*vol[cr[face]]+v[cr[face]].y*vol[cl[face]])/
          (vol[cl[face]]+vol[cr[face]]) ;
        vect3d secondarySourceTerm=-((laminarViscosity[cl[face]]+eddyViscosity
          [cl[face]]+laminarViscosity[cr[face]]+eddyViscosity[cr[face]])/3.0*
          (*thetaParameter)*faceV*area[face].sada)*area[face].n ;
        vSourceTerm[cl[face]]+=secondarySourceTerm ;
        vSourceTerm[cr[face]]-=secondarySourceTerm ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToVelocitySourceTermAxisymmetricInterior>
    registerDiffusiveFluxToVelocitySourceTermAxisymmetricInterior ;

  // Rule to add components of viscous flux that are non-zero only for
  // compressible flows. This rule tacks on the axisymmetric contribution to
  // the divergence term.
  class DiffusiveFluxToVelocitySourceTermAxisymmetricBoundary : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map ci ;
      const_param<real> thetaParameter ;
      const_store<vect3d> v_f ;
      const_store<real> laminarViscosity_f ;
      const_store<real> eddyViscosity_f ;
      const_store<real> noWallFunction ;
      const_store<Area> area ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToVelocitySourceTermAxisymmetricBoundary() {
        name_store("ci",ci) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("v_f",v_f) ;
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        name_store("eddyViscosity_f",eddyViscosity_f) ;
        name_store("noWallFunction",noWallFunction) ;
        name_store("area",area) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("thetaParameter,v_f") ;
        input("laminarViscosity_f,eddyViscosity_f,noWallFunction,area") ;
        output("ci->vSourceTerm") ;
        constraint("boundaryVelocityDiffusion,viscousFlow") ;
        constraint("compressibleFlow,axisymmetric") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        vect3d secondarySourceTerm=-(2.0*(laminarViscosity_f[face]+
          eddyViscosity_f[face])/3.0*(*thetaParameter)*v_f[face].y*
          area[face].sada)*area[face].n*noWallFunction[face] ;
        vSourceTerm[ci[face]]+=secondarySourceTerm ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToVelocitySourceTermAxisymmetricBoundary>
    registerDiffusiveFluxToVelocitySourceTermAxisymmetricBoundary ;

  // Rule to add the axisymmetric volume source.
  class AxisymmetricToVelocitySourceTerm : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_param<real> thetaParameter ;
      const_store<real> cellRadius ;
      const_store<real> vol ;
      const_store<real> rho ;
      const_store<vect3d> v ;
      const_store<real> p ;
      const_store<real> laminarViscosity ;
      const_store<real> eddyViscosity ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      AxisymmetricToVelocitySourceTerm() {
        name_store("thetaParameter",thetaParameter) ;
        name_store("cellRadius",cellRadius) ;
        name_store("vol",vol) ;
        name_store("rho",rho) ;
        name_store("v",v) ;
        name_store("p",p) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("thetaParameter,cellRadius,vol,rho,v,p") ;
        input("laminarViscosity,eddyViscosity") ;
        output("vSourceTerm") ;
        constraint("geom_cells,viscousFlow,axisymmetric") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity cell) {
        real tauThetaTheta=(laminarViscosity[cell]+eddyViscosity[cell])*
          2.0*v[cell].y/cellRadius[cell] ;
        vect3d sourceTerm=vect3d(0.0,(p[cell]-tauThetaTheta)*vol[cell],0.0) ;
        vSourceTerm[cell]+=sourceTerm ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<AxisymmetricToVelocitySourceTerm>
    registerAxisymmetricToVelocitySourceTerm ;

  // Rule to add the axisymmetric volume source for inviscid flow.
  class AxisymmetricToVelocitySourceTermInviscid : public apply_rule
  <store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_param<real> thetaParameter ;
      const_store<real> vol ;
      const_store<real> rho ;
      const_store<vect3d> v ;
      const_store<real> p ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      AxisymmetricToVelocitySourceTermInviscid() {
        name_store("thetaParameter",thetaParameter) ;
        name_store("vol",vol) ;
        name_store("rho",rho) ;
        name_store("v",v) ;
        name_store("p",p) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("thetaParameter,vol,rho,v,p") ;
        output("vSourceTerm") ;
        constraint("geom_cells,inviscidFlow,axisymmetric") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity cell) {
        vect3d sourceTerm=vect3d(0.0,p[cell]*vol[cell],0.0) ;
        vSourceTerm[cell]+=sourceTerm ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<AxisymmetricToVelocitySourceTermInviscid>
    registerAxisymmetricToVelocitySourceTermInviscid ;

  // Rule to add the compressible part of the axisymmetric volume source.
  class AxisymmetricToVelocitySourceTermCompressible : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_param<real> thetaParameter ;
      const_store<real> cellRadius ;
      const_store<real> vol ;
      const_store<vect3d> v ;
      const_store<tens3d> vGradient ;
      const_store<real> laminarViscosity ;
      const_store<real> eddyViscosity ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      AxisymmetricToVelocitySourceTermCompressible() {
        name_store("thetaParameter",thetaParameter) ;
        name_store("cellRadius",cellRadius) ;
        name_store("vol",vol) ;
        name_store("v",v) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("thetaParameter,cellRadius,vol,v,gradv3d(v)") ;
        input("laminarViscosity,eddyViscosity") ;
        output("vSourceTerm") ;
        constraint("geom_cells,viscousFlow,compressibleFlow,axisymmetric") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity cell) {
        vSourceTerm[cell].y+=2.0*(laminarViscosity[cell]+eddyViscosity[cell])*
          (vGradient[cell].x.x+vGradient[cell].y.y+v[cell].y/cellRadius[cell])*
          vol[cell]/3.0 ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<AxisymmetricToVelocitySourceTermCompressible>
    registerAxisymmetricToVelocitySourceTermCompressible ;

//-----------------------------------------------------------------------------
// Rules for pressure correction.

  // Axisymmetric contribution to velocity correction. Only corrects the
  // radial velocity component.
  class CorrectVelocityAxisymmetricInterior : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_param<real> vRelaxationFactor ;
      const_store<real> pPrime ;
      const_store<real> vMainCoefficient ;
      const_store<real> vol ;
      store<vect3d> vCorrected ;
    public:
                                                                                
      // Define input and output.
      CorrectVelocityAxisymmetricInterior() {
        name_store("vRelaxationFactor",vRelaxationFactor) ;
        name_store("pPrime",pPrime) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        name_store("vol",vol) ;
        name_store("vCorrected",vCorrected) ;
        input("vRelaxationFactor,pPrime,vMainCoefficient,vol") ;
        output("vCorrected") ;
        constraint("geom_cells,axisymmetric") ;
      }
                                                                                
      // Correct the velocity for a cell.
      void calculate(Entity cell) {
        vCorrected[cell].y+=pPrime[cell]*vol[cell]/vMainCoefficient[cell]*
          (*vRelaxationFactor) ;
      }
                                                                                
      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<CorrectVelocityAxisymmetricInterior>
    registerCorrectVelocityAxisymmetricInterior ;

//-----------------------------------------------------------------------------
// Rules for the energy equation.

  // Rule to add the axisymmetric viscous stress contribution to the source
  // term for interior faces.
  class ViscousStressToTotalEnthalpySourceTermAxisymmetricInterior : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> vol ;
      const_store<vect3d> v ;
      const_store<real> laminarViscosity,eddyViscosity ;
      const_store<Area> area ;
      store<real> hSourceTerm ;
    public:
                                                                                
      // Define input and output.
      ViscousStressToTotalEnthalpySourceTermAxisymmetricInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("vol",vol) ;
        name_store("v",v) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("area",area) ;
        name_store("hSourceTerm",hSourceTerm) ;
        input("(cl,cr)->(vol,v,laminarViscosity,eddyViscosity),area") ;
        output("(cl,cr)->hSourceTerm") ;
        constraint("internalFaces,viscousFlow,axisymmetric") ;
      }
                                                                                
      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real faceV=0.5*(v[cl[face]].y*vol[cr[face]]+v[cr[face]].y*vol[cl[face]])
          /(vol[cl[face]]+vol[cr[face]]) ;
        real temp=-2.0*faceV/3.0 ;
        real sourceTerm=0.25*area[face].sada*(laminarViscosity[cl[face]]+
          eddyViscosity[cl[face]]+laminarViscosity[cr[face]]+eddyViscosity
          [cr[face]])*dot(temp*(v[cl[face]]+v[cr[face]]),area[face].n) ;
        hSourceTerm[cl[face]]+=sourceTerm ; hSourceTerm[cr[face]]-=sourceTerm ;
      }
                                                                                
      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<ViscousStressToTotalEnthalpySourceTermAxisymmetricInterior>
    registerViscousStressToTotalEnthalpySourceTermAxisymmetricInterior ;

  // Rule to add the axisymmetric viscous stress contribution to the source
  // term for interior faces.
  class ViscousStressToTotalEnthalpySourceTermAxisymmetricBoundary : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> laminarViscosity,eddyViscosity ;
      const_store<vect3d> v_f ;
      const_store<Area> area ;
      store<real> hSourceTerm ;
    public:
                                                                                
      // Define input and output.
      ViscousStressToTotalEnthalpySourceTermAxisymmetricBoundary() {
        name_store("ci",ci) ;
        name_store("v_f",v_f) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("area",area) ;
        name_store("hSourceTerm",hSourceTerm) ;
        input("ci->(laminarViscosity,eddyViscosity),v_f,area") ;
        output("ci->hSourceTerm") ;
        constraint("boundaryEnergyDiffusion,viscousFlow") ;
        constraint("axisymmetric") ;
      }
                                                                                
      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real temp=-2.0*v_f[face].y/3.0 ;
        real sourceTerm=0.5*area[face].sada*(laminarViscosity[ci[face]]+
          eddyViscosity[ci[face]])*dot(temp*v_f[face],area[face].n) ;
        hSourceTerm[ci[face]]+=sourceTerm ;
      }
                                                                                
      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<ViscousStressToTotalEnthalpySourceTermAxisymmetricBoundary>
    registerViscousStressToTotalEnthalpySourceTermAxisymmetricBoundary ;

//-----------------------------------------------------------------------------
// Rules for the turbulence equations.

  // Rule to compute the production term for the k equation for axisymmetric
  // incompressible flow.
  class ComputeKProductionAxisymmetricIncompressible : public pointwise_rule {
    private:
      const_store<real> rho,k,eddyViscosity ;
      const_store<vect3d> v ;
      const_store<tens3d> vGradient ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> kProduction ;
    public:
                                                                                
      // Define input and output.
      ComputeKProductionAxisymmetricIncompressible() {
        name_store("rho",rho) ;
        name_store("k",k) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("v",v) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("vol",vol) ;
        name_store("cellRadius",cellRadius) ;
        name_store("axisymmetric::kProduction",kProduction) ;
        input("rho,k,eddyViscosity,v,gradv3d(v),vol,cellRadius") ;
        output("axisymmetric::kProduction") ;
        constraint("incompressibleFlow,geom_cells,axisymmetric") ;
        constraint("menterBSLSSTTurbulenceModel") ;
      }
                                                                                
      // Add relaxation for a single cell.
      void calculate(Entity cell) {

        // Set the stress tensor, removing all derivatives which are
        // explicitly zero.
        tens3d stress ;
        real temp=2.0*rho[cell]*k[cell]/3.0 ;
        stress.x.x=2.0*eddyViscosity[cell]*vGradient[cell].x.x-temp ;
        stress.y.y=2.0*eddyViscosity[cell]*vGradient[cell].y.y-temp ;
        stress.z.z=2.0*eddyViscosity[cell]*v[cell].y/cellRadius[cell]-temp ;
        stress.x.y=stress.y.x=eddyViscosity[cell]*(vGradient[cell].x.y+
          vGradient[cell].y.x) ;
        stress.x.z=stress.z.x=eddyViscosity[cell]*vGradient[cell].z.x ;
        stress.y.z=stress.z.y=eddyViscosity[cell]*(vGradient[cell].z.y-
          v[cell].z/cellRadius[cell]) ;

        // Set the rate of strain tensor.
        tens3d strain ;
        strain.x.x=vGradient[cell].x.x ;
        strain.y.y=vGradient[cell].y.y ;
        strain.z.z=v[cell].y/cellRadius[cell] ;
        strain.x.y=strain.y.x=0.5*(vGradient[cell].x.y+vGradient[cell].y.x) ;
        strain.x.z=strain.z.x=0.5*vGradient[cell].z.x ;
        strain.y.z=strain.z.y=0.5*(vGradient[cell].z.y-v[cell].z/
          cellRadius[cell]) ;

        // Perform the scalar product.
        kProduction[cell]=ScalarProduct(stress,strain)*vol[cell] ;
      }
                                                                                
      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<ComputeKProductionAxisymmetricIncompressible>
    registerComputeKProductionAxisymmetricIncompressible ;

  // Rule to compute the production term for the k equation for incompressible
  // flow when using Menter's 2003 version of SST.
  class ComputeKProductionSST2003AxisymmetricIncompressible : public pointwise_rule {
    private:
      const_param<real> betaStar ;
      const_store<real> rho,k,omega,eddyViscosity ;
      const_store<vect3d> v ;
      const_store<tens3d> vGradient ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> kProduction ;
    public:

      // Define input and output.
      ComputeKProductionSST2003AxisymmetricIncompressible() {
        name_store("betaStar",betaStar) ;
        name_store("rho",rho) ;
        name_store("k",k) ;
        name_store("omega",omega) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("v",v) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("vol",vol) ;
        name_store("cellRadius",cellRadius) ;
        name_store("axisymmetric::kProduction",kProduction) ;
        input("betaStar,rho,k,omega,eddyViscosity,v,gradv3d(v)") ;
        input("cellRadius,vol") ;
        output("axisymmetric::kProduction") ;
        constraint("incompressibleFlow,geom_cells,axisymmetric") ;
        constraint("menterSST2003TurbulenceModel") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) {

        // Set the stress tensor, removing all derivatives which are
        // explicitly zero.
        tens3d stress ;
        real temp=2.0*rho[cell]*k[cell]/3.0 ;
        stress.x.x=2.0*eddyViscosity[cell]*vGradient[cell].x.x-temp ;
        stress.y.y=2.0*eddyViscosity[cell]*vGradient[cell].y.y-temp ;
        stress.z.z=2.0*eddyViscosity[cell]*v[cell].y/cellRadius[cell]-temp ;
        stress.x.y=stress.y.x=eddyViscosity[cell]*(vGradient[cell].x.y+
          vGradient[cell].y.x) ;
        stress.x.z=stress.z.x=eddyViscosity[cell]*vGradient[cell].z.x ;
        stress.y.z=stress.z.y=eddyViscosity[cell]*(vGradient[cell].z.y-
          v[cell].z/cellRadius[cell]) ;

        // Set the rate of strain tensor.
        tens3d strain ;
        strain.x.x=vGradient[cell].x.x ;
        strain.y.y=vGradient[cell].y.y ;
        strain.z.z=v[cell].y/cellRadius[cell] ;
        strain.x.y=strain.y.x=0.5*(vGradient[cell].x.y+vGradient[cell].y.x) ;
        strain.x.z=strain.z.x=0.5*vGradient[cell].z.x ;
        strain.y.z=strain.z.y=0.5*(vGradient[cell].z.y-v[cell].z/
          cellRadius[cell]) ;

        // Perform the scalar product.
        kProduction[cell]=min(ScalarProduct(stress,strain),10.0*(*betaStar)*
          rho[cell]*omega[cell]*k[cell])*vol[cell] ;
      }

      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeKProductionSST2003AxisymmetricIncompressible>
    registerComputeKProductionSST2003AxisymmetricIncompressible ;

  // Rule to compute the production term for the k equation for axisymmetric
  // compressible flow.
  class ComputeKProductionAxisymmetricCompressible : public pointwise_rule {
    private:
      const_store<real> rho,k,eddyViscosity ;
      const_store<vect3d> v ;
      const_store<tens3d> vGradient ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> kProduction ;
    public:
                                                                                
      // Define input and output.
      ComputeKProductionAxisymmetricCompressible() {
        name_store("rho",rho) ;
        name_store("k",k) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("v",v) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("vol",vol) ;
        name_store("cellRadius",cellRadius) ;
        name_store("axisymmetric::kProduction",kProduction) ;
        input("rho,k,eddyViscosity,v,gradv3d(v),vol,cellRadius") ;
        output("axisymmetric::kProduction") ;
        constraint("compressibleFlow,geom_cells,axisymmetric") ;
        constraint("menterBSLSSTTurbulenceModel") ;
      }
                                                                                
      // Add relaxation for a single cell.
      void calculate(Entity cell) {

        // Set the stress tensor, removing all derivatives which are
        // explicitly zero. Note that in the 'temp' term we precompte by the
        // eddy viscosity so that this cancels when we multiply everything
        // below by eddy viscosity.
        tens3d stress ;
        real temp=2.0*(eddyViscosity[cell]*(vGradient[cell].x.x+
          vGradient[cell].y.y+v[cell].y/cellRadius[cell])+rho[cell]*
          k[cell])/3.0 ;
        stress.x.x=2.0*eddyViscosity[cell]*vGradient[cell].x.x-temp ;
        stress.y.y=2.0*eddyViscosity[cell]*vGradient[cell].y.y-temp ;
        stress.z.z=2.0*eddyViscosity[cell]*v[cell].y/cellRadius[cell]-temp ;
        stress.x.y=stress.y.x=eddyViscosity[cell]*(vGradient[cell].x.y+
          vGradient[cell].y.x) ;
        stress.x.z=stress.z.x=eddyViscosity[cell]*vGradient[cell].z.x ;
        stress.y.z=stress.z.y=eddyViscosity[cell]*(vGradient[cell].z.y-
          v[cell].z/cellRadius[cell]) ;

        // Set the rate of strain tensor.
        tens3d strain ;
        strain.x.x=vGradient[cell].x.x ;
        strain.y.y=vGradient[cell].y.y ;
        strain.z.z=v[cell].y/cellRadius[cell] ;
        strain.x.y=strain.y.x=0.5*(vGradient[cell].x.y+vGradient[cell].y.x) ;
        strain.x.z=strain.z.x=0.5*vGradient[cell].z.x ;
        strain.y.z=strain.z.y=0.5*(vGradient[cell].z.y-v[cell].z/
          cellRadius[cell]) ;

        // Perform the scalar product.
        kProduction[cell]=ScalarProduct(stress,strain)*vol[cell] ;
      }
                                                                                
      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<ComputeKProductionAxisymmetricCompressible>
    registerComputeKProductionAxisymmetricCompressible ;

  // Rule to compute the production term for the k equation for axisymmetric
  // compressible flow.
  class ComputeKProductionSST2003AxisymmetricCompressible : public pointwise_rule {
    private:
      const_param<real> betaStar ;
      const_store<real> rho,k,omega,eddyViscosity ;
      const_store<vect3d> v ;
      const_store<tens3d> vGradient ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> kProduction ;
    public:
                                                                                
      // Define input and output.
      ComputeKProductionSST2003AxisymmetricCompressible() {
        name_store("betaStar",betaStar) ;
        name_store("rho",rho) ;
        name_store("k",k) ;
        name_store("omega",omega) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("v",v) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("vol",vol) ;
        name_store("cellRadius",cellRadius) ;
        name_store("axisymmetric::kProduction",kProduction) ;
        input("betaStar,rho,k,omega,eddyViscosity,v,gradv3d(v)") ;
        input("vol,cellRadius");
        output("axisymmetric::kProduction") ;
        constraint("compressibleFlow,geom_cells,axisymmetric") ;
        constraint("menterSST2003TurbulenceModel") ;
      }
                                                                                
      // Add relaxation for a single cell.
      void calculate(Entity cell) {

        // Set the stress tensor, removing all derivatives which are
        // explicitly zero. Note that in the 'temp' term we precompte by the
        // eddy viscosity so that this cancels when we multiply everything
        // below by eddy viscosity.
        tens3d stress ;
        real temp=2.0*(eddyViscosity[cell]*(vGradient[cell].x.x+
          vGradient[cell].y.y+v[cell].y/cellRadius[cell])+rho[cell]*
          k[cell])/3.0 ;
        stress.x.x=2.0*eddyViscosity[cell]*vGradient[cell].x.x-temp ;
        stress.y.y=2.0*eddyViscosity[cell]*vGradient[cell].y.y-temp ;
        stress.z.z=2.0*eddyViscosity[cell]*v[cell].y/cellRadius[cell]-temp ;
        stress.x.y=stress.y.x=eddyViscosity[cell]*(vGradient[cell].x.y+
          vGradient[cell].y.x) ;
        stress.x.z=stress.z.x=eddyViscosity[cell]*vGradient[cell].z.x ;
        stress.y.z=stress.z.y=eddyViscosity[cell]*(vGradient[cell].z.y-
          v[cell].z/cellRadius[cell]) ;

        // Set the rate of strain tensor.
        tens3d strain ;
        strain.x.x=vGradient[cell].x.x ;
        strain.y.y=vGradient[cell].y.y ;
        strain.z.z=v[cell].y/cellRadius[cell] ;
        strain.x.y=strain.y.x=0.5*(vGradient[cell].x.y+vGradient[cell].y.x) ;
        strain.x.z=strain.z.x=0.5*vGradient[cell].z.x ;
        strain.y.z=strain.z.y=0.5*(vGradient[cell].z.y-v[cell].z/
          cellRadius[cell]) ;

        // Perform the scalar product.
        kProduction[cell]=min(ScalarProduct(stress,strain),10.0*(*betaStar)*
          rho[cell]*k[cell]*omega[cell])*vol[cell] ;
      }
                                                                                
      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<ComputeKProductionSST2003AxisymmetricCompressible>
    registerComputeKProductionSST2003AxisymmetricCompressible ;
}
