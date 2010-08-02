//-----------------------------------------------------------------------------
// Description: This file contains rules for computing both volumetric and
//   surface integrated quantities.
//
// Author: Jeff Wright
//-----------------------------------------------------------------------------

// Standard library includes.
#include<fstream>
using std::ofstream ;
#include<vector>
using std::vector ;

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;

// StreamUns includes.
#include "integrate.h"
#include "sciTypes.h"

// Rules for integrating quantities over volumes.
namespace streamUns {

  // Rule to initialize the total volume.
  class InitializeTotalVolume : public unit_rule {
    private:
      param<real> totalVolume ;
    public:
                                                                                
      // Define input and output.
      InitializeTotalVolume() {
        name_store("totalVolume",totalVolume) ;
        output("totalVolume") ;
        constraint("geom_cells") ;
      }
                                                                                
      // Set the total volume zero.
      virtual void compute(const sequence &seq) { *totalVolume=0.0 ; }
  } ;
                                                                                
  register_rule<InitializeTotalVolume> registerInitializeTotalVolume ;

  // Rule to sum the total volume of all cells.
  class SumTotalVolume : public apply_rule<param<real>,Loci::Summation<real> > {
    private:
      const_store<real> vol ;
      param<real> totalVolume ;
    public:
                                                                                
      // Define input and output.
      SumTotalVolume() {
        name_store("vol",vol) ;
        name_store("totalVolume",totalVolume) ;
        input("vol") ;
        output("totalVolume") ;
        constraint("geom_cells") ;
      }

      // Add the cell volume.
      void calculate(Entity cell) { join(*totalVolume,vol[cell]) ; }
                                                                                
      // Sum over all cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<SumTotalVolume> registerSumTotalVolume ;

  // Rule to initialize the total mass.
  class InitializeTotalMass : public unit_rule {
    private:
      param<real> totalMass ;
    public:
                                                                                
      // Define input and output.
      InitializeTotalMass() {
        name_store("totalMass",totalMass) ;
        output("totalMass") ;
        constraint("geom_cells") ;
      }
                                                                                
      // Set the total mass to zero.
      virtual void compute(const sequence &seq) { *totalMass=0.0 ; }
  } ;
                                                                                
  register_rule<InitializeTotalMass> registerInitializeTotalMass ;

  // Rule to sum the total mass over all cells.
  class SumTotalMass : public apply_rule<param<real>,Loci::Summation<real> > {
    private:
      const_store<real> vol,rho ;
      param<real> totalMass ;
    public:
                                                                                
      // Define input and output.
      SumTotalMass() {
        name_store("vol",vol) ;
        name_store("rho",rho) ;
        name_store("totalMass",totalMass) ;
        input("vol,rho") ;
        output("totalMass") ;
        constraint("geom_cells") ;
      }

      // Add the cell mass.
      void calculate(Entity cell) { join(*totalMass,rho[cell]*vol[cell]) ; }
                                                                                
      // Sum over all cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<SumTotalMass> registerSumTotalMass ;

  // Rule to set zero total energy and enthalpy for incompressible flows.
  class TotalEnergyAndEnthalpyIncompressible : public singleton_rule {
    private:
      param<real> totalEnergy,totalEnthalpy ;
    public:

      // Define input and output.
      TotalEnergyAndEnthalpyIncompressible() {
        name_store("totalEnergy",totalEnergy) ;
        name_store("totalEnthalpy",totalEnthalpy) ;
        output("totalEnergy,totalEnthalpy") ;
        constraint("incompressibleFlow") ;
      }

      // Set the values.
      void compute(const sequence &seq) { *totalEnergy=*totalEnthalpy=0.0 ; }
  } ;

  register_rule<TotalEnergyAndEnthalpyIncompressible>
    registerTotalEnergyAndEnthalpyIncompressible ;

  // Rule to initialize the total energy.
  class InitializeTotalEnergy : public unit_rule {
    private:
      param<real> totalEnergy ;
    public:
                                                                                
      // Define input and output.
      InitializeTotalEnergy() {
        name_store("totalEnergy",totalEnergy) ;
        output("totalEnergy") ;
        constraint("compressibleFlow,geom_cells") ;
      }
                                                                                
      // Set the total energy to zero.
      virtual void compute(const sequence &seq) { *totalEnergy=0.0 ; }
  } ;
                                                                                
  register_rule<InitializeTotalEnergy> registerInitializeTotalEnergy ;

  // Rule to sum the total energy over all cells.
  class SumTotalEnergy : public apply_rule<param<real>,Loci::Summation<real> > {
    private:
      const_store<real> vol,rho,kineticEnergy,p,h ;
      param<real> totalEnergy ;
    public:
                                                                                
      // Define input and output.
      SumTotalEnergy() {
        name_store("vol",vol) ;
        name_store("rho",rho) ;
        name_store("kineticEnergy",kineticEnergy) ;
        name_store("p",p) ;
        name_store("h",h) ;
        name_store("totalEnergy",totalEnergy) ;
        input("vol,rho,kineticEnergy,p,h") ;
        output("totalEnergy") ;
        constraint("compressibleFlow,geom_cells") ;
      }

      // Add the cell energy. Must extract from total enthalpy
      void calculate(Entity cell) {
        join(*totalEnergy,(rho[cell]*(h[cell]-kineticEnergy[cell])-p[cell])*
          vol[cell]) ;
      }
                                                                                
      // Sum over all cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<SumTotalEnergy> registerSumTotalEnergy ;

  // Rule to initialize the total enthalpy.
  class InitializeTotalEnthalpy : public unit_rule {
    private:
      param<real> totalEnthalpy ;
    public:
                                                                                
      // Define input and output.
      InitializeTotalEnthalpy() {
        name_store("totalEnthalpy",totalEnthalpy) ;
        output("totalEnthalpy") ;
        constraint("compressibleFlow,geom_cells") ;
      }
                                                                                
      // Set the total enthalpy to zero.
      virtual void compute(const sequence &seq) { *totalEnthalpy=0.0 ; }
  } ;
                                                                                
  register_rule<InitializeTotalEnthalpy> registerInitializeTotalEnthalpy ;

  // Rule to sum the total enthalpy over all cells.
  class SumTotalEnthalpy : public apply_rule<param<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> vol,rho,kineticEnergy,h ;
      param<real> totalEnthalpy ;
    public:
                                                                                
      // Define input and output.
      SumTotalEnthalpy() {
        name_store("vol",vol) ;
        name_store("rho",rho) ;
        name_store("kineticEnergy",kineticEnergy) ;
        name_store("h",h) ;
        name_store("totalEnthalpy",totalEnthalpy) ;
        input("vol,rho,kineticEnergy,h") ;
        output("totalEnthalpy") ;
        constraint("compressibleFlow,geom_cells") ;
      }

      // Add the cell enthalpy. Must extract from total enthalpy
      void calculate(Entity cell) {
        join(*totalEnthalpy,rho[cell]*(h[cell]-kineticEnergy[cell])*vol[cell]) ;
      }
                                                                                
      // Sum over all cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<SumTotalEnthalpy> registerSumTotalEnthalpy ;

  // Rule to set zero total species masses for cases with no species transport.
  class TotalSpeciesMassNoSpeciesTransport : public singleton_rule {
    private:
      param<StandardVector> totalSpeciesMass ;
    public:

      // Define input and output.
      TotalSpeciesMassNoSpeciesTransport() {
        name_store("totalSpeciesMass",totalSpeciesMass) ;
        output("totalSpeciesMass") ;
        constraint("noSpeciesTransport") ;
      }

      // Set the values.
      void compute(const sequence &seq) {}
  } ;

  register_rule<TotalSpeciesMassNoSpeciesTransport>
    registerTotalSpeciesMassNoSpeciesTransport ;

  // Rule to initialize the total species masses.
  class InitializeTotalSpeciesMass : public unit_rule {
    private:
      const_param<int> numSpecies ;
      param<StandardVector> totalSpeciesMass ;
    public:
                                                                                
      // Define input and output.
      InitializeTotalSpeciesMass() {
        name_store("numSpecies",numSpecies) ;
        name_store("totalSpeciesMass",totalSpeciesMass) ;
        input("numSpecies") ;
        output("totalSpeciesMass") ;
        constraint("geom_cells,speciesTransport") ;
      }
                                                                                
      // Set the total species masses to zero.
      virtual void compute(const sequence &seq) {
        *totalSpeciesMass=StandardVector(*numSpecies) ;
      }
  } ;
                                                                                
  register_rule<InitializeTotalSpeciesMass> registerInitializeTotalSpeciesMass ;

  // Rule to sum the total species masses over all cells.
  class SumTotalSpeciesMass : public apply_rule<param<StandardVector>,
  StandardVectorJoin> {
    private:
      const_param<int> numSpecies ;
      const_store<real> vol,rho ;
      const_storeVec<real> y ;
      param<StandardVector> totalSpeciesMass ;
    public:
                                                                                
      // Define input and output.
      SumTotalSpeciesMass() {
        name_store("numSpecies",numSpecies) ;
        name_store("vol",vol) ;
        name_store("rho",rho) ;
        name_store("y",y) ;
        name_store("totalSpeciesMass",totalSpeciesMass) ;
        input("numSpecies,vol,rho,y") ;
        output("totalSpeciesMass") ;
        constraint("geom_cells,speciesTransport") ;
      }

      // Add the cell species mass.
      void calculate(Entity cell) {
        StandardVector temp(*numSpecies) ;
        for(int i=0;i<(*numSpecies);++i) temp.data[i]=rho[cell]*y[cell][i]*
          vol[cell] ;
        join(*totalSpeciesMass,temp) ;
      }
                                                                                
      // Sum over all cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<SumTotalSpeciesMass> registerSumTotalSpeciesMass ;

  // Rule to write out integrated volumetric quantities.
  class WriteVolumeData : public singleton_rule {
    private:
      const_param<EOS> eos ;
      const_param<real> totalVolume,totalMass,totalEnergy,totalEnthalpy ;
      const_param<StandardVector> totalSpeciesMass ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      WriteVolumeData() {
        name_store("eos{n}",eos) ;
        name_store("totalVolume{n}",totalVolume) ;
        name_store("totalMass{n}",totalMass) ;
        name_store("totalEnergy{n}",totalEnergy) ;
        name_store("totalEnthalpy{n}",totalEnthalpy) ;
        name_store("totalSpeciesMass{n}",totalSpeciesMass) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        input("eos{n}") ;
        input("totalVolume{n},totalMass{n},totalEnergy{n}") ;
        input("totalEnthalpy{n},totalSpeciesMass{n}") ;
        output("OUTPUT{n}") ;
        conditional("do_print{n}") ;
      }

      // Write out the data.
      void compute(const sequence &seq) {
        if(Loci::MPI_rank!=0) return ;
        cout.setf(ios::scientific,ios::floatfield) ; cout.precision(4) ;
        cout << "Integrated Volumetric Data (Complete Domain)" << endl ;
        cout << "  total volume = " << *totalVolume << " m^3" << endl ;
        cout << "  total mass = " << *totalMass << " kg" << endl ;
        cout << "  total energy = " << *totalEnergy << " J" << endl ;
        cout << "  total enthalpy = " << *totalEnthalpy << " J " << endl ;
        unsigned int numSpecies=totalSpeciesMass->Size() ;
        const vector<string> &speciesName=eos->speciesNames() ;
        if(numSpecies!=0){
          cout << "  total species masses: [" ;
          for(unsigned int i=0;i<numSpecies;++i){
            cout << speciesName[i] << " = " << totalSpeciesMass->Data(i) ;
            if(i!=numSpecies-1) cout << "," ;
          }
          cout << "] kg" << endl ;
        }
      }
  } ;

  register_rule<WriteVolumeData> registerWriteVolumeData ;

}

// Rules for integrating quantities over boundary surfaces.
namespace streamUns {

  // Unit rule to compute the total area of a boundary. The parameter specifies
  // the boundary constraint.
  class TotalBoundaryAreaUnit : public unit_rule {
    private:
      param<real> totalArea ;
    public:

      // Define input and output.
      TotalBoundaryAreaUnit() {
        name_store("totalArea(X)",totalArea) ;
        constraint("X") ;
        output("totalArea(X)") ;
      }

      // Initialize the value.
      void compute(const sequence &seq) { *totalArea=0.0 ; }
  } ;

  register_rule<TotalBoundaryAreaUnit> registerTotalBoundaryAreaUnit ;

  // Apply rule to compute the total area of a boundary. The parameter
  // specifies the boundary constraint.
  class TotalBoundaryAreaApply : public apply_rule<param<real>,
  Loci::Summation<real> > {
    private:
      const_store<Area> area ;
      param<real> totalArea ;
    public:

      // Define input and output.
      TotalBoundaryAreaApply() {
        name_store("area",area) ;
        name_store("totalArea(X)",totalArea) ;
        input("area") ;
        output("totalArea(X)") ;
        constraint("X") ;
      }

      // Add the face value to the total.
      void calculate(Entity face) { join(*totalArea,area[face].sada) ; }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TotalBoundaryAreaApply> registerTotalBoundaryAreaApply ;

  // Unit rule to compute the mass transfer through a boundary surface.
  // The parameter specifies the boundary constraint.
  class BoundaryMassTransferUnit : public unit_rule {
    private:
      param<real> massTransfer ;
    public:

      // Define input and output.
      BoundaryMassTransferUnit() {
        name_store("massTransfer(X)",massTransfer) ;
        output("massTransfer(X)") ;
        constraint("UNIVERSE") ;
      }

      // Initialize the value.
      void compute(const sequence &seq) { *massTransfer=0.0 ; }
  } ;

  register_rule<BoundaryMassTransferUnit> registerBoundaryMassTransferUnit ;

  // Apply rule to compute the mass transfer through a boundary surface.
  // The parameter specifies the boundary constraint.
  class BoundaryMassTransferApply : public apply_rule<param<real>,
  Loci::Summation<real> > {
    private:
      const_store<real> massFlux ;
      param<real> massTransfer ;
    public:

      // Define input and output.
      BoundaryMassTransferApply() {
        name_store("massFlux",massFlux) ;
        name_store("massTransfer(X)",massTransfer) ;
        input("massFlux") ;
        output("massTransfer(X)") ;
        constraint("X") ;
      }

      // Add the face value to the total.
      void calculate(Entity face) { join(*massTransfer,massFlux[face]) ; }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryMassTransferApply> registerBoundaryMassTransferApply ;

  // Apply rule for incompressible flows.
  class BoundaryEnergyTransferIncompressible : public singleton_rule {
    private:
      param<real> energyTransfer ;
    public:

      // Define input and output.
      BoundaryEnergyTransferIncompressible() {
        name_store("energyTransfer(X)",energyTransfer) ;
        output("energyTransfer(X)") ;
        constraint("X,incompressibleFlow") ;
      }

      // Loop over the faces.
      void compute(const sequence &seq) { *energyTransfer=0.0 ; }
  } ;

  register_rule<BoundaryEnergyTransferIncompressible>
    registerBoundaryEnergyTransferIncompressible ;

  // Unit rule to compute the energy transfer through a boundary surface.
  // The parameter specifies the boundary constraint.
  class BoundaryEnergyTransferUnit : public unit_rule {
    private:
      param<real> energyTransfer ;
    public:

      // Define input and output.
      BoundaryEnergyTransferUnit() {
        name_store("energyTransfer(X)",energyTransfer) ;
        output("energyTransfer(X)") ;
        constraint("X,compressibleFlow") ;
      }

      // Initialize the value.
      void compute(const sequence &seq) { *energyTransfer=0.0 ; }
  } ;

  register_rule<BoundaryEnergyTransferUnit> registerBoundaryEnergyTransferUnit ;

  // Apply rule to add the convective to the energy transfer.
  class BoundaryEnergyTransferConvectiveApply : public
  apply_rule<param<real>,Loci::Summation<real> > {
    private:
      const_store<real> rho_f,p_f,h_f,kineticEnergy_f ;
      const_store<real> massFlux ;
      param<real> energyTransfer ;
    public:

      // Define input and output.
      BoundaryEnergyTransferConvectiveApply() {
        name_store("rho_f",rho_f) ;
        name_store("p_f",p_f) ;
        name_store("h_f",h_f) ;
        name_store("kineticEnergy_f",kineticEnergy_f) ;
        name_store("massFlux",massFlux) ;
        name_store("energyTransfer(X)",energyTransfer) ;
        input("rho_f,p_f,h_f,kineticEnergy_f,massFlux") ;
        output("energyTransfer(X)") ;
        constraint("X,compressibleFlow") ;
      }

      // Add the face value to the total.
      void calculate(Entity face) {
        real e_f=h_f[face]-kineticEnergy_f[face]-p_f[face]/rho_f[face] ;
        join(*energyTransfer,massFlux[face]*e_f) ;
      }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryEnergyTransferConvectiveApply>
    registerBoundaryEnergyTransferConvectiveApply ;

  // Apply rule to add the diffusive contributions to the energy transfer.
  class BoundaryEnergyTransferDiffusiveApply : public
  apply_rule<param<real>,Loci::Summation<real> > {
    private:
      const_store<real> qWall ;
      const_store<Area> area ;
      param<real> energyTransfer ;
    public:

      // Define input and output.
      BoundaryEnergyTransferDiffusiveApply() {
        name_store("qWall",qWall) ;
        name_store("area",area) ;
        name_store("energyTransfer(X)",energyTransfer) ;
        input("qWall,area") ;
        output("energyTransfer(X)") ;
        constraint("X,compressibleFlow,viscousFlow") ;
      }

      // Add the face value to the total.
      void calculate(Entity face) {
        join(*energyTransfer,qWall[face]*area[face].sada) ;
      }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryEnergyTransferDiffusiveApply>
    registerBoundaryEnergyTransferDiffusiveApply ;

  // Unit rule to compute the total pressure force on the boundary.
  // The parameter specifies the boundary constraint.
  class BoundaryPressureForceUnit : public unit_rule {
    private:
      param<vect3d> boundaryPressureForce ;
    public:

      // Define input and output.
      BoundaryPressureForceUnit() {
        name_store("boundaryPressureForce(X)",boundaryPressureForce) ;
        output("boundaryPressureForce(X)") ;
        constraint("UNIVERSE") ;
      }

      // Initialize the value.
      void compute(const sequence &seq) {
        *boundaryPressureForce=vect3d(0.0,0.0,0.0) ;
      }
  } ;

  register_rule<BoundaryPressureForceUnit>
    registerBoundaryPressureForceUnit ;

  // Apply rule to compute the pressure force on all no-slip boundaries.
  class BoundaryPressureForceApply : public apply_rule<param<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_store<real> p_f ;
      const_store<Area> area ;
      param<vect3d> boundaryPressureForce ;
    public:

      // Define input and output.
      BoundaryPressureForceApply() {
        name_store("p_f",p_f) ;
        name_store("area",area) ;
        name_store("boundaryPressureForce(X)",boundaryPressureForce) ;
        input("p_f,area") ;
        output("boundaryPressureForce(X)") ;
        constraint("X") ;
      }

      // Add the face value to the total.
      void calculate(Entity face) {
        join(*boundaryPressureForce,p_f[face]*area[face].n*area[face].sada) ;
      }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryPressureForceApply> registerBoundaryPressureForceApply ;

  // Unit rule to compute the viscous force on all no-slip boundaries.
  class BoundaryViscousForceUnit : public unit_rule {
    private:
      param<vect3d> boundaryViscousForce ;
    public:

      // Define input and output.
      BoundaryViscousForceUnit() {
        name_store("boundaryViscousForce(X)",boundaryViscousForce) ;
        output("boundaryViscousForce(X)") ;
        constraint("momentumEquationOptions") ;
      }

      // Initialize the value.
      void compute(const sequence &seq) {
        *boundaryViscousForce=vect3d(0.0,0.0,0.0) ;
      }
  } ;

  register_rule<BoundaryViscousForceUnit> registerBoundaryViscousForceUnit ;

  // Bogus rule for inviscid flow.
  class BoundaryViscousForceInviscidApply : public apply_rule<param<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      param<vect3d> boundaryViscousForce ;
    public:

      // Define input and output.
      BoundaryViscousForceInviscidApply() {
        name_store("boundaryViscousForce(X)",boundaryViscousForce) ;
        output("boundaryViscousForce(X)") ;
        constraint("inviscidFlow") ;
      }

      // Add the face value to the total. Note that this is the force on the
      // wall itself, not the force on the fluid.
      void calculate(Entity face) {
        join(*boundaryViscousForce,vect3d(0.0,0.0,0.0)) ;
      }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryViscousForceInviscidApply>
    registerBoundaryViscousForceInviscidApply ;

  // Apply rule to add the contribution to the viscous force for no-slip
  // boundaries where wall functions are in use.
  class BoundaryViscousForceWallFunctionApply : public apply_rule<param<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_store<vect3d> tauWall ;
      const_store<Area> area ;
      param<vect3d> boundaryViscousForce ;
    public:

      // Define input and output.
      BoundaryViscousForceWallFunctionApply() {
        name_store("tauWall",tauWall) ;
        name_store("area",area) ;
        name_store("boundaryViscousForce(X)",boundaryViscousForce) ;
        input("tauWall,area") ;
        output("boundaryViscousForce(X)") ;
        constraint("X,ref->wallFunction_BCoption") ;
      }

      // Add the face value to the total. Note that this is the force on the
      // wall itself, not the force on the fluid.
      void calculate(Entity face) {
        join(*boundaryViscousForce,tauWall[face]*(-area[face].sada)) ;
      }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryViscousForceWallFunctionApply>
    registerBoundaryViscousForceWallFunctionApply ;

  // Apply rule to add the contribution to the viscous force for no-slip
  // boundaries where wall functions are not in use.
  class BoundaryViscousForceNoWallFunctionApply : public
  apply_rule<param<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map ci ;
      const_store<tens3d> vGradient ;
      const_store<real> laminarViscosity_f ;
      const_store<Area> area ;
      const_store<real> noWallFunction ;
      param<vect3d> boundaryViscousForce ;
    public:

      // Define input and output.
      BoundaryViscousForceNoWallFunctionApply() {
        name_store("ci",ci) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        name_store("area",area) ;
        name_store("noWallFunction",noWallFunction) ;
        name_store("boundaryViscousForce(X)",boundaryViscousForce) ;
        input("ci->gradv3d(v),laminarViscosity_f,area,noWallFunction") ;
        output("boundaryViscousForce(X)") ;
        constraint("X") ;
      }

      // Add the face value to the total. Note that this is the force on the
      // wall itself, not the force on the fluid.
      void calculate(Entity face) {
        join(*boundaryViscousForce,noWallFunction[face]*
          dotTemp(vGradient[ci[face]]+Transpose(vGradient[ci[face]]),
          area[face].n*(-area[face].sada)*laminarViscosity_f[face])) ;
      }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryViscousForceNoWallFunctionApply>
    registerBoundaryViscousForceNoWallFunctionApply ;

  // Write out boundary area, as well as mass, momentum and energy transfer
  // data.
  class BoundaryOutput : public pointwise_rule {
    private:
      string boundaryConstraint ;
      const_param<int> nCycle ;
      const_param<real> sTime ;
      const_param<real> totalArea ;
      const_param<real> massTransfer ;
      const_param<real> energyTransfer ;
      const_param<vect3d> boundaryPressureForce ;
      const_param<vect3d> boundaryViscousForce ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      BoundaryOutput(const char *boundaryConstraint) :
      boundaryConstraint(boundaryConstraint) {
        string totalAreaVariable=string("totalArea(")+boundaryConstraint+
          "){n}" ;
        string massTransferVariable=string("massTransfer(")+boundaryConstraint+
          "){n}" ;
        string energyTransferVariable=string("energyTransfer(")+
          boundaryConstraint+"){n}" ;
        string boundaryPressureForceVariable=string("boundaryPressureForce(")+
          boundaryConstraint+"){n}" ;
        string boundaryViscousForceVariable=string("boundaryViscousForce(")+
          boundaryConstraint+"){n}" ;
        name_store("ncycle{n}",nCycle) ;
        name_store("stime{n}",sTime) ;
        name_store(totalAreaVariable,totalArea) ;
        name_store(massTransferVariable,massTransfer) ;
        name_store(energyTransferVariable,energyTransfer) ;
        name_store(boundaryPressureForceVariable,boundaryPressureForce) ;
        name_store(boundaryViscousForceVariable,boundaryViscousForce) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        input("ncycle{n},stime{n}") ;
        input(totalAreaVariable) ;
        input(massTransferVariable) ;
        input(energyTransferVariable) ;
        input(boundaryPressureForceVariable) ;
        input(boundaryViscousForceVariable) ;
        output("OUTPUT{n}") ;
        conditional("do_print{n}") ;
      }

      // Write out the data.
      void compute(const sequence &seq) {
        if(Loci::MPI_rank!=0) return ;

        // Write data to file in output directory.
        string name(boundaryConstraint) ;
        string filename="output/flux_"+name+".dat" ;
        ofstream out  ;
        if(*nCycle==0) out.open(filename.c_str(),ios::out) ;
        else out.open(filename.c_str(),ios::app) ;
        if(!out.fail()){
          out.precision(10) ;
          out << *nCycle << " " << *sTime << " " << *massTransfer << " "
//          << tot_mom.x << ' ' << tot_mom.y << ' ' << tot_mom.z << ' '
            << 0.0 << " " << 0.0 << " " << 0.0 << " " << *energyTransfer << " "
            << *totalArea << endl ;
        }

        // Write data to standard output.
        cout.setf(ios::scientific,ios::floatfield) ; cout.precision(4) ;
        cout << "Integrated Boundary Data (" << boundaryConstraint << ")"
          << endl ;
        cout << "  total area = " << *totalArea << " m^2" << endl ;
        cout << "  mass transfer = " << *massTransfer << " kg/s" << endl ;
        cout << "  energy transfer = " << *energyTransfer << " W" << endl ;
        cout << "  pressure force = " << *boundaryPressureForce << " N"
          << endl ;
        cout << "  viscous force = " << *boundaryViscousForce << " N"
          << endl ;
      }
  } ;

  // Use a macro to make it easy to create new boudary output classes.
  #define OUTPUT_INT_BC(X) class X : public BoundaryOutput {\
                                  public:\
                                  X() : BoundaryOutput(# X ) {}\
                                  };register_rule<X> register_##X;

  // Define individual classes for each boundary constraint.
  OUTPUT_INT_BC(incompressibleInlet_BC) ;
  OUTPUT_INT_BC(subsonicInlet_BC) ;
  OUTPUT_INT_BC(supersonicInlet_BC) ;
  OUTPUT_INT_BC(totalPressureInlet_BC) ;
  OUTPUT_INT_BC(extrapolatedPressureOutlet_BC) ;
  OUTPUT_INT_BC(fixedPressureOutlet_BC) ;
  OUTPUT_INT_BC(noslip_BC) ;
  OUTPUT_INT_BC(symmetry_BC) ;

  // Unit rule to initialize pressure force. Note that we are now attaching
  // the iteration {n,it} to this variable, so that we can explicitly use
  // area{n} for moving grid problems, which is the old area.
  class PressureForceUnit : public unit_rule {
    private:
      param<vect3d> pressureForce ;
    public:

      // Define input and output.
      PressureForceUnit() {
        name_store("pressureForce{n}",pressureForce) ;
        constraint("UNIVERSE{n}") ;
        output("pressureForce{n}") ;
      }

      // Initialize the value.
      void compute(const sequence &seq) { *pressureForce=vect3d(0.0,0.0,0.0) ; }
  } ;

  register_rule<PressureForceUnit> registerPressureForceUnit ;

  // Apply rule to compute the pressure force on all no-slip boundaries.
  class PressureForceApply : public apply_rule<param<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_store<real> p_f ;
      const_store<Area> area ;
      param<vect3d> pressureForce ;
    public:

      // Define input and output.
      PressureForceApply() {
        name_store("p_f{n}",p_f) ;
        name_store("area{n}",area) ;
        name_store("pressureForce{n}",pressureForce) ;
        input("p_f{n},area{n}") ;
        output("pressureForce{n}") ;
        constraint("noslip_BC") ;
      }

      // Add the face value to the total.
      void calculate(Entity face) {
        join(*pressureForce,p_f[face]*area[face].n*area[face].sada) ;
      }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PressureForceApply> registerPressureForceApply ;

  // Apply rule to compute the pressure force on all slip boundaries.
  class PressureForceSlipApply : public apply_rule<param<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_store<real> p_f ;
      const_store<Area> area ;
      param<vect3d> pressureForce ;
    public:

      // Define input and output.
      PressureForceSlipApply() {
        name_store("p_f{n}",p_f) ;
        name_store("area{n}",area) ;
        name_store("pressureForce{n}",pressureForce) ;
        input("p_f{n},area{n}") ;
        output("pressureForce{n}") ;
        constraint("slip_BC") ;
      }

      // Add the face value to the total.
      void calculate(Entity face) {
        join(*pressureForce,p_f[face]*area[face].n*area[face].sada) ;
      }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PressureForceSlipApply> registerPressureForceSlipApply ;

  // Unit rule to compute the viscous force on all no-slip boundaries.
  // Note that we are now attaching the iteration {n,it} to this variable,
  // so that we can explicitly use area{n} and tauWall{n}.
  class ViscousForceUnit : public unit_rule {
    private:
      param<vect3d> viscousForce ;
    public:

      // Define input and output.
      ViscousForceUnit() {
        name_store("viscousForce{n}",viscousForce) ;
        constraint("noslip_BC{n}") ;
        output("viscousForce{n}") ;
      }

      // Initialize the value.
      void compute(const sequence &seq) { *viscousForce=vect3d(0.0,0.0,0.0) ; }
  } ;

  register_rule<ViscousForceUnit> registerViscousForceUnit ;

  // Apply rule to add the contribution to the viscous force for no-slip
  // boundaries where wall functions are in use.
  class ViscousForceWallFunctionApply : public apply_rule<param<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_store<vect3d> tauWall ;
      const_store<Area> area ;
      param<vect3d> viscousForce ;
    public:

      // Define input and output.
      ViscousForceWallFunctionApply() {
        name_store("tauWall{n}",tauWall) ;
        name_store("area{n}",area) ;
        name_store("viscousForce{n}",viscousForce) ;
        input("tauWall{n},area{n}") ;
        output("viscousForce{n}") ;
        constraint("noslip_BC{n},ref->wallFunction_BCoption{n}") ;
      }

      // Add the face value to the total. Note that this is the force on the
      // wall itself, not the force on the fluid.
      void calculate(Entity face) {
        join(*viscousForce,tauWall[face]*(-area[face].sada)) ;
      }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ViscousForceWallFunctionApply>
    registerViscousForceWallFunctionApply ;

  // Apply rule to add the contribution to the viscous force for no-slip
  // boundaries where wall functions are not in use.
  class ViscousForceNoWallFunctionApply : public apply_rule<param<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_Map ci ;
      const_store<tens3d> vGradient ;
      const_store<real> laminarViscosity_f ;
      const_store<Area> area ;
      const_store<real> noWallFunction ;
      param<vect3d> viscousForce ;
    public:

      // Define input and output.
      ViscousForceNoWallFunctionApply() {
        name_store("ci",ci) ;
        name_store("gradv3d(v){n}",vGradient) ;
        name_store("laminarViscosity_f{n}",laminarViscosity_f) ;
        name_store("area{n}",area) ;
        name_store("noWallFunction",noWallFunction) ;
        name_store("viscousForce{n}",viscousForce) ;
        input("ci->gradv3d(v){n},laminarViscosity_f{n}") ;
        input("area{n},noWallFunction") ;
        output("viscousForce{n}") ;
        constraint("noslip_BC") ;
      }

      // Add the face value to the total. Note that this is the force on the
      // wall itself, not the force on the fluid.
      void calculate(Entity face) {
        join(*viscousForce,noWallFunction[face]*dotTemp(vGradient[ci[face]]+
          Transpose(vGradient[ci[face]]),area[face].n*(-area[face].sada)*
          laminarViscosity_f[face])) ;
      }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ViscousForceNoWallFunctionApply>
    registerViscousForceNoWallFunctionApply ;

  // Unit rule to initialize pressure moment.
  class PressureMomentUnit : public unit_rule {
    private:
      param<vect3d> pressureMoment ;
    public:

      // Define input and output.
      PressureMomentUnit() {
        name_store("pressureMoment",pressureMoment) ;
        output("pressureMoment") ;
//      constraint("UNIVERSE") ;
        constraint("ref->momentCenter_BC") ;
      }

      // Initialize the value.
      void compute(const sequence &seq) { *pressureMoment=vect3d(0.0,0.0,0.0) ; }
  } ;

  register_rule<PressureMomentUnit> registerPressureMomentUnit ;

  // Apply rule to compute the pressure moment on all boundaries that
  // have momentCenter specified. Since momentCenter is allowed only
  // for slip and noslip boundaries, this rule will never execute for
  // faces on other boundary types.
  class PressureMomentApply : public apply_rule<param<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_Map ref ;
      const_store<real> p_f ;
      const_store<Area> area ;
      const_store<vect3d> faceCenter,momentCenter_BC ;
      param<vect3d> pressureMoment ;
    public:

      // Define input and output.
      PressureMomentApply() {
        name_store("ref",ref) ;
        name_store("p_f",p_f) ;
        name_store("area",area) ;
        name_store("facecenter",faceCenter) ;
        name_store("momentCenter_BC",momentCenter_BC) ;
        name_store("pressureMoment",pressureMoment) ;
        input("p_f,area,facecenter,ref->momentCenter_BC") ;
        output("pressureMoment") ;
      }

      // Add the face value to the total.
      void calculate(Entity face) {
        join(*pressureMoment,cross(faceCenter[face]-momentCenter_BC[ref[face]],
          p_f[face]*area[face].n*area[face].sada)) ;
      }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PressureMomentApply> registerPressureMomentApply ;

  // Unit rule to initialize the viscous moment.
  class ViscousMomentUnit : public unit_rule {
    private:
      param<vect3d> viscousMoment ;
    public:

      // Define input and output.
      ViscousMomentUnit() {
        name_store("viscousMoment",viscousMoment) ;
        output("viscousMoment") ;
//      constraint("UNIVERSE") ;
        constraint("ref->momentCenter_BC") ;
      }

      // Initialize the value.
      void compute(const sequence &seq) { *viscousMoment=vect3d(0.0,0.0,0.0) ; }
  } ;

  register_rule<ViscousMomentUnit> registerViscousMomentUnit ;

  // Apply rule to add the contribution to the viscous moment for no-slip
  // boundaries where wall functions are in use and momentCenter has been
  // specified.
  class ViscousMomentWallFunctionApply : public apply_rule<param<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_Map ref ;
      const_store<vect3d> tauWall ;
      const_store<Area> area ;
      const_store<vect3d> faceCenter,momentCenter_BC ;
      param<vect3d> viscousMoment ;
    public:

      // Define input and output.
      ViscousMomentWallFunctionApply() {
        name_store("ref",ref) ;
        name_store("tauWall",tauWall) ;
        name_store("area",area) ;
        name_store("facecenter",faceCenter) ;
        name_store("momentCenter_BC",momentCenter_BC) ;
        name_store("viscousMoment",viscousMoment) ;
        input("tauWall,area,facecenter,ref->momentCenter_BC") ;
        output("viscousMoment") ;
        constraint("noslip_BC,ref->(wallFunction_BCoption,momentCenter_BC)") ;
      }

      // Add the face value to the total. Note that this is the moment of
      // the force on the wall itself, not the force on the fluid.
      void calculate(Entity face) {
        join(*viscousMoment,cross(faceCenter[face]-momentCenter_BC[ref[face]],
          tauWall[face]*(-area[face].sada))) ;
      }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ViscousMomentWallFunctionApply>
    registerViscousMomentWallFunctionApply ;

  // Apply rule to add the contribution to the viscous moment for no-slip
  // boundaries where wall functions are not in use and momentCenter has
  // been specified.
  class ViscousMomentNoWallFunctionApply : public apply_rule<param<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_Map ci,ref ;
      const_store<tens3d> vGradient ;
      const_store<real> laminarViscosity_f ;
      const_store<Area> area ;
      const_store<real> noWallFunction ;
      const_store<vect3d> faceCenter,momentCenter_BC ;
      param<vect3d> viscousMoment ;
    public:

      // Define input and output.
      ViscousMomentNoWallFunctionApply() {
        name_store("ci",ci) ;
        name_store("ref",ref) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        name_store("area",area) ;
        name_store("noWallFunction",noWallFunction) ;
        name_store("facecenter",faceCenter) ;
        name_store("momentCenter_BC",momentCenter_BC) ;
        name_store("viscousMoment",viscousMoment) ;
        input("ci->gradv3d(v),laminarViscosity_f,area,noWallFunction") ;
        input("facecenter,ref->momentCenter_BC") ;
        output("viscousMoment") ;
        constraint("noslip_BC,ref->momentCenter_BC") ;
      }

      // Add the face value to the total. Note that this is the moment of the
      // force on the wall itself, not the force on the fluid.
      void calculate(Entity face) {
        join(*viscousMoment,cross(faceCenter[face]-momentCenter_BC[ref[face]],
          noWallFunction[face]*dotTemp(vGradient[ci[face]]+
          Transpose(vGradient[ci[face]]),area[face].n*(-area[face].sada)*
          laminarViscosity_f[face]))) ;
      }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

//register_rule<ViscousMomentNoWallFunctionApply>
//  registerViscousMomentNoWallFunctionApply ;

  // TEMPORARY
  class ViscousMomentNoWallFunctionApplyTemp : public apply_rule<param<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_Map ci,ref ;
      const_store<vect3d> v,v_f;
      const_store<vect3d> cellCenter ;
      const_store<real> laminarViscosity_f ;
      const_store<Area> area ;
      const_store<real> noWallFunction ;
      const_store<vect3d> faceCenter,momentCenter_BC ;
      param<vect3d> viscousMoment ;
    public:

      // Define input and output.
      ViscousMomentNoWallFunctionApplyTemp() {
        name_store("ci",ci) ;
        name_store("ref",ref) ;
        name_store("v",v) ;
        name_store("v_f",v_f) ;
        name_store("cellcenter",cellCenter) ;
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        name_store("area",area) ;
        name_store("noWallFunction",noWallFunction) ;
        name_store("facecenter",faceCenter) ;
        name_store("momentCenter_BC",momentCenter_BC) ;
        name_store("viscousMoment",viscousMoment) ;
        input("ci->(v,cellcenter),v_f") ;
        input("laminarViscosity_f,area,noWallFunction") ;
        input("facecenter,ref->momentCenter_BC") ;
        output("viscousMoment") ;
        constraint("noslip_BC,ref->momentCenter_BC") ;
      }

      // Add the face value to the total. Note that this is the moment of the
      // force on the wall itself, not the force on the fluid.
      void calculate(Entity face) {
        vect3d dV=v[ci[face]]-v_f[face] ;
        vect3d dVTangential=dV-dot(dV,area[face].n)*area[face].n ;
        vect3d t=dVTangential/norm(dVTangential) ;
        vect3d r(0.0,faceCenter[face].y,faceCenter[face].z) ;
        real dX=norm(faceCenter[face]-cellCenter[ci[face]]) ;
        vect3d moment=cross(r,laminarViscosity_f[face]*norm(dVTangential)/dX*
          area[face].sada*t) ;
        join(*viscousMoment,moment) ;
      }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ViscousMomentNoWallFunctionApplyTemp>
    registerViscousMomentNoWallFunctionApplyTemp ;

  // Write out forces and energy transfer for all noslip walls.
  class NoSlipOutput : public pointwise_rule {
    private:
      const_param<vect3d> pressureForce,viscousForce ;
      const_param<real> energyTransfer ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      NoSlipOutput() {
        name_store("pressureForce{n}",pressureForce) ;
        name_store("viscousForce{n}",viscousForce) ;
        name_store("energyTransfer(noslip_BC){n}",energyTransfer) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        input("pressureForce{n},viscousForce{n}") ;
        input("energyTransfer(noslip_BC){n}") ;
        output("OUTPUT{n}") ;
        conditional("do_print{n}") ;
      }

      // Write out the data.
      void compute(const sequence &seq) {
        if(Loci::MPI_rank!=0) return ;
        cout.setf(ios::scientific,ios::floatfield) ; cout.precision(4) ;
        cout << "Force and Energy Transfer Data (noslip_BC)" << endl ;
        cout << "  total force = " << *pressureForce+*viscousForce << " N"
          << endl ;
        cout << "    pressure force = " << *pressureForce << " N" << endl ;
        cout << "    viscous force = " << *viscousForce << " N" << endl ;
        cout << "  energy transfer = " << *energyTransfer << " W" << endl ;
      }
  } ;

  register_rule<NoSlipOutput> registerNoSlipOutput ;

  // Write out moment of forces for noslip walls that have the "momentCenter"
  // spedified.
  class NoSlipMomentOutput : public pointwise_rule {
    private:
      const_param<vect3d> pressureMoment,viscousMoment ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      NoSlipMomentOutput() {
        name_store("pressureMoment{n}",pressureMoment) ;
        name_store("viscousMoment{n}",viscousMoment) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        input("pressureMoment{n},viscousMoment{n}") ;
        output("OUTPUT{n}") ;
        conditional("do_print{n}") ;
      }

      // Write out the data.
      void compute(const sequence &seq) {
        if(Loci::MPI_rank!=0) return ;
        cout.setf(ios::scientific,ios::floatfield) ; cout.precision(4) ;
        cout << "Moment Data (noslip_BC with momentCenter specified)." << endl ;
        cout << "  total moment = " << *pressureMoment+*viscousMoment << " N-m"
          << endl ;
        cout << "    pressure moment = " << *pressureMoment << " N-m" << endl ;
        cout << "    viscous moment = " << *viscousMoment << " N-m" << endl ;
      }
  } ;

  register_rule<NoSlipMomentOutput> registerNoSlipMomentOutput ;

  // Write out pressure force for all slip walls.
  class SlipOutput : public pointwise_rule {
    private:
      const_param<vect3d> pressureForce ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      SlipOutput() {
        name_store("pressureForce{n}",pressureForce) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        input("pressureForce{n}") ;
        output("OUTPUT{n}") ;
        constraint("slip_BC") ;
        conditional("do_print{n}") ;
      }

      // Write out the data.
      void compute(const sequence &seq) {
        if(Loci::MPI_rank!=0) return ;
        cout.setf(ios::scientific,ios::floatfield) ; cout.precision(4) ;
        cout << "Force Data (slip_BC)" << endl ;
        cout << "    pressure force = " << *pressureForce << " N" << endl ;
      }
  } ;

  register_rule<SlipOutput> registerSlipOutput ;
}








