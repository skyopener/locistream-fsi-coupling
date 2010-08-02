//-----------------------------------------------------------------------------
// Description: This file contains rules for wall functions.
//
// Authors: Siddharth Thakur and Jeff Wright
//          Streamline Numerics Inc.
//          Gainesville, Florida 32609
//-----------------------------------------------------------------------------
                                                                                
// Loci includes.
#include <Loci.h>
using Loci::Area ;
                                                                                
// StreamUns includes.
#include "root.h"
#include "sciTypes.h"
                                                                                
// Rules for the Nichols and Nelson wall function. This code was written and
// supplied by Ed Luke. Some recoding has taken place for descriptive purposes
// and to make the rules compatible with our variables.
namespace streamUns {

  // Rule to create the constraint which identifies wall-function cells.
  class WallFunctionCells : public pointwise_rule {
    private:
      const_Map ci ;
      store<bool> wallFunctionCells ;
    public:

      // Define input and output.
      WallFunctionCells() {
        name_store("ci",ci) ;
        name_store("wallFunctionCells",wallFunctionCells) ;
        input("ci") ;
        output("ci->wallFunctionCells") ;
        constraint("ref->wallFunction_BCoption") ;
      }

      // Empty compute method.
      virtual void compute (const sequence &seq) {}
  } ;

  register_rule<WallFunctionCells> registerWallFunctionCells ;


  // Priority rule to set the wall temperature for adiabatic walls when wall
  // functions are in use. The Croco-Busemann equation is used with the near-
  // wall temperature and velocity.
  class BoundaryTemperatureAdiabaticWallFunction : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<vect3d> v ;
      const_store<real> temperature ;
      const_store<real> laminarViscosity ;
      const_store<real> cp ;
      const_store<real> thermalConductivity ;
      const_store<Area> area ;
      store<real> temperature_f ;
    public:

      // Define input and output.
      BoundaryTemperatureAdiabaticWallFunction() {
        name_store("ci",ci) ;
        name_store("v",v) ;
        name_store("temperature",temperature) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("cp",cp) ;
        name_store("kconduct(temperature,p,y)",thermalConductivity) ;
        name_store("area",area) ;
        name_store("wallFunction::temperature_f",temperature_f) ;
        input("ci->(v,temperature,laminarViscosity,cp)") ;
        input("ci->kconduct(temperature,p,y),area") ;
        output("wallFunction::temperature_f") ;
        constraint("noslip_BC,ref->wallFunction_BCoption") ;
        constraint("ref->adiabatic_BCoption") ;
      }

      // Calculate temperature for a single face.
      void calculate(Entity face) {

        // Decompose the cell velocity.
        const real vNormal=dot(v[ci[face]],area[face].n) ;
        const vect3d vTangential=v[ci[face]]-vNormal*area[face].n ;

        // Compute the recovery factor.
        const real recoveryFactor=pow(laminarViscosity[ci[face]]*
          cp[ci[face]]/thermalConductivity[ci[face]],1.0/3.0) ;

        // Compute the temperature at the wall. Note that we are not using
        // the form of the Crocco-Busemann relation in Nichols and Nelson's
        // paper, as it seems strange. The following is what Ed is using.
        // Should check this out since this form apparently violates energy
        // conservation.
        temperature_f[face]=temperature[ci[face]]+0.5*recoveryFactor*
          dot(vTangential,vTangential)/cp[ci[face]] ;
      }

      // Calculate temperature for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryTemperatureAdiabaticWallFunction>
    registerBoundaryTemperatureAdiabaticWallFunction ;

  // Class to compute the implicit functional for the adiabatic, incompressible
  // case.
  class YPlusImplicitAdiabaticIncompressible {
    private:
      real rhoWall,mu,u,y,kappa,expKB ;
    public:
      YPlusImplicitAdiabaticIncompressible(real rhoWall,real mu,real u,real y,
      real kappa,real b) : rhoWall(rhoWall),mu(mu),u(u),y(y),kappa(kappa) {
        expKB=exp(-kappa*b) ; }
    public:
      real operator()(real yPlus) {
        const real uPlus=u*rhoWall*y/(yPlus*mu) ; if(uPlus>100) return -1e30 ;
        return yPlus-uPlus-expKB*(exp(kappa*uPlus)-1.0-kappa*uPlus-pow(kappa*
          uPlus,2)/2.0-pow(kappa*uPlus,3)/6.0) ;
      }
  } ;

  // Class to compute the implicit functional for the adiabatic case.
  class YPlusImplicitAdiabaticCompressible {
    private:
      real cp,tWall,rhoWall,mu,u,y,r,kappa,expKB ;
    public:
      YPlusImplicitAdiabaticCompressible(real cp,real tWall,real rhoWall,real
        mu,real u,real y,real r,real kappa,real b) : cp(cp),tWall(tWall),
        rhoWall(rhoWall),mu(mu),u(u),y(y),r(r),kappa(kappa) {
        expKB=exp(-kappa*b) ; }
    public:
      real operator()(real yPlus) {
        const real uPlus=u*rhoWall*y/(yPlus*mu) ; if(uPlus>100) return -1e30 ;
        const real uTau=u/(uPlus+1e-30) ;
        const real G=r*uTau*uTau/(2.0*cp*tWall) ;
        const real sqrtG=sqrt(G)+1e-30 ;
        const real yPlusWhite=exp((kappa/sqrtG)*asin(sqrtG*uPlus))*expKB ;
        return yPlus-(uPlus+yPlusWhite-expKB*(1.0+kappa*uPlus+pow(kappa*uPlus,
          2)/2.0+pow(kappa*uPlus,3)/6.0)) ;
      }
  } ;

  // Class to compute the implicit functional for the non-adiabatic case.
  class YPlusImplicitNonAdiabaticCompressible {
    private:
      real cp,tWall,tCell,kWall,rhoWall,mu,u,y,r,kappa,expKB ;
    public:
      YPlusImplicitNonAdiabaticCompressible(real cp,real tWall,real tCell,real
        kWall,real rhoWall,real mu,real u,real y,real r,real kappa,real b) :
        cp(cp),tWall(tWall),tCell(tCell),kWall(kWall),rhoWall(rhoWall),mu(mu),
        u(u),y(y),r(r),kappa(kappa) { expKB=exp(-kappa*b) ; }
    public:
      real operator()(real yPlus) {
        const real uPlus=u*rhoWall*y/(yPlus*mu) ; if(uPlus>100) return -1e30 ;
        const real uTau=u/(uPlus+1e-30) ;
        const real uTau2=uTau*uTau ;
        const real g=r*uTau2/(2.0*cp*tWall) ;
        const real qW=uTau2*(tCell/tWall-1.0+r*u*u/(2.0*cp*tWall))*
          (rhoWall*tWall*kWall)/(mu*(u+1e-30)) ;
        const real beta=qW*mu/(rhoWall*tWall*kWall*uTau+1e-30) ;
        const real q=sqrt(beta*beta+4.0*g)+1e-30 ;
        const real phi=asin(-beta/q) ;
        const real sqrtG=sqrt(g)+1e-30 ;
        const real yPlusWhite=exp((kappa/sqrtG)*(asin((2.0*g*uPlus-beta)/q)-
          phi))*expKB ;
        return yPlus-(uPlus+yPlusWhite-expKB*(1.0+kappa*uPlus+pow(kappa*uPlus,
          2)/2.0+pow(kappa*uPlus,3)/6.0)) ;
      }
  } ;

  // Wall function rule for adiabatic incompressible flow.
  class WallLawAdiabaticIncompressible : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<vect3d> cellCenter ;
      const_store<real> rho ;
      const_store<vect3d> v ;
      const_store<real> laminarViscosity ;
      const_store<vect3d> faceCenter ;
      const_store<Area> area ;
      const_store<real> rho_f,laminarViscosity_f ;
      store<vect3d> tauWall ;
      store<real> kWall,omegaWall,eddyViscosityWall,yPlusWall ;
    public:
                                                                                
      // Define input and output.
      WallLawAdiabaticIncompressible() {
        name_store("ci",ci) ;
        name_store("cellcenter",cellCenter) ;
        name_store("rho",rho) ;
        name_store("v",v) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("facecenter",faceCenter) ;
        name_store("area",area) ;
        name_store("rho_f",rho_f) ;
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        name_store("tauWallTemp",tauWall) ;
        name_store("kWall",kWall) ;
        name_store("omegaWall",omegaWall) ;
        name_store("eddyViscosityWall",eddyViscosityWall) ;
        name_store("yPlusWallTemp",yPlusWall) ;
        input("ci->(cellcenter,rho,v),ci->laminarViscosity") ;
        input("facecenter,area,rho_f,laminarViscosity_f") ;
        output("tauWallTemp,kWall,omegaWall,eddyViscosityWall,yPlusWallTemp") ;
        constraint("incompressibleFlow,noslip_BC,ref->wallFunction_BCoption") ;
      }

      // Calculate temperature for a single face.
      void calculate(Entity face) {

        // Decompose the cell velocity and compute the magnitude of the
        // tangential component.
        const real vNormal= dot(v[ci[face]],area[face].n) ;
        const vect3d vTangential=v[ci[face]]-vNormal*area[face].n ;
        const real vTangentialMagnitude=sqrt(dot(vTangential,vTangential)) ;
                                                                                
        // Compute the normal distance from the wall face to the adjacent
        // cell center.
        const real yNormal=dot(faceCenter[face]-cellCenter[ci[face]],
          area[face].n) ;

        // Find yPlus using the wall function. If we are out of the range
        // 0.01<yPlus<5000, use the low Reynolds number model values.
        const real kappa=0.41,B=5.5,expKB=exp(-kappa*B) ;
        YPlusImplicitAdiabaticIncompressible func(rho_f[face],
          laminarViscosity_f[face],vTangentialMagnitude,yNormal,kappa,B) ;
        if(func(0.01)>0.0){
          tauWall[face]=vTangential*(-laminarViscosity[ci[face]]/yNormal) ;
          kWall[face]=0 ; omegaWall[face]=60.0*laminarViscosity[ci[face]]/
            (rho[ci[face]]*0.075*pow(yNormal,2)) ;
          eddyViscosityWall[face]=0.0 ; yPlusWall[face]=0.01 ; return ;
        }
        if(func(5000.0)<0.0){
          tauWall[face]=vTangential*(-laminarViscosity[ci[face]]/yNormal) ;
          kWall[face]=0 ; omegaWall[face]=60.0*laminarViscosity[ci[face]]/
            (rho[ci[face]]*0.075*pow(yNormal,2)) ;
          eddyViscosityWall[face]=0.0 ; yPlusWall[face]=5000.0 ; return ;
        }
        yPlusWall[face]=find_root(func,0.01,5000.0,1.0e-06) ;

        // Compute the wall shear stress.
        const real uPlus=vTangentialMagnitude*rho_f[face]*yNormal/
          (yPlusWall[face]*laminarViscosity_f[face]) ;
        const real frictionVelocity=vTangentialMagnitude/uPlus ;
        const real tauWallMagnitude=rho_f[face]*pow(frictionVelocity,2) ;
        tauWall[face]=-tauWallMagnitude*vTangential/vTangentialMagnitude ;

        // Compute the eddy viscosity for the cell next to the wall.
        eddyViscosityWall[face]=laminarViscosity_f[face]*kappa*expKB*
          (exp(kappa*uPlus)-1-kappa*uPlus-pow(kappa*uPlus,2)/2.0) ;

        // Compute k and omega for the cell next to the wall.
        const real cMu=0.09,omegaI=6.0*laminarViscosity_f[face]/(0.075*
          pow(yNormal,2)*rho_f[face]),omegaO=frictionVelocity/(sqrt(cMu)*
          kappa*yNormal) ;
          omegaWall[face]=sqrt(pow(omegaI,2)+pow(omegaO,2)) ;
          kWall[face]=omegaWall[face]*eddyViscosityWall[face]/rho[ci[face]] ;
      }
                                                                                
      // Calculate temperature for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;
                                                                                
  register_rule<WallLawAdiabaticIncompressible>
    registerWallLawAdiabaticIncompressible ;

  // Wall function rule for adiabatic walls.
  class WallLawAdiabaticCompressible : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<vect3d> cellCenter ;
      const_store<real> rho ;
      const_store<vect3d> v ;
      const_store<real> temperature ;
      const_store<real> laminarViscosity ;
      const_store<real> cp ;
      const_store<real> thermalConductivity ;
      const_store<vect3d> faceCenter ;
      const_store<Area> area ;
      const_store<real> rho_f,temperature_f,laminarViscosity_f ;
      store<vect3d> tauWall ;
      store<real> kWall,omegaWall,eddyViscosityWall,yPlusWall ;
    public:
                                                                                
      // Define input and output.
      WallLawAdiabaticCompressible() {
        name_store("ci",ci) ;
        name_store("cellcenter",cellCenter) ;
        name_store("rho",rho) ;
        name_store("v",v) ;
        name_store("temperature",temperature) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("cp",cp) ;
        name_store("kconduct(temperature,p,y)",thermalConductivity) ;
        name_store("facecenter",faceCenter) ;
        name_store("area",area) ;
        name_store("rho_f",rho_f) ;
        name_store("temperature_f",temperature_f) ;
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        name_store("tauWallTemp",tauWall) ;
        name_store("kWall",kWall) ;
        name_store("omegaWall",omegaWall) ;
        name_store("eddyViscosityWall",eddyViscosityWall) ;
        name_store("yPlusWallTemp",yPlusWall) ;
        input("ci->(cellcenter,rho,v,temperature)") ;
        input("ci->(laminarViscosity,cp)") ;
        input("ci->kconduct(temperature,p,y)") ;
        input("facecenter,area,rho_f,temperature_f,laminarViscosity_f");
        output("tauWallTemp,kWall,omegaWall,eddyViscosityWall,yPlusWallTemp") ;
        constraint("compressibleFlow,noslip_BC") ;
        constraint("ref->adiabatic_BCoption,ref->wallFunction_BCoption") ;
      }

      // Calculate temperature for a single face.
      void calculate(Entity face) {

        // Decompose the cell velocity and compute the magnitude of the
        // tangential component.
        const real vNormal= dot(v[ci[face]],area[face].n) ;
        const vect3d vTangential=v[ci[face]]-vNormal*area[face].n ;
        const real vTangentialMagnitude=sqrt(dot(vTangential,vTangential)) ;
                                                                                
        // Compute the recovery factor.
        const real recoveryFactor=pow(laminarViscosity[ci[face]]*
          cp[ci[face]]/thermalConductivity[ci[face]],1.0/3.0) ;

        // Compute the normal distance from the wall face to the adjacent
        // cell center.
        const real yNormal=dot(faceCenter[face]-cellCenter[ci[face]],
          area[face].n) ;

        // Find yPlus using the wall function. If we are out of the range
        // 0.01<yPlus<5000, use the low Reynolds number model values.
        const real kappa=0.41,B=5.5,expKB=exp(-kappa*B) ;
        YPlusImplicitAdiabaticCompressible func(cp[ci[face]],
          temperature_f[face],rho_f[face],laminarViscosity_f[face],
          vTangentialMagnitude,yNormal,recoveryFactor,kappa,B) ;
        if(func(0.01)>0.0){
          tauWall[face]=vTangential*(-laminarViscosity[ci[face]]/yNormal) ;
          kWall[face]=0 ; omegaWall[face]=60.0*laminarViscosity[ci[face]]/
            (rho[ci[face]]*0.075*pow(yNormal,2)) ;
          eddyViscosityWall[face]=0.0 ; yPlusWall[face]=0.01 ;
          return ;
        }
        if(func(5000.0)<0.0){
          tauWall[face]=vTangential*(-laminarViscosity[ci[face]]/yNormal) ;
          kWall[face]=0 ; omegaWall[face]=60.0*laminarViscosity[ci[face]]/
            (rho[ci[face]]*0.075*pow(yNormal,2)) ;
          eddyViscosityWall[face]=0.0 ; yPlusWall[face]=5000.0 ;
          return ;
        }
        yPlusWall[face]=find_root(func,0.01,5000.0,1.0e-06) ;

        // Compute the wall shear stress.
        const real uPlus=vTangentialMagnitude*rho_f[face]*yNormal/
          (yPlusWall[face]*laminarViscosity_f[face]) ;
        const real frictionVelocity=vTangentialMagnitude/uPlus ;
        const real tauWallMagnitude=rho_f[face]*pow(frictionVelocity,2) ;
        tauWall[face]=-tauWallMagnitude*vTangential/vTangentialMagnitude ;

        // Compute the eddy viscosity for the cell next to the wall.
        eddyViscosityWall[face]=laminarViscosity_f[face]*kappa*expKB*(-1-
          kappa*uPlus-pow(kappa*uPlus,2)/2.0) ;
        const real gamma=recoveryFactor*pow(frictionVelocity,2)/(2.0*
          cp[ci[face]]*temperature_f[face]) ;
        const real yWhite=exp(kappa/sqrt(gamma)*asin(sqrt(gamma)*uPlus))*expKB ;
        const real yWhiteDerivative=yWhite*kappa*sqrt(1.0-gamma*pow(uPlus,2)) ;
        eddyViscosityWall[face]+=(1.0+yWhiteDerivative-laminarViscosity
          [ci[face]]/laminarViscosity_f[face])*laminarViscosity_f[face] ;

        // Compute k and omega for the cell next to the wall.
        const real cMu=0.09,omegaI=6.0*laminarViscosity_f[face]/(0.075*
          pow(yNormal,2)*rho_f[face]),omegaO=frictionVelocity/(sqrt(cMu)*
          kappa*yNormal) ;
          omegaWall[face]=sqrt(pow(omegaI,2)+pow(omegaO,2)) ;
          kWall[face]=omegaWall[face]*eddyViscosityWall[face]/rho[ci[face]] ;
      }
                                                                                
      // Calculate temperature for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;
                                                                                
  register_rule<WallLawAdiabaticCompressible>
    registerWallLawAdiabaticCompressible ;

  // Wall function rule for walls with specified temperature.
  class WallLawSpecifiedTemperatureCompressible : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<vect3d> cellCenter ;
      const_store<real> rho ;
      const_store<vect3d> v ;
      const_store<real> temperature ;
      const_store<real> laminarViscosity ;
      const_store<real> cp ;
      const_store<real> thermalConductivity ;
      const_store<vect3d> faceCenter ;
      const_store<Area> area ;
      const_store<real> rho_f,temperature_f,laminarViscosity_f ;
      store<vect3d> tauWall ;
      store<real> kWall,omegaWall,eddyViscosityWall,yPlusWall,qWall ;
    public:
                                                                                
      // Define input and output.
      WallLawSpecifiedTemperatureCompressible() {
        name_store("ci",ci) ;
        name_store("cellcenter",cellCenter) ;
        name_store("rho",rho) ;
        name_store("v",v) ;
        name_store("temperature",temperature) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("cp",cp) ;
        name_store("kconduct(temperature,p,y)",thermalConductivity) ;
        name_store("facecenter",faceCenter) ;
        name_store("area",area) ;
        name_store("rho_f",rho_f) ;
        name_store("temperature_f",temperature_f) ;
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        name_store("tauWallTemp",tauWall) ;
        name_store("kWall",kWall) ;
        name_store("omegaWall",omegaWall) ;
        name_store("eddyViscosityWall",eddyViscosityWall) ;
        name_store("yPlusWallTemp",yPlusWall) ;
        name_store("qWallTemp",qWall) ;
        input("ci->(cellcenter,rho,v,temperature)") ;
        input("ci->(laminarViscosity,cp)") ;
        input("ci->kconduct(temperature,p,y)") ;
        input("facecenter,area,rho_f,temperature_f,laminarViscosity_f");
        output("tauWallTemp,kWall,omegaWall,eddyViscosityWall,yPlusWallTemp") ;
        output("qWallTemp") ;
        constraint("compressibleFlow,noslip_BC,ref->wallFunction_BCoption") ;
        constraint("specifiedTemperature_BC") ;
      }

      // Calculate temperature for a single face.
      void calculate(Entity face) {
                                                                                
        // Decompose the cell velocity and compute the magnitude of the
        // tangential component.
        const real vNormal=dot(v[ci[face]],area[face].n) ;
        const vect3d vTangential=v[ci[face]]-vNormal*area[face].n ;
        const real vTangentialMagnitude=sqrt(dot(vTangential,vTangential)) ;
                                                                                
        // Compute the recovery factor.
        const real recoveryFactor=pow(laminarViscosity[ci[face]]*
          cp[ci[face]]/thermalConductivity[ci[face]],1.0/3.0) ;
                                                                                
        // Compute the normal distance from the wall face to the adjacent
        // cell center.
        const real yNormal=dot(faceCenter[face]-cellCenter[ci[face]],
          area[face].n) ;

        // Find yPlus using the wall function. If we are out of the range
        // 0.01<yPlus<5000, use the low Reynolds number model values.
        const real kappa=0.41,B=5.5,expKB=exp(-kappa*B) ;
        YPlusImplicitNonAdiabaticCompressible func(cp[ci[face]],
          temperature_f[face],temperature[ci[face]],
          thermalConductivity[ci[face]],rho_f[face],laminarViscosity_f[face],
          vTangentialMagnitude,yNormal,recoveryFactor,kappa,B) ;
        if(func(0.01)>0.0){
          tauWall[face]=vTangential*(-laminarViscosity[ci[face]]/yNormal) ;
          qWall[face]=thermalConductivity[ci[face]]*(temperature[ci[face]]-
            temperature_f[face])/yNormal ;
          kWall[face]=0 ; omegaWall[face]=60.0*laminarViscosity[ci[face]]/
            (rho[ci[face]]*0.075*pow(yNormal,2)) ;
          eddyViscosityWall[face]=0.0 ; yPlusWall[face]=0.01 ;
          return ;
        }
        if(func(5000.0)<0.0){
          tauWall[face]=vTangential*(-laminarViscosity[ci[face]]/yNormal) ;
          qWall[face]=thermalConductivity[ci[face]]*(temperature[ci[face]]-
            temperature_f[face])/yNormal ;
          kWall[face]=0 ; omegaWall[face]=60.0*laminarViscosity[ci[face]]/
            (rho[ci[face]]*0.075*pow(yNormal,2.0)) ;
          eddyViscosityWall[face]=0.0 ; yPlusWall[face]=5000.0 ;
          return ;
        }
        yPlusWall[face]=find_root(func,0.01,5000.0,1.0e-08) ;

        // Compute the wall shear stress.
        real uPlus=vTangentialMagnitude*rho_f[face]*yNormal/
          (yPlusWall[face]*laminarViscosity_f[face]) ;
        const real frictionVelocity=vTangentialMagnitude/uPlus ;
        const real tauWallMagnitude=rho_f[face]*pow(frictionVelocity,2) ;
        tauWall[face]=-tauWallMagnitude*vTangential/vTangentialMagnitude ;

        // Compute the wall heat flux.
        qWall[face]=(temperature[ci[face]]/temperature_f[face]-1.0+
          recoveryFactor*pow(vTangentialMagnitude,2.0)/
          (2.0*cp[ci[face]]*temperature_f[face]))*(rho_f[face]*
          temperature_f[face]*thermalConductivity[ci[face]]*
          pow(frictionVelocity,2.0))/(laminarViscosity[ci[face]]*
          vTangentialMagnitude) ;

        // Compute the eddy viscosity for the cell next to the wall.
        eddyViscosityWall[face]=laminarViscosity_f[face]*kappa*expKB*(-1-
          kappa*uPlus-pow(kappa*uPlus,2.0)/2.0) ;

        // Compute the compressibility and heat transfer constants.
        const real gamma=recoveryFactor*pow(frictionVelocity,2.0)/(2.0*
          cp[ci[face]]*temperature_f[face]) ;
        const real beta=qWall[face]*laminarViscosity_f[face]/(rho_f[face]*
          temperature_f[face]*thermalConductivity[ci[face]]*frictionVelocity) ;
        const real Q=sqrt(beta*beta+4.0*gamma),phi=asin(-beta/Q) ;

        /* Limit the quantitiy within asin(). This sometimes activates itself
        // during intermediate iterations, but should not be active at
        // convergence.
        real temp=sqrt(gamma)*uPlus ;
        if(temp>1.0){
          cout << "  Warning: Temporary limiting in wall function at face: "
            << face << endl ;
          uPlus/=temp ; uPlus*=0.99999 ;
        }*/

        // Compute the wall eddy viscosity.
        const real yWhite=exp(kappa/sqrt(gamma)*(asin((2.0*gamma*uPlus-beta)/
          Q)-phi))*expKB ;
        const real yWhiteDerivative=2.0*yWhite*kappa*sqrt(gamma)/Q*
          sqrt((1.0-pow(2.0*gamma*uPlus-beta,2.0)/pow(Q,2.0))) ;
        eddyViscosityWall[face]+=(1.0+yWhiteDerivative-laminarViscosity
          [ci[face]]/laminarViscosity_f[face])*laminarViscosity_f[face] ;
        eddyViscosityWall[face]=max(eddyViscosityWall[face],0.0) ;

        // Compute k and omega for the cell next to the wall.
        const real cMu=0.09,omegaI=6.0*laminarViscosity_f[face]/(0.075*
          pow(yNormal,2)*rho_f[face]),omegaO=frictionVelocity/(sqrt(cMu)*
          kappa*yNormal) ;
          omegaWall[face]=sqrt(pow(omegaI,2)+pow(omegaO,2)) ;
          kWall[face]=omegaWall[face]*eddyViscosityWall[face]/rho[ci[face]] ;
      }
                                                                                
      // Calculate temperature for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  } ;
                                                                                
  register_rule<WallLawSpecifiedTemperatureCompressible>
    registerWallLawSpecifiedTemperatureCompressible ;
}










