// Loci includes.
#include <Loci.h>
                                                                                
// StreamUns includes.
#include "const.h"
#include "sciTypes.h"

// Fluid physics library includes.
#include "reaction.h"
using fluidPhysics::reaction ;
                                                                                
namespace streamUns {

  // Creates the combustion constraint.
  class CombustionConstraint : public constraint_rule {
    private:
      const_param<reaction> reactor ;
      Constraint combustion ;
    public:

      // Define input and output.
      CombustionConstraint() {
        name_store("reactor",reactor) ;
        name_store("combustion",combustion) ;
        input("reactor") ;
        output("combustion") ;
      }

      // Set up the constraint.
      virtual void compute(const sequence& seq) {
        combustion=(reactor->num_rates()==0)? EMPTY:~EMPTY ;
      }
  } ;

  register_rule<CombustionConstraint> registerCombustionConstraint ;

  // Creates the diagonal dominance constraint.
  class DiagonalDominanceConstraint : public constraint_rule {
    private:
      const_param<string> matrixForm ;
      Constraint diagonalDominance ;
    public:

      // Define input and output.
      DiagonalDominanceConstraint() {
        name_store("matrixForm",matrixForm) ;
        name_store("diagonalDominance",diagonalDominance) ;
        input("matrixForm") ;
        output("diagonalDominance") ;
      }

      // Set up the constraint.
      virtual void compute(const sequence& seq) {
        if(*matrixForm=="natural"){
          diagonalDominance=EMPTY ;
        }else if(*matrixForm=="diagonalDominance"){
          diagonalDominance=~EMPTY ;
        }else{
          cerr << "Bad value for matrixForm in .vars file." << endl ;
          Loci::Abort() ;
        }
      }
  } ;

  register_rule<DiagonalDominanceConstraint>
    registerDiagonalDominanceConstraint ;

  // Creates the flow compressibility constraints.
  class FlowCompressibilityConstraints : public constraint_rule {
    private:
      const_param<string> flowCompressibility ;
      Constraint incompressibleFlow,compressibleFlow ;
    public:

      // Define input and output.
      FlowCompressibilityConstraints() {
        name_store("flowCompressibility",flowCompressibility) ;
        name_store("incompressibleFlow",incompressibleFlow) ;
        name_store("compressibleFlow",compressibleFlow) ;
        input("flowCompressibility") ;
        output("incompressibleFlow,compressibleFlow") ;
      }

      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        if(*flowCompressibility=="incompressible"){
          incompressibleFlow=~EMPTY ; compressibleFlow=EMPTY ;
        }else if(*flowCompressibility=="compressible"){
          incompressibleFlow=EMPTY ; compressibleFlow=~EMPTY ;
        }else{
          cerr << "Bad flowCompressibility in .vars file." << endl ;
          Loci::Abort() ;
        }
      }
  } ;

  register_rule<FlowCompressibilityConstraints>
    registerFlowCompressibilityConstraints ;

  // Creates the flow regime constraints.
  class FlowRegimeConstraints : public constraint_rule {
    private:
      const_param<string> flowRegime,transportModel ;
      Constraint laminarFlow,turbulentFlow ;
    public:

      // Define input and output.
      FlowRegimeConstraints() {
        name_store("flowRegime",flowRegime) ;
        name_store("transport_model",transportModel) ;
        name_store("laminarFlow",laminarFlow) ;
        name_store("turbulentFlow",turbulentFlow) ;
        input("flowRegime,transport_model") ;
        output("laminarFlow,turbulentFlow") ;
      }

      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        if(*transportModel=="none"){
          laminarFlow=EMPTY ; turbulentFlow=EMPTY ;
        }else{
          if(*flowRegime=="laminar"){
            laminarFlow=~EMPTY ; turbulentFlow=EMPTY ;
          }else if(*flowRegime=="turbulent"){
            laminarFlow=EMPTY ; turbulentFlow=~EMPTY ;
          }else{
            cerr << "Bad flowRegime in .vars file." << endl ; Loci::Abort() ;
          }
        }
      }
  } ;

  register_rule<FlowRegimeConstraints> registerFlowRegimeConstraints ;

  // Creates the flow type constraints.
  class FlowTypeConstraints : public constraint_rule {
    private:
      const_param<string> transportModel ;
      Constraint inviscidFlow,viscousFlow ;
    public:

      // Define input and output.
      FlowTypeConstraints() {
        name_store("transport_model",transportModel) ;
        name_store("inviscidFlow",inviscidFlow) ;
        name_store("viscousFlow",viscousFlow) ;
        input("transport_model") ;
        output("inviscidFlow,viscousFlow") ;
      }

      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        if(*transportModel=="none"){
          inviscidFlow=~EMPTY ; viscousFlow=EMPTY ;
        }else{
          inviscidFlow=EMPTY ; viscousFlow=~EMPTY ;
        }
      }
  } ;

  register_rule<FlowTypeConstraints> registerFlowTypeConstraints ;

  // Creates the inviscid flux constraints.
  class InviscidFluxConstraints : public constraint_rule {
    private:
      const_param<string> inviscidFlux ;
      Constraint fouInviscidFlux,fouOrSouInviscidFlux ;
      Constraint souInviscidFlux,souOrRoeInviscidFlux ;
      Constraint cdInviscidFlux,roeInviscidFlux ;
    public:

      // Define input and output.
      InviscidFluxConstraints() {
        name_store("inviscidFlux",inviscidFlux) ;
        name_store("fouInviscidFlux",fouInviscidFlux) ;
        name_store("fouOrSouInviscidFlux",fouOrSouInviscidFlux) ;
        name_store("souInviscidFlux",souInviscidFlux) ;
        name_store("souOrRoeInviscidFlux",souOrRoeInviscidFlux) ;
        name_store("cdInviscidFlux",cdInviscidFlux) ;
        name_store("roeInviscidFlux",roeInviscidFlux) ;
        input("inviscidFlux") ;
        output("fouInviscidFlux,fouOrSouInviscidFlux,souInviscidFlux") ;
        output("souOrRoeInviscidFlux,cdInviscidFlux,roeInviscidFlux") ;
      }

      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        if(*inviscidFlux=="FOU"){
          fouInviscidFlux=~EMPTY ; fouOrSouInviscidFlux=~EMPTY ;
          souInviscidFlux=EMPTY ; souOrRoeInviscidFlux=EMPTY ;
          cdInviscidFlux=EMPTY ; roeInviscidFlux=EMPTY ;
        }else if(*inviscidFlux=="SOU"){
          fouInviscidFlux=EMPTY ; fouOrSouInviscidFlux=~EMPTY ;
          souInviscidFlux=~EMPTY ; souOrRoeInviscidFlux=~EMPTY ;
          cdInviscidFlux=EMPTY ; roeInviscidFlux=EMPTY ;
        }else if(*inviscidFlux=="CD"){
          fouInviscidFlux=~EMPTY ; fouOrSouInviscidFlux=~EMPTY ;
          souInviscidFlux=EMPTY ; souOrRoeInviscidFlux=~EMPTY ;
          cdInviscidFlux=~EMPTY ; roeInviscidFlux=EMPTY ;
        }else if(*inviscidFlux=="Roe"){
          fouInviscidFlux=EMPTY ; fouOrSouInviscidFlux=EMPTY ;
          souInviscidFlux=EMPTY ; souOrRoeInviscidFlux=~EMPTY ;
          cdInviscidFlux=EMPTY ; roeInviscidFlux=~EMPTY ;
        }else{
          cerr << "Bad inviscidFlux in .vars file." << endl ; Loci::Abort() ;
        }
      }
  } ;

  register_rule<InviscidFluxConstraints> registerInviscidFluxConstraints ;

  // Creates the testing constraint which activates case-specific rules in
  // testing.cc .
  class TestingConstraint : public constraint_rule {
    private:
      const_param<string> testing ;
      Constraint testingConstraint ;
    public:

      // Define input and output.
      TestingConstraint() {
        name_store("testing",testing) ;
        name_store(testing->c_str(),testingConstraint) ;
        input("testing") ;
        output(testing->c_str()) ;
      }

      // Set up the constraints.
      virtual void compute(const sequence& seq) { testingConstraint=~EMPTY ; }
  } ;

  //register_rule<TestingConstraint> registerTestingConstraint ;

  // Creates the time integrator constraints.
  class TimeIntegratorConstraints : public constraint_rule {
    private:
      const_param<string> timeIntegrator ;
      Constraint BDFIntegrator,BDF2Integrator,CNIntegrator ;
    public:

      // Define input and output.
      TimeIntegratorConstraints() {
        name_store("timeIntegrator",timeIntegrator) ;
        name_store("BDFIntegrator",BDFIntegrator) ;
        name_store("BDF2Integrator",BDF2Integrator) ;
        name_store("CNIntegrator",CNIntegrator) ;
        input("timeIntegrator") ;
        output("BDFIntegrator,BDF2Integrator,CNIntegrator") ;
      }

      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        if(*timeIntegrator=="BDF"){
          BDFIntegrator=~EMPTY ; BDF2Integrator=EMPTY ; CNIntegrator=EMPTY ;
        }else if(*timeIntegrator=="BDF2"){
          BDFIntegrator=EMPTY ; BDF2Integrator=~EMPTY ; CNIntegrator=EMPTY ;
        }else if(*timeIntegrator=="CN"){
          BDFIntegrator=EMPTY ; BDF2Integrator=EMPTY ; CNIntegrator=~EMPTY ;
        }else{
          cerr << "Bad timeIntegrator in .vars file." << endl ; Loci::Abort() ;
        }
      }
  } ;

  register_rule<TimeIntegratorConstraints> registerTimeIntegratorConstraints ;

  // Creates the turbulence inviscid flux constraints.
  class TurbulenceInviscidFluxConstraints : public constraint_rule {
    private:
      const_param<string> turbulenceInviscidFlux ;
      Constraint fouTurbulence,souTurbulence ;
    public:

      // Define input and output.
      TurbulenceInviscidFluxConstraints() {
        name_store("turbulenceInviscidFlux",turbulenceInviscidFlux) ;
        name_store("fouTurbulence",fouTurbulence) ;
        name_store("souTurbulence",souTurbulence) ;
        input("turbulenceInviscidFlux") ;
        output("fouTurbulence,souTurbulence") ;
      }

      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        if(*turbulenceInviscidFlux!="FOU" && *turbulenceInviscidFlux!="SOU"){
          cerr << "Bad turbulenceInviscidFluxin .vars file." << endl ;
          Loci::Abort() ;
        }
        fouTurbulence=(*turbulenceInviscidFlux=="FOU")? ~EMPTY:EMPTY ;
        souTurbulence=(*turbulenceInviscidFlux=="SOU")? ~EMPTY:EMPTY ;
      }
  } ;

  register_rule<TurbulenceInviscidFluxConstraints>
    registerTurbulenceInviscidFluxConstraints ;

}
