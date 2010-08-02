//-----------------------------------------------------------------------------
// Description: This file contains rules specific to problems for SSC.
//
// Authors: Siddharth Thakur and Jeff Wright
//          Streamline Numerics Inc.
//          Gainesville, Florida 32609
//-----------------------------------------------------------------------------
                                                                                
// Standard library includes.
#include <vector>
using std::vector ;
                                                                                
// Loci includes.
#include <Loci.h>
using Loci::Area ;
                                                                                
// StreamUns includes.
#include "residual.h"
#include "sciTypes.h"
                                                                                
namespace streamUns {

//-----------------------------------------------------------------------------
// Rules for monitoring the mass balance of species for unsteady flow. The

  // Summation join operator for standard library vector<real>.
  class StandardVectorSummation {
    public:
      void operator()(vector<real> &a,const vector<real> &b){
        vector<real>::const_iterator bPtr=b.begin() ;
        for(vector<real>::iterator aPtr=a.begin();aPtr!=a.end();++aPtr,++bPtr)
          *aPtr+=*bPtr ;
      }
  } ;

  // Unit rule to initialize net inlet species mass flow rate.
  class NetInletSpeciesMassFlowRateUnit : public unit_rule {
    private:
      const_param<int> numSpecies ;
      param<vector<real> > netInletSpeciesMassFlowRate ;
    public:
                                                                                
      // Define input and output.
      NetInletSpeciesMassFlowRateUnit() {
        name_store("numSpecies",numSpecies) ;
        name_store("netInletSpeciesMassFlowRate",netInletSpeciesMassFlowRate) ;
        input("numSpecies") ;
        output("netInletSpeciesMassFlowRate") ;
        constraint("speciesTransport") ;
      }
                                                                                
      // Initialize the mass.
      virtual void compute(const sequence &seq) {
        *netInletSpeciesMassFlowRate=vector<real>(*numSpecies,0.0) ;
      }
  } ;
                                                                                
  register_rule<NetInletSpeciesMassFlowRateUnit>
    registerNetInletSpeciesMassFlowRateUnit ;

  // Apply rule for net inlet species mass flow rate.
  class NetInletSpeciesMassFlowRateApply : public
  apply_rule<param<vector<real> >,StandardVectorSummation> {
    private:
      const_param<int> numSpecies ;
      const_store<real> rho_f ;
      const_store<vect3d> v_f ;
      const_storeVec<real> y_f ;
      const_store<Area> area ;
      param<vector<real> > netInletSpeciesMassFlowRate ;
    private:
      vector<real> faceSpeciesMassFlowRate ;
    public:
                                                                                
      // Define input and output.
      NetInletSpeciesMassFlowRateApply() {
        name_store("numSpecies",numSpecies) ;
        name_store("rho_f",rho_f) ;
        name_store("v_f",v_f) ;
        name_store("y_f",y_f) ;
        name_store("area",area) ;
        name_store("netInletSpeciesMassFlowRate",netInletSpeciesMassFlowRate) ;
        input("numSpecies,rho_f,v_f,y_f,area") ;
        output("netInletSpeciesMassFlowRate") ;
        constraint("speciesTransport,inlet_BC") ;
      }
                                                                                
      // Add species mass for a single face.
      void calculate(Entity face) {
        for(int i=0;i<*numSpecies;++i) faceSpeciesMassFlowRate[i]=rho_f[face]*
          y_f[face][i]*dot(v_f[face],area[face].n)*area[face].sada ;
        join(*netInletSpeciesMassFlowRate,faceSpeciesMassFlowRate) ;
      }
                                                                                
      // Add species mass for all faces.
      void compute(const sequence &seq) {
        faceSpeciesMassFlowRate=vector<real>(*numSpecies,0.0) ;
        do_loop(seq,this) ;
      }
  } ;
                                                                                
  register_rule<NetInletSpeciesMassFlowRateApply>
    registerNetInletSpeciesMassFlowRateApply ;

//-----------------------------------------------------------------------------
// Rules for monitoring the mass balance of species for unsteady flow. The
// net mass of each species in the domain as well as the net inlet and outlet
// mass fluxes are tracked in time.

  // Time build rule for net inlet species mass.
  class NetInletSpeciesMassTimeBuild : public singleton_rule {
    private:
      const_param<int> numSpecies ;
      param<vector<real> > netInletSpeciesMass ;
    public:

      // Define input and output.
      NetInletSpeciesMassTimeBuild() {
        name_store("netInletSpeciesMass{n=0}",netInletSpeciesMass) ;
        name_store("numSpecies",numSpecies) ;
        input("numSpecies") ;
        output("netInletSpeciesMass{n=0}") ;
        constraint("speciesTransport,noRestart") ;
      }

      // Initialize.
      virtual void compute(const sequence &seq) {
        *netInletSpeciesMass=vector<real>(*numSpecies,0.0) ;
      }

  } ;

  register_rule<NetInletSpeciesMassTimeBuild>
    registerNetInletSpeciesMassTimeBuild ;

  // Iteration build rule for net inlet species mass.
  class NetInletSpeciesMassIterationBuild : public singleton_rule {
    private:
      const_param<vector<real> > netInletSpeciesMassOld ;
      param<vector<real> > netInletSpeciesMassNew ;
    public:

      // Define input and output.
      NetInletSpeciesMassIterationBuild() {
        name_store("netInletSpeciesMass{n=0}",netInletSpeciesMassOld) ;
        name_store("netInletSpeciesMass{n,it=0}",netInletSpeciesMassNew) ;
        input("netInletSpeciesMass{n=0}") ;
        output("netInletSpeciesMass{n,it=0}") ;
        constraint("speciesTransport") ;
      }

      // Initialize.
      virtual void compute(const sequence &seq) {
        *netInletSpeciesMassNew=*netInletSpeciesMassOld ;
      }

  } ;

  register_rule<NetInletSpeciesMassIterationBuild>
    registerNetInletSpeciesMassIterationBuild ;


  // Time advance rule for net inlet species mass.
  class NetInletSpeciesMassAdvance : public singleton_rule {
    private:
      const_param<int> numSpecies ;
      const_param<real> dt ;
      const_param<vector<real> > netInletSpeciesMassOld ;
      const_param<vector<real> > netInletSpeciesMassFlux ;
      param<vector<real> > netInletSpeciesMassNew ;
    public:

      // Define input and output.
      NetInletSpeciesMassAdvance() {
        name_store("numSpecies{n}",numSpecies) ;
        name_store("dt{n}",dt) ;
        name_store("netInletSpeciesMass{n}",netInletSpeciesMassOld) ;
        name_store("netInletSpeciesMassFlux{n}",netInletSpeciesMassFlux) ;
        name_store("netInletSpeciesMass{n+1}",netInletSpeciesMassNew) ;
        input("numSpecies{n},dt{n}") ;
        input("netInletSpeciesMass{n},netInletSpeciesMassFlux{n}") ;
        output("netInletSpeciesMass{n+1}") ;
      }

      // Initialize.
      virtual void compute(const sequence &seq) {
        *netInletSpeciesMassNew=vector<real>(*numSpecies,0.0) ;
        for(int i=0;i<*numSpecies;++i) (*netInletSpeciesMassNew)[i]=
          (*netInletSpeciesMassOld)[i]+(*netInletSpeciesMassFlux)[i]*
          (*dt) ;
      }

  } ;

  register_rule<NetInletSpeciesMassAdvance>
    registerNetInletSpeciesMassAdvance ;
}
