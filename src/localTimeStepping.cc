// Loci includes.
#include <Loci.h>
                                                                                
// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;

// StreamUns includes.
#include "sciTypes.h"

namespace streamUns {
//-----------------------------------------------------------------------------

  // Set default time-step factor so no local time-stepping happens.
  class DefaultTimeStepFactor : public pointwise_rule {
    private:
      store<real> timeStepFactor ;
    public:
                                                                                
      // Define input and output.
      DefaultTimeStepFactor() {
        name_store("timeStepFactor",timeStepFactor) ;
        output("timeStepFactor") ;
        constraint("geom_cells") ;
      }
                                                                                
      // Set for a cell.
      void calculate(Entity cell) { timeStepFactor[cell]=1.0 ; }
                                                                                
      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<DefaultTimeStepFactor> registerDefaultTimeStepFactor ;

//-----------------------------------------------------------------------------
// Rules to compute the local time-stepping factor by monitoring a specific
// species mass fraction (a radical more than likely).

  // Sets the default species to monitor.
  class DefaultMonitorSpecies : public default_rule {
    private:
      param<string> monitorSpecies ;
    public:

      // Define input and output.
      DefaultMonitorSpecies() {
        name_store("monitorSpecies",monitorSpecies) ;
        output("monitorSpecies") ;
        comments("Sets the default monitor species to OH.") ;
      }

      // Set the default species.
      virtual void compute(const sequence &seq) { *monitorSpecies="OH" ; }
  } ;

  register_rule<DefaultMonitorSpecies> registerDefaultMonitorSpecies ;

  // Sets the power of the monitor species time curve.
  class DefaultMonitorSpeciesPower : public default_rule {
    private:
      param<real> monitorSpeciesPower ;
    public:

      // Define input and output.
      DefaultMonitorSpeciesPower() {
        name_store("monitorSpeciesPower",monitorSpeciesPower) ;
        output("monitorSpeciesPower") ;
        comments("Sets the default monitor species power to 2.0.") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) { *monitorSpeciesPower=2.0 ; }
  } ;

  register_rule<DefaultMonitorSpeciesPower>
    registerDefaultMonitorSpeciesPower ;

  // Sets the lower cut-off value for the monitor species mass fraction. Below
  // this value (which is actually the log of the species mass fraction),
  // dT=timeStep, where timeStep is that specified in the .vars file.
  class DefaultLowerMonitorSpeciesCutoff : public default_rule {
    private:
      param<real> lowerMonitorSpeciesCutoff ;
    public:

      // Define input and output.
      DefaultLowerMonitorSpeciesCutoff() {
        name_store("lowerMonitorSpeciesCutoff",lowerMonitorSpeciesCutoff) ;
        output("lowerMonitorSpeciesCutoff") ;
        comments("Sets the default lower monitor species cut-off to -20.0.") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) {
        *lowerMonitorSpeciesCutoff=-20.0 ;
      }
  } ;

  register_rule<DefaultLowerMonitorSpeciesCutoff>
    registerDefaultLowerMonitorSpeciesCutoff ;

  // Sets the upper cut-off value for the monitor species mass fraction. Above
  // this value dT=minTimeStep, where minTimeStep is that specified in the
  // .vars file.
  class DefaultUpperMonitorSpeciesCutoff : public default_rule {
    private:
      param<real> upperMonitorSpeciesCutoff ;
    public:

      // Define input and output.
      DefaultUpperMonitorSpeciesCutoff() {
        name_store("upperMonitorSpeciesCutoff",upperMonitorSpeciesCutoff) ;
        output("upperMonitorSpeciesCutoff") ;
        comments("Sets the default upper monitor species cut-off to -4.0.") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) {
        *upperMonitorSpeciesCutoff=-4.0 ;
      }
  } ;

  register_rule<DefaultUpperMonitorSpeciesCutoff>
    registerDefaultUpperMonitorSpeciesCutoff ;

  // Minimum time step for the .vars file.
  class OptionalMinTimeStep : public optional_rule {
    private:
      param<real> minTimeStep ;
    public:

      // Define input and output.
      OptionalMinTimeStep() {
        name_store("minTimeStep",minTimeStep) ;
        output("minTimeStep") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<OptionalMinTimeStep> registerOptionalMinTimeStep ;

  // Computes the local time-stepping factor based on the mass fraction of
  // the monitor species.
  class MonitorSpeciesTimeStepFactor : public pointwise_rule {
    private:
      const_param<EOS> eos ;
      const_param<string> monitorSpecies ;
      const_param<real> monitorSpeciesPower ;
      const_param<real> lowerMonitorSpeciesCutoff ;
      const_param<real> upperMonitorSpeciesCutoff ;
      const_param<real> minTimeStep ;
      const_param<real> timeStep ;
      const_storeVec<real> y ;
      store<real> timeStepFactor ;
    private:
      int speciesIndex ;
    public:
                                                                                
      // Define input and output.
      MonitorSpeciesTimeStepFactor() {
        name_store("eos",eos) ;
        name_store("monitorSpecies",monitorSpecies) ;
        name_store("monitorSpeciesPower",monitorSpeciesPower) ;
        name_store("lowerMonitorSpeciesCutoff",lowerMonitorSpeciesCutoff) ;
        name_store("upperMonitorSpeciesCutoff",upperMonitorSpeciesCutoff) ;
        name_store("minTimeStep",minTimeStep) ;
        name_store("timeStep",timeStep) ;
        name_store("y",y) ;
        name_store("monitorSpecies::timeStepFactor",timeStepFactor) ;
        input("eos,monitorSpecies,monitorSpeciesPower") ;
        input("lowerMonitorSpeciesCutoff,upperMonitorSpeciesCutoff") ;
        input("minTimeStep,timeStep,y") ;
        output("monitorSpecies::timeStepFactor") ;
        constraint("geom_cells") ;
      }
                                                                                
      // Compute the value for a cell.
      void calculate(Entity cell) {
        real logY=(y[cell][speciesIndex]<1.0e-30)? -1000.0:
          log10(y[cell][speciesIndex]) ;
        if(logY<(*lowerMonitorSpeciesCutoff)){
          timeStepFactor[cell]=1.0 ;
        }else if(logY>(*upperMonitorSpeciesCutoff)){
          timeStepFactor[cell]=(*minTimeStep)/(*timeStep) ;
        }else{
          timeStepFactor[cell]=((*minTimeStep)/(*timeStep)-1.0)*pow(logY-
            (*lowerMonitorSpeciesCutoff),*monitorSpeciesPower)/
            pow(*upperMonitorSpeciesCutoff-*lowerMonitorSpeciesCutoff,
            *monitorSpeciesPower)+1.0 ;
        }
if(timeStepFactor[cell]<0.0){
  cerr << "ERROR: timeStepFactor is negative!" << endl ;
  cerr << "cell,timeStepFactor: " << cell << " " << timeStepFactor[cell]
    << endl ;
  cerr << "y[cell][speciesIndex]: " << y[cell][speciesIndex] << endl ;
  cerr << "logY: " << logY << endl ;
  Loci::Abort() ;
}
      }
                                                                                
      // Loop over cells.
      virtual void compute(const sequence &seq) {
        speciesIndex=eos->speciesIndex(*monitorSpecies) ;
        if(speciesIndex==-1){
          cerr << "ERROR: Species " << *monitorSpecies << " does not exist in "
            << " MonitorSpeciesTimeStepFactor()." << endl ; Loci::Abort() ;
        }
        do_loop(seq,this) ;
      }
  } ;
                                                                                
  register_rule<MonitorSpeciesTimeStepFactor>
    registerMonitorSpeciesTimeStepFactor ;

}
