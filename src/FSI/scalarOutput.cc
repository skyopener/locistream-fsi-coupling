//-----------------------------------------------------------------------------
// Description: This file contains rules for outputting data when running the
//   grid movement schemes in standalone mode.
//-----------------------------------------------------------------------------
                                                                                
// Standard library includes.
#include <string>
#include <sstream>
#include <stdlib.h>
using std::string ;
                                                                                
// Loci includes.
#include <Loci.h>
                                                                                
// StreamUns includes.
#include "sciTypes.h"
                                                                                
namespace streamUns {

  // Build rule for the solution time.
  class BuildSolutionTimeGridOnly : public singleton_rule {
    private:
      const_param<real> initialSolutionTime ;
      param<real> solutionTime ;
    public:

      // Define input and output.
      BuildSolutionTimeGridOnly() {
        name_store("stime_init",initialSolutionTime) ;
        name_store("stimeGridOnly{nn=0}",solutionTime) ;
        input("stime_init") ;
        output("stimeGridOnly{nn=0}") ;
        constraint("UNIVERSE") ;
      }

      // Initialize the solution time.
      virtual void compute(const sequence &seq) {
        *solutionTime=*initialSolutionTime ;
      }
  } ;

  register_rule<BuildSolutionTimeGridOnly> registerBuildSolutionTimeGridOnly ;

  // Rule to advance the solution time.
  class AdvanceSolutionTimeGridOnly : public singleton_rule {
    private:
      const_param<real> oldSolutionTime ;
      const_param<real> timeStep ;
      param<real> newSolutionTime ;
    public:
                                                                                
      // Define input and output.
      AdvanceSolutionTimeGridOnly() {
        name_store("stimeGridOnly{nn}",oldSolutionTime) ;
        name_store("timeStep{nn}",timeStep) ;
        name_store("stimeGridOnly{nn+1}",newSolutionTime) ;
        input("stimeGridOnly{nn},timeStep{nn}") ;
        output("stimeGridOnly{nn+1}") ;
        constraint("UNIVERSE") ;
      }
                                                                                
      // Increment the solution time.
      virtual void compute(const sequence &seq) {
        *newSolutionTime=*oldSolutionTime+*timeStep ;
      }
  } ;
                                                                                
  register_rule<AdvanceSolutionTimeGridOnly>
    registerAdvanceSolutionTimeGridOnly ;

  // Rule to collapse the solution time.
  class CollapseSolutionTimeGridOnly : public singleton_rule {
    private:
      const_param<real> solutionTime ;
      param<real> finalSolutionTime ;
    public:
                                                                                
      // Define input and output.
      CollapseSolutionTimeGridOnly() {
        name_store("stimeGridOnly{nn}",solutionTime) ;
        name_store("stimeGridOnly",finalSolutionTime) ;
        input("stimeGridOnly{nn}") ;
        output("stimeGridOnly") ;
        conditional("timeSteppingFinished{nn}") ;
        constraint("geom_cells{nn}") ;
      }
                                                                                
      // Set the final solution time value.
      virtual void compute(const sequence &seq) {
        *finalSolutionTime=*solutionTime ;
      }
  } ;
                                                                                
  register_rule<CollapseSolutionTimeGridOnly>
  registerCollapseSolutionTimeGridOnly ;

  // Temporary rule to set initial time-step number.
  class InitialTimeStepNumGridOnly : public singleton_rule {
    private:
      param<int> initialTimeStepNum ;
    public:

      // Define input and output.
      InitialTimeStepNumGridOnly() {
        name_store("ncycleInitGridOnly",initialTimeStepNum) ;
        output("ncycleInitGridOnly") ;
        constraint("UNIVERSE") ;
      }

      // Set the initial time-step number.
      virtual void compute(const sequence &seq) { *initialTimeStepNum=0 ; }
  } ;

  register_rule<InitialTimeStepNumGridOnly>
   registerInitialTimeStepNumGridOnly ;

  // Build rule for the time-step counter.
  class BuildTimeStepCounterGridOnly : public singleton_rule {
    private:
      const_param<int> initialTimeStepNum ;
      param<int> timeStepNum ;
    public:

      // Define input and output.
      BuildTimeStepCounterGridOnly() {
        name_store("ncycleInitGridOnly",initialTimeStepNum) ;
        name_store("ncycleGridOnly{nn=0}",timeStepNum) ;
        input("ncycleInitGridOnly") ;
        output("ncycleGridOnly{nn=0}") ;
        constraint("UNIVERSE") ;
      }

      // Initialize the time step number.
      virtual void compute(const sequence &seq) {
        *timeStepNum=*initialTimeStepNum ;
      }
  } ;

  register_rule<BuildTimeStepCounterGridOnly>
    registerBuildTimeStepCounterGridOnly ;

  // Rule to advance the time-step counter.
  class AdvanceTimeStepCounterGridOnly : public singleton_rule {
    private:
      const_param<int> timeStepNum ;
      param<int> timeStepNumPlusOne ;
    public:

      // Define input and output.
      AdvanceTimeStepCounterGridOnly() {
        name_store("ncycleGridOnly{nn}",timeStepNum) ;
        name_store("ncycleGridOnly{nn+1}",timeStepNumPlusOne) ;
        input("ncycleGridOnly{nn}") ;
        output("ncycleGridOnly{nn+1}") ;
        constraint("UNIVERSE") ;
      }

      // Increment the variable.
      virtual void compute(const sequence &seq) {
        *timeStepNumPlusOne=*timeStepNum+1 ;
      }
  } ;

  register_rule<AdvanceTimeStepCounterGridOnly>
    registerAdvanceTimeStepCounterGridOnly ;

  // Rule to collapse the time-step counter.
  class CollapseTimeStepCounterGridOnly : public singleton_rule {
    private:
      const_param<int> timeStepNum ;
      param<int> finalTimeStepNum ;
    public:
                                                                                
      // Define input and output.
      CollapseTimeStepCounterGridOnly() {
        name_store("ncycleGridOnly{nn}",timeStepNum) ;
        name_store("ncycleGridOnly",finalTimeStepNum) ;
        input("ncycleGridOnly{nn}") ;
        output("ncycleGridOnly") ;
        conditional("timeSteppingFinished{nn}") ;
        constraint("nodes{nn}") ;
      }
                                                                                
      // Set the final time step value.
      virtual void compute(const sequence &seq) {
        *finalTimeStepNum=*timeStepNum ;
      }
  } ;
                                                                                
  register_rule<CollapseTimeStepCounterGridOnly>
    registerCollapseTimeStepCounterGridOnly ;

  // Class to set the plotting flag.
  class DoPlotGridOnly : public singleton_rule {
    private:
      const_param<int> nn ;
      const_param<int> it ;
      const_param<int> gridMoverMaxIterationsPerTimeStep ;
      const_param<int> plotFrequency ;
      param<bool> doPlot ;
    public:
                                                                                
      // Define input and output.
      DoPlotGridOnly() {
        name_store("$nn{nn}",nn) ;
        name_store("$itg{nn,itg}",it) ;
        name_store("gridMoverMaxIterationsPerTimeStep{nn,itg}",
          gridMoverMaxIterationsPerTimeStep) ;
        name_store("plot_freq{nn,itg}",plotFrequency) ;
        name_store("do_plot{nn,itg}",doPlot) ;
        input("$nn{nn},$itg{nn,itg}") ;
        input("gridMoverMaxIterationsPerTimeStep{nn,itg},plot_freq{nn,itg}");
        output("do_plot{nn,itg}") ;
      }
                                                                                
      // Set the output flag. Note that we used to have maxIter-1 here. This
      // is required for the old collapse rules. 5/24/05
      virtual void compute(const sequence &seq) {
        doPlot=((*nn+1)%(*plotFrequency)==0 && *it==
          *gridMoverMaxIterationsPerTimeStep-1);
      }
  } ;
                                                                                
  register_rule<DoPlotGridOnly> registerDoPlotGridOnly ;

  // Hack to get output stuff to work.
  class TimeStepNumberGridOnly : public singleton_rule {
    private:
      const_param<int> ncyc ;
      param<int> timeStepNumber ;
    public:

      TimeStepNumberGridOnly() {
        name_store("$nn{nn}",ncyc) ;
        name_store("timeStepNumber{nn}",timeStepNumber) ;
        input("$nn{nn}") ;
        output("timeStepNumber{nn}") ;
      }

      void compute(const sequence &seq) { *timeStepNumber=*ncyc ; }
  } ;

  register_rule<TimeStepNumberGridOnly> registerTimeStepNumberGridOnly ;

  // Copy of Ed's function with a few minor modifications to make similar to
  // our old dump_scalar().
  void solver_dump_var_local(const sequence &seq,Loci::storeRepP var,
  const_param<int> &ncyc,const_param<int> &ncycle,const_param<int>
  &plot_modulo,const_param<string> &modelName,string type,string sname) {
    if(*ncyc<=-1) return ;
    ostringstream oss ; int cycle = *ncycle+1 ;
    if(*plot_modulo != 0) cycle = cycle % *plot_modulo ;
    oss << "output/" << sname << "_" << type << "." << cycle
      << "_" << *modelName ;
    string filename = oss.str() ;
    if(Loci::MPI_rank == 0) cout << "writing file " << filename << endl ;
     hid_t file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
       H5P_DEFAULT, H5P_DEFAULT) ;
    Loci::writeContainer(file_id,sname,var) ;
    Loci::hdf5CloseFile(file_id) ;
  }

  // Writes out a nodal scalar.
  class scalar_node_output_local : public pointwise_rule {
    private:
      string var_name ;
      string value_name ;
      const_store<float> c2n ;
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
      blackbox<int> blackboxOutput ;
    public:

      // Define input and output.
      scalar_node_output_local(const char *vname,const char *valname) {
        var_name = string(vname) ;
        value_name = string(valname) ;
        string var_name_time = var_name + "{nn,itg}" ;
        string constraint_name = string("scalarOutput_") + value_name ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{nn}",ncyc) ;
        name_store("ncycleGridOnly{nn}",ncycle) ;
        name_store("plot_modulo{nn,itg}",plot_modulo) ;
        name_store("modelName{nn,itg}",modelName) ;
        name_store("OUTPUT{nn,itg}",OUTPUT) ;
        conditional("do_plot{nn,itg}") ;
        constraint("pos{nn,itg}") ;
        input("timeStepNumber{nn},ncycleGridOnly{nn},plot_modulo{nn,itg}") ;
        input("modelName{nn,itg}") ;
        input(var_name_time);
        constraint(constraint_name) ;
        output("OUTPUT{nn,itg}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
        solver_dump_var_local(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "sca",value_name) ;
      }

  } ;

  // Writes out a nodal vector.
  class vector_node_output_local : public pointwise_rule {
    private:
      string var_name ;
      string value_name ;
      const_store<vector3d<real> > c2n ; // Changed to real. 8/3/07 JW
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      vector_node_output_local(const char *vname,const char *valname) {
        var_name = string(vname) ; 
        value_name = string(valname) ;
        string var_name_time = var_name + "{nn,itg}" ;
        string constraint_name = string("scalarOutput_") + value_name ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{nn}",ncyc) ;
        name_store("ncycleGridOnly{nn}",ncycle) ;
        name_store("plot_modulo{nn,itg}",plot_modulo) ;
        name_store("modelName{nn,itg}",modelName) ;
        name_store("OUTPUT{nn,itg}",OUTPUT) ;
        conditional("do_plot{nn,itg}") ;
        constraint("pos{nn,itg}") ;
        input("timeStepNumber{nn},ncycleGridOnly{nn}") ;
        input("plot_modulo{nn,itg},modelName{nn,itg}") ;
        input(var_name_time);
        constraint(constraint_name) ;
        output("OUTPUT{nn,itg}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
        solver_dump_var_local(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "vec",value_name) ;
      }

  } ;

#define OUTPUT_SCALAR_LOCAL(X,Y) class OUT_LOCAL_##Y : public \
	scalar_node_output_local {\
        public:\
        OUT_LOCAL_##Y() : scalar_node_output_local(X,#Y){}\
        }; register_rule<OUT_LOCAL_##Y> register_OUT_LOCAL_##Y ;

#define OUTPUT_VECTOR_LOCAL(X,Y) class OUT_LOCAL_##Y : public \
        vector_node_output_local {\
        public:\
        OUT_LOCAL_##Y() : vector_node_output_local(X,#Y){}\
        }; register_rule<OUT_LOCAL_##Y> register_OUT_LOCAL_##Y ;

  OUTPUT_SCALAR_LOCAL("sMag",sMag) ;
  OUTPUT_SCALAR_LOCAL("node_chi",node_chi) ;
  OUTPUT_VECTOR_LOCAL("node_s",node_s) ;
//OUTPUT_VECTOR_LOCAL("sResidual",sResidual) ;

}
