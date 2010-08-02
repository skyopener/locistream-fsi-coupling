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

  // Class to set the plotting flag.
  class DoPlotGridTimeDependent : public singleton_rule {
    private:
      const_param<int> n,nCycle ;
      const_param<int> it ;
      const_param<int> plotFrequency ;
      param<bool> doPlot ;
    public:

      // Define input and output.
      DoPlotGridTimeDependent() {
        name_store("$n{n}",n) ;
        name_store("ncycle{n}",nCycle) ;
        name_store("$itg{n,itg}",it) ;
        name_store("plot_freq{n,itg}",plotFrequency) ;
        name_store("do_plot_grid{n,itg}",doPlot) ;
        input("$n{n},ncycle{n},$itg{n,itg},plot_freq{n,itg}");
        output("do_plot_grid{n,itg}") ;
        constraint("gridMotionTimeDependent{n,itg}") ;
      }

      // Set the output flag.
      virtual void compute(const sequence &seq) {
        if((*n)==0){ doPlot=false ; return ; }
        doPlot=((*nCycle)%(*plotFrequency)==0 && (*it)==0) ;
      }
  } ;

  register_rule<DoPlotGridTimeDependent> registerDoPlotGridTimeDependent ;

  // Class to set the plotting flag.
  class DoPlotGridSolutionDependent : public singleton_rule {
    private:
      const_param<int> n,nCycle ;
      const_param<int> it ;
      const_param<int> plotFrequency ;
      param<bool> doPlot ;
    public:

      // Define input and output.
      DoPlotGridSolutionDependent() {
        name_store("$n{n}",n) ;
        name_store("ncycle{n}",nCycle) ;
        name_store("$it{n,it}",it) ;
        name_store("plot_freq{n,it}",plotFrequency) ;
        name_store("do_plot_grid{n,it}",doPlot) ;
        input("$n{n},ncycle{n},$it{n,it},plot_freq{n,it}");
        output("do_plot_grid{n,it}") ;
        constraint("gridMotionSolutionDependent{n,it}") ;
      }

      // Set the output flag.
      virtual void compute(const sequence &seq) {
        if((*n)==0){ doPlot=false ; return ; }
        doPlot=((*nCycle)%(*plotFrequency)==0 && (*it)==0) ;
      }
  } ;

  register_rule<DoPlotGridSolutionDependent> registerDoPlotGridSolutionDependent ;

  // Class to set the plotting flag.
  class DoPlotLastIterTimeDependent : public singleton_rule {
    private:
      const_param<int> n,nCycle ;
      const_param<int> it ;
      const_param<int> gridMoverMaxIterationsPerTimeStep ;
      const_param<int> plotFrequency ;
      param<bool> doPlot ;
    public:

      // Define input and output.
      DoPlotLastIterTimeDependent() {
        name_store("$n{n}",n) ;
        name_store("ncycle{n}",nCycle) ;
        name_store("$itg{n,itg}",it) ;
        name_store("gridMoverMaxIterationsPerTimeStep{n,itg}",
          gridMoverMaxIterationsPerTimeStep) ;
        name_store("plot_freq{n,itg}",plotFrequency) ;
        name_store("do_plot_last_iter{n,itg}",doPlot) ;
        input("$n{n},ncycle{n},$itg{n,itg},plot_freq{n,itg}");
        input("gridMoverMaxIterationsPerTimeStep{n,itg}") ;
        output("do_plot_last_iter{n,itg}") ;
        constraint("gridMotionTimeDependent{n,itg}") ;
      }

      // Set the output flag.
      virtual void compute(const sequence &seq) {
        doPlot=((*nCycle)%(*plotFrequency)==0 &&
          (*it)==(*gridMoverMaxIterationsPerTimeStep-1)) ;
      }
  } ;

  register_rule<DoPlotLastIterTimeDependent> registerDoPlotLastIterTimeDependent ;

  // Copy of Ed's function with a few minor modifications to make similar to
  // our old dump_scalar().
  void solver_dump_var_local(const sequence &seq,Loci::storeRepP var,
  const_param<int> &ncyc,const_param<int> &ncycle,const_param<int>
  &plot_modulo,const_param<string> &modelName,string type,string sname) {

    ostringstream oss ; int cycle = *ncycle ;
    if(*plot_modulo != 0) cycle = cycle % *plot_modulo ;

    oss << "output/" << sname << "_" << type << "." << cycle
      << "_" << *modelName ;
    string filename = oss.str() ;
    if(Loci::MPI_rank == 0) cout << "Writing file: " << filename << endl ;
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
        string var_name_time = var_name + "{n,itg}" ;
        string constraint_name = string("scalarOutput_") + value_name ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n,itg}",plot_modulo) ;
        name_store("modelName{n,itg}",modelName) ;
        name_store("OUTPUT{n,itg}",OUTPUT) ;
        conditional("do_plot_grid{n,itg}") ;
        constraint("pos{n,itg}") ;
        input("timeStepNumber{n},ncycle{n},plot_modulo{n,itg}") ;
        input("modelName{n,itg}") ;
        input(var_name_time);
        constraint(constraint_name) ;
        output("OUTPUT{n,itg}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
        solver_dump_var_local(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "sca",value_name) ;
      }

  } ;

  // Writes out a nodal scalar.
  class scalar_node_output_last_iter : public pointwise_rule {
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
      scalar_node_output_last_iter(const char *vname,const char *valname) {
        var_name = string(vname) ;
        value_name = string(valname) ;
        string var_name_time = var_name + "{n,itg}" ;
        string constraint_name = string("scalarOutput_") + value_name ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n,itg}",plot_modulo) ;
        name_store("modelName{n,itg}",modelName) ;
        name_store("OUTPUT{n,itg}",OUTPUT) ;
        conditional("do_plot_last_iter{n,itg}") ;
        constraint("pos{n,itg}") ;
        input("timeStepNumber{n},ncycle{n},plot_modulo{n,itg}") ;
        input("modelName{n,itg}") ;
        input(var_name_time);
        constraint(constraint_name) ;
        output("OUTPUT{n,itg}") ;
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
        string var_name_time = var_name + "{n,itg}" ;
        string constraint_name = string("scalarOutput_") + value_name ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n,itg}",plot_modulo) ;
        name_store("modelName{n,itg}",modelName) ;
        name_store("OUTPUT{n,itg}",OUTPUT) ;
        conditional("do_plot_grid{n,itg}") ;
        constraint("pos{n,itg}") ;
        input("timeStepNumber{n},ncycle{n}") ;
        input("plot_modulo{n,itg},modelName{n,itg}") ;
        input(var_name_time);
        constraint(constraint_name) ;
        output("OUTPUT{n,itg}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
        solver_dump_var_local(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "vec",value_name) ;
      }

  } ;

  // Writes out a nodal vector residual.
  class vector_node_output_last_iter : public pointwise_rule {
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
      vector_node_output_last_iter(const char *vname,const char *valname) {
        var_name = string(vname) ;
        value_name = string(valname) ;
        string var_name_time = var_name + "{n,itg}" ;
        string constraint_name = string("scalarOutput_") + value_name ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n,itg}",plot_modulo) ;
        name_store("modelName{n,itg}",modelName) ;
        name_store("OUTPUT{n,itg}",OUTPUT) ;
        conditional("do_plot_last_iter{n,itg}") ;
        constraint("pos{n,itg}") ;
        input("timeStepNumber{n},ncycle{n}") ;
        input("plot_modulo{n,itg},modelName{n,itg}") ;
        input(var_name_time);
        constraint(constraint_name) ;
        output("OUTPUT{n,itg}") ;
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

#define OUTPUT_VECTOR_LAST_ITER(X,Y) class OUT_LOCAL_##Y : public \
        vector_node_output_last_iter {\
        public:\
        OUT_LOCAL_##Y() : vector_node_output_last_iter(X,#Y){}\
        }; register_rule<OUT_LOCAL_##Y> register_OUT_LOCAL_##Y ;

  OUTPUT_SCALAR_LOCAL("sMag",sMag) ;
  OUTPUT_SCALAR_LOCAL("node_chi",node_chi) ;
  OUTPUT_SCALAR_LOCAL("processNumber",processNumber) ;
  OUTPUT_VECTOR_LOCAL("node_s",node_s) ;
  OUTPUT_VECTOR_LAST_ITER("node_sStar",node_sStar) ;
  OUTPUT_VECTOR_LAST_ITER("sResidual",sResidual) ;

}
