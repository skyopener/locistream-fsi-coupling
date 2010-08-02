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

// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;

// StreamUns includes.
#include "sciTypes.h"

//#define ITERATIONOUTPUT
//#define DOUBLEPRECISIONOUTPUT

namespace streamUns {

  // Temporary rule to set initial time-step number.
  class InitialTimeStepNum : public singleton_rule {
    private:
      param<int> initialTimeStepNum ;
    public:

      // Define input and output.
      InitialTimeStepNum() {
        name_store("ncycle_init",initialTimeStepNum) ;
        output("ncycle_init") ;
        constraint("noRestart") ;
      }

      // Set the initial time-step number.
      virtual void compute(const sequence &seq) { *initialTimeStepNum=0 ; }
  } ;

  register_rule<InitialTimeStepNum> registerInitialTimeStepNum ;

  // Build rule for the time-step counter.
  class BuildTimeStepCounter : public singleton_rule {
    private:
      const_param<int> initialTimeStepNum ;
      param<int> timeStepNum ;
    public:

      // Define input and output.
      BuildTimeStepCounter() {
        name_store("ncycle_init",initialTimeStepNum) ;
        name_store("ncycle{n=0}",timeStepNum) ;
        input("ncycle_init") ;
        output("ncycle{n=0}") ;
//      constraint("UNIVERSE") ;
      }

      // Initialize the time step number.
      virtual void compute(const sequence &seq) {
        *timeStepNum=*initialTimeStepNum ;
      }
  } ;

  register_rule<BuildTimeStepCounter> registerBuildTimeStepCounter ;

  // Rule to advance the time-step counter.
  class AdvanceTimeStepCounter : public singleton_rule {
    private:
      const_param<int> timeStepNum ;
      param<int> timeStepNumPlusOne ;
    public:

      // Define input and output.
      AdvanceTimeStepCounter() {
        name_store("ncycle{n}",timeStepNum) ;
        name_store("ncycle{n+1}",timeStepNumPlusOne) ;
        input("ncycle{n}") ;
        output("ncycle{n+1}") ;
//      constraint("UNIVERSE") ;
      }

      // Increment the variable.
      virtual void compute(const sequence &seq) {
        *timeStepNumPlusOne=*timeStepNum+1 ;
      }
  } ;

  register_rule<AdvanceTimeStepCounter> registerAdvanceTimeStepCounter ;

  // Rule to collapse the time-step counter.
  class CollapseTimeStepCounter : public singleton_rule {
    private:
      const_param<int> timeStepNum ;
      param<int> finalTimeStepNum ;
    public:

      // Define input and output.
      CollapseTimeStepCounter() {
        name_store("ncycle{n}",timeStepNum) ;
        name_store("ncycle",finalTimeStepNum) ;
        input("ncycle{n}") ;
        output("ncycle") ;
        conditional("timeSteppingFinished{n}") ;
//      constraint("geom_cells{n}") ;
      }

      // Set the final time step value.
      virtual void compute(const sequence &seq) {
        *finalTimeStepNum=*timeStepNum ;
      }
  } ;

  register_rule<CollapseTimeStepCounter> registerCollapseTimeStepCounter ;

  // Temporary rule to set initial solution time to zero.
  class InitialSolutionTime : public singleton_rule {
    private:
      param<real> initialSolutionTime ;
    public:

      // Define input and output.
      InitialSolutionTime() {
        name_store("stime_init",initialSolutionTime) ;
        output("stime_init") ;
        constraint("noRestart") ;
      }

      // Set the initial solution time.
      virtual void compute(const sequence &seq) { *initialSolutionTime=0.0 ; }
  } ;

  register_rule<InitialSolutionTime> registerInitialSolutionTime ;

  // Build rule for the solution time.
  class BuildSolutionTime : public singleton_rule {
    private:
      const_param<real> initialSolutionTime ;
      param<real> solutionTime ;
    public:

      // Define input and output.
      BuildSolutionTime() {
        name_store("stime_init",initialSolutionTime) ;
        name_store("stime{n=0}",solutionTime) ;
        input("stime_init") ;
        output("stime{n=0}") ;
//      constraint("UNIVERSE") ;
      }

      // Initialize the solution time.
      virtual void compute(const sequence &seq) {
        *solutionTime=*initialSolutionTime ;
      }
  } ;

  register_rule<BuildSolutionTime> registerBuildSolutionTime ;

  // Rule to advance the solution time.
  class AdvanceSolutionTime : public singleton_rule {
    private:
      const_param<real> oldSolutionTime ;
      const_param<real> dt ;
      param<real> newSolutionTime ;
    public:

      // Define input and output.
      AdvanceSolutionTime() {
        name_store("stime{n}",oldSolutionTime) ;
        name_store("dt{n}",dt) ;
        name_store("stime{n+1}",newSolutionTime) ;
        input("stime{n},dt{n}") ;
        output("stime{n+1}") ;
//      constraint("UNIVERSE") ;
      }

      // Increment the solution time.
      virtual void compute(const sequence &seq) {
        *newSolutionTime=*oldSolutionTime+*dt ;
      }
  } ;

  register_rule<AdvanceSolutionTime> registerAdvanceSolutionTime ;

  // Rule to collapse the solution time.
  class CollapseSolutionTime : public singleton_rule {
    private:
      const_param<real> solutionTime ;
      param<real> finalSolutionTime ;
    public:

      // Define input and output.
      CollapseSolutionTime() {
        name_store("stime{n}",solutionTime) ;
        name_store("stime",finalSolutionTime) ;
        input("stime{n}") ;
        output("stime") ;
        conditional("timeSteppingFinished{n}") ;
//      constraint("geom_cells{n}") ;
      }

      // Set the final solution time value.
      virtual void compute(const sequence &seq) {
        *finalSolutionTime=*solutionTime ;
      }
  } ;

  register_rule<CollapseSolutionTime> registerCollapseSolutionTime ;

  // Class to set the plotting flag.
  class DoPlot : public singleton_rule {
    private:
      const_param<int> n,nCycle ;
      const_param<int> it ;
      const_param<int> plotFrequency ;
      param<bool> doPlot ;
    public:

      // Define input and output.
      DoPlot() {
        name_store("$n{n}",n) ;
        name_store("ncycle{n}",nCycle) ;
        name_store("$it{n,it}",it) ;
        name_store("plot_freq{n,it}",plotFrequency) ;
        name_store("do_plot{n,it}",doPlot) ;
        input("$n{n},ncycle{n},$it{n,it},plot_freq{n,it}") ;
        output("do_plot{n,it}") ;
      }
 
      // Set the output flag.
      virtual void compute(const sequence &seq) {
        if((*n)==0){ doPlot=false ; return ; }
        doPlot=((*nCycle)%(*plotFrequency)==0 && (*it)==0) ;
      }
  } ;

  register_rule<DoPlot> registerDoPlot ;

  // Class to set the plotting flag.
  class DoPlotBoundary : public singleton_rule {
    private:
      const_param<int> n,nCycle ;
      const_param<int> it ;
      const_param<int> plotFrequency,boundaryPlotFrequency ;
      param<bool> doPlot ;
    public:

      // Define input and output.
      DoPlotBoundary() {
        name_store("$n{n}",n) ;
        name_store("ncycle{n}",nCycle) ;
        name_store("$it{n,it}",it) ;
        name_store("plot_freq{n,it}",plotFrequency) ;
        name_store("boundary_plot_freq{n,it}",boundaryPlotFrequency) ;
        name_store("do_plot_boundary{n,it}",doPlot) ;
        input("$n{n},ncycle{n},$it{n,it},plot_freq{n,it}") ;
        input("boundary_plot_freq{n,it}") ;
        output("do_plot_boundary{n,it}") ;
      }
 
      // Set the output flag.
      virtual void compute(const sequence &seq) {
        if((*n)==0){ doPlot=false ; return ; }
        int freq=min(*plotFrequency,*boundaryPlotFrequency) ;
        doPlot=((*nCycle)%freq==0 && (*it)==0) ;
      }
  } ;

  register_rule<DoPlotBoundary> registerDoPlotBoundary ;

  // Class to set the printing flag.
  class DoPrint : public singleton_rule {
    private:
      const_param<int> n,nCycle ;
      const_param<int> it ;
      const_param<int> printFrequency ;
      param<bool> doPrint ;
    public:

      // Define input and output.
      DoPrint() {
        name_store("$n{n}",n) ;
        name_store("ncycle{n}",nCycle) ;
        name_store("$it{n,it}",it) ;
        name_store("print_freq{n,it}",printFrequency) ;
        name_store("do_print{n,it}",doPrint) ;
        input("$n{n},ncycle{n},$it{n,it},print_freq{n,it}") ;
        output("do_print{n,it}") ;
      }
 
      // Set the output flag.
      virtual void compute(const sequence &seq) {
        if((*n)==0){ doPrint=false ; return ; }
        doPrint=((*nCycle)%(*printFrequency)==0 && (*it)==0) ;
      }
  } ;

  register_rule<DoPrint> registerDoPrint ;

  // Rule to put soundSpeed at boundary. Only used for post-processing.
  class SoundSpeedBoundary : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> soundSpeed ;
      store<real> soundSpeed_f ;
    public:

      // Define input and output.
      SoundSpeedBoundary() {
        name_store("ci",ci) ;
        name_store("soundSpeed",soundSpeed) ;
        name_store("soundSpeed_f",soundSpeed_f) ;
        input("ci->soundSpeed") ;
        output("soundSpeed_f") ;
      }

      // Set value for a single face.
      void calculate(Entity face) {
        soundSpeed_f[face]=soundSpeed[ci[face]] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SoundSpeedBoundary> registerSoundSpeedBoundary ;

//-----------------------------------------------------------------------------
// Rule for creating the process number nodal variable.

  class ProcessNumber : public pointwise_rule {
    private:
      store<float> processNumber ;
    public:

      // Define input and output.
      ProcessNumber() {
        name_store("processNumber",processNumber) ;
        output("processNumber") ;
        constraint("nodes") ;
      }

      // Compute for a single node.
      void calculate(Entity node) {
        processNumber[node]=float(Loci::MPI_rank) ;
      }

      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ProcessNumber> registerProcessNumber ;

//-----------------------------------------------------------------------------
// Rules for the computing the vorticity for post-processing analysis.

  class VorticityInterior : public pointwise_rule {
    private:
      const_store<tens3d> vGradient ;
      store<vect3d> vorticity ;
    public:

      // Define input and output.
      VorticityInterior() {
        name_store("gradv3d(v)",vGradient) ;
        name_store("vort",vorticity) ;
        input("gradv3d(v)") ;
        output("vort") ;
        constraint("geom_cells") ;
      }

      // Compute for a single cell.
      void calculate(Entity cell) {
        vorticity[cell]=vect3d(vGradient[cell].z.y-vGradient[cell].y.z,
          vGradient[cell].x.z-vGradient[cell].z.x,vGradient[cell].y.x-
          vGradient[cell].x.y) ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<VorticityInterior> registerVorticityInterior ;

  class VorticityBoundary : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<vect3d> vorticity ;
      store<vect3d> vorticity_f ;
    public:

      // Define input and output.
      VorticityBoundary() {
        name_store("ci",ci) ;
        name_store("vort",vorticity) ;
        name_store("vort_f",vorticity_f) ;
        input("ci->vort") ;
        output("vort_f") ;
        constraint("boundaryFaces") ;
      }

      // Compute for a single face. Simple extrapolation for now.
      void calculate(Entity face) {
        vorticity_f[face]=vorticity[ci[face]] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<VorticityBoundary> registerVorticityBoundary ;

  class VorticityMagnitudeInterior : public pointwise_rule {
    private:
      const_store<tens3d> vGradient ;
      store<real> vorticityMagnitude ;
    public:

      // Define input and output.
      VorticityMagnitudeInterior() {
        name_store("gradv3d(v)",vGradient) ;
        name_store("vort_mag",vorticityMagnitude) ;
        input("gradv3d(v)") ;
        output("vort_mag") ;
        constraint("geom_cells") ;
      }

      // Compute for a single cell.
      void calculate(Entity cell) {
        vorticityMagnitude[cell]=norm(vect3d(vGradient[cell].z.y-
          vGradient[cell].y.z,vGradient[cell].x.z-vGradient[cell].z.x,
          vGradient[cell].y.x-vGradient[cell].x.y)) ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<VorticityMagnitudeInterior> registerVorticityMagnitudeInterior ;

  class VorticityMagnitudeBoundary : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> vorticityMagnitude ;
      store<real> vorticityMagnitude_f ;
    public:

      // Define input and output.
      VorticityMagnitudeBoundary() {
        name_store("ci",ci) ;
        name_store("vort_mag",vorticityMagnitude) ;
        name_store("vort_mag_f",vorticityMagnitude_f) ;
        input("ci->vort_mag") ;
        output("vort_mag_f") ;
        constraint("boundaryFaces") ;
      }

      // Compute for a single face. Simple extrapolation for now.
      void calculate(Entity face) {
        vorticityMagnitude_f[face]=vorticityMagnitude[ci[face]] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<VorticityMagnitudeBoundary> registerVorticityMagnitudeBoundary ;

//-----------------------------------------------------------------------------
// Rules for the computing the thermal conductivity for post-processing
// analysis.

  // Thermal conductivity for cells.
  class ThermalConductivityInterior : public pointwise_rule {
    private:
      const_store<real> k ;
      store<real> thermalConductivity ;
    public:

      // Define input and output.
      ThermalConductivityInterior() {
        name_store("kconduct(temperature,p,y)",k) ;
        name_store("thermalConductivity",thermalConductivity) ;
        input("kconduct(temperature,p,y)") ;
        output("thermalConductivity") ;
        constraint("geom_cells") ;
      }

      // Compute thermal conductivity for a single cell.
      void calculate(Entity cell) { thermalConductivity[cell]=k[cell] ; }

      // Compute thermal conductivity for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ThermalConductivityInterior>
    registerThermalConductivityInterior ;

  // Thermal conductivity for boundary faces.
  class ThermalConductivityBoundary : public pointwise_rule {
    private:
      const_store<real> k_f ;
      store<real> thermalConductivity_f ;
    public:

      // Define input and output.
      ThermalConductivityBoundary() {
        name_store("kconduct(temperature_f,p_f,y_f)",k_f) ;
        name_store("thermalConductivity_f",thermalConductivity_f) ;
        input("kconduct(temperature_f,p_f,y_f)") ;
        output("thermalConductivity_f") ;
        constraint("boundaryFaces") ;
      }

      // Compute thermal conductivity for a single face.
      void calculate(Entity face) { thermalConductivity_f[face]=k_f[face] ; }

      // Compute thermal conductivity for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ThermalConductivityBoundary>
    registerThermalConductivityBoundary ;

//-----------------------------------------------------------------------------

  // Hack to get output stuff to work.
  class TimeStepNumber : public singleton_rule {
    private:
      const_param<int> ncyc ;
      param<int> timeStepNumber ;
    public:

      TimeStepNumber() {
        name_store("$n{n}",ncyc) ;
        name_store("timeStepNumber{n}",timeStepNumber) ;
        input("$n{n}") ;
        output("timeStepNumber{n}") ;
      }

      void compute(const sequence &seq) { *timeStepNumber=*ncyc ; }
  } ;

  register_rule<TimeStepNumber> registerTimeStepNumber ;

  // Copy of Ed's function with a few minor modifications to make similar to
  // our old dump_scalar().
  void solver_dump_var(const sequence &seq,Loci::storeRepP var,const_param<int>
  &ncyc,const_param<int> &ncycle,const_param<int> &plot_modulo,
  const_param<string> &modelName,string type,string sname) {

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
  class scalar_node_output : public pointwise_rule {
    private:
      string var_name ;
      string value_name ;
      const_store<float> c2n ;
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> it ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      scalar_node_output(const char *vname,const char *valname) {
        var_name = string(vname) ;
        value_name = string(valname) ;
        string var_name_time = var_name + "{n,it}" ;
        string constraint_name = string("scalarOutput_") + value_name ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("$it{n,it}",it) ;
        name_store("plot_modulo{n,it}",plot_modulo) ;
        name_store("modelName{n,it}",modelName) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
#ifndef ITERATIONOUTPUT
        conditional("do_plot{n,it}") ;
#endif
        constraint("pos{n,it}") ;
        input("timeStepNumber{n},ncycle{n},$it{n,it}") ;
        input("plot_modulo{n,it},modelName{n,it}") ;
        input(var_name_time);
        constraint(constraint_name) ;
        output("OUTPUT{n,it}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
#ifdef ITERATIONOUTPUT
        solver_dump_var(seq,c2n.Rep(),ncyc,it,plot_modulo,modelName,
          "sca",value_name) ;
#else
        solver_dump_var(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "sca",value_name) ;
#endif
      }	

  } ;

  // Always writes out a nodal scalar.
  class scalar_node_output_always : public pointwise_rule {
    private:
      string var_name ;
      string value_name ;
      const_store<float> c2n ;
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> it ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      scalar_node_output_always(const char *vname,const char *valname) {
        var_name = string(vname) ;
        value_name = string(valname) ;
        string var_name_time = var_name + "{n,it}"  ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("$it{n,it}",it) ;
        name_store("plot_modulo{n,it}",plot_modulo) ;
        name_store("modelName{n,it}",modelName) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
#ifndef ITERATIONOUTPUT
        conditional("do_plot{n,it}") ;
#endif
        constraint("pos{n,it}") ;
        input("timeStepNumber{n},ncycle{n},$it{n,it}") ;
        input("plot_modulo{n,it},modelName{n,it}") ;
        input(var_name_time);
        output("OUTPUT{n,it}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
#ifdef ITERATIONOUTPUT
        solver_dump_var(seq,c2n.Rep(),ncyc,it,plot_modulo,modelName,
          "sca",value_name) ;
#else
        solver_dump_var(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "sca",value_name) ;
#endif
      }	

  } ;

  // Writes out a nodal scalar in double precision.
  class scalard_node_output : public pointwise_rule {
    private:
      string var_name ;
      string value_name ;
      const_store<real> c2n ;
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> it ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
//    blackbox<int> blackboxOutput ;
    public:

      // Define input and output.
      scalard_node_output(const char *vname,const char *valname) {
        var_name = string(vname) ;
        value_name = string(valname) ;
        string var_name_time = var_name + "{n,it}" ;
        string constraint_name = string("scalarOutput_") + value_name ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("$it{n,it}",it) ;
        name_store("plot_modulo{n,it}",plot_modulo) ;
        name_store("modelName{n,it}",modelName) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
#ifndef ITERATIONOUTPUT
        conditional("do_plot{n,it}") ;
#endif
        constraint("pos{n,it}") ;
        input("timeStepNumber{n},ncycle{n},$it{n,it}") ;
        input("plot_modulo{n,it},modelName{n,it}") ;
        input(var_name_time);
        constraint(constraint_name) ;
        output("OUTPUT{n,it}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
#ifdef ITERATIONOUTPUT
        solver_dump_var(seq,c2n.Rep(),ncyc,it,plot_modulo,modelName,
          "sca",value_name) ;
#else
        solver_dump_var(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "sca",value_name) ;
#endif
      }	
  } ;

  // Always writes out a nodal scalar in double precision.
  class scalard_node_output_always : public pointwise_rule {
    private:
      string var_name ;
      string value_name ;
      const_store<real> c2n ;
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> it ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
//    blackbox<int> blackboxOutput ;
    public:

      // Define input and output.
      scalard_node_output_always(const char *vname,const char *valname) {
        var_name = string(vname) ;
        value_name = string(valname) ;
        string var_name_time = var_name + "{n,it}"  ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("$it{n,it}",it) ;
        name_store("plot_modulo{n,it}",plot_modulo) ;
        name_store("modelName{n,it}",modelName) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
#ifndef ITERATIONOUTPUT
        conditional("do_plot{n,it}") ;
#endif
        constraint("pos{n,it}") ;
        input("timeStepNumber{n},ncycle{n},$it{n,it}") ;
        input("plot_modulo{n,it},modelName{n,it}") ;
        input(var_name_time);
        output("OUTPUT{n,it}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
#ifdef ITERATIONOUTPUT
        solver_dump_var(seq,c2n.Rep(),ncyc,it,plot_modulo,modelName,
          "sca",value_name) ;
#else
        solver_dump_var(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "sca",value_name) ;
#endif
      }	

  } ;
  
  // Writes out a nodal vector.
  class vector_node_output : public pointwise_rule {
    private:
      string var_name ;
      string value_name ;
      const_store<vector3d<float> > c2n ;
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> it ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      vector_node_output(const char *vname,const char *valname) {
        var_name = string(vname) ;
        value_name = string(valname) ;
        string var_name_time = var_name + "{n,it}" ;
        string constraint_name = string("scalarOutput_") + value_name ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("$it{n,it}",it) ;
        name_store("plot_modulo{n,it}",plot_modulo) ;
        name_store("modelName{n,it}",modelName) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
#ifndef ITERATIONOUTPUT
        conditional("do_plot{n,it}") ;
#endif
        constraint("pos{n,it}") ;
        input("timeStepNumber{n},ncycle{n},$it{n,it}") ;
        input("plot_modulo{n,it},modelName{n,it}") ;
        input(var_name_time);
        constraint(constraint_name) ;
        output("OUTPUT{n,it}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
#ifdef ITERATIONOUTPUT
        solver_dump_var(seq,c2n.Rep(),ncyc,it,plot_modulo,modelName,
          "vec",value_name) ;
#else
        solver_dump_var(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "vec",value_name) ;
#endif
      }

  } ;
  
  // Always writes out a nodal vector.
  class vector_node_output_always : public pointwise_rule {
    private:
      string var_name ;
      string value_name ;
      const_store<vector3d<float> > c2n ;
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> it ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      vector_node_output_always(const char *vname,const char *valname) {
        var_name = string(vname) ;
        value_name = string(valname) ;
        string var_name_time = var_name + "{n,it}" ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("$it{n,it}",it) ;
        name_store("plot_modulo{n,it}",plot_modulo) ;
        name_store("modelName{n,it}",modelName) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
#ifndef ITERATIONOUTPUT
        conditional("do_plot{n,it}") ;
#endif
        constraint("pos{n,it}") ;
        input("timeStepNumber{n},ncycle{n},$it{n,it}") ;
        input("plot_modulo{n,it},modelName{n,it}") ;
        input(var_name_time);
        output("OUTPUT{n,it}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
#ifdef ITERATIONOUTPUT
        solver_dump_var(seq,c2n.Rep(),ncyc,it,plot_modulo,modelName,
          "vec",value_name) ;
#else
        solver_dump_var(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "vec",value_name) ;
#endif
      }

  } ;

  // Writes out a nodal vector in double precision.
  class vectord_node_output : public pointwise_rule {
    private:
      string var_name ;
      string value_name ;
      const_store<vector3d<real> > c2n ;
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> it ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      vectord_node_output(const char *vname,const char *valname) {
        var_name = string(vname) ;
        value_name = string(valname) ;
        string var_name_time = var_name + "{n,it}" ;
        string constraint_name = string("scalarOutput_") + value_name ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("$it{n,it}",it) ;
        name_store("plot_modulo{n,it}",plot_modulo) ;
        name_store("modelName{n,it}",modelName) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
#ifndef ITERATIONOUTPUT
        conditional("do_plot{n,it}") ;
#endif
        constraint("pos{n,it}") ;
        input("timeStepNumber{n},ncycle{n},$it{n,it}") ;
        input("plot_modulo{n,it},modelName{n,it}") ;
        input(var_name_time);
        constraint(constraint_name) ;
        output("OUTPUT{n,it}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
#ifdef ITERATIONOUTPUT
        solver_dump_var(seq,c2n.Rep(),ncyc,it,plot_modulo,modelName,
          "vec",value_name) ;
#else
        solver_dump_var(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "vec",value_name) ;
#endif
      }

  } ;

  // Always writes out a nodal vector in double precision.
  class vectord_node_output_always : public pointwise_rule {
    private:
      string var_name ;
      string value_name ;
      const_store<vector3d<real> > c2n ;
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> it ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      vectord_node_output_always(const char *vname,const char *valname) {
        var_name = string(vname) ;
        value_name = string(valname) ;
        string var_name_time = var_name + "{n,it}" ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("$it{n,it}",it) ;
        name_store("plot_modulo{n,it}",plot_modulo) ;
        name_store("modelName{n,it}",modelName) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
#ifndef ITERATIONOUTPUT
        conditional("do_plot{n,it}") ;
#endif
        constraint("pos{n,it}") ;
        input("timeStepNumber{n},ncycle{n},$it{n,it}") ;
        input("plot_modulo{n,it},modelName{n,it}") ;
        input(var_name_time);
        output("OUTPUT{n,it}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
#ifdef ITERATIONOUTPUT
        solver_dump_var(seq,c2n.Rep(),ncyc,it,plot_modulo,modelName,
          "vec",value_name) ;
#else
        solver_dump_var(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "vec",value_name) ;
#endif
      }

  } ;

#define OUTPUT_SCALAR(X,Y) class OUT_##Y : public scalar_node_output {\
                           public:\
                           OUT_##Y() : scalar_node_output(X,#Y){}\
                           }; register_rule<OUT_##Y> register_OUT_##Y ;

#define OUTPUT_SCALAR_ALWAYS(X,Y) class OUT_##Y : public \
                           scalar_node_output_always {\
                           public:\
                           OUT_##Y() : scalar_node_output_always(X,#Y){}\
                           }; register_rule<OUT_##Y> register_OUT_##Y

#define OUTPUT_SCALAR_ALWAYS_CONSTRAINT(X,Y,Z) class OUT_##Y : public \
                           scalar_node_output_always {\
                           public:\
                           OUT_##Y() : scalar_node_output_always(X,#Y){ \
                             constraint(Z) ; }\
                           }; register_rule<OUT_##Y> register_OUT_##Y

#define OUTPUT_SCALARD_ALWAYS_CONSTRAINT(X,Y,Z) class DOUT_##Y : public \
                           scalard_node_output_always {\
                           public:\
                           DOUT_##Y() : scalard_node_output_always(X,#Y){ \
                             constraint(Z) ; }\
                           }; register_rule<DOUT_##Y> register_DOUT_##Y

#define OUTPUT_SCALARD(X,Y) class DOUT_##Y : public scalard_node_output {\
                           public:\
                           DOUT_##Y() : scalard_node_output(X,#Y){}\
                           }; register_rule<DOUT_##Y> register_DOUT_##Y

#define OUTPUT_SCALARD_ALWAYS(X,Y) class DOUT_##Y : public \
                           scalard_node_output_always {\
                           public:\
                           DOUT_##Y() : scalard_node_output_always(X,#Y){}\
                           }; register_rule<DOUT_##Y> register_DOUT_##Y

#define OUTPUT_VECTOR(X,Y) class OUT_##Y : public vector_node_output {\
                           public:\
                           OUT_##Y() : vector_node_output(X,#Y){}\
                           }; register_rule<OUT_##Y> register_OUT_##Y ;

#define OUTPUT_VECTOR_ALWAYS(X,Y) class OUT_##Y : public \
                           vector_node_output_always {\
                           public:\
                           OUT_##Y() : vector_node_output_always(X,#Y){}\
                           }; register_rule<OUT_##Y> register_OUT_##Y

#define OUTPUT_VECTORD(X,Y) class DOUT_##Y : public vectord_node_output {\
                           public:\
                           DOUT_##Y() : vectord_node_output(X,#Y){}\
                           }; register_rule<DOUT_##Y> register_DOUT_##Y ;

#define OUTPUT_VECTORD_ALWAYS(X,Y) class DOUT_##Y : public \
                           vectord_node_output_always {\
                           public:\
                           DOUT_##Y() : vectord_node_output_always(X,#Y){}\
                           }; register_rule<DOUT_##Y> register_DOUT_##Y

#ifndef DOUBLEPRECISIONOUTPUT
  OUTPUT_SCALAR_ALWAYS("cell2node(rho)",r) ;
  OUTPUT_SCALAR("cell2node(kineticEnergy)",kineticEnergy) ;
  OUTPUT_SCALAR_ALWAYS("cell2node(p)",pg) ;
  OUTPUT_SCALAR("cell2nodeMaxMag(pResidual)",pResidual) ;
  OUTPUT_SCALAR_ALWAYS("cell2node(temperature)",t) ;
  OUTPUT_SCALAR("cell2node(h)",h) ;
  OUTPUT_SCALAR("cell2nodeMaxMag(hResidual)",hResidual) ;
  OUTPUT_SCALAR("cell2node(k)",k) ;
  OUTPUT_SCALAR("cell2nodeMaxMag(kResidual)",kResidual) ;
  OUTPUT_SCALAR("cell2node(omega)",omega) ;
  OUTPUT_SCALAR("cell2nodeMaxMag(omegaResidual)",omegaResidual) ;
  OUTPUT_SCALAR("cell2node(laminarViscosity)",laminarViscosity) ;
  OUTPUT_SCALAR("cell2node(viscosityRatio)",viscosityRatio) ;
  OUTPUT_SCALAR("cell2node(thermalConductivity)",thermalConductivity) ;
  OUTPUT_SCALAR("cell2node(vort_mag)",vort_mag) ;
  OUTPUT_SCALAR("sMag",sMag) ;
  OUTPUT_SCALAR("cell2nodeMaxMag(timeStepFactor)",timeStepFactor) ;
  OUTPUT_SCALAR_ALWAYS_CONSTRAINT("cell2node(soundSpeed)",a,
    "compressibleFlow") ;
  OUTPUT_SCALAR("cell2node(lDES)",lDES) ;
  OUTPUT_VECTOR_ALWAYS("cell2node_v3d(v)",v) ;
  OUTPUT_VECTOR("cell2node_v3d(vort)",vort) ;
  OUTPUT_VECTOR("cell2nodeMaxMag_v3d(vResidual)",vResidual) ;
  OUTPUT_VECTOR("cell2nodeMaxMag_v3d(dbdSourceTerm)",dbdSourceTerm) ;
  OUTPUT_VECTOR("cell2node_v3d(vAbs)",vAbs) ;
#else
  OUTPUT_SCALARD_ALWAYS("dcell2node(p)",pg) ;
  OUTPUT_VECTORD_ALWAYS("dcell2node_v3d(v)",v) ;

  OUTPUT_SCALARD_ALWAYS("dcell2node(rho)",r) ;
  OUTPUT_SCALARD("dcell2node(kineticEnergy)",kineticEnergy) ;
  OUTPUT_SCALARD_ALWAYS("dcell2node(p)",pg) ;
  OUTPUT_SCALAR("cell2nodeMaxMag(pResidual)",pResidual) ;
  OUTPUT_SCALARD_ALWAYS("dcell2node(temperature)",t) ;
  OUTPUT_SCALARD("dcell2node(h)",h) ;
  OUTPUT_SCALAR("cell2nodeMaxMag(hResidual)",hResidual) ;
  OUTPUT_SCALARD("dcell2node(k)",k) ;
  OUTPUT_SCALAR("cell2nodeMaxMag(kResidual)",kResidual) ;
  OUTPUT_SCALARD("dcell2node(omega)",omega) ;
  OUTPUT_SCALAR("cell2nodeMaxMag(omegaResidual)",omegaResidual) ;
  OUTPUT_SCALARD("dcell2node(laminarViscosity)",laminarViscosity) ;
  OUTPUT_SCALARD("dcell2node(viscosityRatio)",viscosityRatio) ;
  OUTPUT_SCALARD("dcell2node(thermalConductivity)",thermalConductivity) ;
  OUTPUT_SCALARD("dcell2node(vort_mag)",vort_mag) ;
  OUTPUT_SCALARD("sMag",sMag) ;
  OUTPUT_SCALAR("cell2nodeMaxMag(timeStepFactor)",timeStepFactor) ;
  OUTPUT_SCALAR_ALWAYS_CONSTRAINT("cell2node(soundSpeed)",a,
    "compressibleFlow") ;
  OUTPUT_SCALAR("cell2node(lDES)",lDES) ;
  OUTPUT_VECTORD_ALWAYS("dcell2node_v3d(v)",v) ;
  OUTPUT_VECTORD("dcell2node_v3d(vort)",vort) ;
  OUTPUT_VECTOR("cell2nodeMaxMag_v3d(vResidual)",vResidual) ;
#endif

  // Copy of Ed's rule to dump the ambient pressure. Modifications have
  // been made.
  class DumpPAmbient : public singleton_rule {
    private:
      const_param<real> Pambient ;
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> it ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      DumpPAmbient() {
        name_store("Pambient{n,it}", Pambient) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("$it{n,it}",it) ;
        name_store("plot_modulo{n,it}",plot_modulo) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        name_store("modelName{n,it}",modelName) ;
#ifndef ITERATIONOUTPUT
        conditional("do_plot{n,it}") ;
#endif
        constraint("pos{n,it}") ;
        input("modelName{n,it},ncycle{n},plot_modulo{n,it}") ;
        input("Pambient{n,it},timeStepNumber{n},$it{n,it}") ;
        output("OUTPUT{n,it}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
#ifdef ITERATIONOUTPUT
        solver_dump_var(seq,Pambient.Rep(),ncyc,it,plot_modulo,modelName,
          "par","Pambient") ;
#else
        solver_dump_var(seq,Pambient.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "par","Pambient") ;
#endif
      }
  } ;

  register_rule<DumpPAmbient> registerDumpPAmbient ;

  // Copy of Ed's rule to dump the species mass fractions. Modifications have
  // been made.
  class DumpMixture : public pointwise_rule {
    private:
      const_param<int> numSpecies ;
      const_param<EOS> eos ;
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> it ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      const_storeVec<float> mixture ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      DumpMixture() {
        name_store("numSpecies{n,it}",numSpecies) ;
        name_store("eos{n,it}",eos) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("$it{n,it}",it) ;
        name_store("plot_modulo{n,it}",plot_modulo) ;
        name_store("modelName{n,it}",modelName) ;
        name_store("cell2node_v(y){n,it}",mixture) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("cell2node_v(y){n,it},modelName{n,it},ncycle{n},eos{n,it}") ;
        input("plot_modulo{n,it},timeStepNumber{n},numSpecies{n,it},$it{n,it}") ;
        output("OUTPUT{n,it}") ;
#ifndef ITERATIONOUTPUT
        conditional("do_plot{n,it}") ;
#endif
        constraint("pos{n,it},speciesTransport{n,it}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
        ostringstream oss ; int cycle=*ncycle ;
        if(*plot_modulo != 0) cycle = cycle % *plot_modulo ;
#ifndef ITERATIONOUTPUT
        oss << "output/mix." << cycle << "_" << *modelName ;
#else
        oss << "output/mix." << *it << "_" << *modelName ;
#endif
        string filename = oss.str() ;
        if(Loci::MPI_rank == 0) cout << "Writing file: " << filename << endl ;
        hid_t file_id=Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
          H5P_DEFAULT, H5P_DEFAULT) ;
        Loci::writeContainer(file_id,"mixture",mixture.Rep()) ;
        ostringstream species ;
        const vector<string>& name=eos->speciesNames() ;
        for(int i=0;i<*numSpecies;++i) species << name[i] << endl ;
        param<string> species_names ; *species_names=species.str() ;
        Loci::writeContainer(file_id,"species_names",species_names.Rep()) ;
        Loci::hdf5CloseFile(file_id) ;
      }
  } ;

  register_rule<DumpMixture> registerDumpMixture ;

  // Copy of Ed's rule to dump boundary geometry with modifications.
  class DumpBoundaryGeom : public pointwise_rule {
    private:
      const_store<Area> area ;
      const_store<vect3d> facecenter ;
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      DumpBoundaryGeom() {
        name_store("area{n,it}",area) ;
        name_store("facecenter{n,it}",facecenter) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n,it}",plot_modulo) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        name_store("modelName{n,it}",modelName) ;
        input("area{n,it},facecenter{n,it}") ;
        input("modelName{n,it}") ;
        input("ncycle{n},plot_modulo{n,it}") ;
        input("timeStepNumber{n}") ;
        constraint("ci->vol{n,it}") ;
        output("OUTPUT{n,it}") ;
        conditional("do_plot_boundary{n,it}") ;
      }

      // Write the boundary data.
      void compute(const sequence &seq) {
        ostringstream oss ; int cycle=*ncycle ;
        if(*plot_modulo != 0) cycle = cycle % *plot_modulo ;
        oss << "output/bc_geom." << cycle << "_" << *modelName ;
        string filename = oss.str() ;
        entitySet set(seq) ;
        hid_t file_id=createUnorderedFile(filename.c_str(),set) ;
        writeUnorderedStore(file_id,area,set,"area") ;
        writeUnorderedStore(file_id,facecenter,set,"facecenter") ;
        Loci::closeUnorderedFile(file_id) ;
      }

  } ;

  register_rule<DumpBoundaryGeom> registerDumpBoundaryGeom ;

  // Copy of Ed's rule to dump a boundary scalar.
  class DumpBoundaryScalar : public pointwise_rule {
    private:
      string var_name ;
      string value_name ;
      const_store<real> var ;
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      DumpBoundaryScalar(const char *vname,const char *valname) {
        var_name=string(vname)+"{n,it}" ; value_name=string(valname) ;
        name_store(var_name,var) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n,it}",plot_modulo) ;
        name_store("modelName{n,it}",modelName) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input(var_name) ;
        input("modelName{n,it}") ;
        input("ncycle{n},plot_modulo{n,it}") ;
        input("timeStepNumber{n}") ;
        constraint("ci->vol{n,it}") ;
        constraint(var_name) ;
        output("OUTPUT{n,it}") ;
        conditional("do_plot_boundary{n,it}") ;
      }

      // Write the boundary values.
      void compute(const sequence &seq) {
        ostringstream oss ; int cycle=*ncycle ;
        if(*plot_modulo != 0) cycle = cycle % *plot_modulo ;
        oss << "output/" + value_name + "_bnd." << cycle << "_" << *modelName ;
        string filename = oss.str() ;
        entitySet set(seq) ;
        hid_t file_id = createUnorderedFile(filename.c_str(),set) ;
        writeUnorderedStore(file_id,var,set,value_name.c_str()) ;
        Loci::closeUnorderedFile(file_id) ;
      }
  } ;

#define OUTPUT_BNDRY_SCALAR(X,Y,Z) class OUTB_##Y : \
public DumpBoundaryScalar { \
                           public:\
                               OUTB_##Y() : DumpBoundaryScalar(X,#Y){\
                                 constraint(Z);}          \
                           }; register_rule<OUTB_##Y> register_OUTB_##Y

  OUTPUT_BNDRY_SCALAR("qWall",qdot,"noslip_BC")  ;
  OUTPUT_BNDRY_SCALAR("yPlusWall",yplus,"noslip_BC")  ;
  OUTPUT_BNDRY_SCALAR("temperature_f",tw,"noslip_BC")  ;
  OUTPUT_BNDRY_SCALAR("p_f",pw,"noslip_BC")  ;

  // Copy of Ed's rule to dump a boundary vector.
  class DumpBoundaryVector : public pointwise_rule {
    private:
      string var_name ;
      string value_name ;
      const_store<vect3d> var ;
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      DumpBoundaryVector(const char *vname,const char *valname) {
        var_name = string(vname)+"{n,it}" ; value_name = string(valname) ;
        name_store(var_name,var) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n,it}",plot_modulo) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        name_store("modelName{n,it}",modelName) ;
        conditional("do_plot_boundary{n,it}") ;
        input(var_name) ;
        input("modelName{n,it}") ;
        input("ncycle{n},plot_modulo{n,it}") ;
        input("timeStepNumber{n}") ;
        constraint("ci->vol{n,it}") ;
        output("OUTPUT{n,it}") ;
      }

      void compute(const sequence &seq) {
        ostringstream oss ; int cycle=*ncycle ;
        if(*plot_modulo != 0) cycle = cycle % *plot_modulo ;
        oss << "output/" + value_name + "_bndvec." << cycle << "_"
          << *modelName ;
        string filename = oss.str() ;
        entitySet set(seq) ;
        hid_t file_id = createUnorderedFile(filename.c_str(),set) ;
        writeUnorderedStore(file_id,var,set,value_name.c_str()) ;
        Loci::closeUnorderedFile(file_id) ;
      }

  } ;

#define OUTPUT_BNDRY_VECTOR(X,Y,Z) class OUTB_##Y : \
public DumpBoundaryVector { \
                           public:\
                               OUTB_##Y() : DumpBoundaryVector(X,#Y){\
                                 constraint(Z);}                       \
                           }; register_rule<OUTB_##Y> register_OUTB_##Y

  OUTPUT_BNDRY_VECTOR("tauWall",tau,"noslip_BC")  ;


}
