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
      const_param<real> timeStep ;
      param<real> newSolutionTime ;
    public:

      // Define input and output.
      AdvanceSolutionTime() {
        name_store("stime{n}",oldSolutionTime) ;
        name_store("timeStep{n}",timeStep) ;
        name_store("stime{n+1}",newSolutionTime) ;
        input("stime{n},timeStep{n}") ;
        output("stime{n+1}") ;
      }

      // Increment the solution time.
      virtual void compute(const sequence &seq) {
        *newSolutionTime=*oldSolutionTime+*timeStep ;
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
      const_param<int> plotFrequency ;
      param<bool> doPlot ;
    public:

      // Define input and output.
      DoPlot() {
        name_store("$n{n}",n) ;
        name_store("ncycle{n}",nCycle) ;
        name_store("plot_freq{n}",plotFrequency) ;
        name_store("do_plot{n}",doPlot) ;
        input("$n{n},ncycle{n},plot_freq{n}") ;
        output("do_plot{n}") ;
      }
 
      // Set the output flag.
      virtual void compute(const sequence &seq) {
        if((*n)==0){ doPlot=false ; return ; }
        doPlot=((*nCycle)%(*plotFrequency)==0) ;
      }
  } ;

  register_rule<DoPlot> registerDoPlot ;

  // Class to set the printing flag.
  class DoPrint : public singleton_rule {
    private:
      const_param<int> n,nCycle ;
      const_param<int> printFrequency ;
      param<bool> doPrint ;
    public:

      // Define input and output.
      DoPrint() {
        name_store("$n{n}",n) ;
        name_store("ncycle{n}",nCycle) ;
        name_store("print_freq{n}",printFrequency) ;
        name_store("do_print{n}",doPrint) ;
        input("$n{n},ncycle{n},print_freq{n}") ;
        output("do_print{n}") ;
      }
 
      // Set the output flag.
      virtual void compute(const sequence &seq) {
//doPrint=true ; return ;
        if((*n)==0){ doPrint=false ; return ; }
        doPrint=((*nCycle)%(*printFrequency)==0) ;
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
// Rules for the computing the CFL number for post-processing analysis.

  class ComputeCFLIncompressible : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_multiMap upper,lower,boundary_map ;
      const_param<real> timeStep ;
      const_store<vect3d> v,v_f ;
      const_store<Area> area ;
      const_store<real> vol ;
      store<real> CFL ;
    public:

      // Define input and output.
      ComputeCFLIncompressible() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("boundary_map",boundary_map) ;
        name_store("timeStep",timeStep) ;
        name_store("v",v) ;
        name_store("v_f",v_f) ;
        name_store("area",area) ;
        name_store("vol",vol) ;
        name_store("CFL",CFL) ;
        input("(lower,upper)->(cl,cr)->v,(lower,upper)->area") ;
        input("boundary_map->(area,v_f),timeStep,vol") ;
        output("CFL") ;
        constraint("incompressibleFlow,geom_cells") ;
      }

      // Compute CFL using method of Kim and Choi for unstructured meshes.
      void calculate(Entity cell) {
        CFL[cell]=0.0 ;
//cout << "CFLcalc: cell, v{n}: " << cell << " " << v[cell] << endl ;
        for(const int *ui=upper.begin(cell);ui!=upper.end(cell);++ui){
          const real volL=vol[cl[*ui]],volR=vol[cr[*ui]] ;
          const vect3d vFace=(1.0/(volL+volR))*(volL*v[cr[*ui]]+volR*
            v[cl[*ui]]) ;
          CFL[cell]+=abs(dot(vFace,area[*ui].n))*area[*ui].sada ;
        }
        for(const int *li=lower.begin(cell);li!=lower.end(cell);++li){
          const real volL=vol[cl[*li]],volR=vol[cr[*li]] ;
          const vect3d vFace=(1.0/(volL+volR))*(volL*v[cr[*li]]+volR*
            v[cl[*li]]) ;
          CFL[cell]+=abs(dot(vFace,area[*li].n))*area[*li].sada ;
        }
        for(const int *bi=boundary_map.begin(cell);bi!=boundary_map.end(cell);
        ++bi){
          CFL[cell]+=abs(dot(v_f[*bi],area[*bi].n))*area[*bi].sada ;
        }
        CFL[cell]*=0.5*(*timeStep)/vol[cell] ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeCFLIncompressible> registerComputeCFLIncompressible ;

  class ComputeCFLCompressible : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_multiMap upper,lower,boundary_map ;
      const_param<real> timeStep ;
      const_store<real> soundSpeed ;
      const_store<vect3d> v,v_f ;
      const_store<Area> area ;
      const_store<real> vol ;
      store<real> CFL ;
    public:

      // Define input and output.
      ComputeCFLCompressible() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("boundary_map",boundary_map) ;
        name_store("timeStep",timeStep) ;
        name_store("soundSpeed",soundSpeed) ;
        name_store("v",v) ;
        name_store("v_f",v_f) ;
        name_store("area",area) ;
        name_store("vol",vol) ;
        name_store("CFL",CFL) ;
        input("(lower,upper)->(cl,cr)->(v,soundSpeed)") ;
        input("(lower,upper)->area,boundary_map->(area,v_f),timeStep,vol") ;
        output("CFL") ;
        constraint("compressibleFlow,geom_cells") ;
      }

      // Compute CFL using method of Kim and Choi for unstructured meshes.
      void calculate(Entity cell) {
        CFL[cell]=0.0 ;
        for(const int *ui=upper.begin(cell);ui!=upper.end(cell);++ui){
          const real volL=vol[cl[*ui]],volR=vol[cr[*ui]] ;
          const vect3d vFace=(1.0/(volL+volR))*(volL*v[cr[*ui]]+volR*
            v[cl[*ui]]) ;
          const real aFace=(1.0/(volL+volR))*(volL*soundSpeed[cr[*ui]]+
            volR*soundSpeed[cl[*ui]]) ;
          CFL[cell]+=abs(dot(vFace,area[*ui].n)+aFace)*area[*ui].sada ;
        }
        for(const int *li=lower.begin(cell);li!=lower.end(cell);++li){
          const real volL=vol[cl[*li]],volR=vol[cr[*li]] ;
          const vect3d vFace=(1.0/(volL+volR))*(volL*v[cr[*li]]+volR*
            v[cl[*li]]) ;
          const real aFace=(1.0/(volL+volR))*(volL*soundSpeed[cr[*li]]+
            volR*soundSpeed[cl[*li]]) ;
          CFL[cell]+=abs(dot(vFace,area[*li].n)+aFace)*area[*li].sada ;
        }
        for(const int *bi=boundary_map.begin(cell);bi!=boundary_map.end(cell);
        ++bi){
          CFL[cell]+=abs(dot(v_f[*bi],area[*bi].n)+soundSpeed[cell])*
            area[*bi].sada ;
        }
        CFL[cell]*=0.5*(*timeStep)/vol[cell] ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeCFLCompressible> registerComputeCFLCompressible ;

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
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      scalar_node_output(const char *vname,const char *valname) {
        var_name = string(vname) ;
        value_name = string(valname) ;
        string var_name_time = var_name + "{n}" ;
        string constraint_name = string("scalarOutput_") + value_name ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n}",plot_modulo) ;
        name_store("modelName{n}",modelName) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        conditional("do_plot{n}") ;
        constraint("pos{n}") ;
        input("timeStepNumber{n},ncycle{n},plot_modulo{n},modelName{n}") ;
        input(var_name_time);
        constraint(constraint_name) ;
        output("OUTPUT{n}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
        solver_dump_var(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "sca",value_name) ;
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
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      scalar_node_output_always(const char *vname,const char *valname) {
        var_name = string(vname) ;
        value_name = string(valname) ;
        string var_name_time = var_name + "{n}"  ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n}",plot_modulo) ;
        name_store("modelName{n}",modelName) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        conditional("do_plot{n}") ;
        constraint("pos{n}") ;
        input("timeStepNumber{n},ncycle{n},plot_modulo{n},modelName{n}") ;
        input(var_name_time);
        output("OUTPUT{n}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
        solver_dump_var(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "sca",value_name) ;
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
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      scalard_node_output(const char *vname,const char *valname) {
        var_name = string(vname) ;
        value_name = string(valname) ;
        string var_name_time = var_name + "{n}" ;
        string constraint_name = string("scalarOutput_") + value_name ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n}",plot_modulo) ;
        name_store("modelName{n}",modelName) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        conditional("do_plot{n}") ;
        constraint("pos{n}") ;
        input("timeStepNumber{n},ncycle{n},plot_modulo{n},modelName{n}") ;
        input(var_name_time);
        constraint(constraint_name) ;
        output("OUTPUT{n}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
        solver_dump_var(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "sca",value_name) ;
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
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      scalard_node_output_always(const char *vname,const char *valname) {
        var_name = string(vname) ;
        value_name = string(valname) ;
        string var_name_time = var_name + "{n}"  ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n}",plot_modulo) ;
        name_store("modelName{n}",modelName) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        conditional("do_plot{n}") ;
        constraint("pos{n}") ;
        input("timeStepNumber{n},ncycle{n},plot_modulo{n},modelName{n}") ;
        input(var_name_time);
        output("OUTPUT{n}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
        solver_dump_var(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "sca",value_name) ;
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
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      vector_node_output(const char *vname,const char *valname) {
        var_name = string(vname) ;
        value_name = string(valname) ;
        string var_name_time = var_name + "{n}" ;
        string constraint_name = string("scalarOutput_") + value_name ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n}",plot_modulo) ;
        name_store("modelName{n}",modelName) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        conditional("do_plot{n}") ;
        constraint("pos{n}") ;
        input("timeStepNumber{n},ncycle{n},plot_modulo{n},modelName{n}") ;
        input(var_name_time);
        constraint(constraint_name) ;
        output("OUTPUT{n}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
        solver_dump_var(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "vec",value_name) ;
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
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      vector_node_output_always(const char *vname,const char *valname) {
        var_name = string(vname) ;
        value_name = string(valname) ;
        string var_name_time = var_name + "{n}" ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n}",plot_modulo) ;
        name_store("modelName{n}",modelName) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        conditional("do_plot{n}") ;
        constraint("pos{n}") ;
        input("timeStepNumber{n},ncycle{n}") ;
        input("plot_modulo{n},modelName{n}") ;
        input(var_name_time);
        output("OUTPUT{n}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
        solver_dump_var(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "vec",value_name) ;
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
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      vectord_node_output(const char *vname,const char *valname) {
        var_name = string(vname) ;
        value_name = string(valname) ;
        string var_name_time = var_name + "{n}" ;
        string constraint_name = string("scalarOutput_") + value_name ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n}",plot_modulo) ;
        name_store("modelName{n}",modelName) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        conditional("do_plot{n}") ;
        constraint("pos{n}") ;
        input("timeStepNumber{n},ncycle{n},plot_modulo{n},modelName{n}") ;
        input(var_name_time);
        constraint(constraint_name) ;
        output("OUTPUT{n}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
        solver_dump_var(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "vec",value_name) ;
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
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      vectord_node_output_always(const char *vname,const char *valname) {
        var_name = string(vname) ;
        value_name = string(valname) ;
        string var_name_time = var_name + "{n}" ;
        name_store(var_name_time, c2n) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n}",plot_modulo) ;
        name_store("modelName{n}",modelName) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        conditional("do_plot{n}") ;
        constraint("pos{n}") ;
        input("timeStepNumber{n},ncycle{n},plot_modulo{n},modelName{n}") ;
        input(var_name_time);
        output("OUTPUT{n}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
        solver_dump_var(seq,c2n.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "vec",value_name) ;
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
  OUTPUT_SCALAR("cell2nodeMaxMag(pPrimeStarResidual)",pPrimeStarResidual) ;
  OUTPUT_SCALAR_ALWAYS("cell2node(temperature)",t) ;
  OUTPUT_SCALAR("cell2node(h)",h) ;
  OUTPUT_SCALAR("cell2node(k)",k) ;
  OUTPUT_SCALAR("cell2node(omega)",omega) ;
  OUTPUT_SCALAR("cell2node(laminarViscosity)",laminarViscosity) ;
  OUTPUT_SCALAR("cell2node(viscosityRatio)",viscosityRatio) ;
  OUTPUT_SCALAR("cell2node(thermalConductivity)",thermalConductivity) ;
  OUTPUT_SCALAR("cell2node(vort_mag)",vort_mag) ;
  OUTPUT_SCALAR("sMag",sMag) ;
  OUTPUT_SCALAR("cell2nodeMaxMag(timeStepFactor)",timeStepFactor) ;
  OUTPUT_SCALAR_ALWAYS_CONSTRAINT("cell2node(soundSpeed)",a,
    "compressibleFlow") ;
  OUTPUT_SCALAR("cell2node(lDES)",lDES) ;
  OUTPUT_SCALAR("cell2nodeMaxMag(CFL)",CFL) ;
  OUTPUT_VECTOR_ALWAYS("cell2node_v3d(v)",v) ;
  OUTPUT_VECTOR("cell2node_v3d(vort)",vort) ;
#else
  OUTPUT_SCALARD_ALWAYS("dcell2node(p)",pg) ;
  OUTPUT_VECTORD_ALWAYS("dcell2node_v3d(v)",v) ;

  OUTPUT_SCALARD_ALWAYS("dcell2node(rho)",r) ;
  OUTPUT_SCALARD("dcell2node(kineticEnergy)",kineticEnergy) ;
  OUTPUT_SCALARD_ALWAYS("dcell2node(p)",pg) ;
  OUTPUT_SCALAR("cell2nodeMaxMag(pPrimeStarResidual)",pPrimeStarResidual) ;
  OUTPUT_SCALARD_ALWAYS("dcell2node(temperature)",t) ;
  OUTPUT_SCALARD("dcell2node(h)",h) ;
  OUTPUT_SCALARD("dcell2node(k)",k) ;
  OUTPUT_SCALARD("dcell2node(omega)",omega) ;
  OUTPUT_SCALARD("dcell2node(laminarViscosity)",laminarViscosity) ;
  OUTPUT_SCALARD("dcell2node(viscosityRatio)",viscosityRatio) ;
  OUTPUT_SCALARD("dcell2node(thermalConductivity)",thermalConductivity) ;
  OUTPUT_SCALARD("dcell2node(vort_mag)",vort_mag) ;
  OUTPUT_SCALARD("sMag",sMag) ;
  OUTPUT_SCALAR("cell2nodeMaxMag(timeStepFactor)",timeStepFactor) ;
  OUTPUT_SCALAR_ALWAYS_CONSTRAINT("cell2node(soundSpeed)",a,
    "compressibleFlow") ;
  OUTPUT_SCALAR("cell2node(lDES)",lDES) ;
  OUTPUT_SCALAR("cell2nodeMaxMag(CFL)",CFL) ;
  OUTPUT_VECTORD_ALWAYS("dcell2node_v3d(v)",v) ;
  OUTPUT_VECTORD("dcell2node_v3d(vort)",vort) ;
#endif

  // Copy of Ed's rule to dump the ambient pressure. Modifications have
  // been made.
  class DumpPAmbient : public singleton_rule {
    private:
      const_param<real> Pambient ;
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      DumpPAmbient() {
        name_store("Pambient{n}", Pambient) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n}",plot_modulo) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        name_store("modelName{n}",modelName) ;
        conditional("do_plot{n}") ;
        constraint("pos{n}") ;
        input("modelName{n},ncycle{n},plot_modulo{n}") ;
        input("Pambient{n},timeStepNumber{n}") ;
        output("OUTPUT{n}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
        solver_dump_var(seq,Pambient.Rep(),ncyc,ncycle,plot_modulo,modelName,
          "par","Pambient") ;
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
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      const_storeVec<float> mixture ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      DumpMixture() {
        name_store("numSpecies{n}",numSpecies) ;
        name_store("eos{n}",eos) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n}",plot_modulo) ;
        name_store("modelName{n}",modelName) ;
        name_store("cell2node_v(y){n}",mixture) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        input("cell2node_v(y){n},modelName{n},ncycle{n},eos{n}") ;
        input("plot_modulo{n},timeStepNumber{n},numSpecies{n}") ;
        output("OUTPUT{n}") ;
        conditional("do_plot{n}") ;
        constraint("pos{n},speciesTransport{n}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
        ostringstream oss ; int cycle=*ncycle ;
        if(*plot_modulo != 0) cycle = cycle % *plot_modulo ;
        oss << "output/mix." << cycle << "_" << *modelName ;
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
        name_store("area{n}",area) ;
        name_store("facecenter{n}",facecenter) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n}",plot_modulo) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        name_store("modelName{n}",modelName) ;
        input("area{n},facecenter{n}") ;
        input("modelName{n}") ;
        input("ncycle{n},plot_modulo{n}") ;
        input("timeStepNumber{n}") ;
        constraint("ci->vol{n}") ;
        output("OUTPUT{n}") ;
        conditional("do_plot{n}") ;
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
        var_name=string(vname)+"{n}" ; value_name=string(valname) ;
        name_store(var_name,var) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n}",plot_modulo) ;
        name_store("modelName{n}",modelName) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        input(var_name) ;
        input("modelName{n}") ;
        input("ncycle{n},plot_modulo{n}") ;
        input("timeStepNumber{n}") ;
        constraint("ci->vol{n}") ;
        constraint(var_name) ;
        output("OUTPUT{n}") ;
        conditional("do_plot{n}") ;
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
        var_name = string(vname)+"{n}" ; value_name = string(valname) ;
        name_store(var_name,var) ;
        name_store("timeStepNumber{n}",ncyc) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n}",plot_modulo) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        name_store("modelName{n}",modelName) ;
        conditional("do_plot{n}") ;
        input(var_name) ;
        input("modelName{n}") ;
        input("ncycle{n},plot_modulo{n}") ;
        input("timeStepNumber{n}") ;
        constraint("ci->vol{n}") ;
        output("OUTPUT{n}") ;
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
