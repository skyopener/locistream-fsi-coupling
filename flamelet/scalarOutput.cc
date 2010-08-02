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
//    blackbox<int> blackboxOutput ;
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




#define OUTPUT_SCALAR(X,Y) class OUT_##Y : public scalar_node_output {\
                           public:\
                           OUT_##Y() : scalar_node_output(X,#Y){}\
                           }; register_rule<OUT_##Y> register_OUT_##Y ;

  OUTPUT_SCALAR("cell2node(Z)",Z) ;
  OUTPUT_SCALAR("cell2node(Zvar)",Zvar) ;
  OUTPUT_SCALAR("cell2node(Chi)",Chi) ;

}
