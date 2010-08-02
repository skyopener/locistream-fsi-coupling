#include <Loci.h>
#include "sciTypes.h"
#include <string>
#include <sstream>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using std::vector ;

//#define ITERATIONOUTPUT

namespace streamUns {

  // Outputs the grid topology for post-processing purposes.
  class OutputGridTopology : public pointwise_rule {
    private:
      const_param<int> n ;
      const_multiMap upper,lower,boundary_map,face2node ;
      const_store<vect3d> pos ;
      const_Map ref ;
      const_store<string> boundary_names ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      OutputGridTopology() {
        name_store("$n{n}",n) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("boundary_map",boundary_map) ;
        name_store("face2node",face2node) ;
        name_store("ref",ref) ;
        name_store("boundary_names{n}",boundary_names) ;
        name_store("modelName{n}",modelName) ;
        name_store("pos{n}",pos) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        input("$n{n},modelName{n}") ;
        input("(upper,lower,boundary_map)->face2node->pos{n}") ;
        input("boundary_map->ref->boundary_names{n}") ;
        output("OUTPUT{n}") ;
        disable_threading() ;
      }
    
      // Write out topology. Only write out on first time step.
      void compute(const sequence &seq){
        if(*n>0) return ;
        if(Loci::MPI_rank==0) cout << "Writing grid topology." << endl ;
        string filename="output/"+*modelName+".topo" ;
        Loci::parallelWriteGridTopology(filename.c_str(),upper.Rep(),
          lower.Rep(),boundary_map.Rep(),face2node.Rep(),ref.Rep(),
          boundary_names.Rep(),pos.Rep(),entitySet(seq)) ;
      }
  } ;

  register_rule<OutputGridTopology> registerOutputGridTopology ;

  // Function to write node positions.
  void dump_pos(const sequence &seq,const_store<vector3d<double> > &pos,
  const_param<int> &ncycle,const_param<int> &plot_modulo,const_param<string>
  &modelName) {
    ostringstream oss ; int cycle=*ncycle ;
    if(*plot_modulo != 0) cycle = cycle % *plot_modulo ;
    oss << "output/grid_pos." << cycle << "_" << *modelName ;
    string filename=oss.str() ;
    if(Loci::MPI_rank == 0) cout << "writing file " << filename << endl ;
    hid_t file_id=Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
      H5P_DEFAULT, H5P_DEFAULT) ;
    Loci::writeContainer(file_id,"pos",pos.Rep()) ;
    Loci::hdf5CloseFile(file_id) ;
  }

  // Writes out the position of the grid nodes. We used to write out
  // pos{n,it}, but now we write out pos{n}, since this is the proper
  // state of the node position at the beginning of a timestep, before
  // any positions have been advanced. For both time-dependent and
  // solution-dependent grid motion, pos{n,it} contains updated node
  // positions.
  class OutputGridPositions: public pointwise_rule {
    private:
      const_store<vector3d<double> > pos ;
      const_param<int> ncycle ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      OutputGridPositions() {
        name_store("pos{n}",pos) ;
        name_store("ncycle{n}",ncycle) ;
        name_store("plot_modulo{n}",plot_modulo) ;
        name_store("modelName{n}",modelName) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        input("pos{n},ncycle{n},plot_modulo{n},modelName{n}") ;
        output("OUTPUT{n}") ;
        conditional("do_plot{n}") ;
        constraint("pos{n}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
        dump_pos(seq,pos,ncycle,plot_modulo,modelName) ;
      }
  } ;

  register_rule<OutputGridPositions> registerOutputGridPositions ;

  class output_restart_grid: public pointwise_rule {
    const_store<vector3d<double> > pos ;
    const_param<int> ncycle ;
    const_param<int> restart_modulo ;
    const_param<string> modelName ;
    param<bool> OUTPUT ;
  public:
    output_restart_grid() {
      name_store("pos{n}",pos) ;
      name_store("ncycle{n}",ncycle) ;
      name_store("restart_modulo{n}",restart_modulo) ;
      name_store("modelName{n}",modelName) ;
      name_store("OUTPUT{n}",OUTPUT) ;
      
      conditional("do_restart{n}") ;
      constraint("pos{n}") ;
      input("pos{n}") ;
      input("ncycle{n},restart_modulo{n}") ;
      input("modelName{n}") ;
      output("OUTPUT{n}") ;
    }
    void compute(const sequence &seq) {
      ostringstream oss ;
      int cycle = *ncycle ;
      if(*restart_modulo != 0)
        cycle = cycle % *restart_modulo ;
      
      oss << "output/pos_restart." << cycle << "_" << *modelName ;
      string filename = oss.str() ;
      
      if(Loci::MPI_rank == 0)
        cout << "writing file " << filename << endl ;
      
      
      hid_t file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                           H5P_DEFAULT, H5P_DEFAULT) ;
      
      Loci::writeContainer(file_id,"pos",pos.Rep()) ;
      
      Loci::hdf5CloseFile(file_id) ;
    }
  } ;

  //register_rule<output_restart_grid> register_output_restart_grid ;


}

