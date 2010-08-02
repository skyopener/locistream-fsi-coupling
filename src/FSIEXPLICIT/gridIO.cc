#include <Loci.h>
#include "sciTypes.h"
#include <string>
#include <sstream>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using std::vector ;

namespace streamUns {

  // Outputs the grid topology for post-processing purposes.
  class OutputGridTopologyGridOnly : public pointwise_rule {
    private:
      const_param<int> n ;
      const_param<int> it ;
      const_param<int> ncyc ;
      const_multiMap upper,lower,boundary_map,face2node ;
      const_store<vect3d> pos ;
      const_Map ref ;
      const_store<string> boundary_names ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      OutputGridTopologyGridOnly() {
        name_store("$nn{nn}",n) ;
        name_store("$itg{nn,itg}",it) ;
        name_store("upper",upper) ;
        name_store("lower",lower) ;
        name_store("boundary_map",boundary_map) ;
        name_store("face2node",face2node) ;
        name_store("ref",ref) ;
        name_store("boundary_names{nn,itg}",boundary_names) ;
        name_store("timeStepNumber{nn}",ncyc) ;
        name_store("modelName{nn,itg}",modelName) ;
        name_store("pos{nn,itg}",pos) ;
        name_store("OUTPUT{nn,itg}",OUTPUT) ;
        input("$nn{nn},$itg{nn,itg},timeStepNumber{nn},modelName{nn,itg}") ;
        input("(upper,lower,boundary_map)->face2node->pos{nn,itg}") ;
        input("boundary_map->ref->boundary_names{nn,itg}") ;
//      conditional("do_plot{nn,it}") ;
        output("OUTPUT{nn,itg}") ;
        disable_threading() ;
      }
    
      // Write out topology. Only write out on first time step.
      void compute(const sequence &seq){
//      if(*ncyc > 0) return ;
        if(*n>0) return ; if(*it>0) return ;
        if(Loci::MPI_rank==0) cout << "Writing grid topology." << endl ;
        string filename="output/"+*modelName+".topo" ;
        Loci::parallelWriteGridTopology(filename.c_str(),upper.Rep(),
          lower.Rep(),boundary_map.Rep(),face2node.Rep(),ref.Rep(),
          boundary_names.Rep(),pos.Rep(),entitySet(seq)) ;
      }
  } ;

  register_rule<OutputGridTopologyGridOnly>
    registerOutputGridTopologyGridOnly ;

  // Function to write node positions.
  void dump_pos_local(const sequence &seq,const_store<vector3d<double> > &pos,
  const_param<int> &ncycle,const_param<int> &plot_modulo,const_param<string>
  &modelName) {
    if(*ncycle<=-1) return ;
    ostringstream oss ; int cycle=*ncycle+1 ;
    if(*plot_modulo != 0) cycle = cycle % *plot_modulo ;
    oss << "output/grid_pos." << cycle << "_" << *modelName ;
    string filename=oss.str() ;
    if(Loci::MPI_rank == 0) cout << "writing file " << filename << endl ;
    hid_t file_id=Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
      H5P_DEFAULT, H5P_DEFAULT) ;
    Loci::writeContainer(file_id,"pos",pos.Rep()) ;
    Loci::hdf5CloseFile(file_id) ;
  }

  // Writes out the position of the grid nodes.
  class OutputGridPositionsGridOnly : public pointwise_rule {
    private:
      const_store<vector3d<double> > pos ;
      const_param<int> ncycle ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      OutputGridPositionsGridOnly() {
        name_store("pos{nn,itg}",pos) ;
        name_store("ncycleGridOnly{nn}",ncycle) ;
        name_store("plot_modulo{nn,itg}",plot_modulo) ;
        name_store("modelName{nn,itg}",modelName) ;
        name_store("OUTPUT{nn,itg}",OUTPUT) ;
        input("pos{nn,itg},ncycleGridOnly{nn}") ;
        input("plot_modulo{nn,itg},modelName{nn,itg}") ;
        output("OUTPUT{nn,itg}") ;
        conditional("do_plot{nn,itg}") ;
        constraint("pos{nn,itg}") ;
      }

      // Write to file.
      void compute(const sequence &seq) {
        dump_pos_local(seq,pos,ncycle,plot_modulo,modelName) ;
      }
  } ;

  register_rule<OutputGridPositionsGridOnly>
    registerOutputGridPositionsGridOnly ;
}

