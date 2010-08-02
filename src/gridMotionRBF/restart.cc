// Standard library includes.
#include <string>
#include <sstream>
#include <stdlib.h>
using std::string ;
                                                                                
// Loci includes.
#include <Loci.h>
                                                                                
// StreamUns includes.
#include "initialCondition.h"
#include "sciTypes.h"
                                                                                
namespace streamUns {

  // Class for writing vector restart data.
  class RestartVectorNode : public pointwise_rule {
    private:
      string variableName,constraintName ;
      const_param<int> timeStepNum ;
      const_store<vect3d> variable ;
      const_param<int> restart_modulo ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      RestartVectorNode(const char *variableName,const char *constraintName) :
        variableName(variableName),constraintName(constraintName) {
        string fullVariableName=variableName+string("{n,it}") ;
        string fullConstraintName=constraintName+string("{n,it}") ;
        name_store("ncycle{n}",timeStepNum) ;
        name_store(fullVariableName,variable) ;
        name_store("restart_modulo{n,it}",restart_modulo) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input(fullVariableName) ;
        input("ncycle{n},restart_modulo{n,it}") ;
        output("OUTPUT{n,it}") ;
        conditional("do_restart{n,it}") ;
        constraint("nodes{n,it}") ;
        constraint(fullConstraintName) ;
      }

      // Write the cell data to the restart file. Note the use of a temporary
      // store<vect3d> to ensure that we only write out data for geom_cells.
      // This is required since .Rep() method that writes out the data writes
      // out all the data including periodic cells.
      void compute(const sequence &seq) {
        unsigned int fileExtension=(*restart_modulo!=0)? (*timeStepNum)%
          (*restart_modulo):(*timeStepNum) ;
        ostringstream oss ; oss << "restart/" << variableName << "_hdf5."
          << fileExtension ;
        string fileName=oss.str() ;
        if(Loci::MPI_rank==0) cout << "Writing restart file: " << fileName
          << endl ;
        hid_t file_id=Loci::hdf5CreateFile(fileName.c_str(),H5F_ACC_TRUNC,
          H5P_DEFAULT, H5P_DEFAULT) ;
        store<vect3d> tempVariable ; tempVariable.allocate(entitySet(seq)) ;
        for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
          tempVariable[*si]=variable[*si] ;
        Loci::writeContainer(file_id,variableName,tempVariable.Rep()) ;
        Loci::hdf5CloseFile(file_id) ;
      }
  } ;

// Macro which is used to create a restart output class for each vector.
#define RESTART_VECTOR(X,Y) class Restart_##X : public RestartVectorNode {\
                              public:\
                                Restart_##X() : RestartVectorNode(#X,#Y){}\
                            } ;\
                            register_rule<Restart_##X> registerRestart_##X ;

  // Create the restart class for node displacement.
  RESTART_VECTOR(node_s,UNIVERSE) ;

}
