#include "sciTypes.h"
#include "qvi.h"
#include "transport_db.h"
#include "name_var.h"

#include <string>
#include <sstream>
#include <vector>
using std::vector ;
using std::string ;

namespace streamUns {
  using namespace fluidPhysics ;
  
  class readtran : public singleton_rule {
    const_param<name_var> cmodel ;
    const_param<std::string> transport_model ;
    const_param<std::string> diffusion_model ;
    const_param<conservativeVectorInfo> qvi ;
    param<transport_db> tran ;
    void read_transport(istream &input) ;
  public:
    readtran() ;
    virtual void compute( const sequence &seq) ;
  } ;

  readtran::readtran() {
    name_store("chemistry_model",cmodel) ;
    name_store("transport_model",transport_model) ;
    name_store("diffusion_model",diffusion_model) ;
    name_store("qvi",qvi) ;
    name_store("tran",tran) ;
    input("chemistry_model") ;
    input("transport_model") ;
    input("diffusion_model") ;
    input("qvi") ;
    output("tran") ;
  }

  void readtran::read_transport(istream &input) {
    (*tran).Input(input) ;
//    (*tran).Print(cout) ;
  } 
      
  void readtran::compute(const sequence &seq) {
    const char *base = getenv("CHEMISTRY_DATABASE") ;
    if(base == NULL) 
      base = "./" ;
    
    string cmdl = (*cmodel).name ;
    string tmdl = (*transport_model) ;
    string dmdl = (*diffusion_model) ;
    if(tmdl == "chemkin" || dmdl == "chemkin" ) {
      string dbname = base ;
      dbname += "/data_base/chemkin/" ;
      string fname = cmdl + ".tran" ;
      ifstream tranin(fname.c_str(),ios::in) ;
      if(tranin.fail()) {
        tranin.clear() ;
        fname = dbname+cmdl + ".tran" ;
        tranin.open(fname.c_str()) ;
        if(tranin.fail()) {
          cerr << "can't open tran.in, and no model specified." << endl ;
          exit(-1) ;
        }
      }
      if(Loci::MPI_rank == 0)
        cout << "The local file of chemkin transport data is " << fname << endl ;
      read_transport(tranin) ;
    } else {
      cerr << "Transport model " << tmdl << " not implemented or" << endl ;
      cerr << "Diffusion model " << dmdl << " not implemented yet." << endl ;
      exit(-1) ;
    }
    vector<string> species_list ;
    int ns = qvi->numSpecies() ;
    for(int i=0;i<ns;++i)
      species_list.push_back(qvi->speciesName(i)) ;
        
    (*tran).species_reorder(species_list) ;
  }

  register_rule<readtran> register_readtran ;
}
