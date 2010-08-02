// System includes.
#include <rpc/rpc.h>
#include <rpc/xdr.h>

// Loci includes.
#include <fact_db.h>
#include <Tools/parse.h>
using namespace Loci::parse ;

// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;

// StreamUns includes.
#include "cavitationInitialCondition.h"
//#include "cavitation_table.h"
#include "name_var.h"
#include "sciTypes.h"

//extern Cavitation_Table cavitation_table;

namespace streamUns {

  //---------------------------------------------------------------------------
  // Class InitialCondition.

  // Returns the serialized buffer size.
  int CavitationInitialCondition::BufferSize() const {
    return 4 ;
  }

  // Overridden virtual method called by fact database when a fact of type
  // InitialCondition is read.
  istream& CavitationInitialCondition::Input(istream &in) {
    CavitationInitialConditionOptions optionsList ;
    in >> optionsList ; Input(optionsList) ; return in ;
  }

  // Gets the values for the initial conditions from the specified options.
  void CavitationInitialCondition::Input(const options_list &optionsList) {

    if(optionsList.optionExists("Alpha")){
      if(optionsList.getOptionValueType("Alpha")==Loci::REAL) {
        optionsList.getOption("Alpha",Alpha) ; AlphaDefined=true ;
      }else{
        cerr << "Incorrect type for 'Alpha' in cavitationInitialCondition." << endl ;
      }
    }
  }

  // Packs the data into a buffer.
  void CavitationInitialCondition::PackBuffer(real *buffer,int size) {
    int i=0 ;
    buffer[i++]=Alpha ; 
    buffer[i++]=real(AlphaDefined) ;
  }

  // Prints the initial condition.
  ostream& CavitationInitialCondition::Print(ostream &out) const {
    int outputCount=0 ;
    out << "<" ;
    if(AlphaDefined){ out << "Alpha=" << Alpha ; ++outputCount ; }
    out << ">" ; return out ;
  }

  // Unpacks the data from a buffer.
  void CavitationInitialCondition::UnpackBuffer(real *buffer,int size) {
    int i=0 ;
    Alpha=buffer[i++] ; 
    AlphaDefined=bool(buffer[i++]) ;
  }

  //---------------------------------------------------------------------------
  // Rules.

  // Rule to assign the initial condition for Alpha.
  class AlphaInitialConditionCavitation : public pointwise_rule {
    private:
      const_param<CavitationInitialCondition> cavitationInitialCondition ;
      store<real> Alpha_ic ;
    public:

      // Define input and output.
      AlphaInitialConditionCavitation() {
        name_store("cavitationInitialCondition",cavitationInitialCondition) ;
        name_store("Alpha_ic",Alpha_ic) ;
        input("cavitationInitialCondition") ;
        output("Alpha_ic") ;
        constraint("geom_cells,cavitationModel") ;
      }

      // Assign Alpha for a single cell.
      void calculate(Entity cell) {
        Alpha_ic[cell]=cavitationInitialCondition->alpha() ;
      }

      // Assign alpha for a sequence of cells.
      void compute(const sequence &seq) {

        // Check that Alpha was provided by the user.
        if(!cavitationInitialCondition->IsAlphaDefined())
          cerr << "ERROR: Initial condition for Alpha required!" << endl ;

        // Loop through the sequence of cells.
        do_loop(seq,this) ;
      }
  } ;

  register_rule<AlphaInitialConditionCavitation>
    registerAlphaInitialConditionCavitation ;
  
}




