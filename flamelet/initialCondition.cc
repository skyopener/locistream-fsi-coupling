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
#include "flameletInitialCondition.h"
#include "flamelet_table.h"
#include "name_var.h"
#include "sciTypes.h"

extern Flamelet_Table flamelet_table;

namespace streamUns {

  //---------------------------------------------------------------------------
  // Class InitialCondition.

  // Returns the serialized buffer size.
  int FlameletInitialCondition::BufferSize() const {
    return 4 ;
  }

  // Overridden virtual method called by fact database when a fact of type
  // InitialCondition is read.
  istream& FlameletInitialCondition::Input(istream &in) {
    FlameletInitialConditionOptions optionsList ;
    in >> optionsList ; Input(optionsList) ; return in ;
  }

  // Gets the values for the initial conditions from the specified options.
  void FlameletInitialCondition::Input(const options_list &optionsList) {

    if(optionsList.optionExists("Z")){
      if(optionsList.getOptionValueType("Z")==Loci::REAL) {
        optionsList.getOption("Z",Z) ; ZDefined=true ;
      }else{
        cerr << "Incorrect type for 'Z' in flameletInitialCondition." << endl ;
      }
    }

    if(optionsList.optionExists("Zvar")){
      if(optionsList.getOptionValueType("Zvar")==Loci::REAL) {
        optionsList.getOption("Zvar",Zvar) ; ZvarDefined=true ;
      }else{
        cerr << "Incorrect type for 'Zvar' in flameletInitialCondition." << endl ;
      }
    }

  }

  // Packs the data into a buffer.
  void FlameletInitialCondition::PackBuffer(real *buffer,int size) {
    int i=0 ;
    buffer[i++]=Z ; buffer[i++]=Zvar ; 
    buffer[i++]=real(ZDefined) ; buffer[i++]=real(ZvarDefined) ;
  }

  // Prints the initial condition.
  ostream& FlameletInitialCondition::Print(ostream &out) const {
    int outputCount=0 ;
    out << "<" ;
    if(ZDefined){ out << "Z=" << Z ; ++outputCount ; }
	if(ZvarDefined){ out << "Zvar=" << Zvar ; ++outputCount ; }
    out << ">" ; return out ;
  }

  // Unpacks the data from a buffer.
  void FlameletInitialCondition::UnpackBuffer(real *buffer,int size) {
    int i=0 ;
    Z=buffer[i++] ; Zvar=buffer[i++] ; 
    ZDefined=bool(buffer[i++]) ; ZvarDefined=bool(buffer[i++]) ;
  }

  //---------------------------------------------------------------------------
  // Rules.

    class DensityInitialConditionFlamelet : public pointwise_rule {
    private:
          const_store<real> Z_ic,Zvar_ic,omega_ic ;
          store<real> rho_ic ;
    public:

      // Define input and output.
      DensityInitialConditionFlamelet() {
        name_store("Z_ic",Z_ic) ;
        name_store("Zvar_ic",Zvar_ic) ;
        name_store("omega_ic",omega_ic) ;
        name_store("flamelet::rho_ic",rho_ic) ;
        input("Z_ic,Zvar_ic,omega_ic") ;
        output("flamelet::rho_ic") ;
        constraint("geom_cells,flameletModel") ;
      }

      void calculate(Entity cell) {
            double Chi=log10(2.0*0.09*omega_ic[cell]*Zvar_ic[cell]+1.e-30);
            rho_ic[cell]=flamelet_table.get_rho(Z_ic[cell],Zvar_ic[cell],Chi,true);
      }

      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DensityInitialConditionFlamelet>
    registerDensityInitialConditionFlamelet ;

  // Rule to assign the initial condition for Z.
  class ZInitialConditionFlamelet : public pointwise_rule {
    private:
      const_param<FlameletInitialCondition> flameletInitialCondition ;
      store<real> Z_ic ;
    public:

      // Define input and output.
      ZInitialConditionFlamelet() {
        name_store("flameletInitialCondition",flameletInitialCondition) ;
        name_store("Z_ic",Z_ic) ;
        input("flameletInitialCondition") ;
        output("Z_ic") ;
        constraint("geom_cells,flameletModel") ;
      }

      // Assign Z for a single cell.
      void calculate(Entity cell) {
        Z_ic[cell]=flameletInitialCondition->z() ;
      }

      // Assign density for a sequence of cells.
      void compute(const sequence &seq) {

        // Check that Z was provided by the user.
        if(!flameletInitialCondition->IsZDefined())
          cerr << "ERROR: Initial condition for Z required!" << endl ;

        // Loop through the sequence of cells.
        do_loop(seq,this) ;
      }
  } ;

  register_rule<ZInitialConditionFlamelet>
    registerZInitialConditionFlamelet ;
  
  // Rule to assign the initial condition for Zvar.
  class ZvarInitialConditionFlamelet : public pointwise_rule {
    private:
      const_param<FlameletInitialCondition> flameletInitialCondition ;
      store<real> Zvar_ic ;
    public:

      // Define input and output.
      ZvarInitialConditionFlamelet() {
        name_store("flameletInitialCondition",flameletInitialCondition) ;
        name_store("Zvar_ic",Zvar_ic) ;
        input("flameletInitialCondition") ;
        output("Zvar_ic") ;
        constraint("geom_cells,flameletModel") ;
      }

      // Assign Zvar for a single cell.
      void calculate(Entity cell) {
        Zvar_ic[cell]=flameletInitialCondition->zvar() ;
      }

      // Assign density for a sequence of cells.
      void compute(const sequence &seq) {

        // Check that Zvar was provided by the user.
        if(!flameletInitialCondition->IsZvarDefined())
          cerr << "ERROR: Initial condition for Zvar required!" << endl ;

        // Loop through the sequence of cells.
        do_loop(seq,this) ;
      }
  } ;

  register_rule<ZvarInitialConditionFlamelet>
    registerZvarInitialConditionFlamelet ;

}




