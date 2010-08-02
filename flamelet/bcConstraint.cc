//-----------------------------------------------------------------------------
// Description: This file contains rules for creating the additional boundary
//   condition constraints that used to be created in readGrid.cc .
//-----------------------------------------------------------------------------

// Loci includes.
#include <Loci.h>
                                                                                
// StreamUns includes.
#include "sciTypes.h"
                                                                                
namespace streamUns {

  // Add additional constraints for fixed pressure outlets. 
  class FlameletFixedPressureOutletConstraints : public pointwise_rule {
    private:
      store<bool> extrapolatedZ_BC ;
      store<bool> extrapolatedZvar_BC ;
    public:

      // Define input and output.
      FlameletFixedPressureOutletConstraints() {
        name_store("extrapolatedZ_BC",extrapolatedZ_BC) ;
        name_store("extrapolatedZvar_BC",extrapolatedZvar_BC) ;
        output("extrapolatedZ_BC,extrapolatedZvar_BC") ;
        constraint("fixedPressureOutlet_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<FlameletFixedPressureOutletConstraints>
    registerFlameletFixedPressureOutletConstraints ;

 
  // Add additional constraints for no-slip boundaries. 
  class FlameletNoslipConstraints : public pointwise_rule {
    private:
      store<bool> extrapolatedZ_BC,extrapolatedZvar_BC ;
    public:

      // Define input and output.
      FlameletNoslipConstraints() {
        name_store("extrapolatedZ_BC",extrapolatedZ_BC) ;
        name_store("extrapolatedZvar_BC",extrapolatedZvar_BC) ;
        output("extrapolatedZ_BC,extrapolatedZvar_BC") ;
        constraint("noslip_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<FlameletNoslipConstraints> registerFlameletNoslipConstraints ;

  // Add additional constraints for slip boundaries. 
  class FlameletSlipConstraints : public pointwise_rule {
    private:
      store<bool> extrapolatedZ_BC,extrapolatedZvar_BC ;
    public:

      // Define input and output.
      FlameletSlipConstraints() {
        name_store("extrapolatedZ_BC",extrapolatedZ_BC) ;
        name_store("extrapolatedZvar_BC",extrapolatedZvar_BC) ;
        output("extrapolatedZ_BC,extrapolatedZvar_BC") ;
        constraint("slip_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<FlameletSlipConstraints> registerFlameletSlipConstraints ;


  // Add additional constraints for slip boundaries. 
  class FlameletSymmetryConstraints : public pointwise_rule {
    private:
      store<bool> extrapolatedZ_BC,extrapolatedZvar_BC ;
    public:

      // Define input and output.
      FlameletSymmetryConstraints() {
        name_store("extrapolatedZ_BC",extrapolatedZ_BC) ;
        name_store("extrapolatedZvar_BC",extrapolatedZvar_BC) ;
        output("extrapolatedZ_BC,extrapolatedZvar_BC") ;
        constraint("symmetry_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<FlameletSymmetryConstraints> registerFlameletSymmetryConstraints ;

}
