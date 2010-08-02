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
  class CavitationFixedPressureOutletConstraints : public pointwise_rule {
    private:
      store<bool> extrapolatedAlpha_BC ;
    public:

      // Define input and output.
      CavitationFixedPressureOutletConstraints() {
        name_store("extrapolatedAlpha_BC",extrapolatedAlpha_BC) ;
        output("extrapolatedAlpha_BC") ;
        constraint("fixedPressureOutlet_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<CavitationFixedPressureOutletConstraints>
    registerCavitationFixedPressureOutletConstraints ;

 
  // Add additional constraints for no-slip boundaries. 
  class CavitationNoslipConstraints : public pointwise_rule {
    private:
      store<bool> extrapolatedAlpha_BC ;
    public:

      // Define input and output.
      CavitationNoslipConstraints() {
        name_store("extrapolatedAlpha_BC",extrapolatedAlpha_BC) ;
        output("extrapolatedAlpha_BC") ;
        constraint("noslip_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<CavitationNoslipConstraints> registerCavitationNoslipConstraints ;

  // Add additional constraints for slip boundaries. 
  class CavitationSlipConstraints : public pointwise_rule {
    private:
      store<bool> extrapolatedAlpha_BC ;
    public:

      // Define input and output.
      CavitationSlipConstraints() {
        name_store("extrapolatedAlpha_BC",extrapolatedAlpha_BC) ;
        output("extrapolatedAlpha_BC") ;
        constraint("slip_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<CavitationSlipConstraints> registerCavitationSlipConstraints ;


  // Add additional constraints for symmetry boundaries. 
  class CavitationSymmetryConstraints : public pointwise_rule {
    private:
      store<bool> extrapolatedAlpha_BC ;
    public:

      // Define input and output.
      CavitationSymmetryConstraints() {
        name_store("extrapolatedAlpha_BC",extrapolatedAlpha_BC) ;
        output("extrapolatedAlpha_BC") ;
        constraint("symmetry_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<CavitationSymmetryConstraints> registerCavitationSymmetryConstraints ;

}
