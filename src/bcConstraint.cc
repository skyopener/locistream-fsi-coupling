//-----------------------------------------------------------------------------
// Description: This file contains rules for creating the additional boundary
//   condition constraints that used to be created in readGrid.cc .
//-----------------------------------------------------------------------------

// Loci includes.
#include <Loci.h>
                                                                                
// StreamUns includes.
#include "sciTypes.h"
                                                                                
namespace streamUns {

  // Add additional constraints for incompressible inlets.
  class IncompressibleInletConstraints : public pointwise_rule {
    private:
      store<bool> inlet_BC,extrapolatedPressure_BC ;
    public:

      // Define input and output.
      IncompressibleInletConstraints() {
        name_store("inlet_BC",inlet_BC) ;
        name_store("extrapolatedPressure_BC",extrapolatedPressure_BC) ;
        output("inlet_BC,extrapolatedPressure_BC") ;
        constraint("incompressibleInlet_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<IncompressibleInletConstraints>
    registerIncompressibleInletConstraints ;

  // Add additional constraints for subsonic inlets.
  class SubsonicInletConstraints : public pointwise_rule {
    private:
      store<bool> inlet_BC,extrapolatedPressure_BC ;
    public:

      // Define input and output.
      SubsonicInletConstraints() {
        name_store("inlet_BC",inlet_BC) ;
        name_store("extrapolatedPressure_BC",extrapolatedPressure_BC) ;
        output("inlet_BC,extrapolatedPressure_BC") ;
        constraint("subsonicInlet_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<SubsonicInletConstraints> registerSubsonicInletConstraints ;

  // Add additional constraints for supersonic inlets.
  class SupersonicInletConstraints : public pointwise_rule {
    private:
      store<bool> inlet_BC,specifiedPressure_BC ;
    public:

      // Define input and output.
      SupersonicInletConstraints() {
        name_store("inlet_BC",inlet_BC) ;
        name_store("specifiedPressure_BC",specifiedPressure_BC) ;
        output("inlet_BC,specifiedPressure_BC") ;
        constraint("supersonicInlet_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<SupersonicInletConstraints> registerSupersonicInletConstraints ;

  // Add additional constraints for total pressure inlets.
  class TotalPressureInletConstraints : public pointwise_rule {
    private:
      store<bool> inlet_BC,extrapolatedPressure_BC ;
    public:

      // Define input and output.
      TotalPressureInletConstraints() {
        name_store("inlet_BC",inlet_BC) ;
        name_store("extrapolatedPressure_BC",extrapolatedPressure_BC) ;
        output("inlet_BC,extrapolatedPressure_BC") ;
        constraint("totalPressureInlet_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<TotalPressureInletConstraints>
    registerTotalPressureInletConstraints ;

  // Add additional constraints for extrapolated pressure outlets. 
  class ExtrapolatedPressureOutletConstraints : public pointwise_rule {
    private:
      store<bool> outlet_BC,nonSpecifiedMassFlux_BC ;
      store<bool> extrapolatedDensity_BC,extrapolatedPressure_BC ;
      store<bool> extrapolatedTemperature_BC,extrapolatedK_BC ;
      store<bool> extrapolatedOmega_BC,extrapolatedSpeciesMassFraction_BC ;
    public:

      // Define input and output.
      ExtrapolatedPressureOutletConstraints() {
        name_store("outlet_BC",outlet_BC) ;
        name_store("nonSpecifiedMassFlux_BC",nonSpecifiedMassFlux_BC) ;
        name_store("extrapolatedDensity_BC",extrapolatedDensity_BC) ;
        name_store("extrapolatedPressure_BC",extrapolatedPressure_BC) ;
        name_store("extrapolatedTemperature_BC",extrapolatedTemperature_BC) ;
        name_store("extrapolatedK_BC",extrapolatedK_BC) ;
        name_store("extrapolatedOmega_BC",extrapolatedOmega_BC) ;
        name_store("extrapolatedSpeciesMassFraction_BC",
          extrapolatedSpeciesMassFraction_BC) ;
        output("outlet_BC,nonSpecifiedMassFlux_BC,extrapolatedDensity_BC") ;
        output("extrapolatedPressure_BC,extrapolatedTemperature_BC") ;
        output("extrapolatedK_BC,extrapolatedOmega_BC") ;
        output("extrapolatedSpeciesMassFraction_BC") ;
        constraint("extrapolatedPressureOutlet_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<ExtrapolatedPressureOutletConstraints>
    registerExtrapolatedPressureOutletConstraints ;

  // Add additional constraints for fixed pressure outlets. 
  class FixedPressureOutletConstraints : public pointwise_rule {
    private:
      store<bool> outlet_BC,nonSpecifiedMassFlux_BC ;
      store<bool> extrapolatedDensity_BC,specifiedPressure_BC ;
      store<bool> extrapolatedTemperature_BC,extrapolatedK_BC ;
      store<bool> extrapolatedOmega_BC,extrapolatedSpeciesMassFraction_BC ;
    public:

      // Define input and output.
      FixedPressureOutletConstraints() {
        name_store("outlet_BC",outlet_BC) ;
        name_store("nonSpecifiedMassFlux_BC",nonSpecifiedMassFlux_BC) ;
        name_store("extrapolatedDensity_BC",extrapolatedDensity_BC) ;
        name_store("specifiedPressure_BC",specifiedPressure_BC) ;
        name_store("extrapolatedTemperature_BC",extrapolatedTemperature_BC) ;
        name_store("extrapolatedK_BC",extrapolatedK_BC) ;
        name_store("extrapolatedOmega_BC",extrapolatedOmega_BC) ;
        name_store("extrapolatedSpeciesMassFraction_BC",
          extrapolatedSpeciesMassFraction_BC) ;
        output("outlet_BC,nonSpecifiedMassFlux_BC,extrapolatedDensity_BC") ;
        output("specifiedPressure_BC,extrapolatedTemperature_BC") ;
        output("extrapolatedK_BC,extrapolatedOmega_BC") ;
        output("extrapolatedSpeciesMassFraction_BC") ;
        constraint("fixedPressureOutlet_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<FixedPressureOutletConstraints>
    registerFixedPressureOutletConstraints ;

  // Add additional constraints for mean pressure outlets. 
  class MeanPressureOutletConstraints : public pointwise_rule {
    private:
      store<bool> outlet_BC,nonSpecifiedMassFlux_BC ;
      store<bool> extrapolatedDensity_BC,specifiedPressure_BC ;
      store<bool> extrapolatedTemperature_BC,extrapolatedK_BC ;
      store<bool> extrapolatedOmega_BC,extrapolatedSpeciesMassFraction_BC ;
    public:

      // Define input and output.
      MeanPressureOutletConstraints() {
        name_store("outlet_BC",outlet_BC) ;
        name_store("nonSpecifiedMassFlux_BC",nonSpecifiedMassFlux_BC) ;
        name_store("extrapolatedDensity_BC",extrapolatedDensity_BC) ;
        name_store("specifiedPressure_BC",specifiedPressure_BC) ;
        name_store("extrapolatedTemperature_BC",extrapolatedTemperature_BC) ;
        name_store("extrapolatedK_BC",extrapolatedK_BC) ;
        name_store("extrapolatedOmega_BC",extrapolatedOmega_BC) ;
        name_store("extrapolatedSpeciesMassFraction_BC",
          extrapolatedSpeciesMassFraction_BC) ;
        output("outlet_BC,nonSpecifiedMassFlux_BC,extrapolatedDensity_BC") ;
        output("specifiedPressure_BC,extrapolatedTemperature_BC") ;
        output("extrapolatedK_BC,extrapolatedOmega_BC") ;
        output("extrapolatedSpeciesMassFraction_BC") ;
        constraint("meanPressureOutlet_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<MeanPressureOutletConstraints>
    registerMeanPressureOutletConstraints ;

  // Add additional constraints for no-slip boundaries. 
  class NoslipConstraints : public pointwise_rule {
    private:
      store<bool> extrapolatedDensity_BC,extrapolatedPressure_BC ;
      store<bool> extrapolatedSpeciesMassFraction_BC ;
    public:

      // Define input and output.
      NoslipConstraints() {
        name_store("extrapolatedDensity_BC",extrapolatedDensity_BC) ;
        name_store("extrapolatedPressure_BC",extrapolatedPressure_BC) ;
        name_store("extrapolatedSpeciesMassFraction_BC",
          extrapolatedSpeciesMassFraction_BC) ;
        output("extrapolatedDensity_BC,extrapolatedPressure_BC") ;
        output("extrapolatedSpeciesMassFraction_BC") ;
        constraint("noslip_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<NoslipConstraints> registerNoslipConstraints ;

  // Add additional constraints for slip boundaries. 
  class SlipConstraints : public pointwise_rule {
    private:
      store<bool> extrapolatedDensity_BC,extrapolatedPressure_BC ;
      store<bool> extrapolatedTemperature_BC,extrapolatedK_BC ;
      store<bool> extrapolatedOmega_BC,extrapolatedSpeciesMassFraction_BC ;
    public:

      // Define input and output.
      SlipConstraints() {
        name_store("extrapolatedDensity_BC",extrapolatedDensity_BC) ;
        name_store("extrapolatedPressure_BC",extrapolatedPressure_BC) ;
        name_store("extrapolatedTemperature_BC",extrapolatedTemperature_BC) ;
        name_store("extrapolatedK_BC",extrapolatedK_BC) ;
        name_store("extrapolatedOmega_BC",extrapolatedOmega_BC) ;
        name_store("extrapolatedSpeciesMassFraction_BC",
          extrapolatedSpeciesMassFraction_BC) ;
        output("extrapolatedDensity_BC,extrapolatedPressure_BC") ;
        output("extrapolatedTemperature_BC,extrapolatedK_BC") ;
        output("extrapolatedOmega_BC,extrapolatedSpeciesMassFraction_BC") ;
        constraint("slip_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<SlipConstraints> registerSlipConstraints ;

  // Add additional constraints for symmetry boundaries. 
  class SymmetryConstraints : public pointwise_rule {
    private:
      store<bool> extrapolatedDensity_BC,extrapolatedPressure_BC ;
      store<bool> extrapolatedTemperature_BC,extrapolatedK_BC ;
      store<bool> extrapolatedOmega_BC,extrapolatedSpeciesMassFraction_BC ;
    public:

      // Define input and output.
      SymmetryConstraints() {
        name_store("extrapolatedDensity_BC",extrapolatedDensity_BC) ;
        name_store("extrapolatedPressure_BC",extrapolatedPressure_BC) ;
        name_store("extrapolatedTemperature_BC",extrapolatedTemperature_BC) ;
        name_store("extrapolatedK_BC",extrapolatedK_BC) ;
        name_store("extrapolatedOmega_BC",extrapolatedOmega_BC) ;
        name_store("extrapolatedSpeciesMassFraction_BC",
          extrapolatedSpeciesMassFraction_BC) ;
        output("extrapolatedDensity_BC,extrapolatedPressure_BC") ;
        output("extrapolatedTemperature_BC,extrapolatedK_BC") ;
        output("extrapolatedOmega_BC,extrapolatedSpeciesMassFraction_BC") ;
        constraint("symmetry_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<SymmetryConstraints> registerSymmetryConstraints ;

  // Add additional constraints for interfaces.
  class InterfaceConstraints : public pointwise_rule {
    private:
      store<bool> specifiedPressure_BC ;
    public:

      // Define input and output.
      InterfaceConstraints() {
        name_store("specifiedPressure_BC",specifiedPressure_BC) ;
        output("specifiedPressure_BC") ;
        constraint("interface_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<InterfaceConstraints> registerInterfaceConstraints ;

}
