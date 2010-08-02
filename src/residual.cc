//-----------------------------------------------------------------------------
// Description: This file contains rules for writing residual information.
//
// Authors: Siddharth Thakur and Jeff Wright
//          Streamline Numerics Inc.
//          Gainesville, Florida 32609
//-----------------------------------------------------------------------------
                                                                                
// Loci includes.
#include <Loci.h>
                                                                                
// StreamUns includes.
#include "residual.h"
#include "sciTypes.h"
#include "varsFileInputs.h"
                                                                                
namespace streamUns {

  // Class to extract reference scales for residual normalization from the
  // reference values options list.
  class ReferenceScales : public singleton_rule {
    private:
      const_param<ReferenceValue> referenceValue ;
      param<real> lScale,rhoScale,vScale,hScale,kScale,omegaScale ;
    public:

      // Define input and output.
      ReferenceScales() {
        name_store("referenceValue",referenceValue) ;
        name_store("lScale",lScale) ;
        name_store("rhoScale",rhoScale) ;
        name_store("vScale",vScale) ;
        name_store("hScale",hScale) ;
        name_store("kScale",kScale) ;
        name_store("omegaScale",omegaScale) ;
        input("referenceValue") ;
        output("lScale,rhoScale,vScale,hScale,kScale,omegaScale") ;
      }

      void compute(const sequence &seq) {

        // Set default values.
        *lScale=*rhoScale=*vScale=*hScale=*kScale=*omegaScale=1.0 ;

        // Set values for options that have been specified.
        if((*referenceValue).optionExists("L"))
          (*referenceValue).getOption("L",*lScale) ;
        if((*referenceValue).optionExists("rho"))
          (*referenceValue).getOption("rho",*rhoScale) ;
        if((*referenceValue).optionExists("v"))
          (*referenceValue).getOption("v",*vScale) ;
        if((*referenceValue).optionExists("h"))
          (*referenceValue).getOption("h",*hScale) ;
        if((*referenceValue).optionExists("k"))
          (*referenceValue).getOption("k",*kScale) ;
        if((*referenceValue).optionExists("omega"))
          (*referenceValue).getOption("omega",*omegaScale) ;
      }
  } ;

  register_rule<ReferenceScales> registerReferenceScales ;

  // Class to write maximum momentum residual location to debug file.
  class WriteMaximumMomentumResidualLocation : public singleton_rule {
    private:
      const_param<int> n,it ;
      const_param<VectorResidual> vResidualData ;
      const_param<bool> debug ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      WriteMaximumMomentumResidualLocation() {
        name_store("$n{n}",n) ;
        name_store("$it{n,it}",it) ;
        name_store("vResidualData{n,it}",vResidualData) ;
        name_store("debug{n,it}",debug) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("$n{n},$it{n,it},vResidualData{n,it},debug{n,it}") ;
        output("OUTPUT{n,it}") ;
      } ;

      // Write the location.
      void compute(const sequence &seq) {
        if(*debug) Loci::debugout << "n: " << *n << " it: " << *it
          << " vLoc,vMax: " << vResidualData->maxResidualLocation
          << ", " << vResidualData->maxResidual << endl ;
      }
  } ;

  register_rule<WriteMaximumMomentumResidualLocation>
    registerWriteMaximumMomentumResidualLocation ;

  // Class to write maximum pressure correction residual location to debug file.
  class WriteMaximumPressureCorrectionResidualLocation : public singleton_rule {
    private:
      const_param<int> n,it ;
      const_param<ScalarResidual> pPrimeResidualData ;
      const_param<bool> debug ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      WriteMaximumPressureCorrectionResidualLocation() {
        name_store("$n{n}",n) ;
        name_store("$it{n,it}",it) ;
        name_store("pPrimeResidualData{n,it}",pPrimeResidualData) ;
        name_store("debug{n,it}",debug) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("$n{n},$it{n,it},pPrimeResidualData{n,it},debug{n,it}") ;
        output("OUTPUT{n,it}") ;
      } ;

      // Write the location.
      void compute(const sequence &seq) {
        if(*debug) Loci::debugout << "n: " << *n << " it: " << *it
          << " pLoc,pMax: " << pPrimeResidualData->maxResidualLocation
          << ", " << pPrimeResidualData->maxResidual << endl ;
      }
  } ;

  register_rule<WriteMaximumPressureCorrectionResidualLocation>
    registerWriteMaximumPressureCorrectionResidualLocation ;

  // Class to write maximum energy residual location to debug file.
  class WriteMaximumEnergyResidualLocation : public singleton_rule {
    private:
      const_param<int> n,it ;
      const_param<ScalarResidual> hResidualData ;
      const_param<bool> debug ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      WriteMaximumEnergyResidualLocation() {
        name_store("$n{n}",n) ;
        name_store("$it{n,it}",it) ;
        name_store("hResidualData{n,it}",hResidualData) ;
        name_store("debug{n,it}",debug) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("$n{n},$it{n,it},hResidualData{n,it},debug{n,it}") ;
        output("OUTPUT{n,it}") ;
      } ;

      // Write the location.
      void compute(const sequence &seq) {
        if(*debug) Loci::debugout << "n: " << *n << " it: " << *it
          << " hLoc,hMax: " << hResidualData->maxResidualLocation
          << ", " << hResidualData->maxResidual << endl ;
      }
  } ;

  register_rule<WriteMaximumEnergyResidualLocation>
    registerWriteMaximumEnergyResidualLocation ;

  // Class to write maximum k residual location to debug file.
  class WriteMaximumKResidualLocation : public singleton_rule {
    private:
      const_param<int> n,it ;
      const_param<ScalarResidual> kResidualData ;
      const_param<bool> debug ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      WriteMaximumKResidualLocation() {
        name_store("$n{n}",n) ;
        name_store("$it{n,it}",it) ;
        name_store("kResidualData{n,it}",kResidualData) ;
        name_store("debug{n,it}",debug) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("$n{n},$it{n,it},kResidualData{n,it},debug{n,it}") ;
        output("OUTPUT{n,it}") ;
      } ;

      // Write the location.
      void compute(const sequence &seq) {
        if(*debug) Loci::debugout << "n: " << *n << " it: " << *it
          << " kLoc,kMax: " << kResidualData->maxResidualLocation
          << ", " << kResidualData->maxResidual << endl ;
      }
  } ;

  register_rule<WriteMaximumKResidualLocation>
    registerWriteMaximumKResidualLocation ;

  // Class to write maximum omega residual location to debug file.
  class WriteMaximumOmegaResidualLocation : public singleton_rule {
    private:
      const_param<int> n,it ;
      const_param<ScalarResidual> omegaResidualData ;
      const_param<bool> debug ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      WriteMaximumOmegaResidualLocation() {
        name_store("$n{n}",n) ;
        name_store("$it{n,it}",it) ;
        name_store("omegaResidualData{n,it}",omegaResidualData) ;
        name_store("debug{n,it}",debug) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("$n{n},$it{n,it},omegaResidualData{n,it},debug{n,it}") ;
        output("OUTPUT{n,it}") ;
      } ;

      // Write the location.
      void compute(const sequence &seq) {
        if(*debug) Loci::debugout << "n: " << *n << " it: " << *it
          << " omegaLoc,omegaMax: " << omegaResidualData->maxResidualLocation
          << ", " << omegaResidualData->maxResidual << endl ;
      }
  } ;

  register_rule<WriteMaximumOmegaResidualLocation>
    registerWriteMaximumOmegaResidualLocation ;
}
