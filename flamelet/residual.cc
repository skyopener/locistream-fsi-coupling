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

  // Class to write Z and Zvar residuals to screen.
  class PrintZZvarResidual : public singleton_rule {
    private:
      const_param<int> n,nCycle,it ;
      const_param<ScalarResidual> ZResidualData, ZvarResidualData ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      PrintZZvarResidual() {
        name_store("$n{n}",n) ;
	name_store("ncycle{n}",nCycle) ;
        name_store("$it{n,it}",it) ;
        name_store("ZResidualData{n,it}",ZResidualData) ;
        name_store("ZvarResidualData{n,it}",ZvarResidualData) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("$n{n},ncycle{n},$it{n,it},ZResidualData{n,it},ZvarResidualData{n,it}") ;
        output("OUTPUT{n,it}") ;
      } ;

      void compute(const sequence &seq) {

        // Set the output format.
        if(Loci::MPI_rank==0){
          cout.setf(ios::scientific,ios::floatfield) ; cout.precision(6) ;
        }  

        // Write out the residuals.
        if(Loci::MPI_rank==0){
	  cout <<"RZ: " << *nCycle << " " << *it << " " << ZResidualData->
            totalResidual << " " << ZvarResidualData->totalResidual << endl;
        }

      }
  } ;

  register_rule<PrintZZvarResidual>
    registerPrintZZvarResidual ;

  // Class to write maximum Z residual location to debug file.
  class WriteMaximumZResidualLocation : public singleton_rule {
    private:
      const_param<int> n,it ;
      const_param<ScalarResidual> ZResidualData ;
      const_param<bool> debug ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      WriteMaximumZResidualLocation() {
        name_store("$n{n}",n) ;
        name_store("$it{n,it}",it) ;
        name_store("ZResidualData{n,it}",ZResidualData) ;
        name_store("debug{n,it}",debug) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("$n{n},$it{n,it},ZResidualData{n,it},debug{n,it}") ;
        output("OUTPUT{n,it}") ;
      } ;

      // Write the location.
      void compute(const sequence &seq) {
        if(*debug) Loci::debugout << "n: " << *n << " it: " << *it
          << " ZLoc,ZMax: " << ZResidualData->maxResidualLocation
          << ", " << ZResidualData->maxResidual << endl ;
      }
  } ;

  register_rule<WriteMaximumZResidualLocation>
    registerWriteMaximumZResidualLocation ;

 // Class to write maximum Zvar residual location to debug file.
  class WriteMaximumZvarResidualLocation : public singleton_rule {
    private:
      const_param<int> n,it ;
      const_param<ScalarResidual> ZvarResidualData ;
      const_param<bool> debug ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      WriteMaximumZvarResidualLocation() {
        name_store("$n{n}",n) ;
        name_store("$it{n,it}",it) ;
        name_store("ZvarResidualData{n,it}",ZvarResidualData) ;
        name_store("debug{n,it}",debug) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("$n{n},$it{n,it},ZvarResidualData{n,it},debug{n,it}") ;
        output("OUTPUT{n,it}") ;
      } ;

      // Write the location.
      void compute(const sequence &seq) {
        if(*debug) Loci::debugout << "n: " << *n << " it: " << *it
          << " ZvarLoc,ZvarMax: " << ZvarResidualData->maxResidualLocation
          << ", " << ZvarResidualData->maxResidual << endl ;
      }
  } ;

  register_rule<WriteMaximumZvarResidualLocation>
    registerWriteMaximumZvarResidualLocation ;
}
