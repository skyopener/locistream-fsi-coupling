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

  // Class to write Alpha residuals to screen.
  class PrintAlphaResidual : public singleton_rule {
    private:
      const_param<int> n,nCycle,it ;
      const_param<ScalarResidual> AlphaResidualData ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      PrintAlphaResidual() {
        name_store("$n{n}",n) ;
	name_store("ncycle{n}",nCycle) ;
        name_store("$it{n,it}",it) ;
        name_store("AlphaResidualData{n,it}",AlphaResidualData) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("$n{n},ncycle{n},$it{n,it},AlphaResidualData{n,it}") ;
        output("OUTPUT{n,it}") ;
      } ;

      void compute(const sequence &seq) {

        // Set the output format.
        if(Loci::MPI_rank==0){
          cout.setf(ios::scientific,ios::floatfield) ; cout.precision(6) ;
        }  

        // Write out the residuals.
        if(Loci::MPI_rank==0){
	  cout <<"RAlpha: " << *nCycle << " " << *it << " " << AlphaResidualData->
            totalResidual << endl;
        }

      }
  } ;

  register_rule<PrintAlphaResidual>
    registerPrintAlphaResidual ;

  // Class to write maximum Alpha residual location to debug file.
  class WriteMaximumAlphaResidualLocation : public singleton_rule {
    private:
      const_param<int> n,it ;
      const_param<ScalarResidual> AlphaResidualData ;
      const_param<bool> debug ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      WriteMaximumAlphaResidualLocation() {
        name_store("$n{n}",n) ;
        name_store("$it{n,it}",it) ;
        name_store("AlphaResidualData{n,it}",AlphaResidualData) ;
        name_store("debug{n,it}",debug) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("$n{n},$it{n,it},AlphaResidualData{n,it},debug{n,it}") ;
        output("OUTPUT{n,it}") ;
      } ;

      // Write the location.
      void compute(const sequence &seq) {
        if(*debug) Loci::debugout << "n: " << *n << " it: " << *it
          << " AlphaLoc,AlphaMax: " << AlphaResidualData->maxResidualLocation
          << ", " << AlphaResidualData->maxResidual << endl ;
      }
  } ;

  register_rule<WriteMaximumAlphaResidualLocation>
    registerWriteMaximumAlphaResidualLocation ;

}
