//-----------------------------------------------------------------------------
// Description: This file contains rules for reference frames.
//
// Author: Jeff Wright
//-----------------------------------------------------------------------------

// Standard library includes.
#include <vector>
using std::vector ;
                                                                                
// Loci includes.
#include <Loci.h>
                                                                                
// StreamUns includes.
#include "referenceFrame.h"
#include "sciTypes.h"
                                                                                
namespace streamUns {

  // Reads the reference frames from file.
  class ReferenceFrameInput : public singleton_rule {
    private:
      const_param<string> referenceFrameFile ;
      param<vector<ReferenceFrame> > referenceFrame ;
    public:
                                                                                
      // Define input and output.
      ReferenceFrameInput() {
        name_store("referenceFrameFile",referenceFrameFile) ;
        name_store("referenceFrame",referenceFrame) ;
        input("referenceFrameFile") ;
        output("referenceFrame") ;
      }
                                                                                
      // Set paramter values.
      void compute(const sequence &seq) {
        string fileName=*referenceFrameFile ;
        if(Loci::MPI_rank==0) cout << "Reading reference frames from "
          << fileName << endl ;
        ifstream in(fileName.c_str(),ios::in) ;
        if(in.fail()) {
          if(Loci::MPI_rank==0) cerr << "Can't open file: " << fileName
            << endl ; Loci::Abort() ;
        }
        int numReferenceFrame ;
        if(!(in >> numReferenceFrame) || numReferenceFrame<1){
          if(Loci::MPI_rank==0)
            cerr << "Bad number of reference frames in file: " << fileName
              << endl ; Loci::Abort() ;
        }
        *referenceFrame=vector<ReferenceFrame>(numReferenceFrame) ;
        for(int i=0;i<numReferenceFrame;++i){
          if(!(in >> (*referenceFrame)[i])){
            if(Loci::MPI_rank==0) cerr << "Bad data for reference frame "
              << i << " in file: " << fileName << endl ; Loci::Abort() ;
          }
        }
      }
  } ;
                                                                                
  register_rule<ReferenceFrameInput> registerReferenceFrameInput ;
}
