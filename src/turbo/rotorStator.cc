//-----------------------------------------------------------------------------
// Description: This file contains data structures and rules used to define
//   the rotor and stator components of the geometry.
//-----------------------------------------------------------------------------
                                                                                
// Standard library includes.
#include <vector>
using std::vector ;
                                                                                
// Loci includes.
#include <Loci.h>
#include <Tools/parse.h>
using namespace Loci::parse ;
                                                                                
// StreamUns includes.
#include "const.h"
#include "rotorStator.h"
#include "sciTypes.h"
#include "varsFileInputs.h"
                                                                                
namespace streamUns {

  //---------------------------------------------------------------------------
  // Class Rotor.

  // Returns the serialized buffer size.
  int Rotor::BufferSize() const { return 8 ; }

  // Overridden virtual method called by fact database when a fact of type
  // Rotor is read.
  istream& Rotor::Input(istream &in) {
    RotorOptions optionsList ;
    in >> optionsList ; Input(optionsList) ; return in ;
  }

  // Gets the values that define the rotor from the specified options. All
  // values are defined using the same coordinate system as the grid. These
  // values are automatically scaled by Lref here.
  void Rotor::Input(const options_list &optionsList) {

    // Get the reference frame number.
    if(optionsList.optionExists("refFrameNum")){
      double tmp ; optionsList.getOption("refFrameNum",tmp) ;
      referenceFrameNumber=int(tmp) ;
    }else{
      cerr << "ERROR: Must specify a reference frame number for 'rotor'."
        << endl ;
    }

    // Get the center of the cylinder which defines the extent of the rotor.
    if(optionsList.optionExists("center")){
      vect3d val=vect3d(0.,0.,0.) ;
      get_vect3dOption(optionsList,"center","m",val,Lref) ;
      center=val ;
    }else{
      cerr << "ERROR: Must specify center for 'rotor' ." << endl ;
    }

    // Get the cylinder radius.
    if(optionsList.optionExists("radius")){
      double tmp ; optionsList.getOption("radius",tmp) ; radius=tmp*Lref ;
    }else{
      cerr << "ERROR: Must specify a radius for 'rotor'." << endl ;
    }

    // Get the cylinder height.
    if(optionsList.optionExists("length")){
      double tmp ; optionsList.getOption("length",tmp) ; length=tmp*Lref ;
    }else{
      cerr << "ERROR: Must specify a length for 'rotor'." << endl ;
    }

    // Get the axial pressure gradient.
    if(optionsList.optionExists("dpdz")){
      Loci::option_value_type optionValueType=optionsList.getOptionValueType
        ("dpdz") ;
      if(optionValueType==Loci::FUNCTION){
        Loci::options_list::arg_list options ; string name ;
        optionsList.getOption("dpdz",name,options) ;
        if(name=="radial"){
          string validOptions="a:b" ;
          Loci::options_list radialOptions(validOptions) ;
          radialOptions.Input(options) ;
          dpdzA=dpdzB=0.0 ;
          if(radialOptions.optionExists("a")){
            real a ; radialOptions.getOption("a",a) ; dpdzA=a ;
          }
          if(radialOptions.optionExists("b")){
            real b ; radialOptions.getOption("b",b) ; dpdzB=b ;
          }
          cout << "dpdZA,dpdzB: " << dpdzA << " " << dpdzB << endl ;
        }else{
           cerr << "ERROR: Syntax for option dpdz is 'dpdz=radial(a,b)'"
             << endl ;
        }
      }else{
         cerr << "ERROR: Syntax for option dpdz is 'dpdz=radial(a,b)'" << endl ;
      }
    }
  }

  // Checks if the cell center is in the reference frame.
  bool Rotor::InReferenceFrame(const vect3d& cellCenter,const
  vector<ReferenceFrame>& referenceFrame) const {

    // Compute the axis unit vector.
    vect3d axis=referenceFrame[referenceFrameNumber].axisEnd-referenceFrame
      [referenceFrameNumber].axisStart ; axis/=norm(axis) ;

    // Find the radial and axial distance of the cell center from the
    // reference frame center. Use these values to check whether the cell is
    // in the reference frame cylinder defined by 'radius', 'length' and
    // 'center' .
    real deltaR=norm(cross(cellCenter-center,axis)),deltaX=abs(dot(cellCenter-
      center,axis)) ;
    if(deltaR<radius && deltaX<length/2) return true ;
    return false ;
  }
                                                                                
  // Packs the data into a buffer.
  void Rotor::PackBuffer(real *buffer,int size) {
    int i=0 ;
    buffer[i++]=referenceFrameNumber ; buffer[i++]=center.x ;
    buffer[i++]=center.y ; buffer[i++]=center.z ; buffer[i++]=radius ;
    buffer[i++]=length ; buffer[i++]=dpdzA ; buffer[i++]=dpdzB ;
  }
                                                                                
  // Prints the rotor information.
  ostream& Rotor::Print(ostream &out) const {
    out << "<" ;
    out << "refFrameNum=" << referenceFrameNumber << "," ;
    out << "center=[" << center.x << "," << center.y << "," << center.z
      << "]," ;
    out << "radius=" << radius << "," ; out << "length=" << length << "," ;
    out << "dpdzA=" << dpdzA ; out << "dpdzB=" << dpdzB ; out << ">" ;
    return out ;
  }
                                                                                
  // Unpacks the data from a buffer.
  void Rotor::UnpackBuffer(real *buffer,int size) {
    int i=0 ; referenceFrameNumber=(int)buffer[i++] ;
    center.x=buffer[i++] ; center.y=buffer[i++] ; center.z=buffer[i++] ;
    radius=buffer[i++] ; length=buffer[i++] ; dpdzA=buffer[i++] ;
    dpdzB=buffer[i++] ;
  }

  //---------------------------------------------------------------------------
  // Class RotorData.

  // Returns the serialized buffer size.
  int RotorData::BufferSize() const {
    int size=0 ;
    for(unsigned int i=0;i<rotor.size();++i) size+=rotor[i].BufferSize() ;
    return size ;
  }

  // Overridden virtual method called by fact database when a fact of type
  // RotorData is read. This will never be used, so it is empty.
  istream& RotorData::Input(istream &in) { return in ; }

  // Packs the data into a buffer.
  void RotorData::PackBuffer(real *buffer,int size) {
    *buffer=rotor.size() ; ++buffer ;
    for(unsigned int j=0;j<rotor.size();++j){
      rotor[j].PackBuffer(buffer,size) ; buffer+=rotor[j].BufferSize() ;
    }
  }

  // Gets the pressure gradient for a cell given the cell center.
  real RotorData::DPDZA(const vect3d &cellCenter,const vector<ReferenceFrame>
  &referenceFrame) const {
    bool found,foundAtLeastOnce=false ; real dpdzA=0.0 ;
    for(unsigned int i=0;i<rotor.size();++i){
      if((found=rotor[i].InReferenceFrame(cellCenter,referenceFrame)))
        dpdzA=rotor[i].DPDZA() ;
      foundAtLeastOnce|=found ;
    }
    if(!foundAtLeastOnce){
      cerr << "ERROR: Could not find a reference frame for cell located at "
        << cellCenter << " ." << endl ; Loci::Abort() ;
    }
    return dpdzA ;
  }

  // Gets the pressure gradient for a cell given the cell center.
  real RotorData::DPDZB(const vect3d &cellCenter,const vector<ReferenceFrame>
  &referenceFrame) const {
    bool found,foundAtLeastOnce=false ; real dpdzB=0.0 ;
    for(unsigned int i=0;i<rotor.size();++i){
      if((found=rotor[i].InReferenceFrame(cellCenter,referenceFrame)))
        dpdzB=rotor[i].DPDZB() ;
      foundAtLeastOnce|=found ;
    }
    if(!foundAtLeastOnce){
      cerr << "ERROR: Could not find a reference frame for cell located at "
        << cellCenter << " ." << endl ; Loci::Abort() ;
    }
    return dpdzB ;
  }

  // Prints the rotor file information.
  ostream& RotorData::Print(ostream &out) const {
    out << "<" ;
    for(unsigned int i=0;i<rotor.size();++i){
      rotor[i].Print(out) ; if(i!=rotor.size()-1) out << endl << endl ;
    }
    out << ">" ; return out ;
  }

  // Reads the rotor file.
  istream& RotorData::Read(istream &in,const real &Lref) {
    kill_white_space(in) ;
    while(in.good() && !in.eof() && in.peek()!=char_traits<char>::eof()) {
      kill_white_space(in) ;
      if(get_token(in,"rotor")) {
        kill_white_space(in) ;
        if(in.peek()!='=') {
          cerr << "Error: Bad 'rotor' specification." << endl ;
          Loci::Abort() ; return in ;
        }
        in.get() ;
        Rotor r(Lref) ; r.Input(in) ; rotor.push_back(r) ;
        kill_white_space(in) ;
        if(in.peek()!=';') {
          cerr << "Error: Expected ';' after 'rotor' definition ." << endl ;
          Loci::Abort() ; return in ;
        }else in.get() ;
      }else{
        kill_white_space(in) ;
        if(in.peek() != char_traits<char>::eof()) {
          cerr << "Error: Unexpected input in rotor file." << endl ;
          Loci::Abort() ; return in ;
        }
      }
    }
    return in ;
  }

  // Gets the reference frame for a cell given the cell center.
  unsigned int RotorData::ReferenceFrameNumber(const vect3d &cellCenter,const
  vector<ReferenceFrame> &referenceFrame) const {
    bool found,foundAtLeastOnce=false ; int referenceFrameNumber=0 ;
    for(unsigned int i=0;i<rotor.size();++i){
      if((found=rotor[i].InReferenceFrame(cellCenter,referenceFrame)))
        referenceFrameNumber=rotor[i].ReferenceFrameNumber() ;
      foundAtLeastOnce|=found ;
    }
    if(!foundAtLeastOnce){
      cerr << "ERROR: Could not find a reference frame for cell located at "
        << cellCenter << " ." << endl ; Loci::Abort() ;
    }
    return referenceFrameNumber ;
  }

  // Unpacks the data from a buffer.
  void RotorData::UnpackBuffer(real *buffer,int size) {
    int numRotor=int(*buffer) ; ++buffer ; rotor=vector<Rotor>(numRotor) ;
    for(int j=0;j<numRotor;++j){
      rotor[j].UnpackBuffer(buffer,size) ; buffer+=rotor[j].BufferSize() ;
    }
  }

  //---------------------------------------------------------------------------
  // Setup Rules.

  // Rule to declare "rotorFile" as valid .vars file variable. We are making
  // this a default rule so we can always use value to determine whether this
  // variable has been included in the .vars file.
  class DefaultRotorFile : public default_rule {
    private:
      param<string> rotorFile ;
    public:

      // Define input and output.
      DefaultRotorFile() {
        name_store("rotorFile",rotorFile) ;
        output("rotorFile") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) { *rotorFile="" ; }
  } ;

  register_rule<DefaultRotorFile> registerDefaultRotorFile ;

  // Rule to declare "turboMachinery" as a valid .vars file option. The
  // default value is true indicating that source terms due to rotating
  // reference frames are added.
  class DefaultTurboMachinery : public default_rule {
    private:
      param<string> turboMachinery ;
    public:

      // Define input and output.
      DefaultTurboMachinery() {
        name_store("turboMachinery",turboMachinery) ;
        output("turboMachinery") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) { *turboMachinery="no" ; }
  } ;

  register_rule<DefaultTurboMachinery> registerDefaultTurboMachinery ;

  // Creates the turbomachinery constraint. This rule used to use the
  // presence of a rotor file to activate turbomachinery, but now this has
  // been changed since we may use the rotor file even when we don't want
  // the turbomachinery source terms to be active.
  class TurbomachineryConstraint : public constraint_rule {
    private:
//    const_param<string> rotorFile ;
      const_param<string> turboMachinery ;
      Constraint turbomachinery ;
    public:
                                                                                
      // Define input and output.
      TurbomachineryConstraint() {
//      name_store("rotorFile",rotorFile) ;
        name_store("turboMachinery",turboMachinery) ;
        name_store("turbomachinery",turbomachinery) ;
        input("turboMachinery") ;
        output("turbomachinery") ;
      }
                                                                                
      // Set up the constraint.
      virtual void compute(const sequence& seq) {
        if(*turboMachinery=="no") turbomachinery=EMPTY ;
        else turbomachinery=~EMPTY ;
      }
  } ;
                                                                                
  register_rule<TurbomachineryConstraint> registerTurbomachineryConstraint ;

  // Reads the rotor data from file.
  class ReadRotorDataFromFile : public singleton_rule {
    private:
      const_param<real> Lref ;
      const_param<string> rotorFile ;
      param<RotorData> rotorData ;
                                                                                
    public:
                                                                                
      // Define input and output.
      ReadRotorDataFromFile() {
        name_store("Lref",Lref) ;
        name_store("rotorFile",rotorFile) ;
        name_store("rotorData",rotorData) ;
        input("Lref,rotorFile") ;
        output("rotorData") ;
      }
                                                                                
      // Read in the data from the rotor file.
      virtual void compute(const sequence &seq) {
        string fileName=(*rotorFile) ;
        if(Loci::MPI_rank==0) cout << "Reading rotor file: " << fileName
          << endl ;
        ifstream in(fileName.c_str(),ios::in) ;
        if(in.fail()){
          if(Loci::MPI_rank==0) cout << "ERROR: Could not read rotor "
            << "file: " << fileName << endl ;
          Loci::Abort() ;
        }else{
          rotorData->Read(in,*Lref) ;
        }
      }
  } ;
                                                                                
  register_rule<ReadRotorDataFromFile> registerReadRotorDataFromFile ;

  // Rule to determine the reference frame for each cell.
  class CellReferenceFrame : public pointwise_rule {
    private:
      const_param<RotorData> rotorData ;
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<vect3d> cellCenter ;
      store<unsigned int> cellReferenceFrame ;
    public:
                                                                                
      // Define input and output.
      CellReferenceFrame() {
        name_store("rotorData",rotorData) ;
        name_store("referenceFrame",referenceFrame) ;
        name_store("cellcenter",cellCenter) ;
        name_store("cellReferenceFrame",cellReferenceFrame) ;
        input("rotorData,referenceFrame,cellcenter") ;
        output("cellReferenceFrame") ;
        constraint("geom_cells") ;
      }
                                                                                
      // Get the reference frame for a single cell.
      void calculate(Entity cell) {
        cellReferenceFrame[cell]=rotorData->ReferenceFrameNumber
          (cellCenter[cell],*referenceFrame) ;
      }
                                                                                
      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<CellReferenceFrame> registerCellReferenceFrame ;

}
