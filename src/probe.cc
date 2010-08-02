//-----------------------------------------------------------------------------
// Description: This file contains rules for probing solution values at up to
//   ten locations within the problem domain.
//
// Authors: Original coding by Ed Luke. Modified by Jeff Wright
//-----------------------------------------------------------------------------

// Standard library includes.
#include<vector>
using std::vector ;

// Loci includes.
#include <Loci.h>

// StreamUns includes.
#include "sciTypes.h"

/******************************************************************************

  Solution Probes:

  To use:

  Put the following line in the vars file:

  probe: <probe1=[0,0,0],probe2=[1,0,1]>

  This will enable to probes (note the probe names are probe0,probe1,...probe9
  (there are a total of 10 probes available)

  The coordinate vectors are in the grid unit system (as defined by Lref).
  Probes will generate output files for each timestep of the following form
  probeX.dat. For each timestep this file contains a line with the following
  information:

  <timeStepNum> <solutionTime> <rho> <v> <p> <T> <h> <y0,y1,...> <pos> <dist>

  where

    <timeStepNum> is the time-step number in time-dependent mode and the
      iteration number in steady-state mode.
    <solutionTime> is the simulation time
    <rho> is the density
    <v> is the velocity
    <p> is the pressure
    <T> is the temperature
    <h> is the total enthalpy
    <y0,y1,...> are the species mass fractions
    <pos> is the position of the probe
    <dist> is the distance from the probed value and the position given in the
            input file.
******************************************************************************/

 

namespace streamUns {

  class DefaultProbeFrequency : public default_rule {
    private:
      param<int> probe_freq ;
    public:
      DefaultProbeFrequency() {
        name_store("probe_freq",probe_freq) ;
        output("probe_freq") ;
        comments("Frequency of output to probe files.") ;
      }
      void compute(const sequence &seq) { *probe_freq=10 ; }
  } ;

  register_rule<DefaultProbeFrequency> registerDefaultProbeFrequency ;

  class OptionalProbe : public optional_rule {
    private:
      param<options_list> probe ;
    public:
      OptionalProbe() {
        name_store("probe",probe) ;
        output("probe") ;
        comments("Put probes into the solution.  The user specifies a list of probes in the form of < probe1=[x1,y1,z1],probe2=[x2,y2,z2],... > and the solver produces probex.dat files that contains the solution located at that point over time.") ;
                                                                                
      }
      void compute(const sequence &seq) {}
  } ;

  register_rule<OptionalProbe> registerOptionalProbe ;

  // Have now included the cell number as data so we can use this as
  // additional debuggin information.
  class ProbeValues {
    public:
      real cell ;
      real distance ;
      vect3d position ;
      real density,pressure,temperature,totalEnthalpy ;
      vect3d velocity ;
      vector<real> speciesMassFraction ;
    public:
      ProbeValues() : distance(1.0e30) {}

      // Returns the serialized buffer size.
      int BufferSize() const {
        int bufferSize=13+speciesMassFraction.size() ; return bufferSize ;
      }

     // Packs the data into a buffer.
     void PackBuffer(real *buffer,int size) {
       int i=0 ; buffer[i++]=cell ; buffer[i++]=distance ;
       buffer[i++]=position.x ; buffer[i++]=position.y ; buffer[i++]=position.z ;
       buffer[i++]=density ; buffer[i++]=pressure ; buffer[i++]=temperature ;
       buffer[i++]=totalEnthalpy ; buffer[i++]=velocity.x ;
       buffer[i++]=velocity.y ; buffer[i++]=velocity.z ;
       buffer[i++]=speciesMassFraction.size() ;
       for(unsigned int j=0;j<speciesMassFraction.size();++j)
         buffer[i++]=speciesMassFraction[j] ;
     }

     // Unpacks the data from a buffer.
     void UnpackBuffer(real *buffer,int size) {
       int i=0 ; cell=buffer[i++] ; distance=buffer[i++] ; position.x=buffer[i++] ;
       position.y=buffer[i++] ; position.z=buffer[i++] ;
       density=buffer[i++] ; pressure=buffer[i++] ; temperature=buffer[i++] ;
       totalEnthalpy=buffer[i++] ; velocity.x=buffer[i++] ;
       velocity.y=buffer[i++] ; velocity.z=buffer[i++] ;
       int numSpecies=int(buffer[i++]) ;
       speciesMassFraction=vector<real>(numSpecies) ;
       for(int j=0;j<numSpecies;++j) speciesMassFraction[j]=buffer[i++] ;
    }

  } ;

  // Output operator. Required to link.
  inline ostream& operator<<(ostream &out, const ProbeValues &probeValues) {
    return out ;
  }
                                                                                
  // Input operator. Required to link.
  inline istream& operator>>(istream &in,ProbeValues &probeValues) {
    return in ;
  }

  typedef std::vector<ProbeValues> ProbeList ;
}

namespace Loci {

  // Class to serialize ProbeList.
  class ProbeListConverter {
    private:
      streamUns::ProbeList &ref ;
    public:
      explicit ProbeListConverter(streamUns::ProbeList &ref) : ref(ref) {}
    public:
      int getSize() const {
        int size=1 ;
        for(size_t i=0;i<ref.size();++i) size+=ref[i].BufferSize() ;
        return size ;
      }
      void getState(streamUns::real *buf,int &size) {
        size=getSize() ; buf[0]=ref.size() ; ++buf ;
        for(size_t i=0;i<ref.size();++i){
          ref[i].PackBuffer(buf,size) ; buf+=ref[i].BufferSize() ;
        }
      }
      void setState(streamUns::real *buf,int size) {
        ref.clear() ; int sz=int(buf[0]) ; ++buf ; ref.reserve(sz) ;
        for(int i=0;i<sz;++i) {
          streamUns::ProbeValues pV ; pV.UnpackBuffer(buf,size) ;
          ref.push_back(pV) ; buf+=pV.BufferSize() ;
        }
      }
  } ;

  // Schema converter for ProbeList.
  template<> struct data_schema_traits<streamUns::ProbeList> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef streamUns::real Converter_Base_Type ;
    typedef ProbeListConverter Converter_Type ;
  } ;
}

namespace streamUns {

  // Note that this assumes that probes are in order in the .vars file.
  class GetProbePosition : public singleton_rule {
    private:
      const_param<real> Lref ;
      const_param<options_list> probe ;
      param<std::vector<vect3d> > probePos ;
    public:
      GetProbePosition() {
        name_store("Lref",Lref) ;
        name_store("probe",probe) ;
        name_store("probePos",probePos) ;
        input("Lref,probe") ;
        output("probePos") ;
      }

      void compute(const sequence &seq) {
        options_list::option_namelist l = probe->getOptionNameList() ;
        for(options_list::option_namelist::const_iterator li=l.begin();
        li!=l.end();++li) {
          vect3d val=vect3d(0.,0.,0.) ;
          get_vect3dOption(*probe,li->c_str(),"m",val,1.0) ;
          probePos->push_back(val) ;
        }
      }
  } ;
                                                                                
  register_rule<GetProbePosition> registerGetProbePosition ;
 
  struct MinProbeList {
    void operator()(ProbeList &r,const ProbeList &s) {
      fatal(r.size() != s.size()) ;
      for(size_t i=0;i<r.size();++i) {
        if(s[i].distance<r[i].distance) r[i]=s[i] ;
      }
    }
  } ;

  class ProbeUnit : public unit_rule {
    private:
      const_param<std::vector<vect3d> > probePos ;
      param<ProbeList> probeList ;
    public:

      // Define input and output.
      ProbeUnit() {
        name_store("probePos",probePos) ;
        name_store("probeList",probeList) ;
        input("probePos") ;
        output("probeList") ;
      }

      // Initialize probe distance.
    void compute(const sequence &seq) {
      size_t sz=probePos->size() ; ProbeList tmp(sz) ; *probeList=tmp ;
    }
  } ;

  register_rule<ProbeUnit> registerProbeUnit ;

  class ProbeApplyIncompressible : public apply_rule<param<ProbeList>,
  MinProbeList> {
    private:
      const_param<std::vector<vect3d> > probePos ;
      const_param<int> numSpecies ;
      const_store<vect3d> cellcenter ;
      const_store<real> rho,p,T ;
      const_store<vect3d> v ;
      const_storeVec<real> y ;
      param<ProbeList> probeList ;
    private:
      ProbeValues pV ;
    public:

      // Define input and output.
      ProbeApplyIncompressible() {
        name_store("probePos",probePos) ;
        name_store("numSpecies",numSpecies) ;
        name_store("cellcenter",cellcenter) ;
        name_store("rho",rho) ;
        name_store("p",p) ;
        name_store("temperature",T) ;
        name_store("v",v) ;
        name_store("y",y) ;
        name_store("probeList",probeList) ;
        input("probePos,numSpecies,cellcenter,rho,p,temperature,v,y") ;
        input("probeList") ;
        output("probeList") ;
        constraint("incompressibleFlow,geom_cells") ;
      }

      // See if a cell is the current closest cell.
      void calculate(Entity cell) {
        fatal((*probePos).size()!=(*probeList).size()) ;
        for(size_t i=0;i<probeList->size();++i) {
          const vect3d dR=cellcenter[cell]-probePos[cell][i] ;
          const real dist=dot(dR,dR) ;
          if(dist<probeList[cell][i].distance) {
            pV.cell=real(cell) ; pV.distance=dist ; pV.position=cellcenter[cell] ;
            pV.density=rho[cell] ; pV.velocity=v[cell] ; pV.pressure=p[cell] ;
            pV.temperature=T[cell] ; pV.totalEnthalpy=0.0 ;
            if(*numSpecies>0) for(int i=0;i<*numSpecies;++i) pV.
              speciesMassFraction[i]=y[cell][i] ;
            probeList[cell][i]=pV ;
          }
        }
      }

      // Find the nearest cell values.
      virtual void compute(const sequence &seq) {
        if(*numSpecies>0) pV.speciesMassFraction=vector<real>(*numSpecies) ;
        do_loop(seq,this) ;
      }
  } ;

  register_rule<ProbeApplyIncompressible> registerProbeApplyIncompressible ;

  class ProbeApplyCompressible : public apply_rule<param<ProbeList>,
  MinProbeList> {
    private:
      const_param<std::vector<vect3d> > probePos ;
      const_param<int> numSpecies ;
      const_store<vect3d> cellcenter ;
      const_store<real> rho,p,T,h ;
      const_store<vect3d> v ;
      const_storeVec<real> y ;
      param<ProbeList> probeList ;
    private:
      ProbeValues pV ;
    public:

      // Define input and output.
      ProbeApplyCompressible() {
        name_store("probePos",probePos) ;
        name_store("numSpecies",numSpecies) ;
        name_store("cellcenter",cellcenter) ;
        name_store("rho",rho) ;
        name_store("p",p) ;
        name_store("temperature",T) ;
        name_store("h",h) ;
        name_store("v",v) ;
        name_store("y",y) ;
        name_store("probeList",probeList) ;
        input("probePos,numSpecies,cellcenter,rho,p,temperature,h,v,y") ;
        input("probeList") ;
        output("probeList") ;
        constraint("compressibleFlow,geom_cells") ;
      }

      // See if a cell is the current closest cell.
      void calculate(Entity cell) {
        fatal((*probePos).size()!=(*probeList).size()) ;
        for(size_t i=0;i<probeList->size();++i) {
          const vect3d dR=cellcenter[cell]-probePos[cell][i] ;
          const real dist=dot(dR,dR) ;
          if(dist<probeList[cell][i].distance) {
            pV.cell=real(cell) ; pV.distance=dist ; pV.position=cellcenter[cell] ;
            pV.density=rho[cell] ; pV.velocity=v[cell] ; pV.pressure=p[cell] ;
            pV.temperature=T[cell] ; pV.totalEnthalpy=h[cell] ;
            if(*numSpecies>0) for(int i=0;i<*numSpecies;++i) pV.
              speciesMassFraction[i]=y[cell][i] ;
            probeList[cell][i]=pV ;
          }
        }
      }

      // Find the nearest cell values.
      virtual void compute(const sequence &seq) {
        if(*numSpecies>0) pV.speciesMassFraction=vector<real>(*numSpecies) ;
        do_loop(seq,this) ;
      }
  } ;

  register_rule<ProbeApplyCompressible> registerProbeApplyCompressible ;

  // Rule to set the probing flag.
  class DoProbe : public singleton_rule {
    private:
      const_param<int> n ;
      const_param<int> it ;
      const_param<int> numTimeSteps ;
      const_param<int> maxIterationsPerTimeStep ;
      const_param<int> probeFrequency ;
      param<bool> doProbe ;
    public:

      // Define input and output.
      DoProbe() {
        name_store("ncycle{n}",n) ;
        name_store("$it{n,it}",it) ;
        name_store("probe_freq{n,it}",probeFrequency) ;
        name_store("do_probe{n,it}",doProbe) ;
        input("ncycle{n},$it{n,it}") ;
        input("probe_freq{n,it}") ;
        output("do_probe{n,it}") ;
      }

      virtual void compute(const sequence &seq) {
        doProbe=(*n % *probeFrequency==0 && *it==0) ;
      }
  } ;

  register_rule<DoProbe> registerDoProbe ;
  
  // Rule to write out the probe values.
  class OutputProbe : public pointwise_rule {
    private:
      const_param<options_list> probe ;
      const_param<ProbeList> probeList ;
      const_param<int> timeStepNum ;
      const_param<real> solutionTime ;
      const_param<int> numTimeSteps ;
      const_param<int> it ;
      const_param<int> numSpecies ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      OutputProbe() {
        name_store("probe{n,it}",probe) ;
        name_store("probeList{n,it}",probeList) ;
        name_store("ncycle{n}",timeStepNum) ;
        name_store("stime{n}",solutionTime) ;
        name_store("numTimeSteps{n,it}",numTimeSteps) ;
        name_store("$it{n,it}",it) ;
        name_store("numSpecies{n,it}",numSpecies) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("probe{n,it},probeList{n,it}") ;
        input("ncycle{n},stime{n},numTimeSteps{n,it},$it{n,it}") ;
        input("numSpecies{n,it}") ;
        output("OUTPUT{n,it}") ;
        conditional("do_probe{n,it}") ;
      }

      // Write out the probe data. The fixed format with 16 decimal places of
      // precision is used because Jeff West's awk script cannont handle
      // exponential numbers.
      void compute(const sequence &seq) {
        if(Loci::MPI_rank==0){
          size_t i = 0 ;
          options_list::option_namelist l = probe->getOptionNameList() ;
          for(options_list::option_namelist::const_iterator li=l.begin();
          li!=l.end();++li,++i) {
            string filename=*li ; filename += ".dat" ; ofstream out ;
            if(*numTimeSteps==1){
              out.open(filename.c_str(),ios::app) ;
              if(out.fail()) out.open(filename.c_str(),ios::out) ;
            }else if(*timeStepNum==0){
              out.open(filename.c_str(),ios::out) ;
            }else{
              out.open(filename.c_str(),ios::app) ;
            }
            if(out.fail()) {
              cerr << "Probe: can't open "<< filename << endl ; return ;
            }
            out.setf(ios::fixed) ; out.precision(16) ;
            if(*numTimeSteps==1) out << *it ; else out << *timeStepNum ;
            out << " " << *solutionTime << " " << (*probeList)[i].density
              << " " << (*probeList)[i].velocity << " "
              << (*probeList)[i].pressure << " " << (*probeList)[i].temperature
              << " " << (*probeList)[i].totalEnthalpy << " " ;
            if(*numSpecies>0){
              for(int j=0;j<*numSpecies;++j) out << " " << (*probeList)[i].
                speciesMassFraction[j] ;
              out << " " ;
            }
            out << (*probeList)[i].position << " " << sqrt((*probeList)[i].
              distance) << endl ;
            out.close() ;
//          cout << "closest cell to probe " << i << ": "
//            << (*probeList)[i].cell << endl ;
          }
        }
      }
  } ;

  register_rule<OutputProbe> registerOutputProbe ;
}
