#line 1 "timeStep.loci"
// Standard library includes.
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using std::cerr ;
using std::cout ;
using std::endl ;
using std::string ;
using std::vector ;

// Loci includes.
#include <Loci.h>
#line 1 "FVM.lh"
//#############################################################################
//#
//# Copyright 2008, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################

// $type pos store<Loci::vector3d<Loci::real_t> > 
// $type cl Map
// $type cr Map
// $type ci Map
// $type ref Map
// $type pmap Map
// $type face2node multiMap

// $type upper multiMap
// $type lower multiMap
// $type boundary_map multiMap

// $type cellcenter store<Loci::vector3d<Loci::real_t> > 
// $type facecenter store<Loci::vector3d<Loci::real_t> > 
// $type area store<Loci::Area> 
// $type vol store<Loci::real_t> 
// $type grid_vol param<Loci::real_t> 

// $type mn store<Loci::vector3d<Loci::real_t> > 
// $type ln store<Loci::vector3d<Loci::real_t> > 

// $type grads(X0) store<Loci::vector3d<Loci::real_t> > 
// $type gradv(X0) storeVec<Loci::vector3d<Loci::real_t> > 
// $type gradv3d(X0) store<Loci::tensor3d<Loci::real_t> > 

// $type grads_f(X0) store<Loci::vector3d<Loci::real_t> > 
// $type gradv_f(X0) storeVec<Loci::vector3d<Loci::real_t> > 
// $type gradv3d_f(X0) store<Loci::tensor3d<Loci::real_t> > 

// $type limiters(X0) store<Loci::real_t> 
// $type limiterv(X0) storeVec<Loci::real_t> 
// $type limiterv3d(X0) store<Loci::vector3d<Loci::real_t> > 

// $type lefts(X0) store<Loci::real_t> 
// $type rights(X0) store<Loci::real_t> 
// $type leftsP(X0,X1) store<Loci::real_t> 
// $type rightsP(X0,X1) store<Loci::real_t> 
// $type leftvM(X0) storeVec<Loci::real_t> 
// $type rightvM(X0) storeVec<Loci::real_t> 
// $type leftv3d(X0) store<Loci::vector3d<Loci::real_t> > 
// $type rightv3d(X0) store<Loci::vector3d<Loci::real_t> > 

// $type cell2node(X0) store<float> 
// $type cell2node_v(X0) storeVec<float> 
// $type cell2node_v3d(X0) store<Loci::vector3d<float> > 
// $type cell2nodeMax(X0) store<float> 
// $type cell2nodeMin(X0) store<float> 
// $type cell2nodeMaxMag(X0) store<float> 
// $type cell2nodeMaxv3d(X0) store<vector3d<float> > 

// $type BC_options store<Loci::options_list> 

// $type integrateSurface(X0) store<Loci::real_t> 
// $type integrateFlux(X0) store<Loci::real_t> 

// $type petscScalarSolve(X0) store<Loci::real_t> 
// $type petscBlockedSolve(X0) storeVec<Loci::real_t> 
// $type petscBlockedSSolve(X0) storeVec<Loci::real_t> 

// $type L2Norm(X0) param<Loci::real_t> 
// $type L1Norm(X0) param<Loci::real_t> 
// $type LinfNorm(X0) param<Loci::real_t> 
#line 14 "timeStep.loci"

using Loci::Area ;

// StreamUns includes.
#include "sciTypes.h"

namespace streamUns {

  // Class for timestep ramp options.
  class TimeStepRamp : public options_list {
    public:
      TimeStepRamp() : options_list("ramp0:ramp1:ramp2") {}
  } ;

  // Optional vars file input for timestep ramping.
  // $type timeStepRamp param<TimeStepRamp> 
  namespace {class file_timeStep000_1280810698m592 : public Loci::optional_rule {
#line 33 "timeStep.loci"
    Loci::param<TimeStepRamp>  L_timeStepRamp_ ; 
#line 33 "timeStep.loci"
public:
#line 33 "timeStep.loci"
    file_timeStep000_1280810698m592() {
#line 33 "timeStep.loci"
       name_store("timeStepRamp",L_timeStepRamp_) ;
#line 33 "timeStep.loci"
       output("timeStepRamp") ;
#line 33 "timeStep.loci"
       comments("Usage <ramp0=[n=10,dt=1.0e-05],ramp1=[n=20,dt=1.0e-04],...>.") ;
#line 33 "timeStep.loci"
       comments("n: number of timesteps, dt: timestep size") ;
#line 33 "timeStep.loci"
       comments("Three available ramps, ramp0, ramp1 and ramp2.") ;
#line 33 "timeStep.loci"
    }
#line 33 "timeStep.loci"
    void compute(const Loci::sequence &seq) { 
  }} ;
#line 34 "timeStep.loci"
Loci::register_rule<file_timeStep000_1280810698m592> register_file_timeStep000_1280810698m592 ;
#line 34 "timeStep.loci"
}
#line 34 "timeStep.loci"


  // Parse the ramp options and create ramp data.
  // $type rampNumTimeSteps param<vector<int> > 
  // $type rampTimeStep param<vector<real> > 
  namespace {class file_timeStep001_1280810698m593 : public Loci::singleton_rule {
#line 39 "timeStep.loci"
    Loci::const_param<TimeStepRamp>  L_timeStepRamp_ ; 
#line 39 "timeStep.loci"
    Loci::param<vector<int> >  L_rampNumTimeSteps_ ; 
#line 39 "timeStep.loci"
    Loci::param<vector<real> >  L_rampTimeStep_ ; 
#line 39 "timeStep.loci"
public:
#line 39 "timeStep.loci"
    file_timeStep001_1280810698m593() {
#line 39 "timeStep.loci"
       name_store("timeStepRamp",L_timeStepRamp_) ;
#line 39 "timeStep.loci"
       name_store("rampNumTimeSteps",L_rampNumTimeSteps_) ;
#line 39 "timeStep.loci"
       name_store("rampTimeStep",L_rampTimeStep_) ;
#line 39 "timeStep.loci"
       input("timeStepRamp") ;
#line 39 "timeStep.loci"
       output("rampNumTimeSteps") ;
#line 39 "timeStep.loci"
       output("rampTimeStep") ;
#line 39 "timeStep.loci"
    }
#line 39 "timeStep.loci"
    void compute(const Loci::sequence &seq) { 
    (*L_rampNumTimeSteps_).clear() ; (*L_rampTimeStep_).clear() ;
    for(int i=0;i<100;++i){
      ostringstream rampStream ; rampStream << "ramp" << i ;
      string rampName=rampStream.str() ;
      if((*L_timeStepRamp_).optionExists(rampName)){
        if((*L_timeStepRamp_).getOptionValueType(rampName)!=Loci::LIST){
          cerr << "Bad format for " << rampName << endl ; Loci::Abort() ;
        }else{
          bool nSet=false,dtSet=false ;
          Loci::options_list::arg_list rampList ;
          (*L_timeStepRamp_).getOption(rampName,rampList) ;
          for(Loci::options_list::arg_list::iterator p=rampList.begin();
          p!=rampList.end();++p) {
            Loci::option_values::value_list_type rampArg ; p->get_value(rampArg) ;
            if(p->type_of()!=Loci::NAME_ASSIGN || rampArg.size()!=1 ||
            rampArg.front().type_of()!=Loci::REAL){
              cerr << "Bad format for " << rampName << endl ; Loci::Abort() ;
            }
            string argName ; p->get_value(argName) ;
            double argValue ; rampArg.front().get_value(argValue) ;
            if(argName=="n"){
              nSet=true ;
              (*L_rampNumTimeSteps_).push_back(int(argValue+0.01)) ;
            }else if(argName=="dt"){
              dtSet=true ;
              (*L_rampTimeStep_).push_back(argValue) ;
            }else{
              cerr << "ERROR: Argument " << argName << " not valid for ramp."
                << endl ;
            }
          }
          if(!nSet || !dtSet){
            cerr << "ERROR: Must set both n and dt for " << rampName << endl ;
            Loci::Abort() ;
          }
        }
      }else{
        return ;
      }
    }
  }} ;
#line 80 "timeStep.loci"
Loci::register_rule<file_timeStep001_1280810698m593> register_file_timeStep001_1280810698m593 ;
#line 80 "timeStep.loci"
}
#line 80 "timeStep.loci"


  // Default vars file input for timestep ramp type.
  // $type timeStepRampType param<string> 
  namespace {class file_timeStep002_1280810698m594 : public Loci::default_rule {
#line 86 "timeStep.loci"
    Loci::param<string>  L_timeStepRampType_ ; 
#line 86 "timeStep.loci"
public:
#line 86 "timeStep.loci"
    file_timeStep002_1280810698m594() {
#line 86 "timeStep.loci"
       name_store("timeStepRampType",L_timeStepRampType_) ;
#line 86 "timeStep.loci"
       output("timeStepRampType") ;
#line 86 "timeStep.loci"
       comments("Type of interpolation used in timestep ramp.") ;
#line 86 "timeStep.loci"
       comments("Options are 'constant'.") ;
#line 86 "timeStep.loci"
    }
#line 86 "timeStep.loci"
    void compute(const Loci::sequence &seq) { 
    (*L_timeStepRampType_)=string("constant") ;
  }} ;
#line 88 "timeStep.loci"
Loci::register_rule<file_timeStep002_1280810698m594> register_file_timeStep002_1280810698m594 ;
#line 88 "timeStep.loci"
}
#line 88 "timeStep.loci"


  // Build rule for the timestep. Doing it this way, we just need to
  // then have another rule that defines dt{n}, and then we will have
  // access to both dt{n-1} and dt{n}. The bogus value assigned here will
  // never be used, but we need a non-zero value here, otherwise the
  // product dt*Q-P will be zero and we will get an FPE in later rules. JW 04/09/2010
  // $type timeStep param<real> 
  // $type dt param<real> 
  namespace {class file_timeStep003_1280810698m595 : public Loci::singleton_rule {
#line 97 "timeStep.loci"
    Loci::const_param<real>  L_timeStep_ ; 
#line 97 "timeStep.loci"
    Loci::param<real>  L_dt_n_EQ__M_1__ ; 
#line 97 "timeStep.loci"
public:
#line 97 "timeStep.loci"
    file_timeStep003_1280810698m595() {
#line 97 "timeStep.loci"
       name_store("timeStep",L_timeStep_) ;
#line 97 "timeStep.loci"
       name_store("dt{n=-1}",L_dt_n_EQ__M_1__) ;
#line 97 "timeStep.loci"
       input("timeStep") ;
#line 97 "timeStep.loci"
       output("dt{n=-1}") ;
#line 97 "timeStep.loci"
    }
#line 97 "timeStep.loci"
    void compute(const Loci::sequence &seq) { 
    (*L_dt_n_EQ__M_1__)=(*L_timeStep_);
  }} ;
#line 99 "timeStep.loci"
Loci::register_rule<file_timeStep003_1280810698m595> register_file_timeStep003_1280810698m595 ;
#line 99 "timeStep.loci"
}
#line 99 "timeStep.loci"


  // Unit/apply sequence used to get around the Loci restriction that we cannot use
  // priority rules to override params. So, we use a pair<>, where the first value
  // represents the priority and the second value the timestep. The definition of the
  // > and < operators for pair<> hold that if the first value of the pair is not
  // identical, then the first value is used for comparison. If the first values are
  // equal, then the second is used for comparison. We can thus use Loci::Maximum
  // with the first pair value being set as the priority.
  // $type dtTemp param<pair<int,real> > 
  namespace {class file_timeStep004_1280810698m595 : public Loci::unit_rule {
#line 109 "timeStep.loci"
    Loci::param<pair<int,real> >  L_dtTemp_n__ ; 
#line 109 "timeStep.loci"
public:
#line 109 "timeStep.loci"
    file_timeStep004_1280810698m595() {
#line 109 "timeStep.loci"
       name_store("dtTemp{n}",L_dtTemp_n__) ;
#line 109 "timeStep.loci"
       output("dtTemp{n}") ;
#line 109 "timeStep.loci"
       constraint("UNIVERSE{n}") ;
#line 109 "timeStep.loci"
    }
#line 109 "timeStep.loci"
    void compute(const Loci::sequence &seq) { 
    (*L_dtTemp_n__)=pair<int,real>(0,0.0) ;
  }} ;
#line 111 "timeStep.loci"
Loci::register_rule<file_timeStep004_1280810698m595> register_file_timeStep004_1280810698m595 ;
#line 111 "timeStep.loci"
}
#line 111 "timeStep.loci"


  // Apply rule that uses value from the vars file. Priority is 1. Note that
  // we are putting the operation inside the prelude so it is only executed once,
  // as opposed to executing for all entities. Note the semicolon at the end of
  // this rule, which is required when doing this.
  namespace {class file_timeStep005_1280810698m596 : public Loci::apply_rule< param<pair<int,real> > ,Loci::Maximum<pair<int,real> > >  {
#line 117 "timeStep.loci"
    Loci::const_param<real>  L_timeStep_n__ ; 
#line 117 "timeStep.loci"
    Loci::param<pair<int,real> >  L_dtTemp_n__ ; 
#line 117 "timeStep.loci"
public:
#line 117 "timeStep.loci"
    file_timeStep005_1280810698m596() {
#line 117 "timeStep.loci"
       name_store("dtTemp{n}",L_dtTemp_n__) ;
#line 117 "timeStep.loci"
       name_store("timeStep{n}",L_timeStep_n__) ;
#line 117 "timeStep.loci"
       input("timeStep{n}") ;
#line 117 "timeStep.loci"
       output("dtTemp{n}") ;
#line 117 "timeStep.loci"
    }
#line 117 "timeStep.loci"
    void prelude(const Loci::sequence &seq) { 
    join((*L_dtTemp_n__),pair<int,real>(1,(*L_timeStep_n__))) ;
  }    void compute(const Loci::sequence &seq) { 
#line 122 "timeStep.loci"
      prelude(seq) ;
#line 122 "timeStep.loci"
    }
#line 122 "timeStep.loci"
} ;
#line 122 "timeStep.loci"
Loci::register_rule<file_timeStep005_1280810698m596> register_file_timeStep005_1280810698m596 ;
#line 122 "timeStep.loci"
}
#line 122 "timeStep.loci"
// $type ncycle param<int> 
  // $type rampNumTimeSteps param<vector<int> > 
  // $type rampTimeStep param<vector<real> > 
  namespace {class file_timeStep006_1280810698m597 : public Loci::apply_rule< param<pair<int,real> > ,Loci::Maximum<pair<int,real> > >  {
#line 126 "timeStep.loci"
    Loci::const_param<real>  L_timeStep_n__ ; 
#line 126 "timeStep.loci"
    Loci::const_param<int>  L_ncycle_n__ ; 
#line 126 "timeStep.loci"
    Loci::const_param<vector<int> >  L_rampNumTimeSteps_n__ ; 
#line 126 "timeStep.loci"
    Loci::const_param<vector<real> >  L_rampTimeStep_n__ ; 
#line 126 "timeStep.loci"
    Loci::const_param<string>  L_timeStepRampType_n__ ; 
#line 126 "timeStep.loci"
    Loci::param<pair<int,real> >  L_dtTemp_n__ ; 
#line 126 "timeStep.loci"
public:
#line 126 "timeStep.loci"
    file_timeStep006_1280810698m597() {
#line 126 "timeStep.loci"
       name_store("dtTemp{n}",L_dtTemp_n__) ;
#line 126 "timeStep.loci"
       name_store("timeStep{n}",L_timeStep_n__) ;
#line 126 "timeStep.loci"
       name_store("ncycle{n}",L_ncycle_n__) ;
#line 126 "timeStep.loci"
       name_store("rampNumTimeSteps{n}",L_rampNumTimeSteps_n__) ;
#line 126 "timeStep.loci"
       name_store("rampTimeStep{n}",L_rampTimeStep_n__) ;
#line 126 "timeStep.loci"
       name_store("timeStepRampType{n}",L_timeStepRampType_n__) ;
#line 126 "timeStep.loci"
       input("ncycle{n},timeStep{n},rampNumTimeSteps{n},rampTimeStep{n},  timeStepRampType{n}") ;
#line 126 "timeStep.loci"
       output("dtTemp{n}") ;
#line 126 "timeStep.loci"
       constraint("noRestart") ;
#line 126 "timeStep.loci"
    }
#line 126 "timeStep.loci"
    void prelude(const Loci::sequence &seq) { 
    vector<int> interval ; interval.push_back(0) ;
    vector<double> step ; step.push_back(0.0) ;
    for(size_t i=0;i<(*L_rampNumTimeSteps_n__).size();++i){
      int endTimeStep=(*L_rampNumTimeSteps_n__)[i]+interval[i] ;
      interval.push_back(endTimeStep) ;
      step.push_back((*L_rampTimeStep_n__)[i]) ;
    }
    for(size_t i=1;i<interval.size();++i){
      if((*L_ncycle_n__)<interval[i]){
//      int left=interval[i-1],right=interval[i] ;
        if((*L_timeStepRampType_n__)=="constant"){
          join((*L_dtTemp_n__),pair<int,real>(2,step[i])) ;
        }else{
          cerr << "ERROR: Bad value for 'timeStepRampType'." << endl ;
          Loci::Abort() ;
        }
        return ;
      }
    }

    // Once out of the ramp range, use the normal timestep.
    join((*L_dtTemp_n__),pair<int,real>(2,(*L_timeStep_n__))) ;
  }    void compute(const Loci::sequence &seq) { 
#line 152 "timeStep.loci"
      prelude(seq) ;
#line 152 "timeStep.loci"
    }
#line 152 "timeStep.loci"
} ;
#line 152 "timeStep.loci"
Loci::register_rule<file_timeStep006_1280810698m597> register_file_timeStep006_1280810698m597 ;
#line 152 "timeStep.loci"
}
#line 152 "timeStep.loci"
namespace {class file_timeStep007_1280810698m598 : public Loci::singleton_rule {
#line 152 "timeStep.loci"
    Loci::const_param<pair<int,real> >  L_dtTemp_n__ ; 
#line 152 "timeStep.loci"
    Loci::param<real>  L_dt_n__ ; 
#line 152 "timeStep.loci"
public:
#line 152 "timeStep.loci"
    file_timeStep007_1280810698m598() {
#line 152 "timeStep.loci"
       name_store("dt{n}",L_dt_n__) ;
#line 152 "timeStep.loci"
       name_store("dtTemp{n}",L_dtTemp_n__) ;
#line 152 "timeStep.loci"
       input("dtTemp{n}") ;
#line 152 "timeStep.loci"
       output("dt{n}") ;
#line 152 "timeStep.loci"
    }
#line 152 "timeStep.loci"
    void compute(const Loci::sequence &seq) { 
//if(Loci::MPI_rank == 0) cout << "dt: " << $dt{n} << endl ;
    (*L_dt_n__)=(*L_dtTemp_n__).second ;
  }} ;
#line 155 "timeStep.loci"
Loci::register_rule<file_timeStep007_1280810698m598> register_file_timeStep007_1280810698m598 ;
#line 155 "timeStep.loci"
}
#line 155 "timeStep.loci"


}

namespace Loci {

  template<> struct data_schema_traits<streamUns::TimeStepRamp> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<streamUns::TimeStepRamp> Converter_Type ;
  } ;

}
