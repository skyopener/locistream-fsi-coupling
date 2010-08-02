//-----------------------------------------------------------------------------
// Description: This file contains rules for an ignition point, which is used
//   to initiate combustion.
//
// Author: Ed Luke
// Modified by: Siddharth Thakur and Jeff Wright
//              Streamline Numerics Inc.
//              Gainesville, Florida 32609
//-----------------------------------------------------------------------------
                                                                                
// Loci includes.
#include <Loci.h>

// StreamUns includes.
#include "sciTypes.h"
#include "varsFileInputs.h"

namespace streamUns {

//-----------------------------------------------------------------------------
// Rules to process ignition options from the .vars file.

  // Creates the ignition parameters.
  class IgnitionParameters : public singleton_rule {
    private:
      const_param<IgnitionOptions> ignition ;
      const_param<real> Lref ;
      param<real> ignition_Tmax ;
      param<real> ignition_heat ;
      param<vect3d> ignition_pos ;
      param<vect3d> ignition_axis ;
      param<real> ignition_dz,ignition_rmin,ignition_rmax ;
      param<vect3d> ignition_delta ;
      param<int> ignition_nstart ;
      param<int> ignition_nstop ;
    public:
                                                                                
      // Define input and output.
      IgnitionParameters() {
        name_store("ignition",ignition) ;
        name_store("Lref",Lref) ;
        name_store("ignition_Tmax",ignition_Tmax) ;
        name_store("ignition_heat",ignition_heat) ;
        name_store("ignition_pos",ignition_pos) ;
        name_store("ignition_axis",ignition_axis) ;
        name_store("ignition_dz",ignition_dz) ;
        name_store("ignition_rmin",ignition_rmin) ;
        name_store("ignition_rmax",ignition_rmax) ;
        name_store("ignition_delta",ignition_delta) ;
        name_store("ignition_nstart",ignition_nstart) ;
        name_store("ignition_nstop",ignition_nstop) ;
        input("ignition,Lref") ;
        output("ignition_Tmax,ignition_heat,ignition_pos,ignition_delta") ;
        output("ignition_nstart,ignition_nstop") ;
        output("ignition_axis,ignition_dz,ignition_rmin,ignition_rmax") ;
      }
                                                                                
      // Set up the parameters.
      virtual void compute(const sequence& seq) {

        if((*ignition).optionExists("position")){

          // Maximum temperature.
          real Tmax=4000 ;
          Loci::option_value_type ovt=ignition->getOptionValueType("Tmax") ;
          if(ovt==Loci::REAL){ ignition->getOption("Tmax",Tmax) ; }
          if(ovt==Loci::UNIT_VALUE){
            Loci::UNIT_type vu ; ignition->getOption("Tmax",vu) ;
            if(!vu.is_compatible("K")) {
              cerr << "Wrong type of units for ignition Tmax." << endl ;
            }else{
              Tmax=vu.get_value_in("K") ;
            }
          }
          *ignition_Tmax=Tmax ;
                                                                                
          // Heating rate.
          real heat = 1e9;
          ovt=ignition->getOptionValueType("heat") ;
          if(ovt==Loci::REAL){ ignition->getOption("heat",heat) ; }
          if(ovt==Loci::UNIT_VALUE){
            Loci::UNIT_type vu ; ignition->getOption("heat",vu) ;
            if(!vu.is_compatible("watts/m^3")) {
              cerr << "Wrong type of units ignition heat option." << endl ;
            }else{
             heat=vu.get_value_in("watts/m^3") ;
            }
          }
          *ignition_heat=heat ;

          ovt=ignition->getOptionValueType("position") ;
          vect3d val(0.0,0.0,0.0) ;
          if(ovt==Loci::LIST){
            Loci::options_list::arg_list valueList ;
            ignition->getOption("position",valueList) ;
            if(valueList.size()!=3){
              cerr << "Error on reading ignition position: vector input must "
                << "contain 3 terms." << endl ;
            }
            for(int j=0;j<3;++j)
              if(valueList[j].type_of()!=Loci::REAL){
                cerr << "Improper vector specification for ignition position."
                  << endl ;
              }
            double vecVal[3] ;
            for(int j=0;j<3;++j){ valueList[j].get_value(vecVal[j]) ; }
            val=vect3d(vecVal[0],vecVal[1],vecVal[2]) ;
          }else{
            cerr << "Improper coordinate for ignition position." << endl ;
          }
          *ignition_pos=*Lref*val ;

          // Unit axis for annular ignition zone.
          if((*ignition).optionExists("axis")){
            ovt=ignition->getOptionValueType("axis") ;
            val=vect3d(0.0,0.0,0.0) ;
            if(ovt==Loci::LIST){
              Loci::options_list::arg_list valueList ;
              ignition->getOption("axis",valueList) ;
              if(valueList.size()!=3){
                cerr << "Error on reading ignition axis: vector input must "
                  << "contain 3 terms." << endl ;
              }
              for(int j=0;j<3;++j)
                if(valueList[j].type_of()!=Loci::REAL){
                  cerr << "Improper vector specification for ignition axis."
                    << endl ;
                }
              double vecVal[3] ;
              for(int j=0;j<3;++j){ valueList[j].get_value(vecVal[j]) ; }
              val=vect3d(vecVal[0],vecVal[1],vecVal[2]) ;
            }else{
              cerr << "Improper coordinate for ignition axis." << endl ;
            }
            *ignition_axis=val ;
          }

          // Axial extent for annular ignition zone. Include points that lie
          // within +dz and -dz of ignition_pos. Must initialize this even
          // though the option may not be present since we use a non-zero value
          // as a test later in whether we have the standard ignition zone or
          // an annular ignition zone.
          *ignition_dz=0.0 ;
          if((*ignition).optionExists("dz")){
            real dz=0.0 ; ovt=ignition->getOptionValueType("dz") ;
            if(ovt==Loci::REAL){ ignition->getOption("dz",dz) ; }
            if(ovt==Loci::UNIT_VALUE){
              Loci::UNIT_type vu ; ignition->getOption("dz",vu) ;
              if(!vu.is_compatible("m")) {
                cerr << "Wrong type of units for ignition dz." << endl ;
              }else{
                dz=vu.get_value_in("m") ;
              }
            }
            *ignition_dz=dz ;
          }

          // Minimum radial extent for annular ignition zone. Include points that lie
          // more than the minimum radius from ignition_pos.
          if((*ignition).optionExists("rmin")){
            real rmin=0.0 ; ovt=ignition->getOptionValueType("rmin") ;
            if(ovt==Loci::REAL){ ignition->getOption("rmin",rmin) ; }
            if(ovt==Loci::UNIT_VALUE){
              Loci::UNIT_type vu ; ignition->getOption("rmin",vu) ;
              if(!vu.is_compatible("m")) {
                cerr << "Wrong type of units for ignition rmin." << endl ;
              }else{
                rmin=vu.get_value_in("m") ;
              }
            }
            *ignition_rmin=rmin ;
          }

          // Maximum radial extent for annular ignition zone. Include points that lie
          // less than the minimum radius from ignition_pos.
          if((*ignition).optionExists("rmax")){
            real rmax=0.0 ; ovt=ignition->getOptionValueType("rmax") ;
            if(ovt==Loci::REAL){ ignition->getOption("rmax",rmax) ; }
            if(ovt==Loci::UNIT_VALUE){
              Loci::UNIT_type vu ; ignition->getOption("rmax",vu) ;
              if(!vu.is_compatible("m")) {
                cerr << "Wrong type of units for ignition rmax." << endl ;
              }else{
                rmax=vu.get_value_in("m") ;
              }
            }
            *ignition_rmax=rmax ;
          }

          // Delta.
          if((*ignition).optionExists("delta")){
            val=vect3d(0,0,0) ; ovt=ignition->getOptionValueType("delta") ;
            if(ovt==Loci::LIST){
              Loci::options_list::arg_list valueList ;
              ignition->getOption("delta",valueList) ;
              if(valueList.size()!=3) {
                cerr << "Error on reading delta: vector input must contain 3 "
                  << "terms." << endl ;
              }
              for(int j=0;j<3;++j)
                if(valueList[j].type_of()!=Loci::REAL){
                  cerr << "Improper vector specification for ignition delta."
                    << endl ;
                }
              double vecVal[3] ;
              for(int j=0;j<3;++j){ valueList[j].get_value(vecVal[j]) ; }
              val=vect3d(vecVal[0],vecVal[1],vecVal[2]) ;
            }else if(ovt == Loci::REAL){
              real v1 ; ignition->getOption("delta",v1) ; val=vect3d(v1,v1,v1) ;
            }else{
              cerr << "Improper coordinate for ignition delta." << endl ;
            }
            *ignition_delta=*Lref*val ;
          }

          // Starting and ending time step.
          double v1=0 ;
          if(ignition->optionExists("nstart")){
            ignition->getOption("nstart",v1) ; *ignition_nstart=int(v1) ;
          }
          v1=0 ;
          if(ignition->optionExists("nstop")){
            ignition->getOption("nstop",v1) ; *ignition_nstop = int(v1) ;
          }
        }
      }
  } ;
                                                                                
  register_rule<IgnitionParameters> registerIgnitionParamters ;

//-----------------------------------------------------------------------------
// Rule to add ignition source to energy equation. 

  class IgnitionSrc : public apply_rule<storeVec<real>,Loci::Summation
  <Vect<real> > > {
    private:
      const_param<real> ignition_heat ;
      const_param<real> ignition_Tmax ;
      const_param<vect3d> ignition_pos ;
      const_param<vect3d> ignition_delta ;
      const_param<vect3d> ignition_axis ;
      const_param<real> ignition_dz,ignition_rmin,ignition_rmax ;
      const_param<int> ignition_nstart ;
      const_param<int> ignition_nstop ;
      const_param<int> nCycle ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      const_store<vect3d> cellCenter ;
      const_store<real> temperature ;
      store<real> hSourceTerm ;
    public:

      // Define input and output.
      IgnitionSrc() {
        name_store("ignition_Tmax{n,it}",ignition_Tmax) ;
        name_store("ignition_heat{n,it}",ignition_heat) ;
        name_store("ignition_pos{n,it}",ignition_pos) ;
        name_store("ignition_delta{n,it}",ignition_delta) ;
        name_store("ignition_axis{n,it}",ignition_axis) ;
        name_store("ignition_dz{n,it}",ignition_dz) ;
        name_store("ignition_rmin{n,it}",ignition_rmin) ;
        name_store("ignition_rmax{n,it}",ignition_rmax) ;
        name_store("ignition_nstart{n,it}",ignition_nstart) ;
        name_store("ignition_nstop{n,it}",ignition_nstop) ;
        name_store("ncycle{n}",nCycle) ;
        name_store("vol{n,it}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("cellcenter{n,it}",cellCenter) ;
        name_store("temperature{n,it}",temperature) ;
        name_store("hSourceTerm{n,it}",hSourceTerm) ;
        input("ignition_heat{n,it},ignition_Tmax{n,it},ignition_pos{n,it}") ;
        input("ignition_delta{n,it},ignition_nstart{n,it}") ;
        input("ignition_nstop{n,it},ncycle{n}") ;
        input("ignition_axis{n,it},ignition_dz{n,it},ignition_rmin{n,it},ignition_rmax{n,it}") ;
        input("vol{n,it},cellcenter{n,it},temperature{n,it},cellRadius{n,it}") ;
        output("hSourceTerm{n,it}") ;
      }

      // Add volumetric heat at prescribed location.
      void calculate (Entity cell) {
        if(temperature[cell]<*ignition_Tmax) {
          if(ignition_dz[cell]>1.0e-30){
            // Annular ignition zone.
            vect3d p=cellCenter[cell]-ignition_pos[cell] ;
            real dz=abs(dot(p,ignition_axis[cell])) ;
            if(dz<ignition_dz[cell]){
              vect3d r=p-dot(p,ignition_axis[cell])*ignition_axis[cell] ;
              real rMag=norm(r) ;
              if(rMag>ignition_rmin[cell] && rMag<ignition_rmax[cell]){
                hSourceTerm[cell]+=(*ignition_heat)*vol[cell]*cellRadius[cell] ;
              }
            }
          }else{
            // Standard ignition zone.
            real ref_length=pow(vol[cell],0.33333333) ;
            vect3d dp=cellCenter[cell]-ignition_pos[cell] ;
            if(fabs(dp.x)<(ignition_delta[cell].x+ref_length) &&
            fabs(dp.y)<(ignition_delta[cell].y+ref_length) &&
            fabs(dp.z)<(ignition_delta[cell].z+ref_length)) {
              hSourceTerm[cell]+=(*ignition_heat)*vol[cell]*cellRadius[cell] ;
            }
          }
        }
      }

      // Loop through all cells.
      void compute(const sequence &seq) {
        if(*nCycle >= *ignition_nstart && *nCycle <= *ignition_nstop){
          if(Loci::MPI_rank==0) cout << "    INFO: Ignition active" << endl ;
          do_loop(seq,this) ;
        }
      }
  } ;

  register_rule<IgnitionSrc> registerIgnitionSrc ;
}
