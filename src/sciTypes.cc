// Loci includes.
#include <Loci.h>

// StreamUns includes.
#include "const.h"
#include "name_var.h"
#include <readGrid.h>
#include "sciTypes.h"

namespace streamUns {

  vect3d get_vect3d(const options_list &ol,const char *vname,const char
  *units) {
    vect3d vec ;
    Loci::option_value_type ovt= ol.getOptionValueType(vname) ;
    if(ovt == Loci::REAL) {
      double v ;
      ol.getOption(vname,v) ;
      vec = vect3d(v,0,0) ;
    } else if(ol.getOptionValueType(vname) == Loci::UNIT_VALUE) {
      Loci::UNIT_type vu ;
      ol.getOption(vname,vu) ;
      if(!vu.is_compatible(units))
        cerr << "wrong type of units for vector " << vname
             << ": " << vu << endl ;
      else {
        double v ;
        v = vu.get_value_in(units) ;
        vec = vect3d(v,0,0) ;
      }
    } else if(ovt == Loci::LIST) {
      Loci::options_list::arg_list value_list ;
      ol.getOption(vname,value_list) ;
      if(value_list.size() != 3) {
        cerr << "error on reading '" << vname
             <<"': vector input must contain 3 terms"
             << endl ;
      }
      for(int i=0;i<3;++i)
        if(value_list[i].type_of() != Loci::REAL &&
           value_list[i].type_of() != Loci::UNIT_VALUE)
          cerr << "improper vector specification for '"
               << vname << "u' in fluidState"
               << endl ;
      double vecval[3] ;
      for(int i=0;i<3;++i) {
        if(value_list[i].type_of() == Loci::UNIT_VALUE) {
          Loci::UNIT_type vu ;
          value_list[i].get_value(vu) ;
          if(!vu.is_compatible(units))
            cerr << "wrong type of units for vector " << vname
                 << ": " << vu << endl ;
          vecval[i] = vu.get_value_in(units) ;
        } else
          value_list[i].get_value(vecval[i]) ;
      }
      vec.x = vecval[0] ;
      vec.y = vecval[1] ;
      vec.z = vecval[2] ;
    } else if(ovt == Loci::FUNCTION) {
      string name ;
      Loci::options_list::arg_list value_list ;
      ol.getOption(vname,name,value_list) ;
      if(name != "polar") {
        cerr << "don't know coordinate function '" << name
             <<"', defaulting to polar" << endl ;
      }
      if(value_list.size() != 3) {
        cerr << "error on reading '"
             << vname << "': vector input must contain 3 terms"
             << endl ;
      }
      for(int i=0;i<3;++i)
        if(value_list[i].type_of() != Loci::REAL &&
           value_list[i].type_of() != Loci::UNIT_VALUE)
          cerr << "improper vector specification for '"
               << vname << "' in fluidState"
               << endl ;
      real r=1 ,theta=0 ,eta=0 ;
      real conv = M_PI/180.0 ;
      if(value_list[0].type_of() == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        value_list[0].get_value(vu) ;
        if(!vu.is_compatible(units))
          cerr << "wrong type of units for vector " << vname
               << ": " << vu << endl ;
        r = vu.get_value_in(units) ;
      } else
        value_list[0].get_value(r) ;
      if(value_list[1].type_of() == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        value_list[1].get_value(vu) ;
        if(!vu.is_compatible("radians"))
          cerr << "wrong type of units for vector " << vname
               << ": " << vu << endl ;
        theta = vu.get_value_in("radians") ;
      } else {
        value_list[1].get_value(theta) ;
        theta *= conv  ;
      }
      if(value_list[2].type_of() == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        value_list[2].get_value(vu) ;
        if(!vu.is_compatible("radians"))
          cerr << "wrong type of units for vector " << vname
               << ": " << vu << endl ;
        eta = vu.get_value_in("radians") ;
      } else {
        value_list[2].get_value(eta) ;
        eta *= conv  ;
      }

      vec.x = r*cos(theta)*cos(eta) ;
      vec.y = r*sin(theta)*cos(eta) ;
      vec.z = r*sin(eta) ;
    } else {
      cerr << "incorrect type for 'u' in fluidState" << endl ;
    }
    return vec ;
  }

  // Retrieves the value for a vect3d option.
  void get_vect3dOption(const options_list &ol,string vname,string units,
  vect3d &vec,real Lref) {
    Loci::option_value_type ovt= ol.getOptionValueType(vname) ;
    if(ovt == Loci::REAL) {
      double v ;
      ol.getOption(vname,v) ;
      vec = vect3d(v*Lref,0,0) ;
    } else if(ol.getOptionValueType(vname) == Loci::UNIT_VALUE) {
      Loci::UNIT_type vu ;
      ol.getOption(vname,vu) ;
      if(!vu.is_compatible(units)) {
        std::cerr << "wrong type of units for vector " << vname
                  << ": " << vu << std::endl ;
        Loci::Abort() ;
      } else {
        double v ;
        v = vu.get_value_in(units) ;
        vec = vect3d(v,0,0) ;
      }
    } else if(ovt == Loci::LIST) {
      Loci::options_list::arg_list value_list ;
      ol.getOption(vname,value_list) ;
      if(value_list.size() != 3) {
        std::cerr << "error on reading '" << vname
                  <<"': vector input must contain 3 terms"
                  << std::endl ;
        Loci::Abort() ;
      }
      for(int i=0;i<3;++i)
        if(value_list[i].type_of() != Loci::REAL &&
           value_list[i].type_of() != Loci::UNIT_VALUE) {
          std::cerr << "improper vector specification for '"
                    << vname << std::endl ;
          Loci::Abort() ;
        }
      double vecval[3] ;
      for(int i=0;i<3;++i) {
        if(value_list[i].type_of() == Loci::UNIT_VALUE) {
          Loci::UNIT_type vu ;
          value_list[i].get_value(vu) ;
          if(!vu.is_compatible(units)) {
            std::cerr << "wrong type of units for vector " << vname
                      << ": " << vu << std::endl ;
            Loci::Abort() ;
          }
          vecval[i] = vu.get_value_in(units) ;
        } else {
          value_list[i].get_value(vecval[i]) ;
          vecval[i] *= Lref ;
        }
      }
      vec.x = vecval[0] ;
      vec.y = vecval[1] ;
      vec.z = vecval[2] ;
    } else if(ovt == Loci::FUNCTION) {
      string name ;
      Loci::options_list::arg_list value_list ;
      ol.getOption(vname,name,value_list) ;
      if(name != "polar") {
        std::cerr << "don't know coordinate function '" << name
                  <<"', defaulting to polar" << std::endl ;
        Loci::Abort() ;
      }
      if(value_list.size() != 3) {
        std::cerr << "error on reading '"
                  << vname << "': vector input must contain 3 terms"
                  << std::endl ;
        Loci::Abort() ;
      }
      for(int i=0;i<3;++i)
        if(value_list[i].type_of() != Loci::REAL &&
           value_list[i].type_of() != Loci::UNIT_VALUE) {
          std::cerr << "improper vector specification for '"
                    << vname << std::endl ;
          Loci::Abort() ;
        }
      real r=1 ,theta=0 ,eta=0 ;
      real conv = M_PI/180.0 ;
      if(value_list[0].type_of() == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        value_list[0].get_value(vu) ;
        if(!vu.is_compatible(units)) {
          std::cerr << "wrong type of units for vector " << vname
                    << ": " << vu << std::endl ;
          Loci::Abort() ;
        }
        r = vu.get_value_in(units) ;
      } else {
        value_list[0].get_value(r) ;
        r *= Lref ;
      }
      if(value_list[1].type_of() == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        value_list[1].get_value(vu) ;
        if(!vu.is_compatible("radians")) {
          std::cerr << "wrong type of units for vector " << vname
                    << ": " << vu << std::endl ;
          Loci::Abort() ;
        }
        theta = vu.get_value_in("radians") ;
      } else {
        value_list[1].get_value(theta) ;
        theta *= conv  ;
      }
      if(value_list[2].type_of() == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        value_list[2].get_value(vu) ;
        if(!vu.is_compatible("radians")) {
          std::cerr << "wrong type of units for vector " << vname
                    << ": " << vu << std::endl ;
          Loci::Abort() ;
        }
        eta = vu.get_value_in("radians") ;
      } else {
        value_list[2].get_value(eta) ;
        eta *= conv  ;
      }

      vec.x = r*cos(theta)*cos(eta) ;
      vec.y = r*sin(theta)*cos(eta) ;
      vec.z = r*sin(eta) ;
    } else {
      std::cerr << "unable to get vector type!" << std::endl ;
      Loci::Abort() ;
    }
  }

}
