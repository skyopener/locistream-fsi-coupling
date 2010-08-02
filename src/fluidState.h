#ifndef FLUIDSTATE_H
#define FLUIDSTATE_H

#include <map>
#include <string>
#include <iostream>
#include "eos.h"
#include "reaction.h"
using fluidPhysics::EOS ;
using fluidPhysics::reaction ;
#include "qvi.h"
#include "sciTypes.h"

namespace streamUns {

  class fluidState_ol : public options_list {
  public:
    fluidState_ol() :
      options_list("rho:T:p:M:u:equilibrium:mixture:k:e:w:R:nu_t:nu_tm:tmu") {}
  } ;


  
  class fluidState {
    std::map<std::string,real> mixture;
    real p,T ;
    real rho ;
    vect3d u ;
    real k ;
    real e ;
    real w ;
    real R ;
    real nu_t ;
    real nu_tm ;
    real tmu ;	
    bool use_mach_velocity ;
    bool set_density_pressure ;
    bool set_density_temperature ;
    bool set_pressure_temperature ;
    bool set_k ;
    bool set_e ;
    bool set_w ;
    bool set_R ;
    bool set_nu_t ;
    bool set_nu_tm ;
    bool set_tmu ;	
    bool equilibrium_mixture ;
    bool values_defined ;
    void setEquilibriumState(real *q, const conservativeVectorInfo &qvi,
                             const EOS &eos, const reaction &r) const ;
    void assign(const fluidState &f) {
      if(f.values_defined) {
        mixture = f.mixture ;
        p = f.p ;
        T = f.T ;
        rho = f.rho ;
        u = f.u ;
        use_mach_velocity = f.use_mach_velocity ;
        set_density_pressure = f.set_density_pressure ;
        set_density_temperature = f.set_density_temperature ;
        set_pressure_temperature = f.set_pressure_temperature ;
        set_k = f.set_k ;
        set_e = f.set_e ;
        set_w = f.set_w ;
        set_R = f.set_R ;
        set_nu_t = f.set_nu_t ;
        set_nu_tm = f.set_nu_tm ;
        set_tmu = f.set_tmu ;
        equilibrium_mixture = f.equilibrium_mixture ;
        k = f.k ;
        e = f.e ;
        w = f.w ;
        R = f.R ;
        nu_t = f.nu_t ;
        nu_tm = f.nu_tm ;
        tmu = f.tmu ;
        values_defined = true ;
      } else {
        values_defined = false ;
        use_mach_velocity = false ;
        set_density_pressure = false ;
        set_density_temperature = false ;
        set_pressure_temperature = false ;
        set_k = false ;
        set_e = false ;
        set_w = false ;
        set_R = false ;
        set_nu_t = false ;
        set_nu_tm = false ;
        set_tmu = false ;
        equilibrium_mixture = false ;
        rho = 1.0 ;
        p = 1.0 ;
        T = 300.0 ;
        u = vect3d(0,0,0) ;
        k = 0 ;
        e = 0 ;
        w = 0 ;
        R = 0 ;
        nu_t = 0 ;
        nu_tm = 0 ;
        tmu = 0 ;
      }
    }
  public:
    fluidState() ;
    fluidState(const fluidState &f) { values_defined = false ; assign(f) ; }
    fluidState &operator=(const fluidState &f) { assign(f); return *this ; }

    int pack_buf_size() ;
    void pack_buf(real *buf, int sz) ;
    void unpack_buf(real *buf, int sz) ;

    bool defined() const { return values_defined ; } ;
    bool get_k_omega(real &k, real &w) const ;
    bool get_k_epsilon(real &k, real &e) const ;
    bool get_R(real &R) const ;
    bool get_nu_t(real &nu_t) const ;
    bool get_nu_tm(real &nu_tm) const ;
    void setState(real *q, const conservativeVectorInfo &qvi,
                  const EOS &eos, const reaction &r) const ;
    void setPrimitive(real *qp, real Pambient,
                      const conservativeVectorInfo &qvi,
                      const EOS &eos, const reaction &r) const ;

    void Input(const options_list &ol) ;
    void Input(const options_list::arg_list &l) {
      fluidState_ol finput ;
      finput.Input(l) ;
      Input(finput) ;
    }
    std::istream &Input(std::istream &s) {
      fluidState_ol finput ;
      s >> finput ;
      Input(finput) ;
      return s ;
    }
    
    std::ostream &Print(std::ostream &s) const ;

 
    // added for the convenience of initializing fluidState for inflow bc.
    fluidState(const vect3d uu, const real pp, const real tt) :
      p(pp),T(tt),u(uu) {
      values_defined = true ;
      use_mach_velocity = false ;
      set_pressure_temperature = true ;
      equilibrium_mixture = false ; 
      set_density_pressure = false ;
      set_density_temperature = false ;
      set_k = false ;
      set_e = false ;
      set_w = false ;
      set_R = false ;
      set_nu_t = false ;
      set_nu_tm = false ;
      set_tmu = false ;
    }
  
    fluidState(const conservativeVectorInfo &qvi, const real *mf, 
               const vect3d uu, const real pp, const real tt) :
      p(pp),T(tt),u(uu) {
      values_defined = true ;
      use_mach_velocity = false ;
      set_pressure_temperature = true ;
      equilibrium_mixture = false ;
      set_density_pressure = false ;
      set_density_temperature = false ;
      int ns = qvi.numSpecies() ;

      real total = 0.0 ;
      for(int i=0;i<ns;++i) 
        total += mf[i] ;
      for(int i=0;i<ns;++i) 
        mixture[qvi.speciesName(i)] = mf[i]/total ;
      set_k = false ;
      set_e = false ;
      set_w = false ;
      set_R = false ;
      set_nu_t = false ;
      set_nu_tm = false ;
      set_tmu = false ;
    }

    fluidState(const conservativeVectorInfo &qvi, const real *mf, 
               const vect3d uu, const real pp, const real tt, 
	       bool equilibrium ) :
      p(pp),T(tt),u(uu) {
      values_defined = true ;
      use_mach_velocity = false ;
      set_pressure_temperature = true ;
      equilibrium_mixture = equilibrium ;
      set_density_pressure = false ;
      set_density_temperature = false ;
      int ns = qvi.numSpecies() ;

      real total = 0.0 ;
      for(int i=0;i<ns;++i) 
        total += mf[i] ;
      for(int i=0;i<ns;++i) 
	mixture[qvi.speciesName(i)] = mf[i]/total ;
      set_k = false ;
      set_e = false ;
      set_w = false ;
      set_R = false ;
      set_nu_t = false ;
      set_nu_tm = false ;
      set_tmu = false ;
    }

    real getDensity() const {return rho ;} 
    real getPressure() const {return p ;}
    real getTemperature() const {return T ;}
    vect3d getVelocity() const { return u ;}
    real getTenergy() const {return k ;} 
    real getTdissp() const {return e ;}
    real getTR() const {return R ;}
    real getTnu_t() const {return nu_t ;}
    real getTnu_tm() const {return nu_tm ;}
    bool mixture_state() const{ return equilibrium_mixture ; }
    std::map<std::string,real> getMixture() const { return mixture ; }

  } ;

  inline std::ostream & operator<<(std::ostream &s, const fluidState &fs)
  { return fs.Print(s) ; }
  inline std::istream & operator>>(std::istream &s, fluidState &fs)
  { return fs.Input(s) ; }

  vect3d get_vect3d(const options_list &ol,const char *vname, const char *units) ;
}

namespace Loci {
  class fluidStateConverter {
    streamUns::fluidState &ref ;
  public:
    fluidStateConverter(streamUns::fluidState &iref) : ref(iref) {}
    int getSize() const {
      return ref.pack_buf_size() ;
    }
    void getState(streamUns::real *buf, int &size) {
      size = getSize() ;
      ref.pack_buf(buf,size) ;
    }
    void setState(streamUns::real *buf, int size) {
      ref.unpack_buf(buf,size) ;
    }
  } ;
  
  template<> struct data_schema_traits<streamUns::fluidState> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef streamUns::real Converter_Base_Type ;
    typedef fluidStateConverter Converter_Type ;
  } ;
}


#endif
