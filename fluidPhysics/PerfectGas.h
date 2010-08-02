#ifndef PerfectGas_H
#define PerfectGas_H

#include "eos.h"

#include <Tools/tools.h>
#include <Tools/debug.h>

#include <Tools/cptr.h>
using Loci::CPTR ;
using Loci::CPTR_type ;

#include <vector>

#include "chemistry_db.h"
// Thermally perfect gas mixture equation of state
//
// Note: This class is basically composed of two parts. First part
// builds up information such as internal energy, gas constant, mass
// fraction (etc.)  for every species, which is accomplised by
// invoking its member class efunc.  Second part gives the state
// information of mixture based on species data provided by the first
// part. The state information of mixture is not part of data members
// of class ThermallyPerfectEOS, but belongs to its member class
// State.

namespace fluidPhysics {  
  
  class energy_function : public CPTR_type {
  public:
    virtual ~energy_function() {}
    virtual double compute_T(const double *rhoi, double energy,
                             double *ei, double *cvi, double Tguess) const = 0 ;
    virtual double compute_T_from_h(const double *yi,double enthalpy,
                                    double *ei,double *cvi,double Tguess,const std::vector<double> &R) const = 0 ;
    virtual void compute_ei_cvi(double T, double *ei, double *cvi) const = 0  ;
    virtual void compute_dcvdT(double T, double *dcvdT) const = 0 ;
    virtual void compute_Omega(double T,double *O, double *Op) const = 0 ;

  } ;

  class ThermallyPerfectEOS: public equationOfState {
  private:
    // ostream and istream (Input/Output) operators needed to place
    // ThermallyPerfectEOS in a Loci parameter type.
    friend std::ostream &operator<<(std::ostream &s, const ThermallyPerfectEOS &eos) ;
    friend std::istream &operator>>(std::istream &s, ThermallyPerfectEOS &eos) ; 

    std::string thermo_model ;

    CPTR<energy_function> efunc ;
  public:
    ThermallyPerfectEOS() {num_species = 0 ; }
    ThermallyPerfectEOS(const ThermallyPerfectEOS &igm) ;
    ThermallyPerfectEOS &operator=(const ThermallyPerfectEOS &igm) ;
    //  ThermallyPerfectEOS(const species_db &species) {initialize(species) ; }
    ~ThermallyPerfectEOS() ;

    EoSPtr clone() const ;

    void initialize(const species_db &species,
                    std::string thermodynamic_model) ;

    // Calculate the state of mixture 
    State State_from_rho_e(const double *rhoi, double energy, double *mState,
                           const float *hint = 0 ) const;
    State State_from_rho_h(const double *rhoi, double enthalpy, double *mState,
                           const float *hint = 0 ) const;
    State State_from_rho_T(const double *rhoi, double T, double *mState,
                           const float *hint = 0) const ;
    State State_from_p_h(const double *yi,double p,double enthalpy,double *mState,
                         const float *hint=0) const ;
    State State_from_p_e(const double *yi,double p,double energy,double *mState,
                         const float *hint=0) const ;

    // Above two functions are commented in PerfectGas.cc
    // Function State_from_rho_T gets the state from mixture density and pressure
    State State_from_rho_p(const double *rhoi, double p, double *mState,
                           const float *hint = 0) const {
      double rRt = 0.0 ;
      for(int i=0;i<num_species;++i)
        rRt += rhoi[i]*R[i] ;
      return State_from_rho_T(rhoi,p/rRt,mState) ;
    }
    State State_from_mixture_p_T(const double *mixture, double p, double T,
                                 double *mState, const float *hint = 0) const {
      scratch_array<double> rho_temp(num_species) ;

      double Rt = 0.0 ;
      double mixt = 0.0 ;
      for(int i=0;i<num_species;++i) {
        Rt += mixture[i]*R[i] ;
        mixt += mixture[i] ;
      }

      Rt = Rt/mixt ;
      const double rho = p/(Rt*T) ;
      for(int i=0;i<num_species;++i)
        rho_temp[i] = mixture[i]*rho/mixt ;
      State s = State_from_rho_T(rho_temp,T,mState) ;

      return s ;
    }
  
    State State_from_rho_e(const double *rhoi, double energy,
                           const float *hint = 0) const {
      scratch_array<double> mState(mixtureStateSize()) ;
      State s = State_from_rho_e(rhoi,energy,mState,hint) ;
      return s;
    }
  
    State State_from_rho_h(const double *rhoi, double enthalpy,
                           const float *hint = 0) const {
      scratch_array<double> mState(mixtureStateSize()) ;
      State s = State_from_rho_h(rhoi,enthalpy,mState,hint) ;
      return s;
    }

    State State_from_rho_T(const double *rhoi, double T,
                           const float *hint = 0) const {
      scratch_array<double> mState(mixtureStateSize()) ;
      State s =  State_from_rho_T(rhoi,T,mState,hint) ;
      return s;
    }
    
    State State_from_rho_p(const double *rhoi, double p,
                           const float *hint = 0) const {
      scratch_array<double> mState(mixtureStateSize()) ;
      State s =  State_from_rho_p(rhoi,p,mState,hint) ;
      return s;
    }
    
    State State_from_mixture_p_T(const double *mixture, double p, double T,
                                 const float *hint = 0) const {
      scratch_array<double> mState(mixtureStateSize()) ;
      State s =  State_from_mixture_p_T(mixture,p,T,mState,hint) ; 
      return s;
    }


    double rho_from_p_mixture_T(const double p,const double *mixture, double T,
                                const float *hint = 0) const
    {
      double Rt = 0.0 ;
      for(int i=0;i<num_species;++i)
        Rt += mixture[i]*R[i] ;
      return p/(Rt*T) ;
    }
    int mixtureStateSize() const { return num_species*2; }
    int hintSize() const { return 1 ; }

    void getHint(float *hint, const State &s, const double *mixture) const {
      hint[0] = s.temperature() ;
    }
    
    double speciesEnthalpy(int i,const State &s, const double *ms) const {
      const double T = s.temperature() ;
      const double Ri = R[i] ;
      const double ei = ms[i] ;
      const double hi = ei + Ri*T ;
      return hi ;
    }
    void get_ei(const State &s, const double *ms, double *eis) const {
      int ns = num_species ;
      for(int i=0;i<ns;++i)
        eis[i] = ms[i] ;
    }

    void get_cvi(const State &s, const double *ms, double *cvis) const {
      int ns = num_species ;
      for(int i=0;i<ns;++i)
        cvis[i] = ms[i+ns] ;
    }
      
    void dTdri(double *DTDri,const State &s, const double *ms) const {
      for(int i=0;i<num_species;++i) 
        DTDri[i] = -s.temperature()*R[i]/s.get_rRt();
    }
    double dTdP(const State &s, const double *ms) const {
      return 1./s.get_rRt() ;
    }
    void dreidri(double *DreiDri,const State &s, const double *ms) const {
      const double factor = -s.temperature()*dreidT(s,ms)/s.get_rRt() ;
      const int ns = num_species ;
      for(int i=0;i<ns;++i) {
        double ei = ms[i] ;
        DreiDri[i] = ei + factor*R[i] ;
      }
    }
    double dreidT(const State &s, const double *ms) const {
      return s.rho_cvt() ;
    }
    double dreidP(const State &s, const double *ms) const {
      double gm1 = s.Gamma()-1.0 ;
      return 1./gm1 ;
    }
    void dcvdT(const State &s, double *dcvdT) const {
      efunc->compute_dcvdT(s.temperature(),dcvdT) ;
    }
    void getOmegas(double T, double *f, double *fp) const
    { efunc->compute_Omega(T,f,fp); }
    std::ostream &Print(std::ostream &s) const { s << *this ; return s ;}
    std::istream &Input(std::istream &s) { s >> *this ; return s ;} 
  } ;

  std::ostream & operator<<(std::ostream &s, const ThermallyPerfectEOS &eos) ;
  std::istream & operator>>(std::istream &s, ThermallyPerfectEOS &eos) ;
}
  
#endif

