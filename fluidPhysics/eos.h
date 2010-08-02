#ifndef EOS_H
#define EOS_H

#include <Tools/cptr.h>
#include "fluidConst.h"
#include "chemistry_db.h"
#include <vector>
#include <string>
#include <iostream>

#include "scratch_array.h"

namespace fluidPhysics {
  class equationOfState ;
  typedef Loci::CPTR<equationOfState> EoSPtr ;


  class equationOfState: public Loci::CPTR_type {
  protected:
    species_db sdb ;
    int num_species ;
    std::vector<std::string> namelist ;

    std::vector<double> Pref ;      // Reference Pressure (for thermo Kc calc.)
    std::vector<double> Tref ;      // Reference Temperature
    std::vector<double> R ;         // Species Gas Constant
    std::vector<double> m  ;        // Species Molecular Mass
    std::vector<double> mf ;        // Species default mass fraction
    std::vector<double> n ;         // Species rotational/translational energy
    // term
    std::vector<double> href ;      // reference enthalpy
    std::vector<double> sref ;      // reference entropy
  
    int num_elements ;            // number of elements in mixture

    // species_elements is a list of the elemental stochiometric coefficients for
    // each species.  For example, a mixture containing H,O,H2,O2,OH,and H2O
    // this data structure would look something like
    //         | H | O | H2 | O2 | OH | H2O |
    //       ==+===+===+====+====+====+=====|
    //    H  0 | 1 | 0 |  2 |  0 |  1 |   2 |
    //    O  1 | 0 | 1 |  0 |  2 |  1 |   1 |
    //       -------------------------------+

    std::vector<std::vector<int> > species_elements ;
  
  public:
    //class State is a nested class of equations of state, it represents the
    //thermodynamic state of mixture
    struct State {
      double T ;     // temperature of mixture
      double rho ;   // density of mixture
      double P ;
      double e ;
      double cv ;
      double cp ;
      double dfdp ;
      double dpdt ;
      double a2 ;
      State() {}
      State(double T_, double rho_, double rcvt, double rRt, double re_) :
        T(T_),rho(rho_) {
        P = rRt*T ;
        const double rrho = 1./rho ;
        e = re_*rrho ;
        cv = rcvt*rrho ;
        cp = (rcvt+rRt)*rrho ;
        const double gm1 = rRt/rcvt ;
        dfdp = 1./gm1 ;
        dpdt = rRt ;
        a2 = (gm1+1.)*P*rrho ;
      }
      State(double T_, double rho_, double P_, double e_, double cv_, double cp_,
            double dfdp_, double dpdt_, double a2_) :
        T(T_),rho(rho_),P(P_),e(e_),cv(cv_),cp(cp_),
        dfdp(dfdp_),dpdt(dpdt_),a2(a2_) {}
      // Part of following functions are the inspectors of class State */
      double get_rRt() const { return P/T ; }
      double gasConstant() const { return P/(rho*T) ; } 
      double temperature() const { return T ; } 
      double pressure() const    { return P ; } 
      double density() const     { return rho; } 
      double cvt() const         { return cv; } // specific heat
      double rho_cvt() const     { return rho*cv ; }
      double cpt() const         { return cp ; }
      double rho_cpt() const     { return rho*cp ; }
      double Gamma() const       { return cp/cv ;}
      double enthalpy() const    { return e + P/rho ; } 
      double energy() const      { return e ; } 
      double rho_energy() const  { return rho*e ; } // density-internal energy product
      double aSquared() const    { return a2 ; } 
      // Frozen speed of sound
      double soundSpeed() const  { 
        if(a2>0)
          return sqrt(a2) ;
        else
          return -sqrt(-a2) ;
      } 
      
      double dPdT() const {return dpdt ; }
      double dfdP() const {return dfdp ; }
    } ;
    // Calculate the state of mixture
    equationOfState() ;
    virtual ~equationOfState() ;
    void base_initialize(const species_db &species) ;

    virtual EoSPtr clone() const = 0 ;
    virtual void initialize(const species_db &species, std::string info) = 0 ;
  
    virtual State State_from_rho_e(const double *rhoi, double energy, double *mState,
                                   const float *hint = 0 ) const = 0 ;
    virtual State State_from_rho_h(const double *rhoi, double enthalpy, double *mState,
                                   const float *hint = 0 ) const = 0 ;
    virtual State State_from_rho_T(const double *rhoi, double T, double *mState,
                                   const float *hint = 0) const = 0 ;
    virtual State State_from_p_h(const double *yi,double p,double enthalpy,double *mState,
                                 const float *hint=0) const = 0 ;
    virtual State State_from_p_e(const double *yi,double p,double energy,double *mState,
                                 const float *hint=0) const = 0;

    // Above two functions are commented in PerfectGas.cc
    // Function State_from_rho_T gets the state from mixture density and pressure
    virtual State State_from_rho_p(const double *rhoi, double p, double *mState,
                                   const float *hint = 0) const = 0 ;
    virtual State State_from_mixture_p_T(const double *mixture, double p, double T,
                                         double *mState, const float *hint = 0) const =0 ;
    virtual State State_from_rho_e(const double *rhoi, double energy,
                                   const float *hint = 0) const = 0 ;
    virtual State State_from_rho_h(const double *rhoi, double enthalpy,
                                   const float *hint = 0) const = 0 ;
    virtual State State_from_rho_T(const double *rhoi, double T,
                                   const float *hint = 0) const = 0 ;
    virtual State State_from_rho_p(const double *rhoi, double p,
                                   const float *hint = 0) const = 0 ;
    
    virtual State State_from_mixture_p_T(const double *mixture, double p, double T,
                                         const float *hint = 0) const =0 ;

    virtual double rho_from_p_mixture_T(const double p,const double *mixture, double T,
                                        const float *hint = 0) const = 0 ;

    virtual int mixtureStateSize() const = 0 ;
    virtual int hintSize() const = 0;

    virtual void getHint(float *hint, const State &s, const double *mixture) const = 0 ;

    // Mixture property accessors
    virtual double speciesEnthalpy(int i,const State &s, const double *ms) const = 0;
    virtual void get_ei(const State &s, const double *ms, double *eis) const = 0 ;
    virtual void get_cvi(const State &s, const double *ms, double *cvis) const = 0;
    virtual void dTdri(double *DTDri,const State &s, const double *ms) const = 0 ;
    virtual double dTdP(const State &s, const double *ms) const = 0 ;
    virtual void dreidri(double *DreiDri,const State &s, const double *ms) const = 0 ;
    virtual double dreidP(const State &s, const double *ms) const = 0 ;
    //  virtual void dcvdT(const State &s, double *dcvdT) const = 0 ;
    virtual void getOmegas(double T, double *f, double *fp) const = 0 ;

    virtual std::ostream &Print(std::ostream &s) const = 0 ;
    virtual std::istream &Input(std::istream &s) = 0 ;
  
    // Standard utility routines
    int speciesIndex(std::string name) const {
      for(int i=0;i<num_species;++i)
        if(namelist[i] == name)
          return i ;
      return -1 ;
    }
    // Inspectors of EOS (get values of its data members)
    int numSpecies() const { return num_species ; }
    const std::vector<std::string> &speciesNames() const { return namelist ; }
    const std::vector<double> &getMixtureFractions() const { return mf ; }
    int numElements() const { return num_elements ; }
    const std::vector<int> &getSpeciesElements(int idx) const
    { return species_elements[idx] ; }
    double speciesMolecularMass(int i) const { return m[i] ; }
    double speciesR(int i) const { return R[i] ; }
    double speciesPref(int i) const { return Pref[i] ; }

  } ;



  // Equation of state delegation class
  class EOS {
  private:
    EoSPtr eos ;
  public:
    typedef equationOfState::State State ;
    EOS() {}
    EOS(EoSPtr cp) :eos(cp) {}

    EOS clone() { return EOS(eos->clone()) ; }
  
    void initialize(const species_db &species, std::string info) {
      eos->initialize(species,info) ;
    }

    State State_from_rho_e(const double *rhoi, double energy, double *mState,
                           const float *hint = 0 ) const {
      return eos->State_from_rho_e(rhoi,energy,mState,hint) ;
    }
    State State_from_rho_h(const double *rhoi, double enthalpy, double *mState,
                           const float *hint = 0 ) const {
      return eos->State_from_rho_h(rhoi,enthalpy,mState,hint) ;
    }
    State State_from_rho_T(const double *rhoi, double T, double *mState,
                           const float *hint = 0) const {
      return eos->State_from_rho_T(rhoi,T,mState,hint) ;
    }
    State State_from_p_h(const double *yi,double p,double enthalpy,double *mState,
                         const float *hint=0) const {
      return eos->State_from_p_h(yi,p,enthalpy,mState,hint) ;
    }

    State State_from_p_e(const double *yi,double p,double energy,double *mState,
                         const float *hint=0) const {
      return eos->State_from_p_e(yi,p,energy,mState,hint) ;
    }

    State State_from_rho_p(const double *rhoi, double p, double *mState,
                           const float *hint = 0) const {
      return eos->State_from_rho_p(rhoi,p,mState,hint) ;
    }
  
    State State_from_mixture_p_T(const double *mixture, double p, double T,
                                 double *mState, const float *hint = 0) const {
      return eos->State_from_mixture_p_T(mixture,p,T,mState,hint) ;
    }
  
    State State_from_rho_e(const double *rhoi, double energy,
                           const float *hint = 0) const {
      return eos->State_from_rho_e(rhoi,energy,hint) ;
    }
    State State_from_rho_h(const double *rhoi, double enthalpy,
                           const float *hint = 0) const {
      return eos->State_from_rho_h(rhoi,enthalpy,hint) ;
    }
    State State_from_rho_T(const double *rhoi, double T,
                           const float *hint = 0) const {
      return eos->State_from_rho_T(rhoi,T,hint) ;
    }
    State State_from_rho_p(const double *rhoi, double p,
                           const float *hint = 0) const {
      return eos->State_from_rho_p(rhoi,p,hint) ;
    }
    
    State State_from_mixture_p_T(const double *mixture, double p, double T,
                                 const float *hint = 0) const {
      return eos->State_from_mixture_p_T(mixture,p,T,hint) ;
    }

    double rho_from_p_mixture_T(const double p,const double *mixture, double T,
                                const float *hint = 0) const {
      return eos->rho_from_p_mixture_T(p,mixture,T,hint) ;
    }

    int mixtureStateSize() const { return eos->mixtureStateSize() ; }
    
    int hintSize() const { return eos->hintSize(); }

    void getHint(float *hint, const State &s, const double *ms) const {
      eos->getHint(hint,s,ms) ;
    }

    double speciesEnthalpy(int i,const State &s, const double *ms) const {
      return eos->speciesEnthalpy(i,s,ms) ;
    }
    void get_ei(const State &s, const double *ms, double *eis) const {
      eos->get_ei(s,ms,eis) ;
    }
    void get_cvi(const State &s, const double *ms, double *cvis) const {
      eos->get_cvi(s,ms,cvis) ;
    }
    void dTdri(double *DTDri,const State &s, const double *ms) const {
      eos->dTdri(DTDri,s,ms) ;
    }
    double dTdP(const State &s, const double *ms) const {
      return eos->dTdP(s,ms) ;
    }
    void dreidri(double *DreiDri,const State &s, const double *ms) const {
      eos->dreidri(DreiDri,s,ms) ;
    }
    double dreidP(const State &s, const double *ms) const {
      return eos->dreidP(s,ms) ;
    }
    void getOmegas(double T, double *f, double *fp) const {
      eos->getOmegas(T,f,fp) ;
    }

    // Standard utility routines
    int speciesIndex(std::string name) const {
      return eos->speciesIndex(name) ;
    }

    int numSpecies() const { return eos->numSpecies() ; }
    const std::vector<std::string> &speciesNames() const { return eos->speciesNames(); }
    const std::vector<double> &getMixtureFractions() const
    { return eos->getMixtureFractions() ; }
    int numElements() const { return eos->numElements() ; }
    const std::vector<int> &getSpeciesElements(int idx) const
    { return eos->getSpeciesElements(idx) ; }
    double speciesMolecularMass(int i) const
    { return eos->speciesMolecularMass(i) ; }
    double speciesR(int i) const { return eos->speciesR(i) ; }
    double speciesPref(int i) const { return eos->speciesPref(i) ; }

    std::ostream &Print(std::ostream &s) const { return eos->Print(s) ; }
    std::istream &Input(std::istream &s) { return eos->Input(s) ; }
  }  ;

  inline std::ostream &operator<<(std::ostream &s,const EOS &eos) {
    return eos.Print(s) ;
  }
  inline std::istream &operator>>(std::istream &s,EOS &eos) {
    return eos.Input(s) ;
  }

  inline std::ostream & operator<<(std::ostream &s, const equationOfState::State &eos_state)
  {
    s << eos_state.temperature() << ' '
      << eos_state.density() << ' '
      << eos_state.pressure() << ' '
      << eos_state.energy() << ' '
      << eos_state.cvt()  << ' '
      << eos_state.cpt() << ' '
      << eos_state.dfdP() << ' '
      << eos_state.dPdT() << ' '
      << eos_state.aSquared() 
      << std::endl ;
    return s;
  }

  inline std::istream & operator>>(std::istream &s, equationOfState::State &eos_state)
  {
    double T,rho,P, e,cv,cp,dfdp,dpdt,a2 ;
    s >> T >> rho >> P >> e >> cv >>cp >> dfdp >> dpdt >> a2 ;
    eos_state = equationOfState::State(T,rho,P,e,cv,cp,dfdp,dpdt,a2) ;
    return s ;
  }

  class EOSFactory {
    struct EOSInfo {
      int priority ;
      std::string name ;
      EoSPtr eos ;
      bool operator<(const EOSInfo &e) const { return priority<e.priority; }
    } ;
    std::vector<EOSInfo> eos_data ;
  public:
    void insertEOS(EoSPtr eos, std::string name, int priority) {
      EOSInfo entry ;
      entry.eos = eos ;
      entry.name = name ;
      entry.priority = priority ;
      eos_data.push_back(entry) ;
      std::sort(eos_data.begin(),eos_data.end()) ;
    }
    EoSPtr getEOS(std::string name) const {
      for(size_t i=0;i<eos_data.size();++i)
        if(eos_data[i].name == name)
          return eos_data[i].eos->clone() ;
      return EoSPtr(0) ;
    }
    EoSPtr getDefaultEOS() const {
      if(eos_data.size() >= 1)
        return eos_data.back().eos->clone() ;
      return EoSPtr(0) ;
    } 
  } ;
}    

namespace Loci {
  template<> struct data_schema_traits<fluidPhysics::equationOfState::State> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      fluidPhysics::equationOfState::State s ;
      CompoundDatatypeP ct = CompoundFactory(s) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::equationOfState::State,T) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::equationOfState::State,rho) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::equationOfState::State,P) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::equationOfState::State,e) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::equationOfState::State,cv) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::equationOfState::State,cp) ;
      return DatatypeP(ct) ;
    }
  } ;

  template<> struct data_schema_traits<fluidPhysics::EOS> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<fluidPhysics::EOS> Converter_Type ;
  } ;
}


#endif
