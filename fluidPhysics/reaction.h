#ifndef REACTION_H
#define REACTION_H


#include "eos.h"
#include "chemistry_db.h"

#include <vector>
#include <map>

using std::vector ;

namespace fluidPhysics {
  // This is the Gibbs Free energy minimization Kc computation compiler class
  // This class compiles a set of Kc calculation requests and compiles a
  // schedule for computing the collection of Kc's while factoring out common
  // computations
  class thermoKc {
    entitySet speciesSet ;
    
    struct spi {
      double PrefR ;
    } ;
    std::vector<spi> spiv ;
    struct si { // species info
      int species ;
      double nu_diff ;
      bool operator<(const si &s) const {
        if(species<s.species)
          return true ;
        else if(species == s.species && nu_diff < s.nu_diff)
          return true ;
        return false ;
      } 
      bool operator==(const si &s) const
      { return species==s.species && fabs(nu_diff-s.nu_diff)<EPSILON ; }
      bool operator!=(const si &s) const { return !operator==(s); }
    } ;

    struct ri { // reaction info
      std::vector<si> pos ;  //size=num_reaction for each species
      std::vector<si> neg ;
      ri() {}
      ri(const EOS &eos, const reaction_expression &r) ;
      bool operator==(const ri &r) ;
      bool operator<(const ri &r) const {
        if(pos < r.pos)
          return true ;
        else if(pos == r.pos && neg < r.neg)
          return true ;
        return false ;
      }
      bool forwardonly ;  
    } ;
    friend struct ri ;
    friend class KcFunc ;
    std::vector<ri> reactinfo ;
    std::vector<int> res_loc ; // result location

  public:
    thermoKc() {}
    ~thermoKc() {}
    void initialize(const EOS &eos) ;
    int compile_KcFunc(const EOS &eos, const reaction_expression &r,
                       int &f_val) ;
    void calculate(const EOS &eos, double T, std::pair<double,double> *f_val) const ;
  } ;

  // This is the curve fit function compiler class.  It compiles a set of
  // of curve fits (for the moment only Arrhenius curve fits are compiled)
  // factoring out common computations.
  class curveFitFunctions {
    // The func_rep class is a representation of the curve fit function.
    struct func_rep {
      int etai,thetai ;
      double C ;
      func_rep(double Cr,int eta,int theta) : etai(eta),thetai(theta),C(Cr) {}
      bool operator<(const func_rep &t) const {
        if(etai < t.etai)
          return true ;
        else if(etai == t.etai && thetai < t.thetai)
          return true ;
        else if(thetai == t.thetai && C < t.C)
          return true ;
        return false ;
      }      
    } ;

    // These maps are used to find identical eta, and thetas
    std::map<double,int> eta_map,theta_map ;
    // Theta and etas in the database
    std::vector<double>  eta, theta ;  //Note: the elements in vector eta and theta
    //                             //are unrepeatable if you refer to functions 
    //                             //etaindex and thindex appear later

    // sizes of temporary computation space for power and exponential functions
    int num_power, num_exponential ;

    // this is a map of function representations used to find redundant
    // curve fit requests
    std::map<func_rep,int> func_map ;
    std::vector<func_rep> func_list ;
    // a list of pointers for each func_rep in func_list to the place where
    // the resulting function will be placed
    std::vector<int> res_loc ; 

    // Get an integer index corresponding to T^{eta}
    int etaindex(double etaval) {
      std::map<double,int>::iterator cmi ;
      if ((cmi = eta_map.find(etaval)) != eta_map.end()) //if etaval prexists in 
        //                                                   //eta_map
        return cmi->second ;
      int idx = eta.size() ;
      eta.push_back(etaval) ;
      eta_map[etaval] = idx ;
      num_power++ ;
      return idx ;
    }
    // Get an integer index corresponding to temperature theta
    int thindex(double thetaval) {
      std::map<double,int>::iterator cmi ;
      if ((cmi = theta_map.find(thetaval)) != theta_map.end())
        return cmi->second ;
      int idx = theta.size() ;
      theta.push_back(thetaval) ;
      theta_map[thetaval] = idx ;
      num_exponential++ ;
      return idx ;
    }

    // Get an integer index for an Arrhenius curve fit function
    // This index is where calculate will place its results.
    // This function allocates the necessary space for this result
    // in fval.
    int getArrheniusidx(double Cvar,double etavar, double thetavar,
                        int &fval) ;

  public:
    curveFitFunctions() ;

    int compile_Arrhenius(double Ci, double etai, double thetai,
                          int &fval)
    { return getArrheniusidx(Ci,etai,thetai,fval) ; }

    void calculate(double T, std::pair<double,double> *fval) const ;
  } ;

  class reaction {
  public:
    struct sc {
      int species ;
      double nu ;
      double bot ;
      sc() {}
      sc(int s,double n) :species(s),nu(n),bot(0.0){ }
      sc(int s,double n, double b):species(s),nu(n),bot(b){ }
    } ;
    struct sc_int {
      sc_int() {}
      sc_int(int s, int n) { species=s; nu=n; }
      int species ;
      int nu ;
    } ;
    struct scdiff {
      scdiff() {}
      scdiff(int r,double n) { reaction=r; nu_diff = n; }
      int reaction ;
      double nu_diff ;
    } ;
    typedef std::vector<sc_int> stoic_vec_int ;
    typedef std::vector<sc> stoic_vec ;
    typedef std::vector<scdiff> stoic_diff_vec ;
    struct rates {
      double kf ;    // Forward reaction rate 
      double kcr ;   // reciprocal of the equilibrium constant
      double kfp ;   // dkfdT
      double kcp ;   // dkcdT
    } ;
    struct Rox {
      double KA1,KA2,KB1,KB2,KT1,KT2,KZ1,KZ2 ;
      double rho_s, D_s ;
      int Pindex ;
      int reaction_index ;
      Rox() {}
    } ;
    struct Prate {
      int reaction_index ;
      double power ;
      double Pref ;
      Prate(int r, double p, double ref): reaction_index(r),power(p),Pref(ref) {}
      Prate() { reaction_index = -1 ; }
    } ;
//ST------------------
    struct Condrate {
      int reaction_index ;
      double Econd ;
      std::vector<double> metals ;
      std::vector<double> oxides ;
      Condrate(int r, double E, std::vector<double> &M, std::vector<double> &O): reaction_index(r),Econd(E),metals(M),oxides(O)
      {
       std::cout << "Econd in constructor = " << Econd << std::endl ;
       for(size_t i=0; i<oxides.size(); i++)
         std::cout << "oxides in constructor = " << oxides[i] << std::endl ;
       for(size_t i=0; i<O.size(); i++)
         std::cout << "O in constructor = " << O[i] << std::endl ;
       for(size_t i=0; i<metals.size(); i++)
         std::cout << "metals in constructor = " << metals[i] << std::endl ;
       for(size_t i=0; i<M.size(); i++)
         std::cout << "M in constructor = " << M[i] << std::endl ;
      }
      Condrate() { reaction_index = -1 ; }
    } ;
//ST------------------



  private:
    int num_species ;
    int num_reactions ;
    stoic_vec *reactants ; // num_reactions size list of
    stoic_vec_int *reactants_int ; 
    stoic_vec *products ;  // stoichiometric coefficients
    stoic_vec_int *products_int ;

    // num_species size list of product-reactants stoichiometrics
    stoic_diff_vec *pmr ; // products - reactants

    int num_func_results ;

    // num_reactions size list of Kf and Kc function indexes
    std::vector<int> Kfi ; // list of indexes to Kf in func_res
    std::vector<int> Kci ; // list of indexes to Kc in func_res 

    // Thermodynamic Kc compiler/calculator
    thermoKc thermoKcModule ;  
    // curve fit function compiler/calculator
    curveFitFunctions ArrheniusModule ;
    // num_species size list of molecular masses
    vector<double > m ;
    vector<vector<double> > MBwts ; // species weights in the M-Body reactions.
    vector<bool> isMB ; // is M-Body reaction?
    vector<bool> hasMBwts ; // has M-Body weights
    vector<bool> isforwardonly ; // Does the reaction have only forward reaction
    vector<Prate > pressure_modifier ;
    vector<Rox > Rox_modifier ;
//ST--------------
    vector<Condrate > condensation_modifier ;
//ST--------------
  public:
    reaction() { reactants = 0 ; reactants_int = 0 ;
      products = 0 ; products_int = 0; pmr = 0 ;
      num_species = 0 ; num_reactions = 0 ; num_func_results=0 ;}
    void initialize(const EOS &eos, const reaction_db &reactions,
                    const species_db &species) ;

    int num_rates() const { return num_reactions ; }
    void extract_rates(rates *rate_info, const EOS &eos, const double T) const ;
    void compute_w(double *w, const rates *rate_info,
                   const double *mixture, const EOS::State &eos_state,
                   const EOS &eos) const ;
    void compute_w_implicit_explicit(double *a,double *w, const rates *rate_info,
                   const double *mixture, const EOS::State &eos_state,
                   const EOS &eos) const ;
    void compute_dwdrs(double *dwdr, int s, const rates *rate_info,
                       double dTdp, const double *dTdri,
                       const double *mixture, const EOS::State &eos_state,
                       const EOS &eos) const ;
    double compute_dwdp(int s, const rates *rate_info,
                        double dTdp, const double *mixture,
                        const EOS::State &eos_state,const EOS &eos) const ;
  } ;

  inline std::ostream &operator<<(std::ostream &s, const reaction &r)
  {
    return s ;
  }

  inline std::istream &operator>>(std::istream &s, const reaction &r)
  {
    return s ;
  }

  inline std::ostream &operator<<(std::ostream &s, const reaction::rates &r)
  {
    s << r.kf << " " << r.kcr << " " << r.kfp << " " << r.kcp << std::endl ;
    return s ;
  }

  inline std::istream &operator>>(std::istream &s, reaction::rates &r)
  {
    s >> r.kf >> r.kcr >> r.kfp >> r.kcp ;
    return s;
  }
}

namespace Loci {

  template<> struct data_schema_traits<fluidPhysics::reaction::rates> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(fluidPhysics::reaction::rates()) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::reaction::rates,kf) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::reaction::rates,kcr) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::reaction::rates,kfp) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::reaction::rates,kcp) ;
      return DatatypeP(ct) ;
    }
  } ;

  class reactionConverter {
    fluidPhysics::reaction &ref ;
  public:
    reactionConverter(fluidPhysics::reaction &iref) : ref(iref) {}
    int getSize() const {
      return 0 ;
    }
    void getState(char *buf, int &size) {
      size = getSize() ;
      std::cerr << "Reaction converter : getState not implemented" << std::endl ;
    }
    void setState(char *buf, int size) {
      std::cerr << "Reaction converter : setState not implemented" << std::endl ;
    }
  } ;
  
  template<> struct data_schema_traits<fluidPhysics::reaction> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef char Converter_Base_Type ;
    typedef reactionConverter Converter_Type ;
  } ;
}

#endif
