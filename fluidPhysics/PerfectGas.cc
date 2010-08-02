#include "PerfectGas.h"
#include "chemistry_db.h"
#include <vector>
#include <string>
#include <Tools/options_list.h>

using std::vector ;
using std::string ;
using std::cout ;
using std::cerr ;
using std::endl ;
using std::istream ;
using std::ostream ;
namespace fluidPhysics {
  
  class species_energy_func {
  public:
    virtual ~species_energy_func() {}
    virtual void e_cv(double T,double T2, double T3, double T4, double T5,
                      double TR, double TR2, double &e, double &cv) = 0 ;
    virtual void dcvdT(double T, double T2, double T3, double T4, double T5,
                       double TR, double TR2, double &dcvdT) = 0 ;
    virtual void set_energy(double eref, double Tref) = 0;
  } ;

  class species_omega_func {
  public:
    virtual ~species_omega_func() {}
    virtual void Omega(double T, double T2, double T3,double T4, double Tr,
                       double &O, double &Op) = 0 ;
    virtual void set_energy(double eref, double Tref, double Rs) = 0 ;
    virtual void set_entropy(double sref, double Tref,double Rs) = 0 ;
    virtual double entropy(double T, double Rs) = 0 ;
  } ;



  // define the function e(T) for each species in the mixture

  //struct reference contains the array of species reference information 
  struct reference {
    double Tref ;               // Tref: reference temperature
    double Pref ;               // Pref: reference pressure
    double sref ;               // sref: reference entropy 
    double href ;               // href: reference enthalpy
    double R ;                  // R: species gas constant
    reference(double Tref_,double Pref_, double sref_, double href_,double R_)
      : Tref(Tref_),Pref(Pref_),sref(sref_),href(href_),R(R_) {}
  } ;
  

  class species_vibrational_energy_func : public species_energy_func {
    vector<double> theta_v ;
    double niR,hfi,Ri ;
  public:
    void setup(double n, double R, vector<double> &thetas) {
      niR = n*R;
      hfi = 0 ;
      Ri = R ;
      theta_v = thetas ;
    }
    virtual void e_cv(double T, double T2, double T3, double T4, double T5,
                      double TR, double TR2, double &e, double &cv) ;
    virtual void dcvdT(double T, double T2, double T3, double T4, double T5,
                       double TR, double TR2, double &dcvdT) ;
    virtual void set_energy(double eref, double Tref) ;
  } ;


  void species_vibrational_energy_func::e_cv( double T, double T2, double T3, double T4, double T5, double TR, double TR2, double &e, double &cv) {
    cv = niR ;
    e = hfi + niR*T ;
    for(unsigned int i=0;i<theta_v.size();++i) {
      const double tvrT = min(theta_v[i]*TR,100.0) ;
      const double exptv = exp(tvrT) ;
      const double evib = Ri*theta_v[i]/(exptv - 1.) ;
      e += evib ;
      cv += evib*tvrT*exptv/(T*(exptv-1.));
    }
  }
  void species_vibrational_energy_func::dcvdT( double T, double T2, double T3, double T4, double T5, double TR, double TR2, double &dcvdT) {
    dcvdT = 0 ;
    
    for(unsigned int i=0;i<theta_v.size();++i) {
      const double Rv = Ri ;
      const double tv = theta_v[i] ;
      const double tv2 = tv*tv ;
      const double tvrT = min(tv*TR,100.0) ;
      const double tvrT2 = 0.5*tvrT ;
      const double shtv2 = sinh(tvrT2) ;
      const double chtv2 = cosh(tvrT2) ;
        
      dcvdT += -0.25*Rv*tv2*(2.*T*shtv2*shtv2 - tv*chtv2*shtv2)/
        (T4*shtv2*shtv2*shtv2*shtv2) ;
    }
  }

  void species_vibrational_energy_func::set_energy(double eref, double Tref) {
    hfi = 0 ;
    const double T = Tref ;
    const double T2 = T*T ;
    const double T3 = T2*T ;
    const double T4 = T3*T ;
    const double T5 = T4*T ;
    const double TR = 1./T ;
    const double TR2 = TR*TR ;
    double e,cv ;
    e_cv(T,T2,T3,T4,T5,TR,TR2,e,cv) ;
    hfi = eref - e ;
  }

  class species_curve_fit_energy_func : public species_energy_func {
    double cv0,cv1,cv2,cv3,cv4, cvr ;
    double e0, e1, e2, e3, e4, e5, er ;
  public:
    virtual void e_cv(double T,double T2,double T3,double T4,double T5, double TR,double TR2,
                      double &e, double &cv) {
      e =  e0+e1*T+e2*T2+e3*T3+e4*T4+e5*T5 + er*TR;
      cv = cv0+cv1*T+cv2*T2+cv3*T3+cv4*T4 + cvr*TR2;
    }
    virtual void dcvdT(double T,double T2,double T3,double T4,double T5, double TR,double TR2,
                       double &dcvdT) {
      dcvdT = cv1 + 2.*cv2*T + 3.*cv3*T2 + 4.*cv4*T3 - 2.*cvr*(TR2*TR) ;
    }
    double dels(double T, double Rs) {
      return (cv0+Rs)*log(T)+cv1*T+cv2*T*T/2.+cv3*T*T*T/3.+cv4*T*T*T*T
        -cvr/T/T/2. ;
    }
    species_curve_fit_energy_func() {
      cv0 = 0 ; cv1 = 0 ; cv2 = 0; cv3 = 0 ; cv4 = 0 ; cvr = 0 ;
      e0 = 0 ; e1 = 0; e2 = 0; e3 = 0; e4 = 0; e5 = 0 ;
      er = 0 ;
    }
    
    void establish_poly(double a1,double a2, double a3, double a4, double a5) {
      cv0 = a1 ; cv1 = a2 ; cv2 = a3; cv3 = a4 ; cv4 = a5 ; cvr = 0 ;
      e0 = 0 ; e1 = a1; e2 = a2/2.; e3 = a3/3.; e4 = a4/4.; e5 = a5/5. ;
      er = 0 ;
    }
    void establish_shomate(double a1, double a2, double a3, double a4, double ar) {
      cv0 = a1 ; cv1 = a2 ; cv2 = a3; cv3 = a4 ; cv4 = 0 ; cvr = ar ;
      e0 = 0 ; e1 = a1; e2 = a2/2.; e3 = a3/3.; e4 = a4/4.; e5 = 0. ;
      er = -ar ;
    }      
    virtual void set_energy(double eref, double Tref) {
      const double T = Tref ;
      const double T2 = T*T ;
      const double T3 = T2*T ;
      const double T4 = T3*T ;
      const double T5 = T4*T ;
      const double TR = 1./T ;
      e0 = eref -(e1*T+e2*T2+e3*T3+e4*T4+e5*T5 + er*TR) ;
    }
  } ;

  class species_vibrational_omega_func: public species_omega_func {
    double a1,a6,a7 ;
    double Ri ;
    vector<double> theta_v ;
  public:
    void setup(double n, double R, vector<double> &thetas) {
      a6 = 0 ;
      a7 = 0 ;
      a1 = 1+n ;
      Ri = R ;
      theta_v = thetas ;
    }
    virtual void Omega(double T, double T2, double T3, double T4, double Tr,
                       double &O, double &Op) {
      const double T2r = Tr*Tr ;
      const double onemlogT = 1.0 - log(T) ;
      // compute ideal gas contributions
      O = a1*onemlogT + a6*Tr - a7 ;
      Op = -Tr*a1 - a6*T2r ;
      // add in vibrational contributions
      for(unsigned int i=0;i<theta_v.size();++i) {
        const double tv = -theta_v[i] ;
        const double expT = exp(tv*Tr) ;
        O  += log(1.0 - expT) ;
        Op += tv*T2r*expT/(1.0 - expT) ;
      }
    }
    virtual void set_energy(double eref, double Tref, double Rs) {
      double e =  (a1-1.)*Ri*Tref ;
      double TR = 1./Tref ;
      for(unsigned int i=0;i<theta_v.size();++i) {
        const double tvrT = theta_v[i]*TR ;
        const double exptv = exp(min(tvrT,100.0)) ;
        const double evib = Ri*theta_v[i]/(exptv - 1.) ;
        e += evib ;
      }
      double hfi = eref -e ;
      a6 = hfi/Ri ;
    }
    virtual double entropy(double T, double Rs) {
      double s = Ri*(a7 + a1*log(T)) ;
      for(unsigned int i=0;i<theta_v.size();++i) {
        s += Ri*(theta_v[i]/(exp(theta_v[i]/T)-1.0)/T
                 - log(1.0-exp(-theta_v[i]/T))); 
      }
      return s ;
    }
    virtual void set_entropy(double sref, double Tref, double Rs) {
      a7 = 0 ;
      double sT = entropy(Tref,Rs) ;
      a7 = (sref-sT)/Ri ;
    }
  } ;

  class species_curve_fit_omega_func : public species_omega_func {
    double o1, o2, o3, o4, o5, o6, o7 ,er, Ri;
  public:
    virtual void Omega(double T,double T2, double T3, double T4, double Tr,
                       double &O, double &Op) {
      O = (1.+o1)*(1. - log(T)) + o2*T + o3*T2 + o4*T3 + o5*T4 +o6*Tr + o7 ;
      Op = o2 - (1.+o1)*Tr + o3*2.*T + o4*3.*T2 + o5*4.*T3 - o6*Tr*Tr ;
    }
    virtual void set_energy(double eref, double Tref, double Rs) {
      const double e1 = o1*Ri, e2 = -o2*Ri, e3 = -o3*2.*Ri ;
      const double e4 = -o4*3.*Ri, e5 = -o5*4.*Ri ;
      const double T = Tref ;
      const double T2 = T*T ;
      const double T3 = T2*T ;
      const double T4 = T3*T ;
      const double T5 = T4*T ;
      const double TR = 1./T ;
      const double e0 =  eref-(e1*T+e2*T2+e3*T3+e4*T4+e5*T5 + er*TR);
      o6 = e0/Ri ;
    }

    virtual double entropy(double T, double Rs) {
      const double cv0 = o1*Ri , cv1 = -o2*2.*Ri, cv2 = -o3*6.*Ri ;
      const double cv3 = -o4*12.*Ri, cv4 = -o5*20.*Ri, cvr = -er ;
      return (cv0+Ri)*log(T)+cv1*T+cv2*T*T/2.+cv3*T*T*T/3.+cv4*T*T*T*T/4.
        -cvr/(T*T*2.) - Ri*o7 ;
    }
    virtual void set_entropy(double sref, double Tref, double Rs) {
      o7 = 0 ;
      double sT = entropy(Tref,Ri) ;
      o7 = -(sref-sT)/Ri ;
    }
    void establish_poly(double a1,double a2, double a3, double a4, double a5, double Rs) {
      Ri = Rs ;
      o1 = a1/Ri ; o2 = -a2/2./Ri ; o3 = -a3/6./Ri; o4 = -a4/12./Ri;
      o5 = -a5/20./Ri ;
      er = 0 ;
      o6 = 0 ; o7 = 0 ; // ar has no effect on Omega
    }
    void establish_shomate(double a1,double a2, double a3, double a4, double ar, double Rs) {
      Ri = Rs ;
      o1 = a1/Ri ; o2 = -a2/2./Ri ; o3 = -a3/6./Ri; o4 = -a4/12./Ri ; o5 = 0. ;
      er  = -ar ;
      o6 = 0 ; o7 = 0 ; // ar has no effect on Omega
    }
  } ;

  class curve_fit_energy_func : public energy_function {
    int ns ;
    double Tmax ;
    vector<double> temp_ranges ;
    vector<vector<species_energy_func *> > efuncs ;
    vector<vector<species_omega_func *> > ofuncs ;
  public:
    void print_temp_ranges() const {
      cerr << "temp_ranges = " ;
      for(unsigned int k=0;k<temp_ranges.size();++k)
        cerr << temp_ranges[k] << " " ;
      cerr << endl ;
    }
    std::pair<int,int> get_temp_range(double Tl, double Th) ;
    void set_num_species(int numspecies) {
      ns = numspecies ;
      get_temp_range(0.,1e6) ;
      Tmax = 1e6 ;
    }

    int get_temp_interval(double T) const {
      for(int i=1;i<int(temp_ranges.size())-1;++i)
        if(T<temp_ranges[i])
          return i-1 ;
      return int(temp_ranges.size())-2 ;
    }
    void insert_energy_func(int species, double Tl, double Th,
                            species_energy_func *sef,
                            species_omega_func *sof) ;
    void set_energy(int i,double eref, double Tref, double Ri) ;
    void set_entropy(int i,double sref, double Tref, double Ri) ;
    virtual double compute_T(const double *rhoi, double energy,
                             double *ei, double *cvi, double Tguess) const ;
    virtual double compute_T_from_h(const double *yi,double enthalpy,
                                    double *ei,double *cvi,double Tguess,const std::vector<double> &R) const ;
    virtual void compute_ei_cvi(double T, double *ei, double *cvi) const ;

    // compute Omega and (del Omega / del T) for each species
    // The omega function is a partial product that results from a Gibbs free
    // energy minimization.  This is used to compute the equilibrium constant
    // Kc from thermodynamic Parameters
    virtual void compute_dcvdT(double T, double *dcvdT) const ;
    virtual void compute_Omega(double T, double *O, double *Op) const ;

  } ;

  void curve_fit_energy_func::set_energy(int i,double eref,double Tref, double Rs) {
    int t = get_temp_interval(Tref) ;

    warn(Tref < temp_ranges[t]) ;

    if(efuncs[t][i] == 0 || ofuncs[t][i] == 0) {
      cerr << "curve_fit unable to set energy at reference conditions"
           << endl ;
      abort() ;
    }
                            
    efuncs[t][i]->set_energy(eref,Tref) ;
    ofuncs[t][i]->set_energy(eref,Tref,Rs) ;
    for(int k=t+1;k<int(efuncs.size());++k) {
      if(efuncs[k-1][i] == efuncs[k][i])
        continue ;
      double e,cv ;
      const double T = temp_ranges[k]  ;
      const double T2 = T*T ;
      const double T3 = T2*T ;
      const double T4 = T3*T ;
      const double T5 = T4*T ;
      const double Tr = 1./T ;
      const double Tr2 = Tr*Tr ;
      efuncs[k-1][i]->e_cv(T,T2,T3,T4,T5,Tr,Tr2,e,cv) ;

      if(efuncs[k][i] == 0) {
        /* No curve, so initialize it to a linear extrapolation */
        Tmax = min(Tmax,T) ;
        species_curve_fit_energy_func *sef = new species_curve_fit_energy_func;
        species_curve_fit_omega_func *sof = new species_curve_fit_omega_func ;
        sef->establish_poly(cv,0.,0.,0.,0.) ;
        sof->establish_poly(cv,0.,0.,0.,0.,Rs) ;
        efuncs[k][i] = sef ;
        ofuncs[k][i] = sof ;
      }
        
      warn(efuncs[k][i] == 0) ;
      efuncs[k][i]->set_energy(e,T) ;
      ofuncs[k][i]->set_energy(e,T,Rs) ;
    }
    for(int k=t-1;k>=0;--k) {
      if(efuncs[k+1][i] ==efuncs[k][i])
        continue ;
      double e,cv ;
      const double T = temp_ranges[k+1]  ;
      const double T2 = T*T ;
      const double T3 = T2*T ;
      const double T4 = T3*T ;
      const double T5 = T4*T ;
      const double Tr = 1./T ;
      const double Tr2 = Tr*Tr ;

      efuncs[k+1][i]->e_cv(T,T2,T3,T4,T5,Tr,Tr2,e,cv) ;
      warn(efuncs[k][i] == 0) ;
      efuncs[k][i]->set_energy(e,T) ;
      ofuncs[k][i]->set_energy(e,T,Rs) ;
    }
  }

  void curve_fit_energy_func::set_entropy(int i,double sref,double Tref, double Rs) {
    int t = get_temp_interval(Tref) ;

    warn(Tref < temp_ranges[t]) ;
    
    ofuncs[t][i]->set_entropy(sref,Tref,Rs) ;
    for(int k=t+1;k<int(efuncs.size());++k) {
      if(ofuncs[k-1][i] == ofuncs[k][i])
        continue ;
      const double T = temp_ranges[k]  ;
      double s = ofuncs[k-1][i]->entropy(T,Rs) ;
      ofuncs[k][i]->set_entropy(s,T,Rs) ;
    }
    for(int k=t-1;k>=0;--k) {
      if(ofuncs[k+1][i] == ofuncs[k][i])
        continue ;
      const double T = temp_ranges[k+1]  ;
      double s = ofuncs[k+1][i]->entropy(T,Rs) ;
      ofuncs[k][i]->set_entropy(s,T,Rs) ;
    }
  }
  
  std::pair<int,int> curve_fit_energy_func::get_temp_range(double Tl, double Th) {
    warn(temp_ranges.size() == 1) ;
    if(temp_ranges.size() == 0) {
      temp_ranges.push_back(Tl) ;
      temp_ranges.push_back(Th) ;
      efuncs.push_back(vector<species_energy_func *>(ns)) ;
      ofuncs.push_back(vector<species_omega_func *>(ns)) ;
      for(int i=0;i<ns;++i) {
        efuncs[0][i] = 0 ;
        ofuncs[0][i] = 0 ;
      }
      return std::pair<int,int>(0,1) ;
    }
    if(Tl < temp_ranges.front()) {
      temp_ranges.push_back(temp_ranges.back()) ;
      efuncs.push_back(vector<species_energy_func *>(ns)) ;
      ofuncs.push_back(vector<species_omega_func *>(ns)) ;
      for(int i=temp_ranges.size()-1;i>0;--i)
        temp_ranges[i] = temp_ranges[i-1] ;
      for(int i=temp_ranges.size()-2;i>0;--i) {
        efuncs[i] = efuncs[i-1] ;
        ofuncs[i] = ofuncs[i-1] ;
      }
      
      temp_ranges[0] = Tl ;
      for(int i=0;i<ns;++i) {
        efuncs[0][i] = 0 ;
        ofuncs[0][i] = 0 ;
      }
    }

    if(Th > temp_ranges.back()) {
      temp_ranges.push_back(Th) ;
      efuncs.push_back(vector<species_energy_func *>(ns)) ;
      ofuncs.push_back(vector<species_omega_func *>(ns)) ;
      for(int i=0;i<ns;++i) {
        efuncs.back()[i] = 0 ;
        ofuncs.back()[i] = 0 ;
      }
    }
    int rl = -1, rh = -1 ;
    for(int i=1;i<int(temp_ranges.size());++i) {
      if(temp_ranges[i] > Tl) {
        rl = i-1 ;
        break ;
      }
    }
    for(int i=1;i<int(temp_ranges.size());++i) {
      if(temp_ranges[i] >= Th) {
        rh = i ;
        break ;
      }
    }
    if(rl == -1 || rh == -1) {
      cerr << "could not find interval" << endl ;
      abort() ;
    }
    if(temp_ranges[rh] != Th) { // split range at high end
      temp_ranges.push_back(temp_ranges.back()) ;
      for(int i=temp_ranges.size()-1;i>rh;--i)
        temp_ranges[i] = temp_ranges[i-1] ;
      efuncs.push_back(efuncs.back()) ;
      ofuncs.push_back(ofuncs.back()) ;
      for(int i=efuncs.size()-1;i>rh-1;--i) {
        efuncs[i] = efuncs[i-1] ;
        ofuncs[i] = ofuncs[i-1] ;
      }
      temp_ranges[rh] = Th ;
    }
    if(temp_ranges[rl] != Tl) { // split range at low end
      warn(temp_ranges[rl] > Tl) ;
      temp_ranges.push_back(temp_ranges.back()) ;
      for(int i=temp_ranges.size()-1;i>rl+1;--i)
        temp_ranges[i] = temp_ranges[i-1] ;
      efuncs.push_back(efuncs.back()) ;
      ofuncs.push_back(ofuncs.back()) ;
      for(int i=efuncs.size()-1;i>rl;--i) {
        efuncs[i] = efuncs[i-1] ;
        ofuncs[i] = ofuncs[i-1] ;
      }
      temp_ranges[rl+1] = Tl ;
      rl++ ;
      rh++ ;
      cerr << "split range at low end " << endl ;
      print_temp_ranges() ;
    }
    return std::pair<int,int>(rl,rh) ;
  }

  void curve_fit_energy_func::insert_energy_func(int species,
                                                 double Tl, double Th,
                                                 species_energy_func *sef,
                                                 species_omega_func *sof) {
    std::pair<int,int> er = get_temp_range(Tl,Th) ;
    for(int i=er.first;i<er.second;++i) {
      warn(efuncs[i][species] != 0) ;
      warn(ofuncs[i][species] != 0) ;

      efuncs[i][species] = sef ;
      ofuncs[i][species] = sof ;
    }
      
  }
  
  double curve_fit_energy_func::compute_T(const double *rhoi, double energy,
                                          double *ei, double *cvi,
                                          double Tguess) const {
    const int MAX_ITER = 100 ;
    double rho = 0.0 ;

    for(int i=0;i<ns;++i)
      rho+=rhoi[i];

    const double re = rho * energy ;
    // Use Newton Method to find T that corresponds to re
    // eps is error epsilon
    const double eps = 1e-8 ;
    
    double T = Tguess ;
    //    double T = 300.0 ;
    int t=get_temp_interval(T);
    double dT ;
    for(int iter=0;iter<MAX_ITER;++iter) {
      fatal(t<0) ;
      fatal(t>=static_cast<int>(temp_ranges.size())) ;
      fatal(t>=static_cast<int>(efuncs.size())) ;

      const double T2 = T*T ;
      const double T3 = T*T2 ;
      const double T4 = T*T3 ;
      const double T5 = T*T4 ;
      const double Tr = 1./T ;
      const double Tr2 = Tr*Tr ;
      double ref=0.0 ;
      double rcvf=0.0 ;
      //get species internal energy and specific heat
      for(int i=0;i<ns;++i) 
        efuncs[t][i]->e_cv(T,T2,T3,T4,T5,Tr,Tr2,ei[i],cvi[i]) ;
      
      for(int i=0;i<ns;++i) {
        ref += rhoi[i]*ei[i] ;
        rcvf += rhoi[i]*cvi[i] ;
      }
      dT = (re-ref)/rcvf ;
      T += dT ;
      if(abs(dT) < eps)
        break ;
      while(T>temp_ranges[t+1] && t < int(efuncs.size())-1)
        t++ ;
      while(T<temp_ranges[t] && t > 0)
        t-- ;
    }
    if(abs(dT) > eps) {
      // If things didn't converge, try again using bracketing algorithm
      double Tmin = 1 ;
      double Tmax = 15000 ;
      for(int iter=0;iter<MAX_ITER;++iter) {
        double ref=0.0 ;
        double rcvf=0.0 ;
        //get species internal energy and specific heat
        compute_ei_cvi(T,ei,cvi) ; 
        for(int i=0;i<ns;++i) {
          ref += rhoi[i]*ei[i] ;
          rcvf += rhoi[i]*cvi[i] ;
        }
        dT = (re-ref)/rcvf ;
        if(dT > 0)
          Tmin = max(T,Tmin) ;
        else
          Tmax = min(T,Tmax) ;
        
        T += dT ;
        if(T >= Tmax || T <= Tmin ) {
          T = (Tmax+Tmin)/2.0 ;
          dT = Tmax-Tmin ;
        }
        if(abs(dT) < eps)
          break ;
      }
      if(abs(dT) > eps) {
        cerr << "curve_fit_energy_func::getState failed to converge, dT = "
             << dT << ", T = " << T << endl ;
        abort() ;
      }
    }
    if(T<0) {
      cerr << "Warning, EOS converged to a non-physical negative temperature"
           << endl << "T = " << T << endl ;
      Loci::Abort() ;
    }
    if(T>Tmax) {
      cerr << "Warning, EOS converged to a value exceeding curve-fit specifications" << endl ;
      cerr << "Linear extrapolation utilized to compute EOS at T = " << T << endl ;
    }
    return T ;
  }
  
  double curve_fit_energy_func::compute_T_from_h(const double *yi, double enthalpy,
                                                 double *ei,double *cvi,double Tguess,const std::vector<double> &R) const {
    const int MAX_ITER=100 ; double rho = 0.0 ;
    for(int i=0;i<ns;++i) rho+=yi[i] ; const double rh = rho*enthalpy ;
    
    // Use Newton Method to find T that corresponds to rh; eps is error epsilon
    const double eps=1e-8 ; double T=Tguess ; int t=get_temp_interval(T); double dT ;
    const double Tlow=temp_ranges[1],Thigh=temp_ranges[temp_ranges.size()-2] ;
    for(int iter=0;iter<MAX_ITER;++iter) {
      fatal(t<0) ;
      fatal(t>=static_cast<int>(temp_ranges.size())) ; fatal(t>=static_cast<int>
                                                             (efuncs.size())) ;

      // Get species enthalpy and specific heat.
      const double T2=T*T,T3=T*T2,T4=T*T3,T5=T*T4,Tr=1.0/T,Tr2=Tr*Tr ;
      for(int i=0;i<ns;++i) efuncs[t][i]->e_cv(T,T2,T3,T4,T5,Tr,Tr2,ei[i],
                                               cvi[i]) ;
      
      double rhf=0.0,rcpf=0.0 ;
      for(int i=0;i<ns;++i) {
        rhf+=yi[i]*(ei[i]+R[i]*T) ; rcpf+=yi[i]*(cvi[i]+R[i]) ;
      }
      dT=(rh-rhf)/rcpf ; T+=dT ;
      if(T<Tlow) return 0.99*Tlow+0.01*Thigh ; if(T>Thigh) T=0.01*Tlow+0.99*Thigh ;
      if(abs(dT) < eps) break ;
      while(T>temp_ranges[t+1] && t < int(efuncs.size())-1) t++ ;
      while(T<temp_ranges[t] && t > 0) t-- ;
    }

    if(abs(dT) > eps) {
      // If things didn't converge, try again using bracketing algorithm
      double Tmin = 1 ; double Tmax = 15000 ;
      for(int iter=0;iter<MAX_ITER;++iter) {

        // Get species energy and specific heat.
        compute_ei_cvi(T,ei,cvi) ; 

        double rhf=0.0,rcpf=0.0 ;
        for(int i=0;i<ns;++i) {
          rhf+=yi[i]*(ei[i]+R[i]*T) ; rcpf+=yi[i]*(cvi[i]+R[i]) ;
        }
        dT=(rh-rhf)/rcpf ;
        if(dT > 0) Tmin = max(T,Tmin) ; else Tmax = min(T,Tmax) ;
        T += dT ;
        if(T >= Tmax || T <= Tmin ) { T=(Tmax+Tmin)/2.0 ; dT=Tmax-Tmin ; }
        if(abs(dT) < eps) break ;
      }
      if(abs(dT) > eps) {
        cerr << "curve_fit_energy_func::compute_T_from_h failed to converge, "
             << "dT = " << dT << ", T = " << T << endl ; abort() ;
      }
    }
    if(T<Tlow) {
      T=0.99*Tlow+0.01*Thigh ;
      //    cerr << "Warning, EOS converged to a non-physical negative temperature"
      //      << endl << "T = " << T << endl ; Loci::Abort() ;
    }
    if(T>Thigh) {
      T=0.01*Tlow+0.99*Thigh ;
      //    cerr << "Warning, EOS converged to a value exceeding curve-fit "
      //      << "specifications" << endl ;
      //    cerr << "Linear extrapolation utilized to compute EOS at T = "
      //      << T << endl ;
    }
    return T ;
  }

  void curve_fit_energy_func::compute_ei_cvi(double T, double *ei, double *cvi) const
  {
    int t=get_temp_interval(T);
    fatal(t>=static_cast<int>(efuncs.size())) ;
    const double T2 = T*T ;
    const double T3 = T2*T ;
    const double T4 = T3*T ;
    const double T5 = T4*T ;
    const double Tr = 1./T ;
    const double Tr2 = Tr*Tr ;
    for(int i=0;i<ns;++i) {
      fatal(efuncs[t][i] == 0) ;
      efuncs[t][i]->e_cv(T,T2,T3,T4,T5,Tr,Tr2,ei[i],cvi[i]) ;
    }
  }

  void curve_fit_energy_func::compute_dcvdT(double T, double *dcvdT) const
  {
    int t=get_temp_interval(T);
    fatal(t>=static_cast<int>(efuncs.size())) ;
    const double T2 = T*T ;
    const double T3 = T2*T ;
    const double T4 = T3*T ;
    const double T5 = T4*T ;
    const double Tr = 1./T ;
    const double Tr2 = Tr*Tr ;
    for(int i=0;i<ns;++i) {
      fatal(efuncs[t][i] == 0) ;
      efuncs[t][i]->dcvdT(T,T2,T3,T4,T5,Tr,Tr2,dcvdT[i]) ;
    }        
  }

  void curve_fit_energy_func::compute_Omega(double T, double *O, double *Op) const
  { 
    int t=get_temp_interval(T);
    fatal(t >= static_cast<int>(ofuncs.size())) ;
    const double T2 = T*T ;
    const double T3 = T2*T ;
    const double T4 = T3*T ;
    const double Tr = 1./T ;

    for(int i=0;i<ns;++i) {
      fatal(ofuncs[t][i] == 0) ;
      ofuncs[t][i]->Omega(T,T2,T3,T4,Tr,O[i],Op[i]) ;
    }
  }

  class vibrational_energy_func : public energy_function {
    // struct baseline contains the array of elements hf, niR whcih contribute 
    // to the linear part (translational and rotational energy) of species 
    // energy term
    struct baseline {
      double hf ;    // Heat of formation
      double niR ;   // Linear component of energy term
      baseline(double hfi,double niRi): hf(hfi), niR(niRi) {} 
    } ;
    // struct omega_baseline contains the array of elements which contribute 
    // to species Omega term calculated in funtion compute_Omega (Omega term is 
    // used when computing equilibrium constant and represents a partial
    // result of a Gibbs free energy minimization.
    struct omega_baseline {
      double a1 ; // a1 = 1+ni
      double a6 ; // absolute zero enthalapy correction term, e.g. a6 = h0/R
      double a7 ; // entropy constant correction term, e.g. a7 = s0/R
      omega_baseline(double ni,double h0,double s0,double R) {
        a1 = 1.0 + ni ;
        a6 = h0/R ;
        a7 = s0/R ;
      }
    } ;
    // struct vibrational contains the array of vibrational temperatures,
    // theta_v, and species gas constants, R,  which contribute to the
    // vibrational terms of species energy equation
    struct vibrational {
      double theta_v ; // Characteristic vibrational temperature
      double R ;
      int species ;  // species number
      vibrational(double tv,double Ri, int sn) :theta_v(tv),R(Ri),species(sn) {}
    } ;
    // four data members which contains the information described in the above 
    //  four structures
    std::vector<baseline> ideal_part ;
    std::vector<omega_baseline> omega_part ;
    std::vector<reference> ref ;
    std::vector<vibrational> vib_part ;
  public:
    // function add_ideal_part constructs ideal_part, omega_part without 
    // considering vibrational enery. It also contructs ref
    int add_ideal_part(double R, double ni, double hf,
                       double Tref, double Pref, double sref, double href) ;
    // function add_vibrational_part constructs vib_part and omega_part
    // Note: used when a species is diatomic, or polyatomic, and the
    // vibrational terms are considered 
    void add_vibrational_part(double theta_v, int species) ;
    void calculate_ei_cvi(double T, double *ei, double *cvi) const ;

    // compute internal enery and specific heat for each species
    virtual double compute_T(const double *rhoi, double energy,
                             double *ei, double *cvi, double Tguess) const ;
    virtual double compute_T_from_h(const double *yi, double enthalpy,
                                    double *ei, double *cvi, double Tguess,const std::vector<double> &R) const ;

    virtual void compute_ei_cvi(double T, double *ei, double *cvi) const ;
    // compute Omega and (del Omega / del T) for each species
    // The omega function is a partial product that results from a Gibbs free
    // energy minimization.  This is used to compute the equilibrium constant
    // Kc from thermodynamic Parameters
    virtual void compute_dcvdT(double T, double *dcvdT) const ;
    virtual void compute_Omega(double T, double *O, double *Op) const ;

  } ;

  int vibrational_energy_func::add_ideal_part(double R, double ni, double hf,
                                              double Tref, double Pref, double sref, double href) {
    // Note: this function adds a single species to the energy function, returning
    // the allocated address for this species.  ThermallyPerfectEOS::initialize()
    // is responsible for calling this function once for each species in the
    // mixture.  This adds the ideal gas part of the energy equation.  That is,
    // it adds the e(T) = hf + nRT component.  The vibrational components
    // are added using the vibrational_energy_func::add_vibrational_part() method.


    //put species hf, niR in ideal_part which is 
    //related to linear part of the energy
    ideal_part.push_back(baseline(hf,ni*R)) ; 
                                            
    // put species reference parameters in ref 
    ref.push_back(reference(Tref,Pref,sref,href,R)) ; 

    double s0 = sref - R*(1.0+ni)*log(Tref) ;
    double h0 = hf ;
    //put related information which contructs the Omega function in omega_part
    omega_part.push_back(omega_baseline(ni,h0,s0,R)) ; 
                                                 
    // if the number of elements in ideal_part and omega_part are not equal, 
    // give a warning
    warn(ideal_part.size() != omega_part.size()) ; 

    // return the index of last element in ideal_part
    return ideal_part.size() -1 ; 
  }

  void vibrational_energy_func::add_vibrational_part(double theta_v,
                                                     int species) {
    // Note: this function accomplishs the push_back action, resizing the
    // fectors vib_part and omega_part.  The whole vector will  be constructed
    // when calling the method ThermallyPerfectEOS::initialize().
    // The size of vib_part depends on the number of diatomic and polyatomic
    // species, and also on the number of atomic bonds in polyatomic species


    // the index species must have a corresponding component that makes the ideal
    // gas contribution from ideal_part, otherwise the specification is improper.
    // Remember that size of ideal_part is the current species number no matter
    // monatomic, diatomic or polyatomic species since they all contain
    // linear term in their energy equation.
    fatal(species >= static_cast<int>(ideal_part.size()) || species < 0); 

    const double R = ref[species].R ; 
    //put related information in vib_part
    vib_part.push_back(vibrational(theta_v,R,species)) ; 

    // add vibrational contribution to s0
    const double Tref = ref[species].Tref ;
    omega_part[species].a7 -= (theta_v/(exp(theta_v/Tref)-1.0)/Tref
                               - log(1.0-exp(-theta_v/Tref))); 
  }


  // Note: the above two methods are used to provide the species initial,
  // primitive information.  After doing this, we can compute species internal
  // energy, specific heat, Omega etc. defined in the following two methods. 
  // Method vibrational_energy_func::compute_ei_cvi() computes species
  // specific heat and energy, ( e(T), cv(T) )
  inline void vibrational_energy_func::calculate_ei_cvi(double T,
                                                        double *e,
                                                        double *cv) const {
    // Add ideal Gas (not considering vibrational energy) parts
    // (every species contributes this part)
    for(unsigned int i=0;i<ideal_part.size();++i) {
      const baseline &si = ideal_part[i] ;
      cv[i] = si.niR ;
      e[i] = si.hf + si.niR*T ;
    }

    // add in vibrational components ;
    const double Tr = 1.0/T ;
    for(unsigned int i=0;i<vib_part.size();++i) {
      const vibrational &vi = vib_part[i] ;
      const double tvrT = min(vi.theta_v*Tr,100.0) ;
      const double exptv = exp(tvrT) ;
      const double evib = vi.R*vi.theta_v/(exptv - 1.) ;
      e[vi.species] += evib ;
      cv[vi.species] += evib*tvrT*exptv/(T*(exptv-1.));
    }
  }

  // Note: the above two methods are used to provide the species initial,
  // primitive information.  After doing this, we can compute species internal
  // energy, specific heat, Omega etc. defined in the following two methods. 
  // Method vibrational_energy_func::compute_ei_cvi() computes species
  // specific heat and energy, ( e(T), cv(T) ). CHECK THE INTERNALS.

  inline void vibrational_energy_func::compute_dcvdT(double T, double *dcvdT) const {
    for(unsigned int i=0;i<ideal_part.size();++i) 
      dcvdT[i] = 0 ;
    if(vib_part.size() > 0) {
      double Tr = 1./T ;
      double T4 = T*T*T*T ;
      for(unsigned int i=0;i<vib_part.size();++i) {
        const vibrational &vi = vib_part[i] ;
        const double Rv = vi.R ;
        const double tv = vi.theta_v ;
        const double tv2 = tv*tv ;
        const double tvrT = min(tv*Tr,100.0) ;
        const double tvrT2 = 0.5*tvrT ;
        const double shtv2 = sinh(tvrT2) ;
        const double chtv2 = cosh(tvrT2) ;
        
        dcvdT[vi.species] += -0.25*Rv*tv2*(2.*T*shtv2*shtv2 - tv*chtv2*shtv2)/
          (T4*shtv2*shtv2*shtv2*shtv2) ;
      }
    }
  }
  
  inline void vibrational_energy_func::compute_ei_cvi(double T,
                                                      double *e,
                                                      double *cv) const {
    calculate_ei_cvi(T,e,cv) ;
  }
  
  double vibrational_energy_func::compute_T(const double *rhoi, double energy,
                                            double *ei, double *cvi, double Tguess) const {
    const int MAX_ITER = 100 ;
    double rho = 0.0 ;
    const int ns = ideal_part.size() ;
    for(int i=0;i<ns;++i)
      rho+=rhoi[i];

    const double re = rho * energy ;
    // Use Newton Method to find T that corresponds to re
    // eps is error epsilon
    const double eps = 1e-8 ;
    
    double T = Tguess ;
    //    double T = 300 ;
    double dT ;
    for(int iter=0;iter<MAX_ITER;++iter) {
      double ref=0.0 ;
      double rcvf=0.0 ;
      //get species internal energy and specific heat
      calculate_ei_cvi(T,ei,cvi) ; 
      for(int i=0;i<ns;++i) {
        ref += rhoi[i]*ei[i] ;
        rcvf += rhoi[i]*cvi[i] ;
      }
      dT = (re-ref)/rcvf ;
      T += dT ;
      if(abs(dT) < eps)
        break ;
    }
    if(abs(dT) > eps) {
      // If things didn't converge, try again using bracketing algorithm
      cerr << "problems converging EOS:  Trying bracketing algorithm" << endl ;
      cerr << "T = " << T <<", dT = " << dT << endl ;
      double Tmin = 1 ;
      double Tmax = 15000 ;
      T = 300 ;
      for(int iter=0;iter<MAX_ITER;++iter) {
        double ref=0.0 ;
        double rcvf=0.0 ;
        //get species internal energy and specific heat
        calculate_ei_cvi(T,ei,cvi) ; 
        for(int i=0;i<ns;++i) {
          ref += rhoi[i]*ei[i] ;
          rcvf += rhoi[i]*cvi[i] ;
        }
        dT = (re-ref)/rcvf ;

        if(dT > 0)
          Tmin = max(T,Tmin) ;
        else
          Tmax = min(T,Tmax) ;

        T += dT ; 
       
        if(T >= Tmax || T <= Tmin ) {
          T = (Tmax+Tmin)/2.0 ;
          dT = Tmax-Tmin ;
        }
        if(abs(dT) < eps)
          break ;
      }
      if(abs(dT) > eps) {
        cerr << "vibrational_energy_func::getState failed to converge, dT = "
             << dT << ", T = " << T << endl ;
        abort() ;
      }

    }
    if(T<0) {
      cerr << "Warning, EOS converged to a non-physical negative temperature"
           << endl << "T = " << T << endl ;
      Loci::Abort();
    }
    return T ;
  }


  double vibrational_energy_func::compute_T_from_h(const double *yi,double enthalpy,
                                                   double *ei,double *cvi,double Tguess,const std::vector<double> &R) const {
    const int MAX_ITER = 100 ; double rho = 0.0 ;
    const int ns = ideal_part.size() ;
    for(int i=0;i<ns;++i) rho+=yi[i];

    // Use Newton Method to find T that corresponds to re eps is error epsilon
    const double rh=rho*enthalpy ; const double eps=1e-8 ; double T=Tguess ; double dT ;
    for(int iter=0;iter<MAX_ITER;++iter) {
      double rhf=0.0 ; double rcpf=0.0 ;
      calculate_ei_cvi(T,ei,cvi) ; 
      for(int i=0;i<ns;++i) {
        rhf+=yi[i]*(ei[i]+R[i]*T) ; rcpf+=yi[i]*(cvi[i]+R[i]) ;
      }
      dT=(rh-rhf)/rcpf ; T+=dT ; if(abs(dT)<eps) break ;
      if(T<50.0) return 50.0 ;
      if(T>5000.0) return 4999.0 ;
    }
    if(abs(dT) > eps) {
      // If things didn't converge, try again using bracketing algorithm
      //    cerr << "problems converging EOS:  Trying bracketing algorithm" << endl ;
      //    cerr << "T = " << T <<", dT = " << dT << endl ;
      double Tmin = 1 ; double Tmax = 15000 ; T = 300 ;
      for(int iter=0;iter<MAX_ITER;++iter) {
        double rhf=0.0 ; double rcpf=0.0 ;
        calculate_ei_cvi(T,ei,cvi) ; 
        for(int i=0;i<ns;++i) {
          rhf+=yi[i]*(ei[i]+R[i]*T) ; rcpf+=yi[i]*(cvi[i]+R[i]) ;
        }
        dT=(rh-rhf)/rcpf ;
        if(dT > 0) Tmin = max(T,Tmin) ; else Tmax = min(T,Tmax) ;
        T += dT ; 
        if(T >= Tmax || T <= Tmin ) { T = (Tmax+Tmin)/2.0 ; dT = Tmax-Tmin ; }
        if(abs(dT) < eps) break ;
      }
      //    if(abs(dT) > eps) {
      //      cerr << "vibrational_energy_func::compute_T_from_h failed to "
      //        << "converge, dT = " << dT << ", T = " << T << endl ; abort() ;
      //    }
    }
    if(T<50.0) return 50.0 ;
    if(T>5000.0) return 4999.0 ;
    //  if(T<0) {
    //    cerr << "Warning, EOS converged to a non-physical negative temperature"
    //      << endl << "T = " << T << endl ; Loci::Abort();
    //  }
    return T ;
  }

  // Compute the partial product results of the Gibbs free energy minimization
  // used in the calculation of a thermodynamic equilibrium "Constant" Kc
  void vibrational_energy_func::compute_Omega(double T,
                                              double *Omega,
                                              double *Omegap) const {
    double Tr = 1.0/T ;
    double T2r = Tr*Tr ;
    double onemlogT = 1.0 - log(T) ;

    // compute ideal gas contributions
    for(unsigned int i=0;i<omega_part.size();++i) {
      const double a1 = omega_part[i].a1 ;
      const double a6 = omega_part[i].a6 ;
      const double a7 = omega_part[i].a7 ;
      Omega[i] = a1*onemlogT + a6*Tr - a7 ;
      Omegap[i] = -Tr*a1 - a6*T2r ;
    }

    // add in vibrational contributions
    for(unsigned int i=0;i<vib_part.size();++i) {
      const vibrational &vi = vib_part[i] ;
      const double tv = -vi.theta_v ;
      const double expT = exp(tv*Tr) ;
      Omega[vi.species] += log(1.0 - expT) ;
      Omegap[vi.species] += tv*T2r*expT/(1.0 - expT) ;
    }
  }

  /* Constructor of class ThermallyPerfectEOS */
  ThermallyPerfectEOS::ThermallyPerfectEOS(const ThermallyPerfectEOS &igm)
  {
    num_species = igm.num_species ;
    namelist = igm.namelist ;
    Pref = igm.Pref ;
    R = igm.R ;
    m  = igm.m ;
    mf = igm.mf ;

    num_elements = igm.num_elements ;
    species_elements = igm.species_elements ;
    efunc = igm.efunc ;
  }


  // Define operator= for class ThermallyPerfectEOS. 
  // Notes: we can write A=(B=C) if A, B, C are objects of class
  // ThermallyPerfectEOS by defining in the following way. If we
  // do not "return *this" in the end, we can write  A=B, but not A=(B=C) 

  ThermallyPerfectEOS &ThermallyPerfectEOS::operator=(const ThermallyPerfectEOS &igm)
  {
    num_species = igm.num_species ;
    namelist = igm.namelist ;
    Pref = igm.Pref ;
    R = igm.R ;
    m  = igm.m ;
    mf = igm.mf ;

    num_elements = igm.num_elements ;
    species_elements = igm.species_elements ;
    efunc = igm.efunc ;

    //"this" is the built-in pointer which points to the object of class
    return *this ;  
  }

  // Destructor of class ThermallyPerfectEOS 
  ThermallyPerfectEOS::~ThermallyPerfectEOS()
  {
  }

  EoSPtr ThermallyPerfectEOS::clone() const {
    return EoSPtr(new ThermallyPerfectEOS) ;
  }
  // Initialize species information 
  // Note: class species_db is defined in chemistry_db.h. This class provides
  // species data from the user provided chemistry model.
  void ThermallyPerfectEOS::initialize(const species_db &species,
                                       std::string thermodynamic_model) 
  {
    equationOfState::base_initialize(species) ;

    for(int i=0;i<species.numSpecies();++i) {
      const options_list &ol = species.getSpeciesOption(i) ;
      if(ol.optionExists("eos")) {
        cerr << "Thermally perfect EoS does not support generic EoS functions." << endl ;
        cerr << "Check model file or use different equation of state." << endl ;
        Loci::Abort() ;
      }
    }

    thermo_model = thermodynamic_model ;

    bool thermo_model_ve = false ;
    bool thermo_model_cf = false ;
    if(thermodynamic_model=="adaptive") {
      for(int i=0;i<species.numSpecies();++i) {
        const options_list &ol = species.getSpeciesOption(i) ;
        if(ol.optionExists("cp"))
          thermo_model_cf = true ;
      }
      if(thermo_model_cf == false)
        thermo_model_ve = true ;
    }
    if(thermodynamic_model == "vibrational_equilibrium") {
      thermo_model_ve = true ;
    }
    if(thermodynamic_model == "curve_fit") {
      thermo_model_cf = true ;
    }
    if(thermo_model_ve) {
      CPTR<vibrational_energy_func> efuncvib  = new vibrational_energy_func;
      for(int i=0;i<species.numSpecies();++i) {
        string sn = species.getSpeciesName(i) ;
        const options_list &ol = species.getSpeciesOption(i) ; 
      
        vector<double> theta_v ;
        if(ol.optionExists("theta_v")) {
          // If there are vibrational temperatures exist in the specification
          // then add them to the species energy function
          if(ol.getOptionValueType("theta_v") == Loci::REAL) {
            // For diatomic species (one vibrational temperature)
            double tv ;
            ol.getOption("theta_v",tv) ;
            // Note: here only push one theta_v of one species 
            // into vib_part (diatomic species has only one theta_v)
            theta_v.push_back(tv) ;
          } else if(ol.getOptionValueType("theta_v") == Loci::LIST) {
            // For polyatomic species (many vibrational temperatures)
            options_list::arg_list al ;
            ol.getOption("theta_v",al) ;
            double tv ;
            // Note: in the for loop, all theta_v of one species are pushed
            // into vib_part (polyatomic species has more than one theta_v)
            for(unsigned int j=0;j<al.size();++j) { 
              if(al[j].type_of() == Loci::REAL) {
                al[j].get_value(tv) ;
                theta_v.push_back(tv) ;
              } else {
                cerr << "theta_v has invalid type for species"
                     << sn << endl ;
                exit(1) ;
              }
            }
          } else {
            cerr << "theta_v has invalid type for species " << sn << endl ;
            exit(1) ;
          }
        }
  
        double nr,Trefr,srefr,hrefr ;
        Trefr = Tref[i] ;
        hrefr = href[i] ;
        srefr = sref[i] ;
        nr = n[i] ;

        double eref = hrefr - R[i]*Trefr ;
        double hf = eref-(nr*R[i]*Trefr) ;
        for(unsigned int k=0;k<theta_v.size();++k) {
          const double tvrT = theta_v[k]/Trefr ;
          const double evib = R[i]*theta_v[k]/(exp(min(tvrT,100.0))-1.) ;
          hf -= evib ;
        }
        // Add linear contributions to energy equation
        efuncvib->add_ideal_part(R[i],nr,hf,Trefr,Pref[i],
                                 srefr,hrefr) ;

        for(unsigned int t=0;t<theta_v.size();++t) 
          efuncvib->add_vibrational_part(theta_v[t],i) ;
      }
      efunc = CPTR<energy_function>(efuncvib) ;
    } else if(thermo_model_cf == true) {
      CPTR<curve_fit_energy_func> ecurve  = new curve_fit_energy_func;
      const int ns = num_species ;
      ecurve->set_num_species(ns) ;
      for(int i=0;i<species.numSpecies();++i) {
        string sn = species.getSpeciesName(i) ;
        const options_list &ol = species.getSpeciesOption(i) ; 
      
        vector<double> theta_v ;
        if(ol.optionExists("theta_v")) {
          // If there are vibrational temperatures exist in the specification
          // then add them to the species energy function
          if(ol.getOptionValueType("theta_v") == Loci::REAL) {
            // For diatomic species (one vibrational temperature)
            double tv ;
            ol.getOption("theta_v",tv) ;
            // Note: here only push one theta_v of one species 
            // into vib_part (diatomic species has only one theta_v)
            theta_v.push_back(tv) ;
          } else if(ol.getOptionValueType("theta_v") == Loci::LIST) {
            // For polyatomic species (many vibrational temperatures)
            options_list::arg_list al ;
            ol.getOption("theta_v",al) ;
            double tv ;
            // Note: in the for loop, all theta_v of one species are pushed
            // into vib_part (polyatomic species has more than one theta_v)
            for(unsigned int j=0;j<al.size();++j) { 
              if(al[j].type_of() == Loci::REAL) {
                al[j].get_value(tv) ;
                theta_v.push_back(tv) ;
              } else {
                cerr << "theta_v has invalid type for species"
                     << sn << endl ;
                exit(1) ;
              }
            }
          } else {
            cerr << "theta_v has invalid type for species " << sn << endl ;
            exit(1) ;
          }
        }
        double nr,Trefr,srefr,hrefr ;
        Trefr = Tref[i] ;
        hrefr = href[i] ;
        srefr = sref[i] ;
        nr = n[i] ;
        if(!ol.optionExists("cp") || ol.getOptionValueType("cp") != Loci::LIST) {
          cerr << "species " << sn << " does not have a curve fit cp defined"
               << endl ;
          abort() ;
        } else {
          options_list::arg_list al ;
          ol.getOption("cp",al) ;
          double Th=20000, Tl=0 ;
          if(al[0].type_of() == Loci::REAL)
            al[0].get_value(Th) ;
          else {
            cerr << "species " << sn << " error in cp specification" << endl ;
          }
          species_energy_func *sef = 0;
          species_omega_func *sof = 0;

          species_vibrational_energy_func *svf =
            new species_vibrational_energy_func ;
          svf->setup(nr,R[i],theta_v) ;
          sef = svf ;
          species_vibrational_omega_func *svo =
            new species_vibrational_omega_func ;
          svo->setup(nr,R[i],theta_v) ;
          sof = svo ;
          ecurve->insert_energy_func(i,0,Th,sef,sof) ;
        
          // go through list
          warn((al.size() & 1) != 1) ;
          for(unsigned int j=1;j<al.size();++j) {
            if(al[j].type_of() == Loci::FUNCTION) {
              std::string fname ;
              al[j].get_value(fname) ;
              options_list::arg_list fl ;
              al[j].get_value(fl) ;
              if(fname == "shomate") {
                if(fl.size() != 5) {
                  cerr << "shomate curve fit for cp has 5 arguments" << endl ;
                  cerr << "problem occured when parsing species " << sn << endl ;
                  break ;
                }
                for(int k=0;k<5;++k)
                  if(fl[k].type_of() != Loci::REAL) {
                    cerr << "shomate syntax error for species " << sn << endl ;
                  }
                double Ar, Br, Cr, Dr, Er ;
                fl[0].get_value(Ar) ;
                fl[1].get_value(Br) ;
                fl[2].get_value(Cr) ;
                fl[3].get_value(Dr) ;
                fl[4].get_value(Er) ;

                const double coef = 1000.0/m[i] ;
                const double A = Ar*coef -R[i];
                const double B = Br*coef*1.e-3 ;
                const double C = Cr*coef*1.e-6 ;
                const double D = Dr*coef*1.e-9 ;
                const double E = Er*coef*1.e6 ;
                species_curve_fit_energy_func *ef =
                  new species_curve_fit_energy_func ;
                ef->establish_shomate(A,B,C,D,E) ;
                sef = ef ;
                species_curve_fit_omega_func *of =
                  new species_curve_fit_omega_func ;
                of->establish_shomate(A,B,C,D,E,R[i]) ;
                sof = of ;
              } else if(fname == "poly") {
                if(fl.size() != 5) {
                  cerr << "polynomial curve fit for cp has 5 arguments"
                       << endl ;
                  cerr << "problem occured when parsing species " << sn
                       << endl ;
                  break ;
                }
                for(int k=0;k<5;++k)
                  if(fl[k].type_of() != Loci::REAL) {
                    cerr << "poly cp syntax error for species " << sn << endl ;
                  }
                double Ar, Br, Cr, Dr, Er ;
                fl[0].get_value(Ar) ;
                fl[1].get_value(Br) ;
                fl[2].get_value(Cr) ;
                fl[3].get_value(Dr) ;
                fl[4].get_value(Er) ;

                const double coef = 1000.0/m[i] ;
                const double A = Ar*coef -R[i];
                const double B = Br*coef ;
                const double C = Cr*coef ;
                const double D = Dr*coef ;
                const double E = Er*coef ;
                species_curve_fit_energy_func *ef =
                  new species_curve_fit_energy_func ;
                ef->establish_poly(A,B,C,D,E) ;
                sef = ef ;
                species_curve_fit_omega_func *of =
                  new species_curve_fit_omega_func ;
                of->establish_poly(A,B,C,D,E,R[i]) ;
                sof = of ;
              } else {
                cerr << "species " << sn << " has an unknown curve-fit specification function " << fname << endl ;
                break ;
              }
            } else {
              cerr << "species " << sn << "error in cp speicification" << endl ;
              break ;
            }
            ++j ;
            Tl = Th ;
            if(al[j].type_of() == Loci::REAL)
              al[j].get_value(Th) ;
            else {
              cerr << "species " << sn << " error in cp specification" << endl ;
              break ;
            }
            ecurve->insert_energy_func(i,Tl,Th,sef,sof) ;
          }
        }
        double eref = hrefr - R[i]*Trefr ;
        ecurve->set_energy(i,eref,Trefr,R[i]) ;
        ecurve->set_entropy(i,srefr,Trefr,R[i]) ;

      }
      efunc = CPTR<energy_function>(ecurve) ;
    } else {
      cerr << "unknown thermodynamic model `" << thermodynamic_model
           << "', "<< endl
           << "currently valid options are `vibrational_equilibrium'"
           << " and `curve_fit'" << endl ;
      abort() ;
    }
    
  
  }

  // Note: for a non-equilibrium mixture, any combinations of any two state
  // variables and mixture specification can specify a thermodynamic state.
  // If it is assumed that the species densities provide both the density and
  // the mixture specification, then one additional state variable is required
  // to specify the thermodynamic state.  The following three methods get
  // this thermodynamic state for a second state variable of T, P, and e,
  // respectively.

  ThermallyPerfectEOS::State
  ThermallyPerfectEOS::State_from_rho_T(const double *rhoi, double T, double  *ms,
                                        const float *hint) const
  {

    const int ns = num_species ;
    //These two variables are in local scope (they are not the data members rRt,
    // rho of class State, although they have same names.)
    double rRt=0.0,rho=0.0 ; 
                        
    for(int i=0;i<ns;++i) {
      rho += rhoi[i] ;
      rRt += rhoi[i]*R[i] ;
    }
    
    //this re is local variable, same thing as above
    double re  = 0.0 ; 
    double rcvf=0.0 ;
    double *ei = &ms[0] ;
    double *cvi = &ms[ns] ;
    efunc->compute_ei_cvi(T,ei,cvi) ;
    for(int i=0;i<ns;++i) {
      re += rhoi[i]*ei[i] ;
      rcvf += rhoi[i]*cvi[i] ;
    }

    // Copy these local state variables to data memeber of class State and return
    return State(T,rho,rcvf,rRt,re) ; 
  }

  // method State_from_rho_e gets the state from mixture density and
  // internal energy. This involves a newton iteration procedure to
  // solve the non-linear equation (rho*e - Sum_i{rho_i*e_i(T)}) = 0

  ThermallyPerfectEOS::State
  ThermallyPerfectEOS::State_from_rho_e(const double *rhoi, double energy,
                                        double *ms,
                                        const float *hint) const
  {

    const int ns = num_species ;
    //double T = Tguess ;
    double T = 300 ;
    if(hint != 0)
      T = hint[0] ;
    
    double *ei = &ms[0] ;
    double *cvi = &ms[ns] ;
    T = efunc->compute_T(rhoi,energy,ei,cvi,T) ;

    double rcvf = 0, rRt = 0,rho=0 ;

    for(int i=0;i<ns;++i) {
      rho += rhoi[i] ;
      rcvf += rhoi[i]*cvi[i] ;
      rRt += rhoi[i]*R[i] ;
    }
    return State(T,rho,rcvf,rRt,rho*energy) ;
  }

  // method State_from_rho_h gets the state from mixture density and
  // static enthalpy. This involves a newton iteration procedure to
  // solve the non-linear equation (rho*h - Sum_i{rho_i*h_i(T)}) = 0

  ThermallyPerfectEOS::State
  ThermallyPerfectEOS::State_from_rho_h(const double *rhoi, double enthalpy,
                                        double *ms,
                                        const float *hint) const
  {

    const int ns = num_species ;
    //double T = Tguess ;
    double T = 300 ;
    if(hint != 0)
      T = hint[0] ;

    // Note the temporary creation of yi so we can use the existing function
    // compute_T_from_h which takes mass fractions.
    double *ei = &ms[0] ;
    double *cvi = &ms[ns] ;
    double rho=0.0 ; for(int i=0;i<ns;++i) rho+=rhoi[i] ;
    double *yi=new double[ns] ; for(int i=0;i<ns;++i) yi[i]=rhoi[i]/rho ;
    T = efunc->compute_T_from_h(yi,enthalpy,ei,cvi,T,R) ;
    delete [] yi ;

    double rcvf = 0, rRt = 0 ;

    for(int i=0;i<ns;++i) {
      rcvf += rhoi[i]*cvi[i] ; rRt += rhoi[i]*R[i] ;
    }
    return State(T,rho,rcvf,rRt,rho*(enthalpy-rRt*T)) ;
  }

  // method State_from_p_h gets the state from mixture pressure and enthalpy.
  // This involves a newton iteration procedure to solve the non-linear equation
  // (rho*h - Sum_i{rho_i*h_i(T)}) = 0 .
  ThermallyPerfectEOS::State ThermallyPerfectEOS::State_from_p_h(const double
                                                                 *yi,double p,double enthalpy,double *mState,const float *hint) const {
    const int ns = num_species ;
    double T=300 ; if(hint != 0) T=hint[0] ;
    double *ei=&mState[0],*cvi=&mState[ns] ;

    T=efunc->compute_T_from_h(yi,enthalpy,ei,cvi,T,R) ;

    double rcvf=0.0,rRt=0.0 ;
    for(int i=0;i<ns;++i) { rcvf+=yi[i]*cvi[i] ; rRt+=yi[i]*R[i] ; }

    double rho=p/(T*rRt) ;
    return State(T,rho,rho*rcvf,rho*rRt,rho*(enthalpy-rRt*T)) ;
  }

  ThermallyPerfectEOS::State ThermallyPerfectEOS::State_from_p_e(const double
                                                                 *yi,double p,double energy,double *mState,const float *hint) const {
    const int ns = num_species ;
    double T=300 ; if(hint != 0) T=hint[0] ;
    double *ei=&mState[0],*cvi=&mState[ns] ;

    T=efunc->compute_T(yi,energy,ei,cvi,T) ;

    double rcvf=0.0,rRt=0.0 ;
    for(int i=0;i<ns;++i) { rcvf+=yi[i]*cvi[i] ; rRt+=yi[i]*R[i] ; }

    double rho=p/(T*rRt) ;
    return State(T,rho,rho*rcvf,rho*rRt,rho*energy) ;
  }

#ifdef NOTUSED
  double ThermallyPerfectEOS::getEntropy(const State &s, const double *mixture) const{
    return 0 ;
  }
#endif

  ostream & operator<<(ostream &s, const ThermallyPerfectEOS &eos)
  {
    s.precision(16) ;
    s << eos.thermo_model << endl ;
    s << eos.sdb << endl ;
    
    return s ;
  }

  istream & operator>>(istream &s, ThermallyPerfectEOS &eos)
  {
    string therm_model ;
    Loci::parse::kill_white_space(s) ;
    s >> therm_model ;
    species_db sdb ;
    s >> sdb ;
    eos.initialize(sdb,therm_model) ;
    return s ;
  }

}
