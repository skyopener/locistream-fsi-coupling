// This file comes directly from CHEM.
#include "fluidState.h"
#include "reaction.h"

#include "const.h"

#include <vector>
using std::vector ;

using std::map ;
using std::string ;
using std::swap ;

#include "Vector.h"

namespace streamUns {

  EOS::State fluidStateGetState(const real *q,
				const conservativeVectorInfo &qvi,
				const EOS &eos) {
    const int mi = qvi.momentumIndex() ;
    const int ei = qvi.totalEnergyIndex() ;
    const int ns = qvi.numSpecies() ;
    const int lsgki = qvi.lessgkIndex() ;
    const int nsgk = qvi.numLESsgk() ;

    const real &ru = q[mi+0] ;
    const real &rv = q[mi+1] ;
    const real &rw = q[mi+2] ;
    const real &re0 = q[ei] ;
    
    const real *r = q ;
    real rho = 0.0 ;
    for(int i=0;i<ns;++i)
      rho+=r[i];
    real rr = 1.0/rho ;

    real re = re0 - 0.5*(ru*ru+rv*rv+rw*rw)*rr ;
    //for(int i=0;i<nte;++i)
    //re -= q[tei+i] ;
    for(int i=0;i<nsgk;++i)
      re -= q[lsgki+i] ;
    return eos.State_from_rho_e(r,re*rr) ;
  }

  void setTemperature(const EOS &eos, real T, real *q,
                      const conservativeVectorInfo &qvi) {
    const int mi = qvi.momentumIndex() ;
    const real &ru = q[mi+0] ;
    const real &rv = q[mi+1] ;
    const real &rw = q[mi+2] ;
    const int ei = qvi.totalEnergyIndex() ;

    EOS::State s = eos.State_from_rho_T(q,T) ;
    const real rho = s.density() ;
    const real re = s.rho_energy() ;
    q[ei] = re + 0.5*(ru*ru+rv*rv+rw*rw)/rho ;
  }

  void setPressure(const EOS &eos, real p, real *q,
                   const conservativeVectorInfo &qvi) {
    
    const int mi = qvi.momentumIndex() ;
    const real &ru = q[mi+0] ;
    const real &rv = q[mi+1] ;
    const real &rw = q[mi+2] ;
    const int ei = qvi.totalEnergyIndex() ;

    EOS::State s = eos.State_from_rho_p(q,p) ;
    const real rho = s.density() ;
    const real re = s.rho_energy() ;
    q[ei] = re + 0.5*(ru*ru+rv*rv+rw*rw)/rho ;
  }

  void setDensity(const EOS &eos, real rho, real *q,
                  const conservativeVectorInfo &qvi) {
    const vector<real> &mf = eos.getMixtureFractions() ;
    real *r = q ;
    for(int i=0;i<eos.numSpecies();++i)
      r[i] = rho*mf[i] ;
  }

  namespace BLACK_BOX {
    class equilibrium_function {
    protected:
      int ne ;  // number of elements
      int ns ;  // number of species
      int mi ;  // momentum index
      int ei ;  // energy index
      int vs ;  // vector size 
      vector<int> el_species ;  // list,first ne species represent elements
      vector<vector<int> > els ;     // element stoichiometry matrix
      const EOS *eos ;
      const conservativeVectorInfo *qvi ;
      const reaction *r ;
    public:
      virtual void function(real *result, const real *vars) = 0 ;
      virtual void set_ic(real *icval) = 0 ;
      virtual ~equilibrium_function() {}
    } ;

    class rho_P_function : public equilibrium_function {
    protected:
      real T ;
      real P ;
      real rho ;
      real *qtmp,*wtmp ;
      real *ic ;
      real *mcoeff ;
    public:
      rho_P_function(const real *q, const conservativeVectorInfo &s_,
                     const EOS &e, const reaction &rr,
                     const vector<int> & es, const vector<vector<int> > &el,
                     bool useP=true) {
        eos = &e ;
        qvi = &s_ ;
        ne = eos->numElements() ;
        ns = eos->numSpecies() ;
        mi = qvi->momentumIndex() ;
        ei = qvi->totalEnergyIndex() ;
        vs = qvi->vectorSize() ;
        r = &rr ;
        el_species = es ;
        els = el ;
        EOS::State s = fluidStateGetState(q,*qvi,*eos) ;
        P = s.pressure() ;
        T = s.temperature() ;
        rho = s.density() ;
        qtmp = new real[vs] ;
        wtmp = new real[vs] ;
        mcoeff = new real[ne] ;
        ic = new real[vs] ;
        
        for(int e=0;e<ne;++e)
          mcoeff[e] = 0.0 ;
        for(int e=0;e<ne;++e)
          for(int s=0;s<ns;++s)
            mcoeff[e] += els[s][e]*q[s]/(rho*eos->speciesMolecularMass(s)) ;

        // Find reasonable initial guess
        for(int i=0;i<vs;++i)
          qtmp[i] = q[i] ;
        const int MAX_ITER = 20 ;
        for(int it2=0;it2<2;++it2) {
          real cf = 1.0 ;
          int nr = r->num_rates() ;
          tmp_array<reaction::rates> rate_info(nr) ;
          for(int iter=0;iter<MAX_ITER;++iter) {
            // Now begin estimating mass fractions
            if(useP)
              setPressure(*eos,P,qtmp,*qvi) ;
            else
              setTemperature(*eos,T,qtmp,*qvi) ;
            s = fluidStateGetState(qtmp,*qvi,*eos) ;
            rho = s.density() ;
            r->extract_rates(rate_info,*eos,s.temperature()) ;

            tmp_array<real> mixture(ns) ;
	    for(int i=0;i<ns;++i)
	      mixture[i] = qtmp[i]/rho ;
            r->compute_w(wtmp,rate_info,mixture,s,*eos) ;
            for(int s=0;s<ns;++s) {
              if(wtmp[el_species[s]] > 0.0) {
                qtmp[el_species[s]] +=
                  cf*(1.0-qtmp[el_species[s]]/rho)*rho ;
              } else {
                qtmp[el_species[s]] -= cf*qtmp[el_species[s]] ;
              }
            }
            real rho_new = 0.0 ;
            for(int i=0;i<ns;++i)
              rho_new += qtmp[i] ;
            for(int i=0;i<ns;++i)
              qtmp[i] = qtmp[i]*rho/rho_new ;
            cf = cf/2.0 ;
          }
            
        }

        s = fluidStateGetState(qtmp,*qvi,*eos) ;

        for(int i=0;i<ns;++i)
          ic[i] = qtmp[i] ;
        if(useP)
          ic[ns] = s.temperature() ;
        else
          ic[ns] = T ;
        
      }
    
      virtual ~rho_P_function() ;
      virtual void function(real *result, const real *vars) ;
      virtual void set_ic(real *icval) ;
    } ;

    rho_P_function::~rho_P_function()
    {
      delete[] qtmp ;
      delete[] wtmp ;
      delete[] mcoeff ;
      delete[] ic ;
    }

    void rho_P_function::function(real *result, const real *vars)
    {
      for(int s=0;s<ns;++s)
        qtmp[s] = vars[s] ;
      for(int i=0;i<3;++i)
        qtmp[mi+i] = 0.0 ;
      const real T = vars[ns] ;
      setTemperature(*eos,T,qtmp,*qvi) ;
      EOS::State s = fluidStateGetState(qtmp,*qvi,*eos) ;
      int nr = r->num_rates() ;
      tmp_array<reaction::rates> rate_info(nr) ;

      r->extract_rates(rate_info,*eos,s.temperature()) ;

      tmp_array<real> mixture(ns) ;

      for(int i=0;i<ns;++i)
	mixture[i] = qtmp[i]/s.density() ; ;
      r->compute_w(wtmp,rate_info,mixture,s,*eos) ;

      for(int e=0;e<ne;++e) 
        result[e] = 0.0 ;
      for(int e=0;e<ne;++e)
        for(int s=0;s<ns;++s)
          result[e] += els[s][e]*qtmp[s]/(rho*eos->speciesMolecularMass(s)) ;
      for(int e=0;e<ne;++e)
        result[e] -= mcoeff[e] ;

      for(int s=ne;s<ns;++s)
        result[s] = wtmp[el_species[s]] ;
      const real Ptmp = s.pressure() ;
      result[ns] = Ptmp - P ;
    }

    void rho_P_function::set_ic(real *ic_vec)
    {
      for(int i=0;i<ns+1;++i)
        ic_vec[i] = ic[i] ;
    }

    class rho_T_function : public rho_P_function {
    public:
      rho_T_function(const real *q, const conservativeVectorInfo &qvi,
                     const EOS &e, const reaction &rr,
                     const vector<int> & es, const vector<vector<int> > &el) :
        rho_P_function(q,qvi,e,rr,es,el,false) {}
      virtual void function(real *result, const real *vars) ;
    } ;

    void rho_T_function::function(real *result, const real *vars)
    {
      for(int s=0;s<ns;++s)
        qtmp[s] = vars[s] ;
      for(int i=0;i<3;++i)
        qtmp[mi+i] = 0.0 ;
      const real Tval = vars[ns] ;
      setTemperature(*eos,Tval,qtmp,*qvi) ;
      EOS::State s = fluidStateGetState(qtmp,*qvi,*eos) ;
      int nr = r->num_rates() ;
      tmp_array<reaction::rates> rate_info(nr) ;
      r->extract_rates(rate_info,*eos,s.temperature()) ;

      tmp_array<real> mixture(ns) ;
      for(int i=0;i<ns;++i)
	mixture[i] = qtmp[i]/s.density() ;
      r->compute_w(wtmp,rate_info,mixture,s,*eos) ;
      
      for(int e=0;e<ne;++e) 
        result[e] = 0.0 ;
      for(int e=0;e<ne;++e)
        for(int s=0;s<ns;++s)
          result[e] += els[s][e]*qtmp[s]/(rho*eos->speciesMolecularMass(s)) ;
      for(int e=0;e<ne;++e)
        result[e] -= mcoeff[e] ;

      for(int s=ne;s<ns;++s)
        result[s] = wtmp[el_species[s]] ;
      result[ns] = Tval - T ;
    }

    class P_T_function : public equilibrium_function {
    protected:
      real T ;
      real P ;
      real rho ;
      real *qtmp,*wtmp ;
      real *ic ;
      real *mcoeff ;
    public:
      P_T_function(const real *q, const conservativeVectorInfo &s_,
                   const EOS &e, const reaction &rr,
                   const vector<int> & es, const vector<vector<int> > &el,
                   real Pref,real Tref) {
        eos = &e ;
        qvi = &s_ ;
        ne = eos->numElements() ;
        ns = eos->numSpecies() ;
        mi = qvi->momentumIndex() ;
        ei = qvi->totalEnergyIndex() ;
        vs = qvi->vectorSize() ;
        r = &rr ;
        el_species = es ;
        els = el ;
        EOS::State s = fluidStateGetState(q,*qvi,*eos);
        P = Pref ;
        T = Tref ;
        rho = s.density() ;
        qtmp = new real[vs] ;
        wtmp = new real[vs] ;
        mcoeff = new real[ne] ;
        ic = new real[vs] ;
        
        for(int e=0;e<ne;++e)
          mcoeff[e] = 0.0 ;
        for(int e=0;e<ne;++e)
          for(int s=0;s<ns;++s)
            mcoeff[e] += els[s][e]*q[s]/(rho*eos->speciesMolecularMass(s)) ;

        // Find reasonable initial guess
        for(int i=0;i<vs;++i)
          qtmp[i] = q[i] ;
        const int MAX_ITER = 20 ;
        real cf = 1.0 ;
        int nr = r->num_rates() ;
        tmp_array<reaction::rates> rate_info(nr) ;

        for(int iter=0;iter<MAX_ITER;++iter) {
          // First estimate density
          setTemperature(*eos,Tref,qtmp,*qvi) ;
          s = fluidStateGetState(qtmp,*qvi,*eos) ;
          rho = s.density() ;
          const real Rt = s.gasConstant() ;
          real new_rho = Pref/Rt/Tref ;
          for(int i=0;i<ns;++i)
            qtmp[i] = qtmp[i]/rho*new_rho ;
          // Now begin estimating mass fractions
          setTemperature(*eos,Tref,qtmp,*qvi) ;
          s = fluidStateGetState(qtmp,*qvi,*eos) ;
          rho = s.density() ;
          r->extract_rates(rate_info,*eos,s.temperature()) ;

          tmp_array<real> mixture(ns) ;

	  for(int i=0;i<ns;++i)
	    mixture[i] = qtmp[i]/s.density() ;
          r->compute_w(wtmp,rate_info,mixture,s,*eos) ;

          for(int s=0;s<ns;++s) {
            if(wtmp[el_species[s]] > 0.0) {
              qtmp[el_species[s]] +=
                cf*(1.0-qtmp[el_species[s]]/rho)*rho ;
            } else {
              qtmp[el_species[s]] -= cf*qtmp[el_species[s]] ;
            }
          }
          real rho_new = 0.0 ;
          for(int i=0;i<ns;++i)
            rho_new += qtmp[i] ;
          for(int i=0;i<ns;++i)
            qtmp[i] = qtmp[i]*rho/rho_new ;
          cf = cf/2.0 ;
        }

            
        for(int i=0;i<ns;++i)
          ic[i] = qtmp[i] ;
        ic[ns] = rho ;
      }
      virtual ~P_T_function() ;
      virtual void function(real *result, const real *vars) ;
      virtual void set_ic(real *icval) ;
    } ;

    P_T_function::~P_T_function()
    {
      delete[] qtmp ;
      delete[] wtmp ;
      delete[] mcoeff ;
      delete[] ic ;
    }

    void P_T_function::function(real *result, const real *vars)
    {
      const real r2 = vars[ns] ;
    
      for(int s=0;s<ns;++s) 
        qtmp[s] = vars[s] ;
      for(int i=0;i<3;++i)
        qtmp[mi+i] = 0.0 ;
      setTemperature(*eos,T,qtmp,*qvi) ;
      EOS::State s = fluidStateGetState(qtmp,*qvi,*eos) ;
      int nr = r->num_rates() ;
      tmp_array<reaction::rates> rate_info(nr) ;

      r->extract_rates(rate_info,*eos,s.temperature()) ;
      tmp_array<real> mixture(ns) ;

      for(int i=0;i<ns;++i)
	mixture[i] = qtmp[i]/s.density() ;
      r->compute_w(wtmp,rate_info,mixture,s,*eos) ;

      for(int e=0;e<ne;++e) {
        result[e] = 0.0 ;
        for(int s=0;s<ns;++s)
          result[e] += els[s][e]*qtmp[s]/(r2*eos->speciesMolecularMass(s)) ;
        result[e] -= mcoeff[e] ;
      }

      for(int s=ne;s<ns;++s)
        result[s] = wtmp[el_species[s]] ;
      const real Ptmp = s.pressure() ;
      result[ns] = Ptmp - P ;
    }

    void P_T_function::set_ic(real *ic_vec)
    {
      for(int i=0;i<ns+1;++i)
        ic_vec[i] = ic[i] ;
    }

    class find_equilibrium_state {
      int ns ;
      const EOS *eos ;
      equilibrium_function *f ;
      Vector<int> bi ;
      Vector<real> q,x,y,b,h ;
      Matrix<real> A ;
      void compute_jacobian() ;
      void solve_for_x() ;
    public:
      find_equilibrium_state(const EOS &e, equilibrium_function &ef) ;
      void get_equilibrium_densities(real *rho) ;
    } ;

    void find_equilibrium_state::get_equilibrium_densities(real *rho)
    {
      for(int i=0;i<ns;++i)
        rho[i] = q[i] ;
    }
    
    void find_equilibrium_state::compute_jacobian()
    {
      for(int i=0;i<ns+1;++i) {
        real tmp = q[i] ;
        q[i] += h[i] ;
        f->function(x,q) ;
        q[i] = tmp ;
        for(int j=0;j<ns+1;++j) 
          A[i][j] = (x[j]-b[j])/h[i] ;
      }
    }

    void find_equilibrium_state::solve_for_x() {
      const int vs = ns+1 ;

      for(int i=0;i<vs;++i)
        bi[i] = i ;

      // Begin LU factorization
      for(int i=0;i<vs-1;i++) {
        int pivot = i ;
        // find pivot
        for(int j=i;j<vs;++j)
          if(abs(A[i][bi[j]])>abs(A[i][bi[pivot]]))
            pivot = j ;
        if(i != pivot)
          swap(bi[i],bi[pivot]) ;

        // compute l_ij's
        for(int j=i+1;j<vs;++j)
          A[i][bi[j]] = A[i][bi[j]]/A[i][bi[i]] ;
        // compute u_ij's
        for(int j=i+1;j<vs;++j)
          for(int k=i+1;k<vs;++k)
            A[k][bi[j]] -= A[i][bi[j]]*A[k][bi[i]] ;
      }

      // do forward solve Ly = b
      for(int i=0;i<vs;++i) {
        y[i] = b[bi[i]] ;
        for(int j=0;j<i;++j)
          y[i] -= A[j][bi[i]]*y[j] ;
      }

      // do backward solve Ux = y
      for(int i=vs-1;i>=0;--i) {
        x[i] = y[i] ;
        for(int j=i+1;j<vs;++j)
          x[i] -= A[j][bi[i]]*x[j] ;
        x[i] = x[i]/A[i][bi[i]] ;
      }
    }


    find_equilibrium_state::find_equilibrium_state(const EOS &e,
                                                   equilibrium_function &ef) {
      eos = &e ;
      ns = eos->numSpecies() ;
      f = &ef ;
      q.setVecSize(ns+1) ;
      x.setVecSize(ns+1) ;
      y.setVecSize(ns+1) ;
      b.setVecSize(ns+1) ;
      bi.setVecSize(ns+1) ;
      h.setVecSize(ns+1) ;
      A.setMatSize(ns+1) ;
      f->set_ic(q) ;

      const real eps = 1e-7 ;

      real rho =0.0;
      for(int i=0;i<ns;++i)   
        rho += q[i] ;
      for(int i=0;i<ns;++i)
        h[i] = rho*eps ;
      h[ns] = q[ns]*eps ;
      // Begin newton interation
      int iter ;
      const int MAX_ITER = 1000 ;
      for(iter=0;iter<MAX_ITER;++iter) {
        f->function(b,q) ;
        compute_jacobian() ;
        for(int i=0;i<ns+1;++i)
          b[i] = -b[i] ;
        solve_for_x() ;
        int flg=true ;
        for(int j=0;j<ns+1;++j) {
          if(h[j] < abs(x[j]))
            flg = false ;
          q[j] = abs(q[j]+x[j]) ;
        }
        if(flg)
          break ;
      }
      if(iter == MAX_ITER) {
        cerr << "newton method failed to converge in find_equilibrium_state"
             << endl ;
      }
    }
        
  }

  using namespace BLACK_BOX ;

  int fluidState::pack_buf_size() {
    int buf_size = 27 ;
    buf_size += 2*mixture.size() ;
    std::map<std::string,real>::const_iterator mi;
    for(mi=mixture.begin();mi!=mixture.end();++mi)
      buf_size+= mi->first.size() ;
    return buf_size ;
  }
  void fluidState::pack_buf(real *buf, int sz) {
    int i = 0 ;
    buf[i++] = p ;
    buf[i++] = T ;
    buf[i++] = rho ;
    buf[i++] = u.x ;
    buf[i++] = u.y ;
    buf[i++] = u.z ;
    buf[i++] = k ;
    buf[i++] = e ;
    buf[i++] = w ;
    buf[i++] = R ;
    buf[i++] = nu_t ;
    buf[i++] = nu_tm ;
    buf[i++] = tmu ;
    buf[i++] = real(use_mach_velocity) ;
    buf[i++] = real(set_density_pressure) ;
    buf[i++] = real(set_density_temperature) ;
    buf[i++] = real(set_pressure_temperature) ;
    buf[i++] = real(set_k) ;
    buf[i++] = real(set_e) ;
    buf[i++] = real(set_w) ;
    buf[i++] = real(set_R) ;
    buf[i++] = real(set_nu_t) ;
    buf[i++] = real(set_nu_tm) ;
    buf[i++] = real(set_tmu) ;
    buf[i++] = real(equilibrium_mixture) ;
    buf[i++] = real(values_defined) ;
    buf[i++] = real(mixture.size()) ;
    std::map<std::string,real>::const_iterator mi;
    for(mi=mixture.begin();mi!=mixture.end();++mi) {
      string name = mi->first ;
      int ssz = name.size() ;
      buf[i++] = real( ssz ) ;
      for(int j=0;j<ssz;++j)
        buf[i++] = real( name[j] ) ;
      buf[i++] = mi->second ;
    }
    warn(i++!=sz) ;
  }
  void fluidState::unpack_buf(real *buf, int sz) {
    int i = 0 ;
    p  = buf[i++] ;
    T  = buf[i++] ;
    rho  = buf[i++] ;
    u.x  = buf[i++] ;
    u.y  = buf[i++] ;
    u.z  = buf[i++] ;
    k  = buf[i++] ;
    e  = buf[i++] ;
    w  = buf[i++] ;
    R  = buf[i++] ;
    nu_t  = buf[i++] ;
    nu_tm  = buf[i++] ;
    tmu  = buf[i++] ;
    use_mach_velocity = bool(buf[i++]) ;
    set_density_pressure = bool(buf[i++]) ;
    set_density_temperature = bool(buf[i++]) ;
    set_pressure_temperature = bool(buf[i++]) ;
    set_k = bool(buf[i++]) ;
    set_e = bool(buf[i++]) ;
    set_w = bool(buf[i++]) ;
    set_R = bool(buf[i++]) ;
    set_nu_t = bool(buf[i++]) ;
    set_nu_tm = bool(buf[i++]) ;
    set_tmu = bool(buf[i++]) ;
    equilibrium_mixture = bool(buf[i++]) ;
    values_defined = bool(buf[i++]) ;

    int msize = int(buf[i++]) ;

    mixture.clear() ;
    for(int j=0;j<msize;++j) {
      int ssz = int(buf[i++]) ;
      string name ;
      for(int k=0;k<ssz;++k) {
        char c = char(buf[i++]) ;
        name += c ;
      }
      mixture[name] = buf[i++] ;
    }
    warn(i++!=sz) ;
  }
  
  fluidState::fluidState()
  {
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
    use_mach_velocity = false ;
    set_density_pressure = false ;
    set_density_temperature = false ;
    set_pressure_temperature = false ;
    equilibrium_mixture = false ;
    set_k = false ;
    set_e = false ;
    set_w = false ;
    set_R = false ;
    set_nu_t = false ;
    set_nu_tm = false ;
    set_tmu = false ;
    values_defined = false ;
  }

  namespace {
    // returns a vector of species where the first ne elements represent the
    // "elemental" species.
    void find_elemental_species(vector<int> &el_species,
                                const vector<vector<int> > &els,
                                int ns, int ne) {
      // the elemental species are a set of species who have linearly
      // independent stoichiometry vectors. This routine finds these
      // vectors by performing an lu factorization with partial pivoting.
      // The partial pivoting will push the linearly dependent rows to
      // the bottom of the matrix, thus the top ne elements will be
      // linearly independent.

      // Set up stoichiometry matrix
      vector<vector<real> > scratch ;
      for(int s=0;s<ns;++s) {
        vector<real> t ;
        for(int e=0;e<ne;++e)
          t.push_back(real(els[s][e])) ;
        scratch.push_back(t) ;
      }

      // initialize pivot vector
      for(int s=0;s<ns;++s)
        el_species.push_back(s) ;

      // perform a partial lu factorization with pivoting.  We are not
      // interested in the lu factorization, just in the pivot vector.
      // (el_species is the pivot vector)
      for(int i=0;i<ne;++i) {
        int pivot = i ;
        // find pivot
        for(int j=i;j<ns;++j)
          if(abs(scratch[el_species[j]][i])
             >abs(scratch[el_species[pivot]][i]))
            pivot = j ;
        swap(el_species[i],el_species[pivot]) ;
            
        // compute l_ij
        for(int j=i+1;j<ns;++j)
          scratch[el_species[j]][i] /= scratch[el_species[i]][i] ;
        // compute u_ij
        for(int j=i+1;j<ns;++j)
          for(int k=i+1;k<ne;++k)
            scratch[el_species[j]][k] -= scratch[el_species[j]][i]*
              scratch[el_species[i]][k] ;
      }
    }
  }
  void fluidState::setEquilibriumState(real *q, const conservativeVectorInfo &qvi,
                                       const EOS &eos, const reaction &r) const
  {
    vector<vector<int> > els ;
    const int ne = eos.numElements() ;
    const int ns = eos.numSpecies() ;
    for(int s=0;s<ns;++s) {
      vector<int> v = eos.getSpeciesElements(s) ;
      els.push_back(v) ;
    }

    // find species that represent elemental composition
    vector<int> el_species ;
    find_elemental_species(el_species,els,ns,ne) ;

    if(set_density_pressure) {
      rho_P_function frp(q,qvi,eos,r,el_species,els) ;
      find_equilibrium_state e_state(eos,frp) ;
      e_state.get_equilibrium_densities(q) ;
      setPressure(eos,p,q,qvi) ;
    } else if(set_density_temperature) {
      rho_T_function frp(q,qvi,eos,r,el_species,els) ;
      find_equilibrium_state e_state(eos,frp) ;
      e_state.get_equilibrium_densities(q) ;
      setTemperature(eos,T,q,qvi) ;
    } else if(set_pressure_temperature) {
      P_T_function frp(q,qvi,eos,r,el_species,els,p,T) ;
      find_equilibrium_state e_state(eos,frp) ;
      e_state.get_equilibrium_densities(q) ;
      setTemperature(eos,T,q,qvi) ;
    } else
      cerr << "don't know how to set equilibrium state" << endl ;
    
  }

  bool fluidState::get_k_omega(real &kout, real &wout) const {
    const real beta = 1e-2 ;
    real v = norm(u) ;
    if(use_mach_velocity)
      v*= 340 ;
    kout = (3./2.)*(beta*v)*(beta*v) ;
    kout = max(kout,1e-3) ;
    if(set_k)
      kout = k ;
    const real nu_frac = 1e-3 ;
    const real nu = 1e-5 ;
    wout = 0.09*kout/(nu_frac*nu) ;
    if(set_w)
      wout = w ;
    if((!set_k && set_w) || (set_k && !set_w)) {
      cerr << "Both k and w must be set in fluidstate!" << endl ;
    }
    //    if(Loci::MPI_rank == 0) {
      //      cout << "fluidstate: setting k in domain = " << kout << endl ;
      //      cout << "fluidstate: setting w in domain = " << wout << endl ;
    //    }
    return (set_k && set_w) ;
  }

  bool fluidState::get_k_epsilon(real &kout, real &eout) const {
    const real beta = 1e-2 ;
    real v = norm(u) ;
    if(use_mach_velocity)
      v*= 340 ;
    kout = (3./2.)*(beta*v)*(beta*v) ;
    kout = max(kout,1e-3) ;
    if(set_k)
      kout = k ;
    const real nu_frac = 1e-3 ;
    const real nu = 1e-5 ;
    eout = 0.005*0.09*kout*kout/(nu_frac*nu) ;
//    const real len_scale = 1e-3 ;
//    eout = 0.09*pow(kout,1.5)/len_scale ;
    if(set_e)
      eout = e ;
    if((!set_k && set_e) || (set_k && !set_e)) {
      cerr << "Both k and e must be set in fluidstate!" << endl ;
    }
    //    if(Loci::MPI_rank == 0) {
    //      cout << "fluidstate: setting k in domain = " << kout << endl ;
    //      cout << "fluidstate: setting e in domain = " << eout << endl ;
    //    }
    return (set_k && set_e) ;
  }

  bool fluidState::get_R(real &Rout) const {
    Rout = 5e-6 ;
    if(set_R)
      Rout = R ;
    //    if(Loci::MPI_rank == 0) {
    //      cout << "fluidstate: setting R in domain = " << Rout << endl ;
    //    }
    return (set_R) ;
  }

  bool fluidState::get_nu_t(real &nu_tout) const {
    nu_tout = 5e-9 ;
    if(set_nu_t)
      nu_tout = nu_t ;
    //    if(Loci::MPI_rank == 0) {
    //      cout << "fluidstate: setting nu_t in domain = " << nu_tout << endl ;
    //    }
    return (set_nu_t) ;
  }

  bool fluidState::get_nu_tm(real &nu_tmout) const {
    nu_tmout = 5e-9 ;
    if(set_nu_tm)
      nu_tmout = nu_tm ;
    //    if(Loci::MPI_rank == 0) {
    //      cout << "fluidstate: setting nu_tm in domain = " << nu_tmout << endl ;
    //    }
    return (set_nu_tm) ;
  }

  void fluidState::setState(real *q, const conservativeVectorInfo &qvi,
                            const EOS &eos, const reaction &r) const
  {
    if(!values_defined)
      return ;
    const int ns = qvi.numSpecies() ;
    const int mi = qvi.momentumIndex() ;
    const int ei = qvi.totalEnergyIndex() ;
    const int vs = qvi.vectorSize() ;

    for(int i=0;i<vs;++i)
      q[i] = 0.0 ;

    EOS::State s ;

    std::vector<real> mf = eos.getMixtureFractions() ;

    warn(int(mf.size()) != ns) ;
    
    if(mixture.size()) {
      for(int i=0;i<ns;++i)
        mf[i] = 0.0 ;
      map<string,real>::const_iterator mi ;
      for(mi=mixture.begin();mi!=mixture.end();++mi) {
        int s = qvi.speciesIndex(mi->first) ;
        if(s == -1)
          cerr << "species " << mi->first
               << " not in EOS in fluidState::setState" << endl ;
        else
          mf[s] = mi->second ;
      }
      real mft = 0.0 ;
      for(int i=0;i<ns;++i)
        mft += mf[i] ;
      
      for(int i=0;i<ns;++i)
        mf[i] = mf[i]/mft ;
    }
    
    if(set_pressure_temperature) 
      s = eos.State_from_mixture_p_T(&mf[0],p,T) ;
    if(set_density_pressure) {
      for(int i=0;i<ns;++i)
        q[i] = mf[i]*rho ;
      s = eos.State_from_rho_p(q,p) ;
    }
    if(set_density_temperature) {
      for(int i=0;i<ns;++i)
        q[i] = mf[i]*rho ;
      s = eos.State_from_rho_T(q,T) ;
    }

    real srho = s.density() ;

    for(int i=0;i<ns;++i)
      q[i] = srho*mf[i] ;
    q[mi+0] = 0 ;
    q[mi+1] = 0 ;
    q[mi+2] = 0 ;
    q[ei] = s.rho_energy() ;

    if(equilibrium_mixture && ns > 1) {
      setEquilibriumState(q,qvi,eos,r) ;
      srho = 0 ;
      for(int i=0;i<ns;++i)
        srho += q[i] ;
      s = eos.State_from_rho_e(q,q[ei]/srho) ;
    }

    vect3d v = u ;
    if(use_mach_velocity)
      v *= s.soundSpeed() ;

    q[mi+0] = srho*v.x ;
    q[mi+1] = srho*v.y ;
    q[mi+2] = srho*v.z ;
    q[ei] += 0.5*srho*dot(v,v) ;
    
  }
  void fluidState::setPrimitive(real *qp, real Pambient,
                                const conservativeVectorInfo &qvi,
                                const EOS &eos, const reaction &r) const
  {
    if(!values_defined)
      return ;
    const int ns = qvi.numSpecies() ;
    const int mi = qvi.momentumIndex() ;
    const int ei = qvi.totalEnergyIndex() ;
    const int vs = qvi.vectorSize() ;

    for(int i=0;i<vs;++i)
      qp[i] = 0.0 ;

    EOS::State s ;

    std::vector<real> mf = eos.getMixtureFractions() ;

    warn(int(mf.size()) != ns) ;
    
    if(mixture.size()) {
      for(int i=0;i<ns;++i)
        mf[i] = 0.0 ;
      map<string,real>::const_iterator mi ;
      for(mi=mixture.begin();mi!=mixture.end();++mi) {
        int s = qvi.speciesIndex(mi->first) ;
        if(s == -1)
          cerr << "species " << mi->first
               << " not in EOS in fluidState::setState" << endl ;
        else
          mf[s] = mi->second ;
      }
      real mft = 0.0 ;
      for(int i=0;i<ns;++i)
        mft += mf[i] ;
      
      for(int i=0;i<ns;++i)
        mf[i] = mf[i]/mft ;
    }
    
    if(set_pressure_temperature) 
      s = eos.State_from_mixture_p_T(&mf[0],p,T) ;
    if(set_density_pressure) {
      for(int i=0;i<ns;++i)
        qp[i] = mf[i]*rho ;
      s = eos.State_from_rho_p(qp,p) ;
    }
    if(set_density_temperature) {
      for(int i=0;i<ns;++i)
        qp[i] = mf[i]*rho ;
      s = eos.State_from_rho_T(qp,T) ;
    }

    real srho = s.density() ;

    for(int i=0;i<ns;++i)
      qp[i] = srho*mf[i] ;
    qp[mi+0] = 0 ;
    qp[mi+1] = 0 ;
    qp[mi+2] = 0 ;
    qp[ei] = s.pressure() -Pambient ; // Pambient here

    if(equilibrium_mixture && ns > 1) {
      tmp_array<real> q(qvi.vectorSize()) ;
      for(int i=0;i<vs;++i)
        q[i] = qp[i] ;
      q[ei] = s.rho_energy() ;
      setEquilibriumState(q,qvi,eos,r) ;
      for(int i=0;i<ns;++i)
        qp[i] = q[i] ;
    }

    vect3d v = u ;
    if(use_mach_velocity)
      v *= s.soundSpeed() ;

    qp[mi+0] = v.x ;
    qp[mi+1] = v.y ;
    qp[mi+2] = v.z ;
    
  }

  void fluidState::Input(const options_list &finput)
  {

    rho = 1.0 ;
    u = vect3d(0,0,0) ;
    bool rho_set = false ;
    bool p_set = false ;
    bool T_set = false ;
    bool M_set = false ;
    bool u_set = false ;

    if(finput.optionExists("rho")) {
      if(finput.getOptionValueType("rho") == Loci::REAL) {
        finput.getOption("rho",rho) ;
        rho_set = true ;
      } else if(finput.getOptionValueType("rho") == Loci::UNIT_VALUE) {
        Loci::UNIT_type rhou ;
        finput.getOption("rho",rhou) ;
        if(!rhou.is_compatible("kg/m/m/m")) 
          cerr << "wrong type of unit for density: " << rhou << endl ;
        else {
          rho = rhou.get_value_in("kg/m/m/m") ;
          rho_set = true ;
        }
      } else {
        cerr << "incorrect type for 'rho' in fluidState" << endl ;
      }
    }
    if(finput.optionExists("T")) {
      if(finput.getOptionValueType("T") == Loci::REAL) {
        finput.getOption("T",T) ;
        T_set = true ;
      } else if(finput.getOptionValueType("T") == Loci::UNIT_VALUE) {
        Loci::UNIT_type Tu ;
        finput.getOption("T",Tu) ;
        if(!Tu.is_compatible("kelvin")) 
          cerr << "wrong type of unit for temperature: " << Tu << endl ;
        else {
          T = Tu.get_value_in("kelvin") ;
          T_set = true ;
        }
      } else {
        cerr << "incorrect type for 'T' in fluidState" << endl ;
      }
    }
    
    if(finput.optionExists("p")) {
      if(finput.getOptionValueType("p") == Loci::REAL) {
        finput.getOption("p",p) ;
        p_set = true ;
      } else if(finput.getOptionValueType("p") == Loci::UNIT_VALUE) {
        Loci::UNIT_type pu ;
        finput.getOption("p",pu) ;
        if(!pu.is_compatible("Pa")) 
          cerr << "wrong type of unit for pressure: " << pu << endl ;
        else {
          p = pu.get_value_in("Pa") ;
          p_set = true ;
        }
      } else {
        cerr << "incorrect type for 'p' in fluidState" << endl ;
      }
    }

    if(finput.optionExists("k")) {
      finput.getOption("k",k) ;
      set_k = true ;
    }

    if(finput.optionExists("w")) {
      finput.getOption("w",w) ;
      set_w = true ;
    }

    if(finput.optionExists("e")) {
      finput.getOption("e",e) ;
      set_e = true ;
    }

    if(finput.optionExists("R")) {
      finput.getOption("R",R) ;
      set_R = true ;
    }

    if(finput.optionExists("nu_t")) {
      finput.getOption("nu_t",nu_t) ;
      set_nu_t = true ;
    }

    if(finput.optionExists("nu_tm")) {
      finput.getOption("nu_tm",nu_tm) ;
      set_nu_tm = true ;
    }

    if(finput.optionExists("tmu")) {
      finput.getOption("tmu",tmu) ;
      set_tmu = true ;
    }
    
    if(finput.optionExists("M")) {
      u = get_vect3d(finput,"M","") ;
      M_set = true ;
    }
    
    if(finput.optionExists("u")) {
      u = get_vect3d(finput,"u","m/s") ;
      u_set = true ;
    }
    
    equilibrium_mixture = false ;
    if(finput.optionExists("equilibrium")) {
      if(finput.getOptionValueType("equilibrium") != Loci::BOOLEAN) 
        cerr << "equilibrium should be a boolean" << endl ;
      equilibrium_mixture = true ;
    }
    if(finput.optionExists("mixture")) {
      if(finput.getOptionValueType("mixture") != Loci::LIST) {
        cerr << "mixture should be assigned to a species list" << endl ;
      } else {
        Loci::options_list::arg_list species_list ;
        finput.getOption("mixture",species_list) ;
        Loci::options_list::arg_list::iterator li ;
        for(li=species_list.begin();li!=species_list.end();++li) {
          Loci::option_values::value_list_type species_arg ;
          li->get_value(species_arg) ;
          if(li->type_of() != Loci::NAME_ASSIGN || species_arg.size() != 1
             || species_arg.front().type_of() != Loci::REAL) 
            cerr << "error in mixture assignment" << endl ;
          else {
            string sn ;
            double sv ;
            li->get_value(sn) ;
            species_arg.front().get_value(sv) ;
            mixture[sn] = sv ;
          }
        }
      }
    }

    if(rho_set && p_set)
      set_density_pressure = true ;
    if(rho_set && T_set)
      set_density_temperature = true ;
    if(p_set && T_set)
      set_pressure_temperature = true ;
    if(rho_set && p_set && T_set)
      cerr << "overspecification of rho,p, and T in reading fluidState"
           << endl ;
    if(M_set)
      use_mach_velocity = true ;
    else if(u_set) {
      use_mach_velocity = false ;
    } else {
      use_mach_velocity = false ;
      u = vect3d(0,0,0) ;
    }
    if(M_set && u_set)
      cerr << "overspecification of M and u in reading fluidState"
           << endl ;
    values_defined = true ;
  
  }

  ostream &fluidState::Print(ostream &s) const 
  {
    if(values_defined) {
      s << "< " ;
      if(set_density_pressure || set_density_temperature) {
        s << "rho = " << rho ;
        if(set_density_pressure) {
          s << ", p = " << p ;
        } else if(set_density_temperature) {
          s << ", T = " << T ;
        }
      } else if(set_pressure_temperature) {
        s << " p = " << p << ", T = " << T ;
      } else {
        cerr << "fluidState incomprehensible in Print" << endl ;
      }
      
      if(use_mach_velocity)
        s << ", M = " ;
      else
        s << ", u = " ;
      s << "["<<u.x<<","<<u.y<<","<<u.z<<"]" ;
      if(mixture.size() != 0) {
        s << ", mixture = [" ;
        map<string,real>::const_iterator mi = mixture.begin() ;
        s << mi->first << " = " << mi->second ;
        ++mi ;
        for(;mi!=mixture.end();++mi) {
          s << ", " << mi->first << " = " << mi->second ;
        }
        s << "]" ;
      }
      if(set_k)
        s << ",k = " << k ;
      if(set_e)
        s << ",e = " << e ;
      if(set_w)
        s << ",w = " << w ;
      if(set_R)
        s << ",R = " << R ;
      if(set_nu_t)
        s << ",nu_t = " << nu_t ;
      if(set_nu_tm)
        s << ",nu_tm = " << nu_tm ;
      if(set_tmu)
        s << ",tmu = " << tmu ;
      if(equilibrium_mixture)
        s << ", equilibrium" ;
      
      s << " >" << endl ;
            
    } else {
      s << "<>" << endl ;
    }
    return s ;
  }
}
