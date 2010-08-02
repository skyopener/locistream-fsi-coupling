#ifndef QVI_H
#define QVI_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <Tools/debug.h>
#include "sciTypes.h"

namespace streamUns {
  class conservativeVectorInfo {
    friend std::ostream &operator<<(std::ostream &s,
                                    const conservativeVectorInfo &qvi) ;
    friend std::istream &operator>>(std::istream &s,
                                    conservativeVectorInfo &qvi) ;
    int num_species ;
    int momentum_index ;
    int total_energy_index ;
    int num_non_equilibrium_energy ;
    int non_equilibrium_energy_index ;
    int turbulence_scalar_index ;
    int num_turbulence_scalar ;
    int turbulence_energy_index ;
    int num_turbulence_energy ;
    int les_sgk_index ; //index for subgrid kinetic energy k.
    int num_les_sgk ;
    int num_variables ;
    std::vector<std::string> species_names ;
  public:
    //default constructor for class conservativeVectorInfo
    conservativeVectorInfo() {
      std::vector<std::string> nl ;
      nl.push_back("DEFAULT") ;
      initialize(1,0,nl,0,0,0) ;  //ns=1, nee=0
    } 
    conservativeVectorInfo(int ns, int nee, const std::vector<std::string> &sn,
                           int nts, int nte, int nsgk){
      initialize(ns,nee,sn,nts,nte,nsgk) ;
    } 
    //The following function does initialization of conservativeVectorInfo
    void initialize(int ns,int nee,const std::vector<std::string> &sn,
                    int nts,int nte,int nsgk) {
      num_species = ns ;
      momentum_index = num_species ;
      total_energy_index = momentum_index + 3 ;
      non_equilibrium_energy_index = total_energy_index+1 ;
      num_non_equilibrium_energy = nee ;
      turbulence_scalar_index =non_equilibrium_energy_index+nee ;
      num_turbulence_scalar = nts ;
      turbulence_energy_index = turbulence_scalar_index+nts ;
      num_turbulence_energy = nte ;
      les_sgk_index = turbulence_energy_index + nte ;
      num_les_sgk = nsgk ;
      num_variables = num_species + num_non_equilibrium_energy + 
        num_turbulence_scalar + num_turbulence_energy + num_les_sgk + 4 ;
      species_names = sn ;
      warn(static_cast<int>(sn.size()) != ns) ;  //Sanity check
    }
    int numSpecies() const  { return num_species ; }
    int numTurbulenceScalar() const { return num_turbulence_scalar ; }
    int numTurbulenceEnergy() const { return num_turbulence_energy ; }
    int numLESsgk() const { return num_les_sgk ; }
    //The following function returns the species index for a specific species
    int speciesIndex(std::string name) const {
      for(int i=0;i<num_species;++i)
        if(species_names[i] == name)
          return i ;
      return -1 ;
    }
    int momentumIndex() const { return momentum_index ; }
    int totalEnergyIndex() const { return total_energy_index ; }
    int turbulenceScalarIndex() const { return turbulence_scalar_index ; }
    int turbulenceEnergyIndex() const { return turbulence_energy_index ; } 
    int lessgkIndex() const { return les_sgk_index ; } 
    int vectorSize() const { return num_variables ; }
    //Return species name for species i
    const std::string &speciesName(int i) const { return species_names[i] ; }
  } ;

  //Overload ostream operator for class conservativeVectorInfo. This outputs
  //number of species, number of non-equilibrium energy, number of turbulence 
  //variables and all species names
  inline std::ostream &operator<<(std::ostream &s,
                                  const conservativeVectorInfo &qvi)
  {
    s << " " << qvi.num_species << " " << qvi.num_non_equilibrium_energy
      << " " <<qvi.num_turbulence_scalar<<" "<<
      qvi.num_turbulence_energy << " " << qvi.num_les_sgk << std::endl ;
    for(int i=0;i<qvi.num_species;++i)
      s << qvi.species_names[i] << std::endl ;
    return s ;
  }

  //Overload istream operator for class conservativeVectorInfo. 
  inline std::istream &operator>>(std::istream &s, conservativeVectorInfo &qvi)
  {
    int ns,nee,nts,nte,nsgk ;
    s >> ns >> nee >>nts >>nte >> nsgk ; //input number of species, 
    //number of non-equilibrium 
    //energy and number os turbulence variables 
    std::vector<std::string> nl ;
    for(int i=0;i<ns;++i) {
      std::string nm ; 
      s >> nm ;    //input name of the species
      nl.push_back(nm) ;  //put name of the species into the vector nl which 
      //                  //holds the names
    }
    qvi.initialize(ns,nee,nl,nts,nte,nsgk) ;
    return s ;
  }

  namespace conservative_vector {

    
    inline real getRho(const real *q,
		       const conservativeVectorInfo &qvi) {
      const int ns = qvi.numSpecies() ;
      real rho = 0.0 ;
      for(int i=0;i<ns;++i) 
	rho += q[i] ;
      return rho ;
    }

    inline void getMixture(real *mixture,const real *q,
			   const conservativeVectorInfo &qvi) {
      const int num_species = qvi.numSpecies() ;
      
      const real rho = getRho(q,qvi) ;
      for(int i=0;i<num_species;++i)
	mixture[i] = q[i]/rho ;
    }

    //Get state variables from conservative variable q
    template<class EOS>
    inline typename EOS::State getState(const real *q,
                                        const conservativeVectorInfo &qvi,
                                        const EOS &eos, const float *hint = 0) {
      const int mi = qvi.momentumIndex() ;
      const int ei = qvi.totalEnergyIndex() ;
      const int ns = qvi.numSpecies() ;

      const real &ru = q[mi+0] ;
      const real &rv = q[mi+1] ;
      const real &rw = q[mi+2] ;
      const real &re0 = q[ei] ;

      const real rho = getRho(q,qvi) ;
      const real rr = 1./rho ;
      
      tmp_array<real> rhoi(ns) ;

      getMixture(rhoi,q,qvi) ;
      for(int i=0;i<ns;++i)
        rhoi[i] = rho*rhoi[i] ;
      
      real re = re0 - 0.5*(ru*ru+rv*rv+rw*rw)*rr ;

      return eos.State_from_rho_e(rhoi,re*rr,hint) ;
    }

    template<class EOS>
    inline typename EOS::State getState(const real *q,
                                        const conservativeVectorInfo &qvi,
                                        const EOS &eos,
                                        real *ms,
                                        const float *hint=0) {
      const int mi = qvi.momentumIndex() ;
      const int ei = qvi.totalEnergyIndex() ;
      const int ns = qvi.numSpecies() ;

      const real &ru = q[mi+0] ;
      const real &rv = q[mi+1] ;
      const real &rw = q[mi+2] ;
      const real &re0 = q[ei] ;

      const real rho = getRho(q,qvi) ;
      const real rr = 1./rho ;

      tmp_array<real> rhoi(ns) ;

      getMixture(rhoi,q,qvi) ;
      for(int i=0;i<ns;++i)
        rhoi[i] = rho*rhoi[i] ;
      
      real re = re0 - 0.5*(ru*ru+rv*rv+rw*rw)*rr ;

      return eos.State_from_rho_e(rhoi,re*rr,ms,hint) ;
    }

    //convert conservative variables to primitive variables
    //(uses eos_state information that must be computed before calling this
    // routine)

    template<class EOS> inline void
    conserv2primitive(real *qp,real Pambient, const real *q,
                      const conservativeVectorInfo &qvi,
                      const EOS &eos, const typename EOS::State &eos_state) {
      const int ns = qvi.numSpecies() ;
      const int nts = qvi.numTurbulenceScalar() ;
      const int nte = qvi.numTurbulenceEnergy() ;
      const int mi = qvi.momentumIndex() ;
      const int ei = qvi.totalEnergyIndex() ;
      const int tsi = qvi.turbulenceScalarIndex() ;
      const int tei = qvi.turbulenceEnergyIndex() ;
      const int lsgki = qvi.lessgkIndex() ;
      const int nsgk = qvi.numLESsgk() ;

      const real rho = eos_state.density() ;
      const real rhor =  1.0/rho ;
      tmp_array<real> mixture(eos.numSpecies()) ;
      conservative_vector::getMixture(mixture,q,qvi) ;
      const real pg = eos_state.pressure() - Pambient ;
        

      for(int i=0;i<ns;++i)
        qp[i] = rho*mixture[i] ;  //species density
      for(int i=0;i<3;++i)
        qp[mi+i] = q[mi+i]*rhor ;  //velocities in three directions
      qp[ei] = pg ; //pressure
    }

    inline void getDensity(real *density, const real *q,
			   const conservativeVectorInfo &qvi) {
      const int num_species = qvi.numSpecies() ;
 
      for(int i=0; i<num_species;++i) 
	density[i] = q[i] ;
    }

    
    //Temperature derivative with respect to conservative variables (this is 
    //needed when compute chemistry flux Jacobian since chemistry source term is
    //the function of both density and temperature)
    template<class EOS> inline void
    dTdQ(const EOS &eos,
         const typename EOS::State &s, const real *ms,
         const real *q, const conservativeVectorInfo &qvi,
         real *dtdq) {
      const int ns = qvi.numSpecies() ;
      const int mi = qvi.momentumIndex() ;
      const int ei = qvi.totalEnergyIndex() ;
      const real &ru = q[mi+0] ;
      const real &rv = q[mi+1] ;
      const real &rw = q[mi+2] ;
      const real rho = s.density() ;
      const real rhor = 1.0/rho ;
      const real rho_KE = 0.5*(ru*ru+rv*rv+rw*rw)*rhor ;

      const real rrcvt = 1.0/(s.rho_cvt()) ; ;
      const real qh = rho_KE*rhor ;
      tmp_array<real> species_ei(ns) ;
      eos.get_ei(s,ms,species_ei) ;
      for(int i=0;i<ns;++i) 
        dtdq[i] = rrcvt*(qh - species_ei[i]) ; //partial T /partial rho
      for(int i=mi;i<mi+3;++i)
        dtdq[i] = -q[i]*rrcvt*rhor ;  //partial T /partial rho*u
      dtdq[ei] = rrcvt ;             //partial T/ partial rho*e
    }

#ifdef OLD_STUFF

    template<class EOS> inline void
    dpdQ(const EOS &eos,
         const typename EOS::State s, const rael *ms,
         const real *q, const conservativeVectorInfo &qvi, real *dpdq) {
      const int ns = qvi.numSpecies() ;
      const int mi = qvi.momentumIndex() ;
      const int ei = qvi.totalEnergyIndex() ;

      const real rho = s.density() ;
      const real rhor = 1.0/rho ;
      const real &u = q[mi+0]*rhor ;
      const real &v = q[mi+1]*rhor ;
      const real &w = q[mi+2]*rhor ;
      const real V2 = u*u+v*v+w*w ;
      const real gm1 = s.Gamma() - 1.0 ;
      const real T = s.temperature() ;
      tmp_array<real> species_ei(ns) ;
      eos.get_ei(s,ms,species_ei) ;
      for(int i=0;i<ns;++i)
        dpdq[i] = gm1*V2/2.0 + eos.speciesR(i)*T - gm1*species_ei[i] ;
      dpdq[mi+0] = - u * gm1 ;
      dpdq[mi+1] = - v * gm1 ;
      dpdq[mi+2] = - w * gm1 ;
      dpdq[ei] = gm1 ;
    }


    template<class EOS> inline void
    dgamdQ(const EOS &eos, const typename EOS::State s,
           const real *ms,
           const real *species_dcvdT,
           const real *dtdq,
           const real *q,
           const conservativeVectorInfo &qvi, real *dgamdQ) {
      const int ns = qvi.numSpecies() ;
      const int vs = qvi.vectorSize() ;

    
      const real rcv  = s.rho_cvt() ;
      const real rrcv = 1./rcv ;
      const real gam = s.Gamma() ;
      const real gm1 = gam - 1.0 ;

      real rho_dcvdT = 0 ;
      for(int i=0;i<ns;++i) 
        rho_dcvdT += q[i]*species_dcvdT[i] ;

      const real term = -gm1*rho_dcvdT*rrcv ;
      for(int i=0;i<vs;++i) 
        dgamdQ[i] = term*dtdq[i] ;

      for(int i=0;i<ns;++i) {
        const real Ri = eos.speciesR(i) ;
        const real cvi = species_cvi[i] ;
        dgamdQ[i] += (Ri-gm1*cvi)*rrcv ;
      }

    }
    
    template<class EOS> inline void
    dadQ(const EOS &eos, const typename EOS::State s,
         const real *dgamdQ,
         const real *dpdQ,
         const real *q,
         const conservativeVectorInfo &qvi,
         real *dadQ) {
      const real a = s.soundSpeed() ;
      const real gam = s.Gamma() ;
      const real p  = s.pressure() ;
      const real rho = s.density() ;
      const int vs = qvi.vectorSize() ;
      const int ns = qvi.numSpecies() ;
      real coef = .5/(rho*a) ;
      for(int i=0;i<vs;++i) 
        dadQ[i] = coef*(p*dgamdQ[i]+gam*dpdQ[i]) ;
      for(int i=0;i<ns;++i)
        dadQ[i] -= coef*gam*p/rho ;
    }
#endif

  }

  namespace primitive_vector{

    template<class EOS>
    inline void create(real *qp, real pg, real T,
                       const real *mf,
                       vect3d u,
                       real Pambient,
                       const conservativeVectorInfo &qvi,
                       const EOS &eos) {
      const int mi = qvi.momentumIndex() ;
      const int ei = qvi.totalEnergyIndex() ;
      const int ns = qvi.numSpecies() ;
      const int vs = qvi.vectorSize() ;

      for(int i=ei+1;i<vs;++i)
        qp[i] = 0 ;
      qp[mi] = u.x ;
      qp[mi+1] = u.y ;
      qp[mi+2] = u.z ;
      qp[ei] = pg ;
      typename EOS::State s = eos.State_from_mixture_p_T(mf,pg+Pambient,T) ;
      real rho = s.density() ;
      for(int i=0;i<ns;++i)
        qp[i] = mf[i]*rho ;
    }

    inline real getGaugePressure(const real *qp,
                                 const conservativeVectorInfo &qvi) {
      const int ei = qvi.totalEnergyIndex() ;
      return qp[ei] ;
    }      

    inline real getRho(const real *qp,
                       const conservativeVectorInfo &qvi) {
      const int ns =qvi.numSpecies() ;

      real rho =0.0 ;
      for(int i=0;i<ns;++i)
	rho += qp[i] ;
      return rho ;
    }

    inline void getMixture(real *mixture,const real *qp,
                    const conservativeVectorInfo &qvi) {
      const int num_species = qvi.numSpecies() ;
      
      const real rho = getRho(qp,qvi) ;
      for(int i=0;i<num_species;++i)
	mixture[i] = qp[i]/rho ;
    }
    
    //Get state variables from primitive variable qp
    template<class EOS>
    inline typename EOS::State getState(const real *qp,
                                        real Pambient,
                                        const conservativeVectorInfo &qvi,
                                        const EOS &eos,
                                        const float *hint=0) {
      const int ei = qvi.totalEnergyIndex() ;      

      typename EOS::State s = eos.State_from_rho_p(qp,qp[ei]+Pambient,hint) ;
      return s ;
    }

    template<class EOS>
    inline typename EOS::State getState(const real *qp, real Pambient,
                                        const conservativeVectorInfo &qvi,
                                        const EOS &eos,
                                        real *ms,
                                        const float *hint=0) {
      const int ei = qvi.totalEnergyIndex() ;

      return eos.State_from_rho_p(qp,qp[ei]+Pambient,ms,hint) ;
    }

    //convert primitive variables to conservative variables, returns eos_state
    // for the converted vector.  
    template<class EOS> inline typename EOS::State
    primitive2conserv(real *q, const real *qp, real Pambient,
                      const conservativeVectorInfo &qvi,
                      const EOS &eos,
                      real *ms, const float *hint = 0) {

      const int ns = qvi.numSpecies() ;
      typename EOS::State s = primitive_vector::getState(qp,Pambient,
                                                         qvi,eos,ms,hint) ;
      tmp_array<real> mixture(ns) ;
      primitive_vector::getMixture(mixture,qp,qvi) ;
      const real rho = s.density() ;
      for(int i=0;i<ns;++i)
        q[i] = rho*mixture[i] ;
      const int vs = qvi.vectorSize() ;
      for(int i=ns;i<vs;++i)
        q[i] = rho*qp[i] ;
      const int ei = qvi.totalEnergyIndex() ;
      const int mi = qvi.momentumIndex() ;
      vect3d u = vect3d(qp[mi+0],qp[mi+1],qp[mi+2]) ;
      q[ei] = s.rho_energy() + 0.5*dot(u,u)*rho ;
      return s ;
    }

    // New stuff for the new generalized code
    template<class EOS> inline void
    dTdq(const EOS &eos,const typename EOS::State &s,
         const real *ms,
         const conservativeVectorInfo &qvi,
         real *dtdq) {
      const int ns = qvi.numSpecies() ;
      const int vs = qvi.vectorSize() ;
      const int ei = qvi.totalEnergyIndex() ;
      for(int i=ns;i<vs;++i)
        dtdq[i] = 0 ;
      eos.dTdri(dtdq,s,ms) ;
      dtdq[ei] = eos.dTdP(s,ms);
    }

  

    template<class EOS> inline void
    dQdq(Loci::Mat<real_fj> &dQdq,
         vect3d U,
         const real *mixture,
         const typename EOS::State &s,
         const real *ms,
         const EOS &eos,
         const conservativeVectorInfo &qvi) {

    const int ns = qvi.numSpecies() ;
    const int mi = qvi.momentumIndex() ;
    const int ei = qvi.totalEnergyIndex() ;
    const int vs = qvi.vectorSize() ;

    const real rho = s.density() ;
    
    for(int i=0;i<vs;++i)
      for(int j=0;j<vs;++j)
        dQdq[i][j] = 0 ;
    for(int i=0;i<ns;++i)
      dQdq[i][i] = 1. ;

    tmp_array<real> dreidri(ns) ;
    eos.dreidri(dreidri,s,ms) ;
    const real q = 0.5*dot(U,U) ;
    for(int j=0;j<ns;++j) {
      dQdq[ei][j] = q + dreidri[j] ;
      dQdq[mi+0][j] = U.x ;
      dQdq[mi+1][j] = U.y ;
      dQdq[mi+2][j] = U.z ;
    }

    dQdq[mi+0][mi+0] = rho ;
    dQdq[mi+1][mi+1] = rho ;
    dQdq[mi+2][mi+2] = rho ;

    
    dQdq[ei][mi+0] = rho*U.x ;
    dQdq[ei][mi+1] = rho*U.y ;
    dQdq[ei][mi+2] = rho*U.z ;


    const real dreidP = eos.dreidP(s,ms) ;
    dQdq[ei][ei] = dreidP ;

  }
  }
}

namespace Loci {

  template<> struct data_schema_traits<streamUns::conservativeVectorInfo> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<streamUns::conservativeVectorInfo> Converter_Type ;
  } ;
  
}

#endif
