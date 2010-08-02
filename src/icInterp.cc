#include <Loci.h>
#include "sciTypes.h"
#include "eos.h"
using fluidPhysics::EOS ;
#include "fluidState.h"
#include "reaction.h"
#include "initialCondition.h"
#include "interpolateFile.h"
#include <Tools/parse.h>
#include <vector>
#include <string>

using std::vector ;
using std::string ;

namespace streamUns {


  /*class interpolate_ic : public pointwise_rule {
    const_store<vect3d> cellcenter ; //position of nodes
    const_param<conservativeVectorInfo> qvi ;
    const_param<EOS> eos ;
    const_param<reaction> reactor ;
    const_param<string> initialConditionsFile ;
    const_param<real> Pambient ;
    store<real> k_ic ;
    store<real> tmuu_ic ;
    store<real> rho_ic ;
    storeVec<real> q_ic ;
    store<real> pg_ic ;
  public:
    interpolate_ic() ;
    virtual void compute(const sequence &seq) ;
  } ;

  interpolate_ic::interpolate_ic() {
    name_store("initialConditionsFile",initialConditionsFile) ;
    name_store("cellcenter",cellcenter) ;
    name_store("ic_file::q_ic",q_ic) ;
    name_store("ic_file::pg_ic",pg_ic) ;
    name_store("Pambient",Pambient) ;
    name_store("eos",eos) ;
    name_store("qvi",qvi) ;
    name_store("reactor",reactor) ;
    name_store("ic_file::tmuu_ic",tmuu_ic) ;
    name_store("ic_file::k_ic",k_ic) ;
    name_store("ic_file::rho_ic",rho_ic) ;
    input("Pambient") ;
    input("eos,qvi,reactor") ;
    input("cellcenter") ;
    input("initialConditionsFile") ;
    output("ic_file::q_ic") ;
    output("ic_file::k_ic,ic_file::tmuu_ic,ic_file::rho_ic") ;
    output("ic_file::pg_ic") ;
    constraint("geom_cells") ;
    disable_threading();
  }


  //apply initial conditions over sequence (cell by cell)
  void interpolate_ic::compute(const sequence &seq) {
    if(Loci::GLOBAL_AND(seq==EMPTY))
      return ;
    const int qs = qvi->vectorSize() ;
    const int ns = qvi->numSpecies() ;

    q_ic.setVecSize(qs) ;


    store<vect3d> loc,vel ;
    store<real> p,T,k,mu_t ;
    storeVec<real> mix ;
    vector<int> mixids ;

    string filename = *initialConditionsFile ;
    if(Loci::MPI_rank == 0 )
      cout << "interpolating initial conditions from " << filename
           << endl ;
    read_puT_file(filename,*qvi,loc,p,T,vel,k,mu_t,mix,mixids) ;

    int min_mix = mixids[0] ;
    for(size_t i=1;i<mixids.size();++i)
      min_mix = min(min_mix,mixids[i]) ;
    if(min_mix < 0) {
      cerr << "interpolated initial conditions contain invalid species"
           << endl ;
    }


    int npnts = loc.domain().Max()+1 ;
    vector<Loci::kdTree::coord3d> pnts(npnts) ;
    vector<int> pnt_id (npnts) ;
    for(int i=0;i<npnts;++i) {
      pnts[i][0] = loc[i].x ;
      pnts[i][1] = loc[i].y ;
      pnts[i][2] = loc[i].z ;
      pnt_id[i] = i ;
    }

    Loci::kdTree::kd_tree kd(pnts,pnt_id) ;
    tmp_array<real> mixt(mix.vecSize()) ;
    tmp_array<real> mixture(ns) ;

    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si) {
      vector<int> neighbors = get_stencil(kd,cellcenter[*si]) ;

      vector<double> w  ;
      stencil_weights(w,neighbors,loc,cellcenter[*si]) ;

      vect3d ur(0,0,0) ;
      real pr = 0, Tr = 0 ;
      real kr = 0 ;
      real mu_tr = 0 ;
      for(int i=0;i<mix.vecSize();++i)
        mixt[i] = 0 ;

      for(size_t i=0;i<neighbors.size();++i) {
        int pt = neighbors[i] ;

        double weight = w[i] ;
        ur += weight*vel[pt] ;
        pr += weight*p[pt] ;
        Tr += weight*T[pt] ;
        kr += weight*k[pt] ;
        mu_tr += weight*mu_t[pt] ;

        for(int j=0;j<mix.vecSize();++j) {
          mixt[j] += weight*mix[pt][j] ;
        }
      }

      for(int i=0;i<ns;++i)
        mixture[i] = 0.0 ;
      for(int i=0;i<mix.vecSize();++i)
        mixture[mixids[i]] = mixt[i] ;

      // calculate the conservertive fluid state for face id ;
      real *qout = q_ic[*si] ;
      const vect3d U = ur ;

      if(ns == 1)
        mixture[0] = 1.0 ;

      EOS::State s = eos->State_from_mixture_p_T(mixture,pr,Tr) ;

      for(int i=0;i<ns;++i)
        qout[i] = mixture[i]*s.density() ;
      const int mi = qvi->momentumIndex() ;
      qout[mi+0] = U.x*s.density() ;
      qout[mi+1] = U.y*s.density() ;
      qout[mi+2] = U.z*s.density() ;
      const int ei = qvi->totalEnergyIndex() ;
      qout[ei] = s.rho_energy()+.5*dot(U,U)*s.density() ;
      rho_ic[*si] = s.density() ;

      k_ic[*si] = kr ;
      tmuu_ic[*si] = mu_tr ;
      pg_ic[*si] = pr - *Pambient ;
    }

  }

  register_rule<interpolate_ic> register_interpolate_ic ;*/


  /*class interpolate_ic_noslip : public pointwise_rule {
    const_Map min_cell2noslip ;
    const_store<vect3d> facecenter ;
    const_store<Area> area ;
    
    const_store<vect3d> cellcenter ; //position of nodes
    const_param<conservativeVectorInfo> qvi ;
    const_param<EOS> eos ;
    const_param<reaction> reactor ;
    const_param<string> initialConditionsFile ;
    const_param<real> Pambient ;
    store<real> k_ic ;
    store<real> tmuu_ic ;
    store<real> rho_ic ;
    storeVec<real> q_ic ;
    store<real> pg_ic ;
  public:
    interpolate_ic_noslip() ;
    virtual void compute(const sequence &seq) ;
  } ;

  interpolate_ic_noslip::interpolate_ic_noslip() {
    name_store("initialConditionsFile",initialConditionsFile) ;
    name_store("cellcenter",cellcenter) ;
    name_store("noslip::ic_file::q_ic",q_ic) ;
    name_store("noslip::ic_file::pg_ic",pg_ic) ;
    name_store("Pambient",Pambient) ;
    name_store("eos",eos) ;
    name_store("qvi",qvi) ;
    name_store("reactor",reactor) ;
    name_store("noslip::ic_file::tmuu_ic",tmuu_ic) ;
    name_store("noslip::ic_file::k_ic",k_ic) ;
    name_store("noslip::ic_file::rho_ic",rho_ic) 
;
    name_store("min_cell2noslip",min_cell2noslip) ;
    name_store("facecenter",facecenter) ;
    name_store("area",area) ;
    input("min_cell2noslip->(facecenter,area)") ;
    input("Pambient") ;
    input("eos,qvi,reactor") ;
    input("cellcenter") ;
    input("initialConditionsFile") ;
    output("noslip::ic_file::q_ic") ;
    output("noslip::ic_file::k_ic,noslip::ic_file::tmuu_ic,noslip::ic_file::rho_ic") ;
    output("noslip::ic_file::pg_ic") ;
    constraint("geom_cells") ;
    disable_threading();
  }


  //apply initial conditions over sequence (cell by cell)
  void interpolate_ic_noslip::compute(const sequence &seq) {
    if(Loci::GLOBAL_AND(seq==EMPTY))
      return ;
    const int qs = qvi->vectorSize() ;
    const int ns = qvi->numSpecies() ;

    q_ic.setVecSize(qs) ;


    store<vect3d> loc,vel ;
    store<real> p,T,k,mu_t ;
    storeVec<real> mix ;
    vector<int> mixids ;

    string filename = *initialConditionsFile ;
    if(Loci::MPI_rank == 0 )
      cout << "interpolating initial conditions from " << filename
           << endl ;
    read_puT_file(filename,*qvi,loc,p,T,vel,k,mu_t,mix,mixids) ;

    int min_mix = mixids[0] ;
    for(size_t i=1;i<mixids.size();++i)
      min_mix = min(min_mix,mixids[i]) ;
    if(min_mix < 0) {
      cerr << "interpolated initial conditions contain invalid species"
           << endl ;
    }


    int npnts = loc.domain().Max()+1 ;
    vector<Loci::kdTree::coord3d> pnts(npnts) ;
    vector<int> pnt_id (npnts) ;
    for(int i=0;i<npnts;++i) {
      pnts[i][0] = loc[i].x ;
      pnts[i][1] = loc[i].y ;
      pnts[i][2] = loc[i].z ;
      pnt_id[i] = i ;
    }

    Loci::kdTree::kd_tree kd(pnts,pnt_id) ;
    tmp_array<real> mixt(mix.vecSize()) ;
    tmp_array<real> mixture(ns) ;

    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si) {
      vector<int> neighbors = get_stencil(kd,cellcenter[*si]) ;

      // remove any points on the wrong side of a viscous wall
      vect3d fcenter = facecenter[min_cell2noslip[*si]] ;
      vect3d n = area[min_cell2noslip[*si]].n ;
      vector<int> nn ;
      for(size_t i=0;i<neighbors.size();++i) 
        if(dot(n,loc[neighbors[i]]-fcenter) <= 0.0)
          nn.push_back(neighbors[i]) ;
      if(nn.size() != 0)
        neighbors.swap(nn) ;

      // compute stencil weights
      vector<double> w  ;
      stencil_weights(w,neighbors,loc,cellcenter[*si]) ;

      vect3d ur(0,0,0) ;
      real pr = 0, Tr = 0 ;
      real kr = 0 ;
      real mu_tr = 0 ;
      for(int i=0;i<mix.vecSize();++i)
        mixt[i] = 0 ;

      for(size_t i=0;i<neighbors.size();++i) {
        int pt = neighbors[i] ;

        double weight = w[i] ;
        ur += weight*vel[pt] ;
        pr += weight*p[pt] ;
        Tr += weight*T[pt] ;
        kr += weight*k[pt] ;
        mu_tr += weight*mu_t[pt] ;

        for(int j=0;j<mix.vecSize();++j) {
          mixt[j] += weight*mix[pt][j] ;
        }
      }

      for(int i=0;i<ns;++i)
        mixture[i] = 0.0 ;
      for(int i=0;i<mix.vecSize();++i)
        mixture[mixids[i]] = mixt[i] ;

      // calculate the conservertive fluid state for face id ;
      real *qout = q_ic[*si] ;
      const vect3d U = ur ;

      if(ns == 1)
        mixture[0] = 1.0 ;

      EOS::State s = eos->State_from_mixture_p_T(mixture,pr,Tr) ;

      for(int i=0;i<ns;++i)
        qout[i] = mixture[i]*s.density() ;
      const int mi = qvi->momentumIndex() ;
      qout[mi+0] = U.x*s.density() ;
      qout[mi+1] = U.y*s.density() ;
      qout[mi+2] = U.z*s.density() ;
      const int ei = qvi->totalEnergyIndex() ;
      qout[ei] = s.rho_energy()+.5*dot(U,U)*s.density() ;
      rho_ic[*si] = s.density() ;

      k_ic[*si] = kr ;
      tmuu_ic[*si] = mu_tr ;
      pg_ic[*si] = pr - *Pambient ;
    }

  }

  register_rule<interpolate_ic_noslip> register_interpolate_ic_noslip ;*/
  
  /*class output_puT_file : public pointwise_rule {
    const_storeVec<real> q_ic ;
    const_store<real> pg_ic ;
    const_param<real> Pambient ;
    const_store<vect3d> cellcenter ; //position of nodes
    const_store<vec<2> > sst_q_ic ;
    const_param<EOS> eos ;
    const_param<conservativeVectorInfo> qvi ;
    store<real> dummy ;
  public:
    output_puT_file() {
      name_store("q_ic",q_ic) ;
      name_store("pg_ic",pg_ic) ;
      name_store("Pambient",Pambient) ;
      name_store("cellcenter",cellcenter) ;
      name_store("sst_q_ic",sst_q_ic) ;
      name_store("eos",eos) ;
      name_store("qvi",qvi) ;
      name_store("puT_file",dummy) ;
      input("q_ic,pg_ic,Pambient,cellcenter,sst_q_ic,eos,qvi") ;
      output("puT_file") ;
      disable_threading() ;
    }
    void compute(const sequence &seq) ;
  } ;

  void output_puT_file::compute(const sequence &seq) {
    ofstream ofile("puT.dat",ios::out) ;
    ofile.precision(16) ;

    const int nsp = eos->numSpecies() ;
    ofile << nsp << endl ;
    for(int i=0;i<nsp;++i)
      ofile << eos->speciesNames()[i] << endl ;

    int numpnts = seq.size() ;
    ofile << numpnts << endl ;
    const int mi=qvi->momentumIndex() ;
    tmp_array<real> mf(nsp) ;
    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si) {
      Entity e = *si ;

      //      EOS::State s = conservative_vector::getState(q_ic[e],*qvi,*eos) ;
      EOS::State s = eos->State_from_rho_p(q_ic[e],pg_ic[e]+*Pambient) ;
      const real r = s.density() ;
      const real T = s.temperature() ;
      const real P = s.pressure() ;
      const real rr = 1./r ;
      const real k = sst_q_ic[e][0]*rr ;
      const real w = sst_q_ic[e][1]*rr ;
      const real tmuu = r*k/(w+1e-20) ;
      ofile << cellcenter[e].x << ' '
            << cellcenter[e].y << ' '
            << cellcenter[e].z << ' '
            << q_ic[e][mi+0]*rr << ' '
            << q_ic[e][mi+1]*rr << ' '
            << q_ic[e][mi+2]*rr << ' '
            << P << ' '
            << T << ' '
            << k << ' '
            << tmuu << ' ' ;
      conservative_vector::getMixture(mf,q_ic[e],*qvi) ;
      for(int i=0;i<nsp;++i)
        ofile << mf[i] << ' ' ;
      ofile << endl ;
    }
  }

  register_rule<output_puT_file> register_output_puT ;*/

  class scalable_puToutput : public pointwise_rule {
    const_param<string> turbulence_model ;
    const_store<vect3d> cellcenter,u ;
    const_store<real> temperature,p ;
    const_storeVec<real> mixture ;
    const_param<real> Pambient ;
    const_param<EOS> eos ;
    const_param<int> ncyc ;
    const_param<int> ncycle ;
    const_param<int> restart_modulo ;
    const_param<string> modelName ;
    param<bool> OUTPUT ;

  public:
    scalable_puToutput() ;
    virtual void compute(const sequence &seq) ;
  } ;

  scalable_puToutput::scalable_puToutput() {
    name_store("cellcenter{n,it}",cellcenter) ;
    name_store("v{n,it}",u) ;
    name_store("Pambient{n,it}",Pambient) ;
    name_store("y{n,it}",mixture) ;
    name_store("eos{n,it}",eos) ;
    name_store("p{n,it}",p) ;
    name_store("temperature{n,it}",temperature) ;
    name_store("timeStepNumber{n}",ncyc) ;
    name_store("ncycle{n}",ncycle) ;
    name_store("restart_modulo{n,it}",restart_modulo) ;
    name_store("modelName{n,it}",modelName) ;
    name_store("turbulence_model{n,it}",turbulence_model) ;
    input("turbulence_model{n,it}") ;
    name_store("OUTPUT{n,it}",OUTPUT) ;
    conditional("do_restart{n,it}") ;
    input("cellcenter{n,it},v{n,it},Pambient{n,it},p{n,it},temperature{n,it}");
    input("y{n,it}") ;
    input("eos{n,it}") ;
    input("modelName{n,it},ncycle{n},restart_modulo{n,it},timeStepNumber{n}") ;
    constraint("geom_cells{n,it},laminarFlow{n,it}") ;
    output("OUTPUT{n,it}") ;

  }

  void scalable_puToutput::compute(const sequence &seq) {
    ostringstream oss ; int cycle=*ncycle ;
    if(*restart_modulo != 0) cycle = cycle % *restart_modulo ;
    oss << "output/put"<< "." << cycle << "_" << *modelName ;
    string filename = oss.str() ;

    hid_t file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                         H5P_DEFAULT, H5P_DEFAULT) ;

    store<vect3d> vtmp ;
    store<real> stmp ;
    vtmp.allocate(entitySet(seq)) ;
    stmp.allocate(entitySet(seq)) ;

    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      vtmp[*si] = cellcenter[*si] ;
    Loci::writeContainer(file_id,"pos",vtmp.Rep()) ;
    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      vtmp[*si] = u[*si] ;
    Loci::writeContainer(file_id,"u",vtmp.Rep()) ;
    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      stmp[*si] = p[*si]-(*Pambient) ;
    Loci::writeContainer(file_id,"pg",stmp.Rep()) ;
    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      stmp[*si] = temperature[*si] ;
    Loci::writeContainer(file_id,"t",stmp.Rep()) ;

    Loci::writeContainer(file_id,"Pambient",Pambient.Rep()) ;
    Loci::writeContainer(file_id,"turbulence_model",turbulence_model.Rep()) ;

    string speciesNames ;
    speciesNames += eos->speciesNames()[0] ;
    for(int i=1;i<eos->numSpecies();++i) {
      speciesNames += ':' ;
      speciesNames += eos->speciesNames()[i] ;
    }
    param<std::string> sn ;
    *sn = speciesNames ;
    Loci::writeContainer(file_id,"speciesNames",sn.Rep()) ;

    if(eos->numSpecies() > 1) {
      storeVec<real> mtmp ;
      mtmp.allocate(entitySet(seq)) ;
      mtmp.setVecSize(mixture.vecSize()) ;
      for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
        mtmp[*si] = mixture[*si] ;
      Loci::writeContainer(file_id,"mixture",mtmp.Rep()) ;
    }

    Loci::hdf5CloseFile(file_id) ;

  }
  register_rule<scalable_puToutput> register_scalable_puToutput ;

  class scalable_output_kOmega : public pointwise_rule {
    const_param<string> turbulence_model ;
    const_store<vect3d> cellcenter,u ;
    const_store<real> temperature,p,k,w,rho ;
    const_storeVec<real> mixture ;
    const_param<real> Pambient ;
    const_param<EOS> eos ;
    const_param<int> ncyc ;
    const_param<int> ncycle ;
    const_param<int> restart_modulo ;
    const_param<string> modelName ;
    param<bool> OUTPUT ;

  public:
    scalable_output_kOmega() ;
    virtual void compute(const sequence &seq) ;
  } ;

  scalable_output_kOmega::scalable_output_kOmega() {
    name_store("cellcenter{n,it}",cellcenter) ;
    name_store("v{n,it}",u) ;
    name_store("Pambient{n,it}",Pambient) ;
    name_store("y{n,it}",mixture) ;
    name_store("eos{n,it}",eos) ;
    name_store("p{n,it}",p) ;
    name_store("temperature{n,it}",temperature) ;
    name_store("timeStepNumber{n}",ncyc) ;
    name_store("ncycle{n}",ncycle) ;
    name_store("restart_modulo{n,it}",restart_modulo) ;
    name_store("modelName{n,it}",modelName) ;
    name_store("turbulence_model{n,it}",turbulence_model) ;
    name_store("k{n,it}",k) ;
    name_store("omega{n,it}",w) ;
    name_store("rho{n,it}",rho) ;
    input("turbulence_model{n,it},rho{n,it},k{n,it},omega{n,it}") ;
    name_store("OUTPUT{n,it}",OUTPUT) ;
    conditional("do_restart{n,it}") ;
    input("cellcenter{n,it},v{n,it},Pambient{n,it},p{n,it},temperature{n,it}");
    input("y{n,it}") ;
    input("eos{n,it}") ;
    input("modelName{n,it},ncycle{n},restart_modulo{n,it},timeStepNumber{n}") ;
    constraint("geom_cells{n,it},menterTurbulenceModel{n,it}") ;
    output("OUTPUT{n,it}") ;

  }

  void scalable_output_kOmega::compute(const sequence &seq) {
    ostringstream oss ; int cycle=*ncycle ;
    if(*restart_modulo != 0) cycle = cycle % *restart_modulo ;
    oss << "output/put"<< "." << cycle << "_" << *modelName ;
    string filename = oss.str() ;

    hid_t file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                         H5P_DEFAULT, H5P_DEFAULT) ;

    store<vect3d> vtmp ;
    store<real> stmp ;
    vtmp.allocate(entitySet(seq)) ;
    stmp.allocate(entitySet(seq)) ;

    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      vtmp[*si] = cellcenter[*si] ;
    Loci::writeContainer(file_id,"pos",vtmp.Rep()) ;
    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      vtmp[*si] = u[*si] ;
    Loci::writeContainer(file_id,"u",vtmp.Rep()) ;
    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      stmp[*si] = p[*si]-(*Pambient) ;
    Loci::writeContainer(file_id,"pg",stmp.Rep()) ;
    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      stmp[*si] = temperature[*si] ;
    Loci::writeContainer(file_id,"t",stmp.Rep()) ;

    Loci::writeContainer(file_id,"Pambient",Pambient.Rep()) ;
    Loci::writeContainer(file_id,"turbulence_model",turbulence_model.Rep()) ;

    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      stmp[*si] = k[*si] ;
    Loci::writeContainer(file_id,"k",stmp.Rep()) ;

    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      stmp[*si] = rho[*si]*max(k[*si],1e-30)/w[*si] ;
    
    Loci::writeContainer(file_id,"tmuu",stmp.Rep()) ;
    
    string speciesNames ;
    speciesNames += eos->speciesNames()[0] ;
    for(int i=1;i<eos->numSpecies();++i) {
      speciesNames += ':' ;
      speciesNames += eos->speciesNames()[i] ;
    }
    param<std::string> sn ;
    *sn = speciesNames ;
    Loci::writeContainer(file_id,"speciesNames",sn.Rep()) ;

    if(eos->numSpecies() > 1) {
      storeVec<real> mtmp ;
      mtmp.allocate(entitySet(seq)) ;
      mtmp.setVecSize(mixture.vecSize()) ;
      for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
        mtmp[*si] = mixture[*si] ;
      Loci::writeContainer(file_id,"mixture",mtmp.Rep()) ;
    }

    Loci::hdf5CloseFile(file_id) ;

  }
  register_rule<scalable_output_kOmega> register_scalable_output_kOmega ;

  /*class scalable_output_SA : public pointwise_rule {
    const_param<string> turbulence_model ;
    const_store<vect3d> cellcenter,u ;
    const_store<real> temperature,gaugePressure,nu_t,rho ;
    const_storeVec<real> mixture ;
    const_param<real> Pambient ;
    const_param<EOS> eos ;
    const_param<int> ncyc ;
    const_param<int> ncycle ;
    const_param<int> plot_modulo ;
    const_param<string> modelName ;
    param<bool> OUTPUT ;

  public:
    scalable_output_SA() ;
    virtual void compute(const sequence &seq) ;
  } ;

  scalable_output_SA::scalable_output_SA() {
    name_store("cellcenter",cellcenter) ;
    name_store("u",u) ;
    name_store("Pambient",Pambient) ;
    name_store("mixture",mixture) ;
    name_store("eos",eos) ;
    name_store("gaugePressure",gaugePressure) ;
    name_store("temperature",temperature) ;
    name_store("timeStepNumber",ncyc) ;
    name_store("ncycle",ncycle) ;
    name_store("plot_modulo",plot_modulo) ;
    name_store("modelName",modelName) ;
    name_store("turbulence_model",turbulence_model) ;
    name_store("nu_t",nu_t) ;
    name_store("rho",rho) ;
    input("turbulence_model,rho,nu_t") ;
    name_store("OUTPUT",OUTPUT) ;
    conditional("do_plot") ;
    input("cellcenter,u,Pambient,gaugePressure,temperature");
    input("mixture") ;
    input("eos") ;
    input("modelName,ncycle,plot_modulo,timeStepNumber") ;
    constraint("geom_cells,Sp_All") ;
    output("OUTPUT") ;

  }

  void scalable_output_SA::compute(const sequence &seq) {
    if(*ncycle == 0)
      return ;

    ostringstream oss ;
    int cycle = *ncycle ;
    if(*plot_modulo != 0)
      cycle = cycle % *plot_modulo ;

    oss << "output/put"<< "." << cycle
        << "_" << *modelName ;
    string filename = oss.str() ;

    hid_t file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                         H5P_DEFAULT, H5P_DEFAULT) ;

    store<vect3d> vtmp ;
    store<real> stmp ;
    vtmp.allocate(entitySet(seq)) ;
    stmp.allocate(entitySet(seq)) ;

    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      vtmp[*si] = cellcenter[*si] ;
    Loci::writeContainer(file_id,"pos",vtmp.Rep()) ;
    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      vtmp[*si] = u[*si] ;
    Loci::writeContainer(file_id,"u",vtmp.Rep()) ;
    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      stmp[*si] = gaugePressure[*si] ;
    Loci::writeContainer(file_id,"pg",stmp.Rep()) ;
    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      stmp[*si] = temperature[*si] ;
    Loci::writeContainer(file_id,"t",stmp.Rep()) ;

    Loci::writeContainer(file_id,"Pambient",Pambient.Rep()) ;
    Loci::writeContainer(file_id,"turbulence_model",turbulence_model.Rep()) ;

    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      stmp[*si] = 0 ;
    Loci::writeContainer(file_id,"k",stmp.Rep()) ;

    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      stmp[*si] = rho[*si]*nu_t[*si] ;
    
    Loci::writeContainer(file_id,"tmuu",stmp.Rep()) ;
    
    string speciesNames ;
    speciesNames += eos->speciesNames()[0] ;
    for(int i=1;i<eos->numSpecies();++i) {
      speciesNames += ':' ;
      speciesNames += eos->speciesNames()[i] ;
    }
    param<std::string> sn ;
    *sn = speciesNames ;
    Loci::writeContainer(file_id,"speciesNames",sn.Rep()) ;

    if(eos->numSpecies() > 1) {
      storeVec<real> mtmp ;
      mtmp.allocate(entitySet(seq)) ;
      mtmp.setVecSize(mixture.vecSize()) ;
      for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
        mtmp[*si] = mixture[*si] ;
      Loci::writeContainer(file_id,"mixture",mtmp.Rep()) ;
    }

    Loci::hdf5CloseFile(file_id) ;

  }
  register_rule<scalable_output_SA> register_scalable_output_SA ;*/

  // This rule is used when interpolating a new incompressible case from 
  // either an incompressible or compressible puT file. Note that the density
  // is always obtained from InitialCondition since the puT file does not
  // contain density. This for incompressible interpolations, one must also
  // use initialCondition simultaneously.
  class NewInterpolateIncompressibleIC : public pointwise_rule {
    private:
      const_param<EOS> eos ;
      const_param<real> Pambient ;
      const_storeVec<real> data ;
      store<real> rho_ic ;
      store<vect3d> v_ic ;
      store<real> p_ic,T_ic,h_ic ;
      storeVec<real> y_ic ;
    public:

      // Define input and output.
      NewInterpolateIncompressibleIC() {
        name_store("ic_file::v_ic",v_ic) ;
        name_store("ic_file::p_ic",p_ic) ;
        name_store("ic_file::y_ic",y_ic) ;
        name_store("Pambient",Pambient) ;
        name_store("eos",eos) ;
        name_store("interpolateData(interpolateFile(interpolateInitialConditions),cellcenter)",data) ;
        input("interpolateData(interpolateFile(interpolateInitialConditions),cellcenter)") ;
        input("Pambient") ;
        input("eos") ;
        output("ic_file::v_ic,ic_file::p_ic,ic_file::y_ic") ;
        constraint("geom_cells,incompressibleFlow") ;
        disable_threading();
      }

      // Assign cell values.
      void calculate(Entity cc) {
        int base_size = 7 ;
        real pr = data[cc][1] ;
        vect3d ur(data[cc][2],data[cc][3],data[cc][4]) ;
        int ns = eos->numSpecies() ;
        tmp_array<real> mix(ns) ;
        if(ns == 1) mix[0]=y_ic[cc][0]=1 ;
        else for(int i=0;i<ns;++i) mix[i]=y_ic[cc][i]=data[cc][base_size+i] ;
        v_ic[cc]=ur ; p_ic[cc]=pr+(*Pambient) ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) {
        const int ns=eos->numSpecies() ; y_ic.setVecSize(ns) ;

        // Write message for user.
        if(Loci::MPI_rank==0){
          cout << "Restarting incompressible flow from puT file." << endl ;
          cout << "Using density specified in 'initialCondition'." << endl ;
        }

        do_loop(seq,this) ;
      }
  } ;

  register_rule<NewInterpolateIncompressibleIC>
    registerNewInterpolateIncompressibleIC ;

  // This rule is used when interpolating a new compressible case from a
  // compressible puT file.
  class NewInterpolateCompressibleIC : public pointwise_rule {
    private:
      const_param<EOS> eos ;
      const_param<real> Pambient ;
      const_storeVec<real> data ;
      store<real> rho_ic ;
      store<vect3d> v_ic ;
      store<real> p_ic,T_ic,h_ic ;
      storeVec<real> y_ic ;
    public:

      // Define input and output.
      NewInterpolateCompressibleIC() {
        name_store("ic_file::rho_ic",rho_ic) ;
        name_store("ic_file::v_ic",v_ic) ;
        name_store("ic_file::p_ic",p_ic) ;
        name_store("ic_file::T_ic",T_ic) ;
        name_store("ic_file::h_ic",h_ic) ;
        name_store("ic_file::y_ic",y_ic) ;
        name_store("Pambient",Pambient) ;
        name_store("eos",eos) ;
        name_store("interpolateData(interpolateFile(interpolateInitialConditions),cellcenter)",data) ;
        input("interpolateData(interpolateFile(interpolateInitialConditions),cellcenter)") ;
        input("Pambient") ;
        input("eos") ;
        output("ic_file::rho_ic,ic_file::v_ic,ic_file::p_ic,ic_file::T_ic") ;
        output("ic_file::h_ic,ic_file::y_ic") ;
        constraint("geom_cells,compressibleFlow") ;
        disable_threading();
      }

      // Assign cell values.
      void calculate(Entity cc) {
        int base_size = 7 ;
        real Tr = data[cc][0] ;
        real pr = data[cc][1] ;
        vect3d ur(data[cc][2],data[cc][3],data[cc][4]) ;
        int ns = eos->numSpecies() ;
        tmp_array<real> mix(ns) ;
        if(ns == 1) mix[0]=y_ic[cc][0]=1 ;
        else for(int i=0;i<ns;++i) mix[i]=y_ic[cc][i]=data[cc][base_size+i] ;
        EOS::State s = eos->State_from_mixture_p_T(mix,pr,Tr) ;
        v_ic[cc]=ur ; h_ic[cc]=s.enthalpy()+0.5*dot(ur,ur) ;
        rho_ic[cc] = s.density() ; p_ic[cc]=pr+(*Pambient) ; T_ic[cc]=Tr ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) {
        const int ns=eos->numSpecies() ; y_ic.setVecSize(ns) ;

        // Write message for user.
        if(Loci::MPI_rank==0){
          cout << "Restarting compressible flow from puT file." << endl ;
        }

        do_loop(seq,this) ;
      }
  } ;

  register_rule<NewInterpolateCompressibleIC>
    registerNewInterpolateCompressibleIC ;

  // This rule is used when interpolating a new compressible case from an
  // incompressible puT file. Since the incompressible puT file has zero for
  // temperature, we must use initialCondition to provide this. This rule
  // overrides the rula above.
  class NewInterpolateCompressibleFromIncompressibleIC : public pointwise_rule {
    private:
      const_param<EOS> eos ;
      const_param<real> Pambient ;
      const_param<InitialCondition> initialCondition ;
      const_storeVec<real> data ;
      store<real> rho_ic ;
      store<vect3d> v_ic ;
      store<real> p_ic,T_ic,h_ic ;
      storeVec<real> y_ic ;
    public:

      // Define input and output.
      NewInterpolateCompressibleFromIncompressibleIC() {
        name_store("priority::ic_file::rho_ic",rho_ic) ;
        name_store("priority::ic_file::v_ic",v_ic) ;
        name_store("priority::ic_file::p_ic",p_ic) ;
        name_store("priority::ic_file::T_ic",T_ic) ;
        name_store("priority::ic_file::h_ic",h_ic) ;
        name_store("priority::ic_file::y_ic",y_ic) ;
        name_store("Pambient",Pambient) ;
        name_store("eos",eos) ;
        name_store("initialCondition",initialCondition) ;
        name_store("interpolateData(interpolateFile(interpolateInitialConditions),cellcenter)",data) ;
        input("interpolateData(interpolateFile(interpolateInitialConditions),cellcenter)") ;
        input("Pambient,eos,initialCondition") ;
        output("priority::ic_file::rho_ic,priority::ic_file::v_ic") ;
        output("priority::ic_file::p_ic,priority::ic_file::T_ic") ;
        output("priority::ic_file::h_ic,priority::ic_file::y_ic") ;
        constraint("geom_cells,compressibleFlow") ;
        disable_threading();
      }

      // Assign cell values.
      void calculate(Entity cc) {
        int base_size = 7 ;
        real Tr = initialCondition->Temperature() ;
        real pr = data[cc][1] ;
        vect3d ur(data[cc][2],data[cc][3],data[cc][4]) ;
        int ns = eos->numSpecies() ;
        tmp_array<real> mix(ns) ;
        if(ns == 1) mix[0]=y_ic[cc][0]=1 ;
        else for(int i=0;i<ns;++i) mix[i]=y_ic[cc][i]=data[cc][base_size+i] ;
        EOS::State s = eos->State_from_mixture_p_T(mix,pr,Tr) ;
        v_ic[cc]=ur ; h_ic[cc]=s.enthalpy()+0.5*dot(ur,ur) ;
        rho_ic[cc] = s.density() ; p_ic[cc]=pr+(*Pambient) ; T_ic[cc]=Tr ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) {
        const int ns=eos->numSpecies() ; y_ic.setVecSize(ns) ;

        // Write message for user.
        if(Loci::MPI_rank==0){
          cout << "Restarting compressible flow from incompressible puT file."
            << endl ;
          cout << "Using temperature specified in 'initialCondition'." << endl ;
        }

        // Check that temperature was provided by the user.
        if(!initialCondition->IsTemperatureDefined()){
          cerr << "ERROR: Initial condition for temperature required!" << endl ;
          Loci::Abort() ;
        }

        do_loop(seq,this) ;
      }
  } ;

  register_rule<NewInterpolateCompressibleFromIncompressibleIC>
    registerNewInterpolateCompressibleFromIncompressibleIC ;

  /*class new_viscous_interpolate_ic : public pointwise_rule {
    const_param<conservativeVectorInfo> qvi ;
    const_param<EOS> eos ;
    const_param<reaction> reactor ;
    const_param<real> Pambient ;
    const_storeVec<real> data ;
    store<real> rho_ic ;
    storeVec<real> q_ic ;
    store<real> pg_ic ;
  public:
    new_viscous_interpolate_ic() ;
    void calculate(Entity cc) {
      int base_size = 7 ;
      real Tr = data[cc][0] ;
      real pr = data[cc][1] ;
      vect3d ur(data[cc][2],data[cc][3],data[cc][4]) ;
      int ns = eos->numSpecies() ;
      tmp_array<real> mix(ns) ;
      if(ns == 1)
        mix[0] = 1 ;
      else 
        for(int i=0;i<ns;++i)
          mix[i] = data[cc][base_size+i] ;
      EOS::State s = eos->State_from_mixture_p_T(mix,pr,Tr) ;
      for(int i=0;i<ns;++i)
        q_ic[cc][i] = mix[i]*s.density() ;
      const int mi = qvi->momentumIndex() ;
      q_ic[cc][mi+0] = ur.x*s.density() ;
      q_ic[cc][mi+1] = ur.y*s.density() ;
      q_ic[cc][mi+2] = ur.z*s.density() ;
      const int ei = qvi->totalEnergyIndex() ;
      q_ic[cc][ei] = s.rho_energy()+.5*dot(ur,ur)*s.density() ;
      rho_ic[cc] = s.density() ;
      pg_ic[cc] = pr ;
    }
    virtual void compute(const sequence &seq) {
      const int qs = qvi->vectorSize() ;

      q_ic.setVecSize(qs) ;
      do_loop(seq,this) ;
    }
  } ;

  new_viscous_interpolate_ic::new_viscous_interpolate_ic() {
    name_store("ic_file::q_ic",q_ic) ;
    name_store("ic_file::pg_ic",pg_ic) ;
    name_store("Pambient",Pambient) ;
    name_store("eos",eos) ;
    name_store("qvi",qvi) ;
    name_store("reactor",reactor) ;
    name_store("ic_file::rho_ic",rho_ic) ;
    name_store("interpolateDataCell(interpolateFile(interpolateInitialConditions),cellcenter)",data) ;
    input("interpolateDataCell(interpolateFile(interpolateInitialConditions),cellcenter)") ;
    input("Pambient") ;
    input("eos,qvi,reactor") ;
    output("ic_file::q_ic") ;
    output("ic_file::rho_ic") ;
    output("ic_file::pg_ic") ;
    constraint("geom_cells,TurbulentSimulation") ;
    disable_threading();
  }


  register_rule<new_viscous_interpolate_ic> register_new_viscous_interpolate_ic ;*/
  class new_turbulent_interpolate_ic : public pointwise_rule {
    const_storeVec<real> data ;
    store<real> k_ic,tmuu_ic ;
  public:
    new_turbulent_interpolate_ic() ;
    void calculate(Entity cc) {
      k_ic[cc] = data[cc][5] ;
      tmuu_ic[cc] = data[cc][6] ;
    }
    virtual void compute(const sequence &seq) {
      do_loop(seq,this) ;
    }
  } ;

  new_turbulent_interpolate_ic::new_turbulent_interpolate_ic() {
    name_store("ic_file::k_ic",k_ic) ;
    name_store("ic_file::tmuu_ic",tmuu_ic) ;
    name_store("interpolateDataCell(interpolateFile(interpolateInitialConditions),cellcenter)",data) ;
    input("interpolateDataCell(interpolateFile(interpolateInitialConditions),cellcenter)") ;
    output("ic_file::k_ic") ;
    output("ic_file::tmuu_ic") ;
    constraint("geom_cells,turbulentFlow") ;
    disable_threading();
  }


  register_rule<new_turbulent_interpolate_ic> register_new_turbulent_interpolate_ic ;

  // Rule to assign the initial condition for omega. Checked.
  class OmegaInitialConditionPUT : public pointwise_rule {
    private:
      const_store<real> rho_ic ;
      const_store<real> k_ic ;
      const_store<real> tmuu_ic ;
      store<real> omega_ic ;
    public:
                                                                                
      // Define input and output.
      OmegaInitialConditionPUT() {
//      name_store("ic_file::rho_ic",rho_ic) ;
//      name_store("ic_file::k_ic",k_ic) ;
//      name_store("ic_file::tmuu_ic",tmuu_ic) ;
//      name_store("ic_file::omega_ic",omega_ic) ;
//      input("ic_file::rho_ic,ic_file::k_ic,ic_file::tmuu_ic") ;
        name_store("rho_ic",rho_ic) ;
        name_store("k_ic",k_ic) ;
        name_store("tmuu_ic",tmuu_ic) ;
        name_store("ic_file::omega_ic",omega_ic) ;
        input("rho_ic,k_ic,tmuu_ic") ;
        output("ic_file::omega_ic") ;
        constraint("geom_cells") ;
      }
                                                                                
      // Assign omega for a single cell.
      void calculate(Entity cell) {
        omega_ic[cell]=rho_ic[cell]*max(k_ic[cell],1e-30)/max(tmuu_ic[cell],
          1e-300) ;
      }                                                                                
      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<OmegaInitialConditionPUT> registerOmegaInitialConditionPUT ;

  // Creates the turbulence model parameter so we can have compatibility with
  // CHEM on the put stuff. We do not use this anywhere else at the moment.
  class TurbulenceModel : public singleton_rule {
    private:
      const_param<string> flowRegime ;
      param<string> turbulenceModel ;
    public:
                                                                                
      // Define input and output.
      TurbulenceModel() {
        name_store("flowRegime",flowRegime) ;
        name_store("turbulence_model",turbulenceModel) ;
        input("flowRegime") ;
        output("turbulence_model") ;
        constraint("UNIVERSE") ;
      }
                                                                                
      // Set the value to Ed's for k-w since that is all we have.
      virtual void compute(const sequence &seq) {
        if(*flowRegime=="laminar") *turbulenceModel="none" ;
        else *turbulenceModel="SST" ;
      }
  } ;
                                                                                
  register_rule<TurbulenceModel> registerTurbulenceModel ;
}
