#line 1 "initialConditionRegions.loci"
// This is basically a hacked version of Ed's ic.loci file.

// Standard library includes.
#include <map>
#include <string>
using std::map ;
using std::string ;

#include "initialConditionRegions.h"
#line 1 "chem.lh"
#include <Loci>
#include <sciTypes_base.h>
#line 1 "FVM.lh"
//#############################################################################
//#
//# Copyright 2008, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################

// $type pos store<Loci::vector3d<Loci::real_t> > 
// $type cl Map
// $type cr Map
// $type ci Map
// $type ref Map
// $type pmap Map
// $type face2node multiMap

// $type upper multiMap
// $type lower multiMap
// $type boundary_map multiMap

// $type cellcenter store<Loci::vector3d<Loci::real_t> > 
// $type facecenter store<Loci::vector3d<Loci::real_t> > 
// $type area store<Loci::Area> 
// $type vol store<Loci::real_t> 
// $type grid_vol param<Loci::real_t> 

// $type mn store<Loci::vector3d<Loci::real_t> > 
// $type ln store<Loci::vector3d<Loci::real_t> > 

// $type grads(X0) store<Loci::vector3d<Loci::real_t> > 
// $type gradv(X0) storeVec<Loci::vector3d<Loci::real_t> > 
// $type gradv3d(X0) store<Loci::tensor3d<Loci::real_t> > 

// $type grads_f(X0) store<Loci::vector3d<Loci::real_t> > 
// $type gradv_f(X0) storeVec<Loci::vector3d<Loci::real_t> > 
// $type gradv3d_f(X0) store<Loci::tensor3d<Loci::real_t> > 

// $type limiters(X0) store<Loci::real_t> 
// $type limiterv(X0) storeVec<Loci::real_t> 
// $type limiterv3d(X0) store<Loci::vector3d<Loci::real_t> > 

// $type lefts(X0) store<Loci::real_t> 
// $type rights(X0) store<Loci::real_t> 
// $type leftsP(X0,X1) store<Loci::real_t> 
// $type rightsP(X0,X1) store<Loci::real_t> 
// $type leftvM(X0) storeVec<Loci::real_t> 
// $type rightvM(X0) storeVec<Loci::real_t> 
// $type leftv3d(X0) store<Loci::vector3d<Loci::real_t> > 
// $type rightv3d(X0) store<Loci::vector3d<Loci::real_t> > 

// $type cell2node(X0) store<float> 
// $type cell2node_v(X0) storeVec<float> 
// $type cell2node_v3d(X0) store<Loci::vector3d<float> > 
// $type cell2nodeMax(X0) store<float> 
// $type cell2nodeMin(X0) store<float> 
// $type cell2nodeMaxMag(X0) store<float> 
// $type cell2nodeMaxv3d(X0) store<vector3d<float> > 

// $type BC_options store<Loci::options_list> 

// $type integrateSurface(X0) store<Loci::real_t> 
// $type integrateFlux(X0) store<Loci::real_t> 

// $type petscScalarSolve(X0) store<Loci::real_t> 
// $type petscBlockedSolve(X0) storeVec<Loci::real_t> 
// $type petscBlockedSSolve(X0) storeVec<Loci::real_t> 

// $type L2Norm(X0) param<Loci::real_t> 
// $type L1Norm(X0) param<Loci::real_t> 
// $type LinfNorm(X0) param<Loci::real_t> 
#line 3 "chem.lh"


// $type modelName param<string> 

// Transform for periodic bc's (should this be in FVM.lh?)
// $type periodicTransform store<rigid_transform> 

// EOS related variables
// Ambient Pressure used to define gauge pressure
// $type Pambient param<streamUns::real> 
// The equation of state object
// $type eos param<fluidPhysics::EOS> 
// An object that represents the thermodynamic state of a cell
// $type eos_state store<fluidPhysics::EOS::State> 
// An object that represents the thermodynamic state of a mixture at a cell
// $type eos_mixture_state storeVec<streamUns::real> 
// A hint for the equation of state to speed convergence
// $type hint storeVec<float> 
// $type hint_n storeVec<float> 
// Energy of each species 
// $type species_energy storeVec<streamUns::real> 
// cv of each species
// $type species_cv storeVec<streamUns::real> 
// Derivative of internal energy with respect to species density (holding P constant)
// $type dreidri storeVec<streamUns::real> 
// Derivative of internal energy with respect to pressure (holding density constant)
// $type dreidP store<streamUns::real> 
// cell temperature
// $type temperature store<streamUns::real> 
// cell density
// $type rho store<streamUns::real> 
// cell mixture
// $type mixture storeVec<streamUns::real> 
// cell pressure
// $type pressure store<streamUns::real> 
// cell speed of sound
// $type soundSpeed store<streamUns::real> 
// cell cp (and face cp)
// $type cp store<streamUns::real> 
// internal energy for cells
// $type e_internal store<streamUns::real> 
// derivative of temperature with respect to primitive variable
// $type dTdq storeVec<streamUns::real> 
// enthalpy of each species ata cell
// $type speciesEnthalpy storeVec<streamUns::real> 
// molecular mass of each species
// $type speciesMolecularMass param<std::vector<streamUns::real> > 
// Mixture fraction to use when non given
// $type defaultMixtureFraction param<std::vector<streamUns::real> > 

// Other fluid variables
//
// $type gaugePressure store<streamUns::real> 
// fluid velocity at a cell
// $type u store<streamUns::vect3d> 
// grid velocity at face
// $type us_n store<streamUns::real> 
// $type us store<streamUns::vect3d> 

// Face values
// $type gaugePressure_f store<streamUns::real> 
// $type mixture_f storeVec<streamUns::real> 
// $type temperature_f store<streamUns::real> 
// $type rho_f store<streamUns::real> 
// $type soundSpeed_f store<streamUns::real> 
// $type u_f store<streamUns::vect3d> 

// $type boundaryValue(X0) store<streamUns::real> 	
// $type boundaryValueP(X0,X1) store<streamUns::real> 	
// $type boundaryValuevM(X0) storeVec<streamUns::real> 
// $type boundaryValueV3D(X0) store<streamUns::vect3d>  
// $type boundaryState store<EOS::State> 	
// $type boundary_area store<streamUns::real> 
// $type bc_total_area store<streamUns::real> 
// $type bc_total_force store<streamUns::real> 
// $type bc_average_pressure store<streamUns::real> 

// $type flow_direction store<streamUns::vect3d> 
// $type rigid_u store<streamUns::vect3d> 
// $type p0Ref store<streamUns::real> 
// $type T0Ref store<streamUns::real> 
// $type mixtureRef storeVec<streamUns::real> 
// $type massFluxRef store<streamUns::real> 
// $type temperatureRef store<streamUns::real> 
// $type gaugePressureRef store<streamUns::real> 
// $type uRef store<streamUns::vect3d> 

// Epsilon used by numerical jacobians
// $type epsilon param<streamUns::real> 

// Mean Flow Solver Variables
// $type numSpecies param<int> 
// Information on the layout of the conservative and primitive vectors
// $type qvi param<conservativeVectorInfo> 
// primitive solution vector
// $type qp storeVec<streamUns::real> 
// conservative solution vector
// $type q storeVec<streamUns::real> 
// change in primitive variables in newton step
// $type dq storeVec<streamUns::real> 
// All non-temporal terms including source terms, inviscid, and viscous fluxes.
// $type src storeVec<streamUns::real> 
// time dependent residual (right hand side of linear system equations)
// $type rhs storeVec<streamUns::real> 
// inviscid flux
// $type qf storeVec<streamUns::real>  // inviscid flux
// coefficients defining inviscid flux for scalar transport
// $type scalar_l store<streamUns::real> 
// $type scalar_r store<streamUns::real> 
// Preconditioning functions
// $type preconditioning param<streamUns::preconditioning_method> 
// $type pc_param param<streamUns::preconditioning_param> 	
// $type Beta store<streamUns::real> 

// Initial condition support
// $type q_ic storeVec<streamUns::real> 
// $type pg_ic store<streamUns::real> 
// $type pg store<streamUns::real> 

// Time integration parameters
// $type theta param<streamUns::real>  
// $type psi param<streamUns::real>  

// Viscous flux
// $type diff storeVec<streamUns::real> 
// $type diff_f storeVec<streamUns::real> 

// $type vis_flux storeVec<streamUns::real>  // viscous flux
// Fluid Stress
// $type tau store<streamUns::symmetricTensor> 

// Chemistry stuff
// Freezing temperature for chemistry
// $type Tf param<streamUns::TemperatureValue> 
// Object that computes reaction rates
// $type reactor param<fluidPhysics::reaction> 
// reaction rates and their derivatives
// $type rates storeVec<fluidPhysics::reaction::rates> 

// Mean Flow Jacobian
// Diagnal jacobian term
// $type pc_srcJ storeMat<streamUns::real_fj>  
// flux jacobians from left side
// $type pc_fjp storeMat<streamUns::real_fj> 
// flux jacobians from right side
// $type pc_fjm storeMat<streamUns::real_fj> 
// fluid linear system
// $type fluid_B storeVec<streamUns::real_fj> 
// $type fluid_D storeMat<streamUns::real_fj> 
// $type fluid_L storeMat<streamUns::real_fj> 
// $type fluid_U storeMat<streamUns::real_fj> 
// $type fluid_D_inv storeMat<streamUns::real_fj> 
// $type fluid_pivot storeVec<pivot_type> 

// Boundary Condition Parameters
// $type Uwall_BC store<streamUns::vect3d> 
// $type T_BC store<streamUns::real> 
// $type T0_BC store<streamUns::real> 
// $type p0_BC store<streamUns::real> 
// $type rho_BC store<streamUns::real> 
// $type Twall_BC store<streamUns::real> 
// $type mixture_BC storeVec<streamUns::real> 
// $type massFlux_BC store<streamUns::real> 
// $type mdot_BC store<streamUns::real> 
// $type swirlAngle_BC store<streamUns::real> 
// $type swirlCenter_BC store<streamUns::vect3d> 
// $type momentCenter_BC store<streamUns::vect3d> 
// $type M_BC store<streamUns::vect3d> 
// $type u_BC store<streamUns::vect3d> 
// $type p_BC store<streamUns::real> 
// $type BC_inflow_type store<int> 
// $type temperatureRef_BC store<streamUns::real> 
// $type gaugePressureRef_BC store<streamUns::real> 
// $type mixtureRef_BC storeVec<streamUns::real> 
// $type uRef_BC store<streamUns::vect3d> 

// $type angVel_BC store<streamUns::vect3d> 
// $type rotAxis_BC store<streamUns::vect3d> 
// $type rotCenter_BC store<streamUns::vect3d> 
// $type swirlAxis_BC store<streamUns::vect3d> 
// $type flowDir_BC store<streamUns::vect3d> 

// $type Twall store<streamUns::real> 
// $type Twall_prescribed store<streamUns::real> 
// $type T_interface store<streamUns::real> 
// $type wallVelocity store<streamUns::vect3d> 

// $type sst_k_bc store<streamUns::real> 
// $type sst_w_bc store<streamUns::real> 
// $type prescribe_turb_ref store<vec<2> > 
// $type lam_tau_w store<streamUns::real> 

// Wall law parameters
// $type wall_law param<streamUns::wall_func_param> 
// $type tau_wall store<streamUns::vect3d> 
// $type q_wall store<streamUns::real> 
// $type wall_law_k store<streamUns::real> 
// $type wall_law_w store<streamUns::real> 
// $type temp_wlaw store<streamUns::real> 
// $type mixture_f_wlaw storeVec<streamUns::real> 
// $type nonAdiabaticWall store<bool> 
// $type press_wlaw store<streamUns::real> 
// $type wall_law_nut store<streamUns::real> 
// $type wall_cells store<bool> 

// $type hflux store<streamUns::real> 



// $type stime param<streamUns::real>  // simulation time
// $type ncycle param<int>  // simulation iteration number

// $type newton_iter param<int>   // Timestepping controls
// $type urelax param<streamUns::real> 
// $type dtmax param<streamUns::TimeValue> 
// $type cflmax param<streamUns::real> 
// $type cfl param<std::pair<double,double> > 

// $type max_ev store<streamUns::real> 
// preconditioning timestep
// $type dtau store<streamUns::real> 
// physical timestep
// $type dt store<streamUns::real> 
// non-linear limited timestep
// $type dt_urelax store<streamUns::real> 
// local timestep
// $type dtmin store<streamUns::real> 
// delta t/volume
// $type dtvol store<streamUns::real> 

// Generic transport
// $type muu(X0,X1,X2) store<streamUns::real> 
// $type kconduct(X0,X1,X2) store<streamUns::real> 

// Specific transport for cells and faces
// viscosity
// $type muu store<streamUns::real> 
// $type muu_f store<streamUns::real> 

// conductivity
// $type kconduct store<streamUns::real> 
// $type kconduct_f store<streamUns::real> 

// General Turbulence Model
// Map from cell to nearest viscous wall
// $type min_cell2noslip Map

// $type divu store<streamUns::real> 
// vorticity
// $type vortMag store<streamUns::real> 
// $type vort store<streamUns::vect3d> 
// $type vort_f store<streamUns::vect3d> 

// helicity
// $type helicity store<streamUns::real> 
// $type helicity_f store<streamUns::real> 

// distance to nearest viscous wall
// $type dist_noslip store<streamUns::real> 
// $type dist_noslip_f store<streamUns::real> 
// turbulent viscousity
// $type tmuu store<streamUns::real> 
// $type tmuu_f store<streamUns::real> 

// Spalart Allmaras Variables
// $type spalart param<Spa_All_param> 

// conservative integrated variable
// $type rho_nu_t store<streamUns::real> 
// turbulent kinematic viscosity
// $type nu_t store<streamUns::real> 
// $type nu_t_f store<streamUns::real> 
// Chi coefficeint for SA model
// $type xcoeff store<streamUns::real> 
// SA model coefficeitns
// $type fv1 store<streamUns::real> 
// $type fv1_f store<streamUns::real> 
// $type fv2 store<streamUns::real> 
// $type S_tilda store<streamUns::real> 
// $type rcoeff store<streamUns::real> 
// $type rpcoeff store<streamUns::real> 
// $type gcoeff store<streamUns::real> 
// $type gpcoeff store<streamUns::real> 
// $type fw store<streamUns::real> 
// $type fwp store<streamUns::real> 
// $type cw1 store<streamUns::real> 
// convection, diffusion, and source terms
// $type SA_src store<streamUns::real> 
// turmbulent destruction
// $type SA_destruct store<streamUns::real> 
// flux jacobians
// $type SA_fjm store<streamUns::real> 
// $type SA_fjp store<streamUns::real> 
// diagonal jacobian term
// $type SA_srcJ store<streamUns::real> 
// right hand side of newton method
// $type SA_rhs store<streamUns::real> 
// SA matrix
// $type SA_B store<streamUns::real> 
// $type SA_dq store<streamUns::real> 
// $type SA_L store<streamUns::real> 
// $type SA_D store<streamUns::real> 
// $type SA_U store<streamUns::real> 
// $type SA_qi store<streamUns::real> 

// Realizable k-epsilon model
// $type rke_q store<streamUns::vec<2> > 
// $type rke_src store<streamUns::vec<2> > 
// $type rke_rhs store<streamUns::vec<2> > 
// diagonal of jacobian
// $type rke_srcJ storeMat<streamUns::real_fj> 
// flux jacobian terms
// $type rke_fjp storeMat<streamUns::real_fj> 
// $type rke_fjm storeMat<streamUns::real_fj> 
// system matrix for RKE model
// $type rke_B storeVec<streamUns::real_fj> 
// $type rke_D storeMat<streamUns::real_fj> 
// $type rke_L storeMat<streamUns::real_fj> 
// $type rke_U storeMat<streamUns::real_fj> 
// $type e store<streamUns::real> 
// $type e_f store<streamUns::real> 
// $type diff_e store<streamUns::real> 


// K-omega turbulence models
// Conservative variables (rk, rw)
// $type sst_q store<streamUns::vec<2> > 
// $type sst_qi store<streamUns::vec<2> > 
// $type sst_dq storeVec<streamUns::real> 
// turbulent kinetic energy
// $type k store<streamUns::real> 
// $type k_f store<streamUns::real> 
// omega (enstrophy)
// $type w store<streamUns::real> 
// $type w_f store<streamUns::real> 
// SST blending functions
// $type bF1 store<streamUns::real> 
// $type bF2 store<streamUns::real> 
// $type fBetaS store<streamUns::real> 
// $type fBeta store<streamUns::real> 
// Sigmas
// $type sigma_k store<streamUns::real> 
// $type sigmak_f store<streamUns::real> 
// $type sigmae_f store<streamUns::real> 
// SST turbulence model parameters
// $type sst store<streamUns::sst_param> 
// $type sst1 param<sst1_param> 
// $type sst2 param<sst2_param> 
// tke diffusion flux
// $type diff_t store<streamUns::real> 
// omega diffusion flux
// $type diff_w store<streamUns::real> 
//strain
// $type strainRate store<streamUns::real> 
// $type stress_scalar store<streamUns::real> 
// turbulent production
// $type P_k store<streamUns::real> 

// Source terms for k-omega model
// $type sst_src store<streamUns::vec<2> > 
// right hand side of newton, includes time derivative terms
// $type sst_rhs store<streamUns::vec<2> > 

// diagonal of jacobian
// $type sst_srcJ storeMat<streamUns::real_fj> 
// flux jacobian terms
// $type sst_fjp storeMat<streamUns::real_fj> 
// $type sst_fjm storeMat<streamUns::real_fj> 
// system matrix for SST model
// $type sst_B storeVec<streamUns::real_fj> 
// $type sst_D storeMat<streamUns::real_fj> 
// $type sst_L storeMat<streamUns::real_fj> 
// $type sst_U storeMat<streamUns::real_fj> 


// $type plot_modulo param<int> 

// $type restart_freq param<int> 
// $type restart_modulo param<int> 
// $type restart_directory param<string> 
// $type restart_postfix param<string> 

// $type stop_iter param<int> 
// $type plot_freq param<int> 
// $type print_freq param<int> 
// $type time_integration param<std::string> 
// $type pseudo_cflmax param<double> 
// $type gauss_seidel_iter param<int> 

// $type heatf store<streamUns::real> 
// $type wall_stress store<streamUns::vect3d> 
// $type wall_p store<streamUns::real> 
// $type wall_T store<streamUns::real> 
// $type yplus_w store<streamUns::real> 
#line 10 "initialConditionRegions.loci"


namespace streamUns {

  // $type initialConditionRegions param<Loci::options_list> 
  namespace {class file_initialConditionRegions000_1278520192m350 : public Loci::optional_rule {
#line 15 "initialConditionRegions.loci"
    Loci::param<Loci::options_list>  L_initialConditionRegions_ ; 
#line 15 "initialConditionRegions.loci"
public:
#line 15 "initialConditionRegions.loci"
    file_initialConditionRegions000_1278520192m350() {
#line 15 "initialConditionRegions.loci"
       name_store("initialConditionRegions",L_initialConditionRegions_) ;
#line 15 "initialConditionRegions.loci"
       output("initialConditionRegions") ;
#line 15 "initialConditionRegions.loci"
    }
#line 15 "initialConditionRegions.loci"
    void compute(const Loci::sequence &seq) { }} ;
#line 15 "initialConditionRegions.loci"
Loci::register_rule<file_initialConditionRegions000_1278520192m350> register_file_initialConditionRegions000_1278520192m350 ;
#line 15 "initialConditionRegions.loci"
}
#line 15 "initialConditionRegions.loci"


  class leftPlane : public geomTest {
    vect3d pt ;
    vect3d n ;
  public:
    leftPlane(vect3d p1, vect3d n1) : pt(p1),n(n1) {
      n *= 1./norm(n) ;
    }
    leftPlane(const options_list &ol) {
      if(!ol.optionExists("point") || !ol.optionExists("normal")) {
        cerr << "leftPlane needs 'point' and 'normal'" << endl ;
        Loci::Abort() ;
      }
      get_vect3dOption(ol,"point","m",pt) ;
      get_vect3dOption(ol,"normal","",n) ;
      n *= 1./max(norm(n),1e-33) ;
    }
    bool inGeomPt(vect3d pt) const ;
  } ;

  bool leftPlane::inGeomPt(vect3d p) const {
    if(dot(p-pt,n) < 0.0)
      return true ;
    return false ;
  }

  class inBox : public geomTest {
    double xmax,xmin,ymax,ymin,zmax,zmin ;
  public:
    inBox(vect3d p1, vect3d p2) {
      xmax = max(p1.x,p2.x) ;
      xmin = min(p1.x,p2.x) ;

      ymax = max(p1.y,p2.y) ;
      ymin = min(p1.y,p2.y) ;

      zmax = max(p1.z,p2.z) ;
      zmin = min(p1.z,p2.z) ;
    }
    inBox(const options_list &ol) {
      if(!ol.optionExists("p1") || !ol.optionExists("p2")) {
        cerr << "inBox needs two points, 'p1' and 'p2'" << endl ;
        Loci::Abort() ;
      }
      vect3d p1,p2 ;
      get_vect3dOption(ol,"p1","m",p1) ;
      get_vect3dOption(ol,"p2","m",p2) ;
      xmax = max(p1.x,p2.x) ;
      xmin = min(p1.x,p2.x) ;

      ymax = max(p1.y,p2.y) ;
      ymin = min(p1.y,p2.y) ;

      zmax = max(p1.z,p2.z) ;
      zmin = min(p1.z,p2.z) ;
    }

    bool inGeomPt(vect3d pt) const ;
  } ;

  bool inBox::inGeomPt(vect3d p) const {
    if(p.x >= xmin && p.x <= xmax &&
       p.y >= ymin && p.y <= ymax &&
       p.z >= zmin && p.z <= zmax)
      return true ;
    return false ;
  }

  class inSphere : public geomTest {
    double r ;
    double r2 ;
    vect3d center ;
  public:
    inSphere(double ri,  vect3d c) { r=ri; center=c; r2=r*r ; }
    inSphere(const options_list &ol) {
      // Here get sphere information
      if(!ol.optionExists("radius") || !ol.optionExists("center")) {
        cerr << "inSphere needs a radius and center" << endl ;
        Loci::Abort() ;
      }
      ol.getOptionUnits("radius","m",r) ;
      get_vect3dOption(ol,"center","m",center) ;
      r2=r*r ;
    }
    bool inGeomPt(vect3d pt) const ;
  } ;

  bool inSphere::inGeomPt(vect3d pt) const {
    if(dot(pt-center,pt-center) < r2)
      return true ;
    return false ;
  }

  class inCylinder : public geomTest {
    double r ;
    double r2 ;
    vect3d p1,p2 ;
    vect3d n ;
    double mag ;
  public:
    inCylinder(double ri,  vect3d p1i,vect3d p2i) {
      r=ri;r2=r*r ;
      p1=p1i;p2=p2i;
      n = p2-p1 ;
      mag = norm(n) ;
      n *= 1./mag ;
    }
    inCylinder(const options_list &ol) {
      // Here get sphere information
      if(!ol.optionExists("radius") || !ol.optionExists("p1")
         || !ol.optionExists("p2")) {
        cerr << "inCylinder needs a radius and two axis points (p1 and p2)" << endl ;
        Loci::Abort() ;
      }
      ol.getOptionUnits("radius","m",r) ;
      get_vect3dOption(ol,"p1","m",p1) ;
      get_vect3dOption(ol,"p2","m",p2) ;
      n = p2-p1 ;
      mag = norm(n) ;
      n *= 1./mag ;
      r2=r*r ;
    }
    bool inGeomPt(vect3d pt) const ;
  } ;

  bool inCylinder::inGeomPt(vect3d pt) const {
    double v = dot(pt-p1,n) ;
    if(v < 0.0 || v > mag)
      return false ;
    // project point onto axis
    vect3d paxis = p1 + n*v ;
    if(dot(pt-paxis,pt-paxis) > r2) // Check radius
      return false ; // outside
    return true ;// otherwise inside cylinder
  }

  class inCone : public geomTest {
    double r1 ;
    double r2 ;
    vect3d p1,p2 ;
    vect3d n ;
    double mag ;
  public:
    inCone(double r1i, double r2i, vect3d p1i,vect3d p2i) {
      r1=r1i;r2=r2i ;
      p1=p1i;p2=p2i;
      n = p2-p1 ;
      mag = norm(n) ;
      n *= 1./mag ;
    }
    inCone(const options_list &ol) {
      // Here get sphere information
      if(!ol.optionExists("r1") || !ol.optionExists("r2") ||
         !ol.optionExists("p1") || !ol.optionExists("p2")) {
        cerr << "inCone needs a two radii (r1 and r2) and two axis points (p1 and p2)" << endl ;
        Loci::Abort() ;
      }
      ol.getOptionUnits("r1","m",r1) ;
      ol.getOptionUnits("r2","m",r2) ;
      get_vect3dOption(ol,"p1","m",p1) ;
      get_vect3dOption(ol,"p2","m",p2) ;
      n = p2-p1 ;
      mag = norm(n) ;
      n *= 1./mag ;
    }
    bool inGeomPt(vect3d pt) const ;
  } ;

  bool inCone::inGeomPt(vect3d pt) const {
    double v = dot(pt-p1,n) ;
    if(v < 0.0 || v > mag)
      return false ;
    // project point onto axis
    double t = v/mag ;
    double rx = r1*(1.-t) + r2*t ;
    double rx2 = rx*rx ;
    vect3d paxis = p1 + n*v ;
    if(dot(pt-paxis,pt-paxis) > rx2) // Check radius
      return false ; // outside
    return true ;// otherwise inside cone
  }

  Loci::CPTR<geomTest> geomTestFactory(string name, const options_list ol) {
    Loci::CPTR<geomTest> gp ;
    if(name == "inSphere") {
      gp = new inSphere(ol) ;
    } else if(name == "inBox") {
      gp = new inBox(ol) ;
    } else if(name == "inCylinder") {
      gp = new inCylinder(ol) ;
    } else if(name == "inCone") {
      gp = new inCone(ol) ;
    } else if(name == "leftPlane") {
      gp = new leftPlane(ol) ;
    } else {
      cerr << "don't know what to do with '" << name << "'" << endl ;
      return 0 ;
    }
    return gp ;
  }

  // $type icRegionInfo blackbox<ICparsedInitRegion> 

  fluidState getStateFromList(const options_list &ol, string name) {
    fluidState f ;
    using namespace Loci ;
    option_value_type ovt= ol.getOptionValueType(name) ;
    if(ovt != FUNCTION) {
      cerr << "getState failed for variable '" << name
           << "' expecting state function" << endl ;
      Loci::Abort() ;
    }
    options_list::arg_list value_list ;
    string fname ;
    ol.getOption(name,fname,value_list) ;
    if(fname != "state") {
      cerr << "expecting to define function 'state', instead found '"
           << fname << "'" << endl ;
      Loci::Abort() ;
    }
    f.Input(value_list) ;
    return f ;
  }

  namespace {class file_initialConditionRegions001_1278520192m353 : public Loci::singleton_rule {
#line 240 "initialConditionRegions.loci"
    Loci::const_param<streamUns::real>  L_Pambient_ ; 
#line 240 "initialConditionRegions.loci"
    Loci::const_param<fluidPhysics::EOS>  L_eos_ ; 
#line 240 "initialConditionRegions.loci"
    Loci::const_param<conservativeVectorInfo>  L_qvi_ ; 
#line 240 "initialConditionRegions.loci"
    Loci::const_param<fluidPhysics::reaction>  L_reactor_ ; 
#line 240 "initialConditionRegions.loci"
    Loci::const_param<Loci::options_list>  L_initialConditionRegions_ ; 
#line 240 "initialConditionRegions.loci"
    Loci::blackbox<ICparsedInitRegion>  L_icRegionInfo_ ; 
#line 240 "initialConditionRegions.loci"
public:
#line 240 "initialConditionRegions.loci"
    file_initialConditionRegions001_1278520192m353() {
#line 240 "initialConditionRegions.loci"
       name_store("Pambient",L_Pambient_) ;
#line 240 "initialConditionRegions.loci"
       name_store("eos",L_eos_) ;
#line 240 "initialConditionRegions.loci"
       name_store("qvi",L_qvi_) ;
#line 240 "initialConditionRegions.loci"
       name_store("reactor",L_reactor_) ;
#line 240 "initialConditionRegions.loci"
       name_store("initialConditionRegions",L_initialConditionRegions_) ;
#line 240 "initialConditionRegions.loci"
       name_store("icRegionInfo",L_icRegionInfo_) ;
#line 240 "initialConditionRegions.loci"
       input("initialConditionRegions,qvi,eos,Pambient,reactor") ;
#line 240 "initialConditionRegions.loci"
       output("icRegionInfo") ;
#line 240 "initialConditionRegions.loci"
       disable_threading() ;
#line 240 "initialConditionRegions.loci"
    }
#line 240 "initialConditionRegions.loci"
    void compute(const Loci::sequence &seq) { 
    if(!(*L_initialConditionRegions_).optionExists("default")) {
      cerr << "default initial condition must be set in initialConditionRegions" << endl ;
      Loci::Abort() ;
    }
    int vs = (*L_qvi_).vectorSize() ;

    (*L_icRegionInfo_).defaultState = getStateFromList((*L_initialConditionRegions_),
                                                  "default") ;
    if((*L_initialConditionRegions_).optionExists("regions")) {
      using namespace Loci ;
      option_value_type ovt= (*L_initialConditionRegions_).getOptionValueType("regions") ;
      if(ovt != LIST) {
        cerr << "regions in initialConditionRegions should define a list" << endl ;
        Loci::Abort() ;
      }
      options_list::arg_list value_list ;
      (*L_initialConditionRegions_).getOption("regions",value_list) ;
      int sz = value_list.size() ;
      for(int i=0;i<sz;++i) {
        if(value_list[i].type_of() != FUNCTION) {
          cerr << "regions is a list of geometric test functions" << endl ;
          Loci::Abort() ;
        }
        string name ;
        value_list[i].get_value(name) ;
        options_list::arg_list fvalues ;
        value_list[i].get_value(fvalues) ;
        options_list fol ;
        fol.Input(fvalues) ;
        ICstate_info sinfo ;
        string composition = "default" ;
        if(fol.optionExists("composition")) {
          option_value_type ovt= fol.getOptionValueType("composition") ;
          if(ovt != NAME) {
            cerr << "composition should be set to the name of a specified state." << endl ;
          } else {
            fol.getOption("composition",composition) ;
          }
        }

        sinfo.regionState = getStateFromList((*L_initialConditionRegions_),
                                             composition) ;
        { vector<double> tmp(vs) ;
          sinfo.q.swap(tmp) ;
        }
        { vector<double> tmp(vs) ;
          sinfo.qp.swap(tmp) ;
        }
        sinfo.regionState.setState(&sinfo.q[0],
                                   (*L_qvi_),(*L_eos_),(*L_reactor_)) ;
        sinfo.regionState.setPrimitive(&sinfo.qp[0],
                                       (*L_Pambient_),(*L_qvi_),(*L_eos_),(*L_reactor_)) ;
//      print_state($eos,$qvi,&sinfo.qp[0],$Pambient,composition) ;
        sinfo.name = composition ;
        sinfo.geomTestFunc = geomTestFactory(name,fol) ;

        (*L_icRegionInfo_).fluidRegions.push_back(sinfo) ;

      }
    } else {
      if(Loci::MPI_rank == 0)
        cerr << "Warning: No regions defined in initialConditionRegions!" <<endl ;
    }
    { vector<double> tmp(vs) ;
      (*L_icRegionInfo_).default_q.swap(tmp) ;
    }
    { vector<double> tmp(vs) ;
      (*L_icRegionInfo_).default_qp.swap(tmp) ;
    }
    (*L_icRegionInfo_).defaultState.setState(&(*L_icRegionInfo_).default_q[0],
                                        (*L_qvi_),(*L_eos_),(*L_reactor_)) ;
    (*L_icRegionInfo_).defaultState.setPrimitive(&(*L_icRegionInfo_).default_qp[0],
                                            (*L_Pambient_),(*L_qvi_),(*L_eos_),(*L_reactor_)) ;
//  print_state($eos,$qvi,&$icRegionInfo.default_qp[0],$Pambient,"default") ;

  }} ;
#line 316 "initialConditionRegions.loci"
Loci::register_rule<file_initialConditionRegions001_1278520192m353> register_file_initialConditionRegions001_1278520192m353 ;
#line 316 "initialConditionRegions.loci"
}
#line 316 "initialConditionRegions.loci"


  // Initial condition for velocity.
  // $type v_ic store<vect3d> 
  namespace {class file_initialConditionRegions002_1278520192m355 : public Loci::pointwise_rule {
#line 320 "initialConditionRegions.loci"
    Loci::const_store<Loci::vector3d<Loci::real_t> >  L_cellcenter_ ; 
#line 320 "initialConditionRegions.loci"
    Loci::const_blackbox<ICparsedInitRegion>  L_icRegionInfo_ ; 
#line 320 "initialConditionRegions.loci"
    Loci::store<vect3d>  L_v_ic_ ; 
#line 320 "initialConditionRegions.loci"
public:
#line 320 "initialConditionRegions.loci"
    file_initialConditionRegions002_1278520192m355() {
#line 320 "initialConditionRegions.loci"
       name_store("cellcenter",L_cellcenter_) ;
#line 320 "initialConditionRegions.loci"
       name_store("icRegionInfo",L_icRegionInfo_) ;
#line 320 "initialConditionRegions.loci"
       name_store("v_ic",L_v_ic_) ;
#line 320 "initialConditionRegions.loci"
       input("cellcenter,icRegionInfo") ;
#line 320 "initialConditionRegions.loci"
       output("v_ic") ;
#line 320 "initialConditionRegions.loci"
       constraint("geom_cells") ;
#line 320 "initialConditionRegions.loci"
    }
#line 320 "initialConditionRegions.loci"
    void calculate(Entity _e_) { 
#line 321 "initialConditionRegions.loci"

    // Set the default state.
    L_v_ic_[_e_]=L_icRegionInfo_[_e_].defaultState .getVelocity () ;

    // Loop over regions, overwriting velocity with latest region value if the
    // cell center is in the region.
    for (size_t j =0;j <L_icRegionInfo_[_e_].fluidRegions .size ();++j ) {
      if (L_icRegionInfo_[_e_].fluidRegions [j ].geomTestFunc ->inGeomPt (L_cellcenter_[_e_])) {
        L_v_ic_[_e_]=L_icRegionInfo_[_e_].fluidRegions [j ].regionState .getVelocity () ;
      }
    }
  }    void compute(const Loci::sequence &seq) { 
#line 332 "initialConditionRegions.loci"
      do_loop(seq,this) ;
#line 332 "initialConditionRegions.loci"
    }
#line 332 "initialConditionRegions.loci"
} ;
#line 332 "initialConditionRegions.loci"
Loci::register_rule<file_initialConditionRegions002_1278520192m355> register_file_initialConditionRegions002_1278520192m355 ;
#line 332 "initialConditionRegions.loci"
}
#line 332 "initialConditionRegions.loci"


  // Initial condition for pressure.
  // $type p_ic store<real> 
  namespace {class file_initialConditionRegions003_1278520192m356 : public Loci::pointwise_rule {
#line 336 "initialConditionRegions.loci"
    Loci::const_store<Loci::vector3d<Loci::real_t> >  L_cellcenter_ ; 
#line 336 "initialConditionRegions.loci"
    Loci::const_blackbox<ICparsedInitRegion>  L_icRegionInfo_ ; 
#line 336 "initialConditionRegions.loci"
    Loci::store<real>  L_p_ic_ ; 
#line 336 "initialConditionRegions.loci"
public:
#line 336 "initialConditionRegions.loci"
    file_initialConditionRegions003_1278520192m356() {
#line 336 "initialConditionRegions.loci"
       name_store("cellcenter",L_cellcenter_) ;
#line 336 "initialConditionRegions.loci"
       name_store("icRegionInfo",L_icRegionInfo_) ;
#line 336 "initialConditionRegions.loci"
       name_store("p_ic",L_p_ic_) ;
#line 336 "initialConditionRegions.loci"
       input("cellcenter,icRegionInfo") ;
#line 336 "initialConditionRegions.loci"
       output("p_ic") ;
#line 336 "initialConditionRegions.loci"
       constraint("geom_cells") ;
#line 336 "initialConditionRegions.loci"
    }
#line 336 "initialConditionRegions.loci"
    void calculate(Entity _e_) { 
#line 337 "initialConditionRegions.loci"

    // Set the default state.
    L_p_ic_[_e_]=L_icRegionInfo_[_e_].defaultState .getPressure () ;

    // Loop over regions, overwriting pressure with latest region value if the
    // cell center is in the region.
    for (size_t j =0;j <L_icRegionInfo_[_e_].fluidRegions .size ();++j ) {
      if (L_icRegionInfo_[_e_].fluidRegions [j ].geomTestFunc ->inGeomPt (L_cellcenter_[_e_])) {
        L_p_ic_[_e_]=L_icRegionInfo_[_e_].fluidRegions [j ].regionState .getPressure () ;
      }
    }
  }    void compute(const Loci::sequence &seq) { 
#line 348 "initialConditionRegions.loci"
      do_loop(seq,this) ;
#line 348 "initialConditionRegions.loci"
    }
#line 348 "initialConditionRegions.loci"
} ;
#line 348 "initialConditionRegions.loci"
Loci::register_rule<file_initialConditionRegions003_1278520192m356> register_file_initialConditionRegions003_1278520192m356 ;
#line 348 "initialConditionRegions.loci"
}
#line 348 "initialConditionRegions.loci"


  // Initial condition for temperature.
  // $type T_ic store<real> 
  namespace {class file_initialConditionRegions004_1278520192m356 : public Loci::pointwise_rule {
#line 352 "initialConditionRegions.loci"
    Loci::const_store<Loci::vector3d<Loci::real_t> >  L_cellcenter_ ; 
#line 352 "initialConditionRegions.loci"
    Loci::const_blackbox<ICparsedInitRegion>  L_icRegionInfo_ ; 
#line 352 "initialConditionRegions.loci"
    Loci::store<real>  L_T_ic_ ; 
#line 352 "initialConditionRegions.loci"
public:
#line 352 "initialConditionRegions.loci"
    file_initialConditionRegions004_1278520192m356() {
#line 352 "initialConditionRegions.loci"
       name_store("cellcenter",L_cellcenter_) ;
#line 352 "initialConditionRegions.loci"
       name_store("icRegionInfo",L_icRegionInfo_) ;
#line 352 "initialConditionRegions.loci"
       name_store("T_ic",L_T_ic_) ;
#line 352 "initialConditionRegions.loci"
       input("cellcenter,icRegionInfo") ;
#line 352 "initialConditionRegions.loci"
       output("T_ic") ;
#line 352 "initialConditionRegions.loci"
       constraint("geom_cells") ;
#line 352 "initialConditionRegions.loci"
    }
#line 352 "initialConditionRegions.loci"
    void calculate(Entity _e_) { 
#line 353 "initialConditionRegions.loci"

    // Set the default state.
    L_T_ic_[_e_]=L_icRegionInfo_[_e_].defaultState .getTemperature () ;

    // Loop over regions, overwriting temperature with latest region value if the
    // cell center is in the region.
    for (size_t j =0;j <L_icRegionInfo_[_e_].fluidRegions .size ();++j ) {
      if (L_icRegionInfo_[_e_].fluidRegions [j ].geomTestFunc ->inGeomPt (L_cellcenter_[_e_])) {
        L_T_ic_[_e_]=L_icRegionInfo_[_e_].fluidRegions [j ].regionState .getTemperature () ;
      }
    }
  }    void compute(const Loci::sequence &seq) { 
#line 364 "initialConditionRegions.loci"
      do_loop(seq,this) ;
#line 364 "initialConditionRegions.loci"
    }
#line 364 "initialConditionRegions.loci"
} ;
#line 364 "initialConditionRegions.loci"
Loci::register_rule<file_initialConditionRegions004_1278520192m356> register_file_initialConditionRegions004_1278520192m356 ;
#line 364 "initialConditionRegions.loci"
}
#line 364 "initialConditionRegions.loci"


  // Initial condition for k.
  // $type k_ic store<real> 
  // $type omega_ic store<real> 
  namespace {class file_initialConditionRegions005_1278520192m357 : public Loci::pointwise_rule {
#line 369 "initialConditionRegions.loci"
    Loci::const_store<Loci::vector3d<Loci::real_t> >  L_cellcenter_ ; 
#line 369 "initialConditionRegions.loci"
    Loci::const_blackbox<ICparsedInitRegion>  L_icRegionInfo_ ; 
#line 369 "initialConditionRegions.loci"
    Loci::store<real>  L_k_ic_ ; 
#line 369 "initialConditionRegions.loci"
    Loci::store<real>  L_omega_ic_ ; 
#line 369 "initialConditionRegions.loci"
public:
#line 369 "initialConditionRegions.loci"
    file_initialConditionRegions005_1278520192m357() {
#line 369 "initialConditionRegions.loci"
       name_store("cellcenter",L_cellcenter_) ;
#line 369 "initialConditionRegions.loci"
       name_store("icRegionInfo",L_icRegionInfo_) ;
#line 369 "initialConditionRegions.loci"
       name_store("k_ic",L_k_ic_) ;
#line 369 "initialConditionRegions.loci"
       name_store("omega_ic",L_omega_ic_) ;
#line 369 "initialConditionRegions.loci"
       input("cellcenter,icRegionInfo") ;
#line 369 "initialConditionRegions.loci"
       output("k_ic") ;
#line 369 "initialConditionRegions.loci"
       output("omega_ic") ;
#line 369 "initialConditionRegions.loci"
       constraint("geom_cells") ;
#line 369 "initialConditionRegions.loci"
    }
#line 369 "initialConditionRegions.loci"
    void calculate(Entity _e_) { 
#line 370 "initialConditionRegions.loci"

    // Set the default state.
    L_icRegionInfo_[_e_].defaultState .get_k_omega (L_k_ic_[_e_],L_omega_ic_[_e_]) ;

    // Loop over regions, overwriting temperature with latest region value if the
    // cell center is in the region.
    for (size_t j =0;j <L_icRegionInfo_[_e_].fluidRegions .size ();++j ) {
      if (L_icRegionInfo_[_e_].fluidRegions [j ].geomTestFunc ->inGeomPt (L_cellcenter_[_e_])) {
        L_icRegionInfo_[_e_].fluidRegions [j ].regionState .get_k_omega (L_k_ic_[_e_],L_omega_ic_[_e_]) ;
      }
    }
  }    void compute(const Loci::sequence &seq) { 
#line 381 "initialConditionRegions.loci"
      do_loop(seq,this) ;
#line 381 "initialConditionRegions.loci"
    }
#line 381 "initialConditionRegions.loci"
} ;
#line 381 "initialConditionRegions.loci"
Loci::register_rule<file_initialConditionRegions005_1278520192m357> register_file_initialConditionRegions005_1278520192m357 ;
#line 381 "initialConditionRegions.loci"
}
#line 381 "initialConditionRegions.loci"


  // Initial condition for species mass fractions.
  // $type y_ic storeVec<real> 
  namespace {class file_initialConditionRegions006_1278520192m358 : public Loci::pointwise_rule {
#line 386 "initialConditionRegions.loci"
    Loci::const_store<Loci::vector3d<Loci::real_t> >  L_cellcenter_ ; 
#line 386 "initialConditionRegions.loci"
    Loci::const_param<fluidPhysics::EOS>  L_eos_ ; 
#line 386 "initialConditionRegions.loci"
    Loci::const_blackbox<ICparsedInitRegion>  L_icRegionInfo_ ; 
#line 386 "initialConditionRegions.loci"
    Loci::storeVec<real>  L_y_ic_ ; 
#line 386 "initialConditionRegions.loci"
public:
#line 386 "initialConditionRegions.loci"
    file_initialConditionRegions006_1278520192m358() {
#line 386 "initialConditionRegions.loci"
       name_store("cellcenter",L_cellcenter_) ;
#line 386 "initialConditionRegions.loci"
       name_store("eos",L_eos_) ;
#line 386 "initialConditionRegions.loci"
       name_store("icRegionInfo",L_icRegionInfo_) ;
#line 386 "initialConditionRegions.loci"
       name_store("y_ic",L_y_ic_) ;
#line 386 "initialConditionRegions.loci"
       input("eos,cellcenter,icRegionInfo") ;
#line 386 "initialConditionRegions.loci"
       output("y_ic") ;
#line 386 "initialConditionRegions.loci"
       constraint("geom_cells") ;
#line 386 "initialConditionRegions.loci"
    }
#line 386 "initialConditionRegions.loci"
    void prelude(const Loci::sequence &seq) { 
    L_y_ic_.setVecSize(L_eos_->numSpecies()) ;
  }    void calculate(Entity _e_) { 
#line 389 "initialConditionRegions.loci"


    // Set cell value to the default mixture. We must check each species for
    // validity since Ed did not do so during parsing. Initialize all mass
    // fractions to zero so unspecified values will have a default of zero.
    for (int i =0;i <L_eos_[_e_].numSpecies ();++i ) L_y_ic_[_e_][i ]=0.0 ;
    map <string ,real > defaultMixture =L_icRegionInfo_[_e_].defaultState .getMixture () ;
    for (map <string ,real >::const_iterator mi =defaultMixture .begin ();
    mi !=defaultMixture .end ();++mi ) {
      int speciesIndex =L_eos_[_e_].speciesIndex (mi ->first ) ;
      if (speciesIndex ==-1){
        cerr << "ERROR: Species " << mi ->first << " does not exist."
          << endl ; Loci ::Abort () ;
      }
      L_y_ic_[_e_][speciesIndex ]=mi ->second ;
    }

    // Loop over regions, overwriting species mass fractions with latest region
    // value if the cell center is in the region.
    for (size_t j =0;j <L_icRegionInfo_[_e_].fluidRegions .size ();++j ) {
      if (L_icRegionInfo_[_e_].fluidRegions [j ].geomTestFunc ->inGeomPt (L_cellcenter_[_e_])) {
        for (int i =0;i <L_eos_[_e_].numSpecies ();++i ) L_y_ic_[_e_][i ]=0.0 ;
        map <string ,real > mixture =L_icRegionInfo_[_e_].fluidRegions [j ].regionState .getMixture () ;
        for (map <string ,real >::const_iterator mi =mixture .begin ();
        mi !=mixture .end ();++mi ) {
          int speciesIndex =L_eos_[_e_].speciesIndex (mi ->first ) ;
          if (speciesIndex ==-1){
            cerr << "ERROR: Species " << mi ->first << " does not exist."
              << endl ; Loci ::Abort () ;
          }
          L_y_ic_[_e_][speciesIndex ]=mi ->second ;
        }
      }
    }
  }    void compute(const Loci::sequence &seq) { 
#line 423 "initialConditionRegions.loci"
      prelude(seq) ;
#line 423 "initialConditionRegions.loci"
      do_loop(seq,this) ;
#line 423 "initialConditionRegions.loci"
    }
#line 423 "initialConditionRegions.loci"
} ;
#line 423 "initialConditionRegions.loci"
Loci::register_rule<file_initialConditionRegions006_1278520192m358> register_file_initialConditionRegions006_1278520192m358 ;
#line 423 "initialConditionRegions.loci"
}
#line 423 "initialConditionRegions.loci"


}
