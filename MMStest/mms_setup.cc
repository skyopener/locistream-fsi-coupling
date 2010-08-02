#line 1 "mms_setup.loci"
#include <Loci.h>
#include <eos.h>
using fluidPhysics::EOS ;
#include <qvi.h>
#include <readGrid.h>

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
#line 7 "mms_setup.loci"

#line 1 "mms.lh"

// Hacks, these should be in chem.lh
// $type mixture_f storeVec<real> 
// $type numSpecies param<int> 

// mms stuff	
// $type massGeneralSrc(X0) storeVec<real> 
// $type momentumGeneralSrc(X0) store<vect3d> 
// $type energyGeneralSrc(X0) store<real> 

// $type SAGeneralSrc(X0) store<real> 

// $type rkGeneralSrc(X0) store<real> 
// $type rwGeneralSrc(X0) store<real> 

// $type mmsPressure(X0) store<real> 
// $type mmsRho(X0) storeVec<real> 
// $type mmsVelocity(X0) store<vect3d> 

// $type mms_nu_t(X0) store<real> 

// $type mms_k(X0) store<real> 
// $type mms_w(X0) store<real> 
#line 8 "mms_setup.loci"


// Setup for prescibing the mms variables.
namespace streamUns {

  // This bit of code allows us to add a mms option to the boundary conditions
  class mms_boundary_flag_check : public BC_Check {
    std::string error_message ;
  public:
    std::string BoundaryConditions() {
      string s("incompressibleInlet,subsonicInlet,supersonicInlet,") ;
      s+="fixedPressureOutlet,extrapolatedPressureOutlet" ;
      return s ;
    }
    std::string VariablesChecked(fact_db &facts) { return "mms" ; }
    bool CheckOptions(const options_list &bc_options,Loci::fact_db&) {
      error_message = "" ;
      return true ;
    }
    
    std::ostream &ErrorMessage(std::ostream &s) {
      s << error_message << endl ;
      return s ;
    }
  } ;

  register_BC<mms_boundary_flag_check> register_mms_boundary_flag_check ;
  
  // Add source term to continuity equation. Here we are adding all the species
  // continuity equations, which should work for both incompressible/compressible
  // flow with no species transport and compressible flows with species transport.
  // $type pPrime_B store<real> 
  namespace {class file_mms_setup000_1278521289m175 : public Loci::apply_rule< store<real> ,Loci::Summation<real> >  {
#line 40 "mms_setup.loci"
    Loci::const_store<Loci::real_t>  L_vol_ ; 
#line 40 "mms_setup.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 40 "mms_setup.loci"
    Loci::const_storeVec<real>  L_massGeneralSrccellcenter_ ; 
#line 40 "mms_setup.loci"
    Loci::store<real>  L_pPrime_B_ ; 
#line 40 "mms_setup.loci"
public:
#line 40 "mms_setup.loci"
    file_mms_setup000_1278521289m175() {
#line 40 "mms_setup.loci"
       name_store("vol",L_vol_) ;
#line 40 "mms_setup.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 40 "mms_setup.loci"
       name_store("pPrime_B",L_pPrime_B_) ;
#line 40 "mms_setup.loci"
       name_store("massGeneralSrc(cellcenter)",L_massGeneralSrccellcenter_) ;
#line 40 "mms_setup.loci"
       input("massGeneralSrc(cellcenter),numSpecies,vol") ;
#line 40 "mms_setup.loci"
       output("pPrime_B") ;
#line 40 "mms_setup.loci"
    }
#line 40 "mms_setup.loci"
    void calculate(Entity _e_) { 
#line 41 "mms_setup.loci"
    for (int i =0;i <L_numSpecies_[_e_];++i )
      L_pPrime_B_[_e_]+=L_massGeneralSrccellcenter_[_e_][i ]*L_vol_[_e_];
  }    void compute(const Loci::sequence &seq) { 
#line 43 "mms_setup.loci"
      do_loop(seq,this) ;
#line 43 "mms_setup.loci"
    }
#line 43 "mms_setup.loci"
} ;
#line 43 "mms_setup.loci"
Loci::register_rule<file_mms_setup000_1278521289m175> register_file_mms_setup000_1278521289m175 ;
#line 43 "mms_setup.loci"
}
#line 43 "mms_setup.loci"

  
  // Add source term to momentum equation.
  // $type vSourceTerm store<vect3d> 
  namespace {class file_mms_setup001_1278521289m176 : public Loci::apply_rule< store<vect3d> ,Loci::Summation<vect3d> >  {
#line 47 "mms_setup.loci"
    Loci::const_store<Loci::real_t>  L_vol_ ; 
#line 47 "mms_setup.loci"
    Loci::const_store<vect3d>  L_momentumGeneralSrccellcenter_ ; 
#line 47 "mms_setup.loci"
    Loci::store<vect3d>  L_vSourceTerm_ ; 
#line 47 "mms_setup.loci"
public:
#line 47 "mms_setup.loci"
    file_mms_setup001_1278521289m176() {
#line 47 "mms_setup.loci"
       name_store("vol",L_vol_) ;
#line 47 "mms_setup.loci"
       name_store("vSourceTerm",L_vSourceTerm_) ;
#line 47 "mms_setup.loci"
       name_store("momentumGeneralSrc(cellcenter)",L_momentumGeneralSrccellcenter_) ;
#line 47 "mms_setup.loci"
       input("momentumGeneralSrc(cellcenter),vol") ;
#line 47 "mms_setup.loci"
       output("vSourceTerm") ;
#line 47 "mms_setup.loci"
    }
#line 47 "mms_setup.loci"
    void calculate(Entity _e_) { 
#line 48 "mms_setup.loci"
    vect3d tmp =L_vol_[_e_]*L_momentumGeneralSrccellcenter_[_e_];
    L_vSourceTerm_[_e_]+=L_vol_[_e_]*L_momentumGeneralSrccellcenter_[_e_];
  }    void compute(const Loci::sequence &seq) { 
#line 50 "mms_setup.loci"
      do_loop(seq,this) ;
#line 50 "mms_setup.loci"
    }
#line 50 "mms_setup.loci"
} ;
#line 50 "mms_setup.loci"
Loci::register_rule<file_mms_setup001_1278521289m176> register_file_mms_setup001_1278521289m176 ;
#line 50 "mms_setup.loci"
}
#line 50 "mms_setup.loci"


  // Integrate kw turbulence model source term over general cell
  namespace {class file_mms_setup002_1278521289m177 : public Loci::apply_rule< store<streamUns::vec<2> > ,Loci::Summation<streamUns::vec<2> > >  {
#line 54 "mms_setup.loci"
    Loci::const_store<Loci::real_t>  L_vol_ ; 
#line 54 "mms_setup.loci"
    Loci::const_store<real>  L_rkGeneralSrccellcenter_ ; 
#line 54 "mms_setup.loci"
    Loci::const_store<real>  L_rwGeneralSrccellcenter_ ; 
#line 54 "mms_setup.loci"
    Loci::store<streamUns::vec<2> >  L_sst_src_ ; 
#line 54 "mms_setup.loci"
public:
#line 54 "mms_setup.loci"
    file_mms_setup002_1278521289m177() {
#line 54 "mms_setup.loci"
       name_store("vol",L_vol_) ;
#line 54 "mms_setup.loci"
       name_store("sst_src",L_sst_src_) ;
#line 54 "mms_setup.loci"
       name_store("rkGeneralSrc(cellcenter)",L_rkGeneralSrccellcenter_) ;
#line 54 "mms_setup.loci"
       name_store("rwGeneralSrc(cellcenter)",L_rwGeneralSrccellcenter_) ;
#line 54 "mms_setup.loci"
       input("rkGeneralSrc(cellcenter),                       rwGeneralSrc(cellcenter),vol") ;
#line 54 "mms_setup.loci"
       output("sst_src") ;
#line 54 "mms_setup.loci"
    }
#line 54 "mms_setup.loci"
    void calculate(Entity _e_) { 
#line 55 "mms_setup.loci"
    L_sst_src_[_e_][0]+= L_rkGeneralSrccellcenter_[_e_]*L_vol_[_e_];
    L_sst_src_[_e_][1]+= L_rwGeneralSrccellcenter_[_e_]*L_vol_[_e_];
  }    void compute(const Loci::sequence &seq) { 
#line 57 "mms_setup.loci"
      do_loop(seq,this) ;
#line 57 "mms_setup.loci"
    }
#line 57 "mms_setup.loci"
} ;
#line 57 "mms_setup.loci"
Loci::register_rule<file_mms_setup002_1278521289m177> register_file_mms_setup002_1278521289m177 ;
#line 57 "mms_setup.loci"
}
#line 57 "mms_setup.loci"


  // Provide boundary values for pressure for compressible flows.
  // $type p_f store<real> 
  namespace {class file_mms_setup003_1278521289m178 : public Loci::pointwise_rule {
#line 62 "mms_setup.loci"
    Loci::const_store<real>  L_mmsPressurefacecenter_ ; 
#line 62 "mms_setup.loci"
    Loci::store<real>  L_p_f_ ; 
#line 62 "mms_setup.loci"
public:
#line 62 "mms_setup.loci"
    file_mms_setup003_1278521289m178() {
#line 62 "mms_setup.loci"
       name_store("p_f",L_p_f_) ;
#line 62 "mms_setup.loci"
       name_store("mmsPressure(facecenter)",L_mmsPressurefacecenter_) ;
#line 62 "mms_setup.loci"
       input("mmsPressure(facecenter)") ;
#line 62 "mms_setup.loci"
       output("p_f") ;
#line 62 "mms_setup.loci"
       constraint("ref->mms_BCoption,compressibleFlow") ;
#line 62 "mms_setup.loci"
    }
#line 62 "mms_setup.loci"
    void calculate(Entity _e_) { 
#line 63 "mms_setup.loci"
    L_p_f_[_e_]=L_mmsPressurefacecenter_[_e_];
  }    void compute(const Loci::sequence &seq) { 
#line 64 "mms_setup.loci"
      do_loop(seq,this) ;
#line 64 "mms_setup.loci"
    }
#line 64 "mms_setup.loci"
} ;
#line 64 "mms_setup.loci"
Loci::register_rule<file_mms_setup003_1278521289m178> register_file_mms_setup003_1278521289m178 ;
#line 64 "mms_setup.loci"
}
#line 64 "mms_setup.loci"


  // Provide boundary values for velocity.
  // $type v_f store<vect3d> 
  namespace {class file_mms_setup004_1278521289m178 : public Loci::pointwise_rule {
#line 69 "mms_setup.loci"
    Loci::const_store<vect3d>  L_mmsVelocityfacecenter_ ; 
#line 69 "mms_setup.loci"
    Loci::store<vect3d>  L_v_f_ ; 
#line 69 "mms_setup.loci"
public:
#line 69 "mms_setup.loci"
    file_mms_setup004_1278521289m178() {
#line 69 "mms_setup.loci"
       name_store("v_f",L_v_f_) ;
#line 69 "mms_setup.loci"
       name_store("mmsVelocity(facecenter)",L_mmsVelocityfacecenter_) ;
#line 69 "mms_setup.loci"
       input("mmsVelocity(facecenter)") ;
#line 69 "mms_setup.loci"
       output("v_f") ;
#line 69 "mms_setup.loci"
       constraint("ref->mms_BCoption") ;
#line 69 "mms_setup.loci"
    }
#line 69 "mms_setup.loci"
    void calculate(Entity _e_) { 
#line 70 "mms_setup.loci"
    L_v_f_[_e_]=L_mmsVelocityfacecenter_[_e_];
  }    void compute(const Loci::sequence &seq) { 
#line 71 "mms_setup.loci"
      do_loop(seq,this) ;
#line 71 "mms_setup.loci"
    }
#line 71 "mms_setup.loci"
} ;
#line 71 "mms_setup.loci"
Loci::register_rule<file_mms_setup004_1278521289m178> register_file_mms_setup004_1278521289m178 ;
#line 71 "mms_setup.loci"
}
#line 71 "mms_setup.loci"


  // Provide boundary values for temperature for compressible flows.
  namespace {class file_mms_setup005_1278521289m179 : public Loci::pointwise_rule {
#line 75 "mms_setup.loci"
    Loci::const_param<fluidPhysics::EOS>  L_eos_ ; 
#line 75 "mms_setup.loci"
    Loci::const_store<real>  L_mmsPressurefacecenter_ ; 
#line 75 "mms_setup.loci"
    Loci::const_storeVec<real>  L_mmsRhofacecenter_ ; 
#line 75 "mms_setup.loci"
    Loci::store<streamUns::real>  L_temperature_f_ ; 
#line 75 "mms_setup.loci"
public:
#line 75 "mms_setup.loci"
    file_mms_setup005_1278521289m179() {
#line 75 "mms_setup.loci"
       name_store("eos",L_eos_) ;
#line 75 "mms_setup.loci"
       name_store("temperature_f",L_temperature_f_) ;
#line 75 "mms_setup.loci"
       name_store("mmsPressure(facecenter)",L_mmsPressurefacecenter_) ;
#line 75 "mms_setup.loci"
       name_store("mmsRho(facecenter)",L_mmsRhofacecenter_) ;
#line 75 "mms_setup.loci"
       input("mmsPressure(facecenter),mmsRho(facecenter),eos") ;
#line 75 "mms_setup.loci"
       output("temperature_f") ;
#line 75 "mms_setup.loci"
       constraint("ref->mms_BCoption,compressibleFlow") ;
#line 75 "mms_setup.loci"
    }
#line 75 "mms_setup.loci"
    void calculate(Entity _e_) { 
#line 76 "mms_setup.loci"
    EOS ::State s = L_eos_[_e_].State_from_rho_p (L_mmsRhofacecenter_[_e_],
                                         L_mmsPressurefacecenter_[_e_]) ;
    L_temperature_f_[_e_]= s .temperature () ;
  }    void compute(const Loci::sequence &seq) { 
#line 79 "mms_setup.loci"
      do_loop(seq,this) ;
#line 79 "mms_setup.loci"
    }
#line 79 "mms_setup.loci"
} ;
#line 79 "mms_setup.loci"
Loci::register_rule<file_mms_setup005_1278521289m179> register_file_mms_setup005_1278521289m179 ;
#line 79 "mms_setup.loci"
}
#line 79 "mms_setup.loci"


  // Provide boundary values for mixture
  namespace {class file_mms_setup006_1278521289m179 : public Loci::pointwise_rule {
#line 83 "mms_setup.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 83 "mms_setup.loci"
    Loci::const_storeVec<real>  L_mmsRhofacecenter_ ; 
#line 83 "mms_setup.loci"
    Loci::storeVec<real>  L_mmsmixture_f_ ; 
#line 83 "mms_setup.loci"
public:
#line 83 "mms_setup.loci"
    file_mms_setup006_1278521289m179() {
#line 83 "mms_setup.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 83 "mms_setup.loci"
       name_store("mmsRho(facecenter)",L_mmsRhofacecenter_) ;
#line 83 "mms_setup.loci"
       name_store("mms::mixture_f",L_mmsmixture_f_) ;
#line 83 "mms_setup.loci"
       input("mmsRho(facecenter),numSpecies") ;
#line 83 "mms_setup.loci"
       output("mms::mixture_f") ;
#line 83 "mms_setup.loci"
       constraint("ref->mms_BCoption") ;
#line 83 "mms_setup.loci"
    }
#line 83 "mms_setup.loci"
    void prelude(const Loci::sequence &seq) { 
    L_mmsmixture_f_.setVecSize(*L_numSpecies_) ;
  }    void calculate(Entity _e_) { 
#line 86 "mms_setup.loci"
    real rho = 0 ;
    for (int i =0;i <L_numSpecies_[_e_];++i )
      rho += L_mmsRhofacecenter_[_e_][i ];
    for (int i =0;i <L_numSpecies_[_e_];++i )
      L_mmsmixture_f_[_e_][i ]= L_mmsRhofacecenter_[_e_][i ]/rho ;
  }    void compute(const Loci::sequence &seq) { 
#line 91 "mms_setup.loci"
      prelude(seq) ;
#line 91 "mms_setup.loci"
      do_loop(seq,this) ;
#line 91 "mms_setup.loci"
    }
#line 91 "mms_setup.loci"
} ;
#line 91 "mms_setup.loci"
Loci::register_rule<file_mms_setup006_1278521289m179> register_file_mms_setup006_1278521289m179 ;
#line 91 "mms_setup.loci"
}
#line 91 "mms_setup.loci"


  // Provide boundary values for spalart allmaras turbulence variable
  namespace {class file_mms_setup007_1278521289m180 : public Loci::pointwise_rule {
#line 95 "mms_setup.loci"
    Loci::const_store<real>  L_mms_nu_tfacecenter_ ; 
#line 95 "mms_setup.loci"
    Loci::store<streamUns::real>  L_mmsnu_t_f_ ; 
#line 95 "mms_setup.loci"
public:
#line 95 "mms_setup.loci"
    file_mms_setup007_1278521289m180() {
#line 95 "mms_setup.loci"
       name_store("mms_nu_t(facecenter)",L_mms_nu_tfacecenter_) ;
#line 95 "mms_setup.loci"
       name_store("mms::nu_t_f",L_mmsnu_t_f_) ;
#line 95 "mms_setup.loci"
       input("mms_nu_t(facecenter)") ;
#line 95 "mms_setup.loci"
       output("mms::nu_t_f") ;
#line 95 "mms_setup.loci"
       constraint("ref->mms_BCoption") ;
#line 95 "mms_setup.loci"
    }
#line 95 "mms_setup.loci"
    void calculate(Entity _e_) { 
#line 96 "mms_setup.loci"
    L_mmsnu_t_f_[_e_]= L_mms_nu_tfacecenter_[_e_];
  }    void compute(const Loci::sequence &seq) { 
#line 97 "mms_setup.loci"
      do_loop(seq,this) ;
#line 97 "mms_setup.loci"
    }
#line 97 "mms_setup.loci"
} ;
#line 97 "mms_setup.loci"
Loci::register_rule<file_mms_setup007_1278521289m180> register_file_mms_setup007_1278521289m180 ;
#line 97 "mms_setup.loci"
}
#line 97 "mms_setup.loci"


  // Provide turbulent kintentic energy at boundary
  namespace {class file_mms_setup008_1278521289m181 : public Loci::pointwise_rule {
#line 101 "mms_setup.loci"
    Loci::const_store<real>  L_mms_kfacecenter_ ; 
#line 101 "mms_setup.loci"
    Loci::store<streamUns::real>  L_mmsk_f_ ; 
#line 101 "mms_setup.loci"
public:
#line 101 "mms_setup.loci"
    file_mms_setup008_1278521289m181() {
#line 101 "mms_setup.loci"
       name_store("mms_k(facecenter)",L_mms_kfacecenter_) ;
#line 101 "mms_setup.loci"
       name_store("mms::k_f",L_mmsk_f_) ;
#line 101 "mms_setup.loci"
       input("mms_k(facecenter)") ;
#line 101 "mms_setup.loci"
       output("mms::k_f") ;
#line 101 "mms_setup.loci"
       constraint("ref->mms_BCoption") ;
#line 101 "mms_setup.loci"
    }
#line 101 "mms_setup.loci"
    void calculate(Entity _e_) { 
#line 102 "mms_setup.loci"
    L_mmsk_f_[_e_]= L_mms_kfacecenter_[_e_];
  }    void compute(const Loci::sequence &seq) { 
#line 103 "mms_setup.loci"
      do_loop(seq,this) ;
#line 103 "mms_setup.loci"
    }
#line 103 "mms_setup.loci"
} ;
#line 103 "mms_setup.loci"
Loci::register_rule<file_mms_setup008_1278521289m181> register_file_mms_setup008_1278521289m181 ;
#line 103 "mms_setup.loci"
}
#line 103 "mms_setup.loci"


  // Provide turbulent omega at boundary
  namespace {class file_mms_setup009_1278521289m181 : public Loci::pointwise_rule {
#line 107 "mms_setup.loci"
    Loci::const_store<real>  L_mms_wfacecenter_ ; 
#line 107 "mms_setup.loci"
    Loci::store<streamUns::real>  L_mmsw_f_ ; 
#line 107 "mms_setup.loci"
public:
#line 107 "mms_setup.loci"
    file_mms_setup009_1278521289m181() {
#line 107 "mms_setup.loci"
       name_store("mms_w(facecenter)",L_mms_wfacecenter_) ;
#line 107 "mms_setup.loci"
       name_store("mms::w_f",L_mmsw_f_) ;
#line 107 "mms_setup.loci"
       input("mms_w(facecenter)") ;
#line 107 "mms_setup.loci"
       output("mms::w_f") ;
#line 107 "mms_setup.loci"
       constraint("ref->mms_BCoption") ;
#line 107 "mms_setup.loci"
    }
#line 107 "mms_setup.loci"
    void calculate(Entity _e_) { 
#line 108 "mms_setup.loci"
    L_mmsw_f_[_e_]= L_mms_wfacecenter_[_e_];
  }    void compute(const Loci::sequence &seq) { 
#line 109 "mms_setup.loci"
      do_loop(seq,this) ;
#line 109 "mms_setup.loci"
    }
#line 109 "mms_setup.loci"
} ;
#line 109 "mms_setup.loci"
Loci::register_rule<file_mms_setup009_1278521289m181> register_file_mms_setup009_1278521289m181 ;
#line 109 "mms_setup.loci"
}
#line 109 "mms_setup.loci"


  // Compute mms error
  // $type mmsMassResidual store<real> 
  // $type mmsMomentXResidual store<real> 
  // $type mmsMomentYResidual store<real> 
  // $type mmsMomentZResidual store<real> 
  // $type mmsEnergyResidual store<real> 

  // Here we compute the error between the computed conservative variables and
  // the predicted conservative variables
  // $type rho store<real> 
  namespace {class file_mms_setup010_1278521289m182 : public Loci::pointwise_rule {
#line 121 "mms_setup.loci"
    Loci::const_store<real>  L_rho_ ; 
#line 121 "mms_setup.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 121 "mms_setup.loci"
    Loci::const_storeVec<real>  L_mmsRhocellcenter_ ; 
#line 121 "mms_setup.loci"
    Loci::store<real>  L_mmsMassResidual_ ; 
#line 121 "mms_setup.loci"
public:
#line 121 "mms_setup.loci"
    file_mms_setup010_1278521289m182() {
#line 121 "mms_setup.loci"
       name_store("rho",L_rho_) ;
#line 121 "mms_setup.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 121 "mms_setup.loci"
       name_store("mmsMassResidual",L_mmsMassResidual_) ;
#line 121 "mms_setup.loci"
       name_store("mmsRho(cellcenter)",L_mmsRhocellcenter_) ;
#line 121 "mms_setup.loci"
       input("rho,mmsRho(cellcenter),numSpecies") ;
#line 121 "mms_setup.loci"
       output("mmsMassResidual") ;
#line 121 "mms_setup.loci"
    }
#line 121 "mms_setup.loci"
    void calculate(Entity _e_) { 
#line 122 "mms_setup.loci"
    real mmsrho = 0 ;
    for (int i =0;i <L_numSpecies_[_e_];++i )
      mmsrho += L_mmsRhocellcenter_[_e_][i ];
    L_mmsMassResidual_[_e_]= L_rho_[_e_]-mmsrho ;
  }    void compute(const Loci::sequence &seq) { 
#line 126 "mms_setup.loci"
      do_loop(seq,this) ;
#line 126 "mms_setup.loci"
    }
#line 126 "mms_setup.loci"
} ;
#line 126 "mms_setup.loci"
Loci::register_rule<file_mms_setup010_1278521289m182> register_file_mms_setup010_1278521289m182 ;
#line 126 "mms_setup.loci"
}
#line 126 "mms_setup.loci"


  // $type v store<vect3d> 
  namespace {class file_mms_setup011_1278521289m183 : public Loci::pointwise_rule {
#line 130 "mms_setup.loci"
    Loci::const_store<real>  L_rho_ ; 
#line 130 "mms_setup.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 130 "mms_setup.loci"
    Loci::const_storeVec<real>  L_mmsRhocellcenter_ ; 
#line 130 "mms_setup.loci"
    Loci::const_store<vect3d>  L_v_ ; 
#line 130 "mms_setup.loci"
    Loci::const_store<vect3d>  L_mmsVelocitycellcenter_ ; 
#line 130 "mms_setup.loci"
    Loci::store<real>  L_mmsMomentXResidual_ ; 
#line 130 "mms_setup.loci"
    Loci::store<real>  L_mmsMomentYResidual_ ; 
#line 130 "mms_setup.loci"
    Loci::store<real>  L_mmsMomentZResidual_ ; 
#line 130 "mms_setup.loci"
public:
#line 130 "mms_setup.loci"
    file_mms_setup011_1278521289m183() {
#line 130 "mms_setup.loci"
       name_store("rho",L_rho_) ;
#line 130 "mms_setup.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 130 "mms_setup.loci"
       name_store("mmsMomentXResidual",L_mmsMomentXResidual_) ;
#line 130 "mms_setup.loci"
       name_store("mmsMomentYResidual",L_mmsMomentYResidual_) ;
#line 130 "mms_setup.loci"
       name_store("mmsMomentZResidual",L_mmsMomentZResidual_) ;
#line 130 "mms_setup.loci"
       name_store("mmsRho(cellcenter)",L_mmsRhocellcenter_) ;
#line 130 "mms_setup.loci"
       name_store("v",L_v_) ;
#line 130 "mms_setup.loci"
       name_store("mmsVelocity(cellcenter)",L_mmsVelocitycellcenter_) ;
#line 130 "mms_setup.loci"
       input("rho,v,  mmsRho(cellcenter),mmsVelocity(cellcenter),numSpecies") ;
#line 130 "mms_setup.loci"
       output("mmsMomentXResidual") ;
#line 130 "mms_setup.loci"
       output("mmsMomentYResidual") ;
#line 130 "mms_setup.loci"
       output("mmsMomentZResidual") ;
#line 130 "mms_setup.loci"
    }
#line 130 "mms_setup.loci"
    void calculate(Entity _e_) { 
#line 131 "mms_setup.loci"
    real mmsrho = 0 ;
    for (int i =0;i <L_numSpecies_[_e_];++i )
      mmsrho += L_mmsRhocellcenter_[_e_][i ];
    vect3d resid = L_rho_[_e_]*L_v_[_e_]- mmsrho *L_mmsVelocitycellcenter_[_e_];
    L_mmsMomentXResidual_[_e_]= resid .x ;
    L_mmsMomentYResidual_[_e_]= resid .y ;
    L_mmsMomentZResidual_[_e_]= resid .z ;
  }    void compute(const Loci::sequence &seq) { 
#line 138 "mms_setup.loci"
      do_loop(seq,this) ;
#line 138 "mms_setup.loci"
    }
#line 138 "mms_setup.loci"
} ;
#line 138 "mms_setup.loci"
Loci::register_rule<file_mms_setup011_1278521289m183> register_file_mms_setup011_1278521289m183 ;
#line 138 "mms_setup.loci"
}
#line 138 "mms_setup.loci"


  namespace {class file_mms_setup012_1278521289m184 : public Loci::pointwise_rule {
#line 140 "mms_setup.loci"
    Loci::store<real>  L_mmsEnergyResidual_ ; 
#line 140 "mms_setup.loci"
public:
#line 140 "mms_setup.loci"
    file_mms_setup012_1278521289m184() {
#line 140 "mms_setup.loci"
       name_store("mmsEnergyResidual",L_mmsEnergyResidual_) ;
#line 140 "mms_setup.loci"
       output("mmsEnergyResidual") ;
#line 140 "mms_setup.loci"
       constraint("incompressibleFlow,geom_cells") ;
#line 140 "mms_setup.loci"
    }
#line 140 "mms_setup.loci"
    void calculate(Entity _e_) { 
#line 141 "mms_setup.loci"
    L_mmsEnergyResidual_[_e_]= 0 ;
  }    void compute(const Loci::sequence &seq) { 
#line 142 "mms_setup.loci"
      do_loop(seq,this) ;
#line 142 "mms_setup.loci"
    }
#line 142 "mms_setup.loci"
} ;
#line 142 "mms_setup.loci"
Loci::register_rule<file_mms_setup012_1278521289m184> register_file_mms_setup012_1278521289m184 ;
#line 142 "mms_setup.loci"
}
#line 142 "mms_setup.loci"


  namespace {class file_mms_setup013_1278521289m184 : public Loci::pointwise_rule {
#line 145 "mms_setup.loci"
    Loci::const_param<fluidPhysics::EOS>  L_eos_ ; 
#line 145 "mms_setup.loci"
    Loci::const_store<real>  L_rho_ ; 
#line 145 "mms_setup.loci"
    Loci::const_store<streamUns::real>  L_e_internal_ ; 
#line 145 "mms_setup.loci"
    Loci::const_store<streamUns::vect3d>  L_u_ ; 
#line 145 "mms_setup.loci"
    Loci::const_storeVec<real>  L_mmsRhocellcenter_ ; 
#line 145 "mms_setup.loci"
    Loci::const_store<vect3d>  L_mmsVelocitycellcenter_ ; 
#line 145 "mms_setup.loci"
    Loci::const_store<real>  L_mmsPressurecellcenter_ ; 
#line 145 "mms_setup.loci"
    Loci::store<real>  L_mmsEnergyResidual_ ; 
#line 145 "mms_setup.loci"
public:
#line 145 "mms_setup.loci"
    file_mms_setup013_1278521289m184() {
#line 145 "mms_setup.loci"
       name_store("eos",L_eos_) ;
#line 145 "mms_setup.loci"
       name_store("rho",L_rho_) ;
#line 145 "mms_setup.loci"
       name_store("e_internal",L_e_internal_) ;
#line 145 "mms_setup.loci"
       name_store("u",L_u_) ;
#line 145 "mms_setup.loci"
       name_store("mmsEnergyResidual",L_mmsEnergyResidual_) ;
#line 145 "mms_setup.loci"
       name_store("mmsRho(cellcenter)",L_mmsRhocellcenter_) ;
#line 145 "mms_setup.loci"
       name_store("mmsVelocity(cellcenter)",L_mmsVelocitycellcenter_) ;
#line 145 "mms_setup.loci"
       name_store("mmsPressure(cellcenter)",L_mmsPressurecellcenter_) ;
#line 145 "mms_setup.loci"
       input("e_internal,u,rho,mmsRho(cellcenter),mmsPressure(cellcenter),  mmsVelocity(cellcenter),eos") ;
#line 145 "mms_setup.loci"
       output("mmsEnergyResidual") ;
#line 145 "mms_setup.loci"
       constraint("compressibleFlow,geom_cells") ;
#line 145 "mms_setup.loci"
    }
#line 145 "mms_setup.loci"
    void calculate(Entity _e_) { 
#line 146 "mms_setup.loci"
    EOS ::State s = L_eos_[_e_].State_from_rho_p (L_mmsRhocellcenter_[_e_],L_mmsPressurecellcenter_[_e_]) ;
    real re0_mms = s .density ()*(s .energy ()+.5*dot (L_u_[_e_],L_u_[_e_])) ;
    real re0_solver = L_rho_[_e_]*(L_e_internal_[_e_]+.5*dot (L_u_[_e_],L_u_[_e_])) ;
    L_mmsEnergyResidual_[_e_]= re0_solver - re0_mms ;
  }    void compute(const Loci::sequence &seq) { 
#line 150 "mms_setup.loci"
      do_loop(seq,this) ;
#line 150 "mms_setup.loci"
    }
#line 150 "mms_setup.loci"
} ;
#line 150 "mms_setup.loci"
Loci::register_rule<file_mms_setup013_1278521289m184> register_file_mms_setup013_1278521289m184 ;
#line 150 "mms_setup.loci"
}
#line 150 "mms_setup.loci"


  // $type mms_rk_Residual store<real> 
  // $type mms_rw_Residual store<real> 

  namespace {class file_mms_setup014_1278521289m185 : public Loci::pointwise_rule {
#line 155 "mms_setup.loci"
    Loci::const_store<real>  L_rho_ ; 
#line 155 "mms_setup.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 155 "mms_setup.loci"
    Loci::const_store<streamUns::real>  L_k_ ; 
#line 155 "mms_setup.loci"
    Loci::const_storeVec<real>  L_mmsRhocellcenter_ ; 
#line 155 "mms_setup.loci"
    Loci::const_store<real>  L_mms_kcellcenter_ ; 
#line 155 "mms_setup.loci"
    Loci::store<real>  L_mms_rk_Residual_ ; 
#line 155 "mms_setup.loci"
public:
#line 155 "mms_setup.loci"
    file_mms_setup014_1278521289m185() {
#line 155 "mms_setup.loci"
       name_store("rho",L_rho_) ;
#line 155 "mms_setup.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 155 "mms_setup.loci"
       name_store("k",L_k_) ;
#line 155 "mms_setup.loci"
       name_store("mmsRho(cellcenter)",L_mmsRhocellcenter_) ;
#line 155 "mms_setup.loci"
       name_store("mms_rk_Residual",L_mms_rk_Residual_) ;
#line 155 "mms_setup.loci"
       name_store("mms_k(cellcenter)",L_mms_kcellcenter_) ;
#line 155 "mms_setup.loci"
       input("rho,k,mms_k(cellcenter),mmsRho(cellcenter),numSpecies") ;
#line 155 "mms_setup.loci"
       output("mms_rk_Residual") ;
#line 155 "mms_setup.loci"
    }
#line 155 "mms_setup.loci"
    void calculate(Entity _e_) { 
#line 156 "mms_setup.loci"
    real mmsrho = 0 ;
    for (int i =0;i <L_numSpecies_[_e_];++i )
      mmsrho += L_mmsRhocellcenter_[_e_][i ];
    L_mms_rk_Residual_[_e_]= L_rho_[_e_]*L_k_[_e_]- L_mms_kcellcenter_[_e_]*mmsrho ;
  }    void compute(const Loci::sequence &seq) { 
#line 160 "mms_setup.loci"
      do_loop(seq,this) ;
#line 160 "mms_setup.loci"
    }
#line 160 "mms_setup.loci"
} ;
#line 160 "mms_setup.loci"
Loci::register_rule<file_mms_setup014_1278521289m185> register_file_mms_setup014_1278521289m185 ;
#line 160 "mms_setup.loci"
}
#line 160 "mms_setup.loci"


  // $type omega store<real> 
  namespace {class file_mms_setup015_1278521289m186 : public Loci::pointwise_rule {
#line 163 "mms_setup.loci"
    Loci::const_store<real>  L_rho_ ; 
#line 163 "mms_setup.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 163 "mms_setup.loci"
    Loci::const_storeVec<real>  L_mmsRhocellcenter_ ; 
#line 163 "mms_setup.loci"
    Loci::const_store<real>  L_omega_ ; 
#line 163 "mms_setup.loci"
    Loci::const_store<real>  L_mms_wcellcenter_ ; 
#line 163 "mms_setup.loci"
    Loci::store<real>  L_mms_rw_Residual_ ; 
#line 163 "mms_setup.loci"
public:
#line 163 "mms_setup.loci"
    file_mms_setup015_1278521289m186() {
#line 163 "mms_setup.loci"
       name_store("rho",L_rho_) ;
#line 163 "mms_setup.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 163 "mms_setup.loci"
       name_store("mmsRho(cellcenter)",L_mmsRhocellcenter_) ;
#line 163 "mms_setup.loci"
       name_store("mms_rw_Residual",L_mms_rw_Residual_) ;
#line 163 "mms_setup.loci"
       name_store("omega",L_omega_) ;
#line 163 "mms_setup.loci"
       name_store("mms_w(cellcenter)",L_mms_wcellcenter_) ;
#line 163 "mms_setup.loci"
       input("rho,omega,mms_w(cellcenter),mmsRho(cellcenter),numSpecies") ;
#line 163 "mms_setup.loci"
       output("mms_rw_Residual") ;
#line 163 "mms_setup.loci"
    }
#line 163 "mms_setup.loci"
    void calculate(Entity _e_) { 
#line 164 "mms_setup.loci"
    real mmsrho = 0 ;
    for (int i =0;i <L_numSpecies_[_e_];++i )
      mmsrho += L_mmsRhocellcenter_[_e_][i ];
    L_mms_rw_Residual_[_e_]= L_rho_[_e_]*L_omega_[_e_]- L_mms_wcellcenter_[_e_]*mmsrho ;
  }    void compute(const Loci::sequence &seq) { 
#line 168 "mms_setup.loci"
      do_loop(seq,this) ;
#line 168 "mms_setup.loci"
    }
#line 168 "mms_setup.loci"
} ;
#line 168 "mms_setup.loci"
Loci::register_rule<file_mms_setup015_1278521289m186> register_file_mms_setup015_1278521289m186 ;
#line 168 "mms_setup.loci"
}
#line 168 "mms_setup.loci"


  namespace {class file_mms_setup016_1278521289m187 : public Loci::singleton_rule {
#line 172 "mms_setup.loci"
    Loci::const_param<string>  L_modelName_ ; 
#line 172 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_L2NormmmsMassResidual_ ; 
#line 172 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_L2NormmmsMomentXResidual_ ; 
#line 172 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_L2NormmmsMomentYResidual_ ; 
#line 172 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_L2NormmmsMomentZResidual_ ; 
#line 172 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_L2NormmmsEnergyResidual_ ; 
#line 172 "mms_setup.loci"
    Loci::param<bool> L_OUTPUT_ ; 
#line 172 "mms_setup.loci"
public:
#line 172 "mms_setup.loci"
    file_mms_setup016_1278521289m187() {
#line 172 "mms_setup.loci"
       name_store("OUTPUT",L_OUTPUT_) ;
#line 172 "mms_setup.loci"
       name_store("modelName",L_modelName_) ;
#line 172 "mms_setup.loci"
       name_store("L2Norm(mmsMassResidual)",L_L2NormmmsMassResidual_) ;
#line 172 "mms_setup.loci"
       name_store("L2Norm(mmsMomentXResidual)",L_L2NormmmsMomentXResidual_) ;
#line 172 "mms_setup.loci"
       name_store("L2Norm(mmsMomentYResidual)",L_L2NormmmsMomentYResidual_) ;
#line 172 "mms_setup.loci"
       name_store("L2Norm(mmsMomentZResidual)",L_L2NormmmsMomentZResidual_) ;
#line 172 "mms_setup.loci"
       name_store("L2Norm(mmsEnergyResidual)",L_L2NormmmsEnergyResidual_) ;
#line 172 "mms_setup.loci"
       input("L2Norm(mmsMassResidual),L2Norm(mmsMomentXResidual),  L2Norm(mmsMomentYResidual),L2Norm(mmsMomentZResidual),L2Norm(mmsEnergyResidual),  modelName") ;
#line 172 "mms_setup.loci"
       output("OUTPUT") ;
#line 172 "mms_setup.loci"
       conditional("do_plot") ;
#line 172 "mms_setup.loci"
    }
#line 172 "mms_setup.loci"
    void compute(const Loci::sequence &seq) { 
    if(Loci::MPI_rank != 0) return ;

    string fname = (*L_modelName_)+"_meanl2.dat" ;
    ofstream ofile(fname.c_str(),ios::out) ;
    ofile.precision(16) ;
    ofile << (*L_L2NormmmsMassResidual_)<< ' '
          << (*L_L2NormmmsMomentXResidual_)<< ' '
          << (*L_L2NormmmsMomentYResidual_)<< ' '
          << (*L_L2NormmmsMomentZResidual_)<< ' '
          << (*L_L2NormmmsEnergyResidual_)<< endl ;
  }} ;
#line 183 "mms_setup.loci"
Loci::register_rule<file_mms_setup016_1278521289m187> register_file_mms_setup016_1278521289m187 ;
#line 183 "mms_setup.loci"
}
#line 183 "mms_setup.loci"


  namespace {class file_mms_setup017_1278521289m188 : public Loci::singleton_rule {
#line 187 "mms_setup.loci"
    Loci::const_param<string>  L_modelName_ ; 
#line 187 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_L1NormmmsMassResidual_ ; 
#line 187 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_L1NormmmsMomentXResidual_ ; 
#line 187 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_L1NormmmsMomentYResidual_ ; 
#line 187 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_L1NormmmsMomentZResidual_ ; 
#line 187 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_L1NormmmsEnergyResidual_ ; 
#line 187 "mms_setup.loci"
    Loci::param<bool> L_OUTPUT_ ; 
#line 187 "mms_setup.loci"
public:
#line 187 "mms_setup.loci"
    file_mms_setup017_1278521289m188() {
#line 187 "mms_setup.loci"
       name_store("OUTPUT",L_OUTPUT_) ;
#line 187 "mms_setup.loci"
       name_store("modelName",L_modelName_) ;
#line 187 "mms_setup.loci"
       name_store("L1Norm(mmsMassResidual)",L_L1NormmmsMassResidual_) ;
#line 187 "mms_setup.loci"
       name_store("L1Norm(mmsMomentXResidual)",L_L1NormmmsMomentXResidual_) ;
#line 187 "mms_setup.loci"
       name_store("L1Norm(mmsMomentYResidual)",L_L1NormmmsMomentYResidual_) ;
#line 187 "mms_setup.loci"
       name_store("L1Norm(mmsMomentZResidual)",L_L1NormmmsMomentZResidual_) ;
#line 187 "mms_setup.loci"
       name_store("L1Norm(mmsEnergyResidual)",L_L1NormmmsEnergyResidual_) ;
#line 187 "mms_setup.loci"
       input("L1Norm(mmsMassResidual),L1Norm(mmsMomentXResidual),  L1Norm(mmsMomentYResidual),L1Norm(mmsMomentZResidual),L1Norm(mmsEnergyResidual),  modelName") ;
#line 187 "mms_setup.loci"
       output("OUTPUT") ;
#line 187 "mms_setup.loci"
       conditional("do_plot") ;
#line 187 "mms_setup.loci"
    }
#line 187 "mms_setup.loci"
    void compute(const Loci::sequence &seq) { 
    if(Loci::MPI_rank != 0) return ;
    string fname = (*L_modelName_)+"_meanl1.dat" ;
    ofstream ofile(fname.c_str(),ios::out) ;
    ofile.precision(16) ;
    ofile << (*L_L1NormmmsMassResidual_)<< ' '
          << (*L_L1NormmmsMomentXResidual_)<< ' '
          << (*L_L1NormmmsMomentYResidual_)<< ' '
          << (*L_L1NormmmsMomentZResidual_)<< ' '
          << (*L_L1NormmmsEnergyResidual_)<< endl ;
  }} ;
#line 197 "mms_setup.loci"
Loci::register_rule<file_mms_setup017_1278521289m188> register_file_mms_setup017_1278521289m188 ;
#line 197 "mms_setup.loci"
}
#line 197 "mms_setup.loci"


  namespace {class file_mms_setup018_1278521289m189 : public Loci::singleton_rule {
#line 201 "mms_setup.loci"
    Loci::const_param<string>  L_modelName_ ; 
#line 201 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_LinfNormmmsMassResidual_ ; 
#line 201 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_LinfNormmmsMomentXResidual_ ; 
#line 201 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_LinfNormmmsMomentYResidual_ ; 
#line 201 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_LinfNormmmsMomentZResidual_ ; 
#line 201 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_LinfNormmmsEnergyResidual_ ; 
#line 201 "mms_setup.loci"
    Loci::param<bool> L_OUTPUT_ ; 
#line 201 "mms_setup.loci"
public:
#line 201 "mms_setup.loci"
    file_mms_setup018_1278521289m189() {
#line 201 "mms_setup.loci"
       name_store("OUTPUT",L_OUTPUT_) ;
#line 201 "mms_setup.loci"
       name_store("modelName",L_modelName_) ;
#line 201 "mms_setup.loci"
       name_store("LinfNorm(mmsMassResidual)",L_LinfNormmmsMassResidual_) ;
#line 201 "mms_setup.loci"
       name_store("LinfNorm(mmsMomentXResidual)",L_LinfNormmmsMomentXResidual_) ;
#line 201 "mms_setup.loci"
       name_store("LinfNorm(mmsMomentYResidual)",L_LinfNormmmsMomentYResidual_) ;
#line 201 "mms_setup.loci"
       name_store("LinfNorm(mmsMomentZResidual)",L_LinfNormmmsMomentZResidual_) ;
#line 201 "mms_setup.loci"
       name_store("LinfNorm(mmsEnergyResidual)",L_LinfNormmmsEnergyResidual_) ;
#line 201 "mms_setup.loci"
       input("LinfNorm(mmsMassResidual),LinfNorm(mmsMomentXResidual),  LinfNorm(mmsMomentYResidual),LinfNorm(mmsMomentZResidual),LinfNorm(mmsEnergyResidual),  modelName") ;
#line 201 "mms_setup.loci"
       output("OUTPUT") ;
#line 201 "mms_setup.loci"
       conditional("do_plot") ;
#line 201 "mms_setup.loci"
    }
#line 201 "mms_setup.loci"
    void compute(const Loci::sequence &seq) { 
    if(Loci::MPI_rank != 0) return ;

    string fname = (*L_modelName_)+"_meanlinf.dat" ;
    ofstream ofile(fname.c_str(),ios::out) ;
    ofile.precision(16) ;
    ofile << (*L_LinfNormmmsMassResidual_)<< ' '
          << (*L_LinfNormmmsMomentXResidual_)<< ' '
          << (*L_LinfNormmmsMomentYResidual_)<< ' '
          << (*L_LinfNormmmsMomentZResidual_)<< ' '
          << (*L_LinfNormmmsEnergyResidual_)<< endl ;
  }} ;
#line 212 "mms_setup.loci"
Loci::register_rule<file_mms_setup018_1278521289m189> register_file_mms_setup018_1278521289m189 ;
#line 212 "mms_setup.loci"
}
#line 212 "mms_setup.loci"


  namespace {class file_mms_setup019_1278521289m190 : public Loci::singleton_rule {
#line 215 "mms_setup.loci"
    Loci::const_param<string>  L_modelName_ ; 
#line 215 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_L2Normmms_rk_Residual_ ; 
#line 215 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_L2Normmms_rw_Residual_ ; 
#line 215 "mms_setup.loci"
    Loci::param<bool> L_OUTPUT_ ; 
#line 215 "mms_setup.loci"
public:
#line 215 "mms_setup.loci"
    file_mms_setup019_1278521289m190() {
#line 215 "mms_setup.loci"
       name_store("OUTPUT",L_OUTPUT_) ;
#line 215 "mms_setup.loci"
       name_store("modelName",L_modelName_) ;
#line 215 "mms_setup.loci"
       name_store("L2Norm(mms_rk_Residual)",L_L2Normmms_rk_Residual_) ;
#line 215 "mms_setup.loci"
       name_store("L2Norm(mms_rw_Residual)",L_L2Normmms_rw_Residual_) ;
#line 215 "mms_setup.loci"
       input("L2Norm(mms_rk_Residual),L2Norm(mms_rw_Residual),  modelName") ;
#line 215 "mms_setup.loci"
       output("OUTPUT") ;
#line 215 "mms_setup.loci"
       constraint("turbulentFlow,k,omega") ;
#line 215 "mms_setup.loci"
       conditional("do_plot") ;
#line 215 "mms_setup.loci"
    }
#line 215 "mms_setup.loci"
    void compute(const Loci::sequence &seq) { 
    if(Loci::MPI_rank != 0) return ;
    string fname = (*L_modelName_)+"_kwl2.dat" ;
    ofstream ofile(fname.c_str(),ios::out) ;
    ofile.precision(16) ;

    ofile << (*L_L2Normmms_rk_Residual_)<< ' '
          << (*L_L2Normmms_rw_Residual_)<< endl ;
  }} ;
#line 223 "mms_setup.loci"
Loci::register_rule<file_mms_setup019_1278521289m190> register_file_mms_setup019_1278521289m190 ;
#line 223 "mms_setup.loci"
}
#line 223 "mms_setup.loci"


  namespace {class file_mms_setup020_1278521289m190 : public Loci::singleton_rule {
#line 226 "mms_setup.loci"
    Loci::const_param<string>  L_modelName_ ; 
#line 226 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_L1Normmms_rk_Residual_ ; 
#line 226 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_L1Normmms_rw_Residual_ ; 
#line 226 "mms_setup.loci"
    Loci::param<bool> L_OUTPUT_ ; 
#line 226 "mms_setup.loci"
public:
#line 226 "mms_setup.loci"
    file_mms_setup020_1278521289m190() {
#line 226 "mms_setup.loci"
       name_store("OUTPUT",L_OUTPUT_) ;
#line 226 "mms_setup.loci"
       name_store("modelName",L_modelName_) ;
#line 226 "mms_setup.loci"
       name_store("L1Norm(mms_rk_Residual)",L_L1Normmms_rk_Residual_) ;
#line 226 "mms_setup.loci"
       name_store("L1Norm(mms_rw_Residual)",L_L1Normmms_rw_Residual_) ;
#line 226 "mms_setup.loci"
       input("L1Norm(mms_rk_Residual),L1Norm(mms_rw_Residual),  modelName") ;
#line 226 "mms_setup.loci"
       output("OUTPUT") ;
#line 226 "mms_setup.loci"
       constraint("turbulentFlow,k,omega") ;
#line 226 "mms_setup.loci"
       conditional("do_plot") ;
#line 226 "mms_setup.loci"
    }
#line 226 "mms_setup.loci"
    void compute(const Loci::sequence &seq) { 
    if(Loci::MPI_rank != 0) return ;
    string fname = (*L_modelName_)+"_kwl1.dat" ;
    ofstream ofile(fname.c_str(),ios::out) ;
    ofile.precision(16) ;

    ofile << (*L_L1Normmms_rk_Residual_)<< ' '
          << (*L_L1Normmms_rw_Residual_)<< endl ;
  }} ;
#line 234 "mms_setup.loci"
Loci::register_rule<file_mms_setup020_1278521289m190> register_file_mms_setup020_1278521289m190 ;
#line 234 "mms_setup.loci"
}
#line 234 "mms_setup.loci"


  namespace {class file_mms_setup021_1278521289m191 : public Loci::singleton_rule {
#line 237 "mms_setup.loci"
    Loci::const_param<string>  L_modelName_ ; 
#line 237 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_LinfNormmms_rk_Residual_ ; 
#line 237 "mms_setup.loci"
    Loci::const_param<Loci::real_t>  L_LinfNormmms_rw_Residual_ ; 
#line 237 "mms_setup.loci"
    Loci::param<bool> L_OUTPUT_ ; 
#line 237 "mms_setup.loci"
public:
#line 237 "mms_setup.loci"
    file_mms_setup021_1278521289m191() {
#line 237 "mms_setup.loci"
       name_store("OUTPUT",L_OUTPUT_) ;
#line 237 "mms_setup.loci"
       name_store("modelName",L_modelName_) ;
#line 237 "mms_setup.loci"
       name_store("LinfNorm(mms_rk_Residual)",L_LinfNormmms_rk_Residual_) ;
#line 237 "mms_setup.loci"
       name_store("LinfNorm(mms_rw_Residual)",L_LinfNormmms_rw_Residual_) ;
#line 237 "mms_setup.loci"
       input("LinfNorm(mms_rk_Residual),LinfNorm(mms_rw_Residual),  modelName") ;
#line 237 "mms_setup.loci"
       output("OUTPUT") ;
#line 237 "mms_setup.loci"
       constraint("turbulentFlow,k,omega") ;
#line 237 "mms_setup.loci"
       conditional("do_plot") ;
#line 237 "mms_setup.loci"
    }
#line 237 "mms_setup.loci"
    void compute(const Loci::sequence &seq) { 
    if(Loci::MPI_rank != 0) return ;
    string fname = (*L_modelName_)+"_kwlinf.dat" ;
    ofstream ofile(fname.c_str(),ios::out) ;
    ofile.precision(16) ;

    ofile << (*L_LinfNormmms_rk_Residual_)<< ' '
          << (*L_LinfNormmms_rw_Residual_)<< endl ;
  }} ;
#line 245 "mms_setup.loci"
Loci::register_rule<file_mms_setup021_1278521289m191> register_file_mms_setup021_1278521289m191 ;
#line 245 "mms_setup.loci"
}
#line 245 "mms_setup.loci"


  // $type X store<double> 
  // $type ToFloat(X0) store<float> 
  namespace {class file_mms_setup022_1278521289m192 : public Loci::pointwise_rule {
#line 249 "mms_setup.loci"
    Loci::const_store<double>  L_X_ ; 
#line 249 "mms_setup.loci"
    Loci::store<float>  L_ToFloatX_ ; 
#line 249 "mms_setup.loci"
public:
#line 249 "mms_setup.loci"
    file_mms_setup022_1278521289m192() {
#line 249 "mms_setup.loci"
       name_store("X",L_X_) ;
#line 249 "mms_setup.loci"
       name_store("ToFloat(X)",L_ToFloatX_) ;
#line 249 "mms_setup.loci"
       input("X") ;
#line 249 "mms_setup.loci"
       output("ToFloat(X)") ;
#line 249 "mms_setup.loci"
    }
#line 249 "mms_setup.loci"
    void calculate(Entity _e_) { 
#line 250 "mms_setup.loci"
    L_ToFloatX_[_e_]= L_X_[_e_];
  }    void compute(const Loci::sequence &seq) { 
#line 251 "mms_setup.loci"
      do_loop(seq,this) ;
#line 251 "mms_setup.loci"
    }
#line 251 "mms_setup.loci"
} ;
#line 251 "mms_setup.loci"
Loci::register_rule<file_mms_setup022_1278521289m192> register_file_mms_setup022_1278521289m192 ;
#line 251 "mms_setup.loci"
}
#line 251 "mms_setup.loci"


  // $type X store<vector3d<double> > 
  // $type ToFloatV(X0) store<vector3d<float> > 
  namespace {class file_mms_setup023_1278521289m192 : public Loci::pointwise_rule {
#line 255 "mms_setup.loci"
    Loci::const_store<vector3d<double> >  L_X_ ; 
#line 255 "mms_setup.loci"
    Loci::store<vector3d<float> >  L_ToFloatVX_ ; 
#line 255 "mms_setup.loci"
public:
#line 255 "mms_setup.loci"
    file_mms_setup023_1278521289m192() {
#line 255 "mms_setup.loci"
       name_store("X",L_X_) ;
#line 255 "mms_setup.loci"
       name_store("ToFloatV(X)",L_ToFloatVX_) ;
#line 255 "mms_setup.loci"
       input("X") ;
#line 255 "mms_setup.loci"
       output("ToFloatV(X)") ;
#line 255 "mms_setup.loci"
    }
#line 255 "mms_setup.loci"
    void calculate(Entity _e_) { 
#line 256 "mms_setup.loci"
    L_ToFloatVX_[_e_]= vector3d <float >(L_X_[_e_].x ,L_X_[_e_].y ,L_X_[_e_].z ) ;
  }    void compute(const Loci::sequence &seq) { 
#line 257 "mms_setup.loci"
      do_loop(seq,this) ;
#line 257 "mms_setup.loci"
    }
#line 257 "mms_setup.loci"
} ;
#line 257 "mms_setup.loci"
Loci::register_rule<file_mms_setup023_1278521289m192> register_file_mms_setup023_1278521289m192 ;
#line 257 "mms_setup.loci"
}
#line 257 "mms_setup.loci"


  // $type X storeVec<double> 
  // $type ToFloatVec(X0) store<float> 
  namespace {class file_mms_setup024_1278521289m193 : public Loci::pointwise_rule {
#line 261 "mms_setup.loci"
    Loci::const_storeVec<double>  L_X_ ; 
#line 261 "mms_setup.loci"
    Loci::store<float>  L_ToFloatVecX_ ; 
#line 261 "mms_setup.loci"
public:
#line 261 "mms_setup.loci"
    file_mms_setup024_1278521289m193() {
#line 261 "mms_setup.loci"
       name_store("X",L_X_) ;
#line 261 "mms_setup.loci"
       name_store("ToFloatVec(X)",L_ToFloatVecX_) ;
#line 261 "mms_setup.loci"
       input("X") ;
#line 261 "mms_setup.loci"
       output("ToFloatVec(X)") ;
#line 261 "mms_setup.loci"
    }
#line 261 "mms_setup.loci"
    void calculate(Entity _e_) { 
#line 262 "mms_setup.loci"
    int vs = L_X_.vecSize () ;
    double sum = 0 ;
    for (int i =0;i <vs ;++i )
      sum += L_X_[_e_][i ];
    L_ToFloatVecX_[_e_]= sum ;
  }    void compute(const Loci::sequence &seq) { 
#line 267 "mms_setup.loci"
      do_loop(seq,this) ;
#line 267 "mms_setup.loci"
    }
#line 267 "mms_setup.loci"
} ;
#line 267 "mms_setup.loci"
Loci::register_rule<file_mms_setup024_1278521289m193> register_file_mms_setup024_1278521289m193 ;
#line 267 "mms_setup.loci"
}
#line 267 "mms_setup.loci"

  
//OUTPUT_SCALAR("cell2nodeMaxMag(mmsMassResidual)",mmsMassError) ;
//OUTPUT_SCALAR("cell2nodeMaxMag(mmsMomentXResidual)",mmsMomentXError) ;
//OUTPUT_SCALAR("cell2nodeMaxMag(mmsMomentYResidual)",mmsMomentYError) ;
//OUTPUT_SCALAR("cell2nodeMaxMag(mmsMomentZResidual)",mmsMomentZError) ;
//OUTPUT_SCALAR("cell2nodeMaxMag(mmsEnergyResidual)",mmsEnergyError) ;

//OUTPUT_SCALAR("ToFloat(energyGeneralSrc(pos))",mmsEnergySrc) ;
//OUTPUT_SCALAR("ToFloatVec(massGeneralSrc(pos))",mmsMassSrc) ;
//OUTPUT_VECTOR("ToFloatV(momentumGeneralSrc(pos))",mmsMomentSrc) ;
  
}
