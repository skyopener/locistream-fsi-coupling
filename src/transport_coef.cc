#line 1 "transport_coef.loci"
#include <Loci.h> 
#include <Tools/tools.h>
#include "const.h"
#include "transport_db.h"
#include "turb_param.h"
#include "eos.h"
#include "property.h"
#include <Tools/parse.h>
using namespace Loci::parse ;
using namespace fluidPhysics ;
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
#line 11 "transport_coef.loci"


namespace streamUns {

  // $type transport_model param<std::string> 
  namespace {class file_transport_coef000_1280810697m356 : public Loci::default_rule {
#line 18 "transport_coef.loci"
    Loci::param<std::string>  L_transport_model_ ; 
#line 18 "transport_coef.loci"
public:
#line 18 "transport_coef.loci"
    file_transport_coef000_1280810697m356() {
#line 18 "transport_coef.loci"
       name_store("transport_model",L_transport_model_) ;
#line 18 "transport_coef.loci"
       output("transport_model") ;
#line 18 "transport_coef.loci"
       comments("Select the model used for transport of momentum and energy.  This can be 'none' for inviscid flow, 'sutherland' for sutherland's law (with default properties for air), 'powerLaw' for a powerlaw dependence on temperature, 'const_viscosity' for constant viscosity and conductivity, 'chemkin' for multicomponent flows, and 'database' for a species p,T database using Wilke's mixture rule.") ;
#line 18 "transport_coef.loci"
    }
#line 18 "transport_coef.loci"
    void compute(const Loci::sequence &seq) { 
    (*L_transport_model_)= "none" ;
  }} ;
#line 20 "transport_coef.loci"
Loci::register_rule<file_transport_coef000_1280810697m356> register_file_transport_coef000_1280810697m356 ;
#line 20 "transport_coef.loci"
}
#line 20 "transport_coef.loci"
 

  // $type const_vis Constraint
  // $type suther Constraint
  // $type powerLaw Constraint
  // $type chemk Constraint
  // $type database Constraint
  namespace {class file_transport_coef001_1280810697m357 : public Loci::constraint_rule {
#line 27 "transport_coef.loci"
    Loci::const_param<std::string>  L_transport_model_ ; 
#line 27 "transport_coef.loci"
    Loci::Constraint L_const_vis_ ; 
#line 27 "transport_coef.loci"
    Loci::Constraint L_suther_ ; 
#line 27 "transport_coef.loci"
    Loci::Constraint L_powerLaw_ ; 
#line 27 "transport_coef.loci"
    Loci::Constraint L_chemk_ ; 
#line 27 "transport_coef.loci"
    Loci::Constraint L_database_ ; 
#line 27 "transport_coef.loci"
public:
#line 27 "transport_coef.loci"
    file_transport_coef001_1280810697m357() {
#line 27 "transport_coef.loci"
       name_store("transport_model",L_transport_model_) ;
#line 27 "transport_coef.loci"
       name_store("const_vis",L_const_vis_) ;
#line 27 "transport_coef.loci"
       name_store("suther",L_suther_) ;
#line 27 "transport_coef.loci"
       name_store("powerLaw",L_powerLaw_) ;
#line 27 "transport_coef.loci"
       name_store("chemk",L_chemk_) ;
#line 27 "transport_coef.loci"
       name_store("database",L_database_) ;
#line 27 "transport_coef.loci"
       input("transport_model") ;
#line 27 "transport_coef.loci"
       output("const_vis") ;
#line 27 "transport_coef.loci"
       output("suther") ;
#line 27 "transport_coef.loci"
       output("powerLaw") ;
#line 27 "transport_coef.loci"
       output("chemk") ;
#line 27 "transport_coef.loci"
       output("database") ;
#line 27 "transport_coef.loci"
    }
#line 27 "transport_coef.loci"
    void compute(const Loci::sequence &seq) { 
      L_const_vis_= EMPTY ;
      L_suther_= EMPTY ;
      L_powerLaw_= EMPTY ;
      L_chemk_= EMPTY ;
      L_database_= EMPTY ;
      if((*L_transport_model_)== "const_viscosity") {
      L_const_vis_= ~EMPTY ;
      } else if((*L_transport_model_)== "sutherland") {
      L_suther_= ~EMPTY ;
      } else if((*L_transport_model_)== "powerLaw") {
      L_powerLaw_= ~EMPTY ;
      } else if((*L_transport_model_)== "chemkin") {
        L_chemk_= ~EMPTY ;
      } else if((*L_transport_model_)== "database") {
        L_database_= ~EMPTY ;
      } else if((*L_transport_model_)== "none") {
      } else {
        cerr << "Transport model " << (*L_transport_model_)<< " is not available" << endl ;
        cerr << "Available choices are:\n" ;
        cerr << "const_viscosity, sutherland, chemkin, powerLaw, database, and none\n" ;
        Loci::Abort() ;
      }
  }} ;
#line 50 "transport_coef.loci"
Loci::register_rule<file_transport_coef001_1280810697m357> register_file_transport_coef001_1280810697m357 ;
#line 50 "transport_coef.loci"
}
#line 50 "transport_coef.loci"
 

  // $type diffusion_model param<std::string> 
  namespace {class file_transport_coef002_1280810697m358 : public Loci::default_rule {
#line 55 "transport_coef.loci"
    Loci::param<std::string>  L_diffusion_model_ ; 
#line 55 "transport_coef.loci"
public:
#line 55 "transport_coef.loci"
    file_transport_coef002_1280810697m358() {
#line 55 "transport_coef.loci"
       name_store("diffusion_model",L_diffusion_model_) ;
#line 55 "transport_coef.loci"
       output("diffusion_model") ;
#line 55 "transport_coef.loci"
       comments("Mass diffusion model:  This may be 'laminarSchmidt' for a diffusivity based on viscosity and the specified Schmidt number, 'const_diffusion' for a constant diffusivity, and 'chemkin' for the CHEMKIN multicomponent diffusion properties. 'default' will select the diffusion model based on the transport model using 'chemkin' when chemkin is selected for transport model, or laminarSchmidt otherwise.") ;
#line 55 "transport_coef.loci"
    }
#line 55 "transport_coef.loci"
    void compute(const Loci::sequence &seq) { 
    (*L_diffusion_model_)= "default" ;
  }} ;
#line 57 "transport_coef.loci"
Loci::register_rule<file_transport_coef002_1280810697m358> register_file_transport_coef002_1280810697m358 ;
#line 57 "transport_coef.loci"
}
#line 57 "transport_coef.loci"
 

  // $type chemkin_diffusion Constraint
  // $type schmidt_diffusion Constraint
  // $type const_diffusion Constraint
  namespace {class file_transport_coef003_1280810697m359 : public Loci::constraint_rule {
#line 62 "transport_coef.loci"
    Loci::const_param<std::string>  L_transport_model_ ; 
#line 62 "transport_coef.loci"
    Loci::const_param<std::string>  L_diffusion_model_ ; 
#line 62 "transport_coef.loci"
    Loci::Constraint L_chemkin_diffusion_ ; 
#line 62 "transport_coef.loci"
    Loci::Constraint L_schmidt_diffusion_ ; 
#line 62 "transport_coef.loci"
    Loci::Constraint L_const_diffusion_ ; 
#line 62 "transport_coef.loci"
public:
#line 62 "transport_coef.loci"
    file_transport_coef003_1280810697m359() {
#line 62 "transport_coef.loci"
       name_store("transport_model",L_transport_model_) ;
#line 62 "transport_coef.loci"
       name_store("diffusion_model",L_diffusion_model_) ;
#line 62 "transport_coef.loci"
       name_store("chemkin_diffusion",L_chemkin_diffusion_) ;
#line 62 "transport_coef.loci"
       name_store("schmidt_diffusion",L_schmidt_diffusion_) ;
#line 62 "transport_coef.loci"
       name_store("const_diffusion",L_const_diffusion_) ;
#line 62 "transport_coef.loci"
       input("diffusion_model,transport_model") ;
#line 62 "transport_coef.loci"
       output("chemkin_diffusion") ;
#line 62 "transport_coef.loci"
       output("schmidt_diffusion") ;
#line 62 "transport_coef.loci"
       output("const_diffusion") ;
#line 62 "transport_coef.loci"
    }
#line 62 "transport_coef.loci"
    void compute(const Loci::sequence &seq) { 
      L_chemkin_diffusion_= EMPTY ;
      L_schmidt_diffusion_= EMPTY ;
      L_const_diffusion_= EMPTY ;
      if((*L_diffusion_model_)== "chemkin") {
        L_chemkin_diffusion_= ~EMPTY ;
      } else if((*L_diffusion_model_)== "laminarSchmidt") {
        L_schmidt_diffusion_= ~EMPTY ;
      } else if((*L_diffusion_model_)== "const_diffusivity") {
        L_const_diffusion_= ~EMPTY ;
      } else if((*L_diffusion_model_)== "default") {
        if((*L_transport_model_)== "chemkin")
          L_chemkin_diffusion_= ~EMPTY ;
        else if((*L_transport_model_)!= "none") 
          L_schmidt_diffusion_= ~EMPTY ;
      } else if((*L_diffusion_model_)== "none") {
      } else {
        cerr << "Diffusion model " << (*L_diffusion_model_)<< " is not available" << endl ;
        cerr << "Available choices are:\n" ;
        cerr << "chemkin, laminarSchmidt, const_diffusivity, and none\n" ;
        Loci::Abort() ;
      }
  }} ;
#line 84 "transport_coef.loci"
Loci::register_rule<file_transport_coef003_1280810697m359> register_file_transport_coef003_1280810697m359 ;
#line 84 "transport_coef.loci"
}
#line 84 "transport_coef.loci"
 
	
  // $type Sland param<Sutherland_param> 
  namespace {class file_transport_coef004_1280810697m360 : public Loci::default_rule {
#line 89 "transport_coef.loci"
    Loci::param<Sutherland_param>  L_Sland_ ; 
#line 89 "transport_coef.loci"
public:
#line 89 "transport_coef.loci"
    file_transport_coef004_1280810697m360() {
#line 89 "transport_coef.loci"
       name_store("Sland",L_Sland_) ;
#line 89 "transport_coef.loci"
       output("Sland") ;
#line 89 "transport_coef.loci"
       comments("Parameters used by sutherland model to compute conductivity and viscosity. a1-a3 define the viscosity, k1-k3 define the conductivity (if usepr == 0) otherwise conductivity is set by the Prandtl number, pr.") ;
#line 89 "transport_coef.loci"
    }
#line 89 "transport_coef.loci"
    void compute(const Loci::sequence &seq) { 
  }} ;
#line 90 "transport_coef.loci"
Loci::register_rule<file_transport_coef004_1280810697m360> register_file_transport_coef004_1280810697m360 ;
#line 90 "transport_coef.loci"
}
#line 90 "transport_coef.loci"
 

  // $type turbulent_transport param<turbulent_transport> 
  namespace {class file_transport_coef005_1280810697m360 : public Loci::default_rule {
#line 95 "transport_coef.loci"
    Loci::param<turbulent_transport>  L_turbulent_transport_ ; 
#line 95 "transport_coef.loci"
public:
#line 95 "transport_coef.loci"
    file_transport_coef005_1280810697m360() {
#line 95 "transport_coef.loci"
       name_store("turbulent_transport",L_turbulent_transport_) ;
#line 95 "transport_coef.loci"
       output("turbulent_transport") ;
#line 95 "transport_coef.loci"
       comments("Definition of the non-dimensional properties used to define turbulent conduction and diffusion.") ;
#line 95 "transport_coef.loci"
    }
#line 95 "transport_coef.loci"
    void compute(const Loci::sequence &seq) { 
  }} ;
#line 96 "transport_coef.loci"
Loci::register_rule<file_transport_coef005_1280810697m360> register_file_transport_coef005_1280810697m360 ;
#line 96 "transport_coef.loci"
}
#line 96 "transport_coef.loci"
 

  // $type powerLawParam param<powerLaw_param> 
  namespace {class file_transport_coef006_1280810697m361 : public Loci::optional_rule {
#line 101 "transport_coef.loci"
    Loci::param<powerLaw_param>  L_powerLawParam_ ; 
#line 101 "transport_coef.loci"
public:
#line 101 "transport_coef.loci"
    file_transport_coef006_1280810697m361() {
#line 101 "transport_coef.loci"
       name_store("powerLawParam",L_powerLawParam_) ;
#line 101 "transport_coef.loci"
       output("powerLawParam") ;
#line 101 "transport_coef.loci"
       comments("Parameters for power law definition of viscosity and conductivity") ;
#line 101 "transport_coef.loci"
    }
#line 101 "transport_coef.loci"
    void compute(const Loci::sequence &seq) { 
  }} ;
#line 102 "transport_coef.loci"
Loci::register_rule<file_transport_coef006_1280810697m361> register_file_transport_coef006_1280810697m361 ;
#line 102 "transport_coef.loci"
}
#line 102 "transport_coef.loci"
 

  // $type mu param<real> 
  namespace {class file_transport_coef007_1280810697m361 : public Loci::optional_rule {
#line 107 "transport_coef.loci"
    Loci::param<real>  L_mu_ ; 
#line 107 "transport_coef.loci"
public:
#line 107 "transport_coef.loci"
    file_transport_coef007_1280810697m361() {
#line 107 "transport_coef.loci"
       name_store("mu",L_mu_) ;
#line 107 "transport_coef.loci"
       output("mu") ;
#line 107 "transport_coef.loci"
       comments("Viscosity used in const_viscosity model.") ;
#line 107 "transport_coef.loci"
    }
#line 107 "transport_coef.loci"
    void compute(const Loci::sequence &seq) { 
  }} ;
#line 108 "transport_coef.loci"
Loci::register_rule<file_transport_coef007_1280810697m361> register_file_transport_coef007_1280810697m361 ;
#line 108 "transport_coef.loci"
}
#line 108 "transport_coef.loci"
 

  // $type kcond param<real> 

  namespace {class file_transport_coef008_1280810697m361 : public Loci::optional_rule {
#line 114 "transport_coef.loci"
    Loci::param<real>  L_kcond_ ; 
#line 114 "transport_coef.loci"
public:
#line 114 "transport_coef.loci"
    file_transport_coef008_1280810697m361() {
#line 114 "transport_coef.loci"
       name_store("kcond",L_kcond_) ;
#line 114 "transport_coef.loci"
       output("kcond") ;
#line 114 "transport_coef.loci"
       comments("conductivity used in const_viscosity transport model.") ;
#line 114 "transport_coef.loci"
    }
#line 114 "transport_coef.loci"
    void compute(const Loci::sequence &seq) { 
  }} ;
#line 115 "transport_coef.loci"
Loci::register_rule<file_transport_coef008_1280810697m361> register_file_transport_coef008_1280810697m361 ;
#line 115 "transport_coef.loci"
}
#line 115 "transport_coef.loci"
 

  // $type rhod param<real> 
  namespace {class file_transport_coef009_1280810697m362 : public Loci::optional_rule {
#line 120 "transport_coef.loci"
    Loci::param<real>  L_rhod_ ; 
#line 120 "transport_coef.loci"
public:
#line 120 "transport_coef.loci"
    file_transport_coef009_1280810697m362() {
#line 120 "transport_coef.loci"
       name_store("rhod",L_rhod_) ;
#line 120 "transport_coef.loci"
       output("rhod") ;
#line 120 "transport_coef.loci"
       comments("diffusion parameter for constant diffusivity mass diffusion model") ;
#line 120 "transport_coef.loci"
    }
#line 120 "transport_coef.loci"
    void compute(const Loci::sequence &seq) { 
  }} ;
#line 121 "transport_coef.loci"
Loci::register_rule<file_transport_coef009_1280810697m362> register_file_transport_coef009_1280810697m362 ;
#line 121 "transport_coef.loci"
}
#line 121 "transport_coef.loci"
 

  // $type laminarSchmidtNumber param<real> 
  namespace {class file_transport_coef010_1280810697m362 : public Loci::default_rule {
#line 126 "transport_coef.loci"
    Loci::param<real>  L_laminarSchmidtNumber_ ; 
#line 126 "transport_coef.loci"
public:
#line 126 "transport_coef.loci"
    file_transport_coef010_1280810697m362() {
#line 126 "transport_coef.loci"
       name_store("laminarSchmidtNumber",L_laminarSchmidtNumber_) ;
#line 126 "transport_coef.loci"
       output("laminarSchmidtNumber") ;
#line 126 "transport_coef.loci"
       comments("Laminar Schmidt number used in laminarSchmidt mass diffusion model.") ;
#line 126 "transport_coef.loci"
    }
#line 126 "transport_coef.loci"
    void compute(const Loci::sequence &seq) { 
    (*L_laminarSchmidtNumber_)= 0.9 ;
  }} ;
#line 128 "transport_coef.loci"
Loci::register_rule<file_transport_coef010_1280810697m362> register_file_transport_coef010_1280810697m362 ;
#line 128 "transport_coef.loci"
}
#line 128 "transport_coef.loci"
 
  
  namespace {class file_transport_coef011_1280810697m363 : public Loci::pointwise_rule {
#line 130 "transport_coef.loci"
    Loci::const_param<real>  L_mu_ ; 
#line 130 "transport_coef.loci"
    Loci::store<streamUns::real>  L_muuTPMIXTURE_ ; 
#line 130 "transport_coef.loci"
public:
#line 130 "transport_coef.loci"
    file_transport_coef011_1280810697m363() {
#line 130 "transport_coef.loci"
       name_store("muu(T,P,MIXTURE)",L_muuTPMIXTURE_) ;
#line 130 "transport_coef.loci"
       name_store("mu",L_mu_) ;
#line 130 "transport_coef.loci"
       input("mu") ;
#line 130 "transport_coef.loci"
       output("muu(T,P,MIXTURE)") ;
#line 130 "transport_coef.loci"
       constraint("const_vis,T,P,MIXTURE") ;
#line 130 "transport_coef.loci"
    }
#line 130 "transport_coef.loci"
    void calculate(Entity _e_) { 
#line 131 "transport_coef.loci"
    L_muuTPMIXTURE_[_e_]= L_mu_[_e_];
  }    void compute(const Loci::sequence &seq) { 
#line 132 "transport_coef.loci"
      do_loop(seq,this) ;
#line 132 "transport_coef.loci"
    }
#line 132 "transport_coef.loci"
} ;
#line 132 "transport_coef.loci"
Loci::register_rule<file_transport_coef011_1280810697m363> register_file_transport_coef011_1280810697m363 ;
#line 132 "transport_coef.loci"
}
#line 132 "transport_coef.loci"
 

  //specify conductivity as a constant
  namespace {class file_transport_coef012_1280810697m363 : public Loci::pointwise_rule {
#line 135 "transport_coef.loci"
    Loci::const_param<real>  L_kcond_ ; 
#line 135 "transport_coef.loci"
    Loci::store<streamUns::real>  L_kconductTPMIXTURE_ ; 
#line 135 "transport_coef.loci"
public:
#line 135 "transport_coef.loci"
    file_transport_coef012_1280810697m363() {
#line 135 "transport_coef.loci"
       name_store("kconduct(T,P,MIXTURE)",L_kconductTPMIXTURE_) ;
#line 135 "transport_coef.loci"
       name_store("kcond",L_kcond_) ;
#line 135 "transport_coef.loci"
       input("kcond") ;
#line 135 "transport_coef.loci"
       output("kconduct(T,P,MIXTURE)") ;
#line 135 "transport_coef.loci"
       constraint("const_vis,T,P,MIXTURE") ;
#line 135 "transport_coef.loci"
    }
#line 135 "transport_coef.loci"
    void calculate(Entity _e_) { 
#line 136 "transport_coef.loci"
   L_kconductTPMIXTURE_[_e_]= L_kcond_[_e_];
  }    void compute(const Loci::sequence &seq) { 
#line 137 "transport_coef.loci"
      do_loop(seq,this) ;
#line 137 "transport_coef.loci"
    }
#line 137 "transport_coef.loci"
} ;
#line 137 "transport_coef.loci"
Loci::register_rule<file_transport_coef012_1280810697m363> register_file_transport_coef012_1280810697m363 ;
#line 137 "transport_coef.loci"
}
#line 137 "transport_coef.loci"
 

  //calculate viscosity by Sutherland's Law  
  // $type T store<real> 
  namespace {class file_transport_coef013_1280810697m364 : public Loci::pointwise_rule {
#line 141 "transport_coef.loci"
    Loci::const_store<real>  L_T_ ; 
#line 141 "transport_coef.loci"
    Loci::const_param<Sutherland_param>  L_Sland_ ; 
#line 141 "transport_coef.loci"
    Loci::store<streamUns::real>  L_muuTPMIXTURE_ ; 
#line 141 "transport_coef.loci"
public:
#line 141 "transport_coef.loci"
    file_transport_coef013_1280810697m364() {
#line 141 "transport_coef.loci"
       name_store("T",L_T_) ;
#line 141 "transport_coef.loci"
       name_store("muu(T,P,MIXTURE)",L_muuTPMIXTURE_) ;
#line 141 "transport_coef.loci"
       name_store("Sland",L_Sland_) ;
#line 141 "transport_coef.loci"
       input("T,Sland") ;
#line 141 "transport_coef.loci"
       output("muu(T,P,MIXTURE)") ;
#line 141 "transport_coef.loci"
       constraint("suther,T") ;
#line 141 "transport_coef.loci"
    }
#line 141 "transport_coef.loci"
    void calculate(Entity _e_) { 
#line 142 "transport_coef.loci"
      //apply Sutherland's law (refer to "Computational Fluid
      //Mechanics and Heat Transfer", Dale Anderson, pp259)
      real a1 =L_Sland_[_e_].a1 ;
      real a2 =L_Sland_[_e_].a2 ;
      real a3 =L_Sland_[_e_].a3 ;
      L_muuTPMIXTURE_[_e_]= a1 *pow (L_T_[_e_],a2 )/(L_T_[_e_]+a3 ) ;
  }    void compute(const Loci::sequence &seq) { 
#line 148 "transport_coef.loci"
      do_loop(seq,this) ;
#line 148 "transport_coef.loci"
    }
#line 148 "transport_coef.loci"
} ;
#line 148 "transport_coef.loci"
Loci::register_rule<file_transport_coef013_1280810697m364> register_file_transport_coef013_1280810697m364 ;
#line 148 "transport_coef.loci"
}
#line 148 "transport_coef.loci"
 

  //calculate viscosity by using power Law.
  namespace {class file_transport_coef014_1280810697m365 : public Loci::pointwise_rule {
#line 151 "transport_coef.loci"
    Loci::const_store<real>  L_T_ ; 
#line 151 "transport_coef.loci"
    Loci::const_param<powerLaw_param>  L_powerLawParam_ ; 
#line 151 "transport_coef.loci"
    Loci::store<streamUns::real>  L_muuTPMIXTURE_ ; 
#line 151 "transport_coef.loci"
public:
#line 151 "transport_coef.loci"
    file_transport_coef014_1280810697m365() {
#line 151 "transport_coef.loci"
       name_store("T",L_T_) ;
#line 151 "transport_coef.loci"
       name_store("muu(T,P,MIXTURE)",L_muuTPMIXTURE_) ;
#line 151 "transport_coef.loci"
       name_store("powerLawParam",L_powerLawParam_) ;
#line 151 "transport_coef.loci"
       input("T,powerLawParam") ;
#line 151 "transport_coef.loci"
       output("muu(T,P,MIXTURE)") ;
#line 151 "transport_coef.loci"
       constraint("powerLaw,T") ;
#line 151 "transport_coef.loci"
    }
#line 151 "transport_coef.loci"
    void calculate(Entity _e_) { 
#line 152 "transport_coef.loci"
      const real T_ref = L_powerLawParam_[_e_].T_ref ;
      const real mu_ref = L_powerLawParam_[_e_].mu_ref ;
      const real power = L_powerLawParam_[_e_].power ;
      L_muuTPMIXTURE_[_e_]= mu_ref *pow (L_T_[_e_]/T_ref ,power ) ;
  }    void compute(const Loci::sequence &seq) { 
#line 156 "transport_coef.loci"
      do_loop(seq,this) ;
#line 156 "transport_coef.loci"
    }
#line 156 "transport_coef.loci"
} ;
#line 156 "transport_coef.loci"
Loci::register_rule<file_transport_coef014_1280810697m365> register_file_transport_coef014_1280810697m365 ;
#line 156 "transport_coef.loci"
}
#line 156 "transport_coef.loci"
 

  //calculate viscosity and conductivity coefficients 
  //using CHEMKIN curve fit data. 
  // $type MIXTURE storeVec<real> 
  // $type tran param<transport_db> 
  namespace {class file_transport_coef015_1280810697m366 : public Loci::pointwise_rule {
#line 163 "transport_coef.loci"
    Loci::const_param<transport_db>  L_tran_n__ ; 
#line 163 "transport_coef.loci"
    Loci::const_param<int>  L_numSpecies_n__ ; 
#line 163 "transport_coef.loci"
    Loci::const_param<std::vector<streamUns::real> >  L_speciesMolecularMass_n__ ; 
#line 163 "transport_coef.loci"
    Loci::const_storeVec<real>  L_MIXTURE_n__ ; 
#line 163 "transport_coef.loci"
    Loci::const_store<real>  L_T_n__ ; 
#line 163 "transport_coef.loci"
    Loci::store<streamUns::real>  L_kconductTPMIXTURE_n__ ; 
#line 163 "transport_coef.loci"
public:
#line 163 "transport_coef.loci"
    file_transport_coef015_1280810697m366() {
#line 163 "transport_coef.loci"
       name_store("tran{n}",L_tran_n__) ;
#line 163 "transport_coef.loci"
       name_store("numSpecies{n}",L_numSpecies_n__) ;
#line 163 "transport_coef.loci"
       name_store("speciesMolecularMass{n}",L_speciesMolecularMass_n__) ;
#line 163 "transport_coef.loci"
       name_store("MIXTURE{n}",L_MIXTURE_n__) ;
#line 163 "transport_coef.loci"
       name_store("T{n}",L_T_n__) ;
#line 163 "transport_coef.loci"
       name_store("kconduct(T,P,MIXTURE){n}",L_kconductTPMIXTURE_n__) ;
#line 163 "transport_coef.loci"
       input("tran{n},numSpecies{n},  speciesMolecularMass{n},MIXTURE{n},T{n}") ;
#line 163 "transport_coef.loci"
       output("kconduct(T,P,MIXTURE){n}") ;
#line 163 "transport_coef.loci"
       constraint("chemk{n},T{n}") ;
#line 163 "transport_coef.loci"
    }
#line 163 "transport_coef.loci"
    void calculate(Entity _e_) { 
#line 164 "transport_coef.loci"

      int ns =L_numSpecies_n__[_e_];
      tmp_array <real > mf (ns ) ;

      if (ns > 1) {
        // Convert mass fractions to mole fractions
        real sum = 0 ;
        for (int i =0;i <ns ;++i ) {
          mf [i ] = L_MIXTURE_n__[_e_][i ]/L_speciesMolecularMass_n__[_e_][i ]; 
          sum += mf [i ] ;
        }
        real rsum = 1./sum ;
        for (int i =0;i <ns ;++i ) 
          mf [i ] *= rsum ;
      } else mf [0] = 1.0 ;
      L_kconductTPMIXTURE_n__[_e_]= L_tran_n__[_e_].mcacon (L_T_n__[_e_],mf ) ;
  }    void compute(const Loci::sequence &seq) { 
#line 180 "transport_coef.loci"
      do_loop(seq,this) ;
#line 180 "transport_coef.loci"
    }
#line 180 "transport_coef.loci"
} ;
#line 180 "transport_coef.loci"
Loci::register_rule<file_transport_coef015_1280810697m366> register_file_transport_coef015_1280810697m366 ;
#line 180 "transport_coef.loci"
}
#line 180 "transport_coef.loci"
 

  // Force this rule to evaluate only once per timestep. JW 11/24/2008.
  namespace {class file_transport_coef016_1280810697m367 : public Loci::pointwise_rule {
#line 184 "transport_coef.loci"
    Loci::const_param<transport_db>  L_tran_n__ ; 
#line 184 "transport_coef.loci"
    Loci::const_param<int>  L_numSpecies_n__ ; 
#line 184 "transport_coef.loci"
    Loci::const_param<std::vector<streamUns::real> >  L_speciesMolecularMass_n__ ; 
#line 184 "transport_coef.loci"
    Loci::const_storeVec<real>  L_MIXTURE_n__ ; 
#line 184 "transport_coef.loci"
    Loci::const_store<real>  L_T_n__ ; 
#line 184 "transport_coef.loci"
    Loci::store<streamUns::real>  L_muuTPMIXTURE_n__ ; 
#line 184 "transport_coef.loci"
public:
#line 184 "transport_coef.loci"
    file_transport_coef016_1280810697m367() {
#line 184 "transport_coef.loci"
       name_store("tran{n}",L_tran_n__) ;
#line 184 "transport_coef.loci"
       name_store("numSpecies{n}",L_numSpecies_n__) ;
#line 184 "transport_coef.loci"
       name_store("speciesMolecularMass{n}",L_speciesMolecularMass_n__) ;
#line 184 "transport_coef.loci"
       name_store("MIXTURE{n}",L_MIXTURE_n__) ;
#line 184 "transport_coef.loci"
       name_store("T{n}",L_T_n__) ;
#line 184 "transport_coef.loci"
       name_store("muu(T,P,MIXTURE){n}",L_muuTPMIXTURE_n__) ;
#line 184 "transport_coef.loci"
       input("tran{n},speciesMolecularMass{n},  numSpecies{n},MIXTURE{n},T{n}") ;
#line 184 "transport_coef.loci"
       output("muu(T,P,MIXTURE){n}") ;
#line 184 "transport_coef.loci"
       constraint("chemk{n},T{n}") ;
#line 184 "transport_coef.loci"
    }
#line 184 "transport_coef.loci"
    void calculate(Entity _e_) { 
#line 185 "transport_coef.loci"
      int ns =L_numSpecies_n__[_e_];
      tmp_array <real > mf (ns ) ;

      if (ns > 1) {
        // Convert mass fractions to mole fractions
        real sum = 0 ;
        for (int i =0;i <ns ;++i ) {
          mf [i ] = L_MIXTURE_n__[_e_][i ]/L_speciesMolecularMass_n__[_e_][i ]; 
          sum += mf [i ] ;
        }
        real rsum = 1./sum ;
        for (int i =0;i <ns ;++i ) 
          mf [i ] *= rsum ;
      } else mf [0] = 1.0 ;
      L_muuTPMIXTURE_n__[_e_]= L_tran_n__[_e_].mcavis (L_T_n__[_e_],mf ) ; 
  }    void compute(const Loci::sequence &seq) { 
#line 200 "transport_coef.loci"
      do_loop(seq,this) ;
#line 200 "transport_coef.loci"
    }
#line 200 "transport_coef.loci"
} ;
#line 200 "transport_coef.loci"
Loci::register_rule<file_transport_coef016_1280810697m367> register_file_transport_coef016_1280810697m367 ;
#line 200 "transport_coef.loci"
}
#line 200 "transport_coef.loci"


  namespace {class file_transport_coef017_1280810697m369 : public Loci::pointwise_rule {
#line 203 "transport_coef.loci"
    Loci::const_store<streamUns::real>  L_cp_ ; 
#line 203 "transport_coef.loci"
    Loci::const_store<real>  L_T_ ; 
#line 203 "transport_coef.loci"
    Loci::const_store<streamUns::real>  L_muuTPMIXTURE_ ; 
#line 203 "transport_coef.loci"
    Loci::const_param<Sutherland_param>  L_Sland_ ; 
#line 203 "transport_coef.loci"
    Loci::store<streamUns::real>  L_kconductTPMIXTURE_ ; 
#line 203 "transport_coef.loci"
public:
#line 203 "transport_coef.loci"
    file_transport_coef017_1280810697m369() {
#line 203 "transport_coef.loci"
       name_store("cp",L_cp_) ;
#line 203 "transport_coef.loci"
       name_store("T",L_T_) ;
#line 203 "transport_coef.loci"
       name_store("muu(T,P,MIXTURE)",L_muuTPMIXTURE_) ;
#line 203 "transport_coef.loci"
       name_store("kconduct(T,P,MIXTURE)",L_kconductTPMIXTURE_) ;
#line 203 "transport_coef.loci"
       name_store("Sland",L_Sland_) ;
#line 203 "transport_coef.loci"
       input("muu(T,P,MIXTURE),cp,Sland,T") ;
#line 203 "transport_coef.loci"
       output("kconduct(T,P,MIXTURE)") ;
#line 203 "transport_coef.loci"
       constraint("suther,T") ;
#line 203 "transport_coef.loci"
    }
#line 203 "transport_coef.loci"
    void calculate(Entity _e_) { 
#line 204 "transport_coef.loci"
      //apply Sutherland's law for air (refer to "Computational Fluid
      //Mechanics and Heat Transfer", Dale Anderson, pp259)
      //kconduct[fc] = 2.495e-3*pow(T[fc],1.5)/(T[fc]
      //							  +194.0) ;
      //get conductivity by the definition of Prantle number
      if (L_Sland_[_e_].usepr )
        L_kconductTPMIXTURE_[_e_]= L_muuTPMIXTURE_[_e_]*L_cp_[_e_]/L_Sland_[_e_].pr ;
      else L_kconductTPMIXTURE_[_e_]= L_Sland_[_e_].k1 *pow (L_T_[_e_],L_Sland_[_e_].k2 )/(L_T_[_e_]+L_Sland_[_e_].k3 ) ;
  }    void compute(const Loci::sequence &seq) { 
#line 213 "transport_coef.loci"
      do_loop(seq,this) ;
#line 213 "transport_coef.loci"
    }
#line 213 "transport_coef.loci"
} ;
#line 213 "transport_coef.loci"
Loci::register_rule<file_transport_coef017_1280810697m369> register_file_transport_coef017_1280810697m369 ;
#line 213 "transport_coef.loci"
}
#line 213 "transport_coef.loci"
 


  namespace {class file_transport_coef018_1280810697m370 : public Loci::pointwise_rule {
#line 217 "transport_coef.loci"
    Loci::const_store<streamUns::real>  L_cp_ ; 
#line 217 "transport_coef.loci"
    Loci::const_store<real>  L_T_ ; 
#line 217 "transport_coef.loci"
    Loci::const_store<streamUns::real>  L_muuTPMIXTURE_ ; 
#line 217 "transport_coef.loci"
    Loci::const_param<powerLaw_param>  L_powerLawParam_ ; 
#line 217 "transport_coef.loci"
    Loci::store<streamUns::real>  L_kconductTPMIXTURE_ ; 
#line 217 "transport_coef.loci"
public:
#line 217 "transport_coef.loci"
    file_transport_coef018_1280810697m370() {
#line 217 "transport_coef.loci"
       name_store("cp",L_cp_) ;
#line 217 "transport_coef.loci"
       name_store("T",L_T_) ;
#line 217 "transport_coef.loci"
       name_store("muu(T,P,MIXTURE)",L_muuTPMIXTURE_) ;
#line 217 "transport_coef.loci"
       name_store("kconduct(T,P,MIXTURE)",L_kconductTPMIXTURE_) ;
#line 217 "transport_coef.loci"
       name_store("powerLawParam",L_powerLawParam_) ;
#line 217 "transport_coef.loci"
       input("muu(T,P,MIXTURE),cp,powerLawParam,T") ;
#line 217 "transport_coef.loci"
       output("kconduct(T,P,MIXTURE)") ;
#line 217 "transport_coef.loci"
       constraint("powerLaw,T") ;
#line 217 "transport_coef.loci"
    }
#line 217 "transport_coef.loci"
    void calculate(Entity _e_) { 
#line 218 "transport_coef.loci"
      const real Pr = L_powerLawParam_[_e_].Pr ;
      L_kconductTPMIXTURE_[_e_]= L_muuTPMIXTURE_[_e_]*L_cp_[_e_]/Pr ;
  }    void compute(const Loci::sequence &seq) { 
#line 220 "transport_coef.loci"
      do_loop(seq,this) ;
#line 220 "transport_coef.loci"
    }
#line 220 "transport_coef.loci"
} ;
#line 220 "transport_coef.loci"
Loci::register_rule<file_transport_coef018_1280810697m370> register_file_transport_coef018_1280810697m370 ;
#line 220 "transport_coef.loci"
}
#line 220 "transport_coef.loci"
 

 
  namespace {class file_transport_coef019_1280810697m371 : public Loci::pointwise_rule {
#line 224 "transport_coef.loci"
    Loci::store<streamUns::real>  L_kconducttemperature_fgaugePressure_fmixture_f_ ; 
#line 224 "transport_coef.loci"
public:
#line 224 "transport_coef.loci"
    file_transport_coef019_1280810697m371() {
#line 224 "transport_coef.loci"
       name_store("kconduct(temperature_f,gaugePressure_f,mixture_f)",L_kconducttemperature_fgaugePressure_fmixture_f_) ;
#line 224 "transport_coef.loci"
       input("kconduct(temperature_f,gaugePressure_f,mixture_f)") ;
#line 224 "transport_coef.loci"
       output("kconduct_f=kconduct(temperature_f,gaugePressure_f,mixture_f)") ;
#line 224 "transport_coef.loci"
    }
#line 224 "transport_coef.loci"
    void calculate(Entity _e_) { 
#line 225 "transport_coef.loci"
  }    void compute(const Loci::sequence &seq) { 
#line 225 "transport_coef.loci"
      do_loop(seq,this) ;
#line 225 "transport_coef.loci"
    }
#line 225 "transport_coef.loci"
} ;
#line 225 "transport_coef.loci"
Loci::register_rule<file_transport_coef019_1280810697m371> register_file_transport_coef019_1280810697m371 ;
#line 225 "transport_coef.loci"
}
#line 225 "transport_coef.loci"
 


  namespace {class file_transport_coef020_1280810697m371 : public Loci::pointwise_rule {
#line 229 "transport_coef.loci"
    Loci::store<streamUns::real>  L_muutemperature_fgaugePressure_fmixture_f_ ; 
#line 229 "transport_coef.loci"
public:
#line 229 "transport_coef.loci"
    file_transport_coef020_1280810697m371() {
#line 229 "transport_coef.loci"
       name_store("muu(temperature_f,gaugePressure_f,mixture_f)",L_muutemperature_fgaugePressure_fmixture_f_) ;
#line 229 "transport_coef.loci"
       input("muu(temperature_f,gaugePressure_f,mixture_f)") ;
#line 229 "transport_coef.loci"
       output("muu_f=muu(temperature_f,gaugePressure_f,mixture_f)") ;
#line 229 "transport_coef.loci"
    }
#line 229 "transport_coef.loci"
    void calculate(Entity _e_) { 
#line 230 "transport_coef.loci"
  }    void compute(const Loci::sequence &seq) { 
#line 230 "transport_coef.loci"
      do_loop(seq,this) ;
#line 230 "transport_coef.loci"
    }
#line 230 "transport_coef.loci"
} ;
#line 230 "transport_coef.loci"
Loci::register_rule<file_transport_coef020_1280810697m371> register_file_transport_coef020_1280810697m371 ;
#line 230 "transport_coef.loci"
}
#line 230 "transport_coef.loci"
 


  namespace {class file_transport_coef021_1280810697m372 : public Loci::pointwise_rule {
#line 234 "transport_coef.loci"
    Loci::store<streamUns::real>  L_kconducttemperaturegaugePressuremixture_ ; 
#line 234 "transport_coef.loci"
public:
#line 234 "transport_coef.loci"
    file_transport_coef021_1280810697m372() {
#line 234 "transport_coef.loci"
       name_store("kconduct(temperature,gaugePressure,mixture)",L_kconducttemperaturegaugePressuremixture_) ;
#line 234 "transport_coef.loci"
       input("kconduct(temperature,gaugePressure,mixture)") ;
#line 234 "transport_coef.loci"
       output("kconduct=kconduct(temperature,gaugePressure,mixture)") ;
#line 234 "transport_coef.loci"
    }
#line 234 "transport_coef.loci"
    void calculate(Entity _e_) { 
#line 235 "transport_coef.loci"
  }    void compute(const Loci::sequence &seq) { 
#line 235 "transport_coef.loci"
      do_loop(seq,this) ;
#line 235 "transport_coef.loci"
    }
#line 235 "transport_coef.loci"
} ;
#line 235 "transport_coef.loci"
Loci::register_rule<file_transport_coef021_1280810697m372> register_file_transport_coef021_1280810697m372 ;
#line 235 "transport_coef.loci"
}
#line 235 "transport_coef.loci"
 

  namespace {class file_transport_coef022_1280810697m372 : public Loci::pointwise_rule {
#line 238 "transport_coef.loci"
    Loci::store<streamUns::real>  L_muutemperaturegaugePressuremixture_ ; 
#line 238 "transport_coef.loci"
public:
#line 238 "transport_coef.loci"
    file_transport_coef022_1280810697m372() {
#line 238 "transport_coef.loci"
       name_store("muu(temperature,gaugePressure,mixture)",L_muutemperaturegaugePressuremixture_) ;
#line 238 "transport_coef.loci"
       input("muu(temperature,gaugePressure,mixture)") ;
#line 238 "transport_coef.loci"
       output("muu=muu(temperature,gaugePressure,mixture)") ;
#line 238 "transport_coef.loci"
    }
#line 238 "transport_coef.loci"
    void calculate(Entity _e_) { 
#line 239 "transport_coef.loci"
  }    void compute(const Loci::sequence &seq) { 
#line 239 "transport_coef.loci"
      do_loop(seq,this) ;
#line 239 "transport_coef.loci"
    }
#line 239 "transport_coef.loci"
} ;
#line 239 "transport_coef.loci"
Loci::register_rule<file_transport_coef022_1280810697m372> register_file_transport_coef022_1280810697m372 ;
#line 239 "transport_coef.loci"
}
#line 239 "transport_coef.loci"
 
 
// compute diffusion coefficients using CHEMKIN library. 
  namespace {class file_transport_coef023_1280810697m373 : public Loci::pointwise_rule {
#line 244 "transport_coef.loci"
    Loci::const_param<streamUns::real>  L_Pambient_ ; 
#line 244 "transport_coef.loci"
    Loci::const_store<streamUns::real>  L_gaugePressure_f_ ; 
#line 244 "transport_coef.loci"
    Loci::const_storeVec<streamUns::real>  L_mixture_f_ ; 
#line 244 "transport_coef.loci"
    Loci::const_store<streamUns::real>  L_temperature_f_ ; 
#line 244 "transport_coef.loci"
    Loci::const_store<streamUns::real>  L_rho_f_ ; 
#line 244 "transport_coef.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 244 "transport_coef.loci"
    Loci::const_store<streamUns::real>  L_tmuu_f_ ; 
#line 244 "transport_coef.loci"
    Loci::const_param<turbulent_transport>  L_turbulent_transport_ ; 
#line 244 "transport_coef.loci"
    Loci::const_param<transport_db>  L_tran_ ; 
#line 244 "transport_coef.loci"
    Loci::storeVec<streamUns::real>  L_diff_ ; 
#line 244 "transport_coef.loci"
public:
#line 244 "transport_coef.loci"
    file_transport_coef023_1280810697m373() {
#line 244 "transport_coef.loci"
       name_store("Pambient",L_Pambient_) ;
#line 244 "transport_coef.loci"
       name_store("gaugePressure_f",L_gaugePressure_f_) ;
#line 244 "transport_coef.loci"
       name_store("mixture_f",L_mixture_f_) ;
#line 244 "transport_coef.loci"
       name_store("temperature_f",L_temperature_f_) ;
#line 244 "transport_coef.loci"
       name_store("rho_f",L_rho_f_) ;
#line 244 "transport_coef.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 244 "transport_coef.loci"
       name_store("diff",L_diff_) ;
#line 244 "transport_coef.loci"
       name_store("tmuu_f",L_tmuu_f_) ;
#line 244 "transport_coef.loci"
       name_store("turbulent_transport",L_turbulent_transport_) ;
#line 244 "transport_coef.loci"
       name_store("tran",L_tran_) ;
#line 244 "transport_coef.loci"
       input("numSpecies,tran,turbulent_transport,Pambient,temperature_f,mixture_f,gaugePressure_f,rho_f,tmuu_f") ;
#line 244 "transport_coef.loci"
       output("diff") ;
#line 244 "transport_coef.loci"
       constraint("chemkin_diffusion,(cr,cl)->(vol,mixture)") ;
#line 244 "transport_coef.loci"
    }
#line 244 "transport_coef.loci"
    void prelude(const Loci::sequence &seq) { 
    L_diff_.setVecSize(*L_numSpecies_) ;
  }    void calculate(Entity _e_) { 
#line 247 "transport_coef.loci"
    const int ns = L_numSpecies_[_e_];
    if (ns ==1) return ;
    
    tmp_array <real > mf (ns ) ;
    
    const real sct = L_turbulent_transport_[_e_].TurbulentSchmidtNumber ;
    real mutsc = L_tmuu_f_[_e_]/sct ;
    
    real tface = L_temperature_f_[_e_];
    real pface = L_gaugePressure_f_[_e_]+ L_Pambient_[_e_];
    
        real sum = 0.0 ;
        for (int i =0;i <ns ;++i ) {
          mf [i ] = L_mixture_f_[_e_][i ];
          sum += mf [i ]/L_tran_[_e_].MolecularWeight (i ) ;
        }

        for (int i =0;i <ns ;++i )
          mf [i ] = mf [i ]/(L_tran_[_e_].MolecularWeight (i )*sum ) ;
        L_tran_[_e_].mcadif (pface ,tface ,mf ,L_diff_[_e_]) ;

        const real rrf = 1./L_rho_f_[_e_];
        for (int i =0;i <ns ;++i )
          L_diff_[_e_][i ]+= mutsc *rrf ;
    }    void compute(const Loci::sequence &seq) { 
#line 271 "transport_coef.loci"
      prelude(seq) ;
#line 271 "transport_coef.loci"
      do_loop(seq,this) ;
#line 271 "transport_coef.loci"
    }
#line 271 "transport_coef.loci"
} ;
#line 271 "transport_coef.loci"
Loci::register_rule<file_transport_coef023_1280810697m373> register_file_transport_coef023_1280810697m373 ;
#line 271 "transport_coef.loci"
}
#line 271 "transport_coef.loci"


// compute diffusion coefficients on boundary faces using CHEMKIN library. 
// Note: diffusion coefficient was not set to zero in old code
  namespace {class file_transport_coef024_1280810697m374 : public Loci::pointwise_rule {
#line 276 "transport_coef.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 276 "transport_coef.loci"
    Loci::storeVec<streamUns::real>  L_diff_ ; 
#line 276 "transport_coef.loci"
public:
#line 276 "transport_coef.loci"
    file_transport_coef024_1280810697m374() {
#line 276 "transport_coef.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 276 "transport_coef.loci"
       name_store("diff",L_diff_) ;
#line 276 "transport_coef.loci"
       input("numSpecies") ;
#line 276 "transport_coef.loci"
       output("diff") ;
#line 276 "transport_coef.loci"
       constraint("ci->vol,chemkin_diffusion") ;
#line 276 "transport_coef.loci"
    }
#line 276 "transport_coef.loci"
    void prelude(const Loci::sequence &seq) { 
    L_diff_.setVecSize(*L_numSpecies_) ;
	}    void calculate(Entity _e_) { 
#line 279 "transport_coef.loci"
    const int ns = L_numSpecies_[_e_];
      for (int i =0; i <ns ; ++i )
          L_diff_[_e_][i ]= 0 ;
  }    void compute(const Loci::sequence &seq) { 
#line 282 "transport_coef.loci"
      prelude(seq) ;
#line 282 "transport_coef.loci"
      do_loop(seq,this) ;
#line 282 "transport_coef.loci"
    }
#line 282 "transport_coef.loci"
} ;
#line 282 "transport_coef.loci"
Loci::register_rule<file_transport_coef024_1280810697m374> register_file_transport_coef024_1280810697m374 ;
#line 282 "transport_coef.loci"
}
#line 282 "transport_coef.loci"


  namespace {class file_transport_coef025_1280810697m375 : public Loci::pointwise_rule {
#line 286 "transport_coef.loci"
    Loci::const_store<streamUns::real>  L_rho_f_ ; 
#line 286 "transport_coef.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 286 "transport_coef.loci"
    Loci::const_store<streamUns::real>  L_tmuu_f_ ; 
#line 286 "transport_coef.loci"
    Loci::const_param<turbulent_transport>  L_turbulent_transport_ ; 
#line 286 "transport_coef.loci"
    Loci::const_param<real>  L_rhod_ ; 
#line 286 "transport_coef.loci"
    Loci::storeVec<streamUns::real>  L_diff_ ; 
#line 286 "transport_coef.loci"
public:
#line 286 "transport_coef.loci"
    file_transport_coef025_1280810697m375() {
#line 286 "transport_coef.loci"
       name_store("rho_f",L_rho_f_) ;
#line 286 "transport_coef.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 286 "transport_coef.loci"
       name_store("diff",L_diff_) ;
#line 286 "transport_coef.loci"
       name_store("tmuu_f",L_tmuu_f_) ;
#line 286 "transport_coef.loci"
       name_store("turbulent_transport",L_turbulent_transport_) ;
#line 286 "transport_coef.loci"
       name_store("rhod",L_rhod_) ;
#line 286 "transport_coef.loci"
       input("turbulent_transport,tmuu_f,rhod,rho_f,numSpecies") ;
#line 286 "transport_coef.loci"
       output("diff") ;
#line 286 "transport_coef.loci"
       constraint("area,const_diffusion") ;
#line 286 "transport_coef.loci"
    }
#line 286 "transport_coef.loci"
    void prelude(const Loci::sequence &seq) { 
    L_diff_.setVecSize(*L_numSpecies_) ;
	}    void calculate(Entity _e_) { 
#line 289 "transport_coef.loci"
     const int ns = L_numSpecies_[_e_];
      real rho = L_rho_f_[_e_];
      const real sct = L_turbulent_transport_[_e_].TurbulentSchmidtNumber ;
      const real mutsc = L_tmuu_f_[_e_]/sct ;                           
      for (int i =0;i <ns ;++i ) 
        L_diff_[_e_][i ]= (L_rhod_[_e_]+ mutsc ) / rho ;
      for (int i =0;i <ns ;++i )                                 
        L_diff_[_e_][i ]+= mutsc /rho ;                          
  }    void compute(const Loci::sequence &seq) { 
#line 297 "transport_coef.loci"
      prelude(seq) ;
#line 297 "transport_coef.loci"
      do_loop(seq,this) ;
#line 297 "transport_coef.loci"
    }
#line 297 "transport_coef.loci"
} ;
#line 297 "transport_coef.loci"
Loci::register_rule<file_transport_coef025_1280810697m375> register_file_transport_coef025_1280810697m375 ;
#line 297 "transport_coef.loci"
}
#line 297 "transport_coef.loci"
 

  namespace {class file_transport_coef026_1280810697m376 : public Loci::pointwise_rule {
#line 301 "transport_coef.loci"
    Loci::const_store<streamUns::real>  L_rho_f_ ; 
#line 301 "transport_coef.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 301 "transport_coef.loci"
    Loci::const_store<streamUns::real>  L_muu_f_ ; 
#line 301 "transport_coef.loci"
    Loci::const_store<streamUns::real>  L_tmuu_f_ ; 
#line 301 "transport_coef.loci"
    Loci::const_param<turbulent_transport>  L_turbulent_transport_ ; 
#line 301 "transport_coef.loci"
    Loci::const_param<real>  L_laminarSchmidtNumber_ ; 
#line 301 "transport_coef.loci"
    Loci::storeVec<streamUns::real>  L_diff_ ; 
#line 301 "transport_coef.loci"
public:
#line 301 "transport_coef.loci"
    file_transport_coef026_1280810697m376() {
#line 301 "transport_coef.loci"
       name_store("rho_f",L_rho_f_) ;
#line 301 "transport_coef.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 301 "transport_coef.loci"
       name_store("diff",L_diff_) ;
#line 301 "transport_coef.loci"
       name_store("muu_f",L_muu_f_) ;
#line 301 "transport_coef.loci"
       name_store("tmuu_f",L_tmuu_f_) ;
#line 301 "transport_coef.loci"
       name_store("turbulent_transport",L_turbulent_transport_) ;
#line 301 "transport_coef.loci"
       name_store("laminarSchmidtNumber",L_laminarSchmidtNumber_) ;
#line 301 "transport_coef.loci"
       input("numSpecies,turbulent_transport,tmuu_f,muu_f,rho_f,laminarSchmidtNumber") ;
#line 301 "transport_coef.loci"
       output("diff") ;
#line 301 "transport_coef.loci"
       constraint("schmidt_diffusion,muu_f,tmuu_f") ;
#line 301 "transport_coef.loci"
    }
#line 301 "transport_coef.loci"
    void prelude(const Loci::sequence &seq) { 
    L_diff_.setVecSize(*L_numSpecies_) ;
	}    void calculate(Entity _e_) { 
#line 304 "transport_coef.loci"
    const int ns = L_numSpecies_[_e_];
      const real sct = L_turbulent_transport_[_e_].TurbulentSchmidtNumber ;
      const real mutsc = L_tmuu_f_[_e_]/sct ;
      const real sc = L_laminarSchmidtNumber_[_e_];
      const real musc = L_muu_f_[_e_]/sc ;
      const real rrf = 1./L_rho_f_[_e_];
      for (int i =0;i <ns ;++i ) 
        L_diff_[_e_][i ]= (musc + mutsc ) * rrf ;
  }    void compute(const Loci::sequence &seq) { 
#line 312 "transport_coef.loci"
      prelude(seq) ;
#line 312 "transport_coef.loci"
      do_loop(seq,this) ;
#line 312 "transport_coef.loci"
    }
#line 312 "transport_coef.loci"
} ;
#line 312 "transport_coef.loci"
Loci::register_rule<file_transport_coef026_1280810697m376> register_file_transport_coef026_1280810697m376 ;
#line 312 "transport_coef.loci"
}
#line 312 "transport_coef.loci"
 
    
//-----------------------------------------------------------------------------
// Rule to set up the transport property database.
//-----------------------------------------------------------------------------

  class SetupTransportPropertyDatabase : public singleton_rule {
    private:
      const_param<EOS> eos ;
      param<TransportPropertyDatabase> transportPropertyDatabase ;
    public:

      // Define input and output.
      SetupTransportPropertyDatabase() {
        name_store("eos",eos) ;
        name_store("transportPropertyDatabase",transportPropertyDatabase) ;
        input("eos") ;
        output("transportPropertyDatabase") ;
        constraint("database") ;
      }

      // Read in the transport property data.
      void compute(const sequence &seq) {

        // Read the species from the species property files. Look in the
        // current directory if the database enviroment variable is not defined.
        string err("TRANSPORT PROPERTY DATABASE ERROR: ") ;
        const char *propertyBase=getenv("PROPERTY_DATABASE") ;
        char tmp_buf[512] ;
        if(propertyBase==NULL) {
          propertyBase = getenv("CHEMISTRY_DATABASE") ;
          if(propertyBase==NULL) {
            propertyBase="." ;
          } else {
            sprintf(tmp_buf,"%s/data_base/properties",propertyBase) ;
            propertyBase = tmp_buf ;
          }
        }
        const vector<string> &speciesName=eos->speciesNames() ;
        for(int i=0;i<eos->numSpecies();++i){
          cout << "Reading transport properties for species: "
            << speciesName[i] << endl ;
          string fileName=string(propertyBase)+"/"+speciesName[i]+".tran" ;
          ifstream in(fileName.c_str()) ;
          if(!in){
            cerr << err << "Trouble setting up transport property database!\n"
              << "Can't open file: " << fileName << endl ;
            Loci::Abort() ;
          }
          kill_white_space(in) ;
          if(get_token(in,speciesName[i].c_str())){
            kill_white_space(in) ;
            if(in.peek()!='='){
              cerr << err << "Expecting '=' in file " << fileName << endl ;
              Loci::Abort() ;
            }
            in.get() ;
          }else{
            cerr << err << "Expecting species name " << speciesName[i]
              << " in file " << fileName << endl ; Loci::Abort() ;
          }
          kill_white_space(in) ;
          transportPropertyDatabase->Input(in) ;
          transportPropertyDatabase->AddMolecularMass(eos->speciesMolecularMass
            (i)) ;
        }
        transportPropertyDatabase->Setup() ;
      }
  } ;

  register_rule<SetupTransportPropertyDatabase>
    registerSetupTransportPropertyDatabase ;


  // Get thermal conductivity from transport property database.
  class ConductivityFromDatabase : public pointwise_rule {
    private:
      const_param<int> numSpecies ;
      const_param<TransportPropertyDatabase> transportPropertyDatabase ;
      const_store<real> P ;
      const_store<real> T ;
      const_storeVec<real> mixture ;
      store<real> kconduct ;
    public:

      // Define input and output.
      ConductivityFromDatabase() {
        name_store("numSpecies",numSpecies) ;
        name_store("transportPropertyDatabase",transportPropertyDatabase) ;
        name_store("P",P) ;
        name_store("T",T) ;
        name_store("MIXTURE",mixture) ;
        name_store("kconduct(T,P,MIXTURE)",kconduct) ;
        input("numSpecies,transportPropertyDatabase") ;
        input("P,T,MIXTURE") ;
        output("kconduct(T,P,MIXTURE)") ;
        constraint("database,P,T") ;
      }

      // Get a single value from the database.
      void calculate(Entity e) {
        kconduct[e]=transportPropertyDatabase->Conductivity(P[e],T[e],
          mixture[e]) ;
      }

      // Loop over entities.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  // Register the rule.
  register_rule<ConductivityFromDatabase> registerConductivityFromDatabase ;

  // Get viscosity from transport property database.
  class ViscosityFromDatabase : public pointwise_rule {
    private:
      const_param<int> numSpecies ;
      const_param<TransportPropertyDatabase> transportPropertyDatabase ;
      const_param<std::vector<real> > speciesMolecularMass ;
      const_store<real> P ;
      const_store<real> T ;
      const_storeVec<real> mixture ;
      store<real> muu ;
    public:

      // Define input and output.
      ViscosityFromDatabase() {
        name_store("numSpecies",numSpecies) ;
        name_store("transportPropertyDatabase",transportPropertyDatabase) ;
        name_store("speciesMolecularMass",speciesMolecularMass) ;
        name_store("P",P) ;
        name_store("T",T) ;
        name_store("MIXTURE",mixture) ;
        name_store("muu(T,P,MIXTURE)",muu) ;
        input("numSpecies,transportPropertyDatabase,speciesMolecularMass") ;
        input("P,T,MIXTURE") ;
        output("muu(T,P,MIXTURE)") ;
        constraint("database,P,T") ;
      }

      // Get a single value from the database.
      void calculate(Entity e) {
        //int ns=*numSpecies ;
//cout << "before muu" << endl ;
        muu[e]=transportPropertyDatabase->Viscosity(P[e],T[e],mixture[e]) ; 
//cout << "after muu" << endl ;
//cout << "e,P,T,mix,muu: " << e << " " << P[e] << " " << T[e] ;
//for(int i=0;i<ns;++i) cout << " " << mixture[e][i] ;
//cout << " " << muu[e] << " " << endl ;
//Loci::Abort() ;
      }

      // Loop over entities.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  // Register the rule.
  register_rule<ViscosityFromDatabase> registerViscosityFromDatabase ;
}

