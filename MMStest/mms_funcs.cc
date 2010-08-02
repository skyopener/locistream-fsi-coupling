#line 1 "mms_funcs.loci"
#include <Loci.h>
#include <eos.h>
#include <qvi.h>
#include <readGrid.h>
#include <string>
#include <set>
#include <map>

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
#line 9 "mms_funcs.loci"

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
#line 10 "mms_funcs.loci"

using namespace streamUns ;
using std::string ;
using std::set ;
using std::map ;

// $type flowOperator param<Loci::exprP> 
// $type operatorDefs param<Loci::exprP> 

// $type SAOperator param<Loci::exprP> 
// $type tkeEqnOperator param<Loci::exprP> 
// $type omegaEqnOperator param<Loci::exprP> 

// $type wallDistanceFunction param<Loci::exprP> 

namespace {class file_mms_funcs000_1278521288m318 : public Loci::optional_rule {
#line 25 "mms_funcs.loci"
    Loci::param<Loci::exprP>  L_wallDistanceFunction_ ; 
#line 25 "mms_funcs.loci"
public:
#line 25 "mms_funcs.loci"
    file_mms_funcs000_1278521288m318() {
#line 25 "mms_funcs.loci"
       name_store("wallDistanceFunction",L_wallDistanceFunction_) ;
#line 25 "mms_funcs.loci"
       output("wallDistanceFunction") ;
#line 25 "mms_funcs.loci"
    }
#line 25 "mms_funcs.loci"
    void compute(const Loci::sequence &seq) { }} ;
#line 25 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs000_1278521288m318> register_file_mms_funcs000_1278521288m318 ;
#line 25 "mms_funcs.loci"
}
#line 25 "mms_funcs.loci"


namespace {class file_mms_funcs001_1278521288m319 : public Loci::optional_rule {
#line 27 "mms_funcs.loci"
    Loci::param<Loci::exprP>  L_flowOperator_ ; 
#line 27 "mms_funcs.loci"
public:
#line 27 "mms_funcs.loci"
    file_mms_funcs001_1278521288m319() {
#line 27 "mms_funcs.loci"
       name_store("flowOperator",L_flowOperator_) ;
#line 27 "mms_funcs.loci"
       output("flowOperator") ;
#line 27 "mms_funcs.loci"
    }
#line 27 "mms_funcs.loci"
    void compute(const Loci::sequence &seq) { }} ;
#line 27 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs001_1278521288m319> register_file_mms_funcs001_1278521288m319 ;
#line 27 "mms_funcs.loci"
}
#line 27 "mms_funcs.loci"


namespace {class file_mms_funcs002_1278521288m319 : public Loci::optional_rule {
#line 29 "mms_funcs.loci"
    Loci::param<Loci::exprP>  L_operatorDefs_ ; 
#line 29 "mms_funcs.loci"
public:
#line 29 "mms_funcs.loci"
    file_mms_funcs002_1278521288m319() {
#line 29 "mms_funcs.loci"
       name_store("operatorDefs",L_operatorDefs_) ;
#line 29 "mms_funcs.loci"
       output("operatorDefs") ;
#line 29 "mms_funcs.loci"
    }
#line 29 "mms_funcs.loci"
    void compute(const Loci::sequence &seq) { }} ;
#line 29 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs002_1278521288m319> register_file_mms_funcs002_1278521288m319 ;
#line 29 "mms_funcs.loci"
}
#line 29 "mms_funcs.loci"


namespace {class file_mms_funcs003_1278521288m319 : public Loci::optional_rule {
#line 31 "mms_funcs.loci"
    Loci::param<Loci::exprP>  L_SAOperator_ ; 
#line 31 "mms_funcs.loci"
public:
#line 31 "mms_funcs.loci"
    file_mms_funcs003_1278521288m319() {
#line 31 "mms_funcs.loci"
       name_store("SAOperator",L_SAOperator_) ;
#line 31 "mms_funcs.loci"
       output("SAOperator") ;
#line 31 "mms_funcs.loci"
    }
#line 31 "mms_funcs.loci"
    void compute(const Loci::sequence &seq) { }} ;
#line 31 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs003_1278521288m319> register_file_mms_funcs003_1278521288m319 ;
#line 31 "mms_funcs.loci"
}
#line 31 "mms_funcs.loci"


namespace {class file_mms_funcs004_1278521288m319 : public Loci::optional_rule {
#line 33 "mms_funcs.loci"
    Loci::param<Loci::exprP>  L_tkeEqnOperator_ ; 
#line 33 "mms_funcs.loci"
public:
#line 33 "mms_funcs.loci"
    file_mms_funcs004_1278521288m319() {
#line 33 "mms_funcs.loci"
       name_store("tkeEqnOperator",L_tkeEqnOperator_) ;
#line 33 "mms_funcs.loci"
       output("tkeEqnOperator") ;
#line 33 "mms_funcs.loci"
    }
#line 33 "mms_funcs.loci"
    void compute(const Loci::sequence &seq) { }} ;
#line 33 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs004_1278521288m319> register_file_mms_funcs004_1278521288m319 ;
#line 33 "mms_funcs.loci"
}
#line 33 "mms_funcs.loci"


namespace {class file_mms_funcs005_1278521288m320 : public Loci::optional_rule {
#line 35 "mms_funcs.loci"
    Loci::param<Loci::exprP>  L_omegaEqnOperator_ ; 
#line 35 "mms_funcs.loci"
public:
#line 35 "mms_funcs.loci"
    file_mms_funcs005_1278521288m320() {
#line 35 "mms_funcs.loci"
       name_store("omegaEqnOperator",L_omegaEqnOperator_) ;
#line 35 "mms_funcs.loci"
       output("omegaEqnOperator") ;
#line 35 "mms_funcs.loci"
    }
#line 35 "mms_funcs.loci"
    void compute(const Loci::sequence &seq) { }} ;
#line 35 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs005_1278521288m320> register_file_mms_funcs005_1278521288m320 ;
#line 35 "mms_funcs.loci"
}
#line 35 "mms_funcs.loci"



namespace Loci {
  void getVarNames(exprP e, set<string> &namelist) ;
}

// Determine if the prescribed function is time dependent.
// $type MMSTimeDependent Constraint
// $type MMSSteadyState Constraint
namespace {class file_mms_funcs006_1278521288m320 : public Loci::constraint_rule {
#line 45 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_operatorDefs_ ; 
#line 45 "mms_funcs.loci"
    Loci::Constraint L_MMSTimeDependent_ ; 
#line 45 "mms_funcs.loci"
    Loci::Constraint L_MMSSteadyState_ ; 
#line 45 "mms_funcs.loci"
public:
#line 45 "mms_funcs.loci"
    file_mms_funcs006_1278521288m320() {
#line 45 "mms_funcs.loci"
       name_store("operatorDefs",L_operatorDefs_) ;
#line 45 "mms_funcs.loci"
       name_store("MMSTimeDependent",L_MMSTimeDependent_) ;
#line 45 "mms_funcs.loci"
       name_store("MMSSteadyState",L_MMSSteadyState_) ;
#line 45 "mms_funcs.loci"
       input("operatorDefs") ;
#line 45 "mms_funcs.loci"
       output("MMSTimeDependent") ;
#line 45 "mms_funcs.loci"
       output("MMSSteadyState") ;
#line 45 "mms_funcs.loci"
    }
#line 45 "mms_funcs.loci"
    void compute(const Loci::sequence &seq) { 
  Loci::exprP vars = Loci::expression::create(string("(rho,u,v,w,p,T)")) ;
  Loci::exprP p = Loci::substitutionEngine(vars,(*L_operatorDefs_)) ;
  std::set<string> namelist ;
  Loci::getVarNames(p,namelist) ;
  L_MMSTimeDependent_= EMPTY ;
  L_MMSSteadyState_= EMPTY ;
  if(namelist.find(string("t")) != namelist.end())
    L_MMSTimeDependent_= L_operatorDefs_.domain() ;
  else
    L_MMSSteadyState_= L_operatorDefs_.domain() ;
}} ;
#line 56 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs006_1278521288m320> register_file_mms_funcs006_1278521288m320 ;
#line 56 "mms_funcs.loci"
}
#line 56 "mms_funcs.loci"


// $type location store<vect3d> 

// $type MMS_rhoFunc param<Loci::exprP> 
// $type MMS_velFunc param<Loci::exprP> 
// $type MMS_pFunc param<Loci::exprP> 

namespace {class file_mms_funcs007_1278521288m321 : public Loci::singleton_rule {
#line 64 "mms_funcs.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 64 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_operatorDefs_ ; 
#line 64 "mms_funcs.loci"
    Loci::param<Loci::exprP>  L_MMS_rhoFunc_ ; 
#line 64 "mms_funcs.loci"
public:
#line 64 "mms_funcs.loci"
    file_mms_funcs007_1278521288m321() {
#line 64 "mms_funcs.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 64 "mms_funcs.loci"
       name_store("operatorDefs",L_operatorDefs_) ;
#line 64 "mms_funcs.loci"
       name_store("MMS_rhoFunc",L_MMS_rhoFunc_) ;
#line 64 "mms_funcs.loci"
       input("operatorDefs,numSpecies") ;
#line 64 "mms_funcs.loci"
       output("MMS_rhoFunc") ;
#line 64 "mms_funcs.loci"
    }
#line 64 "mms_funcs.loci"
    void compute(const Loci::sequence &seq) { 
  if((*L_numSpecies_)== 1) {
    Loci::exprP vars = Loci::expression::create("rho") ;
    (*L_MMS_rhoFunc_)= Loci::substitutionEngine(vars,(*L_operatorDefs_)) ;
  } else {
    string s ;
    s = "(rho*Y0" ;
    for(int i=1;i<(*L_numSpecies_);++i) {
      s+= ",rho*Y" ;
      s+= char('0'+i) ;
    }
    s+= ")" ;
    Loci::exprP vars = Loci::expression::create(s) ;
    (*L_MMS_rhoFunc_)= Loci::substitutionEngine(vars,(*L_operatorDefs_)) ;
  }
}} ;
#line 79 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs007_1278521288m321> register_file_mms_funcs007_1278521288m321 ;
#line 79 "mms_funcs.loci"
}
#line 79 "mms_funcs.loci"

  
namespace {class file_mms_funcs008_1278521288m321 : public Loci::pointwise_rule {
#line 82 "mms_funcs.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 82 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 82 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_MMS_rhoFunc_ ; 
#line 82 "mms_funcs.loci"
    Loci::storeVec<real>  L_mmsRholocation_ ; 
#line 82 "mms_funcs.loci"
public:
#line 82 "mms_funcs.loci"
    file_mms_funcs008_1278521288m321() {
#line 82 "mms_funcs.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 82 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 82 "mms_funcs.loci"
       name_store("MMS_rhoFunc",L_MMS_rhoFunc_) ;
#line 82 "mms_funcs.loci"
       name_store("mmsRho(location)",L_mmsRholocation_) ;
#line 82 "mms_funcs.loci"
       input("location,numSpecies,MMS_rhoFunc") ;
#line 82 "mms_funcs.loci"
       output("mmsRho(location)") ;
#line 82 "mms_funcs.loci"
       constraint("MMSSteadyState,location") ;
#line 82 "mms_funcs.loci"
    }
#line 82 "mms_funcs.loci"
    void prelude(const Loci::sequence &seq) { 
  L_mmsRholocation_.setVecSize(*L_numSpecies_) ;
}    void calculate(Entity _e_) { 
#line 85 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  if (L_numSpecies_[_e_]== 1) {
    L_mmsRholocation_[_e_][0]= (L_MMS_rhoFunc_[_e_])->evaluate (varmap ) ;
  } else {
    int i =0 ;
    Loci ::expression ::exprList ::const_iterator li ;

    for (li =(L_MMS_rhoFunc_[_e_])->expr_list .begin ();
        li !=(L_MMS_rhoFunc_[_e_])->expr_list .end ();++li ) {
      L_mmsRholocation_[_e_][i ]= (*li )->evaluate (varmap ) ;
    }
  }
}    void compute(const Loci::sequence &seq) { 
#line 100 "mms_funcs.loci"
      prelude(seq) ;
#line 100 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 100 "mms_funcs.loci"
    }
#line 100 "mms_funcs.loci"
} ;
#line 100 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs008_1278521288m321> register_file_mms_funcs008_1278521288m321 ;
#line 100 "mms_funcs.loci"
}
#line 100 "mms_funcs.loci"


namespace {class file_mms_funcs009_1278521288m322 : public Loci::pointwise_rule {
#line 103 "mms_funcs.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 103 "mms_funcs.loci"
    Loci::const_param<streamUns::real>  L_stime_ ; 
#line 103 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 103 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_MMS_rhoFunc_ ; 
#line 103 "mms_funcs.loci"
    Loci::storeVec<real>  L_mmsRholocation_ ; 
#line 103 "mms_funcs.loci"
public:
#line 103 "mms_funcs.loci"
    file_mms_funcs009_1278521288m322() {
#line 103 "mms_funcs.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 103 "mms_funcs.loci"
       name_store("stime",L_stime_) ;
#line 103 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 103 "mms_funcs.loci"
       name_store("MMS_rhoFunc",L_MMS_rhoFunc_) ;
#line 103 "mms_funcs.loci"
       name_store("mmsRho(location)",L_mmsRholocation_) ;
#line 103 "mms_funcs.loci"
       input("location,stime,numSpecies,MMS_rhoFunc") ;
#line 103 "mms_funcs.loci"
       output("mmsRho(location)") ;
#line 103 "mms_funcs.loci"
       constraint("MMSTimeDependent,location") ;
#line 103 "mms_funcs.loci"
    }
#line 103 "mms_funcs.loci"
    void prelude(const Loci::sequence &seq) { 
  L_mmsRholocation_.setVecSize(*L_numSpecies_) ;
}    void calculate(Entity _e_) { 
#line 106 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  varmap ["t"] = L_stime_[_e_];
  if (L_numSpecies_[_e_]== 1) {
    L_mmsRholocation_[_e_][0]= (L_MMS_rhoFunc_[_e_])->evaluate (varmap ) ;
  } else {
    int i =0 ;
    Loci ::expression ::exprList ::const_iterator li ;

    for (li =(L_MMS_rhoFunc_[_e_])->expr_list .begin ();
        li !=(L_MMS_rhoFunc_[_e_])->expr_list .end ();++li ) {
      L_mmsRholocation_[_e_][i ]= (*li )->evaluate (varmap ) ;
    }
  }
}    void compute(const Loci::sequence &seq) { 
#line 122 "mms_funcs.loci"
      prelude(seq) ;
#line 122 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 122 "mms_funcs.loci"
    }
#line 122 "mms_funcs.loci"
} ;
#line 122 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs009_1278521288m322> register_file_mms_funcs009_1278521288m322 ;
#line 122 "mms_funcs.loci"
}
#line 122 "mms_funcs.loci"


namespace {class file_mms_funcs010_1278521288m323 : public Loci::singleton_rule {
#line 124 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_operatorDefs_ ; 
#line 124 "mms_funcs.loci"
    Loci::param<Loci::exprP>  L_MMS_velFunc_ ; 
#line 124 "mms_funcs.loci"
public:
#line 124 "mms_funcs.loci"
    file_mms_funcs010_1278521288m323() {
#line 124 "mms_funcs.loci"
       name_store("operatorDefs",L_operatorDefs_) ;
#line 124 "mms_funcs.loci"
       name_store("MMS_velFunc",L_MMS_velFunc_) ;
#line 124 "mms_funcs.loci"
       input("operatorDefs") ;
#line 124 "mms_funcs.loci"
       output("MMS_velFunc") ;
#line 124 "mms_funcs.loci"
    }
#line 124 "mms_funcs.loci"
    void compute(const Loci::sequence &seq) { 
  Loci::exprP vars = Loci::expression::create("(u,v,w)") ;
  (*L_MMS_velFunc_)= Loci::substitutionEngine(vars,(*L_operatorDefs_)) ;
}} ;
#line 127 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs010_1278521288m323> register_file_mms_funcs010_1278521288m323 ;
#line 127 "mms_funcs.loci"
}
#line 127 "mms_funcs.loci"


namespace {class file_mms_funcs011_1278521288m324 : public Loci::pointwise_rule {
#line 131 "mms_funcs.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 131 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 131 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_MMS_velFunc_ ; 
#line 131 "mms_funcs.loci"
    Loci::store<vect3d>  L_mmsVelocitylocation_ ; 
#line 131 "mms_funcs.loci"
public:
#line 131 "mms_funcs.loci"
    file_mms_funcs011_1278521288m324() {
#line 131 "mms_funcs.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 131 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 131 "mms_funcs.loci"
       name_store("MMS_velFunc",L_MMS_velFunc_) ;
#line 131 "mms_funcs.loci"
       name_store("mmsVelocity(location)",L_mmsVelocitylocation_) ;
#line 131 "mms_funcs.loci"
       input("location,numSpecies,MMS_velFunc") ;
#line 131 "mms_funcs.loci"
       output("mmsVelocity(location)") ;
#line 131 "mms_funcs.loci"
       constraint("MMSSteadyState,location") ;
#line 131 "mms_funcs.loci"
    }
#line 131 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 132 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  Loci ::expression ::exprList ::const_iterator li ;
  li = (L_MMS_velFunc_[_e_])->expr_list .begin () ;
  L_mmsVelocitylocation_[_e_].x = (*li )->evaluate (varmap ) ;
  li ++ ;
  L_mmsVelocitylocation_[_e_].y = (*li )->evaluate (varmap ) ;
  li ++ ;
  L_mmsVelocitylocation_[_e_].z = (*li )->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 143 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 143 "mms_funcs.loci"
    }
#line 143 "mms_funcs.loci"
} ;
#line 143 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs011_1278521288m324> register_file_mms_funcs011_1278521288m324 ;
#line 143 "mms_funcs.loci"
}
#line 143 "mms_funcs.loci"


namespace {class file_mms_funcs012_1278521288m325 : public Loci::pointwise_rule {
#line 147 "mms_funcs.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 147 "mms_funcs.loci"
    Loci::const_param<streamUns::real>  L_stime_ ; 
#line 147 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 147 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_MMS_velFunc_ ; 
#line 147 "mms_funcs.loci"
    Loci::store<vect3d>  L_mmsVelocitylocation_ ; 
#line 147 "mms_funcs.loci"
public:
#line 147 "mms_funcs.loci"
    file_mms_funcs012_1278521288m325() {
#line 147 "mms_funcs.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 147 "mms_funcs.loci"
       name_store("stime",L_stime_) ;
#line 147 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 147 "mms_funcs.loci"
       name_store("MMS_velFunc",L_MMS_velFunc_) ;
#line 147 "mms_funcs.loci"
       name_store("mmsVelocity(location)",L_mmsVelocitylocation_) ;
#line 147 "mms_funcs.loci"
       input("location,stime,numSpecies,MMS_velFunc") ;
#line 147 "mms_funcs.loci"
       output("mmsVelocity(location)") ;
#line 147 "mms_funcs.loci"
       constraint("MMSTimeDependent,location") ;
#line 147 "mms_funcs.loci"
    }
#line 147 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 148 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  varmap ["t"] = L_stime_[_e_];
  Loci ::expression ::exprList ::const_iterator li ;
  li = (L_MMS_velFunc_[_e_])->expr_list .begin () ;
  L_mmsVelocitylocation_[_e_].x = (*li )->evaluate (varmap ) ;
  li ++ ;
  L_mmsVelocitylocation_[_e_].y = (*li )->evaluate (varmap ) ;
  li ++ ;
  L_mmsVelocitylocation_[_e_].z = (*li )->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 160 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 160 "mms_funcs.loci"
    }
#line 160 "mms_funcs.loci"
} ;
#line 160 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs012_1278521288m325> register_file_mms_funcs012_1278521288m325 ;
#line 160 "mms_funcs.loci"
}
#line 160 "mms_funcs.loci"


namespace {class file_mms_funcs013_1278521288m325 : public Loci::singleton_rule {
#line 162 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_operatorDefs_ ; 
#line 162 "mms_funcs.loci"
    Loci::param<Loci::exprP>  L_MMS_pFunc_ ; 
#line 162 "mms_funcs.loci"
public:
#line 162 "mms_funcs.loci"
    file_mms_funcs013_1278521288m325() {
#line 162 "mms_funcs.loci"
       name_store("operatorDefs",L_operatorDefs_) ;
#line 162 "mms_funcs.loci"
       name_store("MMS_pFunc",L_MMS_pFunc_) ;
#line 162 "mms_funcs.loci"
       input("operatorDefs") ;
#line 162 "mms_funcs.loci"
       output("MMS_pFunc") ;
#line 162 "mms_funcs.loci"
    }
#line 162 "mms_funcs.loci"
    void compute(const Loci::sequence &seq) { 
  Loci::exprP vars = Loci::expression::create("p") ;
  (*L_MMS_pFunc_)= Loci::substitutionEngine(vars,(*L_operatorDefs_)) ;
}} ;
#line 165 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs013_1278521288m325> register_file_mms_funcs013_1278521288m325 ;
#line 165 "mms_funcs.loci"
}
#line 165 "mms_funcs.loci"


namespace {class file_mms_funcs014_1278521288m326 : public Loci::pointwise_rule {
#line 169 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 169 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_MMS_pFunc_ ; 
#line 169 "mms_funcs.loci"
    Loci::store<real>  L_mmsPressurelocation_ ; 
#line 169 "mms_funcs.loci"
public:
#line 169 "mms_funcs.loci"
    file_mms_funcs014_1278521288m326() {
#line 169 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 169 "mms_funcs.loci"
       name_store("MMS_pFunc",L_MMS_pFunc_) ;
#line 169 "mms_funcs.loci"
       name_store("mmsPressure(location)",L_mmsPressurelocation_) ;
#line 169 "mms_funcs.loci"
       input("location,MMS_pFunc") ;
#line 169 "mms_funcs.loci"
       output("mmsPressure(location)") ;
#line 169 "mms_funcs.loci"
       constraint("MMSSteadyState,location") ;
#line 169 "mms_funcs.loci"
    }
#line 169 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 170 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  L_mmsPressurelocation_[_e_]= (L_MMS_pFunc_[_e_])->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 175 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 175 "mms_funcs.loci"
    }
#line 175 "mms_funcs.loci"
} ;
#line 175 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs014_1278521288m326> register_file_mms_funcs014_1278521288m326 ;
#line 175 "mms_funcs.loci"
}
#line 175 "mms_funcs.loci"



namespace {class file_mms_funcs015_1278521288m326 : public Loci::pointwise_rule {
#line 180 "mms_funcs.loci"
    Loci::const_param<streamUns::real>  L_stime_ ; 
#line 180 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 180 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_MMS_pFunc_ ; 
#line 180 "mms_funcs.loci"
    Loci::store<real>  L_mmsPressurelocation_ ; 
#line 180 "mms_funcs.loci"
public:
#line 180 "mms_funcs.loci"
    file_mms_funcs015_1278521288m326() {
#line 180 "mms_funcs.loci"
       name_store("stime",L_stime_) ;
#line 180 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 180 "mms_funcs.loci"
       name_store("MMS_pFunc",L_MMS_pFunc_) ;
#line 180 "mms_funcs.loci"
       name_store("mmsPressure(location)",L_mmsPressurelocation_) ;
#line 180 "mms_funcs.loci"
       input("location,stime,MMS_pFunc") ;
#line 180 "mms_funcs.loci"
       output("mmsPressure(location)") ;
#line 180 "mms_funcs.loci"
       constraint("MMSTimeDependent,location") ;
#line 180 "mms_funcs.loci"
    }
#line 180 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 181 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  varmap ["t"] = L_stime_[_e_];
  L_mmsPressurelocation_[_e_]= (L_MMS_pFunc_[_e_])->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 187 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 187 "mms_funcs.loci"
    }
#line 187 "mms_funcs.loci"
} ;
#line 187 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs015_1278521288m326> register_file_mms_funcs015_1278521288m326 ;
#line 187 "mms_funcs.loci"
}
#line 187 "mms_funcs.loci"



// $type meanFlowSrcFunc param<Loci::exprP> 

namespace {class file_mms_funcs016_1278521288m327 : public Loci::singleton_rule {
#line 192 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_flowOperator_ ; 
#line 192 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_operatorDefs_ ; 
#line 192 "mms_funcs.loci"
    Loci::param<Loci::exprP>  L_meanFlowSrcFunc_ ; 
#line 192 "mms_funcs.loci"
public:
#line 192 "mms_funcs.loci"
    file_mms_funcs016_1278521288m327() {
#line 192 "mms_funcs.loci"
       name_store("flowOperator",L_flowOperator_) ;
#line 192 "mms_funcs.loci"
       name_store("operatorDefs",L_operatorDefs_) ;
#line 192 "mms_funcs.loci"
       name_store("meanFlowSrcFunc",L_meanFlowSrcFunc_) ;
#line 192 "mms_funcs.loci"
       input("flowOperator,operatorDefs") ;
#line 192 "mms_funcs.loci"
       output("meanFlowSrcFunc") ;
#line 192 "mms_funcs.loci"
    }
#line 192 "mms_funcs.loci"
    void compute(const Loci::sequence &seq) { 
    
  Loci::exprP sub = Loci::substitutionEngine((*L_flowOperator_),(*L_operatorDefs_)) ;
  sub = sub->simplify() ; // Simplify expression before evaluating 
  // derivative expressions
  sub = sub->symbolic_eval() ; // evaluate symbolic operators (derivatives)
  sub = sub->simplify() ; // simplify result
  (*L_meanFlowSrcFunc_)= sub ;
  if(Loci::MPI_rank == 0) {
    ofstream ofile("flowOperator.dat",ios::out) ;
    
    ofstream defsfile("funcDefs.dat",ios::out) ;
    ofstream flowfile("flowSourceDefs.dat",ios::out) ;
    ofile << "Flow Operator:" << endl ;

    Loci::expression::exprList::const_iterator li ;
    if((*L_flowOperator_)->op != Loci::OP_COMMA) {
      cerr << "flowOperator should be a comma separated list!"<< endl ;
      Loci::Abort() ;
    }
    li = ((*L_flowOperator_))->expr_list.begin() ;
    ofile << "(" << endl ;
    if(li != ((*L_flowOperator_))->expr_list.end()) {// print first item
      ofile << *li ;
    }
    
    for(++li;li != ((*L_flowOperator_))->expr_list.end();++li) {//print rest of list
      ofile << "," << endl << *li ;
    }
    ofile << endl << ")" << endl ;
      
    
    defsfile << "Operator definitions:" << endl ;

    if((*L_operatorDefs_)->op != Loci::OP_COMMA) {
      cerr << "operatorDefs should be a comma separated list!"<< endl ;
      Loci::Abort() ;
    }
    
    li = ((*L_operatorDefs_))->expr_list.begin() ;
    defsfile << "(" << endl ;
    if(li != ((*L_operatorDefs_))->expr_list.end()) { // print first item
      defsfile << *li ;
    }
    for(++li;li != ((*L_operatorDefs_))->expr_list.end();++li) {// print rest of list
      defsfile << "," << endl << *li ;
    }
    defsfile << endl << ")" << endl ;

    flowfile << "mean flow source terms (derived and simplified)" << endl ;

    if(sub->op != Loci::OP_COMMA) {
      cerr << "flow source terms should be a comma separated list!, this shouldn't happen!"<< endl ;
      Loci::Abort() ;
    }
    li = (sub)->expr_list.begin() ;
    flowfile << "(" << endl ;
    if(li != (sub)->expr_list.end()) { // print first item
      flowfile << *li ;
    }
    for(++li;li != (sub)->expr_list.end();++li) { // print rest of list
      flowfile << ","<< endl << *li ;
    }
    flowfile << endl << ")" << endl ;
  }    
}} ;
#line 257 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs016_1278521288m327> register_file_mms_funcs016_1278521288m327 ;
#line 257 "mms_funcs.loci"
}
#line 257 "mms_funcs.loci"


namespace {class file_mms_funcs017_1278521288m327 : public Loci::pointwise_rule {
#line 259 "mms_funcs.loci"
    Loci::const_store<Loci::vector3d<Loci::real_t> >  L_cellcenter_ ; 
#line 259 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_wallDistanceFunction_ ; 
#line 259 "mms_funcs.loci"
    Loci::store<streamUns::real>  L_mmsdist_noslip_ ; 
#line 259 "mms_funcs.loci"
public:
#line 259 "mms_funcs.loci"
    file_mms_funcs017_1278521288m327() {
#line 259 "mms_funcs.loci"
       name_store("cellcenter",L_cellcenter_) ;
#line 259 "mms_funcs.loci"
       name_store("wallDistanceFunction",L_wallDistanceFunction_) ;
#line 259 "mms_funcs.loci"
       name_store("mms::dist_noslip",L_mmsdist_noslip_) ;
#line 259 "mms_funcs.loci"
       input("cellcenter,wallDistanceFunction") ;
#line 259 "mms_funcs.loci"
       output("mms::dist_noslip") ;
#line 259 "mms_funcs.loci"
    }
#line 259 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 260 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_cellcenter_[_e_].x ;
  varmap ["y"] = L_cellcenter_[_e_].y ;
  varmap ["z"] = L_cellcenter_[_e_].z ;
  L_mmsdist_noslip_[_e_]= (L_wallDistanceFunction_[_e_])->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 265 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 265 "mms_funcs.loci"
    }
#line 265 "mms_funcs.loci"
} ;
#line 265 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs017_1278521288m327> register_file_mms_funcs017_1278521288m327 ;
#line 265 "mms_funcs.loci"
}
#line 265 "mms_funcs.loci"

  
namespace {class file_mms_funcs018_1278521288m328 : public Loci::pointwise_rule {
#line 268 "mms_funcs.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 268 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 268 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_meanFlowSrcFunc_ ; 
#line 268 "mms_funcs.loci"
    Loci::storeVec<real>  L_massGeneralSrclocation_ ; 
#line 268 "mms_funcs.loci"
public:
#line 268 "mms_funcs.loci"
    file_mms_funcs018_1278521288m328() {
#line 268 "mms_funcs.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 268 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 268 "mms_funcs.loci"
       name_store("meanFlowSrcFunc",L_meanFlowSrcFunc_) ;
#line 268 "mms_funcs.loci"
       name_store("massGeneralSrc(location)",L_massGeneralSrclocation_) ;
#line 268 "mms_funcs.loci"
       input("location,numSpecies,meanFlowSrcFunc") ;
#line 268 "mms_funcs.loci"
       output("massGeneralSrc(location)") ;
#line 268 "mms_funcs.loci"
       constraint("MMSSteadyState,location") ;
#line 268 "mms_funcs.loci"
    }
#line 268 "mms_funcs.loci"
    void prelude(const Loci::sequence &seq) { 
  L_massGeneralSrclocation_.setVecSize(*L_numSpecies_) ;
}    void calculate(Entity _e_) { 
#line 271 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  Loci ::expression ::exprList ::const_iterator li ;
  li = (L_meanFlowSrcFunc_[_e_])->expr_list .begin () ;
  
  for (int i =0;i <L_numSpecies_[_e_];++i ) {
    L_massGeneralSrclocation_[_e_][i ]= (*li )->evaluate (varmap ) ;
    li ++ ;
  }
}    void compute(const Loci::sequence &seq) { 
#line 282 "mms_funcs.loci"
      prelude(seq) ;
#line 282 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 282 "mms_funcs.loci"
    }
#line 282 "mms_funcs.loci"
} ;
#line 282 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs018_1278521288m328> register_file_mms_funcs018_1278521288m328 ;
#line 282 "mms_funcs.loci"
}
#line 282 "mms_funcs.loci"


namespace {class file_mms_funcs019_1278521288m328 : public Loci::pointwise_rule {
#line 285 "mms_funcs.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 285 "mms_funcs.loci"
    Loci::const_param<streamUns::real>  L_stime_ ; 
#line 285 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 285 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_meanFlowSrcFunc_ ; 
#line 285 "mms_funcs.loci"
    Loci::storeVec<real>  L_massGeneralSrclocation_ ; 
#line 285 "mms_funcs.loci"
public:
#line 285 "mms_funcs.loci"
    file_mms_funcs019_1278521288m328() {
#line 285 "mms_funcs.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 285 "mms_funcs.loci"
       name_store("stime",L_stime_) ;
#line 285 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 285 "mms_funcs.loci"
       name_store("meanFlowSrcFunc",L_meanFlowSrcFunc_) ;
#line 285 "mms_funcs.loci"
       name_store("massGeneralSrc(location)",L_massGeneralSrclocation_) ;
#line 285 "mms_funcs.loci"
       input("location,stime,numSpecies,meanFlowSrcFunc") ;
#line 285 "mms_funcs.loci"
       output("massGeneralSrc(location)") ;
#line 285 "mms_funcs.loci"
       constraint("MMSTimeDependent,location") ;
#line 285 "mms_funcs.loci"
    }
#line 285 "mms_funcs.loci"
    void prelude(const Loci::sequence &seq) { 
  L_massGeneralSrclocation_.setVecSize(*L_numSpecies_) ;
}    void calculate(Entity _e_) { 
#line 288 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  varmap ["t"] = L_stime_[_e_];
  Loci ::expression ::exprList ::const_iterator li ;
  li = (L_meanFlowSrcFunc_[_e_])->expr_list .begin () ;
  
  for (int i =0;i <L_numSpecies_[_e_];++i ) {
    L_massGeneralSrclocation_[_e_][i ]= (*li )->evaluate (varmap ) ;
    li ++ ;
  }
}    void compute(const Loci::sequence &seq) { 
#line 300 "mms_funcs.loci"
      prelude(seq) ;
#line 300 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 300 "mms_funcs.loci"
    }
#line 300 "mms_funcs.loci"
} ;
#line 300 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs019_1278521288m328> register_file_mms_funcs019_1278521288m328 ;
#line 300 "mms_funcs.loci"
}
#line 300 "mms_funcs.loci"


namespace {class file_mms_funcs020_1278521288m329 : public Loci::pointwise_rule {
#line 303 "mms_funcs.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 303 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 303 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_meanFlowSrcFunc_ ; 
#line 303 "mms_funcs.loci"
    Loci::store<vect3d>  L_momentumGeneralSrclocation_ ; 
#line 303 "mms_funcs.loci"
public:
#line 303 "mms_funcs.loci"
    file_mms_funcs020_1278521288m329() {
#line 303 "mms_funcs.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 303 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 303 "mms_funcs.loci"
       name_store("meanFlowSrcFunc",L_meanFlowSrcFunc_) ;
#line 303 "mms_funcs.loci"
       name_store("momentumGeneralSrc(location)",L_momentumGeneralSrclocation_) ;
#line 303 "mms_funcs.loci"
       input("location,numSpecies,meanFlowSrcFunc") ;
#line 303 "mms_funcs.loci"
       output("momentumGeneralSrc(location)") ;
#line 303 "mms_funcs.loci"
       constraint("MMSSteadyState,location") ;
#line 303 "mms_funcs.loci"
    }
#line 303 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 304 "mms_funcs.loci"
  // Momemtum source term  in  Kg/(m^2 sec^2)
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  Loci ::expression ::exprList ::const_iterator li ;
  li = (L_meanFlowSrcFunc_[_e_])->expr_list .begin () ;
  
  // Skipping the species functions to get to momentum functions. JW
  for (int i =0;i <L_numSpecies_[_e_];++i ) {
    li ++ ;
  }

  L_momentumGeneralSrclocation_[_e_].x = (*li )->evaluate (varmap ) ;
  li ++ ;
  L_momentumGeneralSrclocation_[_e_].y = (*li )->evaluate (varmap ) ;
  li ++ ;
  L_momentumGeneralSrclocation_[_e_].z = (*li )->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 322 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 322 "mms_funcs.loci"
    }
#line 322 "mms_funcs.loci"
} ;
#line 322 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs020_1278521288m329> register_file_mms_funcs020_1278521288m329 ;
#line 322 "mms_funcs.loci"
}
#line 322 "mms_funcs.loci"


namespace {class file_mms_funcs021_1278521288m329 : public Loci::pointwise_rule {
#line 325 "mms_funcs.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 325 "mms_funcs.loci"
    Loci::const_param<streamUns::real>  L_stime_ ; 
#line 325 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 325 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_meanFlowSrcFunc_ ; 
#line 325 "mms_funcs.loci"
    Loci::store<vect3d>  L_momentumGeneralSrclocation_ ; 
#line 325 "mms_funcs.loci"
public:
#line 325 "mms_funcs.loci"
    file_mms_funcs021_1278521288m329() {
#line 325 "mms_funcs.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 325 "mms_funcs.loci"
       name_store("stime",L_stime_) ;
#line 325 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 325 "mms_funcs.loci"
       name_store("meanFlowSrcFunc",L_meanFlowSrcFunc_) ;
#line 325 "mms_funcs.loci"
       name_store("momentumGeneralSrc(location)",L_momentumGeneralSrclocation_) ;
#line 325 "mms_funcs.loci"
       input("location,stime,numSpecies,meanFlowSrcFunc") ;
#line 325 "mms_funcs.loci"
       output("momentumGeneralSrc(location)") ;
#line 325 "mms_funcs.loci"
       constraint("MMSTimeDependent,location") ;
#line 325 "mms_funcs.loci"
    }
#line 325 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 326 "mms_funcs.loci"
  // Momemtum source term  in  Kg/(m^2 sec^2)
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  varmap ["t"] = L_stime_[_e_];
  Loci ::expression ::exprList ::const_iterator li ;
  li = (L_meanFlowSrcFunc_[_e_])->expr_list .begin () ;
  
  for (int i =0;i <L_numSpecies_[_e_];++i ) {
    li ++ ;
  }

  L_momentumGeneralSrclocation_[_e_].x = (*li )->evaluate (varmap ) ;
  li ++ ;
  L_momentumGeneralSrclocation_[_e_].y = (*li )->evaluate (varmap ) ;
  li ++ ;
  L_momentumGeneralSrclocation_[_e_].z = (*li )->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 344 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 344 "mms_funcs.loci"
    }
#line 344 "mms_funcs.loci"
} ;
#line 344 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs021_1278521288m329> register_file_mms_funcs021_1278521288m329 ;
#line 344 "mms_funcs.loci"
}
#line 344 "mms_funcs.loci"


namespace {class file_mms_funcs022_1278521288m330 : public Loci::pointwise_rule {
#line 347 "mms_funcs.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 347 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 347 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_meanFlowSrcFunc_ ; 
#line 347 "mms_funcs.loci"
    Loci::store<real>  L_energyGeneralSrclocation_ ; 
#line 347 "mms_funcs.loci"
public:
#line 347 "mms_funcs.loci"
    file_mms_funcs022_1278521288m330() {
#line 347 "mms_funcs.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 347 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 347 "mms_funcs.loci"
       name_store("meanFlowSrcFunc",L_meanFlowSrcFunc_) ;
#line 347 "mms_funcs.loci"
       name_store("energyGeneralSrc(location)",L_energyGeneralSrclocation_) ;
#line 347 "mms_funcs.loci"
       input("location,numSpecies,meanFlowSrcFunc") ;
#line 347 "mms_funcs.loci"
       output("energyGeneralSrc(location)") ;
#line 347 "mms_funcs.loci"
       constraint("MMSSteadyState,location") ;
#line 347 "mms_funcs.loci"
    }
#line 347 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 348 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  Loci ::expression ::exprList ::const_iterator li ;
  li = (L_meanFlowSrcFunc_[_e_])->expr_list .begin () ;
  
  for (int i =0;i <L_numSpecies_[_e_]+3;++i ) {
    li ++ ;
  }

  // energy source term in Kg J/(m^3 sec)
  L_energyGeneralSrclocation_[_e_]= (*li )->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 361 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 361 "mms_funcs.loci"
    }
#line 361 "mms_funcs.loci"
} ;
#line 361 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs022_1278521288m330> register_file_mms_funcs022_1278521288m330 ;
#line 361 "mms_funcs.loci"
}
#line 361 "mms_funcs.loci"


namespace {class file_mms_funcs023_1278521288m330 : public Loci::pointwise_rule {
#line 364 "mms_funcs.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 364 "mms_funcs.loci"
    Loci::const_param<streamUns::real>  L_stime_ ; 
#line 364 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 364 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_meanFlowSrcFunc_ ; 
#line 364 "mms_funcs.loci"
    Loci::store<real>  L_energyGeneralSrclocation_ ; 
#line 364 "mms_funcs.loci"
public:
#line 364 "mms_funcs.loci"
    file_mms_funcs023_1278521288m330() {
#line 364 "mms_funcs.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 364 "mms_funcs.loci"
       name_store("stime",L_stime_) ;
#line 364 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 364 "mms_funcs.loci"
       name_store("meanFlowSrcFunc",L_meanFlowSrcFunc_) ;
#line 364 "mms_funcs.loci"
       name_store("energyGeneralSrc(location)",L_energyGeneralSrclocation_) ;
#line 364 "mms_funcs.loci"
       input("location,stime,numSpecies,meanFlowSrcFunc") ;
#line 364 "mms_funcs.loci"
       output("energyGeneralSrc(location)") ;
#line 364 "mms_funcs.loci"
       constraint("MMSTimeDependent,location") ;
#line 364 "mms_funcs.loci"
    }
#line 364 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 365 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  varmap ["t"] = L_stime_[_e_];
  Loci ::expression ::exprList ::const_iterator li ;
  li = (L_meanFlowSrcFunc_[_e_])->expr_list .begin () ;
  
  for (int i =0;i <L_numSpecies_[_e_]+3;++i ) {
    li ++ ;
  }

  // energy source term in Kg J/(m^3 sec)
  L_energyGeneralSrclocation_[_e_]= (*li )->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 379 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 379 "mms_funcs.loci"
    }
#line 379 "mms_funcs.loci"
} ;
#line 379 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs023_1278521288m330> register_file_mms_funcs023_1278521288m330 ;
#line 379 "mms_funcs.loci"
}
#line 379 "mms_funcs.loci"


// $type MMS_tkeFunc param<Loci::exprP> 
// $type MMS_omegaFunc param<Loci::exprP> 

namespace {class file_mms_funcs024_1278521288m331 : public Loci::singleton_rule {
#line 384 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_operatorDefs_ ; 
#line 384 "mms_funcs.loci"
    Loci::param<Loci::exprP>  L_MMS_tkeFunc_ ; 
#line 384 "mms_funcs.loci"
public:
#line 384 "mms_funcs.loci"
    file_mms_funcs024_1278521288m331() {
#line 384 "mms_funcs.loci"
       name_store("operatorDefs",L_operatorDefs_) ;
#line 384 "mms_funcs.loci"
       name_store("MMS_tkeFunc",L_MMS_tkeFunc_) ;
#line 384 "mms_funcs.loci"
       input("operatorDefs") ;
#line 384 "mms_funcs.loci"
       output("MMS_tkeFunc") ;
#line 384 "mms_funcs.loci"
    }
#line 384 "mms_funcs.loci"
    void compute(const Loci::sequence &seq) { 
  Loci::exprP vars = Loci::expression::create("tke") ;
  (*L_MMS_tkeFunc_)= Loci::substitutionEngine(vars,(*L_operatorDefs_)) ;
}} ;
#line 387 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs024_1278521288m331> register_file_mms_funcs024_1278521288m331 ;
#line 387 "mms_funcs.loci"
}
#line 387 "mms_funcs.loci"


namespace {class file_mms_funcs025_1278521288m331 : public Loci::singleton_rule {
#line 389 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_operatorDefs_ ; 
#line 389 "mms_funcs.loci"
    Loci::param<Loci::exprP>  L_MMS_omegaFunc_ ; 
#line 389 "mms_funcs.loci"
public:
#line 389 "mms_funcs.loci"
    file_mms_funcs025_1278521288m331() {
#line 389 "mms_funcs.loci"
       name_store("operatorDefs",L_operatorDefs_) ;
#line 389 "mms_funcs.loci"
       name_store("MMS_omegaFunc",L_MMS_omegaFunc_) ;
#line 389 "mms_funcs.loci"
       input("operatorDefs") ;
#line 389 "mms_funcs.loci"
       output("MMS_omegaFunc") ;
#line 389 "mms_funcs.loci"
    }
#line 389 "mms_funcs.loci"
    void compute(const Loci::sequence &seq) { 
  Loci::exprP vars = Loci::expression::create("omega") ;
  (*L_MMS_omegaFunc_)= Loci::substitutionEngine(vars,(*L_operatorDefs_)) ;
}} ;
#line 392 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs025_1278521288m331> register_file_mms_funcs025_1278521288m331 ;
#line 392 "mms_funcs.loci"
}
#line 392 "mms_funcs.loci"


namespace {class file_mms_funcs026_1278521288m331 : public Loci::pointwise_rule {
#line 396 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 396 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_MMS_tkeFunc_ ; 
#line 396 "mms_funcs.loci"
    Loci::store<real>  L_mms_klocation_ ; 
#line 396 "mms_funcs.loci"
public:
#line 396 "mms_funcs.loci"
    file_mms_funcs026_1278521288m331() {
#line 396 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 396 "mms_funcs.loci"
       name_store("MMS_tkeFunc",L_MMS_tkeFunc_) ;
#line 396 "mms_funcs.loci"
       name_store("mms_k(location)",L_mms_klocation_) ;
#line 396 "mms_funcs.loci"
       input("location,MMS_tkeFunc") ;
#line 396 "mms_funcs.loci"
       output("mms_k(location)") ;
#line 396 "mms_funcs.loci"
       constraint("MMSSteadyState,location") ;
#line 396 "mms_funcs.loci"
    }
#line 396 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 397 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  L_mms_klocation_[_e_]= (L_MMS_tkeFunc_[_e_])->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 402 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 402 "mms_funcs.loci"
    }
#line 402 "mms_funcs.loci"
} ;
#line 402 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs026_1278521288m331> register_file_mms_funcs026_1278521288m331 ;
#line 402 "mms_funcs.loci"
}
#line 402 "mms_funcs.loci"



namespace {class file_mms_funcs027_1278521288m332 : public Loci::pointwise_rule {
#line 407 "mms_funcs.loci"
    Loci::const_param<streamUns::real>  L_stime_ ; 
#line 407 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 407 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_MMS_tkeFunc_ ; 
#line 407 "mms_funcs.loci"
    Loci::store<real>  L_mms_klocation_ ; 
#line 407 "mms_funcs.loci"
public:
#line 407 "mms_funcs.loci"
    file_mms_funcs027_1278521288m332() {
#line 407 "mms_funcs.loci"
       name_store("stime",L_stime_) ;
#line 407 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 407 "mms_funcs.loci"
       name_store("MMS_tkeFunc",L_MMS_tkeFunc_) ;
#line 407 "mms_funcs.loci"
       name_store("mms_k(location)",L_mms_klocation_) ;
#line 407 "mms_funcs.loci"
       input("location,stime,MMS_tkeFunc") ;
#line 407 "mms_funcs.loci"
       output("mms_k(location)") ;
#line 407 "mms_funcs.loci"
       constraint("MMSTimeDependent,location") ;
#line 407 "mms_funcs.loci"
    }
#line 407 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 408 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  varmap ["t"] = L_stime_[_e_];
  L_mms_klocation_[_e_]= (L_MMS_tkeFunc_[_e_])->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 414 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 414 "mms_funcs.loci"
    }
#line 414 "mms_funcs.loci"
} ;
#line 414 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs027_1278521288m332> register_file_mms_funcs027_1278521288m332 ;
#line 414 "mms_funcs.loci"
}
#line 414 "mms_funcs.loci"


namespace {class file_mms_funcs028_1278521288m332 : public Loci::pointwise_rule {
#line 418 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 418 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_MMS_omegaFunc_ ; 
#line 418 "mms_funcs.loci"
    Loci::store<real>  L_mms_wlocation_ ; 
#line 418 "mms_funcs.loci"
public:
#line 418 "mms_funcs.loci"
    file_mms_funcs028_1278521288m332() {
#line 418 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 418 "mms_funcs.loci"
       name_store("MMS_omegaFunc",L_MMS_omegaFunc_) ;
#line 418 "mms_funcs.loci"
       name_store("mms_w(location)",L_mms_wlocation_) ;
#line 418 "mms_funcs.loci"
       input("location,MMS_omegaFunc") ;
#line 418 "mms_funcs.loci"
       output("mms_w(location)") ;
#line 418 "mms_funcs.loci"
       constraint("MMSSteadyState,location") ;
#line 418 "mms_funcs.loci"
    }
#line 418 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 419 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  L_mms_wlocation_[_e_]= (L_MMS_omegaFunc_[_e_])->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 424 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 424 "mms_funcs.loci"
    }
#line 424 "mms_funcs.loci"
} ;
#line 424 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs028_1278521288m332> register_file_mms_funcs028_1278521288m332 ;
#line 424 "mms_funcs.loci"
}
#line 424 "mms_funcs.loci"



namespace {class file_mms_funcs029_1278521288m333 : public Loci::pointwise_rule {
#line 429 "mms_funcs.loci"
    Loci::const_param<streamUns::real>  L_stime_ ; 
#line 429 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 429 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_MMS_omegaFunc_ ; 
#line 429 "mms_funcs.loci"
    Loci::store<real>  L_mms_wlocation_ ; 
#line 429 "mms_funcs.loci"
public:
#line 429 "mms_funcs.loci"
    file_mms_funcs029_1278521288m333() {
#line 429 "mms_funcs.loci"
       name_store("stime",L_stime_) ;
#line 429 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 429 "mms_funcs.loci"
       name_store("MMS_omegaFunc",L_MMS_omegaFunc_) ;
#line 429 "mms_funcs.loci"
       name_store("mms_w(location)",L_mms_wlocation_) ;
#line 429 "mms_funcs.loci"
       input("location,stime,MMS_omegaFunc") ;
#line 429 "mms_funcs.loci"
       output("mms_w(location)") ;
#line 429 "mms_funcs.loci"
       constraint("MMSTimeDependent,location") ;
#line 429 "mms_funcs.loci"
    }
#line 429 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 430 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  varmap ["t"] = L_stime_[_e_];
  L_mms_wlocation_[_e_]= (L_MMS_omegaFunc_[_e_])->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 436 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 436 "mms_funcs.loci"
    }
#line 436 "mms_funcs.loci"
} ;
#line 436 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs029_1278521288m333> register_file_mms_funcs029_1278521288m333 ;
#line 436 "mms_funcs.loci"
}
#line 436 "mms_funcs.loci"



// $type tkeSrcFunc param<Loci::exprP> 

namespace {class file_mms_funcs030_1278521288m333 : public Loci::singleton_rule {
#line 441 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_operatorDefs_ ; 
#line 441 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_tkeEqnOperator_ ; 
#line 441 "mms_funcs.loci"
    Loci::param<Loci::exprP>  L_tkeSrcFunc_ ; 
#line 441 "mms_funcs.loci"
public:
#line 441 "mms_funcs.loci"
    file_mms_funcs030_1278521288m333() {
#line 441 "mms_funcs.loci"
       name_store("operatorDefs",L_operatorDefs_) ;
#line 441 "mms_funcs.loci"
       name_store("tkeEqnOperator",L_tkeEqnOperator_) ;
#line 441 "mms_funcs.loci"
       name_store("tkeSrcFunc",L_tkeSrcFunc_) ;
#line 441 "mms_funcs.loci"
       input("tkeEqnOperator,operatorDefs") ;
#line 441 "mms_funcs.loci"
       output("tkeSrcFunc") ;
#line 441 "mms_funcs.loci"
    }
#line 441 "mms_funcs.loci"
    void compute(const Loci::sequence &seq) { 
  Loci::exprP sub = Loci::substitutionEngine((*L_tkeEqnOperator_),(*L_operatorDefs_)) ;
  sub = sub->simplify() ; // Simplify expression before evaluating 
  // derivative expressions
  sub = sub->symbolic_eval() ; // evaluate symbolic operators (derivatives)
  sub = sub->simplify() ; // simplify result
  (*L_tkeSrcFunc_)= sub ;
  if(Loci::MPI_rank == 0) {
    ofstream ofile("tkeOperator.dat",ios::out) ;
    ofile << "tke Model Operator:" << endl ;
    ofile << (*L_tkeEqnOperator_)<< endl << endl ;
    ofile << "tke Source Term:" << endl ;
    ofile << sub << endl << endl ;
  }
}} ;
#line 455 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs030_1278521288m333> register_file_mms_funcs030_1278521288m333 ;
#line 455 "mms_funcs.loci"
}
#line 455 "mms_funcs.loci"


namespace {class file_mms_funcs031_1278521288m333 : public Loci::pointwise_rule {
#line 458 "mms_funcs.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 458 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 458 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_tkeSrcFunc_ ; 
#line 458 "mms_funcs.loci"
    Loci::store<real>  L_rkGeneralSrclocation_ ; 
#line 458 "mms_funcs.loci"
public:
#line 458 "mms_funcs.loci"
    file_mms_funcs031_1278521288m333() {
#line 458 "mms_funcs.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 458 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 458 "mms_funcs.loci"
       name_store("tkeSrcFunc",L_tkeSrcFunc_) ;
#line 458 "mms_funcs.loci"
       name_store("rkGeneralSrc(location)",L_rkGeneralSrclocation_) ;
#line 458 "mms_funcs.loci"
       input("location,numSpecies,tkeSrcFunc") ;
#line 458 "mms_funcs.loci"
       output("rkGeneralSrc(location)") ;
#line 458 "mms_funcs.loci"
       constraint("MMSSteadyState,location") ;
#line 458 "mms_funcs.loci"
    }
#line 458 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 459 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  L_rkGeneralSrclocation_[_e_]= (L_tkeSrcFunc_[_e_])->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 464 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 464 "mms_funcs.loci"
    }
#line 464 "mms_funcs.loci"
} ;
#line 464 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs031_1278521288m333> register_file_mms_funcs031_1278521288m333 ;
#line 464 "mms_funcs.loci"
}
#line 464 "mms_funcs.loci"


namespace {class file_mms_funcs032_1278521288m334 : public Loci::pointwise_rule {
#line 467 "mms_funcs.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 467 "mms_funcs.loci"
    Loci::const_param<streamUns::real>  L_stime_ ; 
#line 467 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 467 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_tkeSrcFunc_ ; 
#line 467 "mms_funcs.loci"
    Loci::store<real>  L_rkGeneralSrclocation_ ; 
#line 467 "mms_funcs.loci"
public:
#line 467 "mms_funcs.loci"
    file_mms_funcs032_1278521288m334() {
#line 467 "mms_funcs.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 467 "mms_funcs.loci"
       name_store("stime",L_stime_) ;
#line 467 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 467 "mms_funcs.loci"
       name_store("tkeSrcFunc",L_tkeSrcFunc_) ;
#line 467 "mms_funcs.loci"
       name_store("rkGeneralSrc(location)",L_rkGeneralSrclocation_) ;
#line 467 "mms_funcs.loci"
       input("location,stime,numSpecies,tkeSrcFunc") ;
#line 467 "mms_funcs.loci"
       output("rkGeneralSrc(location)") ;
#line 467 "mms_funcs.loci"
       constraint("MMSTimeDependent,location") ;
#line 467 "mms_funcs.loci"
    }
#line 467 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 468 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  varmap ["t"] = L_stime_[_e_];
  L_rkGeneralSrclocation_[_e_]= (L_tkeSrcFunc_[_e_])->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 474 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 474 "mms_funcs.loci"
    }
#line 474 "mms_funcs.loci"
} ;
#line 474 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs032_1278521288m334> register_file_mms_funcs032_1278521288m334 ;
#line 474 "mms_funcs.loci"
}
#line 474 "mms_funcs.loci"


// $type omegaSrcFunc param<Loci::exprP> 

namespace {class file_mms_funcs033_1278521288m334 : public Loci::singleton_rule {
#line 478 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_operatorDefs_ ; 
#line 478 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_omegaEqnOperator_ ; 
#line 478 "mms_funcs.loci"
    Loci::param<Loci::exprP>  L_omegaSrcFunc_ ; 
#line 478 "mms_funcs.loci"
public:
#line 478 "mms_funcs.loci"
    file_mms_funcs033_1278521288m334() {
#line 478 "mms_funcs.loci"
       name_store("operatorDefs",L_operatorDefs_) ;
#line 478 "mms_funcs.loci"
       name_store("omegaEqnOperator",L_omegaEqnOperator_) ;
#line 478 "mms_funcs.loci"
       name_store("omegaSrcFunc",L_omegaSrcFunc_) ;
#line 478 "mms_funcs.loci"
       input("omegaEqnOperator,operatorDefs") ;
#line 478 "mms_funcs.loci"
       output("omegaSrcFunc") ;
#line 478 "mms_funcs.loci"
    }
#line 478 "mms_funcs.loci"
    void compute(const Loci::sequence &seq) { 
  Loci::exprP sub = Loci::substitutionEngine((*L_omegaEqnOperator_),(*L_operatorDefs_)) ;
  sub = sub->simplify() ; // Simplify expression before evaluating 
  // derivative expressions
  sub = sub->symbolic_eval() ; // evaluate symbolic operators (derivatives)
  sub = sub->simplify() ; // simplify result
  (*L_omegaSrcFunc_)= sub ;
  if(Loci::MPI_rank == 0) {
    ofstream ofile("omegaOperator.dat",ios::out) ;
    ofile << "omega Model Operator:" << endl ;
    ofile << (*L_omegaEqnOperator_)<< endl << endl ;
    ofile << "omega Source Term:" << endl ;
    ofile << sub << endl << endl ;
  }
}} ;
#line 492 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs033_1278521288m334> register_file_mms_funcs033_1278521288m334 ;
#line 492 "mms_funcs.loci"
}
#line 492 "mms_funcs.loci"


namespace {class file_mms_funcs034_1278521288m335 : public Loci::pointwise_rule {
#line 495 "mms_funcs.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 495 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 495 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_omegaSrcFunc_ ; 
#line 495 "mms_funcs.loci"
    Loci::store<real>  L_rwGeneralSrclocation_ ; 
#line 495 "mms_funcs.loci"
public:
#line 495 "mms_funcs.loci"
    file_mms_funcs034_1278521288m335() {
#line 495 "mms_funcs.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 495 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 495 "mms_funcs.loci"
       name_store("omegaSrcFunc",L_omegaSrcFunc_) ;
#line 495 "mms_funcs.loci"
       name_store("rwGeneralSrc(location)",L_rwGeneralSrclocation_) ;
#line 495 "mms_funcs.loci"
       input("location,numSpecies,omegaSrcFunc") ;
#line 495 "mms_funcs.loci"
       output("rwGeneralSrc(location)") ;
#line 495 "mms_funcs.loci"
       constraint("MMSSteadyState,location") ;
#line 495 "mms_funcs.loci"
    }
#line 495 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 496 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  L_rwGeneralSrclocation_[_e_]= (L_omegaSrcFunc_[_e_])->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 501 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 501 "mms_funcs.loci"
    }
#line 501 "mms_funcs.loci"
} ;
#line 501 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs034_1278521288m335> register_file_mms_funcs034_1278521288m335 ;
#line 501 "mms_funcs.loci"
}
#line 501 "mms_funcs.loci"


namespace {class file_mms_funcs035_1278521288m335 : public Loci::pointwise_rule {
#line 504 "mms_funcs.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 504 "mms_funcs.loci"
    Loci::const_param<streamUns::real>  L_stime_ ; 
#line 504 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 504 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_omegaSrcFunc_ ; 
#line 504 "mms_funcs.loci"
    Loci::store<real>  L_rwGeneralSrclocation_ ; 
#line 504 "mms_funcs.loci"
public:
#line 504 "mms_funcs.loci"
    file_mms_funcs035_1278521288m335() {
#line 504 "mms_funcs.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 504 "mms_funcs.loci"
       name_store("stime",L_stime_) ;
#line 504 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 504 "mms_funcs.loci"
       name_store("omegaSrcFunc",L_omegaSrcFunc_) ;
#line 504 "mms_funcs.loci"
       name_store("rwGeneralSrc(location)",L_rwGeneralSrclocation_) ;
#line 504 "mms_funcs.loci"
       input("location,stime,numSpecies,omegaSrcFunc") ;
#line 504 "mms_funcs.loci"
       output("rwGeneralSrc(location)") ;
#line 504 "mms_funcs.loci"
       constraint("MMSTimeDependent,location") ;
#line 504 "mms_funcs.loci"
    }
#line 504 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 505 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  varmap ["t"] = L_stime_[_e_];
  L_rwGeneralSrclocation_[_e_]= (L_omegaSrcFunc_[_e_])->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 511 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 511 "mms_funcs.loci"
    }
#line 511 "mms_funcs.loci"
} ;
#line 511 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs035_1278521288m335> register_file_mms_funcs035_1278521288m335 ;
#line 511 "mms_funcs.loci"
}
#line 511 "mms_funcs.loci"


// $type MMS_nu_tFunc param<Loci::exprP> 

namespace {class file_mms_funcs036_1278521288m335 : public Loci::singleton_rule {
#line 515 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_operatorDefs_ ; 
#line 515 "mms_funcs.loci"
    Loci::param<Loci::exprP>  L_MMS_nu_tFunc_ ; 
#line 515 "mms_funcs.loci"
public:
#line 515 "mms_funcs.loci"
    file_mms_funcs036_1278521288m335() {
#line 515 "mms_funcs.loci"
       name_store("operatorDefs",L_operatorDefs_) ;
#line 515 "mms_funcs.loci"
       name_store("MMS_nu_tFunc",L_MMS_nu_tFunc_) ;
#line 515 "mms_funcs.loci"
       input("operatorDefs") ;
#line 515 "mms_funcs.loci"
       output("MMS_nu_tFunc") ;
#line 515 "mms_funcs.loci"
    }
#line 515 "mms_funcs.loci"
    void compute(const Loci::sequence &seq) { 
  Loci::exprP vars = Loci::expression::create("nu_t") ;
  (*L_MMS_nu_tFunc_)= Loci::substitutionEngine(vars,(*L_operatorDefs_)) ;
}} ;
#line 518 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs036_1278521288m335> register_file_mms_funcs036_1278521288m335 ;
#line 518 "mms_funcs.loci"
}
#line 518 "mms_funcs.loci"


namespace {class file_mms_funcs037_1278521288m336 : public Loci::pointwise_rule {
#line 522 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 522 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_MMS_nu_tFunc_ ; 
#line 522 "mms_funcs.loci"
    Loci::store<real>  L_mms_nu_tlocation_ ; 
#line 522 "mms_funcs.loci"
public:
#line 522 "mms_funcs.loci"
    file_mms_funcs037_1278521288m336() {
#line 522 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 522 "mms_funcs.loci"
       name_store("MMS_nu_tFunc",L_MMS_nu_tFunc_) ;
#line 522 "mms_funcs.loci"
       name_store("mms_nu_t(location)",L_mms_nu_tlocation_) ;
#line 522 "mms_funcs.loci"
       input("location,MMS_nu_tFunc") ;
#line 522 "mms_funcs.loci"
       output("mms_nu_t(location)") ;
#line 522 "mms_funcs.loci"
       constraint("MMSSteadyState,location") ;
#line 522 "mms_funcs.loci"
    }
#line 522 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 523 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  L_mms_nu_tlocation_[_e_]= (L_MMS_nu_tFunc_[_e_])->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 528 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 528 "mms_funcs.loci"
    }
#line 528 "mms_funcs.loci"
} ;
#line 528 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs037_1278521288m336> register_file_mms_funcs037_1278521288m336 ;
#line 528 "mms_funcs.loci"
}
#line 528 "mms_funcs.loci"



namespace {class file_mms_funcs038_1278521288m336 : public Loci::pointwise_rule {
#line 533 "mms_funcs.loci"
    Loci::const_param<streamUns::real>  L_stime_ ; 
#line 533 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 533 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_MMS_nu_tFunc_ ; 
#line 533 "mms_funcs.loci"
    Loci::store<real>  L_mms_nu_tlocation_ ; 
#line 533 "mms_funcs.loci"
public:
#line 533 "mms_funcs.loci"
    file_mms_funcs038_1278521288m336() {
#line 533 "mms_funcs.loci"
       name_store("stime",L_stime_) ;
#line 533 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 533 "mms_funcs.loci"
       name_store("MMS_nu_tFunc",L_MMS_nu_tFunc_) ;
#line 533 "mms_funcs.loci"
       name_store("mms_nu_t(location)",L_mms_nu_tlocation_) ;
#line 533 "mms_funcs.loci"
       input("location,stime,MMS_nu_tFunc") ;
#line 533 "mms_funcs.loci"
       output("mms_nu_t(location)") ;
#line 533 "mms_funcs.loci"
       constraint("MMSTimeDependent,location") ;
#line 533 "mms_funcs.loci"
    }
#line 533 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 534 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  varmap ["t"] = L_stime_[_e_];
  L_mms_nu_tlocation_[_e_]= (L_MMS_nu_tFunc_[_e_])->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 540 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 540 "mms_funcs.loci"
    }
#line 540 "mms_funcs.loci"
} ;
#line 540 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs038_1278521288m336> register_file_mms_funcs038_1278521288m336 ;
#line 540 "mms_funcs.loci"
}
#line 540 "mms_funcs.loci"



// $type SASrcFunc param<Loci::exprP> 

namespace {class file_mms_funcs039_1278521288m337 : public Loci::singleton_rule {
#line 545 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_operatorDefs_ ; 
#line 545 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_SAOperator_ ; 
#line 545 "mms_funcs.loci"
    Loci::param<Loci::exprP>  L_SASrcFunc_ ; 
#line 545 "mms_funcs.loci"
public:
#line 545 "mms_funcs.loci"
    file_mms_funcs039_1278521288m337() {
#line 545 "mms_funcs.loci"
       name_store("operatorDefs",L_operatorDefs_) ;
#line 545 "mms_funcs.loci"
       name_store("SAOperator",L_SAOperator_) ;
#line 545 "mms_funcs.loci"
       name_store("SASrcFunc",L_SASrcFunc_) ;
#line 545 "mms_funcs.loci"
       input("SAOperator,operatorDefs") ;
#line 545 "mms_funcs.loci"
       output("SASrcFunc") ;
#line 545 "mms_funcs.loci"
    }
#line 545 "mms_funcs.loci"
    void compute(const Loci::sequence &seq) { 
  Loci::exprP sub = Loci::substitutionEngine((*L_SAOperator_),(*L_operatorDefs_)) ;
  sub = sub->simplify() ; // Simplify expression before evaluating 
  // derivative expressions
  sub = sub->symbolic_eval() ; // evaluate symbolic operators (derivatives)
  sub = sub->simplify() ; // simplify result
  (*L_SASrcFunc_)= sub ;
  if(Loci::MPI_rank == 0) {
    ofstream ofile("SAOperator.dat",ios::out) ;
    ofile << "SA Model Operator:" << endl ;
    ofile << (*L_SAOperator_)<< endl << endl ;
    ofile << "SA Source Term:" << endl ;
    ofile << sub << endl << endl ;
  }
}} ;
#line 559 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs039_1278521288m337> register_file_mms_funcs039_1278521288m337 ;
#line 559 "mms_funcs.loci"
}
#line 559 "mms_funcs.loci"


namespace {class file_mms_funcs040_1278521288m337 : public Loci::pointwise_rule {
#line 562 "mms_funcs.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 562 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 562 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_SASrcFunc_ ; 
#line 562 "mms_funcs.loci"
    Loci::store<real>  L_SAGeneralSrclocation_ ; 
#line 562 "mms_funcs.loci"
public:
#line 562 "mms_funcs.loci"
    file_mms_funcs040_1278521288m337() {
#line 562 "mms_funcs.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 562 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 562 "mms_funcs.loci"
       name_store("SASrcFunc",L_SASrcFunc_) ;
#line 562 "mms_funcs.loci"
       name_store("SAGeneralSrc(location)",L_SAGeneralSrclocation_) ;
#line 562 "mms_funcs.loci"
       input("location,numSpecies,SASrcFunc") ;
#line 562 "mms_funcs.loci"
       output("SAGeneralSrc(location)") ;
#line 562 "mms_funcs.loci"
       constraint("MMSSteadyState,location") ;
#line 562 "mms_funcs.loci"
    }
#line 562 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 563 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  L_SAGeneralSrclocation_[_e_]= (L_SASrcFunc_[_e_])->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 568 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 568 "mms_funcs.loci"
    }
#line 568 "mms_funcs.loci"
} ;
#line 568 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs040_1278521288m337> register_file_mms_funcs040_1278521288m337 ;
#line 568 "mms_funcs.loci"
}
#line 568 "mms_funcs.loci"


namespace {class file_mms_funcs041_1278521288m337 : public Loci::pointwise_rule {
#line 571 "mms_funcs.loci"
    Loci::const_param<int>  L_numSpecies_ ; 
#line 571 "mms_funcs.loci"
    Loci::const_param<streamUns::real>  L_stime_ ; 
#line 571 "mms_funcs.loci"
    Loci::const_store<vect3d>  L_location_ ; 
#line 571 "mms_funcs.loci"
    Loci::const_param<Loci::exprP>  L_SASrcFunc_ ; 
#line 571 "mms_funcs.loci"
    Loci::store<real>  L_SAGeneralSrclocation_ ; 
#line 571 "mms_funcs.loci"
public:
#line 571 "mms_funcs.loci"
    file_mms_funcs041_1278521288m337() {
#line 571 "mms_funcs.loci"
       name_store("numSpecies",L_numSpecies_) ;
#line 571 "mms_funcs.loci"
       name_store("stime",L_stime_) ;
#line 571 "mms_funcs.loci"
       name_store("location",L_location_) ;
#line 571 "mms_funcs.loci"
       name_store("SASrcFunc",L_SASrcFunc_) ;
#line 571 "mms_funcs.loci"
       name_store("SAGeneralSrc(location)",L_SAGeneralSrclocation_) ;
#line 571 "mms_funcs.loci"
       input("location,stime,numSpecies,SASrcFunc") ;
#line 571 "mms_funcs.loci"
       output("SAGeneralSrc(location)") ;
#line 571 "mms_funcs.loci"
       constraint("MMSTimeDependent,location") ;
#line 571 "mms_funcs.loci"
    }
#line 571 "mms_funcs.loci"
    void calculate(Entity _e_) { 
#line 572 "mms_funcs.loci"
  map <string ,double > varmap ;
  varmap ["x"] = L_location_[_e_].x ;
  varmap ["y"] = L_location_[_e_].y ;
  varmap ["z"] = L_location_[_e_].z ;
  varmap ["t"] = L_stime_[_e_];
  L_SAGeneralSrclocation_[_e_]= (L_SASrcFunc_[_e_])->evaluate (varmap ) ;
}    void compute(const Loci::sequence &seq) { 
#line 578 "mms_funcs.loci"
      do_loop(seq,this) ;
#line 578 "mms_funcs.loci"
    }
#line 578 "mms_funcs.loci"
} ;
#line 578 "mms_funcs.loci"
Loci::register_rule<file_mms_funcs041_1278521288m337> register_file_mms_funcs041_1278521288m337 ;
#line 578 "mms_funcs.loci"
}
#line 578 "mms_funcs.loci"


