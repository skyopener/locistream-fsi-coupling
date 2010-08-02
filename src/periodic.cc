#line 1 "periodic.loci"
// Loci includes.
#include <Loci.h>
using Loci::rigid_transform ;
using Loci::periodic_info ;

// StreamUns includes.
#include "sciTypes.h"
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
#line 8 "periodic.loci"


namespace streamUns {

  // Additional variable types.
  // $type periodicTransform store<rigid_transform> 

  // This rule, in addition to providing the right cell volume for periodic
  // faces, also causes the periodic faces to be included in our contraint
  // "internalFaces" since this constraint is now computed using (cl,cr)->vol
  // instead of (cl,cr)->geom_cells as was done before.
  namespace {class file_periodic000_1278520190m661 : public Loci::pointwise_rule {
#line 19 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 19 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 19 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 19 "periodic.loci"
    Loci::store<Loci::real_t>  L_vol_ ; 
#line 19 "periodic.loci"
public:
#line 19 "periodic.loci"
    file_periodic000_1278520190m661() {
#line 19 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 19 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 19 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 19 "periodic.loci"
       name_store("vol",L_vol_) ;
#line 19 "periodic.loci"
       input("pmap->cl->vol") ;
#line 19 "periodic.loci"
       output("cr->vol") ;
#line 19 "periodic.loci"
    }
#line 19 "periodic.loci"
    void calculate(Entity _e_) { 
#line 20 "periodic.loci"
    L_vol_[L_cr_[_e_]]= L_vol_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 21 "periodic.loci"
      do_loop(seq,this) ;
#line 21 "periodic.loci"
    }
#line 21 "periodic.loci"
} ;
#line 21 "periodic.loci"
Loci::register_rule<file_periodic000_1278520190m661> register_file_periodic000_1278520190m661 ;
#line 21 "periodic.loci"
}
#line 21 "periodic.loci"


  // Cell center.
  namespace {class file_periodic001_1278520190m662 : public Loci::pointwise_rule {
#line 25 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 25 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 25 "periodic.loci"
    Loci::const_Map L_ref_ ; 
#line 25 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 25 "periodic.loci"
    Loci::const_store<rigid_transform>  L_periodicTransform_ ; 
#line 25 "periodic.loci"
    Loci::store<Loci::vector3d<Loci::real_t> >  L_cellcenter_ ; 
#line 25 "periodic.loci"
public:
#line 25 "periodic.loci"
    file_periodic001_1278520190m662() {
#line 25 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 25 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 25 "periodic.loci"
       name_store("ref",L_ref_) ;
#line 25 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 25 "periodic.loci"
       name_store("cellcenter",L_cellcenter_) ;
#line 25 "periodic.loci"
       name_store("periodicTransform",L_periodicTransform_) ;
#line 25 "periodic.loci"
       input("pmap->cl->cellcenter,pmap->ref->  periodicTransform") ;
#line 25 "periodic.loci"
       output("cr->cellcenter") ;
#line 25 "periodic.loci"
       constraint("periodicFaces") ;
#line 25 "periodic.loci"
    }
#line 25 "periodic.loci"
    void calculate(Entity _e_) { 
#line 26 "periodic.loci"
    const rigid_transform &frame = L_periodicTransform_[L_ref_[L_pmap_[_e_]]];
    L_cellcenter_[L_cr_[_e_]]= frame .transform (L_cellcenter_[L_cl_[L_pmap_[_e_]]]) ;
  }    void compute(const Loci::sequence &seq) { 
#line 28 "periodic.loci"
      do_loop(seq,this) ;
#line 28 "periodic.loci"
    }
#line 28 "periodic.loci"
} ;
#line 28 "periodic.loci"
Loci::register_rule<file_periodic001_1278520190m662> register_file_periodic001_1278520190m662 ;
#line 28 "periodic.loci"
}
#line 28 "periodic.loci"


  // Laminar viscosity.
  // $type laminarViscosity store<real> 
  namespace {class file_periodic002_1278520190m663 : public Loci::pointwise_rule {
#line 32 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 32 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 32 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 32 "periodic.loci"
    Loci::store<real>  L_laminarViscosity_ ; 
#line 32 "periodic.loci"
public:
#line 32 "periodic.loci"
    file_periodic002_1278520190m663() {
#line 32 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 32 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 32 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 32 "periodic.loci"
       name_store("laminarViscosity",L_laminarViscosity_) ;
#line 32 "periodic.loci"
       input("pmap->cl->laminarViscosity") ;
#line 32 "periodic.loci"
       output("cr->laminarViscosity") ;
#line 32 "periodic.loci"
    }
#line 32 "periodic.loci"
    void calculate(Entity _e_) { 
#line 33 "periodic.loci"
    L_laminarViscosity_[L_cr_[_e_]]=L_laminarViscosity_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 34 "periodic.loci"
      do_loop(seq,this) ;
#line 34 "periodic.loci"
    }
#line 34 "periodic.loci"
} ;
#line 34 "periodic.loci"
Loci::register_rule<file_periodic002_1278520190m663> register_file_periodic002_1278520190m663 ;
#line 34 "periodic.loci"
}
#line 34 "periodic.loci"


  // Eddy viscosity. Only need this rule for laminar flow, since eddyViscosity
  // is a derived variable for turbulent flows.
  // $type eddyViscosity store<real> 
  namespace {class file_periodic003_1278520190m664 : public Loci::pointwise_rule {
#line 40 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 40 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 40 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 40 "periodic.loci"
    Loci::store<real>  L_eddyViscosity_ ; 
#line 40 "periodic.loci"
public:
#line 40 "periodic.loci"
    file_periodic003_1278520190m664() {
#line 40 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 40 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 40 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 40 "periodic.loci"
       name_store("eddyViscosity",L_eddyViscosity_) ;
#line 40 "periodic.loci"
       input("pmap->cl->eddyViscosity") ;
#line 40 "periodic.loci"
       output("cr->eddyViscosity") ;
#line 40 "periodic.loci"
       constraint("laminarFlow,periodicFaces") ;
#line 40 "periodic.loci"
    }
#line 40 "periodic.loci"
    void calculate(Entity _e_) { 
#line 41 "periodic.loci"
    L_eddyViscosity_[L_cr_[_e_]]= L_eddyViscosity_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 42 "periodic.loci"
      do_loop(seq,this) ;
#line 42 "periodic.loci"
    }
#line 42 "periodic.loci"
} ;
#line 42 "periodic.loci"
Loci::register_rule<file_periodic003_1278520190m664> register_file_periodic003_1278520190m664 ;
#line 42 "periodic.loci"
}
#line 42 "periodic.loci"


  // Speed of sound.
  // $type soundSpeed store<real> 
  namespace {class file_periodic004_1278520190m664 : public Loci::pointwise_rule {
#line 46 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 46 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 46 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 46 "periodic.loci"
    Loci::store<real>  L_soundSpeed_ ; 
#line 46 "periodic.loci"
public:
#line 46 "periodic.loci"
    file_periodic004_1278520190m664() {
#line 46 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 46 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 46 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 46 "periodic.loci"
       name_store("soundSpeed",L_soundSpeed_) ;
#line 46 "periodic.loci"
       input("pmap->cl->soundSpeed") ;
#line 46 "periodic.loci"
       output("cr->soundSpeed") ;
#line 46 "periodic.loci"
    }
#line 46 "periodic.loci"
    void calculate(Entity _e_) { 
#line 47 "periodic.loci"
    L_soundSpeed_[L_cr_[_e_]]=L_soundSpeed_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 48 "periodic.loci"
      do_loop(seq,this) ;
#line 48 "periodic.loci"
    }
#line 48 "periodic.loci"
} ;
#line 48 "periodic.loci"
Loci::register_rule<file_periodic004_1278520190m664> register_file_periodic004_1278520190m664 ;
#line 48 "periodic.loci"
}
#line 48 "periodic.loci"


  // Specific heat.
  // $type cp store<real> 
  namespace {class file_periodic005_1278520190m665 : public Loci::pointwise_rule {
#line 52 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 52 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 52 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 52 "periodic.loci"
    Loci::store<real>  L_cp_ ; 
#line 52 "periodic.loci"
public:
#line 52 "periodic.loci"
    file_periodic005_1278520190m665() {
#line 52 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 52 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 52 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 52 "periodic.loci"
       name_store("cp",L_cp_) ;
#line 52 "periodic.loci"
       input("pmap->cl->cp") ;
#line 52 "periodic.loci"
       output("cr->cp") ;
#line 52 "periodic.loci"
    }
#line 52 "periodic.loci"
    void calculate(Entity _e_) { 
#line 53 "periodic.loci"
    L_cp_[L_cr_[_e_]]=L_cp_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 54 "periodic.loci"
      do_loop(seq,this) ;
#line 54 "periodic.loci"
    }
#line 54 "periodic.loci"
} ;
#line 54 "periodic.loci"
Loci::register_rule<file_periodic005_1278520190m665> register_file_periodic005_1278520190m665 ;
#line 54 "periodic.loci"
}
#line 54 "periodic.loci"


  // Current species enthalpy for cells.
  // $type yCurrEnthalpyCell store<real> 
  namespace {class file_periodic006_1278520190m666 : public Loci::pointwise_rule {
#line 58 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 58 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 58 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 58 "periodic.loci"
    Loci::store<real>  L_yCurrEnthalpyCell_ ; 
#line 58 "periodic.loci"
public:
#line 58 "periodic.loci"
    file_periodic006_1278520190m666() {
#line 58 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 58 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 58 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 58 "periodic.loci"
       name_store("yCurrEnthalpyCell",L_yCurrEnthalpyCell_) ;
#line 58 "periodic.loci"
       input("pmap->cl->yCurrEnthalpyCell") ;
#line 58 "periodic.loci"
       output("cr->yCurrEnthalpyCell") ;
#line 58 "periodic.loci"
    }
#line 58 "periodic.loci"
    void calculate(Entity _e_) { 
#line 59 "periodic.loci"
    L_yCurrEnthalpyCell_[L_cr_[_e_]]=L_yCurrEnthalpyCell_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 60 "periodic.loci"
      do_loop(seq,this) ;
#line 60 "periodic.loci"
    }
#line 60 "periodic.loci"
} ;
#line 60 "periodic.loci"
Loci::register_rule<file_periodic006_1278520190m666> register_file_periodic006_1278520190m666 ;
#line 60 "periodic.loci"
}
#line 60 "periodic.loci"


  // Density initial condition.
  // $type rho_ic store<real> 
  namespace {class file_periodic007_1278520190m666 : public Loci::pointwise_rule {
#line 64 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 64 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 64 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 64 "periodic.loci"
    Loci::store<real>  L_rho_ic_ ; 
#line 64 "periodic.loci"
public:
#line 64 "periodic.loci"
    file_periodic007_1278520190m666() {
#line 64 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 64 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 64 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 64 "periodic.loci"
       name_store("rho_ic",L_rho_ic_) ;
#line 64 "periodic.loci"
       input("pmap->cl->rho_ic") ;
#line 64 "periodic.loci"
       output("cr->rho_ic") ;
#line 64 "periodic.loci"
    }
#line 64 "periodic.loci"
    void calculate(Entity _e_) { 
#line 65 "periodic.loci"
    L_rho_ic_[L_cr_[_e_]]=L_rho_ic_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 66 "periodic.loci"
      do_loop(seq,this) ;
#line 66 "periodic.loci"
    }
#line 66 "periodic.loci"
} ;
#line 66 "periodic.loci"
Loci::register_rule<file_periodic007_1278520190m666> register_file_periodic007_1278520190m666 ;
#line 66 "periodic.loci"
}
#line 66 "periodic.loci"


  // Density.
  // $type rho store<real> 
  namespace {class file_periodic008_1278520190m667 : public Loci::pointwise_rule {
#line 70 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 70 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 70 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 70 "periodic.loci"
    Loci::store<real>  L_rho_ ; 
#line 70 "periodic.loci"
public:
#line 70 "periodic.loci"
    file_periodic008_1278520190m667() {
#line 70 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 70 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 70 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 70 "periodic.loci"
       name_store("rho",L_rho_) ;
#line 70 "periodic.loci"
       input("pmap->cl->rho") ;
#line 70 "periodic.loci"
       output("cr->rho") ;
#line 70 "periodic.loci"
    }
#line 70 "periodic.loci"
    void calculate(Entity _e_) { 
#line 71 "periodic.loci"
    L_rho_[L_cr_[_e_]]=L_rho_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 72 "periodic.loci"
      do_loop(seq,this) ;
#line 72 "periodic.loci"
    }
#line 72 "periodic.loci"
} ;
#line 72 "periodic.loci"
Loci::register_rule<file_periodic008_1278520190m667> register_file_periodic008_1278520190m667 ;
#line 72 "periodic.loci"
}
#line 72 "periodic.loci"


  // Velocity initial condition.
  // $type v_ic store<vect3d> 
  namespace {class file_periodic009_1278520190m667 : public Loci::pointwise_rule {
#line 77 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 77 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 77 "periodic.loci"
    Loci::const_Map L_ref_ ; 
#line 77 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 77 "periodic.loci"
    Loci::const_store<rigid_transform>  L_periodicTransform_ ; 
#line 77 "periodic.loci"
    Loci::store<vect3d>  L_v_ic_ ; 
#line 77 "periodic.loci"
public:
#line 77 "periodic.loci"
    file_periodic009_1278520190m667() {
#line 77 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 77 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 77 "periodic.loci"
       name_store("ref",L_ref_) ;
#line 77 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 77 "periodic.loci"
       name_store("periodicTransform",L_periodicTransform_) ;
#line 77 "periodic.loci"
       name_store("v_ic",L_v_ic_) ;
#line 77 "periodic.loci"
       input("pmap->cl->v_ic,pmap->ref->periodicTransform") ;
#line 77 "periodic.loci"
       output("cr->v_ic") ;
#line 77 "periodic.loci"
       constraint("periodicFaces") ;
#line 77 "periodic.loci"
    }
#line 77 "periodic.loci"
    void calculate(Entity _e_) { 
#line 78 "periodic.loci"
    const rigid_transform &frame = L_periodicTransform_[L_ref_[L_pmap_[_e_]]];
    L_v_ic_[L_cr_[_e_]]=frame .rotate_vec (L_v_ic_[L_cl_[L_pmap_[_e_]]]) ;
  }    void compute(const Loci::sequence &seq) { 
#line 80 "periodic.loci"
      do_loop(seq,this) ;
#line 80 "periodic.loci"
    }
#line 80 "periodic.loci"
} ;
#line 80 "periodic.loci"
Loci::register_rule<file_periodic009_1278520190m667> register_file_periodic009_1278520190m667 ;
#line 80 "periodic.loci"
}
#line 80 "periodic.loci"


  // Velocity.
  // $type v store<vect3d> 
  namespace {class file_periodic010_1278520190m668 : public Loci::pointwise_rule {
#line 85 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 85 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 85 "periodic.loci"
    Loci::const_Map L_ref_ ; 
#line 85 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 85 "periodic.loci"
    Loci::const_store<rigid_transform>  L_periodicTransform_ ; 
#line 85 "periodic.loci"
    Loci::store<vect3d>  L_v_ ; 
#line 85 "periodic.loci"
public:
#line 85 "periodic.loci"
    file_periodic010_1278520190m668() {
#line 85 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 85 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 85 "periodic.loci"
       name_store("ref",L_ref_) ;
#line 85 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 85 "periodic.loci"
       name_store("periodicTransform",L_periodicTransform_) ;
#line 85 "periodic.loci"
       name_store("v",L_v_) ;
#line 85 "periodic.loci"
       input("pmap->cl->v,pmap->ref->periodicTransform") ;
#line 85 "periodic.loci"
       output("cr->v") ;
#line 85 "periodic.loci"
       constraint("periodicFaces") ;
#line 85 "periodic.loci"
    }
#line 85 "periodic.loci"
    void calculate(Entity _e_) { 
#line 86 "periodic.loci"
    const rigid_transform &frame = L_periodicTransform_[L_ref_[L_pmap_[_e_]]];
    L_v_[L_cr_[_e_]]= frame .rotate_vec (L_v_[L_cl_[L_pmap_[_e_]]]) ;
  }    void compute(const Loci::sequence &seq) { 
#line 88 "periodic.loci"
      do_loop(seq,this) ;
#line 88 "periodic.loci"
    }
#line 88 "periodic.loci"
} ;
#line 88 "periodic.loci"
Loci::register_rule<file_periodic010_1278520190m668> register_file_periodic010_1278520190m668 ;
#line 88 "periodic.loci"
}
#line 88 "periodic.loci"


  // New velocity.
  // $type vStar store<vect3d> 
  namespace {class file_periodic011_1278520190m668 : public Loci::pointwise_rule {
#line 93 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 93 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 93 "periodic.loci"
    Loci::const_Map L_ref_ ; 
#line 93 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 93 "periodic.loci"
    Loci::const_store<rigid_transform>  L_periodicTransform_ ; 
#line 93 "periodic.loci"
    Loci::store<vect3d>  L_vStar_ ; 
#line 93 "periodic.loci"
public:
#line 93 "periodic.loci"
    file_periodic011_1278520190m668() {
#line 93 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 93 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 93 "periodic.loci"
       name_store("ref",L_ref_) ;
#line 93 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 93 "periodic.loci"
       name_store("periodicTransform",L_periodicTransform_) ;
#line 93 "periodic.loci"
       name_store("vStar",L_vStar_) ;
#line 93 "periodic.loci"
       input("pmap->cl->vStar,pmap->ref->periodicTransform") ;
#line 93 "periodic.loci"
       output("cr->vStar") ;
#line 93 "periodic.loci"
       constraint("periodicFaces") ;
#line 93 "periodic.loci"
    }
#line 93 "periodic.loci"
    void calculate(Entity _e_) { 
#line 94 "periodic.loci"
    const rigid_transform &frame = L_periodicTransform_[L_ref_[L_pmap_[_e_]]];
    L_vStar_[L_cr_[_e_]]= frame .rotate_vec (L_vStar_[L_cl_[L_pmap_[_e_]]]) ;
  }    void compute(const Loci::sequence &seq) { 
#line 96 "periodic.loci"
      do_loop(seq,this) ;
#line 96 "periodic.loci"
    }
#line 96 "periodic.loci"
} ;
#line 96 "periodic.loci"
Loci::register_rule<file_periodic011_1278520190m668> register_file_periodic011_1278520190m668 ;
#line 96 "periodic.loci"
}
#line 96 "periodic.loci"


  // Pressure initial condition.
  // $type p_ic store<real> 
  namespace {class file_periodic012_1278520190m668 : public Loci::pointwise_rule {
#line 100 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 100 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 100 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 100 "periodic.loci"
    Loci::store<real>  L_p_ic_ ; 
#line 100 "periodic.loci"
public:
#line 100 "periodic.loci"
    file_periodic012_1278520190m668() {
#line 100 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 100 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 100 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 100 "periodic.loci"
       name_store("p_ic",L_p_ic_) ;
#line 100 "periodic.loci"
       input("pmap->cl->p_ic") ;
#line 100 "periodic.loci"
       output("cr->p_ic") ;
#line 100 "periodic.loci"
    }
#line 100 "periodic.loci"
    void calculate(Entity _e_) { 
#line 101 "periodic.loci"
    L_p_ic_[L_cr_[_e_]]=L_p_ic_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 102 "periodic.loci"
      do_loop(seq,this) ;
#line 102 "periodic.loci"
    }
#line 102 "periodic.loci"
} ;
#line 102 "periodic.loci"
Loci::register_rule<file_periodic012_1278520190m668> register_file_periodic012_1278520190m668 ;
#line 102 "periodic.loci"
}
#line 102 "periodic.loci"


  // Pressure.
  // $type p store<real> 
  namespace {class file_periodic013_1278520190m669 : public Loci::pointwise_rule {
#line 106 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 106 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 106 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 106 "periodic.loci"
    Loci::store<real>  L_p_ ; 
#line 106 "periodic.loci"
public:
#line 106 "periodic.loci"
    file_periodic013_1278520190m669() {
#line 106 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 106 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 106 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 106 "periodic.loci"
       name_store("p",L_p_) ;
#line 106 "periodic.loci"
       input("pmap->cl->p") ;
#line 106 "periodic.loci"
       output("cr->p") ;
#line 106 "periodic.loci"
    }
#line 106 "periodic.loci"
    void calculate(Entity _e_) { 
#line 107 "periodic.loci"
    L_p_[L_cr_[_e_]]=L_p_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 108 "periodic.loci"
      do_loop(seq,this) ;
#line 108 "periodic.loci"
    }
#line 108 "periodic.loci"
} ;
#line 108 "periodic.loci"
Loci::register_rule<file_periodic013_1278520190m669> register_file_periodic013_1278520190m669 ;
#line 108 "periodic.loci"
}
#line 108 "periodic.loci"


  // Pressure correction.
  // $type pPrime store<real> 
  namespace {class file_periodic014_1278520190m669 : public Loci::pointwise_rule {
#line 112 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 112 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 112 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 112 "periodic.loci"
    Loci::store<real>  L_pPrime_ ; 
#line 112 "periodic.loci"
public:
#line 112 "periodic.loci"
    file_periodic014_1278520190m669() {
#line 112 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 112 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 112 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 112 "periodic.loci"
       name_store("pPrime",L_pPrime_) ;
#line 112 "periodic.loci"
       input("pmap->cl->pPrime") ;
#line 112 "periodic.loci"
       output("cr->pPrime") ;
#line 112 "periodic.loci"
    }
#line 112 "periodic.loci"
    void calculate(Entity _e_) { 
#line 113 "periodic.loci"
    L_pPrime_[L_cr_[_e_]]=L_pPrime_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 114 "periodic.loci"
      do_loop(seq,this) ;
#line 114 "periodic.loci"
    }
#line 114 "periodic.loci"
} ;
#line 114 "periodic.loci"
Loci::register_rule<file_periodic014_1278520190m669> register_file_periodic014_1278520190m669 ;
#line 114 "periodic.loci"
}
#line 114 "periodic.loci"


  // Temperature initial condition.
  // $type T_ic store<real> 
  namespace {class file_periodic015_1278520190m670 : public Loci::pointwise_rule {
#line 118 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 118 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 118 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 118 "periodic.loci"
    Loci::store<real>  L_T_ic_ ; 
#line 118 "periodic.loci"
public:
#line 118 "periodic.loci"
    file_periodic015_1278520190m670() {
#line 118 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 118 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 118 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 118 "periodic.loci"
       name_store("T_ic",L_T_ic_) ;
#line 118 "periodic.loci"
       input("pmap->cl->T_ic") ;
#line 118 "periodic.loci"
       output("cr->T_ic") ;
#line 118 "periodic.loci"
    }
#line 118 "periodic.loci"
    void calculate(Entity _e_) { 
#line 119 "periodic.loci"
    L_T_ic_[L_cr_[_e_]]=L_T_ic_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 120 "periodic.loci"
      do_loop(seq,this) ;
#line 120 "periodic.loci"
    }
#line 120 "periodic.loci"
} ;
#line 120 "periodic.loci"
Loci::register_rule<file_periodic015_1278520190m670> register_file_periodic015_1278520190m670 ;
#line 120 "periodic.loci"
}
#line 120 "periodic.loci"


  // Temperature.
  // $type temperature store<real> 
  namespace {class file_periodic016_1278520190m670 : public Loci::pointwise_rule {
#line 124 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 124 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 124 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 124 "periodic.loci"
    Loci::store<real>  L_temperature_ ; 
#line 124 "periodic.loci"
public:
#line 124 "periodic.loci"
    file_periodic016_1278520190m670() {
#line 124 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 124 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 124 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 124 "periodic.loci"
       name_store("temperature",L_temperature_) ;
#line 124 "periodic.loci"
       input("pmap->cl->temperature") ;
#line 124 "periodic.loci"
       output("cr->temperature") ;
#line 124 "periodic.loci"
    }
#line 124 "periodic.loci"
    void calculate(Entity _e_) { 
#line 125 "periodic.loci"
    L_temperature_[L_cr_[_e_]]=L_temperature_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 126 "periodic.loci"
      do_loop(seq,this) ;
#line 126 "periodic.loci"
    }
#line 126 "periodic.loci"
} ;
#line 126 "periodic.loci"
Loci::register_rule<file_periodic016_1278520190m670> register_file_periodic016_1278520190m670 ;
#line 126 "periodic.loci"
}
#line 126 "periodic.loci"


  // Total enthalpy.
  // $type h store<real> 
  namespace {class file_periodic017_1278520190m670 : public Loci::pointwise_rule {
#line 130 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 130 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 130 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 130 "periodic.loci"
    Loci::store<real>  L_h_ ; 
#line 130 "periodic.loci"
public:
#line 130 "periodic.loci"
    file_periodic017_1278520190m670() {
#line 130 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 130 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 130 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 130 "periodic.loci"
       name_store("h",L_h_) ;
#line 130 "periodic.loci"
       input("pmap->cl->h") ;
#line 130 "periodic.loci"
       output("cr->h") ;
#line 130 "periodic.loci"
    }
#line 130 "periodic.loci"
    void calculate(Entity _e_) { 
#line 131 "periodic.loci"
    L_h_[L_cr_[_e_]]=L_h_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 132 "periodic.loci"
      do_loop(seq,this) ;
#line 132 "periodic.loci"
    }
#line 132 "periodic.loci"
} ;
#line 132 "periodic.loci"
Loci::register_rule<file_periodic017_1278520190m670> register_file_periodic017_1278520190m670 ;
#line 132 "periodic.loci"
}
#line 132 "periodic.loci"


  // Updated total enthalpy.
  // $type hStar store<real> 
  namespace {class file_periodic018_1278520190m671 : public Loci::pointwise_rule {
#line 136 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 136 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 136 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 136 "periodic.loci"
    Loci::store<real>  L_hStar_ ; 
#line 136 "periodic.loci"
public:
#line 136 "periodic.loci"
    file_periodic018_1278520190m671() {
#line 136 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 136 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 136 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 136 "periodic.loci"
       name_store("hStar",L_hStar_) ;
#line 136 "periodic.loci"
       input("pmap->cl->hStar") ;
#line 136 "periodic.loci"
       output("cr->hStar") ;
#line 136 "periodic.loci"
    }
#line 136 "periodic.loci"
    void calculate(Entity _e_) { 
#line 137 "periodic.loci"
    L_hStar_[L_cr_[_e_]]=L_hStar_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 138 "periodic.loci"
      do_loop(seq,this) ;
#line 138 "periodic.loci"
    }
#line 138 "periodic.loci"
} ;
#line 138 "periodic.loci"
Loci::register_rule<file_periodic018_1278520190m671> register_file_periodic018_1278520190m671 ;
#line 138 "periodic.loci"
}
#line 138 "periodic.loci"


  // K initial condition.
  // $type k_ic store<real> 
  namespace {class file_periodic019_1278520190m671 : public Loci::pointwise_rule {
#line 142 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 142 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 142 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 142 "periodic.loci"
    Loci::store<real>  L_k_ic_ ; 
#line 142 "periodic.loci"
public:
#line 142 "periodic.loci"
    file_periodic019_1278520190m671() {
#line 142 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 142 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 142 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 142 "periodic.loci"
       name_store("k_ic",L_k_ic_) ;
#line 142 "periodic.loci"
       input("pmap->cl->k_ic") ;
#line 142 "periodic.loci"
       output("cr->k_ic") ;
#line 142 "periodic.loci"
    }
#line 142 "periodic.loci"
    void calculate(Entity _e_) { 
#line 143 "periodic.loci"
    L_k_ic_[L_cr_[_e_]]=L_k_ic_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 144 "periodic.loci"
      do_loop(seq,this) ;
#line 144 "periodic.loci"
    }
#line 144 "periodic.loci"
} ;
#line 144 "periodic.loci"
Loci::register_rule<file_periodic019_1278520190m671> register_file_periodic019_1278520190m671 ;
#line 144 "periodic.loci"
}
#line 144 "periodic.loci"


  // K.
  // $type k store<real> 
  namespace {class file_periodic020_1278520190m672 : public Loci::pointwise_rule {
#line 148 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 148 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 148 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 148 "periodic.loci"
    Loci::store<real>  L_k_ ; 
#line 148 "periodic.loci"
public:
#line 148 "periodic.loci"
    file_periodic020_1278520190m672() {
#line 148 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 148 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 148 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 148 "periodic.loci"
       name_store("k",L_k_) ;
#line 148 "periodic.loci"
       input("pmap->cl->k") ;
#line 148 "periodic.loci"
       output("cr->k") ;
#line 148 "periodic.loci"
    }
#line 148 "periodic.loci"
    void calculate(Entity _e_) { 
#line 149 "periodic.loci"
    L_k_[L_cr_[_e_]]=L_k_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 150 "periodic.loci"
      do_loop(seq,this) ;
#line 150 "periodic.loci"
    }
#line 150 "periodic.loci"
} ;
#line 150 "periodic.loci"
Loci::register_rule<file_periodic020_1278520190m672> register_file_periodic020_1278520190m672 ;
#line 150 "periodic.loci"
}
#line 150 "periodic.loci"


  // Updated k.
  // $type kStar store<real> 
  namespace {class file_periodic021_1278520190m672 : public Loci::pointwise_rule {
#line 154 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 154 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 154 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 154 "periodic.loci"
    Loci::store<real>  L_kStar_ ; 
#line 154 "periodic.loci"
public:
#line 154 "periodic.loci"
    file_periodic021_1278520190m672() {
#line 154 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 154 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 154 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 154 "periodic.loci"
       name_store("kStar",L_kStar_) ;
#line 154 "periodic.loci"
       input("pmap->cl->kStar") ;
#line 154 "periodic.loci"
       output("cr->kStar") ;
#line 154 "periodic.loci"
    }
#line 154 "periodic.loci"
    void calculate(Entity _e_) { 
#line 155 "periodic.loci"
    L_kStar_[L_cr_[_e_]]=L_kStar_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 156 "periodic.loci"
      do_loop(seq,this) ;
#line 156 "periodic.loci"
    }
#line 156 "periodic.loci"
} ;
#line 156 "periodic.loci"
Loci::register_rule<file_periodic021_1278520190m672> register_file_periodic021_1278520190m672 ;
#line 156 "periodic.loci"
}
#line 156 "periodic.loci"


  // Omega initial condition.
  // $type omega_ic store<real> 
  namespace {class file_periodic022_1278520190m672 : public Loci::pointwise_rule {
#line 160 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 160 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 160 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 160 "periodic.loci"
    Loci::store<real>  L_omega_ic_ ; 
#line 160 "periodic.loci"
public:
#line 160 "periodic.loci"
    file_periodic022_1278520190m672() {
#line 160 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 160 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 160 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 160 "periodic.loci"
       name_store("omega_ic",L_omega_ic_) ;
#line 160 "periodic.loci"
       input("pmap->cl->omega_ic") ;
#line 160 "periodic.loci"
       output("cr->omega_ic") ;
#line 160 "periodic.loci"
    }
#line 160 "periodic.loci"
    void calculate(Entity _e_) { 
#line 161 "periodic.loci"
    L_omega_ic_[L_cr_[_e_]]=L_omega_ic_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 162 "periodic.loci"
      do_loop(seq,this) ;
#line 162 "periodic.loci"
    }
#line 162 "periodic.loci"
} ;
#line 162 "periodic.loci"
Loci::register_rule<file_periodic022_1278520190m672> register_file_periodic022_1278520190m672 ;
#line 162 "periodic.loci"
}
#line 162 "periodic.loci"


  // Omega.
  // $type omega store<real> 
  namespace {class file_periodic023_1278520190m673 : public Loci::pointwise_rule {
#line 166 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 166 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 166 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 166 "periodic.loci"
    Loci::store<real>  L_omega_ ; 
#line 166 "periodic.loci"
public:
#line 166 "periodic.loci"
    file_periodic023_1278520190m673() {
#line 166 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 166 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 166 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 166 "periodic.loci"
       name_store("omega",L_omega_) ;
#line 166 "periodic.loci"
       input("pmap->cl->omega") ;
#line 166 "periodic.loci"
       output("cr->omega") ;
#line 166 "periodic.loci"
    }
#line 166 "periodic.loci"
    void calculate(Entity _e_) { 
#line 167 "periodic.loci"
    L_omega_[L_cr_[_e_]]=L_omega_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 168 "periodic.loci"
      do_loop(seq,this) ;
#line 168 "periodic.loci"
    }
#line 168 "periodic.loci"
} ;
#line 168 "periodic.loci"
Loci::register_rule<file_periodic023_1278520190m673> register_file_periodic023_1278520190m673 ;
#line 168 "periodic.loci"
}
#line 168 "periodic.loci"


  // Updated omega.
  // $type omegaStar store<real> 
  namespace {class file_periodic024_1278520190m673 : public Loci::pointwise_rule {
#line 172 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 172 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 172 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 172 "periodic.loci"
    Loci::store<real>  L_omegaStar_ ; 
#line 172 "periodic.loci"
public:
#line 172 "periodic.loci"
    file_periodic024_1278520190m673() {
#line 172 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 172 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 172 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 172 "periodic.loci"
       name_store("omegaStar",L_omegaStar_) ;
#line 172 "periodic.loci"
       input("pmap->cl->omegaStar") ;
#line 172 "periodic.loci"
       output("cr->omegaStar") ;
#line 172 "periodic.loci"
    }
#line 172 "periodic.loci"
    void calculate(Entity _e_) { 
#line 173 "periodic.loci"
    L_omegaStar_[L_cr_[_e_]]=L_omegaStar_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 174 "periodic.loci"
      do_loop(seq,this) ;
#line 174 "periodic.loci"
    }
#line 174 "periodic.loci"
} ;
#line 174 "periodic.loci"
Loci::register_rule<file_periodic024_1278520190m673> register_file_periodic024_1278520190m673 ;
#line 174 "periodic.loci"
}
#line 174 "periodic.loci"


  // F1.
  // $type f1 store<real> 
  namespace {class file_periodic025_1278520190m673 : public Loci::pointwise_rule {
#line 178 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 178 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 178 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 178 "periodic.loci"
    Loci::store<real>  L_f1_ ; 
#line 178 "periodic.loci"
public:
#line 178 "periodic.loci"
    file_periodic025_1278520190m673() {
#line 178 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 178 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 178 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 178 "periodic.loci"
       name_store("f1",L_f1_) ;
#line 178 "periodic.loci"
       input("pmap->cl->f1") ;
#line 178 "periodic.loci"
       output("cr->f1") ;
#line 178 "periodic.loci"
    }
#line 178 "periodic.loci"
    void calculate(Entity _e_) { 
#line 179 "periodic.loci"
    L_f1_[L_cr_[_e_]]=L_f1_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 180 "periodic.loci"
      do_loop(seq,this) ;
#line 180 "periodic.loci"
    }
#line 180 "periodic.loci"
} ;
#line 180 "periodic.loci"
Loci::register_rule<file_periodic025_1278520190m673> register_file_periodic025_1278520190m673 ;
#line 180 "periodic.loci"
}
#line 180 "periodic.loci"


  // F2.
  // $type f2 store<real> 
  namespace {class file_periodic026_1278520190m674 : public Loci::pointwise_rule {
#line 184 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 184 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 184 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 184 "periodic.loci"
    Loci::store<real>  L_f2_ ; 
#line 184 "periodic.loci"
public:
#line 184 "periodic.loci"
    file_periodic026_1278520190m674() {
#line 184 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 184 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 184 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 184 "periodic.loci"
       name_store("f2",L_f2_) ;
#line 184 "periodic.loci"
       input("pmap->cl->f2") ;
#line 184 "periodic.loci"
       output("cr->f2") ;
#line 184 "periodic.loci"
    }
#line 184 "periodic.loci"
    void calculate(Entity _e_) { 
#line 185 "periodic.loci"
    L_f2_[L_cr_[_e_]]=L_f2_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 186 "periodic.loci"
      do_loop(seq,this) ;
#line 186 "periodic.loci"
    }
#line 186 "periodic.loci"
} ;
#line 186 "periodic.loci"
Loci::register_rule<file_periodic026_1278520190m674> register_file_periodic026_1278520190m674 ;
#line 186 "periodic.loci"
}
#line 186 "periodic.loci"


  // Species initial condition.
  // $type y_ic storeVec<real> 
  namespace {class file_periodic027_1278520190m674 : public Loci::pointwise_rule {
#line 190 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 190 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 190 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 190 "periodic.loci"
    Loci::storeVec<real>  L_y_ic_ ; 
#line 190 "periodic.loci"
public:
#line 190 "periodic.loci"
    file_periodic027_1278520190m674() {
#line 190 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 190 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 190 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 190 "periodic.loci"
       name_store("y_ic",L_y_ic_) ;
#line 190 "periodic.loci"
       input("pmap->cl->y_ic") ;
#line 190 "periodic.loci"
       output("cr->y_ic") ;
#line 190 "periodic.loci"
    }
#line 190 "periodic.loci"
    void calculate(Entity _e_) { 
#line 191 "periodic.loci"
    L_y_ic_[L_cr_[_e_]]=L_y_ic_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 192 "periodic.loci"
      do_loop(seq,this) ;
#line 192 "periodic.loci"
    }
#line 192 "periodic.loci"
} ;
#line 192 "periodic.loci"
Loci::register_rule<file_periodic027_1278520190m674> register_file_periodic027_1278520190m674 ;
#line 192 "periodic.loci"
}
#line 192 "periodic.loci"


  // Species.
  // $type y storeVec<real> 
  namespace {class file_periodic028_1278520190m674 : public Loci::pointwise_rule {
#line 196 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 196 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 196 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 196 "periodic.loci"
    Loci::storeVec<real>  L_y_ ; 
#line 196 "periodic.loci"
public:
#line 196 "periodic.loci"
    file_periodic028_1278520190m674() {
#line 196 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 196 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 196 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 196 "periodic.loci"
       name_store("y",L_y_) ;
#line 196 "periodic.loci"
       input("pmap->cl->y") ;
#line 196 "periodic.loci"
       output("cr->y") ;
#line 196 "periodic.loci"
    }
#line 196 "periodic.loci"
    void calculate(Entity _e_) { 
#line 197 "periodic.loci"
    L_y_[L_cr_[_e_]]=L_y_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 198 "periodic.loci"
      do_loop(seq,this) ;
#line 198 "periodic.loci"
    }
#line 198 "periodic.loci"
} ;
#line 198 "periodic.loci"
Loci::register_rule<file_periodic028_1278520190m674> register_file_periodic028_1278520190m674 ;
#line 198 "periodic.loci"
}
#line 198 "periodic.loci"


  // Current species.
  // $type yCurr store<real> 
  namespace {class file_periodic029_1278520190m675 : public Loci::pointwise_rule {
#line 202 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 202 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 202 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 202 "periodic.loci"
    Loci::store<real>  L_yCurr_ ; 
#line 202 "periodic.loci"
public:
#line 202 "periodic.loci"
    file_periodic029_1278520190m675() {
#line 202 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 202 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 202 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 202 "periodic.loci"
       name_store("yCurr",L_yCurr_) ;
#line 202 "periodic.loci"
       input("pmap->cl->yCurr") ;
#line 202 "periodic.loci"
       output("cr->yCurr") ;
#line 202 "periodic.loci"
    }
#line 202 "periodic.loci"
    void calculate(Entity _e_) { 
#line 203 "periodic.loci"
    L_yCurr_[L_cr_[_e_]]=L_yCurr_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 204 "periodic.loci"
      do_loop(seq,this) ;
#line 204 "periodic.loci"
    }
#line 204 "periodic.loci"
} ;
#line 204 "periodic.loci"
Loci::register_rule<file_periodic029_1278520190m675> register_file_periodic029_1278520190m675 ;
#line 204 "periodic.loci"
}
#line 204 "periodic.loci"


  // Updated species.
  // $type yCurrStar store<real> 
  namespace {class file_periodic030_1278520190m675 : public Loci::pointwise_rule {
#line 208 "periodic.loci"
    Loci::const_Map L_cl_ ; 
#line 208 "periodic.loci"
    Loci::const_Map L_cr_ ; 
#line 208 "periodic.loci"
    Loci::const_Map L_pmap_ ; 
#line 208 "periodic.loci"
    Loci::store<real>  L_yCurrStar_ ; 
#line 208 "periodic.loci"
public:
#line 208 "periodic.loci"
    file_periodic030_1278520190m675() {
#line 208 "periodic.loci"
       name_store("cl",L_cl_) ;
#line 208 "periodic.loci"
       name_store("cr",L_cr_) ;
#line 208 "periodic.loci"
       name_store("pmap",L_pmap_) ;
#line 208 "periodic.loci"
       name_store("yCurrStar",L_yCurrStar_) ;
#line 208 "periodic.loci"
       input("pmap->cl->yCurrStar") ;
#line 208 "periodic.loci"
       output("cr->yCurrStar") ;
#line 208 "periodic.loci"
    }
#line 208 "periodic.loci"
    void calculate(Entity _e_) { 
#line 209 "periodic.loci"
    L_yCurrStar_[L_cr_[_e_]]=L_yCurrStar_[L_cl_[L_pmap_[_e_]]];
  }    void compute(const Loci::sequence &seq) { 
#line 210 "periodic.loci"
      do_loop(seq,this) ;
#line 210 "periodic.loci"
    }
#line 210 "periodic.loci"
} ;
#line 210 "periodic.loci"
Loci::register_rule<file_periodic030_1278520190m675> register_file_periodic030_1278520190m675 ;
#line 210 "periodic.loci"
}
#line 210 "periodic.loci"

}
