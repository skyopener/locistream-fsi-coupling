#line 1 "eos_interface.loci"
#include <Loci.h>
#include "sciTypes.h"
#include "eos.h"
#include "PerfectGas.h"
#include "reaction.h"
#include "name_var.h"
#include "qvi.h"
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
#line 8 "eos_interface.loci"


namespace streamUns {
  using namespace fluidPhysics ;

  // $type chemistry_model param<name_var> 

  namespace {class file_eos_interface000_1280810721m242 : public Loci::default_rule {
#line 17 "eos_interface.loci"
    Loci::param<name_var>  L_chemistry_model_ ; 
#line 17 "eos_interface.loci"
public:
#line 17 "eos_interface.loci"
    file_eos_interface000_1280810721m242() {
#line 17 "eos_interface.loci"
       name_store("chemistry_model",L_chemistry_model_) ;
#line 17 "eos_interface.loci"
       output("chemistry_model") ;
#line 17 "eos_interface.loci"
       comments("Select the thermodynamic and chemistry model used") ;
#line 17 "eos_interface.loci"
    }
#line 17 "eos_interface.loci"
    void compute(const Loci::sequence &seq) { 
    (*L_chemistry_model_).name = "air_1s0r" ;
  }} ;
#line 19 "eos_interface.loci"
Loci::register_rule<file_eos_interface000_1280810721m242> register_file_eos_interface000_1280810721m242 ;
#line 19 "eos_interface.loci"
}
#line 19 "eos_interface.loci"
 

  // $type eos_type param<name_var> 

  namespace {class file_eos_interface001_1280810721m242 : public Loci::default_rule {
#line 24 "eos_interface.loci"
    Loci::param<name_var>  L_eos_type_ ; 
#line 24 "eos_interface.loci"
public:
#line 24 "eos_interface.loci"
    file_eos_interface001_1280810721m242() {
#line 24 "eos_interface.loci"
       name_store("eos_type",L_eos_type_) ;
#line 24 "eos_interface.loci"
       output("eos_type") ;
#line 24 "eos_interface.loci"
       comments("Select the EoS to use.  Thermally perfect gases use 'gas', other EoS's such as fluid may be available when modules are loaded.  For example, 'fluid' for real fluids.") ;
#line 24 "eos_interface.loci"
    }
#line 24 "eos_interface.loci"
    void compute(const Loci::sequence &seq) { 
    (*L_eos_type_).name = "gas" ;
  }} ;
#line 26 "eos_interface.loci"
Loci::register_rule<file_eos_interface001_1280810721m242> register_file_eos_interface001_1280810721m242 ;
#line 26 "eos_interface.loci"
}
#line 26 "eos_interface.loci"

    
  // $type thermodynamic_model param<name_var> 

  namespace {class file_eos_interface002_1280810721m243 : public Loci::default_rule {
#line 32 "eos_interface.loci"
    Loci::param<name_var>  L_thermodynamic_model_ ; 
#line 32 "eos_interface.loci"
public:
#line 32 "eos_interface.loci"
    file_eos_interface002_1280810721m243() {
#line 32 "eos_interface.loci"
       name_store("thermodynamic_model",L_thermodynamic_model_) ;
#line 32 "eos_interface.loci"
       output("thermodynamic_model") ;
#line 32 "eos_interface.loci"
       comments("Define type of thermodynamic model used.  Currently the choices are curve fit and vibrational equilibrium.  'adaptive' allows the code to select the model based on what is specified in the thermodynamic database.") ;
#line 32 "eos_interface.loci"
    }
#line 32 "eos_interface.loci"
    void compute(const Loci::sequence &seq) { 
    (*L_thermodynamic_model_).name = "adaptive" ;
  }} ;
#line 34 "eos_interface.loci"
Loci::register_rule<file_eos_interface002_1280810721m243> register_file_eos_interface002_1280810721m243 ;
#line 34 "eos_interface.loci"
}
#line 34 "eos_interface.loci"
 


  chemistry_db read_chemistry_db(string filename,const chemistry_db &defaults) {
    string filedata ;
    if(Loci::MPI_rank == 0) {
      ifstream cdf(filename.c_str(),ios::in) ;
      if(!cdf.fail()) {
        char c ;
        while(cdf.get(c)) 
          filedata+= c ;
        cdf.close() ;
      }
    }
    int sz = filedata.size() ;
    MPI_Bcast(&sz,1,MPI_INT,0,MPI_COMM_WORLD) ;
    char *buf = new char[sz] ;
    if(Loci::MPI_rank==0) {
      for(int i=0;i<sz;++i)
        buf[i] = filedata[i] ;
    }
    MPI_Bcast(buf,sz,MPI_CHAR,0,MPI_COMM_WORLD) ;
    if(Loci::MPI_rank!=0) {
      filedata = "" ;
      for(int i=0;i<sz;++i)
        filedata+= buf[i] ;
    }

    istringstream iss(filedata) ;

    chemistry_db cdb ;
    cdb.Input(iss,defaults) ;
    return cdb ;
  }    

  chemistry_db read_chemistry_model(string mdl) {
    string fname = mdl + ".mdl" ;
    string dbase_file ;
    const char *base = 0 ;

    if(Loci::MPI_rank == 0) {
      base = getenv("CHEMISTRY_DATABASE") ;
      if(base == NULL) {
        cerr << "Warning: CHEMISTRY_DATABASE environment variable not defined"
             << endl ;
        base = "./" ;
      }

      // Set up default database path
      char dbname[512] ;
      sprintf(dbname,"%s/data_base/chemistry",base) ;
      dbase_file = string(dbname) ;

      ifstream chemin(fname.c_str(),ios::in) ;
      if(chemin.fail()) {
        fname = string(base)+"/data_base/models/" + mdl + ".mdl";
        ifstream chemin(fname.c_str(),ios::in) ;
        if(chemin.fail()) {
          cerr << "can't find model '"
               << mdl <<"' looked for file '"
               << fname << "' with no luck." << endl ;
          Loci::Abort() ;
        } 
      } else {
        cout << "Using model from local directory, filename = "
             << fname << endl ;
      }
    }

    chemistry_db dummy ; // empty chemistry database for default defaults
    chemistry_db cdb_defaults = read_chemistry_db(dbase_file,dummy) ;
    
    chemistry_db cdb = read_chemistry_db(fname,cdb_defaults) ;

    cdb.removeM() ; // process with m-body reactions
    return cdb ;
  }


  // $type idealGasSimulation Constraint
  // $type imperfectGasSimulation Constraint
  // $type singleSpeciesSimulation Constraint
  // $type multiComponentSimulation Constraint
  // $type reactionMechanismSimulation Constraint

  namespace {class file_eos_interface003_1280810721m245 : public Loci::constraint_rule {
#line 122 "eos_interface.loci"
    Loci::const_param<name_var>  L_chemistry_model_ ; 
#line 122 "eos_interface.loci"
    Loci::Constraint L_idealGasSimulation_ ; 
#line 122 "eos_interface.loci"
    Loci::Constraint L_imperfectGasSimulation_ ; 
#line 122 "eos_interface.loci"
    Loci::Constraint L_singleSpeciesSimulation_ ; 
#line 122 "eos_interface.loci"
    Loci::Constraint L_multiComponentSimulation_ ; 
#line 122 "eos_interface.loci"
    Loci::Constraint L_reactionMechanismSimulation_ ; 
#line 122 "eos_interface.loci"
public:
#line 122 "eos_interface.loci"
    file_eos_interface003_1280810721m245() {
#line 122 "eos_interface.loci"
       name_store("chemistry_model",L_chemistry_model_) ;
#line 122 "eos_interface.loci"
       name_store("idealGasSimulation",L_idealGasSimulation_) ;
#line 122 "eos_interface.loci"
       name_store("imperfectGasSimulation",L_imperfectGasSimulation_) ;
#line 122 "eos_interface.loci"
       name_store("singleSpeciesSimulation",L_singleSpeciesSimulation_) ;
#line 122 "eos_interface.loci"
       name_store("multiComponentSimulation",L_multiComponentSimulation_) ;
#line 122 "eos_interface.loci"
       name_store("reactionMechanismSimulation",L_reactionMechanismSimulation_) ;
#line 122 "eos_interface.loci"
       input("chemistry_model") ;
#line 122 "eos_interface.loci"
       output("idealGasSimulation") ;
#line 122 "eos_interface.loci"
       output("imperfectGasSimulation") ;
#line 122 "eos_interface.loci"
       output("singleSpeciesSimulation") ;
#line 122 "eos_interface.loci"
       output("multiComponentSimulation") ;
#line 122 "eos_interface.loci"
       output("reactionMechanismSimulation") ;
#line 122 "eos_interface.loci"
    }
#line 122 "eos_interface.loci"
    void compute(const Loci::sequence &seq) { 
    L_idealGasSimulation_= EMPTY ;
    L_imperfectGasSimulation_= ~EMPTY ;
    L_singleSpeciesSimulation_= EMPTY ;
    L_multiComponentSimulation_= ~EMPTY ;
    L_reactionMechanismSimulation_= ~EMPTY ;
    string mdl = (*L_chemistry_model_).name ;
    chemistry_db cdb = read_chemistry_model(mdl) ;
    if(cdb.species.numSpecies() == 1) {
      L_singleSpeciesSimulation_= ~EMPTY ;
      L_multiComponentSimulation_= EMPTY ;
      L_reactionMechanismSimulation_= EMPTY ;
      const options_list &ol = cdb.species.getSpeciesOption(0) ;
      if(!ol.optionExists("cp") && !ol.optionExists("theta_v") &&
         !ol.optionExists("eos")) {
        L_idealGasSimulation_= ~EMPTY ;
        L_imperfectGasSimulation_= EMPTY ;
      }
      
    }
    if(cdb.reactions.numReactions() == 0)
      L_reactionMechanismSimulation_= EMPTY ;
  }} ;
#line 144 "eos_interface.loci"
Loci::register_rule<file_eos_interface003_1280810721m245> register_file_eos_interface003_1280810721m245 ;
#line 144 "eos_interface.loci"
}
#line 144 "eos_interface.loci"


  // $type eos_repository blackbox<EOSFactory> 

  namespace {class file_eos_interface004_1280810721m246 : public Loci::unit_rule {
#line 148 "eos_interface.loci"
    Loci::blackbox<EOSFactory>  L_eos_repository_ ; 
#line 148 "eos_interface.loci"
public:
#line 148 "eos_interface.loci"
    file_eos_interface004_1280810721m246() {
#line 148 "eos_interface.loci"
       name_store("eos_repository",L_eos_repository_) ;
#line 148 "eos_interface.loci"
       output("eos_repository") ;
#line 148 "eos_interface.loci"
       constraint("UNIVERSE") ;
#line 148 "eos_interface.loci"
    }
#line 148 "eos_interface.loci"
    void calculate(Entity _e_) { 
#line 149 "eos_interface.loci"
  }    void compute(const Loci::sequence &seq) { 
#line 149 "eos_interface.loci"
      do_loop(seq,this) ;
#line 149 "eos_interface.loci"
    }
#line 149 "eos_interface.loci"
} ;
#line 149 "eos_interface.loci"
Loci::register_rule<file_eos_interface004_1280810721m246> register_file_eos_interface004_1280810721m246 ;
#line 149 "eos_interface.loci"
}
#line 149 "eos_interface.loci"


  // $type thermallyPerfectPriority param<int> 
  namespace {class file_eos_interface005_1280810721m246 : public Loci::singleton_rule {
#line 152 "eos_interface.loci"
    Loci::param<int>  L_thermallyPerfectPriority_ ; 
#line 152 "eos_interface.loci"
public:
#line 152 "eos_interface.loci"
    file_eos_interface005_1280810721m246() {
#line 152 "eos_interface.loci"
       name_store("thermallyPerfectPriority",L_thermallyPerfectPriority_) ;
#line 152 "eos_interface.loci"
       output("thermallyPerfectPriority") ;
#line 152 "eos_interface.loci"
       constraint("UNIVERSE") ;
#line 152 "eos_interface.loci"
    }
#line 152 "eos_interface.loci"
    void compute(const Loci::sequence &seq) { 
    (*L_thermallyPerfectPriority_)= 0 ;
  }} ;
#line 154 "eos_interface.loci"
Loci::register_rule<file_eos_interface005_1280810721m246> register_file_eos_interface005_1280810721m246 ;
#line 154 "eos_interface.loci"
}
#line 154 "eos_interface.loci"

  
  namespace {class file_eos_interface006_1280810721m247 : public Loci::apply_rule< blackbox<EOSFactory> ,Loci::NullOp<EOSFactory> >  {
#line 157 "eos_interface.loci"
    Loci::const_param<int>  L_thermallyPerfectPriority_ ; 
#line 157 "eos_interface.loci"
    Loci::blackbox<EOSFactory>  L_eos_repository_ ; 
#line 157 "eos_interface.loci"
public:
#line 157 "eos_interface.loci"
    file_eos_interface006_1280810721m247() {
#line 157 "eos_interface.loci"
       name_store("eos_repository",L_eos_repository_) ;
#line 157 "eos_interface.loci"
       name_store("thermallyPerfectPriority",L_thermallyPerfectPriority_) ;
#line 157 "eos_interface.loci"
       input("thermallyPerfectPriority") ;
#line 157 "eos_interface.loci"
       output("eos_repository") ;
#line 157 "eos_interface.loci"
    }
#line 157 "eos_interface.loci"
    void prelude(const Loci::sequence &seq) { 
    L_eos_repository_->insertEOS(EoSPtr(new ThermallyPerfectEOS),
                               "gas",*L_thermallyPerfectPriority_) ;
  }    void compute(const Loci::sequence &seq) { 
#line 162 "eos_interface.loci"
      prelude(seq) ;
#line 162 "eos_interface.loci"
    }
#line 162 "eos_interface.loci"
} ;
#line 162 "eos_interface.loci"
Loci::register_rule<file_eos_interface006_1280810721m247> register_file_eos_interface006_1280810721m247 ;
#line 162 "eos_interface.loci"
}
#line 162 "eos_interface.loci"
class read_chem : public singleton_rule {
    const_param<name_var> cmodel ;
    const_param<std::string> tmodel ;
    const_param<name_var> thermo_model ;
    const_param<name_var> eos_type ;
    const_blackbox<EOSFactory> eos_repository ;
    param<EOS> eos ;
    param<reaction> reactor ;
  public:
    read_chem() ;
    virtual void compute( const sequence &seq) ;
  } ;

  read_chem::read_chem() {
    name_store("chemistry_model",cmodel) ;
    name_store("turbulence_model",tmodel) ;
    name_store("thermodynamic_model",thermo_model) ;
    name_store("eos_type",eos_type) ;
    name_store("eos_repository",eos_repository) ;
    name_store("eos",eos) ;
    name_store("reactor",reactor) ;
    input("chemistry_model") ;
    input("turbulence_model") ;
    input("thermodynamic_model") ;
    input("eos_type") ;
    input("eos_repository") ;
    output("eos,reactor") ;
  }

  void read_chem::compute(const sequence &seq) {
    string mdl = (*cmodel).name ;
    chemistry_db cdb = read_chemistry_model(mdl) ;
    if(Loci::MPI_rank==0)
      cdb.Print(cout) ;

    // Old version.
    //EoSPtr e = eos_repository->getEOS(eos_type->name) ;
    // Old version.

    // New version.
    EoSPtr e = 0 ;
    if(eos_type->name == "default")
      e = eos_repository->getDefaultEOS() ;
    else
      e = eos_repository->getEOS(eos_type->name) ;
    // New version.

    if(e == 0) {
      if(Loci::MPI_rank == 0) {
        cerr << "Unable to find EoS named " << eos_type->name
             << ", check eos_type setting in vars file." << endl ;
      }
      Loci::Abort(); 
    }
    *eos = EOS(e) ;
    eos->initialize(cdb.species,thermo_model->name) ;

    reactor->initialize(*eos,cdb.reactions,cdb.species) ;

    if(Loci::MPI_rank==0)
      cout << "Using chemistry model = " << mdl << endl ;

  }

  register_rule<read_chem>   register_read_chem ;

  //The following rule gets state variable (eos_state) and species energy
  // $type eos_state store<EOS::State> 
  // $type eos_mixture_state storeVec<real> 
  // $type hint storeVec<float> 
  // $type hint_n storeVec<float> 
  namespace {class file_eos_interface007_1280810721m249 : public Loci::pointwise_rule {
#line 235 "eos_interface.loci"
    Loci::const_param<fluidPhysics::EOS>  L_eos_nit__ ; 
#line 235 "eos_interface.loci"
    Loci::const_param<conservativeVectorInfo>  L_qvi_nit__ ; 
#line 235 "eos_interface.loci"
    Loci::const_storeVec<streamUns::real>  L_qp_nit__ ; 
#line 235 "eos_interface.loci"
    Loci::const_param<streamUns::real>  L_Pambient_nit__ ; 
#line 235 "eos_interface.loci"
    Loci::storeVec<float>  L_hint_n_nit__ ; 
#line 235 "eos_interface.loci"
    Loci::store<EOS::State>  L_eos_state_nit__ ; 
#line 235 "eos_interface.loci"
    Loci::storeVec<real>  L_eos_mixture_state_nit__ ; 
#line 235 "eos_interface.loci"
public:
#line 235 "eos_interface.loci"
    file_eos_interface007_1280810721m249() {
#line 235 "eos_interface.loci"
       name_store("hint_n{n,it}",L_hint_n_nit__) ;
#line 235 "eos_interface.loci"
       name_store("eos{n,it}",L_eos_nit__) ;
#line 235 "eos_interface.loci"
       name_store("qvi{n,it}",L_qvi_nit__) ;
#line 235 "eos_interface.loci"
       name_store("qp{n,it}",L_qp_nit__) ;
#line 235 "eos_interface.loci"
       name_store("Pambient{n,it}",L_Pambient_nit__) ;
#line 235 "eos_interface.loci"
       name_store("eos_state{n,it}",L_eos_state_nit__) ;
#line 235 "eos_interface.loci"
       name_store("eos_mixture_state{n,it}",L_eos_mixture_state_nit__) ;
#line 235 "eos_interface.loci"
       input("eos{n,it},    qvi{n,it},qp{n,it},hint_n{n,it},Pambient{n,it}") ;
#line 235 "eos_interface.loci"
       output("hint{n,it}=hint_n{n,it}") ;
#line 235 "eos_interface.loci"
       output("eos_state{n,it}") ;
#line 235 "eos_interface.loci"
       output("eos_mixture_state{n,it}") ;
#line 235 "eos_interface.loci"
       constraint("qp{n,it}") ;
#line 235 "eos_interface.loci"
    }
#line 235 "eos_interface.loci"
    void prelude(const Loci::sequence &seq) { 
    L_eos_mixture_state_nit__.setVecSize(L_eos_nit__->mixtureStateSize()) ;
  }    void calculate(Entity _e_) { 
#line 238 "eos_interface.loci"
    //state variables are assigned to eos_state, and species_energy are 
    //assigned with the values since it is passed by the pointer
    L_eos_state_nit__[_e_]= primitive_vector ::getState (L_qp_nit__[_e_], L_Pambient_nit__[_e_],
                                                  L_qvi_nit__[_e_],L_eos_nit__[_e_],
                                                  L_eos_mixture_state_nit__[_e_],
                                                  L_hint_n_nit__[_e_]) ;
    L_eos_nit__[_e_].getHint (L_hint_n_nit__[_e_],L_eos_state_nit__[_e_],L_eos_mixture_state_nit__[_e_]) ;
  }    void compute(const Loci::sequence &seq) { 
#line 245 "eos_interface.loci"
      prelude(seq) ;
#line 245 "eos_interface.loci"
      do_loop(seq,this) ;
#line 245 "eos_interface.loci"
    }
#line 245 "eos_interface.loci"
} ;
#line 245 "eos_interface.loci"
Loci::register_rule<file_eos_interface007_1280810721m249> register_file_eos_interface007_1280810721m249 ;
#line 245 "eos_interface.loci"
}
#line 245 "eos_interface.loci"
   

  // $type species_energy storeVec<real> 
  namespace {class file_eos_interface008_1280810721m250 : public Loci::pointwise_rule {
#line 249 "eos_interface.loci"
    Loci::const_param<fluidPhysics::EOS>  L_eos_ ; 
#line 249 "eos_interface.loci"
    Loci::const_store<EOS::State>  L_eos_state_ ; 
#line 249 "eos_interface.loci"
    Loci::const_storeVec<real>  L_eos_mixture_state_ ; 
#line 249 "eos_interface.loci"
    Loci::storeVec<real>  L_species_energy_ ; 
#line 249 "eos_interface.loci"
public:
#line 249 "eos_interface.loci"
    file_eos_interface008_1280810721m250() {
#line 249 "eos_interface.loci"
       name_store("eos",L_eos_) ;
#line 249 "eos_interface.loci"
       name_store("eos_state",L_eos_state_) ;
#line 249 "eos_interface.loci"
       name_store("eos_mixture_state",L_eos_mixture_state_) ;
#line 249 "eos_interface.loci"
       name_store("species_energy",L_species_energy_) ;
#line 249 "eos_interface.loci"
       input("eos,eos_state,eos_mixture_state") ;
#line 249 "eos_interface.loci"
       output("species_energy") ;
#line 249 "eos_interface.loci"
    }
#line 249 "eos_interface.loci"
    void prelude(const Loci::sequence &seq) { 
    L_species_energy_.setVecSize(L_eos_->numSpecies()) ;
  }    void calculate(Entity _e_) { 
#line 252 "eos_interface.loci"
    L_eos_[_e_].get_ei (L_eos_state_[_e_],L_eos_mixture_state_[_e_],L_species_energy_[_e_]) ;
  }    void compute(const Loci::sequence &seq) { 
#line 253 "eos_interface.loci"
      prelude(seq) ;
#line 253 "eos_interface.loci"
      do_loop(seq,this) ;
#line 253 "eos_interface.loci"
    }
#line 253 "eos_interface.loci"
} ;
#line 253 "eos_interface.loci"
Loci::register_rule<file_eos_interface008_1280810721m250> register_file_eos_interface008_1280810721m250 ;
#line 253 "eos_interface.loci"
}
#line 253 "eos_interface.loci"


  // $type species_cv storeVec<real> 
  namespace {class file_eos_interface009_1280810721m251 : public Loci::pointwise_rule {
#line 257 "eos_interface.loci"
    Loci::const_param<fluidPhysics::EOS>  L_eos_ ; 
#line 257 "eos_interface.loci"
    Loci::const_store<EOS::State>  L_eos_state_ ; 
#line 257 "eos_interface.loci"
    Loci::const_storeVec<real>  L_eos_mixture_state_ ; 
#line 257 "eos_interface.loci"
    Loci::storeVec<real>  L_species_cv_ ; 
#line 257 "eos_interface.loci"
public:
#line 257 "eos_interface.loci"
    file_eos_interface009_1280810721m251() {
#line 257 "eos_interface.loci"
       name_store("eos",L_eos_) ;
#line 257 "eos_interface.loci"
       name_store("eos_state",L_eos_state_) ;
#line 257 "eos_interface.loci"
       name_store("eos_mixture_state",L_eos_mixture_state_) ;
#line 257 "eos_interface.loci"
       name_store("species_cv",L_species_cv_) ;
#line 257 "eos_interface.loci"
       input("eos,eos_state,eos_mixture_state") ;
#line 257 "eos_interface.loci"
       output("species_cv") ;
#line 257 "eos_interface.loci"
    }
#line 257 "eos_interface.loci"
    void prelude(const Loci::sequence &seq) { 
    L_species_cv_.setVecSize(L_eos_->numSpecies()) ;
  }    void calculate(Entity _e_) { 
#line 260 "eos_interface.loci"
    L_eos_[_e_].get_cvi (L_eos_state_[_e_],L_eos_mixture_state_[_e_],L_species_cv_[_e_]) ;
  }    void compute(const Loci::sequence &seq) { 
#line 261 "eos_interface.loci"
      prelude(seq) ;
#line 261 "eos_interface.loci"
      do_loop(seq,this) ;
#line 261 "eos_interface.loci"
    }
#line 261 "eos_interface.loci"
} ;
#line 261 "eos_interface.loci"
Loci::register_rule<file_eos_interface009_1280810721m251> register_file_eos_interface009_1280810721m251 ;
#line 261 "eos_interface.loci"
}
#line 261 "eos_interface.loci"


  // $type dreidri storeVec<real> 
  namespace {class file_eos_interface010_1280810721m251 : public Loci::pointwise_rule {
#line 265 "eos_interface.loci"
    Loci::const_param<fluidPhysics::EOS>  L_eos_ ; 
#line 265 "eos_interface.loci"
    Loci::const_store<EOS::State>  L_eos_state_ ; 
#line 265 "eos_interface.loci"
    Loci::const_storeVec<real>  L_eos_mixture_state_ ; 
#line 265 "eos_interface.loci"
    Loci::storeVec<real>  L_dreidri_ ; 
#line 265 "eos_interface.loci"
public:
#line 265 "eos_interface.loci"
    file_eos_interface010_1280810721m251() {
#line 265 "eos_interface.loci"
       name_store("eos",L_eos_) ;
#line 265 "eos_interface.loci"
       name_store("eos_state",L_eos_state_) ;
#line 265 "eos_interface.loci"
       name_store("eos_mixture_state",L_eos_mixture_state_) ;
#line 265 "eos_interface.loci"
       name_store("dreidri",L_dreidri_) ;
#line 265 "eos_interface.loci"
       input("eos,eos_state,eos_mixture_state") ;
#line 265 "eos_interface.loci"
       output("dreidri") ;
#line 265 "eos_interface.loci"
    }
#line 265 "eos_interface.loci"
    void prelude(const Loci::sequence &seq) { 
    L_dreidri_.setVecSize(L_eos_->numSpecies()) ;
  }    void calculate(Entity _e_) { 
#line 268 "eos_interface.loci"
    L_eos_[_e_].dreidri (L_dreidri_[_e_],L_eos_state_[_e_],L_eos_mixture_state_[_e_]) ;
  }    void compute(const Loci::sequence &seq) { 
#line 269 "eos_interface.loci"
      prelude(seq) ;
#line 269 "eos_interface.loci"
      do_loop(seq,this) ;
#line 269 "eos_interface.loci"
    }
#line 269 "eos_interface.loci"
} ;
#line 269 "eos_interface.loci"
Loci::register_rule<file_eos_interface010_1280810721m251> register_file_eos_interface010_1280810721m251 ;
#line 269 "eos_interface.loci"
}
#line 269 "eos_interface.loci"


  // $type dreidP store<real> 
  namespace {class file_eos_interface011_1280810721m252 : public Loci::pointwise_rule {
#line 273 "eos_interface.loci"
    Loci::const_param<fluidPhysics::EOS>  L_eos_ ; 
#line 273 "eos_interface.loci"
    Loci::const_store<EOS::State>  L_eos_state_ ; 
#line 273 "eos_interface.loci"
    Loci::const_storeVec<real>  L_eos_mixture_state_ ; 
#line 273 "eos_interface.loci"
    Loci::store<real>  L_dreidP_ ; 
#line 273 "eos_interface.loci"
public:
#line 273 "eos_interface.loci"
    file_eos_interface011_1280810721m252() {
#line 273 "eos_interface.loci"
       name_store("eos",L_eos_) ;
#line 273 "eos_interface.loci"
       name_store("eos_state",L_eos_state_) ;
#line 273 "eos_interface.loci"
       name_store("eos_mixture_state",L_eos_mixture_state_) ;
#line 273 "eos_interface.loci"
       name_store("dreidP",L_dreidP_) ;
#line 273 "eos_interface.loci"
       input("eos,eos_state,eos_mixture_state") ;
#line 273 "eos_interface.loci"
       output("dreidP") ;
#line 273 "eos_interface.loci"
    }
#line 273 "eos_interface.loci"
    void calculate(Entity _e_) { 
#line 274 "eos_interface.loci"
    L_dreidP_[_e_]= L_eos_[_e_].dreidP (L_eos_state_[_e_],L_eos_mixture_state_[_e_]) ;
  }    void compute(const Loci::sequence &seq) { 
#line 275 "eos_interface.loci"
      do_loop(seq,this) ;
#line 275 "eos_interface.loci"
    }
#line 275 "eos_interface.loci"
} ;
#line 275 "eos_interface.loci"
Loci::register_rule<file_eos_interface011_1280810721m252> register_file_eos_interface011_1280810721m252 ;
#line 275 "eos_interface.loci"
}
#line 275 "eos_interface.loci"


  // These rules conflict with Loci-Stream, so comment out.

  //get temperature in the cell
  //$rule pointwise(temperature<-eos_state) {
  //  $temperature = $eos_state.temperature() ;
  //}


  //get density in the cell
  //$rule pointwise(rho<-eos_state) {
  //  $rho = $eos_state.density() ;
  //}

  //get pressure in the cell

  //$rule pointwise(pressure<-eos_state) {
  //  $pressure = $eos_state.pressure() ;
  //}

  //get soundSpeed in the cell
  // $type soundSpeed store<real> 
  namespace {class file_eos_interface012_1280810721m252 : public Loci::pointwise_rule {
#line 298 "eos_interface.loci"
    Loci::const_store<EOS::State>  L_eos_state_ ; 
#line 298 "eos_interface.loci"
    Loci::store<real>  L_soundSpeed_ ; 
#line 298 "eos_interface.loci"
public:
#line 298 "eos_interface.loci"
    file_eos_interface012_1280810721m252() {
#line 298 "eos_interface.loci"
       name_store("eos_state",L_eos_state_) ;
#line 298 "eos_interface.loci"
       name_store("soundSpeed",L_soundSpeed_) ;
#line 298 "eos_interface.loci"
       input("eos_state") ;
#line 298 "eos_interface.loci"
       output("soundSpeed") ;
#line 298 "eos_interface.loci"
    }
#line 298 "eos_interface.loci"
    void calculate(Entity _e_) { 
#line 299 "eos_interface.loci"
    L_soundSpeed_[_e_]= L_eos_state_[_e_].soundSpeed () ;
  }    void compute(const Loci::sequence &seq) { 
#line 300 "eos_interface.loci"
      do_loop(seq,this) ;
#line 300 "eos_interface.loci"
    }
#line 300 "eos_interface.loci"
} ;
#line 300 "eos_interface.loci"
Loci::register_rule<file_eos_interface012_1280810721m252> register_file_eos_interface012_1280810721m252 ;
#line 300 "eos_interface.loci"
}
#line 300 "eos_interface.loci"


  // $type cp store<real> 
  namespace {class file_eos_interface013_1280810721m253 : public Loci::pointwise_rule {
#line 303 "eos_interface.loci"
    Loci::const_store<EOS::State>  L_eos_state_ ; 
#line 303 "eos_interface.loci"
    Loci::store<real>  L_cp_ ; 
#line 303 "eos_interface.loci"
public:
#line 303 "eos_interface.loci"
    file_eos_interface013_1280810721m253() {
#line 303 "eos_interface.loci"
       name_store("eos_state",L_eos_state_) ;
#line 303 "eos_interface.loci"
       name_store("cp",L_cp_) ;
#line 303 "eos_interface.loci"
       input("eos_state") ;
#line 303 "eos_interface.loci"
       output("cp") ;
#line 303 "eos_interface.loci"
    }
#line 303 "eos_interface.loci"
    void calculate(Entity _e_) { 
#line 304 "eos_interface.loci"
    L_cp_[_e_]= L_eos_state_[_e_].cpt () ;
  }    void compute(const Loci::sequence &seq) { 
#line 305 "eos_interface.loci"
      do_loop(seq,this) ;
#line 305 "eos_interface.loci"
    }
#line 305 "eos_interface.loci"
} ;
#line 305 "eos_interface.loci"
Loci::register_rule<file_eos_interface013_1280810721m253> register_file_eos_interface013_1280810721m253 ;
#line 305 "eos_interface.loci"
}
#line 305 "eos_interface.loci"


  // $type vol store<real> 
  namespace {class file_eos_interface014_1280810721m253 : public Loci::pointwise_rule {
#line 308 "eos_interface.loci"
    Loci::const_Map L_cl_ ; 
#line 308 "eos_interface.loci"
    Loci::const_Map L_cr_ ; 
#line 308 "eos_interface.loci"
    Loci::const_store<real>  L_vol_ ; 
#line 308 "eos_interface.loci"
    Loci::const_store<EOS::State>  L_eos_state_ ; 
#line 308 "eos_interface.loci"
    Loci::store<real>  L_cp_ ; 
#line 308 "eos_interface.loci"
public:
#line 308 "eos_interface.loci"
    file_eos_interface014_1280810721m253() {
#line 308 "eos_interface.loci"
       name_store("cl",L_cl_) ;
#line 308 "eos_interface.loci"
       name_store("cr",L_cr_) ;
#line 308 "eos_interface.loci"
       name_store("vol",L_vol_) ;
#line 308 "eos_interface.loci"
       name_store("eos_state",L_eos_state_) ;
#line 308 "eos_interface.loci"
       name_store("cp",L_cp_) ;
#line 308 "eos_interface.loci"
       input("(cr,cl)->(eos_state,vol)") ;
#line 308 "eos_interface.loci"
       output("cp") ;
#line 308 "eos_interface.loci"
    }
#line 308 "eos_interface.loci"
    void calculate(Entity _e_) { 
#line 309 "eos_interface.loci"
    const real cptr =L_eos_state_[L_cr_[_e_]].cpt () ;
    const real cptl =L_eos_state_[L_cl_[_e_]].cpt () ;
    L_cp_[_e_]=(L_vol_[L_cl_[_e_]]*cptr +L_vol_[L_cr_[_e_]]*cptl )/
      (L_vol_[L_cl_[_e_]]+L_vol_[L_cr_[_e_]]) ;
  }    void compute(const Loci::sequence &seq) { 
#line 313 "eos_interface.loci"
      do_loop(seq,this) ;
#line 313 "eos_interface.loci"
    }
#line 313 "eos_interface.loci"
} ;
#line 313 "eos_interface.loci"
Loci::register_rule<file_eos_interface014_1280810721m253> register_file_eos_interface014_1280810721m253 ;
#line 313 "eos_interface.loci"
}
#line 313 "eos_interface.loci"


  namespace {class file_eos_interface015_1280810721m254 : public Loci::pointwise_rule {
#line 315 "eos_interface.loci"
    Loci::const_Map L_ci_ ; 
#line 315 "eos_interface.loci"
    Loci::const_store<EOS::State>  L_eos_state_ ; 
#line 315 "eos_interface.loci"
    Loci::store<real>  L_cp_ ; 
#line 315 "eos_interface.loci"
public:
#line 315 "eos_interface.loci"
    file_eos_interface015_1280810721m254() {
#line 315 "eos_interface.loci"
       name_store("ci",L_ci_) ;
#line 315 "eos_interface.loci"
       name_store("eos_state",L_eos_state_) ;
#line 315 "eos_interface.loci"
       name_store("cp",L_cp_) ;
#line 315 "eos_interface.loci"
       input("ci->eos_state") ;
#line 315 "eos_interface.loci"
       output("cp") ;
#line 315 "eos_interface.loci"
    }
#line 315 "eos_interface.loci"
    void calculate(Entity _e_) { 
#line 316 "eos_interface.loci"
    L_cp_[_e_]=L_eos_state_[L_ci_[_e_]].cpt () ;
  }    void compute(const Loci::sequence &seq) { 
#line 317 "eos_interface.loci"
      do_loop(seq,this) ;
#line 317 "eos_interface.loci"
    }
#line 317 "eos_interface.loci"
} ;
#line 317 "eos_interface.loci"
Loci::register_rule<file_eos_interface015_1280810721m254> register_file_eos_interface015_1280810721m254 ;
#line 317 "eos_interface.loci"
}
#line 317 "eos_interface.loci"


  // $type e_internal store<real> 
  namespace {class file_eos_interface016_1280810721m254 : public Loci::pointwise_rule {
#line 320 "eos_interface.loci"
    Loci::const_store<EOS::State>  L_eos_state_ ; 
#line 320 "eos_interface.loci"
    Loci::store<real>  L_e_internal_ ; 
#line 320 "eos_interface.loci"
public:
#line 320 "eos_interface.loci"
    file_eos_interface016_1280810721m254() {
#line 320 "eos_interface.loci"
       name_store("eos_state",L_eos_state_) ;
#line 320 "eos_interface.loci"
       name_store("e_internal",L_e_internal_) ;
#line 320 "eos_interface.loci"
       input("eos_state") ;
#line 320 "eos_interface.loci"
       output("e_internal") ;
#line 320 "eos_interface.loci"
    }
#line 320 "eos_interface.loci"
    void calculate(Entity _e_) { 
#line 321 "eos_interface.loci"
    L_e_internal_[_e_]= L_eos_state_[_e_].energy () ;
  }    void compute(const Loci::sequence &seq) { 
#line 322 "eos_interface.loci"
      do_loop(seq,this) ;
#line 322 "eos_interface.loci"
    }
#line 322 "eos_interface.loci"
} ;
#line 322 "eos_interface.loci"
Loci::register_rule<file_eos_interface016_1280810721m254> register_file_eos_interface016_1280810721m254 ;
#line 322 "eos_interface.loci"
}
#line 322 "eos_interface.loci"


//$rule pointwise(dTdq<-eos,qvi,eos_state,eos_mixture_state),
//  prelude {
//  $dTdq.setVecSize($qvi->vectorSize()) ;
//} compute {
//  primitive_vector::dTdq($eos,$eos_state,$eos_mixture_state,$qvi,$dTdq) ;
//}

  // $type speciesEnthalpy storeVec<real> 
  namespace {class file_eos_interface017_1280810721m254 : public Loci::pointwise_rule {
#line 333 "eos_interface.loci"
    Loci::const_param<fluidPhysics::EOS>  L_eos_ ; 
#line 333 "eos_interface.loci"
    Loci::const_store<EOS::State>  L_eos_state_ ; 
#line 333 "eos_interface.loci"
    Loci::const_storeVec<real>  L_eos_mixture_state_ ; 
#line 333 "eos_interface.loci"
    Loci::storeVec<real>  L_speciesEnthalpy_ ; 
#line 333 "eos_interface.loci"
public:
#line 333 "eos_interface.loci"
    file_eos_interface017_1280810721m254() {
#line 333 "eos_interface.loci"
       name_store("eos",L_eos_) ;
#line 333 "eos_interface.loci"
       name_store("eos_state",L_eos_state_) ;
#line 333 "eos_interface.loci"
       name_store("eos_mixture_state",L_eos_mixture_state_) ;
#line 333 "eos_interface.loci"
       name_store("speciesEnthalpy",L_speciesEnthalpy_) ;
#line 333 "eos_interface.loci"
       input("eos,eos_state,eos_mixture_state") ;
#line 333 "eos_interface.loci"
       output("speciesEnthalpy") ;
#line 333 "eos_interface.loci"
    }
#line 333 "eos_interface.loci"
    void prelude(const Loci::sequence &seq) { 
    L_speciesEnthalpy_.setVecSize(L_eos_->numSpecies()) ;
  }    void calculate(Entity _e_) { 
#line 336 "eos_interface.loci"
    const int ns = L_eos_[_e_].numSpecies () ;
    for (int i =0;i <ns ;++i )
      L_speciesEnthalpy_[_e_][i ]= L_eos_[_e_].speciesEnthalpy (i ,L_eos_state_[_e_],L_eos_mixture_state_[_e_]) ;
  }    void compute(const Loci::sequence &seq) { 
#line 339 "eos_interface.loci"
      prelude(seq) ;
#line 339 "eos_interface.loci"
      do_loop(seq,this) ;
#line 339 "eos_interface.loci"
    }
#line 339 "eos_interface.loci"
} ;
#line 339 "eos_interface.loci"
Loci::register_rule<file_eos_interface017_1280810721m254> register_file_eos_interface017_1280810721m254 ;
#line 339 "eos_interface.loci"
}
#line 339 "eos_interface.loci"


  // $type speciesMolecularMass param<std::vector<real> > 
  namespace {class file_eos_interface018_1280810721m255 : public Loci::singleton_rule {
#line 342 "eos_interface.loci"
    Loci::const_param<fluidPhysics::EOS>  L_eos_ ; 
#line 342 "eos_interface.loci"
    Loci::param<std::vector<real> >  L_speciesMolecularMass_ ; 
#line 342 "eos_interface.loci"
public:
#line 342 "eos_interface.loci"
    file_eos_interface018_1280810721m255() {
#line 342 "eos_interface.loci"
       name_store("eos",L_eos_) ;
#line 342 "eos_interface.loci"
       name_store("speciesMolecularMass",L_speciesMolecularMass_) ;
#line 342 "eos_interface.loci"
       input("eos") ;
#line 342 "eos_interface.loci"
       output("speciesMolecularMass") ;
#line 342 "eos_interface.loci"
    }
#line 342 "eos_interface.loci"
    void compute(const Loci::sequence &seq) { 
    std::vector<real> smm ;
    for(int i=0;i<L_eos_->numSpecies();++i)
      smm.push_back(L_eos_->speciesMolecularMass(i)) ;
    (*L_speciesMolecularMass_)= smm ;
  }} ;
#line 347 "eos_interface.loci"
Loci::register_rule<file_eos_interface018_1280810721m255> register_file_eos_interface018_1280810721m255 ;
#line 347 "eos_interface.loci"
}
#line 347 "eos_interface.loci"


  // $type defaultMixtureFraction param<std::vector<real> > 
  namespace {class file_eos_interface019_1280810721m255 : public Loci::singleton_rule {
#line 350 "eos_interface.loci"
    Loci::const_param<fluidPhysics::EOS>  L_eos_ ; 
#line 350 "eos_interface.loci"
    Loci::param<std::vector<real> >  L_defaultMixtureFraction_ ; 
#line 350 "eos_interface.loci"
public:
#line 350 "eos_interface.loci"
    file_eos_interface019_1280810721m255() {
#line 350 "eos_interface.loci"
       name_store("eos",L_eos_) ;
#line 350 "eos_interface.loci"
       name_store("defaultMixtureFraction",L_defaultMixtureFraction_) ;
#line 350 "eos_interface.loci"
       input("eos") ;
#line 350 "eos_interface.loci"
       output("defaultMixtureFraction") ;
#line 350 "eos_interface.loci"
    }
#line 350 "eos_interface.loci"
    void compute(const Loci::sequence &seq) { 
    (*L_defaultMixtureFraction_)= L_eos_->getMixtureFractions() ;
  }} ;
#line 352 "eos_interface.loci"
Loci::register_rule<file_eos_interface019_1280810721m255> register_file_eos_interface019_1280810721m255 ;
#line 352 "eos_interface.loci"
}
#line 352 "eos_interface.loci"


}
