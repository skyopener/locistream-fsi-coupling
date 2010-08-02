#line 1 "FSI_CSDVarsInput.loci"
//-----------------------------------------------------------------------------
// Description: This file contains some of the basic rules to read in the vars file inputs related to the CSD and FSI.
//-----------------------------------------------------------------------------


// Loci includes.
#include <Loci.h>

// StreamUns includes.
#include "const.h"
#include "sciTypes.h"
#include "varsFileInputs.h"

#line 1 "FSI_CSDvariables.lh"
// $type CSDE1 param<real>  // Young's modulus 1
// $type CSDE2 param<real>  // Young's modulus 2
// $type CSDnu12 param<real>  // Poisson ratio 12
// $type CSDnu21 param<real>  // Poisson ratio 21
// $type CSDG12 param<real>  // Shear modulus 12
// $type CSDrhoStructure param<real>  // structure density
// $type CSDthicknessStructure param<real>  // structure thickness
// $type CSDintegrationScheme param<int>  // integrationScheme for the CSD-> 1:newmark, 2:generalized alpha
// $type CSDdelta param<real>  // CSD convergence criteria for the inner iteration
// $type CSDswitchStiffening param<int>  // 1: on, 2: off -> stiffening tangent stiffness matrix
//$type CSDexcitationType param<int> ; // 0: plunging, 1: flapping <- not used
// $type CSDflappingType param<int>  // 0: switch off, 1: sin, 2: cos, 3: 1-cos, 4: delay ( 1-exp(-4000*t^2))*sin
// $type CSDplungingType param<int>  // 0: switch off, 1: sin, 2: cos, 3: 1-cos
// $type CSDfrequency param<real>  // excitation frequency: ex. sin(2*pi*f*t) -> f (NOT omega!)
// $type CSDgenAlphaCoeff param<real>  // alpha for generalized alpha time integration scheme
// $type CSDnewmarkGammaCoeff param<real>  // Newmark beta coefficient1
// $type CSDnewmarkBetaCoeff param<real>  // Newmark beta coefficient2
// $type CSDdampingCoeff1 param<real>  // damping coefficient for the Rayleigh damping 1
// $type CSDdampingCoeff2 param<real>  // damping coefficient for the Rayleigh damping 2
// $type CSDplungeAmplitudeX param<real>  // plunge amplitude X
// $type CSDplungeAmplitudeY param<real>  // plunge amplitude Y
// $type CSDplungeAmplitudeZ param<real>  // plunge amplitude Z
// $type CSDflappingAmplitudeX param<real>  // flapping amplitude X <- pitch, different coordinate system in CSD, x <- spanwise
// $type CSDflappingAmplitudeY param<real>  // flapping amplitude Y <- flap, different coordinate system in CSD, y <- streamwise
// $type CSDflappingAmplitudeZ param<real>  // flapping amplitude Z <- lag, , different coordinate system in CSD, z <- -gravity 
// $type CSDMeshFilename param<string>  // CSD mesh filename
// $type CSDConnectivityFilename param<string>  // CSD connectivity filename
// $type CSDBCFilename param<string>  // CSD boundary conditions filename
// $type CSDstartingTimeStep param<int>  // 
// $type CSDtipNode param<int>  // Tip Node number inside CSD mesh
// $type CSDdimension param<int>  // CSD solver dimension constraint 2 -> 2D, 3 -> 3D
// $type CSD2dSpanCenter param<real>  // CSD solver dimension == 2D --> set the coordinate for the midspan point
// $type CSDYorigin param<real>  // CSD mesh assumed to be a plate with its Y (height) coordinate being this value. Then when interpolating to&from CFD the difference between the top and bottom CFD surfaces can be eliminated

// $type CFDMaxTotalInnerIterations param<int>  // Maximum allowed CFD iterations per time step

// $type FSICouplingMethod param<string>  // parameter to turn on the FSI coupling
// $type FSIRBFr param<real>  // FSI RBF interpolation r
// $type FSIRBFa param<real>  // FSI RBF interpolation a
// $type FSIRBFnr param<int>  // FSI RBF interpolation function nr
// $type FSIRBFMaxLinearSolverIterations param<int>  // FSI RBF interpolation number of maximum iterations for the linear solver
// $type FSIRBFTolerance param<real>  // FSI RBF interpolation tolerance for the linear solver
// $type FSICSD2CFDRBFr param<real>  // FSI CSD2CFD RBF interpolation r
// $type FSICSD2CFDRBFa param<real>  // FSI CSD2CFD RBF interpolation a
// $type FSICSD2CFDRBFnr param<int>  // FSI CSD2CFD RBF interpolation function nr
// $type FSICSD2CFDRBFMaxLinearSolverIterations param<int>  // FSI CSD2CFD RBF interpolation number of maximum iterations for the linear solver
// $type FSICSD2CFDRBFTolerance param<real>  // FSI RBF interpolation tolerance for the linear solver
// $type FSIIterationTolerance param<real>  // FSI inner iteration tolerance
// $type FSIIterationMinimum param<int>  // FSI minimum number of inner iterations 
// $type maxIterationsPerFSI param<int> 
// $type itfsi param<int>  // FSI inner iteration counter
// $type itfsi_ic param<int>  // FSI inner iteration counter initializer
// $type FSICoupling Constraint
// $type FSINLAMS Constraint
// $type FSIEULERBEAM Constraint
// $type FSI3DCONT Constraint
// $type FSIPRESCRIBED Constraint

  
// $type CSDNumberOfNodes param<int> 
// $type CSDNumberOfElements param<int> 
// $type CSDNumberOfBCs param<int> 
// $type CSDnodes_ic blackbox<ublas::matrix<real,ublas::column_major> >  // CSD initial 
// $type CSDnodes blackbox<ublas::matrix<real,ublas::column_major> >  // CSD current x
// $type CSDnodes_it blackbox<ublas::matrix<real,ublas::column_major> >  // CSD current x
// $type CSDnodesDisp_ic blackbox<ublas::matrix<real,ublas::column_major> >  // CSD initial generalized displacements
// $type CSDnodesDisp blackbox<ublas::matrix<real,ublas::column_major> >  // CSD current generalized displacements	
// $type CSDnodesVel_ic blackbox<ublas::matrix<real,ublas::column_major> >  // CSD initial generalized dot(displacements)
// $type CSDnodesVel blackbox<ublas::matrix<real,ublas::column_major> >  // CSD current generalized dot(displacements)
// $type CSDnodesAcc_ic blackbox<ublas::matrix<real,ublas::column_major> >  // CSD initial generalized ddot(displacements) 
// $type CSDnodesAcc blackbox<ublas::matrix<real,ublas::column_major> >  // CSD current generalized ddot(displacements)
// $type CSDConnectivity blackbox<ublas::matrix<int,ublas::column_major> >  // CSD connectivity matrix
// $type CSDBCdof param<vector<int> > 
// $type CSDBCZeroConstraint param<vector<real> > 

// $type CSDdisplacementsStar blackbox<ublas::matrix<real,ublas::column_major> >  // CSD displacements from NLAMS
// $type CSDdisplacementsOld blackbox<ublas::matrix<real,ublas::column_major> >  // CSD displacements from NLAMS, previous iteration
// $type CSDdisplacementsOldStar blackbox<ublas::matrix<real,ublas::column_major> >  // CSD displacements from NLAMS, previous iteration
// $type CSDdisplacementsOld_ic blackbox<ublas::matrix<real,ublas::column_major> >  // CSD displacements from NLAMS, previous iteration
// $type CSDnodesDispStar blackbox<ublas::matrix<real,ublas::column_major> >  // CSD generalized displacements from NLAMS
// $type CSDnodesVelStar blackbox<ublas::matrix<real,ublas::column_major> >  // CSD generalized dot(displacements) from NLAMS
// $type CSDnodesAccStar blackbox<ublas::matrix<real,ublas::column_major> >  // CSD generalized ddot(displacements) from NLAMS
// $type CSDForce blackbox<ublas::matrix<real,ublas::column_major> >  // Forces from CFD to CSD
// $type CSDForcePreStar blackbox<ublas::matrix<real,ublas::column_major> >  // CFD2CSD Force at the previous fsi iteration
// $type CSDForcePre blackbox<ublas::matrix<real,ublas::column_major> >  // CFD2CSD Force at the previous fsi iteration
// $type CSDForcePre_ic blackbox<ublas::matrix<real,ublas::column_major> >  // CFD2CSD Force at the previous fsi iteration initialization
// $type CSDnodesSysStar blackbox<boost::multi_array<real,4> >  // Nodal coordinates of each elements: from NLAMS
// $type CSDnodesSys blackbox<boost::multi_array<real,4> >  // Nodal coordinates of each elements
// $type CSDnodesSys_ic blackbox<boost::multi_array<real,4> >  // Nodal coordinates of each elements initialization
// $type CSDRBFweights blackbox<ublas::matrix<real,ublas::row_major> > 
	
// EULERBEAM
// $type CSDEulerXstart param<real> 
// $type CSDEulerXend param<real> 
// $type CSDEulerChord param<real> 
// $type CSDEulerXnum param<int> 
// $type CSDEulerAxis param<int>  // beam direction: 0->x, 1->y, 2->z, displacement probably always in y direction
// $type CSDEulerBeamDirection param<int> 
// $type CSDEulerSpanDirection param<int> 
// $type CSDxStar blackbox<ublas::vector<real> > 
// $type CSDxdotStar blackbox<ublas::vector<real> > 
// $type CSDxddotStar blackbox<ublas::vector<real> > 
// $type CSDx_ic blackbox<ublas::vector<real> > 
// $type CSDxdot_ic blackbox<ublas::vector<real> > 
// $type CSDxddot_ic blackbox<ublas::vector<real> > 
// $type CSDx blackbox<ublas::vector<real> > 
// $type CSDxdot blackbox<ublas::vector<real> > 
// $type CSDxddot blackbox<ublas::vector<real> > 

#line 14 "FSI_CSDVarsInput.loci"


namespace streamUns {

  //---------------------------------------------------------------------------
  // Setup Rules.
 
 
 
  namespace {class file_FSI_CSDVarsInput000_1279672461m181 : public Loci::default_rule {
#line 23 "FSI_CSDVarsInput.loci"
    Loci::param<int>  L_CFDMaxTotalInnerIterations_ ; 
#line 23 "FSI_CSDVarsInput.loci"
public:
#line 23 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput000_1279672461m181() {
#line 23 "FSI_CSDVarsInput.loci"
       name_store("CFDMaxTotalInnerIterations",L_CFDMaxTotalInnerIterations_) ;
#line 23 "FSI_CSDVarsInput.loci"
       output("CFDMaxTotalInnerIterations") ;
#line 23 "FSI_CSDVarsInput.loci"
    }
#line 23 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CFDMaxTotalInnerIterations_)= 10000; // aluminum
  }} ;
#line 25 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput000_1279672461m181> register_file_FSI_CSDVarsInput000_1279672461m181 ;
#line 25 "FSI_CSDVarsInput.loci"
}
#line 25 "FSI_CSDVarsInput.loci"
 
 
  namespace {class file_FSI_CSDVarsInput001_1279672461m181 : public Loci::default_rule {
#line 27 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDE1_ ; 
#line 27 "FSI_CSDVarsInput.loci"
public:
#line 27 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput001_1279672461m181() {
#line 27 "FSI_CSDVarsInput.loci"
       name_store("CSDE1",L_CSDE1_) ;
#line 27 "FSI_CSDVarsInput.loci"
       output("CSDE1") ;
#line 27 "FSI_CSDVarsInput.loci"
    }
#line 27 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDE1_)= 70.e9; // aluminum
  }} ;
#line 29 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput001_1279672461m181> register_file_FSI_CSDVarsInput001_1279672461m181 ;
#line 29 "FSI_CSDVarsInput.loci"
}
#line 29 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput002_1279672461m182 : public Loci::default_rule {
#line 31 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDE2_ ; 
#line 31 "FSI_CSDVarsInput.loci"
public:
#line 31 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput002_1279672461m182() {
#line 31 "FSI_CSDVarsInput.loci"
       name_store("CSDE2",L_CSDE2_) ;
#line 31 "FSI_CSDVarsInput.loci"
       output("CSDE2") ;
#line 31 "FSI_CSDVarsInput.loci"
    }
#line 31 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDE2_)= 70.e9; // aluminum
  }} ;
#line 33 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput002_1279672461m182> register_file_FSI_CSDVarsInput002_1279672461m182 ;
#line 33 "FSI_CSDVarsInput.loci"
}
#line 33 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput003_1279672461m182 : public Loci::default_rule {
#line 35 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDnu12_ ; 
#line 35 "FSI_CSDVarsInput.loci"
public:
#line 35 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput003_1279672461m182() {
#line 35 "FSI_CSDVarsInput.loci"
       name_store("CSDnu12",L_CSDnu12_) ;
#line 35 "FSI_CSDVarsInput.loci"
       output("CSDnu12") ;
#line 35 "FSI_CSDVarsInput.loci"
    }
#line 35 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDnu12_)= 0.3; // aluminum
  }} ;
#line 37 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput003_1279672461m182> register_file_FSI_CSDVarsInput003_1279672461m182 ;
#line 37 "FSI_CSDVarsInput.loci"
}
#line 37 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput004_1279672461m182 : public Loci::default_rule {
#line 39 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDnu21_ ; 
#line 39 "FSI_CSDVarsInput.loci"
public:
#line 39 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput004_1279672461m182() {
#line 39 "FSI_CSDVarsInput.loci"
       name_store("CSDnu21",L_CSDnu21_) ;
#line 39 "FSI_CSDVarsInput.loci"
       output("CSDnu21") ;
#line 39 "FSI_CSDVarsInput.loci"
    }
#line 39 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDnu21_)= 0.3; // aluminum
  }} ;
#line 41 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput004_1279672461m182> register_file_FSI_CSDVarsInput004_1279672461m182 ;
#line 41 "FSI_CSDVarsInput.loci"
}
#line 41 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput005_1279672461m182 : public Loci::default_rule {
#line 43 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDG12_ ; 
#line 43 "FSI_CSDVarsInput.loci"
public:
#line 43 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput005_1279672461m182() {
#line 43 "FSI_CSDVarsInput.loci"
       name_store("CSDG12",L_CSDG12_) ;
#line 43 "FSI_CSDVarsInput.loci"
       output("CSDG12") ;
#line 43 "FSI_CSDVarsInput.loci"
    }
#line 43 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDG12_)= 27.e9; // aluminum
  }} ;
#line 45 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput005_1279672461m182> register_file_FSI_CSDVarsInput005_1279672461m182 ;
#line 45 "FSI_CSDVarsInput.loci"
}
#line 45 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput006_1279672461m182 : public Loci::default_rule {
#line 47 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDrhoStructure_ ; 
#line 47 "FSI_CSDVarsInput.loci"
public:
#line 47 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput006_1279672461m182() {
#line 47 "FSI_CSDVarsInput.loci"
       name_store("CSDrhoStructure",L_CSDrhoStructure_) ;
#line 47 "FSI_CSDVarsInput.loci"
       output("CSDrhoStructure") ;
#line 47 "FSI_CSDVarsInput.loci"
    }
#line 47 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDrhoStructure_)= 2700.; // aluminum
  }} ;
#line 49 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput006_1279672461m182> register_file_FSI_CSDVarsInput006_1279672461m182 ;
#line 49 "FSI_CSDVarsInput.loci"
}
#line 49 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput007_1279672461m183 : public Loci::default_rule {
#line 51 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDthicknessStructure_ ; 
#line 51 "FSI_CSDVarsInput.loci"
public:
#line 51 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput007_1279672461m183() {
#line 51 "FSI_CSDVarsInput.loci"
       name_store("CSDthicknessStructure",L_CSDthicknessStructure_) ;
#line 51 "FSI_CSDVarsInput.loci"
       output("CSDthicknessStructure") ;
#line 51 "FSI_CSDVarsInput.loci"
    }
#line 51 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDthicknessStructure_)= 1.e-3; // 0.001
  }} ;
#line 53 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput007_1279672461m183> register_file_FSI_CSDVarsInput007_1279672461m183 ;
#line 53 "FSI_CSDVarsInput.loci"
}
#line 53 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput008_1279672461m183 : public Loci::default_rule {
#line 55 "FSI_CSDVarsInput.loci"
    Loci::param<int>  L_CSDintegrationScheme_ ; 
#line 55 "FSI_CSDVarsInput.loci"
public:
#line 55 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput008_1279672461m183() {
#line 55 "FSI_CSDVarsInput.loci"
       name_store("CSDintegrationScheme",L_CSDintegrationScheme_) ;
#line 55 "FSI_CSDVarsInput.loci"
       output("CSDintegrationScheme") ;
#line 55 "FSI_CSDVarsInput.loci"
    }
#line 55 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDintegrationScheme_)= 1; // newmark
  }} ;
#line 57 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput008_1279672461m183> register_file_FSI_CSDVarsInput008_1279672461m183 ;
#line 57 "FSI_CSDVarsInput.loci"
}
#line 57 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput009_1279672461m183 : public Loci::default_rule {
#line 59 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDdelta_ ; 
#line 59 "FSI_CSDVarsInput.loci"
public:
#line 59 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput009_1279672461m183() {
#line 59 "FSI_CSDVarsInput.loci"
       name_store("CSDdelta",L_CSDdelta_) ;
#line 59 "FSI_CSDVarsInput.loci"
       output("CSDdelta") ;
#line 59 "FSI_CSDVarsInput.loci"
    }
#line 59 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDdelta_)= 1.e-4; // default
  }} ;
#line 61 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput009_1279672461m183> register_file_FSI_CSDVarsInput009_1279672461m183 ;
#line 61 "FSI_CSDVarsInput.loci"
}
#line 61 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput010_1279672461m183 : public Loci::default_rule {
#line 63 "FSI_CSDVarsInput.loci"
    Loci::param<int>  L_CSDswitchStiffening_ ; 
#line 63 "FSI_CSDVarsInput.loci"
public:
#line 63 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput010_1279672461m183() {
#line 63 "FSI_CSDVarsInput.loci"
       name_store("CSDswitchStiffening",L_CSDswitchStiffening_) ;
#line 63 "FSI_CSDVarsInput.loci"
       output("CSDswitchStiffening") ;
#line 63 "FSI_CSDVarsInput.loci"
    }
#line 63 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDswitchStiffening_)= 1; // turned on
  }} ;
#line 65 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput010_1279672461m183> register_file_FSI_CSDVarsInput010_1279672461m183 ;
#line 65 "FSI_CSDVarsInput.loci"
}
#line 65 "FSI_CSDVarsInput.loci"


// Not used    
//  $rule default(CSDexcitationType) {
//  	$CSDexcitationType = 0; // flapping
//  }  
  
  namespace {class file_FSI_CSDVarsInput011_1279672461m183 : public Loci::default_rule {
#line 72 "FSI_CSDVarsInput.loci"
    Loci::param<int>  L_CSDflappingType_ ; 
#line 72 "FSI_CSDVarsInput.loci"
public:
#line 72 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput011_1279672461m183() {
#line 72 "FSI_CSDVarsInput.loci"
       name_store("CSDflappingType",L_CSDflappingType_) ;
#line 72 "FSI_CSDVarsInput.loci"
       output("CSDflappingType") ;
#line 72 "FSI_CSDVarsInput.loci"
    }
#line 72 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDflappingType_)= 1; // sin
  }} ;
#line 74 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput011_1279672461m183> register_file_FSI_CSDVarsInput011_1279672461m183 ;
#line 74 "FSI_CSDVarsInput.loci"
}
#line 74 "FSI_CSDVarsInput.loci"
  
  
  namespace {class file_FSI_CSDVarsInput012_1279672461m184 : public Loci::default_rule {
#line 76 "FSI_CSDVarsInput.loci"
    Loci::param<int>  L_CSDplungingType_ ; 
#line 76 "FSI_CSDVarsInput.loci"
public:
#line 76 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput012_1279672461m184() {
#line 76 "FSI_CSDVarsInput.loci"
       name_store("CSDplungingType",L_CSDplungingType_) ;
#line 76 "FSI_CSDVarsInput.loci"
       output("CSDplungingType") ;
#line 76 "FSI_CSDVarsInput.loci"
    }
#line 76 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDplungingType_)= 1; // sin
  }} ;
#line 78 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput012_1279672461m184> register_file_FSI_CSDVarsInput012_1279672461m184 ;
#line 78 "FSI_CSDVarsInput.loci"
}
#line 78 "FSI_CSDVarsInput.loci"
  
  
  namespace {class file_FSI_CSDVarsInput013_1279672461m184 : public Loci::default_rule {
#line 80 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDgenAlphaCoeff_ ; 
#line 80 "FSI_CSDVarsInput.loci"
public:
#line 80 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput013_1279672461m184() {
#line 80 "FSI_CSDVarsInput.loci"
       name_store("CSDgenAlphaCoeff",L_CSDgenAlphaCoeff_) ;
#line 80 "FSI_CSDVarsInput.loci"
       output("CSDgenAlphaCoeff") ;
#line 80 "FSI_CSDVarsInput.loci"
    }
#line 80 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDgenAlphaCoeff_)= 0.4; // 0.4 
  }} ;
#line 82 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput013_1279672461m184> register_file_FSI_CSDVarsInput013_1279672461m184 ;
#line 82 "FSI_CSDVarsInput.loci"
}
#line 82 "FSI_CSDVarsInput.loci"
  
  
  namespace {class file_FSI_CSDVarsInput014_1279672461m184 : public Loci::default_rule {
#line 84 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDnewmarkGammaCoeff_ ; 
#line 84 "FSI_CSDVarsInput.loci"
public:
#line 84 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput014_1279672461m184() {
#line 84 "FSI_CSDVarsInput.loci"
       name_store("CSDnewmarkGammaCoeff",L_CSDnewmarkGammaCoeff_) ;
#line 84 "FSI_CSDVarsInput.loci"
       output("CSDnewmarkGammaCoeff") ;
#line 84 "FSI_CSDVarsInput.loci"
    }
#line 84 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDnewmarkGammaCoeff_)= 0.5; // linear
  }} ;
#line 86 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput014_1279672461m184> register_file_FSI_CSDVarsInput014_1279672461m184 ;
#line 86 "FSI_CSDVarsInput.loci"
}
#line 86 "FSI_CSDVarsInput.loci"
  
  
  namespace {class file_FSI_CSDVarsInput015_1279672461m184 : public Loci::default_rule {
#line 88 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDnewmarkBetaCoeff_ ; 
#line 88 "FSI_CSDVarsInput.loci"
public:
#line 88 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput015_1279672461m184() {
#line 88 "FSI_CSDVarsInput.loci"
       name_store("CSDnewmarkBetaCoeff",L_CSDnewmarkBetaCoeff_) ;
#line 88 "FSI_CSDVarsInput.loci"
       output("CSDnewmarkBetaCoeff") ;
#line 88 "FSI_CSDVarsInput.loci"
    }
#line 88 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDnewmarkBetaCoeff_)= 0.25; // linear
  }} ;
#line 90 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput015_1279672461m184> register_file_FSI_CSDVarsInput015_1279672461m184 ;
#line 90 "FSI_CSDVarsInput.loci"
}
#line 90 "FSI_CSDVarsInput.loci"
  
  
  namespace {class file_FSI_CSDVarsInput016_1279672461m184 : public Loci::default_rule {
#line 92 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDdampingCoeff1_ ; 
#line 92 "FSI_CSDVarsInput.loci"
public:
#line 92 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput016_1279672461m184() {
#line 92 "FSI_CSDVarsInput.loci"
       name_store("CSDdampingCoeff1",L_CSDdampingCoeff1_) ;
#line 92 "FSI_CSDVarsInput.loci"
       output("CSDdampingCoeff1") ;
#line 92 "FSI_CSDVarsInput.loci"
    }
#line 92 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDdampingCoeff1_)= 0.; // no damping
  }} ;
#line 94 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput016_1279672461m184> register_file_FSI_CSDVarsInput016_1279672461m184 ;
#line 94 "FSI_CSDVarsInput.loci"
}
#line 94 "FSI_CSDVarsInput.loci"


  namespace {class file_FSI_CSDVarsInput017_1279672461m184 : public Loci::default_rule {
#line 96 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDfrequency_ ; 
#line 96 "FSI_CSDVarsInput.loci"
public:
#line 96 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput017_1279672461m184() {
#line 96 "FSI_CSDVarsInput.loci"
       name_store("CSDfrequency",L_CSDfrequency_) ;
#line 96 "FSI_CSDVarsInput.loci"
       output("CSDfrequency") ;
#line 96 "FSI_CSDVarsInput.loci"
    }
#line 96 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDfrequency_)= 1.; // =f, NOT omega!
  }} ;
#line 98 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput017_1279672461m184> register_file_FSI_CSDVarsInput017_1279672461m184 ;
#line 98 "FSI_CSDVarsInput.loci"
}
#line 98 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput018_1279672461m185 : public Loci::default_rule {
#line 100 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDdampingCoeff2_ ; 
#line 100 "FSI_CSDVarsInput.loci"
public:
#line 100 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput018_1279672461m185() {
#line 100 "FSI_CSDVarsInput.loci"
       name_store("CSDdampingCoeff2",L_CSDdampingCoeff2_) ;
#line 100 "FSI_CSDVarsInput.loci"
       output("CSDdampingCoeff2") ;
#line 100 "FSI_CSDVarsInput.loci"
    }
#line 100 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDdampingCoeff2_)= 0.; // no damping
  }} ;
#line 102 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput018_1279672461m185> register_file_FSI_CSDVarsInput018_1279672461m185 ;
#line 102 "FSI_CSDVarsInput.loci"
}
#line 102 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput019_1279672461m185 : public Loci::default_rule {
#line 104 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDplungeAmplitudeX_ ; 
#line 104 "FSI_CSDVarsInput.loci"
public:
#line 104 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput019_1279672461m185() {
#line 104 "FSI_CSDVarsInput.loci"
       name_store("CSDplungeAmplitudeX",L_CSDplungeAmplitudeX_) ;
#line 104 "FSI_CSDVarsInput.loci"
       output("CSDplungeAmplitudeX") ;
#line 104 "FSI_CSDVarsInput.loci"
    }
#line 104 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDplungeAmplitudeX_)= 0.; // no plunging
  }} ;
#line 106 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput019_1279672461m185> register_file_FSI_CSDVarsInput019_1279672461m185 ;
#line 106 "FSI_CSDVarsInput.loci"
}
#line 106 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput020_1279672461m185 : public Loci::default_rule {
#line 108 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDplungeAmplitudeY_ ; 
#line 108 "FSI_CSDVarsInput.loci"
public:
#line 108 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput020_1279672461m185() {
#line 108 "FSI_CSDVarsInput.loci"
       name_store("CSDplungeAmplitudeY",L_CSDplungeAmplitudeY_) ;
#line 108 "FSI_CSDVarsInput.loci"
       output("CSDplungeAmplitudeY") ;
#line 108 "FSI_CSDVarsInput.loci"
    }
#line 108 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDplungeAmplitudeY_)= 0.; // no plunging
  }} ;
#line 110 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput020_1279672461m185> register_file_FSI_CSDVarsInput020_1279672461m185 ;
#line 110 "FSI_CSDVarsInput.loci"
}
#line 110 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput021_1279672461m185 : public Loci::default_rule {
#line 112 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDplungeAmplitudeZ_ ; 
#line 112 "FSI_CSDVarsInput.loci"
public:
#line 112 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput021_1279672461m185() {
#line 112 "FSI_CSDVarsInput.loci"
       name_store("CSDplungeAmplitudeZ",L_CSDplungeAmplitudeZ_) ;
#line 112 "FSI_CSDVarsInput.loci"
       output("CSDplungeAmplitudeZ") ;
#line 112 "FSI_CSDVarsInput.loci"
    }
#line 112 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDplungeAmplitudeZ_)= 0.; // no plunging
  }} ;
#line 114 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput021_1279672461m185> register_file_FSI_CSDVarsInput021_1279672461m185 ;
#line 114 "FSI_CSDVarsInput.loci"
}
#line 114 "FSI_CSDVarsInput.loci"

   
  namespace {class file_FSI_CSDVarsInput022_1279672461m185 : public Loci::default_rule {
#line 116 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDflappingAmplitudeX_ ; 
#line 116 "FSI_CSDVarsInput.loci"
public:
#line 116 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput022_1279672461m185() {
#line 116 "FSI_CSDVarsInput.loci"
       name_store("CSDflappingAmplitudeX",L_CSDflappingAmplitudeX_) ;
#line 116 "FSI_CSDVarsInput.loci"
       output("CSDflappingAmplitudeX") ;
#line 116 "FSI_CSDVarsInput.loci"
    }
#line 116 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDflappingAmplitudeX_)= 15.; // 15 deg
  }} ;
#line 118 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput022_1279672461m185> register_file_FSI_CSDVarsInput022_1279672461m185 ;
#line 118 "FSI_CSDVarsInput.loci"
}
#line 118 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput023_1279672461m186 : public Loci::default_rule {
#line 120 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDflappingAmplitudeY_ ; 
#line 120 "FSI_CSDVarsInput.loci"
public:
#line 120 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput023_1279672461m186() {
#line 120 "FSI_CSDVarsInput.loci"
       name_store("CSDflappingAmplitudeY",L_CSDflappingAmplitudeY_) ;
#line 120 "FSI_CSDVarsInput.loci"
       output("CSDflappingAmplitudeY") ;
#line 120 "FSI_CSDVarsInput.loci"
    }
#line 120 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDflappingAmplitudeY_)= 0.; // 0.
  }} ;
#line 122 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput023_1279672461m186> register_file_FSI_CSDVarsInput023_1279672461m186 ;
#line 122 "FSI_CSDVarsInput.loci"
}
#line 122 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput024_1279672461m186 : public Loci::default_rule {
#line 124 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDflappingAmplitudeZ_ ; 
#line 124 "FSI_CSDVarsInput.loci"
public:
#line 124 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput024_1279672461m186() {
#line 124 "FSI_CSDVarsInput.loci"
       name_store("CSDflappingAmplitudeZ",L_CSDflappingAmplitudeZ_) ;
#line 124 "FSI_CSDVarsInput.loci"
       output("CSDflappingAmplitudeZ") ;
#line 124 "FSI_CSDVarsInput.loci"
    }
#line 124 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDflappingAmplitudeZ_)= 0.; // 0.
  }} ;
#line 126 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput024_1279672461m186> register_file_FSI_CSDVarsInput024_1279672461m186 ;
#line 126 "FSI_CSDVarsInput.loci"
}
#line 126 "FSI_CSDVarsInput.loci"

    
  namespace {class file_FSI_CSDVarsInput025_1279672461m186 : public Loci::default_rule {
#line 128 "FSI_CSDVarsInput.loci"
    Loci::param<int>  L_CSDstartingTimeStep_ ; 
#line 128 "FSI_CSDVarsInput.loci"
public:
#line 128 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput025_1279672461m186() {
#line 128 "FSI_CSDVarsInput.loci"
       name_store("CSDstartingTimeStep",L_CSDstartingTimeStep_) ;
#line 128 "FSI_CSDVarsInput.loci"
       output("CSDstartingTimeStep") ;
#line 128 "FSI_CSDVarsInput.loci"
    }
#line 128 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDstartingTimeStep_)= 0; // 0.
  }} ;
#line 130 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput025_1279672461m186> register_file_FSI_CSDVarsInput025_1279672461m186 ;
#line 130 "FSI_CSDVarsInput.loci"
}
#line 130 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput026_1279672461m186 : public Loci::default_rule {
#line 132 "FSI_CSDVarsInput.loci"
    Loci::param<int>  L_CSDdimension_ ; 
#line 132 "FSI_CSDVarsInput.loci"
public:
#line 132 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput026_1279672461m186() {
#line 132 "FSI_CSDVarsInput.loci"
       name_store("CSDdimension",L_CSDdimension_) ;
#line 132 "FSI_CSDVarsInput.loci"
       output("CSDdimension") ;
#line 132 "FSI_CSDVarsInput.loci"
    }
#line 132 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDdimension_)= 3; // 3 -> 3D, 2 -> 2D
  }} ;
#line 134 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput026_1279672461m186> register_file_FSI_CSDVarsInput026_1279672461m186 ;
#line 134 "FSI_CSDVarsInput.loci"
}
#line 134 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput027_1279672461m186 : public Loci::default_rule {
#line 136 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDYorigin_ ; 
#line 136 "FSI_CSDVarsInput.loci"
public:
#line 136 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput027_1279672461m186() {
#line 136 "FSI_CSDVarsInput.loci"
       name_store("CSDYorigin",L_CSDYorigin_) ;
#line 136 "FSI_CSDVarsInput.loci"
       output("CSDYorigin") ;
#line 136 "FSI_CSDVarsInput.loci"
    }
#line 136 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDYorigin_)= 0.0; // CSD mesh assumed to be a plate with its Y (height) coordinate being this value. Then when interpolating to&from CFD the difference between the top and bottom CFD surfaces can be eliminated
  }} ;
#line 138 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput027_1279672461m186> register_file_FSI_CSDVarsInput027_1279672461m186 ;
#line 138 "FSI_CSDVarsInput.loci"
}
#line 138 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput028_1279672461m187 : public Loci::default_rule {
#line 140 "FSI_CSDVarsInput.loci"
    Loci::param<int>  L_CSDtipNode_ ; 
#line 140 "FSI_CSDVarsInput.loci"
public:
#line 140 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput028_1279672461m187() {
#line 140 "FSI_CSDVarsInput.loci"
       name_store("CSDtipNode",L_CSDtipNode_) ;
#line 140 "FSI_CSDVarsInput.loci"
       output("CSDtipNode") ;
#line 140 "FSI_CSDVarsInput.loci"
    }
#line 140 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDtipNode_)= 1; // 3 -> 3D, 2 -> 2D
  }} ;
#line 142 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput028_1279672461m187> register_file_FSI_CSDVarsInput028_1279672461m187 ;
#line 142 "FSI_CSDVarsInput.loci"
}
#line 142 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput029_1279672461m187 : public Loci::default_rule {
#line 144 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSD2dSpanCenter_ ; 
#line 144 "FSI_CSDVarsInput.loci"
public:
#line 144 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput029_1279672461m187() {
#line 144 "FSI_CSDVarsInput.loci"
       name_store("CSD2dSpanCenter",L_CSD2dSpanCenter_) ;
#line 144 "FSI_CSDVarsInput.loci"
       output("CSD2dSpanCenter") ;
#line 144 "FSI_CSDVarsInput.loci"
    }
#line 144 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	//if ($CSDdimension==2) {
  		(*L_CSD2dSpanCenter_)= 0.5; // 3 -> 3D, 2 -> 2D
 // 	} else {
  //		cerr << "CSDdimension needs to be set to 2" << endl ;
  //		Loci::Abort() ;
  //	}
  }} ;
#line 151 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput029_1279672461m187> register_file_FSI_CSDVarsInput029_1279672461m187 ;
#line 151 "FSI_CSDVarsInput.loci"
}
#line 151 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput030_1279672461m187 : public Loci::default_rule {
#line 153 "FSI_CSDVarsInput.loci"
    Loci::param<string>  L_FSICouplingMethod_ ; 
#line 153 "FSI_CSDVarsInput.loci"
public:
#line 153 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput030_1279672461m187() {
#line 153 "FSI_CSDVarsInput.loci"
       name_store("FSICouplingMethod",L_FSICouplingMethod_) ;
#line 153 "FSI_CSDVarsInput.loci"
       output("FSICouplingMethod") ;
#line 153 "FSI_CSDVarsInput.loci"
    }
#line 153 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) {  // create Constraint for the FSICoupling
  	// nothing now
  	(*L_FSICouplingMethod_)= "none" ;
  	
  }} ;
#line 157 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput030_1279672461m187> register_file_FSI_CSDVarsInput030_1279672461m187 ;
#line 157 "FSI_CSDVarsInput.loci"
}
#line 157 "FSI_CSDVarsInput.loci"

   
  namespace {class file_FSI_CSDVarsInput031_1279672461m187 : public Loci::constraint_rule {
#line 159 "FSI_CSDVarsInput.loci"
    Loci::const_param<string>  L_FSICouplingMethod_ ; 
#line 159 "FSI_CSDVarsInput.loci"
    Loci::Constraint L_FSICoupling_ ; 
#line 159 "FSI_CSDVarsInput.loci"
public:
#line 159 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput031_1279672461m187() {
#line 159 "FSI_CSDVarsInput.loci"
       name_store("FSICouplingMethod",L_FSICouplingMethod_) ;
#line 159 "FSI_CSDVarsInput.loci"
       name_store("FSICoupling",L_FSICoupling_) ;
#line 159 "FSI_CSDVarsInput.loci"
       input("FSICouplingMethod") ;
#line 159 "FSI_CSDVarsInput.loci"
       output("FSICoupling") ;
#line 159 "FSI_CSDVarsInput.loci"
    }
#line 159 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	L_FSICoupling_= EMPTY ;
  	if ((*L_FSICouplingMethod_)== "none") {
  		 if (Loci::MPI_rank==0) cout << "FSI Coupling off: " << (*L_FSICouplingMethod_)<< endl;
  	} else {
  		L_FSICoupling_= ~EMPTY ;
  		if (Loci::MPI_rank==0) cout << "FSI Coupling on: " << (*L_FSICouplingMethod_)<< endl;
  	}		
  }} ;
#line 167 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput031_1279672461m187> register_file_FSI_CSDVarsInput031_1279672461m187 ;
#line 167 "FSI_CSDVarsInput.loci"
}
#line 167 "FSI_CSDVarsInput.loci"

   
  namespace {class file_FSI_CSDVarsInput032_1279672461m188 : public Loci::constraint_rule {
#line 169 "FSI_CSDVarsInput.loci"
    Loci::const_param<string>  L_FSICouplingMethod_ ; 
#line 169 "FSI_CSDVarsInput.loci"
    Loci::Constraint L_FSINLAMS_ ; 
#line 169 "FSI_CSDVarsInput.loci"
    Loci::Constraint L_FSIEULERBEAM_ ; 
#line 169 "FSI_CSDVarsInput.loci"
    Loci::Constraint L_FSI3DCONT_ ; 
#line 169 "FSI_CSDVarsInput.loci"
    Loci::Constraint L_FSIPRESCRIBED_ ; 
#line 169 "FSI_CSDVarsInput.loci"
public:
#line 169 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput032_1279672461m188() {
#line 169 "FSI_CSDVarsInput.loci"
       name_store("FSICouplingMethod",L_FSICouplingMethod_) ;
#line 169 "FSI_CSDVarsInput.loci"
       name_store("FSINLAMS",L_FSINLAMS_) ;
#line 169 "FSI_CSDVarsInput.loci"
       name_store("FSIEULERBEAM",L_FSIEULERBEAM_) ;
#line 169 "FSI_CSDVarsInput.loci"
       name_store("FSI3DCONT",L_FSI3DCONT_) ;
#line 169 "FSI_CSDVarsInput.loci"
       name_store("FSIPRESCRIBED",L_FSIPRESCRIBED_) ;
#line 169 "FSI_CSDVarsInput.loci"
       input("FSICouplingMethod") ;
#line 169 "FSI_CSDVarsInput.loci"
       output("FSINLAMS") ;
#line 169 "FSI_CSDVarsInput.loci"
       output("FSIEULERBEAM") ;
#line 169 "FSI_CSDVarsInput.loci"
       output("FSI3DCONT") ;
#line 169 "FSI_CSDVarsInput.loci"
       output("FSIPRESCRIBED") ;
#line 169 "FSI_CSDVarsInput.loci"
    }
#line 169 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	if ((*L_FSICouplingMethod_)== "none") {
  		 if (Loci::MPI_rank==0) cout << "FSI Coupling off: " << (*L_FSICouplingMethod_)<< endl;
  	} else {
  		if (Loci::MPI_rank==0) cout << "FSI Coupling on: " << (*L_FSICouplingMethod_)<< endl;
			if ((*L_FSICouplingMethod_)== "NLAMS") {
				L_FSINLAMS_= ~EMPTY ;
				L_FSIEULERBEAM_= EMPTY ;
				L_FSI3DCONT_= EMPTY ;
				L_FSIPRESCRIBED_= EMPTY ;
			} else if ((*L_FSICouplingMethod_)== "EULERBEAM") {
				L_FSINLAMS_= EMPTY ;
				L_FSIEULERBEAM_= ~EMPTY ;
				L_FSI3DCONT_= EMPTY ;
				L_FSIPRESCRIBED_= EMPTY ;
			} else if ((*L_FSICouplingMethod_)== "3DCONT") {
				L_FSINLAMS_= EMPTY ;
				L_FSIEULERBEAM_= EMPTY ;
				L_FSI3DCONT_= ~EMPTY ;
				L_FSIPRESCRIBED_= EMPTY ;
			} else if ((*L_FSICouplingMethod_)== "PRESCRIBED") {
				L_FSINLAMS_= EMPTY ;
				L_FSIEULERBEAM_= EMPTY ;
				L_FSI3DCONT_= EMPTY ;
				L_FSIPRESCRIBED_= ~EMPTY ;
			}	else {
				cerr << "FSI Coupling method not recognized" << endl ;
				cerr << "Available FSI Coupling methods are none, NLAMS, EULERBEAM, 3DCONT, PRESCRIBED" << endl ;
				Loci::Abort();  
			}				
		}
  }} ;
#line 200 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput032_1279672461m188> register_file_FSI_CSDVarsInput032_1279672461m188 ;
#line 200 "FSI_CSDVarsInput.loci"
}
#line 200 "FSI_CSDVarsInput.loci"
   
   
  namespace {class file_FSI_CSDVarsInput033_1279672461m188 : public Loci::default_rule {
#line 202 "FSI_CSDVarsInput.loci"
    Loci::param<string>  L_CSDMeshFilename_ ; 
#line 202 "FSI_CSDVarsInput.loci"
public:
#line 202 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput033_1279672461m188() {
#line 202 "FSI_CSDVarsInput.loci"
       name_store("CSDMeshFilename",L_CSDMeshFilename_) ;
#line 202 "FSI_CSDVarsInput.loci"
       output("CSDMeshFilename") ;
#line 202 "FSI_CSDVarsInput.loci"
    }
#line 202 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
 		(*L_CSDMeshFilename_)="trimesh.dat" ;
 	}} ;
#line 204 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput033_1279672461m188> register_file_FSI_CSDVarsInput033_1279672461m188 ;
#line 204 "FSI_CSDVarsInput.loci"
}
#line 204 "FSI_CSDVarsInput.loci"

 	
 	namespace {class file_FSI_CSDVarsInput034_1279672461m188 : public Loci::default_rule {
#line 206 "FSI_CSDVarsInput.loci"
    Loci::param<string>  L_CSDConnectivityFilename_ ; 
#line 206 "FSI_CSDVarsInput.loci"
public:
#line 206 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput034_1279672461m188() {
#line 206 "FSI_CSDVarsInput.loci"
       name_store("CSDConnectivityFilename",L_CSDConnectivityFilename_) ;
#line 206 "FSI_CSDVarsInput.loci"
       output("CSDConnectivityFilename") ;
#line 206 "FSI_CSDVarsInput.loci"
    }
#line 206 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
 		(*L_CSDConnectivityFilename_)="connect.dat" ;
 	}} ;
#line 208 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput034_1279672461m188> register_file_FSI_CSDVarsInput034_1279672461m188 ;
#line 208 "FSI_CSDVarsInput.loci"
}
#line 208 "FSI_CSDVarsInput.loci"

 	
 	namespace {class file_FSI_CSDVarsInput035_1279672461m189 : public Loci::default_rule {
#line 210 "FSI_CSDVarsInput.loci"
    Loci::param<string>  L_CSDBCFilename_ ; 
#line 210 "FSI_CSDVarsInput.loci"
public:
#line 210 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput035_1279672461m189() {
#line 210 "FSI_CSDVarsInput.loci"
       name_store("CSDBCFilename",L_CSDBCFilename_) ;
#line 210 "FSI_CSDVarsInput.loci"
       output("CSDBCFilename") ;
#line 210 "FSI_CSDVarsInput.loci"
    }
#line 210 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
 		(*L_CSDBCFilename_)="bc.dat" ;
 	}} ;
#line 212 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput035_1279672461m189> register_file_FSI_CSDVarsInput035_1279672461m189 ;
#line 212 "FSI_CSDVarsInput.loci"
}
#line 212 "FSI_CSDVarsInput.loci"

 	
  namespace {class file_FSI_CSDVarsInput036_1279672461m189 : public Loci::default_rule {
#line 214 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_FSIRBFr_ ; 
#line 214 "FSI_CSDVarsInput.loci"
public:
#line 214 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput036_1279672461m189() {
#line 214 "FSI_CSDVarsInput.loci"
       name_store("FSIRBFr",L_FSIRBFr_) ;
#line 214 "FSI_CSDVarsInput.loci"
       output("FSIRBFr") ;
#line 214 "FSI_CSDVarsInput.loci"
    }
#line 214 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_FSIRBFr_)= 2.; // 2
  }} ;
#line 216 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput036_1279672461m189> register_file_FSI_CSDVarsInput036_1279672461m189 ;
#line 216 "FSI_CSDVarsInput.loci"
}
#line 216 "FSI_CSDVarsInput.loci"

    
  namespace {class file_FSI_CSDVarsInput037_1279672461m189 : public Loci::default_rule {
#line 218 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_FSIRBFa_ ; 
#line 218 "FSI_CSDVarsInput.loci"
public:
#line 218 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput037_1279672461m189() {
#line 218 "FSI_CSDVarsInput.loci"
       name_store("FSIRBFa",L_FSIRBFa_) ;
#line 218 "FSI_CSDVarsInput.loci"
       output("FSIRBFa") ;
#line 218 "FSI_CSDVarsInput.loci"
    }
#line 218 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_FSIRBFa_)= 0.001; // 2
  }} ;
#line 220 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput037_1279672461m189> register_file_FSI_CSDVarsInput037_1279672461m189 ;
#line 220 "FSI_CSDVarsInput.loci"
}
#line 220 "FSI_CSDVarsInput.loci"

    
  namespace {class file_FSI_CSDVarsInput038_1279672461m189 : public Loci::default_rule {
#line 222 "FSI_CSDVarsInput.loci"
    Loci::param<int>  L_FSIRBFnr_ ; 
#line 222 "FSI_CSDVarsInput.loci"
public:
#line 222 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput038_1279672461m189() {
#line 222 "FSI_CSDVarsInput.loci"
       name_store("FSIRBFnr",L_FSIRBFnr_) ;
#line 222 "FSI_CSDVarsInput.loci"
       output("FSIRBFnr") ;
#line 222 "FSI_CSDVarsInput.loci"
    }
#line 222 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_FSIRBFnr_)= 2; // 2
  }} ;
#line 224 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput038_1279672461m189> register_file_FSI_CSDVarsInput038_1279672461m189 ;
#line 224 "FSI_CSDVarsInput.loci"
}
#line 224 "FSI_CSDVarsInput.loci"

  
    
  namespace {class file_FSI_CSDVarsInput039_1279672461m189 : public Loci::default_rule {
#line 227 "FSI_CSDVarsInput.loci"
    Loci::param<int>  L_FSIRBFMaxLinearSolverIterations_ ; 
#line 227 "FSI_CSDVarsInput.loci"
public:
#line 227 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput039_1279672461m189() {
#line 227 "FSI_CSDVarsInput.loci"
       name_store("FSIRBFMaxLinearSolverIterations",L_FSIRBFMaxLinearSolverIterations_) ;
#line 227 "FSI_CSDVarsInput.loci"
       output("FSIRBFMaxLinearSolverIterations") ;
#line 227 "FSI_CSDVarsInput.loci"
    }
#line 227 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_FSIRBFMaxLinearSolverIterations_)= 10; // 2
  }} ;
#line 229 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput039_1279672461m189> register_file_FSI_CSDVarsInput039_1279672461m189 ;
#line 229 "FSI_CSDVarsInput.loci"
}
#line 229 "FSI_CSDVarsInput.loci"

    
  namespace {class file_FSI_CSDVarsInput040_1279672461m189 : public Loci::default_rule {
#line 231 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_FSIRBFTolerance_ ; 
#line 231 "FSI_CSDVarsInput.loci"
public:
#line 231 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput040_1279672461m189() {
#line 231 "FSI_CSDVarsInput.loci"
       name_store("FSIRBFTolerance",L_FSIRBFTolerance_) ;
#line 231 "FSI_CSDVarsInput.loci"
       output("FSIRBFTolerance") ;
#line 231 "FSI_CSDVarsInput.loci"
    }
#line 231 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_FSIRBFTolerance_)= 1.e-3; // 2
  }} ;
#line 233 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput040_1279672461m189> register_file_FSI_CSDVarsInput040_1279672461m189 ;
#line 233 "FSI_CSDVarsInput.loci"
}
#line 233 "FSI_CSDVarsInput.loci"
      

  namespace {class file_FSI_CSDVarsInput041_1279672461m190 : public Loci::default_rule {
#line 235 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_FSICSD2CFDRBFr_ ; 
#line 235 "FSI_CSDVarsInput.loci"
public:
#line 235 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput041_1279672461m190() {
#line 235 "FSI_CSDVarsInput.loci"
       name_store("FSICSD2CFDRBFr",L_FSICSD2CFDRBFr_) ;
#line 235 "FSI_CSDVarsInput.loci"
       output("FSICSD2CFDRBFr") ;
#line 235 "FSI_CSDVarsInput.loci"
    }
#line 235 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_FSICSD2CFDRBFr_)= 2.; // 2
  }} ;
#line 237 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput041_1279672461m190> register_file_FSI_CSDVarsInput041_1279672461m190 ;
#line 237 "FSI_CSDVarsInput.loci"
}
#line 237 "FSI_CSDVarsInput.loci"

    
  namespace {class file_FSI_CSDVarsInput042_1279672461m190 : public Loci::default_rule {
#line 239 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_FSICSD2CFDRBFa_ ; 
#line 239 "FSI_CSDVarsInput.loci"
public:
#line 239 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput042_1279672461m190() {
#line 239 "FSI_CSDVarsInput.loci"
       name_store("FSICSD2CFDRBFa",L_FSICSD2CFDRBFa_) ;
#line 239 "FSI_CSDVarsInput.loci"
       output("FSICSD2CFDRBFa") ;
#line 239 "FSI_CSDVarsInput.loci"
    }
#line 239 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_FSICSD2CFDRBFa_)= 0.001; // 2
  }} ;
#line 241 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput042_1279672461m190> register_file_FSI_CSDVarsInput042_1279672461m190 ;
#line 241 "FSI_CSDVarsInput.loci"
}
#line 241 "FSI_CSDVarsInput.loci"

    
  namespace {class file_FSI_CSDVarsInput043_1279672461m190 : public Loci::default_rule {
#line 243 "FSI_CSDVarsInput.loci"
    Loci::param<int>  L_FSICSD2CFDRBFnr_ ; 
#line 243 "FSI_CSDVarsInput.loci"
public:
#line 243 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput043_1279672461m190() {
#line 243 "FSI_CSDVarsInput.loci"
       name_store("FSICSD2CFDRBFnr",L_FSICSD2CFDRBFnr_) ;
#line 243 "FSI_CSDVarsInput.loci"
       output("FSICSD2CFDRBFnr") ;
#line 243 "FSI_CSDVarsInput.loci"
    }
#line 243 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_FSICSD2CFDRBFnr_)= 2; // 2
  }} ;
#line 245 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput043_1279672461m190> register_file_FSI_CSDVarsInput043_1279672461m190 ;
#line 245 "FSI_CSDVarsInput.loci"
}
#line 245 "FSI_CSDVarsInput.loci"

  
    
  namespace {class file_FSI_CSDVarsInput044_1279672461m190 : public Loci::default_rule {
#line 248 "FSI_CSDVarsInput.loci"
    Loci::param<int>  L_FSICSD2CFDRBFMaxLinearSolverIterations_ ; 
#line 248 "FSI_CSDVarsInput.loci"
public:
#line 248 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput044_1279672461m190() {
#line 248 "FSI_CSDVarsInput.loci"
       name_store("FSICSD2CFDRBFMaxLinearSolverIterations",L_FSICSD2CFDRBFMaxLinearSolverIterations_) ;
#line 248 "FSI_CSDVarsInput.loci"
       output("FSICSD2CFDRBFMaxLinearSolverIterations") ;
#line 248 "FSI_CSDVarsInput.loci"
    }
#line 248 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_FSICSD2CFDRBFMaxLinearSolverIterations_)= 10; // 2
  }} ;
#line 250 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput044_1279672461m190> register_file_FSI_CSDVarsInput044_1279672461m190 ;
#line 250 "FSI_CSDVarsInput.loci"
}
#line 250 "FSI_CSDVarsInput.loci"

    
  namespace {class file_FSI_CSDVarsInput045_1279672461m190 : public Loci::default_rule {
#line 252 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_FSICSD2CFDRBFTolerance_ ; 
#line 252 "FSI_CSDVarsInput.loci"
public:
#line 252 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput045_1279672461m190() {
#line 252 "FSI_CSDVarsInput.loci"
       name_store("FSICSD2CFDRBFTolerance",L_FSICSD2CFDRBFTolerance_) ;
#line 252 "FSI_CSDVarsInput.loci"
       output("FSICSD2CFDRBFTolerance") ;
#line 252 "FSI_CSDVarsInput.loci"
    }
#line 252 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_FSICSD2CFDRBFTolerance_)= 1.e-3; // 2
  }} ;
#line 254 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput045_1279672461m190> register_file_FSI_CSDVarsInput045_1279672461m190 ;
#line 254 "FSI_CSDVarsInput.loci"
}
#line 254 "FSI_CSDVarsInput.loci"
     
  
  namespace {class file_FSI_CSDVarsInput046_1279672461m191 : public Loci::default_rule {
#line 256 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_FSIIterationTolerance_ ; 
#line 256 "FSI_CSDVarsInput.loci"
public:
#line 256 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput046_1279672461m191() {
#line 256 "FSI_CSDVarsInput.loci"
       name_store("FSIIterationTolerance",L_FSIIterationTolerance_) ;
#line 256 "FSI_CSDVarsInput.loci"
       output("FSIIterationTolerance") ;
#line 256 "FSI_CSDVarsInput.loci"
    }
#line 256 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_FSIIterationTolerance_)= 1.e-5; // 2
  }} ;
#line 258 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput046_1279672461m191> register_file_FSI_CSDVarsInput046_1279672461m191 ;
#line 258 "FSI_CSDVarsInput.loci"
}
#line 258 "FSI_CSDVarsInput.loci"
       
  
  namespace {class file_FSI_CSDVarsInput047_1279672461m191 : public Loci::default_rule {
#line 260 "FSI_CSDVarsInput.loci"
    Loci::param<int>  L_FSIIterationMinimum_ ; 
#line 260 "FSI_CSDVarsInput.loci"
public:
#line 260 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput047_1279672461m191() {
#line 260 "FSI_CSDVarsInput.loci"
       name_store("FSIIterationMinimum",L_FSIIterationMinimum_) ;
#line 260 "FSI_CSDVarsInput.loci"
       output("FSIIterationMinimum") ;
#line 260 "FSI_CSDVarsInput.loci"
    }
#line 260 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_FSIIterationMinimum_)= 5; // 2
  }} ;
#line 262 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput047_1279672461m191> register_file_FSI_CSDVarsInput047_1279672461m191 ;
#line 262 "FSI_CSDVarsInput.loci"
}
#line 262 "FSI_CSDVarsInput.loci"

  
  namespace {class file_FSI_CSDVarsInput048_1279672461m191 : public Loci::default_rule {
#line 264 "FSI_CSDVarsInput.loci"
    Loci::param<int>  L_maxIterationsPerFSI_ ; 
#line 264 "FSI_CSDVarsInput.loci"
public:
#line 264 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput048_1279672461m191() {
#line 264 "FSI_CSDVarsInput.loci"
       name_store("maxIterationsPerFSI",L_maxIterationsPerFSI_) ;
#line 264 "FSI_CSDVarsInput.loci"
       output("maxIterationsPerFSI") ;
#line 264 "FSI_CSDVarsInput.loci"
    }
#line 264 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_maxIterationsPerFSI_)= 5; // 2
  	if ((*L_maxIterationsPerFSI_)<=1) {
	  if (Loci::MPI_rank==0) cerr << "maxIterationsPerFSI should be larger than 1" << endl ;
	}
  }} ;
#line 269 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput048_1279672461m191> register_file_FSI_CSDVarsInput048_1279672461m191 ;
#line 269 "FSI_CSDVarsInput.loci"
}
#line 269 "FSI_CSDVarsInput.loci"
  
    
  namespace {class file_FSI_CSDVarsInput049_1279672461m191 : public Loci::optional_rule {
#line 271 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDEulerXstart_ ; 
#line 271 "FSI_CSDVarsInput.loci"
public:
#line 271 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput049_1279672461m191() {
#line 271 "FSI_CSDVarsInput.loci"
       name_store("CSDEulerXstart",L_CSDEulerXstart_) ;
#line 271 "FSI_CSDVarsInput.loci"
       output("CSDEulerXstart") ;
#line 271 "FSI_CSDVarsInput.loci"
    }
#line 271 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDEulerXstart_)= 0; // 2
  }} ;
#line 273 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput049_1279672461m191> register_file_FSI_CSDVarsInput049_1279672461m191 ;
#line 273 "FSI_CSDVarsInput.loci"
}
#line 273 "FSI_CSDVarsInput.loci"
           

  namespace {class file_FSI_CSDVarsInput050_1279672461m191 : public Loci::optional_rule {
#line 275 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDEulerXend_ ; 
#line 275 "FSI_CSDVarsInput.loci"
public:
#line 275 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput050_1279672461m191() {
#line 275 "FSI_CSDVarsInput.loci"
       name_store("CSDEulerXend",L_CSDEulerXend_) ;
#line 275 "FSI_CSDVarsInput.loci"
       output("CSDEulerXend") ;
#line 275 "FSI_CSDVarsInput.loci"
    }
#line 275 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDEulerXend_)= 1; // 2
  }} ;
#line 277 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput050_1279672461m191> register_file_FSI_CSDVarsInput050_1279672461m191 ;
#line 277 "FSI_CSDVarsInput.loci"
}
#line 277 "FSI_CSDVarsInput.loci"
           
 
  namespace {class file_FSI_CSDVarsInput051_1279672461m192 : public Loci::optional_rule {
#line 279 "FSI_CSDVarsInput.loci"
    Loci::param<int>  L_CSDEulerXnum_ ; 
#line 279 "FSI_CSDVarsInput.loci"
public:
#line 279 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput051_1279672461m192() {
#line 279 "FSI_CSDVarsInput.loci"
       name_store("CSDEulerXnum",L_CSDEulerXnum_) ;
#line 279 "FSI_CSDVarsInput.loci"
       output("CSDEulerXnum") ;
#line 279 "FSI_CSDVarsInput.loci"
    }
#line 279 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDEulerXnum_)= 10; // 2
  }} ;
#line 281 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput051_1279672461m192> register_file_FSI_CSDVarsInput051_1279672461m192 ;
#line 281 "FSI_CSDVarsInput.loci"
}
#line 281 "FSI_CSDVarsInput.loci"
                 
  
  namespace {class file_FSI_CSDVarsInput052_1279672461m192 : public Loci::optional_rule {
#line 283 "FSI_CSDVarsInput.loci"
    Loci::param<real>  L_CSDEulerChord_ ; 
#line 283 "FSI_CSDVarsInput.loci"
public:
#line 283 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput052_1279672461m192() {
#line 283 "FSI_CSDVarsInput.loci"
       name_store("CSDEulerChord",L_CSDEulerChord_) ;
#line 283 "FSI_CSDVarsInput.loci"
       output("CSDEulerChord") ;
#line 283 "FSI_CSDVarsInput.loci"
    }
#line 283 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDEulerChord_)= 1.; // 2
  }} ;
#line 285 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput052_1279672461m192> register_file_FSI_CSDVarsInput052_1279672461m192 ;
#line 285 "FSI_CSDVarsInput.loci"
}
#line 285 "FSI_CSDVarsInput.loci"
   

  namespace {class file_FSI_CSDVarsInput053_1279672461m192 : public Loci::optional_rule {
#line 287 "FSI_CSDVarsInput.loci"
    Loci::param<int>  L_CSDEulerAxis_ ; 
#line 287 "FSI_CSDVarsInput.loci"
public:
#line 287 "FSI_CSDVarsInput.loci"
    file_FSI_CSDVarsInput053_1279672461m192() {
#line 287 "FSI_CSDVarsInput.loci"
       name_store("CSDEulerAxis",L_CSDEulerAxis_) ;
#line 287 "FSI_CSDVarsInput.loci"
       output("CSDEulerAxis") ;
#line 287 "FSI_CSDVarsInput.loci"
    }
#line 287 "FSI_CSDVarsInput.loci"
    void compute(const Loci::sequence &seq) { 
  	(*L_CSDEulerAxis_)= 0; // // beam direction: 0->x, 1->y, 2->z
  	if (Loci::MPI_rank==0) {
	  cout << "CSDEulerAxis=" ;
	  if ((*L_CSDEulerAxis_)==0) cout << "x" << endl ;
	  if ((*L_CSDEulerAxis_)==1) cout << "y" << endl ;
	  if ((*L_CSDEulerAxis_)==2) cout << "z" << endl ;
	}
  }} ;
#line 295 "FSI_CSDVarsInput.loci"
Loci::register_rule<file_FSI_CSDVarsInput053_1279672461m192> register_file_FSI_CSDVarsInput053_1279672461m192 ;
#line 295 "FSI_CSDVarsInput.loci"
}
#line 295 "FSI_CSDVarsInput.loci"

}
   
