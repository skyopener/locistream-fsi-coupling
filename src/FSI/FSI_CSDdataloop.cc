#line 1 "FSI_CSDdataloop.loci"
//-----------------------------------------------------------------------------
// Description: This file contains some of the basic rules common to all
//   grid movement schemes.
//-----------------------------------------------------------------------------

// Standard library includes.
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib> // for exit function
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
using std::vector ;

using namespace std ;

// Loci includes.
#include <Loci.h>

// StreamUns includes.
#include "sciTypes.h"

#define BOOST_DISABLE_ASSERTS // assertion disabled                                      
// boost::ublas includes
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace ublas = boost::numeric::ublas ;

// Including the C-Fortran template for the subroutine call to NLAMS	
	extern "C"{
        void pass_(const int *, const int*, const int*, const double*, const double*, const double*,
		   const int*, const int*,
                   const int*, const int*, const double*,
                   const double*,  const double*,  const double*,  const double*,
                   const int*,  const double*,  const double*,  const double*, 
                   const int*, const double*, const double*, const double*,
                   const double*, const double*, const double*, 
                   const double*, const int*,  const int*,
                   const int*, const int*,
                   const double*, const double*, const double*, const double*,
                   const double*, const double*, const double*, 
                   double*,double*,double*,double*,double*);
   
      	void beam1d_(const int *, const int*, const int*, const double*, const int*, const double*, const double*, const double*, const double*, const double*, const double*,
        const double*, const double*, const double*, const double*, const double*, const double*, const double*, const double*, const int*, 
        double*, double*, double*, double*, double*) ;
        
	}

// Forward declarations using C linkage to avoid mangling the names of
// the Fortran functions that will be called.
// Idea I: Store the data in LOCI which will be required in CSD solver
// Idea II: Excute CSD solver at each step and kill it after get solution.
// 			Only at t=0: readin the mesh, connectivity, bcs, etc.
// Idea III: for constant time, i.e. n=constant, but for advancing inner iteration, i.e. it=changing, we want to keep the unsteady term (CSD) constant.

// CSDvariables
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
// $type CSDForcePreEuler_ic blackbox<ublas::vector<real> > 
// $type CSDForcePreEuler blackbox<ublas::vector<real> > 
// $type CSDForcePreEulerStar blackbox<ublas::vector<real> > 

#line 63 "FSI_CSDdataloop.loci"


namespace streamUns {
            
//-----------------------------------------------------------------------------
// Rules to read in the grids

// -- Read in the nodes, connectivity, bc files for CSD------------------------------------------------------------------------------------- 
// Read in the initial CSD nodes: undeformed

// $type CFDIterationFinished param<bool> 
// $type FSIIterationFinished param<bool> 


namespace {class file_FSI_CSDdataloop000_1280807839m47 : public Loci::blackbox_rule {
#line 78 "FSI_CSDdataloop.loci"
    Loci::const_param<string>  L_CSDMeshFilename_ ; 
#line 78 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 78 "FSI_CSDdataloop.loci"
public:
#line 78 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop000_1280807839m47() {
#line 78 "FSI_CSDdataloop.loci"
       name_store("CSDMeshFilename",L_CSDMeshFilename_) ;
#line 78 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 78 "FSI_CSDdataloop.loci"
       input("CSDMeshFilename") ;
#line 78 "FSI_CSDdataloop.loci"
       output("CSDnodes_ic") ;
#line 78 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 78 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 78 "FSI_CSDdataloop.loci"
    }
#line 78 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
 	string filename = *L_CSDMeshFilename_;
 	int Nb = 0;
	if (Loci::MPI_rank==0) cout << "[I] Reading the CSD mesh" << endl;
	ifstream meshData ;
	meshData.open(filename.c_str(), ios::in) ;
	if (!meshData.is_open()) {
		if (Loci::MPI_rank==0) cerr << "[E] CSD meshfile file " << Nb  << " couldn't be opened" << endl;
		exit(1);
	}
	cout << "[I] ..." << filename << " is successfully opened " << endl;
	
	// initialization
	int counter = 0;
	int tempInt ;
//	vector<double> tempVector(3) ;

	string dumpLine ;
// Count the number of elements
	while (getline(meshData, dumpLine)) 
	{
		counter++;
	}
	
	Nb = counter ;
	
// Resize the CSD nodes
	(*L_CSDnodes_ic_).resize(Nb,3) ;
	
// Clear the stream and go to the first line
	meshData.clear() ;
	meshData.seekg(0) ;
	
// Create nodes by reading the mesh until EOF found
	for (int i=0; i<Nb; ++i) {
		meshData >> tempInt >> (*L_CSDnodes_ic_)(i,2) >> (*L_CSDnodes_ic_)(i,0)>> (*L_CSDnodes_ic_)(i,1); // data CSD coordinates
	}
	
	if (Loci::MPI_rank==0) cout << "[I] ...Number of CSD nodes: " << Nb << endl;

	meshData.close();
}    void compute(const Loci::sequence &seq) { 
#line 122 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 122 "FSI_CSDdataloop.loci"
    }
#line 122 "FSI_CSDdataloop.loci"
} ;
#line 122 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop000_1280807839m47> register_file_FSI_CSDdataloop000_1280807839m47 ;
#line 122 "FSI_CSDdataloop.loci"
}
#line 122 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop001_1280807839m48 : public Loci::blackbox_rule {
#line 123 "FSI_CSDdataloop.loci"
    Loci::const_param<string>  L_CSDConnectivityFilename_ ; 
#line 123 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<int,ublas::column_major> >  L_CSDConnectivity_ ; 
#line 123 "FSI_CSDdataloop.loci"
public:
#line 123 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop001_1280807839m48() {
#line 123 "FSI_CSDdataloop.loci"
       name_store("CSDConnectivityFilename",L_CSDConnectivityFilename_) ;
#line 123 "FSI_CSDdataloop.loci"
       name_store("CSDConnectivity",L_CSDConnectivity_) ;
#line 123 "FSI_CSDdataloop.loci"
       input("CSDConnectivityFilename") ;
#line 123 "FSI_CSDdataloop.loci"
       output("CSDConnectivity") ;
#line 123 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 123 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 123 "FSI_CSDdataloop.loci"
    }
#line 123 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	string filename = *L_CSDConnectivityFilename_;
	int Ne=0;
	
	//Open file 
	if (Loci::MPI_rank==0) cout << "[I] Reading the CSD connectivity" << endl;
	ifstream connectData ;
	connectData.open(filename.c_str(), ios::in) ;
	if (!connectData.is_open()) {
		if (Loci::MPI_rank==0) cerr << "[E] CSD connectivity file " << filename  << " couldn't be opened" << endl;
		exit(1);
	}
	cout << "[I] ..." << filename << " is successfully opened " << endl;
	
	// initialization
	int counter = 0;
	int tempInt ;
//	vector<int> tempVector(3) ;

	string dumpLine ;
// Count the number of elements
	while (getline(connectData, dumpLine)) counter++ ;
	
	Ne = counter ;
	
// Resize the CSD nodes
	(*L_CSDConnectivity_).resize(Ne,3) ;
	
// Clear the stream and go to the first line
	connectData.clear() ;
	connectData.seekg(0) ;
	
// Create nodes by reading the mesh until EOF found
	for (int i=0; i<Ne; ++i) {
		connectData >> tempInt >> tempInt >> (*L_CSDConnectivity_)(i,0) >> (*L_CSDConnectivity_)(i,1) >> (*L_CSDConnectivity_)(i,2) ;
	}
		
	if (Loci::MPI_rank==0) cout << "[I] ...Number of CSD elements: " << Ne << endl;

	connectData.close();
}    void compute(const Loci::sequence &seq) { 
#line 167 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 167 "FSI_CSDdataloop.loci"
    }
#line 167 "FSI_CSDdataloop.loci"
} ;
#line 167 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop001_1280807839m48> register_file_FSI_CSDdataloop001_1280807839m48 ;
#line 167 "FSI_CSDdataloop.loci"
}
#line 167 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop002_1280807839m49 : public Loci::singleton_rule {
#line 167 "FSI_CSDdataloop.loci"
    Loci::const_param<string>  L_CSDBCFilename_ ; 
#line 167 "FSI_CSDdataloop.loci"
    Loci::param<vector<int> >  L_CSDBCdof_ ; 
#line 167 "FSI_CSDdataloop.loci"
    Loci::param<vector<real> >  L_CSDBCZeroConstraint_ ; 
#line 167 "FSI_CSDdataloop.loci"
public:
#line 167 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop002_1280807839m49() {
#line 167 "FSI_CSDdataloop.loci"
       name_store("CSDBCFilename",L_CSDBCFilename_) ;
#line 167 "FSI_CSDdataloop.loci"
       name_store("CSDBCdof",L_CSDBCdof_) ;
#line 167 "FSI_CSDdataloop.loci"
       name_store("CSDBCZeroConstraint",L_CSDBCZeroConstraint_) ;
#line 167 "FSI_CSDdataloop.loci"
       input("CSDBCFilename") ;
#line 167 "FSI_CSDdataloop.loci"
       output("CSDBCdof") ;
#line 167 "FSI_CSDdataloop.loci"
       output("CSDBCZeroConstraint") ;
#line 167 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 167 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 167 "FSI_CSDdataloop.loci"
    }
#line 167 "FSI_CSDdataloop.loci"
    void compute(const Loci::sequence &seq) {  

	string filename = (*L_CSDBCFilename_);
	int Nbc=0;
	//Open file 
	if (Loci::MPI_rank==0) cout << "[I] Reading the CSD boundary conditions" << endl;
	ifstream bcData ;
	bcData.open(filename.c_str(), ios::in) ;
	if (!bcData.is_open()) {
		if (Loci::MPI_rank==0) cerr << "[E] CSD connectivity file " << filename  << " couldn't be opened" << endl;
		exit(1);
	}
	cout << "[I] ..." << filename << " is successfully opened " << endl;
	
	// initialization
	int counter = 0 ;
	int tempInt ;
	real tempReal ;

	string dumpLine ;
// Count the number of elements
	while (getline(bcData, dumpLine)) 
	{
		counter++;
	}
	
	Nbc = counter ;
	
// Resize the CSD nodes
	((*L_CSDBCdof_)).resize(Nbc) ;
	((*L_CSDBCZeroConstraint_)).resize(Nbc) ;
	
// Clear the stream and go to the first line
	bcData.clear() ;
	bcData.seekg(0) ;
	
// Create nodes by reading the mesh until EOF found
	for (int i=0; i<Nbc; ++i) {
	bcData >> tempInt >> tempReal ;
	((*L_CSDBCdof_))[i] = tempInt ;
	((*L_CSDBCZeroConstraint_))[i] = tempReal ;
	}
		
	if (Loci::MPI_rank==0) cout << "[I] ...Number of CSD bcs: " << Nbc<< endl;

	bcData.close();
}} ;
#line 213 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop002_1280807839m49> register_file_FSI_CSDdataloop002_1280807839m49 ;
#line 213 "FSI_CSDdataloop.loci"
}
#line 213 "FSI_CSDdataloop.loci"


// -- Time advancing for the CSDnodes ----------------------------------------------------------------------------------------------------
// Time Build rule for the CSD nodes
namespace {class file_FSI_CSDdataloop003_1280807839m49 : public Loci::blackbox_rule {
#line 218 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 218 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_n_EQ_0__ ; 
#line 218 "FSI_CSDdataloop.loci"
public:
#line 218 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop003_1280807839m49() {
#line 218 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 218 "FSI_CSDdataloop.loci"
       name_store("CSDnodes{n=0}",L_CSDnodes_n_EQ_0__) ;
#line 218 "FSI_CSDdataloop.loci"
       input("CSDnodes_ic") ;
#line 218 "FSI_CSDdataloop.loci"
       output("CSDnodes{n=0}") ;
#line 218 "FSI_CSDdataloop.loci"
       constraint("FSICoupling") ;
#line 218 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 218 "FSI_CSDdataloop.loci"
    }
#line 218 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodes_n_EQ_0__).resize((*L_CSDnodes_ic_).size1(),3) ; std::fill( (*L_CSDnodes_n_EQ_0__).data().begin(), (*L_CSDnodes_n_EQ_0__).data().end(), 0. );  
	(*L_CSDnodes_n_EQ_0__) = (*L_CSDnodes_ic_) ;
}    void compute(const Loci::sequence &seq) { 
#line 224 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 224 "FSI_CSDdataloop.loci"
    }
#line 224 "FSI_CSDdataloop.loci"
} ;
#line 224 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop003_1280807839m49> register_file_FSI_CSDdataloop003_1280807839m49 ;
#line 224 "FSI_CSDdataloop.loci"
}
#line 224 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop004_1280807839m50 : public Loci::blackbox_rule {
#line 225 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_n__ ; 
#line 225 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_it_nit_EQ_0__ ; 
#line 225 "FSI_CSDdataloop.loci"
public:
#line 225 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop004_1280807839m50() {
#line 225 "FSI_CSDdataloop.loci"
       name_store("CSDnodes{n}",L_CSDnodes_n__) ;
#line 225 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_it{n,it=0}",L_CSDnodes_it_nit_EQ_0__) ;
#line 225 "FSI_CSDdataloop.loci"
       input("CSDnodes{n}") ;
#line 225 "FSI_CSDdataloop.loci"
       output("CSDnodes_it{n,it=0}") ;
#line 225 "FSI_CSDdataloop.loci"
       constraint("FSICoupling") ;
#line 225 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 225 "FSI_CSDdataloop.loci"
    }
#line 225 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodes_it_nit_EQ_0__).resize((*L_CSDnodes_n__).size1(),3) ;std::fill( (*L_CSDnodes_it_nit_EQ_0__).data().begin(), (*L_CSDnodes_it_nit_EQ_0__).data().end(), 0. );  
	(*L_CSDnodes_it_nit_EQ_0__) = (*L_CSDnodes_n__) ;
//	if (Loci::MPI_rank==0) cout << "Inside CSDnodes_it{n,it=0}" << endl;
}    void compute(const Loci::sequence &seq) { 
#line 232 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 232 "FSI_CSDdataloop.loci"
    }
#line 232 "FSI_CSDdataloop.loci"
} ;
#line 232 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop004_1280807839m50> register_file_FSI_CSDdataloop004_1280807839m50 ;
#line 232 "FSI_CSDdataloop.loci"
}
#line 232 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop005_1280807839m51 : public Loci::blackbox_rule {
#line 233 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 233 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_it_nit__ ; 
#line 233 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDdisplacementsStar_nit__ ; 
#line 233 "FSI_CSDdataloop.loci"
    Loci::const_param<bool>  L_CFDIterationFinished_nit_M_1__ ; 
#line 233 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_it_nit_P_1__ ; 
#line 233 "FSI_CSDdataloop.loci"
public:
#line 233 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop005_1280807839m51() {
#line 233 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 233 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_it{n,it}",L_CSDnodes_it_nit__) ;
#line 233 "FSI_CSDdataloop.loci"
       name_store("CSDdisplacementsStar{n,it}",L_CSDdisplacementsStar_nit__) ;
#line 233 "FSI_CSDdataloop.loci"
       name_store("CFDIterationFinished{n,it-1}",L_CFDIterationFinished_nit_M_1__) ;
#line 233 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_it{n,it+1}",L_CSDnodes_it_nit_P_1__) ;
#line 233 "FSI_CSDdataloop.loci"
       input("CSDnodes_it{n,it},CSDnodes_ic,CSDdisplacementsStar{n,it},CFDIterationFinished{n,it-1}") ;
#line 233 "FSI_CSDdataloop.loci"
       output("CSDnodes_it{n,it+1}") ;
#line 233 "FSI_CSDdataloop.loci"
       constraint("FSICoupling") ;
#line 233 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 233 "FSI_CSDdataloop.loci"
    }
#line 233 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	
	(*L_CSDnodes_it_nit_P_1__).resize((*L_CSDnodes_ic_).size1(),3) ; std::fill( (*L_CSDnodes_it_nit_P_1__).data().begin(), (*L_CSDnodes_it_nit_P_1__).data().end(), 0. ); 
	if (*L_CFDIterationFinished_nit_M_1__) {	
	//	if (Loci::MPI_rank==0) cout << "Inside CSDnodes_it{n,it+1} priority in" << endl ;
	
//	for(int i=0;i<(*$CSDnodes_it{n,it+1}).size1();++i) {
//		(*$CSDnodes_it{n,it+1})(i,0) = (*$CSDnodes_ic)(i,0) + (*$CSDdisplacementsStar{n,it})(i,2) ; // CSD.x = CFD.z
//		(*$CSDnodes_it{n,it+1})(i,1) = (*$CSDnodes_ic)(i,1) + (*$CSDdisplacementsStar{n,it})(i,0) ; // CSD.y = CFD.x
//		(*$CSDnodes_it{n,it+1})(i,2) = (*$CSDnodes_ic)(i,2) + (*$CSDdisplacementsStar{n,it})(i,1) ; // CSD.z = CFD.y
//	}
		
	(*L_CSDnodes_it_nit_P_1__) = (*L_CSDnodes_ic_) + (*L_CSDdisplacementsStar_nit__); // CHECK!!!!
	//	if (Loci::MPI_rank==0) cout << "Inside CSDnodes_it{n,it+1} priority out" << endl ;
	} else {
		(*L_CSDnodes_it_nit_P_1__) = (*L_CSDnodes_it_nit__); // CHECK!!!!
	}
}    void compute(const Loci::sequence &seq) { 
#line 253 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 253 "FSI_CSDdataloop.loci"
    }
#line 253 "FSI_CSDdataloop.loci"
} ;
#line 253 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop005_1280807839m51> register_file_FSI_CSDdataloop005_1280807839m51 ;
#line 253 "FSI_CSDdataloop.loci"
}
#line 253 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop006_1280807839m53 : public Loci::blackbox_rule {
#line 254 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_it_nit__ ; 
#line 254 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_n_P_1__ ; 
#line 254 "FSI_CSDdataloop.loci"
public:
#line 254 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop006_1280807839m53() {
#line 254 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_it{n,it}",L_CSDnodes_it_nit__) ;
#line 254 "FSI_CSDdataloop.loci"
       name_store("CSDnodes{n+1}",L_CSDnodes_n_P_1__) ;
#line 254 "FSI_CSDdataloop.loci"
       input("CSDnodes_it{n,it}") ;
#line 254 "FSI_CSDdataloop.loci"
       output("CSDnodes{n+1}") ;
#line 254 "FSI_CSDdataloop.loci"
       constraint("FSICoupling") ;
#line 254 "FSI_CSDdataloop.loci"
       conditional("iterationFinished{n,it-1}") ;
#line 254 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 254 "FSI_CSDdataloop.loci"
    }
#line 254 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	//	if (Loci::MPI_rank==0) cout << "Inside CSDnodes{n+1} normal in" << endl ;
	(*L_CSDnodes_n_P_1__).resize((*L_CSDnodes_it_nit__).size1(),3) ; std::fill( (*L_CSDnodes_n_P_1__).data().begin(),(*L_CSDnodes_n_P_1__).data().end(), 0. ); 
	(*L_CSDnodes_n_P_1__) = (*L_CSDnodes_it_nit__)  ;
//	if (Loci::MPI_rank==0) cout << "Inside CSDnodes{n+1} normal out" << endl ;
}    void compute(const Loci::sequence &seq) { 
#line 272 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 272 "FSI_CSDdataloop.loci"
    }
#line 272 "FSI_CSDdataloop.loci"
} ;
#line 272 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop006_1280807839m53> register_file_FSI_CSDdataloop006_1280807839m53 ;
#line 272 "FSI_CSDdataloop.loci"
}
#line 272 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop007_1280807839m54 : public Loci::blackbox_rule {
#line 272 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 272 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDisp_ic_ ; 
#line 272 "FSI_CSDdataloop.loci"
public:
#line 272 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop007_1280807839m54() {
#line 272 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 272 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDisp_ic",L_CSDnodesDisp_ic_) ;
#line 272 "FSI_CSDdataloop.loci"
       input("CSDnodes_ic") ;
#line 272 "FSI_CSDdataloop.loci"
       output("CSDnodesDisp_ic") ;
#line 272 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 272 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 272 "FSI_CSDdataloop.loci"
    }
#line 272 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesDisp_ic_).resize(6*(*L_CSDnodes_ic_).size1(),1) ;
	std::fill( (*L_CSDnodesDisp_ic_).data().begin(),(*L_CSDnodesDisp_ic_).data().end(), 0. ); 
}    void compute(const Loci::sequence &seq) { 
#line 278 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 278 "FSI_CSDdataloop.loci"
    }
#line 278 "FSI_CSDdataloop.loci"
} ;
#line 278 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop007_1280807839m54> register_file_FSI_CSDdataloop007_1280807839m54 ;
#line 278 "FSI_CSDdataloop.loci"
}
#line 278 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop008_1280807839m54 : public Loci::blackbox_rule {
#line 278 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDisp_ic_ ; 
#line 278 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDisp_n_EQ_0__ ; 
#line 278 "FSI_CSDdataloop.loci"
public:
#line 278 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop008_1280807839m54() {
#line 278 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDisp_ic",L_CSDnodesDisp_ic_) ;
#line 278 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDisp{n=0}",L_CSDnodesDisp_n_EQ_0__) ;
#line 278 "FSI_CSDdataloop.loci"
       input("CSDnodesDisp_ic") ;
#line 278 "FSI_CSDdataloop.loci"
       output("CSDnodesDisp{n=0}") ;
#line 278 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 278 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 278 "FSI_CSDdataloop.loci"
    }
#line 278 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesDisp_n_EQ_0__).resize((*L_CSDnodesDisp_ic_).size1(),1) ; std::fill( (*L_CSDnodesDisp_n_EQ_0__).data().begin(),(*L_CSDnodesDisp_n_EQ_0__).data().end(), 0. ); 
	(*L_CSDnodesDisp_n_EQ_0__) = (*L_CSDnodesDisp_ic_) ;
}    void compute(const Loci::sequence &seq) { 
#line 284 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 284 "FSI_CSDdataloop.loci"
    }
#line 284 "FSI_CSDdataloop.loci"
} ;
#line 284 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop008_1280807839m54> register_file_FSI_CSDdataloop008_1280807839m54 ;
#line 284 "FSI_CSDdataloop.loci"
}
#line 284 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop009_1280807839m55 : public Loci::blackbox_rule {
#line 284 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDisp_n__ ; 
#line 284 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDisp_nit_EQ_0__ ; 
#line 284 "FSI_CSDdataloop.loci"
public:
#line 284 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop009_1280807839m55() {
#line 284 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDisp{n}",L_CSDnodesDisp_n__) ;
#line 284 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDisp{n,it=0}",L_CSDnodesDisp_nit_EQ_0__) ;
#line 284 "FSI_CSDdataloop.loci"
       input("CSDnodesDisp{n}") ;
#line 284 "FSI_CSDdataloop.loci"
       output("CSDnodesDisp{n,it=0}") ;
#line 284 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 284 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 284 "FSI_CSDdataloop.loci"
    }
#line 284 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
//	if (Loci::MPI_rank==0) cout << "Inside CSDnodesDisp{n,it=0}" << endl;
	(*L_CSDnodesDisp_nit_EQ_0__).resize((*L_CSDnodesDisp_n__).size1(),1) ;	std::fill( (*L_CSDnodesDisp_nit_EQ_0__).data().begin(),(*L_CSDnodesDisp_nit_EQ_0__).data().end(), 0. ); 
	(*L_CSDnodesDisp_nit_EQ_0__) = (*L_CSDnodesDisp_n__) ;
//	if (Loci::MPI_rank==0) cout << "Inside CSDnodesDisp{n,it=0}" << endl;
}    void compute(const Loci::sequence &seq) { 
#line 292 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 292 "FSI_CSDdataloop.loci"
    }
#line 292 "FSI_CSDdataloop.loci"
} ;
#line 292 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop009_1280807839m55> register_file_FSI_CSDdataloop009_1280807839m55 ;
#line 292 "FSI_CSDdataloop.loci"
}
#line 292 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop010_1280807839m56 : public Loci::blackbox_rule {
#line 292 "FSI_CSDdataloop.loci"
    Loci::const_param<bool>  L_CFDIterationFinished_nit_M_1__ ; 
#line 292 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDisp_nit__ ; 
#line 292 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDispStar_nit__ ; 
#line 292 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDisp_nit_P_1__ ; 
#line 292 "FSI_CSDdataloop.loci"
public:
#line 292 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop010_1280807839m56() {
#line 292 "FSI_CSDdataloop.loci"
       name_store("CFDIterationFinished{n,it-1}",L_CFDIterationFinished_nit_M_1__) ;
#line 292 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDisp{n,it}",L_CSDnodesDisp_nit__) ;
#line 292 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDispStar{n,it}",L_CSDnodesDispStar_nit__) ;
#line 292 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDisp{n,it+1}",L_CSDnodesDisp_nit_P_1__) ;
#line 292 "FSI_CSDdataloop.loci"
       input("CSDnodesDisp{n,it},CSDnodesDispStar{n,it},CFDIterationFinished{n,it-1}") ;
#line 292 "FSI_CSDdataloop.loci"
       output("CSDnodesDisp{n,it+1}") ;
#line 292 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 292 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 292 "FSI_CSDdataloop.loci"
    }
#line 292 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesDisp_nit_P_1__).resize((*L_CSDnodesDispStar_nit__).size1(),1) ;	std::fill( (*L_CSDnodesDisp_nit_P_1__).data().begin(),(*L_CSDnodesDisp_nit_P_1__).data().end(), 0. );
	if (*L_CFDIterationFinished_nit_M_1__) {
		(*L_CSDnodesDisp_nit_P_1__) = (*L_CSDnodesDispStar_nit__) ;
	} else {
		(*L_CSDnodesDisp_nit_P_1__) = (*L_CSDnodesDisp_nit__) ;
	}
}    void compute(const Loci::sequence &seq) { 
#line 307 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 307 "FSI_CSDdataloop.loci"
    }
#line 307 "FSI_CSDdataloop.loci"
} ;
#line 307 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop010_1280807839m56> register_file_FSI_CSDdataloop010_1280807839m56 ;
#line 307 "FSI_CSDdataloop.loci"
}
#line 307 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop011_1280807839m57 : public Loci::blackbox_rule {
#line 307 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDisp_nit__ ; 
#line 307 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDisp_n_P_1__ ; 
#line 307 "FSI_CSDdataloop.loci"
public:
#line 307 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop011_1280807839m57() {
#line 307 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDisp{n,it}",L_CSDnodesDisp_nit__) ;
#line 307 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDisp{n+1}",L_CSDnodesDisp_n_P_1__) ;
#line 307 "FSI_CSDdataloop.loci"
       input("CSDnodesDisp{n,it}") ;
#line 307 "FSI_CSDdataloop.loci"
       output("CSDnodesDisp{n+1}") ;
#line 307 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 307 "FSI_CSDdataloop.loci"
       conditional("iterationFinished{n,it-1}") ;
#line 307 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 307 "FSI_CSDdataloop.loci"
    }
#line 307 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesDisp_n_P_1__).resize((*L_CSDnodesDisp_nit__).size1(),1) ; std::fill( (*L_CSDnodesDisp_n_P_1__).data().begin(),(*L_CSDnodesDisp_n_P_1__).data().end(), 0. );
	(*L_CSDnodesDisp_n_P_1__) = (*L_CSDnodesDisp_nit__) ;
}    void compute(const Loci::sequence &seq) { 
#line 320 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 320 "FSI_CSDdataloop.loci"
    }
#line 320 "FSI_CSDdataloop.loci"
} ;
#line 320 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop011_1280807839m57> register_file_FSI_CSDdataloop011_1280807839m57 ;
#line 320 "FSI_CSDdataloop.loci"
}
#line 320 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop012_1280807839m58 : public Loci::blackbox_rule {
#line 320 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 320 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVel_ic_ ; 
#line 320 "FSI_CSDdataloop.loci"
public:
#line 320 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop012_1280807839m58() {
#line 320 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 320 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVel_ic",L_CSDnodesVel_ic_) ;
#line 320 "FSI_CSDdataloop.loci"
       input("CSDnodes_ic") ;
#line 320 "FSI_CSDdataloop.loci"
       output("CSDnodesVel_ic") ;
#line 320 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 320 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 320 "FSI_CSDdataloop.loci"
    }
#line 320 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesVel_ic_).resize(6*(*L_CSDnodes_ic_).size1(),1) ; std::fill( (*L_CSDnodesVel_ic_).data().begin(),(*L_CSDnodesVel_ic_).data().end(), 0. );
	(*L_CSDnodesVel_ic_).clear() ; // Set all elements to zero
}    void compute(const Loci::sequence &seq) { 
#line 326 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 326 "FSI_CSDdataloop.loci"
    }
#line 326 "FSI_CSDdataloop.loci"
} ;
#line 326 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop012_1280807839m58> register_file_FSI_CSDdataloop012_1280807839m58 ;
#line 326 "FSI_CSDdataloop.loci"
}
#line 326 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop013_1280807839m58 : public Loci::blackbox_rule {
#line 326 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVel_ic_ ; 
#line 326 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVel_n_EQ_0__ ; 
#line 326 "FSI_CSDdataloop.loci"
public:
#line 326 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop013_1280807839m58() {
#line 326 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVel_ic",L_CSDnodesVel_ic_) ;
#line 326 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVel{n=0}",L_CSDnodesVel_n_EQ_0__) ;
#line 326 "FSI_CSDdataloop.loci"
       input("CSDnodesVel_ic") ;
#line 326 "FSI_CSDdataloop.loci"
       output("CSDnodesVel{n=0}") ;
#line 326 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 326 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 326 "FSI_CSDdataloop.loci"
    }
#line 326 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesVel_n_EQ_0__).resize((*L_CSDnodesVel_ic_).size1(),1) ; std::fill( (*L_CSDnodesVel_n_EQ_0__).data().begin(),(*L_CSDnodesVel_n_EQ_0__).data().end(), 0. );
	(*L_CSDnodesVel_n_EQ_0__) = (*L_CSDnodesVel_ic_) ;
}    void compute(const Loci::sequence &seq) { 
#line 333 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 333 "FSI_CSDdataloop.loci"
    }
#line 333 "FSI_CSDdataloop.loci"
} ;
#line 333 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop013_1280807839m58> register_file_FSI_CSDdataloop013_1280807839m58 ;
#line 333 "FSI_CSDdataloop.loci"
}
#line 333 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop014_1280807839m59 : public Loci::blackbox_rule {
#line 333 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVel_n__ ; 
#line 333 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVel_nit_EQ_0__ ; 
#line 333 "FSI_CSDdataloop.loci"
public:
#line 333 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop014_1280807839m59() {
#line 333 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVel{n}",L_CSDnodesVel_n__) ;
#line 333 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVel{n,it=0}",L_CSDnodesVel_nit_EQ_0__) ;
#line 333 "FSI_CSDdataloop.loci"
       input("CSDnodesVel{n}") ;
#line 333 "FSI_CSDdataloop.loci"
       output("CSDnodesVel{n,it=0}") ;
#line 333 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 333 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 333 "FSI_CSDdataloop.loci"
    }
#line 333 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesVel_nit_EQ_0__).resize((*L_CSDnodesVel_n__).size1(),1) ; std::fill( (*L_CSDnodesVel_nit_EQ_0__).data().begin(),(*L_CSDnodesVel_nit_EQ_0__).data().end(), 0. );
	(*L_CSDnodesVel_nit_EQ_0__) = *L_CSDnodesVel_n__ ;
	//if (Loci::MPI_rank==0) cout << "Inside CSDnodesVel{n,it=0}" << endl;
}    void compute(const Loci::sequence &seq) { 
#line 340 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 340 "FSI_CSDdataloop.loci"
    }
#line 340 "FSI_CSDdataloop.loci"
} ;
#line 340 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop014_1280807839m59> register_file_FSI_CSDdataloop014_1280807839m59 ;
#line 340 "FSI_CSDdataloop.loci"
}
#line 340 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop015_1280807839m60 : public Loci::blackbox_rule {
#line 340 "FSI_CSDdataloop.loci"
    Loci::const_param<bool>  L_CFDIterationFinished_nit_M_1__ ; 
#line 340 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVel_nit__ ; 
#line 340 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVelStar_nit__ ; 
#line 340 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVel_nit_P_1__ ; 
#line 340 "FSI_CSDdataloop.loci"
public:
#line 340 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop015_1280807839m60() {
#line 340 "FSI_CSDdataloop.loci"
       name_store("CFDIterationFinished{n,it-1}",L_CFDIterationFinished_nit_M_1__) ;
#line 340 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVel{n,it}",L_CSDnodesVel_nit__) ;
#line 340 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVelStar{n,it}",L_CSDnodesVelStar_nit__) ;
#line 340 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVel{n,it+1}",L_CSDnodesVel_nit_P_1__) ;
#line 340 "FSI_CSDdataloop.loci"
       input("CSDnodesVel{n,it},CSDnodesVelStar{n,it},CFDIterationFinished{n,it-1}") ;
#line 340 "FSI_CSDdataloop.loci"
       output("CSDnodesVel{n,it+1}") ;
#line 340 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 340 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 340 "FSI_CSDdataloop.loci"
    }
#line 340 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesVel_nit_P_1__).resize((*L_CSDnodesVelStar_nit__).size1(),1) ; std::fill( (*L_CSDnodesVel_nit_P_1__).data().begin(),(*L_CSDnodesVel_nit_P_1__).data().end(), 0. );
	
	if (*L_CFDIterationFinished_nit_M_1__) {
		(*L_CSDnodesVel_nit_P_1__) = (*L_CSDnodesVelStar_nit__) ;
	} else {
		(*L_CSDnodesVel_nit_P_1__) = (*L_CSDnodesVel_nit__) ;
	}
}    void compute(const Loci::sequence &seq) { 
#line 356 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 356 "FSI_CSDdataloop.loci"
    }
#line 356 "FSI_CSDdataloop.loci"
} ;
#line 356 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop015_1280807839m60> register_file_FSI_CSDdataloop015_1280807839m60 ;
#line 356 "FSI_CSDdataloop.loci"
}
#line 356 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop016_1280807839m61 : public Loci::blackbox_rule {
#line 356 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVel_nit__ ; 
#line 356 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVel_n_P_1__ ; 
#line 356 "FSI_CSDdataloop.loci"
public:
#line 356 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop016_1280807839m61() {
#line 356 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVel{n,it}",L_CSDnodesVel_nit__) ;
#line 356 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVel{n+1}",L_CSDnodesVel_n_P_1__) ;
#line 356 "FSI_CSDdataloop.loci"
       input("CSDnodesVel{n,it}") ;
#line 356 "FSI_CSDdataloop.loci"
       output("CSDnodesVel{n+1}") ;
#line 356 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 356 "FSI_CSDdataloop.loci"
       conditional("iterationFinished{n,it-1}") ;
#line 356 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 356 "FSI_CSDdataloop.loci"
    }
#line 356 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesVel_n_P_1__).resize((*L_CSDnodesVel_nit__).size1(),1) ;  std::fill( (*L_CSDnodesVel_n_P_1__).data().begin(),(*L_CSDnodesVel_n_P_1__).data().end(), 0. );
	(*L_CSDnodesVel_n_P_1__) = (*L_CSDnodesVel_nit__) ;
}    void compute(const Loci::sequence &seq) { 
#line 370 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 370 "FSI_CSDdataloop.loci"
    }
#line 370 "FSI_CSDdataloop.loci"
} ;
#line 370 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop016_1280807839m61> register_file_FSI_CSDdataloop016_1280807839m61 ;
#line 370 "FSI_CSDdataloop.loci"
}
#line 370 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop017_1280807839m62 : public Loci::blackbox_rule {
#line 370 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 370 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAcc_ic_ ; 
#line 370 "FSI_CSDdataloop.loci"
public:
#line 370 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop017_1280807839m62() {
#line 370 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 370 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAcc_ic",L_CSDnodesAcc_ic_) ;
#line 370 "FSI_CSDdataloop.loci"
       input("CSDnodes_ic") ;
#line 370 "FSI_CSDdataloop.loci"
       output("CSDnodesAcc_ic") ;
#line 370 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 370 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 370 "FSI_CSDdataloop.loci"
    }
#line 370 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesAcc_ic_).resize(6*(*L_CSDnodes_ic_).size1(),1) ; std::fill( (*L_CSDnodesAcc_ic_).data().begin(),(*L_CSDnodesAcc_ic_).data().end(), 0. );
	(*L_CSDnodesAcc_ic_).clear() ;
}    void compute(const Loci::sequence &seq) { 
#line 376 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 376 "FSI_CSDdataloop.loci"
    }
#line 376 "FSI_CSDdataloop.loci"
} ;
#line 376 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop017_1280807839m62> register_file_FSI_CSDdataloop017_1280807839m62 ;
#line 376 "FSI_CSDdataloop.loci"
}
#line 376 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop018_1280807839m62 : public Loci::blackbox_rule {
#line 376 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAcc_ic_ ; 
#line 376 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAcc_n_EQ_0__ ; 
#line 376 "FSI_CSDdataloop.loci"
public:
#line 376 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop018_1280807839m62() {
#line 376 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAcc_ic",L_CSDnodesAcc_ic_) ;
#line 376 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAcc{n=0}",L_CSDnodesAcc_n_EQ_0__) ;
#line 376 "FSI_CSDdataloop.loci"
       input("CSDnodesAcc_ic") ;
#line 376 "FSI_CSDdataloop.loci"
       output("CSDnodesAcc{n=0}") ;
#line 376 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 376 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 376 "FSI_CSDdataloop.loci"
    }
#line 376 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesAcc_n_EQ_0__).resize((*L_CSDnodesAcc_ic_).size1(),1) ;std::fill( (*L_CSDnodesAcc_n_EQ_0__).data().begin(),(*L_CSDnodesAcc_n_EQ_0__).data().end(), 0. );
	(*L_CSDnodesAcc_n_EQ_0__) = (*L_CSDnodesAcc_ic_) ;
}    void compute(const Loci::sequence &seq) { 
#line 383 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 383 "FSI_CSDdataloop.loci"
    }
#line 383 "FSI_CSDdataloop.loci"
} ;
#line 383 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop018_1280807839m62> register_file_FSI_CSDdataloop018_1280807839m62 ;
#line 383 "FSI_CSDdataloop.loci"
}
#line 383 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop019_1280807839m63 : public Loci::blackbox_rule {
#line 383 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAcc_n__ ; 
#line 383 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAcc_nit_EQ_0__ ; 
#line 383 "FSI_CSDdataloop.loci"
public:
#line 383 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop019_1280807839m63() {
#line 383 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAcc{n}",L_CSDnodesAcc_n__) ;
#line 383 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAcc{n,it=0}",L_CSDnodesAcc_nit_EQ_0__) ;
#line 383 "FSI_CSDdataloop.loci"
       input("CSDnodesAcc{n}") ;
#line 383 "FSI_CSDdataloop.loci"
       output("CSDnodesAcc{n,it=0}") ;
#line 383 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 383 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 383 "FSI_CSDdataloop.loci"
    }
#line 383 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesAcc_nit_EQ_0__).resize((*L_CSDnodesAcc_n__).size1(),1) ;std::fill( (*L_CSDnodesAcc_nit_EQ_0__).data().begin(),(*L_CSDnodesAcc_nit_EQ_0__).data().end(), 0. );
	(*L_CSDnodesAcc_nit_EQ_0__) = (*L_CSDnodesAcc_n__) ;
	//if (Loci::MPI_rank==0) cout << "Inside CSDnodesAcc{n,it=0}" << endl;
}    void compute(const Loci::sequence &seq) { 
#line 390 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 390 "FSI_CSDdataloop.loci"
    }
#line 390 "FSI_CSDdataloop.loci"
} ;
#line 390 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop019_1280807839m63> register_file_FSI_CSDdataloop019_1280807839m63 ;
#line 390 "FSI_CSDdataloop.loci"
}
#line 390 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop020_1280807839m64 : public Loci::blackbox_rule {
#line 390 "FSI_CSDdataloop.loci"
    Loci::const_param<bool>  L_CFDIterationFinished_nit_M_1__ ; 
#line 390 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAcc_nit__ ; 
#line 390 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAccStar_nit__ ; 
#line 390 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAcc_nit_P_1__ ; 
#line 390 "FSI_CSDdataloop.loci"
public:
#line 390 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop020_1280807839m64() {
#line 390 "FSI_CSDdataloop.loci"
       name_store("CFDIterationFinished{n,it-1}",L_CFDIterationFinished_nit_M_1__) ;
#line 390 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAcc{n,it}",L_CSDnodesAcc_nit__) ;
#line 390 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAccStar{n,it}",L_CSDnodesAccStar_nit__) ;
#line 390 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAcc{n,it+1}",L_CSDnodesAcc_nit_P_1__) ;
#line 390 "FSI_CSDdataloop.loci"
       input("CSDnodesAcc{n,it},CSDnodesAccStar{n,it},CFDIterationFinished{n,it-1}") ;
#line 390 "FSI_CSDdataloop.loci"
       output("CSDnodesAcc{n,it+1}") ;
#line 390 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 390 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 390 "FSI_CSDdataloop.loci"
    }
#line 390 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesAcc_nit_P_1__).resize((*L_CSDnodesAccStar_nit__).size1(),1) ;std::fill( (*L_CSDnodesAcc_nit_P_1__).data().begin(),(*L_CSDnodesAcc_nit_P_1__).data().end(), 0. );
	if (*L_CFDIterationFinished_nit_M_1__) {
		(*L_CSDnodesAcc_nit_P_1__) = (*L_CSDnodesAccStar_nit__) ;
	} else {
		(*L_CSDnodesAcc_nit_P_1__) = (*L_CSDnodesAcc_nit__) ;
	}
}    void compute(const Loci::sequence &seq) { 
#line 405 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 405 "FSI_CSDdataloop.loci"
    }
#line 405 "FSI_CSDdataloop.loci"
} ;
#line 405 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop020_1280807839m64> register_file_FSI_CSDdataloop020_1280807839m64 ;
#line 405 "FSI_CSDdataloop.loci"
}
#line 405 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop021_1280807839m65 : public Loci::blackbox_rule {
#line 405 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAcc_nit__ ; 
#line 405 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAcc_n_P_1__ ; 
#line 405 "FSI_CSDdataloop.loci"
public:
#line 405 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop021_1280807839m65() {
#line 405 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAcc{n,it}",L_CSDnodesAcc_nit__) ;
#line 405 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAcc{n+1}",L_CSDnodesAcc_n_P_1__) ;
#line 405 "FSI_CSDdataloop.loci"
       input("CSDnodesAcc{n,it}") ;
#line 405 "FSI_CSDdataloop.loci"
       output("CSDnodesAcc{n+1}") ;
#line 405 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 405 "FSI_CSDdataloop.loci"
       conditional("iterationFinished{n,it-1}") ;
#line 405 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 405 "FSI_CSDdataloop.loci"
    }
#line 405 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesAcc_n_P_1__).resize((*L_CSDnodesAcc_nit__).size1(),1) ; std::fill( (*L_CSDnodesAcc_n_P_1__).data().begin(),(*L_CSDnodesAcc_n_P_1__).data().end(), 0. );
	(*L_CSDnodesAcc_n_P_1__) = (*L_CSDnodesAcc_nit__) ;
}    void compute(const Loci::sequence &seq) { 
#line 420 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 420 "FSI_CSDdataloop.loci"
    }
#line 420 "FSI_CSDdataloop.loci"
} ;
#line 420 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop021_1280807839m65> register_file_FSI_CSDdataloop021_1280807839m65 ;
#line 420 "FSI_CSDdataloop.loci"
}
#line 420 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop022_1280807839m66 : public Loci::blackbox_rule {
#line 420 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 420 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePre_ic_ ; 
#line 420 "FSI_CSDdataloop.loci"
public:
#line 420 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop022_1280807839m66() {
#line 420 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 420 "FSI_CSDdataloop.loci"
       name_store("CSDForcePre_ic",L_CSDForcePre_ic_) ;
#line 420 "FSI_CSDdataloop.loci"
       input("CSDnodes_ic") ;
#line 420 "FSI_CSDdataloop.loci"
       output("CSDForcePre_ic") ;
#line 420 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 420 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 420 "FSI_CSDdataloop.loci"
    }
#line 420 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDForcePre_ic_).resize(6*(*L_CSDnodes_ic_).size1(),1) ; std::fill( (*L_CSDForcePre_ic_).data().begin(),(*L_CSDForcePre_ic_).data().end(), 0. );
	(*L_CSDForcePre_ic_).clear();
}    void compute(const Loci::sequence &seq) { 
#line 426 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 426 "FSI_CSDdataloop.loci"
    }
#line 426 "FSI_CSDdataloop.loci"
} ;
#line 426 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop022_1280807839m66> register_file_FSI_CSDdataloop022_1280807839m66 ;
#line 426 "FSI_CSDdataloop.loci"
}
#line 426 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop023_1280807839m66 : public Loci::blackbox_rule {
#line 426 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePre_ic_ ; 
#line 426 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePre_n_EQ_0__ ; 
#line 426 "FSI_CSDdataloop.loci"
public:
#line 426 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop023_1280807839m66() {
#line 426 "FSI_CSDdataloop.loci"
       name_store("CSDForcePre_ic",L_CSDForcePre_ic_) ;
#line 426 "FSI_CSDdataloop.loci"
       name_store("CSDForcePre{n=0}",L_CSDForcePre_n_EQ_0__) ;
#line 426 "FSI_CSDdataloop.loci"
       input("CSDForcePre_ic") ;
#line 426 "FSI_CSDdataloop.loci"
       output("CSDForcePre{n=0}") ;
#line 426 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 426 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 426 "FSI_CSDdataloop.loci"
    }
#line 426 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDForcePre_n_EQ_0__).resize((*L_CSDForcePre_ic_).size1(),1) ;std::fill((*L_CSDForcePre_n_EQ_0__).data().begin(),(*L_CSDForcePre_n_EQ_0__).data().end(), 0. );
	(*L_CSDForcePre_n_EQ_0__) = (*L_CSDForcePre_ic_) ;
}    void compute(const Loci::sequence &seq) { 
#line 433 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 433 "FSI_CSDdataloop.loci"
    }
#line 433 "FSI_CSDdataloop.loci"
} ;
#line 433 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop023_1280807839m66> register_file_FSI_CSDdataloop023_1280807839m66 ;
#line 433 "FSI_CSDdataloop.loci"
}
#line 433 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop024_1280807839m67 : public Loci::blackbox_rule {
#line 433 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePre_n__ ; 
#line 433 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePre_nit_EQ_0__ ; 
#line 433 "FSI_CSDdataloop.loci"
public:
#line 433 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop024_1280807839m67() {
#line 433 "FSI_CSDdataloop.loci"
       name_store("CSDForcePre{n}",L_CSDForcePre_n__) ;
#line 433 "FSI_CSDdataloop.loci"
       name_store("CSDForcePre{n,it=0}",L_CSDForcePre_nit_EQ_0__) ;
#line 433 "FSI_CSDdataloop.loci"
       input("CSDForcePre{n}") ;
#line 433 "FSI_CSDdataloop.loci"
       output("CSDForcePre{n,it=0}") ;
#line 433 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 433 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 433 "FSI_CSDdataloop.loci"
    }
#line 433 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDForcePre_nit_EQ_0__).resize((*L_CSDForcePre_n__).size1(),1) ;std::fill( (*L_CSDForcePre_nit_EQ_0__).data().begin(),(*L_CSDForcePre_nit_EQ_0__).data().end(), 0. );
	(*L_CSDForcePre_nit_EQ_0__) = (*L_CSDForcePre_n__) ;
//	if (Loci::MPI_rank==0) cout << "Inside CSDForcePre{n,it=0}" << endl;
}    void compute(const Loci::sequence &seq) { 
#line 440 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 440 "FSI_CSDdataloop.loci"
    }
#line 440 "FSI_CSDdataloop.loci"
} ;
#line 440 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop024_1280807839m67> register_file_FSI_CSDdataloop024_1280807839m67 ;
#line 440 "FSI_CSDdataloop.loci"
}
#line 440 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop025_1280807839m68 : public Loci::blackbox_rule {
#line 440 "FSI_CSDdataloop.loci"
    Loci::const_param<bool>  L_CFDIterationFinished_nit_M_1__ ; 
#line 440 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePre_nit__ ; 
#line 440 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePreStar_nit__ ; 
#line 440 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePre_nit_P_1__ ; 
#line 440 "FSI_CSDdataloop.loci"
public:
#line 440 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop025_1280807839m68() {
#line 440 "FSI_CSDdataloop.loci"
       name_store("CFDIterationFinished{n,it-1}",L_CFDIterationFinished_nit_M_1__) ;
#line 440 "FSI_CSDdataloop.loci"
       name_store("CSDForcePre{n,it}",L_CSDForcePre_nit__) ;
#line 440 "FSI_CSDdataloop.loci"
       name_store("CSDForcePreStar{n,it}",L_CSDForcePreStar_nit__) ;
#line 440 "FSI_CSDdataloop.loci"
       name_store("CSDForcePre{n,it+1}",L_CSDForcePre_nit_P_1__) ;
#line 440 "FSI_CSDdataloop.loci"
       input("CSDForcePre{n,it},CSDForcePreStar{n,it},CFDIterationFinished{n,it-1}") ;
#line 440 "FSI_CSDdataloop.loci"
       output("CSDForcePre{n,it+1}") ;
#line 440 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 440 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 440 "FSI_CSDdataloop.loci"
    }
#line 440 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDForcePre_nit_P_1__).resize((*L_CSDForcePreStar_nit__).size1(),1) ;std::fill( (*L_CSDForcePre_nit_P_1__).data().begin(),(*L_CSDForcePre_nit_P_1__).data().end(), 0. );
	
	if (*L_CFDIterationFinished_nit_M_1__) {
		(*L_CSDForcePre_nit_P_1__) = (*L_CSDForcePreStar_nit__) ;
	} else {
		(*L_CSDForcePre_nit_P_1__) = (*L_CSDForcePre_nit__) ;
	}
}    void compute(const Loci::sequence &seq) { 
#line 456 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 456 "FSI_CSDdataloop.loci"
    }
#line 456 "FSI_CSDdataloop.loci"
} ;
#line 456 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop025_1280807839m68> register_file_FSI_CSDdataloop025_1280807839m68 ;
#line 456 "FSI_CSDdataloop.loci"
}
#line 456 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop026_1280807839m69 : public Loci::blackbox_rule {
#line 456 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePre_nit__ ; 
#line 456 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePre_n_P_1__ ; 
#line 456 "FSI_CSDdataloop.loci"
public:
#line 456 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop026_1280807839m69() {
#line 456 "FSI_CSDdataloop.loci"
       name_store("CSDForcePre{n,it}",L_CSDForcePre_nit__) ;
#line 456 "FSI_CSDdataloop.loci"
       name_store("CSDForcePre{n+1}",L_CSDForcePre_n_P_1__) ;
#line 456 "FSI_CSDdataloop.loci"
       input("CSDForcePre{n,it}") ;
#line 456 "FSI_CSDdataloop.loci"
       output("CSDForcePre{n+1}") ;
#line 456 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 456 "FSI_CSDdataloop.loci"
       conditional("iterationFinished{n,it-1}") ;
#line 456 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 456 "FSI_CSDdataloop.loci"
    }
#line 456 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDForcePre_n_P_1__).resize((*L_CSDForcePre_nit__).size1(),1) ;std::fill( (*L_CSDForcePre_n_P_1__).data().begin(),(*L_CSDForcePre_n_P_1__).data().end(), 0. );
	(*L_CSDForcePre_n_P_1__) = (*L_CSDForcePre_nit__) ;
}    void compute(const Loci::sequence &seq) { 
#line 470 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 470 "FSI_CSDdataloop.loci"
    }
#line 470 "FSI_CSDdataloop.loci"
} ;
#line 470 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop026_1280807839m69> register_file_FSI_CSDdataloop026_1280807839m69 ;
#line 470 "FSI_CSDdataloop.loci"
}
#line 470 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop027_1280807839m70 : public Loci::blackbox_rule {
#line 470 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<int,ublas::column_major> >  L_CSDConnectivity_ ; 
#line 470 "FSI_CSDdataloop.loci"
    Loci::blackbox<boost::multi_array<real,4> >  L_CSDnodesSys_ic_ ; 
#line 470 "FSI_CSDdataloop.loci"
public:
#line 470 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop027_1280807839m70() {
#line 470 "FSI_CSDdataloop.loci"
       name_store("CSDConnectivity",L_CSDConnectivity_) ;
#line 470 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSys_ic",L_CSDnodesSys_ic_) ;
#line 470 "FSI_CSDdataloop.loci"
       input("CSDConnectivity") ;
#line 470 "FSI_CSDdataloop.loci"
       output("CSDnodesSys_ic") ;
#line 470 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 470 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 470 "FSI_CSDdataloop.loci"
    }
#line 470 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	 size_t numEl = (*L_CSDConnectivity_).size1() ;
	// if (Loci::MPI_rank==0) cout << "ic: numEL=" << numEl << endl;
	 (*L_CSDnodesSys_ic_).resize(boost::extents[3][3][numEl][3]);
	 //(*$CSDnodesSys_ic).resize(boost::extents[3][3][3][3]);
	 std::fill( (*L_CSDnodesSys_ic_).origin(), (*L_CSDnodesSys_ic_).origin() + (*L_CSDnodesSys_ic_).size(), 0. ); // initialize with zero
	 	}    void compute(const Loci::sequence &seq) { 
#line 479 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 479 "FSI_CSDdataloop.loci"
    }
#line 479 "FSI_CSDdataloop.loci"
} ;
#line 479 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop027_1280807839m70> register_file_FSI_CSDdataloop027_1280807839m70 ;
#line 479 "FSI_CSDdataloop.loci"
}
#line 479 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop028_1280807839m70 : public Loci::blackbox_rule {
#line 479 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<int,ublas::column_major> >  L_CSDConnectivity_ ; 
#line 479 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<boost::multi_array<real,4> >  L_CSDnodesSys_ic_ ; 
#line 479 "FSI_CSDdataloop.loci"
    Loci::blackbox<boost::multi_array<real,4> >  L_CSDnodesSys_n_EQ_0__ ; 
#line 479 "FSI_CSDdataloop.loci"
public:
#line 479 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop028_1280807839m70() {
#line 479 "FSI_CSDdataloop.loci"
       name_store("CSDConnectivity",L_CSDConnectivity_) ;
#line 479 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSys_ic",L_CSDnodesSys_ic_) ;
#line 479 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSys{n=0}",L_CSDnodesSys_n_EQ_0__) ;
#line 479 "FSI_CSDdataloop.loci"
       input("CSDnodesSys_ic,CSDConnectivity") ;
#line 479 "FSI_CSDdataloop.loci"
       output("CSDnodesSys{n=0}") ;
#line 479 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 479 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 479 "FSI_CSDdataloop.loci"
    }
#line 479 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	size_t numEl = (*L_CSDConnectivity_).size1() ;
	//if (Loci::MPI_rank==0) cout << "{n=0} numEL=" << numEl << endl;
	//(*$CSDnodesSys{n=0}).resize(boost::extents[0][0][0][0]) ;
	(*L_CSDnodesSys_n_EQ_0__).resize(boost::extents[3][3][numEl][3]) ;  std::fill( (*L_CSDnodesSys_n_EQ_0__).origin(), (*L_CSDnodesSys_n_EQ_0__).origin() + (*L_CSDnodesSys_n_EQ_0__).size(), 0. );
	(*L_CSDnodesSys_n_EQ_0__) = (*L_CSDnodesSys_ic_) ;
	// std::fill( (*$CSDnodesSys{n=0}).origin(), (*$CSDnodesSys{n=0}).origin() + (*$CSDnodesSys{n=0}).size(), 0. ); // initialize with zero
}    void compute(const Loci::sequence &seq) { 
#line 490 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 490 "FSI_CSDdataloop.loci"
    }
#line 490 "FSI_CSDdataloop.loci"
} ;
#line 490 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop028_1280807839m70> register_file_FSI_CSDdataloop028_1280807839m70 ;
#line 490 "FSI_CSDdataloop.loci"
}
#line 490 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop029_1280807839m71 : public Loci::blackbox_rule {
#line 490 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<boost::multi_array<real,4> >  L_CSDnodesSys_n__ ; 
#line 490 "FSI_CSDdataloop.loci"
    Loci::blackbox<boost::multi_array<real,4> >  L_CSDnodesSys_nit_EQ_0__ ; 
#line 490 "FSI_CSDdataloop.loci"
public:
#line 490 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop029_1280807839m71() {
#line 490 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSys{n}",L_CSDnodesSys_n__) ;
#line 490 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSys{n,it=0}",L_CSDnodesSys_nit_EQ_0__) ;
#line 490 "FSI_CSDdataloop.loci"
       input("CSDnodesSys{n}") ;
#line 490 "FSI_CSDdataloop.loci"
       output("CSDnodesSys{n,it=0}") ;
#line 490 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 490 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 490 "FSI_CSDdataloop.loci"
    }
#line 490 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	size_t numEl = (*L_CSDnodesSys_n__).shape()[2] ;
	//if (Loci::MPI_rank==0) cout << "{n,it=0} numEL=" << numEl << endl;
	(*L_CSDnodesSys_nit_EQ_0__).resize(boost::extents[3][3][numEl][3]) ; std::fill( (*L_CSDnodesSys_nit_EQ_0__).origin(), (*L_CSDnodesSys_nit_EQ_0__).origin() + (*L_CSDnodesSys_nit_EQ_0__).size(), 0. );
	(*L_CSDnodesSys_nit_EQ_0__) = (*L_CSDnodesSys_n__) ;
	//if (Loci::MPI_rank==0) cout << "Inside CSDnodesSys{n,it=0}" << endl;
}    void compute(const Loci::sequence &seq) { 
#line 499 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 499 "FSI_CSDdataloop.loci"
    }
#line 499 "FSI_CSDdataloop.loci"
} ;
#line 499 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop029_1280807839m71> register_file_FSI_CSDdataloop029_1280807839m71 ;
#line 499 "FSI_CSDdataloop.loci"
}
#line 499 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop030_1280807839m72 : public Loci::blackbox_rule {
#line 499 "FSI_CSDdataloop.loci"
    Loci::const_param<bool>  L_CFDIterationFinished_nit_M_1__ ; 
#line 499 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<boost::multi_array<real,4> >  L_CSDnodesSys_nit__ ; 
#line 499 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<boost::multi_array<real,4> >  L_CSDnodesSysStar_nit__ ; 
#line 499 "FSI_CSDdataloop.loci"
    Loci::blackbox<boost::multi_array<real,4> >  L_CSDnodesSys_nit_P_1__ ; 
#line 499 "FSI_CSDdataloop.loci"
public:
#line 499 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop030_1280807839m72() {
#line 499 "FSI_CSDdataloop.loci"
       name_store("CFDIterationFinished{n,it-1}",L_CFDIterationFinished_nit_M_1__) ;
#line 499 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSys{n,it}",L_CSDnodesSys_nit__) ;
#line 499 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSysStar{n,it}",L_CSDnodesSysStar_nit__) ;
#line 499 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSys{n,it+1}",L_CSDnodesSys_nit_P_1__) ;
#line 499 "FSI_CSDdataloop.loci"
       input("CSDnodesSys{n,it},CSDnodesSysStar{n,it},CFDIterationFinished{n,it-1}") ;
#line 499 "FSI_CSDdataloop.loci"
       output("CSDnodesSys{n,it+1}") ;
#line 499 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 499 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 499 "FSI_CSDdataloop.loci"
    }
#line 499 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	size_t numEl = (*L_CSDnodesSysStar_nit__).shape()[2] ;
//	if (Loci::MPI_rank==0) cout << "Star{n,it} numEL=" << numEl << endl;
	(*L_CSDnodesSys_nit_P_1__).resize(boost::extents[3][3][numEl][3]) ; std::fill( (*L_CSDnodesSys_nit_P_1__).origin(), (*L_CSDnodesSys_nit_P_1__).origin() + (*L_CSDnodesSys_nit_P_1__).size(), 0. );
		if (*L_CFDIterationFinished_nit_M_1__) {
			(*L_CSDnodesSys_nit_P_1__) = (*L_CSDnodesSysStar_nit__) ;
		} else {
			(*L_CSDnodesSys_nit_P_1__) = (*L_CSDnodesSys_nit__) ;
		}
}    void compute(const Loci::sequence &seq) { 
#line 518 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 518 "FSI_CSDdataloop.loci"
    }
#line 518 "FSI_CSDdataloop.loci"
} ;
#line 518 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop030_1280807839m72> register_file_FSI_CSDdataloop030_1280807839m72 ;
#line 518 "FSI_CSDdataloop.loci"
}
#line 518 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop031_1280807839m73 : public Loci::blackbox_rule {
#line 518 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<boost::multi_array<real,4> >  L_CSDnodesSys_nit__ ; 
#line 518 "FSI_CSDdataloop.loci"
    Loci::blackbox<boost::multi_array<real,4> >  L_CSDnodesSys_n_P_1__ ; 
#line 518 "FSI_CSDdataloop.loci"
public:
#line 518 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop031_1280807839m73() {
#line 518 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSys{n,it}",L_CSDnodesSys_nit__) ;
#line 518 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSys{n+1}",L_CSDnodesSys_n_P_1__) ;
#line 518 "FSI_CSDdataloop.loci"
       input("CSDnodesSys{n,it}") ;
#line 518 "FSI_CSDdataloop.loci"
       output("CSDnodesSys{n+1}") ;
#line 518 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 518 "FSI_CSDdataloop.loci"
       conditional("iterationFinished{n,it-1}") ;
#line 518 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 518 "FSI_CSDdataloop.loci"
    }
#line 518 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	size_t numEl = (*L_CSDnodesSys_nit__).shape()[2] ;
//	if (Loci::MPI_rank==0) cout << "{n+1} numEL=" << numEl << endl;
	(*L_CSDnodesSys_n_P_1__).resize(boost::extents[3][3][numEl][3]) ;std::fill( (*L_CSDnodesSys_n_P_1__).origin(), (*L_CSDnodesSys_n_P_1__).origin() + (*L_CSDnodesSys_n_P_1__).size(), 0. );
	(*L_CSDnodesSys_n_P_1__) = (*L_CSDnodesSys_nit__) ;
//	if (Loci::MPI_rank==0) cout << "inside CSDnodesSys{n+1}" << endl ; 
}    void compute(const Loci::sequence &seq) { 
#line 537 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 537 "FSI_CSDdataloop.loci"
    }
#line 537 "FSI_CSDdataloop.loci"
} ;
#line 537 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop031_1280807839m73> register_file_FSI_CSDdataloop031_1280807839m73 ;
#line 537 "FSI_CSDdataloop.loci"
}
#line 537 "FSI_CSDdataloop.loci"
// $type CFDIterationFinished param<bool> 
// $type stime param<real> 
// $type ncycle param<int> 

namespace {class file_FSI_CSDdataloop032_1280807839m74 : public Loci::singleton_rule {
#line 541 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_stime_n__ ; 
#line 541 "FSI_CSDdataloop.loci"
    Loci::param<real>  L_stime_nit_EQ_0__ ; 
#line 541 "FSI_CSDdataloop.loci"
public:
#line 541 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop032_1280807839m74() {
#line 541 "FSI_CSDdataloop.loci"
       name_store("stime{n}",L_stime_n__) ;
#line 541 "FSI_CSDdataloop.loci"
       name_store("stime{n,it=0}",L_stime_nit_EQ_0__) ;
#line 541 "FSI_CSDdataloop.loci"
       input("stime{n}") ;
#line 541 "FSI_CSDdataloop.loci"
       output("stime{n,it=0}") ;
#line 541 "FSI_CSDdataloop.loci"
    }
#line 541 "FSI_CSDdataloop.loci"
    void compute(const Loci::sequence &seq) { 	
//	if (Loci::MPI_rank==0) cout << "inside stime{n,it=0} in" << endl ; 
	(*L_stime_nit_EQ_0__)=(*L_stime_n__) ;
//	if (Loci::MPI_rank==0) cout << "inside stime{n,it=0} out" << endl ; 
}} ;
#line 545 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop032_1280807839m74> register_file_FSI_CSDdataloop032_1280807839m74 ;
#line 545 "FSI_CSDdataloop.loci"
}
#line 545 "FSI_CSDdataloop.loci"
	 

namespace {class file_FSI_CSDdataloop033_1280807839m75 : public Loci::singleton_rule {
#line 547 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_stime_nit__ ; 
#line 547 "FSI_CSDdataloop.loci"
    Loci::param<real>  L_stime_nit_P_1__ ; 
#line 547 "FSI_CSDdataloop.loci"
public:
#line 547 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop033_1280807839m75() {
#line 547 "FSI_CSDdataloop.loci"
       name_store("stime{n,it}",L_stime_nit__) ;
#line 547 "FSI_CSDdataloop.loci"
       name_store("stime{n,it+1}",L_stime_nit_P_1__) ;
#line 547 "FSI_CSDdataloop.loci"
       input("stime{n,it}") ;
#line 547 "FSI_CSDdataloop.loci"
       output("stime{n,it+1}") ;
#line 547 "FSI_CSDdataloop.loci"
    }
#line 547 "FSI_CSDdataloop.loci"
    void compute(const Loci::sequence &seq) { 
	(*L_stime_nit_P_1__)=(*L_stime_nit__) ;
}} ;
#line 549 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop033_1280807839m75> register_file_FSI_CSDdataloop033_1280807839m75 ;
#line 549 "FSI_CSDdataloop.loci"
}
#line 549 "FSI_CSDdataloop.loci"
	 

namespace {class file_FSI_CSDdataloop034_1280807839m75 : public Loci::singleton_rule {
#line 551 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_ncycle_n__ ; 
#line 551 "FSI_CSDdataloop.loci"
    Loci::param<int>  L_ncycle_nit_EQ_0__ ; 
#line 551 "FSI_CSDdataloop.loci"
public:
#line 551 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop034_1280807839m75() {
#line 551 "FSI_CSDdataloop.loci"
       name_store("ncycle{n}",L_ncycle_n__) ;
#line 551 "FSI_CSDdataloop.loci"
       name_store("ncycle{n,it=0}",L_ncycle_nit_EQ_0__) ;
#line 551 "FSI_CSDdataloop.loci"
       input("ncycle{n}") ;
#line 551 "FSI_CSDdataloop.loci"
       output("ncycle{n,it=0}") ;
#line 551 "FSI_CSDdataloop.loci"
    }
#line 551 "FSI_CSDdataloop.loci"
    void compute(const Loci::sequence &seq) { 	
//	if (Loci::MPI_rank==0) cout << "inside ncycle{n,it=0} in" << endl ; 
	(*L_ncycle_nit_EQ_0__)=(*L_ncycle_n__) ;
//	if (Loci::MPI_rank==0) cout << "inside ncycle{n,it=0} out" << endl ; 
}} ;
#line 555 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop034_1280807839m75> register_file_FSI_CSDdataloop034_1280807839m75 ;
#line 555 "FSI_CSDdataloop.loci"
}
#line 555 "FSI_CSDdataloop.loci"
	 

namespace {class file_FSI_CSDdataloop035_1280807839m76 : public Loci::singleton_rule {
#line 557 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_ncycle_nit__ ; 
#line 557 "FSI_CSDdataloop.loci"
    Loci::param<int>  L_ncycle_nit_P_1__ ; 
#line 557 "FSI_CSDdataloop.loci"
public:
#line 557 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop035_1280807839m76() {
#line 557 "FSI_CSDdataloop.loci"
       name_store("ncycle{n,it}",L_ncycle_nit__) ;
#line 557 "FSI_CSDdataloop.loci"
       name_store("ncycle{n,it+1}",L_ncycle_nit_P_1__) ;
#line 557 "FSI_CSDdataloop.loci"
       input("ncycle{n,it}") ;
#line 557 "FSI_CSDdataloop.loci"
       output("ncycle{n,it+1}") ;
#line 557 "FSI_CSDdataloop.loci"
    }
#line 557 "FSI_CSDdataloop.loci"
    void compute(const Loci::sequence &seq) { 
	(*L_ncycle_nit_P_1__)=(*L_ncycle_nit__) ;
}} ;
#line 559 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop035_1280807839m76> register_file_FSI_CSDdataloop035_1280807839m76 ;
#line 559 "FSI_CSDdataloop.loci"
}
#line 559 "FSI_CSDdataloop.loci"
	 

//$rule singleton(stime{n+1}<-stime{n,it}),conditional(iterationFinished{n,it-1}) {
//	$stime{n+1}=$stime{n,it} ;
//}	 

// $type timeStep param<real> 
//$rule blackbox(CSDdisplacementsStar{n,it},CSDnodesDispStar,CSDnodesVelStar,CSDnodesAccStar,CSDForcePreStar,CSDnodesSysStar<-
//							CSDnodes_ic,CSDnodes,CSDnodesDisp,CSDnodesVel,CSDnodesAcc,CSDForce,CSDForcePre,CSDConnectivity,CSDBCdof,CSDBCZeroConstraint,CSDnodesSys,
//							stime,ncycle,timeStep,CSDstartingTimeStep,
//							CSDE1,CSDE2,CSDnu12,CSDnu21,CSDG12,CSDthicknessStructure,CSDrhoStructure,
//							CSDintegrationScheme,CSDdelta,CSDswitchStiffening,CSDgenAlphaCoeff,CSDnewmarkGammaCoeff,CSDnewmarkBetaCoeff,CSDdampingCoeff1,CSDdampingCoeff2,
//							CSDexcitationType,CSDflappingType,CSDplungingType,
//							CSDfrequency,CSDplungeAmplitudeX,CSDplungeAmplitudeY,CSDplungeAmplitudeZ,CSDflappingAmplitudeX,CSDflappingAmplitudeY,CSDflappingAmplitudeZ,
//							CFDIterationFinished{n,it-1}),
//						constraint(FSICoupling,FSINLAMS),option(disable_threading), prelude {
// $rule blackbox(CSDdisplacementsStar{n},CSDnodesDispStar{n},CSDnodesVelStar{n},CSDnodesAccStar{n},CSDForcePreStar{n,it},CSDnodesSysStar{n}<-
// 							CSDnodes_ic,CSDnodes{n},CSDnodesDisp{n},CSDnodesVel{n},CSDnodesAcc{n},CSDForce{n},CSDForcePre{n},CSDConnectivity,CSDBCdof,CSDBCZeroConstraint,CSDnodesSys{n},
// 							stime{n},ncycle{n},timeStep{n},CSDstartingTimeStep,CSDtipNode,
// 							CSDE1{n},CSDE2{n},CSDnu12{n},CSDnu21{n},CSDG12{n},CSDthicknessStructure{n},CSDrhoStructure{n},
// 							CSDintegrationScheme{n},CSDdelta{n},CSDswitchStiffening{n},CSDgenAlphaCoeff{n},CSDnewmarkGammaCoeff{n},CSDnewmarkBetaCoeff{n},CSDdampingCoeff1{n},CSDdampingCoeff2{n},
// 							CSDflappingType{n},CSDplungingType{n},
// 							CSDfrequency{n},CSDplungeAmplitudeX{n},CSDplungeAmplitudeY{n},CSDplungeAmplitudeZ{n},CSDflappingAmplitudeX{n},CSDflappingAmplitudeY{n},CSDflappingAmplitudeZ{n},
// 							itfsi{n},CSDdimension),
// 						constraint(FSICoupling{n},FSINLAMS{n}),option(disable_threading), prelude {

  
// NLAMS not yet modified
  namespace {class file_FSI_CSDdataloop036_1280807839m77 : public Loci::blackbox_rule {
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDE1_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDE2_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDnu12_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDnu21_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDG12_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDrhoStructure_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDthicknessStructure_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDintegrationScheme_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDdelta_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDswitchStiffening_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDflappingType_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDplungingType_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDfrequency_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDgenAlphaCoeff_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDnewmarkGammaCoeff_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDnewmarkBetaCoeff_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDdampingCoeff1_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDdampingCoeff2_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDplungeAmplitudeX_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDplungeAmplitudeY_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDplungeAmplitudeZ_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDflappingAmplitudeX_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDflappingAmplitudeY_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDflappingAmplitudeZ_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDstartingTimeStep_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDtipNode_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDdimension_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<int,ublas::column_major> >  L_CSDConnectivity_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<vector<int> >  L_CSDBCdof_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<vector<real> >  L_CSDBCZeroConstraint_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_n__ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<bool>  L_CFDIterationFinished_nit_M_1__ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDisp_n__ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVel_n__ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAcc_n__ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePre_n__ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<boost::multi_array<real,4> >  L_CSDnodesSys_n__ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_stime_n__ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_ncycle_n__ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_timeStep_ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForce_nit__ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_itfsi_nit__ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDdisplacementsStar_nit__ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDispStar_nit__ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVelStar_nit__ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAccStar_nit__ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePreStar_nit__ ; 
#line 595 "FSI_CSDdataloop.loci"
    Loci::blackbox<boost::multi_array<real,4> >  L_CSDnodesSysStar_nit__ ; 
#line 595 "FSI_CSDdataloop.loci"
public:
#line 595 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop036_1280807839m77() {
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDE1",L_CSDE1_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDE2",L_CSDE2_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDnu12",L_CSDnu12_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDnu21",L_CSDnu21_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDG12",L_CSDG12_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDrhoStructure",L_CSDrhoStructure_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDthicknessStructure",L_CSDthicknessStructure_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDintegrationScheme",L_CSDintegrationScheme_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDdelta",L_CSDdelta_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDswitchStiffening",L_CSDswitchStiffening_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDflappingType",L_CSDflappingType_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDplungingType",L_CSDplungingType_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDfrequency",L_CSDfrequency_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDgenAlphaCoeff",L_CSDgenAlphaCoeff_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDnewmarkGammaCoeff",L_CSDnewmarkGammaCoeff_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDnewmarkBetaCoeff",L_CSDnewmarkBetaCoeff_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDdampingCoeff1",L_CSDdampingCoeff1_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDdampingCoeff2",L_CSDdampingCoeff2_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDplungeAmplitudeX",L_CSDplungeAmplitudeX_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDplungeAmplitudeY",L_CSDplungeAmplitudeY_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDplungeAmplitudeZ",L_CSDplungeAmplitudeZ_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDflappingAmplitudeX",L_CSDflappingAmplitudeX_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDflappingAmplitudeY",L_CSDflappingAmplitudeY_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDflappingAmplitudeZ",L_CSDflappingAmplitudeZ_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDstartingTimeStep",L_CSDstartingTimeStep_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDtipNode",L_CSDtipNode_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDdimension",L_CSDdimension_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDConnectivity",L_CSDConnectivity_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDBCdof",L_CSDBCdof_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDBCZeroConstraint",L_CSDBCZeroConstraint_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDnodes{n}",L_CSDnodes_n__) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDdisplacementsStar{n,it}",L_CSDdisplacementsStar_nit__) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CFDIterationFinished{n,it-1}",L_CFDIterationFinished_nit_M_1__) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDisp{n}",L_CSDnodesDisp_n__) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDispStar{n,it}",L_CSDnodesDispStar_nit__) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVel{n}",L_CSDnodesVel_n__) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVelStar{n,it}",L_CSDnodesVelStar_nit__) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAcc{n}",L_CSDnodesAcc_n__) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAccStar{n,it}",L_CSDnodesAccStar_nit__) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDForcePre{n}",L_CSDForcePre_n__) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDForcePreStar{n,it}",L_CSDForcePreStar_nit__) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSys{n}",L_CSDnodesSys_n__) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSysStar{n,it}",L_CSDnodesSysStar_nit__) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("stime{n}",L_stime_n__) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("ncycle{n}",L_ncycle_n__) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("timeStep",L_timeStep_) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("CSDForce{n,it}",L_CSDForce_nit__) ;
#line 595 "FSI_CSDdataloop.loci"
       name_store("itfsi{n,it}",L_itfsi_nit__) ;
#line 595 "FSI_CSDdataloop.loci"
       input("							CSDnodes_ic,CSDnodes{n},CSDnodesDisp{n},CSDnodesVel{n},CSDnodesAcc{n},CSDForce{n,it},CSDForcePre{n},CSDConnectivity,CSDBCdof,CSDBCZeroConstraint,CSDnodesSys{n},							stime{n},ncycle{n},timeStep,CSDstartingTimeStep,CSDtipNode,							CSDE1,CSDE2,CSDnu12,CSDnu21,CSDG12,CSDthicknessStructure,CSDrhoStructure,							CSDintegrationScheme,CSDdelta,CSDswitchStiffening,CSDgenAlphaCoeff,CSDnewmarkGammaCoeff,CSDnewmarkBetaCoeff,CSDdampingCoeff1,CSDdampingCoeff2,							CSDflappingType,CSDplungingType,							CSDfrequency,CSDplungeAmplitudeX,CSDplungeAmplitudeY,CSDplungeAmplitudeZ,CSDflappingAmplitudeX,CSDflappingAmplitudeY,CSDflappingAmplitudeZ,							CFDIterationFinished{n,it-1},itfsi{n,it},CSDdimension") ;
#line 595 "FSI_CSDdataloop.loci"
       output("CSDdisplacementsStar{n,it}") ;
#line 595 "FSI_CSDdataloop.loci"
       output("CSDnodesDispStar{n,it}") ;
#line 595 "FSI_CSDdataloop.loci"
       output("CSDnodesVelStar{n,it}") ;
#line 595 "FSI_CSDdataloop.loci"
       output("CSDnodesAccStar{n,it}") ;
#line 595 "FSI_CSDdataloop.loci"
       output("CSDForcePreStar{n,it}") ;
#line 595 "FSI_CSDdataloop.loci"
       output("CSDnodesSysStar{n,it}") ;
#line 595 "FSI_CSDdataloop.loci"
       constraint("FSICoupling{n,it},FSINLAMS{n,it}") ;
#line 595 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 595 "FSI_CSDdataloop.loci"
    }
#line 595 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) { 
	// call NLAMS
	//communicateWithNLAMS(*CSDdisplacementsStar, ....) ;
	//if (*$CFDIterationFinished) {
	
	const int rank = Loci::MPI_rank ;
	const int CSDNumNodes = (*L_CSDnodes_ic_).size1() ;
	const int CSDNumElems = (*L_CSDConnectivity_).size1() ;
	const int CSDNumBC = (*L_CSDBCdof_).size() ;
	const int CSDAnswerSize = 6 * CSDNumNodes ;
		
	// Initialize
	(*L_CSDdisplacementsStar_nit__).resize(CSDNumNodes,3) ; std::fill( (*L_CSDdisplacementsStar_nit__).data().begin(),(*L_CSDdisplacementsStar_nit__).data().end(), 0. );
//	(*$CSDdisplacementsOldStar).resize(CSDNumNodes,3) ; (*$CSDdisplacementsOldStar).clear();
	(*L_CSDnodesDispStar_nit__).resize(CSDAnswerSize,1) ; std::fill( (*L_CSDnodesDispStar_nit__).data().begin(),(*L_CSDnodesDispStar_nit__).data().end(), 0. );
	(*L_CSDnodesVelStar_nit__).resize(CSDAnswerSize,1) ; std::fill( (*L_CSDnodesVelStar_nit__).data().begin(),(*L_CSDnodesVelStar_nit__).data().end(), 0. );
	(*L_CSDnodesAccStar_nit__).resize(CSDAnswerSize,1) ; std::fill( (*L_CSDnodesAccStar_nit__).data().begin(),(*L_CSDnodesAccStar_nit__).data().end(), 0. );
	(*L_CSDForcePreStar_nit__).resize(CSDAnswerSize,1) ; std::fill( (*L_CSDForcePreStar_nit__).data().begin(),(*L_CSDForcePreStar_nit__).data().end(), 0. );
	(*L_CSDnodesSysStar_nit__).resize(boost::extents[3][3][CSDNumElems][3]); std::fill( (*L_CSDnodesSysStar_nit__).origin(), (*L_CSDnodesSysStar_nit__).origin() + (*L_CSDnodesSysStar_nit__).size(), 0. );
		
	
//	if (*$CFDIterationFinished{n,it-1}) {	
		
	if (*L_ncycle_n__ < *L_CSDstartingTimeStep_) {
		 if (rank==0) cout << "CSDstartingTime not yet reached. Current time step = " << *L_ncycle_n__ << ", CSDstartingTimeStep = " << *L_CSDstartingTimeStep_<< endl ;
	} else {
		
	
	int tempTimeStepNumber = *L_ncycle_n__ - *L_CSDstartingTimeStep_+ 1 ;	
		
	if (rank==0) cout << "[I] Communicating with NLAMS.. " << endl ;
	if (rank==0) cout << "[I]... CSDNumNodes = " << CSDNumNodes << endl ;	
	if (rank==0) cout << "[I]... CSDNumElems = " << CSDNumElems << endl ;	
	if (rank==0) cout << "[I]... CSDNumBC = " << CSDNumBC << endl ;	
	if (rank==0) cout << "[I]... stime = " << *L_stime_n__ << endl ;	
	if (rank==0) cout << "[I]... current time step number [c++] = " << *L_ncycle_n__ << endl ;	
	if (rank==0) cout << "[I]... timestep = " << *L_timeStep_<< endl ;	
	if (rank==0) cout << "[I]... tip node = " << *L_CSDtipNode_<< ": " << (*L_CSDnodes_n__)(*L_CSDtipNode_,0) << ", "  << (*L_CSDnodes_n__)(*L_CSDtipNode_,1) << ", "<< (*L_CSDnodes_n__)(*L_CSDtipNode_,2) << endl ;	
	if (rank==0) cout << "[I]... CSDE1 = " << *L_CSDE1_<< endl ;					
	if (rank==0) cout << "[I]... CSDE2 = " << *L_CSDE2_<< endl ;		
	if (rank==0) cout << "[I]... CSDnu12 = " << *L_CSDnu12_<< endl ;		
	if (rank==0) cout << "[I]... CSDnu21 = " << *L_CSDnu21_<< endl ;		
	if (rank==0) cout << "[I]... CSDG12 = " << *L_CSDG12_<< endl ;		
	if (rank==0) cout << "[I]... CSDthicknessStructure = " << *L_CSDthicknessStructure_<< endl ;		
	if (rank==0) cout << "[I]... CSDrhoStructure = " << *L_CSDrhoStructure_<< endl ;		
	if (rank==0) cout << "[I]... CSDdelta = " << *L_CSDdelta_<< endl ;		
	//if (rank==0) cout << "[I]... CSDexcitationType = " << *$CSDexcitationType{n} << endl ;		
	if (rank==0) cout << "[I]... CSDfrequency = " << *L_CSDfrequency_<< endl ;		

	if (rank==0) {
		std::stringstream filename ;	
		filename << "CSDnodesforces" << tempTimeStepNumber << "it" << *L_itfsi_nit__ << ".dat" ;
		ofstream CSDnodesforces ;
		CSDnodesforces.open(filename.str().c_str(), ios::out) ;
			CSDnodesforces << "CSDnodes.x" << ", " << "CSDnodes.y" << ", " << "CSDnodes.z" << ", " << "CSDForce.x" << ", "  << "CSDForce.y" << ", " << "CSDForce.z" << "CSDdisp.x" << ", " << "CSDdisp.y" << ", " << "CSDdisp.z" << endl ;
		for(int i=0; i<CSDNumNodes;++i) { // in CFD coordinates
			CSDnodesforces << (*L_CSDnodes_n__)(i,0) << ", " << (*L_CSDnodes_n__)(i,1) << ", " << (*L_CSDnodes_n__)(i,2) << ", " << (*L_CSDForce_nit__)(i*6+0,0) << ", "  << (*L_CSDForce_nit__)(i*6+1,0) << ", " << (*L_CSDForce_nit__)(i*6+2,0)  << (*L_CSDnodes_n__)(i,0) - (*L_CSDnodes_ic_)(i,0)<< ", " << (*L_CSDnodes_n__)(i,1) - (*L_CSDnodes_ic_)(i,1) << ", " << (*L_CSDnodes_n__)(i,2) - (*L_CSDnodes_ic_)(i,2) << endl ;
		}
		CSDnodesforces.close();
	}
	
	// Displacements at the previous itfis
	//(CSDdisplacementsOld) = (*$CSDnodes) - (*$CSDnodes_ic) ;
	
//	if (rank==0) cout << "connectivity matrix in loci-stream: " << (*$CSDConnectivity) << endl ;

//	if (rank==0) cout << "force matrix in loci-stream: " << (*$CSDForce) << endl ;

  int itfsiInt = 0;

      for(int i=0; i<CSDNumNodes; ++i) {
		 if (rank==0) cout << "displacementsStarBefore: " << (*L_CSDdisplacementsStar_nit__)(i,0) << endl ;
  		 if (rank==0) cout << "displacementsStarBefore: " << (*L_CSDdisplacementsStar_nit__)(i,1) << endl ;
		 if (rank==0) cout << "displacementsStarBefore: " << (*L_CSDdisplacementsStar_nit__)(i,2) << endl ;
  	}		
	
  pass_(&rank, &itfsiInt, &CSDNumNodes, &(*L_CSDnodes_ic_)(0,0), &(*L_CSDnodes_ic_)(0,1),&(*L_CSDnodes_ic_)(0,2), 
  			&CSDNumElems, &(*L_CSDConnectivity_)(0,0), &CSDNumBC, &(*L_CSDBCdof_)[0], &(*L_CSDBCZeroConstraint_)[0], 
  			&(*L_CSDE1_), &(*L_CSDnu12_), &(*L_CSDrhoStructure_), &(*L_CSDthicknessStructure_), 
  			&(*L_CSDintegrationScheme_), &(*L_CSDnewmarkBetaCoeff_), &(*L_CSDnewmarkGammaCoeff_), &(*L_CSDgenAlphaCoeff_), 
  			&CSDAnswerSize, &(*L_CSDnodesDisp_n__)(0,0), &(*L_CSDnodesVel_n__)(0,0), &(*L_CSDnodesAcc_n__)(0,0), &(*L_CSDForce_nit__)(0,0), &(*L_CSDForcePre_n__)(0,0), (*L_CSDnodesSys_n__).data(),
  			&(*L_timeStep_), &tempTimeStepNumber, &(*L_CSDtipNode_), // &(*$ncycle), 
  			&(*L_CSDflappingType_), &(*L_CSDplungingType_), 
  			&(*L_CSDfrequency_), &(*L_CSDflappingAmplitudeX_), &(*L_CSDflappingAmplitudeY_), &(*L_CSDflappingAmplitudeZ_), &(*L_CSDplungeAmplitudeX_), &(*L_CSDplungeAmplitudeY_), &(*L_CSDplungeAmplitudeZ_), 
  			&(*L_CSDdisplacementsStar_nit__)(0,0), &(*L_CSDnodesDispStar_nit__)(0,0), &(*L_CSDnodesVelStar_nit__)(0,0), &(*L_CSDnodesAccStar_nit__)(0,0), (*L_CSDnodesSysStar_nit__).data() );
  
  if (rank==0) cout << "[I] Communication with NLAMS finished.. " << endl ;
  
  // Updates
  // CSDForce -> CSDForcePre
  (*L_CSDForcePreStar_nit__) = (*L_CSDForce_nit__ ) ;
  
  if (*L_CSDdimension_== 2) { // displacement from NLAMS is sometimes not exactly zero due to floating point errors
  	for(int i=0; i<CSDNumNodes; ++i) {
  		 (*L_CSDdisplacementsStar_nit__)(i,2) = 0. ;
  	}
  }
  //-----------
  for(int i=0; i<CSDNumNodes; ++i) {
//	  (*$CSDdisplacementsStar{n,it})(i,0) = 0.0 ;
//	  (*$CSDdisplacementsStar{n,it})(i,1) = 0.0 ;
	  (*L_CSDdisplacementsStar_nit__)(i,2) = 0.0 ;
	}
  //-----------
  
  			
	for (int i=0; i<3; ++i) {
		for (int n=0; n<CSDNumNodes; ++n) {
			if ( fabs((*L_CSDdisplacementsStar_nit__)(n,i))  < 1.0e-5 ) (*L_CSDdisplacementsStar_nit__)(n,i) = 0.0 ;
		}
	}
  
  
    for(int i=0; i<CSDNumNodes; ++i) {
		 if (rank==0) cout << "displacementsStar: " << (*L_CSDdisplacementsStar_nit__)(i,0) << endl ;
  		 if (rank==0) cout << "displacementsStar: " << (*L_CSDdisplacementsStar_nit__)(i,1) << endl ;
		 if (rank==0) cout << "displacementsStar: " << (*L_CSDdisplacementsStar_nit__)(i,2) << endl ;
  	}		
//  
//		const double PI = 4.0*atan(1.0) ;	
//  		real alpha_amp = 5.0 ; 
//		//real freq = *$CSDfrequency ; // period in sec
//  		//real alpha = alpha_amp * sin(2.*PI * freq * (*$stime{n})) ;  		
//		real tT = (*$ncycle{n}) / 500.0 ;
//		real alpha = alpha_amp * sin( 2.0 * PI * tT );
//
//  		for(int i=0;i<CSDNumNodes;++i) {
//		 	  (*$CSDdisplacementsStar{n,it})(i,0) = (*$CSDnodes_ic)(i,0) * ( cos(alpha*PI/180) - 1.) ;
//		 	  (*$CSDdisplacementsStar{n,it})(i,1) = (*$CSDnodes_ic)(i,0) * ( sin(alpha*PI/180) ) ;
//		 	  (*$CSDdisplacementsStar{n,it})(i,2) = 0.0 ;
//		  }
//	   	if (Loci::MPI_rank==0) cout << "CSDdisplacement:" << (*$stime{n}) << ", alpha = " << alpha << endl ;

  
 // (*$CSDdisplacementsOldStar) = (*$CSDForce ) ;
  
  
  
 } // end if CSDstartingTimeStep
 
//}

}    void compute(const Loci::sequence &seq) { 
#line 740 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 740 "FSI_CSDdataloop.loci"
    }
#line 740 "FSI_CSDdataloop.loci"
} ;
#line 740 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop036_1280807839m77> register_file_FSI_CSDdataloop036_1280807839m77 ;
#line 740 "FSI_CSDdataloop.loci"
}
#line 740 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop037_1280807839m84 : public Loci::singleton_rule {
#line 740 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerAxis_ ; 
#line 740 "FSI_CSDdataloop.loci"
    Loci::param<int>  L_CSDEulerBeamDirection_ ; 
#line 740 "FSI_CSDdataloop.loci"
    Loci::param<int>  L_CSDEulerSpanDirection_ ; 
#line 740 "FSI_CSDdataloop.loci"
public:
#line 740 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop037_1280807839m84() {
#line 740 "FSI_CSDdataloop.loci"
       name_store("CSDEulerAxis",L_CSDEulerAxis_) ;
#line 740 "FSI_CSDdataloop.loci"
       name_store("CSDEulerBeamDirection",L_CSDEulerBeamDirection_) ;
#line 740 "FSI_CSDdataloop.loci"
       name_store("CSDEulerSpanDirection",L_CSDEulerSpanDirection_) ;
#line 740 "FSI_CSDdataloop.loci"
       input("CSDEulerAxis") ;
#line 740 "FSI_CSDdataloop.loci"
       output("CSDEulerBeamDirection") ;
#line 740 "FSI_CSDdataloop.loci"
       output("CSDEulerSpanDirection") ;
#line 740 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 740 "FSI_CSDdataloop.loci"
    }
#line 740 "FSI_CSDdataloop.loci"
    void compute(const Loci::sequence &seq) { 
	(*L_CSDEulerBeamDirection_)= (*L_CSDEulerAxis_);
	(*L_CSDEulerSpanDirection_)= 3 - ((*L_CSDEulerAxis_)+ 1) ; // 3 % 1 [x] = 2 [z], 3 % 3 [z] = 0 [x]
}} ;
#line 743 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop037_1280807839m84> register_file_FSI_CSDdataloop037_1280807839m84 ;
#line 743 "FSI_CSDdataloop.loci"
}
#line 743 "FSI_CSDdataloop.loci"

   
namespace {class file_FSI_CSDdataloop038_1280807839m85 : public Loci::blackbox_rule {
#line 745 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDdimension_ ; 
#line 745 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSD2dSpanCenter_ ; 
#line 745 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDEulerXstart_ ; 
#line 745 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDEulerXend_ ; 
#line 745 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerXnum_ ; 
#line 745 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerBeamDirection_ ; 
#line 745 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerSpanDirection_ ; 
#line 745 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 745 "FSI_CSDdataloop.loci"
public:
#line 745 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop038_1280807839m85() {
#line 745 "FSI_CSDdataloop.loci"
       name_store("CSDdimension",L_CSDdimension_) ;
#line 745 "FSI_CSDdataloop.loci"
       name_store("CSD2dSpanCenter",L_CSD2dSpanCenter_) ;
#line 745 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 745 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXstart",L_CSDEulerXstart_) ;
#line 745 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXend",L_CSDEulerXend_) ;
#line 745 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXnum",L_CSDEulerXnum_) ;
#line 745 "FSI_CSDdataloop.loci"
       name_store("CSDEulerBeamDirection",L_CSDEulerBeamDirection_) ;
#line 745 "FSI_CSDdataloop.loci"
       name_store("CSDEulerSpanDirection",L_CSDEulerSpanDirection_) ;
#line 745 "FSI_CSDdataloop.loci"
       input("CSDEulerXstart,CSDEulerXend,CSDEulerXnum,CSD2dSpanCenter,CSDdimension,CSDEulerBeamDirection,CSDEulerSpanDirection") ;
#line 745 "FSI_CSDdataloop.loci"
       output("CSDnodes_ic") ;
#line 745 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 745 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 745 "FSI_CSDdataloop.loci"
    }
#line 745 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) { 
	int nElements = *L_CSDEulerXnum_- 1;
	real dx = (*L_CSDEulerXend_- *L_CSDEulerXstart_) / nElements ;

	real posz =0.0;
	posz = ((*L_CSDdimension_==2)? (*L_CSD2dSpanCenter_) : 0.0 ) ;

	int beamdirection = (*L_CSDEulerBeamDirection_) ;
	int spandirection = (*L_CSDEulerSpanDirection_) ;
	
	
	(*L_CSDnodes_ic_).resize(*L_CSDEulerXnum_,3) ; std::fill( (*L_CSDnodes_ic_).data().begin(), (*L_CSDnodes_ic_).data().end(), 0. );	    	
	
	(*L_CSDnodes_ic_)(0,0) = *L_CSDEulerXstart_;
	for(int i=0; i<*L_CSDEulerXnum_; ++i) {
		(*L_CSDnodes_ic_)(i,beamdirection) = *L_CSDEulerXstart_+ i * dx;
		(*L_CSDnodes_ic_)(i,spandirection) = posz;
	}
}    void compute(const Loci::sequence &seq) { 
#line 766 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 766 "FSI_CSDdataloop.loci"
    }
#line 766 "FSI_CSDdataloop.loci"
} ;
#line 766 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop038_1280807839m85> register_file_FSI_CSDdataloop038_1280807839m85 ;
#line 766 "FSI_CSDdataloop.loci"
}
#line 766 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop039_1280807839m86 : public Loci::blackbox_rule {
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDE1_ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDrhoStructure_ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDthicknessStructure_ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDplungingType_ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDfrequency_ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDplungeAmplitudeY_ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDflappingAmplitudeX_ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDstartingTimeStep_ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDEulerChord_ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerXnum_ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerBeamDirection_ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerSpanDirection_ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_param<bool>  L_CFDIterationFinished_nit_M_1__ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_stime_n__ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_ncycle_n__ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_timeStep_ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForce_nit__ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDx_n__ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxdot_n__ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxddot_n__ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDForcePreEuler_n__ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_it_nit__ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDdisplacementsStar_nit__ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxStar_nit__ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxdotStar_nit__ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxddotStar_nit__ ; 
#line 771 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDForcePreEulerStar_nit__ ; 
#line 771 "FSI_CSDdataloop.loci"
public:
#line 771 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop039_1280807839m86() {
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDE1",L_CSDE1_) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDrhoStructure",L_CSDrhoStructure_) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDthicknessStructure",L_CSDthicknessStructure_) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDplungingType",L_CSDplungingType_) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDfrequency",L_CSDfrequency_) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDplungeAmplitudeY",L_CSDplungeAmplitudeY_) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDflappingAmplitudeX",L_CSDflappingAmplitudeX_) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDstartingTimeStep",L_CSDstartingTimeStep_) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDEulerChord",L_CSDEulerChord_) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXnum",L_CSDEulerXnum_) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDEulerBeamDirection",L_CSDEulerBeamDirection_) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDEulerSpanDirection",L_CSDEulerSpanDirection_) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDdisplacementsStar{n,it}",L_CSDdisplacementsStar_nit__) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CFDIterationFinished{n,it-1}",L_CFDIterationFinished_nit_M_1__) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("stime{n}",L_stime_n__) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("ncycle{n}",L_ncycle_n__) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("timeStep",L_timeStep_) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDForce{n,it}",L_CSDForce_nit__) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDx{n}",L_CSDx_n__) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDxdot{n}",L_CSDxdot_n__) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDxddot{n}",L_CSDxddot_n__) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDForcePreEuler{n}",L_CSDForcePreEuler_n__) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("$it{n,it}",L_it_nit__) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDxStar{n,it}",L_CSDxStar_nit__) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDxdotStar{n,it}",L_CSDxdotStar_nit__) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDxddotStar{n,it}",L_CSDxddotStar_nit__) ;
#line 771 "FSI_CSDdataloop.loci"
       name_store("CSDForcePreEulerStar{n,it}",L_CSDForcePreEulerStar_nit__) ;
#line 771 "FSI_CSDdataloop.loci"
       input("							CSDEulerXnum,CSDnodes_ic,CSDx{n},CSDxdot{n},CSDxddot{n},CSDForce{n,it},CSDForcePreEuler{n},							stime{n},ncycle{n},timeStep,CSDstartingTimeStep,CSDEulerBeamDirection,CSDEulerSpanDirection,$it{n,it},							CSDE1,CSDthicknessStructure,CSDrhoStructure,CSDEulerChord,														CSDfrequency,CSDplungeAmplitudeY,CSDflappingAmplitudeX,CSDplungingType,CFDIterationFinished{n,it-1}") ;
#line 771 "FSI_CSDdataloop.loci"
       output("CSDdisplacementsStar{n,it}") ;
#line 771 "FSI_CSDdataloop.loci"
       output("CSDxStar{n,it}") ;
#line 771 "FSI_CSDdataloop.loci"
       output("CSDxdotStar{n,it}") ;
#line 771 "FSI_CSDdataloop.loci"
       output("CSDxddotStar{n,it}") ;
#line 771 "FSI_CSDdataloop.loci"
       output("CSDForcePreEulerStar{n,it}") ;
#line 771 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 771 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 771 "FSI_CSDdataloop.loci"
    }
#line 771 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) { 
	// call NLAMS
	//communicateWithNLAMS(*CSDdisplacementsStar, ....) ;
  //  if (*$CFDIterationFinished{n,it-1}) {
	
	const int rank = Loci::MPI_rank ;

	if (rank==0) cout << "[I] Communicating with the Euler-Bernoulli beam solver.. " << endl ;
	// Communicate with NLAMS
	int CSDNumNodes = (*L_CSDEulerXnum_) ;
	int CSDNumBC = 2 ;
	int CSDdofNode = 2 ;
	int CSDdofFree = CSDNumNodes * CSDdofNode - CSDNumBC ;
	int CSDAnswerSize = 6 * CSDNumNodes ;

	int beamdirection = (*L_CSDEulerBeamDirection_) ;
	int spandirection = (*L_CSDEulerSpanDirection_) ;
	
	int itcurrent = *L_it_nit__ ;
	
	// Initialize
	(*L_CSDdisplacementsStar_nit__).resize(CSDNumNodes,3) ; std::fill( (*L_CSDdisplacementsStar_nit__).data().begin(),(*L_CSDdisplacementsStar_nit__).data().end(), 0. );	
	(*L_CSDxStar_nit__).resize(CSDdofFree) ; std::fill( (*L_CSDxStar_nit__).begin(),(*L_CSDxStar_nit__).end(), 0. );
	(*L_CSDxdotStar_nit__).resize(CSDdofFree) ; std::fill( (*L_CSDxdotStar_nit__).begin(),(*L_CSDxdotStar_nit__).end(), 0. );
	(*L_CSDxddotStar_nit__).resize(CSDdofFree) ; std::fill( (*L_CSDxddotStar_nit__).begin(),(*L_CSDxddotStar_nit__).end(), 0. );
	(*L_CSDForcePreEulerStar_nit__).resize(CSDdofFree) ; std::fill( (*L_CSDForcePreEulerStar_nit__).begin(),(*L_CSDForcePreEulerStar_nit__).end(), 0. );


	if (rank==0) cout << "[I] Euler-Bernoulli.. initialized, timestep, itcurrent = " << *L_ncycle_n__ << ", " << *L_it_nit__ << endl ;
	
	if (rank==0) cout << "[I] *CFDIterationFinished{n,it-1} = " << (*L_CFDIterationFinished_nit_M_1__) << endl ;

	if (*L_CFDIterationFinished_nit_M_1__) {	

		if (*L_ncycle_n__ < *L_CSDstartingTimeStep_) {
			 if (rank==0) cout << "CSDstartingTime not yet reached. Current time step = " << *L_ncycle_n__ << ", CSDstartingTimeStep = " << *L_CSDstartingTimeStep_<< endl ;
		} else {	
		  
		  if (itcurrent > 1) { // bypass the first iterations
			int tempTimeStepNumber = *L_ncycle_n__ - *L_CSDstartingTimeStep_+ 1 ;	
				
			ublas::vector<real> CSDbeam1dForce(CSDNumNodes) ;	std::fill(CSDbeam1dForce.begin(), CSDbeam1dForce.end(), 0.0) ;
			ublas::vector<real> disp(CSDNumNodes) ; std::fill(disp.begin(), disp.end(), 0.0) ;
			ublas::vector<real> nodes(CSDNumNodes) ; std::fill(nodes.begin(), nodes.end(), 0.0) ;
			
			if (rank==0) cout << "CSDstartingTime reached. Current time step = " << *L_ncycle_n__ << ", CSDstartingTimeStep = " << *L_CSDstartingTimeStep_<< ", tempTimeStepNumber = " << tempTimeStepNumber << endl ;
			if (rank==0) cout << "[I] Euler-Bernoulli.. communication starting: " << CSDNumNodes << ", " << CSDbeam1dForce.size() << ", " << (*L_CSDForce_nit__).size1() << endl ;	

			for (int i=0; i<CSDNumNodes; ++i) {
	//			if (rank==0) cout << "[I] Euler-Bernoulli.. before CSDforce" << endl ;	
				CSDbeam1dForce[i] = (*L_CSDForce_nit__)(6*i+1,0) ; // y-component only
	//			CSDbeam1dForce[i] = 0.0 ; // y-component only			
	//			if (rank==0) cout << "[I] Euler-Bernoulli.. CSDbeam1dForce[" << i << "] = " << CSDbeam1dForce[i] << endl ;	
				nodes[i] = (*L_CSDnodes_ic_)(i,beamdirection) ;
	//			if (rank==0) cout << "[I] Euler-Bernoulli.. nodes[" << i << "] = " << nodes[i] << endl ;	
			}
			
//			for (int i=0; i<CSDdofFree; ++i) {
//				cout << "rank = " << rank << ", " << "[I] Euler-Bernoulli: i, x, xdot, xddot, disp = " << i << "\t" << (*$CSDx{n})[i] << "\t" << (*$CSDxdot{n})[i] << "\t" << (*$CSDxddot{n})[i] << "\t" << endl ;
//			}
			
		if (rank==0 && itcurrent < 5) {
			std::stringstream filename ;	
			//filename << "CSDnodesforces" << setfill('0') << setw(5) << tempTimeStepNumber << "it" << itcurrent << ".dat" ;
			filename << "CSDnodesforces" << setfill('0') << setw(5) << tempTimeStepNumber << ".dat" ;
			ofstream CSDnodesforces ;
			CSDnodesforces.open(filename.str().c_str(), ios::out) ;
			CSDnodesforces << "CSDnodes.x" << ", " << "CSDnodes.y" << ", " << "CSDnodes.z" << ", " << "CSDForce.x" << ", "  << "CSDForce.y" << ", " << "CSDForce.z" << "CSDdisp.x" << ", " << "CSDdisp.y" << ", " << "CSDdisp.z" << endl ;
			for(int i=0; i<CSDNumNodes;++i) { // in CFD coordinates
				CSDnodesforces << nodes[i] << ", " << (*L_CSDForce_nit__)(i*6+0,0) << ", "  << (*L_CSDForce_nit__)(i*6+1,0) << ", " << (*L_CSDForce_nit__)(i*6+2,0)  << endl ;
			}
			CSDnodesforces.close();
		}
	

			if (rank==0) cout << "[I] Euler-Bernoulli.. communication starting " << endl ;
			
			real currentTime = ( (tempTimeStepNumber - 1) * (*L_timeStep_)) ;

			// communiate	
	     beam1d_(&rank, &CSDNumNodes, &CSDdofFree, &(nodes)[0], &(*L_CSDplungingType_), &(*L_CSDplungeAmplitudeY_), &(*L_CSDflappingAmplitudeX_), &(*L_CSDfrequency_), &(*L_CSDEulerChord_), &(*L_CSDrhoStructure_), &(*L_CSDthicknessStructure_),
	        &*L_CSDE1_, &(*L_CSDx_n__)[0], &(*L_CSDxdot_n__)[0], &(*L_CSDxddot_n__)[0], &(*L_CSDForcePreEuler_n__)[0], &(CSDbeam1dForce)[0], &(*L_timeStep_), &(currentTime), &tempTimeStepNumber, 
	        &(*L_CSDxStar_nit__)[0], &(*L_CSDxdotStar_nit__)[0], &(*L_CSDxddotStar_nit__)[0], &(*L_CSDForcePreEulerStar_nit__)[0], &(disp)[0]) ;

	  
	  
		  // update
		for (int i=0; i<CSDNumNodes; ++i) {
				//(*$CSDdisplacementsStar{n,it})(i,0) = 0.0 ;
				(*L_CSDdisplacementsStar_nit__)(i,1) = disp[i];
				//(*$CSDdisplacementsStar{n,it})(i,2) = 0.0 ;
	//			cout << "rank = " << rank << ", " << "[I] Euler-Bernoulli: i, CSDdisplacementsStar, disp = " << i << "\t" << (*$CSDdisplacementsStar{n,it})(i,1) << "\t" << disp[i] << endl ;				
		}
		
	//	(*$CSDForcePreStar{n,it}) = (*$CSDForce{n,it}) ;
		
		for (int i=0; i<3; ++i) {
			for (int n=0; n<CSDNumNodes; ++n) {
//				if ( fabs((*$CSDdisplacementsStar{n,it})(n,i))  < 1.0e-5 ) (*$CSDdisplacementsStar{n,it})(n,i) = 0.0 ;
			}
		}
		
//		const double PI = 4.0*atan(1.0) ;	
//  		real alpha_amp = 15.0 ; 
//		//real freq = *$CSDfrequency ; // period in sec
//  		//real alpha = alpha_amp * sin(2.*PI * freq * (*$stime{n})) ;  		
//		real tT = (*$ncycle{n}) / 500.0 ;
//		real alpha = alpha_amp * sin( 2.0 * PI * tT );
//
//  		for(int i=0;i<CSDNumNodes;++i) {
//		 	  (*$CSDdisplacementsStar{n,it})(i,0) = (*$CSDnodes_ic)(i,0) * ( cos(alpha*PI/180) - 1.) ;
////		 	  (*$CSDdisplacementsStar{n,it})(i,0) = 0.0 ;
//		 	  (*$CSDdisplacementsStar{n,it})(i,1) = (*$CSDnodes_ic)(i,0) * ( sin(alpha*PI/180) ) ;
//		 	  (*$CSDdisplacementsStar{n,it})(i,2) = 0.0 ;
//		  }
//	   	if (Loci::MPI_rank==0) cout << "CSDdisplacement:" << (*$stime{n}) << ", alpha = " << alpha << endl ;			
//			
	//	for (int i=0; i<CSDdofFree; ++i) {
	//		if (rank==0) cout << "[I] Euler-Bernoulli: i, xStar, xdotStar, xddotStar = " << i << "\t" << (*$CSDxStar{n})[i] << "\t" << (*$CSDxdotStar{n})[i] << "\t" << (*$CSDxddotStar{n})[i] << "\t" << endl ;
	//	}
			
		if (rank==0) cout << "[I] Communication with the Euler-Bernoulli beam solver done. " << endl ;
	  } else { if (rank==0) cout << "[I] Communication with the Euler-Bernoulli beam solver bypass because the current iteration number is = " << itcurrent << endl ; }
	}   
//	} else {
//		if (rank==0) cout << "[I] Communication with the Euler-Bernoulli beam solver bypassed. " << endl ;
    } // CFDIterationFinished
}    void compute(const Loci::sequence &seq) { 
#line 903 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 903 "FSI_CSDdataloop.loci"
    }
#line 903 "FSI_CSDdataloop.loci"
} ;
#line 903 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop039_1280807839m86> register_file_FSI_CSDdataloop039_1280807839m86 ;
#line 903 "FSI_CSDdataloop.loci"
}
#line 903 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop040_1280807839m90 : public Loci::blackbox_rule {
#line 903 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerXnum_ ; 
#line 903 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDx_ic_ ; 
#line 903 "FSI_CSDdataloop.loci"
public:
#line 903 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop040_1280807839m90() {
#line 903 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXnum",L_CSDEulerXnum_) ;
#line 903 "FSI_CSDdataloop.loci"
       name_store("CSDx_ic",L_CSDx_ic_) ;
#line 903 "FSI_CSDdataloop.loci"
       input("CSDEulerXnum") ;
#line 903 "FSI_CSDdataloop.loci"
       output("CSDx_ic") ;
#line 903 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 903 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 903 "FSI_CSDdataloop.loci"
    }
#line 903 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) { 
	int CSDNumBC = 2 ;
	int CSDdofNode = 2 ;
	int CSDdofFree = *L_CSDEulerXnum_* CSDdofNode - CSDNumBC ;
	
	(*L_CSDx_ic_).resize(CSDdofFree) ; std::fill( (*L_CSDx_ic_).begin(),(*L_CSDx_ic_).end(), 0. );
}    void compute(const Loci::sequence &seq) { 
#line 911 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 911 "FSI_CSDdataloop.loci"
    }
#line 911 "FSI_CSDdataloop.loci"
} ;
#line 911 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop040_1280807839m90> register_file_FSI_CSDdataloop040_1280807839m90 ;
#line 911 "FSI_CSDdataloop.loci"
}
#line 911 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop041_1280807839m90 : public Loci::blackbox_rule {
#line 911 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDx_ic_ ; 
#line 911 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDx_n_EQ_0__ ; 
#line 911 "FSI_CSDdataloop.loci"
public:
#line 911 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop041_1280807839m90() {
#line 911 "FSI_CSDdataloop.loci"
       name_store("CSDx_ic",L_CSDx_ic_) ;
#line 911 "FSI_CSDdataloop.loci"
       name_store("CSDx{n=0}",L_CSDx_n_EQ_0__) ;
#line 911 "FSI_CSDdataloop.loci"
       input("CSDx_ic") ;
#line 911 "FSI_CSDdataloop.loci"
       output("CSDx{n=0}") ;
#line 911 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 911 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 911 "FSI_CSDdataloop.loci"
    }
#line 911 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDx_n_EQ_0__).resize((*L_CSDx_ic_).size()) ; std::fill( (*L_CSDx_n_EQ_0__).begin(),(*L_CSDx_n_EQ_0__).end(), 0. );
	(*L_CSDx_n_EQ_0__) = (*L_CSDx_ic_) ;
}    void compute(const Loci::sequence &seq) { 
#line 918 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 918 "FSI_CSDdataloop.loci"
    }
#line 918 "FSI_CSDdataloop.loci"
} ;
#line 918 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop041_1280807839m90> register_file_FSI_CSDdataloop041_1280807839m90 ;
#line 918 "FSI_CSDdataloop.loci"
}
#line 918 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop042_1280807839m91 : public Loci::blackbox_rule {
#line 918 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDx_n__ ; 
#line 918 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDx_nit_EQ_0__ ; 
#line 918 "FSI_CSDdataloop.loci"
public:
#line 918 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop042_1280807839m91() {
#line 918 "FSI_CSDdataloop.loci"
       name_store("CSDx{n}",L_CSDx_n__) ;
#line 918 "FSI_CSDdataloop.loci"
       name_store("CSDx{n,it=0}",L_CSDx_nit_EQ_0__) ;
#line 918 "FSI_CSDdataloop.loci"
       input("CSDx{n}") ;
#line 918 "FSI_CSDdataloop.loci"
       output("CSDx{n,it=0}") ;
#line 918 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 918 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 918 "FSI_CSDdataloop.loci"
    }
#line 918 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDx_nit_EQ_0__).resize((*L_CSDx_n__).size()) ;std::fill( (*L_CSDx_nit_EQ_0__).begin(),(*L_CSDx_nit_EQ_0__).end(), 0. );
	(*L_CSDx_nit_EQ_0__) = (*L_CSDx_n__) ;
	//if (Loci::MPI_rank==0) cout << "Inside CSDnodesAcc{n,it=0}" << endl;
}    void compute(const Loci::sequence &seq) { 
#line 925 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 925 "FSI_CSDdataloop.loci"
    }
#line 925 "FSI_CSDdataloop.loci"
} ;
#line 925 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop042_1280807839m91> register_file_FSI_CSDdataloop042_1280807839m91 ;
#line 925 "FSI_CSDdataloop.loci"
}
#line 925 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop043_1280807839m91 : public Loci::blackbox_rule {
#line 925 "FSI_CSDdataloop.loci"
    Loci::const_param<bool>  L_CFDIterationFinished_nit_M_1__ ; 
#line 925 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxStar_nit__ ; 
#line 925 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDx_nit__ ; 
#line 925 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDx_nit_P_1__ ; 
#line 925 "FSI_CSDdataloop.loci"
public:
#line 925 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop043_1280807839m91() {
#line 925 "FSI_CSDdataloop.loci"
       name_store("CFDIterationFinished{n,it-1}",L_CFDIterationFinished_nit_M_1__) ;
#line 925 "FSI_CSDdataloop.loci"
       name_store("CSDxStar{n,it}",L_CSDxStar_nit__) ;
#line 925 "FSI_CSDdataloop.loci"
       name_store("CSDx{n,it}",L_CSDx_nit__) ;
#line 925 "FSI_CSDdataloop.loci"
       name_store("CSDx{n,it+1}",L_CSDx_nit_P_1__) ;
#line 925 "FSI_CSDdataloop.loci"
       input("CSDx{n,it},CSDxStar{n,it},CFDIterationFinished{n,it-1}") ;
#line 925 "FSI_CSDdataloop.loci"
       output("CSDx{n,it+1}") ;
#line 925 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 925 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 925 "FSI_CSDdataloop.loci"
    }
#line 925 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDx_nit_P_1__).resize((*L_CSDx_nit__).size()) ;std::fill( (*L_CSDx_nit_P_1__).begin(),(*L_CSDx_nit_P_1__).end(), 0. );
	if (*L_CFDIterationFinished_nit_M_1__) {
		(*L_CSDx_nit_P_1__) = (*L_CSDxStar_nit__) ;
	} else {
		(*L_CSDx_nit_P_1__) = (*L_CSDx_nit__) ;
	}
}    void compute(const Loci::sequence &seq) { 
#line 935 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 935 "FSI_CSDdataloop.loci"
    }
#line 935 "FSI_CSDdataloop.loci"
} ;
#line 935 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop043_1280807839m91> register_file_FSI_CSDdataloop043_1280807839m91 ;
#line 935 "FSI_CSDdataloop.loci"
}
#line 935 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop044_1280807839m92 : public Loci::blackbox_rule {
#line 935 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDx_nit__ ; 
#line 935 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDx_n_P_1__ ; 
#line 935 "FSI_CSDdataloop.loci"
public:
#line 935 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop044_1280807839m92() {
#line 935 "FSI_CSDdataloop.loci"
       name_store("CSDx{n,it}",L_CSDx_nit__) ;
#line 935 "FSI_CSDdataloop.loci"
       name_store("CSDx{n+1}",L_CSDx_n_P_1__) ;
#line 935 "FSI_CSDdataloop.loci"
       input("CSDx{n,it}") ;
#line 935 "FSI_CSDdataloop.loci"
       output("CSDx{n+1}") ;
#line 935 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 935 "FSI_CSDdataloop.loci"
       conditional("iterationFinished{n,it-1}") ;
#line 935 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 935 "FSI_CSDdataloop.loci"
    }
#line 935 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDx_n_P_1__).resize((*L_CSDx_nit__).size()) ; std::fill( (*L_CSDx_n_P_1__).begin(),(*L_CSDx_n_P_1__).end(), 0. );
	(*L_CSDx_n_P_1__) = (*L_CSDx_nit__) ;
}    void compute(const Loci::sequence &seq) { 
#line 947 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 947 "FSI_CSDdataloop.loci"
    }
#line 947 "FSI_CSDdataloop.loci"
} ;
#line 947 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop044_1280807839m92> register_file_FSI_CSDdataloop044_1280807839m92 ;
#line 947 "FSI_CSDdataloop.loci"
}
#line 947 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop045_1280807839m92 : public Loci::blackbox_rule {
#line 947 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerXnum_ ; 
#line 947 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxdot_ic_ ; 
#line 947 "FSI_CSDdataloop.loci"
public:
#line 947 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop045_1280807839m92() {
#line 947 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXnum",L_CSDEulerXnum_) ;
#line 947 "FSI_CSDdataloop.loci"
       name_store("CSDxdot_ic",L_CSDxdot_ic_) ;
#line 947 "FSI_CSDdataloop.loci"
       input("CSDEulerXnum") ;
#line 947 "FSI_CSDdataloop.loci"
       output("CSDxdot_ic") ;
#line 947 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 947 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 947 "FSI_CSDdataloop.loci"
    }
#line 947 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	
	int CSDNumBC = 2 ;
	int CSDdofNode = 2 ;
	int CSDdofFree = *L_CSDEulerXnum_* CSDdofNode - CSDNumBC ;
	
	(*L_CSDxdot_ic_).resize(CSDdofFree) ; std::fill( (*L_CSDxdot_ic_).begin(),(*L_CSDxdot_ic_).end(), 0. );
}    void compute(const Loci::sequence &seq) { 
#line 956 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 956 "FSI_CSDdataloop.loci"
    }
#line 956 "FSI_CSDdataloop.loci"
} ;
#line 956 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop045_1280807839m92> register_file_FSI_CSDdataloop045_1280807839m92 ;
#line 956 "FSI_CSDdataloop.loci"
}
#line 956 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop046_1280807839m93 : public Loci::blackbox_rule {
#line 956 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxdot_ic_ ; 
#line 956 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxdot_n_EQ_0__ ; 
#line 956 "FSI_CSDdataloop.loci"
public:
#line 956 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop046_1280807839m93() {
#line 956 "FSI_CSDdataloop.loci"
       name_store("CSDxdot_ic",L_CSDxdot_ic_) ;
#line 956 "FSI_CSDdataloop.loci"
       name_store("CSDxdot{n=0}",L_CSDxdot_n_EQ_0__) ;
#line 956 "FSI_CSDdataloop.loci"
       input("CSDxdot_ic") ;
#line 956 "FSI_CSDdataloop.loci"
       output("CSDxdot{n=0}") ;
#line 956 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 956 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 956 "FSI_CSDdataloop.loci"
    }
#line 956 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDxdot_n_EQ_0__).resize((*L_CSDxdot_ic_).size()) ; std::fill( (*L_CSDxdot_n_EQ_0__).begin(),(*L_CSDxdot_n_EQ_0__).end(), 0. );
	(*L_CSDxdot_n_EQ_0__) = (*L_CSDxdot_ic_) ;
}    void compute(const Loci::sequence &seq) { 
#line 962 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 962 "FSI_CSDdataloop.loci"
    }
#line 962 "FSI_CSDdataloop.loci"
} ;
#line 962 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop046_1280807839m93> register_file_FSI_CSDdataloop046_1280807839m93 ;
#line 962 "FSI_CSDdataloop.loci"
}
#line 962 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop047_1280807839m93 : public Loci::blackbox_rule {
#line 962 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxdot_n__ ; 
#line 962 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxdot_nit_EQ_0__ ; 
#line 962 "FSI_CSDdataloop.loci"
public:
#line 962 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop047_1280807839m93() {
#line 962 "FSI_CSDdataloop.loci"
       name_store("CSDxdot{n}",L_CSDxdot_n__) ;
#line 962 "FSI_CSDdataloop.loci"
       name_store("CSDxdot{n,it=0}",L_CSDxdot_nit_EQ_0__) ;
#line 962 "FSI_CSDdataloop.loci"
       input("CSDxdot{n}") ;
#line 962 "FSI_CSDdataloop.loci"
       output("CSDxdot{n,it=0}") ;
#line 962 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 962 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 962 "FSI_CSDdataloop.loci"
    }
#line 962 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDxdot_nit_EQ_0__).resize((*L_CSDxdot_n__).size()) ;std::fill( (*L_CSDxdot_nit_EQ_0__).begin(),(*L_CSDxdot_nit_EQ_0__).end(), 0. );
	(*L_CSDxdot_nit_EQ_0__) = (*L_CSDxdot_n__) ;
	//if (Loci::MPI_rank==0) cout << "Inside CSDnodesAcc{n,it=0}" << endl;
}    void compute(const Loci::sequence &seq) { 
#line 969 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 969 "FSI_CSDdataloop.loci"
    }
#line 969 "FSI_CSDdataloop.loci"
} ;
#line 969 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop047_1280807839m93> register_file_FSI_CSDdataloop047_1280807839m93 ;
#line 969 "FSI_CSDdataloop.loci"
}
#line 969 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop048_1280807839m94 : public Loci::blackbox_rule {
#line 969 "FSI_CSDdataloop.loci"
    Loci::const_param<bool>  L_CFDIterationFinished_nit_M_1__ ; 
#line 969 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxdotStar_nit__ ; 
#line 969 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxdot_nit__ ; 
#line 969 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxdot_nit_P_1__ ; 
#line 969 "FSI_CSDdataloop.loci"
public:
#line 969 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop048_1280807839m94() {
#line 969 "FSI_CSDdataloop.loci"
       name_store("CFDIterationFinished{n,it-1}",L_CFDIterationFinished_nit_M_1__) ;
#line 969 "FSI_CSDdataloop.loci"
       name_store("CSDxdotStar{n,it}",L_CSDxdotStar_nit__) ;
#line 969 "FSI_CSDdataloop.loci"
       name_store("CSDxdot{n,it}",L_CSDxdot_nit__) ;
#line 969 "FSI_CSDdataloop.loci"
       name_store("CSDxdot{n,it+1}",L_CSDxdot_nit_P_1__) ;
#line 969 "FSI_CSDdataloop.loci"
       input("CSDxdot{n,it},CSDxdotStar{n,it},CFDIterationFinished{n,it-1}") ;
#line 969 "FSI_CSDdataloop.loci"
       output("CSDxdot{n,it+1}") ;
#line 969 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 969 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 969 "FSI_CSDdataloop.loci"
    }
#line 969 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDxdot_nit_P_1__).resize((*L_CSDxdot_nit__).size()) ;std::fill( (*L_CSDxdot_nit_P_1__).begin(),(*L_CSDxdot_nit_P_1__).end(), 0. );
	if (*L_CFDIterationFinished_nit_M_1__) {
		(*L_CSDxdot_nit_P_1__) = (*L_CSDxdotStar_nit__) ;
	} else {
		(*L_CSDxdot_nit_P_1__) = (*L_CSDxdot_nit__) ;
	}
}    void compute(const Loci::sequence &seq) { 
#line 984 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 984 "FSI_CSDdataloop.loci"
    }
#line 984 "FSI_CSDdataloop.loci"
} ;
#line 984 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop048_1280807839m94> register_file_FSI_CSDdataloop048_1280807839m94 ;
#line 984 "FSI_CSDdataloop.loci"
}
#line 984 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop049_1280807839m94 : public Loci::blackbox_rule {
#line 984 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxdot_nit__ ; 
#line 984 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxdot_n_P_1__ ; 
#line 984 "FSI_CSDdataloop.loci"
public:
#line 984 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop049_1280807839m94() {
#line 984 "FSI_CSDdataloop.loci"
       name_store("CSDxdot{n,it}",L_CSDxdot_nit__) ;
#line 984 "FSI_CSDdataloop.loci"
       name_store("CSDxdot{n+1}",L_CSDxdot_n_P_1__) ;
#line 984 "FSI_CSDdataloop.loci"
       input("CSDxdot{n,it}") ;
#line 984 "FSI_CSDdataloop.loci"
       output("CSDxdot{n+1}") ;
#line 984 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 984 "FSI_CSDdataloop.loci"
       conditional("iterationFinished{n,it-1}") ;
#line 984 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 984 "FSI_CSDdataloop.loci"
    }
#line 984 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDxdot_n_P_1__).resize((*L_CSDxdot_nit__).size()) ; std::fill( (*L_CSDxdot_n_P_1__).begin(),(*L_CSDxdot_n_P_1__).end(), 0. );
	(*L_CSDxdot_n_P_1__) = (*L_CSDxdot_nit__) ;
}    void compute(const Loci::sequence &seq) { 
#line 996 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 996 "FSI_CSDdataloop.loci"
    }
#line 996 "FSI_CSDdataloop.loci"
} ;
#line 996 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop049_1280807839m94> register_file_FSI_CSDdataloop049_1280807839m94 ;
#line 996 "FSI_CSDdataloop.loci"
}
#line 996 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop050_1280807839m95 : public Loci::blackbox_rule {
#line 996 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerXnum_ ; 
#line 996 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxddot_ic_ ; 
#line 996 "FSI_CSDdataloop.loci"
public:
#line 996 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop050_1280807839m95() {
#line 996 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXnum",L_CSDEulerXnum_) ;
#line 996 "FSI_CSDdataloop.loci"
       name_store("CSDxddot_ic",L_CSDxddot_ic_) ;
#line 996 "FSI_CSDdataloop.loci"
       input("CSDEulerXnum") ;
#line 996 "FSI_CSDdataloop.loci"
       output("CSDxddot_ic") ;
#line 996 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 996 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 996 "FSI_CSDdataloop.loci"
    }
#line 996 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	
	int CSDNumBC = 2 ;
	int CSDdofNode = 2 ;
	int CSDdofFree = *L_CSDEulerXnum_* CSDdofNode - CSDNumBC ;
	
	(*L_CSDxddot_ic_).resize(CSDdofFree) ; std::fill( (*L_CSDxddot_ic_).begin(),(*L_CSDxddot_ic_).end(), 0. );
}    void compute(const Loci::sequence &seq) { 
#line 1005 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 1005 "FSI_CSDdataloop.loci"
    }
#line 1005 "FSI_CSDdataloop.loci"
} ;
#line 1005 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop050_1280807839m95> register_file_FSI_CSDdataloop050_1280807839m95 ;
#line 1005 "FSI_CSDdataloop.loci"
}
#line 1005 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop051_1280807839m95 : public Loci::blackbox_rule {
#line 1005 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxddot_ic_ ; 
#line 1005 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxddot_n_EQ_0__ ; 
#line 1005 "FSI_CSDdataloop.loci"
public:
#line 1005 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop051_1280807839m95() {
#line 1005 "FSI_CSDdataloop.loci"
       name_store("CSDxddot_ic",L_CSDxddot_ic_) ;
#line 1005 "FSI_CSDdataloop.loci"
       name_store("CSDxddot{n=0}",L_CSDxddot_n_EQ_0__) ;
#line 1005 "FSI_CSDdataloop.loci"
       input("CSDxddot_ic") ;
#line 1005 "FSI_CSDdataloop.loci"
       output("CSDxddot{n=0}") ;
#line 1005 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 1005 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 1005 "FSI_CSDdataloop.loci"
    }
#line 1005 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDxddot_n_EQ_0__).resize((*L_CSDxddot_ic_).size()) ; std::fill( (*L_CSDxddot_n_EQ_0__).begin(),(*L_CSDxddot_n_EQ_0__).end(), 0. );
	(*L_CSDxddot_n_EQ_0__) = (*L_CSDxddot_ic_) ;
}    void compute(const Loci::sequence &seq) { 
#line 1011 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 1011 "FSI_CSDdataloop.loci"
    }
#line 1011 "FSI_CSDdataloop.loci"
} ;
#line 1011 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop051_1280807839m95> register_file_FSI_CSDdataloop051_1280807839m95 ;
#line 1011 "FSI_CSDdataloop.loci"
}
#line 1011 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop052_1280807839m95 : public Loci::blackbox_rule {
#line 1011 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxddot_n__ ; 
#line 1011 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxddot_nit_EQ_0__ ; 
#line 1011 "FSI_CSDdataloop.loci"
public:
#line 1011 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop052_1280807839m95() {
#line 1011 "FSI_CSDdataloop.loci"
       name_store("CSDxddot{n}",L_CSDxddot_n__) ;
#line 1011 "FSI_CSDdataloop.loci"
       name_store("CSDxddot{n,it=0}",L_CSDxddot_nit_EQ_0__) ;
#line 1011 "FSI_CSDdataloop.loci"
       input("CSDxddot{n}") ;
#line 1011 "FSI_CSDdataloop.loci"
       output("CSDxddot{n,it=0}") ;
#line 1011 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 1011 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 1011 "FSI_CSDdataloop.loci"
    }
#line 1011 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDxddot_nit_EQ_0__).resize((*L_CSDxddot_n__).size()) ;std::fill( (*L_CSDxddot_nit_EQ_0__).begin(),(*L_CSDxddot_nit_EQ_0__).end(), 0. );
	(*L_CSDxddot_nit_EQ_0__) = (*L_CSDxddot_n__) ;
	//if (Loci::MPI_rank==0) cout << "Inside CSDnodesAcc{n,it=0}" << endl;
}    void compute(const Loci::sequence &seq) { 
#line 1018 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 1018 "FSI_CSDdataloop.loci"
    }
#line 1018 "FSI_CSDdataloop.loci"
} ;
#line 1018 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop052_1280807839m95> register_file_FSI_CSDdataloop052_1280807839m95 ;
#line 1018 "FSI_CSDdataloop.loci"
}
#line 1018 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop053_1280807839m96 : public Loci::blackbox_rule {
#line 1018 "FSI_CSDdataloop.loci"
    Loci::const_param<bool>  L_CFDIterationFinished_nit_M_1__ ; 
#line 1018 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxddotStar_nit__ ; 
#line 1018 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxddot_nit__ ; 
#line 1018 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxddot_nit_P_1__ ; 
#line 1018 "FSI_CSDdataloop.loci"
public:
#line 1018 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop053_1280807839m96() {
#line 1018 "FSI_CSDdataloop.loci"
       name_store("CFDIterationFinished{n,it-1}",L_CFDIterationFinished_nit_M_1__) ;
#line 1018 "FSI_CSDdataloop.loci"
       name_store("CSDxddotStar{n,it}",L_CSDxddotStar_nit__) ;
#line 1018 "FSI_CSDdataloop.loci"
       name_store("CSDxddot{n,it}",L_CSDxddot_nit__) ;
#line 1018 "FSI_CSDdataloop.loci"
       name_store("CSDxddot{n,it+1}",L_CSDxddot_nit_P_1__) ;
#line 1018 "FSI_CSDdataloop.loci"
       input("CSDxddot{n,it},CSDxddotStar{n,it},CFDIterationFinished{n,it-1}") ;
#line 1018 "FSI_CSDdataloop.loci"
       output("CSDxddot{n,it+1}") ;
#line 1018 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 1018 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 1018 "FSI_CSDdataloop.loci"
    }
#line 1018 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDxddot_nit_P_1__).resize((*L_CSDxddotStar_nit__).size()) ;std::fill( (*L_CSDxddot_nit_P_1__).begin(),(*L_CSDxddot_nit_P_1__).end(), 0. );
	if (*L_CFDIterationFinished_nit_M_1__) {
		(*L_CSDxddot_nit_P_1__) = (*L_CSDxddotStar_nit__) ;
	} else {
		(*L_CSDxddot_nit_P_1__) = (*L_CSDxddot_nit__) ;
	}
}    void compute(const Loci::sequence &seq) { 
#line 1028 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 1028 "FSI_CSDdataloop.loci"
    }
#line 1028 "FSI_CSDdataloop.loci"
} ;
#line 1028 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop053_1280807839m96> register_file_FSI_CSDdataloop053_1280807839m96 ;
#line 1028 "FSI_CSDdataloop.loci"
}
#line 1028 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop054_1280807839m97 : public Loci::blackbox_rule {
#line 1028 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxddot_nit__ ; 
#line 1028 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxddot_n_P_1__ ; 
#line 1028 "FSI_CSDdataloop.loci"
public:
#line 1028 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop054_1280807839m97() {
#line 1028 "FSI_CSDdataloop.loci"
       name_store("CSDxddot{n,it}",L_CSDxddot_nit__) ;
#line 1028 "FSI_CSDdataloop.loci"
       name_store("CSDxddot{n+1}",L_CSDxddot_n_P_1__) ;
#line 1028 "FSI_CSDdataloop.loci"
       input("CSDxddot{n,it}") ;
#line 1028 "FSI_CSDdataloop.loci"
       output("CSDxddot{n+1}") ;
#line 1028 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 1028 "FSI_CSDdataloop.loci"
       conditional("iterationFinished{n,it-1}") ;
#line 1028 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 1028 "FSI_CSDdataloop.loci"
    }
#line 1028 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDxddot_n_P_1__).resize((*L_CSDxddot_nit__).size()) ; std::fill( (*L_CSDxddot_n_P_1__).begin(),(*L_CSDxddot_n_P_1__).end(), 0. );
	(*L_CSDxddot_n_P_1__) = (*L_CSDxddot_nit__) ;
}    void compute(const Loci::sequence &seq) { 
#line 1039 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 1039 "FSI_CSDdataloop.loci"
    }
#line 1039 "FSI_CSDdataloop.loci"
} ;
#line 1039 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop054_1280807839m97> register_file_FSI_CSDdataloop054_1280807839m97 ;
#line 1039 "FSI_CSDdataloop.loci"
}
#line 1039 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop055_1280807839m97 : public Loci::blackbox_rule {
#line 1039 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerXnum_ ; 
#line 1039 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDForcePreEuler_ic_ ; 
#line 1039 "FSI_CSDdataloop.loci"
public:
#line 1039 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop055_1280807839m97() {
#line 1039 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXnum",L_CSDEulerXnum_) ;
#line 1039 "FSI_CSDdataloop.loci"
       name_store("CSDForcePreEuler_ic",L_CSDForcePreEuler_ic_) ;
#line 1039 "FSI_CSDdataloop.loci"
       input("CSDEulerXnum") ;
#line 1039 "FSI_CSDdataloop.loci"
       output("CSDForcePreEuler_ic") ;
#line 1039 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 1039 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 1039 "FSI_CSDdataloop.loci"
    }
#line 1039 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
  
	int CSDNumBC = 2 ;
	int CSDdofNode = 2 ;
	int CSDdofFree = *L_CSDEulerXnum_* CSDdofNode - CSDNumBC ;
	
	(*L_CSDForcePreEuler_ic_).resize(CSDdofFree) ; std::fill( (*L_CSDForcePreEuler_ic_).begin(),(*L_CSDForcePreEuler_ic_).end(), 0. );
}    void compute(const Loci::sequence &seq) { 
#line 1049 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 1049 "FSI_CSDdataloop.loci"
    }
#line 1049 "FSI_CSDdataloop.loci"
} ;
#line 1049 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop055_1280807839m97> register_file_FSI_CSDdataloop055_1280807839m97 ;
#line 1049 "FSI_CSDdataloop.loci"
}
#line 1049 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop056_1280807839m97 : public Loci::blackbox_rule {
#line 1049 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDForcePreEuler_ic_ ; 
#line 1049 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDForcePreEuler_n_EQ_0__ ; 
#line 1049 "FSI_CSDdataloop.loci"
public:
#line 1049 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop056_1280807839m97() {
#line 1049 "FSI_CSDdataloop.loci"
       name_store("CSDForcePreEuler_ic",L_CSDForcePreEuler_ic_) ;
#line 1049 "FSI_CSDdataloop.loci"
       name_store("CSDForcePreEuler{n=0}",L_CSDForcePreEuler_n_EQ_0__) ;
#line 1049 "FSI_CSDdataloop.loci"
       input("CSDForcePreEuler_ic") ;
#line 1049 "FSI_CSDdataloop.loci"
       output("CSDForcePreEuler{n=0}") ;
#line 1049 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 1049 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 1049 "FSI_CSDdataloop.loci"
    }
#line 1049 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDForcePreEuler_n_EQ_0__).resize((*L_CSDForcePreEuler_ic_).size()) ;std::fill((*L_CSDForcePreEuler_n_EQ_0__).begin(),(*L_CSDForcePreEuler_n_EQ_0__).end(), 0. );
	(*L_CSDForcePreEuler_n_EQ_0__) = (*L_CSDForcePreEuler_ic_) ;
}    void compute(const Loci::sequence &seq) { 
#line 1056 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 1056 "FSI_CSDdataloop.loci"
    }
#line 1056 "FSI_CSDdataloop.loci"
} ;
#line 1056 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop056_1280807839m97> register_file_FSI_CSDdataloop056_1280807839m97 ;
#line 1056 "FSI_CSDdataloop.loci"
}
#line 1056 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop057_1280807839m98 : public Loci::blackbox_rule {
#line 1056 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDForcePreEuler_n__ ; 
#line 1056 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDForcePreEuler_nit_EQ_0__ ; 
#line 1056 "FSI_CSDdataloop.loci"
public:
#line 1056 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop057_1280807839m98() {
#line 1056 "FSI_CSDdataloop.loci"
       name_store("CSDForcePreEuler{n}",L_CSDForcePreEuler_n__) ;
#line 1056 "FSI_CSDdataloop.loci"
       name_store("CSDForcePreEuler{n,it=0}",L_CSDForcePreEuler_nit_EQ_0__) ;
#line 1056 "FSI_CSDdataloop.loci"
       input("CSDForcePreEuler{n}") ;
#line 1056 "FSI_CSDdataloop.loci"
       output("CSDForcePreEuler{n,it=0}") ;
#line 1056 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 1056 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 1056 "FSI_CSDdataloop.loci"
    }
#line 1056 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDForcePreEuler_nit_EQ_0__).resize((*L_CSDForcePreEuler_n__).size()) ;std::fill( (*L_CSDForcePreEuler_nit_EQ_0__).begin(),(*L_CSDForcePreEuler_nit_EQ_0__).end(), 0. );
	(*L_CSDForcePreEuler_nit_EQ_0__) = (*L_CSDForcePreEuler_n__) ;
//	if (Loci::MPI_rank==0) cout << "Inside CSDForcePre{n,it=0}" << endl;
}    void compute(const Loci::sequence &seq) { 
#line 1063 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 1063 "FSI_CSDdataloop.loci"
    }
#line 1063 "FSI_CSDdataloop.loci"
} ;
#line 1063 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop057_1280807839m98> register_file_FSI_CSDdataloop057_1280807839m98 ;
#line 1063 "FSI_CSDdataloop.loci"
}
#line 1063 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop058_1280807839m98 : public Loci::blackbox_rule {
#line 1063 "FSI_CSDdataloop.loci"
    Loci::const_param<bool>  L_CFDIterationFinished_nit_M_1__ ; 
#line 1063 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDForcePreEulerStar_nit__ ; 
#line 1063 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDForcePreEuler_nit__ ; 
#line 1063 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDForcePreEuler_nit_P_1__ ; 
#line 1063 "FSI_CSDdataloop.loci"
public:
#line 1063 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop058_1280807839m98() {
#line 1063 "FSI_CSDdataloop.loci"
       name_store("CFDIterationFinished{n,it-1}",L_CFDIterationFinished_nit_M_1__) ;
#line 1063 "FSI_CSDdataloop.loci"
       name_store("CSDForcePreEulerStar{n,it}",L_CSDForcePreEulerStar_nit__) ;
#line 1063 "FSI_CSDdataloop.loci"
       name_store("CSDForcePreEuler{n,it}",L_CSDForcePreEuler_nit__) ;
#line 1063 "FSI_CSDdataloop.loci"
       name_store("CSDForcePreEuler{n,it+1}",L_CSDForcePreEuler_nit_P_1__) ;
#line 1063 "FSI_CSDdataloop.loci"
       input("CSDForcePreEuler{n,it},CSDForcePreEulerStar{n,it},CFDIterationFinished{n,it-1}") ;
#line 1063 "FSI_CSDdataloop.loci"
       output("CSDForcePreEuler{n,it+1}") ;
#line 1063 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 1063 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 1063 "FSI_CSDdataloop.loci"
    }
#line 1063 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDForcePreEuler_nit_P_1__).resize((*L_CSDForcePreEulerStar_nit__).size()) ;std::fill( (*L_CSDForcePreEuler_nit_P_1__).begin(),(*L_CSDForcePreEuler_nit_P_1__).end(), 0. );
	
	if (*L_CFDIterationFinished_nit_M_1__) {
		(*L_CSDForcePreEuler_nit_P_1__) = (*L_CSDForcePreEulerStar_nit__) ;
	} else {
		(*L_CSDForcePreEuler_nit_P_1__) = (*L_CSDForcePreEuler_nit__) ;
	}
}    void compute(const Loci::sequence &seq) { 
#line 1079 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 1079 "FSI_CSDdataloop.loci"
    }
#line 1079 "FSI_CSDdataloop.loci"
} ;
#line 1079 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop058_1280807839m98> register_file_FSI_CSDdataloop058_1280807839m98 ;
#line 1079 "FSI_CSDdataloop.loci"
}
#line 1079 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop059_1280807839m99 : public Loci::blackbox_rule {
#line 1079 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDForcePreEuler_nit__ ; 
#line 1079 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDForcePreEuler_n_P_1__ ; 
#line 1079 "FSI_CSDdataloop.loci"
public:
#line 1079 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop059_1280807839m99() {
#line 1079 "FSI_CSDdataloop.loci"
       name_store("CSDForcePreEuler{n,it}",L_CSDForcePreEuler_nit__) ;
#line 1079 "FSI_CSDdataloop.loci"
       name_store("CSDForcePreEuler{n+1}",L_CSDForcePreEuler_n_P_1__) ;
#line 1079 "FSI_CSDdataloop.loci"
       input("CSDForcePreEuler{n,it}") ;
#line 1079 "FSI_CSDdataloop.loci"
       output("CSDForcePreEuler{n+1}") ;
#line 1079 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 1079 "FSI_CSDdataloop.loci"
       conditional("iterationFinished{n,it-1}") ;
#line 1079 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 1079 "FSI_CSDdataloop.loci"
    }
#line 1079 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDForcePreEuler_n_P_1__).resize((*L_CSDForcePreEuler_nit__).size()) ;std::fill( (*L_CSDForcePreEuler_n_P_1__).begin(),(*L_CSDForcePreEuler_n_P_1__).end(), 0. );
	(*L_CSDForcePreEuler_n_P_1__) = (*L_CSDForcePreEuler_nit__) ;
}    void compute(const Loci::sequence &seq) { 
#line 1084 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 1084 "FSI_CSDdataloop.loci"
    }
#line 1084 "FSI_CSDdataloop.loci"
} ;
#line 1084 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop059_1280807839m99> register_file_FSI_CSDdataloop059_1280807839m99 ;
#line 1084 "FSI_CSDdataloop.loci"
}
#line 1084 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop060_1280807839m99 : public Loci::blackbox_rule {
#line 1084 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDEulerXstart_ ; 
#line 1084 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDEulerXend_ ; 
#line 1084 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerXnum_ ; 
#line 1084 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 1084 "FSI_CSDdataloop.loci"
public:
#line 1084 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop060_1280807839m99() {
#line 1084 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 1084 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXstart",L_CSDEulerXstart_) ;
#line 1084 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXend",L_CSDEulerXend_) ;
#line 1084 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXnum",L_CSDEulerXnum_) ;
#line 1084 "FSI_CSDdataloop.loci"
       input("CSDEulerXstart,CSDEulerXend,CSDEulerXnum") ;
#line 1084 "FSI_CSDdataloop.loci"
       output("CSDnodes_ic") ;
#line 1084 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIPRESCRIBED") ;
#line 1084 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 1084 "FSI_CSDdataloop.loci"
    }
#line 1084 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) { 
	int nElements = *L_CSDEulerXnum_- 1;
	real dx = (*L_CSDEulerXend_- *L_CSDEulerXstart_) / nElements ;
	
	(*L_CSDnodes_ic_).resize(*L_CSDEulerXnum_,3) ; std::fill( (*L_CSDnodes_ic_).data().begin(), (*L_CSDnodes_ic_).data().end(), 0. );	    	
	
	(*L_CSDnodes_ic_)(0,0) = *L_CSDEulerXstart_;
	for(int i=0; i<*L_CSDEulerXnum_; ++i) {
		(*L_CSDnodes_ic_)(i,0) = *L_CSDEulerXstart_+ i * dx;
	}
}    void compute(const Loci::sequence &seq) { 
#line 1096 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 1096 "FSI_CSDdataloop.loci"
    }
#line 1096 "FSI_CSDdataloop.loci"
} ;
#line 1096 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop060_1280807839m99> register_file_FSI_CSDdataloop060_1280807839m99 ;
#line 1096 "FSI_CSDdataloop.loci"
}
#line 1096 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop061_1280807839m100 : public Loci::blackbox_rule {
#line 1100 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDfrequency_ ; 
#line 1100 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDflappingAmplitudeX_ ; 
#line 1100 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDstartingTimeStep_ ; 
#line 1100 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 1100 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerXnum_ ; 
#line 1100 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_stime_n__ ; 
#line 1100 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_ncycle_n__ ; 
#line 1100 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_timeStep_ ; 
#line 1100 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForce_nit__ ; 
#line 1100 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDdisplacementsStar_nit__ ; 
#line 1100 "FSI_CSDdataloop.loci"
public:
#line 1100 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop061_1280807839m100() {
#line 1100 "FSI_CSDdataloop.loci"
       name_store("CSDfrequency",L_CSDfrequency_) ;
#line 1100 "FSI_CSDdataloop.loci"
       name_store("CSDflappingAmplitudeX",L_CSDflappingAmplitudeX_) ;
#line 1100 "FSI_CSDdataloop.loci"
       name_store("CSDstartingTimeStep",L_CSDstartingTimeStep_) ;
#line 1100 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 1100 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXnum",L_CSDEulerXnum_) ;
#line 1100 "FSI_CSDdataloop.loci"
       name_store("CSDdisplacementsStar{n,it}",L_CSDdisplacementsStar_nit__) ;
#line 1100 "FSI_CSDdataloop.loci"
       name_store("stime{n}",L_stime_n__) ;
#line 1100 "FSI_CSDdataloop.loci"
       name_store("ncycle{n}",L_ncycle_n__) ;
#line 1100 "FSI_CSDdataloop.loci"
       name_store("timeStep",L_timeStep_) ;
#line 1100 "FSI_CSDdataloop.loci"
       name_store("CSDForce{n,it}",L_CSDForce_nit__) ;
#line 1100 "FSI_CSDdataloop.loci"
       input("							CSDEulerXnum,CSDnodes_ic,CSDForce{n,it},							stime{n},ncycle{n},timeStep,CSDstartingTimeStep,								CSDfrequency,CSDflappingAmplitudeX") ;
#line 1100 "FSI_CSDdataloop.loci"
       output("CSDdisplacementsStar{n,it}") ;
#line 1100 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIPRESCRIBED") ;
#line 1100 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 1100 "FSI_CSDdataloop.loci"
    }
#line 1100 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) { 
	// call NLAMS
	//communicateWithNLAMS(*CSDdisplacementsStar, ....) ;
	//if (*$CFDIterationFinished) {
	
	const int rank = Loci::MPI_rank ;

	if (rank==0) cout << "[I] Communicating with the prescribed motion rule " << endl ;
	// Communicate with NLAMS
	int CSDNumNodes = (*L_CSDEulerXnum_) ;
	
	// Initialize
	(*L_CSDdisplacementsStar_nit__).resize(CSDNumNodes,3) ; std::fill( (*L_CSDdisplacementsStar_nit__).data().begin(),(*L_CSDdisplacementsStar_nit__).data().end(), 0. );	
	
		if (*L_ncycle_n__ < *L_CSDstartingTimeStep_) {
			 if (rank==0) cout << "CSDstartingTime not yet reached. Current time step = " << *L_ncycle_n__ << ", CSDstartingTimeStep = " << *L_CSDstartingTimeStep_<< endl ;
		} else {	
			int tempTimeStepNumber = *L_ncycle_n__ - *L_CSDstartingTimeStep_+ 1 ;	
				
			std::vector<real> CSDbeam1dForce(CSDNumNodes, 0.0) ;	
			std::vector<real> disp(CSDNumNodes, 0.0) ;
			std::vector<real> nodes(CSDNumNodes, 0.0) ;
				
			const double PI = 4.0*atan(1.0) ;	
  		real alpha_amp = *L_CSDflappingAmplitudeX_; 
		  real freq = *L_CSDfrequency_; // period in sec
  		real alpha = alpha_amp * sin(2.*PI * freq * (*L_stime_n__)) ;  		

  		for(int i=0;i<CSDNumNodes;++i) {
		 	  (*L_CSDdisplacementsStar_nit__)(i,0) = (*L_CSDnodes_ic_)(i,0) * ( cos(alpha*PI/180) - 1.) ;
		 	  (*L_CSDdisplacementsStar_nit__)(i,1) = (*L_CSDnodes_ic_)(i,0) * ( sin(alpha*PI/180) ) ;
		 	  (*L_CSDdisplacementsStar_nit__)(i,2) = 0.0 ;
		  }
	   	if (Loci::MPI_rank==0) cout << "CSDdisplacement:" << (*L_stime_n__) << ", alpha = " << alpha << endl ;
	  	//if (Loci::MPI_rank==0) cout << "CSDnodes_ic: " << (*$CSDnodes_ic) << endl ;
	  	//if (Loci::MPI_rank==0) cout << "CSDdisplacementsStar: " << (*$CSDdisplacementsStar{n,it}) << endl ;	
		  
		  if (rank==0) cout << "[I] Communication with the prescribed motion rule done. " << endl ;
		}   
	
}    void compute(const Loci::sequence &seq) { 
#line 1143 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 1143 "FSI_CSDdataloop.loci"
    }
#line 1143 "FSI_CSDdataloop.loci"
} ;
#line 1143 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop061_1280807839m100> register_file_FSI_CSDdataloop061_1280807839m100 ;
#line 1143 "FSI_CSDdataloop.loci"
}
#line 1143 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop062_1280807839m101 : public Loci::singleton_rule {
#line 1143 "FSI_CSDdataloop.loci"
    Loci::param<int>  L_itfsi_ic_ ; 
#line 1143 "FSI_CSDdataloop.loci"
public:
#line 1143 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop062_1280807839m101() {
#line 1143 "FSI_CSDdataloop.loci"
       name_store("itfsi_ic",L_itfsi_ic_) ;
#line 1143 "FSI_CSDdataloop.loci"
       output("itfsi_ic") ;
#line 1143 "FSI_CSDdataloop.loci"
       constraint("UNIVERSE") ;
#line 1143 "FSI_CSDdataloop.loci"
    }
#line 1143 "FSI_CSDdataloop.loci"
    void compute(const Loci::sequence &seq) { 
	(*L_itfsi_ic_)= 0 ;
}} ;
#line 1145 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop062_1280807839m101> register_file_FSI_CSDdataloop062_1280807839m101 ;
#line 1145 "FSI_CSDdataloop.loci"
}
#line 1145 "FSI_CSDdataloop.loci"


namespace {class file_FSI_CSDdataloop063_1280807839m101 : public Loci::singleton_rule {
#line 1147 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_itfsi_ic_ ; 
#line 1147 "FSI_CSDdataloop.loci"
    Loci::param<int>  L_itfsi_n_EQ_0__ ; 
#line 1147 "FSI_CSDdataloop.loci"
public:
#line 1147 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop063_1280807839m101() {
#line 1147 "FSI_CSDdataloop.loci"
       name_store("itfsi_ic",L_itfsi_ic_) ;
#line 1147 "FSI_CSDdataloop.loci"
       name_store("itfsi{n=0}",L_itfsi_n_EQ_0__) ;
#line 1147 "FSI_CSDdataloop.loci"
       input("itfsi_ic") ;
#line 1147 "FSI_CSDdataloop.loci"
       output("itfsi{n=0}") ;
#line 1147 "FSI_CSDdataloop.loci"
    }
#line 1147 "FSI_CSDdataloop.loci"
    void compute(const Loci::sequence &seq) { 
	(*L_itfsi_n_EQ_0__)=(*L_itfsi_ic_);
}} ;
#line 1149 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop063_1280807839m101> register_file_FSI_CSDdataloop063_1280807839m101 ;
#line 1149 "FSI_CSDdataloop.loci"
}
#line 1149 "FSI_CSDdataloop.loci"


namespace {class file_FSI_CSDdataloop064_1280807839m102 : public Loci::singleton_rule {
#line 1151 "FSI_CSDdataloop.loci"
    Loci::param<int>  L_itfsi_n_P_1__ ; 
#line 1151 "FSI_CSDdataloop.loci"
public:
#line 1151 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop064_1280807839m102() {
#line 1151 "FSI_CSDdataloop.loci"
       name_store("itfsi{n+1}",L_itfsi_n_P_1__) ;
#line 1151 "FSI_CSDdataloop.loci"
       output("itfsi{n+1}") ;
#line 1151 "FSI_CSDdataloop.loci"
       conditional("iterationFinished{n,it-1}") ;
#line 1151 "FSI_CSDdataloop.loci"
    }
#line 1151 "FSI_CSDdataloop.loci"
    void compute(const Loci::sequence &seq) { 
	
	//if (Loci::MPI_rank==0) cout << "inside itfsi{n+1} in" << endl ;
	(*L_itfsi_n_P_1__) = 0 ;
	//if (Loci::MPI_rank==0) cout << "inside itfsi{n+1} out" << endl ;
	
}} ;
#line 1157 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop064_1280807839m102> register_file_FSI_CSDdataloop064_1280807839m102 ;
#line 1157 "FSI_CSDdataloop.loci"
}
#line 1157 "FSI_CSDdataloop.loci"


}
