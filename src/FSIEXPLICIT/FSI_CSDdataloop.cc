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
   
      	void beam1d_(const int *, const int*, const int*, const double*, const double*, const double*, const double*, const double*, const double*, const double*,
        const double*, const double*, const double*, const double*, const double*, const double*, const double*, const int*, 
        double*, double*, double*, double*) ;
        
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
//$type CSDForcePre_ic blackbox<ublas::vector<real> > ;
//$type CSDForcePre blackbox<ublas::vector<real> > ;
//$type CSDForcePreStar blackbox<ublas::vector<real> > ;
#line 63 "FSI_CSDdataloop.loci"


namespace streamUns {
            
//-----------------------------------------------------------------------------
// Rules to read in the grids

// -- Read in the nodes, connectivity, bc files for CSD------------------------------------------------------------------------------------- 
// Read in the initial CSD nodes: undeformed

// $type CFDIterationFinished param<bool> 
// $type FSIIterationFinished param<bool> 


namespace {class file_FSI_CSDdataloop000_1279551148m526 : public Loci::blackbox_rule {
#line 78 "FSI_CSDdataloop.loci"
    Loci::const_param<string>  L_CSDMeshFilename_ ; 
#line 78 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 78 "FSI_CSDdataloop.loci"
public:
#line 78 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop000_1279551148m526() {
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
Loci::register_rule<file_FSI_CSDdataloop000_1279551148m526> register_file_FSI_CSDdataloop000_1279551148m526 ;
#line 122 "FSI_CSDdataloop.loci"
}
#line 122 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop001_1279551148m526 : public Loci::blackbox_rule {
#line 123 "FSI_CSDdataloop.loci"
    Loci::const_param<string>  L_CSDConnectivityFilename_ ; 
#line 123 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<int,ublas::column_major> >  L_CSDConnectivity_ ; 
#line 123 "FSI_CSDdataloop.loci"
public:
#line 123 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop001_1279551148m526() {
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
Loci::register_rule<file_FSI_CSDdataloop001_1279551148m526> register_file_FSI_CSDdataloop001_1279551148m526 ;
#line 167 "FSI_CSDdataloop.loci"
}
#line 167 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop002_1279551148m527 : public Loci::singleton_rule {
#line 167 "FSI_CSDdataloop.loci"
    Loci::const_param<string>  L_CSDBCFilename_ ; 
#line 167 "FSI_CSDdataloop.loci"
    Loci::param<vector<int> >  L_CSDBCdof_ ; 
#line 167 "FSI_CSDdataloop.loci"
    Loci::param<vector<real> >  L_CSDBCZeroConstraint_ ; 
#line 167 "FSI_CSDdataloop.loci"
public:
#line 167 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop002_1279551148m527() {
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
Loci::register_rule<file_FSI_CSDdataloop002_1279551148m527> register_file_FSI_CSDdataloop002_1279551148m527 ;
#line 213 "FSI_CSDdataloop.loci"
}
#line 213 "FSI_CSDdataloop.loci"


// -- Time advancing for the CSDnodes ----------------------------------------------------------------------------------------------------
// Time Build rule for the CSD nodes
namespace {class file_FSI_CSDdataloop003_1279551148m527 : public Loci::blackbox_rule {
#line 218 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 218 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_n_EQ_0__ ; 
#line 218 "FSI_CSDdataloop.loci"
public:
#line 218 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop003_1279551148m527() {
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
Loci::register_rule<file_FSI_CSDdataloop003_1279551148m527> register_file_FSI_CSDdataloop003_1279551148m527 ;
#line 224 "FSI_CSDdataloop.loci"
}
#line 224 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop004_1279551148m528 : public Loci::blackbox_rule {
#line 225 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 225 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDdisplacementsStar_n__ ; 
#line 225 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_n_P_1__ ; 
#line 225 "FSI_CSDdataloop.loci"
public:
#line 225 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop004_1279551148m528() {
#line 225 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 225 "FSI_CSDdataloop.loci"
       name_store("CSDdisplacementsStar{n}",L_CSDdisplacementsStar_n__) ;
#line 225 "FSI_CSDdataloop.loci"
       name_store("CSDnodes{n+1}",L_CSDnodes_n_P_1__) ;
#line 225 "FSI_CSDdataloop.loci"
       input("CSDnodes_ic,CSDdisplacementsStar{n}") ;
#line 225 "FSI_CSDdataloop.loci"
       output("CSDnodes{n+1}") ;
#line 225 "FSI_CSDdataloop.loci"
       constraint("FSICoupling") ;
#line 225 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 225 "FSI_CSDdataloop.loci"
    }
#line 225 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  	
	(*L_CSDnodes_n_P_1__).resize((*L_CSDnodes_ic_).size1(),3) ; std::fill( (*L_CSDnodes_n_P_1__).data().begin(), (*L_CSDnodes_n_P_1__).data().end(), 0. ); 
			
	(*L_CSDnodes_n_P_1__) = (*L_CSDnodes_ic_) + (*L_CSDdisplacementsStar_n__); // CHECK!!!!
	//	if (Loci::MPI_rank==0) cout << "Inside CSDnodes_it{n,it+1} priority out" << endl ;	
}    void compute(const Loci::sequence &seq) { 
#line 234 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 234 "FSI_CSDdataloop.loci"
    }
#line 234 "FSI_CSDdataloop.loci"
} ;
#line 234 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop004_1279551148m528> register_file_FSI_CSDdataloop004_1279551148m528 ;
#line 234 "FSI_CSDdataloop.loci"
}
#line 234 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop005_1279551148m528 : public Loci::blackbox_rule {
#line 234 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 234 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDisp_ic_ ; 
#line 234 "FSI_CSDdataloop.loci"
public:
#line 234 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop005_1279551148m528() {
#line 234 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 234 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDisp_ic",L_CSDnodesDisp_ic_) ;
#line 234 "FSI_CSDdataloop.loci"
       input("CSDnodes_ic") ;
#line 234 "FSI_CSDdataloop.loci"
       output("CSDnodesDisp_ic") ;
#line 234 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 234 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 234 "FSI_CSDdataloop.loci"
    }
#line 234 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesDisp_ic_).resize(6*(*L_CSDnodes_ic_).size1(),1) ;
	std::fill( (*L_CSDnodesDisp_ic_).data().begin(),(*L_CSDnodesDisp_ic_).data().end(), 0. ); 
}    void compute(const Loci::sequence &seq) { 
#line 240 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 240 "FSI_CSDdataloop.loci"
    }
#line 240 "FSI_CSDdataloop.loci"
} ;
#line 240 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop005_1279551148m528> register_file_FSI_CSDdataloop005_1279551148m528 ;
#line 240 "FSI_CSDdataloop.loci"
}
#line 240 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop006_1279551148m529 : public Loci::blackbox_rule {
#line 240 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDisp_ic_ ; 
#line 240 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDisp_n_EQ_0__ ; 
#line 240 "FSI_CSDdataloop.loci"
public:
#line 240 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop006_1279551148m529() {
#line 240 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDisp_ic",L_CSDnodesDisp_ic_) ;
#line 240 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDisp{n=0}",L_CSDnodesDisp_n_EQ_0__) ;
#line 240 "FSI_CSDdataloop.loci"
       input("CSDnodesDisp_ic") ;
#line 240 "FSI_CSDdataloop.loci"
       output("CSDnodesDisp{n=0}") ;
#line 240 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 240 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 240 "FSI_CSDdataloop.loci"
    }
#line 240 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesDisp_n_EQ_0__).resize((*L_CSDnodesDisp_ic_).size1(),1) ; std::fill( (*L_CSDnodesDisp_n_EQ_0__).data().begin(),(*L_CSDnodesDisp_n_EQ_0__).data().end(), 0. ); 
	(*L_CSDnodesDisp_n_EQ_0__) = (*L_CSDnodesDisp_ic_) ;
}    void compute(const Loci::sequence &seq) { 
#line 247 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 247 "FSI_CSDdataloop.loci"
    }
#line 247 "FSI_CSDdataloop.loci"
} ;
#line 247 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop006_1279551148m529> register_file_FSI_CSDdataloop006_1279551148m529 ;
#line 247 "FSI_CSDdataloop.loci"
}
#line 247 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop007_1279551148m529 : public Loci::blackbox_rule {
#line 247 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDisp_n__ ; 
#line 247 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDispStar_n__ ; 
#line 247 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDisp_n_P_1__ ; 
#line 247 "FSI_CSDdataloop.loci"
public:
#line 247 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop007_1279551148m529() {
#line 247 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDisp{n}",L_CSDnodesDisp_n__) ;
#line 247 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDispStar{n}",L_CSDnodesDispStar_n__) ;
#line 247 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDisp{n+1}",L_CSDnodesDisp_n_P_1__) ;
#line 247 "FSI_CSDdataloop.loci"
       input("CSDnodesDisp{n},CSDnodesDispStar{n}") ;
#line 247 "FSI_CSDdataloop.loci"
       output("CSDnodesDisp{n+1}") ;
#line 247 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 247 "FSI_CSDdataloop.loci"
       conditional("iterationFinished{n,it-1}") ;
#line 247 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 247 "FSI_CSDdataloop.loci"
    }
#line 247 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesDisp_n_P_1__).resize((*L_CSDnodesDispStar_n__).size1(),1) ;	std::fill( (*L_CSDnodesDisp_n_P_1__).data().begin(),(*L_CSDnodesDisp_n_P_1__).data().end(), 0. );
	(*L_CSDnodesDisp_n_P_1__) = (*L_CSDnodesDispStar_n__) ;	
}    void compute(const Loci::sequence &seq) { 
#line 254 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 254 "FSI_CSDdataloop.loci"
    }
#line 254 "FSI_CSDdataloop.loci"
} ;
#line 254 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop007_1279551148m529> register_file_FSI_CSDdataloop007_1279551148m529 ;
#line 254 "FSI_CSDdataloop.loci"
}
#line 254 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop008_1279551148m529 : public Loci::blackbox_rule {
#line 254 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 254 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVel_ic_ ; 
#line 254 "FSI_CSDdataloop.loci"
public:
#line 254 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop008_1279551148m529() {
#line 254 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 254 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVel_ic",L_CSDnodesVel_ic_) ;
#line 254 "FSI_CSDdataloop.loci"
       input("CSDnodes_ic") ;
#line 254 "FSI_CSDdataloop.loci"
       output("CSDnodesVel_ic") ;
#line 254 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 254 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 254 "FSI_CSDdataloop.loci"
    }
#line 254 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesVel_ic_).resize(6*(*L_CSDnodes_ic_).size1(),1) ; std::fill( (*L_CSDnodesVel_ic_).data().begin(),(*L_CSDnodesVel_ic_).data().end(), 0. );
	(*L_CSDnodesVel_ic_).clear() ; // Set all elements to zero
}    void compute(const Loci::sequence &seq) { 
#line 260 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 260 "FSI_CSDdataloop.loci"
    }
#line 260 "FSI_CSDdataloop.loci"
} ;
#line 260 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop008_1279551148m529> register_file_FSI_CSDdataloop008_1279551148m529 ;
#line 260 "FSI_CSDdataloop.loci"
}
#line 260 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop009_1279551148m530 : public Loci::blackbox_rule {
#line 260 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVel_ic_ ; 
#line 260 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVel_n_EQ_0__ ; 
#line 260 "FSI_CSDdataloop.loci"
public:
#line 260 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop009_1279551148m530() {
#line 260 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVel_ic",L_CSDnodesVel_ic_) ;
#line 260 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVel{n=0}",L_CSDnodesVel_n_EQ_0__) ;
#line 260 "FSI_CSDdataloop.loci"
       input("CSDnodesVel_ic") ;
#line 260 "FSI_CSDdataloop.loci"
       output("CSDnodesVel{n=0}") ;
#line 260 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 260 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 260 "FSI_CSDdataloop.loci"
    }
#line 260 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesVel_n_EQ_0__).resize((*L_CSDnodesVel_ic_).size1(),1) ; std::fill( (*L_CSDnodesVel_n_EQ_0__).data().begin(),(*L_CSDnodesVel_n_EQ_0__).data().end(), 0. );
	(*L_CSDnodesVel_n_EQ_0__) = (*L_CSDnodesVel_ic_) ;
}    void compute(const Loci::sequence &seq) { 
#line 266 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 266 "FSI_CSDdataloop.loci"
    }
#line 266 "FSI_CSDdataloop.loci"
} ;
#line 266 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop009_1279551148m530> register_file_FSI_CSDdataloop009_1279551148m530 ;
#line 266 "FSI_CSDdataloop.loci"
}
#line 266 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop010_1279551148m530 : public Loci::blackbox_rule {
#line 266 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVel_n__ ; 
#line 266 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVelStar_n__ ; 
#line 266 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVel_n_P_1__ ; 
#line 266 "FSI_CSDdataloop.loci"
public:
#line 266 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop010_1279551148m530() {
#line 266 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVel{n}",L_CSDnodesVel_n__) ;
#line 266 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVelStar{n}",L_CSDnodesVelStar_n__) ;
#line 266 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVel{n+1}",L_CSDnodesVel_n_P_1__) ;
#line 266 "FSI_CSDdataloop.loci"
       input("CSDnodesVel{n},CSDnodesVelStar{n}") ;
#line 266 "FSI_CSDdataloop.loci"
       output("CSDnodesVel{n+1}") ;
#line 266 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 266 "FSI_CSDdataloop.loci"
       conditional("iterationFinished{n,it-1}") ;
#line 266 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 266 "FSI_CSDdataloop.loci"
    }
#line 266 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesVel_n_P_1__).resize((*L_CSDnodesVelStar_n__).size1(),1) ; std::fill( (*L_CSDnodesVel_n_P_1__).data().begin(),(*L_CSDnodesVel_n_P_1__).data().end(), 0. );
	(*L_CSDnodesVel_n_P_1__) = (*L_CSDnodesVelStar_n__) ;	
}    void compute(const Loci::sequence &seq) { 
#line 273 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 273 "FSI_CSDdataloop.loci"
    }
#line 273 "FSI_CSDdataloop.loci"
} ;
#line 273 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop010_1279551148m530> register_file_FSI_CSDdataloop010_1279551148m530 ;
#line 273 "FSI_CSDdataloop.loci"
}
#line 273 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop011_1279551148m531 : public Loci::blackbox_rule {
#line 273 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 273 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAcc_ic_ ; 
#line 273 "FSI_CSDdataloop.loci"
public:
#line 273 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop011_1279551148m531() {
#line 273 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 273 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAcc_ic",L_CSDnodesAcc_ic_) ;
#line 273 "FSI_CSDdataloop.loci"
       input("CSDnodes_ic") ;
#line 273 "FSI_CSDdataloop.loci"
       output("CSDnodesAcc_ic") ;
#line 273 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 273 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 273 "FSI_CSDdataloop.loci"
    }
#line 273 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesAcc_ic_).resize(6*(*L_CSDnodes_ic_).size1(),1) ; std::fill( (*L_CSDnodesAcc_ic_).data().begin(),(*L_CSDnodesAcc_ic_).data().end(), 0. );
	(*L_CSDnodesAcc_ic_).clear() ;
}    void compute(const Loci::sequence &seq) { 
#line 279 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 279 "FSI_CSDdataloop.loci"
    }
#line 279 "FSI_CSDdataloop.loci"
} ;
#line 279 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop011_1279551148m531> register_file_FSI_CSDdataloop011_1279551148m531 ;
#line 279 "FSI_CSDdataloop.loci"
}
#line 279 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop012_1279551148m531 : public Loci::blackbox_rule {
#line 279 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAcc_ic_ ; 
#line 279 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAcc_n_EQ_0__ ; 
#line 279 "FSI_CSDdataloop.loci"
public:
#line 279 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop012_1279551148m531() {
#line 279 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAcc_ic",L_CSDnodesAcc_ic_) ;
#line 279 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAcc{n=0}",L_CSDnodesAcc_n_EQ_0__) ;
#line 279 "FSI_CSDdataloop.loci"
       input("CSDnodesAcc_ic") ;
#line 279 "FSI_CSDdataloop.loci"
       output("CSDnodesAcc{n=0}") ;
#line 279 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 279 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 279 "FSI_CSDdataloop.loci"
    }
#line 279 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesAcc_n_EQ_0__).resize((*L_CSDnodesAcc_ic_).size1(),1) ;std::fill( (*L_CSDnodesAcc_n_EQ_0__).data().begin(),(*L_CSDnodesAcc_n_EQ_0__).data().end(), 0. );
	(*L_CSDnodesAcc_n_EQ_0__) = (*L_CSDnodesAcc_ic_) ;
}    void compute(const Loci::sequence &seq) { 
#line 285 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 285 "FSI_CSDdataloop.loci"
    }
#line 285 "FSI_CSDdataloop.loci"
} ;
#line 285 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop012_1279551148m531> register_file_FSI_CSDdataloop012_1279551148m531 ;
#line 285 "FSI_CSDdataloop.loci"
}
#line 285 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop013_1279551148m531 : public Loci::blackbox_rule {
#line 285 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAcc_n__ ; 
#line 285 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAccStar_n__ ; 
#line 285 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAcc_n_P_1__ ; 
#line 285 "FSI_CSDdataloop.loci"
public:
#line 285 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop013_1279551148m531() {
#line 285 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAcc{n}",L_CSDnodesAcc_n__) ;
#line 285 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAccStar{n}",L_CSDnodesAccStar_n__) ;
#line 285 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAcc{n+1}",L_CSDnodesAcc_n_P_1__) ;
#line 285 "FSI_CSDdataloop.loci"
       input("CSDnodesAcc{n},CSDnodesAccStar{n}") ;
#line 285 "FSI_CSDdataloop.loci"
       output("CSDnodesAcc{n+1}") ;
#line 285 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 285 "FSI_CSDdataloop.loci"
       conditional("iterationFinished{n,it-1}") ;
#line 285 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 285 "FSI_CSDdataloop.loci"
    }
#line 285 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDnodesAcc_n_P_1__).resize((*L_CSDnodesAccStar_n__).size1(),1) ;std::fill( (*L_CSDnodesAcc_n_P_1__).data().begin(),(*L_CSDnodesAcc_n_P_1__).data().end(), 0. );
	(*L_CSDnodesAcc_n_P_1__) = (*L_CSDnodesAccStar_n__) ;	
}    void compute(const Loci::sequence &seq) { 
#line 292 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 292 "FSI_CSDdataloop.loci"
    }
#line 292 "FSI_CSDdataloop.loci"
} ;
#line 292 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop013_1279551148m531> register_file_FSI_CSDdataloop013_1279551148m531 ;
#line 292 "FSI_CSDdataloop.loci"
}
#line 292 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop014_1279551148m532 : public Loci::blackbox_rule {
#line 292 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 292 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePre_ic_ ; 
#line 292 "FSI_CSDdataloop.loci"
public:
#line 292 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop014_1279551148m532() {
#line 292 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 292 "FSI_CSDdataloop.loci"
       name_store("CSDForcePre_ic",L_CSDForcePre_ic_) ;
#line 292 "FSI_CSDdataloop.loci"
       input("CSDnodes_ic") ;
#line 292 "FSI_CSDdataloop.loci"
       output("CSDForcePre_ic") ;
#line 292 "FSI_CSDdataloop.loci"
       constraint("FSICoupling") ;
#line 292 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 292 "FSI_CSDdataloop.loci"
    }
#line 292 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDForcePre_ic_).resize(6*(*L_CSDnodes_ic_).size1(),1) ; std::fill( (*L_CSDForcePre_ic_).data().begin(),(*L_CSDForcePre_ic_).data().end(), 0. );
	(*L_CSDForcePre_ic_).clear();
}    void compute(const Loci::sequence &seq) { 
#line 298 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 298 "FSI_CSDdataloop.loci"
    }
#line 298 "FSI_CSDdataloop.loci"
} ;
#line 298 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop014_1279551148m532> register_file_FSI_CSDdataloop014_1279551148m532 ;
#line 298 "FSI_CSDdataloop.loci"
}
#line 298 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop015_1279551148m532 : public Loci::blackbox_rule {
#line 298 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePre_ic_ ; 
#line 298 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePre_n_EQ_0__ ; 
#line 298 "FSI_CSDdataloop.loci"
public:
#line 298 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop015_1279551148m532() {
#line 298 "FSI_CSDdataloop.loci"
       name_store("CSDForcePre_ic",L_CSDForcePre_ic_) ;
#line 298 "FSI_CSDdataloop.loci"
       name_store("CSDForcePre{n=0}",L_CSDForcePre_n_EQ_0__) ;
#line 298 "FSI_CSDdataloop.loci"
       input("CSDForcePre_ic") ;
#line 298 "FSI_CSDdataloop.loci"
       output("CSDForcePre{n=0}") ;
#line 298 "FSI_CSDdataloop.loci"
       constraint("FSICoupling") ;
#line 298 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 298 "FSI_CSDdataloop.loci"
    }
#line 298 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDForcePre_n_EQ_0__).resize((*L_CSDForcePre_ic_).size1(),1) ;std::fill((*L_CSDForcePre_n_EQ_0__).data().begin(),(*L_CSDForcePre_n_EQ_0__).data().end(), 0. );
	(*L_CSDForcePre_n_EQ_0__) = (*L_CSDForcePre_ic_) ;
}    void compute(const Loci::sequence &seq) { 
#line 304 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 304 "FSI_CSDdataloop.loci"
    }
#line 304 "FSI_CSDdataloop.loci"
} ;
#line 304 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop015_1279551148m532> register_file_FSI_CSDdataloop015_1279551148m532 ;
#line 304 "FSI_CSDdataloop.loci"
}
#line 304 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop016_1279551148m533 : public Loci::blackbox_rule {
#line 304 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePre_n__ ; 
#line 304 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePreStar_n__ ; 
#line 304 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePre_n_P_1__ ; 
#line 304 "FSI_CSDdataloop.loci"
public:
#line 304 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop016_1279551148m533() {
#line 304 "FSI_CSDdataloop.loci"
       name_store("CSDForcePre{n}",L_CSDForcePre_n__) ;
#line 304 "FSI_CSDdataloop.loci"
       name_store("CSDForcePreStar{n}",L_CSDForcePreStar_n__) ;
#line 304 "FSI_CSDdataloop.loci"
       name_store("CSDForcePre{n+1}",L_CSDForcePre_n_P_1__) ;
#line 304 "FSI_CSDdataloop.loci"
       input("CSDForcePre{n},CSDForcePreStar{n}") ;
#line 304 "FSI_CSDdataloop.loci"
       output("CSDForcePre{n+1}") ;
#line 304 "FSI_CSDdataloop.loci"
       constraint("FSICoupling") ;
#line 304 "FSI_CSDdataloop.loci"
       conditional("iterationFinished{n,it-1}") ;
#line 304 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 304 "FSI_CSDdataloop.loci"
    }
#line 304 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDForcePre_n_P_1__).resize((*L_CSDForcePreStar_n__).size1(),1) ;std::fill( (*L_CSDForcePre_n_P_1__).data().begin(),(*L_CSDForcePre_n_P_1__).data().end(), 0. );
	(*L_CSDForcePre_n_P_1__) = (*L_CSDForcePreStar_n__) ;
	
}    void compute(const Loci::sequence &seq) { 
#line 312 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 312 "FSI_CSDdataloop.loci"
    }
#line 312 "FSI_CSDdataloop.loci"
} ;
#line 312 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop016_1279551148m533> register_file_FSI_CSDdataloop016_1279551148m533 ;
#line 312 "FSI_CSDdataloop.loci"
}
#line 312 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop017_1279551148m533 : public Loci::blackbox_rule {
#line 312 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<int,ublas::column_major> >  L_CSDConnectivity_ ; 
#line 312 "FSI_CSDdataloop.loci"
    Loci::blackbox<boost::multi_array<real,4> >  L_CSDnodesSys_ic_ ; 
#line 312 "FSI_CSDdataloop.loci"
public:
#line 312 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop017_1279551148m533() {
#line 312 "FSI_CSDdataloop.loci"
       name_store("CSDConnectivity",L_CSDConnectivity_) ;
#line 312 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSys_ic",L_CSDnodesSys_ic_) ;
#line 312 "FSI_CSDdataloop.loci"
       input("CSDConnectivity") ;
#line 312 "FSI_CSDdataloop.loci"
       output("CSDnodesSys_ic") ;
#line 312 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 312 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 312 "FSI_CSDdataloop.loci"
    }
#line 312 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	 size_t numEl = (*L_CSDConnectivity_).size1() ;
	// if (Loci::MPI_rank==0) cout << "ic: numEL=" << numEl << endl;
	 (*L_CSDnodesSys_ic_).resize(boost::extents[3][3][numEl][3]);
	 //(*$CSDnodesSys_ic).resize(boost::extents[3][3][3][3]);
	 std::fill( (*L_CSDnodesSys_ic_).origin(), (*L_CSDnodesSys_ic_).origin() + (*L_CSDnodesSys_ic_).size(), 0. ); // initialize with zero
	 	}    void compute(const Loci::sequence &seq) { 
#line 321 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 321 "FSI_CSDdataloop.loci"
    }
#line 321 "FSI_CSDdataloop.loci"
} ;
#line 321 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop017_1279551148m533> register_file_FSI_CSDdataloop017_1279551148m533 ;
#line 321 "FSI_CSDdataloop.loci"
}
#line 321 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop018_1279551148m534 : public Loci::blackbox_rule {
#line 321 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<int,ublas::column_major> >  L_CSDConnectivity_ ; 
#line 321 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<boost::multi_array<real,4> >  L_CSDnodesSys_ic_ ; 
#line 321 "FSI_CSDdataloop.loci"
    Loci::blackbox<boost::multi_array<real,4> >  L_CSDnodesSys_n_EQ_0__ ; 
#line 321 "FSI_CSDdataloop.loci"
public:
#line 321 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop018_1279551148m534() {
#line 321 "FSI_CSDdataloop.loci"
       name_store("CSDConnectivity",L_CSDConnectivity_) ;
#line 321 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSys_ic",L_CSDnodesSys_ic_) ;
#line 321 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSys{n=0}",L_CSDnodesSys_n_EQ_0__) ;
#line 321 "FSI_CSDdataloop.loci"
       input("CSDnodesSys_ic,CSDConnectivity") ;
#line 321 "FSI_CSDdataloop.loci"
       output("CSDnodesSys{n=0}") ;
#line 321 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 321 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 321 "FSI_CSDdataloop.loci"
    }
#line 321 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	size_t numEl = (*L_CSDConnectivity_).size1() ;
	//if (Loci::MPI_rank==0) cout << "{n=0} numEL=" << numEl << endl;
	//(*$CSDnodesSys{n=0}).resize(boost::extents[0][0][0][0]) ;
	(*L_CSDnodesSys_n_EQ_0__).resize(boost::extents[3][3][numEl][3]) ;  std::fill( (*L_CSDnodesSys_n_EQ_0__).origin(), (*L_CSDnodesSys_n_EQ_0__).origin() + (*L_CSDnodesSys_n_EQ_0__).size(), 0. );
	(*L_CSDnodesSys_n_EQ_0__) = (*L_CSDnodesSys_ic_) ;
	// std::fill( (*$CSDnodesSys{n=0}).origin(), (*$CSDnodesSys{n=0}).origin() + (*$CSDnodesSys{n=0}).size(), 0. ); // initialize with zero
}    void compute(const Loci::sequence &seq) { 
#line 331 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 331 "FSI_CSDdataloop.loci"
    }
#line 331 "FSI_CSDdataloop.loci"
} ;
#line 331 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop018_1279551148m534> register_file_FSI_CSDdataloop018_1279551148m534 ;
#line 331 "FSI_CSDdataloop.loci"
}
#line 331 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop019_1279551148m534 : public Loci::blackbox_rule {
#line 331 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<boost::multi_array<real,4> >  L_CSDnodesSys_n__ ; 
#line 331 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<boost::multi_array<real,4> >  L_CSDnodesSysStar_n__ ; 
#line 331 "FSI_CSDdataloop.loci"
    Loci::blackbox<boost::multi_array<real,4> >  L_CSDnodesSys_n_P_1__ ; 
#line 331 "FSI_CSDdataloop.loci"
public:
#line 331 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop019_1279551148m534() {
#line 331 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSys{n}",L_CSDnodesSys_n__) ;
#line 331 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSysStar{n}",L_CSDnodesSysStar_n__) ;
#line 331 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSys{n+1}",L_CSDnodesSys_n_P_1__) ;
#line 331 "FSI_CSDdataloop.loci"
       input("CSDnodesSys{n},CSDnodesSysStar{n}") ;
#line 331 "FSI_CSDdataloop.loci"
       output("CSDnodesSys{n+1}") ;
#line 331 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSINLAMS") ;
#line 331 "FSI_CSDdataloop.loci"
       conditional("iterationFinished{n,it-1}") ;
#line 331 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 331 "FSI_CSDdataloop.loci"
    }
#line 331 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	size_t numEl = (*L_CSDnodesSysStar_n__).shape()[2] ;
//	if (Loci::MPI_rank==0) cout << "Star{n,it} numEL=" << numEl << endl;
	(*L_CSDnodesSys_n_P_1__).resize(boost::extents[3][3][numEl][3]) ; std::fill( (*L_CSDnodesSys_n_P_1__).origin(), (*L_CSDnodesSys_n_P_1__).origin() + (*L_CSDnodesSys_n_P_1__).size(), 0. );
	(*L_CSDnodesSys_n_P_1__) = (*L_CSDnodesSysStar_n__) ;		
}    void compute(const Loci::sequence &seq) { 
#line 340 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 340 "FSI_CSDdataloop.loci"
    }
#line 340 "FSI_CSDdataloop.loci"
} ;
#line 340 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop019_1279551148m534> register_file_FSI_CSDdataloop019_1279551148m534 ;
#line 340 "FSI_CSDdataloop.loci"
}
#line 340 "FSI_CSDdataloop.loci"
// $type CFDIterationFinished param<bool> 
// $type stime param<real> 
// $type ncycle param<int> 

namespace {class file_FSI_CSDdataloop020_1279551148m535 : public Loci::singleton_rule {
#line 344 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_stime_n__ ; 
#line 344 "FSI_CSDdataloop.loci"
    Loci::param<real>  L_stime_nit_EQ_0__ ; 
#line 344 "FSI_CSDdataloop.loci"
public:
#line 344 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop020_1279551148m535() {
#line 344 "FSI_CSDdataloop.loci"
       name_store("stime{n}",L_stime_n__) ;
#line 344 "FSI_CSDdataloop.loci"
       name_store("stime{n,it=0}",L_stime_nit_EQ_0__) ;
#line 344 "FSI_CSDdataloop.loci"
       input("stime{n}") ;
#line 344 "FSI_CSDdataloop.loci"
       output("stime{n,it=0}") ;
#line 344 "FSI_CSDdataloop.loci"
    }
#line 344 "FSI_CSDdataloop.loci"
    void compute(const Loci::sequence &seq) { 	
//	if (Loci::MPI_rank==0) cout << "inside stime{n,it=0} in" << endl ; 
	(*L_stime_nit_EQ_0__)=(*L_stime_n__) ;
//	if (Loci::MPI_rank==0) cout << "inside stime{n,it=0} out" << endl ; 
}} ;
#line 348 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop020_1279551148m535> register_file_FSI_CSDdataloop020_1279551148m535 ;
#line 348 "FSI_CSDdataloop.loci"
}
#line 348 "FSI_CSDdataloop.loci"
	 

namespace {class file_FSI_CSDdataloop021_1279551148m535 : public Loci::singleton_rule {
#line 350 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_stime_nit__ ; 
#line 350 "FSI_CSDdataloop.loci"
    Loci::param<real>  L_stime_nit_P_1__ ; 
#line 350 "FSI_CSDdataloop.loci"
public:
#line 350 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop021_1279551148m535() {
#line 350 "FSI_CSDdataloop.loci"
       name_store("stime{n,it}",L_stime_nit__) ;
#line 350 "FSI_CSDdataloop.loci"
       name_store("stime{n,it+1}",L_stime_nit_P_1__) ;
#line 350 "FSI_CSDdataloop.loci"
       input("stime{n,it}") ;
#line 350 "FSI_CSDdataloop.loci"
       output("stime{n,it+1}") ;
#line 350 "FSI_CSDdataloop.loci"
    }
#line 350 "FSI_CSDdataloop.loci"
    void compute(const Loci::sequence &seq) { 
	(*L_stime_nit_P_1__)=(*L_stime_nit__) ;
}} ;
#line 352 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop021_1279551148m535> register_file_FSI_CSDdataloop021_1279551148m535 ;
#line 352 "FSI_CSDdataloop.loci"
}
#line 352 "FSI_CSDdataloop.loci"
	 

namespace {class file_FSI_CSDdataloop022_1279551148m535 : public Loci::singleton_rule {
#line 354 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_ncycle_n__ ; 
#line 354 "FSI_CSDdataloop.loci"
    Loci::param<int>  L_ncycle_nit_EQ_0__ ; 
#line 354 "FSI_CSDdataloop.loci"
public:
#line 354 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop022_1279551148m535() {
#line 354 "FSI_CSDdataloop.loci"
       name_store("ncycle{n}",L_ncycle_n__) ;
#line 354 "FSI_CSDdataloop.loci"
       name_store("ncycle{n,it=0}",L_ncycle_nit_EQ_0__) ;
#line 354 "FSI_CSDdataloop.loci"
       input("ncycle{n}") ;
#line 354 "FSI_CSDdataloop.loci"
       output("ncycle{n,it=0}") ;
#line 354 "FSI_CSDdataloop.loci"
    }
#line 354 "FSI_CSDdataloop.loci"
    void compute(const Loci::sequence &seq) { 	
//	if (Loci::MPI_rank==0) cout << "inside ncycle{n,it=0} in" << endl ; 
	(*L_ncycle_nit_EQ_0__)=(*L_ncycle_n__) ;
//	if (Loci::MPI_rank==0) cout << "inside ncycle{n,it=0} out" << endl ; 
}} ;
#line 358 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop022_1279551148m535> register_file_FSI_CSDdataloop022_1279551148m535 ;
#line 358 "FSI_CSDdataloop.loci"
}
#line 358 "FSI_CSDdataloop.loci"
	 

namespace {class file_FSI_CSDdataloop023_1279551148m536 : public Loci::singleton_rule {
#line 360 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_ncycle_nit__ ; 
#line 360 "FSI_CSDdataloop.loci"
    Loci::param<int>  L_ncycle_nit_P_1__ ; 
#line 360 "FSI_CSDdataloop.loci"
public:
#line 360 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop023_1279551148m536() {
#line 360 "FSI_CSDdataloop.loci"
       name_store("ncycle{n,it}",L_ncycle_nit__) ;
#line 360 "FSI_CSDdataloop.loci"
       name_store("ncycle{n,it+1}",L_ncycle_nit_P_1__) ;
#line 360 "FSI_CSDdataloop.loci"
       input("ncycle{n,it}") ;
#line 360 "FSI_CSDdataloop.loci"
       output("ncycle{n,it+1}") ;
#line 360 "FSI_CSDdataloop.loci"
    }
#line 360 "FSI_CSDdataloop.loci"
    void compute(const Loci::sequence &seq) { 
	(*L_ncycle_nit_P_1__)=(*L_ncycle_nit__) ;
}} ;
#line 362 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop023_1279551148m536> register_file_FSI_CSDdataloop023_1279551148m536 ;
#line 362 "FSI_CSDdataloop.loci"
}
#line 362 "FSI_CSDdataloop.loci"
	 

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
namespace {class file_FSI_CSDdataloop024_1279551148m537 : public Loci::blackbox_rule {
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDstartingTimeStep_ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDtipNode_ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDdimension_ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<int,ublas::column_major> >  L_CSDConnectivity_ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<vector<int> >  L_CSDBCdof_ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<vector<real> >  L_CSDBCZeroConstraint_ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDisp_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVel_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAcc_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePre_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<boost::multi_array<real,4> >  L_CSDnodesSys_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_stime_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_ncycle_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForce_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_timeStep_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDE1_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDE2_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDnu12_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDnu21_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDG12_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDthicknessStructure_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDrhoStructure_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDintegrationScheme_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDdelta_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDswitchStiffening_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDgenAlphaCoeff_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDnewmarkGammaCoeff_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDnewmarkBetaCoeff_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDdampingCoeff1_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDdampingCoeff2_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDflappingType_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDplungingType_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDfrequency_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDplungeAmplitudeX_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDplungeAmplitudeY_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDplungeAmplitudeZ_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDflappingAmplitudeX_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDflappingAmplitudeY_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDflappingAmplitudeZ_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_itfsi_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDdisplacementsStar_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesDispStar_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesVelStar_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodesAccStar_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePreStar_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
    Loci::blackbox<boost::multi_array<real,4> >  L_CSDnodesSysStar_n__ ; 
#line 386 "FSI_CSDdataloop.loci"
public:
#line 386 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop024_1279551148m537() {
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDstartingTimeStep",L_CSDstartingTimeStep_) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDtipNode",L_CSDtipNode_) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDdimension",L_CSDdimension_) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDConnectivity",L_CSDConnectivity_) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDBCdof",L_CSDBCdof_) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDBCZeroConstraint",L_CSDBCZeroConstraint_) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDnodes{n}",L_CSDnodes_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDdisplacementsStar{n}",L_CSDdisplacementsStar_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDisp{n}",L_CSDnodesDisp_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDnodesDispStar{n}",L_CSDnodesDispStar_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVel{n}",L_CSDnodesVel_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDnodesVelStar{n}",L_CSDnodesVelStar_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAcc{n}",L_CSDnodesAcc_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDnodesAccStar{n}",L_CSDnodesAccStar_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDForcePre{n}",L_CSDForcePre_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDForcePreStar{n}",L_CSDForcePreStar_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSys{n}",L_CSDnodesSys_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDnodesSysStar{n}",L_CSDnodesSysStar_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("stime{n}",L_stime_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("ncycle{n}",L_ncycle_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDForce{n}",L_CSDForce_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("timeStep{n}",L_timeStep_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDE1{n}",L_CSDE1_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDE2{n}",L_CSDE2_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDnu12{n}",L_CSDnu12_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDnu21{n}",L_CSDnu21_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDG12{n}",L_CSDG12_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDthicknessStructure{n}",L_CSDthicknessStructure_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDrhoStructure{n}",L_CSDrhoStructure_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDintegrationScheme{n}",L_CSDintegrationScheme_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDdelta{n}",L_CSDdelta_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDswitchStiffening{n}",L_CSDswitchStiffening_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDgenAlphaCoeff{n}",L_CSDgenAlphaCoeff_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDnewmarkGammaCoeff{n}",L_CSDnewmarkGammaCoeff_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDnewmarkBetaCoeff{n}",L_CSDnewmarkBetaCoeff_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDdampingCoeff1{n}",L_CSDdampingCoeff1_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDdampingCoeff2{n}",L_CSDdampingCoeff2_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDflappingType{n}",L_CSDflappingType_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDplungingType{n}",L_CSDplungingType_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDfrequency{n}",L_CSDfrequency_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDplungeAmplitudeX{n}",L_CSDplungeAmplitudeX_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDplungeAmplitudeY{n}",L_CSDplungeAmplitudeY_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDplungeAmplitudeZ{n}",L_CSDplungeAmplitudeZ_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDflappingAmplitudeX{n}",L_CSDflappingAmplitudeX_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDflappingAmplitudeY{n}",L_CSDflappingAmplitudeY_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("CSDflappingAmplitudeZ{n}",L_CSDflappingAmplitudeZ_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       name_store("itfsi{n}",L_itfsi_n__) ;
#line 386 "FSI_CSDdataloop.loci"
       input("							CSDnodes_ic,CSDnodes{n},CSDnodesDisp{n},CSDnodesVel{n},CSDnodesAcc{n},CSDForce{n},CSDForcePre{n},CSDConnectivity,CSDBCdof,CSDBCZeroConstraint,CSDnodesSys{n},							stime{n},ncycle{n},timeStep{n},CSDstartingTimeStep,CSDtipNode,							CSDE1{n},CSDE2{n},CSDnu12{n},CSDnu21{n},CSDG12{n},CSDthicknessStructure{n},CSDrhoStructure{n},							CSDintegrationScheme{n},CSDdelta{n},CSDswitchStiffening{n},CSDgenAlphaCoeff{n},CSDnewmarkGammaCoeff{n},CSDnewmarkBetaCoeff{n},CSDdampingCoeff1{n},CSDdampingCoeff2{n},							CSDflappingType{n},CSDplungingType{n},							CSDfrequency{n},CSDplungeAmplitudeX{n},CSDplungeAmplitudeY{n},CSDplungeAmplitudeZ{n},CSDflappingAmplitudeX{n},CSDflappingAmplitudeY{n},CSDflappingAmplitudeZ{n},							itfsi{n},CSDdimension") ;
#line 386 "FSI_CSDdataloop.loci"
       output("CSDdisplacementsStar{n}") ;
#line 386 "FSI_CSDdataloop.loci"
       output("CSDnodesDispStar{n}") ;
#line 386 "FSI_CSDdataloop.loci"
       output("CSDnodesVelStar{n}") ;
#line 386 "FSI_CSDdataloop.loci"
       output("CSDnodesAccStar{n}") ;
#line 386 "FSI_CSDdataloop.loci"
       output("CSDForcePreStar{n}") ;
#line 386 "FSI_CSDdataloop.loci"
       output("CSDnodesSysStar{n}") ;
#line 386 "FSI_CSDdataloop.loci"
       constraint("FSICoupling{n},FSINLAMS{n}") ;
#line 386 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 386 "FSI_CSDdataloop.loci"
    }
#line 386 "FSI_CSDdataloop.loci"
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
	(*L_CSDdisplacementsStar_n__).resize(CSDNumNodes,3) ; std::fill( (*L_CSDdisplacementsStar_n__).data().begin(),(*L_CSDdisplacementsStar_n__).data().end(), 0. );
//	(*$CSDdisplacementsOldStar).resize(CSDNumNodes,3) ; (*$CSDdisplacementsOldStar).clear();
	(*L_CSDnodesDispStar_n__).resize(CSDAnswerSize,1) ; std::fill( (*L_CSDnodesDispStar_n__).data().begin(),(*L_CSDnodesDispStar_n__).data().end(), 0. );
	(*L_CSDnodesVelStar_n__).resize(CSDAnswerSize,1) ; std::fill( (*L_CSDnodesVelStar_n__).data().begin(),(*L_CSDnodesVelStar_n__).data().end(), 0. );
	(*L_CSDnodesAccStar_n__).resize(CSDAnswerSize,1) ; std::fill( (*L_CSDnodesAccStar_n__).data().begin(),(*L_CSDnodesAccStar_n__).data().end(), 0. );
	(*L_CSDForcePreStar_n__).resize(CSDAnswerSize,1) ; std::fill( (*L_CSDForcePreStar_n__).data().begin(),(*L_CSDForcePreStar_n__).data().end(), 0. );
	(*L_CSDnodesSysStar_n__).resize(boost::extents[3][3][CSDNumElems][3]); std::fill( (*L_CSDnodesSysStar_n__).origin(), (*L_CSDnodesSysStar_n__).origin() + (*L_CSDnodesSysStar_n__).size(), 0. );
		
	
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
	if (rank==0) cout << "[I]... timestep = " << *L_timeStep_n__ << endl ;	
	if (rank==0) cout << "[I]... tip node = " << *L_CSDtipNode_<< ": " << (*L_CSDnodes_n__)(*L_CSDtipNode_,0) << ", "  << (*L_CSDnodes_n__)(*L_CSDtipNode_,1) << ", "<< (*L_CSDnodes_n__)(*L_CSDtipNode_,2) << endl ;	
	if (rank==0) cout << "[I]... CSDE1 = " << *L_CSDE1_n__ << endl ;					
	if (rank==0) cout << "[I]... CSDE2 = " << *L_CSDE2_n__ << endl ;		
	if (rank==0) cout << "[I]... CSDnu12 = " << *L_CSDnu12_n__ << endl ;		
	if (rank==0) cout << "[I]... CSDnu21 = " << *L_CSDnu21_n__ << endl ;		
	if (rank==0) cout << "[I]... CSDG12 = " << *L_CSDG12_n__ << endl ;		
	if (rank==0) cout << "[I]... CSDthicknessStructure = " << *L_CSDthicknessStructure_n__ << endl ;		
	if (rank==0) cout << "[I]... CSDrhoStructure = " << *L_CSDrhoStructure_n__ << endl ;		
	if (rank==0) cout << "[I]... CSDdelta = " << *L_CSDdelta_n__ << endl ;		
	//if (rank==0) cout << "[I]... CSDexcitationType = " << *$CSDexcitationType{n} << endl ;		
	if (rank==0) cout << "[I]... CSDfrequency = " << *L_CSDfrequency_n__ << endl ;		

	if (rank==0) {
		std::stringstream filename ;	
		filename << "CSDnodesforces" << tempTimeStepNumber << "it" << *L_itfsi_n__ << ".dat" ;
		ofstream CSDnodesforces ;
		CSDnodesforces.open(filename.str().c_str(), ios::out) ;
			CSDnodesforces << "CSDnodes.x" << ", " << "CSDnodes.y" << ", " << "CSDnodes.z" << ", " << "CSDForce.x" << ", "  << "CSDForce.y" << ", " << "CSDForce.z" << "CSDdisp.x" << ", " << "CSDdisp.y" << ", " << "CSDdisp.z" << endl ;
		for(int i=0; i<CSDNumNodes;++i) { // in CFD coordinates
			CSDnodesforces << (*L_CSDnodes_n__)(i,0) << ", " << (*L_CSDnodes_n__)(i,1) << ", " << (*L_CSDnodes_n__)(i,2) << ", " << (*L_CSDForce_n__)(i*6+0,0) << ", "  << (*L_CSDForce_n__)(i*6+1,0) << ", " << (*L_CSDForce_n__)(i*6+2,0)  << (*L_CSDnodes_n__)(i,0) - (*L_CSDnodes_ic_)(i,0)<< ", " << (*L_CSDnodes_n__)(i,1) - (*L_CSDnodes_ic_)(i,1) << ", " << (*L_CSDnodes_n__)(i,2) - (*L_CSDnodes_ic_)(i,2) << endl ;
		}
		CSDnodesforces.close();
	}
	
	// Displacements at the previous itfis
	//(CSDdisplacementsOld) = (*$CSDnodes) - (*$CSDnodes_ic) ;
	
//	if (rank==0) cout << "connectivity matrix in loci-stream: " << (*$CSDConnectivity) << endl ;

//	if (rank==0) cout << "force matrix in loci-stream: " << (*$CSDForce) << endl ;

  int itfsiInt = 0;

      for(int i=0; i<CSDNumNodes; ++i) {
		 if (rank==0) cout << "displacementsStarBefore: " << (*L_CSDdisplacementsStar_n__)(i,0) << endl ;
  		 if (rank==0) cout << "displacementsStarBefore: " << (*L_CSDdisplacementsStar_n__)(i,1) << endl ;
		 if (rank==0) cout << "displacementsStarBefore: " << (*L_CSDdisplacementsStar_n__)(i,2) << endl ;
  	}		
	
  pass_(&rank, &itfsiInt, &CSDNumNodes, &(*L_CSDnodes_ic_)(0,0), &(*L_CSDnodes_ic_)(0,1),&(*L_CSDnodes_ic_)(0,2), 
  			&CSDNumElems, &(*L_CSDConnectivity_)(0,0), 
  			&CSDNumBC, &(*L_CSDBCdof_)[0], &(*L_CSDBCZeroConstraint_)[0], 
  			&(*L_CSDE1_n__), &(*L_CSDnu12_n__), &(*L_CSDrhoStructure_n__), &(*L_CSDthicknessStructure_n__), 
  			&(*L_CSDintegrationScheme_n__), &(*L_CSDnewmarkBetaCoeff_n__), &(*L_CSDnewmarkGammaCoeff_n__), &(*L_CSDgenAlphaCoeff_n__), 
  			&CSDAnswerSize, &(*L_CSDnodesDisp_n__)(0,0), &(*L_CSDnodesVel_n__)(0,0), &(*L_CSDnodesAcc_n__)(0,0), &(*L_CSDForce_n__)(0,0), &(*L_CSDForcePre_n__)(0,0), (*L_CSDnodesSys_n__).data(),
  			&(*L_timeStep_n__), &tempTimeStepNumber, &(*L_CSDtipNode_), // &(*$ncycle), 
  			&(*L_CSDflappingType_n__), &(*L_CSDplungingType_n__), 
  			&(*L_CSDfrequency_n__), &(*L_CSDflappingAmplitudeX_n__), &(*L_CSDflappingAmplitudeY_n__), &(*L_CSDflappingAmplitudeZ_n__), &(*L_CSDplungeAmplitudeX_n__), &(*L_CSDplungeAmplitudeY_n__), &(*L_CSDplungeAmplitudeZ_n__), 
  			&(*L_CSDdisplacementsStar_n__)(0,0), &(*L_CSDnodesDispStar_n__)(0,0), &(*L_CSDnodesVelStar_n__)(0,0), &(*L_CSDnodesAccStar_n__)(0,0), (*L_CSDnodesSysStar_n__).data() );
  
  if (rank==0) cout << "[I] Communication with NLAMS finished.. " << endl ;
  
  // Updates
  // CSDForce -> CSDForcePre
  (*L_CSDForcePreStar_n__) = (*L_CSDForce_n__ ) ;
  
  if (*L_CSDdimension_== 2) { // displacement from NLAMS is sometimes not exactly zero due to floating point errors
  	for(int i=0; i<CSDNumNodes; ++i) {
  		 (*L_CSDdisplacementsStar_n__)(i,2) = 0. ;
  	}
  }
  //-----------
  for(int i=0; i<CSDNumNodes; ++i) {
//	  (*$CSDdisplacementsStar{n})(i,0) = 0.0 ;
//	  (*$CSDdisplacementsStar{n})(i,1) = 0.0 ;
	  (*L_CSDdisplacementsStar_n__)(i,2) = 0.0 ;
	}
  //-----------
  
  			
	for (int i=0; i<3; ++i) {
		for (int n=0; n<CSDNumNodes; ++n) {
			if ( fabs((*L_CSDdisplacementsStar_n__)(n,i))  < 1.0e-5 ) (*L_CSDdisplacementsStar_n__)(n,i) = 0.0 ;
		}
	}
  
  
    for(int i=0; i<CSDNumNodes; ++i) {
		 if (rank==0) cout << "displacementsStar: " << (*L_CSDdisplacementsStar_n__)(i,0) << endl ;
  		 if (rank==0) cout << "displacementsStar: " << (*L_CSDdisplacementsStar_n__)(i,1) << endl ;
		 if (rank==0) cout << "displacementsStar: " << (*L_CSDdisplacementsStar_n__)(i,2) << endl ;
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
//		 	  (*$CSDdisplacementsStar{n})(i,0) = (*$CSDnodes_ic)(i,0) * ( cos(alpha*PI/180) - 1.) ;
//		 	  (*$CSDdisplacementsStar{n})(i,1) = (*$CSDnodes_ic)(i,0) * ( sin(alpha*PI/180) ) ;
//		 	  (*$CSDdisplacementsStar{n})(i,2) = 0.0 ;
//		  }
//	   	if (Loci::MPI_rank==0) cout << "CSDdisplacement:" << (*$stime{n}) << ", alpha = " << alpha << endl ;

  
 // (*$CSDdisplacementsOldStar) = (*$CSDForce ) ;
  
  
  
 } // end if CSDstartingTimeStep
 
//}

}    void compute(const Loci::sequence &seq) { 
#line 532 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 532 "FSI_CSDdataloop.loci"
    }
#line 532 "FSI_CSDdataloop.loci"
} ;
#line 532 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop024_1279551148m537> register_file_FSI_CSDdataloop024_1279551148m537 ;
#line 532 "FSI_CSDdataloop.loci"
}
#line 532 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop025_1279551148m540 : public Loci::singleton_rule {
#line 532 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerAxis_ ; 
#line 532 "FSI_CSDdataloop.loci"
    Loci::param<int>  L_CSDEulerBeamDirection_ ; 
#line 532 "FSI_CSDdataloop.loci"
    Loci::param<int>  L_CSDEulerSpanDirection_ ; 
#line 532 "FSI_CSDdataloop.loci"
public:
#line 532 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop025_1279551148m540() {
#line 532 "FSI_CSDdataloop.loci"
       name_store("CSDEulerAxis",L_CSDEulerAxis_) ;
#line 532 "FSI_CSDdataloop.loci"
       name_store("CSDEulerBeamDirection",L_CSDEulerBeamDirection_) ;
#line 532 "FSI_CSDdataloop.loci"
       name_store("CSDEulerSpanDirection",L_CSDEulerSpanDirection_) ;
#line 532 "FSI_CSDdataloop.loci"
       input("CSDEulerAxis") ;
#line 532 "FSI_CSDdataloop.loci"
       output("CSDEulerBeamDirection") ;
#line 532 "FSI_CSDdataloop.loci"
       output("CSDEulerSpanDirection") ;
#line 532 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 532 "FSI_CSDdataloop.loci"
    }
#line 532 "FSI_CSDdataloop.loci"
    void compute(const Loci::sequence &seq) { 
	(*L_CSDEulerBeamDirection_)= (*L_CSDEulerAxis_);
	(*L_CSDEulerSpanDirection_)= 3 - ((*L_CSDEulerAxis_)+ 1) ; // 3 % 1 [x] = 2 [z], 3 % 3 [z] = 0 [x]
//	if (rank==0) cout << "[I] Euler-Bernoulli.. $CSDEulerBeamDirection: " << $CSDEulerBeamDirection << ", $CSDEulerSpanDirection" << $CSDEulerSpanDirection << endl ;
}} ;
#line 536 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop025_1279551148m540> register_file_FSI_CSDdataloop025_1279551148m540 ;
#line 536 "FSI_CSDdataloop.loci"
}
#line 536 "FSI_CSDdataloop.loci"

   
namespace {class file_FSI_CSDdataloop026_1279551148m541 : public Loci::blackbox_rule {
#line 538 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDdimension_ ; 
#line 538 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSD2dSpanCenter_ ; 
#line 538 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDEulerXstart_ ; 
#line 538 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDEulerXend_ ; 
#line 538 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerXnum_ ; 
#line 538 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerBeamDirection_ ; 
#line 538 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerSpanDirection_ ; 
#line 538 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 538 "FSI_CSDdataloop.loci"
public:
#line 538 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop026_1279551148m541() {
#line 538 "FSI_CSDdataloop.loci"
       name_store("CSDdimension",L_CSDdimension_) ;
#line 538 "FSI_CSDdataloop.loci"
       name_store("CSD2dSpanCenter",L_CSD2dSpanCenter_) ;
#line 538 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 538 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXstart",L_CSDEulerXstart_) ;
#line 538 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXend",L_CSDEulerXend_) ;
#line 538 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXnum",L_CSDEulerXnum_) ;
#line 538 "FSI_CSDdataloop.loci"
       name_store("CSDEulerBeamDirection",L_CSDEulerBeamDirection_) ;
#line 538 "FSI_CSDdataloop.loci"
       name_store("CSDEulerSpanDirection",L_CSDEulerSpanDirection_) ;
#line 538 "FSI_CSDdataloop.loci"
       input("CSDEulerXstart,CSDEulerXend,CSDEulerXnum,CSD2dSpanCenter,CSDdimension,CSDEulerBeamDirection,CSDEulerSpanDirection") ;
#line 538 "FSI_CSDdataloop.loci"
       output("CSDnodes_ic") ;
#line 538 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 538 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 538 "FSI_CSDdataloop.loci"
    }
#line 538 "FSI_CSDdataloop.loci"
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
#line 559 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 559 "FSI_CSDdataloop.loci"
    }
#line 559 "FSI_CSDdataloop.loci"
} ;
#line 559 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop026_1279551148m541> register_file_FSI_CSDdataloop026_1279551148m541 ;
#line 559 "FSI_CSDdataloop.loci"
}
#line 559 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop027_1279551148m542 : public Loci::blackbox_rule {
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDE1_ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDrhoStructure_ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDthicknessStructure_ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDintegrationScheme_ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDfrequency_ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDgenAlphaCoeff_ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDnewmarkGammaCoeff_ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDnewmarkBetaCoeff_ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDplungeAmplitudeY_ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDflappingAmplitudeX_ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDEulerChord_ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerXnum_ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerBeamDirection_ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerSpanDirection_ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePre_n__ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_stime_n__ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_ncycle_n__ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_timeStep_ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForce_n__ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDx_n__ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxdot_n__ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxddot_n__ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDstartingTimeStep_n__ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDdisplacementsStar_n__ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForcePreStar_n__ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxStar_n__ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxdotStar_n__ ; 
#line 565 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxddotStar_n__ ; 
#line 565 "FSI_CSDdataloop.loci"
public:
#line 565 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop027_1279551148m542() {
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDE1",L_CSDE1_) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDrhoStructure",L_CSDrhoStructure_) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDthicknessStructure",L_CSDthicknessStructure_) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDintegrationScheme",L_CSDintegrationScheme_) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDfrequency",L_CSDfrequency_) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDgenAlphaCoeff",L_CSDgenAlphaCoeff_) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDnewmarkGammaCoeff",L_CSDnewmarkGammaCoeff_) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDnewmarkBetaCoeff",L_CSDnewmarkBetaCoeff_) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDplungeAmplitudeY",L_CSDplungeAmplitudeY_) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDflappingAmplitudeX",L_CSDflappingAmplitudeX_) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDEulerChord",L_CSDEulerChord_) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXnum",L_CSDEulerXnum_) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDEulerBeamDirection",L_CSDEulerBeamDirection_) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDEulerSpanDirection",L_CSDEulerSpanDirection_) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDdisplacementsStar{n}",L_CSDdisplacementsStar_n__) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDForcePre{n}",L_CSDForcePre_n__) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDForcePreStar{n}",L_CSDForcePreStar_n__) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("stime{n}",L_stime_n__) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("ncycle{n}",L_ncycle_n__) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("timeStep",L_timeStep_) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDForce{n}",L_CSDForce_n__) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDx{n}",L_CSDx_n__) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDxdot{n}",L_CSDxdot_n__) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDxddot{n}",L_CSDxddot_n__) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDstartingTimeStep{n}",L_CSDstartingTimeStep_n__) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDxStar{n}",L_CSDxStar_n__) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDxdotStar{n}",L_CSDxdotStar_n__) ;
#line 565 "FSI_CSDdataloop.loci"
       name_store("CSDxddotStar{n}",L_CSDxddotStar_n__) ;
#line 565 "FSI_CSDdataloop.loci"
       input("							CSDEulerXnum,CSDnodes_ic,CSDx{n},CSDxdot{n},CSDxddot{n},CSDForce{n},CSDForcePre{n},							stime{n},ncycle{n},timeStep,CSDstartingTimeStep{n},CSDEulerBeamDirection,CSDEulerSpanDirection,							CSDE1,CSDthicknessStructure,CSDrhoStructure,CSDEulerChord,							CSDintegrationScheme,CSDgenAlphaCoeff,CSDnewmarkGammaCoeff,CSDnewmarkBetaCoeff,							CSDfrequency,CSDplungeAmplitudeY,CSDflappingAmplitudeX") ;
#line 565 "FSI_CSDdataloop.loci"
       output("CSDdisplacementsStar{n}") ;
#line 565 "FSI_CSDdataloop.loci"
       output("CSDForcePreStar{n}") ;
#line 565 "FSI_CSDdataloop.loci"
       output("CSDxStar{n}") ;
#line 565 "FSI_CSDdataloop.loci"
       output("CSDxdotStar{n}") ;
#line 565 "FSI_CSDdataloop.loci"
       output("CSDxddotStar{n}") ;
#line 565 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 565 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 565 "FSI_CSDdataloop.loci"
    }
#line 565 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) { 
	// call NLAMS
	//communicateWithNLAMS(*CSDdisplacementsStar, ....) ;
	//if (*$CFDIterationFinished) {
	
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
	
	// Initialize
	(*L_CSDdisplacementsStar_n__).resize(CSDNumNodes,3) ; std::fill( (*L_CSDdisplacementsStar_n__).data().begin(),(*L_CSDdisplacementsStar_n__).data().end(), 0. );	
	(*L_CSDxStar_n__).resize(CSDdofFree) ; std::fill( (*L_CSDxStar_n__).begin(),(*L_CSDxStar_n__).end(), 0. );
	(*L_CSDxdotStar_n__).resize(CSDdofFree) ; std::fill( (*L_CSDxdotStar_n__).begin(),(*L_CSDxdotStar_n__).end(), 0. );
	(*L_CSDxddotStar_n__).resize(CSDdofFree) ; std::fill( (*L_CSDxddotStar_n__).begin(),(*L_CSDxddotStar_n__).end(), 0. );
	(*L_CSDForcePreStar_n__).resize(CSDAnswerSize,1) ; std::fill( (*L_CSDForcePreStar_n__).data().begin(),(*L_CSDForcePreStar_n__).data().end(), 0. );


	if (rank==0) cout << "[I] Euler-Bernoulli.. initialized " << endl ;

//	if (*$CFDIterationFinished{n,it-1}) {	

		if (*L_ncycle_n__ < *L_CSDstartingTimeStep_n__) {
			 if (rank==0) cout << "CSDstartingTime not yet reached. Current time step = " << *L_ncycle_n__ << ", CSDstartingTimeStep = " << *L_CSDstartingTimeStep_n__ << endl ;
		} else {	
			int tempTimeStepNumber = *L_ncycle_n__ - *L_CSDstartingTimeStep_n__ + 1 ;	
				
			ublas::vector<real> CSDbeam1dForce(CSDNumNodes) ;	std::fill(CSDbeam1dForce.begin(), CSDbeam1dForce.end(), 0.0) ;
			ublas::vector<real> disp(CSDNumNodes) ; std::fill(disp.begin(), disp.end(), 0.0) ;
			ublas::vector<real> nodes(CSDNumNodes) ; std::fill(nodes.begin(), nodes.end(), 0.0) ;
			
			if (rank==0) cout << "[I] Euler-Bernoulli.. communication starting: " << CSDNumNodes << ", " << CSDbeam1dForce.size() << ", " << (*L_CSDForce_n__).size1() << endl ;	
			if (rank==0) cout << "[I] Euler-Bernoulli.. beamdirection: " << beamdirection << ", spandirection" << spandirection << endl ;	

			for (int i=0; i<CSDNumNodes; ++i) {
	//			if (rank==0) cout << "[I] Euler-Bernoulli.. before CSDforce" << endl ;	
	//			CSDbeam1dForce[i] = (*$CSDForce{n})(6*i+1,0) ; // y-component only
				CSDbeam1dForce[i] = ( 2.0 * ((*L_CSDForce_n__)(6*i+1,0)) - (*L_CSDForcePre_n__)(6*i+1,0) ) ; // y-component only
	//			CSDbeam1dForce[i] = 0.0 ; // y-component only			
	//			if (rank==0) cout << "[I] Euler-Bernoulli.. CSDbeam1dForce[" << i << "] = " << CSDbeam1dForce[i] << endl ;	
				nodes[i] = (*L_CSDnodes_ic_)(i,beamdirection) ;
	//			if (rank==0) cout << "[I] Euler-Bernoulli.. nodes[" << i << "] = " << nodes[i] << endl ;	
			}
			
//			for (int i=0; i<CSDdofFree; ++i) {
//				cout << "rank = " << rank << ", " << "[I] Euler-Bernoulli: i, x, xdot, xddot, disp = " << i << "\t" << (*$CSDx{n})[i] << "\t" << (*$CSDxdot{n})[i] << "\t" << (*$CSDxddot{n})[i] << "\t" << endl ;
//			}
			
		if (rank==0) {
			std::stringstream filename ;	
			filename << "CSDnodesforces" << setfill('0') << setw(5) << tempTimeStepNumber << ".dat" ;
			ofstream CSDnodesforces ;
			CSDnodesforces.open(filename.str().c_str(), ios::out) ;
			CSDnodesforces << "CSDnodes.x" << ", " << "CSDnodes.y" << ", " << "CSDnodes.z" << ", " << "CSDForce.x" << ", "  << "CSDForce.y" << ", " << "CSDForce.z" << "CSDdisp.x" << ", " << "CSDdisp.y" << ", " << "CSDdisp.z" << endl ;
			for(int i=0; i<CSDNumNodes;++i) { // in CFD coordinates
				CSDnodesforces << nodes[i] << ", " << (*L_CSDForce_n__)(i*6+0,0) << ", "  << (*L_CSDForce_n__)(i*6+1,0) << ", " << (*L_CSDForce_n__)(i*6+2,0)  << endl ;
			}
			CSDnodesforces.close();
		}
	

			if (rank==0) cout << "[I] Euler-Bernoulli.. communication starting " << endl ;

			// communiate	
	     beam1d_(&rank, &CSDNumNodes, &CSDdofFree, &(nodes)[0], &(*L_CSDplungeAmplitudeY_), &(*L_CSDflappingAmplitudeX_), &(*L_CSDfrequency_), &(*L_CSDEulerChord_), &(*L_CSDrhoStructure_), &(*L_CSDthicknessStructure_),
	        &*L_CSDE1_, &(*L_CSDx_n__)[0], &(*L_CSDxdot_n__)[0], &(*L_CSDxddot_n__)[0], &(CSDbeam1dForce)[0], &(*L_timeStep_), &(*L_stime_n__), &tempTimeStepNumber, 
	        &(*L_CSDxStar_n__)[0], &(*L_CSDxdotStar_n__)[0], &(*L_CSDxddotStar_n__)[0], &(disp)[0]) ;

	  
	  
		  // update
		  for (int i=0; i<CSDNumNodes; ++i) {
				//(*$CSDdisplacementsStar{n})(i,0) = 0.0 ;
				(*L_CSDdisplacementsStar_n__)(i,1) = disp[i];
				//(*$CSDdisplacementsStar{n})(i,2) = 0.0 ;
	//			cout << "rank = " << rank << ", " << "[I] Euler-Bernoulli: i, CSDdisplacementsStar, disp = " << i << "\t" << (*$CSDdisplacementsStar{n})(i,1) << "\t" << disp[i] << endl ;				
			}
		(*L_CSDForcePreStar_n__) = (*L_CSDForce_n__) ;
			
			for (int i=0; i<3; ++i) {
				for (int n=0; n<CSDNumNodes; ++n) {
	//				if ( fabs((*$CSDdisplacementsStar{n})(n,i))  < 1.0e-5 ) (*$CSDdisplacementsStar{n})(n,i) = 0.0 ;
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
//		 	  (*$CSDdisplacementsStar{n})(i,0) = (*$CSDnodes_ic)(i,0) * ( cos(alpha*PI/180) - 1.) ;
////		 	  (*$CSDdisplacementsStar{n})(i,0) = 0.0 ;
//		 	  (*$CSDdisplacementsStar{n})(i,1) = (*$CSDnodes_ic)(i,0) * ( sin(alpha*PI/180) ) ;
//		 	  (*$CSDdisplacementsStar{n})(i,2) = 0.0 ;
//		  }
//	   	if (Loci::MPI_rank==0) cout << "CSDdisplacement:" << (*$stime{n}) << ", alpha = " << alpha << endl ;			
//			
	//	for (int i=0; i<CSDdofFree; ++i) {
	//		if (rank==0) cout << "[I] Euler-Bernoulli: i, xStar, xdotStar, xddotStar = " << i << "\t" << (*$CSDxStar{n})[i] << "\t" << (*$CSDxdotStar{n})[i] << "\t" << (*$CSDxddotStar{n})[i] << "\t" << endl ;
	//	}
			
		if (rank==0) cout << "[I] Communication with the Euler-Bernoulli beam solver done. " << endl ;
	}   
//	} else {
//		if (rank==0) cout << "[I] Communication with the Euler-Bernoulli beam solver bypassed. " << endl ;
//		}
}    void compute(const Loci::sequence &seq) { 
#line 703 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 703 "FSI_CSDdataloop.loci"
    }
#line 703 "FSI_CSDdataloop.loci"
} ;
#line 703 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop027_1279551148m542> register_file_FSI_CSDdataloop027_1279551148m542 ;
#line 703 "FSI_CSDdataloop.loci"
}
#line 703 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop028_1279551148m544 : public Loci::blackbox_rule {
#line 703 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerXnum_ ; 
#line 703 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDx_ic_ ; 
#line 703 "FSI_CSDdataloop.loci"
public:
#line 703 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop028_1279551148m544() {
#line 703 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXnum",L_CSDEulerXnum_) ;
#line 703 "FSI_CSDdataloop.loci"
       name_store("CSDx_ic",L_CSDx_ic_) ;
#line 703 "FSI_CSDdataloop.loci"
       input("CSDEulerXnum") ;
#line 703 "FSI_CSDdataloop.loci"
       output("CSDx_ic") ;
#line 703 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 703 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 703 "FSI_CSDdataloop.loci"
    }
#line 703 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) { 
	int CSDNumBC = 2 ;
	int CSDdofNode = 2 ;
	int CSDdofFree = *L_CSDEulerXnum_* CSDdofNode - CSDNumBC ;
	
	(*L_CSDx_ic_).resize(CSDdofFree) ; std::fill( (*L_CSDx_ic_).begin(),(*L_CSDx_ic_).end(), 0. );
}    void compute(const Loci::sequence &seq) { 
#line 711 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 711 "FSI_CSDdataloop.loci"
    }
#line 711 "FSI_CSDdataloop.loci"
} ;
#line 711 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop028_1279551148m544> register_file_FSI_CSDdataloop028_1279551148m544 ;
#line 711 "FSI_CSDdataloop.loci"
}
#line 711 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop029_1279551148m544 : public Loci::blackbox_rule {
#line 711 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDx_ic_ ; 
#line 711 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDx_n_EQ_0__ ; 
#line 711 "FSI_CSDdataloop.loci"
public:
#line 711 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop029_1279551148m544() {
#line 711 "FSI_CSDdataloop.loci"
       name_store("CSDx_ic",L_CSDx_ic_) ;
#line 711 "FSI_CSDdataloop.loci"
       name_store("CSDx{n=0}",L_CSDx_n_EQ_0__) ;
#line 711 "FSI_CSDdataloop.loci"
       input("CSDx_ic") ;
#line 711 "FSI_CSDdataloop.loci"
       output("CSDx{n=0}") ;
#line 711 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 711 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 711 "FSI_CSDdataloop.loci"
    }
#line 711 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDx_n_EQ_0__).resize((*L_CSDx_ic_).size()) ; std::fill( (*L_CSDx_n_EQ_0__).begin(),(*L_CSDx_n_EQ_0__).end(), 0. );
	(*L_CSDx_n_EQ_0__) = (*L_CSDx_ic_) ;
}    void compute(const Loci::sequence &seq) { 
#line 717 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 717 "FSI_CSDdataloop.loci"
    }
#line 717 "FSI_CSDdataloop.loci"
} ;
#line 717 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop029_1279551148m544> register_file_FSI_CSDdataloop029_1279551148m544 ;
#line 717 "FSI_CSDdataloop.loci"
}
#line 717 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop030_1279551148m545 : public Loci::blackbox_rule {
#line 717 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDx_n__ ; 
#line 717 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxStar_n__ ; 
#line 717 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDx_n_P_1__ ; 
#line 717 "FSI_CSDdataloop.loci"
public:
#line 717 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop030_1279551148m545() {
#line 717 "FSI_CSDdataloop.loci"
       name_store("CSDx{n}",L_CSDx_n__) ;
#line 717 "FSI_CSDdataloop.loci"
       name_store("CSDxStar{n}",L_CSDxStar_n__) ;
#line 717 "FSI_CSDdataloop.loci"
       name_store("CSDx{n+1}",L_CSDx_n_P_1__) ;
#line 717 "FSI_CSDdataloop.loci"
       input("CSDx{n},CSDxStar{n}") ;
#line 717 "FSI_CSDdataloop.loci"
       output("CSDx{n+1}") ;
#line 717 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 717 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 717 "FSI_CSDdataloop.loci"
    }
#line 717 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDx_n_P_1__).resize((*L_CSDx_n__).size()) ;std::fill( (*L_CSDx_n_P_1__).begin(),(*L_CSDx_n_P_1__).end(), 0. );
	(*L_CSDx_n_P_1__) = (*L_CSDxStar_n__) ;	
}    void compute(const Loci::sequence &seq) { 
#line 723 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 723 "FSI_CSDdataloop.loci"
    }
#line 723 "FSI_CSDdataloop.loci"
} ;
#line 723 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop030_1279551148m545> register_file_FSI_CSDdataloop030_1279551148m545 ;
#line 723 "FSI_CSDdataloop.loci"
}
#line 723 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop031_1279551148m545 : public Loci::blackbox_rule {
#line 723 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerXnum_ ; 
#line 723 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxdot_ic_ ; 
#line 723 "FSI_CSDdataloop.loci"
public:
#line 723 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop031_1279551148m545() {
#line 723 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXnum",L_CSDEulerXnum_) ;
#line 723 "FSI_CSDdataloop.loci"
       name_store("CSDxdot_ic",L_CSDxdot_ic_) ;
#line 723 "FSI_CSDdataloop.loci"
       input("CSDEulerXnum") ;
#line 723 "FSI_CSDdataloop.loci"
       output("CSDxdot_ic") ;
#line 723 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 723 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 723 "FSI_CSDdataloop.loci"
    }
#line 723 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	
	int CSDNumBC = 2 ;
	int CSDdofNode = 2 ;
	int CSDdofFree = *L_CSDEulerXnum_* CSDdofNode - CSDNumBC ;
	
	(*L_CSDxdot_ic_).resize(CSDdofFree) ; std::fill( (*L_CSDxdot_ic_).begin(),(*L_CSDxdot_ic_).end(), 0. );
}    void compute(const Loci::sequence &seq) { 
#line 732 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 732 "FSI_CSDdataloop.loci"
    }
#line 732 "FSI_CSDdataloop.loci"
} ;
#line 732 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop031_1279551148m545> register_file_FSI_CSDdataloop031_1279551148m545 ;
#line 732 "FSI_CSDdataloop.loci"
}
#line 732 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop032_1279551148m545 : public Loci::blackbox_rule {
#line 732 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxdot_ic_ ; 
#line 732 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxdot_n_EQ_0__ ; 
#line 732 "FSI_CSDdataloop.loci"
public:
#line 732 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop032_1279551148m545() {
#line 732 "FSI_CSDdataloop.loci"
       name_store("CSDxdot_ic",L_CSDxdot_ic_) ;
#line 732 "FSI_CSDdataloop.loci"
       name_store("CSDxdot{n=0}",L_CSDxdot_n_EQ_0__) ;
#line 732 "FSI_CSDdataloop.loci"
       input("CSDxdot_ic") ;
#line 732 "FSI_CSDdataloop.loci"
       output("CSDxdot{n=0}") ;
#line 732 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 732 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 732 "FSI_CSDdataloop.loci"
    }
#line 732 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDxdot_n_EQ_0__).resize((*L_CSDxdot_ic_).size()) ; std::fill( (*L_CSDxdot_n_EQ_0__).begin(),(*L_CSDxdot_n_EQ_0__).end(), 0. );
	(*L_CSDxdot_n_EQ_0__) = (*L_CSDxdot_ic_) ;
}    void compute(const Loci::sequence &seq) { 
#line 738 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 738 "FSI_CSDdataloop.loci"
    }
#line 738 "FSI_CSDdataloop.loci"
} ;
#line 738 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop032_1279551148m545> register_file_FSI_CSDdataloop032_1279551148m545 ;
#line 738 "FSI_CSDdataloop.loci"
}
#line 738 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop033_1279551148m546 : public Loci::blackbox_rule {
#line 738 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxdot_n__ ; 
#line 738 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxdotStar_n__ ; 
#line 738 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxdot_n_P_1__ ; 
#line 738 "FSI_CSDdataloop.loci"
public:
#line 738 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop033_1279551148m546() {
#line 738 "FSI_CSDdataloop.loci"
       name_store("CSDxdot{n}",L_CSDxdot_n__) ;
#line 738 "FSI_CSDdataloop.loci"
       name_store("CSDxdotStar{n}",L_CSDxdotStar_n__) ;
#line 738 "FSI_CSDdataloop.loci"
       name_store("CSDxdot{n+1}",L_CSDxdot_n_P_1__) ;
#line 738 "FSI_CSDdataloop.loci"
       input("CSDxdot{n},CSDxdotStar{n}") ;
#line 738 "FSI_CSDdataloop.loci"
       output("CSDxdot{n+1}") ;
#line 738 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 738 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 738 "FSI_CSDdataloop.loci"
    }
#line 738 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDxdot_n_P_1__).resize((*L_CSDxdot_n__).size()) ;std::fill( (*L_CSDxdot_n_P_1__).begin(),(*L_CSDxdot_n_P_1__).end(), 0. );
	(*L_CSDxdot_n_P_1__) = (*L_CSDxdotStar_n__) ;	
}    void compute(const Loci::sequence &seq) { 
#line 744 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 744 "FSI_CSDdataloop.loci"
    }
#line 744 "FSI_CSDdataloop.loci"
} ;
#line 744 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop033_1279551148m546> register_file_FSI_CSDdataloop033_1279551148m546 ;
#line 744 "FSI_CSDdataloop.loci"
}
#line 744 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop034_1279551148m546 : public Loci::blackbox_rule {
#line 744 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerXnum_ ; 
#line 744 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxddot_ic_ ; 
#line 744 "FSI_CSDdataloop.loci"
public:
#line 744 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop034_1279551148m546() {
#line 744 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXnum",L_CSDEulerXnum_) ;
#line 744 "FSI_CSDdataloop.loci"
       name_store("CSDxddot_ic",L_CSDxddot_ic_) ;
#line 744 "FSI_CSDdataloop.loci"
       input("CSDEulerXnum") ;
#line 744 "FSI_CSDdataloop.loci"
       output("CSDxddot_ic") ;
#line 744 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 744 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 744 "FSI_CSDdataloop.loci"
    }
#line 744 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	
	int CSDNumBC = 2 ;
	int CSDdofNode = 2 ;
	int CSDdofFree = *L_CSDEulerXnum_* CSDdofNode - CSDNumBC ;
	
	(*L_CSDxddot_ic_).resize(CSDdofFree) ; std::fill( (*L_CSDxddot_ic_).begin(),(*L_CSDxddot_ic_).end(), 0. );
}    void compute(const Loci::sequence &seq) { 
#line 753 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 753 "FSI_CSDdataloop.loci"
    }
#line 753 "FSI_CSDdataloop.loci"
} ;
#line 753 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop034_1279551148m546> register_file_FSI_CSDdataloop034_1279551148m546 ;
#line 753 "FSI_CSDdataloop.loci"
}
#line 753 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop035_1279551148m547 : public Loci::blackbox_rule {
#line 753 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxddot_ic_ ; 
#line 753 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxddot_n_EQ_0__ ; 
#line 753 "FSI_CSDdataloop.loci"
public:
#line 753 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop035_1279551148m547() {
#line 753 "FSI_CSDdataloop.loci"
       name_store("CSDxddot_ic",L_CSDxddot_ic_) ;
#line 753 "FSI_CSDdataloop.loci"
       name_store("CSDxddot{n=0}",L_CSDxddot_n_EQ_0__) ;
#line 753 "FSI_CSDdataloop.loci"
       input("CSDxddot_ic") ;
#line 753 "FSI_CSDdataloop.loci"
       output("CSDxddot{n=0}") ;
#line 753 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 753 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 753 "FSI_CSDdataloop.loci"
    }
#line 753 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDxddot_n_EQ_0__).resize((*L_CSDxddot_ic_).size()) ; std::fill( (*L_CSDxddot_n_EQ_0__).begin(),(*L_CSDxddot_n_EQ_0__).end(), 0. );
	(*L_CSDxddot_n_EQ_0__) = (*L_CSDxddot_ic_) ;
}    void compute(const Loci::sequence &seq) { 
#line 759 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 759 "FSI_CSDdataloop.loci"
    }
#line 759 "FSI_CSDdataloop.loci"
} ;
#line 759 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop035_1279551148m547> register_file_FSI_CSDdataloop035_1279551148m547 ;
#line 759 "FSI_CSDdataloop.loci"
}
#line 759 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop036_1279551148m547 : public Loci::blackbox_rule {
#line 759 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxddot_n__ ; 
#line 759 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::vector<real> >  L_CSDxddotStar_n__ ; 
#line 759 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::vector<real> >  L_CSDxddot_n_P_1__ ; 
#line 759 "FSI_CSDdataloop.loci"
public:
#line 759 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop036_1279551148m547() {
#line 759 "FSI_CSDdataloop.loci"
       name_store("CSDxddot{n}",L_CSDxddot_n__) ;
#line 759 "FSI_CSDdataloop.loci"
       name_store("CSDxddotStar{n}",L_CSDxddotStar_n__) ;
#line 759 "FSI_CSDdataloop.loci"
       name_store("CSDxddot{n+1}",L_CSDxddot_n_P_1__) ;
#line 759 "FSI_CSDdataloop.loci"
       input("CSDxddot{n},CSDxddotStar{n}") ;
#line 759 "FSI_CSDdataloop.loci"
       output("CSDxddot{n+1}") ;
#line 759 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIEULERBEAM") ;
#line 759 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 759 "FSI_CSDdataloop.loci"
    }
#line 759 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) {  
	(*L_CSDxddot_n_P_1__).resize((*L_CSDxddotStar_n__).size()) ;std::fill( (*L_CSDxddot_n_P_1__).begin(),(*L_CSDxddot_n_P_1__).end(), 0. );
	(*L_CSDxddot_n_P_1__) = (*L_CSDxddotStar_n__) ;	
}    void compute(const Loci::sequence &seq) { 
#line 764 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 764 "FSI_CSDdataloop.loci"
    }
#line 764 "FSI_CSDdataloop.loci"
} ;
#line 764 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop036_1279551148m547> register_file_FSI_CSDdataloop036_1279551148m547 ;
#line 764 "FSI_CSDdataloop.loci"
}
#line 764 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop037_1279551148m547 : public Loci::blackbox_rule {
#line 764 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDEulerXstart_ ; 
#line 764 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDEulerXend_ ; 
#line 764 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerXnum_ ; 
#line 764 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 764 "FSI_CSDdataloop.loci"
public:
#line 764 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop037_1279551148m547() {
#line 764 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 764 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXstart",L_CSDEulerXstart_) ;
#line 764 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXend",L_CSDEulerXend_) ;
#line 764 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXnum",L_CSDEulerXnum_) ;
#line 764 "FSI_CSDdataloop.loci"
       input("CSDEulerXstart,CSDEulerXend,CSDEulerXnum") ;
#line 764 "FSI_CSDdataloop.loci"
       output("CSDnodes_ic") ;
#line 764 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIPRESCRIBED") ;
#line 764 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 764 "FSI_CSDdataloop.loci"
    }
#line 764 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) { 
	int nElements = *L_CSDEulerXnum_- 1;
	real dx = (*L_CSDEulerXend_- *L_CSDEulerXstart_) / nElements ;
	
	(*L_CSDnodes_ic_).resize(*L_CSDEulerXnum_,3) ; std::fill( (*L_CSDnodes_ic_).data().begin(), (*L_CSDnodes_ic_).data().end(), 0. );	    	
	
	(*L_CSDnodes_ic_)(0,0) = *L_CSDEulerXstart_;
	for(int i=0; i<*L_CSDEulerXnum_; ++i) {
		(*L_CSDnodes_ic_)(i,0) = *L_CSDEulerXstart_+ i * dx;
	}
}    void compute(const Loci::sequence &seq) { 
#line 776 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 776 "FSI_CSDdataloop.loci"
    }
#line 776 "FSI_CSDdataloop.loci"
} ;
#line 776 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop037_1279551148m547> register_file_FSI_CSDdataloop037_1279551148m547 ;
#line 776 "FSI_CSDdataloop.loci"
}
#line 776 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop038_1279551148m548 : public Loci::blackbox_rule {
#line 780 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDfrequency_ ; 
#line 780 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_CSDflappingAmplitudeX_ ; 
#line 780 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDstartingTimeStep_ ; 
#line 780 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDnodes_ic_ ; 
#line 780 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_CSDEulerXnum_ ; 
#line 780 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_stime_n__ ; 
#line 780 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_ncycle_n__ ; 
#line 780 "FSI_CSDdataloop.loci"
    Loci::const_param<real>  L_timeStep_ ; 
#line 780 "FSI_CSDdataloop.loci"
    Loci::const_blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDForce_n__ ; 
#line 780 "FSI_CSDdataloop.loci"
    Loci::blackbox<ublas::matrix<real,ublas::column_major> >  L_CSDdisplacementsStar_n__ ; 
#line 780 "FSI_CSDdataloop.loci"
public:
#line 780 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop038_1279551148m548() {
#line 780 "FSI_CSDdataloop.loci"
       name_store("CSDfrequency",L_CSDfrequency_) ;
#line 780 "FSI_CSDdataloop.loci"
       name_store("CSDflappingAmplitudeX",L_CSDflappingAmplitudeX_) ;
#line 780 "FSI_CSDdataloop.loci"
       name_store("CSDstartingTimeStep",L_CSDstartingTimeStep_) ;
#line 780 "FSI_CSDdataloop.loci"
       name_store("CSDnodes_ic",L_CSDnodes_ic_) ;
#line 780 "FSI_CSDdataloop.loci"
       name_store("CSDEulerXnum",L_CSDEulerXnum_) ;
#line 780 "FSI_CSDdataloop.loci"
       name_store("CSDdisplacementsStar{n}",L_CSDdisplacementsStar_n__) ;
#line 780 "FSI_CSDdataloop.loci"
       name_store("stime{n}",L_stime_n__) ;
#line 780 "FSI_CSDdataloop.loci"
       name_store("ncycle{n}",L_ncycle_n__) ;
#line 780 "FSI_CSDdataloop.loci"
       name_store("timeStep",L_timeStep_) ;
#line 780 "FSI_CSDdataloop.loci"
       name_store("CSDForce{n}",L_CSDForce_n__) ;
#line 780 "FSI_CSDdataloop.loci"
       input("							CSDEulerXnum,CSDnodes_ic,CSDForce{n},							stime{n},ncycle{n},timeStep,CSDstartingTimeStep,								CSDfrequency,CSDflappingAmplitudeX") ;
#line 780 "FSI_CSDdataloop.loci"
       output("CSDdisplacementsStar{n}") ;
#line 780 "FSI_CSDdataloop.loci"
       constraint("FSICoupling,FSIPRESCRIBED") ;
#line 780 "FSI_CSDdataloop.loci"
       disable_threading() ;
#line 780 "FSI_CSDdataloop.loci"
    }
#line 780 "FSI_CSDdataloop.loci"
    void prelude(const Loci::sequence &seq) { 
	// call NLAMS
	//communicateWithNLAMS(*CSDdisplacementsStar, ....) ;
	//if (*$CFDIterationFinished) {
	
	const int rank = Loci::MPI_rank ;

	if (rank==0) cout << "[I] Communicating with the prescribed motion rule " << endl ;
	// Communicate with NLAMS
	int CSDNumNodes = (*L_CSDEulerXnum_) ;
	
	// Initialize
	(*L_CSDdisplacementsStar_n__).resize(CSDNumNodes,3) ; std::fill( (*L_CSDdisplacementsStar_n__).data().begin(),(*L_CSDdisplacementsStar_n__).data().end(), 0. );	
	
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
		 	  (*L_CSDdisplacementsStar_n__)(i,0) = (*L_CSDnodes_ic_)(i,0) * ( cos(alpha*PI/180) - 1.) ;
		 	  (*L_CSDdisplacementsStar_n__)(i,1) = (*L_CSDnodes_ic_)(i,0) * ( sin(alpha*PI/180) ) ;
		 	  (*L_CSDdisplacementsStar_n__)(i,2) = 0.0 ;
		  }
	   	if (Loci::MPI_rank==0) cout << "CSDdisplacement:" << (*L_stime_n__) << ", alpha = " << alpha << endl ;
	  	//if (Loci::MPI_rank==0) cout << "CSDnodes_ic: " << (*$CSDnodes_ic) << endl ;
	  	//if (Loci::MPI_rank==0) cout << "CSDdisplacementsStar: " << (*$CSDdisplacementsStar{n}) << endl ;	
		  
		  if (rank==0) cout << "[I] Communication with the prescribed motion rule done. " << endl ;
		}   
	
}    void compute(const Loci::sequence &seq) { 
#line 823 "FSI_CSDdataloop.loci"
      prelude(seq) ;
#line 823 "FSI_CSDdataloop.loci"
    }
#line 823 "FSI_CSDdataloop.loci"
} ;
#line 823 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop038_1279551148m548> register_file_FSI_CSDdataloop038_1279551148m548 ;
#line 823 "FSI_CSDdataloop.loci"
}
#line 823 "FSI_CSDdataloop.loci"
namespace {class file_FSI_CSDdataloop039_1279551148m549 : public Loci::singleton_rule {
#line 823 "FSI_CSDdataloop.loci"
    Loci::param<int>  L_itfsi_ic_ ; 
#line 823 "FSI_CSDdataloop.loci"
public:
#line 823 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop039_1279551148m549() {
#line 823 "FSI_CSDdataloop.loci"
       name_store("itfsi_ic",L_itfsi_ic_) ;
#line 823 "FSI_CSDdataloop.loci"
       output("itfsi_ic") ;
#line 823 "FSI_CSDdataloop.loci"
       constraint("UNIVERSE") ;
#line 823 "FSI_CSDdataloop.loci"
    }
#line 823 "FSI_CSDdataloop.loci"
    void compute(const Loci::sequence &seq) { 
	(*L_itfsi_ic_)= 0 ;
}} ;
#line 825 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop039_1279551148m549> register_file_FSI_CSDdataloop039_1279551148m549 ;
#line 825 "FSI_CSDdataloop.loci"
}
#line 825 "FSI_CSDdataloop.loci"


namespace {class file_FSI_CSDdataloop040_1279551148m549 : public Loci::singleton_rule {
#line 827 "FSI_CSDdataloop.loci"
    Loci::const_param<int>  L_itfsi_ic_ ; 
#line 827 "FSI_CSDdataloop.loci"
    Loci::param<int>  L_itfsi_n_EQ_0__ ; 
#line 827 "FSI_CSDdataloop.loci"
public:
#line 827 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop040_1279551148m549() {
#line 827 "FSI_CSDdataloop.loci"
       name_store("itfsi_ic",L_itfsi_ic_) ;
#line 827 "FSI_CSDdataloop.loci"
       name_store("itfsi{n=0}",L_itfsi_n_EQ_0__) ;
#line 827 "FSI_CSDdataloop.loci"
       input("itfsi_ic") ;
#line 827 "FSI_CSDdataloop.loci"
       output("itfsi{n=0}") ;
#line 827 "FSI_CSDdataloop.loci"
    }
#line 827 "FSI_CSDdataloop.loci"
    void compute(const Loci::sequence &seq) { 
	(*L_itfsi_n_EQ_0__)=(*L_itfsi_ic_);
}} ;
#line 829 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop040_1279551148m549> register_file_FSI_CSDdataloop040_1279551148m549 ;
#line 829 "FSI_CSDdataloop.loci"
}
#line 829 "FSI_CSDdataloop.loci"


namespace {class file_FSI_CSDdataloop041_1279551148m549 : public Loci::singleton_rule {
#line 831 "FSI_CSDdataloop.loci"
    Loci::param<int>  L_itfsi_n_P_1__ ; 
#line 831 "FSI_CSDdataloop.loci"
public:
#line 831 "FSI_CSDdataloop.loci"
    file_FSI_CSDdataloop041_1279551148m549() {
#line 831 "FSI_CSDdataloop.loci"
       name_store("itfsi{n+1}",L_itfsi_n_P_1__) ;
#line 831 "FSI_CSDdataloop.loci"
       output("itfsi{n+1}") ;
#line 831 "FSI_CSDdataloop.loci"
       conditional("iterationFinished{n,it-1}") ;
#line 831 "FSI_CSDdataloop.loci"
    }
#line 831 "FSI_CSDdataloop.loci"
    void compute(const Loci::sequence &seq) { 
	
	//if (Loci::MPI_rank==0) cout << "inside itfsi{n+1} in" << endl ;
	(*L_itfsi_n_P_1__) = 0 ;
	//if (Loci::MPI_rank==0) cout << "inside itfsi{n+1} out" << endl ;
	
}} ;
#line 837 "FSI_CSDdataloop.loci"
Loci::register_rule<file_FSI_CSDdataloop041_1279551148m549> register_file_FSI_CSDdataloop041_1279551148m549 ;
#line 837 "FSI_CSDdataloop.loci"
}
#line 837 "FSI_CSDdataloop.loci"


}
