#line 1 "FSI_iterationFinishedPriority.loci"
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
#include <vector>
using std::vector ;

// Loci includes.
#include <Loci.h>
using Loci::Area ;


// StreamUns includes.
#include "const.h"
#include "residual.h"
#include "sciTypes.h"
#include "varsFileInputs.h"
                                      
// boost::ublas includes
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas ;


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

#line 35 "FSI_iterationFinishedPriority.loci"


namespace streamUns {

// Rule to set the CFDIterationFinished constraint
// $type CFDIterationFinished param<bool> 
// $type CFDIterationFinishedConstraint Constraint


// Class to determine if iteration is finished.
  class CheckIterationFinishedCFD : public singleton_rule {
    private:
      const_param<int> numSpecies ;
      const_param<int> nCycle,it ;
      const_param<int> maxIterationsPerFSI ;
      const_param<real> convergenceTolerance ;
      const_param<VectorResidual> vResidualData ;
      const_param<ScalarResidual> pPrimeResidualData ;
      const_param<ScalarResidual> hResidualData ;
      const_param<ScalarResidual> kResidualData ;
      const_param<ScalarResidual> omegaResidualData ;
      const_param<vector<real> > totalSpeciesResidual ;
      const_param<int> CFDMaxTotalInnerIterations ; //add
      param<bool> CFDIterationFinished ;
    public:

      // Define input and output.
      CheckIterationFinishedCFD() {
        name_store("numSpecies{n,it}",numSpecies) ;
        name_store("ncycle{n}",nCycle) ;
        name_store("$it{n,it}",it) ;
        name_store("maxIterationsPerFSI{n,it}",maxIterationsPerFSI) ;
        name_store("convergenceTolerance{n,it}",convergenceTolerance) ;
        name_store("vResidualData{n,it}",vResidualData) ;
        name_store("pPrimeResidualData{n,it}",pPrimeResidualData) ;
        name_store("hResidualData{n,it}",hResidualData) ;
        name_store("kResidualData{n,it}",kResidualData) ;
        name_store("omegaResidualData{n,it}",omegaResidualData) ;
        name_store("totalSpeciesResidual{n,it}",totalSpeciesResidual) ;
        name_store("CFDIterationFinished{n,it}",CFDIterationFinished) ;
        name_store("CFDMaxTotalInnerIterations", CFDMaxTotalInnerIterations) ;
        input("numSpecies{n,it},ncycle{n},$it{n,it}") ;
        input("maxIterationsPerFSI{n,it},convergenceTolerance{n,it}") ;
        input("vResidualData{n,it},pPrimeResidualData{n,it}") ;
        input("hResidualData{n,it},kResidualData{n,it}") ;
        input("omegaResidualData{n,it},totalSpeciesResidual{n,it}") ;
        input("CFDMaxTotalInnerIterations") ; // add
        output("CFDIterationFinished{n,it}") ;
        constraint("vResidualData{n,it},FSICoupling") ;
      }

      // Check if iteration is finished.
      void compute(const sequence &seq) {

        // Set the output format.
        if(Loci::MPI_rank==0){
          cout.setf(ios::scientific,ios::floatfield) ; cout.precision(6) ;
        }

	
        // Compute the maximum species residual.
//	if (Loci::MPI_rank==0) cout << "CFDIterationFinished..0 it, maxIterationsPerFSI" << (*it) << ", " << (*maxIterationsPerFSI) << endl ;
        *CFDIterationFinished = ((*it)%(*maxIterationsPerFSI-1)) == 0  ;
//	if (Loci::MPI_rank==0) cout << "CFDIterationFinished..1 " << endl ;
	*CFDIterationFinished = (*CFDIterationFinished && (*it)>2) ; // To avoid some residual increase at the first few iteration
//	if (Loci::MPI_rank==0) cout << "CFDIterationFinished..2 " << endl ;
	*CFDIterationFinished = (*CFDIterationFinished || *it==*CFDMaxTotalInnerIterations) ;
//	if (Loci::MPI_rank==0) cout << "CFDIterationFinished..3 " << endl ;
          
	if (*CFDIterationFinished) {
	      if (Loci::MPI_rank==0) cout << "CFDIterationFinished.. " << endl ;
	      } else {
	      if (Loci::MPI_rank==0) cout << "CFDIteration not yet Finished.. " << endl ;
	      }
	}
  } ;

  register_rule<CheckIterationFinishedCFD> registerCheckIterationFinishedCFD ;

 // Build rule.
  class CFDIterationFinishedBuild : public singleton_rule {
    private:
      param<bool> CFDIterationFinished ;
    public:
                                                                                
      // Define input and output.
      CFDIterationFinishedBuild() {
        name_store("CFDIterationFinished{n,it=-1}",CFDIterationFinished) ;
        output("CFDIterationFinished{n,it=-1}") ;
        constraint("FSICoupling{n}") ;
      }
                                                                                
      // Set to false.
      void compute(const sequence &seq) { *CFDIterationFinished=false ; }
  } ;
                                                                                
  register_rule<CFDIterationFinishedBuild> registerCFDIterationFinishedBuild ;


// Class to determine if iteration is finished.
//  class CheckIterationFinishedFSI : public singleton_rule {
//    private:
//      const_param<bool> CFDIterationFinished ;
//      const_param<bool> FSIIterationFinished ;
//      const_param<int> CFDMaxTotalInnerIterations ;
//      const_param<int> nCycle,it ;
//      param<bool> iterationFinished ;
//    public:
//
//      // Define input and output.
//      CheckIterationFinishedFSI() {
//        name_store("FSIIterationFinished{n,it-1}",FSIIterationFinished) ;
//        name_store("CFDIterationFinished{n,it}",CFDIterationFinished) ;
//        name_store("priority::iterationFinished{n,it}",iterationFinished) ;
//        name_store("CFDMaxTotalInnerIterations", CFDMaxTotalInnerIterations) ;
//        name_store("ncycle{n}",nCycle) ;
//        name_store("$it{n,it}",it) ;
//        input("CFDMaxTotalInnerIterations") ;
//        input("ncycle{n},$it{n,it}") ;
//        input("CFDIterationFinished{n,it}") ;
//				input("FSIIterationFinished{n,it-1}") ;
//				output("priority::iterationFinished{n,it}") ;
//        constraint("vResidualData{n,it},FSICoupling") ;
////        conditional("FSIIterationFinished{n,itfsi}") ;
//				disable_threading() ;
//      }
//
//      // Check if iteration is finished.
//      void compute(const sequence &seq) {
//
//        // See if we are converged.
//             	
//        *iterationFinished=(((*FSIIterationFinished) || (*it >= *CFDMaxTotalInnerIterations)) && *CFDIterationFinished ) ;
//        if (*FSIIterationFinished) cout << "FSIIterationFinished" << endl ;
//				if (*CFDIterationFinished) cout << "CFDIterationFinished" << endl ;   
//				if (*iterationFinished) cout << "iterationFinished" << endl ;   					
////      *iterationFinished=(*CFDIterationFinished) ;
////         *iterationFinished= (*FSIIterationFinished) ; // modified
//
//         if (*iterationFinished) {
//          	if (Loci::MPI_rank==0) cout << "FSI IterationFinished.. " << endl ;
//         	} else {
//       		if (Loci::MPI_rank==0) cout << "FSI Iteration not yet Finished.. " << endl ;
//         	}          
//  //    *iterationFinished=false ; // test --> iteration never finishes?
//      }
//  } ;
//
//  register_rule<CheckIterationFinishedFSI> registerCheckIterationFinishedFSI ;

// Class to determine if iteration is finished.
//  class CheckIterationFinishedCFD : public singleton_rule {
//    private:
//      const_param<int> numSpecies ;
//      const_param<int> nCycle,it ;
//      const_param<int> maxIterationsPerTimeStep ;
//      const_param<real> convergenceTolerance ;
//      const_param<VectorResidual> vResidualData ;
//      const_param<ScalarResidual> pPrimeResidualData ;
//      const_param<ScalarResidual> hResidualData ;
//      const_param<ScalarResidual> kResidualData ;
//      const_param<ScalarResidual> omegaResidualData ;
//      const_param<vector<real> > totalSpeciesResidual ;
//      const_param<int> CFDMaxTotalInnerIterations ; //add
//      param<bool> CFDIterationFinished ;
//    public:
//
//      // Define input and output.
//      CheckIterationFinishedCFD() {
//        name_store("numSpecies{n,it}",numSpecies) ;
//        name_store("ncycle{n}",nCycle) ;
//        name_store("$it{n,it}",it) ;
//        name_store("maxIterationsPerTimeStep{n,it}",maxIterationsPerTimeStep) ;
//        name_store("convergenceTolerance{n,it}",convergenceTolerance) ;
//        name_store("vResidualData{n,it}",vResidualData) ;
//        name_store("pPrimeResidualData{n,it}",pPrimeResidualData) ;
//        name_store("hResidualData{n,it}",hResidualData) ;
//        name_store("kResidualData{n,it}",kResidualData) ;
//        name_store("omegaResidualData{n,it}",omegaResidualData) ;
//        name_store("totalSpeciesResidual{n,it}",totalSpeciesResidual) ;
//        name_store("CFDIterationFinished{n,it}",CFDIterationFinished) ;
//        name_store("CFDMaxTotalInnerIterations", CFDMaxTotalInnerIterations) ;
//        input("numSpecies{n,it},ncycle{n},$it{n,it}") ;
//        input("maxIterationsPerTimeStep{n,it},convergenceTolerance{n,it}") ;
//        input("vResidualData{n,it},pPrimeResidualData{n,it}") ;
//        input("hResidualData{n,it},kResidualData{n,it}") ;
//        input("omegaResidualData{n,it},totalSpeciesResidual{n,it}") ;
//        input("CFDMaxTotalInnerIterations") ; // add
//        output("CFDIterationFinished{n,it}") ;
//        constraint("vResidualData{n,it},FSICoupling") ;
//      }
//
//      // Check if iteration is finished.
//      void compute(const sequence &seq) {
//
//        // Set the output format.
//        if(Loci::MPI_rank==0){
//          cout.setf(ios::scientific,ios::floatfield) ; cout.precision(6) ;
//        }
//
//        // Compute the maximum species residual.
//        real maxTotalSpeciesResidual=0.0 ;
//        if(*numSpecies>1){
//          for(int i=0;i<*numSpecies;++i)
//            if((*totalSpeciesResidual)[i]>maxTotalSpeciesResidual)
//              maxTotalSpeciesResidual=(*totalSpeciesResidual)[i] ;
//        }
//
//        // Write out the residuals.
//        /*if(Loci::MPI_rank==0){
//          cout <<"R: " << *nCycle << " " << *it << " " << vResidualData->
//            totalResidual << pPrimeResidualData->totalResidual ;
//          if(hResidualData->totalResidual!=0.0)
//            cout << " " << hResidualData->totalResidual ;
//          if(kResidualData->totalResidual!=0.0)
//            cout << " " << kResidualData->totalResidual << " "
//              << omegaResidualData->totalResidual ;
//          if(maxTotalSpeciesResidual!=0.0)
//            cout << " " << maxTotalSpeciesResidual ;
//          cout << endl ;
//        } */ // Already output inside momentum.cc when checking for iterationFinished
//
//        // See if we are converged.
////         *CFDIterationFinished=(*it==*maxIterationsPerTimeStep-1 ||
////          (vResidualData->totalResidual.x<*convergenceTolerance &&
////          vResidualData->totalResidual.y<*convergenceTolerance &&
////          vResidualData->totalResidual.z<*convergenceTolerance &&
////          pPrimeResidualData->totalResidual<*convergenceTolerance &&
////          hResidualData->totalResidual<*convergenceTolerance &&
////          kResidualData->totalResidual<*convergenceTolerance &&
////          omegaResidualData->totalResidual<*convergenceTolerance &&
////          maxTotalSpeciesResidual<*convergenceTolerance)) ;
//
//        *CFDIterationFinished=(((((*it)%(*maxIterationsPerTimeStep-1))==0 && (*it)>1) ||
//          (vResidualData->totalResidual.x<*convergenceTolerance &&
//          vResidualData->totalResidual.y<*convergenceTolerance &&
//          vResidualData->totalResidual.z<*convergenceTolerance &&
//          pPrimeResidualData->totalResidual<*convergenceTolerance &&
//          hResidualData->totalResidual<*convergenceTolerance &&
//          kResidualData->totalResidual<*convergenceTolerance &&
//          omegaResidualData->totalResidual<*convergenceTolerance &&
//          maxTotalSpeciesResidual<*convergenceTolerance))|| *it==*CFDMaxTotalInnerIterations);
//          
//          if (*CFDIterationFinished) {
//          	if (Loci::MPI_rank==0) cout << "CFDIterationFinished.. " << endl ;
//         	} else {
//      		if (Loci::MPI_rank==0) cout << "CFDIteration not yet Finished.. " << endl ;
//         	}
//      }
//  } ;
//
//  register_rule<CheckIterationFinishedCFD> registerCheckIterationFinishedCFD ;
//
// // Build rule.
//  class CFDIterationFinishedBuild : public singleton_rule {
//    private:
//      param<bool> CFDIterationFinished ;
//    public:
//                                                                                
//      // Define input and output.
//      CFDIterationFinishedBuild() {
//        name_store("CFDIterationFinished{n,it=-1}",CFDIterationFinished) ;
//        output("CFDIterationFinished{n,it=-1}") ;
//        constraint("FSICoupling{n}") ;
//      }
//                                                                                
//      // Set to false.
//      void compute(const sequence &seq) { *CFDIterationFinished=false ; }
//  } ;
//                                                                                
//  register_rule<CFDIterationFinishedBuild> registerCFDIterationFinishedBuild ;

// // Build rule.
//  class FSIIterationFinishedBuild : public singleton_rule {
//    private:
//      param<bool> FSIIterationFinished ;
//    public:
//                                                                                
//      // Define input and output.
//      FSIIterationFinishedBuild() {
//        name_store("FSIIterationFinished{n,it=-1}",FSIIterationFinished) ;
//        output("FSIIterationFinished{n,it=-1}") ;
//        constraint("FSICoupling{n}") ;
//      }
//                                                                                
//      // Set to false.
//      void compute(const sequence &seq) { *FSIIterationFinished=false ; }
//  } ;
//                                                                                
//  register_rule<FSIIterationFinishedBuild> registerFSIIterationFinishedBuild ;

//$type FSIIterationFinished param<bool> ;
//$type FSIIterationFinishedStar param<bool> ;
//
//$rule singleton(FSIIterationFinished{n,it=-1}),constraint(FSICoupling{n}) {
////	if (Loci::MPI_rank==0) cout << "inside itfsi{n,it=0} in" << endl ;
//	$FSIIterationFinished{n,it=-1}=false ;
////	if (Loci::MPI_rank==0) cout << "inside itfsi{n,it=0} out" << endl ; 
//}
//
//$rule singleton(FSIIterationFinished{n,it=0}),constraint(FSICoupling{n}) {
////	if (Loci::MPI_rank==0) cout << "inside itfsi{n,it=0} in" << endl ;
//	$FSIIterationFinished{n,it=0}=false ;
////	if (Loci::MPI_rank==0) cout << "inside itfsi{n,it=0} out" << endl ; 
//}
//
//$rule singleton(FSIIterationFinished{n,it+1}<-FSIIterationFinishedStar{n,it},FSIIterationFinished{n,it}),constraint(FSICoupling{n,it})  { 	
//	$FSIIterationFinished{n,it+1} = ($FSIIterationFinished{n,it} || $FSIIterationFinishedStar{n,it}) ; // initially FSIIterationFinished false, but if FSIIterationFinishedStar becomes true any time, FSIIterationFinished stays true.	
//}
//
//// Rule to check if the FSI inner iteration is finished by convergence of the nodal displacements
//// Class to determine if iteration is finished.
//  class CheckIterationConvergedFSI : public singleton_rule {
//    private:
//      const_param<int> itfsi ;
//      const_param<int> ncycle ;
//      const_param<int> FSIIterationMinimum ;
//      const_blackbox<ublas::matrix<real,ublas::column_major> > CSDnodesCurrent ;
//      const_blackbox<ublas::matrix<real,ublas::column_major> > CSDnodesRef ;	
//      const_blackbox<ublas::matrix<real,ublas::column_major> > CSDdisplacementsStar ;
//      const_param<real> FSIIterationTolerance ;
//      const_param<bool> CFDIterationFinished ;
//      param<bool> FSIIterationFinishedStar ;
//    public:
//
//      // Define input and output.
//      CheckIterationConvergedFSI() {
//        name_store("CSDnodes_it{n,it}",CSDnodesCurrent) ;
//        name_store("CSDnodes_ic",CSDnodesRef) ;
//        name_store("CSDdisplacementsStar{n,it}",CSDdisplacementsStar) ;
//        name_store("FSIIterationFinishedStar{n,it}",FSIIterationFinishedStar) ;
//        name_store("FSIIterationTolerance{n,it}",FSIIterationTolerance) ;
//        name_store("itfsi{n,it}",itfsi) ;        
//        name_store("CFDIterationFinished{n,it-1}",CFDIterationFinished) ;
//        name_store("ncycle{n,it}",ncycle) ;
//        name_store("FSIIterationMinimum", FSIIterationMinimum) ;
//        input("FSIIterationMinimum") ;
//        input("ncycle{n,it}") ;
//        input("itfsi{n,it}") ;
//        input("CSDnodes_it{n,it}") ;
//        input("CSDnodes_ic") ;
//        input("CSDdisplacementsStar{n,it}") ;
//        input("FSIIterationTolerance{n,it}") ;
//        input("CFDIterationFinished{n,it-1}") ;
//        output("FSIIterationFinishedStar{n,it}") ;
//        constraint("FSICoupling{n,it}") ;
//  //      conditional("CFDIterationFinished{n,it-1}") ;
//        disable_threading() ;
//      }
//
//      // Check if iteration is finished.
//      void compute(const sequence &seq) {
//
//			if (Loci::MPI_rank==0) cout << "FSI: here_start" << endl ;
//			if (*CFDIterationFinished) {
////					if (Loci::MPI_rank==0) cout << "FSI: here_start1" << endl ;
//					*FSIIterationFinishedStar = false ;
//	        // See if we are converged only if $it > 1
//	        
//	        if (Loci::MPI_rank==0) cout << "itfsi = " << *itfsi << endl ;
//	        if (*itfsi >= *FSIIterationMinimum) {
//	        	
////	        	if (Loci::MPI_rank==0) cout << "FSI: here_start2" << endl ;
//	        	ublas::matrix<real,ublas::column_major> CSDdisplacementsDifference ;
//	        		
//	// cout << "FSI: here0" << endl ;	
//		        (CSDdisplacementsDifference) = (*CSDdisplacementsStar) - ((*CSDnodesCurrent) - (*CSDnodesRef)) ;	
//		        
//		//cout << "FSI: here" << endl ;
//		        if (Loci::MPI_rank==0) cout << "[I] Computing the relative displacement L2 norms" << endl ;
//		        real dCurrent=0., dOld=0., sumCurrent=0., LNORM=0.;          
//		        for(int i=0;i<(*CSDdisplacementsStar).size1();++i) {	        	
//		        	dCurrent=sqrt(pow((CSDdisplacementsDifference)(i,0), 2.)+pow((CSDdisplacementsDifference)(i,1), 2.)+pow((CSDdisplacementsDifference)(i,2), 2.)); // length of the difference
//		        	//dOld=pow((CSDdisplacementsOld)(i,0), 2.)+pow((CSDdisplacementsOld)(i,1), 2.)+pow((CSDdisplacementsOld)(i,2), 2.);
//		        	LNORM += dCurrent ;
//		        	//sumCurrent += dCurrent ;
//		        }
//		        
////		        cout << "FSI: here2 sumCurrent = " << sumCurrent << endl ;
//		        
//		        real Ncsd = (real) (*CSDdisplacementsStar).size1() ;
//		        
//		        real NORM = LNORM / Ncsd ;
//		        
////		        if (sumCurrent > 0.0) {
////		        	L2RELNORM = sqrt(L2RELNORM / sumCurrent) ;
////		        } else {
////		        	L2RELNORM = 0.0 ;
////		        }
//		        
//		        if (Loci::MPI_rank==0) cout << "[I] The relative displacement L2 norm computed" << endl ;
//		        //	cout << "[I] The relative displacement L2 norm computed" << endl ;
//		        	
//		        if ((fabs(NORM)<(*FSIIterationTolerance)) && (*itfsi >= *FSIIterationMinimum)) {
//		        	if (Loci::MPI_rank==0) cout << "FSI: " << *ncycle << " itfsi = " << *itfsi << " Residual = " << NORM << " FSI inneriteration converged" << ", norm: " << NORM << endl ;
//		        	*FSIIterationFinishedStar = true ;
//		        } else {
//		        	if (Loci::MPI_rank==0) cout << "FSI: " << *ncycle << " itfsi = " << *itfsi << " Residual = " << NORM << " FSI inneriteration not converged" << ", Tolerance was " << *FSIIterationTolerance << ", minimum number of FSI iterations is " << *FSIIterationMinimum << endl ;
//		 //      	if (Loci::MPI_rank==0) cout << "FSI inneriteration not converged: Residual = " << sumDiff << ", Tolerance was " << *FSIIterationTolerance << endl ;     
//		        //	if (Loci::MPI_rank==0) cout << "FSI inneriteration not converged: CSDnodesCurrent = " << (*CSDnodesCurrent) << endl ;	      	
//		        //	if (Loci::MPI_rank==0) cout << "FSI inneriteration not converged: CSDnodesRef = " << (*CSDnodesRef) << endl ;	      	
//		        //	if (Loci::MPI_rank==0) cout << "FSI inneriteration not converged: CSDdisplacementsStar = " << (*CSDdisplacementsStar) << endl ;	      	
//		        	*FSIIterationFinishedStar = false ;
//		       	 //if (Loci::MPI_rank==0) cout << "FSI: " << *ncycle << " itfsi = " << *itfsi << endl ;	
//		        }	
//		      } else {
//		      	if (Loci::MPI_rank==0) cout << "FSI: " << *ncycle << " itfsi = " << *itfsi << endl ;	
//		      }
//      	} else {
//      		*FSIIterationFinishedStar = false ;
//      		if (Loci::MPI_rank==0) cout << "FSI: Bypassing FSIconvergence check because CFD{n,it-1} is not converged: FSIIterationFinished = " << ((*FSIIterationFinishedStar)?"true":"false") << endl ;
//      	}
//      }	      
//  	} ;
//
//  register_rule<CheckIterationConvergedFSI> registerCheckIterationConvergedFSI ;
///*
//class CheckIterationConvergedFSIdefault : public singleton_rule {
//    private:      
//    	const_param<int> itfsi ;
//      param<bool> FSIIterationFinished ;      
//    public:
//
//      // Define input and output.
//      CheckIterationConvergedFSIdefault() {
//        name_store("FSIIterationFinished{n,it}",FSIIterationFinished) ; 
//        name_store("itfsi{n,it}",itfsi) ;       
//        input("itfsi{n,it}") ;
//        output("FSIIterationFinished{n,it}") ;
//        disable_threading() ;
//      }
//
//      // Check if iteration is finished.
//      void compute(const sequence &seq) {
//				
//		//		if (Loci::MPI_rank==0) cout << "FSI IterationFinished false by default" << endl ;
//				*FSIIterationFinished = false ;
//        
//      }      
//  } ;
//
//  register_rule<CheckIterationConvergedFSIdefault> registerCheckIterationConvergedFSIdefault ;
//*/
}

