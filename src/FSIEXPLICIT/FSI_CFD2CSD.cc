#line 1 "FSI_CFD2CSD.loci"
//-----------------------------------------------------------------------------
// Description: This file contains rules for assembling and solving the linear
//   elasticity equations for the movement of interior nodes given the
//   boundary node displacements.
//-----------------------------------------------------------------------------
                                                                                
// Standard library includes.
#include <map>
using std::map ;
#include <vector>
using std::vector ;
                                                                                
// Loci includes.
#include <Loci.h>
using Loci::Area ;
using Loci::singleton_rule ;
using Loci::pointwise_rule ;
using Loci::unit_rule ;
using Loci::apply_rule ;
using Loci::sequence ;
using Loci::entitySet ;
using Loci::EMPTY ;
using Loci::Entity ;
using Loci::register_rule ;
using Loci::blackbox ;
using Loci::const_blackbox ;
using Loci::param ;
using Loci::const_param ;
using Loci::store ;
using Loci::const_store ;
using Loci::const_Map ;
using Loci::const_multiMap ;
using Loci::const_MapVec ;
using Loci::const_storeVec ;
	
#include <mpi.h>
                                                                                
// StreamUns includes.
#include "const.h"
#include "FSI_move.h"
#include "residual.h"
#include "sciTypes.h"


#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas ;
                                                                                
namespace streamUns {


//-----------------------------------------------------------------------------
// Compute the force on the top and bottom CSD faces
//

 // Pointwise rule to compute the pressure force on the flexible no-slip boundaries at the initial level.
//  class BoundaryForceCSDic : public pointwise_rule {
//    private:      
//      store<vect3d> boundaryForce ;
//    public:
//
//      // Define input and output.
//      BoundaryForceCSDic() {
//        name_store("boundaryForceCSD(X){n,it=0}",boundaryForce) ;
//        output("boundaryForceCSD(X){n,it=0}") ;
//        constraint("X{n}") ;
//        disable_threading() ;
//      }
//
//      // Add the face value to the total.
//      void calculate(Entity face) {				
//				boundaryForce[face] = vect3d(0.0,0.0,0.0) ; //+ viscousForce ;
//      }
//
//      // Loop over the faces.
//      void compute(const sequence &seq) { do_loop(seq,this) ; }
//  } ;
//
//  register_rule<BoundaryForceCSDic> registerBoundaryForceCSDic ;

  // Pointwise rule to compute the pressure force on the flexible no-slip boundaries.
  class BoundaryForceCSD : public pointwise_rule {
    private:
      const_store<real> p_f ;
	const_store<real> boundaryPressure ;
  //    const_store<real> pCorrected ;
      //const_store<vect3d> tauWall ;
      const_store<Area> area ;
      const_Map ci ;
      const_store<tens3d> vGradient ;
      const_store<real> laminarViscosity_f ;
      //const_store<Area> area ;
      const_store<real> noWallFunction ;
      store<vect3d> boundaryForce ;
//	  store<real> boundaryForce ;
    public:

      // Define input and output.
      BoundaryForceCSD() {
 //       name_store("boundaryPressure{n}",boundaryPressure) ;
//		name_store("boundaryPressure{n}",p_f) ;
		name_store("p_f{n}",p_f) ;
   //     name_store("pCorrected{n,it}",pCorrected) ;        
  //      name_store("tauWall",tauWall) ;
        name_store("area{n}",area) ;
        name_store("boundaryForceCSD(X){n}",boundaryForce) ;
        name_store("ci",ci) ;
    //    name_store("gradv3d(v){n,it}",vGradient) ;
    //    name_store("laminarViscosity_f{n,it}",laminarViscosity_f) ;
        //name_store("area",area) ;
     //   name_store("noWallFunction",noWallFunction) ;
    // 		input("ci->pCorrected{n,it}") ;  
        input("p_f{n}") ;
//		 input("boundaryPressure{n}") ;
        input("area{n}") ;
  //      input("tauWall"); 
   //     input("ci->gradv3d(v){n,it},laminarViscosity_f{n,it},noWallFunction") ;
        output("boundaryForceCSD(X){n}") ;
        constraint("X{n}") ;
        disable_threading() ;
      }

      // Add the face value to the total.
      void calculate(Entity face) {
      	int counter=0 ;
				vect3d viscousForce ; //, pressureForce ;				
				real pressureForce ;
     //   pressureForce=pCorrected[ci[face]]*area[face].n*area[face].sada ;
     //   pressureForce=p_f[face]*area[face].n*area[face].sada ;
		pressureForce=p_f[face]*area[face].sada ;
	//	pressureForce=boundaryPressure[face]*area[face].sada ;
    //    cout << "inside p-force computation: r, " << Loci::MPI_rank << ", " << pressureForce << endl ;
        viscousForce=noWallFunction[face]*dotTemp(vGradient[ci[face]]+Transpose(vGradient[ci[face]]),area[face].n*(-area[face].sada)*laminarViscosity_f[face]) ;
        //viscousForce=tauWall[face]*(-area[face].sada) ;
				//boundaryForce[face] = pressureForce ; //+ viscousForce ;
				boundaryForce[face].x = 0.0; //((fabs(area[face].n.x)<5.0e-7)? 0.0 : pressureForce.x) ;
		//		boundaryForce[face].y = ((fabs(pressureForce)<5.0e-7)? 0.0 : pressureForce) ;
				boundaryForce[face].y = pressureForce ;
				boundaryForce[face].z = 0.0; //((fabs(area[face].n.z)<5.0e-7)? 0.0 : pressureForce.z) ;
	//			cout << "rank: " << Loci::MPI_rank << "counter = " << counter << ", boundaryForce = " << boundaryForce[face].y <<  endl ;;
	//			cout << ", p_f " << p_f[face] << ", area: " << area[face].sada << ", n= " << area[face].n.x << ", " << area[face].n.y << ", " << area[face].n.z << endl ;
				++counter ;
      }

      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryForceCSD> registerBoundaryForceCSD ;

  // Rule to assemble the RBF RHS matrix in the x-direction 
  // with the "sXZero" option and the specified boundary displacement nodes
  class AssembleRHSFSI : public pointwise_rule {
    private:
      const_store<vect3d> boundaryForceCSD ;
      store<vect3d> B ;
    public:
                                                                                
      // Define input and output.
      AssembleRHSFSI() {
        name_store("FSI_B(X){n}",B) ;
        name_store("boundaryForceCSD(X){n}", boundaryForceCSD) ;
        input("boundaryForceCSD(X){n}") ;
        output("FSI_B(X){n}") ;
        constraint("X{n}") ;
        constraint("gridMotionTimeDependent{n}") ;
        disable_threading() ;
      }
                                                                                
      // Add relaxation for a single node.
      void calculate(Entity face) {
      int counter = 0 ;
		  B[face] = boundaryForceCSD[face] ;
	//	  cout << "inside B-p-force computation:, rank " << Loci::MPI_rank << ", counter = " << counter << B[face] << endl ;
		  ++counter ;
      }
                                                                                
      // Loop over face.
      virtual void compute(const sequence &seq) { 
    //  	cout << "Inside FSI_B : rank = " << Loci::MPI_rank << endl ;
      	do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<AssembleRHSFSI> registerAssembleRHSRFSI;

//-----------------------------------------------------------------------------
// Rules for setting up the PETSC data.

  // Gets the number of boundary face assigned on all processes. This rule is
  // definitely not in the Loci spirit since we are collecting data from
  // other processes.
//  class FSIGetLocalBoundaryFaceDataUnit : public unit_rule {
//    private:
//      param<vector<int> > fsiNumBoundaryFace ;
//    public:
//
//      // Define input and output.
//      FSIGetLocalBoundaryFaceDataUnit() {
//        name_store("fsiNumBoundaryFace(X){n,it}",fsiNumBoundaryFace) ;
//        output("fsiNumBoundaryFace(X){n,it}") ;
//        constraint("X{n,it}") ;
//        disable_threading() ;
//      }
//
//      // Get the number of nodes for each process.
//      virtual void compute(const sequence &seq) {}
//  } ;
//
//  register_rule<FSIGetLocalBoundaryFaceDataUnit> registerFSIGetLocalBoundaryFaceDataUnit ;

  // Gets the number of boundary nodes assigned on all processes. This rule is
  // definitely not in the Loci spirit since we are collecting data from
  // other processes.
  //class FSIGetLocalBoundaryFaceDataApply : public apply_rule<param<vector<int> >,Loci::NullOp<vector<int> > > {
  class FSIGetLocalBoundaryFaceDataApply : public singleton_rule {
    private:
      param<vector<int> > fsiNumBoundaryFace ;
    public:

      // Define input and output.
      FSIGetLocalBoundaryFaceDataApply() {
        name_store("fsiNumBoundaryFace(X){n}",fsiNumBoundaryFace) ;
        output("fsiNumBoundaryFace(X){n}") ;
        constraint("X{n}") ;
//        constraint("RBFxConstraints") ;
        disable_threading() ;
      }

      // Get the number of nodes for each process.
      virtual void compute(const sequence &seq) {

		const int rank=Loci::MPI_rank ;
		const int p=Loci::MPI_processes ;
			
	//		cout << "FSInum: starting, r = " << rank << endl ;

		entitySet fsiLocalBoundaryFace ;
        // Get the collection of entities assigned to this processor
        Loci::storeRepP myEntities=Loci::exec_current_fact_db->get_variable
          ("my_entities") ;
        entitySet localEntities=~EMPTY ;
        if(myEntities!=0) localEntities=(*myEntities).domain() ;

        // Get the local nodes.
        fsiLocalBoundaryFace=entitySet(seq) ;

				(*fsiNumBoundaryFace).resize(p) ;

        // Distribute the number of nodes to all processes.
        *fsiNumBoundaryFace=Loci::all_collect_sizes((fsiLocalBoundaryFace).size()) ;
        	for (int i=0;i<p;++i) {
    //    	cout << "FSInum: r, fsiNumBF" << rank << ", " << (*fsiNumBoundaryFace)[i] << endl ;
        	}
      }
  } ;

  register_rule<FSIGetLocalBoundaryFaceDataApply> registerFSIGetLocalBoundaryFaceDataApply ;

  // Creates the node-to-row map. Each boundary node now has a unique row index
  class FSIBoundaryFaceToRow : public pointwise_rule {
    private:
      const_param<vector<int> > fsiNumBoundaryFace ;
      store<int> fsiBoundaryFaceToRow ;
    public:

      // Define input and output.
      FSIBoundaryFaceToRow() {
        name_store("fsiNumBoundaryFace(X){n}",fsiNumBoundaryFace) ;
        name_store("fsiBoundaryFaceToRow(X){n}",fsiBoundaryFaceToRow) ;
        input("fsiNumBoundaryFace(X){n}") ;
        output("fsiBoundaryFaceToRow(X){n}") ;
        constraint("X{n}") ;
//        constraint("FSIxConstraints") ;
        disable_threading() ;
      }

      // Set the Petsc row for each node.
      virtual void compute(const sequence & seq) {

        // Compute the row offset for this process.
//        int offset=(*fsiLocalStart)[Loci::MPI_rank] ;
        int offset=0 ;
        for(int i=0;i<Loci::MPI_rank;++i){ offset+=(*fsiNumBoundaryFace)[i] ; }

        // Assign row number.
        sequence::const_iterator nodePtr=seq.begin() ;
        for(int i=0;i<(*fsiNumBoundaryFace)[Loci::MPI_rank];++nodePtr,++i){
          fsiBoundaryFaceToRow[*nodePtr]=offset+i ;
        }
      }
  } ;

  register_rule<FSIBoundaryFaceToRow> registerFSIBoundaryFaceToRow ;
  
//-----------------------------Q---------------------------------------------------------
//

  // Sets up the Q vector.
//  class FSIBoundaryFaceAssembleQUnit : public unit_rule { // see interpolateFile.cc
//    private:
//      const_param<vector<int> > fsiNumBoundaryFace ;
//	  	const_store<vect3d> faceCenter ;
//	  	blackbox<vector<real> > fsiQ ;
//    public:
//
//      // Define input and output.
//      	FSIBoundaryFaceAssembleQUnit() {
//        name_store("fsiNumBoundaryFace(X){n,it}",fsiNumBoundaryFace) ;
//				name_store("facecenter{n,it-1}",faceCenter) ;
//        name_store("fsiBoundaryFaceQ(X){n,it}", fsiQ) ;
//        input("fsiNumBoundaryFace(X){n,it}") ;
//        input("facecenter{n,it-1}") ;
//        output("fsiBoundaryFaceQ(X){n,it}") ;
//				constraint("X{n,it}");
//        disable_threading() ;
//      }
//
//      // Do the set-up.
//      void compute(const sequence & seq) {
//      	cout << "Inside Qunit: rank = " << Loci::MPI_rank << endl ;
//      }
//  } ;
 
//  register_rule<FSIBoundaryFaceAssembleQUnit> registerFSIBoundaryFaceAssembleQUnit ;

  // Assemble Q -> unit rule to output blackbox
  //class FSIBoundaryFaceAssembleQApply : public apply_rule<blackbox<vector<real> >,Loci::NullOp<blackbox<vector<real> > > > {
  
  
// $type fsiNumBoundaryFace(X0) param<vector<int> >   
// $type facecenter store<vect3d> 
// $type fsiBoundaryFaceQ(X0) blackbox<vector<real> > 

//$rule blackbox(fsiBoundaryFaceQ(X){n}<-fsiNumBoundaryFace(X){n},facecenter{n}),constraint(X{n}), prelu
    // Sets up the Q vector.
  class FSIBoundaryFaceAssembleQUnit : public unit_rule { // see interpolateFile.cc
    private:
		const_param<vector<int> > fsiNumBoundaryFace ;
		const_store<vect3d> faceCenter ;
		blackbox<vector<real> > fsiQ ;
    public:

      // Define input and output.
        FSIBoundaryFaceAssembleQUnit() {
        name_store("fsiNumBoundaryFace(X){n}",fsiNumBoundaryFace) ;
        name_store("facecenter{n}",faceCenter) ;
        name_store("fsiBoundaryFaceQ(X){n}", fsiQ) ;
        input("fsiNumBoundaryFace(X){n}") ;
        input("facecenter{n}") ;
        output("fsiBoundaryFaceQ(X){n}") ;
        constraint("X{n}");
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {
   //     cout << "Inside Qunit: rank = " << Loci::MPI_rank << endl ;
      }
  } ;
 
  register_rule<FSIBoundaryFaceAssembleQUnit> registerFSIBoundaryFaceAssembleQUnit ;

  
  // Assemble Q -> unit rule to output blackbox
  class FSIBoundaryFaceAssembleQApply : public apply_rule<blackbox<vector<real> >,Loci::NullOp<blackbox<vector<real> > > > {
  //class FSIBoundaryFaceAssembleQApply : public singleton_rule {
    private:
      const_param<vector<int> > fsiNumBoundaryFace ;
	  	const_store<vect3d> faceCenter ;
	  	const_param<bool> CFDIterationFinished ;
	  	blackbox<vector<real> > fsiQ ;
    public:

      // Define input and output.
      FSIBoundaryFaceAssembleQApply() {
        name_store("fsiNumBoundaryFace(X){n}",fsiNumBoundaryFace) ;
				name_store("facecenter{n}",faceCenter) ;
        name_store("fsiBoundaryFaceQ(X){n}", fsiQ) ;
    //    name_store("CFDIterationFinished{n,it-1}",CFDIterationFinished) ;
    //		input("CFDIterationFinished{n,it-1}") ;
        input("fsiNumBoundaryFace(X){n}") ;
        input("facecenter{n}") ;
        output("fsiBoundaryFaceQ(X){n}") ;
				constraint("X{n}");
        disable_threading() ;
      }

      // Do the set-up.
      void compute(const sequence & seq) {

        // Get the number of local and global nodes.
       // if (Loci::MPI_rank==0) cout << "FSIBoundaryFaceAssembleQApply" << endl ;
		const int p = Loci::MPI_processes ;
		const int rank = Loci::MPI_rank ;
        int globalNumFace=0 ;
        for(unsigned int i=0;i<(*fsiNumBoundaryFace).size();++i) globalNumFace+=(*fsiNumBoundaryFace)[i] ;
     
		
		// Compute the row offset for this process.
        int localStart=0 ; 
        for(int i=0;i<rank;++i) localStart+=(*fsiNumBoundaryFace)[i]  ;  // 3 coordinates
		int localNumFace=(*fsiNumBoundaryFace)[rank]; 
			
		//	cout << "FSIQ: r, localStart, globalNumFace=" << rank << ", " << localStart << ", " << globalNumFace << endl ;
					
		vector<int> recvcounts(p,0) ;
		vector<int> displs(p,0) ;
		for(int i=0;i<p;++i) {
			for(int j=0;j<i;++j) displs[i]+=(3*(*fsiNumBoundaryFace)[j]);
			recvcounts[i] = 3*(*fsiNumBoundaryFace)[i] ;
		}
				
		// Allocate fsiQ
		(*fsiQ).resize(3*globalNumFace) ; // number of nodes * 3 coordinates
//		(*fsiQ).resize(3*(*fsiGlobalNumFace)) ; // number of nodes * 3 coordinates

		// Fill in fsiQ
        sequence::const_iterator facePtr=seq.begin() ;
//		for(int i=localStart; i<localStart+localNum;++nodePtr, ++i) {
//		int counter=(*fsiLocalStart)[rank];
		int counter=localStart;
		for(facePtr=seq.begin(); facePtr!=seq.end();++facePtr, ++counter) { // Q[row,col]=Q[3*row+col], Q=[0,1,2;3,4,5;...]
			(*fsiQ)[3*counter+0] = faceCenter[*facePtr].x ;
			(*fsiQ)[3*counter+1] = faceCenter[*facePtr].y ;
			(*fsiQ)[3*counter+2] = faceCenter[*facePtr].z ;
		}

		// Allgatherv
//		for(int i=0; i<(*fsiQ).size(); ++i) {
//			cout << "before gather: p,r,fsiQ[" << i << "]: " << p << "," << rank << ", " << (*fsiQ)[i] << endl ;
//		}
	//	if (rank==0) cout << "Before CFD2CSD Q all gather" << endl ;
		MPI_Allgatherv(&(*fsiQ)[3*(localStart)], 3*localNumFace, MPI_DOUBLE, &(*fsiQ)[0], &recvcounts[0], &displs[0], MPI_DOUBLE, MPI_COMM_WORLD);
	//	if (rank==0) cout << "After CFD2CSD Q all gather" << endl ;
	//if (rank==0) {
	//	for(int i=0;i<globalNumFace;++i) cout << "CFD2CSD Q: r, face centers: " << rank << ", " << (*fsiQ)[3*i+0] << ", " << (*fsiQ)[3*i+1] << ", " << (*fsiQ)[3*i+2] << endl ;
	//}
	//if (rank==0) cout << "FSIQ = " << (*fsiQ) << endl ;
  }
  } ;

  register_rule<FSIBoundaryFaceAssembleQApply> registerFSIBoundaryFaceAssembleQApply ;

}





