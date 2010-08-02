
#include "CSDinput.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib> // for exit function
#include <string>
#include <sstream>
using namespace std;

	
// Loci includes.
#include <Loci.h>

void CSDinput::readMesh(string meshFilename) {
 	
	//Open file 
	if (Loci::MPI_rank==0) cout << "[I] Reading the CSD mesh" << endl;
	ifstream meshData ;
	meshData.open(meshFilename.c_str(), ios::in) ;
	if (!meshData.is_open) {
		if (Loci::MPI_rank==0) cerr << "[E] CSD meshfile file " << meshFilename  << " couldn't be opened" << endl;
		exit(1);
	}
	cout << "[I] ..." << meshFilename << " is successfully opened " << endl;
	
	// initialization
	int counter = 0;
	int tempInt ;
	vector<double, 3> tempVector ;

// Create nodes by reading the mesh until EOF found
	while (meshData >> tempInt >> tempVector[0] >> tempVector[1] >> tempVector[2] )
	{
		csdMesh.push_back(tempVector) ;
		counter++ ;
	}
	
	numberOfNodes = counter ;

	if (numberOfNodes == csdMesh.size() ) {
		if (Loci::MPI_rank==0) cout << "[I] ...Number of nodes  = " << numberOfNodes << endl; 
	} else {
		if (Loci::MPI_rank==0) cout << "[E] ...Something wrong constructing the node size = " << csdMesh.size() << " numofNodes = " << numberOfNodes << endl;
	}
	
	if (Loci::MPI_rank==0) cout << "[I] ...Number of CSD nodes: " << numberOfNodes << endl;

	meshData.close();
	
} // CSDinput::readMesh


void CSDinput::readConnectivity(string connectFilename) {

	//Open file 
	if (Loci::MPI_rank==0) cout << "[I] Reading the CSD connectivity" << endl;
	ifstream connectData ;
	connectData.open(connectFilename.c_str(), ios::in) ;
	if (!connectData.is_open) {
		if (Loci::MPI_rank==0) cerr << "[E] CSD connectivity file " << connectFilename  << " couldn't be opened" << endl;
		exit(1);
	}
	cout << "[I] ..." << connectFilename << " is successfully opened " << endl;
	
	// initialization
	int counter = 0;
	int tempInt ;
	vector<int, 3> tempVector ;

// Create nodes by reading the mesh until EOF found
	while (connectData >> tempInt >> tempInt >> tempVector[0] >> tempVector[1] >> tempVector[2] )
	{
		connectData.push_back(tempVector) ;
		counter++ ;
	}
	
	numberOfElements = counter ;

	if (numberOfElements == connectData.size() ) {
		if (Loci::MPI_rank==0) cout << "[I] ...Number of elements  = " << numberOfElements << endl; 
	} else {
		if (Loci::MPI_rank==0) cout << "[E] ...Something wrong constructing the elements size = " << connectData.size() << " numofElements = " << numberOfElements << endl;
	}
	
	if (Loci::MPI_rank==0) cout << "[I] ...Number of CSD elements: " << numberOfElements << endl;

	connectData.close();
	
} // CSDinput::readConnectivity



void CSDinput::readBC(string bcFilename) {

	//Open file 
	if (Loci::MPI_rank==0) cout << "[I] Reading the CSD boundary conditions" << endl;
	ifstream bcData ;
	bcData.open(bcFilename.c_str(), ios::in) ;
	if (!bcData.is_open) {
		if (Loci::MPI_rank==0) cerr << "[E] CSD connectivity file " << bcFilename  << " couldn't be opened" << endl;
		exit(1);
	}
	cout << "[I] ..." << bcFilename << " is successfully opened " << endl;
	
	// initialization
	int counter = 0 ;
	int tempInt ;
	real tempReal ;

// Create nodes by reading the mesh until EOF found
	while (bcData >> tempInt >> tempReal )
	{
		csdBCdof.push_back(tempInt) ;
		csdBCZeroConstraint.push_back(tempReal) ;
		counter++ ;
	}
	
	numberOfBCs = counter ;

	if (numberOfBCs == csdBCdof.size() ) {
		if (Loci::MPI_rank==0) cout << "[I] ...Number of bcs  = " << numberOfBCs << endl; 
	} else {
		if (Loci::MPI_rank==0) cout << "[E] ...Something wrong constructing the bc size = " << csdBCdof.size() << " numofNodes = " << numberOfBCs << endl;
	}
	
	if (Loci::MPI_rank==0) cout << "[I] ...Number of CSD bcs: " << numberOfBCs << endl;

	bcData.close();
	
} // CSDinput::readBC

