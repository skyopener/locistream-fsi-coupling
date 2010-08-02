// Program to test the assembly of a global 1-D array from local
// 1-D arrays on each process.

// Standard library includes.
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using std::cerr ;
using std::cout ;
using std::endl ;
using std::fstream ;
using std::ifstream ;
using std::ios ;
using std::ostream ;
using std::string ;
using std::vector ;

// MPI includes.
#include "mpi.h"

int main(int argC,char *argV[]) {

  // Variables.
  int numProcess,myRank,n ;
  float *xLocal,*xGlobal ;

  // Initialize MPI.
  MPI_Init(&argC,&argV) ;

  // Get rank and number of processes.
  MPI_Comm_rank(MPI_COMM_WORLD,&myRank) ;
  MPI_Comm_size(MPI_COMM_WORLD,&numProcess) ;

  // Allocate the local array and set the data.
  n=myRank+1 ; xLocal=new float[n] ;
  for(unsigned int i=0;i<n;++i) xLocal[i]=float(i) ;

  // Send/get the size of the local array to/from all other processes.
  int *localArraySize=new int[numProcess] ;
  {
    int *sendBuffer=new int[numProcess] ;
    for(unsigned int i=0;i<numProcess;++i) sendBuffer[i]=n ;
    MPI_Alltoall(sendBuffer,1,MPI_INT,localArraySize,1,MPI_INT,MPI_COMM_WORLD) ;
    delete [] sendBuffer ;
  }

  // Allocate memory for the global array.
  int nGlobal=0 ;
  for(unsigned int i=0;i<numProcess;++i) nGlobal+=localArraySize[i] ;
  cout << "proc,nGobal: " << myRank << " " << nGlobal << endl ;

  // Deallocate memory.
  delete [] xLocal ;

  // Finalize MPI.
  MPI_Finalize() ;

}
