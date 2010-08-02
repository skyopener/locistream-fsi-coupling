#ifndef HYPRE_H
#define HYPRE_H

// Standard library includes.
#include<vector>
using std::vector ;
                                                                                
// Loci includes.
#include<Loci>
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
                                                                                
// HYPRE includes.
#include "HYPRE.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#ifdef HYPRE_RELEASE_NAME // hypre-2.6.0b
#include "_hypre_utilities.h"
#else
#ifdef HYPRE_SMG
#include "_hypre_utilities.h" // hypre-2.0.0
#else
#include "utilities.h" // hypre-1.10.0b
#define HYPRE_1_10_0b
#endif
#endif
                                                                                
// StreamUns includes.
#include "sciTypes.h"
                                                                                
namespace streamUns {
                                                                                
//-----------------------------------------------------------------------------
// Wrapper classes for HYPRE objects. These have been created primarily to
// handle memory management.
                                                                                
  // Wrapper for the HYPRE vector class.
  class HypreVector {
    private:
      bool created ;
      mutable HYPRE_IJVector v ;
      mutable HYPRE_ParVector vPar ;
      int numValue,*index ;
      mutable void * object ;
    public:
      HypreVector() : created(false) {}
      ~HypreVector() {
        if(created){ HYPRE_IJVectorDestroy(v) ; delete [] index ; }
      }
    public:
      void Assemble() const {
        if(HYPRE_IJVectorAssemble(v)){
          cerr << "ERROR: HYPRE_IJVectorAssemble() returned error." << endl ;
          Loci::Abort() ;
        }
        if(HYPRE_IJVectorGetObject(v,&object)){
          cerr << "ERROR: HYPRE_IJVectorGetObject() returned error." << endl ;
          Loci::Abort() ;
        }
        vPar=(HYPRE_ParVector) object ;
      }
      void Create(int rowStart,int rowEnd) {
        if(created) return ;
        if(HYPRE_IJVectorCreate(MPI_COMM_WORLD,rowStart,rowEnd,&v)){
          cerr << "ERROR: HYPRE_IJVectorCreate() returned error." << endl ;
          Loci::Abort() ;
        }
        if(HYPRE_IJVectorSetObjectType(v,HYPRE_PARCSR)){
          cerr << "ERROR: HYPRE_IJVectorSetObjectType() returned error."
            << endl ; Loci::Abort() ;
        }
        numValue=rowEnd-rowStart+1 ; index=new int[numValue] ;
        for(int i=0;i<numValue;++i) index[i]=rowStart+i ;
        created=true ;
      }
      const HYPRE_ParVector& Data() const { return vPar ; }
      void GetValues(double *value) const {
        if(HYPRE_IJVectorGetValues(v,numValue,index,value)){
          cerr << "ERROR: HYPRE_IJVectorGetValues() returned error." << endl ;
          Loci::Abort() ;
        }
        return ;
      }
      void Initialize() const {
        if(HYPRE_IJVectorInitialize(v)){
          cerr << "ERROR: HYPRE_IJVectorInitialize() returned error." << endl ;
          Loci::Abort() ;
        }
      }
      void SetValues(const double *value) const {
        if(HYPRE_IJVectorSetValues(v,numValue,index,value)){
          cerr << "ERROR: HYPRE_IJVectorSetValues() returned error." << endl ;
          Loci::Abort() ;
        }
      }
  } ;

  // Wrapper for the HYPRE matrix class. This class assumes a square matrix.
  class HypreMatrix {
    private:
      bool created ;
      mutable HYPRE_IJMatrix m ;
      mutable HYPRE_ParCSRMatrix mParCSR ;
      mutable void *object ;
    public:
      HypreMatrix() : created(false) {}
      ~HypreMatrix() { if(created) HYPRE_IJMatrixDestroy(m) ; }
    public:
      void Assemble() const {
        if(HYPRE_IJMatrixAssemble(m)){
          cerr << "ERROR: HYPRE_IJMatrixAssemble() returned error." << endl ;
          Loci::Abort() ;
        }
        if(HYPRE_IJMatrixGetObject(m,&object)){
          cerr << "ERROR: HYPRE_IJMatrixGetObject() returned error." << endl ;
          Loci::Abort() ;
        }
        mParCSR=(HYPRE_ParCSRMatrix) object ;
      }
      void Create(int rowStart,int rowEnd,int *numDiagonalNonZero,int
      *numOffDiagonalNonZero) {
        if(created) return ;
        if(HYPRE_IJMatrixCreate(MPI_COMM_WORLD,rowStart,rowEnd,rowStart,
        rowEnd,&m)){
          cerr << "ERROR: HYPRE_IJMatrixCreate() returned error." << endl ;
          Loci::Abort() ;
        }
        if(HYPRE_IJMatrixSetObjectType(m,HYPRE_PARCSR)){
          cerr << "ERROR: HYPRE_IJMatrixSetObjectType() returned error."
            << endl ; Loci::Abort() ;
        }
/* TEMPORARY
        if(HYPRE_IJMatrixSetDiagOffdSizes(m,numDiagonalNonZero,
        numOffDiagonalNonZero)){
          cerr << "ERROR: HYPRE_IJMatrixSetDiagOffdSizes() returned error."
            << endl ; Loci::Abort() ;
        }*/
        created=true ;
      }
      const HYPRE_ParCSRMatrix& Data() const { return mParCSR ; }
      void Initialize() const {
        if(HYPRE_IJMatrixInitialize(m)){
          cerr << "ERROR: HYPRE_IJMatrixInitialize() returned error." << endl ;
          Loci::Abort() ;
        }
      }
      void SetValues(int numRow,int *numColPerRow,int *rowNum,int *colNum,
      double *colValue) const {
/*cout << "numRow,numColPerRow[0]: " << numRow << " " << numColPerRow[0] << endl ;
cout << "rowNum[0]: " << rowNum[0] << endl ;
cout << "colNum" << endl ;
for(int i=0;i<numColPerRow[0];++i)
cout << "i,colNum[i]: " << i << " " << colNum[i] << endl ;
cout << "colValue" << endl ;
for(int i=0;i<numColPerRow[0];++i)
cout << "i,colValue[i]: " << i << " " << colValue[i] << endl ;*/
        if(HYPRE_IJMatrixSetValues(m,numRow,numColPerRow,rowNum,colNum,
        colValue)){
          cerr << "ERROR: HYPRE_IJMatrixSetValues() returned error." << endl ;
          Loci::Abort() ;
        }
      }
  } ;
                                                                                
  // Wrapper for the HYPRE solver. NOTE: For the functions that set the
  // number of grid sweeps and relaxation type, we cannot delete the data
  // after passing to HYPRE since HYPRE holds this as a handle instead of
  // making a copy. Bad design on their part.
  class HypreSolver {
    private:
      bool created ;
      mutable HYPRE_Solver s ;
    public:
      HypreSolver() : created(false) {}
      ~HypreSolver() { if(created) HYPRE_BoomerAMGDestroy(s) ; }
    public:
      void Create() {
        if(created) return ;
        if(HYPRE_BoomerAMGCreate(&s)){
          cerr << "ERROR: HYPRE_BoomerAMGCreate() returned error." << endl ;
          Loci::Abort() ;
        }
//      HYPRE_BoomerAMGSetPrintLevel(s,3);
        created=true ;
      }
      int GetNumIterations() const {
        int num ; HYPRE_BoomerAMGGetNumIterations(s,&num) ; return num ;
      }
      double GetFinalRelativeResidualNorm() const {
        double r ; HYPRE_BoomerAMGGetFinalRelativeResidualNorm(s,&r) ; return r ;
      }
      void SetCoarsenType(int type) const {
        if(HYPRE_BoomerAMGSetCoarsenType(s,type)){
          cerr << "ERROR: HYPRE_BoomerAMGSetCoarsenType() returned error."
            << endl ; Loci::Abort() ;
        }
      }
      void SetCycleNumSweeps(int num0,int num1,int num2,int num3) const {
#ifdef HYPRE_1_10_0b
        HYPRE_BoomerAMGSetCycleNumSweeps(s,num0,0) ;
#endif
        HYPRE_BoomerAMGSetCycleNumSweeps(s,num1,1) ;
        HYPRE_BoomerAMGSetCycleNumSweeps(s,num2,2) ;
        HYPRE_BoomerAMGSetCycleNumSweeps(s,num3,3) ;
      }
      void SetCycleRelaxType(int type0,int type1,int type2,int type3) const {
#ifdef HYPRE_1_10_0b
        HYPRE_BoomerAMGSetCycleRelaxType(s,type0,0) ;
#endif
        HYPRE_BoomerAMGSetCycleRelaxType(s,type1,1) ;
        HYPRE_BoomerAMGSetCycleRelaxType(s,type2,2) ;
        HYPRE_BoomerAMGSetCycleRelaxType(s,type3,3) ;
      }
      void SetMaxIterations(int maxIterations) const {
        if(HYPRE_BoomerAMGSetMaxIter(s,maxIterations)){
          cerr << "ERROR: HYPRE_BoomerAMGSetMaxIter() returned error." << endl ;
          Loci::Abort() ;
        }
      }
      void SetMaxLevels(int maxLevels) const {
        if(HYPRE_BoomerAMGSetMaxLevels(s,maxLevels)){
          cerr << "ERROR: HYPRE_BoomerAMGSetMaxLevels() returned error."
            << endl ; Loci::Abort() ;
        }
      }
      void SetMaxRowSum(real maxRowSum) const {
        if(HYPRE_BoomerAMGSetMaxRowSum(s,maxRowSum)){
          cerr << "ERROR: HYPRE_BoomerAMGSetMaxRowSum() returned error."
            << endl ; Loci::Abort() ;
        }
      }
      void SetNumSweeps(int num) const {
        HYPRE_BoomerAMGSetNumSweeps(s,num) ;
      }
      void SetRelaxType(int type) const {
        HYPRE_BoomerAMGSetRelaxType(s,type) ;
      }
      void SetStrongThreshold(double threshold) const {
        if(HYPRE_BoomerAMGSetStrongThreshold(s,threshold)){
          cerr << "ERROR: HYPRE_BoomerAMGSetStrongThreshold() returned error."
            << endl ; Loci::Abort() ;
        }
      }
      void SetTolerance(real tolerance) const {
        if(HYPRE_BoomerAMGSetTol(s,tolerance)){
          cerr << "ERROR: HYPRE_BoomerAMGSetTol() returned error." << endl ;
          Loci::Abort() ;
        }
      }
      void Solve(const HypreMatrix &A,const HypreVector &b,HypreVector &x)
      const {
        int errorFlag ;
        if(HYPRE_BoomerAMGSetup(s,A.Data(),b.Data(),x.Data())){
          cerr << "ERROR: HYPRE_BoomerAMGSetup() returned error." << endl ;
          Loci::Abort() ;
        }
        if((errorFlag=HYPRE_BoomerAMGSolve(s,A.Data(),b.Data(),x.Data()))){
//        cerr << "ERROR: HYPRE_BoomerAMGSolve() returned error flag: "
//          << errorFlag << endl ;
//        Loci::Abort() ;
        }
//      int numIteration ;
//        HYPRE_BoomerAMGGetNumIterations(s,&numIteration) ;
//cout << "numIteration: " << numIteration << endl ;
      }
  } ;

}

#endif
