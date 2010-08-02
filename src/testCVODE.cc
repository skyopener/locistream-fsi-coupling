//-----------------------------------------------------------------------------
// Description: This file contains rules for testing CVODE with a specific
//   state.
//
// Authors: Siddharth Thakur and Jeff Wright
//          Streamline Numerics Inc.
//          Gainesville, Florida 32609
//-----------------------------------------------------------------------------
                                                                                
// Standard library includes.
#include <vector>
using std::vector ;
                                                                                
// Loci includes.
#include <Loci.h>
                                                                                
// Fluid physics library includes.
#include "eos.h"
#include "reaction.h"
using fluidPhysics::EOS ;
using fluidPhysics::reaction ;
                                                                                
// StreamUns includes.
#include "sciTypes.h"
#include "varsFileInputs.h"
                                                                                
// CVODE includes.
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <cvode/cvode_spgmr.h>
#include <cvode/cvode_spbcgs.h>
#include <cvode/cvode_sptfqmr.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
                                                                                
namespace streamUns {

  // This class is used to pass data to CVODE.
  class CVODEReactionData {
    public:
      int numSpecies ;
      const EOS &eos ;
      real rho,P,T,h ;
      const reaction &reactor ;
      reaction::rates *reactionRate ;
      tmp_array<real> rhoTemp,yTemp,w,dwdr ;
      float *hint ;
    public:
      CVODEReactionData(int numSpecies,const EOS &eos,real rho,real P,real T,
      real h,const reaction &reactor) : numSpecies(numSpecies),eos(eos),
      rho(rho),P(P),T(T),h(h),reactor(reactor),rhoTemp(numSpecies),
      yTemp(numSpecies),w(numSpecies),dwdr(numSpecies) {
        reactionRate=new reaction::rates[reactor.num_rates()] ;
        hint=new float[numSpecies+2] ;
      }
                                                                                
      ~CVODEReactionData() {
        delete [] reactionRate ; delete [] hint ;
      }
  } ;

  // The right-hand side function called by the CVODE integrator.
  static int cvodeRHS(realtype t,N_Vector y,N_Vector ydot,void *f_data) {

    // Recast the reaction data.
    CVODEReactionData *reactionData=(CVODEReactionData*)f_data ;
    int n=reactionData->numSpecies ;

    // The values in the argument "y" of this function are species masses,
    // but we must convert to mass fraction to pass to the EOS inversion
    // function.
    for(int i=0;i<n;++i){
      reactionData->rhoTemp[i]=NV_Ith_S(y,i) ;
      reactionData->yTemp[i]=NV_Ith_S(y,i)/reactionData->rho ;
    }

    // Solve the EOS given rho, h and y to get T and P.
    reactionData->hint[0]=reactionData->T ;
    reactionData->hint[1]=reactionData->P ;
    EOS::State eosState=reactionData->eos.State_from_rho_h(reactionData->
      rhoTemp,reactionData->h,reactionData->hint) ;
    reactionData->P=eosState.pressure() ;
    reactionData->T=eosState.temperature() ;

    // Copute the reaction rates based on the new temperature and pressure.
    reactionData->reactor.extract_rates(reactionData->reactionRate,
      reactionData->eos,reactionData->T) ;

    // Compute the production rates based on the new mass fractions that
    // have been passed into this function.
    reactionData->reactor.compute_w(reactionData->w,reactionData->reactionRate,
      reactionData->yTemp,eosState,reactionData->eos) ;

    // Copy computed production rate to CVODE data structure.
    for(int i=0;i<n;++i) NV_Ith_S(ydot,i)=reactionData->w[i] ;

    // Return success.
    return 0 ;
  }

  // The right-hand side function called by the CVODE integrator.
  static int cvodeJacobian(long int N,DenseMat J,realtype t,N_Vector y,N_Vector
  fy,void *jac_data,N_Vector tmp1,N_Vector tmp2,N_Vector tmp3) {

    // Recast the reaction data.
    CVODEReactionData *reactionData=(CVODEReactionData*)jac_data ;
    int n=reactionData->numSpecies ;

    // The values in the argument "y" of this function are species masses,
    // but we must convert to mass fraction to pass to the EOS inversion
    // function.
    for(int i=0;i<n;++i){
      reactionData->rhoTemp[i]=NV_Ith_S(y,i) ;
      reactionData->yTemp[i]=NV_Ith_S(y,i)/reactionData->rho ;
    }

    // Solve the EOS given rho, h and y to get T and P.
    reactionData->hint[0]=reactionData->T ;
    reactionData->hint[1]=reactionData->P ;
    EOS::State eosState=reactionData->eos.State_from_rho_h(reactionData->
      rhoTemp,reactionData->h,reactionData->hint) ;
    reactionData->P=eosState.pressure() ;
    reactionData->T=eosState.temperature() ;

    // Copute the reaction rates based on the new temperature and pressure.
    reactionData->reactor.extract_rates(reactionData->reactionRate,
      reactionData->eos,reactionData->T) ;

    // Compute the Jacobian entries for each species based on the new mass
    // fractions that have been passed into this function.
    for(int i=0;i<n;++i){
      reactionData->reactor.compute_dwdrs(reactionData->dwdr,i,reactionData->
        reactionRate,0.0,NULL,reactionData->yTemp,eosState,reactionData->eos) ;
      for(int j=0;j<n;++j) DENSE_ELEM(J,i,j)=reactionData->dwdr[j] ;
    }

    // Return success.
    return 0 ;
  }

  // Call the mechanism library to compute the species production rates.
  class TestCVODE : public singleton_rule {
    private:
      const_param<real> cvodeTolerance ;
      const_param<EOS> eos ;
      const_param<reaction> reactor ;
      const_param<int> numSpecies ;
      param<int> testCVODE ;
    private:
      int flag ;
      realtype relativeTolerance ;
      N_Vector absoluteTolerance,rhoNew,rhoOld ;
      realtype time ;
      void *cvodeMemory ;
      reaction::rates *reactionRate ;
    public:

       // Define input and output.
      TestCVODE() {
        name_store("cvodeTolerance",cvodeTolerance) ;
        name_store("eos",eos) ;
        name_store("reactor",reactor) ;
        name_store("numSpecies",numSpecies) ;
        name_store("testCVODE",testCVODE) ;
        input("cvodeTolerance,eos,reactor,numSpecies") ;
        output("testCVODE") ;
      }

      // Compute the species production rates for a single cell.
      void calculate() {

        // Set the state to test.
        real integrationTime=1.0e-04,p=101325.0,temperature=1388.88 ;
        real y[6],rhoi[6] ;
        y[0]=0.00689,y[1]=0.99311,y[2]=y[3]=y[4]=y[5]=0.0 ;

        // Call the EOS to get the density and enthalpy, which will remain
        // constant during the ODE integration.
        EOS::State eosState=eos->State_from_mixture_p_T(y,p,temperature,
          (float *)NULL) ;
        real rho=eosState.density(),h=eosState.enthalpy() ;
        cout << "rho,h: " << rho << " " << h << endl ;

        // Compute the species densities.
        for(int i=0;i<6;++i) rhoi[i]=rho*y[i] ;

        // Initialize the old and new species densities. Set the absolute
        // error tolerance based on the old species mass fractions.
        for(int i=0;i<*numSpecies;++i){
          NV_Ith_S(rhoOld,i)=rho*y[i] ;
          NV_Ith_S(rhoNew,i)=NV_Ith_S(rhoOld,i) ;
          NV_Ith_S(absoluteTolerance,i)=(y[i]==0.0)? 1.0e-20:0.0 ;
//        NV_Ith_S(absoluteTolerance,i)=(y[i]==0.0)? 1.0e-07:0.0 ;
        }

        // Allocate and assign the reaction data to be used by CVODE.
        CVODEReactionData f_data(*numSpecies,*eos,rho,p,temperature,h,
          *reactor) ;

        // Set up CVODE.
        if(cvodeMemory==NULL){

          // Create solver memory and request BDF with Newton iteration.
          cvodeMemory=CVodeCreate(CV_BDF,CV_NEWTON) ;
          if(cvodeMemory==NULL){
            cerr << "ERROR: CVodeCreate() failed!" << endl ; Loci::Abort() ;
          }

          // Allocate CVODE memory.
          flag=CVodeMalloc(cvodeMemory,cvodeRHS,0.0,rhoOld,CV_SV,
            relativeTolerance,absoluteTolerance) ;
          if(flag!=CV_SUCCESS){
            cerr << "ERROR: CVodeMalloc() failed!" << endl ; Loci::Abort() ;
          }
        }else{

          // Reinitialize CVODE solver.
          flag=CVodeReInit(cvodeMemory,cvodeRHS,0.0,rhoOld,CV_SV,
            relativeTolerance,absoluteTolerance) ;
          if(flag!=CV_SUCCESS){
            cerr << "ERROR: CVodeReInit() failed!" << endl ; Loci::Abort() ;
          }
        }

        // Set data to be used by RHS function
        CVodeSetFdata(cvodeMemory,&f_data) ;

        // Choose the linear solver.
        flag=CVDense(cvodeMemory,*numSpecies) ;
        if(flag!=CVDENSE_SUCCESS){
          cerr << "ERROR: CVDense() failed!" << endl ; Loci::Abort() ;
        }
                                                                                
        /* We will provide an analytical Jacobian now.
        flag=CVDenseSetJacFn(cvodeMemory,cvodeJacobian,&f_data) ;
        if(flag!=0){
          cerr << "ERROR: CVDenseSetJacFn() failed!" << endl ; Loci::Abort() ;
        }*/

        // Call CVODE to integrate the system of ODEs to get the new species
        // densites (rhoM*y[i]).
        time=0.0 ;
        while(time<integrationTime){
          flag=CVode(cvodeMemory,integrationTime,rhoNew,&time,CV_ONE_STEP) ;
          cout << "time,P,T: " << time << " " << f_data.P << " " <<
            f_data.T << endl ;
        }
        if(flag==CV_TOO_MUCH_WORK){
          cerr << "ERROR: CVODE could not reach end of itegration time!"
            << endl ;
          cerr << "integrationTime: " << integrationTime << endl ;
        }else if(flag==CV_TOO_MUCH_WORK){
          cerr << "ERROR: CVODE could not satisfy itegration accuracy!"
            << endl ;
        }else if(flag==CV_ERR_FAILURE){
          cerr << "ERROR: Too many error test failures in CVODE!" << endl ;
        }else if(flag==CV_CONV_FAILURE){
          cerr << "ERROR: Too many convergence test failures in CVODE!"
            << endl ;
        }

        // Get output information from the solver.
        long int numSolverSteps,numFCall,numJCall,numNewtonIterations ;
        realtype lastStepSize,currentTime,toleranceScaleFactor ;
        CVodeGetNumSteps(cvodeMemory,&numSolverSteps) ;
        CVodeGetNumRhsEvals(cvodeMemory,&numFCall) ;
        CVDenseGetNumJacEvals(cvodeMemory,&numJCall) ;
        CVodeGetNumNonlinSolvIters(cvodeMemory,&numNewtonIterations) ;
        CVodeGetLastStep(cvodeMemory,&lastStepSize) ;
        CVodeGetCurrentTime(cvodeMemory,&currentTime) ;
        CVodeGetTolScaleFactor(cvodeMemory,&toleranceScaleFactor) ;
        cerr << "integrationTime: " << integrationTime << endl ;
        cerr << "Number of internal solver steps: " << numSolverSteps
          << endl ;
        cerr << "Number of call to f: " << numFCall << endl ;
        cerr << "Number of call to J: " << numJCall << endl ;
        cerr << "Number of Newton iterations: " << numNewtonIterations
          << endl ;
        cerr << "Last internal step size: " << lastStepSize << endl ;
        cerr << "Current time reached: " << currentTime << endl ;
        cerr << "Suggested relative tolerance scaling factor: "
          << toleranceScaleFactor << endl ;
        cerr << "rho,p,t,species: " << rho << " " << p << " " << temperature
          << " " ;
        for(int i=0;i<*numSpecies;++i) cerr << y[i] << " " ;
        cerr << endl ;
      }

      // Assign species mass fractions for a sequence of cells.
      void compute(const sequence &seq) {

        // Set the relative tolerance.
        relativeTolerance=(*cvodeTolerance) ;

        // Initialize the variables required by CVODE.
        absoluteTolerance=rhoNew=rhoOld=NULL ; cvodeMemory=NULL ;
        absoluteTolerance=N_VNew_Serial(*numSpecies) ;
        if(absoluteTolerance==NULL){
          cerr << "ERROR: Inside TestCVODE()" << endl ;
          cerr << "Problem allocating memory for absoluteTolerance" << endl ;
          Loci::Abort() ;
        }
        rhoNew=N_VNew_Serial(*numSpecies) ;
        if(rhoNew==NULL){
          cerr << "ERROR: Inside TestCVODE()" << endl ;
          cerr << "Problem allocating memory for rhoNew" << endl ;
          Loci::Abort() ;
        }
        rhoOld=N_VNew_Serial(*numSpecies) ;
        if(rhoOld==NULL){
          cerr << "ERROR: Inside TestCVODE()" << endl ;
          cerr << "Problem allocating memory for rhoOld" << endl ;
          Loci::Abort() ;
        }

        // Do the test.
        calculate() ;

        // Deallocate memory associated with CVODE.
        N_VDestroy_Serial(absoluteTolerance) ; N_VDestroy_Serial(rhoNew) ;
        N_VDestroy_Serial(rhoOld) ;
        if(cvodeMemory) CVodeFree(&cvodeMemory) ;
      }
  } ;

  register_rule<TestCVODE> registerTestCVODE ;

}
