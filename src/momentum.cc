//-----------------------------------------------------------------------------
// Description: This file contains rules for the momentum equation, including
//    rules for velocity boundary conditions.
//
// Authors: Siddharth Thakur and Jeff Wright
//          Streamline Numerics Inc.
//          Gainesville, Florida 32609
//-----------------------------------------------------------------------------

// Standard library includes.
#include <map>
#include <vector>
using std::map ;
using std::vector ;

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;

// StreamUns includes.
#include "const.h"
#include <readGrid.h>
#include "referenceFrame.h"
#include "residual.h"
#include "sciTypes.h"
#include "varsFileInputs.h"

namespace streamUns {

//-----------------------------------------------------------------------------
// Pre-computed quantities.

  // Momentum total viscosity for the interior faces.
  class ViscosityInterior : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<real> laminarViscosity,eddyViscosity ;
      store<real> viscosity ;
    public:
                                                                                
      // Define input and output.
      ViscosityInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("laminarViscosity",laminarViscosity) ;
        name_store("eddyViscosity",eddyViscosity) ;
        name_store("viscosity",viscosity) ;
        input("(cl,cr)->(laminarViscosity,eddyViscosity)") ;
        output("viscosity") ;
        constraint("internalFaces") ;
      }
                                                                                
      // Compute for each face.
      void calculate(Entity face) {
        viscosity[face]=0.5*(laminarViscosity[cl[face]]+laminarViscosity
          [cr[face]]+eddyViscosity[cl[face]]+eddyViscosity[cr[face]]) ;
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<ViscosityInterior> registerViscosityInterior ;

  // Momentum total viscosity for the boundary faces.
  class ViscosityBoundary : public pointwise_rule {
    private:
      const_store<real> laminarViscosity_f,eddyViscosity_f ;
      store<real> viscosity ;
    public:
                                                                                
      // Define input and output.
      ViscosityBoundary() {
        name_store("laminarViscosity_f",laminarViscosity_f) ;
        name_store("eddyViscosity_f",eddyViscosity_f) ;
        name_store("viscosity",viscosity) ;
        input("laminarViscosity_f,eddyViscosity_f") ;
        output("viscosity") ;
        constraint("boundaryFaces") ;
      }
                                                                                
      // Compute for each face.
      void calculate(Entity face) {
        viscosity[face]=laminarViscosity_f[face]+eddyViscosity_f[face] ;
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<ViscosityBoundary> registerViscosityBoundary ;

//-----------------------------------------------------------------------------
// Rules to set the time integrator factors. These factors are used to provide
// a unified implementation of the BDF, BDF2 and CN schemes as well as to turn
// on the BDF scheme for the first time-step when using BDF2.

  // Creates the time integrator constraints. Note that we have now brought the
  // time step(s) into the factors in order to simplify coding that allows for
  // a variable (time-dependent) time step.
  class TimeIntegratorFactors : public singleton_rule {
    private:
      const_param<int> n ;
      const_param<real> dtOld,dt ;
      const_param<string> timeIntegrator ;
      param<real> timeIntegratorFactor0,timeIntegratorFactor1 ;
      param<real> timeIntegratorFactor2 ;
    public:
                                                                                
      // Define input and output.
      TimeIntegratorFactors() {
        name_store("$n{n}",n) ;
        name_store("dt{n-1}",dtOld) ;
        name_store("dt{n}",dt) ;
        name_store("timeIntegrator{n}",timeIntegrator) ;
        name_store("timeIntegratorFactor0{n}",timeIntegratorFactor0) ;
        name_store("timeIntegratorFactor1{n}",timeIntegratorFactor1) ;
        name_store("timeIntegratorFactor2{n}",timeIntegratorFactor2) ;
        input("$n{n},dt{n-1},dt{n},timeIntegrator{n}") ;
        output("timeIntegratorFactor0{n},timeIntegratorFactor1{n}") ;
        output("timeIntegratorFactor2{n}") ;
      }
                                                                                
      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        if(*timeIntegrator=="BDF" || *timeIntegrator=="CN"){
          *timeIntegratorFactor0=1.0/(*dt) ; *timeIntegratorFactor1=0.0 ;
          *timeIntegratorFactor2=0.0 ;
        }else{
          if((*n)!=0){
            real P=(*dt)+(*dtOld),Q=P*P/((*dt)*(*dt)) ;
            *timeIntegratorFactor0=Q/((*dt)*Q-P) ;
            *timeIntegratorFactor1=1.0/((*dt)*Q-P) ;
            *timeIntegratorFactor2=1.0 ;
          }else{
            *timeIntegratorFactor0=1.0/(*dt) ; *timeIntegratorFactor1=0.0 ;
            *timeIntegratorFactor2=0.0 ;
          }
        }
      }
  } ;
                                                                                
  register_rule<TimeIntegratorFactors> registerTimeIntegratorFactors ;

//-----------------------------------------------------------------------------
// Rule to set the two-dimension factor which is needed to zero out the
// w velocity source term and correction.

  class TwoDimensionFactor : public singleton_rule {
    private:
      const_param<grid_options> grid_file_info ;
      param<real> twoDimensionFactor ;
    public:

      // Define input and output.
      TwoDimensionFactor() {
        name_store("grid_file_info",grid_file_info) ;
        name_store("twoDimensionFactor",twoDimensionFactor) ;
        input("grid_file_info") ;
        output("twoDimensionFactor") ;
      }

      // Set paramter value.
      void compute(const sequence &seq) {
        if((*grid_file_info).optionExists("axisymmetric") ||
        (*grid_file_info).optionExists("pieSlice")){
          *twoDimensionFactor=0.0 ;
        }else{
          *twoDimensionFactor=1.0 ;
        }
      }
  } ;

  register_rule<TwoDimensionFactor> registerTwoDimensionFactor ;

//-----------------------------------------------------------------------------
// Rules to process gravity options from the .vars file.

  // Creates the acceleration and density reference parameters.
  class GravityParameters : public singleton_rule {
    private:
      const_param<GravityOptions> gravity ;
      param<vect3d> g ;
      param<real> rho ;
    public:

      // Define input and output.
      GravityParameters() {
        name_store("gravity",gravity) ;
        name_store("gravityAcceleration",g) ;
        name_store("gravityRhoRef",rho) ;
        input("gravity") ;
        output("gravityAcceleration,gravityRhoRef") ;
      }

      // Set paramter values.
      void compute(const sequence &seq) {
        *g=vect3d(9.8,0.0,0.0) ; *rho=0.0 ;
        if((*gravity).optionExists("g"))
          get_vect3dOption(*gravity,"g","m/s/s",*g) ;
        if((*gravity).optionExists("rhoref"))
          (*gravity).getOption("rhoref",*rho) ;
      }
  } ;

  register_rule<GravityParameters> registerGravityParameters ;

//-----------------------------------------------------------------------------
// Rules to process momentum equation options from the .vars file.

  // Creates the momentum equation solver constraints.
  class MomentumEquationSolverConstraints : public constraint_rule {
    private:
      const_param<MomentumEquationOptions> momentumEquationOptions ;
      Constraint vSGSLinearSolver ;
    public:

      // Define input and output.
      MomentumEquationSolverConstraints() {
        name_store("momentumEquationOptions",momentumEquationOptions) ;
        name_store("vSGSLinearSolver",vSGSLinearSolver) ;
        input("momentumEquationOptions") ;
        output("vSGSLinearSolver") ;
      }

      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        if((*momentumEquationOptions).optionExists("linearSolver")){
          Loci::option_value_type optionValueType=momentumEquationOptions->
            getOptionValueType("linearSolver") ;
          switch(optionValueType){
            case Loci::NAME:
              {
                Loci::option_values optionValues=momentumEquationOptions->
                  getOption("linearSolver") ;
                string name ; optionValues.get_value(name) ;
                if(name=="SGS"){
                  vSGSLinearSolver=~EMPTY ;
                }else{
                  cerr << "Bad linearSolver for momentumEquation." << endl ;
                  Loci::Abort() ;
                }
              }
              break ;
            default:
              cerr << "Bad type for linearSolver in momentumEquation." << endl ;
              Loci::Abort() ;
          }
        }else{
          vSGSLinearSolver=~EMPTY ;
        }
      }
  } ;

  register_rule<MomentumEquationSolverConstraints>
    registerMomentumEquationSolverConstraints ;

  // Creates the momentum equation solver parameters.
  class MomentumEquationSolverParameters : public singleton_rule {
    private:
      const_param<MomentumEquationOptions> momentumEquationOptions ;
      param<int> vMaxIterations ;
      param<real> vRelaxationFactor ;
    public:

      // Define input and output.
      MomentumEquationSolverParameters() {
        name_store("momentumEquationOptions",momentumEquationOptions) ;
        name_store("vMaxIterations",vMaxIterations) ;
        name_store("vRelaxationFactor",vRelaxationFactor) ;
        input("momentumEquationOptions") ;
        output("vMaxIterations,vRelaxationFactor") ;
      }

      // Set up the parameters.
      virtual void compute(const sequence& seq) {

        // Relaxation factor.
        if((*momentumEquationOptions).optionExists("relaxationFactor")){
          Loci::option_value_type optionValueType=momentumEquationOptions->
            getOptionValueType("relaxationFactor") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                momentumEquationOptions->getOption("relaxationFactor",temp) ;
                if(temp<=0.0 || temp>1.0){
                  cerr << "Bad relaxationFactor for momentumEquation." << endl ;
                  Loci::Abort() ;
                }
                *vRelaxationFactor=temp ;
              }
              break ;
            default:
              cerr << "Bad type for relaxationFactor in momentumEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *vRelaxationFactor=0.5 ;
        }

        // Maximum number of iterations.
        if((*momentumEquationOptions).optionExists("maxIterations")){
          Loci::option_value_type optionValueType=momentumEquationOptions->
            getOptionValueType("maxIterations") ;
          switch(optionValueType){
            case Loci::REAL:
              {
                real temp ;
                momentumEquationOptions->getOption("maxIterations",temp) ;
                if(int(temp)<0){
                  cerr << "Bad maxIterations value for momentumEquation."
                    << endl ; Loci::Abort() ;
                }
                *vMaxIterations=int(temp) ;
              }
              break ;
            default:
              cerr << "Bad type for maxIterations in momentumEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          *vMaxIterations=5 ;
        }
      }
  } ;

  register_rule<MomentumEquationSolverParameters>
    registerMomentumEquationSolverParamters ;

//-----------------------------------------------------------------------------
// Velocity boundary condition rules.

  // Assigns boundary velocity from interpolated interface value.
  class BoundaryVelocityInterpolation : public pointwise_rule {
    private:
      const_store<vect3d> vI ;
      store<vect3d> v_f ;
    public:

      // Define input and output.
      BoundaryVelocityInterpolation() {
        name_store("interpolateFace_v3d(v)",vI) ;
        name_store("interpolated::v_f",v_f) ;
        input("interpolateFace_v3d(v)") ;
        output("interpolated::v_f") ;
        constraint("interface_BC") ;
      }

      // Calculate velocity for a single face.
      void calculate(Entity face) { v_f[face]=vI[face] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryVelocityInterpolation>
    registerBoundaryVelocityInterpolation ;

  // Rule for boundary faces with specified velocity. Assigns velocity value
  // to all boundary faces that have the property constantV_BC.
  class BoundaryVelocitySpecification : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<vect3d> constantV_BC ;
      store<vect3d> v_f ;
    public:

      // Define input and output.
      BoundaryVelocitySpecification() {
        name_store("ref",ref) ;
        name_store("constantV_BC",constantV_BC) ;
        name_store("v_f",v_f) ;
        input("ref->constantV_BC") ;
        output("v_f") ;
      }

      // Calculate face velocity for a single face.
      void calculate(Entity face) {
        v_f[face]=constantV_BC[ref[face]] ;
      }

      // Calculate face velocity for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundaryVelocitySpecification>
    registerBoundaryVelocitySpecification ;

  // Rule for specifying velocity with a profile. A single Cartesian coordinate
  // is used for the interpolation.
  class BoundaryVelocityProfileCartesian : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<vect3d> faceCenter ;
      const_store<string> cartesianV_BC ;
      store<vect3d> v_f ;
    public:

      // Define input and output.
      BoundaryVelocityProfileCartesian() {
        name_store("ref",ref) ;
        name_store("facecenter",faceCenter) ;
        name_store("cartesianV_BC",cartesianV_BC) ;
        name_store("v_f",v_f) ;
        input("facecenter,ref->cartesianV_BC") ;
        output("v_f") ;
        disable_threading() ;
      }

      // Calculate velocity for all faces in sequence.
      virtual void compute(const sequence &seq) {

        if(Loci::GLOBAL_OR((seq != EMPTY))){

          // Create a map to organize ref values for the faces.
          std::map<int, Loci::entitySet> bcmap ;
          for(sequence::const_iterator si=seq.begin();si!=seq.end();++si){
            bcmap[ref[*si]]+=*si ;
          }

          // Loop through the map. Each map entry has differnet input file.
          std::map<int,Loci::entitySet>::iterator bci ;
          for(bci=bcmap.begin();bci!=bcmap.end();++bci) {

            // Open the file containing the profile.
            string fileName=cartesianV_BC[bci->first] ;
            ifstream in(fileName.c_str(),ios::in) ;
            if(in.fail()) {
              cerr << "Open failed on " << fileName.c_str() << endl ;
              Loci::Abort() ;
            }

            // Skip spaces.
            while(!in.eof() && isspace(in.peek())) in.get() ;

            // Read in the number of points on the boundary and the variable
            // flag which indicates the coordinate direction to use in the
            // interpolation.
            int np,coordFlag ; in >> np >> coordFlag ;
            if(np<2){
              cerr << "Bad number of data points in bc_v.dat." << endl ;
              Loci::Abort() ;
            }
            if(coordFlag<0 || coordFlag>2){
              cerr << "Bad coordinate flag in bc_v.dat." << endl ;
              Loci::Abort() ;
            }

            // Read in the coordinates and velocity values.
            vect3d *center=new vect3d[np] ; vect3d *v =new vect3d[np] ;
            for(int i=0;i<np;++i){ in >> center[i] ; in >> v[i] ; }

            vect3d first = center[0] ; sequence s=sequence(bci->second) ;
            for(sequence::const_iterator si=s.begin();si!=s.end();++si) {
              int face=*si ;
              const vect3d &fcenter=faceCenter[face] ;
              int id1=0 ; real d1=0.0 ;
              switch(coordFlag){
                case 0: d1=fabs(first.x-fcenter.x) ; break ;
                case 1: d1=fabs(first.y-fcenter.y) ; break ;
                case 2: d1=fabs(first.z-fcenter.z) ; break ;
              }

              // Find the nearest point.
              for(int i=1;i<np;++i) {
                real dis=0.0 ;
                switch(coordFlag){
                  case 0: dis=fabs(fcenter.x-center[i].x) ; break ;
                  case 1: dis=fabs(fcenter.y-center[i].y) ; break ;
                  case 2: dis=fabs(fcenter.z-center[i].z) ; break ;
                }
                if(d1>dis){ d1=dis ; id1=i ; }
              }

              // Find the second nearest point
              int id2=0 ; real dd=0.0 ;
              switch(coordFlag){
                case 0: dd=center[id1].x-fcenter.x ; break ;
                case 1: dd=center[id1].y-fcenter.y ; break ;
                case 2: dd=center[id1].z-fcenter.z ; break ;
              }
              if(fabs(dd)>EPSILON)
                if(dd<0.0)
                  if(id1<np-1){
                    id2=id1+1 ;
                  }else{
                    v_f[face]=v[np-1] ; continue ;
                  }
                else
                  if(id1>0){
                    id2=id1-1 ;
                  }else{
                    v_f[face]=v[0] ; continue ;
                  }
              else {
                v_f[face]=v[id1] ; continue ;
              }
              real d2=0.0 ;
              switch(coordFlag){
                case 0: d2=fabs(center[id2].x-fcenter.x) ; break ;
                case 1: d2=fabs(center[id2].y-fcenter.y) ; break ;
                case 2: d2=fabs(center[id2].z-fcenter.z) ; break ;
              }
              v_f[face]=(d2*v[id1]+d1*v[id2])/(d1+d2) ;
            }
            delete [] center ; delete [] v ;
          }
        }
      }
  } ;

  register_rule<BoundaryVelocityProfileCartesian>
    registerBoundaryVelocityProfileCartesian ;

  // Rule for specifying an axisymmetric velocity from a radial profile.
  class BoundaryVelocityProfileAxisymmetric : public pointwise_rule {
    private:
      const_Map ref ;
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<string> axisymmetricV_BC ;
      const_store<unsigned int> referenceFrame_BC ;
      const_store<vect3d> faceCenter ;
      store<vect3d> v_f ;
    public:

      // Define input and output.
      BoundaryVelocityProfileAxisymmetric() {
        name_store("ref",ref) ;
        name_store("referenceFrame",referenceFrame) ;
        name_store("axisymmetricV_BC",axisymmetricV_BC) ;
        name_store("referenceFrame_BC",referenceFrame_BC) ;
        name_store("facecenter",faceCenter) ;
        name_store("v_f",v_f) ;
        input("referenceFrame,ref->(referenceFrame_BC,axisymmetricV_BC)") ;
        input("facecenter") ;
        output("v_f") ;
        disable_threading() ;
      }

      // Calculate velocity for all faces in sequence.
      virtual void compute(const sequence &seq) {

        if(Loci::GLOBAL_OR((seq != EMPTY))){

          // Create a map to organize ref values for the faces.
          std::map<int, Loci::entitySet> bcmap ;
          for(sequence::const_iterator si=seq.begin();si!=seq.end();++si){
            bcmap[ref[*si]]+=*si ;
          }

          // Loop through the map. Each map entry has differnet input file.
          std::map<int,Loci::entitySet>::iterator bci ;
          for(bci=bcmap.begin();bci!=bcmap.end();++bci) {

            // Open the file containing the profile.
            string fileName=axisymmetricV_BC[bci->first] ;
            ifstream in(fileName.c_str(),ios::in) ;
            if(in.fail()) {
              cerr << "Open failed on " << fileName.c_str() << endl ;
              Loci::Abort() ;
            }

            // Skip spaces.
            while(!in.eof() && isspace(in.peek())) in.get() ;

            // Read in the number of points in the profile.
            int np ; in >> np ;
            if(np<2){
              cerr << "Bad number of data points in bc_v.dat." << endl ;
              Loci::Abort() ;
            }

            // Read in the radius and velocity values. The order of the
            // velocity components is vA, vT and vR.
            real *r=new real[np] ; vect3d *v =new vect3d[np] ;
            for(int i=0;i<np;++i){ in >> r[i] ; in >> v[i] ; }

            real first=r[0] ; sequence s=sequence(bci->second) ;
            for(sequence::const_iterator si=s.begin();si!=s.end();++si) {

              // Compute the axial, radial and theta unit vectors.
              int face=*si ; const vect3d &fcenter=faceCenter[face] ;
              vect3d delta=(*referenceFrame)[referenceFrame_BC[bci->first]].
                axisEnd-(*referenceFrame)[referenceFrame_BC[bci->first]].
                axisStart ;
              vect3d pos=fcenter-(*referenceFrame)[referenceFrame_BC
                [bci->first]].axisStart ;
              vect3d iHatZ=(1.0/norm(delta))*delta ;
              vect3d iR=pos-dot(pos,iHatZ)*iHatZ ;
              vect3d iHatR=(1.0/norm(iR))*iR,iHatT=cross(iHatZ,iHatR) ;

              // Set up the transformation tensor.
              tens3d t ;
              t.x.x=iHatZ.x ; t.x.y=iHatT.x ; t.x.z=iHatR.x ;
              t.y.x=iHatZ.y ; t.y.y=iHatT.y ; t.y.z=iHatR.y ;
              t.z.x=iHatZ.z ; t.z.y=iHatT.z ; t.z.z=iHatR.z ;

              // Find the nearest point.
              int id1=0 ; real d1=fabs(first-norm(iR)) ;
              for(int i=1;i<np;++i) {
                real dis=fabs(norm(iR)-r[i]) ; if(d1>dis){ d1=dis ; id1=i ; }
              }

              // Find the second nearest point
              int id2 ; real dd=r[id1]-norm(iR) ;
              if(fabs(dd)>EPSILON)
                if(dd<0.0)
                  if(id1<np-1){
                    id2=id1+1 ;
                  }else{
                    v_f[face]=dotTemp(t,v[np-1]) ; continue ;
                  }
                else
                  if(id1>0){
                    id2=id1-1 ;
                  }else{
                    v_f[face]=dotTemp(t,v[0]) ; continue ;
                  }
              else {
                v_f[face]=dotTemp(t,v[id1]) ; continue ;
              }
              real d2=fabs(r[id2]-norm(iR)) ;
              v_f[face]=dotTemp(t,(d2*v[id1]+d1*v[id2])/(d1+d2)) ;
            }
            delete [] r ; delete [] v ;
          }
        }
      }
  } ;

  register_rule<BoundaryVelocityProfileAxisymmetric>
    registerBoundaryVelocityProfileAxisymmetric ;

  // Rule for specifying velocity as a function of x, y and z.
  class BoundaryVelocityFunctionSteady : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<vect3d> faceCenter ;
      const_store<string> functionVX_BC,functionVY_BC,functionVZ_BC ;
      store<vect3d> v_f ;
    public:

      // Define input and output.
      BoundaryVelocityFunctionSteady() {
        name_store("ref",ref) ;
        name_store("facecenter",faceCenter) ;
        name_store("functionVX_BC",functionVX_BC) ;
        name_store("functionVY_BC",functionVY_BC) ;
        name_store("functionVZ_BC",functionVZ_BC) ;
        name_store("v_f",v_f) ;
        input("facecenter,ref->(functionVX_BC,functionVY_BC,functionVZ_BC)") ;
        output("v_f") ;
        constraint("ref->v_functionSteady_BC") ;
      }

      // Calculate face velocity for a single face.
      void calculate(Entity face) {
        map<string,real> varMap ; varMap["x"]=faceCenter[face].x ;
        varMap["y"]=faceCenter[face].y ; varMap["z"]=faceCenter[face].z ;
        Loci::exprP pX=Loci::expression::create(functionVX_BC[ref[face]]) ;
        Loci::exprP pY=Loci::expression::create(functionVY_BC[ref[face]]) ;
        Loci::exprP pZ=Loci::expression::create(functionVZ_BC[ref[face]]) ;
        real vX=pX->evaluate(varMap),vY=pY->evaluate(varMap),
          vZ=pZ->evaluate(varMap) ;
        v_f[face]=vect3d(vX,vY,vZ) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryVelocityFunctionSteady>
    registerBoundaryVelocityFunctionSteady ;

  // Rule for specifying velocity as a function of x, y, z and time. Since
  // Loci will not promote part of a variable (part meaning a certain sub-
  // sequence of the total sequence) we must implement this in two parts.
  // This rule outputs a temporary variable which computes at {n} and
  // should promote to {n,it}. The following iteration independent rule
  // will use this temporary value to assign v_f at each iteration. So, the
  // computation is done at {n} as we want it.
  class BoundaryVelocityFunctionUnsteadyCompute : public pointwise_rule {
    private:
      const_Map ref ;
      const_param<real> sTime ;
      const_param<real> dt ;
      const_store<vect3d> faceCenter ;
      const_store<string> functionVX_BC,functionVY_BC,functionVZ_BC ;
      store<vect3d> tempV_f ;
    public:

      // Define input and output.
      BoundaryVelocityFunctionUnsteadyCompute() {
        name_store("ref",ref) ;
        name_store("stime{n}",sTime) ;
        name_store("dt{n}",dt) ;
        name_store("facecenter{n}",faceCenter) ;
        name_store("functionVX_BC{n}",functionVX_BC) ;
        name_store("functionVY_BC{n}",functionVY_BC) ;
        name_store("functionVZ_BC{n}",functionVZ_BC) ;
        name_store("tempV_f{n}",tempV_f) ;
        input("stime{n},dt{n},facecenter{n}") ;
        input("ref->(functionVX_BC{n},functionVY_BC{n},functionVZ_BC{n})") ;
        output("tempV_f{n}") ;
        constraint("ref->v_functionUnsteady_BC") ;
      }

      // Calculate face velocity for a single face.
      void calculate(Entity face) {
        map<string,real> varMap ; varMap["x"]=faceCenter[face].x ;
        varMap["y"]=faceCenter[face].y ; varMap["z"]=faceCenter[face].z ;
        varMap["t"]=(*sTime)+(*dt) ;
        Loci::exprP pX=Loci::expression::create(functionVX_BC[ref[face]]) ;
        Loci::exprP pY=Loci::expression::create(functionVY_BC[ref[face]]) ;
        Loci::exprP pZ=Loci::expression::create(functionVZ_BC[ref[face]]) ;
        real vX=pX->evaluate(varMap),vY=pY->evaluate(varMap),
          vZ=pZ->evaluate(varMap) ;
        tempV_f[face]=vect3d(vX,vY,vZ) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryVelocityFunctionUnsteadyCompute>
    registerBoundaryVelocityFunctionUnsteadyCompute ;

  // Assigns the boundary velocity at each iteration.
  class BoundaryVelocityFunctionUnsteadyAssign : public pointwise_rule {
    private:
      const_store<vect3d> tempV_f ;
      store<vect3d> v_f ;
    public:

      // Define input and output.
      BoundaryVelocityFunctionUnsteadyAssign() {
        name_store("tempV_f",tempV_f) ;
        name_store("v_f",v_f) ;
        input("tempV_f") ;
        output("v_f") ;
        constraint("ref->v_functionUnsteady_BC") ;
      }

      // Calculate face velocity for a single face.
      void calculate(Entity face) { v_f[face]=tempV_f[face] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryVelocityFunctionUnsteadyAssign>
    registerBoundaryVelocityFunctionUnsteadyAssign ;

  // Rule for for assigning velocity on inlet faces with specified mass flux.
  // Bug fix on 02/16/2005. Added negative to produce positive velocities from
  // negative mass flux value. Note that this rule used to be constrained to
  // subsonicInlet_BC, but since we now need this rule for multi-species
  // incompressible flow, this has been removed and the rule renamed. Made
  // this a priority rule so as to override the rule below when the flow
  // direction is specified.
  class BoundaryVelocityFromMassFlux : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<real> massFlux_BC ;
      const_store<vect3d> flowDirection_BC ;
      const_store<real> rho_f ;
      store<vect3d> v_f ;
    public:

      // Define input and output.
      BoundaryVelocityFromMassFlux() {
        name_store("ref",ref) ;
        name_store("massFlux_BC",massFlux_BC) ;
        name_store("flowDirection_BC",flowDirection_BC) ;
        name_store("rho_f",rho_f) ;
        name_store("specifiedFlowDirection::v_f",v_f) ;
        input("ref->massFlux_BC,ref->flowDirection_BC,rho_f") ;
        output("specifiedFlowDirection::v_f") ;
      }

      // Calculate face velocity for a single face.
      void calculate(Entity face) {
        v_f[face]=(-massFlux_BC[ref[face]]/rho_f[face])*
          flowDirection_BC[ref[face]] ;
      }

      // Calculate face velocity for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryVelocityFromMassFlux>
    registerBoundaryVelocityFromMassFlux ;

  // Rule for for assigning velocity on inlet faces with specified mass flux.
  // This rule assigns a velocity normal to the boundary face.
  class BoundaryVelocityFromMassFluxNormal : public pointwise_rule {
    private:
      const_Map ref ;
      const_store<real> massFlux_BC ;
      const_store<real> rho_f ;
      const_store<Area> area ;
      store<vect3d> v_f ;
    public:

      // Define input and output.
      BoundaryVelocityFromMassFluxNormal() {
        name_store("ref",ref) ;
        name_store("massFlux_BC",massFlux_BC) ;
        name_store("rho_f",rho_f) ;
        name_store("area",area) ;
        name_store("v_f",v_f) ;
        input("ref->massFlux_BC,rho_f,area") ;
        output("v_f") ;
      }

      // Calculate face velocity for a single face.
      void calculate(Entity face) {
        v_f[face]=(massFlux_BC[ref[face]]/rho_f[face])*area[face].n ;
      }

      // Calculate face velocity for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryVelocityFromMassFluxNormal>
    registerBoundaryVelocityFromMassFluxNormal ;

  // Rule for assigning velocity on no-slip boundary faces.
  class BoundaryVelocityNoSlip : public pointwise_rule {
    private:
      store<vect3d> v_f ;
    public:

      // Define input and output.
      BoundaryVelocityNoSlip() {
        name_store("v_f",v_f) ;
        output("v_f") ;
        constraint("noslip_BC") ;
      }

      // Calculate face velocity for a single face.
      void calculate(Entity face) { 
	v_f[face]=vect3d(0.0,0.0,0.0) ; 
	cout << "face_v=v_f=momentum" << v_f[face] << endl ;
      }

      // Calculate face velocity for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryVelocityNoSlip> registerBoundaryVelocityNoSlip ;

  // Priority rule for assigning velocity on no-slip boundary faces that are
  // rotating. Note that this rule is for non-turbomachinery calculations
  // where all we want is the ability to assign cross(omega,r) as the velocity
  // on a noslip boundary. One cannot use the 'referenceFrame' option and the
  // 'rotating' option at the same time.
  class BoundaryVelocityNoSlipReferenceFrame : public pointwise_rule {
    private:
      const_Map ref ;
      const_param<vector<ReferenceFrame> > referenceFrame ;
      const_store<unsigned int> referenceFrame_BC ;
      const_store<vect3d> faceCenter ;
      store<vect3d> v_f ;
    public:

      // Define input and output.
      BoundaryVelocityNoSlipReferenceFrame() {
        name_store("ref",ref) ;
        name_store("referenceFrame",referenceFrame) ;
        name_store("referenceFrame_BC",referenceFrame_BC) ;
        name_store("facecenter",faceCenter) ;
        name_store("rotating::v_f",v_f) ;
        input("referenceFrame,ref->referenceFrame_BC,facecenter") ;
        output("rotating::v_f") ;
        constraint("noslip_BC,ref->referenceFrame_BCoption") ;
      }

      // Calculate face velocity for a single face.
      void calculate(Entity face) {
        unsigned int n=ref[face] ;
        vect3d delta=(*referenceFrame)[referenceFrame_BC[n]].axisEnd-
        (*referenceFrame)[referenceFrame_BC[n]].axisStart,omega=
        ((*referenceFrame)[referenceFrame_BC[n]].omega/norm(delta))*delta,
          r=faceCenter[face]-(*referenceFrame)[referenceFrame_BC[n]].axisStart ;
        v_f[face]=cross(omega,r) ;
      }

      // Calculate face velocity for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryVelocityNoSlipReferenceFrame>
    registerBoundaryVelocityNoSlipReferenceFrame ;

  // Rule for for assigning velocity on slip boundary faces. We have replaced
  // the old slip bc that we never used for stability reasons with this new
  // one since we need to distinguish between slip and symmetry for moving
  // boundary problems (i.e. the normal component of the fluid velocity on a
  // slip surface is equal to the face velocity, but zero for a symmetry
  // surface).
  class BoundaryVelocitySlip : public pointwise_rule {
    private:
      const_store<Area> area ;
      const_Map ci ;
      const_store<vect3d> v ;
      store<vect3d> v_f ;
    public:

      // Define input and output.
      BoundaryVelocitySlip() {
        name_store("area",area) ;
        name_store("ci",ci) ;
        name_store("v",v) ;
        name_store("v_f",v_f) ;
        input("area,ci->v") ;
        output("v_f") ;
        constraint("slip_BC") ;
      }

      // Calculate face velocity for a single face.
      void calculate(Entity face) {
        const vect3d vCell=v[ci[face]] ;
        v_f[face]=vCell-dot(vCell,area[face].n)*area[face].n ;
      }

      // Calculate face velocity for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryVelocitySlip> registerBoundaryVelocitySlip ;

  // Rule for for assigning velocity on outlet faces. Right now we are using
  // a simple zero-gradient extrapolation. This rule does not allow
  // entrainment to occur.
  class BoundaryVelocityOutlet : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<Area> area ;
      const_store<vect3d> v ;
      store<vect3d> v_f ;
    public:

      // Define input and output.
      BoundaryVelocityOutlet() {
        name_store("ci",ci) ;
        name_store("area",area) ;
        name_store("v",v) ;
        name_store("v_f",v_f) ;
        input("ci->v,area") ;
        output("v_f") ;
        constraint("outlet_BC") ;
      }

      // Calculate face velocity for a single face.
      void calculate(Entity face) {
        v_f[face]=(dot(v[ci[face]],area[face].n)>=0.0)? v[ci[face]]:
          vect3d(0.0,0.0,0.0) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryVelocityOutlet> registerBoundaryVelocityOutlet ;

  // Rule for for assigning velocity on outlet faces. This priority rule
  // allows for entrainment.
  class BoundaryVelocityOutletEntrainment : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<vect3d> v ;
      store<vect3d> v_f ;
    public:

      // Define input and output.
      BoundaryVelocityOutletEntrainment() {
        name_store("ci",ci) ;
        name_store("v",v) ;
        name_store("entrainment::v_f",v_f) ;
        input("ci->v") ;
        output("entrainment::v_f") ;
        constraint("outlet_BC,ref->entrainment_BCoption") ;
      }

      // Calculate face velocity for a single face.
      void calculate(Entity face) { v_f[face]=v[ci[face]] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryVelocityOutletEntrainment>
    registerBoundaryVelocityOutletEntrainment ;

  // Unit rule for mass transfer calculation. A separate calculation is
  // required since we must use v[ci[face]] at the boundary.
  class ExtrapolatedPressureMassFlowRateUnit : public unit_rule {
    private:
      param<real> extrapolatedPressureMassFlowRate ;
    public:
                                                                                
      // Define input and output.
      ExtrapolatedPressureMassFlowRateUnit() {
        name_store("extrapolatedPressureMassFlowRate",
          extrapolatedPressureMassFlowRate) ;
        output("extrapolatedPressureMassFlowRate") ;
        constraint("extrapolatedPressureOutlet_BC") ;
      }
                                                                                
      // Initialize the value.
      void compute(const sequence &seq) {
        *extrapolatedPressureMassFlowRate=0.0 ;
      }
  } ;
                                                                                
  register_rule<ExtrapolatedPressureMassFlowRateUnit>
    registerExtrapolatedPressureMassFlowRateUnit ;

  // Apply rule for mass transfer calculation. A separate calculation is
  // required since we must use v[ci[face]] at the boundary.
  class ExtrapolatedPressureMassFlowRateApply : public
  apply_rule<param<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<vect3d> v ;
      const_store<real> rho_f ;
      const_store<Area> area ;
      param<real> extrapolatedPressureMassFlowRate ;
    public:
                                                                                
      // Define input and output.
      ExtrapolatedPressureMassFlowRateApply() {
        name_store("ci",ci) ;
        name_store("v",v) ;
        name_store("area",area) ;
        name_store("rho_f",rho_f) ;
        name_store("extrapolatedPressureMassFlowRate",
          extrapolatedPressureMassFlowRate) ;
        input("ci->v,area,rho_f") ;
        output("extrapolatedPressureMassFlowRate") ;
        constraint("extrapolatedPressureOutlet_BC") ;
      }
                                                                                
      // Add the face value to the total.
      void calculate(Entity face) { join(*extrapolatedPressureMassFlowRate,
        rho_f[face]*dot(v[ci[face]],area[face].n)*area[face].sada) ; }
                                                                                
      // Loop over the faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<ExtrapolatedPressureMassFlowRateApply>
    registerExtrapolatedPressureMassFlowRateApply ;

  // Priority rule for global mass conservation. IMPORTANT: This rule will not
  // be activated if there are both fixed pressure outlets and extrapolated
  // pressure outlets. In this case, the global conservation constraint is
  // not created in setupBC.cc .
  class BoundaryVelocityGlobalConservationOutlet : public pointwise_rule {
    private:
      const_Map ci ;
      const_param<real> totalInletMassFlowRate ;
      const_param<real> extrapolatedPressureMassFlowRate ;
      const_param<real> totalOutletArea ;
      const_store<Area> area ;
      const_store<vect3d> v ;
      const_store<real> rho_f ;
      store<vect3d> v_f ;
    private:
      real temp ;
    public:

      // Define input and output.
      BoundaryVelocityGlobalConservationOutlet() {
        name_store("ci",ci) ;
        name_store("massTransfer(incompressibleInlet_BC)",
          totalInletMassFlowRate) ;
        name_store("extrapolatedPressureMassFlowRate",
          extrapolatedPressureMassFlowRate) ;
        name_store("totalArea(extrapolatedPressureOutlet_BC)",
          totalOutletArea) ;
        name_store("area",area) ;
        name_store("rho_f",rho_f) ;
        name_store("v",v) ;
        name_store("globalConservation::entrainment::v_f",v_f) ;
        input("ci->v,area,rho_f") ;
        input("massTransfer(incompressibleInlet_BC)") ;
        input("extrapolatedPressureMassFlowRate") ;
        input("totalArea(extrapolatedPressureOutlet_BC)") ;
        output("globalConservation::entrainment::v_f") ;
        constraint("incompressibleFlow,extrapolatedPressureOutlet_BC") ;
        constraint("globalConservation") ;
      }

      // Calculate face velocity for a single face.
      void calculate(Entity face) {
        v_f[face]=v[ci[face]]-(temp/rho_f[face])*area[face].n ;
      }

      // Calculate face velocity for all faces in sequence.
      virtual void compute(const sequence &seq) {
        temp=(*totalInletMassFlowRate+*extrapolatedPressureMassFlowRate)/
          (*totalOutletArea) ;
        do_loop(seq,this) ;
      }
  } ;

  register_rule<BoundaryVelocityGlobalConservationOutlet>
    registerBoundaryVelocityGlobalConservationOutlet ;

  // Rule for for assigning velocity on symmetry boundary faces.
  class BoundaryVelocitySymmetry : public pointwise_rule {
    private:
      const_store<Area> area ;
      const_Map ci ;
      const_store<vect3d> v ;
      store<vect3d> v_f ;
    public:

      // Define input and output.
      BoundaryVelocitySymmetry() {
        name_store("area",area) ;
        name_store("ci",ci) ;
        name_store("v",v) ;
        name_store("v_f",v_f) ;
        input("area,ci->v") ;
        output("v_f") ;
        constraint("symmetry_BC") ;
      }

      // Calculate face velocity for a single face.
      void calculate(Entity face) {
        const vect3d vCell=v[ci[face]] ;
        v_f[face]=vCell-dot(vCell,area[face].n)*area[face].n ;
      }

      // Calculate face velocity for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryVelocitySymmetry> registerBoundaryVelocitySymmetry ;

  // Rule for for assigning velocity and pressure on total-pressure inlet
  // faces where the temperature is specified.
  class TotalPressureInletSpecifiedTemperature : public pointwise_rule {
    private:
      const_Map ci,ref ;
      const_param<EOS> eos ;
      const_store<real> p ;
      const_store<vect3d> flowDirection_BC ;
      const_store<real> p0_BC,T_BC ;
      const_storeVec<real> y_f ;
      store<vect3d> v_f ;
      store<real> p_f ;
    public:

      // Define input and output.
      TotalPressureInletSpecifiedTemperature() {
        name_store("ci",ci) ;
        name_store("ref",ref) ;
        name_store("eos",eos) ;
        name_store("p",p) ;
        name_store("flowDirection_BC",flowDirection_BC) ;
        name_store("p0_BC",p0_BC) ;
        name_store("T_BC",T_BC) ;
        name_store("y_f",y_f) ;
        name_store("v_f",v_f) ;
        name_store("totalPressureInlet::p_f",p_f) ;
        input("eos,ci->p,ref->(flowDirection_BC,p0_BC,T_BC),y_f") ;
        output("v_f,totalPressureInlet::p_f") ;
        constraint("compressibleFlow,totalPressureInlet_BC,ref->T_BC") ;
      }

      // Calculate face velocity for a single face.
      void calculate(Entity face) {

        // Extrapolate pressure to boundary and compute the boundary velocity.
        // Note that we are now inverting the EOS to get the consistent state
        // at the face rather than using cell values for thermodynamic
        // quantities.
        p_f[face]=p[ci[face]] ;
        EOS::State s=eos->State_from_mixture_p_T(y_f[face],p_f[face],
          T_BC[ref[face]]) ;
        const real gamma=s.cpt()/s.cvt() ;
        const real machNumber=sqrt((2.0/(gamma-1.0))*(pow(p0_BC[ref[face]]/
          p_f[face],(gamma-1.0)/gamma)-1.0)) ;
        const real vMag=machNumber*s.soundSpeed() ;
        v_f[face]=vMag*flowDirection_BC[ref[face]] ;
      }

      // Calculate face velocity for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TotalPressureInletSpecifiedTemperature>
    registerTotalPressureInletSpecifiedTemperature ;

  // Rule for for assigning velocity and pressure on total-pressure inlet
  // faces where the total temperature is specified. Note that we cannot use
  // the eos_state for the adjacent cell since this is obtained at the end of
  // the loop and we need this rule to occur to assemble the governing
  // equations.
  class TotalPressureInletSpecifiedTotalTemperature : public pointwise_rule {
    private:
      const_Map ci,ref ;
      const_param<EOS> eos ;
      const_store<real> p,temperature ;
      const_store<vect3d> flowDirection_BC ;
      const_store<real> p0_BC,T0_BC ;
      const_storeVec<real> y_f ;
      store<vect3d> v_f ;
      store<real> p_f,temperature_f ;
    public:

      // Define input and output.
      TotalPressureInletSpecifiedTotalTemperature() {
        name_store("ci",ci) ;
        name_store("ref",ref) ;
        name_store("eos",eos) ;
        name_store("p",p) ;
        name_store("temperature",temperature) ;
        name_store("flowDirection_BC",flowDirection_BC) ;
        name_store("p0_BC",p0_BC) ;
        name_store("T0_BC",T0_BC) ;
        name_store("y_f",y_f) ;
        name_store("v_f",v_f) ;
        name_store("totalPressureInlet::p_f",p_f) ;
        name_store("temperature_f",temperature_f) ;
        input("eos,ci->(p,temperature)") ;
        input("ref->(flowDirection_BC,p0_BC,T0_BC),y_f") ;
        output("v_f,totalPressureInlet::p_f,temperature_f") ;
        constraint("compressibleFlow,ref->p0_BC,ref->T0_BC") ;
      }

      // Calculate face velocity for a single face.
      void calculate(Entity face) {

        // Extrapolate pressure to boundary and compute the boundary velocity
        // and temperature. This should really be a Newton iteration in order
        // to get thermodynamic properties that are consistent with the
        // boundary temperature, but for now we will use the adjacent cell
        // temperature.
        p_f[face]=p[ci[face]] ;
        EOS::State s=eos->State_from_mixture_p_T(y_f[face],p_f[face],
          temperature[ci[face]]) ;
        const real gamma=s.cpt()/s.cvt() ;
        const real machNumber=sqrt((2.0/(gamma-1.0))*(pow(p0_BC[ref[face]]/
          p_f[face],(gamma-1.0)/gamma)-1.0)) ;
        const real vMag=machNumber*s.soundSpeed() ;
        v_f[face]=vMag*flowDirection_BC[ref[face]] ;
        temperature_f[face]=T0_BC[ref[face]]*pow(p_f[face]/p0_BC[ref[face]],
          (gamma-1.0)/gamma) ;
      }

      // Calculate face velocity for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TotalPressureInletSpecifiedTotalTemperature>
    registerTotalPressureInletSpecifiedTotalTemperature ;

//-----------------------------------------------------------------------------
// Testing: Rules to compute net outward mass flux from each cell. Used to
// prevent pathological conditions in inviscid flow

  // Rule to initialize the net mass flux.
  class InitializeNetMassFlux : public unit_rule {
    private:
      store<real> netMassFlux ;
    public:

      // Define input and output.
      InitializeNetMassFlux() {
        name_store("netMassFlux",netMassFlux) ;
        output("netMassFlux") ;
        constraint("vol") ;
      }

      // Set the main coefficient to zero for a single cell.
      void calculate(Entity cell) { netMassFlux[cell]=0.0 ; }

      // Set the main coefficient to zero for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeNetMassFlux> registerInitializeNetMassFlux ;

  // Rule to add interior face contributions to net mass flux.
  class NetMassFluxInterior : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_store<real> massFlux ;
      store<real> netMassFlux ;
    public:

      // Define input and output.
      NetMassFluxInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("massFlux",massFlux) ;
        name_store("netMassFlux",netMassFlux) ;
        input("massFlux") ;
        output("(cl,cr)->netMassFlux") ;
        constraint("internalFaces") ;
      }

      // Increment the net mass flux for the cells attach to a single face.
      void calculate(Entity face) {
        netMassFlux[cl[face]]+=massFlux[face] ;
        netMassFlux[cr[face]]-=massFlux[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<NetMassFluxInterior> registerNetMassFluxInterior ;

  // Rule to add boundary face contributions to net mass flux.
  class NetMassFluxBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_store<real> massFlux ;
      store<real> netMassFlux ;
    public:

      // Define input and output.
      NetMassFluxBoundary() {
        name_store("ci",ci) ;
        name_store("massFlux",massFlux) ;
        name_store("netMassFlux",netMassFlux) ;
        input("massFlux") ;
        output("ci->netMassFlux") ;
        constraint("boundaryFaces") ;
      }

      // Increment the net mass flux for the cells attach to a single face.
      void calculate(Entity face) {
        netMassFlux[ci[face]]+=massFlux[face] ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<NetMassFluxBoundary> registerNetMassFluxBoundary ;

//-----------------------------------------------------------------------------
// Rules for computing eddy viscosity for laminar flows.

  // Rule to assign eddy viscosity for cells.
  class EddyViscosityLaminarFlowInterior : public pointwise_rule {
    private:
      store<real> eddyViscosity ;
    public:

      // Define input and output.
      EddyViscosityLaminarFlowInterior() {
        name_store("eddyViscosity",eddyViscosity) ;
        output("eddyViscosity") ;
        constraint("geom_cells,laminarFlow") ;
      }

      // Set eddy viscosity to zero for a single cell.
      void calculate(Entity cell) { eddyViscosity[cell]=0.0 ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<EddyViscosityLaminarFlowInterior>
    registerEddyViscosityLaminarFlowInterior ;

  // Rule to assign eddy viscosity for boundary faces.
  class EddyViscosityLaminarFlowBoundary : public pointwise_rule {
    private:
      store<real> eddyViscosity_f ;
    public:

      // Define input and output.
      EddyViscosityLaminarFlowBoundary() {
        name_store("eddyViscosity_f",eddyViscosity_f) ;
        output("eddyViscosity_f") ;
        constraint("boundaryFaces,laminarFlow") ;
      }

      // Set eddy viscosity to zero for a single face.
      void calculate(Entity face) { eddyViscosity_f[face]=0.0 ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<EddyViscosityLaminarFlowBoundary>
    registerEddyViscosityLaminarFlowBoundary ;

//-----------------------------------------------------------------------------
// Rules to create a constraint for boundary faces with non-zero velocity
// diffusion flux.

  // All inlet faces have non-zero diffusion flux.
  class BoundaryVelocityDiffusionInlet : public pointwise_rule {
    private:
      store<bool> boundaryVelocityDiffusion ;
    public:

      // Define input and output.
      BoundaryVelocityDiffusionInlet() {
        name_store("boundaryVelocityDiffusion",boundaryVelocityDiffusion) ;
        output("boundaryVelocityDiffusion") ;
        constraint("inlet_BC") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { boundaryVelocityDiffusion[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryVelocityDiffusionInlet>
    registerBoundaryVelocityDiffusionInlet ;

  // All outlet faces have non-zero diffusion flux.
  class BoundaryVelocityDiffusionOutlet : public pointwise_rule {
    private:
      store<bool> boundaryVelocityDiffusion ;
    public:

      // Define input and output.
      BoundaryVelocityDiffusionOutlet() {
        name_store("boundaryVelocityDiffusion",boundaryVelocityDiffusion) ;
        output("boundaryVelocityDiffusion") ;
        constraint("outlet_BC") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { boundaryVelocityDiffusion[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryVelocityDiffusionOutlet>
    registerBoundaryVelocityDiffusionOutlet ;

  // No-slip faces have non-zero diffusion flux.
  class BoundaryVelocityDiffusionNoSlip : public pointwise_rule {
    private:
      store<bool> boundaryVelocityDiffusion ;
    public:

      // Define input and output.
      BoundaryVelocityDiffusionNoSlip() {
        name_store("boundaryVelocityDiffusion",boundaryVelocityDiffusion) ;
        output("boundaryVelocityDiffusion") ;
        constraint("noslip_BC") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { boundaryVelocityDiffusion[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryVelocityDiffusionNoSlip>
    registerBoundaryVelocityDiffusionNoSlip ;

  // Interface faces have non-zero diffusion flux.
  class BoundaryVelocityDiffusionInterface : public pointwise_rule {
    private:
      store<bool> boundaryVelocityDiffusion ;
    public:

      // Define input and output.
      BoundaryVelocityDiffusionInterface() {
        name_store("boundaryVelocityDiffusion",boundaryVelocityDiffusion) ;
        output("boundaryVelocityDiffusion") ;
        constraint("interface_BC") ;
      }

      // Assign flag for a single boundary face.
      void calculate(Entity face) { boundaryVelocityDiffusion[face]=true ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryVelocityDiffusionInterface>
    registerBoundaryVelocityDiffusionInterface ;

//-----------------------------------------------------------------------------
// Scheme independent rules for assembling the momentum equation.

  // Rule to initialize the main coefficient. Checked.
  class InitializeVelocityMainCoefficient : public unit_rule {
    private:
      store<real> vMainCoefficient ;
    public:

      // Define input and output.
      InitializeVelocityMainCoefficient() {
        name_store("vMainCoefficient",vMainCoefficient) ;
        output("vMainCoefficient") ;
        constraint("vol") ;
      }

      // Set the main coefficient to zero for a single cell.
      void calculate(Entity cell) { vMainCoefficient[cell]=0.0 ; }

      // Set the main coefficient to zero for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeVelocityMainCoefficient>
    registerInitializeVelocityMainCoefficient ;

  // Rule to add the first-order inviscid flux contribution to the main
  // coefficient for interior faces. Checked.
  class FOUInviscidFluxToVelocityMainCoefficientInterior : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_param<real> thetaParameter ;
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      store<real> vMainCoefficient ;
    public:

      // Define input and output.
      FOUInviscidFluxToVelocityMainCoefficientInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("massFlux",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        input("thetaParameter,massFlux,gridMassFlux") ;
        output("(cl,cr)->vMainCoefficient") ;
        constraint("internalFaces") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        if(netMassFlux>0.0){
          vMainCoefficient[cr[face]]+=(*thetaParameter)*netMassFlux ;
        }else{
          vMainCoefficient[cl[face]]-=(*thetaParameter)*netMassFlux ;
        }
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToVelocityMainCoefficientInterior>
    registerFOUInviscidFluxToVelocityMainCoefficientInterior ;

  // Rule to add the first-order inviscid flux contribution to the main
  // coefficient for boundary faces.
  class FOUInviscidFluxToVelocityMainCoefficientBoundary : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_param<real> thetaParameter ;
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      store<real> vMainCoefficient ;
    public:

      // Define input and output.
      FOUInviscidFluxToVelocityMainCoefficientBoundary() {
        name_store("ci",ci) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("massFlux",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        input("thetaParameter,massFlux,gridMassFlux") ;
        output("ci->vMainCoefficient") ;
        constraint("boundaryFaces") ;
      }

      // Increment the main coefficient for the cell attach to a single face.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        if(netMassFlux<0.0) vMainCoefficient[ci[face]]-=(*thetaParameter)*
          netMassFlux ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToVelocityMainCoefficientBoundary>
    registerFOUInviscidFluxToVelocityMainCoefficientBoundary ;

  // Rule to add the diffusive flux contribution to the main coefficient for
  // interior faces. Right now we are assuming that we will be dealing with
  // turbulence models that interface to the flow equations via a turbulent
  // eddy viscosity.
  class DiffusiveFluxToVelocityMainCoefficientInterior : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_param<real> thetaParameter ;
      const_store<real> faceRadius ;
      const_store<real> diffusionProduct ;
      const_store<real> viscosity ;
      store<real> vMainCoefficient ;
    public:

      // Define input and output.
      DiffusiveFluxToVelocityMainCoefficientInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("faceRadius",faceRadius) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("viscosity",viscosity) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        input("thetaParameter,faceRadius,diffusionProduct,viscosity") ;
        output("(cl,cr)->vMainCoefficient") ;
        constraint("internalFaces,viscousFlow") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real temp=viscosity[face]*diffusionProduct[face]*(*thetaParameter)*
          faceRadius[face] ;
        vMainCoefficient[cl[face]]+=temp ; vMainCoefficient[cr[face]]+=temp ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToVelocityMainCoefficientInterior>
    registerDiffusiveFluxToVelocityMainCoefficientInterior ;

  // Rule to add the diffusive flux contribution to the velocity main
  // coefficient for boundary faces. Right now we are assuming that we will be
  // dealing with turbulence models that interface to the flow equations via a
  // turbulent eddy viscosity.
  class DiffusiveFluxToVelocityMainCoefficientBoundary : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_param<real> thetaParameter ;
      const_store<real> faceRadius ;
      const_store<real> diffusionProduct ;
      const_store<real> viscosity ;
      const_store<real> noWallFunction ;
      store<real> vMainCoefficient ;
    public:

      // Define input and output.
      DiffusiveFluxToVelocityMainCoefficientBoundary() {
        name_store("ci",ci) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("faceRadius",faceRadius) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("viscosity",viscosity) ;
        name_store("noWallFunction",noWallFunction) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        input("thetaParameter,faceRadius,noWallFunction") ;
        input("diffusionProduct,viscosity") ;
        output("ci->vMainCoefficient") ;
        constraint("boundaryVelocityDiffusion,viscousFlow") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        vMainCoefficient[ci[face]]+=viscosity[face]*diffusionProduct[face]*
         (*thetaParameter)*faceRadius[face]*noWallFunction[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToVelocityMainCoefficientBoundary>
    registerDiffusiveFluxToVelocityMainCoefficientBoundary ;

  // Rule to ensure main coefficient is not zero for inviscid flows.
  class NetMassFluxToVelocityMainCoefficient : public
  apply_rule<store<real>,Loci::Summation<real> > {
    private:
      const_store<real> netMassFlux ;
      store<real> vMainCoefficient ;
    public:

      // Define input and output.
      NetMassFluxToVelocityMainCoefficient() {
        name_store("netMassFlux",netMassFlux) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        input("netMassFlux") ;
        output("vMainCoefficient") ;
        constraint("inviscidFlow,geom_cells") ;
      }

      // Add net mass flux for a single cell.
      void calculate(Entity cell) {
        if(netMassFlux[cell]<0.0) vMainCoefficient[cell]-=netMassFlux[cell] ;
      }

      // Add net mass flux for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<NetMassFluxToVelocityMainCoefficient>
    registerNetMassFluxToVelocityMainCoefficient ;

  // Rule to initialize the source term. Checked.
  class InitializeVelocitySourceTerm : public unit_rule {
    private:
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      InitializeVelocitySourceTerm() {
        name_store("vSourceTerm",vSourceTerm) ;
        output("vSourceTerm") ;
        constraint("vol") ;
      }

      // Set the source term to zero for a single cell.
      void calculate(Entity cell) {
        vSourceTerm[cell]=vector3d<real>(0.0,0.0,0.0) ;
      }

      // Set the source term to zero for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeVelocitySourceTerm>
    registerInitializeVelocitySourceTerm ;

  // Rule to add the first-order inviscid flux contribution to the source term
  // for boundary faces. Checked.
  class FOUInviscidFluxToVelocitySourceTermBoundary : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map ci ;
      const_param<real> thetaParameter ;
      const_store<vect3d> v ;
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      const_store<vect3d> v_f ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      FOUInviscidFluxToVelocitySourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("v",v) ;
        name_store("massFlux",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("v_f",v_f) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("thetaParameter,ci->v,massFlux,gridMassFlux,v_f") ;
        output("ci->vSourceTerm") ;
        constraint("boundaryFaces,fouInviscidFlux") ;
      }

      // Increment the source term for the cell attach to a single face.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        if(netMassFlux<0.0) vSourceTerm[ci[face]]-=netMassFlux*v_f[face]*
         (*thetaParameter) ;
        else vSourceTerm[ci[face]]-=netMassFlux*(v_f[face]-v[ci[face]])*
          (*thetaParameter) ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToVelocitySourceTermBoundary>
    registerFOUInviscidFluxToVelocitySourceTermBoundary ;

  // Rule to add the second-order convection contribution to the source term
  // for interior faces.
  class SOUInviscidFluxToVelocitySourceTermInterior : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
      const_param<real> thetaParameter ;
      const_store<tens3d> vGradient ;
      const_store<vect3d> vLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      SOUInviscidFluxToVelocitySourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("limiterv3d(v)",vLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("massFlux",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("thetaParameter") ;
        input("(cl,cr)->(gradv3d(v),limiterv3d(v),cellcenter)") ;
        input("facecenter,massFlux,gridMassFlux") ;
        output("(cl,cr)->vSourceTerm") ;
        constraint("internalFaces,souOrRoeInviscidFlux") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        vect3d secondOrderSource=(netMassFlux>0.0)? netMassFlux*
          ComponentProduct(vLimiter[cl[face]],dotTemp(vGradient[cl[face]],
          faceCenter[face]-cellCenter[cl[face]]))*(*thetaParameter):
          netMassFlux*ComponentProduct(vLimiter[cr[face]],
          dotTemp(vGradient[cr[face]],faceCenter[face]-cellCenter[cr[face]]))*
          (*thetaParameter) ;
//if(face>=1320 && face<=1519){
//  cout << "face,lGrad: " << face << " " << vGradient[cl[face]] << endl ;
//  cout << "face,rGrad: " << face << " " << vGradient[cr[face]] << endl ;
//}
        vSourceTerm[cl[face]]-=secondOrderSource ;
        vSourceTerm[cr[face]]+=secondOrderSource ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SOUInviscidFluxToVelocitySourceTermInterior>
    registerSOUInviscidFluxToVelocitySourceTermInterior ;

  // Rule to add the second-order convection contribution to the source term
  // for boundary faces.
  class SOUInviscidFluxToVelocitySourceTermBoundary : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map ci ;
      const_param<real> thetaParameter ;
      const_store<tens3d> vGradient ;
      const_store<vect3d> vLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      const_store<vect3d> v_f ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      SOUInviscidFluxToVelocitySourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("limiterv3d(v)",vLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("massFlux",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("v_f",v_f) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("thetaParameter,ci->(gradv3d(v),limiterv3d(v),cellcenter)") ;
        input("facecenter,massFlux,gridMassFlux,v_f") ;
        output("ci->vSourceTerm") ;
        constraint("boundaryFaces,souOrRoeInviscidFlux") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        vSourceTerm[ci[face]]-=(netMassFlux>0.0)? netMassFlux*
          ComponentProduct(vLimiter[ci[face]],dotTemp(vGradient[ci[face]],
          faceCenter[face]-cellCenter[ci[face]]))*(*thetaParameter):
          netMassFlux*v_f[face]*(*thetaParameter) ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<SOUInviscidFluxToVelocitySourceTermBoundary>
    registerSOUInviscidFluxToVelocitySourceTermBoundary ;

  // Rule to add the defect correction term to get the central-difference
  // scheme for interior faces.
  class CDDefectCorrectionToVelocitySourceTermInterior : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
      const_store<vect3d> v ;
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      CDDefectCorrectionToVelocitySourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("v",v) ;
        name_store("massFlux",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("(cl,cr)->v,massFlux,gridMassFlux") ;
        output("(cl,cr)->vSourceTerm") ;
        constraint("internalFaces,cdInviscidFlux") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        vect3d source=(netMassFlux>0.0)? 0.5*netMassFlux*(v[cr[face]]-
          v[cl[face]]):0.5*netMassFlux*(v[cl[face]]-v[cr[face]]) ;
        vSourceTerm[cl[face]]-=source ; vSourceTerm[cr[face]]+=source ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CDDefectCorrectionToVelocitySourceTermInterior>
    registerCDDefectCorrectionToVelocitySourceTermInterior ;

  // Rule to add the pressure contribution to the source term for interior
  // faces.
  class PressureToVelocitySourceTermInterior : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
      const_param<real> thetaParameter ;
      const_store<real> faceRadius ;
      const_store<real> p ;
      const_store<vect3d> pGradient ;
      const_store<real> pLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      const_store<Area> area ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      PressureToVelocitySourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("faceRadius",faceRadius) ;
        name_store("p",p) ;
        name_store("grads(p)",pGradient) ;
        name_store("limiters(p)",pLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("area",area) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("thetaParameter,faceRadius") ;
        input("(cl,cr)->(p,grads(p),limiters(p),cellcenter)") ;
        input("facecenter,area") ;
        output("(cl,cr)->vSourceTerm") ;
        constraint("internalFaces") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        vect3d pressureSource=0.5*((p[cl[face]]+pLimiter[cl[face]]*
          dot(pGradient[cl[face]],(faceCenter[face]-cellCenter[cl[face]])))+
          (p[cr[face]]+pLimiter[cr[face]]*dot(pGradient[cr[face]],
          (faceCenter[face]-cellCenter[cr[face]]))))*area[face].n*
          area[face].sada*(*thetaParameter)*faceRadius[face] ;
        vSourceTerm[cl[face]]-=pressureSource ;
        vSourceTerm[cr[face]]+=pressureSource ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PressureToVelocitySourceTermInterior>
    registerPressureToVelocitySourceTermInterior ;

  // New rule to compute boundary pressure for momentum equation. Default is
  // high-order extrapolation to boundary.
  class BoundaryPressureDefault : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> p ;
      const_store<vect3d> pGradient ;
      const_store<real> pLimiter ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> faceCenter ;
      store<real> boundaryPressure ;
    public:

      // Define input and output.
      BoundaryPressureDefault() {
        name_store("ci",ci) ;
        name_store("p",p) ;
        name_store("grads(p)",pGradient) ;
        name_store("limiters(p)",pLimiter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("facecenter",faceCenter) ;
        name_store("boundaryPressure",boundaryPressure) ;
        input("ci->(p,grads(p),limiters(p),cellcenter),facecenter") ;
        output("boundaryPressure") ;
      }

      // Calculate pressure for a single face.
      void calculate(Entity face) {
        boundaryPressure[face]=p[ci[face]]+pLimiter[ci[face]]*
          dot(pGradient[ci[face]],(faceCenter[face]-cellCenter[ci[face]])) ;
      }

      // Calculate pressure for all faces in sequence.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryPressureDefault>
    registerBoundaryPressureDefault ;

  // New rule to compute boundary pressure for momentum equation. Priority is
  // for boundaries with specified pressure.
  class BoundaryPressurePriority : public pointwise_rule {
    private:
      const_store<real> p_f ;
      store<real> boundaryPressure ;
    public:

      // Define input and output.
      BoundaryPressurePriority() {
        name_store("p_f",p_f) ;
        name_store("specified::boundaryPressure",boundaryPressure) ;
        input("p_f") ;
        output("specified::boundaryPressure") ;
        constraint("specifiedPressure_BC") ;
      }

      // Calculate pressure for a single face.
      void calculate(Entity face) { boundaryPressure[face]=p_f[face] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BoundaryPressurePriority> registerBoundaryPressurePriority ;

  // Rule to add the pressure contribution to the source term for boundary
  // faces. Checked.
  class PressureToVelocitySourceTermBoundary : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_Map ci ;
      const_param<real> thetaParameter ;
      const_store<real> faceRadius ;
      const_store<real> boundaryPressure ;
      const_store<Area> area ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      PressureToVelocitySourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("faceRadius",faceRadius) ;
        name_store("boundaryPressure",boundaryPressure) ;
        name_store("area",area) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("thetaParameter,faceRadius,boundaryPressure,area") ;
        output("ci->vSourceTerm") ;
        constraint("boundaryFaces") ;
      }

      // Increment the source term for the cell attach to a single face.
      void calculate(Entity face) {
        vect3d pressureSource=boundaryPressure[face]*area[face].n*
          area[face].sada*(*thetaParameter)*faceRadius[face] ;
        vSourceTerm[ci[face]]-=pressureSource ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<PressureToVelocitySourceTermBoundary>
    registerPressureToVelocitySourceTermBoundary ;

  // Rule to add the Roe dissipation contribution to the source term for
  // interior faces.
  class RoeDissipationToVelocitySourceTermInterior : public apply_rule
  <store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
      const_store<real> delP ;
      const_store<Area> area ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      RoeDissipationToVelocitySourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("delP",delP) ;
        name_store("area",area) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("delP,area") ;
        output("(cl,cr)->vSourceTerm") ;
        constraint("internalFaces,roeInviscidFlux") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        vect3d source=0.5*delP[face]*area[face].n*area[face].sada ;
        vSourceTerm[cl[face]]+=source ; vSourceTerm[cr[face]]-=source ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<RoeDissipationToVelocitySourceTermInterior>
    registerRoeDissipationToVelocitySourceTermInterior ;

  // Rule to add the buoyancy contribution to the source term.
  class BuoyancyToVelocitySourceTerm : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_param<vect3d> gravity ;
      const_param<real> rhoReference ;
      const_store<real> rho ;
      const_store<real> vol ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      BuoyancyToVelocitySourceTerm() {
        name_store("gravityAcceleration",gravity) ;
        name_store("gravityRhoRef",rhoReference) ;
        name_store("rho",rho) ;
        name_store("vol",vol) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("gravityAcceleration,gravityRhoRef,rho,vol") ;
        output("vSourceTerm") ;
        constraint("geom_cells") ;
      }

      // Add buoyance force for single cell.
      void calculate(Entity cell) {
        vSourceTerm[cell]+=(rho[cell]-*rhoReference)*(*gravity)*vol[cell] ;
      }

      // Call calculate for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<BuoyancyToVelocitySourceTerm>
    registerBuoyancyToVelocitySourceTerm ;

  // Rule to add the diffusive flux contribution to the velocity source term for
  // interior faces. Right now we are assuming that we will be dealing with
  // turbulence models that interface to the flow equations via a turbulent
  // eddy viscosity.
  class DiffusiveFluxToVelocitySourceTermInterior : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
      const_param<real> thetaParameter ; 
      const_store<real> vol ;
      const_store<tens3d> vGradient ;
      const_store<vect3d> geometryFactor0 ;
      const_store<real> faceRadius ; 
      const_store<real> viscosity ; 
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToVelocitySourceTermInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("vol",vol) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("geometryFactor0",geometryFactor0) ;
        name_store("faceRadius",faceRadius) ;
        name_store("viscosity",viscosity) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("thetaParameter,(cl,cr)->(vol,gradv3d(v))") ;
        input("geometryFactor0,faceRadius,viscosity") ;
        output("(cl,cr)->vSourceTerm") ;
        constraint("internalFaces,viscousFlow") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        tens3d tempTensor=product(vol[cr[face]],vGradient[cl[face]])+
          product(vol[cl[face]],vGradient[cr[face]]) ;
        vect3d secondarySourceTerm=viscosity[face]*dotTemp(tempTensor,
          geometryFactor0[face])*(*thetaParameter)*faceRadius[face]/
          (vol[cl[face]]+vol[cr[face]]) ;
        vSourceTerm[cl[face]]+=secondarySourceTerm ;
        vSourceTerm[cr[face]]-=secondarySourceTerm ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToVelocitySourceTermInterior>
    registerDiffusiveFluxToVelocitySourceTermInterior ;

  // Rule to add the diffusive flux contribution to the velocity source term for
  // boundary faces. Right now we are assuming that we will be dealing with
  // turbulence models that interface to the flow equations via a turbulent
  // eddy viscosity.
  class DiffusiveFluxToVelocitySourceTermBoundary : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map ci ;
      const_param<real> thetaParameter ;
      const_store<vect3d> cellCenter ;
      const_store<tens3d> vGradient ;
      const_store<vect3d> faceCenter ;
      const_store<real> diffusionProduct ;
      const_store<vect3d> v_f ;
      const_store<real> viscosity ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      const_store<real> noWallFunction ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToVelocitySourceTermBoundary() {
        name_store("ci",ci) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("facecenter",faceCenter) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("v_f",v_f) ;
        name_store("viscosity",viscosity) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("noWallFunction",noWallFunction) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("thetaParameter,ci->cellcenter,ci->gradv3d(v)") ;
        input("facecenter,diffusionProduct,v_f,noWallFunction") ;
        input("viscosity,area,faceRadius") ;
        output("ci->vSourceTerm") ;
        constraint("boundaryVelocityDiffusion,viscousFlow") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        vect3d sourceTerm=viscosity[face]*(v_f[face]*diffusionProduct[face]+
          dotTemp(vGradient[ci[face]], (area[face].n*area[face].sada-
          diffusionProduct[face]*(faceCenter[face]-cellCenter[ci[face]]))))*
          (*thetaParameter)*faceRadius[face]*noWallFunction[face] ;
        vSourceTerm[ci[face]]+=sourceTerm ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToVelocitySourceTermBoundary>
    registerDiffusiveFluxToVelocitySourceTermBoundary ;

  // Rule to add components of viscous flux that are non-zero only for
  // compressible flows.
  class DiffusiveFluxToVelocitySourceTermCompressibleInterior : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
      const_param<real> thetaParameter ; 
      const_store<real> vol ;
      const_store<tens3d> vGradient ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      const_store<real> viscosity ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToVelocitySourceTermCompressibleInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("vol",vol) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("viscosity",viscosity) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("thetaParameter,faceRadius,(cl,cr)->(vol,gradv3d(v))") ;
        input("area,viscosity") ;
        output("(cl,cr)->vSourceTerm") ;
        constraint("internalFaces,viscousFlow,compressibleFlow") ;
      }

      // Increment the main coefficient for the cells attach to a single face.
      void calculate(Entity face) {
        real vR=vol[cr[face]]/(vol[cl[face]]+vol[cr[face]]) ;
        real vL=vol[cl[face]]/(vol[cl[face]]+vol[cr[face]]) ;
        tens3d tempTensor=Transpose(product(vR,vGradient[cl[face]])+product(vL,
          vGradient[cr[face]])) ;
        real temp=2.0*Trace(tempTensor)/3.0 ;
        tempTensor.x.x-=temp ; tempTensor.y.y-=temp ; tempTensor.z.z-=temp ;
        vect3d secondarySourceTerm=viscosity[face]*dotTemp(tempTensor,
          area[face].n*area[face].sada)*(*thetaParameter)*faceRadius[face] ;
        vSourceTerm[cl[face]]+=secondarySourceTerm ;
        vSourceTerm[cr[face]]-=secondarySourceTerm ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToVelocitySourceTermCompressibleInterior>
    registerDiffusiveFluxToVelocitySourceTermCompressibleInterior ;

  // Rule to add components of viscous flux that are non-zero only for
  // compressible flows.
  class DiffusiveFluxToVelocitySourceTermCompressibleBoundary : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map ci ;
      const_param<real> thetaParameter ;
      const_store<tens3d> vGradient ;
      const_store<real> viscosity ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      const_store<real> noWallFunction ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToVelocitySourceTermCompressibleBoundary() {
        name_store("ci",ci) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("gradv3d(v)",vGradient) ;
        name_store("viscosity",viscosity) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("noWallFunction",noWallFunction) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("thetaParameter,ci->gradv3d(v),noWallFunction") ;
        input("viscosity,area,faceRadius") ;
        output("ci->vSourceTerm") ;
        constraint("boundaryVelocityDiffusion,viscousFlow") ;
        constraint("compressibleFlow") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        tens3d tempTensor=Transpose(vGradient[ci[face]]) ;
        real temp=2.0*Trace(tempTensor)/3.0 ;
        tempTensor.x.x-=temp ; tempTensor.y.y-=temp ; tempTensor.z.z-=temp ;
        vect3d secondarySourceTerm=viscosity[face]*dotTemp(tempTensor,
          area[face].n*area[face].sada)*(*thetaParameter)*faceRadius[face]*
          noWallFunction[face] ;
        vSourceTerm[ci[face]]+=secondarySourceTerm ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToVelocitySourceTermCompressibleBoundary>
    registerDiffusiveFluxToVelocitySourceTermCompressibleBoundary ;

  // Rule to add the wall shear stress determined from wall functions to the
  // source term for boundary faces. We only need yPlusWall here to trigger
  // output of yPlus.
  class DiffusiveFluxToVelocitySourceTermBoundaryWallFunction : public
  apply_rule<store<vect3d>,Loci::Summation<vect3d> > {
    private:
      const_Map ci ;
      const_param<real> thetaParameter ;
      const_store<vect3d> tauWall ;
      const_store<real> yPlusWall ;
      const_store<Area> area ;
      const_store<real> faceRadius ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      DiffusiveFluxToVelocitySourceTermBoundaryWallFunction() {
        name_store("ci",ci) ;
        name_store("thetaParameter",thetaParameter) ;
        name_store("tauWall",tauWall) ;
        name_store("yPlusWall",yPlusWall) ;
        name_store("area",area) ;
        name_store("faceRadius",faceRadius) ;
        name_store("vSourceTerm",vSourceTerm) ;
        input("thetaParameter,tauWall,yPlusWall,area,faceRadius") ;
        output("ci->vSourceTerm") ;
        constraint("viscousFlow,ref->wallFunction_BCoption") ;
      }

      // Increment the source term for the cells attach to a single face.
      void calculate(Entity face) {
        vSourceTerm[ci[face]]+=tauWall[face]*area[face].sada*(*thetaParameter)*
          faceRadius[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToVelocitySourceTermBoundaryWallFunction>
    registerDiffusiveFluxToVelocitySourceTermBoundaryWallFunction ;

  // Rule to compute the diagonal term for the linear system. Checked.
  class ComputeVelocityMatrixDiagonal : public pointwise_rule {
    private:
      const_param<real> vRelaxationFactor ;
      const_store<real> vMainCoefficient ;
      store<real> D ;
    public:

      // Define input and output.
      ComputeVelocityMatrixDiagonal() {
        name_store("vRelaxationFactor",vRelaxationFactor) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        name_store("vStar_D",D) ;
        input("vRelaxationFactor,vMainCoefficient") ;
        output("vStar_D") ;
        constraint("geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) {
        D[cell]=vMainCoefficient[cell]/(*vRelaxationFactor) ;
      }

      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; } } ;

  register_rule<ComputeVelocityMatrixDiagonal>
    registerComputeVelocityMatrixDiagonal ;

  // Rule to copy vStar_D for periodic faces.
  class ComputeVelocityMatrixDiagonalPeriodic : public pointwise_rule {
    private:
      const_Map cl,cr,pmap ;
      store<real> D ;
    public:

      // Define input and output.
      ComputeVelocityMatrixDiagonalPeriodic() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pmap",pmap) ;
        name_store("vStar_D",D) ;
        input("pmap->cl->vStar_D") ;
        output("cr->vStar_D") ;
      }

      // For a face.
      void calculate(Entity face) { D[cr[face]]=D[cl[pmap[face]]] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; } } ;

  register_rule<ComputeVelocityMatrixDiagonalPeriodic>
    registerComputeVelocityMatrixDiagonalPeriodic ;

  // Rule to initialize the lower terms for the linear system. Checked.
  class InitializeVelocityMatrixLower : public unit_rule {
    private:
      store<real> L ;
    public:

      // Define input and output.
      InitializeVelocityMatrixLower() {
        name_store("vStar_L",L) ;
        output("vStar_L") ;
        constraint("internalFaces") ;
      }

      // Initialize for a single face.
      void calculate(Entity face) { L[face]=0.0 ; }

      // Initialize for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeVelocityMatrixLower>
    registerInitializeVelocityMatrixLower ;

  // Rule to add the first-order inviscid flux contribution to the lower terms
  // for the linear system. Checked.
  class FOUInviscidFluxToVelocityMatrixLower : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> thetaParameter ;
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      store<real> L ;
    public:

      // Define input and output.
      FOUInviscidFluxToVelocityMatrixLower() {
        name_store("thetaParameter",thetaParameter) ;
        name_store("massFlux",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("vStar_L",L) ;
        input("thetaParameter,massFlux,gridMassFlux") ;
        output("vStar_L") ;
        constraint("internalFaces") ;
      }

      // Increment the lower term for a single face. Note that the increment
      // is the negative of the one in streamUns since this coefficient is
      // for a term on the lhs of the equation. In streamUns, the coefficient
      // is associated with a term on the rhs.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        if(netMassFlux>0.0) L[face]-=(*thetaParameter)*netMassFlux ;
      }

      // Increment the upper term for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToVelocityMatrixLower>
    registerFOUInviscidFluxToVelocityMatrixLower ;

  // Rule to add the diffusive flux contribution to the lower terms for the
  // linear system. Right now we are assuming that we will be dealing with
  // turbulence models that interface to the flow equations via a turbulent
  // eddy viscosity.
  class DiffusiveFluxToVelocityMatrixLower : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> thetaParameter ;
      const_store<real> faceRadius ;
      const_store<real> diffusionProduct ;
      const_store<real> viscosity ;
      store<real> L ;
    public:

      // Define input and output.
      DiffusiveFluxToVelocityMatrixLower() {
        name_store("thetaParameter",thetaParameter) ;
        name_store("faceRadius",faceRadius) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("viscosity",viscosity) ;
        name_store("vStar_L",L) ;
        input("thetaParameter,faceRadius,diffusionProduct,viscosity") ;
        output("vStar_L") ;
        constraint("internalFaces,viscousFlow") ;
      }

      // Increment the lower term for a single face.
      void calculate(Entity face) {
        L[face]-=viscosity[face]*diffusionProduct[face]*(*thetaParameter)*
          faceRadius[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToVelocityMatrixLower>
    registerDiffusiveFluxToVelocityMatrixLower ;

  // Rule to initialize the upper terms for the linear system.
  class InitializeVelocityMatrixUpper : public unit_rule {
    private:
      store<real> U ;
    public:

      // Define input and output.
      InitializeVelocityMatrixUpper() {
        name_store("vStar_U",U) ;
        output("vStar_U") ;
        constraint("internalFaces") ;
      }

      // Initialize for a single face.
      void calculate(Entity face) { U[face]=0.0 ; }

      // Initialize for a sequence of faces.
      virtual void compute(const sequence &seq) {
        do_loop(seq,this,&InitializeVelocityMatrixUpper::calculate) ;
      }
  } ;

  register_rule<InitializeVelocityMatrixUpper>
    registerInitializeVelocityMatrixUpper ;

  // Rule to add the first-order inviscid flux contribution to the upper terms
  // for the linear system. Checked.
  class FOUInviscidFluxToVelocityMatrixUpper : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> thetaParameter ;
      const_store<real> massFlux ;
      const_store<real> gridMassFlux ;
      store<real> U ;
    public:

      // Define input and output.
      FOUInviscidFluxToVelocityMatrixUpper() {
        name_store("thetaParameter",thetaParameter) ;
        name_store("massFlux",massFlux) ;
        name_store("gridMassFlux",gridMassFlux) ;
        name_store("vStar_U",U) ;
        input("thetaParameter,massFlux,gridMassFlux") ;
        output("vStar_U") ;
        constraint("internalFaces") ;
      }

      // Increment the upper term for a single face. Contribution is opposite
      // to that in streamUns, as noted above for the lower terms.
      void calculate(Entity face) {
        real netMassFlux=massFlux[face]-gridMassFlux[face] ;
        if(netMassFlux<=0.0) U[face]+=(*thetaParameter)*netMassFlux ;
      }

      // Increment the upper term for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FOUInviscidFluxToVelocityMatrixUpper>
    registerFOUInviscidFluxToVelocityMatrixUpper ;

  // Rule to add the diffusive flux contribution to the upper terms for the
  // linear system. Right now we are assuming that we will be dealing with
  // turbulence models that interface to the flow equations via a turbulent
  // eddy viscosity.
  class DiffusiveFluxToVelocityMatrixUpper : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> thetaParameter ;
      const_store<real> faceRadius ;
      const_store<real> diffusionProduct ;
      const_store<real> viscosity ;
      store<real> U ;
    public:

      // Define input and output.
      DiffusiveFluxToVelocityMatrixUpper() {
        name_store("thetaParameter",thetaParameter) ;
        name_store("faceRadius",faceRadius) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("viscosity",viscosity) ;
        name_store("vStar_U",U) ;
        input("thetaParameter,faceRadius,diffusionProduct,viscosity") ;
        output("vStar_U") ;
        constraint("internalFaces,viscousFlow") ;
      }

      // Increment the lower term for a single face.
      void calculate(Entity face) {
        U[face]-=viscosity[face]*diffusionProduct[face]*(*thetaParameter)*
          faceRadius[face] ;
      }

      // Call calculate for a sequence of faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusiveFluxToVelocityMatrixUpper>
    registerDiffusiveFluxToVelocityMatrixUpper ;

  // Rule to compute the right-hand side for the linear system. Checked.
  class ComputeVelocityRHS : public pointwise_rule {
    private:
      const_param<real> twoDimensionFactor ;
      const_param<real> vRelaxationFactor ;
      const_store<vect3d> v ;
      const_store<real> vMainCoefficient ;
      const_store<vect3d> vSourceTerm ;
      store<vect3d> B ;
    public:

      // Define input and output.
      ComputeVelocityRHS() {
        name_store("twoDimensionFactor",twoDimensionFactor) ;
        name_store("vRelaxationFactor",vRelaxationFactor) ;
        name_store("v",v) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        name_store("vSourceTerm",vSourceTerm) ;
        name_store("vStar_B",B) ;
        input("twoDimensionFactor,vRelaxationFactor,v,vMainCoefficient") ;
        input("vSourceTerm") ;
        output("vStar_B") ;
        constraint("geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) {
        B[cell]=vSourceTerm[cell]+(1.0-(*vRelaxationFactor))*
          vMainCoefficient[cell]*v[cell]/(*vRelaxationFactor) ;
        B[cell].z*=(*twoDimensionFactor) ;
      }

      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeVelocityRHS> registerComputeVelocityRHS ;

  // Priority rule for inviscid flow to compute the right-hand side for the
  // linear system. Includes the net mass flux term which cancels the term
  // added to the main coefficient. Checked.
  class ComputeVelocityRHSInviscid : public pointwise_rule {
    private:
      const_param<real> vRelaxationFactor ;
      const_store<vect3d> v ;
      const_store<real> vMainCoefficient ;
      const_store<vect3d> vSourceTerm ;
      const_store<real> netMassFlux ;
      store<vect3d> B ;
    public:

      // Define input and output.
      ComputeVelocityRHSInviscid() {
        name_store("vRelaxationFactor",vRelaxationFactor) ;
        name_store("v",v) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        name_store("vSourceTerm",vSourceTerm) ;
        name_store("netMassFlux",netMassFlux) ;
        name_store("inviscidFlow::vStar_B",B) ;
        input("vRelaxationFactor,v,vMainCoefficient") ;
        input("vSourceTerm,netMassFlux") ;
        output("inviscidFlow::vStar_B") ;
        constraint("inviscidFlow,geom_cells") ;
      }

      // Add relaxation for a single cell.
      void calculate(Entity cell) {
        B[cell]=vSourceTerm[cell]+(1.0-(*vRelaxationFactor))*
          vMainCoefficient[cell]*v[cell]/(*vRelaxationFactor) ;
        if(netMassFlux[cell]<0.0) B[cell]-=netMassFlux[cell]*v[cell] ;
      }

      // Add relaxation for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeVelocityRHSInviscid> registerComputeVelocityRHSInviscid ;

  // Rule to copy vStar_B for periodic faces.
  class ComputeVelocityRHSPeriodic : public pointwise_rule {
    private:
      const_Map cl,cr,pmap ;
      store<vect3d> B ;
    public:

      // Define input and output.
      ComputeVelocityRHSPeriodic() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("pmap",pmap) ;
        name_store("vStar_B",B) ;
        input("pmap->cl->vStar_B") ;
        output("cr->vStar_B") ;
      }

      // For a face.
      void calculate(Entity face) { B[cr[face]]=B[cl[pmap[face]]] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeVelocityRHSPeriodic> registerComputeVelocityRHSPeriodic ;

//-----------------------------------------------------------------------------
// Rules for computing the residual of the momentum equation.

  // Rule to initialize the momentum residual scale factor.
  class InitializeMomentumResidualScaleFactor : public unit_rule {
    private:
      param<real> vResidualScaleFactor ;
    public:

      // Define input and output.
      InitializeMomentumResidualScaleFactor() {
        name_store("vResidualScaleFactor",vResidualScaleFactor) ;
        output("vResidualScaleFactor") ;
        constraint("UNIVERSE") ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { *vResidualScaleFactor=0.0 ; }
  } ;

//  register_rule<InitializeMomentumResidualScaleFactor>
//  registerInitializeMomentumResidualScaleFactor ;

  // Rule to sum for the momentum residual scale factor.
  class ComputeMomentumResidualScaleFactor : public apply_rule<param<real>,
  Loci::Summation<real> > {
    private:
      const_store<vect3d> v ;
      const_store<real> vMainCoefficient ;
      param<real> vResidualScaleFactor ;
    public:

      // Define input and output.
      ComputeMomentumResidualScaleFactor() {
        name_store("v",v) ;
        name_store("vMainCoefficient",vMainCoefficient) ;
        name_store("vResidualScaleFactor",vResidualScaleFactor) ;
        input("v,vMainCoefficient") ;
        output("vResidualScaleFactor") ;
      }

      // Add the cell value.
      void calculate(Entity cell) {
        join(*vResidualScaleFactor,vMainCoefficient[cell]*norm(v[cell])) ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

//register_rule<ComputeMomentumResidualScaleFactor>
//  registerComputeMomentumResidualScaleFactor ;

  // Rule to initialize the residual.
  class InitializeMomentumResidual : public unit_rule {
    private:
      store<vect3d> vResidual ;
    public:

      // Define input and output.
      InitializeMomentumResidual() {
        name_store("vResidual",vResidual) ;
        output("vResidual") ;
        constraint("vol") ;
      }

      // Initialize the residual for a single cell.
      void calculate(Entity cell) { vResidual[cell]=vect3d(0.0,0.0,0.0) ; }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitializeMomentumResidual> registerInitializeMomentumResidual ;

  // Rule to compute the momentum residual for each cell. Checked.
  class ComputeMomentumResidualOne : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
//    const_param<real> vResidualScaleFactor ;
      const_param<real> rhoScale,vScale,lScale ;
      const_store<real> D ;
      const_store<vect3d> v ;
      const_store<vect3d> B ;
      store<vect3d> vResidual ;
    private:
      real vFactor ;
    public:

      // Define input and output.
      ComputeMomentumResidualOne() {
//      name_store("vResidualScaleFactor",vResidualScaleFactor) ;
        name_store("rhoScale",rhoScale) ;
        name_store("vScale",vScale) ;
        name_store("lScale",lScale) ;
        name_store("vStar_D",D) ;
        name_store("v",v) ;
        name_store("vStar_B",B) ;
        name_store("vResidual",vResidual) ;
//      input("vResidualScaleFactor,vStar_D,v,vStar_B") ;
        input("rhoScale,vScale,lScale,vStar_D,v,vStar_B") ;
        output("vResidual") ;
        constraint("vol") ;
      }

      // Add the source and diagonal terms.
      void calculate(Entity cell) {
//      real temp=((*vResidualScaleFactor)==0.0)? 1.0:(*vResidualScaleFactor) ;
//      vResidual[cell]+=(B[cell]-D[cell]*v[cell])/temp ;
        vResidual[cell]+=(B[cell]-D[cell]*v[cell])/vFactor ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) {
        vFactor=(*rhoScale)*(*vScale)*(*vScale)*(*lScale)*(*lScale) ;
        do_loop(seq,this) ;
      }
  } ;

  register_rule<ComputeMomentumResidualOne> registerComputeMomentumResidualOne ;

  // Rule to compute the momentum residual for each cell. Checked.
  class ComputeMomentumResidualTwo : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
//    const_param<real> vResidualScaleFactor ;
      const_param<real> rhoScale,vScale,lScale ;
      const_store<vect3d> v ;
      const_store<real> L,U ;
      store<vect3d> vResidual ;
    private:
      real vFactor ;
    public:

      // Define input and output.
      ComputeMomentumResidualTwo() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
//      name_store("vResidualScaleFactor",vResidualScaleFactor) ;
        name_store("rhoScale",rhoScale) ;
        name_store("vScale",vScale) ;
        name_store("lScale",lScale) ;
        name_store("v",v) ;
        name_store("vStar_L",L) ;
        name_store("vStar_U",U) ;
        name_store("vResidual",vResidual) ;
//      input("vResidualScaleFactor,(cl,cr)->v,vStar_L,vStar_U") ;
        input("rhoScale,vScale,lScale,(cl,cr)->v,vStar_L,vStar_U") ;
        output("(cl,cr)->vResidual") ;
        constraint("internalFaces") ;
      }

      // Add the neighbor contribution to the residual for each of the two
      // cells on either side of the face.
      void calculate(Entity face) {
//      real temp=((*vResidualScaleFactor)==0.0)? 1.0:(*vResidualScaleFactor) ;
//      vResidual[cl[face]]-=U[face]*v[cr[face]]/temp ;
//      vResidual[cr[face]]-=L[face]*v[cl[face]]/temp ;
        vResidual[cl[face]]-=U[face]*v[cr[face]]/vFactor ;
        vResidual[cr[face]]-=L[face]*v[cl[face]]/vFactor ;
      }

      // Add the neighbor contributions for a sequence of faces.
      virtual void compute(const sequence &seq) {
        vFactor=(*rhoScale)*(*vScale)*(*vScale)*(*lScale)*(*lScale) ;
        do_loop(seq,this) ;
      }
  } ;

  register_rule<ComputeMomentumResidualTwo> registerComputeMomentumResidualTwo ;

  // Rule to initialize the total momentum residual. Checked.
  class InitializeTotalMomentumResidual : public unit_rule {
    private:
      param<VectorResidual> vResidualData ;
    public:

      // Define input and output.
      InitializeTotalMomentumResidual() {
        name_store("vResidualData",vResidualData) ;
        output("vResidualData") ;
        constraint("geom_cells") ;
      }

      // Initialize the residual for a sequence of cells.
      virtual void compute(const sequence &seq) {
        *vResidualData=VectorResidual() ;
      }
  } ;

  register_rule<InitializeTotalMomentumResidual>
    registerInitializeTotalMomentumResidual ;

  // Rule to compute the total momentum residual. Checked.
  class ComputeTotalMomentumResidual : public apply_rule<param<VectorResidual>,
  VectorResidualJoin> {
    private:
      const_store<vect3d> vResidual ;
      const_store<vect3d> cellCenter ;
      param<VectorResidual> vResidualData ;
    public:

      // Define input and output.
      ComputeTotalMomentumResidual() {
        name_store("vResidual",vResidual) ;
        name_store("cellcenter",cellCenter) ;
        name_store("vResidualData",vResidualData) ;
        input("vResidual,cellcenter") ;
        output("vResidualData") ;
        constraint("geom_cells") ;
      }

      // Add the contribution to the residual for a single cell.
      void calculate(Entity cell) {
        VectorResidual temp ; temp.maxResidual=vResidual[cell] ;
        temp.totalResidual=vect3d(abs(vResidual[cell].x),
          abs(vResidual[cell].y),abs(vResidual[cell].z)) ;
        temp.maxResidualLocation=cellCenter[cell] ;
        join(*vResidualData,temp) ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<ComputeTotalMomentumResidual>
    registerComputeTotalMomentumResidual ;

//-----------------------------------------------------------------------------
// Rules for marching the momentum equation.

  // Time build rule for velocity when using BDF2 time integrator. Although
  // this rule sets v{n=-1} from v_ic, these values are not really
  // used since BDF is used on the first timestep for non-restarts. The only
  // purpose for this rule is to let Loci know that there are two previous
  // time-levels that need to be maintained for BDF2.
  class TimeBuildVelocityBDF2 : public pointwise_rule {
    private:
      const_store<vect3d> v_ic ;
      store<vect3d> v ;
    public:

      // Define input and output.
      TimeBuildVelocityBDF2() {
        name_store("v_ic",v_ic) ;
        name_store("v{n=-1}",v) ;
        input("v_ic") ;
        output("v{n=-1}") ;
        constraint("geom_cells") ;
      }

      // Assign velocity for a single cell.
      void calculate(Entity cell) { v[cell]=v_ic[cell] ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildVelocityBDF2> registerTimeBuildVelocityBDF2 ;

  // Time build rule for velocity.
  class TimeBuildVelocity : public pointwise_rule {
    private:
      const_store<vect3d> v_ic ;
      store<vect3d> vTimeStepZero ;
    public:

      // Define input and output.
      TimeBuildVelocity() {
        name_store("v_ic",v_ic) ;
        name_store("v{n=0}",vTimeStepZero) ;
        input("v_ic") ;
        output("v{n=0}") ;
        constraint("geom_cells") ;
      }

      // Assign velocity at time zero for a single cell.
      void calculate(Entity cell) { vTimeStepZero[cell]=v_ic[cell] ; }

      // Assign velocity at time zero for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildVelocity> registerTimeBuildVelocity ;

  // Iteration build rule for velocity.
  class IterationBuildVelocity : public pointwise_rule {
    private:
      const_store<vect3d> vTimeStepN ;
      store<vect3d> vIterationZero ;
    public:

      // Define input and output.
      IterationBuildVelocity() {
        name_store("v{n}",vTimeStepN) ;
        name_store("v{n,it=0}",vIterationZero) ;
        input("v{n}") ;
        output("v{n,it=0}") ;
        constraint("geom_cells{n}") ;
      }

      // Assign velocity at iteration zero for a single cell.
      void calculate(Entity cell) { vIterationZero[cell]=vTimeStepN[cell] ; }

      // Assign velocity at iteration zero for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationBuildVelocity> registerIterationBuildVelocity ;

  // Rule to add temporal component of momentum equation to the main
  // coefficient. The use of old density and volume comes from subtracting off
  // continuity to get diagonal dominance.
  class TemporalToVelocityMainCoefficient: public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_param<real> timeIntegratorFactor0 ;
      const_store<real> timeStepFactor ;
      const_store<real> rho ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<real> vMainCoefficient ;
    public:

      // Define input and output.
      TemporalToVelocityMainCoefficient() {
        name_store("timeIntegratorFactor0{n}",timeIntegratorFactor0) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n}",rho) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("vMainCoefficient{n,it}",vMainCoefficient) ;
        input("rho{n},vol{n},cellRadius{n,it}") ;
        input("timeStepFactor{n},timeIntegratorFactor0{n}") ;
        output("vMainCoefficient{n,it}") ;
        constraint("geom_cells") ;
      }

      // Add temporal component for a single cell.
      void calculate(Entity cell) {
        vMainCoefficient[cell]+=rho[cell]*vol[cell]*cellRadius[cell]*
          (*timeIntegratorFactor0)/timeStepFactor[cell] ;
      }

      // Add temporal component for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToVelocityMainCoefficient>
    registerTemporalToVelocityMainCoefficient ;

  // Rule to add temporal component of momentum equation to the source term.
  // The use of old density and volume comes from subtracting off continuity
  // to get diagonal dominance.
  class TemporalToVelocitySourceTerm : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_param<real> timeIntegratorFactor0 ;
      const_store<real> timeStepFactor ;
      const_store<real> rho ;
      const_store<vect3d> v ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      TemporalToVelocitySourceTerm() {
        name_store("timeIntegratorFactor0{n}",timeIntegratorFactor0) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n}",rho) ;
        name_store("v{n}",v) ;
        name_store("vol{n}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("vSourceTerm{n,it}",vSourceTerm) ;
        input("rho{n},v{n},vol{n},cellRadius{n,it}") ;
        input("timeStepFactor{n},timeIntegratorFactor0{n}") ;
        output("vSourceTerm{n,it}") ;
        constraint("geom_cells") ;
      }

      // Add temporal component to source term for a single cell.
      void calculate(Entity cell) {
        vSourceTerm[cell]+=(rho[cell]*vol[cell]*cellRadius[cell]*
          (*timeIntegratorFactor0)/timeStepFactor[cell])*v[cell] ;
      }

      // Add temporal component to source term for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToVelocitySourceTerm>
    registerTemporalToVelocitySourceTerm ;

  // Rule to add temporal component of momentum equation to the source term
  // for the BDF2 scheme. Note that we are using vol{n,it} in here when we
  // should be using vol{n-1}. This is done so that we do not have to carry
  // aroung vol{n-1} for non-deforming mesh calculations. For deforming meshes,
  // we add another source to recover back vol{n-1} from vol{n,it}.
  class TemporalToVelocitySourceTermBDF2 : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_param<real> timeIntegratorFactor1 ;
      const_store<real> timeStepFactor ;
      const_store<real> rhoOld ;
      const_store<vect3d> vOld,v ;
      const_store<real> vol ;
      const_store<real> cellRadius ;
      store<vect3d> vSourceTerm ;
    public:

      // Define input and output.
      TemporalToVelocitySourceTermBDF2() {
        name_store("timeIntegratorFactor1{n}",timeIntegratorFactor1) ;
        name_store("timeStepFactor{n}",timeStepFactor) ;
        name_store("rho{n-1}",rhoOld) ;
        name_store("v{n-1}",vOld) ;
        name_store("v{n,it}",v) ;
        name_store("vol{n,it}",vol) ;
        name_store("cellRadius{n,it}",cellRadius) ;
        name_store("vSourceTerm{n,it}",vSourceTerm) ;
        input("rho{n-1},v{n-1},v{n,it},vol{n,it}") ;
        input("cellRadius{n,it},timeStepFactor{n},timeIntegratorFactor1{n}") ;
        output("vSourceTerm{n,it}") ;
        constraint("geom_cells,BDF2Integrator") ;
      }

      // Add temporal component to source term for a single cell.
      void calculate(Entity cell) {
        vSourceTerm[cell]+=((*timeIntegratorFactor1)*rhoOld[cell]*
          vol[cell]*cellRadius[cell]/timeStepFactor[cell])*
          (v[cell]-vOld[cell]) ;
      }

      // Add temporal component to source term for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TemporalToVelocitySourceTermBDF2>
    registerTemporalToVelocitySourceTermBDF2 ;

  // Iteration advance rule for velocity.
  class IterationAdvanceVelocity : public pointwise_rule {
    private:
      const_store<vect3d> vCorrected ;
      const_param<bool> iterationFinishedHack ;
      store<vect3d> vIterationPlusOne ;
    public:

      // Define input and output.
      IterationAdvanceVelocity() {
        name_store("vCorrected{n,it}",vCorrected) ;
        name_store("iterationFinished{n,it}",iterationFinishedHack) ;
        name_store("v{n,it+1}",vIterationPlusOne) ;
        input("vCorrected{n,it},iterationFinished{n,it}") ;
        output("v{n,it+1}") ;
        constraint("geom_cells{n,it}") ;
      }

      // Assign velocity at end of iteration for a single cell.
      void calculate(Entity cell) { vIterationPlusOne[cell]=vCorrected[cell] ; }

      // Assign velocity at end of iteration for a sequence of cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationAdvanceVelocity> registerIterationAdvanceVelocity ;

  // Build rule.
  class IterationFinishedBuild : public singleton_rule {
    private:
      param<bool> iterationFinished ;
    public:
                                                                                
      // Define input and output.
      IterationFinishedBuild() {
        name_store("iterationFinished{n,it=-1}",iterationFinished) ;
        output("iterationFinished{n,it=-1}") ;
        constraint("UNIVERSE{n}") ;
      }
                                                                                
      // Set to false.
      void compute(const sequence &seq) { *iterationFinished=false ; }
  } ;
                                                                                
  register_rule<IterationFinishedBuild> registerIterationFinishedBuild ;

  // Class to determine if iteration is finished.
  class CheckIterationFinished : public singleton_rule {
    private:
      const_param<int> numSpecies ;
      const_param<int> nCycle,it ;
      const_param<int> maxIterationsPerTimeStep ;
      const_param<real> convergenceTolerance ;
      const_param<VectorResidual> vResidualData ;
      const_param<ScalarResidual> pPrimeResidualData ;
      const_param<ScalarResidual> hResidualData ;
      const_param<ScalarResidual> kResidualData ;
      const_param<ScalarResidual> omegaResidualData ;
      const_param<vector<real> > totalSpeciesResidual ;
      param<bool> iterationFinished ;
    public:

      // Define input and output.
      CheckIterationFinished() {
        name_store("numSpecies{n,it}",numSpecies) ;
        name_store("ncycle{n}",nCycle) ;
        name_store("$it{n,it}",it) ;
        name_store("maxIterationsPerTimeStep{n,it}",maxIterationsPerTimeStep) ;
        name_store("convergenceTolerance{n,it}",convergenceTolerance) ;
        name_store("vResidualData{n,it}",vResidualData) ;
        name_store("pPrimeResidualData{n,it}",pPrimeResidualData) ;
        name_store("hResidualData{n,it}",hResidualData) ;
        name_store("kResidualData{n,it}",kResidualData) ;
        name_store("omegaResidualData{n,it}",omegaResidualData) ;
        name_store("totalSpeciesResidual{n,it}",totalSpeciesResidual) ;
        name_store("iterationFinished{n,it}",iterationFinished) ;
        input("numSpecies{n,it},ncycle{n},$it{n,it}") ;
        input("maxIterationsPerTimeStep{n,it},convergenceTolerance{n,it}") ;
        input("vResidualData{n,it},pPrimeResidualData{n,it}") ;
        input("hResidualData{n,it},kResidualData{n,it}") ;
        input("omegaResidualData{n,it},totalSpeciesResidual{n,it}") ;
        output("iterationFinished{n,it}") ;
        constraint("vResidualData{n,it}") ;
      }

      // Check if iteration is finished.
      void compute(const sequence &seq) {

        // Set the output format.
        if(Loci::MPI_rank==0){
          cout.setf(ios::scientific,ios::floatfield) ; cout.precision(6) ;
        }

        // Compute the maximum species residual.
        real maxTotalSpeciesResidual=0.0 ;
        if(*numSpecies>1){
//for(int i=0;i<*numSpecies;++i)
//cout << "species, residual: " << i << " " << (*totalSpeciesResidual)[i] << endl ;
          for(int i=0;i<*numSpecies;++i)
            if((*totalSpeciesResidual)[i]>maxTotalSpeciesResidual)
              maxTotalSpeciesResidual=(*totalSpeciesResidual)[i] ;
        }

        // Write out the residuals.
        if(Loci::MPI_rank==0){
          cout <<"R: " << *nCycle << " " << *it << " " << vResidualData->
            totalResidual << pPrimeResidualData->totalResidual ;
          if(hResidualData->totalResidual!=0.0)
            cout << " " << hResidualData->totalResidual ;
          if(kResidualData->totalResidual!=0.0)
            cout << " " << kResidualData->totalResidual << " "
              << omegaResidualData->totalResidual ;
          if(maxTotalSpeciesResidual!=0.0)
            cout << " " << maxTotalSpeciesResidual ;
          cout << endl ;
        }

        // See if we are converged.
        *iterationFinished=(*it==*maxIterationsPerTimeStep-1 ||
          (vResidualData->totalResidual.x<*convergenceTolerance &&
          vResidualData->totalResidual.y<*convergenceTolerance &&
          vResidualData->totalResidual.z<*convergenceTolerance &&
          pPrimeResidualData->totalResidual<*convergenceTolerance &&
          hResidualData->totalResidual<*convergenceTolerance &&
          kResidualData->totalResidual<*convergenceTolerance &&
          omegaResidualData->totalResidual<*convergenceTolerance &&
          maxTotalSpeciesResidual<*convergenceTolerance)) ;
      }
  } ;

  register_rule<CheckIterationFinished> registerCheckIterationFinished ;

  // Iteration collapse rule for velocity. Note that we have changed
  // this rule from "v{n+1}=v{n,it}" to the current form. With the old form
  // we could not do one iteration per time step, because v{n,it} would not
  // get updated and thus the residuals would never change.
  class IterationCollapseVelocity : public pointwise_rule {
    private:
      store<vect3d> v ;
    public:

      // Define input and output.
      IterationCollapseVelocity() {
        name_store("v{n,it}",v) ;
        input("v{n,it}") ;
        output("v{n+1}=v{n,it}") ;
        conditional("iterationFinished{n,it-1}") ;
        constraint("geom_cells{n,it}") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<IterationCollapseVelocity> registerIterationCollapseVelocity ;

  // Class to determine if time-stepping is finished.
  class CheckTimeSteppingFinished: public singleton_rule {
    private:
      const_param<int> n ;
      const_param<int> numTimeSteps ;
      param<bool> timeSteppingFinished ;
    public:

      // Define input and output.
      CheckTimeSteppingFinished() {
        name_store("$n",n) ;
        name_store("numTimeSteps",numTimeSteps) ;
        name_store("timeSteppingFinished",timeSteppingFinished) ;
        input("$n,numTimeSteps") ;
        output("timeSteppingFinished") ;
      }

      // Check if time-stepping is finished.
      void compute(const sequence &seq) {
        *timeSteppingFinished=(*n==*numTimeSteps) ;
      }
  } ;

  register_rule<CheckTimeSteppingFinished> registerCheckTimeSteppingFinished ;

  // Time collapse rule for incompressible flow.
  class TimeCollapseIncompressible : public pointwise_rule {
    private:
      const_store<real> rho,p ;
      const_store<vect3d> v ;
      store<int> solution ;
    public:

      // Define input and output.
      TimeCollapseIncompressible() {
        name_store("rho{n}",rho) ;
        name_store("v{n}",v) ;
        name_store("p{n}",p) ;
        name_store("solution",solution) ;
        input("rho{n},v{n},p{n}") ;
        output("solution") ;
        conditional("timeSteppingFinished{n}") ;
        constraint("incompressibleFlow{n},noSpeciesTransport{n}") ;
        constraint("geom_cells{n}") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<TimeCollapseIncompressible> registerTimeCollapseIncompressible ;

  // Time collapse rule for incompressible flow with species transport.
  class TimeCollapseIncompressibleSpeciesTransport : public pointwise_rule {
    private:
      const_store<real> rho,p ;
      const_store<vect3d> v ;
      const_storeVec<real> y ;
      store<int> solution ;
    public:

      // Define input and output.
      TimeCollapseIncompressibleSpeciesTransport() {
        name_store("rho{n}",rho) ;
        name_store("v{n}",v) ;
        name_store("p{n}",p) ;
        name_store("y{n}",y) ;
        name_store("solution",solution) ;
        input("rho{n},v{n},p{n},y{n}") ;
        output("solution") ;
        conditional("timeSteppingFinished{n}") ;
        constraint("incompressibleFlow{n},speciesTransport{n},geom_cells{n}") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<TimeCollapseIncompressibleSpeciesTransport>
    registerTimeCollapseIncompressibleSpeciesTransport ;

  // Time collapse rule for compressible flow.
  class TimeCollapseCompressible : public pointwise_rule {
    private:
      const_store<real> rho,p,temperature,h ;
      const_store<vect3d> v ;
      const_storeVec<real> y ;
      store<int> solution ;
    public:

      // Define input and output.
      TimeCollapseCompressible() {
        name_store("rho{n}",rho) ;
        name_store("v{n}",v) ;
        name_store("p{n}",p) ;
        name_store("temperature{n}",temperature) ;
        name_store("y{n}",y) ;
        name_store("h{n}",h) ;
        name_store("solution",solution) ;
        input("rho{n},v{n},p{n},temperature{n},h{n},y{n}") ;
        output("solution") ;
        conditional("timeSteppingFinished{n}") ;
        constraint("compressibleFlow,geom_cells{n}") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<TimeCollapseCompressible> registerTimeCollapseCompressible ;
}







