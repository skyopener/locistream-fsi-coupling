//-----------------------------------------------------------------------------
// Description: This file contains some of the basic rules common to all
//   grid movement schemes.
//-----------------------------------------------------------------------------

// Standard library includes.
#include <vector>
using std::vector ;

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// StreamUns includes.
#include "const.h"
#include "move.h"
#include "residual.h"
#include "sciTypes.h"
#include "varsFileInputs.h"

namespace streamUns {

  //---------------------------------------------------------------------------
  // Implementation for RigidBodyDisplacement.

  // Returns the serialized buffer size.
  int RigidBodyDisplacement::BufferSize() const {
    int bufferSize=17 ; return bufferSize ;
  }

  // Packs the data into a buffer.
  void RigidBodyDisplacement::PackBuffer(real *buffer,int size) {
    int i=0 ;
    buffer[i++]=sConst ; buffer[i++]=tMag ; buffer[i++]=tFreq ;
    buffer[i++]=tPhi ; buffer[i++]=rAlphaBar ; buffer[i++]=rMag ;
    buffer[i++]=rFreq ; buffer[i++]=rPhi ;
    buffer[i++]=tDir.x ; buffer[i++]=tDir.y ; buffer[i++]=tDir.z ;
    buffer[i++]=rCenter.x ; buffer[i++]=rCenter.y ; buffer[i++]=rCenter.z ;
    buffer[i++]=rAxis.x ; buffer[i++]=rAxis.y ; buffer[i++]=rAxis.z ;
  }

  // Unpacks the data from a buffer.
  void RigidBodyDisplacement::UnpackBuffer(real *buffer,int size) {
    int i=0 ;
    sConst=buffer[i++] ; tMag=buffer[i++] ; tFreq=buffer[i++] ;
    tPhi=buffer[i++] ; rAlphaBar=buffer[i++] ; rMag=buffer[i++] ;
    rFreq=buffer[i++] ; rPhi=buffer[i++] ;
    tDir.x=buffer[i++] ; tDir.y=buffer[i++] ; tDir.z=buffer[i++] ;
    rCenter.x=buffer[i++] ; rCenter.y=buffer[i++] ; rCenter.z=buffer[i++] ;
    rAxis.x=buffer[i++] ; rAxis.y=buffer[i++] ; rAxis.z=buffer[i++] ;
  }

  vect3d RigidBodyDisplacement::Value(real time,vect3d position) const {

    // Compute the rotation origin, which is the point on the axis of
    // rotation that when connected to this node is perpendicular to the
    // axis of rotation.
    vect3d r=position-rCenter ; real s0=dot(r,rAxis)/dot(rAxis,rAxis) ;
    vect3d rOrigin=rCenter+s0*rAxis ;

    // Compute normalized position vector from rotation center. Compute
    // vector perpendicular to postion vector and rotation axis.
    vect3d rr=position-rOrigin ; real R=norm(rr) ; vect3d rHat=(1.0/R)*rr ;
    vect3d tHat=cross(rHat,rAxis) ;

    // Compute displacement.
    real theta=rAlphaBar+rMag*sin(rFreq*time+rPhi) ;
    return (sConst+tMag*sin(tFreq*time+tPhi))*tDir+R*((cos(theta)-1.0)*rHat-
      sin(theta)*tHat) ;
  }

  //---------------------------------------------------------------------------
  // Implementation for BoundaryDisplacement.

  // Returns the serialized buffer size.
  int BoundaryDisplacement::BufferSize() const {
    int bufferSize=4+rigidBody.BufferSize() ; return bufferSize ;
  }

  // Overridden virtual method called by fact database when a fact of type
  // BoundaryDisplacment is read. Will never be used.
  istream& BoundaryDisplacement::Input(istream &in) {
    return in ;
  }

  // Packs the data into a buffer.
  void BoundaryDisplacement::PackBuffer(real *buffer,int size) {
    int i=0 ;
    buffer[i++]=real(type) ;
    buffer[i++]=constant.x ; buffer[i++]=constant.y ; buffer[i++]=constant.z ;
    rigidBody.PackBuffer(buffer+4,size) ;
  }

  // Prints the values. Will never be used.
  ostream& BoundaryDisplacement::Print(ostream &out) const {
    return out ;
  }

  // Unpacks the data from a buffer.
  void BoundaryDisplacement::UnpackBuffer(real *buffer,int size) {
    int i=0 ;
    type=int(buffer[i++]) ;
    constant.x=buffer[i++] ; constant.y=buffer[i++] ; constant.z=buffer[i++] ;
    rigidBody.UnpackBuffer(buffer+4,size) ;
  }

  //---------------------------------------------------------------------------
  // Setup Rules.

  // Set the grid motion type to time-dependent when we are using rigid-body
  // grid motion.
  class RigidBodyGridMotionType : public constraint_rule {
    private:
      const_param<RigidBodyGridMotionOptions> rigidBodyGridMotion ;
      Constraint gridMotionTimeDependent,gridMotionSolutionDependent ;
    public:

      // Define input and output.
      RigidBodyGridMotionType() {
        name_store("rigidBodyGridMotion",rigidBodyGridMotion) ;
        name_store("gridMotionTimeDependent",gridMotionTimeDependent) ;
        name_store("gridMotionSolutionDependent",gridMotionSolutionDependent) ;
        input("rigidBodyGridMotion") ;
        output("gridMotionTimeDependent,gridMotionSolutionDependent") ;
      }

      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        gridMotionTimeDependent=~EMPTY ; gridMotionSolutionDependent=EMPTY ;
      }
  } ;

  register_rule<RigidBodyGridMotionType> registerRigidBodyGridMotionType ;

  // The following unit/apply sequence of rules identifies if the mesh
  // deformation is time-dependent or solution dependent. The value is 1
  // for time-dependent and 0 for solution-dependent. If there are any
  // boundary conditions that are solution-dependent, then the logial AND
  // operation for the apply rule will cause this value to be false.
  class MeshDeformationTimeDependentUnit : public unit_rule {
    private:
      const_param<string> gridMover ;
      param<int> meshDeformationTimeDependent ;
    public:
                                                                                
      // Define input and output.
      MeshDeformationTimeDependentUnit() {
        name_store("gridMover",gridMover) ;
        name_store("meshDeformationTimeDependent",meshDeformationTimeDependent) ;
        input("gridMover") ;
        output("meshDeformationTimeDependent") ;
      }
                                                                                
      // Set the default value.
      virtual void compute(const sequence &seq) {
        *meshDeformationTimeDependent=1 ;
      }
  } ;
                                                                                
  register_rule<MeshDeformationTimeDependentUnit>
    registerMeshDeformationTimeDependentUnit ;

  // The apply rule that checks for solution-dependent deforming boundaries.
  class MeshDeformationTimeDependentApply : public apply_rule<param<int>,
  Loci::Product<int> > {
    private:
      const_store<BoundaryDisplacement> s_BC ;
      param<int> meshDeformationTimeDependent ;
    public:
                                                                                
      // Define input and output.
      MeshDeformationTimeDependentApply() {
        name_store("s_BC",s_BC) ;
        name_store("meshDeformationTimeDependent",meshDeformationTimeDependent) ;
        input("s_BC") ;
        output("meshDeformationTimeDependent") ;
        constraint("s_BCoption") ;
      }

      // Look at a single boundary.
      void calculate(Entity e) {
        int flag=(s_BC[e].IsSolutionDependent())? 0:1 ;
        join(*meshDeformationTimeDependent,flag) ;
      }
                                                                                
      // Loop over the boundaries.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<MeshDeformationTimeDependentApply>
    registerMeshDeformationTimeDependentApply ;

  // Set the grid motion type for deforming meshes.
  class DeformingMeshGridMotionType : public constraint_rule {
    private:
      const_param<int> meshDeformationTimeDependent ;
      Constraint gridMotionTimeDependent,gridMotionSolutionDependent ;
    public:

      // Define input and output.
      DeformingMeshGridMotionType() {
        name_store("meshDeformationTimeDependent",meshDeformationTimeDependent);
        name_store("gridMotionTimeDependent",gridMotionTimeDependent) ;
        name_store("gridMotionSolutionDependent",gridMotionSolutionDependent) ;
        input("meshDeformationTimeDependent") ;
        output("gridMotionTimeDependent,gridMotionSolutionDependent") ;
      }

      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        gridMotionTimeDependent=(*meshDeformationTimeDependent==1)?
          ~EMPTY:EMPTY ;
        gridMotionSolutionDependent=(*meshDeformationTimeDependent==0)?
          ~EMPTY:EMPTY ;
      }
  } ;

  register_rule<DeformingMeshGridMotionType> registerDeformingMeshGridMotionType ;

  // Rule to declare "chi" as valid .vars file variable. This variable is used
  // in the Stein-Tezduyar algorithm.
  class DefaultGridMoverChi : public default_rule {
    private:
      param<real> gridMoverChi ;
    public:
                                                                                
      // Define input and output.
      DefaultGridMoverChi() {
        name_store("gridMoverChi",gridMoverChi) ;
        output("gridMoverChi") ;
      }
                                                                                
      // Set the default value.
      virtual void compute(const sequence &seq) { *gridMoverChi=1.0 ; }
  } ;
                                                                                
  register_rule<DefaultGridMoverChi> registerDefaultGridMoverChi ;

  // Rule to declare "kappa" as valid .vars file variable. This variable is used
  // as a multiplier for the rigid-body rotation term that has been added to the
  // original Stein-Tezduyar algorithm.
  class DefaultGridMoverKappa : public default_rule {
    private:
      param<real> gridMoverKappa ;
    public:
                                                                                
      // Define input and output.
      DefaultGridMoverKappa() {
        name_store("gridMoverKappa",gridMoverKappa) ;
        output("gridMoverKappa") ;
      }
                                                                                
      // Set the default value.
      virtual void compute(const sequence &seq) { *gridMoverKappa=1.0 ; }
  } ;
                                                                                
  register_rule<DefaultGridMoverKappa> registerDefaultGridMoverKappa ;

  // Rule to declare "gridMoverLinearSolver" as valid .vars file variable.
  class DefaultGridMoverLinearSolver : public default_rule {
    private:
      param<string> gridMoverLinearSolver ;
    public:
                                                                                
      // Define input and output.
      DefaultGridMoverLinearSolver() {
        name_store("gridMoverLinearSolver",gridMoverLinearSolver) ;
        output("gridMoverLinearSolver") ;
      }
                                                                                
      // Set the default value.
      virtual void compute(const sequence &seq) {
        *gridMoverLinearSolver="PETSC" ;
      }
  } ;
                                                                                
  register_rule<DefaultGridMoverLinearSolver>
    registerDefaultGridMoverLinearSolver ;

  // Rule to declare "gridMoverMaxIterationsPerTimeStep" as valid .vars file
  // variable.
  class DefaultGridMoverMaxIterationsPerTimeStep : public default_rule {
    private:
      param<int> gridMoverMaxIterationsPerTimeStep ;
    public:
                                                                                
      // Define input and output.
      DefaultGridMoverMaxIterationsPerTimeStep() {
        name_store("gridMoverMaxIterationsPerTimeStep",
          gridMoverMaxIterationsPerTimeStep) ;
        output("gridMoverMaxIterationsPerTimeStep") ;
      }
                                                                                
      // Set the default value.
      virtual void compute(const sequence &seq) {
        *gridMoverMaxIterationsPerTimeStep=100 ;
      }
  } ;
                                                                                
  register_rule<DefaultGridMoverMaxIterationsPerTimeStep>
    registerDefaultGridMoverMaxIterationsPerTimeStep ;

  // Rule to declare "gridMoverMaxLinearSolverIterations" as valid .vars file
  // variable. This variable is passed to the linear solver.
  class DefaultGridMoverMaxLinearSolverIterations : public default_rule {
    private:
      param<int> gridMoverMaxLinearSolverIterations ;
    public:
                                                                                
      // Define input and output.
      DefaultGridMoverMaxLinearSolverIterations() {
        name_store("gridMoverMaxLinearSolverIterations",
          gridMoverMaxLinearSolverIterations) ;
        output("gridMoverMaxLinearSolverIterations") ;
      }
                                                                                
      // Set the default value.
      virtual void compute(const sequence &seq) {
        *gridMoverMaxLinearSolverIterations=10 ;
      }
  } ;
                                                                                
  register_rule<DefaultGridMoverMaxLinearSolverIterations>
    registerDefaultGridMoverMaxLinearSolverIterations ;

  // Rule to declare "gridMoverNu" as valid .vars file variable. This variable
  // represents Poisson's ratio and is used in the Stein-Tezduyar algorithm.
  class DefaultGridMoverNu : public default_rule {
    private:
      param<real> gridMoverNu ;
    public:

      // Define input and output.
      DefaultGridMoverNu() {
        name_store("gridMoverNu",gridMoverNu) ;
        output("gridMoverNu") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) { *gridMoverNu=0.0 ; }
  } ;

  register_rule<DefaultGridMoverNu> registerDefaultGridMoverNu ;

  // Rule to declare "gridMoverRelaxation" as valid .vars file variable.
  class DefaultGridMoverRelaxation : public default_rule {
    private:
      param<real> gridMoverRelaxation ;
    public:
                                                                                
      // Define input and output.
      DefaultGridMoverRelaxation() {
        name_store("gridMoverRelaxation",gridMoverRelaxation) ;
        output("gridMoverRelaxation") ;
      }
                                                                                
      // Set the default value.
      virtual void compute(const sequence &seq) { *gridMoverRelaxation=0.3 ; }
  } ;
                                                                                
  register_rule<DefaultGridMoverRelaxation> registerDefaultGridMoverRelaxation ;

  // Rule to declare "gridMoverTolerance" as valid .vars file variable.
  // This variable is passed to the linear solver.
  class DefaultGridMoverTolerance : public default_rule {
    private:
      param<real> gridMoverTolerance ;
    public:
                                                                                
      // Define input and output.
      DefaultGridMoverTolerance() {
        name_store("gridMoverTolerance",gridMoverTolerance) ;
        output("gridMoverTolerance") ;
      }
                                                                                
      // Set the default value.
      virtual void compute(const sequence &seq) {
        *gridMoverTolerance=1.0e-03 ;
      }
  } ;
                                                                                
  register_rule<DefaultGridMoverTolerance> registerDefaultGridMoverTolerance ;

  // Rule to declare "gridMover" as valid .vars file variable.
  class OptionalGridMover : public optional_rule {
    private:
      param<string> gridMover ;
    public:
                                                                                
      // Define input and output.
      OptionalGridMover() {
        name_store("gridMover",gridMover) ;
        output("gridMover") ;
      }
                                                                                
      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;
                                                                                
  register_rule<OptionalGridMover> registerOptionalGridMover ;

  // Rule to declare "gridMoverChiMax" as valid .vars file variable. This
  // variable is used in computing a spatially varying chi for each edge.
  class OptionalGridMoverMaxChi : public optional_rule {
    private:
      param<real> gridMoverMaxChi ;
    public:

      // Define input and output.
      OptionalGridMoverMaxChi() {
        name_store("gridMoverMaxChi",gridMoverMaxChi) ;
        output("gridMoverMaxChi") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<OptionalGridMoverMaxChi> registerOptionalGridMoverMaxChi ;

  // Rule to declare "gridMoverChiMin" as valid .vars file variable. This
  // variable is used in computing a spatially varying chi for each edge.
  class OptionalGridMoverMinChi : public optional_rule {
    private:
      param<real> gridMoverMinChi ;
    public:

      // Define input and output.
      OptionalGridMoverMinChi() {
        name_store("gridMoverMinChi",gridMoverMinChi) ;
        output("gridMoverMinChi") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<OptionalGridMoverMinChi> registerOptionalGridMoverMinChi ;
                                                                                
  // Creates the linear solver constraint.
  class GridMoverLinearSolverConstraint : public constraint_rule {
    private:
      const_param<string> gridMoverLinearSolver ;
      Constraint s_JacobiLinearSolver,s_PETSCLinearSolver,s_HYPRELinearSolver ;
    public:
                                                                                
      // Define input and output.
      GridMoverLinearSolverConstraint() {
        name_store("gridMoverLinearSolver",gridMoverLinearSolver) ;
        name_store("sStar_JacobiLinearSolver",s_JacobiLinearSolver) ;
        name_store("sStar_PETSCLinearSolver",s_PETSCLinearSolver) ;
        name_store("sStar_HYPRELinearSolver",s_HYPRELinearSolver) ;
        input("gridMoverLinearSolver") ;
        output("sStar_JacobiLinearSolver,sStar_PETSCLinearSolver") ;
        output("sStar_HYPRELinearSolver") ;
      }
                                                                                
      // Set up the constraint.
      virtual void compute(const sequence& seq) {
        if(*gridMoverLinearSolver=="Jacobi"){
          s_JacobiLinearSolver=~EMPTY ; s_PETSCLinearSolver=EMPTY ;
          s_HYPRELinearSolver=EMPTY ;
        }else if(*gridMoverLinearSolver=="PETSC"){
          s_JacobiLinearSolver=EMPTY ; s_PETSCLinearSolver=~EMPTY ;
          s_HYPRELinearSolver=EMPTY ;
        }else if(*gridMoverLinearSolver=="HYPRE"){
          s_JacobiLinearSolver=EMPTY ; s_PETSCLinearSolver=EMPTY ;
          s_HYPRELinearSolver=~EMPTY ;
        }else{
          cerr << "ERROR: Bad linear solver type for grid mover." << endl ;
          Loci::Abort() ;
        }
      }
  } ;
                                                                                
  register_rule<GridMoverLinearSolverConstraint>
    registerGridMoverLinearSolverConstraint ;

  // Creates the maximum iterations variable used by the linear solvers.
  class LinearSolverMaxIterations : public singleton_rule {
    private:
      const_param<int> gridMoverMaxLinearSolverIterations ;
      param<int> s_maxLinearSolverIterations ;
    public:
                                                                                
      // Define input and output.
      LinearSolverMaxIterations() {
        name_store("gridMoverMaxLinearSolverIterations",
          gridMoverMaxLinearSolverIterations) ;
        name_store("sStar_maxLinearSolverIterations",
          s_maxLinearSolverIterations) ;
        input("gridMoverMaxLinearSolverIterations") ;
        output("sStar_maxLinearSolverIterations") ;
      }
                                                                                
      // Set the value.
      virtual void compute(const sequence &seq) {
        *s_maxLinearSolverIterations=*gridMoverMaxLinearSolverIterations ;
      }
  } ;
                                                                                
  register_rule<LinearSolverMaxIterations> registerLinearSolverMaxIterations ;

  // Creates the tolerance variable used by the linear solvers.
  class LinearSolverTolerance : public singleton_rule {
    private:
      const_param<real> gridMoverTolerance ;
      param<real> s_linearSolverTolerance ;
    public:
                                                                                
      // Define input and output.
      LinearSolverTolerance() {
        name_store("gridMoverTolerance",gridMoverTolerance) ;
        name_store("sStar_linearSolverTolerance",s_linearSolverTolerance) ;
        input("gridMoverTolerance") ;
        output("sStar_linearSolverTolerance") ;
      }
                                                                                
      // Set the value.
      virtual void compute(const sequence &seq) {
        *s_linearSolverTolerance=*gridMoverTolerance ;
      }
  } ;
                                                                                
  register_rule<LinearSolverTolerance> registerLinearSolverTolerance ;

//-----------------------------------------------------------------------------
// Priority rule to compute initial face areas when there is grid movement.

  class InitialFaceCenterGridMotion : public pointwise_rule {
    private:
      const_multiMap face2node ;
      const_store<vect3d> pos ;
      store<vect3d> facecenter ;
    public :

      // Define input and output.
     InitialFaceCenterGridMotion() {
        name_store("face2node",face2node) ;
        name_store("pos_ic",pos) ;
        name_store("facecenter_ic",facecenter) ;
        input("face2node->pos_ic") ;
        output("facecenter_ic") ;
     }

     // Compute center for a single face.
     void calculate(Entity face) {
       vect3d nodesum(0.0,0.0,0.0) ; real lensum=0.0 ;
       for(const int *id=face2node.begin(face);id+1!=face2node.end(face);++id) {
         vect3d pos0=pos[*id],pos1=pos[*(id+1)] ;
         vect3d edge_loc=0.5*(pos0+pos1),edge_vec=pos0-pos1 ;
         real len=sqrt(dot(edge_vec,edge_vec)) ;
         nodesum+=len*edge_loc ; lensum += len ;
       }
       const int *id=face2node.begin(face),*idend=face2node.end(face)-1 ;
       vect3d pos0=pos[*id],pos1=pos[*idend] ;
       vect3d edge_loc=0.5*(pos0+pos1),edge_vec=pos0-pos1 ;
       real len=sqrt(dot(edge_vec,edge_vec)) ;
       nodesum+=len*edge_loc ; lensum += len ;
       facecenter[face]=nodesum/lensum ;
     }

     // Loop through faces.
     virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitialFaceCenterGridMotion>
    registerInitialFaceCenterGridMotion ;

  class InitialAreaGridMotion : public pointwise_rule {
    private:
      const_store<vect3d> pos_ic ;
      const_store<vect3d> facecenter_ic ;
      const_multiMap face2node ;
      store<Area> area_ic ;
    public:
                                                                                
      // Define input and output.
      InitialAreaGridMotion() {
        name_store("pos_ic",pos_ic) ;
        name_store("face2node",face2node) ;
        name_store("facecenter_ic",facecenter_ic) ;
        name_store("gridMotion::area_ic",area_ic) ;
        input("facecenter_ic,face2node->pos_ic") ;
        output("gridMotion::area_ic") ;
      }
                                                                                
      // Area for a single face.
      void calculate(Entity fc) {
        const vect3d center=facecenter_ic[fc] ;
        const int first=*(face2node.begin(fc)) ;
        const int last=*(face2node.end(fc)-1) ;
        vect3d sum(cross(pos_ic[last]-center,pos_ic[first]-center)) ;
        for(const int* ni=face2node.begin(fc)+1;ni!=(face2node.end(fc));ni++) {
          sum+=cross(pos_ic[*(ni-1)]-center,pos_ic[*(ni)]-center) ;
        }
        area_ic[fc].n=sum ;
        const real sada=sqrt(dot(area_ic[fc].n,area_ic[fc].n)) ;
        area_ic[fc].n*=1./(sada+EPSILON) ; area_ic[fc].sada=0.5*sada ;
      }
                                                                                
      // Loop over faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<InitialAreaGridMotion> registerInitialAreaGridMotion ;

//-----------------------------------------------------------------------------
// Rules to allow rigidBodyGridMotion in the .vars file.

  // Options for the rigid-body grid mover.
  class OptionalRigidBodyGridMotion : public optional_rule {
    private:
      param<RigidBodyGridMotionOptions> rigidBodyGridMotion ;
    public:
                                                                                
      // Define input and output.
      OptionalRigidBodyGridMotion() {
        name_store("rigidBodyGridMotion",rigidBodyGridMotion) ;
        output("rigidBodyGridMotion") ;
      }
                                                                                
      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;
                                                                                
  register_rule<OptionalRigidBodyGridMotion>
    registerOptionalRigidBodyGridMotion ;

//-----------------------------------------------------------------------------
// Rules to process rigidBodyGridMotion options in the .vars file.

  // Creates the rigid-body grid motion parameters.
  class RigidBodyGridMotionParameters : public singleton_rule {
    private:
      const_param<RigidBodyGridMotionOptions> rigidBodyGridMotion ;
      param<vect3d> sDir,tDir,rCenter,rAxis ;
      param<real> sMag,tMag,tFreq,tPhi,rAlphaBar,rMag,rFreq,rPhi,func ;
    public:
                                                                                
      // Define input and output.
      RigidBodyGridMotionParameters() {
        name_store("rigidBodyGridMotion",rigidBodyGridMotion) ;
        name_store("rigidBodyGridMotionTDir",tDir) ;
        name_store("rigidBodyGridMotionRCenter",rCenter) ;
        name_store("rigidBodyGridMotionRAxis",rAxis) ;
        name_store("rigidBodyGridMotionTMag",tMag) ;
        name_store("rigidBodyGridMotionTFreq",tFreq) ;
        name_store("rigidBodyGridMotionTPhi",tPhi) ;
        name_store("rigidBodyGridMotionRAlphaBar",rAlphaBar) ;
        name_store("rigidBodyGridMotionRMag",rMag) ;
        name_store("rigidBodyGridMotionRFreq",rFreq) ;
        name_store("rigidBodyGridMotionRPhi",rPhi) ;
        name_store("rigidBodyGridMotionSMag",sMag) ;
        name_store("rigidBodyGridMotionSDir",sDir) ;
	name_store("rigidBodyGridMotionFunction",func) ;
        input("rigidBodyGridMotion") ;
	output("rigidBodyGridMotionFunction") ;
        output("rigidBodyGridMotionTDir,rigidBodyGridMotionRCenter") ;
        output("rigidBodyGridMotionRAxis,rigidBodyGridMotionTMag") ;
        output("rigidBodyGridMotionTFreq,rigidBodyGridMotionTPhi") ;
        output("rigidBodyGridMotionRAlphaBar,rigidBodyGridMotionRMag") ;
        output("rigidBodyGridMotionRFreq,rigidBodyGridMotionRPhi") ;
        output("rigidBodyGridMotionSMag,rigidBodyGridMotionSDir") ;
      }
                                                                                
      // Set paramter values.
      void compute(const sequence &seq) {

        // Set default values.
        *tDir=vect3d(0.0,1.0,0.0) ; *rCenter=vect3d(0.0,0.0,0.0) ;
        *rAxis=vect3d(0.0,0.0,1.0) ; *tMag=0.0 ; *tFreq=0.0 ; *tPhi=0.0 ;
        *rAlphaBar=0.0 ; *rMag=0.0 ; *rFreq=0.0 ; *rPhi=0.0 ;
        *sDir=vect3d(-1.0,0.0,0.0) ; *sMag=0.0 ; *func=1.0 ;

        // Set values of existing options.
        if((*rigidBodyGridMotion).optionExists("tDir"))
          get_vect3dOption(*rigidBodyGridMotion,"tDir","",*tDir) ;
        if((*rigidBodyGridMotion).optionExists("rCenter"))
          get_vect3dOption(*rigidBodyGridMotion,"rCenter","m",*rCenter) ;
        if((*rigidBodyGridMotion).optionExists("rAxis"))
          get_vect3dOption(*rigidBodyGridMotion,"rAxis","",*rAxis) ;
        if((*rigidBodyGridMotion).optionExists("tMag"))
          (*rigidBodyGridMotion).getOptionUnits("tMag","m",*tMag) ;
        if((*rigidBodyGridMotion).optionExists("tFreq"))
          (*rigidBodyGridMotion).getOptionUnits("tFreq","rad/s",*tFreq) ;
        if((*rigidBodyGridMotion).optionExists("tPhi"))
          (*rigidBodyGridMotion).getOptionUnits("tPhi","rad",*tPhi) ;
        if((*rigidBodyGridMotion).optionExists("rAlphaBar"))
          (*rigidBodyGridMotion).getOptionUnits("rAlphaBar","rad",*rAlphaBar) ;
        if((*rigidBodyGridMotion).optionExists("rMag"))
          (*rigidBodyGridMotion).getOptionUnits("rMag","rad",*rMag) ;
        if((*rigidBodyGridMotion).optionExists("rFreq"))
          (*rigidBodyGridMotion).getOptionUnits("rFreq","rad/s",*rFreq) ;
        if((*rigidBodyGridMotion).optionExists("rPhi"))
          (*rigidBodyGridMotion).getOptionUnits("rPhi","rad",*rPhi) ;
        if((*rigidBodyGridMotion).optionExists("sDir"))
          get_vect3dOption(*rigidBodyGridMotion,"sDir","",*sDir) ;
        if((*rigidBodyGridMotion).optionExists("sMag"))
          (*rigidBodyGridMotion).getOptionUnits("sMag","m/s",*sMag) ;
	if((*rigidBodyGridMotion).optionExists("func"))
          (*rigidBodyGridMotion).getOptionUnits("func","",*func) ;

        // Make sure direction vector do not have zero magnitude.
        if(norm(*tDir)==0.0){
          cerr << "ERROR: tDir cannot have zero magnitude!" << endl ;
          Loci::Abort() ;
        }
        if(norm(*rAxis)==0.0){
          cerr << "ERROR: rAxis cannot have zero magnitude!" << endl ;
          Loci::Abort() ;
        }
        if(norm(*sDir)==0.0){
          cerr << "ERROR: sDir cannot have zero magnitude!" << endl ;
          Loci::Abort() ;
        }

        // Normalize the translation direction and rotation axis.
        *tDir*=1.0/norm(*tDir) ; *rAxis*=1.0/norm(*rAxis) ;
        *sDir*=1.0/norm(*sDir) ;
      }
  } ;

  register_rule<RigidBodyGridMotionParameters>
    registerRigidBodyGridMotionParameters ;

//-----------------------------------------------------------------------------
// Initial conditions for node displacement.

  // Default node displacement initial condition.
  class NodeDisplacementInitialCondition : public pointwise_rule {
    private:
      store<vect3d> s_ic ;
    public:

      // Define input and output.
      NodeDisplacementInitialCondition() {
        name_store("node_s_ic",s_ic) ;
        output("node_s_ic") ;
        constraint("pos,noRestart") ;
      }

      // Set displacement for single node.
      void calculate(Entity n) { s_ic[n]=vect3d(0.0,0.0,0.0) ; }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<NodeDisplacementInitialCondition>
    registerNodeDisplacementInitialCondition ;

  // Priority node displacement initial condition for rigid-body motion of grid.
  class NodeDisplacementInitialConditionPriority : public pointwise_rule {
    private:
      const_param<vect3d> tDir,rCenter,rAxis ;
      const_param<real> tMag,tPhi,rAlphaBar,rMag,rPhi,func ;
      const_store<vect3d> pos ;
      store<vect3d> s_ic ;
    private:
      vect3d sTranslation ;
      real theta ;
    public:

      // Define input and output.
      NodeDisplacementInitialConditionPriority() {
        name_store("rigidBodyGridMotionTDir",tDir) ;
        name_store("rigidBodyGridMotionTMag",tMag) ;
        name_store("rigidBodyGridMotionTPhi",tPhi) ;
        name_store("rigidBodyGridMotionRCenter",rCenter) ;
        name_store("rigidBodyGridMotionRAxis",rAxis) ;
        name_store("rigidBodyGridMotionRAlphaBar",rAlphaBar) ;
        name_store("rigidBodyGridMotionRMag",rMag) ;
        name_store("rigidBodyGridMotionRPhi",rPhi) ;
        name_store("pos",pos) ;
        name_store("rigidBodyGridMotion::node_s_ic",s_ic) ;
	name_store("rigidBodyGridMotionFunction",func) ;
	input("rigidBodyGridMotionFunction") ;
        input("rigidBodyGridMotionTDir,rigidBodyGridMotionTMag") ;
        input("rigidBodyGridMotionTPhi") ;
        input("rigidBodyGridMotionRCenter,rigidBodyGridMotionRAxis") ;
        input("rigidBodyGridMotionRAlphaBar,rigidBodyGridMotionRMag") ;
        input("rigidBodyGridMotionRPhi") ;
        input("pos") ;
        output("rigidBodyGridMotion::node_s_ic") ;
        constraint("pos,noRestart,rigidBodyGridMotion") ;
      }

      // Set displacement for single node.
      void calculate(Entity n) {

        // Compute the rotation origin, which is the point on the axis of
        // rotation that when connected to this node is perpendicular to the
        // axis of rotation.
        vect3d r=pos[n]-(*rCenter) ; real s0=dot(r,*rAxis)/dot(*rAxis,*rAxis) ;
        vect3d rOrigin=(*rCenter)+s0*(*rAxis) ;

        // Compute normalized position vector from rotation origin. Compute
        // vector perpendicular to postion vector and rotation axis. If a
        // node falls on the axis, there is no rotational displacement.
        vect3d rr=pos[n]-rOrigin ; real R=norm(rr) ;
        if(R<1.0e-12){ s_ic[n]=sTranslation ; return ; }
        vect3d rHat=(1.0/R)*rr ;
        vect3d tHat=cross(rHat,*rAxis) ;

        // Compute initial displacement.
        s_ic[n]=sTranslation+R*((cos(theta)-1.0)*rHat-sin(theta)*tHat) ;
      }

      // Loop through nodes.
      void compute(const sequence &seq) {
	int functionType=0 ;
	functionType = int(*func) ;
	if (functionType == 1) {
	  // cout << "functionType == 1" << endl ;
	  sTranslation = ((*tMag)*sin(*tPhi))*(*tDir) ;
	} else if (functionType == 3) {
	  // cout << "functionType == 3" << endl ;
	  sTranslation = ((*tMag)*(cos(*tPhi)-1.0))*(*tDir) ;
	}
        //sTranslation=((*tMag)*sin(*tPhi))*(*tDir) ;
        theta=(*rAlphaBar)+(*rMag)*sin(*rPhi) ;
        do_loop(seq,this) ;
      }
  } ;

  register_rule<NodeDisplacementInitialConditionPriority>
    registerNodeDisplacementInitialConditionPriority ;

  // Node displacement initial condition for restart.
  class NodeDisplacementInitialConditionRestart : public pointwise_rule {
    private:
      const_param<string> restartNum ;
      store<vect3d> s_ic ;
    public:

      // Define input and output.
      NodeDisplacementInitialConditionRestart() {
        name_store("restartNum",restartNum) ;
        name_store("node_s_ic",s_ic) ;
        input("restartNum") ;
        output("node_s_ic") ;
        constraint("restart,pos") ;
      }
                                                                                
      // Read the node displacement values from file.
      void compute(const sequence &seq) {
        ostringstream oss ; oss << "restart/node_s_hdf5." << *restartNum ;
        string fileName=oss.str() ;
        if(Loci::MPI_rank==0) cout << "Reading node displacements "
          << "from restart file: " << fileName << endl ;
        entitySet dom=entitySet(seq) ;
        hid_t fileID=Loci::hdf5OpenFile(fileName.c_str(),H5F_ACC_RDONLY,
          H5P_DEFAULT);
        Loci::readContainer(fileID,"node_s",s_ic.Rep(),dom) ;
        Loci::hdf5CloseFile(fileID) ;
      }
  } ;
                                                                                
  register_rule<NodeDisplacementInitialConditionRestart>
    registerNodeDisplacementInitialConditionRestart ;

//-----------------------------------------------------------------------------
// Rules that define the iteration history of the node position.

  // Default node position initial condition.
  class NodePositionInitialCondition : public pointwise_rule {
    private:
      const_store<vect3d> pos ;
      const_store<vect3d> s_ic ;
      store<vect3d> pos_ic ;
    public:

      // Define input and output.
      NodePositionInitialCondition() {
        name_store("pos",pos) ;
        name_store("node_s_ic",s_ic) ;
        name_store("pos_ic",pos_ic) ;
        input("pos,node_s_ic") ;
        output("pos_ic") ;
      }

      // Set position for single node.
      void calculate(Entity n) { pos_ic[n]=pos[n]+s_ic[n] ; }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<NodePositionInitialCondition>
    registerNodePositionInitialCondition ;

  // Time build rule for node position.
  class TimeBuildNodePositionBDF2 : public pointwise_rule {
    private:
      const_store<vect3d> pos_ic ;
      store<vect3d> pos ;
    public:

      // Define input and output.
      TimeBuildNodePositionBDF2() {
        name_store("pos_ic",pos_ic) ;
        name_store("pos{n=-1}",pos) ;
        input("pos_ic") ;
        output("pos{n=-1}") ;
      }

      // Set position for single node.
      void calculate(Entity n) { pos[n]=pos_ic[n] ; }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildNodePositionBDF2> registerTimeBuildNodePositionBDF2 ;

  // Time build rule for node position.
  class TimeBuildNodePosition : public pointwise_rule {
    private:
      const_store<vect3d> pos_ic ;
      store<vect3d> pos ;
    public:

      // Define input and output.
      TimeBuildNodePosition() {
        name_store("pos_ic",pos_ic) ;
        name_store("pos{n=0}",pos) ;
        input("pos_ic") ;
        output("pos{n=0}") ;
      }

      // Set position for single node.
      void calculate(Entity n) { pos[n]=pos_ic[n] ; }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildNodePosition> registerTimeBuildNodePosition ;

  // Iteration build rule for node position for time-dependent grid motion,
  // where the node positions have already been updated before going
  // into the "it" loop.
  class IterationBuildNodePositionGridMotionTimeDependent : public
  pointwise_rule {
    private:
      const_store<vect3d> posStar ;
      store<vect3d> pos ;
    public:

      // Define input and output.
      IterationBuildNodePositionGridMotionTimeDependent() {
        name_store("posStar{n}",posStar) ;
        name_store("pos{n,it=0}",pos) ;
        input("posStar{n}") ;
        output("pos{n,it=0}") ;
        constraint("pos,gridMotionTimeDependent") ;
      }

      // Assign position for a single node.
      void calculate(Entity n) { pos[n]=posStar[n] ; }

      // Loop over nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationBuildNodePositionGridMotionTimeDependent>
    registerIterationBuildNodePositionGridMotionTimeDependent ;

  // Iteration build rule for node position for solution-dependent grid motion,
  // where the node positions have not already been updated before going
  // into the "it" loop. Note that we use a different variable to track
  // the node position in the iterations loop here, called "pos_it". We do this
  // so that the updated position can be called "pos{n,it}", which would
  // not be possible if we used "pos" for the looping.
  class IterationBuildNodePositionGridMotionSolutionDependent : public
  pointwise_rule {
    private:
      const_store<vect3d> pos ;
      store<vect3d> pos_it ;
    public:

      // Define input and output.
      IterationBuildNodePositionGridMotionSolutionDependent() {
        name_store("pos{n}",pos) ;
        name_store("pos_it{n,it=0}",pos_it) ;
        input("pos{n}") ;
        output("pos_it{n,it=0}") ;
        constraint("pos,gridMotionSolutionDependent") ;
      }

      // Assign position for a single node.
      void calculate(Entity n) { pos_it[n]=pos[n] ; }

      // Loop over nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationBuildNodePositionGridMotionSolutionDependent>
    registerIterationBuildNodePositionGridMotionSolutionDependent ;

  // Iteration advance rule for node position for time-dependent grid
  // motion. Position stays constant during the "it" loop.
  class IterationAdvanceNodePositionGridMotionTimeDependent : public
  pointwise_rule {
    private:
      store<vect3d> pos ;
    public:

      // Define input and output.
      IterationAdvanceNodePositionGridMotionTimeDependent() {
        name_store("pos{n,it}",pos) ;
        input("pos{n,it}") ;
        output("pos{n,it+1}=pos{n,it}") ;
        constraint("pos,gridMotionTimeDependent") ;
      }

      // Empty compute method.
      void compute(const sequence &seq) {}
  } ;

  register_rule<IterationAdvanceNodePositionGridMotionTimeDependent>
    registerIterationAdvanceNodePositionGridMotionTimeDependent ;

  // Iteration advance rule for node position for solution-dependent grid
  // motion. Here, "pos{n,it}" is the updated position.
  class IterationAdvanceNodePositionGridMotionSolutionDependent : public
  pointwise_rule {
    private:
      store<vect3d> pos ;
    public:

      // Define input and output.
      IterationAdvanceNodePositionGridMotionSolutionDependent() {
        name_store("pos{n,it}",pos) ;
        input("pos{n,it}") ;
        output("pos_it{n,it+1}=pos{n,it}") ;
        constraint("pos,gridMotionSolutionDependent") ;
      }

      // Empty compute method.
      void compute(const sequence &seq) {}
  } ;

  register_rule<IterationAdvanceNodePositionGridMotionSolutionDependent>
    registerIterationAdvanceNodePositionGridMotionSolutionDependent ;

  // Iteration collapse rule for node position for time-dependent grid motion.
  class IterationCollapseNodePositionGridMotionTimeDependent : public
  pointwise_rule {
    private:
      store<vect3d> pos ;
    public:

      // Define input and output.
      IterationCollapseNodePositionGridMotionTimeDependent() {
        name_store("pos{n,it}",pos) ;
        input("pos{n,it}") ;
        output("pos{n+1}=pos{n,it}") ;
        conditional("iterationFinished{n,it-1}") ;
        constraint("pos,gridMotionTimeDependent") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<IterationCollapseNodePositionGridMotionTimeDependent>
    registerIterationCollapseNodePositionGridMotionTimeDependent ;

  // Iteration collapse rule for node position for solution-dependent
  // grid motion.
  class IterationCollapseNodePositionGridMotionSolutionDependent : public
  pointwise_rule {
    private:
      store<real> pos ;
    public:

      // Define input and output.
      IterationCollapseNodePositionGridMotionSolutionDependent() {
        name_store("pos_it{n,it}",pos) ;
        input("pos_it{n,it}") ;
        output("pos{n+1}=pos_it{n,it}") ;
        conditional("iterationFinished{n,it-1}") ;
        constraint("pos,gridMotionSolutionDependent") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<IterationCollapseNodePositionGridMotionSolutionDependent>
    registerIterationCollapseNodePositionGridMotionSolutionDependent ;

//-----------------------------------------------------------------------------
// Rules that define the iteration history of the node displacement.

  // Time build rule for node displacement.
  class TimeBuildNodeDisplacement : public pointwise_rule {
    private:
      const_store<vect3d> s_ic ;
      store<vect3d> s ;
    public:

      // Define input and output.
      TimeBuildNodeDisplacement() {
        name_store("node_s_ic",s_ic) ;
        name_store("node_s{n=0}",s) ;
        input("node_s_ic") ;
        output("node_s{n=0}") ;
      }

      // Set displacement for single node.
      void calculate(Entity n) { s[n]=s_ic[n] ; }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildNodeDisplacement> registerTimeBuildNodeDisplacement ;

  // Iteration build rule for node displacement when using the pde-based solver
  // for time-depdent grid motion.
  class IterationBuildNodeDisplacementGridMotionTimeDependent : public
  pointwise_rule {
    private:
      const_store<vect3d> s ;
      store<vect3d> sStar ;
    public:

      // Define input and output.
      IterationBuildNodeDisplacementGridMotionTimeDependent() {
        name_store("node_s{n}",s) ;
        name_store("node_sStar{n,itg=0}",sStar) ;
        input("node_s{n}") ;
        output("node_sStar{n,itg=0}") ;
        constraint("pos,gridMover,gridMotionTimeDependent") ;
      }

      // Assign displacement for a single node.
      void calculate(Entity n) { sStar[n]=s[n] ; }

      // Loop over nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationBuildNodeDisplacementGridMotionTimeDependent>
    registerIterationBuildNodeDisplacementGridMotionTimeDependent ;

  // Iteration advance rule for node displacement when using the pde-based
  // solver for time-depdent grid motion.
  class IterationAdvanceNodeDisplacementGridMotionTimeDependent : public
  pointwise_rule {
    private:
      store<vect3d> s ;
    public:

      // Define input and output.
      IterationAdvanceNodeDisplacementGridMotionTimeDependent() {
        name_store("sStar{n,itg}",s) ;
        input("sStar{n,itg}") ;
        output("node_sStar{n,itg+1}=sStar{n,itg}") ;
        constraint("pos,gridMover,gridMotionTimeDependent") ;
      }

      // Empty compute method.
      void compute(const sequence &seq) {}
  } ;

  register_rule<IterationAdvanceNodeDisplacementGridMotionTimeDependent>
    registerIterationAdvanceNodeDisplacementGridMotionTimeDependent ;

  // Iteration collapse rule for node displacement when using the pde-based
  // solver for time-depdent grid motion.
  class IterationCollapseNodeDisplacementGridMotionTimeDependent : public
  pointwise_rule {
    private:
      store<vect3d> s ;
    public:

      // Define input and output.
      IterationCollapseNodeDisplacementGridMotionTimeDependent() {
        name_store("node_sStar{n,itg}",s) ;
        input("node_sStar{n,itg}") ;
        output("node_sStar{n}=node_sStar{n,itg}") ;
        constraint("pos,gridMover,gridMotionTimeDependent") ;
        conditional("nodeIterationFinished{n,itg-1}") ;
      }

      // Empty compute method.
      void compute(const sequence &seq) {}
  } ;

  register_rule<IterationCollapseNodeDisplacementGridMotionTimeDependent>
    registerIterationCollapseNodeDisplacementGridMotionTimeDependent ;

  // Updated node positions for time-dependent grid motion.
  class UpdateNodePositionGridMotionTimeDependent : public pointwise_rule {
    private:
      const_store<vect3d> posRef ;
      const_store<vect3d> node_sStar ;
      store<vect3d> pos ;
    public:

      // Define input and output.
      UpdateNodePositionGridMotionTimeDependent() {
        name_store("pos",posRef) ;
        name_store("node_sStar{n}",node_sStar) ;
        name_store("posStar{n}",pos) ;
        input("pos,node_sStar{n}") ;
        output("posStar{n}") ;
        constraint("pos,gridMotionTimeDependent") ;
      }

      // Position for single node.
      void calculate(Entity n) { pos[n]=posRef[n]+node_sStar[n] ; }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<UpdateNodePositionGridMotionTimeDependent>
    registerUpdateNodePositionGridMotionTimeDependent ;

  // Iteration build rule for node displacement when using the pde-based solver
  // for solution-depdent grid motion.
  class IterationBuildNodeDisplacementGridMotionSolutionDependent : public
  pointwise_rule {
    private:
      const_store<vect3d> sN ;
      store<vect3d> s ;
    public:

      // Define input and output.
      IterationBuildNodeDisplacementGridMotionSolutionDependent() {
        name_store("node_s{n}",sN) ;
        name_store("node_s{n,it=0}",s) ;
        input("node_s{n}") ;
        output("node_s{n,it=0}") ;
        constraint("pos,gridMotionSolutionDependent") ;
      }

      // Assign displacement for a single node.
      void calculate(Entity n) { s[n]=sN[n] ; }

      // Loop over nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<IterationBuildNodeDisplacementGridMotionSolutionDependent>
    registerIterationBuildNodeDisplacementGridMotionSolutionDependent ;

  // Iteration advance rule for node displacement when using the pde-based
  // solver for solution-depdent grid motion.
  class IterationAdvanceNodeDisplacementGridMotionSolutionDependent : public
  pointwise_rule {
    private:
      store<vect3d> s ;
    public:

      // Define input and output.
      IterationAdvanceNodeDisplacementGridMotionSolutionDependent() {
        name_store("sStar{n,it}",s) ;
        input("sStar{n,it}") ;
        output("node_s{n,it+1}=sStar{n,it}") ;
        constraint("pos,gridMotionSolutionDependent") ;
      }

      // Empty compute method.
      void compute(const sequence &seq) {}
  } ;

  register_rule<IterationAdvanceNodeDisplacementGridMotionSolutionDependent>
    registerIterationAdvanceNodeDisplacementGridMotionSolutionDependent ;

  // Iteration collapse/time advance rule for node displacement when using
  // the pde-based solver for solution-depdent grid motion.
  class IterationCollapseNodeDisplacementGridMotionSolutionDependent : public
  pointwise_rule {
    private:
      store<vect3d> s ;
    public:

      // Define input and output.
      IterationCollapseNodeDisplacementGridMotionSolutionDependent() {
        name_store("node_s{n,it}",s) ;
        input("node_s{n,it}") ;
        output("node_s{n+1}=node_s{n,it}") ;
        conditional("iterationFinished{n,it-1}") ;
        constraint("pos,gridMotionSolutionDependent") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<IterationCollapseNodeDisplacementGridMotionSolutionDependent>
    registerIterationCollapseNodeDisplacementGridMotionSolutionDependent ;

  // Time advance rule for node displacement for time-dependent grid motion.
  class TimeAdvanceNodeDisplacementGridMotionTimeDependent : public
  pointwise_rule {
    private:
      store<vect3d> s ;
    public:
                                                                                
      // Define input and output.
      TimeAdvanceNodeDisplacementGridMotionTimeDependent() {
        name_store("node_sStar{n}",s) ;
        input("node_sStar{n}") ;
        output("node_s{n+1}=node_sStar{n}") ;
      }

      // Empty compute.
      virtual void compute(const sequence &seq) {}
  } ;
                                                                                
  register_rule<TimeAdvanceNodeDisplacementGridMotionTimeDependent>
    registerTimeAdvanceNodeDisplacementGridMotionTimeDependent ;

  // Time collapse rule for grid mover.
  class TimeCollapseGridMover : public pointwise_rule {
    private:
      const_store<vect3d> s ;
      store<int> solution ;
    public:

      // Define input and output.
      TimeCollapseGridMover() {
        name_store("node_s{n}",s) ;
        name_store("solution",solution) ;
        input("node_s{n}") ;
        output("solution") ;
        conditional("timeSteppingFinished{n}") ;
        constraint("nodes{n}") ;
      }

      // Empty compute method.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<TimeCollapseGridMover> registerTimeCollapseGridMover ;

//-----------------------------------------------------------------------------
// Rule to move grid with a rigid-body translation and rotation.

  // Compute the node displacement.
  class RigidBodyGridMotionDisplacement : public pointwise_rule {
    private:
      const_param<real> stime,timeStep ;
      const_param<vect3d> sDir,tDir,rCenter,rAxis ;
      const_param<real> sMag,tMag,tFreq,tPhi,rAlphaBar,rMag,rFreq,rPhi, func ;
      const_store<vect3d> pos ;
      store<vect3d> node_sStar ;
    private:
      real t,theta ;
      vect3d sTranslation ;
    public:

      // Define input and output.
      RigidBodyGridMotionDisplacement() {
        name_store("stime{n}",stime) ;
        name_store("timeStep{n}",timeStep) ;
        name_store("rigidBodyGridMotionTDir{n}",tDir) ;
        name_store("rigidBodyGridMotionRCenter{n}",rCenter) ;
        name_store("rigidBodyGridMotionRAxis{n}",rAxis) ;
        name_store("rigidBodyGridMotionTMag{n}",tMag) ;
        name_store("rigidBodyGridMotionTFreq{n}",tFreq) ;
        name_store("rigidBodyGridMotionTPhi{n}",tPhi) ;
        name_store("rigidBodyGridMotionRAlphaBar{n}",rAlphaBar) ;
        name_store("rigidBodyGridMotionRMag{n}",rMag) ;
        name_store("rigidBodyGridMotionRFreq{n}",rFreq) ;
        name_store("rigidBodyGridMotionRPhi{n}",rPhi) ;
        name_store("rigidBodyGridMotionSDir{n}",sDir) ;
        name_store("rigidBodyGridMotionSMag{n}",sMag) ;
        name_store("pos",pos) ;
        name_store("node_sStar{n}",node_sStar) ;
	name_store("rigidBodyGridMotionFunction{n}",func) ;
	input("rigidBodyGridMotionFunction{n}") ;
        input("stime{n},timeStep{n}") ;
        input("rigidBodyGridMotionTDir{n},rigidBodyGridMotionRCenter{n}") ;
        input("rigidBodyGridMotionRAxis{n},rigidBodyGridMotionTMag{n}") ;
        input("rigidBodyGridMotionTFreq{n},rigidBodyGridMotionRAlphaBar{n}") ;
        input("rigidBodyGridMotionRMag{n},rigidBodyGridMotionRFreq{n}") ;
        input("rigidBodyGridMotionRPhi{n},rigidBodyGridMotionTPhi{n}") ;
        input("rigidBodyGridMotionSDir{n},rigidBodyGridMotionSMag{n}") ;
        input("pos") ;
        output("node_sStar{n}") ;
      }

      // Set displacement for single node.
      void calculate(Entity n) {

        // Compute the rotation origin, which is the point on the axis of
        // rotation that when connected to this node is perpendicular to the
        // axis of rotation.
        vect3d r=pos[n]-(*rCenter) ; real s0=dot(r,*rAxis)/dot(*rAxis,*rAxis) ;
        vect3d rOrigin=(*rCenter)+s0*(*rAxis) ;

        // Compute normalized position vector from rotation center. Compute
        // vector perpendicular to postion vector and rotation axis. If a
        // node is on the axis, there is no rotational displacement.
        vect3d rr=pos[n]-rOrigin ; real R=norm(rr) ;
        if(R<1.0e-12){ node_sStar[n]=sTranslation ; return ; }
        vect3d rHat=(1.0/R)*rr ;
        vect3d tHat=cross(rHat,*rAxis) ;

        // Compute displacement.
        node_sStar[n]=sTranslation+R*((cos(theta)-1.0)*rHat-sin(theta)*tHat) ;
      }

      // Loop through nodes.
      void compute(const sequence &seq) {
        t=(*stime)+(*timeStep) ;
	
	const real PI2 = 2.0 * 4.0*atan(1.0) ;
	real rFreqNew = PI2 * (*rFreq) ;
	real tFreqNew = PI2 * (*tFreq) ;

	// Compute displacement.
	
	int functionType=0 ;
	functionType = int(*func) ;
	if (functionType == 1) {
	  // cout << "functionType == 1" << endl ;
	  sTranslation = (*sMag)*(*sDir)*t+((*tMag)*sin((*tFreq)*t+*tPhi))*(*tDir) ;
	} else if (functionType == 3) {
	  // cout << "functionType == 3" << endl ;
	  sTranslation = (*sMag)*(*sDir)*t+((*tMag)*(cos((*tFreq)*t+*tPhi)-1.0))*(*tDir) ;
	}
    
        //sTranslation=(*sMag)*(*sDir)*t+((*tMag)*sin((*tFreq)*t+*tPhi))*(*tDir) ;
        theta=(*rAlphaBar)+(*rMag)*sin((*rFreq)*t+(*rPhi)) ;
        if(Loci::MPI_rank==0) cout << "t,theta: " <<  t << " "
          << theta*180/3.1415926 << endl ;
        do_loop(seq,this) ;
      }
  } ;

  register_rule<RigidBodyGridMotionDisplacement>
    registerRigidBodyGridMotionDisplacement ;

//-----------------------------------------------------------------------------
// Generic rules independent of the grid motion type.

  // Compute old node velocities. This rule assumes that we have a constant
  // timestep.
  class OldNodeVelocity : public pointwise_rule {
    private:
      const_param<real> timeStep ;
      const_store<vect3d> posOld,pos ;
      store<vect3d> node_vOld ;
    public:

      // Define input and output.
      OldNodeVelocity() {
        name_store("timeStep{n}",timeStep) ;
        name_store("pos{n-1}",posOld) ;
        name_store("pos{n}",pos) ;
        name_store("node_vOld{n}",node_vOld) ;
        input("timeStep{n},pos{n-1},pos{n}") ;
        output("node_vOld{n}") ;
      }

      // Position for single node.
      void calculate(Entity n) {
        node_vOld[n]=(pos[n]-posOld[n])/(*timeStep) ;
      }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<OldNodeVelocity> registerOldNodeVelocity ;

  // Compute new node velocities for solution-dependent grid motion. Note
  // that pos{n,it} are the updated positions.
  class NodeVelocity : public pointwise_rule {
    private:
      const_param<real> timeStep ;
      const_store<vect3d> posOld ;
      const_store<vect3d> pos ;
      store<vect3d> node_v ;
    public:

      // Define input and output.
      NodeVelocity() {
        name_store("timeStep{n}",timeStep) ;
        name_store("pos{n}",posOld) ;
        name_store("pos{n,it}",pos) ;
        name_store("node_v{n,it}",node_v) ;
        input("timeStep{n},pos{n},pos{n,it}") ;
        output("node_v{n,it}") ;
        constraint("SIMPLE,pos") ;
      }

      // Position for single node.
      void calculate(Entity n) { node_v[n]=(pos[n]-posOld[n])/(*timeStep) ; }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<NodeVelocity> registerNodeVelocity ;

  // Default rule for face velocities.
  class GridFaceVelocityDefault : public pointwise_rule {
    private:
      store<vect3d> face_v ;
    public:

      // Define input and output.
      GridFaceVelocityDefault() {
        name_store("face_v",face_v) ;
        output("face_v") ;
        constraint("faces") ;
      }

      // Velocity for a single face.
      void calculate(Entity face) { face_v[face]=vect3d(0.0,0.0,0.0) ; }

      // Loop through faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GridFaceVelocityDefault> registerGridFaceVelocityDefault ;

  // Compute grid face velocity.
  class GridFaceVelocityPriority : public pointwise_rule {
    private:
      const_multiMap face2node ;
      const_store<vect3d> node_v ;
      store<vect3d> face_v ;
    public:

      // Define input and output.
      GridFaceVelocityPriority() {
        name_store("face2node",face2node) ;
        name_store("node_v",node_v) ;
        name_store("priority::face_v",face_v) ;
        input("face2node->node_v") ;
        output("priority::face_v") ;
      }

      // Velocity for a single face.
      void calculate(Entity face) {
        vect3d sum(0.0,0.0,0.0) ;
        for(const int *ni=face2node.begin(face);ni!=face2node.end(face);++ni)
          sum+=node_v[*ni] ;
        face_v[face]=(1.0/real(face2node.num_elems(face)))*sum ;
      }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GridFaceVelocityPriority> registerGridFaceVelocityPriority ;

//-----------------------------------------------------------------------------
// Rules for producing output data.

  // Magnitude of node displacement.
  class NodeDisplacementMagnitude : public pointwise_rule {
    private:
      const_store<vect3d> node_s ;
      store<float> sMag ;
    public:

      // Define input and output.
      NodeDisplacementMagnitude() {
        name_store("node_s",node_s) ;
        name_store("sMag",sMag) ;
        input("node_s") ;
        output("sMag") ;
        constraint("pos") ;
      }

      // Set displacement magnitude for single node.
      void calculate(Entity n) { sMag[n]=norm(node_s[n]) ; }

      // Loop through nodes.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<NodeDisplacementMagnitude> registerNodeDisplacementMagnitude ;
}
