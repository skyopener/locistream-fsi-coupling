//-----------------------------------------------------------------------------
// Description: This file contains rules for various test problems.
//
// Authors: Siddharth Thakur and Jeff Wright
//          Streamline Numerics Inc.
//          Gainesville, Florida 32609
//-----------------------------------------------------------------------------

// Standard library includes.
#include <string>
#include <sstream>
#include <stdlib.h>
#include <vector>
using std::string ;
using std::vector ;

// Loci includes.
#include <Loci.h>

// StreamUns includes.
#include <readGrid.h>
#include "sciTypes.h"

namespace streamUns {
//-----------------------------------------------------------------------------
// Testing.

  // Variable for use at {n}.
  class JunkVarTime : public pointwise_rule {
    private:
      store<real> junkVar ;
    public:

      // Define input and output.
      JunkVarTime() {
        name_store("junkVar{n}",junkVar) ;
        output("junkVar{n}") ;
        constraint("geom_cells{n}") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<JunkVarTime> registerJunkVarTime ;

  // Variable for use at {n,it}.
  class JunkVarIteration : public pointwise_rule {
    private:
      store<real> junkVar ;
    public:

      // Define input and output.
      JunkVarIteration() {
        name_store("priority::junkVar{n,it}",junkVar) ;
        output("priority::junkVar{n,it}") ;
        constraint("geom_cells{n,it}") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<JunkVarIteration> registerJunkVarIteration ;

//-----------------------------------------------------------------------------
// Flat plate boundary layer problem.

  // Rule to write out data normal to the flat plate at a certain location.
  class FlatPlateDataCell : public pointwise_rule {
    private:
      const_store<real> rho,laminarViscosity,k,omega,eddyViscosity ;
      const_store<vect3d> v ;
      const_store<vect3d> cellCenter ;
      param<bool> OUTPUT ;
    private:
      real xPosition,yPosition,zPosition,epsilon ;
    public:

      // Define input and output.
      FlatPlateDataCell() {
        name_store("rho{n,it}",rho) ;
        name_store("muu(temperature,p,y){n,it}",laminarViscosity) ;
        name_store("k{n,it}",k) ;
        name_store("omega{n,it}",omega) ;
        name_store("eddyViscosity{n,it}",eddyViscosity) ;
        name_store("v{n,it}",v) ;
        name_store("cellcenter{n,it}",cellCenter) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("rho{n,it},muu(temperature,p,y){n,it},k{n,it},omega{n,it}") ;
        input("eddyViscosity{n,it},v{n,it},cellcenter{n,it}") ;
        output("OUTPUT{n,it}") ;
        conditional("do_plot{n,it}") ;
        constraint("geom_cells{n,it},flatPlate{n,it}") ;
      }

      // Write data for the cell if it is at the desired location.
      void calculate(Entity cell) {
        if(abs(cellCenter[cell].x-xPosition)<epsilon && abs(cellCenter[cell].z-
        zPosition)<epsilon){
          ofstream out("flatPlate.txt",ios::app) ;
          out.setf(ios::scientific,ios::floatfield) ; out.precision(6) ;
          out << cellCenter[cell].y << " " << rho[cell] << " " << v[cell].x
            << " " << laminarViscosity[cell] << " " << k[cell] << " "
            << omega[cell] << " " << eddyViscosity[cell] << endl ;
        }
      }

      // Write the data.
      virtual void compute(const sequence &seq) {
        xPosition=0.95 ; zPosition=2.5e-04 ; epsilon=1.0e-07 ;
        cout << "Writing flatPlate data." << endl ;
        do_loop(seq,this) ;
      }
  } ;

  register_rule<FlatPlateDataCell> registerFlatPlateDataCell ;

//-----------------------------------------------------------------------------
// Cavity problem.

  // Rule to write out the x and y velocity profiles.
  class CavityData : public pointwise_rule {
    private:
      const_param<int> n ;
      const_param<real> dt ;
      const_Map cl,cr ;
      const_store<vect3d> v ;
      const_store<vect3d> faceCenter ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      CavityData() {
        name_store("$n{n}",n) ;
        name_store("dt{n}",dt) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("v{n,it}",v) ;
        name_store("facecenter{n,it}",faceCenter) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("$n{n},dt{n}") ;
        input("(cl,cr)->v{n,it},facecenter{n,it}") ;
        output("OUTPUT{n,it}") ;
        conditional("do_print{n,it}") ;
        constraint("internalFaces{n,it},cavity{n,it}") ;
      }

      // Write the data.
      virtual void compute(const sequence &seq) {
        ostringstream ossVX,ossVY ;
        ossVX << "output/vX_" << *n << ".dat" ;
        ossVY << "output/vY_" << *n << ".dat" ;
        string vXFile=ossVX.str(),vYFile=ossVY.str() ;
        ofstream outVX(vXFile.c_str(),ios::out),outVY(vYFile.c_str(),ios::out) ;
        outVX.setf(ios::scientific,ios::floatfield) ; outVX.precision(6) ;
        outVY.setf(ios::scientific,ios::floatfield) ; outVY.precision(6) ;
        cout << "Writing cavity data." << endl ;
        for(sequence::const_iterator si=seq.begin();si!=seq.end();++si) {
          int face=*si ;
          real xFace=faceCenter[face].x,yFace=faceCenter[face].y ;
          if(xFace>0.49999 && xFace < 0.50001)
            outVX << yFace << " " << 0.5*(v[cl[face]].x+v[cr[face]].x) << endl ;
          if(yFace>0.49999 && yFace < 0.50001)
            outVY << xFace << " " << 0.5*(v[cl[face]].y+v[cr[face]].y) << endl ;
        }
      }
  } ;

  register_rule<CavityData> registerCavityData ;

//-----------------------------------------------------------------------------
// Unsteady cylinder problem.

  // Rule to write out the recirculation length behind a cylinder at Re=40.
  class CylinderData : public pointwise_rule {
    private:
      const_param<real> sTime ;
      const_Map cl,cr ;
      const_store<vect3d> v ;
      const_store<vect3d> cellCenter ;
      param<bool> OUTPUT ;
    private:
      real xMin,xMax,vXMin,vXMax ;
    public:

      // Define input and output.
      CylinderData() {
        name_store("stime{n}",sTime) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("v{n,it}",v) ;
        name_store("cellcenter{n,it}",cellCenter) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("stime{n},(cl,cr)->v{n,it},(cl,cr)->cellcenter{n,it}") ;
        output("OUTPUT{n,it}") ;
        conditional("do_print{n,it}") ;
        constraint("internalFaces{n,it},cylinder{n,it}") ;
      }

      // Write data for the cell if it is at the desired location.
      void calculate(Entity face) {
        real epsilon=1.0e-06 ;
        real xFace=0.5*(cellCenter[cl[face]].x+cellCenter[cr[face]].x) ;
        real yFace=0.5*(cellCenter[cl[face]].y+cellCenter[cr[face]].y) ;
        if(abs(yFace)<epsilon){
          if(xFace>0.5){
            real vX=0.5*(v[cl[face]].x+v[cr[face]].x) ;
            if(vX<0.0 && xFace>xMin){ xMin=xFace ; vXMin=vX ; }
            if(vX>0.0 && xFace<xMax){ xMax=xFace ; vXMax=vX ; }
          }
        }
      }

      // Write the data.
      virtual void compute(const sequence &seq) {
        xMin=0.5 ; xMax=20.0 ; do_loop(seq,this) ;
        real x=(-vXMin*xMax+vXMax*xMin)/(vXMax-vXMin) ;
        ofstream out("recirculationLength.txt",ios::app) ;
        out.setf(ios::scientific,ios::floatfield) ; out.precision(6) ;
        out << *sTime << " " << x-0.5 << endl ;
      }
  } ;

  register_rule<CylinderData> registerCylinderData ;

//-----------------------------------------------------------------------------
// Laminar diffusion problem.

  // Rule to hard-code the intial condition for the injector case.
  class DiffusionInitialCondition : public pointwise_rule {
    private:
      const_param<int> numSpecies ;
      const_store<vect3d> cellCenter ;
      storeVec<real> y_ic;
    public:

      // Define input and output.
      DiffusionInitialCondition() {
        name_store("numSpecies",numSpecies) ;
        name_store("cellcenter",cellCenter) ;
        name_store("hardcode::y_ic",y_ic) ;
        input("numSpecies,cellcenter") ;
        output("hardcode::y_ic") ;
        constraint("geom_cells,diffusion,noRestart") ;
      }

      // Assign initial condition based on location.
      void calculate(Entity cell) {
        y_ic[cell][0]=1.0-pow(cellCenter[cell].x,3.0) ;
        y_ic[cell][1]=pow(cellCenter[cell].x,3.0) ;
      }

      // Set initial condition for all cells.
      virtual void compute(const sequence &seq) {
        y_ic.setVecSize(*numSpecies) ; do_loop(seq,this) ;
      }
  } ;

  register_rule<DiffusionInitialCondition> registerDiffusionInitialCondition ;

//-----------------------------------------------------------------------------
// Channel combustion problem benchmark.

  // Rule to write out data normal to the flat plate at a certain location.
  class ChannelLinePlot : public pointwise_rule {
    private:
      const_store<vect3d> cellCenter ;
      const_store<real> temperature ;
      const_storeVec<real> y ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      ChannelLinePlot() {
        name_store("cellcenter{n,it}",cellCenter) ;
        name_store("temperature{n,it}",temperature) ;
        name_store("y{n,it}",y) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("cellcenter{n,it},temperature{n,it},y{n,it}") ;
        output("OUTPUT{n,it}") ;
        conditional("do_plot{n,it}") ;
        constraint("geom_cells{n,it},channel{n,it}") ;
      }

      // Write data for the cell if it is at the desired location.
      void calculate(Entity cell) {
        if(cellCenter[cell].x>0.0 && cellCenter[cell].y>0.006 && cellCenter
        [cell].y<0.012){
          ofstream out("channel.dat",ios::app) ;
          out.setf(ios::scientific,ios::floatfield) ; out.precision(6) ;
          out << cellCenter[cell].x << " " << temperature[cell] << " "
            << y[cell][0] << " " << y[cell][1] << " " << y[cell][2] << " "
            << y[cell][3] << " " << y[cell][4] << " " << y[cell][5] << endl ;
        }
      }

      // Write the data.
      virtual void compute(const sequence &seq) {
        cout << "Writing channel data." << endl ; do_loop(seq,this) ;
      }
  } ;

  register_rule<ChannelLinePlot> registerChannelLinePlot ;

//-----------------------------------------------------------------------------
// Dumping some things here that we need temporarily in going to Loci-3-1 .

  class cell2nodeMaxMag_v3d : public unit_rule {
    private:
      store<vector3d<float> > node_value ;
    public:
                                                                                
      // Define input and output.
      cell2nodeMaxMag_v3d() {
        name_store("cell2nodeMaxMag_v3d(X)",node_value) ;
        constraint("pos") ;
        output("cell2nodeMaxMag_v3d(X)") ;
      }
                                                                                
      // Initialize the value.
      void calculate(Entity e) { node_value[e]=vector3d<float>(0.0,0.0,0.0) ; }
                                                                                
      // Loop over all entities.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<cell2nodeMaxMag_v3d> register_cell2nodeMaxMag_v3d ;

  struct max_mag_v3d_join {
    void operator()(vector3d<float> &f1, const vector3d<float> &f2) {
      if(fabs(f2.x) > fabs(f1.x)) f1.x=f2.x ;
      if(fabs(f2.y) > fabs(f1.y)) f1.y=f2.y ;
      if(fabs(f2.z) > fabs(f1.z)) f1.z=f2.z ;
    }
  } ;
                                                                                
  class c2nMaxMag_v3d : public apply_rule<store<vector3d<float> >,
  max_mag_v3d_join > {
    const_multiMap upper, lower, boundary_map ;
    const_multiMap face2node ;
    const_store<vector3d<real> > X ;
    store<vector3d<float> > cell2nodeMaxMag_v3d ;
    vector<int> node_list ;
  public:
    c2nMaxMag_v3d() ;
    void calculate(Entity cc) ;
    virtual void compute(const sequence &seq) ;
  } ;
                                                                                
  c2nMaxMag_v3d::c2nMaxMag_v3d() {
    name_store("cell2nodeMaxMag_v3d(X)",cell2nodeMaxMag_v3d) ;
    name_store("face2node",face2node) ;
    name_store("upper",upper) ;
    name_store("lower",lower) ;
    name_store("boundary_map",boundary_map) ;
    name_store("X",X) ;
                                                                                
    constraint("geom_cells") ;
    input("X") ;
    input("(upper,lower,boundary_map)->face2node->cell2nodeMaxMag_v3d(X)") ;
    output("(upper,lower,boundary_map)->face2node->cell2nodeMaxMag_v3d(X)") ;
  }
                                                                                
  void c2nMaxMag_v3d::calculate(Entity cc) {
                                                                                
    node_list.clear() ;
    for(const Entity *fi=upper.begin(cc);fi!=upper.end(cc);++fi)
      for(const Entity *ni=face2node.begin(*fi);ni!=face2node.end(*fi);++ni)
        node_list.push_back(*ni) ;
    for(const Entity *fi=lower.begin(cc);fi!=lower.end(cc);++fi)
      for(const Entity *ni=face2node.begin(*fi);ni!=face2node.end(*fi);++ni)
        node_list.push_back(*ni) ;
    for(const Entity *fi=boundary_map.begin(cc);fi!=boundary_map.end(cc);++fi)
      for(const Entity *ni=face2node.begin(*fi);ni!=face2node.end(*fi);++ni)
        node_list.push_back(*ni) ;
    sort(node_list.begin(),node_list.end()) ;
    vector<int>::iterator ns = node_list.begin() ;
    vector<int>::iterator ne = unique(node_list.begin(),node_list.end()) ;
                                                                                
    for(vector<int>::iterator vi = ns;vi!=ne;++vi) {
      int nd = *vi ;
      join(cell2nodeMaxMag_v3d[nd],vector3d<float>(X[cc].x,X[cc].y,X[cc].z)) ;
    }
                                                                                
  }
                                                                                
  void c2nMaxMag_v3d::compute(const sequence &seq) { do_loop(seq,this) ; }
                                                                                
  register_rule<c2nMaxMag_v3d> register_c2nMaxMag_v3d ;

//-----------------------------------------------------------------------------
// Rules to test nodal gradient rule.

  class NodeFunction : public pointwise_rule {
    private:
      const_store<vect3d> pos ;
      store<vect3d> nodeFunction ;
    public:

      // Define input and output.
      NodeFunction() {
        name_store("pos",pos) ;
        name_store("nodeFunction",nodeFunction) ;
        input("pos") ;
        output("nodeFunction") ;
      }

      // Write data for the cell if it is at the desired location.
      void calculate(Entity node) {
        nodeFunction[node].x=1.0*pos[node].x+2.0*pos[node].y+3.0*pos[node].z ;
        nodeFunction[node].y=4.0*pos[node].x+5.0*pos[node].y+6.0*pos[node].z ;
        nodeFunction[node].z=7.0*pos[node].x+8.0*pos[node].y+9.0*pos[node].z ;
      }

      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<NodeFunction> registerNodeFunction ;

  class NodeGradientTest : public pointwise_rule {
    private:
      const_store<tens3d> grad ;
      store<bool> nodeGradientTest ;
    public:

      // Define input and output.
      NodeGradientTest() {
        name_store("nodeGrad(nodeFunction)",grad) ;
        name_store("nodeGradientTest",nodeGradientTest) ;
        input("nodeGrad(nodeFunction)") ;
        output("nodeGradientTest") ;
        constraint("pos") ;
      }

      // Write data for the node.
      void calculate(Entity node) {
        cout << "node,grad: " << node << " " << grad[node].x.x << " "
          << grad[node].x.y << " " << grad[node].x.z << " " << grad[node].y.x
          << " " << grad[node].y.y << " " << grad[node].y.z << " "
          << grad[node].z.x << " " << grad[node].z.y << " " << grad[node].z.z
          << endl ;
      }

      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<NodeGradientTest> registerNodeGradientTest ;

//-----------------------------------------------------------------------------
// Rules to test new cell volume rule for deforming meshes.

  class CellVolumeTest : public pointwise_rule {
    private:
      const_store<real> vol ;
      store<bool> cellVolumeTest ;
    public:

      // Define input and output.
      CellVolumeTest() {
        name_store("vol",vol) ;
        name_store("cellVolumeTest",cellVolumeTest) ;
        input("vol") ;
        output("cellVolumeTest") ;
        constraint("geom_cells") ;
      }

      // Write data for the node.
      void calculate(Entity cell) {
        cout << "cell,vol: " << cell << " " << vol[cell] << endl ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<CellVolumeTest> registerCellVolumeTest ;

//-----------------------------------------------------------------------------
// Taylor decaying vortex problem.

  // This rule creates an "exact" boundary condition type with the allowable
  // options "taylor" and "rho".
  class CheckExactTaylor : public BC_Check {
    private:
      string errorMessage ;
    public:
      CheckExactTaylor() : errorMessage("") {}
    public:
      string BoundaryConditions() { return "exact" ; }
      bool CheckOptions(const options_list& bc_options,fact_db &facts) {
        return true ;
      }
      ostream &ErrorMessage(std::ostream &s) {
        s << errorMessage << endl ; return s;
      }
      string VariablesChecked(fact_db &facts) {
        string s="taylor,rho" ; return s ;
      }
  } ;

  register_BC<CheckExactTaylor> registerCheckExactTaylor ;

  // Add additional constraints for the exact boundary condition.
  class ExactConstraints : public pointwise_rule {
    private:
      store<bool> extrapolatedPressure_BC ;
    public:

      // Define input and output.
      ExactConstraints() {
        name_store("extrapolatedPressure_BC",extrapolatedPressure_BC) ;
        output("extrapolatedPressure_BC") ;
        constraint("exact_BC") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<ExactConstraints> registerExactConstraints ;

  // Sets the default Reynolds number.
  class DefaultReynoldsNumber : public default_rule {
    private:
      param<real> Re ;
    public:

      // Define input and output.
      DefaultReynoldsNumber() {
        name_store("Re",Re) ;
        output("Re") ;
        comments("Sets the default Reynolds number to 100.0.") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) { *Re=100.0 ; }
  } ;

  register_rule<DefaultReynoldsNumber> registerDefaultReynoldsNumber ;

  // Priority rule to assign exact initial condition.
  class InitialConditionTaylor : public pointwise_rule {
    private:
      const_param<real> Re ;
      const_store<vect3d> cellCenter ;
      store<vect3d> v_ic ;
      store<real> rho_ic,p_ic ;
    public:

      // Define input and output.
      InitialConditionTaylor() {
        name_store("Re",Re) ;
        name_store("cellcenter",cellCenter) ;
        name_store("taylor::rho_ic",rho_ic) ;
        name_store("taylor::v_ic",v_ic) ;
        name_store("taylor::p_ic",p_ic) ;
        input("Re,cellcenter") ;
        constraint("taylor,geom_cells") ;
        output("taylor::rho_ic,taylor::v_ic,taylor::p_ic") ;
      }

      // Set velocity and pressure for a single cell.
      void calculate(Entity cell) {
        real x=cellCenter[cell].x,y=cellCenter[cell].y ;
        rho_ic[cell]=1.0 ;
        v_ic[cell]=vect3d(-cos(x)*sin(y),sin(x)*cos(y),0.0) ;
        p_ic[cell]=-0.25*(cos(2.0*x)+cos(2.0*y)) ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<InitialConditionTaylor> registerInitialConditionTaylor ;

  // Priority rule to assign exact solution for velocity at the boundary.
  class BoundarySpecificationTaylor : public pointwise_rule {
    private:
      const_param<int> it ;
      const_param<real> oldSolutionTime ;
      const_param<real> dt ;
      const_param<real> Re ;
      const_store<vect3d> faceCenter ;
      store<vect3d> v_f ;
      store<real> p_f ;
    public:

      // Define input and output.
      BoundarySpecificationTaylor() {
        name_store("$it{n,it}",it) ;
        name_store("Re{n,it}",Re) ;
        name_store("stime{n}",oldSolutionTime) ;
        name_store("dt{n}",dt) ;
        name_store("facecenter{n,it}",faceCenter) ;
        name_store("taylor::v_f{n,it}",v_f) ;
//      name_store("taylor::p_f{n,it}",p_f) ;
        input("$it{n,it}") ;
        input("Re{n,it},stime{n},dt{n},facecenter{n,it}") ;
        constraint("ref->taylor_BCoption{n,it}") ;
//      output("taylor::v_f{n,it},taylor::p_f{n,it}") ;
        output("taylor::v_f{n,it}") ;
      }

      // Set velocity and pressure for a single face.
      void calculate(Entity face) {
        real x=faceCenter[face].x,y=faceCenter[face].y ;
        real time=*oldSolutionTime+*dt ;
        v_f[face]=vect3d(-cos(x)*sin(y)*exp(-2.0*time/(*Re)),sin(x)*cos(y)*
          exp(-2.0*time/(*Re)),0.0) ;
//      p_f[face]=-0.25*(cos(2.0*x)+cos(2.0*y))*exp(-4.0*time/(*Re)) ;
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

  register_rule<BoundarySpecificationTaylor>
    registerBoundarySpecificationTaylor ;

  // Computes the error in each cell.
  class ErrorTaylor : public pointwise_rule {
    private:
      const_param<real> Re ;
      const_param<real> dt ;
      const_param<real> oldSolutionTime ;
      const_store<vect3d> cellCenter ;
      const_store<vect3d> v ;
      const_store<real> p ;
      store<real> vXError,vYError,pError ;
    public:

      // Define input and output.
      ErrorTaylor() {
        name_store("Re{n,it}",Re) ;
        name_store("dt{n}",dt) ;
        name_store("stime{n}",oldSolutionTime) ;
        name_store("cellcenter{n,it}",cellCenter) ;
        name_store("v{n,it}",v) ;
        name_store("p{n,it}",p) ;
        name_store("vXError{n,it}",vXError) ;
        name_store("vYError{n,it}",vYError) ;
        name_store("pError{n,it}",pError) ;
        input("Re{n,it},dt{n},stime{n},cellcenter{n,it}") ;
        input("v{n,it},p{n,it}") ;
        output("vXError{n,it},vYError{n,it},pError{n,it}") ;
        constraint("taylor,geom_cells") ;
      }

      // Compute error for each cell.
      void calculate(Entity cell) {
        real x=cellCenter[cell].x,y=cellCenter[cell].y ;
        real time=*oldSolutionTime+*dt ;
        vXError[cell]=v[cell].x+cos(x)*sin(y)*exp(-2.0*time/(*Re)) ;
        vYError[cell]=v[cell].y-sin(x)*cos(y)*exp(-2.0*time/(*Re)) ;
        pError[cell]=p[cell]+0.25*(cos(2.0*x)+cos(2.0*y))*
          exp(-4.0*time/(*Re)) ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }

  } ;

//register_rule<ErrorTaylor> registerErrorTaylor ;

  // Rule to write out L1 norm error for the various variables.
  class WriteNormTaylor : public singleton_rule {
    private:
      const_param<int> n ;
      const_param<real> l1NormVX,l1NormVY,l1NormP ;
      const_param<real> l2NormVX,l2NormVY,l2NormP ;
      const_param<real> lInfNormVX,lInfNormVY,lInfNormP ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      WriteNormTaylor() {
        name_store("$n{n}",n) ;
        name_store("L1Norm(vXError){n,it}",l1NormVX) ;
        name_store("L1Norm(vYError){n,it}",l1NormVY) ;
        name_store("L1Norm(pError){n,it}",l1NormP) ;
        name_store("L2Norm(vXError){n,it}",l2NormVX) ;
        name_store("L2Norm(vYError){n,it}",l2NormVY) ;
        name_store("L2Norm(pError){n,it}",l2NormP) ;
        name_store("LinfNorm(vXError){n,it}",lInfNormVX) ;
        name_store("LinfNorm(vYError){n,it}",lInfNormVY) ;
        name_store("LinfNorm(pError){n,it}",lInfNormP) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("$n{n}") ;
        input("L1Norm(vXError){n,it},L1Norm(vYError){n,it}") ;
        input("L1Norm(pError){n,it}") ;
        input("L2Norm(vXError){n,it},L2Norm(vYError){n,it}") ;
        input("L2Norm(pError){n,it}") ;
        input("LinfNorm(vXError){n,it},LinfNorm(vYError){n,it}") ;
        input("LinfNorm(pError){n,it}") ;
        output("OUTPUT{n,it}") ;
//      conditional("do_plot{n,it}") ;
      }

      // Write to output file.
      void compute(const sequence &seq) {
        if(Loci::MPI_rank != 0) return ;
        ostringstream oss ; oss << "output/norm_" << *n << ".dat" ;
        string file=oss.str() ;
        ofstream ofile(file.c_str(),ios::out) ; ofile.precision(16) ;
        ofile << *l1NormVX << ' ' << *l1NormVY << ' ' << *l1NormP << endl ;
        ofile << *l2NormVX << ' ' << *l2NormVY << ' ' << *l2NormP << endl ;
        ofile << *lInfNormVX << ' ' << *lInfNormVY << ' ' << *lInfNormP
          << endl ;
      }
  } ;

//register_rule<WriteNormTaylor> registerWriteNormTaylor ;

  // Rule to write out the x and y velocity profiles.
  class TaylorData : public pointwise_rule {
    private:
      const_param<int> n ;
      const_param<real> dt ;
      const_param<real> solutionTime ;
      const_param<real> Re ;
      const_Map cl,cr ;
      const_store<vect3d> v ;
      const_store<real> p ;
      const_store<vect3d> faceCenter ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      TaylorData() {
        name_store("$n{n}",n) ;
        name_store("dt{n}",dt) ;
        name_store("stime{n}",solutionTime) ;
        name_store("Re{n,it}",Re) ;
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("v{n,it}",v) ;
        name_store("p{n,it}",p) ;
        name_store("facecenter{n,it}",faceCenter) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("$n{n},dt{n},stime{n},Re{n,it}") ;
        input("(cl,cr)->(v{n,it},p{n,it}),facecenter{n,it}") ;
        output("OUTPUT{n,it}") ;
        conditional("do_plot{n,it}") ;
        constraint("internalFaces{n,it},taylor{n,it}") ;
      }

      // Write the data.
      virtual void compute(const sequence &seq) {

        // Write out numerical solution.
        ostringstream ossVX,ossVY,ossP ;
        ossVX << "output/vX_" << *n+1 << ".dat" ;
        ossVY << "output/vY_" << *n+1 << ".dat" ;
        ossP << "output/p_" << *n+1 << ".dat" ;
        string vXFile=ossVX.str(),vYFile=ossVY.str(),pFile=ossP.str() ;
        ofstream outVX(vXFile.c_str(),ios::out),outVY(vYFile.c_str(),ios::out),
          outP(pFile.c_str(),ios::out) ;
        outVX.setf(ios::scientific,ios::floatfield) ; outVX.precision(12) ;
        outVY.setf(ios::scientific,ios::floatfield) ; outVY.precision(12) ;
        outP.setf(ios::scientific,ios::floatfield) ; outP.precision(12) ;
        cout << "Writing numerical taylor solution." << endl ;
        for(sequence::const_iterator si=seq.begin();si!=seq.end();++si) {
          int face=*si ;
          real xFace=faceCenter[face].x,yFace=faceCenter[face].y ;
          if(xFace>3.1415 && xFace < 3.1416){
            outVX << yFace << " " << 0.5*(v[cl[face]].x+v[cr[face]].x) << endl ;
          }
        }
        for(sequence::const_iterator si=seq.begin();si!=seq.end();++si) {
          int face=*si ;
          real xFace=faceCenter[face].x,yFace=faceCenter[face].y ;
          if(yFace>3.1415 && yFace < 3.1416){
            outVY << xFace << " " << 0.5*(v[cl[face]].y+v[cr[face]].y) << endl ;
            outP << xFace << " " << 0.5*(p[cl[face]]+p[cr[face]]) << endl ;
          }
        }

        // Write out exact solution.
        ostringstream ossVXExact,ossVYExact,ossPExact ;
        ossVXExact << "output/vXExact_" << *n+1 << ".dat" ;
        ossVYExact << "output/vYExact_" << *n+1 << ".dat" ;
        ossPExact << "output/pExact_" << *n+1 << ".dat" ;
        string vXExactFile=ossVXExact.str(),vYExactFile=ossVYExact.str(),
          pExactFile=ossPExact.str() ;
        ofstream outVXExact(vXExactFile.c_str(),ios::out),
          outVYExact(vYExactFile.c_str(),ios::out),
          outPExact(pExactFile.c_str(),ios::out) ;
        outVXExact.setf(ios::scientific,ios::floatfield) ;
        outVXExact.precision(12) ;
        outVYExact.setf(ios::scientific,ios::floatfield) ;
        outVYExact.precision(12) ;
        outPExact.setf(ios::scientific,ios::floatfield) ;
        outPExact.precision(12) ;
        cout << "Writing exact taylor solution." << endl ;
        real t=*solutionTime+*dt ;
        for(sequence::const_iterator si=seq.begin();si!=seq.end();++si) {
          int face=*si ;
          real xFace=faceCenter[face].x,yFace=faceCenter[face].y ;
          if(xFace>3.1415 && xFace < 3.1416){
            outVXExact << yFace << " " << -cos(xFace)*sin(yFace)*
              exp(-2.0*t/(*Re)) << endl ;
          }
        }
        for(sequence::const_iterator si=seq.begin();si!=seq.end();++si) {
          int face=*si ;
          real xFace=faceCenter[face].x,yFace=faceCenter[face].y ;
          if(yFace>3.1415 && yFace < 3.1416){
            outVYExact << xFace << " " << sin(xFace)*cos(yFace)*
              exp(-2.0*t/(*Re)) << endl ;
            outP << xFace << " " << -0.25*(cos(2.0*xFace)+cos(2.0*yFace))*
              exp(-4.0*t/(*Re)) << endl ;
          }
        }
      }
  } ;

//register_rule<TaylorData> registerTaylorData ;

  // Rule to write out the total kinetic energy in the domain.
  class TaylorKineticEnergy : public pointwise_rule {
    private:
      const_param<real> dt ;
      const_param<real> solutionTime ;
      const_param<real> Re ;
      const_store<vect3d> v ;
      const_store<real> vol ;
      const_store<vect3d> cellCenter ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      TaylorKineticEnergy() {
        name_store("dt{n}",dt) ;
        name_store("stime{n}",solutionTime) ;
        name_store("Re{n,it}",Re) ;
        name_store("v{n,it}",v) ;
        name_store("vol{n,it}",vol) ;
        name_store("cellcenter{n,it}",cellCenter) ;
        name_store("OUTPUT{n,it}",OUTPUT) ;
        input("dt{n},stime{n},Re{n,it}") ;
        input("v{n,it},vol{n,it},cellcenter{n,it}") ;
        output("OUTPUT{n,it}") ;
        conditional("do_print{n,it}") ;
        constraint("geom_cells{n,it},taylor{n,it}") ;
      }

      // Write the data.
      virtual void compute(const sequence &seq) {

        // Write out numerical solution.
        ostringstream ossKE,ossKEExact ;
        ossKE << "output/KE.dat" ; ossKEExact << "output/KEExact.dat" ;
        string kEFile=ossKE.str(),kEExactFile=ossKEExact.str() ;
        ofstream outKE(kEFile.c_str(),ios::app) ;
        ofstream outKEExact(kEExactFile.c_str(),ios::app) ;
        outKE.setf(ios::scientific,ios::floatfield) ; outKE.precision(12) ;
        outKEExact.setf(ios::scientific,ios::floatfield) ;
        outKEExact.precision(12) ;
        cout << "Writing taylor kinetic energy." << endl ;
        real kineticEnergy=0.0 ;
        real t=*solutionTime+*dt ;
        real kineticEnergyExact=acos(-1.0)*acos(-1.0)*exp(-4.0*t/(*Re)) ;
        for(sequence::const_iterator si=seq.begin();si!=seq.end();++si) {
          Entity cell=*si ;
          kineticEnergy+=0.5*vol[cell]*dot(v[cell],v[cell]) ;
        }
        outKE << t << " " << kineticEnergy << endl ;
        outKEExact << t << " " << kineticEnergyExact << endl ;
      }
  } ;

  register_rule<TaylorKineticEnergy> registerTaylorKineticEnergy ;

}
