//-----------------------------------------------------------------------------
// Rules to implement various forms of Detached-Eddy Simulation turbulence
// models.
//-----------------------------------------------------------------------------

// Standard library includes.
#include <vector>
using std::vector ;

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// StreamUns includes.
#include "sciTypes.h"
#include "varsFileInputs.h"

namespace streamUns {

  // Creates the constraints which activate the DES turbulence models.
  class DESConstraints : public constraint_rule {
    private:
      const_param<TurbulenceEquationOptions> turbulenceEquationOptions ;
      Constraint DES2001TurbulenceModel,DDESTurbulenceModel ;
    public:

      // Define input and output.
      DESConstraints() {
        name_store("turbulenceEquationOptions",turbulenceEquationOptions) ;
        name_store("DES2001TurbulenceModel",DES2001TurbulenceModel) ;
        name_store("DDESTurbulenceModel",DDESTurbulenceModel) ;
        input("turbulenceEquationOptions") ;
        output("DES2001TurbulenceModel,DDESTurbulenceModel") ;
      }

      // Set up the constraints.
      virtual void compute(const sequence& seq) {
        if((*turbulenceEquationOptions).optionExists("des")){
          Loci::option_value_type optionValueType=turbulenceEquationOptions->
            getOptionValueType("des") ;
          switch(optionValueType){
            case Loci::NAME:
              {
                Loci::option_values optionValues=turbulenceEquationOptions->
                  getOption("des") ;
                string name ; optionValues.get_value(name) ;
                if(name=="DES2001"){
                  DES2001TurbulenceModel=~EMPTY ;
                  DDESTurbulenceModel=EMPTY ;
                }else if(name=="DDES"){
                  DES2001TurbulenceModel=EMPTY ;
                  DDESTurbulenceModel=~EMPTY ;
                }else{
                  cerr << "Bad DES model for turbulenceEquation." << endl ;
                  Loci::Abort() ;
                }
              }
              break ;
            default:
              cerr << "Bad type for DES model in turbulenceEquation."
                << endl ; Loci::Abort() ;
          }
        }else{
          DES2001TurbulenceModel=EMPTY ; DDESTurbulenceModel=EMPTY ;
        }
      }
  } ;

  register_rule<DESConstraints> registerDESConstraints ;

  // Rule to initialize the maximum edge length for a cell. For two-dimensional
  // simulations there is an implicit assumption that the relevant coordinates
  // are x and y. Edge lengths then do not involve the z-coordinate.
  class MaxEdgeLengthUnit : public unit_rule {
    private:
      store<real> maxEdgeLength ;
    public:

      // Define input and output.
      MaxEdgeLengthUnit() {
        name_store("maxEdgeLength",maxEdgeLength) ;
        output("maxEdgeLength") ;
        constraint("geom_cells") ;
      }

      // Initialize to zero.
      void calculate(Entity cell) { maxEdgeLength[cell]=0.0 ; }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<MaxEdgeLengthUnit> registerMaxEdgeLengthUnit ;

  // Add contributions from interior faces. This rule is kind of brute
  // force and does twice as many calculations as necessary, but there
  // is no way around it. Each face edge will be tested twice. For
  // rigid-body moving mesh problems, we will later add a different rule
  // that uses pos_ic so that we do not do this calculation every time
  // we rigidly move the mesh, which does not change the maximum edge
  // length.
  class MaxEdgeLengthApplyInterior : public apply_rule<store<real>,
  Loci::Maximum<real> > {
    private:
      const_Map cl,cr ;
      const_multiMap face2edge ;
      const_MapVec<2> edge2node ;
      const_param<real> twoDimensionFactor ;
      const_store<vect3d> pos ;
      store<real> maxEdgeLength ;
    public:

      // Define input and output.
      MaxEdgeLengthApplyInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("face2edge",face2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("twoDimensionFactor",twoDimensionFactor) ;
        name_store("pos",pos) ;
        name_store("maxEdgeLength",maxEdgeLength) ;
        input("twoDimensionFactor,face2edge->edge2node->pos") ;
        output("(cl,cr)->maxEdgeLength") ;
        constraint("internalFaces") ;
      }

      // Test each edge for a face.
      void calculate(Entity face) {
        unsigned int numEdge=face2edge.num_elems(face) ;
        for(unsigned int i=0;i<numEdge;++i){
          unsigned int edgeNum=face2edge[face][i] ;
          unsigned int n0=edge2node[edgeNum][0],n1=edge2node[edgeNum][1] ;
          vect3d pos0=pos[n0],pos1=pos[n1] ;
          real edgeLength=sqrt((pos1.x-pos0.x)*(pos1.x-pos0.x)+(pos1.y-pos0.y)*
            (pos1.y-pos0.y)+(*twoDimensionFactor)*(pos1.z-pos0.z)*
            (pos1.z-pos0.z)) ;
          join(maxEdgeLength[cl[face]],edgeLength) ;
          join(maxEdgeLength[cr[face]],edgeLength) ;
        }
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<MaxEdgeLengthApplyInterior> registerMaxEdgeLengthApplyInterior ;

  // Add contributions from boundary faces.
  class MaxEdgeLengthApplyBoundary : public apply_rule<store<real>,
  Loci::Maximum<real> > {
    private:
      const_Map ci ;
      const_multiMap face2edge ;
      const_MapVec<2> edge2node ;
      const_param<real> twoDimensionFactor ;
      const_store<vect3d> pos ;
      store<real> maxEdgeLength ;
    public:

      // Define input and output.
      MaxEdgeLengthApplyBoundary() {
        name_store("ci",ci) ;
        name_store("face2edge",face2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("twoDimensionFactor",twoDimensionFactor) ;
        name_store("pos",pos) ;
        name_store("maxEdgeLength",maxEdgeLength) ;
        input("twoDimensionFactor,face2edge->edge2node->pos") ;
        output("ci->maxEdgeLength") ;
        constraint("boundaryFaces") ;
      }

      // Test each edge for a face.
      void calculate(Entity face) {
        unsigned int numEdge=face2edge.num_elems(face) ;
        for(unsigned int i=0;i<numEdge;++i){
          unsigned int edgeNum=face2edge[face][i] ;
          unsigned int n0=edge2node[edgeNum][0],n1=edge2node[edgeNum][1] ;
          vect3d pos0=pos[n0],pos1=pos[n1] ;
          real edgeLength=sqrt((pos1.x-pos0.x)*(pos1.x-pos0.x)+(pos1.y-pos0.y)*
            (pos1.y-pos0.y)+(*twoDimensionFactor)*(pos1.z-pos0.z)*
            (pos1.z-pos0.z)) ;
          join(maxEdgeLength[ci[face]],edgeLength) ;
        }
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<MaxEdgeLengthApplyBoundary> registerMaxEdgeLengthApplyBoundary ;

  // Priority rule for k-destruction multiplier when using Strelets' 2001 DES
  // turbulence model.
  class FDESStrelets: public pointwise_rule {
    private:
      const_param<real> betaStar ;
      const_store<real> f1,maxEdgeLength ;
      const_store<real> k,omega ;
      store<real> fDES ;
    public:

      // Define input and output.
      FDESStrelets() {
        name_store("betaStar",betaStar) ;
        name_store("f1",f1) ;
        name_store("maxEdgeLength",maxEdgeLength) ;
        name_store("k",k) ;
        name_store("omega",omega) ;
        name_store("DES2001::fDES",fDES) ;
        input("betaStar,f1,maxEdgeLength,k,omega") ;
        output("DES2001::fDES") ;
        constraint("DES2001TurbulenceModel,geom_cells") ;
      }

      // Compute for each cell.
      void calculate(Entity cell) {
        const real l=sqrt(k[cell])/((*betaStar)*omega[cell]) ;
        const real cDES=(1.0-f1[cell])*0.61+f1[cell]*0.78 ;
        fDES[cell]=max(l/(cDES*maxEdgeLength[cell]),1.0) ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FDESStrelets> registerFDESStrelets ;

  // Priority rule for k-destruction multiplier when using Menter's DDES
  // turbulence model. This model uses F2 as the shielding parameter which
  // Menter proposes as the default.
  class FDESMenter : public pointwise_rule {
    private:
      const_param<real> betaStar ;
      const_store<real> f1,f2,maxEdgeLength ;
      const_store<real> k,omega ;
      store<real> fDES ;
    public:

      // Define input and output.
      FDESMenter() {
        name_store("betaStar",betaStar) ;
        name_store("f1",f1) ;
        name_store("f2",f2) ;
        name_store("maxEdgeLength",maxEdgeLength) ;
        name_store("k",k) ;
        name_store("omega",omega) ;
        name_store("DDES::fDES",fDES) ;
        input("betaStar,f1,f2,maxEdgeLength,k,omega") ;
        output("DDES::fDES") ;
        constraint("DDESTurbulenceModel,geom_cells") ;
      }

      // Compute for each cell.
      void calculate(Entity cell) {
        const real l=sqrt(k[cell])/((*betaStar)*omega[cell]) ;
        const real cDES=(1.0-f1[cell])*0.61+f1[cell]*0.78 ;
        fDES[cell]=max((1.0-f2[cell])*l/(cDES*maxEdgeLength[cell]),1.0) ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FDESMenter> registerFDESMenter ;

  // Rule to compute the ratio of the turbulent length scale to the DES
  // length scale when using Strelets' 2001 DES turbulence model. This
  // variable is only used for grid diagnostic purposes.
  class LDESStrelets: public pointwise_rule {
    private:
      const_param<real> betaStar ;
      const_store<real> f1,maxEdgeLength ;
      const_store<real> k,omega ;
      store<real> lDES ;
    public:

      // Define input and output.
      LDESStrelets() {
        name_store("betaStar",betaStar) ;
        name_store("f1",f1) ;
        name_store("maxEdgeLength",maxEdgeLength) ;
        name_store("k",k) ;
        name_store("omega",omega) ;
        name_store("lDES",lDES) ;
        input("betaStar,f1,maxEdgeLength,k,omega") ;
        output("lDES") ;
        constraint("DES2001TurbulenceModel,geom_cells") ;
      }

      // Compute for each cell.
      void calculate(Entity cell) {
        const real l=sqrt(k[cell])/((*betaStar)*omega[cell]) ;
        const real cDES=(1.0-f1[cell])*0.61+f1[cell]*0.78 ;
        lDES[cell]=l/(cDES*maxEdgeLength[cell]) ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<LDESStrelets> registerLDESStrelets ;

  // Rule to compute the ratio of the turbulent length scale to the DES
  // length scale when using Menter's DDES turbulence model. This
  // variable is only used for grid diagnostic purposes.
  class LDESMenter : public pointwise_rule {
    private:
      const_param<real> betaStar ;
      const_store<real> f1,f2,maxEdgeLength ;
      const_store<real> k,omega ;
      store<real> lDES ;
    public:

      // Define input and output.
      LDESMenter() {
        name_store("betaStar",betaStar) ;
        name_store("f1",f1) ;
        name_store("f2",f2) ;
        name_store("maxEdgeLength",maxEdgeLength) ;
        name_store("k",k) ;
        name_store("omega",omega) ;
        name_store("lDES",lDES) ;
        input("betaStar,f1,f2,maxEdgeLength,k,omega") ;
        output("lDES") ;
        constraint("DDESTurbulenceModel,geom_cells") ;
      }

      // Compute for each cell.
      void calculate(Entity cell) {
        const real l=sqrt(k[cell])/((*betaStar)*omega[cell]) ;
        const real cDES=(1.0-f1[cell])*0.61+f1[cell]*0.78 ;
        lDES[cell]=(1.0-f2[cell])*l/(cDES*maxEdgeLength[cell]) ;
      }

      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<LDESMenter> registerLDESMenter ;

  // This rule provides lDES for the boundaries.
  class LDESBoundary : public pointwise_rule {
    private:
      const_Map ci ;
      const_store<real> lDES ;
      store<real> lDES_f ;
    public:

      // Define input and output.
      LDESBoundary() {
        name_store("ci",ci) ;
        name_store("lDES",lDES) ;
        name_store("lDES_f",lDES_f) ;
        input("ci->lDES") ;
        output("lDES_f") ;
        constraint("boundaryFaces") ;
      }

      // Compute for each face.
      void calculate(Entity face) { lDES_f[face]=lDES[ci[face]] ; }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<LDESBoundary> registerLDESBoundary ;
}
