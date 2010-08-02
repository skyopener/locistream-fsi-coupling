//-----------------------------------------------------------------------------
// Description: This file contains rules for the metrics associated with
// nodal control volumes.
//-----------------------------------------------------------------------------

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// StreamUns includes.
#include "const.h"
#include "sciTypes.h"

namespace streamUns {

  // This rule is the only rule required to provide facecenter{n-1}. We do not
  // have to explicitly define the entire iteration history. The quantities
  // facecenter{n} and facecenter{n,it} will be computed locally at {n} and
  // {n,it}. Since we always start with BDF, the value assigned here is never
  // used for the first timestep, so we can set it to zero.
  class TimeBuildFaceCenterBDF2: public pointwise_rule {
    private:
      store<vect3d> faceCenter ;
    public:

      // Define input and output.
      TimeBuildFaceCenterBDF2() {
        name_store("facecenter{n=-1}",faceCenter) ;
        output("facecenter{n=-1}") ;
        constraint("faces") ;
      }

      // Assign center for a single face.
      void calculate(Entity face) { faceCenter[face]=vect3d(0.0,0.0,0.0) ; }

      // Loop over faces.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildFaceCenterBDF2> registerTimeBuildFaceCenterBDF2 ;

  // This rule is the only rule required to provide vol{n-1}. We do not have
  // to explicitly define the entire iteration history for volume. The
  // quantities vol{n} and vol{n,it} will be computed locally at {n} and
  // {n,it}. Since we always start with BDF, the value assigned here is never
  // used for the first timestep, so we can set it to zero.
  class TimeBuildVolumeBDF2: public pointwise_rule {
    private:
      store<real> vol ;
    public:

      // Define input and output.
      TimeBuildVolumeBDF2() {
        name_store("vol{n=-1}",vol) ;
        output("vol{n=-1}") ;
        constraint("geom_cells") ;
      }

      // Assign volume for a single cell.
      void calculate(Entity cell) { vol[cell]=0.0 ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<TimeBuildVolumeBDF2> registerTimeBuildVolumeBDF2 ;

  // Unit rule to initialize cell volume.
  class GridMotionVolumeUnit : public unit_rule {
    private:
      store<real> vol ;
    public:
                                                                                
      // Define input and output.
      GridMotionVolumeUnit() {
        name_store("volGridMotion",vol) ;
        output("volGridMotion") ;
        constraint("geom_cells,gridMoverRbf") ;
      }
                                                                                
      // Initialize for a cell.
      void calculate(Entity cell) { vol[cell]=0.0 ; }
                                                                                
      // Loop over cells.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<GridMotionVolumeUnit> registerGridMotionVolumeUnit ;

  // Apply rule to add contributions from interior faces.
  class GridMotionVolumeApplyInterior : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_multiMap face2edge ;
      const_MapVec<2> edge2node ;
      const_store<vect3d> pos ;
      const_store<Area> area ;
      const_store<vect3d> faceCenter ;
      store<real> vol ;
    public:

      // Define input and output.
      GridMotionVolumeApplyInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("face2edge",face2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("pos",pos) ;
        name_store("area",area) ;
        name_store("facecenter",faceCenter) ;
        name_store("volGridMotion",vol) ;
        input("face2edge->edge2node->pos,area,facecenter") ;
        output("(cl,cr)->volGridMotion") ;
        constraint("internalFaces,gridMoverRbf") ;
      }

      // Increment the volume for the cells attached to the face.
      void calculate(Entity face) {
        unsigned int numEdge=face2edge.num_elems(face) ;
        for(unsigned int i=0;i<numEdge;++i){
          unsigned int edgeNum=face2edge[face][i] ;
          unsigned int n0=edge2node[edgeNum][0],n1=edge2node[edgeNum][1] ;
          vect3d pos0=pos[n0],pos1=pos[n1] ;
          vect3d dA=cross(faceCenter[face]-pos0,faceCenter[face]-pos1) ;
          real dV=(1.0/6.0)*(faceCenter[face].x+pos0.x+pos1.x)*dA.x ;
          if(dot(dA,area[face].n)>0.0){
            vol[cl[face]]+=dV ; vol[cr[face]]-=dV ;
          }else{
            vol[cl[face]]-=dV ; vol[cr[face]]+=dV ;
          }
        }
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GridMotionVolumeApplyInterior>
    registerGridMotionVolumeApplyInterior ;

  class GridMotionRenameVol : public pointwise_rule {
    private:
	  const_store<real> volGridMotion;
	  store<real> vol ;
    public:
	  GridMotionRenameVol() {	 
	    name_store("volGridMotion", volGridMotion) ;
		name_store("gridMotion::vol", vol) ;
		input("volGridMotion") ;
		output("gridMotion::vol=volGridMotion") ;
	}
  
    // Empty
      virtual void compute(const sequence &seq) {}
  
  };

  register_rule<GridMotionRenameVol>
    registerGridMotionRenameVol ;

  // Apply rule to add contributions from boundary faces.
  class GridMotionVolumeApplyBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_multiMap face2edge ;
      const_MapVec<2> edge2node ;
      const_store<vect3d> pos ;
      const_store<Area> area ;
      const_store<vect3d> faceCenter ;
      store<real> vol ;
    public:

      // Define input and output.
      GridMotionVolumeApplyBoundary() {
        name_store("ci",ci) ;
        name_store("face2edge",face2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("pos",pos) ;
        name_store("area",area) ;
        name_store("facecenter",faceCenter) ;
        name_store("volGridMotion",vol) ;
        input("face2edge->edge2node->pos,area,facecenter") ;
        output("ci->volGridMotion") ;
        constraint("boundaryFaces,gridMoverRbf") ;
      }

      // Increment the volume for the cells attached to the face.
      void calculate(Entity face) {
        unsigned int numEdge=face2edge.num_elems(face) ;
        for(unsigned int i=0;i<numEdge;++i){
          unsigned int edgeNum=face2edge[face][i] ;
          unsigned int n0=edge2node[edgeNum][0],n1=edge2node[edgeNum][1] ;
          vect3d pos0=pos[n0],pos1=pos[n1] ;
          vect3d dA=cross(faceCenter[face]-pos0,faceCenter[face]-pos1) ;
          real dV=(1.0/6.0)*(faceCenter[face].x+pos0.x+pos1.x)*dA.x ;
          vol[ci[face]]+=(dot(dA,area[face].n)>0.0)? dV:-dV ;
        }
      }

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<GridMotionVolumeApplyBoundary>
    registerGridMotionVolumeApplyBoundary ;

//-----------------------------------------------------------------------------
// Rules to compute the area of the faces associated with the dual mesh. The
// face normal is defined to point away from edge node 0.

  // Unit rule.
  class DualAreaUnit : public unit_rule {
    private:
      store<vect3d> dualArea ;
    public:
                                                                                
      // Define input and output.
      DualAreaUnit() {
        name_store("dualArea",dualArea) ;
        output("dualArea") ;
        constraint("edges") ;
      }
                                                                                
      // Initialize for an edge.
      void calculate(Entity edge) { dualArea[edge]=vect3d(0.0,0.0,0.0) ; }
                                                                                
      // Loop over edges.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<DualAreaUnit> registerDualAreaUnit ;

  // Add contributions from interior faces.
  class DualAreaApplyInterior : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_Map cl,cr ;
      const_multiMap face2edge ;
      const_MapVec<2> edge2node ;
      const_store<vect3d> pos ;
      const_store<vect3d> faceCenter ;
      const_store<vect3d> cellCenter ;
      store<vect3d> dualArea ;
    public:
                                                                                
      // Define input and output.
      DualAreaApplyInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("face2edge",face2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("pos",pos) ;
        name_store("facecenter",faceCenter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("dualArea",dualArea) ;
        input("facecenter,face2edge->edge2node->pos,(cl,cr)->cellcenter") ;
        output("face2edge->dualArea") ;
        constraint("internalFaces") ;
      }
                                                                                
      // Increment the dual area for the edges attached to the face.
      void calculate(Entity face) {
        unsigned int numEdge=face2edge.num_elems(face) ;
        for(unsigned int i=0;i<numEdge;++i){
          unsigned int edgeNum=face2edge[face][i] ;
          vect3d pos0=pos[edge2node[edgeNum][0]] ;
          vect3d pos1=pos[edge2node[edgeNum][1]] ;
          vect3d edgeCenter=0.5*(pos0+pos1) ;
          vect3d dA0=0.5*cross(cellCenter[cl[face]]-edgeCenter,faceCenter[face]-
            edgeCenter) ;
          vect3d dA1=0.5*cross(cellCenter[cr[face]]-edgeCenter,faceCenter[face]-
            edgeCenter) ;
          dualArea[edgeNum]+=(dot(dA0,pos1-pos0)>0.0)? dA0:-1.0*dA0 ;
          dualArea[edgeNum]+=(dot(dA1,pos1-pos0)>0.0)? dA1:-1.0*dA1 ;
        }
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<DualAreaApplyInterior> registerDualAreaApplyInterior ;

  // Add contributions from boundary faces.
  class DualAreaApplyBoundary : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_Map ci ;
      const_multiMap face2edge ;
      const_MapVec<2> edge2node ;
      const_store<vect3d> pos ;
      const_store<vect3d> faceCenter ;
      const_store<vect3d> cellCenter ;
      store<vect3d> dualArea ;
    public:
                                                                                
      // Define input and output.
      DualAreaApplyBoundary() {
        name_store("ci",ci) ;
        name_store("face2edge",face2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("pos",pos) ;
        name_store("facecenter",faceCenter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("dualArea",dualArea) ;
        input("facecenter,face2edge->edge2node->pos,ci->cellcenter") ;
        output("face2edge->dualArea") ;
        constraint("boundaryFaces") ;
      }
                                                                                
      // Increment the dual area for the edges attached to the face.
      void calculate(Entity face) {
        unsigned int numEdge=face2edge.num_elems(face) ;
        for(unsigned int i=0;i<numEdge;++i){
          unsigned int edgeNum=face2edge[face][i] ;
          vect3d pos0=pos[edge2node[edgeNum][0]] ;
          vect3d pos1=pos[edge2node[edgeNum][1]] ;
          vect3d edgeCenter=0.5*(pos0+pos1) ;
          vect3d dA0=0.5*cross(cellCenter[ci[face]]-edgeCenter,faceCenter[face]-
            edgeCenter) ;
          dualArea[edgeNum]+=(dot(dA0,pos1-pos0)>0.0)? dA0:-1.0*dA0 ;
        }
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<DualAreaApplyBoundary> registerDualAreaApplyBoundary ;

//-----------------------------------------------------------------------------
// Rules to compute the dual cell volumes which are centered about the nodes.

  // Unit rule.
  class DualVolumeUnit : public unit_rule {
    private:
      store<real> dualVolume ;
    public:
                                                                                
      // Define input and output.
      DualVolumeUnit() {
        name_store("dualVolume",dualVolume) ;
        output("dualVolume") ;
        constraint("UNIVERSE") ;
      }
                                                                                
      // Initialize for a node.
      void calculate(Entity node) { dualVolume[node]=0.0 ; }
                                                                                
      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<DualVolumeUnit> registerDualVolumeUnit ;

  // Add contributions from interior faces.
  class DualVolumeApplyInterior : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_multiMap face2edge ;
      const_MapVec<2> edge2node ;
      const_store<vect3d> pos ;
      const_store<vect3d> faceCenter ;
      const_store<vect3d> cellCenter ;
      store<real> dualVolume ;
    public:
                                                                                
      // Define input and output.
      DualVolumeApplyInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("face2edge",face2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("pos",pos) ;
        name_store("facecenter",faceCenter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("dualVolume",dualVolume) ;
        input("facecenter,face2edge->edge2node->pos,(cl,cr)->cellcenter") ;
        output("face2edge->edge2node->dualVolume") ;
        constraint("internalFaces") ;
      }
                                                                                
      // Increment the dual volume for the nodes attached to the face.
      void calculate(Entity face) {
        unsigned int numEdge=face2edge.num_elems(face) ;
        for(unsigned int i=0;i<numEdge;++i){
          unsigned int edgeNum=face2edge[face][i] ;
          unsigned int n0=edge2node[edgeNum][0],n1=edge2node[edgeNum][1] ;
          vect3d pos0=pos[n0],pos1=pos[n1],edgeCenter=0.5*(pos0+pos1) ;
          vect3d dA0=0.5*cross(cellCenter[cl[face]]-edgeCenter,faceCenter[face]-
            edgeCenter) ;
          real dV0=(1.0/3.0)*(edgeCenter.x+cellCenter[cl[face]].x+
            faceCenter[face].x)*dA0.x ;
          vect3d dA1=0.5*cross(cellCenter[cr[face]]-edgeCenter,faceCenter[face]-
            edgeCenter) ;
          real dV1=(1.0/3.0)*(edgeCenter.x+cellCenter[cr[face]].x+
            faceCenter[face].x)*dA1.x ;
          if(dot(dA0,pos1-pos0)>0.0){
            dualVolume[n0]+=dV0 ; dualVolume[n1]-=dV0 ;
          }else{
            dualVolume[n0]-=dV0 ; dualVolume[n1]+=dV0 ;
          }
          if(dot(dA1,pos1-pos0)>0.0){
            dualVolume[n0]+=dV1 ; dualVolume[n1]-=dV1 ;
          }else{
            dualVolume[n0]-=dV1 ; dualVolume[n1]+=dV1 ;
          }
        }
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<DualVolumeApplyInterior> registerDualVolumeApplyInterior ;

  // Add contributions from boundary faces.
  class DualVolumeApplyBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_multiMap face2edge ;
      const_MapVec<2> edge2node ;
      const_store<vect3d> pos ;
      const_store<vect3d> faceCenter ;
      const_store<vect3d> cellCenter ;
      store<real> dualVolume ;
    public:
                                                                                
      // Define input and output.
      DualVolumeApplyBoundary() {
        name_store("ci",ci) ;
        name_store("face2edge",face2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("pos",pos) ;
        name_store("facecenter",faceCenter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("dualVolume",dualVolume) ;
        input("facecenter,face2edge->edge2node->pos,ci->cellcenter") ;
        output("face2edge->edge2node->dualVolume") ;
        constraint("boundaryFaces") ;
      }
                                                                                
      // Increment the dual volume for the nodes attached to the face.
      void calculate(Entity face) {
        unsigned int numEdge=face2edge.num_elems(face) ;
        for(unsigned int i=0;i<numEdge;++i){
          unsigned int edgeNum=face2edge[face][i] ;
          unsigned int n0=edge2node[edgeNum][0],n1=edge2node[edgeNum][1] ;
          vect3d pos0=pos[n0],pos1=pos[n1],edgeCenter=0.5*(pos0+pos1) ;

          // Contribution from the facet interior to the domain.
          {
            vect3d dA0=0.5*cross(cellCenter[ci[face]]-edgeCenter,
              faceCenter[face]-edgeCenter) ;
            real dV0=(1.0/3.0)*(edgeCenter.x+cellCenter[ci[face]].x+
              faceCenter[face].x)*dA0.x ;
            if(dot(dA0,pos1-pos0)>0.0){
              dualVolume[n0]+=dV0 ; dualVolume[n1]-=dV0 ;
            }else{
              dualVolume[n0]-=dV0 ; dualVolume[n1]+=dV0 ;
            }
          }

          // Contribution from the facet on the boundary of the domain.
          {
            vect3d dA0=0.5*cross(edgeCenter-pos0,faceCenter[face]-pos0) ;
            real dV0=(1.0/3.0)*(edgeCenter.x+pos0.x+faceCenter[face].x)*dA0.x ;
            dualVolume[n0]+=(dot(dA0,cellCenter[ci[face]]-edgeCenter)<0.0)?
              dV0:-dV0 ;
            vect3d dA1=0.5*cross(edgeCenter-pos1,faceCenter[face]-pos1) ;
            real dV1=(1.0/3.0)*(edgeCenter.x+pos1.x+faceCenter[face].x)*dA1.x ;
            dualVolume[n1]+=(dot(dA1,cellCenter[ci[face]]-edgeCenter)<0.0)?
              dV1:-dV1 ;
          }
        }
      }
                                                                                
      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<DualVolumeApplyBoundary> registerDualVolumeApplyBoundary ;

//-----------------------------------------------------------------------------
// Rules to compute the diffusion volumes which are centered about the edges.
// This volume is used in the discretization of the diffusion operator.

  // Unit rule.
  class DiffusionVolumeUnit : public unit_rule {
    private:
      store<real> diffusionVolume ;
    public:
                                                                                
      // Define input and output.
      DiffusionVolumeUnit() {
        name_store("diffusionVolume",diffusionVolume) ;
        output("diffusionVolume") ;
        constraint("UNIVERSE") ;
      }
                                                                                
      // Initialize for a node.
      void calculate(Entity node) { diffusionVolume[node]=0.0 ; }
                                                                                
      // Loop over nodes.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;
                                                                                
  register_rule<DiffusionVolumeUnit> registerDiffusionVolumeUnit ;

  // Add contributions from interior faces.
  class DiffusionVolumeApplyInterior : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map cl,cr ;
      const_multiMap face2edge ;
      const_MapVec<2> edge2node ;
      const_store<vect3d> pos ;
      const_store<vect3d> faceCenter ;
      const_store<vect3d> cellCenter ;
      store<real> diffusionVolume ;
    public:

      // Define input and output.
      DiffusionVolumeApplyInterior() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("face2edge",face2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("pos",pos) ;
        name_store("facecenter",faceCenter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("diffusionVolume",diffusionVolume) ;
        input("facecenter,face2edge->edge2node->pos,(cl,cr)->cellcenter") ;
        output("face2edge->diffusionVolume") ;
        constraint("internalFaces") ;
      }

      // Increment the diffusion volume for the edges attached to the face.
      void calculate(Entity face) {
        unsigned int numEdge=face2edge.num_elems(face) ;
        for(unsigned int i=0;i<numEdge;++i){
          unsigned int edgeNum=face2edge[face][i] ;
          unsigned int n0=edge2node[edgeNum][0],n1=edge2node[edgeNum][1] ;
          vect3d pos0=pos[n0],pos1=pos[n1],edgeCenter=0.5*(pos0+pos1) ;

          // Left cell.
          {
            vect3d dA0=cross(cellCenter[cl[face]]-pos0,faceCenter[face]-pos0) ;
            real dV0=(1.0/6.0)*(pos0.x+cellCenter[cl[face]].x+faceCenter
              [face].x)*dA0.x ;
            diffusionVolume[edgeNum]+=(dot(dA0,faceCenter[face]-
              edgeCenter)>0.0)? dV0:-dV0 ;
            vect3d dA1=cross(cellCenter[cl[face]]-pos1,faceCenter[face]-pos1) ;
            real dV1=(1.0/6.0)*(pos1.x+cellCenter[cl[face]].x+faceCenter
              [face].x)*dA1.x ;
            diffusionVolume[edgeNum]+=(dot(dA1,faceCenter[face]-
              edgeCenter)>0.0)? dV1:-dV1 ;
          }

          // Right cell.
          {
            vect3d dA0=cross(cellCenter[cr[face]]-pos0,faceCenter[face]-pos0) ;
            real dV0=(1.0/6.0)*(pos0.x+cellCenter[cr[face]].x+faceCenter
              [face].x)*dA0.x ;
            diffusionVolume[edgeNum]+=(dot(dA0,faceCenter[face]-
              edgeCenter)>0.0)? dV0:-dV0 ;
            vect3d dA1=cross(cellCenter[cr[face]]-pos1,faceCenter[face]-pos1) ;
            real dV1=(1.0/6.0)*(pos1.x+cellCenter[cr[face]].x+faceCenter
              [face].x)*dA1.x ;
            diffusionVolume[edgeNum]+=(dot(dA1,faceCenter[face]-
              edgeCenter)>0.0)? dV1:-dV1 ;
          }
        }
      } 

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusionVolumeApplyInterior>
    registerDiffusionVolumeApplyInterior ;

  // Add contributions from boundary faces.
  class DiffusionVolumeApplyBoundary : public apply_rule<store<real>,
  Loci::Summation<real> > {
    private:
      const_Map ci ;
      const_multiMap face2edge ;
      const_MapVec<2> edge2node ;
      const_store<vect3d> pos ;
      const_store<vect3d> faceCenter ;
      const_store<vect3d> cellCenter ;
      store<real> diffusionVolume ;
    public:

      // Define input and output.
      DiffusionVolumeApplyBoundary() {
        name_store("ci",ci) ;
        name_store("face2edge",face2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("pos",pos) ;
        name_store("facecenter",faceCenter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("diffusionVolume",diffusionVolume) ;
        input("facecenter,face2edge->edge2node->pos,ci->cellcenter") ;
        output("face2edge->diffusionVolume") ;
        constraint("boundaryFaces") ;
      }

      // Increment the diffusion volume for the edges attached to the face.
      void calculate(Entity face) {
        unsigned int numEdge=face2edge.num_elems(face) ;
        for(unsigned int i=0;i<numEdge;++i){
          unsigned int edgeNum=face2edge[face][i] ;
          unsigned int n0=edge2node[edgeNum][0],n1=edge2node[edgeNum][1] ;
          vect3d pos0=pos[n0],pos1=pos[n1],edgeCenter=0.5*(pos0+pos1) ;

          // Interior contribution.
          {
            vect3d dA0=cross(cellCenter[ci[face]]-pos0,faceCenter[face]-pos0) ;
            real dV0=(1.0/6.0)*(pos0.x+cellCenter[ci[face]].x+faceCenter
              [face].x)*dA0.x ;
            diffusionVolume[edgeNum]+=(dot(dA0,faceCenter[face]-
              edgeCenter)>0.0)? dV0:-dV0 ;
            vect3d dA1=cross(cellCenter[ci[face]]-pos1,faceCenter[face]-pos1) ;
            real dV1=(1.0/6.0)*(pos1.x+cellCenter[ci[face]].x+faceCenter
              [face].x)*dA1.x ;
            diffusionVolume[edgeNum]+=(dot(dA1,faceCenter[face]-
              edgeCenter)>0.0)? dV1:-dV1 ;
          }

          // Boundary contribution.
          {
            vect3d dA0=cross(faceCenter[face]-pos0,edgeCenter-pos0) ;
            real dV0=(1.0/6.0)*(pos0.x+edgeCenter.x+faceCenter[face].x)*dA0.x ;
            diffusionVolume[edgeNum]+=(dot(dA0,faceCenter[face]-
              cellCenter[ci[face]])>0.0)? dV0:-dV0 ;
            vect3d dA1=cross(faceCenter[face]-pos1,edgeCenter-pos1) ;
            real dV1=(1.0/6.0)*(pos1.x+edgeCenter.x+faceCenter[face].x)*dA1.x ;
            diffusionVolume[edgeNum]+=(dot(dA1,faceCenter[face]-
              cellCenter[ci[face]])>0.0)? dV1:-dV1 ;
          }
        }
      } 

      // Loop over faces.
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<DiffusionVolumeApplyBoundary>
    registerDiffusionVolumeApplyBoundary ;
}
