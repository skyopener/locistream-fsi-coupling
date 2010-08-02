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
                                                                                
// StreamUns includes.
#include "const.h"
#include "move.h"
#include "residual.h"
#include "sciTypes.h"
                                                                                
namespace streamUns {

  // Rule for creating constraint for all the nodes.
      NodeConstraint() {
        name_store("face2node",face2node) ;
        name_store("nodes",nodes) ;
        input("face2node") ;
        output("face2node->nodes") ;
        constraint("faces") ;
      }
      BoundaryNodeConstraint() {
        name_store("face2node",face2node) ;
        name_store("boundaryNodes",boundaryNodes) ;
        input("face2node") ;
        output("face2node->boundaryNodes") ;
        constraint("boundaryFaces") ;
      }
      EdgeConstraint() {
        name_store("face2edge",face2edge) ;
        name_store("edges",edges) ;
        input("face2edge") ;
        output("face2edge->edges") ;
        constraint("faces") ;
      }
      BoundaryEdgeConstraint() {
        name_store("face2edge",face2edge) ;
        name_store("boundaryEdges",boundaryEdges) ;
        input("face2edge") ;
        output("face2edge->boundaryEdges") ;
        constraint("boundaryFaces") ;
      }
      NodePosition() {
        name_store("pos",pos) ;
        name_store("nodePosition",nodePosition) ;
        output("nodePosition") ;
        constraint("nodes") ;
      }
      BoundaryNodePosition() {
        name_store("pos",pos) ;
        name_store("boundaryNodePosition",boundaryNodePosition) ;
        output("boundaryNodePosition") ;
        constraint("boundaryNodes") ;
      }
      MinimumDualVolumeUnit() {
        name_store("minDualVolume",minDualVolume) ;
        output("minDualVolume") ;
        constraint("UNIVERSE") ;
      }
      MinimumDualVolumeApply() {
        name_store("dualVolume",dualVolume) ;
        name_store("minDualVolume",minDualVolume) ;
        input("dualVolume") ;
        output("minDualVolume") ;
        constraint("nodes") ;
      }
      MaximumDualVolumeUnit() {
        name_store("maxDualVolume",maxDualVolume) ;
        output("maxDualVolume") ;
        constraint("UNIVERSE") ;
      }
      MaximumDualVolumeApply() {
        name_store("dualVolume",dualVolume) ;
        name_store("maxDualVolume",maxDualVolume) ;
        input("dualVolume") ;
        output("maxDualVolume") ;
        constraint("nodes") ;
      }
      DefaultChi() {
        name_store("gridMoverChi",gridMoverChi) ;
        name_store("chi",chi) ;
        input("gridMoverChi") ;
        output("chi") ;
        constraint("edges") ;
      }
      NodeChiFromDualVolume() {
        name_store("gridMoverMinChi",gridMoverMinChi) ;
        name_store("gridMoverMaxChi",gridMoverMaxChi) ;
        name_store("minDualVolume",minDualVolume) ;
        name_store("maxDualVolume",maxDualVolume) ;
        name_store("dualVolume",dualVolume) ;
        name_store("node_chi",node_chi) ;
        input("gridMoverMinChi,gridMoverMaxChi,minDualVolume,maxDualVolume") ;
        input("dualVolume") ;
        output("node_chi") ;
        constraint("nodes") ;
      }
      EdgeChi() {
        name_store("edge2node",edge2node) ;
        name_store("node_chi",node_chi) ;
        name_store("priority::chi",chi) ;
        input("edge2node->node_chi") ;
        output("priority::chi") ;
        constraint("edges") ;
      }
      GridMoverFactor0() {
        name_store("chi",chi) ;
        name_store("diffusionVolume",diffusionVolume) ;
        name_store("gridMoverFactor0",gridMoverFactor0) ;
        input("chi,diffusionVolume") ;
        output("gridMoverFactor0") ;
        constraint("edges") ;
      }
      ScaledDualArea() {
        name_store("gridMoverNu",nu) ;
        name_store("dualArea",dualArea) ;
        name_store("gridMoverFactor0",gridMoverFactor0) ;
        name_store("scaledDualArea",scaledDualArea) ;
        input("gridMoverNu,dualArea,gridMoverFactor0") ;
        output("scaledDualArea") ;
        constraint("edges") ;
      }
      InternalFacetAreaSum() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("face2edge",face2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("pos",pos) ;
        name_store("facecenter",faceCenter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("diffusionVolume",diffusionVolume) ;
        name_store("dASumLeft",dASumLeft) ;
        name_store("dASumRight",dASumRight) ;
        input("face2edge->edge2node->pos,(cl,cr)->cellcenter,facecenter") ;
        input("face2edge->diffusionVolume") ;
        output("dASumLeft,dASumRight") ;
        constraint("internalFaces") ;
      }
      BoundaryFacetAreaSum() {
        name_store("ci",ci) ;
        name_store("face2edge",face2edge) ;
        name_store("edge2node",edge2node) ;
        name_store("pos",pos) ;
        name_store("facecenter",faceCenter) ;
        name_store("cellcenter",cellCenter) ;
        name_store("diffusionVolume",diffusionVolume) ;
        name_store("dASum",dASum) ;
        name_store("dABoundary",dABoundary) ;
        input("face2edge->edge2node->pos,ci->cellcenter,facecenter") ;
        input("face2edge->diffusionVolume") ;
        output("dASum,dABoundary") ;
        constraint("boundaryFaces") ;
      }
      InitializeDisplacementMainCoefficient() {
        name_store("sMainCoefficient{n}",sMainCoefficient) ;
        output("sMainCoefficient{n}") ;
        constraint("nodes{n}") ;
      }
      DiffusionToDisplacementMainCoefficient() {
        name_store("edge2node{n}",edge2node) ;
        name_store("gridMoverFactor0{n}",gridMoverFactor0) ;
        name_store("dualArea{n}",dualArea) ;
        name_store("diffusionVolume{n}",diffusionVolume) ;
        name_store("sMainCoefficient{n}",sMainCoefficient) ;
        input("gridMoverFactor0{n},dualArea{n},diffusionVolume{n}") ;
        output("edge2node{n}->sMainCoefficient{n}") ;
        constraint("edges{n}") ;
      }
      InitializeDisplacementSourceTerm() {
        name_store("sSourceTerm{n,itg}",sSourceTerm) ;
        output("sSourceTerm{n,itg}") ;
        constraint("nodes{n,itg}") ;
      }
      DilatationToDisplacementSourceTerm() {
        name_store("gridMoverNu{n,itg}",nu) ;
        name_store("edge2node{n,itg}",edge2node) ;
        name_store("gridMoverFactor0{n}",gridMoverFactor0) ;
        name_store("dualArea{n}",dualArea) ;
        name_store("scaledDualArea{n}",scaledDualArea) ;
        name_store("diffusionVolume{n}",diffusionVolume) ;
        name_store("node_sStar{n,itg}",node_s) ;
        name_store("sSourceTerm{n,itg}",sSourceTerm) ;
        input("gridMoverNu{n,itg},gridMoverFactor0{n},dualArea{n}") ;
        input("diffusionVolume{n},edge2node{n,itg}->node_sStar{n,itg}") ;
        output("edge2node{n,itg}->sSourceTerm{n,itg}") ;
        constraint("edges{n,itg}") ;
      }
      CrossDiffusionDilatationApplyInterior() {
        name_store("cl{n,itg}",cl) ;
        name_store("cr{n,itg}",cr) ;
        name_store("face2edge{n,itg}",face2edge) ;
        name_store("edge2node{n,itg}",edge2node) ;
        name_store("gridMoverNu{n,itg}",nu) ;
        name_store("gridMoverFactor0{n}",gridMoverFactor0) ;
        name_store("dualArea{n}",dualArea) ;
        name_store("scaledDualArea{n}",scaledDualArea) ;
        name_store("diffusionVolume{n}",diffusionVolume) ;
        name_store("face(node_sStar){n,itg}",face_s) ;
        name_store("cell(node_sStar){n,itg}",cell_s) ;
        name_store("dASumLeft{n}",dASumLeft) ;
        name_store("dASumRight{n}",dASumRight) ;
        name_store("sSourceTerm{n,itg}",sSourceTerm) ;
        input("gridMoverNu{n,itg},face2edge{n,itg}->gridMoverFactor0{n}") ;
        input("face2edge{n,itg}->(dualArea{n},diffusionVolume{n})") ;
        input("face2edge{n,itg}->scaledDualArea{n}") ;
        input("face(node_sStar){n,itg}") ;
        input("(cl{n,itg},cr{n,itg})->cell(node_sStar){n,itg}") ;
        input("dASumLeft{n},dASumRight{n}") ;
        output("face2edge{n,itg}->edge2node{n,itg}->sSourceTerm{n,itg}") ;
        constraint("internalFaces{n,itg}") ;
      }
      CrossDiffusionDilatationApplyBoundary() {
        name_store("ci{n,itg}",ci) ;
        name_store("face2edge{n,itg}",face2edge) ;
        name_store("edge2node{n,itg}",edge2node) ;
        name_store("gridMoverNu{n,itg}",nu) ;
        name_store("gridMoverFactor0{n}",gridMoverFactor0) ;
        name_store("dualArea{n}",dualArea) ;
        name_store("scaledDualArea{n}",scaledDualArea) ;
        name_store("diffusionVolume{n}",diffusionVolume) ;
        name_store("node_sStar{n,itg}",node_s) ;
        name_store("face(node_sStar){n,itg}",face_s) ;
        name_store("cell(node_sStar){n,itg}",cell_s) ;
        name_store("dASum{n}",dASum) ;
        name_store("dABoundary{n}",dABoundary) ;
        name_store("sSourceTerm{n,itg}",sSourceTerm) ;
        input("gridMoverNu{n,itg}") ;
        input("face2edge{n,itg}->edge2node{n,itg}->node_sStar{n,itg}");
        input("face2edge{n,itg}->gridMoverFactor0{n}") ;
        input("face2edge{n,itg}->(dualArea{n},diffusionVolume{n})") ;
        input("face2edge{n,itg}->scaledDualArea{n}") ;
        input("face(node_sStar){n,itg},ci{n,itg}->cell(node_sStar){n,itg}") ;
        input("dASum{n},dABoundary{n}") ;
        output("face2edge{n,itg}->edge2node{n,itg}->sSourceTerm{n,itg}") ;
        constraint("boundaryFaces{n,itg}") ;
      }
      RigidBodyRotationDiffusionToDisplacementSourceTerm() {
        name_store("gridMoverKappa{n,itg}",gridMoverKappa) ;
        name_store("gridMoverFactor0{n}",gridMoverFactor0) ;
        name_store("edge2node{n,itg}",edge2node) ;
        name_store("dualArea{n}",dualArea) ;
        name_store("diffusionVolume{n}",diffusionVolume) ;
        name_store("nodeGrad(node_sStar){n,itg}",node_sGrad) ;
        name_store("sSourceTerm{n,itg}",sSourceTerm) ;
        input("gridMoverKappa{n,itg},gridMoverFactor0{n},dualArea{n}") ;
        input("diffusionVolume{n}") ;
        input("edge2node{n,itg}->nodeGrad(node_sStar){n,itg}") ;
        output("edge2node{n,itg}->sSourceTerm{n,itg}") ;
        constraint("edges{n,itg}") ;
      }
      DisplacementMatrixDiagonal() {
        name_store("gridMoverRelaxation{n}",gridMoverRelaxation) ;
        name_store("sMainCoefficient{n}",sMainCoefficient) ;
        name_store("sStar_D{n}",D) ;
        input("gridMoverRelaxation{n},sMainCoefficient{n}") ;
        output("sStar_D{n}") ;
        constraint("nodes{n}") ;
      }
      DisplacementMatrixDiagonalSpecified() {
        name_store("priority::sStar_D{n}",D) ;
        output("priority::sStar_D{n}") ;
        constraint("boundaryDisplacement{n}") ;
      }
      EdgeCoefficientMultiplierDefault() {
        name_store("edgeCoefficientMultiplier",edgeCoefficientMultiplier) ;
        output("edgeCoefficientMultiplier") ;
        constraint("nodes") ;
      }
      EdgeCoefficientMultiplierPriority() {
        name_store("boundaryDisplacement::edgeCoefficientMultiplier",
          edgeCoefficientMultiplier) ;
        output("boundaryDisplacement::edgeCoefficientMultiplier") ;
        constraint("boundaryDisplacement") ;
      }
      DisplacementMatrixEdgeCoefficient() {
        name_store("edge2node{n}",edge2node) ;
        name_store("gridMoverFactor0{n}",gridMoverFactor0) ;
        name_store("dualArea{n}",dualArea) ;
        name_store("diffusionVolume{n}",diffusionVolume) ;
        name_store("edgeCoefficientMultiplier{n}",edgeCoefficientMultiplier) ;
        name_store("sStar_E{n}",sStar_E) ;
        input("gridMoverFactor0{n},dualArea{n},diffusionVolume{n}") ;
        input("edge2node{n}->edgeCoefficientMultiplier{n}") ;
        output("sStar_E{n}") ;
        constraint("edges{n}") ;
      }
      DisplacementMatrixEdgeCoefficientNoSymmetry() {
        name_store("boundaryNoSymmetryEdges::sStar_E{n}",sStar_E) ;
        output("boundaryNoSymmetryEdges::sStar_E{n}") ;
        constraint("boundaryNoSymmetryEdges{n}") ;
      }
      DisplacementRHS() {
        name_store("twoDimensionFactor{n,itg}",twoDimensionFactor) ;
        name_store("gridMoverRelaxation{n,itg}",gridMoverRelaxation) ;
        name_store("node_sStar{n,itg}",node_sStar) ;
        name_store("sMainCoefficient{n,itg}",sMainCoefficient) ;
        name_store("sSourceTerm{n,itg}",sSourceTerm) ;
        name_store("sStar_B{n,itg}",B) ;
        input("gridMoverRelaxation{n,itg},node_sStar{n,itg}") ;
        input("sMainCoefficient{n,itg}") ;
        input("sSourceTerm{n,itg},twoDimensionFactor{n,itg}") ;
        output("sStar_B{n,itg}") ;
        constraint("nodes{n,itg}") ;
      }
      DisplacementRHSSXZero() {
        name_store("twoDimensionFactor{n,itg}",twoDimensionFactor) ;
        name_store("gridMoverRelaxation{n,itg}",gridMoverRelaxation) ;
        name_store("node_sStar{n,itg}",node_sStar) ;
        name_store("sMainCoefficient{n,itg}",sMainCoefficient) ;
        name_store("sSourceTerm{n,itg}",sSourceTerm) ;
        name_store("priority::sStar_B{n,itg}",B) ;
        input("gridMoverRelaxation{n,itg},node_sStar{n,itg}") ;
        input("sMainCoefficient{n,itg},sSourceTerm{n,itg}") ;
        input("twoDimensionFactor{n,itg}") ;
        output("priority::sStar_B{n,itg}") ;
        constraint("sXZeroNodes{n,itg}") ;
      }
      DisplacementRHSSYZero() {
        name_store("twoDimensionFactor{n,itg}",twoDimensionFactor) ;
        name_store("gridMoverRelaxation{n,itg}",gridMoverRelaxation) ;
        name_store("node_sStar{n,itg}",node_sStar) ;
        name_store("sMainCoefficient{n,itg}",sMainCoefficient) ;
        name_store("sSourceTerm{n,itg}",sSourceTerm) ;
        name_store("priority::sStar_B{n,itg}",B) ;
        input("gridMoverRelaxation{n,itg},node_sStar{n,itg}") ;
        input("sMainCoefficient{n,itg},sSourceTerm{n,itg}") ;
        input("twoDimensionFactor{n,itg}") ;
        output("priority::sStar_B{n,itg}") ;
        constraint("sYZeroNodes{n,itg}") ;
      }
      DisplacementRHSSZZero() {
        name_store("gridMoverRelaxation{n,itg}",gridMoverRelaxation) ;
        name_store("node_sStar{n,itg}",node_sStar) ;
        name_store("sMainCoefficient{n,itg}",sMainCoefficient) ;
        name_store("sSourceTerm{n,itg}",sSourceTerm) ;
        name_store("priority::sStar_B{n,itg}",B) ;
        input("gridMoverRelaxation{n,itg},node_sStar{n,itg}") ;
        input("sMainCoefficient{n,itg},sSourceTerm{n,itg}") ;
        output("priority::sStar_B{n,itg}") ;
        constraint("sZZeroNodes{n,itg}") ;
      }
      DisplacementRHSSpecified() {
        name_store("node_s_b{n,itg}",node_s_b) ;
        name_store("priority::priority::sStar_B{n,itg}",B) ;
        input("node_s_b{n,itg}") ;
        output("priority::priority::sStar_B{n,itg}") ;
        constraint("boundaryDisplacement{n,itg}") ;
      }
      InitializeDisplacementResidual() {
        name_store("sResidualTemp",sResidualTemp) ;
        output("sResidualTemp") ;
        constraint("nodes") ;
      }
      ComputeDisplacementResidualOne() {
        name_store("sStar_D",D) ;
        name_store("node_sStar",node_s) ;
        name_store("sStar_B",B) ;
        name_store("sResidualTemp",sResidualTemp) ;
        input("sStar_D,node_sStar,sStar_B") ;
        output("sResidualTemp") ;
        constraint("nodes") ;
      }
      ComputeDisplacementResidualTwo() {
        name_store("edge2node",edge2node) ;
        name_store("node_sStar",node_s) ;
        name_store("sStar_E",sStar_E) ;
        name_store("sResidualTemp",sResidualTemp) ;
        input("edge2node->node_sStar,sStar_E") ;
        output("edge2node->sResidualTemp") ;
        constraint("edges") ;
      }
      ComputeFinalDisplacementResidualDefault() {
        name_store("sResidualTemp",sResidualTemp) ;
        name_store("sResidual",sResidual) ;
        input("sResidualTemp") ;
        output("sResidual") ;
        constraint("nodes") ;
      }
      ComputeFinalDisplacementResidualSXZero() {
        name_store("sResidualTemp",sResidualTemp) ;
        name_store("priority::sResidual",sResidual) ;
        input("sResidualTemp") ;
        output("priority::sResidual") ;
        constraint("sXZeroNodes") ;
      }
      ComputeFinalDisplacementResidualSYZero() {
        name_store("sResidualTemp",sResidualTemp) ;
        name_store("priority::sResidual",sResidual) ;
        input("sResidualTemp") ;
        output("priority::sResidual") ;
        constraint("sYZeroNodes") ;
      }
      ComputeFinalDisplacementResidualSZZero() {
        name_store("sResidualTemp",sResidualTemp) ;
        name_store("priority::sResidual",sResidual) ;
        input("sResidualTemp") ;
        output("priority::sResidual") ;
        constraint("sZZeroNodes") ;
      }
      ComputeFinalDisplacementResidualBoundaryDisplacement() {
        name_store("priority::priority::sResidual",sResidual) ;
        output("priority::priority::sResidual") ;
        constraint("boundaryDisplacement") ;
      }
      InitializeTotalDisplacementResidual() {
        name_store("sResidualData",sResidualData) ;
        output("sResidualData") ;
        constraint("nodes") ;
      }
      ComputeTotalDisplacementResidual() {
        name_store("sResidual",sResidual) ;
        name_store("pos",pos) ;
        name_store("sResidualData",sResidualData) ;
        input("sResidual,pos") ;
        output("sResidualData") ;
        constraint("nodes") ;
      }
      NodeIterationFinishedBuild() {
        name_store("nodeIterationFinished{n,itg=-1}",nodeIterationFinished) ;
        output("nodeIterationFinished{n,itg=-1}") ;
//      constraint("UNIVERSE,gridMover,gridMotionTimeDependent") ;
        constraint("UNIVERSE{n},gridMover{n},gridMotionTimeDependent{n}") ;
      }
      CheckNodeIterationFinished() {
        name_store("$itg{n,itg}",itg) ;
        name_store("gridMoverMaxIterationsPerTimeStep{n,itg}",
          gridMoverMaxIterationsPerTimeStep) ;
        name_store("convergenceTolerance{n,itg}",convergenceTolerance) ;
        name_store("sResidualData{n,itg}",sResidualData) ;
        name_store("nodeIterationFinished{n,itg}",nodeIterationFinished) ;
        input("$itg{n,itg},gridMoverMaxIterationsPerTimeStep{n,itg}") ;
        input("convergenceTolerance{n,itg},sResidualData{n,itg}") ;
        output("nodeIterationFinished{n,itg}") ;
        constraint("sResidualData{n,itg},gridMover{n,itg}") ;
        constraint("gridMotionTimeDependent{n,itg}") ;
      }
