// Loci includes.
#include <Loci.h>
#include <Tools/tools.h>
using Loci::Area ;

// StreamUns includes.
#include "const.h"
#include "sciTypes.h"

namespace streamUns {

//-----------------------------------------------------------------------------
// Rules for some common geometric factors.

  class GeometryFactor0 : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<vect3d> cellCenter ;
      const_store<Area> area ;
      const_store<real> diffusionProduct ;
      store<vect3d> geometryFactor0 ;
    public:
                                                                                
      // Define input and output.
      GeometryFactor0() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("cellcenter",cellCenter) ;
        name_store("area",area) ;
        name_store("diffusionProduct",diffusionProduct) ;
        name_store("geometryFactor0",geometryFactor0) ;
        input("(cl,cr)->cellcenter,area,diffusionProduct") ;
        output("geometryFactor0") ;
        constraint("internalFaces") ;
      }

      void calculate(Entity face) {
        geometryFactor0[face]=area[face].n*area[face].sada-
          diffusionProduct[face]*(cellCenter[cr[face]]-cellCenter[cl[face]]) ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ;}
  } ;

  register_rule<GeometryFactor0> registerGeometryFactor0 ;

//-----------------------------------------------------------------------------
// Rules to compute radius at cell and face centers for axisymmetric geometries.

  class CellRadiusDefault : public pointwise_rule {
    private:
      store<real> cellRadius ;
    public:
                                                                                
      // Define input and output.
      CellRadiusDefault() {
        name_store("cellRadius",cellRadius) ;
        output("cellRadius") ;
        constraint("UNIVERSE") ;
      }

      inline void calculate(Entity cell) { cellRadius[cell]=1.0 ; }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ;}
  } ;

  register_rule<CellRadiusDefault> registerCellRadiusDefault ;

  class CellRadiusAxisymmetric : public pointwise_rule {
    private:
      const_store<vect3d> cellCenter ;
      store<real> cellRadius ;
    public:
                                                                                
      // Define input and output.
      CellRadiusAxisymmetric() {
        name_store("cellcenter",cellCenter) ;
        name_store("axisymmetric::cellRadius",cellRadius) ;
        input("cellcenter") ;
        output("axisymmetric::cellRadius") ;
        constraint("axisymmetric,geom_cells") ;
      }

      inline void calculate(Entity cell) {
        cellRadius[cell]=cellCenter[cell].y ;
      }

      // Loop over cells.
      void compute(const sequence &seq) { do_loop(seq,this) ;}
  } ;

  register_rule<CellRadiusAxisymmetric> registerCellRadiusAxisymmetric ;

  class FaceRadiusDefault : public pointwise_rule {
    private:
      store<real> faceRadius ;
    public:
                                                                                
      // Define input and output.
      FaceRadiusDefault() {
        name_store("faceRadius",faceRadius) ;
        output("faceRadius") ;
        constraint("UNIVERSE") ;
      }

      inline void calculate(Entity face) { faceRadius[face]=1.0 ; }

      // Loop over faces.
      void compute(const sequence &seq) { do_loop(seq,this) ;}
  } ;

  register_rule<FaceRadiusDefault> registerFaceRadiusDefault ;

  class FaceRadiusAxisymmetric : public pointwise_rule {
    private:
      const_store<vect3d> faceCenter ;
      store<real> faceRadius ;
    public:
                                                                                
      // Define input and output.
      FaceRadiusAxisymmetric() {
        name_store("facecenter",faceCenter) ;
        name_store("axisymmetric::faceRadius",faceRadius) ;
        input("facecenter") ;
        output("axisymmetric::faceRadius") ;
        constraint("axisymmetric,faces") ;
      }

      inline void calculate(Entity face) {
        faceRadius[face]=faceCenter[face].y ;
      }

      // Loop over faces.
      void compute(const sequence &seq) { do_loop(seq,this) ;}
  } ;

  register_rule<FaceRadiusAxisymmetric> registerFaceRadiusAxisymmetric ;

//-----------------------------------------------------------------------------
// Variable area_ic. For non-moving grids this is just the same as area. For
// moving grids, this rule is overridden by priority rules. Tried this as a
// rename rule, but causes Loci to hang.
                                                                                
  class InitialArea : public pointwise_rule {
    private:
      const_store<Area> area ;
      store<Area> area_ic ;
    public:
                                                                                
      // Define input and output.
      InitialArea() {
        name_store("area",area) ;
        name_store("area_ic",area_ic) ;
        input("area") ;
        output("area_ic") ;
      }

      inline void calculate(Entity face) { area_ic[face]=area[face] ; }

      // Empty compute method.
      void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<InitialArea> registerInitialArea ;

}
