// Standard library includes.
#include <fstream>

// Loci includes.
#include <Loci.h>
using Loci::Area ;

// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;
    
// StreamUns includes.
#include "residual.h"
#include "sciTypes.h"
#include "varsFileInputs.h"
       
namespace streamUns {

  class InteriorFaceDensityCavitationFOU : public pointwise_rule {
    private:
      const_Map cl,cr ;
      const_store<real> rho ;
      const_store<vect3d> v ;
      const_store<real> Alpha ;
      const_store<Area> area ;
      store<real> faceDensity ;
    public:
 
       InteriorFaceDensityCavitationFOU() {
        name_store("cl",cl) ;
        name_store("cr",cr) ;
        name_store("rho",rho) ;
        name_store("v",v) ;
        name_store("Alpha",Alpha) ;
        name_store("area",area) ;
        name_store("cavitation::faceDensity",faceDensity) ;
        input("(cl,cr)->(rho,v,Alpha)") ;
        input("area") ;
        output("cavitation::faceDensity") ;
        constraint("internalFaces") ;
      }

      void calculate(Entity face) {
        const real AlphaFace=0.5*(Alpha[cl[face]]+Alpha[cr[face]]) ;
        if(AlphaFace<1.0){
          vect3d faceVelocity=0.5*(v[cl[face]]+v[cr[face]]) ;
          faceDensity[face]=(dot(faceVelocity,area[face].n)>0.0)? rho[cl[face]]:
            rho[cr[face]] ;
        }else{
          faceDensity[face]=0.5*(rho[cl[face]]+rho[cr[face]]) ;
        }
      }
 
      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
 
  } ;
 
  register_rule<InteriorFaceDensityCavitationFOU>
    registerInteriorFaceDensityCavitationFOU ;

}
