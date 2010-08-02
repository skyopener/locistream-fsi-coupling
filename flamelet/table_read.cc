// Standard library includes.
#include <vector>
using std::vector ;

// Loci includes.      
#include <Loci.h>
using Loci::Area ;

// StreamUns includes.                              
#include "sciTypes.h"       

#include "flamelet_table.h"
extern Flamelet_Table flamelet_table;

namespace streamUns {       

/*==============================DENSITY=====================================*/
  // Rule for reading the density off the flamelet table for interior cells
  class FlameletRhoInterior : public pointwise_rule {
    private:
      const_store<real> Z,Zvar,Chi ;
      store<real> rho ;
    public:

      // Define input and output.
      FlameletRhoInterior() {
        name_store("flamelet::rhoStar{n,it}",rho) ;
        name_store("ZStar{n,it}",Z) ;
        name_store("ZvarStar{n,it}",Zvar) ;
        name_store("Chi{n,it}",Chi) ;
        input("ZStar{n,it},ZvarStar{n,it},Chi{n,it}") ;
        output("flamelet::rhoStar{n,it}") ;
        constraint("geom_cells,flameletModel") ;
      }

      void calculate(Entity cell) {
      	rho[cell]=flamelet_table.get_rho(Z[cell],Zvar[cell],Chi[cell],true);
      }

      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;      

  register_rule<FlameletRhoInterior>
    registerFlameletRhoInterior ;   

  // Rule for reading the density off the flamelet table for boundary faces
  class FlameletRhoBoundary : public pointwise_rule {
    private:
      const_store<real> Z_f ;
      const_store<real> Zvar_f ;
      const_store<real> Chi_f ;
      store<real> rho_f ;
    public:

      // Define input and output.
      FlameletRhoBoundary() {       
        name_store("flamelet::extrapolated::rho_f",rho_f) ;
        name_store("Z_f",Z_f) ;
        name_store("Zvar_f",Zvar_f) ;
        name_store("Chi_f",Chi_f) ;
        input("Z_f,Zvar_f,Chi_f") ;
        output("flamelet::extrapolated::rho_f") ;
        constraint("boundaryFaces,flameletModel") ;
      }

      void calculate(Entity face) {
      	rho_f[face]=flamelet_table.get_rho(Z_f[face],Zvar_f[face],Chi_f[face],true);
      }

      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FlameletRhoBoundary>
    registerFlameletRhoBoundary ;
  
/*==============================TEMPERATURE=====================================*/
  // Rule for reading the temperature off the flamelet table for interior cells
  class FlameletTemperatureInterior : public pointwise_rule {
    private:
      const_store<real> Z,Zvar,Chi ;
      store<real> temperature ;
    public:

      // Define input and output.
      FlameletTemperatureInterior() {
        name_store("flamelet::temperature",temperature) ;
        name_store("Z",Z) ;
        name_store("Zvar",Zvar) ;
        name_store("Chi",Chi) ;
        input("Z,Zvar,Chi") ;
        output("flamelet::temperature") ;
        constraint("geom_cells,flameletModel") ;
      }

      void calculate(Entity cell) {
      	temperature[cell]=flamelet_table.get_temperature(Z[cell],Zvar[cell],Chi[cell],true);
      }

      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;      

  register_rule<FlameletTemperatureInterior>
    registerFlameletTemperatureInterior ;   

  // Rule for reading the temperature off the flamelet table for boundary faces
  class FlameletTemperatureBoundary : public pointwise_rule {
    private:
      const_store<real> Z_f ;
      const_store<real> Zvar_f ;
      const_store<real> Chi_f ;
      store<real> temperature_f ;
    public:

      // Define input and output.
      FlameletTemperatureBoundary() {       
        name_store("flamelet::temperature_f",temperature_f) ;
        name_store("Z_f",Z_f) ;
        name_store("Zvar_f",Zvar_f) ;
        name_store("Chi_f",Chi_f) ;
        input("Z_f,Zvar_f,Chi_f") ;
        output("flamelet::temperature_f") ;
        constraint("boundaryFaces,flameletModel") ;
      }

      void calculate(Entity face) {
      	temperature_f[face]=flamelet_table.get_temperature(Z_f[face],Zvar_f[face],Chi_f[face],true);
      }

      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FlameletTemperatureBoundary>
    registerFlameletTemperatureBoundary ;  
  

/*==============================LaminarViscosity=====================================*/
  // Rule for reading the laminarViscosity off the flamelet table for interior cells
  class FlameletLaminarViscosityInterior : public pointwise_rule {
    private:
      const_store<real> Z,Zvar,Chi ;
      store<real> laminarViscosity ;
    public:

      // Define input and output.
      FlameletLaminarViscosityInterior() {
        name_store("flamelet::laminarViscosity",laminarViscosity) ;
        name_store("Z",Z) ;
        name_store("Zvar",Zvar) ;
        name_store("Chi",Chi) ;
        input("Z,Zvar,Chi") ;
        output("flamelet::laminarViscosity") ;
        constraint("geom_cells,flameletModel") ;
      }

      void calculate(Entity cell) {
      	laminarViscosity[cell]=flamelet_table.get_laminarViscosity(Z[cell],Zvar[cell],Chi[cell],true);
      }

      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;      

  register_rule<FlameletLaminarViscosityInterior>
    registerFlameletLaminarViscosityInterior ;   

  // Rule for reading the laminarViscosity off the flamelet table for boundary faces
  class FlameletLaminarViscosityBoundary : public pointwise_rule {
    private:
      const_store<real> Z_f ;
      const_store<real> Zvar_f ;
      const_store<real> Chi_f ;
      store<real> laminarViscosity_f ;
    public:

      // Define input and output.
      FlameletLaminarViscosityBoundary() {       
        name_store("flamelet::laminarViscosity_f",laminarViscosity_f) ;
        name_store("Z_f",Z_f) ;
        name_store("Zvar_f",Zvar_f) ;
        name_store("Chi_f",Chi_f) ;
        input("Z_f,Zvar_f,Chi_f") ;
        output("flamelet::laminarViscosity_f") ;
        constraint("boundaryFaces,flameletModel") ;
      }

      void calculate(Entity face) {
      	laminarViscosity_f[face]=flamelet_table.get_laminarViscosity(Z_f[face],Zvar_f[face],Chi_f[face],true);
      }

      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FlameletLaminarViscosityBoundary>
    registerFlameletLaminarViscosityBoundary ;  

/*==============================Diffusivity=====================================*/
  // Rule for reading the diffusivity off the flamelet table for interior cells
  class FlameletDiffusivityInterior : public pointwise_rule {
    private:
      const_store<real> Z,Zvar,Chi ;
      store<real> diffusivity ;
    public:

      // Define input and output.
      FlameletDiffusivityInterior() {
        name_store("flamelet::diffusivity",diffusivity) ;
        name_store("Z",Z) ;
        name_store("Zvar",Zvar) ;
        name_store("Chi",Chi) ;
        input("Z,Zvar,Chi") ;
        output("flamelet::diffusivity") ;
        constraint("geom_cells,flameletModel") ;
      }

      void calculate(Entity cell) {
      	diffusivity[cell]=flamelet_table.get_diffusivity(Z[cell],Zvar[cell],Chi[cell],true);
      }

      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;      

  register_rule<FlameletDiffusivityInterior>
    registerFlameletDiffusivityInterior ;   

  // Rule for reading the diffusivity off the flamelet table for boundary faces
  class FlameletDiffusivityBoundary : public pointwise_rule {
    private:
      const_store<real> Z_f ;
      const_store<real> Zvar_f ;
      const_store<real> Chi_f ;
      store<real> diffusivity_f ;
    public:

      // Define input and output.
      FlameletDiffusivityBoundary() {       
        name_store("flamelet::diffusivity_f",diffusivity_f) ;
        name_store("Z_f",Z_f) ;
        name_store("Zvar_f",Zvar_f) ;
        name_store("Chi_f",Chi_f) ;
        input("Z_f,Zvar_f,Chi_f") ;
        output("flamelet::diffusivity_f") ;
        constraint("boundaryFaces,flameletModel") ;
      }

      void calculate(Entity face) {
      	diffusivity_f[face]=flamelet_table.get_diffusivity(Z_f[face],Zvar_f[face],Chi_f[face],true);
      }

      virtual void compute(const sequence &seq) { do_loop(seq,this) ; }
  } ;

  register_rule<FlameletDiffusivityBoundary>
    registerFlameletDiffusivityBoundary ;  

}
