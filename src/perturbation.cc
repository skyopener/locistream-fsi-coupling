//-----------------------------------------------------------------------------
// Description: This file contains rules for a perturbation point, which is
//   used normally to trip a flow into asymmetry.
//-----------------------------------------------------------------------------
                                                                                
// Loci includes.
#include <Loci.h>

// StreamUns includes.
#include "sciTypes.h"
#include "varsFileInputs.h"

namespace streamUns {

//-----------------------------------------------------------------------------
// Rules to process ignition options from the .vars file.

  // Creates the ignition parameters.
  class PerturbationParameters : public singleton_rule {
    private:
      const_param<PerturbationOptions> perturbation ;
      const_param<real> Lref ;
      param<vect3d> perturbation_dv ;
      param<vect3d> perturbation_pos ;
      param<vect3d> perturbation_delta ;
      param<int> perturbation_nstart ;
      param<int> perturbation_nstop ;
    public:
                                                                                
      // Define input and output.
      PerturbationParameters() {
        name_store("perturbation",perturbation) ;
        name_store("Lref",Lref) ;
        name_store("perturbation_dv",perturbation_dv) ;
        name_store("perturbation_pos",perturbation_pos) ;
        name_store("perturbation_delta",perturbation_delta) ;
        name_store("perturbation_nstart",perturbation_nstart) ;
        name_store("perturbation_nstop",perturbation_nstop) ;
        input("perturbation,Lref") ;
        output("perturbation_dv,perturbation_pos,perturbation_delta") ;
        output("perturbation_nstart,perturbation_nstop") ;
      }
                                                                                
      // Set up the parameters.
      virtual void compute(const sequence& seq) {

        if((*perturbation).optionExists("position")){

          // Velocity perturbation.
          Loci::option_value_type ovt=perturbation->getOptionValueType("dv") ;
          vect3d val(0.0,0.0,0.0) ;
          if(ovt==Loci::LIST){
            Loci::options_list::arg_list valueList ;
            perturbation->getOption("dv",valueList) ;
            if(valueList.size()!=3){
              cerr << "Error on reading perturbation dv: vector input must "
                << "contain 3 terms." << endl ;
            }
            for(int j=0;j<3;++j)
              if(valueList[j].type_of()!=Loci::REAL){
                cerr << "Improper vector specification for perturbation dv."
                  << endl ;
              }
            double vecVal[3] ;
            for(int j=0;j<3;++j){ valueList[j].get_value(vecVal[j]) ; }
            val=vect3d(vecVal[0],vecVal[1],vecVal[2]) ;
          }else{
            cerr << "Improper value for perturbation dv." << endl ;
          }
          *perturbation_dv=val ;

          // Position.
          ovt=perturbation->getOptionValueType("position") ;
          if(ovt==Loci::LIST){
            Loci::options_list::arg_list valueList ;
            perturbation->getOption("position",valueList) ;
            if(valueList.size()!=3){
              cerr << "Error on reading perturbation position: vector input "
                << "must contain 3 terms." << endl ;
            }
            for(int j=0;j<3;++j)
              if(valueList[j].type_of()!=Loci::REAL){
                cerr << "Improper vector specification for perturbation "
                  << "position." << endl ;
              }
            double vecVal[3] ;
            for(int j=0;j<3;++j){ valueList[j].get_value(vecVal[j]) ; }
            val=vect3d(vecVal[0],vecVal[1],vecVal[2]) ;
          }else{
            cerr << "Improper coordinate for perturbation position." << endl ;
          }
          *perturbation_pos=*Lref*val ;

          // Delta.
          val=vect3d(0,0,0) ; ovt=perturbation->getOptionValueType("delta") ;
          if(ovt==Loci::LIST){
            Loci::options_list::arg_list valueList ;
            perturbation->getOption("delta",valueList) ;
            if(valueList.size()!=3) {
              cerr << "Error on reading perturbation delta: vector input must "
                "contain 3 terms." << endl ;
            }
            for(int j=0;j<3;++j)
              if(valueList[j].type_of()!=Loci::REAL){
                cerr << "Improper vector specification for perturbation delta."
                  << endl ;
              }
            double vecVal[3] ;
            for(int j=0;j<3;++j){ valueList[j].get_value(vecVal[j]) ; }
            val=vect3d(vecVal[0],vecVal[1],vecVal[2]) ;
          }else if(ovt == Loci::REAL){
            real v1 ; perturbation->getOption("delta",v1) ;
            val=vect3d(v1,v1,v1) ;
          }else{
            cerr << "Improper coordinate for perturbation delta." << endl ;
          }
          *perturbation_delta=*Lref*val ;

          // Starting and ending time step.
          double v1=0 ;
          if(perturbation->optionExists("nstart")){
            perturbation->getOption("nstart",v1) ;
            *perturbation_nstart=int(v1) ;
          }
          v1=0 ;
          if(perturbation->optionExists("nstop")){
            perturbation->getOption("nstop",v1) ;
            *perturbation_nstop = int(v1) ;
          }
        }
      }
  } ;
                                                                                
  register_rule<PerturbationParameters> registerPerturbationParamters ;

//-----------------------------------------------------------------------------
// Rule to add a perturbation to the corrected velocity.

  class CorrectedVelocityFromPerturbation : public apply_rule<store<vect3d>,
  Loci::Summation<vect3d> > {
    private:
      const_param<int> perturbation_nstart,perturbation_nstop,nCycle ;
      const_param<vect3d> perturbation_dv ;
      const_param<vect3d> perturbation_pos ;
      const_param<vect3d> perturbation_delta ;
      const_store<vect3d> cellCenter ;
      const_store<real> vol ;
      store<vect3d> vCorrected ;
    public:

      // Define input and output.
      CorrectedVelocityFromPerturbation() {
        name_store("perturbation_nstart{n,it}",perturbation_nstart) ;
        name_store("perturbation_nstop{n,it}",perturbation_nstop) ;
        name_store("ncycle{n}",nCycle) ;
        name_store("perturbation_dv{n,it}",perturbation_dv) ;
        name_store("perturbation_pos{n,it}",perturbation_pos) ;
        name_store("perturbation_delta{n,it}",perturbation_delta) ;
        name_store("cellcenter{n,it}",cellCenter) ;
        name_store("vol{n,it}",vol) ;
        name_store("vCorrected{n,it}",vCorrected) ;
        input("ncycle{n}") ;
        input("perturbation_nstart{n,it},perturbation_nstop{n,it}") ;
        input("perturbation_dv{n,it},perturbation_pos{n,it}") ;
        input("perturbation_delta{n,it},cellcenter{n,it},vol{n,it}") ;
        output("vCorrected{n,it}") ;
        constraint("geom_cells{n,it}") ;
      }

      // Add random perturbation up to dv in magnitude. The sign is
      // determined by the oddness or evenness of rand().
      void calculate(Entity cell) {
        real ref_length=pow(vol[cell],0.33333333) ;
        vect3d dp=cellCenter[cell]-perturbation_pos[cell] ;
        if(fabs(dp.x)<(perturbation_delta[cell].x+ref_length) &&
        fabs(dp.y)<(perturbation_delta[cell].y+ref_length) &&
        fabs(dp.z)<(perturbation_delta[cell].z+ref_length)) {
          vect3d dV=(real(rand())/real(RAND_MAX))*(*perturbation_dv) ;
          vCorrected[cell]+=((rand()%2)? -1.0*dV:dV) ;
        }
      }

      // Loop over cells.
      void compute(const sequence &seq) {
        if(*nCycle >= *perturbation_nstart && *nCycle <= *perturbation_nstop){
          if(Loci::MPI_rank==0) cout << "    INFO: Perturbation active"
            << endl ;
          do_loop(seq,this) ;
        }
      }
  } ;

  register_rule<CorrectedVelocityFromPerturbation>
    registerCorrectedVelocityFromPerturbation ;

}
