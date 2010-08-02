#include "reaction.h"
#include <Tools/tools.h>
#include "storeVec.h"
#include "fluidConst.h"


using std::map ;
using std::pair ;

using std::vector ;
using std::sort ;
using std::string ;
using std::cerr ;
using std::endl ;
// Limit on exp evaluations not to exceed
#define EXP_LIM (512.0)
// Limit on reciprocal of Kc (To prevent div0)
#define KCR_EPS (1e-307)

namespace fluidPhysics {

  curveFitFunctions::curveFitFunctions() {
    num_power = 0 ;
    num_exponential = 0 ;
    thindex(0.0) ;     //theta[0]
    etaindex(0.0) ;    //eta[0]
    etaindex(1.0) ;    //eta[1]
    etaindex(-1.0) ;   //eta[2]
  }

  // the following function provides the index with unrepeatable
  // (Cvar, etavar, thetavar), which is the preparation for the calculate
  // function to compute f_val
  int curveFitFunctions::getArrheniusidx(double Cvar, double etavar, double thetavar,
                                         int &num_func_results )
  {
    func_rep Ar(Cvar,etaindex(etavar),thindex(thetavar)) ;
    map<func_rep,int>::const_iterator fi ;
    if((fi = func_map.find(Ar)) != func_map.end()) {
      return fi->second ;
    }
    int idx = func_list.size() ;
    func_list.push_back(Ar) ; //func_list in which Ar is unrepeatable
    res_loc.push_back(num_func_results) ;
    num_func_results++ ;
    func_map[Ar] = res_loc[idx] ;
    return res_loc[idx] ;
  }

  // the following function calculates Kf/Kc and derivative of Kf/Kc
  //Note: the first value in the pair &f_val represents Kf/Kc, and second one
  //represents the derivative of Kf/Kc
  void curveFitFunctions::calculate(double T, pair<double,double> *f_val) const {
    const double Tr = 1.0/T ;
    const double Tr2 = Tr*Tr ;

    scratch_array<pair<double,double> > power(num_power) ;
    scratch_array<pair<double,double> > exponential(num_exponential) ;

    // Compute powers
    power[0] = pair<double,double>(1.0,0.0) ;
    exponential[0] = pair<double,double>(1.0,0.0) ;
    power[1] = pair<double,double>(T,1.0) ;
    power[2] = pair<double,double>(Tr,-Tr2) ;
    for(size_t i=3;i<eta.size();++i) {
      const double et = eta[i] ;
      const double p = pow(T,et) ;
      power[i] = pair<double,double>(p,et*Tr*p) ;
    }
    // Compute exponentials
    for(size_t i=1;i<theta.size();++i) {
      const double th = theta[i] ;
      const double e = exp(min(-th*Tr,EXP_LIM)) ;
      exponential[i] = pair<double,double>(e,th*Tr2*e) ;
    }
    // Compute functions + function derivatives
    for(size_t i=0;i<func_list.size();++i) {
      const double Cv = func_list[i].C ;
      const double pw = power[func_list[i].etai].first ;
      const double pwp = power[func_list[i].etai].second ;
      const double e = exponential[func_list[i].thetai].first ;
      const double ep = exponential[func_list[i].thetai].second ;
      f_val[res_loc[i]] = pair<double,double>(Cv*pw*e,Cv*(pw*ep + pwp*e)) ;
    }

  }

  int get_function(string var, const reaction_options ro,
                   curveFitFunctions &cfm,
                   int &fnum)
  {
    if(ro.getOptionValueType(var) != Loci::FUNCTION) {
      cerr << var << " has uninterpretable function" << endl ;
      return cfm.compile_Arrhenius(0.0,0.0,0.0,fnum) ;
    }
    string fname ;
    options_list::arg_list alist ;  //alist is (C, eta, theta) list in database,
    //                              //they may be repeatable
    ro.getOption(var,fname,alist) ;

    if(fname != string("Arrhenius")) {
      cerr << fname << ": unknown function type" << endl ;
      return cfm.compile_Arrhenius(0.0,0.0,0.0,fnum) ;
    }
    if(alist.size() != 3) {
      cerr << "Arrhenius function only has three arguments" << endl ;
      return cfm.compile_Arrhenius(0.0,0.0,0.0,fnum) ;
    }
    options_list::arg_list::iterator ali ;
    for(ali=alist.begin();ali!=alist.end();++ali)
      if((*ali).type_of() != Loci::REAL) {
        cerr << "Arrhenius function takes three REAL arguments" << endl ;
        return cfm.compile_Arrhenius(0.0,0.0,0.0,fnum) ;
      }
    ali = alist.begin() ;
    double C ;
    (*ali).get_value(C) ;  //get C value from alist
    ++ali ;
    double eta ;
    (*ali).get_value(eta) ;  //get eta value from alist
    ++ali ;
    double theta ;
    (*ali).get_value(theta) ; //get theta value from alist
    return cfm.compile_Arrhenius(C,eta,theta,fnum) ; //By calling function
    //Arrhenius, wet get unrepeatable (C, eta, theta) list
  }

  thermoKc::ri::ri(const EOS &eos, const reaction_expression &r)
  {
    const reaction_expression::coef_list &cldiff = r.productsMinusReactants() ;
    for(size_t i=0;i<cldiff.size();++i) {
      si sinfo ;
      sinfo.species = eos.speciesIndex(cldiff[i].speciesName) ; //species index
      if(sinfo.species < 0) {
        if(Loci::MPI_rank == 0)
          cerr << "FATAL ERROR: species " << cldiff[i].speciesName << " not found in EoS!" << endl ;
        Loci::Abort() ;
      }
      //                                              // in the database
      sinfo.nu_diff = cldiff[i].nu ; //productsMinusReactants coefficient for
      //                             //each species
      if(sinfo.nu_diff>0)
        pos.push_back(sinfo) ;   //pos(species_index, nu_diff)
      if(sinfo.nu_diff<0)
        neg.push_back(sinfo) ;
    }
    sort(pos.begin(),pos.end()) ; //sorting based on species index
    sort(neg.begin(),neg.end()) ;

    forwardonly = r.isforwardonly() ;
  }

  bool thermoKc::ri::operator==(const ri &r)
  {
    return pos == r.pos && neg == r.neg ;
  }

  void thermoKc::initialize(const EOS &eos)
  {
    int ns = eos.numSpecies() ;
    //  Omega = new double[ns] ;
    //  Omegap = new double[ns] ;
    for(int i=0;i<ns;++i) {
      spi s ;
      s.PrefR = eos.speciesPref(i)/Rh ;
      spiv.push_back(s) ;
    }
  }

  //The following function builds up reaction information list.
  //This is the preparation for function calculate to compute Kc
  int thermoKc::compile_KcFunc(const EOS &eos, const reaction_expression &r,
                               int &fnum)
  {
    ri reaction(eos,r) ; //one specific reaction, already sorted
    for(size_t i=0;i<reactinfo.size();++i)
      if(reactinfo[i] == reaction) { //this can leave out M body reaction
        //as one reaction
        return res_loc[i] ;
      }
    int rnum = reactinfo.size() ;
    reactinfo.push_back(reaction) ;
    res_loc.push_back(fnum) ;
    fnum++ ;
    for(size_t i=0;i<reaction.pos.size();++i)
      speciesSet += reaction.pos[i].species ;
    for(size_t i=0;i<reaction.neg.size();++i)
      speciesSet += reaction.neg[i].species ;
    return res_loc[rnum] ;
  }

  //The following function calculate Kc, dKc/dT using Gibbs Free energy
  //minimization
  void thermoKc::calculate(const EOS &eos, double T,
                           pair<double,double> *f_val)  const {
    int ns = eos.numSpecies() ;
    scratch_array<double> Omega(ns), Omegap(ns) ;

    eos.getOmegas(T,Omega,Omegap) ;


    for(size_t r=0;r<reactinfo.size();++r) {
      const ri &rinfo = reactinfo[r] ;
      if(rinfo.forwardonly) {
        f_val[res_loc[r]].first = 0 ;
        f_val[res_loc[r]].second = 0 ;
        continue ;
      }
      double prod1=1.0, prod2=1.0, sum=0.0, sump=0.0 ;
      double nusum=0 ;
      for(size_t i=0;i<rinfo.pos.size();++i) {
        const int ss = rinfo.pos[i].species ;
        const double nu_diff = rinfo.pos[i].nu_diff ;

        double pp = spiv[ss].PrefR ;
        nusum+=nu_diff ;
        prod1 *= pow(max(pp,0.0),nu_diff) ;
        sum +=  Omega[ss]*nu_diff ;
        sump += Omegap[ss]*nu_diff ;
      }
      for(size_t i=0;i<rinfo.neg.size();++i) {
        const int ss = rinfo.neg[i].species ;
        const double nu_diff = rinfo.neg[i].nu_diff ;
        double pp = spiv[ss].PrefR ;
        nusum+=nu_diff ;
        prod2 *= pow(max(pp,0.0),-nu_diff) ;
        sum +=  Omega[ss]*nu_diff ;
        sump += Omegap[ss]*nu_diff ;
      }

      const double S = prod1/prod2 ;
      const double Kc = pow(T,-nusum)*S*exp(min(-sum,EXP_LIM)) ;


      f_val[res_loc[r]].first = Kc ;
      f_val[res_loc[r]].second = (-double(nusum)/T - sump)*Kc ;
    }
  }

  vector<double> &get_MBweights(const EOS &eos, reaction_options &ro) {
    int num_species = eos.numSpecies() ;
    vector<double> *wts = new vector<double>(num_species) ;
    for(int i=0;i<num_species;++i)
      (*wts)[i] = 1.0 ;

    if(ro.getOptionValueType("MB") != Loci::LIST) {
      cerr << "M-Body species weights should be assigned to list" << endl ;
      return *wts ;
    } else {
      Loci::options_list::arg_list species_list ;
      ro.getOption("MB",species_list) ;
      Loci::options_list::arg_list::iterator li ;
      for(li=species_list.begin();li!=species_list.end();++li) {
        Loci::option_values::value_list_type species_arg ;
        li->get_value(species_arg) ;
        if(li->type_of() != Loci::NAME_ASSIGN || species_arg.size() != 1
           || species_arg.front().type_of() != Loci::REAL) {
          cerr << "error in MBody weights assignment" << endl ;
          Loci::Abort() ;
        } else {
          string sn ;
          double sv ;
          li->get_value(sn) ;
          species_arg.front().get_value(sv) ;
          int sid = eos.speciesIndex(sn) ;
          if(sid < 0) {
            if(Loci::MPI_rank == 0)
              cerr << "WARNING: species " << sn << " not in EoS but used in Mbody Weights Expression" << endl ;
          } else
            (*wts)[sid] = sv ;
        }
      }
    }
    return *wts ;
  }

  vector<double> get_exp_nu(const EOS &eos, reaction_options &ro) {
    int num_species = eos.numSpecies() ;
    vector<double> exp_nu(num_species) ;
    for(int i=0;i<num_species;++i)
      exp_nu[i] = 0.0 ;

    if(ro.getOptionValueType("exp_nu") != Loci::LIST) {
      cerr << "Exponential nu in reactants should be assigned to list" << endl ;
      return exp_nu ;
    } else {
      Loci::options_list::arg_list species_list ;
      ro.getOption("exp_nu",species_list) ;
      Loci::options_list::arg_list::iterator li ;
      for(li=species_list.begin();li!=species_list.end();++li) {
        Loci::option_values::value_list_type species_arg ;
        li->get_value(species_arg) ;
        if(li->type_of() != Loci::NAME_ASSIGN || species_arg.size() != 1
           || species_arg.front().type_of() != Loci::REAL) {
          cerr << "error in exponential nu assignment" << endl ;
          Loci::Abort() ;
        }
        else {
          string sn ;
          double sv ;
          li->get_value(sn) ;
          species_arg.front().get_value(sv) ;
          int sid = eos.speciesIndex(sn) ;
          if(sid < 0) {
            if(Loci::MPI_rank == 0)
              cerr << "WARNING: Invalid species name " << sn << " in exp_nu expression"
                   << endl ;
          } else
            exp_nu[sid] = sv ;
        }
      }
    }
    return exp_nu ;
  }



  void reaction::initialize(const EOS &eos, const reaction_db &reactions,
                            const species_db  &species)
  {
    thermoKcModule.initialize(eos) ;

    num_species = eos.numSpecies() ;
    int num_elements = species.numElements() ;

    pmr = new vector<scdiff>[num_species] ;
    num_reactions = reactions.numReactions() ;
    reactants = new vector<sc>[num_reactions] ;
    reactants_int = new vector<sc_int>[num_reactions] ;
    products = new vector<sc>[num_reactions] ;
    products_int = new vector<sc_int>[num_reactions] ;

    for(int i=0;i<num_reactions;++i) {
      vector<double> el(num_elements) ;
      vector<sc> reactants_i,products_i ;

      for(int j=0;j<num_elements;++j)
        el[j] = 0 ;
      const reaction_expression &r = reactions.getReactionStoichiometry(i) ;
      const reaction_expression::coef_list &clr = r.reactants() ;
      isforwardonly.push_back(r.isforwardonly()) ;
      bool bogus_species = false ;
      for(size_t j=0;j<clr.size();++j) {
        if(clr[j].speciesName[0] == '_')
          bogus_species = true ;

        int sp = eos.speciesIndex(clr[j].speciesName) ;
        if(sp<0) {
          if(Loci::MPI_rank == 0)
            cerr << "WARNING: species " << clr[j].speciesName
                 << " is defined in reaction, but is not in species list!"
                 << endl ;
          Loci::Abort() ;
        }

        double nu      = clr[j].nu ;
        //building up vector reactants[reaction](species, nu)
        reactants_i.push_back(sc(sp,nu)) ;
        const vector<int> &ev = species.getSpeciesElements(sp) ;
        for(int e=0;e<num_elements;++e)
          el[e] += nu*double(ev[e]) ;
      }
      const reaction_expression::coef_list &clp = r.products() ;
      for(size_t j=0;j<clp.size();++j) {
        if(clp[j].speciesName[0] == '_')
          bogus_species = true ;

        int sp = eos.speciesIndex(clp[j].speciesName) ;
        if(sp < 0) {
          if(Loci::MPI_rank == 0)
            cerr << "FATAL ERROR: species " << clp[j].speciesName << " in reaction mechanism but not in EoS." << endl ;
          Loci::Abort() ;
        }


        double nu      = clp[j].nu ;
        //building up vector products[reaction](species, nu)
        products_i.push_back(sc(sp,nu)) ;
        const vector<int> &ev = species.getSpeciesElements(sp) ;
        for(int e=0;e<num_elements;++e)
          el[e] -= nu*double(ev[e]) ;
      }

      if(!bogus_species)
        for(int e=0;e<num_elements;++e)
          if(fabs(el[e]) > EPSILON) {
            cerr << "warning: reaction is not stoichiometrically balanced!"
                 << endl  ;
            cerr << "reaction = " << r << endl ;
            break ;
          }

      const reaction_expression::coef_list &cldiff = r.productsMinusReactants() ;
      for(size_t j=0;j<cldiff.size();++j) {
        int species = eos.speciesIndex(cldiff[j].speciesName) ;

        if(species>=0 && species<num_species)
          //building up vector pmr[species](reaction, nu_diff)
          pmr[species].push_back(scdiff(i,cldiff[j].nu)) ;
        else {
          if(Loci::MPI_rank == 0)
            cerr << "FATAL ERROR: species " << cldiff[j].speciesName
                 << " in reaction mechanism but not in EoS" << endl ;
          Loci::Abort() ;
        }
      }


      reaction_options ro = reactions.getReactionOption(i) ;
      //building up vector Kfi (indexes to Kf), meanwhile by calling function,
      //get_function, func_list, res_loc and func_map func_resand are created
      //(all sorted)
      Kfi.push_back(get_function("Kf",ro,ArrheniusModule,num_func_results)) ;
      //      if(!isforwardonly[i]) {
      if(ro.optionExists("Kc")) {
        if(ro.getOptionValueType("Kc") == Loci::NAME) {
          if(!ro.checkOption("Kc","Thermodynamic")) {
            cerr << "Does not undertand Kc variable, database set Kc=" ;
            string name ;
            ro.getOption("Kc",name) ;
            cerr << name <<", expected Thermodynamic" << endl ;
            exit(1) ;
          }
          //building up vector Kci (indexes to Kc), meanwhile by calling function
          //compile_KcFunc (for Gibbs Free energy minimization) or get_function
          //(for curve-fit), containers for holding the information to compute
          //Kc are created (all sorted)
          Kci.push_back(thermoKcModule.compile_KcFunc(eos,r,num_func_results)) ;
        } else {
          Kci.push_back(get_function("Kc",ro,ArrheniusModule,num_func_results)) ;
        }
      } else {
        Kci.push_back(thermoKcModule.compile_KcFunc(eos,r,num_func_results)) ;
      }

      if(r.isMbodyReaction()) {
        isMB.push_back(true) ;
        if(ro.optionExists("MB")) {
          hasMBwts.push_back(true) ;
          MBwts.push_back(get_MBweights(eos,ro)) ;
        } else {
          std::vector<double> wts ;
          hasMBwts.push_back(false) ;
          MBwts.push_back(wts) ;
        }
      } else {
        isMB.push_back(false) ;
        std::vector<double> wts ;
        hasMBwts.push_back(false) ;
        MBwts.push_back(wts) ;
      }        

      double MinMF = 1e-10 ;
      if(ro.optionExists("MinMF")) {
        ro.getOption("MinMF",MinMF) ;
      }

//rate modifier ---------------------------------------------------------------
      if(ro.optionExists("rate_modifier")) {
        if(ro.getOptionValueType("rate_modifier") != Loci::FUNCTION) {
          cerr << "rate_modifier has uninterpretable function" << endl ;
          Loci::Abort() ;
        }
        string fname ;
        Loci::options_list::arg_list modifier_list ;
        ro.getOption("rate_modifier",fname,modifier_list) ;
        if(fname == "Rox") {
          Rox rx ;
          rx.Pindex = eos.speciesIndex("O2") ;
          if(rx.Pindex < 0) {
            if(Loci::MPI_rank == 0)
              cerr << "FATAL ERROR: Not able to find species O2 for Rox expression in reaction mechanism" << endl ;
            Loci::Abort() ;
          }

          Loci::options_list::arg_list::iterator li=modifier_list.begin() ;

          li->get_value(rx.KA1) ;
          ++li ;
          li->get_value(rx.KA2) ;
          ++li ;
          li->get_value(rx.KB1) ;
          ++li ;
          li->get_value(rx.KB2) ;
          ++li ;
          li->get_value(rx.KT1) ;
          ++li ;
          li->get_value(rx.KT2) ;
          ++li ;
          li->get_value(rx.KZ1) ;
          ++li ;
          li->get_value(rx.KZ2) ;
          ++li ;
          li->get_value(rx.rho_s) ;
          ++li ;
          li->get_value(rx.D_s) ;
          rx.reaction_index = i ;
          Rox_modifier.push_back(rx) ;
        } else if(fname == "pressure") {
          Loci::options_list::arg_list::iterator li=modifier_list.begin() ;
          double power_factor = 0.0 ;
          double Pref = 1e5 ;
          li->get_value(power_factor) ;
          ++li ;
          li->get_value(Pref) ;
          pressure_modifier.push_back(Prate(i,power_factor,Pref)) ;
//ST------------------------------------------------------------------
//--------------------------------------------------------------------
        } else if(fname == "condensation") {
          if(!r.isforwardonly()) {
            cerr << "condensation only valid for forward reactions!" << endl ;
            Loci::Abort() ;
          }
//Read metals list---------
          vector<double> metals(num_species,0.0) ;
          string validOptions="metals:oxides:Econd" ;
          Loci::options_list functionOptions(validOptions) ;
          functionOptions.Input(modifier_list) ;
          if(functionOptions.optionExists("metals")) {
            Loci::options_list::arg_list species_list ;
            functionOptions.getOption("metals",species_list) ;
            Loci::options_list::arg_list::iterator li ;
            for(li=species_list.begin();li!=species_list.end();++li) {
              Loci::option_values::value_list_type species_arg ;
              li->get_value(species_arg) ;
              if(li->type_of() != Loci::NAME_ASSIGN || species_arg.size() != 1
                || species_arg.front().type_of() != Loci::REAL) {
                cerr << "error in metals list assignment" << endl ;
                Loci::Abort() ;
              }
              else {
                string sn ; double sv ;
                li->get_value(sn) ;
                species_arg.front().get_value(sv) ;
                int sid = eos.speciesIndex(sn) ;
                if(sid < 0) {
                  cerr << "ERROR: Invalid species" << sn << " in metals list" << endl;
                  Loci::Abort() ;
                }
                metals[sid] = sv ;
              }
            }
          }    
//Read oxides list----------
          vector<double> oxides(num_species,0.0) ;
          for(int s=0;s<num_species;++s) oxides[s] = 0.0 ;
          if(functionOptions.optionExists("oxides")) {
            Loci::options_list::arg_list species_list ;
            functionOptions.getOption("oxides",species_list) ;
            Loci::options_list::arg_list::iterator li ;
            for(li=species_list.begin();li!=species_list.end();++li) {
              Loci::option_values::value_list_type species_arg ;
              li->get_value(species_arg) ;
              if(li->type_of() != Loci::NAME_ASSIGN || species_arg.size() != 1
                || species_arg.front().type_of() != Loci::REAL) {
                cerr << "error in oxide list assignment" << endl ;
                Loci::Abort() ;
              }
              else {
                string sn ; double sv ;
                li->get_value(sn) ;
                species_arg.front().get_value(sv) ;
                int sid = eos.speciesIndex(sn) ;
                if(sid < 0) {
                  cerr << "ERROR: Invalid species" << sn << " in oxides list" << endl;
                  Loci::Abort() ;
                }
                oxides[sid] = sv ;
              }
            }
          }
//Read Econd--------------------
          double Econd = 0.0 ;
          if(functionOptions.optionExists("Econd")) {
            functionOptions.getOption("Econd",Econd) ;
          }
std::cout << "step 1" << endl ;
//Condrate cx(i,Econd,metals,oxides) ;
//condensation_modifier.push_back(cx) ;
          condensation_modifier.push_back(Condrate(i,Econd,metals,oxides)) ;
std::cout << "step 2" << endl ;
std::cout << "Econd=" << Econd << endl ;
std::cout << "step 3" << endl ;
for(int s=0;s<num_species;++s) {
  std::cout << "metals[" << s << "]=" << metals[s] << endl ;
//  std::cout << "cond_metals[" << s << "]=" << condensation_modifier[i].metals[s] << endl ;
}
for(int s=0;s<num_species;++s) {
  std::cout << "oxides[" << s << "]" << oxides[s] << endl ;
//  std::cout << "cond_oxides[" << s << "]=" << condensation_modifier[i].oxides[s] << endl ;
}
std::cout << "step 4" << endl ;
//-----------------------------------------------------------------------------------
//ST-----------------------------------------------------------------------------------
        } else {
          cerr << "unable to interpret rate modifier function " << fname << endl;
          Loci::Abort() ;
        }

      }

      if(ro.optionExists("exp_nu")) {
        if(!r.isforwardonly()) {
          cerr << "exp_nu only valid for forward reactions!" << endl ;
          Loci::Abort() ;
        }
        vector<double> exp_nu ;
        exp_nu = get_exp_nu(eos,ro) ;

        for(size_t j=0;j<reactants_i.size();++j) {
          reactants_i[j].nu = exp_nu[reactants_i[j].species] ;
          exp_nu[reactants_i[j].species] = 0.0 ;
        }
        for(size_t j=0;j<exp_nu.size();++j) {
          if(exp_nu[j] != 0.0) {
            cerr << "exp_nu set for species not in reaction for reaction:"
                 << endl << r << endl ;
            Loci::Abort() ;
          }
        }
      }
      for(size_t j=0;j<reactants_i.size();++j) {
        if(reactants_i[j].nu >=0 &&
           reactants_i[j].nu < 5 &&
           (fabs(double(int(reactants_i[j].nu))-reactants_i[j].nu)) < EPSILON) {
          reactants_int[i].push_back(sc_int(reactants_i[j].species,
                                            int(reactants_i[j].nu))) ;
        } else if(reactants_i[j].nu < 0) {
          reactants[i].push_back(sc(reactants_i[j].species,
                                    reactants_i[j].nu, MinMF)) ;
        } else {
          reactants[i].push_back(sc(reactants_i[j].species,
                                    reactants_i[j].nu)) ;
        }
      }
      for(size_t j=0;j<products_i.size();++j) {
        if(products_i[j].nu >=0 &&
           products_i[j].nu < 5 &&
           (fabs(double(int(products_i[j].nu))-products_i[j].nu)) < EPSILON) {
          products_int[i].push_back(sc_int(products_i[j].species,
                                           int(products_i[j].nu))) ;
        } else if(products_i[j].nu < 0) {
          products[i].push_back(sc(products_i[j].species,
                                   products_i[j].nu, MinMF)) ;
        } else {
          products[i].push_back(sc(products_i[j].species,
                                   products_i[j].nu)) ;
        }
      }
    }
    for(int i=0;i<eos.numSpecies();++i)
      m.push_back(eos.speciesMolecularMass(i)) ;

  }

  //the following function gets relevant chemical reaction rates
  void reaction::extract_rates(rates *rate_info, const EOS &eos,
                               const double T) const {
    //  std::vector<std::pair<double,double> > fres = func_res ;
    pair<double,double> func_res[1000] ;
    thermoKcModule.calculate(eos,T,func_res) ;
    ArrheniusModule.calculate(T,func_res) ;

    for(int r=0;r<num_reactions;++r) {
      const double kf = func_res[Kfi[r]].first ;
      const double kfp = func_res[Kfi[r]].second ;
      rate_info[r].kf  = kf ;
      rate_info[r].kfp = kfp ; // dKf/dT
      if(!isforwardonly[r]) {
        const double kc = func_res[Kci[r]].first ;
        const double kcp = func_res[Kci[r]].second ;
        const double kcr = (1./(kc+KCR_EPS)) ;
        rate_info[r].kcr = kcr ;
        rate_info[r].kcp =  kcp ;  // dKc/dT
      } else {
        rate_info[r].kcr = 0 ;
        rate_info[r].kcp = 0 ;
      }
    }
  }

  double Evaluate_Rox(reaction::Rox rx, const double *mixture,
                      const EOS::State &eos_state,
                      const EOS &eos) {
    const double rho = eos_state.density() ;
    const double Rs = eos.speciesR(rx.Pindex) ;
    const double T = eos_state.temperature() ;
    const double Tr = 1./T ;
    const double PO2 = mixture[rx.Pindex]*rho*Rs*T/(1.01325e5) ;
    const double KA = rx.KA1*exp(min(-rx.KA2*Tr,EXP_LIM)) ;
    const double KB = rx.KB1*exp(min(-rx.KB2*Tr,EXP_LIM)) ;
    const double KT = rx.KT1*exp(min(-rx.KT2*Tr,EXP_LIM)) ;
    const double KZ = rx.KZ1*exp(min(-rx.KZ2*Tr,EXP_LIM)) ;
    const double X = 1/(1+PO2*KT/KB) ;
    const double Rox = KA*PO2*X/(1+KZ*PO2) + (1.-X)*KB*PO2 ;
    //    cout << "Rox = " << Rox << endl ;
    return Rox/(rx.rho_s*rx.D_s) ;
  }
  double Evaluate_dRoxdT(reaction::Rox rx, const double *mixture,
                         const EOS::State &eos_state,
                         const EOS &eos) {
    // It may be best to not condsider these derivatives. This is a test.
    //    return 0 ;
    const double rho = eos_state.density() ;
    const double Rs = eos.speciesR(rx.Pindex) ;
    const double T = eos_state.temperature() ;
    const double Tr = 1./T ;
    const double Tr2 = Tr*Tr ;
    const double PO2 = mixture[rx.Pindex]*rho*Rs*T/(1.01325e5) ;
    const double KA = rx.KA1*exp(min(-rx.KA2*Tr,EXP_LIM)) ;
    const double KAP = KA*rx.KA2*Tr2 ;
    const double KB = rx.KB1*exp(min(-rx.KB2*Tr,EXP_LIM)) ;
    const double KBP = KB*rx.KB2*Tr2 ;
    const double KT = rx.KT1*exp(min(-rx.KT2*Tr,EXP_LIM)) ;
    const double KTP = KT*rx.KT2*Tr2 ;
    const double KZ = rx.KZ1*exp(min(-rx.KZ2*Tr,EXP_LIM)) ;
    const double KZP = KZ*rx.KZ2*Tr2 ;
    const double X = 1/(1+PO2*KT/KB) ;
    const double XP = (PO2*KBP*KT/(KB*KB) - PO2*KTP/KB)/(X*X);
    //    const double Rox = (KA*PO2*X/(1+KZ*PO2) + (1.-X)*KB*PO2);
    // neglecting dPO2dT
    const double dRoxdT = PO2*(KAP*X+KA*XP)/(PO2*KZ+1) -PO2*PO2*KA*KZP*X/((PO2*KZ+1)*(PO2*KZ+1)) + PO2*((1-X)*KBP - KB*XP) ;
    return dRoxdT/(rx.rho_s*rx.D_s) ;
  }


  //the following function computes the chemical source term
  void reaction::compute_w(double *w, const rates *rate_info,
                           const double *mixture,
                           const EOS::State &eos_state,
                           const EOS &eos) const
  {
    const double rho = eos_state.density() ;
    const double P = eos_state.pressure() ;
    scratch_array<double> mr(num_species) ;
    double sumM = 0.0 ;
    for(int s=0;s<num_species;++s) {
      mr[s] = rho*mixture[s]/m[s] ;   //rho_i / M_i
      sumM+=mr[s] ;
    }
    vector<double> rate_modifier(num_reactions) ;

    // Compute rate modifier for mbody reactions
    for(int r=0;r<num_reactions;++r) {
      rate_modifier[r] = 1.0 ;
      if(isMB[r]) {
        rate_modifier[r] = sumM ;
        if(hasMBwts[r]) {
          double sumMm = 0.0 ;
          for(int j=0;j<num_species;++j)
            sumMm += MBwts[r][j]*mr[j] ;
          rate_modifier[r] = sumMm ;
        }
      }
    }

    // Add other rate modifications
    for(size_t i=0;i<Rox_modifier.size();++i) {
      int r = Rox_modifier[i].reaction_index ;
      rate_modifier[r] *= Evaluate_Rox(Rox_modifier[i],mixture,eos_state,eos) ;
    }
    for(size_t i=0;i<pressure_modifier.size();++i) {
      int r = pressure_modifier[i].reaction_index ;
      const double exponent = pressure_modifier[i].power ;
      const double Pref = pressure_modifier[i].Pref ;
      rate_modifier[r] *= pow(P/Pref,exponent) ;
    }
//ST------------------------------------------------------------begin
//Condensation modifier
    const double T = eos_state.temperature() ;
    for(size_t i=0;i<condensation_modifier.size();++i) {
      int r = condensation_modifier[i].reaction_index ;
      double T3=T*T*T ;
      double sum_P_suboxides=0.0 ; double sum_P_metals=0.0 ;
      for(int s=0;s<num_species;++s) 
        sum_P_suboxides+=condensation_modifier[i].oxides[s]*mixture[s]/m[s] ;
      for(int s=0;s<num_species;++s) 
        sum_P_metals+=condensation_modifier[i].metals[s]*mixture[s]/m[s] ;
      double S_metal=max((sum_P_suboxides/(sum_P_metals+1.0e-20)+1.0),1.0) ;
      double factor=log(S_metal)*log(S_metal) ;
      double Econd=condensation_modifier[i].Econd ;
      rate_modifier[r] *= exp(Econd/(T3*factor+1.0e-10)) ; 
    }
//ST------------------------------------------------------------end
    for(int s=0;s<num_species;++s) {
      w[s]  = 0.0 ;
      for(size_t i=0;i<pmr[s].size();++i) {
        const int r = pmr[s][i].reaction ;
        double p = 1.0 ;
        for(size_t j=0;j<reactants_int[r].size();++j) {
          const int m = reactants_int[r][j].species ;
          const double pp = mr[m] ;
          const int nu = reactants_int[r][j].nu ;
          for(int k=0;k<nu;++k)
            p *= pp ;
        }
        for(size_t j=0;j<reactants[r].size();++j) {
          const int m = reactants[r][j].species ;
          const double pp = mr[m] ;
          double nu = reactants[r][j].nu ;
          double bot = reactants[r][j].bot ;
          p *= pow(max(pp+bot,0.0),nu) ;
        }
        double wt = p ;
        if(!isforwardonly[r]) {
          p = 1.0 ;
          for(size_t j=0;j<products_int[r].size();++j) {
            const int m = products_int[r][j].species ;
            const double pp = mr[m] ;
            const int nu = products_int[r][j].nu ;
            for(int k=0;k<nu;++k)
              p *= pp ;
          }
          for(size_t j=0;j<products[r].size();++j) {
            const int m = products[r][j].species ;
            const double pp = mr[m] ;
            double nu = products[r][j].nu ;
            double bot = products[r][j].bot ;
            p *= pow(max(pp+bot,0.0),nu) ;
          }
          wt -= rate_info[r].kcr*p ;
        }
        w[s] += (pmr[s][i].nu_diff)*rate_info[r].kf*wt*rate_modifier[r] ;
      }
      w[s] *= m[s] ;
//std::cout << "w=" << w[s] << endl ;   //ST--
    }
  }

  // Computes the implicit and explicit parts of the chemical source term for
  // each species. The implicit part will later be brought into the main
  // coefficient of the species equation. The explicit part remains on the RHS
  // of the species equation. IMPORTANT: This function currently assumes that a
  // species does not exist on both sides of a reaction so we can simplify the
  // logic in determining which part goes to the LHS of the equation. Thus, only
  // terms involving reactants (as opposed to products) are candidates.
  void reaction::compute_w_implicit_explicit(double *a,double *w,const rates
  *rate_info,const double *mixture,const EOS::State &eos_state,
  const EOS &eos) const {
    const double rho = eos_state.density() ;
    const double P = eos_state.pressure() ;
    scratch_array<double> mr(num_species) ;
    double sumM = 0.0 ;
    for(int s=0;s<num_species;++s) {
      mr[s] = rho*mixture[s]/m[s] ;   //rho_i / M_i
      sumM+=mr[s] ;
    }
    vector<double> rate_modifier(num_reactions) ;

    // Compute rate modifier for mbody reactions
    for(int r=0;r<num_reactions;++r) {
      rate_modifier[r] = 1.0 ;
      if(isMB[r]) {
        rate_modifier[r] = sumM ;
        if(hasMBwts[r]) {
          double sumMm = 0.0 ;
          for(int j=0;j<num_species;++j)
            sumMm += MBwts[r][j]*mr[j] ;
          rate_modifier[r] = sumMm ;
        }
      }
    }

    // Add other rate modifications
    for(size_t i=0;i<Rox_modifier.size();++i) {
      int r = Rox_modifier[i].reaction_index ;
      rate_modifier[r] *= Evaluate_Rox(Rox_modifier[i],mixture,eos_state,eos) ;
    }
    for(size_t i=0;i<pressure_modifier.size();++i) {
      int r = pressure_modifier[i].reaction_index ;
      const double exponent = pressure_modifier[i].power ;
      const double Pref = pressure_modifier[i].Pref ;
      rate_modifier[r] *= pow(P/Pref,exponent) ;
    }
//ST------------------------------------------------------------begin
//Condensation modifier
    const double T = eos_state.temperature() ;
    for(size_t i=0;i<condensation_modifier.size();++i) {
      int r = condensation_modifier[i].reaction_index ;
      double T3=T*T*T ;
      double sum_P_suboxides=0.0 ; double sum_P_metals=0.0 ;
      for(int s=0;s<num_species;++s) 
        sum_P_suboxides+=condensation_modifier[i].oxides[s]*mixture[s]/m[s] ;
      for(int s=0;s<num_species;++s) 
        sum_P_metals+=condensation_modifier[i].metals[s]*mixture[s]/m[s] ;
      double S_metal=max((sum_P_suboxides/(sum_P_metals+1.0e-20)+1.0),1.0) ;
      double factor=log(S_metal)*log(S_metal) ;
      double Econd=condensation_modifier[i].Econd ;
      rate_modifier[r] *= exp(Econd/(T3*factor+1.0e-10)) ; 
    }
//ST------------------------------------------------------------end
    for(int s=0;s<num_species;++s) {
      a[s]=w[s]=0.0 ;
      for(size_t i=0;i<pmr[s].size();++i) {
        const int r = pmr[s][i].reaction ;
        double p=1.0 ;
        for(size_t j=0;j<reactants_int[r].size();++j) {
          const int m = reactants_int[r][j].species ;
          const double pp = mr[m] ;
          const int nu = reactants_int[r][j].nu ;
          for(int k=0;k<nu;++k)
            p *= pp ;
        }
        for(size_t j=0;j<reactants[r].size();++j) {
          const int m = reactants[r][j].species ;
          const double pp = mr[m] ;
          double nu = reactants[r][j].nu ;
          double bot = reactants[r][j].bot ;
          p *= pow(max(pp+bot,0.0),nu) ;
        }

        // If species s is a reactant and there is a non-zero mass fraction
        // of s, then the forward reaction term can be made implicit.
        double wt=0.0 ;
        if(pmr[s][i].nu_diff<0 && mixture[s]>1.0e-30){
          a[s]-=(pmr[s][i].nu_diff)*rate_info[r].kf*p*rate_modifier[r]/
            mixture[s] ;
        }else{
          wt=p ;
        }

        if(!isforwardonly[r]) {
          p = 1.0 ;
          for(size_t j=0;j<products_int[r].size();++j) {
            const int m = products_int[r][j].species ;
            const double pp = mr[m] ;
            const int nu = products_int[r][j].nu ;
            for(int k=0;k<nu;++k)
              p *= pp ;
          }
          for(size_t j=0;j<products[r].size();++j) {
            const int m = products[r][j].species ;
            const double pp = mr[m] ;
            double nu = products[r][j].nu ;
            double bot = products[r][j].bot ;
            p *= pow(max(pp+bot,0.0),nu) ;
          }
          wt -= rate_info[r].kcr*p ;
        }
        w[s] += (pmr[s][i].nu_diff)*rate_info[r].kf*wt*rate_modifier[r] ;
      }
      a[s] *= m[s] ; w[s] *= m[s] ;
    }
  }

  void reaction::compute_dwdrs(double *dwdr, int s, const rates *rate_info,
                               double dTdp, const double *dTdri,
                               const double *mixture,
                               const EOS::State &eos_state,
                               const EOS &eos) const {
    vector<double> rate_modifier(num_reactions) ;
    for(int r=0;r<num_reactions;++r)
      rate_modifier[r] = 1 ;

    for(size_t i=0;i<Rox_modifier.size();++i) {
      int r = Rox_modifier[i].reaction_index ;
      rate_modifier[r] *= Evaluate_Rox(Rox_modifier[i],mixture,eos_state,eos) ;
    }
    const double P = eos_state.pressure() ;
    for(size_t i=0;i<pressure_modifier.size();++i) {
      int r = pressure_modifier[i].reaction_index ;
      const double exponent = pressure_modifier[i].power ;
      const double Pref = pressure_modifier[i].Pref ;
      rate_modifier[r] *= pow(P/Pref,exponent) ;
    }

    scratch_array<double> X(num_species) ;

    double sum_Xm = 0 ;
    for(int k=0;k<num_species;++k) {
      X[k] = eos_state.density()*mixture[k]/m[k] ;
      sum_Xm+= X[k] ;
    }

    const int ns = num_species ;
    for(int k=0;k<ns;++k) { // Loop over species k
      dwdr[k] = 0 ;
      for(size_t i=0;i<pmr[s].size();++i) { // Loop over reactions
        const int r = pmr[s][i].reaction ;
        double rp = 1.0 ;
        double rpt = 1.0 ;
        double nu_r = 0 ;

        // Loop over reactants for each equations.  First those with
        // integer coefficients, then double ones, compute exponents
        for(size_t j=0;j<reactants_int[r].size();++j) {
          const int m = reactants_int[r][j].species ;
          const int nu = reactants_int[r][j].nu ;
          if(m==k) nu_r = nu ;
          const int nut = (m==k)?nu-1:nu ;
          const double Xm = X[m] ;
          double pw = 1.0;
          for(int kk=0;kk<nut;++kk)
            pw *= Xm ;
          rpt *= pw ;
          if(nut!=nu)
            pw*=Xm ;
          rp *= pw ;
        }
        for(size_t j=0;j<reactants[r].size();++j) {
          const int m = reactants[r][j].species ;
          double nu = reactants[r][j].nu ;
          if(m==k) nu_r = nu ;
          double nut = (m==k)?nu-1:nu ;
          double bot = reactants[r][j].bot ;
          const double Xm = X[m]+bot ;
          double pw = pow(Xm,nut) ;
          rpt *= pw ;
          rp *= pw ;
          if(m==k) rp *= Xm ;
        }

        double ppt = 0 ;
        double pp  = 0 ;
        double nu_p = 0 ;
        if(!isforwardonly[r]) {
          ppt = 1 ;
          pp = 1 ;
          // Loop over products for each equations.  First those with
          // integer coefficients, then double ones, compute exponents
          for(size_t j=0;j<products_int[r].size();++j) {
            const int m = products_int[r][j].species ;
            const int nu = products_int[r][j].nu ;
            if(m==k) nu_p = nu ;
            const int nut = (m==k)?nu-1:nu ;
            const double Xm = X[m] ;
            double pw = 1.0;
            for(int kk=0;kk<nut;++kk)
              pw *= Xm ;
            ppt *= pw ;
            if(nut!=nu)
              pw*=Xm ;
            pp *= pw ;
          }
          for(size_t j=0;j<products[r].size();++j) {
            const int m = products[r][j].species ;
            double nu = products[r][j].nu ;
            if(m==k) nu_p = nu ;
            double nut = (m==k)?nu-1:nu ;
            double bot = products[r][j].bot ;
            const double Xm = X[m] ;
            double pw = pow(Xm+bot,nut) ;
            ppt *= pw ;
            pp *= pw ;
            if(m==k) pp *= Xm+bot ;
          }
        }
        double sumM = 1.0 ;
        if(isMB[r]) {
          sumM = 0.0 ;
          if(hasMBwts[r])
            for(int j=0;j<num_species;++j)
              sumM += MBwts[r][j]*X[j] ;
          else
            for(int j=0;j<num_species;++j)
              sumM += X[j] ;
        }

        const double rmk = 1./m[k] ;
        const double nu_diff = pmr[s][i].nu_diff ;
        const double Kft = rate_info[r].kf*rate_modifier[r] ;
        const double Kf = Kft*sumM ;
        const double dKfdr = rate_modifier[r]*sumM*rate_info[r].kfp*dTdri[k] ;
        // Note, this should expand further to properly account for
        // more for general Kf
        const double rKc = rate_info[r].kcr ;

        const double dKcdr = rate_info[r].kcp*dTdri[k] ; //+ Kc*kc_coef*rmk ;

        dwdr[k] += nu_diff*(rmk*Kf*(nu_r*rpt-nu_p*rKc*ppt) +
                            dKfdr*(rp-rKc*pp) +
                            (Kf*rKc)*(rKc*dKcdr)*pp ) ;
        if(isMB[r]) {// Mbody contribution
          const double wtmp = nu_diff*Kft*(rp-rKc*pp)  ;
          if(hasMBwts[r])
            dwdr[k]+= wtmp*MBwts[r][k]*rmk ;
          else
            dwdr[k]+= wtmp*rmk ;
        }
      }
      dwdr[k] *= m[s];
    }
  }

  double reaction::compute_dwdp(int s,  const rates *rate_info,
                                double dTdp, const double *mixture,
                                const EOS::State &eos_state,
                                const EOS &eos) const {
    vector<double> rate_modifier(num_reactions) ;
    for(int r=0;r<num_reactions;++r)
      rate_modifier[r] = 1 ;

    for(size_t i=0;i<Rox_modifier.size();++i) {
      int r = Rox_modifier[i].reaction_index ;
      rate_modifier[r] *= Evaluate_Rox(Rox_modifier[i],mixture,eos_state,eos) ;
    }
    const double P = eos_state.pressure() ;
    for(size_t i=0;i<pressure_modifier.size();++i) {
      int r = pressure_modifier[i].reaction_index ;
      const double exponent = pressure_modifier[i].power ;
      const double Pref = pressure_modifier[i].Pref ;
      rate_modifier[r] *= pow(P/Pref,exponent) ;
    }

    scratch_array<double> X(num_species) ;

    for(int k=0;k<num_species;++k)
      X[k] = eos_state.density()*mixture[k]/m[k] ;

    double dwdp = 0 ;

    for(size_t i=0;i<pmr[s].size();++i) { // Loop over reactions
      const int r = pmr[s][i].reaction ;

      double rp = 1.0 ;

      // Loop over reactants for each equations.  First those with
      // integer coefficients, then double ones, compute exponents
      for(size_t j=0;j<reactants_int[r].size();++j) {
        const int m = reactants_int[r][j].species ;
        const int nu = reactants_int[r][j].nu ;
        const double Xm = X[m] ;
        double pw = 1.0;
        for(int kk=0;kk<nu;++kk)
          pw *= Xm ;
        rp *= pw ;
      }
      for(size_t j=0;j<reactants[r].size();++j) {
        const int m = reactants[r][j].species ;
        double nu = reactants[r][j].nu ;
        double bot = reactants[r][j].bot ;
        const double Xm = X[m]+bot ;
        double pw = pow(Xm,nu) ;
        rp *= pw ;
      }

      double pp  = 0 ;

      if(!isforwardonly[r]) {
        pp = 1 ;
        // Loop over products for each equations.  First those with
        // integer coefficients, then double ones, compute exponents
        for(size_t j=0;j<products_int[r].size();++j) {
          const int m = products_int[r][j].species ;
          const int nu = products_int[r][j].nu ;
          const double Xm = X[m] ;
          double pw = 1.0;
          for(int kk=0;kk<nu;++kk)
            pw *= Xm ;
          pp *= pw ;
        }
        for(size_t j=0;j<products[r].size();++j) {
          const int m = products[r][j].species ;
          double nu = products[r][j].nu ;
          double bot = products[r][j].bot ;
          const double Xm = X[m]+bot ;
          double pw = pow(Xm,nu) ;
          pp *= pw ;
        }
      }
      double sumM = 1.0 ;
      if(isMB[r]) {
        sumM = 0.0 ;
        if(hasMBwts[r])
          for(int j=0;j<num_species;++j)
            sumM += MBwts[r][j]*X[j] ;
        else
          for(int j=0;j<num_species;++j)
            sumM += X[j] ;
      }
      const double nu_diff = pmr[s][i].nu_diff ;
      const double Kf = rate_info[r].kf*rate_modifier[r]*sumM ;
      const double dKfdp = rate_info[r].kfp*rate_modifier[r]*sumM*dTdp ;
      const double rKc = rate_info[r].kcr ;
      const double dKcdp = rate_info[r].kcp*dTdp ;


      dwdp += nu_diff*(dKfdp*(rp-rKc*pp)+ (Kf*rKc)*(rKc*dKcdp)*pp) ;
    }
    dwdp *= m[s] ;
    return dwdp ;
  }
}
