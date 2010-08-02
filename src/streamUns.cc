//-----------------------------------------------------------------------------
// Description: This file contains the main function for LOCI-STREAMUNS.
//
// Authors: Siddharth Thakur and Jeff Wright
//          Streamline Numerics Inc.
//          Gainesville, Florida 32609
//-----------------------------------------------------------------------------

// System includes.
#include <sys/stat.h>
#include <unistd.h>

// Standard library includes.
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using std::cerr ;
using std::cout ;
using std::endl ;
using std::fstream ;
using std::ifstream ;
using std::ios ;
using std::ostream ;
using std::string ;
using std::vector ;

// Petsc includes.
#include "petsc.h"

// Loci includes.
#include <Loci.h>
#include <Tools/fpe.h>
using Loci::storeRepP ;
using Loci::variableSet ;
using Loci::variable ;
using Loci::MPI_processes ;
using Loci::MPI_rank ;

// StreamUns includes.
#include <readGrid.h>
#include "name_var.h"
using streamUns::ReadGrid ;
using streamUns::name_var ;

//-----------------------------------------------------------------------------
// Functions to print version and release date.

namespace streamUns {
  char *revisionName = "$Name: rel-1-5-4 $" ;
  char *revisionDate = "$Date: 2010/06/24 18:17:18 $" ;

  string Version() {
    string version(revisionName) ; version.erase(0,6) ;
    version.erase(version.size()-1) ; return version ;
  }

  string Date() {
    string date(revisionDate) ; date.erase(0,6) ; date.erase(date.size()-1) ;
    return date ;
  }

  void PrettyPrintString(string i,string s,ostream &o) {
    o << i  << ": " ;
    size_t initial_space=i.size()+2 ; size_t count=initial_space ;
    string current_word ; size_t loc=0 ;
    do {
      current_word = "" ;
      while(loc<s.size() && s[loc]!=' ') current_word += s[loc++] ;
      if(loc<s.size()) loc++ ;
      if(count+current_word.size()>=79){
        o << endl ; for(size_t t=0;t<initial_space;++t) o << ' ' ;
        count = initial_space ;
      }
      o << current_word << ' ' ; count+=current_word.size()+1 ;
    } while(loc<s.size()) ;
    o << endl ;
  }

  void DescribeInputs(rule_db &rdb) {
    using namespace Loci ; fact_db local ;
    ruleSet special_rules = rdb.get_default_rules() ;

    // Process the default rules first.
    cout << "-------------------------------------------------------------"
      << endl ;
    for(ruleSet::const_iterator ri=special_rules.begin();ri!=special_rules.
    end();++ri) {

      // Create the facts in the fact_db.
      variableSet targets=ri->targets() ;
      rule_implP rp=ri->get_rule_implP() ;
      for(variableSet::const_iterator vi=targets.begin();vi!=targets.end();
      ++vi) {

        // Get the storeRep for this variable.
        storeRepP srp=rp->get_store(*vi) ;
        if(srp==0) {
          cerr << "rule " << *ri << " unable to provide type for " << *vi
            << endl ; exit(-1) ;
        }
        local.create_fact(*vi,srp) ;
      }

      // Call the compute method to set the default value for this variable.
      rp->initialize(local) ; rp->compute(sequence(EMPTY)) ;
      for(variableSet::const_iterator vi=targets.begin();vi!=targets.end();
      ++vi) {
        storeRepP srp=local.get_variable(*vi) ;
        cout << *vi << ": " ; srp->Print(cout) ;
        string comment=rp->get_comments() ;
        if(comment.size()!=0) PrettyPrintString("comment",comment,cout) ;
        cout << "-------------------------------------------------------------"
          << endl ;
      }
    }

    // Process the optional rules.
    special_rules=rdb.get_optional_rules() ;
    for(ruleSet::const_iterator ri=special_rules.begin();ri!=special_rules.
    end();++ri) {

      // Create the facts in the fact_db.
      variableSet targets=ri->targets() ; rule_implP rp=ri->get_rule_implP() ;
      for(variableSet::const_iterator vi=targets.begin();vi!=targets.end();
      ++vi) {
        cout << *vi << ": NO DEFAULT VALUE" << endl ;
        string comment=rp->get_comments() ;
        if(comment.size()!=0) PrettyPrintString("comment",comment,cout) ;
        cout << "-------------------------------------------------------------"
          << endl ;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// The main function.

int main(int argC,char *argV[]) {

  // Initialize Loci.
  Loci::Init(&argC,&argV) ;

  // Set the default query.
  string query = "solution" ;

  // Other variables for CHEM compatability.
  bool use_solid_heat=false ; vector<string> fvm_model_file ;

  // Print version data. If executable invoked with no arguments, exit.
  if(MPI_rank==0){
    cout << "StreamUns version: " << streamUns::Version() << ", Date: "
      << streamUns::Date() << endl ;
    cout << "Loci version: " << Loci::version() << endl ;
  }
  if(argC==1) Loci::Abort() ;

  // Process the input.
  bool inputDescription=false ; vector<string> modulesToLoad ;
  while(argC>=2 && argV[1][0]=='-') {
    if(argC>=2 && !strcmp(argV[1],"-q")){
      query = argV[2] ; argC-=2 ; argV+=2 ;
    }else if(argC>=2 && !strcmp(argV[1],"-v")) {
      if(MPI_rank==0){
        cout << "StreamUns version: " << streamUns::Version() << " Date: "
          << streamUns::Date() << endl ;
        cout << "Loci version: " << Loci::version() << endl ;
      }
      if(argC==2){ Loci::Finalize() ; exit(0) ; } argC-- ; argV++ ;
    }else if(argC>=2 && !strcmp(argV[1],"-inputs")){
      inputDescription=true ; argC-- ; argV++ ;
    }else if(argC>=3 && !strcmp(argV[1],"-load_module")) {
      modulesToLoad.push_back(string(argV[2])) ; argC-=2 ; argV+=2 ;
    }else if(!strcmp(argV[1],"-log_summary") ||
    !strcmp(argV[1],"-preload")) { // PETSC options.
      argC-=2 ; argV+=2 ;
    }else{
      if(MPI_rank==0) cerr << "Argument " << argV[1] << " is not understood."
        << endl ;
      argC-- ; argV++ ; Loci::Abort() ;
    }
  }

  // Create the output directory if one doesn't exist.
  struct stat statBufOutput ;
  if(stat("output",&statBufOutput)){
    mkdir("output",0755) ;
  }else{
    if(!S_ISDIR(statBufOutput.st_mode)) {
      if(MPI_rank==0) cerr << "File 'output' should be a directory!, "
        << "rename 'output' and start again." << endl ;
      Loci::Abort() ;
    }
  }

  // Don't know how this works, but allows the user to "touch stop" to stop
  // a job. There must be more than this though.
  if(!stat("stop",&statBufOutput)) ::unlink("stop") ;

  // Create the restart directory if one doesn't exist.
  struct stat statBufRestart ;
  if(stat("restart",&statBufRestart)){
    mkdir("restart",0755) ;
  }else{
    if(!S_ISDIR(statBufRestart.st_mode)) {
      if(MPI_rank==0) cerr << "File 'restart' should be a directory!, "
        << "rename 'restart' and start again." << endl ;
      Loci::Abort() ;
    }
  }

  // Set up floating point exception handlers.
  set_fpe_abort() ;

  /* Courtesy code for Ed so he doesn't have to delete the '.'.
  int argVOneSize=strlen(argV[1]) ;
cout << "argVOneSize: " << argVOneSize << endl ;
  if(argVOneSize>0 && (argV[1][argVOneSize-1]=='.'))
    argV[1][argVOneSize-1]='\0' ;*/

  // Add rules registered into the global rule list to the rule database.
  rule_db rdb ; rdb.add_rules(global_rule_list) ;

  // Load the finite volume module.
  Loci::load_module("fvm",rdb) ;

  // Load all modules specified in the vars file.
  if(argC > 1) {
    string varsfile = string(argV[1])+string(".vars") ;
    try {
      ifstream ifile(varsfile.c_str(),ios::in) ;
      if(ifile.fail()) {
        cerr<<"can't open " << varsfile << endl ;
        Loci::Abort() ;
      }
      Loci::parse::kill_white_space(ifile) ;
      while(ifile.peek() != '{' && ifile.peek() != EOF) {
        if( Loci::parse::is_name(ifile) ) {
          string name = Loci::parse::get_name(ifile) ;
          Loci::parse::kill_white_space(ifile) ;
          if(ifile.peek() != ':') {
            cerr << "expected ':' after '" << name << "' in file \""
                 << varsfile << '"' << endl ;
            Loci::Abort() ;
          }
          ifile.get() ;
          Loci::parse::kill_white_space(ifile) ;
          string argument ;
          if(Loci::parse::is_string(ifile)) {
            argument = Loci::parse::get_string(ifile) ;
          } else if(Loci::parse::is_name(ifile)) {
            argument = Loci::parse::get_name(ifile) ;
          } else {
            cerr << "unable to parse argument to option '" << name
                 << "' in file \"" << varsfile << '"' << endl ;
            Loci::Abort() ;
          }
          if(name == "loadModule") {
            modulesToLoad.push_back(argument) ;
          } else if(name == "conjugateHeat") {
            use_solid_heat = true ;
            fvm_model_file.push_back(argument) ;
          } else if(name == "query") {
            query = argument ;
          } else {
            cerr << "unable to interpret preamble directive '"
                 << name << "' found in file \"" << varsfile << '"' << endl ;
            Loci::Abort() ;
          }
        } else {
          cerr << "problem parsing preamble of '" << varsfile << "'" << endl ;
          string s ;
          ifile >> s ;
          cerr << "problem near token '" << s << "'" << endl ;
          Loci::Abort() ;
        }
        Loci::parse::kill_white_space(ifile) ;
      }
                                                                                
    }catch(const Loci::BasicException &err) {
      err.Print(cerr) ;
      cerr << "aborted reading \"" << varsfile << "\"" << endl ;
      Loci::Abort() ;
    }
  }

  // Load all modules that are specified on the command line.
  fact_db facts ;
  for(size_t i=0;i<modulesToLoad.size();++i){
    Loci::exprP ep=Loci::expression::create(modulesToLoad[i]) ;
    string mod_name ; string mod_namespace ; string load_file ;
    std::set<std::string> str_set ;
    if(ep->op == Loci::OP_FUNC){
      mod_name = ep->name ;
      if(ep->expr_list.front()->op == Loci::OP_NAME)
        load_file = ep->expr_list.front()->name ;
      else
        cerr << "unable to interpret argument in " << modulesToLoad[i] << endl ;
    }else if(ep->op==Loci::OP_NAME){
      mod_name = ep->name ;
    }else if(ep->op==Loci::OP_SCOPE){
      Loci::exprList::const_iterator li ;
      bool first_time=true ; Loci::exprP last_ep=0 ;
      for(li=ep->expr_list.begin();li!=ep->expr_list.end();++li){
        if(!first_time){
          if(last_ep->op==Loci::OP_NAME)
            if(mod_namespace=="") mod_namespace=last_ep->name ;
            else mod_namespace += string("_") + last_ep->name ;
          else cerr << "-load_module namespace not recognized" << endl ;
        }
        last_ep=*li ; first_time=false ;
      }
      ep=last_ep ;
      if(ep->op==Loci::OP_FUNC){
        mod_name=ep->name ;
        if(ep->expr_list.front()->op==Loci::OP_NAME)
          load_file=ep->expr_list.front()->name ;
        else
          cerr << "unable to interpret argument in " << modulesToLoad[i]
               << endl ;
      }else if(ep->op==Loci::OP_NAME){
        mod_name=ep->name ;
      }else{
        cerr << "unable to interpret " << modulesToLoad[i] << endl ;
      }
    }else{
      cerr << "unable to interpret " << modulesToLoad[i] << endl ;
    }
    Loci::load_module(mod_name,mod_namespace,load_file.c_str(),facts,rdb,
      str_set) ;
  }

  // Describe the inputs.
  if(inputDescription){
    streamUns::DescribeInputs(rdb) ; Loci::Finalize() ; exit(0) ;
  }

  // Read grid, connectivity information and user supplied information into
  // the fact database.
  ReadGrid(facts,rdb,argV[1]) ;

  // Scan for the debug flag.
  param<bool> debug ; *debug=false ;
  if(argC>=3){
    if(!strcmp(argV[2],"debug")){
      *debug=true ;
      if(MPI_rank==0) cout << "Debugging output turned on." << endl ; argV++ ;
    }else{
      *debug=false ;
      if(MPI_rank==0) cout << "Debugging output turned off." << endl ;
    }
  }
  facts.create_fact("debug",debug) ;

  // Add restart data to fact database if user specifies a restart file.
  if((!(*debug) && argC==3) || argC==4){
    if(argV[2][0]<'0' || argV[2][0]>'9'){
      if(MPI_rank==0) cout << "ERROR: Restart argument should be the desired "
        << "restart file number." << endl ;
      Loci::Abort() ;
    }
    param<string> restartNum ;
    *restartNum=argV[2] ;
    facts.create_fact("restartNum",restartNum) ;
    param<bool> restart ; *restart=true ;
    facts.create_fact("restart",restart) ;
  }else{
    param<bool> noRestart ; *noRestart=true ;
    facts.create_fact("noRestart",noRestart) ;
  }

  // Create execution schedule that derives the query variable from the database
  // of facts using the rule database rdb.
  if(!Loci::makeQuery(rdb,facts,query)) {
    if(MPI_rank==0) cerr << "Query for " << query << " failed!" << endl ;
    Loci::Abort() ;
  }

  // Finalize Loci. The exit is somehow required to prevent getting
  // an annoying message about "invalid pointer".
  Loci::Finalize() ; exit(0) ;
}

//void LoadEOSRules() { extern void loadEOSRules() ; loadEOSRules() ; }















