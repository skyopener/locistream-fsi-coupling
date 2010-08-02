#ifndef CHEMISTRY_DB_H
#define CHEMISTRY_DB_H

#include <iostream>
#include <string>
#include <map>
#include <vector>

#include <Loci.h>

namespace fluidPhysics {
  struct species_component {
    int element ;
    int num_element ;
  } ;

  std::vector<species_component> interpret_species(std::string s) ;

    

  class species_options : public options_list {
  public:
    species_options() :
      options_list("m:n:gamma:theta_v:mf:href:sref:Tref:Pref:cp:eos") {}
  } ;

  class species_db {
  public:
    bool speciesExists(const std::string &name) const  {
      return species_map.find(name) != species_map.end() ; }
  private:
    std::map<std::string,int> species_map ;
    std::vector<species_options> species_table ;
    std::vector<std::string> species_names ;
    int species_count ;
    int element_count ;
    std::vector< std::vector<int> > species_elements ;
    
    int speciesIndex(const std::string &name) {
      if(speciesExists(name)) 
        return species_map[name] ;
      else {
        species_map[name] = species_count++ ;
        species_options so ;
        species_table.push_back(so) ;
        species_names.push_back(name) ;
        return species_count-1 ;
      }
    }
  public:
    species_db() {species_count = 0 ; }

    std::istream &Input(std::istream &s) ;
    std::istream &Input(std::istream &s, const species_db &defaults) ;
    std::ostream &Print(std::ostream &s) const ;

    int numSpecies() const { return species_count ; }
    const std::string &getSpeciesName(int idx) const
    { return species_names[idx] ; }
    int speciesIndex(const std::string &name) const {
      std::map<std::string,int>::const_iterator smi ;
      if((smi =species_map.find(name)) != species_map.end())
        return smi->second ;
      else 
        return -1 ;
    }
    const species_options &getSpeciesOption (int idx) const
    { return species_table[idx] ; }
    int numElements() const { return element_count ; }
    const std::vector<int> &getSpeciesElements(int idx) const
    { return species_elements[idx] ; }
    
  } ;

  inline std::ostream &operator<<(std::ostream &s, const species_db &sdb)
  { return sdb.Print(s) ; }

  inline std::istream &operator>>(std::istream &s, species_db &sdb)
  { return sdb.Input(s) ; }

  class reaction_options : public options_list {
  public:
    reaction_options() :
      options_list("Kf:Kc:MB:exp_nu:rate_modifier:MinMF") {}
  } ;

  class reaction_expression {
    typedef std::map<std::string,double> reaction_map ;
    reaction_map lhs,rhs ;
    std::string reaction ;
    bool forwardonly ;
    bool MbodyReaction ;
  public:
    struct stoichiometric_coefficient {
      stoichiometric_coefficient() {nu=0;}
      stoichiometric_coefficient(const std::string &n, double val)
      { speciesName = n ; nu = val ; }
      std::string speciesName ;
      double nu ;
    } ;
    typedef std::vector<stoichiometric_coefficient> coef_list ;
    coef_list reactants() const ;
    coef_list products() const ;
    coef_list productsMinusReactants() const ;
    bool isMbodyReaction() const ;
    bool isforwardonly() const { return forwardonly ;} ;
    reaction_expression replaceMbodySpecies(const std::string &sname) const ;
    void removeM() ;
    std::istream &Input(std::istream &s) ;
    std::ostream &Print(std::ostream &s) const ;
    const std::string &reactionName() const { return reaction; }
  } ;

  inline std::ostream &operator<<(std::ostream &s, const reaction_expression &re)
  { return re.Print(s) ; }

  inline std::istream &operator>>(std::istream &s, reaction_expression &re)
  { return re.Input(s) ; }


  class reaction_db {
  private:
    int reaction_count ;
    std::map<std::string,int> reaction_map ;
    std::vector<std::string> reaction_name ;
    std::vector<reaction_options> reaction_table ;
    std::vector<reaction_expression> reaction_exp ;

    bool reactionExists(const std::string &name) const 
    { return reaction_map.find(name) != reaction_map.end() ; }
    bool reactionExists(const reaction_expression &re) const 
    { return reactionExists(re.reactionName()) ;}
    int reactionIndex(const reaction_expression &re) ;
    int reactionIndex(const reaction_expression &re) const ;
    int createReactionIndex(const reaction_expression &re) ;
      
  public:
    reaction_db() { reaction_count = 0 ; }
    std::istream &Input(std::istream &s) ;
    std::istream &Input(std::istream &s, const reaction_db &defaults) ;
    std::ostream &Print(std::ostream &s) const ;

    int numReactions() const { return reaction_count ; }
    const std::string &getReactionName(int idx) const
    { return reaction_exp[idx].reactionName() ; }
    const reaction_options &getReactionOption(int idx) const
    { return reaction_table[idx] ; }
    const reaction_expression &getReactionStoichiometry(int idx) const
    { return reaction_exp[idx] ; }
    void expandMBodyReactions(const species_db &species) ;
    void removeM() ;
  } ;

  inline std::ostream &operator<<(std::ostream &s, const reaction_db &rdb)
  { return rdb.Print(s) ; }

  inline std::istream &operator>>(std::istream &s, reaction_db &rdb)
  { return rdb.Input(s) ; }


  class chemistry_db {
  public:
    species_db species ;
    reaction_db reactions ;

    std::istream &Input(std::istream &s) ;
    std::istream &Input(std::istream &s, const chemistry_db &defaults) ;
    std::ostream &Print(std::ostream &s) const ;

    void expandMBodyReactions() { reactions.expandMBodyReactions(species);}
    void removeM() { reactions.removeM() ;}
  } ;

  inline std::ostream &operator<<(std::ostream &s, const chemistry_db &cdb)
  { return cdb.Print(s) ; }

  inline std::istream &operator>>(std::istream &s, chemistry_db &cdb)
  { return cdb.Input(s) ; }

}

#endif
