#include "chemistry_db.h"
#include "periodic_table.h"

#include <Loci.h>
#include <Tools/stream.h>
#include <Tools/parse.h>
#include <entitySet.h>

#include <vector>
using namespace Loci::parse ;


using std::vector ;
using std::map ;
#include "eos.h"

namespace fluidPhysics {
  equationOfState::equationOfState() {}
  equationOfState::~equationOfState() {}
  void equationOfState::base_initialize(const species_db &species) {
    sdb = species ;
    num_species = species.numSpecies() ;
    R.resize(num_species) ;
    m.resize(num_species) ;
    mf.resize(num_species) ;
    Pref.resize(num_species) ;
    Tref.resize(num_species) ;
    n.resize(num_species) ;
    href.resize(num_species) ;
    sref.resize(num_species) ;

    // Set total species mass fraction initially zero
    double mft = 0.0 ;

    for(int i=0;i<species.numSpecies();++i) {
      string sn = species.getSpeciesName(i) ;
      // Add species name to namelist
      namelist.push_back(sn) ;  
      // Get species_table
      const options_list &ol = species.getSpeciesOption(i) ; 

      if(!ol.optionExists("m")) {
        // if species_table contain incomplete information, print a diagnostic
        cerr << "species '" << sn << "' has an incomplete specification."
             << endl ;
        cerr << "missing either molecular mass 'm;" << endl ;
        cerr << "check your model file (.mdl) for errors." << endl ;
        Loci::Abort() ;
      }

      double mr ;
      // Assign m in species_table to variable mr
      // m is the atomic mass of this species
      ol.getOption("m",mr) ;  
      m[i] = mr ;
      R[i] = Rh/mr ; // Rh is universal gas constant which is defined in const.h
    
      double Prefr = 101325.0 ;
      if(ol.optionExists("Pref"))
        ol.getOptionUnits("Pref","Pa",Prefr) ;
      Pref[i] = Prefr ;

      double Trefr = 298.0 ;
      if(ol.optionExists("Tref") )
        ol.getOptionUnits("Tref","kelvin",Trefr) ;
      Tref[i] = Trefr ;

      double mfr = 0.0 ;
      if(ol.optionExists("mf"))
        ol.getOption("mf",mfr) ;
      mf[i] = mfr ;
      mft += mfr ;  // Accumulating mass fration


      double nr = 1.5 ;
      if(ol.optionExists("n"))
        ol.getOption("n",nr) ;
      n[i] = nr ;

      double gamr = 1.+1./nr ;
      if(ol.optionExists("gamma"))
        ol.getOption("gamma",gamr) ;
      n[i] = 1./(gamr-1.) ;

      double hrefr = 0 ;
      if(ol.optionExists("href"))
        ol.getOption("href",hrefr) ;
      href[i] = hrefr/m[i] ;

      double srefr = 0 ;
      if(ol.optionExists("sref"))
        ol.getOption("sref",srefr) ;
      sref[i] = srefr/m[i] ;

    }

    num_elements = species.numElements() ;

    for(int s=0;s<num_species;++s) {
      // Get array of species element from species database
      std::vector<int> v = species.getSpeciesElements(s) ;
      species_elements.push_back(v) ; 
    }

    // Check to see that the default mixture fractions make sense.
    if(num_species == 1) {
      mf[0] = 1.0 ;
      mft = 1.0 ;
    } else if(abs(1.0-mft) > 0.001) {
      cerr << "Default mixture fractions not specified correctly," << endl
           << "They should sum to 1.0, instead they sum to " << mft << endl ;
      cerr << "mf = {" <<namelist[0] << "=" << mf[0] ;
      for(int i=1;i<num_species;++i)
        cerr << ","<<namelist[i]<<"="<< mf[i] ;
      cerr << "}" << endl ;
      abort() ;
    }
    // normalize any small differences in mf
    // (enforcing total mass fraction mft=1.0)
    for(int i=0;i<num_species;++i)
      mf[i] = mf[i]/mft ;

  }


  namespace {
    periodic_table pt ;
  }

  vector<species_component> interpret_species(string s)
  {
    vector<species_component> v ;
    istringstream iss(s) ;
    if(s[0] == '_') {
      species_component c ;
      c.element = -1 ;
      c.num_element = 1 ;
      v.push_back(c) ;
      return v ;
    }
    while(iss.peek() != char_traits<char>::eof()) {
      if(iss.peek() == '(') {
        species_component c ;
        c.element = pt.find_table_entry("e") ;
        if(get_token(iss,"(+)")) {
          c.num_element = -1 ;
          v.push_back(c) ;
        } else if(get_token(iss,"(-)")) {
          c.num_element = 1 ;
          v.push_back(c) ;
        }
        if(iss.peek() != char_traits<char>::eof()) {
          cerr << "error reading ion from database" << endl ;
          Loci::Abort() ;
        }
        return v ;
      }
                
      if((iss.peek() != 'e') && !isupper(iss.peek())) {
        cerr << "syntax error in species database" << endl ;
        cerr << "iss.peek = '" << iss.peek() << "'"<< endl ;
        Loci::Abort() ;
        vector<species_component> ve ;
        return ve ;
      }
      string element ;
      element += iss.get() ;
      if(!iss.eof() && islower(iss.peek()))
        element += iss.get() ;
      species_component c ;
      c.element = pt.find_table_entry(element) ;
      if(c.element == -1) {
        cerr << "invalid element name " << element << endl ;
        Loci::Abort() ;
        vector<species_component> ve ;
        return ve ;
      }
      c.num_element = 1 ;
      if(!iss.eof() && isdigit(iss.peek()))
        iss >> c.num_element ;
      v.push_back(c) ;
    }
    return v ;
  }

  istream &species_db::Input(istream &s)
  {
    species_db d ;
    return Input(s,d) ;
  }

  istream &species_db::Input(istream &s, const species_db &defaults)
  {
    kill_white_space(s) ;
    if(s.peek() != '{') {
      cerr << "error reading species_db, '{' is missing" << endl ;
      Loci::Abort() ;
      return s ;
    }
    s.get() ;

    for(;;) {
      kill_white_space(s) ;
      if(s.peek() == '}') {
        s.get() ;
        break ;
      }
      if(!is_name(s)) {
        cerr << "error reading species_db" << endl ;
        Loci::Abort() ;
        break ;
      }
      string name = get_name(s) ;
      if(get_token(s,"(+)")) 
        name += "(+)" ;
      else if(get_token(s,"(-)")) 
        name += "(-)" ;
        
      kill_white_space(s) ;
      vector<species_component> scv = interpret_species(name) ;
        
      if(s.peek() == ':' || s.peek() == '=') {
        bool reuse = s.peek() == ':' ;
        s.get() ;
        kill_white_space(s) ;
        species_options so ;

        if(scv.size() > 0) {
          double m = 0.0  ;
          for(vector<species_component>::iterator scvi=scv.begin();
              scvi!=scv.end();++scvi) {
            m += double(scvi->num_element)*
              pt.molecular_mass(scvi->element) ;
          }
          so.setOption("m",m) ;
//          if(scv.size() == 1 && (*scv.begin()).num_element == 1) {
//            so.setOption("n",1.5) ;
//          }
        }
        if(reuse) {
          if(speciesExists(name))
            so = species_table[speciesIndex(name)] ;
          else if(defaults.speciesExists(name)) 
            so = defaults.species_table[defaults.speciesIndex(name)] ;
        } 
        s >> so ;
        species_table[speciesIndex(name)] = so ;
      } else {
        if(!defaults.speciesExists(name)) {
          if(scv.size() > 0) {
            species_options so ;
            double m = 0.0  ;
            for(vector<species_component>::iterator scvi=scv.begin();
                scvi!=scv.end();++scvi) {
              if(scvi->element != -1)
                m += double(scvi->num_element)*
                  pt.molecular_mass(scvi->element) ;
            }
            so.setOption("m",m) ;
            if(scv.size() == 1 && (*scv.begin()).num_element == 1) {
//              so.setOption("n",1.5) ;
              species_table[speciesIndex(name)] = so ;
            }  else {
              cerr << "species " << name
                   << " has no default specification" << endl ;
              Loci::Abort() ;
            }
                        
          } else{
            cerr << "species " << name << " has no default specification"
                 << endl ;
            Loci::Abort() ;
          }
        } else
          species_table[speciesIndex(name)] =
            defaults.species_table[defaults.speciesIndex(name)] ;
      }
      kill_white_space(s) ;
      if(s.peek() != ';') {
        cerr << "error in statement termination for species "
             << name << endl ;
        cerr << "expecting a ';'" << endl ;
        Loci::Abort() ;
      } else
        s.get() ;
    }

    entitySet elements ;
    int k = -1 ;

    for(int sp=0;sp<species_count;++sp) {
      vector<species_component> scv = interpret_species(species_names[sp]) ;
      for(unsigned int e=0;e<scv.size();++e) {
        if(scv[e].element < 0)
          elements += k-- ;
        else
          elements += scv[e].element ;
      }
    }

    map<int, int> element_map ;
    element_count = 0 ;
    for(entitySet::const_iterator ei=elements.begin();ei!=elements.end();++ei)
      element_map[*ei] = element_count++ ;

    warn(element_count > species_count) ;
    vector<int> vi ;
    for(int i=0;i<element_count;++i)
      vi.push_back(0) ;
    for(int sp=0;sp<species_count;++sp)
      species_elements.push_back(vi) ;

    k = -1 ;
    for(int sp=0;sp<species_count;++sp) {
      vector<species_component> scv = interpret_species(species_names[sp]) ;
      for(unsigned int e=0;e<scv.size();++e) {
        int el ;
        if(scv[e].element < 0) {
          el = k-- ; 
          species_elements[sp][element_map[el]] = 1 ;
        } else {
          el = scv[e].element ;
          species_elements[sp][element_map[el]] = scv[e].num_element ;
        }
      }
    }

    return s ;
  }

  ostream &species_db::Print(ostream &s) const
  {
    s << "{" << endl ;
    for(int i=0;i<species_count;++i)
      s << species_names[i] << ":" << species_table[i] << ";" << endl ;
    s << "}" ;
    return s ;
  }


  reaction_expression::coef_list reaction_expression::reactants() const
  {
    reaction_map::const_iterator rmi ;
    coef_list l ;
    for(rmi=lhs.begin();rmi!=lhs.end();++rmi) {
      stoichiometric_coefficient sc(rmi->first,rmi->second) ;
      l.push_back(sc) ;
    }
    return l ;
  }

  reaction_expression::coef_list reaction_expression::products() const
  {
    reaction_map::const_iterator rmi ;
    coef_list l ;
    for(rmi=rhs.begin();rmi!=rhs.end();++rmi) {
      stoichiometric_coefficient sc(rmi->first,rmi->second) ;
      l.push_back(sc) ;
    }
    return l ;
  }

  reaction_expression::coef_list
  reaction_expression::productsMinusReactants() const
  {
    reaction_map diff ;
    reaction_map::const_iterator rmi ;
    for(rmi=rhs.begin();rmi!=rhs.end();++rmi)
      diff[rmi->first] = rmi->second ;
    for(rmi=lhs.begin();rmi!=lhs.end();++rmi) {
      if(diff.find(rmi->first) == diff.end())
        diff[rmi->first] = - rmi->second ;
      else
        diff[rmi->first] += - rmi->second ;
    }
    
    coef_list l ;
    for(rmi=diff.begin();rmi!=diff.end();++rmi) 
      if(fabs(rmi->second) > EPSILON) {
        stoichiometric_coefficient sc(rmi->first,rmi->second) ;
        l.push_back(sc) ;
      }
    return l ;
  }


  bool reaction_expression::isMbodyReaction() const
  {
    //    reaction_map::const_iterator rmi = lhs.find("M") ;
    //    return rmi != lhs.end() ;
    return MbodyReaction ;
  }

  reaction_expression
  reaction_expression::replaceMbodySpecies(const string &sname) const
  {
    int lm = 0 ;
    int rm = 0 ;
    reaction_expression r ;
    reaction_map::const_iterator rmi ;
    for(rmi=lhs.begin();rmi != lhs.end();++rmi) {
      if(rmi->first != string("M")) 
        r.lhs[rmi->first] = rmi->second ;
      else
        lm++ ;
    }
    for(rmi=rhs.begin();rmi != rhs.end();++rmi) {
      if(rmi->first != string("M"))
        r.rhs[rmi->first] = rmi->second ;
      else
        rm++ ;
    }
    if(lm!=1 || rm !=1) {
      cerr << "Malformed M-body Reaction: " << *this << endl ;
      Loci::Abort() ;
    }
    if(r.lhs.find(sname)==r.lhs.end())
      r.lhs[sname] = 1 ;
    else
      r.lhs[sname] += 1 ;
    if(r.rhs.find(sname)==r.rhs.end())
      r.rhs[sname] = 1 ;
    else
      r.rhs[sname] += 1 ;
    
    // create reaction string
    ostringstream oss ;
    for(rmi=r.lhs.begin();rmi!=r.lhs.end();++rmi) {
      if(rmi!=r.lhs.begin())
        oss << "+" ;
      if(rmi->second != 1)
        oss << rmi->second ;
      oss << rmi->first ;
    }
    oss << "<->" ;
    for(rmi=r.rhs.begin();rmi!=r.rhs.end();++rmi) {
      if(rmi!=r.rhs.begin())
        oss << "+" ;
      if(rmi->second != 1)
        oss << rmi->second ;
      oss << rmi->first ;
    }
    r.reaction = oss.str() ;
    return r ;
  }

  void reaction_expression::removeM()
  {
    int lm = 0 ;
    int rm = 0 ;
    reaction_map lhstmp, rhstmp ;
    reaction_map::const_iterator rmi ;
    for(rmi=lhs.begin();rmi != lhs.end();++rmi) {
      if(rmi->first != string("M")) 
        lhstmp[rmi->first] = rmi->second ;
      else
        lm++ ;
    }
    for(rmi=rhs.begin();rmi != rhs.end();++rmi) {
      if(rmi->first != string("M"))
        rhstmp[rmi->first] = rmi->second ;
      else
        rm++ ;
    }
    if(lm!=1 || rm !=1) {
      cerr << "Malformed M-body Reaction: " << *this << endl ;
      Loci::Abort() ;
    }
    
    // create reaction string
    ostringstream oss ;
    for(rmi=lhstmp.begin();rmi!=lhstmp.end();++rmi) {
      if(rmi!=lhstmp.begin())
        oss << "+" ;
      if(rmi->second != 1)
        oss << rmi->second ;
      oss << rmi->first ;
    }
    if(forwardonly) 
      oss << "->" ;
    else 
      oss << "<->" ;
    for(rmi=rhstmp.begin();rmi!=rhstmp.end();++rmi) {
      if(rmi!=rhstmp.begin())
        oss << "+" ;
      if(rmi->second != 1)
        oss << rmi->second ;
      oss << rmi->first ;
    }
    lhs = lhstmp ;
    rhs = rhstmp ;
    reaction = oss.str() ;
  }

  istream &reaction_expression::Input(istream &s)
  {
    reaction_map *side = &lhs ;
    for(int i=0;;) {
      kill_white_space(s) ;
      double v = 1 ;
      string name ;
      if(is_real(s)) {
        v = get_real(s) ;
      }
      kill_white_space(s) ;
      if(is_name(s)) {
        name = get_name(s) ;
        if(get_token(s,"(+)"))
          name += "(+)" ;
        else if(get_token(s,"(-)"))
          name += "(-)" ;
            
        reaction_map::iterator rmi = side->find(name) ;
        if(rmi == side->end())
          (*side)[name] = v ;
        else
          rmi->second += v ;
        kill_white_space(s) ;
        if(!get_token(s,"+")) {
          if(i==1) {
            ostringstream oss ;
            MbodyReaction = false ;
            for(rmi=lhs.begin();rmi!=lhs.end();++rmi) {
              if(rmi!=lhs.begin())
                oss << "+" ;
              if(rmi->second != 1)
                oss << rmi->second ;
              oss << rmi->first ;
              if(rmi->first =="M")
                MbodyReaction = true ;
            }
            if(forwardonly) 
              oss << "->" ;
            else
              oss << "<->" ;
            for(rmi=rhs.begin();rmi!=rhs.end();++rmi) {
              if(rmi!=rhs.begin())
                oss << "+" ;
              if(rmi->second != 1)
                oss << rmi->second ;
              oss << rmi->first ;
            }
            reaction = oss.str() ;
            return s ;
          }
          if(get_token(s,"->"))
            forwardonly = true ;
          else if(get_token(s,"<->")) 
            forwardonly = false ;
          else {
            cerr << "error reading reaction_expression" << endl ;
            Loci::Abort() ;
            return s ;
          }
          i=1 ;
          side = &rhs ;
        }
      } else {
        cerr << "error reading reaction_expression, expected species name" << endl ;
        string tmp ;

        s >> tmp ;
        cerr << "got " << tmp << " instead." << endl ;
        Loci::Abort() ;
        return s ;
      }
    }
  }

  ostream &reaction_expression::Print(ostream &s) const
  {
    s << reaction ;
    return s ;
  }

  int reaction_db::createReactionIndex(const reaction_expression &re) {
    reaction_map[re.reactionName()] = reaction_count++ ;
    reaction_name.push_back(re.reactionName()) ;
    reaction_options ro ;
    reaction_table.push_back(ro) ;
    reaction_exp.push_back(re) ;
    return reaction_count-1 ;
  }
  
  int reaction_db::reactionIndex(const reaction_expression &re)
  {
    if(reactionExists(re.reactionName()))
      return reaction_map[re.reactionName()] ;
    else {
      return createReactionIndex(re) ;
    }
  }

  int reaction_db::reactionIndex(const reaction_expression &re) const
  {
    std::map<std::string,int>::const_iterator rmi ;
    if((rmi=reaction_map.find(re.reactionName())) != reaction_map.end())
      return rmi->second ;
    else
      return -1 ;
  }

  istream &reaction_db::Input(istream &s)
  {
    reaction_db d ;
    return Input(s,d) ;
  }

  istream &reaction_db::Input(istream &s, const reaction_db &defaults)
  {
    kill_white_space(s) ;
    if(s.peek() != '{') {
      cerr << "error reading reaction_db, '{' is missing" << endl ;
      Loci::Abort() ;
      return s ;
    }
    s.get() ;
    for(;;) {
      kill_white_space(s) ;
      if(s.peek() == '}') {
        s.get() ;
        return s ;
      }
      reaction_expression re ;
      s >> re ;
      kill_white_space(s) ;
      if(s.peek() == ':' || s.peek() == '=') {
        bool reuse = s.peek() == ':' ;
        s.get() ;
        kill_white_space(s) ;
        reaction_options so ;
        if(reuse) {
          if(reactionExists(re))
            so = reaction_table[reactionIndex(re)] ;
          else if(defaults.reactionExists(re))
            so = defaults.reaction_table[defaults.reactionIndex(re)] ;
          s >> so ;
          reaction_table[reactionIndex(re)] = so ;
        } else {
          s >> so ;
          reaction_table[createReactionIndex(re)] = so ;
        }
      } else {
        if(!defaults.reactionExists(re)) {
          cerr << "reaction " << re
               << " has no default specification" << endl ;
          Loci::Abort() ;
        } else
          reaction_table[reactionIndex(re)] =
            defaults.reaction_table[defaults.reactionIndex(re)] ;
      }

      kill_white_space(s) ;
      if(s.peek() != ';') {
        cerr << "error in statement termination for reaction "
             << re << endl ;
        Loci::Abort() ;
      }
      else
        s.get() ;
    }
  }

  ostream &reaction_db::Print(ostream &s) const
  {
    s << '{' << endl ;
    for(int i=0;i<reaction_count;++i)
      s << reaction_exp[i] << '=' << reaction_table[i] << ';' << endl ;
    s << '}' ;
    return s ;
  }

  void reaction_db::expandMBodyReactions(const species_db &species)
  {
    vector<reaction_options> rtable ;
    vector<reaction_expression> rexp ;
    map<string,int> rmap ;
    int rcount = 0 ;
    for(int i=0;i<reaction_count;++i) {
      if(reaction_exp[i].isMbodyReaction()) {
        for(int s=0;s<species.numSpecies();++s) {
          string sname = species.getSpeciesName(s) ;
          reaction_expression
            r = reaction_exp[i].replaceMbodySpecies(sname) ;
          if(!reactionExists(r.reactionName())) {
            rexp.push_back(r) ;
            rtable.push_back(reaction_table[i]) ;
            rmap[r.reactionName()] = rcount ;
            rcount++ ;
          }
        }
      } else {
        rexp.push_back(reaction_exp[i]) ;
        rtable.push_back(reaction_table[i]) ;
        rmap[reaction_exp[i].reactionName()] = rcount ;
        rcount++ ;
      }
    }
    reaction_count = rcount ;
    reaction_exp = rexp ;
    reaction_table = rtable ;
    reaction_map = rmap ;
  }

  void reaction_db::removeM() {  
    for(int i=0;i<reaction_count;++i) 
      if(reaction_exp[i].isMbodyReaction()) 
        reaction_exp[i].removeM() ;
  }
      
  istream &chemistry_db::Input(istream &s)
  {
    chemistry_db d ;
    return Input(s,d) ;
  }

  istream &chemistry_db::Input(istream &s,const chemistry_db &d)
  {
    kill_white_space(s) ;
    while(s.good() && !s.eof() && s.peek() != char_traits<char>::eof()) {
      kill_white_space(s) ;
      if(get_token(s,"species")) {
        kill_white_space(s) ;
        if(s.peek() != '=') {
          cerr << "error reading chemistry_db" << endl ;
          Loci::Abort() ;
          return s ;
        }
        s.get() ;
        species.Input(s,d.species) ;
        kill_white_space(s) ;
        if(s.peek() != ';') {
          cerr << "expected ';' reading species" << endl ;
          Loci::Abort() ;
          return s ;
        } else
          s.get() ;
            
      } else if(get_token(s,"reactions")) {
        kill_white_space(s) ;
        if(s.peek() != '=') {
          cerr << "error reading chemistry_db" << endl ;
          Loci::Abort() ;
          return s ;
        }
        s.get() ;
        reactions.Input(s,d.reactions) ;
        kill_white_space(s) ;
        if(s.peek() != ';') {
          cerr << "expected ';' reading species" << endl ;
          Loci::Abort() ;
          return s ;
        } else
          s.get() ;
      } else {
        kill_white_space(s) ;
        if(s.peek() != char_traits<char>::eof()) {
          cerr << "unexpected input in chemistry_db" << endl ;
          Loci::Abort() ;
          return s ;
        }
      }
    }

    return s ;
  }

  ostream &chemistry_db::Print(ostream &s) const {
    s << "species = " << species << ";" << endl ;
    s << "reactions = " << reactions << ";" << endl ;
    return s ;
  }


}
