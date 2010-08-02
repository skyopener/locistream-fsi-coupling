#ifndef PERIODIC_TABLE_H
#define PERIODIC_TABLE_H

#include <map>
#include <string>

namespace fluidPhysics {
  class periodic_table {
    typedef std::basic_string<char> string ;
    typedef std::map<string,int> element_map ;
    element_map elements ;
  public:
    periodic_table() ;
    ~periodic_table() {}
    bool valid_element(string el) ;
    int find_table_entry(string el) ;
    double molecular_mass(int entry) ;
    double molecular_mass(string entry) {
      return molecular_mass(find_table_entry(entry)) ; }
  } ;

}
#endif
