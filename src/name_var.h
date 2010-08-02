#ifndef NAME_VAR_H
#define NAME_VAR_H

#include <Tools/parse.h>

namespace streamUns {
  class name_var {
  public:
    std::string name ;
  } ;

  inline std::ostream & operator <<(std::ostream &s, const name_var &n)
  {
    s << n.name ;
    return s ;
  }

  inline std::istream &operator>>(std::istream &s, name_var &n)
  {
    Loci::parse::kill_white_space(s) ;
    n.name = Loci::parse::get_name(s) ;
    return s ;
  }

}

namespace Loci {

  class name_var_schema_converter {
    streamUns::name_var &ref ;
  public:
    explicit name_var_schema_converter(streamUns::name_var &iref): ref(iref) {}
    int getSize() const {
      return ref.name.size() ;
    }
    void getState(char *buf, int &size) {
      size = getSize() ;
      for(int i=0;i<size;++i)
        buf[i] = ref.name[i] ;
    }
    void setState(char *buf, int size) {
      ref.name = "" ;
      for(int i=0;i<size;++i)
        ref.name += buf[i] ;
    }
  } ;

  template<> struct data_schema_traits<streamUns::name_var> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef char Converter_Base_Type ;
    typedef name_var_schema_converter Converter_Type ;
  } ;
}

#endif
