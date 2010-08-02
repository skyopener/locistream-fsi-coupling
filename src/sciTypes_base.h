#ifndef SCITYPES_H
#define SCITYPES_H

#include <Loci>

#include <iostream>
#include <string>
#include <sstream>
#include <ctype.h>

#include <Tools/tools.h>
#include <Tools/parse.h>
#include <Tools/unit_type.h>

namespace streamUns {

  typedef Loci::real_t real ;
  typedef float  real_fj ;

  typedef Loci::vector3d<real> vect3d ;
  typedef Loci::tensor3d<real> tens3d ;

  using Loci::vector3d ;
  using Loci::tensor3d ;
  using Loci::norm ;
  using Loci::dot ;
  using Loci::cross ;
  
  template <unsigned int n> class vec : public Loci::Array<real,n> {
  } ;

  struct symmetricTensor {
    real xx, xy, xz, yy, yz, zz ;
  } ;

  inline std::ostream & operator<<(std::ostream &s, const symmetricTensor &ss)
  {
    s << ss.xx << ' ' << ss.xy << ' ' << ss.xz << ' ' <<
      ss.yy << ' ' << ss.yz << ' ' << ss.zz << std::endl ; 
    return s ;
  }

  inline std::istream &operator>>(std::istream &s, symmetricTensor &ss)
    {
      s >> ss.xx >> ss.xy >> ss.xz >> ss.yy >> ss.yz >> ss.zz ;
      return s ;
    }

  template <class T> class tmp_array {
    int sz ;
    T data[25] ;
    T * p ;
    void alloc(int size) {
      sz = size ;
      p = data ;
      if(sz > 25)
        p = new T[sz] ;
    }
    void free() {
      if(sz > 25)
        delete[] p ;
    }
    tmp_array() { alloc(0) ; }
  public:
    tmp_array(int size) {
      alloc(size) ;
    }
    tmp_array(const tmp_array &ta) {
      alloc(ta.sz) ;
      for(int i=0;i<sz;++i)
        p[i] = ta.p[i] ;
    }
    tmp_array &operator=(const tmp_array &ta) {
      free() ;
      alloc(ta.sz) ;
      for(int i=0;i<sz;++i)
        p[i] = ta.p[i] ;
      return *this ;
    }
    ~tmp_array() { free(); }
    T & operator[](int i) { return p[i] ; }
    T & operator[](int i) const { return p[i] ; }
    operator T *() { return p ; }
    operator const T *() const { return p ; }
  } ;
      
  void get_vect3dOption(const Loci::options_list &ol,std::string vname,
                        std::string units, vect3d &vec,real Lref=1.0) ;

  struct Quaternion {
    real x,y,z,w ;
    Quaternion() {}
    Quaternion(real xi, real yi, real zi, real wi):x(xi),y(yi),z(zi),w(wi) {}
    Quaternion(vect3d axis, real angle) {
      real sinAngle; 
      angle *= 0.5; 
      axis *= 1.0/(norm(axis)+1e-30) ; 
      sinAngle = sin(angle);
      x = (axis.x * sinAngle); 
      y = (axis.y * sinAngle); 
      z = (axis.z * sinAngle); 
      w = cos(angle);
    }
    Quaternion operator*(const Quaternion &q) const {
      vect3d vector1(x,y,z), vector2(q.x,q.y,q.z); 

      const real angle = ((w * q.w) - (dot(vector1, vector2))); 
      const vect3d across = cross(vector1, vector2);
      vector1 *= q.w ;
      vector2 *= w ;
      Quaternion result; 
      result.x = (vector1.x + vector2.x + across.x); 
      result.y = (vector1.y + vector2.y + across.y); 
      result.z = (vector1.z + vector2.z + across.z); 
      result.w = angle;
      return result ;
    }
    Quaternion &Normalize() {
      const real length = 1.0 / ((x*x) + (y*y) + (z*z) + (w*w)) ;
      x*=length ;
      y*=length ;
      z*=length ;
      w*=length ;
      return *this ;
    }
    
    Quaternion Inverse() const {
      Quaternion result = *this ;
      result.x *= -1 ; 
      result.y *= -1 ; 
      result.z *= -1 ; 
      return result ;
    }
    vect3d operator*(const vect3d &v) const {
      const Quaternion Q = *this ;
      const Quaternion Qinv = Inverse() ;
      const Quaternion vQ(v.x,v.y,v.z,0) ;
      Quaternion result = vQ*Qinv ;
      result = Q*result ;
      return vect3d(result.x,result.y,result.z) ;
    }
  } ;
}

namespace Loci {

  template<unsigned int n> struct data_schema_traits<streamUns::vec<n> > {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      return getLociType(Loci::Array<streamUns::real,n>()) ;
    }
  } ;

  template<> struct data_schema_traits<streamUns::symmetricTensor> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(streamUns::symmetricTensor()) ;
      LOCI_INSERT_TYPE(ct,streamUns::symmetricTensor,xx) ;
      LOCI_INSERT_TYPE(ct,streamUns::symmetricTensor,xy) ;
      LOCI_INSERT_TYPE(ct,streamUns::symmetricTensor,xz) ;
      LOCI_INSERT_TYPE(ct,streamUns::symmetricTensor,yy) ;
      LOCI_INSERT_TYPE(ct,streamUns::symmetricTensor,yz) ;
      LOCI_INSERT_TYPE(ct,streamUns::symmetricTensor,zz) ;
      return DatatypeP(ct) ;
    }
  } ;

  template<> struct data_schema_traits<streamUns::Quaternion> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(streamUns::Quaternion()) ;
      LOCI_INSERT_TYPE(ct,streamUns::Quaternion,x) ;
      LOCI_INSERT_TYPE(ct,streamUns::Quaternion,y) ;
      LOCI_INSERT_TYPE(ct,streamUns::Quaternion,z) ;
      LOCI_INSERT_TYPE(ct,streamUns::Quaternion,w) ;
      return DatatypeP(ct) ;
    }
  } ;

}


namespace streamUns {
  inline double get_val_in_units(std::istream &s, const char *units) {
    double val ;
    s >> val ;
    while(s.peek() == ' ' || s.peek() == '\t')
      s.get() ;
    if(isalpha(s.peek())) {
      std::string units_in ;
      while(s.peek() != EOF && s.peek() != '\n' && s.peek() !='\r')
        units_in += s.get() ;
      Loci::UNIT_type units_value(Loci::UNIT_type::MKS,"general",val,units_in) ;
      if(!units_value.is_compatible(units)) {
        std::ostringstream oss ;
        oss << "Expecting units compatible with " << units
            << "!" << std::endl
            << "incompatible with input of '"<<units<<"'" << std::endl ;
        throw Loci::StringError(oss.str()) ;
      }
      val = units_value.get_value_in(units) ;
    }
    return val ;
  }

  struct TimeValue {
    double val ;
    operator double() const { return val ; }
    TimeValue &operator=(const TimeValue &v) { val = v.val ; return *this ;}
    TimeValue &operator=(double v) { val = v ; return *this; }
    TimeValue &operator=(float v) { val = v ; return *this; }
    TimeValue &operator=(int v) { val = v ; return *this; }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const TimeValue &ss)
  {
    s << ss.val << " second" << std::endl ;
    return s ;
  }
    
  inline std::istream &operator>>(std::istream &s, TimeValue &ss) {
    ss.val = get_val_in_units(s,"second") ;
    return s ;
  }        

  struct TemperatureValue {
    double val ;
    operator double() const { return val ; }
    TemperatureValue &operator=(const TemperatureValue  &v) { val = v.val ; return *this; }
    TemperatureValue &operator=(double v) { val = v ; return *this; }
    TemperatureValue &operator=(float v) { val = v ; return *this; }
    TemperatureValue &operator=(int v) { val = v ; return *this; }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const TemperatureValue &ss)
  {
    s << ss.val << " kelvin" << std::endl ;
    return s ;
  }
    
  inline std::istream &operator>>(std::istream &s, TemperatureValue &ss) {
    ss.val = get_val_in_units(s,"kelvin") ;
    return s ;
  }        

  struct PressureValue {
    double val ;
    operator double() const { return val ; }
    PressureValue &operator=(const PressureValue &v) { val = v.val ; return *this; }
    PressureValue &operator=(double v) { val = v ; return *this; }
    PressureValue &operator=(float v) { val = v ; return *this; }
    PressureValue &operator=(int v) { val = v ; return *this; }
  } ;

  inline std::ostream & operator<<(std::ostream &s, const PressureValue &ss)
  {
    s << ss.val << " pascal" << std::endl ;
    return s ;
  }
    
  inline std::istream &operator>>(std::istream &s, PressureValue &ss) {
    ss.val = get_val_in_units(s,"pascal") ;
    return s ;
  }        
  
}

namespace Loci {

  template<class T> class ValConverter {
    T & ref ;
  public:
    ValConverter(T &iref) : ref(iref) {}
    int getSize() const {
      return 1 ;
    }
    void getState(double *buf, int &size) {
      buf[0] = ref.val ;
      size = 1 ;
    }
    void setState(double *buf, int size) {
      ref.val = buf[0] ;
    }
  } ;
  
  template<> struct data_schema_traits<streamUns::TimeValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef ValConverter<streamUns::TimeValue> Converter_Type ;
  } ;
  template<> struct data_schema_traits<streamUns::TemperatureValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef ValConverter<streamUns::TemperatureValue> Converter_Type ;
  } ;
  template<> struct data_schema_traits<streamUns::PressureValue> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef double Converter_Base_Type ;
    typedef ValConverter<streamUns::PressureValue> Converter_Type ;
  } ;
}


namespace streamUns {

  class scalar_node_output : public Loci::pointwise_rule {
    private:
      string var_name ;
      string value_name ;
      const_store<float> c2n ;
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> it ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:
      scalar_node_output(const char *vname, const char *valname) ;
      virtual void compute(const Loci::sequence &seq) ;
  } ;

  class scalar_node_output_always : public Loci::pointwise_rule {
    private:
      string var_name ;
      string value_name ;
      const_store<float> c2n ;
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> it ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:
      scalar_node_output_always(const char *vname, const char *valname) ;
      virtual void compute(const Loci::sequence &seq) ;
  } ;

  class vector_node_output : public Loci::pointwise_rule {
    private:
      string var_name ;
      string value_name ;
      const_store<vector3d<float> > c2n ;
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> it ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:
      vector_node_output(const char *vname, const char *valname) ;
      virtual void compute(const Loci::sequence &seq) ;
  } ;

  class vector_node_output_always : public Loci::pointwise_rule {
    private:
      string var_name ;
      string value_name ;
      const_store<vector3d<float> > c2n ;
      const_param<int> ncyc ;
      const_param<int> ncycle ;
      const_param<int> it ;
      const_param<int> plot_modulo ;
      const_param<string> modelName ;
      param<bool> OUTPUT ;
    public:
      vector_node_output_always(const char *vname, const char *valname) ;
      virtual void compute(const Loci::sequence &seq) ;
  } ;
  
  
}
#define OUTPUT_SCALAR(X,Y) class OUT_##Y : public scalar_node_output {\
                           public:\
                           OUT_##Y() : scalar_node_output(X,#Y){}\
                           }; register_rule<OUT_##Y> register_OUT_##Y 

#define OUTPUT_SCALAR_ALWAYS(X,Y) class OUT_##Y : public scalar_node_output_always {\
                           public:\
                           OUT_##Y() : scalar_node_output_always(X,#Y){}\
                           }; register_rule<OUT_##Y> register_OUT_##Y 

#define OUTPUT_VECTOR(X,Y) class VECOUT_##Y : public vector_node_output { \
                           public:\
                           VECOUT_##Y() : vector_node_output(X,#Y){}\
                           }; register_rule<VECOUT_##Y> register_VECOUT_##Y 

#define OUTPUT_VECTOR_ALWAYS(X,Y) class OUT_##Y : public vector_node_output_always{\
                           public:\
                           OUT_##Y() : vector_node_output_always(X,#Y){}\
                           }; register_rule<OUT_##Y> register_OUT_##Y 


#endif
