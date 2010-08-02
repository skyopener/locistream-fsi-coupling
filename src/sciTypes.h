#ifndef SCITYPES_H
#define SCITYPES_H

// Standard library includes.
#include <string>
using std::string ;

// Loci includes.
#include <Loci>
#include <Tools/stream.h>
#include <Tools/tools.h>
using Loci::vector3d ;
using Loci::tensor3d ;
using Loci::Array ;
using Loci::options_list ;

namespace Loci {

  // Added by JW 10/29/2008
  template<class T> inline tensor3d<T> operator-(const tensor3d<T> &v1,const
  tensor3d<T> &v2) {
    return tensor3d<T>(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z) ;
  }

  // Added by JW 7/5/2007
  template<class T> inline tensor3d<T> operator+(const tensor3d<T> &v1,const
  tensor3d<T> &v2) {
    return tensor3d<T>(v1.x+v2.x,v1.y+v2.y,v1.z+v2.z) ;
  }

  template<class T> inline tensor3d<T> product(float f,const tensor3d<T> &t) {
    return tensor3d<T>(f*t.x,f*t.y,f*t.z) ;
  }
}

namespace streamUns {

  typedef double real ;
  typedef float  real_fj ;

  typedef vector3d<real> vect3d ;
  typedef tensor3d<real> tens3d ;

  // Function to perform the component-by-component product of a vect3d and a
  // vect3d.
  inline vect3d ComponentProduct(const vect3d &a,const vect3d &b) {
    return vect3d(a.x*b.x,a.y*b.y,a.z*b.z) ;
  }

  // Function to perform the dot product of a tens3d and a vect3d.
  inline vect3d dotTemp(const tens3d &a,const vect3d &b) {
    return vect3d(a.x.x*b.x+a.x.y*b.y+a.x.z*b.z,a.y.x*b.x+a.y.y*b.y+a.y.z*b.z,
      a.z.x*b.x+a.z.y*b.y+a.z.z*b.z) ;
  }

  // This function produces a tensor from two vectors. The components of the
  // first vector are used row-wise and the components of the second vector are
  // used column-wise. This is often called a dyad (Malvern, p. 36).
  inline tens3d Dyad(const vect3d &a,const vect3d &b) {
    return tens3d(vect3d(a.x*b.x,a.x*b.y,a.x*b.z),vect3d(a.y*b.x,a.y*b.y,
      a.y*b.z),vect3d(a.z*b.x,a.z*b.y,a.z*b.z)) ;
  }

  // Function to perform the scalar product of a tens3d and a tens3d. The
  // scalar product is defined in indicial notation as AijBij.
  inline real ScalarProduct(const tens3d &a,const tens3d &b) {
    return a.x.x*b.x.x+a.x.y*b.x.y+a.x.z*b.x.z+a.y.x*b.y.x+a.y.y*b.y.y+
      a.y.z*b.y.z+a.z.x*b.z.x+a.z.y*b.z.y+a.z.z*b.z.z ;
  }

  // Function to return the trace of a tens3d.
  inline real Trace(const tens3d &a) {
    return a.x.x+a.y.y+a.z.z ;
  }

  // Function to return the transpose of a tens3d.
  inline tens3d Transpose(const tens3d &a) {
    tens3d b ;
    b.x.x=a.x.x ; b.x.y=a.y.x ; b.x.z=a.z.x ;
    b.y.x=a.x.y ; b.y.y=a.y.y ; b.y.z=a.z.y ;
    b.z.x=a.x.z ; b.z.y=a.y.z ; b.z.z=a.z.z ;
    return b ;
  }

  vect3d get_vect3d(const options_list &ol,const char *vname,const char*units) ;

  void get_vect3dOption(const options_list &ol,string vname,string units,
    vect3d &vec,real Lref=1.0) ;

  template <class T, unsigned int n> class vecT : public Array<T,n> {
  } ;


  template <unsigned int n> class vec : public Array<real,n> {
  } ;

  template <unsigned int m, unsigned int n> class mat {
    vec<n> a[m] ;
  public:
    vec<n> &operator[](unsigned int indx) { return a[indx] ; }
    const vec<n> &operator[](unsigned int indx) const { return a[indx] ; }
  } ;

  template <unsigned int m, unsigned int n> ostream &operator<<(ostream &s,
  const mat<m,n> &mt) {
    for(unsigned int i=0;i<n;++i)
      s << mt[i] << endl;
    return s ;
  }

  template <unsigned int m, unsigned int n> istream &operator>>(istream &s,
  mat<m,n> &mt) {
    for(unsigned int i=0;i<n;++i)
      s >> mt[i] ;
    return s ;
  }

  struct real_vect3d {
    real r ;
    vect3d v ;
    real_vect3d() {}
    real_vect3d(real rr, vect3d vv) : r(rr),v(vv) {}
  } ;

  inline real_vect3d operator+=(real_vect3d &target, const real_vect3d &val) {
    target.r += val.r ;
    target.v += val.v ;
    return target ;
  }

  inline ostream & operator<<(ostream &s, const real_vect3d &rv)
  {
    s << rv.r << ' ' << rv.v << ' ' ;
    return s ;
  }

  inline istream &operator>>(istream &s, real_vect3d &rv)
  {
    s >> rv.r >> rv.v ;
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

}


namespace Loci {

  template<class T,unsigned int n> struct data_schema_traits<streamUns::
  vecT<T,n> >{
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      return getLociType(Array<T,n>()) ;
    }
  } ;
  
  template<unsigned int n> struct data_schema_traits<streamUns::vec<n> > {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      return getLociType(Array<streamUns::real,n>()) ;
    }
  } ;

  template<> struct data_schema_traits<streamUns::real_vect3d> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(streamUns::real_vect3d()) ;
      LOCI_INSERT_TYPE(ct,streamUns::real_vect3d,r) ;
      LOCI_INSERT_TYPE(ct,streamUns::real_vect3d,v) ;
      return DatatypeP(ct) ;
    }
  } ;

}
  
  
#endif



