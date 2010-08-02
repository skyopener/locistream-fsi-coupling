#ifndef VEC_H
#define VEC_H
#include <Tools/debug.h>
#include <Tools/stream.h>

namespace streamUns {
template<class T> class Vector {
    int sz ;
    T *vec ;
  public:
    Vector() { sz = -1; vec = 0 ; }
    ~Vector() { if(vec != 0) delete[] vec ; }
    Vector &operator=(const Vector &v) {
        if(v.sz != sz) {
            if(vec!=0) delete[] vec ;
            sz = v.sz ;
            vec = new T[sz] ;
        }
        for(int i=0;i<sz;++i)
          vec[i] = v.vec[i] ;
        return *this ;
    }
    void setVecSize(int nels) {
        if (vec != 0) delete[] vec ;
        sz = nels ; vec = new T[sz] ; }
    int getVecSize() const { return sz ; }

    T &operator[](int idx) {
#ifdef BOUNDS_CHECK
        warn(idx<0 || idx >= sz) ;
#endif
        return vec[idx]; }
    const T &operator[](int idx) const {
#ifdef BOUNDS_CHECK
        warn(idx<0 || idx >= sz) ;
#endif
        return vec[idx]; }
    operator T*() { return vec ; }
    operator const T*() const { return vec ; }
} ;

template<class T> inline ostream & operator<<(ostream &s, const Vector<T> &v) 
{
    int sz = v.getVecSize() ;
    s << " " << sz << " " ;
    for(int i=0;i<sz;++i)
      s << v[i] << " " ;
    s << endl ;
    return s ;
}

template<class T> inline istream & operator>>(istream &s, Vector<T> &v)
{
    int sz ;
    s >> sz ;
    v.setVecSize(sz) ;
    for(int i=0;i<sz;++i)
      s >> v[i] ;
    return s ;
}

template<class T> class Matrix {
    int sz ;
    T *mat ;
  public:
    template<class S> class vect {
        S *ptr ;
#ifdef BOUNDS_CHECK
        int size ;
#endif
      public:
        vect(S *p ,int sz) {
            ptr = p ;
#ifdef BOUNDS_CHECK
            size = sz ;
#endif
        }
        S &operator[](int idx) {
#ifdef BOUNDS_CHECK
            warn(idx >= size || idx < 0) ;
#endif
            return ptr[idx] ;
        }
        const S &operator[](int idx) const {
#ifdef BOUNDS_CHECK
            warn(idx >= size || idx < 0) ;
#endif
            return ptr[idx] ;
        }
        operator S*() {
            return ptr ;
        }
        operator const S *() const {
            return ptr ;
        }
    } ;
        
    Matrix() { sz = -1; mat = 0 ; }
    ~Matrix() { if(mat != 0) delete[] mat ; }
    Matrix &operator=(const Matrix &v) {
        if(v.sz != sz) {
            if(mat!=0) delete[] mat ;
            sz = v.sz ;
            mat = new T[sz*sz] ;
        }
        for(int i=0;i<sz*sz;++i)
          mat[i] = v.mat[i] ;
    }
    void setMatSize(int nels) {
        if (mat != 0) delete[] mat ;
        sz = nels ; mat = new T[sz*sz] ; }
    int getMatSize() const { return sz ; }

    vect<T> operator[](int idx) {
#ifdef BOUNDS_CHECK
        warn(idx<0 || idx >= sz) ;
#endif
        return vect<T>(mat+idx*sz,sz) ; }
    const vect<T> operator[](int idx) const {
#ifdef BOUNDS_CHECK
        warn(idx<0 || idx >= sz) ;
#endif
        return vect<T>(mat+idx*sz,sz) ;
    }
    operator T*() { return mat ; }
    operator const T*() const { return mat ; }
} ;
        
template<class T> inline ostream & operator<<(ostream &s, const Matrix<T> &m) 
{
    int sz = m.getMatSize() ;
    s << " " << sz << endl ;
    for(int j=0;j<sz;++j) {
        for(int i=0;i<sz;++i)
          s << m[j][i] << " " ;
        s << endl ;
    }
    s << endl ;
    return s ;
}

template<class T> inline istream & operator>>(istream &s, Matrix<T> &m)
{
    int sz ;
    s >> sz ;
    m.setMatSize(sz) ;
    for(int j=0;j<sz;++j)
      for(int i=0;i<sz;++i)
        s >> m[j][i] ;
    return s ;
}

}

namespace Loci {

  template<class T> class VectorSchemaConverter {
    // For the schema converter, we always store a reference to the object
    // we are converting schmata for.
    streamUns::Vector<T> &eref ;
  public:
    explicit VectorSchemaConverter(streamUns::Vector<T> &iset): eref(iset) {}
    int getSize() const {
      return eref.getVecSize() ;
    }
    void getState(T *buf, int &size) {
      size = getSize() ;
      for(int i=0;i<size;++i)
        buf[i] = eref[i] ;
    }
    void setState(T *buf, int size) {
      eref.setVecSize(size) ;
      for(int i=0;i<size;++i)
        eref[i] = buf[i] ;
    }
  } ;

  template<class T> struct data_schema_traits<streamUns::Vector<T> > {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef T Converter_Base_Type ;
    typedef VectorSchemaConverter<T> Converter_Type ;
  } ;
  
}

#endif
