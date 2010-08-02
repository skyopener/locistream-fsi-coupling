#ifndef FLUID_PHYSICS_MATRIX_H
#define FLUID_PHYSICS_MATRIX_H
#include <Tools/stream.h>
#include <valarray> 
#include <Loci.h>
#include <Tools/parse.h>


namespace fluidPhysics {
  using std::valarray ;
  using std::slice ;



  template<class T>class Slice_iter {
    valarray<T>* v ;
    slice s ;
    size_t curr ;

    T& ref(size_t i) const { return (*v)[s.start()+i*s.stride()] ; }
  public:
    Slice_iter(valarray<T>* vv, slice ss) : v(vv), s(ss) {}
    Slice_iter end() {
      Slice_iter t = *this ;
      t.curr = s.start() + s.size()*s.stride() ;
      return t ;
    }

    Slice_iter& operator++() { curr++ ; return *this ; }
    Slice_iter operator++(int) { Slice_iter t = *this ; curr++ ; return t ; }

    T& operator[] (size_t i) { return ref(curr=i) ; } 
    T& operator() (size_t i) { return ref(curr=i) ; }
    T& operator*() { return ref(curr) ; }

  } ;

  template<class T> bool operator==(const Slice_iter<T>& p, const Slice_iter<T>& q) {
    return p.curr==q.curr && p.s.stride()==q.s.stride() && p.s.start()==q.s.start() ;
  }

  template<class T> bool operator!=(const Slice_iter<T>& p, const Slice_iter<T>& q) {
    return !(p==q) ;
  }

  template<class T> bool operator<(const Slice_iter<T>& p, const Slice_iter<T>& q) {
    return p.curr<q.curr && p.s.stride()==q.s.stride && p.s.start()==q.s.start() ;
  }

  template<class T>class Cslice_iter {
    valarray<T>* v ;
    slice s ;
    size_t curr ;

    T& ref(size_t i) const { return (*v)[s.start()+i*s.stride()] ; }
  public:
    Cslice_iter(valarray<T>* vv, slice ss) : v(vv), s(ss) {}
    Cslice_iter end() {
      Cslice_iter t = *this ;
      t.curr = s.start() + s.size()*s.stride() ;
      return t ;
    }

    Cslice_iter& operator++() { curr++ ; return *this ; }
    Cslice_iter operator++(int)
    {
      Slice_iter<T> t = *this ;
      curr++ ; return t ;
    }

    const T& operator[] (size_t i) { return ref(curr=i) ; }
    const T& operator() (size_t i) { return ref(curr=i) ; }
    const T& operator*() { return ref(curr) ; }

  } ;

  template<class T> bool operator==(const Cslice_iter<T>& p, const Cslice_iter<T>& q) {
    return p.curr==q.curr && p.s.stride()==q.s.stride() && p.s.start()==q.s.start() ;
  }

  template<class T> bool operator!=(const Cslice_iter<T>& p, const Cslice_iter<T>& q) {
    return !(p==q) ;
  }

  template<class T> bool operator<(const Cslice_iter<T>& p, const Cslice_iter<T>& q) {
    return p.curr<q.curr && p.s.stride()==q.s.stride && p.s.start()==q.s.start() ;}


  class Matrix2 {
    valarray<double> v ;
    size_t d1, d2 ;
  public:
    Matrix2() { d1 = 0 ; d2 = 0 ; } 
    Matrix2(size_t x, size_t y) ;
    ~Matrix2() {} ;
    size_t size() const { return d1*d2 ; } 
    size_t dim1() const { return d1 ; }
    size_t dim2() const { return d2 ; }
    std::ostream& Print(std::ostream& s) const ;
    std::istream& Input(std::istream& s) ;

    Slice_iter<double> row(size_t i) 
    { return Slice_iter<double> (&v,slice(i,d2,d1)) ; }
    Cslice_iter<double> row(size_t i) const 
    { return Cslice_iter<double> (const_cast<valarray<double> *>(&v),slice(i,d2,d1));}
    Slice_iter<double> column(size_t i) 
    {return Slice_iter<double> (&v,slice(i*d1,d1,1)) ;}
    Cslice_iter<double> column(size_t i) const 
    { return Cslice_iter<double> (const_cast<valarray<double> *>(&v),slice(i*d1,d1,1)) ;}

    double& operator() (size_t x, size_t y) ;
    double operator() (size_t x, size_t y) const ;

    Slice_iter<double> operator[] (size_t i) { return column(i) ; }
    Cslice_iter<double> operator[] (size_t i) const { return column(i) ; } 

    Matrix2& operator*=(double) ;
    valarray<double>& array() { return v ; } 

    void resize(size_t x, size_t y) ; 

    Matrix2 &operator=(const Matrix2 &mt) {
      if(this!= &mt) 
        resize(mt.dim1(),mt.dim2()) ;

      for(size_t j=0;j<d2;j++)
        for(size_t i=0;i<d1;i++)
          (*this)(i,j) = mt(i,j) ;

      return *this ;
    }
  
    void ludcmp(int *indx, double &d) ;
    void lubksb(int *indx, double *b) ;
    void inverse(Matrix2 &y) ;
    void QR() ; //QR factorization using modified Gramm-Schmidt method.

  } ;

  inline std::ostream &operator<<(std::ostream& s, const Matrix2& m) {
    return m.Print(s) ;
  }

  inline std::istream &operator>>(std::istream& s, Matrix2& m) {
    return m.Input(s) ;
  }


  class Matrix3 {
    valarray<double> v ;
    size_t d1, d2, d3 ;
    size_t d12 ;
  public:
    Matrix3() { d1 = 0 ; d2 = 0 ; d3 = 0 ; d12 = d1*d2 ;  } 
    Matrix3(size_t x, size_t y, size_t z) ;
    size_t size() const { return d1*d2*d3 ; }
    size_t dim1() const { return d1 ; }
    size_t dim2() const { return d2 ; }
    size_t dim3() const { return d3 ; }
    std::ostream& Print(std::ostream& s) const ;
    std::istream& Input(std::istream& s) ;

    Slice_iter<double> first(size_t i, size_t j) ;
    Cslice_iter<double> first(size_t i, size_t j) const ;
    Slice_iter<double> second(size_t i, size_t j) ;
    Cslice_iter<double> second(size_t i, size_t j) const ;
    Slice_iter<double> third(size_t i, size_t j) ;
    Cslice_iter<double> third(size_t i, size_t j) const ;

    double& operator() (size_t x, size_t y, size_t z) ;
    double operator() (size_t x, size_t y, size_t z) const ;

    //  Slice_iter<double> operator[] (size_t i) { return row(i) ; }
    //  Cslice_iter<double> operator[] (size_t i) const { return row(i) ; }

    Matrix3& operator*=(double) ;
    valarray<double>& array() { return v ; }

    void resize(size_t x, size_t y, size_t z) ;

    Matrix3 &operator=(const Matrix3 &mt) {
      if(this!= &mt) 
        resize(mt.dim1(),mt.dim2(),mt.dim3()) ;

      for(size_t k=0;k<d3;k++)
        for(size_t j=0;j<d2;j++)
          for(size_t i=0;i<d1;i++)
            (*this)(i,j,k) = mt(i,j,k) ;

      return *this ;
    }

  } ;

  inline Slice_iter<double> Matrix3::first(size_t i, size_t j)
  {
    return Slice_iter<double> (&v,slice(i*d1+j*d12,d1,1)) ;
  }

  inline Cslice_iter<double> Matrix3::first(size_t i, size_t j) const
  {
    return Cslice_iter<double> (const_cast<valarray<double> *>(&v),
                                slice(i*d1+j*d12,d1,1)) ;
  }

  inline Slice_iter<double> Matrix3::second(size_t i, size_t j)
  {
    return Slice_iter<double> (&v,slice(i+j*d12,d2,d1)) ;
  }

  inline Cslice_iter<double> Matrix3::second(size_t i, size_t j) const
  {
    return Cslice_iter<double> (const_cast<valarray<double> *>(&v),
                                slice(i+j*d12,d2,d1)) ;
  }

  inline Slice_iter<double> Matrix3::third(size_t i, size_t j)
  {
    return Slice_iter<double> (&v,slice(i+j*d1,d3,d12)) ;
  }

  inline Cslice_iter<double> Matrix3::third(size_t i, size_t j) const
  {
    return Cslice_iter<double> (const_cast<valarray<double> *>(&v),
                                slice(i+j*d12,d2,d1)) ;
  }

  inline std::ostream &operator<<(std::ostream& s, Matrix3& m) {
    return m.Print(s) ;
  }

  inline std::istream &operator>>(std::istream& s, Matrix3& m) {
    return m.Input(s) ;
  }

}
namespace Loci {
  using namespace fluidPhysics ;
  class Matrix2Converter {
    Matrix2 &ref ;
  public:
    Matrix2Converter(Matrix2 &iref) : ref(iref) {}
    int getSize() const {
      return 0 ;
    }
    void getState(char *buf, int &size) {
      size = getSize() ;
      std::cerr << "Matrix2Converter: getState not implemented" << std::endl ;
    }
    void setState(char *buf, int size) {
      std::cerr << "Matrix2Converter: setState not implemented" << std::endl ;
    }
  } ;

  template<> struct data_schema_traits<Matrix2> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef char Converter_Base_Type ;
    typedef Matrix2Converter Converter_Type ;
  } ;

}
#endif


