#include "Matrix.h"
#include <Tools/stream.h>
#include <valarray> 
#include <Loci.h>
#include <Tools/parse.h>
#include "fluidConst.h"

using std::valarray ;
using std::slice ;

namespace fluidPhysics {
  Matrix2::Matrix2(size_t x, size_t y) {
    d1 = x ;
    d2 = y ;
    v.resize(x*y) ;
  }

  void Matrix2::resize(size_t x, size_t y) {
    d1 = x ;
    d2 = y ;
    v.resize(x*y) ;
  }

  double& Matrix2::operator() (size_t x, size_t y) {
    return column(y) [x] ;
  }

  double  Matrix2::operator() (size_t x, size_t y) const {
    return column(y) [x] ;
  }

  double mul(const valarray<double>& v1, const valarray<double>& v2) {
    double res = 0 ;
    for(size_t i=0; i<v1.size(); i++)
      res+=v1[i]*v2[i] ;
    return res ;
  }

  /*
    valarray<double> operator*(const Matrix2&m, const valarray<double>& v) {
    valarray<double> res(m.dim1()) ;
    for(int i=0; i<m.dim1(); i++) res(i) = mul(m.row(i),v) ;
    return res ;
    }
  */

  Matrix2& Matrix2::operator*=(double d) {
    (v) *=d ;
    return *this ;
  }

  std::ostream &Matrix2::Print(std::ostream &s) const {
    s << d1 << " " << d2 << endl ;
    for(size_t i=0;i<d1;++i) {
      for(size_t j=0;j<d2;++j) 
        s << (*this)(i,j) << " " ;
      s << endl ;
    }
    return s ;
  }

  std::istream &Matrix2::Input(std::istream &s) {
    size_t n = size() ;
    size_t i = 0 ;
    while(i<n) {
      double f ;
      Loci::parse::kill_white_space(s) ;
      if(!Loci::parse::is_real(s))
        break ;
      f=Loci::parse::get_real(s) ;
      (v)[i++] = f ;
      //    s>>(v)[i++] ;
    }
    if(i != n) {
      cerr << "too few elements on Matrix input" << endl ;
      cout << "i= " << i << " " << "n= " << n << endl ;
    }
    return s ;
  }

  void Matrix2::ludcmp(int *indx, double &d) {
    if(d1 != d2) {
      cerr << "Matrix cannot be LU decomposed for d1 != d2\n" ;
      exit(-1) ;
    }
    Matrix2 &s = (*this) ;
    int n = d1 ;
    double *vv = new double[n] ;
    d = 1.0 ;
    double big, dum,sum,temp ;
    int imax = 0 ;
    for(int i=0;i<n;++i) {
      big = 0.0 ;
      for(int j=0;j<n;++j)
        if((temp=fabs(s(i,j))) > big) 
          big = temp ;
      if(big==0.0) {
        cerr << "Singular Matrix doesn't have LU decomposition\n" ;
        exit(-1) ;
      }  
      vv[i] = 1.0/big ;
    }
    for(int j=0;j<n;++j) {
      for(int i=0;i<j;++i) {
        sum = s(i,j) ;
        for(int k=0;k<i;++k) 
          sum -= s(i,k)*s(k,j) ;
        s(i,j) = sum ;
      }
      big = 0.0 ;
      for(int i=j;i<n;i++) {
        sum = s(i,j) ;
        for(int k=0;k<j;k++)
          sum -= s(i,k)*s(k,j) ;
        s(i,j) = sum ;
        if((dum=vv[i]*fabs(sum)) >= big) {
          big = dum ;
          imax = i ;
        }
      }   
      if(j!=imax) {
        for(int k=0;k<n;++k) {
          dum = s(imax,k) ;
          s(imax,k) = s(j,k) ;
          s(j,k) = dum ;
        }
        d = -d ;
        vv[imax] = vv[j] ;
      }

      indx[j] = imax ;
      if(fabs(s(j,j))<EPSILON) {
        cerr << "The matrix is singular\n" ;
        s(j,j) = EPSILON ;
      }
      if(j!=(n-1)) {
        dum = 1.0/s(j,j) ;
        for(int i=j+1;i<n;++i)
          s(i,j) *= dum ;
      } 
    }
    delete[] vv ;
  }

  void Matrix2::lubksb(int *indx, double *b) {
    int ii=0,ip ;
    double sum ;
    int n = d1 ;
    Matrix2 & s = (*this) ; 
    for(int i=0;i<n;++i) {
      ip=indx[i] ;
      sum=b[ip] ;
      b[ip]=b[i] ;
      if(ii)
        for(int j=ii;j<i;++j) 
          sum -= s(i,j)*b[j] ;
      else if (sum) ii = i ;
      b[i] = sum ;
    }

    for(int i= n-1;i>=0;i--) {
      sum = b[i] ;
      for(int j=i+1;j<n;++j)
        sum -= s(i,j)*b[j] ;
      b[i] = sum/s(i,i) ;
    }
  }

  void Matrix2::inverse(Matrix2 &y) {
    if(d1 != d2) {
      cerr << "Matrix doesn't have inverse for d1 != d2\n" ;
      exit(-1) ;
    }
   
    int n= d1 ;
    int *indx = new int[n] ;
    double *col = new double[n] ;
    double d ;
    ludcmp(indx,d) ;
    for(int j=0;j<n;++j) {
      for(int i=0;i<n;++i) 
        col[i] = 0.0 ;
      col[j] = 1.0 ;
      lubksb(indx,col) ; 
      for(int i=0;i<n;++i)
        y(i,j) = col[i] ;
    }
  } 

  void Matrix2::QR() {
    if(d1<d2) {
      cerr << "Matrix cannot be QR factorized for d1<d2\n" ;
      exit(-1) ;
    }
    /*
      if(d1==d2) {
      Matrix2 y(d1,d1) ;
      inverse(y) ;
      for(int j=0;j<d1;++j) 
      for(int i=0;i<d1;++i)
      (*this)(i,j) = y(j,i) ;
      return ;
      }
    */
    Matrix2 r(d2,d2) ;
    Matrix2 s(d1,d2) ;
    s = *this ;
    for(size_t k=0;k<d2;++k) {
      r(k,k) = 0.0 ;
      for(size_t i=0;i<d1;++i)
        r(k,k) += s(i,k)*s(i,k) ;
      r(k,k) = sqrt(r(k,k)) ;
      if(r(k,k) < EPSILON) {
        cerr << "Ill defined matrix at column " << k+1 << endl ;
        cout << *this << endl ;
        cout << s << endl ;
      
        for(size_t i=0;i<d1;++i)
          s(i,k) = 0.0 ;
        exit(-1) ;
      } else {
        for(size_t i=0;i<d1;++i)
          s(i,k) /= r(k,k) ;
      }
      for(size_t j=k+1;j<d2;++j) {
        r(k,j) = 0.0 ;
        for(size_t i=0;i<d1;++i) 
          r(k,j) += s(i,k)*s(i,j) ;
        for(size_t i=0;i<d1;++i) 
          s(i,j) -= r(k,j) * s(i,k) ;
      }
    }

    //to get inverse of r matrix
    Matrix2 ri(d2,d2) ;
    for(int j=0;j<int(d2);++j) 
      for(int i=int(d2)-1;i>=0;--i) {
        if(i==j) ri(i,j) = 1.0 ;
        else ri(i,j) = 0.0 ;
        for(size_t k=i+1;k<d2;++k) 
          ri(i,j) -= r(i,k)*ri(k,j) ;
        if(fabs(r(i,i)) < EPSILON) {
          cerr << "Ill defined upper D matrix at " << i+1 << endl ;
          cout << *this << endl ;
          cout << r << endl ;
          exit(-1) ;
          ri(i,j) = 0.0 ;
        } else 
          ri(i,j) /= r(i,i) ;
      }
  
    // put ri*s(T) to *this ;
    for(size_t i=0;i<d2;++i)  
      for(size_t j=0;j<d1;++j) {
        (*this)(j,i) = 0.0 ;
        for(size_t k=0;k<d2;++k) 
          (*this)(j,i) += ri(i,k)*s(j,k) ;
      }
  }  

  Matrix3::Matrix3(size_t x, size_t y, size_t z) {
    d1 = x ;
    d2 = y ;
    d3 = z ;
    d12 = x*y ;
    v.resize(x*y*z) ;
  }

  void Matrix3::resize(size_t x, size_t y, size_t z) {
    d1 = x ;
    d2 = y ;
    d3 = z ;
    d12 = x*y ;
    v.resize(x*y*z) ;
  }

  double& Matrix3::operator() (size_t x, size_t y, size_t z) {
    return first(y,z)[x] ;
  }

  double Matrix3::operator() (size_t x, size_t y, size_t z) const {
    return first(y,z)[x] ;
  }

  Matrix3& Matrix3::operator*=(double d) {
    (v) *=d ;
    return *this ;
  }

  std::ostream &Matrix3::Print(std::ostream &s) const {
    size_t n = size() ;
    size_t i = 0 ;
    s << d1 << " " << d2 << " " << d3 << endl ;
    while(i<n && s<<(v)[i++]) s << " " ;
    return s ;
  }

  std::istream &Matrix3::Input(std::istream &s) {
    size_t n = size() ;
    size_t i = 0 ;
    while(i<n) {
      double f ;
      Loci::parse::kill_white_space(s) ;
      if(!Loci::parse::is_real(s))
        break ;
      f=Loci::parse::get_real(s) ;
      (v)[i++] = f ;
    }
    if(i != n) {
      cerr << "too few elements on Matrix input" << endl ;
      cout << "i= " << i << " " << "n= " << n << endl ;
    }
    return s ;
  }
}
