#include <Loci.h>
#include <Tools/stream.h>
#include <Tools/parse.h>

#include "transport_db.h"
#include "Matrix.h"
#include "periodic_table.h"
#include "fluidConst.h"
#include "scratch_array.h"

#include <vector>

namespace fluidPhysics {
  using namespace Loci::parse ;

  using std::vector ;
  using std::map ;

  void transport_db::species_reorder(const vector<std::string> &species_order) {
    //This function changes the species order within transport model to 
    //that of chemistry model 
    map<int,int> order ;
    if(int(species_order.size()) != species_count) {
      cerr << " Transport species number is incorrect " << endl ;
      cerr << " Unable to use chemkin transort model!" << endl ;
      Loci::Abort() ;
    }

    if(species_count>1) {
      typedef std::map<std::string,int>::iterator IT ;
      for(IT p=species_map.begin();p!=species_map.end();++p) {
        int idx = -1 ;
        for(size_t j=0;j<species_order.size();++j)
          if(species_order[j] == p->first) {
            idx = j ;
            break ;
          }
        if(idx == -1) {
          cerr << "species " << p->first << " is not in chemistry model!"
               << endl ;
          cerr << "unable to use chemkin transport model!" << endl ;
          Loci::Abort() ;
        }
        order[p->second] = idx ;
        p->second = order[p->second] ;
      }

      std::vector<double> tmp1(species_count) ;
      for(int i=0;i<species_count;++i) 
        tmp1[i] = mwt[i] ;
      for(int i=0;i<species_count;++i) 
        mwt[order[i]] = tmp1[i] ;


      Matrix2 tmp(number_order,species_count) ;
      for(int i=0;i<species_count;++i)
        for(int j=0;j<number_order;++j) 
          tmp(j,i) = viscof(j,i) ;
      for(int i=0;i<species_count;++i)
        for(int j=0;j<number_order;++j) 
          viscof(j,order[i]) = tmp(j,i) ;

      for(int i=0;i<species_count;++i)
        for(int j=0;j<number_order;++j) 
          tmp(j,i) = condcof(j,i) ;
      for(int i=0;i<species_count;++i)
        for(int j=0;j<number_order;++j)
          condcof(j,order[i]) = tmp(j,i) ;

      Matrix3 tmpd(number_order,species_count,species_count) ;
      for(int i=0;i<species_count;++i)
        for(int j=0;j<species_count;++j)
          for(int k=0;k<number_order;++k)
            tmpd(k,j,i) = difcof(k,j,i) ;
      for(int i=0;i<species_count;++i)
        for(int j=0;j<species_count;++j)
          for(int k=0;k<number_order;++k)
            difcof(k,order[j],order[i]) = tmpd(k,j,i) ;
    }
  }

  void transport_db::mceval(const double tf, const int nk, const int no, 
                            const Matrix2& cof, double* val) const {
    //This function uses Horners algorithm to evaluate a polynomial fit. 
    double b ;
    int nom1 = no - 1 ;
    for(int k=0;k<nk;++k) {
      b = cof(nom1,k) ;
      for(int i=0;i<nom1;++i) 
        b = cof(nom1-i-1,k) + b*tf ;  
      val[k] = b ; 
    }
  }

  double transport_db::mcavis(double T, double* mf) const {
    //This function computes the mixture viscosity, 
    //given temperature T and species mole fractions mf.
    //It uses modification of the Wilke semi-empirical formulas.

    if(species_count==1 && (species_map.begin()->first)=="_Air") {
      cerr << "Use sutherlands model with 1 species air model!" << endl ;
      cerr << "chemkin not supported for _Air species" << endl ;
      Loci::Abort() ;
      return 1.458e-6*pow(T,1.5)/(T+110.4) ;
    }
    if(T<=0) {
      cerr << "Temperature became negative, T = " << T << endl ;
      return 0 ;
    }
    double alogt = log(T) ;
    scratch_array<double> vis(species_count) ;
  
    mceval(alogt,species_count,number_order,viscof,vis) ;

    if(species_count==1) {
      double retval = 0.1*exp(vis[0]) ;
      return retval ;
    }
    for(int k=0;k<species_count;++k)  vis[k] = exp(vis[k]) ;
    double sumo = 0.0 ;
    double sumi, top, bot, phikj ;
    for(int k=0;k<species_count;++k) {
      sumi = 0.0 ;
      for(int j=0;j<species_count;++j) {
        top = pow(1.0+sqrt(vis[k]/vis[j])*pow(mwt[j]/mwt[k],0.25),2) ;
        bot = sqrt(1.0+mwt[k]/mwt[j]) ;
        phikj = 0.3535534*top/bot ;  
        sumi += mf[j]*phikj ;
      }
      sumo += mf[k]*vis[k]/sumi ;
    }
    return 0.1*sumo ; //to change unit
  }

  double transport_db::mcacon(double T, double* mf) const {
    //This function computes the mixture thermal conductivity,
    //given temperature T and species mole fractions mf.

    if(species_count==1 && (species_map.begin()->first)=="_Air") {
      cerr << "Use sutherlands model with 1 species air model!" << endl ;
      cerr << "chemkin not supported for _Air species" << endl ;
      Loci::Abort() ;
      return 1.967e-3*pow(T,1.5)/(T+110.4) ;
    }
    double alogt = log(T) ;
    scratch_array<double> con(species_count) ;
    
    mceval(alogt,species_count,number_order,condcof,con) ;
    if(species_count==1) {
      double retval = 1e-5*exp(con[0]) ;
      return retval ;
    }

    double sum = 0.0 ;
    double sumr = 0.0 ;
  
    for(int k=0;k<species_count;++k) {
      con[k] = exp(con[k]) ;

      sum += mf[k]*con[k] ;
      sumr += mf[k]/con[k] ;
    }

    return 5e-6*(sum+1.0/sumr) ;
  }

  void transport_db::mcadif(double P, double T, double* mf, double* d) const {
    //This function computes mixture-averaged diffusion coefficients
    //given pressure P, temperature T, and spcies mole fractions mf.

    if(species_count==1) {
      d[0] = 0.0 ;
      return ;
    }
    double alogt = log(T) ;

    Matrix2 djk(species_count,species_count) ;
    scratch_array<double> djkj(species_count) ;

    Matrix2 difcofk(number_order,species_count) ;

    for(int k=0;k<species_count;++k) {
      for(int j=0;j<species_count;++j) 
        for(int i=0;i<number_order;++i) 
          difcofk(i,j) = difcof(i,j,k) ;
      mceval(alogt,species_count,number_order,difcofk,djkj) ; 
      for(int j=0;j<species_count;++j)
        djk(j,k) = djkj[j] ;
    }


    for(int k=0;k<species_count;++k) 
      for(int j=0;j<species_count;++j) 
        djk(j,k) = exp(djk(j,k)) ;

    double wts = 0.0 ;
    for(int k=0;k<species_count;++k) 
      wts += mwt[k]*mf[k] ;
 
    double sumxw, sumxod ;
    for(int k=0;k<species_count;++k) {
      sumxw = 0.0 ;
      sumxod = 0.0 ;
      for(int j=0;j<species_count;++j) 
        if(j!=k) {
          sumxw += mf[j]*mwt[j] ;
          sumxod += mf[j]/djk(j,k) ;
        } ;
      if(sumxod < EPSILON)
        d[k] = 0 ;
      else
        d[k] = sumxw/(wts*sumxod) ;
    }

    for(int k=0;k<species_count;++k)
      d[k] *= (1e-4*patmos/P) ;    
 
    return ;

  }

  istream &transport_db::Input(istream &s)
  {
    //This function reads in the transport coeffecients from a local file.
    kill_white_space(s) ;
    while(s.good() && !s.eof() && s.peek() != char_traits<char>::eof()) {
      kill_white_space(s) ;
      if(get_token(s,"species_count")) {
        kill_white_space(s) ;
        if(s.peek() != '=') {
          cerr << "error reading transport_db" << endl ;
          Loci::Abort() ;
          return s ;
        }
        s.get() ;
        kill_white_space(s) ;
        s >> species_count ;
      } else {
        if(s.peek() != char_traits<char>::eof()) {
          cerr << "unexpected input in transport_db" << endl ;
          Loci::Abort() ;
          return s ;
        }
      }
    
      kill_white_space(s) ;
      if(get_token(s,"species")) {
        kill_white_space(s) ;
        if(s.peek() != ':') {
          cerr << "error reading transport_db" << endl ;
          Loci::Abort() ;
          return s ;
        }
        s.get() ;
        for(int i=0;i<species_count;++i) {
          kill_white_space(s) ;
          if(!is_name(s)) {
            cerr << "error reading transport_db" << endl ;
            Loci::Abort() ;
            break ;
          }
          string name = get_name(s) ;
          species_map[name] = i ;
        }
      } else {
        if(s.peek() != char_traits<char>::eof()) {
          cerr << "unexpected input in transport_db" << endl ;
          Loci::Abort() ;
          return s ;
        }
      }

      if(species_count==1 && (species_map.begin()->first)=="_Air")
        return s ;

      kill_white_space(s) ;
      if(get_token(s,"number_order")) {
        kill_white_space(s) ;
        if(s.peek() != '=') {
          cerr << "error reading transport_db" << endl ;
          Loci::Abort() ;
          return s ;
        }
        s.get() ;
        kill_white_space(s) ;
        s >> number_order ;
      } else {
        if(s.peek() != char_traits<char>::eof()) {
          cerr << "unexpected input in transport_db" << endl ;
          Loci::Abort() ;
          return s ;
        }
      }

      kill_white_space(s) ;
      if(get_token(s,"patmos")) {
        kill_white_space(s) ;
        if(s.peek() != '=') {
          cerr << "error reading transport_db" << endl ;
          Loci::Abort() ;
          return s ;
        }
        s.get() ;
        kill_white_space(s) ;
        s >> patmos ;
      } else {
        if(s.peek() != char_traits<char>::eof()) {
          cerr << "unexpected input in transport_db" << endl ;
          Loci::Abort() ;
          return s ;
        }
      }

      patmos *= 1.01e5 ;

      kill_white_space(s) ;
      mwt.resize(species_count) ;

      for(int i=0;i<species_count;++i) s >> mwt[i] ;
      kill_white_space(s) ;

      viscof.resize(number_order,species_count) ;
      condcof.resize(number_order,species_count) ;
      difcof.resize(number_order,species_count,species_count) ;

      s >> condcof ;
      s >> viscof ;
      s >> difcof ;
      kill_white_space(s) ;
    }

    return s ;
  }

  std::ostream &transport_db::Print(std::ostream &s) const
  {
    //This function prints out the transport coefficients.
    s << "species_count " << species_count << endl ;
    typedef std::map<std::string,int>::const_iterator CI ;
    for(CI p=species_map.begin();p!=species_map.end();++p)
      s << p->first << " " << p->second << endl ;
    if(species_count==1 && (species_map.begin()->first)=="_Air") {
      cout << " Sutherland model is used for ideal air transport property " << endl ;
      return s ;
    }
    s << "patmos= " << patmos << endl ;
    s << "molecular weights" << ":" << endl ;
    for(int i=0;i<species_count;++i) s << mwt[i] << " " ;
    s << endl ;
    s << "viscous coefficient" << ":" << endl ;
    viscof.Print(s) ;
    s << endl ;
    for(int k=0;k<species_count;++k)  {
      for(int i=0;i<number_order;++i) 
        s << viscof(i,k) << " ";
      s << endl ;
    }
    s << endl ;
    s << "conductivity coefficient" << ":" << endl ;
    condcof.Print(s) ;
    s << endl ; 
    for(int k=0;k<species_count;++k)  {
      for(int i=0;i<number_order;++i)
        s << condcof(i,k) << " ";
      s << endl ;
    }
    s << endl ;
    s << "diffusion coefficients" << ":" << endl ;
    difcof.Print(s) ;
    s << endl ;
    for(int k=0;k<species_count;++k)  {
      for(int j=0;j<species_count;++j)  {
        for(int i=0;i<number_order;++i)
          s << difcof(i,j,k) << " ";
        s << endl ;
      }
    }

    s << endl ;
 
    return s ;
  }

  class sutherol : public options_list {
  public:
    sutherol() : options_list("a1:a2:a3:k1:k2:k3:pr") {}
  } ;

  Sutherland_param::Sutherland_param() 
  {
    a1=1.458e-6 ;
    a2=1.5 ;
    a3=110.4 ;
    k1=0. ;
    k2=0. ;
    k3=0. ;
    pr = .75 ;
    usepr = true ;
  }

  istream &Sutherland_param::Input(istream &s) {
    sutherol finput ;
    s >> finput ;
    a1=1.458e-6 ;
    a2=1.5 ;
    a3=110.4 ;

    if(finput.optionExists("a1")) {
      finput.getOption("a1",a1) ;
    } 
    if(finput.optionExists("a2")) {
      finput.getOption("a2",a2) ;
    }
    if(finput.optionExists("a3")) {
      finput.getOption("a3",a3) ;
    }
    if(finput.optionExists("pr")) {
      usepr = true ;
      finput.getOption("pr",pr) ;
    } else {
      usepr = false ;
      k1=2.495e-3 ;
      k2=1.5 ;
      k3=194.0 ;
      if(finput.optionExists("k1"))
        finput.getOption("k1",k1) ;
      if(finput.optionExists("k2"))
        finput.getOption("k2",k2) ;
      if(finput.optionExists("k3"))
        finput.getOption("k3",k3) ;
    }
    return s ;
  }

  ostream &Sutherland_param::Print(ostream &s) const {
    if(usepr)
      s<<"<a1=" <<a1<<", a2=" <<a2<<",a3= "<<a3<< ",pr="<<pr<<">"<<endl ;
    else
      s<<"<a1=" <<a1<<", a2=" <<a2<<",a3= "<<a3<<
        ",k1="<<k1<<",k2="<<k2<<",k3="<<k3<<">"<<endl ;
      
    return s ;
  }


}
