#ifndef TRANSPORT_DB_H
#define TRANSPORT_DB_H

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include "Matrix.h"
#include <Tools/tools.h> 
#include <Loci.h>

namespace fluidPhysics {
  class transport_db {
  public:
    transport_db() { species_count = 0 ; number_order = 0 ; } ;
    std::istream &Input(std::istream &s) ;
    std::ostream &Print(std::ostream &s) const ;

    void species_reorder(const std::vector<std::string> &species_order) ;
    //Function species_reorder changes the species order within transport model to
    //that of chemistry model

    double MolecularWeight(int i) const { return mwt[i] ; }

    double mcavis(double T, double* mf) const;
    //Function mcavis computes the mixture viscosity,
    //given temperature T and species mole fractions mf.
    //It uses modification of the Wilke semi-empirical formulas.

    double mcacon(double T, double* mf) const;
    //Function mcacon computes the mixture thermal conductivity,
    //given temperature T and species mole fractions mf.

    void mcadif(double P, double T, double* mf, double* d) const;
    //Function mcadif computes mixture-averaged diffusion coefficients d
    //given pressure P, temperature T, and spcies mole fractions mf.
   
  private:
    int species_count, number_order ;
    std::map<std::string,int> species_map ;
    std::vector<double> mwt ;
    Matrix2 viscof ;
    Matrix2 condcof ;
    double patmos ;
    Matrix3 difcof ;

    void mceval(const double alogt, const int nk, const int no, 
                const Matrix2& cof, double* val) const ;
    //Function mceval uses Horners algorithm to evaluate a polynomial fit.
  } ;

  inline std::ostream &operator<<(std::ostream &s, const transport_db &tdb)
  { return tdb.Print(s) ; }

  inline std::istream &operator>>(std::istream &s, transport_db &tdb)
  { return tdb.Input(s) ; }

  struct Sutherland_param { 
    double a1,a2,a3 ; //coefficients in Sutherland transport model
    double k1,k2,k3 ;
    double pr ;
    bool usepr ;
    Sutherland_param() ;
    std::istream &Input(std::istream &s) ;
    std::ostream &Print(std::ostream &s) const ;
  } ;

  inline std::ostream & operator<<(std::ostream &s, const Sutherland_param 
                                   &suther) 
  {return suther.Print(s) ;}
  inline std::istream & operator>>(std::istream &s, Sutherland_param &suther) 
  {return suther.Input(s) ;}

  struct powerLaw_param {
    double mu_ref, T_ref;
    double power ;
    double Pr ;
    powerLaw_param() {mu_ref = 1.8e-5; T_ref = 300; power = 0.7; Pr = 0.7; }
  } ;

}
namespace Loci {
  template<> struct data_schema_traits<fluidPhysics::transport_db> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<fluidPhysics::transport_db> Converter_Type ;
  } ;
  
  template<> struct data_schema_traits<fluidPhysics::Sutherland_param> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(fluidPhysics::Sutherland_param()) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::Sutherland_param,a1) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::Sutherland_param,a2) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::Sutherland_param,a3) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::Sutherland_param,k1) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::Sutherland_param,k2) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::Sutherland_param,k3) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::Sutherland_param,pr) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::Sutherland_param,usepr) ;
      return DatatypeP(ct) ;
    }
  } ;

  template<> struct data_schema_traits<fluidPhysics::powerLaw_param> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(fluidPhysics::powerLaw_param()) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::powerLaw_param,mu_ref) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::powerLaw_param,T_ref) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::powerLaw_param,power) ;
      LOCI_INSERT_TYPE(ct,fluidPhysics::powerLaw_param,Pr) ;
      return DatatypeP(ct) ;
    }
  } ;
}
#endif

