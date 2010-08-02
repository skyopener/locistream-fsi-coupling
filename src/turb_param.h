#ifndef TURB_PARAM_H
#define TURB_PARAM_H

#include <iostream>
#include "sciTypes.h"
#include <Tools/tools.h>
#include <Loci.h>

namespace streamUns {
  struct turbulent_transport {
    real TurbulentPrandtlNumber ;
    real TurbulentSchmidtNumber ;
    turbulent_transport() {
      TurbulentPrandtlNumber = 0.9 ;
      TurbulentSchmidtNumber = 0.7 ;
    }
  } ;
      
  struct Spa_All_param {
    real cb1,k_coeff,cv1,cb2,cw2,cw3,sigma ; //coefficients in Spalart_Allmaras
    //                                         //model

    Spa_All_param()   {
      cb1=0.14; //production coeffient
      k_coeff=0.41; // coeffient for constructing vorticity \tilde{S}
      cv1=7.1;  //coeffient for constructing f_v1 which is the ratio of
      //        //eddy viscoisty to molecular viscosity
      cb2=0.622; // coeffient in diffusion term
      cw2=0.3; // coeffient in destruction term
      cw3=2.0; //coeffienct in destruction term
      sigma=0.6667; //coeffienct in diffusion term
    }
  } ;

  struct sst_param { 
    real sigmak, sigmae, beta, gama ;
  } ;

  struct sst1_param {
    real sigmak, sigmae, beta, betas, kappa, a1 ;
    sst1_param() {
      //sigmak = 0.5 ;
      sigmak = 0.85 ;
      sigmae = 0.5 ;
      beta = 0.075 ;
      betas = 0.09 ;
      kappa = 0.41 ;
      a1 = 0.31 ;
    }
  } ;

  struct sst2_param {
    real sigmak, sigmae, beta, betas, kappa ;
    sst2_param() {
      sigmak = 1.0 ;
      sigmae = 0.856 ;
      beta = 0.0828 ;
      betas = 0.09 ;
      kappa = 0.41 ;
    }
  } ;

  struct rke_param {
    real a1,ce1,ce2,sigmak,sigmae,ctau,amu,ae,aet,cs ;
    rke_param() {
      a1 = 1.25 ;
      ce1 = 1.44 ;
      ce2 = 1.92 ;
      sigmak = 1.0 ;
      sigmae = 1.3 ;
      ctau = 1.414213562 ; //sqrt(2.0)
      amu = 0.01 ;
      ae = 0.3 ;
      aet = 0.15 ;
      cs = 0.05 ;
    }
  } ;

}

namespace Loci {
  template<> struct data_schema_traits<streamUns::turbulent_transport> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(streamUns::turbulent_transport()) ;
      LOCI_INSERT_TYPE(ct,streamUns::turbulent_transport,TurbulentPrandtlNumber) ;
      LOCI_INSERT_TYPE(ct,streamUns::turbulent_transport,TurbulentSchmidtNumber) ;
      return DatatypeP(ct) ;
    }
  } ;
                       
  template<> struct data_schema_traits<streamUns::Spa_All_param> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(streamUns::Spa_All_param()) ;
      LOCI_INSERT_TYPE(ct,streamUns::Spa_All_param,cb1) ;
      LOCI_INSERT_TYPE(ct,streamUns::Spa_All_param,k_coeff) ;
      LOCI_INSERT_TYPE(ct,streamUns::Spa_All_param,cv1) ;
      LOCI_INSERT_TYPE(ct,streamUns::Spa_All_param,cb2) ;
      LOCI_INSERT_TYPE(ct,streamUns::Spa_All_param,cw2) ;
      LOCI_INSERT_TYPE(ct,streamUns::Spa_All_param,cw3) ;
      LOCI_INSERT_TYPE(ct,streamUns::Spa_All_param,sigma) ;
      return DatatypeP(ct) ;
    }
  } ;

  template<> struct data_schema_traits<streamUns::sst_param> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(streamUns::sst_param()) ;
      LOCI_INSERT_TYPE(ct,streamUns::sst_param,sigmak) ;
      LOCI_INSERT_TYPE(ct,streamUns::sst_param,sigmae) ;
      LOCI_INSERT_TYPE(ct,streamUns::sst_param,beta) ;
      LOCI_INSERT_TYPE(ct,streamUns::sst_param,gama) ;
      return DatatypeP(ct) ;
    }
  } ;
    

  template<> struct data_schema_traits<streamUns::sst1_param> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(streamUns::sst1_param()) ;
      LOCI_INSERT_TYPE(ct,streamUns::sst1_param,sigmak) ;
      LOCI_INSERT_TYPE(ct,streamUns::sst1_param,sigmae) ;
      LOCI_INSERT_TYPE(ct,streamUns::sst1_param,beta) ;
      LOCI_INSERT_TYPE(ct,streamUns::sst1_param,betas) ;
      LOCI_INSERT_TYPE(ct,streamUns::sst1_param,kappa) ;
      LOCI_INSERT_TYPE(ct,streamUns::sst1_param,a1) ;
      return DatatypeP(ct) ;
    }
  } ;

  template<> struct data_schema_traits<streamUns::sst2_param> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(streamUns::sst2_param()) ;
      LOCI_INSERT_TYPE(ct,streamUns::sst2_param,sigmae) ;
      LOCI_INSERT_TYPE(ct,streamUns::sst2_param,beta) ;
      LOCI_INSERT_TYPE(ct,streamUns::sst2_param,betas) ;
      LOCI_INSERT_TYPE(ct,streamUns::sst2_param,kappa) ;
      return DatatypeP(ct) ;
    }
  } ;
  
}

#endif
