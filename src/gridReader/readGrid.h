#ifndef READ_GRID_H
#define READ_GRID_H

// System includes.
#include <iostream>
#include <string>
using std::ostream ;
using std::string ;

// Loci includes.
#include <Loci.h>

// StreamUns includes.
#include "name_var.h"

namespace streamUns {

  typedef double real_t ; // Used only in CreateEdges().

  void ReadGrid(fact_db &facts,const rule_db &rdb,const char *caseName) ;

  class grid_options : public options_list {
    public:
    grid_options() : options_list("axisymmetric:pieSlice:file_type:Lref:translate:bc_check") {} ;
  } ;

  enum matrix_coloring_type {COLOR_DEFAULT, COLOR_DFS} ;
 
  bool CheckBoundaryConditions(fact_db &facts) ;
  void ColorMatrix(fact_db &facts,matrix_coloring_type mct) ;
  void CreateCellStencil(fact_db &facts) ;
  void CreateGhostCells(fact_db &facts) ;
  void CreateLowerUpper(fact_db &facts) ;
  void CreateEdges(fact_db &facts) ;
  void FindMindNoslip(fact_db &facts) ;
  void MakeFacesConsistent(fact_db &facts) ;
  void SetupBoundaryConditions(fact_db &facts) ;
  void SortFaces(fact_db &facts,bool doSort) ;
}

namespace Loci {
  template<> struct data_schema_traits<streamUns::grid_options> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<streamUns::grid_options> Converter_Type ;
  } ;
}

namespace streamUns {

  class BC_Check : public Loci::CPTR_type {
    public:

      // Returns string containing comma separated list of boundary
      // condition tags that are checked by this checker
      virtual string BoundaryConditions() = 0 ;

      // Returns string containing variables that are tested by this
      // checker.
      virtual string VariablesChecked(fact_db &facts) = 0 ;

      // Check the boundary condition
      virtual bool CheckOptions(const options_list &bc_options,fact_db
        &facts) = 0 ;

      // Emit Error Message
      virtual ostream  &ErrorMessage(ostream &s) = 0 ;
  } ;

  typedef Loci::CPTR<BC_Check> BC_implP ;

  class register_BC_type {
  public:
    virtual ~register_BC_type() {}
    virtual BC_implP get_BC() const = 0 ;
  } ;

  class register_BC_impl_list ;

  class BC_impl_list {
  public:
    class BC_list_iterator ;
    friend class BC_list_iterator ;
    friend class register_BC_impl_list ;
    class BC_list_ent {
    public:
      BC_list_ent(register_BC_type *p, BC_list_ent *nxt) :
        rr(p), next(nxt) {}
      register_BC_type *rr ;
      BC_list_ent *next ;
    } ;
    BC_list_ent *list ;
  public:
    class BC_list_iterator {
      BC_list_ent *p ;
    public:
      BC_list_iterator(BC_list_ent *ptr) : p(ptr) {}
      BC_implP operator*() { return p->rr->get_BC() ; }
      BC_list_iterator &operator++() {
        p = p->next ;
        return *this ;
      }
      BC_list_iterator operator++(int ) {
        BC_list_iterator tmp(p) ;
        p = p->next ;
        return tmp ;
      }
      BC_list_ent* get_p() { return p; }
      bool operator==(const BC_list_iterator &i) { return i.p == p ; }
      bool operator!=(const BC_list_iterator &i) { return i.p != p ; }
    } ;
    typedef BC_list_iterator iterator ;
    BC_impl_list() {list = 0 ; }
    ~BC_impl_list() ;

    void push_BC(register_BC_type *rr) ;
    iterator begin() { return iterator(list) ; }
    iterator end() { return iterator(0) ; }
    void clear() ;
    void copy_BC_list(const BC_impl_list &rl) ;
    void copy_BC_list(const register_BC_impl_list &rl) ;
    BC_impl_list(const BC_impl_list &x) {
      list = 0 ;
      copy_BC_list(x) ;
    }
    BC_impl_list &operator=(const BC_impl_list &x) {
      list = 0 ;
      copy_BC_list(x) ;
      return *this ;
    }
  } ;
  class register_BC_impl_list : public BC_impl_list {
  public:
    static BC_list_ent *global_list ;
    register_BC_impl_list() {}
    ~register_BC_impl_list() ;
    void clear() ;
    bool empty() ;
    void push_BC(register_BC_type *p) ;
    iterator begin() { return iterator(global_list) ; }
    iterator end() { return iterator(0) ; }
  } ;
  extern register_BC_impl_list register_BC_list ;

  template<class T> class register_BC : public register_BC_type {
  public:
    register_BC() { register_BC_list.push_BC(this) ; }
    virtual BC_implP get_BC() const { return BC_implP(new T) ; }
  } ;
  bool check_scalar_units(const options_list &o, std::string option,
                          std::string unit) ;
  bool check_vector_units(const options_list &ol,std::string vname,
                          std::string units) ;
}

#endif
