// System includes.
#include <algorithm>

// Standard library includes.
#include <vector>
#include <list>
using std::list ;
using std::make_pair ;
using std::pair ;
using std::vector ;

// Loci includes.
#include <Loci.h>
using Loci::send_clone_non ;
using Loci::send_global_clone_non ;
using Loci::fill_clone ;
using Loci::send_map ;
using Loci::send_global_map ; 
using Loci::GLOBAL_OR ;
using Loci::GLOBAL_MAX ;

// StreamUns includes.
#include "sciTypes.h"
#include "readGrid.h"
#include "name_var.h"

namespace streamUns {

  // Definition of global BC lists.
  register_BC_impl_list register_BC_list ;

  BC_impl_list::~BC_impl_list() {
    BC_list_ent *p,*v ;
    for(p=list;p!=0;p=v) {
      v = p->next ;
      delete p ;
    }
  }
  void BC_impl_list::clear() {
    BC_list_ent *p,*v ;
    for(p=list;p!=0;p=v) {
      v = p->next ;
      delete p ;
    }
    list = 0 ;
  }
  void BC_impl_list::push_BC(register_BC_type *p) {
    BC_list_ent *flp = new BC_list_ent(p,list) ;
    list = flp ;
  }
                                                                                                                                                      
  void BC_impl_list::copy_BC_list(const BC_impl_list& rl) {
    BC_list_ent *p, *v ;
    for(p = rl.list; p != 0; p=v) {
      push_BC(p->rr) ;
      v = p->next ;
    }
  }
  void BC_impl_list::copy_BC_list(const register_BC_impl_list& rl) {
    BC_list_ent *p, *v ;
    for(p = rl.global_list; p != 0; p=v) {
      push_BC(p->rr) ;
      v = p->next ;
    }
  }

  // Declaration of static variable global_list
  BC_impl_list::BC_list_ent *register_BC_impl_list::global_list = 0 ;

  register_BC_impl_list::~register_BC_impl_list() {
    BC_list_ent *p,*v ;
    for(p=global_list;p!=0;p=v) {
      v = p->next ;
      delete p ;
    }
  }
  void register_BC_impl_list::clear() {
    BC_list_ent *p,*v ;
    for(p=global_list;p!=0;p=v) {
      v = p->next ;
      delete p ;
    }
    global_list = 0 ;
  }
   bool register_BC_impl_list::empty() {
     return (global_list == 0) ;
   }
  void register_BC_impl_list::push_BC(register_BC_type *p) {
    BC_list_ent *flp = new BC_list_ent(p,global_list) ;
    global_list = flp ;
  }

  // Used for suboptions of the boundary condition? For example,
  // BC_5=supersonicInflow(p=10.0,T=5.0). The suboptions are p=10.0 and T=5.0.
  struct BCinfo {
    std::string name ;
    entitySet apply_set ;
    options_list bc_options ;
    BCinfo() {}
    BCinfo(const std::string &n,const entitySet &a, const options_list &o) :
      name(n),apply_set(a),bc_options(o) {}
  } ;

  bool check_scalar_units(const options_list &o,string option,string unit) {
    bool check = false ;
    if(o.getOptionValueType(option) == Loci::REAL)
      check = true ;
    if(o.getOptionValueType(option) == Loci::UNIT_VALUE) {
      Loci::UNIT_type Tu ;
      o.getOption(option,Tu) ;
      if(Tu.is_compatible(unit)) {
        check = true ;
      }
    }
    return check ;
  }

  bool check_vector_units(const options_list &ol,string vname,string units) {
    Loci::option_value_type ovt= ol.getOptionValueType(vname) ;
    if(ovt == Loci::REAL) {
      return true ;
    } else if(ol.getOptionValueType(vname) == Loci::UNIT_VALUE) {
      Loci::UNIT_type vu ;
      ol.getOption(vname,vu) ;
      if(!vu.is_compatible(units)) {
        return false ;
      }
      return true ;
    } else if(ovt == Loci::LIST) {
      Loci::options_list::arg_list value_list ;
      ol.getOption(vname,value_list) ;
      if(value_list.size() != 3) {
        return false ;
      }
      for(int i=0;i<3;++i)
        if(value_list[i].type_of() != Loci::REAL &&
           value_list[i].type_of() != Loci::UNIT_VALUE) {
          return false ;
        }
      for(int i=0;i<3;++i) {
        if(value_list[i].type_of() == Loci::UNIT_VALUE) {
          Loci::UNIT_type vu ;
          value_list[i].get_value(vu) ;
          if(!vu.is_compatible(units)) {
            return false ;
          }
        }
      }
      return true ;
    } else if(ovt == Loci::FUNCTION) {
      string name ;
      Loci::options_list::arg_list value_list ;
      ol.getOption(vname,name,value_list) ;
      if(name != "polar") {
        return false ;
      }
      if(value_list.size() != 3) {
        return false ;
      }
      for(int i=0;i<3;++i)
        if(value_list[i].type_of() != Loci::REAL &&
           value_list[i].type_of() != Loci::UNIT_VALUE) {
          return false ;
        }
      if(value_list[0].type_of() == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        value_list[0].get_value(vu) ;
        if(!vu.is_compatible(units)) {
          return false ;
        }
      }
      if(value_list[1].type_of() == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        value_list[1].get_value(vu) ;
        if(!vu.is_compatible("radians")) {
          return false ;
        }
      }
      if(value_list[2].type_of() == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        value_list[2].get_value(vu) ;
        if(!vu.is_compatible("radians")) {
          return false ;
        }
      }
      return true ;
    } else {
      return false ;
    }
  }

  class periodic_check : public BC_Check {
    string error_message ;
  public:
    std::string BoundaryConditions() { return "periodic" ; }
    std::string VariablesChecked(fact_db &facts) {
      return "rotate,translate,name,center,vector" ;
    }
    bool CheckOptions(const options_list& bc_options,fact_db &facts) {
      error_message = "" ;
      bool check = true ;
      if(!bc_options.optionExists("name")) {
        error_message = "all periodic boundaries must have a name " ;
        check = false ;
      }
      if(bc_options.optionExists("center")) {
        if(!check_vector_units(bc_options,"center","m")) {
          error_message += "center has incorrect units " ;
          check = false ;
        }
      }
      if(bc_options.optionExists("translate")) {
        if(!check_vector_units(bc_options,"translate","m")) {
          error_message += "translate has incorrect units " ;
          check = false ;
        }
      }
      if(bc_options.optionExists("vector")) {
        if(!check_vector_units(bc_options,"vector","")) {
          error_message += "vector has incorrect units " ;
          check = false ;
        }
      }
      if(bc_options.optionExists("rotate")) {
        if(!check_scalar_units(bc_options,"rotate","radians")) {
          error_message += "rotate has incorrect units " ;
          check = false ;
        }
      }
                                                                                                                                                      
      return check ;
    }
    std::ostream &ErrorMessage(std::ostream &s) { return s; }
  } ;
                                                                                                                                                      
  register_BC<periodic_check> register_BC_periodic_check ;

  struct bc_checker_info {
    Loci::variableSet bc_check ;
    Loci::variableSet bc_vars ;
    BC_implP checker ;
  } ;
  
  // Ensure the consistency of the boundary conditions. Code has been added
  // here to check for the simultaneous presence of fixedPressureOulet and
  // extrapolatedPressureOutlet boundaries, in which case, global conservation
  // is turned off for incompressible flows. Note that global conservation is
  // not used for compressible flows.
  bool CheckBoundaryConditions(fact_db &facts) {
    list<bc_checker_info> bclist ;
    for(BC_impl_list::iterator bci(register_BC_list.begin());bci!=
    register_BC_list.end();++bci) {
      bclist.push_back(bc_checker_info()) ;
      Loci::variableSet bcs ;
      if((*bci)->BoundaryConditions() == "*") {
        bcs = ~EMPTY ;
      } else {
        Loci::exprP bc_check =
          Loci::expression::create((*bci)->BoundaryConditions()) ;
        bcs = Loci::variableSet(bc_check) ;
      }
      if((*bci)->VariablesChecked(facts).size() == 0) {
        bclist.back().bc_vars = EMPTY ;
      } else {
        Loci::exprP bc_vars =
          Loci::expression::create((*bci)->VariablesChecked(facts)) ;
        bclist.back().bc_vars = Loci::variableSet(bc_vars) ;
      }
      bclist.back().bc_check = bcs ;
      bclist.back().checker = *bci ;
    }
    bool error = false ; bool doprint = (Loci::MPI_rank == 0) ;
    param<options_list> bc_info ;
    bc_info = facts.get_variable("boundary_conditions") ;
    options_list::option_namelist nl = bc_info->getOptionNameList() ;
    options_list::option_namelist::iterator li;
    bool foundFixedPressureOutlet=false,foundExtrapolatedPressureOutlet=false ;
    for(li=nl.begin();li!=nl.end();++li) {
      string bname = *li ;
      Loci::option_value_type vt =
        bc_info->getOptionValueType(bname);
      Loci::option_values ov = bc_info->getOption(bname) ;
      options_list::arg_list value_list ;
      string name ;
      bool func = false ;
      switch(vt) {
      case Loci::NAME :
        ov.get_value(name) ;
        bc_info->setOption(bname,name) ;
        break ;
      case Loci::FUNCTION:
        func = true ;
        ov.get_value(name) ;
        ov.get_value(value_list) ;
        bc_info->setOption(bname,name,value_list) ;
        break ;
      default:
        cerr << "setup_bc can not interpret value assigned to " << bname
             << " in boundary_conditions" << endl ;
        exit(-1) ;
      }
      if(name=="fixedPressureOutlet") foundFixedPressureOutlet=true ;
      if(name=="extrapolatedPressureOutlet") foundExtrapolatedPressureOutlet=
        true ;
      options_list ol ;
      ol.Input(value_list) ;
      Loci::variable bcv(name) ;
      list<bc_checker_info>::iterator bi ;
      options_list::option_namelist nlb = ol.getOptionNameList() ;
      Loci::variableSet bvars ;
      options_list::option_namelist::iterator lii;
      for(lii=nlb.begin();lii!=nlb.end();++lii)
        bvars += Loci::variable(*lii) ;

      bool found_bcmatch = false;
      bool found_empty_match = false ;
      Loci::variableSet unchecked_variables = bvars ;
      for(bi=bclist.begin();bi!=bclist.end();++bi) {
        if(bi->bc_check.inSet(bcv)) {
          if(bi->bc_check != ~EMPTY) {
            found_bcmatch = true ;
            if(bi->bc_vars == EMPTY) {
              found_empty_match = true ;
            }
          }
          if((bi->bc_vars & bvars) != EMPTY) {
            try {
              if(bi->checker->CheckOptions(ol,facts)) {
                unchecked_variables -= bi->bc_vars ;
              }
            } catch(const Loci::BasicException &err) {
              cerr << "ERROR: Boundary type " << name << " for boundary id "
                   << bname << ":" << endl ;
              err.Print(cerr) ;
              error = true ;
            }
          }
        }
      }

      if(!found_bcmatch) {
        if(doprint)
          cerr << "Boundary type '" << name << "' is unknown for boundary id "
               << bname << endl ;
        error = true ;
      } else if (bvars == EMPTY && !found_empty_match) {
        cerr << "Boundary type '" << name << "' requires argument(s) for "
          << "boundary id " << bname << endl ;
        error = true ;
      } else if(unchecked_variables != EMPTY) {
        if(doprint) {
          bool errorprinted = false ;
          for(bi=bclist.begin();bi!=bclist.end();++bi) {
            if(bi->bc_check.inSet(bcv)) {
              if((bi->bc_vars & unchecked_variables) != EMPTY) {
                cerr << "check failed for boundary condition " << name
                     << " on boundary id " << bname << endl ;
                bi->checker->ErrorMessage(cerr);
                errorprinted = true ;
              }
            }
          }
          if(!errorprinted)
            cerr << "variables " << unchecked_variables << " not compatible "
              << "with boundary condition " << name  << " for boundary id "
              << bname << endl  ;
        }
        error = true ;
      }
    }
    if(!foundFixedPressureOutlet && foundExtrapolatedPressureOutlet){
      param<bool> tmp ; facts.create_fact("globalConservation",tmp) ;
    }
    return error ;
  }

  void FindMindNoslip(fact_db &facts) {
    constraint geom_cells(facts.get_variable("geom_cells")) ;
                                                                                
    entitySet viscous_faces ;
                                                                                
    Loci::storeRepP sp = facts.get_variable("noslip_BC") ;
    if(sp != 0) {
      constraint tmp ;
      tmp = sp ;
      viscous_faces = *tmp ;
    }

    // If no viscous faces, then min_cell2noslip doesn't exist
    if(!GLOBAL_OR(viscous_faces.size()!=0))
      return ;

    multiMap face2node(facts.get_variable("face2node")) ;
    entitySet face_dom = face2node.domain() ;

    // Compute facecenter
    store<vect3d> pos_static(facts.get_variable("pos")) ;
    dstore<vect3d> pos ;

    // Copy static pos we currently have to dynamic pos
    FORALL(pos_static.domain(), pi) {
      pos[pi] = pos_static[pi] ;
    } ENDFORALL ;

    if(Loci::MPI_processes > 1) {
      // Compute nodes accessed by face2node
      entitySet total_dom = pos.domain() ;
      total_dom += Loci::MapRepP(face2node.Rep())->image(face_dom) ;
      entitySet out_dom = total_dom - pos.domain() ;
                                                                                
      // Fill in values of pos that are accessed but we don't own
      std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
      Loci::storeRepP posRep = pos.Rep() ;
      fill_clone(posRep, out_dom, init_ptn) ;
    }
                                                                                
    // Now compute facecenter for all faces we own
    dstore<vect3d> facecenter ;
    dstore<real> area ;
    FORALL(face_dom,fc) {
      vect3d nodesum(0.0,0.0,0.0) ;
      real lensum = 0 ;
      int fsz = face2node[fc].size() ;
                                                                                
      for(int id = 0; id < fsz; ++id) {
        const int nd1 = id ;
        const int nd2 = (id+1 == fsz)?0:(id+1) ;
        const vect3d p1 = pos[face2node[fc][nd1]] ;
        const vect3d p2 = pos[face2node[fc][nd2]] ;
                                                                                
        const vect3d edge_loc = 0.5*(p1 + p2) ;
        const vect3d edge_vec = p1-p2 ;
        real len = norm(edge_vec) ;
        nodesum += len*edge_loc ;
        lensum += len ;
      }
                                                                                
      vect3d center =  nodesum/lensum ;
      facecenter[fc] = center ;
                                                                                
      vect3d sum(0,0,0) ;
      for(int id = 0; id < fsz; ++id) {
        const int nd1 = id ;
        const int nd2 = (id+1 == fsz)?0:(id+1) ;
        const vect3d p1 = pos[face2node[fc][nd1]] ;
        const vect3d p2 = pos[face2node[fc][nd2]] ;
        sum += cross(p1-center,p2-center) ;
      }
      area[fc] = .5*norm(sum) ;
    } ENDFORALL ;
                                                                                
    // Release pos memory
    pos.setRep(dstore<vect3d>().Rep()) ;
                                                                                
    // Now compute cell center as average of face centers
    multiMap upper(facts.get_fact("upper")) ;
    multiMap lower(facts.get_fact("lower")) ;
    multiMap boundary_map(facts.get_fact("boundary_map")) ;
                                                                                
    if(Loci::MPI_processes > 1) {
      // Fill in ghost facecenter
      entitySet face_out = face_dom ;
      face_out += Loci::MapRepP(lower.Rep())->image(*geom_cells) ;
      face_out += Loci::MapRepP(boundary_map.Rep())->image(*geom_cells) ;
      face_out += Loci::MapRepP(upper.Rep())->image(*geom_cells) ;
      // Fill clone for facecenter
      std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
      Loci::storeRepP facecenterRep = facecenter.Rep() ;
      fill_clone(facecenterRep, face_out, init_ptn) ;
      // fill clone for face area ;
      Loci::storeRepP areaRep = area.Rep() ;
      fill_clone(areaRep, face_out, init_ptn) ;
    }
                                                                                
    store<vect3d> cellcenter ;
    cellcenter.allocate(*geom_cells) ;
    FORALL(*geom_cells,cc) {
      vect3d nodesum(0.0,0.0,0.0) ;
      real areasum = 0 ;
      for(Entity *id = upper.begin(cc); id != upper.end(cc); ++id) {
        double w = area[*id] ;
        nodesum += facecenter[*id]*w ;
        areasum += w ;
      }
      for(Entity *id = lower.begin(cc); id != lower.end(cc); ++id) {
        double w = area[*id] ;
        nodesum += facecenter[*id]*w ;
        areasum += w;
      }
      for(Entity *id=boundary_map.begin(cc); id!=boundary_map.end(cc); ++id) {
        double w = area[*id] ;
        nodesum += facecenter[*id]*w ;
        areasum += w ;
      }
      cellcenter[cc] = nodesum/areasum ;
    } ENDFORALL ;
                                                                                
    // Find min_noslip here
    double t1 = MPI_Wtime()  ;
                                                                                
    vector<Loci::kdTree::coord3d> viscous_pts(viscous_faces.size()) ;
    vector<int> viscous_ids(viscous_faces.size()) ;
                                                                                
    int cnt = 0 ;
    FORALL(viscous_faces,fc) {
      viscous_pts[cnt][0] = facecenter[fc].x ;
      viscous_pts[cnt][1] = facecenter[fc].y ;
      viscous_pts[cnt][2] = facecenter[fc].z ;
      viscous_ids[cnt] = fc ;
      cnt++ ;
    } ENDFORALL ;
                                                                                
    vector<Loci::kdTree::coord3d> cell_pts((*geom_cells).size()) ;
    vector<int> closest((*geom_cells).size(),-1) ;
    cnt = 0 ;
    FORALL(*geom_cells,cc) {
      cell_pts[cnt][0] = cellcenter[cc].x ;
      cell_pts[cnt][1] = cellcenter[cc].y ;
      cell_pts[cnt][2] = cellcenter[cc].z ;
      cnt++ ;
    } ENDFORALL ;
                                                                                
    Loci::parallelNearestNeighbors(viscous_pts,viscous_ids,cell_pts,closest,
                                   MPI_COMM_WORLD) ;
                                                                                
    Map min_cell2noslip ;
    min_cell2noslip.allocate(*geom_cells) ;
    cnt = 0 ;
    FORALL(*geom_cells,cc) {
      min_cell2noslip[cc] = closest[cnt] ;
      cnt++ ;
    } ENDFORALL ;
                                                                                
    double t2 = MPI_Wtime()  ;
    Loci::debugout << "Time spent in new min_2noslip is " << t2-t1 << endl ;
                                                                                
    facts.create_fact("min_cell2noslip",min_cell2noslip) ;
  }

}
