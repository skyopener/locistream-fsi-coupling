#line 1 "nodeValues.loci"
// Standard library includes.
#include <algorithm>
using std::find ;
#include <map>
using std::map ;
#include <vector>
using std::vector ;

// Loci includes.
#include <Loci.h>

// StreamUns includes.
#include "sciTypes.h"
#line 1 "FVM.lh"
//#############################################################################
//#
//# Copyright 2008, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################

// $type pos store<Loci::vector3d<Loci::real_t> > 
// $type cl Map
// $type cr Map
// $type ci Map
// $type ref Map
// $type pmap Map
// $type face2node multiMap

// $type upper multiMap
// $type lower multiMap
// $type boundary_map multiMap

// $type cellcenter store<Loci::vector3d<Loci::real_t> > 
// $type facecenter store<Loci::vector3d<Loci::real_t> > 
// $type area store<Loci::Area> 
// $type vol store<Loci::real_t> 
// $type grid_vol param<Loci::real_t> 

// $type mn store<Loci::vector3d<Loci::real_t> > 
// $type ln store<Loci::vector3d<Loci::real_t> > 

// $type grads(X0) store<Loci::vector3d<Loci::real_t> > 
// $type gradv(X0) storeVec<Loci::vector3d<Loci::real_t> > 
// $type gradv3d(X0) store<Loci::tensor3d<Loci::real_t> > 

// $type grads_f(X0) store<Loci::vector3d<Loci::real_t> > 
// $type gradv_f(X0) storeVec<Loci::vector3d<Loci::real_t> > 
// $type gradv3d_f(X0) store<Loci::tensor3d<Loci::real_t> > 

// $type limiters(X0) store<Loci::real_t> 
// $type limiterv(X0) storeVec<Loci::real_t> 
// $type limiterv3d(X0) store<Loci::vector3d<Loci::real_t> > 

// $type lefts(X0) store<Loci::real_t> 
// $type rights(X0) store<Loci::real_t> 
// $type leftsP(X0,X1) store<Loci::real_t> 
// $type rightsP(X0,X1) store<Loci::real_t> 
// $type leftvM(X0) storeVec<Loci::real_t> 
// $type rightvM(X0) storeVec<Loci::real_t> 
// $type leftv3d(X0) store<Loci::vector3d<Loci::real_t> > 
// $type rightv3d(X0) store<Loci::vector3d<Loci::real_t> > 

// $type cell2node(X0) store<float> 
// $type cell2node_v(X0) storeVec<float> 
// $type cell2node_v3d(X0) store<Loci::vector3d<float> > 
// $type cell2nodeMax(X0) store<float> 
// $type cell2nodeMin(X0) store<float> 
// $type cell2nodeMaxMag(X0) store<float> 
// $type cell2nodeMaxv3d(X0) store<vector3d<float> > 

// $type BC_options store<Loci::options_list> 

// $type integrateSurface(X0) store<Loci::real_t> 
// $type integrateFlux(X0) store<Loci::real_t> 

// $type petscScalarSolve(X0) store<Loci::real_t> 
// $type petscBlockedSolve(X0) storeVec<Loci::real_t> 
// $type petscBlockedSSolve(X0) storeVec<Loci::real_t> 

// $type L2Norm(X0) param<Loci::real_t> 
// $type L1Norm(X0) param<Loci::real_t> 
// $type LinfNorm(X0) param<Loci::real_t> 
#line 14 "nodeValues.loci"


namespace streamUns {

  // Copy of Ed's cell2node rules which use double instead of float.

  // $type dnodalw_sum store<real> 

  namespace {class file_nodeValues000_1279563563m128 : public Loci::unit_rule {
#line 22 "nodeValues.loci"
    Loci::store<real>  L_dnodalw_sum_ ; 
#line 22 "nodeValues.loci"
public:
#line 22 "nodeValues.loci"
    file_nodeValues000_1279563563m128() {
#line 22 "nodeValues.loci"
       name_store("dnodalw_sum",L_dnodalw_sum_) ;
#line 22 "nodeValues.loci"
       output("dnodalw_sum") ;
#line 22 "nodeValues.loci"
       constraint("pos") ;
#line 22 "nodeValues.loci"
    }
#line 22 "nodeValues.loci"
    void calculate(Entity _e_) { 
#line 23 "nodeValues.loci"
    L_dnodalw_sum_[_e_]= 0 ;
  }    void compute(const Loci::sequence &seq) { 
#line 24 "nodeValues.loci"
      do_loop(seq,this) ;
#line 24 "nodeValues.loci"
    }
#line 24 "nodeValues.loci"
} ;
#line 24 "nodeValues.loci"
Loci::register_rule<file_nodeValues000_1279563563m128> register_file_nodeValues000_1279563563m128 ;
#line 24 "nodeValues.loci"
}
#line 24 "nodeValues.loci"


  namespace {class file_nodeValues001_1279563563m129 : public Loci::apply_rule< store<real> ,Loci::Summation<real> >  {
#line 28 "nodeValues.loci"
    Loci::const_store<Loci::vector3d<Loci::real_t> >  L_pos_ ; 
#line 28 "nodeValues.loci"
    Loci::const_multiMap L_face2node_ ; 
#line 28 "nodeValues.loci"
    Loci::const_multiMap L_upper_ ; 
#line 28 "nodeValues.loci"
    Loci::const_multiMap L_lower_ ; 
#line 28 "nodeValues.loci"
    Loci::const_multiMap L_boundary_map_ ; 
#line 28 "nodeValues.loci"
    Loci::const_store<Loci::vector3d<Loci::real_t> >  L_cellcenter_ ; 
#line 28 "nodeValues.loci"
    Loci::store<real>  L_dnodalw_sum_ ; 
#line 28 "nodeValues.loci"
public:
#line 28 "nodeValues.loci"
    file_nodeValues001_1279563563m129() {
#line 28 "nodeValues.loci"
       name_store("pos",L_pos_) ;
#line 28 "nodeValues.loci"
       name_store("face2node",L_face2node_) ;
#line 28 "nodeValues.loci"
       name_store("upper",L_upper_) ;
#line 28 "nodeValues.loci"
       name_store("lower",L_lower_) ;
#line 28 "nodeValues.loci"
       name_store("boundary_map",L_boundary_map_) ;
#line 28 "nodeValues.loci"
       name_store("cellcenter",L_cellcenter_) ;
#line 28 "nodeValues.loci"
       name_store("dnodalw_sum",L_dnodalw_sum_) ;
#line 28 "nodeValues.loci"
       input("(upper,lower,  boundary_map)->face2node->dnodalw_sum,cellcenter,(upper,lower,boundary_map)->  face2node->pos") ;
#line 28 "nodeValues.loci"
       output("(boundary_map,lower,upper)->face2node->dnodalw_sum") ;
#line 28 "nodeValues.loci"
       constraint("geom_cells") ;
#line 28 "nodeValues.loci"
    }
#line 28 "nodeValues.loci"
    void calculate(Entity _e_) { 
#line 29 "nodeValues.loci"
    int sztot = 0 ;

    for (const Entity *fi =L_upper_[_e_].begin ();fi !=L_upper_[_e_].end ();++fi )
      sztot += L_face2node_[*fi].size () ;
    for (const Entity *fi =L_lower_[_e_].begin ();fi !=L_lower_[_e_].end ();++fi )
      sztot += L_face2node_[*fi].size () ;
    for (const Entity *fi =L_boundary_map_[_e_].begin ();fi !=L_boundary_map_[_e_].end ();++fi )
      sztot += L_face2node_[*fi].size () ;

    vector <Entity > node_list (sztot ) ;
    int cnt = 0 ;
    for (const Entity *fi =L_upper_[_e_].begin ();fi !=L_upper_[_e_].end ();++fi )
      for (const Entity *ni =L_face2node_[*fi].begin ();ni !=L_face2node_[*fi].end ();++ni )
        node_list [cnt ++] = *ni ;
    for (const Entity *fi =L_lower_[_e_].begin ();fi !=L_lower_[_e_].end ();++fi )
      for (const Entity *ni =L_face2node_[*fi].begin ();ni !=L_face2node_[*fi].end ();++ni )
        node_list [cnt ++] = *ni ;
    for (const Entity *fi =L_boundary_map_[_e_].begin ();fi !=L_boundary_map_[_e_].end ();++fi )
      for (const Entity *ni =L_face2node_[*fi].begin ();ni !=L_face2node_[*fi].end ();++ni )
        node_list [cnt ++] = *ni ;

    sort (node_list .begin (),node_list .end ()) ;
    vector <Entity >::iterator ne = unique (node_list .begin (),node_list .end ()) ;
    vector <Entity >::iterator ns = node_list .begin () ;
    for (vector <Entity >::iterator vi = ns ;vi !=ne ;++vi ) {
      const real weight = 1./norm (L_pos_[*vi]-L_cellcenter_[_e_]) ;
      join (L_dnodalw_sum_[*vi],weight ) ;
    }
  }    void compute(const Loci::sequence &seq) { 
#line 57 "nodeValues.loci"
      do_loop(seq,this) ;
#line 57 "nodeValues.loci"
    }
#line 57 "nodeValues.loci"
} ;
#line 57 "nodeValues.loci"
Loci::register_rule<file_nodeValues001_1279563563m129> register_file_nodeValues001_1279563563m129 ;
#line 57 "nodeValues.loci"
}
#line 57 "nodeValues.loci"


  class dc2n_scalar_unit : public unit_rule {
    store<real> c2n_scalar_sum ;
  public:
    dc2n_scalar_unit() {
      name_store("dc2n_scalar_sum(X)",c2n_scalar_sum) ;
      constraint("pos") ;
      output("dc2n_scalar_sum(X)") ;
    }
    void calculate(Entity nd) {
      c2n_scalar_sum[nd] = 0;
    }
    virtual void compute(const sequence &seq) {
      do_loop(seq,this) ;
    }
  } ;
  register_rule<dc2n_scalar_unit> register_dc2n_scalar_unit ;

  class dc2n_scalar : public apply_rule<store<real>,
  Loci::Summation<real> > {
    const_multiMap upper, lower, boundary_map ;
    const_multiMap face2node ;
    const_store<vect3d> pos ;
    const_store<vect3d> cellcenter ;
    const_store<real> X ;

    store<real> c2n_scalar_sum ;
    vector<int> node_list ;
  public:
    dc2n_scalar() ;
    void calculate(Entity cc) ;
    virtual void compute(const sequence &seq) ;
  } ;

  dc2n_scalar::dc2n_scalar() {
    name_store("dc2n_scalar_sum(X)",c2n_scalar_sum) ;
    name_store("face2node",face2node) ;
    name_store("upper",upper) ;
    name_store("lower",lower) ;
    name_store("boundary_map",boundary_map) ;
    name_store("pos",pos) ;
    name_store("cellcenter",cellcenter) ;
    name_store("X",X) ;

    constraint("geom_cells") ;
    input("X") ;
    input("cellcenter") ;
    input("(upper,lower,boundary_map)->face2node->pos") ;
    input("(upper,lower,boundary_map)->face2node->dc2n_scalar_sum(X)") ;
    output("(upper,lower,boundary_map)->face2node->dc2n_scalar_sum(X)") ;
  }

  void dc2n_scalar::calculate(Entity cc) {

    node_list.clear() ;
    for(const Entity *fi=upper.begin(cc);fi!=upper.end(cc);++fi)
      for(const Entity *ni=face2node.begin(*fi);ni!=face2node.end(*fi);++ni)
        node_list.push_back(*ni) ;
    for(const Entity *fi=lower.begin(cc);fi!=lower.end(cc);++fi)
      for(const Entity *ni=face2node.begin(*fi);ni!=face2node.end(*fi);++ni)
        node_list.push_back(*ni) ;
    for(const Entity *fi=boundary_map.begin(cc);fi!=boundary_map.end(cc);++fi)
      for(const Entity *ni=face2node.begin(*fi);ni!=face2node.end(*fi);++ni)
        node_list.push_back(*ni) ;
    sort(node_list.begin(),node_list.end()) ;
    vector<int>::iterator ns = node_list.begin() ;
    vector<int>::iterator ne = unique(node_list.begin(),node_list.end()) ;



    for(vector<int>::iterator vi = ns;vi!=ne;++vi) {
      int nd = *vi ;
      const real weight = 1./norm(pos[nd]-cellcenter[cc]) ;
      join(c2n_scalar_sum[nd],real(X[cc]*weight)) ;
    }

  }

  void dc2n_scalar::compute(const sequence &seq) {
    do_loop(seq,this) ;
  }

  register_rule<dc2n_scalar> register_dc2n_scalar ;

  class dcompute_scalar_nodal : public pointwise_rule {
    const_store<real> nodal_sum ;
    const_store<real> nodalw_sum ;
    store<real> nodal ;
  public:
    dcompute_scalar_nodal() {
      name_store("dc2n_scalar_sum(X)",nodal_sum) ;
      name_store("dnodalw_sum",nodalw_sum) ;
      name_store("dcell2node(X)",nodal) ;
      input("dc2n_scalar_sum(X),dnodalw_sum") ;
      constraint("pos") ;
      output("dcell2node(X)") ;
    }
    void calculate(Entity nd) {
      real rsum = 1./(real(nodalw_sum[nd])+1e-20) ;
      nodal[nd] = nodal_sum[nd]*rsum ;
    } ;
    virtual void compute(const sequence &seq) {
      do_loop(seq,this) ;
    }
  } ;

  register_rule<dcompute_scalar_nodal> register_dcompute_scalar_nodal ;

  class dc2n_v3d_unit : public unit_rule {
    store<vector3d<real> > c2n_v3d_sum ;
  public:
    dc2n_v3d_unit() {
      name_store("dc2n_v3d_sum(X)",c2n_v3d_sum) ;
      constraint("pos") ;
      output("dc2n_v3d_sum(X)") ;
    }
    void calculate(Entity nd) {
      c2n_v3d_sum[nd] = vector3d<real>(0.,0.,0.);
    }
    virtual void compute(const sequence &seq) {
      do_loop(seq,this) ;
    }
  } ;
  register_rule<dc2n_v3d_unit> register_dc2n_v3d_unit ;

  class dc2n_v3d : public apply_rule<store<vector3d<real> >,
                  Loci::Summation<vector3d<real> > > {
    const_multiMap upper, lower,boundary_map ;
    const_multiMap face2node ;
    const_store<vect3d> pos ;
    const_store<vect3d> cellcenter ;
    const_store<vect3d> X ;

    store<vector3d<real> > c2n_v3d_sum ;
    vector<int> node_list ;
  public:
    dc2n_v3d() ;
    void calculate(Entity cc) ;
    virtual void compute(const sequence &seq) ;
  } ;

  dc2n_v3d::dc2n_v3d() {
    name_store("dc2n_v3d_sum(X)",c2n_v3d_sum) ;
    name_store("face2node",face2node) ;
    name_store("upper",upper) ;
    name_store("lower",lower) ;
    name_store("boundary_map",boundary_map) ;
    name_store("pos",pos) ;
    name_store("cellcenter",cellcenter) ;
    name_store("X",X) ;

    constraint("geom_cells") ;
    input("X") ;
    input("cellcenter") ;
    input("(upper,lower,boundary_map)->face2node->pos") ;
    input("(upper,lower,boundary_map)->face2node->dc2n_v3d_sum(X)") ;
    output("(upper,lower,boundary_map)->face2node->dc2n_v3d_sum(X)") ;
  }

  void dc2n_v3d::calculate(Entity cc) {

    node_list.clear() ;
    for(const Entity *fi=upper.begin(cc);fi!=upper.end(cc);++fi)
      for(const Entity *ni=face2node.begin(*fi);ni!=face2node.end(*fi);++ni)
        node_list.push_back(*ni) ;
    for(const Entity *fi=lower.begin(cc);fi!=lower.end(cc);++fi)
      for(const Entity *ni=face2node.begin(*fi);ni!=face2node.end(*fi);++ni)
        node_list.push_back(*ni) ;
    for(const Entity *fi=boundary_map.begin(cc);fi!=boundary_map.end(cc);++fi)
      for(const Entity *ni=face2node.begin(*fi);ni!=face2node.end(*fi);++ni)
        node_list.push_back(*ni) ;
    sort(node_list.begin(),node_list.end()) ;
    vector<int>::iterator ns = node_list.begin() ;
    vector<int>::iterator ne = unique(node_list.begin(),node_list.end()) ;



    for(vector<int>::iterator vi = ns;vi!=ne;++vi) {
      int nd = *vi ;
      const real weight = 1./norm(pos[nd]-cellcenter[cc]) ;
      vector3d<real> v(X[cc].x*weight,X[cc].y*weight,X[cc].z*weight) ;
      join(c2n_v3d_sum[nd],v) ;
    }

  }

  void dc2n_v3d::compute(const sequence &seq) {
    do_loop(seq,this) ;
  }

  register_rule<dc2n_v3d> register_dc2n_v3d ;

  class dcompute_v3d_nodal : public pointwise_rule {
    const_store<vector3d<real> > nodal_sum ;
    const_store<real> nodalw_sum ;
    store<vector3d<real> > nodal ;
  public:
    dcompute_v3d_nodal() {
      name_store("dc2n_v3d_sum(X)",nodal_sum) ;
      name_store("dnodalw_sum",nodalw_sum) ;
      name_store("dcell2node_v3d(X)",nodal) ;
      input("dc2n_v3d_sum(X),dnodalw_sum") ;
      constraint("pos") ;
      output("dcell2node_v3d(X)") ;
    }
    void calculate(Entity nd) {
      real rsum = 1./(real(nodalw_sum[nd])+1e-20) ;
      nodal[nd] = nodal_sum[nd]*rsum ;
    } ;
    virtual void compute(const sequence &seq) {
      do_loop(seq,this) ;
    }
  } ;

  register_rule<dcompute_v3d_nodal> register_dcompute_v3d_nodal ;

}
