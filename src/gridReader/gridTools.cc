// Standard library includes.
#include <algorithm>
#include <iostream>
#include <vector>
using std::endl ;
using std::sort ;
using std::vector ;

// StreamUns includes.
#include "sciTypes.h"
#include "readGrid.h"
#include "name_var.h"

// Loci includes.
using Loci::all_collect_entitySet ;
using Loci::MPI_rank ;
using Loci::MPI_processes ;

namespace streamUns {

  void CreateCellStencil(fact_db &facts) {
    using std::vector ;
    store<vect3d> pos ;
    pos = facts.get_variable("pos") ;
    multiMap face2node ;
    face2node = facts.get_variable("face2node") ;
    constraint faces ;
    constraint interior_faces ;
    faces = facts.get_variable("faces") ;
    interior_faces = facts.get_variable("interior_faces") ;
    entitySet global_nodes = all_collect_entitySet(pos.domain(),facts) ;
    entitySet global_faces = all_collect_entitySet(*faces,facts) ;
    dmultiMap node2face,dface2node ;
    distributed_inverseMap(node2face,face2node,global_nodes,global_faces,
      facts) ;
    FORALL(face2node.domain(),fc) {
      for(Entity *p = face2node.begin(fc);p!=face2node.end(fc);++p)
        dface2node[fc].push_back(*p) ;
    } ENDFORALL ;
    dmultiMap cell2face ;
    dmultiMap lower,upper,boundary_map ;
    lower= facts.get_variable("lower") ;
    upper= facts.get_variable("upper") ;
    boundary_map = facts.get_variable("boundary_map") ;

    entitySet cells = lower.domain() ;
    FORALL(cells,cc) {
      for(size_t i=0;i<lower[cc].size();++i)
        cell2face[cc].push_back(lower[cc][i]) ;
      for(size_t i=0;i<upper[cc].size();++i)
        cell2face[cc].push_back(upper[cc][i]) ;
      for(size_t i=0;i<boundary_map[cc].size();++i)
        cell2face[cc].push_back(boundary_map[cc][i]) ;
    } ENDFORALL ;
    Map cl,cr ;
    cl = facts.get_variable("cl") ;
    cr = facts.get_variable("cr") ;
    if(Loci::MPI_processes > 1) { // Expand maps
      vector<Loci::storeRepP> lp1,lp2 ;
      lp1.push_back(cell2face.Rep()) ;
      lp1.push_back(dface2node.Rep()) ;
      lp1.push_back(node2face.Rep()) ;
      lp2.push_back(cr.Rep()) ;
      lp2.push_back(cl.Rep()) ;
      std::vector<entitySet> ptn = facts.get_init_ptn() ;
      entitySet locdom = lp1[0]->domain() ;
      for(size_t j=0;j<lp1.size();++j) {
        Loci::storeRepP p = lp1[j] ;
        entitySet tmp_dom = p->domain() ;
        Loci::MapRepP mp =  Loci::MapRepP(p->getRep()) ;
        entitySet glob_dom = all_collect_entitySet(tmp_dom,facts) ;
        //          glob_dom &= dom ;
        entitySet tmp_out = (glob_dom & locdom) - tmp_dom ;
        Loci::storeRepP sp = mp->expand(tmp_out, ptn) ;
        entitySet image =  Loci::MapRepP(sp)->image((sp->domain()) & locdom) ;
        lp1[j] = sp ;
        locdom = image ;
      }
      for(size_t j=0;j<lp2.size();++j) {
        Loci::storeRepP p = lp2[j] ;
        entitySet tmp_dom = p->domain() ;
        Loci::MapRepP mp =  Loci::MapRepP(p->getRep()) ;
        entitySet glob_dom = all_collect_entitySet(tmp_dom,facts) ;
        //          glob_dom &= dom ;
        entitySet tmp_out = (glob_dom & locdom) - tmp_dom ;
        Loci::storeRepP sp = mp->expand(tmp_out, ptn) ;
        lp2[j] = sp ;
      }
      cell2face = lp1[0] ;
      dface2node = lp1[1] ;
      node2face = lp1[2] ;
      cr = lp2[0] ;
      cl = lp2[1] ;
    }
    entitySet global_cells = all_collect_entitySet(cells,facts) ;
    dmultiMap cellStencil ;
    FORALL(cells,cc) {
      entitySet cell_list ;
      for(size_t i=0;i<cell2face[cc].size();++i) {
        const Entity fc = cell2face[cc][i] ;
        for(size_t k = 0;k<dface2node[fc].size();++k) {
          int nd = dface2node[fc][k] ;
          for(size_t j=0;j<node2face[nd].size();++j) {
            Entity fc2 = node2face[nd][j] ;
            if(global_cells.inSet(cl[fc2]))
              cell_list += cl[fc2] ;
            if(global_cells.inSet(cr[fc2]))
              cell_list += cr[fc2] ;
          }
        }
      }
      cell_list -= cc ;
      std::vector<int,Loci::malloc_alloc<int> > cvec(cell_list.size()) ;
      int kk=0 ;
      for(entitySet::const_iterator ei=cell_list.begin();ei!= cell_list.end();
      ++ei) {
        cvec[kk++] = *ei ;
      }
      cellStencil[cc].swap(cvec) ;
    } ENDFORALL ;
    facts.create_fact("cellStencil",cellStencil) ;
  }
}
  
