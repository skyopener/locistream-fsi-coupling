#ifndef INTERPOLATE_FILE_H
#define INTERPOLATE_FILE_H
#include <Loci.h>
using Loci::Area ;
#include "sciTypes.h"
#include <string>
#include <vector>

// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;

namespace streamUns {
  void read_puT_file(std::string filename,
                     const EOS &eos,
                     store<vect3d> &loc,
                     store<real> &p,
                     store<real> &T,
                     store<vect3d> &vel,
                     store<real> &k,
                     store<real> &mu_t,
                     storeVec<real> &mixture,
                     std::vector<int> &mixids) ;
  void read_scalar_file(std::string filename,
                        store<vect3d> &loc,
                        store<real> &val) ;
  void read_puT_file_serial(std::string filename,
                     const EOS &eos,
                     store<vect3d> &loc,
                     store<real> &p,
                     store<real> &T,
                     store<vect3d> &vel,
                     store<real> &k,
                     store<real> &mu_t,
                     storeVec<real> &mixture,
                     std::vector<int> &mixids) ;
  void read_scalar_file_serial(std::string filename,
                        store<vect3d> &loc,
                        store<real> &val) ;
  void broadcast_storeRep(Loci::storeRepP rep) ;

  void stencil_weights(std::vector<double> &w,
                       std::vector<int> &neighbors,
                       const store<vect3d> &loc,
                       vect3d fcenter) ;
  std::vector<int> get_stencil(const Loci::kdTree::kd_tree &kd,vect3d pnt,
                               double delta = 1e100) ;


}
#endif


