#include "interpolateFile.h"

using std::string ;
using std::vector ;
namespace Loci {
  void ORBPartition(const vector<vector3d<float> > &pnts,
                    vector<int> &procid,
                    MPI_Comm comm) ;
}
namespace streamUns {

  // Compute stencil using nearest point in 8 octants
  vector<int> get_stencil(const Loci::kdTree::kd_tree &kd,vect3d pnt,
                          double delta) {
    vector<int> neighbors ;
    Loci::kdTree::coord3d ccenter ;
    ccenter[0] = pnt.x ;
    ccenter[1] = pnt.y ;
    ccenter[2] = pnt.z ;

    double rmin = std::numeric_limits<float>::max() ;

    neighbors.push_back(kd.find_closest(ccenter,rmin)) ;

    if(rmin <= 1e-30)
      return neighbors ;

    // First gather postive z quadrants
    Loci::kdTree::kd_tree::bounds box ;

    for(int i=0;i<3;++i) {
      box.minc[i] = ccenter[i] ;
      box.maxc[i] = ccenter[i]+delta ;
    }
    int id = kd.find_closest_box(ccenter,box) ;
    if(id >=0)
      neighbors.push_back(id) ;

    box.maxc[0] = ccenter[0] ;
    box.minc[0] = ccenter[0]-delta ;
    id = kd.find_closest_box(ccenter,box) ;
    if(id >=0)
      neighbors.push_back(id) ;

    box.maxc[1] = ccenter[1] ;
    box.minc[1] = ccenter[1]-delta ;
    id = kd.find_closest_box(ccenter,box) ;
    if(id >=0)
      neighbors.push_back(id) ;

    box.minc[0] = ccenter[0] ;
    box.maxc[0] = ccenter[0]+delta ;
    id = kd.find_closest_box(ccenter,box) ;
    if(id >=0)
      neighbors.push_back(id) ;

    // Now gather negative z quadrants
    for(int i=0;i<3;++i) {
      box.minc[i] = ccenter[i] ;
      box.maxc[i] = ccenter[i]+delta ;
    }
    box.maxc[2] = ccenter[2] ;
    box.minc[2] = ccenter[2]-delta ;
    id = kd.find_closest_box(ccenter,box) ;
    if(id >=0)
      neighbors.push_back(id) ;
    box.maxc[0] = ccenter[0] ;
    box.minc[0] = ccenter[0]-delta ;
    id = kd.find_closest_box(ccenter,box) ;
    if(id >=0)
      neighbors.push_back(id) ;
    box.maxc[1] = ccenter[1] ;
    box.minc[1] = ccenter[1]-delta ;
    id = kd.find_closest_box(ccenter,box) ;
    if(id >=0)
      neighbors.push_back(id) ;
    box.minc[0] = ccenter[0] ;
    box.maxc[0] = ccenter[0]+delta ;
    id = kd.find_closest_box(ccenter,box) ;
    if(id >=0)
      neighbors.push_back(id) ;


    std::sort(neighbors.begin(),neighbors.end()) ;
    vector<int> tmp ;
    tmp.push_back(neighbors[0]) ;
    int cnt = 0 ;
    for(size_t i=1;i<neighbors.size();++i)
      if(tmp[cnt] != neighbors[i]) {
        tmp.push_back(neighbors[i]) ;
        cnt++ ;
      }

    return tmp ;
  }

  // Compute stencil weights from a set of candidate neighbors
  // The code outputs an updated neighbors list and stencil weights
  // This routine uses the following algorithm to determine the stencil:
  //  1) find the closest point, it is in the stencil.
  //  2) find the closest edge (closest orthogonal projection to the
  //     edge formed by connecting the closest point and some other point
  //     in the stencil.  If no edge contains the point projection, then
  //     the interpolation stencil contains only the closest point
  //  3) Find the closest point that projects onto the face formed from
  //     the closest edge and all remaining points.  If no face contains
  //     the projected point, then the stencil only contains the points of
  //     the closest edge.
  //  4) Find the smallest volume tetrahedra that contains the interpolant
  //     and is formed from the closest face and remaining points in the
  //     stencil.  If no tetrahedra satisfy this condition, then interpolate
  //     using projection on closest face.  Otherwise compute weights based
  //     on barycentric coordinates of the smallest tetrahedra.

  void stencil_weights(std::vector<double> &w,
                       std::vector<int> &neighbors,
                       const store<vect3d> &loc,
                       vect3d ipnt) {

    // Find closest point, set it to n0
    int n0 = 0 ;
    int sz = neighbors.size() ;
    double dist = 1e33 ;
    for(int i=0;i<sz;++i) {
      const vect3d dv = loc[neighbors[i]]-ipnt;

      const double d2 = dot(dv,dv) ;
      if(d2 < dist) {
        n0 = i ;
        dist = d2 ;
      }
    }

    // Vertex 0
    vect3d v0 = loc[neighbors[n0]] ;

    // Vertex of interpolatant point relative to v0
    vect3d vc = ipnt-v0 ;

    // Find closest edge
    int n1 = -1 ;
    dist = 1e33 ;
    for(int i=0;i<sz;++i)
      if(i != n0) { // If not closet point, then form edge
        // Compute normalized edge vector
        vect3d ev = loc[neighbors[i]] - v0 ;
        double evlen = norm(ev) ;
        ev *= 1./evlen ;
        // project vc onto ev making sure it stays between endpoints
        double pnte = dot(vc,ev) ;
        if(pnte >= 0 && pnte <= evlen) { // If projection on edge
          // Compute distance to edge
          vect3d v2e = vc- pnte*ev ;
          double d = dot(v2e,v2e) ;
          if(dist > d) { // If it is closer, use this one
            n1 = i ;
            dist = d ;
          }
        }
      }

    if(n1 == -1) {
      // No close edge, so the only choice
      //is to use closest point as stencil
      vector<int> N(1) ;
      vector<double> W(1) ;
      N[0] = neighbors[n0] ;
      W[0] = 1.0 ;
      w.swap(W) ;
      neighbors.swap(N) ;
      return ;
    }

    // Now find closest triangle
    dist = 1e33 ;
    int n2 = -1 ;
    vect3d v1 = loc[neighbors[n1]]-v0 ;
    double scale = dot(v1,v1) ;
    for(int j=0;j<sz;++j)  // Loop over stencil points
      if(j!= n0 && j != n1) {
        // Form triangle (candidate triangle formed from v2)
        const vect3d v2 = loc[neighbors[j]]-v0 ;
        vect3d n = cross(v1,v2) ; // Compute face normal
        double a = norm(n) ;
        if(a < 1e-10*scale) // degenerate triangle, ignore
          continue ;
        n *= 1./a ;
        double d = dot(vc,n) ; // projected distance to triangle
        double ad = abs(d) ;   // magnitude of distance
        if(ad < dist) {
          // Compute projected point on triangle plane
          vect3d pnt = (vc-d*n) ;
          // Now check to see if pnt is in triangle
          // Compute barycentric coordinates and make sure all are positive
          double c1 = dot(n,cross(v1,pnt)) ;
          double c2 = dot(n,cross(pnt,v2)) ;
          double c3 = dot(n,cross(v1-pnt,v2-pnt)) ;
          if(c1>=0 && c2>=0 && c3>=0) {
            dist = ad ; // inside triangle && closest so update
            n2 = j ;
          }
        }
      }

    if(n2 == -1) {
      // No closest triangle, so we must project onto a line segment
      vector<int> N(2) ;
      vector<double> W(2) ;
      N[0] = neighbors[n0] ;
      N[1] = neighbors[n1] ;
      vect3d ev = v1 ;
      double evlen = norm(ev) ;
      ev *= 1./evlen ;
      // project vc onto ev making sure it stays between endpoints
      double pnte = dot(vc,ev) ;
      W[1] = pnte/evlen ;
      W[0] = 1.-W[1] ;
      w.swap(W) ;
      neighbors.swap(N) ;
      return ;
    }


    // Search for smallest tetrahedra that is formed from the closest face
    // and contains the interpolant point
    int n3 = -1 ;
    vect3d v2 = loc[neighbors[n2]]-v0 ;

    // If projected very close to face, then default to face interpolation
    // otherwise search for smallest tetrahedra that also contains the
    // interpolation point
    if(dist > 1e-10*scale) {
      double minvol = 1e33 ;
      // Now find smallest volume tetrahedra that contains interpolant point
      double vol3 = dot(cross(v1,v2),vc) ;
      // Make sure volumes we compare are positive
      if(vol3 < 0) {
        std::swap(v1,v2) ;
        std::swap(n1,n2) ;
        vol3 = -vol3 ;
      }
      for(int i=0;i<sz;++i)
        if(i!= n0 && i!= n1 && i != n2) {
          vect3d v3 = loc[neighbors[i]]-v0 ;
          double vol = dot(cross(v1,v2),v3) ;
          double vol1 = dot(cross(vc,v2),v3) ;
          double vol2 = dot(cross(v1,vc),v3) ;
          double vol0 = vol-vol1-vol2-vol3 ;

          // Check to see if point is inside tetrahedra
          // It is inside if the all of the sub-volumes are postive
          if(vol0>0 && vol1> 0 && vol2>0 && vol > 0)
            if(vol < minvol) { // if smaller vol, update n3
              minvol = vol ;
              n3 = i ;
            }
        }
    }

    // If we didn't find a tetrahedra then interpolate to the closest face
    if(n3 == -1) { // Project onto triangle
      vector<int> N(3) ;
      vector<double> W(3) ;
      N[0] = neighbors[n0] ;
      N[1] = neighbors[n1] ;
      N[2] = neighbors[n2] ;
      vect3d n = cross(v1,v2) ;
      double a = norm(n) ;
      n *= 1./a ;
      double d = dot(vc,n) ; // projected distance to triangle
      vect3d pnt = (vc-d*n) ;
          // Now check to see if pnt is in triangle
      double w2 = dot(n,cross(v1,pnt)) ;
      double w1 = dot(n,cross(pnt,v2)) ;
      double w0 = dot(n,cross(v1-pnt,v2-pnt)) ;
      double wsr = 1./(w0+w1+w2) ;
      W[0] = w0*wsr ;
      W[1] = w1*wsr ;
      W[2] = w2*wsr ;
      w.swap(W) ;
      neighbors.swap(N) ;
      return ;
    }


    // Compute tetrahedra barycentric coordinates for weights
    // and modify stencil to only include the four points that define
    // the bounding tetrahedra.
    vector<int> N(4) ;
    vector<double> W(4) ;
    N[0] = neighbors[n0] ;
    N[1] = neighbors[n1] ;
    N[2] = neighbors[n2] ;
    N[3] = neighbors[n3] ;

    vect3d v3 = loc[neighbors[n3]]-v0 ;
    double w3 = dot(cross(v1,v2),vc) ;
    double w2 = dot(cross(v1,vc),v3) ;
    double w1 = dot(cross(vc,v2),v3) ;
    double vol = dot(cross(v1,v2),v3) ;
    double w0 = max(vol-w1-w2-w3,0.0) ;
    double rvol = 1./(w0+w1+w2+w3) ;
    W[0] = w0*rvol ;
    W[1] = w1*rvol ;
    W[2] = w2*rvol ;
    W[3] = w3*rvol ;


    w.swap(W) ;
    neighbors.swap(N) ;
    return ;
  }

  void read_puT_file_serial(std::string filename,
                 const EOS &eos,
                 store<vect3d> &loc,
                 store<real> &p,
                 store<real> &T,
                 store<vect3d> &vel,
                 store<real> &k,
                 store<real> &mu_t,
                 storeVec<real> &mixture,
                 vector<int> &mixids) {

    ifstream bf(filename.c_str(),ios::in) ;
      if(bf.fail()) {
        cerr << "open failed on '" << filename <<"'"<< endl ;
        Loci::Abort() ;
      }

      Loci::parse::kill_white_space(bf) ;
      int nsp ;
      bf >> nsp ;
      for(int i=0;i<nsp;++i) {
        Loci::parse::kill_white_space(bf) ;
        string name = Loci::parse::get_name(bf) ;
        if(Loci::parse::get_token(bf,"(+)"))
          name += "(+)" ;
        else if(Loci::parse::get_token(bf,"(-)"))
          name += "(-)" ;
        mixids.push_back(eos.speciesIndex(name)) ;
        if(mixids.back() < 0) {
          Loci::debugout << "species name " << name << " not in simulation!"
                         << endl ;
        }
      }

      Loci::parse::kill_white_space(bf) ;
      int np ;
      bf >> np ; // read in number of points on the inflow boundary.

      entitySet dom = interval(0,np-1) ;

      loc.allocate(dom) ;
      p.allocate(dom) ;
      T.allocate(dom) ;
      vel.allocate(dom) ;
      k.allocate(dom) ;
      mu_t.allocate(dom) ;
      mixture.allocate(dom) ;
      mixture.setVecSize(nsp) ;

      Loci::parse::kill_white_space(bf) ;

      for(int i=0;i<np;++i) {
        bf >> loc[i]  // read in position of sample
           >> vel[i]  // velocity vector
           >> p[i]    // pressure
           >> T[i]    // temperature
           >> k[i]    // Turbulent kinetic energy
           >> mu_t[i]    // omega
          ;
        for(int j=0;j<nsp;++j)
          bf >> mixture[i][j] ;
      }
  }

  void broadcast_storeRep(Loci::storeRepP rep) {
    if(Loci::MPI_processes == 1)
      return ;
    entitySet domain = rep->domain() ;
    int size = rep->pack_size(domain) ;
    unsigned char *my_stuff = new unsigned char[size] ;
    int loc_pack = 0 ;
    if(Loci::MPI_rank == 0) {
      rep->pack(my_stuff,loc_pack,size,domain) ;
      MPI_Bcast(my_stuff,size,MPI_PACKED,0,MPI_COMM_WORLD) ;
    } else {
      MPI_Bcast(my_stuff,size,MPI_PACKED,0,MPI_COMM_WORLD) ;
      rep->unpack(my_stuff,loc_pack,size,domain) ;
    }
    delete[] my_stuff ;
  }


  void read_puT_file(std::string filename,
                     const EOS &eos,
                     store<vect3d> &loc,
                     store<real> &p,
                     store<real> &T,
                     store<vect3d> &vel,
                     store<real> &k,
                     store<real> &mu_t,
                     storeVec<real> &mix,
                     vector<int> &mixids) {
    const int ns = eos.numSpecies() ;
    int *data = new int[2+ns] ;
    if(Loci::MPI_rank == 0) {
      read_puT_file_serial(filename,eos,loc,p,T,vel,k,mu_t,mix,mixids) ;
      if(Loci::MPI_processes > 1) {
        int npnts = loc.domain().size() ;
        data[0] = npnts ;
        data[1] = mixids.size() ;
        for(unsigned int i=0;i<mixids.size();++i)
          data[i+2] = mixids[i] ;

        MPI_Bcast(data,2+mixids.size(),MPI_INT,0,MPI_COMM_WORLD) ;
        broadcast_storeRep(loc.Rep()) ;
        broadcast_storeRep(vel.Rep()) ;
        broadcast_storeRep(p.Rep()) ;
        broadcast_storeRep(T.Rep()) ;
        broadcast_storeRep(k.Rep()) ;
        broadcast_storeRep(mu_t.Rep()) ;
        broadcast_storeRep(mix.Rep()) ;

      }
    } else {
      int npnts ;
      MPI_Bcast(data,2+ns,MPI_INT,0,MPI_COMM_WORLD) ;
      npnts = data[0] ;
      int cnt = data[1] ;
      for(int i=0;i<cnt;++i)
        mixids.push_back(data[2+i]) ;

      entitySet dom = interval(0,npnts-1) ;
      loc.allocate(dom) ;
      vel.allocate(dom) ;
      p.allocate(dom) ;
      T.allocate(dom) ;
      k.allocate(dom) ;
      mu_t.allocate(dom) ;
      mix.allocate(dom) ;
      mix.setVecSize(mixids.size()) ;
      broadcast_storeRep(loc.Rep()) ;
      broadcast_storeRep(vel.Rep()) ;
      broadcast_storeRep(p.Rep()) ;
      broadcast_storeRep(T.Rep()) ;
      broadcast_storeRep(k.Rep()) ;
      broadcast_storeRep(mu_t.Rep()) ;
      broadcast_storeRep(mix.Rep()) ;

    }
    delete[] data ;
  }

  void read_scalar_file_serial(std::string filename,
                               store<vect3d> &loc,
                               store<real> &val) {

    ifstream bf(filename.c_str(),ios::in) ;
      if(bf.fail()) {
        cerr << "open failed on '" << filename <<"'"<< endl ;
        Loci::Abort() ;
      }

      Loci::parse::kill_white_space(bf) ;
      int np ;
      bf >> np ; // read in number of points on the inflow boundary.

      entitySet dom = interval(0,np-1) ;

      loc.allocate(dom) ;
      val.allocate(dom) ;

      Loci::parse::kill_white_space(bf) ;

      for(int i=0;i<np;++i) {
        bf >> loc[i]  // read in position of sample
           >> val[i] ; //  for this point value
      }
  }

  void read_scalar_file(std::string filename,
                        store<vect3d> &loc,
                        store<real> &val) {
    if(Loci::MPI_rank == 0) {
      read_scalar_file_serial(filename,loc,val) ;
      if(Loci::MPI_processes > 1) {
        int npnts = loc.domain().size() ;
        MPI_Bcast(&npnts,1,MPI_INT,0,MPI_COMM_WORLD) ;
        broadcast_storeRep(loc.Rep()) ;
        broadcast_storeRep(val.Rep()) ;
      }
    } else {
      int npnts ;
      MPI_Bcast(&npnts,1,MPI_INT,0,MPI_COMM_WORLD) ;

      entitySet dom = interval(0,npnts-1) ;
      loc.allocate(dom) ;
      val.allocate(dom) ;
      broadcast_storeRep(loc.Rep()) ;
      broadcast_storeRep(val.Rep()) ;
    }
  }

  class interpolate_data {
  public:
    Loci::kdTree::kd_tree *kd ;
    store<vector3d<double> > pos ;
    storeVec<double> data ;
    vector<int> distribution ;
    vector<std::string> species ;
    interpolate_data &operator=(const interpolate_data &in) {
      cerr << "interpolate_data shouldn't be copied!" << endl ;
      kd = 0 ;
      pos.setRep(in.pos.Rep()) ;
      data.setRep(in.data.Rep()) ;
      distribution = in.distribution ;
      return *this ;
    }
    ~interpolate_data() {
      if(kd !=0)
        delete kd ;
      kd = 0 ;
      data.allocate(EMPTY) ;
      pos.allocate(EMPTY) ;
    }
    interpolate_data() {
      kd = 0 ;
    }
  } ;



  class interpolate_file_read : public blackbox_rule {
    const_param<string> filename ;
    const_param<real> Pambient ;
    blackbox<interpolate_data> interp_data ;
  public:
    interpolate_file_read() {
      name_store("Pambient",Pambient) ;
      input("Pambient") ;
      name_store("FILENAME",filename) ;
      name_store("interpolateFile(FILENAME)",interp_data) ;
      input("FILENAME") ;
      output("interpolateFile(FILENAME)") ;
    }
    void compute(const sequence &seq) ;
  } ;

  void readContainerSimple(hid_t file_id, std::string vname,
                           Loci::storeRepP var) {
    Loci::readContainerRAW(file_id,vname,var,MPI_COMM_WORLD) ;
  }
  void interpolate_file_read::compute(const sequence &seq) {
    hid_t file_id = Loci::hdf5OpenFile((*filename).c_str(),
                                       H5F_ACC_RDONLY,H5P_DEFAULT) ;
    store<vect3d> pos ;

    readContainerSimple(file_id,"pos",pos.Rep()) ;
    entitySet dom = pos.domain() ;
    vector<vector3d<float> > spos(dom.size()) ;
    int cnt=0;
    FORALL(dom,nd) {
      spos[cnt++] = vector3d<float>(pos[nd].x,pos[nd].y,pos[nd].z) ;
    } ENDFORALL ;

    const int p = Loci::MPI_processes ;
    // Redistribute data
    vector<int> send_count(p,0) ;

    vector<int> procid ;
    if(p> 1)
      Loci::ORBPartition(spos,procid,MPI_COMM_WORLD) ;
    else {
      vector<int> tmp(spos.size(),0) ;
      procid.swap(tmp) ;
    }

    for(size_t i=0;i<procid.size();++i)
      send_count[procid[i]]++ ;


    vector<int> recv_count(p,0) ;

    MPI_Alltoall(&send_count[0],1,MPI_INT,&recv_count[0],1,MPI_INT,
                 MPI_COMM_WORLD) ;

    int final_size = 0 ;
    for(int i=0;i<p;++i)
      final_size += recv_count[i] ;

    vector<int> dist_sizes(p,0) ;
    MPI_Allgather(&final_size,1,MPI_INT,
                  &dist_sizes[0],1,MPI_INT,MPI_COMM_WORLD) ;


    vector<int> offsets(p+1,0) ;
    for(int i=1;i<p+1;++i)
      offsets[i] = offsets[i-1]+dist_sizes[i-1] ;

    interp_data->distribution = offsets ;
    vector<int> soffsets(p,0) ;
    for(int i=1;i<p;++i)
      soffsets[i] = soffsets[i-1] + send_count[i-1] ;

    int npos = pos.domain().size() ;
    vector<vector3d<double> > dpos(npos) ;
    vector<int> counts = soffsets ;
    cnt = 0 ;
    FORALL(dom,nd) {
      int pid = procid[cnt++] ;
      dpos[counts[pid]++] = pos[nd] ;
    } ENDFORALL ;

    int r = Loci::MPI_rank ;
    entitySet mydom = interval(offsets[r],offsets[r]+final_size-1) ;
    interp_data->pos.allocate(mydom) ;

    vector<int> send_num(p),recv_num(p) ;
    for(int i=0;i<p;++i) {
      send_num[i] = send_count[i] ;
      recv_num[i] = recv_count[i] ;
      send_count[i] *= 3 ;
      recv_count[i] *= 3 ;
    }
    vector<int> send_displacement(p,0) ;
    vector<int> recv_displacement(p,0) ;
    send_displacement[0] = 0;
    recv_displacement[0] = 0 ;
    for(int i=1;i<p;++i) {
      send_displacement[i] = send_displacement[i-1]+send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1]+recv_count[i-1] ;
    }

    int rbase = mydom.Min() ;

    MPI_Alltoallv(&dpos[0],
                  &send_count[0],&send_displacement[0], MPI_DOUBLE,
                  &(interp_data->pos[rbase].x),
                  &recv_count[0],&recv_displacement[0], MPI_DOUBLE,
                  MPI_COMM_WORLD) ;

    store<real> t,pg,k,tmuu ;
    store<vect3d> u ;
    storeVec<real> mixture ;

    readContainerSimple(file_id,"t",t.Rep()) ;
    readContainerSimple(file_id,"pg",pg.Rep()) ;

    readContainerSimple(file_id,"u",u.Rep()) ;
    param<std::string> sn ;
    Loci::readContainer(file_id,"speciesNames",sn.Rep(),EMPTY) ;
    vector<std::string> species ;
    string s = *sn ;
    string sp ;
    for(size_t i=0;i<s.size();++i) {
      if(s[i] == ':') {
        species.push_back(sp) ;
        sp = "" ;
      } else
        sp += s[i] ;
    }
    species.push_back(sp) ;
    interp_data->species = species ;

    if(species.size() > 1) {
      readContainerSimple(file_id,"mixture",mixture.Rep()) ;
    } else {
      mixture.setVecSize(0) ;
    }

    param<string> turbulence_model ;
    Loci::readContainer(file_id,"turbulence_model",turbulence_model.Rep()
                        ,EMPTY) ;
    bool has_turbulence_data = false ;
    if(*turbulence_model != "none") {
      has_turbulence_data = true ;
      readContainerSimple(file_id,"k",k.Rep()) ;
      readContainerSimple(file_id,"tmuu",tmuu.Rep()) ;
    }


    param<real> filePambient ;
    Loci::readContainer(file_id,"Pambient",filePambient.Rep(),EMPTY) ;

    Loci::hdf5CloseFile(file_id) ;

    const int base_size = 7 ;
    int vec_size = base_size ; // Number of values for pg, t, and u(x,y,z)

    if(mixture.vecSize() > 1)
      vec_size += mixture.vecSize() ;

    storeVec<real> data_vec ;
    entitySet local_dom = interval(0,dom.size()-1) ;
    data_vec.allocate(local_dom) ;
    data_vec.setVecSize(vec_size) ;
    interp_data->data.allocate(mydom) ;
    interp_data->data.setVecSize(vec_size) ;


    real padjust = *filePambient - *Pambient ;

    counts = soffsets ;
    cnt = 0 ;
    const int ms = vec_size-base_size ;
    FORALL(dom,ii) {
      int pid = procid[cnt] ;
      int addr = counts[pid]++ ;
      cnt++ ;
      data_vec[addr][0] = t[ii] ;
      data_vec[addr][1] = pg[ii] + padjust ;
      data_vec[addr][2] = u[ii].x ;
      data_vec[addr][3] = u[ii].y ;
      data_vec[addr][4] = u[ii].z ;
      data_vec[addr][5] = 0 ;
      data_vec[addr][6] = 0 ;
      if(has_turbulence_data) {
        data_vec[addr][5] = k[ii] ;
        data_vec[addr][6] = tmuu[ii] ;
      }

      for(int i=0;i<ms;++i)
        data_vec[addr][base_size+i] = mixture[ii][i] ;
    } ENDFORALL ;

    for(int i=0;i<p;++i) {
      send_count[i] = send_num[i]*vec_size ;
      recv_count[i] = recv_num[i]*vec_size ;
    }
    send_displacement[0] = 0;
    recv_displacement[0] = 0 ;
    for(int i=1;i<p;++i) {
      send_displacement[i] = send_displacement[i-1]+send_count[i-1] ;
      recv_displacement[i] = recv_displacement[i-1]+recv_count[i-1] ;
    }

   MPI_Alltoallv(&(data_vec[0][0]),
                 &send_count[0],&send_displacement[0], MPI_DOUBLE,
                 &(interp_data->data[rbase][0]),
                 &recv_count[0],&recv_displacement[0], MPI_DOUBLE,
                 MPI_COMM_WORLD) ;

   // Create kd_tree
    vector<Loci::kdTree::coord3d> pnts(mydom.size()) ;
    vector<int> pnt_id(mydom.size()) ;
    cnt = 0 ;
    FORALL(mydom,ii) {
      pnt_id[cnt] = ii ;
      pnts[cnt][0] = interp_data->pos[ii].x ;
      pnts[cnt][1] = interp_data->pos[ii].y ;
      pnts[cnt][2] = interp_data->pos[ii].z ;
      cnt++ ;
    } ENDFORALL ;
    interp_data->kd = new Loci::kdTree::kd_tree(pnts,pnt_id) ;
  }

  register_rule<interpolate_file_read> register_interpolate_file_read ;

  void getStencilBoundingBox(Loci::kdTree::kd_tree::bounds &bnd,
                             real &delta,
                             const const_store<vect3d> &pnts,
                             entitySet dom) {
    for(int d=0;d<3;++d) {
      bnd.minc[d] = .25*std::numeric_limits<double>::max() ;
      bnd.maxc[d] = -.25*std::numeric_limits<double>::max() ;
    }
    FORALL(dom,cc) {
      bnd.maxc[0] = max(bnd.maxc[0],pnts[cc].x) ;
      bnd.maxc[1] = max(bnd.maxc[1],pnts[cc].y) ;
      bnd.maxc[2] = max(bnd.maxc[2],pnts[cc].z) ;
      bnd.minc[0] = min(bnd.minc[0],pnts[cc].x) ;
      bnd.minc[1] = min(bnd.minc[1],pnts[cc].y) ;
      bnd.minc[2] = min(bnd.minc[2],pnts[cc].z) ;
    } ENDFORALL ;

    Loci::kdTree::kd_tree::bounds bndall ;
    MPI_Allreduce(&bnd.maxc[0],&bndall.maxc[0],3,MPI_DOUBLE,MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&bnd.minc[0],&bndall.minc[0],3,MPI_DOUBLE,MPI_MIN,
                  MPI_COMM_WORLD);
    int domsize = dom.size() ;
    int npnts = 0 ;
    MPI_Allreduce(&domsize,&npnts,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
    // Compute a delta to expand the bounding box
    double d1 = bndall.maxc[0]-bndall.minc[0] ;
    double d2 = bndall.maxc[1]-bndall.minc[1] ;
    double d3 = bndall.maxc[2]-bndall.minc[2] ;
    if(d1<d2)
      std::swap(d1,d2) ;
    if(d1<d3)
      std::swap(d1,d3) ;
    if(d2<d3)
      std::swap(d2,d3) ;

    // Compute mean distance
    real rnpnts = 1./double(npnts) ;
    real dist = max(max(d1*rnpnts,sqrt(d1*d2*rnpnts)),
               pow(d1*d2*d3*rnpnts,0.333333333)) ;

    // Over estimate to try to get the distance required to find stencil pnts
    dist *= 5 ;
    dist = max(dist,delta) ;
    delta = dist ;
    for(int d=0;d<3;++d) {
      bnd.maxc[d] += dist ;
      bnd.minc[d] -= dist ;
    }
  }

  void collectPoints(vector<Loci::kdTree::kd_tree::coord_info> &pout,
                     const Loci::kdTree::kd_tree &kd,
                     Loci::kdTree::kd_tree::bounds bnd) {
    // Communicate bounds request to other processors
    using namespace Loci::kdTree ;
    int p = Loci::MPI_processes ;
    vector<kd_tree::bounds> bnd_req(p) ;
    MPI_Allgather(&bnd,6,MPI_DOUBLE,&bnd_req[0],6,MPI_DOUBLE,MPI_COMM_WORLD) ;


    // Now communicate points that processors need to build stencil
    vector<int> scounts(p,0) ;
    vector<kd_tree::coord_info> pntlist ;

    for(int i=0;i<p;++i) { // Find points in i'th processor bounding box
      int bsz = pntlist.size() ;
      kd.find_box(pntlist,bnd_req[i]) ;
      scounts[i] = pntlist.size()-bsz ;
    }

    for(size_t i=0;i<scounts.size();++i)
      scounts[i]*=sizeof(kd_tree::coord_info) ;

    vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i)
      sdispls[i] = sdispls[i-1]+scounts[i-1] ;

    vector<int> rcounts(p) ;
    MPI_Alltoall(&scounts[0],1,MPI_INT,&rcounts[0],1,MPI_INT,MPI_COMM_WORLD) ;

    vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+rcounts[i-1] ;
    }

    int result_size = (rdispls[p-1]+rcounts[p-1])/sizeof(kd_tree::coord_info) ;

    vector<kd_tree::coord_info> pcollect(result_size) ;

    MPI_Alltoallv(&pntlist[0],&scounts[0],&sdispls[0],MPI_BYTE,
                  &pcollect[0],&rcounts[0],&rdispls[0],MPI_BYTE,
                  MPI_COMM_WORLD) ;
    pout.swap(pcollect) ;
  }

  void getCommSchedFromStencil(vector<int> &send_info_out,
                               vector<int> &req_sizes_out,
                               vector<int> &snd_sizes_out,
                               vector<int> &access_out,
                               const vector<int> &stencil,
                               const store<int> &ids,
                               const vector<int> &distribution) {
    vector<int> stmp = stencil;
    std::sort(stmp.begin(),stmp.end()) ;
    vector<int>::const_iterator se = std::unique(stmp.begin(),stmp.end()) ;
    vector<int> access(se-stmp.begin()) ;
    
    int cnt = 0 ;
    for(vector<int>::const_iterator ii=stmp.begin();ii!=se;++ii)
      access[cnt++] = ids[*ii] ;
    std::sort(access.begin(),access.end()) ;

    const int p = Loci::MPI_processes ;
    // Now communicate the accessed info
    vector<int> req_sizes(p,0) ;
    // Count accesses to each processor
    cnt = 0 ;
    for(size_t i=0;i<access.size();++i) {
      while(access[i] >= distribution[cnt+1]&& cnt<p)
        cnt++ ;
      req_sizes[cnt]++ ;
    }

    vector<int> snd_sizes(p,0) ;
    MPI_Alltoall(&req_sizes[0],1,MPI_INT,&snd_sizes[0],1,MPI_INT,
                 MPI_COMM_WORLD) ;

    int snd_tot_size = 0 ;
    for(int i=0;i<p;++i)
      snd_tot_size += snd_sizes[i] ;
    vector<int> send_info(snd_tot_size) ;

    vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i)
      sdispls[i] = sdispls[i-1]+req_sizes[i-1] ;

    vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+snd_sizes[i-1] ;
    }
    MPI_Alltoallv(&access[0],&req_sizes[0],&sdispls[0],MPI_INT,
                  &send_info[0],&snd_sizes[0],&rdispls[0],MPI_INT,
                  MPI_COMM_WORLD) ;

    send_info_out.swap(send_info) ;
    req_sizes_out.swap(req_sizes) ;
    snd_sizes_out.swap(snd_sizes) ;
    access_out.swap(access) ;
  }

  void remapStencil(vector<int> &stencils,
                    const vector<int> &access,
                    const store<int> &ids) {
    if(stencils.size() == 0)
      return ;
    entitySet locdom = ids.domain() ;

    vector<pair<int,int> > idmap(locdom.size()) ;

    FORALL(locdom,ii) {
      idmap[ii].first = ids[ii] ;
      idmap[ii].second = ii ;
    } ENDFORALL ;

    std::sort(idmap.begin(),idmap.end()) ;
    vector<int> acmap(locdom.size(),-1) ;
    size_t cnt = 0 ;
    for(int i=0;i<locdom.size();++i) {
      if(cnt >= access.size())
        break ;
      if(access[cnt] == idmap[i].first)
        acmap[idmap[i].second] = cnt++ ;
    }
    for(size_t i=0;i<stencils.size();++i) {
      if(acmap[stencils[i]] < 0)
        cerr << "bad stencil for remap" << endl ;
      stencils[i] = acmap[stencils[i]] ;
    }
  }

  void sendStencilData(store<double> &stencilData,
                       const_store<double> &sourceData,
                       const vector<int> &send_info,
                       const vector<int> &req_sizes_in,
                       const vector<int> &snd_sizes_in) {

    vector<double> databuf(send_info.size()) ;
    for(size_t i = 0;i<send_info.size();++i) {
      int id = send_info[i] ;
      databuf[i] = sourceData[id] ;
    }

    int p = Loci::MPI_processes ;
    vector<int> req_sizes(p),snd_sizes(p) ;
    for(int i=0;i<p;++i) {
      req_sizes[i] = req_sizes_in[i] ;
      snd_sizes[i] = snd_sizes_in[i] ;
    }

    vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i)
      sdispls[i] = sdispls[i-1]+req_sizes[i-1] ;

    vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+snd_sizes[i-1] ;
    }

    int loc_size = sdispls[p-1]+req_sizes[p-1] ;
    
    stencilData.allocate(entitySet(interval(0,loc_size-1))) ;

    MPI_Alltoallv(&databuf[0],&snd_sizes[0],&rdispls[0],MPI_DOUBLE,
                  &stencilData[0],&req_sizes[0],&sdispls[0],MPI_DOUBLE,
                  MPI_COMM_WORLD) ;
  }

  void sendStencilData(store<vector3d<double> > &stencilData,
                       const_store<vector3d<double> > &sourceData,
                       const vector<int> &send_info,
                       const vector<int> &req_sizes_in,
                       const vector<int> &snd_sizes_in) {

    vector<double> databuf(send_info.size()*3) ;
    for(size_t i = 0;i<send_info.size();++i) {
      int id = send_info[i] ;
      databuf[i*3+0] = sourceData[id].x ;
      databuf[i*3+1] = sourceData[id].y ;
      databuf[i*3+2] = sourceData[id].z ;
    }

    int p = Loci::MPI_processes ;
    vector<int> req_sizes(p),snd_sizes(p) ;
    for(int i=0;i<p;++i) {
      req_sizes[i] = req_sizes_in[i]*3 ;
      snd_sizes[i] = snd_sizes_in[i]*3 ;
    }

    vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i)
      sdispls[i] = sdispls[i-1]+req_sizes[i-1] ;

    vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+snd_sizes[i-1] ;
    }

    int loc_size = 0 ;
    for(int i=0;i<p;++i)
      loc_size += req_sizes_in[i] ;

    stencilData.allocate(entitySet(interval(0,loc_size-1))) ;

    MPI_Alltoallv(&databuf[0],&snd_sizes[0],&rdispls[0],MPI_DOUBLE,
                  &stencilData[0].x,&req_sizes[0],&sdispls[0],MPI_DOUBLE,
                  MPI_COMM_WORLD) ;
  }

  void sendStencilData(storeVec<double> &stencilData,
                       const_storeVec<double> &sourceData,
                       const vector<int> &send_info,
                       const vector<int> &req_sizes_in,
                       const vector<int> &snd_sizes_in) {

    int vec_size = sourceData.vecSize() ;
    vector<double> databuf(send_info.size()*vec_size) ;
    for(size_t i = 0;i<send_info.size();++i) {
      int id = send_info[i] ;
      for(int j=0;j<vec_size;++j) {
        databuf[i*vec_size+j] = sourceData[id][j] ;
      }
    }

    int p = Loci::MPI_processes ;
    vector<int> req_sizes(p),snd_sizes(p) ;
    for(int i=0;i<p;++i) {
      req_sizes[i] = req_sizes_in[i]*vec_size ;
      snd_sizes[i] = snd_sizes_in[i]*vec_size ;
    }

    vector<int> sdispls(p) ;
    sdispls[0] = 0 ;
    for(int i=1;i<p;++i)
      sdispls[i] = sdispls[i-1]+req_sizes[i-1] ;

    vector<int> rdispls(p) ;
    rdispls[0] = 0 ;
    for(int i=1;i<p;++i) {
      rdispls[i] = rdispls[i-1]+snd_sizes[i-1] ;
    }

    int loc_size = 0 ;
    for(int i=0;i<p;++i)
      loc_size += req_sizes_in[i] ;

    stencilData.allocate(entitySet(interval(0,loc_size-1))) ;
    stencilData.setVecSize(vec_size) ;

    MPI_Alltoallv(&databuf[0],&snd_sizes[0],&rdispls[0],MPI_DOUBLE,
                  &stencilData[0][0],&req_sizes[0],&sdispls[0],MPI_DOUBLE,
                  MPI_COMM_WORLD) ;
  }

  class interpolate_Data : public pointwise_rule {
    const_blackbox<interpolate_data> interp_data ;
    const_store<vect3d> pnts ;
    const_param<EOS> eos ;
    const_param<double> interpolateMinStencilSize ;
    storeVec<real> data ;
  public:
    interpolate_Data() {
      name_store("interpolateMinStencilSize",interpolateMinStencilSize) ;
      input("interpolateMinStencilSize") ;
      name_store("interpolateData(DATA,PNTS)",data) ;
      output("interpolateData(DATA,PNTS)") ;
      name_store("eos",eos) ;
      input("eos") ;
      name_store("PNTS",pnts) ;
      input("PNTS") ;
      name_store("DATA",interp_data) ;
      input("DATA") ;
      disable_threading() ;
    }
    void compute(const sequence &seq) ;

  } ;
  void interpolate_Data::compute(const sequence &seq) {
    using Loci::debugout ;

    entitySet dom = entitySet(seq) ;
    int vec_size = interp_data->data.vecSize() ;
    int base_size = 7 ;
    int out_vec_size = base_size+eos->numSpecies() ;
    data.setVecSize(out_vec_size) ;

    Loci::kdTree::kd_tree::bounds bnd ;

    real dist = *interpolateMinStencilSize ;
    getStencilBoundingBox(bnd,dist,pnts,dom) ;

    // Communicate bounds request to other processors
    using namespace Loci::kdTree ;

    vector<kd_tree::coord_info> pcollect ;

    collectPoints(pcollect,*(interp_data->kd),bnd) ;

    int result_size = pcollect.size() ;

    store<vect3d> loc ;
    store<int> ids ;
    entitySet locdom = interval(0,result_size-1) ;
    loc.allocate(locdom) ;
    ids.allocate(locdom) ;
    for(int i=0;i<result_size;++i) {
      loc[i].x = pcollect[i].coords[0] ;
      loc[i].y = pcollect[i].coords[1] ;
      loc[i].z = pcollect[i].coords[2] ;
      ids[i] = pcollect[i].id ;
      pcollect[i].id = i ;
    }
    kd_tree stree(pcollect) ;

    vector<double> weights ;
    vector<int> stencils ;
    vector<int> stencil_sizes ;
    FORALL(dom,cc) {
      vector<int> neighbors = get_stencil(stree,pnts[cc],dist) ;
#ifdef WALL_CONSTRAINT
      // remove any points on the wrong side of a viscous wall
      vect3d fcenter = facecenter[min_cell2noslip[*si]] ;
      vect3d n = area[min_cell2noslip[*si]].n ;
      vector<int> nn ;
      for(size_t i=0;i<neighbors.size();++i)
        if(dot(n,loc[neighbors[i]]-fcenter) <= 0.0)
          nn.push_back(neighbors[i]) ;
      if(nn.size() != 0)
        neighbors.swap(nn) ;
#endif
      // compute stencil weights
      vector<double> w  ;
      stencil_weights(w,neighbors,loc,pnts[cc]) ;

      stencil_sizes.push_back(neighbors.size()) ;

      for(size_t i=0;i<neighbors.size();++i) {
        stencils.push_back(neighbors[i]) ;
        weights.push_back(w[i]) ;
      }
    } ENDFORALL ;

    vector<int> send_info, req_sizes, snd_sizes, access ;


    getCommSchedFromStencil(send_info,req_sizes,snd_sizes, access,
                            stencils,ids,interp_data->distribution) ;

    remapStencil(stencils,access, ids) ;

    const_storeVec<double> tmp ;
    tmp.setRep(interp_data->data.Rep()) ;
    storeVec<double> valx ;
    sendStencilData(valx,tmp,send_info,req_sizes,snd_sizes) ;

    int cnt = 0 ;
    int c2 = 0 ;
    int ninput_species = vec_size-base_size ;
    vector<int> sid(ninput_species) ;
    for(int i=0;i<ninput_species;++i) {
      sid[i] = eos->speciesIndex(interp_data->species[i]) ;
      if(sid[i] < 0) {
        cerr << "species name " << interp_data->species[i]
             << " not in simulation!" << endl ;
        cerr << "problem with interpolated input file" << endl ;
        Loci::Abort() ;
      }

    }
    FORALL(dom,cc) {
      int sz = stencil_sizes[cnt++] ;
      for(int i=0;i<out_vec_size;++i)
        data[cc][i] = 0 ;
      for(int i=0;i<base_size;++i) {
        double dval = 0 ;
        for(int j=0;j<sz;++j) {
          dval += weights[c2+j]*valx[stencils[c2+j]][i] ;
        }
        data[cc][i] = dval ;
      }
      if(ninput_species <2)
        data[cc][base_size] = 1 ;
      for(int i=0;i<ninput_species;++i) {
        double dval = 0 ;
        for(int j=0;j<sz;++j) {
          dval += weights[c2+j]*valx[stencils[c2+j]][i+base_size] ;
        }
        data[cc][base_size+sid[i]] = dval ;
      }
      c2 += sz ;
    } ENDFORALL ;
  }

  register_rule<interpolate_Data> register_interpolate_Data ;

  class interpolate_DataCell : public pointwise_rule {
    const_blackbox<interpolate_data> interp_data ;
    const_store<vect3d> pnts ;
    const_param<EOS> eos ;

    const_Map min_cell2noslip ;
    const_store<vect3d> facecenter ;
    const_store<Area> area ;

    const_param<double> interpolateMinStencilSize ;

    storeVec<real> data ;
  public:
    interpolate_DataCell() {
      name_store("min_cell2noslip",min_cell2noslip) ;
      name_store("facecenter",facecenter) ;
      name_store("area",area) ;
      input("min_cell2noslip->(facecenter,area)") ;
      name_store("interpolateMinStencilSize",interpolateMinStencilSize) ;
      input("interpolateMinStencilSize") ;
       name_store("interpolateDataCell(DATA,PNTS)",data) ;
       output("interpolateDataCell(DATA,PNTS)") ;
       name_store("eos",eos) ;
       input("eos") ;
       name_store("PNTS",pnts) ;
       input("PNTS") ;
       name_store("DATA",interp_data) ;
       input("DATA") ;
       disable_threading() ;
    }
    void compute(const sequence &seq) ;

  } ;
  void interpolate_DataCell::compute(const sequence &seq) {
    using Loci::debugout ;

    entitySet dom = entitySet(seq) ;
    int vec_size = interp_data->data.vecSize() ;
    int base_size = 7 ;
    int out_vec_size = base_size+eos->numSpecies() ;
    data.setVecSize(out_vec_size) ;

    Loci::kdTree::kd_tree::bounds bnd ;

    real dist = *interpolateMinStencilSize ;
    getStencilBoundingBox(bnd,dist,pnts,dom) ;

    // Communicate bounds request to other processors
    using namespace Loci::kdTree ;

    vector<kd_tree::coord_info> pcollect ;

    collectPoints(pcollect,*(interp_data->kd),bnd) ;

    int result_size = pcollect.size() ;

    store<vect3d> loc ;
    store<int> ids ;
    entitySet locdom = interval(0,result_size-1) ;
    loc.allocate(locdom) ;
    ids.allocate(locdom) ;
    for(int i=0;i<result_size;++i) {
      loc[i].x = pcollect[i].coords[0] ;
      loc[i].y = pcollect[i].coords[1] ;
      loc[i].z = pcollect[i].coords[2] ;
      ids[i] = pcollect[i].id ;
      pcollect[i].id = i ;
    }
    kd_tree stree(pcollect) ;

    vector<double> weights ;
    vector<int> stencils ;
    vector<int> stencil_sizes ;
    FORALL(dom,cc) {
      vector<int> neighbors = get_stencil(stree,pnts[cc],dist) ;

      // remove any points on the wrong side of a viscous wall
      vect3d fcenter = facecenter[min_cell2noslip[cc]] ;
      vect3d n = area[min_cell2noslip[cc]].n ;
      vector<int> nn ;
      for(size_t i=0;i<neighbors.size();++i)
        if(dot(n,loc[neighbors[i]]-fcenter) <= 0.0)
          nn.push_back(neighbors[i]) ;
      if(nn.size() != 0)
        neighbors.swap(nn) ;

      //      debugout << "neighbors.size() = " << neighbors.size() << endl;
      // compute stencil weights
      vector<double> w  ;
      stencil_weights(w,neighbors,loc,pnts[cc]) ;

      stencil_sizes.push_back(neighbors.size()) ;

      for(size_t i=0;i<neighbors.size();++i) {
        stencils.push_back(neighbors[i]) ;
        weights.push_back(w[i]) ;
      }
    } ENDFORALL ;

    vector<int> send_info, req_sizes, snd_sizes, access ;


    getCommSchedFromStencil(send_info,req_sizes,snd_sizes, access,
                            stencils,ids,interp_data->distribution) ;

    remapStencil(stencils,access, ids) ;

    const_storeVec<double> tmp ;
    tmp.setRep(interp_data->data.Rep()) ;
    storeVec<double> valx ;
    sendStencilData(valx,tmp,send_info,req_sizes,snd_sizes) ;

    int cnt = 0 ;
    int c2 = 0 ;
    int ninput_species = vec_size-base_size ;
    vector<int> sid(ninput_species) ;
    for(int i=0;i<ninput_species;++i) {
      sid[i] = eos->speciesIndex(interp_data->species[i]) ;
      if(sid[i] < 0) {
        cerr << "species name " << interp_data->species[i]
             << " not in simulation!" << endl ;
        cerr << "problem with interpolated input file" << endl ;
        Loci::Abort() ;
      }

    }
    FORALL(dom,cc) {
      int sz = stencil_sizes[cnt++] ;
      for(int i=0;i<out_vec_size;++i)
        data[cc][i] = 0 ;
      for(int i=0;i<base_size;++i) {
        double dval = 0 ;
        for(int j=0;j<sz;++j) {
          dval += weights[c2+j]*valx[stencils[c2+j]][i] ;
        }
        data[cc][i] = dval ;
      }
      if(ninput_species <2)
        data[cc][base_size] = 1 ;
      for(int i=0;i<ninput_species;++i) {
        double dval = 0 ;
        for(int j=0;j<sz;++j) {
          dval += weights[c2+j]*valx[stencils[c2+j]][i+base_size] ;
        }
        data[cc][base_size+sid[i]] = dval ;
      }
      c2 += sz ;
    } ENDFORALL ;
  }

  register_rule<interpolate_DataCell> register_interpolate_DataCell ;

  class interpolate_points {
  public:
    Loci::kdTree::kd_tree *kd ;
    store<vector3d<double> > pos ;
    store<int> posid ;
    vector<int> distribution ;
    interpolate_points &operator=(const interpolate_points &in) {
      cerr << "interpolate_data shouldn't be copied!" << endl ;
      kd = 0 ;
      pos.setRep(in.pos.Rep()) ;
      return *this ;
    }
    ~interpolate_points() {
      if(kd !=0)
        delete kd ;
      kd = 0 ;
      pos.allocate(EMPTY) ;
      posid.allocate(EMPTY) ;
    }
    interpolate_points() {
      kd = 0 ;
    }
  } ;

  class interpolatePoints : public blackbox_rule {
    const_store<vect3d> pnts ;
    blackbox<interpolate_points> interp_pnts ;
  public:
    interpolatePoints() {
      name_store("PNTS",pnts) ;
      name_store("interpolatePoints(PNTS)",interp_pnts) ;
      input("PNTS") ;
      output("interpolatePoints(PNTS)") ;
    }
    void compute(const sequence &seq) ;
  } ;

  void interpolatePoints::compute(const sequence &seq) {
    int r = Loci::MPI_rank ;
    int p = Loci::MPI_processes ;
    entitySet dom = pnts.domain() ;
    if(p > 1)
      dom &= Loci::exec_current_fact_db->get_distribute_info()->my_entities ;
    int sz = dom.size() ;
    vector<int> dist_sizes(p,0) ;
    MPI_Allgather(&sz,1,MPI_INT,&dist_sizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
    vector<int> offsets(p+1,0) ;
    for(int i=1;i<p+1;++i)
      offsets[i] = offsets[i-1]+dist_sizes[i-1] ;
    interp_pnts->distribution = offsets ;
    entitySet mydom = interval(offsets[r],offsets[r]+sz-1) ;

    interp_pnts->pos.allocate(mydom) ;
    interp_pnts->posid.allocate(mydom) ;
    int ncnt = offsets[r] ;
    FORALL(dom,i) {
      interp_pnts->pos[ncnt] = pnts[i] ;
      interp_pnts->posid[ncnt] = i ;
      ncnt++ ;
    } ENDFORALL ;
    // Create kd_tree
    vector<Loci::kdTree::coord3d> pnts(mydom.size()) ;
    vector<int> pnt_id(mydom.size()) ;
    int cnt = 0 ;
    FORALL(mydom,ii) {
      pnt_id[cnt] = ii ;
      pnts[cnt][0] = interp_pnts->pos[ii].x ;
      pnts[cnt][1] = interp_pnts->pos[ii].y ;
      pnts[cnt][2] = interp_pnts->pos[ii].z ;
      cnt++ ;
    } ENDFORALL ;
    interp_pnts->kd = new Loci::kdTree::kd_tree(pnts,pnt_id) ;
  }

  register_rule<interpolatePoints> register_interpolatePoints ;

  struct stencil_info {
    vector<double> weights ;
    vector<int> stencils ;
    vector<int> stencil_sizes ;
    vector<int> stencil_offsets ;
    vector<int> send_info, req_sizes, snd_sizes ;
    Loci::storeRepP slookup ;
  } ;

  class interpolateFaceStencil_Unit: public unit_rule {
    blackbox<stencil_info> interpolateFaceStencil ;
  public:
    interpolateFaceStencil_Unit() {
      name_store("interpolateFaceStencil",interpolateFaceStencil) ;
      output("interpolateFaceStencil") ;
      constraint("UNIVERSE") ;
      disable_threading() ;
    }
    void compute(const sequence &seq) {
    }
  } ;

  register_rule<interpolateFaceStencil_Unit> register_interpolateFaceStencil_Unit ;

  class interpolateFaceStencil_Apply :
    public apply_rule<blackbox<stencil_info>, Loci::NullOp<stencil_info> > {
    const_store<vect3d> facecenter,cellcenter ;
    const_store<Area> area ;
    const_Map ci ;
    const_param<double> interpolateMinStencilSize ;
    const_blackbox<interpolate_points> interp_pnts ;
    blackbox<stencil_info> interpolateFaceStencil ;
  public:
    interpolateFaceStencil_Apply() {
      name_store("facecenter",facecenter) ;
      name_store("cellcenter",cellcenter) ;
      name_store("area",area) ;
      name_store("ci",ci) ;
      name_store("interpolateMinStencilSize",interpolateMinStencilSize) ;
      name_store("interpolatePoints(cellcenter)",interp_pnts) ;
      input("facecenter,area,ci->cellcenter") ;
      input("interpolateMinStencilSize") ;
      input("interpolatePoints(cellcenter)") ;
      name_store("interpolateFaceStencil",interpolateFaceStencil) ;
      output("interpolateFaceStencil") ;
      disable_threading() ;
    }
    void compute(const sequence &seq) ;
  } ;

  void interpolateFaceStencil_Apply::compute(const sequence &seq) {
    using Loci::debugout ;

    Loci::stopWatch s ;
    s.start() ;
    entitySet dom = entitySet(seq) ;

    Loci::kdTree::kd_tree::bounds bnd ;

    for(int d=0;d<3;++d) {
      bnd.minc[d] = .25*std::numeric_limits<double>::max() ;
      bnd.maxc[d] = -.25*std::numeric_limits<double>::max() ;
    }
    real dist = 0 ;
    FORALL(dom,fc) {
      real locdist = norm(facecenter[fc]-cellcenter[ci[fc]]) ;
      dist = max(dist,locdist) ;
      bnd.maxc[0] = max(bnd.maxc[0],facecenter[fc].x) ;
      bnd.maxc[1] = max(bnd.maxc[1],facecenter[fc].y) ;
      bnd.maxc[2] = max(bnd.maxc[2],facecenter[fc].z) ;
      bnd.minc[0] = min(bnd.minc[0],facecenter[fc].x) ;
      bnd.minc[1] = min(bnd.minc[1],facecenter[fc].y) ;
      bnd.minc[2] = min(bnd.minc[2],facecenter[fc].z) ;
    } ENDFORALL ;
    real dist_g = 0 ;
    MPI_Allreduce(&dist,&dist_g,1,MPI_DOUBLE,MPI_MAX, MPI_COMM_WORLD);
    dist = dist_g ;
    for(int d=0;d<3;++d) {
      bnd.maxc[d] += dist*6 ;
      bnd.minc[d] -= dist*6 ;
    }

    Loci::debugout << "time to find bounding box: " << s.stop() << endl ;
    s.start() ;
    // Communicate bounds request to other processors
    using namespace Loci::kdTree ;

    vector<kd_tree::coord_info> pcollect ;

    collectPoints(pcollect,*(interp_pnts->kd),bnd) ;

    int result_size = pcollect.size() ;

    Loci::debugout << "time to collect points: " << s.stop() << endl ;

    s.start() ;
    store<vect3d> loc ;
    store<int> ids ;
    entitySet locdom = interval(0,result_size-1) ;
    loc.allocate(locdom) ;
    ids.allocate(locdom) ;
    for(int i=0;i<result_size;++i) {
      loc[i].x = pcollect[i].coords[0] ;
      loc[i].y = pcollect[i].coords[1] ;
      loc[i].z = pcollect[i].coords[2] ;
      ids[i] = pcollect[i].id ;
      pcollect[i].id = i ;
    }
    kd_tree stree(pcollect) ;

    Loci::debugout << "time to build tree: " << s.stop() << endl ;

    s.start() ;
    vector<double> weights ;
    vector<int> stencils ;
    vector<int> stencil_sizes ;
    Map slookup ;
    slookup.allocate(dom) ;
    int scnt = 0 ;
    FORALL(dom,fc) {
      slookup[fc] = scnt ;
      scnt++ ;
      vect3d fcenter = facecenter[fc] ;
      vect3d n = area[fc].n ;
      // compute local stencil distance as either the normal projected
      // distance from this side of the face to the interior cell, or
      // to the closest projected point found from the other side,
      // whichever is largest
      real locdist = dot(fcenter-cellcenter[ci[fc]],n) ;
      if(locdist < 0)
        cerr << "logic error " << endl ;
      // Search points in neighborhood of the face
      {
        vector<int> neighbors = get_stencil(stree,facecenter[fc],dist) ;
        real dist2 = dist ; // First set to largest distance
        for(size_t i=0;i<neighbors.size();++i) {
          real dp = dot(loc[neighbors[i]]-fcenter,n) ;
          if(dp > 0) {
            real d = norm(loc[neighbors[i]]-fcenter) ;
            dist2 = min(dist2,d) ;
          }
        }
        locdist = max(locdist,dist2) ;
      }

      locdist *= 1.1 ;  // Saftey factor to ensure proper upwinding
      vect3d prr = facecenter[fc]+n*(locdist*3.) ;
      vect3d pr = facecenter[fc]+n*locdist ;
      vect3d pl = facecenter[fc]-n*locdist ;
      // right right stencil
      {
        vector<int> neighbors = get_stencil(stree,prr,2*dist) ;

        // compute stencil weights
        vector<double> w  ;
        stencil_weights(w,neighbors,loc,prr) ;

        stencil_sizes.push_back(neighbors.size()) ;

        for(size_t i=0;i<neighbors.size();++i) {
          stencils.push_back(neighbors[i]) ;
          weights.push_back(w[i]) ;
        }
      }
      // right stencil
      {
        vector<int> neighbors = get_stencil(stree,pr,2*dist) ;

        // compute stencil weights
        vector<double> w  ;
        stencil_weights(w,neighbors,loc,pr) ;

        stencil_sizes.push_back(neighbors.size()) ;

        for(size_t i=0;i<neighbors.size();++i) {
          stencils.push_back(neighbors[i]) ;
          weights.push_back(w[i]) ;
        }
      }
      // left stencil
      {
        vector<int> neighbors = get_stencil(stree,pl,2*dist) ;

        // compute stencil weights
        vector<double> w  ;
        stencil_weights(w,neighbors,loc,pl) ;

        stencil_sizes.push_back(neighbors.size()) ;

        for(size_t i=0;i<neighbors.size();++i) {
          stencils.push_back(neighbors[i]) ;
          weights.push_back(w[i]) ;
        }
      }
    } ENDFORALL ;

    Loci::debugout << "time to build stencil: " << s.stop()<< endl ;
    s.start() ;
    
    vector<int> send_info, req_sizes, snd_sizes, access ;


    getCommSchedFromStencil(send_info,req_sizes,snd_sizes, access,
                            stencils,ids,interp_pnts->distribution) ;

    for(size_t i=0;i<send_info.size();++i)
      send_info[i] = interp_pnts->posid[send_info[i]] ;
    remapStencil(stencils,access, ids) ;


    vector<int> stencil_offsets(stencil_sizes.size()) ;
    if(stencil_sizes.size() > 0)
      stencil_offsets[0] = 0 ;
    for(size_t i=1;i<stencil_sizes.size();++i)
      stencil_offsets[i] = stencil_offsets[i-1]+stencil_sizes[i-1] ;

    Loci::debugout << "time to generate comm schedule: " << s.stop() << endl ;
    
    interpolateFaceStencil->weights.swap(weights) ;
    interpolateFaceStencil->stencils.swap(stencils) ;
    interpolateFaceStencil->stencil_sizes.swap(stencil_sizes) ;
    interpolateFaceStencil->stencil_offsets.swap(stencil_offsets) ;
    interpolateFaceStencil->send_info.swap(send_info) ;
    interpolateFaceStencil->snd_sizes.swap(snd_sizes) ;
    interpolateFaceStencil->req_sizes.swap(req_sizes) ;
    interpolateFaceStencil->slookup = slookup.Rep() ;
  }

  register_rule<interpolateFaceStencil_Apply> register_interpolateFaceStencil_Apply ;


  class interpolateFaceScalar : public pointwise_rule {
    const_Map ci ;
    const_store<real> X ;
    const_blackbox<stencil_info> interpolateFaceStencil ;
    store<real> fX ;
  public:
    interpolateFaceScalar() {
      name_store("interpolateFaceStencil",interpolateFaceStencil) ;
      input("interpolateFaceStencil") ;
      name_store("ci",ci) ;
      name_store("X",X) ;
      input("ci->X") ;
      name_store("interpolateFace(X)",fX) ;
      output("interpolateFace(X)") ;
      disable_threading() ;
    }
    void compute(const sequence &seq) {
      entitySet dom = entitySet(seq) ;
      Map slookup ;
      slookup = interpolateFaceStencil->slookup ;
      store<real> Xstencil ;
      sendStencilData(Xstencil,X,
                      interpolateFaceStencil->send_info,
                      interpolateFaceStencil->req_sizes,
                      interpolateFaceStencil->snd_sizes)  ;
      const vector<double> &weights = interpolateFaceStencil->weights ;
      const vector<int> &stencils = interpolateFaceStencil->stencils ;
      const vector<int> &stencil_sizes
        = interpolateFaceStencil->stencil_sizes ;
      const vector<int> &stencil_offsets
        = interpolateFaceStencil->stencil_offsets ;

      FORALL(dom,fc) {
        int cnt = slookup[fc]*3 ;

        int szll = stencil_sizes[cnt+0] ;
        int szl  = stencil_sizes[cnt+1] ;
        int szr  = stencil_sizes[cnt+2] ;

        int cll = stencil_offsets[cnt+0] ;
        int cl =  stencil_offsets[cnt+1] ;
        int cr =  stencil_offsets[cnt+2] ;

        double vll = 0 ;
        for(int j=0;j<szll;++j) {
          vll += weights[cll+j]*Xstencil[stencils[cll+j]] ;
        }
        double vl = 0 ;
        for(int j=0;j<szl;++j) {
          vl += weights[cl+j]*Xstencil[stencils[cl+j]] ;
        }
        double vr = 0 ;
        for(int j=0;j<szr;++j) {
          vr += weights[cr+j]*Xstencil[stencils[cr+j]] ;
        }
        real dv = (vl-vll) ;
        if(dv < 0.)
          dv += -1e-30 ;
        else
          dv += 1e-30 ;

        // Use Van Albada limiter to estimate second order upwind
        // extrapolation to the face
        const real r = max(0.0,(vr-vl)/dv) ;
        const real lim = (r+r*r)/(1.+r*r) ;
        fX[fc] = vl + 0.5*lim*(vl-vll) ;

      } ENDFORALL ;
    }
  } ;
  register_rule<interpolateFaceScalar> register_interpolateFaceScalar ;



  class interpolateFacevect3d : public pointwise_rule {
    const_Map ci ;
    const_store<vect3d> X ;
    const_blackbox<stencil_info> interpolateFaceStencil ;
    store<vect3d> fX ;
  public:
    interpolateFacevect3d() {
      name_store("interpolateFaceStencil",interpolateFaceStencil) ;
      input("interpolateFaceStencil") ;
      name_store("ci",ci) ;
      name_store("X",X) ;
      input("ci->X") ;
      name_store("interpolateFace_v3d(X)",fX) ;
      output("interpolateFace_v3d(X)") ;
      disable_threading() ;
    }
    void compute(const sequence &seq) {
      entitySet dom = entitySet(seq) ;
      Map slookup ;
      slookup = interpolateFaceStencil->slookup ;
      store<vect3d> Xstencil ;
      sendStencilData(Xstencil,X,
                      interpolateFaceStencil->send_info,
                      interpolateFaceStencil->req_sizes,
                      interpolateFaceStencil->snd_sizes)  ;
      const vector<double> &weights = interpolateFaceStencil->weights ;
      const vector<int> &stencils = interpolateFaceStencil->stencils ;
      const vector<int> &stencil_sizes
        = interpolateFaceStencil->stencil_sizes ;
      const vector<int> &stencil_offsets
        = interpolateFaceStencil->stencil_offsets ;

      FORALL(dom,fc) {
        int cnt = slookup[fc]*3 ;

        int szll = stencil_sizes[cnt+0] ;
        int szl  = stencil_sizes[cnt+1] ;
        int szr  = stencil_sizes[cnt+2] ;

        int cll = stencil_offsets[cnt+0] ;
        int cl =  stencil_offsets[cnt+1] ;
        int cr =  stencil_offsets[cnt+2] ;

        vect3d vll = vect3d(0,0,0) ;
        for(int j=0;j<szll;++j) {
          vll += weights[cll+j]*Xstencil[stencils[cll+j]] ;
        }
        vect3d vl = vect3d(0,0,0) ;
        for(int j=0;j<szl;++j) {
          vl += weights[cl+j]*Xstencil[stencils[cl+j]] ;
        }
        vect3d vr = vect3d(0,0,0) ;
        for(int j=0;j<szr;++j) {
          vr += weights[cr+j]*Xstencil[stencils[cr+j]] ;
        }

        vect3d dv = (vl-vll) ;
        if(dv.x < 0.)
          dv.x += -1e-30 ;
        else
          dv.x += 1e-30 ;

        if(dv.y < 0.)
          dv.y += -1e-30 ;
        else
          dv.y += 1e-30 ;

        if(dv.z < 0.)
          dv.z += -1e-30 ;
        else
          dv.z += 1e-30 ;

        // Use Van Albada limiter to estimate second order upwind
        // extrapolation to the face
        const real rx = (vr.x-vl.x)/dv.x ;
        const real limx= (rx+rx*rx)/(1.+rx*rx) ;
        fX[fc].x = vl.x + 0.5*limx*(vl.x-vll.x) ;
        const real ry = (vr.y-vl.y)/dv.y ;
        const real limy= (ry+ry*ry)/(1.+ry*ry) ;
        fX[fc].y = vl.y + 0.5*limy*(vl.y-vll.y) ;
        const real rz = (vr.z-vl.z)/dv.z ;
        const real limz= (rz+rz*rz)/(1.+rz*rz) ;
        fX[fc].z = vl.z + 0.5*limz*(vl.z-vll.z) ;
      } ENDFORALL ;
    }
  } ;
  register_rule<interpolateFacevect3d> register_interpolateFacevect3d ;

  class interpolateFaceMixture : public pointwise_rule {
    const_Map ci ;
    const_storeVec<real> X ;
    const_blackbox<stencil_info> interpolateFaceStencil ;
    storeVec<real> fX ;
  public:
    interpolateFaceMixture() {
      name_store("interpolateFaceStencil",interpolateFaceStencil) ;
      input("interpolateFaceStencil") ;
      name_store("ci",ci) ;
      name_store("X",X) ;
      input("ci->X") ;
      name_store("interpolateFace_M(X)",fX) ;
      output("interpolateFace_M(X)") ;
      disable_threading() ;
    }
    void compute(const sequence &seq) {
      const int vs = X.vecSize() ;
      fX.setVecSize(X.vecSize()) ;
      entitySet dom = entitySet(seq) ;
      if(vs == 1) {
        FORALL(dom,fc) {
          fX[fc][0] = 1.0 ;
        } ENDFORALL ;
        return ;
      }
      
      Map slookup ;
      slookup = interpolateFaceStencil->slookup ;
      storeVec<real> Xstencil ;
      sendStencilData(Xstencil,X,
                      interpolateFaceStencil->send_info,
                      interpolateFaceStencil->req_sizes,
                      interpolateFaceStencil->snd_sizes)  ;
      const vector<double> &weights = interpolateFaceStencil->weights ;
      const vector<int> &stencils = interpolateFaceStencil->stencils ;
      const vector<int> &stencil_sizes
        = interpolateFaceStencil->stencil_sizes ;
      const vector<int> &stencil_offsets
        = interpolateFaceStencil->stencil_offsets ;

      FORALL(dom,fc) {
        int cnt = slookup[fc]*3 ;

        int szll = stencil_sizes[cnt+0] ;
        int szl  = stencil_sizes[cnt+1] ;
        int szr  = stencil_sizes[cnt+2] ;

        int cll = stencil_offsets[cnt+0] ;
        int cl =  stencil_offsets[cnt+1] ;
        int cr =  stencil_offsets[cnt+2] ;

        real sum = 0 ;
        for(int i=0;i<vs;++i) {
          real vll = 0 ;
          for(int j=0;j<szll;++j) 
            vll += weights[cll+j]*Xstencil[stencils[cll+j]][i] ;

          real vl = 0 ;
          for(int j=0;j<szl;++j) 
            vl += weights[cl+j]*Xstencil[stencils[cl+j]][i] ;

          real vr = 0 ;
          for(int j=0;j<szr;++j) 
            vr += weights[cr+j]*Xstencil[stencils[cr+j]][i] ;
          
          real dv = (vl-vll) ;
          if(dv < 0.)
            dv += -1e-30 ;
          else
            dv += 1e-30 ;
          
          // Use Van Albada limiter to estimate second order upwind
          // extrapolation to the face
          const real r = (vr-vl)/dv ;
          const real lim = (r+r*r)/(1.+r*r) ;
          // Limit mixture fractions should be between 0 and 1
          fX[fc][i] = min(max(vl + 0.5*lim*(vl-vll),0.0),1.0) ;
          sum += fX[fc][i] ;
        }
        real rsum = sum==0?1.0:1./sum ;
        // renormalize
        for(int i=0;i<vs;++i) 
          fX[fc][i] *= rsum ;
        
      } ENDFORALL ;
    }
  } ;
  register_rule<interpolateFaceMixture> register_interpolateFaceMixture ;


}

