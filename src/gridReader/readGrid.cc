// Standard library includes.
#include <string>
#include <vector>
using std::string ;
using std::vector ;

// Loci includes.
#include "Tools/expr.h"

// Fluid physics library includes.
#include "chemistry_db.h"
using fluidPhysics::chemistry_db ;

// StreamUns includes.
#include "const.h"
#include "name_var.h"
#include "readGrid.h"
#include "referenceFrame.h"
#include "sciTypes.h"

namespace Loci {
  extern void memSpace(string s) ;
  extern void createEdgesPar(fact_db &facts) ;
  extern void distributed_inverseMap(dmultiMap &result,const dmultiMap
    &input_map,const entitySet &input_image,const entitySet &input_preimage,
    fact_db &facts) ;
}

namespace streamUns {

  using Loci::createEdgesPar ;
  using Loci::memSpace ;
  using Loci::distributed_inverseMap ;

  // This currently only works SINGLE PROCESSOR. Copy of Ed's function with
  // additions for node2edge.
  void CreateEdges(fact_db &facts) {
    using std::vector ;
    using std::pair ;
    using namespace Loci ;

    multiMap face2node ;
    face2node = facts.get_variable("face2node") ;
    entitySet faces = face2node.domain() ;

    // Loop over faces and create list of edges (with duplicates)
    entitySet::const_iterator ei = faces.begin();
    vector<pair<Entity,Entity> > emap ;
    for(ei = faces.begin();ei!=faces.end();++ei) {
      int sz = face2node[*ei].size() ;
      for(int i=0;i<sz-1;++i) {
        Entity e1 = face2node[*ei][i] ;
        Entity e2 = face2node[*ei][i+1] ;
        emap.push_back(pair<Entity,Entity>(min(e1,e2),max(e1,e2))) ;
      }
      Entity e1 = face2node[*ei][0] ;
      Entity e2 = face2node[*ei][sz-1] ;
      emap.push_back(pair<Entity,Entity>(min(e1,e2),max(e1,e2))) ;
    }

    // Sort edges and remove duplicates
    sort(emap.begin(),emap.end()) ;
    vector<pair<Entity,Entity> >::iterator uend=unique(emap.begin(),
      emap.end()) ;

    // Allocate entities for new edges
    int num_edges = uend - emap.begin() ;
    cout << "num_edges = " << num_edges << endl ;
    entitySet edges = facts.get_distributed_alloc(num_edges).first ;

    // Copy edge nodes into a MapVec
    MapVec<2> edge ;
    edge.allocate(edges) ;
    vector<pair<Entity,Entity> >::iterator pi = emap.begin() ;
    for(ei=edges.begin();ei!=edges.end();++ei,++pi) {
      edge[*ei][0] = pi->first ;
      edge[*ei][1] = pi->second ;
    }

    // Add edge2node data structure to fact databse
    facts.create_fact("edge2node",edge) ;

    // Now create face2edge data-structure
    // We need to create a lower node to edge mapping to facilitate the
    // searches.  First get map from edge to lower node
    Map el ; // Lower edge map
    el.allocate(edges) ;
    for(ei=edges.begin();ei!=edges.end();++ei,++pi) {
      el[*ei] = edge[*ei][0] ;
    }

    // Now invert this map to get nodes-> edges that have this as a first entry
    multiMap n2e ;
    // Get nodes
    store<vector3d<real_t> > pos ;
    pos = facts.get_variable("pos") ;
    entitySet nodes = pos.domain() ;
    // Get mapping from nodes to edges from lower numbered node
    //    distributed_inverseMap(n2e, el, nodes, edges, facts) ;
    inverseMap(n2e, el, nodes, edges) ;

    // Section added by JW to create full node2edge map.
#ifdef MESH_DEFORMATION
    multiMap node2edge ;
    inverseMap(node2edge,edge,nodes,edges) ;
    facts.create_fact("node2edge",node2edge) ;
#endif

    // Now create face2edge map with same size as face2node
    multiMap face2edge ;
    store<int> count ;
    count.allocate(faces) ;
    for(ei = faces.begin();ei!=faces.end();++ei)
      count[*ei] = face2node[*ei].size() ;
    face2edge.allocate(count) ;

    // Now loop over faces, for each face search for matching edge and
    // store in the new face2edge structure
    for(ei = faces.begin();ei!=faces.end();++ei) {
      int sz = face2node[*ei].size() ;
      // Loop over edges of the face
      for(int i=0;i<sz-1;++i) {
        Entity t1 = face2node[*ei][i] ;
        Entity t2 = face2node[*ei][i+1] ;
        Entity e1 = min(t1,t2) ;
        Entity e2 = max(t1,t2) ;
        face2edge[*ei][i] = -1 ;
        // search for matching edge
        for(int j=0;j<n2e[e1].size();++j) {
          int e = n2e[e1][j] ;
          if(edge[e][0] == e1 && edge[e][1] == e2)
            face2edge[*ei][i] = e ;
        }
        if(face2edge[*ei][i] == -1)
          cerr << "not able to find edge for face " << *ei << endl ;
      }
      // Work on closing edge
      Entity t1 = face2node[*ei][0] ;
      Entity t2 = face2node[*ei][sz-1] ;
      Entity e1 = min(t1,t2) ;
      Entity e2 = max(t1,t2) ;
      face2edge[*ei][sz-1] = -1 ;
      for(int j=0;j<n2e[e1].size();++j) {
        int e = n2e[e1][j] ;
        if(edge[e][0] == e1 && edge[e][1] == e2)
          face2edge[*ei][sz-1] = e ;
      }
      if(face2edge[*ei][sz-1] == -1)
        cerr << "not able to find edge for face " << *ei << endl ;
    }
    // Add face2edge to the fact database
    facts.create_fact("face2edge",face2edge) ;
  }

  // Creates the node-to-edge map. Right now only good for serial.
  void CreateNodeToEdgePar(fact_db &facts) {

    // Get the edge2node map from the facts database and copy it to a
    // regular multiMap that can be used with distributed_inverseMap().
    MapVec<2> edge2node ; edge2node=facts.get_variable("edge2node") ;
    entitySet edges=edge2node.domain() ;
    multiMap e2n ;
    store<int> count ; count.allocate(edges) ;
    for(entitySet::const_iterator ei=edges.begin();ei!=edges.end();++ei) {
      count[*ei]=2 ;
    }
    e2n.allocate(count) ;
    for(entitySet::const_iterator ei=edges.begin();ei!=edges.end();++ei) {
      e2n[*ei][0]=edge2node[*ei][0] ; e2n[*ei][1]=edge2node[*ei][1] ;
    }

    // Invert the e2n multiMap to get n2e.
#ifdef MESH_DEFORMATION
    multiMap n2e ;
    store<vector3d<real_t> > pos ; pos=facts.get_variable("pos") ;
    entitySet nodes=pos.domain() ;
    entitySet global_nodes=all_collect_entitySet(nodes,facts) ;
    entitySet global_edges=all_collect_entitySet(edges,facts) ;
    distributed_inverseMap(n2e,e2n,global_nodes,global_edges,facts) ;
    facts.create_fact("node2edge",n2e) ;
#endif
  }

  // Peeks into the chemistry model file to see how many species there are.
  unsigned int PeekChemistry(ifstream &input) {
    chemistry_db cdb_defaults,cdb ; char dbname[512] ;
    const char *base = getenv("CHEMISTRY_DATABASE") ;
    if(base==NULL) {
      if(Loci::MPI_rank==0)
        cerr << "Warning: CHEMISTRY_DATABASE environment variable not defined"            << endl ;
      base = "./" ;
    }
    sprintf(dbname,"%s/data_base/chemistry",base) ;
    ifstream cdf(dbname,ios::in) ;
    if(!cdf.fail()){
      cdb_defaults.Input(cdf) ; cdf.close() ;
    }else{
      if(Loci::MPI_rank==0)
        cerr << "warning, can't open '" << dbname << "'" << endl ;
    }
    cdb.Input(input,cdb_defaults) ; cdb.removeM() ;
    return cdb.species.numSpecies() ;
  }

  // Performs some simple consistency checks on the .vars file.
  void CheckVarsFile(fact_db &facts) {

    // Check boundary conditions.
    if(facts.get_fact("boundary_conditions")==0){
      cerr << "ERROR: There are no boundary conditions!" << endl ;
      Loci::Abort() ;
    }

    // Check initial condition.
    {
      int count=0 ;
      if(facts.get_fact("ql")!=0 && facts.get_fact("qr")==0){
        cerr << "ERROR: Missing qr!" << endl ; Loci::Abort() ;
      }
      if(facts.get_fact("qr")!=0 && facts.get_fact("ql")==0){
        cerr << "ERROR: Missing ql!" << endl ; Loci::Abort() ;
      }
      if(facts.get_fact("initialCondition")!=0) ++count ;
      if(facts.get_fact("initialConditionsFile")!=0) ++count ;
      if(facts.get_fact("initialConditionRegions")!=0) ++count ;
      if(facts.get_fact("interpolateInitialConditions")!=0){
        ++count ;

        // Need to specify density, so allow this combination.
        if(facts.get_fact("initialCondition")!=0) --count ;
      }
      if(facts.get_fact("ql")!=0) ++count ;
      if(count==0){
        cerr << "ERROR: Initial conditions missing!" << endl ; Loci::Abort() ;
      }
      if(count>1){
        cerr << "ERROR: Initial conditions overspecified!" << endl ;
        Loci::Abort() ;
      }
    }

    // See if transport properties have been specified.
    if(facts.get_fact("transport_model")==0){
      cerr << "ERROR: transport_model missing!" << endl ; Loci::Abort() ;
    }else{
      param<string> transportModel ;
      transportModel=facts.get_fact("transport_model") ;
      param<string> flowCompressibility ;
      flowCompressibility=facts.get_fact("flowCompressibility") ;
      if(*transportModel=="const_viscosity"){
        if(facts.get_fact("mu")==0){
          cerr << "ERROR: Must specify mu with transport_model: "
            << "const_viscosity!" << endl ; Loci::Abort() ;
        }
        if(*flowCompressibility=="compressible" && facts.get_fact("kcond")==0){
          cerr << "ERROR: Must specify kcond with transport_model: "
            << "const_viscosity!" << endl ; Loci::Abort() ;
        }
      }
    }

    // Make sure energy equation is not specified for incompressible flow and
    // is specified for compressible flows.
    {
      param<string> flowCompressibility ;
      flowCompressibility=facts.get_fact("flowCompressibility") ;
      if(*flowCompressibility=="incompressible" && facts.get_fact
      ("energyEquationOptions")!=0){
        cerr << "ERROR: Remove energyEquationOptions for incompressible "
          << "flow!" << endl ; Loci::Abort() ;
      }else if(*flowCompressibility=="compressible" && facts.get_fact
      ("energyEquationOptions")==0){
        cerr << "ERROR: Must specify energyEquationOptions for compressible "
          << "flow!" << endl ; Loci::Abort() ;
      }
    }
  }

  // Reads the .vars file and the grid file.
  void ReadGrid(fact_db &facts,const rule_db &rdb,const char *caseName) {

    // Grid file options.
    param<grid_options> grid_file_info ;
    facts.create_fact("grid_file_info",grid_file_info) ;

    // Add the case name to the fact database.
    param<string> modelName ; *modelName = string(caseName) ;
    facts.create_fact("modelName",modelName) ;

    // Output variable selection.
    param<string> plot_output ; *plot_output = "" ;
    facts.create_fact("plot_output",plot_output) ;

    // Open and read the .vars file.
    string varsFileName=caseName ; varsFileName+=".vars" ;
    try{
      ifstream in(varsFileName.c_str(),ios::in) ;
      if(in.fail()) {
        cerr << "Can't open " << varsFileName << endl ; Loci::Abort() ;
      }
      while(in.peek() != '{' && in.peek() != EOF) { in.get() ; }
      facts.read_vars(in,rdb) ;
    }catch(const Loci::BasicException &err){
      err.Print(cerr) ;
      cerr << "Aborted reading " << varsFileName << endl ; Loci::Abort() ;
    }

    // Do some simple test on .vars file consistency.
    CheckVarsFile(facts) ;

    // Add testing option to the global facts database.
    constraint testingConstraint ; testingConstraint = ~EMPTY ;
    param<string> testing ; testing=facts.get_fact("testing") ;
    facts.create_fact(testing->c_str(),testingConstraint) ;

    // This is a HACK to create the species transport constraints here by
    // peeking in the chemistry model file. No other way appears possible
    // right now. Have now added check to see if there is a chemistry model
    // specified, so we can run grid movement problems without having to
    // specify a chemistry model.
    if(facts.get_fact("chemistry_model")!=0){
      unsigned int numSpecies ;
      param<name_var> chemistryModel ;
      chemistryModel=facts.get_fact("chemistry_model") ;
      if((*chemistryModel).name==""){
        ifstream chemin("chem.in",ios::in) ;
        if(chemin.fail()){
          cerr << "Can't find chemistry model file: chem.in" << endl ;
          Loci::Abort() ;
        }
        numSpecies=PeekChemistry(chemin) ;
      }else{
        string fname=(*chemistryModel).name+".mdl" ;
        ifstream chemin(fname.c_str(),ios::in) ;
        if(chemin.fail()){
          const char *base=getenv("CHEMISTRY_DATABASE") ;
          if(base==NULL) {
            cerr << "Warning: CHEMISTRY_DATABASE environment variable "
              << "not defined" << endl ;
            base = "./" ;
          }
          fname=string(base)+"/data_base/models/"+(*chemistryModel).name+
            ".mdl" ;
          ifstream chemin(fname.c_str(),ios::in) ;
          if(chemin.fail()) {
            cerr << "Can't find chemistry model file: "
              << (*chemistryModel).name << fname << endl ; Loci::Abort() ;
          }else{
            numSpecies=PeekChemistry(chemin) ;
          }
        }else{
          numSpecies=PeekChemistry(chemin) ;
        }
      }

      if(numSpecies==1){
        constraint noSpeciesTransport ; noSpeciesTransport=~EMPTY ;
        facts.create_fact("noSpeciesTransport",noSpeciesTransport) ;
      }else{
        constraint speciesTransport ; speciesTransport=~EMPTY ;
        facts.create_fact("speciesTransport",speciesTransport) ;
      }
    }

    // Create the output constraints.
    if(*plot_output != ""){
      Loci::expression::exprP ep=Loci::expression::create(*plot_output) ;
      if(ep->op == Loci::OP_NAME){
        constraint x ; x = ~EMPTY ;
        string cname = string("scalarOutput_")+ep->name ;
        facts.create_fact(cname,x) ;
      }else if(ep->op == Loci::OP_COMMA){
        Loci::expression::exprList::const_iterator li ;
        for(li=ep->expr_list.begin();li!=ep->expr_list.end();++li) {
          if((*li)->op == Loci::OP_NAME) {
            constraint x ; x = ~EMPTY ;
            string cname=string("scalarOutput_")+(*li)->name ;
            facts.create_fact(cname,x) ;
          }else{
            cerr << "unable to interpret expression in 'plot_output'" << endl ;
            (*li)->Print(cerr) ; Loci::Abort() ;
          }
        }
      }else{
        cerr << "unable to interpret expression in 'plot_output'" << endl ;
        ep->Print(cerr) ; Loci::Abort() ;
      }
    }

    // Check the boundary conditions.
    if(CheckBoundaryConditions(facts)) {
      cerr << "WARNING: boundary condition errors detected!" << endl ;
      if((*grid_file_info).optionExists("bc_check")){
        string check ; (*grid_file_info).getOption("bc_check",check) ;
        if(check != "relaxed") Loci::Abort() ;
      }else{
        Loci::Abort() ;
      }
    }

    // The file_type is an old option.  Here we just check to see if it is
    // set, if so, print a diagnostic if a nonsupported filetype is selected.
    string extension=".xdr" ;
    if((*grid_file_info).optionExists("file_type")){
      string file_type ; (*grid_file_info).getOption("file_type",file_type) ;
      if(file_type == "NEW_BINARY_FORMATTED" || file_type == "XDR"){
      }else if(file_type == "VOG"){
        extension = ".vog" ;
      }else{
        cerr << "unknown file type " << file_type << endl ; Loci::Abort() ;
      }
    }

    // Read in the XDR file
    string file=string(caseName)+extension ;
    if(Loci::MPI_rank==0)
      cout << "Portable Grid File input, reading file = " << file << endl ;
    if(!Loci::setupFVMGrid(facts,file)){
      cerr << "Yikes! Reading grid file '" << file <<"' failed in grid reader!"
        << endl ;
      cerr << "aborting!" << endl ; Loci::Abort() ;
    }
    if(Loci::MPI_rank==0) cout << "Reading Grid File Complete" << endl ;

    // Scale the grid by Lref, if provided
    real Lref=1.0 ; vect3d translate=vect3d(0.,0.,0.) ;
    if((*grid_file_info).optionExists("translate")){
      get_vect3dOption(*grid_file_info,"translate","m",translate,1.0) ;
    }
    if((*grid_file_info).optionExists("Lref")){
      Loci::option_value_type ovt=(*grid_file_info).getOptionValueType("Lref") ;
      if(ovt==Loci::REAL) (*grid_file_info).getOption("Lref",Lref) ;
      if(ovt==Loci::UNIT_VALUE){
        Loci::UNIT_type vu ; (*grid_file_info).getOption("Lref",vu) ;
        if(!vu.is_compatible("m")){
          cerr << "wrong type of units for reference length Lref" << endl ;
        }else{
          Lref = vu.get_value_in("m") ;
        }
      }
      store<vect3d> pos ; pos=facts.get_fact("pos") ;
      entitySet posdom = pos.domain() ;
      if(Loci::MPI_rank==0) cout << "Scaling grid by " << Lref << endl ;
      for(entitySet::const_iterator ei=posdom.begin();ei!=posdom.end();++ei){
        pos[*ei]=pos[*ei]*Lref+translate ;
      }
    }

    // Add Lref to fact database
    param<real> Lrefparam ; *Lrefparam=Lref ;
    facts.create_fact("Lref",Lrefparam) ;

    // Create the axisymmetric constraint.
    if((*grid_file_info).optionExists("axisymmetric")){
      constraint axisymmetric ; axisymmetric=~EMPTY ;
      facts.create_fact("axisymmetric",axisymmetric) ;
    }

    // Perform additional setup.
    memSpace("Before setupBoundaryConditions().") ;
    setupBoundaryConditions(facts) ;
    memSpace("Before CreateLowerUpper().") ;
    createLowerUpper(facts) ;

    /* Create the noWallFunction constraint.
    constraint boundary_faces=facts.get_variable("boundary_faces") ;
    if(facts.get_variable("wallFunction_BCoption")!=0){
      constraint wallFunction_BC=facts.get_variable("wallFunction_BCoption") ;
      entitySet noWallFunction=*boundary_faces-*wallFunction_BC ;
      constraint noWallFunction_BC ; noWallFunction_BC=noWallFunction ;
      facts.create_fact("noWallFunction_BC",noWallFunction_BC) ;
    }else{
      facts.create_fact("noWallFunction_BC",boundary_faces) ;
    }*/

    // Create additional maps for moving meshes.
#ifdef MESH_DEFORMATION
    memSpace("Before CreateEdges().") ;
//  if(Loci::MPI_rank!=0){
//    cerr << "ERROR: Non-parallel version of edge2node and node2edge used "
//      << "with multiple processes." << endl ;
//    Loci::Abort() ;
//  }
//  CreateEdges(facts) ;
    createEdgesPar(facts) ;
    CreateNodeToEdgePar(facts) ;
#endif
    memSpace("Before FindMindNoslip().") ;
    FindMindNoslip(facts) ;
    memSpace("Returning from grid setup.") ;
  }
}

