// Standard library includes.
#include <string>
#include <sstream>
#include <stdlib.h>
using std::string ;
                                                                                
// Loci includes.
#include <Loci.h>

// Fluid physics library includes.
#include "eos.h"
using fluidPhysics::EOS ;

// StreamUns includes.
#include "initialCondition.h"
#include "sciTypes.h"
                                                                                
namespace streamUns {

  // Class to set the restart flag.
  class DoRestart : public singleton_rule {
    private:
      const_param<int> n,nCycle ;
      const_param<int> restartFrequency ;
      param<bool> doRestart ;
    public:
                                                                                
      // Define input and output.
      DoRestart() {
        name_store("$n{n}",n) ;
        name_store("ncycle{n}",nCycle) ;
        name_store("restart_freq{n}",restartFrequency) ;
        name_store("do_restart{n}",doRestart) ;
        input("$n{n},ncycle{n},restart_freq{n}");
        output("do_restart{n}") ;
      }
                                                                                
      // Set the restart flag.
      virtual void compute(const sequence &seq) {
        if((*n)==0){ doRestart=false ; return ; }
        doRestart=((*nCycle)%(*restartFrequency)==0) ;
      }
  } ;
                                                                                
  register_rule<DoRestart> registerDoRestart ;

  // Class for reading the solution time and time-step number from a
  // restart file.
  class ReadRestartData : public singleton_rule {
    private:
      const_param<string> restartNum ;
      param<int> initialTimeStepNum ;
      param<real> initialSolutionTime ;
    public:

      // Define input and output.
      ReadRestartData() {
        name_store("restartNum",restartNum) ;
        name_store("restart::ncycle_init",initialTimeStepNum) ;
        name_store("restart::stime_init",initialSolutionTime) ;
        input("restartNum") ;
        output("restart::ncycle_init,restart::stime_init") ;
        constraint("restart") ;
        disable_threading() ;
      }

      // Set the initial time-step number.
      virtual void compute(const sequence &seq) {
        string restartFile=string("restart/time_hdf5.")+ *restartNum ;
        if(Loci::MPI_rank==0) cout << "Reading solution time and time-step "
          << "number from restart file: " << restartFile << endl ;
        hid_t fileID=Loci::hdf5OpenFile((restartFile).c_str(),H5F_ACC_RDONLY,
          H5P_DEFAULT) ;
        entitySet dom = ~EMPTY ;
        Loci::readContainer(fileID,"ncycle",initialTimeStepNum.Rep(),dom) ;
        Loci::readContainer(fileID,"stime",initialSolutionTime.Rep(),dom) ;
        Loci::hdf5CloseFile(fileID) ;
        if(Loci::MPI_rank==0) cout << "Time Step Number: "
          << *initialTimeStepNum << "  Solution Time (s): "
          << *initialSolutionTime << endl ;
      }
  } ;

  register_rule<ReadRestartData> registerReadRestartData ;

  // Rule to read the initial condition for density from a restart file. We
  // need this now so that we can construct the species densities which are
  // used in determining eos_state_ic.
  class DensityInitialConditionRestart : public pointwise_rule {
    private:
      const_param<string> restartNum ;
      const_store<EOS::State> eos_state_ic ;
      store<real> rho_ic ;
    public:

      // Define input and output.
      DensityInitialConditionRestart() {
        name_store("restartNum",restartNum) ;
        name_store("eos_state_ic",eos_state_ic) ;
        name_store("restart::rho_ic",rho_ic) ;
        input("restartNum,eos_state_ic") ;
        output("restart::rho_ic") ;
        constraint("restart,geom_cells,compressibleFlow") ;
        disable_threading() ;
      }

      // Assign total enthalpy for a single cell.
      void calculate(Entity cell) {
        rho_ic[cell]=eos_state_ic[cell].density() ;
      }

      // Read the density values from the file restart/rho.hdf5
      void compute(const sequence &seq) {
        ostringstream oss ; oss << "restart/rho_hdf5." << *restartNum ;
        string fileName=oss.str() ;
        ifstream test(fileName.c_str()) ;
        if(test){
          test.close() ;
          if(Loci::MPI_rank==0) cout << "Reading density initial "
            << "condition from restart file: " << fileName << endl ;
          entitySet dom=entitySet(seq) ;
          hid_t fileID=Loci::hdf5OpenFile(fileName.c_str(),H5F_ACC_RDONLY,
            H5P_DEFAULT);
          Loci::readContainer(fileID,"rho",rho_ic.Rep(),dom) ;
          Loci::hdf5CloseFile(fileID) ;
        }else{
          if(Loci::MPI_rank==0) cout << "WARNING: No density restart file. "
            << "Computing density from equation of state." << endl ;

          // Loop through the sequence of cells.
          do_loop(seq,this) ;
        }
      }
  } ;
                                                                                
  register_rule<DensityInitialConditionRestart>
    registerDensityInitialConditionRestart ;

  // Rule to read the initial condition for velocity from a restart file.
  class VelocityInitialConditionRestart : public pointwise_rule {
    private:
      const_param<string> restartNum ;
      store<vect3d> v_ic ;
    public:
                                                                                
      // Define input and output.
      VelocityInitialConditionRestart() {
        name_store("restartNum",restartNum) ;
        name_store("restart::v_ic",v_ic) ;
        input("restartNum") ;
        output("restart::v_ic") ;
        constraint("restart,geom_cells") ;
        disable_threading() ;
      }
                                                                                
      // Read the velocity values from the file restart/v.hdf5
      void compute(const sequence &seq) {
        ostringstream oss ; oss << "restart/v_hdf5." << *restartNum ;
        string fileName=oss.str() ;
        if(Loci::MPI_rank==0) cout << "Reading velocity initial "
          << "condition from restart file: " << fileName << endl ;
        entitySet dom=entitySet(seq) ;
        hid_t fileID=Loci::hdf5OpenFile(fileName.c_str(),H5F_ACC_RDONLY,
          H5P_DEFAULT);
        Loci::readContainer(fileID,"v",v_ic.Rep(),dom) ;
        Loci::hdf5CloseFile(fileID) ;
      }
  } ;
                                                                                
  register_rule<VelocityInitialConditionRestart>
    registerVelocityInitialConditionRestart ;

  // Rule to read the initial condition for pressure from a restart file.
  class PressureInitialConditionRestart : public pointwise_rule {
    private:
      const_param<string> restartNum ;
      store<real> p_ic ;
    public:
                                                                                
      // Define input and output.
      PressureInitialConditionRestart() {
        name_store("restartNum",restartNum) ;
        name_store("restart::p_ic",p_ic) ;
        input("restartNum") ;
        output("restart::p_ic") ;
        constraint("restart,geom_cells") ;
        disable_threading() ;
      }
                                                                                
      // Read the pressure values from the file restart/p.hdf5
      void compute(const sequence &seq) {
        ostringstream oss ; oss << "restart/p_hdf5." << *restartNum ;
        string fileName=oss.str() ;
        if(Loci::MPI_rank==0) cout << "Reading pressure initial "
          << "condition from restart file: " << fileName << endl ;
        entitySet dom=entitySet(seq) ;
        hid_t fileID=Loci::hdf5OpenFile(fileName.c_str(),H5F_ACC_RDONLY,
          H5P_DEFAULT);
        Loci::readContainer(fileID,"p",p_ic.Rep(),dom) ;
        Loci::hdf5CloseFile(fileID) ;
      }
  } ;
                                                                                
  register_rule<PressureInitialConditionRestart>
    registerPressureInitialConditionRestart ;

  // Rule to read the initial condition for temperature from a restart file.
  class TemperatureInitialConditionRestart : public pointwise_rule {
    private:
      const_param<InitialCondition> initialCondition ;
      const_param<string> restartNum ;
      store<real> T_ic ;
    public:
                                                                                
      // Define input and output.
      TemperatureInitialConditionRestart() {
        name_store("initialCondition",initialCondition) ;
        name_store("restartNum",restartNum) ;
        name_store("restart::T_ic",T_ic) ;
        input("initialCondition,restartNum") ;
        output("restart::T_ic") ;
        constraint("restart,geom_cells") ;
        disable_threading() ;
      }

      // Assign temperature for a single cell.
      void calculate(Entity cell) {
        T_ic[cell]=initialCondition->Temperature() ;
      }
                                                                                
      // Read the temperature values from the file restart/T.hdf5
      void compute(const sequence &seq) {
        ostringstream oss ; oss << "restart/temperature_hdf5." << *restartNum ;
        string fileName=oss.str() ;
        ifstream test(fileName.c_str()) ;
        if(test){
          test.close() ;
          if(Loci::MPI_rank==0) cout << "Reading temperature initial "
            << "condition from restart file: " << fileName << endl ;
          entitySet dom=entitySet(seq) ;
          hid_t fileID=Loci::hdf5OpenFile(fileName.c_str(),H5F_ACC_RDONLY,
            H5P_DEFAULT);
          Loci::readContainer(fileID,"temperature",T_ic.Rep(),dom) ;
          Loci::hdf5CloseFile(fileID) ;
        }else{
          if(Loci::MPI_rank==0) cout << "WARNING: No temperature restart file. "
            << "Using temperature initial condition value from .vars file."
            << endl ;

          // Check that temperature was provided by the user.
          if(!initialCondition->IsTemperatureDefined())
            cerr << "ERROR: Initial condition for temperature required!"
              << endl ;
                                                                                
          // Loop through the sequence of cells.
          do_loop(seq,this) ;
        }
      }
  } ;
                                                                                
  register_rule<TemperatureInitialConditionRestart>
    registerTemperatureInitialConditionRestart ;

  // Rule to read the initial condition for total enthalpy from a restart file.
  class TotalEnthalpyInitialConditionRestart : public pointwise_rule {
    private:
      const_param<string> restartNum ;
      const_store<EOS::State> eos_state_ic ;
      const_store<vect3d> v_ic ;
      store<real> h_ic ;
    public:

      // Define input and output.
      TotalEnthalpyInitialConditionRestart() {
        name_store("restartNum",restartNum) ;
        name_store("eos_state_ic",eos_state_ic) ;
        name_store("v_ic",v_ic) ;
        name_store("restart::h_ic",h_ic) ;
        input("restartNum,eos_state_ic,v_ic") ;
        output("restart::h_ic") ;
        constraint("restart,compressibleFlow,geom_cells") ;
        disable_threading() ;
      }

      // Assign total enthalpy for a single cell.
      void calculate(Entity cell) {
        h_ic[cell]=eos_state_ic[cell].enthalpy()+0.5*dot(v_ic[cell],
          v_ic[cell]) ;
      }

      // Read the total enthalpy values from the file restart/h.hdf5
      void compute(const sequence &seq) {
        ostringstream oss ; oss << "restart/h_hdf5." << *restartNum ;
        string fileName=oss.str() ;
        ifstream test(fileName.c_str()) ;
        if(test){
          test.close() ;
          if(Loci::MPI_rank==0) cout << "Reading total enthalpy initial "
            << "condition from restart file: " << fileName << endl ;
          entitySet dom=entitySet(seq) ;
          hid_t fileID=Loci::hdf5OpenFile(fileName.c_str(),H5F_ACC_RDONLY,
            H5P_DEFAULT);
          Loci::readContainer(fileID,"h",h_ic.Rep(),dom) ;
          Loci::hdf5CloseFile(fileID) ;
        }else{
          if(Loci::MPI_rank==0) cout << "WARNING: No total enthalpy restart "
            <<"file. Computing total enthalpy from equation of state." << endl ;

          // Loop through the sequence of cells.
          do_loop(seq,this) ;
        }
      }
  } ;

  register_rule<TotalEnthalpyInitialConditionRestart>
    registerTotalEnthalpyInitialConditionRestart ;

  // Rule to read the initial condition for species mass fraction from a
  // restart file.
  class SpeciesInitialConditionRestart : public pointwise_rule {
    private:
      const_param<string> restartNum ;
      const_param<int> numSpecies ;
      storeVec<real> y_ic ;
    public:

      // Define input and output.
      SpeciesInitialConditionRestart() {
        name_store("restartNum",restartNum) ;
        name_store("numSpecies",numSpecies) ;
        name_store("restart::y_ic",y_ic) ;
        input("restartNum,numSpecies") ;
        output("restart::y_ic") ;
        constraint("restart,speciesTransport,geom_cells") ;
        disable_threading() ;
      }
                                                                                
      // Read the species values from the file restart/y.hdf5
      void compute(const sequence &seq) {
        ostringstream oss ; oss << "restart/y_hdf5." << *restartNum ;
        string fileName=oss.str() ;
        if(Loci::MPI_rank==0) cout << "Reading species initial condition "
          << "from restart file: " << fileName << endl ;
        entitySet dom=entitySet(seq) ; y_ic.setVecSize(*numSpecies) ;
        hid_t fileID=Loci::hdf5OpenFile(fileName.c_str(),H5F_ACC_RDONLY,
          H5P_DEFAULT);
        Loci::readContainer(fileID,"y",y_ic.Rep(),dom) ;
        Loci::hdf5CloseFile(fileID) ;
      }
  } ;
                                                                                
  register_rule<SpeciesInitialConditionRestart>
    registerSpeciesInitialConditionRestart ;

  // Rule to read the initial condition for k from a restart file.
  class KInitialConditionRestart : public pointwise_rule {
    private:
      const_param<InitialCondition> initialCondition ;
      const_param<string> restartNum ;
      store<real> k_ic ;
    public:
                                                                                
      // Define input and output.
      KInitialConditionRestart() {
        name_store("initialCondition",initialCondition) ;
        name_store("restartNum",restartNum) ;
        name_store("restart::k_ic",k_ic) ;
        input("initialCondition,restartNum") ;
        output("restart::k_ic") ;
        constraint("restart,geom_cells") ;
        disable_threading() ;
      }

      // Assign k for a single cell.
      void calculate(Entity cell) { k_ic[cell]=initialCondition->K() ; }
                                                                                
      // Read the k values from the file restart/k.hdf5
      void compute(const sequence &seq) {
        ostringstream oss ; oss << "restart/k_hdf5." << *restartNum ;
        string fileName=oss.str() ;
        ifstream test(fileName.c_str()) ;
        if(test){
          test.close() ;
          if(Loci::MPI_rank==0) cout << "Reading k initial "
            << "condition from restart file: " << fileName << endl ;
          entitySet dom=entitySet(seq) ;
          hid_t fileID=Loci::hdf5OpenFile(fileName.c_str(),H5F_ACC_RDONLY,
            H5P_DEFAULT);
          Loci::readContainer(fileID,"k",k_ic.Rep(),dom) ;
          Loci::hdf5CloseFile(fileID) ;
        }else{
          if(Loci::MPI_rank==0) cout << "WARNING: No k restart file. Using "
            << "k initial condition value from .vars file." << endl ;

          // Check that k was provided by the user.
          if(!initialCondition->IsKDefined())
            cerr << "ERROR: Initial condition for k required!" << endl ;
                                                                                
          // Loop through the sequence of cells.
          do_loop(seq,this) ;
        }
      }
  } ;
                                                                                
  register_rule<KInitialConditionRestart> registerKInitialConditionRestart ;

  // Rule to read the initial condition for omega from a restart file.
  class OmegaInitialConditionRestart : public pointwise_rule {
    private:
      const_param<InitialCondition> initialCondition ;
      const_param<string> restartNum ;
      store<real> omega_ic ;
    public:
                                                                                
      // Define input and output.
      OmegaInitialConditionRestart() {
        name_store("initialCondition",initialCondition) ;
        name_store("restartNum",restartNum) ;
        name_store("restart::omega_ic",omega_ic) ;
        input("initialCondition,restartNum") ;
        output("restart::omega_ic") ;
        constraint("restart,geom_cells") ;
        disable_threading() ;
      }

      // Assign omega for a single cell.
      void calculate(Entity cell) { omega_ic[cell]=initialCondition->Omega() ; }
                                                                                
      // Read the omega values from the file restart/omega.hdf5
      void compute(const sequence &seq) {
        ostringstream oss ; oss << "restart/omega_hdf5." << *restartNum ;
        string fileName=oss.str() ;
        ifstream test(fileName.c_str()) ;
        if(test){
          test.close() ;
          if(Loci::MPI_rank==0) cout << "Reading omega initial "
            << "condition from restart file: " << fileName << endl ;
          entitySet dom=entitySet(seq) ;
          hid_t fileID=Loci::hdf5OpenFile(fileName.c_str(),H5F_ACC_RDONLY,
            H5P_DEFAULT);
          Loci::readContainer(fileID,"omega",omega_ic.Rep(),dom) ;
          Loci::hdf5CloseFile(fileID) ;
        }else{
          if(Loci::MPI_rank==0) cout << "WARNING: No omega restart file. Using "
            << "omega initial condition value from .vars file." << endl ;

          // Check that omega was provided by the user.
          if(!initialCondition->IsOmegaDefined())
            cerr << "ERROR: Initial condition for omega required!" << endl ;
                                                                                
          // Loop through the sequence of cells.
          do_loop(seq,this) ;
        }
      }
  } ;
                                                                                
  register_rule<OmegaInitialConditionRestart>
    registerOmegaInitialConditionRestart ;

  // Rule to read the initial condition for mass flux from a restart file.
  class MassFluxInitialConditionRestart : public pointwise_rule {
    private:
      const_param<string> restartNum ;
      store<real> massFluxTimeStepZero ;
    public:
                                                                                
      // Define input and output.
      MassFluxInitialConditionRestart() {
        name_store("restartNum",restartNum) ;
        name_store("massFlux{n=0}",massFluxTimeStepZero) ;
        input("restartNum") ;
        output("massFlux{n=0}") ;
        constraint("restart,massFluxCorrected") ;
        disable_threading() ;
      }

      // Read the pressure values from the file restart/massFlux.hdf5
      void compute(const sequence &seq) {
        ostringstream oss ; oss << "restart/massFlux_hdf5." << *restartNum ;
        string fileName=oss.str() ;
        if(Loci::MPI_rank==0) cout << "Reading mass flux initial "
          << "condition from restart file: " << fileName << endl ;
        entitySet dom=entitySet(seq) ;
        hid_t fileID=Loci::hdf5OpenFile(fileName.c_str(),H5F_ACC_RDONLY,
          H5P_DEFAULT);
        Loci::readContainer(fileID,"massFlux",massFluxTimeStepZero.Rep(),dom) ;
        Loci::hdf5CloseFile(fileID) ;
      }
  } ;
                                                                                
//register_rule<MassFluxInitialConditionRestart>
//  registerMassFluxInitialConditionRestart ;

  // Class for writing out the solution time and time-step number to a
  // restart file.
  class WriteRestartData : public singleton_rule {
    private:
      const_param<int> timeStepNum ;
      const_param<real> solutionTime ;
      const_param<int> restart_modulo ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      WriteRestartData() {
        name_store("ncycle{n}",timeStepNum) ;
        name_store("stime{n}",solutionTime) ;
        name_store("restart_modulo{n}",restart_modulo) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        input("ncycle{n},stime{n},restart_modulo{n}") ;
        output("OUTPUT{n}") ;
        conditional("do_restart{n}") ;
        constraint("UNIVERSE") ;
      }

      // Set the initial time-step number.
      virtual void compute(const sequence &seq) {
        unsigned int fileExtension=(*restart_modulo!=0)? (*timeStepNum)%
          (*restart_modulo):(*timeStepNum) ;
        ostringstream oss ; oss << "restart/time_hdf5." << fileExtension ;
        string fileName=oss.str() ;
        if(Loci::MPI_rank==0) cout << "Writing restart file: " << fileName
          << endl ;
        hid_t fileID=Loci::hdf5CreateFile(fileName.c_str(),H5F_ACC_TRUNC,
          H5P_DEFAULT, H5P_DEFAULT) ;
        param<int> tempTimeStepNum ; *tempTimeStepNum=*timeStepNum ;
        param<real> tempSolutionTime ; *tempSolutionTime=*solutionTime ;
        Loci::writeContainer(fileID,"ncycle",tempTimeStepNum.Rep()) ;
        Loci::writeContainer(fileID,"stime",tempSolutionTime.Rep()) ;
        Loci::hdf5CloseFile(fileID) ;
      }
  } ;

  register_rule<WriteRestartData> registerWriteRestartData ;

  // Class for writing scalar restart data.
  class RestartScalar : public pointwise_rule {
    private:
      string variableName,constraintName ;
      const_param<int> timeStepNum ;
      const_store<real> variable ;
      const_param<int> restart_modulo ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      RestartScalar(const char *variableName,const char *constraintName) :
        variableName(variableName),constraintName(constraintName) {
        string fullVariableName=variableName+string("{n}") ;
        string fullConstraintName=constraintName+string("{n}") ;
        name_store("ncycle{n}",timeStepNum) ;
        name_store(fullVariableName,variable) ;
        name_store("restart_modulo{n}",restart_modulo) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        input(fullVariableName) ;
        input("ncycle{n},restart_modulo{n}") ;
        output("OUTPUT{n}") ;
        conditional("do_restart{n}") ;
        constraint("geom_cells{n}") ;
        constraint(fullConstraintName) ;
      }

      // Write the cell data to the restart file. Note the use of a temporary
      // store<real> to ensure that we only write out data for geom_cells. This
      // is required since .Rep() method that writes out the data writes out
      // all the data including periodic cells.
      void compute(const sequence &seq) {
        unsigned int fileExtension=(*restart_modulo!=0)? (*timeStepNum)%
          (*restart_modulo):(*timeStepNum) ;
        ostringstream oss ; oss << "restart/" << variableName << "_hdf5."
          << fileExtension ;
        string fileName=oss.str() ;
        if(Loci::MPI_rank==0) cout << "Writing restart file: " << fileName
          << endl ;
        hid_t file_id=Loci::hdf5CreateFile(fileName.c_str(),H5F_ACC_TRUNC,
          H5P_DEFAULT, H5P_DEFAULT) ;
        store<real> tempVariable ; tempVariable.allocate(entitySet(seq)) ;
        for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
          tempVariable[*si]=variable[*si] ;
        Loci::writeContainer(file_id,variableName,tempVariable.Rep()) ;
        Loci::hdf5CloseFile(file_id) ;
      }
  } ;

  // Class for writing vector restart data.
  class RestartVector : public pointwise_rule {
    private:
      string variableName,constraintName ;
      const_param<int> timeStepNum ;
      const_store<vect3d> variable ;
      const_param<int> restart_modulo ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      RestartVector(const char *variableName,const char *constraintName) :
        variableName(variableName),constraintName(constraintName) {
        string fullVariableName=variableName+string("{n}") ;
        string fullConstraintName=constraintName+string("{n}") ;
        name_store("ncycle{n}",timeStepNum) ;
        name_store(fullVariableName,variable) ;
        name_store("restart_modulo{n}",restart_modulo) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        input(fullVariableName) ;
        input("ncycle{n},restart_modulo{n}") ;
        output("OUTPUT{n}") ;
        conditional("do_restart{n}") ;
        constraint("geom_cells{n}") ;
        constraint(fullConstraintName) ;
      }

      // Write the cell data to the restart file. Note the use of a temporary
      // store<vect3d> to ensure that we only write out data for geom_cells.
      // This is required since .Rep() method that writes out the data writes
      // out all the data including periodic cells.
      void compute(const sequence &seq) {
        unsigned int fileExtension=(*restart_modulo!=0)? (*timeStepNum)%
          (*restart_modulo):(*timeStepNum) ;
        ostringstream oss ; oss << "restart/" << variableName << "_hdf5."
          << fileExtension ;
        string fileName=oss.str() ;
        if(Loci::MPI_rank==0) cout << "Writing restart file: " << fileName
          << endl ;
        hid_t file_id=Loci::hdf5CreateFile(fileName.c_str(),H5F_ACC_TRUNC,
          H5P_DEFAULT, H5P_DEFAULT) ;
        store<vect3d> tempVariable ; tempVariable.allocate(entitySet(seq)) ;
        for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
          tempVariable[*si]=variable[*si] ;
        Loci::writeContainer(file_id,variableName,tempVariable.Rep()) ;
        Loci::hdf5CloseFile(file_id) ;
      }
  } ;

  // Class for writing storeVec restart data.
  class RestartStoreVec : public pointwise_rule {
    private:
      string variableName,constraintName ;
      const_param<int> timeStepNum ;
      const_storeVec<real> variable ;
      const_param<int> restart_modulo ;
      param<bool> OUTPUT ;
    public:

      // Define input and output. Note that we must output y{n} and not
      // y{n,it} since the ODE integration may have occurred and we do not
      // want to include that contribution since this will mean that the
      // ODE will be doubly integrated at the timestep after restart.
      RestartStoreVec(const char *variableName,const char *constraintName) :
        variableName(variableName),constraintName(constraintName) {
        string fullVariableName=variableName+string("{n}") ;
        string fullConstraintName=constraintName+string("{n}") ;
        name_store("ncycle{n}",timeStepNum) ;
        name_store(fullVariableName,variable) ;
        name_store("restart_modulo{n}",restart_modulo) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        input(fullVariableName) ;
        input("ncycle{n},restart_modulo{n}") ;
        output("OUTPUT{n}") ;
        conditional("do_restart{n}") ;
        constraint("geom_cells{n}") ;
        constraint(fullConstraintName) ;
      }

      // Write the cell data to the restart file. Note the use of a temporary
      // storeVec<real> to ensure that we only write out data for geom_cells.
      // This is required since .Rep() method that writes out the data writes
      // out all the data including periodic cells.
      void compute(const sequence &seq) {
        unsigned int fileExtension=(*restart_modulo!=0)? (*timeStepNum)%
          (*restart_modulo):(*timeStepNum) ;
        ostringstream oss ; oss << "restart/" << variableName << "_hdf5."
          << fileExtension ;
        string fileName=oss.str() ;
        if(Loci::MPI_rank==0) cout << "Writing restart file: " << fileName
          << endl ;
        hid_t file_id=Loci::hdf5CreateFile(fileName.c_str(),H5F_ACC_TRUNC,
          H5P_DEFAULT, H5P_DEFAULT) ;
        storeVec<real> tempVariable ; tempVariable.allocate(entitySet(seq)) ;
        tempVariable.setVecSize(variable.vecSize()) ;
        for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
          tempVariable[*si]=variable[*si] ;
        Loci::writeContainer(file_id,variableName,tempVariable.Rep()) ;
        Loci::hdf5CloseFile(file_id) ;
      }
  } ;

// Macro which is used to create a restart output class for each scalar.
#define RESTART_SCALAR(X,Y) class Restart_##X : public RestartScalar {\
                              public:\
                                Restart_##X() : RestartScalar(#X,#Y){}\
                            } ;\
                            register_rule<Restart_##X> registerRestart_##X ;

// Macro which is used to create a restart output class for each vector.
#define RESTART_VECTOR(X,Y) class Restart_##X : public RestartVector {\
                              public:\
                                Restart_##X() : RestartVector(#X,#Y){}\
                            } ;\
                            register_rule<Restart_##X> registerRestart_##X ;

// Macro which is used to create a restart output class for each storeVec.
#define RESTART_STOREVEC(X,Y) class Restart_##X : public RestartStoreVec {\
                              public:\
                                Restart_##X() : RestartStoreVec(#X,#Y){}\
                            } ;\
                            register_rule<Restart_##X> registerRestart_##X ;

  // Create the scalar restart classes. For scalars that have no contraint,
  // just pass in UNIVERSE.
  RESTART_SCALAR(rho,compressibleFlow) ;
  RESTART_SCALAR(p,UNIVERSE) ;
  RESTART_SCALAR(temperature,compressibleFlow) ;
  RESTART_SCALAR(h,compressibleFlow) ;
  RESTART_SCALAR(k,kOmegaTurbulenceModel) ;
  RESTART_SCALAR(omega,kOmegaTurbulenceModel) ;

  // Create the vector restart classes. For vectors that have no contraint,
  // just pass in UNIVERSE.
  RESTART_VECTOR(v,UNIVERSE) ;

  // Create the storeVec restart classes. For storeVecs that have no contraint,
  // just pass in UNIVERSE.
  RESTART_STOREVEC(y,speciesTransport) ;

  // Class for writing mass flux restart data.
  class RestartMassFlux : public pointwise_rule {
    private:
      const_param<int> timeStepNum ;
      const_store<real> massFlux ;
      const_param<int> restart_modulo ;
      param<bool> OUTPUT ;
    public:

      // Define input and output.
      RestartMassFlux() {
        name_store("ncycle{n}",timeStepNum) ;
        name_store("massFlux{n}",massFlux) ;
        name_store("restart_modulo{n}",restart_modulo) ;
        name_store("OUTPUT{n}",OUTPUT) ;
        input("ncycle{n},massFlux{n}") ;
        input("restart_modulo{n}") ;
        output("OUTPUT{n}") ;
        conditional("do_restart{n}") ;
        constraint("massFluxCorrected{n}") ;
      }

      // Write the mass flux data to the restart file. Note the use of a
      // temporary store<real> to ensure that we only write out data over
      // "faces" and not any duplicate partition data. This is required since
      // the .Rep() method that writes out the data writes out all the data
      // including periodic cells.
      void compute(const sequence &seq) {
        unsigned int fileExtension=(*restart_modulo!=0)? (*timeStepNum)%
          (*restart_modulo):(*timeStepNum) ;
        ostringstream oss ; oss << "restart/" << "massFlux_hdf5."
          << fileExtension ;
        string fileName=oss.str() ;
        if(Loci::MPI_rank==0) cout << "Writing restart file: " << fileName
          << endl ;
        hid_t file_id=Loci::hdf5CreateFile(fileName.c_str(),H5F_ACC_TRUNC,
          H5P_DEFAULT, H5P_DEFAULT) ;
        Loci::writeContainer(file_id,"massFlux",massFlux.Rep()) ;
        Loci::hdf5CloseFile(file_id) ;
      }
  } ;

//register_rule<RestartMassFlux> registerRestartMassFlux ;

}
