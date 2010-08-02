// Loci includes.
#include <Loci.h>
                                                                                
// StreamUns includes.
#include "const.h"
#include "initialCondition.h"
#include "sciTypes.h"
#include "varsFileInputs.h"
                                                                                
namespace streamUns {

  // Sets the default background pressure.
  class DefaultBackgroundPressure : public default_rule {
    private:
      param<real> pAmbient ;
    public:

      // Define input and output.
      DefaultBackgroundPressure() {
        name_store("Pambient",pAmbient) ;
        output("Pambient") ;
        comments("Sets the default background pressure to 0.0.") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) { *pAmbient=0.0 ; }
  } ;

  register_rule<DefaultBackgroundPressure> registerDefaultBackgroundPressure ;

  // Plot frequency for boundary data. Set to high value so that the regular
  // plot_freq takes precedence unless this value is set lower.
  class DefaultBoundaryPlotFreq : public default_rule {
    private:
      param<int> boundary_plot_freq ;
    public:

      // Define input and output.
      DefaultBoundaryPlotFreq() {
        name_store("boundary_plot_freq",boundary_plot_freq) ;
        output("boundary_plot_freq") ;
        comments("Output boundary plot files whenever the timestep ") ;
        comments("number is divisible by this input.") ;
      }

      // Set the default value.
      virtual void compute(const sequence& seq) {
        *boundary_plot_freq=10000000 ;
      }
  } ;

  register_rule<DefaultBoundaryPlotFreq> registerDefaultBoundaryPlotFreq ;

  // Sets the default convergence tolerance.
  class DefaultConvergenceTolerance : public default_rule {
    private:
      param<real> convergenceTolerance ;
    public:

      // Define input and output.
      DefaultConvergenceTolerance() {
        name_store("convergenceTolerance",convergenceTolerance) ;
        output("convergenceTolerance") ;
      }

      // Set the default value.
      virtual void compute(const sequence& seq) {
        *convergenceTolerance=1.0e-03 ;
      }
  } ;

  register_rule<DefaultConvergenceTolerance>
    registerDefaultConvergenceTolerance ;

  // Sets the default CVODE convergence tolerance.
  class DefaultCVODETolerance : public default_rule {
    private:
      param<real> cvodeTolerance ;
    public:

      // Define input and output.
      DefaultCVODETolerance() {
        name_store("cvodeTolerance",cvodeTolerance) ;
        output("cvodeTolerance") ;
      }

      // Set the default value.
      virtual void compute(const sequence& seq) {
        *cvodeTolerance=1.0e-03 ;
      }
  } ;

  register_rule<DefaultCVODETolerance> registerDefaultCVODETolerance ;

  // Sets the default eddy viscosity limit.
  class DefaultEddyViscosityLimit : public default_rule {
    private:
      param<real> eddyViscosityLimit ;
    public:

      // Define input and output.
      DefaultEddyViscosityLimit() {
        name_store("eddyViscosityLimit",eddyViscosityLimit) ;
        output("eddyViscosityLimit") ;
        comments("Sets the default eddy viscosity limit to 5.0e+05.") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) {
        *eddyViscosityLimit=5.0e+05 ;
      }
  } ;

  register_rule<DefaultEddyViscosityLimit> registerDefaultEddyViscosityLimit ;

  // Sets the default flow compressibility.
  class DefaultFlowCompressibility : public default_rule {
    private:
      param<string> flowCompressibility ;
    public:

      // Define input and output.
      DefaultFlowCompressibility() {
        name_store("flowCompressibility",flowCompressibility) ;
        output("flowCompressibility") ;
        comments("Sets the default flow compressibility to 'compressible'.") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) { *flowCompressibility="compressible" ; }
  } ;

  register_rule<DefaultFlowCompressibility> registerDefaultFlowCompressibility ;

  // Sets the default flow regime.
  class DefaultFlowRegime : public default_rule {
    private:
      param<string> flowRegime ;
    public:

      // Define input and output.
      DefaultFlowRegime() {
        name_store("flowRegime",flowRegime) ;
        output("flowRegime") ;
        comments("Sets the default flow regime to 'turbulent'.") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) { *flowRegime="turbulent" ; }
  } ;

  register_rule<DefaultFlowRegime> registerDefaultFlowRegime ;

  // Sets the default flow type.
  class DefaultFlowType : public default_rule {
    private:
      param<string> flowType ;
    public:

      // Define input and output.
      DefaultFlowType() {
        name_store("flowType",flowType) ;
        output("flowType") ;
        comments("Sets the default flow type to 'viscous'.") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) { *flowType="viscous" ; }
  } ;

  register_rule<DefaultFlowType> registerDefaultFlowType ;

  // Sets the default freezing temperature.
  class DefaultFreezingTemperature : public default_rule {
    private:
      param<double> Tf ;
    public:

      // Define input and output.
      DefaultFreezingTemperature() {
        name_store("Tf",Tf) ;
        output("Tf") ;
        comments("Freezing Temperature:  Below this temperature all ") ;
        comments("reactions cease to be evaluated.") ;
      }

      // Set the default value.
      virtual void compute(const sequence& seq) { *Tf = 0 ; }
  } ;

  register_rule<DefaultFreezingTemperature> registerDefaultFreezingTemperature ;

  // Default stencil size for scalable file interpolation.
  class DefaultInterpolateMinStencilSize : public default_rule {
    private:
      param<double> interpolateMinStencilSize ;
    public:

      // Define input and output.
      DefaultInterpolateMinStencilSize() {
        name_store("interpolateMinStencilSize",interpolateMinStencilSize) ;
        output("interpolateMinStencilSize") ;
        comments("Minimum distance that we should use in excluding points from interpolation stencils. The default is zero, in which case the estimate used from the number of points and volume of space are used. In some cases, this may be used to improve interpolation quality.") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {
        *interpolateMinStencilSize=0.0 ;
      }
  } ;

  register_rule<DefaultInterpolateMinStencilSize>
    registerDefaultInterpolateMinStencilSize ;

  // Sets the default inviscid flux.
  class DefaultInviscidFlux : public default_rule {
    private:
      param<string> inviscidFlux ;
    public:

      // Define input and output.
      DefaultInviscidFlux() {
        name_store("inviscidFlux",inviscidFlux) ;
        output("inviscidFlux") ;
        comments("Sets the default inviscidFlux to 'SOU'.") ;
        comments("Valid options are 'FOU','SOU','CD' and 'Roe'.") ;
      }

      // Set the default value.
      virtual void compute(const sequence& seq) { *inviscidFlux="SOU" ; }
  } ;

  register_rule<DefaultInviscidFlux> registerDefaultInviscidFlux ;

  // Sets the default matrix form, which determines whether or not we
  // make the matrix for the transport equations diagonally dominant
  // or not.
  class DefaultMatrixForm : public default_rule {
    private:
      param<string> matrixForm ;
    public:

      // Define input and output.
      DefaultMatrixForm() {
        name_store("matrixForm",matrixForm) ;
        output("matrixForm") ;
        comments("Sets the default matrix form to 'natural'.") ;
        comments("Can override with 'diagonalDominance'.") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) { *matrixForm="natural" ; }
  } ;

  register_rule<DefaultMatrixForm> registerDefaultMatrixForm ;

  // Sets the default maximum number of iterations per time-step.
  class DefaultMaxIterationsPerTimeStep : public default_rule {
    private:
      param<int> maxIterationsPerTimeStep ;
    public:

      // Define input and output.
      DefaultMaxIterationsPerTimeStep() {
        name_store("maxIterationsPerTimeStep",maxIterationsPerTimeStep) ;
        output("maxIterationsPerTimeStep") ;
      }

      // Set the default value.
      virtual void compute(const sequence& seq) {
        *maxIterationsPerTimeStep=1000 ;
      }
  } ;

  register_rule<DefaultMaxIterationsPerTimeStep>
    registerMaxIterationsPerTimeStep ;

  // Sets the default maximum temperature change (in percent) per time-step.
  class DefaultMaxTemperatureChange : public default_rule {
    private:
      param<real> maxTemperatureChange ;
    public:

      // Define input and output.
      DefaultMaxTemperatureChange() {
        name_store("maxTemperatureChange",maxTemperatureChange) ;
        output("maxTemperatureChange") ;
        comments("Sets the default maximum temperature change to -1.0.") ;
        comments("Any negative value turns off this feature.") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) { *maxTemperatureChange=-1.0 ; }
  } ;

  register_rule<DefaultMaxTemperatureChange>
    registerDefaultMaxTemperatureChange ;

  // Sets the default number of time-steps to one.
  class DefaultNumTimeSteps : public default_rule {
    private:
      param<int> numTimeSteps ;
    public:

      // Define input and output.
      DefaultNumTimeSteps() {
        name_store("numTimeSteps",numTimeSteps) ;
        output("numTimeSteps") ;
      }

      // Set the default value.
      virtual void compute(const sequence& seq) { *numTimeSteps=1 ; }
  } ;

  register_rule<DefaultNumTimeSteps> registerDefaultNumTimeSteps ;

  // Sets the default form of operator splitting.
  class DefaultOperatorSplitting : public default_rule {
    private:
      param<string> operatorSplitting ;
    public:

      // Define input and output.
      DefaultOperatorSplitting() {
        name_store("operatorSplitting",operatorSplitting) ;
        output("operatorSplitting") ;
        comments("Sets the default operator splitting to 'source'.") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) {
        *operatorSplitting="source" ;
      }
  } ;

  register_rule<DefaultOperatorSplitting> registerDefaultOperatorSplitting ;

  // Sets the default operator splitting EOS inversion value, which determines
  // when the EOS is inverted during the solution of the ODE.
  class DefaultOperatorSplittingEOS : public default_rule {
    private:
      param<int> operatorSplittingEOS ;
    public:

      // Define input and output.
      DefaultOperatorSplittingEOS() {
        name_store("operatorSplittingEOS",operatorSplittingEOS) ;
        output("operatorSplittingEOS") ;
        comments("Value to determine when EOS is inverted during ODE.") ;
        comments("solution. 1-never,2-each internal ODE timestep,3-") ;
        comments("each call to cvodeRHS") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) { *operatorSplittingEOS=1 ; }
  } ;

  register_rule<DefaultOperatorSplittingEOS>
    registerDefaultOperatorSplittingEOS ;

  // Sets the default convergence tolerance for the PETSC solver.
  class DefaultPetscConvergenceTolerance : public default_rule {
    private:
      param<real> petscConvergenceTolerance ;
    public:

      // Define input and output.
      DefaultPetscConvergenceTolerance() {
        name_store("petscConvergenceTolerance",petscConvergenceTolerance) ;
        output("petscConvergenceTolerance") ;
      }

      // Set the default value.
      virtual void compute(const sequence& seq) {
        *petscConvergenceTolerance=1.0e-03 ;
      }
  } ;

  register_rule<DefaultPetscConvergenceTolerance>
    registerDefaultPetscConvergenceTolerance ;

  // Sets the default plot modulo.
  class DefaultPlotModulo : public default_rule {
    private:
      param<int> plot_modulo ;
    public:

      // Define input and output.
      DefaultPlotModulo() {
        name_store("plot_modulo",plot_modulo) ;
        output("plot_modulo") ;
        comments("Plot files are output indexed by the iteration number modulo this parameter. ") ;
        comments("A zero indicates that the files will be indexed only by iteration number.") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) { *plot_modulo=0 ; }
  } ;

  register_rule<DefaultPlotModulo> registerDefaultPlotModulo ;

  // Sets the default print frequency.
  class DefaultPrintFreq: public default_rule {
    private:
      param<int> print_freq ;
    public:

      // Define input and output.
      DefaultPrintFreq() {
        name_store("print_freq",print_freq) ;
        output("print_freq") ;
        comments("Print integrated values and CFL conditions whenever the timestep is divisible by this input") ;
      }

      // Set the default value.
      virtual void compute(const sequence& seq) { *print_freq=100 ; }
  } ;

  register_rule<DefaultPrintFreq> registerDefaultPrintFreq ;

  // Scale factors for normalizing residuals.
  class DefaultReferenceValue : public default_rule {
    private:
      param<ReferenceValue> referenceValue ;
    public:

      // Define input and output.
      DefaultReferenceValue() {
        name_store("referenceValue",referenceValue) ;
        output("referenceValue") ;
        comments("Used to set reference values for residual normalization.") ;
      }

      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;

  register_rule<DefaultReferenceValue> registerDefaultReferenceValue ;

  // Sets the default restart modulo.
  class DefaultRestartModulo : public default_rule {
    private:
      param<int> restart_modulo ;
    public:

      // Define input and output.
      DefaultRestartModulo() {
        name_store("restart_modulo",restart_modulo) ;
        output("restart_modulo") ;
        comments("Restart files are output indexed by the iteration number ") ;
        comments("modulo this parameter. A zero indicates that the files ") ;
        comments("will be indexed only by iteration number.") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) { *restart_modulo=0 ; }
  } ;

  register_rule<DefaultRestartModulo> registerDefaultRestartModulo ;

  // Sets the default species production type.
  class DefaultSpeciesProduction : public default_rule {
    private:
      param<string> speciesProduction ;
    public:

      // Define input and output.
      DefaultSpeciesProduction() {
        name_store("speciesProduction",speciesProduction) ;
        output("speciesProduction") ;
        comments("Sets the default species production type.") ;
        comments("Valid values are 'CVODE' and 'instantaneous'.") ;
      }

      // Set the default model.
      virtual void compute(const sequence &seq) { *speciesProduction="CVODE" ; }
  } ;

  register_rule<DefaultSpeciesProduction> registerDefaultSpeciesProduction ;

  // Sets the default testing value. This is used to set up a constraint to
  // enable case-specific rules in testing.cc .
  class DefaultTesting : public default_rule {
    private:
      param<string> testing ;
    public:

      // Define input and output.
      DefaultTesting() {
        name_store("testing",testing) ;
        output("testing") ;
        comments("Sets the default testing value to 'none'. This name is ") ;
        comments("used to active case-specific rules found in testing.cc .") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) { *testing="none" ; }
  } ;

  register_rule<DefaultTesting> registerDefaultTesting ;

  // Sets the default time integrator.
  class DefaultTimeIntegrator : public default_rule {
    private:
      param<string> timeIntegrator ;
    public:

      // Define input and output.
      DefaultTimeIntegrator() {
        name_store("timeIntegrator",timeIntegrator) ;
        output("timeIntegrator") ;
        comments("Sets the default time integrator to 'BDF'. Valid options") ;
        comments(" are 'BDF','BDF2' and 'CN'.") ;
      }

      // Set the default value.
      virtual void compute(const sequence &seq) { *timeIntegrator="BDF" ; }
  } ;

  register_rule<DefaultTimeIntegrator> registerDefaultTimeIntegrator ;

  // Sets the default time-step to a large value thus giving a steady-state
  // run.
  class DefaultTimeStep : public default_rule {
    private:
      param<real> timeStep ;
    public:

      // Define input and output.
      DefaultTimeStep() {
        name_store("timeStep",timeStep) ;
        output("timeStep") ;
      }

      // Set the default value.
      virtual void compute(const sequence& seq) { *timeStep=1.0e30 ; }
  } ;

  register_rule<DefaultTimeStep> registerDefaultTimeStep ;

  // Sets default value for freestream turbulence intensity.
  class DefaultTurbulenceIntensityFreestream : public default_rule {
    private:
      param<real> turbulenceIntensityFreestream ;
    public:

      // Define input and output.
      DefaultTurbulenceIntensityFreestream() {
        name_store("turbulenceIntensityFreestream",
          turbulenceIntensityFreestream) ;
        output("turbulenceIntensityFreestream") ;
      }

      // Set the default value.
      virtual void compute(const sequence& seq) {
        *turbulenceIntensityFreestream=0.02 ;
      }
  } ;

  register_rule<DefaultTurbulenceIntensityFreestream>
    registerDefaultTurbulenceIntensityFreestream ;

  // Sets the default turbulence inviscid flux.
  class DefaultTurbulenceInviscidFlux : public default_rule {
    private:
      param<string> turbulenceInviscidFlux ;
    public:

      // Define input and output.
      DefaultTurbulenceInviscidFlux() {
        name_store("turbulenceInviscidFlux",turbulenceInviscidFlux) ;
        output("turbulenceInviscidFlux") ;
        comments("Sets the default turbulenceInviscidFlux to 'SOU'.") ;
        comments("Valid options are 'FOU' and 'SOU'.") ;
      }

      // Set the default value.
      virtual void compute(const sequence& seq) {
        *turbulenceInviscidFlux="SOU" ;
      }
  } ;

  register_rule<DefaultTurbulenceInviscidFlux>
    registerDefaultTurbulenceInviscidFlux ;

  // Sets the default turbulent Prandtl number.
  class DefaultTurbulentPrandtlNumber : public default_rule {
    private:
      param<real> turbulentPrandtlNumber ;
    public:

      // Define input and output.
      DefaultTurbulentPrandtlNumber() {
        name_store("turbulentPrandtlNumber",turbulentPrandtlNumber) ;
        output("turbulentPrandtlNumber") ;
      }

      // Set the default value.
      virtual void compute(const sequence& seq) {
        *turbulentPrandtlNumber=0.7 ;
      }
  } ;

  register_rule<DefaultTurbulentPrandtlNumber>
    registerDefaultTurbulentPrandtlNumber ;

  // Sets the default turbulent Schmidt number.
  class DefaultTurbulentSchmidtNumber : public default_rule {
    private:
      param<real> turbulentSchmidtNumber ;
    public:

      // Define input and output.
      DefaultTurbulentSchmidtNumber() {
        name_store("turbulentSchmidtNumber",turbulentSchmidtNumber) ;
        output("turbulentSchmidtNumber") ;
      }

      // Set the default value.
      virtual void compute(const sequence& seq) {
        *turbulentSchmidtNumber=0.95 ;
      }
  } ;

  register_rule<DefaultTurbulentSchmidtNumber>
    registerDefaultTurbulentSchmidtNumber ;

  // Sets the default midpoint for 'ql' and 'qr' to 0.0.
  class DefaultXMid : public default_rule {
    private:
      param<real> xmid ;
    public:

      // Define input and output.
      DefaultXMid() {
        name_store("xmid",xmid) ;
        output("xmid") ;
      }

      // Set the default value.
      virtual void compute(const sequence& seq) { *xmid=0.0 ; }
  } ;

  register_rule<DefaultXMid> registerDefaultXMid ;

  // Options for the energy equation.
  class OptionalEnergyEquation : public optional_rule {
    private:
      param<EnergyEquationOptions> energyEquationOptions ;
    public:

      // Define input and output.
      OptionalEnergyEquation() {
        name_store("energyEquationOptions",energyEquationOptions) ;
        output("energyEquationOptions") ;
      }

      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;

  register_rule<OptionalEnergyEquation> registerOptionalEnergyEquation ;

  // Gravity options.
  class OptionalGravity : public optional_rule {
    private:
      param<GravityOptions> gravity ;
    public:

      // Define input and output.
      OptionalGravity() {
        name_store("gravity",gravity) ;
        output("gravity") ;
      }

      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;

  register_rule<OptionalGravity> registerOptionalGravity ;

  // Ignition options.
  class OptionalIgnition : public optional_rule {
    private:
      param<IgnitionOptions> ignition ;
    public:

      // Define input and output.
      OptionalIgnition() {
        name_store("ignition",ignition) ;
        output("ignition") ;
      }

      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;

  register_rule<OptionalIgnition> registerOptionalIgnition ;

  // Initial condition to set flowfield values for the entire domain.
  class OptionalInitialCondition : public optional_rule {
    private:
      param<InitialCondition> initialCondition ;
    public:

      // Define input and output.
      OptionalInitialCondition() {
        name_store("initialCondition",initialCondition) ;
        output("initialCondition") ;
      }

      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;

  register_rule<OptionalInitialCondition> registerOptionalInitialCondition ;

  // Optional filename for scalable binary initial conditions from file.
  class OptionalInterpolateInitialConditions : public optional_rule {
    private:
      param<string> interpolateInitialConditions ;
    public:

      // Define input and output.
      OptionalInterpolateInitialConditions() {
        name_store("interpolateInitialConditions",interpolateInitialConditions);
        output("interpolateInitialConditions") ;
        comments("Scalable binary initial conditions file. Set this variable to the name of the 'put' file that you want to use to interpolate the initial conditions. This file is typically generated from another run of CHEM or STREAM.") ;
      }

      // Do nothing.
      virtual void compute(const sequence &seq) {}
  } ;

  register_rule<OptionalInterpolateInitialConditions>
    registerOptionalInterpolateInitialConditions ;

  // Optional freestream k value.
  class OptionalKFreestream : public optional_rule {
    private:
      param<real> kFreestream ;
    public:

      // Define input and output.
      OptionalKFreestream() {
        name_store("kFreestream",kFreestream) ;
        output("kFreestream") ;
      }

      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;

  register_rule<OptionalKFreestream> registerOptionalKFreestream ;

  // Options for the momentum equation .
  class OptionalMomentumEquation : public optional_rule {
    private:
      param<MomentumEquationOptions> momentumEquationOptions ;
    public:

      // Define input and output.
      OptionalMomentumEquation() {
        name_store("momentumEquationOptions",momentumEquationOptions) ;
        output("momentumEquationOptions") ;
      }

      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;

  register_rule<OptionalMomentumEquation> registerOptionalMomentumEquation ;

  // Optional freestream omega value.
  class OptionalOmegaFreestream : public optional_rule {
    private:
      param<real> omegaFreestream ;
    public:

      // Define input and output.
      OptionalOmegaFreestream() {
        name_store("omegaFreestream",omegaFreestream) ;
        output("omegaFreestream") ;
      }

      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;

  register_rule<OptionalOmegaFreestream> registerOptionalOmegaFreestream ;

  // Perturbation options.
  class OptionalPerturbation : public optional_rule {
    private:
      param<PerturbationOptions> perturbation ;
    public:

      // Define input and output.
      OptionalPerturbation() {
        name_store("perturbation",perturbation) ;
        output("perturbation") ;
      }

      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;

  register_rule<OptionalPerturbation> registerOptionalPerturbation ;

  // Plot frequency.
  class OptionalPlotFreq : public optional_rule {
    private:
      param<int> plot_freq ;
    public:
                                                                                
      // Define input and output.
      OptionalPlotFreq() {
        name_store("plot_freq",plot_freq) ;
        output("plot_freq") ;
        comments("Output plot files whenever the timestep number is divisible by this input.") ;
      }
                                                                                
      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;
                                                                                
  register_rule<OptionalPlotFreq> registerOptionalPlotFreq ;

  // Options for the pressure equation .
  class OptionalPressureEquation : public optional_rule {
    private:
      param<PressureEquationOptions> pressureEquationOptions ;
    public:

      // Define input and output.
      OptionalPressureEquation() {
        name_store("pressureEquationOptions",pressureEquationOptions) ;
        output("pressureEquationOptions") ;
      }

      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;

  register_rule<OptionalPressureEquation> registerOptionalPressureEquation ;

  // Flowfield values for the left half of the domain as determined by xmid.
  class OptionalQL : public optional_rule {
    private:
      param<InitialCondition> qL ;
    public:

      // Define input and output.
      OptionalQL() { name_store("ql",qL) ; output("ql") ; }

      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;

  register_rule<OptionalQL> registerOptionalQL ;

  // Flowfield values for the right half of the domain as determined by xmid.
  class OptionalQR : public optional_rule {
    private:
      param<InitialCondition> qR ;
    public:

      // Define input and output.
      OptionalQR() { name_store("qr",qR) ; output("qr") ; }

      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;

  register_rule<OptionalQR> registerOptionalQR ;

  // File for specifying reference frames.
  class OptionalReferenceFrameFile : public optional_rule {
    private:
      param<string> referenceFrameFile ;
    public:
                                                                                
      // Define input and output.
      OptionalReferenceFrameFile() {
        name_store("referenceFrameFile",referenceFrameFile) ;
        output("referenceFrameFile") ;
        comments("File for specifying reference frames.") ;
      }
                                                                                
      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;
                                                                                
  register_rule<OptionalReferenceFrameFile>
    registerOptionalReferenceFrameFile ;

  // Restart frequency.
  class OptionalRestartFreq : public optional_rule {
    private:
      param<int> restart_freq ;
    public:
                                                                                
      // Define input and output.
      OptionalRestartFreq() {
        name_store("restart_freq",restart_freq) ;
        output("restart_freq") ;
        comments("Output restart file when timestep is divisible by this value") ;
      }
                                                                                
      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;
                                                                                
  register_rule<OptionalRestartFreq> registerOptionalRestartFreq ;

  // Options for the species equations.
  class OptionalSpeciesEquation : public optional_rule {
    private:
      param<SpeciesEquationOptions> speciesEquationOptions ;
    public:

      // Define input and output.
      OptionalSpeciesEquation() {
        name_store("speciesEquationOptions",speciesEquationOptions) ;
        output("speciesEquationOptions") ;
      }

      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;

  register_rule<OptionalSpeciesEquation> registerOptionalSpeciesEquation ;

  // Options for the turbulence equations.
  class OptionalTurbulenceEquation : public optional_rule {
    private:
      param<TurbulenceEquationOptions> turbulenceEquationOptions ;
    public:

      // Define input and output.
      OptionalTurbulenceEquation() {
        name_store("turbulenceEquationOptions",turbulenceEquationOptions) ;
        output("turbulenceEquationOptions") ;
      }

      // Do nothing.
      virtual void compute(const sequence& seq) {}
  } ;

  register_rule<OptionalTurbulenceEquation> registerOptionalTurbulenceEquation ;

}
