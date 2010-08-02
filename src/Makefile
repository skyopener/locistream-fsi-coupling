include ../Stream.conf
include $(LOCI_BASE)/Loci.conf

SUNDIALSINCLUDES = -I$(SUNDIALSDIR)/include

HYPREINCLUDES = -I$(HYPREDIR)/include

# Start developer defines.
ODE_PERFORMANCE =#-DODE_PERFORMANCE # timing of ode integration
PETSC_PERFORMANCE = -DPETSC_USE_LOG # PETSC performance logging
DEVELOPER = $(ODE_PERFORMANCE) $(PETSC_PERFORMANCE)
# End developer defines.

MESHDEFORMATION = -DMESH_DEFORMATION
INCLUDES = -I$(STREAM_BASE)/src -I$(STREAM_BASE)/src/gridReader -I$(STREAM_BASE)/fluidPhysics \
$(SUNDIALSINCLUDES) $(HYPREINCLUDES) $(MESHDEFORMATION) $(DEVELOPER)

STREAM_EXEC = stream
TARGET = streamUns

OBJS = axisymmetric.o bcCheck.o bcConstraint.o bcInput.o bcReconstruct.o \
  chem.o constraints.o crankNicholson.o date.o density.o des.o energy.o \
  eos_interface.o equationOfState.o fluidState.o flux.o gridIO.o gsi.o \
  hypre.o icInterp.o ignition.o initialCondition.o initialConditionRegions.o \
  interpolateFile.o integrate.o kOmegaMenter.o \
  localTimeStepping.o metrics.o misc.o momentum.o perturbation.o periodic.o \
  petsc.o pressure.o pressureCorrection.o probe.o \
  properties.o read_tran.o referenceFrame.o residual.o restart.o roe.o \
  scalarOutput.o sciTypes.o species.o ssc.o streamUns.o testCVODE.o testing.o \
  timeStep.o tmp.o transport_coef.o varsFileInputs.o wallFunction.o wallOutput.o

LOCAL_LIBS = gridReader/GridReader.a -L$(STREAM_INSTALL_DIR)/lib -lfluidPhysicsSTREAM

MODULE_LIBS = $(STREAM_BASE)/src/gridMotion/gridMotion_m.so gridReader/GridReader.a \
  $(STREAM_BASE)/src/turbo/turbo_m.so \
  $(STREAM_BASE)/src/FSIEXPLICIT/FSIEXPLICIT_m.so

SUNDIALSLIBS = $(SUNDIALSDIR)/lib/libsundials_cvode.a \
  $(SUNDIALSDIR)/lib/libsundials_nvecserial.a

#HYPRELIBS = -L$(HYPREDIR)/lib -lHYPRE -lHYPRE_IJ_mv -lHYPRE_parcsr_ls \
# -lHYPRE_parcsr_mv -lHYPRE_utilities -lHYPRE_LSI
HYPRELIBS = -L$(HYPREDIR)/lib -lHYPRE

JUNK = *~ core ti_files ii_files rii_files *.d

LOCI_FILES = $(wildcard *.loci)
LOCI_LPP_FILES = $(LOCI_FILES:.loci=.cc)

install: default
	STREAM_INSTALL_DIR=$(STREAM_INSTALL_DIR) STREAM_EXEC=$(STREAM_EXEC) \
	bash Install.bash

default:
	$(MAKE) $(TARGET) ;
	cd piso && ( $(MAKE) ; cd ..)

$(TARGET): $(OBJS) $(MODULE_LIBS) FRC 
	mkdir -p $(STREAM_INSTALL_DIR)/lib
	$(LD) -o streamUns -Wl,-rpath $(STREAM_INSTALL_DIR)/lib $(OBJS) \
        $(LOCAL_LIBS) $(SUNDIALSLIBS) $(HYPRELIBS) $(LIBS) $(LDFLAGS) \
	-Wl,-rpath $(CHEMMODULEDIR)

$(STREAM_BASE)/src/turbo/turbo_m.so: FRC
	cd turbo ; $(MAKE) turbo_m.so LOCI_BASE="$(LOCI_BASE)" \
	MAIN_INCLUDES="$(INCLUDES)"

$(STREAM_BASE)/src/gridMotion/gridMotion_m.so: FRC
	cd gridMotion ; $(MAKE) gridMotion_m.so LOCI_BASE="$(LOCI_BASE)" \
	MAIN_INCLUDES="$(INCLUDES)"

$(STREAM_BASE)/src/FSIEXPLICIT/FSIEXPLICIT_m.so: FRC
	cd FSIEXPLICIT ; $(MAKE) FSIEXPLICIT_m.so LOCI_BASE="$(LOCI_BASE)" \
	MAIN_INCLUDES="$(INCLUDES)"

gridReader/GridReader.a: FRC
	cd gridReader ; $(MAKE) GridReader.a LOCI_BASE="$(LOCI_BASE)" \
	MAIN_INCLUDES="$(INCLUDES)"

streamUns.o: streamUns.cc
	$(CP) $(COPTLESS) $(EXCEPTIONS) $(DEFINES) $(LOCI_INCLUDES) \
	$(INCLUDES) -c streamUns.cc

FRC : 

allclean: clean
	rm -f make.depend make.depend.bak

clean:
	rm -fr $(OBJS) $(TARGET) $(JUNK) $(LOCI_LPP_FILES)
	cd gridMotion && ( $(MAKE) clean LOCI_BASE="$(LOCI_BASE)" ; cd ..)
	cd gridReader && ( $(MAKE) clean LOCI_BASE="$(LOCI_BASE)" ; cd ..)
	cd turbo && ( $(MAKE) clean LOCI_BASE="$(LOCI_BASE)" ; cd ..)
	cd piso && ( $(MAKE) clean LOCI_BASE="$(LOCI_BASE)" ; cd ..)
	cd FSIEXPLICIT && ( $(MAKE) clean LOCI_BASE="$(LOCI_BASE)" ; cd ..)

ALL_OFILES = $(OBJS)
                                                                                
DEPEND_FILES=$(subst .o,.d,$(ALL_OFILES))
                                                                                
                                                                                
#include automatically generated dependencies
include $(DEPEND_FILES)

#make.depend:
#	touch make.depend 

#depend: make.depend
#	$(MAKEDEPEND) $(DEFINES) $(INCLUDES) $(LOCI_INCLUDES) *.cc > make.depend
#	cd gridReader && ( $(MAKE) depend LOCI_BASE="$(LOCI_BASE)" ; cd .. )

#include make.depend

# DO NOT DELETE
