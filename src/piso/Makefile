include ../../Stream.conf
include $(LOCI_BASE)/Loci.conf

SUNDIALSINCLUDES = -I$(SUNDIALSDIR)/include

HYPREINCLUDES = -I$(HYPREDIR)/include

MESHDEFORMATION = -DMESH_DEFORMATION
INCLUDES = -I.. -I$(STREAM_BASE)/fluidPhysics $(SUNDIALSINCLUDES) \
$(HYPREINCLUDES) $(MESHDEFORMATION)

TARGET = streamUns

STREAM_OBJS = ../bcCheck.o ../bcConstraint.o ../bcInput.o ../bcReconstruct.o \
  ../chem.o ../constraints.o ../eos_interface.o ../fluidState.o ../hypre.o \
  ../icInterp.o ../initialCondition.o ../initialConditionRegions.o \
  ../localTimeStepping.o ../metrics.o \
  ../misc.o ../petsc.o ../read_tran.o ../referenceFrame.o ../residual.o \
  ../sciTypes.o ../streamUns.o ../transport_coef.o ../varsFileInputs.o \
  ../wallFunction.o

OBJS = crankNicholson.o density.o energy.o equationOfState.o flux.o gridIO.o \
  gsi.o integrate.o kOmegaMenter.o momentum.o perturbation.o pressure.o \
  pressureCorrection.o probe.o restart.o scalarOutput.o species.o testing.o

STREAM_LIBS = $(STREAM_BASE)/src/gridReader/GridReader.a -L$(STREAM_INSTALL_DIR)/lib \
-lfluidPhysicsSTREAM

SUNDIALSLIBS = $(SUNDIALSDIR)/lib/libsundials_cvode.a \
  $(SUNDIALSDIR)/lib/libsundials_nvecserial.a

HYPRELIBS = -L$(HYPREDIR)/lib -lHYPRE -lHYPRE_IJ_mv -lHYPRE_parcsr_ls \
  -lHYPRE_parcsr_mv -lHYPRE_utilities -lHYPRE_LSI

JUNK = *~ core ti_files ii_files rii_files *.d

default:
	$(MAKE) $(TARGET)

$(TARGET): $(STREAM_OBJS) $(OBJS) FRC 
	$(LD) -o streamUns -Wl,-rpath $(STREAM_INSTALL_DIR)/lib $(OBJS) \
        $(STREAM_OBJS) $(STREAM_LIBS) $(SUNDIALSLIBS) $(HYPRELIBS) $(LIBS) \
        $(LDFLAGS) -Wl,-rpath $(CHEMMODULEDIR)

streamUns.o: streamUns.cc
	$(CP) $(COPTLESS) $(EXCEPTIONS) $(DEFINES) $(LOCI_INCLUDES) \
	$(INCLUDES) -c streamUns.cc

FRC : 

allclean: clean
	rm -f make.depend make.depend.bak

clean:
	rm -fr $(OBJS) $(JUNK)

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
