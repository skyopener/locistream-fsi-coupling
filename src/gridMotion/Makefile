
LOCI_BASE=/nobackup/kangck/codes/loci/Loci-3.2-beta-12/OBJ
include $(LOCI_BASE)/Loci.conf


TARGET = gridMotion_m.so

#OBJS =	boundaryCondition.o gridIO.o hypre.o metrics.o move.o \
  petsc.o restart.o scalarOutput.o solver.o
OBJS =	boundaryCondition.o deformation.o flux.o hypre.o metrics.o \
  move.o nodeGradient.o nodeValues.o petsc.o piso.o restart.o scalarOutput.o \
  solver.o source.o

JUNK = *~  core ti_files ii_files rii_files *.d

LOCI_FILES = $(wildcard *.loci)
LOCI_LPP_FILES = $(LOCI_FILES:.loci=.cc)

INCLUDES = $(CC_LIB_FLAGS) -I. -I.. $(MAIN_INCLUDES)

default:
	$(MAKE) $(TARGET)

LIB_OBJS=$(OBJS:.o=_lo.o)
gridMotion_m.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) gridMotion_m.so $(LIB_FLAGS) \
	$(LIB_OBJS)

FRC : 

allclean: clean
	rm -f make.depend make.depend.bak

clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(MAIN_OBJS) $(TARGET) $(JUNK) $(LOCI_LPP_FILES)

distclean: clean
	rm -fr $(DEPEND_FILES)

ALL_OFILES = $(OBJS)

DEPEND_FILES=$(subst .o,.d,$(ALL_OFILES))


#include automatically generated dependencies
include $(DEPEND_FILES)
