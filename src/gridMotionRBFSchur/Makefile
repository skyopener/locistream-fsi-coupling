LOCI_BASE=/nobackup/kangck/codes/loci/Loci-3.2-beta-12/OBJ
include $(LOCI_BASE)/Loci.conf

# ".cc" files created from .loci files
LOCI_INTERMEDIATE_FILES = $(subst .loci,.cc, $(wildcard *.loci) )

TARGET = gridMotionRBFSchur_m.so

#OBJS =	boundaryCondition.o gridIO.o hypre.o metrics.o move.o \
  petsc.o restart.o scalarOutput.o solver.o
OBJS =	boundaryCondition.o RBF_deformation.o flux.o gridIO.o metrics.o \
  RBF_move.o nodeGradient.o nodeValues.o RBF_petsc.o restart.o scalarOutput.o \
  solver.o source.o deformation.o

#JUNK = *~  core ti_files ii_files rii_files *.d
# Junk files that are created while editing and running cases
JUNK = $(wildcard *~) $(wildcard crash_dump.*) core debug output $(wildcard .schedule*) ti_files ii_files rii_files *.d

INCLUDES = $(CC_LIB_FLAGS) -I. -I.. $(MAIN_INCLUDES)

default:
	$(MAKE) $(TARGET)

LIB_OBJS=$(OBJS:.o=_lo.o)
gridMotionRBFSchur_m.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) gridMotionRBFSchur_m.so $(LIB_FLAGS) \
	$(LIB_OBJS)

FRC : 

allclean: clean
	rm -f make.depend make.depend.bak

clean:
#	rm -fr $(OBJS) $(LIB_OBJS) $(MAIN_OBJS) $(TARGET) $(JUNK)
	rm -fr $(OBJS) $(LIB_OBJS) $(MAIN_OBJS) $(TARGET) $(JUNK) $(LOCI_INTERMEDIATE_FILES)

distclean: clean
	rm -fr $(DEPEND_FILES) $(JUNK) $(LOCI_INTERMEDIATE_FILES)

ALL_OFILES = $(OBJS)

DEPEND_FILES=$(subst .o,.d,$(ALL_OFILES))


#include automatically generated dependencies
include $(DEPEND_FILES)
