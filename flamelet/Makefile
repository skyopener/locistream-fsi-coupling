include ../Stream.conf
include $(LOCI_BASE)/Loci.conf

INCLUDES = -I$(STREAM_BASE)/src -I$(STREAM_BASE)/src/gridReader -I$(STREAM_BASE)/fluidPhysics

# Put objects in the module here
OBJS =	bcCheck.o bcConstraint.o bcInput.o chi.o flamelet.o flamelet_table.o gsi.o \
  initialCondition.o residual.o scalarOutput.o table_read.o flameletVarsFileInputs.o \
  restart.o 

MODULE_NAME = flamelet

LOCAL_LIBS =

JUNK = *~  core ti_files ii_files rii_files *.d

LOCI_FILES = $(wildcard *.loci)
LOCI_LPP_FILES = $(LOCI_FILES:.loci=.cc)

LIB_OBJS=$(OBJS:.o=_lo.o)
$(MODULE_NAME)_m.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) $(MODULE_NAME)_m.so $(LIB_FLAGS) $(LIB_OBJS)
	mv $(MODULE_NAME)_m.so $(STREAM_INSTALL_DIR)/lib
 
FRC : 

all: $(MODULE_NAME)_m.so
 
clean:
	rm -f $(DEPEND_FILES)
	rm -fr $(OBJS) $(LIB_OBJS) $(MODULE_NAME)_m.so $(JUNK) $(LOCI_LPP_FILES)
 
distclean:
	rm -f $(DEPEND_FILES)
	rm -fr $(OBJS) $(LIB_OBJS) $(MODULE_NAME)_m.so $(JUNK) $(LOCI_LPP_FILES)
 
DEPEND_FILES=$(subst .o,.d,$(OBJS))

#include automatically generated dependencies
include $(DEPEND_FILES)
