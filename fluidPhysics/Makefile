include ../Stream.conf
include $(LOCI_BASE)/Loci.conf

INCLUDES = -I.
#COPT=-O0 -g #-DDEBUG -DBOUNDS_CHECK

# Put objects in the module here
OBJS = PerfectGas.o chemistry_db.o periodic_table.o reaction.o \
	transport_db.o Matrix.o property.o

LIB_NAME = fluidPhysicsSTREAM

LOCAL_LIBS = 

JUNK = *~  core ti_files ii_files rii_files *.d

LIB_OBJS=$(OBJS:.o=_lo.o)

all: lib$(LIB_NAME).so 
	mkdir -p $(STREAM_INSTALL_DIR)/lib
	mv lib$(LIB_NAME).so $(STREAM_INSTALL_DIR)/lib

lib$(LIB_NAME).so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) lib$(LIB_NAME).so $(LIB_FLAGS) $(LIB_OBJS)
FRC : 



clean:
	rm -fr $(OBJS) $(LIB_OBJS) lib$(LIB_NAME).so $(JUNK)

LOCI_FILES = $(wildcard *.loci)
LOCI_LPP_FILES = $(LOCI_FILES:.loci=.cc)

distclean: 
	rm $(DEPEND_FILES)
	rm -fr $(OBJS) $(LIB_OBJS) $(DROPLET_OBJS_LO) droplet $(MODULE_NAME)_m.so $(JUNK) $(LOCI_LPP_FILES)



DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
-include $(DEPEND_FILES)
