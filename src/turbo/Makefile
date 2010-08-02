include $(LOCI_BASE)/Loci.conf

TARGET = turbo_m.so

OBJS =	boundaryCondition.o cellGradient.o energy.o momentum.o rotorStator.o

JUNK = *~  core ti_files ii_files rii_files *.d

INCLUDES = $(CC_LIB_FLAGS) -I. -I.. $(MAIN_INCLUDES)

default:
	$(MAKE) $(TARGET)

LIB_OBJS=$(OBJS:.o=_lo.o)
turbo_m.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) turbo_m.so $(LIB_FLAGS) \
	$(LIB_OBJS)

FRC : 

allclean: clean
	rm -f make.depend make.depend.bak

clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(MAIN_OBJS) $(TARGET) $(JUNK)

distclean: clean
	rm -fr $(DEPEND_FILES)

ALL_OFILES = $(OBJS)

DEPEND_FILES=$(subst .o,.d,$(ALL_OFILES))


#include automatically generated dependencies
include $(DEPEND_FILES)
