include $(LOCI_BASE)/Loci.conf


TARGET = GridReader.a  #GridReader.so 

OBJS =  readGrid.o gridTools.o setupBC.o 


JUNK = *~ core ti_files ii_files *.d


INCLUDES = $(CC_LIB_FLAGS) -I.. $(MAIN_INCLUDES)

default:
	$(MAKE) $(TARGET)

$(TARGET): $(OBJS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)

readGrid.o: readGrid.cc
	$(CP) $(EXCEPTIONS) $(DEFINES) $(LOCI_INCLUDES) $(INCLUDES) -c readGrid.cc

GridReader.so: $(OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) GridReader.so $(LIB_FLAGS) $(OBJS)


FRC : 

allclean: clean
	rm -f make.depend make.depend.bak

clean:
	rm -fr $(OBJS) $(TARGET) $(JUNK)

make.depend: 
	touch make.depend 

depend: make.depend
	$(MAKEDEPEND) $(DEFINES) $(INCLUDES) $(LOCI_INCLUDES) *.cc > make.depend

include make.depend				

# DO NOT DELETE
