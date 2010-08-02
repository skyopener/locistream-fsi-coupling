# Load LOCI related configs
LOCI_BASE=/nobackup/kangck/codes/loci/Loci-3.2-beta-12/OBJ
include $(LOCI_BASE)/Loci.conf

#PETSCSHAREDLIBS=-L/home/kangck/codes/loci/petsc-2.3.3-p11/lib/linux-gnu-c-opt
#HYPRELIBS=-L/home/kangck/codes/loci/petsc-2.3.3-p11/externalpackages/hypre-2.0.0/linux-gnu-c-opt/lib/ -lHYPRE -lHYPRE_IJ_mv -lHYPRE_parcsr_ls -lHYPRE_parcsr_mv -lHYPRE_utilities -lHYPRE_LSI
#HYPREINCLUDES=/home/kangck/codes/loci/petsc-2.3.3-p11/externalpackages/hypre-2.0.0/linux-gnu-c-opt/include
#HYPREDIR=/home/esozer/Loci-Stream/hypre-1.10.0b
#HYPREINCLUDES = -I$(HYPREDIR)/include
#HYPRELIBS = -L$(HYPREDIR)/lib -lHYPRE -lHYPRE_IJ_mv -lHYPRE_parcsr_ls \
#  -lHYPRE_parcsr_mv -lHYPRE_utilities -lHYPRE_LSI

# Fortran compiler for NLAMS
F90=ifort
F90FLAGS = -fpp -warn -fPIC -check bounds -ftrapuv -zero -r8 -ftz -fp-model double

# EULERBEAM1D Object file
OBJSBEAM1D= beam1dmesh.o beam1dgauss_elimination.o beam1dstiff.o beam1dinterface.o

# Loci-STREAM FSI Object files
OBJS =	FSI_boundaryCondition.o FSI_move.o FSI_iterationFinishedPriority.o \
	FSI_CSDVarsInput.o FSI_CSD2CFD.o FSI_CSDdataloop.o \
	FSI_RemeshSetup.o FSI_RemeshInterpolation.o \
	FSI_CFD2CSD.o FSI_CFD2CSDInterpolation.o \
	deformation.o flux.o gridIO.o metrics.o nodeGradient.o nodeValues.o restart.o scalarOutput.o solver.o source.o \
	rbf.o

# NLAMS Object files
OBJSNLAMS= lib_kind.o lib_fem.o lib_sparse.o lib_rodrigues.o \
	nlams_transforms.o lib_tools.o lib_entry.o nlams_LociSTREAMdata.o \
	nlams_mod_shared.o nlams_input.o nlams_shared_time.o nlams_mesh_connect.o \
	nlams_loads.o nlams_init.o interface_lapack_v3.o nlams_elemlocal.o nlams_dkt_shape.o \
	nlams_opt_shape.o nlams_index.o nlams_abbd.o nlams_update.o nlams_shell.o nlams_undeformed.o \
	nlams_stiff.o nlams_shared_header.o nlams_header.o nlams_mass.o nlams_damp.o nlams_interf.o nlams_solver.o \
	nlams_solver_generalalfa.o nlams_NLAMS2LociSTREAM.o nlams_excuteNLAMS.o nlams_communicateLociSTREAM2NLAMS.o

# F90 compilation
%.o: %.f90
	$(F90) -c -o $@ $(F90FLAGS) $<

INCLUDES = $(CC_LIB_FLAGS) -I. -I.. $(MAIN_INCLUDES)
#INCLUDES = -I$(HYPREINCLUDES)

LDFLAGS= -L/usr/caen/intel-11.0/fc/11.0.074/lib/intel64 -lifport -lifcore -limf -Wabi -Wcast-align -ldl
#LAPACKLIBS= -L$(MKL_LINK) -lmkl -lguide -lpthread
LAPACKLIBS= -L$(MKL_LINK) -lmkl_solver_lp64_sequential -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lguide -lpthread
LIB_LOCAL_FLAGS= $(LDFLAGS) $(LAPACKLIBS) $(HYPRELIBS) $(PETSCSHAREDLIBS)

# Junk files that are created while editing and running cases
JUNK = $(wildcard *~) $(wildcard crash_dump.*) core debug output $(wildcard .schedule*) ti_files ii_files rii_files *.d
# ".cc" files created from .loci files
LOCI_INTERMEDIATE_FILES = $(subst .loci,.cc, $(wildcard *.loci) )


# Target library
TARGET = FSI_m.so

default:
	$(MAKE) $(TARGET)

ALL_OBJS = $(OBJS)
LIB_OBJS= $(OBJSNLAMS) $(OBJSBEAM1D) $(OBJS:.o=_lo.o)
#LIB_OBJS= $(OBJSNLAMS) $(OBJS:.o=_lo.o)
FSI_m.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) FSI_m.so $(LIB_FLAGS) $(LIB_LOCAL_FLAGS) \
	$(LIB_OBJS) $(INCLUDES)

FRC : 

allclean: clean
	rm -f make.depend make.depend.bak

clean:
#	rm -fr $(OBJSNLAMS) $(OBJS) $(LIB_OBJS) $(MAIN_OBJS) $(TARGET) $(JUNK) $(LOCI_INTERMEDIATE_FILES)
	rm -fr $(OBJSNLAMS) $(OBJSBEAM1D) $(OBJS) $(LIB_OBJS) $(MAIN_OBJS) $(TARGET) $(JUNK) $(LOCI_INTERMEDIATE_FILES)

distclean: clean
	rm -fr $(DEPEND_FILES) $(JUNK) $(LOCI_INTERMEDIATE_FILES)

ALL_OFILES = $(OBJS) 

DEPEND_FILES=$(subst .o,.d,$(ALL_OFILES))


#include automatically generated dependencies
include $(DEPEND_FILES)
