# Architecture-related flags
# Intel MIC arch, default is undefined
SUARCH_MIC  :=
SUARCH_GPU  :=

# comm mode: async or sync
COMM_MODE   := async
ifeq "$(COMM_MODE)" "sync"
CXXFLGS     += -DPGAMODE_SYNC
endif

SHELL       += -x
#CXX          = gcc
CXX          = mpicc
#CXXFLGS     += -g -Wall -DGSL_SPRNG
#CXXFLGS     += -g -Wall -DGSL_SPRNG -DPGAMODE -DPGA_NONBLOCK_MODE
#CXXFLGS     += -g -Wall -DSPRNG -DPGAMODE -DPGA_NONBLOCK_MODE
CXXFLGS     += -DSPRNG -DPGAMODE -DPGA_NONBLOCK_MODE -DT_PROFILING
#CXXFLGS     += -g -Wall -DSPRNG -DPGAMODE -DPGA_NONBLOCK_MODE -DT_PROFILING -DDEBUG_COMM
# use MPI_Ibsend(), not robust 'cause buffer policy differs on diff MPIs
#CXXFLGS     += -g -Wall -DSPRNG -DPGAMODE -DPGA_NONBLOCK_MODE -DDEBUG_COMM -DPGA_USE_IBSEND
SRCC         = c
SRCH         = h
# use static libraries
#STATICLINK   = -static
STATICLINK   =
# standalone obj items
#TARGETS      = log data addr pga gsl-sprng mysprng
TARGETS      = log data addr pga mysprng randseq
#TARGETS      = log data addr
OBJS         = $(TARGETS:=.o)
DEFTARGETS = myrng
DEFS       = $(DEFTARGETS:=.h)
# executables
MAINS        = ga

# external directories, set by env vars
#EXTDIRS      = $(MPICH_GM_HOME) $(GSL_HOME) $(SPRNG_HOME)
EXTDIRS      = $(MPICH_GM_HOME) $(SPRNG_HOME)
# stampede
#EXTDIRS      = $(MPICH_GM_HOME)/lib64 $(SPRNG_HOME)

# include files
INCPATH  = $(EXTDIRS:%=-I%/include)

# library paths
LIBPATH      = $(EXTDIRS:%=-L%/lib)
#LIBS         = -lgsl -lgslcblas -lsprng
LIBS         = -lsprng
#LIBS_DEFAULT = -lm -lpthread
LIBS_DEFAULT = -lm

all: $(MAINS) Makefile

$(MAINS): % : %.$(SRCC) $(DEFS) $(OBJS)
ifdef SUARCH_MIC
	@$(CXX) -xhost $(CXXFLGS) $(STATICLINK) -I. $(INCPATH) -o $@-$(COMM_MODE)-host $< $(OBJS) $(LIBPATH) $(LIBS) $(LIBS_DEFAULT) 
	@$(CXX) -mmic $(CXXFLGS) $(STATICLINK) -I. $(INCPATH) -o $@-$(COMM_MODE)-host $< $(OBJS) $(LIBPATH) $(LIBS) $(LIBS_DEFAULT) 
else
	@$(CXX) $(CXXFLGS) $(STATICLINK) -I. $(INCPATH) -o $@-$(COMM_MODE) $< $(OBJS) $(LIBPATH) $(LIBS) $(LIBS_DEFAULT) 
endif

# compile
$(DEFS):
$(OBJS): %.o: %.$(SRCC) %.$(SRCH)
	@$(CXX) $(CXXFLGS) -I. $(INCPATH) -o $@ -c $<
# link

# clean
clean:
	@rm -f $(OBJS) 
EXECS_ASYNC      = $(MAINS:=-async)
EXECS_SYNC      = $(MAINS:=-sync)
cleanall:
	@rm -f $(EXECS_ASYNC) $(EXECS_SYNC) $(OBJS) 
