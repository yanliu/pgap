SHELL       += -x 
#CXX          = gcc
CXX          = mpicc
#CXXFLGS     += -g -Wall -DGSL_SPRNG
#CXXFLGS     += -g -Wall -DGSL_SPRNG -DPGAMODE -DPGA_NONBLOCK_MODE
#CXXFLGS     += -g -Wall -DSPRNG -DPGAMODE -DPGA_NONBLOCK_MODE
# sync mode
#CXXFLGS     += -DSPRNG -DPGAMODE -DPGA_NONBLOCK_MODE -DPGAMODE_SYNC
# async mode
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
DEFTARGETS = ga myrng
DEFS       = $(DEFTARGETS:=.h)
# executables
MAINS        = ga-async
MAINOBJS         = $(MAINS:=.o)

# external directories, set by env vars
#EXTDIRS      = $(MPICH_GM_HOME) $(GSL_HOME) $(SPRNG_HOME)
EXTDIRS      = $(MPICH_GM_HOME) $(SPRNG_HOME)

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
	@$(CXX) $(CXXFLGS) $(STATICLINK) -I. $(INCPATH) -o $@ $< $(OBJS) $(LIBPATH) $(LIBS) $(LIBS_DEFAULT)

# compile
$(DEFS):
$(OBJS): %.o: %.$(SRCC) %.$(SRCH)
	@$(CXX) $(CXXFLGS) -I. $(INCPATH) -o $@ -c $<
# link

# clean
clean:
#	@rm -f $(MAINS) $(OBJS) $(MAINOBJS)
	@rm -f $(OBJS) $(MAINOBJS)
