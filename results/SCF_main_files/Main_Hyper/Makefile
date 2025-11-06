#INCLUDE=.:/usr/local/include 
INCLUDE=/usr/include 
DIR_GEOMETRY = ../A_Modules_Geometry/Bulk

DIR_GENERAL = ../A_Routines_General
DIR_POLYMERS = ../A_Modules_Polymer
DIR_INTERACTIONS = ../A_Modules_Interaction
DIR_ANALYSIS = ../A_Analysis

OBJS =	\
  mod_global.o \
  mod_mt.o \
  mod_fourier.o \
  mod_propagate.o \
  mod_iterate.o \
  mod_hyper.o \
  sub_wrap.o \
  mod_hyper_extras.o \
  mod_FloryHuggins_Solution.o 

OBJS_ANL =	\
  mod_global.o \
  mod_mt.o \
  mod_fourier.o \
  mod_propagate.o \
  mod_iterate.o \
  mod_hyper_analysis.o \
  sub_wrap.o \
  mod_hyper_extras_analysis.o \
  mod_FloryHuggins_Solution.o 

LIBS =-lfftw3_omp -lfftw3 -lm

CC = g++
CFLAGS = -O2 -g
FC = gfortran 
FFLAGS = -O
F90 = gfortran
F90FLAGS = -O3 -g -pipe -fno-omit-frame-pointer -fopenmp 
#F90FLAGS = -O3 -pipe -fomit-frame-pointer  
LDFLAGS = -fopenmp 

all: 

main_hyper: $(OBJS) main_hyper.f90
	$(F90) -I$(INCLUDE) $(F90FLAGS) -c $@.f90 
	$(F90) -I$(INCLUDE) $(LDFLAGS)   $(OBJS) $@.o -L/usr/lib:/usr/local/lib:/usr/lib64 $(LIBS) -o $@.exe 

analysis: $(OBJS_ANL) analysis.f90
	$(F90) -I$(INCLUDE) $(F90FLAGS) -c $@.f90 
	$(F90) -I$(INCLUDE) $(LDFLAGS)   $(OBJS_ANL) $@.o -L/usr/lib:/usr/local/lib:/usr/lib64 $(LIBS) -o $@.exe 

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90 
	$(F90)  -I$(INCLUDE) $(F90FLAGS) -c $<

%.o: $(DIR_GENERAL)/%.f90
	$(F90)  -I$(INCLUDE) $(F90FLAGS) -c $<

%.o: $(DIR_GEOMETRY)/%.f90
	$(F90)  -I$(INCLUDE) $(F90FLAGS) -c $<

%.o: $(DIR_POLYMERS)/%.f90
	$(F90)  -I$(INCLUDE) $(F90FLAGS) -c $<

%.o: $(DIR_INTERACTIONS)/%.f90
	$(F90)  -I$(INCLUDE) $(F90FLAGS) -c $<

%.o: $(DIR_ANALYSIS)/%.f90
	$(F90)  -I$(INCLUDE) $(F90FLAGS) -c $<

clean:
	rm -f *.o *.mod fort.*
