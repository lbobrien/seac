program = seac 

DEBUG    = no
OPTIMIZE = no

source = $(wildcard *.f)
objects = lsode.o seac1.1.o 

F90 = gfortran
F90FLAGS := -cpp -fbacktrace
LDFLAGS =

# Debugging
ifeq ($(DEBUG),yes)
  F90FLAGS += -g -Wall -pedantic -std=f2008 -fbounds-check \
              -ffpe-trap=invalid,overflow,underflow
  LDFLAGS  += -g
endif

# Optimization
ifeq ($(OPTIMIZE),yes)
  F90FLAGS += -O3
endif

#===============================================================================
# Targets
#===============================================================================

all: $(program)
$(program): $(objects)
	$(F90) $(objects) $(LDFLAGS) -o $@
distclean: clean
	cd xml; make clean
clean:
	@rm -f *.o $(program)
neat:
	@rm -f *.o *.mod

#===============================================================================
# Rules
#===============================================================================

.SUFFIXES: .F90 .o
.PHONY: all clean neat

%.o: %.f
	$(F90) $(F90FLAGS) -c $<
