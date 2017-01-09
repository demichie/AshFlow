AUTOMAKE_OPTIONS = foreign

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build (the names of the final binaries)
bin_PROGRAMS = AshFlow
bindir=@bindir@

# list of sources for the binary
AshFlow_SOURCES = \
		tools.f90 \
        current.f90 \
        environment.f90 \
        particles.f90 \
        mixture.f90 \
		solver_flow.f90 \        
		inpout.f90 \
        flow.f90 \
		flow_model.f90

FCFLAGS         = -g -O3 -ffast-math
FCDEBUGFLAGS    = -Wall -g

LIBS            = -lz -lm

FC              = gfortran
MPFC            = gfortran
LINKER          = gfortran

