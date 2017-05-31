#$preamble
# compilers
FC = gfortran # fortran

# Flags
# Fortran Debugging
FCLFLAGS = -g -fbounds-check -fbacktrace -Wall
FCFLAGS =  -g -fbounds-check -fbacktrace -Wall -J $(OBJDIR)

# Fortran optimization
# FCLFLAGS = -O2
# FCFLAGS = -O2 -J $(OBJDIR)

# Fortran flags forall (e.g. look for system .mod files, *sometimes?* required in gfortran)
# FCLFLAGS += -I/usr/include

# Directories
BINDIR = bin
OBJDIR = obj
SRCDIR = src
ANIDIR = anim
PROGDIR = progs

# Tests - name your test executables here
TESTS = $(addprefix $(BINDIR)/,read_poly_to_vtk read_poly_to_triangulated_vtk crystalize)
# Programs - name your final executables here
PROGRAMS = $(addprefix $(BINDIR)/,read_poly_to_vtk crystalize)
# All directive. What make will make by default
all:$(TESTS) $(PROGRAMS)

# Linker dependencies. List the object files that need to be linked to
# create the executables.

$(BINDIR)/crystalize: $(addprefix $(OBJDIR)/,crystalize.o qmap_input.o trans_lattice_triangle.o triangle_output.o basictypes.o mFileHandling.o triangle_c_wrap.o triangle_input.o evolvers.o energyCalc.o lattice_output.o) triangle_lib/triangle.o

$(BINDIR)/read_poly_to_vtk: $(addprefix $(OBJDIR)/,read_poly_to_vtk.o triangle_output.o triangle_c_wrap.o triangle_input.o mFileHandling.o math_tools.o) triangle_lib/triangle.o

$(BINDIR)/read_poly_to_triangulated_vtk: $(addprefix $(OBJDIR)/, read_poly_to_triangulated_vtk.o triangle_output.o triangle_c_wrap.o triangle_input.o mFileHandling.o math_tools.o) triangle_lib/triangle.o

#Executable Object dependencies. List object files that depend on modules from other object files here

$(OBJDIR)/crystalize.o: $(addprefix $(OBJDIR)/,basictypes.o qmap_input.o trans_lattice_triangle.o)

$(OBJDIR)/read_poly_to_vtk.o: $(addprefix $(OBJDIR)/,triangle_output.o triangle_c_wrap.o triangle_input.o mFileHandling.o)

#Module Object dependencies. List module object files that depend on modules from other object files here

$(OBJDIR)/basictypes.o: $(addprefix $(OBJDIR)/,math_tools.o triangle_c_wrap.o)

$(OBJDIR)/triangle_c_wrap.o: triangle_lib/triangle.o

$(OBJDIR)/triangle_input.o: $(addprefix $(OBJDIR)/,triangle_c_wrap.o mFileHandling.o math_tools.o)

$(OBJDIR)/triangle_output.o: $(addprefix $(OBJDIR)/,triangle_c_wrap.o mFileHandling.o)

$(OBJDIR)/qmap_input.o: $(addprefix $(OBJDIR)/,basictypes.o mFileHandling.o trans_lattice_triangle.o triangle_c_wrap.o triangle_input.o triangle_output.o evolvers.o lattice_output.o)

$(OBJDIR)/trans_lattice_triangle.o: $(addprefix $(OBJDIR)/,basictypes.o triangle_c_wrap.o energyCalc.o)

$(OBJDIR)/evolvers.o: $(addprefix $(OBJDIR)/,triangle_c_wrap.o basictypes.o energyCalc.o)

$(OBJDIR)/enrgyCalc.o: $(OBJDIR)/triangle_c_wrap.o

$(OBJDIR)/lattice_output.o: $(addprefix $(OBJDIR)/,basictypes.o mFileHandling.o)

# target for Triangle lib
triangle_lib/triangle.o : triangle_lib/triangle.c triangle_lib/triangle.h
	cd triangle_lib && make triangle.o

#Default object file recipes

$(BINDIR)/%: $(OBJDIR)/%.o | $(BINDIR)
	$(FC) $(FCLFLAGS) -o $@ $^

$(OBJDIR)/%.o: $(SRCDIR)/%.f90 | $(OBJDIR)
	$(FC) $(FCFLAGS) -o $@ -c $<

$(OBJDIR)/%.o: $(PROGDIR)/%.f90 | $(OBJDIR)
	$(FC) $(FCFLAGS) -o $@ -c $<

# directory setup

$(BINDIR):
	mkdir $(BINDIR)

$(OBJDIR):
	mkdir $(OBJDIR)

$(SRCDIR):
	mkdir $(SRCDIR)
	mv *.f90 $(SRCDIR)

$(ANIDIR):
	mkdir $(ANIDIR)

$(PROGDIR):
	mkdir $(PROGDIR)

# Utility targets
.PHONY: all cleanobj cleanbin clean cleantriobj

cleanobj:
	rm -f *.mod $(OBJDIR)/*.o $(OBJDIR)/*.mod $(OBJDIR)/*.MOD $(OBJDIR)/*~

cleantriobj:
	rm -f triangle_lib/*.mod triangle_lib/*.o

cleanbin:
	rm -rf $(BINDIR)/*

clean: cleanbin cleanobj cleantriobj
	rm -f *~ $(SRCDIR)/*~ *.mod

tests: $(TESTS)

programs: $(PROGRAMS)
