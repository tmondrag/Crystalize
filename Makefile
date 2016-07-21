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
TESTS = $(addprefix $(BINDIR)/,copy_qmap_test quick_qmap_analysis_test read_qmap_to_vtk read_poly_to_vtk read_poly_to_triangulated_vtk qmap_to_lattice_vtk boundary_smoother_test)
# Programs - name your final executables here
PROGRAMS = $(addprefix $(BINDIR)/,generate_qmap generate_qmap_from_file qmap_to_rect_bounds qmap_to_oct_bounds qmap_to_triangulated_vtk)# process_qmap
# All directive. What make will make by default
all:$(TESTS) $(PROGRAMS)

# Linker dependencies. List the object files that need to be linked to 
# create the executables.

$(BINDIR)/boundary_smoother_test: $(addprefix $(OBJDIR)/,boundary_smoother_test.o boundarysmoother.o boundary_output.o boundaries.o basictypes.o filehandling.o math_tools.o weaver.o boundarycells.o)

$(BINDIR)/copy_qmap_test: $(addprefix $(OBJDIR)/,copy_qmap_test.o output.o input.o basictypes.o filehandling.o math_tools.o)

$(BINDIR)/generate_qmap_from_file: $(addprefix $(OBJDIR)/,generate_qmap_from_file.o input.o output.o generators.o basictypes.o math_tools.o filehandling.o )

$(BINDIR)/generate_qmap: $(addprefix $(OBJDIR)/,generate_qmap.o generators.o basictypes.o math_tools.o output.o filehandling.o)

$(BINDIR)/qmap_to_lattice_vtk: $(addprefix $(OBJDIR)/,qmap_to_lattice_vtk.o input.o basictypes.o output.o weaver.o filehandling.o math_tools.o)

$(BINDIR)/qmap_to_oct_bounds: $(addprefix $(OBJDIR)/,qmap_to_oct_bounds.o weaver.o basictypes.o filehandling.o boundarysmoother.o input.o boundaries.o boundary_output.o boundarysegments.o boundarycells.o math_tools.o)

$(BINDIR)/qmap_to_rect_bounds: $(addprefix $(OBJDIR)/,qmap_to_rect_bounds.o input.o filehandling.o basictypes.o weaver.o boundaries.o boundarysegments.o boundarysmoother.o boundary_output.o boundarycells.o math_tools.o)

$(BINDIR)/quick_qmap_analysis_test: $(addprefix $(OBJDIR)/, quick_qmap_analysis_test.o input.o qmap_utility.o basictypes.o filehandling.o math_tools.o)

$(BINDIR)/qmap_to_triangulated_vtk: $(addprefix $(OBJDIR)/,qmap_to_triangulated_vtk.o boundaries.o input.o basictypes.o filehandling.o weaver.o boundarysegments.o boundarysmoother.o polygons.o triangle_c_wrap.o triangle_output.o math_tools.o boundarycells.o qmap_utility.o) triangle_lib/triangle.o

$(BINDIR)/read_poly_to_triangulated_vtk: $(addprefix $(OBJDIR)/, read_poly_to_triangulated_vtk.o triangle_output.o triangle_c_wrap.o triangle_input.o filehandling.o math_tools.o) triangle_lib/triangle.o

$(BINDIR)/read_poly_to_vtk: $(addprefix $(OBJDIR)/,read_poly_to_vtk.o triangle_output.o triangle_c_wrap.o triangle_input.o filehandling.o math_tools.o) triangle_lib/triangle.o

$(BINDIR)/read_qmap_to_vtk: $(addprefix $(OBJDIR)/, read_qmap_to_vtk.o input.o generators.o qmap_utility.o basictypes.o output.o filehandling.o math_tools.o )

#Executable Object dependencies. List object files that depend on modules from other object files here

$(OBJDIR)/boundary_smoother_test.o: $(addprefix $(OBJDIR)/,boundary_output.o boundaries.o boundarysmoother.o)

$(OBJDIR)/copy_qmap_test.o: $(addprefix $(OBJDIR)/,output.o input.o basictypes.o filehandling.o)

$(OBJDIR)/generate_qmap_from_file.o: $(addprefix $(OBJDIR)/,input.o output.o generators.o basictypes.o filehandling.o)

$(OBJDIR)/generate_qmap.o: $(addprefix $(OBJDIR)/,output.o generators.o math_tools.o filehandling.o basictypes.o)

$(OBJDIR)/qmap_to_lattice_vtk.o: $(addprefix $(OBJDIR)/,input.o basictypes.o output.o weaver.o filehandling.o)

$(OBJDIR)/qmap_to_oct_bounds.o: $(addprefix $(OBJDIR)/,weaver.o basictypes.o filehandling.o boundarysmoother.o input.o boundaries.o boundary_output.o boundarysegments.o)

$(OBJDIR)/qmap_to_rect_bounds.o: $(addprefix $(OBJDIR)/,input.o filehandling.o basictypes.o weaver.o boundaries.o boundarysegments.o boundarysmoother.o boundary_output.o)

$(OBJDIR)/qmap_to_triangulated_vtk.o: $(addprefix $(OBJDIR)/, boundaries.o input.o basictypes.o filehandling.o weaver.o boundarysegments.o boundarysmoother.o polygons.o triangle_output.o triangle_c_wrap.o qmap_utility.o)

$(OBJDIR)/quick_qmap_analysis_test.o: $(addprefix $(OBJDIR)/,qmap_utility.o input.o basictypes.o filehandling.o)

$(OBJDIR)/read_poly_to_triangulated_vtk.o: $(addprefix $(OBJDIR)/,triangle_output.o triangle_c_wrap.o triangle_input.o filehandling.o)

$(OBJDIR)/read_poly_to_vtk.o: $(addprefix $(OBJDIR)/,triangle_output.o triangle_c_wrap.o triangle_input.o filehandling.o)

$(OBJDIR)/read_qmap_to_vtk.o: $(addprefix $(OBJDIR)/,input.o generators.o qmap_utility.o basictypes.o output.o filehandling.o)

#Module Object dependencies. List module object files that depend on modules from other object files here

$(OBJDIR)/basictypes.o: $(OBJDIR)/math_tools.o

$(OBJDIR)/boundaries.o: $(OBJDIR)/basictypes.o

$(OBJDIR)/boundary_output.o: $(addprefix $(OBJDIR)/,boundarysegments.o filehandling.o boundaries.o)

$(OBJDIR)/boundarycells.o: $(addprefix $(OBJDIR)/,basictypes.o boundaries.o filehandling.o)

$(OBJDIR)/boundarypoints.o: $(addprefix $(OBJDIR)/,basictypes.o boundaries.o)

$(OBJDIR)/boundarysegments.o: $(addprefix $(OBJDIR)/,basictypes.o boundaries.o weaver.o boundarycells.o filehandling.o)

$(OBJDIR)/boundarysmoother.o: $(addprefix $(OBJDIR)/,basictypes.o boundaries.o filehandling.o)

$(OBJDIR)/generators.o: $(addprefix $(OBJDIR)/,basictypes.o math_tools.o filehandling.o)

$(OBJDIR)/input.o: $(addprefix $(OBJDIR)/,math_tools.o filehandling.o basictypes.o)

$(OBJDIR)/output.o: $(addprefix $(OBJDIR)/,math_tools.o filehandling.o basictypes.o)

$(OBJDIR)/polygons.o: $(addprefix $(OBJDIR)/, boundaries.o triangle_c_wrap.o)

$(OBJDIR)/qmap_utility.o: $(addprefix $(OBJDIR)/,math_tools.o basictypes.o filehandling.o)

$(OBJDIR)/triangle_c_wrap.o: $(OBJDIR)/filehandling.o triangle_lib/triangle.o

$(OBJDIR)/triangle_input.o: $(addprefix $(OBJDIR)/,triangle_c_wrap.o filehandling.o math_tools.o)

$(OBJDIR)/triangle_output.o: $(addprefix $(OBJDIR)/,triangle_c_wrap.o filehandling.o)

$(OBJDIR)/weaver.o: $(addprefix $(OBJDIR)/,basictypes.o math_tools.o filehandling.o)

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

