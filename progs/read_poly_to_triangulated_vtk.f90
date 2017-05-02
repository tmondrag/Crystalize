! Tomas Mondragon
! Computer Scientist
! USACE ERDC Information Technology Laboratory
! Computational Analysis Branch
! Tomas.A.Mondragon@erdc.dren.mil
!   delivered 24 Apr 2017

! Test of integration of triangle library methods
! USAGE: read_poly_to_vtk flags infile outfile
! This program reads in the planar straight line graph in the file specified by the second argument
! and triangulates the shapes described within, according to the flags in the first argument,
! and outputs the result to a vtk file (name specified in the third argument)

! The planar straight line graph need to be in the same PSLG format as those readable by Triangle
! https://www.cs.cmu.edu/~quake/triangle.poly.html
! and the flags are the same as Triangle's flags:
! -p Triangulates a Planar Straight Line Graph (.poly file).
! -r Refines a previously generated mesh.
! -q Quality mesh generation with no angles smaller than 20 degrees. An alternate minimum angle may be specified after the `q'.
! -a Imposes a maximum triangle area constraint. A fixed area constraint (that applies to every triangle) may be specified after the `a', or varying area constraints may be read from a .poly file or .area file.
! -u Imposes a user-defined constraint on triangle size.
! -A Assigns a regional attribute to each triangle that identifies what segment-bounded region it belongs to.
! -c Encloses the convex hull with segments.
! -D Conforming Delaunay: use this switch if you want all triangles in the mesh to be Delaunay, and not just constrained Delaunay; or if you want to ensure that all Voronoi vertices lie within the triangulation.
! -j Jettisons vertices that are not part of the final triangulation from the output .node file (including duplicate input vertices and vertices ``eaten'' by holes).
! -I Suppresses mesh iteration numbers.
! -O Suppresses holes: ignores the holes in the .poly file.
! -X Suppresses exact arithmetic.
! -z Numbers all items starting from zero (rather than one). Note that this switch is normally overrided by the value used to number the first vertex of the input .node or .poly file. However, this switch is useful when calling Triangle from another program.
! -o2 Generates second-order subparametric elements with six nodes each.
! -Y Prohibits the insertion of Steiner points on the mesh boundary. If specified twice (-YY), it prohibits the insertion of Steiner points on any segment, including internal segments.
! -S Specifies the maximum number of added Steiner points.
! -i Uses the incremental algorithm for Delaunay triangulation, rather than the divide-and-conquer algorithm.
! -F Uses Steven Fortune's sweepline algorithm for Delaunay triangulation, rather than the divide-and-conquer algorithm.
! -l Uses only vertical cuts in the divide-and-conquer algorithm. By default, Triangle uses alternating vertical and horizontal cuts, which usually improve the speed except with vertex sets that are small or short and wide. This switch is primarily of theoretical interest.
! -s Specifies that segments should be forced into the triangulation by recursively splitting them at their midpoints, rather than by generating a constrained Delaunay triangulation. Segment splitting is true to Ruppert's original algorithm, but can create needlessly small triangles. This switch is primarily of theoretical interest.
! -C Check the consistency of the final mesh. Uses exact arithmetic for checking, even if the -X switch is used. Useful if you suspect Triangle is buggy.
! -Q Quiet: Suppresses all explanation of what Triangle is doing, unless an error occurs.
! -V Verbose: Gives detailed information about what Triangle is doing. Add more `V's for increasing amount of detail. `-V' gives information on algorithmic progress and detailed statistics.
! -h Help: Displays complete instructions.

program read_poly_to_traingulated_vtk
  USE,INTRINSIC::ISO_FORTRAN_ENV, only: stderr => ERROR_UNIT
  USE triangle_c_wrap, only: triangulateio,f_triangulateio,ctriangulate,copytriangles_c_to_f,init_outputs,deallocate_triangulateio
  USE triangle_input, only: read_shapes
  USE triangle_output, only: output_triangle_to_vtk
  USE,INTRINSIC:: ISO_C_BINDING, only:C_NULL_CHAR,C_CHAR
  IMPLICIT NONE
  TYPE(triangulateio)       :: c_shape,c_shape_out,c_shape_vorout
  TYPE(f_triangulateio)     :: f_shape,f_shape_out
  CHARACTER(kind=C_CHAR,len=256) :: infile, outfile, flags
  INTEGER                   :: numargs
  LOGICAL                   :: neighflag

  numargs = COMMAND_ARGUMENT_COUNT()
  IF (numargs < 3) THEN
    WRITE(stderr,'(A)') "USAGE: read_poly_to_vtk flags infile outfile"
    STOP 1
  ELSE IF (numargs > 3) THEN
    WRITE(stderr,'(A)') "Only first three arguments used, other arguments ignored."
    WRITE(stderr,'(A)') "USAGE: read_poly_to_vtk flags infile outfile"
  END IF
  CALL GET_COMMAND_ARGUMENT(1,flags)
  CALL GET_COMMAND_ARGUMENT(2,infile)
  CALL GET_COMMAND_ARGUMENT(3,outfile)

  neighflag = (0 .NE. SCAN(flags,"n"))
  flags = TRIM(flags)//C_NULL_CHAR
  CALL read_shapes(infile,c_shape,f_shape)
  CALL init_outputs(c_shape_out)
  CALL init_outputs(c_shape_vorout)
  CALL ctriangulate(flags,c_shape,c_shape_out,c_shape_vorout)
  CALL copytriangles_c_to_f(c_shape_out,f_shape_out)
  CALL output_triangle_to_vtk(f_shape_out,outfile)
  CALL deallocate_triangulateio(f_shape, input=.TRUE., neigh=neighflag, voroni=.FALSE.)
  CALL deallocate_triangulateio(f_shape_out, input=.FALSE., neigh=neighflag, voroni=.FALSE.)
end program read_poly_to_traingulated_vtk
