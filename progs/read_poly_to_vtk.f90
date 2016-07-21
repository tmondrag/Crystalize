PROGRAM read_poly_to_vtk
  USE filehandling, only: stderr
  USE triangle_c_wrap, only: triangulateio,f_triangulateio
  USE triangle_input, only: read_shapes
  USE triangle_output, only: output_triangle_to_vtk
  IMPLICIT NONE
  TYPE(triangulateio)       :: c_shape
  TYPE(f_triangulateio)     :: f_shape
  CHARACTER(len=256)        :: infile, outfile
  INTEGER                   :: numargs

  numargs = COMMAND_ARGUMENT_COUNT()
  IF (numargs < 2) THEN
    WRITE(stderr,'(A)') "USAGE: read_poly_to_vtk infile outfile"
    STOP 1
  ELSE IF (numargs > 2) THEN
    WRITE(stderr,'(A)') "Only first two arguments used, other arguments ignored."
    WRITE(stderr,'(A)') "USAGE: read_poly_to_vtk infile outfile"
  END IF
  CALL GET_COMMAND_ARGUMENT(1,infile)
  CALL GET_COMMAND_ARGUMENT(2,outfile)

  CALL read_shapes(infile,c_shape,f_shape)
  CALL output_triangle_to_vtk(f_shape,outfile)
END PROGRAM read_poly_to_vtk
