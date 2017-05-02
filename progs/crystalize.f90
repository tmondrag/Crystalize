! crystalize.f90
program crystalize

  implicit none

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
end program crystalize
