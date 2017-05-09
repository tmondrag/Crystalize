! crystalize.f90
PROGRAM crystalize
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: stderr => ERROR_UNIT
  USE basictypes, only: Lattice
  USE qmap_input, ONLY: read_input_file
  IMPLICIT NONE
  INTEGER                                           :: numargs
  CHARACTER(LEN=256)                                :: infilename, outfileroot
  TYPE(Lattice)                                     :: latticeData

  numargs = COMMAND_ARGUMENT_COUNT()
  IF (numargs < 2) THEN
    WRITE(stderr,'(A)') "USAGE: crystalize infilename outfileroot"
    STOP 1
  ELSE IF (numargs > 2) THEN
    WRITE(stderr,'(A)') "Only first two arguments used, other arguments ignored."
    WRITE(stderr,'(A)') "USAGE: crystalize infilename outfileroot"
  END IF
  CALL GET_COMMAND_ARGUMENT(1,infilename)
  CALL GET_COMMAND_ARGUMENT(2,outfileroot)
  ! Read in the preamble data file, which will lead to the polydata file, and
  ! write out vtk files for the original polydata and the triangulated polydata
  latticeData = read_input_file(infilename,outfileroot)
END PROGRAM crystalize
