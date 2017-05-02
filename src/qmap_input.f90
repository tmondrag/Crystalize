! Tomas Mondragon
! Computer Scientist
! USACE ERDC Information Technology Laboratory
! Computational Analysis Branch
! Tomas.A.Mondragon@erdc.dren.mil
!   delivered 06 Feb 2016
MODULE qmap_input
  IMPLICIT NONE

CONTAINS

  ! Read input from Q-map data file into a lattice data structure
  ! USAGE: latticeData = read_input_file(filename)
  FUNCTION read_input_file(filename,outfileroot) RESULT(latticeData)
    USE mFileHandling, only: FILE, stderr
    USE basictypes, only: Lattice
    USE triangle_c_wrap, only: triangulateio,f_triangulateio,ctriangulate,copytriangles_c_to_f,init_outputs,deallocate_triangulateio
    USE triangle_input, only: read_shapes
    USE,INTRINSIC:: ISO_C_BINDING, only:C_NULL_CHAR,C_CHAR
    IMPLICIT NONE

    CHARACTER(LEN=256),INTENT(IN)                     :: filename
    CHARACTER(LEN=256),INTENT(IN)                     :: outfileroot
    TYPE(Lattice)                                     :: latticeData
    TYPE(triangulateio)                               :: c_shape,c_shape_out,c_shape_vorout
    TYPE(f_triangulateio)                             :: f_shape,f_shape_out
    CHARACTER(len=256)                                :: polyfilename
    INTEGER,DIMENSION(:),ALLOCATABLE                  :: Q
    TYPE(FILE)                                        :: inputFile
    INTEGER                                           :: ios
    CHARACTER(LEN=8)                                  :: mapType
    CHARACTER(kind=C_CHAR,len=256)                    :: flags
    INTEGER                                           :: num_cols,num_rows,col_i,row_j
    LOGICAL                                           :: qstate_present

    CALL inputFile%openReadOnly(filename)
    CALL skip_comments_input(inputFile)
    ! constant used in grain boundary energy calculations (electronVolts * ?m) Length scale could be nanometers, micrometers, millimeters, as long as it matches the scale used for point coordinates
    READ(inputFile%getFunit(),*,iostat=ios,err=700) latticeData%enrgPerLength
    CALL skip_comments_input(inputFile)
    ! temperature of crystal during evolution (useful during Potts generation)
    READ(inputFile%getFunit(),*,iostat=ios,err=600) latticeData%temperature
    CALL skip_comments_input(inputFile)
    ! number of possible discrete states
    READ(inputFile%getFunit(),*,iostat=ios,err=500) latticeData%numStates
    CALL skip_comments_input(inputFile)
    ! number of Monte Carlo sweeps
    READ(inputFile%getFunit(),*,iostat=ios,err=400) latticeData%mcSweeps
    CALL skip_comments_input(inputFile)
    ! x and y bounds
    READ(inputFile%getFunit(),*,iostat=ios,err=200) latticeData%xBounds,latticeData%yBounds
    CALL skip_comments_input(inputFile)
    ! grid spacing (for rigid rectililinear or hex meshes) first is x spacing and second is y spacing. a ratio of 1:1 will result in a square mesh, 1:0.8660254 is a regular hexagon mesh
    READ(inputFile%getFunit(),*,iostat=ios,err=300) latticeData%gridSpacing
    CALL skip_comments_input(inputFile)
    ! 0 =: aperiodic, 1 =: periodic in y, 2 =: periodic in x, 3 =: periodic everywhere
    READ(inputFile%getFunit(),*,iostat=ios,err=900) latticeData%periodicity
    CALL skip_comments_input(inputFile)
    ! Where should the point data come from?
    READ(inputFile%getFunit(),*,iostat=ios,err=800) mapType
    IF (mapType == "datafile") THEN
      CALL skip_comments_input(inputFile)
      READ(inputFile%getFunit(),*,iostat=ios,err=100) polyfilename
      ! Read in point and shape data
      CALL read_shapes(polyfilename,c_shape,f_shape)
      CALL output_triangle_to_vtk(f_shape_out,TRIM(outfileroot)//"init.vtk")
      flags = TRIM("prqcO")//C_NULL_CHAR
      CALL init_outputs(c_shape_out)
      CALL init_outputs(c_shape_vorout)
      CALL ctriangulate(flags,c_shape,c_shape_out,c_shape_vorout)
      CALL copytriangles_c_to_f(c_shape_out,f_shape_out)
      CALL output_triangle_to_vtk(f_shape_out,TRIM(outfileroot)//"first.vtk")
      CALL deallocate_triangulateio(f_shape, input=.TRUE., neigh=.FALSE., voroni=.FALSE.)
      CALL deallocate_triangulateio(f_shape_out, input=.FALSE., neigh=.FALSE., voroni=.FALSE.)
    ELSE IF (mapType == "random") THEN
    ELSE IF (mapType == "rect") THEN
    ELSE IF (mapType == "hex") THEN
    ELSE
      WRITE(stderr,'(A,A,A)') "ERROR: Misread map type from file ",TRIM(filename),"."
      STOP 1
    END IF
    RETURN
100 WRITE(stderr,'(A,A,A)') "ERROR: Misread point data file path name  in ",TRIM(filename),"."
    STOP 1
150 WRITE(stderr,'(A,A,A)') "ERROR: Misread number of nodes from file ",TRIM(filename),"."
    STOP 1
200 WRITE(stderr,'(A,A,A)') "ERROR: Misread x,y-bounds from file ",TRIM(filename),"."
    STOP 1
300 WRITE(stderr,'(A,A,A)') "ERROR: Misread grid spacing information from file ",TRIM(filename),"."
    STOP 1
400 WRITE(stderr,'(A,A,A)') "ERROR: Misread number of Evolver Monte Carlo Steps from file ",TRIM(filename),"."
    STOP 1
500 WRITE(stderr,'(A,A,A)') "ERROR: Misread number of states from file ",TRIM(filename),"."
    STOP 1
600 WRITE(stderr,'(A,A,A)') "ERROR: Misread temperature from file ",TRIM(filename),"."
    STOP 1
700 WRITE(stderr,'(A,A,A)') "ERROR: Misread energy scale factor from file ",TRIM(filename),"."
    STOP 1
800 WRITE(stderr,'(A,A,A)') "ERROR: Misread map type from file ",TRIM(filename),"."
    STOP 1
900 WRITE(stderr,'(A,A,A)') "ERROR: Misread periodicity from file ",TRIM(filename),"."
    STOP 1
  END FUNCTION read_input_file


  SUBROUTINE skip_comments_input(inputFile)
    USE mFileHandling, only: FILE
    IMPLICIT NONE
    TYPE(FILE),INTENT(INOUT)                          :: inputFile
    CHARACTER(LEN=256)                                :: buffer

    ReadComments: DO
      READ(inputFile%getFunit(), '(A)', END=100) buffer
      IF((buffer(1:1) /= "#") .AND. (TRIM(buffer) /= "")) exit ReadComments
    END DO ReadComments
    BACKSPACE(inputFile%getFunit())
100 CONTINUE
  END SUBROUTINE skip_comments_input
END MODULE qmap_input
