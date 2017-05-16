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
    USE triangle_c_wrap, only: C_REAL,allocate_points_ftoc,allocate_segments,allocate_holes,allocate_regions
    USE triangle_input, only: read_shapes
    USE triangle_output, only: output_triangle_to_vtk, output_triangle_to_poly
    USE,INTRINSIC:: ISO_C_BINDING, only:C_NULL_CHAR,C_CHAR,C_INT
    IMPLICIT NONE

    CHARACTER(LEN=256),INTENT(IN)                     :: filename
    CHARACTER(LEN=256),INTENT(IN)                     :: outfileroot
    CHARACTER(LEN=256)                                :: infileroot
    TYPE(Lattice)                                     :: latticeData
    TYPE(triangulateio)                               :: c_shape,c_shape_out,c_shape_vorout
    TYPE(f_triangulateio)                             :: f_shape,f_shape_out
    CHARACTER(len=256)                                :: polyfilename
    TYPE(FILE)                                        :: inputFile
    INTEGER                                           :: ios
    CHARACTER(LEN=8)                                  :: mapType
    CHARACTER(kind=C_CHAR,len=256)                    :: flags
    INTEGER                                           :: num_x,num_y
    INTEGER                                           :: i,j,k,l
    REAL(C_REAL)                                      :: remainder
    REAL(C_REAL),DIMENSION(:),ALLOCATABLE             :: xPosArray, yPosArray

    CHARACTER(LEN=1),PARAMETER                        :: dirSep = '/'
    INTEGER                                           :: ls

    ls = INDEX(filename,dirSep,.TRUE.)
    infileroot = filename(1:ls)

    f_shape%numberofpoints             = 0
    f_shape%numberofpointattributes    = 0
    f_shape%numberofholes              = 0
    f_shape%numberofregions            = 0
    f_shape%numberofsegments           = 0
    f_shape%numberoftriangles          = 0
    f_shape%numberofcorners            = 0
    f_shape%numberoftriangleattributes = 0
    f_shape%numberofedges              = 0
    CALL allocate_holes(f_shape,c_shape)
    CALL allocate_regions(f_shape,c_shape)

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
    ! Point data is in a seperate pslg format file
    IF (mapType == "datafile") THEN
      CALL skip_comments_input(inputFile)
      READ(inputFile%getFunit(),'(A)',iostat=ios,err=100) polyfilename
      polyfilename = TRIM(infileroot)//TRIM(polyfilename)
      ! Read in point and shape data
      CALL read_shapes(polyfilename,c_shape,f_shape)

    ! Generate point data, no Q-state, points at random positions
    ! for use with the positional evolver and the orientation based boundary finder
    ELSE IF (mapType == "random") THEN
      CALL init_random_seed()
      CALL skip_comments_input(inputFile)
      READ(inputFile%getFunit(),*,iostat=ios,err=150) f_shape%numberofpoints
      f_shape%numberofpoints = f_shape%numberofpoints + 4
      ! attributes: Q-state = 0, grain_id
      f_shape%numberofpointattributes    = 2
      CALL allocate_points_ftoc(f_shape,c_shape)
      f_shape%pointlist(SIZE(f_shape%pointlist)-7) = latticeData%xBounds(1)
      f_shape%pointlist(SIZE(f_shape%pointlist)-6) = latticeData%yBounds(1)
      f_shape%pointlist(SIZE(f_shape%pointlist)-5) = latticeData%xBounds(2)
      f_shape%pointlist(SIZE(f_shape%pointlist)-4) = latticeData%yBounds(1)
      f_shape%pointlist(SIZE(f_shape%pointlist)-3) = latticeData%xBounds(2)
      f_shape%pointlist(SIZE(f_shape%pointlist)-2) = latticeData%yBounds(2)
      f_shape%pointlist(SIZE(f_shape%pointlist)-1) = latticeData%xBounds(1)
      f_shape%pointlist(SIZE(f_shape%pointlist))   = latticeData%yBounds(2)
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-7) = 0
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-6) = SIZE(f_shape%pointattributelist)-3
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-5) = 0
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-4) = SIZE(f_shape%pointattributelist)-2
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-3) = 0
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-2) = SIZE(f_shape%pointattributelist)-1
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-1) = 0
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist))   = SIZE(f_shape%pointattributelist)
      f_shape%pointmarkerlist(SIZE(f_shape%pointmarkerlist)-3) = 1
      f_shape%pointmarkerlist(SIZE(f_shape%pointmarkerlist)-2) = 1
      f_shape%pointmarkerlist(SIZE(f_shape%pointmarkerlist)-1) = 1
      f_shape%pointmarkerlist(SIZE(f_shape%pointmarkerlist))   = 1

      f_shape%numberofsegments = 4
      CALL allocate_segments(f_shape,c_shape)
      f_shape%segmentlist(1) = SIZE(f_shape%pointmarkerlist)-3
      f_shape%segmentlist(2) = SIZE(f_shape%pointmarkerlist)-2
      f_shape%segmentlist(3) = SIZE(f_shape%pointmarkerlist)-2
      f_shape%segmentlist(4) = SIZE(f_shape%pointmarkerlist)-1
      f_shape%segmentlist(5) = SIZE(f_shape%pointmarkerlist)-1
      f_shape%segmentlist(6) = SIZE(f_shape%pointmarkerlist)
      f_shape%segmentlist(7) = SIZE(f_shape%pointmarkerlist)
      f_shape%segmentlist(8) = SIZE(f_shape%pointmarkerlist)-3
      f_shape%segmentmarkerlist(1) = 1
      f_shape%segmentmarkerlist(2) = 1
      f_shape%segmentmarkerlist(3) = 1
      f_shape%segmentmarkerlist(4) = 1

      DO i=1,f_shape%numberofpoints - 4
        f_shape%pointattributelist(2*i-1) = 0                   ! q-state = 0
        f_shape%pointattributelist(2*i)   = i                   ! grain_id = vertex index initially
        CALL RANDOM_NUMBER(remainder)
        f_shape%pointlist(2*i-1) = remainder*(latticeData%xBounds(2)-latticeData%xBounds(1))
        f_shape%pointlist(2*i-1) = f_shape%pointlist(2*i-1) + latticeData%xBounds(1)
        CALL RANDOM_NUMBER(remainder)
        f_shape%pointlist(2*i)   = remainder*(latticeData%yBounds(2)-latticeData%yBounds(1))
        f_shape%pointlist(2*i)   = f_shape%pointlist(2*i) + latticeData%yBounds(1)
        f_shape%pointmarkerlist(i) = 0
      END DO

    ! Generate point data, random Q-state, points arranged in a rectilinear grid
    ! for use with the Q-state evolver and the Q-state based boundary finder
    ELSE IF (mapType == "rect") THEN
      CALL init_random_seed()
      num_x = INT((latticeData%xBounds(2) - latticeData%xBounds(1))/latticeData%gridSpacing(1))
      remainder = (latticeData%xBounds(2) - latticeData%xBounds(1)) - REAL(num_x)*latticeData%gridSpacing(1)
      IF (remainder < latticeData%gridSpacing(1)/2.0) THEN
        ALLOCATE(xPosArray(1:num_x+1),STAT=ios)
        IF ( ios /= 0) THEN
          WRITE(stderr,'(A,i16, A)') "ERROR : Allocation request denied. Array needed ",num_x+1," elements."
          STOP 1
        END IF
        xPosArray(1:num_x) = (/(latticeData%xBounds(1)+latticeData%gridSpacing(1)*REAL(i),i = 0,(num_x-1),1)/)
        xPosArray(num_x+1) = latticeData%xBounds(2)
      ELSE
        ALLOCATE(xPosArray(1:num_x+2),STAT=ios)
        IF ( ios /= 0) THEN
          WRITE(stderr,'(A,i16, A)') "ERROR : Allocation request denied. Array needed ",num_x+2," elements."
          STOP 1
        END IF
        xPosArray(1:num_x+1) = (/(latticeData%xBounds(1)+latticeData%gridSpacing(1)*REAL(i),i = 0,num_x,1)/)
        xPosArray(num_x+2) = latticeData%xBounds(2)
      END IF
      num_y = INT((latticeData%yBounds(2) - latticeData%yBounds(1))/latticeData%gridSpacing(2))
      remainder = (latticeData%yBounds(2) - latticeData%yBounds(1)) - REAL(num_y)*latticeData%gridSpacing(2)
      IF (remainder < latticeData%gridSpacing(2)/2.0) THEN
        ALLOCATE(yPosArray(1:num_y+1),STAT=ios)
        IF ( ios /= 0) THEN
          WRITE(stderr,'(A,i16, A)') "ERROR : Allocation request denied. Array needed ",num_y+1," elements."
          STOP 1
        END IF
        yPosArray(1:num_y) = (/(latticeData%yBounds(1)+latticeData%gridSpacing(2)*REAL(i),i = 0,(num_y-1),1)/)
        yPosArray(num_y+1) = latticeData%yBounds(2)
      ELSE
        ALLOCATE(yPosArray(1:num_y+2),STAT=ios)
        IF ( ios /= 0) THEN
          WRITE(stderr,'(A,i16, A)') "ERROR : Allocation request denied. Array needed ",num_y+2," elements."
          STOP 1
        END IF
        yPosArray(1:num_y+1) = (/(latticeData%yBounds(1)+latticeData%gridSpacing(2)*REAL(i),i = 0,num_y,1)/)
        yPosArray(num_y+2) = latticeData%yBounds(2)
      END IF
      f_shape%numberofpoints             = SIZE(xPosArray)*SIZE(yPosArray) + 4
      f_shape%numberofpointattributes    = 2
      CALL allocate_points_ftoc(f_shape,c_shape)
      DO j = 1, SIZE(yPosArray)
        DO i = 1,SIZE(xPosArray)
          f_shape%pointlist(2*(j-1)*SIZE(xPosArray)+2*i-1) = xPosArray(i)
          f_shape%pointlist(2*(j-1)*SIZE(xPosArray)+2*i) = yPosArray(j)
          CALL RANDOM_NUMBER(remainder)
          remainder = REAL(latticeData%numStates,C_REAL)*remainder
          f_shape%pointattributelist(2*(j-1)*SIZE(xPosArray)+2*i-1) = INT(remainder,C_INT) + 1
          f_shape%pointattributelist(2*(j-1)*SIZE(xPosArray)+2*i) = (j-1)*SIZE(xPosArray) + i
          f_shape%pointmarkerlist((j-1)*SIZE(xPosArray)+i) = 0
        END DO
      END DO
      f_shape%pointlist(SIZE(f_shape%pointlist)-7) = latticeData%xBounds(1)
      f_shape%pointlist(SIZE(f_shape%pointlist)-6) = latticeData%yBounds(1)
      f_shape%pointlist(SIZE(f_shape%pointlist)-5) = latticeData%xBounds(2)
      f_shape%pointlist(SIZE(f_shape%pointlist)-4) = latticeData%yBounds(1)
      f_shape%pointlist(SIZE(f_shape%pointlist)-3) = latticeData%xBounds(2)
      f_shape%pointlist(SIZE(f_shape%pointlist)-2) = latticeData%yBounds(2)
      f_shape%pointlist(SIZE(f_shape%pointlist)-1) = latticeData%xBounds(1)
      f_shape%pointlist(SIZE(f_shape%pointlist))   = latticeData%yBounds(2)
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-7) = 0
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-6) = SIZE(f_shape%pointattributelist)-3
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-5) = 0
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-4) = SIZE(f_shape%pointattributelist)-2
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-3) = 0
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-2) = SIZE(f_shape%pointattributelist)-1
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-1) = 0
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist))   = SIZE(f_shape%pointattributelist)
      f_shape%pointmarkerlist(SIZE(f_shape%pointmarkerlist)-3) = 1
      f_shape%pointmarkerlist(SIZE(f_shape%pointmarkerlist)-2) = 1
      f_shape%pointmarkerlist(SIZE(f_shape%pointmarkerlist)-1) = 1
      f_shape%pointmarkerlist(SIZE(f_shape%pointmarkerlist))   = 1

      f_shape%numberofsegments = 4
      CALL allocate_segments(f_shape,c_shape)
      f_shape%segmentlist(1) = SIZE(f_shape%pointmarkerlist)-3
      f_shape%segmentlist(2) = SIZE(f_shape%pointmarkerlist)-2
      f_shape%segmentlist(3) = SIZE(f_shape%pointmarkerlist)-2
      f_shape%segmentlist(4) = SIZE(f_shape%pointmarkerlist)-1
      f_shape%segmentlist(5) = SIZE(f_shape%pointmarkerlist)-1
      f_shape%segmentlist(6) = SIZE(f_shape%pointmarkerlist)
      f_shape%segmentlist(7) = SIZE(f_shape%pointmarkerlist)
      f_shape%segmentlist(8) = SIZE(f_shape%pointmarkerlist)-3
      f_shape%segmentmarkerlist(1) = 1
      f_shape%segmentmarkerlist(2) = 1
      f_shape%segmentmarkerlist(3) = 1
      f_shape%segmentmarkerlist(4) = 1

      DEALLOCATE(xPosArray)
      DEALLOCATE(yPosArray)

    ! Generate point data, random Q-state, points arranged in a hexagonal grid
    ! for use with the Q-state evolver and the Q-state based boundary finder
    ELSE IF (mapType == "hex") THEN
      CALL init_random_seed()
      num_x = INT((latticeData%xBounds(2) - latticeData%xBounds(1))/latticeData%gridSpacing(1))
      remainder = (latticeData%xBounds(2) - latticeData%xBounds(1)) - REAL(num_x)*latticeData%gridSpacing(1)
      IF (remainder < latticeData%gridSpacing(1)/2.0) THEN
        ALLOCATE(xPosArray(1:num_x+1),STAT=ios)
        IF ( ios /= 0) THEN
          WRITE(stderr,'(A,i16, A)') "ERROR : Allocation request denied. Array needed ",num_x+1," elements."
          STOP 1
        END IF
        xPosArray(1:num_x) = (/(latticeData%xBounds(1)+latticeData%gridSpacing(1)*REAL(i),i = 0,(num_x-1),1)/)
        xPosArray(num_x+1) = latticeData%xBounds(2)
      ELSE
        ALLOCATE(xPosArray(1:num_x+2),STAT=ios)
        IF ( ios /= 0) THEN
          WRITE(stderr,'(A,i16, A)') "ERROR : Allocation request denied. Array needed ",num_x+2," elements."
          STOP 1
        END IF
        xPosArray(1:num_x+1) = (/(latticeData%xBounds(1)+latticeData%gridSpacing(1)*REAL(i),i = 0,num_x,1)/)
        xPosArray(num_x+2) = latticeData%xBounds(2)
      END IF
      num_y = INT((latticeData%yBounds(2) - latticeData%yBounds(1))/latticeData%gridSpacing(2))
      remainder = (latticeData%yBounds(2) - latticeData%yBounds(1)) - REAL(num_y)*latticeData%gridSpacing(2)
      IF (remainder < latticeData%gridSpacing(2)/2.0) THEN
        ALLOCATE(yPosArray(1:num_y+1),STAT=ios)
        IF ( ios /= 0) THEN
          WRITE(stderr,'(A,i16, A)') "ERROR : Allocation request denied. Array needed ",num_y+1," elements."
          STOP 1
        END IF
        yPosArray(1:num_y) = (/(latticeData%yBounds(1)+latticeData%gridSpacing(2)*REAL(i),i = 0,(num_y-1),1)/)
        yPosArray(num_y+1) = latticeData%yBounds(2)
      ELSE
        ALLOCATE(yPosArray(1:num_y+2),STAT=ios)
        IF ( ios /= 0) THEN
          WRITE(stderr,'(A,i16, A)') "ERROR : Allocation request denied. Array needed ",num_y+2," elements."
          STOP 1
        END IF
        yPosArray(1:num_y+1) = (/(latticeData%yBounds(1)+latticeData%gridSpacing(2)*REAL(i),i = 0,num_y,1)/)
        yPosArray(num_y+2) = latticeData%yBounds(2)
      END IF
      f_shape%numberofpointattributes    = 2
      IF(MOD(SIZE(yPosArray),2) == 0 ) THEN
        f_shape%numberofpoints = (SIZE(xPosArray)-1)*SIZE(yPosArray) + SIZE(yPosArray)/2 + 4
      ELSE
        f_shape%numberofpoints = (SIZE(xPosArray)-1)*SIZE(yPosArray) + SIZE(yPosArray)/2 + 5
      END IF
      CALL allocate_points_ftoc(f_shape,c_shape)
      f_shape%pointlist(SIZE(f_shape%pointlist)-7) = latticeData%xBounds(1)
      f_shape%pointlist(SIZE(f_shape%pointlist)-6) = latticeData%yBounds(1)
      f_shape%pointlist(SIZE(f_shape%pointlist)-5) = latticeData%xBounds(2)
      f_shape%pointlist(SIZE(f_shape%pointlist)-4) = latticeData%yBounds(1)
      f_shape%pointlist(SIZE(f_shape%pointlist)-3) = latticeData%xBounds(2)
      f_shape%pointlist(SIZE(f_shape%pointlist)-2) = latticeData%yBounds(2)
      f_shape%pointlist(SIZE(f_shape%pointlist)-1) = latticeData%xBounds(1)
      f_shape%pointlist(SIZE(f_shape%pointlist))   = latticeData%yBounds(2)
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-7) = 0
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-6) = SIZE(f_shape%pointattributelist)-3
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-5) = 0
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-4) = SIZE(f_shape%pointattributelist)-2
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-3) = 0
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-2) = SIZE(f_shape%pointattributelist)-1
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist)-1) = 0
      f_shape%pointattributelist(SIZE(f_shape%pointattributelist))   = SIZE(f_shape%pointattributelist)
      f_shape%pointmarkerlist(SIZE(f_shape%pointmarkerlist)-3) = 1
      f_shape%pointmarkerlist(SIZE(f_shape%pointmarkerlist)-2) = 1
      f_shape%pointmarkerlist(SIZE(f_shape%pointmarkerlist)-1) = 1
      f_shape%pointmarkerlist(SIZE(f_shape%pointmarkerlist))   = 1

      f_shape%numberofsegments = 4
      CALL allocate_segments(f_shape,c_shape)
      f_shape%segmentlist(1) = SIZE(f_shape%pointmarkerlist)-3
      f_shape%segmentlist(2) = SIZE(f_shape%pointmarkerlist)-2
      f_shape%segmentlist(3) = SIZE(f_shape%pointmarkerlist)-2
      f_shape%segmentlist(4) = SIZE(f_shape%pointmarkerlist)-1
      f_shape%segmentlist(5) = SIZE(f_shape%pointmarkerlist)-1
      f_shape%segmentlist(6) = SIZE(f_shape%pointmarkerlist)
      f_shape%segmentlist(7) = SIZE(f_shape%pointmarkerlist)
      f_shape%segmentlist(8) = SIZE(f_shape%pointmarkerlist)-3
      f_shape%segmentmarkerlist(1) = 1
      f_shape%segmentmarkerlist(2) = 1
      f_shape%segmentmarkerlist(3) = 1
      f_shape%segmentmarkerlist(4) = 1

      k = 0
      l = 0
      DO j = 1, SIZE(yPosArray)
        DO i = 1,SIZE(xPosArray)
          IF(MOD(j,2) == 0) THEN
            IF(i == SIZE(xPosArray)) CYCLE
            k = k+1
            f_shape%pointlist(k) = xPosArray(i) + latticeData%gridSpacing(1)/2.0_C_REAL
            k = k+1
            f_shape%pointlist(k) = yPosArray(j)
            l = l+1
            CALL RANDOM_NUMBER(remainder)
            remainder = REAL(latticeData%numStates,C_REAL)*remainder
            f_shape%pointattributelist(2*l-1) = INT(remainder,C_INT) + 1
            f_shape%pointattributelist(2*l)   = l
            f_shape%pointmarkerlist(l) = 0
          ELSE
             k = k+1
             f_shape%pointlist(k) = xPosArray(i)
             k = k+1
             f_shape%pointlist(k) = yPosArray(j)
             l = l+1
             CALL RANDOM_NUMBER(remainder)
             remainder = REAL(latticeData%numStates,C_REAL)*remainder
             f_shape%pointattributelist(2*l-1) = INT(remainder,C_INT) + 1
             f_shape%pointattributelist(2*l)   = l
             f_shape%pointmarkerlist(l) = 0
          END IF
        END DO
      END DO

      DEALLOCATE(xPosArray)
      DEALLOCATE(yPosArray)

    ! Fail!!!
    ELSE
      WRITE(stderr,'(A,A,A)') "ERROR: Misread map type from file ",TRIM(filename),"."
      STOP 1
    END IF

    ! output the bare naked points and outer bounds to vtk for visualization/verification
    CALL output_triangle_to_vtk(f_shape,TRIM(outfileroot)//"_init.vtk")
    IF (mapType /= "datafile") THEN
      ! output the bare naked points and outer bounds to PSLG for later runs
      CALL output_triangle_to_poly(f_shape,TRIM(outfileroot)//"_init")
    END IF
    flags = TRIM("pqcOA")//C_NULL_CHAR
    CALL init_outputs(c_shape_out)
    CALL init_outputs(c_shape_vorout)
    ! triangulate initial points
    CALL ctriangulate(flags,c_shape,c_shape_out,c_shape_vorout)
    CALL copytriangles_c_to_f(c_shape_out,f_shape_out)
    ! output the first triangulation to vtk for visualization/verification
    CALL output_triangle_to_vtk(f_shape_out,TRIM(outfileroot)//"_000000.vtk")
    CALL output_triangle_to_poly(f_shape_out,TRIM(outfileroot)//"_000000")
    CALL deallocate_triangulateio(f_shape, input=.TRUE., neigh=.FALSE., voroni=.FALSE.)
    CALL deallocate_triangulateio(f_shape_out, input=.FALSE., neigh=.FALSE., voroni=.FALSE.)
    RETURN
100 WRITE(stderr,'(A,A,A)') "ERROR: Misread point data file path name  in ",TRIM(filename),"."
    STOP 1
150 WRITE(stderr,'(A,A,A)') "ERROR: Misread number of points from file ",TRIM(filename),"."
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

  SUBROUTINE init_random_seed
    INTEGER                             :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE  :: seed

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)

    DEALLOCATE(seed)
  END SUBROUTINE
END MODULE qmap_input
