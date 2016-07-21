MODULE triangle_input
CONTAINS
  ! Read in data for the triangle.c library. Starting file can be either a *.node or *.poly file
  SUBROUTINE read_shapes(filename,c_shape,f_shape)
    USE triangle_c_wrap, only: triangulateio,f_triangulateio
    USE filehandling, only: stderr
    IMPLICIT NONE
    TYPE(triangulateio),INTENT(OUT)   :: c_shape
    TYPE(f_triangulateio),INTENT(OUT) :: f_shape
    CHARACTER(len=*),INTENT(IN)       :: filename
    CHARACTER(len=256)                :: filename_poly, filename_node, filename_ele, filename_root
    INTEGER                           :: Cind
    LOGICAL                           :: exists

    ! find out whether the input file is a *.node or a *.poly or if input string is just the root
    Cind = INDEX(filename,".poly",back=.TRUE.)
    IF(Cind .EQ. 0) THEN
      Cind = INDEX(filename,".node",back=.TRUE.)
      IF(Cind .EQ. 0) THEN
        filename_root = filename
        filename_poly = TRIM(filename_root)//"poly"
        filename_node = TRIM(filename_root)//"node"
        filename_ele  = TRIM(filename_root)//"ele"
      ELSE
        filename_node = filename
        filename_root = filename(1:Cind)
        filename_poly = TRIM(filename_root)//"poly"
        filename_ele  = TRIM(filename_root)//"ele"
      END IF
    ELSE
      filename_poly = filename
      filename_root = filename(1:Cind)
      filename_node = TRIM(filename_root)//"node"
      filename_ele  = TRIM(filename_root)//"ele"
    ENDIF

    f_shape%numberofpoints             = 0
    f_shape%numberofpointattributes    = 0
    f_shape%numberofholes              = 0
    f_shape%numberofregions            = 0
    f_shape%numberofsegments           = 0
    f_shape%numberoftriangles          = 0
    f_shape%numberofcorners            = 0
    f_shape%numberoftriangleattributes = 0
    f_shape%numberofedges              = 0
    ! at the very least, a *.node file or *.poly should exist but if the *.poly file exists, read that first

    INQUIRE(FILE = filename_poly, EXIST = exists)
    IF(exists) THEN
      CALL read_polyfile(filename_poly,filename_node,c_shape,f_shape)
    ELSE
      INQUIRE(FILE = filename_node, EXIST = exists)
      IF(exists) THEN
        CALL read_nodefile(filename_node,c_shape,f_shape)
      ELSE
        WRITE(stderr,'(a,a,a,a,a)') "ERROR: Neither ",TRIM(filename_poly), " nor ", TRIM(filename_node), " exist."
        STOP 1
      END IF
    END IF
    INQUIRE(FILE = filename_ele, EXIST = exists)
    IF(exists) THEN
      CALL read_elefile(filename_ele,c_shape,f_shape)
    ELSE
      WRITE(stderr,'(a,a,a)') "WARNING: ",TRIM(filename_ele), " does not exist."
    END IF
  END SUBROUTINE read_shapes

  SUBROUTINE read_polyfile(filename_poly,filename_node,c_shape,f_shape)
    USE ISO_C_BINDING, only: C_INT
    USE triangle_c_wrap, only: triangulateio,f_triangulateio, allocate_holes,allocate_points_ftoc,allocate_segments,allocate_regions
    USE filehandling, only:stderr,safeopen_readonly
    IMPLICIT NONE
    TYPE(triangulateio),INTENT(OUT)    :: c_shape
    TYPE(f_triangulateio),INTENT(OUT)  :: f_shape
    CHARACTER(len=*),INTENT(IN)        :: filename_poly, filename_node
    INTEGER                            :: pfd        ! file descriptor/IO unit
    INTEGER(kind=C_INT)                :: num_vertices, num_dimensions, num_attributes, num_bmarkers, i
    INTEGER(kind=C_INT)                :: num_segments, num_holes, num_regions
    LOGICAL                            :: has_bound_flag
    CHARACTER(len=1)                   :: dummy

    pfd = safeopen_readonly(filename_poly)

    ! Read vertices
    CALL skip_comments_tri(pfd)
    READ(pfd,*,ERR=100) num_vertices, num_dimensions, num_attributes, num_bmarkers
    IF(num_vertices == 0) THEN
      CALL read_nodefile(filename_node,c_shape,f_shape)
    ELSE
      has_bound_flag = .TRUE.
      IF (num_bmarkers .EQ. 0) has_bound_flag = .FALSE.
      f_shape%numberofpoints             = num_vertices
      f_shape%numberofpointattributes    = num_attributes

      CALL allocate_points_ftoc(f_shape,c_shape)
      DO i=1,num_vertices
        CALL skip_comments_tri(pfd)
        CALL read_vertex(pfd,i,f_shape,has_bound_flag)
      END DO
    END IF

    ! Read segments
    CALL skip_comments_tri(pfd)
    READ(pfd,*,ERR=110) num_segments, num_bmarkers
    has_bound_flag = .TRUE.
    IF (num_bmarkers .EQ. 0) has_bound_flag = .FALSE.
    f_shape%numberofsegments           = num_segments
    CALL allocate_segments(f_shape,c_shape)
    DO i = 1,num_segments
      CALL skip_comments_tri(pfd)
      CALL read_segment(pfd,i,f_shape,has_bound_flag)
    END DO

    ! Read holes
    CALL skip_comments_tri(pfd)
    READ(pfd,*,ERR=120) num_holes
    has_bound_flag = .TRUE.
    IF (num_bmarkers .EQ. 0) has_bound_flag = .FALSE.
    f_shape%numberofholes           = num_holes
    CALL allocate_holes(f_shape,c_shape)
    DO i = 1,num_holes
      CALL skip_comments_tri(pfd)
      CALL read_hole(pfd,i,f_shape)
    END DO

    ! premptively set num_regions to 0
    c_shape%numberofregions = 0
    ! detects EOF or read regions
    CALL skip_comments_tri(pfd)
    READ(pfd,*,END=90,ERR=90) dummy ! quick exit if there are no regions specified
    BACKSPACE(pfd)
    READ(pfd,*,ERR=130) num_regions
    f_shape%numberofregions = num_regions
    CALL allocate_regions(f_shape,c_shape)
    DO i = 1,num_regions
      CALL skip_comments_tri(pfd)
      CALL read_region(pfd,i,f_shape)
    END DO
 90 CONTINUE
    RETURN
100 WRITE(stderr,'(A,A,A)') "ERROR: Mistake reading vertices from .poly file ",filename_poly,"."
    STOP
110 WRITE(stderr,'(A,A,A)') "ERROR: Mistake reading segments from .poly file ",filename_poly,"."
    STOP
120 WRITE(stderr,'(A,A,A)') "ERROR: Mistake reading holes from .poly file ",filename_poly,"."
    STOP
130 WRITE(stderr,'(A,A,A)') "ERROR: Mistake reading regions from .poly file ",filename_poly,"."
    STOP
  END SUBROUTINE read_polyfile

  SUBROUTINE read_nodefile(filename,c_shape,f_shape)
    USE ISO_C_BINDING, only: C_INT
    USE triangle_c_wrap, only: triangulateio, f_triangulateio, allocate_points_ftoc
    USE filehandling, only: stderr,safeopen_readonly
    IMPLICIT NONE
    TYPE(triangulateio),INTENT(OUT)   :: c_shape
    TYPE(f_triangulateio),INTENT(OUT) :: f_shape
    CHARACTER(len=*),INTENT(IN)       :: filename
    INTEGER                           :: nfd        ! file descriptor/IO unit
    INTEGER(kind=C_INT)               :: num_vertices, num_dimensions, num_attributes, num_bmarkers, i
    LOGICAL                           :: has_bound_flag

    nfd = safeopen_readonly(filename)
    CALL skip_comments_tri(nfd)
    READ(nfd,*,ERR=150) num_vertices, num_dimensions, num_attributes, num_bmarkers
    has_bound_flag = .TRUE.
    IF (num_bmarkers .EQ. 0) has_bound_flag = .FALSE.
    f_shape%numberofpoints             = num_vertices
    f_shape%numberofpointattributes    = num_attributes

    CALL allocate_points_ftoc(f_shape,c_shape)
    DO i=1,num_vertices
      CALL skip_comments_tri(nfd)
      CALL read_vertex(nfd,i,f_shape,has_bound_flag)
    END DO
    CLOSE(nfd)
    RETURN

150 WRITE(stderr,'(A,A,A)') "ERROR: Mistake reading from .node file ",filename,"."
    STOP
  END SUBROUTINE read_nodefile

  ! Read if triangle descriptions from file filename into f_shape
  !
  ! USAGE: read_elefile(filename,c_shape,f_shape)
  SUBROUTINE read_elefile(filename,c_shape,f_shape)
    USE ISO_C_BINDING, only: C_INT
    USE triangle_c_wrap, only: triangulateio, f_triangulateio, allocate_triangles
    USE filehandling, only: stderr,safeopen_readonly
    IMPLICIT NONE
    TYPE(triangulateio),INTENT(OUT)   :: c_shape
    TYPE(f_triangulateio),INTENT(OUT) :: f_shape
    CHARACTER(len=*),INTENT(IN)       :: filename
    INTEGER                           :: efd        ! file descriptor/IO unit
    INTEGER                           :: i
    INTEGER(kind=C_INT)               :: num_tri, nodes_per_tri, num_attributes

    efd = safeopen_readonly(filename)
    CALL skip_comments_tri(efd)
    READ(efd,*,ERR=175) num_tri, nodes_per_tri, num_attributes
    f_shape%numberoftriangles          = num_tri
    f_shape%numberofcorners            = nodes_per_tri
    f_shape%numberoftriangleattributes = num_attributes

    CALL allocate_triangles(f_shape,c_shape)

    DO i=1,num_tri
      CALL skip_comments_tri(efd)
      CALL read_triangle(efd,i,f_shape)
    END DO
    CLOSE(efd)
    RETURN

175 WRITE(stderr,'(A,A,A)') "ERROR: Mistake reading from .ele file ",filename,"."
    STOP
  END SUBROUTINE read_elefile

  ! Read a triangle description line from IO unit fd into the appropriate place in f_shape
  !
  ! USAGE: read_triangle(fd,tri_index,f_shape)
  SUBROUTINE read_triangle(fd,tri_index,f_shape)
    USE ISO_C_BINDING, only: C_INT
    USE kindprecision, only: StrBuffLen
    USE filehandling, only: stderr
    USE triangle_c_wrap, only: f_triangulateio
    IMPLICIT NONE
    TYPE(f_triangulateio),INTENT(INOUT)               :: f_shape
    INTEGER,INTENT(IN)                                :: tri_index
    INTEGER(kind=C_INT),INTENT(IN)                    :: fd   ! File descriptor/Unit number
    INTEGER(kind=C_INT)                               :: check_index
    INTEGER(kind=C_INT)                               :: nodes_per_tri
    INTEGER(kind=C_INT)                               :: num_attributes
    CHARACTER(LEN=StrBuffLen)                         :: buffer

    nodes_per_tri = f_shape%numberofcorners
    num_attributes = f_shape%numberoftriangleattributes

    IF(nodes_per_tri == 3) THEN
      READ(fd,*,ERR=100) check_index,f_shape%trianglelist(tri_index*3-2:tri_index*3),buffer
    ELSE IF(nodes_per_tri == 6) THEN
      READ(fd,*,ERR=100) check_index,f_shape%trianglelist(tri_index*6-5:tri_index*6),buffer
    ELSE
      GOTO 100
    END IF
    ! check_index is unused but it should match tri_index. Do a sanity check if you wish.
    ! what follows is a loop for the attributes
    IF(num_attributes .GT. 0) THEN
      READ(buffer,*,ERR=100) f_shape%triangleattributelist((tri_index-1)*num_attributes+1:tri_index*num_attributes)
    END IF
    RETURN
100 WRITE(stderr,'(A)') "ERROR: triangle information misread."
    STOP
  END SUBROUTINE read_triangle

  SUBROUTINE read_vertex(fd,vert_index, f_shape, has_bound_flag)
    USE ISO_C_BINDING, only: C_INT
    USE kindprecision, only: StrBuffLen
    USE filehandling, only: stderr
    USE triangle_c_wrap, only: f_triangulateio
    IMPLICIT NONE
    TYPE(f_triangulateio),INTENT(INOUT)               :: f_shape
    INTEGER,INTENT(IN)                                :: vert_index
    INTEGER(kind=C_INT),INTENT(IN)                    :: fd   ! File descriptor/Unit number
    INTEGER(kind=C_INT)                               :: check_index
    INTEGER(kind=C_INT)                               :: num_attributes
    LOGICAL,INTENT(IN)                                :: has_bound_flag
    CHARACTER(LEN=StrBuffLen)                         :: buffer

    num_attributes = f_shape%numberofpointattributes
    IF((num_attributes .GT. 0) .or. has_bound_flag) THEN
    READ(fd,*,ERR=100) check_index,f_shape%pointlist(vert_index*2-1:vert_index*2),buffer
    ELSE
    READ(fd,*,ERR=100) check_index,f_shape%pointlist(vert_index*2-1:vert_index*2)
    END IF
    ! check_index is unused but it should match vert_index. Do a sanity check if you wish.
    ! what follows is a loop for the attributes
    IF(num_attributes .GT. 0) THEN
      READ(buffer,*,ERR=100) f_shape%pointattributelist((vert_index-1)*num_attributes+1:vert_index*num_attributes),buffer
    END IF
    IF(has_bound_flag) THEN
      READ(buffer,*,ERR=100) f_shape%pointmarkerlist(vert_index)
    END IF
    RETURN
100 WRITE(stderr,'(A)') "ERROR: vertex information misread."
    STOP
  END SUBROUTINE read_vertex

  SUBROUTINE read_segment(fd,seg_index, f_shape, has_bound_flag)
    USE ISO_C_BINDING, only: C_INT
    USE kindprecision, only: StrBuffLen
    USE filehandling, only: stderr
    USE triangle_c_wrap, only: f_triangulateio
    IMPLICIT NONE
    TYPE(f_triangulateio),INTENT(INOUT)               :: f_shape
    INTEGER,INTENT(IN)                                :: seg_index
    INTEGER(kind=C_INT),INTENT(IN)                    :: fd   ! File descriptor/Unit number
    INTEGER(kind=C_INT)                               :: check_index
    LOGICAL,INTENT(IN)                                :: has_bound_flag
    CHARACTER(LEN=StrBuffLen)                         :: buffer

    READ(fd,*,ERR=100) check_index,f_shape%segmentlist(seg_index*2-1:seg_index*2),buffer
    ! check_index is unused but it should match vert_index. Do a sanity check if you wish.
    IF(has_bound_flag) THEN
      READ(buffer,*,ERR=100) f_shape%segmentmarkerlist(seg_index)
    END IF
    RETURN
100 WRITE(stderr,'(A)') "ERROR: segment information misread."
    STOP
  END SUBROUTINE read_segment

  SUBROUTINE read_hole(fd,hole_index, f_shape)
    USE ISO_C_BINDING, only: C_INT
    USE filehandling, only: stderr
    USE triangle_c_wrap, only: f_triangulateio
    IMPLICIT NONE
    TYPE(f_triangulateio),INTENT(INOUT)               :: f_shape
    INTEGER,INTENT(IN)                                :: hole_index
    INTEGER(kind=C_INT),INTENT(IN)                    :: fd   ! File descriptor/Unit number
    INTEGER(kind=C_INT)                               :: check_index

    READ(fd,*,ERR=100) check_index,f_shape%holelist(hole_index*2-1:hole_index*2)
    ! check_index is unused but it should match vert_index. Do a sanity check if you wish.
    RETURN
100 WRITE(stderr,'(A)') "ERROR: hole information misread."
    STOP
  END SUBROUTINE read_hole

  SUBROUTINE read_region(fd,reg_index, f_shape)
    USE ISO_C_BINDING, only: C_INT
    USE filehandling, only: stderr
    USE triangle_c_wrap, only: f_triangulateio,C_REAL
    IMPLICIT NONE
    TYPE(f_triangulateio),INTENT(INOUT)               :: f_shape
    INTEGER,INTENT(IN)                                :: reg_index
    INTEGER(kind=C_INT),INTENT(IN)                    :: fd   ! File descriptor/Unit number
    REAL(kind=C_REAL),DIMENSION(1:5)                  :: buff_array
    CHARACTER(len=128)                                :: buff_string

    READ(fd,'(A)') buff_string
    buff_array = 0
    READ(buff_string,*,ERR=100,END=50) buff_array(1:4),buff_string
    READ(buff_string,*,ERR=100) buff_array(5)
 50 IF( INT(buff_array(1)) .NE. reg_index) PRINT *, "OOPS! ",INT(buff_array(1)),"!=",reg_index 
    f_shape%regionlist(reg_index*4-3:reg_index*4) = buff_array(2:5)
    !PRINT *, f_shape%regionlist(reg_index*4-3:reg_index*4)
    RETURN
100 WRITE(stderr,'(A)') "ERROR: region information misread."
    STOP
  END SUBROUTINE read_region

  SUBROUTINE skip_comments_tri(fd)
    IMPLICIT NONE
    INTEGER, INTENT(IN)                               :: fd   ! FILE DESCRIPTOR NUMBER
    CHARACTER(LEN=256)                                :: buffer

    ReadComments: DO
      READ(fd, '(A)', END=100) buffer
      IF((buffer(1:1) /= "#") .AND. (TRIM(buffer) /= "")) exit ReadComments
    END DO ReadComments
    BACKSPACE(fd)
100 CONTINUE
  END SUBROUTINE skip_comments_tri
END MODULE triangle_input
