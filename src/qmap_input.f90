! Tomas Mondragon
! Computer Scientist
! USACE ERDC Information Technology Laboratory
! Computational Analysis Branch
! Tomas.A.Mondragon@erdc.dren.mil
!   delivered 06 Feb 2016
MODULE qmap_input
  IMPLICIT NONE

CONTAINS

  ! Read input from Q-map dat file into a lattice data structure
  ! USAGE: lattice_data = read_Qmap(filename)
  FUNCTION read_Qmap(filename) RESULT(lattice_data)
    USE filehandling, only: safeopen_readonly, stderr
    USE basictypes, only: lattice
    IMPLICIT NONE

    CHARACTER(LEN=256),INTENT(IN)                     :: filename
    TYPE(lattice)                                     :: lattice_data
    INTEGER,DIMENSION(:),ALLOCATABLE                  :: Q
    INTEGER                                           :: fd,ios
    INTEGER                                           :: is_grid
    INTEGER                                           :: num_cols,num_rows,col_i,row_j
    LOGICAL                                           :: qstate_present

    fd = safeopen_readonly(filename)
    CALL skip_comments_Qmap(fd)
    ! 0 =: aperiodic, 1 =: periodic in y, 2 =: periodic in x, 3 =: periodic everywhere
    READ(fd,*,iostat=ios,err=900) lattice_data%periodicity
    CALL skip_comments_Qmap(fd)
    ! temperature of crystal during evolution (useful during Potts generation)
    READ(fd,*,iostat=ios,err=600) lattice_data%temperature
    CALL skip_comments_Qmap(fd)
    ! constant used in grain boundary energy calculations (electronVolts per ?m)
    READ(fd,*,iostat=ios,err=700) lattice_data%en_per_length
    CALL skip_comments_Qmap(fd)
    ! number of possible discrete states
    READ(fd,*,iostat=ios,err=500) lattice_data%num_states
    CALL skip_comments_Qmap(fd)
    ! whether or not the qmap data is presented in a predefined grid or not
    READ(fd,*,iostat=ios) is_grid
    IF(is_grid .ne. 0) THEN
      ! Data is presented in a grid - nodes are only denoted by their state. Ask for x and y gridspacing
      CALL skip_comments_Qmap(fd)
      READ(fd,*,iostat=ios,err=200) lattice_data%x_bounds
      CALL skip_comments_Qmap(fd)
      READ(fd,*,iostat=ios,err=250) lattice_data%y_bounds
      CALL skip_comments_Qmap(fd)
      READ(fd,*,iostat=ios,err=300) lattice_data%grid_spacing
      num_cols = FLOOR((lattice_data%x_bounds(2) - lattice_data%x_bounds(1))/lattice_data%grid_spacing(1))
      num_rows = FLOOR((lattice_data%y_bounds(2) - lattice_data%y_bounds(1))/lattice_data%grid_spacing(2))
      IF(lattice_data%x_bounds(2) .gt. lattice_data%x_bounds(1) + lattice_data%grid_spacing(1) * REAL(num_cols)) THEN
        WRITE(stderr,'(A,F9.4,A,F9.4)') "Warning: Grid X upper bound reduced from ",lattice_data%x_bounds(2)," to", &
                     lattice_data%x_bounds(1) + lattice_data%grid_spacing(1) * REAL(num_cols)
        lattice_data%x_bounds(2) = lattice_data%x_bounds(1) + lattice_data%grid_spacing(1) * REAL(num_cols)
      ENDIF
      IF(lattice_data%y_bounds(2) .gt. lattice_data%y_bounds(1) + lattice_data%grid_spacing(2) * REAL(num_rows)) THEN
        WRITE(stderr,'(A,F9.4,A,F9.4)') "Warning: Grid Y upper bound reduced from ",lattice_data%y_bounds(2)," to", &
                     lattice_data%y_bounds(1) + lattice_data%grid_spacing(2) * REAL(num_rows)
        lattice_data%y_bounds(2) = lattice_data%y_bounds(1) + lattice_data%grid_spacing(2) * REAL(num_rows)
      ENDIF
      WRITE(stderr,'(A)') "Grid setup:"
      WRITE(stderr,'(A,F9.4,A,F9.4)') "  x range:",lattice_data%x_bounds(1)," - ",lattice_data%x_bounds(2)
      WRITE(stderr,'(A,F9.4,A,I4)') "  x spacing:",lattice_data%grid_spacing(1)," x levels: ",num_cols
      WRITE(stderr,'(A,F9.4,A,F9.4)') "  y range:",lattice_data%y_bounds(1)," - ",lattice_data%y_bounds(2)
      WRITE(stderr,'(A,F9.4,A,I4)') "  y spacing:",lattice_data%grid_spacing(2)," y levels: ",num_rows
      CALL skip_comments_Qmap(fd)
      ! Number of sweeps used to generate (evolve) this grid. Used only if generating sample data
      READ(fd,*,iostat=ios,err=400) lattice_data%mc_sweeps
      ! now allocate the node array
      ALLOCATE(lattice_data%nodes(1:num_cols*num_rows))
      ALLOCATE(Q(1:num_cols*num_rows))
      lattice_data%num_nodes = num_cols*num_rows
      CALL skip_comments_Qmap(fd)
      READ(fd,*,iostat=ios,err=100) Q
      DO row_j=0,num_rows-1
        DO col_i=1,num_cols
          lattice_data%nodes(col_i+row_j*num_cols)%q_state = Q(col_i+row_j*num_cols)
          lattice_data%nodes(col_i+row_j*num_cols)%location(1) = lattice_data%x_bounds(1)+lattice_data%grid_spacing(1) &
                                                                  *REAL(col_i)
          lattice_data%nodes(col_i+row_j*num_cols)%location(2) = lattice_data%y_bounds(1)+lattice_data%grid_spacing(2) &
                                                                  *REAL(row_j+1)
          lattice_data%nodes(col_i+row_j*num_cols)%is_universe = .false.
          lattice_data%nodes(col_i+row_j*num_cols)%n_index = col_i+row_j*num_cols
        END DO
      END DO
      DEALLOCATE(Q)
    ELSE
      ! Data is freeform - nodes are listed with a location and a state(maybe)
      CALL skip_comments_Qmap(fd)
      ! are the q-states already in the input?
      READ(fd,*,iostat=ios) is_grid
      qstate_present = (is_grid .eq. 1)
      is_grid = 0
      CALL skip_comments_Qmap(fd)
      ! number of nodes
      READ(fd,*,iostat=ios,err=150) lattice_data%num_nodes
      ! now allocate the node array
      ALLOCATE(lattice_data%nodes(1:lattice_data%num_nodes))
      CALL skip_comments_Qmap(fd)
      DO row_j = 1,lattice_data%num_nodes
        IF(qstate_present) THEN
          READ(fd,*,iostat=ios,err=100) lattice_data%nodes(row_j)%location(1), lattice_data%nodes(row_j)%location(2), &
            lattice_data%nodes(row_j)%q_state
        ELSE
          READ(fd,*,iostat=ios,err=100) lattice_data%nodes(row_j)%location(1), lattice_data%nodes(row_j)%location(2)
        END IF
      END DO
    END IF
    RETURN
100 WRITE(stderr,'(A,A,A)') "End of file ",TRIM(filename)," reached before all values of Q-map read."
    WRITE(stderr,'(A,A,A)') "Dimensions specified in ",TRIM(filename)," are wrong."
    STOP 1
150 WRITE(stderr,'(A,A,A)') "Error reading number of nodes from file ",TRIM(filename),"."
    STOP 1
200 WRITE(stderr,'(A,A,A)') "Error reading x-bounds from file ",TRIM(filename),"."
    STOP 1
250 WRITE(stderr,'(A,A,A)') "Error reading y-bounds from file ",TRIM(filename),"."
    STOP 1
300 WRITE(stderr,'(A,A,A)') "Error reading grid spacing information from file ",TRIM(filename),"."
    STOP 1
400 WRITE(stderr,'(A,A,A)') "Error reading number of Evolver Monte Carlo Steps from file ",TRIM(filename),"."
    STOP 1
500 WRITE(stderr,'(A,A,A)') "Error reading number of states from file ",TRIM(filename),"."
    STOP 1
600 WRITE(stderr,'(A,A,A)') "Error reading temperature from file ",TRIM(filename),"."
    STOP 1
700 WRITE(stderr,'(A,A,A)') "Error reading energy per unit length from file ",TRIM(filename),"."
    STOP 1
900 WRITE(stderr,'(A,A,A)') "Error reading periodicity from file ",TRIM(filename),"."
    STOP 1
  END FUNCTION read_Qmap


  SUBROUTINE skip_comments_Qmap(fd)
    IMPLICIT NONE
    INTEGER, INTENT(IN)                               :: fd   ! FILE DESCRIPTOR NUMBER
    CHARACTER(LEN=256)                                :: buffer

    ReadComments: DO
      READ(fd, '(A)', END=100) buffer
      IF((buffer(1:1) /= "#") .AND. (TRIM(buffer) /= "")) exit ReadComments
    END DO ReadComments
    BACKSPACE(fd)
100 CONTINUE
  END SUBROUTINE skip_comments_Qmap
END MODULE qmap_input
