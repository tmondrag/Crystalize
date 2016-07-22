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
    INTEGER                                           :: fd,ios

    fd = safeopen_readonly(filename)
    CALL skip_comments_Qmap(fd)
    READ(fd,*,iostat=ios,err=900) lattice_data%periodicity
    CALL skip_comments(fd)
    READ(fd,*,iostat=ios,err=600) lattice_data%temperature
    CALL skip_comments(fd)
    READ(fd,*,iostat=ios,err=700) lattice_data%en_per_length
    CALL skip_comments(fd)
    READ(fd,*,iostat=ios,err=500) lattice_data%num_states
    RETURN
100 WRITE(stderr,'(A,A,A)') "End of file ",TRIM(filename)," reached before all values of Q-map read."
    WRITE(stderr,'(A,A,A)') "Dimensions specified in ",TRIM(filename)," are wrong."
    STOP 1
200 WRITE(stderr,'(A,A,A)') "Error reading number of columns from file ",TRIM(filename),"."
    STOP 1
300 WRITE(stderr,'(A,A,A)') "Error reading number of rows from file ",TRIM(filename),"."
    STOP 1
400 WRITE(stderr,'(A,A,A)') "Error reading number of Evolver Monte Carlo Steps from file ",TRIM(filename),"."
    STOP 1
500 WRITE(stderr,'(A,A,A)') "Error reading number of states from file ",TRIM(filename),"."
    STOP 1
600 WRITE(stderr,'(A,A,A)') "Error reading temperature from file ",TRIM(filename),"."
    STOP 1
700 WRITE(stderr,'(A,A,A)') "Error reading energy per unit length from file ",TRIM(filename),"."
    STOP 1
800 WRITE(stderr,'(A,A,A)') "Error reading length scale from file ",TRIM(filename),"."
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
