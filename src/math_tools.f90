module kindprecision
  use iso_fortran_env
  implicit none

  ! default kinds
  integer,parameter           ::ikind=(kind(1))                 ! default range for integers
  integer,parameter           ::sikind = selected_int_kind(5)   ! short integers - range = +/- 10^5
  integer,parameter           ::likind = selected_int_kind(10)  ! long integers - range = +/- 10^10
  integer,parameter           ::fkind=(kind(0.e0))              ! default precision level for floats
  integer,parameter           ::dkind=(kind(0.d0))              ! default precision level for doubles

  ! custom kinds - set precision and range here using defaults, iso_fortran_env constants, or selected_real_kind intrinsic
  integer,parameter           :: SP = fkind                     ! custom single precision for float vars and constants
  integer,parameter           :: DP = dkind                     ! custom double precision for double vars and constants
  INTEGER, PARAMETER          :: StrBuffLen = 256

  ! selected_real_kind(p,r) will select the smallest kind of floating point storage that can store a float with
  ! precision of at least p decimal places and a range of at least +/- 10^r
  !
  ! These must be the kind parameter corresponding to MPI_DOUBLE_PRECISION
  ! and MPI_REAL
  !
  integer, parameter  :: psb_dpk_ = kind(1.d0)
  integer, parameter  :: psb_spk_ = kind(1.e0)

  integer, parameter :: RWP = psb_dpk_  ! real precision used thoughout the program quick switch
  integer, parameter :: IWP = ikind     ! integer range used in indices thoughout the program, quick switch incase integer width for indices needs to be increased
  integer, parameter :: KR4 = 4 !selected_real_kind(6,37)       ! single precision real
  integer, parameter :: KR8 = 8  !selected_real_kind(15,307)     ! double precision real
  integer, parameter :: KI4 = 4  !selected_int_kind(9)           ! single precision integer
  integer, parameter :: KI8 = 8  !selected_int_kind(18)          ! double precision integer

  integer, parameter :: METIS_IDX_T = 4    !this is to make portanolity with METIS library since they say
                                             ! in C, int can change with system but int32_t is consta
  integer, parameter :: METIS_REAL_T = 4    !this is to make portanolity with METIS library since they say
                                             ! in C, int can change with system but real_t is constant

end module kindprecision

module physical_constants
  use kindprecision
  implicit none
  REAL(RWP),PARAMETER                             :: k_B = 8.6173303d-5   ! Boltzmann constant in eV/K
end module physical_constants

module tools
  implicit none
contains
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

  SUBROUTINE find_neighbors(I,J,max_i,max_j,IP1,IM1,JP1,JM1)
    ! rows: assumes that min for I is 1, the max is max_i, and the universe is 0
    ! columns: assumes that min for J is 1 and the max is max_j, and the universe is 0
    USE kindprecision, only: IWP
    IMPLICIT NONE
    INTEGER(IWP),INTENT(IN)         :: I,J        ! indices of interest
    INTEGER(IWP),INTENT(IN)         :: max_i,max_j    ! max values for I and J
    INTEGER(IWP),INTENT(OUT)        :: IP1,IM1    ! neighbors for I, one step away
    INTEGER(IWP),INTENT(OUT)        :: JP1,JM1    ! neighbors for J, one step away

    ! I plus one, I minus one, and so on
    ip1 = i+1
    im1 = i-1
    jp1 = j+1
    jm1 = j-1
    ! ACCOUNT FOR NONPERIODIC BOUNDARIES
    IF(I+1 .GT. max_i) THEN
      IP1 = 0
    END IF
    IF(I-1 .LT. 1) THEN
      IM1 = 0
    END IF
    IF(J+1 .GT. max_j) THEN
      JP1 = 0
    END IF
    IF(J-1 .LT. 1) THEN
      JM1 = 0
    END IF
  END SUBROUTINE find_neighbors

  SUBROUTINE find_neighbors_periodic(I,J,max_i,max_j,IP1,IM1,JP1,JM1)
    ! rows: assumes that min for I is 1 and the max is max_i
    ! columns: assumes that min for J is 1 and the max is max_j
    USE kindprecision, only: IWP
    IMPLICIT NONE
    INTEGER(IWP),INTENT(IN)         :: I,J        ! indices of interest
    INTEGER(IWP),INTENT(IN)         :: max_i,max_j    ! max values for I and J
    INTEGER(IWP),INTENT(OUT)        :: IP1,IM1    ! neighbors for I, one step away
    INTEGER(IWP),INTENT(OUT)        :: JP1,JM1    ! neighbors for J, one step away

    ! I plus one, I minus one, and so on
    ip1 = i+1
    im1 = i-1
    jp1 = j+1
    jm1 = j-1
    ! ACCOUNT FOR PERIODIC BOUNDARIES
    IF(I+1 .GT. max_i) THEN
      IP1 = 1
    END IF
    IF(I-1 .LT. 1) THEN
      IM1 = max_i
    END IF
    IF(J+1 .GT. max_j) THEN
      JP1 = 1
    END IF
    IF(J-1 .LT. 1) THEN
      JM1 = max_j
    END IF
  END SUBROUTINE find_neighbors_periodic

  SUBROUTINE find_neighbors_column_periodic(I,J,max_i,max_j,IP1,IM1,JP1,JM1)
    ! columns: assumes that min for I is 1 and the max is max_i
    ! rows: assumes that min for J is 1 and the max is max_j, and the universe is 0
    USE kindprecision, only: IWP
    IMPLICIT NONE
    INTEGER(IWP),INTENT(IN)         :: I,J        ! indices of interest
    INTEGER(IWP),INTENT(IN)         :: max_i,max_j    ! max values for I and J
    INTEGER(IWP),INTENT(OUT)        :: IP1,IM1    ! neighbors for I, one step away
    INTEGER(IWP),INTENT(OUT)        :: JP1,JM1    ! neighbors for J, one step away

    ! I plus one, I minus one, and so on
    ip1 = i+1
    im1 = i-1
    jp1 = j+1
    jm1 = j-1
    ! ACCOUNT FOR PERIODIC BOUNDARIES IN X (I) DIRECTION
    IF(I+1 .GT. max_i) THEN
      IP1 = 1
    END IF
    IF(I-1 .LT. 1) THEN
      IM1 = max_i
    END IF
    ! ACCOUNT FOR NONPERIODIC BOUNDARIES IN Y (J) DIRECTION
    IF(J+1 .GT. max_j) THEN
      JP1 = 0
    END IF
    IF(J-1 .LT. 1) THEN
      JM1 = 0
    END IF
  END SUBROUTINE find_neighbors_column_periodic

  SUBROUTINE find_neighbors_row_periodic(I,J,max_i,max_j,IP1,IM1,JP1,JM1)
    ! columns: assumes that min for I is 1 and the max is max_i, and the universe is 0
    ! rows: assumes that min for J is 1 and the max is max_j
    USE kindprecision, only: IWP
    IMPLICIT NONE
    INTEGER(IWP),INTENT(IN)         :: I,J        ! indices of interest
    INTEGER(IWP),INTENT(IN)         :: max_i,max_j    ! max values for I and J
    INTEGER(IWP),INTENT(OUT)        :: IP1,IM1    ! neighbors for I, one step away
    INTEGER(IWP),INTENT(OUT)        :: JP1,JM1    ! neighbors for J, one step away

    ! I plus one, I minus one, and so on
    ip1 = i+1
    im1 = i-1
    jp1 = j+1
    jm1 = j-1
    ! ACCOUNT FOR NONPERIODIC BOUNDARIES IN X (I) DIRECTION
    IF(I+1 .GT. max_i) THEN
      IP1 = 0
    END IF
    IF(I-1 .LT. 1) THEN
      IM1 = 0
    END IF
    ! ACCOUNT FOR PERIODIC BOUNDARIES IN Y (J) DIRECTION
    IF(J+1 .GT. max_j) THEN
      JP1 = 1
    END IF
    IF(J-1 .LT. 1) THEN
      JM1 = max_j
    END IF
  END SUBROUTINE find_neighbors_row_periodic
end module tools
