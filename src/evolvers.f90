! Tomas Mondragon
! Computer Scientist
! USACE ERDC Information Technology Laboratory
! Computational Analysis Branch
! Tomas.A.Mondragon@erdc.dren.mil
!   delivered 16 Apr 2017

MODULE evolvers
  USE, INTRINSIC :: ISO_C_BINDING, only: C_INT,C_BOOL
  USE            :: triangle_c_wrap, only: C_REAL


CONTAINS
  SUBROUTINE evolve_qstates(inLattice)
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: stderr -> ERROR_UNIT, stdout -> OUTPUT_UNIT
    USE basictypes, only: Lattice
    IMPLICIT NONE
    TYPE(Lattice),INTENT(INOUT)                       :: inLattice
    TYPE(LatticeEdgeSHead),DIMENSION(:),ALLOCATABLE   :: edgeHash
    INTEGER                                           :: i,j,k
    REAL                                              :: tmp

    edgeHash = build_edge_hash(inLattice)
    DO i = 1, inLattice%mcSweeps
      DO j = 1, inLattice%numVertices
        CALL RANDOM_NUMBER(tmp)
        k = CEILING(inLattice%numVertices*tmp)
        CALL flip_qstate_isometric(inLattice,edgeHash,k)
      END DO
      IF(MOD(i,10) == 0) THEN
        WRITE(stdout,'(A,I8,A)') "Monte-Carlo step ",i," elapsed." ! print out to entertain the user during long process
      END IF
    END DO
    WRITE(stdout,'(A,I8,A)') "Q-map generated after ", i, " Monte-Carlo steps."
  END SUBROUTINE evolve_qstates

  SUBROUTINE evolve_positions(inLattice)
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: stderr -> ERROR_UNIT, stdout -> OUTPUT_UNIT
    USE basictypes, only: Lattice
    IMPLICIT NONE
    TYPE(Lattice),INTENT(INOUT)                       :: inLattice
    TYPE(LatticeEdgeSHead),DIMENSION(:),ALLOCATABLE   :: edgeHash
    INTEGER                                           :: i,j,k
    REAL                                              :: tmp

    edgeHash = build_edge_hash(inLattice)
    DO i = 1, inLattice%mcSweeps
      DO j = 1, inLattice%numVertices
        CALL RANDOM_NUMBER(tmp)
        k = CEILING(inLattice%numVertices*tmp)
        CALL flip_position_isometric(inLattice,edgeHash,k)
      END DO
      IF(MOD(i,10) == 0) THEN
        WRITE(stdout,'(A,I8,A)') "Monte-Carlo step ",i," elapsed." ! print out to entertain the user during long process
      END IF
    END DO
    WRITE(stdout,'(A,I8,A)') "Position map generated after ", i, " Monte-Carlo steps."
  END SUBROUTINE evolve_positions

  SUBROUTINE flip_qstate_isometric(inLattice,edgeHash,vertexIndex)
    USE basictypes, only: Lattice
    IMPLICIT NONE
    TYPE(Lattice),INTENT(INOUT)                       :: inLattice
    TYPE(LatticeEdgeSHead),DIMENSION(:),INTENT(IN)    :: edgeHash
    INTEGER,INTENT(IN)                                :: vertexIndex
    REAL                                              :: enrg1,enrg2
    INTEGER                                           :: presentState, randomState
    INTEGER                                           :: numNeighbors,i,edgeIndex,neighborIndex
    INTEGER,DIMENSION(:)                              :: neighborStates
    REAL,DIMENSION(:)                                 :: neighborDistances
    TYPE(LatticeEdgeSItem),POINTER                    :: curr

    presentState = inLattice%vertices(vertexIndex)%qState
    curr => edgeHash(vertexIndex)%listHead
    numNeighbors = 0
    DO
      curr => curr%next
      IF (curr%isHead) EXIT
      numNeighbors = numNeighbors + 1
    END DO
    ALLOCATE(neighborStates(numNeighbors),neighborDistances(numNeighbors))
    curr => edgeHash(vertexIndex)%listHead
    i = 0
    DO
      curr => curr%next
      IF (curr%isHead) EXIT
      i = i + 1
      edgeIndex = curr%edgeIndex
      IF(inLattice%edges(edgeIndex)%primoVertex == vertexIndex) THEN
        neighborIndex = inLattice%edges(edgeIndex)%secundoVertex
      ELSE IF (inLattice%edges(edgeIndex)%secundoVertex == vertexIndex) THEN
        neighborIndex = inLattice%edges(edgeIndex)%primoVertex
      END IF
      neighborStates(i) = inLattice%vertices(neighborIndex)%qState
      neighborDistances(i) = inLattice%edges(edgeIndex)%length
    END DO

    enrg1
  END SUBROUTINE flip_qstate_isometric

  SUBROUTINE shift_position_isometric(inLattice,edgeHash,vertexIndex)
    USE basictypes, only: Lattice
    IMPLICIT NONE
    TYPE(Lattice),INTENT(INOUT)                       :: inLattice
    TYPE(LatticeEdgeSHead),DIMENSION(:),INTENT(IN)    :: edgeHash
    INTEGER,INTENT(IN)                                :: vertexIndex
  END SUBROUTINE shift_position_isometric

END MODULE evolvers
