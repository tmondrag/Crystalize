! Tomas Mondragon
! Computer Scientist
! USACE ERDC Information Technology Laboratory
! Computational Analysis Branch
! Tomas.A.Mondragon@erdc.dren.mil
!   delivered 16 Apr 2017

MODULE evolvers
  USE, INTRINSIC :: ISO_C_BINDING, only: C_INT,C_BOOL
  USE            :: triangle_c_wrap, only: C_REAL

  REAL(C_REAL),PARAMETER    :: k_B = 8.6173303e-5_C_REAL   ! Boltzmann constant in eV/K


CONTAINS
  SUBROUTINE evolve_qstates(inLattice)
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: stderr -> ERROR_UNIT, stdout -> OUTPUT_UNIT
    USE basictypes, only: Lattice
    IMPLICIT NONE
    TYPE(Lattice),INTENT(INOUT)                       :: inLattice
    INTEGER                                           :: i,j,k
    REAL                                              :: tmp

    DO i = 1, inLattice%mcSweeps
      DO j = 1, inLattice%numVertices
        CALL RANDOM_NUMBER(tmp)
        k = CEILING(inLattice%numVertices*tmp)
        CALL flip_qstate_isometric(inLattice,k)
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

  SUBROUTINE flip_qstate_isometric(inLattice,vertexIndex)
    USE basictypes, only: Lattice
    USE energyCalc, only: calculate_edge_energy
    IMPLICIT NONE
    TYPE(Lattice),INTENT(INOUT)                       :: inLattice
    TYPE(LatticeEdgeSHead),DIMENSION(:)               :: edgeHash
    INTEGER,INTENT(IN)                                :: vertexIndex
    REAL                                              :: enrg1,enrg2,tmp, flipprob
    INTEGER                                           :: presentState, randomState
    INTEGER                                           :: numNeighbors,i,edgeIndex,neighborIndex
    INTEGER,DIMENSION(:)                              :: neighborStates
    REAL,DIMENSION(:)                                 :: neighborDistances
    TYPE(LatticeEdgeSItem),POINTER                    :: curr

    edgeHash = inLattice%edgeHash
    presentState = inLattice%vertices(vertexIndex)%qState
    CALL RANDOM_NUMBER(tmp)
    tmp = REAL(latticeData%numStates,C_REAL)*tmp
    randomState = INT(tmp,C_INT) + 1
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
    enrg1 = 0.0
    enrg2 = 0.0
    DO
      curr => curr%next
      IF (curr%isHead) EXIT
      i = i + 1
      edgeIndex = curr%edgeIndex
      neighborIndex = curr%vertexIndex
      neighborStates(i) = inLattice%vertices(neighborIndex)%qState
      neighborDistances(i) = inLattice%edges(edgeIndex)%length
      enrg1 = enrg1 + 0.5*calculate_edge_energy(inLattice%enrgScale,inLattice%restBose, &
                                                inLattice%restFermi,presentState,neighborStates(i), &
                                                neighborDistances(i))
      enrg2 = enrg2 + 0.5*calculate_edge_energy(inLattice%enrgScale,inLattice%restBose, &
                                                inLattice%restFermi,randomState,neighborStates(i), &
                                                neighborDistances(i))
    END DO
    flipprob = 0.5*ERF((enrg1-enrg2)/(k_B*inLattice%temperature)) + 0.5
    CALL RANDOM_NUMBER(tmp)
    IF ( tmp > fliprob ) THEN
      ! no qstate flip
      inLattice%vertices(vertexIndex)%enrgVertex = enrg1
    ELSE !( tmp <= flipprob )
      ! qstate flip
      inLattice%vertices(vertexIndex)%qState = randomState
      inLattice%vertices(vertexIndex)%enrgVertex = enrg2
      ! update edge energies
      curr => edgeHash(vertexIndex)%listHead
      i = 0
      DO
        curr => curr%next
        IF (curr%isHead) EXIT
        i = i + 1
        edgeIndex = curr%edgeIndex
        inLattice%edges(edgeIndex)%enrgEdge = calculate_edge_energy(inLattice%enrgScale,inLattice%restBose, &
                                                  inLattice%restFermi,randomState,neighborStates(i), &
                                                  neighborDistances(i))
      END DO
    END IF
  END SUBROUTINE flip_qstate_isometric

  SUBROUTINE shift_position_isometric(inLattice,edgeHash,vertexIndex)
    USE basictypes, only: Lattice
    IMPLICIT NONE
    TYPE(Lattice),INTENT(INOUT)                       :: inLattice
    TYPE(LatticeEdgeSHead),DIMENSION(:),INTENT(IN)    :: edgeHash
    INTEGER,INTENT(IN)                                :: vertexIndex
  END SUBROUTINE shift_position_isometric

END MODULE evolvers
