! Tomas Mondragon
! Computer Scientist
! USACE ERDC Information Technology Laboratory
! Computational Analysis Branch
! Tomas.A.Mondragon@erdc.dren.mil
!   delivered 16 Apr 2017

MODULE evolvers
  USE, INTRINSIC :: ISO_C_BINDING, only: C_INT,C_BOOL
  USE            :: triangle_c_wrap, only: C_REAL
  ! Lattice facet search data structure list item
  TYPE LatticeEdgeSItem
    LOGICAL(C_BOOL)                             :: isHead
    INTEGER(C_INT)                              :: vertexIndex
    INTEGER(C_INT)                              :: edgeIndex
    TYPE(LatticeEdgeSItem),POINTER              :: next
  CONTAINS
    PROCEDURE :: init   => init_lattice_edge_list
    PROCEDURE :: push   => push_lattice_edge_list
    PROCEDURE :: pop    => pop_lattice_edge_list
    PROCEDURE :: delete => delete_lattice_edge_list
  END TYPE LatticeEdgeSItem

  ! Lattice facet search data structure list head/array item
  TYPE LatticeEdgeSHead
    TYPE(LatticeEdgeSItem)                     :: listHead
  END TYPE LatticeEdgeSHead

CONTAINS
  SUBROUTINE evolve_qstates
  END SUBROUTINE evolve_qstates

  SUBROUTINE evolve_positions
  END SUBROUTINE evolve_positions

  SUBROUTINE flip_qstate
  END SUBROUTINE flip_qstate

  SUBROUTINE shift_position
  END SUBROUTINE shift_position

  FUNCTION build_edge_hash(inLattice) RESULT(edgeHash)
    IMPLICIT NONE
    USE basictypes, only: Lattice
    TYPE(Lattice)                                     :: inLattice
    TYPE(LatticeEdgeSHead),DIMENSION(:),ALLOCATABLE   :: edgeHash

    ALLOCATE(edgeHash(1:inLattice%numVertices))
  END FUNCTION build_edge_hash

  ! Initialize a lattice edge search sublist
  SUBROUTINE init_lattice_edge_list(this,vertexIndex)
    IMPLICIT NONE
    CLASS(LatticeEdgeSItem),INTENT(INOUT),TARGET  :: this
    INTEGER(C_INT),INTENT(IN)                     :: vertexIndex

    this%isHead = .TRUE.
    this%vertexIndex = vertexIndex
    this%next => this
  END SUBROUTINE init_lattice_edge_list

  ! push a lattice edge onto the search sublist
  SUBROUTINE push_lattice_edge_list(this,vertexIndex,edgeIndex)
    IMPLICIT NONE
    CLASS(LatticeEdgeSItem),INTENT(INOUT)         :: this
    INTEGER(C_INT),INTENT(IN)                     :: vertexIndex,edgeIndex
    TYPE(LatticeEdgeSItem),POINTER                :: curr

    ALLOCATE(curr)
    curr%isHead = .FALSE.
    curr%vertexIndex = vertexIndex
    curr%facetIndex = facetIndex
    curr%next => this%next
    this%next => curr
  END SUBROUTINE push_lattice_facet_list

  ! mosly useless, but here for the sake of completeness
  FUNCTION pop_lattice_edge_list(this) RESULT(edgeIndex)
    IMPLICIT NONE
    CLASS(LatticeEdgeSItem)                       :: this
    INTEGER                                       :: edgeIndex
    TYPE(LatticeEdgeSItem),POINTER                :: prev,curr,next

    curr => this%next
    next => curr%next
    IF (this%isHead) THEN
      edgeIndex = curr%edgeIndex
      this%next => next
      curr%next => curr
      DEALLOCATE(curr)
      RETURN
    ELSE
      DO
        prev => curr
        curr => prev%next
        next => prev%next%next
        IF (prev%isHead) THEN
          edgeIndex = curr%edgeIndex
          prev%next => next
          curr%next => curr
          DEALLOCATE(curr)
          RETURN
        END IF
      END DO
    END IF
  END FUNCTION pop_lattice_edge_list

  ! delete and deallocate a edge list the proper way
  SUBROUTINE delete_lattice_edge_list(this)
    IMPLICIT NONE
    CLASS(LatticeEdgeSItem),INTENT(IN),TARGET    :: this
    TYPE(LatticeEdgeSItem),POINTER               :: prev,curr,next

    prev => this
    curr => this%next
    next => this%next%next
    IF (.NOT. this%isHead) THEN
      DO
        prev => curr
        curr => prev%next
        next => prev%next%next
        IF (prev%isHead) EXIT
      END DO
    END IF
    DO
      IF (curr%isHead) EXIT
      curr%next => curr
      DEALLOCATE(curr)
      prev%next => next
      curr => prev%next
      next => prev%next%next
    END DO
  END SUBROUTINE delete_lattice_edge_list
END MODULE evolvers
