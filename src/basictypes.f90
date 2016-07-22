! Tomas Mondragon
! Computer Scientist
! USACE ERDC Information Technology Laboratory
! Computational Analysis Branch
! Tomas.A.Mondragon@erdc.dren.mil
!   delivered 06 Feb 2016

! Subroutines
!   deallocate_lattice_node(which_node)
!   deallocate_lattice(which_lattice)

MODULE basictypes
  USE, INTRINSIC :: ISO_C_BINDING, only: C_INT,C_BOOL
  USE            :: triangle_c_wrap, only: C_REAL
  IMPLICIT NONE

  ! data types
  ! lattice node - measured and quantified points from the input plus some derived info
  TYPE lattice_node
    INTEGER(C_INT)                              :: n_index          ! self-referential index for convenience
    !INTEGER(C_INT),DIMENSION(:),ALLOCATABLE     :: neigh_edges      ! array of indices to connected edges(Not useful, may not be allocated and populalated to save time)
    !INTEGER(C_INT),DIMENSION(:),ALLOCATABLE     :: neigh_cells      ! array of indices to connected cells(Not useful, may not be allocated and populalated to save time)
    LOGICAL(C_BOOL)                             :: is_universe      ! flag for special universal nodes
    LOGICAL(C_BOOL)                             :: is_boundary      ! boundary marker flag
    REAL(C_REAL),DIMENSION(1:2)                 :: location         ! x-y coordinates of node
    INTEGER(C_INT)                              :: q_state          ! discretized state (spin state, orthogonal orientation state)
    REAL(C_REAL)                                :: major_axis       ! non-discretized state (orientation)
    INTEGER(C_INT)                              :: enumer           ! node membership in an enumerated grain
    REAL(C_REAL)                                :: enrg_node        ! partial energy of node
  END TYPE lattice_node

  ! Lattice edge - line segments connecting  the nodes and bounding the cells
  TYPE lattice_edge
    INTEGER(C_INT)                              :: e_index          ! self-referential index for convenience
    INTEGER(C_INT),DIMENSION(1:2)               :: neigh_cells      ! array of indices to connected cells
    INTEGER(C_INT),DIMENSION(1:2)               :: neigh_nodes      ! array of indices to connected nodes
    LOGICAL(C_BOOL)                             :: is_broken        ! flag indicating that boundary intersects edge
    LOGICAL(C_BOOL)                             :: was_traversed    ! flag marking whether boundary finder has already visited or not
    LOGICAL(C_BOOL)                             :: is_universe      ! flag for special universal nodes
    LOGICAL(C_BOOL)                             :: is_boundary      ! boundary marker flag
    REAL(C_REAL)                                :: enrg_edge        ! partial energy of edge
    REAL(C_REAL),DIMENSION(1:2)                 :: center           ! center of edge
    INTEGER(C_INT)                              :: enumer           ! edge membership in a grain
    REAL(C_REAL)                                :: orientation      ! orientation angle (in radians, (0,pi))
    REAL(C_REAL)                                :: length           ! length of edge
  END TYPE lattice_edge

  ! lattice cell - triangular elements bounded by three neighboring nodes and three neighboring edges
  TYPE lattice_cell
    INTEGER(C_INT)                              :: c_index          ! self-referential index for convenience
    INTEGER(C_INT),DIMENSION(1:3)               :: neigh_nodes      ! array of indices to connected nodes
    INTEGER(C_INT),DIMENSION(1:3)               :: neigh_edges      ! array of indices to connected edges
    INTEGER(C_INT),DIMENSION(1:3)               :: broken_edges     ! flags to indicate which edges are broken
    INTEGER(C_INT)                              :: num_broken_edges ! also useful for classifying cell quickly
    INTEGER(C_INT)                              :: num_traversals   ! number of times cell was traversed by boundary finder
    INTEGER(C_INT)                              :: q_state          ! discretized state (spin state, orthogonal orientation state)
    REAL(C_REAL)                                :: major_axis       ! non-discretized state (orientation)
    REAL(C_REAL),DIMENSION(1:2)                 :: center           ! Fermat point of triangular cell
    INTEGER(C_INT)                              :: enumer           ! cell membership in an enumerated grain
    LOGICAL(C_BOOL)                             :: is_universe      ! flag for special universal cells
  END TYPE lattice_cell

  ! lattice - container for raw and simplistic data. Nodes, cells, edges, & attached info
  TYPE lattice
    TYPE(lattice_node),DIMENSION(:),ALLOCATABLE :: nodes            ! array of nodes
    TYPE(lattice_edge),DIMENSION(:),ALLOCATABLE :: edges            ! array of edges
    TYPE(lattice_cell),DIMENSION(:),ALLOCATABLE :: cells            ! array of cells
    INTEGER(C_INT)                              :: num_nodes        ! number of nodes for quick reference
    INTEGER(C_INT)                              :: num_edges        ! number of edges for quick reference
    INTEGER(C_INT)                              :: num_cells        ! number of cells for quick reference
    INTEGER(C_INT)                              :: periodicity      ! 0 =: aperiodic, 1 =: periodic in y, 2 =: periodic in x, 3 =: periodic everywhere
    INTEGER(C_INT)                              :: num_states       ! number of possible discrete states
    INTEGER(C_INT)                              :: mc_states        ! number of Monte-Carlo "sweeps" to perform during generation
    REAL(C_REAL)                                :: temperature      ! temperature of crystal during evolution (useful during Potts generation)
    REAL(C_REAL)                                :: en_per_length    ! constant used in grain boundary energy calculations (electronVolts per ?m)
    REAL(C_REAL),DIMENSION(1:2)                 :: grid_spacing     ! x and y spacing used in rectangular grid in Potts generator
    REAL(C_REAL),DIMENSION(1:2)                 :: x_bounds         ! max and min x coordinate
    REAL(C_REAL),DIMENSION(1:2)                 :: y_bounds         ! max and min y coordinate
    INTEGER(C_INT),DIMENSION(:),ALLOCATABLE     :: grain_size_array ! histogram of number of nodes in each grain, numbered the same as enumer
    INTEGER(C_INT)                              :: max_grain_index  ! Highest integer in enumer. Put here to avoid recalculation.
    INTEGER(C_INT)                              :: num_grains       ! Number of grains in crystal
    INTEGER(C_INT)                              :: max_grain_size   ! number of nodes inside of biggest grain
    REAL(C_REAL)                                :: avg_size         ! Average grain size, if it is useful
    INTEGER(C_INT), DIMENSION(0:7)              :: grain_size_histo ! Histogram of grain sizes, binned by the proportion of map coverage
  END TYPE lattice

CONTAINS

  ! Properly deallocate a lattice_node
  ! USAGE: CALL deallocate_lattice_node(which_node)
  SUBROUTINE deallocate_lattice_node(which_node)
    IMPLICIT NONE
    TYPE(lattice_node),INTENT(INOUT)            :: which_node

    !IF( ALLOCATED(which_node%neigh_edges) ) DEALLOCATE(which_node%neigh_edges)
    !IF( ALLOCATED(which_node%neigh_cells) ) DEALLOCATE(which_node%neigh_cells)
  END SUBROUTINE deallocate_lattice_node

  ! Properly deallocate a lattice
  ! USAGE: CALL deallocate_lattice(which_lattice)
  SUBROUTINE deallocate_lattice(which_lattice)
    USE FileHandling, only: stderr
    IMPLICIT NONE
    TYPE(lattice),INTENT(INOUT)                 :: which_lattice
    INTEGER                                     :: i

    IF( ALLOCATED(which_lattice%nodes)) THEN
      IF( which_lattice%num_nodes .ne. SIZE(which_lattice%nodes) ) THEN
        WRITE(stderr,'(A,I8)') "WARNING: During node deallocation, lattice reported ", which_lattice%num_nodes, " nodes"
        WRITE(stderr,'(A,I8,A)') "   but there were ",SIZE(which_lattice%nodes)," nodes to deallocate."
      END IF
      DO i=LBOUND(which_lattice%nodes),UBOUND(which_lattice%nodes)
        CALL deallocate_lattice_node(which_lattice%nodes(i))
      END DO
      DEALLOCATE(which_lattice%nodes)
    ELSE IF( which_lattice%num_nodes .ne. 0 ) THEN
      WRITE(stderr,'(A,I8)') "WARNING: During node deallocation, lattice reported ", which_lattice%num_nodes, " nodes"
      WRITE(stderr,'(A)') "   but there were 0 nodes to deallocate."
    END IF

    IF( ALLOCATED(which_lattice%edges)) THEN
      IF( which_lattice%num_edges .ne. SIZE(which_lattice%edges) ) THEN
        WRITE(stderr,'(A,I8)') "WARNING: During edge deallocation, lattice reported ", which_lattice%num_edges, " edges"
        WRITE(stderr,'(A,I8,A)') "   but there were ",SIZE(which_lattice%edges)," edges to deallocate."
      END IF
      DEALLOCATE(which_lattice%edges)
    ELSE IF( which_lattice%num_edges .ne. 0 ) THEN
      WRITE(stderr,'(A,I8)') "WARNING: During edge deallocation, lattice reported ", which_lattice%num_edges, " edges"
      WRITE(stderr,'(A)') "   but there were 0 edges to deallocate."
    END IF

    IF( ALLOCATED(which_lattice%cells)) THEN
      IF( which_lattice%num_cells .ne. SIZE(which_lattice%cells) ) THEN
        WRITE(stderr,'(A,I8)') "WARNING: During cell deallocation, lattice reported ", which_lattice%num_cells, " cells"
        WRITE(stderr,'(A,I8,A)') "   but there were ",SIZE(which_lattice%cells)," cells to deallocate."
      END IF
      DEALLOCATE(which_lattice%cells)
    ELSE IF( which_lattice%num_cells .ne. 0 ) THEN
      WRITE(stderr,'(A,I8)') "WARNING: During cell deallocation, lattice reported ", which_lattice%num_cells, " cells"
      WRITE(stderr,'(A)') "   but there were 0 cells to deallocate."
    END IF
    IF( ALLOCATED(which_lattice%grain_size_array) ) DEALLOCATE(which_lattice%grain_size_array)
  END SUBROUTINE deallocate_lattice

END MODULE basictypes
