! Tomas Mondragon
! Computer Scientist
! USACE ERDC Information Technology Laboratory
! Computational Analysis Branch
! Tomas.A.Mondragon@erdc.dren.mil
!   delivered 06 May 2017

! Subroutines
!   init_lattice(which_lattice)
!   deallocate_lattice_node(which_node)
!   deallocate_lattice(which_lattice)

MODULE basictypes
  USE, INTRINSIC :: ISO_C_BINDING, only: C_INT,C_BOOL
  USE            :: triangle_c_wrap, only: C_REAL
  IMPLICIT NONE

  ! data types
  ! Lattice types - DCEL and winded edge data stuctures won't do well, since Triangle's data stucture is point centered rater than edge centered
  ! lattice vertex - measured and quantified points from the input plus some derived info
  TYPE LatticeVertex
    INTEGER(C_INT)                              :: nIndex           ! self-referential index for convenience
    INTEGER(C_INT)                              :: boundary         ! boundary marker flag
    REAL(C_REAL),DIMENSION(1:2)                 :: location         ! x-y coordinates of vertex
    INTEGER(C_INT)                              :: qState           ! discretized state (spin state, orthogonal orientation state)
    INTEGER(C_INT)                              :: enumer           ! vertex membership in an enumerated grain
    REAL(C_REAL)                                :: enrgVertex       ! partial energy of vertex (Half the sum of partial energies of connected edges)
  END TYPE LatticeVertex

  ! Lattice edge - line segments connecting  the vertices and bounding the facets
  TYPE LatticeEdge
    INTEGER(C_INT)                              :: eIndex           ! self-referential index for convenience
    LOGICAL(C_BOOL)                             :: isBroken         ! flag indicating that boundary intersects edge
    LOGICAL(C_BOOL)                             :: wasTraversed     ! flag marking whether boundary finder has already visited or not
    LOGICAL(C_BOOL)                             :: is_fake          ! flag marking edge used to implement holes in fancy polygons
    INTEGER(C_INT)                              :: boundary         ! boundary marker flag
    REAL(C_REAL)                                :: enrgEdge         ! partial energy of edge
    REAL(C_REAL),DIMENSION(1:2)                 :: center           ! center of edge
    REAL(C_REAL)                                :: orientation      ! orientation angle (in radians, (0,pi))
    REAL(C_REAL)                                :: length           ! length of edge
    !turn into ints. cant point into an array
    INTEGER(C_INT)                              :: primoVertex      ! index to one of this edge's endpoints
    INTEGER(C_INT)                              :: secundoVertex    ! index to the other edge endpoints
    INTEGER(C_INT)                              :: dexterFacet      ! index to the facet to the right of the edge
    INTEGER(C_INT)                              :: sinesterFacet    ! index to the facet to the left of the edge
  END TYPE LatticeEdge

  ! Lattice facet - triangular elements bounded by three neighboring nodes and three neighboring edges (or it can have more if you are doing something fancy)
  TYPE LatticeFacet
    INTEGER(C_INT)                              :: fIndex           ! self-referential index for convenience
    INTEGER(C_INT)                              :: numBrokenEdges   ! Count of edges crossed by boundary, also useful for classifying the facet/boundary type
    INTEGER(C_INT)                              :: numTraversals    ! number of times facet was traversed by boundary finder
    REAL(C_REAL),DIMENSION(1:2)                 :: center           ! Fermat point of triangular facet
    REAL(C_REAL)                                :: area             ! area of facet
    REAL(C_REAL)                                :: enrgFacet        ! partial energy for facet (Half the sum of partial energies of connected edges)
    LOGICAL(C_BOOL)                             :: isUniverse       ! flag for special universal facet
    INTEGER(C_INT)                              :: primoVertex      ! index to one of this facet's vertices
    INTEGER(C_INT)                              :: secundoVertex    ! index to the next vertex in clockwise order
    INTEGER(C_INT)                              :: tertioVertex     ! index to the third vertex in clockwise order
    INTEGER(C_INT)                              :: primoEdge        ! index to one of this facet's edges
    INTEGER(C_INT)                              :: secundoEdge      ! index to the next edge in clockwise order
    INTEGER(C_INT)                              :: tertioEdge       ! index to the third edge in clockwise order
  END TYPE latticeFacet

  ! Lattice facet search data structure list item
  TYPE LatticeFacetSItem
    LOGICAL(C_BOOL)                             :: isHead
    INTEGER(C_INT)                              :: vertexIndex
    INTEGER(C_INT)                              :: facetIndex
    TYPE(LatticeFacetSItem),POINTER             :: next
  CONTAINS
    PROCEDURE :: init   => init_lattice_facet_list
    PROCEDURE :: push   => push_lattice_facet_list
    PROCEDURE :: pop    => pop_lattice_facet_list
    PROCEDURE :: delete => delete_lattice_facet_list
    PROCEDURE :: search => search_lattice_facet_list
  END TYPE LatticeFacetSItem

  ! Lattice facet search data structure list head/array item
  TYPE LatticeFacetSHead
    TYPE(LatticeFacetSItem)                     :: listHead
  END TYPE LatticeFacetSHead

  ! Lattice - container for raw and simplistic data. vertices, facets, edges, & attached info
  TYPE Lattice
    TYPE(LatticeVertex),DIMENSION(:),ALLOCATABLE:: vertices         ! array of vertices
    TYPE(LatticeEdge),DIMENSION(:),ALLOCATABLE  :: edges            ! array of edges
    TYPE(LatticeFacet),DIMENSION(:),ALLOCATABLE :: facets           ! array of facets
    TYPE(LatticeFacetSHead),DIMENSION(:),ALLOCATABLE:: facetReference
    INTEGER(C_INT)                              :: numVertices      ! number of vertices for quick reference
    INTEGER(C_INT)                              :: numEdges         ! number of edges for quick reference
    INTEGER(C_INT)                              :: numFacets        ! number of facets for quick reference
    INTEGER(C_INT)                              :: periodicity      ! 0 =: aperiodic, 1 =: periodic in y, 2 =: periodic in x, 3 =: periodic everywhere
    INTEGER(C_INT)                              :: numStates        ! number of possible discrete states
    INTEGER(C_INT)                              :: mcSweeps         ! number of Monte-Carlo "sweeps" to perform during generation
    REAL(C_REAL)                                :: temperature      ! temperature of crystal during evolution (useful during Potts generation)
    REAL(C_REAL)                                :: enrgPerLength    ! constant used in grain boundary energy calculations (electronVolts per ?m)
    REAL(C_REAL),DIMENSION(1:2)                 :: gridSpacing      ! x and y spacing used in rectangular grid in Potts generator
    REAL(C_REAL),DIMENSION(1:2)                 :: xBounds          ! min and max x coordinate
    REAL(C_REAL),DIMENSION(1:2)                 :: yBounds          ! min and max y coordinate
    INTEGER(C_INT),DIMENSION(:),ALLOCATABLE     :: grainSizeArray   ! histogram of number of nodes in each grain, numbered the same as enumer
    INTEGER(C_INT)                              :: maxGrainIndex    ! Highest integer in enumer. Put here to avoid recalculation.
    INTEGER(C_INT)                              :: numGrains        ! Number of grains in crystal
    INTEGER(C_INT)                              :: maxGrainSize     ! number of nodes inside of biggest grain
    REAL(C_REAL)                                :: avgSize          ! Average grain size, if it is useful
    INTEGER(C_INT), DIMENSION(0:7)              :: grainSizeHisto   ! Histogram of grain sizes, binned by the proportion of map coverage
  END TYPE Lattice

  ! Grain facet - one whole grain undivided by triangulation
  TYPE GrainFacet
    INTEGER(C_INT)                              :: gFIndex          ! self-referential index for convenience
    TYPE(LatticeVertex), POINTER                :: samplePoint      ! Sample point within the grain
    TYPE(GrainHalfEdge), POINTER                :: boundingEdge     ! one of the half-edges that bound the grain
    INTEGER(C_INT)                              :: latticeVertexCount
    INTEGER(C_INT)                              :: qState           ! discretized state of grain (spin state, orthogonal orientation state)
    REAL(C_REAL)                                :: majorAxis        ! non-discretized state of grain (major axis orientation angle)
  END TYPE GrainFacet

  ! Grain half-edge - entity representing half of an actual grain boundary segment
  TYPE GrainHalfEdge
    INTEGER(C_INT)                              :: gEIndex          ! self-referential index for convenience
    TYPE(GrainHalfEdge), POINTER                :: twin             ! complimentary half-edge going in the opposite direction
    TYPE(GrainHalfEdge), POINTER                :: next             ! the next half-edge on the border of the neighboring facet
    TYPE(GrainHalfEdge), POINTER                :: prev             ! the previous half-edge on the border of the neighboring facet
    TYPE(GrainBoundary), POINTER                :: parent           ! the grain boundary that contains this border segment
    TYPE(GrainFacet), POINTER                   :: neighborFacet    ! the facet bounded by this half-edge
    TYPE(GrainVertex), POINTER                  :: originVertex     ! the origin of this half-edge
  END TYPE GrainHalfEdge

  ! Grain vertex - the boundary points of a grain
  TYPE GrainVertex
    INTEGER(C_INT)                              :: gVIndex          ! self-referential index for convenience
    TYPE(GrainHalfEdge), POINTER                :: neighborHalfEdge ! one (of possibly many) of the half edges originating at this vertex
    REAL(C_REAL),DIMENSION(1:2)                 :: location         ! x-y coordinates of vertex
    LOGICAL(C_BOOL)                             :: isJunction       ! True if the vertex is the junction of boundaries
  END TYPE GrainVertex

  TYPE GrainBoundary
    INTEGER(C_INT)                              :: gBIndex          ! self-referential index for convenience
    TYPE(GrainHalfEdge), POINTER                :: startSegment     ! a half-edge on the grain boundary
    REAL(C_REAL),DIMENSION(:,:),ALLOCATABLE     :: fittedPoints     ! original points on boundary
    REAL(C_REAL),DIMENSION(:,:),ALLOCATABLE     :: smoothedPoints   ! smoothed out points on the boundary
    INTEGER(C_INT)                              :: numPoints        ! number of points in boundary
    REAL(C_REAL)                                :: a,b,c,d          ! boundary curve parameters
    REAL(C_REAL)                                :: alpha_m,beta_m,gamma_m ! curve fitting parameters
    REAL(C_REAL)                                :: p,q,r,s          ! point fitting parameters
  END TYPE GrainBoundary

CONTAINS

  ! Initialize a lattice facet search sublist
  SUBROUTINE init_lattice_facet_list(this,firstIndex)
    IMPLICIT NONE
    CLASS(LatticeFacetSItem),INTENT(INOUT),TARGET :: this
    INTEGER(C_INT),INTENT(IN)                     :: firstIndex

    this%isHead = .TRUE.
    this%vertexIndex = firstIndex
    this%next => this
  END SUBROUTINE init_lattice_facet_list

  ! push a lattice facet onto the search sublist
  SUBROUTINE push_lattice_facet_list(this,secondIndex,facetIndex)
    IMPLICIT NONE
    CLASS(LatticeFacetSItem),INTENT(INOUT)        :: this
    INTEGER(C_INT),INTENT(IN)                     :: secondIndex,facetIndex
    TYPE(LatticeFacetSItem),POINTER               :: curr

    ALLOCATE(curr)
    curr%isHead = .FALSE.
    curr%vertexIndex = secondIndex
    curr%facetIndex = facetIndex
    curr%next => this%next
    this%next => curr
  END SUBROUTINE push_lattice_facet_list

  ! mosly useless, but here for the sake of completeness
  FUNCTION pop_lattice_facet_list(this) RESULT(facetIndex)
    IMPLICIT NONE
    CLASS(LatticeFacetSItem)                      :: this
    INTEGER                                       :: facetIndex
    TYPE(LatticeFacetSItem),POINTER               :: prev,curr,next

    curr => this%next
    next => curr%next
    IF (this%isHead) THEN
      facetIndex = curr%facetIndex
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
          facetIndex = curr%facetIndex
          prev%next => next
          curr%next => curr
          DEALLOCATE(curr)
          RETURN
        END IF
      END DO
    END IF
  END FUNCTION pop_lattice_facet_list

  ! delete and deallocate a facet list the proper way
  SUBROUTINE delete_lattice_facet_list(this)
    IMPLICIT NONE
    CLASS(LatticeFacetSItem),INTENT(IN),TARGET    :: this
    TYPE(LatticeFacetSItem),POINTER               :: prev,curr,next

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
  END SUBROUTINE delete_lattice_facet_list

  ! searc hth list for a facet to the left of a certain edge
  FUNCTION search_lattice_facet_list(this,secondIndex) RESULT(facetIndex)
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: stderr => ERROR_UNIT
    IMPLICIT NONE
    CLASS(LatticeFacetSItem)                      :: this
    INTEGER                                       :: facetIndex, firstIndex, secondIndex
    TYPE(LatticeFacetSItem),POINTER               :: prev,curr

    curr => this%next
    IF (.NOT. this%isHead) THEN
      DO
        prev => curr
        curr => prev%next
        IF (prev%isHead) THEN
          firstIndex = prev%vertexIndex
          EXIT
        END IF
      END DO
    ELSE
      firstIndex = this%vertexIndex
    END IF
    DO
      IF (curr%isHead) THEN
        WRITE(stderr,"(A,I8,A,I8,A)") "WARNING: Search for a facet with vertices ", firstIndex, " , ", secondIndex, " in", &
             " in counter-clockwise order failed. Assuming line segment is part of outer bounds."
        facetIndex = 0
        RETURN
      END IF
      IF (curr%vertexIndex == secondIndex) THEN
        facetIndex = curr%facetIndex
        RETURN
      END IF
      curr => curr%next
    END DO
  END FUNCTION search_lattice_facet_list

END MODULE basictypes
