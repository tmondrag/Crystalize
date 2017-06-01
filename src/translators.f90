MODULE translators
  IMPLICIT NONE

  TYPE Lattice
    TYPE(LatticeVertex),ALLOCATABLE     :: vertices(:)
    TYPE(LatticeEdge),ALLOCATABLE       :: edges(:)
    TYPE(LatticeHalfEdge),ALLOCATABLE   :: halfEdges(:)
    TYPE(LatticeFacet),ALLOCATABLE      :: facets(:)
  END TYPE Lattice

  TYPE LatticeEdge
    INTEGER                             :: sinesterIndex
    INTEGER                             :: dexterIndex
    REAL                                :: length
    REAL                                :: energy
  END TYPE LatticeEdge

  TYPE LatticeHalfEdge
    INTEGER                             :: vertexIndex
    INTEGER                             :: facetIndex
    INTEGER                             :: twinIndex
    INTEGER                             :: nextIndex
    INTEGER                             :: prevIndex
    INTEGER                             :: edgeInfoIndex
    LOGICAL                             :: isAtBoundary
  END TYPE LatticeHalfEdge

  TYPE LatticeVertex
    REAL,DIMENSION(1:2)                 :: location
    INTEGER                             :: halfEdgeIndex
    LOGICAL                             :: isAtBoundary
  END TYPE LatticeVertex

  TYPE LatticeFacet
    INTEGER                             :: edgeIndex
    LOGICAL                             :: isUniverse
  END TYPE LatticeFacet

CONTAINS
  SUBROUTINE transfer_triangle_to_lattice(inFTriangles,outLattice)
    USE triangle_c_wrap, only: triangulateio,f_triangulateio
    IMPLICIT NONE
    TYPE(f_triangulateio),INTENT(IN)    :: inFTriangles
    TYPE(Lattice),INTENT(INOUT)         :: outLattice
    INTEGER                             :: numVertices,numFacets,numEdges

    numVertices = inFTriangles%numberofpoints
    numFacets = inFTriangles%numberoftriangles + 1
    numEdges = numVertices+numFacets-2 ! Euler Characteristic

    ALLOCATE(outLattice%vertices(1:numVertices))
    ALLOCATE(outLattice%edges(1:numEdges))
    ALLOCATE(outLattice%halfEdges(1:2*numEdges))
    ALLOCATE(outLattice%facets(1:numFacets))
  END SUBROUTINE transfer_triangle_to_lattice


END MODULE translators
