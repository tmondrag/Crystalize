MODULE trans_lattice_triangle
  IMPLICIT NONE

CONTAINS
  ! This subroutine has no use right now
  SUBROUTINE transfer_lattice_to_triangle(inLattice,f_shape,c_shape)
    USE basictypes, only: Lattice
    USE triangle_c_wrap, only: triangulateio,f_triangulateio,allocate_points_ftoc,allocate_segments,allocate_regions,allocate_holes
    USE triangle_c_wrap, only: C_REAL,allocate_triangles
    IMPLICIT NONE
    TYPE(lattice),INTENT(IN)                          :: inLattice
    TYPE(triangulateio),INTENT(OUT)                   :: c_shape
    TYPE(f_triangulateio),INTENT(OUT)                 :: f_shape
    INTEGER                                           :: i

    f_shape%numberofpoints = inLattice%numVertices
    f_shape%numberofpointattributes = 2               ! qState, enumer
    f_shape%numberofsegments = inLattice%numEdges
    f_shape%numberoftriangles = inLattice%numFacets
    f_shape%numberofcorners = 3
    f_shape%numberoftriangleattributes = 0
    f_shape%numberofholes = 0
    f_shape%numberofregions = 0

    CALL allocate_points_ftoc(f_shape,c_shape)
    CALL allocate_segments(f_shape,c_shape)
    CALL allocate_triangles(f_shape,c_shape)
    CALL allocate_regions(f_shape,c_shape)
    CALL allocate_holes(f_shape,c_shape)

    DO i=1,inLattice%numVertices
      f_shape%pointlist(2*i-1) = inLattice%vertices(i)%location(1)
      f_shape%pointlist(2*i) = inLattice%vertices(i)%location(2)
      f_shape%pointAttributelist(2*i-1) = REAL(inLattice%vertices(i)%qState,C_REAL)
      f_shape%pointAttributelist(2*i) = REAL(inLattice%vertices(i)%enumer,C_REAL)
      f_shape%pointmarkerlist = inLattice%vertices(i)%boundary
    END DO

    DO i=1,inLattice%numEdges
      f_shape%segmentlist(2*i-1) = inLattice%edges(i)%primoVertex
      f_shape%segmentlist(2*i) = inLattice%edges(i)%secundoVertex
      f_shape%segmentmarkerlist = inLattice%edges(i)%boundary
    END DO

    ! this assumes triangular facets
    DO i=1,inLattice%numFacets
      f_shape%trianglelist(3*i-2) = inLattice%facets(i)%primoVertex
      f_shape%trianglelist(3*i-1) = inLattice%facets(i)%secundoVertex
      f_shape%trianglelist(3*i) = inLattice%facets(i)%tertioVertex
    END DO
    f_shape%trianglearealist = 0.0
  END SUBROUTINE transfer_lattice_to_triangle

  ! This is useful after the first triangulation in preparation for grain discovery and border finding
  SUBROUTINE transfer_triangle_to_lattice(f_shape,outLattice)
    USE, INTRINSIC :: ISO_FORTRAN_ENV, only: stderr => ERROR_UNIT
    USE, INTRINSIC :: ISO_C_BINDING, only: C_INT
    USE basictypes, only: Lattice, LatticeFacetSItem
    USE triangle_c_wrap, only: f_triangulateio,allocate_points_ftoc,allocate_segments,allocate_regions,allocate_holes
    USE triangle_c_wrap, only: C_REAL
    USE energyCalc, only: calculate_edge_energy(enrgScale,restBose,restFermi,q1,q2,length)
    IMPLICIT NONE
    TYPE(lattice),INTENT(OUT)                         :: outLattice
    TYPE(f_triangulateio),INTENT(IN)                  :: f_shape
    INTEGER                                           :: i, err, facetIndex, num_attr
    INTEGER(C_INT)                                    :: primoVertex, secundoVertex, tertioVertex
    REAL(C_REAL),DIMENSION(1:2)                       :: edgeA,edgeB,edgeC
    REAL(C_REAL)                                      :: lenA,lenB,lenC
    REAL(C_REAL)                                      :: baryA,baryB,baryC
    REAL(C_REAL)                                      :: areaABC
    REAL(C_REAL)                                      :: enrgScale,restBose,restFermi,length
    INTEGER(C_INT)                                    :: q1,q2

    outLattice%numVertices = f_shape%numberofpoints
    outLattice%numEdges = f_shape%numberofsegments
    outLattice%numFacets = f_shape%numberoftriangles

    ALLOCATE(outLattice%vertices(1:outLattice%numVertices), STAT=err)
    IF (err /= 0) THEN
      WRITE(stderr,*) "outLattice%vertices: Allocation request denied"
      STOP
    END IF
    ALLOCATE(outLattice%facetReference(1:outLattice%numVertices), STAT=err)
    IF (err /= 0) THEN
      WRITE(stderr,*) "outLattice%facetReference: Allocation request denied"
      STOP
    END IF
    ALLOCATE(outLattice%edgeHash(1:outLattice%numVertices), STAT=err)
    IF (err /= 0) THEN
      WRITE(stderr,*) "outLattice%edgeHash: Allocation request denied"
      STOP
    END IF
    ALLOCATE(outLattice%edges(1:outLattice%numEdges), STAT=err)
    IF (err /= 0) THEN
      WRITE(stderr,*) "outLattice%edges: Allocation request denied"
      STOP
    END IF
    ALLOCATE(outLattice%facets(0:outLattice%numFacets), STAT=err)
    IF (err /= 0) THEN
      WRITE(stderr,*) "outLattice%Facets: Allocation request denied"
      STOP
    END IF

    DO i = 1,outLattice%numVertices
      outLattice%vertices(i)%nIndex = i
      outLattice%vertices(i)%boundary = f_shape%pointmarkerlist(i)
      outLattice%vertices(i)%location = (/f_shape%pointlist(2*i-1), f_shape%pointlist(2*i)/)
      IF(f_shape%numberofpointattributes > 0) THEN
        num_attr = f_shape%numberofpointattributes
        outLattice%vertices(i)%qState = INT(f_shape%pointattributelist(num_attr*i-num_attr+1),C_INT)
      ELSE
        outLattice%vertices(i)%qState = 0
      END IF
      CALL outLattice%facetReference(i)%listHead%init(i)
      CALL outLattice%edgeHash(i)%listHead%init(i)
    END DO
    DO i = 1,outLattice%numFacets
      outLattice%facets(i)%fIndex = i
      outLattice%facets(i)%isUniverse = .FALSE.
      outLattice%facets(i)%primoVertex = f_shape%trianglelist(3*i-2)
      outLattice%facets(i)%secundoVertex = f_shape%trianglelist(3*i-1)
      outLattice%facets(i)%tertioVertex = f_shape%trianglelist(3*i)
      primoVertex = outLattice%facets(i)%primoVertex
      secundoVertex = outLattice%facets(i)%secundoVertex
      tertioVertex = outLattice%facets(i)%tertioVertex
      edgeA = outLattice%vertices(tertioVertex)%location - outLattice%vertices(secundoVertex)%location
      edgeB = outLattice%vertices(primoVertex)%location - outLattice%vertices(tertioVertex)%location
      edgeC = outLattice%vertices(secundoVertex)%location - outLattice%vertices(primoVertex)%location
      lenA = NORM2(edgeA)
      lenB = NORM2(edgeB)
      lenC = NORM2(edgeC)
      areaABC = (edgeC(1)*edgeB(2)-edgeB(1)*edgeC(2))/2.0_C_REAL
      outLattice%facets(i)%area = areaABC
      IF (lenA**2 .GE. (lenB+lenC)*(lenB+lenC)-(lenB*lenC)) THEN
        outLattice%facets(i)%center = outLattice%vertices(primoVertex)%location
      ELSE IF (lenB**2 .GE. (lenC+lenA)*(lenC+lenA)-(lenC*lenA)) THEN
        outLattice%facets(i)%center = outLattice%vertices(secundoVertex)%location
      ELSE IF (lenC**2 .GE. (lenA+lenB)*(lenA+lenB)-(lenA*lenB)) THEN
        outLattice%facets(i)%center = outLattice%vertices(tertioVertex)%location
      ELSE
        baryA = lenA**4 - 2*(lenB**2 - lenC**2)**2 + lenA**2*(lenB**2 + lenC**2 + 4*SQRT(3._C_REAL)*areaABC)
        baryB = lenB**4 - 2*(lenC**2 - lenA**2)**2 + lenB**2*(lenC**2 + lenA**2 + 4*SQRT(3._C_REAL)*areaABC)
        baryC = lenC**4 - 2*(lenA**2 - lenB**2)**2 + lenC**2*(lenA**2 + lenB**2 + 4*SQRT(3._C_REAL)*areaABC)
        outLattice%facets(i)%center = baryA*outLattice%vertices(primoVertex)%location
        outLattice%facets(i)%center = outLattice%facets(i)%center + baryB*outLattice%vertices(secundoVertex)%location
        outLattice%facets(i)%center = outLattice%facets(i)%center + baryC*outLattice%vertices(tertioVertex)%location
        outLattice%facets(i)%center = outLattice%facets(i)%center/(baryA+baryB+baryC)
      END IF
      CALL outLattice%facetReference(primoVertex)%listHead%push(secundoVertex,i)
      CALL outLattice%facetReference(secundoVertex)%listHead%push(tertioVertex,i)
      CALL outLattice%facetReference(tertioVertex)%listHead%push(primoVertex,i)

    END DO
    outLattice%facets(0)%fIndex = 0
    outLattice%facets(0)%isUniverse = .TRUE.
    outLattice%facets(0)%area = 0.0
    ! let center be some garbage value
    ! vertices will have to be handled later
    enrgScale = outLattice%enrgScale
    restBose  = outLattice%restBose
    restFermi = outLattice%restFermi
    DO i = 1,outLattice%numEdges
      outLattice%edges(i)%eIndex = i
      outLattice%edges(i)%is_fake = .FALSE.
      outLattice%edges(i)%boundary = f_shape%segmentmarkerlist(i)
      outLattice%edges(i)%primoVertex = f_shape%segmentlist(2*i-1)
      outLattice%edges(i)%secundoVertex = f_shape%segmentlist(2*i)
      primoVertex = outLattice%edges(i)%primoVertex
      secundoVertex = outLattice%edges(i)%secundoVertex
      edgeA = outLattice%vertices(secundoVertex)%location - outLattice%vertices(primoVertex)%location
      outLattice%edges(i)%center = outLattice%vertices(primoVertex)%location + 0.5_C_REAL*edgeA
      outLattice%edges(i)%length = NORM2(edgeA)
      outLattice%edges(i)%orientation = ATAN2(edgeA(2),edgeA(1))
      IF (outLattice%edges(i)%orientation < 0.0_C_REAL) THEN
        outLattice%edges(i)%orientation = outLattice%edges(i)%orientation + ATAN2(1.0_C_REAL,0.0_C_REAL) ! adds pi
      END IF
      CALL outLattice%edgeHash(primoVertex)%listHead%push(secundoVertex,i)
      CALL outLattice%edgeHash(secundoVertex)%listHead%push(primoVertex,i)
      !!!! ENERGY CALCULATION !!!! fiddling with this may improve how well this simulates reality!
      length    = outLattice%edges(i)%length
      q1        = outLattice%vertices(primoVertex)%qState
      q2        = outLattice%vertices(secundoVertex)%qState
      outLattice%edges(i)%enrgEdge = calculate_edge_energy(enrgScale,restBose,restFermi,q1,q2,length)

      facetIndex = outLattice%facetReference(primoVertex)%listHead%search(secundoVertex)
      outLattice%edges(i)%sinesterFacet = facetIndex
      IF (primoVertex == outLattice%facets(facetIndex)%primoVertex) THEN
        outLattice%facets(facetIndex)%primoEdge = i
      ELSE IF (primoVertex == outLattice%facets(facetIndex)%secundoVertex) THEN
        outLattice%facets(facetIndex)%secundoEdge = i
      ELSE IF (primoVertex == outLattice%facets(facetIndex)%tertioVertex) THEN
        outLattice%facets(facetIndex)%tertioEdge = i
      END IF
      facetIndex = outLattice%facetReference(secundoVertex)%listHead%search(primoVertex)
      outLattice%edges(i)%dexterFacet = facetIndex
    END DO
  END SUBROUTINE transfer_triangle_to_lattice
END MODULE trans_lattice_triangle
