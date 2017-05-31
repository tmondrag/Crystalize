MODULE lattice_output
  IMPLICIT NONE
CONTAINS
  SUBROUTINE output_lattice_to_vtk(inLattice,filename)
    USE basictypes, only: Lattice
    USE mFileHandling, only: stderr, FILE
    IMPLICIT NONE

    TYPE(Lattice),INTENT(IN)            :: inLattice
    CHARACTER(len=*),INTENT(IN)         :: filename
    TYPE(FILE)                          :: outputfile
    INTEGER                             :: i,j, num_attributes
    INTEGER                             :: totPointCount,totCellCount, cellDataCount

    totPointCount = inLattice%numVertices
    totCellCount  = inLattice%numEdges+inLattice%numFacets!+totPointCount

    CALL outputFile%openWriteNew(filename)
    ! Write Header Information
    WRITE(outputFile%getFUnit(),'(A)') '# vtk DataFile Version 3.0'
    WRITE(outputFile%getFUnit(),'(A)') 'Vtk file produced from Lattice Data'
    WRITE(outputFile%getFUnit(),'(A)') 'ASCII'
    WRITE(outputFile%getFUnit(),'(A)') 'DATASET UNSTRUCTURED_GRID'
    ! Write point header
    WRITE(outputFile%getFUnit(),'(/A,i10,A)') 'POINTS ' ,totPointCount,' DOUBLE'
    ! Write out the vertices
    DO i=1,totPointCount
      WRITE(outputFile%getFUnit(),'(F18.12,1X,F18.12,1X,F18.12)') inLattice%vertices(i)%location(1), &
                                                                  inLattice%vertices(i)%location(2), 0.
      FLUSH(outputFile%getFUnit())
    END DO
    ! Write out the connectivity header.
    cellDataCount = 3*inLattice%numEdges+4*inLattice%numFacets!+2*totPointCount
    WRITE(outputFile%getFUnit(),'(/A,1X,I8,1X,I8)') 'CELLS ',totCellCount, cellDataCount
    FLUSH(outputFile%getFUnit())
    ! Write out segment connectivity. Remember, vtk point indices start with 0
    IF(inLattice%numEdges > 0) THEN
      DO i=1,inLattice%numEdges
        WRITE(outputFile%getFUnit(),'(I2,1X,I8,1X,I8)') 2, inLattice%edges(i)%primoVertex-1, inLattice%edges(i)%secundoVertex-1
      FLUSH(outputFile%getFUnit())
      END DO
    END IF
    ! Write out triangle connectivity. Remember, vtk point indices start with 0
    ! remember vtk likes its points ordered on a quad triangle so that the first edge is opposite the third corner, the second edge is opposite the first corner, and the third edge is opposite the second corner
    IF(inLattice%numFacets > 0) THEN
      DO i=1,inLattice%numFacets
        WRITE(outputFile%getFUnit(),'(I2,1X,I8,1X,I8,1X,I8)') 3, inLattice%facets(i)%primoVertex-1, &
                                                                 inLattice%facets(i)%secundoVertex-1, &
                                                                 inLattice%facets(i)%tertioVertex-1
      END DO
    END IF
    ! Write out the naked points as vtk_vertex cells
    ! DO i = 1, totPointCount
    !   WRITE(outputFile%getFUnit(),'(I2,1X,I8)') 1,i - 1
    ! END DO

    ! Write out the cell types: 3 for straight line segments, 5 for simple triangles, and 22 for quadratic triangles, and 1 for points
    WRITE(outputFile%getFUnit(),'(/A,1X,I8)') 'CELL_TYPES',totCellCount

    IF(inLattice%numEdges > 0) THEN
      DO i=1,inLattice%numEdges
        WRITE(outputFile%getFUnit(),'(I2,1X)',ADVANCE="no") 3
      END DO
    END IF
    IF(inLattice%numFacets > 0) THEN
      DO i=1,inLattice%numFacets
        WRITE(outputFile%getFUnit(),'(I2,1X)',ADVANCE="no") 5
      END DO
    END IF
    ! DO i = 1, totPointCount
    !   WRITE(outputFile%getFUnit(),'(I2,1X)',ADVANCE="no") 1
    ! END DO
    WRITE(outputFile%getFUnit(),'(A)') ""
    ! Write out point attributes
    WRITE(outputFile%getFUnit(),'(/A,I10)') 'POINT_DATA ',totPointCount
    ! Boundaries
    WRITE(outputFile%getFUnit(),'(/A,/A)') 'SCALARS vertex_boundary int','LOOKUP_TABLE default'
    DO i=1,totPointCount
      WRITE(outputFile%getFUnit(),'(I4,1X)',ADVANCE="no") inLattice%vertices(i)%boundary
    END DO
    ! q_state
    WRITE(outputFile%getFUnit(),'(/A,/A)') 'SCALARS q_state int','LOOKUP_TABLE default'
    DO i=1,totPointCount
      WRITE(outputFile%getFUnit(),'(I4,1X)',ADVANCE="no") inLattice%vertices(i)%qState
    END DO
    ! q_state
    WRITE(outputFile%getFUnit(),'(/A,/A)') 'SCALARS grain_id int','LOOKUP_TABLE default'
    DO i=1,totPointCount
      WRITE(outputFile%getFUnit(),'(I4,1X)',ADVANCE="no") inLattice%vertices(i)%enumer
    END DO
    ! Write out segment and triangle attributes
    IF(totCellCount > 0) THEN
      WRITE(outputFile%getFUnit(),'(/A10,I10)') 'CELL_DATA ',totCellCount
      ! segments may also be boundaries
      IF(inLattice%numEdges > 0) THEN
        WRITE(outputFile%getFUnit(),'(/A,/A)') 'SCALARS segment_boundary int','LOOKUP_TABLE default'
        DO i=1,inLattice%numEdges
          WRITE(outputFile%getFUnit(),'(I4,1X)',ADVANCE="no") inLattice%edges(i)%boundary
        END DO
        IF(inLattice%numFacets > 0) THEN
          DO i=1,inLattice%numFacets
            WRITE(outputFile%getFUnit(),'(I4,1X)',ADVANCE="no") 0
          END DO
        END IF
        ! DO i=1,totPointCount
        !   WRITE(outputFile%getFUnit(),'(I4,1X)',ADVANCE="no") 0
        ! END DO
        WRITE(outputFile%getFUnit(),'(A)') ""
      END IF
      ! edges in a lattice have brokeness states based on wheather the grain boudary crosses them
      IF(inLattice%numEdges > 0) THEN
        WRITE(outputFile%getFUnit(),'(/A,/A)') 'SCALARS broken_edges int','LOOKUP_TABLE default'
        DO i=1,inLattice%numEdges
          IF(inLattice%edges(i)%isBroken) THEN
            WRITE(outputFile%getFUnit(),'(I4,1X)',ADVANCE="no") 1
          ELSE
            WRITE(outputFile%getFUnit(),'(I4,1X)',ADVANCE="no") 0
          END IF
        END DO
        IF(inLattice%numFacets > 0) THEN
          DO i=1,inLattice%numFacets
            WRITE(outputFile%getFUnit(),'(I4,1X)',ADVANCE="no") 0
          END DO
        END IF
        ! DO i=1,totPointCount
        !   WRITE(outputFile%getFUnit(),'(I4,1X)',ADVANCE="no") 0
        ! END DO
        WRITE(outputFile%getFUnit(),'(A)') ""
      END IF
      ! the internal crystal enery is due to interation between neighbors. visualize the energy
      IF(inLattice%numEdges > 0) THEN
        WRITE(outputFile%getFUnit(),'(/A,/A)') 'SCALARS edge_energy float','LOOKUP_TABLE default'
        DO i=1,inLattice%numEdges
          WRITE(outputFile%getFUnit(),'(F18.12,1X)',ADVANCE="no") inLattice%edges(i)%enrgEdge
        END DO
        IF(inLattice%numFacets > 0) THEN
          DO i=1,inLattice%numFacets
            WRITE(outputFile%getFUnit(),'(F18.12,1X)',ADVANCE="no") 0.
          END DO
        END IF
        ! DO i=1,totPointCount
        !   WRITE(outputFile%getFUnit(),'(F18.12,1X)',ADVANCE="no") 0.
        ! END DO
        WRITE(outputFile%getFUnit(),'(A)') ""
      END IF
      ! Triangle properties
      IF(inLattice%numFacets > 0) THEN
        WRITE(outputFile%getFUnit(),'(/A,/A)') 'SCALARS broken_edges int','LOOKUP_TABLE default'
        IF(inLattice%numEdges > 0) THEN
          DO i=1,inLattice%numEdges
            WRITE(outputFile%getFUnit(),'(I4,1X)',ADVANCE="no") 0
          END DO
        END IF
        DO i=1,inLattice%numFacets
          WRITE(outputFile%getFUnit(),'(I4,1X)',ADVANCE="no") inLattice%facets(i)%numBrokenEdges
        END DO
        ! DO i=1,totPointCount
        !   WRITE(outputFile%getFUnit(),'(I4,1X)',ADVANCE="no") 0
        ! END DO
        WRITE(outputFile%getFUnit(),'(A)') ""
      END IF
    END IF
  END SUBROUTINE output_lattice_to_vtk
END MODULE lattice_output
