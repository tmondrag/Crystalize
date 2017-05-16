MODULE triangle_output
  IMPLICIT NONE
CONTAINS
  SUBROUTINE output_triangle_to_vtk(f_shapes,filename)
    USE triangle_c_wrap, only: f_triangulateio
    USE mFileHandling, only: stderr, FILE
    IMPLICIT NONE
    TYPE(f_triangulateio),INTENT(IN)    :: f_shapes
    CHARACTER(len=*),INTENT(IN)         :: filename
    TYPE(FILE)                          :: outputfile
    INTEGER                             :: i,j, num_attributes
    INTEGER                             :: totPointCount,totCellCount, cellDataCount

    totPointCount = f_shapes%numberofpoints+f_shapes%numberofholes
    totCellCount  = f_shapes%numberofsegments+f_shapes%numberoftriangles+totPointCount

    CALL outputFile%openWriteNew(filename)
    ! Write Header Information
    WRITE(outputFile%getFUnit(),'(A)') '# vtk DataFile Version 3.0'
    WRITE(outputFile%getFUnit(),'(A)') 'Translation of PSLG format file'
    WRITE(outputFile%getFUnit(),'(A)') 'ASCII'
    WRITE(outputFile%getFUnit(),'(A)') 'DATASET UNSTRUCTURED_GRID'
    ! Write point header
    WRITE(outputFile%getFUnit(),'(/A,i10,A)') 'POINTS ' ,totPointCount,' DOUBLE'
    ! Write out the vertices
    DO i=1,f_shapes%numberofpoints
      WRITE(outputFile%getFUnit(),'(F18.12,1X,F18.12,1X,F18.12)') f_shapes%pointlist(2*i-1), f_shapes%pointlist(2*i), 0.
      FLUSH(outputFile%getFUnit())
    END DO
    ! Write out the hole seeds
    IF(f_shapes%numberofholes > 0) THEN
      DO i=1,f_shapes%numberofholes
        WRITE(outputFile%getFUnit(),'(F18.12,1X,F18.12,1X,F18.12)') f_shapes%holelist(2*i-1), f_shapes%holelist(2*i), 0.
      FLUSH(outputFile%getFUnit())
      END DO
    END IF
    ! Write out the connectivity header.
    cellDataCount = 3*f_shapes%numberofsegments+2*totPointCount
    IF(f_shapes%numberoftriangles > 0) THEN
      IF(f_shapes%numberofcorners == 6) THEN
        cellDataCount = cellDataCount+7*(f_shapes%numberoftriangles)
      ELSE IF(f_shapes%numberofcorners == 3) THEN
        cellDataCount = cellDataCount+4*(f_shapes%numberoftriangles)
      ELSE
       WRITE(stderr,'(A,I2)') 'ERROR: unsupported number of corners on triangle = ', f_shapes%numberofcorners
       STOP
      END IF
    END IF

    WRITE(outputFile%getFUnit(),'(/A,1X,I8,1X,I8)') 'CELLS ',totCellCount, cellDataCount
    FLUSH(outputFile%getFUnit())

    ! Write out segment connectivity. Remember, vtk point indices start with 0
    IF(f_shapes%numberofsegments > 0) THEN
      DO i=1,f_shapes%numberofsegments
        WRITE(outputFile%getFUnit(),'(I2,1X,I8,1X,I8)') 2, f_shapes%segmentlist(2*i-1)-1, f_shapes%segmentlist(2*i)-1
      FLUSH(outputFile%getFUnit())
      END DO
    END IF
    ! Write out triangle connectivity. Remember, vtk point indices start with 0
    ! remember vtk likes its points ordered on a quad triangle so that the first edge is opposite the third corner, the second edge is opposite the first corner, and the third edge is opposite the second corner
    IF(f_shapes%numberoftriangles > 0) THEN
      IF(f_shapes%numberofcorners == 6) THEN
        DO i=1,f_shapes%numberoftriangles
          WRITE(outputFile%getFUnit(),'(I2,1X,I8,1X,I8,1X,I8,1X,I8,1X,I8,1X,I8)') 6, f_shapes%trianglelist(6*i-5)-1, &
            f_shapes%trianglelist(6*i-4)-1, f_shapes%trianglelist(6*i-3)-1, f_shapes%trianglelist(6*i)-1, &
            f_shapes%trianglelist(6*i-2)-1, f_shapes%trianglelist(6*i-1)-1
        END DO
      ELSE
        DO i=1,f_shapes%numberoftriangles
          WRITE(outputFile%getFUnit(),'(I2,1X,I8,1X,I8,1X,I8)') 3, f_shapes%trianglelist(3*i-2)-1, f_shapes%trianglelist(3*i-1)-1, &
            f_shapes%trianglelist(3*i)-1
        END DO
      END IF
    END IF
    ! Write out the naked points as vtk_vertex cells
    DO i = 1, totPointCount
      WRITE(outputFile%getFUnit(),'(I2,1X,I8)') 1,i - 1
    END DO

    ! Write out the cell types: 3 for straight line segments, 5 for simple triangles, and 22 for quadratic triangles, and 1 for points
    WRITE(outputFile%getFUnit(),'(/A,1X,I8)') 'CELL_TYPES',totCellCount

    IF(f_shapes%numberofsegments > 0) THEN
      DO i=1,f_shapes%numberofsegments
        WRITE(outputFile%getFUnit(),'(I2,1X)',ADVANCE="no") 3
      END DO
    END IF
    IF(f_shapes%numberoftriangles > 0) THEN
      IF(f_shapes%numberofcorners == 6) THEN
        DO i=1,f_shapes%numberoftriangles
          WRITE(outputFile%getFUnit(),'(I2,1X)',ADVANCE="no") 22
        END DO
      ELSE
        DO i=1,f_shapes%numberoftriangles
          WRITE(outputFile%getFUnit(),'(I2,1X)',ADVANCE="no") 5
        END DO
      END IF
    END IF
    DO i = 1, totPointCount
      WRITE(outputFile%getFUnit(),'(I2,1X)',ADVANCE="no") 1
    END DO
    WRITE(outputFile%getFUnit(),'(A)') ""

    ! write out point attributes
    ! the number of attributes should be the same for each point
    ! each vertex (hole seeds included) will all have to the same number of scalars defined
    ! define one scalar at least to be the boundary flag
    num_attributes = f_shapes%numberofpointattributes
    WRITE(outputFile%getFUnit(),'(/A,I10)') 'POINT_DATA ',totPointCount
    WRITE(outputFile%getFUnit(),'(/A,/A)') 'SCALARS vertex_boundary int','LOOKUP_TABLE default'
    DO i=1,f_shapes%numberofpoints
      WRITE(outputFile%getFUnit(),'(I4,1X)',ADVANCE="no") f_shapes%pointmarkerlist(i)
    END DO
    IF(f_shapes%numberofholes > 0) THEN
      DO i=1,f_shapes%numberofholes
        WRITE(outputFile%getFUnit(),'(I4,1X)',ADVANCE="no") 0
      END DO
    END IF
    WRITE(outputFile%getFUnit(),'(A)') ""
    ! mark off the hole seeds in the point set
    IF(f_shapes%numberofholes > 0) THEN
      WRITE(outputFile%getFUnit(),'(/A,/A)') 'SCALARS hole_seeds int','LOOKUP_TABLE default'
      DO i=1,f_shapes%numberofpoints
        WRITE(outputFile%getFUnit(),'(I1,1X)',ADVANCE="no") 0
      END DO
      DO i=1,f_shapes%numberofholes
        WRITE(outputFile%getFUnit(),'(I1,1X)',ADVANCE="no") 1
      END DO
      WRITE(outputFile%getFUnit(),'(A)') ""
    END IF
    IF(num_attributes .GT. 0) THEN
      DO j=1,num_attributes
        IF (j==1) THEN
          WRITE(outputFile%getFUnit(),'(/A,/A)') 'SCALARS Q-state double','LOOKUP_TABLE default'
        ELSE IF (j==2) THEN
          WRITE(outputFile%getFUnit(),'(/A,/A)') 'SCALARS grain_id double','LOOKUP_TABLE default'
        ELSE
          WRITE(outputFile%getFUnit(),'(/A,I0.2,A,/A)') 'SCALARS attribute_',j,' double','LOOKUP_TABLE default'
        END IF
        DO i=1,f_shapes%numberofpoints
          WRITE(outputFile%getFUnit(),'(F18.12,1X)',ADVANCE="no") f_shapes%pointattributelist((i-1)+j)
        END DO
        IF(f_shapes%numberofholes > 0) THEN
          ! setting the atrribute value to 0 on hole seeds may mess with some visualizations
          ! but what else could I set it to? any value will introduce false data to the viz
          DO i=1,f_shapes%numberofholes
            WRITE(outputFile%getFUnit(),'(F18.12,1X)',ADVANCE="no") 0.
          END DO
        END IF
        WRITE(outputFile%getFUnit(),'(A)') ""
      END DO
    END IF
    ! Write out segment and triangle attributes
    IF(totCellCount > 0) THEN
      num_attributes = f_shapes%numberoftriangleattributes
      WRITE(outputFile%getFUnit(),'(/A10,I10)') 'CELL_DATA ',totCellCount
      ! segments may also be bounds
      IF(f_shapes%numberofsegments > 0) THEN
        WRITE(outputFile%getFUnit(),'(/A,/A)') 'SCALARS segment_boundary int','LOOKUP_TABLE default'
        DO i=1,f_shapes%numberofsegments
          WRITE(outputFile%getFUnit(),'(I4,1X)',ADVANCE="no") f_shapes%segmentmarkerlist(i)
        END DO
        IF(f_shapes%numberoftriangles > 0) THEN
          DO i=1,f_shapes%numberoftriangles
            WRITE(outputFile%getFUnit(),'(I4,1X)',ADVANCE="no") 0
          END DO
        END IF
        DO i=1,totPointCount
          WRITE(outputFile%getFUnit(),'(I4,1X)',ADVANCE="no") 0
        END DO
        WRITE(outputFile%getFUnit(),'(A)') ""
      END IF
      ! write out triangle element attributes
      IF((num_attributes > 0) .AND. (f_shapes%numberoftriangles > 0)) THEN
        DO j=1,num_attributes
          WRITE(outputFile%getFUnit(),'(/A,I0.2,A,/A)') 'SCALARS cell_attribute_',j,' double','LOOKUP_TABLE default'
          IF (f_shapes%numberofsegments > 0) THEN
            DO i=1,f_shapes%numberofsegments
              WRITE(outputFile%getFUnit(),'(F18.12,1X)',ADVANCE="no") 0.
            END DO
          END IF
          DO i=1,f_shapes%numberoftriangles
            WRITE(outputFile%getFUnit(),'(F18.12,1X)',ADVANCE="no") f_shapes%triangleattributelist((i-1)+j)
          END DO
          DO i=1,totPointCount
            WRITE(outputFile%getFUnit(),'(F18.12,1X)',ADVANCE="no") 0.
          END DO
          WRITE(outputFile%getFUnit(),'(A)') ""
        END DO
      END IF
    END IF
  END SUBROUTINE output_triangle_to_vtk

  SUBROUTINE output_triangle_to_poly(f_shapes,outfile_root)
    USE triangle_c_wrap, only: f_triangulateio
    USE mFileHandling, only: stderr, FILE
    IMPLICIT NONE
    TYPE(f_triangulateio),INTENT(IN)    :: f_shapes
    CHARACTER(len=*),INTENT(IN)         :: outfile_root
    TYPE(FILE)                          :: polyfile, nodefile, elefile
    INTEGER                             :: i,j,num_attr,num_corn
    CHARACTER(len=256)                  :: outbuffer,outfile_poly,outfile_node,outfile_ele

    WRITE(outfile_poly,'(A,A)') TRIM(outfile_root),".poly"
    WRITE(outfile_node,'(A,A)') TRIM(outfile_root),".node"
    WRITE(outfile_ele,'(A,A)') TRIM(outfile_root),".ele"

    CALL polyfile%openWriteNew(outfile_poly)
    CALL nodefile%openWriteNew(outfile_node)
    IF (f_shapes%numberofpoints > 0) THEN
      WRITE(polyfile%getFUnit(),'(A)') '# VERTICES '
      WRITE(polyfile%getFUnit(),'(A,I8,A,I1,A,I8,A)') "# Declare ",f_shapes%numberofpoints," vertices, ", &
                    2, " dimensions, ",f_shapes%numberofpointattributes, " attributes per vertex, and 1 boundary marker per vertex"
      WRITE(polyfile%getFUnit(),'(A,A)') '# Actual vertex data is in ', outfile_node
      WRITE(polyfile%getFUnit(),'(I8,1X,I1,1X,I8,1X,I1)')0,2,f_shapes%numberofpointattributes,1
      WRITE(nodefile%getFUnit(),'(A,I8,A,I1,A,I8,A)') "# Declare ",f_shapes%numberofpoints," vertices, ", &
                    2, " dimensions, ",f_shapes%numberofpointattributes, " attributes per vertex, and 1 boundary marker per vertex"
      WRITE(nodefile%getFUnit(),'(I8,1X,I1,1X,I8,1X,I1)')f_shapes%numberofpoints,2,f_shapes%numberofpointattributes,1
      IF (f_shapes%numberofpointattributes > 0) THEN
        num_attr = f_shapes%numberofpointattributes
        IF (f_shapes%numberofpointattributes == 1) THEN
          WRITE(nodefile%getFUnit(),'(/A)') "# <index> <x> <y> <Q-state> <boundary marker>"
          DO i = 1,f_shapes%numberofpoints
            WRITE(outbuffer,'(I8,1X,F14.8,1X,F14.8)')i,f_shapes%pointlist(2*i-1),f_shapes%pointlist(2*i)
            WRITE(outbuffer,'(A,1X,I8)') TRIM(outbuffer),INT(f_shapes%pointattributelist(i))
            WRITE(outbuffer,'(A,1X,I1)') TRIM(outbuffer),f_shapes%pointmarkerlist(i)
            WRITE(nodefile%getFUnit(),'(A)')outbuffer
          END DO
        ELSEIF (f_shapes%numberofpointattributes == 2) THEN
          WRITE(nodefile%getFUnit(),'(/A)') "# <index> <x> <y> <Q-state> <grain-id> <boundary marker>"
          DO i = 1,f_shapes%numberofpoints
            WRITE(outbuffer,'(I8,1X,F14.8,1X,F14.8)')i,f_shapes%pointlist(2*i-1),f_shapes%pointlist(2*i)
            WRITE(outbuffer,'(A,1X,I8)') TRIM(outbuffer),INT(f_shapes%pointattributelist(2*i-1))
            WRITE(outbuffer,'(A,1X,I8)') TRIM(outbuffer),INT(f_shapes%pointattributelist(2*i))
            WRITE(outbuffer,'(A,1X,I1)') TRIM(outbuffer),f_shapes%pointmarkerlist(i)
            WRITE(nodefile%getFUnit(),'(A)')outbuffer
          END DO
        ELSE IF (f_shapes%numberofpointattributes > 2) THEN
          WRITE(nodefile%getFUnit(),'(/A)') "# <index> <x> <y> <Q-state> <grain-id> <other attributes> <boundary marker>"
          DO i = 1,f_shapes%numberofpoints
            WRITE(outbuffer,'(I8,1X,F14.8,1X,F14.8)')i,f_shapes%pointlist(2*i-1),f_shapes%pointlist(2*i)
            WRITE(outbuffer,'(A,1X,I8)') TRIM(outbuffer),INT(f_shapes%pointattributelist(num_attr*i-num_attr+1))
            WRITE(outbuffer,'(A,1X,I8)') TRIM(outbuffer),INT(f_shapes%pointattributelist(num_attr*i-num_attr+2))
            DO j = 3,num_attr
              WRITE(outbuffer,'(A,1X,F14.8)') TRIM(outbuffer),f_shapes%pointattributelist(num_attr*(i-1)+j)
            END DO
            WRITE(outbuffer,'(A,1X,I1)') TRIM(outbuffer),f_shapes%pointmarkerlist(i)
            WRITE(nodefile%getFUnit(),'(A)')outbuffer
          END DO
        END IF

      ELSE
        DO i = 1,f_shapes%numberofpoints
          WRITE(outbuffer,'(I8,1X,F14.8,1X,F14.8)')i,f_shapes%pointlist(2*i-1),f_shapes%pointlist(2*i)
          WRITE(outbuffer,'(A,1X,I1)') TRIM(outbuffer),f_shapes%pointmarkerlist(i)
          WRITE(nodefile%getFUnit(),'(A)')outbuffer
        END DO
      END IF
      WRITE(polyfile%getFUnit(),'(/A)')"# LINE SEGMENTS"
      IF (f_shapes%numberofsegments > 0) THEN
        WRITE(polyfile%getFUnit(),'(A,I8,A,I1,A)') '# Declare ',f_shapes%numberofsegments,' line segments with ', &
                                                  1, ' boundary marker each'
        WRITE(polyfile%getFUnit(),'(I8,1X,I1)') f_shapes%numberofsegments,1
        WRITE(nodefile%getFUnit(),'(/A)') "# <index> <vertex_id 1> <vertex_id 2> <boundary marker>"
        DO i = 1,f_shapes%numberofsegments
          WRITE(polyfile%getFUnit(),'(I8,1X,I8,1X,I8,1X,I1)') i,f_shapes%segmentlist(2*i-1),f_shapes%segmentlist(2*i),1
        END DO
      ELSE
        WRITE(polyfile%getFUnit(),'(I8,1X,I1)') 0,0
      END IF
      WRITE(polyfile%getFUnit(),'(/A)')"# HOLE SEEDS"
      WRITE(polyfile%getFUnit(),'(A)') "# Declare ", f_shapes%numberofholes, "holes"
      IF (f_shapes%numberofholes > 0) THEN
        WRITE(polyfile%getFUnit(),'(I8)') f_shapes%numberofholes
        WRITE(polyfile%getFUnit(),'(/A)') "# <index> <x> <y>"
        DO i = 1,f_shapes%numberofholes
          WRITE(polyfile%getFUnit(),'(I8,1X,F14.8,1X,F14.8)') i,f_shapes%holelist(2*i-1),f_shapes%holelist(2*i)
        END DO
      ELSE
        WRITE(polyfile%getFUnit(),'(I8)') 0
      END IF
      WRITE(polyfile%getFUnit(),'(/A)')"# REGIONS"
      IF (f_shapes%numberofregions > 0) THEN
        WRITE(polyfile%getFUnit(),'(A,I8,A)')"# Declare ", f_shapes%numberofregions, " regions"
        WRITE(polyfile%getFUnit(),'(I8)') f_shapes%numberofregions
        WRITE(polyfile%getFUnit(),'(A,I8,A)')"# <index> <seed x> <seed y> <attribute> <maximum area>"
        DO i = 1,f_shapes%numberofregions
          WRITE(polyfile%getFUnit(),'(I8,1X,F14.8,1X,F14.8,1X,F14.8,1X,F14.8)') i,f_shapes%regionlist(4*i-3), &
                                                f_shapes%regionlist(4*i-2),f_shapes%regionlist(4*i-1),f_shapes%regionlist(4*i)
        END DO
      END IF
      IF (f_shapes%numberoftriangles > 0) THEN
        CALL elefile%openWriteNew(outfile_ele)
        WRITE(elefile%getFUnit(),'(/A)') "# TRIANGLES"
        WRITE(elefile%getFUnit(),'(A,I8,A,I1,A,I8,A)')"# Declare ",f_shapes%numberoftriangles," triangles made of ",&
                  f_shapes%numberofcorners, " vertices and with ",f_shapes%numberoftriangleattributes, " attributes each"
        WRITE(elefile%getFUnit(),'(I8,1X,I1,1X,I8)')f_shapes%numberoftriangles,f_shapes%numberofcorners, &
                  f_shapes%numberoftriangleattributes
        IF (f_shapes%numberoftriangleattributes > 0) THEN
          WRITE(elefile%getFUnit(),'(/A)') "# <index> <vertex_ids> <attributes>"
          num_attr = f_shapes%numberoftriangleattributes
          num_corn = f_shapes%numberofcorners
          DO i = 1,f_shapes%numberoftriangles
            WRITE(outbuffer,'(I8)')i
            DO j= num_corn - 1,0,-1
              WRITE(outbuffer,'(A,1X,I8)') TRIM(outbuffer),f_shapes%trianglelist(num_corn*i-j)
            END DO
            DO j = 1,num_attr
              WRITE(outbuffer,'(A,1X,F14.8)') TRIM(outbuffer),f_shapes%triangleattributelist(num_attr*(i-1)+j)
            END DO
            WRITE(elefile%getFUnit(),'(A)')outbuffer
          END DO
        ELSE
          WRITE(elefile%getFUnit(),'(/A)') "# <index> <vertex_ids>"
          num_corn = f_shapes%numberofcorners
          DO i = 1,f_shapes%numberoftriangles
            WRITE(outbuffer,'(I8)')i
            DO j= num_corn - 1,0,-1
              WRITE(outbuffer,'(A,1X,I8)') TRIM(outbuffer),f_shapes%trianglelist(num_corn*i-j)
            END DO
            WRITE(elefile%getFUnit(),'(A)')outbuffer
          END DO
        END IF
        CALL elefile%close()
      END IF
    ELSE
      WRITE(stderr,'(A)') "Misformed triangulateio data, nothing printed to file"
    END IF
    CALL polyfile%close()
    CALL nodefile%close()
  END SUBROUTINE output_triangle_to_poly
END MODULE triangle_output
