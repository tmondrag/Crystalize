MODULE triangle_output
  IMPLICIT NONE
CONTAINS
  SUBROUTINE output_triangle_to_vtk(f_shapes,outfile)
    USE triangle_c_wrap, only: f_triangulateio
    USE filehandling, only: stderr, safeopen_writenew
    IMPLICIT NONE
    TYPE(f_triangulateio),INTENT(IN)    :: f_shapes
    CHARACTER(len=*),INTENT(IN)         :: outfile
    INTEGER                             :: fd   !file descriptor
    INTEGER                             :: i,j, num_attributes

    fd = safeopen_writenew(outfile)
    ! Write Header Information
    WRITE(fd,'(A)') '# vtk DataFile Version 3.0'
    WRITE(fd,'(A)') 'Translation of PSLG format file'
    WRITE(fd,'(A)') 'ASCII'
    WRITE(fd,'(A)') 'DATASET UNSTRUCTURED_GRID'
    ! Write point header
    WRITE(fd,'(/A,i10,A)') 'POINTS ' ,f_shapes%numberofpoints+f_shapes%numberofholes,' DOUBLE'
    ! Write out the vertices
    DO i=1,f_shapes%numberofpoints
      WRITE(fd,'(F18.12,1X,F18.12,1X,F18.12)') f_shapes%pointlist(2*i-1), f_shapes%pointlist(2*i), 0.
      FLUSH(fd)
    END DO
    ! Write out the hole seeds
    IF(f_shapes%numberofholes > 0) THEN
      DO i=1,f_shapes%numberofholes
        WRITE(fd,'(F18.12,1X,F18.12,1X,F18.12)') f_shapes%holelist(2*i-1), f_shapes%holelist(2*i), 0.
      FLUSH(fd)
      END DO
    END IF
    ! Write out the connectivity header.
    IF(f_shapes%numberoftriangles > 0) THEN
      IF(f_shapes%numberofcorners == 6) THEN
        WRITE(fd,'(/A,1X,I8,1X,I8)') 'CELLS ',f_shapes%numberofsegments+f_shapes%numberoftriangles, &
                                      3*(f_shapes%numberofsegments)+7*(f_shapes%numberoftriangles)
      FLUSH(fd)
      ELSE IF(f_shapes%numberofcorners == 3) THEN
        WRITE(fd,'(/A,1X,I8,1X,I8)') 'CELLS ',(f_shapes%numberofsegments)+(f_shapes%numberoftriangles), &
                                      3*(f_shapes%numberofsegments)+4*(f_shapes%numberoftriangles)
      FLUSH(fd)
      ELSE
       WRITE(stderr,'(A,I2)') 'ERROR: unsupported number of corners on triangle = ', f_shapes%numberofcorners
      END IF
    ELSE
      IF(f_shapes%numberofsegments > 0) THEN
        WRITE(fd,'(/A,1X,I8,1X,I8)') 'CELLS ',f_shapes%numberofsegments, 3*f_shapes%numberofsegments
      FLUSH(fd)
      END IF
    END IF

    ! Write out segment connectivity. Remember, vtk point indices start with 0
    IF(f_shapes%numberofsegments > 0) THEN
      DO i=1,f_shapes%numberofsegments
        WRITE(fd,'(I2,1X,I8,1X,I8)') 2, f_shapes%segmentlist(2*i-1)-1, f_shapes%segmentlist(2*i)-1
      FLUSH(fd)
      END DO
    END IF
    ! Write out triangle connectivity. Remember, vtk point indices start with 0
    ! remember vtk likes its points ordered on a quad triangle so that the first edge is opposite the third corner, the second edge is opposite the first corner, and the third edge is opposite the second corner
    IF(f_shapes%numberoftriangles > 0) THEN
      IF(f_shapes%numberofcorners == 6) THEN
        DO i=1,f_shapes%numberoftriangles
          WRITE(fd,'(I2,1X,I8,1X,I8,1X,I8,1X,I8,1X,I8,1X,I8)') 6, f_shapes%trianglelist(6*i-5)-1, f_shapes%trianglelist(6*i-4)-1,&
            f_shapes%trianglelist(6*i-3)-1, f_shapes%trianglelist(6*i)-1, f_shapes%trianglelist(6*i-2)-1,&
            f_shapes%trianglelist(6*i-1)-1
        END DO
      ELSE
        DO i=1,f_shapes%numberoftriangles
          WRITE(fd,'(I2,1X,I8,1X,I8,1X,I8)') 3, f_shapes%trianglelist(3*i-2)-1, f_shapes%trianglelist(3*i-1)-1, &
            f_shapes%trianglelist(3*i)-1
        END DO
      END IF
    END IF
    ! Write out the cell types: 3 for straight line segments, 5 for simple triangles, and 22 for quadratic triangles
    IF(f_shapes%numberofsegments+f_shapes%numberoftriangles > 0) THEN
      WRITE(fd,'(/A,1X,I8)') 'CELL_TYPES',f_shapes%numberofsegments+f_shapes%numberoftriangles
    END IF
    IF(f_shapes%numberofsegments > 0) THEN
      DO i=1,f_shapes%numberofsegments
        WRITE(fd,'(I2)') 3
      END DO
    END IF
    IF(f_shapes%numberoftriangles > 0) THEN
      IF(f_shapes%numberofcorners == 6) THEN
        DO i=1,f_shapes%numberoftriangles
          WRITE(fd,'(I2)') 22
        END DO
      ELSE
        DO i=1,f_shapes%numberoftriangles
          WRITE(fd,'(I2)') 5
        END DO
      END IF
    END IF
    ! write out point attributes
    ! the number of attributes should be the same for each point
    ! each vertex (hole seeds included) will all have to the same number of scalars defined
    ! define one scalar at least to be the boundary flag
    num_attributes = f_shapes%numberofpointattributes
    WRITE(fd,'(/A,I10)') 'POINT_DATA ',f_shapes%numberofpoints+f_shapes%numberofholes
    WRITE(fd,'(/A,/A)') 'SCALARS vertex_boundary int','LOOKUP_TABLE default'
    DO i=1,f_shapes%numberofpoints
      WRITE(fd,'(I4)') f_shapes%pointmarkerlist(i)
    END DO
    IF(f_shapes%numberofholes > 0) THEN
      DO i=1,f_shapes%numberofholes
        WRITE(fd,'(I1)') 0
      END DO
    END IF
    ! mark off the hole seeds in the point set
    IF(f_shapes%numberofholes > 0) THEN
      WRITE(fd,'(/A,/A)') 'SCALARS hole_seeds int','LOOKUP_TABLE default'
      DO i=1,f_shapes%numberofpoints
        WRITE(fd,'(I1)') 0
      END DO
      DO i=1,f_shapes%numberofholes
        WRITE(fd,'(I1)') 1
      END DO
    END IF
    IF(num_attributes .GT. 0) THEN
      DO j=1,num_attributes
        WRITE(fd,'(/A,I0.2,A,/A)') 'SCALARS attribute_',j,' double','LOOKUP_TABLE default'
        DO i=1,f_shapes%numberofpoints
          WRITE(fd,'(F18.12)') f_shapes%pointattributelist((i-1)+j)
        END DO
        IF(f_shapes%numberofholes > 0) THEN
          ! setting the atrribute value to 0 on hole seeds may mess with some visualizations
          ! but what else could I set it to? any value will introduce false data to the viz
          DO i=1,f_shapes%numberofholes
            WRITE(fd,'(F18.12)') 0.
          END DO
        END IF
      END DO
    END IF
    ! Write out segment and triangle attributes
    IF(f_shapes%numberofsegments+f_shapes%numberoftriangles > 0) THEN
      num_attributes = f_shapes%numberoftriangleattributes
      WRITE(fd,'(/A10,I10)') 'CELL_DATA ',f_shapes%numberofsegments+f_shapes%numberoftriangles
      ! segments may also be bounds
      IF(f_shapes%numberofsegments > 0) THEN
        WRITE(fd,'(/A,/A)') 'SCALARS segment_boundary int','LOOKUP_TABLE default'
        DO i=1,f_shapes%numberofsegments
          WRITE(fd,'(I4)') f_shapes%segmentmarkerlist(i)
        END DO
        IF(f_shapes%numberoftriangles > 0) THEN
          DO i=1,f_shapes%numberoftriangles
            WRITE(fd,'(I4)') 0
          END DO
        END IF
      END IF
      ! write out triangle element attributes
      IF((num_attributes > 0) .AND. (f_shapes%numberoftriangles > 0)) THEN
        DO j=1,num_attributes
          WRITE(fd,'(/A,I0.2,A,/A)') 'SCALARS cell_attribute_',j,' double','LOOKUP_TABLE default'
          IF (f_shapes%numberofsegments > 0) THEN
            DO i=1,f_shapes%numberofsegments
              WRITE(fd,'(F18.12)') 0.
            END DO
          END IF
          DO i=1,f_shapes%numberoftriangles
            WRITE(fd,'(F18.12)') f_shapes%triangleattributelist((i-1)+j)
          END DO
        END DO
      END IF
    END IF
  END SUBROUTINE output_triangle_to_vtk
  
  SUBROUTINE output_triangle_to_poly(f_shapes,outfile_root)
    USE triangle_c_wrap, only: f_triangulateio
    USE filehandling, only: stderr, safeopen_writenew
    IMPLICIT NONE
    TYPE(f_triangulateio),INTENT(IN)    :: f_shapes
    CHARACTER(len=*),INTENT(IN)         :: outfile_root
    INTEGER                             :: pfd,nfd,efd   !file descriptor ints
    INTEGER                             :: i,j,num_attr,num_corn
    CHARACTER(len=256)                  :: outbuffer,outfile_poly,outfile_node,outfile_ele
    
    WRITE(outfile_poly,'(A,A)') TRIM(outfile_root),".poly"
    WRITE(outfile_node,'(A,A)') TRIM(outfile_root),".node"
    WRITE(outfile_ele,'(A,A)') TRIM(outfile_root),".ele"
    
    pfd = safeopen_writenew(outfile_poly)
    nfd = safeopen_writenew(outfile_node)
    IF (f_shapes%numberofpoints > 0) THEN
      WRITE(pfd,'(I8,1X,I1,1X,I8,1X,I1)')0,2,f_shapes%numberofpointattributes,1
      WRITE(nfd,'(I8,1X,I1,1X,I8,1X,I1)')f_shapes%numberofpoints,2,f_shapes%numberofpointattributes,1
      IF (f_shapes%numberofpointattributes > 0) THEN
      num_attr = f_shapes%numberofpointattributes
        DO i = 1,f_shapes%numberofpoints
          WRITE(outbuffer,'(I8,1X,F14.8,1X,F14.8)')i,f_shapes%pointlist(2*i-1),f_shapes%pointlist(2*i)
          DO j = 1,num_attr
            WRITE(outbuffer,'(A,1X,F14.8)') TRIM(outbuffer),f_shapes%pointattributelist(num_attr*(i-1)+j)
          END DO
          WRITE(outbuffer,'(A,1X,I1)') TRIM(outbuffer),f_shapes%pointmarkerlist(i)
          WRITE(nfd,'(A)')outbuffer
        END DO
      ELSE
        DO i = 1,f_shapes%numberofpoints
          WRITE(outbuffer,'(I8,1X,F14.8,1X,F14.8)')i,f_shapes%pointlist(2*i-1),f_shapes%pointlist(2*i)
          WRITE(outbuffer,'(A,1X,I1)') TRIM(outbuffer),f_shapes%pointmarkerlist(i)
          WRITE(nfd,'(A)')outbuffer
        END DO
      END IF
      WRITE(pfd,'(A)')""
      IF (f_shapes%numberofsegments > 0) THEN
        WRITE(pfd,'(I8,1X,I1)') f_shapes%numberofsegments,1
        DO i = 1,f_shapes%numberofsegments
          WRITE(pfd,'(I8,1X,I8,1X,I8,1X,I1)') i,f_shapes%segmentlist(2*i-1),f_shapes%segmentlist(2*i),1
        END DO
      ELSE
        WRITE(pfd,'(I8,1X,I1)') 0,0
      END IF
      WRITE(pfd,'(A)')""
      IF (f_shapes%numberofholes > 0) THEN
        WRITE(pfd,'(I8)') f_shapes%numberofholes
        DO i = 1,f_shapes%numberofholes
          WRITE(pfd,'(I8,1X,F14.8,1X,F14.8)') i,f_shapes%holelist(2*i-1),f_shapes%holelist(2*i)
        END DO
      ELSE
        WRITE(pfd,'(I8)') 0
      END IF
      WRITE(pfd,'(A)')""
      IF (f_shapes%numberofregions > 0) THEN
        WRITE(pfd,'(I8)') f_shapes%numberofregions
        DO i = 1,f_shapes%numberofregions
          WRITE(pfd,'(I8,1X,F14.8,1X,F14.8,1X,F14.8,1X,F14.8)') i,f_shapes%regionlist(4*i-3),f_shapes%regionlist(4*i-2), &
                                                                 f_shapes%regionlist(4*i-1),f_shapes%regionlist(4*i)
        END DO
      END IF
      IF (f_shapes%numberoftriangles > 0) THEN
        efd = safeopen_writenew(outfile_ele)
        WRITE(efd,'(I8,1X,I1,1X,I8)')f_shapes%numberoftriangles,f_shapes%numberofcorners,f_shapes%numberoftriangleattributes
        IF (f_shapes%numberoftriangleattributes > 0) THEN
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
            WRITE(efd,'(A)')outbuffer
          END DO
        ELSE
          num_corn = f_shapes%numberofcorners
          DO i = 1,f_shapes%numberoftriangles
            WRITE(outbuffer,'(I8)')i
            DO j= num_corn - 1,0,-1
              WRITE(outbuffer,'(A,1X,I8)') TRIM(outbuffer),f_shapes%trianglelist(num_corn*i-j)
            END DO
            WRITE(efd,'(A)')outbuffer
          END DO
        END IF
        CLOSE(efd)
      END IF
    ELSE
      WRITE(stderr,'(A)') "Misformed triangulateio data, nothing printed to file"
    END IF
    CLOSE(pfd)
    CLOSE(nfd)
  END SUBROUTINE output_triangle_to_poly
END MODULE triangle_output
