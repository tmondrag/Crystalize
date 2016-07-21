MODULE triangle_c_wrap
  USE, INTRINSIC:: ISO_C_BINDING, only: C_PTR, C_INT, C_FLOAT, C_CHAR, C_DOUBLE
  IMPLICIT NONE
  
  INTEGER,PARAMETER :: C_REAL = C_DOUBLE       ! triangle.c can switch its REAL type between FLOAT and DOUBLE 
                                               ! based on compilation options. match it here

  ! Fortran mirror of triangulateio struct in triangle_lib/triangle.h
  TYPE, BIND(C) :: triangulateio
    TYPE(C_PTR) pointlist                                               !* In / out C_REAL ARRAY*!
    TYPE(C_PTR) pointattributelist                                      !* In / out C_REAL ARRAY*!
    TYPE(C_PTR) pointmarkerlist                                         !* In / out INT ARRAY  *!
    INTEGER(C_INT) numberofpoints                                       !* In / out            *!
    INTEGER(C_INT) numberofpointattributes                              !* In / out            *!

    TYPE(C_PTR) trianglelist                                            !* In / out INT ARRAY  *!
    TYPE(C_PTR) triangleattributelist                                   !* In / out C_REAL ARRAY*!
    TYPE(C_PTR) trianglearealist                                        !* In only C_REAL ARRAY *!
    TYPE(C_PTR) neighborlist                                            !* Out only INT ARRAY  *!
    INTEGER(C_INT) numberoftriangles                                    !* In / out            *!
    INTEGER(C_INT) numberofcorners                                      !* In / out            *!
    INTEGER(C_INT) numberoftriangleattributes                           !* In / out            *!

    TYPE(C_PTR) segmentlist                                             !* In / out INT ARRAY  *!
    TYPE(C_PTR) segmentmarkerlist                                       !* In / out INT ARRAY  *!
    INTEGER(C_INT) numberofsegments                                     !* In / out            *!

    TYPE(C_PTR) holelist                                                !* In / pointer to C_REAL array copied out *!
    INTEGER(C_INT) numberofholes                                        !* In / copied out     *!

    TYPE(C_PTR) regionlist                                              !* In / pointer to C_REAL array copied out *!
    INTEGER(C_INT) numberofregions                                      !* In / copied out     *!

    TYPE(C_PTR) edgelist                                                !* Out only INT ARRAY  *!
    TYPE(C_PTR) edgemarkerlist                                          !* Not used with Voronoi diagram; out only INT ARRAY   *!
    TYPE(C_PTR) normlist                                                !* Used only with Voronoi diagram; out only C_REAL ARRAY*!
    INTEGER(C_INT) numberofedges                                        !* Out only            *!
  END TYPE triangulateio

  ! version of trianglateio struct that is easier to use in Fortran
  TYPE                     :: f_triangulateio
    REAL(C_REAL),POINTER               :: pointlist(:)                  !* In / out C_REAL ARRAY*!
    REAL(C_REAL),POINTER               :: pointattributelist(:)         !* In / out C_REAL ARRAY*!
    INTEGER(C_INT),POINTER             :: pointmarkerlist(:)            !* In / out INT ARRAY  *!
    INTEGER(C_INT) numberofpoints                                       !* In / out            *!
    INTEGER(C_INT) numberofpointattributes                              !* In / out            *!

    INTEGER(C_INT),POINTER             :: trianglelist(:)               !* In / out INT ARRAY  *!
    REAL(C_REAL),POINTER               :: triangleattributelist(:)      !* In / out C_REAL ARRAY*!
    REAL(C_REAL),POINTER               :: trianglearealist(:)           !* In only C_REAL ARRAY *!
    INTEGER(C_INT),POINTER             :: neighborlist(:)               !* Out only INT ARRAY  *!
    INTEGER(C_INT) numberoftriangles                                    !* In / out            *!
    INTEGER(C_INT) numberofcorners                                      !* In / out            *!
    INTEGER(C_INT) numberoftriangleattributes                           !* In / out            *!

    INTEGER(C_INT),POINTER             :: segmentlist(:)                !* In / out INT ARRAY  *!
    INTEGER(C_INT),POINTER             :: segmentmarkerlist(:)          !* In / out INT ARRAY  *!
    INTEGER(C_INT) numberofsegments                                     !* In / out            *!

    REAL(C_REAL),POINTER              :: holelist(:)                    !* In / pointer to C_REAL array copied out *!
    INTEGER(C_INT) numberofholes                                        !* In / copied out     *!

    REAL(C_REAL),POINTER              :: regionlist(:)                  !* In / pointer to C_REAL array copied out *!
    INTEGER(C_INT) numberofregions                                      !* In / copied out     *!

    INTEGER(C_INT),POINTER             :: edgelist(:)                   !* Out only INT ARRAY  *!
    INTEGER(C_INT),POINTER             :: edgemarkerlist(:)             !* Not used with Voronoi diagram; out only INT ARRAY */
    REAL(C_REAL),POINTER               :: normlist(:)                   !* Used only with Voronoi diagram; out only C_REAL ARRAY*!
    INTEGER(C_INT) numberofedges                                        !* Out only            *!
  END TYPE f_triangulateio

  INTERFACE
    ! Interface to allow Fortran to use the functions defined in triangle.c

    ! triangulate
    ! USAGE: call ctriangulate(triswitches, in, out, vorout)
    ! trifree
    ! USAGE: call ctrifree(memptr)
    ! trimalloc
    ! USAGE: cpointer = ctrimalloc(allocationsize)
    SUBROUTINE ctriangulate(triswitches, in, out, vorout) bind(c, name="triangulate")
      IMPORT :: triangulateio,C_CHAR
      CHARACTER(kind=C_CHAR) triswitches(*)                               !* In char array of flags to use      *!
      TYPE(triangulateio) in                                              !* input shape data structure         *!
      TYPE(triangulateio) out                                             !* output shape data structure        *!
      TYPE(triangulateio) vorout                                          !* output Voroni shape data structure *!
    END SUBROUTINE ctriangulate

    SUBROUTINE ctrifree(memptr) bind(c,name="trifree")
      IMPORT :: C_PTR
      TYPE(C_PTR) memptr
    END SUBROUTINE ctrifree

    TYPE(C_PTR) FUNCTION ctrimalloc(allocationsize) bind(c,name="trimalloc")
      IMPORT :: C_INT, C_PTR
      INTEGER(C_INT) allocationsize
    END FUNCTION ctrimalloc
  END INTERFACE

  CONTAINS
    ! this associates the pointers to arrays in c_shapes to arrays in f_shapes
    !
    ! USAGE: call copytriangles_c_to_f(c_shapes,f_shapes)
    SUBROUTINE copytriangles_c_to_f(c_shapes,f_shapes)
      USE, INTRINSIC:: ISO_C_BINDING, only: c_f_pointer
      USE filehandling, only: stderr
      IMPLICIT NONE
      TYPE(triangulateio),INTENT(IN)    :: c_shapes
      TYPE(f_triangulateio),INTENT(OUT) :: f_shapes
      INTEGER(C_INT)                    :: numberofpoints,numberofpointattributes,numberoftriangles,numberofcorners
      INTEGER(C_INT)                    :: numberoftriangleattributes,numberofsegments,numberofholes,numberofregions
      INTEGER(C_INT)                    :: numberofedges

      numberofpoints             = c_shapes%numberofpoints
      numberofpointattributes    = c_shapes%numberofpointattributes
      numberoftriangles          = c_shapes%numberoftriangles
      numberofcorners            = c_shapes%numberofcorners
      numberoftriangleattributes = c_shapes%numberoftriangleattributes
      numberofsegments           = c_shapes%numberofsegments
      numberofholes              = c_shapes%numberofholes
      numberofregions            = c_shapes%numberofregions
      numberofedges              = c_shapes%numberofedges
      f_shapes%numberofpoints             = numberofpoints
      f_shapes%numberofpointattributes    = numberofpointattributes
      f_shapes%numberoftriangles          = numberoftriangles
      f_shapes%numberofcorners            = numberofcorners
      f_shapes%numberoftriangleattributes = numberoftriangleattributes
      f_shapes%numberofsegments           = numberofsegments
      f_shapes%numberofholes              = numberofholes
      f_shapes%numberofregions            = numberofregions
      f_shapes%numberofedges              = numberofedges

      CALL c_f_pointer(c_shapes%pointlist, f_shapes%pointlist, [2*numberofpoints])
      CALL c_f_pointer(c_shapes%pointmarkerlist, f_shapes%pointmarkerlist, [numberofpoints])
      CALL c_f_pointer(c_shapes%pointattributelist, f_shapes%pointattributelist, [numberofpoints*numberofpointattributes])
      CALL c_f_pointer(c_shapes%trianglelist, f_shapes%trianglelist, [numberoftriangles*numberofcorners])
      CALL c_f_pointer(c_shapes%triangleattributelist, f_shapes%triangleattributelist, &
                                                            [numberoftriangles*numberoftriangleattributes])
      CALL c_f_pointer(c_shapes%trianglearealist, f_shapes%trianglearealist, [numberoftriangles])
      CALL c_f_pointer(c_shapes%neighborlist, f_shapes%neighborlist, [3*numberoftriangles])
      CALL c_f_pointer(c_shapes%segmentlist, f_shapes%segmentlist, [2*numberofsegments])
      CALL c_f_pointer(c_shapes%segmentmarkerlist, f_shapes%segmentmarkerlist, [numberofsegments])
      CALL c_f_pointer(c_shapes%holelist, f_shapes%holelist, [2*numberofholes])
      CALL c_f_pointer(c_shapes%regionlist, f_shapes%regionlist, [4*numberofregions])
      CALL c_f_pointer(c_shapes%edgelist, f_shapes%edgelist, [2*numberofedges])
      CALL c_f_pointer(c_shapes%edgemarkerlist, f_shapes%edgemarkerlist, [numberofedges])
      CALL c_f_pointer(c_shapes%normlist, f_shapes%normlist, [2*numberofedges])
    END SUBROUTINE copytriangles_c_to_f

    ! This takes f_shapes and associates its arrays with the pointers to arrays in c_shapes, allocating those arrays properly.
    ! This should be called before the arrays in f_shapes are allocated, especially before any values are put in those arrays, but
    ! it should only be called after f_shapes has values for its members numberofpoints, numberofpointattributes, numberoftriangles,
    ! numberofcorners, numberoftriangleattributes, numberofsegments, numberofholes, numberofregions, and numberofedges so that properly
    ! sized arrays can be allocated

    ! THIS ALLOCATES THE WHOLE OBJECT, EVEN OUTPUT ONLY ARRAYS. DON'T WASTE MEMORY IF YOU DON'T NEED IT
    !
    ! USAGE: call initializetriangles_f_to_c_all(f_shapes,c_shapes)
    SUBROUTINE initializetriangles_f_to_c_all(f_shapes,c_shapes)
      IMPLICIT NONE
      TYPE(f_triangulateio),INTENT(INOUT) :: f_shapes
      TYPE(triangulateio),INTENT(OUT)     :: c_shapes

      CALL allocate_points_ftoc(f_shapes,c_shapes)
      CALL allocate_triangles(f_shapes,c_shapes)
      ! Neigborlist is only used for output! Don't waste memory if you don't need it !!!!!!!!!!!!!!!!
      CALL allocate_neighborlist(f_shapes,c_shapes)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL allocate_segments(f_shapes,c_shapes)
      CALL allocate_holes(f_shapes,c_shapes)
      CALL allocate_regions(f_shapes,c_shapes)
      ! Voroni edges are only used as outputs! Don't waste memory if you don't need it !!!!!!!!!!!!!!
      CALL allocate_voroni_edges(f_shapes,c_shapes)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END SUBROUTINE initializetriangles_f_to_c_all

    ! This takes f_shapes and associates its arrays with the pointers to arrays in c_shapes, allocating those arrays properly.
    ! This should be called before the arrays in f_shapes are allocated, especially before any values are put in those arrays, but
    ! it should only be called after f_shapes has values for its members numberofpoints, numberofpointattributes, numberoftriangles,
    ! numberofcorners, numberoftriangleattributes, numberofsegments, numberofholes, and numberofregions so that properly
    ! sized arrays can be allocated

    ! THIS ALLOCATES ONLY THE INPUT ARRAYS. The neighborlist and voroni edge arrays will remain unallocated
    !
    ! USAGE: call initializetriangles_f_to_c_inputs(f_shapes,c_shapes)
    SUBROUTINE initializetriangles_f_to_c_inputs(f_shapes,c_shapes)
      IMPLICIT NONE
      TYPE(f_triangulateio),INTENT(INOUT) :: f_shapes
      TYPE(triangulateio),INTENT(OUT)     :: c_shapes

      CALL allocate_points_ftoc(f_shapes,c_shapes)
      CALL allocate_triangles(f_shapes,c_shapes)
      CALL allocate_segments(f_shapes,c_shapes)
      CALL allocate_holes(f_shapes,c_shapes)
      CALL allocate_regions(f_shapes,c_shapes)
    END SUBROUTINE initializetriangles_f_to_c_inputs

    SUBROUTINE allocate_points_ftoc(f_shapes,c_shapes)
      USE, INTRINSIC:: ISO_C_BINDING, only: c_f_pointer
      IMPLICIT NONE
      TYPE(f_triangulateio),INTENT(INOUT) :: f_shapes
      TYPE(triangulateio),INTENT(OUT)     :: c_shapes
      INTEGER(C_INT)                      :: numberofpoints,numberofpointattributes
      INTEGER(C_INT)                      :: dumint
      REAL(C_REAL)                        :: dumreal

      numberofpoints                      = f_shapes%numberofpoints
      numberofpointattributes             = f_shapes%numberofpointattributes
      c_shapes%numberofpoints             = numberofpoints
      c_shapes%numberofpointattributes    = numberofpointattributes

      c_shapes%pointlist = ctrimalloc(INT(2*numberofpoints*SIZEOF(dumreal),C_INT))
      CALL c_f_pointer(c_shapes%pointlist, f_shapes%pointlist, [2*numberofpoints])
      c_shapes%pointmarkerlist = ctrimalloc(INT(numberofpoints*SIZEOF(dumint),C_INT))
      CALL c_f_pointer(c_shapes%pointmarkerlist, f_shapes%pointmarkerlist, [numberofpoints])
      c_shapes%pointattributelist = ctrimalloc(INT(numberofpoints*numberofpointattributes*SIZEOF(dumreal),C_INT))
      CALL c_f_pointer(c_shapes%pointattributelist, f_shapes%pointattributelist, [numberofpoints*numberofpointattributes])
    END SUBROUTINE allocate_points_ftoc

    SUBROUTINE allocate_triangles(f_shapes,c_shapes)
      USE, INTRINSIC:: ISO_C_BINDING, only: c_f_pointer
      IMPLICIT NONE
      TYPE(f_triangulateio),INTENT(INOUT) :: f_shapes
      TYPE(triangulateio),INTENT(OUT)     :: c_shapes
      INTEGER(C_INT)                      :: numberoftriangles,numberofcorners,numberoftriangleattributes
      INTEGER(C_INT)                      :: dumint
      REAL(C_REAL)                        :: dumreal

      numberoftriangles                   = f_shapes%numberoftriangles
      numberofcorners                     = f_shapes%numberofcorners
      numberoftriangleattributes          = f_shapes%numberoftriangleattributes
      c_shapes%numberoftriangles          = numberoftriangles
      c_shapes%numberofcorners            = numberofcorners
      c_shapes%numberoftriangleattributes = numberoftriangleattributes

      c_shapes%trianglelist = ctrimalloc(INT(numberoftriangles*numberofcorners*SIZEOF(dumint),C_INT))
      CALL c_f_pointer(c_shapes%trianglelist, f_shapes%trianglelist, [numberoftriangles*numberofcorners])
      c_shapes%triangleattributelist = ctrimalloc(INT(numberoftriangles*numberoftriangleattributes*SIZEOF(dumreal),C_INT))
      CALL c_f_pointer(c_shapes%triangleattributelist, f_shapes%triangleattributelist, &
                                                                  [numberoftriangles*numberoftriangleattributes])
      c_shapes%trianglearealist = ctrimalloc(INT(numberoftriangles*SIZEOF(dumreal),C_INT))
      CALL c_f_pointer(c_shapes%trianglearealist, f_shapes%trianglearealist, [numberoftriangles])
    END SUBROUTINE allocate_triangles

    SUBROUTINE allocate_segments(f_shapes,c_shapes)
      USE, INTRINSIC:: ISO_C_BINDING, only: c_f_pointer
      IMPLICIT NONE
      TYPE(f_triangulateio),INTENT(INOUT) :: f_shapes
      TYPE(triangulateio),INTENT(OUT)     :: c_shapes
      INTEGER(C_INT)                      :: numberofsegments
      INTEGER(C_INT)                      :: dumint

      numberofsegments                    = f_shapes%numberofsegments
      c_shapes%numberofsegments           = numberofsegments

      c_shapes%segmentlist = ctrimalloc(INT(2*numberofsegments*SIZEOF(dumint),C_INT))
      CALL c_f_pointer(c_shapes%segmentlist, f_shapes%segmentlist, [2*numberofsegments])
      c_shapes%segmentmarkerlist = ctrimalloc(INT(numberofsegments*SIZEOF(dumint),C_INT))
      CALL c_f_pointer(c_shapes%segmentmarkerlist, f_shapes%segmentmarkerlist, [numberofsegments])
    END SUBROUTINE allocate_segments

    SUBROUTINE allocate_holes(f_shapes,c_shapes)
      USE, INTRINSIC:: ISO_C_BINDING, only: c_f_pointer
      IMPLICIT NONE
      TYPE(f_triangulateio),INTENT(INOUT) :: f_shapes
      TYPE(triangulateio),INTENT(OUT)     :: c_shapes
      INTEGER(C_INT)                      :: numberofholes
      REAL(C_REAL)                        :: dumreal

      numberofholes                       = f_shapes%numberofholes
      c_shapes%numberofholes              = numberofholes

      c_shapes%holelist = ctrimalloc(INT(2*numberofholes*SIZEOF(dumreal),C_INT))
      CALL c_f_pointer(c_shapes%holelist, f_shapes%holelist, [2*numberofholes])
    END SUBROUTINE allocate_holes

    SUBROUTINE allocate_regions(f_shapes,c_shapes)
      USE, INTRINSIC:: ISO_C_BINDING, only: c_f_pointer
      IMPLICIT NONE
      TYPE(f_triangulateio),INTENT(INOUT) :: f_shapes
      TYPE(triangulateio),INTENT(OUT)     :: c_shapes
      INTEGER(C_INT)                      :: numberofregions
      REAL(C_REAL)                        :: dumreal

      numberofregions                     = f_shapes%numberofregions
      c_shapes%numberofregions            = numberofregions

      c_shapes%regionlist = ctrimalloc(INT(4*numberofregions*SIZEOF(dumreal),C_INT))
      CALL c_f_pointer(c_shapes%regionlist, f_shapes%regionlist, [4*numberofregions])
    END SUBROUTINE allocate_regions

    SUBROUTINE allocate_voroni_edges(f_shapes,c_shapes)
      USE, INTRINSIC:: ISO_C_BINDING, only: c_f_pointer
      IMPLICIT NONE
      TYPE(f_triangulateio),INTENT(INOUT) :: f_shapes
      TYPE(triangulateio),INTENT(OUT)     :: c_shapes
      INTEGER(C_INT)                      :: numberofedges
      INTEGER(C_INT)                      :: dumint
      REAL(C_REAL)                        :: dumreal

      numberofedges                       = f_shapes%numberofedges
      c_shapes%numberofedges              = numberofedges

      c_shapes%edgelist = ctrimalloc(INT(2*numberofedges*SIZEOF(dumint),C_INT))
      CALL c_f_pointer(c_shapes%edgelist, f_shapes%edgelist, [2*numberofedges])
      c_shapes%edgemarkerlist = ctrimalloc(INT(numberofedges*SIZEOF(dumint),C_INT))
      CALL c_f_pointer(c_shapes%edgemarkerlist, f_shapes%edgemarkerlist, [numberofedges])
      c_shapes%normlist = ctrimalloc(INT(2*numberofedges*SIZEOF(dumreal),C_INT))
      CALL c_f_pointer(c_shapes%normlist, f_shapes%normlist, [2*numberofedges])
    END SUBROUTINE allocate_voroni_edges

    SUBROUTINE allocate_neighborlist(f_shapes,c_shapes)
      USE, INTRINSIC:: ISO_C_BINDING, only: c_f_pointer
      IMPLICIT NONE
      TYPE(f_triangulateio),INTENT(INOUT) :: f_shapes
      TYPE(triangulateio),INTENT(OUT)     :: c_shapes
      INTEGER(C_INT)                      :: numberoftriangles
      INTEGER(C_INT)                      :: dumint

      numberoftriangles                   = f_shapes%numberoftriangles

      c_shapes%neighborlist = ctrimalloc(INT(3*numberoftriangles*SIZEOF(dumint),C_INT))
      CALL c_f_pointer(c_shapes%neighborlist, f_shapes%neighborlist, [3*numberoftriangles])
    END SUBROUTINE allocate_neighborlist
    
    SUBROUTINE init_outputs(c_shape)
      USE, INTRINSIC:: ISO_C_BINDING, only: C_NULL_PTR
      IMPLICIT NONE
      TYPE(triangulateio),INTENT(INOUT)   :: c_shape
      
      c_shape%pointlist = C_NULL_PTR
      c_shape%pointattributelist = C_NULL_PTR
      c_shape%pointmarkerlist = C_NULL_PTR
      c_shape%numberofpoints = 0
      c_shape%numberofpointattributes = 0

      c_shape%trianglelist = C_NULL_PTR
      c_shape%triangleattributelist = C_NULL_PTR
      c_shape%trianglearealist = C_NULL_PTR
      c_shape%neighborlist = C_NULL_PTR
      c_shape%numberoftriangles = 0
      c_shape%numberofcorners = 0
      c_shape%numberoftriangleattributes = 0

      c_shape%segmentlist  = C_NULL_PTR
      c_shape%segmentmarkerlist = C_NULL_PTR
      c_shape%numberofsegments = 0

      c_shape%holelist = C_NULL_PTR
      c_shape%numberofholes = 0

      c_shape%regionlist = C_NULL_PTR
      c_shape%numberofregions = 0

      c_shape%edgelist = C_NULL_PTR
      c_shape%edgemarkerlist = C_NULL_PTR
      c_shape%normlist = C_NULL_PTR
      c_shape%numberofedges = 0
    END SUBROUTINE init_outputs
    
    SUBROUTINE deallocate_triangulateio(f_shape,input,neigh,voroni)
      IMPLICIT NONE
      
      TYPE(f_triangulateio),INTENT(INOUT)   :: f_shape
      LOGICAL                               :: input,neigh,voroni
      
      if(f_shape%numberofpoints .NE. 0)deallocate(f_shape%pointlist)
      if(f_shape%numberofpoints*f_shape%numberofpointattributes .NE. 0)deallocate(f_shape%pointattributelist)
      if(f_shape%numberofpoints .NE. 0)deallocate(f_shape%pointmarkerlist)
      f_shape%numberofpoints = 0
      f_shape%numberofpointattributes = 0

      if(f_shape%numberoftriangles .NE. 0)deallocate(f_shape%trianglelist)
      if(f_shape%numberoftriangles*f_shape%numberoftriangleattributes .NE. 0)deallocate(f_shape%triangleattributelist)
      if((f_shape%numberoftriangles .NE. 0) .AND. input)deallocate(f_shape%trianglearealist)
      if((f_shape%numberoftriangles .NE. 0) .AND. (.NOT. input .AND. neigh))deallocate(f_shape%neighborlist)
      f_shape%numberoftriangles = 0
      f_shape%numberofcorners = 0
      f_shape%numberoftriangleattributes = 0

      if(f_shape%numberofsegments .NE. 0)deallocate(f_shape%segmentlist)
      if(f_shape%numberofsegments .NE. 0)deallocate(f_shape%segmentmarkerlist)
      f_shape%numberofsegments = 0

      if((f_shape%numberofholes .NE. 0) .AND. input)deallocate(f_shape%holelist)
      f_shape%numberofholes = 0

      if((f_shape%numberofregions .NE. 0) .AND. input)deallocate(f_shape%regionlist)
      f_shape%numberofregions = 0

      if(f_shape%numberofedges .NE. 0)deallocate(f_shape%edgelist)
      if((f_shape%numberofedges .NE. 0) .AND. .NOT. voroni)deallocate(f_shape%edgemarkerlist)
      if((f_shape%numberofedges .NE. 0) .AND. voroni)deallocate(f_shape%normlist)
      f_shape%numberofedges = 0
    END SUBROUTINE deallocate_triangulateio
END MODULE triangle_c_wrap
