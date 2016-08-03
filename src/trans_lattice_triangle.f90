MODULE trans_lattice_triangle
  IMPLICIT NONE

CONTAINS

  SUBROUTINE transfer_lattice_to_triangle(which_lattice,f_shape,c_shape)
    USE basictypes, only: lattice
    USE triangle_c_wrap, only: triangulateio,f_triangulateio,allocate_points_ftoc,allocate_segments,allocate_regions,allocate_holes,C_REAL
    IMPLICIT NONE
    TYPE(lattice),INTENT(IN)                          :: which_lattice
    TYPE(triangulateio),INTENT(OUT)                   :: c_shape
    TYPE(f_triangulateio),INTENT(OUT)                 :: f_shape
    INTEGER                                           :: i

    f_shape%numberofpoints = which_lattice%num_nodes
    f_shape%numberofpointattributes = 4               ! q_state, major_axis, enumer, enrg_node
    f_shape%numberofsegments = which_lattice%num_edges
    f_shape%numberoftriangles = which_lattice%num_cells
    f_shape%numberofcorners = 3
    f_shape%numberoftriangleattributes = 5            ! q_state, major_axis, enumer, num_broken_edges, enrg_cell
    f_shape%numberofholes = 0
    f_shape%numberofregions = which_lattice%num_grains

    CALL allocate_points_ftoc(f_shape,c_shape)
    CALL allocate_segments(f_shape,c_shape)
    CALL allocate_triangles(f_shape,c_shape)
    CALL allocate_regions(f_shape,c_shape)
    CALL allocate_holes(f_shape,c_shape)

    DO i=1,which_lattice%num_nodes
      f_shape%pointlist(2*i-1) = which_lattice%nodes(i)%location(1)
      f_shape%pointlist(2*i) = which_lattice%nodes(i)%location(2)
      f_shape%pointAttributelist(4*i-3) = REAL(which_lattice%nodes(i)%q_state,C_REAL)
      f_shape%pointAttributelist(4*i-2) = which_lattice%nodes(i)%major_axis
      f_shape%pointAttributelist(4*i-1) = REAL(which_lattice%nodes(i)%enumer,C_REAL)
      f_shape%pointAttributelist(4*i) = which_lattice%nodes(i)%enrg_node
    END DO
    f_shapes%pointmarkerlist = 0

    DO i=1,which_lattice%num_edges
      f_shape%segmentlist(2*i-1) = which_lattice%edges(i)%neigh_nodes(1)
      f_shape%segmentlist(2*i) = which_lattice%edges(i)%neigh_nodes(2)
    END DO
    f_shapes%segmentmarkerlist = 1

    DO i=1,which_lattice%num_cells
      f_shape%trianglelist(3*i-2) = which_lattice%cells(i)%neigh_nodes(1)
      f_shape%trianglelist(3*i-1) = which_lattice%cells(i)%neigh_nodes(2)
      f_shape%trianglelist(3*i) = which_lattice%cells(i)%neigh_nodes(3)
      f_shape%triangleAttributelist(5*i-4) = REAL(which_lattice%cells(i)%q_state,C_REAL)
      f_shape%triangleAttributelist(5*i-3) = which_lattice%cells(i)%major_axis
      f_shape%triangleAttributelist(5*i-2) = REAL(which_lattice%cells(i)%enumer,C_REAL)
      f_shape%triangleAttributelist(5*i-1) = REAL(which_lattice%cells(i)%num_broken_edges,C_REAL)
      f_shape%triangleAttributelist(5*i) = which_lattice%cells(i)%enrg_node
    END DO
    f_shapes%trianglearealist = 0.0
  END SUBROUTINE transfer_lattice_to_triangle
END MODULE trans_lattice_triangle
