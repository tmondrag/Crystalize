MODULE basictypes
  USE, INTRINSIC :: ISO_C_BINDING, only: C_INT,C_BOOL
  USE            :: triangle_c_wrap, only: C_REAL
  IMPLICIT NONE

  ! data types
  ! lattice node - measued and quantified points from the input plus some derived info
  TYPE lattice_node
    INTEGER(C_INT)                              :: n_index        ! self-referential index for convenience
    REAL(C_REAL),DIMENSION(1:2)                 :: location       ! x-y coordinates of node
    INTEGER(C_INT)                              :: q_state        ! discretized state (spin state, orthogonal orientation state)
    REAL(C_REAL),DIMENSION(1:3)                 :: n_orientation  ! non-discretized orientation state (orientation of neighbor edges)
    INTEGER(C_INT)                              :: major_axis     ! Index into n_orientation of major axis
  END TYPE
  ! lattice - container for raw and simplistic data. Nodes, cells, edges, & attached info
  TYPE lattice
  END TYPE
END MODULE basictypes
