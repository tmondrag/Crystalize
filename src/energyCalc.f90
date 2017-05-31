MODULE energyCalc
CONTAINS
  FUNCTION calculate_edge_energy_positional(enrgScale,restBose,restFermi,q1,q2,length) RESULT(edgeEnergy)
    USE triangle_c_wrap, only: C_REAL
    USE, INTRINSIC :: ISO_C_BINDING, only: C_INT
    IMPLICIT NONE

    REAL(C_REAL)              :: enrgScale,restBose,restFermi,length,edgeEnergy
    INTEGER(C_INT)            :: q1, q2

    IF (q1 == q2) THEN
      edgeEnergy = enrgScale*(length-restBose)**2/length**2
    ELSE
      edgeEnergy = enrgScale*(length-restFermi)**2/length**2
    END IF
  END FUNCTION

  FUNCTION calculate_edge_energy_qstate(enrgScale,q1,q2,length) RESULT(edgeEnergy)
    USE triangle_c_wrap, only: C_REAL
    USE, INTRINSIC :: ISO_C_BINDING, only: C_INT
    IMPLICIT NONE

    REAL(C_REAL)              :: enrgScale,length,edgeEnergy
    INTEGER(C_INT)            :: q1, q2

    IF (q1 == q2) THEN
      edgeEnergy = enrgScale*0.1_C_REAL/length**2
    ELSE
      edgeEnergy = enrgScale*5.0_C_REAL/length**2
    END IF
  END FUNCTION
END MODULE energyCalc
