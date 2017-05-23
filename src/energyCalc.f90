MODULE energyCalc
CONTAINS
  FUNCTION calculate_edge_energy(enrgScale,restBose,restFermi,q1,q2,length) RESULT(edgeEnergy)
    IF (q1 == q2) THEN
      edgeEnergy = enrgScale*(length-restBose)**2
    ELSE
      edgeEnergy = enrgScale*(length-restFermi)**2
    END IF
    edgeEnergy = edgeEnergy + enrgScale/length**2
  END FUNCTION
END MODULE energyCalc
