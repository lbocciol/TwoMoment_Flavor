MODULE IntegrationModule 

  USE KindModule, ONLY: &
    DP, Zero, Half, One
  USE UnitsModule, ONLY: &
    MeV, Centimeter, Second
  USE ProgramHeaderModule, ONLY: &
    nDOFX, nDOFE, &
    nNodesX, nNodesE, & 
    iZ_B0, iZ_E0, iE_B0, iE_E0, &
    iX_B0
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeNeutrinoOpacities_EC_Points
  USE RadiationFieldsModule, ONLY: &
    nSpecies, uPR, iPR_D, nPR, &
    iNuEE, iNuEE_Bar, &
    iNuEX, iNuEX_Bar, &
    iNuE, iNuE_Bar
  USE MeshModule, ONLY: &
    MeshE, NodeCoordinate
  USE UtilitiesModule, ONLY: &
    NodeNumber
  USE FluidFieldsModule, ONLY: &                          
    uPF, uAF, &
    nPF, nAF, iPF_D, iAF_T, iAF_Ye                             
  USE FlavorOpacitiesModule, ONLY: & 
    ComputeEmission, &
    ComputeEmissionAbsorptionAvg, &
    ComputeNu4CollisionTerms
  USE InitializationModule, ONLY: &
    AA, BB5, BB6, nRK, nRKOrder

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: TimeStepping_Decoherence

  INTEGER  :: nE_G, nX_G
  REAL(DP), ALLOCATABLE :: D(:), T(:), Ye(:)
  REAL(DP), ALLOCATABLE :: Chi(:,:,:), E_N(:)
  REAL(DP), ALLOCATABLE :: eta(:,:,:), kappa(:,:,:)
  REAL(DP), ALLOCATABLE :: kappa_avg(:,:,:), eta_avg(:,:,:)
  REAL(DP), ALLOCATABLE :: PF_N(:,:) , AF_N(:,:) 
  REAL(DP), ALLOCATABLE :: uPR_N(:,:,:,:)
  REAL(DP), ALLOCATABLE :: J0(:,:,:), J0_old(:,:,:)
  REAL(DP), ALLOCATABLE :: J0_check(:,:,:)

  REAL(DP), ALLOCATABLE :: CollisionTermNu4(:,:,:)

CONTAINS

  SUBROUTINE TimeStepping_Decoherence(dt, nFlavors)
    
    INTEGER,  INTENT(IN)    :: nFlavors
    REAL(DP), INTENT(INOUT) :: dt
    
    LOGICAL  :: Recycle
    INTEGER  :: iS, iEk, iE, iX
    INTEGER  :: i1,i2,i3
    INTEGER  :: MaxIter = 100
    REAL(DP) :: Error

    CALL InitializeArrays(nFlavors)
    
    !$OMP PARALLEL DO PRIVATE(iS)
    DO iS = 1, nFlavors

      CALL ComputeNeutrinoOpacities_EC_Points &   
             ( 1, nE_G, 1, nX_G, &                
               E_N (:), &                         
               PF_N(:,iPF_D ), &                  
               AF_N(:,iAF_T ), &                  
               AF_N(:,iAF_Ye), &                  
               iS, Chi(:,:,iS) )                  

    END DO
    !$OMP END PARALLEL DO

    CALL ComputeEmission(  &
        E_N(:), Chi(:,:,:), &
        kappa(:,:,:), eta(:,:,:), &
        nFlavors, nE_G, nX_G )
      
    CALL ComputeEmissionAbsorptionAvg(   &
       kappa_avg, eta_avg, kappa, eta, &
       nFlavors, nE_G, nX_G )

    J0_old = uPR_N(:,:,iPR_D,:)

    !$OMP PARALLEL DO PRIVATE(iEk,iX)
    DO iEk = 1,nE_G   
      DO iX = 1,nX_G
        
        CALL ComputeNu4CollisionTerms( &
            J0_old(:,iX,:), &
            CollisionTermNu4(iEk,iX,:), &
            E_N, nE_G, iEk, nFlavors)
      
      END DO
    END DO
    !$OMP END PARALLEL DO

    !CALL Rk_Explicit_Step(dt, nFlavors)
    
    Recycle = .TRUE.
    DO WHILE ( Recycle .EQV. .TRUE. )
      
      CALL SolveThisIteration(dt,J0)
      
      CALL SolveThisIteration(0.8d0 * dt,J0_check)
      
      Error = 0
      !$OMP PARALLEL DO 
      DO iE = 1,nE_G
      DO iX = 1,nX_G
      DO iS = 1,nSpecies
        IF (J0_check(iE,iX,iS) .GT. 1.0d-10 ) &
        Error = MAX(Error, &
            ABS( ( J0_check(iE,iX,iS) - J0(iE,iX,iS) ) / &
                     J0_check(iE,iX,iS) ) )
      END DO
      END DO
      END DO
      !$OMP END PARALLEL DO
      
      IF ( Error .GT. 1.0d-3 ) THEN
        
        dt = 0.8d0 * dt
        Recycle = .TRUE.

      ELSE

        dt = 1.2d0 * dt
        Recycle = .FALSE.
      
      END IF

    END DO

    CALL FinalizeArrays

  END SUBROUTINE TimeStepping_Decoherence

  SUBROUTINE Rk_Explicit_Step( dt, nFlavors)

    INTEGER,  INTENT(IN)    :: nFlavors
    REAL(DP), INTENT(INOUT) :: dt
    
    REAL(DP) :: J0_temp(nE_G,nX_G,nSpecies)
    REAL(DP) :: dJ0dt  (nE_G,nX_G,nSpecies,nRK)
    REAL(DP) :: MaxError, Error
    REAL(DP) :: Accuracy = 1d-3
    REAL(DP) :: Increase = 1.1d0
    REAL(DP) :: dJ0

    LOGICAL :: Recycle
    INTEGER :: k, l, iE, iX, iS, iEk

    Recycle = .TRUE.

    DO WHILE( Recycle )

      DO k=1,nRK
  
        J0_temp = J0_old
  
        DO l = 1, k-1
  
          J0_temp(:,:,:) = J0_temp(:,:,:) + dJ0dt(:,:,:,l) * AA(k,l) * dt
  
        END DO
      
        !$OMP PARALLEL DO PRIVATE(iEk,iX)
        DO iEk = 1,nE_G
          DO iX = 1,nX_G
  
            CALL ComputeNu4CollisionTerms( &
                J0_old(:,nX_G,:), &
                CollisionTermNu4(iEk,iX,:), &
                E_N, nE_G, iEk, nFlavors)
            
            dJ0dt(iEk,iX,:,k) = CollisionTermNu4(iEk,iX,:)
          
          END DO
        END DO
        !$OMP END PARALLEL DO
        
      END DO
      
      MaxError = Zero
      DO iE = 1,nE_G
      DO iX = 1,nX_G
      DO iS = 1,nSpecies
        
        dJ0 = Zero
        
        DO k = 1,NRK
          
          dJ0 = dJ0 + dJ0dt(iE,iX,iS,k) * BB5(k) * dt
          Error = Error + dJ0dt(iE,iX,iS,k) * (BB5(k) - BB6(k)) * dt
        
        END DO
      
        J0(iE,iX,iS) = J0_old(iE,iX,iS) + dJ0
  
        MaxError = MAX( MaxError, ABS(Error)/ABS(J0_old(iE,iX,iS)) )
          
       END DO
       END DO
       END DO

      IF (MaxError > Accuracy) THEN

	dt = dt * 0.9 * (Accuracy/Maxerror)**(One/(NRKOrder-One))
	Recycle = .TRUE.

      ELSE
	
        dt = dt * Increase
	Recycle = .FALSE.

        IF( MaxError > 0 ) &
            dt = dt * MIN( One, ( (Accuracy/MaxError)** &
                 (One/MAX(One,NRKOrder)) ) / Increase) ;
      
       END IF

    END DO !End of DO WHILE

  END SUBROUTINE Rk_Explicit_Step


  SUBROUTINE SolveThisIteration(dt,J0_new)
    
    REAL(DP), INTENT(IN)  :: dt
    REAL(DP), INTENT(OUT) :: J0_new(nE_G,nX_G,nSpecies)

    INTEGER :: iS, iE, iX

    !$OMP PARALLEL DO PRIVATE( iE, iX, iS )
    DO iE = 1,nE_G
      DO iX = 1,nX_G
        DO iS = 1,nSpecies
        
          IF ( iS .EQ. iNuEE ) THEN

            J0_new(iE,iX,iS) = ( J0_old(iE,iX,iS) + dt * &
                    ( eta(iE,iX,iNuE) + CollisionTermNu4(iE,iX,iS) ) ) / &
                ( One + ( eta_avg(iE,iX,iS) + kappa_avg(iE,iX,iS) ) * dt )


          ELSE IF ( iS .EQ. iNuEE_Bar ) THEN

            J0_new(iE,iX,iS) = ( J0_old(iE,iX,iS) + dt * &
                    ( eta(iE,iX,iNuE_Bar) + CollisionTermNu4(iE,iX,iS) ) ) / &
                ( One + ( eta_avg(iE,iX,iS) + kappa_avg(iE,iX,iS) ) * dt )


          ELSE IF ( (iS .EQ. iNuEX) .OR. (iS .EQ. iNuEX_Bar) ) THEN

            J0_new(iE,iX,iS) = ( J0_old(iE,iX,iS) + dt * &
                    CollisionTermNu4(iE,iX,iS) ) / &
                ( One + ( eta_avg(iE,iX,iS) + kappa_avg(iE,iX,iS) ) * dt )

          ELSE

            J0_new(iE,iX,iS) = J0_old(iE,iX,iS) + dt * CollisionTermNu4(iE,iX,iS)

          ENDIF

        END DO
      END DO
    END DO
    !$OMP END PARALLEL DO


  END SUBROUTINE SolveThisIteration

  SUBROUTINE InitializeArrays(nFlavors)

    INTEGER, INTENT(IN) :: nFlavors

    INTEGER :: nZ(4), nX(3)
    INTEGER :: iX1, iX2, iX3, iNodeX
    INTEGER :: iE, iNodeE
    INTEGER :: iN_E, iN_X, iPF, iAF, iPR

    INTEGER :: iNodeX_G, iNodeE_G
    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNode

    nZ = iZ_E0 - iZ_B0 + 1           
    nX = nZ(2:4)                     
    nX_G = nDOFX * PRODUCT( nX )     
    nE_G = nNodesE * nZ(1)        

    ALLOCATE( E_N  (nE_G    ) )
    ALLOCATE( PF_N (nX_G,nPF) )
    ALLOCATE( AF_N (nX_G,nAF) )
    ALLOCATE( Chi  (nE_G,nX_G,nFlavors) )
    ALLOCATE( eta  (nE_G,nX_G,nFlavors) )
    ALLOCATE( kappa(nE_G,nX_G,nFlavors) )
 
    ALLOCATE( eta_avg(nE_G,nX_G,nSpecies) )
    ALLOCATE( kappa_avg(nE_G,nX_G,nSpecies) )

    ALLOCATE( uPR_N     (nE_G,nX_G,nPR,nSpecies) )
    ALLOCATE( J0      (nE_G,nX_G,    nSpecies) )
    ALLOCATE( J0_old  (nE_G,nX_G,    nSpecies) )
    ALLOCATE( J0_check(nE_G,nX_G,    nSpecies) )
    
    ALLOCATE( CollisionTermNu4(nE_G,nX_G,nSpecies) )

    !$OMP PARALLEL DO PRIVATE( iE, iNodeE )
    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nZ(1) ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      E_N(iN_E) = NodeCoordinate( MeshE, iE, iNodeE )
    
    END DO
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE( iX1, iX2, iX3, iNodeX )
    DO iN_X = 1, nX_G

      DO iPF = 1, nPF
        iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
        iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
        iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
        iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1    
                                                                           
        PF_N(iN_X,iPF) = uPF(iNodeX,iX1,iX2,iX3,iPF)                       
      END DO                                                               
    
      DO iAF = 1, nAF
        iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
        iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
        iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
        iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

        AF_N(iN_X,iAF) = uAF(iNodeX,iX1,iX2,iX3,iAF)    
      END DO

    END DO
    !$OMP END PARALLEL DO

    iNodeX_G = 0
    DO iX3 = 1, nX(3)
    DO iX2 = 1, nX(2)
    DO iX1 = 1, nX(1)
      DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
      DO iNodeX1 = 1, nNodesX(1)
 
        iNodeX_G = iNodeX_G + 1
    
        iNodeE_G = 0
        DO iE = iE_B0, iE_E0
          DO iNodeE = 1, nNodesE
    
            iNodeE_G = iNodeE_G + 1
    
            iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )
    
            DO iPR = 1,nPR
    
              uPR_N(iNodeE_G,iNodeX_G,iPR,:) &
                 = uPR(iNode,iE,iX1,iX2,iX3,iPR,:)
    
            END DO
    
          END DO
        END DO
  
      END DO
      END DO
      END DO
    END DO
    END DO
    END DO
  

  END SUBROUTINE InitializeArrays

  SUBROUTINE FinalizeArrays
   
    INTEGER :: iN_X, iN_E
    INTEGER :: iX1, iX2, iX3, iNodeX
    INTEGER :: nZ(4), nX(3)
    
    nZ = iZ_E0 - iZ_B0 + 1
    nX = nZ(2:4)
    nZ = iZ_E0 - iZ_B0 + 1
    
    !$OMP PARALLEL DO PRIVATE( iX1, iX2, iX3, iNodeX )
    DO iN_X = 1,nX_G
      DO iN_E = 1,nE_G

        iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
        iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
        iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1)
        iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1

        uPR(iNodeX,iN_E,iX1,iX2,iX3,iPR_D,:) = J0(iN_E,iN_X,:)

      END DO
    END DO
    !$OMP END PARALLEL DO
      
    DEALLOCATE( E_N, PF_N, AF_N, uPR_N, &
                Chi, eta, kappa,      &
                eta_avg, kappa_avg,   &
                J0, J0_old, J0_check, &
                CollisionTermNu4 )

  END SUBROUTINE

END MODULE IntegrationModule 
