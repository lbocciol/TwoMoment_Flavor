MODULE FlavorOpacitiesModule
    
  USE KindModule, ONLY: &
    DP, Zero, Half, One, &
    Two, Three, Four, Five, &
    Pi, TwoPi, FourPi
  USE UnitsModule, ONLY: &
    BoltzmannConstant, &
    MeV, Second
  USE PhysicalConstantsModule, ONLY: &
    SpeedOfLightCGS
  USE ProgramHeaderModule, ONLY: &
    nNodesE, nDOFE, iE_E0, iE_B0
  USE ProgramHeaderModule, ONLY: &
    nDOF, nDOFX, iX_E0, iX_B0, &
    iE_B0
  USE MeshModule, ONLY: &
    MeshE, NodeCoordinate
  USE FluidFieldsModule, ONLY: &
    uAF, iAF_T, &
    iAF_Me, iAF_Mp, iAF_Mn
  USE RadiationFieldsModule, ONLY: &
    nSpecies, iNuE, iNuE_Bar, &
    iNuEE, iNuXX, iNuEX, &
    iNuEE_Bar, iNuXX_Bar, iNuEX_Bar
  USE ReadProfileModule, ONLY: &
    LinearInterpolation1D, &
    QuadraticInterpolation1D

  PRIVATE
  
  PUBLIC :: ComputeEmission
  PUBLIC :: ComputeEmissionAbsorptionAvg
  PUBLIC :: ComputeNu4CollisionTerms
  PUBLIC :: ComputeNu4Kernels

  REAL(DP), ALLOCATABLE :: KernelNu4Pair(:,:,:), KernelNu4Scat(:,:,:)

CONTAINS

  SUBROUTINE ComputeEmission(E_N, Chi, kappa, eta, &
                             nFlavors, nE_G, nX_G)

    INTEGER,  INTENT(IN)  :: nE_G, nX_G, nFlavors
    REAL(DP), INTENT(IN)  :: E_N  (nE_G)
    REAL(DP), INTENT(IN)  :: Chi  (nE_G,nX_G,nFlavors)
    REAL(DP), INTENT(OUT) :: eta  (nE_G,nX_G,nFlavors)
    REAL(DP), INTENT(OUT) :: kappa(nE_G,nX_G,nFlavors)

    INTEGER  :: iS, iN_E, iN_X
    INTEGER  :: iX1, iX2, iX3, iNodeX
    INTEGER  :: nX(3)
    REAL(DP) :: Mnu, kT, FD

    nX = iX_E0 - iX_B0 + 1

    DO iN_X = 1, nX_G
  
      iX3    = MOD( (iN_X-1) / ( nDOFX * nX(1) * nX(2) ), nX(3) ) + iX_B0(3)
      iX2    = MOD( (iN_X-1) / ( nDOFX * nX(1)         ), nX(2) ) + iX_B0(2)
      iX1    = MOD( (iN_X-1) / ( nDOFX                 ), nX(1) ) + iX_B0(1) 
      iNodeX = MOD( (iN_X-1)                            , nDOFX ) + 1    
    
      DO iS = 1, nFlavors
      
        IF ( iS .EQ. iNuE ) THEN
  
          Mnu = + ( uAF(iNodeX,iX1,iX2,iX3,iAF_Me) &
                      + uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
                      - uAF(iNodeX,iX1,iX2,iX3,iAF_Mn) )
          
        ELSE IF ( iS .EQ. iNuE_Bar ) THEN
            
          Mnu = - ( uAF(iNodeX,iX1,iX2,iX3,iAF_Me) &
                      + uAF(iNodeX,iX1,iX2,iX3,iAF_Mp) &
                      - uAF(iNodeX,iX1,iX2,iX3,iAF_Mn) )
        ELSE
  
          Mnu = Zero
  
        END IF
          
        kT = BoltzmannConstant &
             * uAF(iNodeX,iX1,iX2,iX3,iAF_T)
  
        DO iN_E = 1, nE_G
  
          FD = MAX( One / ( EXP( (E_N(iN_E)-Mnu)/kT ) + One ), 1.0d-99 )
          
          eta(iN_E,iN_X,iS) = Chi(iN_E,iN_X,iS) * FD 
          kappa(iN_E,iN_X,iS) = Chi(iN_E,iN_X,iS) - eta(iN_E,iN_X,iS)
          
        END DO

      END DO
    
    END DO

  END SUBROUTINE ComputeEmission

  SUBROUTINE ComputeEmissionABSorptionAvg ( &
        kappa_avg, eta_avg, kappa, eta,   &
        nFlavors, nE_G, nX_G)
    
    INTEGER,  INTENT(IN)  :: nFlavors
    INTEGER,  INTENT(IN)  :: nE_G, nX_G
    REAL(DP), INTENT(IN)  :: kappa  (nE_G,nX_G,nFlavors)
    REAL(DP), INTENT(IN)  :: eta    (nE_G,nX_G,nFlavors)
    REAL(DP), INTENT(OUT) :: kappa_avg(nE_G,nX_G,nSpecies)
    REAL(DP), INTENT(OUT) :: eta_avg(nE_G,nX_G,nSpecies)

    INTEGER :: iS

    DO iS = 1,nSpecies
      
      IF ( iS .EQ. iNuEE ) THEN

        kappa_avg(:,:,iS) = kappa(:,:,iNuE)
        eta_avg(:,:,iS) = eta(:,:,iNuE)

      ELSE IF ( iS .EQ. iNuEE_Bar ) THEN

        kappa_avg(:,:,iS) = kappa(:,:,iNuE_Bar)
        eta_avg(:,:,iS) = eta(:,:,iNuE_Bar)

      ELSE IF ( iS .EQ. iNuEX ) THEN

        kappa_avg(:,:,iS) = Half * kappa(:,:,iNuE)
        eta_avg(:,:,iS) = Half * eta(:,:,iNuE)

      ELSE IF ( iS .EQ. iNuEX_Bar ) THEN

        kappa_avg(:,:,iS) = Half * kappa(:,:,iNuE_Bar)
        eta_avg(:,:,iS) = Half * eta(:,:,iNuE_Bar)

      ELSE

        kappa_avg(:,:,iS) = Zero
        eta_avg(:,:,iS) = Zero

      END IF

    END DO

  END SUBROUTINE ComputeEmissionAbsorptionAvg


  SUBROUTINE ComputeNu4CollisionTerms( J0, CollisionTerm, &
          E_N, nE_G, iEk, nFlavors)
    
    INTEGER,  INTENT(IN)  :: nE_G, iEk
    REAL(DP), INTENT(IN)  :: E_N(nE_G)
    REAL(DP), INTENT(IN)  :: J0           (nE_G,nSpecies)
    REAL(DP), INTENT(OUT) :: CollisionTerm(nSpecies)
    
    INTEGER  :: i1,i2,i3
    INTEGER  :: imat, iE
    REAL(DP) :: J0_E(2,nE_G,nFlavors,nFlavors)

    imat = nSpecies / 2

    ! 1 is for matter, 2 is for antimatter
    DO iE = 1,nE_G

      J0_E(1,iE,:,:) = PackIntoMatrix( J0(iE,1:imat),  nSpecies, nFlavors )
      J0_E(2,iE,:,:) = PackIntoMatrix( J0(iE,imat+1:), nSpecies, nFlavors )

    END DO

    CALL ComputeNu4CollisionTerm( J0_E, CollisionTerm(1:imat), &
        E_N, nE_G, iEk, nFlavors)

    !Flip them to get the Collision term for antimatter
    DO iE = 1,nE_G

      J0_E(2,iE,:,:) = PackIntoMatrix( J0(iE,1:imat),  nSpecies, nFlavors )
      J0_E(1,iE,:,:) = PackIntoMatrix( J0(iE,imat+1:), nSpecies, nFlavors )

    END DO

    CALL ComputeNu4CollisionTerm( J0_E, CollisionTerm(imat+1:), &
        E_N, nE_G, iEk, nFlavors)

  END SUBROUTINE ComputeNu4CollisionTerms


  SUBROUTINE ComputeNu4CollisionTerm( J0_E, CollisionTerm, &
          E_N, nE_G, iEk, nFlavors)

    INTEGER,  INTENT(IN)  :: nE_G, iEk
    REAL(DP), INTENT(IN)  :: J0_E         (2,nE_G,nFlavors,nFlavors)
    REAL(DP), INTENT(IN)  :: E_N(nE_G)
    REAL(DP), INTENT(OUT) :: CollisionTerm(nSpecies/2)

    REAL(DP) :: FirstTerm      (nFlavors,nFlavors)
    REAL(DP) :: SecondTerm     (nFlavors,nFlavors)
    REAL(DP) :: FirstTermG     (nFlavors,nFlavors)
    REAL(DP) :: SecondTermG    (nFlavors,nFlavors)
    REAL(DP) :: AntiCommutator1(nFlavors,nFlavors)
    REAL(DP) :: AntiCommutator2(nFlavors,nFlavors)
    REAL(DP) :: CollisionMatrix(nFlavors,nFlavors)
    REAL(DP) :: tmp            (nFlavors,nFlavors)
    REAL(DP) :: Identity       (nFlavors,nFlavors)
    REAL(DP) :: J0_E2(2,nFlavors,nFlavors)
    REAL(DP) :: E2

    INTEGER  :: iE, iE1, iE3

    !Consider defining Identity in some module without revaluating it every time
    Identity(:,:) = Zero
    DO iS = 1, nFlavors
      Identity(iS,iS) = One
    END DO
    
    FirstTerm (:,:) = Zero
    SecondTerm(:,:) = Zero

    FirstTermG (:,:) = Zero
    SecondTermG(:,:) = Zero

    DO iE1 = 1,nE_G
    DO iE3 = 1,nE_G
    
      E2 = E_N(iE1) + E_N(iE3) - E_N(iEk)
    
      IF ( (E2 >= E_N(1)) .AND. (E2 < E_N(nE_G)) ) THEN

        DO iS1 = 1,nFlavors
          DO iS2 = 1,nFlavors
  
            CALL QuadraticInterpolation1D( E_N, J0_E(1,:,iS1,iS2), nE_G, &
                  E2, J0_E2(1,iS1,iS2) )
  
            CALL QuadraticInterpolation1D( E_N, J0_E(2,:,iS1,iS2), nE_G, &
                  E2, J0_E2(2,iS1,iS2) )
  
          END DO
        END DO

        ! Scattering on nu, i.e. 3rd line in eq. 96 of Blaschke (2016)
        tmp = MatMul( Identity - J0_E(1,iE1,:,:) , J0_E2(1,:,:) )
        FirstTerm = FirstTerm + MatMul( TraceI( tmp, nFlavors ) + tmp , &
                      Identity - J0_E(1,iE3,:,:) ) * &
                         KernelNu4Scat(iE3,iE1,iEk)
        
        ! Scattering on nubar, i.e. first part of 5th line in eq. 96 of Blaschke (2016)
        tmp = MatMul( J0_E2(2,:,:) , Identity - J0_E(2,iE1,:,:) )
        SecondTerm = SecondTerm + MatMul( TraceI( tmp, nFlavors ) + tmp , &
                        Identity - J0_E(1,iE3,:,:) ) * &
                           KernelNu4Pair(iE3,iE1,iEk)

        ! Scattering on nubar, i.e. second part of 5th line in eq. 96 of Blaschke (2016)
        tmp = MatMul( Identity - J0_E(1,iE3,:,:) , Identity - J0_E(2,iE1,:,:) )
        SecondTerm = SecondTerm + &
              MatMul( TraceI( tmp, nFlavors ) + tmp , J0_E2(2,:,:) ) * &
                KernelNu4Pair(iE3,iE1,iEk)

        ! Gain part
        ! Scattering on nu, i.e. 3rd line in eq. 96 of Blaschke (2016)
        tmp = MatMul( J0_E(1,iE1,:,:) , Identity - J0_E2(1,:,:) )
        FirstTermG = FirstTermG + MatMul( TraceI( tmp, nFlavors ) + tmp , &
                        J0_E(1,iE3,:,:) ) * KernelNu4Scat(iE3,iE1,iEk)

        ! Scattering on nubar, i.e. first part of 5th line in eq. 96 of Blaschke (2016)
        tmp = MatMul( Identity - J0_E2(2,:,:) , J0_E(2,iE1,:,:) )
        SecondTermG = SecondTermG + MatMul( TraceI( tmp, nFlavors ) + tmp , &
                         J0_E(1,iE3,:,:) ) * KernelNu4Pair(iE3,iE1,iEk)

        ! Scattering on nubar, i.e. second part of 5th line in eq. 96 of Blaschke (2016)
        tmp = MatMul( J0_E(1,iE3,:,:) , J0_E(2,iE1,:,:) )
        SecondTermG = SecondTermG + MatMul( TraceI( tmp, nFlavors ) + tmp , &
                         Identity - J0_E2(2,:,:) ) * &
                            KernelNu4Pair(iE3,iE1,iEk)

      ENDIF

    END DO
    END DO
    
    AntiCommutator1 = MatMul( FirstTerm       , J0_E(1,iEk,:,:) ) + &
                      MatMul( J0_E(1,iEk,:,:) , FirstTerm       )

    AntiCommutator2 = MatMul( SecondTerm      , J0_E(1,iEk,:,:) ) + &
                      MatMul( J0_E(1,iEk,:,:) , SecondTerm      )

    CollisionMatrix = - ( AntiCommutator1 + AntiCommutator2 )

    !Gain part
    AntiCommutator1 = MatMul( FirstTermG       , Identity - J0_E(1,iEk,:,:) ) + &
                      MatMul( Identity - J0_E(1,iEk,:,:) , FirstTermG       )

    AntiCommutator2 = MatMul( SecondTermG      , Identity - J0_E(1,iEk,:,:) ) + &
                      MatMul( Identity - J0_E(1,iEk,:,:) , SecondTermG      )

    CollisionMatrix = CollisionMatrix + ( AntiCommutator1 + AntiCommutator2 )

    CollisionTerm   = FlattenMatrix( CollisionMatrix, nFlavors, nSpecies / 2 )

  END SUBROUTINE ComputeNu4CollisionTerm



  SUBROUTINE ComputeNu4Kernels

    REAL(DP), ALLOCATABLE :: E_N(:), E_Left(:), E_Right(:)
    INTEGER  :: nE_G
    
    INTEGER  :: iEk, iE1, iE3
    INTEGER  :: nE, iE, iNodeE
    REAL(DP) :: qk, q1, q2, q3, E2
    REAL(DP) :: V1, V3 ! phase-space volumes (cm^-3)
    REAL(DP) :: hbarc_mevcm = 1.97326966d-11 
    
    nE = iE_E0 - iE_B0 + 1
    nE_G = nNodesE * nE

    ALLOCATE( KernelNu4Pair(nE_G,nE_G,nE_G) )
    ALLOCATE( KernelNu4Scat(nE_G,nE_G,nE_G) )

    ALLOCATE( E_N    (nE_G) )
    ALLOCATE( E_Left (nE_G) )
    ALLOCATE( E_Right(nE_G) )
    
    DO iN_E = 1, nE_G

      iE     = MOD( (iN_E-1) / nDOFE, nE    ) + iE_B0
      iNodeE = MOD( (iN_E-1)        , nDOFE ) + 1

      E_N(iN_E) = NodeCoordinate( MeshE, iE, iNodeE )

      ! GET ENERGIES AT INTERFACE !--- WARING ---! ONLY WORKS WITH 1 NODE
      E_Left (iN_E) = E_N(iN_E) - Half * MeshE % Width(iN_E)
      E_Right(iN_E) = E_N(iN_E) + Half * MeshE % Width(iN_E)

    ENDDO

    KernelNu4Pair(:,:,:) = Zero
    KernelNu4Scat(:,:,:) = Zero

    DO iEk=1,nE_G
       
       qk = E_N(iEk) / MeV
       
       DO iE1=1,nE_G
         
         q1 = E_N(iE1) / MeV
         
         V1 = FourPi * ( E_Right(iE1)**3 - E_Left(iE1)**3 ) / Three  &
                     / ( TwoPi )**3 / MeV**3 / hbarc_mevcm**3      
        
         DO iE3=1,nE_G
           
           q3 = E_N(iE3) / MeV

           V3 = Fourpi * ( E_Right(iE3)**3 - E_Left(iE3)**3 ) / Three &
                       / ( TwoPi )**3 / MeV**3 / hbarc_mevcm**3

           E2 = E_N(iE1) + E_N(iE3) - E_N(iEk)

           ! result is 1/cm
           IF ( E2 > E_N(1) .AND. E2 <= E_N(nE_G) ) THEN
             q2 = E2 / MeV
             
             KernelNu4Pair(iE3,iE1,iEk) = Nu4Pair_Kernel_Single( qk, q1, q2, q3) * V1 * V3
             
             KernelNu4Scat(iE3,iE1,iEk) = Nu4Scat_Kernel_Single( qk, q1, q2, q3) * V1 * V3
             
             !NOW FIX UNITS SINCE SHERWOOD SOLVES 1/c*df/dt
             KernelNu4Pair(iE3,iE1,iEk) = KernelNu4Pair(iE3,iE1,iEk) * SpeedOfLightCGS * One / Second 
             KernelNu4Scat(iE3,iE1,iEk) = KernelNu4Scat(iE3,iE1,iEk) * SpeedOfLightCGS * One / Second
        
             IF( KernelNu4Pair(iE3,iE1,iEk) < 0 ) THEN
                WRITE(*,*) E2, iEk, iE1, iE3, KernelNu4Pair(iE3,iE1,iEk)
                STOP "KernelNu4Pair value less than 0"
             END IF

             IF( KernelNu4Scat(iE3,iE1,iEk) < 0 ) THEN
                WRITE(*,*) E2, iEk, iE1, iE3, KernelNu4Scat(iE3,iE1,iEk)
                STOP "KernelNu4Scat value less than 0"
             END IF
    
           ELSE
 
             KernelNu4Scat(iE3,iE1,iEk) = Zero
             KernelNu4Pair(iE3,iE1,iEk) = Zero
           
           ENDIF
        
         END DO
      END DO
    END DO

  DEALLOCATE( E_N, E_right, E_left )

  END SUBROUTINE ComputeNu4Kernels

  !!=============================!!
  !! nu+nu <--> nu+nu scattering !!
  !!=============================!!
  ! second line in Blaschke & Cirigliano 2016 Equation 96
  ! technically applicable only in isotropic limit
  FUNCTION Nu4Scat_Kernel_Single(qk, q1, q2, q3) result(Kernel)

    REAL(DP), INTENT(IN) :: qk, q1, q2, q3
    REAL(DP) :: Kernel
    REAL(DP) :: Gfermi = 1.16637d-11 !MeV^-2
    REAL(DP) :: hbarc_mevcm = 1.97326966d-11

    IF( ABS( q1-q2+q3-qk ) > ( q1+q2+q3+qk )*1.d-10) THEN
       WRITE(*,*) qk, q1, q2, q3
       STOP "ERROR: passed q2 Does not conserve energy"
    ENDIF

    Kernel = Gfermi**2 * Pi * hbarc_mevcm **5 * &
      (                      nu4_D3( q1,q2,q3,qk )   &
         + q1 * q2 * q3 * qk * nu4_D1( q1,q2,q3,qk )   &
         + q2 * qk           * nu4_D2( q2,qk,q1,q3 )   &
         + q1 * q3           * nu4_D2( q1,q3,q2,qk ) ) &
                  / q1**2 / q3**2 / qk**2
    
  END FUNCTION Nu4Scat_Kernel_Single

  !!=====================================!!
  !! nu+nubar <--> nu+nubar annihilation !!
  !!=====================================!!
  ! fourth line in Blaschke & Cirigliano 2016 Equation 96
  ! technically applicable only in isotropic limit
  FUNCTION Nu4Pair_Kernel_Single(qk, q1, q2, q3) result(Kernel)
    
    REAL(DP), INTENT(IN) :: qk, q1, q2, q3
    REAL(DP) :: Kernel
    REAL(DP) :: Gfermi = 1.16637d-11 !MeV^-2 
    REAL(DP) :: hbarc_mevcm = 1.97326966d-11

    IF( ABS( q1-q2+q3-qk ) > ( q1+q2+q3+qk )*1.d-10) THEN
       WRITE(*,*) qk, q1, q2, q3
       STOP "ERROR: passed q2 Does not conserve energy"
    ENDIF
    
    Kernel = Gfermi**2 * Pi * hbarc_mevcm **5 * &
        (                      nu4_D3( q1,q2,q3,qk )   &
         + q1 * q2 * q3 * qk * nu4_D1( q1,q2,q3,qk )   &
         - q1 * qk           * nu4_D2( q1,qk,q2,q3 )   &
         - q2 * q3           * nu4_D2( q2,q3,q1,qk ) ) &
                  / q1**2 / q3**2 / qk**2

  END FUNCTION Nu4Pair_Kernel_Single

  ! equation D4c in Blaschke+Cirigliano 2016
  FUNCTION nu4_D4c(q1,q2,q3,q4) RESULT(D4c) ! MeV^5
    
    REAL(DP), INTENT(IN) :: q1, q2, q3, q4
    REAL(DP) :: D4c

    D4c = One/60.0 * ( q1**5 - Five * q1**3*q2**2 + Five * q1**2*q2**3 - &
                       q2**5 - Five * q1**3*q3**2 + Five * q2**3*q3**2 + &
                               Five * q1**2*q3**3 + Five * q2**2*q3**3 - &
                       q3**5 - Five * q1**3*q4**2 + Five * q2**3*q4**2 + &
                               Five * q3**3*q4**2 + Five * q1**2*q4**3 + &
                               Five * q2**2*q4**3 + Five * q3**2*q4**3 - &
                       q4**5 )
  
  END FUNCTION nu4_D4c
  
  !!====!!
  !! D1 !!
  !!====!!
  FUNCTION nu4_D1(qin1,qin2,qin3,qin4) RESULT(D1)! MeV
  
    REAL(DP), INTENT(IN) :: qin1, qin2, qin3, qin4
    REAL(DP) :: q1, q2, q3, q4, D1

    ! all RESULTs are symmetric on 1<-->2.AND.3<-->4
    ! use the case where q1>q2.AND.q3>q4
    
    q1 = MAX(qin1,qin2)
    q2 = MIN(qin1,qin2)
    q3 = MAX(qin3,qin4)
    q4 = MIN(qin3,qin4)

    ! Case 1 (Eqns. D4)
    IF( q1+q2 >= q3+q4 .AND. q1+q4 >= q2+q3 ) THEN
      
      IF( q1 <= q2+q3+q4 ) THEN
        D1 = Half  * ( q2+q3+q4-q1 )
      ELSE
        D1 = Zero
      END IF

    ! Case 2 (Eqns. D5)
    ELSE IF( q1+q2 >= q3+q4 .AND. q1+q4 <= q2+q3 ) THEN
      
      D1 = q4

    ! Case 3 (Eqns. D6)
    ELSE IF( q1+q2 <= q3+q4 .AND. q1+q4 <= q2+q3 ) THEN
      
      IF( q3 <= q1+q2+q4 ) THEN
        D1 = Half * ( q1+q2+q4-q3 )
      ELSE
        D1 = Zero
      END IF

    ! Case 4 (Eqns. D7)
    ELSE !(q1+q2<q3+q4 .AND. q1+q4>q2+q3)
      
      D1 = q2
    
    END IF
  
  END FUNCTION nu4_D1
  
  !!====!!
  !! D2 !!
  !!====!!
  FUNCTION nu4_D2(qin1,qin2,qin3,qin4) RESULT(D2)! MeV^3
  
    REAL(DP), INTENT(IN) :: qin1, qin2, qin3, qin4
    REAL(DP) :: q1, q2, q3, q4, D2

    ! all RESULTs are symmetric on 1<-->2.AND.3<-->4
    ! use the case where q1>q2.AND.q3>q4
    q1 = MAX(qin1,qin2)
    q2 = MIN(qin1,qin2)
    q3 = MAX(qin3,qin4)
    q4 = MIN(qin3,qin4)

    ! Case 1 (Eqns. D4)
    IF(q1+q2 >= q3+q4 .AND. q1+q4>=q2+q3) THEN
      
      IF(q1<=q2+q3+q4) THEN
        D2 = One/12.0 * ( (q1-q2)**3 + Two * ( q3**3+q4**3 ) &
                          - Three * (q1-q2) * ( q3**2+q4**2) )
      ELSE
        D2 = Zero
      END IF

    ! Case 2 (Eqns. D5)
    ELSE IF ( q1+q2 >= q3+q4 .AND. q1+q4 <= q2+q3 ) THEN
      
      D2 = q4**3 / Three

    ! Case 3 (Eqns. D6)
    ELSE IF( q1+q2 <= q3+q4 .AND. q1+q4 <= q2+q3 ) THEN
    
      IF( q3<=q1+q2+q4 ) THEN
        D2 = One/12.0 * ( - (q1+q2)**3 - Two * q3**3 + &
            Two * q4**3 + Three * (q1+q2) * (q3**2+q4**2) )
      ELSE
        D2 = Zero
      END IF

    ! Case 4 (Eqns. D7)
    ELSE !(q1+q2<q3+q4 .AND. q1+q4>q2+q3)
    
      D2 = q2/6.0 * ( Three * q3**2 + Three * q4**2 - &
                      Three * q1**2 -         q2**2 )
    
    END IF

  END FUNCTION nu4_D2
  
  !!====!!
  !! D3 !!
  !!====!!
  FUNCTION nu4_D3(qin1,qin2,qin3,qin4) RESULT(D3)! MeV^5
  
    REAL(DP), INTENT(IN) :: qin1, qin2, qin3, qin4
    REAL(DP) :: q1, q2, q3, q4, D3

    ! all RESULTs are symmetric on 1<-->2.AND.3<-->4
    ! use the case where q1>q2.AND.q3>q4
    q1 = MAX(qin1,qin2)
    q2 = MIN(qin1,qin2)
    q3 = MAX(qin3,qin4)
    q4 = MIN(qin3,qin4)

    ! Case 1 (Eqns. D4)
    IF( q1+q2 >= q3+q4 .AND. q1+q4 >= q2+q3 ) THEN
      
      IF( q1 <= q2+q3+q4 ) THEN
        D3 = nu4_D4c(q1,q2,q3,q4)
      ELSE
        D3 = Zero
      END IF

    ! Case 2 (Eqns. D5)
    ELSE IF( q1+q2 >= q3+q4 .AND. q1+q4 <= q2+q3 ) THEN
      
      D3 = q4**3 / 30.0 * ( Five * q1**2 + Five * q2**2 + &
                            Five * q3**2 -        q4**2 )
    
    ! Case 3 (Eqns. D6)
    ELSE IF( q1+q2 <= q3+q4 .AND. q1+q4 <= q2+q3 ) THEN
        
      IF( q3 <= q1+q2+q4 ) THEN
        D3 = nu4_D4c(q3,q4,q1,q2)
      ELSE
        D3 = Zero
      END IF

    ! Case 4 (Eqns. D7)
    ELSE !(q1+q2<q3+q4 .AND. q1+q4>q2+q3)
      
      D3 = q2**3/30.0 * ( Five *q1**2 + Five * q3**2 + &
                          Five *q4**2 -        q2**2 )
    END IF
  
  END FUNCTION nu4_D3


  FUNCTION PackIntoMatrix( Vector, dimV, dimM ) RESULT( Matrix )
    
    INTEGER,  INTENT(IN)  :: dimV, dimM
    REAL(DP), INTENT(IN)  :: Vector( dimV )
    REAL(DP)              :: Matrix( dimM,dimM )
    
    ! THIS WORKS ONLY IF THE ORIGINAL FLAVOR MATRIX WAS NUMBERED
    ! ROW-WISE : i.e. 1 2 3      1 2 
    !                   4 5  AND   3
    !                     6

    ! FILL UPPER PART
    DO i = 1,dimM
      DO j = i, dimM

        IF ( j == 3 ) THEN
          Matrix(i,j) = Vector(i+j)
        ELSE 
          Matrix(i,j) = Vector(i+j-1)
        END IF
      
      END DO
    END DO
    
    ! FILL LOWER PART
    DO i = 1,dimM-1
      DO j = i+1,dimM
        
        IF ( j == 3 ) THEN
          Matrix(j,i) = Vector(i+j)
        ELSE
          Matrix(j,i) = Vector(i+j-1)
        END IF

      END DO
    END DO

  END FUNCTION PackIntoMatrix
    
  FUNCTION FlattenMatrix( Matrix, dimM, dimV ) RESULT( Vector )

    INTEGER,  INTENT(IN)  :: dimV, dimM
    REAL(DP), INTENT(IN)  :: Matrix( dimM,dimM )
    REAL(DP)              :: Vector( dimV )

    ! THIS WORKS ONLY IF THE ORIGINAL FLAVOR MATRIX WAS NUMBERED
    ! ROW-WISE : i.e. 1 2 3      1 2 
    !                   4 5  AND   3
    !                     6

    ! FLATTEN UPPER PART, i.e. the one that is stored
    DO i = 1,dimM
      DO j = i, dimM

        IF ( j == 3 ) THEN
          Vector(i+j) = Matrix(i,j)
        ELSE
          Vector(i+j-1) = Matrix(i,j)
        END IF

      END DO
    END DO

  END FUNCTION FlattenMatrix


  FUNCTION TraceI( Matrix, dimM ) RESULT( Trace )
    
    INTEGER,  INTENT(IN)  :: dimM
    REAL(DP), INTENT(IN)  :: Matrix(dimM,dimM)
    REAL(DP)              :: Trace (dimM,dimM)

    INTEGER :: i, j

    Trace(:,:) = Zero

    DO i = 1,dimM
      DO j = 1,dimM

        Trace(i,i) = Trace(i,i) + Matrix(j,j)

      END DO
    END DO

  END FUNCTION TraceI

  FUNCTION MatrixMul( Matrix1, Matrix2, dimM) RESULT(Matrix3)

    INTEGER,  INTENT(IN) :: dimM
    REAL(DP), INTENT(IN) :: Matrix1(dimM,dimM), Matrix2(dimM,dimM)
    REAL(DP)             :: Matrix3(dimM,dimM)
    
    INTEGER :: i,j,k
    
    DO j = 1,dimM
      DO k = 1,dimM
        DO i = 1,dimM
        
          Matrix3(i,j) = Matrix1(i,j) + Matrix1(i,k) * Matrix2(k,j)

        END DO
      END DO
    END DO
  
  END FUNCTION MatrixMul


END MODULE FlavorOpacitiesModule
