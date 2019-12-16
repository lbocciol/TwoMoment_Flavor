PROGRAM Decoherence

  USE KindModule, ONLY: &
    DP, Zero, One, Two, &
    Pi, TwoPi, &
    Half
  USE UnitsModule, ONLY: &
    Kilometer, MeV, Millisecond, &
    Second
  USE ProgramHeaderModule, ONLY: &
    nZ, nNodesZ, &
    iX_B0, iX_E0, iX_B1, iX_E1, &
    iE_B0, iE_E0, iE_B1, iE_E1, &
    iZ_B0, iZ_E0, iZ_B1, iZ_E1
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE TimersModule, ONLY: &
    InitializeTimers, &
    FinalizeTimers, &
    TimersStart, &
    TimersStop, &
    Timer_InputOutput, &
    Timer_Evolve
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  USE ReferenceElementModuleE, ONLY: &
    InitializeReferenceElementE, &
    FinalizeReferenceElementE
  USE ReferenceElementModuleE_Lagrange, ONLY: &
    InitializeReferenceElementE_Lagrange, &
    FinalizeReferenceElementE_Lagrange
  USE ReferenceElementModule, ONLY: &
    InitializeReferenceElement, &
    FinalizeReferenceElement
  USE ReferenceElementModule_Lagrange, ONLY: &
    InitializeReferenceElement_Lagrange, &
    FinalizeReferenceElement_Lagrange
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE GeometryComputationModuleE, ONLY: &
    ComputeGeometryE
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE GeometryFieldsModuleE, ONLY: &
    uGE
  USE FluidFieldsModule, ONLY: &
    uCF
  USE RadiationFieldsModule, ONLY: &
    uCR
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE OpacityModule_TABLE, ONLY: &
    InitializeOpacities_TABLE, &
    FinalizeOpacities_TABLE
  USE NeutrinoOpacitiesModule, ONLY: &
    CreateNeutrinoOpacities, &
    DestroyNeutrinoOpacities
  USE NeutrinoOpacitiesComputationModule, ONLY: &
    ComputeNeutrinoOpacities  
  USE TwoMoment_ClosureModule, ONLY: &
    InitializeClosure_TwoMoment
  USE IntegrationModule, ONLY: &
    TimeStepping_Decoherence
  USE InitializationModule, ONLY: &
    InitializeFields_Decoherence, &
    InitializeSingleZone_Decoherence
  
  IMPLICIT NONE

  LOGICAL  :: wrt
  INTEGER  :: nFlavors, nOpacities
  INTEGER  :: iCycle, iCycleD, iCycleW
  INTEGER  :: wrt_control
  INTEGER  :: nNodes, nSpecies
  INTEGER  :: nE, bcE, nX(3), bcX(3)
  REAL(DP) :: eL, eR, xL(3), xR(3), ZoomX(3)
  REAL(DP) :: t, dt, t_end
  REAL(DP) :: dt_max

  nNodes = 2
  nSpecies = 6

  nFlavors = 2
  nOpacities = 2

  nX  = [ 128, 1, 1 ]
  xL  = [ 0.0d0 * Kilometer, 0.0_DP, 0.0_DP ]
  xR  = [ 5.0d2 * Kilometer, Pi,     TwoPi  ]
  bcX = [ 1, 0, 0 ]
  ZoomX = [ 1.10851_DP, 1.0_DP, 1.0_DP ]

  nE  = 50
  eL  = 0.5d0 * MeV
  eR  = 200.5d0 * MeV
  bcE = 0

  t       = 0.0_DP
  t_end   = 1.0      * Second
  dt      = 1.0d-15  * Second
  
  iCycleD = 1
  iCycleW = 10

  CALL InitializeProgram &
         ( ProgramName_Option &
             = 'Decoherence', &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 1, 1, 1 ], &
           bcX_Option &
             = bcX, &
           xL_Option &
             = xL, &
           xR_Option &
             = xR, &
           ZoomX_Option &
             = ZoomX, &
           nE_Option &
             = nE, &
           swE_Option &
             = 0, &
           bcE_Option &
             = bcE, &
           eL_Option &
             = eL, &
           eR_Option &
             = eR, &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = 'SPHERICAL', &
           ActivateUnits_Option &
             = .TRUE., &
           nSpecies_Option &
             = nSpecies, &
           BasicInitialization_Option &
             = .TRUE. )

  ! --- Position Space Reference Element and Geometry ---

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF )

  ! --- Energy Space Reference Element and Geometry ---

  CALL InitializeReferenceElementE

  CALL InitializeReferenceElementE_Lagrange

  CALL ComputeGeometryE &
         ( iE_B0, iE_E0, iE_B1, iE_E1, uGE )

  ! --- Phase Space Reference Element ---

  CALL InitializeReferenceElement

  CALL InitializeReferenceElement_Lagrange

  ! --- Initialize Equation of State ---

  CALL InitializeEquationOfState &
         ( EquationOfState_Option &
             = 'TABLE', &
           EquationOfStateTableName_Option &
             = 'EquationOfStateTable.h5' )
  
  ! --- Initialize and create Opacities

  CALL InitializeOpacities_TABLE &                         
         ( OpacityTableName_EmAb_Option &                  
             = 'wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5', &    
           OpacityTableName_Iso_Option  &                  
             = 'wl-Op-SFHo-15-25-50-E40-B85-Iso.h5', &     
!           OpacityTableName_NES_Option &                   
!             = 'wl-Op-SFHo-15-25-50-E40-B85-NES.h5', &     
!           OpacityTableName_Pair_Option &                  
!             = 'wl-Op-SFHo-15-25-50-E40-B85-Pair.h5', &    
           Verbose_Option = .TRUE. )                       

  CALL CreateNeutrinoOpacities( nZ, nNodesZ, nOpacities )

  ! --- Initialize Moment Closure ---

  CALL InitializeClosure_TwoMoment

  ! --- Create Profile and initialize Moments ---
  
  CALL InitializeFields_Decoherence

!  CALL InitializeSingleZone_Decoherence

  CALL WriteFieldsHDF &
         ( Time = 0.0_DP, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE. )
  
  ! --- Evolve ---

  WRITE(*,*)
  WRITE(*,'(A6,A,ES8.2E2,A8,ES8.2E2)') &
    '', 'Evolving from t = ', t / Millisecond, &
    ' to t = ', t_end / Millisecond
  WRITE(*,*)

  iCycle = 0
  wrt     = .FALSE.
  wrt_control = 1
  DO WHILE( t < t_end )
    
    iCycle = iCycle + 1
    
    IF( t + dt > t_end )THEN

      dt = t_end - t

    END IF

    IF ( MOD( iCycle, iCycleW ) == 0 ) THEN
      
      wrt = .TRUE.

    ENDIF

!    IF ( MOD( iCycle, wrt_control ) == 0 ) THEN
!      
!      wrt = .TRUE.
!      wrt_control = wrt_control * 5
!
!    ENDIF

    IF( MOD( iCycle, iCycleD ) == 0 )THEN
      WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A5,A5,ES12.6E2,A4)') &
          '', 'Cycle = ', iCycle, &
          '', 't = ',  t  / Second, &
          ' (s) ', 'dt = ', dt / Second, ' (s)'

    END IF

    CALL TimersStart( Timer_Evolve )

    CALL TimeStepping_Decoherence(dt, nOpacities, nFlavors)
        
    t = t + dt

    CALL TimersStop( Timer_Evolve )

    IF( wrt )THEN

      CALL TimersStart( Timer_InputOutput )
      
      CALL WriteFieldsHDF &
             ( Time = t, &
               WriteGF_Option = .FALSE., &
               WriteFF_Option = .FALSE., &
               WriteRF_Option = .TRUE., &
               WriteOP_Option = .FALSE. )

      CALL TimersStop( Timer_InputOutput )

      wrt = .FALSE.

    END IF

  END DO

  ! --- Write Final Solution ---

  CALL TimersStart( Timer_InputOutput )

  CALL WriteFieldsHDF &
         ( Time = t, &
           WriteGF_Option = .TRUE., &
           WriteFF_Option = .TRUE., &
           WriteRF_Option = .TRUE., &
           WriteOP_Option = .TRUE. )

  CALL TimersStop( Timer_InputOutput )

  WRITE(*,*)
  WRITE(*,'(A6,A,I6.6,A,ES12.6E2,A)') &
    '', 'Finished ', iCycle, ' Cycles in ', Timer_Evolve, ' s'
  WRITE(*,*)

  CALL FinalizeTimers

  CALL FinalizeEquationOfState

  CALL FinalizeOpacities_TABLE
    
  CALL DestroyNeutrinoOpacities

  CALL FinalizeReferenceElementX

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeReferenceElementE

  CALL FinalizeReferenceElementE_Lagrange

  CALL FinalizeReferenceElement

  CALL FinalizeReferenceElement_Lagrange

  CALL FinalizeProgram

END PROGRAM Decoherence
