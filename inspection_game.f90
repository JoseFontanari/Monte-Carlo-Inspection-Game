!*******************************************************************************
!* DATE: 18/10/2025
!* PROJECT: Finite Population Dynamics Resolve the Central Paradox of the Inspection Game
!* AUTHORS: Bianca Ishikawa & José Fontanari
!* JOURNAL: European Physical Journal B (EPJ B)
!*
!* DESCRIPTION:
!* This program simulates the imitation dynamics for the Inspection Game in 
!* finite populations. It explores the stochastic evolution of strategies 
!* between two asymmetric populations: Citizens and Inspectors.
!*
!* DYNAMICS:
!* - Pairwise comparison (Imitation dynamics).
!* - A focal agent imitates a model agent with probability proportional 
!* to the payoff difference.
!*
!* SYSTEM STATE:
!* - Absorbing states (Fixation) are detected to calculate fixation probabilities.
!*******************************************************************************

PROGRAM Inspection_Game
    IMPLICIT NONE

    !---------------------------------------------------------------------------
    ! PARAMETERS
    !---------------------------------------------------------------------------
    INTEGER, PARAMETER :: N = 1000      ! Total number of Citizens
    INTEGER, PARAMETER :: M = 100       ! Total number of Inspectors
    INTEGER, PARAMETER :: Nsample = 10000 ! Number of Monte Carlo samples per parameter set
    INTEGER, PARAMETER :: Ntime = 1000000 ! Maximum time steps before forcing stop

    !---------------------------------------------------------------------------
    ! VARIABLES
    !---------------------------------------------------------------------------
    LOGICAL :: S(N)  ! Citizen Strategies: TRUE = Crime, FALSE = No Crime
    LOGICAL :: R(M)  ! Inspector Strategies: TRUE = Inspect, FALSE = No Inspect
    
    INTEGER :: X, Y  ! Current count of Criminals (X) and Inspecting Inspectors (Y)
    INTEGER :: X0, Y0 ! Initial conditions
    
    ! Payoff Parameters
    REAL :: g, p     ! g = Gain from crime, p = Penalty
    REAL :: rr, kk   ! rr = Reward for catching, kk = Cost of inspection
    REAL :: D        ! Normalization factor for switching probability
    REAL :: gg       ! Loop variable for varying the ratio g/p
    
    ! Random Number Generation
    REAL :: r1, aux
    INTEGER :: is, ir, js, jr ! Indices for focal (i) and model (j) agents
    REAL :: fis, fir, fjs, fjr ! Payoffs for focal and model agents
    
    ! Statistics Counters (Fixation events)
    INTEGER :: ix0ym, ixnym, ix0y0, ixny0
    REAL    :: tx0ym, txnym, tx0y0, txny0 ! Time to fixation accumulators
    
    ! Loop counters
    INTEGER :: isample, itime, jtime, k

    !---------------------------------------------------------------------------
    ! INITIALIZATION
    !---------------------------------------------------------------------------
    ! Initialize Random Number Generator
    CALL init_random_seed()

    ! Define fixed payoff parameters (k=1 is the normalization unit)
    rr = 4.0
    kk = 1.0
    p  = 100.0

    ! Initial Conditions (50% of each population)
    X0 = N / 2
    Y0 = M / 2

    !---------------------------------------------------------------------------
    ! PARAMETER SWEEP LOOP
    ! Here we vary 'g' relative to 'p' to sweep the Inspection Threshold g/p
    !---------------------------------------------------------------------------
    DO gg = 0.01, 0.99, 0.01
        g = p * gg

        ! Calculate Normalization Factor D
        ! D ensures that the transition probability (Delta_payoff / D) <= 1.
        ! D = max(payoff differences)
        IF ((p - g) > (rr - kk)) THEN
            D = p - g
        ELSE
            D = rr - kk
        END IF
        IF (D < g)  D = g
        IF (D < kk) D = kk

        ! Reset Statistics for this parameter set
        ix0ym = 0; ixnym = 0; ix0y0 = 0; ixny0 = 0
        tx0ym = 0.0; txnym = 0.0; tx0y0 = 0.0; txny0 = 0.0

        !-----------------------------------------------------------------------
        ! MONTE CARLO SAMPLING LOOP
        !-----------------------------------------------------------------------
        sample_loop: DO isample = 1, Nsample

            ! Reset Populations to Initial Conditions
            S(1:X0)   = .TRUE.   ! Criminals
            S(X0+1:N) = .FALSE.  ! Law-abiding
            R(1:Y0)   = .TRUE.   ! Inspecting
            R(Y0+1:M) = .FALSE.  ! Not Inspecting

            X = X0
            Y = Y0

            !-------------------------------------------------------------------
            ! TIME EVOLUTION LOOP
            !-------------------------------------------------------------------
            time_loop: DO itime = 1, Ntime
                
                ! Perform N interaction steps per time unit (MCS)
                step_loop: DO jtime = 1, N
                    
                    ! 1. Select Focal Agents (Randomly)
                    CALL RANDOM_NUMBER(r1)
                    is = 1 + INT(N * r1)
                    CALL RANDOM_NUMBER(r1)
                    ir = 1 + INT(M * r1)

                    ! 2. Select Model Agents (Randomly, distinct from focal)
                    DO
                        CALL RANDOM_NUMBER(r1)
                        js = 1 + INT(N * r1)
                        IF (js /= is) EXIT
                    END DO
                    DO
                        CALL RANDOM_NUMBER(r1)
                        jr = 1 + INT(M * r1)
                        IF (jr /= ir) EXIT
                    END DO

                    ! 3. Calculate Payoffs
                    ! Payoff depends on the interaction between Citizen and Inspector
                    
                    ! Focal Citizen Payoff
                    IF (S(is)) THEN ! If Citizen is Criminal
                        IF (R(ir)) THEN ! And Inspector Inspects
                            fis = g - p
                            fir = rr - kk
                        ELSE            ! And Inspector does NOT Inspect
                            fis = g
                            fir = 0.0
                        END IF
                    ELSE            ! If Citizen is Law-abiding
                        IF (R(ir)) THEN
                            fis = 0.0
                            fir = -kk ! Cost of inspection wasted
                        ELSE
                            fis = 0.0
                            fir = 0.0
                        END IF
                    END IF

                    ! Model Citizen Payoff (plays against Model Inspector)
                    IF (S(js)) THEN
                        IF (R(jr)) THEN
                            fjs = g - p
                            fjr = rr - kk
                        ELSE
                            fjs = g
                            fjr = 0.0
                        END IF
                    ELSE
                        IF (R(jr)) THEN
                            fjs = 0.0
                            fjr = -kk
                        ELSE
                            fjs = 0.0
                            fjr = 0.0
                        END IF
                    END IF

                    ! 4. Imitation Dynamics (Update Strategies)
                    
                    ! Citizen Update
                    IF (S(is) .NEQV. S(js)) THEN ! Only if strategies differ
                        IF (fjs > fis) THEN      ! Only if model is more successful
                            aux = (fjs - fis) / D
                            CALL RANDOM_NUMBER(r1)
                            IF (r1 < aux) THEN
                                ! Switch strategy
                                IF (S(is)) THEN
                                    X = X - 1 ! Was Crime, becomes No Crime
                                ELSE
                                    X = X + 1 ! Was No Crime, becomes Crime
                                END IF
                                S(is) = S(js)
                            END IF
                        END IF
                    END IF

                    ! Inspector Update
                    IF (R(ir) .NEQV. R(jr)) THEN
                        IF (fjr > fir) THEN
                            aux = (fjr - fir) / D
                            CALL RANDOM_NUMBER(r1)
                            IF (r1 < aux) THEN
                                IF (R(ir)) THEN
                                    Y = Y - 1 ! Was Inspect, becomes No Inspect
                                ELSE
                                    Y = Y + 1 ! Was No Inspect, becomes Inspect
                                END IF
                                R(ir) = R(jr)
                            END IF
                        END IF
                    END IF

                    ! 5. Check for Fixation (Absorbing States)
                    
                    ! State: All Crime, All Inspect
                    IF (X == N .AND. Y == M) THEN
                        ixnym = ixnym + 1
                        txnym = txnym + REAL(itime) * REAL(jtime) / REAL(N)
                        EXIT time_loop
                    END IF

                    ! State: All Crime, No Inspect
                    IF (X == N .AND. Y == 0) THEN
                        ixny0 = ixny0 + 1
                        txny0 = txny0 + REAL(itime) * REAL(jtime) / REAL(N)
                        EXIT time_loop
                    END IF

                    ! State: No Crime, No Inspect
                    IF (X == 0 .AND. Y == 0) THEN
                        ix0y0 = ix0y0 + 1
                        tx0y0 = tx0y0 + REAL(itime) * REAL(jtime) / REAL(N)
                        EXIT time_loop
                    END IF

                    ! State: No Crime, All Inspect (Extinction of Crime)
                    IF (X == 0 .AND. Y == M) THEN
                        ix0ym = ix0ym + 1
                        tx0ym = tx0ym + REAL(itime) * REAL(jtime) / REAL(N)
                        EXIT time_loop
                    END IF

                END DO step_loop
            END DO time_loop

        END DO sample_loop

        !-----------------------------------------------------------------------
        ! DATA OUTPUT
        ! Appends results to file for each parameter g
        !-----------------------------------------------------------------------
        OPEN(UNIT=10, FILE='simg_N1000_M100_r4_p100.dat', STATUS='UNKNOWN', ACCESS='APPEND')
        WRITE(10, *) rr, ix0y0, ix0ym, ixny0, ixnym, tx0y0, &
                     tx0ym, txny0, txnym, X0, Y0, N, M, g, p, Nsample
        CLOSE(10)

    END DO

END PROGRAM Inspection_Game


!*******************************************************************************
! SUBROUTINE: init_random_seed
! DESCRIPTION: Initializes the random number generator using the system clock
! to ensure different seeds on different runs.
!*******************************************************************************
SUBROUTINE init_random_seed()
    IMPLICIT NONE
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    ! Seed based on clock and index to ensure entropy
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)

    DEALLOCATE(seed)
IEND SUBROUTINE init_random_seed
