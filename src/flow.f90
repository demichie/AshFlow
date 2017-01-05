!********************************************************************************
!> \brief PDC flow model
!> This module calls subroutines required from other modules for solving the 
!> model and describes exit conditions for, and the information output, on run 
!> completion. 
!> It contains the loops required for iterative discretization used to 
!> approximate solutions to equations in solver flow module.
!********************************************************************************

MODULE PDCflow

  IMPLICIT NONE
  SAVE

CONTAINS

  !******************************************************************************
  !> \brief PDC run_name
  !> This subroutine calls all subroutines required from other modules for       
  !> solving the model. 
  !******************************************************************************

  SUBROUTINE pdc(run_name)

    ! External variables
    USE envi_module, ONLY: zmet, alpha
    USE gas_module, ONLY: gasvals
    USE mixture_module, ONLY: beta , beta_old, beta_new
    USE mixture_module, ONLY: Ri
    USE current_module, ONLY: r, r0, r_old, r_new, u, h_old, h_new, h,T
    USE solver_module, ONLY: dr, dr0, dr_old, f, ftemp, rhs, rhstemp
    USE solver_module, ONLY: f_stepold, f_oldold

    USE solver_module, ONLY: f_new
    USE solver_module, ONLY: f_RK
    USE solver_module, ONLY: n_RK
    USE solver_module, ONLY: rhs_RK
    USE solver_module, ONLY: rhs4
    USE solver_module, ONLY: f_err

    USE solver_module, ONLY: eps_rel , eps_abs

    ! External procedures
    USE current_module, ONLY: initialize_current
    USE mixture_module, ONLY: compute_mixture, flow_regime
    USE particles_module, ONLY: initialize_particles
    USE solver_module, ONLY: rate, lump, unlump, marching

    USE inpout, ONLY: write_output , write_dakota

    IMPLICIT NONE

    CHARACTER(LEN=80), INTENT(IN) :: run_name

    CHARACTER(LEN=4) :: idx_string

    REAL*8 :: volume_flux

    REAL*8 :: eps_step
    REAL*8 :: convergence
    
    REAL*8 :: Ri_old
    REAL*8 :: Ri_new

    INTEGER :: i_RK
    LOGICAL :: RK_fail

    INTEGER :: i_sum

    REAL*8 :: a_RK(6,6)
    REAL*8 :: b_RK(6)
    REAL*8 :: c_RK(6)

    REAL*8 :: b4_RK(6)

    REAL*8 :: dr_RK

    REAL*8 :: err_check

    ! define the coefficients for the RK 4-5 order scheme

    a_RK = 0.d0

    a_RK(2,1) = 1.d0 / 5.d0

    a_RK(3,1) = 3.d0 / 40.d0
    a_RK(3,2) = 9.d0 / 40.d0

    a_RK(4,1) = 3.d0 / 10.d0
    a_RK(4,2) = -9.d0 / 10.d0
    a_RK(4,3) = 6.d0 / 5.d0

    a_RK(5,1) = -11.d0 / 54.d0
    a_RK(5,2) = 5.d0 / 2.d0
    a_RK(5,3) = -70.d0 / 27.d0
    a_RK(5,4) = 35.d0 / 27.d0

    a_RK(6,1) = 1631.d0 / 55296.d0
    a_RK(6,2) = 175.d0 / 512.d0
    a_RK(6,3) = 575.d0 / 13824.d0
    a_RK(6,4) = 44275.d0 / 110592.d0
    a_RK(6,5) = 253.d0 / 4096.d0

    c_RK = 0.d0

    c_RK(2) = 1.D0 / 5.d0
    c_RK(3) = 3.d0 / 10.d0
    c_RK(4) = 3.d0 / 5.d0
    c_RK(5) = 1.d0
    c_RK(6) = 7.d0 / 8.d0

    ! coefficients for the 5th order RK

    b_RK = 0.d0

    b_RK(1) = 37.d0 / 378.d0
    b_RK(2) = 0.d0
    b_RK(3) = 250.d0 / 621.d0
    b_RK(4) = 125.d0 / 594.d0
    b_RK(5) = 0.d0
    b_RK(6) = 512.d0 / 1771.d0

    ! coefficient for the 4th order

    b4_RK = 0.d0

    b4_RK(1) = 2825.d0 / 27648.d0
    b4_RK(2) = 0.d0
    b4_RK(3) = 18575.d0 / 48384.d0
    b4_RK(4) = 13525.d0 / 55296.d0
    b4_RK(5) = 277.d0 / 14336.d0
    b4_RK(6) = 1.d0 / 4.d0
 
    !
    ! ... Set initial conditions at initial radius
    !
    CALL initialize_current
    !CALL initialize_particles
    !
    ! ... Get meteorological variables
    CALL gasvals
    CALL zmet
    CALL compute_mixture
    !
    ! ... Lump physical variables 
    CALL lump(f)         
    !
    ! ----------------------------------------------------
    ! ... assign initial stepping length and start marching loop
    ! ----------------------------------------------------
    dr = dr0
    dr_old = 0

    CALL write_output

    f_oldold = f

    main_loop: DO

       ! move one step back to get the terms needed to evaluate dbeta_dr and dh_dr:
       ! dbeta_dr = ( beta_new - beta_old ) / ( r_new - r_old)
       ! dh_dr = ( h_new - h_old ) / ( h_new - h_old )
       !
       ! These terms are evaluated in the subroutine rate

       r = r - dr_old
       CALL unlump(f_oldold)

       h_old = h
       r_old = r
       beta_old = beta 
       Ri_old = Ri

       f_stepold = f

       r = r + dr_old

       CALL unlump(f)        

       h_new = h
       r_new = r
       beta_new = beta 
       Ri_new = ri

       !-------------- evaluate the terms for the RK scheme --------------------!

       ! first step
       CALL zmet
       CALL rate( rhs_RK(1,:) )

       h_old = h_new
       r_old = r_new
       beta_old = beta_new 
       Ri_old = Ri_new

       RK_fail = .FALSE.
       
       RK_loop: DO i_RK=2,n_RK

          rhs(:) = 0.D0

          DO i_sum=1,i_RK-1

             rhs = rhs + a_RK(i_RK,i_sum) * rhs_RK(i_sum,:)
             
          END DO

          CALL marching( f_stepold , f_RK , rhs , dr )

          dr_RK = c_RK(i_RK)*dr
          
          r = r_old + dr_RK

          CALL unlump(f_RK)

          IF ( ( beta < alpha) .OR. ( h < 0 ) .OR. ( u < 0 ) .OR.               &
               ( MINVAL(f) < 0 ) ) THEN

             RK_fail = .TRUE.
             
             exit RK_LOOP

          END IF
    
          ! new values to evaluate the terms dh_dr and d_beta_dr
          h_new = h 
          r_new = r 
          beta_new = beta

          CALL zmet
          CALL rate( rhs_RK(i_RK,:) )
       
          
       END DO RK_loop

       rhs(:) = 0.D0
       rhs4(:) = 0.D0

       DO i_RK=1,n_RK
          
          rhs = rhs + b_RK(i_RK) * rhs_RK(i_RK,:)
          rhs4 = rhs4 + b4_RK(i_RK) * rhs_RK(i_RK,:)
          
       END DO

       ! update the solution with the 5th order RK scheme
       CALL marching(f_stepold,f,rhs,dr)
       
       r = r_old + dr

       CALL unlump(f)

       ! check if the solution cannot be accepted
       IF ( ( beta < alpha) .OR. ( h < 0 ) .OR. ( u < 0 ) .OR.                  &
            ( MINVAL(f) < 0 ) ) THEN

          RK_fail = .TRUE.
          
       END IF

       IF ( RK_fail ) THEN

          ! repeat the RK integration with a smaller integration step
          r = r_old          
          dr = 0.5 * dr
          f = f_stepold
          
       ELSE

          ! accept the solution and advance
          f_oldold = f_stepold

          dr_old = dr

          ! che the error 4th-5th and compute the adaptive step
          f_err(:) = dr * ( rhs(:) - rhs4(:) )
          
          err_check = MAXVAL( ABS(f_err) / ( eps_abs + eps_rel * ABS(f) ) )
          
          dr = MIN( dr0 , dr/err_check**0.2D0 )
          
          CALL write_output
    
          IF ( Ri .gt. 1.d0 .and. Ri_old .lt. 1.d0) THEN
       
          WRITE(*,*) 'hydraulicjump'
       
                flow_regime = 3
         
             CALL write_dakota
             EXIT main_loop
       
          ELSE
       
          IF ( ri_new < 1.d0 ) THEN

               !WRITE(*,*) 'supercritical'
          
               flow_regime = 1
          
          ELSE
          
               !WRITE(*,*) 'subcritical'
          
               flow_regime = 2
          
           END IF
          
       END IF

          IF ( ( (beta - alpha) < 1.D-5) .OR. ( dr .LE. 1.D-11 ) ) THEN

             ! exit condition
             CALL write_dakota
             
             EXIT main_loop
             
          END IF

       END IF
    
    END DO main_loop

    ! ----------------------------------------------------
    ! ... End of marching loop
    ! ----------------------------------------------------
    !

    RETURN

  END SUBROUTINE PDC

END MODULE PDCflow
!----------------------------------------------------------------------
