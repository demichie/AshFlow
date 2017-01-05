!********************************************************************************
!> \brief Solver Module
!> This module contains the differential equations that are solved at each 
!> iteration.
!> \date 01/07/2014
!> @author 
!> Sam & Mattia
!********************************************************************************

MODULE solver_module
  USE particles_module, ONLY: iclass
  USE current_module
  USE mixture_module
  USE envi_module, ONLY: p, T_a, gi
  USE gas_module, ONLY: gas_constair, C_vair
  !
  IMPLICIT NONE
  !
  ! ... Right-Hand Side (rhs) and Known Term (f) of the equation system
  !
  REAL*8, ALLOCATABLE, DIMENSION(:) :: rhstemp
  REAL*8, ALLOCATABLE, DIMENSION(:) :: rhs
  REAL*8, ALLOCATABLE, DIMENSION(:) :: ftemp
  REAL*8, ALLOCATABLE, DIMENSION(:) :: f
  REAL*8, ALLOCATABLE, DIMENSION(:) :: f_stepold
  REAL*8, ALLOCATABLE, DIMENSION(:) :: f_oldold

  REAL*8, ALLOCATABLE, DIMENSION(:) :: f_new
  REAL*8, ALLOCATABLE, DIMENSION(:) :: f_RK
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: rhs_RK
  REAL*8, ALLOCATABLE, DIMENSION(:) :: rhs4
  REAL*8, ALLOCATABLE, DIMENSION(:) :: f_err

  REAL*8 :: eps_rel , eps_abs

  INTEGER :: itotal

  INTEGER :: n_RK                 !< steps for the RK scheme

  REAL*8 :: dh_dr                 !< Change in flow height with radius
  REAL*8 :: dbeta_dr              !< Change in flow density with radius
  REAL*8 :: dr                    !< Change in radius
  REAL*8 :: dr0                   !< Change in radius
  REAL*8 :: dr_old                !< Old radius - used for backwards forwards analysis
  !
  LOGICAL :: vel_equation

  SAVE

CONTAINS

  !******************************************************************************
  !> \brief Allocate matrix
  !> This subroutine defines the matrix over which each of the conditions are 
  !> solved.
  !> \date 01/07/2014
  !> @author 
  !> Sam & Mattia
  !******************************************************************************

  SUBROUTINE allocate_matrix

    IMPLICIT NONE

    itotal = iclass + 5

    n_RK = 6

    ALLOCATE(rhstemp(itotal))
    ALLOCATE(rhs(itotal))
    ALLOCATE(ftemp(itotal)) 
    ALLOCATE(f(itotal)) 
    ALLOCATE(f_stepold(itotal))   

    ALLOCATE(f_new(itotal))   
    ALLOCATE(f_RK(itotal))   
    ALLOCATE(rhs_RK(n_RK,itotal))
    ALLOCATE(rhs4(itotal))
    ALLOCATE(f_err(itotal))

    f = 0.D0
    ftemp = 0.D0
    rhs = 0.D0
    rhstemp = 0.D0

    RETURN
  END SUBROUTINE allocate_matrix


  !*****************************************************************************
  !> \brief Right hand side equations
  !> This subroutine contains the right-hand side of the equations 
  !> that are computed for each particle class
  !> \param[out]    rhs_     right-hand side
  !> \date 01/07/2014
  !> @author 
  !> Sam & Mattia
  !*****************************************************************************

  SUBROUTINE rate(rhs_)

    ! External variables
    USE envi_module, ONLY: alpha, theta, fric
    USE particles_module, ONLY: rhosol , diam, solidmassflux_fract, v_s, S,     &
         sumsed, fracsolid, C_d

    IMPLICIT NONE

    REAL*8, DIMENSION(:), INTENT(OUT) :: rhs_

    REAL*8 :: var , coeff_radial
    
    REAL*8 :: term1 , term2

    REAL*8 :: rhs_2A , rhs_2B


    LOGICAL :: no_entr , const_volume_flux

    INTEGER :: i

    IF ( oneD_model ) THEN

       var = 1.D0

    ELSE

       var = r

    END IF

    solid_mass_flux = (beta * u * h * var) * ( 1.D0 - n )  
    ! provides the mass_flux of a particle of a given size
    solidmassflux_fract(1:iclass) = solid_mass_flux*fracsolid(1:iclass)  

    ! SEDIMENTATION for the particle fraction
    v_s(1:iclass) = DSQRT( (rhosol(1:iclass) * gi * diam(1:iclass)) /           &
         (C_d(1:iclass) * beta) )

    S(1:iclass) = 1.D0 / (h * u) * solidmassflux_fract(1:iclass) * v_s(1:iclass)

    sumsed = SUM(S(1:iclass)) ! sum of sedimentation across all particles sizes

    !---- Mass conservation of the mixture   (Eq.1 Bursik & Woods 1996)

    entrainment_rate = epsilon * alpha * u * var
    rhs_(1) = entrainment_rate - sumsed
    
    !---- Momentum Conservation (Eq.4 Bursik & Woods 1996)

    IF ( r_new .NE. r_old ) THEN

       dbeta_dr = (beta_new - beta_old) / (r_new - r_old)
       dh_dr = (h_new - h_old) / (r_new - r_old)

    ELSE

       dbeta_dr = 0.D0
       dh_dr = 0.D0

    END IF

    coeff_radial = ( var - 1.D0 ) / ( r - 1.D0 )

    no_entr = .FALSE.
    
    const_volume_flux = .FALSE.

    IF ( no_entr ) THEN

       term1 = 0.D0

    ELSE
       
       term1 = - u * epsilon * alpha / ( h * beta )

    END IF

    IF ( const_volume_flux ) THEN

       term2 = 0.D0

    ELSE

       term2 = - ( epsilon * alpha * u * var - sumsed ) / ( beta * h * var ) *  &
            dcos(theta) + u / beta * dcos(theta) * dbeta_dr

    END IF

    ! RHS term for the modified momentum equation 
    rhs_2A = ( term1 + Ri * ( term2 + coeff_radial * u / var * dcos(theta) +    &
         u / h * dsin(theta) - u / ( 2.D0 * ( beta - alpha ) ) * dcos(theta) *  &
         dbeta_dr ) - fric * u / h ) / ( 1.D0 - Ri * cos(theta) )

    ! RHS term for the original equation
    rhs_2B =  - epsilon * alpha * u / ( beta * h ) + Ri * u / h * ( - dh_dr     &
         * dcos(theta) + dsin(theta) ) - Ri * u  / ( 2.D0 * ( beta - alpha ))   &
         * dcos(theta) * dbeta_dr - fric * u / h

    IF ( vel_equation ) THEN

       rhs_(2) = rhs_2A
       
    ELSE

       rhs_(2) = rhs_2B

    END IF

    !---- Thermal Energy Conservation (Eq.5 Bursik & Woods 1996)

    rhs_(3) = epsilon * alpha * u * var * (C_vair * T_a + p/alpha) - sumsed *   &
         C_s * T    

    !---- Mass average specific heat
    rhs_(4) = ((C_vair - C_vmix) / beta * epsilon * alpha * u * var + (C_s -    &
         C_vmix) / beta * (-sumsed))/ (u * h * var)    

    !---- Gas constant of mixture
    rhs_(5) = (gas_constair - gas_constmix) / (beta * n * u * h * var ) *       &
         epsilon * alpha * u * var

    DO i=1, iclass

       rhs_(5+i) = -S(i)

    ENDDO

    RETURN
  END SUBROUTINE rate


  !*****************************************************************************
  !> \brief Lump variables
  !> Variables are grouped together to calculate the lumped variables 
  !> \param[out]    f_     lumped variables
  !> \date 01/07/2014
  !> @author 
  !> Sam & Mattia
  !*****************************************************************************

  SUBROUTINE lump(f_)

    ! External variables
    USE envi_module, ONLY: p
    USE gas_module, ONLY: cpwvapour
    USE mixture_module, ONLY: C_vmix
    USE particles_module, ONLY: fracsolid, iclass

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(OUT) :: f_
    INTEGER :: i
    REAL*8 :: var

    CALL compute_mixture
    !

    IF ( oneD_model ) THEN

       var = 1.D0

    ELSE

       var = r

    END IF


    f_(1) = beta * u * h * var

    f_(2) = u
    f_(3) = beta * u * h * var * (C_vmix * T + p/beta + 0.5 * u**2)
    f_(4) = C_vmix
    f_(5) = gas_constmix

    DO i=1, iclass  
     
       f_(5+i) = (beta * u * h * var) * (1-n) * fracsolid(i) 
   
    ENDDO

    RETURN
  END SUBROUTINE lump


  !******************************************************************************
  !> \brief Marching r one step
  !
  !> This subroutine update the solution of the model from r to r+dr as 
  !> fnew=fold+dr*rate.
  !> \param[in]    fold    old lumped variables
  !> \param[in]    rate     rate of change of the lumped variables
  !> \param[out]   fnew    new lumped variables
  !> \date 01/07/2014
  !> @author 
  !> Sam & Mattia
  !*****************************************************************************

  SUBROUTINE marching(fold,fnew,rat,delta_r)

    IMPLICIT NONE
    REAL*8, DIMENSION(:), INTENT(IN) :: fold, rat
    REAL*8, INTENT(IN) :: delta_r
    REAL*8, DIMENSION(:), INTENT(OUT) :: fnew
    INTEGER :: i

    DO i=1,itotal

       fnew(i) = fold(i) + rat(i) * delta_r

    END DO

    RETURN

  END SUBROUTINE marching


  !*****************************************************************************
  !> \brief Unlumping variable
  !> Lumped variables defined in subroutine 'lump' are unlumped to 
  !> calculate physical variables from lumped variables for each particle class
  !> \date 01/07/2014
  !> @author 
  !> Sam & Mattia
  !> \param[in]   f_   lumped variables 
  !*****************************************************************************

  SUBROUTINE unlump(f_)

    ! External varuables
    USE envi_module, ONLY: alpha , p
    USE mixture_module, ONLY: rhosol_ave
    USE particles_module, ONLY: iclass, fracsolid, rhosol

    IMPLICIT NONE

    REAL*8, DIMENSION(:), INTENT(IN) :: f_

    REAL*8 :: C1 , C2 , C3 , C4

    REAL*8 :: var

    IF ( oneD_model ) THEN

       var = 1.D0

    ELSE

       var = r

    END IF

    u = f_(2)
    n = 1.D0 - SUM(f_(5+1:5+iclass))/f_(1)    
    C_vmix = f_(4)    
    gas_constmix = f_(5)
    fracsolid(1:iclass) = (f_(5+1:5+iclass))/SUM(f_(5+1:5+iclass))

    rhosol_ave = 1.D0 / (SUM(fracsolid/rhosol))

    C1 = ( 1.D0 - n ) / rhosol_ave
    C2 = n * gas_constmix / p
    C3 = ( f_(3)/f_(1) - 0.5D0 * u**2 ) / C_vmix
    C4 = - p / C_vmix
    
    beta = ( 1.D0 - C2*C4 ) / ( C2*C3 + C1 )
    T = C4 / beta + C3
    
    !WRITE(*,*) 'rhogas', 1.0/(C2*T/n)
    !WRITE(*,*) 'n, T, gasconstmix, p',n, T, gas_constmix, p

    h = f_(1) / (beta * u * var)

    Ri = ((beta - alpha) * gi * h) / (beta * u**2.D0) 
    epsilon = 0.075D0 / DSQRT(1.D0 + (718.D0 * Ri**2.4D0))

    RETURN

  END SUBROUTINE unlump

END MODULE solver_module

