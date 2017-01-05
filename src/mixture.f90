!********************************************************************************
!> \brief Mixture module
!> This module contains all the variables required for describing and
!> calculating the characteristics of the mixture.
!> \date 01/07/2014
!> @author 
!> Sam & Mattia
!********************************************************************************

MODULE mixture_module

  USE envi_module, ONLY: alpha, T_a, gi

  IMPLICIT NONE

  !> Specific heat mixture
  REAL*8 :: C_vmix

  !> Entrainment coefficient
  REAL*8 :: epsilon

  !> Volatile mass fraction of erupted material
  REAL*8 :: Nmag

  !> Total gas mass fraction in the mixture, sum of entr. air and erupted vol.
  REAL*8 :: N
  
  !> Initial gas mass fraction
  REAL*8 :: init_N

  !> Mass fraction of entrained air in the ash flow 
  REAL*8 :: lambda

  !> Density of the air in the flow
  REAL*8 :: rhogas

  !> Density of the flow, used in iteration
  REAL*8 :: beta_old

  !> Density of the flow, used in iteration
  REAL*8 :: beta_new

  !> Richardson number of flow
  REAL*8 :: Ri
  
  !> Description of flow regime to identify whether sub, or supercritical flow, or whether a hydraulic jump occurs
  REAL*8 :: flow_regime

  !> Gas constant for the mixture
  REAL*8 :: gas_constmix

  !> Specific heat of solid material in the flow
  REAL*8 :: C_s = 1617.D0	

  !> Mixture density (kg/m3)
  REAL*8 :: beta

  !> Average density of solid phase (kg/m3)
  REAL*8 :: rhosol_ave
  
  REAL*8 :: pi
  
  REAL*8 :: initial_velocity
  
  REAL*8 :: initial_density
  
  REAL*8 :: Cv_magmix
  
  REAL*8 :: solidvolumefraction
  REAL*8 :: gasvolumefraction
  REAL*8 :: final_solid_mass_flux
  REAL*8 :: initial_SMF
  
  REAL*8 :: deltabeta
  
  REAL*8 :: deltabeta_g
  
  REAL*8 :: Pdyn
  
  REAL*8 :: entrainment_rate

  SAVE

CONTAINS

  !******************************************************************************
  !> \brief Computing conditions within mixture
  !> This subroutine calculates the characteristics of the mixture by combining 
  !> particle
  !> characteristics for each class with initial conditions
  !> \date 01/07/2014
  !> @author 
  !> Sam & Mattia
  !******************************************************************************

  SUBROUTINE compute_mixture
    !
    USE gas_module, ONLY: cpwvapour, C_vair, rwvapour, gas_constair
    USE envi_module, ONLY: alpha ,p
    USE current_module, ONLY: u, r, h, T, solid_mass_flux, T0, TVENT_FLAG ,     &
         oneD_model, Initial_MF
    USE particles_module, ONLY: rhosol, diam, iclass, fracsolid , v_s, S, C_d , &
         sumsed ,solidmassflux_fract, Acc_rate
    !      
    IMPLICIT NONE

    INTEGER :: i	
    
    REAL*8 :: var
    
    ! mass fraction of magmatic gas in the mixture (with air and solid)
    REAL*8 :: nmag_mix

    ! mass fraction pf magmatic gas in the gas (volcanic+air)
    REAL*8 :: nmag_gas
    
    pi = 4.D0 * ATAN(1.D0)
    
    init_N = N

    nmag_mix = nmag * ( 1.D0 - N ) / (1.D0 - nmag)

    lambda = (n - nmag) / ( 1.D0 - nmag )

    C_vmix = nmag_mix * cpwvapour + lambda * C_vair                   &
         + ( 1.D0 - n ) * C_s 
         
    Cv_magmix = nmag * cpwvapour + ( 1.D0 - nmag ) * C_s

    nmag_gas = nmag_mix / N

    gas_constmix = rwvapour * nmag_gas + gas_constair * ( 1.D0 - nmag_gas ) 
    
    !WRITE(*,*) 'nmag_mix, nmag_gas, lambda', nmag_mix, nmag_gas, lambda
    !WRITE(*,*) 'C_vmix, Cv_magmix, gas_constmix', C_vmix, Cv_magmix, gas_constmix
    !READ(*,*)
    
    ! Calculating the initial temperature using temperature of erupted material 
    ! if not explicitly defined in input file   

    IF (TVENT_FLAG) THEN

       ! Initial temperature in the flow based on the eruption temp (eq. 14 in 
       ! Bursik and Woods 96)    

       T = ( ( 1.D0 - lambda ) * Cv_magmix * T0 + lambda * C_vair * T_a ) /        &
            ( ( 1.D0 - lambda )* Cv_magmix + lambda * C_vair )
	
	!WRITE(*,*) 'T', T
	!READ(*,*)
	
    ELSE   

       T = T0

    END IF

    rhogas = p / ( gas_constmix * T )
    rhosol_ave = 1/(SUM(fracsolid/rhosol))
    beta = 1.D0 / ( n / rhogas + ( 1.D0 - n ) / rhosol_ave)
    
    !WRITE(*,*) 'Density', beta
    !READ(*,*)
    
    initial_density = beta
  
    ! Flow velocity
    u = SQRT((beta - alpha) * gi * h / (Ri * beta))
    initial_velocity = u

    ! Defining entrainment conditions          
     epsilon = 0.075D0 / DSQRT( 1.0D0 + 718.0D0 * Ri**2.4D0 )

    IF ( oneD_model ) THEN

       var = 1.D0

    ELSE

       var = r

    END IF

    solid_mass_flux = (beta * u * h * var) * ( 1.D0 - n ) 
    initial_SMF = (beta * u * h * var) * ( 1.D0 - n ) 
    
    Initial_MF =  (beta * u * h * var) * 2 * pi
 
    ! provides the mass_flux of a particle of a given size
    DO i=1,iclass

       solidmassflux_fract(i) = solid_mass_flux*fracsolid(i)  

       ! Settling velocity for the particle fraction

       v_s(i) = DSQRT( (rhosol(i) * gi * diam(i)) / (C_d(i) * beta) )  

       ! SEDIMENTATION for the particle fraction

       S(i) = 1.D0 / (h * u) * solidmassflux_fract(i) * v_s(i)

    ENDDO

    sumsed = SUM(S)

    RETURN
  END SUBROUTINE compute_mixture

  !----------------------------------------------------------------------
END MODULE mixture_module
!----------------------------------------------------------------------
