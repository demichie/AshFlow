!********************************************************************************
!> \brief environment module
!> This module contains all the variables related to the background environmental 
!> conditions, e.g. atmospheric conditions where h = 0 and the topography 
!> descriptors-
!********************************************************************************

MODULE envi_module

  IMPLICIT NONE

  ! Topography descriptors:
  REAL*8 :: theta!< Angle of slope   
  REAL*8 :: fric !< Friction param. describing relation between flow and ground 

  ! Atmospheric Conditions:
  REAL*8 :: gi		!< Gravitational acceleration 
  REAL*8 :: alpha       	!< Atmospheric density
  REAL*8 :: T_a      	!< Atmospheric temperature
  REAL*8 :: p	!< Atmospheric pressure        
  !
  SAVE

CONTAINS

!******************************************************************************
!> \brief Meteorological conditions
!> This subroutine is used to calculate atmospheric density, alpha. 
!******************************************************************************

  SUBROUTINE zmet

    ! External variables
    USE gas_module, ONLY: gas_constair
    USE current_module, ONLY: MARS_ATM
       
    IMPLICIT NONE

    IF (MARS_ATM) THEN

       !gi = 3.711D0
       !T_a = 237.D0
       !p = 636.D0
       !alpha = 0.0168D0
       
! Ancient Mars
       gi = 3.711D0
       T_a = 288.D0
       p = 20000.D0      
       alpha = 0.4D0
       
    ELSE   

! Earth conditions
       gi = 9.81D0
       T_a = 288.D0
       p = 100000.D0
       alpha = p/(gas_constair*T_a) !< Calculating atmospheric density
       
       !WRITE(*,*) 'gi, ta, p, alpha, gas_constair',        gi, T_a, p, alpha, gas_constair
       !READ(*,*)
       
    END IF

    RETURN

  END SUBROUTINE zmet

END MODULE envi_module

