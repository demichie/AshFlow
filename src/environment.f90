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
  REAL*8 :: gas_constair !< perfect gas constant of ambient air ( J/(kg K) )
  REAL*8 :: C_vair !< specific heat ambient air ( J/(kg K) )       
  !
  SAVE

CONTAINS

!******************************************************************************
!> \brief Meteorological conditions
!> This subroutine is used to calculate atmospheric density, alpha. 
!******************************************************************************

  SUBROUTINE zmet
       
    IMPLICIT NONE
    
    alpha = p/(gas_constair*T_a)
    
    RETURN

  END SUBROUTINE zmet

END MODULE envi_module

