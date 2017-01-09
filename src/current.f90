!********************************************************************************
!> \brief Current module
!> This module contains descriptors for the initial conditions of the flow, and parameters
!> required to evaluate these descriptors as the flow propagates.
!********************************************************************************

MODULE current_module !< provides source conditions for flow
  !       
  IMPLICIT NONE
  !
    
  !> flag for 1d or radial model
  LOGICAL :: oneD_model

  !> Flow thickness (m) at a given radius (r;m) from source
  REAL*8 :: h 
  
  !> Radius (m), describing distance from source 
  REAL*8 :: r
  
  !> Flow velocity (m/s) at a given radius (r;m) from source
  REAL*8 :: u 
  
  !> Temperature (K) of the flow at a given radius (r;m) from source 
  REAL*8 :: T
  
  !> Parameter used to evaluate flow radius integration step; defines radius at previous integration
  REAL*8 :: r_old
  
  !> Parameter used to evaluate flow radius integration step; defines radius at new integration location
  REAL*8 :: r_new
  
  !> Parameter used to evaluate flow thickness for integration; defines thickness at previous integration 
  REAL*8 :: h_old
  
  !> Parameter used to evaluate flow thickness for integration; defines thickness at new integration
  REAL*8 :: h_new
  
  !> Initial thickness of flow (m) provided in input file 
  REAL*8 :: h0      				
  
  !> Initial flow radius (m) provided in input file  
  REAL*8 :: r0   				!< initial radius of flow  
        
  !> Mass mass flux of solids in the flow
  REAL*8 :: solid_mass_flux 
  
  !> Initial flow temperature (K) provided in input file 
  REAL*8 :: T0
  
  !> Logical flag describing how the initial temperature is determined, true means given in input, if false, calculated in mixture
  LOGICAL :: TVENT_FLAG 

  !> Mass flux
  REAL*8 :: mass_flux
  
  !> Initial velocity for coignimbrite plume calculated using current mass flux
  REAL*8 :: plume_velocity
  
  !> Initial Mass Flux
  Real*8 :: Initial_MF
  
  !> Ratio  of Final Mass Flux to Initial Mass Flux
  Real*8 ::solidMF_ratio
    
  !
  SAVE
  !----------------------------------------------------------------------
CONTAINS

!******************************************************************************
!> \brief Initialising current
!> This subroutine allocates conditions at flow source to those assigned in input file
!******************************************************************************

  SUBROUTINE initialize_current
    IMPLICIT NONE
    !     
    h = h0
    r = r0
    !
    RETURN
  END SUBROUTINE initialize_current
  !----------------------------------------------------------------------
END MODULE current_module
!----------------------------------------------------------------------
