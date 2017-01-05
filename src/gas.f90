!********************************************************************************
!> \brief gas module
!> This module contains values for gas constant and specific heat of ambient air 
!> and water vapour.
!********************************************************************************
!
!MODULE gas_module
!  !
!  IMPLICIT NONE
!  !> perfect gas constant of ambient air ( J/(kg K) )
!  REAL, PARAMETER :: gas_constair = 285.D0
!
!  !> perfect gas constant for water vapour ( J/(kg K) )
!  REAL*8, PARAMETER :: rwvapour=461.D0
!
!  !> specific heat ambient air ( J/(kg K) )
!  REAL, PARAMETER :: C_vair=998.D0
!
!  !> specific heat of water vapour ( J/(kg K) )
!  REAL, PARAMETER :: cpwvapour=1617.D0
!  
!  SAVE
!
!END MODULE gas_module

MODULE gas_module

   IMPLICIT NONE

   !> perfect gas constant of ambient air ( J/(kg K) )
   REAL*8 :: gas_constair
       
   !> specific heat ambient air ( J/(kg K) )
   REAL*8 :: C_vair
       
   !> perfect gas constant for water vapour ( J/(kg K) )
   REAL*8 :: rwvapour
       
   !> specific heat of water vapour ( J/(kg K) )
   REAL*8 :: cpwvapour

   SAVE 
CONTAINS

!******************************************************************************
!> \brief Atmospheric consitions
!> This subroutine is used to calculate assign the gas conditions. 
!*****************************************************************************
      
 SUBROUTINE gasvals 
 
 USE current_module, ONLY: MARS_ATM
 
 IMPLICIT NONE
    
      
      IF (MARS_ATM) THEN
      
      	gas_constair = 191.D0
      	C_vair=734.9D0
      	rwvapour=300.D0
      	cpwvapour=2425.D0

      ELSE   
     
      	gas_constair = 287.D0
      	C_vair=998.D0
      	rwvapour=461.D0
      	cpwvapour=1617.D0
      	
      	!WRITE(*,*) 'gas_constair, C_vair, rwvapour, cpwvapour', gas_constair, C_vair, rwvapour, cpwvapour
      	!READ(*,*)
      	

      END IF

   RETURN

 END SUBROUTINE gasvals
 
END MODULE gas_module

