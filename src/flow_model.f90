!********************************************************************************
!> \mainpage   FlowModel  
!> FlowModel is a numerical code for the steady-state ash flow model, describing 
!> the propagation of a mixture of entrained air and volcanic ash with distance 
!> from source during an eruption. The system of equations is formally the same 
!> of that presented in Bursik and Woods 1996, i.e. the equations we formulate 
!> describe the same conservation principles (conservation of mass, momentum and 
!> energy). In this model, the specific heat and gas constant of the flow are 
!> solved at each iteration.
!>
!> Version 1.0:\n
!> 
!> \authors Samantha Engwell (1,*) & Mattia de' Michieli Vitturi (1) \n
!> (1) Istituto Nazionale di Geofisica e vulcanologia, sezione di Pisa\n
!>     Via della Faggiola, 36\n
!>     I-56126 Pisa, Italy \n
!> (*) E-mail: samantha.engwell@pi.ingv.it \n
!********************************************************************************

!> \brief Main Program 


PROGRAM flow_model

  USE inpout, ONLY: initialize , reainp
  USE inpout, ONLY: run_name , close_units

  USE PDCflow, ONLY: pdc
  USE solver_module, ONLY: allocate_matrix

  IMPLICIT NONE
  !
  !***  Initialize the input variables
  !
  CALL initialize

  !
  !***  Open and Read input data from file cpiuc.inp
  !
  CALL reainp

  !
  !***  Allocate variables for the colum model
  !
  CALL allocate_matrix

  !***  Solve the model
  CALL pdc(run_name)

  CALL close_units

END PROGRAM flow_model
