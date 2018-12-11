!********************************************************************
!> \brief Input/Output module
!
!> This module contains all the input/output subroutines and
!> related variables.
!********************************************************************

MODULE inpout

  IMPLICIT NONE

  CHARACTER(LEN=15) :: finp
  CHARACTER*150 :: fout5    !< Name of output file for backup of input parameters

  CHARACTER(LEN=80) :: run_name, input_name

  CHARACTER(LEN=80) :: out_file
  CHARACTER(LEN=80) :: out_file2
  CHARACTER(LEN=80) :: out_file3
  CHARACTER(LEN=80) :: out_file4


  CHARACTER(LEN=3) :: outtime               !< Output units for time
  REAL*8 :: factime
  CHARACTER(LEN=3) :: outpres_b = 'Pa '     !< Output units for pressure
  REAL*8 :: facpres
  CHARACTER(LEN=3) :: outflux = 'd00'       !< Output units for mass flux
  REAL*8 :: facflux
  CHARACTER(LEN=3) :: outmass = 'Kg '       !< Output units for remaining chamber mass
  REAL*8 :: facmass
  CHARACTER(LEN=3) :: outvelo = 'ms '       !< Output units for velocity
  REAL*8 :: facvelo
  CHARACTER(LEN=3) :: outdist = ' m '       !< Output units for distance
  REAL*8 :: facdist

  CHARACTER(LEN=5) :: outdens = 'Kgm-3'     !< Output units for density
  CHARACTER(LEN=3) :: outpres_t = 'Pa'      !< Output units for pressure

  INTEGER*4 :: iprint                       !< Flag for print

  REAL*8 :: pi
  !
  !**  Input files unit numbers
  !
  INTEGER, PARAMETER :: ninp = 7            !< Input data unit
  !
  !**  Output files unit numbers
  !
  INTEGER, PARAMETER :: nout5=14            !< Backup input unit

  SAVE

CONTAINS


  !*****************************************************************************
  !> \brief Initialize variables
  !
  !> This subroutine initialises the input variables
  !> \date 01/07/2014
  !> @author
  !> Sam & Mattia
  !******************************************************************************

  SUBROUTINE initialize

    USE current_module, ONLY: r

    USE mixture_module, ONLY: epsilon

    USE envi_module, ONLY: gi

    USE solver_module, ONLY: dh_dr

    IMPLICIT NONE

    outtime = 'sec'           ! Output units for time
    factime = 1.0d0
    outpres_b = 'Pa '         ! Output units for pressure
    facpres = 1.0d0
    outflux = 'd00'           ! Output units for mass flux
    facflux = 1.0d0
    outmass = 'Kg '           ! Output units for remaining chamber mass
    facmass = 1.0d0
    outvelo = 'ms '           ! Output units for velocity
    facvelo = 1.0d0
    outdist = ' m '           ! Output units for distance
    facdist = 1.0d0
    outpres_t = 'Pa '         ! Output units for pressure
    outdens = 'Kgm-3 '           ! Output units for density

  END SUBROUTINE initialize

  !*****************************************************************************
  !> \brief Read Input data
  !> This subroutine reads input data from the file flow_model.inp
  !> \date 01/07/2014
  !> @author
  !> Sam & Mattia
  !******************************************************************************

  SUBROUTINE reainp

    USE current_module, ONLY: oneD_model , h0, r0, TVENT_FLAG, T0

    USE envi_module, ONLY: fric, theta, gi, T_a, p, gas_constair, C_vair

    USE solver_module, ONLY: vel_equation , dr0 , eps_rel , eps_abs

    USE mixture_module, ONLY: nmag, N, Ri, cpwvapour, rwvapour

    USE particles_module, ONLY: iclass, diam, rhosol, fracsolid, C_d , mean_phi,&
         sigma_phi , min_phi , max_phi , C_d_input , phi_dist

    USE particles_module, ONLY: diam1 , diam2 , min_rho , max_rho

    USE tools, ONLY : runend

    USE particles_module, ONLY : allocate_particles , initialize_particles

    IMPLICIT NONE

    LOGICAL :: tend1
    CHARACTER(LEN=80) :: card

    INTEGER :: i

    LOGICAL :: eval_gsd_flag

    NAMELIST / control_parameters / run_name

    NAMELIST / flow_parameters /  oneD_model , vel_equation , dr0 , eps_rel , eps_abs

    NAMELIST / env_parameters /  fric, theta, gi, T_a, P, gas_constair, c_vair, rwvapour, cpwvapour

    NAMELIST / initial_values / r0, h0, Ri, TVENT_FLAG, T0, n, nmag , eval_gsd_flag

    NAMELIST / density_parameters / diam1 , diam2 , min_rho , max_rho

    finp = 'flow_model.inp'

    OPEN(ninp,FILE=finp,STATUS='old')

    READ(ninp, control_parameters)

    READ(ninp, flow_parameters)

    READ(ninp, env_parameters)

    pi = 4.D0 * ATAN(1.D0)
    theta = theta*pi/180

    READ(ninp, initial_values)

    READ(ninp, density_parameters)

    tend1 = .FALSE.

    granulometry_search: DO

       READ( ninp , *, END = 300 ) card

       IF( TRIM(card) == 'GRANULOMETRY' ) THEN

          EXIT granulometry_search

       END IF

    END DO granulometry_search

    READ(ninp,*) iclass

    IF (iclass .LT. 1) THEN

    	WRITE(*,*) 'INPUT_ASH_FLOW_MODEL: no particle classes'

    ELSE

       CALL  allocate_particles

    ENDIF

    IF ( eval_gsd_flag ) THEN

       READ(ninp,*) mean_phi, sigma_phi , min_phi , max_phi , C_d_input

       CALL initialize_particles

    ELSE

       DO i = 1, iclass

          READ(ninp,*) diam(i), rhosol(i), fracsolid(i), C_d(i)

       ENDDO

    END IF

    GOTO 310
300 tend1 = .TRUE.
310 CONTINUE

    !
    !***  Close input file
    !

    CLOSE(ninp)

    fout5 = TRIM(run_name)//'.bak'

    OPEN(nout5,file=fout5,status='unknown')

    WRITE(nout5, control_parameters)

    WRITE(nout5, flow_parameters)

    WRITE(nout5, env_parameters)

    WRITE(nout5, initial_values)

    WRITE(nout5, density_parameters)

    IF (( tend1) .OR. (iclass .EQ. 0)) THEN

       WRITE(*,*) 'WARNING: input', 'SAMPLING POINTS not found'

    ELSE

       WRITE(nout5,*) '''GRANULOMETRY'''
       WRITE(nout5,*) iclass

       IF ( eval_gsd_flag ) THEN

          WRITE(nout5,*) mean_phi, sigma_phi , min_phi , max_phi , C_d_input

       ELSE

          DO i = 1, iclass

             WRITE(nout5,*) diam(i), rhosol(i), fracsolid(i), C_d(i)

          END DO

       END IF

    END IF

    CLOSE(nout5)

    out_file = TRIM(run_name)//'.flow'
    out_file2 = TRIM(run_name)//'.part'
    out_file3 = TRIM(run_name)//'.dep'
    out_file4 = TRIM(run_name)//'.ep'

    !
    OPEN(51,FILE=out_file)
    OPEN(52,FILE=out_file2)
    OPEN(53,FILE='dakota_run.dak')
    OPEN(54,FILE=out_file3)
    OPEN(56, FILE=out_file4)

    IF ( eval_gsd_flag ) THEN

       WRITE(52,177) iclass,phi_dist(1:iclass)

    END IF

177 FORMAT(i5,1x,50(e12.5, 1x))

    RETURN

1007 CALL runend(-1, 'Error: cannot OPEN file: '//fout5)

  END SUBROUTINE reainp

  SUBROUTINE write_output

    USE current_module, ONLY : r , h , T , solid_mass_flux , u , oneD_model, mass_flux, plume_velocity, Initial_MF,  solidMF_ratio
    USE particles_module, ONLY : sumsed, fracsolid, v_s, iclass, final, Acc_rate, diam, norm_flux, ipf, fpf, rhosol, diam
    USE particles_module, ONLY : wavelengthfreen, wavelengthshearn, wavelengthfreex, wavelengthshearx, wavelength_height
    USE particles_module, ONLY : Pn, k, ustar, H_bl, ks, eta0, eta, SS0, BV, init_fracsolid, rhosol, diam, phi_dist, no_bins
    USE mixture_module, ONLY : Ri , epsilon , n , beta, Pdyn, final_solid_mass_flux
    USE mixture_module, ONLY : init_N, initial_SMF, entrainment_rate
    USE envi_module, ONLY : gi

    IMPLICIT NONE

    REAL*8 :: volume_flux
    REAL*8 :: mean_grainsize
    REAL*8 :: std_dev

    INTEGER :: i
    INTEGER :: j

    IF ( oneD_model ) THEN

       volume_flux = h * u

    ELSE

       volume_flux = h*u*r

       mass_flux = (beta * u * h * r) * 2 * pi
       final_solid_mass_flux = (beta * u * h * r) * ( 1.D0 - n )
       plume_velocity = mass_flux/(3.14D0 * r**2 * beta)
       solidMF_ratio = final_solid_mass_flux/ initial_SMF

	mean_grainsize = SUM( phi_dist(1:no_bins)*fracsolid(1:no_bins))
	std_dev = sqrt(SUM( fracsolid(1:no_bins)*(phi_dist(1:no_bins)- mean_grainsize)**2 ) )

	!WRITE(*,*) 'mean_grainsize', mean_grainsize
	!WRITE(*,*) 'std_dev', std_dev

       ! Particle Mass Flux
       !DO i=1,iclass

       ! initial particle flux

       !ipf = Initial_MF * (1-init_N) * init_fracsolid
       !fpf = mass_flux * (1-n) * fracsolid
       !norm_flux = (ipf - fpf) / ipf

       !WRITE (*,*) 'ipf, fpf, norm_flux', ipf, fpf, norm_flux

       !ENDDO


    END IF

    !
    ! *** equations for postprocessing: flow dynamic pressure and deposit wavelength
    !

    ! Dynamic Pressure (kPa)

    Pdyn = (0.5 * beta * u **2)* 0.001

    ! Boundary Layer Thickness -- taken from Amanda's code from Allen 1970 book

    H_bl = 0.376D0*r*(0.0001D0/(beta*u*r))**(0.2)

    ! Eq. 14 Valentine 1987
    ! ustar = (u * k) / LOG(30 * H_bl /ks)

    DO i=1,iclass

    	ustar = (u * k) / LOG(50 * H_bl /diam(i))

 		! Accumulation rate

 		 Acc_rate = (v_s(i) *  solid_mass_flux)/( h * u) !Accumulation rate kg/s

        ! Rouse Number -- eq 2 from Valentine (1987)

         Pn = v_s(i) / (0.4D0 * ustar)

        ! Non dimension reference height

         eta0 = H_bl/h

        ! Brunt-Vaisala Frequency

    		DO j = 1,final

        		eta(j) =  0.1D0 * (j)                                                 ! Non-dimensionless height in the flow, between 0 and 1

        		SS0(j) = (eta0 / (1 - eta0)) * ((1 - eta(j)) / eta(j)) ** Pn			! Non-dimensionalised concentration profile: Equation 4 from Valentine 1987

        	    BV(j) = (1/(2 * pi)) * (( gi / h ) * ( Pn / (eta(j) *(1-eta(j))))) ** 0.5   !Brunt-Vaisala Frequency as a function on non-dimensional height in the flow, equation 8 in Valentine 1987


    		ENDDO


    	 wavelengthfreen = u / MAXVAL(BV)	                             ! minimum wavelength using freestream velocity
   		 wavelengthshearn = ustar / MAXVAL(BV)                            !  minimum  wavelength using shear velocity

    	 wavelengthfreex = u / MINVAL(BV)                             ! maximum wavelength using freestream velocity
    	 wavelengthshearx = ustar / MINVAL(BV)                            !  maximum  wavelength using shear velocity

    	 wavelength_height = 2 * pi * h * 0.4                             ! Using Yalin, 1972 wavelength = 2*pi*h

       			!WRITE(*,*) 'wavelengthfreen', wavelengthfreen
    			!WRITE(*,*) 'wavelengthshearn', wavelengthshearn
    			!WRITE(*,*) 'wavelengthfreex', wavelengthfreex
    			!WRITE(*,*) 'wavelengthshearx', wavelengthshearx
    			!WRITE(*,*) 'wavelength_height', wavelength_height
    			!READ(*,*)

    ENDDO

    ! Boundary Layer Thickness

    WRITE(51,100) r , h , u , T , beta , solid_mass_flux , sumsed , Ri ,        &
         epsilon , volume_flux , n, Pdyn

    !WRITE(54,155) r , mean_grainsize, std_dev

    WRITE(56,156) r , entrainment_rate

    WRITE(54,155) r , h , u , v_s , beta, Pn, wavelengthfreen,		&
 			wavelength_height, Acc_rate

    WRITE(52,177) r , fracsolid


100 FORMAT(f12.5,1x,e12.5,1x,e12.5,1x,e12.5,1x,e12.5,1x,e12.5,1x,e12.5,         &
         1x,e12.5,1x,e12.5, 1x,e12.5, 1x,e12.5, 1x,e12.5)

!155 FORMAT(f12.5,1x,e12.5,1x,e12.5,1x)

155 FORMAT(f12.5,1x,e12.5,1x,e12.5,1x,e12.5,1x,e12.5,1x,e12.5,1x,e12.5,         &
         1x,e12.5,1x,e12.5)

156 FORMAT(f12.5,1x,e12.5)

177 FORMAT(f12.5,1x,50(e12.5, 1x))


  END SUBROUTINE write_output

  SUBROUTINE write_dakota

    USE current_module, ONLY : r , h , T , solid_mass_flux , u , oneD_model, mass_flux, plume_velocity, Initial_MF, solidMF_ratio
    USE particles_module, ONLY : sumsed, fracsolid, phi_dist, no_bins, iclass, diam, rhosol, v_s, norm_flux
    USE mixture_module, ONLY : Ri , epsilon , n , beta, flow_regime, gas_constmix, C_vmix, initial_velocity
    USE mixture_module, ONLY : initial_velocity, initial_density, solidvolumefraction, final_solid_mass_flux, initial_SMF
    USE envi_module, ONLY : alpha, gi

    IMPLICIT NONE

	INTEGER :: i
	CHARACTER(len=8) :: x1

    WRITE(53,*) 'initial mass flux', LOG10(Initial_MF)
    WRITE(53,*) 'initial_velocity', initial_velocity
    WRITE(53,*) 'initial_density', initial_density
    WRITE(53,*) 'solidvolumefraction', solidvolumefraction
    WRITE(53,*) 'radius', r
    WRITE(53,*) 'temperature', t
    WRITE(53,*) 'Gas Mass Fraction', N
    WRITE(53,*) 'Sedimentation', sumsed
    WRITE(53,*) 'Flow Velocity', u
    WRITE(53,*) 'Flow Density', beta
    WRITE(53,*) 'Ri', ri
    WRITE(53,*) 'Flow Height', h
    WRITE(53,*) 'Flow Regime', flow_regime
    WRITE(53,*) 'Solid Mass Fraction', 1-N
    WRITE(53,*) 'Final Mass Flux', mass_flux
    WRITE(53,*) 'Initial Plume Velocity', plume_velocity
    WRITE(53,*) 'Gas Mass Constant', gas_constmix
    WRITE(53,*) 'C_vmix', C_vmix
    WRITE(53,*) 'deltabeta1', beta - alpha
    WRITE(53,*) 'deltabeta_g', (beta - alpha) * gi
    WRITE(53,*) 'settling velocity', v_s
    WRITE(53,*) 'final_smflux', final_solid_mass_flux
    WRITE(53,*) 'initial_smflux', initial_SMF
    WRITE(53,*) 'solidMF_ratio', solidMF_ratio
    WRITE(53,*) 'alpha', alpha


	 DO i = 1, iclass


     WRITE(x1,'(I2.2)')i
!     WRITE(53,*) 'Diam'//x1,diam(i)*10e6
     WRITE(53,*) 'Diam'//x1,diam(i)
     WRITE(53,*) 'Frac'//x1,fracsolid(i)
     !WRITE(53,*) 'norm_flux'//x1,norm_flux(i)
     !WRITE(53,*) 'Rhosol'//x1,rhosol(i)
     !WRITE(nout5,*) diam(i), rhosol(i), fracsolid(i), C_d(i)


     END DO



  END SUBROUTINE write_dakota


  SUBROUTINE close_units

    CLOSE(51)
    CLOSE(52)
    CLOSE(53)
    CLOSE(54)
    CLOSE(56)


  END SUBROUTINE close_units



END MODULE inpout
