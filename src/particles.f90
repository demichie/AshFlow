!********************************************************************************
!> \brief Particles module
!> This module contains the procedures and the variables related to the solid
!> particles. 
!********************************************************************************
MODULE particles_module
  !
  IMPLICIT NONE

  !> The maximum particle diameter when defining grainsize distribution   
  REAL*8 :: min_phi

  !> The minimum particle diameter when defining grainsize distribution
  REAL*8 :: max_phi

  REAL*8 :: diam2

  !> The maximum particle density when defining grainsize distribution   
  REAL*8 :: max_rho

  REAL*8 :: diam1

  !> The minimum particle density when defining grainsize distribution
  REAL*8 :: min_rho

  !> Number of bins grainsize distribution should be split into
  INTEGER :: no_bins

  !> Standard deviation used to produce grainsize distribution
  REAL*8 :: sigma_phi

  !> Mean used to produce grainsize distribution
  REAL*8 :: mean_phi

  !> Mass fraction of the particle phases with respect to the total solid
  REAL*8, ALLOCATABLE, DIMENSION(:) :: fracsolid 
  
  !> Mass fraction of the particle phases with respect to the total solid
  REAL*8, ALLOCATABLE, DIMENSION(:) :: init_fracsolid 

  !> Density of each solid class 
  REAL*8, ALLOCATABLE, DIMENSION(:) :: rhosol

  !> Diameter of each solid class  
  REAL*8, ALLOCATABLE, DIMENSION(:) :: diam

  !> Drag coefficient for each solid class 
  REAL*8, ALLOCATABLE, DIMENSION(:) :: C_d

  !> Drag coefficient for each solid class 
  REAL*8 :: C_d_input

  !> Solid mass fraction (over the mixture)    
  REAL*8, ALLOCATABLE, DIMENSION(:) :: sepsolid

  !> The fraction of particles of given size   
  REAL*8, ALLOCATABLE, DIMENSION(:) :: solidmassflux_fract

  !> Settling velocity of particles (metres\second)  
  REAL*8, ALLOCATABLE, DIMENSION(:) :: v_s

  !> Term to describe sedimentation for a given particle class at a location
  REAL*8, ALLOCATABLE, DIMENSION(:) :: S 
  
  !> Accumulation Rate
  REAL*8 :: Acc_rate
  
  !> Total Accumulation Rate
  REAL*8 :: Total_Acc_rate
  
  !> Rouse number
  REAL*8 :: Pn

  !> Von Karmans constant (See Valentine 1987)
  REAL*8 :: k = 0.4D0

  !> shear velocity of flow
  REAL*8 :: ustar !> shear velocity of flow
  
  !> Height of boundary thickness layer
  REAL*8 :: h_bl

  !> Height of roughness elements
  REAL*8 :: ks = 1.D0
  
  !> Dimensionless reference level
  REAL*8 :: eta0
  
  !> Dimensionless reference level
  REAL*8 :: wavelengthfreen 
  
  !> Dimensionless reference level
  REAL*8 :: wavelengthshearn
  
  !> Dimensionless reference level
  REAL*8 :: wavelengthfreex 
  
  !> Dimensionless reference level
  REAL*8 :: wavelengthshearx
  
  !> Dimensionless reference level
  REAL*8 :: wavelength_height 
  
  !> Dimensionless height
  REAL*8, ALLOCATABLE, DIMENSION(:) :: eta
  
  !> Brunt-Vaisala fequency
  REAL*8, ALLOCATABLE, DIMENSION(:) :: BV
  
  !> Concentration profile
  REAL*8, ALLOCATABLE, DIMENSION(:) :: SS0
  
  !> Normalised Flux
  REAL*8, ALLOCATABLE, DIMENSION(:) :: norm_flux
  
  !> Initial Particle Flux
  REAL*8, ALLOCATABLE, DIMENSION(:) :: ipf
  
  !> Final Particle Flux
  REAL*8, ALLOCATABLE, DIMENSION(:) :: fpf

  !> Actual number of solid particles specified in input file   
  INTEGER :: iclass
  INTEGER :: final = 9.D0
   
  !> Total sedimentation from mixture 
  REAL*8 :: sumsed

  REAL*8 :: delta_phi
  REAL*8 :: phi_left
  REAL*8 :: phi_right
  REAL*8, ALLOCATABLE, DIMENSION(:) :: phi_dist
  REAL*8 :: left_tail_fract 
  REAL*8 :: right_tail_fract
  REAL*8 :: rho_one

  !
  SAVE
  !----------------------------------------------------------------------
CONTAINS
  !----------------------------------------------------------------------

  !******************************************************************************
  !> \brief Allocating particles variables
  !> This subroutine allocates the variables defining the particle sedimentation. 
  !******************************************************************************
  SUBROUTINE initialize_particles
    INTEGER :: i

    no_bins = iclass

    ! Define a normal distribution

    delta_phi= (max_phi - min_phi) / no_bins

    phi_left = min_phi 

    DO i=1,no_bins

       phi_dist(i) = phi_left + 0.5D0 * delta_phi

       phi_right = min_phi + (i) * delta_phi

       ! Convert grainsize to m   
       diam(i) = 1.D-3 * 2.0**(-phi_dist(i))
    
       IF ( diam(i) .LE. diam1 ) THEN

          rhosol(i) = min_rho

       ELSEIF ( diam(i) .LE. diam2 ) THEN

          rhosol(i) = min_rho + ( diam(i) - diam1 ) / ( diam2 - diam1 ) *       &
               ( max_rho - min_rho)

       ELSE
          
          rhosol(i) = max_rho
   
       END IF
       
       C_d(i) = C_d_input

       fracsolid(i) = cdf(phi_right) - cdf(phi_left)
       
       phi_left = phi_right

    END DO

    left_tail_fract= cdf(min_phi)
    right_tail_fract= 1 - cdf(max_phi)

    If (left_tail_fract .GT. fracsolid(1)) THEN

       fracsolid(1) = fracsolid(1) + left_tail_fract

    END IF

    If (right_tail_fract .GT. fracsolid(no_bins)) THEN

       fracsolid(no_bins) = fracsolid(no_bins) + right_tail_fract

    END IF

    fracsolid(1:no_bins) = fracsolid(1:no_bins)/sum(fracsolid(1:no_bins))
    
    RETURN
    
  END SUBROUTINE initialize_particles

  !******************************************************************************
  !> \brief Defining cumulative density function of particle size
  !
  !> This function produces a normal grainsize distribution using the parameters in input file
  !> 
  !******************************************************************************

  FUNCTION cdf(phi_in)
    !
    IMPLICIT NONE

    REAL*8 :: cdf
    REAL*8, INTENT(IN) :: phi_in

    cdf = 0.5D0 * (1.D0 + erf((phi_in - mean_phi)/(sigma_phi * SQRT(2.D0))))

  END FUNCTION cdf

  !******************************************************************************
  !> \brief Particle inizialization
  !> This subroutine allocates the variables defining the particle physical characteristics 
  !> for each class.
  !******************************************************************************
  SUBROUTINE allocate_particles
    !
    IMPLICIT NONE
    !
    ALLOCATE(solidmassflux_fract(iclass))
    ALLOCATE(v_s(iclass))
    ALLOCATE(S(iclass))
    ALLOCATE(diam(iclass))
    ALLOCATE(rhosol(iclass))
    ALLOCATE(fracsolid(iclass))
    ALLOCATE(phi_dist(iclass))
    ALLOCATE(C_d(iclass))
    ALLOCATE(norm_flux(iclass))
    ALLOCATE(ipf(iclass))
    ALLOCATE(fpf(iclass))
    ALLOCATE(eta(final))
    ALLOCATE(SS0(final))
    ALLOCATE(BV(final))    
    !
    RETURN
  END SUBROUTINE allocate_particles
  !----------------------------------------------------------------------
END MODULE particles_module
!----------------------------------------------------------------------
