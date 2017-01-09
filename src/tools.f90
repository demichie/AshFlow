!*******************************************************************************
!> \brief Tools
!*******************************************************************************

MODULE tools

  IMPLICIT NONE

CONTAINS

  !*********************************************************************
  !> \brief Program end
  !
  !> This routine ends the program according to the flag iflag:\n
  !> - iflag = -1 => Abnormal termination
  !> - iflag =  1 => Correct  termination
  !> .
  !> \date 23/11/2008
  !> \param   iflag     termination flag      (\b input)
  !> \param   message   termination message   (\b output)
  !********************************************************************

  SUBROUTINE runend(iflag,message)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iflag 
    CHARACTER, INTENT(IN) :: message*(*)

    IF ( iflag .EQ. -1 ) THEN 

       !
       !***  Abnormal termination of the program
       !

       WRITE(*,*) message
       WRITE(*,*) 'CPIUC: *** Abnormal Termination ***'
       WRITE(*,*) ' '
       WRITE(*,*) 'Press return to close the window'
       !       READ(*,*) 

       STOP

    ELSE IF ( iflag .EQ. 1 ) THEN 

       !
       !***  Normal termination of the program
       !

       WRITE(*,*) 'CPIUC: *** The program has finished sucessfully ***'
       WRITE(*,*) ' '
       WRITE(*,*) 'Press return to close the window' 
       !       READ(*,*) 

       STOP 

    END IF

    RETURN

  END SUBROUTINE runend

  !**************************************************************************
  !> \brief 1D Interpolation
  !
  !> This function interpolates values from a discrete function y=f(x). \n
  !> The function is stored at n discrete points in an array f(2,n), where:\n
  !> - f(1,i) = x-coordinate of the point i  (i= 1:n);
  !> - f(2,i) = f(x)         of the point i  (i= 1:n).
  !> .
  !> It returns the value of f(x) at a point x (where f(1,1) < x < f(1,n) ).\n
  !> NOTE: It is assumed that the points are already increasingly ordered in x
  !> \date 13/11/2008
  !> \param  getfx    interpolated value  (\b output)
  !> \param  f        array (x,f(x))      (\b input)
  !> \param  x        interp. point       (\b input)   
  !**************************************************************************

  SUBROUTINE interp1D( getfx,f,x)


    IMPLICIT NONE

    REAL*8, INTENT(OUT) :: getfx

    REAL*8, DIMENSION(:,:), INTENT(IN)  :: f
    REAL*8, INTENT(IN)  :: x

    !*** local variables
    INTEGER :: n, i
    REAL*8 ::  u , x0 , x1 , f0 , f1 , f2 , f_1
    LOGICAL :: search

    n = SIZE( f,2 )

    !
    !***  x is equal or lower than the first point 
    !
    IF ( x .LE. f(1,1) ) THEN

       getfx = f(2,1)
       RETURN

    END IF

    !
    !***  x is equal or greater than the last point 
    !

    IF ( x .GE. f(1,n) ) THEN

       getfx = f(2,n)
       RETURN

    END IF

    !
    !***  x lies between the first and the second point. Interpolate
    !***  using a second order forward
    !

    IF ( ( x .GT. f(1,1) ) .AND. ( x .LE. f(1,2) ) ) THEN

       x0 = f(1,1)
       f0 = f(2,1)
       x1 = f(1,2)
       f1 = f(2,2)
       f2 = f(2,3)
       u = (x-x0)/(x1-x0)
       getfx = f0 + u*(f1-f0)+ 0.5d0*u*(u-1)*(f2-2.d0*f1+f0)

       RETURN

    END IF

    !
    !***  For the rest of points, search and interpolate using a second
    !***  order centred 
    !

    search = .TRUE.

    i = 1

    DO WHILE (search)

       i  = i+1
       x0 = f(1,i  )
       x1 = f(1,i+1)

       IF ( ( x .GT. x0 ) .AND. ( x .LE. x1 ) ) THEN

          f0  = f(2,i  )
          f1  = f(2,i+1)
          f_1 = f(2,i-1)
          search = .FALSE.

       END IF

    END DO

    u     = (x-x0)/(x1-x0)
    getfx = f0 + 0.5d0*u*(f1-f_1)+ 0.5d0*u*u*(f1-2.d0*f0+f_1)

    RETURN

  END SUBROUTINE interp1D

  !*************************************************************************
  !> \brief Function Derivative
  !
  !> This routine computes the discrete derivarive of a function f(x)
  !> at the base points where. The function f(2,n) and its derivative
  !> df(2,n) are stored as:\n
  !> - f(1,i) x-coordinate of the the point i
  !> - f(2,i) f(x) at the point i
  !> - df(1,i) x-coordinate of the the point i
  !> - df(2,i) df(x)/dx at the point i
  !> .    
  !> Derivatives are calculated using a second order approximation\n
  !> NOTE: Points do not need to be equally spaced, but they must be ordered
  !>       on increasing values of x.
  !> \date 13/11/2008
  !> \param   f    function      (\b input)
  !> \param   df   derivative    (\b output)
  !> \param   n    n. of points  (\b input)
  !***************************************************************************  

  SUBROUTINE dfdx(f,df,n)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL*8, INTENT(IN)  :: f(2,n)
    REAL*8, INTENT(OUT) :: df(2,n)


    !*** local variables

    INTEGER :: i

    !
    !***  Set the seme base points
    !

    DO i=1,n

       df(1,i) = f(1,i)

    END DO

    !
    !***  First point (i=1). Second order forward approximation
    !
    df(2,1) = (-f(2,3)+4.d0*f(2,2)-3.0d0*f(2,1))/(f(1,3)-f(1,1))

    !
    !***  Last point  (i=n). Second order backward approximation.
    !
    df(2,n) = (3.d0*f(2,n)-4.d0*f(2,n-1)+f(2,n-2))/(f(1,n)-f(1,n-2))

    !
    !***  Central points. Second order centred.
    !

    DO i = 2,n-1

       df(2,i) = (f(2,i+1)-f(2,i-1))/(f(1,i+1)-f(1,i-1)) 

    END DO

    RETURN

  END SUBROUTINE dfdx

END MODULE tools
