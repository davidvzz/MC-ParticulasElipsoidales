      program ellipses
      Implicit None
      real*8, parameter   :: a=0.75, b=0.2
      real*8, parameter :: pi=4.D0*DATAN(1.D0)
      integer, parameter   :: N=2
      integer  :: i
      real*8, dimension(N)   :: cx, cy, theta
      real*8::   theta1

      !Elegimos aleatoriamente la posición y orientación de las elipses
      call random_number(cx)
      call random_number(cy)
      call random_number(theta)
      
      cx=int(cx*5)
      cy=int(cy*5)
      theta=int(theta*180)

      write(*,*) 'x',cx
      write(*,*) 'y',cy   
      write(*,*) 'ángulo',theta

      !Calculamos el ángulo entre cada elipse
      do i=2, N
         theta1=datan( (cy(i)-cy(1)) / (cx(i)-cx(1)) )
         theta1=theta1*180/pi
      end do
      write(*,*) 'ángulo entre elipses',theta1
      
      open(unit=10, file='datos.dat')
      write(10,*) a, b, theta(1), theta(2), theta1

      end 