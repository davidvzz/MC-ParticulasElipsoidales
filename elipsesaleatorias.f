      program ellipses
      Implicit None
      real, parameter   :: a=1, b=0.5
      
      integer, parameter   :: N=2
      integer  :: i, seed
      real  ::  theta1, theta2, theta3
      real, dimension(N)   :: cx, cy, theta

      !Elegimos aleatoriamente la posición y orientación de las elipses
      call random_number(cx)
      call random_number(cy)
      call random_number(theta)
      
      cx=int(cx*10)
      cy=int(cy*10)
      theta=int(theta*180)

      write(*,*) cx
      write(*,*) cy   
      write(*,*) theta

      !Calculamos el ángulo entre cada elipse
      do i=1, N
         theta1=theta(1)-theta(i)
      end do
      write(*,*) theta1
      end 