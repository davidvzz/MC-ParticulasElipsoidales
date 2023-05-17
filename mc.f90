program main
   implicit none
   integer, parameter :: DP = SELECTED_real_KIND(14)
   integer,parameter :: n=100
   real(KIND = DP)   :: rx(n),ry(n),ra(n)
   real(kind = DP), parameter :: rho=0.01d0, aspect_ratio=3.d0
   real(kind = DP)   :: s

   
   !call conf_inicial(rx, ry, ra, n, rho, aspect_ratio)
   !call monte_carlo(rx,ry,ra,n,s)
end program main

subroutine conf_inicial(rx, ry, ra, n, rho, aspect_ratio)
   implicit none
   integer, parameter :: DP = SELECTED_real_KIND(14)
   integer, intent(in)  :: n
   integer  :: i,j
   real(KIND = DP), intent(inout)   :: rx(n),ry(n),ra(n)
   real(KIND = DP)   :: x,y,rr
   real(KIND = DP)   :: xl,rho,s,aspect_ratio
   logical  :: condition


   !Parámetros caja
   xl=(n/rho)**(0.5)
   s=1/xl

   !Creamos una configuración inicial
   call random_number(rx)
   call random_number(ry)
   call random_number(ra)
   ra=ra*180

   write(*,*) rx(1),ry(1),ra(1)
   i=1
   j=1
   condition=.false.

   !Ciclo para revisar traslapes
   do while (i<=n)
      j=1
      do while ( j<=n .and. .not. condition )
         if (i == j) then
            j=j+1
            cycle
         end if
         write(*,*) "i",i,"j",j

         x=rx(j)-rx(i)
         y=ry(j)-ry(i)

         call imagen_minima(x,y,rr)
         call traslape(x, y, ra(j), ra(i), aspect_ratio, rr, s, condition)

         write(*,*) condition
         write(*,*) "-------------------"
         j=j+1
      end do
      
      !En caso de existir traslape modificamos la posición de la partícula 
      !y volvemos a repetir el paso
      if (condition) then
         call random_number(rx(i))
         call random_number(ry(i))
         call random_number(ra(i))
         ra(i)=ra(i)*180
         condition = .false.
         cycle
      end if

      i=i+1
   end do
   open(unit=2,file="conf_inicial.dat",status="unknown")


   do i = 1,n
      write(2,*) rx(i), ry(i), s*aspect_ratio, s, ra(i)
   end do

END subroutine conf_inicial

subroutine monte_carlo(rx, ry, ra, n, s)
   implicit none
   integer, parameter :: DP = SELECTED_real_KIND(14)
   integer, intent(in)  :: n
   real(KIND = DP), intent(in)   :: s
   real(KIND = DP), intent(inout)   :: rx(n),ry(n),ra(n)
   real(kind = DP)  :: ireal, xnew, ynew, anew
   integer  :: ncount=1, nmc, ndata=50, i, j, k

   
   do while(ncount<nmc)
      !escogemos una partícula al azar 
      call random_number(ireal)
      i=int(ireal*n)+1
      ncount=ncount+1

      !movemos la partícula
      call random_number(xnew)
      call random_number(ynew)
      call random_number(anew)
      xnew=rx(i)+max_displ*(xnew-0.5)
      ynew=ry(i)+max_displ*(ynew-0.5)
      anew=ra(i)+max_displ_ang*(anew-0.5)

      call condicion_periodica(xnew, ynew, anew)
   end do 

end subroutine monte_carlo

subroutine imagen_minima(x,  y, rr)
   implicit none
   integer, parameter :: DP = SELECTED_real_KIND(14)
   real(KIND = DP), intent(inout) :: x, y
   real(KIND = DP), intent(out)  :: rr

   if (x > 0.5) then
      x = x-1
   else if (x < -0.5) then
      x = x+1
   end if

   if (y > 0.5) then
      y = y-1
   else if (y < -0.5) then
      y = y+1
   end if

   rr=x*x+y*y

end subroutine imagen_minima

subroutine traslape(x,y,ra1,ra2,ar,rr,s,condition)
   implicit none
   integer, parameter :: DP = SELECTED_real_KIND(14)
   real(KIND = DP), intent(in) :: ar,x,y,ra1,ra2,s,rr
   logical, intent(out)  ::   condition
   real(KIND = DP)   :: ggg,f1,f2,phi
   real(KIND = DP), parameter   ::  pi=4*datan(1.d0)


   ggg = 2.00 +(AR-(1/AR))**2*(sin((RA1-RA2)*PI/180))**2
   F1 = -((X*cos(RA2*PI/180)+ Y*sin(RA2*PI/180))**2)
   F1 = F1/(AR*S/2)**2
   F1 = F1+1.00 + GGG

   F2 = -((X*cos(RA1*PI/180)+ Y*sin(RA1*PI/180))**2)
   F2 = F2/(AR*S/2)**2
   F2 = F2+1.00 +GGG
   F2 = F2-((Y*cos(RA1*PI/180)- X*sin(RA1*PI/180))**2)/(S*S/4)
   PHI = 4*(F1**2-3*F2)*(F2**2-3*F1)-(9-F1*F2)**2

   condition =(PHI.lt.0.or.(PHI.gt.0.and.F1>=0.and.F2>=0).or.RR.lt.S*S)
end subroutine traslape

subroutine condicion_periodica(x, y, a)
   implicit none
   integer, parameter :: DP = SELECTED_real_KIND(14)
   real(KIND = DP), intent(inout)   :: x,y,a

   if (x>1.d0) then
      x=x-1
   else if (x<0) then 
   end if
      
   
end subroutine