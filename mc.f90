PROGRAM MAIN
   IMPLICIT NONE
   INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
   INTEGER  :: n
   PARAMETER (n=700)
   INTEGER  :: i,j
   REAL(KIND = DP)   :: rx(n),ry(n),ra(n)
   REAL(KIND = DP)   :: x,y,rr,aspect_ratio=3.d0
   REAL(KIND = DP)   :: xl,rho,s
   logical  :: condition

   rho=0.1d0

   !Parámetros caja
   xl=(n/rho)**(0.5)
   s=1/xl

   call random_number(rx)
   call random_number(ry)
   call random_number(ra)  
   ra=ra*180

   write(*,*) rx(1),ry(1),ra(1)
   i=1
   j=1
   condition=.false.
   do while (i<=n)
      j=1
      do while ( j<=n .and. .not.condition )
         if (i == j) then
            j=j+1
            cycle
         end if
         write(*,*) "i",i,"j",j
         
         x=rx(j)-rx(i)
         y=ry(j)-ry(i)
         call IMAGEN_MINIMA(x,y,rr)
         call traslape(x, y, ra(j), ra(i), aspect_ratio, rr, s, condition)
         write(*,*) condition
         write(*,*) "-------------------"
         j=j+1
      end do
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

   write(*,*) ra

   do i = 1,n
      write(2,*) rx(i), ry(i), s*aspect_ratio, s, ra(i)
   end do 

END PROGRAM MAIN

subroutine IMAGEN_MINIMA(x,  y, rr)
   IMPLICIT NONE
   INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
   REAL(KIND = DP), intent(inout) :: x, y
   REAL(KIND = DP), intent(out)  :: rr

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
   
end subroutine IMAGEN_MINIMA

subroutine traslape(x,y,ra1,ra2,ar,rr,s,condition)
   IMPLICIT NONE
   INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
   REAL(KIND = DP), intent(in) :: ar,x,y,ra1,ra2,s,rr
   logical, intent(out)  ::   condition
   REAL(KIND = DP)   :: ggg,f1,f2,phi
   REAL(KIND = DP), PARAMETER   ::  pi=4*datan(1.d0)

   
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

