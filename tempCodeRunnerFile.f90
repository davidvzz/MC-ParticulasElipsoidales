PROGRAM MAIN
   IMPLICIT NONE
   INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
   INTEGER  :: n
   PARAMETER (n=3)
   INTEGER  :: i,j,k
   REAL(KIND = DP)   :: x(n),y(n),z(n)
   call random_seed()
   call random_number(x)

   write(*,*) x
END PROGRAM MAIN

