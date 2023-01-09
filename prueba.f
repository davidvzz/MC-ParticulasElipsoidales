      program prueba
      implicit none
      integer :: i,j, i_old
      j=0
      do 10 i = 1, 10
    8   if (j==9) then
            write(*,*) i,j
        end if
        i_old=j+1
        do 10 j = 1, i
            !if (i)
            if (j==9) GOTO 8
   10 end do 
      end 