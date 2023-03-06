        PROGRAM TESTGRAF
            Implicit None
            DOUBLE PRECISION:: a, b, ANGLE, ANGLE2, ANGLE3, dist, U,k
            DOUBLE PRECISION :: A1,A2,A3,A4,ri(3)
            DOUBLE PRECISION, PARAMETER :: pi=4.D0*DATAN(1.D0)
            INTEGER T,i
            NAMELIST /input/ a,k,b,ANGLE, ANGLE2, ANGLE3, A1, A2, A3, A4
            
            OPEN(UNIT=1,FILE='SS-2-80.dat',STATUS='unknown')     
            OPEN(UNIT=2,FILE='ST-2-80.dat',STATUS='unknown')     
            OPEN(UNIT=3,FILE='TT-2-80.dat',STATUS='unknown')     

            OPEN(UNIT=101, FILE='outSSk2.dat', STATUS='unknown')
            OPEN(UNIT=102, FILE='outSTk2.dat', STATUS='unknown')
            OPEN(UNIT=103, FILE='outTTk2.dat', STATUS='unknown')
            

            ri = (/ 2,3,4 /)
            do i=1, 3
                T=i+100
                IF (i==1) THEN 
                !Para angulos SS
                    READ(unit=1,nml=input) 
                    dist=2*b
                    write(T,*) 'k= ',k,'a=',a,'b=',b,'r/b=',ri(i)
                    write(T,*) 'x   ', 'SS'
                ELSE IF (i==2) THEN
                !PARA CONF ST
                    READ(unit=2,nml=input)
                    dist=a+b
                    write(T,*) 'k= ',k,'a=',a,'b=',b,'r/b=',ri(i)
                    write(T,*) 'x   ', 'ST'
                ELSE 
                    !PARA CONF TT
                    READ(UNIT=3,nml=input)
                    dist=2*a
                    write(T,*) 'k= ',k,'a=',a,'b=',b,'r/b=',ri(i)
                    write(T,*) 'x   ', 'TT'
                END IF


                do while (ri(i)<=10.1)
                    
                    U=-A1*cos(2*ANGLE2+2*ANGLE3)*( b /
     +                             ( ri(i)*b - A2*dist + A3 ))**A4
                    WRITE(T,*) ri(i),U 
                    ri(i)=ri(i)+0.1
                end do

            end do
        END PROGRAM

