        PROGRAM TESTGRAF
            Implicit None
            DOUBLE PRECISION:: a, b, ANGLE, ANGLE2, ANGLE3, dist, U,k,ri
            DOUBLE PRECISION :: A1,A2,A3,A4
            DOUBLE PRECISION, PARAMETER :: pi=4.D0*DATAN(1.D0)
            INTEGER CONF,T
            NAMELIST /input/ k,b,ANGLE, ANGLE2, ANGLE3, A1, A2, A3, A4
            
            OPEN(UNIT=1,FILE='SS-2-80.dat',STATUS='unknown')     
            OPEN(UNIT=2,FILE='ST-2-80.dat',STATUS='unknown')     
            OPEN(UNIT=3,FILE='TT-2-80.dat',STATUS='unknown')     

            WRITE(*,*) 'Configuracion?: 1:SS, 2:ST, 3:TT'
            READ(*,*) CONF
            OPEN(UNIT=101, FILE='outSS.dat', STATUS='unknown')
            OPEN(UNIT=102, FILE='outST.dat', STATUS='unknown')
            OPEN(UNIT=103, FILE='outTT.dat', STATUS='unknown')
            T=CONF+100

            ri=4.0

            IF (CONF==1) THEN 
            !Para angulos SS
                READ(unit=1,nml=input) 
                dist=2*b
                a=k*b  
                write(T,*) 'k= ',k,'a=',a,'b=',b,'r=',ri
                write(T,*) 'x   ', 'SS'
            ELSE IF (CONF==2) THEN
            !PARA CONF ST
                READ(unit=2,nml=input)
                a=k*b
                dist=a+b
                write(T,*) 'k= ',k,'a=',a,'b=',b,'r=',ri
                write(T,*) 'x   ', 'ST'
            ELSE 
                !PARA CONF TT
                READ(UNIT=3,nml=input)
                a=k*b  
                dist=2*a
                write(T,*) 'k= ',k,'a=',a,'b=',b,'r=',ri
                write(T,*) 'x   ', 'TT'
            END IF



            do while (ri<10)
                
               U=-A1*cos(2*ANGLE2+2*ANGLE3)*( b /
     +                             ( ri - A2*dist + A3 ))**A4
                WRITE(T,*) ri,U 
                ri=ri+0.2
                
            end do
        END PROGRAM


        