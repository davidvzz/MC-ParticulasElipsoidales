            PROGRAM TESTGRAF2
            Implicit None
            DOUBLE PRECISION:: a, b, ANGLE, ANGLE2, ANGLE3, dist, U,k,ri
            DOUBLE PRECISION :: A1,A2,A3,A4
            DOUBLE PRECISION, PARAMETER :: pi=4.D0*DATAN(1.D0)
            INTEGER CONF,T
            NAMELIST /input/ k,b,ANGLE, ANGLE2, ANGLE3, A1, A2, A3, A4
            
            OPEN(UNIT=1,FILE='SS-5-30.dat',STATUS='unknown')     
            OPEN(UNIT=2,FILE='ST-5-30.dat',STATUS='unknown')     
            OPEN(UNIT=3,FILE='TT-5-30.dat',STATUS='unknown')     

            WRITE(*,*) 'Configuracion?: 1:SS, 2:ST, 3:TT'
            READ(*,*) CONF
            OPEN(UNIT=101, FILE='outSS5-30.dat', STATUS='unknown')
            OPEN(UNIT=102, FILE='outST5-30.dat', STATUS='unknown')
            OPEN(UNIT=103, FILE='outTT5-30.dat', STATUS='unknown')
            T=CONF+100
            IF (CONF==1) THEN 
            !Para angulos SS
                READ(unit=1,nml=input) 
                dist=2*b
                a=k*b  
            ELSE IF (CONF==2) THEN
            !PARA CONF ST
                READ(unit=2,nml=input)
                a=k*b
                dist=a+b
            ELSE 
                !PARA CONF TT
                READ(UNIT=3,nml=input)
                a=k*b  
                dist=2*a
            END IF

            ri=2.0

            do while (ri.LE.18)
                
               U=-A1*cos(2*ANGLE2+2*ANGLE3)*( b /
     +                             ( ri - A2*dist + A3 ))**A4
                WRITE(T,*) ri,U 
                ri=ri+0.5
                
            end do
        END PROGRAM