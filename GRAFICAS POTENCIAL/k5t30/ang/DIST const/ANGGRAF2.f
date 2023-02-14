        PROGRAM ANGGRAF2
            IMPLICIT NONE
            DOUBLE PRECISION:: a, b, ANGLE, ANGLE2, ANGLE3, dist, U,k,ri
            DOUBLE PRECISION :: A1,A2,A3,A4
            DOUBLE PRECISION :: f, ANGINC
            DOUBLE PRECISION, PARAMETER :: pi=4.D0*DATAN(1.D0)
            INTEGER CONF,T
            NAMELIST /input/ k,b,ANGLE, ANGLE2, ANGLE3, A1, A2, A3, A4
            OPEN(UNIT=1,FILE='SS-5-30.dat',STATUS='unknown')     
            OPEN(UNIT=3,FILE='TT-5-30.dat',STATUS='unknown')   
            
            OPEN(UNIT=106,FILE='r6-5-30.dat',STATUS='unknown') 
            OPEN(UNIT=108,FILE='r8-5-30.dat',STATUS='unknown') !para SS

            OPEN(UNIT=110,FILE='r10-5-30.dat',STATUS='unknown') !para TT
            OPEN(UNIT=112,FILE='r12-5-30.dat',STATUS='unknown')

            WRITE(*,*) 'Configuracion?: 1:SS, 3:TT'
            READ(*,*) CONF

            IF (CONF==1) THEN 
            !PARA SS CON R=3B
                READ(UNIT=1,nml=input)
                dist=2*b

            ELSE IF(CONF.NE.1)THEN
            !Para TT
                READ(UNIT=3,nml=input)
                a=k*b  
                dist=2*a
            END IF

            WRITE(*,*) 'ri?: '
            READ(*,*) ri
            T=ri+100
            ANGINC=0
            f=0
            WRITE(*,*) 'angle3'
            WRITE(*,*) ANGLE3
            do while (f.LE.0.5) 

                U=-A1*cos(2*ANGLE2+2*ANGLE3)*( b /
     +                 ( ri - A2*dist + A3 ))**A4
                
                WRITE(T,*) f,U
                ANGINC=ANGINC+pi/180 !rand 
                ANGLE3=ANGLE3+pi/180 !radianes
                f=ANGINC/pi
                
            end do 
            WRITE(*,*) 'angle3 after'
            WRITE(*,*) ANGLE3
            end PROGRAM