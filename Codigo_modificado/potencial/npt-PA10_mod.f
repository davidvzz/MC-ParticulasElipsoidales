      PROGRAM NPT6 
      !                                 FEBRERO 1999
                                      !Noviembre 2017
                                      !Enero 2023
      !
      !     SIMULACION MONTE CARLO Elipses
         !Graficar Energía promedio vs numero de pasos de monte carlo
         !Ver que se genere archivos de posición finales e iniciales

      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NPART=2000,NACC=20,NG=30000)
      COMMON /BPOSITN/ RX(NPART),RY(NPART),RA(NPART),ACC(NACC),G(NG),AR,
     +                 CLU(NPART,NPART)
      COMMON /BCONSTR/ PI,ETA,RHO,XLAMBDA,XLAM2,XLAM3,SL2,SQL2,
     +                 SL3,SQL3,XN,TEMP,DISPL,DISPLAng, DISPLClu,XHISTG,
     +                 S,SS,SSLL,SL,XC,YC,YCINV,
     +                 XL,YL,Y2,RA3,A1, A2, A3, A4,a,b,tension,b_sigma
      COMMON /BCONSTG/ DR1,DR2,DR3,DR4,DR5,PHI,TAU
      COMMON /BCOSTRX/ SL4,SQL4,SL5,SQL5,XLAM4,XLAM5
      COMMON /BCOSTXX/ DR6,SL6,SQL6,XLAM6,E6
      COMMON /BCONSTI/ N,NGOFR,LGOFR,NMOVE,NMOVE2,NSUB,NGOFR0,ISEED
      COMMON /BCONSTV/ PRESS,VOL,SDISPL,NSET,LRHO
      LOGICAL LGOFR
      PI = 4.D0*DATAN(1.D0)
      OPEN(UNIT=15,FILE='graficas/conf_inicial.dat',STATUS='unknown')
      open(unit=16,file="graficas/conf_final.dat",STATUS="unknown")
      OPEN(UNIT=6,FILE='datos/npt6.dat',STATUS='unknown')

      OPEN(UNIT=8,FILE='datos/npt6.tau',STATUS='unknown')
      !  UNIDAD 10 NPT5.SIG GUARDA MUESTREO EN SIGMA
      OPEN(UNIT=10,FILE='datos/npt6.sigma',STATUS='unknown')
      open(unit=20,file='graficas/EvsN.dat',status='unknown')
      
      CALL START
      CALL MCARLO
      CALL FINISH
      !CALL RADIAL

      CLOSE(UNIT=6)
      STOP
      END
      !****************************************************************
      SUBROUTINE START 

       !Lee variables de sistema y estado (rho,temp) de pozos.in
       !al inicio y pozos.old si es continuacion
       !Genera la configuracion inicial de fcc.

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NPART=2000,NACC=20,NG=30000)
      COMMON /BPOSITN/ RX(NPART),RY(NPART),RA(NPART),ACC(NACC),G(NG),AR,
     +                 CLU(NPART,NPART)
      COMMON /BCONSTR/ PI,ETA,RHO,XLAMBDA,XLAM2,XLAM3,SL2,SQL2,
     +                 SL3,SQL3,XN,TEMP,DISPL,DISPLAng, DISPLClu,XHISTG,
     +                 S,SS,SSLL,SL,XC,YC,YCINV,
     +                 XL,YL,Y2,RA3,A1, A2, A3, A4,a,b,tension,b_sigma
      COMMON /BCONSTG/ DR1,DR2,DR3,DR4,DR5,PHI,TAU
      COMMON /BCOSTRX/ SL4,SQL4,SL5,SQL5,XLAM4,XLAM5
      COMMON /BCOSTXX/ DR6,SL6,SQL6,XLAM6,E6
      COMMON /BCONSTI/ N,NGOFR,LGOFR,NMOVE,NMOVE2,NSUB,NGOFR0,ISEED
      COMMON /BCONSTV/ PRESS,VOL,SDISPL,NSET,LRHO

      LOGICAL LGOFR
      !LOGICAL LRHO

      !Lectura de datos de entrada por medio de una namelist
      NAMELIST /input/ N,RHO,TEMP,PHI,XLAMBDA,XLAM2,XLAM3,XLAM4,XLAM5,
     +                 XLAM6,NRUN,NMOVE,NMOVE2,NSUB,DISPL,DISPLAng ,
     +                 DISPLClu,RA3,IAV,NAC,NCX,NCY,QX,QY,SDISPL,NGOFR,
     +                 LGOFR,XHISTG,AR, A1, A2, A3, A4,a,b,tension


      ! ACUMULADORES A CERO
      DO I=1,NACC
         ACC(I)=0.0D00
      end do

      DO I=1,NG
        G(I)=0.0D00
      end do

      CLU(:,:)=0.0D0
      ! LEE PARAMETROS DE ENTRADA: N,TEMP E.T.C.
      OPEN (UNIT=1,FILE='npt6.in',STATUS='old')
      READ(unit=1, nml=input)

      TEMP=TEMP+273.15d0
      ! ESCRIBE PARAMETROS DE ENTRADA
      LRHO=0
      WRITE(6,100)
      WRITE(6,101) N,RHO,TEMP,NRUN,NMOVE,NSUB,DISPL,DISPLAng,IAV,NAC,
     + NCX,NCY,NGOFR
      WRITE(6,*) 'LAMBDA1=',XLAMBDA
      WRITE(6,*) 'LAMBDA2=',XLAM2
      WRITE(6,*) 'LAMBDA3=',XLAM3
      WRITE(6,*) 'LAMBDA4=',XLAM4
      WRITE(6,*) 'LAMBDA5=',XLAM5
      WRITE(6,*) 'LAMBDA6=',XLAM6
      WRITE(6,*) 'PHI=',PHI

      IF (NRUN.EQ.0.0D00) WRITE(6,102)
      IF (NRUN.NE.0.0D00) WRITE(6,103)
      IF (IAV.EQ.0.0D00)  WRITE(6,104)
      IF (IAV.NE.0.0D00)  WRITE(6,105)
      IF (LGOFR)     WRITE(6,106) XHISTG

      ! CONVERSION A UNIDADES DE PROGRAMA (Lx = 1)
      PRESS=1.1547005384D0*PHI
      XN=DFLOAT(N)

      ! SE CALCULA EL TAMANO DE LA MUESTRA PARA BARRAS DE ERROR EN SIGMA
      XNMOVE=DFLOAT(NMOVE)
      XSET=XNMOVE/100D0
      NSET=NINT(XSET)

      ! CADA NSET MOVIMIENTOS SE MUESTREA SIGMA
      XNCX=DFLOAT(NCX)
      XNCY=DFLOAT(NCY)
      XC=(XNCX*QX)/(XNCX*QX)
      YC=(XNCY*QY)/(XNCX*QX)
      Y2=YC/2.0D00
      YCINV=2.0D00/YC

      IF (NRUN.EQ.0) THEN
         ! CALCULA SIGMA EN UNIDADES DE CAJA
         ETA=RHO*PI/4.0D00
         XL=(XN/(RHO*YC))**(1.0D00/2.0D00)
         YL=XL*YC
         S=1.0D00/XL
         !  CONSTRUYE LATIZ INICIAL CON SIGMA DADA
         ISEED=-123456789

         DO I=1,N
            !Genera una posicion y orientacion para npart (se podría usar allocate para no tener 2 variables distintas)(npart>N)
    9       RX(I)=RAN2(ISEED)
            RY(I)=YC*RAN2(ISEED)
            RA(I)=180*RAN2(ISEED)

            DO J=1,I
               IF (J.NE.I)then !No hacemos nada si i=j (es la misma partícula)
                  !Calculamos distancias entre particulas
                  X=RX(J)-RX(I)
                  Y=RY(J)-RY(I)

                  ! convencion de imagen m!nima
                  IF (X.GT.0.5D00) THEN
                     X=X-1.0D00
                  ELSE IF (X.LT.-0.5D00) THEN
                     X=X+1.0D00
                  END IF

                  IF (Y.GT.Y2) THEN
                     Y=Y-YC
                  ELSE IF (Y.LT.-Y2) THEN
                     Y=Y+YC
                  END IF

                  RR=X*X+Y*Y
                  !!!!CONDICIONES DE TRASLAPE
                  GGG = 2.00 +(AR-(1/AR))**2*
     +             (sin((RA(J)-RA(I))*PI/180))**2
                  F1=-((X*cos(RA(I)*PI/180)+ Y*sin(RA(I)*PI/180))**2)
                  F1=F1/(AR*S/2)**2
                  F1=F1+1.00 + GGG

                  F1=F1-((Y*cos(RA(I)*PI/180)- X*sin(RA(I)*PI/180))**2)
     +             /(S*S/4)

                  F2=-((X*cos(RA(J)*PI/180)+ Y*sin(RA(J)*PI/180))**2)
                  F2=F2/(AR*S/2)**2
                  F2=F2+1.00 +GGG
                  F2=F2-((Y*cos(RA(J)*PI/180)- X*sin(RA(J)*PI/180))**2)
     +             /(S*S/4)
                  FI=4*(F1**2-3*F2)*(F2**2-3*F1)-(9-F1*F2)**2

                  if (FI.lt.0.or.(FI.gt.0.and.F1>=0.and.F2>=0)
     +             .or.RR.lt.S*S) GOTO 9     !True si se translapan, en ese caso se vuelve a elegir una nueva posición
               end if
            end do
         end do

      !! CHECAR si se ocupa else (nrun es cero y no se modifica)
      ELSE
         OPEN(UNIT=3,FILE='npt6.old',STATUS='OLD')
         IF (LRHO.EQ.1) THEN
            ! CALCULA SIGMA CON DENSIDAD DADA Y CONFIGURACION VIEJA
            ETA=RHO*PI/4.0D0
            XL=(XN/(RHO*YC))**(1.0D00/2.0D00)
            YL=XL*YC
            S=1.0D00/XL
            READ(3,*) RX(1)
         ELSE
          ! LEE VALOR DE SIGMA
            READ(3,*) S
            VOL=YC/(S*S)
            RHO=XN/VOL
            ETA=RHO*PI/4.0D00
            XL=1/S
            YL=XL*YC
         END IF

         !LEE CONFIGURACION VIEJA
         DO 11 I=1,N
            !WRITE(3,*) RX(I),RY(I), RA(I)
            READ(3,*) RX(I),RY(I),Ancho,S, RA(I)
   11    CONTINUE

      ! lee acumuladores viejos
         IF (IAV.EQ.0) GOTO 2

         DO 22 I=1,NACC
            READ(3,*)ACC(I)
   22    CONTINUE

         DO 33 I=1,NG
            READ(3,*)G(I)
   33    CONTINUE

         ! escribe acumulador y configuracion de inicio o re-inicio
         WRITE(6,111) ACC(1)
         WRITE(6,112)
   2  END IF  !!GOTO 2


      !CONTINUA CONVIRTIENDO A UNIDADES DE PROGRAMA

      SS=S*S
      SL=S*XLAMBDA
      SSLL=SL*SL
      SL2=S*XLAM2
      SQL2=SL2*SL2
      SL3=S*XLAM3
      SQL3=SL3*SL3
      SL4=S*XLAM4
      SQL4=SL4*SL4
      SL5=S*XLAM5
      SQL5=SL5*SL5
      SL6=S*XLAM6
      SQL6=SL6*SL6
      DISPL=DISPL*S
      VOL=YC/(S*S)

      !ESCRIBE LAS UNIDADES CONVERTIDAS
      WRITE(6,107) RHO
      WRITE(6,108) XC,YC
      WRITE(6,109) XL,YL
      WRITE(6,110) S

      DO 44 I=1,N
         WRITE(15,*) RX(I),RY(I),S*AR,S,RA(I)
   44 CONTINUE

      ! formatos de escritura de START
  100 FORMAT(1X,'*********************************************',/
     +         ,' ***SIMULACION MONTE CARLO DE POZOS CUADRADOS ***',/
     +         ,' *********************************************',/)
  101 FORMAT(1X,'N= ',I4,' RHO= ',F8.6,' TEMP= ',F8.4,/
     +       ,' NRUN=',I2,' NMOVE=',I10,' NSUB=',I9,' DISPL=',F10.8,/
     +       ,' DISPLAng=',F10.8,' IAV=',I2,' NAC,NCX,NCY=',4I4,/
     +       ,' NGOFR=',I9,/,/)
  102 FORMAT(1X,'INICIO DE UNA LATIZ',/)
  103 FORMAT(1X,'INICIO DE CONFIGURACION PREVIA',/)
  104 FORMAT(1X,'ACUMULADOR A CERO',/)
  105 FORMAT(1X,'PROMEDIOS CONTINUAN',/)
  106 FORMAT(1X,'CALCULA G(R)','    XHISTG=',F10.4,/)
  107 FORMAT(1X,'DENSIDAD REDUCIDA=',F10.8,/)
  108 FORMAT(1X,'TAMA�O DE CAJA EN UNIDADES DE PROGRAMA :',
     +         8X,'XC=',F10.8,' YC=',F10.8,/)
  109 FORMAT(1X,'EN UNIDADES DE SIGMA:  ',8X,'XL=',F10.6,' YL=',F10.6,/)
  110 FORMAT(1X,'SIGMA EN UNIDADES DE LA CAJA:',4X,F10.8,/)
  111 FORMAT(1X,'NUMERO PREVIO DE CONFIGURACIONES',F12.0,/)
  112 FORMAT(1X,'CONFIGURACON INICIAL')
  115 FORMAT(1X,2(1X,F10.8))
      RETURN
      END


! ******************************************************************

      SUBROUTINE MCARLO
      ! CALCULA LA ENERGIA POTENCIAL DE LA CONFIGURACION INICIAL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NPART=2000,NACC=20,NG=30000)
      COMMON /BPOSITN/ RX(NPART),RY(NPART),RA(NPART),ACC(NACC),G(NG),AR,
     +                 CLU(NPART,NPART)
      COMMON /BCONSTR/ PI,ETA,RHO,XLAMBDA,XLAM2,XLAM3,SL2,SQL2,
     +                 SL3,SQL3,XN,TEMP,DISPL,DISPLAng, DISPLClu,XHISTG,
     +                 S,SS,SSLL,SL,XC,YC,YCINV,
     +                 XL,YL,Y2,RA3,A1, A2, A3, A4,a,b,tension,b_sigma
      COMMON /BCONSTG/ DR1,DR2,DR3,DR4,DR5,PHI,TAU
      COMMON /BCOSTRX/ SL4,SQL4,SL5,SQL5,XLAM4,XLAM5
      COMMON /BCOSTXX/ DR6,SL6,SQL6,XLAM6,E6
      COMMON /BCONSTI/ N,NGOFR,LGOFR,NMOVE,NMOVE2,NSUB,NGOFR0,ISEED
      COMMON /BCONSTV/ PRESS,VOL,SDISPL,NSET,LRHO
      LOGICAL LGOFR
      LOGICAL ALL
      DOUBLE PRECISION XNEW1(N), YNEW1(N), ANEW1(N)
      DOUBLE PRECISION, PARAMETER:: kb=1.380649D-23
      
      ! inicializa secuencia al azar
      ISEED=-123456794  ! 
         
      a=AR*S/2 !semieje mayor en unidad de caja
      b=S/2    !semieje menor en unidad de caja
      b_sigma=(1d-6/2d0) !semieje mejor en metros (suponiendo que sigma=1micrometro)

      ALL=.true.
      UTOT=0
      X=0
      Y=0
      ANG=0
      ANG2=0

      CALL ENERG(UTOT,ALL,X,Y,ANG,ANG2) ! subrutina energia para calcular la energía de la configuración creada inicialmente

      UPERP=UTOT/XN ! promedio energia 
      ! NTEST=0
      USUBAV=0.0D00
      NCOUNT=0

      NERCONT=0 ! contadores y posiciones en cero 
      NACCPT=0
      RTEST=0.0D00
      NSUB0=NSUB
      NGOFR0=NGOFR
      WRITE(6,101) UPERP
      !comienza la secuencia de Monte Carlo

      OPEN(UNIT=2,FILE='datos/unpt6.dat',STATUS='unknown')
      !OPEN(UNIT=93,FILE='angulos.dat',STATUS='NEW')

      ! escoge particula al azar
    1 I=INT(RAN2(ISEED)*N)+1   !!!GOTO 1
      NCOUNT=NCOUNT+1

      ! NERCONT ES EL CONTADOR QUE AL LLEGAR A NSET INDICA MUESTREAR SIGMA
      NERCONT=NERCONT+1
      RTEST=RTEST+RAN2(ISEED)

      ! desplaza la part!cula
      XNEW=RX(I)+DISPL*(RAN2(ISEED)-0.5D00)
      YNEW=RY(I)+DISPL*(RAN2(ISEED)-0.5D00)

      !angular
      ANEW=RA(I)+DISPLAng*(RAN2(ISEED)-0.5D0)

      !condicion periodica de frontera
      IF(XNEW.GT.1.0D00) THEN
         XNEW=XNEW-1.0D00
      ELSE IF (XNEW.LT.0.0D00) THEN
         XNEW=XNEW+1.0D00
      END IF

      IF (YNEW.GT.YC) THEN
         YNEW=YNEW-YC
      ELSE IF (YNEW.LT.0.0D00) THEN
         YNEW=YNEW+YC
      END IF
      
      IF (ANEW.gt.180) THEN
         ANEW=ANEW-180
      ELSE IF(ANEW.lt.0) then
         ANEW=ANEW+180
      end if

      ! elimina particulas traslapadas
      ! y calcula la nueva energia (provisional)
      UNEW=0.0D00

      DO 2 J=1,N
         IF (J.EQ.I) GOTO 2   !Si es la misma partícula continua con el siguiente numero
         X=RX(J)-XNEW
         Y=RY(J)-YNEW
         ! convencion de imagen m!nima
         IF (X.GT.0.5D00) THEN
            X=X-1.0D00
         ELSE IF (X.LT.-0.5D00) THEN
            X=X+1.0D00
         END IF

         IF (Y.GT.Y2) THEN
            Y=Y-YC
         ELSE IF (Y.LT.-Y2) THEN
            Y=Y+YC
         END IF
         RR=X*X+Y*Y
      

         !!!!CONDICIONES DE TRASLAPE
         GGG = 2.00 +(AR-(1/AR))**2*(sin((RA(J)-ANEW)*PI/180))**2
         F1=1.00 + GGG
         F1=F1-((X*cos(ANEW*PI/180)+ Y*sin(ANEW*PI/180))**2)/(AR*S/2)**2
         F1=F1-((Y*cos(ANEW*PI/180)- X*sin(ANEW*PI/180))**2)/(SS/4)
         F2=1.00 +GGG
         F2=F2-((X*cos(RA(J)*PI/180)+ Y*sin(RA(J)*PI/180))**2)
     +    /(AR*S/2)**2
         F2=F2-((Y*cos(RA(J)*PI/180)- X*sin(RA(J)*PI/180))**2)/(SS/4)
         FI=4*(F1**2-3*F2)*(F2**2-3*F1)-(9-F1*F2)**2

         !FI=0 es tangencia
         !no puntos en com�n es FI positiva y al menos uno de f1, f2 negativos
      
         if (FI.lt.0.or.(FI.gt.0.and.F1>=0.and.F2>=0).or.RR.lt.SS) !TRUE si se translapan
     +    GOTO 3


         if (RR.lt.25*a*a )then !CHECAR condición if
            ALL=.false.
            call ENERG(UNEW, ALL, X, Y, ANEW, RA(J))

         end if
      
    2 CONTINUE   !!!GOTO 2


      UOLD=0.0D00
      DO 4 J=1,N
         IF (J.EQ.I) GOTO 4
         X=RX(J)-RX(I)
         Y=RY(J)-RY(I)
      
         ! convencion de imagen m!nima
         IF (X.GT.0.5D00) THEN
            X=X-1.0D00
         ELSE IF (X.LT.-0.5D00) THEN
            X=X+1.0D00
         END IF
         IF (Y.GT.Y2) THEN
            Y=Y-YC
         ELSE IF (Y.LT.-Y2) THEN
            Y=Y+YC
         END IF
         RR=X*X+Y*Y

         !!!!!!! CAMBIAR
         !!Energia pozos angulares

         if (RR.lt. 25*a*a  )then !CHECAR condición if
            ALL=.false.
            call ENERG(UOLD, ALL, X, Y, RA(I), RA(J))

         end if
    4 CONTINUE
      

      DENERG=tension*(b_sigma*b_sigma)*(UNEW-UOLD)
      

      !Si la energía de la nueva configuración es menor que la energía de laconfiguracioń
      !inicial, entonces aceptamos la nueva configuración
      IF (DENERG.LE.0.0D00) GOTO 5
      
      
      ! compara el factor de BOLTZMANN con n#mero al azar
      RND=RAN2(ISEED)
      IF (RND.GT.EXP(-DENERG/(kb*TEMP))) GOTO 3  !TRUE -> se rechaza la configuración
      
      ! actualiza posicion de la part!cula I
    5 RX(I)=XNEW        !GOTO 5
      RY(I)=YNEW
      RA(I)=ANEW

      ALL=.true.
      call ENERG(UTOT,ALL,X,Y,ANG,ANG2)
      NACCPT=NACCPT+1

      ! acumula promedios
    3 ACC(1)=ACC(1)+1.0D00    !GOTO 3
      ACC(2)=ACC(2)+UTOT
      USUBAV=USUBAV+UTOT
      XNTEST=DFLOAT(NCOUNT)
      
   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   !C        PRUEBA PARA OBTENER BARRA DE ERROR EN DENSIDAD
      IF (NERCONT.EQ.NSET) THEN
         WRITE(10,*)RHO, S, VOL
         NERCONT=0
      END IF
   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if ( NCOUNT==NMOVE/2 ) then
         open(unit=30,file="graficas/conf_intermedia.dat")
         do i=1,N
            write(30,*) RX(i),RY(i),S*AR,S,RA(I)
         end do
      end if

      IF (NCOUNT.EQ.NGOFR.AND.LGOFR) CALL GOFR
      IF (NCOUNT.GE.NMOVE) GOTO 6

      IF (NCOUNT.LT.NSUB) GOTO 1

      ! escribe resultados
      
      UAV=USUBAV/(XN*NCOUNT)
      !      WRITE(6,102) NCOUNT,DFLOAT(NACCPT)/DFLOAT(NCOUNT),RTEST/XNTEST,
      !     + UAV,VACCPT
      WRITE(6,102) NCOUNT,DFLOAT(NACCPT)/DFLOAT(NCOUNT),TAU,
     + UAV,VACCPT,RTEST/XNTEST
      write(20,*) NCOUNT, UAV
      WRITE(8,*)TAU
      
      NSUB=NSUB+NSUB0
      XXMOV = DFLOAT(NMOVE)
      XXX = XNTEST/XXMOV
      ! verifica si es necesario cambiar desplazamiento maximo
      !IF (DFLOAT(NACCPT)/DFLOAT(NCOUNT).GT.0.45) DISPL=DISPL*1.05
      !IF (DFLOAT(NACCPT)/DFLOAT(NCOUNT).LT.0.35) DISPL=DISPL*0.95
      ! verifica si es necesario cambiar CAMBIO DE VOLUMEN  maximo
      !IF (VACCPT.GT.0.45.AND.LVOL) SDISPL=SDISPL*1.05
      !IF (VACCPT.LT.0.35.AND.LVOL) SDISPL=SDISPL*0.95

      WRITE(2,*) XXX,UAV
      !WRITE(23,*) XN, NCOUNT
      
      GOTO 1
      ! fin de la corrida por part�cula


      ! marca clusters
      !! Revisa vecinos para clusters



    6 DO I=1,N          !GOTO 6
         DO J=1,N
            IF (J.NE.I) then
               X=RX(I)-RX(J)
               Y=RY(I)-RY(J)

               ! convencion de imagen m!nima
               IF (X.GT.0.5D00) THEN
                  X=X-1.0D00
               ELSE IF (X.LT.-0.5D00) THEN
                  X=X+1.0D00
               END IF

               IF (Y.GT.Y2) THEN
                  Y=Y-YC
               ELSE IF (Y.LT.-Y2) THEN
                  Y=Y+YC
               END IF
               RR=X*X+Y*Y

               if (RR.lt.(5*a)*(5*a)) then ! checa el rango espacial para formar clusters
                  CLU(I,J)=1.0D0
               else
                  CLU(I,J)=0.0D0
               end if

            end if
         END DO
      END DO

      do i=1,N
         write(16,*) RX(i),RY(i),S*AR,S,RA(I)
      end do


      !!!!!CHECAR
      I=1
      DO while (I.NE.0)
         CALL CLUSTERS(I) 
      end do


      ! inicia corrida por clusters
    7 I=INT(RAN2(ISEED)*N)+1     !!!!! GOTO 7
      NCOUNT=NCOUNT+1

      ! NERCONT ES EL CONTADOR QUE AL LLEGAR A NSET INDICA MUESTREAR SIGMA
      NERCONT=NERCONT+1
      RTEST=RTEST+RAN2(ISEED)
      ! desplaza la part!cula
      XNEW1(1:N)=RX(1:N)
      YNEW1(1:N)=RY(1:N)
      ANEW1(1:N)=RA(1:N)
      ! centro de masa
      K=1
      CX= XNEW1(I)
      CY=YNEW1(I)

      DO J=1,N
         IF (J.NE.I.AND.CLU(I,J).NE.0.0D0) THEN
            CX= CX+XNEW1(J)
            CY=CY+YNEW1(J)
            K=K+1
         end if
      END DO

      CX=CX/DFLOAT(K)
      CY=CY/DFLOAT(K)
      K=0

      IF (ABS(XNEW1(I)-CX)<0.4.AND.ABS(YNEW1(I)-CY)<0.4) then
         K=K+1
      end if

      DO J=1,N
         IF (J.NE.I.AND.CLU(I,J).NE.0.0D0) THEN
            IF (ABS(XNEW1(J)-CX)<0.4.AND.ABS(YNEW1(J)-CY)<0.4) then
               K=K+1
            end if
         end if
      END DO

      IF (K==0) GOTO 8

      ! giro en angulo
      DA=DISPLAng*(RAN2(ISEED)-0.5D0)
      ANEW1(I)=ANEW1(I)+DA
      X=XNEW1(I)-CX
      Y=YNEW1(I)-CY
      XNEW1(I)=X*COS(DA*PI/180)-Y*SIN(DA*PI/180)+CX
      YNEW1(I)=X*SIN(DA*PI/180)+Y*COS(DA*PI/180)+CY
      do J=1,N
         if (J.NE.I.AND.CLU(I,J).NE.0.0D0) then
            ANEW1(J)=ANEW1(J)+DA
            X=XNEW1(J)-CX
            Y=YNEW1(J)-CY
            XNEW1(J)=X*COS(DA*PI/180)-Y*SIN(DA*PI/180)+CX
            YNEW1(J)=X*SIN(DA*PI/180)+Y*COS(DA*PI/180)+CY
         end if
      end do

      ! desplazamiento  X, Y
    8 DX=DISPLClu*(RAN2(ISEED)-0.5D00)
      DY=DISPLClu*(RAN2(ISEED)-0.5D00)
      XNEW1(I)=XNEW1(I)+DX ! agrega desplazamientos a las componentes
      YNEW1(I)=YNEW1(I)+DY
     
      IF(XNEW1(I).GT.1.0D00) THEN
         XNEW1(I)=XNEW1(I)-1.0D00 
      ELSE
         IF (XNEW1(I).LT.0.0D00) THEN 
             XNEW1(I)=XNEW1(I)+1.0D00   
         END IF
      End if
      IF (YNEW1(I).GT.YC) THEN
         YNEW1(I)=YNEW1(I)-YC
      ELSE
         IF (YNEW1(I).LT.0.0D00) THEN
           YNEW1(I)=YNEW1(I)+YC
         END IF
      End if
      if (ANEW1(I).gt.180) then
         ANEW1(I)=ANEW1(I)-180
      end if
      if(ANEW1(I).lt.0) then
         ANEW1(I)=ANEW1(I)+180
      end if

      !!!!CHECAR
      do J=1,N
         if (J.NE.I.AND.CLU(I,J).NE.0.0D0) then ! si los valores diferentes de cero, (clu es un valor de 0 o 1)
            XNEW1(J)=XNEW1(J)+DX ! actualiza el valor de la nueva x sumandole un desplazamiento
            YNEW1(J)=YNEW1(J)+DY ! actualiza con desplazamiento de y 
            
            ! condicion periodica de frontera
            IF(XNEW1(J).GT.1.0D00) THEN 
               XNEW1(J)=XNEW1(J)-1.0D00 ! disminuye su valor si se excede 
            End if
            IF (XNEW1(J).LT.0.0D00) THEN
               XNEW1(J)=XNEW1(J)+1.0D00 ! aumenta el valor si le falta
            END IF
            IF (YNEW1(J).GT.YC) THEN
            YNEW1(J)=YNEW1(J)-YC 
            END IF
            IF (YNEW1(J).LT.0.0D00) THEN
            YNEW1(J)=YNEW1(J)+YC
            END IF
            if (ANEW1(J).gt.180) then
            ANEW1(J)=ANEW1(J)-180
            end if
            if(ANEW1(J).lt.0) then
            ANEW1(J)=ANEW1(J)+180
            end if
         end if
      end do

      ! revisa traslapes
      DO J=1,N
         IF (J.NE.I)then
            X=XNEW1(J)-XNEW1(I)
            Y=YNEW1(J)-YNEW1(I)

            ! convencion de imagen m!nima
            IF (X.GT.0.5D00) THEN
               X=X-1.0D00
            ELSE IF (X.LT.-0.5D00) THEN
               X=X+1.0D00
            END IF
            IF (Y.GT.Y2) THEN
               Y=Y-YC
            ELSE IF (Y.LT.-Y2) THEN
               Y=Y+YC
            END IF
            RR=X*X+Y*Y
            !!!!CONDICIONES DE TRASLAPE
            GGG = 2.00 +(AR-(1/AR))**2*(sin((ANEW1(J)-ANEW1(I))
     +       *PI/180))**2
            F1=-((X*cos(ANEW1(I)*PI/180)+ Y*sin(ANEW1(I)*PI/180))**2)
            F1=F1/(AR*S/2)**2
            F1=F1+1.00 + GGG

            F1=F1-((Y*cos(ANEW1(I)*PI/180)- X*sin(ANEW1(I)*PI/180))**2)
     +       /(SS/4)

            F2=-((X*cos(ANEW1(J)*PI/180)+ Y*sin(ANEW1(J)*PI/180))**2)
            F2=F2/(AR*S/2)**2
            F2=F2+1.00 +GGG
            F2=F2-((Y*cos(ANEW1(J)*PI/180)- X*sin(ANEW1(J)*PI/180))**2)
     +       /(SS/4)
            FI=4*(F1**2-3*F2)*(F2**2-3*F1)-(9-F1*F2)**2

            if (FI.lt.0.or.(FI.gt.0.and.F1>=0.and.F2>=0).or.RR.lt.SS) 
     +       GOTO 9
         end if
      end do

      do J=1,N
         IF(I.NE.J.AND.CLU(I,J).NE.0.0D0) then
            do M=1,N
               if (J.NE.M) then
                  X=XNEW1(M)-XNEW1(J)
                  Y=YNEW1(M)-YNEW1(J)

                  ! convencion de imagen m!nima
                  IF (X.GT.0.5D00) THEN
                     X=X-1.0D00
                  end if
                  
                  IF (X.LT.-0.5D00) THEN
                     X=X+1.0D00
                  END IF

                  IF (Y.GT.Y2) THEN
                     Y=Y-YC
                  end if

                  IF (Y.LT.-Y2) THEN
                     Y=Y+YC
                  END IF
                  RR=X*X+Y*Y

                  !!!!CONDICIONES DE TRASLAPE
                  GGG = 2.00 +(AR-(1/AR))**2*(sin((ANEW1(M)-ANEW1(J))
     +             *PI/180))**2
                  F1=-((X*cos(ANEW1(J)*PI/180)+ Y*sin(ANEW1(J)*PI/180))
     +             **2)
                  F1=F1/(AR*S/2)**2
                  F1=F1+1.00 + GGG
                        
                  F1=F1-((Y*cos(ANEW1(J)*PI/180)- X*sin(ANEW1(J)
     +             *PI/180))**2)/(SS/4)

                  F2=-((X*cos(ANEW1(M)*PI/180)+ Y*sin(ANEW1(M)*PI/180))
     +             **2)
                  F2=F2/(AR*S/2)**2
                  F2=F2+1.00 +GGG
                  F2=F2-((Y*cos(ANEW1(M)*PI/180)- X*sin(ANEW1(M)*PI
     +             /180))**2)/(SS/4)
                  FI=4*(F1**2-3*F2)*(F2**2-3*F1)-(9-F1*F2)**2

                  if (FI.lt.0.or.(FI.gt.0.and.F1>=0.and.F2>=0)
     +             .or.RR.lt.SS) GOTO 9
               end if
            end do
         end if
      end do

   ! actualiza posicion de la part!cula I
      RX(1:N)=XNEW1(1:N)
      RY(1:N)=YNEW1(1:N)
      RA(1:N)=ANEW1(1:N)

      ALL=.true.
      call ENERG(UTOT,ALL,X,Y,ANG,ANG2)

      !! Revisa vecinos para clusters
      DO I=1,N
         DO J=1,N
            IF (J.NE.I) then
               X=RX(I)-RX(J)
               Y=RY(I)-RY(J)

               ! convencion de imagen m!nima
               IF (X.GT.0.5D00) THEN
                  X=X-1.0D00
               ELSE IF (X.LT.-0.5D00) THEN
                  X=X+1.0D00
               END IF
               IF (Y.GT.Y2) THEN
                  Y=Y-YC
               ELSE IF (Y.LT.-Y2) THEN
                  Y=Y+YC
               END IF
               RR=X*X+Y*Y
               if (RR.lt.(5*a)*(5*a)) then
                  CLU(I,J)=1.0D0
               else
                  CLU(I,J)=0.0D0
               end if
            
            end if
         END DO
      END DO
      

      !! ubica clusters
      I=1
      DO while (I.NE.0)
         CALL CLUSTERS(I)
      end do
      NACCPT=NACCPT+1 ! acumulador de pasos 

      ! acumula promedios
    9 ACC(1)=ACC(1)+1.0D00   !!!GOTO 9
      ACC(2)=ACC(2)+UTOT
      USUBAV=USUBAV+UTOT ! lo utiliza para calcular el promedio 
      XNTEST=DFLOAT(NCOUNT)

      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !        PRUEBA PARA OBTENER BARRA DE ERROR EN DENSIDAD
      IF (NERCONT.EQ.NSET) THEN
         WRITE(10,*)RHO, S, VOL
         NERCONT=0
      END IF
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF (NCOUNT.EQ.NGOFR.AND.LGOFR) CALL GOFR 
      IF (NCOUNT.GE.NMOVE+NMOVE2) GOTO 12

      IF (NCOUNT.LT.NSUB) GOTO 7
      ! escribe resultados

      UAV=USUBAV/(XN*NCOUNT) ! aqui calcula el prom (mencionado en linea 891

      !      WRITE(6,102) NCOUNT,DFLOAT(NACCPT)/DFLOAT(NCOUNT),RTEST/XNTEST,
      !     + UAV,VACCPT
      WRITE(6,102) NCOUNT,DFLOAT(NACCPT)/DFLOAT(NCOUNT),TAU,
     + UAV,VACCPT,RTEST/XNTEST
      write(20,*) NCOUNT, UAV
      WRITE(8,*)TAU
      NSUB=NSUB+NSUB0 ! le agrega la parte inicial del contador 
      XXMOV = DFLOAT(NMOVE)  ! movimiento de en x 
      XXX = XNTEST/XXMOV ! razon entre los pasos en x contra el movimiento 
      ! verifica si es necesario cambiar desplazamiento maximo
      !IF (DFLOAT(NACCPT)/DFLOAT(NCOUNT).GT.0.45) DISPL=DISPL*1.05
      !IF (DFLOAT(NACCPT)/DFLOAT(NCOUNT).LT.0.35) DISPL=DISPL*0.95
      ! verifica si es necesario cambiar CAMBIO DE VOLUMEN  maximo
      !IF (VACCPT.GT.0.45.AND.LVOL) SDISPL=SDISPL*1.05
      !IF (VACCPT.LT.0.35.AND.LVOL) SDISPL=SDISPL*0.95
      WRITE(2,*)XXX,UAV
      !WRITE(23,*) XN, NCOUNT
      GOTO 7
      !! termina corrida por clusters


   12 UAV=USUBAV/(XN*NCOUNT)   !!!GOTO 12
      UAVRUN=ACC(2)/(XN*ACC(1))
      OPEN(UNIT=17,FILE='graficas/despues.dat',STATUS='unknown')
      DO I=1,N
         WRITE(17,*) RX(I),RY(I),S*AR,S,RA(I) ! escribe datos 
      END DO

      !WRITE(6,102) NCOUNT,DFLOAT(NACCPT)/DFLOAT(NCOUNT),RTEST/XNTEST,
      !     + UAV,VACCPT
      WRITE(6,102) NCOUNT,DFLOAT(NACCPT)/DFLOAT(NCOUNT),TAU,
     + UAV,VACCPT,RTEST/XNTEST
      write(20,*) NCOUNT, UAV
      WRITE(8,*)TAU
      WRITE(6,107)S
      WRITE(6,108)DISPL/S
      WRITE(6,109)SDISPL
      WRITE(6,103) ACC(1),ACC(3)
      WRITE(6,104) UAVRUN

      RETURN 

  101 FORMAT(/,1X,'inicio de corrrida',/,1X,'energ!a inicial',ES10.4,/)
C  102 FORMAT(1X,'NCOUNT=',I10,' NACCPT RATIO=',F8.4,' RTEST=',F8.4
C     + ,' UAV=',F10.4,' VACCPT RATIO=',F8.4,/)
  102 FORMAT(1X,'NCOUNT=',I10,' NACCPT RATIO=',F8.4,' TAU=',F12.8
     + ,' UAV=',ES12.4,' VACCPT RATIO=',F8.4,' random  RATIO=',F8.4,/)
  103 FORMAT(1X,'Num total de configuraciones desde el inicio',F12.0,/
     +,' N#mero de llamadas a g(r)',F12.0,/)
  104 FORMAT(1X,'Energ!a potencial promedio desde el inicio: ',ES12.4,/)

  107 FORMAT(1X,' S=',G16.8,/)
 108  FORMAT(1X,' DISPL=',G16.8,/)
 109  FORMAT(1X,' SDISPL=',G16.8,/)
      RETURN
      END

C     *********************************************
C     Identifica los clusters en la matriz
      SUBROUTINE CLUSTERS(COUNT)   !!todo bien
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NPART=2000,NACC=20,NG=30000)
      COMMON /BPOSITN/ RX(NPART),RY(NPART),RA(NPART),ACC(NACC),G(NG),AR,
     +                 CLU(NPART,NPART)
      COMMON /BCONSTR/ PI,ETA,RHO,XLAMBDA,XLAM2,XLAM3,SL2,SQL2,
     +                 SL3,SQL3,XN,TEMP,DISPL,DISPLAng, DISPLClu,XHISTG,
     +                 S,SS,SSLL,SL,XC,YC,YCINV,
     +                 XL,YL,Y2,RA3,A1, A2, A3, A4,a,b,tension,b_sigma
      COMMON /BCONSTG/ DR1,DR2,DR3,DR4,DR5,PHI,TAU
      COMMON /BCOSTRX/ SL4,SQL4,SL5,SQL5,XLAM4,XLAM5
      COMMON /BCOSTXX/ DR6,SL6,SQL6,XLAM6,E6
      COMMON /BCONSTI/ N,NGOFR,LGOFR,NMOVE,NMOVE2,NSUB,NGOFR0,ISEED
      COMMON /BCONSTV/ PRESS,VOL,SDISPL,NSET,LRHO
      LOGICAL LGOFR
      integer COUNT
      COUNT=0
      DO I=1, N
         DO J=1,N
            IF (I.NE.J.AND.CLU(I,J).NE.0.0D0) then ! si existe clusters 
               IF(CLU(J,I)==0.0D0) then
                  CLU(J,I)=1.0D0 ! cambia el valor 
                  COUNT=COUNT+1 ! y cuenta las iteraciones
               END IF
               DO M=1,N
                  IF (CLU(J,M).NE.0.0D0.AND.CLU(I,M)==0.0D0.AND.M.NE.J) 
     +             then  ! SOLAMENTE para valores del clu(IM) hacer 1 
                     CLU(I,M)=1.0D0
                     COUNT=COUNT+1 ! y contar 
                  END IF
               END DO
            END IF
         END DO
      END DO
      RETURN
      END
C     ***************************************************************
C     Calculo de energia potencial
      !! U(output)   :: aqui se almacena el valor de la energia 
      !!ALL(input)   :: Booleano que controla si queremos calcular la energía de todo el sistema
      !                 o solamente la energía de una partícula respecto a las demás

      !Las siguientes variables solamente se necesitan en caso de calcular la energía de una partícula respecto a las demás
      !!X(input)     :: Componente en x de la distancia entre centros   
      !!Y(input)     :: Componente en y de la distancia entre centros
      !!ANG(input)   :: Orientación 

      SUBROUTINE ENERG(U,ALL,X,Y,ANG,ANG2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NPART=2000,NACC=20,NG=30000)
      COMMON /BPOSITN/ RX(NPART),RY(NPART),RA(NPART),ACC(NACC),G(NG),AR,
     +                 CLU(NPART,NPART)
      COMMON /BCONSTR/ PI,ETA,RHO,XLAMBDA,XLAM2,XLAM3,SL2,SQL2,
     +                 SL3,SQL3,XN,TEMP,DISPL,DISPLAng, DISPLClu,XHISTG,
     +                 S,SS,SSLL,SL,XC,YC,YCINV,
     +                 XL,YL,Y2,RA3,A1, A2, A3, A4,a,b,tension,b_sigma
      COMMON /BCONSTG/ DR1,DR2,DR3,DR4,DR5,PHI,TAU
      COMMON /BCOSTRX/ SL4,SQL4,SL5,SQL5,XLAM4,XLAM5
      COMMON /BCOSTXX/ DR6,SL6,SQL6,XLAM6,E6
      COMMON /BCONSTI/ N,NGOFR,LGOFR,NMOVE,NMOVE2,NSUB,NGOFR0,ISEED
      COMMON /BCONSTV/ PRESS,VOL,SDISPL,NSET,LRHO
      LOGICAL LGOFR
      LOGICAL ALL 
      double Precision, PARAMETER::   kb=1.380649D-23

      open(unit=25,file="datos/bool.dat")
      write(25,*) ALL
      !Si ALL=TRUE calcula la energía de TODO el sistema
      If (ALL .eqv. .true.) then
         U=0
         DO 4 I=1,N-1
         DO 4 J=I+1,N
            IF (J.EQ.I) GOTO 4
            X=RX(J)-RX(I)
            Y=RY(J)-RY(I)
            ! Convencion de imagen m!nima
            IF (X.GT.0.5D00) THEN
               X=X-1.0D00
            ELSE IF (X.LT.-0.5D00) THEN
               X=X+1.0D00
            END IF
            IF (Y.GT.Y2) THEN
               Y=Y-YC
            ELSE IF (Y.LT.-Y2) THEN
               Y=Y+YC
            END IF
            RR=X*X+Y*Y

            !!Energia 
            if (RR .lt. 25*a*a ) then  !!CHECAR condición if
               dist=0
               ANGLE=RA(I)*PI/180 !orientación elipse 1 en radianes, medida respecto a la horizontal

               ANGLE2=RA(J)*PI/180 !orientación elipse 2 en radianes, medida respecto a la horizontal

               !Angulo entre elipses
               ANGLE3=datan( Y/X )
               call ellipses( a, b, a, b, ANGLE, ANGLE2, ANGLE3, dist )  !Calcula la distancia de máxima aproximación
                                                                         !los ángulos tienen que estar medidos respecto a la horizontal 
                                                                         !y en dirección antihoraria para que funcione la subrutina

               pppp= X*dcos(ANGLE)+Y*dsin(ANGLE) !producto punto entre el vector que une los centros y vector unitario de orientación de elipse1
               pppp=pppp/dsqrt(RR)
               ANGLE=acos(pppp) !angulo de orientación elipse1 respecto al vector que une los centros 

               pppp= X*dcos(ANGLE2)+Y*dsin(ANGLE2)
               pppp=pppp/dsqrt(RR)
               ANGLE2=acos(pppp)  !angulo de orientación elipse2 respecto al vector que une los centros

               U = U - A1*dcos(2*ANGLE2+2*ANGLE)*( b /                     !Potencial SIN UNIDADES
     +                             ( dsqrt(RR) - A2*dist+A3*b ) )**A4      !el potencial funciona con los ángulos tomados respecto al vector que une los centros
            end if
    4    CONTINUE
         
         U=U*b_sigma*b_sigma*tension !unidades SI

         U=U/(kb*TEMP)

      Else !!Calculamos energía a pares, con una partícula fija
         RR=X*X+Y*Y
         
         dist=0
         ANGLE=ANG*PI/180 !orientación elipse 1 en radianes, medida respecto a la horizontal

         ANGLE2=ANG2*PI/180 !orientación elipse 2 en radianes, medida respecto a la horizontal
         !Angulo entre elipses
         ANGLE3=datan( Y/X )
         call ellipses( a, b, a, b, ANGLE, ANGLE2, ANGLE3, dist )

         pppp= X*dcos(ANGLE)+Y*dsin(ANGLE)
         pppp=pppp/dsqrt(RR)

         ANGLE=acos(pppp)

         pppp= X*dcos(ANGLE2)+Y*dsin(ANGLE2)
         pppp=pppp/dsqrt(RR)

         ANGLE2=acos(pppp)

         U = U - A1*dcos(2*ANGLE2+2*ANGLE)*( b /                           !Potencial SIN UNIDADES
     +                              ( dsqrt(RR) - A2*dist+A3*b ) )**A4
      end if

      RETURN
      END
C     ************************************************************
C     Calcula g(r)

      SUBROUTINE GOFR !!todo bien

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NPART=2000,NACC=20,NG=30000)
      COMMON /BPOSITN/ RX(NPART),RY(NPART),RA(NPART),ACC(NACC),G(NG),AR,
     +                 CLU(NPART,NPART)
      COMMON /BCONSTR/ PI,ETA,RHO,XLAMBDA,XLAM2,XLAM3,SL2,SQL2,
     +                 SL3,SQL3,XN,TEMP,DISPL,DISPLAng, DISPLClu,XHISTG,
     +                 S,SS,SSLL,SL,XC,YC,YCINV,
     +                 XL,YL,Y2,RA3,A1, A2, A3, A4,a,b,tension,b_sigma
      COMMON /BCONSTG/ DR1,DR2,DR3,DR4,DR5,PHI,TAU
      COMMON /BCOSTRX/ SL4,SQL4,SL5,SQL5,XLAM4,XLAM5
      COMMON /BCOSTXX/ DR6,SL6,SQL6,XLAM6,E6
      COMMON /BCONSTI/ N,NGOFR,LGOFR,NMOVE,NMOVE2,NSUB,NGOFR0,ISEED
      COMMON /BCONSTV/ PRESS,VOL,SDISPL,NSET,LRHO
      LOGICAL LGOFR
      ACC(3)=ACC(3)+1.0D00 ! agrega una unidad al acumulador pero en su posicion 3 (nunca antes utilizada)
      RMAX=0.5D00*0.5D00 
      DO 20 I=1,N-1
      DO 20 J=I+1,N   
         X=RX(I)-RX(J)
         Y=RY(I)-RY(J)
         ! Convencion de imagen minima
         IF (X.GT.0.5D00) THEN
            X=X-1.0D00
         ELSE IF (X.LT.-0.5D00) THEN
            X=X+1.0D00
         END IF

         IF (Y.GT.Y2) THEN
            Y=Y-YC
         ELSE IF (Y.LT.-Y2) THEN
            Y=Y+YC
         END IF
         RR=X*X+Y*Y !termina CIM

         IF (RR.GT.RMAX) GOTO 20     

         R=SQRT(RR) !calcula la norma de r 
         DR1=(SL-S)/XHISTG  
         DR2=(SL2-SL)/XHISTG
         DR3=(SL3-SL2)/XHISTG
         DR4=(SL4-SL3)/XHISTG
         DR5=(SL5-SL4)/XHISTG
         DR6=(SL6-SL5)/XHISTG

         IF (R.LT.SL) THEN
            L=INT((R-S)/DR1)+1
         ELSE
            IF (R.LT.SL2) THEN
               L=INT(XHISTG+(R-SL)/DR2)+1
            ELSE
               IF (R.LT.SL3)THEN
                  L=INT(2.0*XHISTG+(R-SL2)/DR3)+1
               ELSE
                  IF (R.LT.SL4)THEN
                     L=INT(3.0*XHISTG+(R-SL3)/DR4)+1
                  ELSE
                     IF (R.LT.SL5)THEN
                        L=INT(4.0*XHISTG+(R-SL4)/DR5)+1
                     ELSE
                        L=INT(5.0*XHISTG+(R-SL5)/DR6)+1
                     END IF
                  END IF
               END IF
            END IF
         END IF

         if (L>30000) then
            WRITE(6,*) 'L=', L
         end if
         G(L)=G(L)+1.0D00

   20 CONTINUE
      NGOFR=NGOFR+NGOFR0
      RETURN
      END

      !     *****************************************************************
      !     Limpia para terminar

      SUBROUTINE FINISH   !!todo bien    NADAMAS ESCRIBE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NPART=2000,NACC=20,NG=30000)
      COMMON /BPOSITN/ RX(NPART),RY(NPART),RA(NPART),ACC(NACC),G(NG),AR,
     +                 CLU(NPART,NPART)
      COMMON /BCONSTR/ PI,ETA,RHO,XLAMBDA,XLAM2,XLAM3,SL2,SQL2,
     +                 SL3,SQL3,XN,TEMP,DISPL,DISPLAng, DISPLClu,XHISTG,
     +                 S,SS,SSLL,SL,XC,YC,YCINV,
     +                 XL,YL,Y2,RA3,A1, A2, A3, A4,a,b,tension,b_sigma
      COMMON /BCONSTG/ DR1,DR2,DR3,DR4,DR5,PHI,TAU
      COMMON /BCOSTRX/ SL4,SQL4,SL5,SQL5,XLAM4,XLAM5
      COMMON /BCOSTXX/ DR6,SL6,SQL6,XLAM6,E6
      COMMON /BCONSTI/ N,NGOFR,LGOFR,NMOVE,NMOVE2,NSUB,NGOFR0,ISEED
      COMMON /BCONSTV/ PRESS,VOL,SDISPL,NSET,LRHO
      LOGICAL LGOFR

      ! Escribe configuracion final en p2a.new
      OPEN(UNIT=4,FILE='datos/npt6.new',STATUS='unknown')
      WRITE(4,*) S
      DO 55 I=1,N
         WRITE(4,*) RX(I),RY(I),S*AR,S,RA(I)
   55 CONTINUE
      DO 66 I=1,NACC
         WRITE(4,*) ACC(I)
   66 CONTINUE
      DO 77 I=1,NG
         WRITE(4,*) G(I)
   77 CONTINUE
      !OPEN(UNIT=15,FILE='npt5.pot',STATUS='NEW')
      !WRITE(15,*) 1.0D0, E1
      !WRITE(15,*) XLAMBDA, E1
      !WRITE(15,*) XLAMBDA, E2
      !WRITE(15,*) XLAM2, E2
      !WRITE(15,*) XLAM2, E3
      !WRITE(15,*) XLAM3, E3
      !WRITE(15,*) XLAM3, E4
      !WRITE(15,*) XLAM4, E4
      !WRITE(15,*) XLAM4, E5
      !WRITE(15,*) XLAM5, E5
      !WRITE(15,*) XLAM5, E6
      !WRITE(15,*) XLAM6, E6
      !WRITE(15,*) XLAM6, 0.0D0
      RETURN
      END

      !*****************************************************************
      !Calcula el perfil
      
      SUBROUTINE RADIAL  !!todo bien
      !!!CHECAR que hace la subrutina

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NPART=2000,NACC=20,NG=30000)
      COMMON /BPOSITN/ RX(NPART),RY(NPART),RA(NPART),ACC(NACC),G(NG),AR,
     +                 CLU(NPART,NPART)
      COMMON /BCONSTR/ PI,ETA,RHO,XLAMBDA,XLAM2,XLAM3,SL2,SQL2,
     +                 SL3,SQL3,XN,TEMP,DISPL,DISPLAng, DISPLClu,XHISTG,
     +                 S,SS,SSLL,SL,XC,YC,YCINV,
     +                 XL,YL,Y2,RA3,A1, A2, A3, A4,a,b,tension,b_sigma
      COMMON /BCONSTG/ DR1,DR2,DR3,DR4,DR5,PHI,TAU
      COMMON /BCOSTRX/ SL4,SQL4,SL5,SQL5,XLAM4,XLAM5
      COMMON /BCOSTXX/ DR6,SL6,SQL6,XLAM6,E6
      COMMON /BCONSTI/ N,NGOFR,LGOFR,NMOVE,NMOVE2,NSUB,NGOFR0,ISEED
      COMMON /BCONSTV/ PRESS,VOL,SDISPL,NSET,LRHO
      LOGICAL LGOFR
      NHISTG=nint(XHISTG)
      XNAV=ACC(3)
      
      IF(XNAV.EQ.0.0D00)RETURN

      !DELR=(SL-S)*XL/XHISTG
      VOLME=XL*YL
      CONST=XN*XN*PI*XNAV/(VOLME*2.0D00)
      !NEND=INT((0.5D00*XC-S)*XHISTG/(SL-S))
      NEND=INT(5.0D00*XHISTG+((0.5D00*XC)-SL5)/DR6)

      OPEN(UNIT=7,FILE='datos/gdrnpt5.dat',STATUS='unknown')

      DO 1 L=1,NEND
         X=FLOAT(L)
         XHIST2=2.0D00*XHISTG
         XHIST3=3.0D00*XHISTG
         XHIST4=2*XHIST2
         XHIST5=XHIST2+XHIST3
         IF (X.LE.XHISTG) THEN
            DELR=DR1
            RB=((X-1.0D00)*DELR)*XL + 1
         ELSE
            IF (X.LE.XHIST2) THEN
               DELR=DR2
               RB=XLAMBDA + ((X-XHISTG -1.0D00)*DELR)*XL
            ELSE
               IF (X.LE.XHIST3) THEN
                  DELR=DR3
                  RB=XLAM2 + ((X-XHIST2 -1.0D00)*DELR)*XL
               ELSE
                  IF (X.LE.XHIST4) THEN
                     DELR=DR4
                     RB=XLAM3 + ((X-(XHIST3) -1.0D00)*DELR)*XL
                  ELSE
                     IF (X.LE.XHIST5) THEN
                        DELR=DR5
                        RB=XLAM4 + ((X-(XHIST4) -1.0D00)*DELR)*XL
                     ELSE
                        RB=XLAM5 + ((X-(XHIST5) -1.0D00)*DELR)*XL
                     END IF
                  END IF
               END IF
            END IF
         END IF

         R=RB + 0.5d00*DELR*XL
         !R=(SL-S)*(2.0D00*X-1.0D00)*XL/(XHISTG*2.0D00)
         !RB=(SL-S)*(X-1.0D00)*XL/XHISTG
         RSQ=R*R
         RBSQ=RB*RB
         RC=RB+DELR*XL
         RCSQ=RC*RC
         GR=G(L)/(CONST*(RCSQ-RBSQ))
         WRITE(7,*)R,GR,L
1     CONTINUE
      RETURN
      END

C     ****************************************************************


        FUNCTION RAN2(IDUM)
C       RAN2 OF NUMERICAL RECIPES 2ND ED.
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.0D00/IM1,
     &  IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     &  IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB)
        PARAMETER(EPS=1.2D-14,RNMX=1.0D00-EPS)
        DIMENSION IV(NTAB)
        SAVE IV,IY,IDUM2
        DATA IDUM2/123456789/,IV/NTAB*0/,IY/0/
        IF (IDUM.LE.0)THEN
C        WRITE(6,*)'INIT.',IDUM
          IDUM=MAX(-IDUM,1)
          IDUM2=IDUM
          DO J=NTAB+8,1,-1
            K=IDUM/IQ1
            IDUM=IA1*(IDUM-K*IQ1)-K*IR1
            IF(IDUM.LT.0)IDUM=IDUM+IM1
            IF(J.LE.NTAB)IV(J)=IDUM
          ENDDO
          IY=IV(1)
        ENDIF
        K=IDUM/IQ1
        IDUM=IA1*(IDUM-K*IQ1)-K*IR1
        IF(IDUM.LT.0)IDUM=IDUM+IM1
        K=IDUM2/IQ2
        IDUM2=IA2*(IDUM2-K*IQ2)-K*IR2
        IF(IDUM2.LT.0)IDUM2=IDUM2+IM2
        J=1+IY/NDIV
        IY=IV(J)-IDUM2
        IV(J)=IDUM
        IF(IY.LT.1)IY=IY+IMM1
        RAN2=MIN(AM*IY,RNMX)
        RETURN
        END


      !****************************************************************
      !Calcula la distancia de máximo acercamiento entre elipses
      !codigo obtenido de: http://www.math.kent.edu/~zheng/ellipse.html
      !INPUT:  a1, b1 (a1>b1): length of semiaxis of first ellipse
      !        a2, b2 (a2>b1): length of semiaxis of second ellipse
      !        theta1(0, 2*pi): angle associated with the major axis of first ellipse, meaured from lab frame counterclockwise
      !        theta2(0, 2*pi): angle associated with the major axis of second ellipse
      !        theta3(0, 2*pi): angle associated with the vector joining the centers, pointing to the center of the second ellipse

      !OUTPUT: dist: distance between the centers when two ellipses are externally tangent, i.e., distance of closest approach

        Subroutine ellipses(a1,b1,a2,b2,theta1,theta2,theta3,dist)
         Implicit None
         Double Complex,parameter::in=(1.d0,0.d0)
         Double Precision,Intent(in)::a1,b1,a2,b2,theta1,theta2,theta3
         Double Precision,Intent(out)::dist
         Double Precision::k1k2,k1d,k2d 
         Double Precision::e1,e2,eta,a11,a12,a22
         Double Precision::lambda1,lambda2,a2p,b2p
         Double Precision::kpmp,kpmps,t,deltap,A,B,C,D,E,alpha,beta,
     +    gamma,P,Q,Rc
         Double Precision::sn1,sn2,cs1,cs2,cs3,sn3
         Double Complex::U,y,qq
         Double Complex, External::ccbrt
     
         cs1=dcos(theta1);sn1=dsin(theta1)    ! components of the unit vector along the major axis of the first ellipse, k1
         cs2=dcos(theta2);sn2=dsin(theta2)    ! components of the unit vector along the major axis of the second ellipse, k2
         cs3=dcos(theta3);sn3=dsin(theta3)    ! components of the unit vector joining the centers, d
         k1d=dcos(theta3-theta1)              ! inner product of k1 and d
         k2d=dcos(theta3-theta2)              ! inner product of k2 and d
         k1k2=dcos(theta2-theta1)             ! inner product of k1 and k2
         !eccentricity of ellipses
         e1=1.d0-b1**2/a1**2
         e2=1.d0-b2**2/a2**2
         !component of A', matrix associated with the ellipse E2'
         eta=a1/b1-1.d0
         a11=b1**2/b2**2*(1.d0+0.5d0*(1.d0+k1k2)*(eta*(2.d0+eta)-e2*
     +    (1.d0+eta*k1k2)**2))
         a12=b1**2/b2**2*0.5d0*dsqrt(1.d0-k1k2**2)*(eta*(2.d0+eta)
     +    +e2*(1.d0-eta**2*k1k2**2))
         a22=b1**2/b2**2*(1.d0+0.5d0*(1.d0-k1k2)*(eta*(2.d0+eta)-e2*
     +    (1.d0-eta*k1k2)**2))
   
         !eigenvalues of A'
         lambda1=0.5d0*(a11+a22)+0.5d0*dsqrt((a11-a22)**2+4.d0*a12**2)  ! bigger one
         lambda2=0.5d0*(a11+a22)-0.5d0*dsqrt((a11-a22)**2+4.d0*a12**2)  ! smaller one
         !major and minor axes of transformed ellipse
         b2p=1.d0/dsqrt(lambda1)                                        ! length of minor axes
         a2p=1.d0/dsqrt(lambda2)                                        ! length of major axes
     
         deltap=a2p**2/b2p**2-1.d0
         If(dabs(k1k2) == 1.d0) then                                   ! if k1//k2, or k1//(-k2)
             if(a11>a22)then                                           ! if minor axes k+ is along k1  
                 kpmp=1.d0/dsqrt(1.d0-e1*k1d**2)*b1/a1*k1d                
                 kpmps=1.d0/dsqrt(1.d0-e1*k1d**2)*(sn3*cs1-cs3*sn1)
             else                                                      ! if major axes is k- along k1
                 kpmp=(sn3*cs1-cs3*sn1)/dsqrt(1.d0-e1*k1d**2)
                 kpmps=1.d0/dsqrt(1.d0-e1*k1d**2)*b1/a1*k1d
             end if
         Else                                                         ! if k1 is not parallel to k2
                 kpmp=(a12/dsqrt(1.d0+k1k2)*(b1/a1*k1d+k2d+(b1/a1-1.d0)
     +            *k1d*k1k2)+(lambda1-a11)/dsqrt(1.d0-k1k2)*(b1/a1*k1d-
     +             k2d-(b1/a1-1.d0)*k1d*k1k2))/dsqrt(2.d0*(a12**2+
     +             (lambda1-a11)**2)*(1.d0-e1*k1d**2))
                 kpmps=(-(lambda1-a11)/dsqrt(1.d0+k1k2)*(b1/a1*k1d+k2d+
     +            (b1/a1-1.d0)*k1d*k1k2)+a12/dsqrt(1.d0-k1k2)*(b1/a1*
     +            k1d-k2d-(b1/a1-1.d0)*k1d*k1k2))/dsqrt(2.d0*(a12**2+
     +            (lambda1-a11)**2)*(1.d0-e1*k1d**2))

         End If
         IF(kpmp==0.d0 .or. deltap==0.0d0) Then                       !if the major axes k- is along d' or they are both circles
             Rc=a2p+1.d0
         ELSE
         !coefficients of quartic for q
             t=1.d0/kpmp**2-1.d0
             A=-1.d0/b2p**2*(1.d0+t)
             B=-2.d0/b2p*(1.d0+t+deltap)
             C=-t-(1.d0+deltap)**2+1.d0/b2p**2*(1.d0+t+deltap*t)
             D=2.d0/b2p*(1.d0+t)*(1.d0+deltap)
             E=(1.d0+t+deltap)*(1.d0+deltap)
             
             !solution for quartic 
             alpha=-3.d0/8.d0*(B/A)**2+C/A
             beta=(B/A)**3.d0/8.d0-(B/A)*(C/A)/2.d0+D/A
             gamma=-3.d0/256.d0*(B/A)**4+C/A*(B/A)**2/16.d0-(B/A)*
     +        (D/A)/4.+E/A
         
             If(beta==0.d0) Then
                 qq=-B/4.d0/A+cdsqrt((-alpha+cdsqrt(alpha**2-4.d0*
     +            gamma*in))/2.)        
             Else
                 P=-alpha**2/12.d0-gamma
                 Q=-alpha**3/108.d0+gamma*alpha/3.d0-beta**2/8.d0
                 U=ccbrt(-0.5d0*Q+cdsqrt(Q**2/4.d0+P**3/27.d0*in))
     
                 if(abs(U)/=0.0d0) then
                     y=-5.d0/6d0*alpha+U-P/3.d0/U
                 else
                     y=-5.d0/6.d0*alpha-ccbrt(Q*in)
                 end if
                 
                 qq=-B/4.d0/A+0.5d0*(cdsqrt(alpha+2.d0*y)+cdsqrt(-
     +            (3.d0*alpha+2.d0*y+2.d0*beta/cdsqrt(alpha+2.d0*y))))
                 
         End If
     
         !substitute for R'
             Rc=dreal(cdsqrt((qq**2-1.d0)/deltap*(1.d0+b2p*(1.d0+deltap)
     +        /qq)**2+(1.d0-(qq**2-1.d0)/deltap)*(1.d0+b2p/qq)**2))
             
         END IF
         
           ! The distance of closest approach
     
         dist=Rc*b1/dsqrt(1.d0-e1*k1d**2)
         
      End Subroutine ellipses
      !*************************************
      !Calculate the cubic root of a complex number, return the principal value
      !*************************************
      Complex*16 function ccbrt(z)
      Implicit None
      Complex*16:: z
      Real*8::arg,theta,r
      Intrinsic::datan2
      arg=datan2(dimag(z),dreal(z))
      theta = arg / 3.0d0
      r = (dreal(z)**2+dimag(z)**2)**(1.d0/6.d0)
      ccbrt = dcmplx(r*dcos(theta), r*dsin(theta))
      End function ccbrt