        PROGRAM NPT6
C                                      FEBRERO 1999
                                      !Noviembre 2017
C
C     SIMULACION MONTE CARLO Elipses
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NPART=2000,NACC=20,NG=30000)
      COMMON /BPOSITN/ RX(NPART),RY(NPART),RA(NPART),ACC(NACC),G(NG),AR,
     +                 CLU(NPART,NPART)
      COMMON /BCONSTR/ PI,ETA,RHO,XLAMBDA,XLAM2,XLAM3,SL2,SQL2,
     +                 SL3,SQL3,XN,TEMP,DISPL,DISPLAng,XHISTG,E1,E2,E3,
     +                 S,SS,SSLL,SL,XC,YC,YCINV,
     +                 XL,YL,Y2,EA1,EA2,AAng, RA1,RA2,RA3
      COMMON /BCONSTG/ DR1,DR2,DR3,DR4,DR5,PHI,TAU
      COMMON /BCOSTRX/ SL4,SQL4,SL5,SQL5,XLAM4,XLAM5,E4,E5
      COMMON /BCOSTXX/ DR6,SL6,SQL6,XLAM6,E6
      COMMON /BCONSTI/ N,NGOFR,LGOFR,NMOVE,NMOVE2,NSUB,NGOFR0,ISEED
      COMMON /BCONSTV/ PRESS,VOL,SDISPL,NSET,LRHO
      LOGICAL LGOFR
      PI = 3.141592654D0
      OPEN(UNIT=6,FILE='npt6.dat',STATUS='unknown')
      OPEN(UNIT=8,FILE='npt6.tau',STATUS='unknown')
C     UNIDAD 10 NPT5.SIG GUARDA MUESTREO EN SIGMA
      OPEN(UNIT=10,FILE='npt6.sigma',STATUS='unknown')
      CALL START
      CALL MCARLO
      CALL FINISH
      CALL RADIAL
      CLOSE(UNIT=6)
      STOP
      END
C     ****************************************************************
      SUBROUTINE START

C       Lee variables de sistema y estado (rho,temp) de pozos.in
C       al inicio y pozos.old si es continuaci"n
C       Genera la configuraci"n inicial de fcc.
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NPART=2000,NACC=20,NG=30000)
      COMMON /BPOSITN/ RX(NPART),RY(NPART),RA(NPART),ACC(NACC),G(NG),AR,
     +                 CLU(NPART,NPART)
      COMMON /BCONSTR/ PI,ETA,RHO,XLAMBDA,XLAM2,XLAM3,SL2,SQL2,
     +                 SL3,SQL3,XN,TEMP,DISPL,DISPLAng,XHISTG,E1,E2,E3,
     +                 S,SS,SSLL,SL,XC,YC,YCINV,
     +                 XL,YL,Y2,EA1,EA2,AAng, RA1,RA2,RA3
      COMMON /BCONSTG/ DR1,DR2,DR3,DR4,DR5,PHI,TAU
      COMMON /BCOSTRX/ SL4,SQL4,SL5,SQL5,XLAM4,XLAM5,E4,E5
      COMMON /BCOSTXX/ DR6,SL6,SQL6,XLAM6,E6
      COMMON /BCONSTI/ N,NGOFR,LGOFR,NMOVE,NMOVE2,NSUB,NGOFR0,ISEED
      COMMON /BCONSTV/ PRESS,VOL,SDISPL,NSET,LRHO
      DIMENSION XA(6), YA(6)
      LOGICAL LGOFR
c      LOGICAL LRHO
C ACUMULADORES A CERO
      DO 3 I=1,NACC
    3 ACC(I)=0.0D00
      DO 4 I=1,NG
    4 G(I)=0.0D00
    
       CLU(:,:)=0.0D0
C LEE PARAMETROS DE ENTRADA: N,TEMP E.T.C.
      OPEN (UNIT=1,FILE='npt6.in',STATUS='unknown')
      READ(1,*) N
      READ(1,*) RHO
      READ(1,*) TEMP
      READ(1,*) PHI
      READ(1,*) XLAMBDA
      READ(1,*) XLAM2
      READ(1,*) XLAM3
      READ(1,*) XLAM4
      READ(1,*) XLAM5
      READ(1,*) XLAM6
      READ(1,*) E1
      READ(1,*) E2
      READ(1,*) E3
      READ(1,*) E4
      READ(1,*) E5
      READ(1,*) E6
      READ(1,*) EA1  !pozo en 90�, lado con lado
      READ(1,*) EA2  !pozo en 180�, punta con punta
      READ(1,*) NRUN,NMOVE,NMOVE2,NSUB
      READ(1,*) DISPL
      READ(1,*) DISPLAng
      READ(1,*) DISPLClu
      READ(1,*) AAng              !Ancho pozoz angulares
      READ(1,*) RA1      ! pozo angular LL actua de 1 a RA1 sigmas
      READ(1,*) RA2      !pozo angular PP actua de RA2 a RA3 sigmas
      READ(1,*) RA3
      READ(1,*) IAV
      READ(1,*) NAC,NCX,NCY
      READ(1,*) QX, QY
      READ(1,*) SDISPL
      READ(1,*) NGOFR,LGOFR
      READ(1,*) XHISTG
      READ(1,*)(XA(I),YA(I),I=1,NAC)
      READ(1,*)AR                       !!Aspect ratio
C LRHO DEFINE SI SE PARTE DE UNA DENSIDAD NUEVA
 !     READ(1,*) LRHO
C ESCRIBE PARAMETROS DE ENTRADA

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
      WRITE(6,*) 'E1=',E1
      WRITE(6,*) 'E2=',E2
      WRITE(6,*) 'E3=',E3
      WRITE(6,*) 'E4=',E4
      WRITE(6,*) 'E5=',E5
      WRITE(6,*) 'E6=',E6
      WRITE(6,*) 'PHI=',PHI
      WRITE(6,*)(XA(I),YA(I),I=1,NAC)
      IF (NRUN.EQ.0.0D00) WRITE(6,102)
      IF (NRUN.NE.0.0D00) WRITE(6,103)
      IF (IAV.EQ.0.0D00)  WRITE(6,104)
      IF (IAV.NE.0.0D00)  WRITE(6,105)
      IF (LGOFR)     WRITE(6,106) XHISTG
C CONVERSION A UNIDADES DE PROGRAMA (Lx = 1)
      PRESS=1.1547005384D0*PHI
      XN=DFLOAT(N)
C SE CALCULA EL TAMA~NO DE LA MUESTRA PARA BARRAS DE ERROR EN SIGMA
      XNMOVE=DFLOAT(NMOVE)
      XSET=XNMOVE/100D0
      NSET=NINT(XSET)
C CADA NSET MOVIMIENTOS SE MUESTREA SIGMA
      XNCX=DFLOAT(NCX)
      XNCY=DFLOAT(NCY)
      XC=(XNCX*QX)/(XNCX*QX)
      YC=(XNCY*QY)/(XNCX*QX)
      Y2=YC/2.0D00
      YCINV=2.0D00/YC
      IF (NRUN.EQ.0) THEN
C CALCULA SIGMA EN UNIDADES DE CAJA
       ETA=RHO*PI/4.0D00
       XL=(XN/(RHO*YC))**(1.0D00/2.0D00)
       YL=XL*YC
       S=1.0D00/XL
C  CONSTRUYE LATIZ INICIAL CON SIGMA DADA
             ISEED=-123456789

      DO I=1,N
    9 RX(I)=RAN2(ISEED)
      RY(I)=YC*RAN2(ISEED)
      RA(I)=180*RAN2(ISEED)
      DO J=1,I
      IF (J.NE.I)then
      X=RX(J)-RX(I)
      Y=RY(J)-RY(I)

C convenci"n de imagen m!nima
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
       GGG = 2.00 +(AR-(1/AR))**2*(sin((RA(J)-RA(I))*PI/180))**2
       F1=-((X*cos(RA(I)*PI/180)+ Y*sin(RA(I)*PI/180))**2)
       F1=F1/(AR*S/2)**2
       F1=F1+1.00 + GGG

      F1=F1-((Y*cos(RA(I)*PI/180)- X*sin(RA(I)*PI/180))**2)/(S*S/4)

       F2=-((X*cos(RA(J)*PI/180)+ Y*sin(RA(J)*PI/180))**2)
       F2=F2/(AR*S/2)**2
       F2=F2+1.00 +GGG
      F2=F2-((Y*cos(RA(J)*PI/180)- X*sin(RA(J)*PI/180))**2)/(S*S/4)
      FI=4*(F1**2-3*F2)*(F2**2-3*F1)-(9-F1*F2)**2

      if (FI.lt.0.or.(FI.gt.0.and.F1>=0.and.F2>=0).or.RR.lt.S*S) GOTO 9
      end if
      end do
      end do


        ELSE
        OPEN(UNIT=3,FILE='npt6.old',STATUS='OLD')
        IF (LRHO.EQ.1) THEN
C CALCULA SIGMA CON DENSIDAD DADA Y CONFIGURACION VIEJA
            ETA=RHO*PI/4.0D0
            XL=(XN/(RHO*YC))**(1.0D00/2.0D00)
            YL=XL*YC
            S=1.0D00/XL
            READ(3,*) RX(1)
            ELSE
C   LEE VALOR DE SIGMA
               READ(3,*) S
               VOL=YC/(S*S)
               RHO=XN/VOL
               ETA=RHO*PI/4.0D00
               XL=1/S
               YL=XL*YC
           END IF
C     LEE CONFIGURACION VIEJA
       DO 11 I=1,N
C       WRITE(3,*) RX(I),RY(I), RA(I)
       READ(3,*) RX(I),RY(I),Ancho,S, RA(I)
   11  CONTINUE
C lee acumuladores viejos
       IF (IAV.EQ.0) GOTO 2
       DO 22 I=1,NACC
       READ(3,*)ACC(I)
   22  CONTINUE
       DO 33 I=1,NG
       READ(3,*)G(I)
   33  CONTINUE
C escribe acumulador y configuraci"n de inicio o re-inicio
       WRITE(6,111) ACC(1)
       WRITE(6,112)
    2  END IF
C     CONTINUA CONVIRTIENDO A UNIDADES DE PROGRAMA
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
C ESCRIBE LAS UNIDADES CONVERTIDAS
      WRITE(6,107) RHO
      WRITE(6,108) XC,YC
      WRITE(6,109) XL,YL
      WRITE(6,110) S
      DO 44 I=1,20
      WRITE(6,115) RX(I),RY(I), RA(I)
   44 CONTINUE

C formatos de escritura de START
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

C ******************************************************************

      SUBROUTINE MCARLO
C CALCULA LA ENERGIA POTENCIAL DE LA CONFIGURACION INICIAL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NPART=2000,NACC=20,NG=30000)
      COMMON /BPOSITN/ RX(NPART),RY(NPART),RA(NPART),ACC(NACC),G(NG),AR,
     +                 CLU(NPART,NPART)
      COMMON /BCONSTR/ PI,ETA,RHO,XLAMBDA,XLAM2,XLAM3,SL2,SQL2,
     +                 SL3,SQL3,XN,TEMP,DISPL,DISPLAng,XHISTG,E1,E2,E3,
     +                 S,SS,SSLL,SL,XC,YC,YCINV,
     +                 XL,YL,Y2,EA1,EA2,AAng, RA1,RA2,RA3
      COMMON /BCONSTG/ DR1,DR2,DR3,DR4,DR5,PHI,TAU
      COMMON /BCOSTRX/ SL4,SQL4,SL5,SQL5,XLAM4,XLAM5,E4,E5
      COMMON /BCOSTXX/ DR6,SL6,SQL6,XLAM6,E6
      COMMON /BCONSTI/ N,NGOFR,LGOFR,NMOVE,NMOVE2,NSUB,NGOFR0,ISEED
      COMMON /BCONSTV/ PRESS,VOL,SDISPL,NSET,LRHO
      LOGICAL LGOFR
      DOUBLE PRECISION XNEW1(N), YNEW1(N), ANEW1(N)
C inicializa secuencia al azar
      ISEED=-123456789
      CALL ENERG(UTOT)
      UPERP=UTOT/XN
C     NTEST=0
      USUBAV=0.0D00
      NCOUNT=0

      NERCONT=0
      NACCPT=0
      RTEST=0.0D00
      NSUB0=NSUB
      NGOFR0=NGOFR
      WRITE(6,101) UPERP
C comienza la secuencia de Monte Carlo

      OPEN(UNIT=2,FILE='unpt6.dat',STATUS='NEW')
      !OPEN(UNIT=93,FILE='angulos.dat',STATUS='NEW')

C escoge part!cula al azar
    1 I=INT(RAN2(ISEED)*N)+1
      NCOUNT=NCOUNT+1

C     NERCONT ES EL CONTADOR QUE AL LLEGAR A NSET INDICA MUESTREAR SIGMA
      NERCONT=NERCONT+1
      RTEST=RTEST+RAN2(ISEED)
C desplaza la part!cula
      XNEW=RX(I)+DISPL*(RAN2(ISEED)-0.5D00)
      YNEW=RY(I)+DISPL*(RAN2(ISEED)-0.5D00)
       !!!angular
      ANEW=RA(I)+DISPLAng*(RAN2(ISEED)-0.5D0)
C condici"n peri"dica de frontera
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
      
      if (ANEW.gt.180) then
         ANEW=ANEW-180
      end if
      if(ANEW.lt.0) then
         ANEW=ANEW+180
      end if
C elimina part!culas traslapadas
C y calcula la nueva energ!a (provisional)
      UNEW=0.0D00
      DO 2 J=1,N
      IF (J.EQ.I) GOTO 2
      X=RX(J)-XNEW
      Y=RY(J)-YNEW
C convenci"n de imagen m!nima
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
       F2=F2-((X*cos(RA(J)*PI/180)+ Y*sin(RA(J)*PI/180))**2)/(AR*S/2)**2
       F2=F2-((Y*cos(RA(J)*PI/180)- X*sin(RA(J)*PI/180))**2)/(SS/4)
       FI=4*(F1**2-3*F2)*(F2**2-3*F1)-(9-F1*F2)**2

      !FI=0 es tangencia
      !no puntos en com�n es FI positiva y al menos uno de f1, f2 negativos
      
      if (FI.lt.0.or.(FI.gt.0.and.F1>=0.and.F2>=0).or.RR.lt.SS) GOTO 3

       
       !!Energia pozos angulares
      if (RR.lt.RA3*RA3*SS )then
          pppp= X*cos(ANEW*PI/180)+Y*sin(ANEW*PI/180)
          pppp=pppp/sqrt(RR)
          ANGLE=ACOS(pppp)
          ANGLE=ANGLE*180/PI
          pppp= X*cos(RA(J)*PI/180)+Y*sin(RA(J)*PI/180)
          pppp=pppp/sqrt(RR)
          ANGLE2=ACOS(pppp)
          ANGLE2=ANGLE2*180/PI
      if (ANGLE>(90-AAng/2).and.ANGLE<(90+AAng/2).and.RR<RA1*RA1*SS)then
          if (ANGLE2>(90-AAng/2).and.ANGLE2<(90+AAng/2) )then
             UNEW=UNEW+EA1
          end if
      end if
       if (ANGLE>(180-AAng/2).or.ANGLE<AAng/2.and.RR>RA2*RA2*SS)then
          if (ANGLE2>(180-AAng/2).or.ANGLE2<AAng/2) then
             UNEW=UNEW+EA2
          end if
          end if
      end if
      
      IF (RR.LT.SSLL) THEN

         UNEW=UNEW+ E1+(sqrt(RR)-S)*(E2-E1)/(SL-S)

         ELSE
         IF(RR.LT.SQL2) THEN
         
            UNEW=UNEW+ E2+(sqrt(RR)-SL)*(E3-E2)/(SL2-SL)
            
            ELSE
            IF (RR.LT.SQL3) THEN
               UNEW=UNEW+E3
               ELSE
               IF (RR.LT.SQL4) THEN
                  UNEW=UNEW+E4
                  ELSE
                  IF (RR.LT.SQL5) THEN
                     UNEW=UNEW+E5
                     ELSE
                        IF (RR.LT.SQL6) UNEW=UNEW+E6
                  END IF
               END IF
            END IF
         END IF
      END IF
    2 CONTINUE
      UOLD=0.0D00
      DO 4 J=1,N
      IF (J.EQ.I) GOTO 4
      X=RX(J)-RX(I)
      Y=RY(J)-RY(I)
      
C convenci"n de imagen m!nima
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
       GGG = 2.00 +(AR-(1/AR))**2*(sin((RA(J)-RA(I))*PI/180))**2
       F1=1.00 + GGG
       F1=F1-((X*cos(RA(I)*PI/180)+ Y*sin(RA(I)*PI/180))**2)/(AR*S/2)**2
       F1=F1-((Y*cos(RA(I)*PI/180)- X*sin(RA(I)*PI/180))**2)/(SS/4)
       F2=1.00 +GGG
       F2=F2-((X*cos(RA(J)*PI/180)+ Y*sin(RA(J)*PI/180))**2)/(AR*S/2)**2
       F2=F2-((Y*cos(RA(J)*PI/180)- X*sin(RA(J)*PI/180))**2)/(SS/4)
       FI=4*(F1**2-3*F2)*(F2**2-3*F1)-(9-F1*F2)**2

      if (FI.lt.0.or.(FI.gt.0.and.F1>=0.and.F2>=0).or.RR.lt.SS) GOTO 42



      !!Energia pozos angulares
       if (RR.lt.RA3*RA3*SS) then
          pppp= X*cos(RA(I)*PI/180)+Y*sin(RA(I)*PI/180)
          pppp=pppp/sqrt(RR)
          ANGLE=ACOS(pppp)
          ANGLE=ANGLE*180/PI
          pppp= X*cos(RA(I)*PI/180)+Y*sin(RA(I)*PI/180)
          pppp=pppp/sqrt(RR)
          ANGLE2=ACOS(pppp)
          ANGLE2=ANGLE2*180/PI
      if (ANGLE>(90-AAng/2).and.ANGLE<(90+AAng/2).and.RR<RA1*RA1*SS)then
         if (ANGLE2>(90-AAng/2).and.ANGLE2<(90+AAng/2) )then
             UOLD=UOLD+EA1
           end if
          end if
       if (ANGLE>(180-AAng/2).or.ANGLE<AAng/2.and.RR>RA2*RA2*SS )then
          if (ANGLE2>(180-AAng/2).or.ANGLE2<AAng/2 )then
             UOLD=UOLD+EA2
             end if
          end if
       end if
       
      
      IF (RR.LT.SSLL) THEN
         UOLD=UOLD+E1+(sqrt(RR)-S)*(E2-E1)/(SL-S)
         ELSE
         IF(RR.LT.SQL2) THEN
            UOLD=UOLD+E2+(sqrt(RR)-SL)*(E3-E2)/(SL2-SL)
            ELSE
            IF (RR.LT.SQL3) THEN
               UOLD=UOLD+E3
               ELSE
               IF (RR.LT.SQL4) THEN
                  UOLD=UOLD+E4
                  ELSE
                  IF (RR.LT.SQL5) THEN
                     UOLD=UOLD+E5
                     ELSE
                        IF (RR.LT.SQL6) UOLD=UOLD+E6
                  END IF
               END IF
            END IF
         END IF
      END IF
    4 CONTINUE
      DENERG=UNEW-UOLD
      
      IF (DENERG.LE.0.0D00) GOTO 5
C compara el factor de BOLTZMANN con n#mero al azar
      RND=RAN2(ISEED)
      IF (RND.GT.EXP(-DENERG/TEMP)) GOTO 3
C actualiza posici"n de la part!cula I
    5 RX(I)=XNEW
      RY(I)=YNEW
      RA(I)=ANEW
      call ENERG(UTOT)
      NACCPT=NACCPT+1.0D00
C acumula promedios
    3 ACC(1)=ACC(1)+1.0D00
      ACC(2)=ACC(2)+UTOT
      USUBAV=USUBAV+UTOT
      XNTEST=DFLOAT(NCOUNT)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        PRUEBA PARA OBTENER BARRA DE ERROR EN DENSIDAD
           IF (NERCONT.EQ.NSET) THEN
             WRITE(10,*)RHO, S, VOL
             NERCONT=0
           END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (NCOUNT.EQ.NGOFR.AND.LGOFR) CALL GOFR
      IF (NCOUNT.GE.NMOVE) GOTO 6

      IF (NCOUNT.LT.NSUB) GOTO 1
C escribe resultados

      UAV=USUBAV/(XN*NCOUNT)

c      WRITE(6,102) NCOUNT,DFLOAT(NACCPT)/DFLOAT(NCOUNT),RTEST/XNTEST,
c     + UAV,VACCPT
      WRITE(6,102) NCOUNT,DFLOAT(NACCPT)/DFLOAT(NCOUNT),TAU,
     + UAV,VACCPT,RTEST/XNTEST
      WRITE(8,*)TAU
      NSUB=NSUB+NSUB0
      XXMOV = DFLOAT(NMOVE)
      XXX = XNTEST/XXMOV
C verifica si es necesario cambiar desplazamiento maximo
      !IF (DFLOAT(NACCPT)/DFLOAT(NCOUNT).GT.0.45) DISPL=DISPL*1.05
      !IF (DFLOAT(NACCPT)/DFLOAT(NCOUNT).LT.0.35) DISPL=DISPL*0.95
C verifica si es necesario cambiar CAMBIO DE VOLUMEN  maximo
      !IF (VACCPT.GT.0.45.AND.LVOL) SDISPL=SDISPL*1.05
      !IF (VACCPT.LT.0.35.AND.LVOL) SDISPL=SDISPL*0.95
      WRITE(2,*)XXX,UAV
      !WRITE(23,*) XN, NCOUNT
      GOTO 1

C fin de la corrida por part�cula
C marca clusters

      !! Revisa vecinos para clusters
    6  DO I=1,N
      DO J=1,N
      IF (J.NE.I) then
      X=RX(I)-RX(J)
      Y=RY(I)-RY(J)

C convenci"n de imagen m!nima
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
      if (RR.lt.(RA3+3)*(RA3+3)*SS) then
          CLU(I,J)=1.0D0
       else
           CLU(I,J)=0.0D0
      end if

      end if
      END DO
      END DO
       I=1
       DO while (I.NE.0)
          CALL CLUSTERS(I)
       end do
      OPEN(UNIT=15,FILE='antes.dat',STATUS='unknown')
      DO I=1,N
      WRITE(15,*)RX(I),RY(I),S*AR,S,RA(I)
      END DO

C inicia corrida por clusters
    7 I=INT(RAN2(ISEED)*N)+1
      NCOUNT=NCOUNT+1

C     NERCONT ES EL CONTADOR QUE AL LLEGAR A NSET INDICA MUESTREAR SIGMA
      NERCONT=NERCONT+1
      RTEST=RTEST+RAN2(ISEED)
C desplaza la part!cula
      XNEW1(1:N)=RX(1:N)
      YNEW1(1:N)=RY(1:N)
      ANEW1(1:N)=RA(1:N)
C centro de masa
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
C giro en angulo
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
C desplazamiento  X, Y
    8 DX=DISPLClu*(RAN2(ISEED)-0.5D00)
      DY=DISPLClu*(RAN2(ISEED)-0.5D00)
      XNEW1(I)=XNEW1(I)+DX
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
      do J=1,N
         if (J.NE.I.AND.CLU(I,J).NE.0.0D0) then
             XNEW1(J)=XNEW1(J)+DX
             YNEW1(J)=YNEW1(J)+DY
C condici"n peri"dica de frontera
             IF(XNEW1(J).GT.1.0D00) THEN
                XNEW1(J)=XNEW1(J)-1.0D00
             ENd if
             IF (XNEW1(J).LT.0.0D00) THEN
                XNEW1(J)=XNEW1(J)+1.0D00
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
C revisa traslapes
      DO J=1,N
      IF (J.NE.I)then
      X=XNEW1(J)-XNEW1(I)
      Y=YNEW1(J)-YNEW1(I)

C convenci"n de imagen m!nima
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
       GGG = 2.00 +(AR-(1/AR))**2*(sin((ANEW1(J)-ANEW1(I))*PI/180))**2
       F1=-((X*cos(ANEW1(I)*PI/180)+ Y*sin(ANEW1(I)*PI/180))**2)
       F1=F1/(AR*S/2)**2
       F1=F1+1.00 + GGG

      F1=F1-((Y*cos(ANEW1(I)*PI/180)- X*sin(ANEW1(I)*PI/180))**2)/(SS/4)

       F2=-((X*cos(ANEW1(J)*PI/180)+ Y*sin(ANEW1(J)*PI/180))**2)
       F2=F2/(AR*S/2)**2
       F2=F2+1.00 +GGG
      F2=F2-((Y*cos(ANEW1(J)*PI/180)- X*sin(ANEW1(J)*PI/180))**2)/(SS/4)
      FI=4*(F1**2-3*F2)*(F2**2-3*F1)-(9-F1*F2)**2

      if (FI.lt.0.or.(FI.gt.0.and.F1>=0.and.F2>=0).or.RR.lt.SS) GOTO 9
      end if
      end do
      
      do J=1,N
         IF(I.NE.J.AND.CLU(I,J).NE.0.0D0) then
             do M=1,N
              if (J.NE.M) then
                X=XNEW1(M)-XNEW1(J)
                Y=YNEW1(M)-YNEW1(J)

C convenci"n de imagen m!nima
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
       GGG = 2.00 +(AR-(1/AR))**2*(sin((ANEW1(M)-ANEW1(J))*PI/180))**2
       F1=-((X*cos(ANEW1(J)*PI/180)+ Y*sin(ANEW1(J)*PI/180))**2)
       F1=F1/(AR*S/2)**2
       F1=F1+1.00 + GGG

      F1=F1-((Y*cos(ANEW1(J)*PI/180)- X*sin(ANEW1(J)*PI/180))**2)/(SS/4)

       F2=-((X*cos(ANEW1(M)*PI/180)+ Y*sin(ANEW1(M)*PI/180))**2)
       F2=F2/(AR*S/2)**2
       F2=F2+1.00 +GGG
      F2=F2-((Y*cos(ANEW1(M)*PI/180)- X*sin(ANEW1(M)*PI/180))**2)/(SS/4)
      FI=4*(F1**2-3*F2)*(F2**2-3*F1)-(9-F1*F2)**2

      if (FI.lt.0.or.(FI.gt.0.and.F1>=0.and.F2>=0).or.RR.lt.SS) GOTO 9
             end if
             end do
         end if
      end do

C actualiza posici"n de la part!cula I
   11 RX(1:N)=XNEW1(1:N)
      RY(1:N)=YNEW1(1:N)
      RA(1:N)=ANEW1(1:N)


      call ENERG(UTOT)
      !! Revisa vecinos para clusters
      DO I=1,N
      DO J=1,N
      IF (J.NE.I) then
      X=RX(I)-RX(J)
      Y=RY(I)-RY(J)

C convenci"n de imagen m!nima
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
      if (RR.lt.(RA3+3)*(RA3+3)*SS) then
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
      NACCPT=NACCPT+1.0D00
C acumula promedios
    9 ACC(1)=ACC(1)+1.0D00
      ACC(2)=ACC(2)+UTOT
      USUBAV=USUBAV+UTOT
      XNTEST=DFLOAT(NCOUNT)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        PRUEBA PARA OBTENER BARRA DE ERROR EN DENSIDAD
           IF (NERCONT.EQ.NSET) THEN
             WRITE(10,*)RHO, S, VOL
             NERCONT=0
           END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (NCOUNT.EQ.NGOFR.AND.LGOFR) CALL GOFR
      IF (NCOUNT.GE.NMOVE+NMOVE2) GOTO 12

      IF (NCOUNT.LT.NSUB) GOTO 7
C escribe resultados

      UAV=USUBAV/(XN*NCOUNT)

c      WRITE(6,102) NCOUNT,DFLOAT(NACCPT)/DFLOAT(NCOUNT),RTEST/XNTEST,
c     + UAV,VACCPT
      WRITE(6,102) NCOUNT,DFLOAT(NACCPT)/DFLOAT(NCOUNT),TAU,
     + UAV,VACCPT,RTEST/XNTEST
      WRITE(8,*)TAU
      NSUB=NSUB+NSUB0
      XXMOV = DFLOAT(NMOVE)
      XXX = XNTEST/XXMOV
C verifica si es necesario cambiar desplazamiento maximo
      !IF (DFLOAT(NACCPT)/DFLOAT(NCOUNT).GT.0.45) DISPL=DISPL*1.05
      !IF (DFLOAT(NACCPT)/DFLOAT(NCOUNT).LT.0.35) DISPL=DISPL*0.95
C verifica si es necesario cambiar CAMBIO DE VOLUMEN  maximo
      !IF (VACCPT.GT.0.45.AND.LVOL) SDISPL=SDISPL*1.05
      !IF (VACCPT.LT.0.35.AND.LVOL) SDISPL=SDISPL*0.95
      WRITE(2,*)XXX,UAV
      !WRITE(23,*) XN, NCOUNT
      GOTO 7
C termina corrida por clusters
   12  UAV=USUBAV/(XN*NCOUNT)
       UAVRUN=ACC(2)/(XN*ACC(1))
       OPEN(UNIT=17,FILE='despues.dat',STATUS='unknown')
      DO I=1,N
      WRITE(17,*)RX(I),RY(I),S*AR,S,RA(I)
      END DO

c      WRITE(6,102) NCOUNT,DFLOAT(NACCPT)/DFLOAT(NCOUNT),RTEST/XNTEST,
c     + UAV,VACCPT
      WRITE(6,102) NCOUNT,DFLOAT(NACCPT)/DFLOAT(NCOUNT),TAU,
     + UAV,VACCPT,RTEST/XNTEST
      WRITE(8,*)TAU
      WRITE(6,107)S
      WRITE(6,108)DISPL/S
      WRITE(6,109)SDISPL
       WRITE(6,103) ACC(1),ACC(3)
       WRITE(6,104) UAVRUN
       RETURN
C part!culas traslapadas
 42   WRITE(6,105) I,J,RR
      WRITE(6,106) SQRT(RR),S
      WRITE(6,*) RX(I),RY(I)
      WRITE(6,*) RX(J),RY(J)
  101 FORMAT(/,1X,'inicio de corrrida',/,1X,'energ!a inicial',G16.8,/)
C  102 FORMAT(1X,'NCOUNT=',I10,' NACCPT RATIO=',F8.4,' RTEST=',F8.4
C     + ,' UAV=',F10.4,' VACCPT RATIO=',F8.4,/)
  102 FORMAT(1X,'NCOUNT=',I10,' NACCPT RATIO=',F8.4,' TAU=',F12.8
     + ,' UAV=',F10.4,' VACCPT RATIO=',F8.4,' random  RATIO=',F8.4,/)
  103 FORMAT(1X,'Num total de configuraciones desde el inicio',F12.0,/
     +,' N#mero de llamadas a g(r)',F12.0,/)
  104 FORMAT(1X,'Energ!a potencial promedio desde el inicio',F10.4,/)
  105 FORMAT(1X,'*****part!culas traslapadas*****',' I=',I4,' J=',I4,
     + ' R**2=',G16.8,/)
  106 FORMAT(1X,'R=',G16.8,' S=',G16.8,/)
  107 FORMAT(1X,' S=',G16.8,/)
 108  FORMAT(1X,' DISPL=',G16.8,/)
 109  FORMAT(1X,' SDISPL=',G16.8,/)
      RETURN
      END

C     *********************************************
C      Identifica los clusters en la matriz
       SUBROUTINE CLUSTERS(COUNT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NPART=2000,NACC=20,NG=30000)
      COMMON /BPOSITN/ RX(NPART),RY(NPART),RA(NPART),ACC(NACC),G(NG),AR,
     +                 CLU(NPART,NPART)
      COMMON /BCONSTR/ PI,ETA,RHO,XLAMBDA,XLAM2,XLAM3,SL2,SQL2,
     +                 SL3,SQL3,XN,TEMP,DISPL,DISPLAng,XHISTG,E1,E2,E3,
     +                 S,SS,SSLL,SL,XC,YC,YCINV,
     +                 XL,YL,Y2,EA1,EA2,AAng, RA1,RA2,RA3
      COMMON /BCONSTG/ DR1,DR2,DR3,DR4,DR5,PHI,TAU
      COMMON /BCOSTRX/ SL4,SQL4,SL5,SQL5,XLAM4,XLAM5,E4,E5
      COMMON /BCOSTXX/ DR6,SL6,SQL6,XLAM6,E6
      COMMON /BCONSTI/ N,NGOFR,LGOFR,NMOVE,NMOVE2,NSUB,NGOFR0,ISEED
      COMMON /BCONSTV/ PRESS,VOL,SDISPL,NSET,LRHO
      LOGICAL LGOFR
      integer COUNT
       COUNT=0
       DO I=1, N
          DO J=1,N
             IF (I.NE.J.AND.CLU(I,J).NE.0.0D0) then
                IF(CLU(J,I)==0.0D0) then
                    CLU(J,I)=1.0D0
                    COUNT=COUNT+1
                END IF
                DO M=1,N
        IF (CLU(J,M).NE.0.0D0.AND.CLU(I,M)==0.0D0.AND.M.NE.J) then
                       CLU(I,M)=1.0D0
                       COUNT=COUNT+1
                    END IF
                END DO
             END IF
           END DO
       END DO
       RETURN
       END
C     ***************************************************************
C     C lculo de energ!a potencial

      SUBROUTINE ENERG(UTOT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NPART=2000,NACC=20,NG=30000)
      COMMON /BPOSITN/ RX(NPART),RY(NPART),RA(NPART),ACC(NACC),G(NG),AR,
     +                 CLU(NPART,NPART)
      COMMON /BCONSTR/ PI,ETA,RHO,XLAMBDA,XLAM2,XLAM3,SL2,SQL2,
     +                 SL3,SQL3,XN,TEMP,DISPL,DISPLAng,XHISTG,E1,E2,E3,
     +                 S,SS,SSLL,SL,XC,YC,YCINV,
     +                 XL,YL,Y2,EA1,EA2,AAng, RA1,RA2,RA3
      COMMON /BCONSTG/ DR1,DR2,DR3,DR4,DR5,PHI,TAU
      COMMON /BCOSTRX/ SL4,SQL4,SL5,SQL5,XLAM4,XLAM5,E4,E5
      COMMON /BCOSTXX/ DR6,SL6,SQL6,XLAM6,E6
      COMMON /BCONSTI/ N,NGOFR,LGOFR,NMOVE,NMOVE2,NSUB,NGOFR0,ISEED
      COMMON /BCONSTV/ PRESS,VOL,SDISPL,NSET,LRHO
      LOGICAL LGOFR
      UTOT=0.0D00
      DO 4 I=1,N-1
      DO 4 J=I+1,N
      IF (J.EQ.I) GOTO 4
      X=RX(J)-RX(I)
      Y=RY(J)-RY(I)
C Convenci"n de imagen m!nima
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
      
      !!Energia pozos angulares
       if (RR.lt.RA3*RA3*SS )then
          pppp= X*cos(RA(I)*PI/180)+Y*sin(RA(I)*PI/180)
          pppp=pppp/sqrt(RR)
          ANGLE=ACOS(pppp)
          ANGLE=ANGLE*180/PI
          pppp= X*cos(RA(J)*PI/180)+Y*sin(RA(J)*PI/180)
          pppp=pppp/sqrt(RR)
          ANGLE2=ACOS(pppp)
          ANGLE2=ANGLE2*180/PI
      if (ANGLE>(90-AAng/2).and.ANGLE<(90+AAng/2).and.RR<RA1*RA1*SS)then
        if (ANGLE2>(90-AAng/2).and.ANGLE2<(90+AAng/2)) then
             UTOT=UTOT+EA1
            end if
          end if
       if (ANGLE>(180-AAng/2).or.ANGLE<AAng/2.and.RR>RA2*RA2*SS) then
          if (ANGLE2>(180-AAng/2).or.ANGLE2<AAng/2)then
             UTOT=UTOT+EA2
          end if
          end if
       end if
      
      IF (RR.LT.SSLL) THEN
         UTOT=UTOT+E1+(sqrt(RR)-S)*(E2-E1)/(SL-S)
         ELSE
         IF(RR.LT.SQL2) THEN
            UTOT=UTOT+E2+(sqrt(RR)-SL)*(E3-E2)/(SL2-SL)
            ELSE
            IF (RR.LT.SQL3) THEN
               UTOT=UTOT+E3
               ELSE
               IF (RR.LT.SQL4) THEN
                  UTOT=UTOT+E4
                  ELSE
                  IF (RR.LT.SQL5) THEN
                     UTOT=UTOT+E5
                     ELSE
                        IF (RR.LT.SQL6) UTOT=UTOT+E6
                  END IF
               END IF
            END IF
         END IF
      END IF
    4 CONTINUE
      RETURN
      END


C     ************************************************************
C     Calcula g(r)

      SUBROUTINE GOFR
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NPART=2000,NACC=20,NG=30000)
      COMMON /BPOSITN/ RX(NPART),RY(NPART),RA(NPART),ACC(NACC),G(NG),AR,
     +                 CLU(NPART,NPART)
      COMMON /BCONSTR/ PI,ETA,RHO,XLAMBDA,XLAM2,XLAM3,SL2,SQL2,
     +                 SL3,SQL3,XN,TEMP,DISPL,DISPLAng,XHISTG,E1,E2,E3,
     +                 S,SS,SSLL,SL,XC,YC,YCINV,
     +                 XL,YL,Y2,EA1,EA2,AAng, RA1,RA2,RA3
      COMMON /BCONSTG/ DR1,DR2,DR3,DR4,DR5,PHI,TAU
      COMMON /BCOSTRX/ SL4,SQL4,SL5,SQL5,XLAM4,XLAM5,E4,E5
      COMMON /BCOSTXX/ DR6,SL6,SQL6,XLAM6,E6
      COMMON /BCONSTI/ N,NGOFR,LGOFR,NMOVE,NMOVE2,NSUB,NGOFR0,ISEED
      COMMON /BCONSTV/ PRESS,VOL,SDISPL,NSET,LRHO
      LOGICAL LGOFR
      ACC(3)=ACC(3)+1.0D00
      RMAX=0.5D00*0.5D00
      DO 1 I=1,N-1
      DO 1 J=I+1,N
      X=RX(I)-RX(J)
      Y=RY(I)-RY(J)
C Convenci�n de imagen m�nima
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
      IF (RR.GT.RMAX) GOTO 1
      R=SQRT(RR)
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

    1 CONTINUE
      NGOFR=NGOFR+NGOFR0
      RETURN
      END

C     *****************************************************************
C     Limpia para terminar

      SUBROUTINE FINISH
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NPART=2000,NACC=20,NG=30000)
      COMMON /BPOSITN/ RX(NPART),RY(NPART),RA(NPART),ACC(NACC),G(NG),AR,
     +                 CLU(NPART,NPART)
      COMMON /BCONSTR/ PI,ETA,RHO,XLAMBDA,XLAM2,XLAM3,SL2,SQL2,
     +                 SL3,SQL3,XN,TEMP,DISPL,DISPLAng,XHISTG,E1,E2,E3,
     +                 S,SS,SSLL,SL,XC,YC,YCINV,
     +                 XL,YL,Y2,EA1,EA2,AAng, RA1,RA2,RA3
      COMMON /BCONSTG/ DR1,DR2,DR3,DR4,DR5,PHI,TAU
      COMMON /BCOSTRX/ SL4,SQL4,SL5,SQL5,XLAM4,XLAM5,E4,E5
      COMMON /BCOSTXX/ DR6,SL6,SQL6,XLAM6,E6
      COMMON /BCONSTI/ N,NGOFR,LGOFR,NMOVE,NMOVE2,NSUB,NGOFR0,ISEED
      COMMON /BCONSTV/ PRESS,VOL,SDISPL,NSET,LRHO
      LOGICAL LGOFR
C Escribe configuracion final en p2a.new
      OPEN(UNIT=4,FILE='npt6.new',STATUS='NEW')
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
C      OPEN(UNIT=15,FILE='npt5.pot',STATUS='NEW')
C      WRITE(15,*) 1.0D0, E1
C      WRITE(15,*) XLAMBDA, E1
C      WRITE(15,*) XLAMBDA, E2
C      WRITE(15,*) XLAM2, E2
C      WRITE(15,*) XLAM2, E3
C      WRITE(15,*) XLAM3, E3
C      WRITE(15,*) XLAM3, E4
C      WRITE(15,*) XLAM4, E4
C      WRITE(15,*) XLAM4, E5
C      WRITE(15,*) XLAM5, E5
C      WRITE(15,*) XLAM5, E6
C      WRITE(15,*) XLAM6, E6
C      WRITE(15,*) XLAM6, 0.0D0
      RETURN
      END

C     *****************************************************************
C     Calcula el perfil

      SUBROUTINE RADIAL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NPART=2000,NACC=20,NG=30000)
      COMMON /BPOSITN/ RX(NPART),RY(NPART),RA(NPART),ACC(NACC),G(NG),AR,
     +                 CLU(NPART,NPART)
      COMMON /BCONSTR/ PI,ETA,RHO,XLAMBDA,XLAM2,XLAM3,SL2,SQL2,
     +                 SL3,SQL3,XN,TEMP,DISPL,DISPLAng,XHISTG,E1,E2,E3,
     +                 S,SS,SSLL,SL,XC,YC,YCINV,
     +                 XL,YL,Y2,EA1,EA2,AAng, RA1,RA2,RA3
      COMMON /BCONSTG/ DR1,DR2,DR3,DR4,DR5,PHI,TAU
      COMMON /BCOSTRX/ SL4,SQL4,SL5,SQL5,XLAM4,XLAM5,E4,E5
      COMMON /BCOSTXX/ DR6,SL6,SQL6,XLAM6,E6
      COMMON /BCONSTI/ N,NGOFR,LGOFR,NMOVE,NMOVE2,NSUB,NGOFR0,ISEED
      COMMON /BCONSTV/ PRESS,VOL,SDISPL,NSET,LRHO
      LOGICAL LGOFR
      NHISTG=nint(XHISTG)
      XNAV=ACC(3)
      IF(XNAV.EQ.0.0D00)RETURN
C     DELR=(SL-S)*XL/XHISTG
      VOLME=XL*YL
      CONST=XN*XN*PI*XNAV/(VOLME*2.0D00)
C      NEND=INT((0.5D00*XC-S)*XHISTG/(SL-S))
      NEND=INT(5.0D00*XHISTG+((0.5D00*XC)-SL5)/DR6)

      OPEN(UNIT=7,FILE='gdrnpt5.dat',STATUS='NEW')

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
C      R=(SL-S)*(2.0D00*X-1.0D00)*XL/(XHISTG*2.0D00)
C      RB=(SL-S)*(X-1.0D00)*XL/XHISTG
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
