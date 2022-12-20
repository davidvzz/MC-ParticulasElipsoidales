SUBROUTINE START

    !Lee variables de sistema y estado (rho,temp) de pozos.in
    !al inicio y pozos.old si es continuaci"n
    !Genera la configuraci"n inicial de fcc.
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
    !      LOGICAL LRHO
    ! ACUMULADORES A CERO
      DO 3 I=1,NACC
      3 ACC(I)=0.0D00
      DO 4 I=1,NG
      4 G(I)=0.0D00
        
      CLU(:,:)=0.0D0

    ! LEE PARAMETROS DE ENTRADA: N,TEMP E.T.C.
      OPEN (UNIT=1,FILE='npt6.in',STATUS='unknown')
      READ(1,*) N         !1
      READ(1,*) RHO       !2
      READ(1,*) TEMP      !3
      READ(1,*) PHI       !4
      READ(1,*) XLAMBDA   !5
      READ(1,*) XLAM2     !6
      READ(1,*) XLAM3     !7
      READ(1,*) XLAM4     !8
      READ(1,*) XLAM5     !9
      READ(1,*) XLAM6     !10
      READ(1,*) E1        !11
      READ(1,*) E2        !12
      READ(1,*) E3        !13
      READ(1,*) E4        !14
      READ(1,*) E5        !15
      READ(1,*) E6        !16
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
    ! LRHO DEFINE SI SE PARTE DE UNA DENSIDAD NUEVA
     !     READ(1,*) LRHO
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
      WRITE(6,*) 'E1=',E1
      WRITE(6,*) 'E2=',E2
      WRITE(6,*) 'E3=',E3
      WRITE(6,*) 'E4=',E4
      WRITE(6,*) 'E5=',E5
      WRITE(6,*) 'E6=',E6
      WRITE(6,*) 'PHI=',PHI
      WRITE(6,*)(XA(I),YA(I),I=1,NAC)

      IF (NRUN == 0.0D00) WRITE(6,102)
      IF (NRUN.NE.0.0D00) WRITE(6,103)
      IF (IAV == 0.0D00)  WRITE(6,104)
      IF (IAV.NE.0.0D00)  WRITE(6,105)
      IF (LGOFR)     WRITE(6,106) XHISTG

    ! CONVERSION A UNIDADES DE PROGRAMA (Lx = 1)
      PRESS=1.1547005384D0*PHI
      XN=DFLOAT(N)

    ! SE CALCULA EL TAMA~NO DE LA MUESTRA PARA BARRAS DE ERROR EN SIGMA
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

      IF (NRUN == 0) THEN
       ! CALCULA SIGMA EN UNIDADES DE CAJA
            ETA=RHO*PI/4.0D00
            XL=(XN/(RHO*YC))**(1.0D00/2.0D00)
            YL=XL*YC
            S=1.0D00/XL
       !  CONSTRUYE LATIZ INICIAL CON SIGMA DADA
            ISEED=-123456789
      
            DO I=1,N
                  9 RX(I)=RAN2(ISEED)
                  RY(I)=YC*RAN2(ISEED)
                  RA(I)=180*RAN2(ISEED)
                  DO J=1,I
                        IF (J .NE. I)then
                              X=RX(J)-RX(I)
                              Y=RY(J)-RY(I)
                  
                              ! convencion de imagen minima
                              IF (X > 0.5D00) THEN
                                    X=X-1.0D00
                              ELSE IF (X < -0.5D00) THEN
                                    X=X+1.0D00
                              END IF
                              IF (Y > Y2) THEN
                                    Y=Y-YC
                              ELSE IF (Y < -Y2) THEN
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
                        
                              if ( FI < 0 .or. (FI > 0 .and. F1 >= 0 .and. F2 >= 0) .or. RR < S*S ) GOTO 9
                        
                        end if
                  end do
            end do
      ELSE
            OPEN(UNIT=3,FILE='npt6.old',STATUS='OLD')
            IF (LRHO == 1) THEN
            ! CALCULA SIGMA CON DENSIDAD DADA Y CONFIGURACION VIEJA
                  ETA=RHO*PI/4.0D0
                  XL=(XN/(RHO*YC))**(1.0D00/2.0D00)
                  YL=XL*YC
                  S=1.0D00/XL
                  READ(3,*) RX(1)
            ELSE
            !   LEE VALOR DE SIGMA
                  READ(3,*) S
                  VOL=YC/(S*S)
                  RHO=XN/VOL
                  ETA=RHO*PI/4.0D00
                  XL=1/S
                  YL=XL*YC
            END IF
       !     LEE CONFIGURACION VIEJA
            DO 11 I=1,N
       !       WRITE(3,*) RX(I),RY(I), RA(I)
                  READ(3,*) RX(I),RY(I),Ancho,S, RA(I)
        11  CONTINUE
            ! lee acumuladores viejos
            IF (IAV == 0) GOTO 2

            DO 22 I=1,NACC
                  READ(3,*)ACC(I)
            22  CONTINUE
            DO 33 I=1,NG
                  READ(3,*)G(I)
            33  CONTINUE
            ! escribe acumulador y configuraci"n de inicio o re-inicio
            WRITE(6,111) ACC(1)
            WRITE(6,112)
2     END IF

    !     CONTINUA CONVIRTIENDO A UNIDADES DE PROGRAMA
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

    ! ESCRIBE LAS UNIDADES CONVERTIDAS
      WRITE(6,107) RHO
      WRITE(6,108) XC,YC
      WRITE(6,109) XL,YL
      WRITE(6,110) S
      DO 44 I=1,20
            WRITE(6,115) RX(I),RY(I), RA(I)
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

END SUBROUTINE