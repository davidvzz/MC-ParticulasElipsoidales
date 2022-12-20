SUBROUTINE MCARLO
! CALCULA LA ENERGIA POTENCIAL DE LA CONFIGURACION INICIAL
   IMPLICIT DOUBLE PRECISION( A-H,O-Z )
   PARAMETER ( NPART=2000,NACC=20,NG=30000 )
   COMMON /BPOSITN/ RX( NPART ),RY( NPART ),RA( NPART ),ACC( NACC ),G( NG ),AR,
   +                 CLU( NPART,NPART )
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
   DOUBLE PRECISION XNEW1( N ), YNEW1( N ), ANEW1( N )

! inicializa secuencia al azar
   ISEED=-123456789
   CALL ENERG( UTOT )
   UPERP=UTOT/XN
!  NTEST=0
   USUBAV=0.0D00
   NCOUNT=0

   NERCONT=0
   NACCPT=0
   RTEST=0.0D00
   NSUB0=NSUB
   NGOFR0=NGOFR
   WRITE( 6,101 ) UPERP

! comienza la secuencia de Monte Carlo

   OPEN( UNIT=2,FILE='unpt6.dat',STATUS='NEW' )
   !OPEN( UNIT=93,FILE='angulos.dat',STATUS='NEW' )

! escoge part!cula al azar
   1 I=INT( RAN2( ISEED )*N )+1
   NCOUNT=NCOUNT+1

!  NERCONT ES EL CONTADOR QUE AL LLEGAR A NSET INDICA MUESTREAR SIGMA
   NERCONT=NERCONT+1
   RTEST=RTEST+RAN2( ISEED )
! desplaza la part!cula
   XNEW=RX( I )+DISPL*( RAN2( ISEED )-0.5D00 )
   YNEW=RY( I )+DISPL*( RAN2( ISEED )-0.5D00 )
   !!!angular
   ANEW=RA( I )+DISPLAng*( RAN2( ISEED )-0.5D0 )
! condici"n peri"dica de frontera
   IF( XNEW > 1.0D00 ) THEN
      XNEW=XNEW-1.0D00
   ELSE IF ( XNEW < 0.0D00 ) THEN
      XNEW=XNEW+1.0D00
   END IF

   IF ( YNEW > YC ) THEN
      YNEW=YNEW-YC
   ELSE IF ( YNEW < 0.0D00 ) THEN
      YNEW=YNEW+YC
   END IF
   
   if ( ANEW > 180 ) then
      ANEW=ANEW-180
   end if
   if( ANEW < 0 ) then
      ANEW=ANEW+180
   end if

! elimina part!culas traslapadas
! y calcula la nueva energ!a ( provisional )
   !****************
   UNEW=0.0D00
   DO 2 J=1,N

      IF ( J == I ) GOTO 2

      X=RX( J )-XNEW
      Y=RY( J )-YNEW
      !************

    ! convenci"n de imagen m!nima
      IF ( X > 0.5D00 ) THEN
         X=X-1.0D00
      ELSE IF ( X < -0.5D00 ) THEN
         X=X+1.0D00
      END IF

      IF ( Y > Y2 ) THEN
         Y=Y-YC
      ELSE IF ( Y < -Y2 ) THEN
         Y=Y+YC
      END IF
      RR=X*X+Y*Y

      !!!!CONDICIONES DE TRASLAPE
      GGG = 2.00 +( AR-( 1/AR ) )**2*( sin( ( RA( J )-ANEW )*PI/180 ) )**2
      F1=1.00 + GGG
      F1=F1-( ( X*cos( ANEW*PI/180 )+ Y*sin( ANEW*PI/180 ) )**2 )/( AR*S/2 )**2
      F1=F1-( ( Y*cos( ANEW*PI/180 )- X*sin( ANEW*PI/180 ) )**2 )/( SS/4 )
      F2=1.00 +GGG
      F2=F2-( ( X*cos( RA( J )*PI/180 )+ Y*sin( RA( J )*PI/180 ) )**2 )/( AR*S/2 )**2
      F2=F2-( ( Y*cos( RA( J )*PI/180 )- X*sin( RA( J )*PI/180 ) )**2 )/( SS/4 )
      FI=4*( F1**2-3*F2 )*( F2**2-3*F1 )-( 9-F1*F2 )**2

      !FI=0 es tangencia
      !no puntos en com�n es FI positiva y al menos uno de f1, f2 negativos
         
      if ( FI < 0 .or. ( FI > 0 .and. F1 > =0 .and. F2 > =0 ) .or. RR < SS ) GOTO 3

         
         !!Energia pozos angulares
      if ( RR < RA3*RA3*SS  )then
         pppp= X*cos( ANEW*PI/180 )+Y*sin( ANEW*PI/180 )
         pppp=pppp/sqrt( RR )
         ANGLE=ACOS( pppp )
         ANGLE=ANGLE*180/PI
         pppp= X*cos( RA( J )*PI/180 )+Y*sin( RA( J )*PI/180 )
         pppp=pppp/sqrt( RR )
         ANGLE2=ACOS( pppp )
         ANGLE2=ANGLE2*180/PI

         if ( ANGLE > ( 90-AAng/2 ) .and. ANGLE < ( 90+AAng/2 ) .and. RR < RA1*RA1*SS )then
            if ( ANGLE2 > ( 90-AAng/2 ) .and. ANGLE2 < ( 90+AAng/2 )  )then
               UNEW=UNEW+EA1
            end if
         end if

         if ( ANGLE > ( 180-AAng/2 ) .or. ANGLE < AAng/2 .and. RR > RA2*RA2*SS )then
            if ( ANGLE2 > ( 180-AAng/2 ) .or. ANGLE2 < AAng/2 ) then
               UNEW=UNEW+EA2
            end if
         end if

      end if
         

      IF ( RR < SSLL ) THEN
         UNEW=UNEW+ E1+( sqrt( RR )-S )*( E2-E1 )/( SL-S )

      ELSE
         IF( RR < SQL2 ) THEN
            UNEW=UNEW+ E2+( sqrt( RR )-SL )*( E3-E2 )/( SL2-SL )

         ELSE
            IF ( RR < SQL3 ) THEN
               UNEW=UNEW+E3

            ELSE
               IF ( RR < SQL4 ) THEN
                  UNEW=UNEW+E4

               ELSE
                  IF ( RR < SQL5 ) THEN
                     UNEW=UNEW+E5

                  ELSE
                     IF ( RR < SQL6 ) UNEW=UNEW+E6

                  END IF
               END IF
            END IF
         END IF
      END IF

   2 CONTINUE
   UOLD=0.0D00

    !*********
   DO 4 J=1,N
      
      IF ( J == I ) GOTO 4
      X=RX( J )-RX( I )
      Y=RY( J )-RY( I )
    !***************

    ! convenci"n de imagen m!nima
      IF ( X > 0.5D00 ) THEN
         X=X-1.0D00
      ELSE IF ( X < -0.5D00 ) THEN
         X=X+1.0D00
      END IF

      IF ( Y > Y2 ) THEN
         Y=Y-YC
      ELSE IF ( Y < -Y2 ) THEN
         Y=Y+YC
      END IF

      RR=X*X+Y*Y
      !!!!CONDICIONES DE TRASLAPE
      GGG = 2.00 +( AR-( 1/AR ) )**2*( sin( ( RA( J )-RA( I ) )*PI/180 ) )**2
      F1=1.00 + GGG
      F1=F1-( ( X*cos( RA( I )*PI/180 )+ Y*sin( RA( I )*PI/180 ) )**2 )/( AR*S/2 )**2
      F1=F1-( ( Y*cos( RA( I )*PI/180 )- X*sin( RA( I )*PI/180 ) )**2 )/( SS/4 )
      F2=1.00 +GGG
      F2=F2-( ( X*cos( RA( J )*PI/180 )+ Y*sin( RA( J )*PI/180 ) )**2 )/( AR*S/2 )**2
      F2=F2-( ( Y*cos( RA( J )*PI/180 )- X*sin( RA( J )*PI/180 ) )**2 )/( SS/4 )
      FI=4*( F1**2-3*F2 )*( F2**2-3*F1 )-( 9-F1*F2 )**2
      if ( FI < 0 .or. ( FI > 0 .and. F1 > =0 .and. F2 > =0 ) .or. RR < SS ) GOTO 42

      !!Energia pozos angulares
      if ( RR < RA3*RA3*SS ) then
         pppp= X*cos( RA( I )*PI/180 )+Y*sin( RA( I )*PI/180 )
         pppp=pppp/sqrt( RR )
         ANGLE=ACOS( pppp )
         ANGLE=ANGLE*180/PI
         pppp= X*cos( RA( I )*PI/180 )+Y*sin( RA( I )*PI/180 )
         pppp=pppp/sqrt( RR )
         ANGLE2=ACOS( pppp )
         ANGLE2=ANGLE2*180/PI

         if ( ANGLE > ( 90-AAng/2 ) .and. ANGLE < ( 90+AAng/2 ) .and. RR < RA1*RA1*SS )then
            if ( ANGLE2 > ( 90-AAng/2 ) .and. ANGLE2 < ( 90+AAng/2 )  )then
               UOLD=UOLD+EA1
            end if
         end if

         if ( ANGLE > ( 180-AAng/2 ) .or. ANGLE < AAng/2 .and. RR > RA2*RA2*SS  )then
            if ( ANGLE2 > ( 180-AAng/2 ) .or. ANGLE2 < AAng/2  )then
               UOLD=UOLD+EA2
            end if
         end if

      end if
         
      IF ( RR < SSLL ) THEN
         UOLD=UOLD+E1+( sqrt( RR )-S )*( E2-E1 )/( SL-S )

      ELSE
         IF( RR < SQL2 ) THEN
            UOLD=UOLD+E2+( sqrt( RR )-SL )*( E3-E2 )/( SL2-SL )

         ELSE
            IF ( RR < SQL3 ) THEN
               UOLD=UOLD+E3

            ELSE
               IF ( RR < SQL4 ) THEN
                  UOLD=UOLD+E4

               ELSE
                  IF ( RR < SQL5 ) THEN
                     UOLD=UOLD+E5

                  ELSE
                        IF ( RR < SQL6 ) UOLD=UOLD+E6
                  END IF
               END IF
            END IF
         END IF
      END IF
   4 CONTINUE
   DENERG=UNEW-UOLD
   
   IF ( DENERG <= 0.0D00 ) GOTO 5
   
! compara el factor de BOLTZMANN con n#mero al azar
   RND=RAN2( ISEED )
   IF ( RND > EXP( -DENERG/TEMP ) ) GOTO 3

! actualiza posici"n de la part!cula I
   5 RX( I )=XNEW
     RY( I )=YNEW
     RA( I )=ANEW
     call ENERG( UTOT )
     NACCPT=NACCPT+1.0D00

! acumula promedios
   3 ACC( 1 )=ACC( 1 )+1.0D00
     ACC( 2 )=ACC( 2 )+UTOT
     USUBAV=USUBAV+UTOT
     XNTEST=DFLOAT( NCOUNT )

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!        PRUEBA PARA OBTENER BARRA DE ERROR EN DENSIDAD
   IF ( NERCONT == NSET ) THEN
      WRITE( 10,* )RHO, S, VOL
      NERCONT=0
   END IF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   IF ( NCOUNT == NGOFR .and. LGOFR ) CALL GOFR
   IF ( NCOUNT >= NMOVE ) GOTO 6

   IF ( NCOUNT < NSUB ) GOTO 1
! escribe resultados

   UAV=USUBAV/( XN*NCOUNT )

!      WRITE( 6,102 ) NCOUNT,DFLOAT( NACCPT )/DFLOAT( NCOUNT ),RTEST/XNTEST,
!     + UAV,VACCPT
   WRITE( 6,102 ) NCOUNT,DFLOAT( NACCPT )/DFLOAT( NCOUNT ),TAU,
   + UAV,VACCPT,RTEST/XNTEST
   WRITE( 8,* )TAU
   NSUB=NSUB+NSUB0
   XXMOV = DFLOAT( NMOVE )
   XXX = XNTEST/XXMOV
! verifica si es necesario cambiar desplazamiento maximo
   !IF ( DFLOAT( NACCPT )/DFLOAT( NCOUNT ) > 0.45 ) DISPL=DISPL*1.05
   !IF ( DFLOAT( NACCPT )/DFLOAT( NCOUNT ) < 0.35 ) DISPL=DISPL*0.95
! verifica si es necesario cambiar CAMBIO DE VOLUMEN  maximo
   !IF ( VACCPT > 0.45 .and. LVOL ) SDISPL=SDISPL*1.05
   !IF ( VACCPT < 0.35 .and. LVOL ) SDISPL=SDISPL*0.95
   WRITE( 2,* )XXX,UAV
   !WRITE( 23,* ) XN, NCOUNT
   GOTO 1

! fin de la corrida por part�cula
! marca clusters

   !! Revisa vecinos para clusters
   6  DO I=1,N
         DO J=1,N
            IF ( J.NE.I ) then
               X=RX( I )-RX( J )
               Y=RY( I )-RY( J )

               ! convenci"n de imagen m!nima
               IF ( X > 0.5D00 ) THEN
                  X=X-1.0D00
               ELSE IF ( X < -0.5D00 ) THEN
                  X=X+1.0D00
               END IF

               IF ( Y > Y2 ) THEN
                  Y=Y-YC
               ELSE IF ( Y < -Y2 ) THEN
                  Y=Y+YC
               END IF

               RR=X*X+Y*Y

               if ( RR < ( RA3+3 )*( RA3+3 )*SS ) then
                  CLU( I,J )=1.0D0
               else
                  CLU( I,J )=0.0D0
               end if
            end if
         END DO
      END DO

      I=1
      DO while ( I.NE.0 )
         CALL CLUSTERS( I )
      end do

      OPEN( UNIT=15,FILE='antes.dat',STATUS='unknown' )
      DO I=1,N
         WRITE( 15,* )RX( I ),RY( I ),S*AR,S,RA( I )
      END DO

! inicia corrida por clusters
   7  I=INT( RAN2( ISEED )*N )+1
      NCOUNT=NCOUNT+1

      !NERCONT ES EL CONTADOR QUE AL LLEGAR A NSET INDICA MUESTREAR SIGMA
      NERCONT=NERCONT+1
      RTEST=RTEST+RAN2( ISEED )
      ! desplaza la part!cula
      XNEW1( 1:N )=RX( 1:N )
      YNEW1( 1:N )=RY( 1:N )
      ANEW1( 1:N )=RA( 1:N )
      ! centro de masa
      K=1
      CX= XNEW1( I )
      CY=YNEW1( I )
      DO J=1,N
         IF ( J.NE.I .and. CLU( I,J ).NE.0.0D0 ) THEN
            CX= CX+XNEW1( J )
            CY=CY+YNEW1( J )
            K=K+1
         end if
      END DO

      CX=CX/DFLOAT( K )
      CY=CY/DFLOAT( K )
      K=0

      IF ( ABS( XNEW1( I )-CX ) < 0.4 .and. ABS( YNEW1( I )-CY ) < 0.4 ) then
         K=K+1
      end if

      DO J=1,N
         IF ( J.NE.I .and. CLU( I,J ).NE.0.0D0 ) THEN
            IF ( ABS( XNEW1( J )-CX ) < 0.4 .and. ABS( YNEW1( J )-CY ) < 0.4 ) then
               K=K+1
            end if
         end if
      END DO

      IF ( K==0 ) GOTO 8

      ! giro en angulo
      DA=DISPLAng*( RAN2( ISEED )-0.5D0 )
      ANEW1( I )=ANEW1( I )+DA
      X=XNEW1( I )-CX
      Y=YNEW1( I )-CY
      XNEW1( I )=X*COS( DA*PI/180 )-Y*SIN( DA*PI/180 )+CX
      YNEW1( I )=X*SIN( DA*PI/180 )+Y*COS( DA*PI/180 )+CY

      do J=1,N
         if ( J.NE.I .and. CLU( I,J ).NE.0.0D0 ) then
            ANEW1( J )=ANEW1( J )+DA
            X=XNEW1( J )-CX
            Y=YNEW1( J )-CY
            XNEW1( J )=X*COS( DA*PI/180 )-Y*SIN( DA*PI/180 )+CX
            YNEW1( J )=X*SIN( DA*PI/180 )+Y*COS( DA*PI/180 )+CY
         end if
      end do

! desplazamiento  X, Y
   8  DX=DISPLClu*( RAN2( ISEED )-0.5D00 )
      DY=DISPLClu*( RAN2( ISEED )-0.5D00 )
      XNEW1( I )=XNEW1( I )+DX
      YNEW1( I )=YNEW1( I )+DY
      IF( XNEW1( I ) > 1.0D00 ) THEN
         XNEW1( I )=XNEW1( I )-1.0D00
      ELSE
         IF ( XNEW1( I ) < 0.0D00 ) THEN
            XNEW1( I )=XNEW1( I )+1.0D00
         END IF
      End if

      IF ( YNEW1( I ) > YC ) THEN
        YNEW1( I )=YNEW1( I )-YC
      ELSE
         IF ( YNEW1( I ) < 0.0D00 ) THEN
            YNEW1( I )=YNEW1( I )+YC
         END IF
      End if
      
      if ( ANEW1( I ) > 180 ) then
         ANEW1( I )=ANEW1( I )-180
      end if

      if( ANEW1( I ) < 0 ) then
         ANEW1( I )=ANEW1( I )+180
      end if

      do J            ! condici"n peri"dica de frontera
         IF( XNEW1( J ) > 1.0D00 ) THEN
            XNEW1( J )=XNEW1( J )-1.0D00
         ENd if
         IF ( XNEW1( J ) < 0.0D00 ) THEN
            XNEW1( J )=XNEW1( J )+1.0D00
         END IF
         IF ( YNEW1( J ) > YC ) THEN
            YNEW1( J )=YNEW1( J )-YC
         END IF
            IF ( YNEW1( J ) < 0.0D00 ) THEN
         YNEW1( J )=YNEW1( J )+YC
         END IF
         if ( ANEW1( J ) > 180 ) then
            ANEW1( J )=ANEW1( J )-180
         end if
         if( ANEW1( J ) < 0 ) then
            ANEW1( J )=ANEW1( J )+180
         end if=1,N

         if ( J.NE.I .and. CLU( I,J ).NE.0.0D0 ) then
            XNEW1( J )=XNEW1( J )+DX
            YNEW1( J )=YNEW1( J )+DY
            ! condici"n peri"dica de frontera
            IF( XNEW1( J ) > 1.0D00 ) THEN
               XNEW1( J )=XNEW1( J )-1.0D00
            ENd if
            IF ( XNEW1( J ) < 0.0D00 ) THEN
               XNEW1( J )=XNEW1( J )+1.0D00
            END IF
            IF ( YNEW1( J ) > YC ) THEN
               YNEW1( J )=YNEW1( J )-YC
            END IF
            IF ( YNEW1( J ) < 0.0D00 ) THEN
               YNEW1( J )=YNEW1( J )+YC
            END IF
            if ( ANEW1( J ) > 180 ) then
               ANEW1( J )=ANEW1( J )-180
            end if
            if( ANEW1( J ) < 0 ) then
               ANEW1( J )=ANEW1( J )+180
            end if
         end if
      end do

      ! revisa traslapes
      DO J=1,N
         IF ( J.NE.I )then
            X=XNEW1( J )-XNEW1( I )
            Y=YNEW1( J )-YNEW1( I )

            !      convenci"n de imagen m!nima
            IF ( X > 0.5D00 ) THEN
               X=X-1.0D00
            ELSE IF ( X < -0.5D00 ) THEN
               X=X+1.0D00
            END IF

            IF ( Y > Y2 ) THEN
               Y=Y-YC
            ELSE IF ( Y < -Y2 ) THEN
               Y=Y+YC
            END IF

            RR=X*X+Y*Y
            !!!!CONDICIONES DE TRASLAPE
            GGG = 2.00 +( AR-( 1/AR ) )**2*( sin( ( ANEW1( J )-ANEW1( I ) )*PI/180 ) )**2
            F1=-( ( X*cos( ANEW1( I )*PI/180 )+ Y*sin( ANEW1( I )*PI/180 ) )**2 )
            F1=F1/( AR*S/2 )**2
            F1=F1+1.00 + GGG

            F1=F1-( ( Y*cos( ANEW1( I )*PI/180 )- X*sin( ANEW1( I )*PI/180 ) )**2 )/( SS/4 )

            F2=-( ( X*cos( ANEW1( J )*PI/180 )+ Y*sin( ANEW1( J )*PI/180 ) )**2 )
            F2=F2/( AR*S/2 )**2
            F2=F2+1.00 +GGG
            F2=F2-( ( Y*cos( ANEW1( J )*PI/180 )- X*sin( ANEW1( J )*PI/180 ) )**2 )/( SS/4 )
            FI=4*( F1**2-3*F2 )*( F2**2-3*F1 )-( 9-F1*F2 )**2

            if ( FI < 0 .or. ( FI > 0 .and. F1 > =0 .and. F2 > =0 ) .or. RR < SS ) GOTO 9
         end if
      end do
      
      do J=1,N
         IF( I.NE.J .and. CLU( I,J ).NE.0.0D0 ) then
            do M=1,N
               if ( J.NE.M ) then
                  X=XNEW1( M )-XNEW1( J )
                  Y=YNEW1( M )-YNEW1( J )

                  ! convenci"n de imagen m!nima
                  IF ( X > 0.5D00 ) THEN
                     X=X-1.0D00
                  end if
                  IF ( X < -0.5D00 ) THEN
                     X=X+1.0D00
                  END IF
                  IF ( Y > Y2 ) THEN
                     Y=Y-YC
                  end if
                  IF ( Y < -Y2 ) THEN
                     Y=Y+YC
                  END IF
                  RR=X*X+Y*Y
                  !!!!CONDICIONES DE TRASLAPE
                  GGG = 2.00 +( AR-( 1/AR ) )**2*( sin( ( ANEW1( M )-ANEW1( J ) )*PI/180 ) )**2
                  F1=-( ( X*cos( ANEW1( J )*PI/180 )+ Y*sin( ANEW1( J )*PI/180 ) )**2 )
                  F1=F1/( AR*S/2 )**2
                  F1=F1+1.00 + GGG

                  F1=F1-( ( Y*cos( ANEW1( J )*PI/180 )- X*sin( ANEW1( J )*PI/180 ) )**2 )/( SS/4 )

                  F2=-( ( X*cos( ANEW1( M )*PI/180 )+ Y*sin( ANEW1( M )*PI/180 ) )**2 )
                  F2=F2/( AR*S/2 )**2
                  F2=F2+1.00 +GGG
                  F2=F2-( ( Y*cos( ANEW1( M )*PI/180 )- X*sin( ANEW1( M )*PI/180 ) )**2 )/( SS/4 )
                  FI=4*( F1**2-3*F2 )*( F2**2-3*F1 )-( 9-F1*F2 )**2

                  if ( FI < 0 .or. ( FI > 0 .and. F1 > =0 .and. F2 > =0 ) .or. RR < SS ) GOTO 9
               end if
            end do
         end if
      end do


! actualiza posici"n de la part!cula I
   11 RX( 1:N )=XNEW1( 1:N )
      RY( 1:N )=YNEW1( 1:N )
      RA( 1:N )=ANEW1( 1:N )


      call ENERG( UTOT )
      !! Revisa vecinos para clusters
      DO I=1,N
         DO J=1,N
            IF ( J.NE.I ) then
               X=RX( I )-RX( J )
               Y=RY( I )-RY( J )

               ! convenci"n de imagen m!nima
               IF ( X > 0.5D00 ) THEN
                  X=X-1.0D00
               ELSE IF ( X < -0.5D00 ) THEN
                  X=X+1.0D00
               END IF

               IF ( Y > Y2 ) THEN
                  Y=Y-YC
               ELSE IF ( Y < -Y2 ) THEN
                  Y=Y+YC
               END IF

               RR=X*X+Y*Y
               if ( RR < ( RA3+3 )*( RA3+3 )*SS ) then
                  CLU( I,J )=1.0D0
               else
                  CLU( I,J )=0.0D0
               end if
            
            end if
         END DO
      END DO
      

      !! ubica clusters
      I=1
      DO while ( I.NE.0 )
         CALL CLUSTERS( I )
      end do
      NACCPT=NACCPT+1.0D00
   

   ! acumula promedios
    9 ACC( 1 )=ACC( 1 )+1.0D00
      ACC( 2 )=ACC( 2 )+UTOT
      USUBAV=USUBAV+UTOT
      XNTEST=DFLOAT( NCOUNT )

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!        PRUEBA PARA OBTENER BARRA DE ERROR EN DENSIDAD
           IF ( NERCONT == NSET ) THEN
             WRITE( 10,* )RHO, S, VOL
             NERCONT=0
           END IF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF ( NCOUNT == NGOFR .and. LGOFR ) CALL GOFR
      IF ( NCOUNT >= NMOVE+NMOVE2 ) GOTO 12

      IF ( NCOUNT < NSUB ) GOTO 7
! escribe resultados

      UAV=USUBAV/( XN*NCOUNT )

!      WRITE( 6,102 ) NCOUNT,DFLOAT( NACCPT )/DFLOAT( NCOUNT ),RTEST/XNTEST,
!     + UAV,VACCPT
      WRITE( 6,102 ) NCOUNT,DFLOAT( NACCPT )/DFLOAT( NCOUNT ),TAU,
     + UAV,VACCPT,RTEST/XNTEST
      WRITE( 8,* )TAU
      NSUB=NSUB+NSUB0
      XXMOV = DFLOAT( NMOVE )
      XXX = XNTEST/XXMOV
! verifica si es necesario cambiar desplazamiento maximo
      !IF ( DFLOAT( NACCPT )/DFLOAT( NCOUNT ) > 0.45 ) DISPL=DISPL*1.05
      !IF ( DFLOAT( NACCPT )/DFLOAT( NCOUNT ) < 0.35 ) DISPL=DISPL*0.95
! verifica si es necesario cambiar CAMBIO DE VOLUMEN  maximo
      !IF ( VACCPT > 0.45 .and. LVOL ) SDISPL=SDISPL*1.05
      !IF ( VACCPT < 0.35 .and. LVOL ) SDISPL=SDISPL*0.95
      WRITE( 2,* )XXX,UAV
      !WRITE( 23,* ) XN, NCOUNT
      GOTO 7


! termina corrida por clusters
   12  UAV=USUBAV/( XN*NCOUNT )
       UAVRUN=ACC( 2 )/( XN*ACC( 1 ) )
       OPEN( UNIT=17,FILE='despues.dat',STATUS='unknown' )
      DO I=1,N
      WRITE( 17,* )RX( I ),RY( I ),S*AR,S,RA( I )
      END DO

!      WRITE( 6,102 ) NCOUNT,DFLOAT( NACCPT )/DFLOAT( NCOUNT ),RTEST/XNTEST,
!     + UAV,VACCPT
      WRITE( 6,102 ) NCOUNT,DFLOAT( NACCPT )/DFLOAT( NCOUNT ),TAU,
     + UAV,VACCPT,RTEST/XNTEST
      WRITE( 8,* )TAU
      WRITE( 6,107 )S
      WRITE( 6,108 )DISPL/S
      WRITE( 6,109 )SDISPL
       WRITE( 6,103 ) ACC( 1 ),ACC( 3 )
       WRITE( 6,104 ) UAVRUN
       RETURN
! part!culas traslapadas
 42   WRITE( 6,105 ) I,J,RR
      WRITE( 6,106 ) SQRT( RR ),S
      WRITE( 6,* ) RX( I ),RY( I )
      WRITE( 6,* ) RX( J ),RY( J )
  101 FORMAT( /,1X,'inicio de corrrida',/,1X,'energ!a inicial',G16.8,/ )
!  102 FORMAT( 1X,'NCOUNT=',I10,' NACCPT RATIO=',F8.4,' RTEST=',F8.4
!     + ,' UAV=',F10.4,' VACCPT RATIO=',F8.4,/ )
  102 FORMAT( 1X,'NCOUNT=',I10,' NACCPT RATIO=',F8.4,' TAU=',F12.8
     + ,' UAV=',F10.4,' VACCPT RATIO=',F8.4,' random  RATIO=',F8.4,/ )
  103 FORMAT( 1X,'Num total de configuraciones desde el inicio',F12.0,/
     +,' N#mero de llamadas a g( r )',F12.0,/ )
  104 FORMAT( 1X,'Energ!a potencial promedio desde el inicio',F10.4,/ )
  105 FORMAT( 1X,'*****part!culas traslapadas*****',' I=',I4,' J=',I4,
     + ' R**2=',G16.8,/ )
  106 FORMAT( 1X,'R=',G16.8,' S=',G16.8,/ )
  107 FORMAT( 1X,' S=',G16.8,/ )
 108  FORMAT( 1X,' DISPL=',G16.8,/ )
 109  FORMAT( 1X,' SDISPL=',G16.8,/ )
      RETURN
      END