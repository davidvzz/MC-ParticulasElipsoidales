! Calculo de energ!a potencial

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
! Convenci"n de imagen m!nima
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
