        PROGRAM ANGGRAF
            IMPLICIT NONE
            DOUBLE PRECISION:: a, b, ANGLE, ANGLE2, ANGLE3, dist, U,k
            DOUBLE PRECISION :: A1,A2,A3,A4
            DOUBLE PRECISION :: f, ANGINC,ANGLE2_0
            DOUBLE PRECISION, PARAMETER :: pi=4.D0*DATAN(1.D0)
            INTEGER:: T,i,j
            CHARACTER   ::  rstr*2
            INTEGER :: ri(2,2)
            NAMELIST /input/ a,k,b,ANGLE, ANGLE2, ANGLE3, A1, A2, A3, A4
            
            OPEN(UNIT=1,FILE='SS-5-30.dat',STATUS='unknown')     
            OPEN(UNIT=3,FILE='TT-5-30.dat',STATUS='unknown')   
            
            OPEN(UNIT=106,FILE='r6k5.dat',STATUS='unknown') 
            OPEN(UNIT=108,FILE='r8k5.dat',STATUS='unknown') !para SS

            OPEN(UNIT=110,FILE='r10k5.dat',STATUS='unknown') !para TT
            OPEN(UNIT=112,FILE='r12k5.dat',STATUS='unknown')

            ri(1,:) = (/ 6, 8 /)    !configuración 1, valores de distancia 
            ri(2,:) = (/ 10, 12 /)    !configuración 2, valores de distancia
            
            do i=1, 2
                do j=1,2
                    T=ri(i,j)+100
                    ANGINC=0
                    f=0
                    if (ri(i,j)<10) then
                        write(rstr, '(I1)') ri(i,j)
                    else
                        write(rstr, '(I2)') ri(i,j)
                    end if

                    IF (i==1) THEN 
                    !PARA SS CONF R=3B
                        if (j==1) THEN 
                            READ(UNIT=1,nml=input)
                            ANGLE2_0=ANGLE2
                        end if
                        ANGLE2=ANGLE2_0
                        dist=2*b
                        write(T,*) 'k= ',k,'a=',a,'b=',b,'r/b=',ri(i,j)
                        WRITE(T,*) 'x   ','SS-'//rstr//'b'
                    ELSE IF(i.NE.1)THEN
                    !Para TT
                        if (j==1) THEN 
                            READ(UNIT=3,nml=input)
                            ANGLE2_0=ANGLE2
                        end if
                        ANGLE2=ANGLE2_0
                        dist=2*a
                        write(T,*) 'k= ',k,'a=',a,'b=',b,'r/b=',ri(i,j)
                        WRITE(T,*) 'x   ','TT-'//rstr//'b'
                    END IF

                    do while (f.LE.0.5) 
                        call ellipses(a,b,a,b,ANGLE, ANGLE2, ANGLE3
     +                   ,dist)
                        
                        U=-A1*cos(2*ANGLE+2*ANGLE2)*( b /
     +                 ( ri(i,j)*b - A2*dist + A3*b ))**A4
                        WRITE(T,*) f,U
                        ANGINC=ANGINC+pi/180 !rand 
                        ANGLE2=ANGLE2+pi/180 !radianes
                        f=ANGINC/pi    
                    end do 


                end do
            end do

         end PROGRAM

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
     +           (lambda1-a11)**2)*(1.d0-e1*k1d**2))
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
     +          gamma*in))/2.)        
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
            Rc=cdsqrt((qq**2-1.d0)/deltap*(1.d0+b2p*(1.d0+deltap)/qq)
     +        **2+(1.d0-(qq**2-1.d0)/deltap)*(1.d0+b2p/qq)**2)
                    
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