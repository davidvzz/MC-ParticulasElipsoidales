!INPUT:  a1, b1 (a1>b1): length of semiaxis of first ellipse
    !        a2, b2 (a2>b1): length of semiaxis of second ellipse
    !        theta1(0, 2*pi): angle associated with the major axis of first ellipse, meaured from lab frame counterclockwise
    !        theta2(0, 2*pi): angle associated with the major axis of second ellipse
    !        theta3(0, 2*pi): angle associated with the vector joining the centers, pointing to the center of the second ellipse
    !        xc1, yc1:        x and y coordinates of the center of the first ellipse
    !OUTPUT: dist: distance between the centers when two ellipses are externally tangent, i.e., distance of closest approach
    !        xt, yt:  x and y coordinates of the contact points in the global lab frame     
    Subroutine ellipses(a1,b1,a2,b2,theta1,theta2,theta3,xc1,yc1,dist,xt,yt)
    Implicit None
    Double Complex,parameter::in=(1.d0,0.d0)
    Double Precision,Intent(in)::a1,b1,a2,b2,theta1,theta2,theta3,xc1,yc1
    Double Precision,Intent(out)::dist,xt,yt
    Double Precision::k1k2,k1d,k2d 
    Double Precision::e1,e2,eta,a11,a12,a22
    Double Precision::lambda1,lambda2,a2p,b2p
    Double Precision::kpmp,kpmps,t,deltap,A,B,C,D,E,alpha,beta,gamma,P,Q,Rc
    Double Precision::sn1,sn2,cs1,cs2,cs3,sn3,cps,sps,rp,rm,ami,bmi,ama,bma
    Double Complex::U,y,qq,z,mysqrt
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
    a11=b1**2/b2**2*(1.d0+0.5d0*(1.d0+k1k2)*(eta*(2.d0+eta)-e2*(1.d0+eta*k1k2)**2))
    a12=b1**2/b2**2*0.5d0*dsqrt(1.d0-k1k2**2)*(eta*(2.d0+eta)+e2*(1.d0-eta**2*k1k2**2))
    a22=b1**2/b2**2*(1.d0+0.5d0*(1.d0-k1k2)*(eta*(2.d0+eta)-e2*(1.d0-eta*k1k2)**2))
 
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
            kpmp=(a12/dsqrt(1.d0+k1k2)*(b1/a1*k1d+k2d+(b1/a1-1.d0)*k1d*k1k2)+ &
            (lambda1-a11)/dsqrt(1.d0-k1k2)*(b1/a1*k1d-k2d-(b1/a1-1.d0)*k1d*k1k2))&
            /dsqrt(2.d0*(a12**2+(lambda1-a11)**2)*(1.d0-e1*k1d**2))
            kpmps=(-(lambda1-a11)/dsqrt(1.d0+k1k2)*(b1/a1*k1d+k2d+(b1/a1-1.d0)*k1d*k1k2)+ &
            a12/dsqrt(1.d0-k1k2)*(b1/a1*k1d-k2d-(b1/a1-1.d0)*k1d*k1k2))&
            /dsqrt(2.d0*(a12**2+(lambda1-a11)**2)*(1.d0-e1*k1d**2))
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
        gamma=-3.d0/256.d0*(B/A)**4+C/A*(B/A)**2/16.d0-(B/A)*(D/A)/4.+E/A
       
        If(beta==0.d0) Then
            qq=-B/4.d0/A+cdsqrt((-alpha+cdsqrt(alpha**2-4.d0*gamma*in))/2.)        
        Else
            P=-alpha**2/12.d0-gamma
            Q=-alpha**3/108.d0+gamma*alpha/3.d0-beta**2/8.d0
            U=ccbrt(-0.5d0*Q+cdsqrt(Q**2/4.d0+P**3/27.d0*in))

            if(abs(U)/=0.0d0) then
                y=-5.d0/6d0*alpha+U-P/3.d0/U
            else
                y=-5.d0/6.d0*alpha-ccbrt(Q)
            end if
            
            qq=-B/4.d0/A+0.5d0*(cdsqrt(alpha+2.d0*y)+cdsqrt(-(3.d0*alpha+2.d0*y+2.d0*beta/cdsqrt(alpha+2.d0*y))))
               
       End If

    !substitute for R'
        Rc=cdsqrt((qq**2-1.d0)/deltap*(1.d0+b2p*(1.d0+deltap)/qq)**2+(1.d0-(qq**2-1.d0)/deltap)*(1.d0+b2p/qq)**2)
        
    END IF
    
! The distance of closest approach

    dist=Rc*b1/dsqrt(1.d0-e1*k1d**2)
        
!	-------now we calculate the coordinates of the point of contact 
	IF((deltap .EQ. 0.D0) .OR. (kpmp .EQ. 0.D0) ) THEN
		xt=xc1+b1/DSQRT(1.D0-e1*cs1**2)*cs3
		yt=yc1+b1/DSQRT(1.D0-e1*cs1**2)*sn3

	ELSE
	sps = CDSQRT((qq*qq-1.D0)/deltap)
	cps = CDSQRT((1.D0-(qq*qq-1.D0)/deltap))   	
!	now make sure that the sign of psi is same as of phi 
 	sps = DSIGN(sps,kpmps)
	cps = DSIGN(cps,kpmp)

     If(DABS(k1k2).EQ. 1.0D0) Then
		if (A11>A22) then
             xt = xc1+a1*cps*cs1-b1*sps*sn1
             yt = yc1+a1*cps*sn1+b1*sps*cs1
		else
	         xt = xc1-b1*cps*sn1+sps*a1*cs1
             yt = yc1+ b1*cps*cs1+sps*a1*sn1
		end if
	 Else
     !general case
	 rp = DSQRT(1.D0+k1k2)
	 rm = DSQRT(1.D0-k1k2)

	 ami = (1.D0/DSQRT(2.D0*(a12*a12+(lambda1-a11)**2)))*a12
	 bmi = (1.D0/DSQRT(2.D0*(a12*a12+(lambda1-a11)**2)))*(lambda1-a11)

	 ama = cps*(ami/rp+bmi/rm) + sps*(ami/rm-bmi/rp)
	 bma = cps*(ami/rp-bmi/rm) - sps*(ami/rm+bmi/rp)

     xt = xc1+ama*a1*cs1+bma*(a1-b1)*k1k2*cs1+bma*b1*cs2 
     yt = yc1+ama*a1*sn1+bma*(a1-b1)*k1k2*sn1+bma*b1*sn2
     End If
     
    END IF
	
    End Subroutine ellipses