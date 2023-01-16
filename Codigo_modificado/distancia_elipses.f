      Program distancia
        Implicit None
        real*8  :: a,b, k1k2, k1d, k2d ,dist=0.0, ang1, ang2, ang3 
        real*8,parameter :: pi=3.141592654D0

        a=2d0
        b=0.5d0

        ang1=0d0*pi/180d0
        ang2=90d0*pi/180d0
        ang3=270d0*pi/180d0
        
        k1k2=dcos(ang1-ang2)
        k1d=dcos(ang1-ang3)
        k2d=dcos(ang2-ang3)

        call ellipses(a, b, a, b, k1k2, k1d, k2d, dist )
        write(*,*) dist      
      End

      !**************************************
      !subroutine to calculate the distance of closest approach of two ellipses
      !***************************************
      Subroutine ellipses(a1,b1,a2,b2,k1k2,k1d,k2d,dist)
      Implicit None
      complex*16,parameter::in=(1.d0,0.d0)
      !k1k2, k1d,k2d are dot product between the vectors, k1, k2 and d
      real*8,Intent(in)::a1,b1,a2,b2,k1k2,k1d,k2d 
      real*8,Intent(out)::dist
      real*8::e1,e2,eta,a11,a12,a22
      real*8::lambda1,lambda2,a2p,b2p
      real*8::kpmp,t,deltap,A,B,C,D,E,alpha,beta,gamma,P,Q,Rp
      Complex*16::U,y,qq
      Complex*16, External::ccbrt
      
      !eccentricity of ellipses
      e1=1.d0-b1**2/a1**2
      e2=1.d0-b2**2/a2**2
      !component of A'
      eta=a1/b1-1.d0
      a11=b1**2/b2**2*(1.d0+0.5d0*(1.d0+k1k2)*(eta*(2.d0+eta)-e2*(1.d0+
     + eta*k1k2)**2))
      a12=b1**2/b2**2*0.5d0*dsqrt(1.d0-k1k2**2)*(eta*(2.d0+eta)+e2*(1.d0-eta
     + **2*k1k2**2))
      a22=b1**2/b2**2*(1.d0+0.5d0*(1.d0-k1k2)*(eta*(2.d0+eta)-e2*(1.d0-eta
     + *k1k2)**2))
      !eigenvalues of A'
      lambda1=0.5d0*(a11+a22)+0.5d0*dsqrt((a11-a22)**2+4.d0*a12**2)
      lambda2=0.5d0*(a11+a22)-0.5d0*dsqrt((a11-a22)**2+4.d0*a12**2)
      !major and minor axes of transformed ellipse
      b2p=1.d0/dsqrt(lambda1)
      a2p=1.d0/dsqrt(lambda2)
  
      deltap=a2p**2/b2p**2-1.d0
      If(dabs(k1k2)==1.d0) then
          if(a11>a22)then
               kpmp=1.d0/dsqrt(1.d0-e1*k1d**2)*b1/a1*k1d
          elseif(a11<a22) then
              kpmp=dsqrt(1.d0-k1d**2)/dsqrt(1.d0-e1*k1d**2)
          end if
       Else
        kpmp=(a12/dsqrt(1.d0+k1k2)*(b1/a1*k1d+k2d+(b1/a1-1.d0)*k1d*k1k2) 
     +  +(lambda1-a11)/dsqrt(1.d0-k1k2)*(b1/a1*k1d-k2d-(b1/a1-1.d0)*k1d
     +  *k1k2))/dsqrt(2.d0*(a12**2+(lambda1-a11)**2)*(1.d0-e1*k1d**2))
       End If
  
      IF(kpmp==0.d0 .or. deltap==0.0d0) Then
          Rp=a2p+1.d0
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
          gamma=-3.d0/256.d0*(B/A)**4+C/A*(B/A)**2/16.d0-(B/A)*(D/A)
     +     /4.0d0+E/A
         
          If(beta==0.) Then
              qq=-B/4.d0/A+cdsqrt((-alpha+cdsqrt(alpha**2-4.d0*gamma
     +         *in))/2.d0)        
          Else
              P=-alpha**2/12.d0-gamma
              Q=-alpha**3/108.d0+gamma*alpha/3.d0-beta**2/8.d0
              U=ccbrt(-0.5d0*Q+cdsqrt(Q**2/4.d0+P**3/27.d0*in))
               
  
              if(cdabs(U)/=0.0d0) then
                  y=-5.d0/6d0*alpha+U-P/3.d0/U
              else
                  y=-5.d0/6.d0*alpha-ccbrt(Q*in)
              end if
              
              qq=-B/4.d0/A+0.5d0*(cdsqrt(alpha+2.d0*y)+cdsqrt(-(3.d0
     +         *alpha+2.d0*y+2.d0*beta/cdsqrt(alpha+2.d0*y))))
                  
                 
         End If
  
      !substitute for R'
          Rp=cdsqrt((qq**2-1.d0)/deltap*(1.d0+b2p*(1.d0+deltap)/qq)
     +     **2+(1.d0-(qq**2-1.d0)/deltap)*(1.d0+b2p/qq)**2)
          
      END IF
      
      dist=Rp*b1/dsqrt(1.d0-e1*k1d**2)
          
          
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