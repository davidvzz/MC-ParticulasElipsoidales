    !**************************************
    !subroutine to calculate the distance of closest approach of two ellipses
    !***************************************
    Subroutine ellipses(a1,b1,a2,b2,k1k2,k1d,k2d,dist)
      Implicit None
      complex*32,parameter::in=(1.q0,0.q0)
      !k1k2, k1d,k2d are dot product between the vectors, k1, k2 and d
      real*16,Intent(in)::a1,b1,a2,b2,k1k2,k1d,k2d 
      real*16,Intent(out)::dist
      real*16::e1,e2,eta,a11,a12,a22
      real*16::lambda1,lambda2,a2p,b2p
      real*16::kpmp,t,deltap,A,B,C,D,E,alpha,beta,gamma,P,Q,Rp
      Complex*32::U,y,qq,z,mysqrt
      Complex*32, External::ccbrt
      
      !eccentricity of ellipses
      e1=1.q0-b1**2/a1**2
      e2=1.q0-b2**2/a2**2
      !component of A'
      eta=a1/b1-1.q0
      a11=b1**2/b2**2*(1.q0+0.5q0*(1.q0+k1k2)*(eta*(2.q0+eta)-e2*(1.q0+eta*k1k2)**2))
      a12=b1**2/b2**2*0.5q0*qsqrt(1.q0-k1k2**2)*(eta*(2.q0+eta)+e2*(1.q0-eta**2*k1k2**2))
      a22=b1**2/b2**2*(1.q0+0.5q0*(1.q0-k1k2)*(eta*(2.q0+eta)-e2*(1.q0-eta*k1k2)**2))
      !eigenvalues of A'
      lambda1=0.5q0*(a11+a22)+0.5q0*qsqrt((a11-a22)**2+4.q0*a12**2)
      lambda2=0.5q0*(a11+a22)-0.5q0*qsqrt((a11-a22)**2+4.q0*a12**2)
      !major and minor axes of transformed ellipse
      b2p=1.q0/qsqrt(lambda1)
      a2p=1.q0/qsqrt(lambda2)
  
      deltap=a2p**2/b2p**2-1.q0
      If(qabs(k1k2)==1.q0) then
          if(a11>a22)then
               kpmp=1.q0/qsqrt(1.q0-e1*k1d**2)*b1/a1*k1d
          elseif(a11<a22) then
              kpmp=qsqrt(1.q0-k1d**2)/qsqrt(1.q0-e1*k1d**2)
          end if
       Else
              kpmp=(a12/qsqrt(1.q0+k1k2)*(b1/a1*k1d+k2d+(b1/a1-1.q0)*k1d*k1k2)+ &
              (lambda1-a11)/qsqrt(1.q0-k1k2)*(b1/a1*k1d-k2d-(b1/a1-1.q0)*k1d*k1k2))&
              /qsqrt(2.q0*(a12**2+(lambda1-a11)**2)*(1.q0-e1*k1d**2))
       End If
  
      IF(kpmp==0.q0 .or. deltap==0.0q0) Then
          Rp=a2p+1.q0
      ELSE
    !coefficients of quartic for q
          t=1.q0/kpmp**2-1.q0
          A=-1.q0/b2p**2*(1.q0+t)
          B=-2.q0/b2p*(1.q0+t+deltap)
          C=-t-(1.q0+deltap)**2+1.q0/b2p**2*(1.q0+t+deltap*t)
          D=2.q0/b2p*(1.q0+t)*(1.q0+deltap)
          E=(1.q0+t+deltap)*(1.q0+deltap)
          
          !solution for quartic 
          alpha=-3.q0/8.q0*(B/A)**2+C/A
          beta=(B/A)**3.q0/8.q0-(B/A)*(C/A)/2.q0+D/A
          gamma=-3.q0/256.q0*(B/A)**4+C/A*(B/A)**2/16.q0-(B/A)*(D/A)/4.+E/A
         
          If(beta==0.) Then
              qq=-B/4.q0/A+cqsqrt((-alpha+cqsqrt(alpha**2-4.q0*gamma*in))/2.)        
          Else
              P=-alpha**2/12.q0-gamma
              Q=-alpha**3/108.q0+gamma*alpha/3.q0-beta**2/8.q0
              U=ccbrt(-0.5q0*Q+cqsqrt(Q**2/4.q0+P**3/27.q0*in))
               
  
              if(cqabs(U)/=0.0q0) then
                  y=-5.q0/6q0*alpha+U-P/3.q0/U
              else
                  y=-5.q0/6.q0*alpha-ccbrt(Q*in)
              end if
              
              qq=-B/4.q0/A+0.5q0*(cqsqrt(alpha+2.q0*y)+cqsqrt(-(3.q0*alpha+2.q0*y+2.q0*beta/cqsqrt(alpha+2.q0*y))))
                  
                 
         End If
  
      !substitute for R'
          Rp=cqsqrt((qq**2-1.q0)/deltap*(1.q0+b2p*(1.q0+deltap)/qq)**2+(1.q0-(qq**2-1.q0)/deltap)*(1.q0+b2p/qq)**2)
          
      END IF
      
          dist=Rp*b1/qsqrt(1.q0-e1*k1d**2)
          
          
      End Subroutine ellipses
     !*************************************
     !Calculate the cubic root of a complex number, return the principal value
     !*************************************
      Complex*32 function ccbrt(z)
      Implicit None
      Complex*32:: z
      Real*16::arg,theta,r
      Real*16, Intrinsic::qatan2
      arg=qatan2(qimag(z),qreal(z))
      theta = arg / 3.0q0
      r = (qreal(z)**2+qimag(z)**2)**(1.q0/6.q0)
      ccbrt = qcmplx(r*qcos(theta), r*qsin(theta))
      End function ccbrt
  