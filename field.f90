      Subroutine field(j)

      !	      Field point calculation
		!			场点计算

	      Common x(2, 100), y(2, 100), p(2, 100), t(2, 100) !两行100列矩阵
	      Common xslp(100), yslp(100), amp(100), theta(100), php(100), ams(100), phs(100), asass(100), dasdx(100), islp
	      Common xis(21, 26), yis(21, 26), w(21), tau(26), nsonic, nangle
  
	      Common xsonic(26), ysonic(26), psonic(26), tsonic(26), isonic
	      Common xcone(100), ycone(100), pcone(100), tcone(100), icone
  
	      Common xshd(100), yshd(100), dysdx(100), nshd
	      Common xbdy(100), ybdy(100), dybdx(100), nbdy
				
              Common wtfl, hshp, tos, top, gams, gamp, fung, amr, angr, apref, assaps, xprim, yprim, aprim, dprim, cfl, cvl, cona,&
              & dbdy, dshd, edn, pamb, yratio, pi, conva, convr, fdim, ndata, nstop, solve, choke, change, charge, typ, point, stag
  
	      Common pts(25), area(25), wleak(25), title(18), niter, try
        Dimension a(16), b(16), d(6)
				
        ave(x1, x2) = (x1+x2)/2.0
				
		!      Subroutine field
                                
        write(*,*)"call field"
        k=j-1
        iter=0
        error=0.1
        x(2,j)=ave(x(1,j),x(2,k))
        y(2,j)=ave(y(1,j),y(2,k))
        p(2,j)=ave(p(1,j),p(2,k))
        t(2,j)=ave(t(1,j),t(2,k))
10		iter=1+iter
        pref=p(2,j)
        Call coave(xave,yave,pave,tave,2, k, 2, j)
        Call coeff(xave,yave,pave,tave,a,gamp)
        Call coave(xave,yave,pave,tave,1, j, 2, j)
        Call coeff(xave,yave,pave,tave,b,gamp)
        d(1)=a(1)
        d(2)=b(2)
        d(3)=1.0/(a(3)*a(5)*a(7))
        d(4)=1.0/(b(4)*b(6)*b(8))
        d(5)=fdim*a(9)*a(11)/(a(13)*a(15))
        d(6)=fdim*b(10)*b(12)/(b(14)*b(16))
        x(2,j)=(y(1,j)-y(2,k)+d(1)*x(2,k)-d(2)*x(1,j))/(d(1)-d(2))
        y(2,j)=y(1,j)+d(2)*(x(2,j)-x(1,j))
        p(2,j)=(d(3)*p(2,k)+d(4)*p(1,j)+t(2,k)-t(1,j)-d(5)*(y(2,j)-y(2,k))-d(6)*(y(2,j)-y(1,j)))/(d(3)+d(4))
        t(2,j)=t(1,j)+d(4)*(p(2,j)-p(1,j))+d(6)*(y(2,j)-y(1,j))
        test=abs(pref-p(2,j))
        If(iter<95) Goto 12
        If(iter==95) Write(7,600)
        Write(7,602) iter,x(2,j),y(2,j),p(2,j),t(2,j),test
        If(iter==100) Call exit
12		If(test>error*p(2,j)) Goto 10
write(*,*)"end field"
        Return
				
600		Format('1', //40X,'unable to obtain cunvergence in subroutine field',////) 	
602		Format(8X,'ITER=',I3,4X,'X(2,j)=',F8.5,4X,'Y(2,j)=',F8.5,4X,'P(2,j)=',F8.5,4X,'T(2,j)=',F8.5,4X,'TEST=',E12.5)

      End Subroutine field
