      Subroutine fonics(k, j, angle)
	
      !	Construction of first point along isocline
		!		�ص�б�ߵĵ�һ����Ĺ���
		  Common x(2, 100), y(2, 100), p(2, 100), t(2, 100) !����100�о���
	      Common xslp(100), yslp(100), amp(100), theta(100), php(100), ams(100), phs(100), asass(100), dasdx(100), islp
	      Common xis(21, 26), yis(21, 26), w(21), tau(26), nsonic, nangle			
				
	      Common xsonic(26), ysonic(26), psonic(26), tsonic(26), isonic
	      Common xcone(100), ycone(100), pcone(100), tcone(100), icone
				
	      Common xshd(100), yshd(100), dysdx(100), nshd
	      Common xbdy(100), ybdy(100), dybdx(100), nbdy
		
         Common wtfl, hshp, tos, top, gams, gamp, fung, amr, angr, apref, assaps, xprim, yprim, aprim, dprim, cfl, cvl, cona,&
         & dbdy, dshd, edn, pamb, yratio, pi, conva, convr, fdim, ndata, nstop, solve, choke, change, charge, typ, point, stag
  
	      Common pts(25), area(25), wleak(25), title(18), niter, try
		
        Dimension a(16), b(16), c(16), d(8)
        Real mach, muave
				
        funmu(am) = asin(1.0/am)
        funm(g, ph) = sqrt(2.0/(g-1.0)*(ph**(-(g-1.0)/g)-1.0))
        funp(g, am) = (1.0+(g-1.0)/2.0*am*am)**(-g/(g-1.0))
        ave(x1, x2) = (x1+x2)/2.0

      !      Subroutine fonics

      !  write(*,*)"call fonics"
        t(2, 1) = tau(k)
        p(2, 1) = funp(gamp, amr)
        mach = funm(gamp, p(1,j))
        msonic = nsonic - 1
        tave = ave(t(2,1), t(1,j))
        muave = ave(funmu(amr), funmu(mach))
        d(2) = tan(tave-muave)
        Do i = 1, msonic
          dydx = (yis(i+1,k)-yis(i,k))/(xis(i+1,k)-xis(i,k))
          ! write(*,*)"dxdy=",dxdy
          ! write(*,*)"d(2)=",d(2)
          ! write(*,*)"yis(i,k)=",yis(i,k)
          ! write(*,*)"y(1,j)=",y(1,j)
          ! write(*,*)"x(1,j)=",x(1,j)
          ! write(*,*)"xis(i,k)=",xis(i,k)
          x(2, 1) = (y(1,j)-yis(i,k)-d(2)*x(1,j)+dydx*xis(i,k))/(dydx-d(2))
          ! write(*,*)"x(2,1)=",x(2,1)
          y(2, 1) = y(1, j) + d(2)*(x(2,1)-x(1,j))
          If (k == nangle) Goto 12
          If (x(2, 1)>=xis(i+1,k) .And. x(2, 1)<=xis(i,k)) Goto 12
          ! write(*,*)"x(2,1)=",x(2,1)
          ! write(*,*)"xis(i+1,k)=",xis(i+1,k)
          ! write(*,*)"xis(i,k)=",xis(i,k)
          ! write(*,*)"k=",k
          ! write(*,*)"nangle=",nangle
        End Do
        Write (7, 600)
        Write (7, 602) xis(nsonic, k), x(2, 1), xis(1, k)
        Call outsnp
        Call exit
12		Call coave(xave, yave, pave, tave, 1, j, 2, 1)
        Call coeff(xave, yave, pave, tave, c, gamp)
        d(2) = c(2)
        d(4) = 1.0/(c(4)*c(6)*c(8))
        d(6) = fdim*c(10)*c(12)/(c(14)*c(16))
        angle = t(1,j) + d(4)*(p(2,1))-(p(1,j))+d(6)*(y(2,1)-y(1,j))
!14	 write(*,*)"end fonics"
14	Return
		
600     format('1',//40X,'Sonic point (XP) lies outside range of isoclines',//)				
602		format (//24X,'xis(nsonic,k)=',f9.5,15x,'xp=',f8.5,15x,'xis(1,k)=',f8.5,//)	
      End Subroutine fonics