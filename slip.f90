      Subroutine slip(j)
	      !      Slipline calculation
			!		�����߼���
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
        Dimension c(16), d(6)
			
		funp(g,am)=(1.0+(g-1.0)/2.0*am*am)**(-g/(g-1.0))
		ave(x1,x2)=(x1+x2)/2.0
			  
      !      Subroutine slip
        
        write(*,*)"enter slip"
        k = j - 1
        iter = 0
        error =0.0001
        Call coeff(x(2,k), y(2,k), p(2,k), t(2,k), c, gamp)
        dydx = tan(t(1,k))
10		x(2, j) = (y(2,k)-y(1,k)+dydx*x(1,k)-c(1)*x(2,k))/(dydx-c(1))
        y(2, j) = y(2, k) + c(1)*(x(2,j)-x(2,k))
        p(2, j) = ave(p(1,k), p(2,k))
        t(2, j) = ave(t(1,k), t(2,k))
        If (p(2,j)>=hshp) p(2, j) = php(islp-1)
        Call store(j)
12		iter = 1 + iter
        yref = y(2, j)
        pref = p(2, j)
14		Call coave(xave, yave, pave, tave, 2, j, 2, k)
        Call coeff(xave, yave, pave, tave, c, gamp)
        d(1) = c(1)
        d(3) = 1.0/(c(3)*c(5)*c(7))
        d(5) = fdim*c(9)*c(11)/(c(13)*c(15))
        If (change==1.0) Goto 16
        x(2, j) = (y(2,k)-y(1,k)+dydx*x(1,k)-d(1)*x(2,k))/(dydx-d(1))
        y(2, j) = y(2, k) + d(1)*(x(2,j)-x(2,k))
        Call flow(j)
        If (point == -1.0) Goto 24
        If (change == 1.0) Goto 18
        t(2, j) = t(2, k) - d(3)*(p(2,j)-p(2,k)) - d(5)*(y(2,j)-y(2,k))
        dydx = ave(tan(t(1,k)), tan(t(2,j)))
        Call store(j)
        Goto 20
16		Call flow(j)
        If (point == -1.0) Goto 24
18		y(2, j) = yslp(islp)
        If (abs(yref-y(2,j))<=error*y(2,j)) y(2, j) = yref
        x(2, j) = x(2, k) + (y(2,j)-y(2,k))/d(1)
        dydx = (y(2,j)-y(1,k))/(x(2,j)-x(1,k))
        dzdx = tan(t(1,k))
        t(2, j) = atan(ave(dydx,dzdx))
        p(2, j) = (d(3)*p(2,k)-(t(2,j)-t(2,k))-d(5)*(y(2,j)-y(2,k)))/d(3)
        pest = p(2, j)/hshp
        If (pest>0.0) pest = pest**(-(gams-1.0)/gams) - 1.0
        If (pest>0.0) Call store(j)
        If (pest>0.0) Goto 20
        point = -1.0
        asass(islp) = 0.5
        Goto 24
20		testp = abs(p(2,j)-pref)
        testy = abs(y(2,j)-yref)
        If (iter<95) Goto 22
        If (iter==95) Write (7, 600)
        If (wtfl==0.0) Write (7, 602) iter, x(2, j), y(2, j), p(2, j), t(2, j), testy
        If (wtfl>0.0) Write (7, 602) iter, x(2, j), y(2, j), p(2, j), t(2, j), testp
        If (iter==100) Goto 23
22		If (iter<1) Goto 12
		If (testp>error*p(2,j)) Goto 12
        If (testy>error*y(2,j)) Goto 12
23		If (stag>=0.0) Goto 24
        Call finde(xslp(islp), yp, dypdx, 2.0)
        delta = yp - yslp(islp)
        If (delta <= 0.0) point = -1.0
24		write(*,*)"end slip"
Return

600		format('1', //40X,'unable to obtain convergence in subroutine slip',////) 	
602		format(8X,'ITER=',I3,4X,'X(2,j)=',F8.5,4X,'Y(2,j)=',F8.5,4X,'P(2,j)=',F8.5,4X,'T(2,j)=',F8.5,4X,'TEST=',1PE12.5)
604		format('j')		
      End Subroutine slip