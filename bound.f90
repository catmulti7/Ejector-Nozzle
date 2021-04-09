      Subroutine bound(j)

      !      Centerline or solid boundary calculaiion

	    Common x(2, 100), y(2, 100), p(2, 100), t(2, 100) !����100�о���
	    Common xslp(100), yslp(100), amp(100), theta(100), php(100), ams(100), phs(100), asass(100), dasdx(100), islp
	    Common xis(21, 26), yis(21, 26), w(21), tau(26), nsonic, nangle
  
	    Common xsonic(26), ysonic(26), psonic(26), tsonic(26), isonic
	    Common xcone(100), ycone(100), pcone(100), tcone(100), icone
  
	    Common xshd(100), yshd(100), dysdx(100), nshd
	    Common xbdy(100), ybdy(100), dybdx(100), nbdy
  
            Common wtfl, hshp, tos, top, gams, gamp, fung, amr, angr, apref, assaps, xprim, yprim, aprim, dprim, cfl, cvl, cona,&
             &dbdy,dshd, edn, pamb, yratio, pi, conva, convr, fdim, ndata, nstop, solve, choke, change, charge, typ, point, stag
  
	    Common pts(25), area(25), wleak(25), title(18), niter, try
  
            Common /bnlyr/xsum(100), ysum(100), delshd(100), thetas(100), cfshd(100), xcne(100), ycne(100), delcne(100), &
            &thetac(100),cfcne(100), reyprm, pop, aop, vop, reysec, pos, aos, vos, pex, xscale
            Common /cplot/iplot, xstart, ystart, xspan, yspan, scale, span, axis, xorgn, yorgn, xshft, yshft, kkk(14), &
            &pp(14), xdown(100),yacros(100)
	  
		  Dimension c(16), d(6)
      
			ave(x1,x2)=(x1+x2)/2.0
      !     Subroutine buund

        write(*,*)"call bound"
        iter=0
        error=0.0001
        x(2, j) = x(1, j)
        p(2, j) = p(1, j)
        If (icone == 1 .And. nbdy > 0) Call cnlyr(ddsdx)
        icone = icone + 1
        delcne(icone) = delcne(icone - 1)
        If (nbdy > 0) Goto 10
        ddsdx = 0.0
        delcne(icone) = 0.0
10		iter = 1 + iter
        pref = p(2, j)
        Call finde(x(2, j), ycne(icone) ,dycdx, 1.0)
        tcne = abs(atan(dycdx))
        dypdx = dycdx + ddsdx
        y(2, j) = ycne(icone) + delcne(icone)/cos(tcne)
        t(2, j) = atan(dypdx)
        Call coave(xave, yave, pave, tave, 1, j, 2, j)
        Call coeff(xave, yave, pave, tave, c, gamp)
        d(2) = c(2)
        d(4) = 1.0/(c(4)*c(6)*c(8))
        d(6) = fdim*c(10)*c(12)/(c(14)*c(16))
        x(2, j) = (y(1, j) - y(2, j) - d(2)*x(1, j) + dypdx*x(2, j))/(dypdx - d(2))
        p(2, j) = (d(4)*p(1, j) + t(2, j) - t(1, j) - d(6)*(y(2,j)-y(1,j)))/d(4)
        xcone(icone) = x(2, j) - delcne(icone)*sin(tcne)
        Call finde(xcone(icone), ycone(icone), dycdx, 1.0)
        pcone(icone) = p(2, j)
        tcone(icone) = atan(dycdx)
        If (nbdy > 0 .And. reyprm > 0.0) Call cnlyr(ddsdx)
12		test = abs(p(2, j) - pref)
        If (iter < 45) Goto 14
        If (iter == 45) Write (6, 600)
        Write (6, 602) iter, x(2, j), y(2, j), p(2, j), t(2, j), test
        If (iter < 45) Goto 14
        Call outslp
        If (choke == -1.0) Call outlyr
        If (iplot == 1) Call plotl
        Call exit
14		If (test > error*p(2, j)) Goto 10
16	write(*,*)"call bound"
        Return
				
600		format('1',//,'unable to obtain cunvergence in subroutine line',//)				
602		format(8X,'ITER=',I3,4X,'X(2,j)=',F8.5,4X,'Y(2,j)=',F8.5,4X,'P(2,j)=',F8.5,4X,'T(2,j)=',F8.5,4X,'TEST=',1PE12.5)
604		format(19X,7F12.5)
606		format('j')
				
      End Subroutine bound
