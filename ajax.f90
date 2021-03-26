      Subroutine ajax(xp, yp, alpha, asec, dadx)

      !     Calculation of flow area in secondary passage
      !         二次通道流面积的计算
	
	      Common x(2, 100), y(2, 100), p(2, 100), t(2, 100) !两行100列矩阵
	      Common xslp(100), yslp(100), amp(100), theta(100), php(100), ams(100), phs(100), asass(100), dasdx(100), islp
	      Common xis(21, 26), yis(21, 26), w(21), tau(26), nsonic, nangle
  
	      Common xsonic(26), ysonic(26), psonic(26), tsonic(26), isonic
	      Common xcone(100), ycone(100), pcone(100), tcone(100), icone
  
	      Common xshd(100), yshd(100), dysdx(100), nshd
	      Common xbdy(100), ybdy(100), dybdx(100), nbdy

Common wtfl, hshp, tos, top, gams, gamp, fung, amr, angr, apref, assaps, xprim, yprim, aprim, dprim, cfl, cvl, cona, dbdy, dshd,&
	  &edn, pamb, yratio, pi, conva, convr, fdim, ndata, nstop, solve, choke, change, charge, typ, point, stag
  
	  Common pts(25), area(25), wleak(25), title(18), niter, try
  
	  Common /bnlyr/xsum(100), ysum(100), delshd(100), thetas(100), cfshd(100), xcne(100), ycne(100), delcne(100), thetac(100),&
	  &cfcne(100), reyprm, pop, aop, vop, reysec, pos, aos, vos, pex, xscale
	  Common /prflr/xprf(6), yprf(6, 25), qprf(6, 25), nprf(6), max, nmax
	  Common /cplot/iplot, xstart, ystart, xspan, yspan, scale, span, axis, xorgn, yorgn, xshft, yshft, kkk(14), pp(14), xdown(100),&
	  &yacros(100)
	  Common /outpn/ipnch, iprnt, icomp
        ave(x1, x2) = (x1+x2)/2.0
        i = islp
        iter = 0
        alpha = 0.0
        error = 0.05
        xp = xslp(i)
        dypdx = tan(theta(i))
        Call finde(xp, yp, a, 2.0)
        If (islp==1) pest = 0.0        
10		iter = 1 + iter
        xref = xp
        yref = yp
        b = -1.0/ave(a, dypdx)
        xp = (yslp(i)-yref+a*xref-b*xslp(i))/(a-b)
        Call finde(xp, yp, a, 2.0)
        test = abs(yp-yref)
        If (iter<100) Goto 12
        If (pest==0.0) Write (7, 600)(title(k), k=1, 18)
        Write (7, 602) xslp(islp), yslp(islp)
        pest = 1.0
        Goto 14        
12		If (test>error*yp) Goto 10        
14		Call shlyr(ddsdx)
        yp = yp - delshd(islp)
        a = a - ddsdx
        alpha = atan(ave(a,dypdx))
        asec = (yp*yp**fdim-yslp(i)*yslp(i)**fdim)/cos(alpha)
        dadx = (1.0+fdim)*(yp**fdim*a-yslp(i)**fdim*dypdx)/cos(alpha)
        Return
		
600		Format ('1', //28X, 18A4, //)        
602		Format (//25X, 'unable to obtain convergence in subroutine ajax,xslp=', F8.5,', yslp=',F8.5)
      end Subroutine 
	
