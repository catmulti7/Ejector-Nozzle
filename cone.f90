      Subroutine cone
      !
      !     Calculation of conic plug contour
	    !      Բ׶�������ļ���
      !
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
				
        nbdy = 45
        div = 44.0
        b = -tan(convr*cona)
        a = dbdy/2.0
        x0 = -dprim/4.0
        xr = -a/b
        delx = (xr-x0)/div
        xbdy(1) = x0
        ybdy(1) = a + b*x0
        Do i = 2, nbdy
          xbdy(i) = xbdy(i-1) + delx
          ybdy(i) = a + b*xbdy(i)
        End Do
        Return
      End Subroutine cone