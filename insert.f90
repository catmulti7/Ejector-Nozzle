      Subroutine insert(npt)
      !
      !     Insertion of field points
      !			�����ֵ
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
      !
      !     subroutine insert
      !
        npts = 1 + npt
        xnpts = npts
        dx = (x(1,2)-x(1,1))/xnpts
        dy = (y(1,2)-y(1,1))/xnpts
        dt = (t(1,2)-t(1,1))/xnpts
        Do 12,j = 1, 100
		  x(2, j) = x(1, j)
		  y(2, j) = y(1, j)
		  p(2, j) = p(1, j)
		  t(2, j) = t(1, j)
		  Call coave(xave, yave, pave, tave, 1, 1, 1, 2)
		  Call coeff(xavf, yave, pave, tave, c, gamp)
		  d(1) = c(1)
		  d(3) = 1.0/(c(3)*c(5)*c(7))
		  d(5) = fdim*c(9)*c(11)/(c(13)*c(15))
12		Continue
        Do j = 2, npts
          x(1, j) = x(1, j-1) + dx
          y(1, j) = y(1, j-1) + d(1)*(x(1,j)-x(1,j-1))
          t(1, j) = t(1, j-1) + dt
          p(1, j) = (d(3)*p(1,j-1)-(t(1,j)-t(1,j-1))-d(s)*(y(1,j)-y(1,j-1)))/d(3)
        End Do
        i = npts
        Do j = 2, 100
          If (p(2,j)==0.0) Goto 18
          i = 1 + i
          x(1, i) = x(2, j)
          y(1, i) = y(2, j)
          p(1, i) = p(2, j)
          t(1, i) = t(2, j)
16		End Do
18		Do j = 1, 100
          x(2, j) = 0.0
          y(2, j) = 0.0
          p(2, j) = 0.0
          t(2, j) = 0.0
        End Do
        Return
      End Subroutine insert