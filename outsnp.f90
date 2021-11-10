	Subroutine outsnp
	!
	!     write-out ofconditions along sonic line
	!
	  Common x(2, 100), y(2, 100), p(2, 100), t(2, 100) !����100�о���
	  Common xslp(100), yslp(100), amp(100), theta(100), php(100), ams(100), phs(100), asass(100), dasdx(100), islp
	  Common xis(21, 26), yis(21, 26), w(21), tau(26), nsonic, nangle
  
	  Common xsonic(26), ysonic(26), psonic(26), tsonic(26), isonic
	  Common xcone(100), ycone(100), pcone(100), tcone(100), icone
  
	  Common xshd(100), yshd(100), dysdx(100), nshd
	  Common xbdy(100), ybdy(100), dybdx(100), nbdy
  
	  Common wtfl, hshp, tos, top, gams, gamp, fung, amr, angr, apref, assaps, xprim, yprim, aprim, dprim, cfl, cvl, cona, dbdy,&
	  & dshd, edn, pamb, yratio, pi, conva, convr, fdim, ndata, nstop, solve, choke, change, charge, typ, point, stag
  
	  Common pts(25), area(25), wleak(25), title(18), niter, try
  
	  Common /bnlyr/xsum(100), ysum(100), delshd(100), thetas(100), cfshd(100), xcne(100), ycne(100), delcne(100), &
	  &thetac(100), cfcne(100), reyprm, pop, aop, vop, reysec, pos, aos, vos, pex, xscale	
	  
	  Real k1, k2	  
	  ave(x1, x2) = (x1+x2)/2.0
	  
	  write(*,*)"call outsnp"
	  pex = 7.0
	  cfl = 0.0
	  cvl = 0.0
	  aprm = 0.0
	  pave = psonic(1)
	  Write (7, 600)(title(k), k=1, 18)
	  Do i = 1, isonic
		alpha = conva*tsonic(i)!��i��1�ģ�
		Write (7, 604) xsonic(i), ysonic(i), alpha
		If (i==1) Goto 10
		yave = ave(ysonic(i), ysonic(i-1))
		tave = ave(tsonic(i), tsonic(i-1))
		dx = xsonic(i-1) - xsonic(i)
		dy = ysonic(i-1) - ysonic(i)
		aprm = aprm + (1.0+fdim)*yave**fdim*dy
		cvl = cvl + (1.0+fdim)*yave**fdim*cos(tave)*(cos(tave)*dy-sin(tave)*dx)
		cfl = cfl + (1.0+fdim)*yave**fdim*(cos(tave)*dy-sin(tave)*dx)	  
10	  End Do
	  cvl = cvl/cfl
	  cfl = cfl/aprm
	  If (reyprm==0.0) Goto 12
	  k1 = 4.0*0.370/(1.0+pex)/cos(abs(angr))
	  k2 = 4.0*0.370*pex/((1.0+pex)*(2.0+pex))
	  cfl = (1.0-k1*reyprm**(-0.20))*cfl
	  cvl = (1.0-k2*reyprm**(-0.20))*cvl
12	 write(*,*)"call outsnp"
 Return
	  
600	  Format ('1', //33X, 18A4, //)	  
604	  Format (35X, 'xsonic=',f9.5, 8X, 'ysonic=',f8.5, 8X, 'tsonic=',f10.5)
	End Subroutine outsnp