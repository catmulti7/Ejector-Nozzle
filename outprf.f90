	Subroutine outprf
	!
	!     write-out of mach number profiles
	!			д������������
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
	  
	  Common /prflr/xprf(6), yprf(6, 25), qprf(6, 25), nprf(6), max, nmax
	  Common /outpn/ipnch, iprnt, icomp

	  funp(g, am) = (1.0+(g-1.0)/2.0*am*am)**(-g/(g-1.0))
	  funt(g, am) = (1.0+(g-1.0)/2.0*am*am)**(-1)
	  iskip = 0
	  Write (7, 600)(title(j), j=1, 12)
	  !If (kpnch>0) punch 700,(title(i), i=1, 12)
	  Do i = 1, max
		If (nprf(i)==0) Goto 14
		kprf = nprf(i)
		iskip = 1 + iskip
		If (iskip==2) Write (7, 604)
		Do k = 1, kprf
		  pres = funp(gamp, qprf(i,k))
		  temp = funt(gamp, qprf(i,k))
		  Write (7, 602) xprf(i), yprf(i, k), qprf(i, k), pres, temp
		End Do
		If (iskip==1) Goto 12
		If (i==max) Goto 12
		iskip = 0
		Write (7, 600)(title(j), j=1, 12)
12		If (kpnch==0) Goto 14
		!punch 702, xprf(i)                                         
		!punch 704, kprf
		!punch 706, (yprf(i,k), oprf(i,k), k=1, kprf)
14	  End Do
	  Return
	  
600	  Format ('1', //27X, 12A6, ///)	  
602	  Format (19X,'xp=',F9.5,4X,'yp=',F8.5,4X,'mach no.=',f8.5,4X,'P/Pop=',F8.5,4X,'t/top=',F8.5)	  
604	  Format (/'j')	  

700	  Format (12A6)	  
702	  Format ('ejector mach no. flow field at x/r=',F6.3)	  
704	  Format (I3)	  
706	  Format (8F10.5)
	End Subroutine outprf