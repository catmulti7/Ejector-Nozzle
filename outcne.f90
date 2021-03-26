	Subroutine outcne
	!
	!     write-out of pertinent information along the ejector centerbody
	!			��������������������Ϣ���

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

	Real mach
	funm(g, ph) = sqrt(2.0/(g-1.0)*(ph**(-(g-1.0)/g)-1.0))

	lchk = 45
	line = 0
	Write (7, 600)(title(k), k=1, 18)
	Do i = 1, icone
		mach = funm(gamp, pcone(i))
		Write (7, 602) xcone(i), ycone(i), mach, tcone(i), pcone(i)
		line = 1 + line
		If (line<=lchk) Goto 10
		If (i<icone) Write (7, 600)(title(k), k=1, 12)	  
10	End Do
	Return
	
600	Format('1',//28X,18A4,//)	  	     
602	Format(21X,'xc=',f8.5,4X,'yc=',f8.5,4x,'mach no.=',f8.5,4x,'theta=',f8.5,4X,'p/hp=',f8.5)
	End Subroutine outcne
	