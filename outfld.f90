	Subroutine outfld(i)
	!
	!     write-out of pertinent information in primary flow field
	!			�������������Ϣ���
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
	  Common /cplot/iplot, xstart, ystart, xspan, yspan, scale, span, axis, xorgn, yorgn, xshft, yshft, kkk(14), pp(14),&
	  & xdown(100), yacros(100)

	  Real mach
	  funm(g, ph) = sqrt(2.0/(g-1.0)*(ph**(-(g-1.0)/g)-1.0))

!	  WRITE(*,*)"enter outfld"
	  If(iplot==1) Goto 12
	  Write(7,600)
	  Do k=1, 100
		If(p(i,k)==0.0) Goto 12
		mach=funm(gamp,p(i,k))
		Write(7,602) x(i,k), y(i,k), mach, p(i,k), t(i,k)
	  End Do
12	  If(iplot==1) Call plotf(i)
! WRITE(*,*)"quit outfld"
	  Return
	  
600	  Format('1',//48X,'CONDITIONS IN THE PRIMARY FLOW FIELD',//)	  	  
602	  format(24X, 'xp=',F9.5, 4X, 'yp=',F8.5, 4X, 'mach no.=',F8.5,4X, 'p/hp=',f8.5,4X,'THETA=',F9.5)
	End Subroutine outfld