			Subroutine coave(xave, yave, pave, tave, i1, j1, i2, j2)
	!		Averaging subroutine
	!			ƽ���ӳ���
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
			
			Common /clpot/xpen, ypen, nx, ny, ipen, xlabel(10), ylabel(10)
			Common /specl/test, orgset, spaset !ע��˴�����������������
			
			
			ave(x1,x2)=(x1+x2)/2.0
			funp(g, am)=(1.0+(g-1.0)/2.0*am*am)**(-g/(g-1.0))
			
			!write(*,*)"call coave"
			pmax=funp(gamp,1.0)
			xave=ave(x(i1,j1),x(i2,j2))
			yave=ave(y(i1,j1),y(i2,j2))
			pave=ave(p(i1,j1),p(i2,j2))
			tave=ave(t(i1,j1),t(i2,j2))
			if(pave<pmax) Goto 10
			write(7,600)
			write(7,602) xave, yave, pave, tave
			Call outslp
			if(choke==-1.0) Call outlyr
			Call exit

!10		write(*,*)"end coave"
	10	 Return		
600		format('1',//33X,'conditions of sonic flow have been reached in subroutine coave',///)
602		format(29X,'xave=',F9.5,4X,'yave=',F9.5,4X,'pave=',F8.5,4X,'tave=',F8.5,4X)
			
	End subroutine coave
	