	Subroutine outlyr
	!
	!     Write-out of pertinent boundary layer information
	!			д����صı߽����Ϣ
	  Common x(2, 100), y(2, 100), p(2, 100), t(2, 100) !����100�о���
	  Common xslp(100), yslp(100), amp(100), theta(100), php(100), ams(100), phs(100), asass(100), dasdx(100), islp
	  Common xis(21, 26), yis(21, 26), w(21), tau(26), nsonic, nangle
  
	  Common xsonic(26), ysonic(26), psonic(26), tsonic(26), isonic
	  Common xcone(100), ycone(100), pcone(100), tcone(100), icone
  
	  Common xshd(100), yshd(100), dysdx(100), nshd
	  Common xbdy(100), ybdy(100), dybdx(100), nbdy
  
	  Common wtfl, hshp, tos, top, gams, gamp, fung, amr, angr, apref, assaps, xprim, yprim, aprim, dprim, cfl, cvl, cona, dbdy, &
	  &dshd, edn, pamb, yratio, pi, conva, convr, fdim, ndata, nstop, solve, choke, change, charge, typ, point, stag
  
	  Common pts(25), area(25), wleak(25), title(18), niter, try
  
	  Common /bnlyr/xsum(100), ysum(100), delshd(100), thetas(100), cfshd(100), xcne(100), ycne(100), delcne(100), &
	  &thetac(100), cfcne(100), reyprm, pop, aop, vop, reysec, pos, aos, vos, pex, xscale
	  Common /outpn/ipnch, iprnt, icomp

	  Dimension zslp(100), zcone(100), amcne(100)
	  funm(g, ph) = sqrt(2.0/(g-1.0)*(ph**(-(g-1.0)/g)-1.0))
	!
	!     Subroutine outlyr
	!
	  write(*,*)"call outlyr"
	  If (reyprm==0.0) Goto 18
	  Write (7, 600)(title(k), k=1, 18)
	  Write (7, 602) wtfl, pts(niter), reyprm, reysec
	  Write (7, 604)
	  Do i = 1, islp
		zslp(i) = xslp(i) + xsum(1)/xscale
		Write (7, 606) xslp(i), ysum(i), ams(i), delshd(i), thetas(i), cfshd(i)
		If (i/=45) Goto 10
		Write (7, 600)(title(k), k=1, 12)
		Write (7, 602) wtfl, pts(niter), reyprm, reysec
		Write (7, 604)	  
10	  End Do
	  If (nbdy==0) Goto 14
	  istop = icone - 1
	  Write (7, 600)(title(k), k=1, 18)
	  Write (7, 602) wtfl, pts(niter), reyprm, reysec
	  Write (7, 608)
	  Do i = 1, istop
		zcone(i) = xcone(i) + xcne(i)/xscale  !�ڶ���xcne()��1��i���������i
		amcne(i) = funm(gamp, pcone(i))
		Write (7, 606) xcone(i), ycone(i), amcne(i), delcne(i), thetac(i), cfcne(i)
		If(i/=46) Goto 12
		Write(7,600) (title(k),k=1,18)
		Write(7,602) wtfl, pts(niter), reyprm, reysec
		Write(7,608)	  
12	  End Do	    
14	  If(ipnch==0) Goto 18
	!
	!     Punch cards for sasman cresci turbulent boundary layer program
	!				��˹��cresci�����߽�����Ĵ�׿�
	  pr=0.705
	  recov=0.950
	  dim=(2.0+fdim)
	  scale=dprim/2.0
	  div=4.0
	  iprof=0
	  icard=0
	  iplot=1
	  If(ipnch==2) Goto 16
	  p0s=pos/144.0
	  theta1=scale*thetas(1)
	  h1=delshd(1)/thetas(1)
	  !punch 500, (title(k) , k=1,9)
	  !punch 502, pos, tos, amhu0, gams, pr, recov
	  !punch 502, dim, scale, ams(1) theta1, h1, div
	  !punch 504, islp, islp, tprof, icard, iplot
	  !punch 506, (zslp(i),ysum(i),ams(i),phs(i),delshd(i),thetas(i),cfshd(i),i=1,islp)	  
16	  If(ipnch==1 .Or. nbdy==0) Goto 18
	  p0p=pop/144.0
	  theta1=scale*thetac(1)
	  h1=delcne(1)/thetac(1)
	  !punch 500, (title(k),k=1,9)
	  !punch 502, p0p, top, amhu0, gamp, pr, recov
	  !punch 502, dim, scale, amcne(1) theta1, h1, div
	  !punch 504, istop, istop, iprof, icard, iplot
	  !punch 506, (zcone(i),ycone(i),amcne(i),pcone(i),delcne(i),thetac(i),cfcne(i),i=1,istop)	  
18	  write(*,*)"end outlyr"
	Return	
	  
600	Format('1',//28X,18A4,//)	  
602	Format(22X,'wtfl=',f9.6,5X,'pts/pts=',f9.6,5X,'reyprm=',E10.3,5X,'reysec=',E10.3)	  
604	Format(//32X,'xshd',8X,'yshd', 8X, 'ams', 8X, 'delshd', 6X, 'thetas', 6X, 'cfshd')	  
606	Format (25X, 6F12.5)	  
608	Format (//32X, 'xcne', 8X, 'ycne', 8X, 'amc', 8X, 'delcne', 6X, 'thetac', 6X, 'cfcne')	  
610	Format ('j')
	  
500	Format (9A6)	  
502	Format (6F12.6)	  
504	Format (5I6)	  
506	Format (2F10.6, 10X, 5F10.6)
	End Subroutine outlyr