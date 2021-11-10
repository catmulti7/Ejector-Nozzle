      Subroutine finde(xp, yp, dypdx, surf)
		!
		!     Location of shroud contour point
		!        ����������λ��
		Common x(2, 100), y(2, 100), p(2, 100), t(2, 100) !����100�о���
		Common xslp(100), yslp(100), amp(100), theta(100), php(100), ams(100), phs(100), asass(100), dasdx(100), islp
		Common xis(21, 26), yis(21, 26), w(21), tau(26), nsonic, nangle
  
		Common xsonic(26), ysonic(26), psonic(26), tsonic(26), isonic
		Common xcone(100), ycone(100), pcone(100), tcone(100), icone
  
		Common xshd(100), yshd(100), dysdx(100), nshd
		Common xbdy(100), ybdy(100), dybdx(100), nbdy
  
                Common wtfl, hshp, tos, top, gams, gamp, fung, amr, angr, apref, assaps, xprim, yprim, aprim, dprim, cfl, cvl, &
                &cona, dbdy, dshd, edn, pamb, yratio, pi, conva, convr, fdim, ndata, nstop, solve, choke, change, charge,&
                & typ, point, stag
  
	Common pts(25), area(25), wleak(25), title(18), niter, try
				
	If (surf==1.0) Goto 14
        If (nshd>1) Goto 10
        yp = yshd(1)
        dypdx = 0.0
        Goto 20
        
10   If (xp>=xshd(1) .And. xp<=xshd(nshd)) Goto 12
        Write (7, 600)
        Write (7, 602) xshd(1), xp, xshd(nshd)
        Call exit
        
12   Call sintp(xshd, yshd, nshd, xp, yp)
        Call sintp(xshd, dysdx, nshd, xp, dypdx)
        Goto 20
        
14		If (nbdy>0) Goto 16
        yp = 0.0
        dypdx = 0.0
        Goto 20
        
16		If (xp>=xbdy(1) .And. xp<=xbdy(nbdy)) Goto 18
        Write (7, 604)
        Write (7, 606) xbdy(1), xp, xbdy(nbdy)
        Call exit        
18		Call sintp(xbdy, ybdy, nbdy, xp, yp)
        Call sintp(xbdy, dybdx, nbdy, xp, dypdx)        
20		Return      
600	Format ('1', //37X, 'shroud point (xp) lies outside range of input contour',//)       
602  Format (//24X, 'xshd(1)=',F8.5, 15X, 'xp=',F8.5, 15X, 'xshd(nshd)=',f8.5,//)        
604  Format ('1', //38X, 'body point (xp) lies outside range of input contour',//)       
606  Format (//24X, 'xbdy(1)=',F8.5, 15X, 'xp=',F8.5, 15X, 'xbdy(nbdy)=',F8.5,//)
	End Subroutine finde
	