      Subroutine cnlyr(ddsdx)
	
	!     Boundary layer along centerbody surface
	
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
  
              Common /bnlyr/xsum(100), ysum(100), delshd(100), thetas(100), cfshd(100), xcne(100), ycne(100), delcne(100), &
              &thetac(100), cfcne(100), reyprm, pop, aop, vop, reysec, pos, aos, vos, pex, xscale

        ave(x1, x2) = (x1+x2)/2.0
        pr(am) = (am/(1.0+0.2*am*am))**4
        funm(g, ph) = sqrt(2.0/(g-1.0)*(ph**(-(g-1.0)/g)-1.0))
        funt(g, am) = (1.0+(g-1.0)/2.0*am*am)**(-1.0)
        funp(g, am) = (1.0+(g-1.0)/2.0*am*am)**(-g/(g-1.0))
        funu(temp) = 2.27*32.17*temp**1.5/(198.6+temp)*1.E-8

        write(*,*)"call cnlyr"
        If (reyprm==0.0) Goto 20
        s = 0.702
        i = icone
        beta = 0.0
        twto = 0.950
        If (fdim==1.0) beta = 1.25
        amcne = funm(gamp, pcone(i))
        If (i>1) amcne = funm(gamp, ave(pcone(i),pcone(i-1)))
        pp = pop*funp(gamp, amcne)
        tp = top*funt(gamp, amcne)
        ap = 49.02*sqrt(tp)
        tpto = funt(gamp, amcne)
        taveto = 0.50*twto + 0.22*s**(1.0/3.0) + (0.50-0.22*s**(1.0/3.0))*tpto
        tptave = tpto/taveto
        tave = tp/tptave
        vave = 53.3*tp*funu(tave)/pp
        If (icone>1) Goto 10
        delcne(1) = delshd(1)
        xrex = xscale*delcne(i)/(0.046*(1.0+0.80*amcne**2)**(0.44))
        rex = (aop/vop*amcne*(1.0+0.20*amcne**2)**(-2.25))**(-0.20)
		xdim = (xrex/rex)**(5.0/4.0)
        xcne(1) = xdim*ycone(i)**beta*pr(amcne)
        Goto 12
		
10		yave = ave(ycone(i), ycone(i-1))
        pave = ave(pcone(i), pcone(i-1))
        amave = funm(gamp, pave)
        delx = xcone(i) - xcone(i-1)
        xcne(i) = xcne(i-1) + yave**beta*pr(amave)*xscale*delx
12		xdim = 1.0/(ycone(icone)**beta*pr(amcne))*xcne(i)
        rex = aop/vop*xdim*amcne*(1.0+0.20*amcne**2)**(-2.25)
        If (rex<=0.0) Goto 18
        If (rex>=1.E7)Goto 14
		delcne(i) = 0.046*xdim/xscale*(1.0+0.80*amcne**2)**(0.44)*rex**(-1.0/5.0)
        thetac(i) = 0.036*xdim/xscale*(1.0+0.10*amcne**2)**(-0.70)*rex**(-1.0/5.0)
		ddsdx = 0.0368*(1.0+0.80*amcne**2)**(0.44)*rex**(-1.0/5.0)
		Goto 16
14		delcne(i) = 0.028*xdim/xscale*(1.0+0.80*amcne**2)**(0.44)*rex**(-1.0/6.0)
        thetac(i) = 0.022*xdim/xscale*(1.0+0.10*amcne**2)**(-0.70)*rex**(-1.0/6.0)
		ddsdx = 0.0233*(1.0+0.80*amcne**2)**(0.44)*rex**(-1.0/6.0)

16		rey=ap/vave*thetac(i)*xscale*amcne
        hc=delcne(i)/thetac(i)
        hi=(hc/twto-0.20*amcne**2)*tpto
        cfcne(i)=0.246*exp(-1.561*hi)*rey**(-0.268)*(tptave)**1.268
        Goto 20
18		icone=icone-1
        Call outlyr
        Call exit
20		 write(*,*)"end cnlyr"
        Return

      !      boundary layer along shroud wall

600		Format('1',//28X,18A4,//)       
602		Format(//32X,'xcne',8X,'ycne',8X,'amc',8X,'delcne',6X,'thetac',6X,'thetac',6X,'cfcne',7X,'ddsdx'//)        
604		Format('j')        
606		Format(25X,7F12.5)
			
	End Subroutine cnlyr
	