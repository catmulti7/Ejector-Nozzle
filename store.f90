      Subroutine store(j)
      !
      !     Storage of pertinent information along slipline
      !				�ش��������Ϣ�Ĵ洢
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
				
        funa(g, am) = ((g+1.0)/2.0)**(-(g+1.0)/(2.0*(g-1.0)))*1.0/am*(1.0+(g-1.0)/2.0*am*am)**((g+1.0)/(2.0*(g-1.0)))
        funm(g, ph) = sqrt(2.0/(g-1.0)*(ph**(-(g-1.0)/g)-1.0))
        xslp(islp) = x(2, j)
        yslp(islp) = y(2, j)
        php(islp) = p(2, j)
        theta(islp) = t(2, j)
        amp(islp) = funm(gamp, php(islp))
        If (stag>=0.0) Goto 10
        phs(islp) = 1.0
        ams(islp) = 0.0
        asass(islp) = 500.0
        Goto 12        
10		phs(islp) = php(islp)/hshp
		ams(islp) = funm(gams, phs(islp))    
		asass(islp) = funa(gams, ams(islp))
12		Return
      End Subroutine store