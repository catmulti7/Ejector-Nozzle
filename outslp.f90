Subroutine outslp
      !
      !     write-out of pertinent information in secondary flow field
      !			�����������е������Ϣ���
	
		
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
	Common /outpn/ipnch, iprnt, icomp

 
	funa(g, am) = ((g+1.0)/2.0)**(-(g+1.0)/(2.0*(g-1.0)))*1.0/am*(1.0+(g-1.0)/2.0*am*am)**((g+1.0)/(2.0*(g-1.0)))
	funm(g, ph) = sqrt(2.0/(g-1.0)*(ph**(-(g-1.0)/g)-1.0))
      !
      !     Subroutine outslp
	  !
	write(*,*)"call outslp"
	
	write(*,*)"xslp(islp)=",xslp(islp)
	write(*,*)"edn=",edn
	If (xslp(islp)<=edn) Goto 8
	ratio = php(islp)/phs(islp)
	dx = xslp(islp) - xslp(islp-1)
	dy = yslp(islp) - yslp(islp-1)
	dp = php(islp) - php(islp-1)
	dt = theta(islp) - theta(islp-1)
	xslp(islp) = edn
	yslp(islp) = yslp(islp-1) + (dy/dx)*(xslp(islp)-xslp(islp-1))
	php(islp) = php(islp-1) + (dp/dx)*(xslp(islp)-xslp(islp-1))
	theta(islp) = theta(islp-1) + (dt/dx)*(xslp(islp)-xslp(islp-1))
	amp(islp) = funm(gamp, php(islp))
	phs(islp) = php(islp)/ratio
	ams(islp) = funm(gams, phs(islp))
	asass(islp) = funa(gams, ams(islp))
	Call shlyr(ddsdx)        
8	Write (7, 600)(title(k), k=1, 18)
	Write (7, 602) wtfl, wleak(niter), pts(niter), area(niter)
	Write (7, 604)
	Do i = 1, islp
		If (i<=45 .Or. i>=47) Goto 10
		Write (7, 600)(title(k), k=1, 18)
		Write (7, 602) wtfl, wleak(niter), pts(niter), area(niter)
		Write (7, 604)
10		angle = conva*theta(i)
		Write (7, 606) xslp(i), yslp(i), amp(i), angle, php(i), ams(i), phs(i), asass(i)
	End Do
	If (choke==-1.0) Goto 20
	If (icomp==2) Call outprf
	Call outlyr
18	If (nbdy>0) Call outcne
	Call outsnp
	hshp = pts(niter+1)
	assaps = fung*wtfl/pts(niter+1)
	Write (7, 600)(title(k), k=1, 18)
	Write (7, 608) wtfl
	Write (7, 610) pts(niter+1)
	Write (7, 612) assaps
	Write (7, 614) cfl
	Write (7, 616) cvl
	pts(niter) = pts(niter+1)
	If (nshd==1) Call perf
	If (nshd>1 .And. xslp(islp)==edn) Call perf        
20	write(*,*)"end outslp"
		Return
      !     
      !     Format statements
      !        
600		Format ('1', //28X, 18A4, //)        
602		Format (25X, 'wtfl=',f9.6, 5X, 'wleak=',f9.6, 5X, 'pts/ptp=',f9.6, 5X, 'as/as*=',f9.6)        
604		Format (//20X, 'xslp', 8X, 'yslp', 8X, 'amp', 8X, 'theta',7X, 'p/ptp', 8X, 'ams', 8X, 'p/pts', 7X, 'as/as*')       
606		Format (13X, 8F12.5)       
608		Format (//37X, 'secondary corrected weight flow ratio,wtfl=',F9.6)       
610		Format (//37X, 'secondary total pressure ratio,pts/ptp=',f9.6)       
612		Format (//37X, 'secondary critical area patio,as*/ap*=',f9.6)        
614		Format (//37X, 'primary nozzle flow coefficient,cfl=',f8.5)        
616		Format (//37X, 'primary nozzle velocity coefficient,cvl=',f8.5)
      End Subroutine outslp