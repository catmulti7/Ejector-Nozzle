			Subroutine shlyr(ddsdx)
	!			Boundary layer along shroud wall
	!				�����ֱ߽��
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
  
			Common /bnlyr/xsum(100), ysum(100), delshd(100), thetas(100), cfshd(100), xcne(100), ycne(100), delcne(100), thetac(100), &
			&cfcne(100), reyprm, pop, aop, vop, reysec, pos, aos, vos, pex, xscale
			
			ave(x1,x2)=(x1+x2)/2.0
			pr(am)=(am/(1.0+0.2*am*am))**4
			funm(g,ph)=sqrt(2.0/(g-1.0)*(ph**(-(g-1.0)/g)-1.0))
			funt(g,am)=(1.0+(g-1.0)/2.0*am*am)**(-1.0)
			funp(g,am)=(1.0+(g-1.0)/2.0*am*am)**(-g/(g-1.0))
			funu(temp)=2.27*32.17*temp**1.5/(198.6+temp)*1.e-8
			
		!	Subroutine layer
			
			DDSDX=0.0
			IF(reyprm==0.0) Goto 30
			i=islp
			s=0.702
			beta=0.0
			start=0.0
			twto=0.950
			if(fdim==1.0) beta=1.25
			if(wtfl>0.0 .And. islp==1) Goto 12
			if(wtfl>0.0 .And. islp>1)	Goto 14
			if(wtfl==0.0 ) then
			if(ams(i-1)>0.0) Goto 14
			end if
			if(islp==1)	Goto 12
10		xsum(i)=0.0
			delshd(i)=0.0
			thetas(i)=0.0
			cfshd(i)=0.0
			reysec=0.0
12		tp=top*funt(gamp,amr)
			ap=49.02*sqrt(tp)
			vp=2.0*ap/reyprm*xscale*amr
			pp=53.3*tp*funu(tp)/vp
			pop=pp/funp(gamp,amr)
			pos=pop*pts(niter)
			aop=49.02*sqrt(top)
			aos=49.02*sqrt(tos)
			vop=53.3*top*funu(top)/pop
			vos=53.3*top*funu(tos)/pos
14		Call finde(xslp(i), ysum(i), dydx, 2.0)
			if(i>=1 .And. ams(i)==0) Goto 30
			if(i>1 .And. ams(i-1)==0.0) Goto 30
			if(delshd(i)==0.0) start=1.0
			amshd=ams(i)
			if(i>1) amshd=ave(ams(i),ams(i-1))	!i��1����
			ps=pos*funp(gams,amshd)
			ts=tos*funt(gams,amshd)
			as=49.02*sqrt(ts)
			vs=53.3*ts*funu(ts)/ps
			tsto=1.0/(1.0+0.20*amshd**2)
			taveto=0.50*twto+0.22*s**(1.0/3.0)+(0.50-0.22*s**(1.0/3.0))*tsto
			tstave=tsto/taveto
			tave=ts/tstave
			vave=53.3*ts*funu(tave)/ps
			if(islp>1) Goto 16
			reysec=2.0*as/vs*xscale*ams(1)!****ע��ams()����1����i
			xrex=xscale*delshd(i)/(0.046*(1.0+0.80*ams(i)**2)**(0.44))
			rex=(aos/vos*ams(i)*(1.0+0.20*ams(i)**2)**(-2.25))**(-0.20)
			xdim=(xrex/rex)**(5.0/4.0)
			xsum(i)=xdim*ysum(i)**beta*pr(ams(i))
			Goto 18
16		yave=ave(ysum(i),ysum(i-1))
			amave=ave(ams(i),ams(i-1))
			delx=xslp(i)-xslp(i-1)
			xsum(i)=xsum(i-1)+yave**beta*pr(amave)*xscale*delx
18		xdim=1.0/(ysum(i)**beta*pr(amshd))*xsum(i)	
			rex=aos/vos*xdim*amshd*(1.0+0.20*amshd**2)**(-2.25)
			if(rex<0.0) Goto 22
			delshd(i)=0.046*xdim/xscale*(1.0+0.80*amshd**2)**(0.44)*rex**(-1.0/5.0)
			thetas(i)=0.036*xdim/xscale*(1.0+0.10*amshd**2)**(-0.70)*rex**(-1.0/5.0)
			ddsdx=0.0368*(1.0+0.80*amshd**2)**(0.44)*rex**(-1.0/5.0)
			rey=as/vave*thetas(i)*xscale*amshd
			hc=delshd(i)/thetas(i)
			hi=(hc/twto-0.20*amshd**2)*tsto
20		cfshd(i)=0.246*exp(-1.561*hi)*rey**(-0.268)*(tstave)**1.268	
			if(islp==1) Goto 30
			if(delshd(i-1)==0.0) Goto 30
			if(start==1.0) delshd(i)=delshd(i-1)
			if(delshd(i)>5.0*delshd(i-1)) delshd(i)=delshd(i-1)
			Goto 30
22		islp=islp-1
			Call outlyr
			Call exit
				
30		Return 	
			
600		format('1',//28X,18A4,//)
602		format(22X,'WTFL=',F9.6,5X,'PTS/PTS=',F9.6,5X,'REYPRM=',1PE10.3,5X,'REYSEC',1PE10.3)
604		format(//26X,'XSHD',8X,'YSHD',8X,'AMS',8X,'DELSHD',6X,'THETAS',6X,'CFSHD',6X,'DDSDX',//)
606		format(19X,7F12.5)
610		format('j')
			
	End subroutine shlyr
	