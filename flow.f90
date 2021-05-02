      Subroutine flow(j)



      !      Calculation of flow conditions along the slipline

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
        funp(g, am) = (1.0+(g-1.0)/2.0*am*am)**(-g/(g-1.0))
        ave(x1, x2) = (x1+x2)/2.0
                        
        
        write(*,*)"call flow"
        i=islp
        iter=0
        !atol=1.05
        error=0.0001
        asass(i)=1.0
        If(islp==1) amin=10.0
        If(stag>=0.0) Goto 10
        p(2,j)=hshp
        php(i)=hshp
        amp(i)=funm(gamp,php(i))
        ams(i)=0.0
        phs(i)=1.0
        asass(i)=500.0
        Call shlyr(ddsdx)
        If(i==1) theta(i)=pmer(amr,angr,amp(i),gamp) 
        Goto 20
10		iter=1+iter
        asave=asass(i)
        write(*,*)"amp(islp)=",amp(i)
        If(i==1) theta(i)=pmer(amr,angr,amp(i),gamp)
        If(change==1.0) Goto 12
        Call ajax(xp,yp,alpha,asec,dadx)
        asass(i)=asec/aprim*apref/assaps !(36)
        If(i==1) Goto 14
        If(asass(i)>1.05 .And. asass(i-1)>1.07) Goto 14
        ams(i)=ams(i-1)
12		change=1.0
        Call ajax(xp,yp,alpha,asec,dadx)
        asass(i)=ave(asass(i-1),funa(gams,ams(i)))
        asec=asass(i)*aprim/apref*assaps
        yslp(i)=yp*yp**fdim-asec*cos(alpha)
        If(yslp(i)<=0.0) asass(i)=0.5
        If(yslp(i)>0.0 .And. fdim==1.0) yslp(i)=sqrt(yslp(i))
        y(2,j)=yslp(i)
        dzdx=tan(t(1,j-1))
        dydx=(y(2,j)-y(1,j-1))/(x(2,j)-x(1,j-1))
        theta(i)=atan(ave(dydx,dzdx))
        If(asass(i)>1.0) Goto 16
14		If(asass(i)<1.0) Goto 20
        Call astar(ams(i),asass(i),gams)
        phs(i)=funp(gams,ams(i))
        php(i)=phs(i)*hshp
        amp(i)=funm(gamp,php(i))
        p(2,j)=php(i)
        If(i>1) Goto 20
16		test=abs(asass(i)-asave)
        If(iter<25) Goto 18
        Write(7,600)
        Call exit
18		If(test>error*asass(i)) Goto 10
20		point=0.0
        If(asass(i)<1.0) point=-1.0
        If(islp>1 .Or. point==0.0) Goto 22
        asass(i)=0.50
        If(niter==50) Call exit
 22 write(*,*)"end flow"
         Return
				
600		Format('1',//40X,'unable to obtain convengence in subroutine flow',////)
602		Format('j')
				
      End Subroutine flow