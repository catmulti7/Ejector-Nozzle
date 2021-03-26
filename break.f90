      Subroutine break(j)
	
	!     Calculation of wake imlpingment on shroud wall
	!		波撞击到外罩壁面上的计算
	      Common x(2, 100), y(2, 100), p(2, 100), t(2, 100) !两行100列矩阵
	      Common xslp(100), yslp(100), amp(100), theta(100), php(100), ams(100), phs(100), asass(100), dasdx(100), islp
	      Common xis(21, 26), yis(21, 26), w(21), tau(26), nsonic, nangle
  
	      Common xsonic(26), ysonic(26), psonic(26), tsonic(26), isonic
	      Common xcone(100), ycone(100), pcone(100), tcone(100), icone
  
	      Common xshd(100), yshd(100), dysdx(100), nshd
	      Common xbdy(100), ybdy(100), dybdx(100), nbdy
  
              Common wtfl, hshp, tos, top, gams, gamp, fung, amr, angr, apref, assaps, xprim, yprim, aprim, dprim, cfl, cvl, cona,&
              & dbdy,dshd, edn, pamb, yratio, pi, conva, convr, fdim, ndata, nstop, solve, choke, change, charge, typ, point, stag
  
	      Common pts(25), area(25), wleak(25), title(18), niter, try
				
        Real kappa
				
	!Function statement
				
        funa(g, am) = ((g+1.0)/2.0)**(-(g+1.0)/(2.0*(g-1.0)))*1.0/am*(1.0+(g-1.0)/2.0*am*am)**((g+1.0)/(2.0*(g-1.0)))
        funm(g, ph) = sqrt(2.0/(g-1.0)*(ph**(-(g-1.0)/g)-1.0))
        funq(g, am) = sqrt((g+1.0)/2.0*am*am/(1.0+(g-1.0)/2.0*am*am))
        !recv(g, p) = (((g+1.0)*p+(g-1.0))/((g-1.0)*p+(g+1.0)))**(g/(g-1.0))*p**(-1.0/(g-1.0))
        escape(g, pr) = sqrt((g+1.0)/(g-1.0)*(1.0-pr**(-(g-1.0)/g)))
        crocco(g, am) = am*am/(2.0/(g-1.0)+am*am)
        test(g, pr) = 1.0-pr**(-(g-1.0)/g)
        ave(x1, x2) = (x1+x2)/2.0
				
        kappa=1.000
        wtol=0.0005
        ptol=0.0005
10		i=islp-1
        dydx=tan(theta(i))
        dtdx=(theta(i)-theta(i-1))/(xslp(i)-xslp(i-1))
        Call finde(xslp(i),yp,dypdx,2.0)
        xref=xslp(i)+(yp-yslp(i))/(dydx-dypdx)
        Call finde(xref,yref,dyrdx,2.0)
        thetap=theta(i)+dtdx*(xref-xslp(i))
        pref=hshp
        amref=funm(gamp,pref)
        thetar=pmer(amp(i),thetap,amref,gamp)
12		delta=thetar-atan(dyrdx)
        If(delta<=0.0) Goto 18
        If(delta>0.0) pratio=psi(gamp,amref,delta)
        recomp=1.0+kappa*(pratio-1.0)
        pest=test(gamp,recomp)
        If(pest<0.0) Call exit
        ama=funq(gamp,amref)
        amd=escape(gamp,recomp)
        b=(1.0-tos/top)*(amd/ama)**2/2.0
	c=tos/top*(amd/ama)**2
        phid=b+sqrt(b*b+c)
        tostop=tos/top+(1.0-tos/top)*phid
        phib=0.0
				
        gami=(gamp+wtfl*gams)/(1.0+wtfl)
        siisi=(1.0+phib)/(1.0-phib)
        caiisq=crocco(gamp,amref)
        caisq=caiisq*(1.0-phib)**2/(caiisq*(1.0-phib)**2+(1.0-caiisq))
        sigi=12.0+2.758*sqrt(caisq)/(sqrt((1.0-caisq)*(gami-1.0)/2.0))
        sigii=siisi*sigi
        aratio=funa(gamp,amr)/funa(gamp,amref)
        Call mix(tos/top,caiisq,phib,phid,sigvb)
        wleak(niter)=2.0/sigii*xref*yref*aratio*sigvb
				
        xslp(islp)=xslp(islp-1)
        yslp(islp)=yslp(islp-1)
        theta(islp)=atan(dypdx)
        phs(islp)=1.0/recomp
        ams(islp)=funm(gams,phs(islp))
        asass(islp)=funa(gams,ams(islp))
        php(islp)=recomp*pref
        amp(islp)=funm(gamp,php(islp))
        area(niter)=amin1(area(niter),asass(islp))
        point=-1.0
        nhigh=0
        Do 14, i=1,niter
          If(wleak(i)==0.0) Goto 14
          If(wleak(i)>wtfl) nhigh=1+nhigh
14		Continue
        wiff=abs(wtfl-wleak(niter))
        If(niter==1) piff=1.0
        If(niter>1) piff=abs(pts(niter)-pts(niter-1))
        If(nhigh>0 .And. piff<=ptol*pts(niter)) Goto 16
        If(wiff>wtol) Goto 18
16		point=0.0
        typ=0.0
        stag=1.0
        change=1.0
        pts(niter+1)=hshp
        wleak(niter)=wtfl
        asec=yp*yp**fdim-yslp(islp)*yslp(islp)**fdim
        If(abs(dypdx)>0.0) Call ajax(xp,yp,alpha,asec,dadx)
        aspaps=asec/aprim*apref/asass(islp)
        If(wtfl>0.0) assaps=amin1(assaps,aspaps)
        If(wtfl==0.0) assaps=aspaps
        hshp=php(islp)/phs(islp)
        p(1,j-1)=ave(php(islp),php(islp-1))
        t(1,j-1)=ave(theta(islp),theta(islp-1))
        Call field(j-1)
17			islp=islp+1
        Call slip(j)
        Call slid(j)
18		Return
        

600		Format('1'//35X,'wake impinges on shroud wall at xslp=',f9.6,'yslp =',F9.6,//)        
602		Format(24X,'wtfl =',F9.6,4X,'wleak =',F9.6,4X,'pratio =',F9.6,4X,'tos/top =',F9.6)        
604		Format('1')
				
	End Subroutine break

	