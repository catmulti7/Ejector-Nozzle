      Subroutine slid(j)

      !      Calculation of wake effect on flow field
	!			��������Ӱ��ļ���

	      Common x(2, 100), y(2, 100), p(2, 100), t(2, 100) !����100�о���
	      Common xslp(100), yslp(100), amp(100), theta(100), php(100), ams(100), phs(100), asass(100), dasdx(100), islp
	      Common xis(21, 26), yis(21, 26), w(21), tau(26), nsonic, nangle
		
	      Common xsonic(26), ysonic(26), psonic(26), tsonic(26), isonic
	      Common xcone(100), ycone(100), pcone(100), tcone(100), icone
  
	      Common xshd(100), yshd(100), dysdx(100), nshd
	      Common xbdy(100), ybdy(100), dybdx(100), nbdy
  
              Common wtfl, hshp, tos, top, gams, gamp, fung, amr, angr, apref, assaps, xprim, yprim, aprim, dprim, cfl, cvl, cona, &
              &dbdy, dshd, edn, pamb, yratio, pi, conva, convr, fdim, ndata, nstop, solve, choke, change, charge, typ, point, stag
         
	      Common pts(25), area(25), wleak(25), title(18), niter, try
				
        funm(g, ph) = sqrt(2.0/(g-1.0)*(ph**(-(g-1.0)/g)-1.0))
        funp(g, am) = (1.0+(g-1.0)/2.0*am*am)**(-g/(g-1.0))
        ave(x1, x2) = (x1+x2)/2.0
				
        iter=0
        error=0.1
        If(islp==2) npoint=0
        If(solve==0.0 .Or. stag==-1.0) Goto 18
        If(charge==1.0 .Or. change==-1.0) Goto 18
        amsave=amp(islp)
        psave=php(islp)
        tsave=theta(islp)
        ams(islp)=ams(islp-1)
        Call ajax(xp,yp,alpha,asec,dadx)
10  iter=iter+1
        asave=asass(islp)
        wsave=wsec	!ǰ�������ޡ�wsec���ñ������Ƿ���asec?
        If(iter==2) asass(islp)=wsec/wtfl*asass(islp)
        If(iter>2) asass(islp)=asass(islp)+dadw*(wtfl-wsec)
        If(asass(islp)<1.0) Goto 16
        Call astar(ams(islp),asass(islp),gams)
        phs(islp)=funp(gams,ams(islp))
        php(islp)=phs(islp)*hshp
        amp(islp)=funm(gamp,php(islp))
        Call wake(delw,dela)
        aspref=(asec-dela)/aprim*apref/asass(islp)
        ws=1.0/fung*hshp*aspref
        wsec=ws+delw
        wratio=ws/wtfl
        wmix=delw/wtfl
        test=abs(wsec-wtfl)
        If(dela/asec>=0.50) Goto 16
        If(iter>1) dadw=(asave-asass(islp))/(wsave-wsec)
        If(iter<45) Goto 12
        If(iter==45) Write(7,600)
        Write(7,602) iter,asass(islp),wsec,delw,test
        If(iter==50) Call exit
12   If(test>0.10*error) Goto 10
        asave=asass(islp)
        If(area(1)<1.20) asave=area(niter)
        If(aspref>assaps) Goto 16
        If(asave>area(niter) .And. area(niter)<1.20) Goto 16
        assaps=aspref
        p(2,j)=php(islp)
        t(2,j)=pmer(amsave,tsave,amp(islp),gamp)
14  Call store(j)
        Goto 18
16  charge=1.0
        Call slip(j)
18   write(*,*)"point=",point
write(*,*)"npoint=",npoint
If(point/=-1.0) Goto 22
        If(stag==1.0) Goto 20
        If(npoint<=1) Goto 24
20  stag=2.0
        point=1.0
        islp=islp-1
        hshp=pts(niter)
        Goto 24
22  If(asass(islp)<=1.05) npoint=1+npoint
        Call ajax(xp,yp,alpha,asec,dasdx(islp))
        dadx=-1.0
        point=0.0
        dela=0.0
        If(wtfl==0.0) Goto 21
        If(islp<3) Goto 24
        delp=asass(islp)-asass(islp-1)
        delq=asass(islp-1)-asass(islp-2)
        dela=-abs(delp)
        If(delp>0.0 .And. delq>0.0) dela=delp
21  If(dasdx(islp)>0.0 .And. dasdx(islp-1)>0.0) dadx=1.0
        If(dela>=0.0 .And. dadx>0.0) point=1.0
        If(try==1.0 .And. islp<=ndata) point=0.0
        If(xslp(islp)<edn .And. area(niter)<=1.05) point=0.0
        If(xslp(islp)>=edn) point=1.0
        If(solve<=1.0 .Or. stag==1.0) Goto 24
        If(dela>0.0 .And. dadx>0.0) point=1.0
        If(try==1.0 .And. islp<=ndata) point=0.0
24		Return
				
600		Format('1', //43X,'unable to obtain convergence in subroutine slid',////)
602		Format(22X,'ITER=',I6,4X,'asass=',F8.5,4X,'wsec=',F8.5,4X,'delw=',F8.5,4X,'TEST=',1PE12.5)
604		Format('1',//28X,12A6,////)
606		Format(22X,'ISLP=',I6,4X,'XSLP=',F8.5,4X,'WSEC=',F8.5,4X,'WMIX=',F8.5,4X,'assaps',F8.5)
608		Format('j')				
	End Subroutine slid
	