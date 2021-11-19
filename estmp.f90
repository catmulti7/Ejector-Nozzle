      Subroutine estmp
      !
      !     Estimation of total pressure ratio for choked secondary flow
      !			����������ѹ�ȹ���
	
      !     if choke=-1.0 solution has not been found
      !     if choke=0.0 solution has been found,secondary flow choked
      !     if choke=1.0 solution has been found,secondary flow unchoked
      !
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
        
        ! WRITE(*,*)"call estmp"
        ptoL=0.0001
        atol=1.05
        delp=1.02
        aref=1.001
        hshp=pts(niter)
        area(niter)=asass(islp)
        Do i=1, islp
          area(niter)=amin1(area(niter),asass(i))
        End Do
        If(point/=1.0) Goto 10
        If(area(niter)>=1.0 .And. area(niter)<=atol) Goto 22
10		If(niter>1) Goto 12
        If(point==-1.0) hshp=hshp*delp
        If(point==1.0) hshp=hshp/delp
        Goto 26
12		nlow=0
        nhigh=0
        plow=0.0
        pmin=1.0
        Do i=1, niter
          If(area(i)>=1.0) Goto 14
          nlow=1+nlow
          plow=amax1(plow,pts(i)) !amax1����?
          Goto 16
14		  nhigh=1+nhigh
          pmin=amin1(pmin,pts(i))
          If(pmin==pts(i)) amin=area(i)
16		End Do
        If(nhigh>1) Goto 18
        If(abs(pmin-plow)<=ptol*pmin) Goto 24
        If(nhigh==0 .And. nlow>=2) delp=delp*delp
        If(nhigh==0) hshp=plow*delp
		If(nhigh==1) hshp=pmin/delp
        If(hshp<=plow) hshp=pmin-(pmin-plow)/4.0
        Goto 26
18		If(abs(pmin-plow)<=ptol*pmin) Goto 24
        pmax=1.0
        Do 20, i=1, niter
          If(area(i)<1.0) Goto 20
          If(pts(i)==pmin) Goto 20
          pmax=amin1(pts(i),pmax)
          if(pmax==pts(i)) amax=area(i)
20		Continue
          dmin=amin*(1.0+amin/2.0)
          dmax=amax*(1.0+amax/2.0)
          alpha=(pmin*dmax-pmax*dmin)/(dmax-dmin)
          gamma=(pmax-pmin)/(dmax-dmin)
          hshp=alpha+gamma*aref*(1.0+aref/2.0)
          hlow=(1.0+ptol)*plow
          hmin=(1.0-ptol)*pmin
          If(hshp <= hlow .Or. hshp>=hmin) hshp=pmin-(pmin-plow)/4.0
          Goto 26
22		  hshp=hshp
          choke=0.0
          pts(niter+1)=hshp
          Goto 26
24		  hshp=pmin
          choke=1.0
          pts(niter+1)=hshp
! 26   WRITE(*,*)"end estmp" 
 26 Return
600		  Format('1')
					
        End Subroutine estmp