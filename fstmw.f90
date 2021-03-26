      Subroutine estmw
      !
      !     Estimation of total pressure ratio for wleak=wtfl
      !			wleak=wtflʱ��ѹ�ȹ���
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
				
        ave(x1,x2)=(x1+x2)/2.0
				
        nlow=0
        nmax=0
        nhigh=0
        plow=0.0
        pmax=1.0
        phigh=1.0
        delp=1.02
        Do i=1, niter
          If(abs(wleak(i))>0.0) Goto 10
          nmax=1+nmax
          pmax=amin1(pts(i),pmax)		!�˴�amini�Ƿ���Fortran�ڲ�����amin1()
10		End Do
        Do i=1, niter
          If(wleak(i)==0.0) Goto 14
          If(wleak(i)>wtfl) Goto 12
          If(pts(i)>pmax) Goto 14
          nlow=1+nlow
          plow=amax1(pts(i),plow)
          If(plow==pts(i)) wlow=wleak(i)
          Goto 14          
12		  If(pts(i)>pmax) Goto 14
          nhigh=1+nhigh
          phigh=amin1(pts(i),phigh)		!�˴�amini�Ƿ���Fortran�ڲ�����amin1()
          If(phigh==pts(i)) whigh=wleak(i)        
14		End Do
16		If(nlow>0 .And. nhigh>0) Goto 18
        If(nmax>=0 .And. nhigh==0) hshp=pts(niter)*delp
        If(nmax>=0 .And. nlow==0) hshp=pts(niter)*delp
        If(nmax>0 .And. nlow>0) hshp=pmax-(pmax-plow)/2.0
        Goto 22        
18		If(abs(whigh-wlow)<=0.001) Goto 20
        hshp=phigh-(phigh-plow)/4.0
        Goto 22        
20		dpdw=(phigh-plow)/(whigh-wlow)
        hshp=plow+dpdw*(wtfl-wlow)
22		Return
      End Subroutine estmw