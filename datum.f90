    Subroutine datum
    !
    !     calculation of entering flow conditions
    !			������������
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
			
		Dimension d(2), amq(21)
      Real muave
	  
      funmu(am) = asin(1.0/am)
      funm(g, ph) = sqrt(2.0/(g-1.0)*(ph**(-(g-1.0)/g)-1.0))
      funp(g, am) = (1.0+(g-1.0)/2.0*am*am)**(-g/(g-1.0))
      funq(g, am) = sqrt(am*am/(1.0+(g-1.0)/(g+1.0)*(am*am-1.0)))
      funr(g, q) = sqrt(((1.0-(g-1.0)/(g+1.0))*q*q)/(1.0-(g-1.0)/(g+1.0)*q*q))
      funt(g, am) = (1.0+(g-1.0)/2.0*am*am)**(-1.0)
      funu(temp) = 2.27*32.17*temp**1.5/(198.6+temp)*1.e-8
      ave(x1, x2) = (x1+x2)/2.0
    !
    !     subroutine datum
    !
      xp = 0.0
      div = ndata - 1
      isonic = ndata
      nangle = ndata
      Call finde(xp, yp, dypdx, 1.0)
      dely = (1.0-yp)/div
      angle = abs(atan(dypdx))
      icone = 1
      xref = 0.0
      tcne = abs(atan(dypdx))
      delcne(icone) = 0.0
      If (reyprm==0.0) Goto 6
			
      tp = top*funt(gamp, amr)
      ap = 49.02*sqrt(tp)
      vp = 2.0*ap/reyprm*xscale*amr
      pp = 53.3*tp*funu(tp)/vp
      pop = pp/funp(gamp, amr)
      aop = 49.02*sqrt(top)
      vop = 53.3*top*funu(top)/pop
			
6		If (angle==0.0) Goto 8
      ycone(icone) = yp
      pcone(icone) = funp(gamp, amr)
      Call cnlyr(ddsdx)
      yp = yp + delcne(icone)
      dypdx = dypdx + ddsdx
      xref = xp - yp/dypdx
      angle = abs(atan(dypdx))
8	  Do i = 1, ndata
        xsonic(i) = 0.0
        ysonic(i) = 1.0
        psonic(i) = funp(gamp, amr)
        tsonic(i) = atan(dypdx)
        amq(i) = funm(gamp, psonic(i))
      End Do
      If (fdim==0.0 .Or. nbdy==0) Goto 12
      qp = funq(gamp, amr)!qpΪȫ�ٶȴ�С��up��vp�ֱ�Ϊˮƽ�ʹ�ֱ�����ٶ�ֵ
      up = qp*cos(angle)
      vp = qp*sin(angle)
      vup = -1.0/tan(angle)
      Call conic(xref, ysonic(1), up, vp, vup, gamp)
      qp = sqrt(up*up+vp*vp)
      amq(1) = funr(gamp, qp)
      psonic(1) = funp(gamp, amq(1))
      tsonic(1) = -atan(vp/up)
      angr = tsonic(1)
12    Do 18, i = 2, ndata
        If (i==ndata) Goto 16
        xsonic(i) = xsonic(i-1)
        ysonic(i) = ysonic(i-1) - dely
        dyqdx = 0.0
        If (xref>0.0) dyqdx = -ysonic(i)/xref
        If (fdim==0.0 .Or. nbdy==0) Goto 14
        qp = funq(gamp, amr)
        up = qp*cos(angle)
        vp = qp*sin(angle)
        vup = -1.0/tan(angle)
        Call conic(xref, ysonic(i), up, vp, vup, gamp)
        qp = sqrt(up*up+vp*vp)
        amq(i) = funr(gamp, qp)
        psonic(i) = funp(gamp, amp(i))
        tsonic(i) = -atan(vp/up)
14			tave = ave(tsonic(i),tsonic(i-1))
        muave = ave(funmu(amq(i)), funmu(amq(i-1)))
        d(2) = tan(tave-muave)
        xsonic(i) = (ysonic(i-1)-ysonic(i)-d(2)*xsonic(i-1))/(dyqdx-d(2))
        ysonic(i) = ysonic(i-1)+d(2)*(xsonic(i)-xsonic(i-1))
        Goto 18
16			xsonic(i) = xsonic(i-1)
        Call finde(xsonic(i), ysonic(i), dyqdx, 1.0)
        dyqdx = dyqdx +	ddsdx
        ysonic(i) = ysonic(i) + delcne(icone)
        tave = ave(tsonic(i), tsonic(i-1))
        muave = ave(funmu(amq(i)), funmu(amq(i-1)))
        d(2) = tan(tave-muave)
        xsonic(i) = (ysonic(i-1)-ysonic(i)-d(2)*xsonic(i-1)+dyqdx*xsonic(i))/(dyqdx-d(2))
18    Continue
      xcone(icone) = xsonic(ndata) - delcne(icone)*sin(tcne)
      Call finde(xcone(icone), ycone(icone), dypdx, 1.0)            
      pcone(icone) = psonic(ndata)
      tcone(icone) = atan(dypdx)
      amr = funm(gamp, psonic(1))
      Return
600		Format ('1', //28X, 12A6, ///)
    End Subroutine datum