      Subroutine profle
      !
      !     Computation of initial velocity profile
      !         ��ʼ�ٶ����ߵļ���
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
		
    Common /bnlyr/xsum(100), ysum(100), delshd(100), thetas(100), cfshd(100), xcne(100), ycne(100), delcne(100), &
    &thetac(100), cfcne(100), reyprm, pop, aop, vop, reysec, pos, aos, vos, pex, xscale
		Common /prflr/xprf(6), yprf(6, 25), qprf(6, 25), nprf(6), max, nmax
		Common /outpn/ipnch, iprnt, icomp

        Dimension c(16), d(6), ywall(10, 2)
        funm(g, ph) = sqrt(2.0/(g-1.0)*(ph**(-(g-1.0)/g)-1.0))

        write(*,*)"call profle"
        max = 6
        nmax = 25
        step = max - 1
        If (xprf(6)==0.0) xprf(6) = edn
        delx = (xprf(6)-xprf(1))/step
        If (islp>2) Goto 12
        Do i = 1, max
          nprf(i) = 0
          If (i>1) xprf(i) = xprf(i-1) + delx
		End Do
12		Do i = 1, max
          If (xprf(i)==0.0 .Or. nprf(i)==nmax) Goto 26
          If (xprf(i)<xslp(islp-1) .Or. xprf(i)>xslp(islp)) Goto 14
          If (nprf(i)>0) Goto 16
          nprf(i) = 3
          Call finde(xprf(i), ywall(i,1), dywdx, 1.0)
          Call finde(xprf(i), ywall(i,2), dywdx, 2.0)
          dydx = (yslp(islp)-yslp(islp-1))/(xslp(islp)-xslp(islp-1))
          dpdx = (php(islp)-php(islp-1))/(xslp(islp)-xslp(islp-1))
          dtdx = (theta(islp) - theta(islp-1))/(xslp(islp) - xslp(islp-1))
          yprfl = yslp(islp - 1) + dydx*(xprf(i) - xslp(islp-1))
					
          pprf = php(islp - 1) + dpdx*(xprf(i) - xslp(islp-1))
          tprf = theta(islp - 1) + dtdx*(xprf(i) - xslp(islp-1))
          yprf(i, 1) = 1.0
          yprf(i, 2) = (yprfl - ywall(i, 1))/(ywall(i, 2) - ywall(i, 1))
          yprf(i, 3) = yprf(1, 2)
          qprf(i, 1) = funm(gams, pprf/pts(niter))
          qprf(i, 2) = qprf(i, 1)
          qprf(i, 3) = funm(gamp, pprf)          
14		  If (icone == 1) Goto 16
          If (xcone(icone-1) <= xprf(1) .And. xcone(icone) >= xprf(i)) Goto 22
          If (xcone(icone-1) > xprf(i)) Goto 26          
16		  If (xprf(i) < x(1, 1)) Goto 26
          Do j = 2, 100
            If (p(1, j) == 0.0) Goto 26
            If (x(1, j-1) <= xprf(i) .And. x(1, j) >= xprf(i)) Goto 20
          End Do
          Goto 26          
20		  nprf(i) = 1 + nprf(i)
          iprf = nprf(i)
          dydx = (y(1, j) - y(1, j-1))/(x(1, j) - x(1, j-1))
          dtdx = (t(1, j) - t(1, j-1))/(x(1, j) - x(1, j-1))
          Call coave(xave, yave, pave, tave, 1, j - 1, 1, j)
          Call coeff(xave, yave, pave, tave, c, gamp)
          d(3) = 1.0/(c(3)*c(5)*c(7))
          d(5) = fdim*c(9)*c(11)/(c(13)*c(15))
          yprfl = y(1, j - 1) + dydx*(xprf(i) - x(1, j-1))
          tprf = t(1, j - 1) + dtdx*(xprf(i) - x(1, j-1))
          pprf = (d(3)*p(i, j-1) - (tprf-t(1,j-1))-d(5)*(yprfl-y(1,j-1)))/d(3)
          yprf(i, iprf) = (yprfl - ywall(i, 1))/(ywall(i,2)-ywall(i,1))
          qprf(i, iprf) = funm(gamp, pprf)
          Goto 26          
22		  If (icone>2) Goto 24
          xcone(1) = xsonic(ndata)
          ycone(1) = ysonic(ndata)
          pcone(1) = psonic(ndata)
          tcone(1) = tsonic(ndata)
          If (xprf(i)<xcone(1)) Goto 26          
24		  nprf(i) = 1 + nprf(i)
          iprf = nprf(i)
          dpdx = (pcone(icone)-pcone(icone-1))/(xcone(icone)-xcone(icone-1))
          dtdx = (tcone(icone)-tcone(icone-1))/(xcone(icone)-xcone(icone-1))
          pprf = pcone(icone-1) + dpdx*(xprf(i)-xcone(icone-1))
          tprf = tcone(icone-1) + dtdx*(xprf(i)-xcone(icone-1))
          yprf(i, iprf) = 0.0
          qprf(i, iprf) = funm(gamp, pprf)        
26			End Do
write(*,*)"call profle"
        Return
        
600		Format ('1', //27X, 18A4, //)        
602		Format (20X, 'i =',I3, 5X, 'j =',I3, 5X, 'xprf(i) =',F8.5, 5X, 'yprf(i,j)=',F8.5, 5X, 'qprf(i,j) =',F8.5)
	End Subroutine profle
	