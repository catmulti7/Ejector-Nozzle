      Subroutine comp(k4, skip)
	
      !     Completion of flow field for an ejector with a centerbody
      !     带中心体的引射流场完成结果
	
	      Common x(2, 100), y(2, 100), p(2, 100), t(2, 100) !两行100列矩阵
	      Common xslp(100), yslp(100), amp(100), theta(100), php(100), ams(100), phs(100), asass(100), dasdx(100), islp
	      Common xis(21, 26), yis(21, 26), w(21), tau(26), nsonic, nangle
  
	      Common xsonic(26), ysonic(26), psonic(26), tsonic(26), isonic
	      Common xcone(100), ycone(100), pcone(100), tcone(100), icone
  
	      Common xshd(100), yshd(100), dysdx(100), nshd
	      Common xbdy(100), ybdy(100), dybdx(100), nbdy
  
        Common wtfl, hshp, tos, top, gams, gamp, fung, amr, angr, apref, assaps, xprim, yprim, aprim, dprim, cfl, cvl, cona,&
        & dbdy, dshd, edn, pamb, yratio, pi, conva, convr, fdim, ndata, nstop, solve, choke, change, charge, typ, point, stag
  
	      Common pts(25), area(25), wleak(25), title(18), niter, try
  
	      Common /outpn/ipnch, iprnt, icomp
				
        Real k4
        Dimension c(16), d(6)

      !     Subroutine comp

        write(*,*)"call comp"
        Do k = 1, 100
          x(2, k) = 0.0
          y(2, k) = 0.0
          p(2, k) = 0.0
          t(2, k) = 0.0
          If (p(1,k)>0.0) j = k
        End Do
        If (xslp(islp)<edn) Goto 10
        delx = edn - xslp(islp-1)
        dx = xslp(islp) - xslp(islp-1)
        dy = yslp(islp) - yslp(islp-1)
        dp = php(islp) - php(islp-1)
        dt = theta(islp) - theta(islp-1)	!NASA附录有可能打错了，thvta(islp)应为theat(islp)
        x(1, j) = xslp(islp-1) + delx
        y(1, j) = yslp(islp-1) + dy/dx*delx
        p(1, j) = php(islp-1) + dp/dx*delx
        t(1, j) = theta(islp-1) + dt/dx*delx        
10		xcomp = xslp(islp)
		
        If (nbdy>0) xcomp = xbdy(nbdy)
        Do i = 1, 100
          Call bound(2)
          Do j = 3, 100
            If (p(1,j)==0.0) Goto 14
            Call field(j)
            Call check(j, shock)
            If (x(2,j)>xcomp .Or. shock==1.0) Goto 14
		  End Do         
14		  Call clear(1, j)
          radius = sqrt((x(1,2)-x(1,1))**2+(y(1,2)-y(1,1))**2)
          nsert = radius/k4
          If (p(1,2)==0.0) nsert = 0
          If (nsert>=2) Call insert(nsert-1)
          If (skip==1.0) Call outfld(1)
          If (icomp==2) Call profle
          xtest = x(1, 1) + 2.0*(x(1,2)-x(1,1))
          If (x(1,1)>xcomp) Goto 18
          If (xtest>xcomp) Goto 18
          If (x(1,2)>=xcomp .Or. p(1,2)==0.0) Goto 18
		End Do
        
18		If (nbdy==0) Goto 20
        icone = icone + 1
        xcone(icone) = xbdy(nbdy)
        ycone(icone) = ybdy(nbdy)
        pcone(icone) = pcone(icone-1) + dp/dx*(xcone(icone)-xcone(icone-1))
        tcone(icone) = atan(dybdx(nbdy))
        
20	 write(*,*)"end comp"
	Return
	End Subroutine comp	