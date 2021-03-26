			Subroutine sonpt(k)
	
	!		Construction of sonic point 
	!			���ٵ�Ĺ���
	
			Common x(2, 100), y(2, 100), p(2, 100), t(2, 100) !����100�о���
			Common xslp(100), yslp(100), amp(100), theta(100), php(100), ams(100), phs(100), asass(100), dasdx(100), islp
			Common xis(21, 26), yis(21, 26), w(21), tau(26), nsonic, nangle
			
			Common xsonic(26), ysonic(26), psonic(26), tsonic(26), isonic
			Common xcone(100), ycone(100), pcone(100), tcone(100), icone
  
			Common xshd(100), yshd(100), dysdx(100), nshd
			Common xbdy(100), ybdy(100), dybdx(100), nbdy
			
                  Common wtfl, hshp, tos, top, gams, gamp, fung, amr, angr, apref, assaps, xprim, yprim, aprim, dprim, cfl, &
                  &cvl, cona, dbdy, dshd, edn, pamb, yratio, pi, conva, convr, fdim, ndata, nstop, solve, choke, change,&
                  & charge, typ, point, stag
  
			Common pts(25), area(25), wleak(25), title(18), niter, try
			
			Dimension c(18), d(16), ang(2)

	!		Subroutine snopt
      
      ! WRITE(*,*)"call sonpt"
      iter = 0
      error = 0.15
      Do j = 1,100
        If (p(1,j)==0.0) Goto 12
        Call fonics(k, j, ang(2))
        If (j>1 .And. ang(2)>=t(2,1)) Goto 14
        ang(1) = ang(2)
	 End Do
12	 Write (7, 604)
      Write (*,*)"failed "
      Call outslp
      Call outsnp
      Call exit
14	i = 0
      jmin = j - 1
      Do 18, j = 1, 100
        If (p(1,j)==0.0) Goto 20
        If (j<jmin) Goto 16
        i = 1 + i
        x(1, i) = x(1, j)
        y(1, i) = y(1, j)
        p(1, i) = p(1, j)
        t(1, i) = t(1, j)
        If (i==j) Goto 18
16		x(1, j) = 0.0
        y(1, j) = 0.0
        p(1, j) = 0.0
        t(1, j) = 0.0
18		Continue
20		iter = 1 + iter
      Call coave(xave, yave, pave, tave, 1, 1, 1, 2)
      Call coeff(xave, yave, pave, tave, c, gamp)
      d(1) = c(1)
      d(3) = 1.0/(c(3)*c(5)*c(7))
      d(5) = fdim*c(9)*c(11)/(c(13)*c(15))
      dxda = (x(1,1)-x(1,2))/(ang(1)-ang(2))
      dtda = (t(1,1)-t(1,2))/(ang(1)-ang(2))
      x(1, 1) = x(1, 2) + dxda*(t(2,1)-ang(2))
      y(1, 1) = y(1, 2) + d(1)*(x(1,1)-x(1,2))
      t(1, 1) = t(1, 2) + dtda*(t(2,1)-ang(2))
      p(1, 1) = (d(3)*p(1,2)-(t(1,1)-t(1,2))-d(5)*(y(1,1)-y(1,2)))/d(3)
      ! WRITE(*,*)"p(1,1)",p(1,1)
      ! WRITE(*,*)"psonic(isonic)",psonic(isonic)
      If (p(1,1)>psonic(isonic)) Goto 12
      Call fonics(k, 1, ang(1))
      test = abs(t(2,1)-ang(1))
    !  WRITE(*,*)"test=",test
      If (iter<45) Goto 22
      If (iter==45) Write (7, 600)
      Write (7, 602) iter, x(2, 1), y(2, 1), p(2, 1), ang(1), test
      If (iter==50) Call exit
22		If (test>error) Goto 20
      isonic = 1 + isonic
      xsonic(isonic) = x(2, 1)
      ysonic(isonic) = y(2, 1)
      psonic(isonic) = p(2, 1)
      tsonic(isonic) = t(2, 1)
!24	 WRITE(*,*)"end sonpt"
 24 Return
	  
600		format ('1',//40X,'unable to obtain convergence in subroutine sonpt',////)
602		format (8X,'iter =', I3,4X,'x(2,1)=', f8.5,4X,'y(2,1)=', f8.5,4X,'p(2,1)=', f8.5,4X,'t(2,1)=', f8.5,4X,'test=', 1PE12.5)!ע��IPE12.5			
604		format ('1',//40X,'unable to obtain a solution in subroutine sonpt',////) 
	End Subroutine sonpt