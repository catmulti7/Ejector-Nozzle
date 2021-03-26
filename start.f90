
	Subroutine start
!
!     Calculation of pertinent ejector parameters
!			�й������������ļ���
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
	  
	Real k1, k2
	
	funa(g,am)=((g+1.0)/2.0)**(-(g+1.0)/(2.0*(g-1.0)))*1.0/am*(1.0+(g-1.0)/2.0*am*am)**((g+1.0)/(2.0*(g-1.0)))
	fgam(g)=sqrt(g)*(g+1.0/2.0)**(-(g+1.0)/(2.0*(g-1.0)))

!
!     Subroutine start
!
  edn=edn-xprim
  If(angr>0.0) angr=0.0
  
  Write(7,600)(title(i),i=1,18)
  Write(7,607) dprim
  Write(7,602) amr
  Write(7,608) angr
  Write(7,613) yratio
  Write(7,611) xprim
  If(nshd==1) Write(7,610) dshd
  Write(7,612) edn
  Write(7,615) reyprm
  Write(7,616) tos
  Write(7,617) gams
  Write(7,618) top
  Write(7,619) gamp
  
	If(fdim==0.0) Write(7,620)
	If(fdim==1.0) Write(7,622)
	Write(7,624) nshd
	If(abs(cona)==0.0) Write(7,626) nbdy
	If(abs(cona)>0.0) Write(7,628)
	Write(7,630) ndata
	lchk=45
	pi=3.1415927
	convr=0.01745329
	conva=1.0/convr
	angr=convr*angr
	dpref=dprim
	yprim=1.0
	If(nshd==1) xshd(1)=0.0
	If(nshd==1) yshd(1)=dshd/2.0
  Do 10 i=1,nshd
    xshd(i)=2.0*(xshd(i)-xprim)/dpref
    yshd(i)=2.0*yshd(i)/dpref
10 continue
	If(nshd<=1) Goto 14
	Call spline(xshd,yshd,nshd,dysdx,d2fdx2) 
	line=0
	Write(7,632)
  Do 12, i=1,nshd
    line=1+line
    angle=conva*atan(dysdx(i))
    Write(7,634) i,xshd(i),yshd(i),dysdx(i),angle
    If(line<lchk) Goto 12
    line=0
    If(i<nshd) Write(7,632)
12 Continue
14 If(abs(cona)>0.0) Call cone
  If(nbdy==0) Goto 20
  Do i=1,nbdy
    xbdy(i)=2.0*(xbdy(i)-xprim)/dpref
    ybdy(i)=2.0*ybdy(i)/dpref
  End Do
  Call spline(xbdy,ybdy,nbdy,dybdx,d2fdx2)
  line=0
  Write(7,636)
  Do 18, i=1,nbdy
    line=1+line
    angle=conva*atan(dybdx(i))
    Write(7,638) i,xbdy(i),ybdy(i),dybdx(i),angle
    If(line<lchk) Goto 18
    line=0
    If(i<nbdy) Write(7,636)
18 Continue
20 xprim=0.0
  Call finde(xprim,ybody,slope,1.0)
  aprim=(yprim+ybody)*(yprim-ybody)**fdim*cos(atan(slope))
  edn=2.0*edn/dpref
  apref=funa(gamp,amr)
  fung=fgam(gamp)/fgam(gams)
  xscale=dprim/24.0
  If(dprim==2.0) xscale=6.0/24.0
  cfl=1.0
  cvl=1.0
22  Return

!
!     Format statements
!	
	
600 Format('1',//27X,18A4) !fortran77�е�H������ʾ�ַ������*H��������,//��������
602 Format(//35X,'primary nozzle exit mach no., amr = ',F9.5)
607 Format(//35X,'primary nozzle exit diameter, dprim = ',F9.5)
608 Format(//35X,'primary nozzle lip angle, angr = ',F10.5)
610 Format(//35X,'diameter of the cyclindrical shroud ,dshd = ',F9.5)
611 Format(//35X,'location of primary nozzle, xprim = ',F10.5)
  
612 format(//35X,'ejector length measured from primary nozzle, edn= ', F9.5)
613 Format(//35X,'primary nozzle radius ratio, yratio = ',F8.5)
615 Format(//35X,'primary nozzle reynolds number, reyprm = ',E10.3)
616 Format(//35X,'total temperature of secondary flow, tos = ',F9.3)
  
617 format(//35X,'ratio of specific heats for secondary flow, gams= ',F8.5) 
618 Format(//35X,'total tempeprature of primary flow, top = ' ,F9.3)
	
619 format (//35X,'ratio of specific heats for primary flow, gamp =' ,F8.5) 
620 format (//35X,'flow in ejector is two dimensional, fdim = 0.0')	
622 format (//35X,'flow in ejector is axiisymmetric, fdim= 1.0')	
		
624 Format(//35X,'mnumber of points specifing,shroud contour, nshd = ',I3)
626 Format(//35X,'number of points specifing plug contour, nbdy = ' ,I3)
628 Format(//35X,'number of pints specifing plug contour, nbdy =45 ')       !*****
630 Format(//35X,'number of points specifing enterence conditions, ndata = ',I3)
632 Format('1',//51X,'shroud contour specifcations'//)
		
634 Format(19X,'i =',I3,4X,'xshd(i) =',F9.5,4X,'yshd(i) =', F8.5,4X,'dysdx(i) =', F9.5,4X,'angle(i) =', F9.4)
636 Format('1',//52X,'plug contour specifications'//)
638 Format(19X,'i=', I3,4X,'xbdy(i) =', F8.5,4X,'ybdy(i) =', F8.5,4X,'dybdx(i) =', F9.5,4X,'angle(i) =', F9.4)

	End Subroutine start