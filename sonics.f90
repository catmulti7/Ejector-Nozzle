      Subroutine sonics

       !	Construction of isoclines for primary flow field
		!		Ϊ��Ҫ��������ȸ���
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
				
        Real match, mamb, maprch

        funm(g, v) = sqrt(2.0/(g+1.0)*v*v/(1.0-(g-1.0)/(g+1.0)*v*v))
        funp(g, am) = (1.0+(g-1.0)/2.0*am*am)**(-g/(g-1.0))
        funr(g, am) = (1.0+(g-1.0)/2.0*am*am)**(-1.0/(g-1.0))
        funq(g,pa)=sqrt((g+1.0)/(g-1.0)*(1.0-pa**((g-1.0)/g)))
        funv(g, am) = sqrt((g+1.0)/2.0*am*am/(1.0+(g-1.0)/2.0*am*am))
        omega(v, am) = alog(2.0*v*am/sqrt(1.0-am*am)/(1.0+sqrt(1.0+v*v*am*am/(1.0-am*am))))
        ave(x1,x2)=(x1+x2)/2.0

      !
      !     Subroutine sonics
      !
        nmax=26
        inner=0
        nsonic=21
        wtflw=0.80
        nangle=ndata
        angle=conva*abs(angr)
        error=1
        match=0.99-0.001*(angle-10.0)
        vmatch=funv(gamp,match)
        vamb=funq(gamp,pamb)
        mamb=funm(gamp,vamb)
        vmax=vamb/vmatch
10      write(*,*)"10"
		inner=1+inner
        iter=0
        velapr=0.20
        rhovel=wtflw/(yratio*yratio**fdim)
        ! write(*,*)"wtflw=",wtflw
        ! write(*,*)"yratio=",yratio
        ! write(*,*)"fdim=",fdim
12	write(*,*)"12"	
iter=1+iter
        vaprch=velapr
        velstr=vaprch*vmatch
        maprch=funm(gamp,velstr)
        velapr=rhovel*funr(gamp,match)/funr(gamp,maprch)
        ! write(*,*)"rhovel changed=",rhovel
        test=abs(velapr-vaprch)
        If(iter<50) Goto 14
        Write(7,600)
        Call exit
14      write(*,*)"14"
! write(*,*)"error*velapr=",error*velapr
!         write(*,*)"test=",test
	If(test>error*velapr) Goto 12
        vmin=0.750
        If(nangle>nmax) nangle=nmax
        If(nangle<=nstop) nangle=1+nstop
        xnoa=nangle-nstop
        xnow=nsonic-1
        wj=omega(vmax,match)
        delw=wj-omega(velapr,match)
        ! write(*,*)"delw=",delw
        ! write(*,*)"velapr=",velapr
        ! write(*,*)"match=",match
        delv=(vmax-vmin)/xnow
        xis(1,1)=0.0
        yis(1,1)=1.0
        w(1)=omega(vmax,match)
        tau(1)=angr
        alpha=abs(angr)
        dt=alpha/xnoa
		!vratio=vmax
        Do i=2,nsonic
          vratio=vratio-delv  !vratio�����״γ���,ǰ��δ������
          If(vratio<1.0) vratio=1.0
          w(i)=omega(vratio,match)
          dw=w(i)-w(i-1)
          Call dzdxdy(w(i),wj,alpha,tau(1),delw,dw,0.0,dxs,dys)
          xis(i,1)=xis(i-1,1)+dxs
          yis(i,1)=yis(i-1,1)+dys
          If(vratio==1.0) Goto 18
		End Do
18		write(*,*)"18"
        wtflw=0.0
        nchnge=nangle-nstop
        Do 20, j=2,nangle
          If(j<=nchnge) dt=dt
          If(j>nchnge) dt=dt/2.0
          If(j==nangle) dt=-tau(j-1)
          tau(j)=tau(j-1)+dt
          avetau=ave(tau(j),tau(j-1))
          Call dzdxdy(w(i),wj,alpha,tau(j),delw,0.0,dt,dxs,dys)
          xis(i,j)=xis(i,j-1)+dxs
          yis(i,j)=yis(i,j-1)+dys
          wtflw=wtflw+(dys*cos(avetau)-dxs*sin(avetau))
20      Continue
        scale=1.0-yis(i,j)
        
        wtflw=-wtflw/scale
        test=abs(rhovel-wtflw/(yratio*yratio**fdim))
        If(inner<25) Goto 22
        Write(6,600)
        Call exit
22		write(*,*)"22"
        If(test>error*rhovel) Goto 10
        vratio=vmax
        Do i=2,nsonic
          vratio=vratio-delv
          w(i)=omega(vratio,match)
          dw=w(i)-w(i-1)
        !   write(*,*)"w(i)=",w(i)
        !   write(*,*)"wj=",wj
        !   write(*,*)"alpha=",alpha
        !   write(*,*)"tau(1)=",tau(1)
        !   write(*,*)"delw=",delw
        !   write(*,*)"dw=",dw
          Call dzdxdy(w(i),wj,alpha,tau(1),delw,dw,0.0,dxs,dys)!ע��tau()����1����i

          xis(i,1)=xis(i-1,1)+dxs/scale
          yis(i,1)=yis(i-1,1)+dys/scale
        !   write(*,*)"xis(i,1)=",xis(i,1)
        !   write(*,*)"yis(i,1)=",yis(i,1)
        End Do

        Do i=1,nsonic
          Do j=2,nangle
			dt=tau(j)-tau(j-1)
            If(i==1 .And. j==nangle) Goto 25
        !     write(*,*)"w(i)=",w(i)
        !     write(*,*)"wj=",wj
        !     write(*,*)"alpha=",alpha
        !     write(*,*)"tau(1)=",tau(1)
        !     write(*,*)"delw=",delw
        !     write(*,*)"dw=",dw
            Call dzdxdy(w(i),wj,alpha,tau(j),delw,0.0,dt,dxs,dys)
25	write(*,*)"25"
        xis(i,j)=xis(i,j-1)+dxs/scale
           yis(i,j)=yis(i,j-1)+dys/scale
        !    write(*,*)"scale",scale
        !    write(*,*)"dxs(i,k)=",dxs
        !   write(*,*)"dys(i,k)=",dys
            If(j==nangle) yis(i,j)=0.0
          End Do
        End Do
        isonic=1
        xsonic(isonic)=xprim
	ysonic(isonic)=yprim
        psonic(isonic)=funp(gamp,amr)
        tsonic(isonic)=tau(isonic)
        Return
		
600		Format ('1',//40x,'unable to obtain convergence in subroutine sonic',////)		
		!�޷����ӳ���sonic����������		
      End Subroutine sonics
	