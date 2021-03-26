      Subroutine perf
      !
      !     Calculation of pertinent ejector parameters
      !       �й������������ļ���(����P37��������������������)
	
		Common x(2, 100), y(2, 100), p(2, 100), t(2, 100) !����100�о���
		Common xslp(100), yslp(100), amp(100), theta(100), php(100), ams(100), phs(100), asass(100), dasdx(100), islp
		Common xis(21, 26), yis(21, 26), w(21), tau(26), nsonic, nangle
  
		Common xsonic(26), ysonic(26), psonic(26), tsonic(26), isonic
		Common xcone(100), ycone(100), pcone(100), tcone(100), icone
  
		Common xshd(100), yshd(100), dysdx(100), nshd
		Common xbdy(100), ybdy(100), dybdx(100), nbdy
  
    Common wtfl, hshp, tos, top, gams, gamp, fung, amr, angr, apref, assaps, xprim, yprim, aprim, dprim, cfl, cvl, cona, dbdy, &
    &dshd, edn, pamb, yratio, pi, conva, convr, fdim, ndata, nstop, solve, choke, change, charge, typ, point, stag
  
		Common pts(25), area(25), wleak(25), title(18), niter, try
  
    Common /bnlyr/xsum(100), ysum(100), delshd(100), thetas(100), cfshd(100), xcne(100), ycne(100), delcne(100), thetac(100), &
    &cfcne(100), reyprm, pop, aop, vop, reysec, pos, aos, vos, pex, xscale

        funa(g, am) = ((g+1.0)/2.0)**(-(g+1.0)/(2.0*(g-1.0)))*1.0/am*(1.0+(g-1.0)/2.0*am*am)**((g+1.0)/(2.0*(g-1.0)))
        funm(g, ph) = sqrt(2.0/(g-1.0)*(ph**(-(g-1.0)/g)-1.0))
        funp(g, am) = (1.0+(g-1.0)/2.0*am*am)**(-g/(g-1.0))
        funv(g, am) = sqrt(am*am/(1.0+(g-1.0)/2.0*am*am))
        fideal(g, ph) = sqrt(2.0*g*g/(g-1.0)*(2.0/(g+1.0))**((g+1.0)/(g-1.0))*(1.0-ph**((g-1.0)/g)))
        ave(x1, x2) = (x1+x2)/2.0
      !
      !     Subroutine perf
      !
        cpl = 1.0
        iref = islp
        cfi = 0.003
        pres = funp(gamp, amr)
        fp = pres*(1.0+gamp*cfl*cvl*amr*amr)
        islp = 1
        Call finde(xslp(1), ywall, dywdx, 2.0)
        Call ajax(xsec, ysec, alpha, aref, dadx)
        If (reyprm==0.0) ysum(1) = ywall
        ywall = ywall - delshd(1)
        asec = (ywall*ywall**fdim-yslp(1)*yslp(1)**fdim)/aprim
        If (assaps>0.0) asass(1) = asec/aprim*apref/assaps
        If (wtfl>0.0) Call astar(ams(1), asass(1), gams)
        php(1)=pts(niter)*funp(gams,ams(1))
        angsec=ave(theta(1), atan(dywdx))
        amsfc=ams(1)*cos(angsec)
        fs=php(1)*(1.0+gams*amsec*amsec)*cpl*asec
        fshd=0.0
        fbdy=0.0
        fdrag=0.0
        Do 10 i=2, iref
          If(reyprm==0.0) Call finde(xslp(1),ysum(i),dywdx,2.0) !xslp()��i��1
          xave=ave(xslp(i),xslp(i-1))
          yave=ave(ysum(i),ysum(i-1))
          If(xslp(1)<=xsec) pave=ave(php(1),php(2))     !!!!�ֲ����i,1
          If(xslp(i)>xsec) pave=ave(php(i),php(i-1))
          cfave=ave(cfshd(i),cfshd(i-1))
          amave=ave(ams(i),ams(i-1))
          amsq=amave*amave
          dx=xslp(i)-xslp(i-1)
          dy=ysum(i)-ysum(i-1)
          da=(ysum(i)*ysum(i)**fdim-ysum(i-1)*ysum(i-1))/aprim
          ds=yave**fdim*sqrt(dx*dx+dy*dy)/aprim
          pda=pave*da
          drag=cfave*gams*pave*amsq*ds
          fshd=fshd+pda
          fdrag=fdrag+drag					
10		End do 
12		If(nbdy==0) Goto 16
        Do i=2,icone
          yave=ave(ycone(i),ycone(i-1))
          pave=ave(pcone(i),pcone(i-1))
          amave=funm(gamp,pave)
          cfave=ave(cfcne(i),cfcne(i-1))
          amsq=amave*amave
          dx=xcone(i-1)-xcone(i)
          dy=ycone(i-1)-ycone(i)
          da=(ycone(i-1)*ycone(i-1)**fdim-ycone(i)*ycone(i)**fdim)/aprim
          ds=yave**fdim*sqrt(dx*dx+dy*dy)/aprim
          pda=pave*da
          drag=cfave*gamp*pave*amsq*ds
          fbdy=fbdy+pda
          fdrag=fdrag+drag
		End Do        
16		ftotal=fs+fp+fshd+fbdy-fdrag
        Call finde(edn,yp,dypdx,2.0)
        aexit=yp*yp**fdim/aprim
        fgross=ftotal-php(iref)*aexit
        hppo=1.0/php(iref)
        Write(7,600) hppo
        Write(7,602) fp
        Write(7,604) fs
        Write(7,606) fshd
        Write(7,608) fbdy
        Write(7,609) fdrag
        Write(7,610) ftotal
        Write(7,612) fgross
        Write(7,618)
        Write(7,620)
        hppo=1.0
        Do i=2,40
          hppo=1.0+hppo
          pohp=1.0/hppo
          pohs=pohp/pts(niter)
          fip=cfl*fideal(gamp,pohp)/apref
          If(pohs>=1.0) fis=0.0
          If(pohs<1.0) fis=pts(niter)*assaps*fideal(gams,pohs)/apref
          fgross=ftotal-pohp*aexit
          If(fgross<0.0) Goto 18
          cvi=fgross/(fip+fis)
          cvp=fgross/fip
          Write(7,622) hppo,fgross,fip,fis,cvp,cvi        
18		End Do
		Return
      !
      !     Format statements
      !        
600		Format(//37X,'nozzle pressure ratio,ptp/p0=',f10.5)        
602		Format(//37X,'primary stream thrust,fp/(ptp*ap)=',f9.6)        
604		Format(//37X,'secondary stream thrust,fs/(ptp*ap)=',f9.6)        
606		Format(//37X,'pressure force on shroud,fshd/(ptp*ap)=',f10.6)        
608		Format(//37X,'pressure force on body,fbdy/(ptp*ap)=',f9.6)        
609		Format(//37X,'skin friction drag,fdrag/(ptp*ap)=',f9.6)       
610		Format(//37X,'total stream thrust,ft/(ptp*ap)=',f9.6)        
612		Format(//37X,'gross stream thrust,fgross/(ptp*ap)=',f9.6)        
618		Format('1',//49X,'ejector throst characteristics'/49X,'******************************')        
620		Format(//26X,'ptp/p0',8X,'fgross',10X,'fip',11X,'fis',11X,'cvp',11X,'cv')       
622		Format(19X,6F14.5)
			
	End Subroutine perf
	