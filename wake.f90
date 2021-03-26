      Subroutine wake(delw, dela)

      !      Calculation of viscous mixing region
		!		ճ�Ի�����ļ���
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


        funa(g, am) = ((g+1.0)/2.0)**(-(g+1.0)/(2.0*(g-1.0)))*1.0/am*(1.0+(g-1.0)/2.0*am*am)**((g+1.0)/(2.0*(g-1.0)))
        funq(g, am) = sqrt((g+1.0)/2.0*am*am/(1.0+(g-1.0)/2.0*am*am))			
        crocco(g, am) = am*am/(2.0/(g-1.0)+am*am)
        funcp(ca, phi) = (1.0-ca)*phi/(1.0-ca*phi*phi)
        ave(x1, x2) = (x1+x2)/2.0
				
        phib = funq(gams, ams(islp))/funq(gamp, amp(islp))*sqrt(tos/top)
        gami = (gamp+wtfl*gams)/(1.0+wtfl)
        siisi = (1.0+phib)/(1.0-phib)
        caiisq = crocco(gamp,amp(islp))                           
        caisq = caiisq*(1.0-phib)**2/(caiisq*(1.0-phib)**2+(1.0-caiisq))
        sigi = 12.0+2.758*sqrt(caisq)/(sqrt((1.0-caisq)*(gami-1.0)/2.0))
        sigii = siisi*sigi
        aratio=funa(gamp,amr)/funa(gamp,amp(islp))
        Call mix(tos/top,caiisq,phib,0.0,sigvb)                   
        dela = 2.0/sigii*xslp(islp)*yslp(islp)*sigvb
        delw = 2.0/sigii*xslp(islp)*yslp(islp)*aratio*funcp(caisq,phib)*sigvb
        Return
      End Subroutine wake
