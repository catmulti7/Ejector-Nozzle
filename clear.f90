      Subroutine clear(kshift, jref)

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
  
      Common /bnlyr/xsum(100), ysum(100), delshd(100), thetas(100), cfshd(100), xcne(100), ycne(100), delcne(100),&
      & thetac(100), cfcne(100), reyprm, pop, aop, vop, reysec, pos, aos, vos, pex, xscale

      write(*,*)"call clear"
        error=0.0001
        xref=x(2,jref)
        If(islp>1) Goto 14
        Do i=1, 2
          Do j=1,100
            x(i,j)=0.0
            y(i,j)=0.0
            p(i,j)=0.0
            t(i,j)=0.0
          End Do
        End Do
        Do i=1, 100
          If(i>1) delshd(i)=0.0
          thetas(i)=0.0
          cfshd(i)=0.0		!ע����������i����1
          delcne(i)=0.0
          thetac(i)=0.0
          cfcne(i)=0.0
        End Do
        Goto 26
14		j=0
        kmax=100-kshift
        Do k=1, jref
          i=k+kshift
          test=abs(x(2,i)-x(2,i+1))
          If(x(2,i)==x(2,i+1)) Goto 16
          If(kshift>0 .And. test<=error) Goto 16
          j=i+j ! i+j or 1+j????
          x(1,j)=x(2,i)
          y(1,j)=y(2,i)
          p(1,j)=p(2,i)
          t(1,j)=t(2,i)
16		  x(2,i)=0.0
          y(2,i)=0.0
          p(2,i)=0.0
          t(2,i)=0.0
18		End Do
20		kref=1+j
        Do k=kref, 100
          Do i=1, 2
            x(i,k)=0.0
            y(i,k)=0.0
            p(i,k)=0.0
            t(i,k)=0.0
          End Do
        End Do
26		write(*,*)"end clear"
      Return
      End Subroutine clear