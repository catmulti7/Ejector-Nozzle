      Subroutine check(j, shock)
      !
      !     Check for coalescence of characteristics
      !			��������ĺϲ�
	      Common x(2, 100), y(2, 100), p(2, 100), t(2, 100) !����100�о���
	      Common xslp(100), yslp(100), amp(100), theta(100), php(100), ams(100), phs(100), asass(100), dasdx(100), islp
	      Common xis(21, 26), yis(21, 26), w(21), tau(26), nsonic, nangle
  
	      Common xsonic(26), ysonic(26), psonic(26), tsonic(26), isonic
	      Common xcone(100), ycone(100), pcone(100), tcone(100), icone
  
	      Common xshd(100), yshd(100), dysdx(100), nshd
	      Common xbdy(100), ybdy(100), dybdx(100), nbdy
  
              Common wtfl, hshp, tos, top, gams, gamp, fung, amr, angr, apref, assaps, xprim, yprim, aprim, dprim, cfl, cvl, cona,&
              &dbdy, dshd, edn, pamb, yratio, pi, conva, convr, fdim, ndata, nstop, solve, choke, change, charge, typ, point, stag
  
	      Common pts(25), area(25), wleak(25), title(18), niter, try
				
        ave(x1,x2)=(x1+x2)/2.0
                
        write(*,*)"call check"
        shock=0.0
        Lref=j-1
        Do L=1, Lref			
          k=j-L
          If(x(2,j)>x(2,2)) Goto 6
          If(y(2,j)>y(2,2)) Goto 6
          x(2,j)=x(2,2)
          y(2,j)=y(2,2)
          p(2,j)=p(2,2)
          t(2,j)=t(2,2)
          Goto 8
6		  If(x(2,j)>x(2,k)) Goto 12
          If(y(2,j)>y(2,k)) Goto 12
8		  Do n=k, Lref
            x(2,n)=x(2,j)
            y(2,n)=y(2,j)
            p(2,n)=p(2,j)
            t(2,n)=t(2,j)
		  End Do
12		End Do
14		If(x(2,j)>x(1,j)) Goto 20
        If(y(2,j)<y(1,j)) Goto 20
        Do i=j, 100
          If(p(1,i)==0.0) Goto 18
          x(2,i)=x(1,i)
          y(2,i)=y(1,i)
          p(2,i)=p(1,i)
          t(2,i)=t(1,i)
		End Do
18		j=i-1
        shock=1.0
20		write(*,*)"end check"
        Return
600		Format(//35X,'Coalesence has occurred at xp='F8.5,',yp=',f8.5)
      End Subroutine check