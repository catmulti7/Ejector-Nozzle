      Subroutine mix(tob, ca2, phib, phid, sigvb)
	
      !     Calculation of pertinent paraneters for mixing solution
		!			混合解相关参数的计算
	
        Dimension eta(410), phi(410), tor(410), ei1(410), ei2(410), ei3(410)
        pi = sqrt(3.1415927)
        etarb = -5.0
        deeta = 0.020
        eta(1) = etarb
        phi(1) = phib
        tor(1) = tob
        a1 = phib*(1.0-ca2)/(tob-ca2*phib**2)
        b1 = phib*a1
        c1 = tob*a1
        ei1(1) = a1*etarb
        ei2(1) = b1*etarb
        ei3(1) = c1*etarb
        Do 12, i = 1, 400
          eta(i+1) = eta(i) + deeta
          aveta = (eta(i+1)+eta(i))/2.0
          phi(i+1) = phi(i) + (1.0-phib)*exp(-(aveta**2))*deeta/pi
          tor(i+1) = (tob*(1.0-phi(i+1))+phi(i+1)-phib)/(1.0-phib)
          a2 = phi(i+1)*(1.0-ca2)/(tor(i+1)-ca2*phi(i+1)**2)
          b2 = phi(i+1)*a2
          c2 = tor(i+1)*a2
          ei1(i+1) = ei1(i) + (a1+a2)*deeta/2.0
          ei2(i+1) = ei2(i) + (b1+b2)*deeta/2.0
          ei3(i+1) = ei3(i) + (c1+c2)*deeta/2.0
          a1 = a2
          b1 = b2
          c1 = c2
          j = i + 1
          If (eta(i+1)) 12, 12, 10
10		  a = ei1(i+1) - ei2(i+1)
          b = ei1(i) - ei2(i)
          If (abs(a-b)-1.0E-05) 14, 14, 12
12		Continue
14		ei1j = (ei1(j)-ei2(j))/(1.0-phib)
        etam = eta(j) - (ei2(j)-phib*ei1(j))/(1.0-phib)                  
        Do 16, k = 2, 400
          If (ei1(k)-ei1j) 16, 18, 18                               
16		Continue
18		deta = (ei1j-ei1(k-1))/(ei1(k)-ei1(k-1))*deeta
        etaj = eta(k-1) + deta
        phij = phi(k-1) + deta/deeta*(phi(k)-phi(k-1))
        ei2j = ei2(k-1) + deta/deeta*(ei2(k)-ei2(k-1))
        ei3j = ei3(k-1) + deta/deeta*(ei3(k)-ei3(k-1))
        stnsg = (ei2j-phib*ei1j)/(1.0-phib)
        vbsig = ei1j*(tob-ca2*phib**2)/(1.0-ca2) - etam*phib
        If (phib==0.0) Goto 20
        sigvb = vbsig/phib
        Goto 26
20		Do i = 1, 400
          If (phi(i+1)>phid) Goto 24
		End Do
24		delp = phid - phi(i)
		dphi = phi(i+1)-phi(i)
        deta = eta(i+1) - eta(i)
        dei1 = ei1(i+1) - ei1(i)
        etad = eta(i) + deta/dphi*delp
        ei1d = ei1(i) + dei1/dphi*delp
        sigvb = ei1j - ei1d
26		Return
	End Subroutine mix
	