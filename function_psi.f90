      Function psi(g, am1, del)
      !
      !     Series expansion for pressure ratio across an oblique shock
	    !     斜激波多次膨胀的压比
      !
        exp = 7.0/2.0
        am2 = am1*am1
        am4 = am2*am2
        am6 = am2*am4
        am8 = am4*am4
        beta = am2 - 1.0
        a = 1.0
        b = g*am2/sqrt(beta)
        c = g*am2/(4.0*beta**2)*((g+1.0)*am4-4.0*beta)
        d = g*am2/beta**exp*((g+1.0)**2/32.0*am8-(7.0+12.0*g-3.0*g*g)/24.0*am6+3.0/4.0*(g+1.0)*am4-am2+2.0/3.0)
        psi = a + b*del + c*del**2 + d*del**3
        Return
      End Function psi