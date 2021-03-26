      Subroutine coeff(xp, yp, pps, angp, cfs, gam)
      !
      !     Coefficients to the characteristic equations
      !				特征方程的系数
        Dimension cfs(16)

        funm(g, ph) = sqrt(2.0/(g-1.0)*(ph**(-(g-1.0)/g)-1.0))
				
        amp = funm(gam, pps)
        amup = asin(1.0/amp)
        angpu = angp + amup
        angmu = angp - amup
        cfs(1) = tan(angpu)
        cfs(2) = tan(angmu)
        cfs(3) = gam*pps
        cfs(4) = cfs(3)
        cfs(5) = amp*amp
        cfs(6) = cfs(5)
        cfs(7) = tan(amup)
        cfs(8) = cfs(7)
        cfs(9) = sin(angp)
        cfs(10) = cfs(9)
        cfs(11) = sin(amup)
        cfs(12) = cfs(11)
        cfs(13) = yp
        cfs(14) = cfs(13)
        cfs(15) = sin(angpu)
        cfs(16) = sin(angmu)
        Return
      End Subroutine coeff