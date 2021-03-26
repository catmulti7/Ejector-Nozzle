      Subroutine runge(up, vp, vup, du, sigsq)
      !
      !     runge kutta integration of taylor-maccoll equation   
		!     taylor-maccoll方程的runge kutta积分  四阶龙格-库塔格式
      !
        Real k1, k2, k3, k4
        k1 = du*funv(up, vp, vup, sigsq)
        k2 = du*funv(up+du/2.0, vp+du/2.0*vup+k1*du/8.0, vup+k1/8.0, sigsq)
        k3 = du*funv(up+du/2.0, vp+du/2.0*vup+k1*du/8.0, vup+k2/8.0, sigsq)
        k4 = du*funv(up+du, vp+du*vup+k3*du/2.0, vup+k3, sigsq)
        up = up + du
        vp = vp + du*(vup+(k1+k2+k3)/6.0)
        vup = vup + (k1+2.0*k2+2.0*k3+k4)/6.0
        Return
      End Subroutine runge
