      Subroutine dfdqb(w, wj, alpha, theta, delw, dfdz)
      !
      !     Evaluation of transformation derivative 转换导数的评估
      !
        Complex p, q, r, dfdz
        pi = 3.1415927
        a = pi/alpha*(w-wj)
        b = pi/alpha*theta
        c = pi/alpha*delw
        p = cmplx(sinh(a)*cos(b), -cosh(a)*sin(b))
        q = cmplx(cosh(a)*cos(b)-cosh(c), -sinh(a)*sin(b))
        r = cmplx(cosh(a)*cos(b)-1.0, -sinh(a)*sin(b))
        If (delw==0.0) dfdz = -pi/alpha*p/r
        If (delw/=0.0) dfdz = pi/alpha*(p/q-p/r)
        Return
      End Subroutine dfdqb