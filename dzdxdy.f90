      Subroutine dzdxdy(w, wj, alpha, theta, delw, dw, dt, x, y)
      !
      !     Transformation integration 转型整合
      !
        Complex f, g, h, df, dz, dfdz
        Dimension a(3), b(3), dx(3), dy(3)
        a(1) = w
        a(2) = w - dw/2.0
        a(3) = w - dw
        b(1) = theta
        b(2) = theta - dt/2.0
        b(3) = theta - dt
        Do i = 1, 3
          dz = cmplx(dw, -dt)
          Call dfdqb(a(i), wj, alpha, b(i), delw, dfdz)
          f = cmplx(cos(b(i)),sin(b(i)))
		  g=exp(-a(i))*f
          h = exp(a(i))*conjg(f)
          df = g*dfdz*dz - conjg(h*dfdz*dz/4.0) !conjg()函数:求复数的共轭(即实部不变,虚部取相反数)
          !write(*,*)"dfdz=",dfdz
          dx(i) = real(df)
          dy(i) = aimag(df)   !aimag(df)函数:获取df的虚部
        End Do
        x = (dx(1)+4.0*dx(2)+dx(3))/6.0
        y = (dy(1)+4.0*dy(2)+dy(3))/6.0
        Return
      End Subroutine dzdxdy
