      Subroutine spline(x, y, n, slope, dummy)
      !    	Calculation of first and second derivatives
	    !         计算一二阶导数
        Dimension x(100), y(100), s(100), a(100), b(100), c(100), f(100), w(100), sb(100), g(100), em(100), slope(100)
        Do i = 2, n
          s(i) = x(i) - x(i-1)
        End Do
        no = n - 1
        Do i = 2, no
          a(i) = s(i)/6.0
          b(i) = (s(i)+s(i+1))/3.0
          c(i) = s(i+1)/6.0
          f(i) = (y(i+1)-y(i))/s(i+1) - (y(i)-y(i-1))/s(i)
        End Do
        a(n) = -0.5
        b(1) = 1.0
        b(n) = 1.0
        c(1) = -0.5
        f(1) = 0.0
        f(n) = 0.0
        w(1) = b(1)
        sb(1) = c(1)/w(1)
        g(1) = 0.0
        Do i = 2, n
          w(i) = b(i) - a(i)*sb(i-1)
          sb(i) = c(i)/w(i)
          g(i) = (f(i)-a(i)*g(i-1))/w(i)
        End Do
        em(n) = g(n)
        Do i = 2, n
          k = n + 1 - i
          em(k) = g(k) - sb(k)*em(k+1)
        End Do
        slope(1) = -s(2)/6.0*(2.0*em(1)+em(2)) + (y(2)-y(1))/s(2)
        Do i = 2, n
          slope(i) = s(i)/6.0*(2.0*em(i)+em(i-1)) + (y(i)-y(i-1))/s(i)
		End Do       
		Return		
	End Subroutine spline	
	
