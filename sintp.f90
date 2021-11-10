      Subroutine sintp(x, y, n, x1, y1)
      !
      !   Interpolation subroutine  ��ֵ�ӳ���
      !
        Dimension x(100), y(100)
        error = 0.0001
        xmin = 100.0
        Do i = 1, n
          delx = x1 - x(i)
          If (abs(delx)<error) delx = 0.0
          If (delx<0.0) Goto 10
          xmin = amin1(delx, xmin)
          If (xmin==delx) k = i        
10		End Do
        If (xmin==0.0) Goto 12
        If (k==1) k = k + 1
        If (k==n) k = k - 1
        delx = x1 - x(k)
        a = y(k)
        d = (y(k)-y(k-1))/(x(k)-x(k-1))
        b = (d+(y(k+1)-y(k))/(x(k+1)-x(k)))/2.0
        c = ((y(k+1)-y(k))-b*(x(k+1)-x(k)))/((x(k+1)-x(k))*(x(k+1)-x(k-1)))
        If (delx < 0.0) b = d
        y1 = a+b*(x1-x(k))+c*(x1-x(k))*(x1-x(k-1))
        Goto 14        
12		y1 = y(k)        
14		Return
	End Subroutine sintp
	