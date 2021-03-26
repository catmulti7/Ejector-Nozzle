      Subroutine conic(xp, yp, up, vp, vup, gam)
      !
      !     Calculation of conditions in a conical flow field   锥形流场中的条件计算
      !
        delup = 0.0001
        error = 0.0001
        sigsq = (gam-1.0)/(gam+1.0)
        delta = yp/xp
        test = abs(delta) - abs(1.0/vup)
        If (abs(test)<=error*abs(delta)) Goto 12
        signt = test/abs(test)
        delup = signt*delup       

10		Call runge(up, vp, vup, delup, sigsq)
        pest = test
        signp = pest/abs(pest)
        test = abs(delta) - abs(1.0/vup)
        If (abs(test)<=error*abs(delta)) Goto 12
        signt = test/abs(test)
        If (signp==signt) Goto 10
        delup = -0.50*((1.0/delta)+vup)/funv(up, vp, vup, sigsq)
        If (delup==0.0) Goto 12
        Goto 10

12		Return
	End Subroutine conic
	