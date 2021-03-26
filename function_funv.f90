      Function funv(up, vp, vup, sigsq)
      !
      !     Taylor-maccolli equation
      !
        fun1 = (1.0+vup**2)
        fun2 = (1.0-sigsq)*(up+vp*vup)**2
        fun3 = (1.0-sigsq*(up**2+vp**2))
        funv = (1.0/vp)*(fun1-fun2/fun3)
        Return
      End Function funv
