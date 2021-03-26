      Function pmer(amp, angp, amq, gam)
      !
      !     Prandtl meyer expansion angle 普朗特-迈耶膨胀角
      !
        ak = sqrt((gam-1.0)/(gam+1.0))
        famp = sqrt(amp*amp-1.0)
        thetap = atan(ak*famp)/ak - atan(famp)
        famq = sqrt(amq*amq-1.0)
        thetaq = atan(ak*famq)/ak - atan(famq)
        pmer = angp + (thetaq-thetap)
        Return
      End Function pmer