      Subroutine astar(mach, area, gam)
      !
      !     Calculation of mach number as a function of a/a*
      !
        Real mach				
        funa(g, am) = (2.0/(g+1.0)*(1.0+(g-1.0)/2.0*am*am))**((g+1.0)/(2.0*(g-1.0)))
        funb(g, am) = am*(2.0/(g+1.0)*(1.0+(g-1.0)/2.0*am*am))**((3.0-g)/(2.0*(g-1.0)))
        
        write(*,*)"call astar"
        iter=0
        error=0.00001
10		iter=1+iter
        f=mach-funa(gam,mach)/area
        dfdm=1.0-funb(gam,mach)/area
        mach=mach-f/dfdm
        If(iter<100) Goto 12
        Call exit
        
12		If(abs(f/dfdm/mach)>error) Goto 10
14	write(*,*)"end astar"	
        Return				
	End Subroutine astar
	