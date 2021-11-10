! 	超声速引射喷管，等熵流场计算(特征线法)
!	作者：曹飞飞
!	版本：Ejector C-D Nozzle V1.1
!	日期：2021年1月19日
!	描述：
!		//用于详细说明此程序文件完成的主要功能，与其他模块或函数的接口，
!		//输出值、取值范围、含义及参数间的控制、顺序、独立或依赖等关系
!	其他：
!		//其他内容的说明
!	修改历史记录列表：
!		//修改日期：修改者：修改内容描述：
!		//主要函数列表，每条修改记录应包含修改日期、修改者及修改内容简述
	
	
	  Program  CADCP!(input,output,tape5=input,tape6=output,tape7=output)
!     Eject nozzle computer program 引射喷管计算程序
!     input parmeters 输入参数

!     wtfl=ratio of secondary to primary corrected weight flow 二次流与主流换算质量流量比
!     hshp=ratio of secondary to primary total pressure 次主流总压比
!     tos=total temperature of secondary flow K(R°)二次流总温
!	  gams=ratio of cp/cv of secondary flow 次流比热比
!	  top=total temperature of primary flow 主流总温
!     gamp=ratio of cp/cv of primary flow 主流比热比
!     amr=initial primary mach number 主流流场计算的初始马赫数
!	  对于塞式喷管amr=1.003；对于其他喷管amr=1.001,format(6e12.0)
!     angr=angle of primary lip(degrees) 主喷管唇口角度
      !对于平面声速线angr=0；对于其他喷管angr<0
!     yratio=primary nozzle radius ratio 主喷管半径比
!     xprim=primry flow station !主喷嘴出口相对于护罩轮廓点坐标系统的位置，厘米（英寸）
!     dprim=diamieter of primary nozzle exit 主喷管出口直径
!     dshd=diameter of shroud 外罩直径，将NSHD设置为1时，其值为DSHD直径值,厘米（英寸）的圆柱形导流罩喷射器的性能
!     dbdy=diameter of body at primary flow station 塞式喷管主喷嘴出口直径cm(in)
!     cona=cone centerbody angle (degrees)塞锥半角
!     end=eject length 喷管出口相对轴向位置

!     Command parameters

!     set fdim=0.0 for two dimensional flow 对于二维平面流动，此值设为0.0
!     set fdim=1.0 for axisymmetric flow 对于二维旋转轴对称流动，此值设为1.0
!     set solve=0.0 for non-mixing solution 对于非混合解设置，solve=0.0
!     set solve=1.0 for mixing solution 对于计算混合解设置，solve=1.0 for 换算流量wtfl＞0.04
!     set solve=2.0 for impingement solution 对于计算冲击碰撞解设置，solve=2.for 0≤换算流量wtfl<0.04
!     set print=0.0 for no prtintout of primary flow field 此print在程序代码中改为prt代替;print=0.0 for 不输出主流流场设置
!     set print=1.0 for print-out of final primary flow field 输出主流最终流场;设置print=1.0 for 要输出主流流场
!     set print=2.0 for print-out of every iteration primary flow field 输出主流每步迭代的流场；
!     set iplot=0 for no calcomp plot 不输出主流流场
!     set iplot=1 for calcomp plot on every iteration solution
!     set iplot=2 for calcomp plot on final solution

!					Card 6 program constants and options 程序约束和选项
!		set ndata=21 沿声速线的离散点数,ndata=21 for angr<0
!		set ndata=8 沿声速线的离散点数,ndata=21 for angr=0
!		nshd 从输入文件读进来的外罩轮廓点数
!		nbdy 从输入文件读进来的中心体轮廓点数
!		niter 重新启动选项之前完成的迭代次数
!		iplot parameter which controls plotting routines
!		ipnch set ipnch=0
!		iprnt set iprnt=0
!		icomp set icomp=0
	
	!	Card 7 niter>0	
!	Pts secondary total pressure ratio (hshp) at the I th iteration ,where I=1,niter
!	area(I)		minimum computed secondary flow area ratio As/As* at I th iteration , for solve=0.0 ,1.0, where I=1,niter
!	wleak(I)	computed leakage secondary weight flow ratio at the I th iteration for solve =2.0, where I=1, niter
!	REJECT中提供了一个重启选项，用于在卡7上读取的值之间迭代求解。
!	每次迭代后，打印出泵送特性的PTS(1)、ARea(1)和WLEAK(1)的值。
	
	!	Card 8 shroud geometry (for NSHD>1)	外罩几何
!	xshd(i) axial location of shroud coordiate point, cm(in.)
!	yshd(i) radial location of shroud coordiate point, where I=1 ,NSHD ,cm(in.)	
!	xbdy(i) axial location of plug coordiate point, cm(in.)
!	ybdy(i) radial location of plug coordiate point, cm(in.)

!	    Description of output	输出列表描述 Paper 31-32
!	xslp jet boundary axial position
!	yslp radial location of jet boundary
!	amp primary stream mach number
!	theta flow angle of jet boundary 流向角
!	P/Ptp ratio of static to primary total pressure
!	ams secondary stream mach number
!	P/Pts ratio of static to secondary total pressure
!	as/as* ratio of secondary flow area to secondary critical area
	
		!	other 第五页
	!	Reyprm	主流雷诺数
	!	Reysec	二次流雷诺数
	!	delshd	位移厚度
	!	thetas	动量厚度
	!	cfshd	壁面摩擦阻力系数
	
	!	第六页 音速线求解计算
	!	xsnoic 
	!	ysnoic 
	!	tsnoic 
	
	!	第七页 重要的性能参数
!	喷管性能重要性能参数
	
		!	第八页 不同主流落压比下引射喷管性能参数	
!		FGROSS gross stream thrust，(F-PeAe)/(PpAp) 总推力
!   FIP ideal thrust of primary nozzle 主流理想推力
!   FIS ideal thrust of secondary stream 二次流理想推力
!   CVP gross thrust coefficient 总推力系数
!		CV nozzle efficiency 喷管效率
	
!	    Appendix D sample output listing
!	primary nozzle exit diameter;			dprim=2.00000   主喷管出口直径
!	primary nozzle exit mach NO.;			amr=1.00100     主喷管出口马赫数
!	primary nozzle lip angle;					angr=-16.00000  主喷管收敛角
!	primary nozzle radius ratio;			yratio=1.40000  主喷管收缩比：主喷管7截面直径/8截面直径
!   location of primary nozzle;     xprim=0.95860   外罩最小截面x周向位置
!   ejector length measured from primary nozzle;	end=3.19592    引射喷管长度终点位置****后文用edn变量表示
!   primary nozzle reynolds number;								reyprm=4.000E+06 主喷管雷诺数
!	total temperature of secondary flow;						Tts=560.000     次流总温
!	ratio of specific heats for secondary flow;			gams=1.40000    次流比热比
!	total temperature of primary flow;							Ttp=560.000     主流总温
!	ratio of specific heats for primary flow;				gamp=1.40000    主流比热比
!	flow in ejector is axisymetric;									flow=1.0        轴对称流动
!	number of points specifing shroud contour;			nshd=86     指定引射外罩离散点数
!	number of points specifing plug contour;				nbdy=0      指定引射塞式离散点数
!	number of points specifing enterence conditions; ndata=21   指定入口离散点数


!	全局变量区
	! x,y 坐标 p,t 对应坐标上的温度与压强
	! 最终计算结果存放于此，是流场中离散点的坐标与该点的参数
	! 两行用于进行迭代计算
	  Common x(2, 100), y(2, 100), p(2, 100), t(2, 100) !两行100列矩阵


	  Common xslp(100), yslp(100), amp(100), theta(100), php(100), ams(100), phs(100), asass(100), dasdx(100), islp
	  Common xis(21, 26), yis(21, 26), w(21), tau(26), nsonic, nangle
  
	  Common xsonic(26), ysonic(26), psonic(26), tsonic(26), isonic


	  Common xcone(100), ycone(100), pcone(100), tcone(100), icone
		
	! 外罩坐标与斜率，nshd=总的外罩轮廓点数目
	  Common xshd(100), yshd(100), dysdx(100), nshd

	
	  Common xbdy(100), ybdy(100), dybdx(100), nbdy
	
	! 命名规律：s=secondary 次流
	! 命名规律：p=primary 主流
	! 标？意为不确定
	!     wtfl=ratio of secondary to primary corrected weight flow 二次流与主流换算质量流量比
	!     hshp=ratio of secondary to primary total pressure 次主流总压比
	!     tos=total temperature of secondary flow K(R°)二次流总温
	!	  top=total temperature of primary flow 主流总温
	!	  gams=ratio of cp/cv of secondary flow 次流比热比
	!     gamp=ratio of cp/cv of primary flow 主流比热比
	!     fung 公式（33）中参数
	!     amr=initial primary mach number 主流流场计算的初始马赫数
	!	  对于塞式喷管amr=1.003；对于其他喷管amr=1.001,format(6e12.0)
	!     angr=angle of primary lip(degrees) 主喷管唇口角度
			  !对于平面声速线angr=0；对于其他喷管angr<0
	!	  ？ apref 主流参考面积
	!	 assaps=as*/ap*(s=star,即为*)
	!     xprim=primry flow station !主喷嘴出口相对于护罩轮廓点坐标系统的位置，厘米（英寸）
	!     dprim=diamieter of primary nozzle exit 主喷管出口直径
	!	  yprim和aprim的定义不清楚
	!     dshd=diameter of shroud 外罩直径，将NSHD设置为1时，其值为DSHD直径值,厘米（英寸）的圆柱形导流罩喷射器的性能
	!     dbdy=diameter of body at primary flow station 塞式喷管主喷嘴出口直径cm(in)
	!     cona=cone centerbody angle (degrees)塞锥半角
	!     end(edn)=eject length 喷管出口相对轴向位置
	!     ？pamb 某压强
	!     yratio=primary nozzle radius ratio 主喷管半径比
	!	  conva, convr 角度转弧度和弧度转角度参数
	!	slove 控制参数
	  !     set solve=0.0 for non-mixing solution 对于非混合解设置，solve=0.0
	  !     set solve=1.0 for mixing solution 对于计算混合解设置，solve=1.0 for 换算流量wtfl＞0.04
	  !     set solve=2.0 for impingement solution 对于计算冲击碰撞解设置，solve=2.for 0≤换算流量wtfl<0.04
	!	choke=-1/0/1 表示次流是否拥塞（详见estmp程序)
	!	？change 子程序flow中使用，类似一个标志量
	!	？charge 子程序slid中使用，类似一个标志量
	  Common wtfl, hshp, tos, top, gams, gamp, fung, amr, angr, apref, assaps, xprim, yprim, aprim, dprim, cfl, cvl, cona, dbdy&
	  &, dshd,edn, pamb, yratio, pi, conva, convr, fdim, ndata, nstop, solve, choke, change, charge, typ, point, stag
  
	  Common pts(25), area(25), wleak(25), title(18), niter, try
  
	  Common /bnlyr/xsum(100), ysum(100), delshd(100), thetas(100), cfshd(100), xcne(100), ycne(100), delcne(100), thetac(100),&
	  &cfcne(100), reyprm, pop, aop, vop, reysec, pos, aos, vos, pex, xscale
	  Common /prflr/xprf(6), yprf(6, 25), qprf(6, 25), nprf(6), max, nmax
	  Common /cplot/iplot, xstart, ystart, xspan, yspan, scale, span, axis, xorgn, yorgn, xshft, yshft, kkk(14), pp(14), xdown(100),&
	  &yacros(100)
	  Common /outpn/ipnch, iprnt, icomp
  
	  Real k1, k2, k3, k4, mach
!
!     function statements
!
	  funm(g, ph) = sqrt(2.0/(g-1.0)*(ph**(-(g-1.0)/g)-1.0))
	  funp(g, am) = (1.0+(g-1.0)/2.0*am*am)**(-g/(g-1.0))
	  funq(g, am) = sqrt((g-1.0)/2.0*am*am/(1.0+(g-1.0)/2.0*am*am))
	  funw(g, vel) = sqrt(2.0/(g-1.0)*vel*vel/(1.0-vel*vel))

!   输入
CHARACTER(len=80):: paramfilename ='../params.txt'
CHARACTER(len=80):: resultfilename ='../result.txt'
open(unit=1, file=paramfilename)
open(unit=7, file=resultfilename)


600 Format (10XE12.0,/,10XE12.0,/,10XE12.0,/,10XE12.0,/,10XE12.0,/,10XE12.0)
602 Format (10XE12.0,/,10XE12.0,/,10XE12.0,/,10XE12.0,/,10XE12.0,/,10XE12.0)
604 Format(10XI6,/,10XI6,/,10XI6,/,10XI6,/,10XI6,/,10XI6,/,10XI6,/,10XI6)
606 Format(10XE12.0,/,10XE12.0,/,10XE12.0)
608 Format(10XE12.0,/,10XE12.0)


! 标题
5 read (1,"(18A4)")(title(i), i=1, 18) !Card 1 Title Card 任何字母数字特征,暂注释掉 

! 初始参数
read(1,600)wtfl, hshp, tos, gams, top, gamp
Write(*,*)"wtfl, hshp, tos, gams, top, gamp ="
Write(*,*)wtfl, hshp, tos, gams, top, gamp

Read (1,600) amr, angr, yratio, xprim, dprim, dshd
Write(*,*)"amr, angr, yratio, xprim, dprim, dshd="
Write(*,*)amr, angr, yratio, xprim, dprim, dshd

Read (1,602) dbdy, cona, edn, reyprm, delshd(1), fdim
Write(*,*)"dbdy, cona, edn, reyprm, delshd(1), fdim="
Write(*,*)dbdy, cona, edn, reyprm, delshd(1), fdim



Read (1,600) k1, k2, k3, k4, solve, prt
Write(*,*)"k1, k2, k3, k4, solve, prt="
Write(*,*)k1, k2, k3, k4, solve, prt

Read (1,604) ndata, nshd, nbdy, niter, iplot,ipnch, iprnt, icomp
Write(*,*)"ndata, nshd, nbdy, niter, iplot, ipnch, iprnt, icomp="
Write(*,*)ndata, nshd, nbdy, niter, iplot, ipnch, iprnt, icomp

If (niter>0)  Read (1, 606)(pts(i), area(i), wleak(i), i=1, niter)
If (niter>0)  write (*,*)"prev data"
If (niter>0)  write (*,*)(pts(i), area(i), wleak(i), i=1, niter)

If (icomp==2) Read (1, 608) xprf(1), xprf(6)
If (icomp==2) write (*,*) xprf(1), xprf(6)

! 几何参数，外罩坐标
If (nshd>1)  Read (1, 608)(xshd(i), yshd(i), i=1, nshd)	  !Card 8 Shoud Geometry (for NSHD>1)
If (nshd>1)  write (*,*)"shoud contour points data"
If (nshd>1)  write (*,*)(xshd(i), yshd(i), i=1, nshd)

! 几何参数，中心锥体坐标
If (nbdy>0)  Read (1, 608)(xbdy(i), ybdy(i), i=1, nbdy)
If (nbdy>0)  write (*,*)"centerbody contour points data"
If (nbdy>0)  write (*,*)(xbdy(i), ybdy(i), i=1, nbdy)

Write(*,*)"file read done"
close(1)

! 输入结束

! !     program eject
! 	  Write(*,*)"title"
! !5	  Read (*,*)(title(i), i=1, 18) !Card 1 Title Card 任何字母数字特征,暂注释掉 
! 	  5	  Read (*,*)(title(i), i=1, 2) !Card 1 Title Card 任何字母数字特征,暂注释掉 
! 		Write(*,*)"card 1 done"
! 	  !If (eof(5)) 7, 9, 7 !********************************
! !7	  Stop
! !9	  Continue
! 	  Read (*,*) wtfl, hshp, tos, gams, top, gamp					!Card 2 Data Variables 
! 	  Write(*,*)"card 2 done"
! 	  Write(*,*) wtfl, hshp, tos, gams, top, gamp
! 	  Read (*,*) amr, angr, yratio, xprim, dprim, dshd			!Card 3 Data Variables
! 	  Write(*,*)"card 3 done"
! 	  Read (*,*) dbdy, cona, edn, reyprm, delshd(1), fdim	!Card 4 Data Variables	****end用edn代替
! 	  Write(*,*)"card 4 done"
! 	  Read (*,*) k1, k2, k3, k4, solve, prt								!Card 5 Diagram Constains and Options	****print用prt代替
! 	  Write(*,*)"card 5 done"
! 	  Read (*,*) ndata, nshd, nbdy, niter, iplot, ipnch, iprnt, icomp		!Card 6 Program Constants and Options
! 	  Write(*,*)"card 6 done"

! 	  If (niter>0)	Read (5, 502)(pts(i), area(i), wleak(i), i=1, niter)		!Card 7 Niter>0
		
! 	  If (iplot>0)	Read (5, 502) xstart, ystart, xspan, yspan, axis, scale	!P26页：本报告画图程序已被删除
! 	  If (icomp==2) Read (5, 506) xprf(1), xprf(6)
		
! 	  If (nshd>1)		Read (5, 502)(xshd(i), yshd(i), i=1, nshd)	!Card 8 Shoud Geometry (for NSHD>1)
! 	  If (nbdy>0)		Read (5, 502)(xbdy(i), ybdy(i), i=1, nbdy)

	  nstop = 6
	  try = 0.0
	  cas = 0.0			!****case用cas代替
	  choke = -1.0
	  If (icomp==0 .And. nbdy>0) icomp = 1
	  If (angr<0.0) try = 1.0
	  Call start
	  If (try==0.0) Call datum
	  If (iplot>0) prt = -1.0
		
10	niter = 1 + niter
	  islp = 1
	  skip = 0.0
	  typ = 1.0		!****type用typ代替
	  stag = 0.0   ! in original program stag=0
	  icone = 1
	  If (wtfl==0.0) stag = -1.0
	  change = 0.0
	  charge = 0.0
	  If (iplot==1) Call plotc
	  If (iplot==1) skip = 1.0
	  If (prt==1.0 .And. choke/=-1.0) skip = 1.0
	  If (prt==2.0) skip = 1.0
	  pts(niter) = hshp
	  wleak(niter) = 0.0
	  Call clear(0,1)
	  assaps = fung*wtfl/hshp		!公式（33）
	  xslp(islp) = xprim
	  yslp(islp) = yprim
	  ams(islp) = 0.200
	  amp(islp) = 1.500
	  Call flow(islp)			!P53页flow()函数
	!   php(1)=0.28564
	!   phs(1)=0.97090
	!   amp(1)=1.46709
	!   ams(1)=0.20582
	!   asass(1)=2.88378
	  area(niter) = asass(islp)
	  If (point==-1.0) Goto 10
		

	  pamb = php(islp)
	!   write(*,*)"pamb=",pamb
	!   write(*,*)"stag=",stag
	  If (try==1.0) Call sonics
	  vel = funq(gamp, amr)
	  write(*,*)"x,y,p,t(1,j) init"
		Do j = 1, 100
			write(*,*)"j=",j
			x(1, j) = xprim
			y(1, j) = yprim
			If (j==1) delv = 0.0
			If (j==2 .Or. j==3) delv = k1*(1.0+k2)**2/3.0
			If (j>=4) delv = k1*(1.0+k2)**(j-2)
			If (delv>k3) delv = k3
			vel = vel+delv
			mach = funw(gamp, vel)
			If (mach>amp(islp)) mach = amp(islp)
			p(1, j) = funp(gamp, mach)
			t(1, j) = pmer(amr, angr, mach, gamp) !后续function pmer(amp,angp,amq,gam)
			If (j==3) vel = funq(gamp, amr)
			If (mach==amp(islp)) Goto 14
		End Do
	  Call exit

14	If (skip==1.0) Call outfld(1)

		Do 22 i = 2, ndata
			x(2, i) = xsonic(i)
			y(2, i) = ysonic(i)
			p(2, i) = psonic(i)
			t(2, i) = tsonic(i)
			If (try==1.0) Call sonpt(i)
			Do j = 2, 100
				! write(*,*)"bp5"
				! write(*,*)"p(1,j)=",p(1,j)
				If (p(1,j)==0.0) Goto 18
				Call field(j)
			End Do
18		islp = 1 + islp
			Call slip(j)
			Call slid(j)
			If (solve==2.0 .And. point==-1.0) Call break(j)
			! write(*,*)"stag=",stag
20		If (stag==2.0) Goto 34
			area(niter) = amin1(area(niter), asass(islp))
			write(*,*)"area(niter)",area(niter)
			Call clear(0, j)
			radius = sqrt((x(1,2)-x(1,1))**2+(y(1,2)-y(1,1))**2)
			nsert = radius/k4
			If (i==ndata .And. nsert>=2) Call insert(nsert-1)
			If (skip==1.0) Call outfld(1)
			If (icomp==2) Call profle
			If (point/=0.0) Goto 34
22 continue

24	lslp = islp + 1
		Do jslp = lslp, 100		!注意此处jslp可能有错
			Call bound(2)
			Do j = 3, 100
				If (p(i,j)==0.0) Goto 28
				Call field(j)
				Call check(j, shock)
				If (shock==1.0) Goto 30
			End Do
28		islp = 1 + islp
			Call slip(j)
			Call slid(j)
			If (solve==2.0 .And. point==-1.0) Call break(j)
30		If (stag==2.0) Goto 34
			area(niter) = amin1(area(niter), asass(islp))
			Call clear(1, j)
			radius = sqrt((x(1,2)-x(1,1))**2+(y(1,2)-y(1,1))**2)
			nsert = radius/k4
			If (nsert>=2) Call insert(nsert-1)
			If (skip==1.0) Call outfld(1)
			If (icomp==2) Call profle
			If (point/=0.0) Goto 34
		End Do

! 34	 write(*,*)"bp6"
	34	If (cas==1.0) Goto 38
		If (solve<=1.0) Call estmp
		If (solve==2.0 .And. typ==0.0) choke = 0.0
		If (solve==2.0 .And. typ==1.0) Call estmw
36	If (icomp>0 .And. choke==0.0) Call comp(k4, skip)
		
		Call outslp
		! write(*,*)"bp7"
		If (iplot==1) Call plotl
		If (iplot==1 .And. choke/=-1.0) Goto 40
		If (prt==0.0 .And. choke/=-1.0) Goto 40
		If (prt==1.0 .And. cas==1.0) Goto 38
		If (prt==2.0 .And. choke/=-1.0) Goto 40
		If (choke/=-1.0) cas = 1.0
		If (cas==1.0 .And. iplot==2) iplot = 1
		Goto 10
38	If (icomp>0 .And. prt==0.0) Call comp(k4, skip)
		If (cas==1.0 .And. iplot==1) Call plotl
40	Goto 5

write(*,*)"program end normally"
		
500	  Format (18a4)
502	  Format (6E12.0)
504	  Format (8I6)
506	  Format (2E12.0)
		End Program CADCP
