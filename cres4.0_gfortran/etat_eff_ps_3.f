 
c****************************************************************************
 
	subroutine etat_eff_ps_3(p,t,xchim,deriv,
     1       ro,drop,drot,drox,drott,drotp,drotx,
     2       u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
 
c	Si P < (1 + al) *  Prad alors Ps =Prad(1+al)exp(P/(1+al)Prad - 1)**3
c	derivees secondes numeriques dans ce cas car on ne sort
c	pas toutes les derivees secondes
 
c	on utilise etat_eff
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version 3
 
c entree :
c	p : pression
c	t : temperature
c	xchim : composition chimique
c	deriv=.true. : calcul des derivees secondes
 
c sortie :
c	ro : densite et derivees
c	u : energie interne et derivees
c	nh1, nhe1, nhe2 : taux d'ionisation
c	lamb: degenerescence
 
 
	implicit none
 
	include 'ctephy.common'
 
       real*8	p,t,xchim(1),ps,al, unpal,dx,unpdx,s,s1,pr,dsdp,
     + stor,stor0,
     1 ro,drop,drot,drox,drott,drotp,drotx,dsdt,dsdpr,dprdt,ds1ds,
     2 u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb,ds1dpr,dstor,
     3 ro1,drop1,drot1,drox1,drott1,drotp1,drotx1,dpsdpr,dpsdt,
     4 u1,dup1,dut1,dux1,dutt1,dutp1,dutx1,nh11,nhe11,nhe21,lamb1,
     5	dpsdp,ds1dt
 
	logical init,deriv
	data init/.true./
 
2000	format(1x,1p8d10.3)
 
	if(init)then
	 init=.false.
	 al=.5
	 unpal=1.+al
	 dx=1.d-6
	 unpdx=1.+dx
c	 write(6,2000)aradia
	 write(2,1)unpal
	 write(6,1)unpal
1	 format(//,5x,'-------------------------------------------',//,
     1	1x,'equation d''etat EFF modifiee si P < Prad*',1pe10.3,//)
	endif
 
	pr=aradia/3.*t**4
	if(p .ge. unpal*pr)then	
	 call etat_eff_3(p,t,xchim,deriv,
     1       ro,drop,drot,drox,drott,drotp,drotx,
     2       u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
	else
 
	 dprdt=pr/t*4.
	 s=p/unpal/pr-1
	 dsdp=(s+1)/p
	 dsdpr=-(s+1)/pr
	 dsdt=dsdpr*dprdt
	 s1=pr*al*exp(s**3)
	 ds1ds=3.*s**2*s1
	 ds1dpr=s1/pr+ds1ds*dsdpr
	 ds1dt=ds1dpr*dprdt
	 ps=pr+s1
	 dpsdpr=1.+ds1dpr
	 dpsdp=ds1ds*dsdp
	 dpsdt=dpsdpr*dprdt
 
	 call etat_eff_3(ps,t,xchim,.false.,
     1       ro,drop,drot,drox,drott,drotp,drotx,
     2       u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
	 drot=drot+drop*dpsdt
	 drop=drop*dpsdp
	 dut=dut+dup*dpsdt
	 dup=dup*dpsdp
 
	 if(deriv)then		!derivees secondes numeriques
 
c	  write(6,*)'derivees par rapport a T'
	  stor0=t
	  stor=stor0*unpdx
	  if(stor .eq. 0.)stor=dx
	  dstor=stor-stor0
	  t=stor
 
	  pr=aradia/3*t**4
	  dprdt=pr/t*4.
	  s=p/unpal/pr-1
	  dsdp=(s+1)/p
	  dsdpr=-(s+1)/pr
	  dsdt=dsdpr*dprdt
	  s1=pr*al*exp(s**3)
	  ds1ds=3.*s**2*s1
	  ds1dpr=s1/pr+ds1ds*dsdpr
	  ds1dt=ds1dpr*dprdt
	  ps=pr+s1
	  dpsdpr=1.+ds1dpr
	  dpsdp=ds1ds*dsdp
	  dpsdt=dpsdpr*dprdt
 
	  call etat_eff_3(ps,t,xchim,.false.,
     1       ro1,drop1,drot1,drox1,drott1,drotp1,drotx1,
     2       u1,dup1,dut1,dux1,dutt1,dutp1,dutx1,nh11,nhe11,nhe21,lamb1)
	  drot1=drot1+drop1*dpsdt
	  dut1=dut1+dup1*dpsdt
	  drop1=drop1*dpsdp
	  dup1=dup1*dpsdp
 
	  drott=(drot1-drot)/dstor
	  drotp=(drop1-drop)/dstor
	  drotx=(drox1-drox)/dstor
	
	  dutt=(dut1-dut)/dstor
	  dutp=(dup1-dup)/dstor
	  dutx=(dux1-dux)/dstor
 
	  t=stor0
	 endif
	endif
 
	return
 
	end
