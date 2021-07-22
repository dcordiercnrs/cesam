c****************************************************************************
 
	subroutine etat_eff_dc(p,t,xchim,deriv,
     1       ro,drop,drot,drox,drott,drotp,drotx,
     2       u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
 
c	Si P < (1 + al) *  Prad alors Ps =Prad(1+al)exp(P/(1+al)Prad - 1)**3
c	derivees secondes numeriques dans ce cas car on ne sort
c	pas toutes les derivees secondes
 
c	on utilise etat_eff
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	        D. Cordier, DASGAL
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
	 call equation_etat(p,t,xchim,deriv,
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
 
	 call equation_etat(ps,t,xchim,.false.,
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
 
	  call equation_etat(ps,t,xchim,.false.,
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
 
c****************************************************************************
 
	subroutine equation_etat(p,t,xchim,deriv,
     1       ro,drop,drot,drox,drott,drotp,drotx,
     2       u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
c
c       adaptation du logiciel eff ecrit par jorgen christensen_dalsgaard
c       a partir de p et t, il calcule ro et u et toutes les derivees
c       jusqu'au deuxieme ordre
c       annie baglin le 1 12 1990
c
c       modification: suppression de implicit real*8 (a-h,o-z)
c       et declaration explicite des variables.
c       michel juillet 1990.
 
c	ajout de lamb=5 P. Morel decembre 1990
c	suppression des variables inutiles
 
c	ajout d'un appel a etat_gong2 en cas d'echec de convergence dans eqstp
c	P. Morel octobre 91
 
c	lamb=psi, Daniel Cordier, Novembre 1997
c	suppression de la possibilite de faire Z_CTE=.true.
 
c	CESAM version 3
 
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
 
	include	'cesam_3.parametres'
	include 'modele_3.common'
	include 'evol_chim_3.common'
	
	real*8 lamb
	real*8 ro,drop,drot,drox,drott,drotp,drotx
	real*8 dropl,dropn,drotn,droxn,drofn
	real*8 f,dfpn,dftn,dfxn,dfpl,dftl,dfxl,dptn
	real*8 droffn,droftn,drottn,drotxn,drofxn,drotpn,droxxn
	real*8 dpffn,dptxn,dpfxn,betgam,betbet,dpttn,dpxn
	real*8 xchim(1),p10,t10,p,t
	real*8 nh1,nhe1,nhe2
	real*8 unsro,u,ro2,ro3,psro2,psro3,h
	real*8 dutn,duxn,dupn,dutt,dutx,dutp,dup,dut,dux,
     1       duttn,dutpn,dutxn
	real*8 dhtn,dhxn,dhtpn,dhtxn,dhttn,
     1	dhffn,dhftn,dhfxn,dhpn,dhfn,dpftn
 
	logical ok
	common /marche/ok
 
	logical deriv,init
c	data init/.true./
 
	real*8 droppn
	common/bidon/droppn
	
	real*8 av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     1       ct,crho,cpe,che,cao3,caa,ckh,car,ergev,exhm
	common/consts/av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     1       ct,crho,cpe,che,cao3,caa,ckh,car,ergev,exhm
	
	real*8 xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     1       dlt(4),gm1,tprh,trph,rhxp
	common/eqstd/xii1,ane,rho,ht,pt,cp,dad,
     1       dlt,gm1,tprh,trph,rhxp
	
	integer ihvz
	real*8 anh0
	common/eqscnt/anh0,ihvz
	
	real*8 amm,amm2,amm3
	common/ln10/amm,amm2,amm3
 
        real*8 psi
        common/eqsaux/psi
	
	external bleqst
 
	data init/.true./
 
2000	format(1x,1p8d10.3)
 
	if(init)then
	 call setcns		!initialisation des ctes
	 anh0=.5		!initialisation ???????
	 init=.false.
	 write(2,1)
	 write(6,1)
1	 format(//,5x,'------------------------------------------------',//,
     1	1x,'equation d''etat EFF',//)
	endif
 
	ok=.true.
	p10=log10(p)
	t10=log10(t)
c	print*,ihe4
c	write(6,2000)p10,t10,xchim(1),xchim(ihe4),x0,y0,z0
	if(ihe4 .gt. 1)then
	 if(iz .gt. 1)then	!ihe4-1: indice de He3
c	  write(6,2000)xchim(1),xchim(ihe4),xchim(iz)
c	  pause'eff'	
	  call theffp(p10,t10,xchim(1),1.d0-xchim(1)-xchim(iz),f)
	 elseif(.false.)then     ! Modif. Z_CTE ---> .false. Daniel
	  call theffp(p10,t10,xchim(1),1.d0-xchim(1)-z0,f)	  	
	 else
	  call theffp(p10,t10,xchim(1),xchim(ihe4)+xchim(ihe4-1),f)
	 endif
	else
	 call theffp(p10,t10,xchim(1),1.d0-xchim(1)-z0,f)
        endif
 
	if(ok)then
 
c        sorties
 
c        ro et ses derivees
 
         ro=rho(1)
 
c        derivees premieres de rho
 
	 drofn=rho(2)*rho(1)/f
         dropl=rho(2)/pt(2)
         dropn=(rho(1)/pt(1))*dropl
         drotn=(rho(1)/t)*(rho(3)-dropl*pt(3))
         droxn=rho(1)*(rho(4)-dropl*pt(4))*amm
         drop=dropn
         drot=drotn
         drox=droxn
 
c        derivees secondes de rho et p
 
c        1. derivees secondes naturelles dans le systeme ftx
 
	 dpffn=(pt(1)/(f*f))*((pt(5)/amm)+pt(2)*(pt(2)-1.d0))
         droffn=(rho(1)/(f*f))*((rho(5)/amm)+rho(2)*(rho(2)-1.d0))
         drottn=(rho(1)/(t*t))*((rho(8)/amm)+rho(3)*(rho(3)-1.d0))
	 dpttn=(pt(1)/(t*t))*((pt(8)/amm+pt(3)*(pt(3)-1.d0)))
         droxxn=rho(1)*(rho(10)+rho(4)*rho(4)*amm)*amm
         droftn=(rho(1)/(f*t))*(rho(2)*rho(3)+rho(6)/amm)
         drotxn=(rho(1)/t)*(rho(9)+rho(4)*rho(3)*amm)
	 dptxn=(pt(1)/t)*(pt(9)+pt(4)*pt(3)*amm)
         drofxn=(rho(1)/f)*(rho(7)+rho(4)*rho(2)*amm)
	 dpfxn=(pt(1)/f)*(pt(7)+pt(4)*pt(2)*amm)
         dpftn=(pt(1)/(f*t))*(pt(2)*pt(3)+pt(6)/amm)
 
c        2. derivees premieres naturelles de f dans le systeme ptx
 
         dfpl=1.d0/pt(2)
         dftl=-pt(3)/pt(2)
         dfxl=-pt(4)/pt(2)
         dfpn=(f/pt(1))*dfpl
         dftn=(f/t)*dftl
	 dpxn=pt(4)*pt(1)*amm
	 dfxn=-dfpn*dpxn
         dfxn=f*dfxl
	 dptn=pt(3)*pt(1)/t
 
c        3. derivees secondes naturelles dans le systeme ptx
 
         droppn=dfpn*dfpn*(droffn-dpffn*drofn*dfpn)
         drotpn=dfpn*(droftn-dfpn*dptn*droffn+drofn*dfpn*
     1	(dpffn*dptn*dfpn-dpftn))
	 betbet=-(dpttn+2.d0*dftn*dpftn+dftn*dftn*dpffn)
	 betgam=-(dptxn+dfxn*dpftn+dftn*dpfxn+dftn*dfxn*dpffn)
	 drotxn=drotxn+dftn*drofxn+dfxn*droftn+dftn*dfxn*droffn+betgam*dropn
	 drottn=drottn+2.d0*dftn*droftn+dftn*dftn*droffn+betbet*dropn
         drott=drottn
         drotp=drotpn
         drotx=drotxn
 
c        energie interne  u   ! attention jorgen calcule l'enthalpie
 
         unsro=1.d0/ro
         h=ht(1)
	 u=ht(1)-pt(1)*unsro
 
c        derivees premieres de u
 
         ro2=ro*ro
         psro2=pt(1)/ro2
	 dhfn=ht(2)/(f*amm)
         dhpn=dhfn*dfpn
	 dupn=dhpn-unsro+psro2*dropn
	 dhtn=(ht(3)-ht(2)*pt(3)/pt(2))/(t*amm)
	 dutn=dhtn+psro2*drotn
	 dhxn=ht(4)-ht(2)*pt(4)/pt(2)
	 duxn=dhxn+psro2*droxn
	 dux=duxn
	 dut=dutn
	 dup=dupn	
 
c        derivees secondes de u
 
c	 1- derivees naturelles de h dans le systeme ftx
 
	 dhttn=((ht(8)/amm)-ht(3))/(t*t*amm)
	 dhtxn=ht(9)/(t*amm)
	 dhffn=((ht(5)/amm)-ht(2))/(f*f*amm)
	 dhftn=ht(6)/(f*t*amm2)
	 dhfxn=ht(7)/(f*amm)
 
c	 2- derivees naturelles de h dans le systeme ptx
 
	 dhttn=dhttn+2.d0*dftn*dhftn+dftn*dftn*dhffn+betbet*dhpn
	 dhtpn=dfpn*(dhftn-dfpn*dptn*dhffn+dhfn*dfpn*
     1	(dpffn*dptn*dfpn-dpftn))
	 dhtxn=dhtxn+dftn*dhfxn+dfxn*dhftn+dftn*dfxn*dhffn+betgam*dhpn
 
c	 3- derivees de u dans le systeme ptx
 
         ro3=rho(1)*ro2
	 psro3=2.d0*psro2/rho(1)
	 dutpn=dhtpn+drotn/ro2-psro3*dropn*drotn+psro2*drotpn
	 duttn=dhttn-psro3*drotn*drotn+psro2*drottn
	 dutxn=dhtxn-psro3*drotn*droxn+psro2*drotxn+drotn*dptn/ro2
	 dutt=duttn
	 dutp=dutpn
	 dutx=dutxn
 
c        degres d'ionisation
 
         nh1=xii1(1)
         nhe1=xii1(2)
         nhe2=xii1(3)
 
	 lamb=psi
 
	else	!appel a GONG2 car il y a pb avec EFF pour les cas extremes
	 write(6,*)'appel a ETAT_GONG2 car defaillance de EFF'
	 call etat_gong2_3(p,t,xchim,deriv,
     1       ro,drop,drot,drox,drott,drotp,drotx,
     2       u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
	endif
 
        return
 
        end
 
 
