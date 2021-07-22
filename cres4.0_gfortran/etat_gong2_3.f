 
c******************************************************************
 
	subroutine etat_gong2_3(p,t,xchim,deriv,
     1	ro,drop,drot,drox,drott,drotp,drotx,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
 
c	equation d'etat pour GONG etape 2
c	d'apres la note de J. Christensen-Dalsgard 1 / 3 / 88
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 3
 
c entree :
c	p : pression
c	t : temperature
c	xchim : composition chimique en grammes
c	deriv=.true. : calcul des derivees
 
c sortie :
c	ro : densite et derivees
c	u : energie interne et derivees
c	nh1, nhe1, nhe2 : taux d'ionisation
c	lamb: degenerescence
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'ctephy.common'
	include 'evol_chim_3.common'
 
	integer ntour,kder
 
	real*8 p,t,xchim(1),x,ro,drop,drot,u,dup,dut,granr,zisai,zai,unsai,
     1	zizaiz,mum1,corr,stor,stor0,dx,unpdx,dstor,ro0,u0,
     2	nh1,nhe1,nhe2,y,lamb,
     3	cte1,cte2,cte15,cte16,cte17,cte18,cte14,cte19,cte20,
     4	cte21,cte22,cte23,cte24,cte25,cte26,
     5	drott,drotp,dutt,dutp,drotx,dutx,
     6	dmunx,ah2ahe,dmux,dgx,dnex,dpx,drox0,dux0,fifi
 
	real*8 aa(2,2),bb(2),tpre,ahmu,ahemu,zih,zihz,cf,a03,
     1	ysahe,kihkt,kihekt,kihe1kt,kt,kihk,kihek,
     3	kihe1k,psi,dnh1ne,dnhene,muem1,bid,t32,drox,dux,
     4	deltp,deltpne,deltpro,munm1,ne00,zbar,zbar2,zbar3,dmu,
     5	xsah,dmuem1ne,dnero,dune,
     6	kihah,kihehe,kihe1he,kihx,kihey,kihe1y,ropre,nepre,
     7	dgne,dgro,dgt,dnet,dpro,dpt,dpne,duro,drop0,drot0,dup0,dut0
 
	real*8	ne0,nel,phi1,phi2,phi3,nemne,dphi1t,dphi2t,dphi3t,npart
 
	logical init,deriv
	data init/.false./
 
	if(.not. init)then	!initialisations
	 init=.true.	!Z proportion d'elements lourds doit avoir ete
			!initialise par appel fictif a opa
	 dx=1.d-8
	 unpdx=1.d0+dx
 
	 granr=kbol/amu
	 ahmu=ah*amu
	 ahemu=ahe4*amu
 
	 zisai=.5d0	![ Zi / Ai ]h		page 9
	 zizaiz=zisai*z0/amu	![ Zi /Ai ]h Z /mu
	 zih=8.d0			!< Zi >h	page 9
	 zihz=zih*z0
	 cf=15.d0			!Cf		  "
 
	 unsai=.0625d0	![ 1 / Ai ]h
	 zai=unsai*z0/amu	![ Zi / Ai ]h / mu
 
	 a03=(hpl/4.d0/pi/pi/me*hpl/echarg/echarg)**3		!2.14
 
	 kihk=kih*eve/kbol
	 kihek=kihe*eve/kbol
	 kihe1k=kihe1*eve/kbol
	 kihah=kih*eve/ah/amu		!pour 2.9
	 kihehe=kihe*eve/ahe4/amu
	 kihe1he=(kihe+kihe1)*eve/ahe4/amu
 
	 tpre=0.d0		!initialisation de "t ro nel appel precedent"
	 ropre=0.d0
	 nepre=1.d0
 
	 cte1=3.d0/2.d0*kbol
	 cte2=sqrt(hpl/me*hpl/kbol/2.d0/pi)**3/exp(1.d0)**2/2.d0 !f/4=cte2*nel/T p.8
	 cte14=cf*a03*20.d0*kih*eve
	 cte15=cf*a03			!2.13
	 cte16=20.d0*kih*eve		!2.13
	 cte17=10.d0*cf*kih*eve*a03	!2.17
	 cte18=cf/2.d0*a03			!2.16
 
	 dmunx=1.d0/ah-1.d0/ahe4
	 ah2ahe=(1.d0/ah-2.d0/ahe4)/amu
 
	 write(2,1)
	 write(6,1)
1	 format(//,1x,'equation d''etat de GONG2 : equation de Saha evitant ',
     1	'la recombinaison, ',/,5x,
     2	'pas de pression de radiation, pas de degenerescence',//)
	endif
 
c	print*,nbelem
c	write(6,2000)(xchim(kder),kder=1,nbelem)
 
	x=xchim(1)
	y=1.d0-x-z0
	xsah=x/ahmu
	ysahe=y/ahemu
 
	kihx=kihah*x		!pour 2.9
	kihey=kihehe*y
	kihe1y=kihe1he*y
c	write(6,*)'kihx,kihey,kihe1y'
c	write(6,2000)kihx,kihey,kihe1y
 
	munm1=xsah+ysahe+zai			!2.7/amu
	ne00=xsah+2.d0*ysahe+zizaiz		!2.18/amu
	zbar=x+2.d0*y+zihz
	zbar2=zbar**2
	zbar3=zbar*zbar2
	cte19=cte15/zbar3
	cte20=cte16*zbar2
	cte21=-cte14/zbar
	cte22=cte19/2.d0
	cte26=cte20*cte19
	cte23=cte26/zbar
	cte24=cte23*2.d0
	cte25=cte21/2.d0
 
c	solution de l'equation de SAHA par la methode de Mihalas
c	Stellar Atmosphere chapitre 5 formules 5-17 et 5-20 et Newton Raphson
c	2 variables : ro et nel
 
c	write(6,*)'p,t'
c	write(6,2000)p,t
 
	if(abs(t-tpre)/t .gt. .2)then	!initialisation de ro et nel si T a
	 if(t .ge. 1.d5)then	!varie de plus de 20% depuis l'appel precedent
	  ro=p/kbol/t/(ne00+munm1)
	  nel=ro*ne00			!ionisation totale
	 else
	  ro=p/kbol/t/munm1		!ionisation partielle
	  nel=ro*zizaiz
	 endif
	else
	 ro=ropre
	 nel=nepre
	endif
c	write(6,*)'ro,nel provisoires'
c	write(6,2000)ro,nel,zizaiz
 
19	kder=1	!indice pour derivees 1(p,t), 2dt
20	kihkt=kihk/t		!ki / k /T
	kihekt=kihek/t
	kihe1kt=kihe1k/t
	kt=kbol*t
	t32=sqrt(t)**3
	bid=cte19*(kt+cte20)
 
	ntour=0
c	write(6,*)'kder / ro,nel',kder
c	write(6,2000)ro,nel
21	ntour=ntour+1
	if(ntour .gt. 30)then
	 write(6,*)'pas de conv. dans saha pour etat_gong2/P,T,epsi,nel,Tp'
	 write(6,2000)p,t,corr,nel,tpre
	 write(6,*)'x,y,z0'
	 write(6,2000)x,y,z0
	 x=-10.d2
	 y=log(x)	!ln(-10) pour sortir avec la trace
	endif
	
	if(nel .le. 0.d0)then
	 write(6,*)'ntour,nel,p,t,xchim(1), ro ini, nel ini',ntour
	 write(6,2000)nel,p,t,xchim(1),p/kbol/t/munm1,ro*zizaiz
	 stop
	endif
 
	psi=2.d0+log(cte2*nel/t32)		!2.19
	dmu=cte19*(kt+cte20)*nel/kt	!2.13
	phi1=psi+kihkt-dmu
	phi1=exp(phi1)*2.d0	! 1 / 2.22
	nh1=1.d0/(1.d0+phi1)
	phi2=psi+kihekt-dmu
	phi2=exp(phi2)/2.d0	! 1 / 2.23
	phi3=psi+kihe1kt-dmu
	phi3=exp(phi3)*2.d0	! 1 / 2.24
	nhe2=1.d0/(phi3*(phi2+1.d0)+1.d0)
	nhe1=phi3*nhe2
 
	dnh1ne=-nh1**2*phi1
	dnhene=-nhe2**2*phi3*(phi2*(phi3+4.d0)+1.d0)
	muem1=nh1*xsah+(nhe1+2.d0*nhe2)*ysahe+zizaiz	!a 1/amu pres
	mum1=munm1+muem1
	dmuem1ne=(xsah*dnh1ne+dnhene*ysahe)*(1.d0-dmu)/nel
 
	ne0=ne00*ro
	nemne=(ne0-nel)*(ne0+nel)
	deltp=bid*nemne/2.d0
	deltpne=-bid*nel
	deltpro= bid*ne00*ne0
 
	bb(1)=kt*ro*mum1+deltp-p	!2.4  1/amu dans mum1
	bb(2)=ro*muem1-nel		!2.8
 
	aa(1,1)=kt*ro*dmuem1ne+deltpne	!jacobien	 d 1 /dne
	aa(1,2)=kt*mum1+deltpro			!d 1 /dro
	aa(2,1)=ro*dmuem1ne-1.d0			!d 2 /dne
	aa(2,2)=muem1					!d 2 /dro
 
	call simq(aa,bb,2)		!solution
 
	corr=1.
27	if(corr*bb(1) .gt. .6*nel .or. corr*bb(2) .gt. .6*ro)then
	 corr=corr/2.d0
	 goto 27
	endif
	nel=nel-bb(1)*corr
	ro=ro-bb(2)*corr
	corr=max(abs(bb(1))/dble(nel),abs(bb(2)/ro))
 
c	write(6,*)'ntour kder/ ro,nel,bb(1),bb(2),corr',ntour,kder
c	write(6,2000)ro,nel,bb(1),bb(2),corr
	if(corr .gt. 1.d-13)goto 21
 
c	write(6,*)  'nel,ne0,ne0-nel,ro,nh1,nhe1,nhe2'
c	write(6,2000)nel,ne0,ne0-nel,ro,nh1,nhe1,nhe2
 
c	calcul de delta, cp, gradad
 
	dmux=3.d0*dmu/zbar-cte24*nel/kt
	dgne=dmuem1ne*ro
	dgro=muem1
	dphi1t=-phi1*(cte21*nel/kt+1.5d0+kihkt)/t
	dphi2t=-phi2*(cte21*nel/kt+1.5d0+kihekt)/t
	dphi3t=-phi3*(cte21*nel/kt+1.5d0+kihe1kt)/t
	dgt=-ro*(nh1**2*dphi1t*xsah+nhe2**2*(dphi2t*phi3*(phi3+2.d0)+
     1	dphi3t*(2.d0*phi2+1.d0))*ysahe)
	dgx=ro*(nh1/ahmu-(nhe1+2.d0*nhe2)/ahemu+dmux*(nh1**2*phi1*xsah+
     1	nhe2**2*phi3*(phi2*(phi3+4.d0)+1.)*ysahe))
	dnero=dgro/(1.d0-dgne)
	dnet =dgt /(1.d0-dgne)
	dnex =dgx /(1.d0-dgne)
 
	npart=munm1*ro+nel	!nombre de particules libres
	dpro=kt*munm1+deltpro
	dpt=kbol*npart+cte22*nemne*kbol
	dpne=kt-bid*nel
	dpx=granr*ro*t*dmunx+bid*ne0*ah2ahe*ro+3.d0*deltp/zbar-nemne*cte23
 
	dpro=dpro+dpne*dnero
	dpt= dpt +dpne*dnet
	dpx= dpx +dpne*dnex
 
	drop=1.d0/dpro
	drot=-drop*dpt
	drox=-drop*dpx
 
c        10        20        30        40        50        60        70
	fifi=nh1**2*phi1*kihx+
     1 nhe2**2*(phi3*(phi3*phi2-1.d0)*kihey+phi3*(2.d0*phi2+1.d0)
     2 *kihe1y)
	u=(1.5d0*kt*npart-cte25*nemne)/ro+kihx*nh1+nhe1*kihey+nhe2
     1 *kihe1y
	duro=-(1.5d0*kt*nel+cte25*(ne0**2+nel**2))/ro/ro
	dut=cte1*npart/ro-(nh1**2*dphi1t*kihx+nhe2**2*((dphi2t*phi3**2-
     1	dphi3t)*kihey+(dphi2t*phi3+dphi3t*(phi2+1.d0))*kihe1y))
	dune=(1.5d0*kt+cte21*nel)/ro-(1.d0-dmu)/nel*fifi
	dux=1.5d0*granr*t*dmunx-cte25/zbar*nemne/ro+cte26*ah2ahe*ne0+
     1	nh1*kihah-nhe1*kihehe-nhe2*kihe1he+dmux*fifi
 
c	write(6,*)'drox,dux,dune,dnex'
c	write(6,2000)drox,dux,dune,dnex
 
	duro=duro+dune*dnero
	dut= dut +dune*dnet
	dux= dux +dune*dnex
 
	dup=     duro*drop
	dut=dut +duro*drot
	dux=dux +duro*drox
 
	if(.not.deriv)return
 
	goto(30,31),kder
30	tpre=t
	ropre=ro
	nepre=nel
 
	ro0=ro
	drop0=drop
	drot0=drot
	drox0=drox
	u0=u
	dup0=dup
	dut0=dut
	dux0=dux
 
	stor0=t
	stor=stor0*unpdx
	dstor=stor-stor0
	t=stor
	kder=2
	goto20
 
31	t=stor0		!derivee/t
	drott=(drot-drot0)/dstor
	drotp=(drop-drop0)/dstor
	drotx=(drox-drox0)/dstor
	dutt =(dut -dut0 )/dstor
	dutp =(dup -dup0 )/dstor
	dutx =(dux -dux0 )/dstor
 
c	write(6,*)'d.,d.-d.0'
c	write(6,2000)drot,drot-drot0,drop,drop-drop0,dut,dut-dut0,
c	1	dup,dup-dup0
 
	ro=ro0
	drop=drop0
	drot=drot0
	drox=drox0
	u=u0
	dup=dup0
	dut=dut0
	dux=dux0
 
	lamb=5		!degenerescence non calculee
 
2000	format((1x,1p8d10.3))
 
	return
 
	end
