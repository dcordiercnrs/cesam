c*********************************************************************
 
	subroutine thermo_3(p,t,xchim,m,l,r,deriv,dxchim,mstar,
     1	ro,drop,drot,drox,u,dup,dut,dux,lamb,r_zc,r_ov,lim,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,
     3	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
     4	delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
     5	hp,gradad,gradrad,gradconv,d_grad,w,dw,nh1,nhe1,nhe2,
     6	etat,opa,conv,nuc)
 
c	calcul de la thermodynamique
 
c	Auteur: P. Morel, Departement J.D. Cassini,
c	O.C.A., Observatoire de Nice
c	version: 3
 
c	MODIF:
c	09 10 96 introduction de la gravite effective et de w, dw en entree
c		 dont on tient compte dans Hp
 
c entree :
c	p : pression cgs
c	t : temperature K
c	xchim : composition chimique ATTENTION en 1/mole : *ah pour avoir X
c	xchim(iw)=omega (si iw > 1 ie. rotation non solide)
c	dxchim : d xchim/d m pour critere de Ledoux uniquement
c	m : masse/msol
c	l : luminosite/lsol
c	r : rayon / rsol
c	deriv=.true. : calcul des derivees d'ordre superieur
c	lim: nombre de limites ZR/ZC
c	mstar: masse au temps du calcul, avec perte de masse
c	w, dw : rotation et derivee
 
c sortie : (drop : d ro /dp etc..)
c	ro : densite cgs
c	u : energie interne
c	gradient : gradient pour les calculs
c	epsilon : energie nucleaire cgs
c	kap : opacite Rosseland cm2/gr
c	delta,cp,hp,gradad,gradrad : notation evidentes
c	dcomp : d Xi / dt pour elements chimiques
c	nh1, nhe1, nhe2 : taux d'ionisation
c	lamb: degenerescence
c	d_grad: difference des gradients (rad-ad) pour test de convection
 
c routines externes :
c	etat, opa, conv, nuc
 
c-----------------------------------------------------------------------
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'ctephy.common'
	include 'modele_3.common'
	include 'evol_chim_3.common'
 
	integer i,lim
 
	real*8	p,t,xchim(1),m,l,r,xchimm(pnelem),r_zc(1),r_ov(1),gradconv,
     1	ro,drop,drot,gradient,dgradp,dgradt,dgradl,dgradm,dgradr,
     2	epsilon(5),depsp,depst,u,dup,dut,drox,dux,dkapx,dgradx,
     3	kap,cte8,cte9,cte1,cte13,gradrad,dkapp,dkapt,nh1,nhe1,nhe2,
     4	drott,drotp,drotx,dutt,dutp,dutx,dcomp(1),mstar,cte2,
     5	jac(pnelem*pnelem),depsro,hh,be7,b8,n13,o15,f17,beta,
     6	dtaurp,dtaurt,dtaurx,dtaurr,dtaurm,depsx(1)
 
	real*8	delta,deltap,deltat,deltax,cpp,cp,dcpp,dcpt,dcpx,gradad,
     1	dgradadp,dgradadt,dgradadx,krad,dkradp,dkradt,dkradx,
     2	gravite,dgravr,dgravm,hp,dhpp,dhpt,dhpx,dhpr,dhpm,z,
     3	dgradkra,dgradgra,dgradel,dgradcp,dgradro,dgradhp,
     4	dgradgrad,dgradgad,dkapro,aradias3,dxchim(1),lamb,
     5	gradmu,d_grad,cte7,mu,taur,dgradtaur,rot,w,dw,
     6	nuc_m(pnelem),dxchimm(pnelem)
 
	logical init,deriv,ovsht
c	data init/.true./
 
	external etat,opa,conv,nuc
 
	data init/.true./
 
	save init,cte7,cte8,cte9,cte1,cte13,aradias3,cte2
 
2000	format((1x,1p8d10.3))
 
c	write(6,2000)p,t,xchim(1),m,l,r,dxchim(1),mstar
c	print*,deriv
c	pause'entree thermo_3'
c	les elements lourds	
c	Z0: abondance initiale en elements lourds Z0=1-x0-Y0
c	on impose nucleo(14)=Z0
c	Z est eventuellement diffuse pour tenir compte + correctement de
c	l'abondance des metaux dans le calcul de l'opaciteles opacites
 
	if(init)then	!initialisations
	 init=.false.
 
	 cte7=rsol/msol/g
	 cte8=lsol/4./pi/rsol/rsol	!de 5.9
	 cte9=3./16./pi/aradia/clight/g
c	 cte8=cte9*lsol/msol		!calcul de gradrad sans grv. eff.
	 cte1=4./3.*aradia*clight
	 cte13=g*msol/rsol/rsol
	 aradias3=aradia/3.
	 cte2=2./3.*rsol
 
	 if(ledoux)then
	  write(6,*)' '
	  write(6,*)'convection: critere de LEDOUX'
	  write(6,*)' '
	  write(2,*)' '
	  write(2,*)'convection: critere de LEDOUX'
	  write(2,*)' '
	 else
	  write(6,*)' '
	  write(6,*)'convection: critere de SCHWARZSCHILD'
	  write(6,*)' '
	  write(2,*)' '
	  write(2,*)'convection: critere de SCHWARZSCHILD'
	  write(2,*)' '
	 endif
	endif		!initialisation
 
	do i=1,nbelem
	 xchimm(i)=abs(xchim(i))	!composition chimique /gr
	 dxchimm(i)=dxchim(i)	
	enddo
	call chim_gram_3(xchimm,dxchimm,nuc_m)
c	write(6,*)'deriv/p,t,xchim,nucleo,xchimm',deriv
c	write(6,2000)p,t,xchim(1),nucleo(1),xchimm(1)
c	pause'3'
 
	call etat(p,t,xchimm,deriv,	!deriv=.true. calcul derivees secondes
     1	ro,drop,drot,drox,drott,drotp,drotx,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
 
	drox=drox/nuc_m(1)	!par mole
	drotx=drotx/nuc_m(1)
	dux=dux/nuc_m(1)
	dutx=dutx/nuc_m(1)
 
c	write(6,*)'ro,drop,drot,drox,drott,drotp,drotx,
c	1	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2'
c	write(6,2000)ro,drop,drot,drox,drott,drotp,drotx,
c	1	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2
c	pause'4'
 
	delta=-t/ro*drot
 
	cpp=p/ro/t*delta
	cp=dut+cpp
 
	gradad=p/ro*delta/cp/t			!gradient adiabatique
 
	if(t .gt. t_inf)then
	 call nuc(t,ro,xchim,dcomp,jac,deriv,3,
     1	epsilon,depst,depsro,depsx,hh,be7,b8,n13,o15,f17)
	else
	 epsilon(1)=0.	!total
	 epsilon(2)=0.	!pp
	 epsilon(3)=0.	!cno
	 epsilon(4)=0.	!3 alpha
	 depst=0.
	 depsro=0.
	 do i=1,nchim
	  depsx(i)=0.
	 enddo
	endif
	
	call opa(xchimm,t,ro,kap,dkapt,dkapro,dkapx)
	dkapx=dkapx/nuc_m(1)	!par mole
 
	if(deriv)then
	 deltap=delta*(-drop/ro+drotp/drot)
	 deltat=delta*(-drot/ro+drott/drot+1.d0/t)
	 deltax=delta*(-drox/ro+drotx/drot)
 
	 dcpp=dutp+cpp*(-drop/ro+deltap/delta+1.d0/p)
	 dcpt=dutt+cpp*(-drot/ro+deltat/delta-1.d0/t)
	 dcpx=dutx+cpp*(-drox/ro+deltax/delta)
 
	 dgradadp=gradad*(-drop/ro+deltap/delta-dcpp/cp+1.d0/p)
	 dgradadt=gradad*(-drot/ro+deltat/delta-dcpt/cp-1.d0/t)
	 dgradadx=gradad*(-drox/ro+deltax/delta-dcpx/cp)
 
c	 write(6,*)'delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
c	1	gradad,dgradadp,dgradadt,dgradadx'
c	 write(6,2000)delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
c	1	gradad,dgradadp,dgradadt,dgradadx
c	 pause'5'
 
	 depsp=depsro*drop
	 depst=depsro*drot+depst
	 depsx(1)=depsro*drox+depsx(1)
 
c	 write(6,*)'epsilon,depsp,depst,depsx,p,t,r,l,m'
c	 write(6,2000)epsilon(1),depsp,depst,depsx,p,t,r,l,m
 
	 dkapp=dkapro*drop
	 dkapt=dkapro*drot+dkapt
	 dkapx=dkapro*drox+dkapx
 
c	 write(6,*)'kap,dkapp,dkapt,dkapx'
c	 write(6,2000)kap,dkapp,dkapt,dkapx
c	 pause'6'
 
	endif
	
	krad=cte1/kap/ro*t**3		!5.1 conductivite radiative
	if(m*l*r .ne. 0.)then		!gradient radiatif
	 gravite=cte13*m/r**2		!5.8
	 dgravr=-2.*gravite/r
	 dgravm= gravite/m
	
	 rot=-cte2*w**2*r	!gravite effective avec rotation
	 gravite=gravite+rot	!remarque de N. Audard
	 dgravr=dgravr-3.*rot/r
	 dgravm=dgravm-2.*cte2*w*r*dw
	
	 hp=p/gravite/ro		!5.7
	 gradrad=cte8*l*hp/r**2/krad/t	!5.9	
c	 gradrad=cte8*l*p*kap/m/t**4	!sans gravite effective
	 taur=kap*ro*alpha*hp		!epaisseur optique de la bulle
 
	 if(deriv)then
	  dkradp=krad*(-drop/ro-dkapp/kap)
	  dkradt=krad*(-drot/ro-dkapt/kap+3./t)
	  dkradx=krad*(-drox/ro-dkapx/kap)
 
	  dhpp=hp*(-drop/ro+1.d0/p)
	  dhpt=-hp*drot/ro
	  dhpx=-hp*drox/ro
	  dhpr=-hp*dgravr/gravite
	  dhpm=-hp*dgravm/gravite
 
c	  dgradp=gradrad*(dkapp/kap+1/p)	!sans la gravite effective
c	  dgradt=gradrad*(dkapt/kap-4./t)
c	  dgradx=gradrad*dkapx/kap
c	  dgradm=-gradrad/m
c	  dgradr=0.
c	  dgradl=gradrad/l
	
	  dgradp=gradrad*(dhpp/hp-dkradp/krad)
	  dgradt=gradrad*(dhpt/hp-dkradt/krad-1./t)
	  dgradx=gradrad*(dhpx/hp-dkradx/krad)
	  dgradm=gradrad*dhpm/hp
	  dgradr=gradrad*(dhpr/hp-2./r)
	  dgradl=gradrad/l
 
	  dtaurp=taur*(dkapp/kap+drop/ro+dhpp/hp)
	  dtaurt=taur*(dkapt/kap+drot/ro+dhpt/hp)
	  dtaurx=taur*(dkapx/kap+drox/ro+dhpx/hp)
	  dtaurr=taur*dhpr/hp
	  dtaurm=taur*dhpm/hp
 
c	  write(6,*)'krad,dkradp,dkradt,dkradx,gravite,dgravr,dgravm,
c	1	hp,dhpp,dhpt,dhpx,dhpr,dhpm,gradrad,dgradp,dgradt,dgradx,
c	2	dgradm,dgradr,dgradl (en dehors du centre)'
c	  write(6,2000)krad,dkradp,dkradt,dkradx,gravite,dgravr,dgravm,
c	1	hp,dhpp,dhpt,dhpx,dhpr,dhpm,gradrad,dgradp,dgradt,dgradx,
c	2	dgradm,dgradr,dgradl
	 endif
 
	else	!au centre, mais approximativement
	 gradrad=cte9*kap*epsilon(1)*p/t**4	!au centre l/m ~ epsilon
	 hp=0.d0
	 gravite=0.d0
	 if(deriv)then
	  if(gradrad .ne. 0.d0)then
	   dgradp=gradrad*(dkapp/kap+depsp/epsilon(1)+1.d0/p)
	   dgradt=gradrad*(dkapt/kap+depst/epsilon(1)-4.d0/t)
	   dgradx=gradrad*(dkapx/kap+depsx(1)/epsilon(1))
	  else
	   dgradp=0.
	   dgradt=0.
	   dgradx=0.
	  endif
	  dgradm=0.
	  dgradr=0.
	  dgradl=0.
c	  write(6,*)'gradrad,dgradp,dgradt,dgradx (au centre)'
c	  write(6,2000)gradrad,dgradp,dgradt,dgradx
	 endif
	endif
 
	gradient=gradrad
	gradconv=0
 
	d_grad=gradrad-gradad	!critere de Schwarzschild pour la convection
 
c	pause'7'
 
c	critere de Ledoux pour la convection
c	comme il n'y a pas de gradient de composition chimique dans le milieu
c	froid on suppose la totale ionisation (xchim est en mole**-1)
c	z est le Z local s'il est diffuse ou Z0
 
	if(ledoux)then
	 beta=aradias3*t**4/p		!x=X/nucleo
	 if(iz .gt. 1)then
	  z=abs(xchim(iz))
	 else
	  z=z0
	 endif
	 mu=2.*abs(xchim(1))+z/2.		!1/mu=1/(2X+3/4*Y+Z/2)
	 gradmu=2.*dxchim(1)		!d 1/mu / dm
	 if(ihe4 .gt. 1)mu=mu+3.*abs(xchim(ihe4))
	 if(m .ne. 0.d0)then
	  gradmu=mu*gradmu*cte7*p*r**4/m
	 else
	  gradmu=0.
	 endif
	 d_grad=(4.d0-3.d0*beta)/beta*d_grad-gradmu	!X (-1) le critere
	endif
 
c	write(6,*)'jpz,lim/xchim(1),d_grad,r,r_zc,r_ov',jpz,lim
c	write(6,2000)xchim(1),d_grad,r,(r_zc(i),r_ov(i),i=1,lim)
c	write(6,2000)d_grad,beta,gradrad,gradad,gradmu,dxchim(1),z,mu
 
	if(d_grad .le. 0.d0)then	!zone radiative grad_rad<grad_ad
	 gradconv=0
	 ovsht=.false.		!.true. : il y a overshoot
	 do i=1,lim		!est-on dans une zone d'overshoot?
	
c	  si r_zc a ete fixe et que la limite ZR/ZC s'est legerement deplacee
c	  r est considere comme en dehors de la zone d'overshoot
c	  on ajoute/retire arbitrairement 0.01
	
	  if(r_ov(i) .ge. 0.d0)then
	   if(ovshts .gt. 0.d0)ovsht=
     1	ovsht .or. (r .le. r_ov(i) .and. r .ge. r_zc(i)-0.01)
	   if(ovshti .gt. 0.d0)ovsht=
     1	ovsht .or. (r .ge. r_ov(i) .and. r .le. r_zc(i)+0.01)
	  endif
	 enddo		!lim
c	 ovsht=.false.		!pour avoir grad=grad_rad dans zone ovsht
c	 write(6,*)ovsht,lim
c	 write(6,2000)d_grad,r,(r_zc(i),r_ov(i),i=1,lim)
	 if(ovsht)then		!.true. : il y a penetration convective
	  gradient=gradad
	  dgradp=dgradadp
	  dgradt=dgradadt
	  dgradx=dgradadx
	  dgradm=0.
	  dgradr=0.
	  dgradl=0.
	 else		!pc
	  gradient=gradrad		!pas de penetration convective
	 endif		!pc
	else				!zone convective
	 if(m*l*r .eq. 0.d0 .or. t .gt. 5.d5)then !ZC a m=0, ou profonde
	  gradconv=gradad
	  if(deriv)then
	   dgradp=dgradadp
	   dgradt=dgradadt
	   dgradx=dgradadx
	   dgradm=0.
	   dgradr=0.
	   dgradl=0.
c	   write(6,*)'gradconv,dgradp,dgradt,dgradx (ZC au centre)'
c	   write(6,2000)gradconv,dgradp,dgradt,dgradx
	  endif
	 else			!zone convective pour m .ne. 0
	  call conv(krad,gravite,delta,cp,ro,hp,taur,gradrad,gradad,deriv,
     1	gradconv,dgradkra,dgradgra,dgradel,dgradcp,dgradro,
     2	dgradtaur,dgradhp,dgradgrad,dgradgad)
	  if(deriv)then
	   dgradp=dgradkra*dkradp+dgradel*deltap+dgradcp*dcpp+dgradro*drop+
     1	dgradhp*dhpp+dgradgrad*dgradp+dgradgad*dgradadp
     2	+dgradtaur*dtaurp
	   dgradt=dgradkra*dkradt+dgradel*deltat+dgradcp*dcpt+dgradro*drot+
     1	dgradhp*dhpt+dgradgrad*dgradt+dgradgad*dgradadt
     2	+dgradtaur*dtaurt
	   dgradx=dgradkra*dkradx+dgradel*deltax+dgradcp*dcpx+dgradro*drox+
     1	dgradhp*dhpx+dgradgrad*dgradx+dgradgad*dgradadx
     2	+dgradtaur*dtaurx
	   dgradm=dgradgra*dgravm+dgradhp*dhpm+dgradgrad*dgradm
     1	+dgradtaur*dtaurm
	   dgradr=dgradgra*dgravr+dgradhp*dhpr+dgradtaur*dtaurr
	   dgradl=dgradgrad*dgradl
c	   write(6,*)'zone convective'
c	   write(6,*)'gradconv,dgradkra,dgradgra,dgradel,dgradcp,dgradro,
c	1	dgradhp,dgradgrad,dgradgad,dgradp,dgradt,dgradx,dgradm,
c	2	dgradr,dgradl'
c	   write(6,2000)gradconv,dgradkra,dgradgra,dgradel,dgradcp,dgradro,
c	1	dgradhp,dgradgrad,dgradgad,dgradp,dgradt,dgradx,dgradm,
c	2	dgradr,dgradl
	  endif		!deriv
	 endif		!en m=0
	 gradient=gradconv
	endif
c	write(6,*)'jpz,lim/xchim(1),d_grad,r,r_zc,r_ov',jpz,lim
c	write(6,2000)gradient,gradad,gradrad,d_grad,r,(r_zc(i),r_ov(i),i=1,lim)
 
c	pause'8'
 
	return
 
	end
