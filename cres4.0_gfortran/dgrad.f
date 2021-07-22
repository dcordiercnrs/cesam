c****************************************************************************
 
	function dgrad(p,t,xchim,m,l,r,dxchim,mstar,w,etat,opa,nuc)
	
c	calcul de la difference des gradients pour le critere de convection
 
c	Auteur: P. Morel + N. Audard, Departement J.D. Cassini,
c	O.C.A., Observatoire de Nice
 
c	version: 3
 
c	Modif
c	09 10 96 introduction de w
 
c entree :
c	p : pression cgs
c	t : temperature K
c	xchim : composition chimique ATTENTION en 1/mole : *ah pour avoir X
c	dxchim : d xchim/d m pour critere de Ledoux uniquement
c	m : masse/msol
c	l : luminosite/lsol
c	r : rayon / rsol
c	mstar: masse au temps du calcul, avec perte de masse
c	w : rotation
 
c routines externes :
c	etat, opa, conv, nuc
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'ctephy.common'
	include 'modele_3.common'
	include 'evol_chim_3.common'
 
	integer i
 
	real*8	p,t,xchim(1),m,l,r,xchimm(pnelem),dgrad,
     1	ro,drop,drot,dkapro,beta,depsx(1),w,krad,gravite,
     2	epsilon(5),depst,u,dup,dut,drox,dux,dkapx,
     3	kap,cte1,cte2,cte8,cte9,cte13,gradrad,dkapt,nh1,nhe1,nhe2,
     4	drott,drotp,drotx,dutt,dutp,dutx,dcomp(1),mstar,
     5	jac(pnelem*pnelem),depsro,hh,be7,b8,n13,o15,f17,
     6	delta,cp,gradad,z,aradias3,dxchim(1),lamb,
     7	gradmu,mu,nuc_m(pnelem),dxchimm(pnelem)
 
	logical init
c	data init/.true./
 
	external etat,opa,nuc
 
	data init/.true./
 
	save init,cte8,cte9,aradias3
 
2000	format((1x,1p8d10.3))
 
	if(init)then	!initialisations
	 init=.false.
	 cte1=4./3.*aradia*clight
	 cte13=g*msol/rsol/rsol	 	
	 cte2=2./3.*rsol
	 cte8=lsol/4./pi/rsol/rsol	!de 5.9
	 cte9=3./16./pi/aradia/clight/g	 	
	 aradias3=aradia/3.
	endif		!initialisation
 
	do i=1,nbelem
	 xchimm(i)=abs(xchim(i))	!composition chimique /gr
	 dxchimm(i)=dxchim(i)	
	enddo
	call chim_gram_3(xchimm,dxchimm,nuc_m)
 
	call etat(p,t,xchimm,.false.,	!deriv=.true. calcul derivees secondes
     1	ro,drop,drot,drox,drott,drotp,drotx,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
 
	delta=-t/ro*drot
	cp=dut+p/ro/t*delta
	gradad=p/ro*delta/cp/t			!gradient adiabatique
 
	call opa(xchimm,t,ro,kap,dkapt,dkapro,dkapx)
	krad=cte1/kap/ro*t**3		!5.1 conductivite radiative
	
	if(m*l*r .ne. 0.)then		!gradient radiatif
	 gravite=cte13*m/r**2-cte2*w**2*r !gravite effective avec rotation
	 gradrad=cte8*l*p/gravite/ro/r**2/krad/t	!5.9		
	else		!au centre
	 if(t .gt. t_inf)then
	  call nuc(t,ro,xchim,dcomp,jac,.false.,3,
     1	epsilon,depst,depsro,depsx,hh,be7,b8,n13,o15,f17)
	 else
	  epsilon(1)=0.	!total
	 endif
	 gradrad=cte9*kap*epsilon(1)*p/t**4	!au centre l/m ~ epsilon
	endif
	dgrad=gradrad-gradad	!critere de Schwarzschild pour la convection
 
c	critere de Ledoux pour la convection
c	comme il n'y a pas de gradient de composition chimique dans le milieu
c	froid on suppose la totale ionisation (xchim est en mole**-1)
 
	if(ledoux)then
	 beta=aradias3*t**4/p		!x=X/nucleo
	 if(iz .gt. 1)then
	  z=abs(xchim(iz))
	 else
	  z=z0
	 endif
	 mu=2.*abs(xchim(1))+z/2.	!1/mu=1/(2X+3/4*Y+Z/2)
	 gradmu=2.*dxchim(1)		!d 1/mu / dm
	 if(ihe4 .gt. 1)mu=mu+3.*abs(xchim(ihe4))
	 dgrad=(4.d0-3.d0*beta)/beta*dgrad-gradmu	!X (-1) le critere
	endif
	
	return
	
	end
