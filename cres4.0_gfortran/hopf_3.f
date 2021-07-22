 
c****************************************************************
 
	subroutine hopf_3(tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
     1	ro_ext,dro_grav,dro_teff)
 
c	loi t(tau) de hopf mihalas stellar atmospheres (3.16) p. 55 et 72
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version 3
 
c entree :
c	tau : profondeur optique Rosseland
c	teff : temperature effective
c	grav : gravite
 
c sortie :
c	t : temperature
c	ro_ext: densite externe
c	dtsd* : derivees t/ tau, teff, grav
c	dro_** : derivees ro_ext/ teff, grav
c	tau_ext: profondeur optique externe
 
	implicit none
 
	include 'atmosphere_3.common'
 
	integer pm,pm_
	parameter (pm=20, pm_=4)
 
	integer knot,l,nx,mx,nrx
 
	real*8 tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
     1	cte1,taus(pm),qs(pm),taut(pm+2*pm_),dqdtau,q,taur(pm),
     2	stau(pm_*pm),ro_ext,dro_grav,dro_teff
 
        logical init
 
	data taus/.0,.01,.03,.05,.1,.2,.3,.4,.5,.6,
     1	  .8,1.,1.5,2.,2.5,3.,3.5,4.,5.,100./
	data qs/.577351,.588236,.601242,.610758,.627919,
     1	.649550,.663365,.673090,.680240,.685801,
     2	.693534,.698540,.705130,.707916,.709191,
     3	.709806,.710120,.710270,.710398,.710446/
 
c	logical init
	data init/.true./
 
	save taus,taut,taur,nx,nrx,mx,knot,qs,stau,cte1,init
 
	if(init)then
	 l=1
	 mx=4
	 nx=20
	 cte1=(3./4.)**.25
	 init=.false.
	 call pp1dn(1,taus,taut,taur,nx,nrx,mx,tau,knot,l,dqdtau,qs,stau,
     1	q,.false.)
	 write(2,*)' '
	 write(2,*)'loi t(tau,teff,grav) de Hopf, tau_min=1.d-4, ro_ext=1.d-9'
	 write(2,*)' '
	 write(6,*)' '
	 write(6,*)'loi t(tau,teff,grav) de hopf, tau_min=1.d-4, ro_ext=1.d-9'
	 write(6,*)' '
	 if(tau_max .gt. 100.)then
	  write(6,*)'avec Hopf CESAM prend tau_max .le. 100'
	  tau_max=100
	 elseif(tau_max .lt. 0.8)then
	  write(6,*)'avec Hopf CESAM prend tau_max .ge. 0.8'
	  tau_max=0.8
	 endif
	 tau_min=1.d-4
	endif
 
c	write(6,*)'hopf'
c	write(6,*)tau,teff,grav
c	pause'hopf'
 
	call pp1dn(1,taus,taut,taur,nx,nrx,mx,tau,knot,l,dqdtau,qs,stau,
     1	q,.true.)
c	write(6,*)q,dqdtau
	t=teff*cte1*(tau+q)**.25
	dtsdtau=t/4./(tau+q)*(1.+dqdtau)
	dtsdteff=t/teff
	dtsdg=0.
	ro_ext=3.24d-9
	dro_grav=0.
	dro_teff=0.
 
c	write(6,*)t,dtsdtau,dtsdteff,dtsdg,cte1
c	pause'hopfs'
 
	return
 
	end
