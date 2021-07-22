c****************************************************************
 
	subroutine edding_TESTm20(tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
     1	ro_ext,dro_grav,dro_teff)

c==========================================================
c BUT : voir l'influence de la condition limite externe en faisant
c       Tedding donnée par eddington - 20 %
c      
c D.C., 8 mai 2002
c==========================================================

c	loi t(tau) de edding mihalas stellar atmospheres (3.16) p. 55 et 72
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 3
 
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
	
	real*8 tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,cte1,
     1	ro_ext,dro_grav,dro_teff, pour
 
	logical init
	data init/.true./
 
	save cte1,init
 
	if(init)then
	 init=.false.
	 cte1=(3./4.)**.25
	 write(2,*)' '
	 write(2,*)'loi t(tau,teff,grav) de Edding, tau_min=1.d-4, ro_ext=1.d-9'
	 write(2,*)' '
	 write(6,*)' '
	 write(6,*)'loi t(tau,teff,grav) de edding, tau_min=1.d-4, ro_ext=1.d-9'
	 write(6,*)' '
	 if(tau_max .gt. 100.)then
	  write(6,*)'avec Edding CESAM prend tau_max .le. 100'
	  tau_max=100
	 elseif(tau_max .lt. 0.8)then
	  write(6,*)'avec Edding CESAM prend tau_max .ge. 0.8'
	  tau_max=0.8
	 endif
	 if(tau_min .ne. 1.d-4)then
	  write(6,*)'avec Edding CESAM prend tau_min = 1.d-4'
	  write(6,*)'pour modifier, changer tau_min et ro_ext dans Edding,'
	 endif
	 tau_min=1.d-4
	endif
 
c	write(6,*)'edding'
c	write(6,*)tau,teff,grav
c	pause'edding'

	pour = 0.80d0

	t=teff*cte1*(tau+2./3.)**.25 * pour
	dtsdtau=t/4./(tau+2./3.) * pour
	dtsdteff=t/teff
	dtsdg=0.
	ro_ext=1.d-9
	dro_grav=0.
	dro_teff=0.
 
c	write(6,*)t,dtsdtau,dtsdteff,dtsdg,cte1
c	pause'edding'
 
	return
 
	end
