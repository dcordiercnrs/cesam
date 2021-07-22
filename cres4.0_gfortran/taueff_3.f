 
c**********************************************************************
 
	subroutine taueff_3(teff,grav,t_detau,tau)
 
c	on determine la valeur de tau tel que T=T(tau,teff,grav)
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version 3
 
c entrees
c	teff : temperature effective
c	grav : gravite
c	t_detau : loi T(tau)
 
c sortie
c	tau : tel que t=T(tau,teff,grav)
 
	implicit none
 
	integer tour
 
	real*8 tau,teff,epsi,grav,dtsdtau,t,dtsdteff,dtsdg,
     1	ro_ext,dro_grav,dro_teff
 
	logical ini
c	data ini/.true./
c	save ini
 
	external t_detau
 
	data ini/.true./
	save ini
 
2000	format((1x,1p8d10.3))
 
c	write(6,2000)teff,grav
	tau=2./3.
	epsi=1.d3
	tour=0
	do while (abs(epsi) .gt. 1.d-4 .and. tour .lt. 30)
	 tour=tour+1
	 call t_detau(tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
     1	ro_ext,dro_grav,dro_teff)
	 epsi=(t-teff)/dtsdtau
	 tau=tau-epsi
	 teff=t
c	 write(6,2000)tau,t,epsi,teff
	enddo
c	pause'taueff'
 
 
	if(abs(epsi) .gt. 1.d-4 .or. tau .lt. .1)then
	 if(ini)then
	  ini=.false.
	  write(6,*)'dans taueff, non CV ou tau erronne: on fixe taueff=2/3'
	 endif
	 tau=2./3.
	endif
 
	return
 
	end
