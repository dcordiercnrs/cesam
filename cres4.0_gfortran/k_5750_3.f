 
c****************************************************************
 
	subroutine k_5750_3(tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
     1	ro_ext,dro_grav,dro_teff)
 
c	loi t(tau) de C van Veer issue du modele Kurucz 5750, 4.47, X=.7, Z=.02
c	mis sous la forme t=Teff(3/4(tau+q(tau)))**1/4, l/Hp=1.78
c	source: claude2.dat
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version 3 de CESAM
 
c entree :
c	tau : profondeur optique Rosseland
c	teff : temperature effective
c	grav : gravite
 
c sortie :
c	t : temperature
c	dtsd* : derivees t/ tau, teff, grav
c	dtsd* : derivees t/ tau, teff, grav
c	dro_** : derivees ro_ext/ teff, grav
 
	implicit none
 
	include 'atmosphere_3.common'
 
	integer pm,pmm
	parameter (pm=64, pmm=4)
        logical init
	integer knot,l,nx,mx,nrx,i
 
	real*8 tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,qs(pm),teffc,cte1,q,
     1	taus(pm),ts(pm),taut(pm+2*pmm),taur(pm),stau(pmm*pm),dqdtau,
     2	ro_ext,dro_grav,dro_teff
 
	data ts
     1/ 3.903D+03, 3.927D+03, 3.951D+03, 3.974D+03, 3.998D+03,
     + 4.021D+03,
     2  4.044D+03, 4.067D+03, 4.090D+03, 4.114D+03, 4.138D+03,
     + 4.163D+03,
     3  4.188D+03, 4.214D+03, 4.240D+03, 4.267D+03, 4.294D+03,
     + 4.322D+03,
     4  4.350D+03, 4.379D+03, 4.407D+03, 4.436D+03, 4.464D+03,
     + 4.492D+03,
     5  4.521D+03, 4.550D+03, 4.578D+03, 4.607D+03, 4.635D+03,
     + 4.664D+03,
     6  4.694D+03, 4.725D+03, 4.758D+03, 4.792D+03, 4.829D+03,
     + 4.869D+03,
     7  4.913D+03, 4.963D+03, 5.020D+03, 5.088D+03, 5.165D+03,
     + 5.256D+03,
     8  5.362D+03, 5.487D+03, 5.635D+03, 5.804D+03, 6.007D+03,
     + 6.259D+03,
     9  6.564D+03, 6.907D+03, 7.192D+03, 7.451D+03, 7.690D+03,
     + 7.904D+03,
     1  8.112D+03, 8.300D+03, 8.491D+03, 8.663D+03, 8.845D+03,
     + 9.008D+03,
     2  9.186D+03, 9.345D+03, 9.522D+03, 9.681D+03/
	data taus
     1/ 1.330D-06, 1.780D-06, 2.370D-06, 3.160D-06, 4.220D-06,
     + 5.620D-06,
     2  7.500D-06, 1.000D-05, 1.330D-05, 1.780D-05, 2.370D-05,
     + 3.160D-05,
     3  4.220D-05, 5.620D-05, 7.500D-05, 1.000D-04, 1.330D-04,
     + 1.780D-04,
     4  2.370D-04, 3.160D-04, 4.220D-04, 5.620D-04, 7.500D-04,
     + 1.000D-03,
     5  1.330D-03, 1.780D-03, 2.370D-03, 3.160D-03, 4.220D-03,
     + 5.620D-03,
     6  7.500D-03, 1.000D-02, 1.330D-02, 1.780D-02, 2.370D-02,
     + 3.160D-02,
     7  4.220D-02, 5.620D-02, 7.500D-02, 1.000D-01, 1.330D-01,
     + 1.780D-01,
     8  2.370D-01, 3.160D-01, 4.220D-01, 5.620D-01, 7.500D-01,
     + 1.000D+00,
     9  1.330D+00, 1.780D+00, 2.370D+00, 3.160D+00, 4.220D+00,
     + 5.620D+00,
     1  7.500D+00, 1.000D+01, 1.330D+01, 1.780D+01, 2.370D+01,
     + 3.160D+01,
     2  4.220D+01, 5.620D+01, 7.500D+01, 1.000D+02/
 
c	logical init
	data init/.true./
 
	save init,mx,cte1,taus,taut,taur,qs,stau,nx,nrx,knot
 
	if(init)then
	 l=1
	 mx=4
	 nx=pm
	 init=.false.
	 teffc=5750.d0		!Teff du calcul
	 cte1=(3./4.)**.25
	 do i=1,pm
	  qs(i)=(ts(i)/teffc)**4*4./3.-taus(i)
	 enddo
	 call pp1dn(1,taus,taut,taur,nx,nrx,mx,tau,knot,l,dqdtau,qs,stau,
     1	q,.false.)
	 write(2,*)
	 write(2,*)'loi t(tau,teff,grav), k_5750, log g=4.47, l/Hp=1.78'
	 write(2,*)'tau_min=1.d-4, ro_ext=3.29d-9'
	 write(2,*)' '
	 write(6,*)
	 write(6,*)'loi t(tau,teff,grav), k_5750, log g=4.47, l/Hp=1.78'
	 write(6,*)'tau_min=1.d-4, ro_ext=3.29d-9'
	 write(6,*)' '
	 if(tau_max .gt. 100.)then
	  write(6,*)'avec k_5750 CESAM prend tau_max .le. 100'
	  tau_max=100
	 elseif(tau_max .lt. 0.8)then
	  write(6,*)'avec k_5750 CESAM prend tau_max .ge. 0.8'
	  tau_max=0.8
	 endif
	 if(tau_min .ne. 1.d-4)then
	  write(6,*)'avec k_5750 CESAM prend tau_min = 1.d-4'
	  write(6,*)'pour modifier, changer tau_min et ro_ext dans k_5750,'
	 endif
	 tau_min=1.d-4
	endif
 
	call pp1dn(1,taus,taut,taur,nx,nrx,mx,tau,knot,l,dqdtau,qs,stau,
     1	q,.true.)
	t=teff*cte1*(tau+q)**.25
	dtsdtau=t/4./(tau+q)*(1.+dqdtau)
	dtsdteff=t/teff
	dtsdg=0.d0
	ro_ext=3.29d-9
	dro_grav=0.
	dro_teff=0.
 
	return
 
	end
