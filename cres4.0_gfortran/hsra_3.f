 
c****************************************************************
 
	subroutine hsra_3(tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,
     1	ro_ext,dro_grav,dro_teff)
 
c	loi t(tau) HSRA 5780, 4.44, He/h by number .1 .Sol. Phys. 18,347
c	mis sous la forme t=Teff(3/4(tau+q(tau)))**1/4
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version 3
 
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
	parameter (pm=55, pmm=4)
 
	integer knot,l,nx,mx,nrx,i
 
	logical init
 
	real*8 tau,teff,grav,t,dtsdtau,dtsdteff,dtsdg,qs(pm),teffc,cte1,q,
     1	taus(pm),ts(pm),taut(pm+2*pmm),taur(pm),stau(pmm*pm),
     2	dqdtau,ro_ext,dro_grav,dro_teff
 
	data ts
     1/4.170d3, 4.175d3, 4.190d3, 4.205d3, 4.225d3, 4.250d3, 4.280d3,
     2 4.305d3, 4.330d3, 4.355d3, 4.380d3, 4.405d3, 4.430d3, 4.460d3,
     3 4.490d3, 4.525d3, 4.550d3, 4.575d3, 4.600d3, 4.630d3, 4.660d3,
     4 4.690d3, 4.720d3, 4.750d3, 4.790d3, 4.840d3, 4.895d3, 4.950d3,
     5 5.010d3, 5.080d3, 5.160d3, 5.240d3, 5.330d3, 5.430d3, 5.540d3,
     6 5.650d3, 5.765d3, 5.890d3, 6.035d3, 6.200d3, 6.390d3, 6.610d3,
     7 6.860d3, 7.140d3, 7.440d3, 7.750d3, 8.030d3, 8.290d3, 8.520d3,
     8 8.710d3, 8.880d3, 9.050d3, 9.220d3, 9.390d3, 9.560d3/
 
 	data taus
     1/1.000d-4, 1.259d-4, 1.585d-4, 1.995d-4, 2.512d-4, 3.162d-4, 3.981d-4
     +,
     2 5.012d-4, 6.310d-4, 7.943d-4, 1.000d-3, 1.259d-3, 1.585d-3, 1.995d-3
     +,
     3 2.512d-3, 3.162d-3, 3.981d-3, 5.012d-3, 6.310d-3, 7.943d-3, 1.000d-2
     +,
     4 1.259d-2, 1.585d-2, 1.995d-2, 2.512d-2, 3.162d-2, 3.981d-2, 5.012d-2
     +,
     5 6.310d-2, 7.943d-2, 1.000d-1, 1.259d-1, 1.585d-1, 1.995d-1, 2.512d-1
     +,
     6 3.162d-1, 3.981d-1, 5.012d-1, 6.310d-1, 7.943d-1, 1.000d+0, 1.259d+0
     +,
     7 1.585d+0, 1.995d+0, 2.512d+0, 3.162d+0, 3.981d+0, 5.012d+0, 6.310d+0
     +,
     8 7.943d+0, 1.000d+1, 1.259d+1, 1.585d+1, 1.995d+1, 2.512d+1/
 
	data init/.true./
c	logical init
 
	save init,mx,nx,taus,taut,taur,nrx,knot,qs,stau
 
	if(init)then
	 l=1
	 mx=4
	 nx=pm
	 init=.false.
	 teffc=5780.		!Teff du calcul
	 cte1=(3./4.)**.25
	 do i=1,pm
	  qs(i)=(ts(i)/teffc)**4*4./3.-taus(i)
	 enddo
	 call pp1dn(1,taus,taut,taur,nx,nrx,mx,tau,knot,l,dqdtau,qs,stau,
     1	q,.false.)
	 tau_min=1.d-4	
	 ro_ext=3.24d-9
	 write(2,*)' '
	 write(2,1)tau_min,ro_ext
1	 format(1x,'loi t(tau,teff,grav) HSRA 5780, 4.44, He/h by number .1',
     1	' Sol. Phys. 18,347, forme: t=Teff(3/4(tau+q(tau)))**1/4',/,
     2	1x,'ep. opt. 1-iere couche=',1pd10.3,' densite ext.=',1pd10.3)
	 write(6,*)' '
	 write(6,1)tau_min,ro_ext
	 if(tau_max .gt. 100.)then
	  write(6,*)'avec HSRA CESAM prend tau_max .le. 100'
	  tau_max=100
	 elseif(tau_max .lt. 0.8)then
	  write(6,*)'avec HSRA CESAM prend tau_max .ge. 0.8'
	  tau_max=0.8
	 endif
	endif
 
	call pp1dn(1,taus,taut,taur,nx,nrx,mx,tau,knot,l,dqdtau,qs,stau,
     1	q,.true.)
	t=teff*cte1*(tau+q)**.25
	dtsdtau=t/4./(tau+q)*(1.+dqdtau)
	dtsdteff=t/teff
	dtsdg=0.
	ro_ext=3.24d-9
	dro_grav=0.
	dro_teff=0.
 
	return
 
	end
