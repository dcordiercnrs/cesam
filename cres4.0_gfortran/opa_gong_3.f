 
c*****************************************************************
 
	subroutine opa_gong_3(xchim,t,ro,kap,dkapt,dkapro,dkapx)
 
c	opacite du projet GONG
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 3
 
c entree :
c	xchim : composition chimique
c	t : temperature
c	ro : densite
 
c sortie :
c	kap : opacite
c	dkapt : d kap / dt
c	dkapro : d kap / dro
c	dkapx : d kap /dx
 
	implicit none
 
	real*8 t,ro,kap,dkapt,dkapro,dkapx,ce,kapme,kapne,ci,kapmi,kapni,
     1	dkero,dket,dkiro,dkit,xchim(1)
 
	real*8 t16,ke,ki
 
	logical init
	data init /.false./
 
	if(.not. init)then
	 init=.true.
 
	 if(ro .le. 0.d0)ro=1.d0
	 if(t  .le. 0.d0)t=1.d6
 
	 ce=1.6236784d-33		!ce	page 10
	 kapme=0.407895d0			!me	  "
	 kapne=9.28289d0			!ne	  "
	 ci=7.1548412d13		!ci	  "
	 kapmi=0.138316d0			!mi	  "
	 kapni=-1.97541d0			!ni	  "
	 write(6,1)
	 write(2,1)
1	 format(//,1x,'opacite de GONG : formule analytique, Z=0.02',//)
	endif
 
	t16=t
	ke=ro**kapme*ce*t16**kapne	!opacite 3.2
	dkero=ke*kapme/ro
	dket= ke*kapne/t
	ki=ci*ro**kapmi*t**kapni	!3.3
	dkiro=ki*kapmi/ro
	dkit =ki*kapni/t
	kap=ke*ki/(ke+ki)		!3.1
	dkapro=(ke**2*dkiro+ki**2*dkero)/(ki+ke)**2
	dkapt= (ke**2*dkit +ki**2*dket )/(ki+ke)**2
	dkapx=0.d0
 
	return
 
	end
