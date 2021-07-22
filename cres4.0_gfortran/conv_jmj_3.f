c****************************************************************
 
	subroutine conv_jmj_3(krad,gravite,delta,cp,ro,hp,taur,gradrad,gradad,
     1	deriv,gradconv,dgradkra,dgradgra,dgradel,dgradcp,dgradro,
     2	dgradtaur,dgradhp,dgradgrad,dgradgad)
 
c	calcul du gradient convectif
 
c	formulation de:
c	J. Provost, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	M. J. Goupil DASGAL Observatoire de Meudon
c	Correction erreur sur b S. Brun 30 08 96
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 3
 
c entree :
c	krad : conductivite radiative 4 ac T**3 / 3 kap ro
c	gravite : G m / r**2
c	delta : - ( d ln ro / d ln T )p
c	cp : chaleur specifique
c	ro : densite
c	hp : echelle de hauteur de pression
c	gradrad : gradient radiatif
c	gradad : gradient adiabatique
c	taur : ep. opt. de la bulle
c	deriv : on calcule les derivees
 
c	alpha : longueur de melange (dans le common /modele/)
 
c sortie :
c	gradconv : gradient convectif
c	dgrad.. : derivees
 
	implicit none
 
	include 'modele_3.common'
 
	integer ntour
 
	real*8 krad,gravite,delta,cp,ro,hp,gradrad,gradad,taur,dgradtaur,
     1	gradconv,dgradkra,dgradgra,dgradel,dgradcp,dgradro,dgradhp,
     2	dgradgrad,dgradgad,ksi,alpha4,phi,b,dbgra,dbkra,dbhp,dbdel,
     3	dbro,dbcp,gradmd,a0,da0gra,da0kra,da0hp,da0del,da0ro,da0cp,
     4	da0rad,da0ad,gam,gam1,gam2,corr,da0taur,dbtaur,dal4taur,
     5	dksitaur,ksi0,dvsltaur,vsal32,dphitaur,phi0,dvst,dvstaur
 
	logical init,deriv
	data init/.false./
	
	save ksi0,phi0
 
	if(.not. init)then	!vsal: V/Al
	 init=.true.
	 vsal=2./9.
	 ksi0=1./72.*(3.*vsal)**2			!ksi=1/72 (3V/Al)**2
	 phi0=3./2./(3.*vsal)				!phi=3/2 / (3V/Al)
	 write(2,1)vsal,ksi0,phi0
	 write(6,1)vsal,ksi0,phi0
1	 format(//,1x,'-------------------------------------------',//,
     1	1x,'gradient convectif calcule par longueur de melange: ',/,
     2	1x,'formulation tenant compte de l''ep. op. de la bulle',
     3	/,1x,' V/Al= ',1pd10.3,' ksi0=',1pd10.3,' phi0=',1pd10.3,//)
	endif
 
	dvst=1.+2./3./vsal/taur**2
	dvstaur=-2.*(dvst-1.)/taur
	phi=phi0/dvst
	dphitaur=-phi/dvst*dvstaur
 
	vsal32=dvst**2
	dvsltaur=2.*dvst*dvstaur
	ksi=ksi0*vsal32
	dksitaur=ksi0*dvsltaur
 
	alpha4=ksi*alpha**4
	dal4taur=alpha4/ksi*dksitaur
 
	b=alpha4*gravite/krad**2*hp**3*delta*(ro*cp)**2
	dbgra=b/gravite
	dbkra=-b*2./krad
	dbhp=b*3./hp
	dbdel=b/delta
	dbro=b*2./ro
	dbcp=b*2./cp
	dbtaur=b*dal4taur/alpha4
 
c	solution de 5-4 par methode iterative
 
c	write(6,*)'ksi,gravite,delta,ro,cp,krad,hp,b,
c	1	gradrad,gradad'
c	write(6,2000)ksi,gravite,delta,ro,cp,krad,hp,b,
c	1	gradrad,gradad
2000	format((1x,1p8d10.3))
	
	gradmd=gradrad-gradad
	a0=b*gradmd			!second membre de 5.4
	da0gra=dbgra*gradmd
	da0kra=dbkra*gradmd
	da0hp=dbhp*gradmd
	da0del=dbdel*gradmd
	da0ro=dbro*gradmd
	da0cp=dbcp*gradmd
	da0rad=b
	da0ad=-b
	da0taur=dbtaur*gradmd
 
	gam=abs(a0/phi)**(1./3.)			!5.12 initialisation
	ntour=0
18	ntour=ntour+1
	if(ntour .gt. 10)then
	 write(6,*)'pas de convergence pour 5.4'
	 stop
	endif
 
	corr=(((phi*gam+1.)*gam+1.)*gam-a0)/
     1	((3.*phi*gam+2.)*gam+1.)!correction de Newton Raphson
 
c	write(6,*)'gam,corr,t16,corr/gam'
c	write(6,2000)gam,corr,t16,corr/gam
 
	gam=gam-corr	 	!algorithme de Newton Raphson
	if(gam .ne. 0.d0)then
	 if(abs(corr/gam) .gt. 1.d-7)goto18		!test arret
	 gam1=gam*(gam+1.)/b
	 gam2=(2.*gam+1.)/((3.*phi*gam+2.)*gam+1.)
 
	 gradconv=gradad+gam1			!5.3
	 dgradkra=(gam2*da0kra-gam1*dbkra)/b
	 dgradgra=(gam2*da0gra-gam1*dbgra)/b
	 dgradel= (gam2*da0del-gam1*dbdel)/b
	 dgradhp= (gam2*da0hp -gam1*dbhp )/b
	 dgradro= (gam2*da0ro -gam1*dbro )/b
	 dgradcp= (gam2*da0cp -gam1*dbcp )/b
	 dgradgrad=  gam2*da0rad/b
	 dgradgad=1.+gam2*da0ad/b
	 dgradtaur=(gam2*(da0taur-gam**3*dphitaur) -gam1*dbtaur)/b
	else
	 gradconv=gradad			!5.3
	 dgradkra=0.
	 dgradgra=0.
	 dgradel=0.
	 dgradhp=0.
	 dgradro=0.
	 dgradcp=0.
	 dgradgrad=0.
	 dgradgad=1.
	 dgradtaur=0.
	endif
 
	return
 
	end
