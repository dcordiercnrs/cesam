 
c****************************************************************
 
	subroutine conv_gong_3(krad,gravite,delta,cp,ro,hp,taur,gradrad,gradad,
     1	deriv,gradconv,dgradkra,dgradgra,dgradel,dgradcp,dgradro,
     2	dgradtaur,dgradhp,dgradgrad,dgradgad)
 
c	calcul du gradient convectif pour gong
 
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
     4	da0rad,da0ad,gam,gam1,gam2,corr
 
	logical init,deriv
	data init/.false./
 
	if(.not. init)then
	 init=.true.
	 ksi=6.172840271d-03	!1./162.			!ksi	page 13
	 alpha4=ksi*alpha**4
	 phi=9.d0/4.d0
	 write(2,1)
1	 format(//,1x,'gradient convectif par longueur de melange, ',
     1	' formulation de GONG',//)
	endif
 
	b=alpha4*gravite/krad**2*hp**3*delta*(ro*cp)**2	!5.5
	dbgra=b/gravite
	dbkra=-b*2.d0/krad
	dbhp=b*3.d0/hp
	dbdel=b/delta
	dbro=b*2.d0/ro
	dbcp=b*2.d0/cp
 
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
 
	gam=abs(a0/phi)**(1.d0/3.d0)			!5.12 initialisation
	ntour=0
18	ntour=ntour+1
	if(ntour .gt. 10)then
	 write(6,*)'pas de convergence pour 5.4'
	 stop
	endif
 
	corr=(((phi*gam+1.d0)*gam+1.d0)*gam-a0)/
     1	((3.d0*phi*gam+2.d0)*gam+1.d0)!correction de Newton Raphson
 
c	write(6,*)'gam,corr,t16,corr/gam'
c	write(6,2000)gam,corr,t16,corr/gam
 
	gam=gam-corr	 	!algorithme de Newton Raphson
	if(gam .ne. 0.d0)then
	 if(abs(corr/gam) .gt. 1.d-7)goto18		!test arret
	 gam1=gam*(gam+1.d0)/b
	 gam2=(2.d0*gam+1.d0)/((3.d0*phi*gam+2.d0)*gam+1.d0)
 
	 gradconv=gradad+gam1			!5.3
	 dgradkra=(gam2*da0kra-gam1*dbkra)/b
	 dgradgra=(gam2*da0gra-gam1*dbgra)/b
	 dgradel= (gam2*da0del-gam1*dbdel)/b
	 dgradhp= (gam2*da0hp -gam1*dbhp )/b
	 dgradro= (gam2*da0ro -gam1*dbro )/b
	 dgradcp= (gam2*da0cp -gam1*dbcp )/b
	 dgradgrad=  gam2*da0rad/b
	 dgradgad=1.d0+gam2*da0ad/b
	else
	 gradconv=gradad			!5.3
	 dgradkra=0.d0
	 dgradgra=0.d0
	 dgradel=0.d0
	 dgradhp=0.d0
	 dgradro=0.d0
	 dgradcp=0.d0
	 dgradgrad=0.d0
	 dgradgad=1.d0
	endif
	dgradtaur=0.d0
	
	return
 
	end
