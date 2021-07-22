 
c******************************************************************
 
	subroutine iben_3(t,ro,comp,dcomp,jac,deriv,fait,
     1	epsilon,depst,depsro,depsx,hhe,be7e,b8e,n13e,o15e,f17e)
 
c	fausses reactions thermonucleaires pour modele initial de
c	pre main sequence: epsilon=cT d'apres Iben
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	CESAM version 3.2
 
c entree :
c	t : temperature cgs
c	ro : densite cgs
c	comp : abondances
c	deriv=.true. : on calcule les derivees dont le jacobien
c	fait=1 : initialisation de la composition chimique
c	    =2 : derivee premiere et jacobien si deriv
c	    =3 : energie nucleaire et derivees / t
c	    =4 : production de neutrinos
c	    =5 : calcul des derivee/t et energie
c sorties
c	dcomp : derivee temporelle (unite de temps : 10**6 ans)
c	jac : jacobien (unite de temps : 10**6 ans)
c	e, et, ero, ex : energie thermonucleaire (unite de temps : s)
c			   : et derivees /t, ro ,X
c	hhe, be7e, b8e, n13e, o15e, f17e : nombre de neutrinos g/s
c	pour les reactions designees par le symbole (e pour electron +)
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'evol_chim_3.common'
 
	integer fait,i
 
	real*8 t,ro,dcomp(1),jac(1),comp(1),
     1	epsilon(5),depst,depsro,depsx(1),hhe,be7e,b8e,n13e,o15e,f17e
 
	logical deriv
 
	real*8 c_iben
	common/premain/c_iben
 
2000	format((1x,1p8d10.3))
 
	goto(100,200,300,400,200),fait
 
c	initialisations
 
100	 write(6,*)' '
	 write(6,*)'methode de Iben pour la recherche du modele initial PMS'
	 write(6,*)' '
	 t_inf=1.
c	 pause'iben'
	return
 
200	dcomp(1)=0
	if(deriv)jac(1)=0
	if(fait .ne. 5)return
 
300	epsilon(1)=c_iben*t
	epsilon(2)=0.
	do i=3,4
	 epsilon(i)=0.
	enddo
c	write(6,*)'Iben : epsilon(1),c_iben,t'
c	write(6,2000)epsilon(1),c_iben,t
 
	if(deriv)then
	 depsro=0
	 depst=c_iben
	 depsx(1)=0
	endif
 
	return
 
400	hhe=0			!neutrino
	be7e=0
	b8e=0
	n13e=0
	o15e=0
	f17e=0
 
	return
 
	end
