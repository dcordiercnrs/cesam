 
c******************************************************************
 
	subroutine nuc_gong_3(t,ro,comp,dcomp,jac,deriv,fait,
     1	epsilon,depst,depsro,depsx,hhe,be7e,b8e,n13e,o15e,f17e)
 
c	reactions thermonucleaires pour gong
c	avec vecteur de comp. chim. generalisee
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version 19 07 91
 
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
 
c initialisation de COMMON
c	nbelem : nombre d'elements variables /modele/
c	nreac : nombre de reactions thermonucleaires utilises
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'evol_chim_3.common'
	include 'ctephy.common'
 
	integer fait,i
 
	real*8 t,ro,a11,b,re,dcomp(1),jac(1),comp(1),t913,
     1	epsilon(5),depst,depsro,depsx(1),hhe,be7e,b8e,n13e,o15e,f17e,
     2	cte1,cte2,cte3,cte5,cte6,qe,t9
 
	logical init,deriv
	data init/.false./
 
	save init,re,a11,b,cte1,cte2,qe,cte3,cte5,cte6
 
	if(.not. init)then
	 init=.true.
	 re=6.5d0
	 a11=4.21d-15
         b=3.6d0
	 cte1=-re*a11*secon6/2.d0*ah	!ah pour coherence avec comp
	 cte2=2.d0*cte1
	 qe=6.5d-5			!Qe	  "
	 cte3=1.d-9
	 cte5=1.d0/3.d0
	 cte6=qe*a11/2.d0/amu*ah**2
	 write(2,1)
	 write(6,1)
1	 format(//,1x,'reactions nucleaires simplifiees de GONG',/)
	endif
 
	goto(100,200,300,400,200),fait
 
c	initialisations
 
100	nchim=1		!1 seul element
	nom_elem(1)=' H '
	ab_min(1)=1.d-3
	nucleo(1)=ah
	comp(1)=x0/nucleo(1)
	ab_ini(1)=x0
	
	nreac=1	
	t_inf=1.d6
	ihe4=-100
 
	write(2,*)'Reaction thermonucleaire du cycle PP simplifie'
	write(2,*)' '
	write(2,*)'nombre de reactions : ',nreac
	write(2,*)'nombre d''elements chimiques : ',nchim
	write(2,*)' '
	write(2,3)(comp(i),i=1,nchim)
3	format(3x,'Abondance initiale H : ',1pd10.3)
2	format(1x,'abondance negligeable H :', 1pd10.3)
	write(2,2)(ab_min(i),i=1,nchim)
	write(2,*)' '
	write(2,*)'pour l''evolution temporelle, test de precision sur H'
	write(2,*)' '
	
	write(6,*)'Reaction thermonucleaire du cycle PP simplifie'
	write(6,*)' '
	write(6,*)'nombre de reactions : ',nreac
	write(6,*)'nombre d''elements chimiques : ',nchim
	write(6,*)' '
	write(6,3)(comp(i),i=1,nchim)
	write(6,2)(ab_min(i),i=1,nchim)
	write(6,*)' '
	write(6,*)'pour l''evolution temporelle, test de precision sur H'
	write(6,*)' '
 
	return
 
200	if(t .lt. t_inf)then	!si t<t_inf
	 do i=1,nbelem
	  dcomp(i)=0.d0
	  depsx(i)=0.d0
	 enddo
	 do i=1,nbelem*nbelem
	  jac(i)=0.d0
	 enddo
	 do i=1,4
	  epsilon(i)=0.d0
	 enddo
	 depst=0.d0
	 depsro=0.d0
	 return
	endif
	
	t913=(t*1.d-9)**(1.d0/3.d0)			!evolution chimique
	dcomp(1)=cte1*comp(1)**2*ro/t913**2*exp(-b/t913)	!f=aX**2
	if(iw .gt. 0)dcomp(iw)=0.d0	!pour MA
	if(iz .gt. 0)dcomp(iz)=0.d0	!pour Z		
	if(deriv)then
	 do i=1,nbelem*nbelem
	  jac(i)=0.d0
	 enddo		
	 jac(1)=dcomp(1)*2.d0/comp(1)		!jac=d y'/dX=2*aX
	endif
		
c	pour Mw ou Z jac=0
	
	if(fait .ne. 5)return
300	if(t .le. t_inf)then
	 do i=1,4
	  epsilon(i)=0.d0
	 enddo
	 depst=0.d0
	 depsro=0.d0
	 depsx(1)=0.d0
	else
	 t9=t*cte3			!energie nucleaire
	 t913=t9**cte5
	 epsilon(1)=0.d0
	 epsilon(2)=cte6*comp(1)**2*ro/t913/t913*exp(-b/t913)	!4.1
	 do i=3,4
	  epsilon(i)=0.d0
	 enddo
	 do i=2,4
	  epsilon(1)=epsilon(1)+epsilon(i)
	 enddo
	 if(deriv)then
	  depsro=epsilon(1)/ro
	  depst=epsilon(1)*((b/t913-2.d0)/3.d0/t9*cte3)
	  depsx(1)=epsilon(1)*2.d0/comp(1)
	 endif
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
