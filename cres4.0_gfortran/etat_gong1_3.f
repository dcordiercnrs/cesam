 
c******************************************************************
 
	subroutine etat_gong1_3(p,t,xchim,deriv,
     1	ro,drop,drot,drox,drott,drotp,drotx,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
 
c	equation d'etat pour GONG etape 1
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version 3
 
c entree :
c	p : pression
c	t : temperature
c	xchim : composition chimique
c	deriv=.true. : calcul des derivees secondes
 
c sortie :
c	ro : densite et derivees
c	u : energie interne et derivees
c	nh1, nhe1, nhe2 : taux d'ionisation
c	lamb: degenerescence
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'ctephy.common'
 
	real*8 p,t,xchim(1),cte1,x,mum1,dmum1x,y,lamb,
     1	ro,drop,drot,drox,drott,drotp,drotx,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,
     3	granr,zisai,zai,unsai,zizaiz
 
	logical init,deriv
	data init/.true./
	
	save init,granr,zisai,zai,unsai,zizaiz,cte1
 
2000	format((1x,1p8d10.3))
 
c	write(6,*)'entree etat_gong1_3 p,t,xchim(1),ah,ahe4,zai',deriv
c	write(6,2000)p,t,xchim(1),ah,ahe4,zai
 
	if(init)then	!initialisations
	 init=.false.	!Z proportion en elements lourds doit avoir
			!ete initialise par appel fictif a opa
	 granr=kbol/amu
	 cte1=3./2.*granr
 
	 zisai=.5	![ Zi / Ai ]h		page 9
	 zizaiz=zisai*z0
 
	 unsai=.0625	![ 1 / Ai ]h
	 zai=(zisai+unsai)*z0	![ Zi + 1 / Ai ]h
 
	 write(2,1)
	 write(6,1)
1	 format(//,1x,'equation d''etat de GONG1 : ionisation totale, ',
     1	'pas de pression de radiation, pas de degenerescence',//)
	endif
 
	x=xchim(1)
	y=1-x-z0	!on oublie que Z est diffusable	
 
	mum1=2.*x/ah+3.*y/ahe4+zai		!2.2
	dmum1x=2./ah-3./ahe4
 
	ro=p/granr/mum1/t	!2.1
	drop= ro/p
	drot=-ro/t
	drox=-ro*dmum1x/mum1
c	write(6,2000)ro,p,mum1,t,x,ah,y,ahe4
c	write(6,2000)zai,granr,cte1
 
	u=cte1*mum1*t		!energie interne par gramme 2.3
	dup=0.
	dut=u/t
	dux=u*dmum1x/mum1
 
	if(deriv)then
	 drott=-2.*drot/t
	 drotp=-drop/t
	 drotx=-drox/t
 
	 dutp=0.
	 dutt=0.
	 dutx=dux/t
	endif
 
	nh1=1
	nhe1=0
	nhe2=1
 
	lamb=5.		!degenerescence non calculee
 
c	write(6,*)'sortie gong1 ro,drop,drot,drox,drott,drotp,drotx/u'
c	write(6,2000)ro,drop,drot,drox,drott,drotp,drotx
c	write(6,2000)u,dup,dut,dux,dutt,dutp,dutx
 
	return
 
	end
