 
c****************************************************************
 
	subroutine lim_gong1_3(list,l,r,xchim,p,t,dpl,dpr,dtl,dtr,teff,rtot,
     1	m,dml,dmr,p_atm,t_atm,m_atm,tau,r_atm,mstar,
     2	tdetau,etat,opa)
 
c	calcul de la condition limite pour gong cas 1
c	relations Pext(l,r), Text(l,r), Mext(l,r) et derivees
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	CESAM Version 3
 
c entree
c	list=.true. : on calcule p_atm,t_atm,r_atm,tau,m_atm	
c	l :luminosite
c	r :rayon
c	xchim : composition chimique par gramme
c	mstar: masse avec perte de masse
 
c sortie
c	p : pression
c	t : temperature
c	m : masse
c	dpl, dpr, dtl, dtr, dml, dmr : derivees p,t,m / l, r
c	teff : temperature effective
c	rtot : rayon externe
c	p_atm,t_atm,r_atm,tau,m_atm : sans objet
c routines externes : tdetau,etat,opa
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'atmosphere_3.common'
	include 'ctephy.common'
	include 'evol_chim_3.common'
 
	integer ntour
 
	real*8 l,r,xchim(1),p,t,dpl,dpr,dtl,dtr,cte1,cte2,pprec,cte20,
     1	ro,nh1,nhe1,nhe2,corr,lambda,beta,teff,rtot,lamb,w,cte3,
     2	dux,drox,dkapx,u,kap,dkapp,dkapt,drop,drot,dup,dut,
     3	drott,drotp,drotx,dutt,dutp,dutx,dkapro,mstar,mstarp,
     4	m,dml,dmr,p_atm(1),t_atm(1),m_atm(1),tau(1),r_atm(1)
c	data mstarp/-1.d0/
 
	real*8 beta_g,lambda_g
	common/gong/beta_g,lambda_g
 
	logical init,list
c	data init /.false./
 
	external tdetau,etat,opa
        data init /.false./
        data mstarp/-1.d0/
 
	save init,cte1,cte2,pprec,cte20
 
2000	format((1x,1p8d10.3))
 
c	write(6,*)'entree lim_gong1'
 
	if(.not. init)then
	 init=.true.
	 lambda=6.
c	 write(6,*)'entrer beta'
c	 read(5,*)beta
	 beta=7.22
	 write(6,1)lambda,beta
1	 format(1x,'lambda=',1pd10.3,' beta=',1pd10.3)
	 write(2,2)lambda,beta
2	 format(/,t10,'pour GONG 1  lambda=',1pd10.3,' beta=',1pd10.3,//)
	 cte1=lsol/rsol/rsol/lambda/aradia/clight*4./pi/4.
	 cte20=beta*g*msol/rsol/rsol
	 cte3=2./3.*rsol*beta	
	 pprec=1.d7	!valeur provisoire de pext
	 n_atm=0
	 beta_g=beta
	 lambda_g=lambda	
	endif
 
c	print*,'lim_gong1_3',mstar
	if(mstar .ne. mstarp)then
	 mstarp=mstar
	 cte2=cte20*mstar
	endif
 
c	resolution de kap p r**2 = G m beta par iteration newton-raphson
 
	t=(cte1*l/r**2)**(1./4.)
 
c	write(6,*)'dans lim_gong1 pprec,r,l,xchim(1)'
c	write(6,2000)pprec,r,l,xchim(1)
c	write(6,*)'dans lim_gong1 p,corr,t,r,l,mtot,xchim(1),kap,dkapp'
 
	if(iw .gt. 0)then		!la rotation
	 w=xchim(iw)/r**2
	else
	 w=w_rot
	endif
	w=w**2
	p=pprec	!initialisation
	ntour=0
10	ntour=ntour+1
	if(ntour .gt. 10)then
	 write(6,*)'pas de convergence pour pext'
	 stop
	endif
 
	call etat(p,t,xchim,.true.,
     1	ro,drop,drot,drox,drott,drotp,drotx,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
	call opa(xchim,t,ro,kap,dkapt,dkapro,dkapx)
	dkapp=dkapro*drop
 
	corr=(kap*p-cte2/r**2+cte3*w*r)/(dkapp*p+kap)
	p=p-corr
c	write(6,2000)p,corr,t,r,l,mtot,xchim(1),kap,dkapp
	if(abs(corr/p) .gt. 1.d-8)goto10
 
	dtr=-t/r/2.
	dtl=t/l/4.
	dpr=-(cte2*2./r**3+cte3*w
     1   +p*(dkapt+dkapro*drot)*dtr)/kap/(1.+p*dkapp/kap)
	dpl=-p*(dkapt+dkapro*drot)*dtl /kap/(1.+p*dkapp/kap)
 
	m=mstar		!le raccord en masse est fait a l'exterieur
	dmr=0
	dml=0
 
	rtot=r
	teff=t
 
	pprec=p		!pour initialisation appel suivant
 
	return
 
	end
