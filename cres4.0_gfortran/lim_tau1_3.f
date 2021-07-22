 
c******************************************************************
 
	subroutine lim_tau1_3(list,l,r,xchim,p,t,dpl,dpr,dtl,dtr,teff,rtot,
     1	m,dml,dmr,p_atm,t_atm,m_atm,tau,r_atm,mstar,
     2	tdetau,etat,opa)
 
c	calcul de la condition limite en tau=1
c	relations Pext(l,r), Text(l,r), Mext(l,r) et derivees
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version 3
 
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
c	n_atm : mis a 0
c routines externes : tdetau,etat,opa
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'atmosphere_3.common'
	include 'ctephy.common'
	include 'evol_chim_3.common'
 
	integer ntour
 
	real*8 l,r,xchim(1),p,t,dpl,dpr,dtl,dtr,cte1,cte2,pprec,
     1	ro,nh1,nhe1,nhe2,corr,teff,rtot,lamb,mstar,cte3,
     2	dux,drox,dkapx,u,kap,dkapp,dkapt,drop,drot,dup,dut,w,
     3	drott,drotp,drotx,dutt,dutp,dutx,dkapro,tau_b,
     4	m,dml,dmr,p_atm(1),t_atm(1),m_atm(1),tau(1),r_atm(1)
 
	logical init,list
c	data init /.true./
 
	external tdetau,opa,etat
        data init /.true./
 
	save init,cte1,cte2,pprec
 
2000	format((1x,1p8d10.3))
 
	if(init)then
	 init=.false.
	 tau_b=1.d0
	 write(6,2)tau_b
	 write(2,2)tau_b
2	 format(//,1x,'LIMITE EXTERNE sans calcul d''atmosphere ',
     1	'd P / d tau = Pext = gravite / kappa',
     2	' en tau=',1pd10.3,' : m=Mstar, r=Rstar, t=Teff',//)
	 cte1=lsol/rsol/rsol/aradia/clight*4.d0/pi/4.d0
	 cte2=g*msol/rsol/rsol*tau_b
	 cte3=2.d0/3.d0*rsol*tau_b
	 pprec=1.d5	!initialisation de la pression
	 n_atm=0
 
	endif
 
 
c	resolution de kap p r**2 = G m par iteration newton-raphson
 
	t=(cte1*l/r**2)**(1.d0/4.d0)
 
 
c	write(6,*)'dans lim_tau1 p,corr,t,r,l,mtot,xchim(1),kap,dkapp'
c	write(6,2000)p,corr,t,r,l,mtot,xchim(1),kap,dkapp
	p=pprec		!valeur provisoire de pext
	
	if(iw .gt. 0)then		!la rotation
	 w=xchim(iw)/r**2
	else
	 w=w_rot
	endif
	w=w**2
	
	ntour=0
 
10	ntour=ntour+1
	if(ntour .gt. 30)then
	 write(6,*)'pas de convergence pour pext dans lim_tau1'
	 stop
	endif
 
	call etat(p,t,xchim,.true.,
     1	ro,drop,drot,drox,drott,drotp,drotx,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
	call opa(xchim,t,ro,kap,dkapt,dkapro,dkapx)
	dkapp=dkapro*drop
 
	corr=(kap*p-cte2*mstar/r**2+cte3*w*r)/(dkapp*p+kap)
	p=p-corr
c	write(6,2000)p,corr,t,r,l,mtot,xchim(1),kap,dkapp
	if(abs(corr/p) .gt. 1.d-8)goto10
 
	dtr=-t/r/2.d0
	dtl=t/l/4.d0
	dpr=-(cte2*mstar*2./r**3+cte3*w
     1   +p*(dkapt+dkapro*drot)*dtr)/kap/(1.d0+p*dkapp/kap)
	dpl=-p*(dkapt+dkapro*drot)*dtl /kap/(1.d0+p*dkapp/kap)
 
	rtot=r
	teff=t
 
	m=mstar		!raccord en masse a l'exterieur
c	print*,'mstar',mstar
 
	dml=0
	dmr=0
 
	pprec=p		!pour initialisation appel suivant
	
 
	tau(1)=tau_b
 
	return
 
	end
