 
c**********************************************************************
 
	subroutine perte_tot(bp,q,n,qt,knot,chim,mc,nc,mct,knotc,mstar,
     1	dt,age,old_m23,new_m23,nm,new_m23t,knotm,etat,nuc)
 
 
c	routine d'interpolation m(t+dt)--->m(t) en tenant compte 	
c	des variations de masse dues a E=mc**2
c	de la perte de masse (mdot>0 : gain de masse, mdot<0 : perte de masse)
c	la perte de masse est concentree dans la couche nm-1 nm
 
c	on entre avec mstar=masse au temps t'
c	on sort avec mstar=masse au temps t+dt
 
c	on integre les variations dues a E=mc**2
c	d'ou la variation de masse totale et la nouvelle masse totale mstar
c	calcule la masse qui corrspondait au temps t a la masse en chaque point
c	tient compte de la perte de masse dans la couche [nm-1 nm]
c	tabule l'ancienne masse en fonction de la nouvelle (en m**2/3)
c	on neglige la variation de masse de l'atmosphere entre les
c	instants t+dt et t
 
c	utilisation par sbsp1dn (m**2/3 ---> m**2/3 ancien)
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
c	version 3
 
c entrees
c	bp,q,n,qt,knot,chim,mc,nc,mct,knotc : modele au temps t
c	dt : pas temporel
c	age : age
 
c entree/sortie
c	mstar : masse totale avec perte de masse
c	en entree au temps t, en sortie au temps t+dt
 
c sortie
c	old_m23,new_m23,nm,new_m23t,knotm : interpolation de l'ancienne masse
c	en fonction de la nouvelle (en m**2/3) normalise (Mstar x Msol)
 
	include 'cesam_3.parametres'
	include 'ctephy.common'
	include 'evol_chim_3.common'
	include 'modele_3.common'
	
	integer i,j,nm,knotm,n,knot,nc,knotc,l
	
	real*8	bp(1),q(1),qt(1),chim(1),mc(1),mct(1),p,t,dm(pn),
     1	dt,old_m23(1),new_m23(1),new_m23t(1),f(pne),df(pne),cte2,
     2	m(pn),cte1,xchim(pnelem),dxchim(pnelem),ro,drop,drot,drox,
     3	drott,drotp,drotx,u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,
     4	nhe2,lamb,jac(1),epsilon(5),depst,depsro,mt(pn+6),age,mstar,
     5	depsx(pnelem),hh,be7,b8,n13,o15,f17,xchimm(pnelem),
     6	mev,ehh,ebe7,eb8,en13,eo15,ef17,nuc_m(pnelem)
	
	logical init
c	data init/.true./
 
	save init
	
	external etat,nuc
	
	data init/.true./
 
2000	format(1x,1p8d10.3)
	
	if(init)then
	 init=.false.
	 cte1=secon6/clight/clight
	 mev=1.d6*eve	!energie des neutrinos selon Clayton p.380 et 392
	 ehh=mev*0.263
	 ebe7=mev*0.80
	 eb8=mev*7.2
	 en13=mev*0.71
	 eo15=mev*1.
	 ef17=mev*0.94
	endif
	cte2=cte1*dt
 
c	extraction des masses
c	calcul des reactions thermonucleaires puis des pertes par neutrinos
c	m(i)=delta m du a e=mc2 pendant dt en fraction de Mstar x Msol
	
	dm(1)=0.
	m(1)=0	
	new_m23(1)=0
	do i=2,n
	 call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q(i),l,f,df)
c	 write(6,2000)f
	 p=exp(f(1))
	 t=exp(f(2))
	 new_m23(i)=f(5)
	 m(i)=sqrt(abs(f(5)))**3	!m est en Msol
c	 write(6,2000)m(i),f(5)
	 call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,min(f(5),mc(nc)),l,xchim,dxchim)
	 if(t .gt. t_inf)then
	  do j=1,nbelem
	   xchimm(j)=xchim(j)
	  enddo
	  call chim_gram_3(xchimm,dxchim,nuc_m)
	  call etat(p,t,xchimm,.false.,
     1	ro,drop,drot,drox,drott,drotp,drotx,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
	  call nuc(t,ro,xchim,dxchim,jac,.false.,3,	!3 : epsilon
     1	epsilon,depst,depsro,depsx,hh,be7,b8,n13,o15,f17)
	
c	  hh,be7,b8,n13,o15,f17 : nombre de neutrinos par unite de masse
c	  et de temps pour les diverses reactions
	
	  call nuc(t,ro,xchim,dxchim,jac,.false.,4,	!4 : neutrinos
     1	epsilon,depst,depsro,depsx,hh,be7,b8,n13,o15,f17)
	  dm(i)=(epsilon(1)+hh*ehh+be7*ebe7+b8*eb8+n13*en13+
     1	o15*eo15+f17*ef17)*cte2
	 else
	  dm(i)=0.
	 endif	
c	 write(6,2000)m(i),f(5),dm(i)
	enddo
	
c	write(6,2000)dt,cte2,mstar
c	write(6,2000)(dm(i),i=1,n)
c	write(6,2000)(m(i),i=1,n)
 
c	suppression des inversions, elles sont anormales avec perte de masse
c	mais peuvent se produire avec gain de masse
 
	nm=n
	i=nm-1
	do while(i .gt. 1)
	 if(m(i) .ge. m(i+1))then
c	  print*,'perte_gu, inversion en',i,nm,m(i),m(i+1)
c	  pause
	  nm=nm-1
	  do l=i,nm
	   m(l)=m(l+1)
	   dm(l)=dm(l+1)
	   new_m23(l)=new_m23(l+1)
	  enddo
	  i=min(i,nm-1)
	 else
	  i=i-1
	 endif
	enddo	
		
c	spline de dm en fonction de m pour integration
	
	call sbsp1dn(1,dm,m,mt,nm,2,
     1	knotm,.false.,m(1),l,xchim,dxchim)
	
c	integration mdot<0, pour perte de masse	
c	mstar(t+dt)=mstar(t)+m_dot*dt-defaut de masse
c	mstar contient l'atmosphere alors que m(nm) ne la contient pas
	
	call sum_n(1,dm,mt,2,knotm,.false.,m(1),m(nm),f)
	
	mstar=mstar+mdot*1.d6*dt-f(1)	!mstar au temps t en entree
	
c	on enleve la masse perdue par perte de masse due a e=mc2
 
	do i=1,nm
	 call sum_n(1,dm,mt,2,knotm,.true.,m(1),m(i),f)
	 old_m23(i)=m(i)+f(1)
c	 write(6,2000)f(1)
	enddo
	
c	la perte de masse externe
	
	old_m23(nm)=old_m23(nm)-mdot*1.d6*dt
	
c	write(6,2000)(m(i),i=1,n)
c	write(6,2000)(old_m23(i),i=1,n)
	
c	on repasse en m**2/3
 
	do i=1,nm
	 old_m23(i)=old_m23(i)**(2.d0/3.d0)
	enddo
	
c	suppression des masses negatives anormales
 
	old_m23(1)=0.d0
	new_m23(1)=0.d0
	i=2
	do while(i .le. nm)
	 if(old_m23(i) .le. 0.d0 .or. new_m23(i) .le. 0.d0)then
	  nm=nm-1
	  do l=i,nm
	   old_m23(l)=old_m23(l+1)
	   new_m23(l)=new_m23(l+1)
	  enddo
	 else
	  i=i+1
	 endif
	enddo
	
c	suppression des inversions, elles sont anormales avec perte de masse
c	mais peuvent se produire avec gain de masse
 
	i=nm-1
	do while(i .gt. 1)
	 if(old_m23(i) .ge. old_m23(i+1) .or. new_m23(i) .ge. new_m23(i+1))then
c	  print*,'perte_ext, inversion en',i,nm,
c	1	old_m23(i),old_m23(i+1),new_m23(i),new_m23(i+1)
c	  pause
	  nm=nm-1
	  do l=i,nm
	   old_m23(l)=old_m23(l+1)
	   new_m23(l)=new_m23(l+1)
	  enddo
	  i=min(i,nm-1)
	 else
	  i=i-1
	 endif
	enddo	
	
c	tabulation
 
	call sbsp1dn(1,old_m23,new_m23,new_m23t,nm,2,knotm,.false.,
     1	new_m23(1),l,f,df)
c	print*,mstar,dt,age
 
c	print*,'mtot-mstar/mtot,mstar'	
c	write(6,2000)mtot-mstar
c	print*,mtot,mstar
c	pause'sortie perte_ext'
	
	return
 
	end
 
 
