 
c**********************************************************************
 
	subroutine perte_ext(bp,q,n,qt,knot,chim,mc,nc,mct,knotc,mstar,
     1	dt,age,old_m23,new_m23,nm,new_m23t,knotm,etat,nuc)
 
 
c	routine d'interpolation m(t+dt)--->m(t) en tenant compte 	
c	de la perte de masse (mdot>0 : gain de masse, mdot<0 : perte de masse)
c	la perte de masse est concentree dans la couche nm-1 nm
 
c	on entre avec mstar=masse au temps t'
c	on sort avec mstar=masse au temps t+dt
 
c	on tient compte de la perte de masse dans la couche [n-1 n]
c	tabule l'ancienne masse en fonction de la nouvelle (en m**2/3)
 
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
	
	integer i,nm,knotm,n,knot,nc,knotc,l
	
	real*8	bp(1),q(1),qt(1),chim(1),mc(1),mct(1),age,mstar,
     1	dt,old_m23(1),new_m23(1),new_m23t(1),f(pne),df(pne)
	
	external etat,nuc
	
2000	format(1x,1p8d10.3)
	
c	extraction des masses au temps t+dt
	
	do i=1,n
	 call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q(i),l,f,df)
c	 write(6,2000)f
	 new_m23(i)=f(5)
	 old_m23(i)=sqrt(abs(f(5)))**3	!m est en fraction de Msol
c	 write(6,2000)m(i),f(5)
	enddo
	
c	calcul de mstar(t+dt), en entree mstar est la masse au temps t
 
	mstar=mstar+mdot*1.d6*dt
	
c	la perte de masse
 
	old_m23(n)=old_m23(n)-mdot*1.d6*dt
 
c	write(6,2000)new_mstar
c	write(6,2000)(old_m23(i),i=1,n)
 
c	pause'old_m23'	
 
c	on repasse en m**2/3
 
	do i=1,n
	 old_m23(i)=old_m23(i)**(2.d0/3.d0)
	enddo
	
c	suppression des masses negatives anormales
 
	old_m23(1)=0.d0
	new_m23(1)=0.d0
	nm=n
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
 
c	write(6,2000)mtot-mstar
c	print*,mtot,mstar
	
c	pause
	
	return
 
	end
