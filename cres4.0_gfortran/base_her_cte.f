 
c********************************************************************
 
	subroutine base_her_cte(x,z,l,her,h,h2,h3)
 
c	determination des 2 valeurs des 2 fonctions de base
c	non identiquement nulles Hermite et de
c	leurs derivees 1iere et 2nde au point z compris entre x(l) et x(l+1)
c	x(l) .le. z .lt. x(l+1)
c	le pas h est constant
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
c entree
c	x : abscisses
c	z : point de calcul
c	l : indice de couche
c	h, h2, h3 : pas spatial constant **1, **2, **3
 
c sortie
c	her(derivee,fonction phi ou psi,partie droite ou gauche)
 
c	her(0,1,1)=phi1=1 en xi, 0 en xi+1	(derivee d'ordre 0)  * fi
c	her(0,1,2)=phi2=0 en xi, 1 en xi+1	(derivee d'ordre 0)  * fi+1
c	her(1,1,1)=phi'1=0 en x, i 0 en xi+1	(derivee d'ordre 1)
c	her(1,1,2)=phi'2=0 en xi, 0 en xi+1	(derivee d'ordre 1)
 
c	her(0,2,1)=psi1=0 en xi, 0 en xi+1	(derivee d'ordre 0)  * f'i
c	her(0,2,2)=psi2=0 en xi, 0 en xi+1	(derivee d'ordre 0)  * f'i+1
c	her(1,2,1)=psi'1=1 en xi, 0 en xi+1	(derivee d'ordre 1)
c	her(1,2,2)=psi'2=0 en xi, 1 en xi+1	(derivee d'ordre 1)
 
	implicit none
 
	integer l
 
	real*8 her(0:2,2,2),x(1),z,a(0:3),p(0:3),h,h2,h3,y


2000	format((1x,1p8d10.3))
 
	y=z-x(l)
	a(1)=0.	
 
	a(0)=1.
	a(2)=-3.*h2		!hermite 1.1 en z (*fonction en l)
	a(3)=2.*h3
c	write(6,*)'l/a(0,3),h,y',l
c	write(6,2000)(a(i),i=0,3),h,y
	call polyder(a,3,2,y,p)
	her(0,1,1)=p(0)
	her(1,1,1)=p(1)
	her(2,1,1)=p(2)
 
	a(0)=0.
	a(2)=3.*h2		!hermite 1.2 en z (*fonction en l+1)
	a(3)=-2.*h3
	call polyder(a,3,2,y,p)
	her(0,1,2)=p(0)
	her(1,1,2)=p(1)
	her(2,1,2)=p(2)
 
	a(1)=1.	
	a(2)=-2.*h		!hermite 2.1 en z (*derivee en l)
	a(3)=h2
	call polyder(a,3,2,y,p)
	her(0,2,1)=p(0)
	her(1,2,1)=p(1)
	her(2,2,1)=p(2)
 
	a(1)=0.	
	a(2)=-h			!hermite 2.2 en z (*derivee en l)
	a(3)=h2
	call polyder(a,3,2,y,p)
	her(0,2,2)=p(0)
	her(1,2,2)=p(1)
	her(2,2,2)=p(2)
 
	return
 
	end
