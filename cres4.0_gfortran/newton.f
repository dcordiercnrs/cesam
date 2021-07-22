 
c********************************************************************
 
	subroutine newton(l,k,f,x,z,p,der)
 
c	polynome de newton
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
c entree
c	l : premier indice
c	k : degre du polynome
c	f : table des f
c	x : abscisses
c	z : point d'interpolation
c	der : ordre de la derivee
 
c sortie
c	p : valeur des derivees
 
	implicit none
 
	integer pd
	parameter (pd=100)
 
	integer l,k,j,der
 
	real*8 z,f(0:*),x(0:*),p(0:*),a(0:pd),xh(0:pd),ah(0:pd)
 
c	write(6,*)l,k,der,z
 
	call difdiv(l,k,f,x,a)
 
	do j=0,k
	 xh(j)=x(l+j)
	 ah(j)=a(j)
	enddo
c	write(6,*)(f(j+l),j=0,k)
c	write(6,*)(xh(j),j=0,k)
c	write(6,*)(ah(j),j=0,k)
	call horner(k,der,z,ah,xh,p)
 
	return
 
	end
