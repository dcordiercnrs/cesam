 
c******************************************************************
 
	subroutine difdiv(l,k,f,x,a)
 
c	differences divisees
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
c entree
c	l : premier indice
c	k : degre du polynome
c	f : fonction
c	x : abscisses
 
c sortie
c	a : differences divisees
 
	implicit none
 
	integer pd
	parameter (pd=100)
 
	integer l,k,i,j
 
	real*8 f(0:*),x(0:*),a(0:*),d(0:pd)
 
c	write(6,*)(f(l+j),j=0,k)
c	write(6,*)(x(l+j),j=0,k)
 
	do j=0,k
	 d(j)=f(l+j)
	 do i=j-1,0,-1
	  d(i)=(d(i+1)-d(i))/(x(l+j)-x(l+i))
	 enddo
	 a(j)=d(0)
	enddo	
c	write(6,*)'les a'
c	write(6,*)(a(j),j=0,k)
 
	return
 
	end
