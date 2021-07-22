 
c********************************************************************
 
	subroutine horner(k,der,z,a,x,p)
 
c	algorithme de horner
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
c entrees
c	k : degre du polynome
c	der : ordre de derivation
c	z : point du calcul
 
c entrees/sorties
c	a : coefficients / nouveaux coefficients
c	x : points de centrage / nouveaux points
 
c sorties
c	p : polynome et derivees
 
	implicit none
 
	integer k,der,i,j
 
	real*8 z,a(0:*),x(0:*),p(0:*),fac
 
2000	format((1x,1p8d10.3))
 
c	write(6,*)'k,der,z/a/x/p',k,der,z
c	write(6,2000)(a(i),i=0,k)
c	write(6,2000)(x(i),i=0,k)
 
	fac=1.
	do j=0,der
	 if(j .ne. 0)then
	  do i=k-1,0,-1
	   x(i+1)=x(i)
	  enddo
	  x(0)=z
	  fac=fac*j
	 endif
	 do i=k-1,0,-1
	  a(i)=a(i+1)*(z-x(i))+a(i)
	 enddo
	 p(j)=fac*a(j)
	enddo
c	write(6,2000)(p(i),i=0,k)
 
	return
 
	end
