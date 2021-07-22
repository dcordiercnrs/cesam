 
c*********************************************************************
 
	subroutine strans(a,n,b)
 
c	b = transposee de a(n,n)
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
	implicit none
 
	integer*4 n
	integer*4 i,j
 
	real*8 a(1),b(1)
 
	do i=1,n
		do j=1,n
		b(n*(j-1)+i)=a(n*(i-1)+j)
		enddo
	enddo
 
	return
 
	end
