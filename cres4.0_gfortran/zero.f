 
c********************************************************************
 
	subroutine  zero(a,nn)
 
c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P. Morel
 
	implicit real*8 (a-h,o-z)
 
c	sets a(n)=0., n=1,nn
 
	integer n,nn
 
	real*8 a(1)
 
	do n=1,nn
	 a(n)=0.
	enddo
	
	return
 
	end
