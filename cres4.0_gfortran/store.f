 
c*************************************************************************
 
	subroutine store(a,b,n)
 
c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P. Morel
 
	implicit real*8 (a-h,o-z)
 
c     stores first n elements of a into b
c
      dimension a(1),b(1)
	integer i,n
c
   10 do 11 i=1,n
   11 b(i)=a(i)
      return
 
	end
