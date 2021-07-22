 
c********************************************************************
 
	function omega(i,iz)
 
c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P. Morel
 
      implicit real*8 (a-h,o-z)
 
c      implicit real*16 (a-h,o-z)
c  calculates statistical weight of i-th ionization stage of element
c  with number iz
      dimension iom(26),iom1(20)
 
	real*8 omega
	integer i,iz,iom,iom1
 
      data iom/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,10,21,28,25,6,25,
     .   28,21/
      data iom1/2,1,1,2,10,15,21,28,28,25,7,6,6,7,25,30,28,21,21,10/
c
      if(i.le.1.and.iz.ge.19) go to 20
      if(i.eq.iz) go to 15
      omega=iom(iz-i)
      return
   15 omega=15
      return
   20 omega=iom1(2*(iz-19)+i+1)
      return
 
	end
