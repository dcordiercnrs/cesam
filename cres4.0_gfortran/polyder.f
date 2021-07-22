 
c************************************************************************
			
	subroutine polyder(a,n,der,x,p)
 
c	calcul, au point x, des derivees jusqu'a l'ordre der du polynome
c	a0 + a1 x + ... +an x**n algorithme de horner
 
c	le tableau a est conserve
 
c	derivees dans p(0), p(1), ... , p(der)
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
	implicit none
	
	integer n, der, i, j
 
	real*8 a(0:*), p(0:*), fac, x
 
	do i=0,n
	 p(i)=a(i)
	enddo	!i
 
	do j=0,der
	 do i=n-1,j,-1
	  p(i)=p(i+1)*x+p(i)
	 enddo	!i
	enddo	!j
 
	if(der .le. 1)return
 
	fac=1.d0
	do i=2,der
	 fac=fac*i
	 p(i)=p(i)*fac
	enddo	!i
 
	return
 
	end
