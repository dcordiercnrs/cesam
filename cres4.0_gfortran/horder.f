 
c**********************************************************
 
	subroutine horder(m,l,a,x,p,d)
 
c	calcul de p(x) = a(l) * x**(m-1) + a(l+1) * x**(m-2) + ... + a(l+m-1)
c	et de d=dp/dx algorithme de Horner cf. burlirsh et stoer p. 270
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
	implicit none
 
	integer l,i,m
 
	real*8 d,p,x,a(1)
 
	p=a(l)
	d=p
	do i=1,m-2
	 p=p*x+a(l+i)
	 d=d*x+p
	enddo	!i
	p=p*x+a(l+m-1)
 
	return
 
	end
