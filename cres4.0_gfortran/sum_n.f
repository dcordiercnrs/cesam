 
c*************************************************************************
 
	subroutine sum_n(n,f,xt,m,knot,init,c,d,sum)
 
c	somme de c a d des n fonctions f mise sous forme de spline
c	(par sbsp1dn par exemple) aux n points
c	x(1)<...<x(n), c avec d dans [x(1), x(n)]
 
c	!!! ATTENTION le tableau des splines f est modifie !!!!
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
c entrees
c	n: nombre de fonctions dans f
c	xt, m, knot: vecteur nodal, ordre des splines et dimension
c	c, d: intervalle d'integration
c	init=.true.: le tableau des f est deja adapte a l'integration
 
c entrees/sorties
c	f(n,knot-m): coefficients des splines modifie si init=.false.
 
c sorties
c	sum: vecteur des integrales
 
	implicit none
 
	integer n,m,knot,l,i,j
 
	real*8 f(1),sum(1),c,d,xt(1),val(20),bid
 
	logical init
 
2000	format((1x,1p8d10.3))
 
c	l'integrale de la spline est une spline d'ordre m+1
c	calcul des ci !!!! mis dans f !!!! algorithme 5.19 de Schumaker
 
	if(.not.init)then	!calcul de la spline de l'integrale
	 do j=1,n
	  bid=0.
	  do i=1,knot-m
	   bid=bid+(xt(i+m)-xt(i))*f(n*(i-1)+j)
	   f(n*(i-1)+j)=bid
	  enddo
	 enddo
	endif
 
c	calcul de somme de c -> d de f(x) dx
 
	call slinf(d,xt,knot,l)
	call schu58_n(n,f,d,xt,m+1,l,sum)
	call slinf(c,xt,knot,l)
	call schu58_n(n,f,c,xt,m+1,l,val)
	do i=1,n
	 sum(i)=(sum(i)-val(i))/dfloat(m)
	enddo
 
	return
 
	end
