 
c**********************************************************
 
	function schu58(xx,xt,m,l,c)
 
c	interpolation par l'algorithme 5-8 p.194 de schumaker
c	au point xx xt(l) .le. x. lt. xt(l+1) par spline d'ordre m<10
c	dont les coefficients sont dans le tableau c
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
	implicit none
 
c entree
c	xx: point d'interpolation
c	xt: vecteur nodal
c	m: ordre des splines
c	l: localisation
c	c: fonction en spline
 
	integer m,l,i,j
 
	real*8 xt(1),c(1),xx,schu58,cx(11),a1,a2
 
	do j=1,m
	 if(j+l-m .gt. 0)then
	  cx(j)=c(j+l-m)
	 else
	  cx(j)=0.d0
	 endif
	enddo
 
	do j=2,m
	 do i=m,j,-1
	  a1=(xx-xt(i+l-m))/(xt(i+l-j+1)-xt(i+l-m))
	  a2=1.d0-a1
	  cx(i)=a1*cx(i)+a2*cx(i-1)
	 enddo
	enddo
 
	schu58=cx(m)
 
	return
 
	end
