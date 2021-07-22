 
c**********************************************************
 
	subroutine schu58_n(n,f,xx,xt,m,l,val)
 
c	interpolation pour n fonctions  par l'algorithme 5-8 p.194 de schumaker
c	au point xx, xt(l) .le. xx. lt. xt(l+1) par spline d'ordre m
c	dont les coefficients sont dans le tableau f
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
c entrees
c	n, f:nombre de fonctions f en spline
c	xx, xt: abscisse d'interpolation et vecteur nodal
c	m: ordre des splines
c	l: localisation
 
c sortie
c	val: vecteur des interpolees
 
	implicit none
 
	integer m,l,i,j,n,k
 
	real*8 xt(1),f(1),xx,cx(100),a1,a2,val(1)
 
2000	format((1x,1p8d10.3))
 
	do k=1,n
	 do j=1,m
	  cx(n*(j-1)+k)=f(n*(j+l-m-1)+k)
	  if(j+l-m .gt. 0)then
	   cx(n*(j-1)+k)=f(n*(j+l-m-1)+k)
	  else
	   cx(n*(j-1)+k)=0.d0
	  endif
	 enddo
	enddo
 
	do j=2,m
	 do i=m,j,-1
	  a1=(xx-xt(i+l-m))/(xt(i+l-j+1)-xt(i+l-m))
	  a2=1.d0-a1
	  do k=1,n
	   cx(n*(i-1)+k)=a1*cx(n*(i-1)+k)+a2*cx(n*(i-2)+k)
	  enddo
	 enddo
	enddo
 
	do k=1,n
	 val(k)=cx(n*(m-1)+k)
	enddo
 
	return
 
	end
