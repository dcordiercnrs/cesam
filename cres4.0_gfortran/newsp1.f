 
c*******************************************************************
 
	subroutine newsp1(x,x1t,kno1x,m1x,x2t,kno2x,m2x,c1,c2,q,ax)
 
c	transforme la spline c1 d'ordre m1x sur le reseau x1t\kno1x,
c	en la spline c2 d'ordre m2x sur le reseau x2t\kno2x
c	on n'utilise pas la base duale, ax : matrice de travail
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
	implicit none
 
	integer k,i,kno1x,kno2x,m1x,m2x,n1x,n2x,l1x,l2x,ix
 
	real*8 ax(1),c1(1),c2(1),q(1),x1t(1),x2t(1),x(1),bspint,p
 
c	calcul des tau selon de-Boor p.123
 
	n1x=kno1x-m1x	!dimension des splines
	n2x=kno2x-m2x
	l1x=m1x
	l2x=m2x
	do i=1,n2x
	 do ix=1,2*m2x-1
	  ax(n2x*(ix-1)+i)=0.
	 enddo	!ix
	enddo	!i
	do ix=1,n2x
	 p=0.d0
	 do i=1,m2x-1
	  p=p+x2t(ix+i)
	 enddo	!i
	 p=p/dble(m2x-1)
	 c2(ix)=bspint(c1,x,x1t,n1x,m1x,kno1x,.true.,p,l1x)
	 call slinf(p,x2t,kno2x,l2x)
	 call bval(p,x2t,m2x,l2x,q)
	 do k=1,m2x
	  ax(n2x*(l2x-ix+k-1)+ix)=q(k)
	 enddo	!k
	enddo	!ix
 
	call sgauss(ax,c2,n2x,m2x)
 
	return
 
	end
