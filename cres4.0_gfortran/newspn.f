 
c*******************************************************************
 
	subroutine newspn(n,x,x1t,kno1,m1,x2t,kno2,m2,s1,s2,a,indpc)
 
c	pour n fonctions
c	transforme la spline s1 d'ordre m1 sur le reseau x1t\kno1,
c	en la spline s2 d'ordre m2 sur le reseau x2t\kno2
c	on n'utilise pas la base duale
 
c	le cas ou la spline s2 a des discontinuites est envisage:
c	on on se deplace de dx legerement cote gauche, eventuellement ce
c	deplacement doit etre adapte
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 04 05 93
 
c entrees
c	n: nombre de fonctions
c	x, x1t, x2t, kno1, kno2: abcisses et vecteurs nodaux
c	m1, m2: ordre des splines
c	s1: premiere spline
c	a(knot2*2*m2),indpc(knot2): vecteurs de travail
 
 
 
	implicit none
 
	integer pn
	parameter (pn=12)
 
	integer i,kno1,kno2,m1,m2,n1,n2,l,j,n,indpc(1)
 
	real*8 a(1),s1(1),s2(1),q(10),x1t(1),x2t(1),x(1),p,previous,
     1	dx,fx(pn),dfxdx(pn),p1
	data previous/-1.d30/
 
2000	format((1x,1p8d10.3))
 
	if(n .gt. pn)then
	 write(6,*)'dans newsp1 mettre le parametre pn a',n
	 stop
	endif
 
c	dx=1.d-6		!deplacement vers la droite a adapter
	dx=1.d-12
 
c	calcul des tau selon de-Boor p.123
 
	n1=kno1-m1	!dimension des splines
	n2=kno2-m2
 
c	write(6,*)n,m1,kno1
c	do i=1,20
c	 call sbsp1dn(n,s1,x,x1t,j,m1,kno1,.true.,x(i),l,fx,dfxdx)
c	 write(6,2000)x(i),fx(1),x(i)**3
c	enddo
c	write(6,2000)(x2t(i),i=1,kno2)
c	pause
c	write(6,*)n,m1,kno1,m2,kno2,n1,n2,n2*m2,n2*n
 
	do i=1,n2*m2
	 a(i)=0.d0
	enddo	!i
	i=1
	do while(i .le. n2)
	 p=0.
	 do j=1,m2-1
	  p=p+x2t(j+i)
	 enddo	!j
	 p=p/(m2-1)
c	 write(6,*)i,p,previous
	 if(abs(p-previous) .lt. dx/2.)then
	  p1=p+dx
c	  i=i-1
	  call sbsp1dn(n,s1,x,x1t,j,m1,kno1,.true.,p1,l,fx,dfxdx)
	  do j=1,n
	   s2(n*(i-1)+j)=fx(j)
	  enddo	!j
	  call slinf(p1,x2t,kno2,l)
	  call bval(p1,x2t,m2,l,q)
	  do j=1,m2
	   a(n2*(j-1)+i)=q(j)
	  enddo	!j
	  indpc(i)=l-m2+1
c	  write(6,*)p,p1,i,l
c	  write(6,*)x2t(l-1),x2t(l),p,x2t(l+1)
 
	  p1=p+dx
	  i=i+1
	  call sbsp1dn(n,s1,x,x1t,j,m1,kno1,.true.,p1,l,fx,dfxdx)
	  do j=1,n
	   s2(n*(i-1)+j)=fx(j)
	  enddo	!j
	  call slinf(p1,x2t,kno2,l)
	  call bval(p1,x2t,m2,l,q)
	  do j=1,m2
	   a(n2*(j-1)+i)=q(j)
	  enddo	!j
	  indpc(i)=l-m2+1
c	  write(6,*)p,p1,i,l
c	  write(6,*)x2t(l-1),x2t(l),p,x2t(l+1)
c	  pause
	 else
	  previous=p
	  call sbsp1dn(n,s1,x,x1t,j,m1,kno1,.true.,p,l,fx,dfxdx)
	  do j=1,n
	   s2(n*(i-1)+j)=fx(j)
	  enddo	!j
	  call slinf(p,x2t,kno2,l)
	  call bval(p,x2t,m2,l,q)
	  do j=1,m2
	   a(n2*(j-1)+i)=q(j)
	  enddo	!j
	  indpc(i)=l-m2+1
	 endif
	 i=i+1
	enddo	!i
c	write(6,*)(indpc(i),i=1,n2)
c	do i=1,n2
c	 write(6,2000)(a(n2*(j-1)+i),j=1,m2),(s2(n*(i-1)+j),j=1,n)
c	enddo		!while
c	pause'newspn'
 
	call gausdn(a,s2,indpc,n2,m2,n)
 
	return
 
	end
 
