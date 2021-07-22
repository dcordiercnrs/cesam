 
c**********************************************************
 
	subroutine sbsp1d(f,x,xt,n,m,knot,init,xx,l,fx,dfxdx)
 
c	interpolation dans le tableau f(n) avec une spline polynomiale
c	d'ordre m>1, au point xx : x(1) .le. xx .le. x(n)
c	et	xt(l) .le. xx .lt. xt(l+1)
c	en entree prendre l quelconque : 1 .le. l .le. n+2m
c	au premier appel, init=.false. il y a initialisation :
c	calcul de knot=2m+n, formation du tableau xt(knot) des points de table
c	et de f(n) coefficients des B-splines
c	!!le tableau des donnees f est donc modifie!!
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
	implicit none
 
	integer*4 n,m,l,knot,i,j,k
 
	real*8 f(1),xt(1),x(1),xx,fx,dfxdx,q(11),d(10),a(4500),xx1
 
	logical init
 
2000	format((1x,1p8d10.3))
 
	if(init)goto 100
 
c	initialisations
 
	if(n*(2*m-1).gt.4500)then
	 write(6,201)m,n,m+1,m,n*(2*m-1)
 201	 format(1x,i2,'=m',i3,'=n dans sbsp1d : augmenter les tableaux',
     1	' q(',i2,'), d(',i2,'), a(',i4,')')
	 write(6,*)'changer la limite de la boucle d initialisation de a'
	 stop
	endif
	
	do i=1,4500
	a(i)=0.
	end do
 
c	formation du tableau xt
	
	call snoein(x,xt,n,m,knot)
 
c	initialisation de l
 
	l=m
 
c	determination des coefficients mis dans f
 
	do j=1,n
	call slinf(x(j),xt,knot,l)
	call sbval(x(j),xt,m,l,q,d)
	 do k=1,m
	  a(n*(l-j+k-1)+j)=q(k)
	 enddo
	enddo
	call sgauss(a,f,n,m)
 
c	localisation de xx
 
 100	if(xx .lt. xt(1))then
	 write(6,*)'dans sbsp1d le point sort de la grille xx < xt(1), n, m'
	 write(6,*)xx,xt(1),n,m
c	 write(6,2000)(x(i),i=1,n)
	 call slinf(x(1),xt,knot,l)
	 xx1=x(1)
	elseif(xx .gt. xt(knot))then
	 write(6,*)'dans sbsp1d le point sort de la grille xt(knot) < xx, n, m'
	 write(6,*)xx,xt(knot),n,m
	 call slinf(x(n),xt,knot,l)
	 xx1=x(n)
	else
	 call slinf(xx,xt,knot,l)
	 xx1=xx
	endif
	call sbval(xx1,xt,m,l,q,d)
 
c	interpolation par B-splines
 
	fx=0.
	dfxdx=0.
	do i=1,m
	fx=fx      +q(i)*f(l-m+i)
	dfxdx=dfxdx+d(i)*f(l-m+i)
	enddo
 
	return
 
	end
