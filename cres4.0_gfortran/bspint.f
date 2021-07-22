 
c***************************************************************
 
	function bspint(f,x,xt,n,m,knot,init,xx,l)
 
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
 
	integer n,m,l,knot
	integer i,j,k, dim_dc

        parameter ( dim_dc = 10000 )
 
	real*8 f(1),xt(1),x(1),xx,bspint
	real*8 q(11),a(dim_dc),schu58	!utiliser des parametres
 
	logical init
 
	if(init)goto 100
 
c	initialisations
 
	if(n*(2*m-1) .gt. dim_dc)then
	 write(6,201)m,n,m+1,n*(2*m-1)
201	 format(1x,i2,'=m ',i4,'=n dans bspint : augmenter les tableaux',
     1	'q(',i2,'), a(',i4,')')
	 write(6,*)'changer la limite de la boucle d initialisation de a'
	 stop
	endif
	
	do i=1,dim_dc
	 a(i)=0.
	enddo
 
c	formation du tableau xt
	
	call snoein(x,xt,n,m,knot)
 
c	initialisation de l
 
	l=m
 
c	determination des coefficients mis dans f
 
	do j=1,n
	 call slinf(x(j),xt,knot,l)
	 call bval(x(j),xt,m,l,q)
	 do k=1,m
	  a(n*(l-j+k-1)+j)=q(k)
	 enddo
	enddo	!plus simple en utilisant gausdp cf. pp1d
	call sgauss(a,f,n,m)
 
c	localisation de xx
 
 100	if(xx .lt. xt(1))then
	 write(6,*)'dans bspint le point sort de la grille xx < x(1)'
	 write(6,*)xx,xt(1)
	 l=1
	elseif(xx .gt. xt(knot))then
	 write(6,*)'dans bspint le point sort de la grille x(n) < xx'
	 write(6,*)xt(knot),xx
	 l=knot-m
	else
	 call slinf(xx,xt,knot,l)
	endif
 
c	interpolation par B-splines
 
	bspint=schu58(xx,xt,m,l,f)
 
	return
	
	end
