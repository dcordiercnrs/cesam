c**********************************************************
 
	subroutine sbsp1dn(nf,f,x,xt,n,m,knot,init,xx,l,fx,dfxdx)
 
c	interpolation dans les tableaux f(nf,n) avec une spline polynomiale
c	d'ordre m>1, au point xx : x(1) .le. xx .le. x(n)
c	et	xt(l) .le. xx .lt. xt(l+1)
 
c	en entree prendre l quelconque : 1 .le. l .le. n+2m
c	au premier appel, init=.false. il y a initialisation :
c	calcul de knot=2m+n, formation du tableau xt(knot) des points de table
c	et de f(n) coefficients des B-splines
c	!!le tableau des donnees f est donc modifie!!
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 07 12 92
 
c entrees
c	nf: nombre de fonctions
c	f(nf,n): fonctions a interpoler
c	x: abscisses
c	n: nombre  de points
c	m: ordre des splines
c	knot: nb de noeuds
c	xx: point d'interpolation
 
c entrees/sorties
c	f(nf,n): fonctions a interpoler/coefficients des splines
c	xt: points de raccord
c	l: localisation
 
c sorties
c	f, dfdx: fonctions, derivees 1-ieres
 
	implicit none
 
	integer pn,pm
	parameter (pn=4000, pm=4)
 
	integer n,m,l,knot,i,j,nf,indpc(pn)
 
	real*8 f(1),xt(1),x(1),xx,fx(1),dfxdx(1),q(pm+1),d(pm),a(pn*pm),xx1
 
	logical init
 
2000	format((1x,1p8d10.3))
 
c	initialisations
 
	if(.not.init)then		!les coefficient sont a calculer
	 if(n .gt. pn .or. m .gt. pm)then	!test de dimensions
	  write(6,*)'dans sbsp1dn ajustement de parametres'
	  write(6,*)'mettre le parametre pm a ',m
	  write(6,*)'mettre le parametre pn a ',n
	  stop
	 endif
 
c	 calcul des B-splines
 
	 call snoein(x,xt,n,m,knot)
	 l=m
	 do i=1,n
	  call slinf(x(i),xt,knot,l)
	  call bval(x(i),xt,m,l,q)
	  do j=1,m
	   a(n*(j-1)+i)=q(j)
	  enddo	!j
	  indpc(i)=l-m+1
c	  write(6,*)indpc(i)
c	  write(6,2000)(a(n*(j-1)+i),j=1,m)
	 enddo	!i
	 call gausdn(a,f,indpc,n,m,nf)
	endif
 
c	localisation de xx
 
	if(xx .lt. xt(1))then
	 write(6,*)'dans sbsp1dn le point sort de la grille'
	 print*,'xx=',xx,' < xt(1)=',xt(1),' n, m, knot',n,m,knot
	 call slinf(x(1),xt,knot,l)
	 xx1=x(1)
	elseif(xx .gt. xt(knot))then
	 write(6,*)'dans sbsp1dn le point sort de la grille'
	 print*,'xt(knot)=',xt(knot),' < xx=',xx,' n, m, knot',n,m,knot
	 call slinf(x(n),xt,knot,l)
	 xx1=x(n)
	else
	 call slinf(xx,xt,knot,l)
	 xx1=xx
	endif
	call sbval(xx1,xt,m,l,q,d)
 
c	interpolation par B-splines
 
	do j=1,nf
	 fx(j)=0.
	 dfxdx(j)=0.
	 do i=1,m
	  fx(j)=fx(j)      +q(i)*f(nf*(l-m+i-1)+j)
	  dfxdx(j)=dfxdx(j)+d(i)*f(nf*(l-m+i-1)+j)
	 enddo
	enddo
 
	return
 
	end
