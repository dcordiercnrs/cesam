 
c***********************************************************************
 
	subroutine zone_3(n,m,p,newn,newm)
 
c	determine les abscisses mnew de facon a repartir la fonction p
c	a peu pres uniformement
 
c	on utilise une interpolation lineaire inverse
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 3
 
c entree
c	n : nombre de points
c	m : abscisses
c	p : fonction de repartition
c	newn : nombre de points de la nouvelle repartition
 
c sortie
c	newm : abscisses assurant une repartition uniforme de p
 
	implicit none
 
	include 'cesam_3.parametres'
 
	integer n,i,newn,knot,l
 
	real*8 m(1),p(1),newm(1),y(pn),x(pn),bspint,xt(pqt),pas,xi
 
	logical init
 
	do i=1,n
	 y(i)=m(i)
	 x(i)=p(i)
	enddo
c	write(6,*)'x,y'
c	write(6,10)(x(i),i=1,n)
c	write(6,10)(y(i),i=1,n)
10	format((1x,1p8d10.3))
c	write(6,*)'newn',newn
	pas=(x(n)-x(1))/(newn-1)
	newm(1)=m(1)
	newm(newn)=m(n)
	init=.false.
	l=2
	do i=2,newn-1
	 xi=x(1)+pas*(i-1)
	 newm(i)=bspint(y,x,xt,n,2,knot,init,xi,l)
	 init=.true.
c	 write(6,10)xi,newm(i)
	enddo
 
	return
 
	end
