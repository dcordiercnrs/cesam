 
c******************************************************************
 
	subroutine snoein(x,y,n,m,knot)
 
c	determine la sequence de noeuds de raccord y(knot)
c	pour une interpolation "optimale"
c	par B-splines d'ordre m sur la suite strictement
c	croissante x de n points de donnee, cf. de Boor p.219 formule (10)
c	aux limites le polynome d'interpolation s'appuie sur m points de donnee
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
	implicit none
 
	integer m,knot,n,i,j
 
	real*8 x(1),y(1),mm1,eps
	data eps/1.d-12/
 
c	verification de la stricte croissance de la suite des x
 
	do i=1,n-1
	 if(x(i) .ge. x(i+1))then
	  write(6,*)'dans snoein la suite des abscisses n est pas',
     1' strictement croissante en i=',i
	  write(6,*)'nombre de points : ',n
	  write(6,*)'abscisses x=',x(i)
	  write(6,*)(x(knot),knot=1,n)
	  stop
	 endif
	enddo
 
c	pour l'interpolation spline il faut n.ge.m
 
	if(n .lt. m)then
	 write(6,11)n,m
11	 format(1x,'dans snoein n=',i3,'.lt.',i3,'=m')
	 stop
	endif
 
	mm1=m-1
 
c	formation du tableau des y
 
	knot=0
 
	do i=1,m
	 knot=knot+1
	 y(knot)=x(1)-eps
	enddo
 
	do i=1,n-m
	 knot=knot+1
	 y(knot)=0.
	 do j=i+1,i+m-1
	  y(knot)=y(knot)+x(j)
	 enddo
	 y(knot)=y(knot)/dble(mm1)
	enddo
 
	do i=1,m
	 knot=knot+1
	 y(knot)=x(n)+eps
	enddo
 
	return
 
	end
