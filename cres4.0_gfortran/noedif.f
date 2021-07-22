c**********************************************************
 
	subroutine noedif(x,y,p,m,r,knot)
 
c	forme le tableau y des knot points de table a partir du
c	tableau x de p points de raccord pour l'integration d'une
c	equation differentielle	d'ordre r par collocation avec
c	des B-splines d'ordre m
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91, modif 01 04 93 si suite des x non strict. croissante
 
	implicit none
 
	integer m,knot,p,r,i,j
 
	real*8 x(1),y(1)
 
2000	format((1x,1p8d10.3))
 
c	verification de la stricte croissance de la suite des x
 
	do i=1,p-1
	 if(x(i) .ge. x(i+1))then
	  write(6,*)'dans noedif la suite des abscisses n''est pas',
     1	' strictement croissante en i=',i
	  write(6,2000)(x(j),j=1,i-1)
	  write(6,*)x(i-2),x(i),x(i+1)
	  write(6,2000)(x(j),j=i+1,p)
	  stop
	 endif
	enddo
 
c	formation du tableau des y
 
	knot=0
 
	do i=1,m+r
	 knot=knot+1
	 y(knot)=x(1)
	enddo
 
	do i=2,p-1
	 do j=1,m
	  knot=knot+1
	  y(knot)=x(i)
	 enddo
	enddo
 
	do i=1,m+r
	 knot=knot+1
	 y(knot)=x(p)
	enddo
 
	return
 
	end	
