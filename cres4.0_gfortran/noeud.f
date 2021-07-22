 
c**************************************************************
 
	subroutine noeud(x,y,mi,n,knot)
 
c	formation de la suite croissante des noeuds y(knot) a partir
c	de la suite strictement croissante des points de raccordement
c	x(n)
c	mi(n) : multiplicite des noeuds definie par :
c		mi(1):=m, mi(n):=m
c	au point x(i) 1<i<n, la derivee d'ordre j=m-(mi(i)+1) est continue
c	mi(i)=m -> j=-1 (discontinuite)	* mi(1)=1 -> j=m-2=2 (pour m=4)
c	     =m-1    0   		*       2      m-3=1
c	      m-2    1	(derivee 1iere)	*       3      m-4=0
c	      m-3    2  (derivee 2de)   *       4      m-5=-1
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
	implicit none
 
	integer mi(1),knot,n,i,k
 
	real*8 x(1),y(1)
 
	knot=0
 
	do k=1,mi(1)
	 knot=knot+1
	 y(knot)=x(1)
	enddo
 
	do i=2,n-1
	 do k=1,mi(i)
	  knot=knot+1
	  y(knot)=x(i)
	 enddo
	enddo
 
	do k=1,mi(n)
	 knot=knot+1
	 y(knot)=x(n)
	enddo
 
2000	format(1x,1p8d10.3)
c	print*,'noeud',n,knot
c	print*,(x(i),i=1,n)
c	print*,(mi(i),i=1,n)
c	write(6,2000)(y(i),i=1,knot)
c	print*,(y(i),i=1,knot)
c	pause
 
	return
 
	end
