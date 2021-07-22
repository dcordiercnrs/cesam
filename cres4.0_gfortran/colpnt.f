 
c************************************************************************
 
	function colpnt(x,i,l,m,rho,init)
 
c	colpnt = abscisse du i-ieme point de collocation compris
c	entre x(l)<x<x(l+1) pour la resolution d'une equation
c	differentielle d'ordre r avec m points de collocation
c			cf. de Boor p.293
 
c	m : ordre des splines
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
	implicit none
 
	integer i,l,m,j
 
	real*8 x(1),colpnt,rho(1),xp,xm
 
	logical init
 
	if(init)goto 100
 
c	initialisations
 
	init=.true.
 
	if(m.gt.8)goto 90
	goto(10,20,30,40,50,60,70,80),m
10	rho(1)=0.d0
	goto 100
20	rho(2)=.57735 02691 89626 d0
	rho(1)=-rho(2)
	goto 100
30	rho(3)=.77459 66692 41483 d0
	rho(1)=-rho(3)
	rho(2)=0.d0
	goto 100
40	rho(3)=.33998 10435 84856 d0
	rho(2)=-rho(3)
	rho(4)=.86113 63115 94053 d0
	rho(1)=-rho(4)
	goto 100
50	rho(4)=.53846 93101 05683 d0
	rho(2)=-rho(4)
	rho(5)=.90617 98459 38664 d0
	rho(1)=-rho(5)
	rho(3)=0.d0
	goto 100
60	rho(4)=.23861 91860 83197 d0
	rho(3)=-rho(4)
	rho(5)=.66120 93864 66265 d0
	rho(2)=-rho(5)
	rho(6)=.93246 95142 03152 d0
	rho(1)=-rho(6)
	goto 100
70	rho(5)=.40584 51513 77397 d0
	rho(3)=-rho(5)
	rho(6)=.74153 11855 99394 d0
	rho(2)=-rho(6)
	rho(7)=.94910 79123 42759 d0
	rho(1)=-rho(7)
	rho(4)=0.d0
	goto 100
80	rho(5)=.18343 46424 95650 d0
	rho(4)=-rho(5)
	rho(6)=.52553 24099 16329 d0
	rho(3)=-rho(6)
	rho(7)=.79666 64774 13627 d0
	rho(2)=-rho(7)
	rho(8)=.96028 98564 97536 d0
	rho(1)=-rho(8)
	goto 100
90	do j=1,m
	rho(j)=float(j-1)/float(m-1)/2.d0-1.d0
	enddo
 
100	if(x(l+1).le.x(l))then
			write(6,*)'dans colpnt x(l+1)=x(l)',l,x(l)
			stop
	endif
	xp=(x(l+1)+x(l))/2.
	xm=(x(l+1)-x(l))/2.
	colpnt=xp+rho(i)*xm
 
	return
 
	end
